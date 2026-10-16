//! Nightly CI ledger writer (foundations CI/CD bead; golden-ledger
//! convention, plan §16.2): records a nightly gate run as an op + metrics in
//! an fs-ledger database that CI retains as an artifact — the nightly run is
//! itself provenance, queryable like any other study.
//!
//! Usage:
//!   nightly_ledger <db-path> <outcome: ok|error> <suite> <value>
//!
//! The Five Explicits come from the CI invocation. `GITHUB_SHA` and
//! `RUNNER_OS` are self-reported process-environment provenance: a missing
//! variable is represented explicitly as unavailable, while a present but
//! malformed value is refused. No provenance value is fabricated.
//!
//! Exit code 0 on success; 1 with a structured JSON error line on failure.

use std::ffi::{OsStr, OsString};
use std::fmt;

use fs_ledger::{
    EventRow, FiveExplicits, Ledger, LedgerError, MAX_OP_CAPABILITY_BYTES, MAX_OP_FIELD_BYTES,
    MAX_OP_IR_BYTES, MAX_OP_VERSIONS_BYTES, OpOutcome, now_wall_ns,
};

const MAX_DB_PATH_BYTES: usize = 4 * 1024;
const MAX_OUTCOME_BYTES: usize = "error".len();
const MAX_SUITE_BYTES: usize = 64 * 1024;
const MAX_VALUE_TEXT_BYTES: usize = 256;
const MAX_PROVENANCE_VALUE_BYTES: usize = 128 * 1024;
const MAX_EVENT_PAYLOAD_BYTES: usize = MAX_OP_FIELD_BYTES;

#[derive(Debug, Clone, PartialEq, Eq)]
enum NightlyError {
    Input {
        field: &'static str,
        problem: String,
    },
    Ledger(LedgerError),
    TransactionCleanup {
        primary: LedgerError,
        rollback: LedgerError,
    },
}

impl NightlyError {
    fn input(field: &'static str, problem: impl Into<String>) -> Self {
        Self::Input {
            field,
            problem: problem.into(),
        }
    }
}

impl From<LedgerError> for NightlyError {
    fn from(value: LedgerError) -> Self {
        Self::Ledger(value)
    }
}

impl fmt::Display for NightlyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Input { field, problem } => {
                write!(f, "NightlyInput: field '{field}' rejected: {problem}")
            }
            Self::Ledger(error) => write!(f, "{error}"),
            Self::TransactionCleanup { primary, rollback } => write!(
                f,
                "NightlyTransactionCleanup: primary[{}]={primary}; rollback[{}]={rollback}",
                primary.code(),
                rollback.code()
            ),
        }
    }
}

fn escape_json_string(value: &str) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";

    let mut escaped = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\u{0008}' => escaped.push_str("\\b"),
            '\u{000c}' => escaped.push_str("\\f"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            '\u{0000}'..='\u{001f}' => {
                let byte = ch as u8;
                escaped.push_str("\\u00");
                escaped.push(char::from(HEX[usize::from(byte >> 4)]));
                escaped.push(char::from(HEX[usize::from(byte & 0x0f)]));
            }
            _ => escaped.push(ch),
        }
    }
    escaped
}

fn error_json(error: &NightlyError) -> String {
    format!(
        "{{\"error\":\"NightlyLedger\",\"detail\":\"{}\"}}",
        escape_json_string(&error.to_string())
    )
}

fn fail(error: &NightlyError) -> std::process::ExitCode {
    eprintln!("{}", error_json(error));
    std::process::ExitCode::FAILURE
}

fn utf8<'a>(field: &'static str, value: &'a OsStr) -> Result<&'a str, NightlyError> {
    value
        .to_str()
        .ok_or_else(|| NightlyError::input(field, "must be valid UTF-8"))
}

fn require_byte_range(field: &'static str, value: &str, max: usize) -> Result<(), NightlyError> {
    if value.is_empty() || value.len() > max {
        return Err(NightlyError::input(
            field,
            format!(
                "UTF-8 encoding must contain 1..={max} bytes; observed {}",
                value.len()
            ),
        ));
    }
    Ok(())
}

fn require_json_bound(field: &'static str, value: &str, max: usize) -> Result<(), NightlyError> {
    if value.len() > max {
        return Err(NightlyError::input(
            field,
            format!(
                "encoded JSON exceeds the {max}-byte limit; observed {}",
                value.len()
            ),
        ));
    }
    Ok(())
}

fn provenance_json(
    field: &'static str,
    variable: &'static str,
    raw: Option<OsString>,
) -> Result<String, NightlyError> {
    match raw {
        Some(raw) => {
            let value = utf8(field, &raw)?;
            require_byte_range(field, value, MAX_PROVENANCE_VALUE_BYTES)?;
            Ok(format!(
                "{{\"availability\":\"reported\",\"source\":\"process_environment\",\
                 \"variable\":\"{variable}\",\"value\":\"{}\"}}",
                escape_json_string(value)
            ))
        }
        None => Ok(format!(
            "{{\"availability\":\"unavailable\",\"source\":\"process_environment\",\
             \"variable\":\"{variable}\"}}"
        )),
    }
}

#[derive(Debug)]
struct PreparedRequest {
    db: String,
    outcome_text: String,
    outcome: OpOutcome,
    suite: String,
    value: f64,
    versions: String,
    capability: String,
    ir: String,
    event_payload: String,
}

fn prepare_request(
    args: &[OsString],
    github_sha: Option<OsString>,
    runner_os: Option<OsString>,
) -> Result<PreparedRequest, NightlyError> {
    let [_, db, outcome_text, suite, value_text] = args else {
        return Err(NightlyError::input(
            "argv",
            "usage: nightly_ledger <db-path> <ok|error> <suite> <value>",
        ));
    };

    // Validation order is part of the agent-facing refusal contract. Keep all
    // pure admission complete before Ledger::open can create or alter a file.
    let db = utf8("db-path", db)?;
    require_byte_range("db-path", db, MAX_DB_PATH_BYTES)?;
    if db.contains('\0') {
        return Err(NightlyError::input("db-path", "must not contain NUL"));
    }

    let outcome_text = utf8("outcome", outcome_text)?;
    require_byte_range("outcome", outcome_text, MAX_OUTCOME_BYTES)?;
    let outcome = match outcome_text {
        "ok" => OpOutcome::Ok,
        "error" => OpOutcome::Error,
        _ => {
            return Err(NightlyError::input(
                "outcome",
                "must be exactly 'ok' or 'error'",
            ));
        }
    };

    let suite = utf8("suite", suite)?;
    require_byte_range("suite", suite, MAX_SUITE_BYTES)?;

    let value_text = utf8("value", value_text)?;
    require_byte_range("value", value_text, MAX_VALUE_TEXT_BYTES)?;
    let value = value_text
        .parse::<f64>()
        .map_err(|_| NightlyError::input("value", "must be a base-10 floating-point number"))?;
    if !value.is_finite() {
        return Err(NightlyError::input("value", "must be finite"));
    }

    let sha = provenance_json("GITHUB_SHA", "GITHUB_SHA", github_sha)?;
    let runner = provenance_json("RUNNER_OS", "RUNNER_OS", runner_os)?;
    let versions = format!("{{\"frankensim\":{sha}}}");
    let capability = format!("{{\"ops\":[\"ci.nightly\"],\"runner\":{runner}}}");
    let suite_json = escape_json_string(suite);
    let ir = format!("{{\"op\":\"ci.nightly\",\"suite\":\"{suite_json}\"}}");
    let event_payload =
        format!("{{\"kernel\":\"ci.nightly\",\"metric\":\"{suite_json}\",\"value\":{value}}}");

    require_json_bound("ir", &ir, MAX_OP_IR_BYTES)?;
    require_json_bound("versions", &versions, MAX_OP_VERSIONS_BYTES)?;
    require_json_bound("capability", &capability, MAX_OP_CAPABILITY_BYTES)?;
    require_json_bound("event-payload", &event_payload, MAX_EVENT_PAYLOAD_BYTES)?;

    Ok(PreparedRequest {
        db: db.to_owned(),
        outcome_text: outcome_text.to_owned(),
        outcome,
        suite: suite.to_owned(),
        value,
        versions,
        capability,
        ir,
        event_payload,
    })
}

fn rollback_after_error<T>(
    primary: LedgerError,
    rollback: impl FnOnce() -> Result<(), LedgerError>,
) -> Result<T, NightlyError> {
    match rollback() {
        Ok(()) => Err(NightlyError::Ledger(primary)),
        Err(rollback) => Err(NightlyError::TransactionCleanup { primary, rollback }),
    }
}

fn settle_transaction<T>(
    write_result: Result<T, LedgerError>,
    commit: impl FnOnce() -> Result<(), LedgerError>,
    rollback: impl FnOnce() -> Result<(), LedgerError>,
) -> Result<T, NightlyError> {
    match write_result {
        Ok(value) => match commit() {
            Ok(()) => Ok(value),
            Err(primary) => rollback_after_error(primary, rollback),
        },
        Err(primary) => rollback_after_error(primary, rollback),
    }
}

fn write_request(request: &PreparedRequest) -> Result<i64, NightlyError> {
    let ledger = Ledger::open(&request.db)?;
    ledger.begin()?;
    let explicits = FiveExplicits {
        seed: b"nightly",
        versions: &request.versions,
        budget: "{\"wall_s\":5400}",
        capability: &request.capability,
    };
    let write_result = (|| -> Result<i64, LedgerError> {
        let t0 = now_wall_ns();
        let op = ledger.begin_op(Some(b"ci-nightly"), &request.ir, &explicits, t0)?;
        ledger.record_metric(op, 0, &request.suite, request.value)?;
        ledger.append_event(&EventRow {
            session: Some(b"ci-nightly"),
            t: t0,
            kind: "benchmark_result",
            payload: Some(&request.event_payload),
        })?;
        ledger.finish_op(op, request.outcome, None, now_wall_ns())?;
        Ok(op)
    })();
    settle_transaction(write_result, || ledger.commit(), || ledger.rollback())
}

fn main() -> std::process::ExitCode {
    let args: Vec<OsString> = std::env::args_os().collect();
    let request = match prepare_request(
        &args,
        std::env::var_os("GITHUB_SHA"),
        std::env::var_os("RUNNER_OS"),
    ) {
        Ok(request) => request,
        Err(error) => return fail(&error),
    };

    match write_request(&request) {
        Ok(op) => {
            println!(
                "{{\"suite\":\"ci.nightly\",\"op\":{op},\"outcome\":\"{}\",\
                 \"metric\":\"{}\",\"value\":{},\"db\":\"{}\"}}",
                escape_json_string(&request.outcome_text),
                escape_json_string(&request.suite),
                request.value,
                escape_json_string(&request.db)
            );
            std::process::ExitCode::SUCCESS
        }
        Err(error) => fail(&error),
    }
}

#[cfg(test)]
mod tests {
    use std::cell::Cell;

    use super::*;

    fn args(db: &str, outcome: &str, suite: &str, value: &str) -> Vec<OsString> {
        ["nightly_ledger", db, outcome, suite, value]
            .into_iter()
            .map(OsString::from)
            .collect()
    }

    #[test]
    fn g0_json_string_escaping_is_exhaustive_for_ascii_controls() {
        const EXPECTED: [&str; 32] = [
            "\\u0000", "\\u0001", "\\u0002", "\\u0003", "\\u0004", "\\u0005", "\\u0006", "\\u0007",
            "\\b", "\\t", "\\n", "\\u000b", "\\f", "\\r", "\\u000e", "\\u000f", "\\u0010",
            "\\u0011", "\\u0012", "\\u0013", "\\u0014", "\\u0015", "\\u0016", "\\u0017", "\\u0018",
            "\\u0019", "\\u001a", "\\u001b", "\\u001c", "\\u001d", "\\u001e", "\\u001f",
        ];
        for (byte, expected) in (0_u8..=0x1f).zip(EXPECTED) {
            assert_eq!(escape_json_string(&char::from(byte).to_string()), expected);
        }
        assert_eq!(escape_json_string("\"\\/é🙂"), "\\\"\\\\/é🙂");

        let hostile = NightlyError::input("x\n\u{0001}", "bad\r\t\"\\\u{0002}");
        let json = error_json(&hostile);
        assert!(!json.chars().any(|ch| ch < ' '));
        assert!(json.contains("\\n\\u0001"));
        assert!(json.contains("\\r\\t\\\"\\\\\\u0002"));
    }

    #[test]
    fn g0_admission_order_and_typed_unavailable_provenance_are_exact() {
        let error = prepare_request(&args("", "bogus", "", "NaN"), None, None).unwrap_err();
        assert_eq!(
            error.to_string(),
            "NightlyInput: field 'db-path' rejected: UTF-8 encoding must contain 1..=4096 bytes; observed 0"
        );

        let error =
            prepare_request(&args("ledger.db", "bogus", "", "NaN"), None, None).unwrap_err();
        assert_eq!(
            error.to_string(),
            "NightlyInput: field 'outcome' rejected: must be exactly 'ok' or 'error'"
        );

        let error = prepare_request(
            &args("ledger.db", "ok", "suite", "NaN"),
            Some(OsString::new()),
            Some(OsString::new()),
        )
        .unwrap_err();
        assert_eq!(
            error.to_string(),
            "NightlyInput: field 'value' rejected: must be finite"
        );
        let error = prepare_request(
            &args("ledger.db", "ok", "suite", "1"),
            Some(OsString::new()),
            Some(OsString::new()),
        )
        .unwrap_err();
        assert_eq!(
            error.to_string(),
            "NightlyInput: field 'GITHUB_SHA' rejected: UTF-8 encoding must contain 1..=131072 bytes; observed 0"
        );
        let error = prepare_request(
            &args("ledger.db", "ok", "suite", "1"),
            Some(OsString::from("sha")),
            Some(OsString::new()),
        )
        .unwrap_err();
        assert_eq!(
            error.to_string(),
            "NightlyInput: field 'RUNNER_OS' rejected: UTF-8 encoding must contain 1..=131072 bytes; observed 0"
        );

        let prepared = prepare_request(&args("ledger.db", "ok", "suite", "1"), None, None)
            .expect("missing provenance is explicit, not fabricated");
        assert_eq!(
            prepared.versions,
            "{\"frankensim\":{\"availability\":\"unavailable\",\"source\":\"process_environment\",\"variable\":\"GITHUB_SHA\"}}"
        );
        assert_eq!(
            prepared.capability,
            "{\"ops\":[\"ci.nightly\"],\"runner\":{\"availability\":\"unavailable\",\"source\":\"process_environment\",\"variable\":\"RUNNER_OS\"}}"
        );
        assert!(!prepared.versions.contains("local"));
        assert!(!prepared.capability.contains("local"));

        let reported = prepare_request(
            &args("ledger.db", "error", "suite", "-0"),
            Some(OsString::from("sha\"\\\n")),
            Some(OsString::from("runner\r\t")),
        )
        .expect("valid reported provenance");
        assert!(reported.versions.contains(
            "{\"availability\":\"reported\",\"source\":\"process_environment\",\"variable\":\"GITHUB_SHA\",\"value\":\"sha\\\"\\\\\\n\"}"
        ));
        assert!(reported.capability.contains(
            "{\"availability\":\"reported\",\"source\":\"process_environment\",\"variable\":\"RUNNER_OS\",\"value\":\"runner\\r\\t\"}"
        ));
    }

    #[test]
    fn g0_every_public_input_cap_has_an_exact_boundary() {
        let exact_db = "d".repeat(MAX_DB_PATH_BYTES);
        prepare_request(&args(&exact_db, "ok", "s", "1"), None, None).expect("exact db-path cap");
        let over_db = "d".repeat(MAX_DB_PATH_BYTES + 1);
        assert!(
            prepare_request(&args(&over_db, "ok", "s", "1"), None, None)
                .unwrap_err()
                .to_string()
                .contains("observed 4097")
        );

        prepare_request(&args("d", "error", "s", "1"), None, None).expect("longest outcome token");
        assert_eq!(
            prepare_request(&args("d", "errors", "s", "1"), None, None)
                .unwrap_err()
                .to_string(),
            "NightlyInput: field 'outcome' rejected: UTF-8 encoding must contain 1..=5 bytes; observed 6"
        );

        let exact_suite = "s".repeat(MAX_SUITE_BYTES);
        prepare_request(&args("d", "ok", &exact_suite, "1"), None, None).expect("exact suite cap");
        let over_suite = "s".repeat(MAX_SUITE_BYTES + 1);
        assert!(
            prepare_request(&args("d", "ok", &over_suite, "1"), None, None)
                .unwrap_err()
                .to_string()
                .contains("observed 65537")
        );

        let exact_value = format!("{}1", "0".repeat(MAX_VALUE_TEXT_BYTES - 1));
        prepare_request(&args("d", "ok", "s", &exact_value), None, None)
            .expect("exact value-text cap");
        let over_value = format!("{}1", "0".repeat(MAX_VALUE_TEXT_BYTES));
        assert!(
            prepare_request(&args("d", "ok", "s", &over_value), None, None)
                .unwrap_err()
                .to_string()
                .contains("observed 257")
        );

        let exact_sha = "a".repeat(MAX_PROVENANCE_VALUE_BYTES);
        prepare_request(
            &args("d", "ok", "s", "1"),
            Some(OsString::from(exact_sha)),
            None,
        )
        .expect("exact provenance cap");
        let over_sha = "a".repeat(MAX_PROVENANCE_VALUE_BYTES + 1);
        assert!(
            prepare_request(
                &args("d", "ok", "s", "1"),
                Some(OsString::from(over_sha)),
                None,
            )
            .unwrap_err()
            .to_string()
            .contains("observed 131073")
        );

        let over_runner = "r".repeat(MAX_PROVENANCE_VALUE_BYTES + 1);
        assert!(
            prepare_request(
                &args("d", "ok", "s", "1"),
                Some(OsString::from("sha")),
                Some(OsString::from(over_runner)),
            )
            .unwrap_err()
            .to_string()
            .contains("field 'RUNNER_OS' rejected")
        );
    }

    #[test]
    fn g4_primary_and_rollback_failures_are_both_retained_in_order() {
        let primary = LedgerError::Invalid {
            field: "metric".to_owned(),
            problem: "primary failed".to_owned(),
        };
        let rollback = LedgerError::Sql {
            context: "rollback".to_owned(),
            detail: "cleanup failed".to_owned(),
        };
        let commit_called = Cell::new(false);
        let error = settle_transaction::<i64>(
            Err(primary.clone()),
            || {
                commit_called.set(true);
                Ok(())
            },
            || Err(rollback.clone()),
        )
        .unwrap_err();
        assert!(!commit_called.get());
        assert_eq!(
            error,
            NightlyError::TransactionCleanup {
                primary: primary.clone(),
                rollback: rollback.clone(),
            }
        );
        assert_eq!(
            error.to_string(),
            "NightlyTransactionCleanup: primary[LedgerInvalid]=LedgerInvalid: field 'metric' rejected: primary failed; rollback[LedgerSql]=LedgerSql: rollback: cleanup failed"
        );

        let preserved =
            settle_transaction::<i64>(Err(primary.clone()), || Ok(()), || Ok(())).unwrap_err();
        assert_eq!(preserved, NightlyError::Ledger(primary));

        let commit = LedgerError::Busy {
            context: "commit".to_owned(),
            detail: "write conflict".to_owned(),
        };
        let error = settle_transaction(Ok(7_i64), || Err(commit.clone()), || Err(rollback.clone()))
            .unwrap_err();
        assert_eq!(
            error,
            NightlyError::TransactionCleanup {
                primary: commit,
                rollback,
            }
        );
    }
}
