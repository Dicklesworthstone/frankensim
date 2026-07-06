//! Nightly CI ledger writer (foundations CI/CD bead; golden-ledger
//! convention, plan §16.2): records a nightly gate run as an op + metrics in
//! an fs-ledger database that CI retains as an artifact — the nightly run is
//! itself provenance, queryable like any other study.
//!
//! Usage:
//!   nightly_ledger <db-path> <outcome: ok|error> <suite> <value>
//!
//! The Five Explicits come from the CI environment (deterministic-friendly:
//! all inputs are explicit argv/env, no hidden state):
//!   seed        — not applicable to a gate run; the fixed literal b"nightly"
//!   versions    — {"frankensim": $GITHUB_SHA or "local", "constellation": lock hash}
//!   budget      — {"wall_s": 5400} (the job timeout)
//!   capability  — {"ops": ["ci.nightly"], "runner": $RUNNER_OS or "local"}
//!
//! Exit code 0 on success; 1 with a structured JSON error line on failure.

use fs_ledger::{EventRow, FiveExplicits, Ledger, OpOutcome, now_wall_ns};

fn fail(detail: &str) -> std::process::ExitCode {
    eprintln!("{{\"error\":\"NightlyLedger\",\"detail\":\"{detail}\"}}");
    std::process::ExitCode::FAILURE
}

fn main() -> std::process::ExitCode {
    let args: Vec<String> = std::env::args().collect();
    let [_, db, outcome_text, suite, value_text] = args.as_slice() else {
        return fail("usage: nightly_ledger <db-path> <ok|error> <suite> <value>");
    };
    let outcome = match outcome_text.as_str() {
        "ok" => OpOutcome::Ok,
        "error" => OpOutcome::Error,
        other => return fail(&format!("outcome must be ok|error, got {other}")),
    };
    let Ok(value) = value_text.parse::<f64>() else {
        return fail(&format!("value must be a float, got {value_text}"));
    };

    let sha = std::env::var("GITHUB_SHA").unwrap_or_else(|_| "local".to_string());
    let runner = std::env::var("RUNNER_OS").unwrap_or_else(|_| "local".to_string());
    let versions = format!("{{\"frankensim\":\"{sha}\"}}");
    let capability = format!("{{\"ops\":[\"ci.nightly\"],\"runner\":\"{runner}\"}}");
    let explicits = FiveExplicits {
        seed: b"nightly",
        versions: &versions,
        budget: "{\"wall_s\":5400}",
        capability: &capability,
    };

    let result = (|| -> Result<i64, fs_ledger::LedgerError> {
        let ledger = Ledger::open(db)?;
        let t0 = now_wall_ns();
        let ir = format!("{{\"op\":\"ci.nightly\",\"suite\":\"{suite}\"}}");
        let op = ledger.begin_op(Some(b"ci-nightly"), &ir, &explicits, t0)?;
        ledger.record_metric(op, 0, suite, value)?;
        ledger.append_event(&EventRow {
            session: Some(b"ci-nightly"),
            t: t0,
            kind: "benchmark_result",
            payload: Some(&format!(
                "{{\"kernel\":\"ci.nightly\",\"metric\":\"{suite}\",\"value\":{value}}}"
            )),
        })?;
        ledger.finish_op(op, outcome, None, now_wall_ns())?;
        Ok(op)
    })();

    match result {
        Ok(op) => {
            println!(
                "{{\"suite\":\"ci.nightly\",\"op\":{op},\"outcome\":\"{outcome_text}\",\
                 \"metric\":\"{suite}\",\"value\":{value},\"db\":\"{db}\"}}"
            );
            std::process::ExitCode::SUCCESS
        }
        Err(e) => fail(&e.to_string().replace('"', "'")),
    }
}
