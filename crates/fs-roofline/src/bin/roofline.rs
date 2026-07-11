//! Roofline harness CLI (plan §14.4 nightly lane).
//!
//! Usage:
//!   roofline [--n <elements>] [--warmup <k>] [--reps <k>] [--ledger <db>]
//!            [--baseline <jsonl>] [--firmware <identity>]
//!
//! Probes the machine axes, runs the default kernel registry, prints one
//! JSON line per kernel (plus the axes line and the §14.1 coverage table),
//! and — when `--ledger` is given — records the run as ledger provenance
//! and reports staleness for every registered kernel.

use fs_roofline::kernels::production_registry_with_ledger;
use fs_roofline::{
    AxisBaselinePolicy, BaselineIdentity, BaselineStore, MachineAxes, SECTION_14_1_TARGETS,
    STALENESS_MAX_AGE_NS, days_since_epoch_now, finalize_registry_tuning, run_registry, staleness,
};
use std::io::Read;

fn json_escape(value: &str) -> String {
    use core::fmt::Write as _;

    let mut escaped = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            c if c.is_control() => {
                let _ = write!(escaped, "\\u{:04x}", u32::from(c));
            }
            c => escaped.push(c),
        }
    }
    escaped
}

fn fail(detail: &str) -> std::process::ExitCode {
    eprintln!(
        "{{\"error\":\"Roofline\",\"detail\":\"{}\"}}",
        json_escape(detail)
    );
    std::process::ExitCode::FAILURE
}

#[derive(Debug, PartialEq, Eq)]
struct CliArgs {
    n: usize,
    warmup: usize,
    reps: usize,
    ledger_path: Option<String>,
    baseline_path: Option<String>,
    firmware: Option<String>,
}

impl Default for CliArgs {
    fn default() -> Self {
        Self {
            n: 1 << 22, // 32 MiB per f64 buffer: streams past every L2/L3
            warmup: 2,
            reps: 9,
            ledger_path: None,
            baseline_path: None,
            firmware: None,
        }
    }
}

fn positive_usize(flag: &str, value: &str) -> Result<usize, String> {
    value
        .parse::<usize>()
        .ok()
        .filter(|value| *value > 0)
        .ok_or_else(|| format!("{flag} must be a positive integer"))
}

fn parse_args(args: &[String]) -> Result<CliArgs, String> {
    let mut parsed = CliArgs::default();
    let mut seen = std::collections::BTreeSet::new();
    let mut args = args.iter().skip(1);
    while let Some(flag) = args.next() {
        if !matches!(
            flag.as_str(),
            "--n" | "--warmup" | "--reps" | "--ledger" | "--baseline" | "--firmware"
        ) {
            return Err(format!("unknown roofline argument {flag:?}"));
        }
        if !seen.insert(flag.as_str()) {
            return Err(format!("duplicate roofline argument {flag:?}"));
        }
        let value = args
            .next()
            .filter(|value| !value.starts_with("--"))
            .ok_or_else(|| format!("{flag} requires a value"))?;
        if value.is_empty() {
            return Err(format!("{flag} requires a non-empty value"));
        }
        match flag.as_str() {
            "--n" => parsed.n = positive_usize(flag, value)?,
            "--warmup" => parsed.warmup = positive_usize(flag, value)?,
            "--reps" => parsed.reps = positive_usize(flag, value)?,
            "--ledger" => parsed.ledger_path = Some(value.clone()),
            "--baseline" => parsed.baseline_path = Some(value.clone()),
            "--firmware" => parsed.firmware = Some(value.clone()),
            _ => return Err(format!("unknown roofline argument {flag:?}")),
        }
    }
    if parsed.baseline_path.is_some() && parsed.firmware.is_none() {
        return Err("--firmware is required when --baseline is supplied".to_string());
    }
    Ok(parsed)
}

struct BaselineInputs {
    store: Option<BaselineStore>,
    identity: BaselineIdentity,
    now_day: u64,
}

impl BaselineInputs {
    fn policy(&self, fingerprint: u64) -> AxisBaselinePolicy<'_> {
        AxisBaselinePolicy::new(
            self.store
                .as_ref()
                .and_then(|store| store.for_fingerprint(fingerprint)),
            &self.identity,
            self.now_day,
        )
    }
}

fn parse_bounded_baseline_store(reader: impl Read, source: &str) -> Result<BaselineStore, String> {
    let limit = fs_roofline::baseline::MAX_BASELINE_STORE_BYTES;
    let bounded_bytes = limit
        .checked_add(1)
        .ok_or_else(|| "baseline-store read bound overflows usize".to_string())?;
    let read_limit = u64::try_from(bounded_bytes)
        .map_err(|_| "baseline-store read bound does not fit u64".to_string())?;
    let mut bytes = Vec::with_capacity(bounded_bytes);
    reader
        .take(read_limit)
        .read_to_end(&mut bytes)
        .map_err(|error| format!("cannot read baseline store {source:?}: {error}"))?;
    if bytes.len() > limit {
        return Err(format!(
            "baseline store {source:?} exceeds the {limit}-byte bound"
        ));
    }
    let text =
        String::from_utf8(bytes).map_err(|_| format!("baseline store {source:?} is not UTF-8"))?;
    BaselineStore::from_jsonl(&text).map_err(|error| error.to_string())
}

fn load_baseline_inputs(args: &CliArgs, axes: &MachineAxes) -> Result<BaselineInputs, String> {
    let identity = BaselineIdentity::current(
        axes,
        args.firmware.as_deref().unwrap_or("unbaselined-candidate"),
    )
    .map_err(|error| error.to_string())?;
    let now_day = days_since_epoch_now().map_err(|error| error.to_string())?;
    let store = match args.baseline_path.as_deref() {
        Some(path) => {
            let file = std::fs::File::open(path)
                .map_err(|error| format!("cannot read baseline store {path:?}: {error}"))?;
            Some(parse_bounded_baseline_store(file, path)?)
        }
        None => None,
    };
    Ok(BaselineInputs {
        store,
        identity,
        now_day,
    })
}

fn main() -> std::process::ExitCode {
    let raw_args: Vec<String> = std::env::args().collect();
    let args = match parse_args(&raw_args) {
        Ok(args) => args,
        Err(error) => return fail(&error),
    };

    let tune_ledger = match args.ledger_path.as_deref() {
        Some(path) => match fs_ledger::Ledger::open(path) {
            Ok(ledger) => Some(ledger),
            Err(error) => return fail(&error.to_string()),
        },
        None => None,
    };

    let axes = MachineAxes::probe();
    println!("{}", axes.to_jsonl());

    let baseline_inputs = match load_baseline_inputs(&args, &axes) {
        Ok(inputs) => inputs,
        Err(error) => return fail(&error),
    };
    let baseline_policy = baseline_inputs.policy(axes.fingerprint);

    let mut registry = production_registry_with_ledger(args.n, &axes, tune_ledger);
    let mut results = run_registry(&mut registry, args.warmup, args.reps, &axes);
    let post_axes = MachineAxes::probe();
    println!("{}", post_axes.to_jsonl());
    println!("{}", baseline_policy.receipt_json(&axes, &post_axes));
    let mut finalized =
        match finalize_registry_tuning(&mut registry, &axes, &post_axes, baseline_policy, &results)
        {
            Ok(finalized) => finalized,
            Err(error) => return fail(&error),
        };
    let citable = finalized.admitted();
    for r in &results {
        println!("{}", r.to_jsonl());
    }
    for row in SECTION_14_1_TARGETS {
        println!(
            "{{\"target\":\"{}\",\"statement\":\"{}\",\"landed\":{}}}",
            json_escape(row.kernel),
            json_escape(row.statement),
            row.landed
        );
    }

    // The registry owns fsqlite's deliberately !Send tune connection. Drop it
    // before reopening the same database for the atomic evidence transaction;
    // run_registry is synchronous, so no kernel work survives this point.
    drop(registry);
    if let Some(db) = args.ledger_path {
        let ledger = match fs_ledger::Ledger::open(&db) {
            Ok(l) => l,
            Err(e) => return fail(&e.to_string()),
        };
        match fs_roofline::record_run(
            &ledger,
            &axes,
            &post_axes,
            baseline_policy,
            &mut finalized,
            &mut results,
        ) {
            Ok(op) => {
                println!(
                    "{{\"ledgered\":true,\"citable\":{citable},\"op\":{op},\"db\":\"{}\"}}",
                    json_escape(&db)
                );
            }
            Err(e) => return fail(&e.to_string()),
        }
        for r in &results {
            match staleness(
                &ledger,
                &r.kernel,
                &r.version,
                axes.fingerprint,
                baseline_policy.baseline_hash(),
            ) {
                Ok(s) => println!(
                    "{{\"kernel\":\"{}\",\"staleness\":\"{s:?}\",\"max_age_ns\":{STALENESS_MAX_AGE_NS}}}",
                    json_escape(&r.kernel),
                ),
                Err(e) => return fail(&e.to_string()),
            }
        }
    }
    std::process::ExitCode::SUCCESS
}

#[cfg(test)]
mod tests {
    use super::{json_escape, load_baseline_inputs, parse_args, parse_bounded_baseline_store};
    use fs_roofline::MachineAxes;
    use std::io::Cursor;

    fn args(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| (*value).to_string()).collect()
    }

    #[test]
    fn manual_json_fields_escape_hostile_paths_and_diagnostics() {
        assert_eq!(
            json_escape("ledger\\\"row\n\t\u{0001}.db"),
            "ledger\\\\\\\"row\\n\\t\\u0001.db"
        );
    }

    #[test]
    fn baseline_store_requires_explicit_firmware_identity() {
        let axes = MachineAxes {
            fingerprint: 1,
            cpu_brand: "synthetic".to_string(),
            logical_cpus: 1,
            bandwidth_single_gbs: 10.0,
            bandwidth_all_core_gbs: 10.0,
            peak_single_gflops: 10.0,
            peak_all_core_gflops: 10.0,
        };
        let error = parse_args(&args(&["roofline", "--baseline", "x"]))
            .expect_err("firmware omission must fail before file access");
        assert!(error.contains("--firmware"));

        let parsed = parse_args(&args(&["roofline"])).expect("default invocation parses");
        let candidate =
            load_baseline_inputs(&parsed, &axes).expect("report-only invocation remains available");
        assert!(
            !candidate
                .policy(axes.fingerprint)
                .verdict(&axes, &axes)
                .trusted()
        );
    }

    #[test]
    fn parser_rejects_unknown_duplicate_missing_and_invalid_values() {
        for (case, expected) in [
            (vec!["roofline", "--unknown", "x"], "unknown"),
            (vec!["roofline", "--n", "1", "--n", "2"], "duplicate"),
            (vec!["roofline", "--ledger"], "requires a value"),
            (vec!["roofline", "--ledger", "--n", "1"], "requires a value"),
            (vec!["roofline", "--n", "0"], "positive integer"),
            (vec!["roofline", "--reps", "nope"], "positive integer"),
        ] {
            let error = parse_args(&args(&case)).expect_err("malformed argv must fail");
            assert!(
                error.contains(expected),
                "{error:?} did not contain {expected:?}"
            );
        }
    }

    #[test]
    fn parser_accepts_every_flag_once_and_preserves_report_only_default() {
        let defaults = parse_args(&args(&["roofline"])).expect("defaults");
        assert!(defaults.baseline_path.is_none());
        assert!(defaults.ledger_path.is_none());

        let parsed = parse_args(&args(&[
            "roofline",
            "--n",
            "8",
            "--warmup",
            "1",
            "--reps",
            "2",
            "--ledger",
            "run.db",
            "--baseline",
            "axes.jsonl",
            "--firmware",
            "os-build-1",
        ]))
        .expect("complete argv");
        assert_eq!(parsed.n, 8);
        assert_eq!(parsed.warmup, 1);
        assert_eq!(parsed.reps, 2);
        assert_eq!(parsed.ledger_path.as_deref(), Some("run.db"));
        assert_eq!(parsed.baseline_path.as_deref(), Some("axes.jsonl"));
        assert_eq!(parsed.firmware.as_deref(), Some("os-build-1"));
    }

    #[test]
    fn baseline_reader_stops_at_the_store_bound_plus_one_byte() {
        let oversized = vec![b'x'; fs_roofline::baseline::MAX_BASELINE_STORE_BYTES + 1];
        let error = parse_bounded_baseline_store(Cursor::new(oversized), "oversized.jsonl")
            .err()
            .expect("oversized input must fail before parsing");
        assert!(error.contains("exceeds"));
    }
}
