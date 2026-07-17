//! Cross-ISA divergence report generator (bead 6nb.6): classify two
//! per-ISA artifact ledgers under a declared divergence policy and emit
//! the documentation-of-record markdown to stdout.
//!
//! Usage:
//!   cross-isa-report <isa-a-label> <ledger-a.jsonl> <isa-b-label> \
//!                    <ledger-b.jsonl> <policy.tsv>
//!
//! Ledger lines are the `emit_isa_ledger` test's output:
//!   {"detaudit_ledger":"<artifact>","hash":"<16hex>","value_bits":"<16hex>"|null}
//! Policy lines: `<artifact>\tfma` or `<artifact>\tulp:<max>`; artifacts
//! absent from the policy must match bit-for-bit.
//!
//! Exit status is non-zero when the report is not clean (any unclassified
//! row is a build failure).

use fs_detaudit::{DivergenceClass, DivergencePolicy, IsaLedger, LedgerRow, classify_cross_isa};
use std::collections::BTreeMap;

fn field<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    let tag = format!("\"{key}\":");
    let start = line.find(&tag)? + tag.len();
    let rest = &line[start..];
    if let Some(stripped) = rest.strip_prefix('"') {
        stripped.split('"').next()
    } else {
        rest.split([',', '}']).next()
    }
}

fn parse_ledger(isa: &str, text: &str) -> IsaLedger {
    let mut rows = BTreeMap::new();
    for line in text.lines() {
        let Some(artifact) = field(line, "detaudit_ledger") else {
            continue;
        };
        let Some(hash_hex) = field(line, "hash") else {
            continue;
        };
        let Ok(hash) = u64::from_str_radix(hash_hex, 16) else {
            continue;
        };
        let value_bits = match field(line, "value_bits") {
            Some("null") | None => None,
            Some(hex) => u64::from_str_radix(hex, 16).ok(),
        };
        rows.insert(artifact.to_owned(), LedgerRow { hash, value_bits });
    }
    IsaLedger {
        isa: isa.to_owned(),
        rows,
    }
}

fn parse_policy(text: &str) -> DivergencePolicy {
    let mut declared = BTreeMap::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let Some((artifact, class)) = line.split_once('\t') else {
            continue;
        };
        let class = if class == "fma" {
            DivergenceClass::FmaContraction
        } else if let Some(max) = class.strip_prefix("ulp:") {
            match max.parse::<u32>() {
                Ok(max_ulps) => DivergenceClass::LibmUlp { max_ulps },
                Err(_) => continue,
            }
        } else {
            continue;
        };
        declared.insert(artifact.to_owned(), class);
    }
    DivergencePolicy { declared }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 6 {
        eprintln!(
            "usage: cross-isa-report <isa-a> <ledger-a.jsonl> <isa-b> <ledger-b.jsonl> <policy.tsv>"
        );
        std::process::exit(2);
    }
    let read = |path: &str| -> String {
        std::fs::read_to_string(path).unwrap_or_else(|e| {
            eprintln!("cannot read {path}: {e}");
            std::process::exit(2);
        })
    };
    let a = parse_ledger(&args[1], &read(&args[2]));
    let b = parse_ledger(&args[3], &read(&args[4]));
    let policy = parse_policy(&read(&args[5]));
    let report = classify_cross_isa(&a, &b, &policy);
    print!("{}", report.render_markdown());
    if !report.clean() {
        std::process::exit(1);
    }
}
