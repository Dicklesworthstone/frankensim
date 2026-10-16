//! E10.1 battery (bead wf-root-guzez.11.2): the harness EXECUTES the
//! pinned E4.9a batch (per-row class oracles, never totals-only) and
//! ingests the E4.9b dense reference; correction tables obey §4.2 —
//! calibration pins only, MANDATORY holdout evaluation carried in
//! the table, domain-bound, registry-bound. The DONE-WHEN hostile
//! twins are EXECUTED: stale application refuses, out-of-domain
//! application refuses (AT the bounds admits). Caps at cap AND
//! cap+1; determinism golden.
//! Repro: cargo test -p fs-flyer --test refereeharness_battery

use fs_flyer::referee::DiscrepancyReceipt;
use fs_flyer::refereeharness::{
    MAX_CORRECTION_KNOTS, PINNED_ALPHAS, apply_correction, build_correction_table, run_e49a_batch,
    run_harness,
};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-refereeharness\",\"case\":\"{case}\",{payload}}}");
}

#[test]
fn harness_executes_the_pinned_batch_with_per_row_oracles() {
    let receipt = run_harness().unwrap();
    assert_eq!(receipt.rows.len(), 5, "3 mono alphas + biplane + momentum");
    for row in &receipt.rows {
        println!("{}", row.to_jsonl());
        match row.comparison_class {
            "formulation-band-0.15" => assert!(
                row.rel_discrepancy.abs() < 0.15,
                "{}: {}",
                row.case_id,
                row.rel_discrepancy
            ),
            "below-ideal-within-0.5" => {
                assert!(
                    row.production < row.referee,
                    "{}: BEMT above ideal",
                    row.case_id
                );
                assert!(row.rel_discrepancy.abs() < 0.5, "{}", row.case_id);
            }
            "reported-only" => assert!(row.rel_discrepancy.abs() < 0.5, "{}", row.case_id),
            other => panic!("undeclared comparison class {other}"),
        }
    }
    // The E4.9b ingestion: the dense referee's Wagner-class ratio in
    // its CONTRACT band (shallow deficiency, ~0.9 class).
    assert!(
        receipt.wakeref_wagner_ratio > 0.8 && receipt.wakeref_wagner_ratio < 1.0,
        "wagner {}",
        receipt.wakeref_wagner_ratio
    );
    assert_eq!(receipt.wakeref_series_digest.len(), 64);
    // Determinism.
    let again = run_harness().unwrap();
    assert_eq!(
        again.receipt_digest, receipt.receipt_digest,
        "bit-identical twice"
    );
    jlog(
        "harness",
        &format!(
            "\"worst_abs_rel\":{},\"wagner\":{},\"digest\":\"{}\"",
            receipt.worst_abs_rel, receipt.wakeref_wagner_ratio, receipt.receipt_digest
        ),
    );
}

fn mono_pairs() -> Vec<(f64, DiscrepancyReceipt)> {
    let rows = run_e49a_batch().unwrap();
    PINNED_ALPHAS
        .iter()
        .copied()
        .zip(rows.into_iter().take(3))
        .collect()
}

const REG_ID: &str = "registry-v1-test-id";

#[test]
fn correction_table_calibrates_pins_and_evaluates_holdout() {
    let pairs = mono_pairs();
    // Calibrate on the outer pins, hold out the middle.
    let cal = vec![pairs[0].clone(), pairs[2].clone()];
    let holdout = vec![pairs[1].clone()];
    let raw_holdout_rel = pairs[1].1.rel_discrepancy.abs();
    let table = build_correction_table("wing_lift", &cal, &holdout, REG_ID).unwrap();
    assert_eq!(
        table.calibration_ids,
        vec!["e101-mono-a03", "e101-mono-a07"]
    );
    // The holdout error is IN the table (honest), and the correction
    // genuinely helps on the held-out pin.
    assert!(
        table.holdout_rel_worst < raw_holdout_rel,
        "correction must beat raw on holdout: {} vs {raw_holdout_rel}",
        table.holdout_rel_worst
    );
    // Applying at a calibration pin reproduces the referee exactly.
    let corrected = apply_correction(&table, pairs[0].0, pairs[0].1.production, REG_ID).unwrap();
    assert!(
        (corrected / pairs[0].1.referee - 1.0).abs() < 1e-12,
        "pin correction exact"
    );
    jlog(
        "correction",
        &format!(
            "\"holdout_rel\":{},\"raw_rel\":{raw_holdout_rel}",
            table.holdout_rel_worst
        ),
    );
}

#[test]
fn hostile_twins_stale_and_out_of_domain_refuse() {
    let pairs = mono_pairs();
    let table = build_correction_table(
        "wing_lift",
        &[pairs[0].clone(), pairs[2].clone()],
        &[pairs[1].clone()],
        REG_ID,
    )
    .unwrap();
    // STALE: applied under a different registry id.
    let stale = apply_correction(&table, 0.05, 1000.0, "registry-v2-other").unwrap_err();
    assert_eq!(stale.code, "correction-stale");
    // OUT-OF-DOMAIN: AT the bounds admits; a hair beyond refuses.
    assert!(apply_correction(&table, table.domain_lo, 1000.0, REG_ID).is_ok());
    assert!(apply_correction(&table, table.domain_hi, 1000.0, REG_ID).is_ok());
    let below = apply_correction(&table, table.domain_lo - 1e-12, 1000.0, REG_ID).unwrap_err();
    assert_eq!(below.code, "correction-out-of-domain");
    let above = apply_correction(&table, table.domain_hi + 1e-12, 1000.0, REG_ID).unwrap_err();
    assert_eq!(above.code, "correction-out-of-domain");
    jlog(
        "hostile-twins",
        "\"stale\":\"refused\",\"out_of_domain\":\"refused\"",
    );
}

#[test]
fn correction_caps_and_holdout_mandatory() {
    let pairs = mono_pairs();
    let mk = |x: f64, id_row: &DiscrepancyReceipt| (x, id_row.clone());
    // Holdout is MANDATORY: a table without one never exists.
    let no_holdout =
        build_correction_table("wing_lift", &[pairs[0].clone()], &[], REG_ID).unwrap_err();
    assert_eq!(no_holdout.code, "correction-holdout-missing");
    // Caps: 32 knots admits, 33 refuses; unordered refuses; empty refuses.
    let many: Vec<(f64, DiscrepancyReceipt)> = (0..MAX_CORRECTION_KNOTS)
        .map(|i| mk(0.01 + i as f64 * 2e-3, &pairs[0].1))
        .collect();
    assert!(build_correction_table("wing_lift", &many, &[pairs[1].clone()], REG_ID).is_ok());
    let over: Vec<(f64, DiscrepancyReceipt)> = (0..=MAX_CORRECTION_KNOTS)
        .map(|i| mk(0.01 + i as f64 * 2e-3, &pairs[0].1))
        .collect();
    assert_eq!(
        build_correction_table("wing_lift", &over, &[pairs[1].clone()], REG_ID)
            .unwrap_err()
            .code,
        "correction-calibration-invalid"
    );
    let unordered = vec![mk(0.07, &pairs[0].1), mk(0.03, &pairs[0].1)];
    assert_eq!(
        build_correction_table("wing_lift", &unordered, &[pairs[1].clone()], REG_ID)
            .unwrap_err()
            .code,
        "correction-calibration-invalid"
    );
    assert_eq!(
        build_correction_table("wing_lift", &[], &[pairs[1].clone()], REG_ID)
            .unwrap_err()
            .code,
        "correction-calibration-invalid"
    );
    jlog("caps", &format!("\"max_knots\":{MAX_CORRECTION_KNOTS}"));
}
