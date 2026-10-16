//! Shared generated gate metadata for ASCENT optimizer-budget trend rows.

#![allow(dead_code)] // Each integration-test crate consumes a different shared subset.

use core::fmt::Write as _;
use std::collections::BTreeSet;

pub(super) const BUDGET_TREND_SCHEMA: &str = "frankensim-ascent-budget-trend-v1";
pub(super) const BBOB_COMPONENT: &str = "fs-ascent/bbob-budget";
pub(super) const GRADIENT_COMPONENT: &str = "fs-ascent/gradient-budget";
pub(super) const MACHINE_INDEPENDENT: u64 = 0;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct BudgetTrendRow {
    pub(super) suite: &'static str,
    pub(super) kernel: &'static str,
    pub(super) observed_metric: &'static str,
    pub(super) unit: &'static str,
    pub(super) ceiling: usize,
    pub(super) sanity_floor_exclusive: usize,
    pub(super) attempts: usize,
    pub(super) minimum_successes: usize,
}

/// Canonical row order is `(suite, kernel)` lexical order.
pub(super) const BUDGET_TREND_ROWS: [BudgetTrendRow; 14] = [
    BudgetTrendRow {
        suite: BBOB_COMPONENT,
        kernel: "de/rastrigin2",
        observed_metric: "ert_nfev",
        unit: "objective_evaluations",
        ceiling: 2_500,
        sanity_floor_exclusive: 10,
        attempts: 5,
        minimum_successes: 3,
    },
    BudgetTrendRow {
        suite: BBOB_COMPONENT,
        kernel: "de/rosen5",
        observed_metric: "ert_nfev",
        unit: "objective_evaluations",
        ceiling: 28_500,
        sanity_floor_exclusive: 10,
        attempts: 5,
        minimum_successes: 5,
    },
    BudgetTrendRow {
        suite: BBOB_COMPONENT,
        kernel: "de/sphere5",
        observed_metric: "ert_nfev",
        unit: "objective_evaluations",
        ceiling: 8_700,
        sanity_floor_exclusive: 10,
        attempts: 5,
        minimum_successes: 5,
    },
    BudgetTrendRow {
        suite: BBOB_COMPONENT,
        kernel: "nelder-mead/rastrigin2-local",
        observed_metric: "nfev_at_convergence",
        unit: "objective_evaluations",
        ceiling: 400,
        sanity_floor_exclusive: 10,
        attempts: 1,
        minimum_successes: 1,
    },
    BudgetTrendRow {
        suite: BBOB_COMPONENT,
        kernel: "nelder-mead/rosen5",
        observed_metric: "nfev_at_convergence",
        unit: "objective_evaluations",
        ceiling: 2_600,
        sanity_floor_exclusive: 10,
        attempts: 1,
        minimum_successes: 1,
    },
    BudgetTrendRow {
        suite: BBOB_COMPONENT,
        kernel: "nelder-mead/sphere5",
        observed_metric: "nfev_at_convergence",
        unit: "objective_evaluations",
        ceiling: 1_100,
        sanity_floor_exclusive: 10,
        attempts: 1,
        minimum_successes: 1,
    },
    BudgetTrendRow {
        suite: BBOB_COMPONENT,
        kernel: "powell/rosen5",
        observed_metric: "nfev_at_convergence",
        unit: "objective_evaluations",
        ceiling: 2_900,
        sanity_floor_exclusive: 10,
        attempts: 1,
        minimum_successes: 1,
    },
    BudgetTrendRow {
        suite: BBOB_COMPONENT,
        kernel: "powell/sphere5",
        observed_metric: "nfev_at_convergence",
        unit: "objective_evaluations",
        ceiling: 700,
        sanity_floor_exclusive: 10,
        attempts: 1,
        minimum_successes: 1,
    },
    BudgetTrendRow {
        suite: GRADIENT_COMPONENT,
        kernel: "auglag/shared-constrained",
        observed_metric: "nfev_at_convergence",
        unit: "objective_evaluations",
        ceiling: 450,
        sanity_floor_exclusive: 10,
        attempts: 1,
        minimum_successes: 1,
    },
    BudgetTrendRow {
        suite: GRADIENT_COMPONENT,
        kernel: "interior-point/shared-constrained",
        observed_metric: "nfev_at_convergence",
        unit: "objective_evaluations",
        ceiling: 4_850,
        sanity_floor_exclusive: 10,
        attempts: 1,
        minimum_successes: 1,
    },
    BudgetTrendRow {
        suite: GRADIENT_COMPONENT,
        kernel: "lbfgs/branin",
        observed_metric: "nfev_at_convergence",
        unit: "objective_evaluations",
        ceiling: 20,
        sanity_floor_exclusive: 10,
        attempts: 1,
        minimum_successes: 1,
    },
    BudgetTrendRow {
        suite: GRADIENT_COMPONENT,
        kernel: "lbfgs/hartmann3",
        observed_metric: "nfev_at_convergence",
        unit: "objective_evaluations",
        ceiling: 40,
        sanity_floor_exclusive: 10,
        attempts: 1,
        minimum_successes: 1,
    },
    BudgetTrendRow {
        suite: GRADIENT_COMPONENT,
        kernel: "lbfgs/rosen4",
        observed_metric: "nfev_at_convergence",
        unit: "objective_evaluations",
        ceiling: 45,
        sanity_floor_exclusive: 10,
        attempts: 1,
        minimum_successes: 1,
    },
    BudgetTrendRow {
        suite: GRADIENT_COMPONENT,
        kernel: "sqp/shared-constrained",
        observed_metric: "nfev_at_convergence",
        unit: "objective_evaluations",
        ceiling: 18,
        // A vacuous 0-step stop at x0 uses 3 objective evaluations (initial
        // validation + fg(x0) + kkt(x0)); any 1-step accepted iterate + KKT
        // verification requires >= 6 objective evaluations.
        sanity_floor_exclusive: 5,
        attempts: 1,
        minimum_successes: 1,
    },
];

#[must_use]
pub(super) fn audit_budget_trend_manifest(rows: &[BudgetTrendRow]) -> Vec<String> {
    let mut diagnostics = Vec::new();
    let mut keys = BTreeSet::new();
    for row in rows {
        let key = (row.suite, row.kernel);
        if !keys.insert(key) {
            diagnostics.push(format!(
                "duplicate budget trend row {}/{}",
                row.suite, row.kernel
            ));
        }
        if row.suite.trim().is_empty()
            || row.kernel.trim().is_empty()
            || row.observed_metric.trim().is_empty()
            || row.unit.trim().is_empty()
        {
            diagnostics.push(format!(
                "budget trend row {}/{} has a blank identity field",
                row.suite, row.kernel
            ));
        }
        if row.ceiling <= row.sanity_floor_exclusive {
            diagnostics.push(format!(
                "budget trend row {}/{} ceiling {} does not exceed sanity floor {}",
                row.suite, row.kernel, row.ceiling, row.sanity_floor_exclusive
            ));
        }
        if row.attempts == 0 || row.minimum_successes == 0 || row.minimum_successes > row.attempts {
            diagnostics.push(format!(
                "budget trend row {}/{} has invalid success gate {}/{}",
                row.suite, row.kernel, row.minimum_successes, row.attempts
            ));
        }
    }
    for pair in rows.windows(2) {
        let left = (pair[0].suite, pair[0].kernel);
        let right = (pair[1].suite, pair[1].kernel);
        if left >= right {
            diagnostics.push(format!(
                "budget trend rows are not in canonical (suite, kernel) order at {}/{} then {}/{}",
                pair[0].suite, pair[0].kernel, pair[1].suite, pair[1].kernel
            ));
        }
    }

    for expected in BUDGET_TREND_ROWS {
        match rows
            .iter()
            .find(|row| row.suite == expected.suite && row.kernel == expected.kernel)
        {
            Some(actual) if *actual != expected => diagnostics.push(format!(
                "budget trend row {}/{} metadata drifted from the canonical gate",
                expected.suite, expected.kernel
            )),
            None => diagnostics.push(format!(
                "canonical budget trend row {}/{} is missing",
                expected.suite, expected.kernel
            )),
            Some(_) => {}
        }
    }
    for row in rows {
        if !BUDGET_TREND_ROWS
            .iter()
            .any(|expected| expected.suite == row.suite && expected.kernel == row.kernel)
        {
            diagnostics.push(format!(
                "unexpected budget trend row {}/{}",
                row.suite, row.kernel
            ));
        }
    }
    diagnostics.sort();
    diagnostics.dedup();
    diagnostics
}

#[must_use]
pub(super) fn canonical_budget_trend_manifest_json() -> String {
    let mut out = format!(
        "{{\"schema\":\"{BUDGET_TREND_SCHEMA}\",\"authority\":\"regression-gate-declaration\",\"machine\":{MACHINE_INDEPENDENT},\"rows\":["
    );
    for (index, row) in BUDGET_TREND_ROWS.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        write!(
            out,
            "{{\"suite\":\"{}\",\"kernel\":\"{}\",\"observed_metric\":\"{}\",\"unit\":\"{}\",\"ceiling\":{},\"sanity_floor_exclusive\":{},\"attempts\":{},\"minimum_successes\":{}}}",
            row.suite,
            row.kernel,
            row.observed_metric,
            row.unit,
            row.ceiling,
            row.sanity_floor_exclusive,
            row.attempts,
            row.minimum_successes,
        )
        .expect("writing to a String is infallible");
    }
    out.push_str("]}");
    out
}

pub(super) fn gate_and_emit_budget_observation(
    emitter: &mut fs_obs::Emitter,
    suite: &str,
    kernel: &str,
    budget: usize,
    successes: usize,
    attempts: usize,
) {
    let row = BUDGET_TREND_ROWS
        .iter()
        .find(|row| row.suite == suite && row.kernel == kernel)
        .unwrap_or_else(|| panic!("{suite}/{kernel}: missing canonical budget trend row"));
    assert_eq!(
        attempts, row.attempts,
        "{suite}/{kernel}: attempt count drifted from the canonical fixture"
    );
    assert!(
        successes <= attempts,
        "{suite}/{kernel}: successes {successes} exceed attempts {attempts}"
    );
    assert!(
        successes >= row.minimum_successes,
        "{suite}/{kernel}: {successes}/{attempts} successes is below the canonical gate of {}",
        row.minimum_successes
    );
    assert!(
        budget > row.sanity_floor_exclusive,
        "{suite}/{kernel}: {budget} evaluations is vacuous at floor {}",
        row.sanity_floor_exclusive
    );
    assert!(
        budget <= row.ceiling,
        "{suite}/{kernel}: budget {budget} exceeds canonical ceiling {}",
        row.ceiling
    );

    let success_rate = successes as f64 / attempts as f64;
    let minimum_success_rate = row.minimum_successes as f64 / row.attempts as f64;
    for (metric, value) in [
        (row.observed_metric, budget as f64),
        ("budget_ceiling_nfev", row.ceiling as f64),
        ("success_rate", success_rate),
        ("minimum_success_rate", minimum_success_rate),
    ] {
        let event = emitter.emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::BenchmarkResult {
                kernel: kernel.to_string(),
                metric: metric.to_string(),
                value,
                machine: MACHINE_INDEPENDENT,
            },
            None,
        );
        let line = event.to_jsonl();
        fs_obs::validate_line(&line).expect("budget trend rows stay wire-valid");
        println!("{line}");
    }
}
