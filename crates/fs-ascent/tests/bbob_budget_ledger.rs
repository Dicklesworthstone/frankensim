//! Bead 7tv.21.5: BBOB-class ERT-style budget ledger.
//!
//! nfev is MACHINE-INDEPENDENT given seeds and options, so these rows
//! CI-gate optimizer regressions with no perf lane: "an optimizer
//! regression is a perf regression". Stochastic families use the COCO
//! ERT estimator (total nfev across seeded runs ÷ number of successes),
//! so failed seeds honestly inflate the budget instead of vanishing;
//! deterministic DFO engines report single-run nfev.
//!
//! Every kernel asserts (1) its success gate (all seeds for the
//! separable/unimodal fixtures; ≥ 3/5 documented for DE on Rastrigin —
//! see the measured seed-1 finding below), (2) ERT/nfev at or under a
//! pinned ceiling (measured value + ~30% headroom, recorded per row),
//! and (3) an nfev sanity floor (vacuous-success guard). Rows are
//! emitted as fs-obs `BenchmarkResult` events (`machine: 0` =
//! machine-independent) and wire-validated — the ledger seed for trend
//! tracking.
//!
//! MEASURED FINDINGS kept on record:
//! - fsci `differential_evolution` under default options on Rastrigin-2
//!   over its CANONICAL domain [-5.12, 5.12]² converges to the
//!   neighboring local well (f ≈ 0.99496 at (0, −0.995)) on seed 1 and
//!   stays there at any maxiter; seeds 2–5 reach the origin (4/5).
//!   Shrinking the domain to [-5, 5]² drops that to 2/5 — the domain is
//!   part of the fixture identity, so the canonical bounds are pinned.
//! - fsci Powell on Rosenbrock-5 (`N = 5`): for `N >= 4`, extended
//!   Rosenbrock has a secondary strict local minimum near
//!   `(-0.96205, 0.93574, 0.88071, 0.77788, 0.60509)` with `f ≈ 3.930839`
//!   (Kok & Sandrock 2009); starting Powell with `x[1] = -1.5 < 0` sends
//!   the first coordinate line search into that negative-`x[0]` basin, so
//!   `powell/rosen5` uses the global-basin start `[2.0, 1.5, 1.0, 0.5, 2.5]`
//!   and reaches `f ≈ 1.08e-20 < 1e-8` in `nfev = 2505`.

use fsci_opt::{
    DifferentialEvolutionOptions, MinimizeOptions, OptimizeMethod, differential_evolution,
    minimize, rosen,
};

#[path = "support/budget_trend.rs"]
mod budget_trend;

use budget_trend::{BBOB_COMPONENT, BUDGET_TREND_SCHEMA, gate_and_emit_budget_observation};

fn sphere5(x: &[f64]) -> f64 {
    x.iter().map(|v| v * v).sum()
}

fn rastrigin2(x: &[f64]) -> f64 {
    20.0 + x
        .iter()
        .map(|&v| v * v - 10.0 * (2.0 * core::f64::consts::PI * v).cos())
        .sum::<f64>()
}

/// One ledger row: pinned ceiling + measured budget, emitted and gated.
struct Row {
    kernel: &'static str,
    /// ERT for stochastic kernels; single-run nfev for deterministic.
    budget: usize,
    successes: usize,
    attempts: usize,
}

fn ledger_and_gate(rows: &[Row]) {
    let mut em = fs_obs::Emitter::new(BBOB_COMPONENT, BUDGET_TREND_SCHEMA);
    for row in rows {
        gate_and_emit_budget_observation(
            &mut em,
            BBOB_COMPONENT,
            row.kernel,
            row.budget,
            row.successes,
            row.attempts,
        );
    }
}

#[test]
fn de_family_ert_rows_hold_their_ceilings() {
    // Seeded DE, DEFAULT options, canonical per-fixture domains, seeds
    // 1..=5. ERT = total nfev / successes (COCO estimator). Ceilings are
    // the measured ERT + ~30%: sphere5 6645 -> 8700; rosen5 21810 ->
    // 28500; rastrigin2 1898 (7590/4, seed-1 local trap on record) ->
    // 2500.
    let mut rows = Vec::new();
    for (kernel, f, n, bound, target) in [
        (
            "de/sphere5",
            &sphere5 as &(dyn Fn(&[f64]) -> f64),
            5usize,
            5.0f64,
            1e-8,
        ),
        ("de/rosen5", &(|x: &[f64]| rosen(x)), 5, 5.0, 1e-8),
        ("de/rastrigin2", &rastrigin2, 2, 5.12, 1e-6),
    ] {
        let bounds = vec![(-bound, bound); n];
        let mut total_nfev = 0usize;
        let mut successes = 0usize;
        for seed in 1u64..=5 {
            let res = differential_evolution(
                f,
                &bounds,
                DifferentialEvolutionOptions {
                    seed: Some(seed),
                    ..DifferentialEvolutionOptions::default()
                },
            )
            .expect("DE runs");
            total_nfev += res.nfev;
            if f(&res.x) < target {
                successes += 1;
            }
        }
        rows.push(Row {
            kernel,
            budget: total_nfev / successes.max(1),
            successes,
            attempts: 5,
        });
    }
    ledger_and_gate(&rows);
}

#[test]
fn dfo_family_budget_rows_hold_their_ceilings() {
    // Deterministic DFO engines from fixed starts; single-run nfev.
    // Success is judged on f(res.x) directly (some engines return
    // fun: None). Ceilings = measured nfev + ~30%, recorded per row.
    let start5 = [2.0f64, -1.5, 1.0, -0.5, 2.5];
    // For N >= 4, Rosenbrock has a secondary local minimum near
    // (-0.96205, 0.93574, 0.88071, 0.77788, 0.60509) with f ≈ 3.930839
    // (Kok & Sandrock 2009); starting Powell with x[1] = -1.5 < 0 sends the
    // first coordinate line search into that negative-x[0] basin.
    let start5_rosen_global = [2.0f64, 1.5, 1.0, 0.5, 2.5];
    let start2 = [0.4f64, -0.3]; // rastrigin: global-basin start (DFO is local)
    let mut rows = Vec::new();
    for (kernel, f, x0, target) in [
        (
            "nelder-mead/sphere5",
            &sphere5 as &(dyn Fn(&[f64]) -> f64),
            &start5[..],
            1e-8,
        ),
        (
            "nelder-mead/rosen5",
            &(|x: &[f64]| rosen(x)),
            &start5[..],
            1e-8,
        ),
        (
            "nelder-mead/rastrigin2-local",
            &rastrigin2,
            &start2[..],
            1e-6,
        ),
        ("powell/sphere5", &sphere5, &start5[..], 1e-8),
        ("powell/rosen5", &(|x: &[f64]| rosen(x)), &start5_rosen_global[..], 1e-8),
    ] {
        let method = if kernel.starts_with("powell") {
            OptimizeMethod::Powell
        } else {
            OptimizeMethod::NelderMead
        };
        let res = minimize(
            f,
            x0,
            MinimizeOptions {
                method: Some(method),
                tol: Some(1e-12),
                maxiter: Some(20_000),
                ..MinimizeOptions::default()
            },
        )
        .expect("DFO engine runs");
        let reached = f(&res.x) < target;
        rows.push(Row {
            kernel,
            budget: res.nfev,
            successes: usize::from(reached),
            attempts: 1,
        });
    }
    ledger_and_gate(&rows);
}
