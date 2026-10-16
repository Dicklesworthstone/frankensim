//! Bead 7tv.21.7: gradient-stack budget ledger — the machine-independent
//! optimizer-regression gate of `bbob_budget_ledger.rs`, extended to the
//! gradient-consuming engines.
//!
//! Engines do not all report evaluation counts, so nfev is
//! CLOSURE-COUNTED: the objective callback is ours, and every call is
//! tallied — engine-agnostic and honest. Each row asserts (1) its
//! success gate (the known optimum reached to documented precision),
//! (2) nfev at or under a pinned ceiling (measured value + ~30%
//! headroom, which also absorbs minor cross-ISA float drift in
//! iteration counts), and (3) a sanity floor. Rows are emitted as
//! fs-obs `BenchmarkResult` events (`machine: 0`) and wire-validated.

use core::cell::Cell;

use fs_ascent::auglag::ConstrainedProblem;
use fs_ascent::{LbfgsState, StopRule, augmented_lagrangian, interior_point, sqp};
use fsci_opt::{rosen, rosen_der};

#[path = "support/budget_trend.rs"]
mod budget_trend;

use budget_trend::{BUDGET_TREND_SCHEMA, GRADIENT_COMPONENT, gate_and_emit_budget_observation};

/// Branin (as in fsci_bo_cutest_oracle.rs); f* = 0.39788735772973816.
fn branin(x: &[f64]) -> (f64, Vec<f64>) {
    let (x1, x2) = (x[0], x[1]);
    let b = 5.1 / (4.0 * core::f64::consts::PI * core::f64::consts::PI);
    let c = 5.0 / core::f64::consts::PI;
    let t = 1.0 / (8.0 * core::f64::consts::PI);
    let inner = x2 - b * x1 * x1 + c * x1 - 6.0;
    (
        inner * inner + 10.0 * (1.0 - t) * x1.cos() + 10.0,
        vec![
            2.0 * inner * (-2.0 * b * x1 + c) - 10.0 * (1.0 - t) * x1.sin(),
            2.0 * inner,
        ],
    )
}

/// Hartmann-3 (as in fsci_bo_cutest_oracle.rs); optimum ≈ -3.86277978.
fn hartmann3(x: &[f64]) -> (f64, Vec<f64>) {
    const ALPHA: [f64; 4] = [1.0, 1.2, 3.0, 3.2];
    const A: [[f64; 3]; 4] = [
        [3.0, 10.0, 30.0],
        [0.1, 10.0, 35.0],
        [3.0, 10.0, 30.0],
        [0.1, 10.0, 35.0],
    ];
    const P: [[f64; 3]; 4] = [
        [0.3689, 0.1170, 0.2673],
        [0.4699, 0.4387, 0.7470],
        [0.1091, 0.8732, 0.5547],
        [0.0381, 0.5743, 0.8828],
    ];
    let mut f = 0.0;
    let mut g = vec![0.0; 3];
    for i in 0..4 {
        let mut inner = 0.0;
        for j in 0..3 {
            let d = x[j] - P[i][j];
            inner += A[i][j] * d * d;
        }
        let e = ALPHA[i] * (-inner).exp();
        f -= e;
        for (j, gj) in g.iter_mut().enumerate() {
            *gj += e * 2.0 * A[i][j] * (x[j] - P[i][j]);
        }
    }
    (f, g)
}

#[test]
fn lbfgs_budget_rows_hold_their_ceilings() {
    let mut em = fs_obs::Emitter::new(GRADIENT_COMPONENT, BUDGET_TREND_SCHEMA);
    // (kernel, fg, start, known optimum check). Canonical ceilings and
    // success gates live in the generated trend manifest.
    let fixtures: [(
        &str,
        &dyn Fn(&[f64]) -> (f64, Vec<f64>),
        Vec<f64>,
        &dyn Fn(f64) -> bool,
    ); 3] = [
        (
            "lbfgs/rosen4",
            &|x: &[f64]| (rosen(x), rosen_der(x)),
            vec![0.9, 0.9, 0.9, 0.9],
            &|f: f64| f < 1e-12,
        ),
        ("lbfgs/branin", &branin, vec![3.0, 2.0], &|f: f64| {
            (f - 0.397_887_357_729_738_16).abs() < 1e-8
        }),
        (
            "lbfgs/hartmann3",
            &hartmann3,
            vec![0.2, 0.5, 0.8],
            &|f: f64| (f - (-3.862_779_78)).abs() < 1e-5,
        ),
    ];
    for (kernel, fg, start, reached) in fixtures {
        let count = Cell::new(0usize);
        let mut counted = |x: &[f64]| {
            count.set(count.get() + 1);
            fg(x)
        };
        let mut st = LbfgsState::new(&start, 10, &mut counted);
        st.run(&mut counted, &StopRule::GradNorm(1e-10), 4000);
        let f = fg(&st.x).0;
        assert!(reached(f), "{kernel}: optimum not reached (f = {f:.6e})");
        gate_and_emit_budget_observation(&mut em, GRADIENT_COMPONENT, kernel, count.get(), 1, 1);
    }
}

#[test]
fn constrained_stack_budget_rows_hold_their_ceilings() {
    // Shared fixture from fsci_oracle.rs: min (x−2)² + (y−1)²
    // s.t. x + y = 2, x ≤ 1.2 — analytic optimum (1.2, 0.8).
    let mut em = fs_obs::Emitter::new(GRADIENT_COMPONENT, BUDGET_TREND_SCHEMA);
    let ce = |x: &[f64]| vec![x[0] + x[1] - 2.0];
    let ce_jt = |_: &[f64], w: &[f64]| vec![w[0], w[0]];
    let ci = |x: &[f64]| vec![x[0] - 1.2];
    let ci_jt = |_: &[f64], w: &[f64]| vec![w[0], 0.0];

    // Canonical ceilings and success gates live in the generated trend manifest.
    type Runner<'a> = &'a dyn Fn(&mut ConstrainedProblem<'_>) -> Vec<f64>;
    let runners: [(&str, Runner<'_>); 3] = [
        (
            "auglag/shared-constrained",
            &|p: &mut ConstrainedProblem<'_>| augmented_lagrangian(p, &[0.0, 0.0], 1e-9, 60).x,
        ),
        (
            "interior-point/shared-constrained",
            &|p: &mut ConstrainedProblem<'_>| interior_point(p, &[0.0, 0.0], 1e-8, 80).x,
        ),
        ("sqp/shared-constrained", &|p: &mut ConstrainedProblem<
            '_,
        >| {
            sqp(p, &[0.0, 0.0], 1e-9, 80).x
        }),
    ];
    for (kernel, run) in runners {
        let count = Cell::new(0usize);
        let mut fg = |x: &[f64]| {
            count.set(count.get() + 1);
            (
                (x[0] - 2.0).powi(2) + (x[1] - 1.0).powi(2),
                vec![2.0 * (x[0] - 2.0), 2.0 * (x[1] - 1.0)],
            )
        };
        let mut problem = ConstrainedProblem {
            fg: &mut fg,
            ce: &ce,
            ce_jt: &ce_jt,
            ci: &ci,
            ci_jt: &ci_jt,
        };
        let x = run(&mut problem);
        let dev = (x[0] - 1.2).abs().max((x[1] - 0.8).abs());
        assert!(
            dev < 1e-5,
            "{kernel}: analytic optimum (1.2, 0.8) missed by {dev:.3e} at {x:?}"
        );
        gate_and_emit_budget_observation(&mut em, GRADIENT_COMPONENT, kernel, count.get(), 1, 1);
    }
}
