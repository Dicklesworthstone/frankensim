//! Gradient-based Pareto TRACING: scalarization sweeps that produce
//! fronts of certificate-grade points — warm-started L-BFGS along a
//! weighted-sum schedule (convex fronts), and warm-started augmented-
//! Lagrangian ε-constraint continuation (the form that also covers
//! CONCAVE fronts, where weighted sums provably collapse to the
//! extremes — exhibited, not cited, in the battery). Every
//! ε-constraint point carries its KKT residual certificate.

use crate::auglag::{
    ConstrainedProblem, KktResidual, augmented_lagrangian, augmented_lagrangian_warm,
};
use crate::lbfgs::LbfgsState;
use crate::stop::StopRule;

/// One traced Pareto point.
#[derive(Debug, Clone)]
pub struct ParetoPoint {
    /// Decision vector.
    pub x: Vec<f64>,
    /// Objective values (f₁, f₂).
    pub f: [f64; 2],
    /// KKT certificate (ε-constraint path; None for weighted sums,
    /// whose certificate is the scalarized gradient norm).
    pub kkt: Option<KktResidual>,
    /// Scalarized gradient ∞-norm at the solution (weighted-sum path).
    pub grad_norm: f64,
}

/// Objective callback: x ↦ (f, ∇f). `Fn` (not FnMut) so the sweep can
/// wrap it for the constrained solver's split borrows.
pub type Objective<'a> = &'a dyn Fn(&[f64]) -> (f64, Vec<f64>);

fn assert_decision(x: &[f64]) {
    assert!(
        !x.is_empty() && x.iter().all(|v| v.is_finite()),
        "Pareto decision vectors must be non-empty and finite"
    );
}

fn eval_objective(label: &str, f: Objective<'_>, x: &[f64]) -> (f64, Vec<f64>) {
    assert_decision(x);
    let (value, grad) = f(x);
    assert!(value.is_finite(), "{label} objective value must be finite");
    assert_eq!(
        grad.len(),
        x.len(),
        "{label} gradient length must match the decision dimension"
    );
    assert!(
        grad.iter().all(|v| v.is_finite()),
        "{label} gradient entries must be finite"
    );
    (value, grad)
}

/// Weighted-sum sweep: for each w in `weights` (processed in order),
/// minimize w·f₁ + (1−w)·f₂ by L-BFGS, WARM-STARTED from the previous
/// solution (continuation along the front). Exact on convex fronts;
/// on concave fronts this collapses to extremes — use
/// [`epsilon_constraint_sweep`] there.
#[must_use]
pub fn weighted_sum_sweep(
    f1: Objective<'_>,
    f2: Objective<'_>,
    weights: &[f64],
    x0: &[f64],
) -> Vec<ParetoPoint> {
    assert_decision(x0);
    assert!(
        weights
            .iter()
            .all(|w| w.is_finite() && (0.0..=1.0).contains(w)),
        "Pareto weights must be finite and inside [0, 1]"
    );
    let mut x = x0.to_vec();
    let mut out = Vec::with_capacity(weights.len());
    for &w in weights {
        let mut fg = |xv: &[f64]| -> (f64, Vec<f64>) {
            let (a, ga) = eval_objective("f1", f1, xv);
            let (b, gb) = eval_objective("f2", f2, xv);
            let val = w.mul_add(a, (1.0 - w) * b);
            let g: Vec<f64> = ga
                .iter()
                .zip(&gb)
                .map(|(p, q)| w.mul_add(*p, (1.0 - w) * q))
                .collect();
            (val, g)
        };
        let mut st = LbfgsState::new(&x, 10, &mut fg);
        let rep = st.run(&mut fg, &StopRule::GradNorm(1e-10), 500);
        x = st.x.clone();
        let (a, _) = eval_objective("f1", f1, &x);
        let (b, _) = eval_objective("f2", f2, &x);
        out.push(ParetoPoint {
            x: x.clone(),
            f: [a, b],
            kkt: None,
            grad_norm: rep.grad_norm,
        });
    }
    out
}

/// ε-constraint sweep: for each ε in `epsilons` (in order), solve
/// min f₂ s.t. f₁ ≤ ε by the augmented Lagrangian, WARM-STARTED from
/// the previous solution. Covers concave fronts; every point carries
/// its KKT certificate.
#[must_use]
pub fn epsilon_constraint_sweep(
    f1: Objective<'_>,
    f2: Objective<'_>,
    epsilons: &[f64],
    x0: &[f64],
    tol: f64,
) -> Vec<ParetoPoint> {
    assert_decision(x0);
    assert!(
        epsilons.iter().all(|eps| eps.is_finite()),
        "Pareto epsilon constraints must be finite"
    );
    assert!(
        tol.is_finite() && tol > 0.0,
        "Pareto epsilon-constraint tolerance must be finite and positive"
    );
    let mut x = x0.to_vec();
    let mut out = Vec::with_capacity(epsilons.len());
    for &eps in epsilons {
        let mut fg = |xv: &[f64]| eval_objective("f2", f2, xv);
        let ci = |xv: &[f64]| -> Vec<f64> {
            let (a, _) = eval_objective("f1", f1, xv);
            vec![a - eps]
        };
        let ci_jt = |xv: &[f64], wv: &[f64]| -> Vec<f64> {
            assert_eq!(
                wv.len(),
                1,
                "epsilon-constraint multiplier dimension must be one"
            );
            let (_, ga) = eval_objective("f1", f1, xv);
            ga.iter().map(|g| g * wv[0]).collect()
        };
        let ce = |_: &[f64]| Vec::new();
        let ce_jt = |xv: &[f64], _: &[f64]| vec![0.0f64; xv.len()];
        let mut problem = ConstrainedProblem {
            fg: &mut fg,
            ce: &ce,
            ce_jt: &ce_jt,
            ci: &ci,
            ci_jt: &ci_jt,
        };
        let rep = augmented_lagrangian(&mut problem, &x, tol, 40);
        x.clone_from(&rep.x);
        let (a, _) = eval_objective("f1", f1, &x);
        let (b, _) = eval_objective("f2", f2, &x);
        out.push(ParetoPoint {
            x: x.clone(),
            f: [a, b],
            grad_norm: rep.kkt.stationarity,
            kkt: Some(rep.kkt),
        });
    }
    out
}

/// One traced tri-objective Pareto point (7tv.21.8).
#[derive(Debug, Clone)]
pub struct ParetoPoint3 {
    /// Decision vector.
    pub x: Vec<f64>,
    /// Objective values (f₁, f₂, f₃).
    pub f: [f64; 3],
    /// KKT certificate from the ε-constraint solve.
    pub kkt: KktResidual,
}

/// Tri-objective ε-constraint sweep (7tv.21.8): for each (ε₁, ε₂) pair
/// (in order), solve min f₃ s.t. f₁ ≤ ε₁, f₂ ≤ ε₂, and
/// lower ≤ x ≤ upper by the augmented Lagrangian, WARM-STARTED (both
/// iterate and multipliers) from the previous solution.
///
/// The box is a NATIVE part of the formulation, not a chart: standard
/// MOO fixtures (ZDT/DTLZ/WFG) are box-bounded, and smooth ℝ→(0,1)
/// charts necessarily have vanishing derivatives toward the box faces —
/// minimizing an objective whose unconstrained descent points at a face
/// (DTLZ2's f₃, measured) escapes into the flat asymptote where no
/// gradient signal can recover it. Explicit box rows keep full-strength
/// gradients everywhere. Every point carries its KKT certificate.
#[must_use]
#[expect(
    clippy::too_many_arguments,
    reason = "sweep = problem + schedule; grouping would obscure the call sites"
)]
pub fn epsilon_constraint_sweep3(
    f1: Objective<'_>,
    f2: Objective<'_>,
    f3: Objective<'_>,
    epsilon_pairs: &[(f64, f64)],
    lower: &[f64],
    upper: &[f64],
    x0: &[f64],
    tol: f64,
) -> Vec<ParetoPoint3> {
    assert_decision(x0);
    let n = x0.len();
    assert!(
        lower.len() == n && upper.len() == n,
        "box bounds must match the decision dimension"
    );
    assert!(
        lower
            .iter()
            .zip(upper)
            .zip(x0)
            .all(|((lo, hi), x)| lo.is_finite() && hi.is_finite() && lo < hi && lo <= x && x <= hi),
        "box bounds must be finite, ordered, and contain the start"
    );
    assert!(
        epsilon_pairs
            .iter()
            .all(|(e1, e2)| e1.is_finite() && e2.is_finite()),
        "Pareto epsilon constraints must be finite"
    );
    assert!(
        tol.is_finite() && tol > 0.0,
        "Pareto epsilon-constraint tolerance must be finite and positive"
    );
    let mut x = x0.to_vec();
    let mut out = Vec::with_capacity(epsilon_pairs.len());
    // Warm-started multiplier continuation: cold penalties activate only
    // after the iterate has already moved, which is too late when the
    // objective's descent path crosses a constraint (measured on the
    // DTLZ2 sphere octant).
    let mut nu_warm = vec![0.0f64; 2 + 2 * n];
    for &(e1, e2) in epsilon_pairs {
        let mut fg = |xv: &[f64]| eval_objective("f3", f3, xv);
        let ci = |xv: &[f64]| -> Vec<f64> {
            let (a, _) = eval_objective("f1", f1, xv);
            let (b, _) = eval_objective("f2", f2, xv);
            let mut rows = Vec::with_capacity(2 + 2 * xv.len());
            rows.push(a - e1);
            rows.push(b - e2);
            for ((lo, hi), xi) in lower.iter().zip(upper).zip(xv) {
                rows.push(lo - xi);
                rows.push(xi - hi);
            }
            rows
        };
        let ci_jt = |xv: &[f64], wv: &[f64]| -> Vec<f64> {
            assert_eq!(
                wv.len(),
                2 + 2 * xv.len(),
                "tri-objective epsilon-constraint multiplier dimension must be 2 + 2n"
            );
            let (_, ga) = eval_objective("f1", f1, xv);
            let (_, gb) = eval_objective("f2", f2, xv);
            let mut pull: Vec<f64> = ga
                .iter()
                .zip(&gb)
                .map(|(p, q)| p * wv[0] + q * wv[1])
                .collect();
            for (i, entry) in pull.iter_mut().enumerate() {
                // Box rows: ∇(lo−x_i) = −e_i, ∇(x_i−hi) = +e_i.
                *entry += wv[2 + 2 * i + 1] - wv[2 + 2 * i];
            }
            pull
        };
        let ce = |_: &[f64]| Vec::new();
        let ce_jt = |xv: &[f64], _: &[f64]| vec![0.0f64; xv.len()];
        let mut problem = ConstrainedProblem {
            fg: &mut fg,
            ce: &ce,
            ce_jt: &ce_jt,
            ci: &ci,
            ci_jt: &ci_jt,
        };
        let rep = augmented_lagrangian_warm(&mut problem, &x, tol, 40, &[], &nu_warm);
        x.clone_from(&rep.x);
        nu_warm.clone_from(&rep.nu);
        let (a, _) = eval_objective("f1", f1, &x);
        let (b, _) = eval_objective("f2", f2, &x);
        let (c, _) = eval_objective("f3", f3, &x);
        out.push(ParetoPoint3 {
            x: x.clone(),
            f: [a, b, c],
            kkt: rep.kkt,
        });
    }
    out
}
