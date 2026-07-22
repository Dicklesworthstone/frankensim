//! Second-order adjoints (adjoint-of-adjoint): the EXACT Hessian-vector
//! FORMULA through the implicit function theorem for the density-linear
//! problem class K(ρ)u = b with quadratic misfit J = ½‖u − u*‖². Exact
//! means no finite differencing and no dropped second-order term; the
//! four linear solves that evaluate the formula are iterative, so the
//! delivered accuracy is bounded by the tolerance they are run at (and
//! that tolerance is validated, not trusted).
//!
//! The general IFT Hessian needs four linear solves per Hv (tangent,
//! two second-derivative contractions, second adjoint). For this
//! class R is LINEAR in ρ (∂²R/∂ρ² = 0, ∂²R/∂ρ∂u·v = K(v)) and
//! ∂²J/∂u² = I, so the structure REDUCES to two extra solves per Hv —
//! stated, not hidden; the nonlinear-residual contractions land with
//! their first consumer:
//!
//!   1. tangent:        K·w  = −K(v)·u          (w = ∂u/∂ρ · v)
//!   2. second adjoint: K·λ̇ = w − K(v)·λ
//!   3. (Hv)_t = −(λ̇ᵀK_t u + λᵀK_t w)
//!
//! where K(v) = Σ_t v_t·K_t and λ is the FIRST-order adjoint at the
//! base point. Every solve shares the primal operator — the house
//! IFT discipline (never through Krylov iterations).

use crate::ift::{DensityOp, DensityPoisson};
use fs_solver::CgState;
use fs_sparse::precond::IdentityPrecond;

fn solve(op: &DensityOp<'_, '_>, b: &[f64], tol: f64) -> Vec<f64> {
    let mut st = CgState::new(op, &IdentityPrecond, b);
    let rep = st.run(op, &IdentityPrecond, tol, 50_000);
    assert!(rep.converged, "Hessian inner solve failed: {rep:?}");
    st.x
}

/// Misfit Hessian-vector product for the density-Poisson family: given
/// the base density `rho` (through `problem`), the right-hand side `b`,
/// the target `u_star`, and a density-space direction `v`, returns H·v
/// where H = ∇²_ρ J(ρ), J = ½‖u(ρ) − u*‖². Costs the primal + adjoint
/// (recomputed here for a self-contained call) + TWO extra solves.
///
/// EXACT refers to the FORMULA, not the arithmetic: the second-order
/// terms are the closed-form IFT contractions of a density-linear
/// operator (no finite differencing anywhere), but the four linear
/// solves that evaluate them are iterative, so the accuracy of the
/// returned vector is bounded by the residual each solve achieved at
/// `tol`.
///
/// # Panics
/// If `tol` is not a usable relative-residual tolerance, or an inner
/// solve fails to reach it. `tol` must satisfy `0 < tol < 1`: CG starts
/// at x = 0 with a relative residual of exactly 1, so any `tol ≥ 1`
/// exits before the first iteration and reports `converged = true` with
/// x = 0. All four solves would then return zeros and this function
/// would publish a silent-zero Hessian behind an assertion that passed
/// vacuously — a converged flag the caller's own tolerance can make
/// trivially true is not evidence.
#[must_use]
pub fn density_misfit_hvp(
    problem: &DensityPoisson<'_>,
    b: &[f64],
    u_star: &[f64],
    v: &[f64],
    tol: f64,
) -> Vec<f64> {
    assert!(
        tol.is_finite() && tol > 0.0 && tol < 1.0,
        "Hessian solve tolerance must satisfy 0 < tol < 1 (got {tol}): the initial relative \
         residual is 1, so a looser tolerance admits the zero iterate as converged"
    );
    let op = DensityOp::new(problem);
    // Base state and first-order adjoint.
    let u = solve(&op, b, tol);
    let misfit: Vec<f64> = u.iter().zip(u_star).map(|(a, s)| a - s).collect();
    let lambda = solve(&op, &misfit, tol); // K symmetric
    // K(v)·u and K(v)·λ via the perturbed-density apply: the operator
    // family is linear in ρ, so K(v) IS DensityOp at density v (with
    // zero floor disabled — v may be signed; apply_k handles signs).
    let kv_u = problem.apply_density(v, &u);
    let kv_lambda = problem.apply_density(v, &lambda);
    // Tangent solve: K·w = −K(v)·u.
    let neg_kvu: Vec<f64> = kv_u.iter().map(|x| -x).collect();
    let w = solve(&op, &neg_kvu, tol);
    // Second adjoint: K·λ̇ = w − K(v)·λ  (∂²J/∂u² = I ⇒ the w term).
    let rhs2: Vec<f64> = w.iter().zip(&kv_lambda).map(|(a, c)| a - c).collect();
    let lambda_dot = solve(&op, &rhs2, tol);
    // (Hv)_t = −(λ̇ᵀ K_t u + λᵀ K_t w): two mixed pullbacks.
    let p1 = problem.density_pullback(&lambda_dot, &u);
    let p2 = problem.density_pullback(&lambda, &w);
    p1.iter().zip(&p2).map(|(a, c)| -(a + c)).collect()
}
