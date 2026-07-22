//! Claim-integrity regressions for fs-solver (E02 sweep, beads
//! `frankensim-extreal-program-f85xj.2.24` and `.2.25`).
//!
//! Both defects are about a solve that REPORTS more than it established:
//! a `converged` flag decided on a residual in the wrong norm behind a
//! guard that admitted the very state it names, and a V-cycle whose
//! coarse-solve receipts were dropped on the floor so its fixity
//! assumption had no evidence at all.

use fs_solver::{
    CgState, GmresState, MaskedTensorOp, PMultigrid, PminresState, ResidualClaim, dot, norm2,
};
use fs_sparse::precond::{IdentityPrecond, Precond};

fn log(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-solver/claim-integrity\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

/// A DIAGONAL operator (symmetric, so P-MINRES accepts it).
struct Diag(Vec<f64>);

impl fs_solver::LinearOp for Diag {
    fn n(&self) -> usize {
        self.0.len()
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        for (i, yi) in y.iter_mut().enumerate() {
            *yi = self.0[i] * x[i];
        }
    }
}

/// A diagonal preconditioner. `diag = [1.0, 0.0]` is positive
/// SEMIdefinite — legal to construct, and exactly the input the old
/// positivity guard admitted.
struct DiagPrecond(Vec<f64>);

impl Precond for DiagPrecond {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        for (i, zi) in z.iter_mut().enumerate() {
            *zi = self.0[i] * r[i];
        }
    }
}

#[test]
#[should_panic(expected = "preconditioner lost positivity")]
fn pminres_refuses_a_semidefinite_preconditioner() {
    // frankensim-extreal-program-f85xj.2.24(b). The guard was
    // `assert!(vz >= -1e-30)`, which ADMITS vz == 0 — precisely where
    // positivity is lost. beta_next then collapsed to the smallest
    // denormal (2.2e-308), s_k to ~0, eta to ~0, and the next tolerance
    // check reported converged:true with rel_residual ~1e-308 for an
    // ARBITRARY iterate. A guard whose message is "preconditioner lost
    // positivity" must not pass exactly there.
    let a = Diag(vec![1.0, 2.0]);
    let m = DiagPrecond(vec![1.0, 0.0]);
    let b = [1.0, 1.0];
    let mut st = PminresState::new(&a, &m, &b);
    let _ = st.run(&a, &m, 1e-10, 50);
}

#[test]
fn pminres_happy_breakdown_still_solves_exactly() {
    // The other state that shares `⟨p, Mp⟩ ≈ 0`: p == 0 is a Lanczos
    // HAPPY breakdown, where the final Givens step with beta = 0 makes
    // the iterate exact. Refusing it would be over-refusal, so the two
    // are discriminated on ‖p‖ rather than conflated.
    let a = Diag(vec![1.0, 1.0]);
    let m = DiagPrecond(vec![1.0, 1e-12]);
    let b = [1.0, 0.0];
    let mut st = PminresState::new(&a, &m, &b);
    let report = st.run(&a, &m, 1e-10, 50);
    assert!(report.converged, "{report:?}");
    let mut ax = vec![0.0; 2];
    fs_solver::LinearOp::apply(&a, &st.x, &mut ax);
    let residual: Vec<f64> = b.iter().zip(&ax).map(|(bi, ai)| bi - ai).collect();
    assert!(
        norm2(&residual) < 1e-14,
        "the happy-breakdown iterate must be the exact solution: x = {:?}",
        st.x
    );
    log(
        "pminres-happy-breakdown",
        "p == 0 terminates with the exact iterate; only a NONZERO direction with a \
         non-positive M-inner product is refused",
    );
}

#[test]
fn residual_claims_name_the_norm_each_solver_reports() {
    // frankensim-extreal-program-f85xj.2.24(a). `SolveReport.rel_residual`
    // is documented as "‖r‖/‖b‖" but carries three different quantities
    // by producer. The quantities are now NAMED, so a driver can stop
    // reading `converged` as a Euclidean statement.
    let a = Diag(vec![2.0, 3.0]);
    let b = [1.0, 1.0];

    let cg = CgState::new(&a, &IdentityPrecond, &b);
    assert!(
        matches!(cg.residual_claim(), ResidualClaim::RecursiveEstimate(_)),
        "CG never recomputes b - Ax: {:?}",
        cg.residual_claim()
    );
    assert!(!cg.residual_claim().is_true_euclidean());

    let gmres = GmresState::new(&b, 4);
    assert!(
        gmres.residual_claim().is_true_euclidean(),
        "GMRES recomputes the true residual at every cycle end: {:?}",
        gmres.residual_claim()
    );

    let m = DiagPrecond(vec![1.0, 1e-12]);
    let pminres = PminresState::new(&a, &m, &b);
    let claim = pminres.residual_claim();
    assert!(
        matches!(claim, ResidualClaim::MNormEstimate(_)),
        "P-MINRES reports the M-norm, not the Euclidean norm: {claim:?}"
    );
    assert_eq!(
        claim.value().to_bits(),
        pminres.rel_residual().to_bits(),
        "the claim carries the number the report publishes"
    );

    // The bead's numeric witness for WHY the distinction is load-bearing:
    // with the SPD (hence legal) preconditioner M = diag(1, 1e-12), a
    // residual r = (0,1) against b = (1,0) is a 100% Euclidean relative
    // residual and a 1e-6 M-norm one. At tol = 1e-5 the M-norm number
    // reports `converged` for a completely unsolved system.
    let m_diag = [1.0f64, 1e-12];
    let residual = [0.0f64, 1.0];
    let rhs = [1.0f64, 0.0];
    let m_norm = |v: &[f64]| -> f64 {
        let scaled: Vec<f64> = v.iter().zip(&m_diag).map(|(vi, mi)| mi * vi).collect();
        fs_math::det::sqrt(dot(v, &scaled))
    };
    let m_relative = m_norm(&residual) / m_norm(&rhs);
    let euclidean_relative = norm2(&residual) / norm2(&rhs);
    assert!(
        (m_relative - 1e-6).abs() < 1e-18,
        "M-norm relative residual: {m_relative}"
    );
    assert!(
        (euclidean_relative - 1.0).abs() < 1e-18,
        "Euclidean relative residual: {euclidean_relative}"
    );
    assert!(
        m_relative < 1e-5 && euclidean_relative > 1e-5,
        "the two norms straddle a plausible tolerance: {m_relative} vs {euclidean_relative}"
    );
    log(
        "residual-claim-provenance",
        "CG/MINRES report a recursive estimate, P-MINRES an M-norm estimate, GMRES/FGMRES \
         the true Euclidean residual; the M-norm gap reaches 1e-6 vs 1.0 on a legal SPD \
         preconditioner",
    );
}

#[test]
fn pmg_retains_its_coarse_solve_receipts() {
    // frankensim-extreal-program-f85xj.2.25. Both `pcg(...)` call sites
    // discarded their PcgReport with `let _ =`, so if the r = 1 coarse
    // solve ever exited on its 2000-iteration cap instead of at 1e-13
    // the V-cycle became an inexact, application-dependent (hence
    // VARYING) preconditioner — and no surface could say it happened.
    // The CONTRACT invariant "the V-cycle preconditioner is symmetric
    // (… near-exact coarse)" had literally no evidence behind it.
    let op = MaskedTensorOp::new(3, 2);
    let mut b = vec![0.0f64; op.space().ndof()];
    for (i, (bi, &mk)) in b.iter_mut().zip(op.mask()).enumerate() {
        if mk {
            #[allow(clippy::cast_precision_loss)]
            let value = 1.0 + (i % 7) as f64 / 7.0;
            *bi = value;
        }
    }
    let pmg = PMultigrid::new(3, 2, 3);
    // Construction already measures lambda_max through the smoother, so
    // scope the evidence to the solve under test.
    pmg.reset_coarse_evidence();
    assert_eq!(pmg.coarse_solves(), 0);
    assert!(pmg.coarse_solves_converged());

    let mut st = CgState::new(&op, &pmg, &b);
    let report = st.run(&op, &pmg, 1e-10, 100);
    assert!(report.converged, "pMG-CG failed: {report:?}");

    assert!(
        pmg.coarse_solves() > 0,
        "the V-cycle must record the coarse solves it performed"
    );
    assert!(
        pmg.coarse_solves_converged(),
        "every coarse solve must have met its 1e-13 request; worst = {}",
        pmg.worst_coarse_rel_residual()
    );
    assert!(
        pmg.worst_coarse_rel_residual() < 1e-13,
        "the retained worst coarse residual is the evidence for V-cycle fixity: {}",
        pmg.worst_coarse_rel_residual()
    );
    log(
        "pmg-coarse-receipts",
        &format!(
            "{} coarse solves retained, worst relative residual {:.3e}",
            pmg.coarse_solves(),
            pmg.worst_coarse_rel_residual()
        ),
    );
}
