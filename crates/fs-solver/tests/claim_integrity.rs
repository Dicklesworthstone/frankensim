//! Claim-integrity regressions for fs-solver (E02 sweep, beads
//! `frankensim-extreal-program-f85xj.2.24` and `.2.25`).
//!
//! Both defects are about a solve that REPORTS more than it established:
//! a `converged` flag decided on a residual in the wrong norm behind a
//! guard that admitted the very state it names, and a V-cycle whose
//! coarse-solve receipts were dropped on the floor so its fixity
//! assumption had no evidence at all.

use fs_solver::{
    CgState, GmresState, MaskedTensorOp, PMultigrid, PminresState, ResidualClaim, SolveReport,
    StallDiagnosis, dot, norm2,
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
fn a_report_alone_names_its_residual_and_refuses_the_euclidean_reading() {
    // frankensim-extreal-program-f85xj.2.24, REMAINING half. Naming the
    // three quantities on the STATES was not enough: `SolveReport` is
    // what crosses the crate boundary, and a driver holding only a
    // report could not tell which quantity it had. Now the report
    // carries the claim, and the Euclidean reading is a REFUSAL when the
    // solve never established that number.
    //
    // Against the pre-fix shape this test does not compile: `SolveReport`
    // had no `residual_claim()`, no `euclidean_rel_residual()`, no
    // `converged_euclidean()`, and no constructor — it was an open struct
    // literal any crate could fill in without provenance.
    let a = Diag(vec![1.0, 1.0]);
    let b = [1.0, 1.0];

    // (1) P-MINRES: the M-norm estimate. `converged` is true, the number
    // is small, and the Euclidean reading REFUSES rather than hand the
    // M-norm number over under the Euclidean name.
    let m = DiagPrecond(vec![1.0, 1e-12]);
    let mut pminres = PminresState::new(&a, &m, &b);
    let pm_report = pminres.run(&a, &m, 1e-10, 50);
    assert!(
        matches!(pm_report.residual_claim(), ResidualClaim::MNormEstimate(_)),
        "the report itself must name the M-norm: {pm_report:?}"
    );
    assert_eq!(
        pm_report.euclidean_rel_residual(),
        None,
        "an M-norm estimate must not be readable as ‖b − Ax‖₂/‖b‖₂"
    );
    assert!(
        !pm_report.converged_euclidean(),
        "converged_euclidean() is the Euclidean reading and must fail closed \
         on an M-norm claim: {pm_report:?}"
    );
    assert_eq!(
        pm_report.rel_residual.to_bits(),
        pm_report.residual_claim().value().to_bits(),
        "the published magnitude IS the claim's magnitude — they cannot drift"
    );

    // (2) CG: the recursive estimate. Same refusal, different provenance.
    let mut cg = CgState::new(&a, &IdentityPrecond, &b);
    let cg_report = cg.run(&a, &IdentityPrecond, 1e-12, 50);
    assert!(
        cg_report.converged,
        "the CG fixture is meant to converge: {cg_report:?}"
    );
    assert!(
        matches!(
            cg_report.residual_claim(),
            ResidualClaim::RecursiveEstimate(_)
        ),
        "CG never recomputes b − Ax: {cg_report:?}"
    );
    assert_eq!(
        cg_report.euclidean_rel_residual(),
        None,
        "a recursively propagated estimate must not be readable as the true residual"
    );
    assert!(
        !cg_report.converged_euclidean(),
        "CG's `converged` is not a Euclidean correctness statement: {cg_report:?}"
    );
    assert_ne!(
        cg_report.residual_provenance(),
        pm_report.residual_provenance(),
        "the two estimates must not describe themselves identically"
    );

    // (3) GMRES: the one producer that recomputed it. The Euclidean
    // reading is granted, and it is the same number.
    let mut gmres = GmresState::new(&b, 2);
    let gm_report = gmres.run(&a, &b, 1e-12, 8, false);
    assert!(gm_report.converged, "GMRES fixture: {gm_report:?}");
    assert_eq!(
        gm_report.euclidean_rel_residual().map(f64::to_bits),
        Some(gm_report.rel_residual.to_bits()),
        "a recomputed true residual IS readable as the Euclidean residual"
    );
    assert!(
        gm_report.converged_euclidean(),
        "GMRES recomputes ‖b − Ax‖₂ at every cycle end: {gm_report:?}"
    );

    // (4) The gap the refusal exists for: the M-norm report says
    // `converged` while the TRUE Euclidean residual of its own iterate is
    // the number a driver would actually care about. The typed accessor
    // is what stops `pm_report.converged` from being read as that.
    let mut ax = vec![0.0; 2];
    fs_solver::LinearOp::apply(&a, &pminres.x, &mut ax);
    let residual: Vec<f64> = b.iter().zip(&ax).map(|(bi, ai)| bi - ai).collect();
    let true_euclidean = norm2(&residual) / norm2(&b);
    assert!(
        pm_report.converged && pm_report.euclidean_rel_residual().is_none(),
        "the report claims convergence in ITS measure and refuses the other one"
    );
    log(
        "report-carries-its-claim",
        &format!(
            "P-MINRES report: converged={} rel_residual={:.3e} ({}) euclidean=None; \
             its iterate's true Euclidean relative residual is {:.3e}",
            pm_report.converged,
            pm_report.rel_residual,
            pm_report.residual_provenance(),
            true_euclidean
        ),
    );
}

#[test]
fn every_report_is_built_from_a_claim_and_cannot_contradict_it() {
    // The constructor is the ONLY way in (the claim field is private and
    // the struct is #[non_exhaustive]), so `rel_residual`, `converged`,
    // and the claim are one decision rather than three fields a producer
    // could fill in inconsistently. Before the fix, fs-bem built a
    // `SolveReport` by struct literal with no provenance at all.
    let claim = ResidualClaim::MNormEstimate(1e-6);
    let report = SolveReport::from_claim(7, claim, 1e-5, vec![1.0, 1e-6]);
    assert_eq!(report.rel_residual.to_bits(), 1e-6_f64.to_bits());
    assert!(report.converged, "1e-6 < 1e-5 in the claim's own measure");
    assert!(
        !report.converged_euclidean(),
        "but that convergence is NOT a Euclidean statement"
    );
    assert_eq!(report.euclidean_rel_residual(), None);
    assert_eq!(report.iters, 7);
    assert!(report.diagnosis.is_none());

    // A diagnosis supplied by the producer is recorded only when the
    // claim did NOT meet tol: a success cannot be dressed with a failure
    // story, and a failure cannot be laundered into a success.
    let unresolved = SolveReport::from_claim_with_diagnosis(
        3,
        ResidualClaim::TrueEuclidean(0.5),
        1e-8,
        vec![1.0, 0.5],
        StallDiagnosis::Breakdown,
    );
    assert!(!unresolved.converged);
    assert_eq!(unresolved.diagnosis, Some(StallDiagnosis::Breakdown));
    let converged = SolveReport::from_claim_with_diagnosis(
        3,
        ResidualClaim::TrueEuclidean(1e-12),
        1e-8,
        vec![1.0, 1e-12],
        StallDiagnosis::Breakdown,
    );
    assert!(converged.converged);
    assert_eq!(
        converged.diagnosis, None,
        "a converged report carries no stall diagnosis, whatever the producer passes"
    );
    assert_eq!(
        converged.euclidean_rel_residual().map(f64::to_bits),
        Some(1e-12_f64.to_bits())
    );
    log(
        "report-constructor-is-total",
        "SolveReport::from_claim derives rel_residual/converged from the typed claim; \
         no producer in any crate can publish a magnitude without its provenance",
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
