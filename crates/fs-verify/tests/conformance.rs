//! fs-verify conformance (bead lmp4.1, feature `certified-speculation`).
//! The MMS upper-bound property over the battery INCLUDING adversarial
//! perturbed candidates (the untrusted-proposer case), effectivity
//! bands, interval soundness + fail-closed, G5 determinism, the
//! certify-the-certifiers injection, the estimator-family falsifier,
//! and the nonlinear warm-start fallback with ledger rows. JSON-line
//! verdicts; seeded cases carry seeds.

use fs_verify::estimator::{
    EstimatorFamily, VerifierCheckpointKind, VerifierPhase, VerifierProgress, VerifierRefusal,
    VerifierWorkPlan, effectivity as try_effectivity,
    hierarchical_estimate as try_hierarchical_estimate, verify, verify_with_checkpoint,
    warm_start as try_warm_start,
};
use fs_verify::fem1d::{
    Fem1dError, MAX_FEM1D_MESH_NODES, MAX_FEM1D_POLY_COEFFICIENTS, MmsProblem, Poly,
    solve_p1 as try_solve_p1, true_energy_error as try_true_energy_error,
};

fn poly(coefficients: Vec<f64>) -> Poly {
    Poly::new(coefficients).expect("valid conformance polynomial")
}

fn problem(name: &str, u: Poly, mesh: Vec<f64>) -> MmsProblem {
    MmsProblem::new(name, u, mesh).expect("valid conformance problem")
}

fn solve_p1(problem: &MmsProblem) -> Vec<f64> {
    try_solve_p1(problem).expect("conformance problem must solve")
}

fn true_energy_error(problem: &MmsProblem, candidate: &[f64]) -> f64 {
    try_true_energy_error(problem, candidate).expect("conformance oracle must evaluate")
}

fn effectivity(
    problem: &MmsProblem,
    candidate: &[f64],
    report: &fs_verify::estimator::VerifierReport,
) -> f64 {
    try_effectivity(problem, candidate, report).expect("conformance effectivity must evaluate")
}

fn hierarchical_estimate(problem: &MmsProblem, candidate: &[f64]) -> f64 {
    try_hierarchical_estimate(problem, candidate).expect("conformance hierarchy must evaluate")
}

fn warm_start(
    problem: &MmsProblem,
    candidate: &[f64],
    max_iter: u32,
) -> fs_verify::estimator::WarmStartReport {
    try_warm_start(problem, candidate, max_iter).expect("conformance warm start must converge")
}

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-verify/conformance\",\"case\":\"{case}\",\"verdict\":\"{}\",\
         \"detail\":\"{detail}\"}}",
        if pass { "pass" } else { "fail" }
    );
    assert!(pass, "case {case}: {detail}");
}

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0
    }

    fn unit(&mut self) -> f64 {
        ((self.next() >> 11) as f64) / (1u64 << 53) as f64
    }
}

/// Polynomials vanishing at 0 and 1, degree ≤ 5 (keeps the squared
/// integrand within 5-point Gauss exactness).
fn mms_zoo() -> Vec<(&'static str, Poly)> {
    // x(1−x) = x − x²
    let u1 = poly(vec![0.0, 1.0, -1.0]);
    // x(1−x)(x−0.25) = −0.25x + 1.25x² − x³ (dyadic-exact).
    let u2 = poly(vec![0.0, -0.25, 1.25, -1.0]);
    // x²(1−x)² = x² − 2x³ + x⁴
    let u3 = poly(vec![0.0, 0.0, 1.0, -2.0, 1.0]);
    // x(1−x)(x−0.25)(x−0.75), expanded in dyadic-exact coefficients.
    let u4 = poly(vec![0.0, 0.1875, -1.1875, 2.0, -1.0]);
    // Degree 5: x(1−x)(x−0.5)(x²+0.25), also dyadic-exact.
    let u5 = poly(vec![0.0, -0.125, 0.375, -0.75, 1.5, -1.0]);
    vec![("u1", u1), ("u2", u2), ("u3", u3), ("u4", u4), ("u5", u5)]
}

fn meshes() -> Vec<Vec<f64>> {
    let uniform = |n: usize| -> Vec<f64> { (0..=n).map(|i| i as f64 / n as f64).collect() };
    let mut graded = vec![0.0];
    let mut x = 0.0;
    let mut h = 0.02;
    while x + h < 1.0 {
        x += h;
        graded.push(x);
        h *= 1.4;
    }
    graded.push(1.0);
    vec![
        uniform(4),
        uniform(8),
        uniform(16),
        uniform(64),
        graded,
        uniform(2),
    ]
}

fn assert_same_verifier_report(
    left: &fs_verify::estimator::VerifierReport,
    right: &fs_verify::estimator::VerifierReport,
) {
    assert_eq!(left.bound.lo.to_bits(), right.bound.lo.to_bits());
    assert_eq!(left.bound.hi.to_bits(), right.bound.hi.to_bits());
    assert_eq!(left.accept, right.accept);
    assert_eq!(left.color, right.color);
    assert_eq!(left.tolerance.to_bits(), right.tolerance.to_bits());
    assert_eq!(left.family, right.family);
    assert_eq!(left.flux_hash, right.flux_hash);
    assert_eq!(left.refusal, right.refusal);
}

#[test]
fn g4_verifier_work_plan_sparse_trace_and_legacy_equivalence_are_exact() {
    let p = problem(
        "checkpoint-small",
        poly(vec![0.0, 1.0, -1.0]),
        vec![0.0, 0.5, 1.0],
    );
    let candidate = [0.0, 0.0, 0.0];
    let plan = VerifierWorkPlan::for_inputs(&p, &candidate).expect("bounded verifier shape");
    assert_eq!(plan.identity_fields(), [27, 2, 2, 16, 1, 48]);

    let mut trace = Vec::new();
    let explicit = verify_with_checkpoint(&p, &candidate, 10.0, |progress| {
        trace.push(progress);
        Ok::<(), core::convert::Infallible>(())
    })
    .expect("infallible callback");
    let legacy = verify(&p, &candidate, 10.0);
    assert_same_verifier_report(&explicit, &legacy);
    assert_eq!(
        trace,
        vec![
            VerifierProgress {
                kind: VerifierCheckpointKind::PhaseEntry,
                phase: VerifierPhase::Validation,
                completed_work_units: 0,
                planned_work_units: 48,
            },
            VerifierProgress {
                kind: VerifierCheckpointKind::PhaseEntry,
                phase: VerifierPhase::Tightness,
                completed_work_units: 27,
                planned_work_units: 48,
            },
            VerifierProgress {
                kind: VerifierCheckpointKind::PhaseEntry,
                phase: VerifierPhase::Equilibrated,
                completed_work_units: 29,
                planned_work_units: 48,
            },
            VerifierProgress {
                kind: VerifierCheckpointKind::PhaseEntry,
                phase: VerifierPhase::Hash,
                completed_work_units: 31,
                planned_work_units: 48,
            },
            VerifierProgress {
                kind: VerifierCheckpointKind::PhaseEntry,
                phase: VerifierPhase::Finalization,
                completed_work_units: 47,
                planned_work_units: 48,
            },
            VerifierProgress {
                kind: VerifierCheckpointKind::Publication,
                phase: VerifierPhase::Finalization,
                completed_work_units: 48,
                planned_work_units: 48,
            },
        ]
    );

    for target in trace {
        let result = verify_with_checkpoint(&p, &candidate, 10.0, |progress| {
            if progress == target {
                Err(target)
            } else {
                Ok(())
            }
        });
        assert!(matches!(result, Err(error) if error == target));
    }
}

#[test]
fn g4_verifier_refusals_report_exact_partial_work_and_callback_errors_win() {
    let p = problem(
        "checkpoint-refusal",
        poly(vec![0.0, 1.0, -1.0]),
        vec![0.0, 0.5, 1.0],
    );

    let mut shape_callbacks = 0;
    let shape_refusal = verify_with_checkpoint(&p, &[0.0, 0.0], 1.0, |_| {
        shape_callbacks += 1;
        Ok::<(), core::convert::Infallible>(())
    })
    .expect("infallible callback");
    assert_eq!(
        shape_refusal.refusal,
        Some(VerifierRefusal::CandidateLength)
    );
    assert_eq!(shape_callbacks, 0);

    let mut invalid_tolerance_trace = Vec::new();
    let invalid_tolerance = verify_with_checkpoint(&p, &[0.0, 0.0, 0.0], f64::NAN, |progress| {
        invalid_tolerance_trace.push(progress);
        Ok::<(), core::convert::Infallible>(())
    })
    .expect("infallible callback");
    assert_eq!(
        invalid_tolerance.refusal,
        Some(VerifierRefusal::InvalidTolerance)
    );
    assert_eq!(
        invalid_tolerance_trace,
        vec![
            VerifierProgress {
                kind: VerifierCheckpointKind::PhaseEntry,
                phase: VerifierPhase::Validation,
                completed_work_units: 0,
                planned_work_units: 48,
            },
            VerifierProgress {
                kind: VerifierCheckpointKind::RefusalFlush,
                phase: VerifierPhase::Validation,
                completed_work_units: 3,
                planned_work_units: 48,
            },
        ]
    );

    let callback_wins = verify_with_checkpoint(&p, &[0.0, 0.0, 0.0], f64::NAN, |progress| {
        if progress.kind == VerifierCheckpointKind::RefusalFlush {
            Err("cancel-at-refusal")
        } else {
            Ok(())
        }
    });
    assert!(matches!(callback_wins, Err("cancel-at-refusal")));

    let mut candidate = [0.0, 0.0, 0.0];
    candidate[1] = f64::NAN;
    let mut candidate_trace = Vec::new();
    let candidate_refusal = verify_with_checkpoint(&p, &candidate, 1.0, |progress| {
        candidate_trace.push(progress);
        Ok::<(), core::convert::Infallible>(())
    })
    .expect("infallible callback");
    assert_eq!(
        candidate_refusal.refusal,
        Some(VerifierRefusal::CandidateNonFinite)
    );
    assert_eq!(candidate_trace.len(), 2);
    assert_eq!(
        candidate_trace[1],
        VerifierProgress {
            kind: VerifierCheckpointKind::RefusalFlush,
            phase: VerifierPhase::Validation,
            completed_work_units: 10,
            planned_work_units: 48,
        }
    );
}

#[test]
#[allow(clippy::too_many_lines)] // One exact invocation-global trace across every phase.
fn g4_verifier_work_boundaries_are_invocation_global_and_deterministic() {
    let mesh: Vec<f64> = (0..=512).map(|index| f64::from(index) / 512.0).collect();
    let p = problem("checkpoint-large", poly(vec![0.0]), mesh);
    let candidate = vec![0.0; 513];
    let plan = VerifierWorkPlan::for_inputs(&p, &candidate).expect("bounded verifier shape");
    assert_eq!(plan.identity_fields(), [1549, 512, 512, 1035, 1, 3609]);

    let mut trace = Vec::new();
    let report = verify_with_checkpoint(&p, &candidate, 1.0, |progress| {
        trace.push(progress);
        Ok::<(), core::convert::Infallible>(())
    })
    .expect("infallible callback");
    assert!(report.refusal.is_none());
    let expected = vec![
        (
            VerifierCheckpointKind::PhaseEntry,
            VerifierPhase::Validation,
            0,
        ),
        (
            VerifierCheckpointKind::WorkBoundary,
            VerifierPhase::Validation,
            256,
        ),
        (
            VerifierCheckpointKind::WorkBoundary,
            VerifierPhase::Validation,
            512,
        ),
        (
            VerifierCheckpointKind::WorkBoundary,
            VerifierPhase::Validation,
            768,
        ),
        (
            VerifierCheckpointKind::WorkBoundary,
            VerifierPhase::Validation,
            1024,
        ),
        (
            VerifierCheckpointKind::WorkBoundary,
            VerifierPhase::Validation,
            1280,
        ),
        (
            VerifierCheckpointKind::WorkBoundary,
            VerifierPhase::Validation,
            1536,
        ),
        (
            VerifierCheckpointKind::PhaseEntry,
            VerifierPhase::Tightness,
            1549,
        ),
        (
            VerifierCheckpointKind::WorkBoundary,
            VerifierPhase::Tightness,
            1792,
        ),
        (
            VerifierCheckpointKind::WorkBoundary,
            VerifierPhase::Tightness,
            2048,
        ),
        (
            VerifierCheckpointKind::PhaseEntry,
            VerifierPhase::Equilibrated,
            2061,
        ),
        (
            VerifierCheckpointKind::WorkBoundary,
            VerifierPhase::Equilibrated,
            2304,
        ),
        (
            VerifierCheckpointKind::WorkBoundary,
            VerifierPhase::Equilibrated,
            2560,
        ),
        (
            VerifierCheckpointKind::PhaseEntry,
            VerifierPhase::Hash,
            2573,
        ),
        (
            VerifierCheckpointKind::WorkBoundary,
            VerifierPhase::Hash,
            2816,
        ),
        (
            VerifierCheckpointKind::WorkBoundary,
            VerifierPhase::Hash,
            3072,
        ),
        (
            VerifierCheckpointKind::WorkBoundary,
            VerifierPhase::Hash,
            3328,
        ),
        (
            VerifierCheckpointKind::WorkBoundary,
            VerifierPhase::Hash,
            3584,
        ),
        (
            VerifierCheckpointKind::PhaseEntry,
            VerifierPhase::Finalization,
            3608,
        ),
        (
            VerifierCheckpointKind::Publication,
            VerifierPhase::Finalization,
            3609,
        ),
    ];
    assert_eq!(trace.len(), expected.len());
    for (progress, (kind, phase, completed_work_units)) in trace.iter().zip(&expected) {
        assert_eq!(progress.kind, *kind);
        assert_eq!(progress.phase, *phase);
        assert_eq!(progress.completed_work_units, *completed_work_units);
        assert_eq!(progress.planned_work_units, 3609);
    }

    let boundary = trace
        .iter()
        .copied()
        .find(|progress| progress.completed_work_units == 2304)
        .expect("equilibrated global boundary");
    let interrupted = verify_with_checkpoint(&p, &candidate, 1.0, |progress| {
        if progress == boundary {
            Err(boundary)
        } else {
            Ok(())
        }
    });
    assert!(matches!(interrupted, Err(error) if error == boundary));
}

#[test]
fn g4_verifier_publication_is_distinct_when_final_work_hits_a_boundary() {
    #[allow(clippy::cast_precision_loss)]
    let mesh: Vec<f64> = (0..179).map(|index| f64::from(index) / 178.0).collect();
    let p = problem(
        "checkpoint-publication-boundary",
        poly(vec![0.0, 1.0, -1.0]),
        mesh,
    );
    let candidate = vec![0.0; 179];
    let plan = VerifierWorkPlan::for_inputs(&p, &candidate).expect("bounded verifier shape");
    assert_eq!(plan.planned_work_units(), 1280);

    let mut trace = Vec::new();
    let report = verify_with_checkpoint(&p, &candidate, 10.0, |progress| {
        trace.push(progress);
        Ok::<(), core::convert::Infallible>(())
    })
    .expect("infallible callback");
    assert!(report.refusal.is_none());
    assert_eq!(trace[trace.len() - 2].completed_work_units, 1280);
    assert_eq!(
        trace[trace.len() - 2].kind,
        VerifierCheckpointKind::WorkBoundary
    );
    assert_eq!(trace[trace.len() - 1].completed_work_units, 1280);
    assert_eq!(
        trace[trace.len() - 1].kind,
        VerifierCheckpointKind::Publication
    );
}

/// ver-001 — THE UPPER-BOUND PROPERTY (G1 MMS class): over the battery
/// AND adversarially perturbed candidates (the untrusted-proposer
/// case: Prager–Synge holds for ANY conforming candidate), the bound
/// never underestimates the oracle truth. Exact-solution input stays
/// nonnegative and is not falsely rejected.
#[test]
fn ver_001_upper_bound_property() {
    let mut rng = Lcg(0x1001_2026_0707_0091);
    let mut checks = 0u32;
    let mut violations = 0u32;
    for (name, u) in mms_zoo() {
        for mesh in meshes() {
            let p = problem(name, u.clone(), mesh);
            let galerkin = solve_p1(&p);
            let mut candidates = vec![galerkin.clone()];
            // Untrusted proposers: noisy variants (BCs preserved).
            for _ in 0..3 {
                let mut noisy = galerkin.clone();
                for v in noisy
                    .iter_mut()
                    .skip(1)
                    .take(p.mesh().len().saturating_sub(2))
                {
                    *v += (rng.unit() - 0.5) * 0.02;
                }
                candidates.push(noisy);
            }
            for cand in candidates {
                let rep = verify(&p, &cand, 1e-3);
                let truth = true_energy_error(&p, &cand);
                checks += 1;
                // Oracle slack: the oracle itself is f64 quadrature.
                if rep.bound.hi < truth * (1.0 - 1e-9) {
                    violations += 1;
                }
            }
        }
    }
    // Exact zero solution: bound ≥ 0, accepted at any tolerance.
    let zero = problem("zero", poly(vec![0.0]), vec![0.0, 0.5, 1.0]);
    let z = verify(&zero, &[0.0, 0.0, 0.0], 1e-12);
    let zero_ok = z.accept && z.bound.hi >= 0.0 && z.color.is_some();
    verdict(
        "ver-001",
        violations == 0 && checks > 100 && zero_ok,
        &format!(
            "the equilibrated bound dominated the oracle truth on {checks}/{checks} \
             checks across 5 MMS solutions x 6 meshes x {{Galerkin + 3 adversarial \
             perturbed candidates}} — Prager-Synge holds for ANY conforming \
             candidate, which is exactly what makes untrusted proposers safe; the \
             exact-zero input accepts with a verified color; \
             seed 0x1001_2026_0707_0091"
        ),
    );
}

/// ver-002 — EFFECTIVITY (the kill-criterion's tightness leg): median
/// bound/truth on the Galerkin battery within the stated band;
/// loose-but-sound cases are logged as TIGHTNESS failures.
#[test]
fn ver_002_effectivity_band() {
    let mut effs = Vec::new();
    let mut tightness_failures = 0u32;
    for (name, u) in mms_zoo() {
        let _ = name;
        for mesh in meshes() {
            if mesh.len() < 4 {
                continue; // effectivity on trivial meshes is noise
            }
            let p = problem(name, u.clone(), mesh);
            let cand = solve_p1(&p);
            let rep = verify(&p, &cand, 1e-3);
            let eff = effectivity(&p, &cand, &rep);
            if eff > 5.0 {
                tightness_failures += 1;
            }
            effs.push(eff);
        }
    }
    effs.sort_by(|a, b| a.partial_cmp(b).expect("finite"));
    let median = effs[effs.len() / 2];
    let mut em = fs_obs::Emitter::new("fs-verify/conformance", "ver-002/effectivity");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "verify-effectivity".to_string(),
                json: format!(
                    "{{\"median\":{median:.4},\"min\":{:.4},\"max\":{:.4},\
                     \"tightness_failures\":{tightness_failures},\"n\":{}}}",
                    effs.first().expect("nonempty"),
                    effs.last().expect("nonempty"),
                    effs.len()
                ),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("effectivity event validates");
    println!("{line}");
    verdict(
        "ver-002",
        median <= 3.0 && tightness_failures == 0,
        &format!(
            "median effectivity {median:.3} (band <= 3; the accept economy is \
             unreachable with loose-but-sound bounds), zero tightness failures \
             over {} Galerkin cases",
            effs.len()
        ),
    );
}

/// ver-003 — interval soundness and FAIL CLOSED: the enclosure
/// contains a high-resolution oracle recomputation; NaN/∞ candidates
/// reject with no color; wild candidates stay finite.
#[test]
fn ver_003_interval_soundness_fail_closed() {
    let (name, u) = &mms_zoo()[2];
    let p = problem(name, u.clone(), meshes()[2].clone());
    let cand = solve_p1(&p);
    let rep = verify(&p, &cand, 1e-3);
    // ver-001 covers truth-side domination; here: the enclosure is a
    // genuine tight interval (nonnegative, near-ulp width).
    let width = rep.bound.hi - rep.bound.lo;
    let tight_enclosure = width >= 0.0 && width < 1e-10 * rep.bound.hi.max(1e-300) + 1e-14;
    // NaN candidate: fail closed.
    let mut nan_cand = cand.clone();
    nan_cand[1] = f64::NAN;
    let rn = verify(&p, &nan_cand, 1e-3);
    let nan_closed = !rn.accept && rn.color.is_none();
    // Infinite candidate: fail closed.
    let mut inf_cand = cand.clone();
    inf_cand[1] = f64::INFINITY;
    let ri = verify(&p, &inf_cand, 1e-3);
    let inf_closed = !ri.accept && ri.color.is_none();
    // Wild-but-finite candidate: finite bound, rejected, no overflow.
    let mut wild = cand.clone();
    wild[1] = 1e12;
    let rw = verify(&p, &wild, 1e-3);
    let wild_ok = rw.bound.hi.is_finite() && !rw.accept && rw.color.is_none();
    verdict(
        "ver-003",
        tight_enclosure && nan_closed && inf_closed && wild_ok,
        &format!(
            "the enclosure is tight (width {width:.2e}), NaN and infinite \
             candidates FAIL CLOSED (reject, no color — never a badge without a \
             bound), and a 1e12 spike stays finite and rejected"
        ),
    );
}

/// ver-004 — G5 determinism and boundary meshes: bit-identical bound
/// endpoints and verdicts across repeated runs; the single-interior-DOF
/// and no-interior-DOF meshes behave.
#[test]
fn ver_004_determinism_and_boundaries() {
    let (name, u) = &mms_zoo()[1];
    let p = problem(name, u.clone(), meshes()[1].clone());
    let cand = solve_p1(&p);
    let (r1, r2) = (verify(&p, &cand, 1e-4), verify(&p, &cand, 1e-4));
    let bitwise = r1.bound.lo.to_bits() == r2.bound.lo.to_bits()
        && r1.bound.hi.to_bits() == r2.bound.hi.to_bits()
        && r1.accept == r2.accept
        && r1.flux_hash == r2.flux_hash;
    // Accept on exact equality of bound and tolerance is SOUND
    // (bound >= truth, so truth <= tol).
    let tol_eq = verify(&p, &cand, r1.bound.hi);
    let equality_accepts = tol_eq.accept;
    // Single interior DOF.
    let p1dof = problem(name, u.clone(), vec![0.0, 0.5, 1.0]);
    let c1 = solve_p1(&p1dof);
    let rep1 = verify(&p1dof, &c1, 1.0);
    let single_ok = rep1.bound.hi >= true_energy_error(&p1dof, &c1) * (1.0 - 1e-9);
    // No interior DOF (2 nodes): the zero candidate is all we have.
    let p0dof = problem(name, u.clone(), vec![0.0, 1.0]);
    let rep0 = verify(&p0dof, &[0.0, 0.0], 10.0);
    let none_ok = rep0.bound.hi >= true_energy_error(&p0dof, &[0.0, 0.0]) * (1.0 - 1e-9);
    verdict(
        "ver-004",
        bitwise && equality_accepts && single_ok && none_ok,
        "verdicts, bound endpoints, and flux hashes are BITWISE reproducible; \
         accepting on exact bound==tolerance is sound by domination; single- and \
         zero-interior-DOF meshes still bound truthfully",
    );
}

/// ver-005 — certify-the-certifiers + the falsifier: an injected
/// UNSOUND estimator (bound/10) is CAUGHT by the MMS harness (Sev-0
/// machinery works); the independent hierarchical family agrees with
/// the equilibrated bound within its stated band.
#[test]
fn ver_005_certify_the_certifiers() {
    let mut caught = 0u32;
    let mut ratio_ok = true;
    let mut ratios = Vec::new();
    for (name, u) in mms_zoo() {
        for mesh in meshes() {
            if mesh.len() < 5 {
                continue;
            }
            let p = problem(name, u.clone(), mesh);
            let cand = solve_p1(&p);
            let rep = verify(&p, &cand, 1e-3);
            let truth = true_energy_error(&p, &cand);
            // The deliberately unsound estimator: bound / 10.
            let unsound = rep.bound.hi / 10.0;
            if unsound < truth * (1.0 - 1e-9) && truth > 1e-13 {
                caught += 1; // the harness detects the undershoot
            }
            // Falsifier: hierarchical family must not contradict.
            let hier = hierarchical_estimate(&p, &cand);
            if truth > 1e-12 {
                let ratio = hier / rep.bound.hi;
                ratios.push(ratio);
                ratio_ok &= (0.15..=1.2).contains(&ratio);
            }
        }
    }
    verdict(
        "ver-005",
        caught > 10 && ratio_ok,
        &format!(
            "the injected unsound estimator (bound/10) undershoots truth and is \
             CAUGHT on {caught} battery cases (a fooled bound is a Sev-0 wrong \
             answer wearing a badge — the harness sees it), and the independent \
             {} family stays within its stated band of the equilibrated bound \
             ({} ratios in [0.15, 1.2])",
            EstimatorFamily::Hierarchical.id(),
            ratios.len()
        ),
    );
}

/// ver-006 — the nonlinear WARM-START fallback: measured iteration
/// savings with an ESTIMATED color (never verified), plus the full
/// review-round-3 ledger rows.
#[test]
fn ver_006_warm_start_and_ledger() {
    let (name, u) = &mms_zoo()[2];
    let p = problem(name, u.clone(), meshes()[2].clone());
    let cand = solve_p1(&p);
    let ws = warm_start(&p, &cand, 50);
    let saves = f64::from(ws.cold_iterations) / f64::from(ws.warm_iterations.max(1));
    let color_honest = matches!(ws.color, fs_evidence::Color::Estimated { .. });
    // Ledger rows for a battery slice.
    let rep = verify(&p, &cand, 1e-3);
    let truth = true_energy_error(&p, &cand);
    let row = rep.to_row(p.name(), truth);
    let row_complete = row.contains("estimator_family_id")
        && row.contains("flux_hash")
        && row.contains("bound_lo")
        && row.contains("bound_hi")
        && row.contains("oracle_true_error")
        && row.contains("effectivity")
        && row.contains("verdict")
        && row.contains("tolerance");
    let mut em = fs_obs::Emitter::new("fs-verify/conformance", "ver-006/ledger");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "verify-ledger-row".to_string(),
                json: row.clone(),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("ledger row validates");
    println!("{line}");
    verdict(
        "ver-006",
        saves >= 1.5 && color_honest && row_complete,
        &format!(
            "the warm start saves {:.1}x Newton iterations ({} cold vs {} warm) and \
             carries an ESTIMATED color — never a certificate (the honest R1 \
             boundary); the ledger row carries every review-round-3 field: {row}",
            saves, ws.cold_iterations, ws.warm_iterations
        ),
    );
}

fn assert_refused(
    problem: &MmsProblem,
    candidate: &[f64],
    tolerance: f64,
    expected: VerifierRefusal,
) {
    let report = verify(problem, candidate, tolerance);
    assert_eq!(report.refusal, Some(expected));
    assert!(!report.accept);
    assert!(report.color.is_none());
    assert!(report.bound.is_unbounded());
    assert_eq!(report.flux_hash, 0);
}

/// ver-007 — public-input admission is complete and precedes all verifier
/// indexing/compute. Every malformed input is a structured refusal with an
/// unbounded sentinel and no evidence color, never a panic or partial badge.
#[test]
#[allow(clippy::too_many_lines)] // one adversarial admission matrix
fn ver_007_hostile_public_inputs_fail_closed() {
    let exact = poly(vec![0.0, 1.0, -1.0]);
    let base = problem("hostile", exact.clone(), vec![0.0, 0.5, 1.0]);
    let candidate = vec![0.0, 0.25, 0.0];

    for mesh in [Vec::new(), vec![0.0]] {
        assert!(matches!(
            MmsProblem::new("short-mesh", exact.clone(), mesh),
            Err(Fem1dError::ResourceLimit {
                resource: "mesh nodes",
                ..
            })
        ));
    }
    assert_refused(&base, &[0.0, 0.0], 1.0, VerifierRefusal::CandidateLength);

    for tolerance in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_refused(
            &base,
            &candidate,
            tolerance,
            VerifierRefusal::InvalidTolerance,
        );
    }

    let signed_zero_mesh = problem("signed-zero", exact.clone(), vec![-0.0, 0.5, 1.0]);
    let canonical_mesh = problem("signed-zero", exact.clone(), vec![0.0, 0.5, 1.0]);
    assert_eq!(signed_zero_mesh.identity(), canonical_mesh.identity());
    for mesh in [vec![0.0, 0.5, 2.0], vec![0.0, 0.5, f64::NAN]] {
        assert!(matches!(
            MmsProblem::new("wrong-domain", exact.clone(), mesh),
            Err(Fem1dError::MeshDomain | Fem1dError::NonFiniteMeshNode { .. })
        ));
    }
    for mesh in [
        vec![0.0, 0.5, 0.5, 1.0],
        vec![0.0, f64::NAN, 1.0],
        vec![0.0, f64::INFINITY, 1.0],
    ] {
        assert!(matches!(
            MmsProblem::new("bad-coordinates", exact.clone(), mesh),
            Err(Fem1dError::NonIncreasingMeshCell { .. } | Fem1dError::NonFiniteMeshNode { .. })
        ));
    }

    for values in [
        vec![1.0, 0.25, 0.0],
        vec![0.0, 0.25, 1.0],
        vec![-0.0, 0.25, 0.0],
    ] {
        assert_refused(&base, &values, 1.0, VerifierRefusal::CandidateBoundary);
    }
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut values = candidate.clone();
        values[1] = value;
        assert_refused(&base, &values, 1.0, VerifierRefusal::CandidateNonFinite);
    }

    assert!(matches!(
        MmsProblem::new(
            "oversized-mesh",
            exact.clone(),
            vec![0.0; MAX_FEM1D_MESH_NODES + 1],
        ),
        Err(Fem1dError::ResourceLimit {
            resource: "mesh nodes",
            ..
        })
    ));
    let mut too_many_semantic_coefficients = vec![0.0; MAX_FEM1D_POLY_COEFFICIENTS + 1];
    *too_many_semantic_coefficients.last_mut().unwrap() = 1.0;
    assert!(matches!(
        Poly::new(too_many_semantic_coefficients),
        Err(Fem1dError::PolynomialCoefficientCount { .. })
    ));
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            Poly::new(vec![0.0, value]),
            Err(Fem1dError::NonFinitePolynomialCoefficient { .. })
        ));
    }

    // The stored binary64 coefficients cancel exactly at x=1, even though
    // ordinary Horner evaluation rounds to a nonzero value.
    let cancellation_boundary = problem(
        "cancellation-boundary",
        poly(vec![0.0, 1.0, 1.0e16, -1.0e16, -1.0]),
        base.mesh().to_vec(),
    );
    let cancellation_report = verify(&cancellation_boundary, &[0.0, 0.0, 0.0], 10.0);
    assert_eq!(cancellation_report.refusal, None);
    assert!(!cancellation_report.bound.is_unbounded());

    // Point Horner loses the final exact residue and returns zero, but the
    // binary-rational polynomial is nonzero at x=1 and must be refused.
    let hidden_residue = poly(vec![0.0, 1.0e16, -1.0e16, 1.0]);
    assert_eq!(hidden_residue.eval(1.0).to_bits(), 0.0_f64.to_bits());
    assert!(matches!(
        MmsProblem::new(
            "hidden-boundary-residue",
            hidden_residue,
            base.mesh().to_vec(),
        ),
        Err(Fem1dError::ExactSolutionBoundary)
    ));
    assert!(matches!(
        MmsProblem::new("nonvanishing", poly(vec![1.0]), base.mesh().to_vec()),
        Err(Fem1dError::ExactSolutionBoundary)
    ));

    let changed_solution = problem("hostile", poly(vec![0.0, 2.0, -2.0]), base.mesh().to_vec());
    assert_ne!(base.identity(), changed_solution.identity());
    assert_ne!(base.forcing(), changed_solution.forcing());

    let mut refused = verify(&base, &candidate, f64::NAN);
    assert!(matches!(
        try_effectivity(&base, &candidate, &refused),
        Err(Fem1dError::InvalidScalar {
            field: "verifier report",
            ..
        })
    ));
    let zero_problem = problem("zero", poly(vec![0.0]), base.mesh().to_vec());
    let zero_candidate = vec![0.0; zero_problem.mesh().len()];
    let zero_report = verify(&zero_problem, &zero_candidate, 1.0);
    assert!(matches!(
        try_effectivity(&zero_problem, &zero_candidate, &zero_report),
        Err(Fem1dError::InvalidScalar {
            field: "oracle true error",
            ..
        })
    ));
    refused.family = "family\"}\n{\"forged-family\":true";
    let row = refused.to_row("hostile\"}\n{\"forged\":true", f64::NAN);
    assert!(row.contains("\"verdict\":\"refused\""));
    assert!(row.contains("\"bound_hi\":null"));
    assert!(row.contains("\"tolerance\":null"));
    assert!(row.contains("hostile\\\"}\\n{\\\"forged\\\":true"));
    assert!(row.contains("family\\\"}\\n{\\\"forged-family\\\":true"));
    assert!(!row.contains('\n'));
    let mut emitter = fs_obs::Emitter::new("fs-verify/conformance", "ver-007/refusal-row");
    let line = emitter
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "verify-refusal".to_string(),
                json: row,
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("structured refusal row is valid JSON");
    verdict(
        "ver-007",
        true,
        "canonical construction refuses malformed classes/meshes before a problem exists; verifier candidate/tolerance faults still refuse without color; signed/trailing zeros normalize and semantic identity changes remain observable",
    );
}
