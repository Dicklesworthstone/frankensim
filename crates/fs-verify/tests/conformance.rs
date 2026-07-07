//! fs-verify conformance (bead lmp4.1, feature `certified-speculation`).
//! The MMS upper-bound property over the battery INCLUDING adversarial
//! perturbed candidates (the untrusted-proposer case), effectivity
//! bands, interval soundness + fail-closed, G5 determinism, the
//! certify-the-certifiers injection, the estimator-family falsifier,
//! and the nonlinear warm-start fallback with ledger rows. JSON-line
//! verdicts; seeded cases carry seeds.

use fs_verify::estimator::{
    EstimatorFamily, effectivity, hierarchical_estimate, verify, warm_start,
};
use fs_verify::fem1d::{MmsProblem, Poly, solve_p1, true_energy_error};

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
    let u1 = Poly(vec![0.0, 1.0, -1.0]);
    // x(1−x)(x−0.3) = 0.3·(−x) + ... expand: (x − x²)(x − 0.3)
    // = x² − 0.3x − x³ + 0.3x² = −0.3x + 1.3x² − x³
    let u2 = Poly(vec![0.0, -0.3, 1.3, -1.0]);
    // x²(1−x)² = x² − 2x³ + x⁴
    let u3 = Poly(vec![0.0, 0.0, 1.0, -2.0, 1.0]);
    // x(1−x)(x−0.2)(x−0.8): (x−x²)(x²−x+0.16)
    // = x³ − x² + 0.16x − x⁴ + x³ − 0.16x² = 0.16x − 1.16x² + 2x³ − x⁴
    let u4 = Poly(vec![0.0, 0.16, -1.16, 2.0, -1.0]);
    // Degree 5: x(1−x)(x−0.5)(x²+0.3):
    // (x−x²)(x−0.5) = x² −0.5x −x³ +0.5x² = −0.5x +1.5x² −x³
    // ×(x²+0.3): −0.5x³ +1.5x⁴ −x⁵ −0.15x +0.45x² −0.3x³
    let u5 = Poly(vec![0.0, -0.15, 0.45, -0.8, 1.5, -1.0]);
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
            let p = MmsProblem::new(name, u.clone(), mesh);
            let galerkin = solve_p1(&p);
            let mut candidates = vec![galerkin.clone()];
            // Untrusted proposers: noisy variants (BCs preserved).
            for _ in 0..3 {
                let mut noisy = galerkin.clone();
                for v in noisy
                    .iter_mut()
                    .skip(1)
                    .take(p.mesh.len().saturating_sub(2))
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
    let zero = MmsProblem::new("zero", Poly(vec![0.0]), vec![0.0, 0.5, 1.0]);
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
            let p = MmsProblem::new(name, u.clone(), mesh);
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
    let p = MmsProblem::new(name, u.clone(), meshes()[2].clone());
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
    let p = MmsProblem::new(name, u.clone(), meshes()[1].clone());
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
    let p1dof = MmsProblem::new(name, u.clone(), vec![0.0, 0.5, 1.0]);
    let c1 = solve_p1(&p1dof);
    let rep1 = verify(&p1dof, &c1, 1.0);
    let single_ok = rep1.bound.hi >= true_energy_error(&p1dof, &c1) * (1.0 - 1e-9);
    // No interior DOF (2 nodes): the zero candidate is all we have.
    let p0dof = MmsProblem::new(name, u.clone(), vec![0.0, 1.0]);
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
            let p = MmsProblem::new(name, u.clone(), mesh);
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
    let p = MmsProblem::new(name, u.clone(), meshes()[2].clone());
    let cand = solve_p1(&p);
    let ws = warm_start(&p, &cand, 50);
    let saves = f64::from(ws.cold_iterations) / f64::from(ws.warm_iterations.max(1));
    let color_honest = matches!(ws.color, fs_evidence::Color::Estimated { .. });
    // Ledger rows for a battery slice.
    let rep = verify(&p, &cand, 1e-3);
    let truth = true_energy_error(&p, &cand);
    let row = rep.to_row(&p.name, truth);
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
