//! Proposer-zoo conformance (bead lmp4.2, feature
//! `certified-speculation`). The type-level safety invariant, advisory
//! confidence in both directions, neighbor extrapolation with and
//! without warm adjoints (plus the equidistant tie-break), coarse-rung
//! prolongation with the fp16 precision-discipline demo, the
//! adversarial-surrogate falsifier with zero incorrect accepts and
//! auto-demotion, and the end-to-end economics loop with ledger rows.
//! Completed aggregates emit canonical fs-obs verdicts. The two randomized
//! campaigns carry their literal root seeds; fixed inputs use zero, and this
//! suite has no execution seed. Assertions and expectations reached before a
//! verdict remain ordinary Rust test diagnostics.

use fs_math::eft::two_sum;
use fs_verify::estimator::verify;
use fs_verify::fem1d::{
    Fem1dError, MmsProblem, Poly, solve_p1 as try_solve_p1,
    true_energy_error as try_true_energy_error,
};
use fs_verify::zoo::{
    CoarseRungProlongation, NeighborExtrapolation, Outcome, Proposal, Proposer, Registry,
    SpeculationQuery, ZooTelemetry, quantize_f16, speculate as try_speculate,
};

fn solve_p1(problem: &MmsProblem) -> Vec<f64> {
    try_solve_p1(problem).expect("zoo problem must solve")
}

fn poly(coefficients: Vec<f64>) -> Poly {
    Poly::new(coefficients).expect("valid zoo polynomial")
}

fn problem(name: &str, u: Poly, mesh: Vec<f64>) -> MmsProblem {
    MmsProblem::new(name, u, mesh).expect("valid zoo problem")
}

fn true_energy_error(problem: &MmsProblem, candidate: &[f64]) -> f64 {
    try_true_energy_error(problem, candidate).expect("zoo oracle must evaluate")
}

fn speculate(
    query: &SpeculationQuery,
    registry: &Registry,
    telemetry: &mut ZooTelemetry,
) -> Outcome {
    try_speculate(query, registry, telemetry).expect("zoo speculation must execute")
}

const SUITE: &str = "fs-verify/zoo";
const FIXED_INPUT_SEED: u64 = 0;
const ZOO_004_INPUT_SEED: u64 = 0x1001_2026_0707_00A4;
const ZOO_005_INPUT_SEED: u64 = 0x1001_2026_0707_00A5;

fn verdict(case: &str, pass: bool, detail: &str, seed: u64) {
    let mut emitter = fs_obs::Emitter::new(SUITE, case);
    let event = emitter.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: case.to_string(),
            pass,
            detail: detail.to_string(),
            seed,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("proposer-zoo verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("proposer-zoo verdict must use the fs-obs wire schema");
    println!("{line}");
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

/// The parameterized design family: u(x; θ) = x(1−x)(x−θ).
fn family(theta: f64) -> Poly {
    // `(1 + theta)` may round. Preserve its error term on x^4 so the stored
    // binary64 polynomial, not only the ideal formula, vanishes exactly at 1.
    let (middle, correction) = two_sum(1.0, theta);
    poly(vec![0.0, -theta, middle, -1.0, correction])
}

fn uniform(n: usize) -> Vec<f64> {
    (0..=n).map(|i| i as f64 / n as f64).collect()
}

fn query(theta: f64, n: usize, tol: f64) -> SpeculationQuery {
    SpeculationQuery {
        problem: problem("family", family(theta), uniform(n)),
        theta,
        tolerance: tol,
        regime: "wedge-v0".to_string(),
    }
}

/// A deliberately adversarial surrogate: garbage with max confidence.
struct AdversarialSurrogate;

impl Proposer for AdversarialSurrogate {
    fn name(&self) -> &'static str {
        "adversarial-surrogate"
    }

    fn propose(&self, q: &SpeculationQuery) -> Result<Option<Proposal>, Fem1dError> {
        let mut candidate: Vec<f64> = q.problem.mesh().iter().map(|x| (x * 941.0).sin()).collect();
        candidate[0] = 0.0;
        let last = candidate.len() - 1;
        candidate[last] = 0.0;
        Ok(Some(Proposal {
            candidate,
            confidence: 1.0,
        }))
    }
}

/// A good proposer that self-reports NaN confidence (the advisory
/// property: it must still be tried and accepted when it verifies).
struct HumbleGood;

impl Proposer for HumbleGood {
    fn name(&self) -> &'static str {
        "humble-good"
    }

    fn propose(&self, q: &SpeculationQuery) -> Result<Option<Proposal>, Fem1dError> {
        Ok(Some(Proposal {
            candidate: try_solve_p1(&q.problem)?,
            confidence: f64::NAN,
        }))
    }
}

struct QuotedIdentity;

impl Proposer for QuotedIdentity {
    fn name(&self) -> &'static str {
        "quoted\"\\proposer"
    }

    fn propose(&self, q: &SpeculationQuery) -> Result<Option<Proposal>, Fem1dError> {
        Ok(Some(Proposal {
            candidate: vec![0.0; q.problem.mesh().len()],
            confidence: 0.5,
        }))
    }
}

/// zoo-001 — the interface + THE SAFETY INVARIANT: answers exist only
/// through the verifier; empty registries are honest misses;
/// confidence is advisory in BOTH directions.
#[test]
fn zoo_001_interface_and_safety() {
    let mut telemetry = ZooTelemetry::default();
    // Empty registry: honest NoCandidates.
    let empty = Registry::new();
    let q = query(0.45, 16, 2e-1);
    let none = matches!(speculate(&q, &empty, &mut telemetry), Outcome::NoCandidates);
    // Register/deregister round-trip.
    let mut reg = Registry::new();
    reg.register(Box::new(CoarseRungProlongation))
        .expect("register coarse proposer");
    reg.register(Box::new(AdversarialSurrogate))
        .expect("register adversarial proposer");
    let has_two = reg.names().len() == 2;
    reg.deregister("adversarial-surrogate");
    let has_one = reg.names() == vec!["coarse-rung-prolongation"];
    // The accepted answer carries a VERIFIED color and its bound; the
    // only constructor is the verifier's yes.
    let out = speculate(&q, &reg, &mut telemetry);
    let sound = match &out {
        Outcome::Accepted(ans) => {
            ans.report().accept
                && matches!(
                    ans.report().color,
                    Some(fs_evidence::Color::Verified { .. })
                )
                && ans.report().bound.hi <= q.tolerance
        }
        _ => false,
    };
    // Advisory-NaN: the humble-good proposer (NaN confidence) is tried
    // LAST but still accepted when the noisy one fails.
    let mut reg2 = Registry::new();
    reg2.register(Box::new(HumbleGood))
        .expect("register humble proposer");
    reg2.register(Box::new(AdversarialSurrogate))
        .expect("register adversarial proposer");
    let tight = query(0.45, 16, 5e-2);
    let out2 = speculate(&tight, &reg2, &mut telemetry);
    let nan_still_wins = matches!(&out2, Outcome::Accepted(ans) if ans.proposer() == "humble-good");
    verdict(
        "zoo-001",
        none && has_two && has_one && sound && nan_still_wins,
        "empty registries miss honestly, registration hot-swaps, accepted answers \
         carry the verifier's bound and VERIFIED color (no other constructor \
         exists), and a NaN-confidence good proposer is ordered last yet still \
         accepted — confidence is advisory in both directions",
        FIXED_INPUT_SEED,
    );
}

/// zoo-002 — neighbor extrapolation: warm adjoints beat zeroth-order
/// measurably; equidistant neighbors tie-break deterministically to
/// the smaller θ; verified accepts at honest tolerances.
#[test]
fn zoo_002_neighbor_extrapolation() {
    let n = 32;
    // Certified cache at θ ∈ {0.2, 0.5, 0.8} with FD sensitivities.
    let solved = |th: f64| solve_p1(&problem("f", family(th), uniform(n)));
    let sens = |th: f64| -> Vec<f64> {
        let h = 1e-4;
        let (up, dn) = (solved(th + h), solved(th - h));
        up.iter()
            .zip(&dn)
            .map(|(a, b)| (a - b) / (2.0 * h))
            .collect()
    };
    let cache_warm: Vec<(f64, Vec<f64>, Option<Vec<f64>>)> = [0.2, 0.5, 0.8]
        .iter()
        .map(|&t| (t, solved(t), Some(sens(t))))
        .collect();
    let cache_cold: Vec<(f64, Vec<f64>, Option<Vec<f64>>)> = cache_warm
        .iter()
        .map(|(t, u, _)| (*t, u.clone(), None))
        .collect();
    let q = query(0.45, n, 1e-1);
    let warm = NeighborExtrapolation { cache: cache_warm }
        .propose(&q)
        .expect("warm executes")
        .expect("warm");
    let cold = NeighborExtrapolation { cache: cache_cold }
        .propose(&q)
        .expect("cold executes")
        .expect("cold");
    let warm_err = true_energy_error(&q.problem, &warm.candidate);
    let cold_err = true_energy_error(&q.problem, &cold.candidate);
    let warm_wins = warm_err < 0.5 * cold_err;
    // Both still pass the verifier at a loose tolerance (graceful
    // degradation to zeroth-order remains USEFUL).
    let warm_accepts = verify(&q.problem, &warm.candidate, 1e-1).accept;
    let cold_accepts = verify(&q.problem, &cold.candidate, 2e-1).accept;
    // Equidistant tie: θ = 0.35 sits exactly between 0.2 and 0.5; the
    // rule picks the SMALLER θ, deterministically.
    let qt = query(0.35, n, 1e-1);
    let cache2: Vec<(f64, Vec<f64>, Option<Vec<f64>>)> =
        [0.2, 0.5].iter().map(|&t| (t, solved(t), None)).collect();
    let pick1 = NeighborExtrapolation {
        cache: cache2.clone(),
    }
    .propose(&qt)
    .expect("tie 1 executes")
    .expect("tie 1");
    let pick2 = NeighborExtrapolation { cache: cache2 }
        .propose(&qt)
        .expect("tie 2 executes")
        .expect("tie 2");
    let tie_deterministic = pick1.candidate == pick2.candidate && pick1.candidate == solved(0.2); // zeroth-order from θ=0.2
    verdict(
        "zoo-002",
        warm_wins && warm_accepts && cold_accepts && tie_deterministic,
        &format!(
            "the warm adjoint cuts extrapolation error to {warm_err:.1e} vs \
             zeroth-order {cold_err:.1e} (>2x), both degrade gracefully into \
             verified accepts at honest tolerances, and the equidistant tie \
             resolves deterministically to the smaller theta"
        ),
        FIXED_INPUT_SEED,
    );
}

/// zoo-003 — coarse-rung prolongation + the PRECISION DISCIPLINE:
/// accepts at loose tolerance, rejects honestly at tight tolerance,
/// and an fp16-quantized candidate still verifies — speculate LOW,
/// verify HIGH.
#[test]
fn zoo_003_coarse_rung_and_precision() {
    let q_loose = query(0.3, 32, 1e-1);
    let prop = CoarseRungProlongation
        .propose(&q_loose)
        .expect("coarse executes")
        .expect("coarse");
    let loose = verify(&q_loose.problem, &prop.candidate, 1e-1);
    let tight = verify(&q_loose.problem, &prop.candidate, 1e-6);
    // fp16 quantization: the proposer's precision is nobody's business.
    let quantized: Vec<f64> = prop.candidate.iter().map(|&v| quantize_f16(v)).collect();
    let q_accept = verify(&q_loose.problem, &quantized, 1e-1);
    // Tiny meshes have no coarser rung: honest decline.
    let q_small = query(0.3, 3, 1e-1);
    let declines = CoarseRungProlongation
        .propose(&q_small)
        .expect("small coarse proposal executes")
        .is_none();
    verdict(
        "zoo-003",
        loose.accept && !tight.accept && q_accept.accept && declines,
        &format!(
            "the prolongated coarse solve accepts at 1e-1 (bound {:.2e}) and rejects \
             honestly at 1e-6; the fp16-QUANTIZED candidate still accepts (speculate \
             low, verify high — the certificate inherits the VERIFIER's precision); \
             mesh-too-small declines honestly",
            loose.bound.hi
        ),
        FIXED_INPUT_SEED,
    );
}

/// zoo-004 — THE FALSIFIER: an adversarial surrogate never lands a
/// single incorrect accept over the battery, its accept-rate collapse
/// AUTO-DEMOTES it in the regime, and demoted proposers stop being
/// consulted.
#[test]
fn zoo_004_adversarial_falsifier() {
    let mut telemetry = ZooTelemetry::default();
    let mut reg = Registry::new();
    reg.register(Box::new(AdversarialSurrogate))
        .expect("register adversarial proposer");
    reg.register(Box::new(CoarseRungProlongation))
        .expect("register coarse proposer");
    let mut rng = Lcg(ZOO_004_INPUT_SEED);
    let mut incorrect_accepts = 0u32;
    for _ in 0..25 {
        let theta = 0.2 + 0.6 * rng.unit();
        let q = query(theta, 32, 5e-2);
        match speculate(&q, &reg, &mut telemetry) {
            Outcome::Accepted(ans) => {
                // Cross-check the accept against the oracle: an
                // accepted bound must dominate the true error.
                let truth = true_energy_error(&q.problem, ans.candidate());
                if truth > ans.report().bound.hi * (1.0 + 1e-9) {
                    incorrect_accepts += 1;
                }
                if ans.proposer() == "adversarial-surrogate" {
                    incorrect_accepts += 1; // garbage must never verify
                }
            }
            Outcome::AllRejected(_) | Outcome::NoCandidates => {}
        }
    }
    let adv_rate = telemetry
        .accept_rate("adversarial-surrogate", "wedge-v0")
        .expect("tried");
    let demotions = telemetry
        .demote_collapsed(0.05, 10)
        .expect("valid demotion policy");
    let demoted = telemetry.is_demoted("adversarial-surrogate", "wedge-v0");
    // After demotion the adversary is not consulted (tries frozen).
    let tries_before = 25;
    let q = query(0.5, 32, 5e-2);
    let _ = speculate(&q, &reg, &mut telemetry);
    let frozen = telemetry
        .accept_rate("adversarial-surrogate", "wedge-v0")
        .is_some()
        && telemetry.rows().iter().any(|r| {
            r.contains("adversarial-surrogate") && r.contains(&format!("\"tries\":{tries_before}"))
        });
    verdict(
        "zoo-004",
        incorrect_accepts == 0 && adv_rate == 0.0 && !demotions.is_empty() && demoted && frozen,
        &format!(
            "zero incorrect accepts over 25 adversarial-first speculations (the \
             verifier gates everything), the adversary's accept rate is {adv_rate:.2}, \
             the collapse AUTO-DEMOTES it in the regime, and demoted proposers stop \
             being consulted; seed {ZOO_004_INPUT_SEED:#x}"
        ),
        ZOO_004_INPUT_SEED,
    );
}

/// zoo-005 — the economics loop end-to-end: a mixed registry over a
/// seeded query stream stays sound, telemetry rows ship to the ledger,
/// and the accept-rate ordering matches proposer quality.
#[test]
fn zoo_005_economics_loop() {
    let n = 32;
    let solved = |th: f64| solve_p1(&problem("f", family(th), uniform(n)));
    let cache: Vec<(f64, Vec<f64>, Option<Vec<f64>>)> = [0.2, 0.4, 0.6, 0.8]
        .iter()
        .map(|&t| (t, solved(t), None))
        .collect();
    let mut reg = Registry::new();
    reg.register(Box::new(NeighborExtrapolation { cache }))
        .expect("register neighbor proposer");
    reg.register(Box::new(CoarseRungProlongation))
        .expect("register coarse proposer");
    reg.register(Box::new(AdversarialSurrogate))
        .expect("register adversarial proposer");
    let mut telemetry = ZooTelemetry::default();
    let mut rng = Lcg(ZOO_005_INPUT_SEED);
    let mut accepted = 0u32;
    for _ in 0..30 {
        let theta = 0.25 + 0.5 * rng.unit();
        let q = query(theta, n, 8e-2);
        match speculate(&q, &reg, &mut telemetry) {
            Outcome::Accepted(ans) => {
                accepted += 1;
                assert!(ans.report().accept, "type invariant");
            }
            Outcome::AllRejected(_) | Outcome::NoCandidates => {}
        }
    }
    let rows = telemetry.rows();
    let mut em = fs_obs::Emitter::new("fs-verify/zoo", "zoo-005/economics");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "speculation-accept-rates".to_string(),
                json: format!(
                    "{{\"rows\":[{}],\"input_seed\":{ZOO_005_INPUT_SEED}}}",
                    rows.join(",")
                ),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("economics rows validate");
    println!("{line}");
    let adv = telemetry
        .accept_rate("adversarial-surrogate", "wedge-v0")
        .unwrap_or(1.0);
    let good_beats_bad = telemetry
        .accept_rate("neighbor-extrapolation", "wedge-v0")
        .into_iter()
        .chain(telemetry.accept_rate("coarse-rung-prolongation", "wedge-v0"))
        .any(|r| r > adv);
    verdict(
        "zoo-005",
        accepted > 15 && adv == 0.0 && good_beats_bad,
        &format!(
            "{accepted}/30 speculations accepted with certified bounds, the \
             adversary landed nothing, honest proposers out-rate it, and the \
             per-proposer-per-regime rows ship to the ledger; \
             seed {ZOO_005_INPUT_SEED:#x}"
        ),
        ZOO_005_INPUT_SEED,
    );
}

/// zoo-006 — invalid queries and built-in proposer state propagate structured
/// errors before proposal work; no panic or ordinary miss hides the failure.
#[test]
fn zoo_006_invalid_queries_refuse_before_proposal_work() {
    let mut telemetry = ZooTelemetry::default();
    let mut registry = Registry::new();
    registry
        .register(Box::new(AdversarialSurrogate))
        .expect("register adversarial proposer");

    let mut nonfinite_theta = query(0.4, 8, 1e-2);
    nonfinite_theta.theta = f64::NAN;
    assert!(matches!(
        try_speculate(&nonfinite_theta, &registry, &mut telemetry),
        Err(Fem1dError::InvalidScalar { field: "theta", .. })
    ));

    assert!(matches!(
        MmsProblem::new("malformed", family(0.4), Vec::new()),
        Err(Fem1dError::ResourceLimit {
            resource: "mesh nodes",
            ..
        })
    ));

    let q = query(0.4, 8, 1e-2);
    let poisoned_neighbor = NeighborExtrapolation {
        cache: vec![(f64::NAN, vec![0.0; q.problem.mesh().len()], None)],
    };
    assert!(matches!(
        poisoned_neighbor.propose(&q),
        Err(Fem1dError::InvalidScalar {
            field: "neighbor cache theta",
            ..
        })
    ));

    let mut control_regime = query(0.4, 8, 1e-2);
    control_regime.regime = "bad\nregime".to_string();
    assert!(matches!(
        try_speculate(&control_regime, &registry, &mut telemetry),
        Err(Fem1dError::InvalidScalar {
            field: "regime",
            ..
        })
    ));
    verdict(
        "zoo-006",
        true,
        "non-finite query coordinates and poisoned neighbor identities return structured errors before ordering, solving, or verification; malformed meshes cannot construct a query problem",
        FIXED_INPUT_SEED,
    );
}

/// zoo-007 — valid identity punctuation is deterministically JSON-escaped in
/// telemetry rows; it cannot inject a sibling field or line.
#[test]
fn zoo_007_telemetry_rows_escape_identities() {
    let mut registry = Registry::new();
    registry
        .register(Box::new(QuotedIdentity))
        .expect("register quoted identity");
    let mut telemetry = ZooTelemetry::default();
    let mut q = query(0.4, 8, 1e-12);
    q.regime = "quoted\"\\regime".to_string();
    let outcome = try_speculate(&q, &registry, &mut telemetry)
        .expect("quoted identities must remain valid data");
    assert!(matches!(outcome, Outcome::AllRejected(_)));
    let rows = telemetry.rows();
    assert_eq!(rows.len(), 1);
    assert!(rows[0].contains("quoted\\\"\\\\proposer"));
    assert!(rows[0].contains("quoted\\\"\\\\regime"));
    let mut emitter = fs_obs::Emitter::new("fs-verify/zoo", "zoo-007/identity-json");
    let line = emitter
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "quoted-identity-row".to_string(),
                json: rows[0].clone(),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("escaped identity row is valid JSON");
    println!("{line}");
    verdict(
        "zoo-007",
        true,
        "quoted and backslashed proposer/regime identities remain one valid JSON object with no field or line injection",
        FIXED_INPUT_SEED,
    );
}

/// zoo-008 — coarse-rung routing remains linear and correct on a large,
/// strongly nonuniform mesh; every injected coarse node is reproduced exactly.
#[test]
fn zoo_008_coarse_rung_routes_large_nonuniform_mesh() {
    let cells = 4_096usize;
    let mesh: Vec<f64> = (0..=cells)
        .map(|index| {
            let coordinate = index as f64 / cells as f64;
            coordinate * coordinate
        })
        .collect();
    let fine_problem = problem("nonuniform-coarse-rung", family(0.3), mesh.clone());
    let query = SpeculationQuery {
        problem: fine_problem,
        theta: 0.3,
        tolerance: 1.0,
        regime: "nonuniform-routing".to_string(),
    };
    let proposal = CoarseRungProlongation
        .propose(&query)
        .expect("large nonuniform proposal must execute")
        .expect("large nonuniform mesh has a coarse rung");

    let coarse_mesh: Vec<f64> = mesh.iter().step_by(2).copied().collect();
    let coarse_problem = problem("nonuniform-coarse-rung", family(0.3), coarse_mesh);
    let coarse_solution = solve_p1(&coarse_problem);
    assert_eq!(proposal.candidate.len(), mesh.len());
    assert!(proposal.candidate.iter().all(|value| value.is_finite()));
    for (coarse_index, value) in coarse_solution.iter().enumerate() {
        assert_eq!(
            proposal.candidate[coarse_index * 2].to_bits(),
            value.to_bits(),
            "coarse node {coarse_index} changed during prolongation"
        );
    }
    verdict(
        "zoo-008",
        true,
        "the monotone segment cursor prolongated 4,097 strongly nonuniform nodes and preserved every injected coarse-node value bitwise",
        FIXED_INPUT_SEED,
    );
}
