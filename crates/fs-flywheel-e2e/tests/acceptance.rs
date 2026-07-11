//! EPISTEMIC-ENGINE ACCEPTANCE (bead xpck.8): the TOP-LEVEL runnable
//! proof that the whole addendum works FOR A USER — a declarative
//! query becomes a colored, priced, auditable answer that a third
//! party re-verifies WITHOUT trusting the vendor. "The product is
//! justified belief at minimum cost", executed.
//!
//! The path: admission (typed, teaching refusals) → flywheel discharge
//! (planner + cache) → anytime colored answer (+ the VoI-priced hint)
//! → signed evidence package → SOLVER-FREE third-party re-check →
//! G5 whole-path replay → the laundering invariant at every hop.
#![cfg(feature = "flywheel-e2e")]

use fs_evidence::{Color, IntervalOp, compose};
use fs_ir::planner::{CostTable, MemCache, ProblemFamily};
use fs_ir::{admission, sexpr};
use fs_package::{Claim, EvidencePackage, Provenance};

/// A deliberately corrupted content root (one byte flipped): the v4
/// 32-byte replacement for the old `root ^ 0xdead` tamper idiom.
fn flip(root: fs_package::ContentHash) -> fs_package::ContentHash {
    let mut bytes = *root.as_bytes();
    bytes[0] ^= 0xde;
    fs_package::ContentHash(bytes)
}
use fs_verify::fem1d::Poly;
use std::collections::BTreeMap;

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-flywheel-e2e/acceptance\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

/// The user's wedge query as a declarative study (the CHT wedge with a
/// priced budget — "answer for under the budget or teach me why not").
const WEDGE: &str = r#"(study "cht-wedge-acceptance"
  (seed 0x5EED0008) (versions (constellation :lock "2026-07"))
  (budget (wall 1h) (mem 32GiB) (qoi-rel-error 1e-2))
  (let wedge (frep (revolve (cheb-profile "wedge.chb"))))
  (let field (flux.free-surface-lbm wedge
               (fluid :model (carreau :mu0 0.12Pa*s :n 0.8) :sigma 0.061N/m)
               (schedule :rate 0.5L/s :tilt (ramp 0deg 65deg 3s))))
  (let J (min (perturbation-growth field :at lip :modes (1 .. 4))))
  (ascent.optimize J :over wedge :method (lbfgs :m 7)
    :until (any (grad-norm 1e-4) (budget-exhausted))
    :emit (ledger report)))"#;

/// The ill-posed variant: the same study demanding wall and memory far
/// beyond the session capability token (capability infeasibility — the
/// check that needs no cost model).
const WEDGE_ILL: &str = r#"(study "cht-wedge-illposed"
  (seed 0x5EED0009) (versions (constellation :lock "2026-07"))
  (capability :wall 100h :mem 512GiB)
  (budget (wall 1h) (mem 32GiB) (qoi-rel-error 1e-2))
  (let wedge (frep (revolve (cheb-profile "wedge.chb"))))
  (let field (flux.free-surface-lbm wedge
               (fluid :model (carreau :mu0 0.12Pa*s :n 0.8) :sigma 0.061N/m)
               (schedule :rate 0.5L/s :tilt (ramp 0deg 65deg 3s))))
  (let J (min (perturbation-growth field :at lip :modes (1 .. 4))))
  (ascent.optimize J :over wedge :method (lbfgs :m 7)
    :until (any (grad-norm 1e-4) (budget-exhausted))
    :emit (ledger report)))"#;

fn admission_cx() -> admission::AdmissionContext<'static> {
    let token = admission::SessionCapability {
        ops: vec![
            "flux.*".to_string(),
            "ascent.*".to_string(),
            "frep".to_string(),
            "xform.*".to_string(),
        ],
        cores: 32.0,
        mem_bytes: 64.0 * 1024.0 * 1024.0 * 1024.0,
        wall_s: 7200.0,
    };
    admission::AdmissionContext {
        router: None,
        chart_requirements: Vec::new(),
        cost_models: BTreeMap::new(),
        capability: Some(token),
        regime: None,
        regime_policy: admission::RegimePolicy::Warn,
    }
}

fn steep_family() -> ProblemFamily {
    let mut c = vec![0.0; 11];
    c[1] = 0.2;
    c[2] = -0.2;
    c[9] = 1.0;
    c[10] = -1.0;
    ProblemFamily {
        base: Poly(c),
        kernel: "cht-wedge-acceptance".to_string(),
    }
}

const RUNGS: [usize; 4] = [12, 24, 48, 96];

#[test]
fn ac_001_admission_teaches_in_milliseconds() {
    let cx = admission_cx();
    let node = sexpr::parse(WEDGE).expect("well-formed study");
    let start = std::time::Instant::now();
    let report = admission::admit(&node, &cx);
    let ok_ms = start.elapsed().as_millis();
    assert!(report.admitted, "the well-posed wedge query is admitted");
    // The ILL-POSED variant is refused fast, with ranked teaching fixes.
    let bad = sexpr::parse(WEDGE_ILL).expect("parses");
    let start = std::time::Instant::now();
    let refusal = admission::admit(&bad, &cx);
    let bad_ms = start.elapsed().as_millis();
    assert!(!refusal.admitted, "the infeasible budget is refused");
    let has_fix = refusal.findings.iter().any(|f| !f.fixes.is_empty());
    assert!(has_fix, "the refusal carries ranked fixes (teaching)");
    println!(
        "{{\"metric\":\"admission\",\"ok_ms\":{ok_ms},\"refusal_ms\":{bad_ms},\
         \"findings\":{}}}",
        refusal.findings.len()
    );
    assert!(bad_ms < 100, "refusal in milliseconds: {bad_ms} ms");
    verdict(
        "ac-001",
        "the wedge query admits; the infeasible variant refuses in milliseconds with \
         ranked teaching fixes",
    );
}

#[test]
fn ac_002_flywheel_discharge_and_anytime_answer() {
    use fs_ir::anytime::run_anytime;
    let family = steep_family();
    let tol = 6e-3;
    let ladder = [30.0, 90.0, 400.0];
    let mut cache = MemCache::default();
    let mut costs = CostTable::new(200.0);
    let report = run_anytime(&family, 1.0, tol, &ladder, &RUNGS, &mut cache, &mut costs);
    // ANYTIME: an immediate colored interval that tightens.
    assert!(!report.trajectory.is_empty(), "an immediate answer exists");
    for step in &report.trajectory {
        assert!(
            matches!(step.color, Color::Verified { .. }),
            "every step is a CERTIFIED interval"
        );
    }
    for w in report.trajectory.windows(2) {
        assert!(w[1].bound <= w[0].bound + 1e-12, "monotone tightening");
    }
    // The flywheel reuse: the repeat query discharges from cache at the
    // smallest budget (the cheap-query loop).
    let again = run_anytime(&family, 1.0, tol, &[5.0], &RUNGS, &mut cache, &mut costs);
    assert!(
        again.refusal.is_none() && again.trajectory.last().expect("step").discharged,
        "the repeat query is a cache hit within a 5-cell budget"
    );
    // TEACHING REFUSAL: an impossible tolerance returns the achieved
    // interval, the priced gap, and the no-point-estimate clause.
    let mut cache2 = MemCache::default();
    let mut costs2 = CostTable::new(200.0);
    let refused = run_anytime(
        &family,
        1.0,
        1e-9,
        &[60.0],
        &RUNGS,
        &mut cache2,
        &mut costs2,
    );
    let note = refused.refusal.expect("the refusal note");
    assert!(
        note.contains("achieved a certified") && note.contains("No best-effort point estimate"),
        "the refusal teaches: {note}"
    );
    println!(
        "{{\"metric\":\"anytime\",\"steps\":{},\"final_bound\":{:.3e},\
         \"cache_hit_budget\":5}}",
        report.trajectory.len(),
        report.trajectory.last().expect("step").bound
    );
    verdict(
        "ac-002",
        "immediate certified interval, monotone tightening, cache-hit repeat at a 5-cell \
         budget, and the teaching refusal with the no-point-estimate clause",
    );
}

#[test]
fn ac_003_package_recheck_solver_free_and_voi_hint() {
    use fs_ir::anytime::run_anytime;
    use fs_plan::voi::{LiveDecision, Probe, ProbeKind, UncertaintyNode, rank_purchases};
    // Discharge the query, wrap the answer, and let the STANDALONE
    // checker re-verify it — certificates, composition, and content
    // address — with zero solver dependency.
    let family = steep_family();
    let mut cache = MemCache::default();
    let mut costs = CostTable::new(200.0);
    let report = run_anytime(&family, 1.0, 6e-3, &[400.0], &RUNGS, &mut cache, &mut costs);
    let last = report.trajectory.last().expect("answer");
    let bound = last.bound;
    // The VoI-priced hint (Proposal C riding the answer).
    let margin = move |v: &[f64]| v[0] - 5e-3;
    let decision = LiveDecision {
        margin: &margin,
        arity: 1,
    };
    let nodes = vec![UncertaintyNode {
        name: "qoi-bound".to_string(),
        lo: 0.0,
        hi: bound.max(1e-6) * 2.0,
        nominal: bound,
    }];
    let menu = vec![Probe {
        name: "climb-final-rung".to_string(),
        target: "qoi-bound".to_string(),
        cost: 12.0,
        shrink: 0.25,
        kind: ProbeKind::Computational,
    }];
    let ranked = rank_purchases(&decision, &nodes, &menu, 32);
    let hint = fs_plan::voi::hint_for_query(&ranked);
    // The package: colored claims, signed, Merkle-rooted.
    let pkg = EvidencePackage::new(Provenance::new("acceptance-e2e", "Cargo.lock"))
        .with_claim({
            let Color::Verified { lo, hi } = last.color else {
                panic!("the wedge trajectory ends verified");
            };
            Claim::from_certificate(
                "wedge-qoi-interval",
                format!("certified half-width {bound:.3e} at tol 6e-3"),
                lo,
                hi,
                "fs-wedge/dwr-certifier",
                "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef",
            )
        })
        .with_claim(Claim::estimated("voi-hint", hint.clone(), "voi-myopic", 1.0))
        .signed("acceptance-gate");
    // THIRD-PARTY RE-VERIFICATION: fs-checker has no solver deps —
    // it checks certificates, composition, signature, and the root.
    let check = fs_checker::check(&pkg);
    assert!(check.passed(), "the package re-verifies solver-free");
    let root = pkg.merkle_root();
    assert!(
        fs_checker::check_against_root(&pkg, root).passed(),
        "the content address matches"
    );
    assert!(
        !fs_checker::check_against_root(&pkg, flip(root)).passed(),
        "a tampered root FAILS the third-party check"
    );
    let pie = check.render_pie();
    println!(
        "{{\"metric\":\"package\",\"root\":\"{root}\",\"hint\":\"{hint}\",\
         \"pie\":\"{}\"}}",
        pie.replace('"', "'").replace('\n', " | ")
    );
    verdict(
        "ac-003",
        "the colored answer ships as a signed Merkle-rooted package carrying its \
         VoI-priced hint; the standalone checker re-verifies it solver-free and \
         catches a tampered root",
    );
}

#[test]
fn ac_004_g5_whole_path_replay() {
    use fs_ir::anytime::run_anytime;
    let run = || -> (Vec<u64>, fs_package::ContentHash) {
        let family = steep_family();
        let mut cache = MemCache::default();
        let mut costs = CostTable::new(200.0);
        let report = run_anytime(
            &family,
            1.0,
            6e-3,
            &[30.0, 400.0],
            &RUNGS,
            &mut cache,
            &mut costs,
        );
        let bits: Vec<u64> = report
            .trajectory
            .iter()
            .map(|s| s.bound.to_bits())
            .collect();
        let pkg = EvidencePackage::new(Provenance::new("acceptance-e2e", "Cargo.lock"))
            .with_claim({
                let Color::Verified { lo, hi } =
                    report.trajectory.last().expect("step").color.clone()
                else {
                    panic!("the replay trajectory ends verified");
                };
                Claim::from_certificate(
                    "wedge-qoi-interval",
                    "replay claim",
                    lo,
                    hi,
                    "fs-wedge/dwr-certifier",
                    "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef",
                )
            })
            .signed("acceptance-gate");
        (bits, pkg.merkle_root())
    };
    let (bits_a, root_a) = run();
    let (bits_b, root_b) = run();
    assert_eq!(bits_a, bits_b, "the trajectory replays bit-exact");
    assert_eq!(root_a, root_b, "the artifact hash replays exactly");
    verdict(
        "ac-004",
        "the whole path — discharge, trajectory, package root — replays bit-exact (G5)",
    );
}

#[test]
fn ac_005_laundering_invariant_across_the_path() {
    // An ESTIMATED intermediate anywhere in the composition can never
    // surface as VERIFIED — checked at the color algebra AND at the
    // package layer.
    let estimated = Color::Estimated {
        estimator: "dwr-guess".to_string(),
        dispersion: 0.1,
    };
    let verified = Color::Verified { lo: 1.0, hi: 1.1 };
    let composed = compose(&verified, &estimated, IntervalOp::Add);
    assert!(
        !matches!(composed, Color::Verified { .. }),
        "weakest-input: verified x estimated is NOT verified: {composed:?}"
    );
    // At the package layer the breakdown keeps them apart — an audit
    // sees exactly how much of the answer is estimated.
    let pkg = EvidencePackage::new(Provenance::new("acceptance-e2e", "Cargo.lock"))
        .with_claim(Claim::from_certificate(
            "hard",
            "certified part",
            1.0,
            1.1,
            "test-solver/cert", "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef",
        ))
        .with_claim(Claim::estimated("soft", "estimated part", "dwr-guess", 0.1))
        .signed("acceptance-gate");
    let breakdown = pkg.color_breakdown();
    assert!(
        breakdown.verified == 1 && breakdown.estimated == 1,
        "the package cannot blur colors: {breakdown:?}"
    );
    verdict(
        "ac-005",
        "estimated inputs never launder to verified — enforced by the compose algebra \
         and visible in the package breakdown",
    );
}
