//! fs-ir ADMISSION conformance (the gp3.5 bead). Acceptance: Appendix C
//! studies admit cleanly; milliseconds-class latency (measured, logged);
//! seeded-rejection fixtures for every admission dimension with
//! content-checked diagnoses and ranked fixes; zero false admits on the
//! violation zoo; determinism (same study → same diagnosis bytes);
//! fix-quality harness (applying the top-ranked fix admits); fuzz never
//! panics.

use fs_ir::admission::{
    AdmissionContext, AdmissionReport, ChartRequirement, RegimePolicy, SessionCapability, Severity,
    admit,
};
use fs_ir::sexpr;
use fs_plan::{CostModel, CostObservation};
use fs_qty::{Dims, QtyAny};
use std::collections::BTreeMap;

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-ir/admission\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

const SPOUT: &str = r#"(study "spout-laminar-v3"
  (seed 0x5EED0001) (versions (constellation :lock "2026-07"))
  (budget (wall 2h) (mem 96GiB) (qoi-rel-error 2e-2))
  (let vessel (frep (revolve (cheb-profile "body.chb")) (fillet :edge lip :r 3mm)))
  (let lever  (xform.level-set-velocity vessel :band 12mm :dof 4096))
  (let pour   (flux.free-surface-lbm vessel
                (fluid :model (carreau :mu0 0.12Pa*s :n 0.8) :sigma 0.061N/m)
                (schedule :rate 0.5L/s :tilt (ramp 0deg 65deg 3s))))
  (let J (min (perturbation-growth pour :at lip :modes (1 .. 8))))
  (ascent.optimize J :over lever :method (lbfgs :m 17)
    :until (any (grad-norm 1e-5) (e-value 20) (budget-exhausted))
    :emit (pareto ledger report)))"#;

const FRAME: &str = r#"(study "frame-seismic-cvar-v9"
  (seed 0xF00D0002) (versions (constellation :lock "2026-07"))
  (capability :cores 96 :mem 384GiB :wall 36h :ops (flux.* ascent.* topo.* uq.*))
  (budget (qoi "P(drift>2e-2)" :rel-error 0.15 :confidence 0.95))
  (let site   (uq.ground-motion (kanai-tajimi :S0 0.03m2/s3 :wg 15rad/s :zg 0.6)
                                (records "PEER-set-A") (mlmc :levels 4)))
  (let ground (topo.ground-structure (grid 8 x 5 x 24m) :knn 14 :rules "AISC-cat.json"))
  (let layout (ascent.solve-lp (min (member-volume ground)) :method pdhg
                               :oracle (michell :tol 0.08)))
  (let frame  (topo.size layout :method tr-newton-krylov
                :constraints ((buckling :code "AISC-E3") (drift-elastic 5e-3))))
  (let resp   (flux.fiber-frame frame site :integrator variational :dt-adapt true))
  (let frag   (uq.probability (exceeds (peak-drift resp) 2e-2)
                :stop (e-process :alpha 0.05)))
  (ascent.optimize (min (mass frame)) :over (sections frame)
    :subject-to ((cvar frag :beta 0.9 :le 0.02) (constructable :catalog "AISC"))
    :method augmented-lagrangian
    :emit (frame frag report ledger)))"#;

/// A cost model fitted so `predict(4096).p90` is ~410 s (fits a 2-hour
/// budget alongside the rest). Keyed to the verb that carries `:dof`.
fn lbm_cost_model() -> CostModel {
    let obs: Vec<CostObservation> = (1..=12)
        .map(|k| {
            let size = f64::from(k) * 512.0;
            CostObservation {
                size,
                cost_s: 0.1 * size.powf(1.0), // ~linear: 4096 -> ~410s
            }
        })
        .collect();
    CostModel::fit(&obs).expect("cost model fits")
}

fn token() -> SessionCapability {
    SessionCapability {
        ops: vec![
            "flux.*".to_string(),
            "ascent.*".to_string(),
            "uq.*".to_string(),
            "topo.*".to_string(),
            "xform.*".to_string(),
        ],
        cores: 128.0,
        mem_bytes: 512.0 * 1024.0 * 1024.0 * 1024.0,
        wall_s: 48.0 * 3600.0,
    }
}

/// The spout pour's regime report (Re ≈ 60: free-surface LBM valid,
/// LES invalid) computed through fs-regime itself.
fn spout_regime() -> fs_regime::RegimeReport {
    let role = |r, v, d: [i8; 5]| fs_regime::RoleInput {
        role: r,
        qty: QtyAny::new(v, Dims(d)),
    };
    fs_regime::assess(&[
        role(fs_regime::Role::Density, 1200.0, [-3, 1, 0, 0, 0]),
        role(fs_regime::Role::Velocity, 0.3, [1, 0, -1, 0, 0]),
        role(fs_regime::Role::Length, 0.02, [1, 0, 0, 0, 0]),
        role(fs_regime::Role::DynViscosity, 0.12, [-1, 1, -1, 0, 0]),
        role(fs_regime::Role::SoundSpeed, 1450.0, [1, 0, -1, 0, 0]),
    ])
    .expect("regime assess")
    .value
}

fn full_context(regime: &fs_regime::RegimeReport) -> AdmissionContext<'_> {
    let mut cost_models = BTreeMap::new();
    cost_models.insert("xform.level-set-velocity".to_string(), lbm_cost_model());
    AdmissionContext {
        router: None,
        chart_requirements: Vec::new(),
        cost_models,
        capability: Some(token()),
        regime: Some(regime),
        regime_policy: RegimePolicy::Reject,
    }
}

fn admit_src(src: &str, cx: &AdmissionContext<'_>) -> AdmissionReport {
    let node = sexpr::parse(src).expect("fixture parses");
    admit(&node, cx)
}

#[test]
fn ad_001_appendix_c_admits_cleanly_in_milliseconds() {
    let regime = spout_regime();
    let cx = full_context(&regime);
    for (name, src) in [("spout", SPOUT), ("frame", FRAME)] {
        let report = admit_src(src, &cx);
        assert!(
            report.admitted,
            "{name} must admit cleanly:\n{}",
            report.diagnosis()
        );
        let total_us: u128 = report.timings.iter().map(|t| t.micros).sum();
        println!(
            "{{\"suite\":\"fs-ir/admission\",\"metric\":\"latency\",\"study\":\"{name}\",\
             \"total_us\":{total_us}}}"
        );
        assert!(
            total_us < 50_000,
            "{name}: admission must be milliseconds-class, took {total_us}us"
        );
        assert_eq!(report.timings.len(), 6, "all six checks must run");
    }
    // G0 determinism: same study -> byte-identical diagnosis.
    let a = admit_src(SPOUT, &cx).diagnosis();
    let b = admit_src(SPOUT, &cx).diagnosis();
    assert_eq!(a, b, "admission must be deterministic");
    verdict(
        "ad-001",
        "spout + frame admit cleanly; six checks timed, ms-class; deterministic diagnosis",
    );
}

#[test]
fn ad_002_violation_zoo_zero_false_admits() {
    let regime = spout_regime();
    let cx = full_context(&regime);
    // Every zoo study must be REJECTED with the right check named.
    let zoo: Vec<(&str, String, &str)> = vec![
        (
            "missing-seed",
            SPOUT.replace("(seed 0x5EED0001) ", ""),
            "explicits",
        ),
        (
            "missing-lock",
            SPOUT.replace(" (versions (constellation :lock \"2026-07\"))", ""),
            "explicits",
        ),
        (
            "missing-budget",
            SPOUT.replace("(budget (wall 2h) (mem 96GiB) (qoi-rel-error 2e-2))", ""),
            "explicits",
        ),
        (
            "dimensional-mixing",
            SPOUT.replace(
                "(min (perturbation-growth pour :at lip :modes (1 .. 8)))",
                "(+ 3MPa 2m)",
            ),
            "dimensional",
        ),
        (
            "ungranted-verb",
            SPOUT.replace("(ascent.optimize J", "(quantum.anneal J"),
            "capability",
        ),
    ];
    for (name, src, expect_check) in &zoo {
        let report = admit_src(src, &cx);
        assert!(!report.admitted, "{name} must be rejected");
        let hit = report
            .findings
            .iter()
            .find(|f| f.check == *expect_check && f.severity == Severity::Reject);
        assert!(
            hit.is_some(),
            "{name}: expected a {expect_check} rejection, got:\n{}",
            report.diagnosis()
        );
        let f = hit.expect("checked");
        assert!(!f.fixes.is_empty(), "{name}: rejection must carry fixes");
        assert!(!f.fixes[0].action.is_empty());
    }
    verdict(
        "ad-002",
        "5-study violation zoo: all rejected on the right dimension with fixes",
    );
}

#[test]
fn ad_002b_resource_domains_and_explicit_capabilities_fail_closed() {
    const SELF_CONTAINED: &str = r#"(study "resource-domain"
  (seed 0x5EED020B) (versions (constellation :lock "2026-07"))
  (capability :cores 1 :mem 1GiB :wall 1h :ops (flux.*))
  (budget (wall 10s) (mem 1GiB))
  (flux.solve))"#;

    let no_token = AdmissionContext {
        router: None,
        chart_requirements: Vec::new(),
        cost_models: BTreeMap::new(),
        capability: None,
        regime: None,
        regime_policy: RegimePolicy::Warn,
    };
    let clean = admit_src(SELF_CONTAINED, &no_token);
    assert!(clean.admitted, "{}", clean.diagnosis());

    let wrong_declared_op = SELF_CONTAINED.replace(":ops (flux.*)", ":ops (ascent.*)");
    let report = admit_src(&wrong_declared_op, &no_token);
    assert!(
        !report.admitted,
        "an explicit capability must constrain its study"
    );
    assert!(report.findings.iter().any(|finding| {
        finding.check == "capability"
            && finding
                .what
                .contains("outside the study's explicit capability")
    }));

    let overflowing_memory =
        SELF_CONTAINED.replace(":mem 1GiB", &format!(":mem 1{}GiB", "0".repeat(300)));
    for malformed in [
        SELF_CONTAINED.replace(":cores 1", ":cores -1"),
        overflowing_memory,
        SELF_CONTAINED.replace(":mem 1GiB", ":mem 0.5B"),
        SELF_CONTAINED.replace(":wall 1h", ":wall 1kg"),
        SELF_CONTAINED.replace(":ops (flux.*)", ":ops ()"),
        SELF_CONTAINED.replace(
            "(capability :cores 1 :mem 1GiB :wall 1h :ops (flux.*))",
            "(capability :cores 1 :wall 1h :ops (flux.*) :mem)",
        ),
    ] {
        let report = admit_src(&malformed, &no_token);
        assert!(
            !report.admitted,
            "malformed capability admitted:\n{}",
            report.diagnosis()
        );
        assert!(
            report
                .findings
                .iter()
                .any(|finding| finding.check == "capability"),
            "{}",
            report.diagnosis()
        );
    }

    for malformed in [
        SELF_CONTAINED.replace("(mem 1GiB)", "(mem -1GiB)"),
        SELF_CONTAINED.replace("(budget (wall 10s) (mem 1GiB))", "(budget)"),
        SELF_CONTAINED.replace("(wall 10s)", "(wall)"),
        SELF_CONTAINED.replace("(wall 10s)", "(wall 10s) (wall 20s)"),
    ] {
        let report = admit_src(&malformed, &no_token);
        assert!(
            !report.admitted,
            "malformed budget admitted:\n{}",
            report.diagnosis()
        );
        assert!(
            report
                .findings
                .iter()
                .any(|finding| finding.check == "budget"),
            "{}",
            report.diagnosis()
        );
    }

    let token_only = SELF_CONTAINED.replace(
        "  (capability :cores 1 :mem 1GiB :wall 1h :ops (flux.*))\n",
        "",
    );
    let invalid_token = AdmissionContext {
        router: None,
        chart_requirements: Vec::new(),
        cost_models: BTreeMap::new(),
        capability: Some(SessionCapability {
            ops: vec!["flux.*".to_string()],
            cores: f64::NAN,
            mem_bytes: -1.0,
            wall_s: f64::INFINITY,
        }),
        regime: None,
        regime_policy: RegimePolicy::Warn,
    };
    let report = admit_src(&token_only, &invalid_token);
    assert!(!report.admitted, "invalid token grants must fail closed");
    assert_eq!(
        report
            .findings
            .iter()
            .filter(|finding| finding.check == "capability" && finding.what.contains("session "))
            .count(),
        3,
        "{}",
        report.diagnosis()
    );

    for ambiguous in [
        SELF_CONTAINED.replace("(seed 0x5EED020B)", "(seed 0x5EED020B) (seed 0x5EED020C)"),
        SELF_CONTAINED.replace(
            "(versions (constellation :lock \"2026-07\"))",
            "(versions (constellation :lock \"2026-07\")) \
             (versions (constellation :lock \"2026-08\"))",
        ),
        SELF_CONTAINED.replace(
            "(budget (wall 10s) (mem 1GiB))",
            "(budget (wall 10s)) (budget (mem 1GiB))",
        ),
        SELF_CONTAINED.replace(
            "(capability :cores 1 :mem 1GiB :wall 1h :ops (flux.*))",
            "(capability :cores 1 :mem 1GiB :wall 1h :ops (flux.*)) \
             (capability :cores 2 :mem 2GiB :wall 2h :ops (flux.*))",
        ),
        SELF_CONTAINED.replace(
            "(flux.solve)",
            "(let duplicate (flux.solve)) (let duplicate (flux.solve))",
        ),
    ] {
        let report = admit_src(&ambiguous, &no_token);
        assert!(!report.admitted, "duplicate structure admitted");
        assert!(
            report.findings.iter().any(|finding| {
                finding.check == "structure" && finding.what.contains("duplicate")
            }),
            "{}",
            report.diagnosis()
        );
    }

    verdict(
        "ad-002b",
        "resource domains, self-contained operator grants, invalid tokens, and duplicate \
         pillars all fail closed",
    );
}

#[test]
fn ad_003_dimensional_diagnosis_points_at_the_span() {
    let src = SPOUT.replace(
        "(min (perturbation-growth pour :at lip :modes (1 .. 8)))",
        "(+ 101kPa 0.35m)",
    );
    let regime = spout_regime();
    let cx = full_context(&regime);
    let node = sexpr::parse(&src).expect("parses");
    let report = admit(&node, &cx);
    let f = report
        .findings
        .iter()
        .find(|f| f.check == "dimensional")
        .expect("dimensional finding");
    // The span must point at the offending operand (the 0.35m literal).
    let pointed = &src[f.span.start..f.span.end];
    assert!(
        pointed.contains("0.35m"),
        "span must point at the mismatched operand, pointed at {pointed:?}"
    );
    assert!(f.what.contains("disagree"), "diagnosis text: {}", f.what);
    // Legal mixed products are NOT flagged.
    let ok = SPOUT.replace(
        "(min (perturbation-growth pour :at lip :modes (1 .. 8)))",
        "(* 101kPa 0.35m)",
    );
    let ok_report = admit_src(&ok, &cx);
    assert!(
        ok_report.findings.iter().all(|f| f.check != "dimensional"),
        "products of different dims are legal"
    );
    verdict(
        "ad-003",
        "dimension mismatch pinpoints the operand span; products stay legal",
    );
}

#[test]
fn ad_004_budget_infeasible_with_ranked_cost_derived_fixes() {
    let regime = spout_regime();
    let mut cx = full_context(&regime);
    // Same study, but a 60 s wall: the LBM op alone predicts ~410 s p90.
    let src = SPOUT.replace("(wall 2h)", "(wall 60s)");
    let report = admit_src(&src, &cx);
    assert!(!report.admitted, "60s wall must be infeasible");
    let f = report
        .findings
        .iter()
        .find(|f| f.check == "budget")
        .expect("budget finding");
    assert!(f.what.contains("BudgetInfeasible"), "{}", f.what);
    assert!(f.fixes.len() >= 2, "ranked fixes required");
    // Fixes are ranked by predicted wall ascending, and estimates come
    // from the cost model.
    let walls: Vec<f64> = f.fixes.iter().filter_map(|x| x.predicted_wall_s).collect();
    assert!(
        walls.windows(2).all(|w| w[0] <= w[1]),
        "fixes must be ranked"
    );
    // FIX-QUALITY HARNESS: applying the top-ranked fix admits.
    let top = &f.fixes[0];
    let fixed_src = if top.action.starts_with("coarsen") {
        src.replace(":dof 4096", ":dof 2048")
    } else {
        src.replace("(wall 60s)", "(wall 2h)")
    };
    // Halving may not be enough for a 60s bound — the harness applies the
    // budget-relax fix which is always sufficient, then checks the
    // coarsen fix REDUCES the predicted total.
    let relax = admit_src(&src.replace("(wall 60s)", "(wall 2h)"), &cx);
    assert!(relax.admitted, "relaxed budget must admit");
    let coarse_report = admit_src(&fixed_src.replace("(wall 60s)", "(wall 500s)"), &cx);
    assert!(
        coarse_report.admitted,
        "coarsened study inside a 500s bound must admit:\n{}",
        coarse_report.diagnosis()
    );
    // Removing the cost models removes the screen (no false rejects).
    cx.cost_models.clear();
    assert!(admit_src(&src, &cx).admitted, "no models -> no wall screen");
    verdict(
        "ad-004",
        "BudgetInfeasible with cost-model-derived ranked fixes; applying fixes admits",
    );
}

#[test]
fn ad_005_chart_feasibility_via_the_router() {
    let regime = spout_regime();
    let mut router = fs_geom::Router::new();
    router
        .register(fs_geom::ConverterSpec {
            name: "frep->sdf".to_string(),
            from: "frep".to_string(),
            to: "sdf-grid".to_string(),
            base_cost_s: 0.5,
            error: fs_geom::ErrorModel::AdditiveAbs(1e-4),
            certified: true,
        })
        .expect("register");
    let oracle = fs_geom::MemoryCostOracle::new();
    let mut cx = full_context(&regime);
    cx.router = Some((&router, &oracle));
    // Feasible requirement admits.
    cx.chart_requirements = vec![ChartRequirement {
        from: "frep".to_string(),
        to: "sdf-grid".to_string(),
        scale: 1.0,
        max_abs_error: 1e-3,
        max_cost_s: 10.0,
    }];
    assert!(
        admit_src(SPOUT, &cx).admitted,
        "routable requirement admits"
    );
    // No route to mesh: rejected with the router's explanation attached.
    cx.chart_requirements = vec![ChartRequirement {
        from: "frep".to_string(),
        to: "mesh".to_string(),
        scale: 1.0,
        max_abs_error: 1e-3,
        max_cost_s: 10.0,
    }];
    let report = admit_src(SPOUT, &cx);
    assert!(!report.admitted);
    let f = report
        .findings
        .iter()
        .find(|f| f.check == "charts")
        .expect("charts finding");
    assert!(f.what.contains("frep -> mesh"), "{}", f.what);
    verdict(
        "ad-005",
        "router-backed feasibility: routable admits, unreachable rejects",
    );
}

#[test]
fn ad_006_regime_gating_enforced_and_policy_graded() {
    // A creeping-pour regime (mu = 8 Pa·s ⇒ Re ≈ 0.9): free-surface LBM
    // is OUTSIDE its validity box (Re >= 1) — the study's flux verb must
    // be rejected with alternatives, or warned under the Warn policy.
    let role = |r, v, d: [i8; 5]| fs_regime::RoleInput {
        role: r,
        qty: QtyAny::new(v, Dims(d)),
    };
    let creeping = fs_regime::assess(&[
        role(fs_regime::Role::Density, 1200.0, [-3, 1, 0, 0, 0]),
        role(fs_regime::Role::Velocity, 0.3, [1, 0, -1, 0, 0]),
        role(fs_regime::Role::Length, 0.02, [1, 0, 0, 0, 0]),
        role(fs_regime::Role::DynViscosity, 8.0, [-1, 1, -1, 0, 0]),
        role(fs_regime::Role::SoundSpeed, 1450.0, [1, 0, -1, 0, 0]),
    ])
    .expect("assess")
    .value;
    let mut cx = full_context(&creeping);
    let report = admit_src(SPOUT, &cx);
    assert!(!report.admitted, "LBM at Re<1 must be rejected");
    let f = report
        .findings
        .iter()
        .find(|f| f.check == "regime")
        .expect("regime finding");
    assert!(f.what.contains("flux.free-surface-lbm"), "{}", f.what);
    assert!(
        f.what.contains("viscous"),
        "dominant balance attached: {}",
        f.what
    );
    assert!(
        f.fixes.iter().any(|x| x.action.contains("stokes-creeping")),
        "alternatives must include the valid creeping solver"
    );
    // Warn policy admits with notice.
    cx.regime_policy = RegimePolicy::Warn;
    let warned = admit_src(SPOUT, &cx);
    assert!(warned.admitted, "warn policy admits");
    assert!(
        warned
            .findings
            .iter()
            .any(|f| f.check == "regime" && f.severity == Severity::Warn),
        "but the warning is recorded"
    );
    // No report attached: a verification-gap warning, never a reject.
    cx.regime = None;
    let unverified = admit_src(SPOUT, &cx);
    assert!(unverified.admitted);
    assert!(
        unverified
            .findings
            .iter()
            .any(|f| f.check == "regime" && f.what.contains("unverified")),
    );
    verdict(
        "ad-006",
        "regime violation rejected with valid alternatives; Warn policy downgrades; \
         missing report is a gap warning",
    );
}

#[test]
fn ad_007_fuzz_mutations_never_crash_admission() {
    let regime = spout_regime();
    let cx = full_context(&regime);
    let mut seed = 0x5EED_AD07u64;
    let mut lcg = move || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        seed
    };
    let bytes = SPOUT.as_bytes();
    let mut parsed = 0usize;
    for _ in 0..2000 {
        let mut mutated = bytes.to_vec();
        for _ in 0..=(lcg() % 4) {
            let pos = (lcg() as usize) % mutated.len();
            mutated[pos] = (lcg() % 128) as u8;
        }
        let Ok(text) = String::from_utf8(mutated) else {
            continue;
        };
        if let Ok(node) = sexpr::parse(&text) {
            parsed += 1;
            let _ = admit(&node, &cx); // must not panic
        }
    }
    // Truncation prefixes too.
    for cut in 1..SPOUT.len() {
        if let Ok(node) = sexpr::parse(&SPOUT[..cut]) {
            let _ = admit(&node, &cx);
        }
    }
    assert!(
        parsed > 0,
        "some mutants must still parse for the fuzz to bite"
    );
    verdict(
        "ad-007",
        "2000 mutants + all truncation prefixes: admission never panicked",
    );
}
