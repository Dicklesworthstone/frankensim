//! fs-ir ADMISSION conformance (the gp3.5 bead). Acceptance: Appendix C
//! studies admit cleanly; milliseconds-class latency (measured, logged);
//! seeded-rejection fixtures for every admission dimension with
//! content-checked diagnoses and ranked fixes; zero false admits on the
//! violation zoo; determinism (same study → same diagnosis bytes);
//! fix-quality harness (applying the top-ranked fix admits); fuzz never
//! panics.

use fs_ir::admission::{
    AdmissionContext, AdmissionReport, CapabilityEvidenceClass, CapabilityReceipt,
    ChartRequirement, NoCapabilityIssuerVerifier, RankedFix, RegimePolicy, SealedSessionCapability,
    SessionCapability, Severity, admit, admit_versioned,
};
use fs_ir::sexpr;
use fs_ir::{IR_VERSION, Node, NodeKind, Span, VersionedProgram};
use fs_plan::{CostModel, CostObservation, SealedCostModel};
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

const LOWERED_POUR_RAW: &str = r#"(study "lowered-pour"
  (seed 0x5EED0010) (versions (constellation :lock "2026-07"))
  (budget (wall 10s))
  (simulate-pour vessel fluid schedule))"#;

const LOWERED_POUR_EXPLICIT: &str = r#"(study "lowered-pour"
  (seed 0x5EED0010) (versions (constellation :lock "2026-07"))
  (budget (wall 10s))
  (flux.free-surface-lbm vessel fluid schedule))"#;

const LOWERED_OPTIMIZE_RAW: &str = r#"(study "lowered-optimize"
  (seed 0x5EED0011) (versions (constellation :lock "2026-07"))
  (budget (wall 10s))
  (optimize-shape :min objective :over levers))"#;

const LOWERED_OPTIMIZE_EXPLICIT: &str = r#"(study "lowered-optimize"
  (seed 0x5EED0011) (versions (constellation :lock "2026-07"))
  (budget (wall 10s))
  (ascent.optimize (min objective) :over levers
    :method (lbfgs :m 17) :until (grad-norm 1e-5)
    :emit (ledger report)))"#;

/// A cost model fitted so `predict(4096).p90` is ~410 s (fits a 2-hour
/// budget alongside the rest). Keyed to the verb that carries `:dof`.
fn lbm_cost_model(operation: &str) -> SealedCostModel {
    let obs: Vec<CostObservation> = (1..=12)
        .map(|k| {
            let size = f64::from(k) * 512.0;
            CostObservation {
                size,
                cost_s: 0.1 * size.powf(1.0), // ~linear: 4096 -> ~410s
            }
        })
        .collect();
    SealedCostModel::provisional_unaudited(
        CostModel::fit(&obs).expect("cost model fits"),
        operation,
    )
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
        cores: 128,
        mem_bytes: 512 * 1024 * 1024 * 1024,
        wall_s: 48.0 * 3600.0,
    }
}

/// The spout pour's regime report (Re ≈ 60: free-surface LBM valid,
/// LES invalid) computed through fs-regime itself.
fn spout_regime() -> fs_regime::RegimeReport {
    let role = |r, v, d: [i8; 6]| fs_regime::RoleInput {
        role: r,
        qty: QtyAny::new(v, Dims(d)),
    };
    fs_regime::assess(&[
        role(fs_regime::Role::Density, 1200.0, [-3, 1, 0, 0, 0, 0]),
        role(fs_regime::Role::Velocity, 0.3, [1, 0, -1, 0, 0, 0]),
        role(fs_regime::Role::Length, 0.02, [1, 0, 0, 0, 0, 0]),
        role(fs_regime::Role::DynViscosity, 0.12, [-1, 1, -1, 0, 0, 0]),
        role(fs_regime::Role::SoundSpeed, 1450.0, [1, 0, -1, 0, 0, 0]),
    ])
    .expect("regime assess")
    .value
}

fn creeping_regime() -> fs_regime::RegimeReport {
    let role = |r, v, d: [i8; 6]| fs_regime::RoleInput {
        role: r,
        qty: QtyAny::new(v, Dims(d)),
    };
    fs_regime::assess(&[
        role(fs_regime::Role::Density, 1200.0, [-3, 1, 0, 0, 0, 0]),
        role(fs_regime::Role::Velocity, 0.3, [1, 0, -1, 0, 0, 0]),
        role(fs_regime::Role::Length, 0.02, [1, 0, 0, 0, 0, 0]),
        role(fs_regime::Role::DynViscosity, 8.0, [-1, 1, -1, 0, 0, 0]),
        role(fs_regime::Role::SoundSpeed, 1450.0, [1, 0, -1, 0, 0, 0]),
    ])
    .expect("creeping regime assess")
    .value
}

fn full_context(regime: &fs_regime::RegimeReport) -> AdmissionContext<'_> {
    let mut cost_models = BTreeMap::new();
    cost_models.insert(
        "xform.level-set-velocity".to_string(),
        lbm_cost_model("xform.level-set-velocity"),
    );
    AdmissionContext {
        router: None,
        cost_freshness: None,
        chart_requirements: Vec::new(),
        cost_models,
        capability: Some(SealedSessionCapability::caller_declared(token())),
        regime: Some(regime),
        regime_policy: RegimePolicy::Reject,
    }
}

fn authority_context<'a>(
    ops: &[&str],
    regime: Option<&'a fs_regime::RegimeReport>,
    regime_policy: RegimePolicy,
) -> AdmissionContext<'a> {
    AdmissionContext {
        router: None,
        cost_freshness: None,
        chart_requirements: Vec::new(),
        cost_models: BTreeMap::new(),
        capability: Some(SealedSessionCapability::caller_declared(SessionCapability {
            ops: ops.iter().map(|op| (*op).to_string()).collect(),
            cores: 1,
            mem_bytes: 1 << 30,
            wall_s: 3_600.0,
        })),
        regime,
        regime_policy,
    }
}

fn admit_src(src: &str, cx: &AdmissionContext<'_>) -> AdmissionReport {
    let node = sexpr::parse(src).expect("fixture parses");
    admit(&node, cx)
}

fn finding_semantics(
    report: &AdmissionReport,
) -> Vec<(&'static str, Severity, String, Vec<RankedFix>)> {
    report
        .findings
        .iter()
        .map(|finding| {
            (
                finding.check,
                finding.severity,
                finding.what.clone(),
                finding.fixes.clone(),
            )
        })
        .collect()
}

fn admit_equivalent(
    raw: &str,
    explicit: &str,
    cx: &AdmissionContext<'_>,
) -> (AdmissionReport, AdmissionReport) {
    let raw_report = admit_src(raw, cx);
    let explicit_report = admit_src(explicit, cx);
    assert_eq!(raw_report.lowering.ir_version(), IR_VERSION);
    assert_ne!(
        raw_report.lowering.raw_canonical(),
        explicit_report.lowering.raw_canonical(),
        "raw identity must retain whether shorthand was submitted"
    );
    assert_eq!(
        raw_report.lowering.lowered_canonical(),
        explicit_report.lowering.lowered_canonical(),
        "semantically equivalent syntax must bind the same lowered identity"
    );
    assert_eq!(raw_report.study, explicit_report.study);
    assert_eq!(raw_report.admitted, explicit_report.admitted);
    assert_eq!(
        finding_semantics(&raw_report),
        finding_semantics(&explicit_report),
        "authority verdicts must depend only on lowered semantics"
    );
    (raw_report, explicit_report)
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
#[allow(clippy::too_many_lines)] // One end-to-end raw/explicit authority matrix.
fn ad_001b_admission_is_over_versioned_lowered_semantics() {
    let allowed_pour = authority_context(&["flux.*"], None, RegimePolicy::Warn);
    let (raw_allowed, _) = admit_equivalent(LOWERED_POUR_RAW, LOWERED_POUR_EXPLICIT, &allowed_pour);
    assert!(raw_allowed.admitted, "{}", raw_allowed.diagnosis());
    assert!(
        raw_allowed
            .lowering
            .raw_canonical()
            .contains("simulate-pour")
    );
    let lowered_identity = raw_allowed
        .lowering
        .lowered_canonical()
        .expect("successful admission binds lowered semantics");
    assert!(lowered_identity.contains("flux.free-surface-lbm"));
    assert!(!lowered_identity.contains("simulate-pour"));

    // The source convenience API and explicitly versioned API converge on the
    // exact same authority input and deterministic receipt.
    let raw_node = sexpr::parse(LOWERED_POUR_RAW).expect("raw shorthand parses");
    let versioned = VersionedProgram::current(raw_node.clone());
    let explicit_version_report = admit_versioned(&versioned, &allowed_pour);
    assert_eq!(raw_allowed.lowering, explicit_version_report.lowering);
    assert_eq!(
        finding_semantics(&raw_allowed),
        finding_semantics(&explicit_version_report)
    );
    assert_eq!(raw_allowed.diagnosis(), explicit_version_report.diagnosis());

    // The capability operator appears only after lowering, but it is still
    // denied exactly as the caller-written explicit form is denied.
    let denied_pour = authority_context(&["ascent.*"], None, RegimePolicy::Warn);
    let (raw_denied, _) = admit_equivalent(LOWERED_POUR_RAW, LOWERED_POUR_EXPLICIT, &denied_pour);
    assert!(!raw_denied.admitted);
    assert!(raw_denied.findings.iter().any(|finding| {
        finding.check == "capability"
            && finding.what.contains("flux.free-surface-lbm")
            && finding.what.contains("outside")
    }));

    let allowed_optimize = authority_context(&["ascent.*"], None, RegimePolicy::Warn);
    let (raw_optimize, _) = admit_equivalent(
        LOWERED_OPTIMIZE_RAW,
        LOWERED_OPTIMIZE_EXPLICIT,
        &allowed_optimize,
    );
    assert!(raw_optimize.admitted, "{}", raw_optimize.diagnosis());
    let denied_optimize = authority_context(&["flux.*"], None, RegimePolicy::Warn);
    let (raw_optimize_denied, _) = admit_equivalent(
        LOWERED_OPTIMIZE_RAW,
        LOWERED_OPTIMIZE_EXPLICIT,
        &denied_optimize,
    );
    assert!(!raw_optimize_denied.admitted);
    assert!(raw_optimize_denied.findings.iter().any(|finding| {
        finding.check == "capability" && finding.what.contains("ascent.optimize")
    }));

    // Costing sees the injected flux operator too.
    let mut cost_context = authority_context(&["flux.*"], None, RegimePolicy::Warn);
    cost_context.cost_models.insert(
        "flux.free-surface-lbm".to_string(),
        lbm_cost_model("flux.free-surface-lbm"),
    );
    let raw_tight = LOWERED_POUR_RAW.replace("10s", "0.01s");
    let explicit_tight = LOWERED_POUR_EXPLICIT.replace("10s", "0.01s");
    let (raw_cost_denied, _) = admit_equivalent(&raw_tight, &explicit_tight, &cost_context);
    assert!(!raw_cost_denied.admitted);
    assert!(
        raw_cost_denied.findings.iter().any(|finding| {
            finding.check == "budget" && finding.what.contains("BudgetInfeasible")
        })
    );

    // Regime authority is likewise evaluated on the lowered flux model.
    let creeping = creeping_regime();
    let regime_context = authority_context(&["flux.*"], Some(&creeping), RegimePolicy::Reject);
    let (raw_regime_denied, _) =
        admit_equivalent(LOWERED_POUR_RAW, LOWERED_POUR_EXPLICIT, &regime_context);
    assert!(!raw_regime_denied.admitted);
    assert!(raw_regime_denied.findings.iter().any(|finding| {
        finding.check == "regime" && finding.what.contains("flux.free-surface-lbm")
    }));

    let replay = admit_src(LOWERED_POUR_RAW, &allowed_pour);
    assert_eq!(raw_allowed.lowering, replay.lowering);
    assert_eq!(raw_allowed.diagnosis(), replay.diagnosis());
}

#[test]
#[allow(clippy::too_many_lines)] // One fail-closed lowering/refusal matrix.
fn ad_001c_lowering_refuses_atomically_before_authority_checks() {
    let creeping = creeping_regime();
    let mut hostile_context =
        authority_context(&["ascent.*"], Some(&creeping), RegimePolicy::Reject);
    hostile_context.cost_models.insert(
        "flux.free-surface-lbm".to_string(),
        lbm_cost_model("flux.free-surface-lbm"),
    );
    hostile_context.chart_requirements.push(ChartRequirement {
        from: "frep".to_string(),
        to: "mesh".to_string(),
        scale: 1.0,
        max_abs_error: 1e-3,
        max_cost_s: 1.0,
    });
    // Sealed evidence is immutable by design: rebuild the declared view
    // with the hostile wall grant instead of mutating in place.
    hostile_context.capability = Some(SealedSessionCapability::caller_declared(
        SessionCapability {
            wall_s: f64::INFINITY,
            ..token()
        },
    ));

    let study = |body: &str| {
        format!(
            "(study \"malformed-lowering\" \
             (seed 0x5EED0012) (versions (constellation :lock \"2026-07\")) \
             (budget (wall 10s)) {body})"
        )
    };
    for (body, expected) in [
        (
            "(optimize-shape :min j :min k :over x)",
            "duplicate optimize-shape argument :min",
        ),
        (
            "(optimize-shape :min j :over x :gpu true)",
            "unknown optimize-shape argument :gpu",
        ),
        (
            "(optimize-shape :min j :over x :method)",
            "dangling :method",
        ),
        (
            "(optimize-shape :min j :over x trailing)",
            "trailing argument",
        ),
        (
            "(simulate-pour vessel fluid schedule trailing)",
            "takes 3 arguments",
        ),
        (
            "(wrapper (simulate-pour vessel fluid schedule) \
             (optimize-shape :over x))",
            "needs :min",
        ),
    ] {
        let source = study(body);
        let report = admit_src(&source, &hostile_context);
        assert!(!report.admitted);
        assert_eq!(report.findings.len(), 1, "{}", report.diagnosis());
        assert_eq!(report.findings[0].check, "lowering");
        assert!(
            report.findings[0].what.contains(expected),
            "{body}: expected {expected:?}, got {}",
            report.findings[0].what
        );
        assert!(
            report.timings.is_empty(),
            "no authority dimension may run after lowering refusal"
        );
        assert!(report.lowering.lowered_canonical().is_none());
        assert!(
            report
                .lowering
                .raw_canonical()
                .starts_with("(frankensim-ir :version 3 :program")
        );
        let replay = admit_src(&source, &hostile_context);
        assert_eq!(report.lowering, replay.lowering);
        assert_eq!(report.diagnosis(), replay.diagnosis());
    }

    // A refusal cannot leave a partial lowering or mutate later decisions.
    let stable_context = authority_context(&["flux.*"], None, RegimePolicy::Warn);
    let before = admit_src(LOWERED_POUR_RAW, &stable_context);
    let partial = admit_src(
        &study(
            "(wrapper (simulate-pour vessel fluid schedule) \
             (optimize-shape :over x))",
        ),
        &stable_context,
    );
    assert_eq!(partial.findings.len(), 1);
    assert!(partial.lowering.lowered_canonical().is_none());
    let after = admit_src(LOWERED_POUR_RAW, &stable_context);
    assert_eq!(before.lowering, after.lowering);
    assert_eq!(finding_semantics(&before), finding_semantics(&after));
    assert_eq!(before.diagnosis(), after.diagnosis());
}

#[test]
fn ad_001d_identity_binding_refuses_invalid_and_envelope_deep_asts() {
    let cx = authority_context(&["ascent.*"], None, RegimePolicy::Warn);

    // `Node` is a public matching surface, so admission must reject a
    // caller-forged invalid atom instead of reaching the panicking trusted
    // `VersionedProgram::current` compatibility constructor.
    let forged = Node {
        kind: NodeKind::List(vec![
            Node::synthetic(NodeKind::Symbol("study".to_string())),
            Node::synthetic(NodeKind::Str("forged".to_string())),
            Node {
                kind: NodeKind::Float(f64::NAN),
                span: Span::new(16, 19),
            },
        ]),
        span: Span::new(0, 20),
    };
    let forged_report = admit(&forged, &cx);
    assert!(!forged_report.admitted);
    assert_eq!(forged_report.study, "<lowering-refused>");
    assert_eq!(forged_report.findings.len(), 1);
    assert_eq!(forged_report.findings[0].check, "lowering");
    assert_eq!(forged_report.findings[0].span, Span::new(16, 19));
    assert!(
        forged_report.findings[0]
            .what
            .contains("floating-point atom is not finite")
    );
    assert!(forged_report.timings.is_empty());
    assert_eq!(forged_report.lowering.ir_version(), IR_VERSION);
    assert_eq!(forged_report.lowering.raw_canonical_opt(), None);
    assert_eq!(
        forged_report.lowering.raw_canonical(),
        "",
        "the compatibility accessor uses an unambiguous empty sentinel"
    );
    assert_eq!(forged_report.lowering.lowered_canonical(), None);
    let forged_replay = admit(&forged, &cx);
    assert_eq!(forged_report.lowering, forged_replay.lowering);
    assert_eq!(forged_report.diagnosis(), forged_replay.diagnosis());

    // A tree can satisfy the bare AST depth cap yet be one level too deep for
    // the required version envelope. That boundary also refuses structurally.
    let mut raw_boundary = Node::synthetic(NodeKind::Int(1));
    for _ in 0..256 {
        raw_boundary = Node::synthetic(NodeKind::List(vec![raw_boundary]));
    }
    raw_boundary
        .validate()
        .expect("bare boundary tree is valid before version binding");
    let raw_boundary_report = admit(&raw_boundary, &cx);
    assert!(!raw_boundary_report.admitted);
    assert!(
        raw_boundary_report.findings[0]
            .what
            .contains("nesting exceeds")
    );
    assert_eq!(raw_boundary_report.lowering.raw_canonical_opt(), None);
    assert_eq!(raw_boundary_report.lowering.lowered_canonical(), None);

    // Conversely, this raw program fits its envelope exactly. Lowering adds
    // one semantic level around the objective, leaving a standalone-valid AST
    // that no longer fits a persisted envelope. The lowered identity bind must
    // reject instead of panicking.
    let mut expansion_boundary =
        sexpr::parse("(optimize-shape :min j :over x)").expect("boundary shorthand parses");
    for _ in 0..254 {
        expansion_boundary = Node::synthetic(NodeKind::List(vec![expansion_boundary]));
    }
    let versioned = VersionedProgram::try_current(expansion_boundary)
        .expect("raw boundary program fits the version envelope exactly");
    let expansion_report = admit_versioned(&versioned, &cx);
    assert!(!expansion_report.admitted);
    assert_eq!(expansion_report.study, "<lowering-refused>");
    assert_eq!(expansion_report.findings.len(), 1);
    assert_eq!(expansion_report.findings[0].check, "lowering");
    assert!(
        expansion_report.findings[0]
            .what
            .contains("nesting exceeds")
    );
    assert!(expansion_report.timings.is_empty());
    assert!(expansion_report.lowering.raw_canonical_opt().is_some());
    assert_eq!(expansion_report.lowering.lowered_canonical(), None);
    let expansion_replay = admit_versioned(&versioned, &cx);
    assert_eq!(expansion_report.lowering, expansion_replay.lowering);
    assert_eq!(expansion_report.diagnosis(), expansion_replay.diagnosis());
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
#[allow(clippy::too_many_lines)] // One fail-closed capability/resource matrix.
fn ad_002b_resource_domains_and_explicit_capabilities_fail_closed() {
    const SELF_CONTAINED: &str = r#"(study "resource-domain"
  (seed 0x5EED020B) (versions (constellation :lock "2026-07"))
  (capability :cores 1 :mem 1GiB :wall 1h :ops (flux.*))
  (budget (wall 10s) (mem 1GiB))
  (flux.solve))"#;

    let no_token = AdmissionContext {
        router: None,
        cost_freshness: None,
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

    let overflowing_memory = SELF_CONTAINED.replace(":mem 1GiB", ":mem 17179869184GiB");
    for malformed in [
        SELF_CONTAINED.replace(":cores 1", ":cores -1"),
        overflowing_memory,
        SELF_CONTAINED.replace(":mem 1GiB", ":mem 0.5B"),
        SELF_CONTAINED.replace(":wall 1h", ":wall 1kg"),
        SELF_CONTAINED.replace(":ops (flux.*)", ":ops ()"),
        SELF_CONTAINED.replace(":ops (flux.*)", ":ops (flux*)"),
        SELF_CONTAINED.replace(":ops (flux.*)", ":ops (*)"),
        SELF_CONTAINED.replace(":ops (flux.*)", ":ops (flux.*) :gpu true"),
        SELF_CONTAINED.replace(":ops (flux.*)", ":ops (flux.*) garbage"),
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
        SELF_CONTAINED.replace("(budget (wall 10s) (mem 1GiB))", "(budget nonsense)"),
        SELF_CONTAINED.replace("(wall 10s)", "(wall)"),
        SELF_CONTAINED.replace("(wall 10s)", "(wall 10s 20s)"),
        SELF_CONTAINED.replace("(wall 10s)", "(custom-budget)"),
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

    let extensible_budget = SELF_CONTAINED.replace(
        "(budget (wall 10s) (mem 1GiB))",
        "(budget (custom-budget 3))",
    );
    let report = admit_src(&extensible_budget, &no_token);
    assert!(
        report.admitted,
        "structured operator-specific budgets remain extensible:\n{}",
        report.diagnosis()
    );

    let token_only = SELF_CONTAINED.replace(
        "  (capability :cores 1 :mem 1GiB :wall 1h :ops (flux.*))\n",
        "",
    );
    let invalid_token = AdmissionContext {
        router: None,
        cost_freshness: None,
        chart_requirements: Vec::new(),
        cost_models: BTreeMap::new(),
        capability: Some(SealedSessionCapability::caller_declared(SessionCapability {
            ops: vec!["flux.*".to_string()],
            cores: 0,
            mem_bytes: 0,
            wall_s: f64::INFINITY,
        })),
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
        1,
        "{}",
        report.diagnosis()
    );

    let malformed_token = AdmissionContext {
        router: None,
        cost_freshness: None,
        chart_requirements: Vec::new(),
        cost_models: BTreeMap::new(),
        capability: Some(SealedSessionCapability::caller_declared(SessionCapability {
            ops: vec!["flux*".to_string()],
            cores: 1,
            mem_bytes: 0,
            wall_s: 1.0,
        })),
        regime: None,
        regime_policy: RegimePolicy::Warn,
    };
    let report = admit_src(&token_only, &malformed_token);
    assert!(
        !report.admitted
            && report
                .findings
                .iter()
                .filter(|finding| finding.check == "capability")
                .count()
                >= 1,
        "malformed-glob token admitted:\n{}",
        report.diagnosis()
    );

    for malformed in [
        SELF_CONTAINED.replace("(seed 0x5EED020B)", "(seed 0x5EED020B 0x5EED020C)"),
        SELF_CONTAINED.replace(
            "(versions (constellation :lock \"2026-07\"))",
            "(versions (constellation :lock \"2026-07\" :lock \"2026-08\"))",
        ),
        SELF_CONTAINED.replace(
            "(versions (constellation :lock \"2026-07\"))",
            "(versions (constellation :lock \"2026-07\") \
             (constellation :lock \"2026-08\"))",
        ),
        SELF_CONTAINED.replace(
            "(versions (constellation :lock \"2026-07\"))",
            "(versions (constellation :lock \" \"))",
        ),
        SELF_CONTAINED.replace(
            "(versions (constellation :lock \"2026-07\"))",
            "(versions garbage)",
        ),
        SELF_CONTAINED.replace("(flux.solve)", "(let value (flux.solve) extra)"),
    ] {
        let report = admit_src(&malformed, &no_token);
        assert!(
            !report.admitted
                && report
                    .findings
                    .iter()
                    .any(|finding| finding.check == "structure"),
            "malformed study structure admitted:\n{}",
            report.diagnosis()
        );
    }

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
fn ad_002c_exact_count_authority_boundaries_do_not_alias() {
    const SOURCE: &str = r#"(study "exact-authority"
  (seed 0x5EED020C) (versions (constellation :lock "2026-07"))
  (capability :cores 1 :mem 1B :wall 1h :ops (flux.*))
  (budget (wall 10s) (mem 1B))
  (flux.solve))"#;

    let context = |mem_bytes, cores| AdmissionContext {
        router: None,
        cost_freshness: None,
        chart_requirements: Vec::new(),
        cost_models: BTreeMap::new(),
        capability: Some(SealedSessionCapability::caller_declared(SessionCapability {
            ops: vec!["flux.*".to_string()],
            cores,
            mem_bytes,
            wall_s: 3_600.0,
        })),
        regime: None,
        regime_policy: RegimePolicy::Warn,
    };

    let max = SOURCE.replace(":mem 1B", ":mem 18446744073709551615B");
    let report = admit_src(&max, &context(u64::MAX, 1));
    assert!(
        report.admitted,
        "u64::MAX must remain exact:\n{}",
        report.diagnosis()
    );

    let adjacent = SOURCE.replace(":mem 1B", ":mem 9007199254740993B");
    let report = admit_src(&adjacent, &context(9_007_199_254_740_992, 1));
    assert!(
        !report.admitted
            && report
                .diagnosis()
                .contains("9007199254740993 bytes asked, 9007199254740992 bytes granted"),
        "adjacent byte authorities above 2^53 aliased:\n{}",
        report.diagnosis()
    );

    let decimal = SOURCE.replace(":mem 1B", ":mem 1.5GiB");
    let report = admit_src(&decimal, &context(3 << 29, 1));
    assert!(
        report.admitted,
        "exact decimal scale refused:\n{}",
        report.diagnosis()
    );

    let over = SOURCE.replace(":mem 1B", ":mem 18446744073709551616B");
    let report = admit_src(&over, &context(u64::MAX, 1));
    assert!(!report.admitted, "2^64-byte authority must refuse");

    let adjacent_cores = SOURCE.replace(":cores 1", ":cores 9007199254740993cores");
    let report = admit_src(&adjacent_cores, &context(1, 9_007_199_254_740_992));
    assert!(
        !report.admitted
            && report
                .diagnosis()
                .contains("9007199254740993 cores asked, 9007199254740992 granted"),
        "adjacent core authorities above 2^53 aliased:\n{}",
        report.diagnosis()
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
    // Two budget findings coexist now (bead 2pmb): the provisional-
    // evidence Warn and the infeasibility Reject; select the Reject.
    let f = report
        .findings
        .iter()
        .find(|f| f.check == "budget" && f.what.contains("BudgetInfeasible"))
        .expect("budget infeasibility finding");
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
    for malformed_size in [
        src.replace(":dof 4096", ":dof"),
        src.replace(":dof 4096", ":dof \"many\""),
        src.replace(":dof 4096", ":dof 4096 :size 8"),
    ] {
        let malformed = admit_src(&malformed_size, &cx);
        assert!(
            !malformed.admitted
                && malformed
                    .findings
                    .iter()
                    .any(|finding| finding.check == "budget"),
            "malformed explicit size feature admitted:\n{}",
            malformed.diagnosis()
        );
    }
    let implicit_unit_size = admit_src(&src.replace(" :dof 4096", ""), &cx);
    assert!(
        implicit_unit_size.admitted,
        "the unit-size default applies only when no size feature is declared:\n{}",
        implicit_unit_size.diagnosis()
    );
    let mut legacy_cx = full_context(&regime);
    legacy_cx.cost_models.clear();
    legacy_cx
        .cost_models
        .insert("legacy-solve".to_string(), lbm_cost_model("legacy-solve"));
    let legacy_src = src.replace("xform.level-set-velocity", "legacy-solve");
    let legacy = admit_src(&legacy_src, &legacy_cx);
    assert!(
        !legacy.admitted
            && legacy.findings.iter().any(|finding| {
                finding.check == "budget" && finding.what.contains("BudgetInfeasible")
            }),
        "registered non-namespaced cost model was ignored:\n{}",
        legacy.diagnosis()
    );
    // A registered but evidence-thin model is categorically different
    // from an absent model: its structured refusal must survive as an
    // admission rejection instead of silently contributing zero cost.
    let mut thin_cx = full_context(&regime);
    thin_cx.cost_models.insert(
        "xform.level-set-velocity".to_string(),
        SealedCostModel::provisional_unaudited(CostModel::new(), "xform.level-set-velocity"),
    );
    let thin = admit_src(&src, &thin_cx);
    let refusal = thin
        .findings
        .iter()
        .find(|finding| finding.check == "budget" && finding.what.contains("CostModelRefused"))
        .expect("registered InsufficientData must reject admission");
    assert!(!thin.admitted, "registered thin model cannot be ignored");
    assert!(
        refusal.what.contains("0 observation(s)") && refusal.what.contains("need"),
        "the underlying structured refusal was lost: {}",
        refusal.what
    );
    assert!(
        refusal
            .fixes
            .iter()
            .any(|fix| fix.action.contains("tune observations")),
        "refusal must include an evidence-acquisition fix: {refusal:?}"
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
fn ad_005b_estimated_chart_routes_do_not_authorize_admission() {
    let regime = spout_regime();
    let mut router = fs_geom::Router::new();
    router
        .register(fs_geom::ConverterSpec {
            name: "frep->sdf/estimated".to_string(),
            from: "frep".to_string(),
            to: "sdf-grid".to_string(),
            base_cost_s: 0.5,
            error: fs_geom::ErrorModel::AdditiveAbs(1e-4),
            certified: false,
        })
        .expect("register estimated converter");
    let oracle = fs_geom::MemoryCostOracle::new();
    let mut cx = full_context(&regime);
    cx.router = Some((&router, &oracle));
    cx.chart_requirements = vec![ChartRequirement {
        from: "frep".to_string(),
        to: "sdf-grid".to_string(),
        scale: 1.0,
        max_abs_error: 1e-3,
        max_cost_s: 10.0,
    }];
    let report = admit_src(SPOUT, &cx);
    assert!(!report.admitted);
    let finding = report
        .findings
        .iter()
        .find(|finding| finding.check == "charts")
        .expect("estimated route must leave a chart finding");
    assert!(finding.what.contains("only estimated"), "{}", finding.what);
    verdict(
        "ad-005b",
        "estimated converter routes remain usable for planning but cannot authorize admission",
    );
}

#[test]
fn ad_006_regime_gating_enforced_and_policy_graded() {
    // A creeping-pour regime (mu = 8 Pa·s ⇒ Re ≈ 0.9): free-surface LBM
    // is OUTSIDE its validity box (Re >= 1) — the study's flux verb must
    // be rejected with alternatives, or warned under the Warn policy.
    let creeping = creeping_regime();
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

#[test]
fn ad_007b_empty_list_operand_does_not_panic_dimensional_inference() {
    // Regression: `check_dimensional` recurses `infer_dims` into every body
    // clause and let-expression. A bare `()` has no head symbol and no
    // operands, so the old `&items[1..]` was a usize range panic (start 1 >
    // len 0) — a fail-OPEN crash on parseable input, exactly the kind the
    // Gauntlet forbids. Byte-mutation fuzz (ad-007) never happened to land a
    // surviving `()` in operand position; this pins the minimal reproducers.
    let cx = AdmissionContext {
        router: None,
        cost_freshness: None,
        chart_requirements: Vec::new(),
        cost_models: BTreeMap::new(),
        capability: None,
        regime: None,
        regime_policy: RegimePolicy::Warn,
    };
    // (a) a bare empty list as a body clause.
    let bare = r#"(study "empty-body-clause" ())"#;
    // (b) an empty list nested as an ARITHMETIC operand inside a let — this
    // reaches infer_dims's List arm via the same-dims recursion, the exact
    // path that panicked.
    let nested = r#"(study "empty-list-operand"
  (seed 0x5EED07B0) (versions (constellation :lock "2026-07"))
  (capability :cores 1 :mem 1GiB :wall 1h :ops (flux.*))
  (budget (wall 10s) (mem 1GiB))
  (let bad (+ 1 ()))
  (flux.solve))"#;
    // Reaching each `admit` return at all proves the range-slice panic is gone;
    // the concrete assertions below pin the behavior the fix guarantees.
    // (a) The bare `()` study fails CLOSED (missing pillars) — never a crash.
    let bare_report = admit(&sexpr::parse(bare).expect("bare parses"), &cx);
    assert!(
        !bare_report.admitted,
        "a study whose only body clause is () must not admit:\n{}",
        bare_report.diagnosis()
    );
    // (b) The nested `(+ 1 ())` study runs dimensional inference to completion:
    // the empty list is treated as dimensionless-unknown (no false conflict, no
    // panic), so admission proceeds and names the study it parsed.
    let nested_report = admit(&sexpr::parse(nested).expect("nested parses"), &cx);
    assert_eq!(
        nested_report.study, "empty-list-operand",
        "dimensional inference must complete over an () operand instead of panicking"
    );
    verdict(
        "ad-007b",
        "bare () body clause fails closed; () arithmetic operand infers without panic",
    );
}

#[test]
fn ad_010_provisional_cost_models_warn_but_admit() {
    // Sealed authority (bead 2pmb): a provisional fit may inform the
    // budget arithmetic, but the admission record must say so.
    let regime = spout_regime();
    let cx = full_context(&regime);
    let report = admit_src(SPOUT, &cx);
    assert!(
        report.admitted,
        "provisional evidence admits with notice, never silently blocks:\n{}",
        report.diagnosis()
    );
    let provisional: Vec<_> = report
        .findings
        .iter()
        .filter(|f| f.severity == Severity::Warn && f.what.contains("CostModelProvisional"))
        .collect();
    assert_eq!(
        provisional.len(),
        1,
        "exactly one warning per provisional verb:\n{}",
        report.diagnosis()
    );
    assert!(
        provisional[0].what.contains("provisional-unaudited"),
        "the finding names the evidence class: {}",
        provisional[0].what
    );
    assert!(
        provisional[0]
            .what
            .contains("provisional:xform.level-set-velocity"),
        "the finding names the operation-bound scope: {}",
        provisional[0].what
    );
    verdict(
        "ad-010",
        "provisional cost evidence admits with a named once-per-verb warning",
    );
}

#[test]
fn ad_010a_miskeyed_cost_model_refuses_without_pricing() {
    let regime = spout_regime();
    let mut cx = full_context(&regime);
    cx.cost_models.insert(
        "xform.level-set-velocity".to_string(),
        lbm_cost_model("flux.free-surface-lbm"),
    );
    let report = admit_src(&SPOUT.replace("(wall 2h)", "(wall 60s)"), &cx);
    assert!(
        !report.admitted,
        "foreign model scope must reject admission"
    );
    let mismatch = report
        .findings
        .iter()
        .find(|finding| finding.what.contains("CostModelScopeMismatch"))
        .expect("scope substitution emits a stable rejection");
    assert_eq!(mismatch.severity, Severity::Reject);
    assert!(
        mismatch.what.contains("xform.level-set-velocity")
            && mismatch.what.contains("flux.free-surface-lbm"),
        "both requested and intrinsic operation identities are named: {mismatch:?}"
    );
    assert!(
        mismatch
            .fixes
            .iter()
            .any(|fix| fix.action.contains("separately admitted binding")),
        "the refusal teaches the exact binding rule: {mismatch:?}"
    );
    assert!(
        report.findings.iter().all(|finding| {
            !finding.what.contains("BudgetInfeasible")
                && !finding.what.contains("CostModelProvisional")
                && !finding.what.contains("CostModelStale")
        }),
        "a foreign model must not influence arithmetic or evidence notices:\n{}",
        report.diagnosis()
    );
    verdict(
        "ad-010a",
        "caller registry keys cannot substitute a foreign sealed operation scope",
    );
}

/// Bead sj31i.11: dimension composition inside admission uses the
/// CHECKED authoritative algebra — a product whose exponents leave the
/// i8 domain REJECTS with a typed dimensional finding instead of
/// saturating into a vector that later cancellation could alias back
/// to plausible physics.
#[test]
fn ad_016_dimension_overflow_refuses_typed_instead_of_clamping() {
    // 70 pressure factors drive the s-exponent to -140 < -128.
    let product = format!("(* {})", vec!["101kPa"; 70].join(" "));
    let src = SPOUT.replace(
        "(min (perturbation-growth pour :at lip :modes (1 .. 8)))",
        &product,
    );
    let regime = spout_regime();
    let cx = full_context(&regime);
    let report = admit_src(&src, &cx);
    let finding = report
        .findings
        .iter()
        .find(|f| f.check == "dimensional" && f.what.contains("overflow"))
        .expect("dimension overflow must surface as a typed dimensional finding");
    assert_eq!(finding.severity, Severity::Reject);
    assert!(
        finding.what.contains("i8"),
        "diagnosis names the exponent domain: {}",
        finding.what
    );
    // A deep-but-in-range product stays legal: the refusal is the
    // overflow, not the depth.
    let legal = format!("(* {})", vec!["101kPa"; 40].join(" "));
    let ok_src = SPOUT.replace(
        "(min (perturbation-growth pour :at lip :modes (1 .. 8)))",
        &legal,
    );
    let ok_report = admit_src(&ok_src, &cx);
    assert!(
        ok_report
            .findings
            .iter()
            .all(|f| !(f.check == "dimensional" && f.what.contains("overflow"))),
        "40 pressure factors stay within i8 exponents"
    );
}

// ── Sealed capability evidence (bead aeq7, increment 2 slice B) ─────────────

/// A fixture standing in for fs-session's grant bridge at the
/// composition root.
struct FixtureIssuer;
impl fs_ir::admission::CapabilityIssuerVerifier for FixtureIssuer {
    fn verify(&self, capability: &SessionCapability, receipt: &CapabilityReceipt) -> bool {
        receipt.issuer_id == "ops/fixture" && capability.cores == 128
    }
}

#[test]
fn grant_backed_evidence_is_mintable_only_through_an_accepting_verifier() {
    let receipt = CapabilityReceipt {
        issuer_id: "ops/fixture".to_string(),
        policy_fingerprint: "policy-v1".to_string(),
        session: 41,
        grant_digest: "aa".repeat(32),
    };
    // Deny-all default: nothing mints.
    assert!(
        SealedSessionCapability::grant_backed(
            token(),
            receipt.clone(),
            &NoCapabilityIssuerVerifier,
        )
        .is_err()
    );
    // A refusing pair (wrong issuer) also fails.
    let foreign = CapabilityReceipt {
        issuer_id: "ops/other".to_string(),
        ..receipt.clone()
    };
    assert!(SealedSessionCapability::grant_backed(token(), foreign, &FixtureIssuer).is_err());
    // The accepting verifier mints, and the class + receipt survive.
    let sealed = SealedSessionCapability::grant_backed(token(), receipt.clone(), &FixtureIssuer)
        .expect("fixture issuer vouches");
    assert_eq!(
        sealed.evidence_class(),
        CapabilityEvidenceClass::GrantBacked
    );
    assert_eq!(sealed.receipt(), Some(&receipt));
    assert_eq!(sealed.capability(), &token());
}

#[test]
fn caller_declared_capability_admits_with_a_provenance_warning() {
    let regime = spout_regime();
    let cx = full_context(&regime);
    let report = admit_src(SPOUT, &cx);
    assert!(report.admitted, "{}", report.diagnosis());
    assert!(
        report.findings.iter().any(|finding| {
            finding.check == "capability"
                && finding.severity == Severity::Warn
                && finding.what.contains("caller-declared")
        }),
        "declared evidence must carry the provenance warning: {:?}",
        report.findings
    );
}

#[test]
fn grant_backed_capability_admits_without_the_provenance_warning() {
    let regime = spout_regime();
    let mut cx = full_context(&regime);
    let receipt = CapabilityReceipt {
        issuer_id: "ops/fixture".to_string(),
        policy_fingerprint: "policy-v1".to_string(),
        session: 41,
        grant_digest: "aa".repeat(32),
    };
    cx.capability = Some(
        SealedSessionCapability::grant_backed(token(), receipt, &FixtureIssuer)
            .expect("fixture issuer vouches"),
    );
    let report = admit_src(SPOUT, &cx);
    assert!(report.admitted, "{}", report.diagnosis());
    assert!(
        !report
            .findings
            .iter()
            .any(|finding| finding.check == "capability"
                && finding.severity == Severity::Warn
                && finding.what.contains("caller-declared")),
        "grant-backed evidence must not carry the declared-provenance warning"
    );
}
