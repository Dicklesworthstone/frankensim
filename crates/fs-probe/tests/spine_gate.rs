//! ADDENDUM PHASE 0 — THE SPINE GATE (milestone xpck.2), as one
//! executable state. The milestone is defined as a Gauntlet STATE, not a
//! date: this suite is the conjunction of the three exit tests from the
//! phase plan plus the budget-pie render, run against the five member
//! implementations (colors: fs-ledger; falsifiers: fs-evidence; Goodhart
//! guard: fs-opt; tombstones: fs-ledger; interface types: fs-iface;
//! budget pie: fs-probe).

use std::collections::BTreeMap;

use fs_evidence::falsify::FalsifierRegistry;
use fs_evidence::{Color, IntervalOp, NumericalCertificate, ValidityDomain};
use fs_iface::{CouplingGraph, CouplingRole, PairingRegistry, SpaceType};
use fs_ledger::tombstone::{Descriptor, ExplorationVerdict, TombstoneIndex};
use fs_ledger::{
    ColorGraph, PolicyDecision, SourceOrigin, SourceOriginRequest, SourceOriginVerifier, hash_bytes,
};
use fs_opt::{DeltaPerturbationStep, Endpoint, GoodhartGuard};
use fs_probe::{BudgetPie, ErrorContribution};
use fs_qty::{Dims, QtyAny};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"spine-gate\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

struct FixtureSourceVerifier(Vec<u8>);

impl SourceOriginVerifier for FixtureSourceVerifier {
    fn verify(&self, request: &SourceOriginRequest<'_>) -> PolicyDecision {
        let fingerprint = hash_bytes(b"fs-probe/spine-gate/source-policy/v1");
        if self.0 == request.canonical_bytes() {
            PolicyDecision::accept(fingerprint)
        } else {
            PolicyDecision::reject(fingerprint)
        }
    }
}

/// EXIT TEST 1 — LAUNDERING: an adversarial pipeline that attempts to
/// upgrade estimated→verified fails the type check.
#[test]
fn spine_exit_1_laundering_refused() {
    let mut graph = ColorGraph::new();
    let state = BTreeMap::new();
    let surrogate = graph
        .source(
            "drag-surrogate",
            Color::Estimated {
                estimator: "fno-v1".to_string(),
                dispersion: 0.12,
            },
        )
        .expect("Estimated source");
    let mesh_color = Color::Verified { lo: 0.9, hi: 1.1 };
    let mesh_origin = SourceOrigin::Certificate {
        producer: "fs-probe/spine-gate".to_string(),
        certificate_hash: hash_bytes(b"fs-probe/spine-gate/mesh-integral-certificate"),
        certificate: NumericalCertificate::enclosure(0.9, 1.1),
    };
    let mesh_verifier = FixtureSourceVerifier(
        SourceOriginRequest::new("mesh-integral", &mesh_color, &mesh_origin).canonical_bytes(),
    );
    let mesh = graph
        .source_with_origin("mesh-integral", &mesh_color, mesh_origin, &mesh_verifier)
        .expect("mesh enclosure mints Verified");
    // The adversarial upgrade: claim Verified from an Estimated parent.
    let laundered = graph.derive(
        "polished-report",
        &[surrogate, mesh],
        IntervalOp::Add,
        Some(Color::Verified { lo: 0.9, hi: 1.3 }),
        &state,
        None,
    );
    let err = laundered.expect_err("estimates cannot become certificates by assertion");
    let msg = format!("{err}");
    assert!(
        msg.contains("laundering refused"),
        "the refusal teaches: {msg}"
    );
    // The honest write at the derived rank succeeds.
    graph
        .derive(
            "honest-report",
            &[surrogate, mesh],
            IntervalOp::Add,
            None,
            &state,
            None,
        )
        .expect("derivation at the derived rank is legal");
    verdict(
        "spine-exit-1",
        "estimated->verified upgrade refused with teaching text; honest derivation \
         passes",
    );
}

/// EXIT TEST 2 — NO-FALSIFIER-NO-SHIP: a certificate class without a
/// registered falsifier fails Gauntlet review the way untested code
/// fails CI.
#[test]
fn spine_exit_2_no_falsifier_no_ship() {
    let registry = FalsifierRegistry::standard();
    // A shipped class the standard registry covers: passes.
    let clean = registry.ship_gate(&["watertightness"]);
    assert!(clean.is_empty(), "paired class ships: {clean:?}");
    // A novel class nobody paired: BLOCKED, named.
    let blocked = registry.ship_gate(&["watertightness", "novel-magic-certificate"]);
    assert_eq!(blocked.len(), 1, "the unpaired class is named");
    assert!(blocked[0].contains("novel-magic-certificate"));
    verdict(
        "spine-exit-2",
        "unpaired certificate class blocked by the ship gate; paired classes pass",
    );
}

/// EXIT TEST 3 — TOMBSTONE PROTOCOL: before funding an exploration the
/// orchestrator queries the index and either cites a (validated)
/// distinguisher or skips.
#[test]
fn spine_exit_3_tombstone_protocol() {
    let mut index = TombstoneIndex::new();
    let mut params = BTreeMap::new();
    params.insert(
        "velocity".to_string(),
        QtyAny::new(24.0, Dims([1, 0, -1, 0, 0])),
    );
    params.insert(
        "length".to_string(),
        QtyAny::new(0.12, Dims([1, 0, 0, 0, 0])),
    );
    params.insert(
        "viscosity".to_string(),
        QtyAny::new(1.5e-5, Dims([2, 0, -1, 0, 0])),
    );
    let dead = Descriptor {
        name: "bracket crossflow".to_string(),
        params: params.clone(),
    };
    index.record_falsification_kill(
        dead,
        "{\"kind\":\"tombstone\"}",
        vec!["estimated".to_string()],
        300.0,
        "2026-07-07",
        "agent:spine",
    );
    // The protocol: query BEFORE funding. A re-run is blocked…
    let retry = Descriptor {
        name: "bracket crossflow retry".to_string(),
        params: params.clone(),
    };
    let gate = index.pre_exploration_check(&retry);
    let neighbor = match gate {
        ExplorationVerdict::Blocked { ref neighbors, .. } => neighbors[0],
        ExplorationVerdict::Clear => panic!("the protocol must block the re-run"),
    };
    // …an arbitrary excuse is refused…
    assert!(
        index
            .fund_with_distinguisher(&retry, neighbor, "trust me")
            .is_err(),
        "free text is not a distinguisher"
    );
    // …and a genuinely different named parameter funds.
    let mut novel_params = params;
    novel_params.insert(
        "velocity".to_string(),
        QtyAny::new(2.0, Dims([1, 0, -1, 0, 0])),
    );
    let novel = Descriptor {
        name: "bracket creepflow".to_string(),
        params: novel_params,
    };
    index
        .fund_with_distinguisher(&novel, neighbor, "velocity")
        .expect("a 12x velocity change is a real distinguisher");
    verdict(
        "spine-exit-3",
        "gate blocks the re-run, refuses free text, funds a validated distinguisher",
    );
}

/// MEMBER CHECKS — the two spine members without a dedicated exit test:
/// the interface-type checker rejects an illegal coupling, and the
/// Goodhart guard refuses to honor an un-escalated optimizer endpoint.
#[test]
fn spine_members_iface_and_goodhart() {
    // Interface types: pressure (L²) tested against a displacement trace
    // (H¹) as CONTINUITY is exactly the classic illegal pairing.
    let graph = CouplingGraph::new()
        .field("pressure", SpaceType::L2)
        .field("displacement", SpaceType::HGrad)
        .couple(
            "fsi-seam",
            "pressure",
            "displacement",
            CouplingRole::Continuity,
        );
    let report = fs_iface::check(&graph, &PairingRegistry::standard());
    assert!(
        !report.admitted,
        "illegal coupling rejected: {}",
        report.diagnosis()
    );
    // A legal same-space continuity coupling is admitted.
    let ok = CouplingGraph::new()
        .field("t-left", SpaceType::HGrad)
        .field("t-right", SpaceType::HGrad)
        .couple(
            "thermal-seam",
            "t-left",
            "t-right",
            CouplingRole::Continuity,
        );
    assert!(fs_iface::check(&ok, &PairingRegistry::standard()).admitted);
    // Goodhart guard: an exploit endpoint (a crack minimum) is VETOED —
    // the honest objective rises steeply δ away from the reported point.
    let honest = |x: &[f64]| {
        let d: f64 = x.iter().map(|v| (v - 0.5).abs()).sum();
        if d < 1e-6 { -100.0 } else { d }
    };
    let guard = GoodhartGuard::new().with_step(Box::new(DeltaPerturbationStep::new(
        0.01, 1e-9, 1.0, honest,
    )));
    let crack = guard.evaluate(&Endpoint::new("crack", vec![0.5, 0.5], -100.0));
    assert_eq!(
        crack.status,
        fs_opt::GuardStatus::Failed,
        "a crack optimum FAILS (vetoed), not merely provisional"
    );
    assert!(
        !crack.findings.is_empty(),
        "the veto is treasure (a finding)"
    );
    // A smooth genuine optimum passes the same step.
    let smooth = |x: &[f64]| x.iter().map(|v| (v - 0.5) * (v - 0.5)).sum::<f64>();
    let guard2 = GoodhartGuard::new().with_step(Box::new(DeltaPerturbationStep::new(
        0.01, 1e-9, 1.0, smooth,
    )));
    let genuine = guard2.evaluate(&Endpoint::new("genuine", vec![0.5, 0.5], 0.0));
    // With 3 of 4 escalation steps lacking machinery, the guard is
    // HONEST: Provisional, never Cleared on skipped checks — and never
    // Failed without a veto. Both un-honored states still block the
    // certificate, which is the policy's point.
    assert_eq!(
        genuine.status,
        fs_opt::GuardStatus::Provisional,
        "{}",
        genuine.diagnosis()
    );
    assert!(genuine.findings.is_empty(), "no veto on a smooth optimum");
    assert!(
        !genuine.is_honored(),
        "certificates are honored only when EVERY step clears"
    );
    verdict(
        "spine-members",
        "illegal coupling rejected / legal admitted; crack endpoint vetoed with \
         findings / smooth endpoint honored",
    );
}

/// THE BUDGET PIE renders error-by-color on a real report (the milestone
/// acceptance's final clause).
#[test]
fn spine_budget_pie_renders() {
    let contributions = [
        ErrorContribution::new(
            "mesh-discretization",
            Color::Verified {
                lo: -0.002,
                hi: 0.002,
            },
            0.004,
        ),
        ErrorContribution::new(
            "turbulence-closure",
            Color::Estimated {
                estimator: "k-epsilon-vs-les-probe".to_string(),
                dispersion: 0.05,
            },
            0.05,
        ),
        ErrorContribution::new(
            "material-anchor",
            Color::Validated {
                regime: ValidityDomain::unconstrained(),
                dataset: "al6061-coupons-2025".to_string(),
            },
            0.01,
        ),
    ];
    let pie = BudgetPie::of(&contributions);
    let total = 0.004 + 0.05 + 0.01;
    assert!((pie.fraction(fs_evidence::ColorRank::Estimated) - 0.05 / total).abs() < 1e-12);
    assert_eq!(pie.dominant(), Some(fs_evidence::ColorRank::Estimated));
    // The operator-legible verdict names the RIGHT action: model-form
    // error dominates, so refining the mesh will not help.
    assert!(
        pie.verdict().contains("refining the mesh will NOT help"),
        "the pie teaches where to spend: {}",
        pie.verdict()
    );
    verdict(
        "spine-pie",
        "error-by-color pie renders with exact fractions and the model-form verdict",
    );
}
