//! fs-contact Stage-1 battery (bead tqag, increment 1).
//!
//! - ct-001 G0/G1: analytic screw-motion broad phase — an approach
//!   window yields the pair, a retreat window provably prunes it, both
//!   against hand-computed enclosure geometry.
//! - ct-002 G5: identical inputs replay identical reports.
//! - ct-003 G0: budget exhaustion lists the exact unresolved pairs;
//!   the resolved prefix is never presented as complete.
//! - ct-004 G0: capability refusal names the pair; the Convex×Convex
//!   route contains the analytic distance at a frozen time.

use asupersync::types::Budget;
use fs_contact::{
    BroadPhaseReport, ContactError, NarrowRoute, NarrowVerdict, SpacetimeBody, narrow_phase,
    spacetime_candidates,
};
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_ga::Motor;
use fs_geom::{Aabb, Point3};
use fs_ivl::Interval;
use fs_motion::{CertifiedMotorTube, ScrewParams, screw_tube};
use fs_query::{ConvexSphere, QueryError};

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-contact\",\"case\":\"{case}\",\"verdict\":\"{}\",\
         \"detail\":\"{detail}\"}}",
        if pass { "pass" } else { "fail" }
    );
    assert!(pass, "case {case}: {detail}");
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0xC0A7,
                kernel_id: 21,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

/// A pure-translation tube along `axis` at `speed` (screw with zero
/// angular rate: the enclosure is exact translation plus Taylor slack).
fn translation_tube(axis: [f64; 3], speed: f64, domain: Interval) -> CertifiedMotorTube {
    screw_tube(
        &ScrewParams {
            axis,
            center: [0.0, 0.0, 0.0],
            omega: 0.0,
            axial_velocity: speed,
            base_pose: Motor::identity(),
        },
        domain,
        4,
        8,
    )
    .expect("analytic translation tube")
}

fn unit_box() -> Aabb {
    Aabb::new(Point3::new(-0.5, -0.5, -0.5), Point3::new(0.5, 0.5, 0.5))
}

#[test]
fn ct_001_screw_motion_broad_phase_matches_analytic_geometry() {
    // Body A rides +x at 1 m/s from x=0; body B rides -x at 1 m/s
    // from x=6 (base pose folded into its support box instead of the
    // motor, keeping both tubes purely analytic screws about origin).
    let domain = Interval::new(0.0, 4.0);
    let tube_a = translation_tube([1.0, 0.0, 0.0], 1.0, domain);
    let tube_b = translation_tube([-1.0, 0.0, 0.0], 1.0, domain);
    let support_b = Aabb::new(Point3::new(5.5, -0.5, -0.5), Point3::new(6.5, 0.5, 0.5));
    let bodies = [
        SpacetimeBody::new(unit_box(), &tube_a).expect("body a"),
        SpacetimeBody::new(support_b, &tube_b).expect("body b"),
    ];
    // Early window [0, 1]: A spans ⊆ [-0.5, 1.5], B spans ⊆ [4.5, 6.5]
    // — a certified gap; the pair MUST be pruned.
    let early = with_cx(|cx| spacetime_candidates(&bodies, Interval::new(0.0, 1.0), 16, cx))
        .expect("early window");
    assert!(
        early.pairs.is_empty(),
        "a certified 3-unit gap cannot be a candidate: {:?}",
        early.pairs
    );
    // Full window [0, 4]: at t=3 the boxes provably meet (A reaches
    // [2.5, 3.5], B reaches [2.5, 3.5]); the pair MUST appear.
    let full = with_cx(|cx| spacetime_candidates(&bodies, domain, 16, cx)).expect("full window");
    assert_eq!(full.pairs, vec![(0, 1)], "approach window finds the pair");
    assert!(full.max_defect.is_finite() && full.max_defect >= 0.0);
    verdict(
        "ct-001",
        true,
        &format!(
            "early window pruned ({} checked, {} pruned); full window found {:?}, \
             defect {:.3e}",
            early.checked_pairs, early.pruned_pairs, full.pairs, full.max_defect
        ),
    );
}

#[test]
fn ct_002_reports_replay_identically() {
    let domain = Interval::new(0.0, 2.0);
    let tube = translation_tube([0.0, 0.0, 1.0], 0.25, domain);
    let bodies = [
        SpacetimeBody::new(unit_box(), &tube).expect("a"),
        SpacetimeBody::new(
            Aabb::new(Point3::new(0.25, 0.25, 0.25), Point3::new(1.25, 1.25, 1.25)),
            &tube,
        )
        .expect("b"),
        SpacetimeBody::new(
            Aabb::new(Point3::new(4.0, 4.0, 4.0), Point3::new(5.0, 5.0, 5.0)),
            &tube,
        )
        .expect("c"),
    ];
    let (first, second): (BroadPhaseReport, BroadPhaseReport) = with_cx(|cx| {
        (
            spacetime_candidates(&bodies, domain, 16, cx).expect("first"),
            spacetime_candidates(&bodies, domain, 16, cx).expect("second"),
        )
    });
    assert_eq!(first, second, "broad phase is a pure function");
    assert_eq!(first.pairs, vec![(0, 1)]);
    verdict("ct-002", true, "identical inputs, identical reports");
}

#[test]
fn ct_003_budget_exhaustion_lists_unresolved_pairs() {
    let domain = Interval::new(0.0, 1.0);
    let tube = translation_tube([1.0, 0.0, 0.0], 0.0, domain);
    // Four co-located bodies: 6 overlapping pairs against a budget of 2.
    let bodies: Vec<SpacetimeBody<'_>> = (0..4)
        .map(|_| SpacetimeBody::new(unit_box(), &tube).expect("body"))
        .collect();
    let refused = with_cx(|cx| spacetime_candidates(&bodies, domain, 2, cx));
    match refused {
        Err(ContactError::CandidateBudgetExhausted {
            max_pairs,
            unresolved,
        }) => {
            assert_eq!(max_pairs, 2);
            assert_eq!(
                unresolved.len(),
                4,
                "6 overlapping pairs minus the 2 budgeted must be listed"
            );
            verdict(
                "ct-003",
                true,
                &format!("budget 2 exhausted; unresolved {unresolved:?} listed"),
            );
        }
        other => panic!("expected budget exhaustion, got {other:?}"),
    }
}

#[test]
fn ct_004_capability_routing_and_convex_containment() {
    // Missing capability refuses by name.
    let sphere_a = ConvexSphere::new(Point3::new(-1.0, 0.0, 0.0), 0.5).expect("a");
    let refused = with_cx(|cx| {
        narrow_phase(
            (3, 7),
            &NarrowRoute::Convex(&sphere_a),
            &NarrowRoute::Undeclared,
            256,
            cx,
        )
    });
    match refused {
        Err(ContactError::MissingCapability {
            body_a,
            body_b,
            capability,
        }) => {
            assert_eq!((body_a, body_b), (3, 7));
            assert_eq!(capability, "convex-support-map");
        }
        other => panic!("expected capability refusal, got {other:?}"),
    }
    // Convex route: analytic distance 1.0 between the spheres.
    let sphere_b = ConvexSphere::new(Point3::new(1.0, 0.0, 0.0), 0.5).expect("b");
    let sep = with_cx(|cx| {
        narrow_phase(
            (0, 1),
            &NarrowRoute::Convex(&sphere_a),
            &NarrowRoute::Convex(&sphere_b),
            256,
            cx,
        )
    })
    .expect("convex route");
    let NarrowVerdict::Convex(separation) = sep;
    assert!(
        separation.lo <= 1.0 && 1.0 <= separation.hi,
        "convex verdict [{}, {}] must contain the analytic 1.0",
        separation.lo,
        separation.hi
    );
    assert!(separation.separation_proven);
    // Query refusals pass through typed (zero iteration budget).
    let passthrough = with_cx(|cx| {
        narrow_phase(
            (0, 1),
            &NarrowRoute::Convex(&sphere_a),
            &NarrowRoute::Convex(&sphere_b),
            0,
            cx,
        )
    });
    assert!(matches!(
        passthrough,
        Err(ContactError::Query(QueryError::ConvexInvalidShape { .. }))
    ));
    verdict(
        "ct-004",
        true,
        &format!(
            "capability refusal named (3,7); convex [{:.6}, {:.6}] ∋ 1.0; \
             query refusals pass through",
            separation.lo, separation.hi
        ),
    );
}
