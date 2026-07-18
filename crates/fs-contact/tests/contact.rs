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

// ── Increment 2: certified CCD battery (ct-005..ct-008) ─────────────────────

use fs_contact::{CcdVerdict, certified_ccd};

fn thin_plate() -> Aabb {
    // A wall in the yz-plane: 2 cm thick, 4 m tall/wide.
    Aabb::new(Point3::new(-0.01, -2.0, -2.0), Point3::new(0.01, 2.0, 2.0))
}

fn bullet() -> Aabb {
    // A 2 cm cube.
    Aabb::new(
        Point3::new(-0.01, -0.01, -0.01),
        Point3::new(0.01, 0.01, 0.01),
    )
}

/// ct-005 Sev-0 G3: a fast bullet fully crosses a thin static plate
/// INSIDE the window — endpoint sampling provably misses it (both
/// endpoint enclosures are disjoint from the plate), yet certified CCD
/// must report a possible-contact window containing the true crossing.
#[test]
fn ct_005_thin_fast_crossing_is_never_reported_clear() {
    with_cx(|cx| {
        let domain = Interval::new(0.0, 1.0);
        // Plate static at the origin; bullet starts at x=-50 riding
        // +x at 100 m/s: true wall crossing at t* = 0.5 (x sweeps
        // [-0.02, 0.02] against the plate's [-0.01, 0.01] near
        // t ∈ [0.4997, 0.5003]).
        let plate_tube = translation_tube([1.0, 0.0, 0.0], 0.0, domain);
        let mut bullet_support = bullet();
        bullet_support.min.x -= 50.0;
        bullet_support.max.x -= 50.0;
        let bullet_tube = translation_tube([1.0, 0.0, 0.0], 100.0, domain);
        let plate = SpacetimeBody::new(thin_plate(), &plate_tube).expect("plate body");
        let shot = SpacetimeBody::new(bullet_support, &bullet_tube).expect("bullet body");

        // The sampling counterexample: at both window endpoints the
        // enclosures are far apart.
        for t in [domain.lo(), domain.hi()] {
            let inst = Interval::point(t);
            let pb = plate_tube
                .box_action_over(&thin_plate(), inst, cx)
                .expect("plate endpoint enclosure");
            let bb = bullet_tube
                .box_action_over(
                    &{
                        let mut s = bullet();
                        s.min.x -= 50.0;
                        s.max.x -= 50.0;
                        s
                    },
                    inst,
                    cx,
                )
                .expect("bullet endpoint enclosure");
            let disjoint_x = bb.bounds.max.x < pb.bounds.min.x || pb.bounds.max.x < bb.bounds.min.x;
            verdict(
                "ct-005-endpoints",
                disjoint_x,
                "endpoint sampling must see disjoint bodies (that is the trap)",
            );
        }

        let report =
            certified_ccd(&plate, &shot, domain, 1e-4, 1 << 16, cx).expect("ccd completes");
        let CcdVerdict::PossibleContact { windows } = &report.verdict else {
            verdict(
                "ct-005",
                false,
                "a real crossing must never be reported ClearWindow",
            );
            unreachable!()
        };
        let crossing_covered = windows.iter().any(|w| w.lo() <= 0.5 && 0.5 <= w.hi());
        verdict(
            "ct-005",
            crossing_covered,
            "some possible-contact window must contain the true crossing t*=0.5",
        );
        // The report localizes the event: unresolved time is a sliver of
        // the window, not a give-up.
        let total: f64 = windows.iter().map(|w| w.width()).sum();
        verdict(
            "ct-005-localized",
            total < 0.01,
            "possible windows must localize the crossing to well under 1% of the window",
        );
    });
}

/// ct-006 G0/G5: a margin-separated parallel pass is PROVEN clear with a
/// positive certified gap, and the report replays bit-identically.
#[test]
fn ct_006_grazing_pass_is_proven_clear_and_replays() {
    with_cx(|cx| {
        let domain = Interval::new(0.0, 4.0);
        let tube_a = translation_tube([1.0, 0.0, 0.0], 1.0, domain);
        let tube_b = translation_tube([1.0, 0.0, 0.0], -1.0, domain);
        let a_support = unit_box();
        let mut b_support = unit_box();
        // Parallel tracks offset by 3 in y: closest approach keeps a
        // >= 2 m gap forever.
        b_support.min.y += 3.0;
        b_support.max.y += 3.0;
        let a = SpacetimeBody::new(a_support, &tube_a).expect("body a");
        let b = SpacetimeBody::new(b_support, &tube_b).expect("body b");
        let first = certified_ccd(&a, &b, domain, 1e-3, 1 << 12, cx).expect("ccd completes");
        let CcdVerdict::ClearWindow { min_gap } = first.verdict else {
            verdict(
                "ct-006",
                false,
                "a margin-separated pass must be proven clear",
            );
            unreachable!()
        };
        verdict(
            "ct-006",
            min_gap > 1.0,
            "the certified gap must reflect the >=2m analytic margin (allowing enclosure slack)",
        );
        let replay = certified_ccd(&a, &b, domain, 1e-3, 1 << 12, cx).expect("ccd replays");
        verdict(
            "ct-006-replay",
            replay == first,
            "identical inputs must replay an identical report",
        );
    });
}

/// ct-007 G3 (the global-root-guard refusal, executable): two bodies
/// overlapping for the WHOLE window have no separation sign change for
/// a root guard to find — certified CCD must stay honest and report
/// possible contact covering essentially the entire window, never
/// ClearWindow.
#[test]
fn ct_007_persistent_contact_is_never_cleared() {
    with_cx(|cx| {
        let domain = Interval::new(0.0, 1.0);
        let tube = translation_tube([1.0, 0.0, 0.0], 0.0, domain);
        let a = SpacetimeBody::new(unit_box(), &tube).expect("body a");
        let mut shifted = unit_box();
        shifted.min.x += 0.25;
        shifted.max.x += 0.25;
        let b = SpacetimeBody::new(shifted, &tube).expect("body b");
        let report = certified_ccd(&a, &b, domain, 1e-2, 1 << 12, cx).expect("ccd completes");
        let CcdVerdict::PossibleContact { windows } = &report.verdict else {
            verdict(
                "ct-007",
                false,
                "persistently overlapping bodies must never be proven clear",
            );
            unreachable!()
        };
        let covered: f64 = windows.iter().map(|w| w.width()).sum();
        verdict(
            "ct-007",
            (covered - domain.width()).abs() < 1e-9 && windows.len() == 1,
            "persistent contact must surface as one window covering the whole domain",
        );
    });
}

/// ct-008 G0: budget exhaustion returns the exact partial state and is
/// a refusal — the resolved prefix is never presented as a verdict.
#[test]
fn ct_008_ccd_budget_exhaustion_lists_partial_state() {
    with_cx(|cx| {
        let domain = Interval::new(0.0, 1.0);
        let plate_tube = translation_tube([1.0, 0.0, 0.0], 0.0, domain);
        let mut bullet_support = bullet();
        bullet_support.min.x -= 50.0;
        bullet_support.max.x -= 50.0;
        let bullet_tube = translation_tube([1.0, 0.0, 0.0], 100.0, domain);
        let plate = SpacetimeBody::new(thin_plate(), &plate_tube).expect("plate body");
        let shot = SpacetimeBody::new(bullet_support, &bullet_tube).expect("bullet body");
        match certified_ccd(&plate, &shot, domain, 1e-4, 5, cx) {
            Err(ContactError::CcdBudgetExhausted {
                max_windows,
                examined,
                pending,
                ..
            }) => {
                verdict(
                    "ct-008",
                    max_windows == 5 && examined == 5 && !pending.is_empty(),
                    "exhaustion must report the exact budget, count, and pending windows",
                );
                let ascending = pending.windows(2).all(|p| p[0].hi() <= p[1].lo() + 1e-12);
                verdict(
                    "ct-008-order",
                    ascending,
                    "pending windows must be reported in ascending time order",
                );
            }
            other => verdict(
                "ct-008",
                false,
                &format!("a starved budget must refuse, got {other:?}"),
            ),
        }
    });
}

// ── Increment 3: swept-vertex-hull refinement (ct-009, ct-010) ──────────────

use fs_contact::{RefinedWindow, refine_possible_windows};

/// Vertices of a cube with half-extent `h`, rotated 45° about z, offset
/// by `(dx, dy)` — in body frame.
fn rotated_cube(h: f64, dx: f64, dy: f64) -> Vec<Point3> {
    let r = h * std::f64::consts::SQRT_2;
    vec![
        Point3::new(dx + r, dy, -h),
        Point3::new(dx, dy + r, -h),
        Point3::new(dx - r, dy, -h),
        Point3::new(dx, dy - r, -h),
        Point3::new(dx + r, dy, h),
        Point3::new(dx, dy + r, h),
        Point3::new(dx - r, dy, h),
        Point3::new(dx, dy - r, h),
    ]
}

/// ct-009 G1: two static 45°-rotated cubes passing corner-to-corner —
/// their axis-aligned boxes overlap at EVERY instant (the box verdict
/// can never clear this window at any tolerance), yet the actual bodies
/// keep an analytic gap. The hull refinement must prune with a certified
/// gap close to the analytic value.
#[test]
fn ct_009_rotated_near_miss_boxes_cannot_clear_but_hulls_prune() {
    with_cx(|cx| {
        let domain = Interval::new(0.0, 1.0);
        let tube = translation_tube([1.0, 0.0, 0.0], 0.0, domain);
        // Two diamonds (45°-rotated cubes, |x|+|y| <= sqrt(2)) with B's
        // center at (c, c) on the diagonal, sqrt(2) < c < 2*sqrt(2): the
        // facing EDGES are parallel with certified gap sqrt(2)*(c-sqrt(2)),
        // while both AABBs (side 2*sqrt(2)) still overlap in x AND y — the
        // structural trap the box route can never clear at any tolerance.
        let c = 2.0;
        let a_verts = rotated_cube(1.0, 0.0, 0.0);
        let b_verts = rotated_cube(1.0, c, c);

        // Establish the trap: the BOX route retains the whole window.
        let a_aabb = Aabb::new(
            Point3::new(-std::f64::consts::SQRT_2, -std::f64::consts::SQRT_2, -1.0),
            Point3::new(std::f64::consts::SQRT_2, std::f64::consts::SQRT_2, 1.0),
        );
        let b_aabb = Aabb::new(
            Point3::new(
                c - std::f64::consts::SQRT_2,
                c - std::f64::consts::SQRT_2,
                -1.0,
            ),
            Point3::new(
                c + std::f64::consts::SQRT_2,
                c + std::f64::consts::SQRT_2,
                1.0,
            ),
        );
        let a_body = SpacetimeBody::new(a_aabb, &tube).expect("body a");
        let b_body = SpacetimeBody::new(b_aabb, &tube).expect("body b");
        let boxed =
            certified_ccd(&a_body, &b_body, domain, 1e-2, 1 << 12, cx).expect("box ccd completes");
        let CcdVerdict::PossibleContact { windows } = &boxed.verdict else {
            verdict(
                "ct-009-trap",
                false,
                "overlapping AABBs must retain the window under the box route",
            );
            unreachable!()
        };

        // The refinement prunes every retained window with a certified gap.
        let refined = refine_possible_windows(&a_verts, &tube, &b_verts, &tube, windows, 256, cx)
            .expect("refinement completes");
        // Analytic gap sqrt(2)*(2-sqrt(2)) ≈ 0.828; allow enclosure slack.
        let all_pruned = refined
            .windows
            .iter()
            .all(|w| matches!(w, RefinedWindow::Pruned { gap, .. } if *gap > 0.5));
        verdict(
            "ct-009",
            all_pruned,
            "the hull route must prune the rotated near-miss with a certified gap",
        );
    });
}

/// ct-010 G0 (refinement soundness): the ct-005 bullet's true crossing
/// window must SURVIVE refinement — a certified-separation prune can
/// never drop a window containing real contact.
#[test]
fn ct_010_refinement_never_drops_a_true_crossing() {
    with_cx(|cx| {
        let domain = Interval::new(0.0, 1.0);
        let plate_tube = translation_tube([1.0, 0.0, 0.0], 0.0, domain);
        let bullet_tube = translation_tube([1.0, 0.0, 0.0], 100.0, domain);
        let plate_verts = vec![
            Point3::new(-0.01, -2.0, -2.0),
            Point3::new(0.01, -2.0, -2.0),
            Point3::new(-0.01, 2.0, -2.0),
            Point3::new(0.01, 2.0, -2.0),
            Point3::new(-0.01, -2.0, 2.0),
            Point3::new(0.01, -2.0, 2.0),
            Point3::new(-0.01, 2.0, 2.0),
            Point3::new(0.01, 2.0, 2.0),
        ];
        let bullet_verts: Vec<Point3> = (0..8)
            .map(|i| {
                Point3::new(
                    -50.0 + if i & 1 == 0 { -0.01 } else { 0.01 },
                    if i & 2 == 0 { -0.01 } else { 0.01 },
                    if i & 4 == 0 { -0.01 } else { 0.01 },
                )
            })
            .collect();
        let mut bullet_support = bullet();
        bullet_support.min.x -= 50.0;
        bullet_support.max.x -= 50.0;
        let plate = SpacetimeBody::new(thin_plate(), &plate_tube).expect("plate body");
        let shot = SpacetimeBody::new(bullet_support, &bullet_tube).expect("bullet body");
        let report =
            certified_ccd(&plate, &shot, domain, 1e-4, 1 << 16, cx).expect("ccd completes");
        let CcdVerdict::PossibleContact { windows } = &report.verdict else {
            unreachable!("ct-005 already proves this arm");
        };
        let refined = refine_possible_windows(
            &plate_verts,
            &plate_tube,
            &bullet_verts,
            &bullet_tube,
            windows,
            256,
            cx,
        )
        .expect("refinement completes");
        let crossing_retained = refined.windows.iter().any(|w| {
            matches!(w, RefinedWindow::Retained { window } if window.lo() <= 0.5 && 0.5 <= window.hi())
        });
        verdict(
            "ct-010",
            crossing_retained,
            "the true crossing window must survive refinement as Retained",
        );
    });
}

// ── SDF-obstacle route (ct-011, ct-012) ─────────────────────────────────────

use fs_contact::refine_windows_against_sdf;
use fs_geom::Chart as _;
use fs_geom::fixtures::{LyingSphereChart, SphereChart};

/// ct-011 G1: a bullet passing the sphere's AABB-corner region — the box
/// route retains (AABBs overlap) but the ball-Lipschitz bound prunes
/// fine windows against the actual curved surface with a certified gap;
/// and a weaker-claim chart refuses at entry.
#[test]
fn ct_011_sdf_route_prunes_near_miss_and_gates_the_claim() {
    with_cx(|cx| {
        let domain = Interval::new(0.0, 1.0);
        // Sphere radius 1 at origin. Bullet flies parallel to x at
        // y = z = 0.9: it cuts through the sphere's AABB corner region
        // (|y|,|z| < 1) while staying sqrt(0.81+0.81) ≈ 1.273 > 1 from
        // the center — a real miss the AABB test cannot see.
        let sphere = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let sphere_aabb = sphere.support();
        let static_tube = translation_tube([1.0, 0.0, 0.0], 0.0, domain);
        let bullet_tube = translation_tube([1.0, 0.0, 0.0], 8.0, domain);
        let bullet_verts: Vec<Point3> = (0..8)
            .map(|i| {
                Point3::new(
                    -4.0 + if i & 1 == 0 { -0.01 } else { 0.01 },
                    0.9 + if i & 2 == 0 { -0.01 } else { 0.01 },
                    0.9 + if i & 4 == 0 { -0.01 } else { 0.01 },
                )
            })
            .collect();
        let mut bullet_support = bullet();
        for p in [&mut bullet_support.min, &mut bullet_support.max] {
            p.x -= 4.0;
            p.y += 0.9;
            p.z += 0.9;
        }
        let obstacle_body = SpacetimeBody::new(sphere_aabb, &static_tube).expect("obstacle");
        let shot = SpacetimeBody::new(bullet_support, &bullet_tube).expect("bullet");
        let report = certified_ccd(&obstacle_body, &shot, domain, 1e-3, 1 << 14, cx).expect("ccd");
        let CcdVerdict::PossibleContact { windows } = &report.verdict else {
            verdict(
                "ct-011-trap",
                false,
                "the corner pass must be retained by the AABB route",
            );
            unreachable!()
        };
        let refined = refine_windows_against_sdf(
            &bullet_verts,
            &bullet_tube,
            &sphere,
            windows,
            1e-3,
            1 << 14,
            cx,
        )
        .expect("sdf refinement");
        let all_pruned = refined
            .windows
            .iter()
            .all(|w| matches!(w, RefinedWindow::Pruned { gap, .. } if *gap > 0.0));
        verdict(
            "ct-011",
            all_pruned,
            "sub-tolerance windows must prune against the curved surface",
        );

        // Claim gating: a chart with a weaker trace claim refuses at entry.
        let lying = LyingSphereChart {
            sphere: SphereChart {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 1.0,
            },
            bias: 0.5,
        };
        verdict(
            "ct-011-gate",
            matches!(
                refine_windows_against_sdf(
                    &bullet_verts,
                    &bullet_tube,
                    &lying,
                    windows,
                    1e-3,
                    1 << 14,
                    cx
                ),
                Err(ContactError::MissingCapability {
                    capability: "exact-distance-chart",
                    ..
                })
            ),
            "a weaker-claim chart must refuse by capability name at entry",
        );
    });
}

/// ct-012 G0 (soundness arm): a bullet aimed THROUGH the sphere keeps
/// its crossing windows Retained — the ball-Lipschitz prune can never
/// clear real contact.
#[test]
fn ct_012_sdf_route_never_prunes_a_real_hit() {
    with_cx(|cx| {
        let domain = Interval::new(0.0, 1.0);
        let sphere = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let static_tube = translation_tube([1.0, 0.0, 0.0], 0.0, domain);
        let bullet_tube = translation_tube([1.0, 0.0, 0.0], 8.0, domain);
        let bullet_verts: Vec<Point3> = (0..8)
            .map(|i| {
                Point3::new(
                    -4.0 + if i & 1 == 0 { -0.01 } else { 0.01 },
                    if i & 2 == 0 { -0.01 } else { 0.01 },
                    if i & 4 == 0 { -0.01 } else { 0.01 },
                )
            })
            .collect();
        let mut bullet_support = bullet();
        bullet_support.min.x -= 4.0;
        bullet_support.max.x -= 4.0;
        let obstacle_body = SpacetimeBody::new(sphere.support(), &static_tube).expect("obstacle");
        let shot = SpacetimeBody::new(bullet_support, &bullet_tube).expect("bullet");
        let report = certified_ccd(&obstacle_body, &shot, domain, 1e-3, 1 << 14, cx).expect("ccd");
        let CcdVerdict::PossibleContact { windows } = &report.verdict else {
            unreachable!("a through-shot must be retained by the AABB route");
        };
        let refined = refine_windows_against_sdf(
            &bullet_verts,
            &bullet_tube,
            &sphere,
            windows,
            1e-3,
            1 << 14,
            cx,
        )
        .expect("sdf refinement");
        // Entry time of the bullet center into the sphere: x(t) = -4 + 8t
        // reaches -1 at t = 0.375.
        let entry_retained = refined.windows.iter().any(|w| {
            matches!(w, RefinedWindow::Retained { window }
                if window.lo() <= 0.375 && 0.375 <= window.hi())
        });
        verdict(
            "ct-012",
            entry_retained,
            "the true sphere-entry window must survive SDF refinement as Retained",
        );
    });
}
