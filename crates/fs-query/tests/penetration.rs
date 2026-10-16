//! Certified convex penetration conformance (bead hk8f5).
//!
//! G0 analytic containment, G3 budget-prefix monotonicity, G4 fail-closed
//! cancellation/refusals, and G5 deterministic replay.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::fixtures::SphereChart;
use fs_geom::{Aabb, Point3, Vec3};
use fs_query::{
    CONVEX_PENETRATION_MAX_ITERATIONS, ConvexBox, ConvexOverlapWitness, ConvexPenetration,
    ConvexSphere, ConvexSupportMap, ImplicitGapOracle, QueryError, convex_overlap_witness,
    convex_penetration_depth,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0xEFA0_0001,
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn sphere(center: [f64; 3], radius: f64) -> ConvexSphere {
    ConvexSphere::new(Point3::new(center[0], center[1], center[2]), radius)
        .expect("valid convex sphere")
}

fn boxx(min: [f64; 3], max: [f64; 3]) -> ConvexBox {
    ConvexBox::new(Aabb::new(
        Point3::new(min[0], min[1], min[2]),
        Point3::new(max[0], max[1], max[2]),
    ))
    .expect("valid convex box")
}

fn sphere_fixture() -> (ConvexSphere, ConvexSphere, ConvexOverlapWitness, f64) {
    let a = sphere([0.0, 0.0, 0.0], 1.5);
    let b = sphere([0.75, 1.0, 0.0], 1.0);
    let center = Point3::new(0.75, 1.0, 0.0);
    let witness = convex_overlap_witness(&a, &b, center).expect("positive common ball");
    // Center distance is exactly 1.25, so depth = 1.5 + 1.0 - 1.25.
    (a, b, witness, 1.25)
}

fn box_fixture() -> (ConvexBox, ConvexBox, ConvexOverlapWitness, f64) {
    let a = boxx([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
    let b = boxx([0.25, -0.5, -0.25], [1.5, 0.75, 0.5]);
    let witness =
        convex_overlap_witness(&a, &b, Point3::new(0.5, 0.0, 0.0)).expect("positive common ball");
    (a, b, witness, 0.75)
}

fn assert_bracket(depth: ConvexPenetration, truth: f64, cap: u32) {
    assert!(depth.lo.is_finite() && depth.hi.is_finite());
    assert!(depth.lo > 0.0 && depth.lo <= depth.hi);
    assert!(
        depth.lo <= truth && truth <= depth.hi,
        "certified bracket [{}, {}] must contain analytic depth {truth}",
        depth.lo,
        depth.hi
    );
    assert!(depth.iterations <= cap.min(CONVEX_PENETRATION_MAX_ITERATIONS));
    assert_eq!(depth.iterations, depth.support_evaluations);
}

struct MalformedSupportMap {
    inner: ConvexSphere,
    bad_point: bool,
    slack: f64,
}

impl ConvexSupportMap for MalformedSupportMap {
    fn support_point(&self, direction: Vec3) -> Point3 {
        if self.bad_point {
            Point3::new(f64::NAN, 0.0, 0.0)
        } else {
            self.inner.support_point(direction)
        }
    }

    fn interior_point(&self) -> Point3 {
        self.inner.interior_point()
    }

    fn support_slack(&self) -> f64 {
        self.slack
    }

    fn contained_ball_radius(&self, center: Point3) -> Option<f64> {
        self.inner.contained_ball_radius(center)
    }

    fn name(&self) -> &'static str {
        "convex/malformed-test"
    }
}

#[test]
fn gp_001_sphere_and_box_depths_are_contained() {
    let (sphere_a, sphere_b, sphere_witness, sphere_truth) = sphere_fixture();
    let sphere_depth =
        with_cx(|cx| convex_penetration_depth(&sphere_a, &sphere_b, &sphere_witness, 256, cx))
            .expect("sphere penetration bracket");
    assert_bracket(sphere_depth, sphere_truth, 256);

    let (box_a, box_b, box_witness, box_truth) = box_fixture();
    let box_depth = with_cx(|cx| convex_penetration_depth(&box_a, &box_b, &box_witness, 256, cx))
        .expect("box penetration bracket");
    assert_bracket(box_depth, box_truth, 256);
}

#[test]
fn gp_002_gap_witness_is_sealed_and_revalidated_against_the_pair() {
    let chart_a = SphereChart {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.5,
    };
    let chart_b = SphereChart {
        center: Point3::new(0.75, 1.0, 0.0),
        radius: 1.0,
    };
    let probe = chart_b.center;
    let gap_witness = with_cx(|cx| {
        ImplicitGapOracle::new(&chart_a, &chart_b)
            .expect("exact-distance pair")
            .gap_at(probe, cx)
            .expect("gap sample")
            .overlap_witness()
            .expect("sealed common-ball witness")
    });
    let convex_a = sphere([0.0, 0.0, 0.0], 1.5);
    let convex_b = sphere([0.75, 1.0, 0.0], 1.0);
    let depth = with_cx(|cx| convex_penetration_depth(&convex_a, &convex_b, &gap_witness, 256, cx))
        .expect("matching support maps revalidate the gap witness");
    assert_bracket(depth, 1.25, 256);

    let unrelated = sphere([4.0, 0.0, 0.0], 0.5);
    let refused =
        with_cx(|cx| convex_penetration_depth(&convex_a, &unrelated, &gap_witness, 64, cx));
    assert!(matches!(
        refused,
        Err(QueryError::ConvexOverlapUnproven { .. })
    ));
}

#[test]
fn gp_003_touching_and_missing_positive_ball_refuse() {
    let a = sphere([-1.0, 0.0, 0.0], 1.0);
    let b = sphere([1.0, 0.0, 0.0], 1.0);
    let touching = convex_overlap_witness(&a, &b, Point3::new(0.0, 0.0, 0.0));
    assert!(matches!(
        touching,
        Err(QueryError::ConvexOverlapUnproven { .. })
    ));

    let (overlap_a, overlap_b, witness, _) = sphere_fixture();
    let zero_budget =
        with_cx(|cx| convex_penetration_depth(&overlap_a, &overlap_b, &witness, 0, cx));
    assert!(matches!(
        zero_budget,
        Err(QueryError::ConvexInvalidShape { .. })
    ));
}

#[test]
fn gp_004_budget_prefixes_tighten_monotonically() {
    let (a, b, witness, truth) = sphere_fixture();
    let mut previous: Option<ConvexPenetration> = None;
    for cap in [1, 2, 4, 8, 16, 32, 64, 128] {
        let depth = with_cx(|cx| convex_penetration_depth(&a, &b, &witness, cap, cx))
            .expect("bounded prefix");
        assert_bracket(depth, truth, cap);
        if let Some(previous) = previous {
            assert!(
                depth.lo >= previous.lo,
                "lower bound regressed at cap {cap}"
            );
            assert!(
                depth.hi <= previous.hi,
                "upper bound regressed at cap {cap}"
            );
        }
        previous = Some(depth);
    }
}

#[test]
fn gp_005_replay_and_cancellation_are_fail_closed() {
    let (a, b, witness, truth) = box_fixture();
    let first = with_cx(|cx| convex_penetration_depth(&a, &b, &witness, 128, cx))
        .expect("first deterministic run");
    let second = with_cx(|cx| convex_penetration_depth(&a, &b, &witness, 128, cx))
        .expect("second deterministic run");
    assert_bracket(first, truth, 128);
    assert_eq!(first.lo.to_bits(), second.lo.to_bits());
    assert_eq!(first.hi.to_bits(), second.hi.to_bits());
    assert_eq!(
        first.normal.map(f64::to_bits),
        second.normal.map(f64::to_bits)
    );
    assert_eq!(first.iterations, second.iterations);
    assert_eq!(first.support_evaluations, second.support_evaluations);
    assert_eq!(first.converged, second.converged);

    let gate = CancelGate::new();
    gate.request();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let cancelled = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0xEFA0_0001,
                kernel_id: 30,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        convex_penetration_depth(&a, &b, &witness, 128, &cx)
    });
    assert!(matches!(cancelled, Err(QueryError::Cancelled)));
}

#[test]
fn gp_006_malformed_support_and_slack_refuse() {
    let a = sphere([0.0, 0.0, 0.0], 1.0);
    let b = sphere([0.25, 0.0, 0.0], 1.0);
    let witness =
        convex_overlap_witness(&a, &b, Point3::new(0.125, 0.0, 0.0)).expect("positive common ball");

    let bad_slack = MalformedSupportMap {
        inner: a,
        bad_point: false,
        slack: f64::NAN,
    };
    let slack_refusal = with_cx(|cx| convex_penetration_depth(&bad_slack, &b, &witness, 16, cx));
    assert!(matches!(
        slack_refusal,
        Err(QueryError::ConvexInvalidSupport { .. })
    ));

    let bad_point = MalformedSupportMap {
        inner: a,
        bad_point: true,
        slack: 0.0,
    };
    let point_refusal = with_cx(|cx| convex_penetration_depth(&bad_point, &b, &witness, 16, cx));
    assert!(matches!(
        point_refusal,
        Err(QueryError::ConvexInvalidSupport { .. })
    ));
}
