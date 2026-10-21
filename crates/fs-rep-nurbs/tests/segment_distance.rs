//! Continuous segment and swept-ball clearance, including exact trim coupling.
#![cfg(feature = "nurbs-sdf")]

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_ivl::Interval;
use fs_rep_nurbs::{
    KnotVector, NurbsCurve, NurbsSurface, Rat, SegmentDistance, SegmentDistanceBudget,
    SegmentDistanceRun, SegmentDistanceStop, SurfaceDistanceBudget, SweptBallClearance, TrimLoop,
    TrimmedPatch,
};

fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 41,
                kernel_id: 59,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn budget(splits: u32) -> SegmentDistanceBudget {
    SegmentDistanceBudget {
        max_splits: splits,
        point_budget: SurfaceDistanceBudget {
            max_splits: 0,
            max_work_units: 16_000_000,
            max_retained_bytes: 96 * 1024 * 1024,
        },
        max_work_units: 256_000_000,
        max_retained_bytes: 100 * 1024 * 1024,
    }
}

fn plane() -> NurbsSurface<f64> {
    let k = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap();
    NurbsSurface::new(
        k.clone(),
        k,
        &[
            vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
        ],
        &[vec![1.0; 2], vec![1.0; 2]],
    )
    .unwrap()
}

fn rectangle(lo: [i128; 2], hi: [i128; 2], den: i128) -> TrimLoop {
    let lo = lo.map(|n| Rat::new(n, den));
    let hi = hi.map(|n| Rat::new(n, den));
    let k = KnotVector::new([0, 0, 1, 2, 3, 4, 4].into_iter().map(Rat::int).collect(), 1).unwrap();
    TrimLoop::new(
        NurbsCurve::new(
            k,
            &[lo, [hi[0], lo[1]], hi, [lo[0], hi[1]], lo],
            &[Rat::int(1); 5],
        )
        .unwrap(),
    )
    .unwrap()
}

fn nonempty<'s, 't>(run: SegmentDistanceRun<'s, 't>) -> SegmentDistance<'s, 't> {
    match run {
        SegmentDistanceRun::Complete { distance }
        | SegmentDistanceRun::Unresolved { distance, .. } => distance,
        other => panic!("expected distance, got {other:?}"),
    }
}

#[test]
fn g1_parallel_toolpath_encloses_continuous_clearance_and_decides_swept_radius() {
    let surface = plane();
    let source = surface.admit().unwrap();
    let endpoints = [[0.25, 0.5, 1.0], [0.75, 0.5, 1.0]];
    with_cx(false, |cx| {
        let zero = nonempty(
            source
                .segment_distance_enclosure_with_cx(endpoints, None, 0.02, budget(0), cx)
                .unwrap(),
        );
        let SegmentDistanceRun::Complete { distance } = source
            .segment_distance_enclosure_with_cx(endpoints, None, 0.02, budget(128), cx)
            .unwrap()
        else {
            panic!("planar toolpath should close its width");
        };
        let bounds = distance.bounds().unwrap();
        assert!(bounds.contains(1.0) && bounds.width() <= 0.02, "{bounds:?}");
        assert!(bounds.lo() >= zero.bounds().unwrap().lo());
        assert!(bounds.hi() <= zero.bounds().unwrap().hi());
        assert!(std::ptr::eq(distance.source(), &surface));
        assert!(distance.trim().is_none());
        assert_eq!(distance.endpoints(), endpoints);
        assert_eq!(
            distance.swept_ball_clearance(0.9).unwrap(),
            SweptBallClearance::Separated
        );
        assert_eq!(
            distance.swept_ball_clearance(1.1).unwrap(),
            SweptBallClearance::WithinRadius
        );
        assert_eq!(
            distance.swept_ball_clearance(1.0).unwrap(),
            SweptBallClearance::Undetermined
        );
        let witness = distance.witness().unwrap();
        assert!(witness.segment_position[0].contains(0.25 + 0.5 * witness.segment_parameter));
        assert!(witness.segment_position[2].contains(1.0));
        assert!(witness.surface_position[2].contains(0.0));
        let work = distance.work();
        assert_eq!(
            1 + work.splits as usize,
            work.remaining_cells + work.excluded_cells
        );
        assert!(work.charged_work_units <= budget(128).max_work_units);
        assert!(work.retained_payload_bytes <= budget(128).max_retained_bytes);
    });
}

#[test]
fn g1_crossing_endpoint_minimum_and_degenerate_segment_keep_their_true_distance() {
    let surface = plane();
    let source = surface.admit().unwrap();
    with_cx(false, |cx| {
        for (endpoints, truth) in [
            ([[0.5, 0.5, -1.0], [0.5, 0.5, 1.0]], 0.0),
            ([[0.5, 0.5, 1.0], [0.5, 0.5, 3.0]], 1.0),
            ([[0.5, 0.5, 2.0]; 2], 2.0),
        ] {
            let SegmentDistanceRun::Complete { distance } = source
                .segment_distance_enclosure_with_cx(endpoints, None, 0.01, budget(128), cx)
                .unwrap()
            else {
                panic!("analytic segment should close its width");
            };
            assert!(distance.bounds().unwrap().contains(truth));
            if truth == 0.0 {
                assert_eq!(
                    distance.swept_ball_clearance(0.0).unwrap(),
                    SweptBallClearance::Undetermined
                );
                assert_eq!(
                    distance.swept_ball_clearance(0.01).unwrap(),
                    SweptBallClearance::WithinRadius
                );
            }
        }
    });
}

#[test]
fn g3_rounded_center_off_the_true_segment_never_becomes_an_exact_pair() {
    let k = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap();
    let surface = NurbsSurface::new(
        k.clone(),
        k,
        &[vec![[1.0, 1.0, 0.0]; 2], vec![[1.0, 1.0, 0.0]; 2]],
        &[vec![1.0; 2], vec![1.0; 2]],
    )
    .unwrap();
    let eps = f64::EPSILON;
    // At t=1/2, x is exactly 1 + eps/2, which f64 cannot represent.
    // Rounding that point onto the constant surface must not prove zero distance.
    let endpoints = [[1.0, 0.0, 0.0], [1.0 + eps, 2.0, 0.0]];
    let truth = Interval::point(eps)
        / (Interval::point(4.0) + Interval::point(eps) * Interval::point(eps)).sqrt();
    with_cx(false, |cx| {
        let distance = nonempty(
            surface
                .admit()
                .unwrap()
                .segment_distance_enclosure_with_cx(endpoints, None, 1e-30, budget(0), cx)
                .unwrap(),
        );
        assert!(distance.bounds().unwrap().encloses(truth));
        let witness = distance.witness().unwrap();
        assert_eq!(witness.segment_parameter, 0.5);
        assert!(witness.segment_position[0].lo() <= 1.0);
        assert!(witness.segment_position[0].hi() >= 1.0 + eps);
        assert!(distance.bounds().unwrap().hi() > 0.0);
    });
}

#[test]
fn g2_trim_removes_the_intersected_plane_and_preserves_empty_and_unwitnessed_outcomes() {
    let surface = plane();
    let source = surface.admit().unwrap();
    let trim = TrimmedPatch::new(vec![rectangle([2, 0], [4, 4], 4)]);
    let endpoints = [[0.125, 0.5, 0.0], [0.25, 0.5, 0.0]];
    with_cx(false, |cx| {
        let mut limits = budget(32);
        limits.point_budget.max_splits = 256;
        let distance = nonempty(
            source
                .segment_distance_enclosure_with_cx(
                    endpoints,
                    Some(trim.admit().unwrap()),
                    0.04,
                    limits,
                    cx,
                )
                .unwrap(),
        );
        let bounds = distance.bounds().unwrap();
        assert!(bounds.contains(0.25), "{bounds:?}");
        assert!(
            bounds.lo() > 0.1 && bounds.hi() < 0.4,
            "useful kept clearance: {bounds:?}"
        );
        assert!(std::ptr::eq(distance.trim().unwrap(), &trim));
        assert!(distance.witness().unwrap().surface_parameter[0] > 0.5);
        assert_eq!(
            distance.swept_ball_clearance(0.1).unwrap(),
            SweptBallClearance::Separated
        );

        let empty = TrimmedPatch::new(vec![]);
        let SegmentDistanceRun::Empty { distance } = source
            .segment_distance_enclosure_with_cx(
                endpoints,
                Some(empty.admit().unwrap()),
                0.01,
                budget(0),
                cx,
            )
            .unwrap()
        else {
            panic!("zero winding must prove an empty target");
        };
        assert!(distance.bounds().is_none() && distance.witness().is_none());
        assert_eq!(
            distance.swept_ball_clearance(1.0).unwrap(),
            SweptBallClearance::EmptyTarget
        );

        let island = TrimmedPatch::new(vec![rectangle([0, 0], [1, 1], 8)]);
        let distance = nonempty(
            source
                .segment_distance_enclosure_with_cx(
                    endpoints,
                    Some(island.admit().unwrap()),
                    0.01,
                    budget(0),
                    cx,
                )
                .unwrap(),
        );
        assert_eq!(distance.bounds().unwrap().hi(), f64::INFINITY);
        assert!(distance.witness().is_none());
        assert_eq!(
            distance.swept_ball_clearance(0.1).unwrap(),
            SweptBallClearance::Undetermined
        );
    });
}

#[test]
fn g4_shared_work_exhaustion_after_one_child_retains_the_whole_parent() {
    let surface = plane();
    let source = surface.admit().unwrap();
    let endpoints = [[0.25, 0.5, 1.0], [0.75, 0.5, 1.0]];
    with_cx(false, |cx| {
        let mut limits = budget(4);
        let seed = nonempty(
            source
                .segment_distance_enclosure_with_cx(endpoints, None, 10.0, limits, cx)
                .unwrap(),
        );
        limits.max_work_units = 2 * seed.work().charged_work_units;
        let SegmentDistanceRun::Unresolved {
            distance,
            stop: SegmentDistanceStop::WorkLimit,
        } = source
            .segment_distance_enclosure_with_cx(endpoints, None, 1e-9, limits, cx)
            .unwrap()
        else {
            panic!("shared work must stop before the second child");
        };
        assert_eq!(distance.work().point_queries, 2);
        assert_eq!(distance.work().splits, 0);
        assert_eq!(distance.work().remaining_cells, 1);
        assert!(distance.bounds().unwrap().contains(1.0));
        assert!(distance.work().charged_work_units <= limits.max_work_units);
    });
}

#[test]
fn g4_request_refusals_and_cancellation_publish_no_clearance() {
    let surface = plane();
    let source = surface.admit().unwrap();
    let endpoints = [[0.0, 0.0, 1.0], [1.0, 1.0, 1.0]];
    with_cx(false, |cx| {
        let mut no_memory = budget(1);
        no_memory.max_retained_bytes = 1;
        assert!(
            source
                .segment_distance_enclosure_with_cx(endpoints, None, 0.1, no_memory, cx)
                .is_err()
        );
        let mut no_work = budget(0);
        no_work.max_work_units = 1;
        assert!(
            source
                .segment_distance_enclosure_with_cx(endpoints, None, 0.1, no_work, cx)
                .is_err()
        );
        assert!(
            source
                .segment_distance_enclosure_with_cx([[f64::NAN; 3]; 2], None, 0.1, budget(0), cx)
                .is_err()
        );
        assert!(
            source
                .segment_distance_enclosure_with_cx(endpoints, None, 0.0, budget(0), cx)
                .is_err()
        );
        let distance = nonempty(
            source
                .segment_distance_enclosure_with_cx(endpoints, None, 0.1, budget(0), cx)
                .unwrap(),
        );
        for radius in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(distance.swept_ball_clearance(radius).is_err());
        }
    });
    with_cx(true, |cx| {
        assert!(matches!(
            source
                .segment_distance_enclosure_with_cx(endpoints, None, 0.1, budget(0), cx,)
                .unwrap(),
            SegmentDistanceRun::Cancelled
        ))
    });
}
