//! Analytic kept-region distance and adversarial trim/source boundaries.
#![cfg(feature = "nurbs-sdf")]

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_rep_nurbs::{
    KnotVector, NurbsCurve, NurbsSurface, Rat, SurfaceDistanceBudget, TrimLoop, TrimmedPatch,
    TrimmedSurfaceDistance, TrimmedSurfaceDistanceRun, TrimmedSurfaceDistanceStop,
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
                seed: 31,
                kernel_id: 59,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn budget(splits: u32) -> SurfaceDistanceBudget {
    SurfaceDistanceBudget {
        max_splits: splits,
        max_work_units: 256_000_000,
        max_retained_bytes: 96 * 1024 * 1024,
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

fn rectangle(min: [Rat; 2], max: [Rat; 2]) -> TrimLoop {
    let k = KnotVector::new([0, 0, 1, 2, 3, 4, 4].into_iter().map(Rat::int).collect(), 1).unwrap();
    TrimLoop::new(
        NurbsCurve::new(
            k,
            &[min, [max[0], min[1]], max, [min[0], max[1]], min],
            &[Rat::int(1); 5],
        )
        .unwrap(),
    )
    .unwrap()
}

fn fractional_rectangle(lo: [i128; 2], hi: [i128; 2], den: i128) -> TrimLoop {
    rectangle(lo.map(|n| Rat::new(n, den)), hi.map(|n| Rat::new(n, den)))
}

fn nonempty<'s, 't>(run: TrimmedSurfaceDistanceRun<'s, 't>) -> TrimmedSurfaceDistance<'s, 't> {
    match run {
        TrimmedSurfaceDistanceRun::Complete { distance }
        | TrimmedSurfaceDistanceRun::Unresolved { distance, .. } => distance,
        other => panic!("expected a nonempty-or-unresolved query, got {other:?}"),
    }
}

#[test]
fn g1_hole_removes_the_untrimmed_minimum_and_bounds_converge_from_kept_witnesses() {
    let surface = plane();
    let trim = TrimmedPatch::new(vec![
        fractional_rectangle([0, 0], [4, 4], 4),
        fractional_rectangle([1, 1], [3, 3], 4)
            .reversed_for_hole()
            .unwrap(),
    ]);
    with_cx(false, |cx| {
        let mut previous = fs_ivl::Interval::new(0.0, f64::INFINITY);
        for splits in [0, 64, 1024] {
            let distance = nonempty(
                surface
                    .admit()
                    .unwrap()
                    .trimmed_distance_enclosure_with_cx(
                        trim.admit().unwrap(),
                        [0.5, 0.5, 0.0],
                        0.01,
                        budget(splits),
                        cx,
                    )
                    .unwrap(),
            );
            let bounds = distance.bounds().unwrap();
            assert!(bounds.contains(0.25), "{splits}: {bounds:?}");
            assert!(bounds.lo() >= previous.lo() && bounds.hi() <= previous.hi());
            assert!(std::ptr::eq(distance.source(), &surface));
            assert!(std::ptr::eq(distance.trim(), &trim));
            assert_eq!(distance.query(), [0.5, 0.5, 0.0]);
            assert!(distance.classification_refusal().is_none());
            let work = distance.work();
            assert_eq!(
                work.search.original_span_cells + work.search.splits as usize,
                work.search.remaining_cells + work.search.excluded_cells + work.outside_cells,
                "every leaf must remain covered or have an exclusion reason"
            );
            if let Some([u, v]) = distance.witness_param() {
                assert!(u < 0.25 || u > 0.75 || v < 0.25 || v > 0.75);
            } else {
                assert_eq!(bounds.hi(), f64::INFINITY);
            }
            assert!(distance.work().search.charged_work_units <= budget(splits).max_work_units);
            assert!(
                distance.work().search.retained_payload_bytes <= budget(splits).max_retained_bytes
            );
            previous = bounds;
        }
        assert!(
            previous.width() < 0.04,
            "useful trim-boundary refinement: {previous:?}"
        );
    });
}

#[test]
fn g2_disconnected_islands_need_no_kept_seed_and_replay_deterministically() {
    let surface = plane();
    let trim = TrimmedPatch::new(vec![
        fractional_rectangle([0, 3], [1, 5], 8),
        fractional_rectangle([7, 3], [8, 5], 8),
    ]);
    with_cx(false, |cx| {
        let run = |splits| {
            surface
                .admit()
                .unwrap()
                .trimmed_distance_enclosure_with_cx(
                    trim.admit().unwrap(),
                    [0.5, 0.5, 0.0],
                    0.01,
                    budget(splits),
                    cx,
                )
                .unwrap()
        };
        let TrimmedSurfaceDistanceRun::Unresolved { distance: zero, .. } = run(0) else {
            panic!("the midpoint lies between the kept components")
        };
        assert!(zero.witness_param().is_none());
        assert_eq!(zero.bounds().unwrap().hi(), f64::INFINITY);
        let a = nonempty(run(1024));
        let b = nonempty(run(1024));
        assert!(a.bounds().unwrap().contains(0.375));
        assert!(a.bounds().unwrap().width() < 0.04);
        assert!(a.work().outside_cells > 0);
        assert_eq!(a.bounds(), b.bounds());
        assert_eq!(a.witness_param(), b.witness_param());
        assert_eq!(a.witness_position(), b.witness_position());
        assert_eq!(a.work(), b.work());
    });
}

#[test]
fn g1_rational_curved_trim_boundary_is_approached_only_from_inside() {
    // The rational quadratic arch has top exactly (1/2,1), with middle
    // Cartesian control (1/2,3/2) and weight 2. A quadratic bottom closes it.
    let k = KnotVector::new(
        [0, 0, 0, 1, 1, 2, 2, 2].into_iter().map(Rat::int).collect(),
        2,
    )
    .unwrap();
    let r = |n, d| Rat::new(n, d);
    let curve = NurbsCurve::new(
        k,
        &[
            [r(0, 1), r(0, 1)],
            [r(1, 2), r(3, 2)],
            [r(1, 1), r(0, 1)],
            [r(1, 2), r(0, 1)],
            [r(0, 1), r(0, 1)],
        ],
        &[r(1, 1), r(2, 1), r(1, 1), r(1, 1), r(1, 1)],
    )
    .unwrap();
    let trim = TrimmedPatch::new(vec![TrimLoop::new(curve).unwrap()]);
    let surface = plane();
    with_cx(false, |cx| {
        let TrimmedSurfaceDistanceRun::Complete { distance } = surface
            .admit()
            .unwrap()
            .trimmed_distance_enclosure_with_cx(
                trim.admit().unwrap(),
                [0.5, 1.0, 0.0],
                0.05,
                budget(1024),
                cx,
            )
            .unwrap()
        else {
            panic!("the rational arch boundary should close to the requested width")
        };
        assert!(distance.bounds().unwrap().contains(0.0));
        assert!(distance.bounds().unwrap().width() <= 0.05);
        assert!(distance.witness_param().unwrap()[1] < 1.0);
    });
}

#[test]
fn g2_thinner_than_legacy_grid_trim_uses_the_actual_binary_witness() {
    let den = 1i128 << 28;
    let lo = 0.5;
    let hi = 0.5 + 2.0f64.powi(-28);
    let ku = KnotVector::new(vec![0.0, 0.0, lo, hi, 1.0, 1.0], 1).unwrap();
    let kv = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap();
    let points: Vec<_> = [0.0, lo, hi, 1.0]
        .into_iter()
        .map(|u| vec![[u, 0.0, 0.0], [u, 1.0, 0.0]])
        .collect();
    let surface = NurbsSurface::new(ku, kv, &points, &vec![vec![1.0; 2]; 4]).unwrap();
    let trim = TrimmedPatch::new(vec![fractional_rectangle(
        [den / 2, 0],
        [den / 2 + 1, den],
        den,
    )]);
    with_cx(false, |cx| {
        let TrimmedSurfaceDistanceRun::Complete { distance } = surface
            .admit()
            .unwrap()
            .trimmed_distance_enclosure_with_cx(
                trim.admit().unwrap(),
                [f64::midpoint(lo, hi), 0.5, 1.0],
                1e-10,
                budget(0),
                cx,
            )
            .unwrap()
        else {
            panic!("an exact kept binary midpoint must survive a sub-grid trim")
        };
        assert!(distance.bounds().unwrap().contains(1.0));
        let [u, v] = distance.witness_param().unwrap();
        assert!(u > lo && u < hi && v > 0.0 && v < 1.0);
        assert_eq!(distance.work().search.original_span_cells, 3);
        assert_eq!(distance.work().search.splits, 0);
    });
}

#[test]
fn g2_full_surface_break_never_turns_an_unkept_left_limit_into_a_witness() {
    let ku = KnotVector::new(vec![0.0, 0.0, 0.5, 0.5, 1.0, 1.0], 1).unwrap();
    let kv = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap();
    let points: Vec<_> = [0.0, 1.0, 10.0, 11.0]
        .into_iter()
        .map(|x| vec![[x, 0.0, 0.0], [x, 1.0, 0.0]])
        .collect();
    let surface = NurbsSurface::new(ku, kv, &points, &vec![vec![1.0; 2]; 4]).unwrap();
    let trim = TrimmedPatch::new(vec![fractional_rectangle([1, 0], [2, 2], 2)]);
    with_cx(false, |cx| {
        let TrimmedSurfaceDistanceRun::Unresolved { distance, .. } = surface
            .admit()
            .unwrap()
            .trimmed_distance_enclosure_with_cx(
                trim.admit().unwrap(),
                [1.0, 0.5, 0.0],
                0.01,
                budget(128),
                cx,
            )
            .unwrap()
        else {
            panic!("a retained unkept span limit must not manufacture width closure")
        };
        assert!(distance.bounds().unwrap().contains(9.0));
        assert!(distance.bounds().unwrap().hi() >= 9.0);
        assert!(distance.witness_param().unwrap()[0] > 0.5);
        assert_eq!(distance.work().search.original_span_cells, 2);
    });
}

#[test]
fn g4_empty_unwitnessed_refused_and_cancelled_results_remain_distinct() {
    let surface = plane();
    let outside = TrimmedPatch::new(vec![fractional_rectangle([2, 2], [3, 3], 1)]);
    let kept = TrimmedPatch::new(vec![fractional_rectangle([3, 0], [4, 4], 4)]);
    with_cx(false, |cx| {
        let TrimmedSurfaceDistanceRun::Empty { distance } = surface
            .admit()
            .unwrap()
            .trimmed_distance_enclosure_with_cx(
                outside.admit().unwrap(),
                [0.5, 0.5, 1.0],
                1e-5,
                budget(0),
                cx,
            )
            .unwrap()
        else {
            panic!("all original cells are exactly outside")
        };
        assert!(distance.bounds().is_none() && distance.witness_param().is_none());
        assert_eq!(distance.work().outside_cells, 1);
        assert_eq!(distance.work().search.remaining_cells, 0);

        let small_work = SurfaceDistanceBudget {
            max_work_units: 2000,
            ..budget(8)
        };
        let TrimmedSurfaceDistanceRun::Unresolved {
            distance,
            stop: TrimmedSurfaceDistanceStop::ClassificationRefused,
        } = surface
            .admit()
            .unwrap()
            .trimmed_distance_enclosure_with_cx(
                kept.admit().unwrap(),
                [0.5, 0.5, 1.0],
                1e-5,
                small_work,
                cx,
            )
            .unwrap()
        else {
            panic!("classifier exhaustion must retain the complete source coverage")
        };
        assert!(distance.classification_refusal().is_some());
        assert!(distance.witness_param().is_none());
        assert_eq!(distance.bounds().unwrap().hi(), f64::INFINITY);
        assert_eq!(distance.work().search.original_span_cells, 1);
        assert_eq!(distance.work().search.remaining_cells, 1);
        assert!(distance.work().search.charged_work_units <= 2000);
        assert!(
            surface
                .admit()
                .unwrap()
                .trimmed_distance_enclosure_with_cx(
                    kept.admit().unwrap(),
                    [0.5, 0.5, 1.0],
                    1e-5,
                    SurfaceDistanceBudget {
                        max_retained_bytes: 1024,
                        ..budget(0)
                    },
                    cx
                )
                .is_err()
        );
    });
    with_cx(true, |cx| {
        assert!(matches!(
            surface
                .admit()
                .unwrap()
                .trimmed_distance_enclosure_with_cx(
                    kept.admit().unwrap(),
                    [0.5, 0.5, 1.0],
                    1e-5,
                    budget(8),
                    cx
                )
                .unwrap(),
            TrimmedSurfaceDistanceRun::Cancelled
        ));
        assert!(
            surface
                .admit()
                .unwrap()
                .trimmed_distance_enclosure_with_cx(
                    kept.admit().unwrap(),
                    [f64::NAN, 0.0, 0.0],
                    1e-5,
                    budget(8),
                    cx
                )
                .is_err()
        );
    });
}

#[test]
fn g4_unrepresentable_parameter_membership_retains_coverage_without_a_witness() {
    let k = KnotVector::new(vec![1e100, 1e100, 2e100, 2e100], 1).unwrap();
    let surface = NurbsSurface::new(
        k.clone(),
        k,
        &vec![vec![[0.0; 3]; 2]; 2],
        &vec![vec![1.0; 2]; 2],
    )
    .unwrap();
    let trim = TrimmedPatch::new(vec![fractional_rectangle([0, 0], [1, 1], 1)]);
    with_cx(false, |cx| {
        let TrimmedSurfaceDistanceRun::Unresolved { distance, .. } = surface
            .admit()
            .unwrap()
            .trimmed_distance_enclosure_with_cx(
                trim.admit().unwrap(),
                [1.0, 0.0, 0.0],
                1e-4,
                budget(0),
                cx,
            )
            .unwrap()
        else {
            panic!("unsupported parameter conversion proves neither membership nor emptiness")
        };
        assert!(distance.witness_param().is_none());
        assert_eq!(distance.bounds().unwrap().hi(), f64::INFINITY);
        assert!(distance.work().unsupported_parameter_queries >= 2);
        assert_eq!(distance.work().search.original_span_cells, 1);
        assert_eq!(distance.work().search.remaining_cells, 1);
        assert!(distance.classification_refusal().is_none());
    });
}
