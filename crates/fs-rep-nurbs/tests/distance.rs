//! G1 analytic distance, G3 representation/scale, and G4 budget/cancellation
//! checks for the source-bound interval lane. This is separate from measured SDF.
#![cfg(feature = "nurbs-sdf")]

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_ivl::Interval;
use fs_rep_nurbs::{
    KnotVector, NurbsSurface, SurfaceDistanceBudget, SurfaceDistanceEnclosure, SurfaceDistanceRun,
    SurfaceDistanceStop,
};

fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    if cancelled {
        gate.request();
    }
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 17,
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
        max_work_units: 32_000_000,
        max_retained_bytes: 8 * 1024 * 1024,
    }
}

fn result(run: SurfaceDistanceRun<'_>) -> SurfaceDistanceEnclosure<'_> {
    match run {
        SurfaceDistanceRun::Complete { enclosure }
        | SurfaceDistanceRun::Unresolved { enclosure, .. } => enclosure,
        SurfaceDistanceRun::Cancelled => panic!("active operation cancelled"),
    }
}

fn plane(weights: [[f64; 2]; 2]) -> NurbsSurface<f64> {
    let knots = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap();
    NurbsSurface::new(
        knots.clone(),
        knots,
        &[
            vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
        ],
        &[weights[0].to_vec(), weights[1].to_vec()],
    )
    .unwrap()
}

/// The homogeneous polynomials are exactly (1-u², 2u, v(1+u²), 1+u²).
/// Integer binary controls therefore represent a TRUE unit quarter cylinder;
/// there is no irrational circle weight rounded before the enclosure query.
fn quarter_cylinder() -> NurbsSurface<f64> {
    let ku = KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2).unwrap();
    let kv = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap();
    NurbsSurface::from_homogeneous(
        ku,
        kv,
        vec![
            vec![[1.0, 0.0, 0.0, 1.0], [1.0, 0.0, 1.0, 1.0]],
            vec![[1.0, 1.0, 0.0, 1.0], [1.0, 1.0, 1.0, 1.0]],
            vec![[0.0, 2.0, 0.0, 2.0], [0.0, 2.0, 2.0, 2.0]],
        ],
    )
    .unwrap()
}

#[test]
fn g1_plane_zero_split_enclosure_and_rational_parameterization_close_width() {
    for weights in [[[1.0; 2]; 2], [[1.0, 2.0], [4.0, 8.0]]] {
        let surface = plane(weights);
        let admitted = surface.admit().unwrap();
        with_cx(false, |cx| {
            let q = [0.25, 0.75, 1.0];
            let zero = result(
                admitted
                    .distance_enclosure_with_cx(q, 1e-4, budget(0), cx)
                    .unwrap(),
            );
            assert!(zero.bounds().contains(1.0), "{:?}", zero.bounds());
            assert_eq!(zero.work().original_span_cells, 1);
            assert_eq!(zero.work().splits, 0);
            assert!(std::ptr::eq(zero.source(), &surface));
            assert_eq!(zero.query(), q);
            let SurfaceDistanceRun::Complete { enclosure } = admitted
                .distance_enclosure_with_cx(q, 1e-4, budget(2048), cx)
                .unwrap()
            else {
                panic!("planar distance must close the admitted width")
            };
            assert!(enclosure.bounds().contains(1.0));
            assert!(enclosure.bounds().width() <= enclosure.tolerance());
            assert!(enclosure.bounds().lo() >= zero.bounds().lo());
            assert!(enclosure.bounds().hi() <= zero.bounds().hi());
            let p = admitted
                .eval(enclosure.witness_param()[0], enclosure.witness_param()[1])
                .unwrap();
            // Ordinary evaluation is only a cross-check of the retained witness,
            // not the oracle authorizing its interval or the distance bounds.
            for (actual, interval) in p.into_iter().zip(enclosure.witness_position()) {
                assert!(actual >= interval.lo() - 2e-14 && actual <= interval.hi() + 2e-14);
            }
        });
    }
}

#[test]
fn g1_exact_rational_cylinder_keeps_analytic_distance_and_refines_monotonically() {
    let surface = quarter_cylinder();
    let admitted = surface.admit().unwrap();
    let truth = Interval::point(2.0).sqrt() - Interval::point(1.0);
    with_cx(false, |cx| {
        let mut previous = Interval::new(0.0, f64::INFINITY);
        for splits in [0, 16, 128, 1024] {
            let run = admitted
                .distance_enclosure_with_cx([1.0, 1.0, 0.5], 1e-6, budget(splits), cx)
                .unwrap();
            let enclosure = result(run);
            assert!(
                enclosure.bounds().encloses(truth),
                "{splits}: {:?} vs {truth:?}",
                enclosure.bounds()
            );
            assert!(enclosure.bounds().lo() >= previous.lo());
            assert!(enclosure.bounds().hi() <= previous.hi());
            assert!(enclosure.work().splits <= splits);
            assert!(enclosure.work().charged_work_units <= budget(splits).max_work_units);
            assert!(enclosure.work().retained_payload_bytes <= budget(splits).max_retained_bytes);
            previous = enclosure.bounds();
        }
        assert!(
            previous.width() < 0.04,
            "refinement should be useful: {previous:?}"
        );
        let a = result(
            admitted
                .distance_enclosure_with_cx([1.0, 1.0, 0.5], 1e-6, budget(128), cx)
                .unwrap(),
        );
        let b = result(
            admitted
                .distance_enclosure_with_cx([1.0, 1.0, 0.5], 1e-6, budget(128), cx)
                .unwrap(),
        );
        assert_eq!(a.bounds(), b.bounds());
        assert_eq!(a.witness_param(), b.witness_param());
        assert_eq!(a.witness_position(), b.witness_position());
        assert_eq!(a.work(), b.work());
    });
}

#[test]
fn g3_repeated_knots_cover_both_components_and_never_publish_a_left_limit_witness() {
    let ku = KnotVector::new(vec![0.0, 0.0, 0.5, 0.5, 1.0, 1.0], 1).unwrap();
    let kv = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap();
    let points: Vec<_> = [0.0, 1.0, 10.0, 11.0]
        .into_iter()
        .map(|x| vec![[x, 0.0, 0.0], [x, 1.0, 0.0]])
        .collect();
    let surface = NurbsSurface::new(ku, kv, &points, &vec![vec![1.0; 2]; 4]).unwrap();
    let admitted = surface.admit().unwrap();
    with_cx(false, |cx| {
        for (q, distance, first_component) in
            [([5.0, 0.5, 0.0], 4.0, true), ([10.5, 0.5, 1.0], 1.0, false)]
        {
            let enclosure = result(
                admitted
                    .distance_enclosure_with_cx(q, 1e-5, budget(1024), cx)
                    .unwrap(),
            );
            assert_eq!(enclosure.work().original_span_cells, 2);
            assert!(
                enclosure.bounds().contains(distance),
                "{:?}",
                enclosure.bounds()
            );
            if first_component {
                assert!(enclosure.witness_param()[0] < 0.5);
                assert!(enclosure.bounds().width() < 0.05);
            } else {
                assert!(enclosure.witness_param()[0] >= 0.5);
            }
            assert!((0.0..=1.0).contains(&enclosure.witness_param()[1]));
        }
    });
}

#[test]
fn g3_point_degeneracy_extreme_weights_and_norm_scales_keep_outward_bounds() {
    let knots = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap();
    let surface = NurbsSurface::new(
        knots.clone(),
        knots,
        &vec![vec![[0.0; 3]; 2]; 2],
        &[
            vec![f64::MIN_POSITIVE, 1e308],
            vec![1e308, f64::MIN_POSITIVE],
        ],
    )
    .unwrap();
    with_cx(false, |cx| {
        for magnitude in [2.0f64.powi(-600), 1.0, 2.0f64.powi(600)] {
            let SurfaceDistanceRun::Complete { enclosure } = surface
                .admit()
                .unwrap()
                .distance_enclosure_with_cx([magnitude, 0.0, 0.0], magnitude * 1e-12, budget(0), cx)
                .unwrap()
            else {
                panic!("constant-point distance should need no refinement")
            };
            assert!(
                enclosure.bounds().contains(magnitude),
                "{:?}",
                enclosure.bounds()
            );
            assert!(enclosure.witness_position().iter().all(|i| i.contains(0.0)));
        }
    });
}

#[test]
fn g4_complete_coverage_survives_refinement_limits_and_pre_cancellation_publishes_nothing() {
    let surface = plane([[1.0; 2]; 2]);
    let admitted = surface.admit().unwrap();
    let q = [0.13, 0.22, 1.0];
    with_cx(false, |cx| {
        let SurfaceDistanceRun::Unresolved {
            enclosure: zero,
            stop: SurfaceDistanceStop::SplitLimit,
        } = admitted
            .distance_enclosure_with_cx(q, 1e-12, budget(0), cx)
            .unwrap()
        else {
            panic!("zero splits must report unresolved width")
        };
        let memory = SurfaceDistanceBudget {
            max_retained_bytes: zero.work().retained_payload_bytes,
            ..budget(64)
        };
        let SurfaceDistanceRun::Unresolved {
            enclosure,
            stop: SurfaceDistanceStop::MemoryLimit,
        } = admitted
            .distance_enclosure_with_cx(q, 1e-12, memory, cx)
            .unwrap()
        else {
            panic!("limited frontier must keep complete initial coverage")
        };
        assert!(enclosure.bounds().contains(1.0));
        // A borrowed source retains the actual allocations supplied by its
        // constructor, including spare Vec capacity outside the control count.
        let knots = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap();
        let mut padded_net = Vec::with_capacity(32);
        for x in [0.0, 1.0] {
            let mut row = Vec::with_capacity(4096);
            row.extend([[x, 0.0, 0.0, 1.0], [x, 1.0, 0.0, 1.0]]);
            padded_net.push(row);
        }
        let padded = NurbsSurface::from_homogeneous(knots.clone(), knots, padded_net).unwrap();
        assert!(
            padded
                .admit()
                .unwrap()
                .distance_enclosure_with_cx(q, 1e-12, memory, cx)
                .is_err(),
            "source storage must be charged by retained capacity, not live length"
        );
        let initial_cost = result(
            admitted
                .distance_enclosure_with_cx([0.5, 0.5, 1.0], 1e-10, budget(64), cx)
                .unwrap(),
        )
        .work()
        .charged_work_units;
        let work = SurfaceDistanceBudget {
            max_work_units: initial_cost,
            ..budget(64)
        };
        let SurfaceDistanceRun::Unresolved {
            enclosure,
            stop: SurfaceDistanceStop::WorkLimit,
        } = admitted
            .distance_enclosure_with_cx(q, 1e-12, work, cx)
            .unwrap()
        else {
            panic!("work exhaustion must not assert requested-width closure")
        };
        assert!(enclosure.bounds().contains(1.0));
        assert_eq!(enclosure.work().charged_work_units, initial_cost);
        assert!(
            admitted
                .distance_enclosure_with_cx(
                    q,
                    1e-12,
                    SurfaceDistanceBudget {
                        max_work_units: 0,
                        ..budget(0)
                    },
                    cx
                )
                .is_err()
        );
        assert!(
            admitted
                .distance_enclosure_with_cx(
                    q,
                    1e-12,
                    SurfaceDistanceBudget {
                        max_retained_bytes: 1,
                        ..budget(0)
                    },
                    cx
                )
                .is_err()
        );
    });
    with_cx(true, |cx| {
        assert!(matches!(
            admitted
                .distance_enclosure_with_cx(q, 1e-12, budget(64), cx)
                .unwrap(),
            SurfaceDistanceRun::Cancelled
        ));
        assert!(
            admitted
                .distance_enclosure_with_cx([f64::NAN, 0.0, 0.0], 1e-12, budget(0), cx)
                .is_err()
        );
        assert!(
            admitted
                .distance_enclosure_with_cx(q, 0.0, budget(0), cx)
                .is_err()
        );
    });
}

#[test]
fn g4_parameter_resolution_and_finite_range_refusals_are_explicit() {
    let hi = fs_math::next_up(1.0);
    let knots = KnotVector::new(vec![1.0, 1.0, hi, hi], 1).unwrap();
    let surface = NurbsSurface::new(
        knots.clone(),
        knots,
        &[
            vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
        ],
        &vec![vec![1.0; 2]; 2],
    )
    .unwrap();
    with_cx(false, |cx| {
        let SurfaceDistanceRun::Unresolved {
            enclosure,
            stop: SurfaceDistanceStop::ParameterResolution,
        } = surface
            .admit()
            .unwrap()
            .distance_enclosure_with_cx([0.25, 0.25, 0.5], 1e-12, budget(64), cx)
            .unwrap()
        else {
            panic!("no representable interior split must not close the width")
        };
        assert!(enclosure.bounds().contains(0.5));
        assert_eq!(enclosure.witness_param(), [1.0, 1.0]);
        let knots = KnotVector::new(vec![0.0, 1.0], 0).unwrap();
        let extreme = NurbsSurface::from_homogeneous(
            knots.clone(),
            knots,
            vec![vec![[f64::MAX, 0.0, 0.0, 1.0]]],
        )
        .unwrap();
        assert!(
            extreme
                .admit()
                .unwrap()
                .distance_enclosure_with_cx([0.0; 3], 1.0, budget(0), cx)
                .is_err()
        );
    });
}
