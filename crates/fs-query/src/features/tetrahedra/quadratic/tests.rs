use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

fn run(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &gate,
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 17,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
        );
    });
}

fn nodes(map: impl Fn([f64; 3]) -> [f64; 3]) -> [[f64; 3]; 10] {
    // Independently specify the native 01,12,02,03,13,23 midside order.
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.5, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.0],
        [0.0, 0.0, 0.5],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ]
    .map(map)
}

fn lower(report: QuadraticTetPathReport) -> f64 {
    match report.verdict {
        QuadraticTetPathVerdict::Positive {
            minimum_scaled_determinant,
        } => minimum_scaled_determinant,
        other => panic!("positive closed-domain path required, got {other:?}"),
    }
}

#[test]
fn isochoric_curved_shear_keeps_every_material_point_positive() {
    run(|_, cx| {
        let start = nodes(|p| p);
        let end = nodes(|[x, y, z]| [x - 8.0 * y * z, y, z]);
        assert_eq!(start[..4], end[..4]);
        let report = quadratic_tetrahedron_orientation(start, end, Default::default(), cx).unwrap();
        let bound = lower(report);
        let exact = report.normalization_length.powi(-3);
        assert!(bound > 0.99 * exact && bound <= exact);
        assert_eq!(report.examined_cells, 1);
    });
}

#[test]
fn bernstein_dependencies_admit_correlated_curvature_without_subdivision() {
    run(|_, cx| {
        // Output rotation/scaling of the isochoric shear X -> X + 4tY².
        // det(J)=2 identically; independent entry ranges include negatives.
        let start = nodes(|[x, y, z]| [x + y, -x + y, z]);
        let end = nodes(|[x, y, z]| [x + y + 4.0 * y * y, -x + y - 4.0 * y * y, z]);
        let (jacobians, _) = normalized_jacobians(start, end).unwrap();
        assert!(determinant(&matrix_hull(&jacobians)).lo() < 0.0);
        let report = quadratic_tetrahedron_orientation(start, end, Default::default(), cx).unwrap();
        let exact = 2.0 * report.normalization_length.powi(-3);
        let bound = lower(report);
        assert!(bound > 0.99 * exact && bound <= exact);
        assert_eq!(report.examined_cells, 1);
    });
}

#[test]
fn unchanged_positive_corner_tetrahedron_does_not_hide_a_quadratic_fold() {
    run(|_, cx| {
        let straight = nodes(|p| p);
        let folded = nodes(|[x, y, z]| [x + 1.5 * x * (1.0 - x), y, z]);
        assert_eq!(straight[..4], folded[..4]);
        let report =
            quadratic_tetrahedron_orientation(folded, folded, Default::default(), cx).unwrap();
        let QuadraticTetPathVerdict::Unresolved {
            fraction_start,
            fraction_end,
            barycentric_bounds,
            scaled_determinant_bounds,
        } = report.verdict
        else {
            panic!("quadratic map reverses orientation next to corner 1");
        };
        assert_eq!([fraction_start, fraction_end], [0.0, 0.0]);
        assert_eq!(barycentric_bounds[1], [1.0, 1.0]);
        assert!(scaled_determinant_bounds[1] < 0.0);
        let exact = -0.5 * report.normalization_length.powi(-3);
        assert!(scaled_determinant_bounds[0] <= exact);
        assert!(scaled_determinant_bounds[1] >= exact);
    });
}

#[test]
fn positive_endpoints_do_not_admit_an_intermediate_collapse() {
    run(|_, cx| {
        let start = nodes(|p| p);
        let end = nodes(|[x, y, z]| [-x, -y, z]);
        let controls = QuadraticTetPathControls {
            minimum_scaled_determinant: 0.0,
            ..Default::default()
        };
        assert!(
            lower(quadratic_tetrahedron_orientation(start, start, controls, cx).unwrap()) > 0.0
        );
        assert!(lower(quadratic_tetrahedron_orientation(end, end, controls, cx).unwrap()) > 0.0);
        let report = quadratic_tetrahedron_orientation(start, end, controls, cx).unwrap();
        let QuadraticTetPathVerdict::Unresolved {
            fraction_start,
            fraction_end,
            scaled_determinant_bounds,
            ..
        } = report.verdict
        else {
            panic!("det(J)=(1-2t)^2 vanishes at t=1/2");
        };
        assert!(fraction_start <= 0.5 && fraction_end >= 0.5);
        assert!(fraction_end - fraction_start <= controls.fraction_tolerance);
        assert!(scaled_determinant_bounds[0] <= 0.0);
    });
}

#[test]
fn translation_and_dyadic_size_do_not_change_the_normalized_geometry() {
    run(|_, cx| {
        let start = nodes(|p| p);
        let end = nodes(|[x, y, z]| [x - 8.0 * y * z, y, z]);
        let baseline =
            quadratic_tetrahedron_orientation(start, end, Default::default(), cx).unwrap();
        for scale in [2.0_f64.powi(-500), 2.0_f64.powi(500)] {
            let scaled = quadratic_tetrahedron_orientation(
                start.map(|p| p.map(|x| x * scale)),
                end.map(|p| p.map(|x| x * scale)),
                Default::default(),
                cx,
            )
            .unwrap();
            assert_eq!(
                scaled.normalization_length,
                baseline.normalization_length * scale
            );
            assert!((lower(scaled) - lower(baseline)).abs() < 1.0e-12);
        }
        let translated = quadratic_tetrahedron_orientation(
            start.map(|p| [p[0] + 1_048_576.0, p[1] - 2_097_152.0, p[2] + 4_194_304.0]),
            end.map(|p| [p[0] + 1_048_576.0, p[1] - 2_097_152.0, p[2] + 4_194_304.0]),
            Default::default(),
            cx,
        )
        .unwrap();
        assert_eq!(
            translated.normalization_length,
            baseline.normalization_length
        );
        assert_eq!(lower(translated), lower(baseline));
        // An independently moving origin is removed at each endpoint too.
        let moving_origin = quadratic_tetrahedron_orientation(
            start,
            end.map(|p| [p[0] + 1024.0, p[1] - 2048.0, p[2] + 4096.0]),
            Default::default(),
            cx,
        )
        .unwrap();
        assert_eq!(lower(moving_origin), lower(baseline));
    });
}

#[test]
fn original_edge_order_and_all_quadratic_derivatives_are_retained() {
    let curved = nodes(|[x, y, z]| {
        [
            x + 0.75 * x * y - 0.5 * y * z,
            y + 0.25 * x * z - 0.5 * y * y,
            z + 0.5 * x * z + 0.125 * x * x,
        ]
    });
    let (jacobians, scale) = normalized_jacobians(curved, curved).unwrap();
    let lambda = [0.125, 0.25, 0.125, 0.5];
    let sample = restrict(
        &jacobians,
        &Cell {
            vertices: [lambda; 4],
            fraction: [0.0, 1.0],
            spatial_depth: 0,
            time_depth: 0,
            serial: 0,
        },
    );
    let [_, x, y, z] = lambda;
    let exact = [
        [1.0 + 0.75 * y, 0.75 * x - 0.5 * z, -0.5 * y],
        [0.25 * z, 1.0 - y, 0.25 * x],
        [0.5 * z + 0.25 * x, 0.0, 1.0 + 0.5 * x],
    ];
    for row in 0..3 {
        for column in 0..3 {
            assert!(sample[0][0][row][column].contains(exact[row][column] / scale));
        }
    }
}

#[test]
fn bounded_work_storage_cancellation_and_invalid_data_refuse_admission() {
    run(|gate, cx| {
        let start = nodes(|p| p);
        let end = nodes(|[x, y, z]| [-x, -y, z]);
        let limited = QuadraticTetPathControls {
            maximum_cells: 1,
            ..Default::default()
        };
        assert_eq!(
            quadratic_tetrahedron_orientation(start, end, limited, cx),
            Err(QuadraticTetPathError::CellBudgetExhausted { maximum_cells: 1 })
        );
        let limited = QuadraticTetPathControls {
            maximum_active_cells: 1,
            ..Default::default()
        };
        assert_eq!(
            quadratic_tetrahedron_orientation(start, end, limited, cx),
            Err(QuadraticTetPathError::ActiveCellBudgetExhausted {
                maximum_active_cells: 1
            })
        );
        assert!(matches!(
            quadratic_tetrahedron_orientation([[f64::NAN; 3]; 10], end, Default::default(), cx),
            Err(QuadraticTetPathError::Invalid(_))
        ));
        assert_eq!(
            quadratic_tetrahedron_orientation(
                [[0.0; 3]; 10],
                [[0.0; 3]; 10],
                Default::default(),
                cx
            ),
            Err(QuadraticTetPathError::DegenerateGeometry)
        );
        gate.request();
        assert_eq!(
            quadratic_tetrahedron_orientation(start, end, Default::default(), cx),
            Err(QuadraticTetPathError::Cancelled)
        );
    });
}

#[test]
fn power_of_two_normalization_handles_subnormal_and_overflow_limits() {
    assert_eq!(
        power_of_two_ceiling(f64::from_bits(3)).unwrap(),
        f64::from_bits(4)
    );
    assert_eq!(power_of_two_ceiling(1.0).unwrap(), 1.0);
    assert_eq!(power_of_two_ceiling(1.0_f64.next_up()).unwrap(), 2.0);
    assert_eq!(
        power_of_two_ceiling(f64::MAX),
        Err(QuadraticTetPathError::NumericalRange)
    );
}
