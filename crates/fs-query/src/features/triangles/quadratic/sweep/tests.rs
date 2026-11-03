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
                    seed: 7,
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

const FLAT: [[f64; 3]; 6] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.5, 0.0, 0.0],
    [0.5, 0.5, 0.0],
    [0.0, 0.5, 0.0],
];

fn paraboloid(scale: f64) -> [[f64; 3]; 6] {
    FLAT.map(|[x, y, _]| [x, y, scale * (x * x + y * y)])
}

fn controls() -> QuadraticSweepControls {
    QuadraticSweepControls {
        time: LinearSweepControls {
            maximum_windows: 1024,
            fraction_tolerance: 1e-3,
        },
        distance: QuadraticQueryControls {
            distance_tolerance: 1e-3,
            witness_tolerance: 1.0,
            maximum_cells: 32768,
            maximum_active_cells: 8192,
        },
        maximum_distance_cells: 1 << 20,
    }
}

#[test]
fn affine_reduction_preserves_exact_touch_and_refuses_one_ulp_penetration() {
    run(|_, cx| {
        let mut limits = controls();
        limits.time.maximum_windows = 1;
        for height in [0.125, 0.25] {
            let report = linear_sphere_quadratic_triangle(
                [0.25, 0.25, 0.125],
                [0.25, 0.25, height],
                0.125,
                FLAT,
                FLAT,
                limits,
                cx,
            )
            .unwrap();
            assert_eq!(
                report.path.verdict,
                LinearSweepVerdict::Clear {
                    minimum_clearance: 0.0
                }
            );
            assert_eq!(report.path.examined_windows, 1);
            assert_eq!(report.examined_distance_cells, 0);
        }
        let penetrating = linear_sphere_quadratic_triangle(
            [0.25, 0.25, 0.125_f64.next_down()],
            [0.25, 0.25, 0.25],
            0.125,
            FLAT,
            FLAT,
            controls(),
            cx,
        )
        .unwrap();
        assert!(matches!(
            penetrating.path.verdict,
            LinearSweepVerdict::PossibleContact {
                fraction_start: 0.0,
                ..
            }
        ));
    });
}

#[test]
fn moving_midside_nodes_cannot_hide_behind_stationary_corners() {
    run(|_, cx| {
        // At fraction 1/2 the bowed face contains (.25,.25,1.25).
        // Both endpoint patches are separated from this sphere, while the
        // corner triangle remains at z=0 throughout the motion.
        let mut bowed = FLAT;
        for node in &mut bowed[3..] {
            node[2] = 2.0;
        }
        let center = [0.25, 0.25, 1.25];
        let linear = super::super::super::sweep::linear_sphere_triangle(
            center,
            center,
            0.0625,
            [FLAT[0], FLAT[1], FLAT[2]],
            [bowed[0], bowed[1], bowed[2]],
            controls().time,
            cx,
        )
        .unwrap();
        assert!(matches!(linear.verdict, LinearSweepVerdict::Clear { .. }));
        let curved = linear_sphere_quadratic_triangle(
            center,
            center,
            0.0625,
            FLAT,
            bowed,
            controls(),
            cx,
        )
        .unwrap();
        match curved.path.verdict {
            LinearSweepVerdict::PossibleContact {
                fraction_start,
                fraction_end,
            } => {
                assert!(fraction_start > 0.0 && fraction_end <= 0.5);
                assert!(fraction_end - fraction_start <= controls().time.fraction_tolerance);
            }
            other => panic!("missed midside crossing: {other:?}"),
        }
        assert!(curved.examined_distance_cells > 0);
    });
}

#[test]
fn curved_distance_bound_proves_clearance_inside_the_control_hull() {
    run(|_, cx| {
        // S(u,v)=(u,v,u²+v²), P=S(1/4,1/4), C=P+(1/8)(-1/2,-1/2,1).
        // |S-C|²-|P-C|² = (3/4)|delta_uv|²+(S_z-P_z)² >= 0.
        // Thus the exact closest distance is sqrt(3/2)/8. C lies inside the
        // original Bernstein hull, so root-hull planes alone cannot prove it.
        let analytic_clearance = 1.5_f64.sqrt() / 8.0 - 0.125;
        for origin in [
            [0.0; 3],
            [1073741824.0, -1073741824.0, 1073741824.0],
        ] {
            let start = paraboloid(1.0)
                .map(|node| core::array::from_fn(|axis| origin[axis] + node[axis]));
            let center =
                core::array::from_fn(|axis| origin[axis] + [0.1875, 0.1875, 0.25][axis]);
            for translation in [[0.0; 3], [4.0, -2.0, 1.0]] {
                let end = start
                    .map(|node| core::array::from_fn(|axis| node[axis] + translation[axis]));
                let center_end = core::array::from_fn(|axis| center[axis] + translation[axis]);
                let report = linear_sphere_quadratic_triangle(
                    center,
                    center_end,
                    0.125,
                    start,
                    end,
                    controls(),
                    cx,
                )
                .unwrap();
                match report.path.verdict {
                    LinearSweepVerdict::Clear { minimum_clearance } => {
                        assert!(minimum_clearance > 0.025);
                        assert!(minimum_clearance <= analytic_clearance);
                    }
                    other => panic!("failed curved common-translation clearance: {other:?}"),
                }
                assert_eq!(report.path.examined_windows, 1);
                assert!(report.examined_distance_cells > 0);
            }
        }
    });
}

#[test]
fn non_axis_curved_touch_remains_unresolved_without_a_proof() {
    run(|_, cx| {
        // S=(u,v,4(u²+v²)), P=S(1/4,1/4), C=P+(1/16)(-2,-2,1).
        // The normal has length 3, so the exact distance is 3/16.
        // |S-C|²-|P-C|²=(1/2)|delta_uv|²+(S_z-P_z)² proves globality.
        let patch = paraboloid(4.0);
        let mut limits = controls();
        limits.time.fraction_tolerance = 0.0625;
        for radius in [0.1875, 0.1875_f64.next_up()] {
            let report = linear_sphere_quadratic_triangle(
                [0.125, 0.125, 0.5625],
                [0.125, 0.125, 0.5625],
                radius,
                patch,
                patch,
                limits,
                cx,
            )
            .unwrap();
            assert_eq!(
                report.path.verdict,
                LinearSweepVerdict::PossibleContact {
                    fraction_start: 0.0,
                    fraction_end: 0.0625,
                }
            );
        }
    });
}

#[test]
fn time_distance_and_cancellation_limits_refuse_without_clearance() {
    run(|gate, cx| {
        let patch = paraboloid(1.0);
        let center = [0.1875, 0.1875, 0.25];
        let mut limits = controls();
        limits.maximum_distance_cells = 1;
        assert!(matches!(
            linear_sphere_quadratic_triangle(center, center, 0.125, patch, patch, limits, cx),
            Err(QuadraticQueryError::CellBudgetExhausted { .. })
                | Err(QuadraticQueryError::ActiveCellBudgetExhausted { .. })
        ));
        limits = controls();
        limits.time.maximum_windows = 1;
        assert!(matches!(
            linear_sphere_quadratic_triangle(
                [0.25, 0.25, 1.0],
                [0.25, 0.25, -1.0],
                0.125,
                FLAT,
                FLAT,
                limits,
                cx,
            ),
            Err(QuadraticQueryError::Geometry(
                TriangleQueryError::WorkBudgetExhausted { .. }
            ))
        ));
        gate.request();
        assert!(matches!(
            linear_sphere_quadratic_triangle(center, center, 0.125, patch, patch, controls(), cx),
            Err(QuadraticQueryError::Geometry(TriangleQueryError::Cancelled))
        ));
    });
}
