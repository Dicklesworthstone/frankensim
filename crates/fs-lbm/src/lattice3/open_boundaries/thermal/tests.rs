use super::*;
use crate::d3q19::{BoundarySpec3, CollisionModel3};
use crate::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Equilibrium27};
use crate::lattice3::{D3q19, Grid3, GridBudget3};
use crate::thermal::boundary::ThermalWall;
use fs_exec::CancelGate;

mod tensor;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 0x0007_03A1,
                kernel_id: 227,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        f(&cx)
    })
}
fn limits(maximum: f64) -> ThermalOpenLimits3 {
    ThermalOpenLimits3 {
        maximum_cells: 512,
        maximum_storage_bytes: 1024 * 1024,
        temperature_range: [-2.0, maximum],
        balance_tolerance: 1e-10,
    }
}
fn thermal_faces(axis: usize, temperature: f64) -> ThermalOpenBoundaries3 {
    let mut faces = [ThermalOpenFace3::Periodic; 6];
    faces[2 * axis] = ThermalOpenFace3::Reservoir { temperature };
    faces[2 * axis + 1] = ThermalOpenFace3::Reservoir { temperature };
    ThermalOpenBoundaries3 {
        faces,
        internal_wall: ThermalWall::Adiabatic,
    }
}
fn native_faces(axis: usize, velocity: [f64; 3]) -> BoundarySpec3 {
    let mut faces = [FaceBoundary3::Periodic; 6];
    faces[2 * axis] = FaceBoundary3::Velocity { velocity };
    faces[2 * axis + 1] = FaceBoundary3::Pressure { density: 1.2 };
    BoundarySpec3::new(faces)
}
fn fixture<L: EquilibriumLattice3<Q>, const Q: usize>(
    cx: &Cx<'_>,
    collision: L::Collision,
    axis: usize,
    velocity: [f64; 3],
    force: [f64; 3],
) -> OpenBoundaryGrid3<L, Q> {
    let mut periodic = [true; 3];
    periodic[axis] = false;
    let grid = Grid3::try_new_equilibrium(
        [5, 4, 3],
        collision,
        force,
        periodic,
        GridBudget3 {
            maximum_cells: 512,
            maximum_mach: 0.2,
        },
        |_| false,
        |_| (1.2, velocity),
    )
    .unwrap();
    OpenBoundaryGrid3::try_new(grid, cx, native_faces(axis, velocity)).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

fn uniform<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        for axis in 0..3 {
            for sign in [-1.0, 1.0] {
                let mut velocity = [0.011, -0.013, 0.009];
                velocity[axis] = 0.035 * sign;
                for force in [[0.0; 3], [0.0001, -0.00007, 0.00008]] {
                    let mut flow = fixture::<L, Q>(cx, collision, axis, velocity, force);
                    let cells = flow.populations().len();
                    let mut heat = flow
                        .prepare_thermal_field(
                            cx,
                            &vec![0.8; cells],
                            0.84,
                            thermal_faces(axis, 0.8),
                            limits(3.0),
                        )
                        .unwrap();
                    let result = flow
                        .step_serial_with_thermal(
                            &mut heat,
                            cx,
                            None,
                            &vec![0.0; cells],
                            None,
                            [0.0; 3],
                        )
                        .unwrap();
                    let area = cells / flow.dimensions()[axis];
                    near(
                        result.heat.reservoirs[2 * axis].heat_in,
                        area as f64 * velocity[axis] * 0.8,
                        2e-13,
                    );
                    near(
                        result.heat.reservoirs[2 * axis + 1].heat_in,
                        -(area as f64) * velocity[axis] * 0.8,
                        2e-13,
                    );
                    near(result.heat.defect, 0.0, 2e-13);
                    for cell in 0..cells {
                        near(heat.temperature(cell).unwrap(), 0.8, 2e-14);
                    }
                    assert_eq!(flow.completed_steps(), 1);
                    assert_eq!(heat.completed_steps(), 1);
                    assert_eq!(heat.last_balance(), Some(result.heat));
                }
            }
        }
    });
}

#[test]
fn d3q19_and_d3q27_transport_heat_from_native_force_corrected_velocity() {
    uniform::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    uniform::<D3q27, 27>(CollisionConfig27 {
        model: CollisionModel27::Bgk {
            tau: 0.9,
            equilibrium: Equilibrium27::Factorized,
        },
        maximum_mach: 0.2,
    });
}

#[test]
fn final_source_refusal_rolls_back_native_flow_and_thermal_history_then_retries_exactly() {
    with_cx(|cx| {
        let make = || {
            fixture::<D3q19, 19>(
                cx,
                CollisionModel3::Bgk { tau: 0.9 },
                0,
                [0.025, 0.0, 0.0],
                [0.0001, 0.0, 0.0],
            )
        };
        let mut actual = make();
        let mut expected = make();
        let cells = actual.populations().len();
        let mut heat = actual
            .prepare_thermal_field(
                cx,
                &vec![0.4; cells],
                0.82,
                thermal_faces(0, 0.4),
                limits(0.7),
            )
            .unwrap();
        let mut control = expected
            .prepare_thermal_field(
                cx,
                &vec![0.4; cells],
                0.82,
                thermal_faces(0, 0.4),
                limits(0.7),
            )
            .unwrap();
        let original_flow = actual.populations().to_vec();
        let original_heat = heat.populations().to_vec();
        let refusal =
            actual.step_serial_with_thermal(&mut heat, cx, None, &vec![0.4; cells], None, [0.0; 3]);
        assert!(matches!(
            refusal,
            Err(OpenThermalError3::Thermal(
                ThermalOpenError3::Temperature { .. }
            ))
        ));
        assert_eq!(actual.populations(), original_flow);
        assert_eq!(heat.populations(), original_heat);
        assert_eq!(actual.completed_steps(), 0);
        assert_eq!(heat.completed_steps(), 0);
        assert_eq!(heat.last_balance(), None);
        let source = vec![0.002; cells];
        let result = actual
            .step_serial_with_thermal(&mut heat, cx, None, &source, None, [0.0; 3])
            .unwrap();
        let reference = expected
            .step_serial_with_thermal(&mut control, cx, None, &source, None, [0.0; 3])
            .unwrap();
        assert_eq!(actual.populations(), expected.populations());
        assert_eq!(heat.populations(), control.populations());
        assert_eq!(result.momentum, reference.momentum);
        assert_eq!(result.heat, reference.heat);
    });
}

#[test]
fn momentum_refusal_and_skipped_native_time_cannot_publish_thermal_state() {
    with_cx(|cx| {
        let mut flow = fixture::<D3q19, 19>(
            cx,
            CollisionModel3::Bgk { tau: 0.9 },
            0,
            [0.02, 0.0, 0.0],
            [0.0; 3],
        );
        let cells = flow.populations().len();
        let mut thermal = flow
            .prepare_thermal_field(
                cx,
                &vec![0.4; cells],
                0.9,
                thermal_faces(0, 0.4),
                limits(2.0),
            )
            .unwrap();
        let before = thermal.populations().to_vec();
        let mut laws = vec![CollisionModel3::Bgk { tau: 0.9 }; cells];
        laws[cells - 1] = CollisionModel3::Bgk { tau: 0.0 };
        assert!(matches!(
            flow.step_serial_with_thermal(
                &mut thermal,
                cx,
                Some(&laws),
                &vec![0.01; cells],
                None,
                [0.0; 3]
            ),
            Err(OpenThermalError3::Flow(_))
        ));
        assert_eq!(thermal.populations(), before);
        assert_eq!(flow.completed_steps(), 0);
        flow.step_serial(cx, None, [0.0; 3]).unwrap();
        assert!(
            flow.step_serial_with_thermal(
                &mut thermal,
                cx,
                None,
                &vec![0.0; cells],
                None,
                [0.0; 3]
            )
            .is_err()
        );
        assert_eq!(thermal.populations(), before);
        assert_eq!(thermal.completed_steps(), 0);
    });
}

#[test]
fn every_native_open_face_requires_explicit_backflow_temperature() {
    with_cx(|cx| {
        let flow = fixture::<D3q19, 19>(
            cx,
            CollisionModel3::Bgk { tau: 0.9 },
            0,
            [0.02, 0.0, 0.0],
            [0.0; 3],
        );
        let cells = flow.populations().len();
        let mut boundaries = thermal_faces(0, 0.4);
        boundaries.faces[1] = ThermalOpenFace3::Wall(ThermalWall::Adiabatic);
        assert!(
            flow.prepare_thermal_field(cx, &vec![0.4; cells], 0.9, boundaries, limits(2.0))
                .is_err()
        );
        let mut budget = limits(2.0);
        budget.maximum_storage_bytes = 1;
        assert!(
            flow.prepare_thermal_field(cx, &vec![0.4; cells], 0.9, thermal_faces(0, 0.4), budget)
                .is_err()
        );
    });
}
