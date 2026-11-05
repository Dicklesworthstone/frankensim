//! The same typed histories reach the original D3Q27 source and reservoir step.

use super::{WORLD, boundary, fluid, near};
use fs_couple::scenario_velocity::ScenarioVelocityError;
use fs_couple::scenario_velocity::lattice3::{
    Error, rebind_velocity_face3, rebind_velocity_faces3,
};
use fs_exec::{CancelGate, Cx};
use fs_ga::{Quat, Vec3};
use fs_lbm::d3q19::{BoundarySpec3, Face3, FaceBoundary3};
use fs_lbm::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Equilibrium27};
use fs_lbm::lattice3::open_boundaries::OpenBoundaryGrid3;
use fs_lbm::lattice3::{Grid3, GridBudget3, GridError3};
use fs_lbm::material::plan_scaling_from_material;
use fs_scenario::vector::VectorEvaluationError;
use fs_scenario::{Frame, FrameId, FrameMotion, FrameTree};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            fs_exec::StreamKey {
                seed: 0xCE10,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}

fn grid(
    cx: &Cx<'_>,
    tau: f64,
    maximum_mach: f64,
    paired_velocity: bool,
) -> OpenBoundaryGrid3<D3q27, 27> {
    let lattice = Grid3::<D3q27, 27>::try_new_equilibrium(
        [4, 4, 5],
        CollisionConfig27 {
            model: CollisionModel27::Bgk {
                tau,
                equilibrium: Equilibrium27::Factorized,
            },
            maximum_mach,
        },
        [2e-5, -1e-5, 3e-5],
        [true, true, false],
        GridBudget3 {
            maximum_cells: 80,
            maximum_mach: 0.2,
        },
        |_| false,
        |_| (1.0, [0.0; 3]),
    )
    .unwrap();
    OpenBoundaryGrid3::try_new(
        lattice,
        cx,
        BoundarySpec3::new([
            FaceBoundary3::Periodic,
            FaceBoundary3::Periodic,
            FaceBoundary3::Periodic,
            FaceBoundary3::Periodic,
            FaceBoundary3::Velocity { velocity: [0.0; 3] },
            if paired_velocity {
                FaceBoundary3::Velocity { velocity: [0.0; 3] }
            } else {
                FaceBoundary3::Pressure { density: 1.0 }
            },
        ]),
    )
    .unwrap()
}

#[test]
fn signed_rotated_history_reaches_force_corrected_d3q27_inlet() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let mut frames = FrameTree::new();
    frames.add(Frame {
        id: FrameId(1),
        name: "inlet_axes".into(),
        parent: WORLD,
        motion: FrameMotion::Fixed {
            orientation: Quat::from_axis_angle(
                Vec3::new(0.0, 0.0, 1.0),
                std::f64::consts::FRAC_PI_2,
            ),
            translation: Vec3::new(1e6, -2e6, 3e6),
        },
    });
    let inlet = boundary([0.4, -0.2, 0.6], [-0.2, 0.4, -0.6], FrameId(1));
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut flow = grid(cx, scaling.flow().tau, 0.2, false);
        let force = flow.force_density();
        for (step, (time_s, expected_si)) in [
            (0.0, [0.2, 0.4, 0.6]),
            (0.5, [-0.1, 0.1, 0.0]),
            (1.0, [-0.4, -0.2, -0.6]),
        ]
        .into_iter()
        .enumerate()
        {
            let before = flow.populations().to_vec();
            let sample = rebind_velocity_face3(
                &mut flow,
                Face3::ZMin,
                &inlet,
                &frames,
                time_s,
                &scaling,
                &gate,
            )
            .unwrap();
            near(sample.world_m_s, expected_si, 1e-14);
            near(sample.lattice, expected_si.map(|u| u / 20.0), 1e-15);
            assert_eq!(flow.populations(), before);
            assert_eq!(flow.completed_steps(), step);
            assert_eq!(flow.force_density(), force);
            flow.step_serial(cx, None, [0.0; 3]).unwrap();
            assert_eq!(flow.completed_steps(), step + 1);
            for y in 0..4 {
                for x in 0..4 {
                    near(
                        flow.macroscopic(y * 4 + x).unwrap().1,
                        sample.lattice,
                        5e-14,
                    );
                    let outlet = flow.macroscopic((4 * 4 + y) * 4 + x).unwrap();
                    assert!((outlet.0 - 1.0).abs() < 5e-14);
                }
            }
            assert_eq!(
                flow.boundaries().face(Face3::ZMax),
                FaceBoundary3::Pressure { density: 1.0 }
            );
        }
    });
}

#[test]
fn later_face_or_native_collision_bound_refuses_without_publication() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let frames = FrameTree::new();
    let gate = CancelGate::new_clock_free();
    let inlet = boundary([0.0, 0.0, 0.5], [0.0, 0.0, 0.5], WORLD);
    let too_fast = boundary([0.7; 3], [0.7; 3], WORLD);
    with_cx(&gate, |cx| {
        let mut flow = grid(cx, scaling.flow().tau, 0.2, true);
        let before = flow.boundaries();
        let populations = flow.populations().to_vec();
        assert!(matches!(
            rebind_velocity_faces3(
                &mut flow,
                [None, None, None, None, Some(&inlet), Some(&too_fast)],
                &frames,
                0.5,
                &scaling,
                &gate,
            ),
            Err(Error::Scenario(ScenarioVelocityError::SpeedBound { .. }))
        ));
        assert_eq!(flow.boundaries(), before);
        assert_eq!(flow.populations(), populations);
        assert_eq!(flow.completed_steps(), 0);

        let mut restricted = grid(cx, scaling.flow().tau, 0.04, false);
        let before = restricted.boundaries();
        let populations = restricted.populations().to_vec();
        let admitted_by_scaling = boundary([0.0, 0.0, 0.8], [0.0, 0.0, 0.8], WORLD);
        assert!(matches!(
            rebind_velocity_face3(
                &mut restricted,
                Face3::ZMin,
                &admitted_by_scaling,
                &frames,
                0.5,
                &scaling,
                &gate,
            ),
            Err(Error::Boundary(GridError3::Mach { maximum, .. })) if maximum == 0.04
        ));
        assert_eq!(restricted.boundaries(), before);
        assert_eq!(restricted.populations(), populations);
        assert_eq!(restricted.completed_steps(), 0);
    });
}

#[test]
fn generic_native_kind_error_and_cancellation_keep_all_controls_and_state() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let frames = FrameTree::new();
    let gate = CancelGate::new_clock_free();
    let inlet = boundary([0.0, 0.0, 0.5], [0.0, 0.0, 0.5], WORLD);
    let mut flow = with_cx(&gate, |cx| grid(cx, scaling.flow().tau, 0.2, false));
    let before = flow.boundaries();
    let populations = flow.populations().to_vec();
    let force = flow.force_density();
    assert!(matches!(
        rebind_velocity_faces3(
            &mut flow,
            [None, None, None, None, Some(&inlet), Some(&inlet)],
            &frames,
            0.5,
            &scaling,
            &gate,
        ),
        Err(Error::Boundary(GridError3::Input(_)))
    ));
    gate.request();
    assert!(matches!(
        rebind_velocity_face3(
            &mut flow,
            Face3::ZMin,
            &inlet,
            &frames,
            0.5,
            &scaling,
            &gate,
        ),
        Err(Error::Scenario(ScenarioVelocityError::Scenario(
            VectorEvaluationError::Cancelled
        )))
    ));
    assert_eq!(flow.boundaries(), before);
    assert_eq!(flow.populations(), populations);
    assert_eq!(flow.completed_steps(), 0);
    assert_eq!(flow.force_density(), force);
}
