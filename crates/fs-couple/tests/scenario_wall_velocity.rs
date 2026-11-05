//! G0/G1/G3/G4: physical scenario drives use the original moving-wall step.

use fs_couple::scenario_wall_velocity::{
    FrameWallStepError, FrameWallVelocity3, FrameWallVelocityError,
};
use fs_exec::{CancelGate, Cx};
use fs_ga::Vec3;
use fs_lbm::d3q19::CollisionModel3;
use fs_lbm::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Equilibrium27};
use fs_lbm::lattice3::{D3q19, EquilibriumLattice3, Grid3, GridBudget3};
use fs_lbm::material::plan_scaling_from_material;
use fs_material::fluid_transport::FluidTransportState;
use fs_qty::{Dims, QtyAny};
use fs_scenario::frame::FrameKinematicsBudget;
use fs_scenario::{Frame, FrameId, FrameMotion, FrameTree, TimeSignal, WORLD};
use std::ops::ControlFlow;

const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 7391,
                kernel_id: 7391,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}

fn frames(motion: FrameMotion) -> FrameTree {
    let mut result = FrameTree::new();
    result.add(Frame {
        id: FrameId(1),
        parent: WORLD,
        name: "wall-drive".into(),
        motion,
    });
    result
}

fn translation(direction: [f64; 3], distance: f64) -> FrameMotion {
    FrameMotion::Translating {
        direction,
        displacement: TimeSignal::Ramp {
            t_start: 0.0,
            t_end: 1.0,
            from: QtyAny::new(0.0, LENGTH),
            to: QtyAny::new(distance, LENGTH),
        },
    }
}

fn budget() -> FrameKinematicsBudget {
    FrameKinematicsBudget {
        max_frames: 8,
        max_history_scalars: 64,
    }
}

fn fluid() -> FluidTransportState {
    FluidTransportState::try_new(1.0, 1000.0, 0.03, 0.02).unwrap()
}

fn grid<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) -> Grid3<L, Q> {
    Grid3::try_new_equilibrium(
        [4, 3, 5],
        collision,
        [0.0; 3],
        [true, true, false],
        GridBudget3 {
            maximum_cells: 60,
            maximum_mach: 0.2,
        },
        |_| false,
        |_| (1.0, [0.0; 3]),
    )
    .unwrap()
}

fn q27(tau: f64) -> CollisionConfig27 {
    CollisionConfig27 {
        model: CollisionModel27::Bgk {
            tau,
            equilibrium: Equilibrium27::Factorized,
        },
        maximum_mach: 0.2,
    }
}

fn uniform_translation_matches_native<L: EquilibriumLattice3<Q>, const Q: usize>(
    collision: L::Collision,
) {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let frames = frames(translation([1.0, 0.0, 0.0], 0.2));
    let prepared = frames
        .prepare_kinematics(FrameId(1), budget(), &mut || ControlFlow::Continue(()))
        .unwrap();
    let drive = FrameWallVelocity3::try_new(&prepared, &scaling, [1.0, -2.0, 0.1]).unwrap();
    let mut actual = grid::<L, Q>(collision);
    let mut oracle = grid::<L, Q>(collision);
    with_cx(|cx| {
        let probe = actual
            .wall_momentum_probe(cx, [2.0, 1.5, 0.0], |link| link.exterior_faces[4])
            .unwrap();
        let expected_velocity = [scaling.velocity_lattice(0.2).unwrap(), 0.0, 0.0];
        let velocities = vec![expected_velocity; probe.measured_links()];
        let mut scratch = vec![[0.0; 3]; probe.measured_links()];
        for time in [0.25, 0.75] {
            let expected = oracle
                .step_serial_with_moving_wall_momentum(cx, &probe, &velocities, None)
                .unwrap();
            let observed = drive
                .step_serial(&mut actual, cx, &probe, &mut scratch, time, None)
                .unwrap();
            assert_eq!(scratch, velocities);
            assert_eq!(observed, expected);
            assert_eq!(actual.populations(), oracle.populations());
            assert_eq!(actual.completed_steps(), oracle.completed_steps());
        }
    });
}

#[test]
fn translating_drive_matches_native_momentum_and_work_on_both_lattices() {
    let material = fluid();
    let tau = plan_scaling_from_material(&material, 0.02, 0.001, 1.0)
        .unwrap()
        .flow()
        .tau;
    uniform_translation_matches_native::<D3q19, 19>(CollisionModel3::Bgk { tau });
    uniform_translation_matches_native::<D3q27, 27>(q27(tau));
}

#[test]
fn rotating_drive_samples_fixed_world_links_instead_of_advected_material_points() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let origin = [1.0, -2.0, 0.1];
    let center = Vec3::new(1.06, -1.94, 0.1);
    let omega = 1.25;
    let frames = frames(FrameMotion::Rotating {
        axis: [0.0, 0.0, 1.0],
        center,
        rate: QtyAny::new(omega, Dims([0, 0, -1, 0, 0, 0])),
    });
    let prepared = frames
        .prepare_kinematics(FrameId(1), budget(), &mut || ControlFlow::Continue(()))
        .unwrap();
    let drive = FrameWallVelocity3::try_new(&prepared, &scaling, origin).unwrap();
    let flow = grid::<D3q27, 27>(q27(scaling.flow().tau));
    with_cx(|cx| {
        let probe = flow
            .wall_momentum_probe(cx, [0.0; 3], |link| link.exterior_faces[4])
            .unwrap();
        let mut scratch = vec![[0.0; 3]; probe.measured_links()];
        for time in [0.25, 0.75] {
            drive
                .sample_velocities(
                    time,
                    &probe,
                    &mut scratch,
                    &mut || ControlFlow::Continue(()),
                )
                .unwrap();
            for (link, observed) in probe.wall_links().zip(&scratch) {
                let x = origin[0] + 0.02 * link.position[0];
                let y = origin[1] + 0.02 * link.position[1];
                let expected = [
                    -omega * (y - center.y) / 20.0,
                    omega * (x - center.x) / 20.0,
                    0.0,
                ];
                for axis in 0..3 {
                    assert!((observed[axis] - expected[axis]).abs() < 2e-14);
                }
            }
        }
        assert_eq!(flow.completed_steps(), 0);
    });
}

#[test]
fn forbidden_normal_motion_and_excessive_speed_do_not_publish_fluid_steps() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let mut flow = grid::<D3q27, 27>(q27(scaling.flow().tau));
    let before = flow.populations().to_vec();
    with_cx(|cx| {
        let probe = flow
            .wall_momentum_probe(cx, [0.0; 3], |link| link.exterior_faces[4])
            .unwrap();
        let mut scratch = vec![[0.0; 3]; probe.measured_links()];
        for (direction, speed, native_refusal) in
            [([0.0, 0.0, 1.0], 0.2, true), ([1.0, 0.0, 0.0], 2.0, false)]
        {
            let frames = frames(translation(direction, speed));
            let prepared = frames
                .prepare_kinematics(FrameId(1), budget(), &mut || ControlFlow::Continue(()))
                .unwrap();
            let drive = FrameWallVelocity3::try_new(&prepared, &scaling, [0.0; 3]).unwrap();
            let result = drive.step_serial(&mut flow, cx, &probe, &mut scratch, 0.5, None);
            if native_refusal {
                assert!(matches!(result, Err(FrameWallStepError::Native(_))));
            } else {
                assert!(matches!(
                    result,
                    Err(FrameWallStepError::Velocity(
                        FrameWallVelocityError::SpeedBound { .. }
                    ))
                ));
            }
            assert_eq!(flow.populations(), before);
            assert_eq!(flow.completed_steps(), 0);
        }
    });
}

#[test]
fn cancelled_sampling_and_wrong_scratch_length_never_touch_the_grid() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let frames = frames(translation([1.0, 0.0, 0.0], 0.2));
    let prepared = frames
        .prepare_kinematics(FrameId(1), budget(), &mut || ControlFlow::Continue(()))
        .unwrap();
    let drive = FrameWallVelocity3::try_new(&prepared, &scaling, [0.0; 3]).unwrap();
    let flow = grid::<D3q19, 19>(CollisionModel3::Bgk {
        tau: scaling.flow().tau,
    });
    let before = flow.populations().to_vec();
    with_cx(|cx| {
        let probe = flow
            .wall_momentum_probe(cx, [0.0; 3], |link| link.exterior_faces[4])
            .unwrap();
        let mut scratch = vec![[0.0; 3]; probe.measured_links()];
        assert!(matches!(
            drive.sample_velocities(0.5, &probe, &mut scratch, &mut || ControlFlow::Break(())),
            Err(FrameWallVelocityError::Cancelled)
        ));
        assert!(
            drive
                .sample_velocities(0.5, &probe, &mut [], &mut || ControlFlow::Continue(()))
                .is_err()
        );
        assert_eq!(flow.populations(), before);
        assert_eq!(flow.completed_steps(), 0);
    });
}
