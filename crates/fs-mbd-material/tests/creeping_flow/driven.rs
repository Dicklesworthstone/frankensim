//! G1/G3: external drives and drag share one native velocity/work solve.
use super::*;
use fs_mbd::{Pose, RigidBodyState, UnitQuaternion, Vec3 as NativeVec3, Wrench};
use fs_mbd_material::creeping_flow::finite_reynolds::{FiniteReynoldsLimits, FiniteReynoldsSphere};

fn at_velocity(
    body: &fs_mbd_material::principal::PrincipalBody<'_>,
    velocity: NativeVec3,
) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(NativeVec3::ZERO, UnitQuaternion::IDENTITY).unwrap(),
        velocity.scale(body.properties().mass()),
        NativeVec3::ZERO,
    )
    .unwrap()
}

fn norm(v: NativeVec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}

#[test]
fn driven_stokes_translation_and_rotation_converge_to_independent_exponentials() {
    let density = 20000.0;
    let source = sphere(density);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 1.0);
    let ambient = NativeVec3::new(0.03, -0.02, 0.04);
    let flow = CreepingSphere::new(&material, RADIUS, Gravity::ZERO, ambient, limits()).unwrap();
    let terminal_relative = NativeVec3::new(0.2, -0.15, 0.1);
    let terminal_angular = NativeVec3::new(1.0, -0.5, 0.25);
    let c = 6.0 * core::f64::consts::PI * RADIUS;
    let c_rot = 8.0 * core::f64::consts::PI * RADIUS.powi(3);
    let applied = Wrench {
        force_world: terminal_relative.scale(c),
        torque_body: terminal_angular.scale(c_rot),
    };
    // Independent continuous-time solutions to m*v'=F-c*(v-U), I*w'=T-c_rot*w.
    let time = 0.8;
    let tau_translation = 2.0 * density * RADIUS.powi(2) / 9.0;
    let tau_rotation = density * RADIUS.powi(2) / 15.0;
    let exact_velocity = ambient.add(terminal_relative.scale(-(-time / tau_translation).exp_m1()));
    let exact_angular = terminal_angular.scale(-(-time / tau_rotation).exp_m1());
    let exact_position = ambient
        .scale(time)
        .add(terminal_relative.scale(time + tau_translation * (-time / tau_translation).exp_m1()));
    let mut errors = Vec::new();
    for steps in [40, 80, 160] {
        let dt = time / f64::from(steps);
        let mut state = at_velocity(&body, ambient);
        let mut applied_work = 0.0;
        let mut background_work = 0.0;
        let mut heat = 0.0;
        let initial_energy = RigidBodyIntegrator::new(Gravity::ZERO)
            .diagnostics(state, body.properties())
            .unwrap()
            .mechanical_energy;
        for _ in 0..steps {
            let step = flow
                .step_with_wrench(&body, state, applied, controls(dt), || false)
                .unwrap();
            let before_v = state
                .center_of_mass_velocity_world(body.properties())
                .unwrap();
            let after_v = step
                .dynamics
                .state_after
                .center_of_mass_velocity_world(body.properties())
                .unwrap();
            let before_w = body
                .properties()
                .angular_velocity_body_checked(state.angular_momentum_body())
                .unwrap();
            let after_w = body
                .properties()
                .angular_velocity_body_checked(step.dynamics.state_after.angular_momentum_body())
                .unwrap();
            close(
                step.applied_work_j,
                dt * (applied.force_world.dot(before_v.add(after_v).scale(0.5))
                    + applied.torque_body.dot(before_w.add(after_w).scale(0.5))),
            );
            assert!(step.viscous_heat_j >= 0.0);
            assert!(step.energy_defect_j.abs() < 1e-12);
            applied_work += step.applied_work_j;
            background_work += step.background_work_j;
            heat += step.viscous_heat_j;
            state = step.dynamics.state_after;
        }
        let final_energy = RigidBodyIntegrator::new(Gravity::ZERO)
            .diagnostics(state, body.properties())
            .unwrap()
            .mechanical_energy;
        assert!(
            (final_energy - initial_energy + heat - background_work - applied_work).abs() < 1e-12
        );
        errors.push([
            norm(
                state
                    .center_of_mass_velocity_world(body.properties())
                    .unwrap()
                    .sub(exact_velocity),
            ),
            norm(
                body.properties()
                    .angular_velocity_body_checked(state.angular_momentum_body())
                    .unwrap()
                    .sub(exact_angular),
            ),
            norm(state.pose().position_world().sub(exact_position)),
        ]);
    }
    for pair in errors.windows(2) {
        for axis in 0..3 {
            let ratio = pair[0][axis] / pair[1][axis];
            assert!(
                (3.9..4.1).contains(&ratio),
                "second-order error ratio {ratio}: {errors:?}"
            );
        }
    }
}

#[test]
fn driven_finite_reynolds_terminal_balance_separates_drive_flow_and_heat() {
    let source = sphere(20000.0);
    let body = source.principal_body(1e-12, || false).unwrap();
    let material = fluid(1.0, 0.001);
    let slip = NativeVec3::new(3.0, 4.0, 0.0);
    // An independently evaluated empirical coefficient supplies the drive;
    // it is not taken from the production force being tested.
    let resistance =
        6.0 * core::f64::consts::PI * 0.001 * RADIUS * (1.0 + 0.15 * 100.0_f64.powf(0.687));
    let applied = Wrench {
        force_world: slip.scale(resistance),
        torque_body: NativeVec3::ZERO,
    };
    let mut heats = Vec::new();
    for ambient in [NativeVec3::ZERO, NativeVec3::new(2.0, -1.0, 0.5)] {
        let flow = FiniteReynoldsSphere::new(
            &material,
            RADIUS,
            Gravity::ZERO,
            ambient,
            FiniteReynoldsLimits {
                maximum_reynolds: 1000.0,
                maximum_diffusion_ratio: 0.05,
            },
        )
        .unwrap();
        let velocity = ambient.add(slip);
        let initial = at_velocity(&body, velocity);
        let dt = 0.2;
        let step = flow
            .step_with_wrench(&body, initial, applied, controls(dt), || false)
            .unwrap();
        let final_velocity = step
            .dynamics
            .state_after
            .center_of_mass_velocity_world(body.properties())
            .unwrap();
        assert!(norm(final_velocity.sub(velocity)) < 1e-13);
        assert!(
            norm(
                step.dynamics
                    .state_after
                    .pose()
                    .position_world()
                    .sub(velocity.scale(dt))
            ) < 1e-13
        );
        close(step.load_midpoint.reynolds, 100.0);
        close(step.applied_work_j, dt * applied.force_world.dot(velocity));
        close(
            step.background_work_j,
            -dt * applied.force_world.dot(ambient),
        );
        close(step.viscous_heat_j, dt * resistance * 25.0);
        assert!((step.viscous_heat_j - step.applied_work_j - step.background_work_j).abs() < 1e-12);
        heats.push(step.viscous_heat_j);
        assert_eq!(
            flow.step(&body, initial, controls(dt), || false).unwrap(),
            flow.step_with_wrench(&body, initial, Wrench::ZERO, controls(dt), || false)
                .unwrap()
        );
        for invalid in [
            Wrench {
                force_world: NativeVec3::new(f64::NAN, 0.0, 0.0),
                ..Wrench::ZERO
            },
            Wrench {
                torque_body: NativeVec3::new(1e-12, 0.0, 0.0),
                ..Wrench::ZERO
            },
            Wrench {
                force_world: NativeVec3::new(1e8, 0.0, 0.0),
                ..Wrench::ZERO
            },
        ] {
            assert!(
                flow.step_with_wrench(&body, initial, invalid, controls(dt), || false)
                    .is_err()
            );
        }
        assert_eq!(
            flow.step_with_wrench(&body, initial, applied, controls(dt), || true),
            Err(CreepingFlowError::Cancelled)
        );
    }
    close(heats[0], heats[1]);
}
