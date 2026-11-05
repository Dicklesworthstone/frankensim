//! Public rigid-body stages: retained angular impulses must also drive attitude.

use fs_mbd::spatial_loads::PreparedSpatialLoads;
use fs_mbd::{
    DynamicsError, Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState,
    UnitQuaternion, Vec3, Wrench,
};

const LARGE: f64 = 9_007_199_254_740_992.0; // 2^53; a unit half-kick rounds away.
const AXES: [Vec3; 3] = [
    Vec3::new(1.0, 0.0, 0.0),
    Vec3::new(0.0, 1.0, 0.0),
    Vec3::new(0.0, 0.0, 1.0),
];

fn engine() -> RigidBodyIntegrator {
    RigidBodyIntegrator::new(Gravity::ZERO)
}

fn torque(value: Vec3) -> Wrench {
    Wrench {
        force_world: Vec3::ZERO,
        torque_body: value,
    }
}

fn spinning(orientation: UnitQuaternion, angular_momentum: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(Vec3::ZERO, orientation).unwrap(),
        Vec3::ZERO,
        angular_momentum,
    )
    .unwrap()
}

fn reversed(state: RigidBodyState) -> RigidBodyState {
    RigidBodyState::new(
        state.pose(),
        state.linear_momentum_world().scale(-1.0),
        state.angular_momentum_body().scale(-1.0),
    )
    .unwrap()
}

fn world_spin(state: RigidBodyState) -> Vec3 {
    state
        .pose()
        .orientation()
        .rotate_body_to_world(state.angular_momentum_body())
}

fn vector_close(actual: Vec3, expected: Vec3, tolerance: f64) {
    let error = actual.sub(expected);
    assert!(
        error.x.abs().max(error.y.abs()).max(error.z.abs()) <= tolerance,
        "{actual:?} != {expected:?}; tolerance={tolerance:e}"
    );
}

fn rotation_close(actual: UnitQuaternion, expected: UnitQuaternion, tolerance: f64) {
    for axis in AXES {
        vector_close(
            actual.rotate_body_to_world(axis),
            expected.rotate_body_to_world(axis),
            tolerance,
        );
    }
}

fn kicked_axial_orientation(axis: Vec3, sign: f64) -> UnitQuaternion {
    // Independent scalar-trig oracle for theta = sign*(2^53 + 1).
    // Do not round that sum before taking trig: its half-angle is the
    // unevaluated pair (sign*2^52, sign/2). No compensated time-owner helper
    // enters this oracle. Losing the unit kick changes attitude by one radian.
    let (sin_high, cos_high) = (sign * (0.5 * LARGE)).sin_cos();
    let (sin_low, cos_low) = (0.5 * sign).sin_cos();
    let sine = sin_high.mul_add(cos_low, cos_high * sin_low);
    let cosine = cos_high.mul_add(cos_low, -sin_high * sin_low);
    UnitQuaternion::new(cosine, axis.x * sine, axis.y * sine, axis.z * sine).unwrap()
}

#[test]
fn spherical_stages_retain_the_kick_in_attitude_and_signed_endpoint_impulses() {
    let properties = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    let axis = AXES[2];
    for sign in [-1.0, 1.0] {
        let initial = spinning(UnitQuaternion::IDENTITY, axis.scale(sign * LARGE));
        let first = torque(axis.scale(sign * 2.0));
        let expected_rotation = kicked_axial_orientation(axis, sign);
        for endpoint in [-6.0, -2.0, 2.0, 6.0] {
            let stage = engine()
                .begin_verlet(initial, properties, first, 1.0)
                .unwrap();
            // The callback keeps its original f64 representation, while the
            // retained unit midpoint tail must already have driven the pose.
            assert_eq!(
                stage.drifted_state().angular_momentum_body(),
                initial.angular_momentum_body()
            );
            rotation_close(
                stage.drifted_state().pose().orientation(),
                expected_rotation,
                4e-14,
            );
            let after = stage
                .finish(torque(axis.scale(sign * endpoint)))
                .unwrap()
                .state_after;
            let impulse = axis.scale(sign * (1.0 + 0.5 * endpoint));
            assert_eq!(
                after.angular_momentum_body(),
                initial.angular_momentum_body().add(impulse)
            );
            rotation_close(after.pose().orientation(), expected_rotation, 4e-14);
            assert_eq!(after.linear_momentum_world(), Vec3::ZERO);
        }

        // Both held-load public paths must use the same complete angular
        // kick/drift treatment. Constant axial torque has an exact phase here.
        let expected_momentum = axis.scale(sign * (LARGE + 2.0));
        let symplectic = engine()
            .step_symplectic(initial, properties, first, 1.0)
            .unwrap()
            .state_after;
        let loads = PreparedSpatialLoads::new(first, Vec3::ZERO, &[], || false).unwrap();
        let spatial = engine()
            .step_spatial(initial, properties, &loads, 1.0)
            .unwrap()
            .step
            .state_after;
        for after in [symplectic, spatial] {
            assert_eq!(after.angular_momentum_body(), expected_momentum);
            rotation_close(after.pose().orientation(), expected_rotation, 4e-14);
        }
    }
}

#[test]
fn rotated_full_inertia_retains_the_exact_principal_axis_phase() {
    // This tensor has a nonzero yz entry, while x remains an exact principal
    // axis. It exercises the cached basis and increment conjugation without
    // rounding an authored large spin into a different physical direction.
    let basis = UnitQuaternion::from_axis_angle(AXES[0], 0.37).unwrap();
    let properties =
        MassProperties::from_principal_frame(1.0, Vec3::ZERO, Vec3::new(1.0, 1.5, 2.0), basis)
            .unwrap();
    assert!(properties.inertia_tensor_body().m[5].abs() > 0.05);
    for sign in [-1.0, 1.0] {
        let initial = spinning(UnitQuaternion::IDENTITY, AXES[0].scale(sign * LARGE));
        let held = torque(AXES[0].scale(sign * 2.0));
        let (step, ()) = engine()
            .step_verlet_with(initial, properties, held, 1.0, |_| {
                Ok::<_, DynamicsError>((held, ()))
            })
            .unwrap();
        assert_eq!(
            step.state_after.angular_momentum_body(),
            AXES[0].scale(sign * (LARGE + 2.0))
        );
        rotation_close(
            step.state_after.pose().orientation(),
            kicked_axial_orientation(AXES[0], sign),
            6e-14,
        );
    }
}

#[test]
fn full_kick_retains_the_inertia_quotient_tail_in_its_rotation_phase() {
    let properties = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(3.0, 3.0, 3.0)).unwrap();
    let initial = spinning(UnitQuaternion::IDENTITY, AXES[2].scale(LARGE));
    let after = engine()
        .step_kick_drift(initial, properties, torque(AXES[2].scale(2.0)), 1.0)
        .unwrap()
        .state_after;

    // The full kicked momentum is representable, but its exact quaternion
    // half-angle (2^53 + 2)/6 is not. Rounding angular velocity first changes
    // this half-angle by 1/12 radian. Build the independent scalar-trig oracle
    // from the unevaluated pair 1501199875790165.5 + 1/6 instead.
    let (sin_high, cos_high) = 1_501_199_875_790_165.5_f64.sin_cos();
    let (sin_low, cos_low) = (1.0_f64 / 6.0).sin_cos();
    let sine = sin_high.mul_add(cos_low, cos_high * sin_low);
    let cosine = cos_high.mul_add(cos_low, -sin_high * sin_low);
    let expected = UnitQuaternion::new(cosine, 0.0, 0.0, sine).unwrap();

    assert_eq!(after.angular_momentum_body(), AXES[2].scale(LARGE + 2.0));
    assert_eq!(after.linear_momentum_world(), Vec3::ZERO);
    assert_eq!(after.pose().position_world(), Vec3::ZERO);
    rotation_close(after.pose().orientation(), expected, 4e-14);
}

#[test]
fn held_symplectic_retains_the_full_impulse_at_the_smallest_duration() {
    let properties = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    let initial = spinning(UnitQuaternion::IDENTITY, Vec3::ZERO);
    let duration = f64::from_bits(1);
    assert_eq!(0.5 * duration, 0.0);

    for magnitude in [1.0, 1.0e308] {
        let force = Vec3::new(magnitude, 0.0, 0.0);
        let expected_impulse = force.scale(duration);
        assert!(expected_impulse.x > 0.0);

        // The complete impulses are one minimum subnormal and about 4.94e-16.
        // Half the duration underflows in both cases; with unit force, even
        // the half-impulse is unrepresentable. Preserve the full held impulse.
        let after = engine()
            .step_symplectic(
                initial,
                properties,
                Wrench {
                    force_world: force,
                    torque_body: Vec3::ZERO,
                },
                duration,
            )
            .unwrap()
            .state_after;
        assert_eq!(after.linear_momentum_world(), expected_impulse);
        assert_eq!(after.angular_momentum_body(), Vec3::ZERO);
        assert_eq!(after.pose().position_world(), Vec3::ZERO);
        assert_eq!(after.pose().orientation(), UnitQuaternion::IDENTITY);
    }
}

#[test]
fn public_rigid_steps_retain_the_linear_kick_through_position_cancellation() {
    let properties = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    for sign in [-1.0, 1.0] {
        let initial = RigidBodyState::new(
            Pose::new(Vec3::new(-sign * LARGE, 0.0, 0.0), UnitQuaternion::IDENTITY).unwrap(),
            Vec3::new(sign * LARGE, 0.0, 0.0),
            Vec3::ZERO,
        )
        .unwrap();
        let held = Wrench {
            force_world: Vec3::new(sign * 2.0, 0.0, 0.0),
            torque_body: Vec3::ZERO,
        };
        let expected_momentum = Vec3::new(sign * (LARGE + 2.0), 0.0, 0.0);
        for step in [
            engine().step(initial, properties, held, 1.0).unwrap(),
            engine()
                .step_cayley(initial, properties, held, 1.0)
                .unwrap(),
            engine()
                .step_symplectic(initial, properties, held, 1.0)
                .unwrap(),
        ] {
            // x0 cancels the large free displacement exactly; the midpoint
            // impulse must survive that contraction as a unit displacement.
            assert_eq!(
                step.state_after.pose().position_world(),
                Vec3::new(sign, 0.0, 0.0)
            );
            assert_eq!(step.state_after.linear_momentum_world(), expected_momentum);
        }
        let kicked = engine()
            .step_kick_drift(initial, properties, held, 1.0)
            .unwrap()
            .state_after;
        assert_eq!(
            kicked.pose().position_world(),
            Vec3::new(sign * 2.0, 0.0, 0.0)
        );
        assert_eq!(kicked.linear_momentum_world(), expected_momentum);
    }
}

fn asymmetric_body() -> (MassProperties, RigidBodyState) {
    let basis = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 0.7), 0.61).unwrap();
    let properties =
        MassProperties::from_principal_frame(2.3, Vec3::ZERO, Vec3::new(1.2, 1.8, 2.4), basis)
            .unwrap();
    let orientation = UnitQuaternion::from_axis_angle(Vec3::new(-0.4, 0.3, 1.0), 0.43).unwrap();
    (properties, spinning(orientation, Vec3::new(0.7, -0.5, 1.3)))
}

#[test]
fn differing_endpoint_torques_close_world_spin_and_reverse_in_their_own_frames() {
    let (properties, initial) = asymmetric_body();
    let first = Vec3::new(0.31, -0.27, 0.14);
    let last = Vec3::new(-0.21, 0.16, 0.35);
    let dt = 0.08;
    let (step, ()) = engine()
        .step_verlet_with(initial, properties, torque(first), dt, |_| {
            Ok::<_, DynamicsError>((torque(last), ()))
        })
        .unwrap();
    let after = step.state_after;
    let expected_impulse = initial
        .pose()
        .orientation()
        .rotate_body_to_world(first)
        .add(after.pose().orientation().rotate_body_to_world(last))
        .scale(0.5 * dt);
    vector_close(
        world_spin(after).sub(world_spin(initial)),
        expected_impulse,
        2e-12,
    );

    // Reversing momentum and exchanging the endpoint loads reverses the
    // symmetric discrete map; no exact forced-energy claim is needed.
    let (back, ()) = engine()
        .step_verlet_with(reversed(after), properties, torque(last), dt, |_| {
            Ok::<_, DynamicsError>((torque(first), ()))
        })
        .unwrap();
    rotation_close(
        back.state_after.pose().orientation(),
        initial.pose().orientation(),
        2e-12,
    );
    vector_close(
        back.state_after.angular_momentum_body(),
        initial.angular_momentum_body().scale(-1.0),
        2e-12,
    );
}

#[test]
fn full_inertia_spatial_couples_close_their_impulse_and_retrace_the_same_field() {
    let (properties, initial) = asymmetric_body();
    let body_couple = Vec3::new(0.17, -0.23, 0.19);
    let world_couple = Vec3::new(0.2, 0.1, -0.3);
    let loads =
        PreparedSpatialLoads::new(torque(body_couple), world_couple, &[], || false).unwrap();
    let dt = 0.08;
    let forward = engine()
        .step_spatial(initial, properties, &loads, dt)
        .unwrap();
    let after = forward.step.state_after;
    let expected_impulse = world_couple.scale(dt).add(
        initial
            .pose()
            .orientation()
            .rotate_body_to_world(body_couple)
            .add(after.pose().orientation().rotate_body_to_world(body_couple))
            .scale(0.5 * dt),
    );
    vector_close(forward.angular_impulse_world, expected_impulse, 2e-12);
    vector_close(
        world_spin(after).sub(world_spin(initial)),
        expected_impulse,
        2e-12,
    );
    vector_close(forward.angular_momentum_residual_world, Vec3::ZERO, 2e-12);

    // The world couple must be pulled back again at each reversed endpoint;
    // freezing the first body torque does not retrace this physical field.
    let back = engine()
        .step_spatial(reversed(after), properties, &loads, dt)
        .unwrap()
        .step
        .state_after;
    rotation_close(
        back.pose().orientation(),
        initial.pose().orientation(),
        2e-12,
    );
    vector_close(
        back.angular_momentum_body(),
        initial.angular_momentum_body().scale(-1.0),
        2e-12,
    );
}
