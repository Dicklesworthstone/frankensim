//! G1/G2/G3: analytic gravity, Kepler refinement and conservative tidal loads.
use fs_mbd::gravity::{
    CentralGravity, GravityError, RigidGravityField, UniformGravity, step_in_field,
};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

const FRAME: u64 = 17;
fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e} (tolerance {tolerance:e})"
    );
}
fn body(mass: f64) -> MassProperties {
    MassProperties::new(mass, Vec3::ZERO, Vec3::new(0.2, 0.3, 0.4)).unwrap()
}
fn state(
    position: Vec3,
    orientation: UnitQuaternion,
    momentum: Vec3,
    spin: Vec3,
) -> RigidBodyState {
    RigidBodyState::new(Pose::new(position, orientation).unwrap(), momentum, spin).unwrap()
}

#[test]
fn point_field_matches_inverse_square_potential_and_extreme_radius_scaling() {
    let field = CentralGravity::new(FRAME, Vec3::new(1.0, -2.0, 3.0), 100.0).unwrap();
    let position = Vec3::new(3.0, -2.0, 3.0);
    let (a, u) = field.sample(FRAME, position).unwrap();
    assert_eq!(a, Vec3::new(-25.0, 0.0, 0.0));
    assert_eq!(u, -50.0);
    assert_eq!(
        field.sample(FRAME, Vec3::new(5.0, -2.0, 3.0)).unwrap().0.x,
        -6.25
    );
    let h = 1e-5;
    let up = field
        .sample(FRAME, position.add(Vec3::new(h, 0.0, 0.0)))
        .unwrap()
        .1;
    let um = field
        .sample(FRAME, position.sub(Vec3::new(h, 0.0, 0.0)))
        .unwrap()
        .1;
    close(-(up - um) / (2.0 * h), a.x, 2e-9);
    // Naive r*r or r*r*r intermediates would respectively overflow/underflow.
    for (mu, radius, expected) in [(1e300, 1e200, -1e-100), (1e-300, 1e-200, -1e100)] {
        let huge = CentralGravity::new(FRAME, Vec3::ZERO, mu).unwrap();
        close(
            huge.sample(FRAME, Vec3::new(radius, 0.0, 0.0)).unwrap().0.x / expected,
            1.0,
            2e-15,
        );
    }
}

#[test]
fn uniform_field_reproduces_native_ballistics_and_datum_shift() {
    let acceleration = Vec3::new(1.1, -0.7, -9.0);
    let gravity = Gravity::new(acceleration).unwrap();
    let datum = Vec3::new(2.0, -1.0, 0.5);
    let field = UniformGravity::new(FRAME, gravity, datum).unwrap();
    let properties = body(3.0);
    let initial = state(
        Vec3::new(2.0, 3.0, 4.0),
        UnitQuaternion::IDENTITY,
        Vec3::new(2.0, -1.0, 3.0),
        Vec3::new(0.1, -0.2, 0.3),
    );
    let expected = RigidBodyIntegrator::new(gravity)
        .step_symplectic(initial, properties, Wrench::ZERO, 0.03)
        .unwrap();
    let actual = step_in_field(
        &field,
        FRAME,
        initial,
        properties,
        Wrench::ZERO,
        0.03,
        || false,
    )
    .unwrap();
    assert_eq!(actual.state_after.pose(), expected.state_after.pose());
    close(
        norm(
            actual
                .state_after
                .linear_momentum_world()
                .sub(expected.state_after.linear_momentum_world()),
        ),
        0.0,
        1e-15,
    );
    assert_eq!(
        actual.state_after.angular_momentum_body(),
        expected.state_after.angular_momentum_body()
    );
    close(
        actual.diagnostics_before.gravitational_potential_energy
            - expected.diagnostics_before.gravitational_potential_energy,
        properties.mass() * acceleration.dot(datum),
        2e-14,
    );
}

fn kepler(steps: usize, periods: usize) -> (f64, f64, f64) {
    let field = CentralGravity::new(FRAME, Vec3::ZERO, 1.0).unwrap();
    let properties = body(2.0);
    let initial = state(
        Vec3::new(0.6, 0.0, 0.0),
        UnitQuaternion::IDENTITY,
        Vec3::new(0.0, 2.0 * (1.4_f64 / 0.6).sqrt(), 0.0),
        Vec3::ZERO,
    );
    let mut current = initial;
    let dt = core::f64::consts::TAU / steps as f64;
    let mut max_energy_error = 0.0_f64;
    let mut max_momentum_error = 0.0_f64;
    let angular = initial
        .pose()
        .position_world()
        .cross(initial.linear_momentum_world());
    for _ in 0..steps * periods {
        let step = step_in_field(&field, FRAME, current, properties, Wrench::ZERO, dt, || {
            false
        })
        .unwrap();
        max_energy_error =
            max_energy_error.max((step.diagnostics_after.mechanical_energy + 1.0).abs());
        max_momentum_error = max_momentum_error.max(norm(
            step.diagnostics_after.angular_momentum_world.sub(angular),
        ));
        current = step.state_after;
    }
    (
        norm(
            current
                .pose()
                .position_world()
                .sub(initial.pose().position_world()),
        ),
        max_energy_error,
        max_momentum_error,
    )
}

#[test]
fn eccentric_kepler_orbit_refines_at_second_order_without_secular_energy_growth() {
    let coarse = kepler(400, 1);
    let fine = kepler(800, 1);
    assert!(
        coarse.0 / fine.0 > 3.8 && coarse.0 / fine.0 < 4.2,
        "{coarse:?} {fine:?}"
    );
    assert!(coarse.1 / fine.1 > 3.8);
    let long = kepler(800, 40);
    assert!(
        long.1 < 1.01 * fine.1,
        "energy defect grew: {fine:?} {long:?}"
    );
    assert!(long.2 < 2e-12, "angular momentum drift: {long:?}");
}

#[test]
fn quadrupole_force_and_torque_are_gradients_of_the_same_potential() {
    let field = CentralGravity::new(FRAME, Vec3::ZERO, 100.0)
        .unwrap()
        .quadrupole(0.8, 0.15)
        .unwrap();
    let principal_frame = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.8).unwrap();
    let properties = MassProperties::from_principal_frame(
        2.0,
        Vec3::ZERO,
        Vec3::new(0.2, 0.3, 0.4),
        principal_frame,
    )
    .unwrap();
    let pose = Pose::new(
        Vec3::new(8.0, -3.0, 5.0),
        UnitQuaternion::from_axis_angle(Vec3::new(2.0, 1.0, -1.0), 0.4).unwrap(),
    )
    .unwrap();
    let load = field.load(FRAME, pose, properties).unwrap();
    for axis in [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ] {
        let h = 1e-4;
        let plus = Pose::new(pose.position_world().add(axis.scale(h)), pose.orientation()).unwrap();
        let minus =
            Pose::new(pose.position_world().sub(axis.scale(h)), pose.orientation()).unwrap();
        let derivative = (field
            .load(FRAME, plus, properties)
            .unwrap()
            .potential_energy_j
            - field
                .load(FRAME, minus, properties)
                .unwrap()
                .potential_energy_j)
            / (2.0 * h);
        close(-derivative, load.wrench.force_world.dot(axis), 5e-10);
        let h = 1e-5;
        let plus = Pose::new(
            pose.position_world(),
            pose.orientation().right_exp(axis.scale(h)).unwrap(),
        )
        .unwrap();
        let minus = Pose::new(
            pose.position_world(),
            pose.orientation().right_exp(axis.scale(-h)).unwrap(),
        )
        .unwrap();
        let derivative = (field
            .load(FRAME, plus, properties)
            .unwrap()
            .potential_energy_j
            - field
                .load(FRAME, minus, properties)
                .unwrap()
                .potential_energy_j)
            / (2.0 * h);
        close(-derivative, load.wrench.torque_body.dot(axis), 4e-10);
    }
    // The translational quadrupole correction balances the orientation torque.
    let net = pose.position_world().cross(load.wrench.force_world).add(
        pose.orientation()
            .rotate_body_to_world(load.wrench.torque_body),
    );
    assert!(norm(net) < 1e-14, "net central torque: {net:?}");
}

#[test]
fn tidal_model_matches_independent_two_mass_potential_with_fourth_order_remainder() {
    let source = CentralGravity::new(FRAME, Vec3::ZERO, 100.0).unwrap();
    let field = source.quadrupole(1.02, 0.2).unwrap();
    // Two unit-mass spherical bobs of radius .01 at body x = +/-1.
    let properties =
        MassProperties::new(2.0, Vec3::ZERO, Vec3::new(0.00008, 2.00008, 2.00008)).unwrap();
    let orientation = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), 0.4).unwrap();
    let arm = orientation.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
    let error = |radius| {
        let position = Vec3::new(radius, 0.0, 0.0);
        let exact = source.sample(FRAME, position.add(arm)).unwrap().1
            + source.sample(FRAME, position.sub(arm)).unwrap().1;
        let approximate = field
            .load(FRAME, Pose::new(position, orientation).unwrap(), properties)
            .unwrap()
            .potential_energy_j;
        (approximate / exact - 1.0).abs()
    };
    assert!(error(10.0) / error(20.0) > 15.0);
    assert!(error(20.0) < 3e-6);
}

#[test]
fn tidal_motion_preserves_total_angular_momentum() {
    let field = CentralGravity::new(FRAME, Vec3::ZERO, 100.0)
        .unwrap()
        .quadrupole(0.8, 0.15)
        .unwrap();
    let properties = body(2.0);
    let initial = state(
        Vec3::new(10.0, 0.0, 0.0),
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.7).unwrap(),
        Vec3::new(0.0, 2.0 * 10.0_f64.sqrt(), 0.0),
        Vec3::new(0.03, -0.02, 0.06),
    );
    let initial_l = RigidBodyIntegrator::new(Gravity::ZERO)
        .diagnostics(initial, properties)
        .unwrap()
        .angular_momentum_world;
    let mut current = initial;
    for _ in 0..2000 {
        let step = step_in_field(
            &field,
            FRAME,
            current,
            properties,
            Wrench::ZERO,
            0.01,
            || false,
        )
        .unwrap();
        assert!(norm(step.diagnostics_after.angular_momentum_world.sub(initial_l)) < 2e-11);
        current = step.state_after;
    }
}

#[test]
fn singular_frames_extent_and_late_cancellation_refuse_without_advancing_state() {
    assert!(CentralGravity::new(FRAME, Vec3::ZERO, 0.0).is_err());
    let source = CentralGravity::new(FRAME, Vec3::ZERO, 1.0).unwrap();
    assert!(source.sample(FRAME, Vec3::ZERO).is_err());
    assert_eq!(
        source.sample(FRAME + 1, Vec3::new(1.0, 0.0, 0.0)),
        Err(GravityError::FrameMismatch)
    );
    assert!(source.sample(FRAME, Vec3::new(f64::NAN, 0.0, 0.0)).is_err());
    let initial = state(
        Vec3::new(2.0, 0.0, 0.0),
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Vec3::ZERO,
    );
    let properties = body(2.0);
    let field = source.quadrupole(0.8, 0.1).unwrap();
    assert_eq!(
        field.load(FRAME, initial.pose(), properties),
        Err(GravityError::OutsideExpansionDomain)
    );
    let invalid_extent = source.quadrupole(0.001, 0.1).unwrap();
    assert!(
        invalid_extent
            .load(FRAME, initial.pose(), properties)
            .is_err()
    );
    for stop in 1..=4 {
        let mut polls = 0;
        let saved = initial;
        assert_eq!(
            step_in_field(
                &source,
                FRAME,
                initial,
                properties,
                Wrench::ZERO,
                0.01,
                || {
                    polls += 1;
                    polls == stop
                }
            ),
            Err(GravityError::Cancelled)
        );
        assert_eq!(initial, saved);
    }
    assert!(
        step_in_field(
            &source,
            FRAME,
            initial,
            properties,
            Wrench::ZERO,
            0.01,
            || false
        )
        .is_ok()
    );
}
