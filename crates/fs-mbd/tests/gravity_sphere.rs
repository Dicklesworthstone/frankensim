//! G1/G2/G3: exact finite-source gravity, conservative work and native motion.
use fs_mbd::gravity::{
    CentralGravity, GravityAdvanceOptions, GravityError, PointGravityField, RigidGravityField,
    UniformSphereGravity, advance_in_field,
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3, Wrench};

const FRAME: u64 = 113;

fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn vector_close(actual: Vec3, expected: Vec3, tolerance: f64) {
    assert!(
        norm(actual.sub(expected)) <= tolerance,
        "actual={actual:?}, expected={expected:?}, tolerance={tolerance:e}"
    );
}

fn axis(index: usize, value: f64) -> Vec3 {
    let mut xyz = [0.0; 3];
    xyz[index] = value;
    Vec3::new(xyz[0], xyz[1], xyz[2])
}

fn values(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

fn mass(value: f64) -> MassProperties {
    MassProperties::new(value, Vec3::ZERO, Vec3::new(0.1, 0.1, 0.1)).unwrap()
}

#[test]
fn finite_center_and_interior_retain_poisson_trace_and_original_rigid_mass() {
    let center = Vec3::new(0.5, -0.25, 1.0);
    // mu/R^3 = 2 s^-2; phi(0) = -12 J/kg.
    let field = UniformSphereGravity::new(FRAME, center, 16.0, 2.0).unwrap();
    for displacement in [
        Vec3::ZERO,
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.5, -0.25, 0.75),
    ] {
        let position = center.add(displacement);
        let sample = field.sample_point(FRAME, position).unwrap();
        let radius_squared = displacement.dot(displacement);
        vector_close(
            sample.acceleration_world_m_s2,
            displacement.scale(-2.0),
            2e-15,
        );
        close(
            sample.specific_potential_j_kg,
            -12.0 + radius_squared,
            3e-15,
        );
        assert_eq!(
            sample.acceleration_gradient_s2.m,
            [-2.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, -2.0]
        );
        let d = sample.acceleration_gradient_s2.m;
        close(d[0] + d[4] + d[8], -6.0, 1e-15);

        let properties = mass(3.5);
        let pose = Pose::new(
            position,
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 3.0), 0.7).unwrap(),
        )
        .unwrap();
        let load = field.load(FRAME, pose, properties).unwrap();
        vector_close(
            load.wrench.force_world,
            sample.acceleration_world_m_s2.scale(3.5),
            2e-14,
        );
        assert_eq!(load.wrench.torque_body, Vec3::ZERO);
        close(
            load.potential_energy_j,
            3.5 * sample.specific_potential_j_kg,
            2e-14,
        );
    }
}

#[test]
fn exterior_is_native_central_gravity_and_derivatives_match_the_same_potential() {
    let center = Vec3::new(0.5, -0.25, 1.0);
    let field = UniformSphereGravity::new(FRAME, center, 16.0, 2.0).unwrap();
    let point_source = CentralGravity::new(FRAME, center, 16.0).unwrap();
    for displacement in [Vec3::new(2.0, 0.0, 0.0), Vec3::new(3.0, -4.0, 0.0)] {
        assert_eq!(
            field.sample_point(FRAME, center.add(displacement)).unwrap(),
            point_source
                .sample_point(FRAME, center.add(displacement))
                .unwrap(),
        );
    }
    // Avoid the density jump, where no unique two-sided gravity tangent exists.
    for displacement in [Vec3::new(0.7, -0.4, 0.3), Vec3::new(3.1, 0.7, -0.2)] {
        let position = center.add(displacement);
        let sample = field.sample_point(FRAME, position).unwrap();
        let h = 2e-5;
        for column in 0..3 {
            let plus = field
                .sample_point(FRAME, position.add(axis(column, h)))
                .unwrap();
            let minus = field
                .sample_point(FRAME, position.sub(axis(column, h)))
                .unwrap();
            close(
                -(plus.specific_potential_j_kg - minus.specific_potential_j_kg) / (2.0 * h),
                values(sample.acceleration_world_m_s2)[column],
                8e-10,
            );
            let derivative = values(
                plus.acceleration_world_m_s2
                    .sub(minus.acceleration_world_m_s2)
                    .scale(0.5 / h),
            );
            for row in 0..3 {
                close(
                    derivative[row],
                    sample.acceleration_gradient_s2.m[3 * row + column],
                    8e-10,
                );
            }
        }
    }
}

#[test]
fn surface_continuity_and_tiny_interior_crossing_and_moving_source_work() {
    let field = UniformSphereGravity::new(FRAME, Vec3::ZERO, 1.0, 1.0).unwrap();
    let surface = field.sample_point(FRAME, Vec3::new(1.0, 0.0, 0.0)).unwrap();
    close(surface.specific_potential_j_kg, -1.0, 0.0);
    close(surface.acceleration_world_m_s2.x, -1.0, 0.0);
    close(surface.acceleration_gradient_s2.m[0], 2.0, 0.0);
    let inner = field
        .sample_point(FRAME, Vec3::new(1.0 - 1e-8, 0.0, 0.0))
        .unwrap();
    let outer = field
        .sample_point(FRAME, Vec3::new(1.0 + 1e-8, 0.0, 0.0))
        .unwrap();
    close(
        inner.specific_potential_j_kg,
        outer.specific_potential_j_kg,
        3e-8,
    );
    vector_close(
        inner.acceleration_world_m_s2,
        outer.acceleration_world_m_s2,
        3e-8,
    );
    close(inner.acceleration_gradient_s2.m[0], -1.0, 0.0);

    let transverse = 2.0_f64.powi(-30);
    let before = Vec3::new(0.25, 0.0, 0.0);
    let after = Vec3::new(0.25, transverse, 0.0);
    assert_eq!(
        field
            .sample_point(FRAME, before)
            .unwrap()
            .specific_potential_j_kg,
        field
            .sample_point(FRAME, after)
            .unwrap()
            .specific_potential_j_kg,
        "absolute potentials intentionally cannot represent the small work",
    );
    let expected = 0.5 * transverse * transverse;
    let change = field
        .specific_potential_increment(FRAME, before, after)
        .unwrap()
        .unwrap();
    close(change / expected, 1.0, 3e-15);
    let properties = mass(3.0);
    let pose = |p| Pose::new(p, UnitQuaternion::IDENTITY).unwrap();
    let rigid_change = field
        .potential_energy_increment(FRAME, pose(before), pose(after), properties)
        .unwrap()
        .unwrap();
    close(rigid_change / (3.0 * expected), 1.0, 3e-15);

    let just_inside = Vec3::new(1.0_f64.next_down(), 0.0, 0.0);
    let just_outside = Vec3::new(1.0_f64.next_up(), 0.0, 0.0);
    let inward = 1.0 - just_inside.x;
    let outward = just_outside.x - 1.0;
    let expected_crossing = inward * (1.0 - 0.5 * inward) + outward / (1.0 + outward);
    let crossing = field
        .specific_potential_increment(FRAME, just_inside, just_outside)
        .unwrap()
        .unwrap();
    close(crossing / expected_crossing, 1.0, 4e-15);
    close(
        field
            .specific_potential_increment(FRAME, just_outside, just_inside)
            .unwrap()
            .unwrap(),
        -crossing,
        0.0,
    );

    let source_before =
        UniformSphereGravity::new(FRAME, Vec3::new(1e10, 0.0, 0.0), 1.0, 1.0).unwrap();
    let source_after =
        UniformSphereGravity::new(FRAME, Vec3::new(1e10 + 1.0, 0.0, 0.0), 1.0, 1.0).unwrap();
    let moved_change = source_before
        .specific_potential_increment_to(
            &source_after,
            FRAME,
            Vec3::new(1e10 + 0.25, 0.0, 0.0),
            Vec3::new(1e10 + 1.25, transverse, 0.0),
        )
        .unwrap();
    close(moved_change / expected, 1.0, 3e-15);
}

fn options(tolerance: f64) -> GravityAdvanceOptions {
    GravityAdvanceOptions {
        initial_step_seconds: 0.5,
        minimum_step_seconds: 1e-10,
        maximum_step_seconds: 0.75,
        maximum_trials: 20_000,
        absolute_position_tolerance_m: tolerance,
        absolute_velocity_tolerance_m_s: tolerance,
        absolute_orientation_tolerance_rad: tolerance,
        absolute_angular_velocity_tolerance_rad_s: tolerance,
        relative_tolerance: 0.0,
    }
}

#[test]
fn native_adaptive_motion_crosses_the_center_and_tracks_the_harmonic_solution() {
    // mu/R^3 = 1. Starting at the center with v_x=0.75 gives x=0.75 sin(t).
    // The whole path stays inside the sphere and returns through its center.
    let field = UniformSphereGravity::new(FRAME, Vec3::ZERO, 8.0, 2.0).unwrap();
    let properties = mass(2.0);
    let initial = RigidBodyState::new(
        Pose::new(Vec3::ZERO, UnitQuaternion::IDENTITY).unwrap(),
        Vec3::new(1.5, 0.0, 0.0),
        Vec3::ZERO,
    )
    .unwrap();
    let duration = std::f64::consts::PI + 0.25;
    let solve = |tolerance| {
        advance_in_field(
            &field,
            FRAME,
            initial,
            properties,
            Wrench::ZERO,
            duration,
            options(tolerance),
            || false,
        )
        .unwrap()
    };
    let loose = solve(1e-5);
    let tight = solve(1e-8);
    let expected_position = Vec3::new(0.75 * duration.sin(), 0.0, 0.0);
    let expected_momentum = Vec3::new(1.5 * duration.cos(), 0.0, 0.0);
    let loose_error = norm(
        loose.states_after[0]
            .pose()
            .position_world()
            .sub(expected_position),
    );
    let tight_error = norm(
        tight.states_after[0]
            .pose()
            .position_world()
            .sub(expected_position),
    );
    assert!(
        tight_error < loose_error / 10.0,
        "loose={loose_error:e}, tight={tight_error:e}"
    );
    vector_close(
        tight.states_after[0].pose().position_world(),
        expected_position,
        3e-5,
    );
    vector_close(
        tight.states_after[0].linear_momentum_world(),
        expected_momentum,
        1e-4,
    );
    assert!(tight.stats.accepted_trials > loose.stats.accepted_trials);
    assert!(tight.stats.maximum_accepted_error_ratio <= 1.0);
    close(
        tight.diagnostics_after.mechanical_energy,
        tight.diagnostics_before.mechanical_energy,
        2e-5,
    );
    close(
        tight.potential_energy_change_j,
        tight.diagnostics_after.gravitational_potential_energy
            - tight.diagnostics_before.gravitational_potential_energy,
        3e-13,
    );
    assert_eq!(tight, solve(1e-8));
}

#[test]
fn scaled_sphere_arithmetic_and_source_frame_admission() {
    // All expected outputs are finite even when naive R^3 would overflow or underflow.
    for (radius, mu) in [(1e100, 1e100), (1e-100, 1e-300), (1e200, 1e300)] {
        let field = UniformSphereGravity::new(FRAME, Vec3::ZERO, mu, radius).unwrap();
        let sample = field
            .sample_point(FRAME, Vec3::new(0.5 * radius, 0.0, 0.0))
            .unwrap();
        let strength = mu / radius;
        close(
            sample.specific_potential_j_kg / (-1.375 * strength),
            1.0,
            5e-15,
        );
        close(
            sample.acceleration_world_m_s2.x / (-0.5 * strength / radius),
            1.0,
            5e-15,
        );
        close(
            sample.acceleration_gradient_s2.m[0] / (-strength / radius / radius),
            1.0,
            5e-15,
        );
    }
    for (mu, radius) in [
        (0.0, 1.0),
        (-1.0, 1.0),
        (f64::NAN, 1.0),
        (1.0, 0.0),
        (1.0, f64::INFINITY),
    ] {
        assert!(UniformSphereGravity::new(FRAME, Vec3::ZERO, mu, radius).is_err());
    }
    assert!(UniformSphereGravity::new(FRAME, Vec3::new(f64::NAN, 0.0, 0.0), 1.0, 1.0).is_err());
    let field = UniformSphereGravity::new(FRAME, Vec3::ZERO, 1.0, 1.0).unwrap();
    assert_eq!(
        field.sample_point(FRAME + 1, Vec3::ZERO),
        Err(GravityError::FrameMismatch)
    );
    assert!(
        field
            .sample_point(FRAME, Vec3::new(f64::INFINITY, 0.0, 0.0))
            .is_err()
    );
    let changed = UniformSphereGravity::new(FRAME, Vec3::ZERO, 2.0, 1.0).unwrap();
    assert!(
        field
            .specific_potential_increment_to(&changed, FRAME, Vec3::ZERO, Vec3::ZERO)
            .is_err()
    );
}
