//! Independent potential derivatives and closed-system tidal momentum exchange.
use fs_mbd::gravity::GravityError;
use fs_mbd::gravity::mutual::MutualGravity;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

const FRAME: u64 = 73;
const RADII: [f64; 2] = [0.6, 0.7];
const AXES: [Vec3; 3] = [
    Vec3::new(1.0, 0.0, 0.0),
    Vec3::new(0.0, 1.0, 0.0),
    Vec3::new(0.0, 0.0, 1.0),
];

fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}
fn vector_close(actual: Vec3, expected: Vec3, tolerance: f64) {
    close(norm(actual.sub(expected)), 0.0, tolerance);
}
fn fixture() -> ([RigidBodyState; 2], [MassProperties; 2]) {
    let body_axes = [
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 3.0), 0.6).unwrap(),
        UnitQuaternion::from_axis_angle(Vec3::new(-2.0, 1.0, 0.7), -0.4).unwrap(),
    ];
    let properties = [
        MassProperties::from_principal_frame(
            2.0,
            Vec3::ZERO,
            Vec3::new(0.07, 0.09, 0.12),
            body_axes[0],
        )
        .unwrap(),
        MassProperties::from_principal_frame(
            3.0,
            Vec3::ZERO,
            Vec3::new(0.10, 0.16, 0.20),
            body_axes[1],
        )
        .unwrap(),
    ];
    let states = [
        RigidBodyState::new(
            Pose::new(
                Vec3::new(-3.0, 0.2, -0.4),
                UnitQuaternion::from_axis_angle(Vec3::new(0.2, 1.0, -0.3), 0.5).unwrap(),
            )
            .unwrap(),
            Vec3::new(0.1, 0.3, -0.2),
            Vec3::new(0.04, -0.02, 0.07),
        )
        .unwrap(),
        RigidBodyState::new(
            Pose::new(
                Vec3::new(3.0, -0.4, 0.7),
                UnitQuaternion::from_axis_angle(Vec3::new(-0.5, 0.8, 1.0), -0.7).unwrap(),
            )
            .unwrap(),
            Vec3::new(-0.2, -0.1, 0.3),
            Vec3::new(-0.03, 0.01, 0.05),
        )
        .unwrap(),
    ];
    (states, properties)
}
fn with_pose(state: RigidBodyState, pose: Pose) -> RigidBodyState {
    RigidBodyState::new(
        pose,
        state.linear_momentum_world(),
        state.angular_momentum_body(),
    )
    .unwrap()
}
fn spin_world(states: &[RigidBodyState]) -> Vec3 {
    states.iter().fold(Vec3::ZERO, |sum, state| {
        sum.add(
            state
                .pose()
                .orientation()
                .rotate_body_to_world(state.angular_momentum_body()),
        )
    })
}

#[test]
fn both_forces_and_full_tensor_torques_differentiate_one_symmetric_potential() {
    let (states, properties) = fixture();
    let model = MutualGravity::new(FRAME, 0.8, 0.01, 1).unwrap()
        .quadrupole(&RADII, 0.25).unwrap();
    let load = model.load(FRAME, &states, &properties, || false).unwrap();
    assert_eq!(load.pair_count, 1);
    assert_eq!(load.wrenches[1].force_world, load.wrenches[0].force_world.scale(-1.0));
    let h = 2e-5;
    for body in 0..2 {
        assert!(norm(load.wrenches[body].torque_body) > 1e-5);
        for axis in AXES {
            let mut plus = states;
            let mut minus = states;
            for (candidate, sign) in [(&mut plus, 1.0), (&mut minus, -1.0)] {
                candidate[body] = with_pose(states[body], Pose::new(
                    states[body].pose().position_world().add(axis.scale(sign * h)),
                    states[body].pose().orientation(),
                ).unwrap());
            }
            let upper = model.load(FRAME, &plus, &properties, || false).unwrap().potential_energy_j;
            let lower = model.load(FRAME, &minus, &properties, || false).unwrap().potential_energy_j;
            close(-(upper - lower) / (2.0 * h), load.wrenches[body].force_world.dot(axis), 3e-9);

            for (candidate, sign) in [(&mut plus, 1.0), (&mut minus, -1.0)] {
                candidate[body] = with_pose(states[body], Pose::new(
                    states[body].pose().position_world(),
                    states[body].pose().orientation().right_exp(axis.scale(sign * h)).unwrap(),
                ).unwrap());
            }
            let upper = model.load(FRAME, &plus, &properties, || false).unwrap().potential_energy_j;
            let lower = model.load(FRAME, &minus, &properties, || false).unwrap().potential_energy_j;
            close(-(upper - lower) / (2.0 * h), load.wrenches[body].torque_body.dot(axis), 3e-9);
        }
    }
    let net_moment = states.iter().zip(&load.wrenches).fold(Vec3::ZERO, |sum, (s, w)| {
        sum.add(s.pose().position_world().cross(w.force_world))
            .add(s.pose().orientation().rotate_body_to_world(w.torque_body))
    });
    vector_close(net_moment, Vec3::ZERO, 3e-14);

    let reversed = [states[1], states[0]];
    let reversed_properties = [properties[1], properties[0]];
    let reversed_radii = [RADII[1], RADII[0]];
    let reordered = MutualGravity::new(FRAME, 0.8, 0.01, 1).unwrap()
        .quadrupole(&reversed_radii, 0.25).unwrap()
        .load(FRAME, &reversed, &reversed_properties, || false).unwrap();
    close(reordered.potential_energy_j, load.potential_energy_j, 2e-15);
    for i in 0..2 {
        vector_close(reordered.wrenches[i].force_world, load.wrenches[1 - i].force_world, 2e-15);
        vector_close(reordered.wrenches[i].torque_body, load.wrenches[1 - i].torque_body, 2e-15);
    }
}

#[test]
fn native_recoiling_steps_exchange_spin_and_orbital_momentum_and_reverse() {
    let (initial, properties) = fixture();
    let model = MutualGravity::new(FRAME, 0.8, 0.01, 1).unwrap()
        .quadrupole(&RADII, 0.25).unwrap();
    let dt = 0.005;
    let first = model.step(FRAME, &initial, &properties, dt, || false).unwrap();
    let before = first.diagnostics_before;
    let mut states = first.states_after;
    let mut after = first.diagnostics_after;
    for _ in 1..80 {
        let step = model.step(FRAME, &states, &properties, dt, || false).unwrap();
        states = step.states_after;
        after = step.diagnostics_after;
    }
    vector_close(after.linear_momentum_world, before.linear_momentum_world, 3e-12);
    vector_close(after.angular_momentum_world, before.angular_momentum_world, 3e-11);
    assert!(norm(spin_world(&states).sub(spin_world(&initial))) > 1e-6);
    assert!(states[0].linear_momentum_world() != initial[0].linear_momentum_world());
    assert!(states[1].linear_momentum_world() != initial[1].linear_momentum_world());
    let mut reversed: Vec<_> = states.iter().map(|s| {
        RigidBodyState::new(s.pose(), s.linear_momentum_world().scale(-1.0),
            s.angular_momentum_body().scale(-1.0)).unwrap()
    }).collect();
    for _ in 0..80 {
        reversed = model.step(FRAME, &reversed, &properties, dt, || false).unwrap().states_after;
    }
    for (actual, original) in reversed.iter().zip(initial) {
        vector_close(actual.pose().position_world(), original.pose().position_world(), 2e-10);
        vector_close(actual.linear_momentum_world().scale(-1.0), original.linear_momentum_world(), 2e-10);
        vector_close(actual.angular_momentum_body().scale(-1.0), original.angular_momentum_body(), 2e-10);
        for axis in AXES {
            vector_close(actual.pose().orientation().rotate_body_to_world(axis),
                original.pose().orientation().rotate_body_to_world(axis), 2e-10);
        }
    }
}

// Two independently resolved centrosymmetric clouds. Six equal point masses
// at +/- each declared principal axis determine the exact inertia and all 36
// pair potentials. Their leading omitted contribution is fourth mass moment.
fn cloud_potential_error(separation: f64) -> f64 {
    let (mut states, base) = fixture();
    let dimensions = [Vec3::new(0.4, 0.3, 0.2), Vec3::new(0.5, 0.2, 0.4)];
    let mut properties = base;
    let mut clouds = [[Vec3::ZERO; 6]; 2];
    for body in 0..2 {
        let a = dimensions[body];
        let m = base[body].mass();
        properties[body] = MassProperties::from_principal_frame(
            m, Vec3::ZERO,
            Vec3::new(a.y * a.y + a.z * a.z, a.x * a.x + a.z * a.z, a.x * a.x + a.y * a.y)
                .scale(m / 3.0),
            base[body].principal_axes_to_body(),
        ).unwrap();
        states[body] = with_pose(states[body], Pose::new(
            Vec3::new(if body == 0 { -0.5 * separation } else { 0.5 * separation }, 0.0, 0.0),
            states[body].pose().orientation(),
        ).unwrap());
        for (point, principal) in clouds[body].iter_mut().zip([
            Vec3::new(a.x, 0.0, 0.0), Vec3::new(-a.x, 0.0, 0.0),
            Vec3::new(0.0, a.y, 0.0), Vec3::new(0.0, -a.y, 0.0),
            Vec3::new(0.0, 0.0, a.z), Vec3::new(0.0, 0.0, -a.z),
        ]) {
            let body_point = base[body].principal_axes_to_body().rotate_body_to_world(principal);
            *point = states[body].pose().point_world_from_body(body_point).unwrap();
        }
    }
    let model = MutualGravity::new(FRAME, 0.8, 0.0, 1).unwrap();
    let monopole = model.load(FRAME, &states, &properties, || false).unwrap().potential_energy_j;
    let quadrupole = model.quadrupole(&RADII, 0.25).unwrap()
        .load(FRAME, &states, &properties, || false).unwrap().potential_energy_j;
    let mut exact = 0.0;
    for first in clouds[0] {
        for second in clouds[1] {
            exact -= 0.8 * properties[0].mass() * properties[1].mass()
                / (36.0 * norm(second.sub(first)));
        }
    }
    let error = (quadrupole - exact).abs();
    assert!(error < 0.02 * (monopole - exact).abs());
    error
}

#[test]
fn mutual_second_moments_match_independent_resolved_mass_clouds() {
    let near = cloud_potential_error(8.0);
    let far = cloud_potential_error(16.0);
    assert!((28.0..36.0).contains(&(near / far)), "{near:e}, {far:e}");
}

#[test]
fn combined_extent_drift_budget_and_late_cancellation_refuse_atomically() {
    let (states, properties) = fixture();
    let base = MutualGravity::new(FRAME, 0.8, 0.01, 1).unwrap();
    let model = base.quadrupole(&RADII, 0.25).unwrap();
    for ratio in [0.0, 1.0, -1.0, f64::NAN] {
        assert!(base.quadrupole(&RADII, ratio).is_err());
    }
    assert!(matches!(model.load(FRAME + 1, &states, &properties, || false), Err(GravityError::FrameMismatch)));
    assert!(base.quadrupole(&RADII[..1], 0.25).unwrap()
        .load(FRAME, &states, &properties, || false).is_err());
    assert!(MutualGravity::new(FRAME, 0.8, 0.01, 0).unwrap()
        .quadrupole(&RADII, 0.25).unwrap().load(FRAME, &states, &properties, || false).is_err());
    for radius in [0.0, -1.0, 0.001, f64::NAN] {
        assert!(base.quadrupole(&[radius, RADII[1]], 0.25).unwrap()
            .load(FRAME, &states, &properties, || false).is_err());
    }

    // Each radius/distance is below 0.8, but the declared spheres overlap.
    let close_states = [
        with_pose(states[0], Pose::new(Vec3::ZERO, states[0].pose().orientation()).unwrap()),
        with_pose(states[1], Pose::new(Vec3::new(1.0, 0.0, 0.0), states[1].pose().orientation()).unwrap()),
    ];
    assert!(matches!(base.quadrupole(&RADII, 0.8).unwrap()
        .load(FRAME, &close_states, &properties, || false), Err(GravityError::OutsideExpansionDomain)));

    // Separated endpoints alone would pass. The split drift passes at distance
    // one, above the monopole cutoff but outside the second-moment far field.
    let mut crossing = states;
    for i in 0..2 {
        let sign = if i == 0 { 1.0 } else { -1.0 };
        crossing[i] = RigidBodyState::new(
            Pose::new(Vec3::new(-4.0 * sign, 0.5 * sign, 0.0), states[i].pose().orientation()).unwrap(),
            Vec3::new(8.0 * sign * properties[i].mass(), 0.0, 0.0),
            Vec3::ZERO,
        ).unwrap();
    }
    let weak = MutualGravity::new(FRAME, 1e-12, 0.01, 1).unwrap();
    assert!(weak.step(FRAME, &crossing, &properties, 1.0, || false).is_ok());
    assert!(matches!(weak.quadrupole(&RADII, 0.25).unwrap()
        .step(FRAME, &crossing, &properties, 1.0, || false), Err(GravityError::OutsideExpansionDomain)));

    let mut polls = 0;
    model.step(FRAME, &states, &properties, 0.005, || { polls += 1; false }).unwrap();
    for stop in [1, polls / 2, polls] {
        let mut seen = 0;
        assert!(matches!(model.step(FRAME, &states, &properties, 0.005, || {
            seen += 1; seen == stop
        }), Err(GravityError::Cancelled)));
        assert_eq!(states, fixture().0);
    }
}
