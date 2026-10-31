//! G1/G2/G3: recoiling Newtonian systems, analytic Kepler motion and exact refusal.
use fs_mbd::gravity::GravityError;
use fs_mbd::gravity::mutual::MutualGravity;
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

const FRAME: u64 = 41;

fn norm(value: Vec3) -> f64 {
    value.x.hypot(value.y).hypot(value.z)
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}
fn vector_close(actual: Vec3, expected: Vec3, tolerance: f64) {
    close(norm(actual.sub(expected)), 0.0, tolerance);
}
fn mass(value: f64) -> MassProperties {
    MassProperties::new(value, Vec3::ZERO, Vec3::new(0.2, 0.3, 0.4)).unwrap()
}
fn state(position: Vec3, momentum: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(position, UnitQuaternion::IDENTITY).unwrap(),
        momentum,
        Vec3::ZERO,
    )
    .unwrap()
}
fn three_bodies() -> ([RigidBodyState; 3], [MassProperties; 3]) {
    (
        [
            state(Vec3::new(-2.0, 0.0, 0.3), Vec3::new(0.2, -0.1, 0.4)),
            state(Vec3::new(1.2, -0.6, 0.7), Vec3::new(-0.3, 0.2, 0.1)),
            state(Vec3::new(0.3, 2.0, -0.4), Vec3::new(0.1, -0.2, -0.5)),
        ],
        [mass(1.0), mass(2.0), mass(4.0)],
    )
}

#[test]
fn finite_sources_recoil_and_pair_potential_is_counted_once() {
    let model = MutualGravity::new(FRAME, 3.0, 0.0, 1).unwrap();
    let states = [
        state(Vec3::new(-1.0, 0.0, 0.0), Vec3::ZERO),
        state(Vec3::new(2.0, 0.0, 0.0), Vec3::ZERO),
    ];
    let properties = [mass(2.0), mass(5.0)];
    let load = model.load(FRAME, &states, &properties, || false).unwrap();
    assert_eq!(load.pair_count, 1);
    close(load.potential_energy_j, -10.0, 1e-14);
    vector_close(load.forces_world[0], Vec3::new(10.0 / 3.0, 0.0, 0.0), 1e-14);
    assert_eq!(load.forces_world[1], load.forces_world[0].scale(-1.0));
    close(
        load.forces_world[0].x / properties[0].mass(),
        5.0 / 3.0,
        1e-14,
    );
    close(
        load.forces_world[1].x / properties[1].mass(),
        -2.0 / 3.0,
        1e-14,
    );

    let receipt = model
        .step(FRAME, &states, &properties, 0.01, || false)
        .unwrap();
    assert!(receipt.states_after[0].linear_momentum_world().x > 0.0);
    assert!(receipt.states_after[1].linear_momentum_world().x < 0.0);
    vector_close(
        receipt.diagnostics_after.linear_momentum_world,
        Vec3::ZERO,
        1e-15,
    );
    close(
        receipt.diagnostics_before.gravitational_potential_energy,
        -10.0,
        1e-14,
    );
}

#[test]
fn three_body_force_differentiates_one_potential_and_has_no_internal_moment() {
    let (states, properties) = three_bodies();
    let model = MutualGravity::new(FRAME, 0.7, 0.01, 3).unwrap();
    let load = model.load(FRAME, &states, &properties, || false).unwrap();
    let axes = [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ];
    let h = 1e-5;
    for i in 0..states.len() {
        for axis in axes {
            let mut plus = states;
            let mut minus = states;
            plus[i] = state(
                states[i].pose().position_world().add(axis.scale(h)),
                states[i].linear_momentum_world(),
            );
            minus[i] = state(
                states[i].pose().position_world().sub(axis.scale(h)),
                states[i].linear_momentum_world(),
            );
            let upper = model
                .load(FRAME, &plus, &properties, || false)
                .unwrap()
                .potential_energy_j;
            let lower = model
                .load(FRAME, &minus, &properties, || false)
                .unwrap()
                .potential_energy_j;
            close(
                -(upper - lower) / (2.0 * h),
                load.forces_world[i].dot(axis),
                7e-11,
            );
        }
    }
    let force = load
        .forces_world
        .iter()
        .fold(Vec3::ZERO, |sum, &f| sum.add(f));
    let moment = states
        .iter()
        .zip(&load.forces_world)
        .fold(Vec3::ZERO, |sum, (s, &f)| {
            sum.add(s.pose().position_world().cross(f))
        });
    vector_close(force, Vec3::ZERO, 2e-15);
    vector_close(moment, Vec3::ZERO, 2e-15);
}

// The two-body relative orbit obeys r'' = -G(m1+m2) r/|r|³, not -G*m1 r/|r|³.
// Analytic apocenter and period therefore test source recoil as well as order.
fn eccentric_binary() -> (MutualGravity, [RigidBodyState; 2], [MassProperties; 2], f64) {
    let properties = [mass(2.0), mass(5.0)];
    let model = MutualGravity::new(FRAME, 0.3, 0.1, 1).unwrap();
    let total = 7.0;
    let semi_major: f64 = 2.0;
    let eccentricity: f64 = 0.35;
    let radius = semi_major * (1.0 - eccentricity);
    let speed = (0.3 * total * (1.0 + eccentricity) / radius).sqrt();
    let reduced_mass = 10.0 / total;
    let states = [
        state(
            Vec3::new(-5.0 / total * radius, 0.0, 0.0),
            Vec3::new(0.0, -reduced_mass * speed, 0.0),
        ),
        state(
            Vec3::new(2.0 / total * radius, 0.0, 0.0),
            Vec3::new(0.0, reduced_mass * speed, 0.0),
        ),
    ];
    let period = 2.0 * core::f64::consts::PI * (semi_major.powi(3) / (0.3 * total)).sqrt();
    (model, states, properties, period)
}

#[test]
fn binary_kepler_apocenter_converges_at_second_order() {
    let (model, initial, properties, period) = eccentric_binary();
    let solve = |steps: usize| {
        let mut states = initial.to_vec();
        for _ in 0..steps {
            states = model
                .step(
                    FRAME,
                    &states,
                    &properties,
                    0.5 * period / steps as f64,
                    || false,
                )
                .unwrap()
                .states_after;
        }
        let relative = states[1]
            .pose()
            .position_world()
            .sub(states[0].pose().position_world());
        let barycenter = states[0]
            .pose()
            .position_world()
            .scale(2.0 / 7.0)
            .add(states[1].pose().position_world().scale(5.0 / 7.0));
        vector_close(barycenter, Vec3::ZERO, 3e-14);
        norm(relative.sub(Vec3::new(-2.7, 0.0, 0.0)))
    };
    let coarse = solve(160);
    let fine = solve(320);
    eprintln!(
        "binary apocenter errors: {coarse:e}, {fine:e}; ratio {}",
        coarse / fine
    );
    assert!((3.8..4.2).contains(&(coarse / fine)));
    assert!(fine < 0.001);
}

#[test]
fn resolved_binary_energy_stays_bounded_and_total_momenta_close() {
    let (model, initial, properties, period) = eccentric_binary();
    let dt = period / 320.0;
    let first = model
        .step(FRAME, &initial, &properties, dt, || false)
        .unwrap();
    let reference = first.diagnostics_before;
    close(reference.mechanical_energy, -0.75, 3e-15);
    let mut states = initial.to_vec();
    let mut largest_relative_energy_error: f64 = 0.0;
    for _ in 0..(20 * 320) {
        let receipt = model
            .step(FRAME, &states, &properties, dt, || false)
            .unwrap();
        largest_relative_energy_error = largest_relative_energy_error.max(
            (receipt.diagnostics_after.mechanical_energy / reference.mechanical_energy - 1.0).abs(),
        );
        vector_close(
            receipt.diagnostics_after.linear_momentum_world,
            reference.linear_momentum_world,
            4e-14,
        );
        vector_close(
            receipt.diagnostics_after.angular_momentum_world,
            reference.angular_momentum_world,
            5e-13,
        );
        states = receipt.states_after;
    }
    eprintln!("20-orbit maximum relative energy error: {largest_relative_energy_error:e}");
    assert!(largest_relative_energy_error < 0.001);
}

#[test]
fn three_body_motion_is_covariant_reversible_and_preserves_native_spin() {
    let (mut initial, mut properties) = three_bodies();
    let axes = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -0.7), 0.6).unwrap();
    for i in 0..3 {
        properties[i] = MassProperties::from_principal_frame(
            properties[i].mass(),
            Vec3::ZERO,
            Vec3::new(0.2, 0.3, 0.4),
            axes,
        )
        .unwrap();
        initial[i] = RigidBodyState::new(
            initial[i].pose(),
            initial[i].linear_momentum_world(),
            Vec3::new(0.03, -0.07, 0.05),
        )
        .unwrap();
    }
    let model = MutualGravity::new(FRAME, 0.7, 0.01, 3).unwrap();
    let dt = 0.03;
    let reference = model
        .step(FRAME, &initial, &properties, dt, || false)
        .unwrap();
    vector_close(
        reference.diagnostics_after.linear_momentum_world,
        reference.diagnostics_before.linear_momentum_world,
        2e-15,
    );
    vector_close(
        reference.diagnostics_after.angular_momentum_world,
        reference.diagnostics_before.angular_momentum_world,
        3e-15,
    );
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(-0.2, 0.5, 0.7), 0.8).unwrap();
    let offset = Vec3::new(4.0, -3.0, 0.7);
    let boost = Vec3::new(0.3, -0.4, 0.2);
    let transformed: Vec<_> = initial
        .iter()
        .zip(properties)
        .map(|(s, m)| {
            RigidBodyState::new(
                Pose::new(
                    rotation
                        .rotate_body_to_world(s.pose().position_world())
                        .add(offset),
                    rotation,
                )
                .unwrap(),
                rotation
                    .rotate_body_to_world(s.linear_momentum_world())
                    .add(boost.scale(m.mass())),
                s.angular_momentum_body(),
            )
            .unwrap()
        })
        .collect();
    let moved = model
        .step(FRAME, &transformed, &properties, dt, || false)
        .unwrap();
    for i in 0..3 {
        let expected = reference.states_after[i];
        let actual = moved.states_after[i];
        vector_close(
            actual.pose().position_world(),
            rotation
                .rotate_body_to_world(expected.pose().position_world())
                .add(offset)
                .add(boost.scale(dt)),
            3e-15,
        );
        vector_close(
            actual.linear_momentum_world(),
            rotation
                .rotate_body_to_world(expected.linear_momentum_world())
                .add(boost.scale(properties[i].mass())),
            2e-15,
        );
        vector_close(
            actual.angular_momentum_body(),
            expected.angular_momentum_body(),
            2e-15,
        );
        vector_close(
            actual
                .pose()
                .orientation()
                .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
            rotation.rotate_body_to_world(
                expected
                    .pose()
                    .orientation()
                    .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
            ),
            2e-15,
        );
        let free = RigidBodyIntegrator::new(Gravity::ZERO)
            .step_symplectic(initial[i], properties[i], Wrench::ZERO, dt)
            .unwrap()
            .state_after;
        assert_eq!(free.pose().orientation(), expected.pose().orientation());
        assert_eq!(
            free.angular_momentum_body(),
            expected.angular_momentum_body()
        );
    }
    let permutation = [2, 0, 1];
    let permuted = permutation.map(|i| initial[i]);
    let permuted_properties = permutation.map(|i| properties[i]);
    let reordered = model
        .step(FRAME, &permuted, &permuted_properties, dt, || false)
        .unwrap();
    for (slot, &original) in permutation.iter().enumerate() {
        vector_close(
            reordered.states_after[slot].pose().position_world(),
            reference.states_after[original].pose().position_world(),
            1e-15,
        );
        vector_close(
            reordered.states_after[slot].linear_momentum_world(),
            reference.states_after[original].linear_momentum_world(),
            1e-15,
        );
    }
    let reversed: Vec<_> = reference
        .states_after
        .iter()
        .map(|s| {
            RigidBodyState::new(
                s.pose(),
                s.linear_momentum_world().scale(-1.0),
                s.angular_momentum_body().scale(-1.0),
            )
            .unwrap()
        })
        .collect();
    let returned = model
        .step(FRAME, &reversed, &properties, dt, || false)
        .unwrap();
    for (actual, expected) in returned.states_after.iter().zip(initial) {
        vector_close(
            actual.pose().position_world(),
            expected.pose().position_world(),
            1e-15,
        );
        vector_close(
            actual.linear_momentum_world().scale(-1.0),
            expected.linear_momentum_world(),
            1e-15,
        );
        vector_close(
            actual.angular_momentum_body().scale(-1.0),
            expected.angular_momentum_body(),
            1e-15,
        );
        vector_close(
            actual
                .pose()
                .orientation()
                .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
            Vec3::new(1.0, 0.0, 0.0),
            2e-15,
        );
    }
}

#[test]
fn empty_single_body_and_extreme_radius_queries_use_original_owners() {
    let model = MutualGravity::new(FRAME, 1.0, 0.0, 0).unwrap();
    let empty = model.step(FRAME, &[], &[], 0.1, || false).unwrap();
    assert!(empty.states_after.is_empty());
    assert_eq!(empty.pair_count, 0);
    assert_eq!(empty.diagnostics_after.mechanical_energy, 0.0);
    let initial = state(Vec3::new(1.0, 2.0, 3.0), Vec3::new(0.3, -0.2, 0.1));
    let expected = RigidBodyIntegrator::new(Gravity::ZERO)
        .step_symplectic(initial, mass(2.0), Wrench::ZERO, 0.1)
        .unwrap();
    let actual = model
        .step(FRAME, &[initial], &[mass(2.0)], 0.1, || false)
        .unwrap();
    assert_eq!(actual.states_after, vec![expected.state_after]);
    assert_eq!(actual.diagnostics_after, expected.diagnostics_after);
    let model = MutualGravity::new(FRAME, 1.0, 0.0, 1).unwrap();
    for (body_mass, radius, expected_force) in [(1e100, 1e200, 1e-200), (1e-100, 1e-200, 1e200)] {
        let bodies = [
            state(Vec3::ZERO, Vec3::ZERO),
            state(Vec3::new(radius, 0.0, 0.0), Vec3::ZERO),
        ];
        let load = model
            .load(FRAME, &bodies, &[mass(body_mass); 2], || false)
            .unwrap();
        close(load.potential_energy_j, -1.0, 3e-15);
        close(load.forces_world[0].x / expected_force, 1.0, 3e-15);
    }

    // Valid forces/energies must not disappear merely because an intermediate
    // G*m underflows or overflows. Reordering the same bodies cannot alter
    // admission or which physical mass receives the reaction.
    let bodies = [
        state(Vec3::ZERO, Vec3::ZERO),
        state(Vec3::new(1.0, 0.0, 0.0), Vec3::ZERO),
    ];
    for (g, first_mass, second_mass, expected) in [
        (1e-200, 1e-200, 1e200, 1e-200),
        (1e200, 1e200, 1e-200, 1e200),
        (f64::from_bits(1), 1.0, 1.0, f64::from_bits(1)),
    ] {
        let model = MutualGravity::new(FRAME, g, 0.0, 1).unwrap();
        let properties = [mass(first_mass), mass(second_mass)];
        let load = model.load(FRAME, &bodies, &properties, || false).unwrap();
        let reversed = model
            .load(
                FRAME,
                &[bodies[1], bodies[0]],
                &[properties[1], properties[0]],
                || false,
            )
            .unwrap();
        close(load.potential_energy_j / -expected, 1.0, 3e-15);
        close(load.forces_world[0].x / expected, 1.0, 3e-15);
        assert_eq!(load.potential_energy_j, reversed.potential_energy_j);
        assert_eq!(load.forces_world[0], reversed.forces_world[1]);
        assert_eq!(load.forces_world[1], reversed.forces_world[0]);
    }

    let model = MutualGravity::new(FRAME, 1e100, 0.0, 1).unwrap();
    let mixed = [
        state(Vec3::ZERO, Vec3::ZERO),
        state(Vec3::new(1e200, 1e-200, 0.0), Vec3::ZERO),
    ];
    let load = model
        .load(FRAME, &mixed, &[mass(1e200); 2], || false)
        .unwrap();
    close(load.forces_world[0].x / 1e100, 1.0, 3e-15);
    close(load.forces_world[0].y / 1e-300, 1.0, 3e-15);
    assert_eq!(load.forces_world[1], load.forces_world[0].scale(-1.0));
}

#[test]
fn separation_guard_catches_interior_source_crossing_with_separated_endpoints() {
    let model = MutualGravity::new(FRAME, 1e-8, 1e-6, 1).unwrap();
    let initial = [
        state(Vec3::new(-1.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)),
        state(Vec3::new(1.0, 0.0, 0.0), Vec3::new(-2.0, 0.0, 0.0)),
    ];
    let properties = [mass(1.0); 2];
    assert!(model.load(FRAME, &initial, &properties, || false).is_ok());
    // The free-drift endpoints swap sides and are again separated by 2 m.
    assert!(
        norm(
            initial[1]
                .pose()
                .position_world()
                .add(initial[1].linear_momentum_world())
                .sub(
                    initial[0]
                        .pose()
                        .position_world()
                        .add(initial[0].linear_momentum_world())
                )
        ) > 1.9
    );
    assert_eq!(
        model.step(FRAME, &initial, &properties, 1.0, || false),
        Err(GravityError::Invalid(
            "mutual gravity separation cutoff reached"
        ))
    );

    // The closest interpolation at t=1/3 leaves a false 1.1e-16 gap unless
    // the interior line-distance calculation preserves exact collinearity.
    let unsoftened = MutualGravity::new(FRAME, 1e-30, 0.0, 1).unwrap();
    let crossing = [
        state(Vec3::ZERO, Vec3::ZERO),
        state(Vec3::new(1.0, 0.0, 0.0), Vec3::new(-3.0, 0.0, 0.0)),
    ];
    assert_eq!(
        unsoftened.step(FRAME, &crossing, &properties, 1.0, || false),
        Err(GravityError::Invalid(
            "mutual gravity separation cutoff reached"
        ))
    );
}

#[test]
fn budgets_frames_singularities_arithmetic_and_late_cancellation_refuse_atomically() {
    let (states, properties) = three_bodies();
    let saved = states;
    let model = MutualGravity::new(FRAME, 0.7, 0.01, 3).unwrap();
    assert_eq!(
        model.load(FRAME + 1, &states, &properties, || false),
        Err(GravityError::FrameMismatch)
    );
    assert!(
        model
            .load(FRAME, &states, &properties[..2], || false)
            .is_err()
    );
    assert!(
        MutualGravity::new(FRAME, 0.7, 0.01, 2)
            .unwrap()
            .load(FRAME, &states, &properties, || false)
            .is_err()
    );
    assert!(
        model
            .load(FRAME, &[states[0]; 2], &properties[..2], || false)
            .is_err()
    );
    assert!(
        MutualGravity::new(FRAME, f64::MAX, 0.0, 3)
            .unwrap()
            .load(FRAME, &states, &properties, || false)
            .is_err()
    );
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(MutualGravity::new(FRAME, invalid, 0.01, 3).is_err());
        assert!(
            model
                .step(FRAME, &states, &properties, invalid, || false)
                .is_err()
        );
    }
    assert!(MutualGravity::new(FRAME, 1.0, -0.1, 3).is_err());
    let mut polls = 0;
    let expected = model
        .step(FRAME, &states, &properties, 0.01, || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut calls = 0;
        let refused = model.step(FRAME, &states, &properties, 0.01, || {
            calls += 1;
            calls == stop
        });
        assert_eq!(refused, Err(GravityError::Cancelled));
        assert_eq!(states, saved);
    }
    assert_eq!(
        model
            .step(FRAME, &states, &properties, 0.01, || false)
            .unwrap(),
        expected
    );
}
