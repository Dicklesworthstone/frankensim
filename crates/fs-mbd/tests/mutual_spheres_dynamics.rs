//! Finite-sphere recoil through overlap, native work and resolved trajectories.

use fs_mbd::gravity::mutual::MutualGravity;
use fs_mbd::gravity::{GravityAdvanceError, GravityAdvanceOptions, GravityError};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

const FRAME: u64 = 127;
const CONTAINED_RADII: [f64; 2] = [2.0, 0.25];

fn norm(value: Vec3) -> f64 {
    value.x.hypot(value.y).hypot(value.z)
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}",
    );
}

fn vector_close(actual: Vec3, expected: Vec3, tolerance: f64) {
    close(norm(actual.sub(expected)), 0.0, tolerance);
}

fn sphere_mass(mass: f64, radius: f64) -> MassProperties {
    let inertia = 0.4 * mass * radius * radius;
    MassProperties::new(mass, Vec3::ZERO, Vec3::new(inertia, inertia, inertia)).unwrap()
}

fn state(position: Vec3, momentum: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(position, UnitQuaternion::IDENTITY).unwrap(),
        momentum,
        Vec3::ZERO,
    )
    .unwrap()
}

fn contained_pair(relative_velocity: Vec3) -> ([RigidBodyState; 2], [MassProperties; 2]) {
    // G=1, m1+m2=8, A=2: relative motion obeys r''=-r while |r|<7/4.
    // The reduced mass is 15/8 and the barycenter is stationary at the origin.
    (
        [
            state(
                Vec3::new(-5.0 / 8.0, 0.0, 0.0),
                relative_velocity.scale(-15.0 / 8.0),
            ),
            state(
                Vec3::new(3.0 / 8.0, 0.0, 0.0),
                relative_velocity.scale(15.0 / 8.0),
            ),
        ],
        [sphere_mass(3.0, 2.0), sphere_mass(5.0, 0.25)],
    )
}

fn relative_motion(states: &[RigidBodyState]) -> (Vec3, Vec3) {
    (
        states[1]
            .pose()
            .position_world()
            .sub(states[0].pose().position_world()),
        states[1]
            .linear_momentum_world()
            .scale(1.0 / 5.0)
            .sub(states[0].linear_momentum_world().scale(1.0 / 3.0)),
    )
}

fn harmonic_error(states: &[RigidBodyState], initial_velocity: Vec3, time: f64) -> f64 {
    let (position, velocity) = relative_motion(states);
    let initial_position = Vec3::new(1.0, 0.0, 0.0);
    let expected_position = initial_position
        .scale(time.cos())
        .add(initial_velocity.scale(time.sin()));
    let expected_velocity = initial_velocity
        .scale(time.cos())
        .sub(initial_position.scale(time.sin()));
    norm(position.sub(expected_position)).hypot(norm(velocity.sub(expected_velocity)))
}

fn reverse(states: &[RigidBodyState]) -> Vec<RigidBodyState> {
    states
        .iter()
        .map(|state| {
            RigidBodyState::new(
                state.pose(),
                state.linear_momentum_world().scale(-1.0),
                state.angular_momentum_body().scale(-1.0),
            )
            .unwrap()
        })
        .collect()
}

fn options(duration: f64, tolerance: f64) -> GravityAdvanceOptions {
    GravityAdvanceOptions {
        initial_step_seconds: duration / 2.0,
        minimum_step_seconds: duration * 1e-9,
        maximum_step_seconds: duration / 2.0,
        maximum_trials: 10_000,
        absolute_position_tolerance_m: tolerance,
        absolute_velocity_tolerance_m_s: tolerance,
        absolute_orientation_tolerance_rad: tolerance,
        absolute_angular_velocity_tolerance_rad_s: tolerance,
        relative_tolerance: 0.0,
    }
}

#[test]
fn native_recoil_passes_through_coincident_centers_and_retains_explicit_cutoffs() {
    let (initial, properties) = contained_pair(Vec3::new(-0.5, 0.0, 0.0));
    let base = MutualGravity::new(FRAME, 1.0, 0.0, 1).unwrap();
    let model = base.uniform_spheres(&CONTAINED_RADII).unwrap();

    // Every translation and impulse here is an exact binary fraction. Both
    // finite centers coincide after the first kick/drift, where the force is zero.
    let first = model
        .step(FRAME, &initial, &properties, 1.0, || false)
        .unwrap();
    assert_eq!(first.pair_count, 1);
    for body in &first.states_after {
        assert_eq!(body.pose().position_world(), Vec3::ZERO);
    }
    assert_eq!(
        first.states_after[0].linear_momentum_world(),
        Vec3::new(15.0 / 8.0, 0.0, 0.0),
    );
    assert_eq!(
        first.states_after[1].linear_momentum_world(),
        Vec3::new(-15.0 / 8.0, 0.0, 0.0),
    );
    close(first.potential_energy_change_j, -15.0 / 16.0, 3e-15);
    assert_eq!(first.diagnostics_after.linear_momentum_world, Vec3::ZERO);

    let second = model
        .step(FRAME, &first.states_after, &properties, 1.0, || false)
        .unwrap();
    for (actual, before) in second.states_after.iter().zip(initial) {
        assert_eq!(
            actual.pose().position_world(),
            before.pose().position_world().scale(-1.0),
        );
        assert_eq!(
            actual.linear_momentum_world(),
            before.linear_momentum_world()
        );
    }
    close(second.potential_energy_change_j, 15.0 / 16.0, 3e-15);

    // A declared exclusion distance still binds even though the selected
    // source physics itself has no center singularity.
    let guarded = MutualGravity::new(FRAME, 1.0, 0.01, 1)
        .unwrap()
        .uniform_spheres(&CONTAINED_RADII)
        .unwrap();
    assert!(matches!(
        guarded.step(FRAME, &initial, &properties, 1.0, || false),
        Err(GravityError::SeparationCutoff { .. }),
    ));
    assert_eq!(initial, contained_pair(Vec3::new(-0.5, 0.0, 0.0)).0);
}

#[test]
fn contained_orbit_converges_at_second_order_and_preserves_total_momenta() {
    let initial_velocity = Vec3::new(0.0, 0.8, 0.2);
    let (initial, properties) = contained_pair(initial_velocity);
    let model = MutualGravity::new(FRAME, 1.0, 0.0, 1)
        .unwrap()
        .uniform_spheres(&CONTAINED_RADII)
        .unwrap();
    let mut errors = Vec::new();
    for steps in [32, 64] {
        let dt = 2.0 / f64::from(steps);
        let first = model
            .step(FRAME, &initial, &properties, dt, || false)
            .unwrap();
        let reference = first.diagnostics_before;
        let mut states = initial.to_vec();
        let mut largest_energy_error: f64 = 0.0;
        for _ in 0..steps {
            let receipt = model
                .step(FRAME, &states, &properties, dt, || false)
                .unwrap();
            vector_close(
                receipt.diagnostics_after.linear_momentum_world,
                reference.linear_momentum_world,
                2e-14,
            );
            vector_close(
                receipt.diagnostics_after.angular_momentum_world,
                reference.angular_momentum_world,
                3e-14,
            );
            largest_energy_error = largest_energy_error.max(
                (receipt.diagnostics_after.mechanical_energy - reference.mechanical_energy).abs(),
            );
            states = receipt.states_after;
        }
        assert!(largest_energy_error < 0.002);
        errors.push(harmonic_error(&states, initial_velocity, 2.0));
    }
    eprintln!("contained harmonic orbit errors: {errors:?}");
    assert!((3.9..4.1).contains(&(errors[0] / errors[1])));
    assert!(errors[1] < 0.001);
}

#[test]
fn direct_overlap_work_survives_unchanged_binding_energies_in_both_regimes() {
    let initial = [
        state(Vec3::new(-5e7, 0.0, 0.0), Vec3::new(0.0, -7500.0, 0.0)),
        state(Vec3::new(5e7, 0.0, 0.0), Vec3::new(0.0, 7500.0, 0.0)),
    ];
    // The relative transverse drift is exactly 5/8 m while both rounded
    // x coordinates and the rounded center distance remain unchanged.
    // Containment has exactly quadratic potential. At equal radii and d=R,
    // F/(G*m1*m2/R^2)=15/32; the first omitted transverse term is <1e-29 J.
    for (radii, expected) in [([2e8, 1e7], 1.46484375e-13), ([1e8, 1e8], 5.4931640625e-13)] {
        let base = MutualGravity::new(FRAME, 1e12, 0.0, 1).unwrap();
        let model = base.uniform_spheres(&radii).unwrap();
        let properties = [sphere_mass(2.0, radii[0]), sphere_mass(3.0, radii[1])];
        for shift in [Vec3::ZERO, Vec3::new(2e9, 1024.0, -512.0)] {
            let states = initial.map(|body| {
                state(
                    body.pose().position_world().add(shift),
                    body.linear_momentum_world(),
                )
            });
            let receipt = model
                .step(FRAME, &states, &properties, 1e-4, || false)
                .unwrap();
            assert_eq!(
                receipt.diagnostics_before.gravitational_potential_energy,
                receipt.diagnostics_after.gravitational_potential_energy,
            );
            close(receipt.potential_energy_change_j, expected, 2e-26);
            for (actual, before) in receipt.states_after.iter().zip(states) {
                assert_eq!(
                    actual.pose().position_world().x,
                    before.pose().position_world().x,
                );
            }
            let reversed_radii = [radii[1], radii[0]];
            let reordered = base
                .uniform_spheres(&reversed_radii)
                .unwrap()
                .step(
                    FRAME,
                    &[states[1], states[0]],
                    &[properties[1], properties[0]],
                    1e-4,
                    || false,
                )
                .unwrap();
            assert_eq!(
                receipt.potential_energy_change_j,
                reordered.potential_energy_change_j
            );
            assert_eq!(receipt.states_after[0], reordered.states_after[1]);
            assert_eq!(receipt.states_after[1], reordered.states_after[0]);

            let returned = model
                .step(
                    FRAME,
                    &reverse(&receipt.states_after),
                    &properties,
                    1e-4,
                    || false,
                )
                .unwrap();
            close(returned.potential_energy_change_j, -expected, 2e-26);
        }
    }
}

#[test]
fn direct_work_crosses_containment_partial_overlap_and_exterior_boundaries() {
    let model = MutualGravity::new(FRAME, 1.0, 0.0, 1)
        .unwrap()
        .uniform_spheres(&CONTAINED_RADII)
        .unwrap();
    for (speed, distance_after) in [(1.5, 2.0), (2.5, 3.0)] {
        let (initial, properties) = contained_pair(Vec3::new(speed, 0.0, 0.0));
        let receipt = model
            .step(FRAME, &initial, &properties, 1.0, || false)
            .unwrap();
        assert_eq!(relative_motion(&receipt.states_after).0.x, distance_after);
        // These changes are well above roundoff. The static law has a separate
        // shell-volume oracle; its endpoint difference checks every intervening
        // segment of the direct-work path without reusing the work algorithm.
        let endpoint_change = receipt.diagnostics_after.gravitational_potential_energy
            - receipt.diagnostics_before.gravitational_potential_energy;
        close(receipt.potential_energy_change_j, endpoint_change, 3e-14);
        if distance_after == 3.0 {
            // U_before=-2631/256 and U_after=-G*m1*m2/3=-5 exactly.
            close(receipt.potential_energy_change_j, 5.27734375, 3e-14);
        }
        let returned = model
            .step(
                FRAME,
                &reverse(&receipt.states_after),
                &properties,
                1.0,
                || false,
            )
            .unwrap();
        close(returned.potential_energy_change_j, -endpoint_change, 3e-14);
        for (actual, expected) in returned.states_after.iter().zip(initial) {
            vector_close(
                actual.pose().position_world(),
                expected.pose().position_world(),
                3e-15,
            );
            vector_close(
                actual.linear_momentum_world().scale(-1.0),
                expected.linear_momentum_world(),
                3e-15,
            );
        }
    }
}

#[test]
fn partially_overlapping_three_body_motion_keeps_native_spin_and_is_reversible() {
    let radii = [1.0, 0.7, 0.5];
    let properties = [
        sphere_mass(2.0, 1.0),
        sphere_mass(3.0, 0.7),
        sphere_mass(5.0, 0.5),
    ];
    let states = [
        state(Vec3::new(-0.8, 0.1, 0.2), Vec3::new(0.2, -0.1, 0.3)),
        state(Vec3::new(0.5, -0.2, -0.3), Vec3::new(-0.3, 0.2, 0.1)),
        state(Vec3::new(0.1, 0.8, 0.2), Vec3::new(0.1, -0.1, -0.4)),
    ]
    .map(|body| {
        RigidBodyState::new(
            body.pose(),
            body.linear_momentum_world(),
            Vec3::new(0.03, -0.02, 0.01),
        )
        .unwrap()
    });
    let model = MutualGravity::new(FRAME, 0.7, 0.0, 3)
        .unwrap()
        .uniform_spheres(&radii)
        .unwrap();
    let receipt = model
        .step(FRAME, &states, &properties, 0.03, || false)
        .unwrap();
    assert_eq!(receipt.pair_count, 3);
    vector_close(
        receipt.diagnostics_after.linear_momentum_world,
        receipt.diagnostics_before.linear_momentum_world,
        2e-15,
    );
    vector_close(
        receipt.diagnostics_after.angular_momentum_world,
        receipt.diagnostics_before.angular_momentum_world,
        3e-15,
    );
    for i in 0..states.len() {
        let free = RigidBodyIntegrator::new(Gravity::ZERO)
            .step_symplectic(states[i], properties[i], Wrench::ZERO, 0.03)
            .unwrap()
            .state_after;
        assert_eq!(
            receipt.states_after[i].pose().orientation(),
            free.pose().orientation()
        );
        assert_eq!(
            receipt.states_after[i].angular_momentum_body(),
            free.angular_momentum_body()
        );
    }
    let returned = model
        .step(
            FRAME,
            &reverse(&receipt.states_after),
            &properties,
            0.03,
            || false,
        )
        .unwrap();
    for (actual, expected) in returned.states_after.iter().zip(states) {
        vector_close(
            actual.pose().position_world(),
            expected.pose().position_world(),
            3e-15,
        );
        vector_close(
            actual.linear_momentum_world().scale(-1.0),
            expected.linear_momentum_world(),
            3e-15,
        );
        vector_close(
            actual.angular_momentum_body().scale(-1.0),
            expected.angular_momentum_body(),
            3e-15,
        );
        vector_close(
            actual
                .pose()
                .orientation()
                .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
            Vec3::new(1.0, 0.0, 0.0),
            3e-15,
        );
    }
}

#[test]
fn adaptive_containment_resolves_center_crossing_with_refining_tolerances() {
    let velocity = Vec3::new(-0.5, 0.0, 0.0);
    let (states, properties) = contained_pair(velocity);
    let model = MutualGravity::new(FRAME, 1.0, 0.0, 1)
        .unwrap()
        .uniform_spheres(&CONTAINED_RADII)
        .unwrap();
    let mut errors = Vec::new();
    let mut accepted = Vec::new();
    for tolerance in [1e-4, 1e-6, 1e-8] {
        let result = model
            .advance(
                FRAME,
                &states,
                &properties,
                2.0,
                options(2.0, tolerance),
                || false,
            )
            .unwrap();
        assert_eq!(result.duration_seconds, 2.0);
        assert!(result.stats.rejected_trials > 0);
        assert_eq!(
            result.stats.trials,
            result.stats.accepted_trials + result.stats.rejected_trials
        );
        assert_eq!(result.stats.native_steps, 3 * result.stats.trials);
        assert!(result.stats.maximum_accepted_error_ratio <= 1.0);
        vector_close(
            result.diagnostics_after.linear_momentum_world,
            Vec3::ZERO,
            3e-14,
        );
        let (position, _) = relative_motion(&result.states_after);
        close(
            result.potential_energy_change_j,
            15.0 / 16.0 * (position.dot(position) - 1.0),
            8e-14,
        );
        errors.push(harmonic_error(&result.states_after, velocity, 2.0));
        accepted.push(result.stats.accepted_trials);
    }
    eprintln!("adaptive containment errors: {errors:?}; accepted trials: {accepted:?}");
    assert!(errors[1] < errors[0] / 8.0);
    assert!(errors[2] < errors[1] / 8.0);
    assert!(errors[2] < 2e-5);
    assert!(accepted[0] < accepted[1] && accepted[1] < accepted[2]);
}

#[test]
fn interval_work_limits_and_late_cancellation_publish_no_partial_candidate() {
    let (states, properties) = contained_pair(Vec3::new(-0.5, 0.0, 0.0));
    let original = states;
    let model = MutualGravity::new(FRAME, 1.0, 0.0, 1)
        .unwrap()
        .uniform_spheres(&CONTAINED_RADII)
        .unwrap();
    let mut limited = options(2.0, 1e-12);
    limited.maximum_trials = 1;
    assert_eq!(
        model.advance(FRAME, &states, &properties, 2.0, limited, || false),
        Err(GravityAdvanceError::TrialLimit),
    );
    let mut too_coarse = options(2.0, 1e-12);
    too_coarse.minimum_step_seconds = 1.0;
    assert_eq!(
        model.advance(FRAME, &states, &properties, 2.0, too_coarse, || false),
        Err(GravityAdvanceError::MinimumStepReached),
    );
    let mut polls = 0;
    model
        .advance(FRAME, &states, &properties, 0.5, options(0.5, 1e-5), || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 1);
    let mut seen = 0;
    assert_eq!(
        model.advance(FRAME, &states, &properties, 0.5, options(0.5, 1e-5), || {
            seen += 1;
            seen == polls
        }),
        Err(GravityAdvanceError::Gravity(GravityError::Cancelled)),
    );
    assert_eq!(states, original);
}
