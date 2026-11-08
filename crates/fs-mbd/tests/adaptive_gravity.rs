//! Adaptive native gravity against Kepler motion and resolved tidal exchange.
use fs_mbd::gravity::mutual::MutualGravity;
use fs_mbd::gravity::{GravityAdvanceError, GravityAdvanceOptions, GravityError};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

const FRAME: u64 = 97;
const AXES: [Vec3; 3] = [
    Vec3::new(1.0, 0.0, 0.0),
    Vec3::new(0.0, 1.0, 0.0),
    Vec3::new(0.0, 0.0, 1.0),
];

fn norm(value: Vec3) -> f64 {
    value.x.hypot(value.y).hypot(value.z)
}

fn vector_close(actual: Vec3, expected: Vec3, tolerance: f64) {
    assert!(
        norm(actual.sub(expected)) <= tolerance,
        "{actual:?} != {expected:?}; tolerance {tolerance:e}"
    );
}

fn options(duration: f64, tolerance: f64) -> GravityAdvanceOptions {
    GravityAdvanceOptions {
        initial_step_seconds: duration / 4.0,
        minimum_step_seconds: duration * 1e-10,
        maximum_step_seconds: duration / 4.0,
        maximum_trials: 50_000,
        absolute_position_tolerance_m: tolerance,
        absolute_velocity_tolerance_m_s: tolerance,
        absolute_orientation_tolerance_rad: tolerance,
        absolute_angular_velocity_tolerance_rad_s: tolerance,
        relative_tolerance: 0.0,
    }
}

fn body(position: Vec3, momentum: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(position, UnitQuaternion::IDENTITY).unwrap(),
        momentum,
        Vec3::ZERO,
    )
    .unwrap()
}

fn eccentric_binary() -> (MutualGravity, [RigidBodyState; 2], [MassProperties; 2], f64) {
    let masses = [2.0, 5.0];
    let properties =
        masses.map(|m| MassProperties::new(m, Vec3::ZERO, Vec3::new(0.2, 0.3, 0.4)).unwrap());
    let semimajor: f64 = 2.0;
    let eccentricity: f64 = 0.65;
    let mu: f64 = 0.3 * 7.0;
    let radius = semimajor * (1.0 - eccentricity);
    let speed = (mu * (1.0 + eccentricity) / radius).sqrt();
    let momentum = (10.0 / 7.0) * speed;
    let states = [
        body(
            Vec3::new(-5.0 / 7.0 * radius, 0.0, 0.0),
            Vec3::new(0.0, -momentum, 0.0),
        ),
        body(
            Vec3::new(2.0 / 7.0 * radius, 0.0, 0.0),
            Vec3::new(0.0, momentum, 0.0),
        ),
    ];
    let half_period = core::f64::consts::PI * (semimajor.powi(3) / mu).sqrt();
    (
        MutualGravity::new(FRAME, 0.3, 0.1, 1).unwrap(),
        states,
        properties,
        half_period,
    )
}

#[test]
fn eccentric_binary_reaches_analytic_apocenter_with_refining_tolerances() {
    let (model, initial, properties, duration) = eccentric_binary();
    let mut errors = Vec::new();
    let mut accepted = Vec::new();
    for tolerance in [1e-4, 1e-6, 1e-8] {
        let result = model
            .advance(
                FRAME,
                &initial,
                &properties,
                duration,
                options(duration, tolerance),
                || false,
            )
            .unwrap();
        assert_eq!(result.duration_seconds, duration);
        assert!(result.stats.rejected_trials > 0);
        assert!(result.stats.maximum_accepted_error_ratio <= 1.0);
        assert_eq!(
            result.stats.trials,
            result.stats.accepted_trials + result.stats.rejected_trials
        );
        let relative = result.states_after[1]
            .pose()
            .position_world()
            .sub(result.states_after[0].pose().position_world());
        errors.push(norm(relative.sub(Vec3::new(-3.3, 0.0, 0.0))));
        accepted.push(result.stats.accepted_trials);
        let barycenter = result.states_after[0]
            .pose()
            .position_world()
            .scale(2.0 / 7.0)
            .add(
                result.states_after[1]
                    .pose()
                    .position_world()
                    .scale(5.0 / 7.0),
            );
        vector_close(barycenter, Vec3::ZERO, 2e-12);
        vector_close(
            result.diagnostics_after.linear_momentum_world,
            result.diagnostics_before.linear_momentum_world,
            2e-12,
        );
        vector_close(
            result.diagnostics_after.angular_momentum_world,
            result.diagnostics_before.angular_momentum_world,
            2e-11,
        );
        let direct_change = result.diagnostics_after.gravitational_potential_energy
            - result.diagnostics_before.gravitational_potential_energy;
        assert!((result.potential_energy_change_j - direct_change).abs() < 2e-11);
    }
    eprintln!("adaptive Kepler apocenter errors {errors:?}; accepted trials {accepted:?}");
    assert!(errors[1] < errors[0] / 5.0 && errors[2] < errors[1] / 5.0);
    assert!(errors[2] < 2e-4);
    assert!(accepted[0] < accepted[1] && accepted[1] < accepted[2]);
}

fn tidal_binary() -> ([RigidBodyState; 2], [MassProperties; 2]) {
    let principal_frames = [
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 3.0), 0.6).unwrap(),
        UnitQuaternion::from_axis_angle(Vec3::new(-2.0, 1.0, 0.7), -0.4).unwrap(),
    ];
    let properties = [
        MassProperties::from_principal_frame(
            2.0,
            Vec3::ZERO,
            Vec3::new(0.07, 0.09, 0.12),
            principal_frames[0],
        )
        .unwrap(),
        MassProperties::from_principal_frame(
            3.0,
            Vec3::ZERO,
            Vec3::new(0.10, 0.16, 0.20),
            principal_frames[1],
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

fn rotation_error(actual: &[RigidBodyState], reference: &[RigidBodyState]) -> f64 {
    actual
        .iter()
        .zip(reference)
        .flat_map(|(a, b)| {
            AXES.map(|axis| {
                norm(
                    a.pose()
                        .orientation()
                        .rotate_body_to_world(axis)
                        .sub(b.pose().orientation().rotate_body_to_world(axis)),
                )
            })
        })
        .fold(0.0, f64::max)
}

#[test]
fn tidal_rotation_controls_refinement_and_keeps_spin_orbit_exchange() {
    let (initial, properties) = tidal_binary();
    let radii = [0.6, 0.7];
    let model = MutualGravity::new(FRAME, 0.8, 0.01, 1)
        .unwrap()
        .quadrupole(&radii, 0.25)
        .unwrap();
    let duration = 0.8;
    let mut reference = initial.to_vec();
    for _ in 0..2048 {
        reference = model
            .step(FRAME, &reference, &properties, duration / 2048.0, || false)
            .unwrap()
            .states_after;
    }
    let mut orientation_errors = Vec::new();
    let mut spin_errors = Vec::new();
    let mut accepted = Vec::new();
    for tolerance in [1e-4, 1e-8] {
        let mut controls = options(duration, tolerance);
        // Translation cannot cause these refinements: only rotation and spin
        // accuracy are tightened, with the actual full body inertia retained.
        controls.absolute_position_tolerance_m = 1.0;
        controls.absolute_velocity_tolerance_m_s = 1.0;
        let result = model
            .advance(FRAME, &initial, &properties, duration, controls, || false)
            .unwrap();
        orientation_errors.push(rotation_error(&result.states_after, &reference));
        spin_errors.push(
            result
                .states_after
                .iter()
                .zip(&reference)
                .map(|(a, b)| norm(a.angular_momentum_body().sub(b.angular_momentum_body())))
                .fold(0.0, f64::max),
        );
        accepted.push(result.stats.accepted_trials);
        vector_close(
            result.diagnostics_after.linear_momentum_world,
            result.diagnostics_before.linear_momentum_world,
            3e-12,
        );
        vector_close(
            result.diagnostics_after.angular_momentum_world,
            result.diagnostics_before.angular_momentum_world,
            3e-11,
        );
        let spin_world = |states: &[RigidBodyState]| {
            states.iter().fold(Vec3::ZERO, |sum, s| {
                sum.add(
                    s.pose()
                        .orientation()
                        .rotate_body_to_world(s.angular_momentum_body()),
                )
            })
        };
        assert!(norm(spin_world(&result.states_after).sub(spin_world(&initial))) > 1e-6);
    }
    eprintln!(
        "adaptive tidal orientation errors {orientation_errors:?}; spin errors {spin_errors:?}"
    );
    assert!(accepted[1] > accepted[0]);
    assert!(orientation_errors[1] < orientation_errors[0] / 10.0);
    assert!(spin_errors[1] < spin_errors[0] / 10.0);
    assert!(orientation_errors[1] < 2e-5);
}

#[test]
fn cancellation_after_native_work_and_budget_exhaustion_publish_no_endpoint() {
    let (model, initial, properties, duration) = eccentric_binary();
    let snapshot = initial;
    for stop_at in [1, 128] {
        let mut polls = 0;
        let result = model.advance(
            FRAME,
            &initial,
            &properties,
            duration,
            options(duration, 1e-8),
            || {
                polls += 1;
                polls == stop_at
            },
        );
        assert!(matches!(
            result,
            Err(GravityAdvanceError::Gravity(GravityError::Cancelled))
        ));
        assert_eq!(polls, stop_at);
        assert_eq!(initial, snapshot);
    }
    let mut controls = options(duration, 1e-8);
    controls.maximum_trials = 1;
    assert!(matches!(
        model.advance(FRAME, &initial, &properties, duration, controls, || false),
        Err(GravityAdvanceError::TrialLimit)
    ));
    let mut controls = options(duration, 1e-14);
    controls.initial_step_seconds = duration;
    controls.minimum_step_seconds = duration;
    controls.maximum_step_seconds = duration;
    assert!(matches!(
        model.advance(FRAME, &initial, &properties, duration, controls, || false),
        Err(GravityAdvanceError::MinimumStepReached)
    ));
    assert_eq!(initial, snapshot);
}

#[test]
fn invalid_controls_frames_and_pair_domains_are_not_reported_as_convergence() {
    let (model, initial, properties, duration) = eccentric_binary();
    let mut controls = options(duration, 1e-6);
    controls.absolute_orientation_tolerance_rad = f64::NAN;
    assert!(matches!(
        model.advance(FRAME, &initial, &properties, duration, controls, || false),
        Err(GravityAdvanceError::InvalidOptions(_))
    ));
    assert!(matches!(
        model.advance(
            FRAME + 1,
            &initial,
            &properties,
            duration,
            options(duration, 1e-6),
            || false
        ),
        Err(GravityAdvanceError::Gravity(GravityError::FrameMismatch))
    ));
    assert!(matches!(
        model.advance(
            FRAME,
            &initial,
            &properties[..1],
            duration,
            options(duration, 1e-6),
            || false
        ),
        Err(GravityAdvanceError::Gravity(_))
    ));
    let insufficient = MutualGravity::new(FRAME, 0.3, 0.1, 0).unwrap();
    assert!(matches!(
        insufficient.advance(
            FRAME,
            &initial,
            &properties,
            duration,
            options(duration, 1e-6),
            || false
        ),
        Err(GravityAdvanceError::Gravity(_))
    ));
    let overlapping = [initial[0], initial[0]];
    assert!(matches!(
        model.advance(
            FRAME,
            &overlapping,
            &properties,
            duration,
            options(duration, 1e-6),
            || false
        ),
        Err(GravityAdvanceError::Gravity(_))
    ));
}

#[test]
fn zero_pair_systems_keep_exact_free_motion_and_deterministic_replay() {
    let model = MutualGravity::new(FRAME, 0.3, 0.0, 0).unwrap();
    let duration = 2.0;
    let empty = model
        .advance(FRAME, &[], &[], duration, options(duration, 1e-8), || false)
        .unwrap();
    assert!(empty.states_after.is_empty());
    assert_eq!(empty.potential_energy_change_j, 0.0);
    let properties = [MassProperties::new(2.0, Vec3::ZERO, Vec3::new(0.4, 0.4, 0.4)).unwrap()];
    let initial = [RigidBodyState::new(
        Pose::new(Vec3::new(1.0, -2.0, 0.5), UnitQuaternion::IDENTITY).unwrap(),
        Vec3::new(0.5, -0.25, 0.75),
        Vec3::new(0.0, 0.0, 0.2),
    )
    .unwrap()];
    let first = model
        .advance(
            FRAME,
            &initial,
            &properties,
            duration,
            options(duration, 1e-8),
            || false,
        )
        .unwrap();
    let replay = model
        .advance(
            FRAME,
            &initial,
            &properties,
            duration,
            options(duration, 1e-8),
            || false,
        )
        .unwrap();
    assert_eq!(first.states_after, replay.states_after);
    assert_eq!(first.suggested_step_seconds, replay.suggested_step_seconds);
    assert_eq!(first.stats.trials, replay.stats.trials);
    let actual = first.states_after[0];
    vector_close(
        actual.pose().position_world(),
        Vec3::new(1.5, -2.25, 1.25),
        2e-14,
    );
    vector_close(
        actual.linear_momentum_world(),
        initial[0].linear_momentum_world(),
        0.0,
    );
    let expected_orientation =
        UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), 1.0).unwrap();
    for axis in AXES {
        vector_close(
            actual.pose().orientation().rotate_body_to_world(axis),
            expected_orientation.rotate_body_to_world(axis),
            2e-13,
        );
    }
}
