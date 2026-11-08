//! Physical domain refusals caused by an unresolved split drift can be refined.
use fs_mbd::gravity::mutual::MutualGravity;
use fs_mbd::gravity::{
    GravityAdvanceError, GravityAdvanceOptions, GravityAdvanceReceipt, GravityError,
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

const FRAME: u64 = 101;
const DURATION: f64 = 4.0;

fn norm(value: Vec3) -> f64 {
    value.x.hypot(value.y).hypot(value.z)
}

fn circular_binary() -> ([RigidBodyState; 2], [MassProperties; 2]) {
    let body = |x, py| {
        RigidBodyState::new(
            Pose::new(Vec3::new(x, 0.0, 0.0), UnitQuaternion::IDENTITY).unwrap(),
            Vec3::new(0.0, py, 0.0),
            Vec3::ZERO,
        )
        .unwrap()
    };
    (
        [body(-0.5, -0.5), body(0.5, 0.5)],
        [MassProperties::new(1.0, Vec3::ZERO, Vec3::new(0.004, 0.004, 0.004)).unwrap(); 2],
    )
}

fn options() -> GravityAdvanceOptions {
    GravityAdvanceOptions {
        initial_step_seconds: DURATION,
        minimum_step_seconds: 1e-8,
        maximum_step_seconds: DURATION,
        maximum_trials: 10_000,
        absolute_position_tolerance_m: 1e-8,
        absolute_velocity_tolerance_m_s: 1e-8,
        absolute_orientation_tolerance_rad: 1e-8,
        absolute_angular_velocity_tolerance_rad_s: 1e-8,
        relative_tolerance: 0.0,
    }
}

fn assert_resolved_circle(result: &GravityAdvanceReceipt) {
    // G(m1+m2)=1 and the initial relative radius and speed are both 1.
    // The exact relative orbit is (cos t, sin t), entirely outside r=0.5.
    let exact_position = Vec3::new(DURATION.cos(), DURATION.sin(), 0.0);
    let exact_velocity = Vec3::new(-DURATION.sin(), DURATION.cos(), 0.0);
    let relative_position = result.states_after[1]
        .pose()
        .position_world()
        .sub(result.states_after[0].pose().position_world());
    let relative_velocity = result.states_after[1]
        .linear_momentum_world()
        .sub(result.states_after[0].linear_momentum_world());
    assert_eq!(result.duration_seconds, DURATION);
    assert!(result.stats.rejected_trials > 0);
    assert!(result.stats.maximum_accepted_error_ratio <= 1.0);
    assert!(norm(relative_position.sub(exact_position)) < 2e-4);
    assert!(norm(relative_velocity.sub(exact_velocity)) < 2e-4);
    assert!(norm(result.diagnostics_after.linear_momentum_world) < 1e-13);
    assert!(
        norm(
            result
                .diagnostics_after
                .angular_momentum_world
                .sub(Vec3::new(0.0, 0.0, 0.5))
        ) < 1e-11
    );
    assert!((result.diagnostics_after.mechanical_energy + 0.25).abs() < 1e-6);
}

#[test]
fn false_coarse_cutoff_crossing_refines_to_the_physical_circular_orbit() {
    let (initial, properties) = circular_binary();
    let model = MutualGravity::new(FRAME, 0.5, 0.5, 1).unwrap();
    // One native drift joins relative positions (1,0) and (-7,4), passing
    // at 1/sqrt(5) < 0.5. The true circular orbit never enters this cutoff.
    assert!(matches!(
        model.step(FRAME, &initial, &properties, DURATION, || false),
        Err(GravityError::SeparationCutoff { .. })
    ));
    let result = model
        .advance(FRAME, &initial, &properties, DURATION, options(), || false)
        .unwrap();
    assert_resolved_circle(&result);
}

#[test]
fn false_coarse_quadrupole_domain_crossing_refines_without_changing_the_model() {
    let (initial, properties) = circular_binary();
    let radii = [0.1, 0.1];
    // Isotropic moments make the quadrupole correction vanish. The declared
    // combined-radius ratio 0.2/r <= 0.4 admits precisely r >= 0.5.
    let model = MutualGravity::new(FRAME, 0.5, 0.0, 1)
        .unwrap()
        .quadrupole(&radii, 0.4)
        .unwrap();
    assert!(matches!(
        model.step(FRAME, &initial, &properties, DURATION, || false),
        Err(GravityError::OutsideExpansionDomain)
    ));
    let result = model
        .advance(FRAME, &initial, &properties, DURATION, options(), || false)
        .unwrap();
    assert_resolved_circle(&result);
}

#[test]
fn initial_domain_errors_and_pair_budget_errors_are_returned_without_refinement() {
    let (initial, properties) = circular_binary();
    let mut controls = options();
    controls.maximum_trials = 1;
    let excluded = MutualGravity::new(FRAME, 0.5, 1.0, 1).unwrap();
    assert!(matches!(
        excluded.advance(FRAME, &initial, &properties, DURATION, controls, || false),
        Err(GravityAdvanceError::Gravity(
            GravityError::SeparationCutoff { .. }
        ))
    ));
    let radii = [0.1, 0.1];
    let excluded_tidal = MutualGravity::new(FRAME, 0.5, 0.0, 1)
        .unwrap()
        .quadrupole(&radii, 0.1)
        .unwrap();
    assert!(matches!(
        excluded_tidal.advance(FRAME, &initial, &properties, DURATION, controls, || false),
        Err(GravityAdvanceError::Gravity(
            GravityError::OutsideExpansionDomain
        ))
    ));
    let insufficient = MutualGravity::new(FRAME, 0.5, 0.0, 0).unwrap();
    assert!(matches!(
        insufficient.advance(FRAME, &initial, &properties, DURATION, controls, || false),
        Err(GravityAdvanceError::Gravity(GravityError::Invalid(
            "mutual gravity pair budget exceeded"
        )))
    ));
}

#[test]
fn a_real_radial_collision_cannot_publish_a_completed_interval() {
    let (circular, properties) = circular_binary();
    let initial =
        circular.map(|state| RigidBodyState::new(state.pose(), Vec3::ZERO, Vec3::ZERO).unwrap());
    let snapshot = initial;
    let model = MutualGravity::new(FRAME, 0.5, 0.5, 1).unwrap();
    assert!(model.load(FRAME, &initial, &properties, || false).is_ok());
    // From rest at r=1 with mu=1, radial free fall reaches r=0.5 at
    // (pi/4 + 1/2)/sqrt(2) seconds and the singularity at pi/(2sqrt(2)).
    // No valid conservative trajectory covers the requested two seconds.
    let result = model.advance(FRAME, &initial, &properties, 2.0, options(), || false);
    assert!(
        result.is_err(),
        "a true domain crossing escaped as {result:?}"
    );
    assert_eq!(initial, snapshot);
}
