//! G1/G2/G4: adaptive intervals retain the original field and rigid-body owners.
use std::cell::Cell;

use fs_mbd::gravity::{
    CentralGravity, GravityAdvanceError, GravityAdvanceOptions, GravityError, GravityLoad,
    RigidGravityField, SuperposedGravity, UniformGravity, advance_in_field,
};
use fs_mbd::{Gravity, MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3, Wrench};

const FRAME: u64 = 71;

fn options(tolerance: f64) -> GravityAdvanceOptions {
    GravityAdvanceOptions {
        initial_step_seconds: 0.5,
        minimum_step_seconds: 1e-10,
        maximum_step_seconds: 1.0,
        maximum_trials: 20_000,
        absolute_position_tolerance_m: tolerance,
        absolute_velocity_tolerance_m_s: tolerance,
        absolute_orientation_tolerance_rad: tolerance,
        absolute_angular_velocity_tolerance_rad_s: tolerance,
        relative_tolerance: 0.0,
    }
}

fn mass(value: f64) -> MassProperties {
    MassProperties::new(value, Vec3::ZERO, Vec3::new(0.1, 0.1, 0.1)).unwrap()
}

fn state(position: Vec3, momentum: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(position, UnitQuaternion::IDENTITY).unwrap(),
        momentum,
        Vec3::ZERO,
    )
    .unwrap()
}

fn distance(a: Vec3, b: Vec3) -> f64 {
    let d = a.sub(b);
    d.x.hypot(d.y).hypot(d.z)
}

#[test]
fn superposed_ballistics_retains_small_work_under_a_large_potential_datum() {
    let vertical = UniformGravity::new(
        FRAME,
        Gravity::new(Vec3::new(0.0, 0.0, -8.0)).unwrap(),
        Vec3::new(0.0, 0.0, -1e30),
    )
    .unwrap();
    let horizontal = UniformGravity::new(
        FRAME,
        Gravity::new(Vec3::new(2.0, 0.0, 0.0)).unwrap(),
        Vec3::ZERO,
    )
    .unwrap();
    let sources: [&dyn RigidGravityField; 2] = [&vertical, &horizontal];
    let field = SuperposedGravity::new(FRAME, &sources, 2).unwrap();
    let initial = state(Vec3::new(1.0, -2.0, 2.0), Vec3::new(3.0, 0.0, 0.0));
    let external = Wrench {
        force_world: Vec3::new(0.0, 12.0, 0.0),
        torque_body: Vec3::ZERO,
    };
    let result = advance_in_field(
        &field,
        FRAME,
        initial,
        mass(3.0),
        external,
        0.25,
        options(1e-11),
        || false,
    )
    .unwrap();
    let after = result.states_after[0];
    assert!(
        distance(
            after.pose().position_world(),
            Vec3::new(1.3125, -1.875, 1.75)
        ) < 2e-14
    );
    assert!(distance(after.linear_momentum_world(), Vec3::new(4.5, 3.0, -6.0)) < 2e-14);
    assert_eq!(
        result.diagnostics_before.gravitational_potential_energy,
        result.diagnostics_after.gravitational_potential_energy,
        "the absolute datum intentionally hides representable physical work"
    );
    assert!((result.potential_energy_change_j + 7.875).abs() < 2e-13);
    assert_eq!(result.stats.accepted_trials, 1);
    assert_eq!(result.stats.rejected_trials, 0);
    assert_eq!(result.stats.native_steps, 3);
}

#[test]
fn eccentric_fixed_source_orbit_resolves_pericenter_when_tolerance_tightens() {
    // a = 1, e = 0.7, mu = 1; apocenter to pericenter takes exactly pi seconds.
    let e: f64 = 0.7;
    let initial = state(
        Vec3::new(1.0 + e, 0.0, 0.0),
        Vec3::new(0.0, ((1.0 - e) / (1.0 + e)).sqrt(), 0.0),
    );
    let field = CentralGravity::new(FRAME, Vec3::ZERO, 1.0).unwrap();
    let solve = |tolerance| {
        advance_in_field(
            &field,
            FRAME,
            initial,
            mass(1.0),
            Wrench::ZERO,
            std::f64::consts::PI,
            options(tolerance),
            || false,
        )
        .unwrap()
    };
    let loose = solve(1e-4);
    let tight = solve(1e-7);
    let expected_position = Vec3::new(-(1.0 - e), 0.0, 0.0);
    let expected_velocity = Vec3::new(0.0, -((1.0 + e) / (1.0 - e)).sqrt(), 0.0);
    let loose_error = distance(
        loose.states_after[0].pose().position_world(),
        expected_position,
    );
    let tight_error = distance(
        tight.states_after[0].pose().position_world(),
        expected_position,
    );
    assert!(
        tight_error < loose_error / 20.0,
        "loose={loose_error:e}, tight={tight_error:e}"
    );
    assert!(
        tight_error < 5e-4,
        "pericenter position error {tight_error:e}"
    );
    assert!(
        distance(
            tight.states_after[0].linear_momentum_world(),
            expected_velocity
        ) < 3e-3
    );
    assert!(tight.stats.accepted_trials > loose.stats.accepted_trials);
    assert!(tight.stats.rejected_trials > 0);
    assert!(tight.stats.maximum_accepted_error_ratio <= 1.0);
    assert!(
        distance(
            tight.diagnostics_after.angular_momentum_world,
            tight.diagnostics_before.angular_momentum_world,
        ) < 1e-12
    );
    assert!(
        (tight.potential_energy_change_j
            - (tight.diagnostics_after.gravitational_potential_energy
                - tight.diagnostics_before.gravitational_potential_energy))
            .abs()
            < 2e-12
    );
    assert_eq!(
        tight,
        solve(1e-7),
        "identical inputs must replay identically"
    );
}

struct ObservedField {
    field: UniformGravity,
    calls: Cell<usize>,
    refuse_at: usize,
}

impl RigidGravityField for ObservedField {
    fn load(
        &self,
        frame: u64,
        pose: Pose,
        properties: MassProperties,
    ) -> Result<GravityLoad, GravityError> {
        let calls = self.calls.get() + 1;
        self.calls.set(calls);
        if calls == self.refuse_at {
            return Err(GravityError::Invalid("declared source refusal"));
        }
        self.field.load(frame, pose, properties)
    }
}

#[test]
fn late_field_refusal_cancellation_and_invalid_options_publish_no_interval() {
    let source = || ObservedField {
        field: UniformGravity::new(FRAME, Gravity::ZERO, Vec3::ZERO).unwrap(),
        calls: Cell::new(0),
        refuse_at: 10,
    };
    let initial = state(Vec3::ZERO, Vec3::new(1.0, 0.0, 0.0));
    let run = |field: &ObservedField, opts, stop: &mut dyn FnMut() -> bool| {
        advance_in_field(
            field,
            FRAME,
            initial,
            mass(1.0),
            Wrench::ZERO,
            2.0,
            opts,
            stop,
        )
    };

    let field = source();
    assert_eq!(
        run(&field, options(1e-8), &mut || false),
        Err(GravityAdvanceError::Gravity(GravityError::Invalid(
            "declared source refusal"
        ))),
    );
    assert_eq!(
        field.calls.get(),
        10,
        "an arbitrary source error must not be retried"
    );

    let field = ObservedField {
        refuse_at: usize::MAX,
        ..source()
    };
    assert_eq!(
        run(&field, options(1e-8), &mut || field.calls.get() >= 9),
        Err(GravityAdvanceError::Gravity(GravityError::Cancelled)),
    );

    let field = source();
    let mut invalid = options(1e-8);
    invalid.maximum_trials = 0;
    assert!(matches!(
        run(&field, invalid, &mut || false),
        Err(GravityAdvanceError::InvalidOptions(_)),
    ));
    assert_eq!(
        field.calls.get(),
        0,
        "options are admitted before source callbacks"
    );
}
