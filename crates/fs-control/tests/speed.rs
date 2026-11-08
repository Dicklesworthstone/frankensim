//! Pure sampled policy tests; no plant or fs-phs dependency is imported here.
use fs_control::ControlError;
use fs_control::speed::{
    ShaftPowerLimits, SpeedCommand, SpeedController, SpeedControllerParameters,
    SpeedControllerState,
};

fn parameters() -> SpeedControllerParameters {
    SpeedControllerParameters {
        sample_period_s: 0.1,
        proportional_gain_nm_s_per_rad: 2.0,
        integral_gain_nm_per_rad: 3.0,
        anti_windup_time_s: 0.3,
    }
}
fn limits() -> ShaftPowerLimits {
    ShaftPowerLimits {
        positive_torque_limit_nm: 12.0,
        negative_torque_limit_nm: 8.0,
        motoring_power_w: 60.0,
        regeneration_power_w: 15.0,
    }
}
fn command(speed: f64, torque: f64) -> SpeedCommand {
    SpeedCommand {
        relative_speed_rad_s: speed,
        feedforward_torque_nm: torque,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn torque_and_power_limits_follow_all_four_relative_speed_power_quadrants() {
    let controller = SpeedController::new(SpeedControllerParameters {
        proportional_gain_nm_s_per_rad: 0.0,
        integral_gain_nm_per_rad: 0.0,
        ..parameters()
    })
    .unwrap();
    for (speed, torque, expected, bounds, power_limited) in [
        (5.0, 100.0, 12.0, [-3.0, 12.0], false),
        (5.0, -100.0, -3.0, [-3.0, 12.0], true),
        (-5.0, 100.0, 3.0, [-8.0, 3.0], true),
        (-5.0, -100.0, -8.0, [-8.0, 3.0], false),
    ] {
        let demand = controller
            .demand_controlled(
                &SpeedControllerState::default(),
                [speed, 0.0],
                command(speed, torque),
                limits(),
                || false,
            )
            .unwrap();
        assert_eq!(demand.torque_bounds_nm, bounds);
        assert_eq!(demand.requested_torque_nm, expected);
        assert_eq!(demand.requested_shaft_power_w, speed * expected);
        assert!(demand.torque_limited);
        assert_eq!(demand.power_limited, power_limited);
        assert!(demand.requested_shaft_power_w <= limits().motoring_power_w);
        assert!(demand.requested_shaft_power_w >= -limits().regeneration_power_w);
        let shifted = controller
            .demand_controlled(
                &SpeedControllerState::default(),
                [speed + 17.0, 17.0],
                command(speed, torque),
                limits(),
                || false,
            )
            .unwrap();
        assert_eq!(
            shifted, demand,
            "a common rotation cannot change reciprocal shaft power"
        );
        let proposal = demand.finish_controlled(expected, || false).unwrap();
        assert_eq!(proposal.achievable_shaft_power_w, speed * expected);
    }
}

#[test]
fn actual_allocator_torque_owns_anti_windup_and_the_exact_sampled_pi_increment() {
    let controller = SpeedController::new(parameters()).unwrap();
    let state = SpeedControllerState {
        integral_torque_nm: 0.4,
    };
    let capability = ShaftPowerLimits {
        motoring_power_w: 6.0,
        ..limits()
    };
    let demand = controller
        .demand_controlled(&state, [4.0, 1.0], command(5.0, 0.5), capability, || false)
        .unwrap();
    close(demand.speed_error_rad_s, 2.0, 0.0);
    close(demand.unlimited_torque_nm, 5.5, 1e-15);
    close(demand.requested_torque_nm, 2.0, 0.0);
    assert!(demand.power_limited && !demand.torque_limited);
    let proposal = demand.finish_controlled(1.5, || false).unwrap();
    close(proposal.next_state.integral_torque_nm, 0.0, 3e-16);
    assert!(proposal.allocation_limited);
    close(proposal.achievable_shaft_power_w, 4.5, 0.0);
    assert_eq!(state.integral_torque_nm, 0.4);
    let unconstrained = controller
        .demand_controlled(&state, [4.0, 1.0], command(5.0, 0.5), limits(), || false)
        .unwrap();
    let accepted = unconstrained
        .finish_controlled(unconstrained.requested_torque_nm, || false)
        .unwrap();
    close(accepted.next_state.integral_torque_nm, 1.0, 3e-16);
}

#[test]
fn saturated_tracking_matches_its_analytic_recurrence_and_zero_time_is_immediate() {
    let controller = SpeedController::new(SpeedControllerParameters {
        proportional_gain_nm_s_per_rad: 4.0,
        integral_gain_nm_per_rad: 2.0,
        anti_windup_time_s: 0.2,
        ..parameters()
    })
    .unwrap();
    let mut state = SpeedControllerState::default();
    let gamma = 1.0 / 3.0;
    let equilibrium = -2.6;
    for n in 1..=60 {
        let demand = controller
            .demand_controlled(&state, [1.0, 0.0], command(2.0, 0.0), limits(), || false)
            .unwrap();
        state = demand.finish_controlled(1.0, || false).unwrap().next_state;
        close(
            state.integral_torque_nm,
            equilibrium * (1.0 - (1.0_f64 - gamma).powi(n)),
            4e-15,
        );
    }
    let immediate = SpeedController::new(SpeedControllerParameters {
        anti_windup_time_s: 0.0,
        ..parameters()
    })
    .unwrap();
    let demand = immediate
        .demand_controlled(&state, [1.0, 0.0], command(4.0, 0.5), limits(), || false)
        .unwrap();
    close(
        demand
            .finish_controlled(1.0, || false)
            .unwrap()
            .next_state
            .integral_torque_nm,
        1.0 - 0.5 - 2.0 * 3.0,
        2e-15,
    );
}

#[test]
fn rest_uses_finite_torque_and_zero_capabilities_never_invent_power_or_isolation() {
    let controller = SpeedController::new(parameters()).unwrap();
    let state = SpeedControllerState::default();
    let zero_power = ShaftPowerLimits {
        motoring_power_w: 0.0,
        regeneration_power_w: 0.0,
        ..limits()
    };
    let resting = controller
        .demand_controlled(&state, [4.0, 4.0], command(10.0, 0.0), zero_power, || false)
        .unwrap();
    assert_eq!(resting.requested_torque_nm, 12.0);
    assert_eq!(resting.requested_shaft_power_w, 0.0);
    let rotating = controller
        .demand_controlled(&state, [4.0, 3.0], command(10.0, 0.0), zero_power, || false)
        .unwrap();
    assert_eq!(rotating.requested_torque_nm, 0.0);
    let disabled = ShaftPowerLimits {
        positive_torque_limit_nm: 0.0,
        negative_torque_limit_nm: 0.0,
        ..zero_power
    };
    let stopped = controller
        .demand_controlled(&state, [0.0; 2], command(10.0, 0.0), disabled, || false)
        .unwrap();
    assert_eq!(stopped.requested_torque_nm, 0.0);
    assert!(stopped.finish_controlled(1e-20, || false).is_err());
    // An unrelated large reverse capability must not excuse positive power
    // when the positive side is disabled.
    let one_way = ShaftPowerLimits {
        positive_torque_limit_nm: 0.0,
        negative_torque_limit_nm: 1e20,
        ..limits()
    };
    let demand = controller
        .demand_controlled(&state, [0.0; 2], command(0.0, 0.0), one_way, || false)
        .unwrap();
    assert!(demand.finish_controlled(1.0, || false).is_err());
}

#[test]
fn invalid_samples_allocator_escape_and_cancellation_leave_accepted_state_unchanged() {
    for bad in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(
            SpeedController::new(SpeedControllerParameters {
                sample_period_s: bad,
                ..parameters()
            })
            .is_err()
        );
    }
    let controller = SpeedController::new(parameters()).unwrap();
    let state = SpeedControllerState {
        integral_torque_nm: 0.25,
    };
    assert!(
        controller
            .demand_controlled(&state, [f64::NAN, 0.0], command(1.0, 0.0), limits(), || {
                false
            })
            .is_err()
    );
    assert!(
        controller
            .demand_controlled(
                &state,
                [0.0; 2],
                command(1.0, 0.0),
                ShaftPowerLimits {
                    regeneration_power_w: -1.0,
                    ..limits()
                },
                || false
            )
            .is_err()
    );
    let mut polls = 0;
    assert_eq!(
        controller.demand_controlled(&state, [0.0; 2], command(1.0, 0.0), limits(), || {
            polls += 1;
            polls == 2
        }),
        Err(ControlError::Cancelled)
    );
    let demand = controller
        .demand_controlled(&state, [0.0; 2], command(1.0, 0.0), limits(), || false)
        .unwrap();
    assert!(demand.finish_controlled(20.0, || false).is_err());
    assert!(demand.finish_controlled(f64::NAN, || false).is_err());
    let mut polls = 0;
    assert_eq!(
        demand.finish_controlled(1.0, || {
            polls += 1;
            polls == 2
        }),
        Err(ControlError::Cancelled)
    );
    let first = demand.finish_controlled(1.0, || false).unwrap();
    let replay = controller
        .demand_controlled(&state, [0.0; 2], command(1.0, 0.0), limits(), || false)
        .unwrap()
        .finish_controlled(1.0, || false)
        .unwrap();
    assert_eq!(first, replay);
    assert_eq!(state.integral_torque_nm, 0.25);
}
