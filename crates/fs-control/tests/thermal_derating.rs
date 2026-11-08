use fs_control::ControlError;
use fs_control::speed::{
    ShaftPowerLimits, SpeedCommand, SpeedController, SpeedControllerParameters,
    SpeedControllerState,
};
use fs_control::thermal::{
    MAX_THERMAL_DERATING_CHANNELS, ThermalDeratingBand, ThermalDeratingPolicy, ThermalDeratingState,
};

const LIMITS: ShaftPowerLimits = ShaftPowerLimits {
    positive_torque_limit_nm: 12.0,
    negative_torque_limit_nm: 8.0,
    motoring_power_w: 300.0,
    regeneration_power_w: 100.0,
};
const BANDS: [ThermalDeratingBand; 2] = [
    ThermalDeratingBand {
        full_output_temperature_k: 350.0,
        trip_temperature_k: 400.0,
        reset_temperature_k: 330.0,
    },
    ThermalDeratingBand {
        full_output_temperature_k: 310.0,
        trip_temperature_k: 330.0,
        reset_temperature_k: 300.0,
    },
];

#[test]
fn each_component_derates_its_original_directional_capabilities_without_averaging() {
    let policy = ThermalDeratingPolicy::new(&BANDS).unwrap();
    let state = ThermalDeratingState::default();
    for (temperatures, fraction, channel) in [
        ([300.0, 300.0], 1.0, None),
        ([375.0, 300.0], 0.5, Some(0)),
        ([300.0, 325.0], 0.25, Some(1)),
        ([375.0, 320.0], 0.5, Some(0)),
    ] {
        let p = policy
            .propose_controlled(&state, &temperatures, LIMITS, false, || false)
            .unwrap();
        assert_eq!(p.output_fraction, fraction);
        assert_eq!(p.limiting_channel, channel);
        assert_eq!(p.limits.positive_torque_limit_nm, 12.0 * fraction);
        assert_eq!(p.limits.negative_torque_limit_nm, 8.0 * fraction);
        assert_eq!(p.limits.motoring_power_w, 300.0 * fraction);
        assert_eq!(p.limits.regeneration_power_w, 100.0 * fraction);
        assert_eq!(p.next_state, state);
    }
}

#[test]
fn overtemperature_latches_until_every_component_allows_an_explicit_reset() {
    let policy = ThermalDeratingPolicy::new(&BANDS).unwrap();
    let original = ThermalDeratingState::default();
    let trip = policy
        .propose_controlled(&original, &[400.0, 300.0], LIMITS, true, || false)
        .unwrap();
    assert!(trip.trip_started && trip.next_state.tripped && trip.reset_blocked);
    assert!(!trip.reset_accepted && !original.tripped);
    assert_eq!(trip.output_fraction, 0.0);
    let latched = trip.next_state;
    let cool_without_reset = policy
        .propose_controlled(&latched, &[300.0, 290.0], LIMITS, false, || false)
        .unwrap();
    assert!(cool_without_reset.next_state.tripped);
    assert_eq!(cool_without_reset.limiting_channel, None);
    let blocked = policy
        .propose_controlled(&latched, &[330.0, 301.0], LIMITS, true, || false)
        .unwrap();
    assert!(blocked.reset_blocked && blocked.next_state.tripped);
    let reset = policy
        .propose_controlled(&latched, &[330.0, 300.0], LIMITS, true, || false)
        .unwrap();
    assert!(reset.reset_accepted && !reset.next_state.tripped);
    assert_eq!(reset.limits, LIMITS);
}

#[test]
fn derating_changes_real_speed_demand_in_both_power_directions_and_tracks_the_stop() {
    let policy = ThermalDeratingPolicy::new(&BANDS).unwrap();
    let state = ThermalDeratingState::default();
    let speed = SpeedController::new(SpeedControllerParameters {
        sample_period_s: 0.01,
        proportional_gain_nm_s_per_rad: 1.0,
        integral_gain_nm_per_rad: 0.0,
        anti_windup_time_s: 0.0,
    })
    .unwrap();
    let warm = policy
        .propose_controlled(&state, &[375.0, 300.0], LIMITS, false, || false)
        .unwrap();
    for (shafts, reference, expected_torque, expected_power) in [
        ([60.0, 10.0], 150.0, 3.0, 150.0),
        ([60.0, 10.0], -50.0, -1.0, -50.0),
        ([-60.0, -10.0], -150.0, -3.0, 150.0),
        ([-60.0, -10.0], 50.0, 1.0, -50.0),
    ] {
        let d = speed
            .demand_controlled(
                &SpeedControllerState::default(),
                shafts,
                SpeedCommand {
                    relative_speed_rad_s: reference,
                    feedforward_torque_nm: 0.0,
                },
                warm.limits,
                || false,
            )
            .unwrap();
        assert_eq!(d.requested_torque_nm, expected_torque);
        assert_eq!(d.requested_shaft_power_w, expected_power);
    }
    let stopped = policy
        .propose_controlled(&state, &[400.0, 300.0], LIMITS, false, || false)
        .unwrap();
    let d = speed
        .demand_controlled(
            &SpeedControllerState::default(),
            [60.0, 10.0],
            SpeedCommand {
                relative_speed_rad_s: 150.0,
                feedforward_torque_nm: 0.0,
            },
            stopped.limits,
            || false,
        )
        .unwrap();
    assert_eq!(d.requested_torque_nm, 0.0);
    let p = d.finish_controlled(0.0, || false).unwrap();
    assert_eq!(p.next_state.integral_torque_nm, -100.0);
    assert_eq!(p.achievable_shaft_power_w, 0.0);
}

#[test]
fn cancellation_and_bad_sensor_samples_leave_the_accepted_latch_retryable() {
    let policy = ThermalDeratingPolicy::new(&BANDS).unwrap();
    let state = ThermalDeratingState::default();
    for temperatures in [&[300.0][..], &[300.0, f64::NAN], &[300.0, 0.0]] {
        assert!(
            policy
                .propose_controlled(&state, temperatures, LIMITS, false, || false)
                .is_err()
        );
        assert_eq!(state, ThermalDeratingState::default());
    }
    let mut polls = 0;
    assert_eq!(
        policy.propose_controlled(&state, &[400.0, 320.0], LIMITS, false, || {
            polls += 1;
            polls == 3
        }),
        Err(ControlError::Cancelled)
    );
    assert!(!state.tripped);
    let first = policy
        .propose_controlled(&state, &[400.0, 320.0], LIMITS, false, || false)
        .unwrap();
    let retry = policy
        .propose_controlled(&state, &[400.0, 320.0], LIMITS, false, || false)
        .unwrap();
    assert_eq!(first, retry);
}

#[test]
fn channel_budget_and_kelvin_threshold_order_are_admitted_before_sampling() {
    assert!(ThermalDeratingPolicy::new(&[]).is_err());
    let maximum = [BANDS[0]; MAX_THERMAL_DERATING_CHANNELS];
    assert!(ThermalDeratingPolicy::new(&maximum).is_ok());
    let too_many = [BANDS[0]; MAX_THERMAL_DERATING_CHANNELS + 1];
    assert!(ThermalDeratingPolicy::new(&too_many).is_err());
    for invalid in [
        ThermalDeratingBand {
            reset_temperature_k: 0.0,
            ..BANDS[0]
        },
        ThermalDeratingBand {
            reset_temperature_k: 351.0,
            ..BANDS[0]
        },
        ThermalDeratingBand {
            trip_temperature_k: 350.0,
            ..BANDS[0]
        },
        ThermalDeratingBand {
            trip_temperature_k: f64::INFINITY,
            ..BANDS[0]
        },
    ] {
        assert!(ThermalDeratingPolicy::new(&[invalid]).is_err());
    }
}
