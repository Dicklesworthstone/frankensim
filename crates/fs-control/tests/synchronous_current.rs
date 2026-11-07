use fs_control::ControlError;
use fs_control::synchronous::{
    CurrentControllerParameters, CurrentControllerState, SynchronousCurrentController,
    SynchronousMeasurement,
};

fn parameters() -> CurrentControllerParameters {
    CurrentControllerParameters {
        sample_period_s: 0.002,
        pole_pairs: 3,
        proportional_gain_ohm: [2.0, 3.0],
        integral_gain_ohm_per_s: [5.0, 7.0],
        anti_windup_time_s: [0.005, 0.007],
        winding_resistance_ohm: 0.4,
        current_limit_a: 10.0,
    }
}

fn measurement() -> SynchronousMeasurement {
    SynchronousMeasurement {
        current_alpha_beta_a: [0.0; 2],
        flux_alpha_beta_wb: [0.0; 2],
        relative_mechanical_angle_rad: 0.0,
        shaft_speed_rad_s: [0.0; 2],
        dc_bus_voltage_v: 100.0,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.17e}, expected {expected:.17e}, tolerance {tolerance:.3e}"
    );
}

fn close2(actual: [f64; 2], expected: [f64; 2], tolerance: f64) {
    for axis in 0..2 {
        close(actual[axis], expected[axis], tolerance);
    }
}

// Independent inverse Park oracle. This uses the platform trigonometric call,
// while production uses fs-math's deterministic range reduction and kernels.
fn stator_vector(dq: [f64; 2], phi: f64) -> [f64; 2] {
    let (s, c) = phi.sin_cos();
    [c * dq[0] - s * dq[1], s * dq[0] + c * dq[1]]
}

#[test]
fn actual_flux_relative_speed_and_sampled_pi_own_the_command() {
    let p = parameters();
    let controller = SynchronousCurrentController::new(p).unwrap();
    let state = CurrentControllerState {
        integral_voltage_dq_v: [0.2, -0.1],
    };
    let theta = 0.4;
    let phi = f64::from(p.pole_pairs) * theta;
    let mut sample = SynchronousMeasurement {
        current_alpha_beta_a: stator_vector([0.4, -0.9], phi),
        // These are deliberately supplied independently of nominal L*i+psi;
        // there is no inductance in the controller's runtime API.
        flux_alpha_beta_wb: stator_vector([0.12, -0.18], phi),
        relative_mechanical_angle_rad: theta,
        shaft_speed_rad_s: [3.0, -1.0],
        dc_bus_voltage_v: 100.0,
    };
    let proposal = controller
        .propose_controlled(&state, &sample, [0.7, -0.2], || false)
        .unwrap();
    close2(proposal.measured_current_dq_a, [0.4, -0.9], 3e-15);
    close2(proposal.actual_flux_dq_wb, [0.12, -0.18], 8e-16);
    close2(proposal.speed_voltage_dq_v, [2.16, 1.44], 8e-15);
    close2(proposal.unlimited_voltage_dq_v, [3.243, 3.3698], 2e-14);
    close2(
        proposal.commanded_voltage_dq_v,
        proposal.unlimited_voltage_dq_v,
        3e-14,
    );
    close2(
        proposal.next_state.integral_voltage_dq_v,
        [0.203, -0.0902],
        2e-14,
    );
    assert!(!proposal.current_limited && !proposal.voltage_limited);

    // Adding the same absolute speed to BOTH bodies must change nothing.
    sample.shaft_speed_rad_s = [10.0, 6.0];
    assert_eq!(
        proposal,
        controller
            .propose_controlled(&state, &sample, [0.7, -0.2], || false)
            .unwrap()
    );
    sample.shaft_speed_rad_s = [8.0, 8.0];
    let co_rotating = controller
        .propose_controlled(&state, &sample, [0.7, -0.2], || false)
        .unwrap();
    close2(co_rotating.speed_voltage_dq_v, [0.0, 0.0], 0.0);
}

#[test]
fn exact_hexagon_uses_physical_phase_span_and_preserves_signed_power() {
    let mut p = parameters();
    p.proportional_gain_ohm = [0.0; 2];
    p.integral_gain_ohm_per_s = [0.0; 2];
    p.winding_resistance_ohm = 0.0;
    let controller = SynchronousCurrentController::new(p).unwrap();
    let sample = SynchronousMeasurement {
        dc_bus_voltage_v: 10.0,
        ..measurement()
    };
    let physical_current = [0.9, -0.2, -0.7];
    let current_ab = [
        (2.0_f64 / 3.0).sqrt()
            * (physical_current[0] - (physical_current[1] + physical_current[2]) / 2.0),
        (physical_current[1] - physical_current[2]) / 2.0_f64.sqrt(),
    ];
    for raw in [
        [7.8, 0.0], // Above the inscribed-circle radius, inside the true hexagon.
        [12.0, 5.0],
        [-4.0, 17.0],
        [-8.0, -3.0],
        [0.0, 0.0],
    ] {
        let proposal = controller
            .propose_controlled(
                &CurrentControllerState {
                    integral_voltage_dq_v: raw,
                },
                &sample,
                [0.0; 2],
                || false,
            )
            .unwrap();
        // Direct balanced physical-phase oracle, independent of production
        // duty reconstruction and of any plant conversion helper.
        let phase = [
            (2.0_f64 / 3.0).sqrt() * raw[0],
            -raw[0] / 6.0_f64.sqrt() + raw[1] / 2.0_f64.sqrt(),
            -raw[0] / 6.0_f64.sqrt() - raw[1] / 2.0_f64.sqrt(),
        ];
        let low = phase.into_iter().fold(f64::INFINITY, f64::min);
        let high = phase.into_iter().fold(f64::NEG_INFINITY, f64::max);
        let span = high - low;
        let scale = if span > 10.0 { 10.0 / span } else { 1.0 };
        close2(
            proposal.commanded_voltage_alpha_beta_v,
            raw.map(|v| v * scale),
            5e-15,
        );
        assert_eq!(proposal.voltage_limited, span > 10.0);
        assert!(
            proposal
                .duty_high_fraction
                .iter()
                .all(|d| (0.0..=1.0).contains(d))
        );
        let dmin = proposal
            .duty_high_fraction
            .into_iter()
            .fold(f64::INFINITY, f64::min);
        let dmax = proposal
            .duty_high_fraction
            .into_iter()
            .fold(f64::NEG_INFINITY, f64::max);
        close(dmin + dmax, 1.0, 3e-16);
        if proposal.voltage_limited {
            close(dmin, 0.0, 0.0);
            close(dmax, 1.0, 0.0);
        }
        for sign in [-1.0, 1.0] {
            let bus_current: f64 = proposal
                .duty_high_fraction
                .iter()
                .zip(physical_current)
                .map(|(d, i)| d * sign * i)
                .sum();
            let electrical_power =
                proposal.commanded_voltage_alpha_beta_v[0] * sign * current_ab[0]
                    + proposal.commanded_voltage_alpha_beta_v[1] * sign * current_ab[1];
            close(10.0 * bus_current, electrical_power, 4e-15);
        }
    }
}

#[test]
fn current_command_limit_is_radial_and_does_not_hide_measured_overcurrent() {
    let mut p = parameters();
    p.current_limit_a = 2.0;
    let controller = SynchronousCurrentController::new(p).unwrap();
    let sample = SynchronousMeasurement {
        current_alpha_beta_a: [100.0, -100.0],
        ..measurement()
    };
    let proposal = controller
        .propose_controlled(
            &CurrentControllerState::default(),
            &sample,
            [3.0, 4.0],
            || false,
        )
        .unwrap();
    assert_eq!(proposal.requested_current_dq_a, [3.0, 4.0]);
    close2(proposal.current_reference_dq_a, [1.2, 1.6], 3e-16);
    close(
        proposal.current_reference_dq_a[0].hypot(proposal.current_reference_dq_a[1]),
        2.0,
        3e-16,
    );
    assert!(proposal.current_limited);
    assert_eq!(proposal.measured_current_dq_a, [100.0, -100.0]);
    assert!(proposal.voltage_limited);
}

#[test]
fn implicit_tracking_stays_bounded_for_a_large_sample_to_tracking_ratio() {
    let p = CurrentControllerParameters {
        sample_period_s: 0.1,
        proportional_gain_ohm: [2.0, 0.0],
        integral_gain_ohm_per_s: [5.0, 0.0],
        anti_windup_time_s: [0.001; 2],
        winding_resistance_ohm: 0.0,
        ..parameters()
    };
    let controller = SynchronousCurrentController::new(p).unwrap();
    let mut state = CurrentControllerState::default();
    let mut sample = SynchronousMeasurement {
        dc_bus_voltage_v: 1.0,
        ..measurement()
    };
    let voltage_limit = (2.0_f64 / 3.0).sqrt();
    let fixed_state = voltage_limit - 2.0 + p.anti_windup_time_s[0] * 5.0;
    let tracking_pole = p.anti_windup_time_s[0] / (p.anti_windup_time_s[0] + p.sample_period_s);
    for n in 1..=12 {
        let proposal = controller
            .propose_controlled(&state, &sample, [1.0, 0.0], || false)
            .unwrap();
        assert!(proposal.voltage_limited);
        close(proposal.commanded_voltage_dq_v[0], voltage_limit, 2e-16);
        state = proposal.next_state;
        // Closed form of the affine digital tracking recurrence. An explicit
        // Ts/Taw correction has pole -99 here and fails this bounded sequence.
        close(
            state.integral_voltage_dq_v[0],
            fixed_state * (1.0 - tracking_pole.powi(n)),
            2e-15,
        );
    }
    // Loss of bus voltage is a real actuator saturation. It must unwind the
    // same accepted state without a P/V division, invented voltage or NaNs.
    sample.dc_bus_voltage_v = 0.0;
    let before = state.integral_voltage_dq_v[0];
    let zero_bus_fixed_state = -2.0 + p.anti_windup_time_s[0] * 5.0;
    for n in 1..=8 {
        let proposal = controller
            .propose_controlled(&state, &sample, [1.0, 0.0], || false)
            .unwrap();
        assert_eq!(proposal.duty_high_fraction, [0.5; 3]);
        assert_eq!(proposal.commanded_voltage_alpha_beta_v, [0.0; 2]);
        assert!(proposal.voltage_limited);
        state = proposal.next_state;
        close(
            state.integral_voltage_dq_v[0],
            zero_bus_fixed_state + tracking_pole.powi(n) * (before - zero_bus_fixed_state),
            2e-15,
        );
    }
}

#[test]
fn invalid_cancelled_and_retried_samples_preserve_caller_state() {
    let p = parameters();
    for bad in [
        CurrentControllerParameters {
            sample_period_s: 0.0,
            ..p
        },
        CurrentControllerParameters { pole_pairs: 0, ..p },
        CurrentControllerParameters {
            integral_gain_ohm_per_s: [-1.0, 1.0],
            ..p
        },
        CurrentControllerParameters {
            current_limit_a: f64::INFINITY,
            ..p
        },
    ] {
        assert!(matches!(
            SynchronousCurrentController::new(bad),
            Err(ControlError::InvalidParameter { .. })
        ));
    }
    let controller = SynchronousCurrentController::new(p).unwrap();
    let state = CurrentControllerState {
        integral_voltage_dq_v: [0.3, -0.1],
    };
    let accepted = state;
    let sample = measurement();
    let expected = controller
        .propose_controlled(&state, &sample, [0.7, 0.2], || false)
        .unwrap();
    for cancel_at in [1, 2] {
        let mut polls = 0;
        assert_eq!(
            controller.propose_controlled(&state, &sample, [0.7, 0.2], || {
                polls += 1;
                polls == cancel_at
            }),
            Err(ControlError::Cancelled)
        );
        assert_eq!(polls, cancel_at);
        assert_eq!(state, accepted);
        assert_eq!(
            controller
                .propose_controlled(&state, &sample, [0.7, 0.2], || false)
                .unwrap(),
            expected
        );
    }
    for invalid in [
        SynchronousMeasurement {
            dc_bus_voltage_v: -0.01,
            ..sample
        },
        SynchronousMeasurement {
            flux_alpha_beta_wb: [f64::NAN, 0.0],
            ..sample
        },
        SynchronousMeasurement {
            shaft_speed_rad_s: [f64::INFINITY, 0.0],
            ..sample
        },
    ] {
        assert!(matches!(
            controller.propose_controlled(&state, &invalid, [0.7, 0.2], || false),
            Err(ControlError::InvalidSample { .. })
        ));
        assert_eq!(state, accepted);
    }
    let immediate = SynchronousCurrentController::new(CurrentControllerParameters {
        anti_windup_time_s: [0.0; 2],
        ..p
    })
    .unwrap();
    let zero_bus = SynchronousMeasurement {
        dc_bus_voltage_v: 0.0,
        ..sample
    };
    let proposal = immediate
        .propose_controlled(&state, &zero_bus, [0.7, 0.2], || false)
        .unwrap();
    close2(proposal.commanded_voltage_alpha_beta_v, [0.0; 2], 0.0);
    // Immediate tracking removes the entire predictor's unsupported voltage.
    close2(
        proposal.next_state.integral_voltage_dq_v,
        [
            -0.7 * (p.winding_resistance_ohm + p.proportional_gain_ohm[0]),
            -0.2 * (p.winding_resistance_ohm + p.proportional_gain_ohm[1]),
        ],
        4e-16,
    );
}
