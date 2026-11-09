//! G0 command timing and native device integration for a manufactured PWM drive.
use fs_math::dd::Dd;
use fs_phs::SolveBudget;
use fs_phs::circuit::converter::commutation::pwm::{CenterAlignedPwm, MAX_PWM_INTERVALS};
use fs_phs::circuit::converter::commutation::{
    BodyDiode, CommutationInterval, CommutationParameters, DeviceGateState, DeviceOperatingBounds,
    PreparedCommutatedSwitchCycle, SwitchDevice,
};
use fs_phs::circuit::converter::power_load::PowerLoadControls;
use fs_phs::circuit::converter::switching::SwitchCycleControls;

const GATES: [[bool; 4]; 9] = [
    [false, true, false, true],
    [false, false, false, true],
    [true, false, false, true],
    [true, false, false, false],
    [true, false, true, false],
    [true, false, false, false],
    [true, false, false, true],
    [false, false, false, true],
    [false, true, false, true],
];
const DURATIONS: [f64; 9] = [3.5, 1.0, 1.0, 1.0, 3.0, 1.0, 1.0, 1.0, 3.5];

#[test]
fn both_legs_retain_the_authored_center_and_each_full_dead_interval() {
    let request = CenterAlignedPwm {
        period_s: 16.0,
        high_fraction: [0.5, 0.25],
        dead_time_s: [1.0, 1.0],
    };
    let program = request.prepare().unwrap();
    assert_eq!(program.intervals().len(), MAX_PWM_INTERVALS);
    for (i, interval) in program.intervals().iter().enumerate() {
        assert_eq!(interval.duration_s, DURATIONS[i]);
        assert_eq!(interval.gates.on, GATES[i]);
    }
    assert_eq!(program.timing().channel_on_s, [7.0, 7.0, 3.0, 11.0]);
    assert_eq!(program.timing().both_off_s, [2.0, 2.0]);
    assert_eq!(program.command(), request);
    assert_eq!(program.boundary_gates().on, GATES[0]);
    assert_eq!(program, request.prepare().unwrap());
    program.admit_after(&program).unwrap();
}

#[test]
fn coincident_edges_and_zero_dead_time_apply_atomic_complementary_commands() {
    for dead in [0.0, 1.0] {
        let program = CenterAlignedPwm {
            period_s: 16.0,
            high_fraction: [0.5, 0.5],
            dead_time_s: [dead; 2],
        }
        .prepare()
        .unwrap();
        assert_eq!(program.intervals().len(), if dead == 0.0 { 3 } else { 5 });
        for interval in program.intervals() {
            assert!(interval.duration_s > 0.0);
            assert_eq!(&interval.gates.on[..2], &interval.gates.on[2..]);
            assert!(!(interval.gates.on[0] && interval.gates.on[1]));
        }
        assert_eq!(program.timing().channel_on_s, [8.0 - dead; 4]);
        assert_eq!(program.timing().both_off_s, [2.0 * dead; 2]);
        if dead == 0.0 {
            assert_eq!(
                program
                    .intervals()
                    .iter()
                    .map(|x| x.duration_s)
                    .collect::<Vec<_>>(),
                [4.0, 8.0, 4.0]
            );
        }
    }
}

#[test]
fn held_endpoints_have_no_synthetic_edges_and_boundary_gate_changes_refuse() {
    let request = CenterAlignedPwm {
        period_s: 0.001,
        high_fraction: [0.0, 1.0],
        dead_time_s: [1e-6, 2e-6],
    };
    let held = request.prepare().unwrap();
    assert_eq!(held.intervals().len(), 1);
    assert_eq!(held.intervals()[0].gates.on, [false, true, true, false]);
    assert_eq!(
        held.timing().channel_on_s,
        [0.0, request.period_s, request.period_s, 0.0]
    );
    assert_eq!(held.timing().both_off_s, [0.0; 2]);
    let changing = CenterAlignedPwm {
        high_fraction: [0.4, 1.0],
        ..request
    }
    .prepare()
    .unwrap();
    changing.admit_after(&held).unwrap();
    held.admit_after(&changing).unwrap();
    let opposite = CenterAlignedPwm {
        high_fraction: [1.0, 1.0],
        ..request
    }
    .prepare()
    .unwrap();
    assert!(held.admit_after(&opposite).is_err());
    assert!(opposite.admit_after(&held).is_err());
}

#[test]
fn pulse_erasure_invalid_data_and_unrepresentable_blank_intervals_are_not_clipped() {
    let valid = CenterAlignedPwm {
        period_s: 1.0,
        high_fraction: [0.5; 2],
        dead_time_s: [0.1; 2],
    };
    for period_s in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(CenterAlignedPwm { period_s, ..valid }.prepare().is_err());
    }
    for duty in [-0.01, 1.01, f64::NAN, f64::INFINITY] {
        assert!(
            CenterAlignedPwm {
                high_fraction: [duty, 0.5],
                ..valid
            }
            .prepare()
            .is_err()
        );
    }
    for dead in [-0.01, 0.5, 0.6, f64::NAN, f64::INFINITY] {
        assert!(
            CenterAlignedPwm {
                dead_time_s: [dead, 0.1],
                ..valid
            }
            .prepare()
            .is_err()
        );
    }
    assert!(
        CenterAlignedPwm {
            high_fraction: [0.1, 0.5],
            ..valid
        }
        .prepare()
        .is_err()
    );
    assert!(
        CenterAlignedPwm {
            dead_time_s: [f64::EPSILON / 32.0, 0.0],
            ..valid
        }
        .prepare()
        .is_err()
    );
    assert!(
        CenterAlignedPwm {
            high_fraction: [f64::from_bits(1), 0.5],
            dead_time_s: [0.0; 2],
            ..valid
        }
        .prepare()
        .is_err()
    );
    valid.prepare().unwrap();
}

#[test]
fn triangular_carrier_oracle_and_time_rescaling_retain_physical_waveform() {
    let reference = CenterAlignedPwm {
        period_s: 0.0013,
        high_fraction: [0.31, 0.67],
        dead_time_s: [0.000013, 0.000021],
    }
    .prepare()
    .unwrap();
    for scale in [1e-9, 1.0, 1e9] {
        let request = CenterAlignedPwm {
            period_s: reference.period_s() * scale,
            high_fraction: reference.command().high_fraction,
            dead_time_s: reference.command().dead_time_s.map(|x| x * scale),
        };
        let program = request.prepare().unwrap();
        assert_eq!(program.intervals().len(), reference.intervals().len());
        let mut elapsed = Dd::ZERO;
        for (interval, original) in program.intervals().iter().zip(reference.intervals()) {
            assert_eq!(interval.gates, original.gates);
            assert!((interval.duration_s / scale - original.duration_s).abs() < 2e-18);
            for fraction in [0.25, 0.5, 0.75] {
                let t = elapsed.to_f64() + fraction * interval.duration_s;
                // Compare the physical waveform to an independent triangular
                // carrier and two separated thresholds, not to sorted edges.
                let carrier = (2.0 * t / request.period_s - 1.0).abs();
                for leg in 0..2 {
                    let blank_fraction = request.dead_time_s[leg] / request.period_s;
                    let expected_high = carrier < request.high_fraction[leg] - blank_fraction;
                    let expected_low = carrier > request.high_fraction[leg] + blank_fraction;
                    assert_eq!(interval.gates.on[2 * leg], expected_high);
                    assert_eq!(interval.gates.on[2 * leg + 1], expected_low);
                }
            }
            elapsed = elapsed + Dd::from_f64(interval.duration_s);
        }
        assert_eq!(elapsed.to_f64(), request.period_s);
    }
}

#[test]
fn generated_program_drives_the_native_diode_model_like_the_physical_gate_deck() {
    let period_s = 1.0 / 2048.0;
    let program = CenterAlignedPwm {
        period_s,
        high_fraction: [0.5, 0.25],
        dead_time_s: [period_s / 16.0; 2],
    }
    .prepare()
    .unwrap();
    let authored: [CommutationInterval; 9] = core::array::from_fn(|i| CommutationInterval {
        gates: DeviceGateState { on: GATES[i] },
        duration_s: DURATIONS[i] * period_s / 16.0,
    });
    let device = SwitchDevice {
        on_resistance_ohm: 0.2,
        off_resistance_ohm: 1000.0,
        body_diode: BodyDiode {
            forward_voltage_v: 0.7,
            differential_resistance_ohm: 0.1,
        },
        bounds: DeviceOperatingBounds {
            maximum_reverse_voltage_v: 100.0,
            maximum_channel_current_a: 100.0,
            maximum_diode_current_a: 100.0,
        },
    };
    let parameters = CommutationParameters {
        inductance_h: 0.02,
        bus_capacitance_f: [0.2, 0.3],
        winding_resistance_ohm: 0.05,
        devices: [device; 4],
    };
    let electrical = PowerLoadControls {
        solve: SolveBudget {
            relative: 1e-13,
            absolute: 1e-12,
            max_iterations: 40,
        },
        backtracking_halvings: 8,
        absolute_energy_tolerance_j: 1e-10,
        relative_energy_tolerance: 1e-12,
    };
    let controls = SwitchCycleControls {
        max_intervals: 64,
        absolute_energy_tolerance_j: 3e-9,
        relative_energy_tolerance: 1e-11,
    };
    let initial = parameters.state(0.7, [5.0, 2.0]).unwrap();
    let run = |intervals: &[CommutationInterval]| {
        let mut native = PreparedCommutatedSwitchCycle::new_controlled(
            parameters,
            [10.0; 2],
            intervals,
            4,
            electrical,
            controls,
            || false,
        )
        .unwrap();
        assert_eq!(native.period_s(), period_s);
        let mut next = [0.0; 3];
        let mut mean = [0.0; 2];
        let report = native
            .step_controlled(
                &initial,
                &[5.0, 2.0],
                &[None; 2],
                &mut next,
                &mut mean,
                || false,
            )
            .unwrap();
        (next, mean, report)
    };
    let (actual, mean, report) = run(program.intervals());
    let (expected, expected_mean, expected_report) = run(&authored);
    assert_eq!(actual, expected);
    assert_eq!(mean, expected_mean);
    assert_eq!(report.channel_heat_j, expected_report.channel_heat_j);
    assert_eq!(report.diode_heat_j, expected_report.diode_heat_j);
    assert!(report.diode_heat_j.iter().sum::<f64>() > 0.0);
    assert_eq!(report.integration_steps, 36);
}
