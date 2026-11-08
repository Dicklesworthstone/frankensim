//! Physical bridge states, reciprocal DC/phase power, losses and warm duty changes.
use fs_phs::circuit::synchronous::{
    SynchronousMotor, SynchronousMotorParameters,
    inverter::{AveragedSynchronousInverter, InverterParameters, PreparedSynchronousInverter},
};
use fs_phs::{PreparedStepError, SolveBudget};

#[path = "circuit_synchronous_inverter/battery.rs"]
mod battery;
#[path = "circuit_synchronous_inverter/trajectory.rs"]
mod trajectory;

fn machine() -> SynchronousMotorParameters {
    SynchronousMotorParameters {
        d_inductance_h: 0.012,
        q_inductance_h: 0.019,
        winding_resistance_ohm: 0.35,
        pm_flux_linkage_wb: 0.08,
        pole_pairs: 3,
        rotor_inertia_kg_m2: 0.025,
        stator_inertia_kg_m2: 0.07,
        bearing_drag_nm_s_per_rad: 0.002,
    }
}

fn ideal() -> InverterParameters {
    InverterParameters {
        dc_capacitance_f: 0.004,
        high_on_resistance_ohm: [0.0; 3],
        low_on_resistance_ohm: [0.0; 3],
    }
}

fn budget() -> SolveBudget {
    SolveBudget {
        relative: 2e-11,
        absolute: 2e-12,
        max_iterations: 30,
    }
}

fn prepared(
    p: SynchronousMotorParameters,
    bridge: InverterParameters,
    current: [f64; 2],
    voltage: f64,
) -> (PreparedSynchronousInverter, [f64; 6]) {
    let motor = SynchronousMotor::new(p).unwrap();
    let state = motor.state(current, 0.31, [2.0, -0.5]).unwrap();
    let inverter = AveragedSynchronousInverter::new(motor, bridge).unwrap();
    let state = inverter.state(voltage, &state).unwrap();
    (
        inverter.prepare_controlled(budget(), 8, || false).unwrap(),
        state,
    )
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; tolerance={tolerance:e}"
    );
}

fn heat(channels: &[[f64; 2]; 3]) -> f64 {
    channels.iter().flatten().sum()
}

// Snapshot the published control and all original operators/source factors,
// rather than only the caller's endpoint buffers.
fn stamp(
    model: &PreparedSynchronousInverter,
) -> ([f64; 3], Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let (j, r, g) = model.system().structure();
    (
        std::array::from_fn(|i| model.high_fraction()[i]),
        j.to_vec(),
        r.to_vec(),
        g.to_vec(),
        model.system().resistance_witness().unwrap().rows().to_vec(),
    )
}

#[test]
fn eight_physical_switch_vectors_match_dc_current_and_phase_power() {
    let p = SynchronousMotorParameters {
        winding_resistance_ohm: 0.0,
        bearing_drag_nm_s_per_rad: 0.0,
        ..machine()
    };
    let (mut model, initial) = prepared(p, ideal(), [-0.7, 2.3], 24.0);
    let dt = 1e-4;
    for mask in 0..8 {
        let duty = std::array::from_fn(|i| f64::from((mask >> i) & 1));
        let (mut next, mut ports) = ([0.0; 6], [0.0; 3]);
        let step = model
            .step_controlled(&initial, duty, &[0.0; 3], dt, &mut next, &mut ports, || {
                false
            })
            .unwrap();
        let [ia, ib, ic] = step.phase_current_a;
        // Physical upper-switch current identities, using a,b,c as bits 0,1,2.
        let expected_current = match mask {
            0 | 7 => 0.0,
            1 => ia,
            2 => ib,
            3 => -ic,
            4 => ic,
            5 => -ib,
            6 => -ia,
            _ => unreachable!(),
        };
        close(ia + ib + ic, 0.0, 2e-14);
        close(step.dc_bridge_current_a, expected_current, 2e-13);
        close(next[0] - initial[0], -dt * expected_current, 3e-12);
        close(
            ports[0],
            0.5 * (initial[0] + next[0]) / ideal().dc_capacitance_f,
            2e-12,
        );
        let neutral = ports[0] * duty.iter().sum::<f64>() / 3.0;
        for phase in 0..3 {
            close(
                step.phase_voltage_v[phase],
                ports[0] * duty[phase] - neutral,
                2e-12,
            );
        }
        let phase_power: f64 = step
            .phase_voltage_v
            .iter()
            .zip(step.phase_current_a)
            .map(|(v, i)| v * i)
            .sum();
        close(phase_power, ports[0] * expected_current, 2e-11);
        close(step.bridge_work_j, dt * ports[0] * expected_current, 2e-13);
        close(step.machine_electrical_work_j, dt * phase_power, 2e-13);
        close(step.native.balance_residual(), 0.0, 2e-10);
        assert_eq!(heat(&step.channel_heat_j), 0.0);
    }
}

#[test]
fn common_mode_preserves_ideal_or_matched_channel_dynamics_and_zero_vectors() {
    for matched in [[0.0; 3], [0.02, 0.05, 0.09]] {
        let bridge = InverterParameters {
            high_on_resistance_ohm: matched,
            low_on_resistance_ohm: matched,
            ..ideal()
        };
        let (mut model, initial) = prepared(machine(), bridge, [0.9, -1.6], 24.0);
        for (first, second) in [([0.1, 0.4, 0.6], [0.3, 0.6, 0.8]), ([0.0; 3], [1.0; 3])] {
            let (mut a, mut b, mut ya, mut yb) = ([0.0; 6], [0.0; 6], [0.0; 3], [0.0; 3]);
            let sa = model
                .step_controlled(&initial, first, &[0.0; 3], 2e-4, &mut a, &mut ya, || false)
                .unwrap();
            let sb = model
                .step_controlled(&initial, second, &[0.0; 3], 2e-4, &mut b, &mut yb, || false)
                .unwrap();
            for i in 0..6 {
                close(a[i], b[i], 3e-12);
            }
            for i in 0..3 {
                close(ya[i], yb[i], 3e-12);
                close(sa.phase_voltage_v[i], sb.phase_voltage_v[i], 3e-12);
            }
            close(sa.dc_bridge_current_a, sb.dc_bridge_current_a, 3e-13);
            close(heat(&sa.channel_heat_j), heat(&sb.channel_heat_j), 2e-13);
            if first == [0.0; 3] {
                close(sa.dc_bridge_current_a, 0.0, 2e-14);
                close(a[0], initial[0], 2e-14);
                for phase in 0..3 {
                    assert_eq!(sa.channel_heat_j[phase][0], 0.0);
                    assert_eq!(sb.channel_heat_j[phase][1], 0.0);
                }
            }
        }
    }
}

#[test]
fn asymmetric_channels_have_their_own_heat_and_actual_floating_neutral_voltage() {
    let p = machine();
    let bridge = InverterParameters {
        high_on_resistance_ohm: [0.11, 0.08, 0.17],
        low_on_resistance_ohm: [0.01, 0.03, 0.02],
        ..ideal()
    };
    let (mut model, initial) = prepared(p, bridge, [0.9, -1.6], 24.0);
    let dt = 2e-4;
    let mut totals = Vec::new();
    for duty in [[0.1, 0.4, 0.6], [0.3, 0.6, 0.8]] {
        let inputs = [0.4, -0.03, 0.01];
        let (mut next, mut ports) = ([0.0; 6], [0.0; 3]);
        let step = model
            .step_controlled(&initial, duty, &inputs, dt, &mut next, &mut ports, || false)
            .unwrap();
        let terminal: [f64; 3] = std::array::from_fn(|i| {
            let r = duty[i] * bridge.high_on_resistance_ohm[i]
                + (1.0 - duty[i]) * bridge.low_on_resistance_ohm[i];
            ports[0] * duty[i] - r * step.phase_current_a[i]
        });
        let neutral = terminal.iter().sum::<f64>() / 3.0;
        for i in 0..3 {
            let square = step.phase_current_a[i] * step.phase_current_a[i];
            close(
                step.channel_heat_j[i][0],
                dt * duty[i] * bridge.high_on_resistance_ohm[i] * square,
                3e-14,
            );
            close(
                step.channel_heat_j[i][1],
                dt * (1.0 - duty[i]) * bridge.low_on_resistance_ohm[i] * square,
                3e-14,
            );
            close(step.phase_voltage_v[i], terminal[i] - neutral, 3e-12);
        }
        close(step.phase_voltage_v.iter().sum(), 0.0, 3e-12);
        close(
            step.copper_heat_j,
            dt * p.winding_resistance_ohm * step.phase_current_a.iter().map(|i| i * i).sum::<f64>(),
            3e-13,
        );
        let dc_current: f64 = duty
            .iter()
            .zip(step.phase_current_a)
            .map(|(d, i)| d * i)
            .sum();
        close(step.dc_bridge_current_a, dc_current, 3e-13);
        close(step.bridge_work_j, dt * ports[0] * dc_current, 3e-13);
        close(
            step.machine_electrical_work_j,
            step.bridge_work_j - heat(&step.channel_heat_j),
            3e-13,
        );
        close(
            step.native.dissipated,
            step.copper_heat_j + step.bearing_heat_j + heat(&step.channel_heat_j),
            3e-13,
        );
        close(
            step.native.supplied,
            step.source_work_j + step.shaft_work_j.iter().sum::<f64>(),
            3e-13,
        );
        close(step.native.balance_residual(), 0.0, 2e-10);
        totals.push(heat(&step.channel_heat_j));
    }
    // Common-mode switching changes high/low occupation and hence losses when
    // the physical channels differ; ideal common-mode invariance must not leak in.
    assert!((totals[1] - totals[0]).abs() > 1e-8);
}

#[test]
fn changing_duty_on_one_warm_workspace_matches_a_fresh_native_owner() {
    let bridge = InverterParameters {
        high_on_resistance_ohm: [0.06, 0.09, 0.12],
        low_on_resistance_ohm: [0.03, 0.15, 0.05],
        ..ideal()
    };
    let (mut warm, initial) = prepared(machine(), bridge, [1.1, -1.4], 24.0);
    let (mut start, mut unused) = ([0.0; 6], [0.0; 3]);
    warm.step_controlled(
        &initial,
        [0.0; 3],
        &[0.8, -0.02, 0.01],
        3e-4,
        &mut start,
        &mut unused,
        || false,
    )
    .unwrap();
    let (mut fresh, _) = prepared(machine(), bridge, [0.0; 2], 24.0);
    let duty = [0.95, 0.1, 0.65];
    let inputs = [0.2, -0.15, 0.0];
    let (mut a, mut b, mut ya, mut yb) = ([0.0; 6], [0.0; 6], [0.0; 3], [0.0; 3]);
    let wa = warm
        .step_controlled(&start, duty, &inputs, 7e-4, &mut a, &mut ya, || false)
        .unwrap();
    let wb = fresh
        .step_controlled(&start, duty, &inputs, 7e-4, &mut b, &mut yb, || false)
        .unwrap();
    for i in 0..6 {
        close(a[i], b[i], 3e-12);
    }
    for i in 0..3 {
        close(ya[i], yb[i], 3e-12);
    }
    close(wa.bridge_work_j, wb.bridge_work_j, 3e-13);
    close(heat(&wa.channel_heat_j), heat(&wb.channel_heat_j), 3e-13);
    assert_eq!(stamp(&warm), stamp(&fresh));
}

#[test]
fn zero_and_rank_one_phase_losses_retain_the_original_positive_source_factor() {
    let p = SynchronousMotorParameters {
        winding_resistance_ohm: 0.0,
        bearing_drag_nm_s_per_rad: 0.0,
        ..machine()
    };
    let bridge = InverterParameters {
        high_on_resistance_ohm: [0.0, 0.2, 0.0],
        ..ideal()
    };
    let (mut model, initial) = prepared(p, bridge, [-0.7, 2.3], 24.0);
    for duty in [[0.0; 3], [0.0, 1.0, 0.0]] {
        let (mut next, mut ports) = ([0.0; 6], [0.0; 3]);
        let step = model
            .step_controlled(
                &initial,
                duty,
                &[0.0; 3],
                1e-4,
                &mut next,
                &mut ports,
                || false,
            )
            .unwrap();
        let (_, resistance, _) = model.system().structure();
        let root = model.system().resistance_witness().unwrap().rows();
        for i in 0..6 {
            for j in 0..6 {
                let gram: f64 = (0..6).map(|k| root[6 * k + i] * root[6 * k + j]).sum();
                close(gram, resistance[6 * i + j], 3e-14);
            }
        }
        assert_eq!(step.copper_heat_j, 0.0);
        assert_eq!(step.bearing_heat_j, 0.0);
        close(step.native.dissipated, heat(&step.channel_heat_j), 3e-14);
        if duty[1] == 0.0 {
            assert_eq!(heat(&step.channel_heat_j), 0.0);
            assert!(resistance.iter().all(|x| *x == 0.0));
        } else {
            assert!(step.channel_heat_j[1][0] > 0.0);
            close(heat(&step.channel_heat_j), step.channel_heat_j[1][0], 1e-15);
        }
    }
}

#[test]
fn invalid_controls_and_late_cancellation_preserve_state_operator_and_source_witness() {
    let bridge = InverterParameters {
        high_on_resistance_ohm: [0.06, 0.09, 0.12],
        low_on_resistance_ohm: [0.03, 0.15, 0.05],
        ..ideal()
    };
    let (mut model, initial) = prepared(machine(), bridge, [1.1, -1.4], 24.0);
    let (mut accepted, mut output) = ([0.0; 6], [0.0; 3]);
    model
        .step_controlled(
            &initial,
            [0.2, 0.6, 0.4],
            &[0.0; 3],
            2e-4,
            &mut accepted,
            &mut output,
            || false,
        )
        .unwrap();
    let original_stamp = stamp(&model);
    for bad in [[-0.01, 0.5, 0.5], [0.5, 1.01, 0.5], [0.5, 0.5, f64::NAN]] {
        let (mut next, mut ports) = ([71.0; 6], [83.0; 3]);
        assert!(
            model
                .step_controlled(
                    &accepted,
                    bad,
                    &[0.0; 3],
                    2e-4,
                    &mut next,
                    &mut ports,
                    || false
                )
                .is_err()
        );
        assert_eq!(next, [71.0; 6]);
        assert_eq!(ports, [83.0; 3]);
        assert_eq!(stamp(&model), original_stamp);
    }
    let duty = [0.85, 0.15, 0.7];
    let (mut reference, _) = prepared(machine(), bridge, [0.0; 2], 24.0);
    let (mut expected, mut expected_ports) = ([0.0; 6], [0.0; 3]);
    let mut calls = 0;
    reference
        .step_controlled(
            &accepted,
            duty,
            &[0.0; 3],
            2e-4,
            &mut expected,
            &mut expected_ports,
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    let (mut next, mut ports) = ([71.0; 6], [83.0; 3]);
    let mut repeated = 0;
    assert!(matches!(
        model.step_controlled(
            &accepted,
            duty,
            &[0.0; 3],
            2e-4,
            &mut next,
            &mut ports,
            || {
                repeated += 1;
                repeated == calls
            }
        ),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!(next, [71.0; 6]);
    assert_eq!(ports, [83.0; 3]);
    assert_eq!(stamp(&model), original_stamp);
    model
        .step_controlled(
            &accepted,
            duty,
            &[0.0; 3],
            2e-4,
            &mut next,
            &mut ports,
            || false,
        )
        .unwrap();
    for i in 0..6 {
        close(next[i], expected[i], 3e-12);
    }
    for i in 0..3 {
        close(ports[i], expected_ports[i], 3e-12);
    }
}

#[test]
fn zero_voltage_bus_can_charge_but_a_negative_candidate_refuses_without_clipping() {
    let p = SynchronousMotorParameters {
        pm_flux_linkage_wb: 0.0,
        winding_resistance_ohm: 0.0,
        bearing_drag_nm_s_per_rad: 0.0,
        ..machine()
    };
    let (mut model, initial) = prepared(p, ideal(), [0.0; 2], 0.0);
    let (mut next, mut ports) = ([0.0; 6], [0.0; 3]);
    model
        .step_controlled(
            &initial,
            [0.2, 0.6, 0.4],
            &[0.0; 3],
            1e-4,
            &mut next,
            &mut ports,
            || false,
        )
        .unwrap();
    let original = stamp(&model);
    next = [71.0; 6];
    ports = [83.0; 3];
    assert!(
        model
            .step_controlled(
                &initial,
                [1.0; 3],
                &[-1.0, 0.0, 0.0],
                1e-3,
                &mut next,
                &mut ports,
                || false
            )
            .is_err()
    );
    assert_eq!(next, [71.0; 6]);
    assert_eq!(ports, [83.0; 3]);
    assert_eq!(stamp(&model), original);
    let step = model
        .step_controlled(
            &initial,
            [1.0; 3],
            &[1.0, 0.0, 0.0],
            1e-3,
            &mut next,
            &mut ports,
            || false,
        )
        .unwrap();
    close(next[0], 1e-3, 3e-14);
    close(
        step.source_work_j,
        0.5 * 1e-6 / ideal().dc_capacitance_f,
        3e-14,
    );
    close(step.dc_bridge_current_a, 0.0, 3e-14);
    assert!(
        AveragedSynchronousInverter::new(
            SynchronousMotor::new(p).unwrap(),
            InverterParameters {
                dc_capacitance_f: 0.0,
                ..ideal()
            }
        )
        .is_err()
    );
    assert!(
        AveragedSynchronousInverter::new(
            SynchronousMotor::new(p).unwrap(),
            InverterParameters {
                high_on_resistance_ohm: [0.0, -0.1, 0.0],
                ..ideal()
            }
        )
        .is_err()
    );
}
