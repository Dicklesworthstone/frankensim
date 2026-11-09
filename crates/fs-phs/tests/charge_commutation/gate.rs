//! Independent gate-driver, reciprocal Miller-storage and complete-drive fixtures.

use super::*;
use fs_phs::circuit::converter::commutation::capacitance::PreparedGateCommutatedSwitchCycle;
use fs_phs::circuit::converter::commutation::capacitance::battery_motor::PreparedGateBatteryMotorCycle;
use fs_phs::circuit::converter::commutation::capacitance::gate::GateCommutationParameters;
use fs_phs::circuit::converter::commutation::gate::{
    GateChannelResponse, GateConductanceCurve, GateDriver,
};

#[path = "gate/current.rs"]
mod current;

const GATE_CAP: [f64; 4] = [0.003, 0.004, 0.005, 0.006];
const MILLER_CAP: [f64; 4] = [0.001, 0.002, 0.0015, 0.0025];
const DRIVER: GateDriver = GateDriver {
    on_voltage_v: 3.0,
    off_voltage_v: 0.0,
    on_resistance_ohm: 2.0,
    off_resistance_ohm: 3.0,
    maximum_current_a: 100.0,
};

fn conductance() -> GateConductanceCurve {
    GateConductanceCurve::try_new(vec![
        (-4.0, 0.025),
        (0.0, 0.025),
        (1.0, 0.025),
        (2.0, 5.0),
        (4.0, 5.0),
    ])
    .unwrap()
}

fn gate_parameters(cgs: [f64; 4], cgd: [f64; 4]) -> GateCommutationParameters {
    GateCommutationParameters {
        converter: drive_parameters().converter,
        intrinsic_drain_charge_curves: DRIVE_CAP.map(curve),
        gate_source_capacitance_f: cgs,
        gate_drain_capacitance_f: cgd,
        channel_conductance: core::array::from_fn(|_| conductance()),
        gate_drivers: [DRIVER; 4],
        maximum_cells_per_chord: 64,
    }
}

// Independent inversion of this fixture's nine original source cells per leg.
// Eliminating a gate at fixed charge adds Cgs*Cgd/(Cgs+Cgd) to its drain
// differential capacitance and shifts the two electrode charges reciprocally.
fn gate_physical<C: GateChannelResponse + Clone + 'static>(
    p: &GateCommutationParameters<C>,
    x: &[f64; 9],
) -> [f64; 9] {
    let mut out = [0.0; 9];
    out[0] = x[0] / p.converter.inductance_h;
    for k in 0..2 {
        let h = 2 * k;
        let l = h + 1;
        let sum = [
            p.gate_source_capacitance_f[h] + p.gate_drain_capacitance_f[h],
            p.gate_source_capacitance_f[l] + p.gate_drain_capacitance_f[l],
        ];
        let ratio = [
            p.gate_drain_capacitance_f[h] / sum[0],
            p.gate_drain_capacitance_f[l] / sum[1],
        ];
        let shifted = [
            x[1 + k] + ratio[0] * x[5 + h],
            x[3 + k] + ratio[1] * x[5 + l] - ratio[0] * x[5 + h],
        ];
        let mut found = false;
        for high in segments(DRIVE_CAP[h]) {
            for low in segments(DRIVE_CAP[l]) {
                let mut high = high;
                let mut low = low;
                high.c += p.gate_source_capacitance_f[h] * ratio[0];
                low.c += p.gate_source_capacitance_f[l] * ratio[1];
                let v = in_cell(
                    shifted,
                    high,
                    low,
                    inverse(p.converter.bus_capacitance_f[k], high, low),
                );
                let vd = v[0] - v[1];
                if vd >= high.low - 1e-11
                    && vd <= high.high + 1e-11
                    && v[1] >= low.low - 1e-11
                    && v[1] <= low.high + 1e-11
                {
                    out[1 + k] = v[0];
                    out[3 + k] = v[1];
                    out[5 + h] = (x[5 + h] + p.gate_drain_capacitance_f[h] * vd) / sum[0];
                    out[5 + l] = (x[5 + l] + p.gate_drain_capacitance_f[l] * v[1]) / sum[1];
                    found = true;
                    break;
                }
            }
            if found {
                break;
            }
        }
        assert!(
            found,
            "fixture charge must be supported by an original cell"
        );
    }
    out
}

fn gate_parts<C: GateChannelResponse + Clone + 'static>(
    p: &GateCommutationParameters<C>,
    x: &[f64; 9],
) -> [f64; 15] {
    let v = gate_physical(p, x);
    let vd = [v[1] - v[3], v[3], v[2] - v[4], v[4]];
    let mut parts = [0.0; 15];
    parts[0] = 0.5 * p.converter.inductance_h * v[0] * v[0];
    for i in 0..2 {
        parts[1 + i] = 0.5 * p.converter.bus_capacitance_f[i] * v[1 + i] * v[1 + i];
    }
    for i in 0..4 {
        parts[3 + i] = energy(DRIVE_CAP[i], vd[i]);
        parts[7 + i] = 0.5 * p.gate_source_capacitance_f[i] * v[5 + i] * v[5 + i];
        parts[11 + i] = 0.5 * p.gate_drain_capacitance_f[i] * (v[5 + i] - vd[i]).powi(2);
    }
    parts
}

fn prepare_gate(
    p: GateCommutationParameters,
    gates: [bool; 4],
    dt: f64,
    n: usize,
) -> PreparedGateCommutatedSwitchCycle {
    PreparedGateCommutatedSwitchCycle::new_gate_controlled(
        p,
        [1.0; 2],
        &[interval(gates, dt)],
        n,
        electrical(),
        controls(),
        || false,
    )
    .unwrap()
}

#[test]
fn four_finite_gate_drivers_recover_rc_source_work_heat_and_physical_refinement() {
    let law = conductance();
    let point = law.evaluate(1.5).unwrap();
    close(point.conductance_s, 2.5125, 1e-14);
    close(point.slope_s_per_v, 4.975, 1e-14);
    let epsilon = 1e-6;
    close(
        (law.evaluate(1.5 + epsilon).unwrap().conductance_s
            - law.evaluate(1.5 - epsilon).unwrap().conductance_s)
            / (2.0 * epsilon),
        point.slope_s_per_v,
        2e-9,
    );
    assert!(law.evaluate(4.01).is_err());
    for on in [true, false] {
        let p = DRIVER.evaluate(on, 0.7).unwrap();
        close(
            p.source_power_w - p.terminal_power_w,
            p.resistor_power_w,
            2e-15,
        );
        assert_eq!(p.current_a.is_sign_positive(), on);
    }
    // The second scale has one-nanofarad gates and a four-nanosecond interval.
    // Both scales compare time refinement to the continuous exponential, not
    // merely to a discrete energy identity.
    for scale in [1.0, 1e-7] {
        let c = 0.01 * scale;
        let tau = DRIVER.on_resistance_ohm * c;
        let duration = 2.0 * tau;
        let exact = DRIVER.on_voltage_v * (1.0 - (-2.0f64).exp());
        let mut errors = Vec::new();
        for n in [8, 16, 32] {
            let p = gate_parameters([c; 4], [0.0; 4]);
            let mut solve = electrical();
            solve.solve.absolute *= scale;
            solve.absolute_energy_tolerance_j *= scale;
            let mut cycle = controls();
            cycle.absolute_energy_tolerance_j *= scale;
            let mut native = PreparedGateCommutatedSwitchCycle::new_gate_controlled(
                p.clone(),
                [1.0; 2],
                &[interval([true; 4], duration)],
                n,
                solve,
                cycle,
                || false,
            )
            .unwrap();
            let initial = p.state(0.0, [0.0; 2], [0.0; 2], [0.0; 4]).unwrap();
            let mut next = [0.0; 9];
            let step = native
                .step_controlled(
                    &initial,
                    &[0.0; 2],
                    &[None; 2],
                    &mut next,
                    &mut [0.0; 2],
                    || false,
                )
                .unwrap();
            let v = gate_physical(&p, &next);
            let g = step.converter.gate.as_ref().unwrap();
            let ratio = (2.0 * tau - duration / n as f64) / (2.0 * tau + duration / n as f64);
            let discrete = DRIVER.on_voltage_v * (1.0 - ratio.powi(n as i32));
            for i in 0..4 {
                close(v[5 + i], discrete, 2e-11);
                close(
                    g.source_work_j[i],
                    DRIVER.on_voltage_v * c * v[5 + i],
                    2e-11 * scale,
                );
                close(
                    g.gate_source_energy_change_j[i],
                    0.5 * c * v[5 + i] * v[5 + i],
                    2e-11 * scale,
                );
                close(
                    g.resistor_heat_j[i],
                    c * (DRIVER.on_voltage_v * v[5 + i] - 0.5 * v[5 + i] * v[5 + i]),
                    2e-11 * scale,
                );
                close(g.gate_drain_energy_change_j[i], 0.0, 2e-13 * scale);
            }
            for value in step.converter.heat_joules() {
                close(value, 0.0, 2e-13 * scale);
            }
            assert_eq!(step.converter.all_heat_joules().len(), 13);
            errors.push((v[5] - exact).abs());
        }
        for pair in errors.windows(2) {
            assert!(pair[0] / pair[1] > 3.8 && pair[0] / pair[1] < 4.2);
        }
        assert!(errors[2] < 3e-4);
    }
}

#[test]
fn a_command_edge_charges_the_gate_before_channel_conduction_changes() {
    let p = gate_parameters([0.01; 4], [0.0; 4]);
    let dt = 0.005;
    let initial = p.state(0.2, [8.0, 5.0], [4.0, 2.0], [0.0; 4]).unwrap();
    let command = PowerSwitchCommand {
        source_voltage_v: [8.0, 5.0],
        loads: [None; 2],
    };
    let mut actual = prepare_gate(p.clone(), [true; 4], dt, 1);
    let mut next = [0.0; 9];
    let step = actual
        .step_controlled(
            &initial,
            &command.source_voltage_v,
            &command.loads,
            &mut next,
            &mut [0.0; 2],
            || false,
        )
        .unwrap();
    let gate = step.converter.gate.as_ref().unwrap();
    // All four commands are on. The actual gate trajectory is still below
    // the original channel law's threshold, so each channel retains Goff.
    for v in gate.mean_voltage_v {
        assert!(v > 0.0 && v < 1.0);
    }
    for v in &gate_physical(&p, &next)[5..] {
        assert!(*v > 0.0 && *v < 1.0);
    }
    let baseline = ChargeCommutationParameters {
        converter: p.converter,
        output_charge_curves: p.intrinsic_drain_charge_curves.clone(),
        maximum_cells_per_chord: 64,
    };
    let mut old = PreparedChargeCommutatedSwitchCycle::new_charge_controlled(
        baseline,
        [1.0; 2],
        &[interval([false; 4], dt)],
        1,
        electrical(),
        controls(),
        || false,
    )
    .unwrap();
    let old_initial: [f64; 5] = initial[..5].try_into().unwrap();
    let mut old_next = [0.0; 5];
    let old_step = old
        .step_controlled(
            &old_initial,
            &command.source_voltage_v,
            &command.loads,
            &mut old_next,
            &mut [0.0; 2],
            || false,
        )
        .unwrap();
    for i in 0..5 {
        close(next[i], old_next[i], 3e-11);
    }
    for (a, b) in step
        .converter
        .heat_joules()
        .into_iter()
        .zip(old_step.converter.heat_joules())
    {
        close(a, b, 3e-11);
    }
    // Prescribing the following off edge continues all nine physical charges.
    let deck = [interval([true; 4], dt), interval([false; 4], dt)];
    let mut whole = PreparedGateCommutatedSwitchCycle::new_gate_controlled(
        p.clone(),
        [1.0; 2],
        &deck,
        1,
        electrical(),
        controls(),
        || false,
    )
    .unwrap();
    let mut second = prepare_gate(p.clone(), [false; 4], dt, 1);
    let mut split = [0.0; 9];
    second
        .step_controlled(
            &next,
            &command.source_voltage_v,
            &command.loads,
            &mut split,
            &mut [0.0; 2],
            || false,
        )
        .unwrap();
    let mut joined = [0.0; 9];
    whole
        .step_interval_commands_controlled(
            &initial,
            &[command; 2],
            &mut joined,
            &mut [0.0; 2],
            || false,
        )
        .unwrap();
    assert_eq!(joined.map(f64::to_bits), split.map(f64::to_bits));
    assert!(gate_physical(&p, &split)[5] > 0.0);
}

fn gate_drive_state<C: GateChannelResponse + Clone + 'static>(
    native: &PreparedGateBatteryMotorCycle<C>,
    p: &GateCommutationParameters<C>,
    y: [f64; 11],
    gates: [f64; 4],
) -> [f64; 15] {
    let converter = p.state(y[3], [y[4], y[5]], [y[6], y[7]], gates).unwrap();
    native
        .state(
            &[
                y[0],
                BATTERY.polarization_capacitance_f * y[1],
                BATTERY.bus_capacitance_f * y[2],
            ],
            &converter,
            &[
                MOTOR.inductance_h * y[8],
                MOTOR.rotor_inertia_kg_m2 * y[9],
                MOTOR.stator_inertia_kg_m2 * y[10],
            ],
        )
        .unwrap()
}

fn gate_drive_parts<C: GateChannelResponse + Clone + 'static>(
    p: &GateCommutationParameters<C>,
    x: &[f64; 15],
) -> [f64; 7] {
    let converter: [f64; 9] = x[3..12].try_into().unwrap();
    [
        12.0 * x[0],
        0.5 * x[1] * x[1] / BATTERY.polarization_capacitance_f,
        0.5 * x[2] * x[2] / BATTERY.bus_capacitance_f,
        gate_parts(p, &converter).iter().sum(),
        0.5 * x[12] * x[12] / MOTOR.inductance_h,
        0.5 * x[13] * x[13] / MOTOR.rotor_inertia_kg_m2,
        0.5 * x[14] * x[14] / MOTOR.stator_inertia_kg_m2,
    ]
}

#[test]
fn fifteen_state_miller_drive_uses_one_effort_for_all_kcl_storage_work_and_eighteen_heats() {
    let p = gate_parameters(GATE_CAP, MILLER_CAP);
    let dt = 1.0 / 1024.0;
    let gates = [true, false, true, false];
    let mut native = PreparedGateBatteryMotorCycle::new_gate_controlled(
        drive_model(),
        p.clone(),
        &[interval(gates, dt)],
        1,
        electrical(),
        controls(),
        || false,
    )
    .unwrap();
    let before = gate_drive_state(
        &native,
        &p,
        [52.0, 0.1, 12.1, 0.5, 11.8, 6.4, 10.75, 4.0, 0.1, 4.0, -0.4],
        [1.4, 0.2, 1.3, 0.1],
    );
    let command = BatteryMotorSwitchCommand {
        shaft_torque_nm: [-0.2, 0.03],
        auxiliary_load: None,
    };
    let mut next = [0.0; 15];
    let step = native
        .step_controlled(&before, command, &mut next, || false)
        .unwrap();
    let c = &step.converter;
    let g = c.gate.as_ref().unwrap();
    let i = c.mean_current_a;
    let vi = c.mean_bus_voltage_v[0];
    let vo = c.mean_bus_voltage_v[1];
    let left = c.mean_node_voltage_v[0];
    let right = c.mean_node_voltage_v[1];
    let vb = step.mean_bus_voltage_v[0];
    let pol = f64::midpoint(before[1], next[1]) / BATTERY.polarization_capacitance_f;
    let cell = (12.0 - pol - vb) / BATTERY.series_resistance_ohm;
    let lead = (vb - vi) / LEAD;
    let im = step.mean_motor_current_a;
    let relative = step.mean_shaft_speed_rad_s[0] - step.mean_shaft_speed_rad_s[1];
    let vd = [vi - left, left, vo - right, right];
    // This is the explicit authored G(Vgs) law, independent of its evaluator.
    let gs = g.mean_voltage_v.map(|v| {
        if v <= 1.0 {
            0.025
        } else if v >= 2.0 {
            5.0
        } else {
            0.025 + 4.975 * (v - 1.0)
        }
    });
    let channel = core::array::from_fn::<_, 4, _>(|k| gs[k] * vd[k]);
    let diode = vd.map(|v| ((-v - 0.5) / 0.1).max(0.0));
    let current = core::array::from_fn::<_, 4, _>(|k| channel[k] - diode[k]);
    let driver = core::array::from_fn::<_, 4, _>(|k| {
        let u = if gates[k] {
            DRIVER.on_voltage_v
        } else {
            DRIVER.off_voltage_v
        };
        let r = if gates[k] {
            DRIVER.on_resistance_ohm
        } else {
            DRIVER.off_resistance_ohm
        };
        (u - g.mean_voltage_v[k]) / r
    });
    let torque = MOTOR.coupling_nm_per_a * im - MOTOR.bearing_drag_nm_s_per_rad * relative;
    let rates = [
        -cell,
        cell - pol / BATTERY.polarization_resistance_ohm,
        cell - lead,
        left - right - DRIVE_CONVERTER.winding_resistance_ohm * i,
        lead - current[0],
        -im - current[2],
        current[0] - current[1] - i,
        current[2] - current[3] + i,
        driver[0],
        driver[1],
        driver[2],
        driver[3],
        vo - MOTOR.winding_resistance_ohm * im - MOTOR.coupling_nm_per_a * relative,
        torque + command.shaft_torque_nm[0],
        -torque + command.shaft_torque_nm[1],
    ];
    for k in 0..15 {
        close(next[k] - before[k], dt * rates[k], 4e-10);
    }
    close(step.mean_lead_current_a, lead, 2e-12);
    let mut power = [0.0; 18];
    power[0] = DRIVE_CONVERTER.winding_resistance_ohm * i * i;
    for k in 0..4 {
        power[1 + k] = vd[k] * channel[k];
        power[10 + k] = -vd[k] * diode[k];
        let r = if gates[k] {
            DRIVER.on_resistance_ohm
        } else {
            DRIVER.off_resistance_ohm
        };
        let u = if gates[k] {
            DRIVER.on_voltage_v
        } else {
            DRIVER.off_voltage_v
        };
        power[14 + k] = r * driver[k] * driver[k];
        close(g.source_work_j[k], dt * u * driver[k], 2e-12);
        close(g.mean_current_a[k], driver[k], 2e-12);
    }
    power[5] = BATTERY.series_resistance_ohm * cell * cell;
    power[6] = pol * pol / BATTERY.polarization_resistance_ohm;
    power[7] = LEAD * lead * lead;
    power[8] = MOTOR.winding_resistance_ohm * im * im;
    power[9] = MOTOR.bearing_drag_nm_s_per_rad * relative * relative;
    let all_heat = step.all_heat_joules();
    assert_eq!(all_heat.len(), 18);
    for (q, power) in all_heat.iter().zip(power) {
        close(*q, dt * power, 3e-10);
    }
    assert!(all_heat[14..].iter().all(|q| *q > 0.0));
    let a = gate_drive_parts(&p, &before);
    let b = gate_drive_parts(&p, &next);
    for k in 0..7 {
        close(step.storage_change_j[k], b[k] - a[k], 3e-10);
    }
    let ca: [f64; 9] = before[3..12].try_into().unwrap();
    let cb: [f64; 9] = next[3..12].try_into().unwrap();
    let before_parts = gate_parts(&p, &ca);
    let after_parts = gate_parts(&p, &cb);
    for k in 0..4 {
        close(
            c.device_capacitor_energy_change_j[k],
            after_parts[3 + k] - before_parts[3 + k],
            3e-11,
        );
        close(
            g.gate_source_energy_change_j[k],
            after_parts[7 + k] - before_parts[7 + k],
            3e-11,
        );
        close(
            g.gate_drain_energy_change_j[k],
            after_parts[11 + k] - before_parts[11 + k],
            3e-11,
        );
    }
    close(step.battery_terminal_work_j, -dt * vb * lead, 3e-10);
    close(c.terminal_work_j[0], dt * vi * lead, 3e-10);
    close(c.terminal_work_j[1], -dt * vo * im, 3e-10);
    close(step.motor_electrical_work_j, dt * vo * im, 3e-10);
    close(
        step.electromagnetic_work_j,
        dt * MOTOR.coupling_nm_per_a * relative * im,
        3e-10,
    );
    for k in 0..2 {
        close(
            step.shaft_work_j[k],
            dt * command.shaft_torque_nm[k] * step.mean_shaft_speed_rad_s[k],
            3e-10,
        );
    }
    let independent_change: f64 = b.into_iter().zip(a).map(|(b, a)| b - a).sum();
    close(
        independent_change + all_heat.iter().sum::<f64>()
            - g.source_work_j.iter().sum::<f64>()
            - step.shaft_work_j.iter().sum::<f64>(),
        0.0,
        6e-10,
    );
    close(step.native.delta_h, independent_change, 6e-10);
    close(step.native.dissipated, all_heat.iter().sum(), 3e-10);
}

#[test]
fn fifteen_state_steady_regeneration_recharges_the_finite_battery_without_resetting_gates() {
    let p = gate_parameters(GATE_CAP, MILLER_CAP);
    let cell = -2.0;
    let pol = BATTERY.polarization_resistance_ohm * cell;
    let vb = 12.0 - pol - BATTERY.series_resistance_ohm * cell;
    let vi = vb - LEAD * cell;
    let left = vi - DEVICE.on_resistance_ohm * cell;
    let ic = cell - left / DEVICE.off_resistance_ohm;
    let right = left - DRIVE_CONVERTER.winding_resistance_ohm * ic;
    let output_high = -ic + right / DEVICE.off_resistance_ohm;
    let vo = right + DEVICE.on_resistance_ohm * output_high;
    let im = -output_high;
    let relative = (vo - MOTOR.winding_resistance_ohm * im) / MOTOR.coupling_nm_per_a;
    let reaction = -MOTOR.coupling_nm_per_a * im + MOTOR.bearing_drag_nm_s_per_rad * relative;
    let ws = -2.0;
    let wr = ws + relative;
    let dt = 1.0 / 128.0;
    let mut native = PreparedGateBatteryMotorCycle::new_gate_controlled(
        drive_model(),
        p.clone(),
        &[interval([true, false, true, false], dt)],
        4,
        electrical(),
        controls(),
        || false,
    )
    .unwrap();
    let before = gate_drive_state(
        &native,
        &p,
        [40.0, pol, vb, ic, vi, vo, left, right, im, wr, ws],
        [3.0, 0.0, 3.0, 0.0],
    );
    let command = BatteryMotorSwitchCommand {
        shaft_torque_nm: [reaction, -reaction],
        auxiliary_load: None,
    };
    let mut next = [0.0; 15];
    let step = native
        .step_controlled(&before, command, &mut next, || false)
        .unwrap();
    close(next[0] - before[0], -cell * dt, 2e-11);
    for k in 1..15 {
        close(next[k], before[k], 2e-11);
    }
    let a = gate_drive_parts(&p, &before);
    let b = gate_drive_parts(&p, &next);
    for k in 0..7 {
        close(step.storage_change_j[k], b[k] - a[k], 3e-10);
    }
    let power = [
        DRIVE_CONVERTER.winding_resistance_ohm * ic * ic,
        DEVICE.on_resistance_ohm * cell * cell,
        left * left / DEVICE.off_resistance_ohm,
        DEVICE.on_resistance_ohm * output_high * output_high,
        right * right / DEVICE.off_resistance_ohm,
        BATTERY.series_resistance_ohm * cell * cell,
        pol * pol / BATTERY.polarization_resistance_ohm,
        LEAD * cell * cell,
        MOTOR.winding_resistance_ohm * im * im,
        MOTOR.bearing_drag_nm_s_per_rad * relative * relative,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ];
    for (q, power) in step.all_heat_joules().into_iter().zip(power) {
        close(q, dt * power, 3e-10);
    }
    let g = step.converter.gate.as_ref().unwrap();
    for k in 0..4 {
        close(g.source_work_j[k], 0.0, 2e-12);
        close(g.gate_source_energy_change_j[k], 0.0, 2e-12);
        close(g.gate_drain_energy_change_j[k], 0.0, 2e-12);
    }
    assert!(step.storage_change_j[0] > 0.0 && step.motor_electrical_work_j < 0.0);
    assert!(step.shaft_work_j.iter().all(|w| *w > 0.0));
    close(step.balance_residual_j(), 0.0, 6e-10);
}

#[test]
fn finite_gate_support_current_limits_and_final_cancellation_are_transactional() {
    let p = gate_parameters([0.01; 4], [0.0; 4]);
    for count in [0, usize::MAX] {
        assert!(
            PreparedGateCommutatedSwitchCycle::new_gate_controlled(
                p.clone(),
                [1.0; 2],
                &[interval([true; 4], 0.005)],
                count,
                electrical(),
                controls(),
                || false,
            )
            .is_err()
        );
    }
    let initial = p.state(0.0, [0.0; 2], [0.0; 2], [0.0; 4]).unwrap();
    let deck = [interval([false; 4], 0.001), interval([true; 4], 0.04)];
    let sentinel = [73.0; 9];
    for limited_current in [true, false] {
        let mut refused = p.clone();
        if limited_current {
            refused.gate_drivers[0].maximum_current_a = 0.01;
        } else {
            refused.channel_conductance[0] =
                GateConductanceCurve::try_new(vec![(-1.0, 0.025), (0.5, 0.025)]).unwrap();
        }
        let mut native = PreparedGateCommutatedSwitchCycle::new_gate_controlled(
            refused,
            [1.0; 2],
            &deck,
            1,
            electrical(),
            controls(),
            || false,
        )
        .unwrap();
        let mut next = sentinel;
        let mut bus = [73.0; 2];
        assert!(
            native
                .step_controlled(&initial, &[0.0; 2], &[None; 2], &mut next, &mut bus, || {
                    false
                })
                .is_err()
        );
        assert_eq!(next, sentinel);
        assert_eq!(bus, [73.0; 2]);
    }
    let mut native = prepare_gate(p, [true; 4], 0.005, 2);
    let polls = Cell::new(0usize);
    let mut expected = [0.0; 9];
    native
        .step_controlled(
            &initial,
            &[0.0; 2],
            &[None; 2],
            &mut expected,
            &mut [0.0; 2],
            || {
                polls.set(polls.get() + 1);
                false
            },
        )
        .unwrap();
    let limit = polls.get() - 1;
    polls.set(0);
    let mut next = sentinel;
    let mut bus = [73.0; 2];
    assert!(matches!(
        native.step_controlled(&initial, &[0.0; 2], &[None; 2], &mut next, &mut bus, || {
            polls.set(polls.get() + 1);
            polls.get() >= limit
        }),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!(next, sentinel);
    assert_eq!(bus, [73.0; 2]);
    native
        .step_controlled(
            &initial,
            &[0.0; 2],
            &[None; 2],
            &mut next,
            &mut [0.0; 2],
            || false,
        )
        .unwrap();
    assert_eq!(next.map(f64::to_bits), expected.map(f64::to_bits));
    let p = gate_parameters(GATE_CAP, MILLER_CAP);
    let mut drive = PreparedGateBatteryMotorCycle::new_gate_controlled(
        drive_model(),
        p.clone(),
        &[interval([true, false, true, false], 1.0 / 1024.0)],
        1,
        electrical(),
        controls(),
        || false,
    )
    .unwrap();
    let before = gate_drive_state(
        &drive,
        &p,
        [52.0, 0.1, 12.1, 0.5, 11.8, 6.4, 10.75, 4.0, 0.1, 4.0, -0.4],
        [1.4, 0.2, 1.3, 0.1],
    );
    let command = BatteryMotorSwitchCommand {
        shaft_torque_nm: [-0.2, 0.03],
        auxiliary_load: None,
    };
    let mut accepted = [0.0; 15];
    polls.set(0);
    drive
        .step_controlled(&before, command, &mut accepted, || {
            polls.set(polls.get() + 1);
            false
        })
        .unwrap();
    let limit = polls.get() - 1;
    polls.set(0);
    let mut candidate = [73.0; 15];
    assert!(matches!(
        drive.step_controlled(&before, command, &mut candidate, || {
            polls.set(polls.get() + 1);
            polls.get() >= limit
        }),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!(candidate, [73.0; 15]);
    drive
        .step_controlled(&before, command, &mut candidate, || false)
        .unwrap();
    assert_eq!(candidate.map(f64::to_bits), accepted.map(f64::to_bits));
}
