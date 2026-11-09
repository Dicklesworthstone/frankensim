//! Independent current-source saturation, physical node KCL and complete power.

use super::*;
use fs_phs::PhsError;
use fs_phs::circuit::converter::commutation::gate::{GateChannelLaw, GateChannelPoint};

type Parameters = GateCommutationParameters<GateChannelLaw>;
type Native = PreparedGateBatteryMotorCycle<GateChannelLaw>;

#[derive(Clone, Debug)]
struct SaturatingChannel {
    gate_max: f64,
    active_above: Option<f64>,
}
impl GateChannelResponse for SaturatingChannel {
    fn gate_voltage_range_v(&self) -> [f64; 2] {
        [-1.0, self.gate_max]
    }
    fn evaluate(&self, vd: f64, vg: f64) -> Result<GateChannelPoint, PhsError> {
        if !vd.is_finite() || !vg.is_finite() || vd.abs() > 20.0 || vg < -1.0 || vg > self.gate_max
        {
            return Err(PhsError::Dimension {
                what: "authored square-law source domain",
            });
        }
        let u = (vg - 1.0).max(0.0);
        let a = vd.abs();
        let k = if vd < 0.0 { 8.0 } else { 20.0 };
        let (i, dd, dg) = if a < u {
            (k * (2.0 * u * a - a * a), 2.0 * k * (u - a), 2.0 * k * a)
        } else {
            (k * u * u, 0.0, 2.0 * k * u)
        };
        let sign = if self.active_above.is_some_and(|threshold| vg > threshold) {
            -1.0
        } else {
            1.0
        };
        Ok(GateChannelPoint {
            current_a: sign * vd.signum() * i,
            drain_derivative_s: sign * dd,
            gate_derivative_s: sign * vd.signum() * dg,
        })
    }
}
fn channel() -> GateChannelLaw {
    GateChannelLaw::new(SaturatingChannel {
        gate_max: 4.0,
        active_above: None,
    })
    .unwrap()
}
fn parameters() -> Parameters {
    let old = gate_parameters(GATE_CAP, MILLER_CAP);
    GateCommutationParameters {
        converter: old.converter,
        intrinsic_drain_charge_curves: old.intrinsic_drain_charge_curves,
        gate_source_capacitance_f: old.gate_source_capacitance_f,
        gate_drain_capacitance_f: old.gate_drain_capacitance_f,
        channel_conductance: core::array::from_fn(|_| channel()),
        gate_drivers: old.gate_drivers,
        maximum_cells_per_chord: old.maximum_cells_per_chord,
    }
}
fn prepare(p: Parameters, gates: [bool; 4], dt: f64, n: usize) -> Native {
    Native::new_gate_controlled(
        drive_model(),
        p,
        &[interval(gates, dt)],
        n,
        electrical(),
        controls(),
        || false,
    )
    .unwrap()
}
// Independent algebraic expression for the authored current, without calling
// the production source or converting I/V into an invented conductance.
fn authored_current(vd: f64, vg: f64) -> f64 {
    let u = (vg - 1.0).max(0.0);
    let scale = if vd < 0.0 { 8.0 } else { 20.0 };
    vd.signum() * scale * (u * u - (u - vd.abs()).max(0.0).powi(2))
}

struct PhysicalAccount {
    rate: [f64; 15],
    heat: [f64; 18],
    work: [f64; 12],
}
fn physical_account(
    p: &Parameters,
    gates: [bool; 4],
    v: [f64; 9],
    battery: [f64; 2],
    motor: [f64; 3],
    torques: [f64; 2],
) -> PhysicalAccount {
    let [pol, vb] = battery;
    let [im, wr, ws] = motor;
    let [i, vi, vo, left, right, _, _, _, _] = v;
    let relative = wr - ws;
    let cell = (12.0 - pol - vb) / BATTERY.series_resistance_ohm;
    let lead = (vb - vi) / LEAD;
    let vd = [vi - left, left, vo - right, right];
    let channels = core::array::from_fn::<_, 4, _>(|k| authored_current(vd[k], v[5 + k]));
    let diode = core::array::from_fn::<_, 4, _>(|k| ((-vd[k] - 0.5) / 0.1).max(0.0));
    let current = core::array::from_fn::<_, 4, _>(|k| channels[k] - diode[k]);
    let driver = core::array::from_fn::<_, 4, _>(|k| {
        let d = p.gate_drivers[k];
        let u = if gates[k] {
            d.on_voltage_v
        } else {
            d.off_voltage_v
        };
        let r = if gates[k] {
            d.on_resistance_ohm
        } else {
            d.off_resistance_ohm
        };
        (u - v[5 + k]) / r
    });
    let torque = MOTOR.coupling_nm_per_a * im - MOTOR.bearing_drag_nm_s_per_rad * relative;
    let rate = [
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
        torque + torques[0],
        -torque + torques[1],
    ];
    let mut heat = [0.0; 18];
    heat[0] = DRIVE_CONVERTER.winding_resistance_ohm * i * i;
    heat[5] = BATTERY.series_resistance_ohm * cell * cell;
    heat[6] = pol * pol / BATTERY.polarization_resistance_ohm;
    heat[7] = LEAD * lead * lead;
    heat[8] = MOTOR.winding_resistance_ohm * im * im;
    heat[9] = MOTOR.bearing_drag_nm_s_per_rad * relative * relative;
    let mut work = [
        -vb * lead,
        vi * lead,
        -vo * im,
        vo * im,
        MOTOR.coupling_nm_per_a * relative * im,
        torques[0] * wr,
        torques[1] * ws,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ];
    for k in 0..4 {
        let d = p.gate_drivers[k];
        let u = if gates[k] {
            d.on_voltage_v
        } else {
            d.off_voltage_v
        };
        let r = if gates[k] {
            d.on_resistance_ohm
        } else {
            d.off_resistance_ohm
        };
        heat[1 + k] = vd[k] * channels[k];
        heat[10 + k] = -vd[k] * diode[k];
        heat[14 + k] = r * driver[k] * driver[k];
        work[8 + k] = u * driver[k];
    }
    PhysicalAccount { rate, heat, work }
}
fn receipt_work(s: &CapacitiveBatteryMotorCycleStep) -> [f64; 12] {
    let gate = s.converter.gate.as_ref().unwrap();
    [
        s.battery_terminal_work_j,
        s.converter.terminal_work_j[0],
        s.converter.terminal_work_j[1],
        s.motor_electrical_work_j,
        s.electromagnetic_work_j,
        s.shaft_work_j[0],
        s.shaft_work_j[1],
        s.load_work_j,
        gate.source_work_j[0],
        gate.source_work_j[1],
        gate.source_work_j[2],
        gate.source_work_j[3],
    ]
}
fn rhs(p: &Parameters, x: &[f64; 45], gates: [bool; 4], torques: [f64; 2]) -> [f64; 45] {
    let converter: [f64; 9] = x[3..12].try_into().unwrap();
    let a = physical_account(
        p,
        gates,
        gate_physical(p, &converter),
        [
            x[1] / BATTERY.polarization_capacitance_f,
            x[2] / BATTERY.bus_capacitance_f,
        ],
        [
            x[12] / MOTOR.inductance_h,
            x[13] / MOTOR.rotor_inertia_kg_m2,
            x[14] / MOTOR.stator_inertia_kg_m2,
        ],
        torques,
    );
    let mut out = [0.0; 45];
    out[..15].copy_from_slice(&a.rate);
    out[15..33].copy_from_slice(&a.heat);
    out[33..].copy_from_slice(&a.work);
    out
}
fn rk4(
    p: &Parameters,
    initial: [f64; 15],
    gates: [bool; 4],
    torques: [f64; 2],
    dt: f64,
    n: usize,
) -> [f64; 45] {
    let h = dt / n as f64;
    let mut x = [0.0; 45];
    x[..15].copy_from_slice(&initial);
    for _ in 0..n {
        let a = rhs(p, &x, gates, torques);
        let b = rhs(
            p,
            &core::array::from_fn(|k| x[k] + 0.5 * h * a[k]),
            gates,
            torques,
        );
        let c = rhs(
            p,
            &core::array::from_fn(|k| x[k] + 0.5 * h * b[k]),
            gates,
            torques,
        );
        let d = rhs(
            p,
            &core::array::from_fn(|k| x[k] + h * c[k]),
            gates,
            torques,
        );
        for k in 0..45 {
            x[k] += h * (a[k] + 2.0 * b[k] + 2.0 * c[k] + d[k]) / 6.0;
        }
    }
    x
}
fn dynamic_state(native: &Native, p: &Parameters) -> [f64; 15] {
    gate_drive_state(
        native,
        p,
        [52.0, 0.1, 12.1, 0.5, 11.8, 6.4, 10.5, 7.1, 0.1, 4.0, -0.4],
        [1.4, 0.2, 1.3, 0.1],
    )
}

#[test]
fn saturated_forward_and_asymmetric_reverse_channels_keep_exact_steady_motoring_and_regeneration() {
    for cell in [1.0, -1.0] {
        let mut p = parameters();
        let high = if cell > 0.0 {
            [1.0 + (1.0f64 / 20.0).sqrt(), 1.0 + (1.0f64 / 8.0).sqrt()]
        } else {
            [1.0 + (1.0f64 / 8.0).sqrt(), 1.0 + (1.0f64 / 20.0).sqrt()]
        };
        p.gate_drivers[0].on_voltage_v = high[0];
        p.gate_drivers[2].on_voltage_v = high[1];
        let pol = BATTERY.polarization_resistance_ohm * cell;
        let vb = 12.0 - pol - BATTERY.series_resistance_ohm * cell;
        let vi = vb - LEAD * cell;
        let left = vi - 0.4 * cell;
        let right = left - DRIVE_CONVERTER.winding_resistance_ohm * cell;
        let vo = right - 0.4 * cell;
        let relative = (vo - MOTOR.winding_resistance_ohm * cell) / MOTOR.coupling_nm_per_a;
        let ws = -2.0;
        let wr = ws + relative;
        let reaction = -MOTOR.coupling_nm_per_a * cell + MOTOR.bearing_drag_nm_s_per_rad * relative;
        let command = BatteryMotorSwitchCommand {
            shaft_torque_nm: [reaction, -reaction],
            auxiliary_load: None,
        };
        let dt = 1.0 / 128.0;
        let gates = [true, false, true, false];
        let mut native = prepare(p.clone(), gates, dt, 4);
        let before = gate_drive_state(
            &native,
            &p,
            [40.0, pol, vb, cell, vi, vo, left, right, cell, wr, ws],
            [high[0], 0.0, high[1], 0.0],
        );
        let mut next = [0.0; 15];
        let step = native
            .step_controlled(&before, command, &mut next, || false)
            .unwrap();
        close(next[0] - before[0], -dt * cell, 3e-11);
        for k in 1..15 {
            close(next[k], before[k], 3e-11);
        }
        let a = physical_account(
            &p,
            gates,
            [cell, vi, vo, left, right, high[0], 0.0, high[1], 0.0],
            [pol, vb],
            [cell, wr, ws],
            command.shaft_torque_nm,
        );
        for (q, power) in step.all_heat_joules().into_iter().zip(a.heat) {
            close(q, dt * power, 3e-10);
        }
        for (w, power) in receipt_work(&step).into_iter().zip(a.work) {
            close(w, dt * power, 3e-10);
        }
        close(step.converter.channel_heat_j[0], 0.4 * dt, 2e-11);
        close(step.converter.channel_heat_j[2], 0.4 * dt, 2e-11);
        for k in [1, 3] {
            close(step.converter.channel_heat_j[k], 0.0, 2e-13);
        }
        assert!(step.converter.diode_heat_j.iter().all(|q| q.abs() < 2e-13));
        for k in [0, 2] {
            let vd = if k == 0 { 0.4 * cell } else { -0.4 * cell };
            let point = p.channel_conductance[k].evaluate(vd, high[k / 2]).unwrap();
            assert_eq!(point.drain_derivative_s, 0.0);
        }
        let initial_energy = gate_drive_parts(&p, &before);
        let final_energy = gate_drive_parts(&p, &next);
        for k in 0..7 {
            close(
                step.storage_change_j[k],
                final_energy[k] - initial_energy[k],
                3e-10,
            );
        }
        close(step.storage_change_j[0], -12.0 * cell * dt, 3e-10);
        assert_eq!(step.storage_change_j[0] > 0.0, cell < 0.0);
        assert_eq!(step.motor_electrical_work_j < 0.0, cell < 0.0);
        close(step.balance_residual_j(), 0.0, 6e-10);
    }
}

#[test]
fn actual_current_source_and_miller_efforts_own_all_fifteen_kcl_rates_and_eighteen_heats() {
    let p = parameters();
    let gates = [true, false, true, false];
    let dt = 1.0 / 2048.0;
    let mut native = prepare(p.clone(), gates, dt, 1);
    let before = dynamic_state(&native, &p);
    let command = BatteryMotorSwitchCommand {
        shaft_torque_nm: [-0.2, 0.03],
        auxiliary_load: None,
    };
    let mut next = [0.0; 15];
    let s = native
        .step_controlled(&before, command, &mut next, || false)
        .unwrap();
    let c = &s.converter;
    let g = c.gate.as_ref().unwrap();
    let v = [
        c.mean_current_a,
        c.mean_bus_voltage_v[0],
        c.mean_bus_voltage_v[1],
        c.mean_node_voltage_v[0],
        c.mean_node_voltage_v[1],
        g.mean_voltage_v[0],
        g.mean_voltage_v[1],
        g.mean_voltage_v[2],
        g.mean_voltage_v[3],
    ];
    let a = physical_account(
        &p,
        gates,
        v,
        [
            f64::midpoint(before[1], next[1]) / BATTERY.polarization_capacitance_f,
            s.mean_bus_voltage_v[0],
        ],
        [
            s.mean_motor_current_a,
            s.mean_shaft_speed_rad_s[0],
            s.mean_shaft_speed_rad_s[1],
        ],
        command.shaft_torque_nm,
    );
    for k in 0..15 {
        close(next[k] - before[k], dt * a.rate[k], 4e-10);
    }
    for (q, power) in s.all_heat_joules().into_iter().zip(a.heat) {
        close(q, dt * power, 3e-10);
    }
    for (w, power) in receipt_work(&s).into_iter().zip(a.work) {
        close(w, dt * power, 3e-10);
    }
    assert!(s.all_heat_joules()[12] > 0.0); // Output high body diode has its own work.
    assert!(g.resistor_heat_j.iter().all(|q| *q > 0.0));
    let ea = gate_drive_parts(&p, &before);
    let eb = gate_drive_parts(&p, &next);
    for k in 0..7 {
        close(s.storage_change_j[k], eb[k] - ea[k], 3e-10);
    }
    let ca: [f64; 9] = before[3..12].try_into().unwrap();
    let cb: [f64; 9] = next[3..12].try_into().unwrap();
    let qa = gate_parts(&p, &ca);
    let qb = gate_parts(&p, &cb);
    for k in 0..4 {
        close(
            c.device_capacitor_energy_change_j[k],
            qb[3 + k] - qa[3 + k],
            3e-11,
        );
        close(
            g.gate_source_energy_change_j[k],
            qb[7 + k] - qa[7 + k],
            3e-11,
        );
        close(
            g.gate_drain_energy_change_j[k],
            qb[11 + k] - qa[11 + k],
            3e-11,
        );
    }
    let dh: f64 = eb.into_iter().zip(ea).map(|(b, a)| b - a).sum();
    close(
        dh + s.all_heat_joules().iter().sum::<f64>()
            - g.source_work_j.iter().sum::<f64>()
            - s.shaft_work_j.iter().sum::<f64>(),
        0.0,
        6e-10,
    );
    close(s.native.delta_h, dh, 6e-10);
    close(s.native.dissipated, s.all_heat_joules().iter().sum(), 3e-10);
}

#[test]
fn saturated_channel_drive_states_and_separate_work_integrals_refine_to_independent_si_ode() {
    let p = parameters();
    let gates = [true, false, true, false];
    let dt = 1.0 / 1024.0;
    let command = BatteryMotorSwitchCommand {
        shaft_torque_nm: [-0.2, 0.03],
        auxiliary_load: None,
    };
    let before = dynamic_state(&prepare(p.clone(), gates, dt, 1), &p);
    let reference = rk4(&p, before, gates, command.shaft_torque_nm, dt, 4096);
    let fine_reference = rk4(&p, before, gates, command.shaft_torque_nm, dt, 8192);
    for k in 0..45 {
        close(reference[k], fine_reference[k], 2e-11);
    }
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let mut native = prepare(p.clone(), gates, dt, n);
        let mut next = [0.0; 15];
        let s = native
            .step_controlled(&before, command, &mut next, || false)
            .unwrap();
        let mut actual = [0.0; 45];
        actual[..15].copy_from_slice(&next);
        actual[15..33].copy_from_slice(&s.all_heat_joules());
        actual[33..].copy_from_slice(&receipt_work(&s));
        let error = core::array::from_fn::<_, 3, _>(|group| {
            let range = match group {
                0 => 0..15,
                1 => 15..33,
                _ => 33..45,
            };
            range
                .map(|k| (actual[k] - reference[k]).abs())
                .fold(0.0, f64::max)
        });
        errors.push(error);
        for k in 0..45 {
            close(actual[k], reference[k], 2e-5);
        }
    }
    for pair in errors.windows(2) {
        for group in 0..3 {
            let ratio = pair[0][group] / pair[1][group];
            assert!(
                ratio > 3.7 && ratio < 4.3,
                "state/heat/work group {group} did not refine quadratically: {errors:?}"
            );
        }
    }
}

#[test]
fn zero_current_cutoff_and_late_source_device_or_cancellation_refusals_preserve_charges() {
    let mut p = parameters();
    p.gate_source_capacitance_f = [0.01; 4];
    p.gate_drain_capacitance_f = [0.0; 4];
    let before = p.state(0.0, [2.0, 2.0], [1.0, 1.0], [0.0; 4]).unwrap();
    let dt = 0.002;
    let mut cutoff = PreparedGateCommutatedSwitchCycle::new_gate_controlled(
        p.clone(),
        [1.0; 2],
        &[interval([true; 4], dt)],
        1,
        electrical(),
        controls(),
        || false,
    )
    .unwrap();
    let mut next = [0.0; 9];
    let s = cutoff
        .step_controlled(
            &before,
            &[2.0; 2],
            &[None; 2],
            &mut next,
            &mut [0.0; 2],
            || false,
        )
        .unwrap();
    assert!(s.converter.channel_heat_j.iter().all(|q| *q == 0.0));
    assert!(s.converter.diode_heat_j.iter().all(|q| *q == 0.0));
    for k in 0..5 {
        close(next[k], before[k], 2e-12);
    }
    let physical = gate_physical(&p, &next);
    assert!(physical[5..].iter().all(|vg| *vg > 0.0 && *vg < 1.0));
    let deck = [interval([false; 4], 0.0001), interval([true; 4], 0.04)];
    for reason in 0..3 {
        let mut limited = p.clone();
        if reason == 0 {
            limited.channel_conductance[0] = GateChannelLaw::new(SaturatingChannel {
                gate_max: 1.1,
                active_above: None,
            })
            .unwrap();
        }
        if reason == 1 {
            limited.converter.devices[0]
                .bounds
                .maximum_channel_current_a = 0.05;
        }
        if reason == 2 {
            limited.channel_conductance[0] = GateChannelLaw::new(SaturatingChannel {
                gate_max: 4.0,
                active_above: Some(1.0),
            })
            .unwrap();
        }
        let mut native = PreparedGateCommutatedSwitchCycle::new_gate_controlled(
            limited,
            [1.0; 2],
            &deck,
            2,
            electrical(),
            controls(),
            || false,
        )
        .unwrap();
        let mut output = [73.0; 9];
        let mut bus = [73.0; 2];
        assert!(
            native
                .step_controlled(
                    &before,
                    &[2.0; 2],
                    &[None; 2],
                    &mut output,
                    &mut bus,
                    || false
                )
                .is_err(),
            "reason {reason}"
        );
        assert_eq!(output, [73.0; 9]);
        assert_eq!(bus, [73.0; 2]);
    }
    let p = parameters();
    let gates = [true, false, true, false];
    let dt = 1.0 / 1024.0;
    let mut native = prepare(p.clone(), gates, dt, 4);
    let before = dynamic_state(&native, &p);
    let command = BatteryMotorSwitchCommand {
        shaft_torque_nm: [-0.2, 0.03],
        auxiliary_load: None,
    };
    let polls = Cell::new(0usize);
    let mut accepted = [0.0; 15];
    native
        .step_controlled(&before, command, &mut accepted, || {
            polls.set(polls.get() + 1);
            false
        })
        .unwrap();
    let limit = polls.get() - 1;
    polls.set(0);
    let mut output = [73.0; 15];
    assert!(matches!(
        native.step_controlled(&before, command, &mut output, || {
            polls.set(polls.get() + 1);
            polls.get() >= limit
        }),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!(output, [73.0; 15]);
    native
        .step_controlled(&before, command, &mut output, || false)
        .unwrap();
    assert_eq!(output.map(f64::to_bits), accepted.map(f64::to_bits));
}
