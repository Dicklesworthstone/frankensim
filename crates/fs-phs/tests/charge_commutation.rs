//! G1/G3/G4: original nonlinear device Q(V), exact physical stored energy,
//! work-conjugate electrode-charge chords, and finite commutation dynamics.
//! The synthetic source has Cdiff=3c below |V|=1 and c above it. Its energy
//! above the knee is c*(1+V^2/2), rather than Cdiff(V)*V^2/2.

use fs_phs::circuit::battery::{BatteryOcvCurve, BatteryParameters, BatterySource};
use fs_phs::circuit::capacitance::CapacitorChargeCurve;
use fs_phs::circuit::converter::battery_fed::motor::BatteryMotorDrive;
use fs_phs::circuit::converter::commutation::capacitance::battery_motor::{
    CapacitiveBatteryMotorCycleStep, PreparedChargeBatteryMotorCycle,
};
use fs_phs::circuit::converter::commutation::capacitance::charge::ChargeCommutationParameters;
use fs_phs::circuit::converter::commutation::capacitance::{
    CapacitiveCommutatedStep, CapacitiveCommutationParameters,
    PreparedCapacitiveCommutatedSwitchCycle, PreparedChargeCommutatedSwitchCycle,
};
use fs_phs::circuit::converter::commutation::{
    BodyDiode, CommutationInterval, CommutationParameters, DeviceGateState, DeviceOperatingBounds,
    SwitchDevice,
};
use fs_phs::circuit::converter::power_load::{PowerLoadControls, RegulatedPowerLoad};
use fs_phs::circuit::converter::switching::SwitchCycleControls;
use fs_phs::circuit::converter::switching::battery_motor::BatteryMotorSwitchCommand;
use fs_phs::circuit::converter::switching::power_load::PowerSwitchCommand;
use fs_phs::circuit::converter::{AveragedBuckBoost, BridgeModulation, BuckBoostParameters};
use fs_phs::circuit::motor::{DcMotor, DcMotorParameters};
use fs_phs::{PreparedStepError, SolveBudget};
use std::cell::Cell;

const DEVICE: SwitchDevice = SwitchDevice {
    on_resistance_ohm: 0.2,
    off_resistance_ohm: 40.0,
    body_diode: BodyDiode {
        forward_voltage_v: 0.5,
        differential_resistance_ohm: 0.1,
    },
    bounds: DeviceOperatingBounds {
        maximum_reverse_voltage_v: 100.0,
        maximum_channel_current_a: 100.0,
        maximum_diode_current_a: 100.0,
    },
};
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; tolerance {tolerance:e}"
    );
}
fn curve(c: f64) -> CapacitorChargeCurve {
    CapacitorChargeCurve::try_new(vec![
        (-22.0 * c, -20.0),
        (-3.0 * c, -1.0),
        (0.0, 0.0),
        (3.0 * c, 1.0),
        (22.0 * c, 20.0),
    ])
    .unwrap()
}
fn parameters(bus: [f64; 2], cap: [f64; 4]) -> ChargeCommutationParameters {
    ChargeCommutationParameters {
        converter: CommutationParameters {
            inductance_h: 0.4,
            bus_capacitance_f: bus,
            winding_resistance_ohm: 0.1,
            devices: [DEVICE; 4],
        },
        output_charge_curves: cap.map(curve),
        maximum_cells_per_chord: 64,
    }
}
fn electrical() -> PowerLoadControls {
    PowerLoadControls {
        solve: SolveBudget {
            relative: 1e-13,
            absolute: 1e-13,
            max_iterations: 40,
        },
        backtracking_halvings: 10,
        absolute_energy_tolerance_j: 1e-10,
        relative_energy_tolerance: 1e-11,
    }
}
fn controls() -> SwitchCycleControls {
    SwitchCycleControls {
        absolute_energy_tolerance_j: 4e-9,
        relative_energy_tolerance: 2e-11,
        ..SwitchCycleControls::default()
    }
}
fn interval(on: [bool; 4], duration_s: f64) -> CommutationInterval {
    CommutationInterval {
        gates: DeviceGateState { on },
        duration_s,
    }
}
fn prepare(
    p: ChargeCommutationParameters,
    r: [f64; 2],
    deck: &[CommutationInterval],
    n: usize,
) -> PreparedChargeCommutatedSwitchCycle {
    PreparedChargeCommutatedSwitchCycle::new_charge_controlled(
        p,
        r,
        deck,
        n,
        electrical(),
        controls(),
        || false,
    )
    .unwrap()
}

// Independent source integrals in PHYSICAL voltage coordinates.
fn charge(c: f64, v: f64) -> f64 {
    if v.abs() <= 1.0 {
        3.0 * c * v
    } else {
        c * (v + 2.0 * v.signum())
    }
}
fn energy(c: f64, v: f64) -> f64 {
    if v.abs() <= 1.0 {
        1.5 * c * v * v
    } else {
        c * (1.0 + 0.5 * v * v)
    }
}
fn expected_state(
    p: &ChargeCommutationParameters,
    c: [f64; 4],
    i: f64,
    b: [f64; 2],
    n: [f64; 2],
) -> [f64; 5] {
    let h = [charge(c[0], b[0] - n[0]), charge(c[2], b[1] - n[1])];
    [
        p.converter.inductance_h * i,
        p.converter.bus_capacitance_f[0] * b[0] + h[0],
        p.converter.bus_capacitance_f[1] * b[1] + h[1],
        charge(c[1], n[0]) - h[0],
        charge(c[3], n[1]) - h[1],
    ]
}
#[derive(Clone, Copy)]
struct Segment {
    low: f64,
    high: f64,
    c: f64,
    offset: f64,
}
fn segments(c: f64) -> [Segment; 3] {
    [
        Segment {
            low: -20.0,
            high: -1.0,
            c,
            offset: -2.0 * c,
        },
        Segment {
            low: -1.0,
            high: 1.0,
            c: 3.0 * c,
            offset: 0.0,
        },
        Segment {
            low: 1.0,
            high: 20.0,
            c,
            offset: 2.0 * c,
        },
    ]
}
fn inverse(cb: f64, h: Segment, l: Segment) -> [f64; 3] {
    let det = cb * (h.c + l.c) + h.c * l.c;
    [(h.c + l.c) / det, h.c / det, (cb + h.c) / det]
}
fn in_cell(q: [f64; 2], h: Segment, l: Segment, inverse: [f64; 3]) -> [f64; 2] {
    let b = q[0] - h.offset;
    let n = q[1] - l.offset + h.offset;
    [
        inverse[0] * b + inverse[1] * n,
        inverse[1] * b + inverse[2] * n,
    ]
}
fn leg_voltage(cb: f64, ch: f64, cl: f64, q: [f64; 2]) -> [f64; 2] {
    // Enumerate the nine analytical voltage cells of this FIXED fixture.
    // Production instead traces original source knots in device charge.
    for high in segments(ch) {
        for low in segments(cl) {
            let v = in_cell(q, high, low, inverse(cb, high, low));
            let h = v[0] - v[1];
            if h >= high.low - 1e-12
                && h <= high.high + 1e-12
                && v[1] >= low.low - 1e-12
                && v[1] <= low.high + 1e-12
            {
                return v;
            }
        }
    }
    panic!("fixture electrode charge is outside its explicit voltage support");
}
fn physical(p: &ChargeCommutationParameters, c: [f64; 4], x: &[f64; 5]) -> [f64; 5] {
    let mut out = [0.0; 5];
    out[0] = x[0] / p.converter.inductance_h;
    for k in 0..2 {
        let v = leg_voltage(
            p.converter.bus_capacitance_f[k],
            c[2 * k],
            c[2 * k + 1],
            [x[1 + k], x[3 + k]],
        );
        out[1 + k] = v[0];
        out[3 + k] = v[1];
    }
    out
}
fn parts(p: &ChargeCommutationParameters, c: [f64; 4], x: &[f64; 5]) -> [f64; 7] {
    let v = physical(p, c, x);
    [
        0.5 * p.converter.inductance_h * v[0] * v[0],
        0.5 * p.converter.bus_capacitance_f[0] * v[1] * v[1],
        0.5 * p.converter.bus_capacitance_f[1] * v[2] * v[2],
        energy(c[0], v[1] - v[3]),
        energy(c[1], v[3]),
        energy(c[2], v[2] - v[4]),
        energy(c[3], v[4]),
    ]
}
fn heat(step: &CapacitiveCommutatedStep) -> [f64; 11] {
    let c = &step.converter;
    [
        c.winding_heat_j,
        c.channel_heat_j[0],
        c.channel_heat_j[1],
        c.channel_heat_j[2],
        c.channel_heat_j[3],
        c.diode_heat_j[0],
        c.diode_heat_j[1],
        c.diode_heat_j[2],
        c.diode_heat_j[3],
        step.terminal_resistor_heat_j[0],
        step.terminal_resistor_heat_j[1],
    ]
}
fn account(
    p: &ChargeCommutationParameters,
    c: [f64; 4],
    a: &[f64; 5],
    b: &[f64; 5],
    step: &CapacitiveCommutatedStep,
    scale: f64,
) {
    let start = parts(p, c, a);
    let end = parts(p, c, b);
    let delta = core::array::from_fn::<_, 7, _>(|i| end[i] - start[i]);
    close(
        (step.converter.inductor_energy_change_j - delta[0]) / scale,
        0.0,
        3e-9,
    );
    for i in 0..2 {
        close(
            (step.converter.bus_capacitor_energy_change_j[i] - delta[1 + i]) / scale,
            0.0,
            3e-9,
        );
        close(
            (step.source_work_j[i]
                - step.terminal_resistor_heat_j[i]
                - step.load_work_j[i]
                - step.converter.terminal_work_j[i])
                / scale,
            0.0,
            3e-9,
        );
    }
    for i in 0..4 {
        close(
            (step.converter.device_capacitor_energy_change_j[i] - delta[3 + i]) / scale,
            0.0,
            3e-9,
        );
    }
    let q = heat(step);
    assert!(q.iter().all(|v| v.is_finite() && *v >= 0.0));
    close(
        (step.converter.storage_change_j - delta.iter().sum::<f64>()) / scale,
        0.0,
        3e-9,
    );
    close(
        (step.native.delta_h - delta.iter().sum::<f64>()) / scale,
        0.0,
        3e-9,
    );
    close(
        (step.native.dissipated - q[..9].iter().sum::<f64>()) / scale,
        0.0,
        3e-9,
    );
    close(
        (delta.iter().sum::<f64>() + q.iter().sum::<f64>() + step.load_work_j.iter().sum::<f64>()
            - step.source_work_j.iter().sum::<f64>())
            / scale,
        0.0,
        3e-9,
    );
}

#[test]
fn linear_charge_laws_recover_constant_capacitor_states_and_all_original_energy_owners() {
    let c = [0.02, 0.03, 0.025, 0.015];
    let mut p = parameters([0.3, 0.4], c);
    p.output_charge_curves = c.map(|cap| {
        CapacitorChargeCurve::try_new(vec![(-20.0 * cap, -20.0), (0.0, 0.0), (20.0 * cap, 20.0)])
            .unwrap()
    });
    let constant = CapacitiveCommutationParameters {
        converter: p.converter,
        output_capacitance_f: c,
    };
    let dt = 1.0 / 1024.0;
    let deck = [
        interval([true, false, true, false], dt),
        interval([false; 4], dt),
        interval([true; 4], dt),
    ];
    let r = [0.5, 0.7];
    let before = constant.state(0.5, [7.0, 3.0], [2.0, 1.0]).unwrap();
    let authored = p.state(0.5, [7.0, 3.0], [2.0, 1.0]).unwrap();
    for (a, b) in before.into_iter().zip(authored) {
        close(a, b, 3e-14);
    }
    let mut original = PreparedCapacitiveCommutatedSwitchCycle::new_controlled(
        constant,
        r,
        &deck,
        3,
        electrical(),
        controls(),
        || false,
    )
    .unwrap();
    let mut nonlinear = prepare(p, r, &deck, 3);
    let loads = [
        None,
        Some(RegulatedPowerLoad::try_new(2.0, 1.0, 0.5, 30.0).unwrap()),
    ];
    let (mut a, mut b, mut va, mut vb) = ([0.0; 5], [0.0; 5], [0.0; 2], [0.0; 2]);
    let sa = original
        .step_controlled(&before, &[8.0, 3.0], &loads, &mut a, &mut va, || false)
        .unwrap();
    let sb = nonlinear
        .step_controlled(&before, &[8.0, 3.0], &loads, &mut b, &mut vb, || false)
        .unwrap();
    for (a, b) in a.into_iter().zip(b) {
        close(a, b, 2e-10);
    }
    for (a, b) in va.into_iter().zip(vb) {
        close(a, b, 2e-10);
    }
    for (a, b) in heat(&sa).into_iter().zip(heat(&sb)) {
        close(a, b, 2e-10);
    }
    for i in 0..2 {
        close(sa.source_work_j[i], sb.source_work_j[i], 2e-10);
        close(sa.load_work_j[i], sb.load_work_j[i], 2e-10);
        close(
            sa.converter.terminal_work_j[i],
            sb.converter.terminal_work_j[i],
            2e-10,
        );
        close(
            sa.converter.bus_capacitor_energy_change_j[i],
            sb.converter.bus_capacitor_energy_change_j[i],
            2e-10,
        );
    }
    for i in 0..4 {
        close(
            sa.converter.device_capacitor_energy_change_j[i],
            sb.converter.device_capacitor_energy_change_j[i],
            2e-10,
        );
    }
    close(sa.native.delta_h, sb.native.delta_h, 2e-10);
    close(sa.native.dissipated, sb.native.dissipated, 2e-10);
    close(sa.native.supplied, sb.native.supplied, 2e-10);
}

#[test]
fn true_nonlinear_rc_charge_and_resistor_heat_refine_across_a_source_knee_in_si_and_nanoseconds() {
    for scale in [1.0, 1e-9] {
        let c = 0.25 * scale;
        let r = 2.0;
        let mut p = parameters([c; 2], [c; 4]);
        for device in &mut p.converter.devices {
            device.on_resistance_ohm = r;
            device.off_resistance_ohm = r;
        }
        let crossing = r * c * 2.0_f64.ln();
        let deck = [
            interval([true, false, true, false], crossing),
            interval([false, true, false, true], r * c),
        ];
        let before = p.state(0.0, [4.0; 2], [2.0; 2]).unwrap();
        for (a, b) in before.into_iter().zip([0.0, 8.0 * c, 8.0 * c, 0.0, 0.0]) {
            close((a - b) / scale, 0.0, 3e-13);
        }
        // q_node=0 forces Vbus=2*Vnode. Above the knee, the exact node
        // voltage decays with rate1/(Rc); below it, the rate is3/(5Rc).
        // Authored intervals straddle the analytical crossing time, without
        // resetting any charge or otherwise forcing the numerical crossing.
        let v = (-0.6_f64).exp();
        let final_charge = [0.0, 5.0 * c * v, 5.0 * c * v, 0.0, 0.0];
        let release = 14.0 * c - 5.0 * c * v * v;
        let mut errors = Vec::new();
        for n in [8, 16, 32] {
            let mut solve = electrical();
            let mut cycle = controls();
            // Scale dimensionful tolerances with the C/time similarity.
            // Relative accuracy and the independent voltage/normalized-heat
            // assertions remain identical in the nanofarad/nanosecond case.
            solve.solve.absolute *= scale;
            solve.absolute_energy_tolerance_j *= scale;
            cycle.absolute_energy_tolerance_j *= scale;
            let mut native = PreparedChargeCommutatedSwitchCycle::new_charge_controlled(
                p.clone(),
                [r; 2],
                &deck,
                n,
                solve,
                cycle,
                || false,
            )
            .unwrap();
            let (mut next, mut mean) = ([0.0; 5], [0.0; 2]);
            let step = native
                .step_controlled(&before, &[0.0; 2], &[None; 2], &mut next, &mut mean, || {
                    false
                })
                .unwrap();
            account(&p, [c; 4], &before, &next, &step, scale);
            let actual = physical(&p, [c; 4], &next);
            let mut error: f64 = 0.0;
            for i in 0..5 {
                error = error.max((next[i] - final_charge[i]).abs() / scale);
            }
            for leg in 0..2 {
                error = error
                    .max((actual[1 + leg] - 2.0 * v).abs())
                    .max((actual[3 + leg] - v).abs());
                error = error.max(
                    (step.terminal_resistor_heat_j[leg] - (2.0 / 3.0) * release).abs() / scale,
                );
                close(step.source_work_j[leg], 0.0, 0.0);
            }
            for device in 0..4 {
                error = error
                    .max((step.converter.channel_heat_j[device] - release / 6.0).abs() / scale);
                error = error.max(
                    (step.converter.device_capacitor_energy_change_j[device]
                        - (1.5 * c * v * v - 3.0 * c))
                        .abs()
                        / scale,
                );
                assert_eq!(step.converter.diode_heat_j[device], 0.0);
            }
            assert!(actual[3] < 1.0 && actual[4] < 1.0);
            errors.push(error);
        }
        assert!(errors[2] < 4e-4, "scale{scale:e}: {errors:?}");
        for pair in errors.windows(2) {
            assert!(
                (3.0..5.0).contains(&(pair[0] / pair[1])),
                "scale{scale:e}: {errors:?}"
            );
        }
    }
}

#[test]
fn charge_history_gate_commands_and_bounded_refusals_are_transactional() {
    let c = [0.2, 0.1, 0.3, 0.2];
    let p = parameters([0.3, 0.4], c);
    let dt = 1.0 / 1024.0;
    let deck = [
        interval([true, false, true, false], dt),
        interval([false; 4], dt),
        interval([true; 4], dt),
    ];
    let commands = [
        PowerSwitchCommand {
            source_voltage_v: [8.0, 3.0],
            loads: [None; 2],
        },
        PowerSwitchCommand {
            source_voltage_v: [7.0, 4.0],
            loads: [
                None,
                Some(RegulatedPowerLoad::try_new(2.0, 1.0, 0.5, 30.0).unwrap()),
            ],
        },
        PowerSwitchCommand {
            source_voltage_v: [6.0, 2.0],
            loads: [None; 2],
        },
    ];
    let before = p.state(1.0, [7.0, 3.0], [2.0, 1.0]).unwrap();
    let mut native = prepare(p.clone(), [0.5, 0.7], &deck, 3);
    let (mut complete, mut means) = ([0.0; 5], [0.0; 2]);
    let polls = Cell::new(0usize);
    let good = native
        .step_interval_commands_controlled(&before, &commands, &mut complete, &mut means, || {
            polls.set(polls.get() + 1);
            false
        })
        .unwrap();
    account(&p, c, &before, &complete, &good, 1.0);
    let mut x = before;
    let mut heats = [0.0; 11];
    let mut work = [0.0; 2];
    for (edge, command) in deck.into_iter().zip(commands) {
        let mut split = prepare(p.clone(), [0.5, 0.7], &[edge], 3);
        let (mut next, mut voltage) = ([0.0; 5], [0.0; 2]);
        let step = split
            .step_controlled(
                &x,
                &command.source_voltage_v,
                &command.loads,
                &mut next,
                &mut voltage,
                || false,
            )
            .unwrap();
        for (total, value) in heats.iter_mut().zip(heat(&step)) {
            *total += value;
        }
        for i in 0..2 {
            work[i] += step.source_work_j[i];
        }
        x = next;
    }
    assert_eq!(x.map(f64::to_bits), complete.map(f64::to_bits));
    for (a, b) in heat(&good).into_iter().zip(heats) {
        close(a, b, 2e-12);
    }
    for i in 0..2 {
        close(good.source_work_j[i], work[i], 2e-12);
    }
    let stop = polls.get() - 1;
    polls.set(0);
    let (mut next, mut voltage) = ([77.0; 5], [88.0; 2]);
    assert!(matches!(
        native.step_interval_commands_controlled(
            &before,
            &commands,
            &mut next,
            &mut voltage,
            || {
                polls.set(polls.get() + 1);
                polls.get() >= stop
            }
        ),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!(next, [77.0; 5]);
    assert_eq!(voltage, [88.0; 2]);
    let mut bad = commands;
    bad[2].loads[1] = Some(RegulatedPowerLoad::try_new(1.0, 1.0, 0.25, 0.5).unwrap());
    assert!(
        native
            .step_interval_commands_controlled(&before, &bad, &mut next, &mut voltage, || false)
            .is_err()
    );
    assert_eq!(next, [77.0; 5]);
    assert_eq!(voltage, [88.0; 2]);
    let retry = native
        .step_interval_commands_controlled(&before, &commands, &mut next, &mut voltage, || false)
        .unwrap();
    assert_eq!(next.map(f64::to_bits), complete.map(f64::to_bits));
    assert_eq!(heat(&retry), heat(&good));
    assert!(p.state(0.0, [21.0, 4.0], [0.0, 2.0]).is_err());
    let mut limited = parameters([0.25; 2], [0.25; 4]);
    limited.maximum_cells_per_chord = 1;
    for d in &mut limited.converter.devices {
        d.on_resistance_ohm = 2.0;
        d.off_resistance_ohm = 2.0;
    }
    let start = limited.state(0.0, [4.0; 2], [2.0; 2]).unwrap();
    let mut limited = prepare(
        limited,
        [2.0; 2],
        &[interval([true; 4], 0.5 * (2.0_f64.ln() + 1.0))],
        1,
    );
    next = [77.0; 5];
    voltage = [88.0; 2];
    assert!(
        limited
            .step_controlled(
                &start,
                &[0.0; 2],
                &[None; 2],
                &mut next,
                &mut voltage,
                || false
            )
            .is_err()
    );
    assert_eq!(next, [77.0; 5]);
    assert_eq!(voltage, [88.0; 2]);
}

const DRIVE_CAP: [f64; 4] = [0.02, 0.03, 0.025, 0.015];
const BATTERY: BatteryParameters = BatteryParameters {
    series_resistance_ohm: 0.1,
    polarization_resistance_ohm: 0.2,
    polarization_capacitance_f: 0.5,
    bus_capacitance_f: 0.03,
};
const DRIVE_CONVERTER: BuckBoostParameters = BuckBoostParameters {
    inductance_h: 0.04,
    bus_capacitance_f: [0.04, 0.05],
    winding_resistance_ohm: 0.07,
    input_on_resistance_ohm: [0.031, 0.047],
    output_on_resistance_ohm: [0.053, 0.071],
};
const MOTOR: DcMotorParameters = DcMotorParameters {
    inductance_h: 0.02,
    winding_resistance_ohm: 0.2,
    rotor_inertia_kg_m2: 0.08,
    stator_inertia_kg_m2: 0.4,
    coupling_nm_per_a: 0.4,
    bearing_drag_nm_s_per_rad: 0.01,
};
const LEAD: f64 = 0.08;
fn drive_model() -> BatteryMotorDrive {
    AveragedBuckBoost::new(
        DRIVE_CONVERTER,
        BridgeModulation {
            input_high_fraction: 0.2,
            output_high_fraction: 0.8,
        },
    )
    .unwrap()
    .with_battery_source(
        BatterySource::new(
            BatteryOcvCurve::new(80.0, &[(0.0, 12.0), (1.0, 12.0)]).unwrap(),
            BATTERY,
        )
        .unwrap(),
        LEAD,
    )
    .unwrap()
    .with_dc_motor(DcMotor::new(MOTOR).unwrap())
    .unwrap()
}
fn drive_parameters() -> ChargeCommutationParameters {
    let mut p = parameters(DRIVE_CONVERTER.bus_capacitance_f, DRIVE_CAP);
    p.converter.inductance_h = DRIVE_CONVERTER.inductance_h;
    p.converter.winding_resistance_ohm = DRIVE_CONVERTER.winding_resistance_ohm;
    p
}
fn drive_prepare(gates: [bool; 4], dt: f64, n: usize) -> PreparedChargeBatteryMotorCycle {
    PreparedChargeBatteryMotorCycle::new_charge_controlled(
        drive_model(),
        [DEVICE; 4],
        DRIVE_CAP.map(curve),
        64,
        &[interval(gates, dt)],
        n,
        electrical(),
        controls(),
        || false,
    )
    .unwrap()
}
fn drive_state(native: &PreparedChargeBatteryMotorCycle, y: [f64; 11]) -> [f64; 11] {
    let p = drive_parameters();
    let c = p.state(y[3], [y[4], y[5]], [y[6], y[7]]).unwrap();
    let authored = expected_state(&p, DRIVE_CAP, y[3], [y[4], y[5]], [y[6], y[7]]);
    for (a, b) in c.into_iter().zip(authored) {
        close(a, b, 2e-13);
    }
    native
        .state(
            &[
                y[0],
                BATTERY.polarization_capacitance_f * y[1],
                BATTERY.bus_capacitance_f * y[2],
            ],
            &c,
            &[
                MOTOR.inductance_h * y[8],
                MOTOR.rotor_inertia_kg_m2 * y[9],
                MOTOR.stator_inertia_kg_m2 * y[10],
            ],
        )
        .unwrap()
}
fn drive_parts(x: &[f64; 11]) -> [f64; 7] {
    let c = core::array::from_fn(|i| x[3 + i]);
    [
        12.0 * x[0],
        0.5 * x[1] * x[1] / BATTERY.polarization_capacitance_f,
        0.5 * x[2] * x[2] / BATTERY.bus_capacitance_f,
        parts(&drive_parameters(), DRIVE_CAP, &c).iter().sum(),
        0.5 * x[8] * x[8] / MOTOR.inductance_h,
        0.5 * x[9] * x[9] / MOTOR.rotor_inertia_kg_m2,
        0.5 * x[10] * x[10] / MOTOR.stator_inertia_kg_m2,
    ]
}
fn drive_account(a: &[f64; 11], b: &[f64; 11], step: &CapacitiveBatteryMotorCycleStep) {
    let ea = drive_parts(a);
    let eb = drive_parts(b);
    let change = core::array::from_fn::<_, 7, _>(|i| eb[i] - ea[i]);
    let q = step.heat_joules();
    let w = [
        step.battery_terminal_work_j,
        step.converter.terminal_work_j[0],
        step.converter.terminal_work_j[1],
        step.motor_electrical_work_j,
        step.electromagnetic_work_j,
        step.shaft_work_j[0],
        step.shaft_work_j[1],
        step.load_work_j,
    ];
    for i in 0..7 {
        close(step.storage_change_j[i], change[i], 3e-9);
    }
    let ca = parts(
        &drive_parameters(),
        DRIVE_CAP,
        &core::array::from_fn(|i| a[3 + i]),
    );
    let cb = parts(
        &drive_parameters(),
        DRIVE_CAP,
        &core::array::from_fn(|i| b[3 + i]),
    );
    for i in 0..4 {
        close(
            step.converter.device_capacitor_energy_change_j[i],
            cb[3 + i] - ca[3 + i],
            3e-9,
        );
    }
    assert!(q.iter().all(|v| v.is_finite() && *v >= 0.0));
    close(change[..3].iter().sum::<f64>() + q[5] + q[6], w[0], 3e-9);
    close(w[0] + q[7] + w[1], 0.0, 3e-9);
    close(
        change[3] + q[..5].iter().sum::<f64>() + q[10..].iter().sum::<f64>(),
        w[1] + w[2],
        3e-9,
    );
    close(w[2] + w[3] + w[7], 0.0, 3e-9);
    close(change[4] + q[8] + w[4], w[3], 3e-9);
    close(change[5] + change[6] + q[9], w[4] + w[5] + w[6], 3e-9);
    close(step.native.delta_h, change.iter().sum(), 3e-9);
    close(step.native.dissipated, q.iter().sum(), 3e-9);
    close(step.native.supplied, w[5] + w[6] - w[7], 3e-9);
    close(step.balance_residual_j(), 0.0, 3e-9);
    close(step.angular_momentum_residual_nms, 0.0, 3e-11);
}

#[test]
fn finite_battery_regeneration_retains_original_nonlinear_capacitor_charge_and_both_shafts() {
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
    let mut native = drive_prepare([true, false, true, false], dt, 4);
    let before = drive_state(
        &native,
        [40.0, pol, vb, ic, vi, vo, left, right, im, wr, ws],
    );
    let mut next = [0.0; 11];
    let step = native
        .step_controlled(
            &before,
            BatteryMotorSwitchCommand {
                shaft_torque_nm: [reaction, -reaction],
                auxiliary_load: None,
            },
            &mut next,
            || false,
        )
        .unwrap();
    close(next[0], before[0] - cell * dt, 2e-11);
    for i in 1..11 {
        close(next[i], before[i], 2e-11);
    }
    drive_account(&before, &next, &step);
    close(step.storage_change_j[0], -12.0 * cell * dt, 2e-10);
    for q in step.converter.device_capacitor_energy_change_j {
        close(q, 0.0, 2e-12);
    }
    // All branch efforts and currents are exactly steady in this manufactured
    // trajectory. Only here can a cycle-mean current legitimately reproduce I²R.
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
    ];
    for (actual, power) in step.heat_joules().into_iter().zip(power) {
        close(actual, dt * power, 2e-10);
    }
    close(step.battery_terminal_work_j, -vb * cell * dt, 2e-10);
    close(step.converter.terminal_work_j[0], vi * cell * dt, 2e-10);
    close(step.converter.terminal_work_j[1], -vo * im * dt, 2e-10);
    close(step.motor_electrical_work_j, vo * im * dt, 2e-10);
    close(
        step.electromagnetic_work_j,
        MOTOR.coupling_nm_per_a * relative * im * dt,
        2e-10,
    );
    close(step.shaft_work_j[0], reaction * wr * dt, 2e-10);
    close(step.shaft_work_j[1], -reaction * ws * dt, 2e-10);
    assert!(step.storage_change_j[0] > 0.0 && step.motor_electrical_work_j < 0.0);
    assert!(step.shaft_work_j.iter().all(|v| *v > 0.0));
}

#[test]
fn eleven_state_crossing_and_dead_time_use_the_same_avf_efforts_for_kcl_and_fourteen_heats() {
    let dt = 1.0 / 512.0;
    for on in [true, false] {
        let gates = if on {
            [true, false, true, false]
        } else {
            [false; 4]
        };
        let mut native = drive_prepare(gates, dt, 1);
        let nodes = if on { [10.75, 4.0] } else { [-1.1, 7.5] };
        let before = drive_state(
            &native,
            [
                52.0, 0.1, 12.1, 0.5, 11.8, 6.4, nodes[0], nodes[1], 0.1, 4.0, -0.4,
            ],
        );
        let command = BatteryMotorSwitchCommand {
            shaft_torque_nm: [-0.2, 0.03],
            auxiliary_load: None,
        };
        let mut next = [0.0; 11];
        let step = native
            .step_controlled(&before, command, &mut next, || false)
            .unwrap();
        drive_account(&before, &next, &step);
        let c = &step.converter;
        let ic = c.mean_current_a;
        let vi = c.mean_bus_voltage_v[0];
        let vo = c.mean_bus_voltage_v[1];
        let left = c.mean_node_voltage_v[0];
        let right = c.mean_node_voltage_v[1];
        let vb = step.mean_bus_voltage_v[0];
        let lead = step.mean_lead_current_a;
        let pol = f64::midpoint(before[1], next[1]) / BATTERY.polarization_capacitance_f;
        let cell = (12.0 - pol - vb) / BATTERY.series_resistance_ohm;
        let im = step.mean_motor_current_a;
        let relative = step.mean_shaft_speed_rad_s[0] - step.mean_shaft_speed_rad_s[1];
        let vds = [vi - left, left, vo - right, right];
        let channel = core::array::from_fn::<_, 4, _>(|i| {
            vds[i]
                / if gates[i] {
                    DEVICE.on_resistance_ohm
                } else {
                    DEVICE.off_resistance_ohm
                }
        });
        let diode = vds.map(|v| {
            ((-v - DEVICE.body_diode.forward_voltage_v)
                / DEVICE.body_diode.differential_resistance_ohm)
                .max(0.0)
        });
        let current = core::array::from_fn::<_, 4, _>(|i| channel[i] - diode[i]);
        let torque = MOTOR.coupling_nm_per_a * im - MOTOR.bearing_drag_nm_s_per_rad * relative;
        let rates = [
            -cell,
            cell - pol / BATTERY.polarization_resistance_ohm,
            cell - lead,
            left - right - DRIVE_CONVERTER.winding_resistance_ohm * ic,
            lead - current[0],
            -im - current[2],
            current[0] - current[1] - ic,
            current[2] - current[3] + ic,
            vo - MOTOR.winding_resistance_ohm * im - MOTOR.coupling_nm_per_a * relative,
            torque + command.shaft_torque_nm[0],
            -torque + command.shaft_torque_nm[1],
        ];
        for i in 0..11 {
            close(next[i] - before[i], dt * rates[i], 3e-10);
        }
        close(lead, (vb - vi) / LEAD, 3e-12);
        let mut power = [0.0; 14];
        power[0] = DRIVE_CONVERTER.winding_resistance_ohm * ic * ic;
        for i in 0..4 {
            power[1 + i] = vds[i] * channel[i];
            power[10 + i] = -vds[i] * diode[i];
        }
        power[5] = BATTERY.series_resistance_ohm * cell * cell;
        power[6] = pol * pol / BATTERY.polarization_resistance_ohm;
        power[7] = LEAD * lead * lead;
        power[8] = MOTOR.winding_resistance_ohm * im * im;
        power[9] = MOTOR.bearing_drag_nm_s_per_rad * relative * relative;
        for (heat, power) in step.heat_joules().into_iter().zip(power) {
            close(heat, dt * power, 2e-10);
        }
        assert!(
            c.device_capacitor_energy_change_j
                .iter()
                .map(|x| x.abs())
                .sum::<f64>()
                > 1e-5
        );
        let endpoint = physical(
            &drive_parameters(),
            DRIVE_CAP,
            &core::array::from_fn(|i| next[3 + i]),
        );
        if on {
            assert!(
                endpoint[1] - endpoint[3] < 1.0,
                "input-high original source knee must be crossed"
            );
        } else {
            assert!(step.heat_joules()[11] > 0.0 && step.heat_joules()[12] > 0.0);
        }
        let sentinel = [77.0; 11];
        let mut trial = sentinel;
        let calls = Cell::new(0usize);
        native
            .step_controlled(&before, command, &mut trial, || {
                calls.set(calls.get() + 1);
                false
            })
            .unwrap();
        let limit = calls.get() - 1;
        calls.set(0);
        trial = sentinel;
        assert!(matches!(
            native.step_controlled(&before, command, &mut trial, || {
                calls.set(calls.get() + 1);
                calls.get() >= limit
            }),
            Err(PreparedStepError::Cancelled)
        ));
        assert_eq!(trial, sentinel);
    }
}
