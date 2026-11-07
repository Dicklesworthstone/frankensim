//! Finite supply storage, original material feedback and seven heat destinations.
use super::*;
use fs_phs::circuit::converter::capacitor_terminals::CapacitorTerminals;

fn capacitor_budget() -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: 5,
        ..budget()
    }
}
fn observations() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 10,
        max_dense_bytes: 1_000_000,
        max_multiply_terms: 100_000_000,
    }
}

#[test]
fn finite_capacitor_charge_sharing_heats_shared_resistance_with_actual_temperature_feedback() {
    let alpha = 0.0039;
    let cp = [1.0, 1.0, 1.0, 1.0, 1.0, 2.0];
    let (cards, curves) = sources([R0[0], R0[1], R0[2], R0[3], R0[4], 0.5], cp, alpha);
    let context = point(300.0);
    let bindings = bindings(&cards, &curves, &context);
    let resistance = key();
    let mut receivers = ALL_IDS;
    receivers[6] = receivers[5];
    let channels = paths(receivers, &resistance);
    let network = VoltageDrivenConverterThermalNetwork {
        cells: &bindings,
        links: &[],
        channels: &channels,
    };
    let thermal = body(&curves, cp);
    let terminals = CapacitorTerminals {
        capacitance_f: [0.4, 0.6],
        resistance_ohm: [0.5; 2],
    };
    let mode = modulation(0.0, 0.0);
    let native = AveragedBuckBoost::new(parameters(), mode)
        .unwrap()
        .with_capacitor_terminals(terminals)
        .unwrap();
    assert_eq!(native.modulation(), mode);
    let dt = 0.05;
    let intervals = [SwitchInterval {
        state: BridgeSwitchState {
            input_high: false,
            output_high: false,
        },
        duration_s: dt,
    }];
    for (bus, capacitor) in [([0.0, 4.0], [6.0, 2.0]), ([6.0, 4.0], [2.0, 1.0])] {
        let initial = native.state(0.0, bus, capacitor).unwrap();
        let converter = network
            .bind_capacitor(&thermal, &native, initial, || false)
            .unwrap();
        // Independent isolated two-capacitor RC law, with only the shared
        // resistor body's final temperature determining each frozen R.
        let heat_at = |theta: f64| -> [f64; 2] {
            core::array::from_fn(|p| {
                let cb = parameters().bus_capacitance_f[p];
                let cs = terminals.capacitance_f[p];
                let equivalent = cb * cs / (cb + cs);
                let difference = capacitor[p] - bus[p];
                let r = 0.5 * (1.0 + alpha * theta);
                0.5 * equivalent
                    * difference
                    * difference
                    * -(-2.0 * dt / (r * equivalent)).exp_m1()
            })
        };
        let (mut low, mut high) = (0.0, 10.0);
        for _ in 0..80 {
            let theta = 0.5 * (low + high);
            if cp[5] * theta > heat_at(theta).iter().sum::<f64>() {
                high = theta;
            } else {
                low = theta;
            }
        }
        let theta = 0.5 * (low + high);
        let expected_heat = heat_at(theta);
        for drive in [
            ConverterElectricalDrive::Averaged(mode),
            ConverterElectricalDrive::Switching {
                intervals: &intervals,
                controls: SwitchCycleControls::default(),
            },
        ] {
            let result = network
                .step_capacitor_candidate(
                    &thermal,
                    dt,
                    &converter,
                    drive,
                    &[0.0; 2],
                    capacitor_budget(),
                    observations(),
                    controls(),
                    &[0.0; 6],
                    || false,
                )
                .unwrap();
            let step = &result.report().electrical;
            for id in &ALL_IDS[..5] {
                near(result.thermal().temperature_k(*id).unwrap(), 300.0);
            }
            near(
                result.thermal().temperature_k(ALL_IDS[5]).unwrap(),
                300.0 + theta,
            );
            near(result.converter().state()[0], 0.0);
            assert_eq!(step.external_work_j, [0.0; 2]);
            assert_eq!(
                result.converter().terminals().capacitance_f,
                terminals.capacitance_f
            );
            for p in 0..2 {
                let cb = parameters().bus_capacitance_f[p];
                let cs = terminals.capacitance_f[p];
                let equivalent = cb * cs / (cb + cs);
                let r = 0.5 * (1.0 + alpha * theta);
                let transfer =
                    equivalent * (capacitor[p] - bus[p]) * -(-dt / (r * equivalent)).exp_m1();
                let bus_work = bus[p] * transfer + transfer * transfer / (2.0 * cb);
                let cap_outflow = capacitor[p] * transfer - transfer * transfer / (2.0 * cs);
                near(result.converter().state()[p + 1], cb * bus[p] + transfer);
                near(
                    result.converter().state()[p + 3],
                    cs * capacitor[p] - transfer,
                );
                near(result.converter().terminals().resistance_ohm[p], r);
                near(step.terminal_resistor_heat_j[p], expected_heat[p]);
                near(step.converter.terminal_work_j[p], bus_work);
                near(step.terminal_outflow_work_j[p], cap_outflow);
                near(step.terminal_storage_change_j[p], -cap_outflow);
                near(
                    result.converter().state()[p + 1] + result.converter().state()[p + 3],
                    initial[p + 1] + initial[p + 3],
                );
            }
            assert_eq!(step.terminal_outflow_work_j[0] > 0.0, capacitor[0] > bus[0]);
            for heat in step.converter.conduction.joules {
                near(heat, 0.0);
            }
            near(step.total.delta_h, -cp[5] * theta);
            near(step.total.dissipated, cp[5] * theta);
            near(result.report().deposited_joule_j, cp[5] * theta);
            near(
                result.thermal().energy_j() - thermal.energy_j(),
                cp[5] * theta,
            );
            near(result.report().total_energy_residual_j, 0.0);
        }
        assert_eq!(converter.state(), &initial);
    }
    assert_eq!(thermal.time_s(), 0.0);
}

fn rk4<const N: usize>(x: &mut [f64; N], h: f64, rhs: impl Fn(&[f64; N]) -> [f64; N]) {
    let k1 = rhs(x);
    let k2 = rhs(&core::array::from_fn(|i| x[i] + 0.5 * h * k1[i]));
    let k3 = rhs(&core::array::from_fn(|i| x[i] + 0.5 * h * k2[i]));
    let k4 = rhs(&core::array::from_fn(|i| x[i] + h * k3[i]));
    for i in 0..N {
        x[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
}

#[test]
fn finite_supply_pwm_uses_hot_materials_and_publishes_all_histories_atomically() {
    let r0 = [R0[0], R0[1], R0[2], R0[3], R0[4], 0.25, 0.8];
    let cp = [0.05, 0.05, 0.05, 0.05, 0.05, 0.2, 0.4];
    let alpha = 0.0039;
    let (cards, curves) = sources(r0, cp, alpha);
    let context = point(300.0);
    let bindings = bindings(&cards, &curves, &context);
    let resistance = key();
    let channels = paths(ALL_IDS, &resistance);
    let network = VoltageDrivenConverterThermalNetwork {
        cells: &bindings,
        links: &[],
        channels: &channels,
    };
    let thermal = body(&curves, cp);
    let original_cells = thermal.cells().clone();
    let dt = 0.01;
    let mode = modulation(0.3, 0.7);
    let intervals = centered_pwm(dt, mode).unwrap();
    let drive = ConverterElectricalDrive::Switching {
        intervals: &intervals,
        controls: SwitchCycleControls::default(),
    };
    let terminal = CapacitorTerminals {
        capacitance_f: [0.4, 0.6],
        resistance_ohm: [r0[5], r0[6]],
    };
    let native = AveragedBuckBoost::new(parameters(), mode)
        .unwrap()
        .with_capacitor_terminals(terminal)
        .unwrap();
    let initial = native.state(1.5, [24.0, 12.0], [26.0, 5.0]).unwrap();
    let converter = network
        .bind_capacitor(&thermal, &native, initial, || false)
        .unwrap();
    let source = [0.4, -0.2];
    let calls = Cell::new(0usize);
    let result = network
        .step_capacitor_candidate(
            &thermal,
            dt,
            &converter,
            drive,
            &source,
            capacitor_budget(),
            observations(),
            controls(),
            &[0.0; 7],
            || {
                calls.set(calls.get() + 1);
                false
            },
        )
        .unwrap();
    let step = &result.report().electrical;
    let heat: [f64; 7] = core::array::from_fn(|i| {
        if i < 5 {
            step.converter.conduction.joules[i]
        } else {
            step.terminal_resistor_heat_j[i - 5]
        }
    });
    for i in 0..7 {
        assert!(heat[i] > 0.0);
        let temperature = result.thermal().temperature_k(ALL_IDS[i]).unwrap();
        near(temperature, 300.0 + heat[i] / cp[i]);
        near(
            result.converter().channels()[i].resistance_ohm(),
            r0[i] * (1.0 + alpha * (temperature - 300.0)),
        );
        assert_eq!(
            converter.channels()[i].material().unwrap().card_identity(),
            result.converter().channels()[i]
                .material()
                .unwrap()
                .card_identity()
        );
    }
    assert!(
        result.converter().state()[3] < initial[3],
        "input supply discharges"
    );
    assert!(
        result.converter().state()[4] > initial[4],
        "receiving capacitor charges"
    );
    assert!(step.external_work_j[0] > 0.0 && step.external_work_j[1] < 0.0);
    let total_heat: f64 = heat.iter().sum();
    near(result.report().deposited_joule_j, total_heat);
    near(result.thermal().energy_j() - thermal.energy_j(), total_heat);
    let storage = |x: &[f64; 5]| {
        0.5 * (x[0] * x[0] / parameters().inductance_h
            + x[1] * x[1] / parameters().bus_capacitance_f[0]
            + x[2] * x[2] / parameters().bus_capacitance_f[1]
            + x[3] * x[3] / terminal.capacitance_f[0]
            + x[4] * x[4] / terminal.capacitance_f[1])
    };
    near(
        step.total.delta_h,
        storage(result.converter().state()) - storage(&initial),
    );
    near(
        step.total.delta_h + total_heat,
        step.external_work_j.iter().sum(),
    );
    near(result.report().total_energy_residual_j, 0.0);
    assert_eq!(result.report().circuit_trials % intervals.len(), 0);

    // Independent current/voltage ODE at actual accepted material R(T), with
    // its seven original heat integrands and signed charging-current powers.
    let p = result.converter().parameters();
    let t = result.converter().terminals();
    let mut oracle = [0.0; 14];
    oracle[..5].copy_from_slice(&[1.5, 24.0, 12.0, 26.0, 5.0]);
    for gate in &intervals {
        let a = if gate.state.input_high { 1.0 } else { 0.0 };
        let b = if gate.state.output_high { 1.0 } else { 0.0 };
        let r = [
            p.winding_resistance_ohm,
            a * p.input_on_resistance_ohm[0],
            (1.0 - a) * p.input_on_resistance_ohm[1],
            b * p.output_on_resistance_ohm[0],
            (1.0 - b) * p.output_on_resistance_ohm[1],
        ];
        for _ in 0..256 {
            rk4(&mut oracle, gate.duration_s / 256.0, |x| {
                let flows = [
                    (x[3] - x[1]) / t.resistance_ohm[0],
                    (x[4] - x[2]) / t.resistance_ohm[1],
                ];
                let mut rhs = [0.0; 14];
                rhs[0] = (a * x[1] - b * x[2] - r.iter().sum::<f64>() * x[0]) / p.inductance_h;
                rhs[1] = (flows[0] - a * x[0]) / p.bus_capacitance_f[0];
                rhs[2] = (flows[1] + b * x[0]) / p.bus_capacitance_f[1];
                rhs[3] = (source[0] - flows[0]) / t.capacitance_f[0];
                rhs[4] = (source[1] - flows[1]) / t.capacitance_f[1];
                for k in 0..5 {
                    rhs[5 + k] = r[k] * x[0] * x[0];
                }
                for k in 0..2 {
                    rhs[10 + k] = t.resistance_ohm[k] * flows[k] * flows[k];
                    rhs[12 + k] = source[k] * x[3 + k];
                }
                rhs
            });
        }
    }
    for (i, coefficient) in [
        p.inductance_h,
        p.bus_capacitance_f[0],
        p.bus_capacitance_f[1],
        t.capacitance_f[0],
        t.capacitance_f[1],
    ]
    .into_iter()
    .enumerate()
    {
        near(result.converter().state()[i], coefficient * oracle[i]);
    }
    for i in 0..7 {
        near(heat[i], oracle[5 + i]);
    }
    for i in 0..2 {
        near(step.external_work_j[i], oracle[12 + i]);
    }

    let last = calls.get();
    calls.set(0);
    assert!(matches!(
        network.step_capacitor_candidate(
            &thermal,
            dt,
            &converter,
            drive,
            &source,
            capacitor_budget(),
            observations(),
            controls(),
            &[0.0; 7],
            || {
                calls.set(calls.get() + 1);
                calls.get() >= last
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    let mut limited = controls();
    limited.maximum_circuit_trials = intervals.len() - 1;
    assert!(matches!(
        network.step_capacitor_candidate(
            &thermal,
            dt,
            &converter,
            drive,
            &source,
            capacitor_budget(),
            observations(),
            limited,
            &[0.0; 7],
            || false
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    for (target, current, energy) in [
        (2.0 * dt, source, [0.0; 7]),
        (dt, [f64::NAN, 0.0], [0.0; 7]),
        (dt, source, [10000.0; 7]),
    ] {
        assert!(
            network
                .step_capacitor_candidate(
                    &thermal,
                    target,
                    &converter,
                    drive,
                    &current,
                    capacitor_budget(),
                    observations(),
                    controls(),
                    &energy,
                    || false
                )
                .is_err()
        );
    }
    let mut changed = paths(ALL_IDS, &resistance);
    changed[6].cell = ALL_IDS[5];
    assert!(
        VoltageDrivenConverterThermalNetwork {
            cells: &bindings,
            links: &[],
            channels: &changed
        }
        .step_capacitor_candidate(
            &thermal,
            dt,
            &converter,
            drive,
            &source,
            capacitor_budget(),
            observations(),
            controls(),
            &[0.0; 7],
            || false
        )
        .is_err()
    );
    assert_eq!(thermal.cells(), &original_cells);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(converter.state(), &initial);
    let retry = network
        .step_capacitor_candidate(
            &thermal,
            dt,
            &converter,
            drive,
            &source,
            capacitor_budget(),
            observations(),
            controls(),
            &[0.0; 7],
            || false,
        )
        .unwrap();
    assert_eq!(
        (*retry.converter().state()).map(f64::to_bits),
        (*result.converter().state()).map(f64::to_bits)
    );

    let (advanced_thermal, advanced_converter, _) = result.into_parts();
    let second = network
        .step_capacitor_candidate(
            &advanced_thermal,
            2.0 * dt,
            &advanced_converter,
            drive,
            &source,
            capacitor_budget(),
            observations(),
            controls(),
            &[0.0; 7],
            || false,
        )
        .unwrap();
    near(
        second.report().electrical.total.delta_h,
        storage(second.converter().state()) - storage(advanced_converter.state()),
    );
    near(second.report().total_energy_residual_j, 0.0);
    assert_eq!(second.converter().time_s(), 2.0 * dt);
    assert!(
        network
            .step_capacitor_candidate(
                &advanced_thermal,
                2.0 * dt,
                &converter,
                drive,
                &source,
                capacitor_budget(),
                observations(),
                controls(),
                &[0.0; 7],
                || false
            )
            .is_err()
    );
}
