//! Independent RC/enthalpy balance plus the original loaded PWM trajectory.
#[path = "voltage/capacitor.rs"]
mod capacitor;
#[path = "voltage/power_load.rs"]
mod power_load;
use super::*;
use fs_couple::circuit_enthalpy::converter::voltage::*;
use fs_phs::circuit::converter::switching::voltage_driven::PreparedVoltageDrivenSwitchCycle;
use fs_phs::linear_exact::DenseExponentialBudget;

const ALL_IDS: [CellId; 7] = [
    IDS[0],
    IDS[1],
    IDS[2],
    IDS[3],
    IDS[4],
    CellId(106),
    CellId(117),
];
fn heat_budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 6,
        max_dense_bytes: 16_384,
        max_multiply_terms: 500_000,
    }
}
fn sources<const N: usize>(
    r: [f64; N],
    cp: [f64; N],
    alpha: f64,
) -> ([MaterialCard; N], [ResolvedHeatCapacityCurve; N]) {
    let pairs: [_; N] = std::array::from_fn(|i| material(r[i], cp[i], alpha));
    (
        pairs.each_ref().map(|(card, _)| card.clone()),
        pairs.map(|(_, curve)| curve),
    )
}
fn body<'a, const N: usize>(
    curves: &'a [ResolvedHeatCapacityCurve; N],
    cp: [f64; N],
) -> DepositionSimulation<'a> {
    let events = (0..N)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: ALL_IDS[i],
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: 1.0,
                specific_h_j_kg: 100.0 * cp[i],
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    let laws = curves
        .iter()
        .enumerate()
        .map(|(i, curve)| (MaterialId(i as u64), curve.curve() as &dyn ThermalLaw))
        .collect::<Vec<_>>();
    DepositionSimulation::new(
        events,
        &laws,
        DepositionBudget {
            max_cells: 8,
            max_links: 8,
            max_events: 16,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn bindings<'a, const N: usize>(
    cards: &'a [MaterialCard; N],
    curves: &'a [ResolvedHeatCapacityCurve; N],
    context: &'a QueryPoint,
) -> Vec<MaterialCell<'a>> {
    (0..N)
        .map(|i| MaterialCell {
            cell: ALL_IDS[i],
            material_id: MaterialId(i as u64),
            curve: curves[i].curve(),
            point: ThermalMaterialPoint {
                material: &cards[i],
                context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        })
        .collect()
}
fn paths(ids: [CellId; 7], key: &PropertyKey) -> [ConverterThermalChannel<'_>; 7] {
    std::array::from_fn(|i| ConverterThermalChannel {
        cell: ids[i],
        resistance: if i == 0 {
            ConverterResistance::Bulk {
                basis: ElectricalPropertyBasis::Resistivity,
                selection: MaterialPropertySelection::SingleClaimOnly,
                length_m: 1.0,
                area_m2: 1.0,
            }
        } else {
            ConverterResistance::Resistance {
                key,
                selection: MaterialPropertySelection::SingleClaimOnly,
            }
        },
    })
}

#[test]
fn terminal_resistor_feedback_and_shared_heat_follow_independent_rc_balances() {
    let alpha = 0.0039;
    let cp = [1.0, 1.0, 1.0, 1.0, 1.0, 2.0];
    let (cards, curves) = sources([R0[0], R0[1], R0[2], R0[3], R0[4], 0.5], cp, alpha);
    let context = point(300.0);
    let bindings = bindings(&cards, &curves, &context);
    let resistance = key();
    let mut receivers = ALL_IDS;
    receivers[6] = receivers[5]; // Both terminal resistors explicitly share one body.
    let channels = paths(receivers, &resistance);
    let network = VoltageDrivenConverterThermalNetwork {
        cells: &bindings,
        links: &[],
        channels: &channels,
    };
    let thermal = body(&curves, cp);
    let native = AveragedBuckBoost::new(parameters(), modulation(0.0, 0.0)).unwrap();
    let dt = 0.05;
    let intervals = [SwitchInterval {
        state: BridgeSwitchState {
            input_high: false,
            output_high: false,
        },
        duration_s: dt,
    }];
    for (source, initial_bus) in [([6.0, 0.0], [0.0, 4.0]), ([2.0, 0.0], [6.0, 4.0])] {
        let initial = native.state(0.0, initial_bus).unwrap();
        let converter = network
            .bind(&thermal, &native, [0.5; 2], initial, || false)
            .unwrap();
        let heat_at = |theta: f64| -> [f64; 2] {
            std::array::from_fn(|port| {
                let capacitance = parameters().bus_capacitance_f[port];
                let resistance = 0.5 * (1.0 + alpha * theta);
                let difference = initial_bus[port] - source[port];
                0.5 * capacitance
                    * difference
                    * difference
                    * -(-2.0 * dt / (resistance * capacitance)).exp_m1()
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
            ConverterElectricalDrive::Averaged(modulation(0.0, 0.0)),
            ConverterElectricalDrive::Switching {
                intervals: &intervals,
                controls: SwitchCycleControls::default(),
            },
        ] {
            let result = network
                .step_candidate(
                    &thermal,
                    dt,
                    &converter,
                    drive,
                    &source,
                    budget(),
                    heat_budget(),
                    controls(),
                    &[0.0; 6],
                    || false,
                )
                .unwrap();
            for id in &ALL_IDS[..5] {
                near(result.thermal().temperature_k(*id).unwrap(), 300.0);
            }
            near(
                result.thermal().temperature_k(ALL_IDS[5]).unwrap(),
                300.0 + theta,
            );
            near(result.converter().state()[0], 0.0);
            let electrical = &result.report().electrical;
            for port in 0..2 {
                let capacitance = parameters().bus_capacitance_f[port];
                let r = 0.5 * (1.0 + alpha * theta);
                let voltage = source[port]
                    + (initial_bus[port] - source[port]) * (-dt / (r * capacitance)).exp();
                let work = source[port] * capacitance * (voltage - initial_bus[port]);
                near(result.converter().state()[port + 1], capacitance * voltage);
                near(result.converter().terminal_resistance_ohm()[port], r);
                near(
                    electrical.terminal_resistor_heat_j[port],
                    expected_heat[port],
                );
                near(electrical.source_work_j[port], work);
                near(
                    electrical.converter.terminal_work_j[port],
                    work - expected_heat[port],
                );
            }
            assert_eq!(
                electrical.source_work_j[0] > 0.0,
                initial_bus[0] < source[0]
            );
            near(electrical.source_work_j[1], 0.0);
            for heat in electrical.converter.conduction.joules {
                near(heat, 0.0);
            }
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

#[test]
fn loaded_pwm_deposits_all_seven_native_losses_and_preserves_failed_histories() {
    let r0 = [R0[0], R0[1], R0[2], R0[3], R0[4], 0.25, 8.0];
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
    let modulation = modulation(0.3, 0.7);
    let intervals = centered_pwm(dt, modulation).unwrap();
    let drive = ConverterElectricalDrive::Switching {
        intervals: &intervals,
        controls: SwitchCycleControls::default(),
    };
    let native = AveragedBuckBoost::new(parameters(), modulation).unwrap();
    for (current, bus, source) in [
        (1.5, [24.0, 12.0], [26.0, 0.0]),
        (-1.5, [28.0, 18.0], [24.0, 0.0]),
    ] {
        let initial = native.state(current, bus).unwrap();
        let converter = network
            .bind(&thermal, &native, [r0[5], r0[6]], initial, || false)
            .unwrap();
        let calls = Cell::new(0usize);
        let result = network
            .step_candidate(
                &thermal,
                dt,
                &converter,
                drive,
                &source,
                budget(),
                heat_budget(),
                controls(),
                &[0.0; 7],
                || {
                    calls.set(calls.get() + 1);
                    false
                },
            )
            .unwrap();
        let electrical = &result.report().electrical;
        let heat: [f64; 7] = std::array::from_fn(|i| {
            if i < 5 {
                electrical.converter.conduction.joules[i]
            } else {
                electrical.terminal_resistor_heat_j[i - 5]
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
        let total_heat = heat.iter().sum::<f64>();
        near(result.report().deposited_joule_j, total_heat);
        near(result.thermal().energy_j() - thermal.energy_j(), total_heat);
        let storage = |x: &[f64; 3]| {
            0.5 * (x[0] * x[0] / parameters().inductance_h
                + x[1] * x[1] / parameters().bus_capacitance_f[0]
                + x[2] * x[2] / parameters().bus_capacitance_f[1])
        };
        near(
            electrical.converter.native.delta_h,
            storage(result.converter().state()) - storage(&initial),
        );
        near(
            electrical.converter.native.delta_h + total_heat,
            electrical.source_work_j.iter().sum(),
        );
        near(electrical.source_work_j[1], 0.0);
        assert_eq!(electrical.source_work_j[0] > 0.0, current > 0.0);
        near(result.report().total_energy_residual_j, 0.0);
        assert_eq!(result.report().circuit_trials % intervals.len(), 0);

        // Reuse the original physical solver at accepted material R(T), not a
        // frozen input current or a squared mean current for terminal heating.
        let mut replay = PreparedVoltageDrivenSwitchCycle::new_controlled(
            result.converter().parameters(),
            result.converter().terminal_resistance_ohm(),
            &intervals,
            budget(),
            heat_budget(),
            SwitchCycleControls::default(),
            || false,
        )
        .unwrap();
        let (mut state, mut voltage) = ([0.0; 3], [0.0; 2]);
        let exact = replay
            .step_controlled(&initial, &source, &mut state, &mut voltage, || false)
            .unwrap();
        for i in 0..3 {
            near(result.converter().state()[i], state[i]);
        }
        for port in 0..2 {
            near(
                electrical.terminal_resistor_heat_j[port],
                exact.terminal_resistor_heat_j[port],
            );
            near(electrical.source_work_j[port], exact.source_work_j[port]);
        }
        let last = calls.get();
        calls.set(0);
        assert!(matches!(
            network.step_candidate(
                &thermal,
                dt,
                &converter,
                drive,
                &source,
                budget(),
                heat_budget(),
                controls(),
                &[0.0; 7],
                || {
                    calls.set(calls.get() + 1);
                    calls.get() >= last
                },
            ),
            Err(CircuitEnthalpyError::Cancelled)
        ));
        let mut limited = controls();
        limited.maximum_circuit_trials = intervals.len() - 1;
        assert!(matches!(
            network.step_candidate(
                &thermal,
                dt,
                &converter,
                drive,
                &source,
                budget(),
                heat_budget(),
                limited,
                &[0.0; 7],
                || false,
            ),
            Err(CircuitEnthalpyError::CircuitBudget)
        ));
        assert!(
            network
                .step_candidate(
                    &thermal,
                    2.0 * dt,
                    &converter,
                    drive,
                    &source,
                    budget(),
                    heat_budget(),
                    controls(),
                    &[0.0; 7],
                    || false,
                )
                .is_err()
        );
        assert!(
            network
                .step_candidate(
                    &thermal,
                    dt,
                    &converter,
                    drive,
                    &source,
                    budget(),
                    heat_budget(),
                    controls(),
                    &[10000.0; 7],
                    || false,
                )
                .is_err()
        );
        assert!(
            network
                .step_candidate(
                    result.thermal(),
                    2.0 * dt,
                    &converter,
                    drive,
                    &source,
                    budget(),
                    heat_budget(),
                    controls(),
                    &[0.0; 7],
                    || false,
                )
                .is_err()
        );
        assert_eq!(converter.state(), &initial);
        assert_eq!(thermal.cells(), &original_cells);
        assert_eq!(thermal.time_s(), 0.0);
    }
}

#[test]
fn terminal_bindings_refuse_ideal_clamps_and_wrong_physical_source_units() {
    let cp = [1.0; 7];
    let (cards, curves) = sources([R0[0], R0[1], R0[2], R0[3], R0[4], 0.25, 8.0], cp, 0.0);
    let context = point(300.0);
    let bindings = bindings(&cards, &curves, &context);
    let resistance = key();
    let thermal = body(&curves, cp);
    let native = AveragedBuckBoost::new(parameters(), modulation(0.0, 0.0)).unwrap();
    let mut channels = paths(ALL_IDS, &resistance);
    channels[5].resistance = ConverterResistance::Fixed;
    let network = VoltageDrivenConverterThermalNetwork {
        cells: &bindings,
        links: &[],
        channels: &channels,
    };
    for terminal in [[0.0, 8.0], [f64::from_bits(1), 8.0]] {
        assert!(matches!(
            network.bind(&thermal, &native, terminal, [0.0; 3], || false),
            Err(CircuitEnthalpyError::Source(
                ConverterThermalError::Invalid(_)
            ))
        ));
    }
    let rho = PropertyKey::new(
        ElectricalPropertyBasis::Resistivity.property_name(),
        ElectricalPropertyBasis::Resistivity.dimensions(),
    );
    channels[5].resistance = ConverterResistance::Resistance {
        key: &rho,
        selection: MaterialPropertySelection::SingleClaimOnly,
    };
    assert!(matches!(
        VoltageDrivenConverterThermalNetwork {
            cells: &bindings,
            links: &[],
            channels: &channels
        }
        .bind(&thermal, &native, [0.25, 8.0], [0.0; 3], || false),
        Err(CircuitEnthalpyError::Source(
            ConverterThermalError::Invalid(_)
        ))
    ));
    assert_eq!(thermal.time_s(), 0.0);
}
