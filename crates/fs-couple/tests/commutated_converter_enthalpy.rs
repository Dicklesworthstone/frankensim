//! G1/G3/G4: independent diode calorimetry, source feedback and joint refusal.
#[path = "commutated_converter_enthalpy/battery_motor.rs"]
mod battery_motor;
#[path = "circuit_enthalpy/common.rs"]
mod common;

use fs_couple::circuit_enthalpy::{
    CircuitEnthalpySolve,
    converter::{ConverterResistance, ConverterThermalChannel, commutation::*},
    material::MaterialCell,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::{
    phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase},
    state_point::{MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS},
    thermal::source::ThermalMaterialPoint,
};
use fs_phs::{
    SolveBudget,
    circuit::converter::{
        commutation::{
            BodyDiode, CommutationInterval, CommutationParameters, DeviceGateState,
            DeviceOperatingBounds, SwitchDevice,
        },
        power_load::{PowerLoadControls, RegulatedPowerLoad},
        switching::{SwitchCycleControls, power_load::PowerSwitchCommand},
    },
};
use fs_qty::{Current, Density, QuantitySpec, Resistance, Temperature, Voltage};
use fs_thermal::enthalpy::*;
use std::cell::Cell;

const IDS: [CellId; 7] = [
    CellId(10),
    CellId(20),
    CellId(30),
    CellId(40),
    CellId(50),
    CellId(60),
    CellId(70),
];

struct Fixture {
    card: MaterialCard,
    caloric: ResolvedHeatCapacityCurve,
    keys: [PropertyKey; 7],
    cp: f64,
}

fn fixture(cp: f64, threshold_slope: f64, upper_k: f64, channel_limit: f64) -> Fixture {
    let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let keys = [
        PropertyKey::new("channel_on_resistance", Resistance::DIMS),
        PropertyKey::new("channel_off_resistance", Resistance::DIMS),
        PropertyKey::new("body_diode_forward_voltage", Voltage::DIMS),
        PropertyKey::new("body_diode_differential_resistance", Resistance::DIMS),
        PropertyKey::new("device_reverse_voltage_limit", Voltage::DIMS),
        PropertyKey::new("channel_current_limit", Current::DIMS),
        PropertyKey::new("diode_current_limit", Current::DIMS),
    ];
    let values = [0.1, 100.0, 0.5, 0.2, 100.0, channel_limit, 100.0];
    let mut claims = ClaimSet::new();
    let mut properties = vec![(heat.clone(), cp, 0.0), (density.clone(), 1.0, 0.0)];
    properties.extend(
        keys.iter()
            .cloned()
            .zip(values)
            .enumerate()
            .map(|(i, (key, value))| (key, value, if i == 2 { threshold_slope } else { 0.0 })),
    );
    for (key, value, slope) in properties {
        let dims = key.quantity().dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![
                        (299.0, value - slope),
                        (upper_k, value + slope * (upper_k - 300.0)),
                    ],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        299.0,
                        upper_k,
                    ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic conduction and caloric fixture, not qualified device data"
                        .into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic-device".into(),
            phase: "solid".into(),
            process: "authored fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &heat,
            density: &density,
            temperature_axis: "T",
            lower: &common::point(299.0),
            upper: &common::point(upper_k),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    Fixture {
        card,
        caloric,
        keys,
        cp,
    }
}

impl Fixture {
    fn thermal(&self) -> DepositionSimulation<'_> {
        self.thermal_at(300.0)
    }
    /// The seven bodies born at one common temperature (the caloric curve starts at 299 K), time 0.
    fn thermal_at(&self, kelvin: f64) -> DepositionSimulation<'_> {
        let events = IDS
            .into_iter()
            .enumerate()
            .map(|(i, id)| TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id,
                    segment_id: i as u64,
                    body_id: i as u64,
                    material_id: MaterialId(0),
                    mass_kg: 1.0,
                    specific_h_j_kg: self.cp * (kelvin - 299.0),
                    geometry: CellGeometry {
                        start_m: [i as f64, 0.0, 0.0],
                        end_m: [i as f64 + 1.0, 0.0, 0.0],
                        volume_m3: 1.0,
                        surface_area_m2: 6.0,
                    },
                    boundary: common::boundary(),
                }),
            })
            .collect();
        DepositionSimulation::new(
            events,
            &[(MaterialId(0), self.caloric.curve() as &dyn ThermalLaw)],
            DepositionBudget {
                max_cells: 7,
                max_links: 0,
                max_events: 7,
                max_substeps: 20,
            },
        )
        .unwrap()
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 1e-10,
            },
            || Ok(()),
        )
        .unwrap()
        .0
    }
    fn cells<'a>(&'a self, context: &'a QueryPoint) -> [MaterialCell<'a>; 7] {
        IDS.map(|cell| MaterialCell {
            cell,
            material_id: MaterialId(0),
            curve: self.caloric.curve(),
            point: ThermalMaterialPoint {
                material: &self.card,
                context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        })
    }
    fn devices(&self) -> [DeviceConductionBinding<'_>; 4] {
        core::array::from_fn(|i| DeviceConductionBinding {
            cell: IDS[i + 1],
            keys: DeviceConductionKeys {
                on_resistance: &self.keys[0],
                off_resistance: &self.keys[1],
                diode_forward_voltage: &self.keys[2],
                diode_differential_resistance: &self.keys[3],
                maximum_reverse_voltage: &self.keys[4],
                maximum_channel_current: &self.keys[5],
                maximum_diode_current: &self.keys[6],
            },
            selection: MaterialPropertySelection::SingleClaimOnly,
        })
    }
}

fn resistors() -> [ConverterThermalChannel<'static>; 3] {
    [IDS[0], IDS[5], IDS[6]].map(|cell| ConverterThermalChannel {
        cell,
        resistance: ConverterResistance::Fixed,
    })
}
fn parameters() -> CommutationParameters {
    CommutationParameters {
        inductance_h: 0.3,
        bus_capacitance_f: [0.2, 0.4],
        winding_resistance_ohm: 0.1,
        devices: [SwitchDevice {
            on_resistance_ohm: 0.1,
            off_resistance_ohm: 100.0,
            body_diode: BodyDiode {
                forward_voltage_v: 0.5,
                differential_resistance_ohm: 0.2,
            },
            bounds: DeviceOperatingBounds {
                maximum_reverse_voltage_v: 100.0,
                maximum_channel_current_a: 100.0,
                maximum_diode_current_a: 100.0,
            },
        }; 4],
    }
}
fn power_controls() -> PowerLoadControls {
    PowerLoadControls {
        solve: SolveBudget {
            relative: 1e-13,
            absolute: 1e-12,
            max_iterations: 40,
        },
        backtracking_halvings: 10,
        absolute_energy_tolerance_j: 1e-10,
        relative_energy_tolerance: 1e-12,
    }
}
fn cycle_controls() -> SwitchCycleControls {
    SwitchCycleControls {
        max_intervals: 64,
        absolute_energy_tolerance_j: 1e-9,
        relative_energy_tolerance: 1e-11,
    }
}
fn controls() -> CircuitEnthalpySolve {
    let mut c = common::config();
    c.difference_step_j = 1e-7;
    c.newton.absolute_tolerance = 1e-11;
    c.cell_energy_tolerance_j = 1e-9;
    c.thermal_energy_tolerance_j = 1e-9;
    c.total_energy_tolerance_j = 1e-9;
    c
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-7, "{a:e} != {b:e}");
}
fn interval(dt: f64) -> [CommutationInterval; 1] {
    [CommutationInterval {
        gates: DeviceGateState {
            on: [true, false, false, false],
        },
        duration_s: dt,
    }]
}

// Independent steady KCL/KVL: input high channel on; output both off with
// its high-side body diode conducting. Explicit off leakage is retained.
fn steady(power_w: f64) -> ([f64; 3], PowerSwitchCommand, [f64; 4], [f64; 2]) {
    let (current, output_v) = (1.0, 2.0);
    let (g_on, g_off, g_diode) = (10.0, 0.01, 5.0);
    let right_node =
        ((g_off + g_diode) * output_v + g_diode * 0.5 + current) / (2.0 * g_off + g_diode);
    let left_node = right_node + 0.1 * current;
    let input_v = ((g_on + g_off) * left_node + current) / g_on;
    let source_i = [
        current + g_off * left_node,
        g_off * right_node - current + power_w / output_v,
    ];
    let command = PowerSwitchCommand {
        source_voltage_v: [input_v + source_i[0], output_v + source_i[1]],
        loads: [
            None,
            Some(RegulatedPowerLoad::try_new(power_w, 10.0, 0.25, 50.0).unwrap()),
        ],
    };
    (
        parameters().state(current, [input_v, output_v]).unwrap(),
        command,
        [
            input_v - left_node,
            left_node,
            output_v - right_node,
            right_node,
        ],
        source_i,
    )
}

#[test]
fn independent_dead_time_calorimetry_keeps_channel_diode_source_and_load_work() {
    let f = fixture(1.0, 0.0, 350.0, 100.0);
    let thermal = f.thermal();
    let context = common::point(300.0);
    let cells = f.cells(&context);
    let devices = f.devices();
    let resistors = resistors();
    let network = CommutatedConverterThermalNetwork {
        cells: &cells,
        links: &[],
        resistors: &resistors,
        devices: &devices,
    };
    let (initial, command, vds, source_i) = steady(0.5);
    let converter = network
        .bind(&thermal, parameters(), [1.0; 2], initial, || false)
        .unwrap();
    let dt = 1.0 / 64.0;
    let intervals = interval(dt);
    let run = |receiver| {
        network
            .step_candidate(
                &thermal,
                dt,
                &converter,
                &intervals,
                &[command],
                4,
                [None, receiver],
                power_controls(),
                cycle_controls(),
                controls(),
                &[0.0; 7],
                || false,
            )
            .unwrap()
    };
    let exported = run(None);
    let heated = run(Some(IDS[3]));
    let step = &heated.report().electrical;
    near(step.native.delta_h, 0.0);
    near(step.winding_heat_j, dt * 0.1);
    let expected_channel = core::array::from_fn::<_, 4, _>(|i| {
        dt * vds[i] * vds[i] / if i == 0 { 0.1 } else { 100.0 }
    });
    let diode_current = (-vds[2] - 0.5) / 0.2;
    let diode_heat = dt * (-vds[2]) * diode_current;
    for i in 0..4 {
        near(step.channel_heat_j[i], expected_channel[i]);
        near(step.diode_heat_j[i], if i == 2 { diode_heat } else { 0.0 });
        let shared_load = if i == 2 { 0.5 * dt } else { 0.0 };
        near(
            heated.thermal().temperature_k(IDS[i + 1]).unwrap() - 300.0,
            expected_channel[i] + step.diode_heat_j[i] + shared_load,
        );
    }
    for port in 0..2 {
        near(
            step.terminal_resistor_heat_j[port],
            dt * source_i[port] * source_i[port],
        );
        near(
            step.source_work_j[port],
            dt * command.source_voltage_v[port] * source_i[port],
        );
        near(
            heated.thermal().temperature_k(IDS[5 + port]).unwrap() - 300.0,
            step.terminal_resistor_heat_j[port],
        );
    }
    near(step.load_work_j[1], 0.5 * dt);
    assert!(step.source_work_j[1] < 0.0);
    assert_eq!(step.integration_steps, 4);
    assert_eq!(heated.report().circuit_trials % 4, 0);
    near(
        heated.report().deposited_joule_j - exported.report().deposited_joule_j,
        0.5 * dt,
    );
    near(
        heated.thermal().energy_j() - thermal.energy_j(),
        step.source_work_j.iter().sum(),
    );
    near(heated.report().total_energy_residual_j, 0.0);
    for i in 0..3 {
        near(heated.converter().state()[i], initial[i]);
    }
}

#[test]
fn temperature_dependent_diode_threshold_refreshes_when_every_resistance_is_constant() {
    let f = fixture(0.001, 0.02, 350.0, 100.0);
    let thermal = f.thermal();
    let context = common::point(300.0);
    let cells = f.cells(&context);
    let devices = f.devices();
    let resistors = resistors();
    let network = CommutatedConverterThermalNetwork {
        cells: &cells,
        links: &[],
        resistors: &resistors,
        devices: &devices,
    };
    let (initial, command, _, _) = steady(0.0);
    let converter = network
        .bind(&thermal, parameters(), [1.0; 2], initial, || false)
        .unwrap();
    let dt = 1.0 / 1024.0;
    let candidate = network
        .step_candidate(
            &thermal,
            dt,
            &converter,
            &interval(dt),
            &[command],
            1,
            [None; 2],
            power_controls(),
            cycle_controls(),
            controls(),
            &[0.0; 7],
            || false,
        )
        .unwrap();
    let step = &candidate.report().electrical;
    let device = candidate.converter().devices()[2].parameters();
    let theta = candidate.thermal().temperature_k(IDS[3]).unwrap() - 300.0;
    assert!(theta > 0.1);
    near(device.on_resistance_ohm, 0.1);
    near(device.off_resistance_ohm, 100.0);
    near(device.body_diode.forward_voltage_v, 0.5 + 0.02 * theta);
    near(f.cp * theta, step.channel_heat_j[2] + step.diode_heat_j[2]);
    let point = step.last_devices[2];
    let current = ((-point.drain_source_voltage_v - device.body_diode.forward_voltage_v)
        / device.body_diode.differential_resistance_ohm)
        .max(0.0);
    near(point.diode_current_a, current);
    near(
        step.diode_heat_j[2],
        dt * -point.drain_source_voltage_v * current,
    );
    near(candidate.report().total_energy_residual_j, 0.0);
    near(
        converter.devices()[2]
            .parameters()
            .body_diode
            .forward_voltage_v,
        0.5,
    );
}

#[test]
fn source_limits_quantity_domains_budget_and_cancellation_preserve_joint_state() {
    for (upper, limit) in [(300.001, 100.0), (350.0, 0.01), (350.0, 100.0)] {
        let f = fixture(0.001, 0.0, upper, limit);
        let thermal = f.thermal();
        let energy = thermal.energy_j();
        let context = common::point(300.0);
        let cells = f.cells(&context);
        let devices = f.devices();
        let resistors = resistors();
        let network = CommutatedConverterThermalNetwork {
            cells: &cells,
            links: &[],
            resistors: &resistors,
            devices: &devices,
        };
        let (initial, command, _, _) = steady(0.0);
        let converter = network
            .bind(&thermal, parameters(), [1.0; 2], initial, || false)
            .unwrap();
        let dt = 1.0 / 1024.0;
        let intervals = interval(dt);
        let run = |target, solve, cancel: &mut dyn FnMut() -> bool| {
            network.step_candidate(
                &thermal,
                target,
                &converter,
                &intervals,
                &[command],
                4,
                [None; 2],
                power_controls(),
                cycle_controls(),
                solve,
                &[0.0; 7],
                cancel,
            )
        };
        if upper < 301.0 || limit < 1.0 {
            assert!(run(dt, controls(), &mut || false).is_err());
        } else {
            let calls = Cell::new(0usize);
            let expected = run(dt, controls(), &mut || {
                calls.set(calls.get() + 1);
                false
            })
            .unwrap();
            let stop = calls.get() / 2;
            calls.set(0);
            assert!(
                run(dt, controls(), &mut || {
                    calls.set(calls.get() + 1);
                    calls.get() >= stop
                })
                .is_err()
            );
            let retry = run(dt, controls(), &mut || false).unwrap();
            assert_eq!(retry.converter().state(), expected.converter().state());
            assert_eq!(retry.thermal().energy_j(), expected.thermal().energy_j());
            let mut too_small = controls();
            too_small.maximum_circuit_trials = 3;
            assert!(run(dt, too_small, &mut || false).is_err());
            assert!(run(2.0 * dt, controls(), &mut || false).is_err());
            let mut wrong = f.devices();
            wrong[0].keys.on_resistance = &f.keys[2];
            let incompatible = CommutatedConverterThermalNetwork {
                devices: &wrong,
                ..network
            };
            assert!(
                incompatible
                    .bind(&thermal, parameters(), [1.0; 2], initial, || false)
                    .is_err()
            );
        }
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(thermal.energy_j(), energy);
        assert_eq!(converter.state(), &initial);
    }
}

// The accepted-state check compares the sources a step resolves at the accepted heat state with
// the bound ones. The bound ones come from the binding's selection (SingleClaimOnly), the
// resolved ones from a refresh pinned to the claim already selected; equal sources must pass.
#[test]
fn a_converter_bound_under_single_claim_only_takes_its_first_step_and_another_claim_is_another_source() {
    let f = fixture(1.0, 0.0, 350.0, 100.0);
    let thermal = f.thermal();
    let context = common::point(300.0);
    let cells = f.cells(&context);
    let devices = f.devices();
    let resistors = resistors();
    let network = CommutatedConverterThermalNetwork {
        cells: &cells,
        links: &[],
        resistors: &resistors,
        devices: &devices,
    };
    let (initial, command, _, _) = steady(0.0);
    let converter = network
        .bind(&thermal, parameters(), [1.0; 2], initial, || false)
        .unwrap();
    let dt = 1.0 / 1024.0;
    match network.step_candidate(
        &thermal,
        dt,
        &converter,
        &interval(dt),
        &[command],
        1,
        [None; 2],
        power_controls(),
        cycle_controls(),
        controls(),
        &[0.0; 7],
        || false,
    ) {
        Ok(_) => {}
        Err(error) => panic!("a step from the bound state is refused: {error:?}"),
    }
    // The comparison is not vacuous: another claim value is another source.
    // (Sources refreshed by the step sit at the stepped state, so their query point differs from the
    // bound one by design; the accepted step above is the positive.)
    let other = fixture(1.0, 0.0, 350.0, 50.0);
    let other_thermal = other.thermal();
    let other_cells = other.cells(&context);
    let other_devices = other.devices();
    let other_network = CommutatedConverterThermalNetwork {
        cells: &other_cells,
        links: &[],
        resistors: &resistors,
        devices: &other_devices,
    };
    let bound_other = other_network
        .bind(&other_thermal, parameters(), [1.0; 2], initial, || false)
        .unwrap();
    assert_ne!(bound_other.devices(), converter.devices());
}

// Planted negative: the converter stays bound at 300 K and the accepted heat state is another one at
// the same time. Temperature is a coordinate of every device and resistor source's query point, so
// the sources resolved at the accepted state differ from the bound ones even with a constant diode
// threshold. The step must still refuse with the accepted-state message. (The threshold is constant
// so that the step from the bound state is the dead-time test's own, which fs-phs accepts; a
// temperature-dependent threshold meets the solve-budget refusal of
// temperature_dependent_diode_threshold_refreshes_when_every_resistance_is_constant.)
#[test]
fn a_converter_bound_at_one_state_refuses_a_step_from_another_accepted_state_at_the_same_time() {
    let f = fixture(1.0, 0.0, 350.0, 100.0);
    let thermal = f.thermal();
    let hot = f.thermal_at(320.0);
    assert!((hot.temperature_k(IDS[1]).unwrap() - 320.0).abs() < 1e-5);
    assert_eq!(hot.time_s().to_bits(), thermal.time_s().to_bits());
    let context = common::point(300.0);
    let cells = f.cells(&context);
    let devices = f.devices();
    let resistors = resistors();
    let network = CommutatedConverterThermalNetwork {
        cells: &cells,
        links: &[],
        resistors: &resistors,
        devices: &devices,
    };
    let (initial, command, _, _) = steady(0.5);
    let converter = network
        .bind(&thermal, parameters(), [1.0; 2], initial, || false)
        .unwrap();
    let dt = 1.0 / 64.0;
    let attempt = |at: &DepositionSimulation<'_>| {
        network
            .step_candidate(
                at,
                dt,
                &converter,
                &interval(dt),
                &[command],
                4,
                [None; 2],
                power_controls(),
                cycle_controls(),
                controls(),
                &[0.0; 7],
                || false,
            )
            .map(|_| ())
            .map_err(|error| format!("{error:?}"))
    };
    assert_eq!(attempt(&thermal), Ok(()));
    let Err(text) = attempt(&hot) else {
        panic!("a step from another accepted state was taken");
    };
    assert!(
        text.contains(
            "all commutated device and resistor sources must match the accepted heat state"
        ),
        "{text}"
    );
}
