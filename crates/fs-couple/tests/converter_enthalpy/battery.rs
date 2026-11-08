//! Finite battery, regulated load and all eight material heat receivers.
#[path = "battery/motor.rs"]
mod motor;

use super::*;
use fs_couple::circuit_enthalpy::converter::battery::*;
use fs_phs::{
    SolveBudget,
    circuit::{
        battery::{BatteryOcvCurve, BatteryParameters, BatterySource},
        converter::{
            battery_fed::BatteryFedBuckBoost,
            power_load::{PowerLoadControls, PowerLoadMode, RegulatedPowerLoad},
        },
    },
};

const ALL_IDS: [CellId; 8] = [
    IDS[0],
    IDS[1],
    IDS[2],
    IDS[3],
    IDS[4],
    CellId(106),
    CellId(117),
    CellId(128),
];
const RESISTANCES: [f64; 8] = [R0[0], R0[1], R0[2], R0[3], R0[4], 0.2, 0.3, 0.1];

struct Fixture {
    cards: [MaterialCard; 8],
    curves: [ResolvedHeatCapacityCurve; 8],
    cp: [f64; 8],
}
impl Fixture {
    fn new(cp: [f64; 8], alpha: f64) -> Self {
        let pairs: [_; 8] = core::array::from_fn(|i| material(RESISTANCES[i], cp[i], alpha));
        Self {
            cards: pairs.each_ref().map(|(card, _)| card.clone()),
            curves: pairs.map(|(_, curve)| curve),
            cp,
        }
    }
    fn body(&self) -> DepositionSimulation<'_> {
        let events = (0..8)
            .map(|i| TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: ALL_IDS[i],
                    segment_id: i as u64,
                    body_id: i as u64,
                    material_id: MaterialId(i as u64),
                    mass_kg: 1.0,
                    specific_h_j_kg: 100.0 * self.cp[i],
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
        let laws = self
            .curves
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
    fn cells<'a>(&'a self, context: &'a QueryPoint) -> [MaterialCell<'a>; 8] {
        core::array::from_fn(|i| MaterialCell {
            cell: ALL_IDS[i],
            material_id: MaterialId(i as u64),
            curve: self.curves[i].curve(),
            point: ThermalMaterialPoint {
                material: &self.cards[i],
                context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        })
    }
}
fn paths(key: &PropertyKey) -> [ConverterThermalChannel<'_>; 8] {
    core::array::from_fn(|i| ConverterThermalChannel {
        cell: ALL_IDS[i],
        resistance: ConverterResistance::Resistance {
            key,
            selection: MaterialPropertySelection::SingleClaimOnly,
        },
    })
}
fn electrical_controls() -> PowerLoadControls {
    PowerLoadControls {
        solve: SolveBudget {
            relative: 1e-13,
            absolute: 1e-12,
            max_iterations: 40,
        },
        backtracking_halvings: 8,
        absolute_energy_tolerance_j: 1e-10,
        relative_energy_tolerance: 1e-12,
    }
}
fn steady() -> (BatteryFedBuckBoost, [f64; 6], RegulatedPowerLoad, [f64; 8]) {
    let battery_parameters = BatteryParameters {
        series_resistance_ohm: RESISTANCES[5],
        polarization_resistance_ohm: RESISTANCES[6],
        polarization_capacitance_f: 2.0,
        bus_capacitance_f: 0.5,
    };
    let battery = BatterySource::new(
        BatteryOcvCurve::new(100.0, &[(0.0, 12.0), (1.0, 12.0)]).unwrap(),
        battery_parameters,
    )
    .unwrap();
    let native = AveragedBuckBoost::new(parameters(), modulation(0.7, 0.6))
        .unwrap()
        .with_battery_source(battery, RESISTANCES[7])
        .unwrap();
    let (a, b, current) = (0.7, 0.6, 1.0);
    let battery_current = a * current;
    let polarization = battery_current * RESISTANCES[6];
    let battery_bus = 12.0 - polarization - battery_current * RESISTANCES[5];
    let input_bus = battery_bus - battery_current * RESISTANCES[7];
    let effective = [
        R0[0],
        a * R0[1],
        (1.0 - a) * R0[2],
        b * R0[3],
        (1.0 - b) * R0[4],
    ];
    let output_bus = (a * input_bus - effective.iter().sum::<f64>() * current) / b;
    let state = native
        .state(
            0.8,
            polarization,
            battery_bus,
            current,
            [input_bus, output_bus],
        )
        .unwrap();
    let power = output_bus * b * current;
    let load = RegulatedPowerLoad::try_new(power, 10.0, 1.0, 50.0).unwrap();
    let heat = core::array::from_fn(|i| {
        if i < 5 {
            effective[i] * current * current
        } else {
            RESISTANCES[i] * battery_current * battery_current
        }
    });
    (native, state, load, heat)
}

#[test]
fn eight_channel_calorimetry_keeps_regulated_load_work_outside_the_thermal_body() {
    let fixture = Fixture::new([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0], 0.0);
    let context = point(300.0);
    let cells = fixture.cells(&context);
    let resistance = key();
    let channels = paths(&resistance);
    let network = BatteryFedConverterThermalNetwork {
        cells: &cells,
        links: &[],
        channels: &channels,
    };
    let body = fixture.body();
    let (native, before, load, heat_w) = steady();
    let converter = network.bind(&body, &native, before, || false).unwrap();
    let dt = 0.02;
    let result = network
        .step_candidate(
            &body,
            dt,
            &converter,
            native.modulation(),
            0.0,
            Some(load),
            electrical_controls(),
            controls(),
            &[0.0; 8],
            || false,
        )
        .unwrap();
    let report = result.report();
    assert_eq!(report.electrical.load.mode, PowerLoadMode::ConstantPower);
    near(load.requested_power_w(), 8.031);
    near(heat_w.iter().sum(), 0.369);
    for i in 0..8 {
        let expected = dt * heat_w[i];
        let caloric = fixture.cp[i] * (result.thermal().temperature_k(ALL_IDS[i]).unwrap() - 300.0);
        assert!(
            (caloric - expected).abs() < 2e-9,
            "heat channel {i}: {caloric} != {expected}"
        );
        near(report.electrical.heat_joules()[i], expected);
    }
    let expected_heat = dt * heat_w.iter().sum::<f64>();
    near(report.deposited_joule_j, expected_heat);
    near(result.thermal().energy_j() - body.energy_j(), expected_heat);
    near(report.electrical.load_work_j, dt * load.requested_power_w());
    near(report.electrical.external_work_j, 0.0);
    near(report.non_electrical_heat_j, 0.0);
    near(
        report.electrical.battery.chemical_energy_change_j,
        -dt * 12.0 * 0.7,
    );
    near(report.total_energy_residual_j, 0.0);
    assert!(report.electrical.load_work_j > 20.0 * report.deposited_joule_j);
}

#[test]
fn all_eight_endpoint_resistances_feed_back_into_the_same_six_state_solution() {
    let (_, _, _, reference_heat) = steady();
    let cp = reference_heat.map(|power| 0.02 * power);
    let mut states = Vec::new();
    for alpha in [0.0, 0.0039] {
        let fixture = Fixture::new(cp, alpha);
        let context = point(300.0);
        let cells = fixture.cells(&context);
        let resistance = key();
        let channels = paths(&resistance);
        let network = BatteryFedConverterThermalNetwork {
            cells: &cells,
            links: &[],
            channels: &channels,
        };
        let body = fixture.body();
        let (native, before, load, _) = steady();
        let converter = network.bind(&body, &native, before, || false).unwrap();
        let dt = 0.1;
        let result = network
            .step_candidate(
                &body,
                dt,
                &converter,
                native.modulation(),
                0.0,
                Some(load),
                electrical_controls(),
                controls(),
                &[0.0; 8],
                || false,
            )
            .unwrap();
        let next = *result.converter().state();
        let actual_r =
            core::array::from_fn::<_, 8, _>(|i| result.converter().channels()[i].resistance_ohm());
        let i_converter = f64::midpoint(before[3], next[3]) / parameters().inductance_h;
        let i_battery = -(next[0] - before[0]) / dt;
        let vp = f64::midpoint(before[1], next[1])
            / native.battery_parameters().polarization_capacitance_f;
        let vb = f64::midpoint(before[2], next[2]) / native.battery_parameters().bus_capacitance_f;
        let vi = f64::midpoint(before[4], next[4]) / parameters().bus_capacitance_f[0];
        let i_lead = (vb - vi) / actual_r[7];
        let factors = [1.0, 0.7, 0.3, 0.6, 0.4];
        for i in 0..8 {
            let temperature = result.thermal().temperature_k(ALL_IDS[i]).unwrap();
            near(
                actual_r[i],
                RESISTANCES[i] * (1.0 + alpha * (temperature - 300.0)),
            );
            let power = match i {
                0..=4 => factors[i] * actual_r[i] * i_converter * i_converter,
                5 => actual_r[i] * i_battery * i_battery,
                6 => vp * vp / actual_r[i],
                _ => actual_r[i] * i_lead * i_lead,
            };
            assert!(
                (cp[i] * (temperature - 300.0) - dt * power).abs() < 2e-9,
                "alpha={alpha}, channel={i}: caloric/source balance"
            );
            assert_eq!(
                converter.channels()[i].material().unwrap().card_identity(),
                result.converter().channels()[i]
                    .material()
                    .unwrap()
                    .card_identity()
            );
        }
        near(result.report().total_energy_residual_j, 0.0);
        states.push(next);
    }
    assert!(
        (states[1][3] - states[0][3]).abs() > 1e-6,
        "thermal resistance must alter native current/flux"
    );
}

#[test]
fn unsupported_load_and_cancelled_trials_leave_both_histories_unchanged() {
    let fixture = Fixture::new([1.0; 8], 0.0);
    let context = point(300.0);
    let cells = fixture.cells(&context);
    let resistance = key();
    let channels = paths(&resistance);
    let network = BatteryFedConverterThermalNetwork {
        cells: &cells,
        links: &[],
        channels: &channels,
    };
    let body = fixture.body();
    let (native, before, load, _) = steady();
    let converter = network.bind(&body, &native, before, || false).unwrap();
    let energy = body.energy_j();
    let unsupported =
        RegulatedPowerLoad::try_new(load.requested_power_w(), 10.0, 1.0, 5.0).unwrap();
    assert!(
        network
            .step_candidate(
                &body,
                0.02,
                &converter,
                native.modulation(),
                0.0,
                Some(unsupported),
                electrical_controls(),
                controls(),
                &[0.0; 8],
                || false
            )
            .is_err()
    );
    assert!(matches!(
        network.step_candidate(
            &body,
            0.02,
            &converter,
            native.modulation(),
            0.0,
            Some(load),
            electrical_controls(),
            controls(),
            &[0.0; 8],
            || true
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert_eq!(converter.state(), &before);
    assert_eq!(converter.time_s(), 0.0);
    assert_eq!(body.time_s(), 0.0);
    assert_eq!(body.energy_j(), energy);
    let accepted = network
        .step_candidate(
            &body,
            0.02,
            &converter,
            native.modulation(),
            0.0,
            Some(load),
            electrical_controls(),
            controls(),
            &[0.0; 8],
            || false,
        )
        .unwrap();
    near(accepted.report().total_energy_residual_j, 0.0);
}
