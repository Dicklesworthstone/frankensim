//! Native averaged/switched conversion heats the original five material bodies.
#[path = "circuit_enthalpy/common.rs"]
mod common;
use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpySolve, converter::*, material::MaterialCell,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::{
    conductor::ElectricalPropertyBasis,
    phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase},
    state_point::{MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS},
    thermal::source::ThermalMaterialPoint,
};
use fs_phs::circuit::converter::{
    AveragedBuckBoost, BridgeModulation, BuckBoostParameters,
    switching::{
        BridgeSwitchState, PreparedSwitchCycle, SwitchCycleControls, SwitchInterval, centered_pwm,
    },
};
use fs_phs::linear_exact::LinearSplitBudget;
use fs_qty::{Density, QuantitySpec, Resistance, Temperature};
use fs_thermal::enthalpy::*;
use std::cell::Cell;

const IDS: [CellId; 5] = [CellId(91), CellId(12), CellId(73), CellId(24), CellId(55)];
const R0: [f64; 5] = [0.03, 0.02, 0.04, 0.015, 0.025];
fn key() -> PropertyKey {
    PropertyKey::new("linear-on-state-resistance", Resistance::DIMS)
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
}
fn parameters() -> BuckBoostParameters {
    BuckBoostParameters {
        inductance_h: 0.02,
        bus_capacitance_f: [0.2, 0.3],
        winding_resistance_ohm: R0[0],
        input_on_resistance_ohm: [R0[1], R0[2]],
        output_on_resistance_ohm: [R0[3], R0[4]],
    }
}
fn modulation(a: f64, b: f64) -> BridgeModulation {
    BridgeModulation {
        input_high_fraction: a,
        output_high_fraction: b,
    }
}
fn budget() -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: 3,
        max_omega_dt: 0.04,
    }
}
fn controls() -> CircuitEnthalpySolve {
    let mut c = common::config();
    c.difference_step_j = 1e-5;
    c.newton.absolute_tolerance = 1e-11;
    c.cell_energy_tolerance_j = 1e-9;
    c.thermal_energy_tolerance_j = 1e-8;
    c.total_energy_tolerance_j = 1e-8;
    c
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-7, "{a:.15e} != {b:.15e}");
}

// Manufactured linear R(T) and caloric laws, not semiconductor device data.
// Unit length/area makes bulk rho numerically equal to the winding resistance.
fn material(r: f64, cp: f64, alpha: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let basis = ElectricalPropertyBasis::Resistivity;
    let mut claims = ClaimSet::new();
    for (property, low, high) in [
        (heat.clone(), cp, cp),
        (density.clone(), 1.0, 1.0),
        (
            PropertyKey::new(basis.property_name(), basis.dimensions()),
            r * (1.0 - 100.0 * alpha),
            r * (1.0 + 1700.0 * alpha),
        ),
        (key(), r * (1.0 - 100.0 * alpha), r * (1.0 + 1700.0 * alpha)),
    ] {
        let dims = property.quantity().dims();
        claims
            .insert_claim(PropertyClaim {
                key: property,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, low), (2000.0, high)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    200.0,
                    2000.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic converter thermal fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "authored integration fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &heat,
            density: &density,
            temperature_axis: "T",
            lower: &point(200.0),
            upper: &point(2000.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-8,
    )
    .unwrap();
    (card, curve)
}
fn materials(cp: [f64; 5], alpha: f64) -> ([MaterialCard; 5], [ResolvedHeatCapacityCurve; 5]) {
    let pairs: [_; 5] = std::array::from_fn(|i| material(R0[i], cp[i], alpha));
    (
        pairs.each_ref().map(|(card, _)| card.clone()),
        pairs.map(|(_, curve)| curve),
    )
}
fn thermal<'a>(
    curves: &'a [ResolvedHeatCapacityCurve; 5],
    cp: [f64; 5],
) -> DepositionSimulation<'a> {
    let events = (0..5)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: IDS[i],
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
        .map(|(i, c)| (MaterialId(i as u64), c.curve() as &dyn ThermalLaw))
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
fn cells<'a>(
    cards: &'a [MaterialCard; 5],
    curves: &'a [ResolvedHeatCapacityCurve; 5],
    context: &'a QueryPoint,
) -> Vec<MaterialCell<'a>> {
    (0..5)
        .map(|i| MaterialCell {
            cell: IDS[i],
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
fn channels(key: &PropertyKey) -> [ConverterThermalChannel<'_>; 5] {
    std::array::from_fn(|i| ConverterThermalChannel {
        cell: IDS[i],
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
fn source_resistance_changes_actual_freewheel_current_and_five_finite_temperatures() {
    let active_r = R0[0] + R0[2] + R0[4];
    let cp = [
        0.01 * R0[0] / active_r,
        0.01,
        0.01 * R0[2] / active_r,
        0.01,
        0.01 * R0[4] / active_r,
    ];
    let alpha = 0.0039;
    let (cards, curves) = materials(cp, alpha);
    let context = point(300.0);
    let bindings = cells(&cards, &curves, &context);
    let resistance = key();
    let paths = channels(&resistance);
    let network = ConverterThermalNetwork {
        cells: &bindings,
        links: &[],
        channels: &paths,
    };
    let mut body = thermal(&curves, cp);
    let native = AveragedBuckBoost::new(parameters(), modulation(0.0, 0.0)).unwrap();
    let mut converter = network
        .bind(&body, &native, native.state(3.0, [0.0; 2]).unwrap(), || {
            false
        })
        .unwrap();
    let (mut expected_i, mut theta) = (3.0_f64, 0.0_f64);
    let dt = 0.1;
    for index in 1..=2 {
        let before = body.energy_j();
        let result = network
            .step_candidate(
                &body,
                index as f64 * dt,
                &converter,
                ConverterElectricalDrive::Averaged(modulation(0.0, 0.0)),
                &[0.0; 2],
                budget(),
                controls(),
                &[0.0; 5],
                || false,
            )
            .unwrap();
        // Independent scalar energy elimination. Equal alpha and Cj/Rj make
        // the three active temperature rises equal; electrical energy is
        // E0 exp(-2 R0_total (1+alpha theta_end) dt/L).
        let energy = 0.5 * parameters().inductance_h * expected_i * expected_i;
        let (mut low, mut high) = (theta, theta + energy / 0.01);
        for _ in 0..80 {
            let end = 0.5 * (low + high);
            let heat = energy
                * (1.0
                    - (-2.0 * active_r * (1.0 + alpha * end) * dt / parameters().inductance_h)
                        .exp());
            if 0.01 * (end - theta) > heat {
                high = end;
            } else {
                low = end;
            }
        }
        let next_theta = 0.5 * (low + high);
        let expected_heat = 0.01 * (next_theta - theta);
        expected_i *=
            (-active_r * (1.0 + alpha * next_theta) * dt / parameters().inductance_h).exp();
        for j in 0..5 {
            let active = j == 0 || j == 2 || j == 4;
            let temperature = 300.0 + if active { next_theta } else { 0.0 };
            near(result.thermal().temperature_k(IDS[j]).unwrap(), temperature);
            near(
                result.report().electrical.conduction.joules[j],
                if active {
                    expected_heat * R0[j] / active_r
                } else {
                    0.0
                },
            );
            near(
                result.converter().channels()[j].resistance_ohm(),
                R0[j] * (1.0 + alpha * (temperature - 300.0)),
            );
            assert_eq!(
                converter.channels()[j].material().unwrap().card_identity(),
                result.converter().channels()[j]
                    .material()
                    .unwrap()
                    .card_identity()
            );
        }
        near(
            result.converter().state()[0],
            parameters().inductance_h * expected_i,
        );
        near(result.converter().state()[1], 0.0);
        near(result.converter().state()[2], 0.0);
        near(result.report().electrical.native.delta_h, -expected_heat);
        near(result.report().electrical.native.supplied, 0.0);
        near(result.report().deposited_joule_j, expected_heat);
        near(result.thermal().energy_j() - before, expected_heat);
        near(result.report().total_energy_residual_j, 0.0);
        if index == 1 {
            near(next_theta, 5.5906860425769);
            near(expected_i, 1.846432765475921);
            let frozen = energy * (1.0 - (-2.0 * active_r * dt / parameters().inductance_h).exp());
            assert!((result.report().deposited_joule_j - frozen).abs() > 5e-4);
        }
        theta = next_theta;
        let (next_body, next_converter, _) = result.into_parts();
        body = next_body;
        converter = next_converter;
    }
}

#[test]
fn all_channels_route_once_during_forward_and_reverse_steady_conversion() {
    let cp = [1.0, 2.0, 3.0, 4.0, 5.0];
    let (cards, curves) = materials(cp, 0.0);
    let context = point(300.0);
    let bindings = cells(&cards, &curves, &context);
    let resistance = key();
    let mut paths = channels(&resistance);
    paths[4].resistance = ConverterResistance::Fixed;
    let network = ConverterThermalNetwork {
        cells: &bindings,
        links: &[],
        channels: &paths,
    };
    let body = thermal(&curves, cp);
    let (a, b, dt) = (0.7, 0.6, 0.02);
    let effective = [
        R0[0],
        a * R0[1],
        (1.0 - a) * R0[2],
        b * R0[3],
        (1.0 - b) * R0[4],
    ];
    let total: f64 = effective.iter().sum();
    for current in [3.0, -3.0] {
        let volts = [24.0, (a * 24.0 - total * current) / b];
        let native = AveragedBuckBoost::new(parameters(), modulation(a, b)).unwrap();
        let state = native.state(current, volts).unwrap();
        let converter = network.bind(&body, &native, state, || false).unwrap();
        let result = network
            .step_candidate(
                &body,
                dt,
                &converter,
                ConverterElectricalDrive::Averaged(modulation(a, b)),
                &[a * current, -b * current],
                budget(),
                controls(),
                &[0.0; 5],
                || false,
            )
            .unwrap();
        for j in 0..5 {
            let heat = dt * current * current * effective[j];
            near(result.report().electrical.conduction.joules[j], heat);
            near(
                result.thermal().temperature_k(IDS[j]).unwrap(),
                300.0 + heat / cp[j],
            );
        }
        for j in 0..3 {
            near(result.converter().state()[j], state[j]);
        }
        near(
            result.report().electrical.terminal_work_j[0],
            dt * a * current * volts[0],
        );
        near(
            result.report().electrical.terminal_work_j[1],
            -dt * b * current * volts[1],
        );
        near(result.report().electrical.native.delta_h, 0.0);
        near(
            result.report().deposited_joule_j,
            dt * current * current * total,
        );
        near(result.report().total_energy_residual_j, 0.0);
        assert!(result.converter().channels()[4].material().is_none());
        assert_eq!(
            result.report().electrical.terminal_work_j[0] > 0.0,
            current > 0.0
        );
        assert_eq!(
            result.report().electrical.terminal_work_j[1] > 0.0,
            current < 0.0
        );
    }
}

#[test]
fn switched_ripple_and_endpoint_sources_share_heat_state_and_atomic_refusals() {
    let cp = [0.05; 5];
    let alpha = 0.0039;
    let (cards, curves) = materials(cp, alpha);
    let context = point(300.0);
    let bindings = cells(&cards, &curves, &context);
    let resistance = key();
    let paths = channels(&resistance);
    let network = ConverterThermalNetwork {
        cells: &bindings,
        links: &[],
        channels: &paths,
    };
    let body = thermal(&curves, cp);
    let dt = 0.02;
    let modulation = modulation(0.3, 0.7);
    let intervals = centered_pwm(dt, modulation).unwrap();
    let native = AveragedBuckBoost::new(parameters(), modulation).unwrap();
    let initial = native.state(1.5, [24.0, 12.0]).unwrap();
    let converter = network.bind(&body, &native, initial, || false).unwrap();
    let current = [0.2, -0.1];
    let drive = ConverterElectricalDrive::Switching {
        intervals: &intervals,
        controls: SwitchCycleControls::default(),
    };
    let calls = Cell::new(0usize);
    let result = network
        .step_candidate(
            &body,
            dt,
            &converter,
            drive,
            &current,
            budget(),
            controls(),
            &[0.0; 5],
            || {
                calls.set(calls.get() + 1);
                false
            },
        )
        .unwrap();
    for j in 0..5 {
        let t = result.thermal().temperature_k(IDS[j]).unwrap();
        near(
            t,
            300.0 + result.report().electrical.conduction.joules[j] / cp[j],
        );
        near(
            result.converter().channels()[j].resistance_ohm(),
            R0[j] * (1.0 + alpha * (t - 300.0)),
        );
        assert!(result.report().electrical.conduction.joules[j] > 0.0);
    }
    // Independently replay the original native recipe at the *accepted* R(T).
    // This is a consistency check for the coupled owner; native switched tests
    // separately compare its trajectory to physical-current/voltage ODEs.
    let mut switched = PreparedSwitchCycle::new_controlled(
        result.converter().parameters(),
        &intervals,
        budget(),
        SwitchCycleControls::default(),
        || false,
    )
    .unwrap();
    let (mut state, mut voltage) = ([0.0; 3], [0.0; 2]);
    let exact = switched
        .step_controlled(&initial, &current, &mut state, &mut voltage, || false)
        .unwrap();
    for j in 0..3 {
        near(result.converter().state()[j], state[j]);
    }
    for j in 0..5 {
        near(
            result.report().electrical.conduction.joules[j],
            exact.conduction.joules[j],
        );
    }
    near(result.report().total_energy_residual_j, 0.0);
    assert_eq!(result.report().circuit_trials % intervals.len(), 0);
    let mut averaged = AveragedBuckBoost::new(result.converter().parameters(), modulation)
        .unwrap()
        .prepare_controlled(dt, budget(), || false)
        .unwrap();
    let average = averaged
        .step_controlled(&initial, &current, &mut state, &mut voltage, || false)
        .unwrap();
    assert!((average.native.dissipated - result.report().deposited_joule_j).abs() > 1e-6);
    let last = calls.get();
    calls.set(0);
    assert!(matches!(
        network.step_candidate(
            &body,
            dt,
            &converter,
            drive,
            &current,
            budget(),
            controls(),
            &[0.0; 5],
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
        network.step_candidate(
            &body,
            dt,
            &converter,
            drive,
            &current,
            budget(),
            limited,
            &[0.0; 5],
            || false
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    assert!(
        network
            .step_candidate(
                &body,
                2.0 * dt,
                &converter,
                drive,
                &current,
                budget(),
                controls(),
                &[0.0; 5],
                || false
            )
            .is_err()
    );
    assert!(
        network
            .step_candidate(
                &body,
                dt,
                &converter,
                drive,
                &current,
                budget(),
                controls(),
                &[10000.0; 5],
                || false
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
                &current,
                budget(),
                controls(),
                &[0.0; 5],
                || false
            )
            .is_err()
    );
    assert_eq!(converter.state(), &initial);
    near(body.time_s(), 0.0);
    for id in IDS {
        near(body.temperature_k(id).unwrap(), 300.0);
    }
}

#[test]
fn binding_refuses_wrong_source_units_and_ambiguous_channel_receivers() {
    let cp = [1.0; 5];
    let (cards, curves) = materials(cp, 0.0039);
    let context = point(300.0);
    let bindings = cells(&cards, &curves, &context);
    let resistance = key();
    let mut paths = channels(&resistance);
    let body = thermal(&curves, cp);
    let native = AveragedBuckBoost::new(parameters(), modulation(0.0, 0.0)).unwrap();
    let state = native.state(3.0, [0.0; 2]).unwrap();
    paths[4].cell = IDS[0];
    assert!(
        ConverterThermalNetwork {
            cells: &bindings,
            links: &[],
            channels: &paths
        }
        .bind(&body, &native, state, || false)
        .is_err()
    );
    paths[4].cell = IDS[4];
    let rho = PropertyKey::new(
        ElectricalPropertyBasis::Resistivity.property_name(),
        ElectricalPropertyBasis::Resistivity.dimensions(),
    );
    paths[1].resistance = ConverterResistance::Resistance {
        key: &rho,
        selection: MaterialPropertySelection::SingleClaimOnly,
    };
    // This high-side channel is off, but its unsupported source still refuses.
    assert!(
        ConverterThermalNetwork {
            cells: &bindings,
            links: &[],
            channels: &paths
        }
        .bind(&body, &native, state, || false)
        .is_err()
    );
}

#[test]
fn short_switch_period_continues_at_a_nonzero_clock_without_inventing_ticks() {
    let cp = [1.0; 5];
    let (cards, curves) = materials(cp, 0.0);
    let context = point(300.0);
    let bindings = cells(&cards, &curves, &context);
    let resistance = key();
    let paths = channels(&resistance);
    let network = ConverterThermalNetwork {
        cells: &bindings,
        links: &[],
        channels: &paths,
    };
    let body = thermal(&curves, cp)
        .advance_to(
            1.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 1e-8,
            },
            || Ok(()),
        )
        .unwrap()
        .0;
    let native = AveragedBuckBoost::new(parameters(), modulation(0.0, 0.0)).unwrap();
    let converter = network
        .bind(&body, &native, native.state(3.0, [0.0; 2]).unwrap(), || {
            false
        })
        .unwrap();
    let period = 1e-4;
    let intervals = [SwitchInterval {
        state: BridgeSwitchState {
            input_high: false,
            output_high: false,
        },
        duration_s: period,
    }];
    let drive = ConverterElectricalDrive::Switching {
        intervals: &intervals,
        controls: SwitchCycleControls::default(),
    };
    let result = network
        .step_candidate(
            &body,
            1.0 + period,
            &converter,
            drive,
            &[0.0; 2],
            budget(),
            controls(),
            &[0.0; 5],
            || false,
        )
        .unwrap();
    let resistance = R0[0] + R0[2] + R0[4];
    let end_i = 3.0 * (-resistance * period / parameters().inductance_h).exp();
    near(
        result.converter().state()[0],
        parameters().inductance_h * end_i,
    );
    near(
        result.report().deposited_joule_j,
        0.5 * parameters().inductance_h * (9.0 - end_i * end_i),
    );
    near(result.report().total_energy_residual_j, 0.0);
    assert_eq!(result.converter().time_s(), 1.0 + period);
    // Absolute clock tolerance must not inflate a cycle smaller than one
    // representable forward clock increment into a complete thermal step.
    let tiny = [SwitchInterval {
        duration_s: 0.25 * (1.0_f64.next_up() - 1.0),
        ..intervals[0]
    }];
    assert!(
        network
            .step_candidate(
                &body,
                1.0_f64.next_up(),
                &converter,
                ConverterElectricalDrive::Switching {
                    intervals: &tiny,
                    controls: SwitchCycleControls::default()
                },
                &[0.0; 2],
                budget(),
                controls(),
                &[0.0; 5],
                || false
            )
            .is_err()
    );
    assert!(
        network
            .step_candidate(
                &body,
                1.0 + 2.0 * period,
                &converter,
                drive,
                &[0.0; 2],
                budget(),
                controls(),
                &[0.0; 5],
                || false
            )
            .is_err()
    );
}
