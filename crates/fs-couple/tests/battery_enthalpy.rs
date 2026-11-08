//! Native finite battery work heats original material cells and feeds back R(T).
#[path = "circuit_enthalpy/common.rs"]
mod common;

use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpySolve,
    battery::*,
    material::{ConductionLink, MaterialCell},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::{
    phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase},
    state_point::{
        MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS,
    },
    thermal::source::ThermalMaterialPoint,
};
use fs_phs::{
    SolveBudget,
    circuit::battery::{BatteryOcvCurve, BatteryParameters, BatterySource},
};
use fs_qty::{Density, QuantitySpec, Resistance, Temperature};
use fs_thermal::enthalpy::{CellId, MaterialId};

const IDS: [CellId; 2] = [CellId(30), CellId(10)];
const R: [f64; 2] = [0.2, 0.3];
const LINK: [ConductionLink; 1] = [ConductionLink {
    link: 7,
    distance_a_m: 0.5,
    distance_b_m: 0.5,
}];
fn key() -> PropertyKey {
    PropertyKey::new("battery-equivalent-resistance", Resistance::DIMS)
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
}
fn parameters() -> BatteryParameters {
    BatteryParameters {
        series_resistance_ohm: R[0],
        polarization_resistance_ohm: R[1],
        polarization_capacitance_f: 2.0,
        bus_capacitance_f: 0.5,
    }
}
fn native(capacity_c: f64, knots: &[(f64, f64)]) -> BatterySource {
    BatterySource::new(
        BatteryOcvCurve::new(capacity_c, knots).unwrap(),
        parameters(),
    )
    .unwrap()
}
fn electrical() -> SolveBudget {
    SolveBudget {
        relative: 2e-12,
        absolute: 2e-13,
        max_iterations: 30,
    }
}
fn controls() -> CircuitEnthalpySolve {
    let mut controls = common::config();
    controls.difference_step_j = 1e-5;
    controls.newton.absolute_tolerance = 1e-11;
    controls.cell_energy_tolerance_j = 1e-9;
    controls.thermal_energy_tolerance_j = 1e-9;
    controls.total_energy_tolerance_j = 1e-8;
    controls
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-8, "{a:.15e} != {b:.15e}");
}

// Authored equivalent-circuit/caloric laws, not measured battery chemistry.
fn material(
    r: f64,
    cp: f64,
    alpha: f64,
    resistance_upper_k: f64,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let mut claims = ClaimSet::new();
    for (property, low, high, upper) in [
        (heat.clone(), cp, cp, 2000.0),
        (density.clone(), 1.0, 1.0, 2000.0),
        (
            PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
            0.0,
            0.0,
            2000.0,
        ),
        (
            key(),
            r * (1.0 - 100.0 * alpha),
            r * (1.0 + 1700.0 * alpha),
            resistance_upper_k,
        ),
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
                    upper,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic battery/enthalpy fixture".into(),
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
            process: "authored battery fixture".into(),
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
fn materials(
    cp: [f64; 2],
    alpha: f64,
    upper: f64,
) -> ([MaterialCard; 2], [ResolvedHeatCapacityCurve; 2]) {
    let pairs: [_; 2] = core::array::from_fn(|i| material(R[i], cp[i], alpha, upper));
    (
        pairs.each_ref().map(|(card, _)| card.clone()),
        pairs.map(|(_, curve)| curve),
    )
}
fn bindings<'a>(
    cards: &'a [MaterialCard; 2],
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    context: &'a QueryPoint,
    thermal_key: &'a PropertyKey,
) -> [MaterialCell<'a>; 2] {
    core::array::from_fn(|i| MaterialCell {
        cell: IDS[i],
        material_id: MaterialId(i as u64),
        curve: curves[i].curve(),
        point: ThermalMaterialPoint {
            material: &cards[i],
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: Some((thermal_key, MaterialPropertySelection::SingleClaimOnly)),
    })
}
fn sourced(key: &PropertyKey) -> [BatteryThermalChannel<'_>; 2] {
    core::array::from_fn(|i| BatteryThermalChannel {
        cell: IDS[i],
        resistance: BatteryResistance::Source {
            key,
            selection: MaterialPropertySelection::SingleClaimOnly,
        },
    })
}

#[test]
fn crossing_an_ocv_kink_deposits_both_native_losses_with_independent_calorimetry() {
    let cp = [2.0, 3.0];
    let (cards, curves) = materials(cp, 0.0, 2000.0);
    let context = point(300.0);
    let thermal_key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = bindings(&cards, &curves, &context, &thermal_key);
    let paths = core::array::from_fn(|i| BatteryThermalChannel {
        cell: IDS[i],
        resistance: BatteryResistance::Fixed,
    });
    let network = BatteryThermalNetwork {
        cells: &cells,
        links: &LINK,
        channels: &paths,
    };
    let body = common::initial(&curves, 0.0, None);
    let native = native(1.0, &[(0.0, 3.0), (0.5, 3.0), (1.0, 4.2)]);
    let before = native.state(0.6, 0.2, 0.4).unwrap();
    let battery = network.bind(&body, &native, before, || false).unwrap();
    let dt = 0.03;
    let result = network
        .step_candidate(
            &body,
            dt,
            &battery,
            0.0,
            electrical(),
            16,
            controls(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    let after = *result.battery().state();
    assert!(before[0] > 0.5 && after[0] < 0.5);
    let chemical_energy = |q: f64| 3.0 * q + 1.2 * (q - 0.5).max(0.0).powi(2);
    let chemical_delta = chemical_energy(after[0]) - chemical_energy(before[0]);
    let current = -(after[0] - before[0]) / dt;
    let vp = f64::midpoint(before[1], after[1]) / parameters().polarization_capacitance_f;
    let heat = [dt * R[0] * current * current, dt * vp * vp / R[1]];
    let report = result.report();
    near(
        report.electrical.battery.chemical_energy_change_j,
        chemical_delta,
    );
    near(report.electrical.battery.series_heat_j, heat[0]);
    near(report.electrical.battery.polarization_heat_j, heat[1]);
    let midpoint_ocv = native
        .curve()
        .voltage_v(f64::midpoint(before[0], after[0]))
        .unwrap();
    assert!((report.electrical.battery.discrete_ocv_v - midpoint_ocv).abs() > 1e-3);
    for i in 0..2 {
        near(
            result.thermal().temperature_k(IDS[i]).unwrap(),
            300.0 + heat[i] / cp[i],
        );
    }
    near(
        result.thermal().energy_j() - body.energy_j(),
        heat[0] + heat[1],
    );
    near(report.deposited_joule_j, heat[0] + heat[1]);
    near(report.non_electrical_heat_j, 0.0);
    near(report.total_energy_residual_j, 0.0);
    assert_eq!(battery.state(), &before);
    assert_eq!(body.time_s(), 0.0);
}

#[test]
fn endpoint_material_resistance_changes_the_actual_charge_and_heat_solution() {
    let cp = [0.02, 0.03];
    let alpha = 0.004;
    let (cards, curves) = materials(cp, alpha, 2000.0);
    let context = point(300.0);
    let thermal_key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = bindings(&cards, &curves, &context, &thermal_key);
    let resistance = key();
    let paths = sourced(&resistance);
    let network = BatteryThermalNetwork {
        cells: &cells,
        links: &LINK,
        channels: &paths,
    };
    let body = common::initial(&curves, 0.0, None);
    let native = native(10.0, &[(0.0, 4.0), (1.0, 4.0)]);
    let before = native.state(0.8, 0.3, 3.0).unwrap();
    let battery = network.bind(&body, &native, before, || false).unwrap();
    let (dt, input) = (0.2, -1.0);
    let result = network
        .step_candidate(
            &body,
            dt,
            &battery,
            input,
            electrical(),
            16,
            controls(),
            &[0.0; 2],
            || false,
        )
        .unwrap();

    // Independent elimination of the two RC midpoint equations at fixed R.
    // A small contraction then closes the two authored R(T_end) laws; this
    // oracle uses no battery, thermal residual or material-resolution solver.
    let trajectory = |r: [f64; 2]| {
        let p = parameters();
        let a = 2.0 * p.polarization_capacitance_f + dt / r[1];
        let current = (4.0
            - 2.0 * p.polarization_capacitance_f * 0.3 / a
            - 3.0
            - dt * input / (2.0 * p.bus_capacitance_f))
            / (r[0] + dt / a + dt / (2.0 * p.bus_capacitance_f));
        let vp = (2.0 * p.polarization_capacitance_f * 0.3 + dt * current) / a;
        let vb = 3.0 + dt * (current + input) / (2.0 * p.bus_capacitance_f);
        (
            current,
            vp,
            vb,
            [dt * r[0] * current * current, dt * vp * vp / r[1]],
        )
    };
    let mut expected_r = R;
    for _ in 0..80 {
        let (_, _, _, heat) = trajectory(expected_r);
        expected_r = core::array::from_fn(|i| R[i] * (1.0 + alpha * heat[i] / cp[i]));
    }
    let (current, vp, vb, heat) = trajectory(expected_r);
    let expected_state = [
        before[0] - dt * current,
        2.0 * parameters().polarization_capacitance_f * vp - before[1],
        2.0 * parameters().bus_capacitance_f * vb - before[2],
    ];
    for i in 0..3 {
        near(result.battery().state()[i], expected_state[i]);
    }
    for i in 0..2 {
        near(
            result.thermal().temperature_k(IDS[i]).unwrap(),
            300.0 + heat[i] / cp[i],
        );
        near(
            result.battery().channels()[i].resistance_ohm(),
            expected_r[i],
        );
        assert_eq!(
            battery.channels()[i].material().unwrap().card_identity(),
            result.battery().channels()[i]
                .material()
                .unwrap()
                .card_identity()
        );
    }
    assert!((heat[0] - trajectory(R).3[0]).abs() > 1e-5);
    near(result.report().deposited_joule_j, heat.iter().sum());
    near(result.report().total_energy_residual_j, 0.0);
    assert_eq!(
        result.battery().parameters().polarization_capacitance_f,
        battery.parameters().polarization_capacitance_f
    );
    assert_eq!(
        result.battery().parameters().bus_capacitance_f,
        battery.parameters().bus_capacitance_f
    );

    // The accepted source snapshots can continue on the same coupled timeline.
    let (next_body, next_battery, _) = result.into_parts();
    let continued = network
        .step_candidate(
            &next_body,
            2.0 * dt,
            &next_battery,
            input,
            electrical(),
            16,
            controls(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    assert!(continued.battery().state()[0] < next_battery.state()[0]);
}

#[test]
fn discharge_and_regeneration_route_both_losses_once_to_an_explicit_common_receiver() {
    let cp = [2.0, 3.0];
    let (cards, curves) = materials(cp, 0.0, 2000.0);
    let context = point(300.0);
    let thermal_key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = bindings(&cards, &curves, &context, &thermal_key);
    let paths = core::array::from_fn(|_| BatteryThermalChannel {
        cell: IDS[0],
        resistance: BatteryResistance::Fixed,
    });
    let network = BatteryThermalNetwork {
        cells: &cells,
        links: &LINK,
        channels: &paths,
    };
    let body = common::initial(&curves, 0.0, None);
    let native = native(10.0, &[(0.0, 4.0), (1.0, 4.0)]);
    let dt = 0.2;
    for input in [-1.0, 1.0] {
        // Constant OCV gives exact stationary polarization and bus voltage,
        // while chemical charge still decreases/increases at one coulomb/s.
        let before = native
            .state(0.5, -input * R[1], 4.0 + input * (R[0] + R[1]))
            .unwrap();
        let battery = network.bind(&body, &native, before, || false).unwrap();
        let result = network
            .step_candidate(
                &body,
                dt,
                &battery,
                input,
                electrical(),
                16,
                controls(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        near(result.battery().state()[0], before[0] + dt * input);
        near(result.battery().state()[1], before[1]);
        near(result.battery().state()[2], before[2]);
        near(result.report().electrical.battery.series_heat_j, dt * R[0]);
        near(
            result.report().electrical.battery.polarization_heat_j,
            dt * R[1],
        );
        near(
            result.report().electrical.battery.terminal_work_j,
            dt * input * (4.0 + input * (R[0] + R[1])),
        );
        near(
            result.thermal().temperature_k(IDS[0]).unwrap(),
            300.0 + dt * (R[0] + R[1]) / cp[0],
        );
        near(result.thermal().temperature_k(IDS[1]).unwrap(), 300.0);
        near(result.report().deposited_joule_j, dt * (R[0] + R[1]));
        near(result.report().total_energy_residual_j, 0.0);
        assert!(
            result
                .battery()
                .channels()
                .iter()
                .all(|channel| channel.material().is_none())
        );
    }
}

#[test]
fn unsupported_sources_exhaustion_events_and_cancellation_publish_neither_history() {
    let (cards, curves) = materials([0.02, 0.03], 0.004, 301.0);
    let context = point(300.0);
    let thermal_key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = bindings(&cards, &curves, &context, &thermal_key);
    let resistance = key();
    let paths = sourced(&resistance);
    let network = BatteryThermalNetwork {
        cells: &cells,
        links: &LINK,
        channels: &paths,
    };
    let body = common::initial(&curves, 0.0, None);
    let native = native(10.0, &[(0.0, 4.0), (1.0, 4.0)]);
    let before = native.state(0.8, 0.3, 3.0).unwrap();
    let battery = network.bind(&body, &native, before, || false).unwrap();
    let initial_heat = body.energy_j();
    // The physical heating solution exceeds the selected 301 K resistance
    // support, although the original caloric owner remains valid to 2000 K.
    assert!(
        network
            .step_candidate(
                &body,
                0.2,
                &battery,
                -1.0,
                electrical(),
                16,
                controls(),
                &[0.0; 2],
                || false
            )
            .is_err()
    );
    assert!(matches!(
        network.step_candidate(
            &body,
            0.01,
            &battery,
            -1.0,
            electrical(),
            16,
            controls(),
            &[0.0; 2],
            || true
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    let mut budget = controls();
    budget.maximum_circuit_trials = 1;
    assert!(matches!(
        network.step_candidate(
            &body,
            0.01,
            &battery,
            -1.0,
            electrical(),
            16,
            budget,
            &[0.0; 2],
            || false
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    let accepted = network
        .step_candidate(
            &body,
            0.01,
            &battery,
            -1.0,
            electrical(),
            16,
            controls(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    assert!(
        network
            .step_candidate(
                accepted.thermal(),
                0.02,
                &battery,
                -1.0,
                electrical(),
                16,
                controls(),
                &[0.0; 2],
                || false
            )
            .is_err(),
        "stale battery clock"
    );
    let event_body = common::initial(&curves, 0.0, Some(0.005));
    let event_battery = network
        .bind(&event_body, &native, before, || false)
        .unwrap();
    assert!(
        network
            .step_candidate(
                &event_body,
                0.01,
                &event_battery,
                -1.0,
                electrical(),
                16,
                controls(),
                &[0.0; 2],
                || false
            )
            .is_err()
    );
    let empty = network
        .bind(
            &body,
            &native,
            native.state(0.0001, 0.3, 3.5).unwrap(),
            || false,
        )
        .unwrap();
    assert!(
        network
            .step_candidate(
                &body,
                0.01,
                &empty,
                -1.0,
                electrical(),
                16,
                controls(),
                &[0.0; 2],
                || false
            )
            .is_err(),
        "finite SOC must not be clipped"
    );
    assert_eq!(body.time_s(), 0.0);
    assert_eq!(body.energy_j(), initial_heat);
    assert_eq!(battery.state(), &before);
}
