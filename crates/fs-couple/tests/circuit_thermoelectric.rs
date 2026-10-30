//! Actual thermocouple circuit + common material network + finite thermal stores.
#[path = "circuit_enthalpy/common.rs"]
mod common;
use common::{boundary, config, near, point};
use fs_couple::circuit_enthalpy::{CircuitEnthalpyError, ResistorCell};
use fs_couple::circuit_enthalpy::material::{CellResistor, ConductionLink, MaterialCell, MaterialNetwork};
use fs_couple::circuit_enthalpy::thermoelectric::{self, ThermoelectricLeg};
use fs_couple::material_circuit::MaterialCircuit;
use fs_evidence::ValidityDomain;
use fs_material::conductor::ElectricalPropertyBasis;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_material::thermoelectric::{ABSOLUTE_SEEBECK_PROPERTY, SEEBECK_DIMS};
use fs_matdb::{InterpolationPolicy, MaterialCard, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{Density, QuantitySpec, Temperature};
use fs_thermal::enthalpy::*;

#[path = "circuit_thermoelectric/motion.rs"]
mod motion_tests;
#[path = "circuit_thermoelectric/graph.rs"]
mod graph_tests;
#[path = "circuit_thermoelectric/heat.rs"]
mod heat_tests;

fn source(s: f64, resistance: f64, k: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    source_value(PropertyValue::Scalar { value: s, dims: SEEBECK_DIMS },
        InterpolationPolicy::ConstantWithinValidity, resistance, k)
}
fn source_value(value: PropertyValue, interpolation: InterpolationPolicy,
    resistance: f64, k: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (base, _) = common::material_with_conductivity(resistance, 10.0, 0.0, k);
    let mut claims = base.claims().clone();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(ABSOLUTE_SEEBECK_PROPERTY, SEEBECK_DIMS),
        value,
        validity: ValidityDomain::unconstrained()
            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 200.0, 2000.0)
            .with("lot", 7.0, 7.0),
        uncertainty: UncertaintyModel::Unstated, interpolation,
        observations: vec![], provenance: Provenance { source: "authored thermocouple fixture, not measured alloy data".into(),
            license: "CC0-1.0".into(), artifact: None },
    }).unwrap();
    let card = MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
        density: &PropertyKey::new("density", Density::DIMS), temperature_axis: "T",
        lower: &point(200.0), upper: &point(2000.0), selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 1e-6).unwrap();
    (card, caloric)
}
fn fixture(s1: f64, s2: f64, r2: f64, k: f64) -> ([MaterialCard; 4], [ResolvedHeatCapacityCurve; 4]) {
    let (a, ca) = source(0.0, 1.0, 1.0);
    let (b, cb) = source(0.0, 1.0, 1.0);
    let (c, cc) = source(s1, 2.0, k);
    let (d, cd) = source(s2, r2, k);
    ([a, b, c, d], [ca, cb, cc, cd])
}
fn initial(curves: &[ResolvedHeatCapacityCurve; 4], temperatures: [f64; 4]) -> DepositionSimulation<'_> {
    initial_with_links(curves, temperatures, &[])
}
fn initial_with_links<'a>(curves: &'a [ResolvedHeatCapacityCurve; 4], temperatures: [f64; 4],
    extra: &[[u64; 2]]) -> DepositionSimulation<'a> {
    let mut events = (0..4).map(|i| TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
        id: CellId(i as u64), segment_id: i as u64, body_id: i as u64, material_id: MaterialId(i as u64),
        mass_kg: 1.0, specific_h_j_kg: 10.0 * (temperatures[i] - 200.0),
        geometry: CellGeometry { start_m: [i as f64, 0.0, 0.0], end_m: [i as f64 + 1.0, 0.0, 0.0],
            volume_m3: 1.0, surface_area_m2: 6.0 }, boundary: boundary(),
    }) }).collect::<Vec<_>>();
    for (id, [a, b]) in [[0, 2], [2, 1], [1, 3], [3, 0]].into_iter()
        .chain(extra.iter().copied()).enumerate() {
        events.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
            id: id as u64, a: CellId(a), b: CellId(b), area_m2: 1.0, conductance_w_k: 999.0,
        }) });
    }
    let laws = (0..4).map(|i| (MaterialId(i as u64), curves[i].curve() as &dyn ThermalLaw)).collect::<Vec<_>>();
    DepositionSimulation::new(events, &laws, DepositionBudget {
        max_cells: 4, max_links: 4 + extra.len(), max_events: 8 + extra.len(), max_substeps: 100,
    }).unwrap().advance_to(0.0, StepConfig { max_dt_s: 0.1, energy_tolerance_j: 1e-7 }, || Ok(())).unwrap().0
}
fn cells<'a>(cards: &'a [MaterialCard; 4], curves: &'a [ResolvedHeatCapacityCurve; 4],
    context: &'a QueryPoint, key: &'a PropertyKey) -> [MaterialCell<'a>; 4] {
    [0, 1, 2, 3].map(|i| MaterialCell {
        cell: CellId(i as u64), material_id: MaterialId(i as u64), curve: curves[i].curve(),
        point: ThermalMaterialPoint { material: &cards[i], context, temperature_axis: "T", liquid_fraction_axis: None },
        conductivity: Some((key, MaterialPropertySelection::SingleClaimOnly)),
    })
}
fn resistors() -> [CellResistor; 2] {
    [(1, 2), (3, 3)].map(|(branch, cell)| CellResistor {
        destination: ResistorCell { branch, cell: CellId(cell) }, selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 1.0, area_m2: 1.0,
    })
}
fn links() -> [ConductionLink; 4] {
    [0, 1, 2, 3].map(|link| ConductionLink { link, distance_a_m: 0.5, distance_b_m: 0.5 })
}
fn legs() -> [ThermoelectricLeg; 2] {
    [ThermoelectricLeg { source_branch: 0, resistor_branch: 1, terminals: [CellId(0), CellId(1)],
        thermal_links: [0, 1], selection: MaterialPropertySelection::SingleClaimOnly },
     ThermoelectricLeg { source_branch: 2, resistor_branch: 3, terminals: [CellId(1), CellId(0)],
        thermal_links: [2, 3], selection: MaterialPropertySelection::SingleClaimOnly }]
}
const BASES: [(usize, ElectricalPropertyBasis); 2] = [
    (1, ElectricalPropertyBasis::Resistivity), (3, ElectricalPropertyBasis::Resistivity),
];
fn circuit(network: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'_>, current_driven: bool, drive: f64) -> MaterialCircuit {
    MaterialCircuit::new(CircuitGraph { node_count: 5, branches: vec![
        (1, 2, Branch::VoltageSource { port: 0 }), (2, 3, Branch::Resistor { ohms: 999.0 }),
        (3, 4, Branch::VoltageSource { port: 1 }), (4, 0, Branch::Resistor { ohms: 999.0 }),
        (1, 0, if current_driven { Branch::CurrentSource { port: 2 } } else { Branch::VoltageSource { port: 2 } }),
    ], transformers: vec![] }, network.resolve(thermal, || false).unwrap().resistors,
        vec![], vec![], &[0.0, 0.0, drive]).unwrap()
}

#[test]
fn open_circuit_two_material_voltage_uses_actual_conducting_junction_temperatures() {
    let (cards, curves) = fixture(0.02, -0.01, 3.0, 1.0);
    let context = point(300.0); let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key); let resistors = resistors(); let links = links();
    let material = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
    let thermal = initial(&curves, [400.0, 300.0, 350.0, 350.0]);
    let electrical = circuit(&material, &thermal, true, 0.0);
    let next = thermoelectric::step_candidate(&material, &thermal, 0.1, &electrical, &[0.0; 3], &legs(),
        &BASES, config(), &[0.0; 4], || false).unwrap();
    // Independent BE solution: the two 350 K cores remain at their symmetric
    // temperature and both unit-conductance links cool/heat the junctions.
    let hot = (10.0 * 400.0 + 0.2 * 350.0) / 10.2;
    let cold = (10.0 * 300.0 + 0.2 * 350.0) / 10.2;
    near(next.thermal().temperature_k(CellId(0)).unwrap(), hot);
    near(next.thermal().temperature_k(CellId(1)).unwrap(), cold);
    near(next.report().coupled.electrical.native.y[2], (0.02 - -0.01) * (cold - hot));
    near(next.report().coupled.deposited_joule_j, 0.0);
    near(next.report().external_electrical_work_j, 0.0);
    near(next.thermal().energy_j(), thermal.energy_j());
    for leg in &next.report().coupled.source.legs {
        near(leg.transfer.current_a, 0.0);
        for power in leg.transfer.terminal_heat_w { near(power, 0.0); }
        assert_eq!(leg.material.material()[0].card_identity(), cards[leg.resistor_branch / 2 + 2].content_hash());
    }
}

#[test]
fn current_reversal_reverses_peltier_heat_while_joule_loss_stays_positive() {
    let (cards, curves) = fixture(0.02, -0.01, 3.0, 1.0);
    let context = point(300.0); let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key); let resistors = resistors(); let links = links();
    let material = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
    let thermal = initial(&curves, [300.0; 4]);
    for current in [-0.5, 0.5] {
        let electrical = circuit(&material, &thermal, true, current);
        let next = thermoelectric::step_candidate(&material, &thermal, 0.2, &electrical,
            &[0.0, 0.0, current], &legs(), &BASES, config(), &[0.0; 4], || false).unwrap();
        let r = next.report(); let [a, b] = [CellId(0), CellId(1)].map(|id| next.thermal().temperature_k(id).unwrap());
        let l = &r.coupled.source.legs;
        near(l[0].transfer.current_a, current); near(l[1].transfer.current_a, current);
        near(l[0].transfer.terminal_heat_w[0] + l[1].transfer.terminal_heat_w[1], -(0.02 - -0.01) * a * current);
        near(l[0].transfer.terminal_heat_w[1] + l[1].transfer.terminal_heat_w[0], (0.02 - -0.01) * b * current);
        near(r.coupled.deposited_joule_j, 0.2 * current * current * (2.0 + 3.0));
        near(r.external_electrical_work_j, 0.2 * current * (current * 5.0 + 0.03 * (b - a)));
        near(next.thermal().energy_j() - thermal.energy_j(), r.external_electrical_work_j);
        near(r.external_heat_j, 0.0); near(r.total_energy_residual_j, 0.0);
        assert!((b - a) * current > 0.0);
        for sample in l {
            let [ta, tb] = sample.material.temperature_k();
            let [qa, qb] = sample.transfer.terminal_heat_w;
            near(qa / ta + qb / tb, 0.0); // Reversible entropy transfer.
        }
    }
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn material_substitution_changes_real_generator_current_and_internal_work_stays_internal() {
    let run = |s, r2, k| {
        let (cards, curves) = fixture(s, -0.01, r2, k);
        let context = point(300.0); let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = cells(&cards, &curves, &context, &key); let resistors = resistors(); let links = links();
        let material = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
        let thermal = initial(&curves, [400.0, 300.0, 350.0, 350.0]);
        let electrical = circuit(&material, &thermal, false, 0.0);
        let next = thermoelectric::step_candidate(&material, &thermal, 0.1, &electrical, &[0.0; 3],
            &legs(), &BASES, config(), &[0.0; 4], || false).unwrap();
        let [a, b] = [CellId(0), CellId(1)].map(|id| next.thermal().temperature_k(id).unwrap());
        let current = next.report().coupled.source.legs[0].transfer.current_a;
        near(current, (s - -0.01) * (a - b) / (2.0 + r2));
        near(next.report().coupled.source.network.links[0].conductance_w_k, 1.0 / (0.5 + 0.5 / k));
        near(next.report().external_electrical_work_j, 0.0);
        near(next.report().external_heat_j, 0.0);
        near(next.thermal().energy_j(), thermal.energy_j());
        near(next.report().coupled.source.reversible_electrical_work_j, next.report().coupled.deposited_joule_j);
        (current, a - b)
    };
    let reference = run(0.02, 3.0, 1.0);
    assert!(run(0.04, 3.0, 1.0).0 > reference.0);
    assert!(run(0.02, 6.0, 1.0).0 < reference.0);
    assert!(run(0.02, 3.0, 3.0).1 < reference.1);
}

#[test]
fn temperature_dependent_source_drives_bulk_heat_and_the_same_native_energy_balance() {
    let run = |slope: f64, current: f64| {
        let (mut cards, mut curves) = fixture(0.02, -0.01, 3.0, 1.0);
        (cards[2], curves[2]) = source_value(PropertyValue::Curve {
            abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0, 0.02 - 150.0*slope), (2000.0, 0.02 + 1650.0*slope)],
            dims: SEEBECK_DIMS,
        }, InterpolationPolicy::LinearInside, 2.0, 1.0);
        let context = point(300.0);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = cells(&cards, &curves, &context, &key);
        let resistors = resistors(); let links = links();
        let material = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
        let thermal = initial(&curves, [400.0, 300.0, 350.0, 350.0]);
        let electrical = circuit(&material, &thermal, true, current);
        let next = thermoelectric::step_candidate(&material, &thermal, 0.1, &electrical,
            &[0.0, 0.0, current], &legs(), &BASES, config(), &[0.0; 4], || false).unwrap();
        let report = next.report();
        let sample = &report.coupled.source.legs[0];
        let [a, b] = sample.material.temperature_k();
        let middle = next.thermal().temperature_k(CellId(2)).unwrap();
        let q = sample.transfer;
        // Independent affine-S integrals, including the interior-node transfer.
        let left = 0.02 + slope * (0.5*(a + middle) - 350.0);
        let right = 0.02 + slope * (0.5*(middle + b) - 350.0);
        near(q.voltage_drop_v, (0.02 + slope*(0.5*(a+b)-350.0))*(b-a));
        near(q.current_a, current);
        near(q.terminal_heat_w[0], -current*left*a);
        near(q.terminal_heat_w[1], current*right*b);
        near(q.bulk_heat_w, 0.5*current*slope*(a-b)*middle);
        near(q.terminal_heat_w[0]/a + q.bulk_heat_w/middle + q.terminal_heat_w[1]/b, 0.0);
        // The actual native bulk enthalpy receives conduction, Joule heat once,
        // and the same signed reversible transfer. No diagnostic-only heat.
        near(10.0*(middle-350.0), 0.1*(a+b-2.0*middle + 2.0*current*current + q.bulk_heat_w));
        near(report.coupled.deposited_joule_j, 0.1*5.0*current*current);
        near(report.total_energy_residual_j, 0.0);
        near(next.thermal().energy_j()-thermal.energy_j(), report.external_electrical_work_j);
        near(report.external_heat_j, 0.0);
        assert_eq!(sample.bulk, CellId(2));
        assert_eq!(sample.material.bulk_material().card_identity(), cards[2].content_hash());
        assert_eq!(thermal.time_s(), 0.0);
        middle
    };
    for current in [-0.5, 0.5] {
        let constant = run(0.0, current);
        let varying = run(2e-4, current);
        assert!((varying-constant)*current > 0.005);
    }
}

#[test]
fn seebeck_curve_support_includes_the_actual_bulk_and_refuses_atomically() {
    let (mut cards, mut curves) = fixture(0.02, -0.01, 3.0, 1.0);
    (cards[2], curves[2]) = source_value(PropertyValue::Curve {
        abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, 0.01), (400.0, 0.03)], dims: SEEBECK_DIMS,
    }, InterpolationPolicy::LinearInside, 2.0, 1.0);
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let resistors = resistors(); let links = links();
    let material = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
    let thermal = initial(&curves, [400.0, 300.0, 450.0, 350.0]);
    let electrical = circuit(&material, &thermal, true, 0.5);
    let saved_cells = thermal.cells().clone(); let saved_circuit = electrical.clone();
    assert!(matches!(thermoelectric::step_candidate(&material, &thermal, 0.1, &electrical,
        &[0.0, 0.0, 0.5], &legs(), &BASES, config(), &[0.0; 4], || false),
        Err(CircuitEnthalpyError::Source(thermoelectric::ThermoelectricNetworkError::Seebeck(_)))));
    assert_eq!(thermal.cells(), &saved_cells);
    assert_eq!(electrical, saved_circuit);
}

#[test]
fn geometry_routing_and_late_cancellation_refuse_without_changing_any_history() {
    let (cards, curves) = fixture(0.02, -0.01, 3.0, 1.0);
    let context = point(300.0); let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key); let resistors = resistors(); let mut links = links();
    let thermal = initial(&curves, [400.0, 300.0, 350.0, 350.0]);
    let saved = thermal.cells().clone();
    let electrical = circuit(&MaterialNetwork { cells: &cells, resistors: &resistors, links: &links }, &thermal, false, 0.0);
    let original = electrical.clone();
    links[0].distance_b_m = 0.6;
    assert!(thermoelectric::step_candidate(&MaterialNetwork { cells: &cells, resistors: &resistors, links: &links },
        &thermal, 0.1, &electrical, &[0.0; 3], &legs(), &BASES, config(), &[0.0; 4], || false).is_err());
    links[0].distance_b_m = 0.5;
    let material = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
    let mut wrong = legs(); wrong[0].thermal_links.swap(0, 1);
    assert!(thermoelectric::step_candidate(&material, &thermal, 0.1, &electrical, &[0.0; 3], &wrong,
        &BASES, config(), &[0.0; 4], || false).is_err());
    assert!(thermoelectric::step_candidate(&material, &thermal, 0.1, &electrical, &[1.0, 0.0, 0.0], &legs(),
        &BASES, config(), &[0.0; 4], || false).is_err());
    let mut polls = 0;
    thermoelectric::step_candidate(&material, &thermal, 0.1, &electrical, &[0.0; 3], &legs(),
        &BASES, config(), &[0.0; 4], || { polls += 1; false }).unwrap();
    let mut actual = 0;
    assert!(matches!(thermoelectric::step_candidate(&material, &thermal, 0.1, &electrical, &[0.0; 3], &legs(),
        &BASES, config(), &[0.0; 4], || { actual += 1; actual == polls }), Err(CircuitEnthalpyError::Cancelled)));
    assert_eq!(thermal.cells(), &saved); assert_eq!(electrical, original);
}
