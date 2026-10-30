//! Manufactured principal-source tests; no measured material qualification.
use super::*;
use fs_blake3::ContentHash;
use fs_material::phase::{HeatCapacitySourceQuery, SolidLiquidPhase};
use fs_material::state_point::SPECIFIC_HEAT_CAPACITY_DIMS;
use fs_matdb::{ClaimSet, PropertyValue, SelectionPolicy};
use fs_qty::{Density, Dims, QuantitySpec, Temperature};
use std::collections::BTreeMap;

const ROLES: [&str; 3] = ["k_longitudinal", "k_radial", "k_tangential"];
const FIELD: Dims = Dims([1, 1, -3, 0, -1, 0]);

fn principal_material(cp: f64, competing: bool, field: bool)
    -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (original, _) = material(1.0, cp, 0.0);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        if claim.key.name() != "thermal_conductivity" {
            claims.insert_claim(claim.clone()).unwrap();
            continue;
        }
        for (i, value) in [1.0, 4.0, 9.0].into_iter().enumerate() {
            let mut principal = claim.clone();
            principal.key = PropertyKey::new(ROLES[i], THERMAL_CONDUCTIVITY_DIMS);
            principal.value = if field && i == 0 {
                let dims = FIELD;
                principal.validity = principal.validity.with_quantity("E", QuantitySpec::dimensional(dims), 0.0, 20.0);
                PropertyValue::Curve { abscissa: "E".into(), abscissa_dims: dims,
                    knots: vec![(0.0, 1.0), (20.0, 11.0)], dims: THERMAL_CONDUCTIVITY_DIMS }
            } else {
                PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, 0.8 * value), (2000.0, 4.4 * value)], dims: THERMAL_CONDUCTIVITY_DIMS }
            };
            if competing && i == 0 {
                let mut alternate = principal.clone();
                principal.validity = principal.validity.with_quantity("T",
                    QuantitySpec::dimensional(Temperature::DIMS), 200.0, 300.05);
                alternate.validity = alternate.validity.with_quantity("T",
                    QuantitySpec::dimensional(Temperature::DIMS), 300.01, 2000.0);
                claims.insert_claim(alternate).unwrap();
            }
            claims.insert_claim(principal).unwrap();
        }
    }
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
        density: &PropertyKey::new("density", Density::DIMS), temperature_axis: "T",
        lower: &point(200.0), upper: &point(2000.0), selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 1e-6).unwrap();
    (card, curve)
}

fn cells<'a>(cards: &'a [MaterialCard; 2], curves: &'a [ResolvedHeatCapacityCurve; 2],
    contexts: &'a [QueryPoint; 2]) -> [MaterialCell<'a>; 2] {
    [0, 1].map(|i| MaterialCell {
        cell: DESTINATIONS[i].cell, material_id: MaterialId(i as u64), curve: curves[i].curve(),
        point: ThermalMaterialPoint { material: &cards[i], context: &contexts[i],
            temperature_axis: "T", liquid_fraction_axis: None }, conductivity: None,
    })
}

fn orientation(a: [f64; 3], b: [f64; 3]) -> [DirectionalConductionLink; 1] {
    [DirectionalConductionLink { link: 7, sides: [a, b].map(|direction| DirectionalConductionSide {
        principal_properties: ROLES.map(|name| PropertyKey::new(name, THERMAL_CONDUCTIVITY_DIMS)),
        selection: MaterialPropertySelection::SingleClaimOnly,
        direction_in_material_axes: direction, material_frame: ContentHash([17; 32]),
    }) }]
}

fn initial_links<'a>(curves: &'a [ResolvedHeatCapacityCurve; 2], extra: bool, reverse: bool)
    -> DepositionSimulation<'a> {
    let mut events = [0, 1].map(|i| TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
        id: DESTINATIONS[i].cell, segment_id: i as u64, body_id: i as u64,
        material_id: MaterialId(i as u64), mass_kg: 1.0,
        specific_h_j_kg: (curves[i].curve().knots().last().unwrap().specific_enthalpy_j_kg
            - curves[i].curve().knots()[0].specific_enthalpy_j_kg) / 18.0,
        geometry: CellGeometry { start_m: [i as f64, 0.0, 0.0], end_m: [i as f64 + 1.0, 0.0, 0.0],
            volume_m3: 1.0, surface_area_m2: 6.0 }, boundary: boundary(),
    }) }).to_vec();
    for id in 7..=if extra { 8 } else { 7 } {
        events.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
            id, a: if reverse { CellId(10) } else { CellId(30) },
            b: if reverse { CellId(30) } else { CellId(10) }, area_m2: 1.0, conductance_w_k: 0.0,
        }) });
    }
    DepositionSimulation::new(events, &[(MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
        (MaterialId(1), curves[1].curve() as &dyn ThermalLaw)], DepositionBudget {
        max_cells: 2, max_links: 2, max_events: 4, max_substeps: 100,
    }).unwrap().advance_to(0.0, StepConfig { max_dt_s: 0.1, energy_tolerance_j: 1e-7 }, || Ok(())).unwrap().0
}

#[test]
fn directional_sources_project_before_series_lowering_and_compose_external_links() {
    let (a, ca) = principal_material(1.0, false, false);
    let (b, cb) = principal_material(2.0, false, false);
    let cards = [a, b]; let curves = [ca, cb]; let contexts = [point(300.0), point(300.0)];
    assert!(cards.iter().all(|card| card.claims_for("thermal_conductivity").is_empty()));
    let cells = cells(&cards, &curves, &contexts); let links = links();
    let network = MaterialNetwork { cells: &cells, resistors: &[], links: &links };
    let thermal = initial_links(&curves, true, false);
    let diagonal = std::f64::consts::FRAC_1_SQRT_2;
    let directions = orientation([1.0, 0.0, 0.0], [0.0, diagonal, diagonal]);
    let additional = BTreeMap::from([(8, 0.75)]);
    let prepared = network.prepare_with_directional_conduction_and_link_conductances(
        &thermal, &[], &directions, &additional, || false).unwrap();
    let sample = prepared.resolve_with_link_conductances(&thermal, &additional, || false).unwrap();
    near(sample.source.links[0].conductivity_w_m_k[0], 1.0);
    near(sample.source.links[0].conductivity_w_m_k[1], 6.5);
    near(sample.conductances_w_k[&7], 1.0 / (0.5 + 0.5 / 6.5));
    assert_eq!(sample.conductances_w_k[&8], 0.75);
    for cell in &sample.source.cells {
        let snapshot = cell.conductivity.as_ref().unwrap();
        assert_eq!(snapshot.properties().len(), 3);
        for role in ROLES { assert!(snapshot.property(role).is_some()); }
    }
    let refreshed = prepared.resolve_with_link_conductances(&thermal, &BTreeMap::from([(8, 2.0)]), || false).unwrap();
    assert_eq!(refreshed.source, sample.source);
    assert_eq!(refreshed.conductances_w_k[&8], 2.0);
    // Reversing an authored gradient preserves its quadratic projection exactly.
    let mut reversed = directions.clone();
    for side in &mut reversed[0].sides { for x in &mut side.direction_in_material_axes { *x = -*x; } }
    let reverse = network.prepare_with_directional_conduction_and_link_conductances(
        &thermal, &[], &reversed, &additional, || false).unwrap();
    assert_eq!(reverse.resolve_with_link_conductances(&thermal, &additional, || false).unwrap().source, sample.source);
    // Rotation changes only data, not the original link geometry or material code.
    let rotated = orientation([0.0, 0.0, 1.0], [0.0, diagonal, diagonal]);
    let rotation = network.prepare_with_directional_conduction_and_link_conductances(
        &thermal, &[], &rotated, &additional, || false).unwrap();
    near(rotation.resolve_with_link_conductances(&thermal, &additional, || false).unwrap()
        .conductances_w_k[&7], 1.0 / (0.5 / 9.0 + 0.5 / 6.5));
    assert!(prepared.resolve(&thermal, || false).is_err());
    assert!(prepared.resolve_with_link_conductances(&thermal, &BTreeMap::from([(7, 2.0)]), || false).is_err());
    let mut stolen = directions.clone(); stolen[0].link = 8;
    assert!(network.prepare_with_directional_conduction_and_link_conductances(
        &thermal, &[], &stolen, &additional, || false).is_err());
    let mut unframed = directions.clone(); unframed[0].sides[0].material_frame = ContentHash([0; 32]);
    assert!(network.prepare_with_directional_conduction_and_link_conductances(
        &thermal, &[], &unframed, &additional, || false).is_err());
    let mut nonunit = directions.clone(); nonunit[0].sides[0].direction_in_material_axes = [2.0, 0.0, 0.0];
    assert!(network.prepare_with_directional_conduction_and_link_conductances(
        &thermal, &[], &nonunit, &additional, || false).is_err());
}

#[test]
fn directional_continuation_retains_every_claim_and_refuses_expiry_or_endpoint_migration() {
    let (a, ca) = principal_material(1.0, true, false);
    let (b, cb) = principal_material(2.0, false, false);
    let cards = [a, b]; let curves = [ca, cb]; let contexts = [point(300.0), point(300.0)];
    let cells = cells(&cards, &curves, &contexts); let links = links();
    let network = MaterialNetwork { cells: &cells, resistors: &[], links: &links };
    let thermal = initial(&curves, 0.0, Some(0.0));
    let directions = orientation([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    let prepared = network.prepare_with_directional_conduction(&thermal, &[], &directions, || false).unwrap();
    let before = prepared.resolve(&thermal, || false).unwrap().source;
    let observer = prepared.observe(()).unwrap();
    let cfg = StepConfig { max_dt_s: 0.005, energy_tolerance_j: 1e-7 };
    let (warm, _, audit) = thermal.advance_observed(&observer, 0.01, cfg, || Ok(())).unwrap();
    near(audit.stored_energy_change_j, 0.04);
    assert!(warm.temperature_k(CellId(30)).unwrap() > 300.01);
    assert!(network.prepare_with_directional_conduction(&warm, &[], &directions, || false).is_err());
    let mut polls = 0;
    let continued = prepared.extend_thermal(&network, &warm, || { polls += 1; false }).unwrap();
    let now = continued.resolve(&warm, || false).unwrap().source;
    for old in &before.cells {
        let current = now.cells.iter().find(|s| s.cell == old.cell).unwrap();
        near(current.phase.temperature_k(), warm.temperature_k(current.cell).unwrap());
        for role in ROLES {
            assert_eq!(old.conductivity.as_ref().unwrap().property(role).unwrap().answer().receipt.selected,
                current.conductivity.as_ref().unwrap().property(role).unwrap().answer().receipt.selected);
        }
    }
    assert!(now.links[0].conductance_w_k > before.links[0].conductance_w_k);
    let mut called = 0;
    assert!(matches!(prepared.extend_thermal(&network, &warm, || { called += 1; called == polls }),
        Err(MaterialNetworkError::Thermal(DepositionError::Cancelled))));
    // The thermal owner can construct a hotter state where only another claim
    // is eligible; extending the retained physical source must still refuse.
    let hot = thermal.advance_to(0.02, cfg, || Ok(())).unwrap().0;
    assert!(network.prepare_with_directional_conduction(&hot, &[], &directions, || false).is_ok());
    assert!(continued.extend_thermal(&network, &hot, || false).is_err());
    assert!(continued.resolve(&hot, || false).is_err());
    let migrated = initial_links(&curves, false, true);
    assert!(prepared.resolve(&migrated, || false).is_err());
    assert_eq!(prepared.resolve(&thermal, || false).unwrap().source, before);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn principal_conductivity_uses_the_actual_solved_electrical_field_bundle() {
    use fs_couple::material_circuit::{MaterialCircuit, conductivity::{ConductionAxis, ConductionCoordinate, ConductivitySolve}};
    use fs_phs::circuit::{Branch, CircuitGraph};
    let (a, ca) = principal_material(1.0, false, true);
    let (b, cb) = principal_material(2.0, false, false);
    let cards = [a, b]; let curves = [ca, cb];
    let contexts = [point(300.0).with_quantity("E", QuantitySpec::dimensional(FIELD), 0.0).unwrap(), point(300.0)];
    let cells = cells(&cards, &curves, &contexts); let links = links(); let resistors = branches();
    let network = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
    let thermal = initial(&curves, 999.0, None);
    let directions = orientation([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    let axes = [ConductionAxis { name: "E", coordinate: ConductionCoordinate::ElectricFieldMagnitude }];
    let responses = [CellConductionResponse { branch: 1, basis: ElectricalPropertyBasis::Resistivity, axes: &axes },
        CellConductionResponse { branch: 2, basis: ElectricalPropertyBasis::Resistivity, axes: &[] }];
    assert!(network.prepare_with_directional_conduction(&thermal,
        &[(1, ElectricalPropertyBasis::Resistivity), (2, ElectricalPropertyBasis::Resistivity)], &directions, || false).is_err());
    let prepared = network.prepare_conductivity_with_directional_conduction(&thermal, &responses, &directions, || false).unwrap();
    let seed = prepared.resolve_conductivity_seeds(&thermal, &responses, || false).unwrap();
    let circuit = MaterialCircuit::new(CircuitGraph { node_count: 3, branches: vec![
        (1, 0, Branch::Capacitor { farads: 1.0 }), (1, 2, Branch::Resistor { ohms: 999.0 }),
        (2, 0, Branch::Resistor { ohms: 999.0 })], transformers: vec![] },
        seed.resistors, vec![], vec![10.0], &[]).unwrap();
    let run = || prepared.conductivity_candidate(&thermal, 0.1, &circuit, &[], &responses,
        config(), ConductivitySolve { maximum_iterations: 64,
            absolute_resistance_tolerance_ohm: 1e-12, relative_resistance_tolerance: 1e-12 }, &[0.0; 2], || false).unwrap();
    let next = run(); let replay = run();
    assert_eq!(next.circuit(), replay.circuit());
    assert_eq!(next.thermal().cells(), replay.thermal().cells());
    let report = next.report();
    let source = report.source.cells.iter().find(|cell| cell.cell == CellId(30)).unwrap();
    let material = source.conductivity.as_ref().unwrap();
    let field = material.query().axes()["E"];
    assert!(field > 1.0);
    near(report.source.links[0].conductivity_w_m_k[0], 1.0 + 0.5 * field);
    near(material.property("k_longitudinal").unwrap().value_si(), 1.0 + 0.5 * field);
    near(material.query().axes()["T"], next.thermal().temperature_k(CellId(30)).unwrap());
    assert_eq!(material.properties().len(), 4);
    assert_eq!(material, report.electrical.evaluated[0].conductor.material());
    near(next.thermal().energy_j() + report.electrical.electrical.stored_after_j, thermal.energy_j() + 50.0);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(circuit.charges_c(), &[10.0]);
}
