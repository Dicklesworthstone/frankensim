use super::*;
#[path = "../../../fs-conduction/tests/radiation_material/common.rs"]
mod radiation_common;
use fs_conduction::{AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as CONTACT_RD,
    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as CONTACT_RP, EMISSIVITY_DIMS,
    SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity};
use fs_conduction::transient::deposition::{
    contact::{ContactCell, MaterialContactLink}, radiation::{MaterialRadiation, RadiatingCell}};
use fs_couple::circuit_enthalpy::{CircuitThermalCoefficients,
    contacts::{ContactThermalCoefficients, heat::{HeatAssembly, HeatNetworkError,
        HeatNetworkSample, MaterialHeatNetwork, PreparedHeatNetwork}}};
use fs_exec::{CancelGate, Cx};
use fs_matdb::{ClaimSelection, InterfaceSystemCard};
use fs_couple::circuit_enthalpy::motion::support::{SupportThermalBinding, SupportThermalError,
    step_heat_assembly_candidate};
use fs_couple::material_circuit::motion::support::AxialSupport;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyValue,
    Provenance, SelectionPolicy, UncertaintyModel};
use fs_material::conductor::{ELECTRICAL_RESISTIVITY_DIMS, ELECTRICAL_RESISTIVITY_PROPERTY};
use fs_material::phase::{HeatCapacitySourceQuery, SolidLiquidPhase};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::section_temperature::{ConcentricSectionTemperatureLaw,
    thermoelastic::ConcentricThermoelasticLaw};
use fs_material::state_point::{DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as ALPHA, ScalarAdmissibility,
    ScalarPropertyRequirement, SPECIFIC_HEAT_CAPACITY_DIMS, YOUNG_MODULUS_PROPERTY,
    resolve_material_state_point};
use fs_qty::{Pressure, QuantitySpec, Temperature};
use fs_thermal::enthalpy::{CellGeometry, DepositedCell, DepositionBudget, DepositionSimulation,
    Event, StepConfig, ThermalLaw, ThermalLink, TimedEvent};

fn support_card(resistivity: f64, expansion: f64)
    -> (MaterialCard, ResolvedHeatCapacityCurve, Vec<ScalarPropertyRequirement>)
{
    support_card_with_capacity(resistivity, expansion, 10.0)
}

fn support_card_with_capacity(resistivity: f64, expansion: f64, capacity: f64)
    -> (MaterialCard, ResolvedHeatCapacityCurve, Vec<ScalarPropertyRequirement>)
{
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value, admissibility) in [
        (DENSITY_PROPERTY, Density::DIMS, 1.0, ScalarAdmissibility::StrictlyPositive),
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1000.0, ScalarAdmissibility::StrictlyPositive),
        (ALPHA, INVERSE_TEMPERATURE_DIMS, expansion, ScalarAdmissibility::Finite),
        ("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS, capacity, ScalarAdmissibility::StrictlyPositive),
        (ELECTRICAL_RESISTIVITY_PROPERTY, ELECTRICAL_RESISTIVITY_DIMS, resistivity,
            ScalarAdmissibility::StrictlyPositive),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 0.0, ScalarAdmissibility::NonNegative),
    ] {
        let affine = name == YOUNG_MODULUS_PROPERTY;
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: if affine { PropertyValue::Curve { abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS, knots: vec![(280.0, 1010.0), (400.0, 950.0)], dims }
            } else { PropertyValue::Scalar { value, dims } },
            validity: ValidityDomain::unconstrained().with_quantity("T",
                QuantitySpec::dimensional(Temperature::DIMS), 280.0, 400.0).with("lot", 7.0, 7.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: if affine { InterpolationPolicy::LinearInside }
                else { InterpolationPolicy::ConstantWithinValidity }, observations: vec![],
            provenance: Provenance { source: "synthetic thermally live actuator support".into(),
                license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name, dims, admissibility).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(),
        phase: "solid".into(), process: "authored thermoelastic support".into(), revision: 0 },
        claims, vec![]).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
        density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS), temperature_axis: "T",
        lower: &point(280.0), upper: &point(400.0), selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 1e-6).unwrap();
    (card, caloric, requirements)
}

fn law<'a>(cards: &'a [MaterialCard; 2], requirements: &[Vec<ScalarPropertyRequirement>; 2])
    -> ConcentricThermoelasticLaw<'a>
{
    let layers = (0..2).map(|i| ConcentricMaterialLayer {
        outer_radius_m: ((i + 1) as f64 / core::f64::consts::PI).sqrt(),
        material: resolve_material_state_point(&cards[i], &point(300.0), &requirements[i],
            MaterialPropertySelection::SingleClaimOnly).unwrap(),
    }).collect();
    ConcentricThermoelasticLaw::try_new(ConcentricSectionTemperatureLaw::try_new(
        ConcentricElasticSection::try_new(0.0, layers).unwrap(), &[&cards[0], &cards[1]], &["T", "T"],
    ).unwrap()).unwrap()
}

fn thermal<'a>(curves: &'a [ResolvedHeatCapacityCurve; 2], law: &ConcentricThermoelasticLaw<'_>)
    -> DepositionSimulation<'a>
{
    thermal_with_contact(curves, law, false)
}

fn thermal_with_contact<'a>(curves: &'a [ResolvedHeatCapacityCurve; 2],
    law: &ConcentricThermoelasticLaw<'_>, contact: bool) -> DepositionSimulation<'a>
{
    let mut events: Vec<_> = (0..2).map(|i| {
        let layer = &law.reference().layer_integrals()[i];
        TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
            id: DESTINATIONS[i].cell, segment_id: i as u64, body_id: 0,
            material_id: MaterialId(i as u64), mass_kg: layer.line_density_kg_m(),
            specific_h_j_kg: curves[i].state_at_temperature(300.0).unwrap().specific_enthalpy_j_kg(),
            geometry: CellGeometry { start_m: [0.0; 3], end_m: [1.0, 0.0, 0.0],
                volume_m3: layer.area_m2(), surface_area_m2: 6.0 }, boundary: boundary(),
        }) }
    }).collect();
    if contact {
        events.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
            id: 7, a: DESTINATIONS[0].cell, b: DESTINATIONS[1].cell,
            area_m2: 1.0, conductance_w_k: 0.0,
        }) });
    }
    DepositionSimulation::new(events, &[
        (MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
        (MaterialId(1), curves[1].curve() as &dyn ThermalLaw),
    ], DepositionBudget { max_cells: 8, max_links: 8, max_events: 16, max_substeps: 1000 }).unwrap()
        .advance_to(0.0, StepConfig { max_dt_s: 0.1, energy_tolerance_j: 1e-8 }, || Ok(())).unwrap().0
}

fn driven_circuit(resistors: Vec<MaterialResistor>) -> MaterialCircuit {
    // Independent voltage heats both layers while zero voltage holds magnetic
    // flux linkage. Only the common material support couples these subsystems.
    MaterialCircuit::new_with_storage(CircuitGraph { node_count: 3, branches: vec![
        (1, 0, Branch::Inductor { henries: f64::NAN }),
        (2, 0, Branch::Resistor { ohms: f64::NAN }),
        (2, 0, Branch::Resistor { ohms: f64::NAN }),
        (1, 0, Branch::VoltageSource { port: 0 }),
        (2, 0, Branch::VoltageSource { port: 1 })], transformers: vec![] },
        resistors, vec![gap_storage()], vec![0.4], vec![], &[0.0, 3.0]).unwrap()
}

fn total_internal(thermal: &DepositionSimulation<'_>, motion: &MaterialFieldMotion, alpha: f64) -> f64 {
    total_internal_layers(thermal, motion, [alpha; 2])
}

fn total_internal_layers(thermal: &DepositionSimulation<'_>, motion: &MaterialFieldMotion,
    alpha: [f64; 2]) -> f64 {
    let support = motion.support().unwrap();
    let strain = (motion.gap_m() - support.reference_gap_m()) / support.length_m();
    let member: f64 = support.section().section().layer_integrals().iter().enumerate().map(|(i, layer)| {
        let t = thermal.temperature_k(DESTINATIONS[i].cell).unwrap();
        let e = strain - alpha[i] * (t - 300.0);
        let modulus = 1000.0 - 0.5 * (t - 300.0);
        let psi = 0.5 * modulus * e * e;
        let entropy = modulus * alpha[i] * e + 0.25 * e * e;
        layer.area_m2() * support.length_m() * (psi + t * entropy)
    }).sum();
    thermal.energy_j() + 0.5 * motion.circuit().fluxes_wb()[0].powi(2) * motion.gap_m()
        + 0.5 * motion.generalized_mass_kg() * motion.velocity_m_s().powi(2) + member
}

#[test]
fn unequal_layer_heat_changes_the_same_actuator_and_closes_physical_internal_energy() {
    let run = |alpha| {
        let (a, ca, ra) = support_card(1.0, alpha); let (b, cb, rb) = support_card(3.0, alpha);
        let cards = [a, b]; let curves = [ca, cb]; let context = point(300.0);
        let binding = SupportThermalBinding { law: law(&cards, &[ra, rb]),
            layer_cells: DESTINATIONS.iter().map(|d| d.cell).collect(), volume_tolerance_m3: 1e-12 };
        let mut thermal = thermal(&curves, &binding.law);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = cells(&cards, &curves, &context, &key);
        let resistors = [0, 1].map(|i| CellResistor { area_m2: binding.law.reference().layer_integrals()[i].area_m2(),
            length_m: 1.0, destination: DESTINATIONS[i], selection: MaterialPropertySelection::SingleClaimOnly });
        let network = MaterialNetwork { cells: &cells, resistors: &resistors, links: &[] };
        let support = AxialSupport::from_thermoelastic_law(&binding.law, &[300.0; 2], 1.0, 1.0, 0.2).unwrap();
        let circuit = driven_circuit(network.resolve(&thermal, || false).unwrap().resistors);
        let mut motion = MaterialFieldMotion::new_with_support(circuit, 0, body(100.0), support,
            0.005, [0.9, 1.1]).unwrap();
        let initial = total_internal(&thermal, &motion, alpha);
        let mass = motion.generalized_mass_kg();
        let mut supplied = 0.0;
        for n in 1..=4 {
            let old = motion.clone();
            let next = network.field_motion_support_candidate(&thermal, f64::from(n) * 0.05,
                &motion, &[0.0, 3.0], 0.0, &bases(), &binding, config(), &[0.0; 2], || false).unwrap();
            let report = next.report();
            supplied += report.coupled.electrical.electrical_supply_j;
            near(total_internal(next.thermal(), next.motion(), alpha), initial + supplied);
            near(report.total_energy_residual_j, 0.0);
            near(report.external_heat_j, 0.0);
            assert_eq!(next.motion().generalized_mass_kg(), mass);
            near(next.motion().circuit().fluxes_wb()[0], 0.4);
            let temperatures = DESTINATIONS.map(|d| next.thermal().temperature_k(d.cell).unwrap());
            assert!(temperatures[0] > temperatures[1]);
            let expected_stiffness: f64 = binding.law.reference().layer_integrals().iter().zip(temperatures)
                .map(|(layer, t)| layer.area_m2() * (1000.0 - 0.5 * (t - 300.0))).sum();
            near(next.motion().support().unwrap().stiffness_n_m(), expected_stiffness);
            assert_ne!(next.motion().support().unwrap().section(), old.support().unwrap().section());
            for layer in &report.coupled.source.layers {
                near(layer.reversible_heat_j + layer.parameter_work_j
                    + layer.temperature_k[1] * layer.entropy_j_k[1]
                    - layer.temperature_k[0] * layer.entropy_j_k[0], 0.0);
            }
            let (next_thermal, next_motion, _) = next.into_parts(); thermal = next_thermal; motion = next_motion;
        }
        (motion.gap_m(), motion.velocity_m_s())
    };
    let expands = run(1e-3);
    let no_expansion = run(0.0);
    assert!((expands.0 - no_expansion.0).abs() > 1e-6);
    assert!((expands.1 - no_expansion.1).abs() > 1e-5);
}

#[test]
fn live_support_refuses_frozen_stepping_bad_inventory_and_cancellation_atomically() {
    let (a, ca, ra) = support_card(1.0, 1e-3); let (b, cb, rb) = support_card(3.0, 1e-3);
    let cards = [a, b]; let curves = [ca, cb]; let context = point(300.0);
    let binding = SupportThermalBinding { law: law(&cards, &[ra, rb]),
        layer_cells: DESTINATIONS.iter().map(|d| d.cell).collect(), volume_tolerance_m3: 1e-12 };
    let thermal = thermal(&curves, &binding.law);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let resistors = [0, 1].map(|i| CellResistor { area_m2: binding.law.reference().layer_integrals()[i].area_m2(),
        length_m: 1.0, destination: DESTINATIONS[i], selection: MaterialPropertySelection::SingleClaimOnly });
    let network = MaterialNetwork { cells: &cells, resistors: &resistors, links: &[] };
    let support = AxialSupport::from_thermoelastic_law(&binding.law, &[300.0; 2], 1.0, 1.0, 0.2).unwrap();
    let model = MaterialFieldMotion::new_with_support(
        driven_circuit(network.resolve(&thermal, || false).unwrap().resistors), 0, body(100.0),
        support, 0.005, [0.9, 1.1]).unwrap();
    let saved_cells = thermal.cells().clone(); let saved_model = model.clone();
    assert!(model.step_candidate(model.circuit().resistors().to_vec(), FieldMotionInput {
        electrical_ports: &[0.0, 3.0], external_force_n: 0.0, duration_s: 0.05 }, config().circuit, || false).is_err());
    let mut polls = 0;
    let normal = network.field_motion_support_candidate(&thermal, 0.05, &model, &[0.0, 3.0], 0.0,
        &bases(), &binding, config(), &[0.0; 2], || { polls += 1; false }).unwrap();
    let mut actual = 0;
    assert!(network.field_motion_support_candidate(&thermal, 0.05, &model, &[0.0, 3.0], 0.0,
        &bases(), &binding, config(), &[0.0; 2], || { actual += 1; actual == polls }).is_err());
    let invalid = SupportThermalBinding { law: binding.law.clone(), layer_cells: vec![CellId(30); 2],
        volume_tolerance_m3: binding.volume_tolerance_m3 };
    assert!(network.field_motion_support_candidate(&thermal, 0.05, &model, &[0.0, 3.0], 0.0,
        &bases(), &invalid, config(), &[0.0; 2], || false).is_err());
    assert_eq!(thermal.cells(), &saved_cells); assert_eq!(thermal.time_s(), 0.0); assert_eq!(model, saved_model);
    let retry = network.field_motion_support_candidate(&thermal, 0.05, &model, &[0.0, 3.0], 0.0,
        &bases(), &binding, config(), &[0.0; 2], || false).unwrap();
    assert_eq!(retry.motion(), normal.motion());
    assert_eq!(retry.thermal().cells(), normal.thermal().cells());
}

#[test]
fn nonpositive_combined_capacity_refuses_before_either_history_advances() {
    let (a, ca, ra) = support_card_with_capacity(1.0, 1e-3, 0.01);
    let (b, cb, rb) = support_card_with_capacity(3.0, 1e-3, 0.01);
    let cards = [a, b]; let curves = [ca, cb]; let context = point(300.0);
    let binding = SupportThermalBinding { law: law(&cards, &[ra, rb]),
        layer_cells: DESTINATIONS.iter().map(|d| d.cell).collect(), volume_tolerance_m3: 1e-12 };
    let thermal = thermal(&curves, &binding.law);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let resistors = [0, 1].map(|i| CellResistor { area_m2: binding.law.reference().layer_integrals()[i].area_m2(),
        length_m: 1.0, destination: DESTINATIONS[i], selection: MaterialPropertySelection::SingleClaimOnly });
    let network = MaterialNetwork { cells: &cells, resistors: &resistors, links: &[] };
    let support = AxialSupport::from_thermoelastic_law(&binding.law, &[300.0; 2], 1.0, 1.0, 0.2).unwrap();
    let model = MaterialFieldMotion::new_with_support(
        driven_circuit(network.resolve(&thermal, || false).unwrap().resistors), 0, body(100.0),
        support, 0.005, [0.9, 1.1]).unwrap();
    let saved_cells = thermal.cells().clone(); let saved_model = model.clone();
    assert!(matches!(network.field_motion_support_candidate(&thermal, 0.05, &model, &[0.0, 3.0], 0.0,
        &bases(), &binding, config(), &[0.0; 2], || false),
        Err(CircuitEnthalpyError::Source(SupportThermalError::Invalid(
            "caloric plus thermoelastic capacity must stay positive over the admitted support envelope")))));
    assert_eq!(thermal.cells(), &saved_cells); assert_eq!(thermal.time_s(), 0.0); assert_eq!(model, saved_model);
}

#[test]
fn contact_heat_drives_live_support_expansion_with_one_closed_energy_account() {
    let run = |contact_resistance| {
        let alpha = 1e-3;
        let (a, ca, ra) = support_card(1.0, alpha);
        let (b, cb, rb) = support_card(3.0, alpha);
        let cards = [a, b]; let curves = [ca, cb]; let context = point(300.0);
        let binding = SupportThermalBinding { law: law(&cards, &[ra, rb]),
            layer_cells: DESTINATIONS.iter().map(|d| d.cell).collect(), volume_tolerance_m3: 1e-12 };
        let mut thermal = thermal_with_contact(&curves, &binding.law, true);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = cells(&cards, &curves, &context, &key);
        let resistors = [0, 1].map(|i| CellResistor {
            area_m2: binding.law.reference().layer_integrals()[i].area_m2(), length_m: 1.0,
            destination: DESTINATIONS[i], selection: MaterialPropertySelection::SingleClaimOnly });
        let links = links(); let bases = bases();
        let bulk = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
        let card = source::card(&cards[0], &cards[1], contact_resistance, 0.005);
        let interface_context = source::point();
        let bindings = [MaterialElectricalContact { link: 7, branch: 1,
            point: ThermalInterfacePoint { interface: &card, materials: [&cards[0], &cards[1]],
                context: &interface_context, temperature_axes: ["Ta", "Tb"],
                liquid_fraction_axes: [Some("melt_a"), Some("melt_b")], mean_temperature_axis: None },
            selection: MaterialPropertySelection::SingleClaimOnly, heat_fraction_to_first: 0.25 }];
        let network = ElectricalContactNetwork { bulk: &bulk, contacts: &bindings, electrical_bases: &bases };
        let original = network.resolve(&thermal, || false).unwrap();
        let circuit = driven_circuit(original.bulk.resistors)
            .rebind_contacts_candidate(original.contacts, &[0.0, 3.0]).unwrap();
        let support = AxialSupport::from_thermoelastic_law(&binding.law, &[300.0; 2], 1.0, 1.0, 0.2).unwrap();
        let mut motion = MaterialFieldMotion::new_with_support(circuit, 0, body(100.0),
            support, 0.005, [0.9, 1.1]).unwrap();
        let initial = total_internal(&thermal, &motion, alpha);
        let mut supplied = 0.0;
        for n in 1..=4 {
            let target = f64::from(n) * 0.05;
            let dt = target - thermal.time_s();
            let candidate = network.field_motion_support_candidate(&thermal, target, &motion,
                &[0.0, 3.0], 0.0, &binding, config(), &[0.0; 2], || false).unwrap();
            let report = candidate.report();
            let electrical = &report.coupled.electrical;
            let contact = &electrical.contacts[0];
            let temperature = candidate.thermal().temperature_k(DESTINATIONS[0].cell).unwrap();
            let resistance = contact_resistance * (1.0 + 0.005 * (temperature - 300.0));
            near(contact.contact.resistance_ohm(), resistance);
            // Independent series branch law at the actual thermal endpoint.
            let bulk_resistance = electrical.resistors[0].material.conductor.resistance_ohm();
            let current = 3.0 / (bulk_resistance + resistance);
            near(contact.current_a.abs(), current);
            near(contact.heat_j, dt * current * current * resistance);
            for (i, destination) in DESTINATIONS.iter().enumerate() {
                let old = &thermal.cells()[&destination.cell];
                let new = &candidate.thermal().cells()[&destination.cell];
                let heat = old.mass_kg * (new.specific_h_j_kg - old.specific_h_j_kg);
                let fraction = if i == 0 { 0.25 } else { 0.75 };
                let bulk_heat = electrical.resistors.iter()
                    .find(|r| r.material.branch == destination.branch).unwrap().heat_j;
                let reversible = report.coupled.source.layers[i].reversible_heat_j;
                near(heat, bulk_heat + fraction * contact.heat_j + reversible);
            }
            near(report.coupled.deposited_joule_j, electrical.heat_j);
            near(report.external_heat_j, 0.0);
            near(report.total_energy_residual_j, 0.0);
            assert_eq!(report.coupled.source.contact_destinations[0].cells,
                [DESTINATIONS[0].cell, DESTINATIONS[1].cell]);
            near(report.coupled.source.contact_destinations[0].heat_fraction_to_first, 0.25);
            supplied += electrical.electrical_supply_j;
            near(total_internal(candidate.thermal(), candidate.motion(), alpha), initial + supplied);
            let (next_thermal, next_motion, _) = candidate.into_parts();
            thermal = next_thermal; motion = next_motion;
        }
        (motion.gap_m(), motion.velocity_m_s(), thermal.temperature_k(DESTINATIONS[0].cell).unwrap())
    };
    let low = run(0.2);
    let high = run(2.0);
    assert!(high.2 < low.2, "series contact resistance must change actual layer heating");
    assert!((low.0 - high.0).abs() > 1e-6, "contact heating must reach support geometry");
    assert!((low.1 - high.1).abs() > 1e-5, "support response must reach native momentum");
}

#[test]
fn live_support_contact_source_budget_and_late_cancellation_preserve_all_owners() {
    let (a, ca, ra) = support_card(1.0, 1e-3);
    let (b, cb, rb) = support_card(3.0, 1e-3);
    let cards = [a, b]; let curves = [ca, cb]; let context = point(300.0);
    let binding = SupportThermalBinding { law: law(&cards, &[ra, rb]),
        layer_cells: DESTINATIONS.iter().map(|d| d.cell).collect(), volume_tolerance_m3: 1e-12 };
    let thermal = thermal_with_contact(&curves, &binding.law, true);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let resistors = [0, 1].map(|i| CellResistor {
        area_m2: binding.law.reference().layer_integrals()[i].area_m2(), length_m: 1.0,
        destination: DESTINATIONS[i], selection: MaterialPropertySelection::SingleClaimOnly });
    let links = links(); let bases = bases();
    let bulk = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
    let card = source::card(&cards[0], &cards[1], 0.2, 0.005);
    let interface_context = source::point();
    let bindings = [MaterialElectricalContact { link: 7, branch: 1,
        point: ThermalInterfacePoint { interface: &card, materials: [&cards[0], &cards[1]],
            context: &interface_context, temperature_axes: ["Ta", "Tb"],
            liquid_fraction_axes: [Some("melt_a"), Some("melt_b")], mean_temperature_axis: None },
        selection: MaterialPropertySelection::SingleClaimOnly, heat_fraction_to_first: 0.25 }];
    let network = ElectricalContactNetwork { bulk: &bulk, contacts: &bindings, electrical_bases: &bases };
    let original = network.resolve(&thermal, || false).unwrap();
    let circuit = driven_circuit(original.bulk.resistors)
        .rebind_contacts_candidate(original.contacts, &[0.0, 3.0]).unwrap();
    let support = AxialSupport::from_thermoelastic_law(&binding.law, &[300.0; 2], 1.0, 1.0, 0.2).unwrap();
    let model = MaterialFieldMotion::new_with_support(circuit, 0, body(100.0),
        support, 0.005, [0.9, 1.1]).unwrap();
    let saved_cells = thermal.cells().clone(); let saved_model = model.clone();
    let mut polls = 0;
    let normal = network.field_motion_support_candidate(&thermal, 0.05, &model,
        &[0.0, 3.0], 0.0, &binding, config(), &[0.0; 2], || { polls += 1; false }).unwrap();
    let mut actual = 0;
    assert!(matches!(network.field_motion_support_candidate(&thermal, 0.05, &model,
        &[0.0, 3.0], 0.0, &binding, config(), &[0.0; 2], || { actual += 1; actual == polls }),
        Err(CircuitEnthalpyError::Cancelled)));
    let limited = CircuitEnthalpySolve { maximum_circuit_trials: normal.report().coupled.circuit_trials - 1,
        ..config() };
    assert!(matches!(network.field_motion_support_candidate(&thermal, 0.05, &model,
        &[0.0, 3.0], 0.0, &binding, limited, &[0.0; 2], || false),
        Err(CircuitEnthalpyError::CircuitBudget)));
    let foreign_card = source::card(&cards[0], &cards[1], 0.4, 0.005);
    let foreign_bindings = [MaterialElectricalContact { link: 7, branch: 1,
        point: ThermalInterfacePoint { interface: &foreign_card, materials: [&cards[0], &cards[1]],
            context: &interface_context, temperature_axes: ["Ta", "Tb"],
            liquid_fraction_axes: [Some("melt_a"), Some("melt_b")], mean_temperature_axis: None },
        selection: MaterialPropertySelection::SingleClaimOnly, heat_fraction_to_first: 0.25 }];
    let foreign = ElectricalContactNetwork { bulk: &bulk, contacts: &foreign_bindings, electrical_bases: &bases };
    assert!(matches!(foreign.field_motion_support_candidate(&thermal, 0.05, &model,
        &[0.0, 3.0], 0.0, &binding, config(), &[0.0; 2], || false),
        Err(CircuitEnthalpyError::Invalid("electrical contact source changed during thermal feedback"))));
    assert!(bulk.field_motion_support_candidate(&thermal, 0.05, &model,
        &[0.0, 3.0], 0.0, &bases, &binding, config(), &[0.0; 2], || false).is_err());
    assert_eq!(thermal.cells(), &saved_cells); assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(model, saved_model);
    let retry = network.field_motion_support_candidate(&thermal, 0.05, &model,
        &[0.0, 3.0], 0.0, &binding, config(), &[0.0; 2], || false).unwrap();
    assert_eq!(retry.motion(), normal.motion());
    assert_eq!(retry.thermal().cells(), normal.thermal().cells());
    assert_eq!(retry.report().coupled.source, normal.report().coupled.source);
}

// The existing support fixture gains sourced transport, without changing the
// original mechanics, caloric law or reference inventory.
fn transport_card(resistivity: f64, expansion: f64, emissivity: f64)
    -> (MaterialCard, ResolvedHeatCapacityCurve, Vec<ScalarPropertyRequirement>)
{
    let (original, _, requirements) = support_card(resistivity, expansion);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == "thermal_conductivity" {
            claim.value = PropertyValue::Scalar { value: 4.0, dims: THERMAL_CONDUCTIVITY_DIMS };
        }
        claims.insert_claim(claim).unwrap();
    }
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
        value: PropertyValue::Scalar { value: emissivity, dims: EMISSIVITY_DIMS },
        validity: ValidityDomain::unconstrained().with_quantity("T",
            QuantitySpec::dimensional(Temperature::DIMS), 280.0, 400.0).with("lot", 7.0, 7.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
        provenance: Provenance { source: "synthetic support emissivity".into(),
            license: "CC0-1.0".into(), artifact: None },
    }).unwrap();
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
        density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS), temperature_axis: "T",
        lower: &point(280.0), upper: &point(400.0), selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 1e-6).unwrap();
    (card, caloric, requirements)
}

fn support_interface(cards: &[MaterialCard; 2], resistance: f64) -> InterfaceSystemCard {
    let original = source::card(&cards[0], &cards[1], 0.2, 0.005);
    let mut claims = original.claims().clone();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(CONTACT_RP, CONTACT_RD),
        value: PropertyValue::Curve { abscissa: "Ta".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(280.0, 0.98 * resistance), (400.0, 1.1 * resistance)], dims: CONTACT_RD },
        validity: ValidityDomain::unconstrained().with("Ta", 280.0, 400.0)
            .with("Tb", 280.0, 400.0).with("pressure", 100000.0, 100000.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: Provenance { source: "synthetic support thermal interface".into(),
            license: "CC0-1.0".into(), artifact: None },
    }).unwrap();
    InterfaceSystemCard::assemble(original.surface_a().clone(), original.surface_b().clone(),
        original.context().clone(), claims, vec![]).unwrap()
}

fn with_heat_support<R>(emissivity: f64, resistance: f64,
    run: impl FnOnce(&Cx<'_>, &CancelGate, &MaterialHeatNetwork<'_>, &SupportThermalBinding<'_>,
        DepositionSimulation<'_>, MaterialFieldMotion) -> R) -> R
{
    let (a, ca, ra) = transport_card(1.0, 1e-3, emissivity);
    let (b, cb, rb) = transport_card(3.0, 2e-3, 0.6);
    let cards = [a, b]; let curves = [ca, cb]; let context = point(300.0);
    let binding = SupportThermalBinding { law: law(&cards, &[ra, rb]),
        layer_cells: DESTINATIONS.iter().map(|d| d.cell).collect(), volume_tolerance_m3: 1e-12 };
    let original = thermal_with_contact(&curves, &binding.law, true);
    let mut events = original.cells().values().map(|cell| {
        let mut cell = cell.clone();
        // Prescribed removal belongs to the original caloric boundary owner.
        // It is distinct from reciprocal conduction and closed radiation.
        cell.boundary.source_w = -1.0;
        TimedEvent { time_s: 0.0, event: Event::Birth(cell) }
    }).collect::<Vec<_>>();
    events.extend(original.links().values().map(|link| TimedEvent {
        time_s: 0.0, event: Event::Contact(link.clone()) }));
    let thermal = DepositionSimulation::new(events, &[
        (MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
        (MaterialId(1), curves[1].curve() as &dyn ThermalLaw),
    ], DepositionBudget { max_cells: 8, max_links: 8, max_events: 16, max_substeps: 1000 }).unwrap()
        .advance_to(0.0, StepConfig { max_dt_s: 0.1, energy_tolerance_j: 1e-8 }, || Ok(())).unwrap().0;
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let resistors = [0, 1].map(|i| CellResistor {
        area_m2: binding.law.reference().layer_integrals()[i].area_m2(), length_m: 1.0,
        destination: DESTINATIONS[i], selection: MaterialPropertySelection::SingleClaimOnly });
    let links = links(); let bases = bases();
    let bulk = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
    let card = support_interface(&cards, resistance);
    let interface_context = source::point();
    let interface = ThermalInterfacePoint { interface: &card, materials: [&cards[0], &cards[1]],
        context: &interface_context, temperature_axes: ["Ta", "Tb"],
        liquid_fraction_axes: [Some("melt_a"), Some("melt_b")], mean_temperature_axis: None };
    let electrical_contacts = [MaterialElectricalContact { link: 7, branch: 1, point: interface,
        selection: MaterialPropertySelection::SingleClaimOnly, heat_fraction_to_first: 0.25 }];
    let electrical = ElectricalContactNetwork { bulk: &bulk, contacts: &electrical_contacts,
        electrical_bases: &bases };
    let thermal_contacts = [MaterialContactLink { link: 7, point: interface,
        sides: [0, 1].map(|i| ContactCell { cell: DESTINATIONS[i].cell,
            material_id: MaterialId(i as u64), curve: curves[i].curve() }),
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly) }];
    let surfaces = [0, 1].map(|i| RadiatingCell { surface: i, cell: DESTINATIONS[i].cell,
        material_id: MaterialId(i as u64), curve: curves[i].curve(), material: cells[i].point,
        selection: MaterialPropertySelection::SingleClaimOnly });
    let enclosure = radiation_common::enclosure([0, 1].map(|i|
        SurfaceEmissivity::from_thermal_state("support enclosure", cells[i].point,
            curves[i].state_at_temperature(300.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly).unwrap().0));
    let radiation = MaterialRadiation { enclosure: &enclosure, cells: &surfaces,
        power_tolerance_w: 1e-8, radiosity_tolerance_w_m2: 1e-8 };
    let assembly = MaterialHeatNetwork { electrical: &electrical,
        thermal_contacts: &thermal_contacts, radiation: Some(&radiation) };
    let resolved = electrical.resolve(&thermal, || false).unwrap();
    let circuit = driven_circuit(resolved.bulk.resistors)
        .rebind_contacts_candidate(resolved.contacts, &[0.0, 3.0]).unwrap();
    let support = AxialSupport::from_thermoelastic_law(&binding.law, &[300.0; 2],
        1.0, 1.0, 0.2).unwrap();
    let motion = MaterialFieldMotion::new_with_support(circuit, 0, body(100.0),
        support, 0.005, [0.9, 1.1]).unwrap();
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena, fs_exec::StreamKey { seed: 73, kernel_id: 51, tile: 0, iteration: 0 },
            fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic);
        run(&cx, &gate, &assembly, &binding, thermal, motion)
    })
}

#[test]
fn complete_heat_assembly_changes_live_support_with_one_physical_energy_account() {
    let run = |emissivity, resistance| with_heat_support(emissivity, resistance,
        |cx, _, assembly, binding, mut thermal, mut motion| {
        let prepared = assembly.prepare(cx, &thermal).unwrap();
        let initial = total_internal_layers(&thermal, &motion, [1e-3, 2e-3]);
        let mut supplied = 0.0;
        for n in 1..=4 {
            let target = f64::from(n) * 0.05;
            let dt = target - thermal.time_s();
            let next = prepared.field_motion_support_candidate(cx, &thermal, target,
                &motion, &[0.0, 3.0], 0.0, binding, config(), &[0.0; 2]).unwrap();
            let report = next.report();
            let network = &report.coupled.source.network;
            let contact = &network.thermal_contacts[0];
            let temperatures = DESTINATIONS.map(|d| next.thermal().temperature_k(d.cell).unwrap());
            near(contact.interface.value_m2_k_per_w(),
                resistance * (1.0 + 0.001 * (temperatures[0] - 300.0)));
            near(contact.conductance_w_k, 1.0 / (0.25
                + contact.interface.value_m2_k_per_w() / contact.area_m2));
            let radiation = network.radiation.as_ref().unwrap();
            let plate_power = 5.670374419e-8 * (temperatures[0].powi(4) - temperatures[1].powi(4))
                / (1.0 / emissivity + 1.0 / 0.6 - 1.0);
            near(radiation.sources_w()[&DESTINATIONS[0].cell], -plate_power);
            near(radiation.sources_w()[&DESTINATIONS[1].cell], plate_power);
            let electrical = &report.coupled.electrical;
            for (i, destination) in DESTINATIONS.iter().enumerate() {
                let old = &thermal.cells()[&destination.cell];
                let endpoint = &next.thermal().cells()[&destination.cell];
                let caloric = old.mass_kg * (endpoint.specific_h_j_kg - old.specific_h_j_kg);
                let bulk_heat = electrical.resistors.iter()
                    .find(|r| r.material.branch == destination.branch).unwrap().heat_j;
                let interface_heat = electrical.contacts[0].heat_j * if i == 0 { 0.25 } else { 0.75 };
                let transport = contact.conductance_w_k * (temperatures[1 - i] - temperatures[i])
                    + radiation.sources_w()[&destination.cell] - 1.0;
                near(caloric, bulk_heat + interface_heat
                    + report.coupled.source.layers[i].reversible_heat_j + dt * transport);
                near(contact.phases[i].temperature_k(), temperatures[i]);
            }
            near(report.coupled.deposited_joule_j, electrical.heat_j);
            near(report.external_heat_j, -2.0 * dt);
            near(report.external_heat_j, report.coupled.non_electrical_heat_j
                - report.coupled.source.reversible_heat_j - dt * network.radiation_net_power_w());
            near(report.total_energy_residual_j, 0.0);
            supplied += electrical.electrical_supply_j + report.external_heat_j;
            near(total_internal_layers(next.thermal(), next.motion(), [1e-3, 2e-3]),
                initial + supplied);
            let (next_thermal, next_motion, _) = next.into_parts();
            thermal = next_thermal; motion = next_motion;
        }
        [thermal.temperature_k(DESTINATIONS[0].cell).unwrap(), motion.gap_m(), motion.velocity_m_s()]
    });
    let baseline = run(0.2, 0.2);
    for changed in [run(0.8, 0.2), run(0.2, 2.0)] {
        for i in 0..3 {
            assert!((changed[i] - baseline[i]).abs() > 1e-9,
                "sourced thermal transport must reach temperature, gap and momentum");
        }
    }
}

#[test]
fn empty_heat_assembly_reproduces_the_existing_contact_support_candidate() {
    with_heat_support(0.2, 0.2, |cx, _, assembly, binding, thermal, motion| {
        let empty = MaterialHeatNetwork { electrical: assembly.electrical,
            thermal_contacts: &[], radiation: None };
        let prepared = empty.prepare(cx, &thermal).unwrap();
        let actual = prepared.field_motion_support_candidate(cx, &thermal, 0.05,
            &motion, &[0.0, 3.0], 0.0, binding, config(), &[0.0; 2]).unwrap();
        let original = assembly.electrical.field_motion_support_candidate(&thermal, 0.05,
            &motion, &[0.0, 3.0], 0.0, binding, config(), &[0.0; 2], || false).unwrap();
        assert_eq!(actual.thermal().cells(), original.thermal().cells());
        assert_eq!(actual.motion(), original.motion());
        assert_eq!(actual.report().coupled.source.network.material, original.report().coupled.source.network);
        assert_eq!(actual.report().coupled.source.layers, original.report().coupled.source.layers);
        assert_eq!(actual.report().external_heat_j, original.report().external_heat_j);
    });
}

#[derive(Debug)]
enum SupportAssemblyFailure {
    Heat(HeatNetworkError),
    Convection(fs_convection::CorrelationError),
    Refused,
}
impl From<HeatNetworkError> for SupportAssemblyFailure {
    fn from(error: HeatNetworkError) -> Self { Self::Heat(error) }
}
impl core::fmt::Display for SupportAssemblyFailure {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Heat(error) => write!(f, "{error}"),
            Self::Convection(error) => write!(f, "{error}"),
            Self::Refused => write!(f, "injected late complete-assembly refusal"),
        }
    }
}
struct SupportAssemblySample {
    heat: HeatNetworkSample,
    temperature_k: [f64; 2],
    reservoir_conductance_w_k: f64,
    reservoir_heat_w: f64,
}
impl AsRef<HeatNetworkSample> for SupportAssemblySample {
    fn as_ref(&self) -> &HeatNetworkSample { &self.heat }
}

// Compose the same provider seam used by downstream fs-airflow without adding
// a dependency cycle. Native fs-convection evaluates the prescribed-fluid
// coefficient on a disjoint exposed patch; the original heat assembly resolves
// every other path. Instrumentation advances no physical history.
struct SupportAssemblyProbe<'p, 'n, 'a> {
    prepared: &'p PreparedHeatNetwork<'n, 'a>,
    gate: &'p CancelGate,
    calls: std::cell::Cell<usize>,
    refuse_at: std::cell::Cell<Option<usize>>,
    cancel_at: std::cell::Cell<Option<usize>>,
}
impl HeatAssembly for SupportAssemblyProbe<'_, '_, '_> {
    type Sample = SupportAssemblySample;
    type Error = SupportAssemblyFailure;
    fn prepared_heat_network(&self) -> &PreparedHeatNetwork<'_, '_> { self.prepared }
    fn resolve_heat(&self, cx: &Cx<'_>, thermal: &DepositionSimulation<'_>)
        -> Result<ContactThermalCoefficients<Self::Sample>, CircuitEnthalpyError<Self::Error>>
    {
        let mut resolved = self.prepared.resolve(cx, thermal).map_err(|error| match error {
            HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
            other => CircuitEnthalpyError::Source(SupportAssemblyFailure::Heat(other)),
        })?;
        let call = self.calls.get() + 1;
        self.calls.set(call);
        if self.refuse_at.get() == Some(call) {
            return Err(CircuitEnthalpyError::Source(SupportAssemblyFailure::Refused));
        }
        let temperature_k = [thermal.temperature_k(DESTINATIONS[0].cell)?,
            thermal.temperature_k(DESTINATIONS[1].cell)?];
        let nu = fs_convection::evaluate(fs_convection::CorrelationId::FlatPlateLaminarAverage,
            fs_convection::CorrelationInputs::forced(10_000.0, 20.0 / 3.0))
            .map_err(|error| CircuitEnthalpyError::Source(SupportAssemblyFailure::Convection(error)))?;
        let reservoir_conductance_w_k = 1e-4 * nu.heat_transfer_coefficient(
            fs_convection::ThermalConductivity::new(0.6), fs_qty::Length::new(0.1))
            .map_err(|error| CircuitEnthalpyError::Source(SupportAssemblyFailure::Convection(error)))?
            .value.value();
        let reservoir_heat_w = reservoir_conductance_w_k * (290.0 - temperature_k[0]);
        *resolved.bulk.sources_w.entry(DESTINATIONS[0].cell).or_default() += reservoir_heat_w;
        if self.cancel_at.get() == Some(call) { self.gate.request(); }
        Ok(ContactThermalCoefficients { contacts: resolved.contacts,
            bulk: CircuitThermalCoefficients { resistors: resolved.bulk.resistors,
                sources_w: resolved.bulk.sources_w, conductances_w_k: resolved.bulk.conductances_w_k,
                source: SupportAssemblySample { heat: resolved.bulk.source, temperature_k,
                    reservoir_conductance_w_k, reservoir_heat_w } } })
    }
}

#[test]
fn complete_support_assembly_retains_typed_samples_and_refuses_atomically() {
    with_heat_support(0.2, 0.2, |cx, gate, assembly, binding, thermal, motion| {
        let prepared = assembly.prepare(cx, &thermal).unwrap();
        let probe = SupportAssemblyProbe { prepared: &prepared, gate,
            calls: std::cell::Cell::new(0), refuse_at: std::cell::Cell::new(None),
            cancel_at: std::cell::Cell::new(None) };
        let run = |controls| step_heat_assembly_candidate(&probe, cx, &thermal, 0.05,
            &motion, &[0.0, 3.0], 0.0, binding, controls, &[0.0; 2]);
        let saved_cells = thermal.cells().clone();
        let saved_motion = motion.clone();
        let normal = run(config()).unwrap();
        let sample = &normal.report().coupled.source.network;
        assert!(sample.reservoir_heat_w < 0.0);
        near(sample.reservoir_heat_w,
            sample.reservoir_conductance_w_k * (290.0 - sample.temperature_k[0]));
        near(normal.report().external_heat_j, 0.05 * (sample.reservoir_heat_w - 2.0));
        near(normal.report().external_heat_j, normal.report().coupled.non_electrical_heat_j
            - normal.report().coupled.source.reversible_heat_j
            - 0.05 * sample.heat.radiation_net_power_w());
        near(total_internal_layers(normal.thermal(), normal.motion(), [1e-3, 2e-3]),
            total_internal_layers(&thermal, &motion, [1e-3, 2e-3])
                + normal.report().coupled.electrical.electrical_supply_j
                + normal.report().external_heat_j);
        let calls = probe.calls.get();
        assert!(calls > 1);
        assert_eq!(normal.report().coupled.source.network.temperature_k,
            DESTINATIONS.map(|d| normal.thermal().temperature_k(d.cell).unwrap()));
        probe.calls.set(0);
        probe.refuse_at.set(Some(calls));
        assert!(matches!(run(config()),
            Err(CircuitEnthalpyError::Source(SupportThermalError::Heat(SupportAssemblyFailure::Refused)))));
        probe.calls.set(0);
        probe.refuse_at.set(None);
        assert!(matches!(run(CircuitEnthalpySolve {
            maximum_circuit_trials: normal.report().coupled.circuit_trials - 1, ..config() }),
            Err(CircuitEnthalpyError::CircuitBudget)));
        probe.calls.set(0);
        let retry = run(config()).unwrap();
        assert_eq!(retry.thermal().cells(), normal.thermal().cells());
        assert_eq!(retry.motion(), normal.motion());
        probe.calls.set(0);
        probe.cancel_at.set(Some(calls));
        assert!(matches!(run(config()), Err(CircuitEnthalpyError::Cancelled)));
        assert_eq!(thermal.cells(), &saved_cells);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(motion, saved_motion);
    });
}
