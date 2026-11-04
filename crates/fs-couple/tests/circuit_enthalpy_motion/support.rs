use super::*;
use fs_couple::circuit_enthalpy::motion::support::{SupportThermalBinding, SupportThermalError};
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
    Event, StepConfig, ThermalLaw, TimedEvent};

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
    let events = (0..2).map(|i| {
        let layer = &law.reference().layer_integrals()[i];
        TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
            id: DESTINATIONS[i].cell, segment_id: i as u64, body_id: 0,
            material_id: MaterialId(i as u64), mass_kg: layer.line_density_kg_m(),
            specific_h_j_kg: curves[i].state_at_temperature(300.0).unwrap().specific_enthalpy_j_kg(),
            geometry: CellGeometry { start_m: [0.0; 3], end_m: [1.0, 0.0, 0.0],
                volume_m3: layer.area_m2(), surface_area_m2: 6.0 }, boundary: boundary(),
        }) }
    }).collect();
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
    let support = motion.support().unwrap();
    let strain = (motion.gap_m() - support.reference_gap_m()) / support.length_m();
    let member: f64 = support.section().section().layer_integrals().iter().enumerate().map(|(i, layer)| {
        let t = thermal.temperature_k(DESTINATIONS[i].cell).unwrap();
        let e = strain - alpha * (t - 300.0);
        let modulus = 1000.0 - 0.5 * (t - 300.0);
        let psi = 0.5 * modulus * e * e;
        let entropy = modulus * alpha * e + 0.25 * e * e;
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
