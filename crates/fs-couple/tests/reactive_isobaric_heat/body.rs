use super::*;
#[path = "body/adaptive.rs"]
mod adaptive;
use fs_blake3::ContentHash;
use fs_couple::cavity_phs::piston::reactive::implicit::coupled::isobaric::thermal::body::IsobaricBodySolve;
use fs_couple::cavity_phs::piston::thermal::layers::MaterialConductionLayers;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve,
    HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase, UniformEnthalpyStepInput};
use fs_material::state_point::{MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS};
use fs_matdb::{ClaimSet, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, InterpolationPolicy, UncertaintyModel, SelectionPolicy};
use fs_evidence::ValidityDomain;
use fs_qty::{Density, Temperature};

fn limits() -> IsobaricBodySolve {
    IsobaricBodySolve { heat: exchange(), entropy_tolerance_j_k: 1e-8 }
}
fn sensible() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(ContentHash([0x41; 32]),
        SolidLiquidPhase::Solid, vec![
            EnthalpyPhaseKnot { specific_enthalpy_j_kg: 0.0, temperature_k: 200.0,
                liquid_mass_fraction: 0.0, bulk_density_kg_m3: 4000.0 },
            EnthalpyPhaseKnot { specific_enthalpy_j_kg: 400000.0, temperature_k: 1000.0,
                liquid_mass_fraction: 0.0, bulk_density_kg_m3: 4000.0 },
        ]).unwrap()
}
fn body(curve: &EquilibriumEnthalpyPhaseCurve, h: f64, dt: f64) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput { curve, initial: curve.state_at_specific_enthalpy(h).unwrap(),
        mass_kg: 0.01, volume_m3: 0.01/4000.0, surface_area_m2: 0.02,
        internal_heat_j: 0.0, duration_s: dt }
}

#[test]
fn finite_body_and_gas_match_independent_capacity_heat_work_and_entropy() {
    let law = law(0.0);
    let old = gas(&law, [0.5, 0.5], 400.0);
    let curve = sensible();
    let input = UniformEnthalpyStepInput { internal_heat_j: 3.0,
        ..body(&curve, 250000.0, 0.1) };
    let old_body = input.initial;
    let result = law.isobaric_thermal_body_candidate(&old, 0.1, 2.0, input, controls(), limits(),
        |gas, state| Ok::<_, Infallible>(ConductanceSample { conductance_w_k: 0.4,
            source: (gas.state().temperature, state) }), || false).unwrap();
    let cg = old.mass_kg()/0.028*3.5*R;
    let cb = input.mass_kg*500.0;
    let q = 0.1*0.4*(700.0-400.0+3.0/cb-2.0/cg)/(1.0+0.1*0.4*(1.0/cg+1.0/cb));
    let tg = 400.0+(2.0+q)/cg;
    let tb = 700.0+(3.0-q)/cb;
    let report = result.step();
    near(report.boundary.exchanged_heat_j, q, 2e-8);
    near(result.gas().state().state().temperature, tg, 2e-7);
    near(result.thermal().state.temperature_k(), tb, 2e-7);
    assert_eq!(result.gas().state().state().pressure, 100000.0);
    assert_eq!(result.gas().mass_kg(), old.mass_kg());
    assert_eq!(result.thermal().external_heat_j, -report.boundary.exchanged_heat_j);
    assert_eq!(result.thermal().report, (result.gas().state().state().temperature, result.thermal().state));
    near(report.thermal_energy_change_j, 3.0-q, 1e-7);
    near(report.boundary.chemistry.boundary.boundary.work_by_gas_j,
        old.mass_kg()/0.028*R*(tg-400.0), 1e-7);
    assert!(report.total_enthalpy_residual_j.abs() <= limits().heat.exchange.energy_tolerance_j);
    assert!(report.total_energy_residual_j.abs() <= limits().heat.exchange.energy_tolerance_j);
    near(report.thermal_entropy_change_j_k, cb*(tb/700.0).ln(), 2e-9);
    let entropy = cg*(tg/400.0).ln()+cb*(tb/700.0).ln()-2.0/tg-3.0/tb;
    near(report.total_entropy_residual_j_k, entropy, 2e-9);
    assert!(report.total_entropy_residual_j_k > 0.0);
    assert_eq!(input.initial, old_body);
}

#[test]
fn latent_wall_changes_fraction_without_an_apparent_heat_capacity() {
    let law = law(0.0);
    let old = gas(&law, [0.5, 0.5], 400.0);
    let knot = |h, t, f| EnthalpyPhaseKnot { specific_enthalpy_j_kg: h,
        temperature_k: t, liquid_mass_fraction: f, bulk_density_kg_m3: 4000.0 };
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x42; 32]),
        vec![knot(0.0, 300.0, 0.0), knot(50000.0, 600.0, 0.0),
            knot(150000.0, 600.0, 1.0), knot(250000.0, 1000.0, 1.0)]).unwrap();
    let input = body(&curve, 100000.0, 0.1);
    let result = law.isobaric_thermal_body_candidate(&old, 0.1, 0.0, input, controls(), limits(),
        |_, state| Ok::<_, Infallible>(ConductanceSample { conductance_w_k: 0.1,
            source: state.liquid_mass_fraction() }), || false).unwrap();
    let cg = old.mass_kg()/0.028*3.5*R;
    let q = 0.1*0.1*(600.0-400.0)/(1.0+0.1*0.1/cg);
    near(result.step().boundary.exchanged_heat_j, q, 2e-8);
    assert_eq!(result.thermal().state.temperature_k(), 600.0);
    near(result.thermal().state.liquid_mass_fraction(), 0.5-q/(input.mass_kg*100000.0), 2e-10);
    assert_eq!(result.thermal().report, result.thermal().state.liquid_mass_fraction());
    near(result.step().thermal_entropy_change_j_k, -q/600.0, 2e-10);
    near(result.step().thermal_energy_change_j, -q, 1e-7);
    assert!(result.step().total_energy_residual_j.abs() <= limits().heat.exchange.energy_tolerance_j);
    assert!(result.step().total_entropy_residual_j_k > 0.0);
}

#[test]
fn finite_wall_feedback_keeps_one_actual_exothermic_endpoint() {
    let law = law(6000.0);
    let old = gas(&law, [1.0, 0.0], 400.0);
    let curve = sensible();
    let input = body(&curve, 75000.0, 0.1);
    let result = law.isobaric_thermal_body_candidate(&old, 0.1, 0.0, input, controls(), limits(),
        |gas, state| Ok::<_, Infallible>(ConductanceSample {
            conductance_w_k: 0.03*(gas.state().temperature/400.0)*(state.temperature_k()/350.0),
            source: (gas.state().temperature, state.temperature_k()) }), || false).unwrap();
    let report = result.step();
    let q = report.boundary.exchanged_heat_j;
    assert!(q < 0.0);
    assert!(result.thermal().state.temperature_k() > input.initial.temperature_k());
    let tg = result.gas().state().state().temperature;
    let tb = result.thermal().state.temperature_k();
    near(q, 0.1*0.03*(tg/400.0)*(tb/350.0)*(tb-tg), exchange().exchange.source_tolerance_j);
    let chemical = &report.boundary.chemistry.boundary.reaction;
    near(chemical.extent_mol,
        0.1*result.gas().volume_m3()*chemical.rate.net_progress_rate().unwrap().value(), 1e-11);
    near(tg, 400.0+(6000.0*chemical.extent_mol+q)/(old.mass_kg()/0.028*3.5*R), 2e-7);
    let original = law.isobaric_implicit_candidate(&old, 0.1, q, controls(), || false).unwrap();
    assert_eq!(result.gas(), original.gas());
    assert_eq!(&report.boundary.chemistry, original.step());
    assert!(report.total_energy_residual_j.abs() <= exchange().exchange.energy_tolerance_j);
    assert!(report.total_entropy_residual_j_k > 0.0);
}

#[test]
fn finite_capacity_trajectory_refines_against_the_independent_exact_exchange() {
    let law = law(0.0);
    let initial = gas(&law, [0.5, 0.5], 400.0);
    let curve = sensible();
    let cg = initial.mass_kg()/0.028*3.5*R;
    let cb = 5.0;
    let equilibrium = (cg*400.0+cb*700.0)/(cg+cb);
    let exact = equilibrium-cb/(cg+cb)*300.0*(-0.4*(1.0/cg+1.0/cb)*0.5).exp();
    let mut errors = Vec::new();
    for count in [8, 16, 32] {
        let dt = 0.5/f64::from(count);
        let mut next_gas = initial.clone();
        let mut input = body(&curve, 250000.0, dt);
        let mut work = 0.0;
        for _ in 0..count {
            let candidate = law.isobaric_thermal_body_candidate(&next_gas, dt, 0.0,
                input, controls(), limits(), |_, _| Ok::<_, Infallible>(ConductanceSample {
                    conductance_w_k: 0.4, source: () }), || false).unwrap();
            work += candidate.step().boundary.chemistry.boundary.boundary.work_by_gas_j;
            assert!(candidate.step().total_entropy_residual_j_k >= -limits().entropy_tolerance_j_k);
            let (gas, thermal, _) = candidate.into_parts();
            next_gas = gas;
            input.initial = thermal.state;
        }
        let tg = next_gas.state().state().temperature;
        let tb = input.initial.temperature_k();
        errors.push((tg-exact).abs());
        near(cg*(tg-400.0)+cb*(tb-700.0), 0.0, 5e-6);
        near(work, 100000.0*(next_gas.volume_m3()-initial.volume_m3()), 1e-10);
        assert_eq!(next_gas.mass_kg(), initial.mass_kg());
    }
    for pair in errors.windows(2) { assert!((1.8..2.2).contains(&(pair[0]/pair[1]))); }
}

fn material_wall() -> (MaterialCard, PropertyKey, ResolvedHeatCapacityCurve, QueryPoint) {
    let quantity = QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static));
    let point = |t| QueryPoint::new().with_quantity("T", quantity, t).unwrap().with("lot", 7.0).unwrap();
    let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cp = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let mut claims = ClaimSet::new();
    for (key, lo, hi) in [(conductivity.clone(), 1.0, 3.0),
        (cp.clone(), 500.0, 500.0), (density.clone(), 4000.0, 4000.0)] {
        let dims = key.quantity().dims();
        claims.insert_claim(PropertyClaim { key,
            value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(200.0, lo), (1000.0, hi)], dims },
            validity: ValidityDomain::unconstrained().with_quantity("T", quantity, 200.0, 1000.0)
                .with("lot", 7.0, 7.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
            observations: vec![], provenance: Provenance { source: "manufactured wall conductivity and storage".into(),
                license: "CC0-1.0".into(), artifact: None } }).unwrap();
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "manufactured wall".into(),
        phase: "solid".into(), process: "test".into(), revision: 0 }, claims, vec![]).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
        heat_capacity: &cp, density: &density, temperature_axis: "T",
        lower: &point(200.0), upper: &point(1000.0), selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 1e-5).unwrap();
    (card, conductivity, caloric, point(700.0))
}

#[test]
fn original_material_layers_use_endpoint_conductivity_and_reject_foreign_storage() {
    let law = law(0.0);
    let gas = gas(&law, [0.5, 0.5], 400.0);
    let (card, key, caloric, context) = material_wall();
    let input = body(caloric.curve(), 250000.0, 0.1);
    let layers = || MaterialConductionLayers { material: &card, context: &context,
        temperature_axis: "T", liquid_fraction_axis: None, conductivity: &key,
        selection: MaterialPropertySelection::SingleClaimOnly,
        area_m2: 0.001, gas_distance_m: 0.001, body_distance_m: 0.002 };
    let result = law.isobaric_thermal_layers_candidate(&gas, 0.1, 0.0,
        input, controls(), limits(), layers(), || false).unwrap();
    let exact = layers().resolve(result.gas().state(), result.thermal().state).unwrap();
    assert_eq!(result.thermal().report, exact.source);
    assert_eq!(result.step().conductance_w_k, exact.conductance_w_k);
    let kb = 1.0+2.0*(result.thermal().state.temperature_k()-200.0)/800.0;
    let kg = result.gas().state().state().thermal_conductivity;
    near(result.step().conductance_w_k, 0.001/(0.001/kg+0.002/kb), 1e-14);
    near(result.thermal().report.property("thermal_conductivity").unwrap().value_si(), kb, 1e-13);
    let foreign = sensible();
    assert!(matches!(law.isobaric_thermal_layers_candidate(&gas, 0.1, 0.0,
        body(&foreign, 250000.0, 0.1), controls(), limits(), layers(), || false),
        Err(IsobaricThermalError::Thermal(ThermalPistonError::Interface(_)))));
}

#[test]
fn finite_wall_source_domain_work_and_late_cancellation_fail_atomically() {
    let law = law(0.0);
    let gas = gas(&law, [0.5, 0.5], 400.0);
    let saved = gas.clone();
    let curve = sensible();
    let input = body(&curve, 250000.0, 0.1);
    let saved_body = input.initial;
    let solve = |input, limits, cancel: &mut dyn FnMut() -> bool| law.isobaric_thermal_body_candidate(
        &gas, 0.1, 0.0, input, controls(), limits,
        |_, _| Ok::<_, Infallible>(ConductanceSample { conductance_w_k: 0.4, source: () }), cancel);
    for bad in [UniformEnthalpyStepInput { mass_kg: 0.0, ..input },
        UniformEnthalpyStepInput { duration_s: 0.2, ..input },
        UniformEnthalpyStepInput { internal_heat_j: f64::NAN, ..input },
        UniformEnthalpyStepInput { mass_kg: 1e-12, volume_m3: 1e-12/4000.0, ..input }] {
        assert!(solve(bad, limits(), &mut || false).is_err());
    }
    let limited = IsobaricBodySolve { heat: IsobaricHeatSolve {
        maximum_chemical_evaluations: 1, ..exchange() }, ..limits() };
    assert!(matches!(solve(input, limited, &mut || false),
        Err(IsobaricThermalError::ChemicalEvaluationLimit { .. })));
    assert!(matches!(law.isobaric_thermal_body_candidate(&gas, 0.1, 0.0,
        input, controls(), limits(), |_, _| Err::<ConductanceSample<()>, _>("missing wall source"), || false),
        Err(IsobaricThermalError::Thermal(ThermalPistonError::Interface("missing wall source")))));
    let mut count = 0;
    solve(input, limits(), &mut || { count += 1; false }).unwrap();
    for stop in [1, count/2, count] {
        let mut polls = 0;
        assert!(matches!(solve(input, limits(), &mut || { polls += 1; polls == stop }),
            Err(IsobaricThermalError::Cancelled)));
    }
    assert_eq!(gas, saved);
    assert_eq!(input.initial, saved_body);
}
