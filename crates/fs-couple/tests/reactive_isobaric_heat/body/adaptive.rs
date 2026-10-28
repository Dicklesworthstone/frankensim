use super::*;
use fs_couple::cavity_phs::piston::reactive::implicit::coupled::isobaric::advance::{
    IsobaricAdvanceError, IsobaricAdvanceSolve,
    body::{IsobaricBodyAdvanceError, IsobaricBodyAdvanceSolve},
};

fn accuracy() -> IsobaricBodyAdvanceSolve {
    IsobaricBodyAdvanceSolve {
        trajectory: IsobaricAdvanceSolve { absolute_amount_mol: 2e-5,
            absolute_temperature_k: 0.3, absolute_volume_m3: 1e-6, relative_tolerance: 0.0,
            maximum_substeps: 256, maximum_endpoint_evaluations: 511,
            maximum_chemical_evaluations: 1_048_576 },
        absolute_body_enthalpy_j_kg: 100.0,
    }
}

#[test]
fn adaptive_hot_and_cold_finite_bodies_match_exact_exchange_and_combined_entropy() {
    let law = law(0.0);
    let curve = sensible();
    let mut solve = limits();
    // Half a second at the initial 120 W contact needs wider integrated
    // heat bounds than the original 0.1 s endpoint fixture.
    solve.heat.exchange.heat_bracket_j = [-80.0, 80.0];
    for (tg0, tb0) in [(400.0, 700.0), (700.0, 400.0)] {
        let initial = gas(&law, [0.5, 0.5], tg0);
        let input = body(&curve, 500.0 * (tb0 - 200.0), 0.5);
        let result = law.isobaric_thermal_body_adaptive_candidate(&initial, 0.5, 0.0,
            input, controls(), solve, accuracy(), |gas, state| Ok::<_, Infallible>(
                ConductanceSample { conductance_w_k: 0.4,
                    source: (gas.state().temperature, state.temperature_k()) }), || false).unwrap();
        let (cg, cb) = (initial.mass_kg() / 0.028 * 3.5 * R, 5.0);
        let equilibrium = (cg * tg0 + cb * tb0) / (cg + cb);
        let difference = (tg0 - tb0) * (-0.4 * (1.0 / cg + 1.0 / cb) * 0.5).exp();
        let tg = equilibrium + cb / (cg + cb) * difference;
        let tb = equilibrium - cg / (cg + cb) * difference;
        near(result.gas().state().state().temperature, tg, 1.1 * accuracy().trajectory.absolute_temperature_k);
        near(result.thermal().state.temperature_k(), tb, 1.1 * accuracy().trajectory.absolute_temperature_k);
        let report = result.report();
        assert!(report.trajectory.steps.len() > 2);
        assert_eq!(result.gas().mass_kg(), initial.mass_kg());
        assert_eq!(result.gas().state().state().pressure, 100000.0);
        assert_eq!(result.thermal().external_heat_j, -report.exchanged_heat_j);
        near(report.thermal_energy_change_j, -report.exchanged_heat_j, 1e-7);
        near(report.thermal_energy_change_j, input.mass_kg
            * (result.thermal().state.specific_enthalpy_j_kg() - input.initial.specific_enthalpy_j_kg()), 1e-7);
        near(report.trajectory.enthalpy_change_j, report.exchanged_heat_j, 1e-7);
        near(report.trajectory.stored_energy_change_j + report.trajectory.work_by_gas_j
            + report.thermal_energy_change_j, 0.0, 1e-7);
        let actual_tg = result.gas().state().state().temperature;
        let actual_tb = result.thermal().state.temperature_k();
        let exact_entropy = cg * (actual_tg / tg0).ln() + cb * (actual_tb / tb0).ln();
        near(report.total_entropy_residual_j_k, exact_entropy, 1e-8);
        assert!(report.total_entropy_residual_j_k > 0.0);
        assert!(report.applied_contact_entropy_j_k > 0.0);
        for interval in &report.trajectory.steps {
            assert_eq!(interval.thermal.report,
                (interval.gas_temperature_k, interval.thermal.state.temperature_k()));
            assert_eq!(interval.thermal.external_heat_j, -interval.step.boundary.exchanged_heat_j);
        }
    }
}

fn nonlinear_law() -> GasReaction {
    let original = law(6000.0);
    let models = original.catalog().thermodynamics().components().iter().map(|model| {
        let mut card = model.regions()[0].card().clone();
        // Reuse the manufactured slope from fs-material's original nonlinear
        // isobaric fixture. These are synthetic isomers, not fitted substances.
        card.parameters.get_mut("a3").unwrap().value = 0.001;
        Nasa9StandardStateModelV1::new(model.species().clone(), model.molar_mass(),
            model.convention().clone(), vec![Nasa9RegionV1::from_card(card).unwrap()]).unwrap()
    }).collect::<Vec<_>>();
    let old = original.kinetics().reverse_rate_closure().equilibrium();
    let equilibrium = IdealGasReactionEquilibriumV1::new(old.stoichiometric_matrix().clone(),
        old.conservation_certificate(), old.reaction().clone(), models.clone()).unwrap();
    let catalog = ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(models,
        original.catalog().thermodynamics().mole_composition().clone()).unwrap(),
        original.catalog().transport_materials().to_vec(), None).unwrap();
    GasReaction::new(catalog, DeclaredStoichiometricMassActionV1::new(
        ThermodynamicReverseRateClosureV1::new(equilibrium,
            original.kinetics().reverse_rate_closure().forward_progress_rate_scale()),
        MassActionKineticOrderConventionV1::CallerDeclaredStoichiometricCoefficients)).unwrap()
        .with_implicit_mechanism(original.implicit_mechanism().unwrap().clone()).unwrap()
}

fn reference(n: f64, duration: f64) -> [f64; 3] {
    // Independent continuous extent/Tgas/Tbody equations. The same NASA-9
    // slope appears in Cp(T); formation release is counted once through dxi.
    let rhs = |state: [f64; 3]| {
        let [xi, tg, tb] = state;
        let extent = 20.0 * R * tg / 100000.0 * (n - xi - (-6000.0 / (R * tg)).exp() * xi);
        let q = 0.4 * (tb - tg);
        [extent, (6000.0 * extent + q + 2.0 / duration) / (n * R * (3.5 + 0.001 * tg)),
            (3.0 / duration - q) / 5.0]
    };
    let add = |a: [f64; 3], b: [f64; 3], scale: f64| std::array::from_fn(|i| a[i] + scale * b[i]);
    let dt = duration / 8192.0;
    let mut state = [0.0, 400.0, 350.0];
    for _ in 0..8192 {
        let a = rhs(state);
        let b = rhs(add(state, a, 0.5 * dt));
        let c = rhs(add(state, b, 0.5 * dt));
        let d = rhs(add(state, c, dt));
        state = std::array::from_fn(|i| state[i] + dt * (a[i] + 2.0*b[i] + 2.0*c[i] + d[i]) / 6.0);
    }
    state
}

#[test]
fn nonlinear_cp_reaction_and_finite_body_match_independent_coupled_ode() {
    let law = nonlinear_law();
    let initial = gas(&law, [1.0, 0.0], 400.0);
    let curve = sensible();
    let duration = 0.25;
    let input = UniformEnthalpyStepInput { internal_heat_j: 3.0,
        ..body(&curve, 75000.0, duration) };
    let result = law.isobaric_thermal_body_adaptive_candidate(&initial, duration, 2.0,
        input, controls(), limits(), accuracy(), |gas, state| Ok::<_, Infallible>(
            ConductanceSample { conductance_w_k: 0.4,
                source: (gas.state().temperature, state.temperature_k()) }), || false).unwrap();
    let report = result.report();
    let n = initial.mass_kg() / 0.028;
    let [xi, tg, tb] = reference(n, duration);
    near(report.trajectory.reaction_extents_mol[0].1, xi, 1.1 * accuracy().trajectory.absolute_amount_mol);
    near(result.gas().state().state().temperature, tg, 1.1 * accuracy().trajectory.absolute_temperature_k);
    near(result.thermal().state.temperature_k(), tb, 1.1 * accuracy().trajectory.absolute_temperature_k);
    let actual_t = result.gas().state().state().temperature;
    let actual_xi = report.trajectory.reaction_extents_mol[0].1;
    let gas_h = n * R * (actual_t - 400.0) * (3.5 + 0.0005 * (actual_t + 400.0)) - 6000.0 * actual_xi;
    near(gas_h, 2.0 + report.exchanged_heat_j, 1e-7);
    near(gas_h + input.mass_kg * (result.thermal().state.specific_enthalpy_j_kg()
        - input.initial.specific_enthalpy_j_kg()), 5.0, 1e-7);
    assert!(report.exchanged_heat_j < 0.0);
    assert!(report.total_enthalpy_residual_j.abs() <= limits().heat.exchange.energy_tolerance_j);
    assert!(report.total_energy_residual_j.abs() <= limits().heat.exchange.energy_tolerance_j);
    let mut replay_gas = initial.clone();
    let mut replay_body = input.initial;
    let (mut gas_heat, mut body_heat, mut external_entropy) = (0.0, 0.0, 0.0);
    let count = report.trajectory.steps.len() as f64;
    let mut inner = limits();
    inner.heat.exchange.source_tolerance_j /= count;
    inner.heat.exchange.heat_bracket_j = inner.heat.exchange.heat_bracket_j.map(|q| q/count);
    for interval in &report.trajectory.steps {
        let dt = interval.end_s - interval.start_s;
        let candidate = law.isobaric_thermal_body_candidate(&replay_gas, dt,
            interval.prescribed_gas_heat_j, UniformEnthalpyStepInput {
                initial: replay_body, duration_s: dt, internal_heat_j: interval.prescribed_body_heat_j, ..input },
            controls(), inner, |gas, state| Ok::<_, Infallible>(ConductanceSample {
                conductance_w_k: 0.4, source: (gas.state().temperature, state.temperature_k()) }), || false).unwrap();
        assert_eq!(candidate.step(), &interval.step);
        assert_eq!(candidate.thermal(), &interval.thermal);
        gas_heat += interval.prescribed_gas_heat_j;
        body_heat += interval.prescribed_body_heat_j;
        external_entropy += interval.prescribed_gas_heat_j / interval.gas_temperature_k
            + interval.prescribed_body_heat_j / interval.thermal.state.temperature_k();
        let (gas, thermal, _) = candidate.into_parts();
        replay_gas = gas; replay_body = thermal.state;
    }
    assert_eq!(&replay_gas, result.gas());
    assert_eq!(replay_body, result.thermal().state);
    near(gas_heat, 2.0, 1e-14); near(body_heat, 3.0, 1e-14);
    near(report.prescribed_heat_entropy_j_k, external_entropy, 1e-14);
    assert!(report.total_entropy_residual_j_k >= -limits().entropy_tolerance_j_k);
}

#[test]
fn latent_entry_and_flat_temperature_are_controlled_by_body_enthalpy() {
    let law = law(0.0);
    let initial = gas(&law, [0.5, 0.5], 400.0);
    let knot = |h, t, f| EnthalpyPhaseKnot { specific_enthalpy_j_kg: h,
        temperature_k: t, liquid_mass_fraction: f, bulk_density_kg_m3: 4000.0 };
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x42; 32]),
        vec![knot(0.0, 300.0, 0.0), knot(50000.0, 600.0, 0.0),
            knot(150000.0, 600.0, 1.0), knot(250000.0, 1000.0, 1.0)]).unwrap();
    let accuracy = IsobaricBodyAdvanceSolve { absolute_body_enthalpy_j_kg: 50.0,
        trajectory: IsobaricAdvanceSolve { absolute_temperature_k: 1e6,
            absolute_volume_m3: 1.0, absolute_amount_mol: 1.0, ..accuracy().trajectory } };
    for h0 in [100000.0, 150100.0] {
        let input = body(&curve, h0, 0.5);
        let result = law.isobaric_thermal_body_adaptive_candidate(&initial, 0.5, 0.0,
            input, controls(), limits(), accuracy, |_, state| Ok::<_, Infallible>(ConductanceSample {
                conductance_w_k: 0.4, source: state.liquid_mass_fraction() }), || false).unwrap();
        let cg = initial.mass_kg() / 0.028 * 3.5 * R;
        let (crossing_time, crossing_tg, crossing_h) = if h0 <= 150000.0 { (0.0, 400.0, h0) } else {
            let cb = input.mass_kg * 250.0;
            let tb0 = input.initial.temperature_k();
            let equilibrium = (cg*400.0 + cb*tb0)/(cg+cb);
            let rate = 0.4*(1.0/cg + 1.0/cb);
            let time = -((600.0-equilibrium)/(tb0-equilibrium)).ln()/rate;
            (time, 400.0 + input.mass_kg*(h0-150000.0)/cg, 150000.0)
        };
        let exact_tg = 600.0 - (600.0-crossing_tg)*(-0.4*(0.5-crossing_time)/cg).exp();
        let exact_h = crossing_h - cg*(exact_tg-crossing_tg)/input.mass_kg;
        near(result.thermal().state.specific_enthalpy_j_kg(), exact_h,
            1.2 * accuracy.absolute_body_enthalpy_j_kg);
        assert_eq!(result.thermal().state.temperature_k(), 600.0);
        assert!(result.report().trajectory.steps.len() > 2);
        assert!(result.thermal().state.liquid_mass_fraction() < input.initial.liquid_mass_fraction());
        near(result.report().thermal_energy_change_j, -result.report().exchanged_heat_j, 1e-7);
        near(result.report().thermal_energy_change_j,
            input.mass_kg * (result.thermal().state.specific_enthalpy_j_kg() - h0), 1e-7);
        if h0 < 150000.0 {
            assert!(result.report().trajectory.steps.iter().all(|step| step.thermal.state.temperature_k() == 600.0));
            near(result.report().thermal_entropy_change_j_k, -result.report().exchanged_heat_j / 600.0, 1e-9);
        }
    }
}

#[test]
fn material_bound_conductivity_refreshes_on_every_accepted_phase_endpoint() {
    let law = law(6000.0);
    let gas = gas(&law, [1.0, 0.0], 400.0);
    let (card, key, caloric, context) = material_wall();
    let input = body(caloric.curve(), 250000.0, 0.2);
    let layers = || MaterialConductionLayers { material: &card, context: &context,
        temperature_axis: "T", liquid_fraction_axis: None, conductivity: &key,
        selection: MaterialPropertySelection::SingleClaimOnly,
        area_m2: 0.001, gas_distance_m: 0.001, body_distance_m: 0.002 };
    let result = law.isobaric_thermal_layers_adaptive_candidate(&gas, 0.2, 0.0,
        input, controls(), limits(), accuracy(), layers(), || false).unwrap();
    for interval in &result.report().trajectory.steps {
        let source = &interval.thermal.report;
        assert_eq!(source.card_identity(), card.content_hash());
        assert_eq!(source.query_point().iter().find(|(axis, _)| axis == "lot")
            .map(|(_, value)| *value), Some(7.0));
        near(source.property("thermal_conductivity").unwrap().value_si(),
            1.0 + 2.0*(interval.thermal.state.temperature_k()-200.0)/800.0, 1e-12);
    }
    let exact = layers().resolve(result.gas().state(), result.thermal().state).unwrap();
    assert_eq!(result.report().trajectory.steps.last().unwrap().thermal.report, exact.source);
    let foreign = sensible();
    assert!(matches!(law.isobaric_thermal_layers_adaptive_candidate(&gas, 0.2, 0.0,
        body(&foreign, 250000.0, 0.2), controls(), limits(), accuracy(), layers(), || false),
        Err(IsobaricBodyAdvanceError::Thermal(IsobaricThermalError::Thermal(ThermalPistonError::Interface(_))))));
}

#[test]
fn failed_heat_root_spends_chemical_work_and_late_cancellation_keeps_both_owners() {
    let law = law(0.0);
    let initial = gas(&law, [0.5, 0.5], 400.0);
    let saved = initial.clone();
    let curve = sensible();
    let input = body(&curve, 250000.0, 0.1);
    let mut bounded = limits();
    bounded.heat.exchange.maximum_iterations = 3;
    bounded.heat.exchange.source_tolerance_j = 1e-6;
    bounded.heat.exchange.energy_tolerance_j = 1e-6;
    let sample = |gas: &fs_material::thermochemical_gas::ResolvedThermochemicalGasState,
        state: fs_material::phase::EquilibriumPhaseState| Ok::<_, Infallible>(ConductanceSample {
            conductance_w_k: 0.03*(gas.state().temperature/400.0)*(state.temperature_k()/350.0), source: () });
    let mut callbacks = 0;
    assert!(matches!(law.isobaric_thermal_body_candidate(&initial, 0.1, 0.0,
        input, controls(), bounded, |gas, state| { callbacks += 1; sample(gas,state) }, || false),
        Err(IsobaricThermalError::Thermal(ThermalPistonError::NotConverged { .. }))));
    // Equilibrated identical isomers require exactly one original chemical
    // evaluation per thermal proposal, so source calls independently count work.
    assert!(callbacks >= 3);
    let capped = IsobaricBodyAdvanceSolve { trajectory: IsobaricAdvanceSolve {
        maximum_chemical_evaluations: callbacks, ..accuracy().trajectory }, ..accuracy() };
    assert_eq!(law.isobaric_thermal_body_adaptive_candidate(&initial, 0.1, 0.0,
        input, controls(), bounded, capped, sample, || false),
        Err(IsobaricBodyAdvanceError::Trajectory(IsobaricAdvanceError::ChemicalEvaluationLimit {
            evaluations: callbacks, maximum: callbacks, endpoint_evaluations: 1 })));
    let result = law.isobaric_thermal_body_adaptive_candidate(&initial, 0.1, 0.0,
        input, controls(), bounded, accuracy(), sample, || false).unwrap();
    assert!(result.report().trajectory.refused_meshes > 0);
    assert!(result.report().trajectory.chemical_evaluations >
        result.report().trajectory.steps.iter().map(|step| step.step.boundary.chemical_evaluations).sum());
    let mut count = 0;
    law.isobaric_thermal_body_adaptive_candidate(&initial, 0.1, 0.0, input,
        controls(), limits(), accuracy(), sample, || { count += 1; false }).unwrap();
    let mut polls = 0;
    assert_eq!(law.isobaric_thermal_body_adaptive_candidate(&initial, 0.1, 0.0, input,
        controls(), limits(), accuracy(), sample, || { polls += 1; polls == count }),
        Err(IsobaricBodyAdvanceError::Trajectory(IsobaricAdvanceError::Cancelled)));
    assert_eq!(initial, saved);
    assert_eq!(input.initial, curve.state_at_specific_enthalpy(250000.0).unwrap());
}
