//! Reacting thermal storage on original source and finite-heat owners.
use super::*;
use super::conduction::GasEquilibriumConductionLaw;
use crate::thermal::conduction::ConductivityCurve;
use crate::thermochemical_gas::energy::equilibrium::tests::{fixture, isomers};
use crate::thermochemical_gas::tests::{
    point, species, transport, transport_in_domain,
};
use fs_evidence::ValidityDomain;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget,
    DepositionSimulation, Event, MaterialId, StepConfig, TimedEvent,
};
use fs_thermal::enthalpy::conduction::implicit::ImplicitConductionLaw;
use fs_thermochem::equilibrium::composition::{
    EquilibriumEnthalpySolveV1, EquilibriumSolveV1,
};
use fs_thermochem::{
    Nasa9RegionV1, Nasa9StandardStateModelV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
};

fn controls() -> GasEquilibriumEnthalpySolve {
    let mut equilibrium = EquilibriumSolveV1::default();
    equilibrium.relative_balance_tolerance = 1e-13;
    equilibrium.dimensionless_potential_tolerance = 1e-12;
    GasEquilibriumEnthalpySolve {
        equilibrium: EquilibriumEnthalpySolveV1 {
            equilibrium, minimum_temperature_k: 200.0, maximum_temperature_k: 1500.0,
            enthalpy_tolerance_j: 1e-7, maximum_temperature_iterations: 128,
        },
        mass_tolerance_kg: 1e-12, enthalpy_tolerance_j: 1e-6,
    }
}
fn feed() -> [(SpeciesId, Amount); 1] {
    [(SpeciesId::new("A").unwrap(), Amount::new(2.0))]
}
fn close(a: f64, b: f64, allowance: f64) {
    assert!((a-b).abs() <= allowance, "{a:.17e} != {b:.17e}");
}
fn fraction(t: f64) -> f64 { 1.0 / (1.0 + (-24000.0 / (R*t)).exp()) }
fn h(t: f64) -> f64 {
    (4.0*R*(t-200.0) - 24000.0*(fraction(t)-fraction(200.0))) / 0.01
}
fn cp(t: f64) -> f64 {
    let y = fraction(t);
    (4.0*R + 24000.0_f64.powi(2)*y*(1.0-y)/(R*t*t)) / 0.01
}
fn entropy(t: f64) -> f64 {
    let y = fraction(t);
    R/0.01*(4.0*t.ln() - y*y.ln() - (1.0-y)*(1.0-y).ln())
}

#[test]
fn g1_reacting_thermal_law_retains_sources_mixing_entropy_and_capacity_bound() {
    let (bank, chemistry) = isomers();
    let original = bank.clone();
    let context = point(900.0, 1e5).with("sample_id", 7.0).unwrap();
    let law = GasEquilibriumThermalStorage::new(
        &bank, &chemistry, &context, &feed(), controls(), || false,
    ).unwrap();
    assert_eq!(law.storage_basis(), ThermalStorageBasis::Enthalpy);
    assert_eq!(law.feed_mass_kg(), 0.02);
    assert_eq!(law.reference_state().parcel().state().state().temperature, 200.0);
    assert_eq!(law.reference_state().parcel().point().axes()["sample_id"], 7.0);
    let frozen_bound = 0.01/(4.0*R);
    assert!(law.maximum_temperature_slope() >= frozen_bound);
    assert!(law.maximum_temperature_slope() <= frozen_bound*(1.0+1e-10));
    for t in [200.0, 400.0, 900.0, 1300.0, 1500.0] {
        let (actual_h, actual_cp) =
            law.enthalpy_and_heat_capacity_at_temperature(Temperature::new(t), || false).unwrap();
        close(actual_h, h(t), 2e-5);
        close(actual_cp, cp(t), 2e-7);
        close(law.temperature_k(actual_h).unwrap(), t, 2e-8);
        close(law.bulk_density_kg_m3(actual_h).unwrap().unwrap(), 1e5*0.01/(R*t), 2e-12);
        assert!(1.0/actual_cp <= law.maximum_temperature_slope());
    }
    let state = law.resolve_specific_enthalpy(h(900.0), || false).unwrap();
    assert_eq!(state.parcel().mass_kg(), law.feed_mass_kg());
    assert_eq!(state.parcel().model().thermodynamics().components(), chemistry.models());
    assert_eq!(state.parcel().model().transport_materials(), bank.transport_materials());
    assert!(state.native_sample().response.heat_capacity_p_j_k
        > state.parcel().mass_kg()*state.parcel().state().state().specific_heat_cp);
    close(law.thermal_entropy_change_j_kg_k(h(400.0), h(1300.0), || false).unwrap(),
        entropy(1300.0)-entropy(400.0), 2e-7);
    close(law.thermal_entropy_change_j_kg_k(h(1300.0), h(400.0), || false).unwrap(),
        entropy(400.0)-entropy(1300.0), 2e-7);
    assert!(law.advected_composition().is_none());
    assert_eq!(bank, original);
}

#[test]
fn g1_finite_isobaric_heat_uses_existing_timeline_with_independent_cell_mass() {
    let (bank, chemistry) = isomers();
    let law = GasEquilibriumThermalStorage::new(
        &bank, &chemistry, &point(900.0, 1e5), &feed(), controls(), || false,
    ).unwrap();
    let cell_mass = 0.03; // Deliberately differs from the chemistry normalization.
    let heat = cell_mass*(h(1200.0)-h(400.0));
    let boundary = |source_w| Boundary {
        ambient_k: 400.0, htc_w_m2_k: 0.0, reservoir: None, source_w,
    };
    let geometry = CellGeometry {
        start_m: [0.0,0.0,0.0], end_m: [1.0,0.0,0.0],
        volume_m3: cell_mass/(1e5*0.01/(R*400.0)), surface_area_m2: 1.0,
    };
    let original = DepositionSimulation::new(
        vec![
            TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
                id: CellId(1), segment_id: 1, body_id: 1, material_id: MaterialId(1),
                mass_kg: cell_mass, specific_h_j_kg: h(400.0), geometry,
                boundary: boundary(heat),
            }) },
            TimedEvent { time_s: 1.0, event: Event::Boundary(CellId(1), boundary(-heat)) },
        ],
        &[(MaterialId(1), &law as &dyn ThermalLaw)],
        DepositionBudget { max_cells: 1, max_links: 0, max_events: 2, max_substeps: 128 },
    ).unwrap();
    let config = StepConfig { max_dt_s: 0.125, energy_tolerance_j: 1e-6 };
    let (initial, _) = original.advance_to(0.0, config, || Ok(())).unwrap();
    let (hot, audit) = initial.advance_to(1.0, config, || Ok(())).unwrap();
    close(hot.temperature_k(CellId(1)).unwrap(), 1200.0, 2e-7);
    close(audit.external_heat_j, heat, 1e-7);
    close(audit.stored_energy_change_j, heat, 1e-7);
    assert!(audit.residual_j.abs() <= config.energy_tolerance_j);
    assert_eq!(hot.mass_kg(), cell_mass);
    assert_eq!(hot.cells()[&CellId(1)].geometry, geometry);
    // Reference geometry is not a rigid-volume constraint; equilibrium
    // density supplies the current volume under the fixed-pressure law.
    close(cell_mass/law.bulk_density_kg_m3(h(1200.0)).unwrap().unwrap(),
        3.0*geometry.volume_m3, 1e-11);
    let (cooled, audit) = hot.advance_to(2.0, config, || Ok(())).unwrap();
    close(cooled.temperature_k(CellId(1)).unwrap(), 400.0, 2e-7);
    close(audit.external_heat_j, -heat, 1e-7);
    assert!(audit.residual_j.abs() <= config.energy_tolerance_j);
    assert_eq!(initial.time_s(), 0.0);
    close(initial.temperature_k(CellId(1)).unwrap(), 400.0, 2e-8);
    assert!(initial.advance_to(1.0, config, || Err(DepositionError::Cancelled)).is_err());
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn g1_conduction_derivatives_include_reacting_cp_and_molecular_count() {
    let mut models = Vec::new();
    for (name, mass, cp_r) in [("A",0.001,2.5), ("A2",0.002,3.5)] {
        let source = species(name, mass, cp_r, 0.0);
        let mut card = source.regions()[0].card().clone();
        // Both g/RT vanish at 1000 K, fixing independent Kp=1 there.
        card.parameters.get_mut("a8").unwrap().value = cp_r*(1.0-1000.0_f64.ln());
        models.push(Nasa9StandardStateModelV1::new(
            source.species().clone(), source.molar_mass(), source.convention().clone(),
            vec![Nasa9RegionV1::from_card(card).unwrap()],
        ).unwrap());
    }
    let (bank, chemistry) = fixture(models,
        vec![transport("A",1.7e-5,0.025), transport("A2",2.0e-5,0.15)],
        &["H"], vec![vec![1,2]]);
    let pressure = 30000.0;
    let law = GasEquilibriumThermalStorage::new(
        &bank, &chemistry, &point(900.0, pressure),
        &[(SpeciesId::new("A2").unwrap(), Amount::new(1.0))], controls(), || false,
    ).unwrap();
    let conductor = GasEquilibriumConductionLaw::try_new(&law,
        ConductivityCurve::new(vec![[200.0,0.1],[800.0,0.22],[1500.0,0.5]]).unwrap(),
    ).unwrap();
    let (stored_h, _) = law.enthalpy_and_heat_capacity_at_temperature(
        Temperature::new(1000.0), || false).unwrap();
    let state = conductor.state_at_specific_enthalpy(stored_h).unwrap();
    let alpha = (1.0_f64/(1.0+4.0*pressure/1e5)).sqrt();
    let d_alpha = 0.75*alpha*(1.0-alpha*alpha)/1000.0;
    let capacity = R*(3.5+1.5*alpha+1500.0*d_alpha);
    close(state.temperature_k, 1000.0, 2e-8);
    close(state.temperature_slope_kg_k_j, 0.002/capacity, 2e-13);
    close(state.conductivity_w_m_k, 0.3, 1e-11);
    close(state.conductivity_slope_kg_w_m_k_j, 0.0004*0.002/capacity, 2e-14);
    let volume_derivative = conductor
        .specific_volume_derivative_at_specific_enthalpy(stored_h).unwrap().unwrap();
    let expected = R*(1.0+alpha+1000.0*d_alpha)/(pressure*capacity);
    close(volume_derivative, expected, 2e-12);
    assert!(volume_derivative > R*(1.0+alpha)/(pressure*capacity));
    assert_eq!(conductor.enthalpy_bounds_j_kg(), law.specific_enthalpy_domain_j_kg());
    close(conductor.caloric_law().temperature_k(stored_h).unwrap(), 1000.0, 2e-8);
    assert!(GasEquilibriumConductionLaw::try_new(&law,
        ConductivityCurve::new(vec![[300.0,0.1],[1500.0,0.5]]).unwrap()).is_err());
}

#[test]
fn g4_storage_refuses_unsupported_material_endpoints_and_discontinuous_sources() {
    let (_, chemistry) = isomers();
    let (partial, _) = fixture(chemistry.models().to_vec(), vec![
        transport_in_domain("A",1.7e-5,0.025,[200.0,950.0],[1000.0,1e7]),
        transport_in_domain("B",2.0e-5,0.15,[200.0,950.0],[1000.0,1e7]),
    ], &["H"], vec![vec![1,1]]);
    assert!(GasEquilibriumThermalStorage::new(
        &partial, &chemistry, &point(900.0,1e5), &feed(), controls(), || false,
    ).is_err());
    let source = species("A",0.01,4.0,0.0);
    let mut lower = source.regions()[0].card().clone();
    lower.validity = ValidityDomain::unconstrained().with("T",200.0,1000.0);
    let mut upper = lower.clone();
    upper.validity = ValidityDomain::unconstrained().with("T",1000.0,2000.0);
    upper.parameters.get_mut("a7").unwrap().value -= 4000.0;
    let discontinuous = Nasa9StandardStateModelV1::new(
        source.species().clone(), source.molar_mass(), source.convention().clone(),
        vec![Nasa9RegionV1::from_card(lower).unwrap(), Nasa9RegionV1::from_card(upper).unwrap()],
    ).unwrap();
    let (bank, chemistry) = fixture(vec![discontinuous],
        vec![transport("A",1.7e-5,0.025)], &["H"], vec![vec![1]]);
    assert!(GasEquilibriumThermalStorage::new(
        &bank, &chemistry, &point(900.0,1e5), &feed(), controls(), || false,
    ).is_err());
}

#[test]
fn g4_refused_and_cancelled_queries_preserve_the_immutable_material_law() {
    let (bank, chemistry) = isomers();
    let original = bank.clone();
    let make = || GasEquilibriumThermalStorage::new(
        &bank, &chemistry, &point(900.0,1e5), &feed(), controls(), || false,
    ).unwrap();
    let law = make();
    let reference = law.reference_state().clone();
    let expected = law.caloric_state_at_specific_enthalpy(h(900.0), || false).unwrap();
    for value in [-1.0, f64::NAN, f64::INFINITY, law.specific_enthalpy_domain_j_kg()[1].next_up()] {
        assert!(law.temperature_k(value).is_err());
        assert_eq!(law.caloric_state_at_specific_enthalpy(h(900.0), || false).unwrap(), expected);
    }
    let mut count = 0;
    law.resolve_specific_enthalpy(h(900.0), || { count += 1; false }).unwrap();
    let cancelled = make();
    let mut polls = 0;
    assert!(matches!(cancelled.resolve_specific_enthalpy(h(900.0), || {
        polls += 1; polls == count
    }), Err(Error::Cancelled)));
    assert_eq!(polls, count);
    assert_eq!(cancelled.reference_state(), &reference);
    assert_eq!(cancelled.caloric_state_at_specific_enthalpy(h(900.0), || false).unwrap(), expected);
    assert_eq!(bank, original);
}
