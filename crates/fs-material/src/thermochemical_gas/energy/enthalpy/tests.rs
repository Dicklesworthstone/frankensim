//! Synthetic NASA-9 coefficients provide an independent caloric oracle.
use super::*;
use crate::thermochemical_gas::tests::{model, point};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation,
    Event, MaterialId, StepConfig, TimedEvent,
};

#[path = "tests/selective.rs"]
mod selective;

fn controls() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 280.0, maximum_temperature_k: 360.0,
        energy_tolerance_j_kg: 1e-8, maximum_iterations: 64 }
}
fn mole(fractions: Vec<f64>) -> Composition {
    Composition::new(CompositionBasis::MoleFraction, fractions).unwrap()
}

fn thermal_config() -> StepConfig {
    StepConfig { max_dt_s: 1.0, energy_tolerance_j: 1e-8 }
}

fn thermal_timeline<'a>(
    materials: &[(MaterialId, &'a dyn ThermalLaw)], h: f64, source_w: f64,
) -> DepositionSimulation<'a> {
    let cell = DepositedCell {
        id: CellId(1), segment_id: 1, body_id: 1, material_id: MaterialId(1),
        mass_kg: 1.0, specific_h_j_kg: h,
        geometry: CellGeometry { start_m: [0.0; 3], end_m: [1.0, 0.0, 0.0],
            volume_m3: 1.0, surface_area_m2: 6.0 },
        boundary: Boundary { ambient_k: 320.0, htc_w_m2_k: 0.0,
            reservoir: None, source_w },
    };
    DepositionSimulation::new(vec![TimedEvent { time_s: 0.0, event: Event::Birth(cell) }],
        materials, DepositionBudget { max_cells: materials.len().max(2), max_links: 0,
            max_events: 1, max_substeps: 8 })
        .unwrap().advance_to(0.0, thermal_config(), || Ok(())).unwrap().0
}

#[test]
fn pressure_enthalpy_inverse_uses_nasa_curvature_and_returns_one_coherent_medium() {
    let source = model(vec![0.75, 0.25], 0.001);
    let law = GasEnthalpyLaw::from_mole(&source, mole(vec![0.75, 0.25]),
        point(300.0, 140000.0), controls()).unwrap();
    let t: f64 = 337.0;
    let molar_mass = 0.75 * 0.028 + 0.25 * 0.004;
    // NASA-9 fixture h/R = x_A*(3.5*T+0.001*T^2/2)+x_B*2.5*T.
    let expected_h = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
        * (3.25 * t + 0.75 * 0.001 * t.powi(2) / 2.0) / molar_mass;
    let solved = law.resolve_enthalpy(expected_h, || false).unwrap();
    assert!((solved.state().state().temperature - t).abs() < 1e-9);
    assert!(solved.residual_j_kg().abs() <= controls().energy_tolerance_j_kg);
    assert_eq!(solved.state().state().pressure, 140000.0);
    solved.state().state().validate_coherence().unwrap();
    let expected_density = 140000.0 * molar_mass / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * t);
    assert!((solved.state().state().density / expected_density - 1.0).abs() < 1e-11);
    for t in [280.0, 300.0, 337.0, 360.0] {
        let h = law.enthalpy_at_temperature(t).unwrap();
        assert!((law.temperature_k(h).unwrap() - t).abs() < 1e-9);
        let exact_cp = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
            * (3.25 + 0.75 * 0.001 * t) / molar_mass;
        assert!(law.maximum_temperature_slope() >= 1.0 / exact_cp);
    }
    // A frozen initial cp would miss the quadratic enthalpy by a finite amount.
    let h300 = law.enthalpy_at_temperature(300.0).unwrap();
    let cp300 = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * (3.25 + 0.00075 * 300.0) / molar_mass;
    assert!((h300 + cp300 * (t - 300.0) - expected_h).abs() > 100.0);
}

#[test]
fn relative_caloric_sample_integrates_mixture_cp_and_preserves_its_bracket() {
    let source = model(vec![0.75, 0.25], 0.001);
    let law = GasEnthalpyLaw::from_mole(&source, mole(vec![0.75, 0.25]),
        point(300.0, 140000.0), controls()).unwrap();
    let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / (0.75 * 0.028 + 0.25 * 0.004);
    for (reference, temperature) in [(280.0, 360.0), (360.0, 280.0), (337.0, 300.0),
        (320.0, 320.0), (320.0, 320.0_f64.next_up())]
    {
        let (change, capacity) = law.relative_enthalpy_and_heat_capacity(reference, temperature).unwrap();
        // The independent integral of Cp/R = 3.25 + 0.00075*T retains the
        // small temperature increment before either absolute h is evaluated.
        let expected_change = r * (temperature - reference)
            * (3.25 + 0.000375 * (reference + temperature));
        let expected_capacity = r * (3.25 + 0.00075 * temperature);
        assert!((change - expected_change).abs() <= 2e-14 * expected_change.abs());
        assert!((capacity - expected_capacity).abs() <= 2e-14 * expected_capacity);
        assert_eq!(capacity, law.model().thermodynamics().evaluate_properties(
            Temperature::new(temperature), fs_qty::Pressure::new(140000.0),
        ).unwrap().mass_specific_cp().value());
        assert_eq!(law.relative_enthalpy_and_heat_capacity(temperature, reference).unwrap().0, -change);
    }
    for invalid in [279.0, 361.0, 0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for (reference, temperature) in [(invalid, 320.0), (320.0, invalid), (invalid, invalid)] {
            assert_eq!(law.relative_enthalpy_and_heat_capacity(reference, temperature),
                Err(GasEnergyError::Invalid("temperature outside the declared caloric bracket")));
        }
    }
}

#[test]
fn relative_caloric_sample_preserves_small_heat_under_large_formation_references() {
    use crate::thermochemical_gas::tests::{species_with_energy_reference, transport};

    let reference = 320.0;
    let temperature = reference + 2.0_f64.powi(-16);
    let mut original = None;
    for formation in [0.0, 2.0_f64.powi(70), -2.0_f64.powi(70)] {
        let source = ThermochemicalGasModel::new(
            fs_thermochem::mixture::FrozenIdealGasMixtureModelV1::new(
                vec![species_with_energy_reference("A", 0.028, 3.5, 0.001, formation)],
                mole(vec![1.0]),
            ).unwrap(), vec![transport("A", 1.7e-5, 0.025)], None,
        ).unwrap();
        let law = GasEnthalpyLaw::from_mole(&source, mole(vec![1.0]),
            point(reference, 100000.0), controls()).unwrap();
        let sample = law.relative_enthalpy_and_heat_capacity(reference, temperature).unwrap();
        let expected = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.028
            * (temperature - reference) * (3.5 + 0.0005 * (reference + temperature));
        assert!(sample.0 > 0.0);
        assert!((sample.0 - expected).abs() < 2e-14 * expected);
        if let Some(original) = original {
            assert_eq!(sample, original);
        } else {
            original = Some(sample);
        }
        if formation != 0.0 {
            let rounded_change = law.enthalpy_at_temperature(temperature).unwrap()
                - law.enthalpy_at_temperature(reference).unwrap();
            assert!((rounded_change - expected).abs() > 0.5 * expected);
        }
    }
}

#[test]
fn thermal_entropy_adapter_integrates_original_variable_cp_at_fixed_pressure() {
    let source = model(vec![0.75,0.25],0.001);
    let law = GasEnthalpyLaw::from_mole(&source,mole(vec![0.75,0.25]),
        point(300.0,140000.0),controls()).unwrap();
    let h0 = law.enthalpy_at_temperature(300.0).unwrap();
    let h1 = law.enthalpy_at_temperature(337.0).unwrap();
    let entropy = law.specific_entropy_change_j_kg_k(h0,h1).unwrap().unwrap();
    let molar_mass = 0.75*0.028 + 0.25*0.004;
    let exact = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / molar_mass
        * (3.25*(337.0_f64/300.0).ln() + 0.00075*37.0);
    assert!((entropy-exact).abs() < 1e-9);
    assert_eq!(law.specific_entropy_change_j_kg_k(h1,h0).unwrap().unwrap(),-entropy);
    assert_eq!(law.specific_entropy_change_j_kg_k(h0,h0).unwrap(),Some(0.0));
    assert!(law.specific_entropy_change_j_kg_k(h0-1e9,h1).is_err());
}

fn entropy_shifted_source(offset: f64) -> ThermochemicalGasModel {
    let source = model(vec![0.5, 0.5], 0.001);
    let components = source.thermodynamics().components().iter().enumerate()
        .map(|(index, original)| {
            let mut card = original.regions()[0].card().clone();
            card.parameters.get_mut("a8").unwrap().value = if index == 0 { offset } else { -offset };
            Nasa9StandardStateModelV1::new(original.species().clone(), original.molar_mass(),
                original.convention().clone(),
                vec![fs_thermochem::Nasa9RegionV1::from_card(card).unwrap()]).unwrap()
        }).collect();
    ThermochemicalGasModel::new(
        fs_thermochem::mixture::FrozenIdealGasMixtureModelV1::new(components,
            mole(vec![0.5, 0.5])).unwrap(), source.transport_materials().to_vec(), None,
    ).unwrap()
}

#[test]
fn source_entropy_anchor_retains_ideal_mixing_pressure_and_compatible_reference() {
    let source = model(vec![0.75, 0.25], 0.001);
    let prepare = |pressure, controls| GasEnthalpyLaw::from_mole(&source,
        mole(vec![0.75, 0.25]), point(320.0, pressure), controls).unwrap();
    let standard = prepare(100000.0, controls());
    let compressed = prepare(140000.0, controls());
    let a = standard.advected_entropy_anchor().unwrap();
    let b = compressed.advected_entropy_anchor().unwrap();
    let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / (0.75 * 0.028 + 0.25 * 0.004);
    let mixing = -r * (0.75 * 0.75_f64.ln() + 0.25 * 0.25_f64.ln());
    assert!((a.reduced_specific_entropy_j_kg_k - mixing).abs() < 1e-11);
    assert!((b.reduced_specific_entropy_j_kg_k - mixing + r * 1.4_f64.ln()).abs() < 1e-11);
    assert_eq!(a.entropy_reference, b.entropy_reference);
    assert_eq!(a.specific_h_j_kg, standard.enthalpy_at_temperature(280.0).unwrap());
    let other_anchor = prepare(100000.0,
        GasEnergySolve { minimum_temperature_k: 281.0, ..controls() });
    assert_ne!(a.entropy_reference, other_anchor.advected_entropy_anchor().unwrap().entropy_reference);
}

#[test]
fn native_stationary_mixture_entropy_uses_source_calorics_under_species_reference_shifts() {
    use fs_thermal::enthalpy::advection::{Advected, MassBoundary, species::InletSpecies};
    let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
    for offset in [0.0, 2.0_f64.powi(50), -2.0_f64.powi(50)] {
        let source = entropy_shifted_source(offset);
        let law = |fractions| GasEnthalpyLaw::new(&source,
            Composition::new(CompositionBasis::MassFraction, fractions).unwrap(),
            point(320.0, 100000.0), controls()).unwrap();
        let mixed = law(vec![0.5, 0.5]);
        let a = law(vec![1.0, 0.0]);
        let b = law(vec![0.0, 1.0]);
        for (ta, tb) in [(320.0, 320.0), (300.0, 340.0)] {
            let ha = a.enthalpy_at_temperature(ta).unwrap();
            let hb = b.enthalpy_at_temperature(tb).unwrap();
            let h = 0.5 * ha + 0.5 * hb;
            let temperature = mixed.temperature_k(h).unwrap();
            let initial = thermal_timeline(&[
                (MaterialId(1), &mixed), (MaterialId(2), &a), (MaterialId(3), &b),
            ], h, 0.0);
            let flow = Advected::new((), vec![], 0).unwrap().with_boundaries(vec![
                MassBoundary::Inlet { id: 1, cell: CellId(1), kg_s: 0.1,
                    material_id: MaterialId(2), specific_h_j_kg: ha },
                MassBoundary::Inlet { id: 2, cell: CellId(1), kg_s: 0.1,
                    material_id: MaterialId(3), specific_h_j_kg: hb },
                MassBoundary::Outlet { id: 3, cell: CellId(1), kg_s: 0.2 },
            ], 3).unwrap().with_stationary_species(vec!["A".into(), "B".into()], vec![
                InletSpecies { boundary_id: 1, mass_fractions: vec![1.0, 0.0] },
                InletSpecies { boundary_id: 2, mass_fractions: vec![0.0, 1.0] },
            ], 1, 2).unwrap();
            let (next, observed, energy, entropy) = initial.advance_advected_with_entropy(
                &flow, 0.25, thermal_config(), || Ok(())).unwrap();
            // Equal mass inflows give x_A=1/8, x_B=7/8. Exact NASA Cp/T
            // integrals and ideal mixing independently determine the production.
            let expected_rate = r * (0.1 / 0.028
                * (3.5 * (temperature / ta).ln() + 0.001 * (temperature - ta) - 0.125_f64.ln())
                + 0.1 / 0.004 * (2.5 * (temperature / tb).ln() - 0.875_f64.ln()));
            let advection = entropy.advection.unwrap();
            assert!(expected_rate > 0.0);
            assert!((advection.upwind_mixing_j_k - 0.25 * expected_rate).abs() < 2e-9);
            assert!((advection.boundary_entropy_input_j_k + 0.25 * expected_rate).abs() < 2e-9);
            assert!(entropy.thermal_storage_change_j_k.abs() < 1e-10);
            assert!(entropy.temporal_defect_j_k.abs() < 1e-10);
            assert!(energy.external_heat_j.abs() < 1e-8);
            assert!((next.temperature_k(CellId(1)).unwrap() - temperature).abs() < 1e-9);
            assert!(observed.entropy_audit().is_some());
            assert!(!entropy.has_unclassified_sources);
            assert!(matches!(initial.advance_advected_with_entropy(&flow, 0.25,
                thermal_config(), || Err(DepositionError::Cancelled)), Err(DepositionError::Cancelled)));
            assert_eq!(initial.cells()[&CellId(1)].specific_h_j_kg, h);
        }
    }
}

#[test]
fn thermal_entropy_retains_heat_smaller_than_the_caloric_inverse_tolerance() {
    let source = model(vec![0.75, 0.25], 0.001);
    for tolerance in [1e-8, 1e-3] {
        let law = GasEnthalpyLaw::from_mole(&source, mole(vec![0.75, 0.25]),
            point(320.0, 140000.0), GasEnergySolve {
                energy_tolerance_j_kg: tolerance, ..controls()
            }).unwrap();
        let h0 = law.enthalpy_at_temperature(320.0).unwrap();
        let h1 = h0 + 0.25 * tolerance;
        assert!(h1 > h0);
        // Both inversions legitimately accept the same bracket midpoint.
        assert_eq!(law.temperature_k(h0).unwrap(), 320.0);
        assert_eq!(law.temperature_k(h1).unwrap(), 320.0);
        let entropy = law.specific_entropy_change_j_kg_k(h0, h1).unwrap().unwrap();
        assert!(entropy > 0.0);
        assert_eq!(entropy, (h1 - h0) / 320.0);
        assert_eq!(law.specific_entropy_change_j_kg_k(h1, h0).unwrap(), Some(-entropy));
        assert_eq!(law.specific_entropy_change_j_kg_k(h0, h0).unwrap(), Some(0.0));
    }
}

#[test]
fn gas_heating_timeline_accounts_for_heat_below_temperature_reconstruction_tolerance() {
    let source = model(vec![0.75, 0.25], 0.001);
    let law = GasEnthalpyLaw::from_mole(&source, mole(vec![0.75, 0.25]),
        point(320.0, 140000.0), GasEnergySolve {
            energy_tolerance_j_kg: 0.01, ..controls()
        }).unwrap();
    let h0 = law.enthalpy_at_temperature(320.0).unwrap();
    // One millijoule over one second in a one-kilogram isobaric cell.
    let initial = thermal_timeline(&[(MaterialId(1), &law)], h0, 0.001);
    let (heated, energy, entropy) = initial.advance_to_with_entropy(
        1.0, thermal_config(), || Ok(())).unwrap();
    let stored_heat = heated.cells()[&CellId(1)].specific_h_j_kg - h0;
    assert_eq!(energy.substeps, 1);
    assert_eq!(energy.external_heat_j, 0.001);
    assert!(energy.residual_j.abs() < thermal_config().energy_tolerance_j);
    assert!(stored_heat > 0.0);
    assert_eq!(law.temperature_k(heated.cells()[&CellId(1)].specific_h_j_kg).unwrap(), 320.0);
    assert_eq!(entropy.thermal_storage_change_j_k, stored_heat / 320.0);
    assert_eq!(entropy.unclassified_source_at_cell_j_k, 0.001 / 320.0);
    assert!(entropy.has_unclassified_sources);
    assert_eq!(entropy.conduction_generation_j_k, 0.0);
    assert!(entropy.temporal_defect_j_k.abs() < 3e-13);
    assert_eq!(initial.cells()[&CellId(1)].specific_h_j_kg, h0);
}

#[test]
fn zero_species_and_shared_reference_survive_caloric_basis_conversion() {
    let source = model(vec![0.75, 0.25], 0.0);
    let a = GasEnthalpyLaw::from_mole(&source, mole(vec![1.0, 0.0]),
        point(300.0, 101325.0), controls()).unwrap();
    let b = GasEnthalpyLaw::from_mole(&source, mole(vec![0.0, 1.0]),
        point(300.0, 101325.0), controls()).unwrap();
    assert_eq!(a.mole_composition().fractions(), &[1.0, 0.0]);
    assert_eq!(a.advected_composition().unwrap().mass_fractions, &[1.0, 0.0]);
    assert_eq!(a.advected_composition().unwrap().species, b.advected_composition().unwrap().species);
    assert_eq!(a.advected_composition().unwrap().enthalpy_reference,
        b.advected_composition().unwrap().enthalpy_reference);
    let h = a.enthalpy_at_temperature(320.0).unwrap();
    assert!(b.temperature_k(h).is_err());
    let mixed = GasEnthalpyLaw::new(&source,
        Composition::new(CompositionBasis::MassFraction, vec![0.25, 0.75]).unwrap(),
        point(300.0, 101325.0), controls()).unwrap();
    let actual = mixed.model().evaluate(&point(320.0, 101325.0)).unwrap();
    let expected_h = 320.0 * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
        * (0.25 * 3.5 / 0.028 + 0.75 * 2.5 / 0.004);
    assert!((actual.specific_enthalpy_j_kg() - expected_h).abs() < 1e-8);
}

#[test]
fn caloric_domain_work_and_cancellation_refuse_without_publishing_a_state() {
    let source = model(vec![0.75, 0.25], 0.001);
    let law = GasEnthalpyLaw::from_mole(&source, mole(vec![0.75, 0.25]),
        point(300.0, 101325.0), controls()).unwrap();
    let h = law.enthalpy_at_temperature(337.0).unwrap();
    assert!(matches!(law.resolve_enthalpy(h, || true), Err(GasEnergyError::Cancelled)));
    assert!(law.resolve_enthalpy(law.enthalpy_at_temperature(360.0).unwrap() + 1.0,
        || false).is_err());
    assert!(GasEnthalpyLaw::from_mole(&source, mole(vec![0.75, 0.25]),
        point(300.0, 101325.0),
        GasEnergySolve { minimum_temperature_k: 199.0, ..controls() }).is_err());
    assert!(GasEnthalpyLaw::from_mole(&source, mole(vec![0.75, 0.25]),
        point(300.0, 101325.0),
        GasEnergySolve { maximum_iterations: 0, ..controls() }).is_err());
    assert!(law.enthalpy_at_temperature(279.0).is_err());
}

fn coefficient_species(coefficients: [f64; 9], temperatures: [f64; 2]) -> Nasa9StandardStateModelV1 {
    let original = crate::thermochemical_gas::tests::species_in_domain(
        "N2", 0.0280134, 3.5, 0.0, 0.0, temperatures,
    );
    let mut card = original.regions()[0].card().clone();
    for (index, value) in coefficients.into_iter().enumerate() {
        card.parameters.get_mut(&format!("a{index}")).unwrap().value = value;
    }
    Nasa9StandardStateModelV1::new(original.species().clone(), original.molar_mass(),
        original.convention().clone(), vec![fs_thermochem::Nasa9RegionV1::from_card(card).unwrap()]).unwrap()
}

// Read the actual checked-in NASA source coefficients, not a second fitted law.
fn nitrogen_coefficients(region: &str) -> [f64; 9] {
    let text = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../data/matdb/seed-v1/nitrogen-nasa9-mcbride2002/nasa9.tsv"));
    let mut coefficients = [0.0; 9];
    let mut count = 0;
    for line in text.lines() {
        let fields: Vec<_> = line.split('\t').collect();
        if fields.first() == Some(&"coefficient") && fields[2] == region {
            let index: usize = fields[3].strip_prefix('a').unwrap().parse().unwrap();
            coefficients[index] = fields[4].parse().unwrap();
            count += 1;
        }
    }
    assert_eq!(count, 9);
    coefficients
}

fn nitrogen_model(region: &str) -> ThermochemicalGasModel {
    let temperatures = if region == "low" { [200.0, 1000.0] } else { [1000.0, 6000.0] };
    ThermochemicalGasModel::new(
        fs_thermochem::mixture::FrozenIdealGasMixtureModelV1::new(
            vec![coefficient_species(nitrogen_coefficients(region), temperatures)],
            mole(vec![1.0]),
        ).unwrap(),
        vec![crate::thermochemical_gas::tests::transport_in_domain(
            "N2", 1.8e-5, 0.03, [200.0, 6000.0], [1000.0, 1.0e7],
        )], None,
    ).unwrap()
}

fn hot_controls(lo: f64, hi: f64) -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: lo, maximum_temperature_k: hi,
        energy_tolerance_j_kg: 1e-6, maximum_iterations: 128 }
}

#[test]
fn nasa9_nitrogen_capacity_bounds_cover_real_cancelling_polynomials() {
    for (region, lo, hi) in [("low", 280.0, 900.0), ("high", 1000.0, 6000.0)] {
        let gas = nitrogen_model(region);
        let species = &gas.thermodynamics().components()[0];
        let controls = hot_controls(lo, hi);
        // This is the previous failing whole-bracket coefficient enclosure.
        assert!(capacity_interval(&species.regions()[0], species.molar_mass().value(), lo, hi).lo < 0.0);
        let cp_lower = minimum_cp(species, controls).unwrap();
        let law = GasEnthalpyLaw::from_mole(&gas, mole(vec![1.0]),
            point(lo, 100000.0), controls).unwrap();
        let isochoric_slope = gas.maximum_isochoric_temperature_slope(controls).unwrap();
        let mut sampled_min = f64::INFINITY;
        for index in 0..=64 {
            let t = lo + (hi - lo) * index as f64 / 64.0;
            let state = gas.evaluate(&point(t, 100000.0)).unwrap();
            let properties = state.thermodynamics().properties();
            let cp = properties.mass_specific_cp().value();
            let cv = properties.mass_specific_cv().value();
            sampled_min = sampled_min.min(cp);
            assert!(cp_lower <= cp);
            assert!(law.maximum_temperature_slope() >= 1.0 / cp);
            assert!(isochoric_slope >= 1.0 / cv);
        }
        // A useful bound avoids the old false refusal and arbitrary tiny steps.
        assert!(cp_lower > 0.8 * sampled_min);
        assert_eq!(cp_lower.to_bits(), minimum_cp(species, controls).unwrap().to_bits());
    }
}

fn nitrogen_specific_enthalpy(t: f64) -> f64 {
    let a = nitrogen_coefficients("low");
    // Independent integrated polynomial, rather than T times h/(RT).
    let h_over_r = -a[0] / t + a[1] * t.ln() + a[2] * t
        + a[3] * t * t / 2.0 + a[4] * t * t * t / 3.0
        + a[5] * t * t * t * t / 4.0 + a[6] * t * t * t * t * t / 5.0 + a[7];
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * h_over_r / 0.0280134
}

#[test]
fn nitrogen_heating_closes_integrated_energy_with_admitted_capacity() {
    let gas = nitrogen_model("low");
    let controls = hot_controls(280.0, 900.0);
    gas.maximum_isochoric_temperature_slope(controls).unwrap();
    let law = GasEnthalpyLaw::from_mole(&gas, mole(vec![1.0]),
        point(300.0, 100000.0), controls).unwrap();
    let h300 = nitrogen_specific_enthalpy(300.0);
    let h850 = nitrogen_specific_enthalpy(850.0);
    let solved = law.resolve_enthalpy(h850, || false).unwrap();
    assert!((solved.state().state().temperature - 850.0).abs() < 1e-8);
    let cp300 = gas.evaluate(&point(300.0, 100000.0)).unwrap().state().specific_heat_cp;
    assert!((h850 - h300 - cp300 * 550.0).abs() > 1000.0);
    let mut parcel = super::super::IsochoricGasParcel::new(gas, point(300.0, 100000.0), 0.01).unwrap();
    let mass = parcel.mass_kg();
    let heat = mass * (h850 - h300 - UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.0280134 * 550.0);
    let report = parcel.add_heat(heat, controls, || false).unwrap();
    assert!((parcel.state().state().temperature - 850.0).abs() < 1e-8);
    assert!((parcel.state().state().pressure - 100000.0 * 850.0 / 300.0).abs() < 1e-5);
    assert!(report.residual_j.abs() <= mass * controls.energy_tolerance_j_kg);
    assert!((report.stored_energy_change_j - heat).abs() <= mass * controls.energy_tolerance_j_kg);
}

#[test]
fn capacity_bound_refuses_nonpositive_cv_including_interior_dips() {
    for cp_over_r in [0.5, 1.0] {
        let mut coefficients = [0.0; 9];
        coefficients[2] = cp_over_r;
        let species = coefficient_species(coefficients, [200.0, 1000.0]);
        assert!(minimum_cp(&species, hot_controls(300.0, 900.0)).is_err());
    }
    // Cp/R = 0.9 + 1e-5*(T-600)^2: both endpoints have positive cv,
    // while the interior violates the ideal-gas caloric stability condition.
    let mut coefficients = [0.0; 9];
    coefficients[2] = 4.5;
    coefficients[3] = -0.012;
    coefficients[4] = 1e-5;
    let species = coefficient_species(coefficients, [200.0, 1000.0]);
    assert!(minimum_cp(&species, hot_controls(300.0, 900.0)).is_err());
}

fn multiregion_model(low_coefficients: [f64; 9], high_coefficients: [f64; 9]) -> ThermochemicalGasModel {
    let low = coefficient_species(low_coefficients, [200.0, 1000.0]);
    let high = coefficient_species(high_coefficients, [1000.0, 6000.0]);
    let species = Nasa9StandardStateModelV1::new(low.species().clone(), low.molar_mass(),
        low.convention().clone(), vec![low.regions()[0].clone(), high.regions()[0].clone()]).unwrap();
    ThermochemicalGasModel::new(
        fs_thermochem::mixture::FrozenIdealGasMixtureModelV1::new(vec![species], mole(vec![1.0])).unwrap(),
        vec![crate::thermochemical_gas::tests::transport_in_domain(
            "N2", 1.8e-5, 0.03, [200.0, 6000.0], [1000.0, 1e7],
        )], None,
    ).unwrap()
}

#[test]
fn relative_caloric_sample_integrates_reconciled_regions_with_the_endpoint_cp() {
    let mut low = [0.0; 9]; low[2] = 3.5;
    let mut high = [0.0; 9]; high[2] = 4.5;
    let raw = multiregion_model(low, high);
    let gas = raw.with_caloric_continuity(Temperature::new(300.0),
        fs_thermochem::MolarEnergyQuantityV1::new(100000.0),
        fs_thermochem::MolarThermalQuantityV1::new(1000.0)).unwrap();
    let law = GasEnthalpyLaw::from_mole(&gas, mole(vec![1.0]),
        point(300.0, 100000.0), hot_controls(280.0, 3000.0)).unwrap();
    let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.0280134;
    for (reference, temperature, integral, cp_over_r) in [
        (900.0, 1000.0, 350.0, 4.5), (1000.0, 1100.0, 450.0, 4.5),
        (900.0, 1100.0, 800.0, 4.5), (1100.0, 900.0, -800.0, 3.5),
        (1000.0, 1000.0, 0.0, 4.5),
    ] {
        let (change, capacity) = law.relative_enthalpy_and_heat_capacity(reference, temperature).unwrap();
        assert!((change - r * integral).abs() <= 2e-14 * (r * integral).abs());
        assert!((capacity - r * cp_over_r).abs() <= 2e-14 * r * cp_over_r);
    }
}

#[test]
fn reconciled_nitrogen_heating_crosses_source_regions_with_one_energy_owner() {
    let raw = multiregion_model(nitrogen_coefficients("low"), nitrogen_coefficients("high"));
    let controls = hot_controls(280.0, 3000.0);
    assert!(GasEnthalpyLaw::from_mole(&raw, mole(vec![1.0]),
        point(300.0, 100000.0), controls).is_err());
    let gas = raw.with_caloric_continuity(Temperature::new(298.15),
        fs_thermochem::MolarEnergyQuantityV1::new(1e-3),
        fs_thermochem::MolarThermalQuantityV1::new(1e-5)).unwrap();
    let law = GasEnthalpyLaw::from_mole(&gas, mole(vec![1.0]),
        point(300.0, 100000.0), controls).unwrap();
    let cv_slope = gas.maximum_isochoric_temperature_slope(controls).unwrap();
    for t in [280.0, 900.0, 999.999, 1000.0, 1000.001, 2500.0, 3000.0] {
        let h = law.enthalpy_at_temperature(t).unwrap();
        let solved = law.resolve_enthalpy(h, || false).unwrap();
        assert!((solved.state().state().temperature - t).abs() < 1e-8);
        assert!(solved.residual_j_kg().abs() <= controls.energy_tolerance_j_kg);
        let properties = solved.state().thermodynamics().properties();
        assert!(law.maximum_temperature_slope() >= 1.0 / properties.mass_specific_cp().value());
        assert!(cv_slope >= 1.0 / properties.mass_specific_cv().value());
        assert!(solved.state().thermodynamics().receipt().component_receipts()[0]
            .caloric_reference().is_some());
    }
    let h0 = law.enthalpy_at_temperature(300.0).unwrap();
    let h1 = law.enthalpy_at_temperature(2500.0).unwrap();
    let entropy = law.specific_entropy_change_j_kg_k(h0, h1).unwrap().unwrap();
    let expected_entropy = gas.thermodynamics().specific_entropy_change(
        Temperature::new(300.0), fs_qty::Pressure::new(100000.0),
        Temperature::new(2500.0), fs_qty::Pressure::new(100000.0)).unwrap().value();
    assert!((entropy - expected_entropy).abs() < 1e-8);
    let mut parcel = super::super::IsochoricGasParcel::new(gas, point(300.0, 100000.0), 0.01).unwrap();
    let mass = parcel.mass_kg();
    let heat = mass * (h1 - h0 - UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.0280134 * 2200.0);
    let report = parcel.add_heat(heat, controls, || false).unwrap();
    assert!((parcel.state().state().temperature - 2500.0).abs() < 1e-8);
    assert!((parcel.state().state().pressure - 100000.0 * 2500.0 / 300.0).abs() < 1e-4);
    assert!(report.residual_j.abs() <= mass * controls.energy_tolerance_j_kg);
    assert!((report.stored_energy_change_j - heat).abs() <= mass * controls.energy_tolerance_j_kg);
}

#[test]
fn caloric_transport_reference_binds_policy_but_preserves_the_full_species_catalog() {
    let raw = model(vec![0.75, 0.25], 0.001);
    let corrected = raw.with_caloric_continuity(Temperature::new(300.0),
        fs_thermochem::MolarEnergyQuantityV1::new(0.0),
        fs_thermochem::MolarThermalQuantityV1::new(0.0)).unwrap();
    let law = |source: &ThermochemicalGasModel, fractions: Vec<f64>| {
        GasEnthalpyLaw::from_mole(source, mole(fractions), point(300.0, 100000.0), controls()).unwrap()
    };
    let original = law(&raw, vec![1.0, 0.0]);
    let a = law(&corrected, vec![1.0, 0.0]);
    let b = law(&corrected, vec![0.0, 1.0]);
    assert_eq!(a.enthalpy_at_temperature(330.0).unwrap(), original.enthalpy_at_temperature(330.0).unwrap());
    assert_ne!(a.advected_composition().unwrap().enthalpy_reference,
        original.advected_composition().unwrap().enthalpy_reference);
    assert_eq!(a.advected_composition().unwrap().enthalpy_reference,
        b.advected_composition().unwrap().enthalpy_reference);
}

#[test]
fn stationary_advection_refuses_unreconciled_source_formation_reference_mismatch() {
    use fs_thermal::enthalpy::advection::{Advected, MassBoundary, species::InletSpecies};
    use crate::thermochemical_gas::tests::{species_with_energy_reference, transport};

    let source = |reference_j_kg| ThermochemicalGasModel::new(
        fs_thermochem::mixture::FrozenIdealGasMixtureModelV1::new(
            vec![species_with_energy_reference("A", 0.028, 3.5, 0.001, reference_j_kg)],
            mole(vec![1.0]),
        ).unwrap(), vec![transport("A", 1.7e-5, 0.025)], None,
    ).unwrap();
    let cell_source = source(0.0);
    let cell_law = GasEnthalpyLaw::from_mole(&cell_source, mole(vec![1.0]),
        point(320.0, 100000.0), controls()).unwrap();
    let h0 = cell_law.enthalpy_at_temperature(320.0).unwrap();
    for reference_j_kg in [0.0, 1000.0] {
        let inlet_source = source(reference_j_kg);
        let inlet_law = GasEnthalpyLaw::from_mole(&inlet_source, mole(vec![1.0]),
            point(320.0, 100000.0), controls()).unwrap();
        // The species, composition and named convention agree; only the
        // source's formation-energy coefficient changes. Without source
        // binding the 1000 J/kg offset acts as spurious sensible inlet heat.
        assert_eq!(cell_source.thermodynamics().components()[0].convention(),
            inlet_source.thermodynamics().components()[0].convention());
        assert!(inlet_source.thermodynamics().components()[0].caloric_reference().is_none());
        let inlet_h = inlet_law.enthalpy_at_temperature(320.0).unwrap();
        assert!((inlet_h - h0 - reference_j_kg).abs() < 1e-8);
        let initial = thermal_timeline(
            &[(MaterialId(1), &cell_law), (MaterialId(2), &inlet_law)], h0, 0.0);
        let flow = Advected::new((), vec![], 0).unwrap().with_boundaries(vec![
            MassBoundary::Inlet { id: 1, cell: CellId(1), kg_s: 0.001,
                material_id: MaterialId(2), specific_h_j_kg: inlet_h },
            MassBoundary::Outlet { id: 2, cell: CellId(1), kg_s: 0.001 },
        ], 2).unwrap().with_stationary_species(vec!["A".into()],
            vec![InletSpecies { boundary_id: 1, mass_fractions: vec![1.0] }], 1, 1).unwrap();
        let result = initial.advance_observed(&flow, 1.0, thermal_config(), || Ok(()));
        if reference_j_kg == 0.0 {
            let (advanced, _, audit) = result.unwrap();
            assert_eq!(advanced.cells()[&CellId(1)].specific_h_j_kg, h0);
            assert_eq!(audit.residual_j, 0.0);
        } else {
            assert_eq!(result.err(), Some(DepositionError::InvalidInput(
                "stationary species caloric laws have different enthalpy references")));
        }
        assert_eq!(initial.cells()[&CellId(1)].specific_h_j_kg, h0);
    }
}

#[test]
fn reconciled_regions_still_refuse_nonpositive_cv_at_the_selected_seam_endpoint() {
    let mut low = [0.0; 9]; low[2] = 3.5;
    let mut high = [0.0; 9]; high[2] = 0.5;
    // Match the raw h/s reference at 1000 K so this exercises Cv admission.
    high[7] = 3000.0;
    high[8] = 3.0 * fs_math::det::ln(1000.0);
    let raw = multiregion_model(low, high);
    let gas = raw.with_caloric_continuity(Temperature::new(300.0),
        fs_thermochem::MolarEnergyQuantityV1::new(100000.0),
        fs_thermochem::MolarThermalQuantityV1::new(1000.0)).unwrap();
    assert!(gas.maximum_isochoric_temperature_slope(hot_controls(280.0, 999.0)).is_ok());
    assert!(gas.maximum_isochoric_temperature_slope(hot_controls(280.0, 1000.0)).is_err());
    assert!(GasEnthalpyLaw::from_mole(&gas, mole(vec![1.0]), point(300.0, 100000.0),
        hot_controls(280.0, 3000.0)).is_err());
}
