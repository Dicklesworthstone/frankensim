//! Synthetic NASA-9 coefficients provide an independent caloric oracle.
use super::*;
use crate::thermochemical_gas::tests::{model, point};

fn controls() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 280.0, maximum_temperature_k: 360.0,
        energy_tolerance_j_kg: 1e-8, maximum_iterations: 64 }
}
fn mole(fractions: Vec<f64>) -> Composition {
    Composition::new(CompositionBasis::MoleFraction, fractions).unwrap()
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
