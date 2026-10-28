//! Synthetic analytic fixtures, not empirical accuracy claims for named gases.
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy,
    LawId, LawParameter, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, UncertaintyModel};
use fs_qty::MolarMass;
use fs_thermochem::{Composition, CompositionBasis, ElementalReferenceIdV1,
    NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1, NASA9_STATE_SCHEMA_VERSION_V1,
    Nasa9RegionV1, Nasa9StandardStateModelV1, ReferenceEquationOfStateV1,
    StandardStateConventionV1, StandardStatePhaseV1, required_nasa9_parameter_dims,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

#[test]
fn generated_mole_amounts_preserve_sources_and_extensive_scaling() {
    let bank = model(vec![0.5, 0.5], 0.001);
    let expected = bank.with_composition(
        Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
    ).unwrap();
    for scale in [1.0e-280, 1.0, 1.0e280] {
        let gas = bank.with_mole_amounts(&[
            fs_qty::Amount::new(3.0 * scale), fs_qty::Amount::new(scale),
        ]).unwrap();
        assert_eq!(gas, expected);
    }
    // A normalized constitutive composition does not need a representable
    // total amount; the largest-first scaling avoids an overflowing sum.
    assert_eq!(bank.with_mole_amounts(&[
        fs_qty::Amount::new(f64::MAX), fs_qty::Amount::new(f64::MAX),
    ]).unwrap(), bank);
    let pure = bank.with_mole_amounts(&[
        fs_qty::Amount::new(0.0), fs_qty::Amount::new(2.0),
    ]).unwrap();
    assert_eq!(pure.thermodynamics().components(), &bank.thermodynamics().components()[1..]);
    assert_eq!(pure.transport_materials(), &bank.transport_materials()[1..]);
}

#[test]
fn generated_mole_amounts_refuse_lost_traces_and_invalid_inventories() {
    let bank = model(vec![0.5, 0.5], 0.001);
    for amounts in [vec![], vec![1.0], vec![0.0, 0.0], vec![-1.0, 2.0],
        vec![f64::NAN, 1.0], vec![f64::INFINITY, 1.0], vec![f64::from_bits(1), f64::MAX]]
    {
        let amounts: Vec<_> = amounts.into_iter().map(fs_qty::Amount::new).collect();
        assert!(bank.with_mole_amounts(&amounts).is_err());
    }
}

pub(crate) fn species(name: &str, mass: f64, cp_r: f64, slope: f64) -> Nasa9StandardStateModelV1 {
    species_with_energy_reference(name, mass, cp_r, slope, 0.0)
}

pub(crate) fn species_with_energy_reference(
    name: &str,
    mass: f64,
    cp_r: f64,
    slope: f64,
    energy_reference_j_kg: f64,
) -> Nasa9StandardStateModelV1 {
    species_in_domain(
        name,
        mass,
        cp_r,
        slope,
        energy_reference_j_kg,
        [200.0, 2000.0],
    )
}

pub(crate) fn species_in_domain(
    name: &str,
    mass: f64,
    cp_r: f64,
    slope: f64,
    energy_reference_j_kg: f64,
    temperatures: [f64; 2],
) -> Nasa9StandardStateModelV1 {
    let mut coefficients = [0.0; 9];
    coefficients[2] = cp_r;
    coefficients[3] = slope;
    // NASA-9's eighth coefficient adds R*a8 to the molar enthalpy.
    coefficients[7] = mass * energy_reference_j_kg / R;
    let mut parameters = std::collections::BTreeMap::new();
    for ((key, dims), value) in required_nasa9_parameter_dims().into_iter().take(9).zip(coefficients) {
        parameters.insert(key.to_owned(), LawParameter { value, dims });
    }
    parameters.insert("reference_pressure".into(), LawParameter { value: 100000.0, dims: Pressure::DIMS });
    let card = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()), law_version: NASA9_LAW_VERSION_V1,
        parameters, state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", temperatures[0], temperatures[1]),
        sources: vec![], provenance: Provenance { source: "synthetic NASA-9 gas".into(),
            license: "CC0-1.0".into(), artifact: None },
    };
    Nasa9StandardStateModelV1::new(SpeciesId::new(name).unwrap(), MolarMass::new(mass),
        StandardStateConventionV1::new(StandardStatePhaseV1::Gas, ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0), ElementalReferenceIdV1::new("synthetic-reference").unwrap()).unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()]).unwrap()
}

pub(crate) fn transport(name: &str, mu: f64, k: f64) -> SpeciesTransportMaterial {
    transport_in_domain(name, mu, k, [200.0, 2000.0], [1000.0, 1.0e7])
}

pub(crate) fn transport_in_domain(
    name: &str,
    mu: f64,
    k: f64,
    temperatures: [f64; 2],
    pressures: [f64; 2],
) -> SpeciesTransportMaterial {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [("dynamic_viscosity", DynViscosity::DIMS, mu),
        (THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, k)] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve { abscissa: "temperature".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(temperatures[0], value), (temperatures[1], 2.0 * value)], dims },
            validity: ValidityDomain::unconstrained()
                .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), temperatures[0], temperatures[1])
                .with_quantity("pressure", typed(QuantityKind::Pressure), pressures[0], pressures[1]),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside, observations: vec![],
            provenance: Provenance { source: "synthetic pure transport curve".into(),
                license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
    }
    SpeciesTransportMaterial { species: SpeciesId::new(name).unwrap(),
        material: MaterialCard::assemble(MaterialStateId { chemistry: name.into(),
            phase: "gas".into(), process: "synthetic".into(), revision: 0 }, claims, vec![]).unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly }
}

pub(crate) fn model(fractions: Vec<f64>, slope: f64) -> ThermochemicalGasModel {
    model_in_domains(
        fractions,
        slope,
        [200.0, 2000.0],
        [200.0, 2000.0],
        [1000.0, 1.0e7],
    )
}

pub(crate) fn model_in_domains(
    fractions: Vec<f64>,
    slope: f64,
    thermodynamic_temperatures: [f64; 2],
    transport_temperatures: [f64; 2],
    pressures: [f64; 2],
) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(
        vec![species_in_domain("A", 0.028, 3.5, slope, 0.0, thermodynamic_temperatures),
            species_in_domain("B", 0.004, 2.5, 0.0, 0.0, thermodynamic_temperatures)],
        Composition::new(CompositionBasis::MoleFraction, fractions).unwrap()).unwrap(),
        vec![transport_in_domain("B", 2.0e-5, 0.15, transport_temperatures, pressures),
            transport_in_domain("A", 1.7e-5, 0.025, transport_temperatures, pressures)], None).unwrap()
}
pub(crate) fn point(t: f64, p: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), t).unwrap()
        .with_quantity("pressure", typed(QuantityKind::Pressure), p).unwrap()
}
pub(super) fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 2e-12 * a.abs().max(b.abs()).max(1e-100), "{a} != {b}");
}

#[test]
fn g0_gamma_cp_coherence() {
    let model = model(vec![0.75, 0.25], 0.001);
    for t in [200.0, 293.15, 308.15, 600.0, 1000.0, 2000.0] {
        let resolved = model.evaluate(&point(t, 101325.0)).unwrap();
        let gas = resolved.state();
        let properties = resolved.thermodynamics().properties();
        let cp = properties.mass_specific_cp().value();
        let cv = properties.mass_specific_cv().value();
        assert!(((cp - cv) / gas.specific_gas_constant - 1.0).abs() <= 1e-14);
        assert!((gas.gamma / (cp / cv) - 1.0).abs() <= 1e-14);
        gas.validate_coherence().unwrap();
    }
}

#[test]
fn source_domains_govern_hot_gas_acoustics_and_transport_together() {
    // Manufactured NASA-9/transport data, not a hot-air or real-gas accuracy claim.
    let model = model_in_domains(
        vec![0.75, 0.25],
        0.001,
        [200.0, 6000.0],
        [200.0, 5000.0],
        [1000.0, 3.0e7],
    );
    let resolved = model.evaluate(&point(3200.0, 2.0e7)).unwrap();
    let gas = resolved.state();
    gas.validate_coherence().unwrap();
    assert_eq!(
        gas.transport_state().unwrap(),
        resolved.transport_state().unwrap()
    );
    let molar_mass = 0.75 * 0.028 + 0.25 * 0.004;
    let cp_over_r = 0.75 * (3.5 + 0.001 * 3200.0) + 0.25 * 2.5;
    close(gas.specific_heat_cp, cp_over_r * R / molar_mass);
    close(gas.gamma, cp_over_r / (cp_over_r - 1.0));
    close(gas.density, 2.0e7 * molar_mass / (R * 3200.0));
    close(
        gas.sound_speed,
        (gas.gamma * R * 3200.0 / molar_mass).sqrt(),
    );
    close(gas.characteristic_impedance, gas.density * gas.sound_speed);

    // The same data-driven composition substitution updates all three views.
    let changed = model
        .with_composition(
            Composition::new(CompositionBasis::MoleFraction, vec![0.25, 0.75]).unwrap(),
        )
        .unwrap()
        .evaluate(&point(3200.0, 2.0e7))
        .unwrap();
    changed.state().validate_coherence().unwrap();
    assert!(changed.state().density < gas.density);
    assert!(changed.state().sound_speed > gas.sound_speed);
    assert!(
        changed
            .transport_state()
            .unwrap()
            .kinematic_viscosity_m2_s()
            > resolved
                .transport_state()
                .unwrap()
                .kinematic_viscosity_m2_s()
    );

    // Each original source still refuses beyond its own support; no extrapolation.
    let thermodynamic_error = model.evaluate(&point(6000.01, 100000.0)).unwrap_err();
    assert!(thermodynamic_error.to_string().contains("thermodynamics"));
    for unsupported in [point(5000.01, 100000.0), point(3200.0, 3.0e7 + 1.0)] {
        let error = model.evaluate(&unsupported).unwrap_err();
        assert!(error.to_string().contains("transport material"));
    }
    // The unrelated legacy Sutherland/USSA model keeps its declared limits.
    let legacy = crate::gas::GasSpec::dry_air_ussa1976();
    assert!(GasState::try_new(&legacy, 3200.0, 100000.0).is_err());
    assert!(GasState::try_new(&legacy, 300.0, 2.0e7).is_err());
}

#[test]
fn refuses_missing_transport_card() {
    let model = model(vec![0.75, 0.25], 0.0);
    let error = ThermochemicalGasModel::new(
        model.thermodynamics().clone(),
        vec![transport("A", 1.7e-5, 0.025)],
        None,
    )
    .unwrap_err();
    assert!(error.to_string().contains("exactly one transport material"));
}

#[test]
fn refuses_temperature_outside_nasa9_region() {
    let model = model(vec![0.75, 0.25], 0.0);
    assert!(model.evaluate(&point(200.0, 101325.0)).is_ok());
    let error = model.evaluate(&point(200.0 - 1e-9, 101325.0)).unwrap_err();
    assert!(error.to_string().contains("thermodynamics"));
}

#[test]
fn refuses_supersaturated_breath() {
    let t = 308.15;
    let p = 101325.0;
    let saturated = crate::gas::saturation_pressure_water_pa(t).unwrap() / p;
    let water = SpeciesId::new("B").unwrap();
    let make = |fraction| {
        let original = model(vec![1.0 - fraction, fraction], 0.001);
        ThermochemicalGasModel::new(
            original.thermodynamics().clone(),
            original.transport_materials().to_vec(),
            Some(&water),
        )
        .unwrap()
    };
    // The water identity and the liquid-water phase condition are caller
    // declarations, independent of this synthetic species' spelling.
    assert!(
        make(saturated)
            .evaluate_with_liquid_water_limit(&point(t, p))
            .is_ok()
    );
    let error = make(saturated * 1.01)
        .evaluate_with_liquid_water_limit(&point(t, p))
        .unwrap_err();
    assert!(error.to_string().contains("supersaturated"));
}

#[test]
fn component_cap_and_cap_plus_one() {
    use fs_thermochem::mixture::{
        MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1 as CAP, FrozenIdealGasMixtureErrorV1,
    };
    assert_eq!(CAP, 128);
    let make = |count| {
        let models = (0..count)
            .map(|i| species(&format!("component-{i:03}"), 0.028, 3.5, 0.0))
            .collect();
        FrozenIdealGasMixtureModelV1::new(
            models,
            Composition::new(
                CompositionBasis::MoleFraction,
                vec![1.0 / count as f64; count],
            )
            .unwrap(),
        )
    };
    let thermodynamics = make(CAP).unwrap();
    let transport = (0..CAP)
        .map(|i| transport(&format!("component-{i:03}"), 1.7e-5, 0.025))
        .collect();
    let model = ThermochemicalGasModel::new(thermodynamics, transport, None).unwrap();
    model
        .evaluate(&point(293.15, 101325.0))
        .unwrap()
        .state()
        .validate_coherence()
        .unwrap();
    assert!(
        matches!(make(CAP + 1), Err(FrozenIdealGasMixtureErrorV1::TooManyComponents {
        offered, limit,
    }) if offered == CAP + 1 && limit == CAP)
    );
}

#[test]
fn species_thermodynamics_and_mechanical_medium_share_one_gas_constant() {
    let model = model(vec![0.75, 0.25], 0.0);
    let evaluated = model.evaluate(&point(300.0, 100000.0)).unwrap();
    let state = evaluated.state();
    let mass = 0.75 * 0.028 + 0.25 * 0.004;
    let cp_r = 0.75 * 3.5 + 0.25 * 2.5;
    close(state.specific_gas_constant, R / mass);
    close(state.specific_heat_cp, R * cp_r / mass);
    close(state.gamma, cp_r / (cp_r - 1.0));
    close(state.density, 100000.0 * mass / (R * 300.0));
    close(state.sound_speed.powi(2), state.gamma * state.pressure / state.density);
    close(evaluated.specific_enthalpy_j_kg() - evaluated.specific_internal_energy_j_kg(),
        state.pressure / state.density);
    let transport = evaluated.transport_state().unwrap();
    close(transport.prandtl(), state.dynamic_viscosity * state.specific_heat_cp / state.thermal_conductivity);
    assert_eq!(transport.density_kg_m3(), state.density);
    assert_eq!(evaluated, model.evaluate(&point(300.0, 100000.0)).unwrap());
    assert_eq!(evaluated.components().len(), 2);
}

#[test]
fn heating_reevaluates_species_cp_and_every_transport_coefficient() {
    let model = model(vec![0.75, 0.25], 0.001);
    let cold = model.evaluate(&point(300.0, 100000.0)).unwrap();
    let hot = model.evaluate(&point(900.0, 100000.0)).unwrap();
    let (a, b) = (cold.state(), hot.state());
    assert!(b.specific_heat_cp > a.specific_heat_cp);
    assert!(b.gamma < a.gamma);
    assert!(b.dynamic_viscosity > a.dynamic_viscosity);
    assert!(b.thermal_conductivity > a.thermal_conductivity);
    assert!(b.sound_speed > a.sound_speed);
    close(b.density, a.density / 3.0);
    // A frozen-gamma sqrt(T) rescaling is demonstrably not this result.
    assert!((b.sound_speed / a.sound_speed - 3.0f64.sqrt()).abs() > 0.01);
    assert_ne!(cold.components(), hot.components());
}

#[test]
fn composition_and_pressure_change_physical_loading_not_an_audio_gain() {
    let a = model(vec![0.75, 0.25], 0.0).evaluate(&point(300.0, 100000.0)).unwrap();
    let b = model(vec![0.25, 0.75], 0.0).evaluate(&point(300.0, 100000.0)).unwrap();
    assert!(b.state().sound_speed > a.state().sound_speed);
    assert!(b.state().density < a.state().density);
    assert_ne!(b.state().dynamic_viscosity, a.state().dynamic_viscosity);
    assert_ne!(b.state().thermal_conductivity, a.state().thermal_conductivity);
    let high = model(vec![0.75, 0.25], 0.0).evaluate(&point(300.0, 200000.0)).unwrap();
    close(high.state().sound_speed, a.state().sound_speed);
    close(high.state().density, 2.0 * a.state().density);
    close(high.state().characteristic_impedance, 2.0 * a.state().characteristic_impedance);
    close(high.transport_state().unwrap().thermal_diffusivity_m2_s(),
        0.5 * a.transport_state().unwrap().thermal_diffusivity_m2_s());
}

#[test]
fn conductivity_uses_viscosity_interactions_not_conductivity_ratios() {
    let evaluated = model(vec![0.75, 0.25], 0.0).evaluate(&point(200.0, 100000.0)).unwrap();
    let x = [0.75, 0.25]; let mu = [1.7e-5_f64, 2.0e-5];
    let mass = [0.028_f64, 0.004]; let k = [0.025, 0.15];
    let mut expected = 0.0;
    for i in 0..2 {
        let mut d = 0.0;
        for j in 0..2 {
            let phi = (1.0 + (mu[i] / mu[j]).sqrt() * (mass[j] / mass[i]).powf(0.25)).powi(2)
                / (8.0 * (1.0 + mass[i] / mass[j])).sqrt();
            d += x[j] * phi;
        }
        expected += x[i] * k[i] / d;
    }
    close(evaluated.state().thermal_conductivity, expected);
}

#[test]
fn mismatched_species_invalid_coordinates_and_unsupported_sources_refuse() {
    let valid = model(vec![0.75, 0.25], 0.0);
    let thermo = valid.thermodynamics().clone();
    assert!(ThermochemicalGasModel::new(thermo.clone(), vec![transport("A", 1e-5, 0.1)], None).is_err());
    assert!(ThermochemicalGasModel::new(thermo.clone(),
        vec![transport("A", 1e-5, 0.1), transport("A", 1e-5, 0.1)], None).is_err());
    assert!(ThermochemicalGasModel::new(thermo, valid.transport_materials().to_vec(),
        Some(&SpeciesId::new("absent-water").unwrap())).is_err());
    for (t, p) in [(100.0, 100000.0), (2100.0, 100000.0), (300.0, 500.0), (300.0, 2e7)] {
        assert!(valid.evaluate(&point(t, p)).is_err());
    }
    let untyped = QueryPoint::new().with("temperature", 300.0).unwrap().with("pressure", 100000.0).unwrap();
    assert!(valid.evaluate(&untyped).is_err());
    let unstable = ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(
        vec![species("A", 0.028, 0.5, 0.0)],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap()).unwrap(),
        vec![transport("A", 1e-5, 0.1)], None).unwrap();
    assert!(unstable.evaluate(&point(300.0, 100000.0)).is_err());
}
