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

pub(super) fn species(name: &str, mass: f64, cp_r: f64, slope: f64) -> Nasa9StandardStateModelV1 {
    let mut coefficients = [0.0; 9];
    coefficients[2] = cp_r;
    coefficients[3] = slope;
    let mut parameters = std::collections::BTreeMap::new();
    for ((key, dims), value) in required_nasa9_parameter_dims().into_iter().take(9).zip(coefficients) {
        parameters.insert(key.to_owned(), LawParameter { value, dims });
    }
    parameters.insert("reference_pressure".into(), LawParameter { value: 100000.0, dims: Pressure::DIMS });
    let card = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()), law_version: NASA9_LAW_VERSION_V1,
        parameters, state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 2000.0),
        sources: vec![], provenance: Provenance { source: "synthetic NASA-9 gas".into(),
            license: "CC0-1.0".into(), artifact: None },
    };
    Nasa9StandardStateModelV1::new(SpeciesId::new(name).unwrap(), MolarMass::new(mass),
        StandardStateConventionV1::new(StandardStatePhaseV1::Gas, ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0), ElementalReferenceIdV1::new("synthetic-reference").unwrap()).unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()]).unwrap()
}

pub(super) fn transport(name: &str, mu: f64, k: f64) -> SpeciesTransportMaterial {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [("dynamic_viscosity", DynViscosity::DIMS, mu),
        (THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, k)] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve { abscissa: "temperature".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(200.0, value), (2000.0, 2.0 * value)], dims },
            validity: ValidityDomain::unconstrained()
                .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), 200.0, 2000.0)
                .with_quantity("pressure", typed(QuantityKind::Pressure), 1000.0, 1.0e7),
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

pub(super) fn model(fractions: Vec<f64>, slope: f64) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(
        vec![species("A", 0.028, 3.5, slope), species("B", 0.004, 2.5, 0.0)],
        Composition::new(CompositionBasis::MoleFraction, fractions).unwrap()).unwrap(),
        vec![transport("B", 2.0e-5, 0.15), transport("A", 1.7e-5, 0.025)], None).unwrap()
}
pub(super) fn point(t: f64, p: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), t).unwrap()
        .with_quantity("pressure", typed(QuantityKind::Pressure), p).unwrap()
}
pub(super) fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 2e-12 * a.abs().max(b.abs()).max(1e-100), "{a} != {b}");
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
