//! Shared analytic/source fixtures for material-to-field ionization tests.
#![allow(dead_code)]
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::ionized_gas::{
    CarrierMomentumRelaxation, IonizedGasConductivity, IonizedGasDcLimits,
    MOMENTUM_RELAXATION_RATE_DIMS,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Amount, MolarMass, Pressure, QuantitySpec, Temperature};
use fs_thermochem::equilibrium::composition::{
    EquilibriumSolveV1, IdealGasEquilibriumSolutionV1, IdealGasEquilibriumV1,
};
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, NASA9_LAW_ID_V1,
    NASA9_LAW_VERSION_V1, NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1, StandardStatePhaseV1,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K, required_nasa9_parameter_dims,
};
use std::collections::BTreeMap;

pub const NA: f64 = 6.022_140_76e23;
pub const KB: f64 = 1.380_649e-23;
pub const E: f64 = 1.602_176_634e-19;
pub const ME: f64 = 9.109_383_713_9e-31;
pub const MH: f64 = 1.673_557_5e-27;

pub fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn standard_state(name: &str, mass: f64, weight: f64, formation: f64) -> Nasa9StandardStateModelV1 {
    // Sackur-Tetrode ideal translational partition function, with explicit
    // ground-state degeneracy and ionization reference energy; not a fitted gas.
    let h: f64 = 6.626_070_15e-34;
    let entropy_constant = 2.5
        + (weight * (2.0 * core::f64::consts::PI * mass * KB / (h * h)).powf(1.5) * KB / 100_000.0)
            .ln();
    let c = [
        0.0,
        0.0,
        2.5,
        0.0,
        0.0,
        0.0,
        0.0,
        formation / UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
        entropy_constant,
    ];
    let mut parameters = BTreeMap::new();
    for ((key, dims), value) in required_nasa9_parameter_dims().into_iter().take(9).zip(c) {
        parameters.insert(key.to_owned(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100_000.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()), law_version: NASA9_LAW_VERSION_V1,
        parameters, state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 2000.0, 20000.0),
        sources: vec![], provenance: Provenance {
            source: "analytic monatomic H/H+/e Sackur-Tetrode fixture; NIST H ionization 13.598434599702 eV".into(),
            license: "analytic test fixture".into(), artifact: None },
    };
    Nasa9StandardStateModelV1::new(
        id(name),
        MolarMass::new(mass * NA),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100_000.0),
            ElementalReferenceIdV1::new("monatomic-H-neutral-ground-reference").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}
pub fn chemistry() -> IdealGasEquilibriumV1 {
    let models = vec![
        standard_state("H", MH, 2.0, 0.0),
        standard_state("H+", MH - ME, 1.0, 13.598_434_599_702 * E * NA),
        standard_state("e-", ME, 2.0, 0.0),
    ];
    let species: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    let matrix = ElementalMatrix::new(
        vec![ElementId::new("H").unwrap()],
        species.clone(),
        vec![vec![1, 1, 0]],
    )
    .unwrap();
    let charges = ChargeVector::new(species, vec![0, 1, -1]).unwrap();
    IdealGasEquilibriumV1::new(matrix, charges, models).unwrap()
}
pub fn equilibrium(owner: &IdealGasEquilibriumV1, t: f64, p: f64) -> IdealGasEquilibriumSolutionV1 {
    owner
        .equilibrate_tp(
            Temperature::new(t),
            Pressure::new(p),
            &[(id("H"), Amount::new(1.0))],
            EquilibriumSolveV1::default(),
            || false,
        )
        .unwrap()
}
pub fn source(species: &str, t: f64, p: f64, rate: f64) -> CarrierMomentumRelaxation {
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let pressure =
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
    let key = PropertyKey::new(
        "effective_momentum_relaxation_rate",
        MOMENTUM_RELAXATION_RATE_DIMS,
    );
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Scalar {
                value: rate,
                dims: MOMENTUM_RELAXATION_RATE_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("T", temperature, 2000.0, 20000.0)
                .with_quantity("p", pressure, 1.0, 1e8),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: format!("synthetic {species}/fixed-neutral-background momentum relaxation"),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: species.into(),
            phase: "gas".into(),
            process: "synthetic collision fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let query = QueryPoint::new()
        .with_quantity("T", temperature, t)
        .unwrap()
        .with_quantity("p", pressure, p)
        .unwrap();
    let resolved = resolve_material_state_point(
        &card,
        &query,
        &[
            ScalarPropertyRequirement::try_with_key(&key, ScalarAdmissibility::StrictlyPositive)
                .unwrap(),
        ],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    CarrierMomentumRelaxation::new(id(species), resolved, key).unwrap()
}
pub fn limits() -> IonizedGasDcLimits {
    IonizedGasDcLimits {
        maximum_charged_mole_fraction: 0.2,
        maximum_coulomb_coupling: 0.1,
        maximum_thermal_phase_space_density: 0.1,
        maximum_drift_thermal_ratio: 0.01,
    }
}
pub fn rates(t: f64, p: f64, multiplier: f64) -> Vec<CarrierMomentumRelaxation> {
    vec![
        source("e-", t, p, multiplier * 1e10),
        source("H+", t, p, multiplier * 1e9),
    ]
}
pub fn material(t: f64, p: f64) -> IonizedGasConductivity {
    let owner = chemistry();
    let state = equilibrium(&owner, t, p);
    IonizedGasConductivity::new(&owner, state, rates(t, p, 1.0), limits(), || false).unwrap()
}
pub fn relative(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance * b.abs(), "{a:e} != {b:e}");
}
