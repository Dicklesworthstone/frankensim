//! Synthetic public-API NASA/transport fixtures with original methane PR data.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_material::thermochemical_gas::nonideal::equilibrium::PureFluidEquilibriumControls;
use fs_material::thermochemical_gas::nonideal::{
    PengRobinsonSpeciesBinding, ThermochemicalPengRobinsonModel,
};
use fs_material::thermochemical_gas::{SpeciesTransportMaterial, ThermochemicalGasModel};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{DynViscosity, MolarMass, Pressure, QuantitySpec, Temperature};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1,
    NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1, StandardStatePhaseV1,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, required_nasa9_parameter_dims,
};

pub const MOLAR_MASS: f64 = 0.016_042_46;
fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn provenance() -> Provenance {
    Provenance {
        source: "synthetic equilibrium caloric fixture".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}
pub fn model(datum_j_kg: f64, domain: [f64; 2]) -> ThermochemicalPengRobinsonModel {
    let mut coefficients = [0.0; 9];
    coefficients[2] = 3.5;
    coefficients[7] = MOLAR_MASS * datum_j_kg / R;
    let mut parameters = std::collections::BTreeMap::new();
    for ((key, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(key.to_owned(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100_000.0,
            dims: Pressure::DIMS,
        },
    );
    let nasa = Nasa9StandardStateModelV1::new(
        SpeciesId::new("methane").unwrap(),
        MolarMass::new(MOLAR_MASS),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100_000.0),
            ElementalReferenceIdV1::new("synthetic-reference").unwrap(),
        )
        .unwrap(),
        vec![
            Nasa9RegionV1::from_card(ConstitutiveModelCard {
                law: LawId(NASA9_LAW_ID_V1.into()),
                law_version: NASA9_LAW_VERSION_V1,
                parameters,
                state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
                initial_state: InitialStatePolicy::ZeroInternalState,
                validity: ValidityDomain::unconstrained().with("T", domain[0], domain[1]),
                sources: vec![],
                provenance: provenance(),
            })
            .unwrap(),
        ],
    )
    .unwrap();
    let mut claims = ClaimSet::new();
    for (property, dims, value) in [
        ("dynamic_viscosity", DynViscosity::DIMS, 1e-5),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            0.03,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(property, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "temperature",
                        typed(QuantityKind::AbsoluteTemperature),
                        200.0,
                        2000.0,
                    )
                    .with_quantity("pressure", typed(QuantityKind::Pressure), 1000.0, 1e7),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
    }
    let source = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![nasa],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        vec![SpeciesTransportMaterial {
            species: SpeciesId::new("methane").unwrap(),
            material: MaterialCard::assemble(
                MaterialStateId {
                    chemistry: "methane".into(),
                    phase: "gas".into(),
                    process: "synthetic".into(),
                    revision: 0,
                },
                claims,
                vec![],
            )
            .unwrap(),
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        None,
    )
    .unwrap();
    ThermochemicalPengRobinsonModel::new(
        source,
        vec![PengRobinsonSpeciesBinding {
            species: SpeciesId::new("methane").unwrap(),
            critical_temperature: Temperature::new(190.564),
            critical_pressure: Pressure::new(4.5992e6),
            acentric_factor: 0.01142,
        }],
        vec![],
    )
    .unwrap()
}

pub fn controls() -> PureFluidEquilibriumControls {
    PureFluidEquilibriumControls::new(
        [Temperature::new(130.0), Temperature::new(180.0)],
        [Pressure::new(1e5), Pressure::new(4.4e6)],
    )
}
