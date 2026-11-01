//! Synthetic binary calorics with original methane/ethane PR constants.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_material::thermochemical_gas::nonideal::flash::enthalpy::FlashEnthalpyControls;
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

pub const MIXTURE_MASS: f64 = (0.016_042_46 + 0.030_069_04) / 2.0;
pub const METHANE: usize = 1; // Canonical species order is ethane, methane.

fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn provenance() -> Provenance {
    Provenance {
        source: "synthetic binary flash calorics".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}

/// Input datum arrays use methane, ethane order, independently of canonical order.
pub fn model(
    enthalpy_datums_j_kg: [f64; 2],
    entropy_datums_j_kg_k: [f64; 2],
    domain: [f64; 2],
) -> ThermochemicalPengRobinsonModel {
    let mut models = Vec::new();
    let mut transport = Vec::new();
    let mut bindings = Vec::new();
    for (index, (name, mass, cp_r, tc, pc, omega)) in [
        ("methane", 0.016_042_46, 3.5, 190.564, 4.5992e6, 0.01142),
        ("ethane", 0.030_069_04, 4.5, 305.322, 4.8722e6, 0.0995),
    ]
    .into_iter()
    .enumerate()
    {
        let species = SpeciesId::new(name).unwrap();
        let mut coefficients = [0.0; 9];
        coefficients[2] = cp_r;
        coefficients[7] = mass * enthalpy_datums_j_kg[index] / R;
        coefficients[8] = mass * entropy_datums_j_kg_k[index] / R;
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
                value: 1e5,
                dims: Pressure::DIMS,
            },
        );
        models.push(
            Nasa9StandardStateModelV1::new(
                species.clone(),
                MolarMass::new(mass),
                StandardStateConventionV1::new(
                    StandardStatePhaseV1::Gas,
                    ReferenceEquationOfStateV1::IdealGas,
                    Pressure::new(1e5),
                    ElementalReferenceIdV1::new("synthetic-binary-reference").unwrap(),
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
            .unwrap(),
        );
        let mut claims = ClaimSet::new();
        for (name, dims, value) in [
            ("dynamic_viscosity", DynViscosity::DIMS, 1e-5),
            (
                THERMAL_CONDUCTIVITY_PROPERTY,
                THERMAL_CONDUCTIVITY_DIMS,
                0.03,
            ),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: PropertyValue::Scalar { value, dims },
                    // Deliberately excludes all cryogenic flash states. Calorics
                    // must not require or extrapolate dilute gas transport cards.
                    validity: ValidityDomain::unconstrained()
                        .with_quantity(
                            "temperature",
                            typed(QuantityKind::AbsoluteTemperature),
                            300.0,
                            2000.0,
                        )
                        .with_quantity("pressure", typed(QuantityKind::Pressure), 1e3, 1e7),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    observations: vec![],
                    provenance: provenance(),
                })
                .unwrap();
        }
        transport.push(SpeciesTransportMaterial {
            species: species.clone(),
            material: MaterialCard::assemble(
                MaterialStateId {
                    chemistry: name.into(),
                    phase: "gas".into(),
                    process: "synthetic".into(),
                    revision: 0,
                },
                claims,
                vec![],
            )
            .unwrap(),
            selection: MaterialPropertySelection::SingleClaimOnly,
        });
        bindings.push(PengRobinsonSpeciesBinding {
            species,
            critical_temperature: Temperature::new(tc),
            critical_pressure: Pressure::new(pc),
            acentric_factor: omega,
        });
    }
    let source = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            models,
            Composition::new(CompositionBasis::MoleFraction, vec![0.5, 0.5]).unwrap(),
        )
        .unwrap(),
        transport,
        None,
    )
    .unwrap();
    ThermochemicalPengRobinsonModel::new(source, bindings, vec![]).unwrap()
}

pub fn source() -> ThermochemicalPengRobinsonModel {
    model([0.0; 2], [0.0; 2], [100.0, 450.0])
}
pub fn controls() -> FlashEnthalpyControls {
    FlashEnthalpyControls::new([Temperature::new(210.0), Temperature::new(230.0)])
}
