//! Explicit synthetic pack fixtures; no species name supplies physical data.
use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId,
    LawParameter, MaterialStateId, ModelNormalizationReceipt, ModelNormalizationTarget,
    NormalizedMaterialCardPack, NormalizedModelPack, NormalizedPack, NormalizedSpeciesPack,
    ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance, SpeciesAssociation,
    SpeciesNormalizationReceipt, SpeciesNormalizationTarget, UncertaintyModel,
    ValidityBoundSide, MODEL_PACK_TARGET_BASIS, SPECIES_MOLAR_MASS_DIMS,
    SPECIES_PACK_TARGET_BASIS, SPECIES_REFERENCE_PRESSURE_DIMS,
};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermochemical_gas::packs::SpeciesPackRef;
use fs_qty::{DynViscosity, Pressure, Temperature};
use fs_qty::chemistry::SpeciesId;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::QuantitySpec;
use fs_thermochem::{
    NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1, NASA9_STATE_SCHEMA_VERSION_V1,
    required_nasa9_parameter_dims,
};

pub(super) struct Packs {
    pub species: NormalizedSpeciesPack,
    pub nasa9: NormalizedModelPack,
    pub transport: NormalizedMaterialCardPack,
}
impl Packs {
    pub fn binding(&self) -> SpeciesPackRef<'_> {
        SpeciesPackRef {
            species: &self.species, nasa9: Some(&self.nasa9),
            transport: Some((&self.transport, MaterialPropertySelection::SingleClaimOnly)),
        }
    }
}
fn hash(text: &str) -> ContentHash {
    hash_domain("org.frankensim.tests.reed-field-input.v1", text.as_bytes())
}
fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

pub(super) fn packs(name: &str, mass: f64, cp_r: f64, minimum_temperature: f64) -> Packs {
    let artifact = hash(name);
    let provenance = Provenance {
        source: format!("synthetic {name} source; no measured gas validation"),
        license: "CC0-1.0".into(), artifact: Some(artifact),
    };
    let association = SpeciesAssociation::new(
        SpeciesId::new(name).unwrap(), mass, "gas", "ideal-gas", 100000.0,
        "synthetic-common-reference", vec![artifact], provenance.clone(),
    ).unwrap();
    let species = NormalizedSpeciesPack::new(
        name, "synthetic-field-source-v1", artifact, "CC0-1.0", association,
        vec![
            SpeciesNormalizationReceipt::new(
                SpeciesNormalizationTarget::MolarMass, hash(&format!("{mass} kg/mol")),
                SPECIES_MOLAR_MASS_DIMS, 1.0, 0.0, "kg/mol", SPECIES_PACK_TARGET_BASIS,
            ),
            SpeciesNormalizationReceipt::new(
                SpeciesNormalizationTarget::ReferencePressure, hash("100000 Pa"),
                SPECIES_REFERENCE_PRESSURE_DIMS, 1.0, 0.0, "Pa", SPECIES_PACK_TARGET_BASIS,
            ),
        ],
    ).unwrap();
    let mut coefficients = [0.0;9];
    coefficients[2] = cp_r;
    coefficients[3] = 0.001;
    let mut parameters = std::collections::BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims().into_iter().take(9).zip(coefficients) {
        parameters.insert(name.to_owned(), LawParameter { value, dims });
    }
    parameters.insert("reference_pressure".into(), LawParameter { value: 100000.0, dims: Pressure::DIMS });
    let model = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()), law_version: NASA9_LAW_VERSION_V1,
        parameters, state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", minimum_temperature, 1000.0),
        sources: vec![artifact], provenance: provenance.clone(),
    };
    let identity = model.content_hash();
    let mut normalizations = Vec::new();
    for (name, parameter) in &model.parameters {
        normalizations.push(ModelNormalizationReceipt::new(
            ModelNormalizationTarget::Parameter { model: identity, parameter: name.clone() },
            hash(&format!("{name} {}", parameter.value)), parameter.dims, 1.0, 0.0,
            "authored SI", MODEL_PACK_TARGET_BASIS, None, None,
        ));
    }
    for (side, value) in [(ValidityBoundSide::Lower, minimum_temperature), (ValidityBoundSide::Upper, 1000.0)] {
        normalizations.push(ModelNormalizationReceipt::new(
            ModelNormalizationTarget::ValidityBound { model: identity, axis: "T".into(), side },
            hash(&format!("{value} K")), Temperature::DIMS, 1.0, 0.0, "K",
            MODEL_PACK_TARGET_BASIS, None, None,
        ));
    }
    let nasa9 = NormalizedModelPack::new(
        name, "synthetic-field-model-v1", artifact, "CC0-1.0", vec![model], normalizations,
    ).unwrap();
    let mut claims = ClaimSet::new();
    let observation = claims.register_observation(ObservationDataset {
        specimen: format!("synthetic {name} transport"), method: "manufactured physical regression".into(),
        artifact, caveats: "not a measured gas".into(), provenance: provenance.clone(),
    }).unwrap();
    for (name, dims, value) in [
        ("dynamic_viscosity", DynViscosity::DIMS, 1.8e-5),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 0.026),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained()
                .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), minimum_temperature, 1000.0)
                .with_quantity("pressure", typed(QuantityKind::Pressure), 1000.0, 1e7),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation.clone()], provenance: provenance.clone(),
        }).unwrap();
    }
    let claims = NormalizedPack::new(
        name, "synthetic-field-transport-v1", artifact, "CC0-1.0", claims, vec![], vec![],
    ).unwrap();
    let transport = NormalizedMaterialCardPack::new(MaterialStateId {
        chemistry: name.into(), phase: "gas".into(), process: "explicit synthetic association".into(),
        revision: 0,
    }, claims).unwrap();
    Packs { species, nasa9, transport }
}
