//! Authored synthetic NASA9/transport packs; no empirical-gas claim.
use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::ValidityDomain;
use fs_matdb::*;
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermochemical_gas::{
    ThermochemicalGasModel,
    packs::{SpeciesPackRef, from_packs},
};
use fs_qty::chemistry::SpeciesId;
use fs_qty::semantic::{Composition, CompositionBasis, QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, DynViscosity, Pressure, QuantitySpec, Temperature};

fn hash(text: &str) -> ContentHash {
    hash_domain(
        "org.frankensim.tests.nsf-thermochemical.v1",
        text.as_bytes(),
    )
}
fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
pub fn point(t: f64, p: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), t)
        .unwrap()
        .with_quantity("pressure", typed(QuantityKind::Pressure), p)
        .unwrap()
}

/// cp/R = 3 + 0.001 T; a7 is an arbitrary formation-energy offset [K].
/// Both transport properties vary with pressure; k_scale changes data only.
pub fn model(a7: f64, pressure_max: f64, k_scale: f64) -> ThermochemicalGasModel {
    model_with_temperature_limit(a7, pressure_max, k_scale, 700.0)
}

pub fn model_with_temperature_limit(
    a7: f64,
    pressure_max: f64,
    k_scale: f64,
    temperature_max: f64,
) -> ThermochemicalGasModel {
    let (species, nasa9, transport) = source_packs(
        "manufactured",
        0.028,
        a7,
        pressure_max,
        k_scale,
        temperature_max,
    );
    from_packs(
        &[SpeciesPackRef {
            species: &species,
            nasa9: Some(&nasa9),
            transport: Some((&transport, MaterialPropertySelection::SingleClaimOnly)),
        }],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        None,
    )
    .unwrap()
}

/// The same authored source construction with explicit water and carrier IDs.
/// Water's thermal coefficients remain a synthetic coupling fixture.
#[allow(dead_code)]
pub fn moist_model(water_fraction: f64, a7: f64) -> ThermochemicalGasModel {
    let carrier = source_packs("manufactured", 0.028, a7, 250000.0, 1.0, 700.0);
    let water = source_packs("water", 0.01801528, a7, 250000.0, 1.0, 700.0);
    from_packs(
        &[
            SpeciesPackRef {
                species: &carrier.0,
                nasa9: Some(&carrier.1),
                transport: Some((&carrier.2, MaterialPropertySelection::SingleClaimOnly)),
            },
            SpeciesPackRef {
                species: &water.0,
                nasa9: Some(&water.1),
                transport: Some((&water.2, MaterialPropertySelection::SingleClaimOnly)),
            },
        ],
        Composition::new(
            CompositionBasis::MoleFraction,
            vec![1.0 - water_fraction, water_fraction],
        )
        .unwrap(),
        Some(&SpeciesId::new("water").unwrap()),
    )
    .unwrap()
}

fn source_packs(
    name: &str,
    molar_mass: f64,
    a7: f64,
    pressure_max: f64,
    k_scale: f64,
    temperature_max: f64,
) -> (
    NormalizedSpeciesPack,
    NormalizedModelPack,
    NormalizedMaterialCardPack,
) {
    let artifact = hash("synthetic gas");
    let provenance = Provenance {
        source: "authored synthetic NSF fixture; not measured".into(),
        license: "CC0-1.0".into(),
        artifact: Some(artifact),
    };
    let species_id = SpeciesId::new(name).unwrap();
    let association = SpeciesAssociation::new(
        species_id,
        molar_mass,
        "gas",
        "ideal-gas",
        100000.0,
        "synthetic-reference",
        vec![artifact],
        provenance.clone(),
    )
    .unwrap();
    let species = NormalizedSpeciesPack::new(
        name,
        "nsf-v1",
        artifact,
        "CC0-1.0",
        association,
        vec![
            SpeciesNormalizationReceipt::new(
                SpeciesNormalizationTarget::MolarMass,
                hash("mass"),
                SPECIES_MOLAR_MASS_DIMS,
                1.0,
                0.0,
                "kg/mol",
                SPECIES_PACK_TARGET_BASIS,
            ),
            SpeciesNormalizationReceipt::new(
                SpeciesNormalizationTarget::ReferencePressure,
                hash("p0"),
                SPECIES_REFERENCE_PRESSURE_DIMS,
                1.0,
                0.0,
                "Pa",
                SPECIES_PACK_TARGET_BASIS,
            ),
        ],
    )
    .unwrap();
    let mut parameters = std::collections::BTreeMap::new();
    for (i, power) in [2, 1, 0, -1, -2, -3, -4, 1, 0].into_iter().enumerate() {
        parameters.insert(
            format!("a{i}"),
            LawParameter {
                value: match i {
                    2 => 3.0,
                    3 => 0.001,
                    7 => a7,
                    _ => 0.0,
                },
                dims: Dims([0, 0, 0, power, 0, 0]),
            },
        );
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100000.0,
            dims: Pressure::DIMS,
        },
    );
    let law = ConstitutiveModelCard {
        law: LawId("nasa9-standard-state".into()),
        law_version: 1,
        parameters,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, temperature_max),
        sources: vec![artifact],
        provenance: provenance.clone(),
    };
    let identity = law.content_hash();
    let mut normalized = law
        .parameters
        .iter()
        .map(|(name, p)| {
            ModelNormalizationReceipt::new(
                ModelNormalizationTarget::Parameter {
                    model: identity,
                    parameter: name.clone(),
                },
                hash(&format!("{name} {}", p.value)),
                p.dims,
                1.0,
                0.0,
                "authored SI",
                MODEL_PACK_TARGET_BASIS,
                None,
                None,
            )
        })
        .collect::<Vec<_>>();
    for (side, value) in [
        (ValidityBoundSide::Lower, 200.0),
        (ValidityBoundSide::Upper, temperature_max),
    ] {
        normalized.push(ModelNormalizationReceipt::new(
            ModelNormalizationTarget::ValidityBound {
                model: identity,
                axis: "T".into(),
                side,
            },
            hash(&format!("{value} K")),
            Temperature::DIMS,
            1.0,
            0.0,
            "K",
            MODEL_PACK_TARGET_BASIS,
            None,
            None,
        ));
    }
    let nasa9 =
        NormalizedModelPack::new(name, "nsf-v1", artifact, "CC0-1.0", vec![law], normalized)
            .unwrap();
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "synthetic gas".into(),
            method: "authored algebra".into(),
            artifact,
            caveats: "not empirical".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    for (name, dims, base) in [
        ("dynamic_viscosity", DynViscosity::DIMS, 1.8e-5),
        (
            "thermal_conductivity",
            THERMAL_CONDUCTIVITY_DIMS,
            0.026 * k_scale,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "pressure".into(),
                    abscissa_dims: Pressure::DIMS,
                    knots: vec![(1000.0, 0.5 * base), (250000.0, 1.5 * base)],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "temperature",
                        typed(QuantityKind::AbsoluteTemperature),
                        200.0,
                        temperature_max,
                    )
                    .with_quantity(
                        "pressure",
                        typed(QuantityKind::Pressure),
                        1000.0,
                        pressure_max,
                    ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![observation.clone()],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    let pack = NormalizedPack::new(
        "synthetic",
        "nsf-v1",
        artifact,
        "CC0-1.0",
        claims,
        vec![],
        vec![],
    )
    .unwrap();
    let transport = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: name.into(),
            phase: "gas".into(),
            process: "authored".into(),
            revision: 0,
        },
        pack,
    )
    .unwrap();
    (species, nasa9, transport)
}
