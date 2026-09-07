//! G0/G3 conformance for the normalized material-card pack boundary.

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MATERIAL_CARD_MODEL_PACK_SCHEMA_VERSION, MATERIAL_CARD_PACK_SCHEMA_VERSION,
    MODEL_PACK_TARGET_BASIS, MaterialStateId, ModelNormalizationReceipt, ModelNormalizationTarget,
    NormalizedMaterialCardPack, NormalizedModelPack, NormalizedPack, ObservationDataset, PackError,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
    UncertaintyModel, ValidityBoundSide,
};
use fs_qty::Dims;

const SOURCE_DOMAIN: &str = "org.frankensim.tests.material-card-pack.source.v1";
const THERMAL_CONDUCTIVITY_DIMS: Dims = Dims([1, 1, -3, -1, 0, 0]);

fn provenance() -> Provenance {
    Provenance {
        source: "guarded-hot-plate campaign GHP-7".to_string(),
        license: "CC-BY-4.0; redistribution permitted with attribution".to_string(),
        artifact: Some(hash_domain(SOURCE_DOMAIN, b"fixture-table")),
    }
}

fn aluminum_state() -> MaterialStateId {
    MaterialStateId {
        chemistry: "AA6061".to_string(),
        phase: "wrought".to_string(),
        process: "T6".to_string(),
        revision: 0,
    }
}

fn claims_pack() -> NormalizedPack {
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "AA6061-T6 guarded-hot-plate coupon".to_string(),
            method: "GHP-7 steady conduction campaign".to_string(),
            artifact: hash_domain(SOURCE_DOMAIN, b"raw-observation"),
            caveats: "fixture value; not a seed-dataset authority".to_string(),
            provenance: provenance(),
        })
        .expect("licensed observation inserts");
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("thermal-conductivity", THERMAL_CONDUCTIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: 167.0,
                dims: THERMAL_CONDUCTIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 273.15, 373.15),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation],
            provenance: provenance(),
        })
        .expect("conductivity claim inserts");
    NormalizedPack::new(
        "fixture-aa6061-t6-thermal",
        "frankensim-material-card-pack-compiler-v1",
        hash_domain(SOURCE_DOMAIN, b"source-envelope"),
        "CC-BY-4.0: redistribution permitted with attribution",
        claims,
        Vec::new(),
        Vec::new(),
    )
    .expect("claim pack admits")
}

fn sample_pack() -> NormalizedMaterialCardPack {
    NormalizedMaterialCardPack::new(aluminum_state(), claims_pack())
        .expect("material-card pack admits")
}

fn model_pack(conductivities: &[f64], compiler: &str) -> NormalizedModelPack {
    let mut models = Vec::new();
    let mut receipts = Vec::new();
    for &value in conductivities {
        let model = ConstitutiveModelCard {
            law: LawId("synthetic-fourier".into()),
            law_version: 3,
            parameters: std::collections::BTreeMap::from([(
                "conductivity".into(),
                LawParameter {
                    value,
                    dims: THERMAL_CONDUCTIVITY_DIMS,
                },
            )]),
            state_schema_version: 0,
            initial_state: InitialStatePolicy::ZeroInternalState,
            validity: ValidityDomain::unconstrained().with("T", 273.15, 373.15),
            sources: vec![hash_domain(SOURCE_DOMAIN, b"synthetic-model-source")],
            provenance: Provenance {
                source: "synthetic transport fixture; no measured calibration".into(),
                ..provenance()
            },
        };
        let hash = model.content_hash();
        for (target, dims, literal) in [
            (
                ModelNormalizationTarget::Parameter {
                    model: hash,
                    parameter: "conductivity".into(),
                },
                THERMAL_CONDUCTIVITY_DIMS,
                format!("{value} W/(m K)"),
            ),
            (
                ModelNormalizationTarget::ValidityBound {
                    model: hash,
                    axis: "T".into(),
                    side: ValidityBoundSide::Lower,
                },
                Dims([0, 0, 0, 1, 0, 0]),
                "273.15 K".into(),
            ),
            (
                ModelNormalizationTarget::ValidityBound {
                    model: hash,
                    axis: "T".into(),
                    side: ValidityBoundSide::Upper,
                },
                Dims([0, 0, 0, 1, 0, 0]),
                "373.15 K".into(),
            ),
        ] {
            receipts.push(ModelNormalizationReceipt::new(
                target,
                hash_domain(SOURCE_DOMAIN, literal.as_bytes()),
                dims,
                1.0,
                0.0,
                "authored SI",
                MODEL_PACK_TARGET_BASIS,
                None,
                None,
            ));
        }
        models.push(model);
    }
    NormalizedModelPack::new(
        "explicitly-associated-models",
        compiler,
        hash_domain(SOURCE_DOMAIN, b"model-envelope"),
        "synthetic fixture redistribution",
        models,
        receipts,
    )
    .unwrap()
}

#[test]
fn g0_model_association_preserves_membership_and_complete_nested_evidence() {
    let models = model_pack(&[167.0, 170.0], "test-v1");
    let pack = NormalizedMaterialCardPack::new_with_models(
        aluminum_state(),
        claims_pack(),
        models.clone(),
    )
    .unwrap();
    assert_eq!(
        pack.schema_version(),
        MATERIAL_CARD_MODEL_PACK_SCHEMA_VERSION
    );
    let bytes = pack.to_bytes();
    let decoded =
        NormalizedMaterialCardPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap();
    assert_eq!(decoded, pack);
    assert_eq!(decoded.model_pack(), Some(&models));
    assert_eq!(decoded.card().models(), models.models());
    assert_eq!(
        decoded
            .card()
            .models_for(&LawId("synthetic-fourier".into()))
            .len(),
        2
    );
    assert_eq!(decoded.model_pack().unwrap().normalizations().len(), 6);
    assert_eq!(decoded.claims_pack(), &claims_pack());
    let reversed = NormalizedMaterialCardPack::new_with_models(
        aluminum_state(),
        claims_pack(),
        model_pack(&[170.0, 167.0], "test-v1"),
    )
    .unwrap();
    assert_eq!(reversed.to_bytes(), bytes);

    let changed_source = NormalizedMaterialCardPack::new_with_models(
        aluminum_state(),
        claims_pack(),
        model_pack(&[167.0, 170.0], "test-v2"),
    )
    .unwrap();
    assert_eq!(
        changed_source.card().content_hash(),
        pack.card().content_hash()
    );
    assert_ne!(changed_source.content_hash(), pack.content_hash());
    let changed_law = NormalizedMaterialCardPack::new_with_models(
        aluminum_state(),
        claims_pack(),
        model_pack(&[168.0, 170.0], "test-v1"),
    )
    .unwrap();
    assert_ne!(
        changed_law.card().content_hash(),
        pack.card().content_hash()
    );
}

#[test]
fn g0_model_association_refuses_substitution_and_silent_downgrade() {
    let plain = sample_pack();
    assert_eq!(plain.schema_version(), 1);
    assert_eq!(plain.model_pack(), None);
    assert_eq!(
        plain.content_hash(),
        hash_domain(
            "org.frankensim.fs-matdb.normalized-material-card-pack.v1",
            &plain.to_bytes()
        )
    );
    let models = model_pack(&[167.0], "test-v1");
    let pack = NormalizedMaterialCardPack::new_with_models(
        aluminum_state(),
        claims_pack(),
        models.clone(),
    )
    .unwrap();
    let bytes = pack.to_bytes();
    let model_offset = plain.to_bytes().len();
    assert_eq!(
        &bytes[model_offset..model_offset + 32],
        &models.content_hash().0
    );
    assert_eq!(
        pack.content_hash(),
        hash_domain(
            "org.frankensim.fs-matdb.normalized-material-card-pack.v2",
            &bytes
        )
    );
    let mut changed = bytes.clone();
    changed[model_offset + 36] ^= 1;
    assert!(matches!(
        NormalizedMaterialCardPack::from_bytes(&changed),
        Err(PackError::IdentityMismatch {
            kind: "model pack",
            ..
        })
    ));
    assert!(matches!(
        NormalizedMaterialCardPack::from_bytes_verified(pack.content_hash(), &changed),
        Err(PackError::IdentityMismatch {
            kind: "material_card_pack",
            ..
        })
    ));

    // Removing the nested models and relabeling as v1 cannot retain the card
    // identity that declared those members, even with a recomputed outer hash.
    let mut downgraded = bytes[..model_offset].to_vec();
    downgraded[8..12].copy_from_slice(&1_u32.to_le_bytes());
    assert!(matches!(
        NormalizedMaterialCardPack::from_bytes(&downgraded),
        Err(PackError::IdentityMismatch {
            kind: "material_card",
            ..
        })
    ));
    let mut upgraded = plain.to_bytes();
    upgraded[8..12].copy_from_slice(&2_u32.to_le_bytes());
    assert!(matches!(
        NormalizedMaterialCardPack::from_bytes(&upgraded),
        Err(PackError::Malformed { .. })
    ));
}

#[test]
fn material_card_pack_round_trips_deterministically() {
    let pack = sample_pack();
    let first = pack.to_bytes();
    let second = sample_pack().to_bytes();
    // Independent frozen v1 grammar: no optional-model marker or extra field
    // may enter an existing model-free artifact.
    let claims = claims_pack();
    let card =
        fs_matdb::MaterialCard::assemble(aluminum_state(), claims.claims().clone(), Vec::new())
            .unwrap();
    let mut legacy = b"FSMCDPK\0".to_vec();
    legacy.extend_from_slice(&1_u32.to_le_bytes());
    for text in ["AA6061", "wrought", "T6"] {
        legacy.extend_from_slice(&(text.len() as u32).to_le_bytes());
        legacy.extend_from_slice(text.as_bytes());
    }
    legacy.extend_from_slice(&0_u32.to_le_bytes());
    legacy.extend_from_slice(&card.content_hash().0);
    legacy.extend_from_slice(&claims.content_hash().0);
    let payload = claims.to_bytes();
    legacy.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    legacy.extend_from_slice(&payload);
    assert_eq!(first, legacy, "model-free v1 wire layout is frozen");
    assert_eq!(first, second, "canonical material-card bytes moved");
    assert_eq!(&first[..8], b"FSMCDPK\0");
    assert_eq!(
        u32::from_le_bytes(first[8..12].try_into().expect("version width")),
        MATERIAL_CARD_PACK_SCHEMA_VERSION
    );

    let decoded = NormalizedMaterialCardPack::from_bytes(&first).expect("pack decodes");
    assert_eq!(decoded, pack);
    assert_eq!(decoded.pack_id(), "fixture-aa6061-t6-thermal");
    assert_eq!(
        decoded.compiler(),
        "frankensim-material-card-pack-compiler-v1"
    );
    assert_eq!(decoded.card().id(), &aluminum_state());
    assert_eq!(decoded.card().claims_for("thermal-conductivity").len(), 1);
    assert!(decoded.card().models().is_empty(), "v1 carries no models");
    assert_eq!(
        decoded.claims_pack().content_hash(),
        pack.claims_pack().content_hash()
    );
    assert_eq!(decoded.card().content_hash(), pack.card().content_hash());
    assert_eq!(decoded.content_hash(), pack.content_hash());
    assert_eq!(
        NormalizedMaterialCardPack::from_bytes_verified(pack.content_hash(), &first)
            .expect("whole pack identity verifies"),
        pack
    );
}

#[test]
fn reconstructed_card_answers_a_binding_style_query_with_a_receipt() {
    // The card must be consumable exactly the way fs-project's binding
    // resolution consumes it: a policy-driven claim-set query at a typed
    // point, leaving a replayable usage receipt.
    let pack = sample_pack();
    let decoded = NormalizedMaterialCardPack::from_bytes(&pack.to_bytes()).expect("pack decodes");
    let point = QueryPoint::new().with("T", 300.0).expect("finite point");
    let answer = decoded
        .card()
        .claims()
        .query(
            "thermal-conductivity",
            &point,
            SelectionPolicy::SingleClaimOnly,
        )
        .expect("in-domain conductivity query resolves");
    assert_eq!(answer.evidence.value.value.to_bits(), 167.0_f64.to_bits());
    decoded
        .card()
        .claims()
        .verify_receipt(&answer.receipt)
        .expect("usage receipt must verify against the reconstructed claim set");
}

#[test]
fn declared_state_moves_the_pack_identity() {
    let baseline = sample_pack();
    let annealed = NormalizedMaterialCardPack::new(
        MaterialStateId {
            process: "O-annealed".to_string(),
            ..aluminum_state()
        },
        claims_pack(),
    )
    .expect("annealed state admits");

    assert_ne!(
        baseline.card().content_hash(),
        annealed.card().content_hash()
    );
    assert_ne!(baseline.content_hash(), annealed.content_hash());
}

#[test]
fn malformed_or_unpinned_material_artifacts_refuse() {
    let pack = sample_pack();
    let bytes = pack.to_bytes();

    let mut bad_magic = bytes.clone();
    bad_magic[0] ^= 0xff;
    assert!(matches!(
        NormalizedMaterialCardPack::from_bytes(&bad_magic),
        Err(PackError::Malformed { .. })
    ));

    let mut bad_version = bytes.clone();
    bad_version[8..12].copy_from_slice(&99_u32.to_le_bytes());
    assert!(matches!(
        NormalizedMaterialCardPack::from_bytes(&bad_version),
        Err(PackError::Malformed { .. })
    ));

    assert!(matches!(
        NormalizedMaterialCardPack::from_bytes(&bytes[..bytes.len() - 1]),
        Err(PackError::Malformed { .. })
    ));

    let mut trailing = bytes.clone();
    trailing.push(0);
    assert!(matches!(
        NormalizedMaterialCardPack::from_bytes(&trailing),
        Err(PackError::Malformed { .. })
    ));

    let wrong_hash = hash_domain(SOURCE_DOMAIN, b"wrong-whole-pack");
    assert!(matches!(
        NormalizedMaterialCardPack::from_bytes_verified(wrong_hash, &bytes),
        Err(PackError::IdentityMismatch {
            kind: "material_card_pack",
            ..
        })
    ));
}

#[test]
fn incomplete_or_nonzero_revision_state_refuses_before_publication() {
    let mut blank_chemistry = aluminum_state();
    blank_chemistry.chemistry = " ".to_string();
    assert!(matches!(
        NormalizedMaterialCardPack::new(blank_chemistry, claims_pack()),
        Err(PackError::InvalidField {
            field: "material_state",
            ..
        })
    ));

    let mut revised = aluminum_state();
    revised.revision = 3;
    assert!(matches!(
        NormalizedMaterialCardPack::new(revised, claims_pack()),
        Err(PackError::MatDb(_))
    ));
}
