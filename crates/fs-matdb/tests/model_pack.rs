//! G0/G3 conformance for normalized constitutive-model-card packs.

use std::collections::BTreeMap;

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, MODEL_PACK_TARGET_BASIS,
    ModelNormalizationReceipt, ModelNormalizationTarget, NormalizedModelPack, PackError,
    Provenance, ValidityBoundSide,
};
use fs_qty::Dims;

const FIXTURE_DOMAIN: &str = "org.frankensim.fs-matdb.model-pack-test.v1";

fn fixture_hash(name: &str) -> ContentHash {
    hash_domain(FIXTURE_DOMAIN, name.as_bytes())
}

fn nasa9_card(source_name: &str) -> ConstitutiveModelCard {
    let coefficient_dims = [
        Dims([0, 0, 0, 2, 0, 0]),
        Dims([0, 0, 0, 1, 0, 0]),
        Dims::NONE,
        Dims([0, 0, 0, -1, 0, 0]),
        Dims([0, 0, 0, -2, 0, 0]),
        Dims([0, 0, 0, -3, 0, 0]),
        Dims([0, 0, 0, -4, 0, 0]),
        Dims([0, 0, 0, 1, 0, 0]),
        Dims::NONE,
    ];
    let mut parameters = BTreeMap::new();
    for (index, dims) in coefficient_dims.into_iter().enumerate() {
        parameters.insert(
            format!("a{index}"),
            LawParameter {
                value: f64::from(index as u32 + 1),
                dims,
            },
        );
    }
    parameters.insert(
        "reference_pressure".to_string(),
        LawParameter {
            value: 100_000.0,
            dims: Dims([-1, 1, -2, 0, 0, 0]),
        },
    );
    let source = fixture_hash(source_name);
    ConstitutiveModelCard {
        law: LawId("nasa9-standard-state".to_string()),
        law_version: 1,
        parameters,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 6_000.0),
        sources: vec![source],
        provenance: Provenance {
            source: format!("fixture NASA-9 table {source_name}"),
            license: "CC-BY-4.0".to_string(),
            artifact: Some(source),
        },
    }
}

fn kinetics_card(source_name: &str) -> ConstitutiveModelCard {
    let source = fixture_hash(source_name);
    ConstitutiveModelCard {
        law: LawId("arrhenius-rate".to_string()),
        law_version: 1,
        parameters: BTreeMap::from([
            (
                "activation_temperature".to_string(),
                LawParameter {
                    value: 12_000.0,
                    dims: Dims([0, 0, 0, 1, 0, 0]),
                },
            ),
            (
                "pre_exponential".to_string(),
                LawParameter {
                    value: 2.5e7,
                    dims: Dims([0, 0, -1, 0, 0, 0]),
                },
            ),
        ]),
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 300.0, 2_000.0),
        sources: vec![source],
        provenance: Provenance {
            source: format!("fixture kinetics table {source_name}"),
            license: "CC-BY-4.0".to_string(),
            artifact: Some(source),
        },
    }
}

fn receipts(card: &ConstitutiveModelCard) -> Vec<ModelNormalizationReceipt> {
    let model = card.content_hash();
    let mut receipts = Vec::new();
    for (name, parameter) in &card.parameters {
        receipts.push(ModelNormalizationReceipt::new(
            ModelNormalizationTarget::Parameter {
                model,
                parameter: name.clone(),
            },
            fixture_hash(&format!("{model}:{name}:literal")),
            parameter.dims,
            1.0,
            0.0,
            "fixture-source-basis",
            MODEL_PACK_TARGET_BASIS,
            None,
            None,
        ));
    }
    for axis in card.validity.bounds().keys() {
        for side in [ValidityBoundSide::Lower, ValidityBoundSide::Upper] {
            receipts.push(ModelNormalizationReceipt::new(
                ModelNormalizationTarget::ValidityBound {
                    model,
                    axis: axis.clone(),
                    side,
                },
                fixture_hash(&format!("{model}:{axis}:{side:?}:literal")),
                Dims([0, 0, 0, 1, 0, 0]),
                1.0,
                0.0,
                "K",
                MODEL_PACK_TARGET_BASIS,
                None,
                None,
            ));
        }
    }
    receipts
}

fn pack(
    cards: Vec<ConstitutiveModelCard>,
    receipts: Vec<ModelNormalizationReceipt>,
) -> Result<NormalizedModelPack, PackError> {
    NormalizedModelPack::new(
        "fixture-models",
        "fixture-compiler-v1",
        fixture_hash("raw-source-envelope"),
        "redistribution permitted with attribution",
        cards,
        receipts,
    )
}

#[test]
fn g0_model_cards_round_trip_and_canonicalize_input_order() {
    let nasa = nasa9_card("nasa");
    let kinetics = kinetics_card("kinetics");
    let mut all_receipts = receipts(&nasa);
    all_receipts.extend(receipts(&kinetics));
    let first =
        pack(vec![nasa.clone(), kinetics.clone()], all_receipts.clone()).expect("first model pack");
    all_receipts.reverse();
    let second = pack(vec![kinetics, nasa], all_receipts).expect("permuted model pack");

    assert_eq!(first, second);
    assert_eq!(first.to_bytes(), second.to_bytes());
    assert_eq!(first.content_hash(), second.content_hash());
    assert_eq!(first.models().len(), 2);
    assert!(
        first
            .models()
            .windows(2)
            .all(|pair| pair[0].content_hash() < pair[1].content_hash())
    );
    assert_eq!(first.normalizations().len(), 16);

    let bytes = first.to_bytes();
    let decoded = NormalizedModelPack::from_bytes_verified(first.content_hash(), &bytes)
        .expect("verified model pack decodes");
    assert_eq!(decoded, first);
    assert_eq!(decoded.to_bytes(), bytes);
}

#[test]
fn g3_every_model_numeric_field_requires_one_exact_receipt() {
    let card = nasa9_card("coverage");
    let mut complete = receipts(&card);
    let missing = complete.pop().expect("one receipt");
    let error = pack(vec![card.clone()], complete.clone()).expect_err("missing receipt refuses");
    assert!(matches!(
        error,
        PackError::InvalidField {
            field: "model_normalizations",
            ..
        }
    ));

    complete.push(missing.clone());
    complete.push(missing);
    let error = pack(vec![card.clone()], complete).expect_err("duplicate receipt refuses");
    assert!(matches!(
        error,
        PackError::InvalidField {
            field: "model_normalizations",
            ..
        }
    ));

    let mut wrong_dims = receipts(&card);
    let target = wrong_dims[0].target().clone();
    wrong_dims[0] = ModelNormalizationReceipt::new(
        target,
        fixture_hash("wrong-dims"),
        Dims([1, 0, 0, 0, 0, 0]),
        1.0,
        0.0,
        "fixture-source-basis",
        MODEL_PACK_TARGET_BASIS,
        None,
        None,
    );
    let error = pack(vec![card], wrong_dims).expect_err("dimension mismatch refuses");
    assert!(matches!(
        error,
        PackError::InvalidField {
            field: "model_normalization.dims",
            ..
        }
    ));
}

#[test]
fn g3_validity_endpoint_receipts_must_share_one_transform() {
    let card = kinetics_card("coherence");
    let mut all = receipts(&card);
    let upper_index = all
        .iter()
        .position(|receipt| {
            matches!(
                receipt.target(),
                ModelNormalizationTarget::ValidityBound {
                    side: ValidityBoundSide::Upper,
                    ..
                }
            )
        })
        .expect("upper validity receipt");
    let target = all[upper_index].target().clone();
    all[upper_index] = ModelNormalizationReceipt::new(
        target,
        fixture_hash("incoherent-upper"),
        Dims([0, 0, 0, 1, 0, 0]),
        2.0,
        0.0,
        "K",
        MODEL_PACK_TARGET_BASIS,
        None,
        None,
    );

    let error = pack(vec![card], all).expect_err("incoherent endpoints refuse");
    assert!(matches!(
        error,
        PackError::InvalidField {
            field: "model_normalizations",
            ..
        }
    ));
}

#[test]
fn g3_portable_profile_refuses_noncanonical_cards() {
    let mut zero_version = kinetics_card("zero-version");
    zero_version.law_version = 0;
    let error = pack(vec![zero_version], Vec::new()).expect_err("zero law version refuses");
    assert!(matches!(
        error,
        PackError::InvalidField {
            field: "model.law_version",
            ..
        }
    ));

    let mut negative_zero = kinetics_card("negative-zero");
    negative_zero
        .parameters
        .get_mut("pre_exponential")
        .expect("fixture parameter")
        .value = -0.0;
    let error = pack(vec![negative_zero], Vec::new()).expect_err("negative zero refuses");
    assert!(matches!(
        error,
        PackError::InvalidField {
            field: "model.parameter.value",
            ..
        }
    ));

    let mut missing_artifact = kinetics_card("missing-artifact");
    missing_artifact.provenance.artifact = None;
    let error = pack(vec![missing_artifact], Vec::new()).expect_err("artifact-less card refuses");
    assert!(matches!(
        error,
        PackError::InvalidField {
            field: "model.provenance.artifact",
            ..
        }
    ));

    let mut unsorted_sources = kinetics_card("unsorted-sources");
    unsorted_sources.sources = vec![fixture_hash("z"), fixture_hash("a")];
    if unsorted_sources.sources[0] < unsorted_sources.sources[1] {
        unsorted_sources.sources.reverse();
    }
    let error = pack(vec![unsorted_sources], Vec::new()).expect_err("source order refuses");
    assert!(matches!(
        error,
        PackError::InvalidField {
            field: "model.sources",
            ..
        }
    ));
}

#[test]
fn g3_whole_pack_and_nested_card_identities_fail_closed() {
    let card = nasa9_card("tamper");
    let admitted = pack(vec![card.clone()], receipts(&card)).expect("fixture pack");
    let expected = admitted.content_hash();
    let bytes = admitted.to_bytes();

    let mut top_level_tamper = bytes.clone();
    let pack_id = b"fixture-models";
    let pack_id_offset = top_level_tamper
        .windows(pack_id.len())
        .position(|window| window == pack_id)
        .expect("pack id bytes");
    top_level_tamper[pack_id_offset] = b'F';
    assert!(NormalizedModelPack::from_bytes(&top_level_tamper).is_ok());
    assert!(matches!(
        NormalizedModelPack::from_bytes_verified(expected, &top_level_tamper),
        Err(PackError::IdentityMismatch {
            kind: "model pack",
            ..
        })
    ));

    let mut nested_tamper = bytes.clone();
    let model_identity = card.content_hash();
    let model_offset = nested_tamper
        .windows(model_identity.as_bytes().len())
        .position(|window| window == model_identity.as_bytes())
        .expect("nested model identity");
    nested_tamper[model_offset] ^= 1;
    assert!(matches!(
        NormalizedModelPack::from_bytes(&nested_tamper),
        Err(PackError::IdentityMismatch {
            kind: "model card",
            ..
        })
    ));

    let mut trailing = bytes;
    trailing.push(0);
    assert!(matches!(
        NormalizedModelPack::from_bytes(&trailing),
        Err(PackError::Malformed { .. })
    ));
}
