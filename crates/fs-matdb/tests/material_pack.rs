//! G0/G3 conformance for the normalized material-card pack boundary.

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MATERIAL_CARD_PACK_SCHEMA_VERSION, MaterialStateId,
    NormalizedMaterialCardPack, NormalizedPack, ObservationDataset, PackError, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
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

#[test]
fn material_card_pack_round_trips_deterministically() {
    let pack = sample_pack();
    let first = pack.to_bytes();
    let second = sample_pack().to_bytes();
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
    assert_eq!(answer.evidence.value.value, 167.0);
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
    bad_version[8..12].copy_from_slice(&(MATERIAL_CARD_PACK_SCHEMA_VERSION + 1).to_le_bytes());
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
