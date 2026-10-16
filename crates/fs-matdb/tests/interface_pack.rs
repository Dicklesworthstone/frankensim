//! G0/G3 conformance for the normalized interface-system pack boundary.

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, INTERFACE_PACK_SCHEMA_VERSION, InterpolationPolicy, MaterialStateId,
    NormalizedInterfacePack, NormalizedPack, ObservationDataset, PackError, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_qty::Dims;

const SOURCE_DOMAIN: &str = "org.frankensim.tests.interface-pack.source.v1";

fn provenance() -> Provenance {
    Provenance {
        source: "pin-on-disk campaign POD-19".to_string(),
        license: "CC-BY-4.0; redistribution permitted with attribution".to_string(),
        artifact: Some(hash_domain(SOURCE_DOMAIN, b"fixture-table")),
    }
}

fn steel(texture_frame: &str) -> SurfaceSpec {
    SurfaceSpec {
        material: MaterialStateId {
            chemistry: "AISI-52100".to_string(),
            phase: "tempered-martensite".to_string(),
            process: "hardened-60HRC".to_string(),
            revision: 2,
        },
        texture_frame: texture_frame.to_string(),
    }
}

fn bronze() -> SurfaceSpec {
    SurfaceSpec {
        material: MaterialStateId {
            chemistry: "C93200".to_string(),
            phase: "cast-bearing-bronze".to_string(),
            process: "machined-bore".to_string(),
            revision: 1,
        },
        texture_frame: "bore-honed-frame-8".to_string(),
    }
}

fn context(history: &str) -> SystemContext {
    SystemContext {
        medium: "oil-film".to_string(),
        third_body: Some("named-reference-oil-lot-4".to_string()),
        environment: "laboratory-air".to_string(),
        history: history.to_string(),
    }
}

fn claims_pack() -> NormalizedPack {
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "ordered steel journal on bronze bearing coupon".to_string(),
            method: "POD-19 run-in campaign".to_string(),
            artifact: hash_domain(SOURCE_DOMAIN, b"raw-observation"),
            caveats: "fixture value; not a seed-dataset authority".to_string(),
            provenance: provenance(),
        })
        .expect("licensed observation inserts");
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("kinetic_friction_coefficient", Dims::NONE),
            value: PropertyValue::Scalar {
                value: 0.08,
                dims: Dims::NONE,
            },
            validity: ValidityDomain::unconstrained()
                .with("temperature", 293.15, 313.15)
                .with("normal_pressure", 1.0e5, 5.0e6),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation],
            provenance: provenance(),
        })
        .expect("interface claim inserts");
    NormalizedPack::new(
        "fixture-steel-bronze-journal-interface",
        "frankensim-interface-pack-compiler-v1",
        hash_domain(SOURCE_DOMAIN, b"source-envelope"),
        "CC-BY-4.0: redistribution permitted with attribution",
        claims,
        Vec::new(),
        Vec::new(),
    )
    .expect("claim pack admits")
}

fn sample_pack() -> NormalizedInterfacePack {
    NormalizedInterfacePack::new(
        steel("journal-ground-frame-3"),
        bronze(),
        context("run-in-1000-cycles"),
        claims_pack(),
    )
    .expect("interface pack admits")
}

#[test]
fn ordered_interface_pack_round_trips_deterministically() {
    let pack = sample_pack();
    let first = pack.to_bytes();
    let second = sample_pack().to_bytes();
    assert_eq!(first, second, "canonical interface bytes moved");
    assert_eq!(&first[..8], b"FSINTPK\0");
    assert_eq!(
        u32::from_le_bytes(first[8..12].try_into().expect("version width")),
        INTERFACE_PACK_SCHEMA_VERSION
    );

    let decoded = NormalizedInterfacePack::from_bytes(&first).expect("interface pack decodes");
    assert_eq!(decoded, pack);
    assert_eq!(decoded.pack_id(), "fixture-steel-bronze-journal-interface");
    assert_eq!(decoded.compiler(), "frankensim-interface-pack-compiler-v1");
    assert_eq!(decoded.card().surface_a(), pack.card().surface_a());
    assert_eq!(decoded.card().surface_b(), pack.card().surface_b());
    assert_eq!(decoded.card().context(), pack.card().context());
    assert_eq!(
        decoded
            .card()
            .claims_for("kinetic_friction_coefficient")
            .len(),
        1
    );
    assert!(decoded.card().models().is_empty(), "v1 carries no models");
    assert_eq!(
        decoded.claims_pack().content_hash(),
        pack.claims_pack().content_hash()
    );
    assert_eq!(decoded.card().content_hash(), pack.card().content_hash());
    assert_eq!(decoded.content_hash(), pack.content_hash());
    assert_eq!(
        NormalizedInterfacePack::from_bytes_verified(pack.content_hash(), &first)
            .expect("whole interface identity verifies"),
        pack
    );
}

#[test]
fn surface_order_and_history_move_pack_identity() {
    let forward = sample_pack();
    let reversed = NormalizedInterfacePack::new(
        bronze(),
        steel("journal-ground-frame-3"),
        context("run-in-1000-cycles"),
        claims_pack(),
    )
    .expect("reversed interface admits");
    let virgin = NormalizedInterfacePack::new(
        steel("journal-ground-frame-3"),
        bronze(),
        context("virgin"),
        claims_pack(),
    )
    .expect("virgin interface admits");

    assert_ne!(
        forward.card().content_hash(),
        reversed.card().content_hash()
    );
    assert_ne!(forward.content_hash(), reversed.content_hash());
    assert_ne!(forward.card().content_hash(), virgin.card().content_hash());
    assert_ne!(forward.content_hash(), virgin.content_hash());
}

#[test]
fn malformed_or_unpinned_interface_artifacts_refuse() {
    let pack = sample_pack();
    let bytes = pack.to_bytes();

    let mut bad_magic = bytes.clone();
    bad_magic[0] ^= 0xff;
    assert!(matches!(
        NormalizedInterfacePack::from_bytes(&bad_magic),
        Err(PackError::Malformed { .. })
    ));

    let mut bad_version = bytes.clone();
    bad_version[8..12].copy_from_slice(&(INTERFACE_PACK_SCHEMA_VERSION + 1).to_le_bytes());
    assert!(matches!(
        NormalizedInterfacePack::from_bytes(&bad_version),
        Err(PackError::Malformed { .. })
    ));

    assert!(matches!(
        NormalizedInterfacePack::from_bytes(&bytes[..bytes.len() - 1]),
        Err(PackError::Malformed { .. })
    ));

    let mut trailing = bytes.clone();
    trailing.push(0);
    assert!(matches!(
        NormalizedInterfacePack::from_bytes(&trailing),
        Err(PackError::Malformed { .. })
    ));

    let wrong_hash = hash_domain(SOURCE_DOMAIN, b"wrong-whole-pack");
    assert!(matches!(
        NormalizedInterfacePack::from_bytes_verified(wrong_hash, &bytes),
        Err(PackError::IdentityMismatch {
            kind: "interface_pack",
            ..
        })
    ));
}

#[test]
fn incomplete_surface_identity_refuses_before_publication() {
    let mut blank_chemistry = steel("journal-ground-frame-3");
    blank_chemistry.material.chemistry = " ".to_string();
    assert!(matches!(
        NormalizedInterfacePack::new(blank_chemistry, bronze(), context("virgin"), claims_pack()),
        Err(PackError::InvalidField {
            field: "interface_surface_material",
            ..
        })
    ));

    let mut blank_third_body = context("virgin");
    blank_third_body.third_body = Some(String::new());
    assert!(matches!(
        NormalizedInterfacePack::new(
            steel("journal-ground-frame-3"),
            bronze(),
            blank_third_body,
            claims_pack()
        ),
        Err(PackError::InvalidField {
            field: "interface_context",
            ..
        })
    ));
}
