//! G0/G3 conformance for normalized species-association packs.

use fs_blake3::{ContentHash, hash_domain};
use fs_matdb::{
    NormalizedSpeciesPack, PackError, Provenance, SPECIES_MOLAR_MASS_DIMS,
    SPECIES_PACK_TARGET_BASIS, SPECIES_REFERENCE_PRESSURE_DIMS, SpeciesAssociation,
    SpeciesNormalizationReceipt, SpeciesNormalizationTarget,
};
use fs_qty::Dims;
use fs_qty::chemistry::SpeciesId;

const FIXTURE_DOMAIN: &str = "org.frankensim.fs-matdb.species-pack-test.v1";

fn fixture_hash(name: &str) -> ContentHash {
    hash_domain(FIXTURE_DOMAIN, name.as_bytes())
}

fn association(sources: Vec<ContentHash>) -> Result<SpeciesAssociation, PackError> {
    SpeciesAssociation::new(
        SpeciesId::new("N2").expect("canonical species"),
        0.028_013_4,
        "gas",
        "ideal-gas",
        100_000.0,
        "NASA-TP-2002-211556",
        sources,
        Provenance {
            source: "licensed nitrogen species metadata fixture".to_string(),
            license: "CC-BY-4.0".to_string(),
            artifact: Some(fixture_hash("source-a")),
        },
    )
}

fn receipts() -> Vec<SpeciesNormalizationReceipt> {
    vec![
        SpeciesNormalizationReceipt::new(
            SpeciesNormalizationTarget::ReferencePressure,
            fixture_hash("pressure-literal"),
            SPECIES_REFERENCE_PRESSURE_DIMS,
            1.0,
            0.0,
            "Pa",
            SPECIES_PACK_TARGET_BASIS,
        ),
        SpeciesNormalizationReceipt::new(
            SpeciesNormalizationTarget::MolarMass,
            fixture_hash("molar-mass-literal"),
            SPECIES_MOLAR_MASS_DIMS,
            0.001,
            0.0,
            "g/mol",
            SPECIES_PACK_TARGET_BASIS,
        ),
    ]
}

fn pack(
    association: SpeciesAssociation,
    receipts: Vec<SpeciesNormalizationReceipt>,
) -> Result<NormalizedSpeciesPack, PackError> {
    NormalizedSpeciesPack::new(
        "N2",
        "fixture-species-compiler-v1",
        fixture_hash("source-envelope"),
        "redistribution permitted with attribution",
        association,
        receipts,
    )
}

#[test]
fn g0_species_association_round_trips_and_canonicalizes_sources_and_receipts() {
    let source_a = fixture_hash("source-a");
    let source_b = fixture_hash("source-b");
    let first = pack(
        association(vec![source_b, source_a]).expect("association"),
        receipts(),
    )
    .expect("species pack");
    let mut reversed_receipts = receipts();
    reversed_receipts.reverse();
    let second = pack(
        association(vec![source_a, source_b]).expect("association"),
        reversed_receipts,
    )
    .expect("permuted species pack");

    assert_eq!(first, second);
    assert_eq!(first.to_bytes(), second.to_bytes());
    assert_eq!(first.content_hash(), second.content_hash());
    assert_eq!(first.pack_id(), "N2");
    assert_eq!(first.association().species().as_str(), "N2");
    assert_eq!(
        first.association().molar_mass().to_bits(),
        0.028_013_4f64.to_bits()
    );
    assert_eq!(first.association().standard_state_phase(), "gas");
    assert_eq!(first.association().reference_eos(), "ideal-gas");
    assert_eq!(
        first.association().reference_pressure().to_bits(),
        100_000.0f64.to_bits()
    );
    assert_eq!(
        first.association().elemental_reference(),
        "NASA-TP-2002-211556"
    );
    let mut expected_sources = vec![source_a, source_b];
    expected_sources.sort_unstable();
    assert_eq!(first.association().sources(), expected_sources.as_slice());
    assert_eq!(first.normalizations().len(), 2);
    assert_eq!(
        first.normalizations()[0].target(),
        SpeciesNormalizationTarget::MolarMass
    );
    assert_eq!(
        first.normalizations()[1].target(),
        SpeciesNormalizationTarget::ReferencePressure
    );

    let bytes = first.to_bytes();
    let decoded = NormalizedSpeciesPack::from_bytes_verified(first.content_hash(), &bytes)
        .expect("verified species pack decodes");
    assert_eq!(decoded, first);
    assert_eq!(decoded.to_bytes(), bytes);
}

#[test]
#[allow(clippy::too_many_lines)] // One table keeps related runtime refusal gates visible together.
fn g3_species_association_and_receipt_admission_fail_closed() {
    let source_a = fixture_hash("source-a");
    for (phase, eos, molar_mass, pressure, field) in [
        (
            "liquid",
            "ideal-gas",
            0.028,
            100_000.0,
            "species.standard_state_phase",
        ),
        ("gas", "cubic", 0.028, 100_000.0, "species.reference_eos"),
        ("gas", "ideal-gas", 0.0, 100_000.0, "species.molar_mass"),
        (
            "gas",
            "ideal-gas",
            0.028,
            f64::NAN,
            "species.reference_pressure",
        ),
    ] {
        let error = SpeciesAssociation::new(
            SpeciesId::new("N2").expect("species"),
            molar_mass,
            phase,
            eos,
            pressure,
            "NASA-ref",
            vec![source_a],
            Provenance {
                source: "fixture".to_string(),
                license: "CC-BY-4.0".to_string(),
                artifact: Some(source_a),
            },
        )
        .expect_err("invalid association refuses");
        assert!(matches!(error, PackError::InvalidField { field: found, .. } if found == field));
    }

    let missing_artifact = SpeciesAssociation::new(
        SpeciesId::new("N2").expect("species"),
        0.028,
        "gas",
        "ideal-gas",
        100_000.0,
        "NASA-ref",
        vec![source_a],
        Provenance {
            source: "fixture".to_string(),
            license: "CC-BY-4.0".to_string(),
            artifact: None,
        },
    )
    .expect_err("missing provenance artifact refuses");
    assert!(matches!(
        missing_artifact,
        PackError::InvalidField {
            field: "species.provenance.artifact",
            ..
        }
    ));

    let duplicate_source =
        association(vec![source_a, source_a]).expect_err("duplicate source artifacts refuse");
    assert!(matches!(
        duplicate_source,
        PackError::InvalidField {
            field: "species.sources",
            ..
        }
    ));

    let source_b = fixture_hash("source-b");
    let unretained_artifact = SpeciesAssociation::new(
        SpeciesId::new("N2").expect("species"),
        0.028,
        "gas",
        "ideal-gas",
        100_000.0,
        "NASA-ref",
        vec![source_b],
        Provenance {
            source: "fixture".to_string(),
            license: "CC-BY-4.0".to_string(),
            artifact: Some(source_a),
        },
    )
    .expect_err("unretained provenance artifact refuses");
    assert!(matches!(
        unretained_artifact,
        PackError::InvalidField {
            field: "species.provenance.artifact",
            ..
        }
    ));

    let valid_association = association(vec![source_a]).expect("valid association");
    let mut one_receipt = receipts();
    one_receipt.truncate(1);
    assert!(matches!(
        pack(valid_association.clone(), one_receipt),
        Err(PackError::InvalidField {
            field: "species_normalizations",
            ..
        })
    ));

    let mut duplicated = receipts();
    duplicated[1] = duplicated[0].clone();
    assert!(matches!(
        pack(valid_association.clone(), duplicated),
        Err(PackError::InvalidField {
            field: "species_normalizations",
            ..
        })
    ));

    let mut wrong_dims = receipts();
    wrong_dims[0] = SpeciesNormalizationReceipt::new(
        wrong_dims[0].target(),
        fixture_hash("wrong-dims"),
        Dims::NONE,
        1.0,
        0.0,
        "Pa",
        SPECIES_PACK_TARGET_BASIS,
    );
    assert!(matches!(
        pack(valid_association.clone(), wrong_dims),
        Err(PackError::InvalidField {
            field: "species_normalization.dims",
            ..
        })
    ));

    for (offset, name) in [(273.15, "affine"), (-0.0, "negative-zero")] {
        let mut translated = receipts();
        translated[0] = SpeciesNormalizationReceipt::new(
            translated[0].target(),
            fixture_hash(name),
            translated[0].dims(),
            1.0,
            offset,
            "translated-pressure",
            SPECIES_PACK_TARGET_BASIS,
        );
        assert!(matches!(
            pack(valid_association.clone(), translated),
            Err(PackError::InvalidField {
                field: "species_normalization.offset",
                ..
            })
        ));
    }

    let wrong_pack_id = NormalizedSpeciesPack::new(
        "O2",
        "fixture-species-compiler-v1",
        fixture_hash("source-envelope"),
        "permitted",
        valid_association,
        receipts(),
    );
    assert!(matches!(
        wrong_pack_id,
        Err(PackError::InvalidField {
            field: "pack_id",
            ..
        })
    ));
}

#[test]
fn g3_species_pack_identity_resource_and_canonical_byte_barriers_hold() {
    let admitted = pack(
        association(vec![fixture_hash("source-a")]).expect("association"),
        receipts(),
    )
    .expect("pack");
    let expected = admitted.content_hash();
    let bytes = admitted.to_bytes();

    let mut top_level_tamper = bytes.clone();
    // Magic + version + pack-id length + "N2" + compiler length.
    let compiler_first_byte = 8 + 4 + 4 + 2 + 4;
    top_level_tamper[compiler_first_byte] ^= 1;
    assert!(NormalizedSpeciesPack::from_bytes(&top_level_tamper).is_ok());
    assert!(matches!(
        NormalizedSpeciesPack::from_bytes_verified(expected, &top_level_tamper),
        Err(PackError::IdentityMismatch {
            kind: "species pack",
            ..
        })
    ));

    let mut oversized_string = bytes.clone();
    oversized_string[12..16].copy_from_slice(&u32::MAX.to_le_bytes());
    assert!(matches!(
        NormalizedSpeciesPack::from_bytes(&oversized_string),
        Err(PackError::ResourceLimit {
            resource: "string_bytes",
            ..
        })
    ));

    let mut trailing = bytes;
    trailing.push(0);
    assert!(matches!(
        NormalizedSpeciesPack::from_bytes(&trailing),
        Err(PackError::Malformed { .. })
    ));
}
