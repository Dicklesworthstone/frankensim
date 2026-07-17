//! G0 conformance for dependency-neutral, exact receipt-schema metadata.

use fs_blake3::{ContentHash, hash_domain};
use fs_package::{
    MAX_RECEIPT_FAMILY_ID_BYTES, MAX_RECEIPT_IDENTITY_DOMAIN_BYTES,
    MAX_RECEIPT_SCHEMA_CATALOG_BYTES, MAX_RECEIPT_SCHEMA_ENTRIES, MAX_RECEIPT_TRANSPORT_BYTES,
    RECEIPT_SCHEMA_CATALOG_IDENTITY_DOMAIN, RECEIPT_SCHEMA_CATALOG_IDENTITY_VERSION,
    RECEIPT_SCHEMA_CATALOG_VERSION, RECEIPT_SCHEMA_DESCRIPTOR_IDENTITY_DOMAIN,
    RECEIPT_SCHEMA_DESCRIPTOR_IDENTITY_VERSION, ReceiptSchemaCatalog, ReceiptSchemaCatalogError,
    ReceiptSchemaDescriptor, ReceiptTransportProfile,
};

const CATALOG_MAGIC: &[u8; 8] = b"FSPRCAT\0";
const FIELD_FAMILY: u8 = 1;
const FIELD_WIRE_SCHEMA_VERSION: u8 = 2;
const FIELD_OWNER_IDENTITY_VERSION: u8 = 3;
const FIELD_OWNER_IDENTITY_DOMAIN: u8 = 4;
const FIELD_TRANSPORT: u8 = 5;
const FIELD_OWNER_SCHEMA_FINGERPRINT: u8 = 6;
const FIELD_DESCRIPTOR_HASH: u8 = 7;
const TRANSPORT_DIGEST_ONLY: u8 = 1;
const TRANSPORT_CANONICAL_BYTES: u8 = 2;

fn fingerprint(byte: u8) -> ContentHash {
    ContentHash([byte; 32])
}

fn descriptor(
    family: &str,
    wire_schema_version: u32,
    owner_identity_version: u32,
    owner_identity_domain: &str,
    transport: ReceiptTransportProfile,
    fingerprint_byte: u8,
) -> ReceiptSchemaDescriptor {
    ReceiptSchemaDescriptor::try_new(
        family,
        wire_schema_version,
        owner_identity_version,
        owner_identity_domain,
        transport,
        fingerprint(fingerprint_byte),
    )
    .expect("fixture descriptor is canonical")
}

fn matdb_v2() -> ReceiptSchemaDescriptor {
    descriptor(
        "fs-matdb:property-usage-receipt",
        2,
        2,
        "org.frankensim.fs-matdb.property-usage-receipt.v2",
        ReceiptTransportProfile::CanonicalBytes {
            maximum_bytes: 1024 * 1024,
        },
        0x51,
    )
}

fn matdb_v1() -> ReceiptSchemaDescriptor {
    descriptor(
        "fs-matdb:property-usage-receipt",
        1,
        1,
        "org.frankensim.fs-matdb.property-usage-receipt.v1",
        ReceiptTransportProfile::DigestOnly,
        0x41,
    )
}

fn ledger_checkpoint_v1() -> ReceiptSchemaDescriptor {
    descriptor(
        "fs-ledger:state-checkpoint-receipt",
        1,
        1,
        "org.frankensim.fs-ledger.state-checkpoint-receipt.v1",
        ReceiptTransportProfile::DigestOnly,
        0x61,
    )
}

fn push_u32(output: &mut Vec<u8>, value: u32) {
    output.extend_from_slice(&value.to_le_bytes());
}

fn push_u64(output: &mut Vec<u8>, value: u64) {
    output.extend_from_slice(&value.to_le_bytes());
}

fn push_string(output: &mut Vec<u8>, value: &str) {
    push_u64(
        output,
        u64::try_from(value.len()).expect("fixture string length fits u64"),
    );
    output.extend_from_slice(value.as_bytes());
}

fn independent_descriptor_preimage(descriptor: &ReceiptSchemaDescriptor) -> Vec<u8> {
    let mut output = Vec::new();
    output.push(FIELD_FAMILY);
    push_string(&mut output, descriptor.family_id());
    output.push(FIELD_WIRE_SCHEMA_VERSION);
    push_u32(&mut output, descriptor.wire_schema_version());
    output.push(FIELD_OWNER_IDENTITY_VERSION);
    push_u32(&mut output, descriptor.owner_identity_version());
    output.push(FIELD_OWNER_IDENTITY_DOMAIN);
    push_string(&mut output, descriptor.owner_identity_domain());
    output.push(FIELD_TRANSPORT);
    match descriptor.transport() {
        ReceiptTransportProfile::DigestOnly => {
            output.push(TRANSPORT_DIGEST_ONLY);
            push_u64(&mut output, 0);
        }
        ReceiptTransportProfile::CanonicalBytes { maximum_bytes } => {
            output.push(TRANSPORT_CANONICAL_BYTES);
            push_u64(&mut output, maximum_bytes);
        }
    }
    output.push(FIELD_OWNER_SCHEMA_FINGERPRINT);
    output.extend_from_slice(descriptor.owner_schema_fingerprint().as_bytes());
    output
}

fn independent_descriptor_hash(descriptor: &ReceiptSchemaDescriptor) -> ContentHash {
    let fields = independent_descriptor_preimage(descriptor);
    let mut preimage = Vec::with_capacity(4 + fields.len());
    push_u32(&mut preimage, RECEIPT_SCHEMA_DESCRIPTOR_IDENTITY_VERSION);
    preimage.extend_from_slice(&fields);
    hash_domain(RECEIPT_SCHEMA_DESCRIPTOR_IDENTITY_DOMAIN, &preimage)
}

fn independent_descriptor_row(descriptor: &ReceiptSchemaDescriptor) -> Vec<u8> {
    let mut output = independent_descriptor_preimage(descriptor);
    output.push(FIELD_DESCRIPTOR_HASH);
    output.extend_from_slice(independent_descriptor_hash(descriptor).as_bytes());
    output
}

fn independent_catalog_wire(entries: &[ReceiptSchemaDescriptor]) -> Vec<u8> {
    let mut output = Vec::new();
    output.extend_from_slice(CATALOG_MAGIC);
    push_u32(&mut output, RECEIPT_SCHEMA_CATALOG_VERSION);
    push_u64(
        &mut output,
        u64::try_from(entries.len()).expect("fixture entry count fits u64"),
    );
    for descriptor in entries {
        output.extend_from_slice(&independent_descriptor_row(descriptor));
    }
    output
}

fn independent_catalog_hash(wire: &[u8]) -> ContentHash {
    let mut preimage = Vec::with_capacity(4 + wire.len());
    push_u32(&mut preimage, RECEIPT_SCHEMA_CATALOG_IDENTITY_VERSION);
    preimage.extend_from_slice(wire);
    hash_domain(RECEIPT_SCHEMA_CATALOG_IDENTITY_DOMAIN, &preimage)
}

fn read_usize(bytes: &[u8], offset: usize) -> usize {
    let raw: [u8; 8] = bytes[offset..offset + 8]
        .try_into()
        .expect("fixture contains a complete u64");
    usize::try_from(u64::from_le_bytes(raw)).expect("fixture u64 fits usize")
}

fn first_transport_tag_offset(bytes: &[u8]) -> usize {
    let mut cursor = 20;
    assert_eq!(bytes[cursor], FIELD_FAMILY);
    cursor += 1;
    let family_len = read_usize(bytes, cursor);
    cursor += 8 + family_len;
    assert_eq!(bytes[cursor], FIELD_WIRE_SCHEMA_VERSION);
    cursor += 1 + 4;
    assert_eq!(bytes[cursor], FIELD_OWNER_IDENTITY_VERSION);
    cursor += 1 + 4;
    assert_eq!(bytes[cursor], FIELD_OWNER_IDENTITY_DOMAIN);
    cursor += 1;
    let domain_len = read_usize(bytes, cursor);
    cursor += 8 + domain_len;
    assert_eq!(bytes[cursor], FIELD_TRANSPORT);
    cursor + 1
}

fn assert_descriptor_identity_moves(
    baseline: &ReceiptSchemaDescriptor,
    changed: &ReceiptSchemaDescriptor,
) {
    assert_ne!(baseline.content_hash(), changed.content_hash());
}

#[test]
#[allow(clippy::too_many_lines)]
fn receipt_schema_descriptor_identity_binds_every_field() {
    let baseline = matdb_v2();
    assert_eq!(
        baseline.content_hash(),
        independent_descriptor_hash(&baseline)
    );

    assert_descriptor_identity_moves(
        &baseline,
        &descriptor(
            "fs-matdb:property-usage-recordx",
            2,
            2,
            baseline.owner_identity_domain(),
            baseline.transport(),
            0x51,
        ),
    );
    assert_descriptor_identity_moves(
        &baseline,
        &descriptor(
            "fs-matdb:property-usage-receipt-extended",
            2,
            2,
            baseline.owner_identity_domain(),
            baseline.transport(),
            0x51,
        ),
    );
    assert_descriptor_identity_moves(
        &baseline,
        &descriptor(
            baseline.family_id(),
            2,
            2,
            "org.frankensim.fs-matdb.property-usage-receipt-extended.v2",
            baseline.transport(),
            0x51,
        ),
    );
    assert_descriptor_identity_moves(
        &baseline,
        &descriptor(
            baseline.family_id(),
            3,
            2,
            baseline.owner_identity_domain(),
            baseline.transport(),
            0x51,
        ),
    );
    assert_descriptor_identity_moves(
        &baseline,
        &descriptor(
            baseline.family_id(),
            2,
            3,
            baseline.owner_identity_domain(),
            baseline.transport(),
            0x51,
        ),
    );
    assert_descriptor_identity_moves(
        &baseline,
        &descriptor(
            baseline.family_id(),
            2,
            2,
            "org.frankensim.fs-matdb.property-usage-receipt.v3",
            baseline.transport(),
            0x51,
        ),
    );
    assert_descriptor_identity_moves(
        &baseline,
        &descriptor(
            baseline.family_id(),
            2,
            2,
            baseline.owner_identity_domain(),
            ReceiptTransportProfile::DigestOnly,
            0x51,
        ),
    );
    assert_descriptor_identity_moves(
        &baseline,
        &descriptor(
            baseline.family_id(),
            2,
            2,
            baseline.owner_identity_domain(),
            ReceiptTransportProfile::CanonicalBytes {
                maximum_bytes: 1024 * 1024 + 1,
            },
            0x51,
        ),
    );
    assert_descriptor_identity_moves(
        &baseline,
        &descriptor(
            baseline.family_id(),
            2,
            2,
            baseline.owner_identity_domain(),
            baseline.transport(),
            0x52,
        ),
    );

    let retained = baseline.content_hash();
    assert_eq!(
        ReceiptSchemaDescriptor::admit_retained_content_hash(
            RECEIPT_SCHEMA_DESCRIPTOR_IDENTITY_VERSION,
            retained.as_bytes(),
        ),
        Some(retained)
    );
    assert_eq!(
        ReceiptSchemaDescriptor::admit_retained_content_hash(
            RECEIPT_SCHEMA_DESCRIPTOR_IDENTITY_VERSION + 1,
            retained.as_bytes(),
        ),
        None
    );
    assert_eq!(
        ReceiptSchemaDescriptor::admit_retained_content_hash(
            RECEIPT_SCHEMA_DESCRIPTOR_IDENTITY_VERSION,
            &retained.as_bytes()[..31],
        ),
        None
    );
}

#[test]
fn receipt_schema_catalog_round_trips_and_locks_the_independent_preimage() {
    let catalog = ReceiptSchemaCatalog::try_new(vec![matdb_v2(), ledger_checkpoint_v1()])
        .expect("valid catalog");
    let wire = independent_catalog_wire(catalog.entries());
    let expected_hash = independent_catalog_hash(&wire);
    let mut expected_bytes = wire.clone();
    expected_bytes.extend_from_slice(expected_hash.as_bytes());

    assert_eq!(catalog.catalog_version(), RECEIPT_SCHEMA_CATALOG_VERSION);
    assert_eq!(catalog.content_hash(), expected_hash);
    assert_eq!(catalog.to_bytes(), expected_bytes);
    assert_eq!(
        &wire[12..20],
        &u64::try_from(catalog.entries().len())
            .expect("entry count fits u64")
            .to_le_bytes()
    );
    for descriptor in catalog.entries() {
        assert_eq!(
            descriptor.content_hash(),
            independent_descriptor_hash(descriptor)
        );
    }

    let decoded = ReceiptSchemaCatalog::from_bytes(&expected_bytes).expect("round trip");
    assert_eq!(decoded, catalog);
    assert_eq!(decoded.to_bytes(), expected_bytes);
    assert_eq!(
        ReceiptSchemaCatalog::from_bytes_verified(&expected_bytes, expected_hash)
            .expect("external pin"),
        catalog
    );
    assert_eq!(
        ReceiptSchemaCatalog::admit_retained_content_hash(
            RECEIPT_SCHEMA_CATALOG_IDENTITY_VERSION,
            expected_hash.as_bytes(),
        ),
        Some(expected_hash)
    );
    assert_eq!(
        ReceiptSchemaCatalog::admit_retained_content_hash(
            RECEIPT_SCHEMA_CATALOG_IDENTITY_VERSION + 1,
            expected_hash.as_bytes(),
        ),
        None
    );
    assert_eq!(
        ReceiptSchemaCatalog::admit_retained_content_hash(
            RECEIPT_SCHEMA_CATALOG_IDENTITY_VERSION,
            &expected_hash.as_bytes()[..31],
        ),
        None
    );
    let mut extended_hash = expected_hash.as_bytes().to_vec();
    extended_hash.push(0);
    assert_eq!(
        ReceiptSchemaCatalog::admit_retained_content_hash(
            RECEIPT_SCHEMA_CATALOG_IDENTITY_VERSION,
            &extended_hash,
        ),
        None
    );

    let mut changed_count_wire = wire.clone();
    changed_count_wire[12] ^= 1;
    assert_ne!(
        independent_catalog_hash(&changed_count_wire),
        expected_hash,
        "entry count must enter the catalog identity"
    );

    let descriptor = &catalog.entries()[0];
    let fields = independent_descriptor_preimage(descriptor);
    let mut changed_length_fields = fields.clone();
    changed_length_fields[1] ^= 1;
    let mut changed_length_preimage = Vec::new();
    push_u32(
        &mut changed_length_preimage,
        RECEIPT_SCHEMA_DESCRIPTOR_IDENTITY_VERSION,
    );
    changed_length_preimage.extend_from_slice(&changed_length_fields);
    assert_ne!(
        hash_domain(
            RECEIPT_SCHEMA_DESCRIPTOR_IDENTITY_DOMAIN,
            &changed_length_preimage,
        ),
        descriptor.content_hash(),
        "string byte counts must enter descriptor identity"
    );
}

#[test]
fn receipt_schema_catalog_input_order_is_nonsemantic_and_lookup_is_exact() {
    let v1 = matdb_v1();
    let v2 = matdb_v2();
    let forward =
        ReceiptSchemaCatalog::try_new(vec![v1.clone(), v2.clone()]).expect("forward catalog");
    let reverse =
        ReceiptSchemaCatalog::try_new(vec![v2.clone(), v1.clone()]).expect("reverse catalog");

    assert_eq!(forward, reverse);
    assert_eq!(forward.content_hash(), reverse.content_hash());
    assert_eq!(forward.to_bytes(), reverse.to_bytes());
    assert_eq!(forward.entries()[0].wire_schema_version(), 1);
    assert_eq!(forward.entries()[1].wire_schema_version(), 2);

    let exact = forward
        .require_exact(v2.family_id(), 2, v2.content_hash())
        .expect("exact row");
    assert_eq!(exact, &v2);
    assert!(matches!(
        forward.require_exact("fs-absent:receipt", 1, fingerprint(0x99)),
        Err(ReceiptSchemaCatalogError::UnknownFamily { .. })
    ));
    assert!(matches!(
        forward.require_exact(v2.family_id(), 3, fingerprint(0x99)),
        Err(ReceiptSchemaCatalogError::UnsupportedWireSchema { found: 3, .. })
    ));
    assert!(matches!(
        forward.require_exact(v2.family_id(), 2, v1.content_hash()),
        Err(ReceiptSchemaCatalogError::DescriptorMismatch {
            wire_schema_version: 2,
            ..
        })
    ));
}

#[test]
#[allow(clippy::too_many_lines)]
fn receipt_schema_catalog_decoder_refuses_hostile_and_noncanonical_bytes() {
    let canonical = ReceiptSchemaCatalog::try_new(vec![matdb_v2()]).expect("catalog");
    let bytes = canonical.to_bytes();

    for end in 0..bytes.len() {
        assert!(
            ReceiptSchemaCatalog::from_bytes(&bytes[..end]).is_err(),
            "truncation at {end} must refuse"
        );
    }

    let mut wrong_magic = bytes.clone();
    wrong_magic[0] ^= 1;
    assert!(matches!(
        ReceiptSchemaCatalog::from_bytes(&wrong_magic),
        Err(ReceiptSchemaCatalogError::Malformed { .. })
    ));

    let mut wrong_version = bytes.clone();
    wrong_version[8..12].copy_from_slice(&(RECEIPT_SCHEMA_CATALOG_VERSION + 1).to_le_bytes());
    assert!(matches!(
        ReceiptSchemaCatalog::from_bytes(&wrong_version),
        Err(ReceiptSchemaCatalogError::UnsupportedCatalogVersion { .. })
    ));

    let mut hostile_count = bytes.clone();
    hostile_count[12..20].copy_from_slice(&u64::MAX.to_le_bytes());
    assert!(matches!(
        ReceiptSchemaCatalog::from_bytes(&hostile_count),
        Err(ReceiptSchemaCatalogError::ResourceLimit {
            resource: "schema-entries",
            ..
        })
    ));

    let mut wrong_field_tag = bytes.clone();
    wrong_field_tag[20] = 0xff;
    assert!(matches!(
        ReceiptSchemaCatalog::from_bytes(&wrong_field_tag),
        Err(ReceiptSchemaCatalogError::Malformed { .. })
    ));

    let mut hostile_family_length = bytes.clone();
    hostile_family_length[21..29].copy_from_slice(&u64::MAX.to_le_bytes());
    assert!(matches!(
        ReceiptSchemaCatalog::from_bytes(&hostile_family_length),
        Err(ReceiptSchemaCatalogError::ResourceLimit {
            resource: "family-id",
            ..
        })
    ));

    let mut non_utf8_family = bytes.clone();
    non_utf8_family[29] = 0xff;
    assert!(matches!(
        ReceiptSchemaCatalog::from_bytes(&non_utf8_family),
        Err(ReceiptSchemaCatalogError::Malformed { .. })
    ));

    let transport_tag = first_transport_tag_offset(&bytes);
    let mut unknown_transport = bytes.clone();
    unknown_transport[transport_tag] = 0xff;
    assert!(matches!(
        ReceiptSchemaCatalog::from_bytes(&unknown_transport),
        Err(ReceiptSchemaCatalogError::UnknownTransportTag { tag: 0xff, .. })
    ));

    let digest_catalog = ReceiptSchemaCatalog::try_new(vec![matdb_v1()]).expect("digest catalog");
    let mut nonzero_digest_limit = digest_catalog.to_bytes();
    let digest_transport_tag = first_transport_tag_offset(&nonzero_digest_limit);
    nonzero_digest_limit[digest_transport_tag + 1..digest_transport_tag + 9]
        .copy_from_slice(&1_u64.to_le_bytes());
    assert!(matches!(
        ReceiptSchemaCatalog::from_bytes(&nonzero_digest_limit),
        Err(ReceiptSchemaCatalogError::Malformed { .. })
    ));

    let mut descriptor_tamper = bytes.clone();
    let descriptor_hash_offset = descriptor_tamper.len() - 64;
    descriptor_tamper[descriptor_hash_offset] ^= 1;
    assert!(matches!(
        ReceiptSchemaCatalog::from_bytes(&descriptor_tamper),
        Err(ReceiptSchemaCatalogError::IdentityMismatch {
            scope: "descriptor",
            ..
        })
    ));

    let mut catalog_tamper = bytes.clone();
    let catalog_hash_offset = catalog_tamper.len() - 32;
    catalog_tamper[catalog_hash_offset] ^= 1;
    assert!(matches!(
        ReceiptSchemaCatalog::from_bytes(&catalog_tamper),
        Err(ReceiptSchemaCatalogError::IdentityMismatch {
            scope: "catalog",
            ..
        })
    ));

    let mut trailing = bytes.clone();
    trailing.push(0);
    assert!(matches!(
        ReceiptSchemaCatalog::from_bytes(&trailing),
        Err(ReceiptSchemaCatalogError::Malformed { .. })
    ));

    assert!(matches!(
        ReceiptSchemaCatalog::from_bytes_verified(&bytes, fingerprint(0xa5)),
        Err(ReceiptSchemaCatalogError::ExternalIdentityMismatch { .. })
    ));

    let substituted = ReceiptSchemaCatalog::try_new(vec![ledger_checkpoint_v1()])
        .expect("self-consistent alternate catalog");
    let substituted_bytes = substituted.to_bytes();
    assert_eq!(
        ReceiptSchemaCatalog::from_bytes(&substituted_bytes).expect("self-consistent catalog"),
        substituted
    );
    assert!(matches!(
        ReceiptSchemaCatalog::from_bytes_verified(&substituted_bytes, canonical.content_hash()),
        Err(ReceiptSchemaCatalogError::ExternalIdentityMismatch { .. })
    ));

    let first = matdb_v2();
    let second = ledger_checkpoint_v1();
    assert!(second.family_id().as_bytes() < first.family_id().as_bytes());
    let mut reversed_wire = independent_catalog_wire(&[first, second]);
    let reversed_hash = independent_catalog_hash(&reversed_wire);
    reversed_wire.extend_from_slice(reversed_hash.as_bytes());
    assert!(matches!(
        ReceiptSchemaCatalog::from_bytes(&reversed_wire),
        Err(ReceiptSchemaCatalogError::NonCanonicalOrder { .. })
    ));

    let over_envelope = vec![0_u8; MAX_RECEIPT_SCHEMA_CATALOG_BYTES + 1];
    assert!(matches!(
        ReceiptSchemaCatalog::from_bytes(&over_envelope),
        Err(ReceiptSchemaCatalogError::ResourceLimit {
            resource: "catalog-bytes",
            ..
        })
    ));

    let mut insufficient_rows = bytes.clone();
    insufficient_rows[12..20].copy_from_slice(
        &u64::try_from(MAX_RECEIPT_SCHEMA_ENTRIES)
            .expect("entry cap fits u64")
            .to_le_bytes(),
    );
    assert!(matches!(
        ReceiptSchemaCatalog::from_bytes(&insufficient_rows),
        Err(ReceiptSchemaCatalogError::Malformed { .. })
    ));

    let mut zero_canonical_limit = bytes.clone();
    zero_canonical_limit[transport_tag + 1..transport_tag + 9]
        .copy_from_slice(&0_u64.to_le_bytes());
    assert!(matches!(
        ReceiptSchemaCatalog::from_bytes(&zero_canonical_limit),
        Err(ReceiptSchemaCatalogError::InvalidField {
            field: "maximum-transport-bytes",
            ..
        })
    ));

    let mut excessive_canonical_limit = bytes;
    excessive_canonical_limit[transport_tag + 1..transport_tag + 9]
        .copy_from_slice(&(MAX_RECEIPT_TRANSPORT_BYTES + 1).to_le_bytes());
    assert!(matches!(
        ReceiptSchemaCatalog::from_bytes(&excessive_canonical_limit),
        Err(ReceiptSchemaCatalogError::ResourceLimit {
            resource: "maximum-transport-bytes",
            ..
        })
    ));
}

#[test]
#[allow(clippy::too_many_lines)]
fn receipt_schema_catalog_refuses_invalid_identities_limits_and_aliases() {
    for family in [
        "",
        " fs-matdb:receipt",
        "fs-matdb:receipt ",
        "FS-matdb:receipt",
        "unqualified",
        "fs-matdb:todo",
        "fs-matdb::receipt",
        "fs-matdb:récépissé",
    ] {
        assert!(matches!(
            ReceiptSchemaDescriptor::try_new(
                family,
                1,
                1,
                "org.frankensim.fixture.receipt.v1",
                ReceiptTransportProfile::DigestOnly,
                fingerprint(1),
            ),
            Err(ReceiptSchemaCatalogError::InvalidField {
                field: "family-id",
                ..
            })
        ));
    }

    for domain in [
        "",
        " org.frankensim.receipt.v1",
        "org.frankensim.receipt.v1 ",
        "Org.frankensim.receipt.v1",
        "unqualified",
        "org.frankensim.pending.v1",
        "org..frankensim.receipt.v1",
        ".org.frankensim.receipt.v1",
    ] {
        assert!(matches!(
            ReceiptSchemaDescriptor::try_new(
                "fs-fixture:receipt",
                1,
                1,
                domain,
                ReceiptTransportProfile::DigestOnly,
                fingerprint(1),
            ),
            Err(ReceiptSchemaCatalogError::InvalidField {
                field: "owner-identity-domain",
                ..
            })
        ));
    }

    assert!(matches!(
        ReceiptSchemaDescriptor::try_new(
            "fs-fixture:receipt",
            0,
            1,
            "org.frankensim.fixture.receipt.v1",
            ReceiptTransportProfile::DigestOnly,
            fingerprint(1),
        ),
        Err(ReceiptSchemaCatalogError::InvalidField {
            field: "wire-schema-version",
            ..
        })
    ));
    assert!(matches!(
        ReceiptSchemaDescriptor::try_new(
            "fs-fixture:receipt",
            1,
            0,
            "org.frankensim.fixture.receipt.v1",
            ReceiptTransportProfile::DigestOnly,
            fingerprint(1),
        ),
        Err(ReceiptSchemaCatalogError::InvalidField {
            field: "owner-identity-version",
            ..
        })
    ));
    assert!(matches!(
        ReceiptSchemaDescriptor::try_new(
            "fs-fixture:receipt",
            1,
            1,
            "org.frankensim.fixture.receipt.v1",
            ReceiptTransportProfile::DigestOnly,
            ContentHash([0; 32]),
        ),
        Err(ReceiptSchemaCatalogError::InvalidField {
            field: "owner-schema-fingerprint",
            ..
        })
    ));

    for maximum_bytes in [0, MAX_RECEIPT_TRANSPORT_BYTES + 1] {
        assert!(
            ReceiptSchemaDescriptor::try_new(
                "fs-fixture:receipt",
                1,
                1,
                "org.frankensim.fixture.receipt.v1",
                ReceiptTransportProfile::CanonicalBytes { maximum_bytes },
                fingerprint(1),
            )
            .is_err()
        );
    }

    let maximum_family = format!(
        "fs:{}",
        "a".repeat(MAX_RECEIPT_FAMILY_ID_BYTES - "fs:".len())
    );
    let maximum_domain = format!(
        "d.{}",
        "a".repeat(MAX_RECEIPT_IDENTITY_DOMAIN_BYTES - "d.".len())
    );
    let boundary = ReceiptSchemaDescriptor::try_new(
        maximum_family,
        1,
        1,
        maximum_domain,
        ReceiptTransportProfile::CanonicalBytes {
            maximum_bytes: MAX_RECEIPT_TRANSPORT_BYTES,
        },
        fingerprint(1),
    )
    .expect("inclusive descriptor caps");
    assert_eq!(
        boundary.transport().maximum_bytes(),
        Some(MAX_RECEIPT_TRANSPORT_BYTES)
    );

    let oversized_family = format!(
        "fs:{}",
        "a".repeat(MAX_RECEIPT_FAMILY_ID_BYTES - "fs:".len() + 1)
    );
    assert!(matches!(
        ReceiptSchemaDescriptor::try_new(
            oversized_family,
            1,
            1,
            "org.frankensim.fixture.receipt.v1",
            ReceiptTransportProfile::DigestOnly,
            fingerprint(1),
        ),
        Err(ReceiptSchemaCatalogError::ResourceLimit {
            resource: "family-id",
            ..
        })
    ));
    let oversized_domain = format!(
        "d.{}",
        "a".repeat(MAX_RECEIPT_IDENTITY_DOMAIN_BYTES - "d.".len() + 1)
    );
    assert!(matches!(
        ReceiptSchemaDescriptor::try_new(
            "fs-fixture:receipt",
            1,
            1,
            oversized_domain,
            ReceiptTransportProfile::DigestOnly,
            fingerprint(1),
        ),
        Err(ReceiptSchemaCatalogError::ResourceLimit {
            resource: "owner-identity-domain",
            ..
        })
    ));

    let base = matdb_v2();
    let oversized_lookup = format!(
        "fs:{}",
        "a".repeat(MAX_RECEIPT_FAMILY_ID_BYTES - "fs:".len() + 1)
    );
    let one = ReceiptSchemaCatalog::try_new(vec![base.clone()]).expect("lookup catalog");
    assert!(matches!(
        one.require_exact(&oversized_lookup, 1, fingerprint(1)),
        Err(ReceiptSchemaCatalogError::ResourceLimit {
            resource: "family-id",
            ..
        })
    ));
    assert!(matches!(
        ReceiptSchemaCatalog::try_new(vec![base.clone(), base.clone()]),
        Err(ReceiptSchemaCatalogError::DuplicateSchema { .. })
    ));
    let reused_domain = descriptor(
        "fs-other:receipt",
        1,
        1,
        base.owner_identity_domain(),
        ReceiptTransportProfile::DigestOnly,
        0x71,
    );
    assert!(matches!(
        ReceiptSchemaCatalog::try_new(vec![base.clone(), reused_domain]),
        Err(ReceiptSchemaCatalogError::ReusedOwnerIdentityDomain { .. })
    ));
    let reused_fingerprint = ReceiptSchemaDescriptor::try_new(
        "fs-other:receipt",
        1,
        1,
        "org.frankensim.fs-other.receipt.v1",
        ReceiptTransportProfile::DigestOnly,
        base.owner_schema_fingerprint(),
    )
    .expect("individually valid alias fixture");
    assert!(matches!(
        ReceiptSchemaCatalog::try_new(vec![base.clone(), reused_fingerprint]),
        Err(ReceiptSchemaCatalogError::ReusedOwnerSchemaFingerprint { .. })
    ));

    assert!(matches!(
        ReceiptSchemaCatalog::try_new(vec![base; MAX_RECEIPT_SCHEMA_ENTRIES + 1]),
        Err(ReceiptSchemaCatalogError::ResourceLimit {
            resource: "schema-entries",
            ..
        })
    ));
}
