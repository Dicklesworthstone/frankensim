//! G0 identity and refusal battery for resumable long-job request envelopes.

use fs_blake3::{ContentHash, hash_domain};
use fs_package::{
    ReceiptSchemaCatalog, ReceiptSchemaCatalogError, ReceiptSchemaDescriptor,
    ReceiptTransportProfile,
};
use fs_session::long_job::{
    LONG_JOB_REQUEST_IDENTITY_DOMAIN, LONG_JOB_REQUEST_IDENTITY_VERSION, LongJobBudget,
    LongJobKind, LongJobRequest, LongJobRequestError, ResumableModelIdentity,
};

fn hash(byte: u8) -> ContentHash {
    ContentHash([byte; 32])
}

fn descriptor(
    family: &str,
    wire_schema_version: u32,
    owner_identity_domain: &str,
    transport: ReceiptTransportProfile,
    fingerprint: u8,
) -> ReceiptSchemaDescriptor {
    ReceiptSchemaDescriptor::try_new(
        family,
        wire_schema_version,
        wire_schema_version,
        owner_identity_domain,
        transport,
        hash(fingerprint),
    )
    .expect("valid descriptor fixture")
}

fn baseline_descriptor() -> ReceiptSchemaDescriptor {
    descriptor(
        "fs-ir:hybrid-machine-checkpoint",
        1,
        "org.frankensim.fs-ir.hybrid-machine-checkpoint.v1",
        ReceiptTransportProfile::CanonicalBytes {
            maximum_bytes: 1_048_576,
        },
        0x31,
    )
}

#[derive(Clone)]
struct RequestSpec {
    kind: LongJobKind,
    operator: String,
    core_nanoseconds: u64,
    memory_bytes: u64,
    wall_nanoseconds: u64,
    max_parallel_cores: u64,
    canonical_program_hash: ContentHash,
    model_family: String,
    model_version: u32,
    state_schema_version: u32,
    model_instance_hash: ContentHash,
    contract_hash: ContentHash,
    code_hash: ContentHash,
    resume_descriptor: ReceiptSchemaDescriptor,
    additional_descriptors: Vec<ReceiptSchemaDescriptor>,
}

fn baseline_spec() -> RequestSpec {
    RequestSpec {
        kind: LongJobKind::HybridMachine,
        operator: "machine.run".to_string(),
        core_nanoseconds: 12_000_000_000,
        memory_bytes: 8 * 1024 * 1024 * 1024,
        wall_nanoseconds: 7_000_000_000,
        max_parallel_cores: 4,
        canonical_program_hash: hash(0x41),
        model_family: "machine.hybrid".to_string(),
        model_version: 7,
        state_schema_version: 3,
        model_instance_hash: hash(0x42),
        contract_hash: hash(0x43),
        code_hash: hash(0x44),
        resume_descriptor: baseline_descriptor(),
        additional_descriptors: Vec::new(),
    }
}

fn build(spec: &RequestSpec) -> LongJobRequest {
    let mut descriptors = vec![spec.resume_descriptor.clone()];
    descriptors.extend(spec.additional_descriptors.clone());
    let catalog = ReceiptSchemaCatalog::try_new(descriptors).expect("catalog fixture");
    let model = ResumableModelIdentity::try_new(
        &spec.model_family,
        spec.model_version,
        spec.state_schema_version,
        spec.model_instance_hash,
        spec.contract_hash,
        spec.code_hash,
        &spec.resume_descriptor,
    )
    .expect("model fixture");
    LongJobRequest::try_new(
        spec.kind,
        &spec.operator,
        LongJobBudget::try_new(
            spec.core_nanoseconds,
            spec.memory_bytes,
            spec.wall_nanoseconds,
            spec.max_parallel_cores,
        )
        .expect("budget fixture"),
        spec.canonical_program_hash,
        model,
        &catalog.to_bytes(),
        catalog.content_hash(),
    )
    .expect("request fixture")
}

fn append_u8(bytes: &mut Vec<u8>, value: u8) {
    bytes.push(value);
}

fn append_u32(bytes: &mut Vec<u8>, value: u32) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn append_u64(bytes: &mut Vec<u8>, value: u64) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn append_string(bytes: &mut Vec<u8>, value: &str) {
    append_u64(bytes, value.len() as u64);
    bytes.extend_from_slice(value.as_bytes());
}

fn append_hash(bytes: &mut Vec<u8>, value: ContentHash) {
    bytes.extend_from_slice(value.as_bytes());
}

fn independent_identity(request: &LongJobRequest) -> ContentHash {
    let mut fields = Vec::new();
    append_u8(&mut fields, 1);
    append_u8(
        &mut fields,
        match request.kind() {
            LongJobKind::HybridMachine => 1,
            LongJobKind::TheoremCheck => 2,
        },
    );
    append_u8(&mut fields, 2);
    append_string(&mut fields, request.operator());
    append_u8(&mut fields, 3);
    append_u64(&mut fields, request.budget().core_nanoseconds());
    append_u8(&mut fields, 4);
    append_u64(&mut fields, request.budget().memory_bytes());
    append_u8(&mut fields, 5);
    append_u64(&mut fields, request.budget().wall_nanoseconds());
    append_u8(&mut fields, 6);
    append_u64(&mut fields, request.budget().max_parallel_cores());
    append_u8(&mut fields, 7);
    append_hash(&mut fields, request.canonical_program_hash());
    append_u8(&mut fields, 8);
    append_string(&mut fields, request.model().family());
    append_u8(&mut fields, 9);
    append_u32(&mut fields, request.model().model_version());
    append_u8(&mut fields, 10);
    append_u32(&mut fields, request.model().state_schema_version());
    append_u8(&mut fields, 11);
    append_hash(&mut fields, request.model().model_instance_hash());
    append_u8(&mut fields, 12);
    append_hash(&mut fields, request.model().contract_hash());
    append_u8(&mut fields, 13);
    append_hash(&mut fields, request.model().code_hash());
    append_u8(&mut fields, 14);
    append_hash(&mut fields, request.receipt_catalog_hash());
    append_u8(&mut fields, 15);
    append_string(&mut fields, request.model().resume_schema().family_id());
    append_u8(&mut fields, 16);
    append_u32(
        &mut fields,
        request.model().resume_schema().wire_schema_version(),
    );
    append_u8(&mut fields, 17);
    append_hash(
        &mut fields,
        request.model().resume_schema().descriptor_hash(),
    );

    let mut preimage = Vec::with_capacity(4 + fields.len());
    preimage.extend_from_slice(&LONG_JOB_REQUEST_IDENTITY_VERSION.to_le_bytes());
    preimage.extend_from_slice(&fields);
    hash_domain(LONG_JOB_REQUEST_IDENTITY_DOMAIN, &preimage)
}

#[test]
fn long_job_request_identity_matches_independent_preimage() {
    let request = build(&baseline_spec());
    assert_eq!(request.content_hash(), independent_identity(&request));
    assert_eq!(request, build(&baseline_spec()));
    assert_eq!(request.kind(), LongJobKind::HybridMachine);
    assert_eq!(request.operator(), "machine.run");
    assert_eq!(request.model().family(), "machine.hybrid");
    assert_eq!(
        request.model().resume_schema().family_id(),
        "fs-ir:hybrid-machine-checkpoint"
    );
}

fn assert_moves(label: &str, spec: &RequestSpec, baseline: ContentHash) {
    let moved = build(spec).content_hash();
    assert_ne!(baseline, moved, "{label} must move request identity");
}

#[test]
fn long_job_kind_and_every_field_move_identity() {
    let baseline = build(&baseline_spec()).content_hash();

    let mut spec = baseline_spec();
    spec.kind = LongJobKind::TheoremCheck;
    assert_moves("job kind", &spec, baseline);

    let mut spec = baseline_spec();
    spec.operator = "theorem.run".to_string();
    assert_moves("operator bytes", &spec, baseline);
    let mut spec = baseline_spec();
    spec.operator = "machine.execute".to_string();
    assert_moves("operator byte count", &spec, baseline);

    let mut spec = baseline_spec();
    spec.core_nanoseconds += 1;
    assert_moves("core nanoseconds", &spec, baseline);
    let mut spec = baseline_spec();
    spec.memory_bytes += 1;
    assert_moves("memory bytes", &spec, baseline);
    let mut spec = baseline_spec();
    spec.wall_nanoseconds += 1;
    assert_moves("wall nanoseconds", &spec, baseline);
    let mut spec = baseline_spec();
    spec.max_parallel_cores += 1;
    assert_moves("parallel cores", &spec, baseline);

    let mut spec = baseline_spec();
    spec.canonical_program_hash = hash(0x51);
    assert_moves("program hash", &spec, baseline);

    let mut spec = baseline_spec();
    spec.model_family = "theorem.check".to_string();
    assert_moves("model family bytes", &spec, baseline);
    let mut spec = baseline_spec();
    spec.model_family = "machine.hybrid.extended".to_string();
    assert_moves("model family byte count", &spec, baseline);
    let mut spec = baseline_spec();
    spec.model_version += 1;
    assert_moves("model version", &spec, baseline);
    let mut spec = baseline_spec();
    spec.state_schema_version += 1;
    assert_moves("state schema version", &spec, baseline);
    let mut spec = baseline_spec();
    spec.model_instance_hash = hash(0x52);
    assert_moves("model instance", &spec, baseline);
    let mut spec = baseline_spec();
    spec.contract_hash = hash(0x53);
    assert_moves("contract", &spec, baseline);
    let mut spec = baseline_spec();
    spec.code_hash = hash(0x54);
    assert_moves("code", &spec, baseline);

    let mut spec = baseline_spec();
    spec.additional_descriptors.push(descriptor(
        "fs-ir:theorem-frontier",
        1,
        "org.frankensim.fs-ir.theorem-frontier.v1",
        ReceiptTransportProfile::CanonicalBytes {
            maximum_bytes: 2_097_152,
        },
        0x61,
    ));
    assert_moves("catalog hash", &spec, baseline);

    let mut spec = baseline_spec();
    spec.resume_descriptor = descriptor(
        "fs-ir:hybrid-machine-checkpoint-next",
        1,
        "org.frankensim.fs-ir.hybrid-machine-checkpoint-next.v1",
        ReceiptTransportProfile::CanonicalBytes {
            maximum_bytes: 1_048_576,
        },
        0x62,
    );
    assert_moves("resume family", &spec, baseline);

    let mut spec = baseline_spec();
    spec.resume_descriptor = descriptor(
        "fs-ir:hybrid-machine-checkpoint",
        2,
        "org.frankensim.fs-ir.hybrid-machine-checkpoint.v2",
        ReceiptTransportProfile::CanonicalBytes {
            maximum_bytes: 1_048_576,
        },
        0x63,
    );
    assert_moves("resume wire schema", &spec, baseline);

    let mut spec = baseline_spec();
    spec.resume_descriptor = descriptor(
        "fs-ir:hybrid-machine-checkpoint",
        1,
        "org.frankensim.fs-ir.hybrid-machine-checkpoint-alternate.v1",
        ReceiptTransportProfile::CanonicalBytes {
            maximum_bytes: 1_048_576,
        },
        0x64,
    );
    assert_moves("resume descriptor", &spec, baseline);
}

#[test]
fn long_job_request_refuses_invalid_numeric_and_name_domains() {
    for invalid in [
        LongJobBudget::try_new(0, 2, 3, 4),
        LongJobBudget::try_new(1, 0, 3, 4),
        LongJobBudget::try_new(1, 2, 0, 4),
        LongJobBudget::try_new(1, 2, 3, 0),
    ] {
        assert!(matches!(
            invalid,
            Err(LongJobRequestError::InvalidField { .. })
        ));
    }

    let resume = baseline_descriptor();
    for family in [
        "",
        "hybrid",
        "machine.*",
        ".machine",
        "Machine.hybrid",
        "machine..hybrid",
        "machine.unknown",
    ] {
        assert!(matches!(
            ResumableModelIdentity::try_new(family, 1, 1, hash(1), hash(2), hash(3), &resume,),
            Err(LongJobRequestError::InvalidField {
                field: "model-family",
                ..
            })
        ));
    }
    let oversized_family = format!("machine.{}", "x".repeat(129));
    assert!(matches!(
        ResumableModelIdentity::try_new(
            &oversized_family,
            1,
            1,
            hash(1),
            hash(2),
            hash(3),
            &resume,
        ),
        Err(LongJobRequestError::ResourceLimit {
            resource: "model-family",
            ..
        })
    ));

    for (model_version, state_schema_version, instance, contract, code) in [
        (0, 1, hash(1), hash(2), hash(3)),
        (1, 0, hash(1), hash(2), hash(3)),
        (1, 1, hash(0), hash(2), hash(3)),
        (1, 1, hash(1), hash(0), hash(3)),
        (1, 1, hash(1), hash(2), hash(0)),
    ] {
        assert!(matches!(
            ResumableModelIdentity::try_new(
                "machine.hybrid",
                model_version,
                state_schema_version,
                instance,
                contract,
                code,
                &resume,
            ),
            Err(LongJobRequestError::InvalidField { .. })
        ));
    }

    let catalog = ReceiptSchemaCatalog::try_new(vec![resume.clone()]).expect("catalog");
    let model =
        ResumableModelIdentity::try_new("machine.hybrid", 1, 1, hash(1), hash(2), hash(3), &resume)
            .expect("model");
    let budget = LongJobBudget::try_new(1, 2, 3, 4).expect("budget");
    for operator in ["", "machine.*", ".machine", "machine run"] {
        assert!(matches!(
            LongJobRequest::try_new(
                LongJobKind::HybridMachine,
                operator,
                budget,
                hash(4),
                model.clone(),
                &catalog.to_bytes(),
                catalog.content_hash(),
            ),
            Err(LongJobRequestError::InvalidField {
                field: "operator",
                ..
            })
        ));
    }
    let oversized_operator = "x".repeat(129);
    assert!(matches!(
        LongJobRequest::try_new(
            LongJobKind::HybridMachine,
            &oversized_operator,
            budget,
            hash(4),
            model.clone(),
            &catalog.to_bytes(),
            catalog.content_hash(),
        ),
        Err(LongJobRequestError::ResourceLimit {
            resource: "operator",
            ..
        })
    ));
    assert!(matches!(
        LongJobRequest::try_new(
            LongJobKind::HybridMachine,
            "machine.run",
            budget,
            hash(0),
            model.clone(),
            &catalog.to_bytes(),
            catalog.content_hash(),
        ),
        Err(LongJobRequestError::InvalidField {
            field: "canonical-program-hash",
            ..
        })
    ));
    assert!(matches!(
        LongJobRequest::try_new(
            LongJobKind::HybridMachine,
            "machine.run",
            budget,
            hash(4),
            model.clone(),
            &catalog.to_bytes(),
            hash(0),
        ),
        Err(LongJobRequestError::InvalidField {
            field: "receipt-catalog-hash",
            ..
        })
    ));
    let substituted_descriptor = descriptor(
        "fs-ir:theorem-frontier",
        1,
        "org.frankensim.fs-ir.theorem-frontier.v1",
        ReceiptTransportProfile::CanonicalBytes { maximum_bytes: 1 },
        0x75,
    );
    let substituted_catalog =
        ReceiptSchemaCatalog::try_new(vec![resume, substituted_descriptor]).expect("substitute");
    assert!(matches!(
        LongJobRequest::try_new(
            LongJobKind::HybridMachine,
            "machine.run",
            budget,
            hash(4),
            model,
            &substituted_catalog.to_bytes(),
            catalog.content_hash(),
        ),
        Err(LongJobRequestError::ReceiptCatalog {
            source: ReceiptSchemaCatalogError::ExternalIdentityMismatch { .. },
            ..
        })
    ));
}

#[test]
fn catalog_and_resume_schema_are_exactly_validated() {
    let cataloged = baseline_descriptor();
    let catalog = ReceiptSchemaCatalog::try_new(vec![cataloged.clone()]).expect("catalog");
    let budget = LongJobBudget::try_new(1, 2_000_000, 3, 4).expect("budget");

    let absent_family = descriptor(
        "fs-ir:absent-checkpoint-family",
        1,
        "org.frankensim.fs-ir.absent-checkpoint-family.v1",
        ReceiptTransportProfile::CanonicalBytes { maximum_bytes: 64 },
        0x71,
    );
    let absent_family_model = ResumableModelIdentity::try_new(
        "machine.hybrid",
        1,
        1,
        hash(1),
        hash(2),
        hash(3),
        &absent_family,
    )
    .expect("model");
    let expected_catalog_hash = catalog.content_hash();
    let expected_descriptor_hash = absent_family.content_hash();
    let error = LongJobRequest::try_new(
        LongJobKind::HybridMachine,
        "machine.run",
        budget,
        hash(4),
        absent_family_model,
        &catalog.to_bytes(),
        expected_catalog_hash,
    )
    .expect_err("absent family must refuse");
    let message = error.to_string();
    assert!(message.contains(&expected_catalog_hash.to_hex()));
    assert!(message.contains("model machine.hybrid state schema v1"));
    assert!(message.contains("resume row fs-ir:absent-checkpoint-family wire v1"));
    assert!(message.contains(&expected_descriptor_hash.to_hex()));
    assert!(matches!(
        error,
        LongJobRequestError::ReceiptCatalog {
            catalog_hash,
            ref model_family,
            state_schema_version: 1,
            ref resume_family,
            resume_wire_schema_version: 1,
            resume_descriptor_hash,
            source: ReceiptSchemaCatalogError::UnknownFamily { .. },
        } if catalog_hash == expected_catalog_hash
            && model_family == "machine.hybrid"
            && resume_family == "fs-ir:absent-checkpoint-family"
            && resume_descriptor_hash == expected_descriptor_hash
    ));

    let absent_version = descriptor(
        cataloged.family_id(),
        cataloged.wire_schema_version() + 1,
        "org.frankensim.fs-ir.hybrid-machine-checkpoint.v2",
        ReceiptTransportProfile::CanonicalBytes { maximum_bytes: 64 },
        0x72,
    );
    let absent_version_model = ResumableModelIdentity::try_new(
        "machine.hybrid",
        1,
        2,
        hash(1),
        hash(2),
        hash(3),
        &absent_version,
    )
    .expect("model");
    assert!(matches!(
        LongJobRequest::try_new(
            LongJobKind::HybridMachine,
            "machine.run",
            budget,
            hash(4),
            absent_version_model,
            &catalog.to_bytes(),
            catalog.content_hash(),
        ),
        Err(LongJobRequestError::ReceiptCatalog {
            source: ReceiptSchemaCatalogError::UnsupportedWireSchema { .. },
            ..
        })
    ));

    let mismatched = descriptor(
        cataloged.family_id(),
        cataloged.wire_schema_version(),
        "org.frankensim.fs-ir.hybrid-machine-checkpoint-alternate.v1",
        ReceiptTransportProfile::CanonicalBytes { maximum_bytes: 64 },
        0x73,
    );
    let mismatched_model = ResumableModelIdentity::try_new(
        "machine.hybrid",
        1,
        1,
        hash(1),
        hash(2),
        hash(3),
        &mismatched,
    )
    .expect("model");
    assert!(matches!(
        LongJobRequest::try_new(
            LongJobKind::HybridMachine,
            "machine.run",
            budget,
            hash(4),
            mismatched_model,
            &catalog.to_bytes(),
            catalog.content_hash(),
        ),
        Err(LongJobRequestError::ReceiptCatalog {
            source: ReceiptSchemaCatalogError::DescriptorMismatch { .. },
            ..
        })
    ));

    let digest_only = descriptor(
        "fs-ir:theorem-frontier",
        1,
        "org.frankensim.fs-ir.theorem-frontier.v1",
        ReceiptTransportProfile::DigestOnly,
        0x74,
    );
    let digest_catalog =
        ReceiptSchemaCatalog::try_new(vec![digest_only.clone()]).expect("digest catalog");
    let digest_model = ResumableModelIdentity::try_new(
        "theorem.check",
        1,
        1,
        hash(1),
        hash(2),
        hash(3),
        &digest_only,
    )
    .expect("model");
    assert!(matches!(
        LongJobRequest::try_new(
            LongJobKind::TheoremCheck,
            "theorem.check",
            budget,
            hash(4),
            digest_model,
            &digest_catalog.to_bytes(),
            digest_catalog.content_hash(),
        ),
        Err(LongJobRequestError::ResumeTransportUnavailable { .. })
    ));

    let bounded = descriptor(
        "fs-ir:bounded-frontier",
        1,
        "org.frankensim.fs-ir.bounded-frontier.v1",
        ReceiptTransportProfile::CanonicalBytes { maximum_bytes: 65 },
        0x76,
    );
    let bounded_catalog =
        ReceiptSchemaCatalog::try_new(vec![bounded.clone()]).expect("bounded catalog");
    let bounded_model =
        ResumableModelIdentity::try_new("theorem.check", 1, 1, hash(1), hash(2), hash(3), &bounded)
            .expect("bounded model");
    assert!(matches!(
        LongJobRequest::try_new(
            LongJobKind::TheoremCheck,
            "theorem.check",
            LongJobBudget::try_new(1, 64, 3, 4).expect("small memory"),
            hash(4),
            bounded_model.clone(),
            &bounded_catalog.to_bytes(),
            bounded_catalog.content_hash(),
        ),
        Err(LongJobRequestError::ResumeTransportExceedsMemoryBudget {
            maximum_bytes: 65,
            requested_memory_bytes: 64,
            ..
        })
    ));
    let theorem_request = LongJobRequest::try_new(
        LongJobKind::TheoremCheck,
        "theorem.check",
        LongJobBudget::try_new(1, 65, 3, 4).expect("boundary memory"),
        hash(4),
        bounded_model,
        &bounded_catalog.to_bytes(),
        bounded_catalog.content_hash(),
    )
    .expect("accepted theorem request");
    assert_eq!(theorem_request.kind(), LongJobKind::TheoremCheck);

    assert!(
        LongJobRequest::try_new(
            LongJobKind::HybridMachine,
            "machine.run",
            budget,
            hash(4),
            ResumableModelIdentity::try_new(
                "machine.hybrid",
                1,
                1,
                hash(1),
                hash(2),
                hash(3),
                &cataloged,
            )
            .expect("cataloged model"),
            &catalog.to_bytes(),
            catalog.content_hash(),
        )
        .is_ok()
    );
}
