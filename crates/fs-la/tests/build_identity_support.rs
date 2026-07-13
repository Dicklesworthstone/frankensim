//! G5 canonical-input tests for the GEMM build fingerprint.

#[path = "../build_identity_support.rs"]
mod support;

use support::{
    ASUPERSYNC_NON_SRC_INPUTS, EXECUTABLE_CONTEXT, GEMM_BUILD_FINGERPRINT_IDENTITY_DOMAIN,
    GEMM_BUILD_FINGERPRINT_IDENTITY_VERSION, GEMM_BUILD_PAYLOAD_SCHEMA, GemmBuildIdentityInput,
    append_executable_identity, append_external_identity, append_source_fields,
    gemm_build_fingerprint, gemm_build_fingerprint_identity_version_is_supported,
    gemm_build_fingerprint_with_domain, push_field, push_optional_field,
};

fn fingerprint(payload: &[u8]) -> fs_blake3::ContentHash {
    gemm_build_fingerprint(&GemmBuildIdentityInput {
        canonical_payload: payload,
    })
}

#[test]
fn gemm_build_identity_domain_moves_fingerprint() {
    let input = GemmBuildIdentityInput {
        canonical_payload: b"canonical-build-payload",
    };
    assert_ne!(
        gemm_build_fingerprint(&input),
        gemm_build_fingerprint_with_domain(
            "frankensim.fs-la.gemm-build-fingerprint.v1.alternate",
            &input,
        )
    );
}

#[test]
fn embedded_payload_schema_moves_gemm_build_fingerprint() {
    let mut current = Vec::new();
    push_field(&mut current, "schema", GEMM_BUILD_PAYLOAD_SCHEMA.as_bytes());
    let mut changed = Vec::new();
    push_field(&mut changed, "schema", b"fs-la-gemm-codegen-v3");
    assert_ne!(fingerprint(&current), fingerprint(&changed));
}

#[test]
fn depgraph_evidence_moves_gemm_build_fingerprint() {
    let mut current = Vec::new();
    push_field(
        &mut current,
        "depgraph-receipt",
        br#"{"schema":"fs-la-depgraph-receipt-v1","packages":["a"]}"#,
    );
    let mut changed = Vec::new();
    push_field(
        &mut changed,
        "depgraph-receipt",
        br#"{"schema":"fs-la-depgraph-receipt-v1","packages":["b"]}"#,
    );
    assert_ne!(fingerprint(&current), fingerprint(&changed));
}

#[test]
fn gemm_build_fingerprint_identity_version_fails_closed() {
    assert_eq!(GEMM_BUILD_FINGERPRINT_IDENTITY_VERSION, 1);
    assert_eq!(
        GEMM_BUILD_FINGERPRINT_IDENTITY_DOMAIN,
        "frankensim.fs-la.gemm-build-fingerprint.v1"
    );
    assert_eq!(GEMM_BUILD_PAYLOAD_SCHEMA, "fs-la-gemm-codegen-v2");
    assert!(gemm_build_fingerprint_identity_version_is_supported(1));
    assert!(!gemm_build_fingerprint_identity_version_is_supported(0));
    assert!(!gemm_build_fingerprint_identity_version_is_supported(2));
}

#[test]
fn source_fields_are_order_independent_and_content_sensitive() {
    let fields = vec![
        ("crates/fs-la/src/gemm.rs".to_string(), b"alpha".to_vec()),
        ("crates/fs-exec/src/pool.rs".to_string(), b"beta".to_vec()),
    ];
    let mut forward = Vec::new();
    append_source_fields(&mut forward, fields.clone());
    let mut reverse = Vec::new();
    append_source_fields(&mut reverse, fields.into_iter().rev().collect());
    assert_eq!(forward, reverse, "directory order must be irrelevant");
    assert_eq!(
        fingerprint(&forward),
        fingerprint(&reverse),
        "directory order must not move the outer build fingerprint"
    );

    let mut changed = Vec::new();
    append_source_fields(
        &mut changed,
        vec![
            ("crates/fs-la/src/gemm.rs".to_string(), b"alphb".to_vec()),
            ("crates/fs-exec/src/pool.rs".to_string(), b"beta".to_vec()),
        ],
    );
    assert_ne!(forward, changed, "one source byte must change the payload");
    assert_ne!(
        fingerprint(&forward),
        fingerprint(&changed),
        "one source byte must change the outer build fingerprint"
    );
}

#[test]
fn field_framing_separates_names_and_values() {
    assert!(!GEMM_BUILD_FINGERPRINT_IDENTITY_DOMAIN.is_empty());
    assert!(!EXECUTABLE_CONTEXT.is_empty());
    let mut left = Vec::new();
    push_field(&mut left, "ab", b"c");
    let mut right = Vec::new();
    push_field(&mut right, "a", b"bc");
    assert_ne!(left, right);
}

#[test]
fn executable_identity_binds_resolved_path_and_content() {
    fn payload(path: &str, bytes: &[u8]) -> Vec<u8> {
        let mut payload = Vec::new();
        append_executable_identity(&mut payload, "RUSTC", path, bytes);
        payload
    }

    let baseline = payload("/toolchain/bin/rustc", b"compiler-a");
    assert_ne!(baseline, payload("/other/bin/rustc", b"compiler-a"));
    assert_ne!(baseline, payload("/toolchain/bin/rustc", b"compiler-b"));
    assert_eq!(baseline, payload("/toolchain/bin/rustc", b"compiler-a"));
    assert_ne!(
        fingerprint(&baseline),
        fingerprint(&payload("/toolchain/bin/rustc", b"compiler-b"))
    );
}

#[test]
fn optional_fields_separate_absent_empty_and_literal_sentinel() {
    fn payload(value: Option<&[u8]>) -> Vec<u8> {
        let mut payload = Vec::new();
        push_optional_field(&mut payload, "SALT", value);
        payload
    }

    assert_ne!(payload(None), payload(Some(b"")));
    assert_ne!(payload(None), payload(Some(b"<unset>")));
    assert_ne!(payload(Some(b"")), payload(Some(b"<unset>")));
    assert_ne!(
        fingerprint(&payload(None)),
        fingerprint(&payload(Some(b"")))
    );
}

#[test]
fn external_identity_binds_lock_head_source_and_include_inputs() {
    fn payload(lock: &[u8], head: &str, source: &[u8], dashboard: &[u8]) -> Vec<u8> {
        let mut payload = Vec::new();
        append_external_identity(
            &mut payload,
            lock,
            head,
            vec![
                (
                    "external/asupersync/src/lib.rs".to_string(),
                    source.to_vec(),
                ),
                (
                    "external/asupersync/assets/dashboard.html".to_string(),
                    dashboard.to_vec(),
                ),
            ],
        );
        payload
    }

    assert_eq!(ASUPERSYNC_NON_SRC_INPUTS, &["assets/dashboard.html"]);
    let baseline = payload(
        b"lock-a",
        "1111111111111111111111111111111111111111",
        b"src-a",
        b"dashboard-a",
    );
    assert_ne!(
        baseline,
        payload(
            b"lock-b",
            "1111111111111111111111111111111111111111",
            b"src-a",
            b"dashboard-a"
        )
    );
    assert_ne!(
        baseline,
        payload(
            b"lock-a",
            "2222222222222222222222222222222222222222",
            b"src-a",
            b"dashboard-a"
        )
    );
    assert_ne!(
        baseline,
        payload(
            b"lock-a",
            "1111111111111111111111111111111111111111",
            b"src-b",
            b"dashboard-a"
        )
    );
    assert_ne!(
        baseline,
        payload(
            b"lock-a",
            "1111111111111111111111111111111111111111",
            b"src-a",
            b"dashboard-b"
        ),
        "an included non-source compiler input must change the payload"
    );
    assert_ne!(
        fingerprint(&baseline),
        fingerprint(&payload(
            b"lock-a",
            "1111111111111111111111111111111111111111",
            b"src-a",
            b"dashboard-b"
        )),
        "an included external input must change the outer build fingerprint"
    );
}
