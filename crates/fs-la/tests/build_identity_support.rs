//! G5 canonical-input tests for the GEMM build fingerprint.

#[path = "../build_identity_support.rs"]
mod support;

use support::{
    ASUPERSYNC_NON_SRC_INPUTS, EXECUTABLE_CONTEXT, FINGERPRINT_CONTEXT,
    append_executable_identity, append_external_identity, append_source_fields, push_field,
    push_optional_field,
};

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

    let mut changed = Vec::new();
    append_source_fields(
        &mut changed,
        vec![
            ("crates/fs-la/src/gemm.rs".to_string(), b"alphb".to_vec()),
            ("crates/fs-exec/src/pool.rs".to_string(), b"beta".to_vec()),
        ],
    );
    assert_ne!(forward, changed, "one source byte must change the payload");
}

#[test]
fn field_framing_separates_names_and_values() {
    assert!(!FINGERPRINT_CONTEXT.is_empty());
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
}
