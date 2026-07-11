//! Pure canonical framing helpers shared by the build script and its tests.

pub const FINGERPRINT_CONTEXT: &str = "frankensim.fs-la.gemm-build-fingerprint.v1";

/// Domain for compiler and wrapper executable content identities.
pub const EXECUTABLE_CONTEXT: &str = "frankensim.fs-la.gemm-build-executable.v1";

/// Compiler inputs included by the normal asupersync graph from outside its
/// package source directories. Keep this explicit so each addition is audited
/// against a concrete `include_*` site rather than sweeping unrelated assets.
pub const ASUPERSYNC_NON_SRC_INPUTS: &[&str] = &["assets/dashboard.html"];

pub fn push_field(payload: &mut Vec<u8>, name: &str, value: &[u8]) {
    payload.extend_from_slice(&(name.len() as u64).to_le_bytes());
    payload.extend_from_slice(name.as_bytes());
    payload.extend_from_slice(&(value.len() as u64).to_le_bytes());
    payload.extend_from_slice(value);
}

pub fn push_optional_field(payload: &mut Vec<u8>, name: &str, value: Option<&[u8]>) {
    push_field(
        payload,
        &format!("{name}:presence"),
        if value.is_some() { b"present" } else { b"absent" },
    );
    if let Some(value) = value {
        push_field(payload, name, value);
    }
}

pub fn append_source_fields(payload: &mut Vec<u8>, mut fields: Vec<(String, Vec<u8>)>) {
    fields.sort_unstable_by(|left, right| left.0.cmp(&right.0));
    assert!(
        fields
            .windows(2)
            .all(|pair| pair[0].0.as_str() != pair[1].0.as_str()),
        "duplicate path in GEMM build-identity source closure"
    );
    for (relative, bytes) in fields {
        push_field(payload, &format!("source:{relative}"), &bytes);
    }
}

pub fn append_external_identity(
    payload: &mut Vec<u8>,
    constellation_lock: &[u8],
    git_head: &str,
    fields: Vec<(String, Vec<u8>)>,
) {
    push_field(payload, "constellation.lock", constellation_lock);
    push_field(payload, "asupersync-git-head", git_head.as_bytes());
    append_source_fields(payload, fields);
}

pub fn append_executable_identity(
    payload: &mut Vec<u8>,
    label: &str,
    resolved_path: &str,
    bytes: &[u8],
) {
    push_field(payload, &format!("executable:{label}:path"), resolved_path.as_bytes());
    let digest = fs_blake3::hash_domain(EXECUTABLE_CONTEXT, bytes);
    push_field(
        payload,
        &format!("executable:{label}:content"),
        digest.as_bytes(),
    );
}
