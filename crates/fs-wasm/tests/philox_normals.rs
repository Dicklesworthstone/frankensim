//! Envelope tests for the fs-wasm `philox_normals` re-export.

use fs_wasm::{philox_normals, philox_normals_envelope_json};

#[test]
fn envelope_ok_is_not_empty_and_names_the_export() {
    let json = philox_normals_envelope_json(1, 0, 0, 0, 4);
    assert!(json.starts_with("{\"ok\":"), "{json}");
    assert!(json.contains("\"export\":\"philox_normals\""));
    let values = philox_normals(1, 0, 0, 0, 4).expect("admitted");
    assert_eq!(values.len(), 4);
}

#[test]
fn envelope_zero_count_is_refusal_not_empty_ok() {
    let json = philox_normals_envelope_json(1, 0, 0, 0, 0);
    assert!(json.contains("\"refusal\""), "{json}");
    assert!(json.contains("invalid-parameter"), "{json}");
    assert!(philox_normals(1, 0, 0, 0, 0).is_err());
}

#[test]
fn envelope_overflow_is_stream_index_overflow() {
    let json = philox_normals_envelope_json(0, 0, 0, 18_446_744_073_709_551_614, 1);
    assert!(json.contains("stream-index-overflow"), "{json}");
}
