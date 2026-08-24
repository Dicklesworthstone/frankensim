//! Cinematic export battery (bead `frankensim-wf-root-guzez.11.8`, E10.4).

use fs_flyer::cinematic_export::{
    export_hero_clip, CinematicClipManifest, QUARANTINED_MUX_ADAPTER,
};
use fs_flyer::replay::{AppliedEvent, InputTrace};

fn test_trace() -> InputTrace {
    InputTrace {
        end_tick_exclusive: 120,
        events: vec![
            AppliedEvent {
                channel: 0,
                applied_tick: 10,
                ordinal_within_tick: 0,
                quantized_value: 0.1,
            },
            AppliedEvent {
                channel: 1,
                applied_tick: 50,
                ordinal_within_tick: 0,
                quantized_value: -0.2,
            },
        ],
    }
}

#[test]
fn hero_clip_export_binds_run_id_and_manifest() {
    let trace = test_trace();
    let run_id = "run-wright-1903-dec17-001";
    let receipt = export_hero_clip(run_id, &trace, 120).expect("hero clip export succeeds");

    assert_eq!(receipt.run_id, run_id);
    assert_eq!(receipt.frames_rendered, 120);
    assert_eq!(receipt.quarantined_adapter, QUARANTINED_MUX_ADAPTER);
    assert!(receipt.tamper_detection_verified);
    assert!(!receipt.manifest_digest.is_empty());
    assert!(!receipt.mux_receipt_id.is_empty());
}

#[test]
fn hostile_twin_tampered_run_id_refuses() {
    let manifest = CinematicClipManifest::new(
        "authoritative-run-001",
        "trace-120",
        60,
        60,
        [1920, 1080],
    )
    .expect("valid manifest");

    let err = manifest
        .verify_origin("tampered-run-999")
        .expect_err("must refuse tampered run_id");
    assert_eq!(err.code, "cinematic-run-id-tampered");
}

#[test]
fn hostile_twin_tampered_digest_refuses() {
    let mut manifest = CinematicClipManifest::new(
        "authoritative-run-001",
        "trace-120",
        60,
        60,
        [1920, 1080],
    )
    .expect("valid manifest");

    // Tamper frame count without recomputing digest
    manifest.frame_count = 120;

    let err = manifest
        .verify_origin("authoritative-run-001")
        .expect_err("must refuse tampered payload");
    assert_eq!(err.code, "cinematic-manifest-tampered");
}
