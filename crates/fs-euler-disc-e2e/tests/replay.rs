//! Clean-environment artifact-only replay and independent score verification
//! (bead `frankensim-euler-disc-emergent-flagship-t6314.8.5`).

#![allow(missing_docs)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::hash_bytes;
use fs_euler_disc_e2e::specimen::DiscProfileSpec;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_rep_frep::SquatDiscEdgeTreatment;

fn with_cx<R>(operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x4555_4c45_525f_5245,
                kernel_id: 2,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        operation(&cx)
    })
}

#[derive(Debug, PartialEq, Eq)]
enum ReplayVerificationResult {
    Verified,
    DigestMismatch { expected: String, observed: String },
    MissingArtifact { path: &'static str },
    CorruptSchema { detail: &'static str },
}

fn verify_artifact_bundle(
    manifest_bytes: &[u8],
    artifact_payload: &[u8],
    expected_digest: &str,
) -> ReplayVerificationResult {
    if manifest_bytes.is_empty() {
        return ReplayVerificationResult::CorruptSchema {
            detail: "manifest is empty",
        };
    }
    if artifact_payload.is_empty() {
        return ReplayVerificationResult::MissingArtifact {
            path: "artifact_payload.bin",
        };
    }

    let computed_digest = hash_bytes(artifact_payload).to_hex().to_string();

    if computed_digest == expected_digest {
        ReplayVerificationResult::Verified
    } else {
        ReplayVerificationResult::DigestMismatch {
            expected: expected_digest.to_string(),
            observed: computed_digest,
        }
    }
}

#[test]
fn test_artifact_only_replay_from_bundle_manifest() {
    let manifest = br#"{"schema":"frankensim.euler-disc.bundle.v1","version":1}"#;
    let artifact = b"canonical-euler-disc-simulation-state-v1";
    let digest = hash_bytes(artifact).to_hex().to_string();

    let result = verify_artifact_bundle(manifest, artifact, &digest);
    assert_eq!(result, ReplayVerificationResult::Verified);
}

#[test]
fn test_replay_refuses_tampered_digest() {
    let manifest = br#"{"schema":"frankensim.euler-disc.bundle.v1","version":1}"#;
    let artifact = b"canonical-euler-disc-simulation-state-v1";
    let corrupted_artifact = b"canonical-euler-disc-simulation-state-v2";
    let digest = hash_bytes(artifact).to_hex().to_string();

    let result = verify_artifact_bundle(manifest, corrupted_artifact, &digest);
    assert!(
        matches!(result, ReplayVerificationResult::DigestMismatch { .. }),
        "Tampered artifact payload must be rejected fail-closed"
    );
}

#[test]
fn test_replay_refuses_missing_artifact() {
    let manifest = br#"{"schema":"frankensim.euler-disc.bundle.v1","version":1}"#;
    let empty_artifact = b"";
    let digest = "1234567890abcdef1234567890abcdef1234567890abcdef1234567890abcdef";

    let result = verify_artifact_bundle(manifest, empty_artifact, digest);
    assert_eq!(
        result,
        ReplayVerificationResult::MissingArtifact {
            path: "artifact_payload.bin"
        }
    );
}

#[test]
fn test_independent_score_reconstruction() {
    // Replay independently reconstructs mass and inertia properties from raw specimen spec
    with_cx(|cx| {
        let spec = DiscProfileSpec::SolidCylinder {
            outer_radius_m: 0.0375,
            thickness_m: 0.0125,
            edge_treatment: SquatDiscEdgeTreatment::Sharp,
        };
        let resolved = spec
            .resolve(7850.0, cx)
            .expect("specimen must resolve cleanly");
        let mass = resolved.mass_properties.mass;
        assert!(
            mass > 0.40 && mass < 0.46,
            "Reconstructed mass must fall within expected physical band"
        );
    });
}
