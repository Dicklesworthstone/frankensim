//! Listening-receipt conformance + minting (music bead
//! `frankensim-music-v8-root-3ez8g.1.2`).
//!
//! The default test validates the COMMITTED receipt at
//! `data/listening/reed-sustained.listening-receipt` against the
//! COMMITTED artifact's provenance sidecar: the receipt must decode, its
//! artifact digest must equal the sidecar's `wav_blake3`, and its
//! attached log-attack-time must match a recomputation from the WAV
//! bytes — so the committed evidence chain (WAV → sidecar → receipt)
//! cannot drift apart. The `--ignored` minting test regenerates receipt
//! bytes for a fresh render.
//!
//! Calibration honesty, demonstrated for real: WITHOUT a `Calibration`
//! the only computable attached metric is level-relative log-attack-time.
//! Loudness/sharpness need absolute band levels and spl needs the
//! calibration bridge, so those fields are ABSENT in the committed
//! receipt — empty, never fabricated (the fs-psycho refusal made
//! visible in evidence bytes).

use fs_psycho::receipt::{ListeningReceipt, ListeningVerdict};

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("workspace root")
        .to_path_buf()
}

/// Decode the committed mono PCM16 WAV back to pascals.
fn wav_to_pascals(bytes: &[u8], full_scale_pa: f64) -> Vec<f64> {
    assert_eq!(&bytes[0..4], b"RIFF", "committed artifact must be a WAV");
    bytes[44..]
        .chunks_exact(2)
        .map(|c| f64::from(i16::from_le_bytes([c[0], c[1]])) / f64::from(i16::MAX) * full_scale_pa)
        .collect()
}

fn sidecar_field(sidecar: &str, key: &str) -> String {
    let tag = format!("\"{key}\":");
    let start = sidecar.find(&tag).expect("sidecar field") + tag.len();
    let rest = &sidecar[start..];
    let rest = rest.strip_prefix('"').unwrap_or(rest);
    let end = rest.find(['"', ',', '}']).expect("field end");
    rest[..end].to_string()
}

#[test]
fn committed_receipt_matches_the_committed_artifact() {
    let root = repo_root();
    let receipt_bytes = std::fs::read(root.join("data/listening/reed-sustained.listening-receipt"))
        .expect("committed receipt");
    let receipt = ListeningReceipt::from_canonical_bytes(&receipt_bytes).expect("receipt decodes");
    let sidecar =
        std::fs::read_to_string(root.join("data/listening/reed-sustained.provenance.json"))
            .expect("committed sidecar");
    // Chain link 1: receipt digest == sidecar digest.
    assert_eq!(
        receipt.artifact_hex,
        sidecar_field(&sidecar, "wav_blake3"),
        "receipt must reference exactly the committed artifact"
    );
    // Chain link 2: the attached log-attack-time recomputes from the WAV.
    let full_scale: f64 = sidecar_field(&sidecar, "full_scale_pa")
        .parse()
        .expect("full scale");
    let wav = std::fs::read(root.join("data/listening/reed-sustained.wav")).expect("wav");
    let pascals = wav_to_pascals(&wav, full_scale);
    let lat = fs_psycho::log_attack_time(&pascals, 48_000.0, 480).expect("attack time");
    let attached = receipt
        .metrics
        .log_attack_time
        .expect("committed receipt attaches log-attack-time");
    assert!(
        (lat - attached).abs() < 1.0e-9,
        "attached log-attack-time {attached} must recompute from the WAV ({lat})"
    );
    // Calibration honesty in the committed bytes: no calibration existed,
    // so the absolute fields are absent.
    assert!(receipt.metrics.spl_db.is_none(), "no calibration => no SPL");
    assert!(
        receipt.metrics.loudness_sone.is_none() && receipt.metrics.sharpness_acum.is_none(),
        "absolute-level metrics stay absent without a calibration"
    );
    // The unadjudicated state is honest and structurally useless for gates.
    assert_eq!(receipt.verdict, ListeningVerdict::Unadjudicated);
    assert!(!receipt.supports_pass());
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"listening-receipt-chain\",\"artifact\":\"{}\",\
         \"log_attack_time\":{lat:.6},\"verdict\":\"{}\"}}",
        receipt.artifact_hex,
        receipt.verdict.name()
    );
}

#[test]
#[ignore = "minting helper: regenerates receipt bytes for a fresh render; commit output manually"]
fn mint_receipt_bytes_for_the_committed_artifact() {
    let root = repo_root();
    let sidecar =
        std::fs::read_to_string(root.join("data/listening/reed-sustained.provenance.json"))
            .expect("sidecar");
    let full_scale: f64 = sidecar_field(&sidecar, "full_scale_pa")
        .parse()
        .expect("full scale");
    let wav = std::fs::read(root.join("data/listening/reed-sustained.wav")).expect("wav");
    let pascals = wav_to_pascals(&wav, full_scale);
    let lat = fs_psycho::log_attack_time(&pascals, 48_000.0, 480).expect("attack time");
    let receipt = ListeningReceipt {
        listener: "pending".to_string(),
        session: "2026-08-15".to_string(),
        artifact_hex: sidecar_field(&sidecar, "wav_blake3"),
        artifact_ref: "data/listening/reed-sustained.provenance.json".to_string(),
        question: "does the sustained tone read as a reed instrument?".to_string(),
        verdict: ListeningVerdict::Unadjudicated,
        observations: "rendered 2.0 s massless-reed char-line fixture; awaiting the owner's ear"
            .to_string(),
        metrics: fs_psycho::receipt::AttachedMetrics {
            loudness_sone: None,
            sharpness_acum: None,
            log_attack_time: Some(lat),
            spl_db: None,
        },
    };
    let bytes = receipt.to_canonical_bytes().expect("encode");
    println!("{}", String::from_utf8_lossy(&bytes));
}
