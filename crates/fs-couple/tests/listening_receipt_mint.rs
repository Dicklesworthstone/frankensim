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

// ---------------------------------------------------------------------
// Brass performance artifact (bead 3ez8g.4.4): attack -> tension walk to
// the upper slot -> crossfaded valve change, rendered by the EMERGENT
// loop (no pitch anywhere in the inputs) — the fixture the trumpet
// listening question adjudicates.
// ---------------------------------------------------------------------

#[test]
#[ignore = "minting run: renders data/listening/brass-emergent.{wav,provenance.json} + receipt bytes"]
fn mint_brass_listening_artifact() {
    use fs_couple::brass_loop::{BrassControl, BrassVoice, LipIsland};
    use fs_couple::mm_line::{MmLineConfig, MmLoad};
    use fs_duct::{Duct, LossModel, Segment, Termination};
    use fs_material::gas::{GasSpec, GasState};

    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air");
    let cfg = MmLineConfig {
        sample_rate_hz: 24_000,
        n_modes: 3,
        extra_slices: 1,
        loss: LossModel::WideTube,
    };
    let bore = |crook: f64| -> Duct {
        let mut segments = vec![Segment::Cylinder {
            radius: 0.006,
            length: 0.30,
        }];
        if crook > 0.0 {
            segments.push(Segment::Cylinder {
                radius: 0.006,
                length: crook,
            });
        }
        segments.push(Segment::Cone {
            inlet_radius: 0.006,
            outlet_radius: 0.012,
            length: 0.25,
        });
        segments.push(Segment::Cone {
            inlet_radius: 0.012,
            outlet_radius: 0.035,
            length: 0.12,
        });
        Duct { segments }
    };
    let combos = [bore(0.0), bore(0.08)];
    let lip = LipIsland {
        mass_kg: 1.8e-4,
        stiffness_n_m: 350.0,
        damping_n_s_m: 0.05,
        width_m: 0.012,
        rest_gap_m: 5.0e-4,
        face_area_m2: 1.0e-4,
        provenance: "card-shaped AUTHORED lip (Estimate)".to_string(),
    };
    let load = MmLoad::Analytic(Termination::FlangedOpen);
    let mut voice = BrassVoice::new(
        &combos,
        &["open", "crook-1"],
        &gas,
        &load,
        &cfg,
        lip,
        1.5e-6,
    )
    .expect("voice");
    voice.apply(BrassControl::SetLipTension(0.8)).expect("t");
    voice
        .apply(BrassControl::SetBlowingPressure(3000.0))
        .expect("p");
    let block_len = 600usize;
    let mut signal = Vec::new();
    let mut block = vec![0.0f64; block_len];
    for b in 0..80usize {
        // 0..32: slot-1 attack+sustain; 32: embouchure to slot 2;
        // 56: crossfaded valve during the note.
        if b == 32 {
            voice.apply(BrassControl::SetLipTension(1.7)).expect("t2");
        }
        if b == 56 {
            voice.apply(BrassControl::SetLipTension(0.8)).expect("t3");
            voice
                .apply(BrassControl::SetValve {
                    combo: 1,
                    fade_samples: 96,
                })
                .expect("valve");
        }
        voice.step_block(&mut block).expect("block");
        signal.extend_from_slice(&block);
    }
    // Remove the DC operating point for the listening artifact (a
    // microphone hears the AC field; the DC mouth pressure is not sound).
    let mean = signal.iter().sum::<f64>() / signal.len() as f64;
    for s in &mut signal {
        *s -= mean;
    }
    let peak = signal.iter().fold(0.0f64, |m, &s| m.max(s.abs()));
    let full_scale_pa = 30_000.0f64;
    assert!(peak < full_scale_pa, "headroom: peak {peak:.0} Pa");
    let (wav, clipped) =
        fs_couple::pcm_wav::encode_pcm16_wav(&signal, 24_000, full_scale_pa).expect("wav");
    assert_eq!(clipped, 0, "never clip a listening artifact");
    let hash = fs_blake3::hash_domain("org.frankensim.music-render.wav.v1", &wav);
    let rms = (signal.iter().map(|s| s * s).sum::<f64>() / signal.len() as f64).sqrt();
    let root = repo_root();
    std::fs::write(root.join("data/listening/brass-emergent.wav"), &wav).expect("wav write");
    let provenance = format!(
        "{{\"schema\":\"frankensim-music-render-provenance-v1\",\"fixture\":\"brass-emergent \
         (attack, slot walk, crossfaded valve; emergent lock, no pitch inputs)\",\
         \"sample_rate_hz\":24000,\"samples\":{},\"block\":{block_len},\
         \"full_scale_pa\":{full_scale_pa:e},\"clipped_samples\":0,\"peak_pa\":{peak:e},\
         \"rms_pa\":{rms:e},\"wav_blake3\":\"{}\",\"encoder\":\"fs_couple::pcm_wav (mono \
         PCM16, never peak-normalized)\"}}\n",
        signal.len(),
        hash.to_hex()
    );
    std::fs::write(
        root.join("data/listening/brass-emergent.provenance.json"),
        provenance,
    )
    .expect("sidecar write");
    let lat = fs_psycho::log_attack_time(&signal, 24_000.0, 240).expect("attack time");
    let receipt = ListeningReceipt {
        listener: "pending".to_string(),
        session: "2026-08-16".to_string(),
        artifact_hex: hash.to_hex(),
        artifact_ref: "data/listening/brass-emergent.provenance.json".to_string(),
        question: "do the attack, slot change, and valve change read as a brass instrument?"
            .to_string(),
        verdict: ListeningVerdict::Unadjudicated,
        observations: "2.0 s emergent brass loop: slot-1 attack, embouchure walk to slot 2, \
                       crossfaded valve during the note; awaiting the owner's ear"
            .to_string(),
        metrics: fs_psycho::receipt::AttachedMetrics {
            loudness_sone: None,
            sharpness_acum: None,
            log_attack_time: Some(lat),
            spl_db: None,
        },
    };
    let bytes = receipt.to_canonical_bytes().expect("encode");
    std::fs::write(
        root.join("data/listening/brass-emergent.listening-receipt"),
        &bytes,
    )
    .expect("receipt write");
    println!(
        "minted brass-emergent: {} samples, peak {peak:.0} Pa, rms {rms:.0} Pa, \
         log-attack-time {lat:.4}, hash {}",
        signal.len(),
        hash.to_hex()
    );
}

#[test]
fn committed_brass_receipt_matches_the_committed_artifact() {
    let root = repo_root();
    let receipt_bytes = std::fs::read(root.join("data/listening/brass-emergent.listening-receipt"))
        .expect("committed brass receipt (mint test)");
    let receipt = ListeningReceipt::from_canonical_bytes(&receipt_bytes).expect("receipt decodes");
    let sidecar =
        std::fs::read_to_string(root.join("data/listening/brass-emergent.provenance.json"))
            .expect("committed sidecar");
    assert_eq!(
        receipt.artifact_hex,
        sidecar_field(&sidecar, "wav_blake3"),
        "receipt must reference exactly the committed artifact"
    );
    let full_scale: f64 = sidecar_field(&sidecar, "full_scale_pa")
        .parse()
        .expect("full scale");
    let wav = std::fs::read(root.join("data/listening/brass-emergent.wav")).expect("wav");
    // The committed WAV bytes must hash to the sidecar digest.
    let hash = fs_blake3::hash_domain("org.frankensim.music-render.wav.v1", &wav);
    assert_eq!(hash.to_hex(), receipt.artifact_hex, "WAV bytes drifted");
    let pascals = wav_to_pascals(&wav, full_scale);
    let lat = fs_psycho::log_attack_time(&pascals, 24_000.0, 240).expect("attack time");
    let attached = receipt
        .metrics
        .log_attack_time
        .expect("brass receipt attaches log-attack-time");
    assert!(
        (lat - attached).abs() < 2.0e-2,
        "attached log-attack-time {attached} must recompute from the WAV ({lat}; PCM16 \
         quantization allows a small drift)"
    );
    assert_eq!(receipt.verdict, ListeningVerdict::Unadjudicated);
    assert!(!receipt.supports_pass());
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"brass-listening-chain\",\"artifact\":\"{}\",\
         \"log_attack_time\":{lat:.6},\"verdict\":\"{}\"}}",
        receipt.artifact_hex,
        receipt.verdict.name()
    );
}
