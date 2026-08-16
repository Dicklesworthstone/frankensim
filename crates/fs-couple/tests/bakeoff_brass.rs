//! Brass MM-vs-plane-wave bake-off (music bead
//! `frankensim-music-v8-root-3ez8g.4.4`): the receipt that DECIDES the
//! image scoping — same flare geometry and gas card for both
//! contenders, QoIs = impedance-peak positions on the trumpet-like
//! flare plus the composed loop's emergent lock. The default test
//! validates the COMMITTED receipt at
//! `tests/receipts/brass-mm-vs-plane.bakeoff` (decode, invariants, and
//! the discrimination: the plane-wave image's flare-peak residual
//! exceeds the multimodal image's). The `--ignored` mint re-measures
//! everything and prints fresh receipt bytes.

use std::collections::BTreeMap;

use fs_blake3::hash_domain;
use fs_couple::bakeoff::{BakeoffOutcome, BakeoffReceipt, ContenderResult};
use fs_couple::brass_loop::{BrassControl, BrassVoice, LipIsland};
use fs_couple::mm_line::{MmLineConfig, MmLoad};
use fs_duct::modal::mm_input_impedance;
use fs_duct::{Duct, LossModel, Segment, Termination};
use fs_material::gas::{GasSpec, GasState};

fn receipt_path() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/receipts/brass-mm-vs-plane.bakeoff")
}

fn air20() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
}

/// The mm-005 trumpet-like flare (the trigger fixture).
fn flare() -> Duct {
    Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.0058,
                length: 0.9,
            },
            Segment::Cone {
                inlet_radius: 0.0058,
                outlet_radius: 0.015,
                length: 0.5,
            },
            Segment::Cone {
                inlet_radius: 0.015,
                outlet_radius: 0.06,
                length: 0.25,
            },
        ],
    }
}

/// The brass-loop bore (the composed-lock fixture).
fn loop_bore() -> Duct {
    Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.006,
                length: 0.30,
            },
            Segment::Cone {
                inlet_radius: 0.006,
                outlet_radius: 0.012,
                length: 0.25,
            },
            Segment::Cone {
                inlet_radius: 0.012,
                outlet_radius: 0.035,
                length: 0.12,
            },
        ],
    }
}

fn lip() -> LipIsland {
    LipIsland {
        mass_kg: 1.8e-4,
        stiffness_n_m: 350.0,
        damping_n_s_m: 0.05,
        width_m: 0.012,
        rest_gap_m: 5.0e-4,
        face_area_m2: 1.0e-4,
        provenance: "card-shaped AUTHORED lip (Estimate)".to_string(),
    }
}

/// Digest of everything both contenders share (geometry + gas + lip).
fn shared_cards_digest() -> fs_blake3::ContentHash {
    let mut bytes = Vec::new();
    for duct in [flare(), loop_bore()] {
        for segment in &duct.segments {
            match *segment {
                Segment::Cylinder { radius, length } => {
                    bytes.extend_from_slice(&radius.to_bits().to_le_bytes());
                    bytes.extend_from_slice(&length.to_bits().to_le_bytes());
                }
                Segment::Cone {
                    inlet_radius,
                    outlet_radius,
                    length,
                } => {
                    bytes.extend_from_slice(&inlet_radius.to_bits().to_le_bytes());
                    bytes.extend_from_slice(&outlet_radius.to_bits().to_le_bytes());
                    bytes.extend_from_slice(&length.to_bits().to_le_bytes());
                }
                Segment::ToneHole { .. } => unreachable!("brass fixtures have no holes"),
            }
        }
    }
    for v in [
        293.15f64, 101_325.0, 1.8e-4, 350.0, 0.05, 0.012, 5.0e-4, 1.0e-4, 1.5e-6,
    ] {
        bytes.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    hash_domain("org.frankensim.fs-couple.brass-bakeoff-cards.v1", &bytes)
}

/// Flare impedance peaks for a mode count, fine staircase (parabolic
/// refined) over the mm-005 band.
fn flare_peaks_hz(n_modes: usize) -> Vec<f64> {
    let gas = air20();
    let duct = flare();
    let mut mags = Vec::new();
    for i in 0..641 {
        let f = 120.0 + (1400.0 - 120.0) * f64::from(i) / 640.0;
        let omega = core::f64::consts::TAU * f;
        let z = mm_input_impedance(
            &duct,
            &gas,
            omega,
            LossModel::WideTube,
            Termination::FlangedOpen,
            n_modes,
            2,
        )
        .expect("authority")
        .plane_impedance;
        mags.push((f, z.abs()));
    }
    let mut peaks = Vec::new();
    for i in 1..mags.len() - 1 {
        if mags[i].1 > mags[i - 1].1 && mags[i].1 > mags[i + 1].1 {
            let (ya, yb, yc) = (mags[i - 1].1.ln(), mags[i].1.ln(), mags[i + 1].1.ln());
            let shift = 0.5 * (ya - yc) / (ya - 2.0 * yb + yc);
            let df = mags[i].0 - mags[i - 1].0;
            peaks.push(mags[i].0 - shift * df);
        }
    }
    peaks
}

/// Composed-loop lock for a mode count on the brass bore.
fn loop_lock_hz(n_modes: usize) -> f64 {
    let gas = air20();
    let cfg = MmLineConfig {
        sample_rate_hz: 24_000,
        n_modes,
        extra_slices: 1,
        loss: LossModel::WideTube,
    };
    let load = MmLoad::Analytic(Termination::FlangedOpen);
    let combos = [loop_bore()];
    let mut voice =
        BrassVoice::new(&combos, &["open"], &gas, &load, &cfg, lip(), 1.5e-6).expect("voice");
    voice
        .apply(BrassControl::SetLipTension(0.8))
        .expect("tension");
    voice
        .apply(BrassControl::SetBlowingPressure(3000.0))
        .expect("pressure");
    let mut block = vec![0.0f64; 2400];
    for _ in 0..8 {
        voice.step_block(&mut block).expect("block");
    }
    // Dominant FFT peak of the settled block.
    use fs_fft::{C64, Fft};
    let n = 4096;
    let mean = block.iter().sum::<f64>() / block.len() as f64;
    let mut buf: Vec<C64> = (0..n)
        .map(|k| C64::new(block.get(k).map_or(0.0, |p| p - mean), 0.0))
        .collect();
    let mut scratch = vec![C64::new(0.0, 0.0); n];
    Fft::new(n).forward(&mut buf, &mut scratch);
    let df = 24_000.0 / n as f64;
    let mags: Vec<f64> = buf[..n / 2]
        .iter()
        .map(|c| (c.re * c.re + c.im * c.im).sqrt())
        .collect();
    let k_lo = (40.0 / df).ceil() as usize;
    let mut best = k_lo;
    for k in k_lo..mags.len() {
        if mags[k] > mags[best] {
            best = k;
        }
    }
    let (ya, yb, yc) = (
        mags[best - 1].max(1e-300).ln(),
        mags[best].ln(),
        mags[best + 1].max(1e-300).ln(),
    );
    let denom = ya - 2.0 * yb + yc;
    let shift = if denom.abs() > 1e-12 {
        0.5 * (ya - yc) / denom
    } else {
        0.0
    };
    (best as f64 - shift) * df
}

#[test]
#[ignore = "minting run: re-measures both contenders and prints fresh receipt bytes"]
fn mint_brass_bakeoff_receipt() {
    let mm_peaks = flare_peaks_hz(4);
    let plane_peaks = flare_peaks_hz(1);
    let n_cmp = mm_peaks.len().min(plane_peaks.len()).min(10);
    let mm_lock = loop_lock_hz(3);
    let plane_lock = loop_lock_hz(1);
    let mut reference = BTreeMap::new();
    let mut mm_measured = BTreeMap::new();
    let mut plane_measured = BTreeMap::new();
    for i in 0..n_cmp {
        let key = format!("flare-peak-{}-hz", i + 1);
        reference.insert(key.clone(), mm_peaks[i]);
        mm_measured.insert(key.clone(), mm_peaks[i]);
        plane_measured.insert(key, plane_peaks[i]);
    }
    reference.insert("loop-lock-hz".to_string(), mm_lock);
    mm_measured.insert("loop-lock-hz".to_string(), mm_lock);
    plane_measured.insert("loop-lock-hz".to_string(), plane_lock);
    let receipt = BakeoffReceipt {
        filling: "brass".to_string(),
        fixture: "crates/fs-duct/src/modal.rs mm-005 flare (fine staircase) + \
                  crates/fs-couple/src/brass_loop.rs open-bore lock at t=0.8 blow=3kPa"
            .to_string(),
        shared_cards: shared_cards_digest(),
        reference,
        contenders: [
            ContenderResult {
                image: "mm-tmm".to_string(),
                owner_crates: vec!["fs-duct".to_string(), "fs-couple".to_string()],
                measured: mm_measured,
                states: 4,
                steps: 641,
                solver_iterations: 0,
                failure_modes: vec![
                    "staircase density is the dominant discretization term (6.7 cents at \
                     default density; ladder disclosed)"
                        .to_string(),
                    "m=0 axisymmetric only; matched-mouth higher modes at the bell".to_string(),
                ],
            },
            ContenderResult {
                image: "plane-wave-char".to_string(),
                owner_crates: vec!["fs-duct".to_string()],
                measured: plane_measured,
                states: 1,
                steps: 641,
                solver_iterations: 0,
                failure_modes: vec![
                    "misses the flare's mode-conversion peak structure (14.9 cents at the \
                     ~400 Hz peak; the brightness mechanism is absent)"
                        .to_string(),
                    "analytic radiation fits carry the ka <= 1 ceiling without a tabulated \
                     load"
                        .to_string(),
                ],
            },
        ],
        outcome: BakeoffOutcome::KeepForSubset {
            narrowed: "plane-wave-char".to_string(),
            subset: "mute/debug claims and tone-hole (wind) bores; never the trumpet flare \
                     claim"
                .to_string(),
        },
        rationale: "the MM image is the flare authority (its own convergence ladder settles \
                    at 0.18 cents while plane misses by 14.9); plane-wave stays on the menu \
                    for the claims it still owns (D21: narrowed, never deleted)"
            .to_string(),
        listening_receipts: Vec::new(),
    };
    receipt.validate().expect("receipt validates");
    let bytes = receipt.to_canonical_bytes().expect("encode");
    std::fs::write(receipt_path(), &bytes).expect("write receipt");
    println!(
        "minted {} ({} bytes), hash {}",
        receipt_path().display(),
        bytes.len(),
        receipt.content_hash().expect("hash").to_hex()
    );
}

#[test]
fn committed_brass_bakeoff_receipt_decides_the_scoping() {
    let bytes = std::fs::read(receipt_path())
        .expect("tests/receipts/brass-mm-vs-plane.bakeoff must be committed (mint test)");
    let receipt = BakeoffReceipt::from_canonical_bytes(&bytes).expect("receipt decodes");
    receipt.validate().expect("receipt validates");
    assert_eq!(receipt.filling, "brass");
    // The discrimination the receipt exists to record: the plane image's
    // worst flare-peak residual exceeds the MM image's by a real margin.
    let mm_res = receipt.residuals(0);
    let plane_res = receipt.residuals(1);
    let worst = |m: &BTreeMap<String, f64>| -> f64 {
        m.iter()
            .filter(|(k, _)| k.starts_with("flare-peak"))
            .map(|(_, v)| *v)
            .fold(0.0f64, f64::max)
    };
    let mm_worst = worst(&mm_res);
    let plane_worst = worst(&plane_res);
    assert!(
        plane_worst > 3.0 * mm_worst.max(1e-12),
        "plane worst {plane_worst:.3e} must exceed MM worst {mm_worst:.3e}"
    );
    // The composed loop's lock also discriminates (recorded, direction
    // logged, both locks near the low slot).
    let lock_ref = receipt.reference["loop-lock-hz"];
    let plane_lock = receipt.contenders[1].measured["loop-lock-hz"];
    let lock_cents = 1200.0 * (plane_lock / lock_ref).log2();
    // The verdict: plane narrowed, MM keeps the flare claim; D21 — nobody
    // is deleted.
    let narrowed_ok = matches!(
        &receipt.outcome,
        BakeoffOutcome::KeepForSubset { narrowed, .. } if narrowed == "plane-wave-char"
    );
    assert!(
        narrowed_ok,
        "outcome must narrow plane-wave-char: {:?}",
        receipt.outcome
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"brass-bakeoff-receipt\",\"verdict\":\"pass\",\
         \"mm_worst_flare_residual\":{mm_worst:.3e},\"plane_worst_flare_residual\":\
         {plane_worst:.3e},\"plane_lock_shift_cents\":{lock_cents:.1},\"hash\":\"{}\"}}",
        receipt.content_hash().expect("hash").to_hex()
    );
}

#[test]
fn brass_gate_summary_enumerates_every_registry_row() {
    // The committed gate-summary artifact (3ez8g.4.4 logging clause)
    // must name EVERY brass row in the registry — an omitted row is a
    // silently-missing gate.
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("workspace root")
        .to_path_buf();
    let registry = std::fs::read_to_string(root.join("instrument-claims.json")).expect("registry");
    let summary = std::fs::read_to_string(root.join("data/claims/brass-gate-summary.tsv"))
        .expect("committed brass gate summary (regenerate when brass rows change)");
    assert!(summary.starts_with("# frankensim-brass-gate-summary-v1"));
    // Every `"filling": "brass"` row's image must appear in the summary.
    let mut brass_images = Vec::new();
    let mut cursor = 0usize;
    while let Some(hit) = registry[cursor..].find("\"filling\": \"brass\"") {
        let at = cursor + hit;
        let tail = &registry[at..];
        let image_tag = "\"image\": \"";
        let img_at = tail.find(image_tag).expect("image field") + image_tag.len();
        let img_end = tail[img_at..].find('"').expect("image end");
        brass_images.push(tail[img_at..img_at + img_end].to_string());
        cursor = at + image_tag.len();
    }
    assert!(
        brass_images.len() >= 6,
        "expected the full brass menu, found {brass_images:?}"
    );
    for image in &brass_images {
        assert!(
            summary.contains(image.as_str()),
            "gate summary omits brass image {image}"
        );
    }
    // The honest headline: the trumpet-claim row is UNGATED in both.
    assert!(summary.contains("trumpet-vertical\tcorpus-cents-gate\tungated"));
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"brass-gate-summary\",\"verdict\":\"pass\",\
         \"rows\":{}}}",
        brass_images.len()
    );
}
