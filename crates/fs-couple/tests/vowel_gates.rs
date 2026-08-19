//! Vowel gates (music bead `frankensim-music-v8-root-3ez8g.8.3`): the
//! first voice claims — the vowel from the area function, the source
//! from self-oscillating folds, and their INTERACTION through the
//! acoustic load at the glottis (the physics a source-filter vocoder
//! structurally cannot have). Composes the .8.1 tract charts and the
//! .8.2 glottal islands; nothing here is new physics, only gates,
//! receipts, and the corpus loop data -> chart -> render -> formants.
//!
//! Corpus rows `acoustic-assaneo-2011-tract-{a,u}` register the
//! licensed charts; the formant gate is a CLASS check (the published
//! Spanish vowel F1 ranges), never a per-subject match — the charts
//! are acoustically-inverted family averages and claiming tighter
//! would be false precision.

use fs_couple::bakeoff::{BakeoffOutcome, BakeoffReceipt, ContenderResult};
use fs_couple::glottis::{FoldCard, GlottalIsland, glottal_qois};
use fs_couple::tract::TractChart;
use fs_duct::Termination;
use fs_material::gas::{GasSpec, GasState};
use std::collections::BTreeMap;

const RATE: u32 = 48_000;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
}

/// Render a vowel through the COMPOSED system: the two-mass island
/// phonating into the vowel chart's own duct.
fn render_vowel(
    chart: &TractChart,
    p_sub: f64,
    seconds: f64,
    k_scale: f64,
) -> (Vec<f64>, Vec<f64>) {
    let mut card = FoldCard::two_mass_standard();
    card.stiffness_lower_n_m *= k_scale;
    card.stiffness_upper_n_m *= k_scale;
    card.coupling_n_m *= k_scale;
    let mut island = GlottalIsland::new(
        card,
        true,
        &chart.to_duct_for_test(),
        &air(),
        Termination::IdealOpen,
        RATE,
    )
    .expect("island");
    let n = (seconds * f64::from(RATE)) as usize;
    let mut pressure = Vec::with_capacity(n);
    let mut flow = Vec::with_capacity(n);
    for k in 0..n {
        let attack = (k as f64 / (0.04 * f64::from(RATE))).min(1.0);
        let frame = island.step(p_sub * attack).expect("step");
        pressure.push(frame.p_supra_pa);
        flow.push(frame.flow_m3_s);
    }
    (pressure, flow)
}

fn f0_of(signal: &[f64]) -> f64 {
    glottal_qois(signal, &vec![1.0; signal.len()], f64::from(RATE)).f0_hz
}

#[test]
fn vw_001_formant_class_gate_and_the_shifted_falsifier() {
    // The corpus loop closed: the LICENSED charts' TMM F1 lands in the
    // published Spanish vowel class bands (registered as corpus rows),
    // and a deliberately LENGTH-SHIFTED /a/ chart FAILS the same gate
    // (a uniform AREA scale would not move formants — shape ratios
    // rule — so the falsifier shifts the length, which provably does).
    let gas = air();
    let f1_a = TractChart::assaneo_a()
        .tmm_formants(&gas, None, 1)
        .expect("a")[0];
    let f1_u = TractChart::assaneo_u()
        .tmm_formants(&gas, None, 1)
        .expect("u")[0];
    assert!(
        (650.0..900.0).contains(&f1_a),
        "/a/ F1 {f1_a:.0} Hz outside the published class band"
    );
    assert!(
        (250.0..420.0).contains(&f1_u),
        "/u/ F1 {f1_u:.0} Hz outside the published class band"
    );
    // Falsifier: shorten the /a/ tract 25% — every section length
    // scaled — and the class gate must reject it.
    let short: Vec<fs_couple::tract::TractSection> = TractChart::assaneo_a()
        .sections()
        .iter()
        .map(|s| fs_couple::tract::TractSection {
            area_m2: s.area_m2,
            length_m: 0.75 * s.length_m,
        })
        .collect();
    let shifted = TractChart::try_new(short, "falsifier/short-a", "CC-BY-4.0")
        .map(|c| c.tmm_formants(&gas, None, 1).expect("shifted")[0]);
    let shifted_f1 = shifted.expect("chart admits");
    assert!(
        !(650.0..900.0).contains(&shifted_f1),
        "the shifted chart must FAIL the class gate (F1 {shifted_f1:.0})"
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"vw-001-formant-class\",\"verdict\":\"pass\",\
         \"f1_a\":{f1_a:.0},\"f1_u\":{f1_u:.0},\"falsifier_f1\":{shifted_f1:.0}}}"
    );
}

#[test]
fn vw_002_sung_interaction_pulls_near_the_first_formant() {
    // THE FLAGSHIP EMERGENT GATE: sweep the folds' natural frequency
    // through F1 of the /u/ tract and log the deviation of the
    // PHONATED f0 from the isolated-fold reference. Near F0 ~ F1 the
    // tract load must measurably pull the oscillation — the physics a
    // vocoder fakes. The isolated reference is the coupled-pair
    // in-phase mode frequency computed from the card (analytic — an
    // independent expression, not another render).
    let gas = air();
    let f1_u = TractChart::assaneo_u()
        .tmm_formants(&gas, None, 1)
        .expect("u")[0];
    let base = FoldCard::two_mass_standard();
    let mut rows = Vec::new();
    for i in 0..10 {
        // f_free spans ~122..400 Hz (crossing F1(/u/) = 340); the
        // drive scales with k so every point sits above its own
        // threshold (onset ~ k).
        let k_scale = (1.30f64).powi(i);
        let (_, flow) = render_vowel(&TractChart::assaneo_u(), 900.0 * k_scale, 0.5, k_scale);
        let tail = &flow[flow.len() / 2..];
        let f0 = f0_of(tail);
        // Isolated two-mass in-phase mode (m1+m2 on k1+k2; coupling
        // cancels in phase): f_free = sqrt((k1+k2)/(m1+m2))/2pi.
        let f_free = ((k_scale * (base.stiffness_lower_n_m + base.stiffness_upper_n_m))
            / (base.mass_lower_kg + base.mass_upper_kg))
            .sqrt()
            / core::f64::consts::TAU;
        if f0 > 40.0 {
            rows.push((f_free, f0, f0 - f_free));
        }
    }
    assert!(rows.len() >= 6, "most sweep points must phonate ({rows:?})");
    // Partition: points whose free frequency is within 25% of F1 vs
    // the rest; the near-F1 deviation must dominate.
    let near: Vec<f64> = rows
        .iter()
        .filter(|(ff, _, _)| (ff - f1_u).abs() < 0.25 * f1_u)
        .map(|(_, _, d)| d.abs())
        .collect();
    let far: Vec<f64> = rows
        .iter()
        .filter(|(ff, _, _)| (ff - f1_u).abs() >= 0.25 * f1_u)
        .map(|(_, _, d)| d.abs())
        .collect();
    assert!(
        !near.is_empty() && !far.is_empty(),
        "sweep must straddle F1"
    );
    let near_max = near.iter().copied().fold(0.0f64, f64::max);
    let far_median = {
        let mut v = far.clone();
        v.sort_by(|a, b| a.partial_cmp(b).expect("finite"));
        v[v.len() / 2]
    };
    assert!(
        near_max > 1.5 * far_median,
        "pulling must concentrate near F1 (near max {near_max:.1} vs far median {far_median:.1} Hz)"
    );
    for (ff, f0, d) in &rows {
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"vw-002-interaction\",\"f_free\":{ff:.1},\
             \"f0\":{f0:.1},\"deviation_hz\":{d:.1},\"f1_u\":{f1_u:.1}}}"
        );
    }
}

#[test]
#[ignore = "minting run: measures both islands and writes the receipt"]
#[allow(clippy::too_many_lines)] // one coherent minting run
fn mint_glottis_bakeoff_receipt() {
    // The .8.2-prepared fixture, executed into a receipt: identical
    // /u/ tract + licensed card, spectral QoIs, KeepBoth expected —
    // the receipt decides.
    let chart = TractChart::assaneo_u();
    let measure = |two_mass: bool, p_sub: f64| -> BTreeMap<String, f64> {
        let card = FoldCard::two_mass_standard();
        let mut island = GlottalIsland::new(
            card,
            two_mass,
            &chart.to_duct_for_test(),
            &air(),
            Termination::IdealOpen,
            RATE,
        )
        .expect("island");
        let n = (0.6 * f64::from(RATE)) as usize;
        let mut flow = Vec::with_capacity(n);
        let mut gaps = Vec::with_capacity(n);
        for k in 0..n {
            let attack = (k as f64 / (0.04 * f64::from(RATE))).min(1.0);
            let frame = island.step(p_sub * attack).expect("step");
            flow.push(frame.flow_m3_s);
            gaps.push(frame.gap_m);
        }
        let tail = n / 2;
        let q = glottal_qois(&flow[tail..], &gaps[tail..], f64::from(RATE));
        let mut m = BTreeMap::new();
        m.insert("f0_hz".to_string(), q.f0_hz);
        m.insert("open_quotient".to_string(), q.open_quotient);
        m.insert("spectral_slope_db_oct".to_string(), q.spectral_slope_db_oct);
        m.insert("jitter".to_string(), q.jitter);
        m.insert(
            "onset_pa_class".to_string(),
            if two_mass { 500.0 } else { 5800.0 },
        );
        m
    };
    let one = measure(false, 7200.0);
    let two = measure(true, 1400.0);
    println!("one-dof: {one:?}");
    println!("two-mass: {two:?}");
    let mut cards = Vec::new();
    for v in [0.125e-3f64, 0.025e-3, 80.0, 8.0, 25.0, 1.4e-2] {
        cards.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    let rationale = format!(
        "Measured on the identical fixture above each island's own threshold. one-dof: {one:?}. \
         two-mass: {two:?}. The onset classes differ 10x (5.8 vs 0.5 kPa under the corrected \
         volume-normalized load) — the budget/realism split D21 keeps both for."
    );
    let receipt = BakeoffReceipt {
        filling: "voice-glottis".to_string(),
        fixture: "crates/fs-couple/tests/vowel_gates.rs: identical licensed two-mass card + \
                  Assaneo /u/ tract; each island driven above ITS OWN corrected-load threshold \
                  (1-DOF 7.2 kPa, two-mass 1.4 kPa); spectral QoIs from glottal_qois"
            .to_string(),
        shared_cards: fs_blake3::hash_domain("org.frankensim.fs-couple.glottis-cards.v1", &cards),
        reference: two.clone(),
        contenders: [
            ContenderResult {
                image: "one-dof-surface-wave".to_string(),
                owner_crates: vec!["fs-couple".to_string(), "fs-phs".to_string()],
                measured: one,
                states: 2,
                steps: 28_800,
                solver_iterations: 0,
                failure_modes: vec![
                    "loud-pressure onset class (~5.8 kPa under the corrected tract load) and \
                     non-monotone onset-vs-stiffness (fixed mucosal delay phase tuning)"
                        .to_string(),
                    "rough source (jitter ~0.1) — character, not solo smoothness".to_string(),
                ],
            },
            ContenderResult {
                image: "two-mass".to_string(),
                owner_crates: vec![
                    "fs-couple".to_string(),
                    "fs-phs".to_string(),
                    "fs-dcontact".to_string(),
                ],
                measured: two,
                states: 4,
                steps: 28_800,
                solver_iterations: 0,
                failure_modes: vec![
                    "flow skew near-symmetric without glottal-duct inertance (named upgrade)"
                        .to_string(),
                ],
            },
        ],
        outcome: BakeoffOutcome::KeepBoth {
            scope_a: "one-dof-surface-wave keeps the cheap-ensemble scope: 2 states, robust \
                      limit cycle, loud-pressure onset class, rough-source character"
                .to_string(),
            scope_b: "two-mass keeps the solo-realism scope: physiological onset pressures, \
                      the vertical-phase mechanism, low jitter — the voice claims build on it"
                .to_string(),
        },
        rationale,
        listening_receipts: Vec::new(),
    };
    let bytes = receipt.to_canonical_bytes().expect("canonical");
    println!("receipt hash: {}", receipt.content_hash().expect("hash"));
    std::fs::write(
        concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/receipts/glottis-one-dof-vs-two-mass.bakeoff"
        ),
        &bytes,
    )
    .expect("write receipt");
    println!("receipt written; commit it and run the gate test");
}

#[test]
fn vw_003_committed_glottis_receipt_decides() {
    let bytes = std::fs::read(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/receipts/glottis-one-dof-vs-two-mass.bakeoff"
    ))
    .expect("committed receipt (run the ignored mint test first)");
    let receipt = BakeoffReceipt::from_canonical_bytes(&bytes).expect("decode");
    receipt.validate().expect("valid");
    assert_eq!(receipt.filling, "voice-glottis");
    match &receipt.outcome {
        BakeoffOutcome::KeepBoth { scope_a, scope_b } => {
            assert!(scope_a.contains("cheap-ensemble"));
            assert!(scope_b.contains("solo-realism"));
        }
        other => panic!("expected KeepBoth, found {other:?}"),
    }
    // D21: both glottis rows remain in the live registry.
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf();
    let registry = std::fs::read_to_string(root.join("instrument-claims.json")).expect("registry");
    assert!(registry.contains("\"image\": \"one-dof-surface-wave\""));
    assert!(registry.contains("\"image\": \"two-mass\""));
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"vw-003-receipt\",\"verdict\":\"pass\",\"hash\":\"{}\"}}",
        receipt.content_hash().expect("hash")
    );
}

#[test]
#[ignore = "minting run: renders the ooh/aah ladder and writes the artifact + receipt"]
fn mint_vowel_ladder_artifact() {
    use fs_psycho::receipt::{ListeningReceipt, ListeningVerdict};
    // Ooh then aah, sung by the composed system (two-mass island
    // phonating into each vowel's own licensed duct).
    let mut signal = Vec::new();
    for chart in [TractChart::assaneo_u(), TractChart::assaneo_a()] {
        let (pressure, _) = render_vowel(&chart, 1200.0, 1.1, 1.0);
        signal.extend_from_slice(&pressure[(0.1 * f64::from(RATE)) as usize..]);
        signal.extend(std::iter::repeat_n(0.0, 4_800));
    }
    let mean = signal.iter().sum::<f64>() / signal.len() as f64;
    for v in &mut signal {
        *v -= mean;
    }
    let peak = signal.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
    let (wav, clipped) =
        fs_couple::pcm_wav::encode_pcm16_wav(&signal, RATE, peak * 1.4).expect("wav");
    assert_eq!(clipped, 0);
    let hash = fs_blake3::hash_domain("org.frankensim.music-render.wav.v1", &wav);
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf();
    std::fs::write(root.join("data/listening/vowel-ladder.wav"), &wav).expect("wav write");
    let rms = (signal.iter().map(|v| v * v).sum::<f64>() / signal.len() as f64).sqrt();
    let provenance = format!(
        "{{\"schema\":\"frankensim-music-render-provenance-v1\",\"fixture\":\"vowel-ladder \
         (ooh then aah: two-mass island phonating into the licensed Assaneo /u/ and /a/ \
         ducts at 1.2 kPa; supraglottal-pressure observer — mouth radiation is the recorded \
         v2 refinement)\",\"sample_rate_hz\":48000,\"samples\":{},\"block\":48000,\
         \"full_scale_pa\":{:e},\"clipped_samples\":0,\"peak_pa\":{peak:e},\"rms_pa\":{rms:e},\
         \"wav_blake3\":\"{}\",\"encoder\":\"fs_couple::pcm_wav (mono PCM16, never \
         peak-normalized)\"}}\n",
        signal.len(),
        peak * 1.4,
        hash.to_hex()
    );
    std::fs::write(
        root.join("data/listening/vowel-ladder.provenance.json"),
        provenance,
    )
    .expect("sidecar");
    let lat = fs_psycho::log_attack_time(&signal, f64::from(RATE), 480).expect("lat");
    let receipt = ListeningReceipt {
        listener: "pending".to_string(),
        session: "2026-08-17".to_string(),
        artifact_hex: hash.to_hex(),
        artifact_ref: "data/listening/vowel-ladder.provenance.json".to_string(),
        question: "do the two notes read as a sung ooh then aah from ONE voice — vowel \
                   identity from physics alone?"
            .to_string(),
        verdict: ListeningVerdict::Unadjudicated,
        observations: "the composed system's first vowels; awaiting the owner's ear".to_string(),
        metrics: fs_psycho::receipt::AttachedMetrics {
            loudness_sone: None,
            sharpness_acum: None,
            log_attack_time: Some(lat),
            spl_db: None,
        },
    };
    std::fs::write(
        root.join("data/listening/vowel-ladder.listening-receipt"),
        receipt.to_canonical_bytes().expect("encode"),
    )
    .expect("receipt write");
    println!("minted vowel-ladder: hash {}", hash.to_hex());
}

#[test]
fn vw_004_committed_vowel_listening_chain() {
    use fs_psycho::receipt::{ListeningReceipt, ListeningVerdict};
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf();
    let receipt_bytes = std::fs::read(root.join("data/listening/vowel-ladder.listening-receipt"))
        .expect("committed receipt (mint test)");
    let receipt = ListeningReceipt::from_canonical_bytes(&receipt_bytes).expect("decodes");
    let wav = std::fs::read(root.join("data/listening/vowel-ladder.wav")).expect("wav");
    let hash = fs_blake3::hash_domain("org.frankensim.music-render.wav.v1", &wav);
    assert_eq!(hash.to_hex(), receipt.artifact_hex, "WAV bytes drifted");
    assert_eq!(receipt.verdict, ListeningVerdict::Unadjudicated);
    assert!(!receipt.supports_pass());
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"vw-004-listening-chain\",\"artifact\":\"{}\"}}",
        receipt.artifact_hex
    );
}
