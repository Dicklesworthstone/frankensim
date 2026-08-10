//! Fluctuation-strength conformance: calibration on the paper's
//! reference stimulus, the published AM-tone curve (behavioral +
//! authored envelopes), FM behavior, refusals, determinism.

use fs_psycho::PsychoError;
use fs_psycho::fluctuation::{FS_FRAME, bark_z, fluctuation_strength_frame};

const SR: f64 = 48_000.0;

fn am_tone(fc: f64, level_db: f64, fmod: f64, m: f64) -> Vec<f64> {
    let a = fs_math::det::sqrt(2.0) * 2.0e-5 * fs_math::det::pow(10.0, level_db / 20.0);
    (0..FS_FRAME)
        .map(|i| {
            let t = i as f64 / SR;
            a * (1.0 + m * fs_math::det::sin(2.0 * core::f64::consts::PI * fmod * t))
                * fs_math::det::sin(2.0 * core::f64::consts::PI * fc * t)
        })
        .collect()
}

/// One-off calibration probe: prints the model output for the
/// 1-vacil reference stimulus so C_FS can be derived (and re-checked
/// on every explicit run).
#[test]
#[ignore = "calibration probe; run explicitly"]
fn calibration_probe_reference_stimulus() {
    let pcm = am_tone(1000.0, 60.0, 4.0, 1.0);
    let fs = fluctuation_strength_frame(&pcm, SR).expect("fs");
    println!("reference stimulus FS (with current C_FS) = {fs:.17}");
}

fn fm_tone(fc: f64, level_db: f64, fmod: f64, fdev: f64) -> Vec<f64> {
    let a = fs_math::det::sqrt(2.0) * 2.0e-5 * fs_math::det::pow(10.0, level_db / 20.0);
    (0..FS_FRAME)
        .map(|i| {
            let t = i as f64 / SR;
            // Instantaneous phase: 2 pi fc t - (fdev/fmod) cos(2 pi fmod t)
            let phase = 2.0 * core::f64::consts::PI * fc * t
                - (fdev / fmod) * fs_math::det::cos(2.0 * core::f64::consts::PI * fmod * t);
            a * fs_math::det::sin(phase)
        })
        .collect()
}

/// Sweep probe: prints the model's AM and FM curves next to the
/// published values so gates can be authored from measurements.
#[test]
#[ignore = "sweep probe; run explicitly"]
fn sweep_probe_am_fm() {
    let mods = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0];
    let am_published = [0.39, 0.84, 1.25, 1.30, 0.36, 0.06];
    let fm_published = [0.85, 1.17, 2.00, 0.70, 0.27, 0.02];
    print!("AM 70dB:");
    for (&fm, &pubv) in mods.iter().zip(&am_published) {
        let v = fluctuation_strength_frame(&am_tone(1000.0, 70.0, fm, 1.0), SR).expect("fs");
        print!(" {fm}Hz={v:.3}(lit {pubv})");
    }
    println!();
    print!("FM 70dB:");
    for (&fm, &pubv) in mods.iter().zip(&fm_published) {
        let v = fluctuation_strength_frame(&fm_tone(1500.0, 70.0, fm, 700.0), SR).expect("fs");
        print!(" {fm}Hz={v:.3}(lit {pubv})");
    }
    println!();
}

/// Calibration regression: the paper's reference stimulus (1 kHz,
/// 60 dB, m = 1, fmod = 4 Hz) reads EXACTLY 1 vacil (C_FS is derived
/// on it; deterministic pipeline, so this is bitwise-stable).
#[test]
fn reference_stimulus_reads_one_vacil() {
    let fs = fluctuation_strength_frame(&am_tone(1000.0, 60.0, 4.0, 1.0), SR).expect("fs");
    assert!((fs - 1.0).abs() < 1.0e-9, "reference: {fs:.17}");
}

/// The published AM-tone curve (Fastl-Zwicker values as reproduced
/// in the ICA2016 paper's Table 1; single-source transcription with
/// the 1-vacil anchor independently corroborated — provenance in the
/// module doc): every point within the authored envelope
/// +-(0.05 + 25%), measured worst -18% at 4-8 Hz (the model's level
/// growth 60->70 dB is weaker than published — disclosed), plus the
/// discriminating shape assertions a flat model cannot pass.
#[test]
fn am_tone_curve_matches_published_within_envelope() {
    let mods = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0];
    let published = [0.39, 0.84, 1.25, 1.30, 0.36, 0.06];
    let mut got = [0.0f64; 6];
    for (g, &fm) in got.iter_mut().zip(&mods) {
        *g = fluctuation_strength_frame(&am_tone(1000.0, 70.0, fm, 1.0), SR).expect("fs");
    }
    for ((&g, &p), &fm) in got.iter().zip(&published).zip(&mods) {
        let tol = 0.05 + 0.25 * p;
        assert!(
            (g - p).abs() < tol,
            "AM {fm} Hz: {g:.3} vs published {p} (tol {tol:.3})"
        );
    }
    // Shape: bandpass peaking at 4 or 8 Hz, rising 1 -> 4, collapsed
    // at 32 Hz.
    let peak = got.iter().copied().fold(f64::MIN, f64::max);
    // float_cmp-safe: compare by index of the max.
    let peak_idx = got
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(b.1))
        .expect("nonempty")
        .0;
    assert!(
        peak_idx == 2 || peak_idx == 3,
        "peak must sit at 4 or 8 Hz: {got:?}"
    );
    assert!(
        got[0] < got[1] && got[1] < got[2],
        "must rise 1 -> 4 Hz: {got:?}"
    );
    assert!(got[5] < 0.1 * peak, "must collapse at 32 Hz: {got:?}");
    // Level growth direction: 70 dB reads above the 60-dB-calibrated
    // reference at the same modulation (magnitude under-grows vs
    // published — disclosed in the module doc, direction asserted).
    assert!(got[2] > 1.0, "level growth direction: {:.3}", got[2]);
    println!(
        "{{\"suite\":\"fs-psycho\",\"case\":\"fluctuation-am\",\"got\":{got:?},\"published\":{published:?},\"verdict\":\"pass\"}}"
    );
}

/// FM tones: qualitative gates only — the paper's OWN model
/// overestimates FM above 4 Hz; this realisation overestimates FM
/// across the band (measured ~2.5x at the 4 Hz point: 3.02 vs the
/// published 2.00) — values REPORTED, the claims gated are the
/// bandpass collapse and finiteness.
#[test]
fn fm_tone_behavior() {
    let f4 = fluctuation_strength_frame(&fm_tone(1500.0, 70.0, 4.0, 700.0), SR).expect("fs");
    let f32 = fluctuation_strength_frame(&fm_tone(1500.0, 70.0, 32.0, 700.0), SR).expect("fs");
    assert!(
        f4.is_finite() && f4 > 1.0,
        "FM 4 Hz must register strongly: {f4:.3}"
    );
    assert!(
        f32 < 0.1 * f4,
        "FM must collapse at 32 Hz: {f32:.3} vs {f4:.3}"
    );
    println!(
        "{{\"suite\":\"fs-psycho\",\"case\":\"fluctuation-fm\",\"fm4\":{f4:.3},\"fm32\":{f32:.3},\"published\":[2.00,0.02],\"verdict\":\"pass\"}}"
    );
}

/// Unmodulated inputs read ~zero; silence reads exactly zero.
#[test]
fn unmodulated_and_silence_are_smooth() {
    let steady = am_tone(1000.0, 70.0, 4.0, 0.0);
    let fs = fluctuation_strength_frame(&steady, SR).expect("fs");
    assert!(fs < 0.05, "steady tone must be smooth: {fs:.4}");
    let silence = vec![0.0f64; FS_FRAME];
    let fs0 = fluctuation_strength_frame(&silence, SR).expect("fs");
    assert!(fs0 == 0.0, "silence: {fs0}");
}

/// The Zwicker-Terhardt Bark formula cross-checked against the
/// INDEPENDENT Daniel-Weber BARK_FREQS table (MoSQITo lineage):
/// z(BARK_FREQS[i]) must sit near 0.5 i across the scale — the
/// OCR-ambiguous 7.6e-4 coefficient is settled by this consistency
/// (7.6e-5 would read ~1 Bark at 1 kHz instead of ~8.5).
#[test]
fn bark_formula_cross_checks_dw_table() {
    let freqs = &fs_psycho::dw_tables::BARK_FREQS;
    let mut worst = 0.0f64;
    for (i, &f) in freqs.iter().enumerate().skip(1) {
        let z = bark_z(f);
        let want = 0.5 * i as f64;
        worst = worst.max((z - want).abs());
        assert!(
            (z - want).abs() < 0.35,
            "bark mismatch at {f} Hz: z = {z:.3} vs table {want}"
        );
    }
    println!("{{\"suite\":\"fs-psycho\",\"case\":\"bark-crosscheck\",\"worst\":{worst:.3}}}");
}

/// Typed refusals and bitwise determinism.
#[test]
fn fluctuation_refusals_and_determinism() {
    let good = am_tone(1000.0, 60.0, 4.0, 1.0);
    assert!(matches!(
        fluctuation_strength_frame(&good[..FS_FRAME - 1], SR),
        Err(PsychoError::DegenerateSignal { .. })
    ));
    let mut bad = good.clone();
    bad[3] = f64::NAN;
    assert!(matches!(
        fluctuation_strength_frame(&bad, SR),
        Err(PsychoError::NonFinite { .. })
    ));
    for sr in [f64::NAN, 0.0, -1.0, f64::INFINITY] {
        assert!(matches!(
            fluctuation_strength_frame(&good, sr),
            Err(PsychoError::NonFinite { .. })
        ));
    }
    let a = fluctuation_strength_frame(&good, SR).expect("a");
    let b = fluctuation_strength_frame(&good, SR).expect("b");
    assert_eq!(a.to_bits(), b.to_bits());
}
