//! fs-psycho conformance battery: ISO 532-1 Annex B reference values
//! (test signal 1 exact-path; tone/noise signals cross-path), DIN
//! sharpness behavior, calibration refusals, monotonicity and scale
//! properties, determinism, the dropped-table mutation, and the
//! listening-law pin.

use fs_psycho::{
    Calibration, LISTENING_LAW, N_THIRD_OCTAVE_BANDS, PsychoError, SoundField, log_attack_time,
    loudness_stationary, sharpness_din, spl_from_pcm_rms,
};

/// ISO 532-1 Annex B.2 test signal 1: the 28 third-octave levels
/// verbatim from the standard's distributed file
/// "Annex B.2/Test signal 1.txt".
const TEST_SIGNAL_1: [f64; 28] = [
    -60.0, -60.0, 78.0, 79.0, 89.0, 72.0, 80.0, 89.0, 75.0, 87.0, 85.0, 79.0, 86.0, 80.0, 71.0,
    70.0, 72.0, 71.0, 72.0, 74.0, 69.0, 65.0, 67.0, 77.0, 68.0, 58.0, 45.0, 30.0,
];

/// Reference: ISO 532-1 Annex B.2 expected loudness for test signal 1
/// (free field), as distributed in the standard's results workbook and
/// reproduced by the MoSQITo validation suite: 83.296 sone, with the
/// standard's stationary compliance criterion (5% relative,
/// 0.1 sone absolute).
#[test]
fn iso_test_signal_1_reference_sones() {
    let out = loudness_stationary(&TEST_SIGNAL_1, SoundField::Free).expect("loudness");
    // The port is statement-for-statement on the SAME path as the
    // reference (compiled reference binary at full precision:
    // 83.29566042436214), so the gate is EXACTNESS at 1e-9 relative —
    // hiding a port bug inside the standard's 5% compliance tolerance
    // is exactly what caught the -.25f transcription mangling, and a
    // 1e-3 gate is what hid the float32-literal table rounding
    // (C's `-0.6f` widened into a double is NOT 0.6) until the
    // signal-path exactness work surfaced it.
    let reference = 83.295_660_424_362_14;
    assert!(
        ((out.sones - reference) / reference).abs() <= 1.0e-9,
        "test signal 1: {:.17} sone vs reference {reference:.17}",
        out.sones
    );
    println!(
        "{{\"suite\":\"fs-psycho\",\"case\":\"iso-signal-1\",\"sones\":{:.3},\"reference\":{reference},\"verdict\":\"pass\"}}",
        out.sones
    );
}

/// Single-band level vector: `level` dB in the band whose center is
/// `freq_hz`, quiet (-60 dB) elsewhere.
fn tone_levels(freq_hz: f64, level: f64) -> [f64; 28] {
    let centers = [
        25.0, 31.5, 40.0, 50.0, 63.0, 80.0, 100.0, 125.0, 160.0, 200.0, 250.0, 315.0, 400.0, 500.0,
        630.0, 800.0, 1000.0, 1250.0, 1600.0, 2000.0, 2500.0, 3150.0, 4000.0, 5000.0, 6300.0,
        8000.0, 10000.0, 12500.0,
    ];
    let mut v = [-60.0; 28];
    let idx = centers
        .iter()
        .position(|&c| (c - freq_hz).abs() / freq_hz < 0.05)
        .expect("band center");
    v[idx] = level;
    v
}

/// Cross-path references (ISO Annex B.3 wav signals, reference sones
/// from the standard's workbook via the MoSQITo validation suite):
/// a single-band level vector under-represents the wav path, whose
/// third-octave filterbank leaks ~-20 dB into neighbor bands that
/// loudness compression turns into REAL extra sones (measured: the
/// 1 kHz tone reads 13% under the wav reference). The authored
/// cross-path gate is therefore 15% + 0.1 sone; the EXACT-path pin at
/// the standard tolerance is test signal 1.
#[test]
fn iso_tone_and_reference_sones_cross_path() {
    for (freq, level, reference) in [
        (250.0, 80.0, 14.655),
        (1000.0, 60.0, 4.019),
        (4000.0, 40.0, 1.549),
    ] {
        let out =
            loudness_stationary(&tone_levels(freq, level), SoundField::Free).expect("loudness");
        let tol = 0.15 * reference + 0.1;
        assert!(
            (out.sones - reference).abs() <= tol,
            "{freq} Hz {level} dB: {} sone vs {reference}",
            out.sones
        );
    }
    println!("{{\"suite\":\"fs-psycho\",\"case\":\"iso-tones\",\"verdict\":\"pass\"}}");
}

#[test]
fn loudness_monotone_in_level_and_1khz_40db_anchor() {
    // The sone scale's anchor: 1 kHz at 40 dB = 1 sone (definition;
    // the reference implementation reproduces it on the levels path).
    let out = loudness_stationary(&tone_levels(1000.0, 40.0), SoundField::Free).expect("l");
    assert!(
        (out.sones - 1.0).abs() <= 0.15,
        "1 kHz 40 dB anchor: {} sone",
        out.sones
    );
    // Monotone in level; doubling per 10 dB in the plateau region
    // (Zwicker's law) within 15%.
    let mut prev = 0.0;
    for level in [40.0, 50.0, 60.0, 70.0, 80.0] {
        let n = loudness_stationary(&tone_levels(1000.0, level), SoundField::Free)
            .expect("l")
            .sones;
        assert!(n > prev, "loudness must be monotone in level");
        if prev > 0.0 {
            let ratio = n / prev;
            assert!(
                (1.7..2.4).contains(&ratio),
                "10 dB should roughly double sones: ratio {ratio:.3} at {level} dB"
            );
        }
        prev = n;
    }
}

#[test]
fn sharpness_increases_with_frequency_and_refuses_silence() {
    // Sharpness is driven by high-Bark content: a 4 kHz tone must be
    // sharper than a 250 Hz tone at the same level; the DIN weighting
    // kicks in above 15.8 Bark.
    let s_low = sharpness_din(
        &loudness_stationary(&tone_levels(250.0, 70.0), SoundField::Free)
            .expect("l")
            .specific,
    )
    .expect("s");
    let s_high = sharpness_din(
        &loudness_stationary(&tone_levels(4000.0, 70.0), SoundField::Free)
            .expect("l")
            .specific,
    )
    .expect("s");
    assert!(
        s_high > 2.0 * s_low,
        "sharpness must rise with frequency: {s_high:.3} vs {s_low:.3}"
    );
    // Approximate acum scale sanity: DIN's reference point is a 1 kHz
    // narrowband at 60 dB = 1 acum; a 1 kHz TONE on the levels path
    // lands near it (authored envelope, cross-signal-class).
    let s_ref = sharpness_din(
        &loudness_stationary(&tone_levels(1000.0, 60.0), SoundField::Free)
            .expect("l")
            .specific,
    )
    .expect("s");
    assert!(
        (0.7..1.3).contains(&s_ref),
        "1 kHz 60 dB sharpness {s_ref:.3} far from the 1-acum anchor"
    );
    // Silence refuses BY NAME, never fabricates.
    let silent = vec![0.0; 240];
    assert!(matches!(
        sharpness_din(&silent),
        Err(PsychoError::DegenerateSignal { .. })
    ));
}

#[test]
fn dropped_critical_band_mutation_exceeds_tolerance() {
    // MUTATION: zeroing the low-frequency critical-band aggregation
    // (the DLL/LCB machinery) by feeding the low bands as silence
    // shifts test signal 1 far outside the standard tolerance — the
    // low-band tables are load-bearing, not decorative.
    let mut mutated = TEST_SIGNAL_1;
    for v in mutated.iter_mut().take(11) {
        *v = -60.0;
    }
    let full = loudness_stationary(&TEST_SIGNAL_1, SoundField::Free).expect("l");
    let cut = loudness_stationary(&mutated, SoundField::Free).expect("l");
    let reference = 83.296;
    let tol = 0.05 * reference + 0.1;
    assert!(
        (cut.sones - reference).abs() > 3.0 * tol,
        "dropping the low bands must be visible: {} vs {}",
        cut.sones,
        full.sones
    );
}

#[test]
fn uncalibrated_absolute_refuses_and_calibrated_works() {
    let pcm: Vec<f64> = (0i32..4800)
        .map(|i| {
            0.5 * fs_math::det::sin(2.0 * core::f64::consts::PI * 1000.0 * f64::from(i) / 48000.0)
        })
        .collect();
    assert!(matches!(
        spl_from_pcm_rms(&pcm, None),
        Err(PsychoError::UncalibratedAbsolute)
    ));
    let spl = spl_from_pcm_rms(
        &pcm,
        Some(Calibration {
            db_spl_at_full_scale: 94.0,
        }),
    )
    .expect("calibrated");
    // Half-scale sine = -6.02 dB re full scale.
    assert!(
        (spl - (94.0 - 6.0206)).abs() < 0.01,
        "calibrated SPL {spl:.3}"
    );
}

#[test]
fn log_attack_time_orders_fast_vs_slow_attacks() {
    let sr = 48_000.0;
    let make = |attack_s: f64| -> Vec<f64> {
        (0i32..48_000)
            .map(|i| {
                let t = f64::from(i) / sr;
                let env = (t / attack_s).min(1.0);
                env * fs_math::det::sin(2.0 * core::f64::consts::PI * 440.0 * t)
            })
            .collect()
    };
    let fast = log_attack_time(&make(0.005), sr, 128).expect("fast");
    let slow = log_attack_time(&make(0.2), sr, 128).expect("slow");
    assert!(
        slow > fast + 1.0,
        "attack ordering: slow {slow:.3} vs fast {fast:.3} (log10 s)"
    );
    // The measured 10-90% rise of a linear ramp is 0.8 * attack time.
    let expected_slow = fs_math::det::ln(0.8 * 0.2) / fs_math::det::ln(10.0);
    assert!(
        (slow - expected_slow).abs() < 0.15,
        "slow LAT {slow:.3} vs expected {expected_slow:.3}"
    );
    // Silence refuses.
    assert!(matches!(
        log_attack_time(&vec![0.0; 1000], sr, 64),
        Err(PsychoError::DegenerateSignal { .. })
    ));
}

#[test]
fn determinism_bitwise_and_typed_refusals() {
    let a = loudness_stationary(&TEST_SIGNAL_1, SoundField::Free).expect("a");
    let b = loudness_stationary(&TEST_SIGNAL_1, SoundField::Free).expect("b");
    assert_eq!(a.sones.to_bits(), b.sones.to_bits());
    for (x, y) in a.specific.iter().zip(&b.specific) {
        assert_eq!(x.to_bits(), y.to_bits());
    }
    // Typed refusals: wrong length, NaN.
    assert!(matches!(
        loudness_stationary(&[0.0; 27], SoundField::Free),
        Err(PsychoError::Shape { .. })
    ));
    let mut bad = TEST_SIGNAL_1;
    bad[5] = f64::NAN;
    assert!(matches!(
        loudness_stationary(&bad, SoundField::Free),
        Err(PsychoError::NonFinite { .. })
    ));
    assert!(matches!(
        sharpness_din(&[0.0; 10]),
        Err(PsychoError::Shape { .. })
    ));
    let _ = N_THIRD_OCTAVE_BANDS;
}

#[test]
fn diffuse_field_differs_from_free_field() {
    // The DDF table must actually act: diffuse-field loudness of a
    // high-frequency tone differs measurably from free-field.
    let free = loudness_stationary(&tone_levels(4000.0, 70.0), SoundField::Free).expect("f");
    let diff = loudness_stationary(&tone_levels(4000.0, 70.0), SoundField::Diffuse).expect("d");
    assert!(
        (free.sones - diff.sones).abs() / free.sones > 0.02,
        "DDF must be live: {} vs {}",
        free.sones,
        diff.sones
    );
}

#[test]
fn listening_law_is_pinned() {
    // The not-a-substitute-for-listening statement is data, asserted
    // here so removing it breaks the build's tests, and present in
    // the crate docs.
    assert!(LISTENING_LAW.contains("never a substitute for human listening"));
}

/// 100% AM tone at 1 kHz carrier, 60 dB SPL: the Daniel-Weber
/// band-pass character must PEAK near 70 Hz modulation (the published
/// signature) at ~1 asper (the model's calibration anchor), falling
/// off on both sides; an unmodulated tone is smooth (near-zero
/// roughness).
#[test]
fn roughness_am_tone_peaks_near_70_hz() {
    use fs_psycho::roughness::{DW_BLOCK, roughness_dw_block};
    let sr = 48_000.0;
    let level_db = 60.0;
    let p_amp = fs_math::det::sqrt(2.0) * 2.0e-5 * fs_math::det::pow(10.0, level_db / 20.0);
    let make = |fmod: f64| -> Vec<f64> {
        (0..DW_BLOCK)
            .map(|i| {
                let t = i as f64 / sr;
                let carrier = fs_math::det::sin(2.0 * core::f64::consts::PI * 1000.0 * t);
                let am = 1.0 + fs_math::det::sin(2.0 * core::f64::consts::PI * fmod * t);
                // 100% AM, normalized so the CARRIER keeps 60 dB rms.
                p_amp * am * carrier
            })
            .collect()
    };
    let mods = [20.0, 40.0, 55.0, 70.0, 90.0, 120.0, 160.0];
    let r: Vec<f64> = mods
        .iter()
        .map(|&fm| roughness_dw_block(&make(fm), sr).expect("roughness"))
        .collect();
    let peak_idx = r
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(b.1))
        .expect("nonempty")
        .0;
    let peak_mod = mods[peak_idx];
    assert!(
        (55.0..=90.0).contains(&peak_mod),
        "roughness must peak near 70 Hz: peaked at {peak_mod} Hz, curve {r:?}"
    );
    // EXACTNESS pin (the loudness lesson: never hide a port bug in a
    // behavioral envelope): the Apache-2.0 reference implementation
    // (hw._H_weighting + main_calc), re-run standalone at this exact
    // block length on this exact signal, gives these seven values —
    // pinned at 1e-12 RELATIVE (an adversarial review caught the
    // first port reading only 7-9 matching digits: its H-weighting
    // support ignored the reference's floor() bin truncation; a loose
    // 1e-3 pin had hidden it). R(70 Hz) also sits at the published
    // ~1-asper anchor.
    let reference = [
        0.197_327_796_978_835_65,
        0.660_982_441_631_880_9,
        0.949_126_590_618_143_1,
        1.044_791_270_606_281_9,
        0.940_671_646_615_856_8,
        0.641_113_651_282_113_5,
        0.335_153_389_958_774_8,
    ];
    for ((&got, &want), &fm) in r.iter().zip(&reference).zip(&mods) {
        assert!(
            ((got - want) / want).abs() < 1.0e-12,
            "R({fm} Hz) = {got:.17} vs reference {want:.17}"
        );
    }
    // Falloff on both sides of the peak region.
    assert!(r[0] < 0.7 * r[peak_idx], "low-side falloff: {r:?}");
    assert!(
        *r.last().expect("nonempty") < 0.7 * r[peak_idx],
        "high-side falloff: {r:?}"
    );
    // Unmodulated tone: nearly smooth.
    let steady: Vec<f64> = (0..DW_BLOCK)
        .map(|i| p_amp * fs_math::det::sin(2.0 * core::f64::consts::PI * 1000.0 * i as f64 / sr))
        .collect();
    let r0 = roughness_dw_block(&steady, sr).expect("steady");
    assert!(
        r0 < 0.2 * r[peak_idx],
        "unmodulated tone must be far smoother: {r0} vs {}",
        r[peak_idx]
    );
    println!(
        "{{\"suite\":\"fs-psycho\",\"case\":\"roughness-am\",\"mods_hz\":{mods:?},\"asper\":{r:?},\"verdict\":\"pass\"}}"
    );
}

/// Roughness refusal paths are TYPED, by name (review-caught: the
/// happy path alone left DegenerateSignal/NonFinite unexecuted, and
/// an infinite sample rate slipped a NaN-only guard to fabricate a
/// value).
#[test]
fn roughness_refusals_are_typed() {
    use fs_psycho::PsychoError;
    use fs_psycho::roughness::{DW_BLOCK, roughness_dw_block};
    let good = vec![0.01f64; DW_BLOCK];
    // Short input.
    assert!(matches!(
        roughness_dw_block(&good[..DW_BLOCK - 1], 48_000.0),
        Err(PsychoError::DegenerateSignal { .. })
    ));
    // Non-finite sample.
    let mut bad = good.clone();
    bad[7] = f64::NAN;
    assert!(matches!(
        roughness_dw_block(&bad, 48_000.0),
        Err(PsychoError::NonFinite { .. })
    ));
    // Bad sample rates: NaN, zero, negative, and INFINITE (which
    // once produced Ok(garbage) through an all-inf frequency axis).
    for sr in [f64::NAN, 0.0, -48_000.0, f64::INFINITY] {
        assert!(
            matches!(
                roughness_dw_block(&good, sr),
                Err(PsychoError::NonFinite { .. })
            ),
            "sample rate {sr} must refuse"
        );
    }
}

/// Deterministic PCM fixtures shared by the signal-path tests and,
/// via `dump_reference_signals`, by the ISO reference binary run that
/// establishes the pinned ground-truth values (bit-identical input on
/// both sides).
mod sigfix {
    pub const SR: f64 = 48_000.0;
    pub const LEN: usize = 72_000; // 1.5 s

    fn amp(level_db: f64) -> f64 {
        fs_math::det::sqrt(2.0) * 2.0e-5 * fs_math::det::pow(10.0, level_db / 20.0)
    }

    /// Steady tone.
    pub fn tone(freq: f64, level_db: f64) -> Vec<f64> {
        let a = amp(level_db);
        (0..LEN)
            .map(|i| a * fs_math::det::sin(2.0 * core::f64::consts::PI * freq * i as f64 / SR))
            .collect()
    }

    /// Deterministic uniform noise in (-0.5, 0.5) from a 32-bit LCG
    /// (Numerical Recipes constants) — reproducible bit-exactly in
    /// any language, so the reference run consumes identical bytes.
    pub fn lcg_noise(n: usize, seed: u32) -> Vec<f64> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                f64::from(x) / 4_294_967_296.0 - 0.5
            })
            .collect()
    }

    /// Tonality fixture length (power of two for fs-fft).
    pub const TONAL_LEN: usize = 65_536;

    /// Tonal fixture A: 1 kHz 60 dB tone + 3 kHz at half amplitude +
    /// broadband LCG noise (~40 dB band level).
    pub fn tonal_a() -> Vec<f64> {
        let a = amp(60.0);
        let noise = lcg_noise(TONAL_LEN, 12_345);
        (0..TONAL_LEN)
            .map(|i| {
                let t = i as f64 / SR;
                a * fs_math::det::sin(2.0 * core::f64::consts::PI * 1000.0 * t)
                    + 0.5 * a * fs_math::det::sin(2.0 * core::f64::consts::PI * 3000.0 * t)
                    + 0.0007 * noise[i]
            })
            .collect()
    }

    /// Noise-only fixture B (no deterministic tonal component).
    pub fn noise_b() -> Vec<f64> {
        lcg_noise(TONAL_LEN, 777)
            .into_iter()
            .map(|v| 0.028 * v)
            .collect()
    }

    /// Low-SNR fixture C: 500 Hz tone barely above the noise.
    pub fn lowsnr_c() -> Vec<f64> {
        let noise = lcg_noise(TONAL_LEN, 999);
        (0..TONAL_LEN)
            .map(|i| {
                let t = i as f64 / SR;
                0.0089 * fs_math::det::sin(2.0 * core::f64::consts::PI * 500.0 * t)
                    + 0.007 * noise[i]
            })
            .collect()
    }

    /// 1 kHz 70 dB tone pulse, on during 0.25..0.75 s.
    pub fn pulse() -> Vec<f64> {
        let a = amp(70.0);
        (0..LEN)
            .map(|i| {
                let t = i as f64 / SR;
                if (0.25..0.75).contains(&t) {
                    a * fs_math::det::sin(2.0 * core::f64::consts::PI * 1000.0 * t)
                } else {
                    0.0
                }
            })
            .collect()
    }
}

/// Provenance tool, not a gate: dumps the fixture signals as raw
/// little-endian f64 for the compiled ISO reference binary
/// (`FS_PSYCHO_DUMP_DIR` names the target directory). The pinned
/// values in the signal-path tests were produced by running that
/// binary on exactly these bytes.
#[test]
#[ignore = "provenance tool: set FS_PSYCHO_DUMP_DIR and run explicitly"]
fn dump_reference_signals() {
    let dir = std::env::var("FS_PSYCHO_DUMP_DIR").expect("FS_PSYCHO_DUMP_DIR");
    let write = |name: &str, data: &[f64]| {
        let bytes: Vec<u8> = data.iter().flat_map(|v| v.to_le_bytes()).collect();
        std::fs::write(format!("{dir}/{name}.f64"), bytes).expect("write");
    };
    write("tone1k60", &sigfix::tone(1000.0, 60.0));
    write("tone250_80", &sigfix::tone(250.0, 80.0));
    write("pulse1k70", &sigfix::pulse());
    write("tonal_a", &sigfix::tonal_a());
    write("noise_b", &sigfix::noise_b());
    write("lowsnr_c", &sigfix::lowsnr_c());
}

/// Stationary loudness through the 48 kHz PCM filterbank path,
/// EXACTNESS-pinned against the compiled ISO reference binary
/// (harness over `f_loudness_from_signal`, LoudnessMethodStationary,
/// time_skip 0.5 s) on bit-identical fixture PCM (see
/// `dump_reference_signals`). This path also lands on the Annex B.3
/// PUBLISHED values (4.019 / 14.655 sone) — evidence the level-vector
/// path's 13% tone gap is filterbank leakage, not a port bug.
#[test]
// Pins are the reference binary's %.17g output verbatim; some carry
// one digit beyond f64 (clippy excessive_precision) — kept verbatim
// so provenance diffs are textual.
#[allow(clippy::excessive_precision)]
fn signal_path_stationary_matches_reference() {
    use fs_psycho::signal::loudness_stationary_from_pcm;
    let cases: [(&str, Vec<f64>, f64, f64); 2] = [
        (
            "1k60",
            sigfix::tone(1000.0, 60.0),
            4.017_710_193_205_207_5,
            4.019,
        ),
        (
            "250_80",
            sigfix::tone(250.0, 80.0),
            14.654_231_012_082_258,
            14.655,
        ),
    ];
    for (name, pcm, reference, published) in cases {
        let out = loudness_stationary_from_pcm(&pcm, sigfix::SR, 0.5, SoundField::Free)
            .expect("loudness");
        assert!(
            ((out.sones - reference) / reference).abs() < 1.0e-9,
            "{name}: {:.17} vs reference {reference:.17}",
            out.sones
        );
        assert!(
            (out.sones - published).abs() < 0.05,
            "{name}: {} vs Annex B.3 published {published}",
            out.sones
        );
    }
    println!("{{\"suite\":\"fs-psycho\",\"case\":\"signal-stationary\",\"verdict\":\"pass\"}}");
}

/// Time-varying loudness of a steady 1 kHz 60 dB tone: monotone rise
/// to a plateau that matches the reference run exactly (Nmax and N5
/// pins plus probe frames), and the plateau sits just above the
/// stationary value (temporal weighting settles to it).
#[test]
// Pins are the reference binary's %.17g output verbatim; some carry
// one digit beyond f64 (clippy excessive_precision) — kept verbatim
// so provenance diffs are textual.
#[allow(clippy::excessive_precision)]
fn signal_path_time_varying_steady_tone() {
    use fs_psycho::signal::loudness_time_varying;
    let out = loudness_time_varying(&sigfix::tone(1000.0, 60.0), sigfix::SR, SoundField::Free)
        .expect("tv");
    assert_eq!(out.loudness.len(), 3000);
    // Reference binary pins (bit-identical input, full precision).
    let pins = [
        (187usize, 3.488_182_192_159_547_5),
        (748, 4.009_182_414_423_507_1),
        (1496, 4.018_800_724_193_021_7),
        (2992, 4.018_846_940_720_695_9),
    ];
    for (frame, want) in pins {
        let got = out.loudness[frame];
        assert!(
            ((got - want) / want).abs() < 1.0e-9,
            "frame {frame}: {got:.17} vs {want:.17}"
        );
    }
    let nmax_ref = 4.018_846_940_772_244;
    let n5_ref = 4.018_846_938_852_739_5;
    assert!(((out.n_max - nmax_ref) / nmax_ref).abs() < 1.0e-9);
    assert!(((out.n5 - n5_ref) / n5_ref).abs() < 1.0e-9);
    // Behavioral: monotone rise onto the plateau (sampled coarsely).
    for w in out.loudness.chunks(300).collect::<Vec<_>>().windows(2) {
        assert!(w[1][0] >= w[0][0] - 1e-9, "rise must be monotone");
    }
    println!("{{\"suite\":\"fs-psycho\",\"case\":\"signal-tv-steady\",\"verdict\":\"pass\"}}");
}

/// Time-varying loudness of a 500 ms 1 kHz 70 dB tone pulse:
/// reference-exact Nmax/N5 and probe frames through rise and decay,
/// plus the temporal asymmetry the nonlinear decay exists to model
/// (loudness persists after offset: still >15% of max 90 ms later,
/// while the pre-onset region is silent).
#[test]
// Pins are the reference binary's %.17g output verbatim; some carry
// one digit beyond f64 (clippy excessive_precision) — kept verbatim
// so provenance diffs are textual.
#[allow(clippy::excessive_precision)]
fn signal_path_time_varying_pulse_decay() {
    use fs_psycho::signal::loudness_time_varying;
    let out = loudness_time_varying(&sigfix::pulse(), sigfix::SR, SoundField::Free).expect("tv");
    let nmax_ref = 8.808_581_173_222_430_4;
    let n5_ref = 8.077_504_421_840_171;
    assert!(((out.n_max - nmax_ref) / nmax_ref).abs() < 1.0e-9);
    assert!(((out.n5 - n5_ref) / n5_ref).abs() < 1.0e-9);
    // The decay TAIL carries a one-time branch-flip offset: near the
    // pulse offset a sub-step of the nonlinear decay element flips a
    // case boundary (its |Ui - UoLast| < 1e-5 equality band) under
    // ulp-level libm differences. Full-series diff vs the reference
    // (review-measured): divergence starts ~frame 1519, peaks at
    // -7.3e-9 relative around frames 1539-1550, then settles to a
    // constant -2.9e-9 for the whole tail (zero frames exceed 1e-8;
    // the steady tone holds 9.3e-15 everywhere). Tail gates are
    // therefore 1e-7; rise/plateau pins stay at 1e-9.
    let pins = [
        (561usize, 5.769_132_334_749_696_3, 1.0e-9),
        (1496, 8.082_105_897_289_144_6, 1.0e-9),
        (1683, 1.699_043_551_606_589_4, 1.0e-7),
        (2244, 0.030_630_403_130_731_315, 1.0e-7),
    ];
    for (frame, want, tol) in pins {
        let got = out.loudness[frame];
        assert!(
            ((got - want) / want).abs() < tol,
            "frame {frame}: {got:.17} vs {want:.17}"
        );
    }
    // Silence before onset; persistence after offset.
    assert!(out.loudness[300] < 1e-6, "pre-onset must be silent");
    assert!(
        out.loudness[1683] > 0.15 * out.n_max,
        "post-offset loudness must persist (nonlinear decay)"
    );
    println!("{{\"suite\":\"fs-psycho\",\"case\":\"signal-tv-pulse\",\"verdict\":\"pass\"}}");
}

/// The verified phon conversion (the reference's own f_sone_to_phon)
/// and the signal-path refusals, all typed by name.
#[test]
fn phon_conversion_and_signal_refusals() {
    use fs_psycho::signal::{loudness_stationary_from_pcm, loudness_time_varying, phon_from_sone};
    // Exact by construction: 1 sone = 40 phon, each doubling +10.
    assert!((phon_from_sone(1.0).expect("phon") - 40.0).abs() < 1e-12);
    assert!((phon_from_sone(2.0).expect("phon") - 50.0).abs() < 1e-12);
    assert!((phon_from_sone(4.0).expect("phon") - 60.0).abs() < 1e-12);
    // Below 1 sone: the reference's 40 (N + 0.0005)^0.35 branch with
    // its 3-phon floor.
    let half = phon_from_sone(0.5).expect("phon");
    assert!(
        (half - 40.0 * fs_math::det::pow(0.5005, 0.35)).abs() < 1e-12,
        "sub-sone branch: {half}"
    );
    assert!((phon_from_sone(0.0).expect("phon") - 3.0).abs() < 1e-12);
    assert!(phon_from_sone(f64::NAN).is_err());
    assert!(matches!(
        phon_from_sone(-0.1),
        Err(PsychoError::DegenerateSignal { .. })
    ));
    // Signal-path refusals.
    let good = sigfix::tone(1000.0, 60.0);
    assert!(matches!(
        loudness_time_varying(&good, 44_100.0, SoundField::Free),
        Err(PsychoError::UnsupportedRate { .. })
    ));
    let mut bad = good.clone();
    bad[5] = f64::INFINITY;
    assert!(matches!(
        loudness_time_varying(&bad, sigfix::SR, SoundField::Free),
        Err(PsychoError::NonFinite { .. })
    ));
    // Huge FINITE samples overflow the squaring stage — must refuse,
    // not return Ok(inf/NaN) (review-executed hole, closed).
    let huge = vec![1.0e200f64; sigfix::LEN];
    assert!(matches!(
        loudness_time_varying(&huge, sigfix::SR, SoundField::Free),
        Err(PsychoError::NonFinite { .. })
    ));
    assert!(matches!(
        loudness_time_varying(&good[..40], sigfix::SR, SoundField::Free),
        Err(PsychoError::DegenerateSignal { .. })
    ));
    assert!(matches!(
        loudness_stationary_from_pcm(&good, sigfix::SR, 1.5, SoundField::Free),
        Err(PsychoError::DegenerateSignal { .. })
    ));
    assert!(matches!(
        loudness_stationary_from_pcm(&good[..100], sigfix::SR, 0.5, SoundField::Free),
        Err(PsychoError::DegenerateSignal { .. })
    ));
}

/// The batch Pareto API is aggregation-exact: every field is BITWISE
/// equal to the corresponding standalone call on the same input (the
/// wiring-mistake contract), refusals propagate typed, and the
/// roughness mean covers exactly the whole blocks.
#[test]
fn pareto_batch_is_aggregation_exact() {
    use fs_psycho::roughness::{DW_BLOCK, roughness_dw_block};
    use fs_psycho::signal::{
        loudness_stationary_from_pcm, loudness_time_varying, pareto_metrics, phon_from_sone,
    };
    // An AM tone so roughness, sharpness, and attack all have signal.
    let pcm: Vec<f64> = (0..sigfix::LEN)
        .map(|i| {
            let t = i as f64 / sigfix::SR;
            let a = fs_math::det::sqrt(2.0) * 2.0e-5 * fs_math::det::pow(10.0, 65.0 / 20.0);
            a * (1.0 + fs_math::det::sin(2.0 * core::f64::consts::PI * 70.0 * t))
                * fs_math::det::sin(2.0 * core::f64::consts::PI * 1000.0 * t)
        })
        .collect();
    let m = pareto_metrics(&pcm, sigfix::SR, 0.5, SoundField::Free, 480).expect("batch");
    let stationary =
        loudness_stationary_from_pcm(&pcm, sigfix::SR, 0.5, SoundField::Free).expect("stat");
    let tv = loudness_time_varying(&pcm, sigfix::SR, SoundField::Free).expect("tv");
    assert_eq!(m.sones_stationary.to_bits(), stationary.sones.to_bits());
    assert_eq!(
        m.phon_stationary.to_bits(),
        phon_from_sone(stationary.sones).expect("phon").to_bits()
    );
    assert_eq!(m.n5.to_bits(), tv.n5.to_bits());
    assert_eq!(m.n_max.to_bits(), tv.n_max.to_bits());
    assert_eq!(
        m.sharpness_acum.to_bits(),
        sharpness_din(&stationary.specific).expect("s").to_bits()
    );
    assert_eq!(m.roughness_blocks, sigfix::LEN / DW_BLOCK);
    let mut r_sum = 0.0;
    for b in 0..m.roughness_blocks {
        r_sum += roughness_dw_block(&pcm[b * DW_BLOCK..(b + 1) * DW_BLOCK], sigfix::SR)
            .expect("roughness");
    }
    assert_eq!(
        m.roughness_asper_mean.to_bits(),
        (r_sum / m.roughness_blocks as f64).to_bits()
    );
    assert_eq!(
        m.log_attack_time.to_bits(),
        fs_psycho::log_attack_time(&pcm, sigfix::SR, 480)
            .expect("lat")
            .to_bits()
    );
    // Sanity: this fixture's roughness is strong (70 Hz AM near the
    // anchor) — the batch value must reflect it, not a stub.
    assert!(
        m.roughness_asper_mean > 0.5,
        "AM fixture roughness {}",
        m.roughness_asper_mean
    );
    // Refusal propagation: too short for a roughness block.
    assert!(matches!(
        pareto_metrics(&pcm[..DW_BLOCK - 1], sigfix::SR, 0.0, SoundField::Free, 48),
        Err(PsychoError::DegenerateSignal { .. })
    ));
    // Refusal propagation: unsupported rate comes from the loudness
    // path before anything else runs.
    assert!(matches!(
        pareto_metrics(&pcm, 44_100.0, 0.5, SoundField::Free, 480),
        Err(PsychoError::UnsupportedRate { .. })
    ));
    println!("{{\"suite\":\"fs-psycho\",\"case\":\"pareto-batch\",\"verdict\":\"pass\"}}");
}

/// ECMA tonality (TNR + PR), EXACTNESS-pinned against the extracted
/// Apache-2.0 MoSQITo reference run on bit-identical fixture PCM
/// (dump_reference_signals + tnrref/drive_tnr.py, full precision).
#[test]
// float_cmp: the zero totals are the reference's EXACT no-prominent
// convention, not a computed value near zero.
#[allow(clippy::excessive_precision, clippy::float_cmp)]
fn tonality_matches_reference() {
    use fs_psycho::tonality::{prominence_ratio_ecma, tone_to_noise_ecma};
    let sr = sigfix::SR;
    // Fixture A: two deterministic tones over LCG noise.
    let a = sigfix::tonal_a();
    let tnr = tone_to_noise_ecma(&a, sr).expect("tnr");
    assert_eq!(tnr.tones.len(), 2, "{:?}", tnr.tones);
    let pins_tnr: [(f64, f64); 2] = [
        (1_000.488_281_25, 12.434_966_853_609_048),
        (3_000.732_421_875, 51.331_313_916_320_724),
    ];
    for (tone, (f_ref, v_ref)) in tnr.tones.iter().zip(pins_tnr) {
        assert_eq!(tone.frequency_hz.to_bits(), f_ref.to_bits(), "tone freq");
        assert!(
            ((tone.ratio_db - v_ref) / v_ref).abs() < 1.0e-9,
            "TNR {:.17} vs {v_ref:.17}",
            tone.ratio_db
        );
        assert!(tone.prominent);
    }
    let t_tnr_ref = 51.331_873_83;
    assert!(
        (tnr.total_db - t_tnr_ref).abs() < 1.0e-6,
        "{}",
        tnr.total_db
    );
    let pr = prominence_ratio_ecma(&a, sr).expect("pr");
    assert_eq!(pr.tones.len(), 2);
    let pins_pr: [(f64, f64); 2] = [
        (1_000.488_281_25, 61.485_567_666_731_83),
        (3_000.732_421_875, 50.696_977_402_016_664),
    ];
    for (tone, (f_ref, v_ref)) in pr.tones.iter().zip(pins_pr) {
        assert_eq!(tone.frequency_hz.to_bits(), f_ref.to_bits());
        assert!(
            ((tone.ratio_db - v_ref) / v_ref).abs() < 1.0e-9,
            "PR {:.17} vs {v_ref:.17}",
            tone.ratio_db
        );
        assert!(tone.prominent);
    }
    assert!(
        (pr.total_db - 61.833_436_68).abs() < 1.0e-6,
        "{}",
        pr.total_db
    );
    // Fixture B: noise only — TNR finds NO tones; PR finds exactly 10
    // weak candidates, none prominent, total 0 (reference-exact).
    let b = sigfix::noise_b();
    let tnr_b = tone_to_noise_ecma(&b, sr).expect("tnr");
    assert!(tnr_b.tones.is_empty(), "{:?}", tnr_b.tones);
    assert_eq!(tnr_b.total_db, 0.0);
    let pr_b = prominence_ratio_ecma(&b, sr).expect("pr");
    assert_eq!(pr_b.tones.len(), 10, "{:?}", pr_b.tones);
    assert!(pr_b.tones.iter().all(|t| !t.prominent));
    assert_eq!(pr_b.total_db, 0.0);
    // One PR candidate pinned mid-list (953.6 Hz, 0.5146 dB).
    let probe = &pr_b.tones[2];
    assert_eq!(probe.frequency_hz.to_bits(), 953.613_281_25_f64.to_bits());
    assert!(
        ((probe.ratio_db - 0.514_649_934_651_174_6) / 0.514_649_934_651_174_6).abs() < 1.0e-6,
        "{:.17}",
        probe.ratio_db
    );
    // Fixture C: low-SNR 500 Hz tone — detected, prominent, pinned.
    let c = sigfix::lowsnr_c();
    let tnr_c = tone_to_noise_ecma(&c, sr).expect("tnr");
    assert_eq!(tnr_c.tones.len(), 1);
    assert_eq!(
        tnr_c.tones[0].frequency_hz.to_bits(),
        500.976_562_5_f64.to_bits()
    );
    assert!(
        ((tnr_c.tones[0].ratio_db - 12.384_273_273_044_364) / 12.384_273_273_044_364).abs()
            < 1.0e-9
    );
    assert!(tnr_c.tones[0].prominent);
    let pr_c = prominence_ratio_ecma(&c, sr).expect("pr");
    assert_eq!(pr_c.tones.len(), 6, "{:?}", pr_c.tones);
    assert!(
        ((pr_c.tones[0].ratio_db - 32.775_056_728_327_13) / 32.775_056_728_327_13).abs() < 1.0e-9
    );
    assert!(pr_c.tones[0].prominent);
    assert!(pr_c.tones[1..].iter().all(|t| !t.prominent));
    println!("{{\"suite\":\"fs-psycho\",\"case\":\"tonality\",\"verdict\":\"pass\"}}");
}

/// Tonality refusals, typed by name.
#[test]
fn tonality_refusals_are_typed() {
    use fs_psycho::tonality::tone_to_noise_ecma;
    let good = sigfix::tonal_a();
    // Non-power-of-two length.
    assert!(matches!(
        tone_to_noise_ecma(&good[..60_000], sigfix::SR),
        Err(PsychoError::Shape { .. })
    ));
    // Too short.
    assert!(matches!(
        tone_to_noise_ecma(&good[..2048], sigfix::SR),
        Err(PsychoError::Shape { .. })
    ));
    // NaN sample.
    let mut bad = good.clone();
    bad[9] = f64::NAN;
    assert!(matches!(
        tone_to_noise_ecma(&bad, sigfix::SR),
        Err(PsychoError::NonFinite { .. })
    ));
    // Bad rates.
    for sr in [f64::NAN, 0.0, -1.0, f64::INFINITY] {
        assert!(matches!(
            tone_to_noise_ecma(&good, sr),
            Err(PsychoError::NonFinite { .. })
        ));
    }
    // A rate so low the detection band is unresolvable.
    assert!(matches!(
        tone_to_noise_ecma(&good[..4096], 100.0),
        Err(PsychoError::DegenerateSignal { .. })
    ));
}
