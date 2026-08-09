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
    // reference (verified against the compiled reference binary:
    // 83.295660), so the gate is EXACTNESS at 0.001 sone — hiding a
    // port bug inside the standard's 5% compliance tolerance is
    // exactly what caught the -.25f transcription mangling.
    let reference = 83.29566;
    assert!(
        (out.sones - reference).abs() <= 1.0e-3,
        "test signal 1: {} sone vs reference {reference}",
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
