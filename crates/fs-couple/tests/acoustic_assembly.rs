//! Description → waveform: prestressed string + viscothermal duct,
//! gas/material parameter motion, deterministic PCM16, typed refusals.

use fs_couple::acoustic_realize::{AcousticRealizeError, assembly_wav, realize_assembly};
use fs_couple::pcm_wav::{WavError, encode_pcm16_wav};
use fs_scenario::{
    AcousticAssembly, AmbientGas, BeatingReed, BowStroke, CylinderSegment, Listener, Pluck,
    PrestressedString, RadiatingPlate, ToneHole, ViscothermalDuct, VolumeVelocityPulse,
    WaveguideEnd,
};

fn empty_base() -> AcousticAssembly {
    AcousticAssembly {
        ambient: AmbientGas::sea_level(),
        string: None,
        duct: None,
        pluck: None,
        bow: None,
        blow: None,
        reed: None,
        soundboard: None,
        listener: Listener { distance_m: 1.0 },
        sample_rate_hz: 8_000,
        duration_s: 0.06,
    }
}

fn nylon_like(tension_n: f64, lin_density_kg_m: f64) -> PrestressedString {
    PrestressedString {
        length_m: 0.65,
        tension_n,
        lin_density_kg_m,
        axial_stiffness_n: 0.0,
        width_m: 0.01,
        n_modes: 4,
        damping_ratio: 0.004,
        rayleigh: None,
    }
}

fn plucked(tension_n: f64, lin_density_kg_m: f64, height_m: f64) -> AcousticAssembly {
    AcousticAssembly {
        string: Some(nylon_like(tension_n, lin_density_kg_m)),
        pluck: Some(Pluck {
            station_frac: 0.25,
            height_m,
        }),
        ..empty_base()
    }
}

fn closed_duct(temperature_k: f64) -> AcousticAssembly {
    AcousticAssembly {
        ambient: AmbientGas {
            temperature_k,
            pressure_pa: 101_325.0,
        },
        duct: Some(ViscothermalDuct {
            segments: vec![CylinderSegment {
                radius_m: 0.012,
                length_m: 0.34,
            }],
            tone_holes: vec![],
            termination: WaveguideEnd::Closed,
        }),
        blow: Some(VolumeVelocityPulse {
            peak_m3_s: 2.0e-5,
            duration_s: 0.002,
        }),
        duration_s: 0.08,
        ..empty_base()
    }
}

/// Mean period from interpolated falling zero crossings, in samples.
fn zero_cross_period(x: &[f64]) -> f64 {
    let mut prev = x[0];
    let mut times = Vec::new();
    for (i, &s) in x.iter().enumerate().skip(1) {
        if prev > 0.0 && s <= 0.0 {
            let frac = prev / (prev - s);
            times.push(i as f64 - 1.0 + frac);
        }
        prev = s;
    }
    assert!(times.len() >= 3, "need crossings, got {}", times.len());
    (times[times.len() - 1] - times[0]) / (times.len() - 1) as f64
}

/// Lag of the strongest positive autocorrelation in `[min_lag, max_lag]`.
fn dominant_period_samples(x: &[f64], min_lag: usize, max_lag: usize) -> usize {
    let mut best_lag = min_lag;
    let mut best = f64::NEG_INFINITY;
    for lag in min_lag..=max_lag {
        let mut acc = 0.0;
        for i in 0..x.len() - lag {
            acc += x[i] * x[i + lag];
        }
        if acc > best {
            best = acc;
            best_lag = lag;
        }
    }
    best_lag
}

fn peak_abs(x: &[f64]) -> f64 {
    x.iter().fold(0.0_f64, |m, &v| m.max(v.abs()))
}

#[test]
fn plucked_string_fundamental_moves_with_mu_and_tension() {
    let base = realize_assembly(&plucked(80.0, 0.006, 0.003)).expect("base");
    let heavy = realize_assembly(&plucked(80.0, 0.012, 0.003)).expect("heavy");
    let taut = realize_assembly(&plucked(160.0, 0.006, 0.003)).expect("taut");
    let p_base = dominant_period_samples(&base.pressure_pa, 40, 200);
    let p_heavy = dominant_period_samples(&heavy.pressure_pa, 40, 220);
    let p_taut = dominant_period_samples(&taut.pressure_pa, 30, 160);
    let mu_ratio = p_heavy as f64 / p_base as f64;
    let t_ratio = p_base as f64 / p_taut as f64;
    assert!(
        (mu_ratio - 2.0_f64.sqrt()).abs() < 0.08,
        "doubling μ should lengthen the period by √2, got {mu_ratio} ({p_heavy}/{p_base})"
    );
    assert!(
        (t_ratio - 2.0_f64.sqrt()).abs() < 0.08,
        "doubling T0 should shorten the period by √2, got {t_ratio} ({p_base}/{p_taut})"
    );
}

#[test]
fn plucked_height_scales_pressure_and_does_not_peak_normalize() {
    let soft = realize_assembly(&plucked(80.0, 0.006, 0.002)).expect("soft");
    let hard = realize_assembly(&plucked(80.0, 0.006, 0.004)).expect("hard");
    let a = peak_abs(&soft.pressure_pa);
    let b = peak_abs(&hard.pressure_pa);
    assert!(a > 0.0 && b > 0.0, "radiated pressure must be live");
    assert!(
        (b / a - 2.0).abs() < 0.05,
        "linear radiation must track pluck height, got {} / {}",
        b,
        a
    );
}

#[test]
fn duct_period_tracks_sound_speed_from_ambient_temperature() {
    let cold = realize_assembly(&closed_duct(288.15)).expect("cold");
    let hot = realize_assembly(&closed_duct(330.0)).expect("hot");
    assert!(
        (hot.gas.sound_speed - cold.gas.sound_speed).abs() > 10.0,
        "GasState must move c with T"
    );
    let p_cold = dominant_period_samples(&cold.pressure_pa, 8, 40);
    let p_hot = dominant_period_samples(&hot.pressure_pa, 8, 40);
    let measured = p_hot as f64 / p_cold as f64;
    let expected = cold.gas.sound_speed / hot.gas.sound_speed;
    assert!(
        (measured - expected).abs() < 0.08,
        "closed-pipe period should track 1/c(T): measured {measured} expected {expected} ({p_hot}/{p_cold})"
    );
}

#[test]
fn pcm16_wav_is_deterministic_and_not_peak_normalized() {
    let realized = realize_assembly(&plucked(80.0, 0.006, 0.003)).expect("pluck");
    let (a, clips_a) = assembly_wav(&realized, 50.0).expect("wav a");
    let (b, clips_b) = assembly_wav(&realized, 50.0).expect("wav b");
    assert_eq!(a, b);
    assert_eq!(clips_a, clips_b);
    assert_eq!(&a[0..4], b"RIFF");
    assert_eq!(&a[8..12], b"WAVE");
    assert_eq!(&a[12..16], b"fmt ");
    assert_eq!(&a[36..40], b"data");
    let rate = u32::from_le_bytes(a[24..28].try_into().expect("rate"));
    assert_eq!(rate, realized.sample_rate_hz);
    let bits = u16::from_le_bytes(a[34..36].try_into().expect("bits"));
    assert_eq!(bits, 16);
    let (wide, _) = encode_pcm16_wav(&realized.pressure_pa, realized.sample_rate_hz, 100.0)
        .expect("wide scale");
    let (narrow, _) = encode_pcm16_wav(&realized.pressure_pa, realized.sample_rate_hz, 50.0)
        .expect("narrow scale");
    assert_ne!(
        wide, narrow,
        "a smaller full-scale must change PCM codes; peak-normalize would hide it"
    );
}

#[test]
fn empty_assembly_and_bad_events_refuse() {
    let err = realize_assembly(&empty_base()).expect_err("empty");
    assert!(matches!(
        err,
        AcousticRealizeError::InvalidDescription { what } if what.contains("neither")
    ));

    let mut no_pluck = plucked(80.0, 0.006, 0.003);
    no_pluck.pluck = None;
    assert!(matches!(
        realize_assembly(&no_pluck),
        Err(AcousticRealizeError::InvalidDescription { what }) if what.contains("pluck")
    ));

    let mut bad_station = plucked(80.0, 0.006, 0.003);
    bad_station.pluck = Some(Pluck {
        station_frac: 0.0,
        height_m: 0.003,
    });
    assert!(matches!(
        realize_assembly(&bad_station),
        Err(AcousticRealizeError::InvalidDescription { what }) if what.contains("station")
    ));

    let mut no_blow = closed_duct(288.15);
    no_blow.blow = None;
    assert!(matches!(
        realize_assembly(&no_blow),
        Err(AcousticRealizeError::InvalidDescription { what }) if what.contains("volume-velocity")
    ));
}

#[test]
fn open_duct_refuses_nyquist_ka_above_the_radiation_fit() {
    let mut open = closed_duct(288.15);
    open.sample_rate_hz = 44_100;
    if let Some(duct) = open.duct.as_mut() {
        duct.termination = WaveguideEnd::UnflangedOpen;
    }
    let err = realize_assembly(&open).expect_err("ka");
    assert!(
        matches!(
            err,
            AcousticRealizeError::Duct(fs_duct::DuctError::RadiationKaTooLarge { .. })
        ),
        "got {err:?}"
    );
}

#[test]
fn empty_or_invalid_wav_refuses() {
    assert!(matches!(
        encode_pcm16_wav(&[], 8_000, 1.0),
        Err(WavError::InvalidInput { .. })
    ));
    assert!(matches!(
        encode_pcm16_wav(&[0.1], 0, 1.0),
        Err(WavError::InvalidInput { .. })
    ));
    assert!(matches!(
        encode_pcm16_wav(&[0.1], 8_000, 0.0),
        Err(WavError::InvalidInput { .. })
    ));
}

#[test]
fn kirchhoff_carrier_loud_pluck_raises_pitch() {
    let mut soft = plucked(70.0, 0.005, 4.0e-4);
    if let Some(s) = soft.string.as_mut() {
        s.axial_stiffness_n = 4.0e4;
        s.n_modes = 3;
    }
    let mut loud = soft.clone();
    loud.pluck = Some(Pluck {
        station_frac: 0.25,
        height_m: 3.5e-3,
    });
    let p_soft = realize_assembly(&soft).expect("soft KC");
    let p_loud = realize_assembly(&loud).expect("loud KC");
    let t_soft = zero_cross_period(&p_soft.pressure_pa);
    let t_loud = zero_cross_period(&p_loud.pressure_pa);
    assert!(
        t_loud < t_soft * 0.997,
        "KC pitch glide must raise f0 at larger amplitude ({t_loud:.3} vs {t_soft:.3} samples)"
    );
}

#[test]
fn open_tone_hole_shortens_bore_period() {
    let mut closed = empty_base();
    closed.duration_s = 0.08;
    closed.blow = Some(VolumeVelocityPulse {
        peak_m3_s: 2.0e-5,
        duration_s: 0.002,
    });
    closed.duct = Some(ViscothermalDuct {
        segments: vec![
            CylinderSegment {
                radius_m: 0.012,
                length_m: 0.08,
            },
            CylinderSegment {
                radius_m: 0.012,
                length_m: 0.26,
            },
        ],
        tone_holes: vec![],
        termination: WaveguideEnd::Closed,
    });
    let mut vented = closed.clone();
    if let Some(duct) = vented.duct.as_mut() {
        duct.tone_holes = vec![ToneHole {
            after_segment: 0,
            radius_m: 0.003,
            chimney_m: 0.003,
            open: true,
        }];
    }
    let a = realize_assembly(&closed).expect("plain");
    let b = realize_assembly(&vented).expect("hole");
    let ta = zero_cross_period(&a.pressure_pa);
    let tb = zero_cross_period(&b.pressure_pa);
    assert!(
        tb < ta * 0.98,
        "an open tone hole must raise the ringing frequency ({tb:.2} vs {ta:.2})"
    );
}

#[test]
fn beating_reed_locks_near_the_quarter_wave() {
    let mut a = empty_base();
    a.duration_s = 0.10;
    a.duct = Some(ViscothermalDuct {
        segments: vec![CylinderSegment {
            radius_m: 0.0075,
            length_m: 0.50,
        }],
        tone_holes: vec![],
        termination: WaveguideEnd::UnflangedOpen,
    });
    a.reed = Some(BeatingReed {
        rest_opening_m: 4.0e-4,
        width_m: 0.013,
        closing_pressure_pa: 6_000.0,
        blowing_pressure_pa: 2_000.0,
        attack_s: 0.008,
    });
    let out = realize_assembly(&a).expect("reed");
    let tail = &out.pressure_pa[out.pressure_pa.len() / 2..];
    let mean = tail.iter().sum::<f64>() / tail.len() as f64;
    let ac: Vec<f64> = tail.iter().map(|p| p - mean).collect();
    let rms: f64 = (ac.iter().map(|p| p * p).sum::<f64>() / ac.len() as f64).sqrt();
    assert!(rms > 1.0, "reed-bore must self-oscillate, rms={rms}");
    let period = zero_cross_period(&ac);
    assert!(
        (4.0..120.0).contains(&period),
        "reed-bore must be periodic, period={period:.2} samples"
    );
    let mut silent = a.clone();
    if let Some(reed) = silent.reed.as_mut() {
        reed.blowing_pressure_pa = 0.0;
    }
    let quiet = realize_assembly(&silent).expect("silent reed");
    assert!(
        peak_abs(&quiet.pressure_pa) < 0.05 * peak_abs(&out.pressure_pa),
        "zero blowing pressure must not speak"
    );
}

#[test]
fn soundboard_adds_body_radiation() {
    let bare = plucked(80.0, 0.006, 0.003);
    let mut body = bare.clone();
    body.soundboard = Some(RadiatingPlate {
        area_m2: 0.12,
        mass_kg: 0.18,
        frequency_hz: 110.0,
        damping_ratio: 0.03,
    });
    let a = realize_assembly(&bare).expect("bare");
    let b = realize_assembly(&body).expect("body");
    assert!(
        peak_abs(&b.pressure_pa) > peak_abs(&a.pressure_pa) * 1.05,
        "a driven soundboard must add observer pressure"
    );
}

#[test]
fn bow_stroke_is_a_live_excitation() {
    let mut a = empty_base();
    a.string = Some(nylon_like(80.0, 0.006));
    a.bow = Some(BowStroke {
        station_frac: 0.12,
        normal_force_n: 0.8,
        velocity_m_s: 0.4,
        mu_static: 0.8,
        mu_dynamic: 0.3,
        stribeck_m_s: 0.05,
    });
    let out = realize_assembly(&a).expect("bow");
    assert!(
        peak_abs(&out.pressure_pa) > 1.0e-4,
        "a bow stroke must drive the string"
    );
}
