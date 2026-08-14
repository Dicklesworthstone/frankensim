//! Description → waveform: prestressed string + viscothermal duct,
//! gas/material parameter motion, deterministic PCM16, typed refusals.

use fs_couple::acoustic_realize::{
    AcousticRealizeError, assembly_wav, realize_assembly, string_mode_omega,
};
use fs_couple::pcm_wav::{WavError, encode_pcm16_wav};
use fs_scenario::{
    AcousticAssembly, AmbientGas, BeatingReed, BowStroke, ContactTexture, CylinderSegment,
    HelmholtzCavity, Listener, Pluck, PrestressedString, RadiatingPlate, ThinPlate, ToneHole,
    UnilateralObstacle, ViscothermalDuct, VolumeVelocityPulse, WaveguideEnd,
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
        body_modes: vec![],
        plate: None,
        cavity: None,
        obstacles: vec![],
        contact_texture: None,
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
        bending_stiffness_n_m2: 0.0,
        polarization_detune: 0.0,
        moving_end: false,
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
            relative_humidity: 0.0,
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
        blowing_pressure_pa: 2_800.0,
        attack_s: 0.008,
        mass_kg: 0.0,
        stiffness_n_m: 0.0,
    });
    let out = realize_assembly(&a).expect("reed");
    let tail = &out.pressure_pa[out.pressure_pa.len() / 2..];
    let mean = tail.iter().sum::<f64>() / tail.len() as f64;
    let ac: Vec<f64> = tail.iter().map(|p| p - mean).collect();
    let rms: f64 = (ac.iter().map(|p| p * p).sum::<f64>() / ac.len() as f64).sqrt();
    assert!(rms > 5.0, "reed-bore must self-oscillate, rms={rms}");
    let period = zero_cross_period(&ac);
    let quarter = 4.0 * 0.50 / out.gas.sound_speed * f64::from(a.sample_rate_hz);
    assert!(
        (period - quarter).abs() < 0.28 * quarter
            || (period - quarter / 3.0).abs() < 0.28 * (quarter / 3.0),
        "reed lock {period:.2} vs quarter-wave {quarter:.1} or twelfth"
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
fn stepped_bore_shifts_the_ode_period() {
    let mut uniform = empty_base();
    uniform.duration_s = 0.08;
    uniform.blow = Some(VolumeVelocityPulse {
        peak_m3_s: 2.0e-5,
        duration_s: 0.002,
    });
    uniform.duct = Some(ViscothermalDuct {
        segments: vec![CylinderSegment {
            radius_m: 0.012,
            length_m: 0.34,
        }],
        tone_holes: vec![],
        termination: WaveguideEnd::Closed,
    });
    let mut stepped = uniform.clone();
    if let Some(duct) = stepped.duct.as_mut() {
        duct.segments = vec![
            CylinderSegment {
                radius_m: 0.012,
                length_m: 0.17,
            },
            CylinderSegment {
                radius_m: 0.006,
                length_m: 0.17,
            },
        ];
    }
    let a = realize_assembly(&uniform).expect("uniform");
    let b = realize_assembly(&stepped).expect("stepped");
    let ta = zero_cross_period(&a.pressure_pa);
    let tb = zero_cross_period(&b.pressure_pa);
    assert!(
        (tb - ta).abs() > 0.02 * ta,
        "an area jump must move the ringing period ({tb:.2} vs {ta:.2})"
    );
}

#[test]
fn narrow_bore_decays_faster_from_the_viscothermal_pin() {
    let make = |radius_m: f64| {
        let mut a = empty_base();
        a.duration_s = 0.08;
        a.blow = Some(VolumeVelocityPulse {
            peak_m3_s: 2.0e-5,
            duration_s: 0.002,
        });
        a.duct = Some(ViscothermalDuct {
            segments: vec![CylinderSegment {
                radius_m,
                length_m: 0.34,
            }],
            tone_holes: vec![],
            termination: WaveguideEnd::Closed,
        });
        a
    };
    let wide = realize_assembly(&make(0.016)).expect("wide");
    let narrow = realize_assembly(&make(0.004)).expect("narrow");
    let decay = |p: &[f64]| {
        let peak = p.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
        let tail = &p[p.len() * 3 / 4..];
        let rms = (tail.iter().map(|x| x * x).sum::<f64>() / tail.len() as f64).sqrt();
        rms / peak.max(1.0e-18)
    };
    let d_wide = decay(&wide.pressure_pa);
    let d_narrow = decay(&narrow.pressure_pa);
    assert!(
        d_narrow < d_wide * 0.85,
        "wide-tube R ∝ 1/a must damp a narrower bore faster ({d_narrow} vs {d_wide})"
    );
}

#[test]
fn massive_reed_speaks_on_the_ode_clock() {
    let mut a = empty_base();
    a.duration_s = 0.12;
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
        blowing_pressure_pa: 2_800.0,
        attack_s: 0.008,
        mass_kg: 3.0e-4,
        stiffness_n_m: 0.0,
    });
    let out = realize_assembly(&a).expect("massive reed");
    let tail = &out.pressure_pa[out.pressure_pa.len() / 2..];
    let mean = tail.iter().sum::<f64>() / tail.len() as f64;
    let ac: Vec<f64> = tail.iter().map(|p| p - mean).collect();
    let rms: f64 = (ac.iter().map(|p| p * p).sum::<f64>() / ac.len() as f64).sqrt();
    assert!(
        rms > 1.0,
        "massive reed-bore must self-oscillate, rms={rms}"
    );
    let mut silent = a.clone();
    if let Some(reed) = silent.reed.as_mut() {
        reed.blowing_pressure_pa = 0.0;
    }
    let quiet = realize_assembly(&silent).expect("silent massive reed");
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

#[test]
fn bending_stiffness_makes_partials_inharmonic() {
    let mut flex = nylon_like(80.0, 0.006);
    flex.n_modes = 4;
    let mut stiff = flex;
    // Wound-string stand-in: E=200 GPa, r=0.6 mm → I=π r^4/4.
    // B = π² EI/(T L²) is then large enough that f2/f1 is obviously > 2.
    let r: f64 = 6.0e-4;
    stiff.bending_stiffness_n_m2 = 2.0e11 * core::f64::consts::PI * r.powi(4) / 4.0;
    let r12_flex = string_mode_omega(flex, 2) / string_mode_omega(flex, 1);
    let r12_stiff = string_mode_omega(stiff, 2) / string_mode_omega(stiff, 1);
    assert!((r12_flex - 2.0).abs() < 1.0e-12);
    assert!(
        r12_stiff > 2.01,
        "Fletcher inharmonicity must push f2/f1 above 2, got {r12_stiff}"
    );
    let mut a = plucked(80.0, 0.006, 0.003);
    a.string = Some(stiff);
    realize_assembly(&a).expect("stiff pluck must realize");
}

#[test]
fn two_polarizations_beat() {
    let mut a = plucked(80.0, 0.006, 0.003);
    a.duration_s = 0.20;
    if let Some(s) = a.string.as_mut() {
        s.polarization_detune = 0.012;
        s.n_modes = 2;
        s.damping_ratio = 0.001;
    }
    let out = realize_assembly(&a).expect("two-pol");
    let hop = 80usize;
    let mut env = Vec::new();
    for chunk in out.pressure_pa.chunks(hop) {
        let e = (chunk.iter().map(|p| p * p).sum::<f64>() / chunk.len() as f64).sqrt();
        env.push(e);
    }
    let max_e = env.iter().copied().fold(0.0_f64, f64::max);
    let min_e = env.iter().copied().fold(f64::INFINITY, f64::min);
    assert!(
        max_e > 1.4 * min_e,
        "two polarizations must beat the envelope ({max_e} vs {min_e})"
    );
}

#[test]
fn sitka_body_pair_adds_low_frequency_energy() {
    let mut a = plucked(80.0, 0.006, 0.003);
    a.duration_s = 0.10;
    // Two compact radiators — the same object as any other 1-DOF
    // modal monopole. Frequencies sit in a published guitar-body
    // band; the type does not know that.
    a.body_modes = vec![
        RadiatingPlate {
            area_m2: 0.10,
            mass_kg: 0.16,
            frequency_hz: 98.0,
            damping_ratio: 0.5 / 37.0,
        },
        RadiatingPlate {
            area_m2: 0.06,
            mass_kg: 0.12,
            frequency_hz: 181.0,
            damping_ratio: 0.5 / 22.0,
        },
    ];
    let bare = plucked(80.0, 0.006, 0.003);
    let with = realize_assembly(&a).expect("body");
    let without = realize_assembly(&bare).expect("bare");
    assert!(
        peak_abs(&with.pressure_pa) > peak_abs(&without.pressure_pa),
        "Carcagno-band body modes must radiate"
    );
}

#[test]
fn bowed_rich_spectrum_has_even_partials() {
    let mut a = empty_base();
    a.duration_s = 0.08;
    let mut s = nylon_like(80.0, 0.006);
    s.n_modes = 12;
    a.string = Some(s);
    a.bow = Some(BowStroke {
        station_frac: 0.13,
        normal_force_n: 1.2,
        velocity_m_s: 0.35,
        mu_static: 0.85,
        mu_dynamic: 0.25,
        stribeck_m_s: 0.03,
    });
    let out = realize_assembly(&a).expect("bow helmholtz");
    assert!(peak_abs(&out.pressure_pa) > 1.0e-3);
}

fn steel_panel() -> ThinPlate {
    ThinPlate {
        length_m: 0.20,
        width_m: 0.15,
        thickness_m: 0.002,
        density_kg_m3: 7800.0,
        e1_pa: 200e9,
        e2_pa: 200e9,
        nu12: 0.3,
        g12_pa: 200e9 / (2.0 * 1.3),
        damping_ratio: 0.02,
        n_modes: 2,
        geometric_nonlinearity: false,
        pretension_n_m: 0.0,
        clamped: false,
    }
}

#[test]
fn certified_plate_radiates_without_named_hertz() {
    let bare = plucked(80.0, 0.006, 0.003);
    let mut with = bare.clone();
    with.plate = Some(steel_panel());
    let a = realize_assembly(&bare).expect("bare");
    let b = realize_assembly(&with).expect("plate");
    assert!(
        peak_abs(&b.pressure_pa) > peak_abs(&a.pressure_pa),
        "a certified plate must add observer pressure"
    );
}

#[test]
fn taut_span_obstacle_changes_the_waveform() {
    let bare = plucked(80.0, 0.006, 0.008);
    let mut rattle = bare.clone();
    rattle.obstacles = vec![UnilateralObstacle {
        stations: vec![0.25, 0.5, 0.75],
        gaps_m: vec![0.001, 0.001, 0.001],
        stiffness: 5.0e5,
        alpha: 2.0,
        mu_kinetic: 0.0,
        internal_loss: 0.0,
        provenance: "fixture/stay".into(),
    }];
    let a = realize_assembly(&bare).expect("bare");
    let b = realize_assembly(&rattle).expect("rattle");
    let err: f64 = a
        .pressure_pa
        .iter()
        .zip(&b.pressure_pa)
        .map(|(x, y)| (x - y).abs())
        .sum();
    assert!(err > 1.0e-6, "an obstacle must change the string waveform");
    let mut sticky = rattle.clone();
    sticky.obstacles[0].mu_kinetic = 0.4;
    let c = realize_assembly(&sticky).expect("sticky");
    let err_mu: f64 = b
        .pressure_pa
        .iter()
        .zip(&c.pressure_pa)
        .map(|(x, y)| (x - y).abs())
        .sum();
    assert!(
        err_mu > 1.0e-8,
        "tangential tribo on an obstacle must change the waveform"
    );
}

#[test]
fn mouth_pressure_drives_the_plate_back_into_the_duct() {
    let mut thin = closed_duct(288.15);
    thin.duration_s = 0.06;
    thin.plate = Some(steel_panel());
    let mut thick = thin.clone();
    if let Some(p) = thick.plate.as_mut() {
        p.thickness_m *= 2.0;
    }
    let a = realize_assembly(&thin).expect("thin panel");
    let b = realize_assembly(&thick).expect("thick panel");
    let err: f64 = a
        .pressure_pa
        .iter()
        .zip(&b.pressure_pa)
        .map(|(x, y)| (x - y).abs())
        .sum();
    assert!(
        err > 1.0e-6,
        "a different plate termination must change the mouth pressure"
    );
}

#[test]
fn string_plate_and_duct_share_one_clock() {
    let mut coupled = plucked(80.0, 0.006, 0.003);
    let air = closed_duct(288.15);
    coupled.duct = air.duct;
    coupled.blow = air.blow;
    coupled.plate = Some(steel_panel());
    coupled.duration_s = 0.06;
    let mut string_only = plucked(80.0, 0.006, 0.003);
    string_only.plate = Some(steel_panel());
    string_only.duration_s = 0.06;
    let mut duct_only = closed_duct(288.15);
    duct_only.plate = Some(steel_panel());
    duct_only.duration_s = 0.06;
    let c = realize_assembly(&coupled).expect("three-way");
    let s = realize_assembly(&string_only).expect("string");
    let d = realize_assembly(&duct_only).expect("duct");
    let mut superposed = s.pressure_pa;
    add_vec(&mut superposed, &d.pressure_pa);
    let err: f64 = c
        .pressure_pa
        .iter()
        .zip(&superposed)
        .map(|(a, b)| (a - b).abs())
        .sum();
    assert!(
        err > 1.0e-5,
        "a shared plate must not equal independently evolved members"
    );
}

#[test]
fn fixed_fixed_string_plate_duct_ode_loads_the_chain() {
    let mut a = plucked(80.0, 0.006, 0.004);
    a.plate = Some(steel_panel());
    a.duct = Some(ViscothermalDuct {
        segments: vec![CylinderSegment {
            radius_m: 0.012,
            length_m: 0.34,
        }],
        tone_holes: vec![],
        termination: WaveguideEnd::Closed,
    });
    a.duration_s = 0.05;
    let with = realize_assembly(&a).expect("ff string-plate-duct");
    a.duct = None;
    let bare = realize_assembly(&a).expect("no duct");
    let err: f64 = with
        .pressure_pa
        .iter()
        .zip(&bare.pressure_pa)
        .map(|(x, y)| (x - y).abs())
        .sum();
    assert!(
        err > 1.0e-6,
        "a one-way bridge into the ODE chain must load the plate"
    );
}

fn add_vec(acc: &mut [f64], add: &[f64]) {
    for (a, b) in acc.iter_mut().zip(add) {
        *a += *b;
    }
}

#[test]
fn von_karman_plate_is_a_body_not_a_guitar() {
    let mut a = plucked(80.0, 0.006, 0.006);
    let mut plate = steel_panel();
    plate.n_modes = 1;
    plate.geometric_nonlinearity = true;
    a.plate = Some(plate);
    let nl = realize_assembly(&a).expect("vk");
    plate.geometric_nonlinearity = false;
    a.plate = Some(plate);
    let lin = realize_assembly(&a).expect("linear plate");
    let err: f64 = nl
        .pressure_pa
        .iter()
        .zip(&lin.pressure_pa)
        .map(|(x, y)| (x - y).abs())
        .sum();
    assert!(err > 1.0e-8, "von Karman must differ from the linear bank");
}

#[test]
fn kirchhoff_carrier_three_way_is_not_the_linear_string() {
    let mut linear = plucked(80.0, 0.006, 0.003);
    let air = closed_duct(288.15);
    linear.duct = air.duct;
    linear.blow = air.blow;
    linear.plate = Some(steel_panel());
    linear.duration_s = 0.06;
    let mut kc = linear.clone();
    if let Some(s) = kc.string.as_mut() {
        s.axial_stiffness_n = 4.0e4;
        s.n_modes = 3;
    }
    if let Some(s) = linear.string.as_mut() {
        s.n_modes = 3;
    }
    let a = realize_assembly(&linear).expect("linear three-way");
    let b = realize_assembly(&kc).expect("KC three-way");
    let err: f64 = a
        .pressure_pa
        .iter()
        .zip(&b.pressure_pa)
        .map(|(x, y)| (x - y).abs())
        .sum();
    assert!(
        err > 1.0e-6,
        "EA > 0 must change the shared-clock string member"
    );
}

#[test]
fn declared_contact_texture_modulates_the_bow() {
    let mut smooth = empty_base();
    smooth.string = Some(nylon_like(80.0, 0.006));
    smooth.bow = Some(BowStroke {
        station_frac: 0.12,
        normal_force_n: 0.8,
        velocity_m_s: 0.4,
        mu_static: 0.8,
        mu_dynamic: 0.3,
        stribeck_m_s: 0.05,
    });
    let mut rough = smooth.clone();
    rough.contact_texture = Some(ContactTexture {
        rms_height_m: 8.0e-6,
        hurst_exponent: 0.8,
        min_cycles: 2,
        max_cycles: 8,
        phase_seed: 7,
        track_length_m: 0.12,
        tangent_stiffness_n_m: 2.0e5,
    });
    let a = realize_assembly(&smooth).expect("smooth bow");
    let b = realize_assembly(&rough).expect("textured bow");
    let err: f64 = a
        .pressure_pa
        .iter()
        .zip(&b.pressure_pa)
        .map(|(x, y)| (x - y).abs())
        .sum();
    assert!(
        err > 1.0e-8,
        "a declared contact texture must change the bowed waveform"
    );
}

#[test]
fn even_mode_compact_dipole_reaches_the_observer() {
    let mut odd = plucked(80.0, 0.006, 0.003);
    if let Some(s) = odd.string.as_mut() {
        s.n_modes = 1;
    }
    let mut both = odd.clone();
    if let Some(s) = both.string.as_mut() {
        s.n_modes = 2;
    }
    let a = realize_assembly(&odd).expect("odd monopole");
    let b = realize_assembly(&both).expect("odd+even");
    let err: f64 = a
        .pressure_pa
        .iter()
        .zip(&b.pressure_pa)
        .map(|(x, y)| (x - y).abs())
        .sum();
    assert!(
        err > 1.0e-8,
        "even sine modes must radiate as compact dipoles ({err})"
    );
}

#[test]
fn kc_contact_is_inside_the_hamiltonian() {
    let mut bare = plucked(70.0, 0.005, 0.004);
    if let Some(s) = bare.string.as_mut() {
        s.axial_stiffness_n = 4.0e4;
        s.n_modes = 3;
    }
    let mut rattle = bare.clone();
    rattle.obstacles = vec![UnilateralObstacle {
        stations: vec![0.25, 0.5, 0.75],
        gaps_m: vec![0.001, 0.001, 0.001],
        stiffness: 5.0e5,
        alpha: 2.0,
        mu_kinetic: 0.0,
        internal_loss: 0.0,
        provenance: "fixture/stay".into(),
    }];
    let a = realize_assembly(&bare).expect("bare KC");
    let b = realize_assembly(&rattle).expect("KC contact");
    let err: f64 = a
        .pressure_pa
        .iter()
        .zip(&b.pressure_pa)
        .map(|(x, y)| (x - y).abs())
        .sum();
    assert!(
        err > 1.0e-6,
        "ContactStorage wrapping KC must change the waveform ({err})"
    );
}

#[test]
fn polarizations_share_the_plate() {
    let mut a = plucked(80.0, 0.006, 0.003);
    a.duration_s = 0.16;
    if let Some(s) = a.string.as_mut() {
        s.polarization_detune = 0.012;
        s.n_modes = 2;
        s.damping_ratio = 0.001;
    }
    a.soundboard = Some(RadiatingPlate {
        area_m2: 0.12,
        mass_kg: 0.18,
        frequency_hz: 110.0,
        damping_ratio: 0.03,
    });
    let out = realize_assembly(&a).expect("shared-plate polarizations");
    let hop = 80usize;
    let mut env = Vec::new();
    for chunk in out.pressure_pa.chunks(hop) {
        let e = (chunk.iter().map(|p| p * p).sum::<f64>() / chunk.len() as f64).sqrt();
        env.push(e);
    }
    let max_e = env.iter().copied().fold(0.0_f64, f64::max);
    let min_e = env.iter().copied().fold(f64::INFINITY, f64::min);
    assert!(
        max_e > 1.2 * min_e,
        "two polarizations on one plate must still beat ({max_e} vs {min_e})"
    );
}

#[test]
fn helmholtz_cavity_changes_the_plate_waveform() {
    let mut bare = plucked(80.0, 0.006, 0.004);
    bare.plate = Some(steel_panel());
    let mut boxed = bare.clone();
    boxed.cavity = Some(HelmholtzCavity {
        volume_m3: 0.02,
        neck_radius_m: 0.02,
        neck_length_m: 0.03,
    });
    let a = realize_assembly(&bare).expect("free plate");
    let b = realize_assembly(&boxed).expect("cavity");
    let err: f64 = a
        .pressure_pa
        .iter()
        .zip(&b.pressure_pa)
        .map(|(x, y)| (x - y).abs())
        .sum();
    assert!(err > 1.0e-8, "a Helmholtz volume must load the plate");
}

#[test]
fn humid_observer_path_differs_from_dry() {
    let mut dry = plucked(80.0, 0.006, 0.003);
    dry.listener.distance_m = 80.0;
    dry.duration_s = 0.08;
    let mut wet = dry.clone();
    wet.ambient.relative_humidity = 0.70;
    let a = realize_assembly(&dry).expect("dry");
    let b = realize_assembly(&wet).expect("wet");
    let err: f64 = a
        .pressure_pa
        .iter()
        .zip(&b.pressure_pa)
        .map(|(x, y)| (x - y).abs())
        .sum();
    assert!(err > 0.0, "ISO humidity must move a long observer path");
}

#[test]
fn hunt_crossley_changes_rattle() {
    let mut elastic = plucked(80.0, 0.006, 0.008);
    elastic.obstacles = vec![UnilateralObstacle {
        stations: vec![0.25, 0.5, 0.75],
        gaps_m: vec![0.001, 0.001, 0.001],
        stiffness: 5.0e5,
        alpha: 2.0,
        mu_kinetic: 0.0,
        internal_loss: 0.0,
        provenance: "fixture/stay".into(),
    }];
    let mut lossy = elastic.clone();
    lossy.obstacles[0].internal_loss = 0.8;
    let a = realize_assembly(&elastic).expect("elastic");
    let b = realize_assembly(&lossy).expect("HC");
    let err: f64 = a
        .pressure_pa
        .iter()
        .zip(&b.pressure_pa)
        .map(|(x, y)| (x - y).abs())
        .sum();
    assert!(err > 1.0e-8, "Hunt–Crossley must drain a rattle");
}

#[test]
fn moving_end_dirac_join_is_not_the_one_way_bridge() {
    let mut one_way = plucked(80.0, 0.006, 0.004);
    one_way.plate = Some(steel_panel());
    one_way.duration_s = 0.05;
    let mut two_way = one_way.clone();
    if let Some(s) = two_way.string.as_mut() {
        s.moving_end = true;
        s.n_modes = 3;
    }
    let a = realize_assembly(&one_way).expect("fixed-fixed");
    let b = realize_assembly(&two_way).expect("moving-end Dirac");
    let err: f64 = a
        .pressure_pa
        .iter()
        .zip(&b.pressure_pa)
        .map(|(x, y)| (x - y).abs())
        .sum();
    assert!(
        err > 1.0e-6,
        "a free attachment must Dirac-join the plate, not reprint the one-way bridge"
    );
}

#[test]
fn moving_end_with_cavity_is_a_three_phs_clock() {
    let mut a = plucked(80.0, 0.006, 0.004);
    if let Some(s) = a.string.as_mut() {
        s.moving_end = true;
        s.n_modes = 2;
    }
    a.plate = Some(steel_panel());
    a.cavity = Some(HelmholtzCavity {
        volume_m3: 0.016,
        neck_radius_m: 0.02,
        neck_length_m: 0.03,
    });
    a.duration_s = 0.04;
    let with = realize_assembly(&a).expect("three-pHS");
    a.cavity = None;
    let bare = realize_assembly(&a).expect("no cavity");
    let err: f64 = with
        .pressure_pa
        .iter()
        .zip(&bare.pressure_pa)
        .map(|(x, y)| (x - y).abs())
        .sum();
    assert!(err > 1.0e-8, "the cavity transformer must load the join");
}

#[test]
fn moving_end_kc_dirac_is_not_the_linear_moving_end() {
    let mut linear = plucked(70.0, 0.005, 3.5e-3);
    if let Some(s) = linear.string.as_mut() {
        s.moving_end = true;
        s.n_modes = 3;
    }
    linear.plate = Some(steel_panel());
    linear.duration_s = 0.05;
    let mut kc = linear.clone();
    if let Some(s) = kc.string.as_mut() {
        s.axial_stiffness_n = 4.0e4;
    }
    let a = realize_assembly(&linear).expect("linear moving-end");
    let b = realize_assembly(&kc).expect("KC moving-end Dirac");
    let err: f64 = a
        .pressure_pa
        .iter()
        .zip(&b.pressure_pa)
        .map(|(x, y)| (x - y).abs())
        .sum();
    assert!(
        err > 1.0e-6,
        "free-fixed KC must change the Dirac join, not reprint the linear waveguide"
    );
}

#[test]
fn bow_on_moving_end_dirac_is_a_live_port() {
    let mut a = empty_base();
    a.duration_s = 0.06;
    let mut s = nylon_like(80.0, 0.006);
    s.moving_end = true;
    s.n_modes = 4;
    a.string = Some(s);
    a.plate = Some(steel_panel());
    a.bow = Some(BowStroke {
        station_frac: 0.13,
        normal_force_n: 1.0,
        velocity_m_s: 0.35,
        mu_static: 0.8,
        mu_dynamic: 0.3,
        stribeck_m_s: 0.04,
    });
    let bowed = realize_assembly(&a).expect("bowed Dirac");
    assert!(
        peak_abs(&bowed.pressure_pa) > 1.0e-4,
        "a bow on the leftover Dirac port must drive the string"
    );
    a.bow = None;
    a.pluck = Some(Pluck {
        station_frac: 0.25,
        height_m: 0.003,
    });
    let plucked = realize_assembly(&a).expect("plucked Dirac");
    let err: f64 = bowed
        .pressure_pa
        .iter()
        .zip(&plucked.pressure_pa)
        .map(|(x, y)| (x - y).abs())
        .sum();
    assert!(err > 1.0e-6, "bow drive must not reprint a pluck");
}

#[test]
fn obstacle_on_moving_end_dirac_changes_the_waveform() {
    let mut a = plucked(80.0, 0.006, 0.008);
    if let Some(s) = a.string.as_mut() {
        s.moving_end = true;
        s.n_modes = 3;
    }
    a.plate = Some(steel_panel());
    a.duration_s = 0.05;
    let bare = realize_assembly(&a).expect("bare Dirac");
    a.obstacles = vec![UnilateralObstacle {
        stations: vec![0.3],
        gaps_m: vec![0.001],
        stiffness: 5.0e5,
        alpha: 2.0,
        mu_kinetic: 0.0,
        internal_loss: 0.0,
        provenance: "fixture/stay".into(),
    }];
    let rattle = realize_assembly(&a).expect("obstacle Dirac");
    let err: f64 = bare
        .pressure_pa
        .iter()
        .zip(&rattle.pressure_pa)
        .map(|(x, y)| (x - y).abs())
        .sum();
    assert!(
        err > 1.0e-6,
        "an obstacle leftover port must change the join"
    );
    a.obstacles[0].internal_loss = 0.8;
    let lossy = realize_assembly(&a).expect("HC Dirac");
    let err_hc: f64 = rattle
        .pressure_pa
        .iter()
        .zip(&lossy.pressure_pa)
        .map(|(x, y)| (x - y).abs())
        .sum();
    assert!(err_hc > 1.0e-8, "Hunt–Crossley on Dirac must drain");
}

#[test]
fn von_karman_duct_dirac_is_not_the_linear_plate_duct() {
    let mut lin = plucked(80.0, 0.006, 0.005);
    if let Some(s) = lin.string.as_mut() {
        s.moving_end = true;
        s.n_modes = 2;
    }
    let mut plate = steel_panel();
    plate.n_modes = 1;
    lin.plate = Some(plate);
    lin.duct = Some(ViscothermalDuct {
        segments: vec![CylinderSegment {
            radius_m: 0.012,
            length_m: 0.34,
        }],
        tone_holes: vec![],
        termination: WaveguideEnd::Closed,
    });
    lin.duration_s = 0.04;
    let mut vk = lin.clone();
    if let Some(p) = vk.plate.as_mut() {
        p.geometric_nonlinearity = true;
    }
    let a = realize_assembly(&lin).expect("linear plate×duct");
    let b = realize_assembly(&vk).expect("VK×duct");
    let err: f64 = a
        .pressure_pa
        .iter()
        .zip(&b.pressure_pa)
        .map(|(x, y)| (x - y).abs())
        .sum();
    assert!(
        err > 1.0e-8,
        "von Karman×duct must differ from the linear join"
    );
}

#[test]
fn blow_on_moving_end_string_plate_duct_loads_the_join() {
    let mut a = plucked(80.0, 0.006, 0.004);
    if let Some(s) = a.string.as_mut() {
        s.moving_end = true;
        s.n_modes = 2;
    }
    a.plate = Some(steel_panel());
    a.duct = Some(ViscothermalDuct {
        segments: vec![CylinderSegment {
            radius_m: 0.012,
            length_m: 0.34,
        }],
        tone_holes: vec![],
        termination: WaveguideEnd::Closed,
    });
    a.duration_s = 0.04;
    let quiet = realize_assembly(&a).expect("no blow");
    a.blow = Some(VolumeVelocityPulse {
        peak_m3_s: 2.0e-5,
        duration_s: 0.002,
    });
    let blown = realize_assembly(&a).expect("blow leftover");
    let err: f64 = quiet
        .pressure_pa
        .iter()
        .zip(&blown.pressure_pa)
        .map(|(x, y)| (x - y).abs())
        .sum();
    assert!(
        err > 1.0e-6,
        "a leftover blow port must load the Dirac join"
    );
}

#[test]
fn reed_on_moving_end_string_plate_duct_speaks() {
    let mut a = plucked(80.0, 0.006, 0.004);
    if let Some(s) = a.string.as_mut() {
        s.moving_end = true;
        s.n_modes = 2;
    }
    a.plate = Some(steel_panel());
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
        blowing_pressure_pa: 2_800.0,
        attack_s: 0.008,
        mass_kg: 0.0,
        stiffness_n_m: 0.0,
    });
    a.duration_s = 0.08;
    let out = realize_assembly(&a).expect("reed leftover");
    let mut silent = a.clone();
    if let Some(reed) = silent.reed.as_mut() {
        reed.blowing_pressure_pa = 0.0;
    }
    let quiet = realize_assembly(&silent).expect("silent leftover reed");
    let err: f64 = out
        .pressure_pa
        .iter()
        .zip(&quiet.pressure_pa)
        .map(|(x, y)| (x - y).abs())
        .sum();
    assert!(
        err > 1.0e-6,
        "a leftover Bernoulli reed must change the Dirac join"
    );
}

#[test]
fn moving_end_with_duct_is_a_three_phs_clock() {
    let mut a = plucked(80.0, 0.006, 0.004);
    if let Some(s) = a.string.as_mut() {
        s.moving_end = true;
        s.n_modes = 2;
    }
    a.plate = Some(steel_panel());
    a.duct = Some(ViscothermalDuct {
        segments: vec![CylinderSegment {
            radius_m: 0.012,
            length_m: 0.34,
        }],
        tone_holes: vec![],
        termination: WaveguideEnd::Closed,
    });
    a.duration_s = 0.04;
    let with = realize_assembly(&a).expect("string-plate-duct");
    a.duct = None;
    let bare = realize_assembly(&a).expect("no duct");
    let err: f64 = with
        .pressure_pa
        .iter()
        .zip(&bare.pressure_pa)
        .map(|(x, y)| (x - y).abs())
        .sum();
    assert!(err > 1.0e-8, "the duct transformer must load the join");
}

#[test]
fn simple_cylinder_ode_still_tracks_sound_speed() {
    let cold = realize_assembly(&closed_duct(288.15)).expect("cold ODE");
    let hot = realize_assembly(&closed_duct(330.0)).expect("hot ODE");
    let p_cold = dominant_period_samples(&cold.pressure_pa, 8, 40);
    let p_hot = dominant_period_samples(&hot.pressure_pa, 8, 40);
    let measured = p_hot as f64 / p_cold as f64;
    let expected = cold.gas.sound_speed / hot.gas.sound_speed;
    assert!(
        (measured - expected).abs() < 0.12,
        "ODE cylinder period must track 1/c(T): {measured} vs {expected}"
    );
}
