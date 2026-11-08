//! Continuous-band transport checked against independent dimensional Benney
//! coefficients, analytic extrema, and the existing single-frequency owner.

use fs_lbm::film::{FilmError, FilmMaterial, LongWaveLimits};
use fs_vessel::{
    LaminarFilmError, LaminarFilmSection, LaminarFilmTransitBandStudy, LaminarFilmTransitStudy,
    LaminarPourConditions, LaminarTransitSection,
};

const MATERIAL: FilmMaterial = FilmMaterial {
    density_kg_m3: 1000.0,
    kinematic_viscosity_m2_s: 1e-5,
    surface_tension_n_m: 0.072,
};
const FIRST: LaminarFilmSection = LaminarFilmSection {
    width_m: 0.04,
    tangential_acceleration_m_s2: 4.0,
    normal_acceleration_m_s2: 3.0,
};
const SECOND: LaminarFilmSection = LaminarFilmSection {
    width_m: 0.08,
    ..FIRST
};
const LIMITS: LongWaveLimits = LongWaveLimits {
    max_kh: 0.1,
    max_re_kh: 0.1,
};

fn conditions() -> LaminarPourConditions {
    LaminarPourConditions {
        volume_flow_m3_s: 4e-7,
        material: MATERIAL,
    }
}

fn study(sections: &[LaminarTransitSection]) -> LaminarFilmTransitBandStudy<'_> {
    LaminarFilmTransitBandStudy {
        sections,
        angular_frequencies_rad_s: [1.0, 12.0],
        limits: LIMITS,
        max_growth_frequency_ratio: 0.1,
    }
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5e-12 * expected.abs().max(1e-20),
        "{actual:e} != {expected:e}"
    );
}

struct Reference {
    thickness: f64,
    speed: f64,
    time: f64,
    growth_quadratic: f64,
    growth_quartic: f64,
    gain_quadratic: f64,
    gain_quartic: f64,
}

// Independent dimensional coefficients from the film thickness. The native
// dispersion, coefficient accessor and band maximizer are not used here.
fn reference(segment: LaminarTransitSection, state: LaminarPourConditions) -> Reference {
    let q = state.volume_flow_m3_s / segment.film.width_m;
    let nu = state.material.kinematic_viscosity_m2_s;
    let gt = segment.film.tangential_acceleration_m_s2;
    let gn = segment.film.normal_acceleration_m_s2;
    let h = (3.0 * nu * q / gt).cbrt();
    let speed = gt * h * h / nu;
    let a = 2.0 * gt * gt * h.powi(6) / (15.0 * nu.powi(3)) - gn * h.powi(3) / (3.0 * nu);
    let b =
        state.material.surface_tension_n_m * h.powi(3) / (3.0 * state.material.density_kg_m3 * nu);
    let time = segment.length_m / speed;
    Reference {
        thickness: h,
        speed,
        time,
        growth_quadratic: a,
        growth_quartic: b,
        gain_quadratic: time * a / speed.powi(2),
        gain_quartic: time * b / speed.powi(4),
    }
}

#[test]
fn g1_band_transport_has_distinct_interior_outlet_and_upstream_optima() {
    let path = [
        LaminarTransitSection {
            film: FIRST,
            length_m: 0.1,
        },
        LaminarTransitSection {
            film: SECOND,
            length_m: 0.02,
        },
    ];
    let first = reference(path[0], conditions());
    let second = reference(path[1], conditions());
    let a = first.gain_quadratic + second.gain_quadratic;
    let b = first.gain_quartic + second.gain_quartic;
    let outlet_frequency = (a / (2.0 * b)).sqrt();
    let upstream_frequency = (first.gain_quadratic / (2.0 * first.gain_quartic)).sqrt();
    assert!(outlet_frequency > 1.0 && outlet_frequency < upstream_frequency);
    assert!(upstream_frequency < 12.0);
    let report = study(&path).evaluate(1.0, conditions()).unwrap();
    near(report.outlet_angular_frequency_rad_s, outlet_frequency);
    near(report.max_outlet_log_amplification, a * a / (4.0 * b));
    near(
        report.peak_angular_frequency_rad_s.unwrap(),
        upstream_frequency,
    );
    near(
        report.peak_log_amplification,
        first.gain_quadratic.powi(2) / (4.0 * first.gain_quartic),
    );
    assert_eq!(report.peak_after_sections, 1);
    assert!(report.peak_log_amplification > report.max_outlet_log_amplification);
    for omega in [1.0_f64, 12.0] {
        assert!(report.max_outlet_log_amplification > a * omega.powi(2) - b * omega.powi(4));
    }
    near(report.wave_transit_time_s, first.time + second.time);
    near(
        report.liquid_residence_time_s,
        3.0 * (first.time + second.time),
    );
    near(report.wavenumber_range_rad_m[0], 1.0 / first.speed);
    near(report.wavenumber_range_rad_m[1], 12.0 / second.speed);
    near(report.max_kh, 12.0 * second.thickness / second.speed);
    near(report.reynolds_range[0], 0.5);
    near(report.reynolds_range[1], 1.0);
}

#[test]
fn g2_weak_growth_admission_includes_the_hidden_cubic_maximum_and_damping_tail() {
    let path = [LaminarTransitSection {
        film: FIRST,
        length_m: 0.1,
    }];
    let r = reference(path[0], conditions());
    let cutoff = r.speed * (r.growth_quadratic / r.growth_quartic).sqrt();
    let cubic_frequency = cutoff / 3.0_f64.sqrt();
    let alpha = r.growth_quadratic / r.speed.powi(2);
    let beta = r.growth_quartic / r.speed.powi(4);
    let expected = alpha * cubic_frequency - beta * cubic_frequency.powi(3);
    let band = LaminarFilmTransitBandStudy {
        angular_frequencies_rad_s: [0.0, cutoff],
        ..study(&path)
    };
    let admitted = band.evaluate(1.0, conditions()).unwrap();
    near(admitted.max_growth_frequency_ratio, expected);
    let strict = LaminarFilmTransitBandStudy {
        max_growth_frequency_ratio: 0.5 * expected,
        ..band
    };
    // Both endpoints are neutral, yet the interior convective growth violates
    // the approximation budget. Checking only the endpoints would accept it.
    for omega in [0.0, cutoff] {
        assert!(
            LaminarFilmTransitBandStudy {
                angular_frequencies_rad_s: [omega, omega],
                ..strict
            }
            .evaluate(1.0, conditions())
            .is_ok()
        );
    }
    assert_eq!(
        strict.evaluate(1.0, conditions()),
        Err(LaminarFilmError::OutsideSpatialGrowthRegime { section: 0 })
    );
    assert!(
        LaminarFilmTransitBandStudy {
            max_growth_frequency_ratio: admitted.max_growth_frequency_ratio,
            ..band
        }
        .evaluate(1.0, conditions())
        .is_ok()
    );

    // Strong damping also invalidates the temporal-to-spatial conversion.
    let upper = 1.2 * cutoff;
    let tail = (alpha * upper - beta * upper.powi(3)).abs();
    assert!(tail > expected);
    let extended = LaminarFilmTransitBandStudy {
        angular_frequencies_rad_s: [0.0, upper],
        ..band
    };
    near(
        extended
            .evaluate(1.0, conditions())
            .unwrap()
            .max_growth_frequency_ratio,
        tail,
    );
    assert_eq!(
        LaminarFilmTransitBandStudy {
            max_growth_frequency_ratio: 0.5 * (tail + expected),
            ..extended
        }
        .evaluate(1.0, conditions()),
        Err(LaminarFilmError::OutsideSpatialGrowthRegime { section: 0 })
    );
}

#[test]
fn g2_band_transport_refuses_the_complete_unsupported_tail_and_invalid_inputs() {
    let path = [LaminarTransitSection {
        film: FIRST,
        length_m: 0.1,
    }];
    let limited = LaminarFilmTransitBandStudy {
        limits: LongWaveLimits {
            max_kh: 0.05,
            ..LIMITS
        },
        ..study(&path)
    };
    let r = reference(path[0], conditions());
    let optimum = (r.gain_quadratic / (2.0 * r.gain_quartic)).sqrt();
    assert!(
        LaminarFilmTransitBandStudy {
            angular_frequencies_rad_s: [optimum, optimum],
            ..limited
        }
        .evaluate(1.0, conditions())
        .is_ok()
    );
    assert!(matches!(
        limited.evaluate(1.0, conditions()),
        Err(LaminarFilmError::Band {
            section: 0,
            source: FilmError::OutsideLongWaveRegime
        })
    ));
    for band in [
        [-1.0, 1.0],
        [2.0, 1.0],
        [0.0, f64::INFINITY],
        [f64::NAN, 1.0],
    ] {
        assert_eq!(
            LaminarFilmTransitBandStudy {
                angular_frequencies_rad_s: band,
                ..study(&path)
            }
            .evaluate(1.0, conditions()),
            Err(LaminarFilmError::InvalidTransitInput)
        );
    }
    for max_growth_frequency_ratio in [0.0, 1.0, f64::NAN] {
        assert_eq!(
            LaminarFilmTransitBandStudy {
                max_growth_frequency_ratio,
                ..study(&path)
            }
            .evaluate(1.0, conditions()),
            Err(LaminarFilmError::InvalidTransitInput)
        );
    }
    assert_eq!(
        study(&[]).evaluate(1.0, conditions()),
        Err(LaminarFilmError::EmptyPath)
    );
    let invalid = [LaminarTransitSection {
        length_m: 0.0,
        ..path[0]
    }];
    assert_eq!(
        study(&invalid).evaluate(1.0, conditions()),
        Err(LaminarFilmError::InvalidLength { section: 0 })
    );
}

#[test]
fn g1_decaying_and_neutral_bands_leave_the_inlet_peak_without_a_frequency() {
    let path = [LaminarTransitSection {
        film: SECOND,
        length_m: 0.5,
    }];
    let band = LaminarFilmTransitBandStudy {
        angular_frequencies_rad_s: [1.0, 8.0],
        ..study(&path)
    };
    let r = reference(path[0], conditions());
    let report = band.evaluate(1.0, conditions()).unwrap();
    assert!(report.max_outlet_log_amplification < 0.0);
    near(
        report.max_outlet_log_amplification,
        r.gain_quadratic - r.gain_quartic,
    );
    assert_eq!(report.outlet_angular_frequency_rad_s, 1.0);
    assert_eq!(report.peak_log_amplification, 0.0);
    assert_eq!(report.peak_after_sections, 0);
    assert_eq!(report.peak_angular_frequency_rad_s, None);
    for frequencies in [[0.0, 8.0], [0.0, 0.0]] {
        let report = LaminarFilmTransitBandStudy {
            angular_frequencies_rad_s: frequencies,
            ..band
        }
        .evaluate(1.0, conditions())
        .unwrap();
        assert_eq!(report.max_outlet_log_amplification, 0.0);
        assert_eq!(report.outlet_angular_frequency_rad_s, 0.0);
        assert_eq!(report.peak_log_amplification, 0.0);
        assert_eq!(report.peak_after_sections, 0);
        assert_eq!(report.peak_angular_frequency_rad_s, None);
        near(report.wave_transit_time_s, r.time);
        if frequencies[1] == 0.0 {
            assert_eq!(report.max_growth_frequency_ratio, 0.0);
            assert_eq!(report.wavenumber_range_rad_m, [0.0, 0.0]);
        }
    }
}

#[test]
fn g1_zero_capillarity_retains_the_increasing_band_endpoint() {
    let path = [LaminarTransitSection {
        film: FIRST,
        length_m: 0.1,
    }];
    let state = LaminarPourConditions {
        material: FilmMaterial {
            surface_tension_n_m: 0.0,
            ..MATERIAL
        },
        ..conditions()
    };
    let r = reference(path[0], state);
    let report = study(&path).evaluate(1.0, state).unwrap();
    assert_eq!(report.outlet_angular_frequency_rad_s, 12.0);
    assert_eq!(report.peak_angular_frequency_rad_s, Some(12.0));
    assert_eq!(report.peak_after_sections, 1);
    near(
        report.max_outlet_log_amplification,
        r.gain_quadratic * 12.0_f64.powi(2),
    );
    near(
        report.peak_log_amplification,
        report.max_outlet_log_amplification,
    );
}

#[test]
fn g3_point_bands_and_section_partition_recover_the_original_transit_model() {
    let path = [
        LaminarTransitSection {
            film: FIRST,
            length_m: 0.1,
        },
        LaminarTransitSection {
            film: SECOND,
            length_m: 0.5,
        },
    ];
    let band = LaminarFilmTransitBandStudy {
        angular_frequencies_rad_s: [0.5, 0.5],
        ..study(&path)
    };
    let one = LaminarFilmTransitStudy {
        sections: &path,
        angular_frequency_rad_s: 0.5,
        limits: band.limits,
        max_growth_frequency_ratio: band.max_growth_frequency_ratio,
    }
    .evaluate(1.0, conditions())
    .unwrap();
    let whole_band = band.evaluate(1.0, conditions()).unwrap();
    near(
        whole_band.max_outlet_log_amplification,
        one.log_amplification,
    );
    near(
        whole_band.peak_log_amplification,
        one.peak_log_amplification,
    );
    assert_eq!(whole_band.peak_after_sections, one.peak_after_sections);
    assert_eq!(whole_band.outlet_angular_frequency_rad_s, 0.5);
    assert_eq!(whole_band.peak_angular_frequency_rad_s, Some(0.5));
    near(whole_band.wave_transit_time_s, one.wave_transit_time_s);
    near(
        whole_band.max_growth_frequency_ratio,
        one.max_growth_frequency_ratio,
    );

    let whole = [LaminarTransitSection {
        film: FIRST,
        length_m: 0.3,
    }];
    let split = [LaminarTransitSection {
        film: FIRST,
        length_m: 0.1,
    }; 3];
    let a = study(&whole).evaluate(1.0, conditions()).unwrap();
    let b = study(&split).evaluate(1.0, conditions()).unwrap();
    near(
        a.max_outlet_log_amplification,
        b.max_outlet_log_amplification,
    );
    near(a.peak_log_amplification, b.peak_log_amplification);
    near(
        a.outlet_angular_frequency_rad_s,
        b.outlet_angular_frequency_rad_s,
    );
    near(
        a.peak_angular_frequency_rad_s.unwrap(),
        b.peak_angular_frequency_rad_s.unwrap(),
    );
    near(a.wave_transit_time_s, b.wave_transit_time_s);
    near(a.liquid_residence_time_s, b.liquid_residence_time_s);
    assert_eq!(a.peak_after_sections, 1);
    assert_eq!(b.peak_after_sections, 3);
}
