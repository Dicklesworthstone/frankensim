//! Physical film integration: independent dimensional references and refusals.

use fs_lbm::film::{FilmError, FilmMaterial, LongWaveLimits};
use fs_vessel::race::{ScreenError, race_base_losses, screen_laminar_lips};
use fs_vessel::robust::{LaminarRobustError, robustify_laminar};
use fs_vessel::stability::laminar::{
    LaminarFilmBandStudy, LaminarFilmError, LaminarFilmSection, LaminarFilmStudy,
    LaminarPourConditions,
};

const MATERIAL: FilmMaterial = FilmMaterial {
    density_kg_m3: 1000.0,
    kinematic_viscosity_m2_s: 1e-5,
    surface_tension_n_m: 0.072,
};
const SECTIONS: [LaminarFilmSection; 2] = [
    LaminarFilmSection {
        width_m: 0.04,
        tangential_acceleration_m_s2: 4.0,
        normal_acceleration_m_s2: 3.0,
    },
    LaminarFilmSection {
        width_m: 0.08,
        tangential_acceleration_m_s2: 4.0,
        normal_acceleration_m_s2: 3.0,
    },
];
const WAVES: [f64; 3] = [2.0, 4.0, 6.0];
const WIDTHS: [f64; 5] = [0.7, 1.0, 1.5, 2.0, 3.0];

#[test]
fn continuous_band_reaches_existing_cvar_and_race_owners() {
    let study = LaminarFilmBandStudy {
        sections: &SECTIONS,
        wavenumbers_rad_m: [2.0, 200.0],
        limits: LongWaveLimits {
            max_kh: 0.3,
            max_re_kh: 0.5,
        },
    };
    let nominal = conditions(1.0, 1.0);
    let band = [
        conditions(0.8, 1.6),
        conditions(1.3, 0.7),
        conditions(1.1, 1.2),
    ];
    let risk = robustify_laminar(&study, &WIDTHS, nominal, &band, 0.5).unwrap();
    assert!(risk.robust_band_cvar_per_s <= risk.nominal_band_cvar_per_s);
    let mut nominal_losses = Vec::new();
    let mut risks = Vec::new();
    for &width in &WIDTHS {
        let result = study.evaluate(width, nominal).unwrap();
        nominal_losses.push(result.max_growth_per_s);
        let independent_max = SECTIONS
            .iter()
            .map(|section| {
                let q = nominal.volume_flow_m3_s / (width * section.width_m);
                let a = q * (1.2 * q / MATERIAL.kinematic_viscosity_m2_s - 0.75);
                let b = q * MATERIAL.surface_tension_n_m / (MATERIAL.density_kg_m3 * 4.0);
                let k = if a > 0.0 {
                    (a / (2.0 * b)).sqrt().clamp(2.0, 200.0)
                } else {
                    2.0
                };
                a * k * k - b * k.powi(4)
            })
            .fold(f64::NEG_INFINITY, f64::max);
        assert!((result.max_growth_per_s - independent_max).abs() < 2e-13);
        let losses: Vec<_> = band
            .iter()
            .map(|&fluid| study.evaluate(width, fluid).unwrap().max_growth_per_s)
            .collect();
        risks.push(fs_vessel::robust::cvar(&losses, 0.5).unwrap());
    }
    assert_eq!(
        risk.robust_candidate,
        risks
            .iter()
            .enumerate()
            .min_by(|a, b| a.1.total_cmp(b.1))
            .unwrap()
            .0
    );
    assert_eq!(
        screen_laminar_lips(&study, &WIDTHS, nominal, 1234).unwrap(),
        race_base_losses(&nominal_losses, 1234).unwrap()
    );
    let refused = LaminarFilmBandStudy {
        wavenumbers_rad_m: [0.0, 1e9],
        ..study
    };
    assert!(matches!(
        robustify_laminar(&refused, &WIDTHS, nominal, &band, 0.5),
        Err(LaminarRobustError::Film {
            candidate: 0,
            fluid: None,
            source: LaminarFilmError::Band {
                section: 0,
                source: FilmError::OutsideLongWaveRegime
            }
        })
    ));
}

fn study() -> LaminarFilmStudy<'static> {
    LaminarFilmStudy {
        sections: &SECTIONS,
        wavenumbers_rad_m: &WAVES,
        limits: LongWaveLimits {
            max_kh: 0.1,
            max_re_kh: 0.1,
        },
    }
}

fn conditions(rate: f64, viscosity: f64) -> LaminarPourConditions {
    LaminarPourConditions {
        volume_flow_m3_s: 4e-7 * rate,
        material: FilmMaterial {
            kinematic_viscosity_m2_s: MATERIAL.kinematic_viscosity_m2_s * viscosity,
            ..MATERIAL
        },
    }
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1e-20),
        "{actual:e} != {expected:e}"
    );
}

/// Independent dimensional coefficients in powers of h, not the native
/// factorized Reynolds formula or any call to its film/mode evaluator.
fn reference(width_scale: f64, state: LaminarPourConditions) -> f64 {
    let nu = state.material.kinematic_viscosity_m2_s;
    let rho = state.material.density_kg_m3;
    let sigma = state.material.surface_tension_n_m;
    let mut worst = f64::NEG_INFINITY;
    for section in SECTIONS {
        let q = state.volume_flow_m3_s / (section.width_m * width_scale);
        let gt = section.tangential_acceleration_m_s2;
        let gn = section.normal_acceleration_m_s2;
        let h = (3.0 * nu * q / gt).cbrt();
        let inertial = 2.0 * gt * gt * h.powi(6) / (15.0 * nu.powi(3));
        let hydrostatic = gn * h.powi(3) / (3.0 * nu);
        let capillary = sigma * h.powi(3) / (3.0 * rho * nu);
        for k in WAVES {
            worst = worst.max((inertial - hydrostatic) * k * k - capillary * k.powi(4));
        }
    }
    worst
}

#[test]
fn laminar_path_uses_conserved_si_flux_and_free_surface_modes() {
    let study = study();
    for width in WIDTHS {
        for state in [conditions(1.0, 1.0), conditions(1.3, 0.7)] {
            let result = study.evaluate(width, state).unwrap();
            near(result.max_growth_per_s, reference(width, state));
            near(
                result.reynolds_range[1],
                state.volume_flow_m3_s
                    / (SECTIONS[0].width_m * width * state.material.kinematic_viscosity_m2_s),
            );
            near(result.reynolds_range[0], 0.5 * result.reynolds_range[1]);
            assert!(result.max_kh < 0.1 && result.max_re_kh < 0.1);
        }
    }
    let nominal = study.evaluate(1.0, conditions(1.0, 1.0)).unwrap();
    assert_eq!((nominal.section, nominal.mode), (0, 2));
    assert!(nominal.max_growth_per_s > 0.0);
    let wider = study.evaluate(2.0, conditions(1.0, 1.0)).unwrap();
    let viscous = study.evaluate(1.0, conditions(1.0, 2.0)).unwrap();
    near(wider.reynolds_range[1], nominal.reynolds_range[1] / 2.0);
    near(viscous.reynolds_range[1], nominal.reynolds_range[1] / 2.0);
    assert!(wider.max_growth_per_s < 0.0 && viscous.max_growth_per_s < 0.0);
}

#[test]
fn laminar_domain_failure_never_becomes_a_stable_score_or_skipped_station() {
    let study = study();
    let nominal = conditions(1.0, 1.0);
    for width in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            study.evaluate(width, nominal),
            Err(LaminarFilmError::InvalidWidthScale)
        );
    }
    assert_eq!(
        study.evaluate(1.0, conditions(0.0, 1.0)),
        Err(LaminarFilmError::InvalidVolumeFlow)
    );
    assert_eq!(
        LaminarFilmStudy {
            sections: &[],
            ..study
        }
        .evaluate(1.0, nominal),
        Err(LaminarFilmError::EmptyPath)
    );
    assert_eq!(
        LaminarFilmStudy {
            wavenumbers_rad_m: &[],
            ..study
        }
        .evaluate(1.0, nominal),
        Err(LaminarFilmError::EmptyModes)
    );
    let mut sections = SECTIONS;
    sections[1].tangential_acceleration_m_s2 = 0.0;
    assert!(matches!(
        LaminarFilmStudy {
            sections: &sections,
            ..study
        }
        .evaluate(1.0, nominal),
        Err(LaminarFilmError::Film { section: 1, .. })
    ));
    assert!(matches!(
        LaminarFilmStudy {
            wavenumbers_rad_m: &[2.0, 1e6],
            ..study
        }
        .evaluate(1.0, nominal),
        Err(LaminarFilmError::Mode {
            section: 0,
            mode: 1,
            source: FilmError::OutsideLongWaveRegime,
        })
    ));
}

#[test]
fn laminar_cvar_selects_physical_candidates_and_allows_a_shared_boundary_optimum() {
    let study = study();
    let nominal = conditions(1.0, 1.0);
    let band = [
        conditions(0.8, 1.6),
        conditions(1.3, 0.7),
        conditions(1.1, 1.2),
    ];
    let report = robustify_laminar(&study, &WIDTHS, nominal, &band, 0.5).unwrap();
    let mut best_nominal = (usize::MAX, f64::INFINITY);
    let mut best_robust = (usize::MAX, f64::INFINITY);
    for (candidate, &width) in WIDTHS.iter().enumerate() {
        let value = reference(width, nominal);
        if value < best_nominal.1 {
            best_nominal = (candidate, value);
        }
        let mut losses = band.map(|state| reference(width, state));
        losses.sort_by(f64::total_cmp);
        // Three equally weighted samples, alpha=1/2: one full largest
        // sample and half the middle sample fill the upper 1.5-sample tail.
        let risk = (losses[2] + 0.5 * losses[1]) / 1.5;
        if risk < best_robust.1 {
            best_robust = (candidate, risk);
        }
    }
    assert_eq!(report.nominal_candidate, best_nominal.0);
    assert_eq!(report.robust_candidate, best_robust.0);
    near(report.robust_band_cvar_per_s, best_robust.1);
    assert!(report.robust_band_cvar_per_s <= report.nominal_band_cvar_per_s);
    // Here physics chooses the widest authored candidate for both objectives.
    assert_eq!(report.nominal_width_scale, 3.0);
    assert_eq!(report.robust_width_scale, 3.0);
    let invalid_band = [nominal, conditions(0.0, 1.0)];
    assert!(matches!(
        robustify_laminar(&study, &WIDTHS, nominal, &invalid_band, 0.5),
        Err(LaminarRobustError::Film {
            candidate: 0,
            fluid: Some(1),
            ..
        })
    ));
    assert!(matches!(
        robustify_laminar(&study, &WIDTHS, nominal, &band, f64::NAN),
        Err(LaminarRobustError::Risk(_))
    ));
}

#[test]
fn laminar_screen_drives_the_original_race_with_native_growth_and_refuses_bad_modes() {
    let study = study();
    let nominal = conditions(1.0, 1.0);
    let losses: Vec<_> = WIDTHS
        .iter()
        .map(|&width| study.evaluate(width, nominal).unwrap().max_growth_per_s)
        .collect();
    let expected = race_base_losses(&losses, 0x7E55).unwrap();
    let actual = screen_laminar_lips(&study, &WIDTHS, nominal, 0x7E55).unwrap();
    assert_eq!(actual, expected);
    assert!(matches!(
        screen_laminar_lips(&study, &[1.0, 0.0], nominal, 0x7E55),
        Err(ScreenError::Film {
            candidate: 1,
            source: LaminarFilmError::InvalidWidthScale
        })
    ));
    assert!(matches!(
        screen_laminar_lips(&study, &[1.0], nominal, 0x7E55),
        Err(ScreenError::TooFewCandidates { count: 1 })
    ));
}
