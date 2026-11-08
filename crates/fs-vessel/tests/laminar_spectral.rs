use fs_cheb::orr_sommerfeld::film::{FilmSpectrumBudget, FilmSpectrumError};
use fs_lbm::film::FilmMaterial;
use fs_vessel::{
    LaminarFilmSection, LaminarFilmSpectralStudy, LaminarGrowthObjective, LaminarPourConditions,
    LaminarSpectralError,
};

fn budget() -> FilmSpectrumBudget {
    FilmSpectrumBudget {
        coarse_order: 16,
        fine_order: 24,
        modes: 1,
        max_qr_sweeps: 256,
        max_svd_sweeps: 64,
        relative_residual_tolerance: 1.0e-8,
        absolute_resolution_tolerance: 1.0e-8,
        relative_resolution_tolerance: 1.0e-7,
    }
}

fn conditions() -> LaminarPourConditions {
    // h=.001 m, Ubar=.001 m/s, Re=1, S=10; width=.1 m.
    LaminarPourConditions {
        volume_flow_m3_s: 1.0e-7,
        material: FilmMaterial {
            density_kg_m3: 1000.0,
            kinematic_viscosity_m2_s: 1.0e-6,
            surface_tension_n_m: 1.0e-5,
        },
    }
}

fn section() -> LaminarFilmSection {
    LaminarFilmSection {
        width_m: 0.1,
        tangential_acceleration_m_s2: 0.003,
        normal_acceleration_m_s2: 0.0,
    }
}

#[test]
fn actual_si_film_retains_finite_wave_decay_and_correct_clock_scaling() {
    let sections = [section()];
    let study = LaminarFilmSpectralStudy {
        sections: &sections,
        wavenumbers_rad_m: &[500.0],
        budget: budget(),
        max_thickness_width_ratio: 0.02,
    };
    let report = study.evaluate(1.0, conditions()).unwrap();
    assert!((report.max_growth_per_s + 0.037_603_889_986_04).abs() < 2.0e-8);
    assert!((report.angular_frequency_rad_s - 1.257_004_847_164_1).abs() < 2.0e-8);
    assert!((report.reynolds_range[0] - 1.0).abs() < 1.0e-12);
    assert!((report.max_kh - 0.5).abs() < 1.0e-12);
    assert!((report.max_thickness_width_ratio - 0.01).abs() < 1.0e-12);

    // h,width -> 10x; Ubar -> Ubar/10; accelerations -> /1000;
    // sigma -> /10. Re, kh, G and S stay fixed; rates scale by /100.
    let scaled_sections = [LaminarFilmSection {
        width_m: 1.0,
        tangential_acceleration_m_s2: 0.000_003,
        normal_acceleration_m_s2: 0.0,
    }];
    let mut scaled_conditions = conditions();
    scaled_conditions.volume_flow_m3_s *= 10.0;
    scaled_conditions.material.surface_tension_n_m /= 10.0;
    let scaled = LaminarFilmSpectralStudy {
        sections: &scaled_sections,
        wavenumbers_rad_m: &[50.0],
        ..study
    }
    .evaluate(1.0, scaled_conditions)
    .unwrap();
    assert!((100.0 * scaled.max_growth_per_s - report.max_growth_per_s).abs() < 2.0e-8);
    assert!(
        (100.0 * scaled.angular_frequency_rad_s - report.angular_frequency_rad_s).abs() < 2.0e-8
    );
}

#[test]
fn spectral_objective_plugs_into_physical_risk_and_refuses_unsupported_samples() {
    let sections = [section()];
    let study = LaminarFilmSpectralStudy {
        sections: &sections,
        wavenumbers_rad_m: &[500.0],
        budget: budget(),
        max_thickness_width_ratio: 0.02,
    };
    fn objective<T: LaminarGrowthObjective>(study: &T) -> f64 {
        study.growth_per_s(1.0, conditions()).unwrap()
    }
    assert!(objective(&study) < 0.0);
    assert!(matches!(
        study.evaluate_with_cancel(1.0, conditions(), || true),
        Err(LaminarSpectralError::Spectrum {
            source: FilmSpectrumError::Cancelled,
            ..
        })
    ));
    assert!(
        LaminarFilmSpectralStudy {
            max_thickness_width_ratio: 0.005,
            ..study
        }
        .evaluate(1.0, conditions())
        .is_err()
    );
    assert!(
        LaminarFilmSpectralStudy {
            wavenumbers_rad_m: &[0.0],
            ..study
        }
        .evaluate(1.0, conditions())
        .is_err()
    );
}
