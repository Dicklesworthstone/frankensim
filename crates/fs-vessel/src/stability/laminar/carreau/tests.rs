use super::*;
use crate::stability::laminar::{LaminarFilmSpectralStudy, LaminarPourConditions};
use fs_lbm::film::FilmMaterial;
use fs_material::fluid_transport::rheology::{CarreauYasuda, ViscosityShiftMode};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

fn law(parameters: [f64; 5]) -> ThermalCarreauYasuda {
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(
            parameters[0],
            parameters[1],
            parameters[2],
            parameters[3],
            parameters[4],
            [0.0, 1.0e5],
        )
        .unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 12_000.0,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
        ViscosityShiftMode::VerticalAndRate,
    )
}
fn budget() -> FilmSpectrumBudget {
    FilmSpectrumBudget {
        coarse_order: 16,
        fine_order: 24,
        modes: 1,
        max_qr_sweeps: 256,
        max_svd_sweeps: 64,
        relative_residual_tolerance: 1.0e-8,
        absolute_resolution_tolerance: 1.0e-7,
        relative_resolution_tolerance: 1.0e-6,
    }
}
fn base_controls() -> CarreauFilmControls {
    let mut controls = CarreauFilmControls::default();
    controls.slit.local.relative_constitutive_tolerance = 1.0e-12;
    controls.slit.local.relative_quadrature_tolerance = 1.0e-10;
    controls.relative_flow_tolerance = 1.0e-10;
    controls
}
fn conditions() -> CarreauPourConditions {
    CarreauPourConditions {
        volume_flow_m3_s: 1.0e-6,
        density_kg_m3: 1000.0,
        surface_tension_n_m: 0.05,
        law: law([0.02, 0.001, 0.01, 0.5, 2.0]),
        temperature_k: 300.0,
        maximum_laminar_shear_rate_per_s: 1.0e5,
    }
}
fn sections() -> [LaminarFilmSection; 1] {
    [LaminarFilmSection {
        width_m: 0.1,
        tangential_acceleration_m_s2: 3.0,
        normal_acceleration_m_s2: 3.0,
    }]
}

#[test]
fn exact_constant_source_uses_its_actual_temperature_and_matches_newtonian_si() {
    let geometry = [LaminarFilmSection {
        width_m: 0.1,
        tangential_acceleration_m_s2: 0.003,
        normal_acceleration_m_s2: 0.0,
    }];
    let data = CarreauPourConditions {
        volume_flow_m3_s: 1.0e-7,
        density_kg_m3: 1000.0,
        surface_tension_n_m: 1.0e-5,
        law: law([0.001, 0.0001, 0.01, 1.0, 1.5]),
        temperature_k: 330.0,
        maximum_laminar_shear_rate_per_s: 1.0e5,
    };
    let study = CarreauFilmSpectralStudy {
        sections: &geometry,
        wavenumbers_rad_m: &[500.0],
        budget: budget(),
        base_controls: base_controls(),
        max_thickness_width_ratio: 0.02,
    };
    let actual = study.evaluate(1.0, data).unwrap();
    let mu = data
        .law
        .evaluate(data.temperature_k, 0.0)
        .unwrap()
        .dynamic_viscosity_pa_s();
    let original = LaminarFilmSpectralStudy {
        sections: &geometry,
        wavenumbers_rad_m: &[500.0],
        budget: budget(),
        max_thickness_width_ratio: 0.02,
    }
    .evaluate(
        1.0,
        LaminarPourConditions {
            volume_flow_m3_s: data.volume_flow_m3_s,
            material: FilmMaterial {
                density_kg_m3: data.density_kg_m3,
                kinematic_viscosity_m2_s: mu / data.density_kg_m3,
                surface_tension_n_m: data.surface_tension_n_m,
            },
        },
    )
    .unwrap();
    assert!((actual.spectrum.max_growth_per_s - original.max_growth_per_s).abs() < 2.0e-7);
    assert!(
        (actual.spectrum.angular_frequency_rad_s - original.angular_frequency_rad_s).abs() < 2.0e-7
    );
    assert!((actual.spectrum.reynolds_range[0] - original.reynolds_range[0]).abs() < 1.0e-10);
}

#[test]
fn original_nonlinear_profile_and_si_clock_match_an_independent_full_pencil() {
    let geometry = sections();
    let study = CarreauFilmSpectralStudy {
        sections: &geometry,
        wavenumbers_rad_m: &[707.322_721_555_680_3],
        budget: budget(),
        base_controls: base_controls(),
        max_thickness_width_ratio: 0.02,
    };
    let report = study.evaluate(1.0, conditions()).unwrap();
    // Independent scalar constitutive inversion and adaptive flux quadrature,
    // then generalized QZ retaining eta explicitly and analytic mu/M derivatives.
    let h = 0.000_565_512_725_393_923_4;
    let clock = 31.269_089_526_808_408;
    assert!((report.thickness_range_m[0] - h).abs() < 2.0e-11);
    assert!((report.spectrum.max_growth_per_s - (-1.361_176_954_124_344 * clock)).abs() < 2.0e-5);
    assert!(
        (report.spectrum.angular_frequency_rad_s - 1.768_299_360_222_106_4 * clock).abs() < 2.0e-5
    );
    assert!((report.spectrum.reynolds_range[0] - 0.5).abs() < 1.0e-9);
    assert!((report.spectrum.max_kh - 0.4).abs() < 1.0e-8);
    assert!(report.max_relative_flow_residual <= base_controls().relative_flow_tolerance);
    assert!(
        report.max_profile_quadrature_relative_change
            <= base_controls().slit.local.relative_quadrature_tolerance
    );
    // h,width -> 10x; U -> /10; g -> /1000; sigma -> /10.
    // Original constitutive time -> 100x because physical shear -> /100.
    // This preserves the complete nonlinear normalized base and spectral pencil.
    let scaled_geometry = [LaminarFilmSection {
        width_m: 1.0,
        tangential_acceleration_m_s2: 0.003,
        normal_acceleration_m_s2: 0.003,
    }];
    let scaled_data = CarreauPourConditions {
        volume_flow_m3_s: 1.0e-5,
        surface_tension_n_m: 0.005,
        law: law([0.02, 0.001, 1.0, 0.5, 2.0]),
        maximum_laminar_shear_rate_per_s: 1000.0,
        ..conditions()
    };
    let scaled = CarreauFilmSpectralStudy {
        sections: &scaled_geometry,
        wavenumbers_rad_m: &[70.732_272_155_568_03],
        ..study
    }
    .evaluate(1.0, scaled_data)
    .unwrap();
    assert!((scaled.thickness_range_m[0] / 10.0 - h).abs() < 2.0e-11);
    assert!(
        (100.0 * scaled.spectrum.max_growth_per_s - report.spectrum.max_growth_per_s).abs()
            < 2.0e-5
    );
    assert!(
        (100.0 * scaled.spectrum.angular_frequency_rad_s - report.spectrum.angular_frequency_rad_s)
            .abs()
            < 2.0e-5
    );
}

#[test]
fn original_rate_support_surface_regularity_and_cancellation_remain_authoritative() {
    let geometry = sections();
    let study = CarreauFilmSpectralStudy {
        sections: &geometry,
        wavenumbers_rad_m: &[700.0],
        budget: budget(),
        base_controls: base_controls(),
        max_thickness_width_ratio: 0.02,
    };
    let singular = CarreauPourConditions {
        law: law([0.02, 0.001, 0.01, 0.5, 1.5]),
        ..conditions()
    };
    assert!(matches!(
        study.evaluate(1.0, singular),
        Err(CarreauFilmSpectralError::UnsupportedSurfaceRegularity { .. })
    ));
    let limited = CarreauPourConditions {
        maximum_laminar_shear_rate_per_s: 5.0,
        ..conditions()
    };
    assert!(matches!(
        study.evaluate(1.0, limited),
        Err(CarreauFilmSpectralError::Base {
            source: CarreauFilmError::OutsideRateSupport,
            ..
        })
    ));
    assert!(matches!(
        study.evaluate_with_cancel(1.0, conditions(), || true),
        Err(CarreauFilmSpectralError::Spectrum {
            source: FilmSpectrumError::Cancelled,
            ..
        })
    ));
}
