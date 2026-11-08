use super::*;

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

#[test]
fn long_waves_recover_nusselt_speed_and_yih_mean_reynolds_threshold() {
    let solver = PreparedFilmSpectrum::new(budget(), || false).unwrap();
    let alpha = 0.001;
    // g_n/g_t=1 implies G=3/Re. Re_c=5/6 in mean-speed convention.
    for re in [0.7, 1.0] {
        let parameters = FilmParameters {
            reynolds: re,
            kh: alpha,
            normal_gravity: 3.0 / re,
            capillarity: 0.0,
        };
        let lambda = solver.solve(parameters, || false).unwrap()[0].eigenvalue;
        let leading_growth = alpha * alpha * (1.2 * re - 1.0);
        assert!((lambda.re - leading_growth).abs() < 1.0e-3 * leading_growth.abs());
        assert!((lambda.im / alpha + 3.0).abs() < 2.0e-5);
    }
}

#[test]
fn finite_wavelength_and_finite_reynolds_match_independent_generalized_qz() {
    let solver = PreparedFilmSpectrum::new(budget(), || false).unwrap();
    // Independent generalized-QZ evaluation of the original generalized pencil
    // at degrees 24,32,40. These are numerical references, not experiment.
    for (re, alpha, capillarity, expected) in [
        (
            1.0,
            0.5,
            10.0,
            C64::new(-0.037_603_889_986_04, -1.257_004_847_164_1),
        ),
        (
            100.0,
            0.5,
            1.0,
            C64::new(0.030_564_725_685_8, -0.902_424_917_853_4),
        ),
    ] {
        let actual = solver
            .solve(
                FilmParameters {
                    reynolds: re,
                    kh: alpha,
                    normal_gravity: 0.0,
                    capillarity,
                },
                || false,
            )
            .unwrap()[0]
            .eigenvalue;
        assert!(
            (actual - expected).abs() < 2.0e-8,
            "{actual:?} != {expected:?}"
        );
    }
    // At Re=1, kh=.5, S=10 the actual mode decays, although extending the
    // Benney polynomial to this wavelength gives positive growth.
    let extrapolated_benney = 0.5 * 0.5 * (1.2 - 10.0 * 0.5 * 0.5 / 3.0);
    assert!(extrapolated_benney > 0.0);
}

#[test]
fn viscous_shear_modes_recover_the_free_surface_half_integer_spectrum() {
    let mut options = budget();
    options.modes = 3;
    let solver = PreparedFilmSpectrum::new(options, || false).unwrap();
    let modes = solver
        .solve(
            FilmParameters {
                reynolds: 1.0,
                kh: 0.0001,
                normal_gravity: 0.0,
                capillarity: 0.0,
            },
            || false,
        )
        .unwrap();
    // At k -> 0 the decaying velocity modes have u(0)=0, Du(1)=0:
    // lambda_j = -[(j+1/2)*pi]^2/Re, alongside the neutral surface mode.
    let pi2 = core::f64::consts::PI * core::f64::consts::PI;
    assert!((modes[1].eigenvalue.re + 0.25 * pi2).abs() < 1.0e-5);
    assert!((modes[2].eigenvalue.re + 2.25 * pi2).abs() < 1.0e-5);
}

#[test]
fn underresolution_invalid_inputs_and_cancellation_return_no_spectrum() {
    assert!(matches!(
        PreparedFilmSpectrum::new(budget(), || true),
        Err(FilmSpectrumError::Cancelled)
    ));
    let solver = PreparedFilmSpectrum::new(budget(), || false).unwrap();
    let data = FilmParameters {
        reynolds: 100.0,
        kh: 0.5,
        normal_gravity: 0.0,
        capillarity: 1.0,
    };
    assert_eq!(
        solver.solve(data, || true),
        Err(FilmSpectrumError::Cancelled)
    );
    assert!(matches!(
        solver.solve(FilmParameters { kh: 0.0, ..data }, || false),
        Err(FilmSpectrumError::Invalid(_))
    ));
    let mut unresolved = budget();
    unresolved.coarse_order = 8;
    unresolved.fine_order = 12;
    unresolved.absolute_resolution_tolerance = 1.0e-12;
    unresolved.relative_resolution_tolerance = 0.0;
    assert!(matches!(
        PreparedFilmSpectrum::new(unresolved, || false)
            .unwrap()
            .solve(data, || false),
        Err(FilmSpectrumError::Resolution { .. })
    ));
}
