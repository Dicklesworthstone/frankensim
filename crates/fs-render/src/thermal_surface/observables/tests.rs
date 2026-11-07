use super::*;
use crate::thermal_surface::{
    GrayDiffuseThermalSurface, SpectralDiffuseThermalSurface, SpectralEmissivity,
};

fn gray(temperature: f64, emissivity: f64) -> GrayDiffuseThermalSurface {
    GrayDiffuseThermalSurface::try_new(temperature, emissivity, ContentHash([62; 32])).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual / expected - 1.0).abs() < tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn g1_absolute_blackbody_outputs_match_independent_si_quadrature() {
    // 50-decimal-digit SI Planck integration in metres, Simpson integration
    // with 64 subdivisions per observer bin. The published CIE data remain
    // piecewise constant as declared by observer::xyz_at_bin. These are
    // numerical reference values, not measured filament calibration.
    for (temperature, visible, xyz) in [
        (
            1000.0,
            0.631_247_649_997_230_5,
            [
                5.094_096_889_857_176,
                2.688_840_515_182_997_5,
                0.021_800_178_124_525_768,
            ],
        ),
        (
            2856.0,
            397_161.907_685_263_43,
            [
                21_693_588.659_813_337,
                19_749_837.905_032_407,
                7_030_513.209_178_548,
            ],
        ),
        (
            6500.0,
            49_341_835.356_637_28,
            [
                2_976_262_683.576_671_6,
                3_072_480_379.702_308,
                3_443_299_160.074_424,
            ],
        ),
    ] {
        let source = gray(temperature, 1.0);
        let report = DiffuseThermalSurface::from(source)
            .visible_observables()
            .unwrap();
        close(report.visible_exitance_w_m2(), visible, 3.0e-11);
        for (actual, expected) in report.xyz().into_iter().zip(xyz) {
            close(actual, expected, 3.0e-11);
        }
        assert_eq!(report.luminance_cd_m2(), report.xyz()[1]);
        assert_eq!(
            report.luminous_exitance_lm_m2(),
            PI * report.luminance_cd_m2()
        );
        assert!(report.visible_exitance_w_m2() < source.bolometric_exitance_w_m2());
        let sum: f64 = xyz.iter().sum();
        let chromaticity = report.chromaticity_xy().unwrap();
        close(chromaticity[0], xyz[0] / sum, 3.0e-11);
        close(chromaticity[1], xyz[1] / sum, 3.0e-11);
        assert_eq!(report.surface_identity(), source.identity());
        assert_eq!(report.observer_id(), observer::OBSERVER_ID);
    }
}

#[test]
fn g3_emissivity_scales_power_without_renormalizing_color() {
    let dim = DiffuseThermalSurface::from(gray(2400.0, 0.2))
        .visible_observables()
        .unwrap();
    let bright = DiffuseThermalSurface::from(gray(2400.0, 0.6))
        .visible_observables()
        .unwrap();
    close(
        bright.visible_exitance_w_m2(),
        3.0 * dim.visible_exitance_w_m2(),
        2.0e-14,
    );
    close(
        bright.luminous_exitance_lm_m2(),
        3.0 * dim.luminous_exitance_lm_m2(),
        2.0e-14,
    );
    for (actual, expected) in bright.xyz().into_iter().zip(dim.xyz()) {
        close(actual, 3.0 * expected, 2.0e-14);
    }
    for (actual, expected) in bright
        .chromaticity_xy()
        .unwrap()
        .into_iter()
        .zip(dim.chromaticity_xy().unwrap())
    {
        close(actual, expected, 2.0e-14);
    }
    assert_ne!(bright.surface_identity(), dim.surface_identity());
}

#[test]
fn g3_borrowed_spectral_law_matches_gray_and_ignores_declared_outside_band() {
    let expected = DiffuseThermalSurface::from(gray(3000.0, 0.7))
        .visible_observables()
        .unwrap();
    for curve in [
        SpectralEmissivity::constant(0.7).unwrap(),
        SpectralEmissivity::table_in_band(&[
            (200.0, 0.1),
            (380.0, 0.7),
            (501.3, 0.7),
            (780.0, 0.7),
            (2000.0, 1.0),
        ])
        .unwrap(),
    ] {
        let surface =
            SpectralDiffuseThermalSurface::try_new(3000.0, curve, ContentHash([63; 32])).unwrap();
        let before = surface;
        let view = DiffuseThermalSurface::from(&surface);
        let report = view.visible_observables().unwrap();
        close(
            report.visible_exitance_w_m2(),
            expected.visible_exitance_w_m2(),
            5.0e-14,
        );
        for (actual, expected) in report.xyz().into_iter().zip(expected.xyz()) {
            close(actual, expected, 5.0e-14);
        }
        assert_eq!(surface, before);
        assert_eq!(report.surface_identity(), surface.identity());
        assert_eq!(report, view.visible_observables().unwrap());
    }
}

#[test]
fn g1_narrow_emissivity_peak_retains_physical_power_between_observer_nodes() {
    let lower = 552.999;
    let peak = 553.0;
    let upper = 553.001;
    let curve = SpectralEmissivity::table(&[
        (380.0, 0.0),
        (lower, 0.0),
        (peak, 1.0),
        (upper, 0.0),
        (780.0, 0.0),
    ])
    .unwrap();
    let surface =
        SpectralDiffuseThermalSurface::try_new(2400.0, curve, ContentHash([64; 32])).unwrap();
    for bin in 0..QUAD_BINS {
        let wavelength = LAMBDA_MIN + (bin as f64 + 0.5) * 5.0;
        assert_eq!(surface.radiance_w_m2_sr_nm(wavelength), Ok(0.0));
    }
    let report = DiffuseThermalSurface::from(&surface)
        .visible_observables()
        .unwrap();
    // Over this symmetric 0.002 nm triangle, the smooth Planck factor varies
    // only at second order in the integrated result. Its actual source value
    // times the triangle area is an independent narrow-band limit.
    let radiance_integral = surface.radiance_w_m2_sr_nm(peak).unwrap() * (upper - lower) / 2.0;
    close(
        report.visible_exitance_w_m2(),
        PI * radiance_integral,
        2.0e-9,
    );
    let response = observer::xyz_at_bin(34);
    for (actual, weight) in report.xyz().into_iter().zip(response) {
        close(
            actual,
            LUMINOUS_EFFICACY_LM_W * radiance_integral * weight,
            2.0e-9,
        );
    }
}

#[test]
fn g0_dark_sources_have_no_invented_chromaticity_and_unknown_bands_refuse() {
    for source in [gray(2400.0, 0.0), gray(1.0, 1.0)] {
        let report = DiffuseThermalSurface::from(source)
            .visible_observables()
            .unwrap();
        assert_eq!(report.visible_exitance_w_m2(), 0.0);
        assert_eq!(report.luminous_exitance_lm_m2(), 0.0);
        assert_eq!(report.xyz(), [0.0; 3]);
        assert_eq!(report.chromaticity_xy(), None);
    }
    for band in [(381.0, 780.0), (380.0, 779.0), (8000.0, 14000.0)] {
        for emissivity in [0.0, 0.7] {
            let surface = SpectralDiffuseThermalSurface::try_new(
                2400.0,
                SpectralEmissivity::constant_in_band(emissivity, band).unwrap(),
                ContentHash([65; 32]),
            )
            .unwrap();
            assert_eq!(
                DiffuseThermalSurface::from(&surface).visible_observables(),
                Err(ThermalSurfaceError::Wavelength)
            );
        }
    }
}

#[test]
fn g0_unrepresentable_observables_refuse_instead_of_infinite_luminance() {
    // The finite-band spectral source is representable; its integrated
    // absolute photometric output is not. No broadband T^4 assumption enters.
    let surface = SpectralDiffuseThermalSurface::try_new(
        1.0e303,
        SpectralEmissivity::constant(1.0).unwrap(),
        ContentHash([66; 32]),
    )
    .unwrap();
    assert_eq!(
        DiffuseThermalSurface::from(&surface).visible_observables(),
        Err(ThermalSurfaceError::Unrepresentable)
    );
}
