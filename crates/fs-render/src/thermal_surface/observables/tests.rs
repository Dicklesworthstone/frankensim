use super::*;
use crate::motion::{ShotTimeBounds, ShutterConvention, ShutterDistribution};
use crate::thermal_surface::{
    GrayDiffuseThermalSurface, SpectralDiffuseThermalSurface, SpectralEmissivity,
};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 62,
                kernel_id: 0x5448_5649,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn shutter(start: f64, duration: f64) -> ShutterInterval {
    ShutterInterval::resolve(
        start,
        duration,
        ShutterConvention::FrontLoaded,
        ShutterDistribution::UniformCounterV1,
        ShotTimeBounds::try_new(start, start + duration).unwrap(),
    )
    .unwrap()
}

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

#[test]
fn g1_nonuniform_patches_keep_absolute_area_power_and_hot_spot_color() {
    with_cx(|cx, _| {
        let patches = [
            ThermalEmissionPatch {
                surface: gray(1000.0, 0.8).into(),
                emitting_area_m2: 1e-4,
            },
            ThermalEmissionPatch {
                surface: gray(2856.0, 0.8).into(),
                emitting_area_m2: 3e-6,
            },
        ];
        let observed = ThermalVisibleEmission::observe(cx, &patches, 2).unwrap();
        // Independent SI quadrature values retained by the single-source G1
        // fixture above. Each source contributes before any spatial averaging.
        close(
            observed.visible_power_w(),
            0.8 * (1e-4 * 0.631_247_649_997_230_5 + 3e-6 * 397_161.907_685_263_43),
            3e-11,
        );
        let cold_xyz = [
            5.094_096_889_857_176,
            2.688_840_515_182_997_5,
            0.021_800_178_124_525_768,
        ];
        let hot_xyz = [
            21_693_588.659_813_337,
            19_749_837.905_032_407,
            7_030_513.209_178_548,
        ];
        let expected = core::array::from_fn::<_, 3, _>(|i| {
            PI * 0.8 * (1e-4 * cold_xyz[i] + 3e-6 * hot_xyz[i])
        });
        for (actual, expected) in observed.xyz_flux_lm().into_iter().zip(expected) {
            close(actual, expected, 3e-11);
        }
        assert_eq!(observed.luminous_flux_lm(), observed.xyz_flux_lm()[1]);
        let averaged_temperature = (1e-4 * 1000.0 + 3e-6 * 2856.0) / (1e-4 + 3e-6);
        let averaged = ThermalVisibleEmission::observe(
            cx,
            &[ThermalEmissionPatch {
                surface: gray(averaged_temperature, 0.8).into(),
                emitting_area_m2: 1e-4 + 3e-6,
            }],
            1,
        )
        .unwrap();
        assert!(observed.visible_power_w() > 100.0 * averaged.visible_power_w());
        assert_ne!(observed.chromaticity_xy(), averaged.chromaticity_xy());
        let both_sides = patches.map(|mut patch| {
            patch.emitting_area_m2 *= 2.0;
            patch
        });
        let doubled = ThermalVisibleEmission::observe(cx, &both_sides, 2).unwrap();
        close(
            doubled.visible_power_w(),
            2.0 * observed.visible_power_w(),
            1e-14,
        );
        close(
            doubled.luminous_flux_lm(),
            2.0 * observed.luminous_flux_lm(),
            1e-14,
        );
        assert_eq!(doubled.chromaticity_xy(), observed.chromaticity_xy());
        assert_ne!(doubled.identity(), observed.identity());
        let split_hot = ThermalEmissionPatch {
            emitting_area_m2: 1.5e-6,
            ..patches[1]
        };
        let split =
            ThermalVisibleEmission::observe(cx, &[patches[0], split_hot, split_hot], 3).unwrap();
        close(split.visible_power_w(), observed.visible_power_w(), 1e-14);
        assert_eq!(observed.observer_id(), observer::OBSERVER_ID);
        let empty = ThermalVisibleEmission::observe(cx, &[], 0).unwrap();
        assert_eq!(empty.visible_power_w(), 0.0);
        assert_eq!(empty.chromaticity_xy(), None);
    });
}

#[test]
fn g1_finite_exposure_clips_actual_intervals_and_preserves_pulse_energy() {
    with_cx(|cx, _| {
        let observe = |temperature| {
            ThermalVisibleEmission::observe(
                cx,
                &[ThermalEmissionPatch {
                    surface: gray(temperature, 0.8).into(),
                    emitting_area_m2: 1e-4,
                }],
                1,
            )
            .unwrap()
        };
        let cold = observe(1000.0);
        let hot = observe(2856.0);
        // Manufactured held observations on an absolute clock; this fixture
        // verifies observation integration, not instantaneous physical heating.
        let history = [
            HeldThermalEmissionInterval {
                start_s: 12.0,
                end_s: 13.0,
                emission: cold,
            },
            HeldThermalEmissionInterval {
                start_s: 13.0,
                end_s: 13.25,
                emission: hot,
            },
            HeldThermalEmissionInterval {
                start_s: 13.25,
                end_s: 14.0,
                emission: cold,
            },
        ];
        let exposure =
            ThermalVisibleExposure::integrate(cx, shutter(12.5, 1.0), &history, 3).unwrap();
        close(
            exposure.visible_energy_j(),
            0.75 * cold.visible_power_w() + 0.25 * hot.visible_power_w(),
            1e-14,
        );
        close(
            exposure.luminous_energy_lm_s(),
            0.75 * cold.luminous_flux_lm() + 0.25 * hot.luminous_flux_lm(),
            1e-14,
        );
        assert_eq!(exposure.mean_visible_power_w(), exposure.visible_energy_j());
        assert_eq!(
            exposure.mean_luminous_flux_lm(),
            exposure.luminous_energy_lm_s()
        );
        assert_eq!(exposure.shutter().open_s(), 12.5);
        let expected_xyz = core::array::from_fn::<_, 3, _>(|i| {
            0.75 * cold.xyz_flux_lm()[i] + 0.25 * hot.xyz_flux_lm()[i]
        });
        for (actual, expected) in exposure.xyz_lm_s().into_iter().zip(expected_xyz) {
            close(actual, expected, 1e-14);
        }
        let left = ThermalVisibleExposure::integrate(cx, shutter(12.5, 0.5), &history, 3).unwrap();
        let right = ThermalVisibleExposure::integrate(cx, shutter(13.0, 0.5), &history, 3).unwrap();
        close(
            left.visible_energy_j() + right.visible_energy_j(),
            exposure.visible_energy_j(),
            1e-14,
        );
        let split = [
            history[0],
            HeldThermalEmissionInterval {
                end_s: 13.125,
                ..history[1]
            },
            HeldThermalEmissionInterval {
                start_s: 13.125,
                ..history[1]
            },
            history[2],
        ];
        let partitioned =
            ThermalVisibleExposure::integrate(cx, shutter(12.5, 1.0), &split, 4).unwrap();
        close(
            partitioned.visible_energy_j(),
            exposure.visible_energy_j(),
            1e-14,
        );
        for (actual, expected) in partitioned
            .chromaticity_xy()
            .unwrap()
            .into_iter()
            .zip(exposure.chromaticity_xy().unwrap())
        {
            close(actual, expected, 1e-14);
        }
        let short =
            ThermalVisibleExposure::integrate(cx, shutter(13.0, 0.125), &history, 3).unwrap();
        close(short.mean_luminous_flux_lm(), hot.luminous_flux_lm(), 1e-14);
        assert_eq!(short.observer_id(), observer::OBSERVER_ID);
        assert_eq!(hot, history[1].emission);
    });
}

#[test]
fn g0_emission_observation_refuses_unknown_coverage_bad_work_and_cancellation() {
    with_cx(|cx, gate| {
        let patch = ThermalEmissionPatch {
            surface: gray(2856.0, 0.8).into(),
            emitting_area_m2: 1e-4,
        };
        let accepted = ThermalVisibleEmission::observe(cx, &[patch], 1).unwrap();
        assert_eq!(
            ThermalVisibleEmission::observe(cx, &[patch], 0),
            Err(ThermalEmissionError::WorkLimit)
        );
        for area in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert_eq!(
                ThermalVisibleEmission::observe(
                    cx,
                    &[ThermalEmissionPatch {
                        emitting_area_m2: area,
                        ..patch
                    }],
                    1
                ),
                Err(ThermalEmissionError::InvalidArea)
            );
        }
        assert_eq!(
            ThermalVisibleEmission::observe(
                cx,
                &[ThermalEmissionPatch {
                    emitting_area_m2: f64::MAX,
                    ..patch
                }],
                1
            ),
            Err(ThermalEmissionError::NumericalRange)
        );
        let incomplete = SpectralDiffuseThermalSurface::try_new(
            2856.0,
            SpectralEmissivity::constant_in_band(0.0, (400.0, 700.0)).unwrap(),
            ContentHash([67; 32]),
        )
        .unwrap();
        assert_eq!(
            ThermalVisibleEmission::observe(
                cx,
                &[ThermalEmissionPatch {
                    surface: (&incomplete).into(),
                    emitting_area_m2: 1e-4,
                }],
                1
            ),
            Err(ThermalEmissionError::Surface(
                ThermalSurfaceError::Wavelength
            ))
        );
        let row = HeldThermalEmissionInterval {
            start_s: 0.0,
            end_s: 1.0,
            emission: accepted,
        };
        let history = [
            row,
            HeldThermalEmissionInterval {
                start_s: 1.0,
                end_s: 2.0,
                emission: accepted,
            },
        ];
        for second_start in [0.5, 1.5, f64::NAN] {
            let bad = [
                row,
                HeldThermalEmissionInterval {
                    start_s: second_start,
                    ..history[1]
                },
            ];
            assert_eq!(
                ThermalVisibleExposure::integrate(cx, shutter(0.0, 2.0), &bad, 2),
                Err(ThermalEmissionError::InvalidHistory)
            );
        }
        assert_eq!(
            ThermalVisibleExposure::integrate(cx, shutter(0.0, 2.0), &[row], 1),
            Err(ThermalEmissionError::IncompleteExposure)
        );
        assert_eq!(
            ThermalVisibleExposure::integrate(cx, shutter(0.0, 2.0), &history, 1),
            Err(ThermalEmissionError::WorkLimit)
        );
        assert_eq!(
            ThermalVisibleExposure::integrate(cx, shutter(0.0, 0.0), &history, 2),
            Err(ThermalEmissionError::InvalidHistory)
        );
        gate.request();
        assert_eq!(
            ThermalVisibleEmission::observe(cx, &[patch], 1),
            Err(ThermalEmissionError::Cancelled)
        );
        assert_eq!(
            ThermalVisibleExposure::integrate(cx, shutter(0.0, 2.0), &history, 2),
            Err(ThermalEmissionError::Cancelled)
        );
        assert_eq!(history[0].emission, accepted);
    });
}
