//! G0/G1/G3/G4: analytic cell transport and explicit admission boundaries.
use super::*;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 3,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    })
}

fn budget() -> AbsorbingGridBudget {
    AbsorbingGridBudget {
        max_grid_cells: 4096,
        max_steps_per_segment: 100,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1e-100),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_heterogeneous_beer_lambert_and_partial_inside_segment() {
    with_cx(|cx, _| {
        let density = [2.0, 3.0, 1.0];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 3], &density, [0.0; 3], [1.0, 1.0, 0.2]),
            MassAbsorptionSpectrum::constant(0.7).unwrap(),
            None,
            budget(),
        )
        .unwrap();
        assert!(core::ptr::eq(medium.grid.data.as_ptr(), density.as_ptr()));
        for (z, direction) in [(-1.0, 1.0), (1.0, -1.0)] {
            let t = medium
                .transmittance(
                    cx,
                    Point3::new(0.25, 0.25, z),
                    Vec3::new(0.0, 0.0, direction),
                    f64::INFINITY,
                    &[380.0, 550.0, 780.0],
                )
                .unwrap();
            for value in t {
                close(value, (-0.84_f64).exp());
            }
        }
        let t = medium
            .transmittance(
                cx,
                Point3::new(0.25, 0.25, 0.1),
                Vec3::new(0.0, 0.0, 1.0),
                0.35,
                &[550.0],
            )
            .unwrap();
        close(t[0], (-0.595_f64).exp());
    });
}

#[test]
fn g1_emission_uses_physical_planck_with_temperature_ordering() {
    with_cx(|cx, _| {
        let density = [1.0, 1.0];
        let temperatures = [1400.0, 2000.0];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 2], &density, [0.0; 3], [1.0; 3]),
            MassAbsorptionSpectrum::constant(0.5).unwrap(),
            Some(&temperatures),
            budget(),
        )
        .unwrap();
        assert!(medium.has_thermal_emission());
        let t = (-0.5_f64).exp();
        let a = planck_radiance_w_m2_sr_nm(1400.0, 550.0).unwrap();
        let b = planck_radiance_w_m2_sr_nm(2000.0, 550.0).unwrap();
        for (z, d, near, far) in [(-1.0, 1.0, a, b), (3.0, -1.0, b, a)] {
            let result = medium
                .transport(
                    cx,
                    Point3::new(0.25, 0.25, z),
                    Vec3::new(0.0, 0.0, d),
                    f64::INFINITY,
                    &[550.0],
                )
                .unwrap();
            close(result.transmittance[0], t * t);
            close(result.emission_w_m2_sr_nm[0], (1.0 - t) * (near + t * far));
        }
    });
}

#[test]
fn g3_density_length_scaling_and_cell_refinement() {
    with_cx(|cx, _| {
        for (cells, spacing, rho) in [(1, 1.0, 2.0), (4, 0.25, 2.0), (4, 0.5, 1.0)] {
            let density = vec![rho; cells * cells * cells];
            let medium = AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([cells; 3], &density, [0.0; 3], [spacing; 3]),
                MassAbsorptionSpectrum::constant(0.3).unwrap(),
                None,
                budget(),
            )
            .unwrap();
            let d = Vec3::new(1.0, 1.0, 1.0).scale(1.0 / 3.0_f64.sqrt());
            let result = medium
                .transmittance(
                    cx,
                    Point3::new(-1.0, -1.0, -1.0),
                    d,
                    f64::INFINITY,
                    &[550.0],
                )
                .unwrap();
            close(result[0], (-0.6 * 3.0_f64.sqrt()).exp());
        }
    });
}

#[test]
fn g0_boundaries_parallel_rays_and_spectral_interpolation() {
    with_cx(|cx, _| {
        let density = [1.0; 8];
        let spectrum =
            MassAbsorptionSpectrum::table(&[(380.0, 0.0), (580.0, 1.0), (780.0, 3.0)]).unwrap();
        assert_eq!(spectrum.at(480.0).unwrap(), 0.5);
        assert_eq!(spectrum.at(680.0).unwrap(), 2.0);
        assert_eq!(spectrum.at(800.0), Err(VolumeTransportError::Wavelength));
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([2; 3], &density, [0.0; 3], [0.5; 3]),
            spectrum,
            None,
            budget(),
        )
        .unwrap();
        for (z, d, optical_length) in [
            (0.0, 1.0, 1.0_f64),
            (1.0, -1.0, 1.0),
            (0.5, -1.0, 0.5),
            (0.5, 1.0, 0.5),
            (1.0, 1.0, 0.0),
            (0.0, -1.0, 0.0),
        ] {
            let result = medium
                .transmittance(
                    cx,
                    Point3::new(0.5, 0.5, z),
                    Vec3::new(0.0, 0.0, d),
                    f64::INFINITY,
                    &[580.0],
                )
                .unwrap();
            close(result[0], (-optical_length).exp());
        }
        for x in [-0.1, 1.0] {
            assert_eq!(
                medium
                    .transmittance(
                        cx,
                        Point3::new(x, 0.5, -1.0),
                        Vec3::new(0.0, 0.0, 1.0),
                        f64::INFINITY,
                        &[580.0]
                    )
                    .unwrap(),
                [1.0]
            );
        }
    });
}

#[test]
fn g1_dda_matches_independent_cell_intersection_integrals() {
    with_cx(|cx, _| {
        let dims = [3, 2, 4];
        let lo = [-0.4, 0.7, -1.2];
        let cell = [0.3, 0.45, 0.22];
        let density: Vec<_> = (0..24).map(|i| 0.2 * (i + 1) as f64).collect();
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new(dims, &density, lo, cell),
            MassAbsorptionSpectrum::constant(0.3).unwrap(),
            None,
            budget(),
        )
        .unwrap();
        let center = [0.05, 1.15, -0.76];
        for ray in 0..64 {
            let a = ray as f64 * 0.371;
            let radius = if ray % 3 == 0 { 0.15 } else { 2.0 };
            let o = [
                center[0] + radius * a.cos(),
                center[1] + radius * a.sin(),
                center[2] + radius * (0.723 * a).sin(),
            ];
            let raw = Vec3::new(center[0] - o[0], center[1] - o[1], center[2] - o[2]);
            let dir = raw.scale(1.0 / raw.norm());
            let d = [dir.x, dir.y, dir.z];
            let distance = if ray % 4 == 0 {
                raw.norm()
            } else {
                f64::INFINITY
            };
            let mut column_density = 0.0;
            for k in 0..dims[2] {
                for j in 0..dims[1] {
                    for i in 0..dims[0] {
                        let idx = [i, j, k];
                        let (mut enter, mut exit) = (0.0_f64, distance);
                        for axis in 0..3 {
                            let lower = lo[axis] + idx[axis] as f64 * cell[axis];
                            let upper = lo[axis] + (idx[axis] + 1) as f64 * cell[axis];
                            if d[axis] == 0.0 {
                                if o[axis] < lower || o[axis] >= upper {
                                    exit = 0.0;
                                }
                            } else {
                                let x = (lower - o[axis]) / d[axis];
                                let y = (upper - o[axis]) / d[axis];
                                enter = enter.max(x.min(y));
                                exit = exit.min(x.max(y));
                            }
                        }
                        column_density +=
                            (exit - enter).max(0.0) * density[i + dims[0] * (j + dims[1] * k)];
                    }
                }
            }
            let result = medium
                .transmittance(cx, Point3::new(o[0], o[1], o[2]), dir, distance, &[550.0])
                .unwrap();
            close(result[0], (-0.3 * column_density).exp());
        }
    });
}

#[test]
fn g1_optically_thin_emission_and_zero_absorption_identity() {
    with_cx(|cx, _| {
        let density = [1.0];
        let temperature = [1800.0];
        for coefficient in [0.0, 1e-20] {
            let medium = AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [0.0; 3], [1.0; 3]),
                MassAbsorptionSpectrum::constant(coefficient).unwrap(),
                Some(&temperature),
                budget(),
            )
            .unwrap();
            let result = medium
                .transport(
                    cx,
                    Point3::new(0.5, 0.5, -1.0),
                    Vec3::new(0.0, 0.0, 1.0),
                    f64::INFINITY,
                    &[550.0],
                )
                .unwrap();
            assert_eq!(medium.is_vacuum(), coefficient == 0.0);
            if coefficient == 0.0 {
                assert_eq!(
                    result,
                    SegmentTransport {
                        transmittance: [1.0],
                        emission_w_m2_sr_nm: [0.0]
                    }
                );
            } else {
                close(
                    result.emission_w_m2_sr_nm[0],
                    coefficient * planck_radiance_w_m2_sr_nm(1800.0, 550.0).unwrap(),
                );
            }
        }
    });
}

#[test]
fn g0_g4_admission_identity_budgets_and_cancellation() {
    with_cx(|cx, gate| {
        let absorption = MassAbsorptionSpectrum::constant(1.0).unwrap();
        let data = [1.0, 2.0];
        let make = |data: &'static [f64], temperature, steps| {
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1, 1, 2], data, [0.0; 3], [1.0; 3]),
                absorption,
                temperature,
                AbsorbingGridBudget {
                    max_steps_per_segment: steps,
                    ..budget()
                },
            )
        };
        assert!(matches!(
            make(&[-1.0, 2.0], None, 3),
            Err(VolumeTransportError::Density { index: 0 })
        ));
        assert!(matches!(
            make(&[1.0, 2.0], Some(&[300.0]), 3),
            Err(VolumeTransportError::Temperature)
        ));
        let medium = make(&[1.0, 2.0], None, 1).unwrap();
        assert_ne!(
            medium.content_identity(),
            make(&[1.0, 3.0], None, 1).unwrap().content_identity()
        );
        assert_eq!(
            medium.content_identity(),
            make(&[1.0, 2.0], None, 5).unwrap().content_identity()
        );
        let trace = || {
            medium.transport(
                cx,
                Point3::new(0.5, 0.5, -1.0),
                Vec3::new(0.0, 0.0, 1.0),
                f64::INFINITY,
                &[550.0],
            )
        };
        assert_eq!(trace(), Err(VolumeTransportError::Budget));
        assert!(matches!(
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1, 1, 2], &data, [0.0; 3], [1.0; 3]),
                absorption,
                None,
                AbsorbingGridBudget {
                    max_grid_cells: 1,
                    ..budget()
                }
            ),
            Err(VolumeTransportError::Budget)
        ));
        gate.request();
        assert_eq!(trace(), Err(VolumeTransportError::Cancelled));
    });
}

#[test]
fn g0_broadband_sources_retain_exact_support_and_source_identity() {
    let knots = [(250.0, 2.0), (1064.0, 0.2), (10_000.0, 1.0)];
    let source = ContentHash([7; 32]);
    let spectrum = MassAbsorptionSpectrum::source_table(&knots, source).unwrap();
    assert_eq!(spectrum.wavelength_bounds_nm(), (250.0, 10_000.0));
    assert_eq!(spectrum.source_identity(), Some(source));
    assert_eq!(spectrum.at(1064.0), Ok(0.2));
    close(spectrum.at(5532.0).unwrap(), 0.6);
    for wavelength in [0.0, -1.0, 249.0, 10_001.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            spectrum.at(wavelength),
            Err(VolumeTransportError::Wavelength)
        );
    }
    assert_eq!(
        MassAbsorptionSpectrum::source_table(&knots, ContentHash([0; 32])),
        Err(VolumeTransportError::SourceIdentity)
    );
    for invalid in [
        [(0.0, 1.0), (1000.0, 1.0)],
        [(1000.0, 1.0), (250.0, 1.0)],
        [(250.0, -1.0), (1000.0, 1.0)],
    ] {
        assert_eq!(
            MassAbsorptionSpectrum::source_table(&invalid, source),
            Err(VolumeTransportError::AbsorptionSpectrum)
        );
    }
    assert_eq!(
        MassAbsorptionSpectrum::table(&knots),
        Err(VolumeTransportError::AbsorptionSpectrum)
    );
    assert_eq!(
        MassAbsorptionSpectrum::constant(1.0).unwrap().at(1064.0),
        Err(VolumeTransportError::Wavelength)
    );
    with_cx(|cx, _| {
        let density = [1.0];
        let make = |identity| {
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [0.0; 3], [1.0; 3]),
                MassAbsorptionSpectrum::source_table(&knots, identity).unwrap(),
                None,
                budget(),
            )
            .unwrap()
        };
        assert_ne!(
            make(source).content_identity(),
            make(ContentHash([8; 32])).content_identity()
        );
        let medium = make(source);
        // Even a zero-power beam cannot manufacture support outside its source.
        assert_eq!(
            medium.deposit_beam_power(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                3.0,
                &[10_001.0],
                &[0.0],
                0,
            ),
            Err(VolumeTransportError::Wavelength)
        );
    });
}

#[test]
fn g1_broadband_absorption_and_cold_infrared_emission_follow_si_laws() {
    with_cx(|cx, _| {
        let density = [2.0];
        let temperature = [20.0];
        let wavelengths = [250.0, 1064.0, 150_000.0];
        let coefficients: [f64; 3] = [1.0, 0.1, 0.5];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [0.0; 3], [0.5; 3]),
            MassAbsorptionSpectrum::source_table(
                &[(250.0, 1.0), (1064.0, 0.1), (150_000.0, 0.5)],
                ContentHash([5; 32]),
            )
            .unwrap(),
            Some(&temperature),
            budget(),
        )
        .unwrap();
        // A 20 K body's visible Planck source underflows; its admitted far-IR
        // source still has a direction proposal and physical spectral power.
        assert!(medium.has_emission_proposal());
        let transport = medium
            .transport(
                cx,
                Point3::new(-1.0, 0.25, 0.25),
                Vec3::new(1.0, 0.0, 0.0),
                2.0,
                &wavelengths,
            )
            .unwrap();
        for lane in 0..3 {
            let transmitted = (-coefficients[lane]).exp();
            close(transport.transmittance[lane], transmitted);
            let lambda_m = wavelengths[lane] * 1e-9;
            // Independent direct SI Planck law, converted from per metre to per nm.
            let h = 6.626_070_15e-34;
            let c = 299_792_458.0;
            let k = 1.380_649e-23;
            let b = 2.0 * h * c * c * 1e-9
                / (lambda_m.powi(5) * (h * c / (lambda_m * k * temperature[0])).exp_m1());
            close(transport.emission_w_m2_sr_nm[lane], b * (1.0 - transmitted));
        }
        assert!(transport.emission_w_m2_sr_nm[2] > 0.0);
        let deposition = medium
            .deposit_beam_power(
                cx,
                Point3::new(-1.0, 0.25, 0.25),
                Vec3::new(1.0, 0.0, 0.0),
                2.0,
                &wavelengths,
                &[10.0; 3],
                1,
            )
            .unwrap();
        assert_eq!(deposition.wavelength_bounds_nm, (250.0, 150_000.0));
        for lane in 0..3 {
            close(
                deposition.absorbed_power_w[lane],
                10.0 * (1.0 - (-coefficients[lane]).exp()),
            );
            close(
                deposition.transmitted_power_w[lane] + deposition.absorbed_power_w[lane],
                10.0,
            );
        }
    });
}

#[test]
fn g1_broadband_scattering_separates_heat_and_refuses_missing_optical_support() {
    use super::scattering::{MassScatteringSpectrum, PhaseFunction, ScatteringGrid};
    with_cx(|cx, _| {
        let density = [1.0];
        let source = ContentHash([9; 32]);
        let make = || {
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [0.0; 3], [1.0; 3]),
                MassAbsorptionSpectrum::source_table(&[(250.0, 1.0), (10_000.0, 1.0)], source)
                    .unwrap(),
                None,
                budget(),
            )
            .unwrap()
        };
        let scattering =
            MassScatteringSpectrum::source_table(&[(250.0, 3.0), (10_000.0, 3.0)], source).unwrap();
        assert!(matches!(
            ScatteringGrid::try_new(
                cx,
                make(),
                MassScatteringSpectrum::constant(0.0).unwrap(),
                PhaseFunction::Isotropic
            ),
            Err(VolumeTransportError::Wavelength)
        ));
        let visible_phase = PhaseFunction::spectral_table(&[
            (380.0, PhaseFunction::Isotropic),
            (780.0, PhaseFunction::Rayleigh),
        ])
        .unwrap();
        assert!(matches!(
            ScatteringGrid::try_new(cx, make(), scattering, visible_phase),
            Err(VolumeTransportError::Wavelength)
        ));
        let medium =
            ScatteringGrid::try_new(cx, make(), scattering, PhaseFunction::Isotropic).unwrap();
        let deposition = medium
            .deposit_beam_power(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                3.0,
                &[250.0, 10_000.0],
                &[8.0; 2],
                1,
            )
            .unwrap();
        for lane in 0..2 {
            close(
                deposition.absorbed_power_w[lane],
                2.0 * (1.0 - (-4.0_f64).exp()),
            );
            close(
                deposition.scattered_power_w[lane],
                6.0 * (1.0 - (-4.0_f64).exp()),
            );
            close(deposition.transmitted_power_w[lane], 8.0 * (-4.0_f64).exp());
        }
    });
}

#[test]
fn g1_observed_cell_transport_retains_depth_work_and_ordered_thermal_radiance() {
    with_cx(|cx, _| {
        let density = [1.0, 2.0, 1.0];
        let temperature = [1200.0, 1700.0, 900.0];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([3, 1, 1], &density, [0.0; 3], [1.0; 3]),
            MassAbsorptionSpectrum::table(&[(380.0, 0.25), (780.0, 0.5)]).unwrap(),
            Some(&temperature),
            budget(),
        )
        .unwrap();
        let wavelengths = [380.0, 780.0];
        for (x, sign, cells) in [(-1.0, 1.0, [0, 1, 2]), (4.0, -1.0, [2, 1, 0])] {
            let origin = Point3::new(x, 0.5, 0.5);
            let direction = Vec3::new(sign, 0.0, 0.0);
            let result = medium
                .transport_observed::<2, true>(
                    cx,
                    origin,
                    direction,
                    f64::INFINITY,
                    &wavelengths,
                    3,
                )
                .unwrap();
            let original = medium
                .transport(cx, origin, direction, f64::INFINITY, &wavelengths)
                .unwrap();
            assert_eq!(result.transport, original);
            assert_eq!(result.cell_visits, 3);
            assert_eq!(result.optical_depth, [1.0, 2.0]);
            let shadow = medium
                .transport_observed::<2, false>(
                    cx,
                    origin,
                    direction,
                    f64::INFINITY,
                    &wavelengths,
                    3,
                )
                .unwrap();
            assert_eq!(shadow.optical_depth, result.optical_depth);
            assert_eq!(
                shadow.transport.transmittance,
                result.transport.transmittance
            );
            assert_eq!(shadow.transport.emission_w_m2_sr_nm, [0.0; 2]);
            for lane in 0..2 {
                let kappa = [0.25, 0.5][lane];
                let mut prefix = 1.0;
                let mut expected = 0.0;
                for cell in cells {
                    let tau = density[cell] * kappa;
                    let source =
                        planck_radiance_w_m2_sr_nm(temperature[cell], wavelengths[lane]).unwrap();
                    expected += prefix * -(-tau).exp_m1() * source;
                    prefix *= (-tau).exp();
                }
                close(result.transport.emission_w_m2_sr_nm[lane], expected);
            }
        }
    });
}

#[test]
fn g4_observed_cells_respect_native_and_remaining_budgets_and_cancellation() {
    with_cx(|cx, gate| {
        let density = [1.0; 2];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([2, 1, 1], &density, [0.0; 3], [1.0; 3]),
            MassAbsorptionSpectrum::constant(1.0).unwrap(),
            None,
            AbsorbingGridBudget {
                max_grid_cells: 2,
                max_steps_per_segment: 1,
            },
        )
        .unwrap();
        let origin = Point3::new(-1.0, 0.5, 0.5);
        let direction = Vec3::new(1.0, 0.0, 0.0);
        for cap in [0, 1, 2, u32::MAX] {
            assert_eq!(
                medium.transport_observed::<1, false>(
                    cx,
                    origin,
                    direction,
                    f64::INFINITY,
                    &[550.0],
                    cap,
                ),
                Err(VolumeTransportError::Budget)
            );
        }
        let short = medium
            .transport_observed::<1, false>(cx, origin, direction, 1.5, &[550.0], 1)
            .unwrap();
        assert_eq!(short.cell_visits, 1);
        assert_eq!(short.optical_depth, [0.5]);
        for (point, distance) in [(Point3::new(-1.0, 2.0, 0.5), f64::INFINITY), (origin, 0.0)] {
            let miss = medium
                .transport_observed::<1, true>(cx, point, direction, distance, &[550.0], 0)
                .unwrap();
            assert_eq!(miss.cell_visits, 0);
            assert_eq!(miss.optical_depth, [0.0]);
            assert_eq!(miss.transport.transmittance, [1.0]);
        }
        gate.request();
        assert_eq!(
            medium.transport_observed::<1, false>(cx, origin, direction, 0.0, &[550.0], 0,),
            Err(VolumeTransportError::Cancelled)
        );
    });
}

#[test]
fn g0_observed_finite_columns_keep_log_survival_and_refuse_range_loss() {
    with_cx(|cx, _| {
        let direction = Vec3::new(1.0, 0.0, 0.0);
        let origin = Point3::new(-1.0, 0.5, 0.5);
        for (density, coefficient, expected) in [
            (1.0, 1000.0, Some(1000.0)),
            (1e300, 1e300, None),
            (1e-300, 1e-300, None),
        ] {
            let values = [density];
            let medium = AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &values, [0.0; 3], [1.0; 3]),
                MassAbsorptionSpectrum::constant(coefficient).unwrap(),
                None,
                budget(),
            )
            .unwrap();
            let result = medium.transport_observed::<1, false>(
                cx,
                origin,
                direction,
                f64::INFINITY,
                &[550.0],
                1,
            );
            if let Some(depth) = expected {
                let result = result.unwrap();
                assert_eq!(result.optical_depth, [depth]);
                assert_eq!(result.transport.transmittance, [0.0]);
                assert_eq!(result.cell_visits, 1);
            } else {
                assert_eq!(result, Err(VolumeTransportError::Unrepresentable));
            }
        }
    });
}

#[test]
fn g0_bright_thermal_source_survives_underflowed_foreground_weight() {
    with_cx(|cx, _| {
        let temperatures = [1.0, 1e8];
        let wavelength = 550.0;
        let source = planck_radiance_w_m2_sr_nm(temperatures[1], wavelength).unwrap();
        assert_eq!(
            planck_radiance_w_m2_sr_nm(temperatures[0], wavelength).unwrap(),
            0.0
        );
        let origin = Point3::new(-1.0, 0.5, 0.5);
        let direction = Vec3::new(1.0, 0.0, 0.0);
        for foreground_depth in [750.0, 1000.0] {
            let density = [foreground_depth, 1.0];
            let medium = AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([2, 1, 1], &density, [0.0; 3], [1.0; 3]),
                MassAbsorptionSpectrum::constant(1.0).unwrap(),
                Some(&temperatures),
                budget(),
            )
            .unwrap();
            assert_eq!((-foreground_depth).exp(), 0.0);
            let public = medium
                .transport(cx, origin, direction, f64::INFINITY, &[wavelength])
                .unwrap();
            let result = medium.transport_observed::<1, true>(
                cx,
                origin,
                direction,
                f64::INFINITY,
                &[wavelength],
                2,
            );
            let expected = (source.ln() + (-(-1.0_f64).exp_m1()).ln() - foreground_depth).exp();
            if expected > 0.0 {
                let result = result.unwrap();
                let actual = result.transport.emission_w_m2_sr_nm[0];
                assert!(actual > 0.0);
                assert!((actual - expected).abs() <= 8.0 * f64::from_bits(1));
                assert_eq!(public, result.transport);
                assert_eq!(result.optical_depth, [foreground_depth + 1.0]);
                assert_eq!(result.cell_visits, 2);
            } else {
                assert_eq!(public.emission_w_m2_sr_nm, [0.0]);
                assert_eq!(result, Err(VolumeTransportError::Unrepresentable));
            }
        }
    });
}

#[test]
fn g0_subnormal_transfer_weights_do_not_quantize_finite_thermal_radiance() {
    with_cx(|cx, _| {
        let temperatures = [1.0, 1e100];
        let wavelength = 550.0;
        let source = planck_radiance_w_m2_sr_nm(temperatures[1], wavelength).unwrap();
        let origin = Point3::new(-1.0, 0.5, 0.5);
        let direction = Vec3::new(1.0, 0.0, 0.0);
        // The first survival is subnormal; the second survival is normal but
        // multiplying it by a thin cell's absorption underflows to zero.
        for (foreground_depth, emitting_depth) in [(740.0_f64, 1.0_f64), (700.0, 1e-20)] {
            let density = [foreground_depth, emitting_depth];
            let medium = AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([2, 1, 1], &density, [0.0; 3], [1.0; 3]),
                MassAbsorptionSpectrum::constant(1.0).unwrap(),
                Some(&temperatures),
                budget(),
            )
            .unwrap();
            let absorbed = -(-emitting_depth).exp_m1();
            assert!((-foreground_depth).exp() * absorbed < f64::MIN_POSITIVE);
            let expected = (source.ln() + absorbed.ln() - foreground_depth).exp();
            assert!(expected > f64::MIN_POSITIVE);
            let public = medium
                .transport(cx, origin, direction, f64::INFINITY, &[wavelength])
                .unwrap();
            let observed = medium
                .transport_observed::<1, true>(
                    cx,
                    origin,
                    direction,
                    f64::INFINITY,
                    &[wavelength],
                    2,
                )
                .unwrap();
            assert_eq!(public, observed.transport);
            assert!((public.emission_w_m2_sr_nm[0] / expected - 1.0).abs() < 2e-12);
        }
    });
}

#[test]
fn g0_tiny_cell_absorption_retains_source_scaled_emission() {
    with_cx(|cx, _| {
        let density = [1e-200];
        let temperatures = [1e100];
        let wavelength = 550.0;
        let source = planck_radiance_w_m2_sr_nm(temperatures[0], wavelength).unwrap();
        let origin = Point3::new(-1.0, 0.5, 0.5);
        let direction = Vec3::new(1.0, 0.0, 0.0);
        for coefficient in [1e-120, 1e-200] {
            let medium = AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [0.0; 3], [1.0; 3]),
                MassAbsorptionSpectrum::constant(coefficient).unwrap(),
                Some(&temperatures),
                budget(),
            )
            .unwrap();
            // For these depths, 1-exp(-tau) equals tau to far beyond binary64
            // accuracy. This independent multiplication order stays normal.
            let expected = (source * density[0]) * coefficient;
            assert!(expected > f64::MIN_POSITIVE);
            let public = medium
                .transport(cx, origin, direction, f64::INFINITY, &[wavelength])
                .unwrap();
            assert_eq!(public.transmittance, [1.0]);
            assert!((public.emission_w_m2_sr_nm[0] / expected - 1.0).abs() < 2e-12);
            let observed = medium.transport_observed::<1, true>(
                cx,
                origin,
                direction,
                f64::INFINITY,
                &[wavelength],
                1,
            );
            if density[0] * coefficient == 0.0 {
                assert_eq!(observed, Err(VolumeTransportError::Unrepresentable));
            } else {
                assert_eq!(public, observed.unwrap().transport);
            }
        }
    });
}

#[test]
fn g3_public_shadow_transport_retains_many_weak_layers_in_either_direction() {
    with_cx(|cx, _| {
        const THIN_CELLS: usize = 8192;
        let thin = 2.0_f64.powi(-60);
        let mut density = vec![thin; THIN_CELLS + 1];
        density[0] = 1.0;
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([density.len(), 1, 1], &density, [0.0; 3], [1.0; 3]),
            MassAbsorptionSpectrum::constant(1.0).unwrap(),
            None,
            AbsorbingGridBudget {
                max_grid_cells: density.len(),
                max_steps_per_segment: density.len() as u32,
            },
        )
        .unwrap();
        let depth = 1.0 + THIN_CELLS as f64 * thin;
        let expected = (-depth).exp();
        let mut forward = None;
        for (x, direction) in [(-1.0, 1.0), (density.len() as f64 + 1.0, -1.0)] {
            let origin = Point3::new(x, 0.5, 0.5);
            let direction = Vec3::new(direction, 0.0, 0.0);
            let public = medium
                .transmittance(cx, origin, direction, f64::INFINITY, &[550.0])
                .unwrap();
            let observed = medium
                .transport_observed::<1, false>(
                    cx,
                    origin,
                    direction,
                    f64::INFINITY,
                    &[550.0],
                    density.len() as u32,
                )
                .unwrap();
            assert_eq!(observed.optical_depth, [depth]);
            assert_eq!(public, observed.transport.transmittance);
            assert!((public[0] - expected).abs() <= 2.0 * f64::EPSILON * expected);
            if let Some(forward) = forward {
                assert_eq!(public, forward);
            } else {
                forward = Some(public);
            }
        }
    });
}

#[test]
fn g0_public_opaque_limit_does_not_poison_compensated_depth_or_emission() {
    with_cx(|cx, _| {
        let density = [1e300, 1.0];
        let temperatures = [2500.0, 1e8];
        let wavelength = 550.0;
        let source = planck_radiance_w_m2_sr_nm(temperatures[0], wavelength).unwrap();
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([2, 1, 1], &density, [0.0; 3], [1.0; 3]),
            MassAbsorptionSpectrum::constant(1e300).unwrap(),
            Some(&temperatures),
            budget(),
        )
        .unwrap();
        let origin = Point3::new(-1.0, 0.5, 0.5);
        let direction = Vec3::new(1.0, 0.0, 0.0);
        let public = medium
            .transport(cx, origin, direction, f64::INFINITY, &[wavelength])
            .unwrap();
        assert_eq!(public.transmittance, [0.0]);
        assert_eq!(public.emission_w_m2_sr_nm, [source]);
        assert_eq!(
            medium
                .transmittance(cx, origin, direction, f64::INFINITY, &[wavelength])
                .unwrap(),
            public.transmittance,
        );
        assert_eq!(
            medium.transport_observed::<1, true>(
                cx,
                origin,
                direction,
                f64::INFINITY,
                &[wavelength],
                2,
            ),
            Err(VolumeTransportError::Unrepresentable),
        );
    });
}
