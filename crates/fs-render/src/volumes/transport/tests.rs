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
