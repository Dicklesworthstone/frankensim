//! Analytic wave and work regressions for the original full pressure owner.

use super::*;
use crate::fixtures::masked_cube_grid;
use fs_rep_mesh::TetComplex;

const UNIFORM: Medium = Medium {
    density_kg_m3: 1.3,
    sound_speed_m_s: 2.0,
};

fn brick(cells: [usize; 3], lengths: [f64; 3]) -> (TetComplex, Vec<[f64; 3]>) {
    let (complex, mut positions) = masked_cube_grid(cells[0], cells[1], cells[2], &|_, _, _| true);
    for point in &mut positions {
        for axis in 0..3 {
            point[axis] *= lengths[axis] / cells[axis] as f64;
        }
    }
    (complex, positions)
}

fn plane(fe: &HelmholtzFe<'_>, axis: usize, coordinate: f64) -> Vec<[u32; 3]> {
    fe.boundary_faces()
        .into_iter()
        .filter(|face| {
            face.iter()
                .all(|&v| (fe.positions()[v as usize][axis] - coordinate).abs() < 1e-12)
        })
        .collect()
}

fn constant_velocity(faces: &[[u32; 3]], value: C64) -> Vec<VelocityTrace> {
    faces.iter().map(|&face| (face, [value; 3])).collect()
}

fn wave(phase: f64) -> C64 {
    C64::new(det::cos(phase), -det::sin(phase))
}

fn uniform<'f, 'c>(
    geometry: &'f HelmholtzFe<'c>,
    medium: Medium,
) -> HeterogeneousHelmholtzFe<'f, 'c> {
    HeterogeneousHelmholtzFe::new(
        geometry,
        geometry.order() + 2,
        &mut |_, _| Ok::<_, AcousticFeError>(medium),
        &mut || false,
    )
    .unwrap()
}

fn max_field_error(
    model: &HeterogeneousHelmholtzFe<'_, '_>,
    field: &HeterogeneousPressureField,
    exact: impl Fn([f64; 3]) -> C64,
) -> f64 {
    let geometry = model.geometry();
    let mut error: f64 = 0.0;
    for (i, &point) in geometry.positions().iter().enumerate() {
        error = error.max((field.pressure_pa[i] - exact(point)).abs());
    }
    for (element, tet) in geometry.tetrahedra().iter().enumerate() {
        let mut center = [0.0; 3];
        for &vertex in tet {
            for axis in 0..3 {
                center[axis] += 0.25 * geometry.positions()[vertex as usize][axis];
            }
        }
        let sampled = model
            .pressure_at(&field.pressure_pa, element, [0.25; 4])
            .unwrap();
        error = error.max((sampled - exact(center)).abs());
    }
    error
}

fn mean_pressure(geometry: &HelmholtzFe<'_>, pressure: &[C64], faces: &[[u32; 3]]) -> C64 {
    let (load, area) = geometry.boundary_load(faces).unwrap();
    pressure
        .iter()
        .zip(load)
        .fold(C64::ZERO, |sum, (&p, b)| sum + p.scale(b / area))
}

#[test]
fn uniform_rigid_modes_recover_the_original_homogeneous_owner() {
    let (complex, positions) = brick([2, 1, 1], [1.0, 0.35, 0.2]);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, UNIFORM).unwrap();
    let model = uniform(&geometry, UNIFORM);
    let options = SliceOptions {
        mass_diagonal_equilibration: true,
        ..SliceOptions::default()
    };
    let old = geometry.modes((0.5, 1.5), &options).unwrap();
    let actual = model.modes((0.5, 1.5), &options, &mut || false).unwrap();
    assert_eq!(actual.frequency_hz.len(), 1);
    assert_eq!(
        actual.below_high - actual.below_low,
        actual.frequency_hz.len()
    );
    assert_eq!(old.frequency_hz.len(), actual.frequency_hz.len());
    assert!((actual.frequency_hz[0] - 1.0).abs() < 0.01);
    for i in 0..old.frequency_hz.len() {
        assert!((actual.frequency_hz[i] / old.frequency_hz[i] - 1.0).abs() < 2e-8);
        let expected = old.lambda[i] * UNIFORM.sound_speed_m_s.powi(2);
        assert!((actual.omega_squared_rad2_s2[i] / expected - 1.0).abs() < 4e-8);
    }
}

#[test]
fn uniform_full_pressure_recovers_traveling_wave_work_and_energy() {
    let lengths = [1.0, 0.2, 0.2];
    let (complex, positions) = brick([4, 1, 1], lengths);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, UNIFORM).unwrap();
    let model = uniform(&geometry, UNIFORM);
    let inlet = plane(&geometry, 0, 0.0);
    let outlet = plane(&geometry, 0, 1.0);
    let impedance = UNIFORM.density_kg_m3 * UNIFORM.sound_speed_m_s;
    let omega = 1.6;
    let prepared = model
        .prepare(
            omega / core::f64::consts::TAU,
            &[BoundaryAdmittance {
                faces: outlet,
                admittance_m_pa_s: C64::from_re(1.0 / impedance),
            }],
            HarmonicLimits::default(),
            &mut || false,
        )
        .unwrap();
    let field = prepared
        .solve(&constant_velocity(&inlet, C64::ONE), &mut || false)
        .unwrap();
    let error = max_field_error(&model, &field, |point| {
        wave(omega * point[0] / UNIFORM.sound_speed_m_s).scale(impedance)
    }) / impedance;
    let area = lengths[1] * lengths[2];
    let expected_power = 0.5 * impedance * area;
    let expected_energy = 0.25 * UNIFORM.density_kg_m3 * lengths.iter().product::<f64>();
    assert!(error < 2e-4, "traveling pressure error {error:e}");
    assert!((field.input_power_w / expected_power - 1.0).abs() < 4e-4);
    assert!(field.power_balance_residual_w().abs() < 1e-10 * expected_power);
    assert!((field.kinetic_energy_j / expected_energy - 1.0).abs() < 4e-4);
    assert!((field.compressional_energy_j / expected_energy - 1.0).abs() < 4e-4);
}

#[test]
fn a_density_jump_reflects_even_when_sound_speed_is_unchanged() {
    let left = Medium {
        density_kg_m3: 1.2,
        sound_speed_m_s: 2.0,
    };
    let right = Medium {
        density_kg_m3: 4.8,
        sound_speed_m_s: 2.0,
    };
    let (complex, positions) = brick([8, 1, 1], [1.0, 0.2, 0.2]);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, UNIFORM).unwrap();
    let model = HeterogeneousHelmholtzFe::new(
        &geometry,
        4,
        &mut |_, point| Ok::<_, AcousticFeError>(if point[0] < 0.5 { left } else { right }),
        &mut || false,
    )
    .unwrap();
    let inlet = plane(&geometry, 0, 0.0);
    let outlet = plane(&geometry, 0, 1.0);
    let z1 = left.density_kg_m3 * left.sound_speed_m_s;
    let z2 = right.density_kg_m3 * right.sound_speed_m_s;
    let reflection = (z2 - z1) / (z2 + z1);
    let omega = 2.0;
    let k = omega / left.sound_speed_m_s;
    let prepared = model
        .prepare(
            omega / core::f64::consts::TAU,
            &[
                BoundaryAdmittance {
                    faces: inlet.clone(),
                    admittance_m_pa_s: C64::from_re(1.0 / z1),
                },
                BoundaryAdmittance {
                    faces: outlet,
                    admittance_m_pa_s: C64::from_re(1.0 / z2),
                },
            ],
            HarmonicLimits::default(),
            &mut || false,
        )
        .unwrap();
    // Matched source boundary: incident pressure is exactly one pascal.
    let field = prepared
        .solve(
            &constant_velocity(&inlet, C64::from_re(2.0 / z1)),
            &mut || false,
        )
        .unwrap();
    let error = max_field_error(&model, &field, |point| {
        if point[0] <= 0.5 {
            wave(k * point[0]) + wave(k * (1.0 - point[0])).scale(reflection)
        } else {
            wave(k * point[0]).scale(1.0 + reflection)
        }
    });
    let recovered = (mean_pressure(&geometry, &field.pressure_pa, &inlet) - C64::ONE) * wave(-k);
    assert!((reflection - 0.6).abs() < 1e-15);
    assert!(error < 5e-5, "density-interface field error {error:e}");
    assert!((recovered - C64::from_re(reflection)).abs() < 3e-5);
    assert!(field.power_balance_residual_w().abs() < 1e-11);
}

#[test]
fn continuous_sound_speed_recovers_exact_travel_time_on_a_mesh_ladder() {
    // rho = exp(-a x), c = c0 exp(a x) has CONSTANT characteristic impedance.
    // Its exact wave is exp(-i omega integral_0^x dx/c), with no reflection.
    let alpha = 1.2;
    let c0 = 1.5;
    let omega = 2.0;
    let exact = |point: [f64; 3]| {
        let delay = -det::expm1(-alpha * point[0]) / (alpha * c0);
        wave(omega * delay)
    };
    let mut errors = Vec::new();
    for nx in [2, 4, 8] {
        let (complex, positions) = brick([nx, 1, 1], [1.0, 0.2, 0.2]);
        let geometry = HelmholtzFe::new(&complex, &positions, 2, UNIFORM).unwrap();
        let model = HeterogeneousHelmholtzFe::new(
            &geometry,
            4,
            &mut |_, point| {
                Ok::<_, AcousticFeError>(Medium {
                    density_kg_m3: det::exp(-alpha * point[0]),
                    sound_speed_m_s: c0 * det::exp(alpha * point[0]),
                })
            },
            &mut || false,
        )
        .unwrap();
        let inlet = plane(&geometry, 0, 0.0);
        let prepared = model
            .prepare(
                omega / core::f64::consts::TAU,
                &[BoundaryAdmittance {
                    faces: plane(&geometry, 0, 1.0),
                    admittance_m_pa_s: C64::from_re(1.0 / c0),
                }],
                HarmonicLimits::default(),
                &mut || false,
            )
            .unwrap();
        let field = prepared
            .solve(
                &constant_velocity(&inlet, C64::from_re(1.0 / c0)),
                &mut || false,
            )
            .unwrap();
        let error = max_field_error(&model, &field, exact);
        println!("heterogeneous continuous sound-speed nx={nx} error={error:e}");
        errors.push(error);
        assert!(field.power_balance_residual_w().abs() < 1e-11);
    }
    assert!((wave(omega / c0) - exact([1.0, 0.0, 0.0])).abs() > 0.25);
    assert!(errors[2] < 5e-4, "{errors:?}");
    assert!(
        errors[0] > 10.0 * errors[2] && errors[1] > 2.5 * errors[2],
        "{errors:?}"
    );
}

fn forward(normal_k: C64, distance: f64) -> C64 {
    wave(normal_k.re * distance).scale(det::exp(normal_k.im * distance))
}

// Exact Snell interface at x=1/2, with incident unit pressure at the interface.
// The outgoing complex root has Im(k_n) <= 0 in exp(+i omega t).
fn interface_wave(point: [f64; 3], right: bool, tangential_k: f64) -> (C64, [C64; 3]) {
    let omega: f64 = 2.0;
    let kn1 = C64::from_re(det::sqrt(omega * omega - tangential_k * tangential_k));
    let q2 = 1.0 - tangential_k * tangential_k; // c2=2
    let kn2 = if q2 >= 0.0 {
        C64::from_re(det::sqrt(q2))
    } else {
        C64::new(0.0, -det::sqrt(-q2))
    };
    let y1 = kn1.scale(1.0 / omega); // rho1=1
    let y2 = kn2.scale(1.0 / (2.0 * omega)); // rho2=2
    let reflection = (y1 - y2) / (y1 + y2);
    let transmission = C64::ONE + reflection;
    let distance = point[0] - 0.5;
    let tangent = wave(tangential_k * point[1]);
    let (p, dx) = if right {
        let p = transmission * forward(kn2, distance);
        (p, C64::new(kn2.im, -kn2.re) * p)
    } else {
        let incident = forward(kn1, distance);
        let reflected = reflection * forward(kn1, -distance);
        (
            incident + reflected,
            C64::new(kn1.im, -kn1.re) * incident + C64::new(-kn1.im, kn1.re) * reflected,
        )
    };
    let p = p * tangent;
    (
        p,
        [dx * tangent, C64::new(0.0, -tangential_k) * p, C64::ZERO],
    )
}

#[test]
fn oblique_refraction_and_evanescent_pressure_converge_to_snell_fields() {
    for tangential_k in [0.6, 1.4] {
        let mut errors = Vec::new();
        let mut intensity_errors = Vec::new();
        for refinement in [1, 2, 4] {
            let (complex, positions) = brick([2 * refinement, refinement, 1], [1.0, 0.5, 0.2]);
            let geometry = HelmholtzFe::new(&complex, &positions, 2, UNIFORM).unwrap();
            let mut source = |_: usize, point: [f64; 3]| {
                Ok::<_, AcousticFeError>(if point[0] < 0.5 {
                    Medium {
                        density_kg_m3: 1.0,
                        sound_speed_m_s: 1.0,
                    }
                } else {
                    Medium {
                        density_kg_m3: 2.0,
                        sound_speed_m_s: 2.0,
                    }
                })
            };
            let model = HeterogeneousHelmholtzFe::new(
                &geometry,
                4,
                &mut source,
                &mut || false,
            )
            .unwrap();
            // Impose the analytic normal velocity on the outer box only.
            // The internal material interface is entirely solved by the FE weak form.
            let traces = geometry
                .boundary_faces()
                .into_iter()
                .map(|face| {
                    let (vertices, normal, _) = geometry.boundary_geometry(face).unwrap();
                    let right = vertices.iter().map(|p| p[0]).sum::<f64>() / 3.0 > 0.5;
                    let density = if right { 2.0 } else { 1.0 };
                    let values = vertices.map(|point| {
                        let (_, gradient) = interface_wave(point, right, tangential_k);
                        let derivative = gradient
                            .iter()
                            .zip(normal)
                            .fold(C64::ZERO, |sum, (&g, n)| sum + g.scale(n));
                        // Inward velocity = -i (outward pressure derivative)/(omega rho).
                        C64::new(derivative.im, -derivative.re).scale(1.0 / (2.0 * density))
                    });
                    (face, values)
                })
                .collect::<Vec<_>>();
            let prepared = model
                .prepare(
                    2.0 / core::f64::consts::TAU,
                    &[],
                    HarmonicLimits::default(),
                    &mut || false,
                )
                .unwrap();
            let field = prepared.solve(&traces, &mut || false).unwrap();
            let error = max_field_error(&model, &field, |point| {
                interface_wave(point, point[0] > 0.5, tangential_k).0
            });
            let normal_k = if tangential_k < 1.0 {
                C64::from_re(det::sqrt(1.0 - tangential_k * tangential_k))
            } else {
                C64::new(0.0, -det::sqrt(tangential_k * tangential_k - 1.0))
            };
            let mut intensity_error: f64 = 0.0;
            for (element, tet) in geometry.tetrahedra().iter().enumerate() {
                let x = tet.iter().map(|&v| positions[v as usize][0]).sum::<f64>() / 4.0;
                if x <= 0.5 {
                    continue;
                }
                let probe = model.probe(element, [0.25; 4], &mut source).unwrap();
                let sample = probe.harmonic(&field.pressure_pa, prepared.frequency_hz()).unwrap();
                let exact_pressure = interface_wave(probe.position_m(), true, tangential_k).0;
                // Transmitted complex intensity is |p|² conj(k)/(2 omega rho).
                // Its evanescent normal component is reactive, while tangential
                // active flux survives. Zero components use the full-vector scale.
                let base = exact_pressure.norm_sq() / 8.0;
                let active = [normal_k.re * base, tangential_k * base, 0.0];
                let reactive = [-normal_k.im * base, 0.0, 0.0];
                let defect: f64 = (0..3).map(|axis| {
                    (sample.active_intensity_w_m2[axis] - active[axis]).powi(2)
                        + (sample.reactive_intensity_w_m2[axis] - reactive[axis]).powi(2)
                }).sum();
                let scale = base * det::sqrt(normal_k.norm_sq() + tangential_k * tangential_k);
                intensity_error = intensity_error.max(det::sqrt(defect) / scale);
                if refinement == 4 {
                    if tangential_k < 1.0 {
                        let [ix, iy, _] = sample.active_intensity_w_m2;
                        let along = normal_k.re * ix + tangential_k * iy;
                        let across = normal_k.re * iy - tangential_k * ix;
                        assert!(along > 0.0 && across.abs() < 0.12 * along,
                            "transported power must follow the refracted direction");
                    } else {
                        assert!(sample.reactive_intensity_w_m2[0] > 0.0);
                        assert!(sample.active_intensity_w_m2[1] > 0.0);
                        assert!(sample.active_intensity_w_m2[0].abs() < 0.2 * reactive[0],
                            "evanescent normal flow must approach zero active power");
                    }
                }
            }
            println!(
                "heterogeneous Snell kt={tangential_k} refinement={refinement} pressure={error:e} intensity={intensity_error:e}"
            );
            errors.push(error);
            intensity_errors.push(intensity_error);
            assert_eq!(field.boundary_dissipation_w, 0.0);
            assert!(field.input_power_w.abs() < 1e-10);
            if tangential_k > 1.0 && refinement == 4 {
                let vertex = |x: f64| {
                    positions
                        .iter()
                        .position(|p| (p[0] - x).abs() < 1e-12 && p[1] == 0.0 && p[2] == 0.0)
                        .unwrap()
                };
                let ratio =
                    field.pressure_pa[vertex(1.0)].abs() / field.pressure_pa[vertex(0.5)].abs();
                let expected = det::exp(-0.5 * det::sqrt(tangential_k * tangential_k - 1.0));
                assert!(
                    (ratio - expected).abs() < 0.02,
                    "evanescent ratio {ratio} vs {expected}"
                );
            }
        }
        assert!(errors[2] < 0.015, "kt={tangential_k}: {errors:?}");
        assert!(
            errors[0] > 5.0 * errors[2] && errors[1] > 2.5 * errors[2],
            "kt={tangential_k}: {errors:?}"
        );
        assert!(intensity_errors[2] < 0.1, "kt={tangential_k}: {intensity_errors:?}");
        assert!(intensity_errors[2] < intensity_errors[1] && intensity_errors[1] < intensity_errors[0],
            "gradient/intensity convergence is independent of pressure: kt={tangential_k}: {intensity_errors:?}");
    }
}

fn real_form(matrix: &Csr, field: &[f64]) -> f64 {
    let mut applied = vec![0.0; field.len()];
    matrix.spmv(field, &mut applied);
    field.iter().zip(applied).map(|(p, q)| p * q).sum()
}

#[test]
fn continuous_coefficients_have_a_live_independent_quadrature_axis() {
    let (complex, positions) = brick([1, 1, 1], [1.0; 3]);
    let geometry = HelmholtzFe::new(&complex, &positions, 1, UNIFORM).unwrap();
    let x = positions.iter().map(|p| p[0]).collect::<Vec<_>>();
    let one = vec![1.0; geometry.ndof()];
    let exact_k = -det::expm1(-5.0) / 5.0;
    let exact_m = -det::expm1(-6.6) / 6.6;
    let mut errors = Vec::new();
    for quadrature in [3, 5, 8] {
        let mut calls = 0;
        let model = HeterogeneousHelmholtzFe::new(
            &geometry,
            quadrature,
            &mut |_, point| {
                calls += 1;
                Ok::<_, AcousticFeError>(Medium {
                    density_kg_m3: det::exp(5.0 * point[0]),
                    sound_speed_m_s: det::exp(0.8 * point[0]),
                })
            },
            &mut || false,
        )
        .unwrap();
        assert_eq!(calls, complex.tets.len() * quadrature.pow(3));
        let error = (real_form(model.stiffness(), &x) / exact_k - 1.0)
            .abs()
            .max((real_form(model.mass(), &one) / exact_m - 1.0).abs());
        errors.push(error);
    }
    assert!(errors[0] > 1e-8 && errors[2] < 3e-12, "{errors:?}");
    assert!(errors[1] < 0.01 * errors[0], "{errors:?}");
}

#[test]
fn prepared_spatial_drives_are_reciprocal_and_vertex_order_independent() {
    let (complex, positions) = brick([2, 1, 1], [1.0, 0.3, 0.2]);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, UNIFORM).unwrap();
    let model = HeterogeneousHelmholtzFe::new(
        &geometry,
        5,
        &mut |_, point| {
            Ok::<_, AcousticFeError>(Medium {
                density_kg_m3: 1.0 + point[0],
                sound_speed_m_s: 1.5 + 0.5 * point[0],
            })
        },
        &mut || false,
    )
    .unwrap();
    let left = plane(&geometry, 0, 0.0);
    let right = plane(&geometry, 0, 1.0);
    let prepared = model
        .prepare(
            0.3,
            &[
                BoundaryAdmittance {
                    faces: left.clone(),
                    admittance_m_pa_s: C64::new(0.7, 0.2),
                },
                BoundaryAdmittance {
                    faces: right.clone(),
                    admittance_m_pa_s: C64::new(0.4, -0.1),
                },
            ],
            HarmonicLimits::default(),
            &mut || false,
        )
        .unwrap();
    let a = (left[0], [1.0, 0.2, 1.8]);
    let b = (right[0], [0.5, 1.1, 0.7]);
    let drive = |trace: ([u32; 3], [f64; 3])| [(trace.0, trace.1.map(C64::from_re))];
    let pa = prepared.solve(&drive(a), &mut || false).unwrap();
    let pb = prepared.solve(&drive(b), &mut || false).unwrap();
    let la = geometry.boundary_linear_load(&[a]).unwrap();
    let lb = geometry.boundary_linear_load(&[b]).unwrap();
    let response = |pressure: &[C64], load: &[f64]| {
        pressure
            .iter()
            .zip(load)
            .fold(C64::ZERO, |sum, (&p, &w)| sum + p.scale(w))
    };
    let ab = response(&pa.pressure_pa, &lb);
    let ba = response(&pb.pressure_pa, &la);
    assert!((ab - ba).abs() < 3e-10 * ab.abs().max(ba.abs()));
    assert!(pa.power_balance_residual_w().abs() < 1e-11);
    assert!(pb.power_balance_residual_w().abs() < 1e-11);
    let reversed = ([a.0[2], a.0[0], a.0[1]], [a.1[2], a.1[0], a.1[1]]);
    let pr = prepared.solve(&drive(reversed), &mut || false).unwrap();
    assert!(
        pa.pressure_pa
            .iter()
            .zip(&pr.pressure_pa)
            .all(|(&x, &y)| (x - y).abs() < 1e-10)
    );
    let zero = prepared.solve(&[], &mut || false).unwrap();
    assert!(zero.pressure_pa.iter().all(|&p| p == C64::ZERO));
    assert_eq!(zero.input_power_w, 0.0);
    assert_eq!(zero.kinetic_energy_j, 0.0);
    assert_eq!(zero.compressional_energy_j, 0.0);
}

#[test]
fn original_source_errors_and_physical_preparation_limits_are_preserved() {
    #[derive(Debug, PartialEq)]
    enum SourceError {
        Original,
        Interior(AcousticFeError),
    }
    impl From<AcousticFeError> for SourceError {
        fn from(value: AcousticFeError) -> Self {
            Self::Interior(value)
        }
    }
    let (complex, positions) = brick([1, 1, 1], [1.0; 3]);
    let geometry = HelmholtzFe::new(&complex, &positions, 1, UNIFORM).unwrap();
    let result = HeterogeneousHelmholtzFe::new(
        &geometry,
        3,
        &mut |_, _| Err::<Medium, _>(SourceError::Original),
        &mut || false,
    );
    assert!(matches!(result, Err(SourceError::Original)));
    let mut calls = 0;
    let cancelled = HeterogeneousHelmholtzFe::new(
        &geometry,
        3,
        &mut |_, _| {
            calls += 1;
            Ok::<_, AcousticFeError>(UNIFORM)
        },
        &mut || true,
    );
    assert!(cancelled.is_err());
    assert_eq!(calls, 0);
    let model = uniform(&geometry, UNIFORM);
    assert!(matches!(
        model.prepare(
            1.0,
            &[],
            HarmonicLimits {
                max_dofs: 1,
                ..HarmonicLimits::default()
            },
            &mut || false,
        ),
        Err(HeterogeneousError::Limit { .. })
    ));
    assert!(matches!(
        model.prepare(1.0, &[], HarmonicLimits::default(), &mut || true,),
        Err(HeterogeneousError::Cancelled)
    ));
    let active = [BoundaryAdmittance {
        faces: plane(&geometry, 0, 1.0),
        admittance_m_pa_s: C64::from_re(-1.0),
    }];
    assert!(
        model
            .prepare(1.0, &active, HarmonicLimits::default(), &mut || false)
            .is_err()
    );
}
