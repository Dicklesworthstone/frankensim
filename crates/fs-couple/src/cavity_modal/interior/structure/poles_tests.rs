//! Free-mode checks against analytic original-coordinate models and the
//! independent existing time propagator. No forced peak is used as a pole oracle.

use super::super::poles::parameter::{
    PoleResponseOptions, PoleSensitivityOptions, StructuralInteriorPoleDerivative,
    StructuralInteriorPoleDirection,
};
use super::super::poles::{PoleOptions, StructuralInteriorPole, StructuralInteriorPoleError};
use super::super::time::StructuralInteriorTimeLimits;
use super::*;

#[path = "poles/relaxation/residue_tests.rs"]
mod relaxation_residue_tests;
#[path = "poles/relaxation_tests.rs"]
mod relaxation_tests;

fn options() -> PoleOptions {
    PoleOptions {
        spectral_scale: 4.0,
        step_tolerance: 2e-13,
        residual_tolerance: 2e-12,
        max_iterations: 40,
        max_dimension: 32,
        max_dense_bytes: 16 * 1024 * 1024,
        max_work: 100_000_000,
    }
}

fn sensitivity_options() -> PoleSensitivityOptions {
    PoleSensitivityOptions {
        numerical: PoleResponseOptions {
            spectral_scale: 4.0,
            pole_residual_tolerance: 3e-12,
            pole_step_tolerance: 3e-12,
            linear_residual_tolerance: 3e-11,
            max_dimension: 32,
            max_vectors: 8,
            max_dense_bytes: 16 * 1024 * 1024,
            max_work: 100_000_000,
        },
        max_observations: 0,
        max_physical_bytes: 16 * 1024 * 1024,
        max_physical_work: 100_000_000,
    }
}

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual:?} != {expected:?}"
    );
}

fn original_rows_and_work(model: &StructuralInterior<'_, '_, '_>, pole: &StructuralInteriorPole) {
    let ns = model.structural_dofs;
    let n = ns + model.pressure_dofs();
    let q: Vec<_> = pole
        .structural_displacement
        .iter()
        .chain(&pole.pressure_potential_pa_s)
        .copied()
        .collect();
    let v: Vec<_> = pole
        .structural_velocity
        .iter()
        .chain(&pole.pressure_pa)
        .copied()
        .collect();
    let s = pole.laplace_pole_per_s;
    let mut mass_norm = C64::ZERO;
    for i in 0..n {
        close(v[i], s * q[i], 3e-12);
        let mut row = C64::ZERO;
        let mut scale = 0.0;
        for j in 0..n {
            let m = model.mass[i * n + j];
            let c = model.damping[i * n + j];
            let k = model.stiffness[i * n + j];
            let term = (s * s).scale(m) + s.scale(c) + C64::from_re(k);
            row = row + term * q[j];
            scale += (s.norm_sq() * m.abs() + s.abs() * c.abs() + k.abs()) * q[j].abs();
            mass_norm = mass_norm + q[i].conj() * q[j].scale(m);
        }
        assert!(
            row.abs() <= 2e-10 * scale.max(1.0),
            "original row {i}: {row:?}"
        );
    }
    close(mass_norm, C64::ONE, 3e-12);
    let scale = (s.abs() * pole.energy.total_j()).max(1.0);
    assert!(pole.decay_power_residual_w.abs() < 3e-10 * scale);
    let structure = 2.0 * s.re * pole.energy.structural_j() + pole.structural_dissipation_w
        - pole.fluid_to_structure_power_w;
    let acoustic = 2.0 * s.re * pole.energy.acoustic_j()
        + pole.wall_dissipation_w
        + pole.fluid_to_structure_power_w;
    assert!(
        structure.abs() < 3e-10 * scale,
        "structural envelope: {structure}"
    );
    assert!(
        acoustic.abs() < 3e-10 * scale,
        "acoustic envelope: {acoustic}"
    );
    assert!(pole.relative_residual < options().residual_tolerance);
}

fn dilation(geometry: &HelmholtzFe<'_>) -> Vec<NormalDisplacementTrace> {
    // P1 affine dilation u(x)=q*(x-centroid) has B_i=3V/4=1/8.
    // B is proportional to M_acoustic * 1, so constant compression is an
    // EXACT invariant subspace of this original FE model, not a low-ka fit.
    geometry
        .boundary_faces()
        .into_iter()
        .map(|face| {
            let (vertices, normal, _) = geometry.boundary_geometry(face).unwrap();
            (
                face,
                vertices.map(|point| (0..3).map(|axis| (point[axis] - 0.25) * normal[axis]).sum()),
            )
        })
        .collect()
}

#[test]
fn scalar_damped_overdamped_and_growing_modes_retain_the_original_poles() {
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 1.2,
        sound_speed_m_s: 3.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 1, medium).unwrap();
    let interior = uniform(&geometry, medium);
    let mass = dense(&[2.0], 1);
    let stiffness = dense(&[18.0], 1);
    for drag in [0.4_f64, 16.0, -0.4] {
        let damping = dense(&[drag], 1);
        let discriminant = drag * drag - 144.0;
        let root = if discriminant >= 0.0 {
            C64::from_re(discriminant.sqrt() / 4.0)
        } else {
            C64::new(0.0, (-discriminant).sqrt() / 4.0)
        };
        with_cx(|cx| {
            let model = StructuralInterior::new(
                &interior,
                StructuralMatrices {
                    mass: &mass,
                    stiffness: &stiffness,
                    damping: Some(&damping),
                },
                &[Vec::new()],
                &[],
                controls(),
                cx,
            )
            .unwrap();
            let prepared = model.prepare_resonance_modes(options(), cx).unwrap();
            for sign in [-1.0, 1.0] {
                let expected = C64::from_re(-drag / 4.0) + root.scale(sign);
                let pole = prepared
                    .refine(expected.scale(1.05), &[C64::ONE], &[C64::ZERO; 4], cx)
                    .unwrap();
                close(pole.laplace_pole_per_s, expected, 2e-12);
                close(
                    pole.structural_displacement[0],
                    C64::from_re(1.0 / 2.0_f64.sqrt()),
                    2e-12,
                );
                assert!(pole.pressure_pa.iter().all(|&p| p == C64::ZERO));
                assert!(pole.pressure_potential_pa_s.iter().all(|&p| p == C64::ZERO));
                assert_eq!(pole.wall_dissipation_w, 0.0);
                assert_eq!(pole.fluid_to_structure_power_w, 0.0);
                if drag < 0.0 {
                    assert!(pole.amplitude_decay_rate_per_s < 0.0);
                    assert!(pole.structural_dissipation_w < 0.0);
                }
                if discriminant > 0.0 {
                    assert_eq!(pole.oscillation_frequency_rad_s, 0.0);
                }
                original_rows_and_work(&model, &pole);
            }
        });
    }
}

#[test]
fn gyroscopic_frequency_split_preserves_two_original_coordinates_without_heat() {
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 1.2,
        sound_speed_m_s: 3.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 1, medium).unwrap();
    let interior = uniform(&geometry, medium);
    let (m, k, g) = (2.0_f64, 18.0_f64, 1.3_f64);
    let mass = dense(&[m, 0.0, 0.0, m], 2);
    let stiffness = dense(&[k, 0.0, 0.0, k], 2);
    let damping = dense(&[0.0, g, -g, 0.0], 2);
    with_cx(|cx| {
        let model = StructuralInterior::new(
            &interior,
            StructuralMatrices {
                mass: &mass,
                stiffness: &stiffness,
                damping: Some(&damping),
            },
            &[Vec::new(), Vec::new()],
            &[],
            controls(),
            cx,
        )
        .unwrap();
        let prepared = model.prepare_resonance_modes(options(), cx).unwrap();
        for sign in [-1.0, 1.0] {
            let omega = ((g * g + 4.0 * m * k).sqrt() + sign * g) / (2.0 * m);
            let expected = C64::new(0.0, -omega);
            let pole = prepared
                .refine(
                    expected.scale(1.03),
                    &[C64::ONE, C64::new(0.0, sign)],
                    &[C64::ZERO; 4],
                    cx,
                )
                .unwrap();
            close(pole.laplace_pole_per_s, expected, 3e-12);
            close(
                pole.structural_displacement[1] / pole.structural_displacement[0],
                C64::new(0.0, sign),
                3e-12,
            );
            assert_eq!(pole.structural_dissipation_w, 0.0);
            assert_eq!(pole.wall_dissipation_w, 0.0);
            original_rows_and_work(&model, &pole);
        }
    });
}

#[test]
fn closed_cavity_free_mode_matches_exact_loading_and_the_original_time_evolution() {
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 1.2,
        sound_speed_m_s: 3.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 1, medium).unwrap();
    let interior = uniform(&geometry, medium);
    let (m, k, c) = (2.0_f64, 3.0_f64, 0.4_f64);
    let mass = dense(&[m], 1);
    let stiffness = dense(&[k], 1);
    let damping = dense(&[c], 1);
    let bulk = medium.density_kg_m3 * medium.sound_speed_m_s.powi(2);
    let loaded_k = k + 9.0 * (1.0 / 6.0) * bulk;
    let decay = c / (2.0 * m);
    let expected = C64::new(-decay, -(loaded_k / m - decay * decay).sqrt());
    with_cx(|cx| {
        let model = StructuralInterior::new(
            &interior,
            StructuralMatrices {
                mass: &mass,
                stiffness: &stiffness,
                damping: Some(&damping),
            },
            &[dilation(&geometry)],
            &[],
            controls(),
            cx,
        )
        .unwrap();
        let prepared = model.prepare_resonance_modes(options(), cx).unwrap();
        let seed = C64::from_re(-3.0 * bulk) / expected;
        let pole = prepared
            .refine(expected.scale(1.04), &[C64::ONE], &[seed; 4], cx)
            .unwrap();
        close(pole.laplace_pole_per_s, expected, 5e-12);
        assert!(pole.oscillation_frequency_rad_s > (k / m).sqrt());
        for a in 0..4 {
            close(
                pole.pressure_pa[a],
                pole.structural_displacement[0].scale(-3.0 * bulk),
                5e-12,
            );
            close(
                pole.displaced_flow_m3_s[a],
                pole.structural_velocity[0].scale(0.125),
                5e-12,
            );
            let wrong = C64::new(0.0, expected.im) * pole.pressure_potential_pa_s[a];
            assert!(
                (wrong - pole.pressure_pa[a]).abs() > 0.01,
                "dropping the real decay contribution must change the physical pressure"
            );
        }
        close(pole.fluid_force[0], pole.pressure_pa[0].scale(0.5), 5e-12);
        original_rows_and_work(&model, &pole);

        // The exact held-input owner independently propagates REAL original
        // coordinates; a free complex mode must reproduce the same trajectory.
        let duration = 0.17;
        let exact = model
            .prepare_exact_time(duration, StructuralInteriorTimeLimits::default(), cx)
            .unwrap();
        let real = |values: &[C64]| values.iter().map(|z| z.re).collect::<Vec<_>>();
        let initial = exact
            .state(
                0.0,
                &real(&pole.structural_displacement),
                &real(&pole.structural_velocity),
                &real(&pole.pressure_potential_pa_s),
                &real(&pole.pressure_pa),
                cx,
            )
            .unwrap();
        let step = exact.step(&initial, &[0.0], &[0.0; 4], cx).unwrap();
        let phase = expected.im * duration;
        let advance = C64::new(phase.cos(), phase.sin()).scale((expected.re * duration).exp());
        for (actual, shape) in [
            (
                step.state.structural_displacement(),
                pole.structural_displacement.as_slice(),
            ),
            (
                step.state.structural_velocity(),
                pole.structural_velocity.as_slice(),
            ),
            (
                step.state.pressure_potential_pa_s(),
                pole.pressure_potential_pa_s.as_slice(),
            ),
            (step.state.pressure_pa(), pole.pressure_pa.as_slice()),
        ] {
            for (&x, &phi) in actual.iter().zip(shape) {
                close(C64::from_re(x), C64::from_re((advance * phi).re), 2e-10);
            }
        }
        assert!(step.frame.work.balance_residual_j().abs() < 2e-10);
    });
}

#[test]
fn original_spatial_wall_loss_survives_refinement_and_budget_or_seed_refusal() {
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 1.2,
        sound_speed_m_s: 3.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 1, medium).unwrap();
    let interior = uniform(&geometry, medium);
    let mass = dense(&[2.0], 1);
    let stiffness = dense(&[3.0], 1);
    let damping = dense(&[0.4], 1);
    with_cx(|cx| {
        let model = StructuralInterior::new(
            &interior,
            StructuralMatrices {
                mass: &mass,
                stiffness: &stiffness,
                damping: Some(&damping),
            },
            &[dilation(&geometry)],
            &[ResistiveBoundary {
                faces: vec![[0, 1, 2]],
                admittance_m_pa_s: 0.002,
            }],
            controls(),
            cx,
        )
        .unwrap();
        for constrained in [
            PoleOptions {
                max_dense_bytes: 0,
                ..options()
            },
            PoleOptions {
                max_work: 0,
                ..options()
            },
            PoleOptions {
                max_dimension: 4,
                ..options()
            },
        ] {
            assert!(model.prepare_resonance_modes(constrained, cx).is_err());
        }
        let prepared = model.prepare_resonance_modes(options(), cx).unwrap();
        let seed = C64::new(-0.1, -3.1);
        let potential_seed = [C64::from_re(-32.4) / seed; 4];
        let pole = prepared
            .refine(seed, &[C64::ONE], &potential_seed, cx)
            .unwrap();
        assert!(pole.wall_dissipation_w > 0.0);
        assert!(pole.structural_dissipation_w > 0.0);
        assert!(pole.amplitude_decay_rate_per_s > 0.1);
        original_rows_and_work(&model, &pole);
        // This finite wall breaks uniform compression; its actual spatial
        // eigenvector and pressure/flow phase remain in the returned mode.
        assert!(
            pole.pressure_pa
                .iter()
                .skip(1)
                .any(|p| (*p - pole.pressure_pa[0]).abs() > 1e-5)
        );
        assert!(matches!(
            prepared.refine(seed, &[], &potential_seed, cx),
            Err(StructuralInteriorPoleError::Physical(
                StructuralInteriorError::Invalid(_)
            ))
        ));
        let repeated = prepared
            .refine(seed, &[C64::ONE], &potential_seed, cx)
            .unwrap();
        assert_eq!(pole.laplace_pole_per_s, repeated.laplace_pole_per_s);
        assert_eq!(pole.pressure_pa, repeated.pressure_pa);
        assert_eq!(
            pole.structural_displacement,
            repeated.structural_displacement
        );

        let gate = CancelGate::new();
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cancelled_cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 81,
                    kernel_id: 19,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            gate.request();
            assert!(matches!(
                prepared.refine(seed, &[C64::ONE], &potential_seed, &cancelled_cx),
                Err(StructuralInteriorPoleError::Physical(
                    StructuralInteriorError::Harmonic(HarmonicResponseError::Cancelled)
                ))
            ));
            assert!(
                model
                    .prepare_resonance_modes(options(), &cancelled_cx)
                    .is_err()
            );
        });
        let after_cancel = prepared
            .refine(seed, &[C64::ONE], &potential_seed, cx)
            .unwrap();
        assert_eq!(pole.laplace_pole_per_s, after_cancel.laplace_pole_per_s);
        assert_eq!(pole.pressure_pa, after_cancel.pressure_pa);
    });
}

fn parameter_rows(
    model: &StructuralInterior<'_, '_, '_>,
    mode: &StructuralInteriorPole,
    direction: StructuralInteriorPoleDirection<'_>,
    derivative: &StructuralInteriorPoleDerivative,
) {
    let ns = model.structural_dofs;
    let na = model.pressure_dofs();
    let n = ns + na;
    let q: Vec<_> = mode
        .structural_displacement
        .iter()
        .chain(&mode.pressure_potential_pa_s)
        .copied()
        .collect();
    let dq: Vec<_> = derivative
        .structural_displacement_derivative
        .iter()
        .chain(&derivative.pressure_potential_derivative_pa_s)
        .copied()
        .collect();
    let dv: Vec<_> = derivative
        .structural_velocity_derivative
        .iter()
        .chain(&derivative.pressure_derivative_pa)
        .copied()
        .collect();
    let s = mode.laplace_pole_per_s;
    let ds = derivative.laplace_pole_derivative_per_s;
    let mut normalization = C64::ZERO;
    for i in 0..n {
        close(dv[i], ds * q[i] + s * dq[i], 2e-10);
        let mut residual = C64::ZERO;
        for j in 0..n {
            let (m, c, k) = (
                model.mass[i * n + j],
                model.damping[i * n + j],
                model.stiffness[i * n + j],
            );
            let (dm, dc, dk) = if i < ns && j < ns {
                (
                    direction.structural_mass[i * ns + j],
                    direction.structural_damping[i * ns + j],
                    direction.structural_stiffness[i * ns + j],
                )
            } else if i >= ns && j >= ns {
                let a = (i - ns) * na + j - ns;
                (
                    direction.acoustic_mass[a],
                    direction.wall_damping[a],
                    direction.acoustic_stiffness[a],
                )
            } else if i < ns {
                (0.0, -direction.coupling[(j - ns) * ns + i], 0.0)
            } else {
                (0.0, direction.coupling[(i - ns) * ns + j], 0.0)
            };
            let operator = (s * s).scale(m) + s.scale(c) + C64::from_re(k);
            let frequency = s.scale(2.0 * m) + C64::from_re(c);
            let parameter = (s * s).scale(dm) + s.scale(dc) + C64::from_re(dk);
            residual = residual + operator * dq[j] + (frequency * ds + parameter) * q[j];
            normalization = normalization
                + (q[i].conj() * dq[j]).scale(2.0 * m)
                + (q[i].conj() * q[j]).scale(dm);
        }
        assert!(
            residual.abs() < 5e-9,
            "differentiated original row {i}: {residual:?}"
        );
    }
    assert!(normalization.re.abs() < 3e-10);
    for i in 0..ns {
        let expected = (0..na).fold(C64::ZERO, |sum, a| {
            sum + derivative.pressure_derivative_pa[a].scale(model.coupling[a * ns + i])
                + mode.pressure_pa[a].scale(direction.coupling[a * ns + i])
        });
        close(derivative.fluid_force_derivative[i], expected, 2e-10);
    }
    for a in 0..na {
        let expected = (0..ns).fold(C64::ZERO, |sum, i| {
            sum + derivative.structural_velocity_derivative[i].scale(model.coupling[a * ns + i])
                + mode.structural_velocity[i].scale(direction.coupling[a * ns + i])
        });
        close(
            derivative.displaced_flow_derivative_m3_s[a],
            expected,
            2e-10,
        );
    }
    assert!(derivative.linear_residual < 3e-11);
    assert!(derivative.polynomial_relative_residual < 3e-11);
    assert!(derivative.normalization_relative_residual < 3e-11);
}

#[test]
fn affine_cavity_parameter_gradients_keep_mass_normalization_and_reciprocal_interface_terms() {
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 1.2,
        sound_speed_m_s: 3.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 1, medium).unwrap();
    let interior = uniform(&geometry, medium);
    let (m, c, k, volume) = (2.0_f64, 0.4_f64, 3.0_f64, 1.0 / 6.0);
    let bulk = medium.density_kg_m3 * medium.sound_speed_m_s.powi(2);
    let mass = dense(&[m], 1);
    let stiffness = dense(&[k], 1);
    let damping = dense(&[c], 1);
    let expected = C64::new(
        -c / (2.0 * m),
        -((k + 9.0 * volume * bulk) / m - (c / (2.0 * m)).powi(2)).sqrt(),
    );
    with_cx(|cx| {
        let model = StructuralInterior::new(
            &interior,
            StructuralMatrices {
                mass: &mass,
                stiffness: &stiffness,
                damping: Some(&damping),
            },
            &[dilation(&geometry)],
            &[],
            controls(),
            cx,
        )
        .unwrap();
        let modes = model.prepare_resonance_modes(options(), cx).unwrap();
        let seed = C64::from_re(-3.0 * bulk) / expected;
        let prepared = modes
            .refine_with_sensitivities(
                expected.scale(1.04),
                &[C64::ONE],
                &[seed; 4],
                sensitivity_options(),
                cx,
            )
            .unwrap();
        let mode = prepared.mode();
        close(mode.laplace_pole_per_s, expected, 5e-12);
        let acoustic_mass: Vec<_> = (0..4)
            .flat_map(|i| (0..4).map(move |j| (i, j)))
            .map(|(i, j)| interior.mass().get(i, j))
            .collect();
        let acoustic_stiffness: Vec<_> = (0..4)
            .flat_map(|i| (0..4).map(move |j| (i, j)))
            .map(|(i, j)| interior.stiffness().get(i, j))
            .collect();
        let d_compliance: Vec<_> = acoustic_mass.iter().map(|&x| -x / bulk).collect();
        let zero = StructuralInteriorPoleDirection {
            structural_mass: &[0.0],
            structural_damping: &[0.0],
            structural_stiffness: &[0.0],
            acoustic_mass: &[0.0; 16],
            acoustic_stiffness: &[0.0; 16],
            wall_damping: &[0.0; 16],
            coupling: &[0.0; 4],
        };
        let directions = [
            StructuralInteriorPoleDirection {
                structural_mass: &[1.0],
                ..zero
            },
            StructuralInteriorPoleDirection {
                structural_damping: &[1.0],
                ..zero
            },
            StructuralInteriorPoleDirection {
                structural_stiffness: &[1.0],
                ..zero
            },
            StructuralInteriorPoleDirection {
                acoustic_mass: &d_compliance,
                ..zero
            },
            StructuralInteriorPoleDirection {
                coupling: model.coupling(),
                ..zero
            },
            StructuralInteriorPoleDirection {
                structural_mass: &[m],
                structural_damping: &[c],
                structural_stiffness: &[k],
                acoustic_mass: &acoustic_mass,
                acoustic_stiffness: &acoustic_stiffness,
                coupling: model.coupling(),
                ..zero
            },
            zero,
        ];
        let derivatives = prepared.parameter_derivatives(&directions, cx).unwrap();
        let denominator = expected.scale(2.0 * m) + C64::from_re(c);
        let expected_derivatives = [
            -(expected * expected) / denominator,
            -expected / denominator,
            -C64::ONE / denominator,
            C64::from_re(-9.0 * volume) / denominator,
            C64::from_re(-18.0 * volume * bulk) / denominator,
            C64::ZERO,
            C64::ZERO,
        ];
        for (index, ((direction, derivative), expected_ds)) in directions
            .iter()
            .zip(&derivatives)
            .zip(expected_derivatives)
            .enumerate()
        {
            close(derivative.laplace_pole_derivative_per_s, expected_ds, 2e-10);
            assert_eq!(
                derivative.amplitude_decay_rate_derivative_per_s,
                -derivative.laplace_pole_derivative_per_s.re
            );
            parameter_rows(&model, mode, *direction, derivative);
            let u = mode.structural_displacement[0];
            let du = derivative.structural_displacement_derivative[0];
            for a in 0..4 {
                let ratio_slope =
                    (derivative.pressure_derivative_pa[a] * u - mode.pressure_pa[a] * du) / (u * u);
                let expected_ratio = if index == 3 {
                    -3.0
                } else if index == 4 {
                    -3.0 * bulk
                } else {
                    0.0
                };
                close(ratio_slope, C64::from_re(expected_ratio), 3e-9);
            }
            let phase = prepared.phase_reference_coordinate();
            let fixed = if phase == 0 {
                derivative.structural_displacement_derivative[0]
            } else {
                derivative.pressure_potential_derivative_pa_s[phase - 1]
            };
            assert_eq!(fixed.im, 0.0);
        }
        // Scaling all original coefficients leaves s unchanged but changes
        // the physical unit-mass shape. dB additionally changes force and flow.
        let scale = &derivatives[5];
        close(
            scale.structural_displacement_derivative[0],
            mode.structural_displacement[0].scale(-0.5),
            2e-10,
        );
        close(
            scale.fluid_force_derivative[0],
            mode.fluid_force[0].scale(0.5),
            2e-10,
        );
        for a in 0..4 {
            close(
                scale.pressure_derivative_pa[a],
                mode.pressure_pa[a].scale(-0.5),
                2e-10,
            );
            close(
                scale.displaced_flow_derivative_m3_s[a],
                mode.displaced_flow_m3_s[a].scale(0.5),
                2e-10,
            );
        }
        assert_eq!(derivatives[6].laplace_pole_derivative_per_s, C64::ZERO);

        // An independently rebuilt bulk-modulus pair approaches the analytic
        // slope; the derivative implementation itself performs no differencing.
        let mut previous_error = f64::INFINITY;
        for relative_step in [1e-2, 1e-3] {
            let mut pair = Vec::new();
            for sign in [-1.0, 1.0] {
                let changed = Medium {
                    density_kg_m3: medium.density_kg_m3,
                    sound_speed_m_s: (bulk * (1.0 + sign * relative_step) / medium.density_kg_m3)
                        .sqrt(),
                };
                let fluid = uniform(&geometry, changed);
                let perturbed = StructuralInterior::new(
                    &fluid,
                    StructuralMatrices {
                        mass: &mass,
                        stiffness: &stiffness,
                        damping: Some(&damping),
                    },
                    &[dilation(&geometry)],
                    &[],
                    controls(),
                    cx,
                )
                .unwrap();
                pair.push(
                    perturbed
                        .prepare_resonance_modes(options(), cx)
                        .unwrap()
                        .refine(
                            expected,
                            &mode.structural_displacement,
                            &mode.pressure_potential_pa_s,
                            cx,
                        )
                        .unwrap()
                        .laplace_pole_per_s,
                );
            }
            let fd = (pair[1] - pair[0]).scale(1.0 / (2.0 * bulk * relative_step));
            let error = (fd - derivatives[3].laplace_pole_derivative_per_s).abs();
            assert!(error < previous_error * 0.03);
            previous_error = error;
        }
        assert!(previous_error < 1e-7);

        let bad = StructuralInteriorPoleDirection {
            coupling: &[0.0; 3],
            ..zero
        };
        assert!(prepared.parameter_derivatives(&[bad], cx).is_err());
        assert!(prepared.parameter_derivatives(&[], cx).is_err());
        assert!(prepared.parameter_derivatives(&[zero; 9], cx).is_err());
        let mut nonsymmetric = [0.0; 16];
        nonsymmetric[1] = 1.0;
        assert!(
            prepared
                .parameter_derivatives(
                    &[StructuralInteriorPoleDirection {
                        acoustic_mass: &nonsymmetric,
                        ..zero
                    }],
                    cx
                )
                .is_err()
        );
        for limited in [
            PoleSensitivityOptions {
                max_physical_bytes: 0,
                ..sensitivity_options()
            },
            PoleSensitivityOptions {
                max_physical_work: 0,
                ..sensitivity_options()
            },
        ] {
            assert!(
                modes
                    .refine_with_sensitivities(expected, &[C64::ONE], &[seed; 4], limited, cx)
                    .is_err()
            );
        }
        let gate = CancelGate::new();
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cancelled = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 81,
                    kernel_id: 19,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            gate.request();
            assert!(
                prepared
                    .parameter_derivatives(&directions, &cancelled)
                    .is_err()
            );
            assert!(
                modes
                    .refine_with_sensitivities(
                        expected,
                        &[C64::ONE],
                        &[seed; 4],
                        sensitivity_options(),
                        &cancelled
                    )
                    .is_err()
            );
        });
        assert_eq!(
            prepared.parameter_derivatives(&directions, cx).unwrap(),
            derivatives
        );
    });
}

#[test]
fn original_gyroscopic_direction_has_the_analytic_split_frequency_tangent() {
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 1.2,
        sound_speed_m_s: 3.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 1, medium).unwrap();
    let interior = uniform(&geometry, medium);
    let (m, k, g) = (2.0_f64, 18.0_f64, 1.3_f64);
    let mass = dense(&[m, 0.0, 0.0, m], 2);
    let stiffness = dense(&[k, 0.0, 0.0, k], 2);
    let damping = dense(&[0.0, g, -g, 0.0], 2);
    with_cx(|cx| {
        let model = StructuralInterior::new(
            &interior,
            StructuralMatrices {
                mass: &mass,
                stiffness: &stiffness,
                damping: Some(&damping),
            },
            &[Vec::new(), Vec::new()],
            &[],
            controls(),
            cx,
        )
        .unwrap();
        let modes = model.prepare_resonance_modes(options(), cx).unwrap();
        let direction = StructuralInteriorPoleDirection {
            structural_mass: &[0.0; 4],
            structural_stiffness: &[0.0; 4],
            structural_damping: &[0.0, 1.0, -1.0, 0.0],
            acoustic_mass: &[0.0; 16],
            acoustic_stiffness: &[0.0; 16],
            wall_damping: &[0.0; 16],
            coupling: &[0.0; 8],
        };
        for sign in [-1.0, 1.0] {
            let discriminant = (g * g + 4.0 * m * k).sqrt();
            let omega = (discriminant + sign * g) / (2.0 * m);
            let prepared = modes
                .refine_with_sensitivities(
                    C64::new(0.0, -omega),
                    &[C64::ONE, C64::new(0.0, sign)],
                    &[C64::ZERO; 4],
                    sensitivity_options(),
                    cx,
                )
                .unwrap();
            let derivative = prepared
                .parameter_derivatives(&[direction], cx)
                .unwrap()
                .remove(0);
            close(
                derivative.laplace_pole_derivative_per_s,
                C64::new(0.0, -(g / discriminant + sign) / (2.0 * m)),
                3e-11,
            );
            for &shape in &derivative.structural_displacement_derivative {
                close(shape, C64::ZERO, 3e-11);
            }
            parameter_rows(&model, prepared.mode(), direction, &derivative);
            assert_eq!(prepared.mode().structural_dissipation_w, 0.0);
        }
    });
}
