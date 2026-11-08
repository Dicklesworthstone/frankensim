//! Free-mode checks against analytic original-coordinate models and the
//! independent existing time propagator. No forced peak is used as a pole oracle.

use super::super::poles::{PoleOptions, StructuralInteriorPole, StructuralInteriorPoleError};
use super::super::time::StructuralInteriorTimeLimits;
use super::*;

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
