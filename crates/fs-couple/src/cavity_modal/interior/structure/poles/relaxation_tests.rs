//! Original-field relaxation checks against independent polynomial and work laws.

use super::*;
use crate::cavity_modal::interior::structure::poles::relaxation::{
    MaxwellDampingBranch, StructuralInteriorRelaxationError, StructuralInteriorRelaxationPole,
};
use fs_modal::relaxation_modes::MaxwellModeError;

fn original_equations_and_envelope(
    model: &StructuralInterior<'_, '_, '_>,
    branches: &[MaxwellDampingBranch],
    pole: &StructuralInteriorRelaxationPole,
) {
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
            for term in [s * v[j].scale(m), v[j].scale(c), q[j].scale(k)] {
                row = row + term;
                scale += term.abs();
            }
            mass_norm = mass_norm + q[i].conj() * q[j].scale(m);
        }
        for (branch, memory) in branches.iter().zip(&pole.memory) {
            let term = memory.stress.scale(branch.strain[i]);
            row = row + term;
            scale += term.abs();
        }
        assert!(row.abs() <= 3e-10 * scale.max(1.0), "row {i}: {row:?}");
    }
    close(mass_norm, C64::ONE, 3e-12);
    let (mut stored, mut loss, mut structure_work, mut acoustic_work) = (0.0, 0.0, 0.0, 0.0);
    assert_eq!(branches.len(), pole.memory.len());
    for (branch, memory) in branches.iter().zip(&pole.memory) {
        let mut strain_rate = C64::ZERO;
        for (i, &velocity) in v.iter().enumerate() {
            let rate = velocity.scale(branch.strain[i]);
            strain_rate = strain_rate + rate;
            let work = (memory.stress.conj() * rate).re;
            if i < ns {
                structure_work += work;
            } else {
                acoustic_work += work;
            }
        }
        close(
            (C64::ONE + s.scale(branch.relaxation_time_s)) * memory.stress,
            strain_rate.scale(branch.stiffness * branch.relaxation_time_s),
            3e-11,
        );
        let expected_energy = memory.stress.norm_sq() / (2.0 * branch.stiffness);
        let expected_loss = memory.stress.norm_sq() / (branch.stiffness * branch.relaxation_time_s);
        assert!((memory.energy_j - expected_energy).abs() < 1e-13 * expected_energy.max(1.0));
        assert!((memory.dissipation_w - expected_loss).abs() < 1e-13 * expected_loss.max(1.0));
        let balance = 2.0 * s.re * memory.energy_j + memory.dissipation_w
            - (memory.stress.conj() * strain_rate).re;
        assert!(balance.abs() < 3e-10 * expected_loss.max(1.0));
        stored += expected_energy;
        loss += expected_loss;
    }
    assert!((pole.relaxation_energy_j - stored).abs() < 2e-13 * stored.max(1.0));
    assert!((pole.relaxation_dissipation_w - loss).abs() < 2e-13 * loss.max(1.0));
    let scale = (s.abs() * (pole.energy.total_j() + stored)).max(1.0);
    let structure =
        2.0 * s.re * pole.energy.structural_j() + pole.structural_dissipation_w + structure_work
            - pole.fluid_to_structure_power_w;
    let acoustic = 2.0 * s.re * pole.energy.acoustic_j()
        + pole.wall_dissipation_w
        + acoustic_work
        + pole.fluid_to_structure_power_w;
    assert!(
        structure.abs() < 3e-10 * scale,
        "structure balance: {structure}"
    );
    assert!(
        acoustic.abs() < 3e-10 * scale,
        "acoustic balance: {acoustic}"
    );
    assert!(pole.decay_power_residual_w.abs() < 3e-10 * scale);
    for a in 0..model.pressure_dofs() {
        let flow = (0..ns).fold(C64::ZERO, |sum, j| {
            sum + v[j].scale(model.coupling[a * ns + j])
        });
        close(pole.displaced_flow_m3_s[a], flow, 1e-13);
    }
    for j in 0..ns {
        let force = (0..model.pressure_dofs()).fold(C64::ZERO, |sum, a| {
            sum + pole.pressure_pa[a].scale(model.coupling[a * ns + j])
        });
        close(pole.fluid_force[j], force, 1e-13);
    }
    assert!(pole.mechanical_relative_residual <= options().residual_tolerance);
    assert!(pole.relaxation_relative_residual <= options().residual_tolerance);
}

#[test]
fn fluid_loaded_relaxation_has_three_exact_cubic_poles_and_independent_memory_heat() {
    // Complete original P1 tetrahedron. Affine dilation gives exact loading
    // 9 V bulk=.3, so (s²+1)(1+s)+(5/4)s =
    // (s+1/2)(s²+s/2+2), with one real and two conjugate physical roots.
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 0.2,
        sound_speed_m_s: 1.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 1, medium).unwrap();
    let interior = uniform(&geometry, medium);
    let mass = dense(&[1.0], 1);
    let stiffness = dense(&[0.7], 1);
    let branches = [MaxwellDampingBranch {
        strain: vec![1.0, 0.0, 0.0, 0.0, 0.0],
        stiffness: 1.25,
        relaxation_time_s: 1.0,
    }];
    with_cx(|cx| {
        let model = StructuralInterior::new(
            &interior,
            StructuralMatrices {
                mass: &mass,
                stiffness: &stiffness,
                damping: None,
            },
            &[dilation(&geometry)],
            &[],
            controls(),
            cx,
        )
        .unwrap();
        let prepared = model
            .prepare_relaxation_modes(&branches, options(), cx)
            .unwrap();
        for expected in [
            C64::from_re(-0.5),
            C64::new(-0.25, -31.0_f64.sqrt() / 4.0),
            C64::new(-0.25, 31.0_f64.sqrt() / 4.0),
        ] {
            let potential = C64::from_re(-0.6) / expected;
            let pole = prepared
                .refine(expected.scale(1.03), &[C64::ONE], &[potential; 4], cx)
                .unwrap();
            close(pole.laplace_pole_per_s, expected, 3e-11);
            let u = pole.structural_displacement[0];
            let amplitude = 1.0 / (1.0 + 0.3 / expected.norm_sq()).sqrt();
            close(u, C64::from_re(amplitude), 3e-11);
            for &pressure in &pole.pressure_pa {
                close(pressure, u.scale(-0.6), 3e-11);
            }
            close(
                pole.memory[0].stress / u,
                expected.scale(1.25) / (C64::ONE + expected),
                3e-11,
            );
            let stored_per_u2 = if expected.im == 0.0 { 0.625 } else { 0.5 };
            assert!((pole.relaxation_energy_j / u.norm_sq() - stored_per_u2).abs() < 3e-11);
            assert!((pole.relaxation_dissipation_w - 2.0 * pole.relaxation_energy_j).abs() < 1e-13);
            assert_eq!(pole.structural_dissipation_w, 0.0);
            assert_eq!(pole.wall_dissipation_w, 0.0);
            original_equations_and_envelope(&model, &branches, &pole);
        }
        assert!(matches!(
            prepared.refine(C64::from_re(-1.0), &[C64::ONE], &[C64::ONE; 4], cx),
            Err(StructuralInteriorRelaxationError::Mode(
                MaxwellModeError::MaterialPole { branch: 0, .. }
            ))
        ));
        assert!(
            prepared
                .refine(C64::from_re(-0.5), &[], &[C64::ONE; 4], cx)
                .is_err()
        );
        for limited in [
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
            assert!(
                model
                    .prepare_relaxation_modes(&branches, limited, cx)
                    .is_err()
            );
        }
        let bad = [MaxwellDampingBranch {
            strain: vec![1.0],
            ..branches[0].clone()
        }];
        assert!(model.prepare_relaxation_modes(&bad, options(), cx).is_err());
        let initial = prepared
            .refine(
                C64::from_re(-0.51),
                &[C64::ONE],
                &[C64::from_re(1.2); 4],
                cx,
            )
            .unwrap();
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
                    .refine(
                        C64::from_re(-0.51),
                        &[C64::ONE],
                        &[C64::from_re(1.2); 4],
                        &cancelled
                    )
                    .is_err()
            );
        });
        let repeated = prepared
            .refine(
                C64::from_re(-0.51),
                &[C64::ONE],
                &[C64::from_re(1.2); 4],
                cx,
            )
            .unwrap();
        assert_eq!(initial.laplace_pole_per_s, repeated.laplace_pole_per_s);
        assert_eq!(initial.pressure_pa, repeated.pressure_pa);
        assert_eq!(initial.memory, repeated.memory);
    });
}

#[test]
fn pressure_potential_memory_is_the_declared_causal_wall_admittance() {
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 0.2,
        sound_speed_m_s: 1.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 1, medium).unwrap();
    let interior = uniform(&geometry, medium);
    let mass = dense(&[1.0], 1);
    let stiffness = dense(&[4.0], 1);
    let branches = [MaxwellDampingBranch {
        strain: vec![0.0, 0.25, 0.25, 0.25, 0.25],
        stiffness: 5.0 / 6.0,
        relaxation_time_s: 1.0,
    }];
    with_cx(|cx| {
        let model = StructuralInterior::new(
            &interior,
            StructuralMatrices {
                mass: &mass,
                stiffness: &stiffness,
                damping: None,
            },
            &[Vec::new()],
            &[],
            controls(),
            cx,
        )
        .unwrap();
        let prepared = model
            .prepare_relaxation_modes(&branches, options(), cx)
            .unwrap();
        // Uniform compression: M_total=V/bulk=5/6, unit summed trace.
        // The pressure admittance is Y(s)=(5/6)/(1+s), hence s²+s+1=0.
        for sign in [-1.0, 1.0] {
            let expected = C64::new(-0.5, sign * 3.0_f64.sqrt() / 2.0);
            let pole = prepared
                .refine(expected.scale(1.02), &[C64::ZERO], &[C64::ONE; 4], cx)
                .unwrap();
            close(pole.laplace_pole_per_s, expected, 3e-11);
            close(pole.structural_displacement[0], C64::ZERO, 3e-12);
            for &pressure in &pole.pressure_pa {
                close(pressure, pole.pressure_pa[0], 3e-11);
            }
            close(
                pole.memory[0].stress / pole.pressure_pa[0],
                C64::from_re(5.0 / 6.0) / (C64::ONE + expected),
                3e-11,
            );
            assert_eq!(pole.wall_dissipation_w, 0.0);
            assert!(pole.relaxation_dissipation_w > 0.0);
            original_equations_and_envelope(&model, &branches, &pole);
        }
    });
}

#[test]
fn mutual_relaxation_channels_retain_gyroscopic_handedness_without_gyro_heat() {
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 0.2,
        sound_speed_m_s: 1.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 1, medium).unwrap();
    let interior = uniform(&geometry, medium);
    let g = 31.0_f64.sqrt() / 12.0;
    let mass = dense(&[1.0, 0.0, 0.0, 1.0], 2);
    let stiffness = dense(&[1.5, 0.0, 0.0, 1.5], 2);
    let damping = dense(&[0.0, g, -g, 0.0], 2);
    let branches: Vec<_> = [[0.6, 0.8], [-0.8, 0.6]]
        .into_iter()
        .map(|b| MaxwellDampingBranch {
            strain: vec![b[0], b[1], 0.0, 0.0, 0.0, 0.0],
            stiffness: 35.0 / 24.0,
            relaxation_time_s: 1.0,
        })
        .collect();
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
        let prepared = model
            .prepare_relaxation_modes(&branches, options(), cx)
            .unwrap();
        // On q=[1,i], Cq=i*g*q. These original coefficients give the exact
        // root -.25+i*sqrt(31)/4; the conjugate has the opposite handedness.
        for sign in [-1.0, 1.0] {
            let expected = C64::new(-0.25, sign * 31.0_f64.sqrt() / 4.0);
            let pole = prepared
                .refine(
                    expected.scale(1.02),
                    &[C64::ONE, C64::new(0.0, sign)],
                    &[C64::ZERO; 4],
                    cx,
                )
                .unwrap();
            close(pole.laplace_pole_per_s, expected, 3e-11);
            close(
                pole.structural_displacement[1] / pole.structural_displacement[0],
                C64::new(0.0, sign),
                3e-11,
            );
            assert_eq!(pole.structural_dissipation_w, 0.0);
            assert!(pole.memory.iter().all(|m| m.dissipation_w > 0.0));
            original_equations_and_envelope(&model, &branches, &pole);
        }
    });
}

#[test]
fn empty_memory_matches_the_original_full_quadratic_cavity_pole() {
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 0.2,
        sound_speed_m_s: 1.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 1, medium).unwrap();
    let interior = uniform(&geometry, medium);
    let mass = dense(&[1.0], 1);
    let stiffness = dense(&[0.7], 1);
    let damping = dense(&[0.2], 1);
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
        let expected = C64::new(-0.1, -0.99_f64.sqrt());
        let potential = C64::from_re(-0.6) / expected;
        let reference = model
            .prepare_resonance_modes(options(), cx)
            .unwrap()
            .refine(expected.scale(1.03), &[C64::ONE], &[potential; 4], cx)
            .unwrap();
        let mode = model
            .prepare_relaxation_modes(&[], options(), cx)
            .unwrap()
            .refine(expected.scale(1.03), &[C64::ONE], &[potential; 4], cx)
            .unwrap();
        close(mode.laplace_pole_per_s, reference.laplace_pole_per_s, 3e-11);
        for (a, b) in mode
            .structural_displacement
            .iter()
            .zip(&reference.structural_displacement)
        {
            close(*a, *b, 3e-11);
        }
        for (a, b) in mode.pressure_pa.iter().zip(&reference.pressure_pa) {
            close(*a, *b, 3e-11);
        }
        assert!(mode.memory.is_empty());
        assert_eq!(mode.relaxation_energy_j, 0.0);
        assert_eq!(mode.relaxation_dissipation_w, 0.0);
        original_equations_and_envelope(&model, &[], &mode);
    });
}
