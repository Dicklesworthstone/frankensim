//! Local Maxwell modes checked against original mechanical/material equations.

use fs_material::visco::GeneralizedMaxwell;
use fs_math::c64::C64;
use fs_modal::relaxation_modes::{
    MaxwellDampingBranch, MaxwellMode, MaxwellModeError, PoleOptions, PreparedMaxwellModes,
};

const I: C64 = C64::new(0.0, 1.0);

fn options(rate: f64) -> PoleOptions {
    PoleOptions {
        spectral_scale: rate,
        step_tolerance: 2e-13,
        residual_tolerance: 2e-12,
        max_iterations: 40,
        max_dimension: 8,
        max_dense_bytes: 16 * 1024 * 1024,
        max_work: 200_000_000,
    }
}
fn branch(k: f64, tau: f64, b: &[f64]) -> MaxwellDampingBranch {
    MaxwellDampingBranch {
        stiffness: k,
        relaxation_time_s: tau,
        strain: b.to_vec(),
    }
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-280),
        "{actual:?} != {expected:?}"
    );
}
fn check_source(
    mode: &MaxwellMode,
    mass: &[f64],
    damping: &[f64],
    stiffness: &[f64],
    branches: &[MaxwellDampingBranch],
) {
    let n = mode.displacement.len();
    let s = mode.laplace_frequency_per_s;
    let mut mass_norm = C64::ZERO;
    for i in 0..n {
        let mut residual = C64::ZERO;
        let mut envelope = 0.0;
        for j in 0..n {
            mass_norm = mass_norm
                + mode.displacement[i].conj() * mode.displacement[j].scale(mass[i * n + j]);
            for term in [
                s * mode.velocity[j].scale(mass[i * n + j]),
                mode.velocity[j].scale(damping[i * n + j]),
                mode.displacement[j].scale(stiffness[i * n + j]),
            ] {
                residual = residual + term;
                envelope += term.abs();
            }
        }
        for (b, &stress) in branches.iter().zip(&mode.branch_stress) {
            let term = stress.scale(b.strain[i]);
            residual = residual + term;
            envelope += term.abs();
        }
        assert!(
            residual.abs() <= 2e-11 * envelope,
            "mechanical row {i}: {residual:?} / {envelope}"
        );
        close(mode.velocity[i], s * mode.displacement[i], 3e-15);
    }
    for (b, &stress) in branches.iter().zip(&mode.branch_stress) {
        let mut strain = C64::ZERO;
        let mut strain_rate = C64::ZERO;
        for j in 0..n {
            strain = strain + mode.displacement[j].scale(b.strain[j]);
            strain_rate = strain_rate + mode.velocity[j].scale(b.strain[j]);
        }
        let expected = (s.scale(b.relaxation_time_s) / (C64::ONE + s.scale(b.relaxation_time_s)))
            * strain.scale(b.stiffness);
        if expected == C64::ZERO {
            assert_eq!(stress, C64::ZERO);
        } else {
            close(stress, expected, 8e-12);
        }
        let stored = (C64::ONE + s.scale(b.relaxation_time_s)) * stress;
        let driven = strain_rate.scale(b.stiffness * b.relaxation_time_s);
        assert!((stored - driven).abs() <= 2e-11 * (stored.abs() + driven.abs()));
    }
    close(mass_norm, C64::ONE, 1e-13);
    assert!(mode.mechanical_relative_residual <= 2e-12);
    assert!(mode.relaxation_relative_residual <= 2e-12);
    assert!(mode.solver_relative_residual <= 2e-12);
    assert!(mode.relative_correction <= 2e-13);
}

#[test]
fn source_material_recovers_all_three_known_mechanical_decay_shapes() {
    // V=4, strain=q, added tip mass=4. The original constitutive owner gives
    // det(T)*(1+s/5) proportional to (s+1)(s+2)(s+4).
    let material = GeneralizedMaxwell::new(1.6, vec![(2.4, 0.2)]).unwrap();
    let mass = [4.0];
    let damping = [8.0];
    let stiffness = [4.0 * material.e_inf];
    let branches: Vec<_> = material
        .terms
        .iter()
        .map(|&(k, tau)| branch(4.0 * k, tau, &[1.0]))
        .collect();
    let prepared = PreparedMaxwellModes::new(
        &mass,
        &damping,
        &stiffness,
        1,
        &branches,
        options(5.0),
        || false,
    )
    .unwrap();
    for (root, stress) in [(-1.0, -1.2), (-2.0, -3.2), (-4.0, -19.2)] {
        let mode = prepared
            .refine(C64::from_re(root - 0.03), &[C64::new(7.0, -2.0)], || false)
            .unwrap();
        close(mode.laplace_frequency_per_s, C64::from_re(root), 2e-12);
        close(mode.displacement[0], C64::from_re(0.5), 2e-15);
        close(mode.branch_stress[0], C64::from_re(stress), 4e-12);
        check_source(&mode, &mass, &damping, &stiffness, &branches);
    }
}

fn coupled_source(gyro: f64) -> ([f64; 4], [f64; 4], [f64; 4], Vec<MaxwellDampingBranch>) {
    // Start from an independent quadratic pencil whose mode at -1+2i has
    // q proportional to [1,i]. Compensate the original real C and K for the
    // exact branch values at that one pole; their derivatives remain full.
    let mass = [1.0, 0.0, 0.0, 1.0];
    let mut damping = [1.0, 0.5, 0.5, 3.0];
    let mut stiffness = [5.0, 2.5, 2.5, 5.0];
    let branches = vec![branch(0.8, 0.4, &[1.0, 0.5]), branch(0.3, 0.7, &[0.3, 1.0])];
    let s = C64::new(-1.0, 2.0);
    for b in &branches {
        let h = (s.scale(b.relaxation_time_s) / (C64::ONE + s.scale(b.relaxation_time_s)))
            .scale(b.stiffness);
        let dc = -h.im / s.im;
        let dk = -h.re - s.re * dc;
        for i in 0..2 {
            for j in 0..2 {
                damping[2 * i + j] += dc * b.strain[i] * b.strain[j];
                stiffness[2 * i + j] += dk * b.strain[i] * b.strain[j];
            }
        }
    }
    // G*[1,i]=-i*gyro*[1,i], retaining the same planted physical mode.
    damping[1] -= gyro;
    damping[2] += gyro;
    for i in [0, 3] {
        damping[i] -= gyro / 2.0;
        stiffness[i] -= 2.5 * gyro;
    }
    (mass, damping, stiffness, branches)
}

#[test]
fn coupled_branches_and_gyroscopic_sign_retain_physical_complex_motion() {
    for gyro in [-0.3, 0.3] {
        let (m, c, k, branches) = coupled_source(gyro);
        let prepared =
            PreparedMaxwellModes::new(&m, &c, &k, 2, &branches, options(3.0), || false).unwrap();
        let seed = [C64::ONE, C64::new(0.02, 0.9)];
        let mode = prepared
            .refine(C64::new(-0.93, 2.08), &seed, || false)
            .unwrap();
        close(mode.laplace_frequency_per_s, C64::new(-1.0, 2.0), 4e-12);
        close(mode.displacement[1] / mode.displacement[0], I, 5e-12);
        assert_eq!(mode.branch_stress.len(), 2);
        check_source(&mode, &m, &c, &k, &branches);
        let replay = prepared
            .refine(C64::new(-0.93, 2.08), &seed, || false)
            .unwrap();
        assert_eq!(mode.laplace_frequency_per_s, replay.laplace_frequency_per_s);
        assert_eq!(mode.displacement, replay.displacement);
        assert_eq!(mode.branch_stress, replay.branch_stress);
    }
}

#[test]
fn mutual_mass_and_original_strain_projection_survive_congruence() {
    let (m, c, k, branches) = coupled_source(0.2);
    let change = [[1.0, 0.5], [0.0, 2.0]];
    let mut matrices = [[0.0; 4]; 3];
    for (order, a) in [m, c, k].iter().enumerate() {
        for i in 0..2 {
            for j in 0..2 {
                for r in 0..2 {
                    for col in 0..2 {
                        matrices[order][2 * i + j] +=
                            change[r][i] * a[2 * r + col] * change[col][j];
                    }
                }
            }
        }
    }
    let transformed: Vec<_> = branches
        .iter()
        .map(|b| {
            branch(
                b.stiffness,
                b.relaxation_time_s,
                &[b.strain[0], 0.5 * b.strain[0] + 2.0 * b.strain[1]],
            )
        })
        .collect();
    let prepared = PreparedMaxwellModes::new(
        &matrices[0],
        &matrices[1],
        &matrices[2],
        2,
        &transformed,
        options(3.0),
        || false,
    )
    .unwrap();
    let seed = [C64::new(1.0, -0.25), I.scale(0.5)];
    let mode = prepared
        .refine(C64::new(-0.95, 2.03), &seed, || false)
        .unwrap();
    close(mode.laplace_frequency_per_s, C64::new(-1.0, 2.0), 4e-12);
    close(
        mode.displacement[1] / mode.displacement[0],
        seed[1] / seed[0],
        5e-12,
    );
    check_source(
        &mode,
        &matrices[0],
        &matrices[1],
        &matrices[2],
        &transformed,
    );
}

#[test]
fn coordinate_units_and_material_rate_scales_preserve_modes_and_branch_stress() {
    let (m, c, k, branches) = coupled_source(-0.2);
    for (coordinates, rate) in [
        ([2.0_f64.powi(300), 2.0_f64.powi(-300)], 1.0),
        ([1.0, 1.0], 2.0_f64.powi(400)),
        ([1.0, 1.0], 2.0_f64.powi(-400)),
    ] {
        let mut matrices = [[0.0; 4]; 3];
        for (order, a) in [m, c, k].iter().enumerate() {
            let factor = match order {
                0 => 1.0,
                1 => rate,
                _ => rate * rate,
            };
            for i in 0..2 {
                for j in 0..2 {
                    matrices[order][2 * i + j] =
                        a[2 * i + j] * coordinates[i] * coordinates[j] * factor;
                }
            }
        }
        let transformed: Vec<_> = branches
            .iter()
            .map(|b| {
                branch(
                    (b.stiffness * rate) * rate,
                    b.relaxation_time_s / rate,
                    &[b.strain[0] * coordinates[0], b.strain[1] * coordinates[1]],
                )
            })
            .collect();
        let prepared = PreparedMaxwellModes::new(
            &matrices[0],
            &matrices[1],
            &matrices[2],
            2,
            &transformed,
            options(3.0 * rate),
            || false,
        )
        .unwrap();
        let seed = [
            C64::from_re(1.0 / coordinates[0]),
            I.scale(1.0 / coordinates[1]),
        ];
        let mode = prepared
            .refine(C64::new(-0.95, 2.03).scale(rate), &seed, || false)
            .unwrap();
        close(
            mode.laplace_frequency_per_s.scale(1.0 / rate),
            C64::new(-1.0, 2.0),
            4e-12,
        );
        let mapped = [
            mode.displacement[0].scale(coordinates[0]),
            mode.displacement[1].scale(coordinates[1]),
        ];
        close(mapped[1] / mapped[0], I, 5e-12);
        check_source(
            &mode,
            &matrices[0],
            &matrices[1],
            &matrices[2],
            &transformed,
        );
    }
}

#[test]
fn nearby_weak_mechanical_decay_is_not_replaced_by_the_rigid_zero() {
    for (stiffness, tau) in [(1e-20, 1.0), (1.0, 1e-20)] {
        let branches = [branch(stiffness, tau, &[1.0])];
        let prepared =
            PreparedMaxwellModes::new(&[1.0], &[0.0], &[0.0], 1, &branches, options(1e-20), || {
                false
            })
            .unwrap();
        let mode = prepared
            .refine(C64::from_re(-1.03e-20), &[C64::ONE], || false)
            .unwrap();
        assert!((mode.laplace_frequency_per_s.re / (-1e-20) - 1.0).abs() < 5e-13);
        assert_eq!(mode.laplace_frequency_per_s.im, 0.0);
        close(mode.branch_stress[0], C64::from_re(-1e-40), 5e-13);
        check_source(&mode, &[1.0], &[0.0], &[0.0], &branches);
    }
}

#[test]
fn active_viscosity_preserves_growing_mechanical_modes_and_material_stress() {
    // The exact characteristic polynomial is s*(s-1)*(s-2); retaining the
    // signed C is essential even though the Maxwell branch itself is passive.
    let branches = [branch(6.0, 1.0, &[1.0])];
    let prepared =
        PreparedMaxwellModes::new(&[1.0], &[-4.0], &[0.0], 1, &branches, options(2.0), || {
            false
        })
        .unwrap();
    for (root, stress) in [(1.0, 3.0), (2.0, 4.0)] {
        let mode = prepared
            .refine(C64::from_re(root + 0.02), &[C64::ONE], || false)
            .unwrap();
        close(mode.laplace_frequency_per_s, C64::from_re(root), 3e-12);
        close(mode.branch_stress[0], C64::from_re(stress), 4e-12);
        check_source(&mode, &[1.0], &[-4.0], &[0.0], &branches);
    }
}

#[test]
fn material_singularities_dark_modes_and_defective_mechanics_refuse() {
    for strain in [0.0, 1.0] {
        let branches = [branch(1.0, 0.5, &[strain])];
        let prepared =
            PreparedMaxwellModes::new(&[1.0], &[0.2], &[4.0], 1, &branches, options(2.0), || false)
                .unwrap();
        assert!(matches!(
            prepared.refine(C64::from_re(-2.0), &[C64::ONE], || false),
            Err(MaxwellModeError::MaterialPole { branch: 0, .. })
        ));
        assert!(matches!(
            prepared.refine(C64::new(-0.1, 2.0), &[C64::ZERO], || false),
            Err(MaxwellModeError::Invalid(_))
        ));
        if strain == 0.0 {
            let mode = prepared
                .refine(C64::new(-0.09, 2.02), &[C64::ONE], || false)
                .unwrap();
            close(
                mode.laplace_frequency_per_s,
                C64::new(-0.1, 3.99_f64.sqrt()),
                3e-12,
            );
            assert_eq!(mode.branch_stress, [C64::ZERO]);
            check_source(&mode, &[1.0], &[0.2], &[4.0], &branches);
        }
    }
    let repeated =
        PreparedMaxwellModes::new(&[1.0], &[2.0], &[1.0], 1, &[], options(1.0), || false).unwrap();
    assert!(matches!(
        repeated.refine(C64::from_re(-1.0), &[C64::ONE], || false),
        Err(MaxwellModeError::Solve(_))
    ));
}

#[test]
fn invalid_sources_resource_caps_and_cancellation_publish_no_partial_mode() {
    let (m, c, k, branches) = coupled_source(0.0);
    for constrained in [
        PoleOptions {
            max_dimension: 1,
            ..options(3.0)
        },
        PoleOptions {
            max_dense_bytes: 1,
            ..options(3.0)
        },
        PoleOptions {
            max_work: 1,
            ..options(3.0)
        },
    ] {
        assert!(
            PreparedMaxwellModes::new(&m, &c, &k, 2, &branches, constrained, || false).is_err()
        );
    }
    for invalid in [
        branch(1.0, 0.0, &[1.0, 0.0]),
        branch(-1.0, 0.5, &[1.0, 0.0]),
        branch(1.0, 0.5, &[1.0]),
        branch(1.0, 0.5, &[f64::NAN, 1.0]),
    ] {
        assert!(matches!(
            PreparedMaxwellModes::new(&m, &c, &k, 2, &[invalid], options(3.0), || false),
            Err(MaxwellModeError::Invalid(_))
        ));
    }
    assert!(matches!(
        PreparedMaxwellModes::new(&m, &c, &k, 2, &branches, options(3.0), || true),
        Err(MaxwellModeError::Cancelled)
    ));
    let prepared =
        PreparedMaxwellModes::new(&m, &c, &k, 2, &branches, options(3.0), || false).unwrap();
    let mut polls = 0;
    assert!(matches!(
        prepared.refine(C64::new(-0.95, 2.03), &[C64::ONE, I], || {
            polls += 1;
            polls > 10
        }),
        Err(MaxwellModeError::Cancelled)
    ));
    let mode = prepared
        .refine(C64::new(-0.95, 2.03), &[C64::ONE, I], || false)
        .unwrap();
    check_source(&mode, &m, &c, &k, &branches);
}
