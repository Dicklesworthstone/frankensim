//! G0/G1/G3 oracles for physical incremental motion under a held preload.
use fs_nlmodal::{PreloadedModalStorage, SosModalStorage, StressChannel};
use fs_phs::{Storage, discrete_gradient};

fn scalar(omega: f64, beta: f64, reference: f64) -> PreloadedModalStorage {
    PreloadedModalStorage::new(
        SosModalStorage {
            omegas: vec![omega],
            channels: vec![StressChannel {
                coefficient: beta,
                coupling: vec![1.0],
            }],
        },
        &[reference],
    )
    .unwrap()
}

fn coupled(reference: &[f64]) -> PreloadedModalStorage {
    PreloadedModalStorage::new(
        SosModalStorage {
            omegas: vec![3.0, 4.6],
            channels: vec![
                StressChannel {
                    coefficient: 7.0,
                    coupling: vec![1.0, 0.25, 0.25, -0.5],
                },
                StressChannel {
                    coefficient: 3.0,
                    coupling: vec![0.0, 0.7, 0.7, 0.0],
                },
            ],
        },
        reference,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, relative: f64) {
    let scale = expected.abs().max(1e-290);
    assert!(
        (actual - expected).abs() <= relative * scale,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn g1_scalar_duffing_has_the_physical_preload_tangent_and_cubic_energy() {
    let (omega, beta, reference) = (3.0_f64, 20.0_f64, 0.4_f64);
    let storage = scalar(omega, beta, reference);
    let stiffness = omega * omega + 3.0 * beta * reference * reference;
    close(
        storage.preload_forces()[0],
        omega * omega * reference + beta * reference.powi(3),
        1e-15,
    );
    close(storage.local_stiffness()[0], stiffness, 1e-15);
    assert_eq!(storage.mode_count(), 1);
    assert_eq!(storage.equilibrium_displacements(), &[reference]);
    for q in [-0.3_f64, -0.02, 0.0, 0.04, 0.5] {
        let p = -0.23;
        let expected_energy = 0.5 * p * p
            + 0.5 * stiffness * q * q
            + beta * reference * q.powi(3)
            + 0.25 * beta * q.powi(4);
        close(storage.hamiltonian(&[q, p]), expected_energy, 2e-15);
        let mut effort = [0.0; 2];
        storage.gradient(&[q, p], &mut effort);
        let expected = stiffness * q + 3.0 * beta * reference * q * q + beta * q.powi(3);
        if expected == 0.0 {
            assert_eq!(effort[0], 0.0);
        } else {
            close(effort[0], expected, 3e-15);
        }
        assert_eq!(effort[1], p);
    }
    assert_eq!(storage.hamiltonian(&[0.0, 0.0]), 0.0);
}

#[test]
fn g1_bias_generates_dc_and_even_harmonics_from_the_same_restoring_force() {
    let (omega, beta, reference, amplitude) = (3.0, 20.0, 0.4, 0.13);
    let storage = scalar(omega, beta, reference);
    let mut coefficients = [0.0; 4];
    let samples = 1024;
    for k in 0..samples {
        let phase = core::f64::consts::TAU * k as f64 / samples as f64;
        let q = amplitude * phase.cos();
        let mut effort = [0.0; 2];
        storage.gradient(&[q, 0.0], &mut effort);
        coefficients[0] += effort[0] / samples as f64;
        for h in 1..=3 {
            coefficients[h] += 2.0 * effort[0] * (h as f64 * phase).cos() / samples as f64;
        }
    }
    // cos^2=(1+cos 2t)/2 and cos^3=(3 cos t+cos 3t)/4.
    let even = 1.5 * beta * reference * amplitude * amplitude;
    close(coefficients[0], even, 3e-14);
    close(coefficients[2], even, 3e-14);
    close(
        coefficients[1],
        storage.local_stiffness()[0] * amplitude + 0.75 * beta * amplitude.powi(3),
        3e-14,
    );
    close(coefficients[3], 0.25 * beta * amplitude.powi(3), 2e-12);

    let centered = scalar(omega, beta, 0.0);
    let mut even_centered = 0.0;
    for k in 0..samples {
        let phase = core::f64::consts::TAU * k as f64 / samples as f64;
        let mut effort = [0.0; 2];
        centered.gradient(&[amplitude * phase.cos(), 0.0], &mut effort);
        even_centered += 2.0 * effort[0] * (2.0 * phase).cos() / samples as f64;
    }
    assert!(even_centered.abs() < 1e-14);
}

#[test]
fn g1_translated_dense_law_retains_off_diagonal_tangent_and_source_physics() {
    let reference = [0.2, -0.13];
    let storage = coupled(&reference);
    assert!(storage.local_stiffness()[1].abs() > 0.1);
    assert_eq!(storage.local_stiffness()[1], storage.local_stiffness()[2]);
    let mut base_effort = [0.0; 4];
    let base = [reference[0], 0.0, reference[1], 0.0];
    storage.source().gradient(&base, &mut base_effort);
    for i in 0..2 {
        close(storage.preload_forces()[i], base_effort[2 * i], 1e-15);
    }
    let x = [-0.07, 0.3, 0.09, -0.2];
    let absolute = [reference[0] + x[0], x[1], reference[1] + x[2], x[3]];
    let translated_energy = storage.source().hamiltonian(&absolute)
        - storage.source().hamiltonian(&base)
        - storage.preload_forces()[0] * x[0]
        - storage.preload_forces()[1] * x[2];
    close(storage.hamiltonian(&x), translated_energy, 5e-15);
    let mut actual = [0.0; 4];
    let mut expected = [0.0; 4];
    storage.gradient(&x, &mut actual);
    storage.source().gradient(&absolute, &mut expected);
    for i in 0..2 {
        expected[2 * i] -= base_effort[2 * i];
        close(actual[2 * i], expected[2 * i], 4e-15);
        assert_eq!(actual[2 * i + 1], expected[2 * i + 1]);
    }
    let direction = [0.11, -0.2, 0.3, 0.4];
    storage.hessian_action(&x, &direction, &mut actual).unwrap();
    storage
        .source()
        .hessian_action(&absolute, &direction, &mut expected)
        .unwrap();
    for (a, e) in actual.iter().zip(expected) {
        close(*a, e, 4e-15);
    }
    let q = storage.quadratic_part().unwrap();
    assert!(!q.exact);
    assert_eq!(q.matrix[2], storage.local_stiffness()[1]);
    assert_eq!(q.matrix[8], storage.local_stiffness()[2]);
    storage
        .hessian_action(&[0.0; 4], &direction, &mut actual)
        .unwrap();
    for i in 0..4 {
        let value: f64 = (0..4).map(|j| q.matrix[4 * i + j] * direction[j]).sum();
        close(actual[i], value, 3e-15);
    }
    let boxed: Box<dyn Storage> = Box::new(storage);
    assert!(
        boxed
            .continuous_tangent()
            .unwrap()
            .apply(&x, &direction, &mut actual)
    );
    for (a, e) in actual.iter().zip(expected) {
        close(*a, e, 4e-15);
    }
}

#[test]
fn g3_tiny_increment_beside_a_large_reference_is_not_subtracted_away() {
    let reference = 1e8;
    let beta = 1e-16;
    let increment = 1e-12;
    assert_eq!(reference + increment, reference);
    let storage = scalar(2.0, beta, reference);
    let stiffness = 7.0;
    let mut effort = [0.0; 2];
    storage.gradient(&[increment, 0.0], &mut effort);
    close(effort[0], stiffness * increment, 5e-16);
    close(
        storage.hamiltonian(&[increment, 0.0]),
        0.5 * stiffness * increment * increment,
        5e-16,
    );
    let remainder = storage
        .midpoint_energy_remainder(&[0.0, 0.0], &[increment, 0.0], &[0.5 * increment, 0.0])
        .unwrap();
    close(
        remainder,
        0.25 * beta * (reference + 0.5 * increment) * increment.powi(3),
        8e-16,
    );
    assert!(remainder > 0.0);
}

#[test]
fn g1_midpoint_remainder_closes_work_against_independent_gauss_integration() {
    let storage = coupled(&[0.2, -0.13]);
    let a = [-0.11, 0.23, 0.07, -0.19];
    let b = [0.17, -0.08, -0.15, 0.13];
    let delta: [f64; 4] = core::array::from_fn(|i| b[i] - a[i]);
    let mid: [f64; 4] = core::array::from_fn(|i| f64::midpoint(a[i], b[i]));
    // Two-point Gauss quadrature integrates the cubic source effort exactly.
    let mut independently_integrated_work = 0.0;
    for t in [0.5 - 0.5 / 3.0_f64.sqrt(), 0.5 + 0.5 / 3.0_f64.sqrt()] {
        let mut absolute: [f64; 4] = core::array::from_fn(|i| a[i] + t * delta[i]);
        for i in 0..2 {
            absolute[2 * i] += storage.equilibrium_displacements()[i];
        }
        let mut effort = [0.0; 4];
        storage.source().gradient(&absolute, &mut effort);
        for i in 0..2 {
            effort[2 * i] -= storage.preload_forces()[i];
        }
        independently_integrated_work += 0.5 * (0..4).map(|i| effort[i] * delta[i]).sum::<f64>();
    }
    let discrete = discrete_gradient(&storage, &a, &b);
    let work: f64 = discrete.iter().zip(delta).map(|(e, d)| e * d).sum();
    close(work, independently_integrated_work, 2e-14);
    close(
        work,
        storage.hamiltonian(&b) - storage.hamiltonian(&a),
        2e-14,
    );
    let reverse = discrete_gradient(&storage, &b, &a);
    for (forward, backward) in discrete.iter().zip(reverse) {
        close(*forward, backward, 2e-15);
    }
    let mut effort = [0.0; 4];
    storage.gradient(&mid, &mut effort);
    let midpoint_work: f64 = effort.iter().zip(delta).map(|(e, d)| e * d).sum();
    let remainder = storage.midpoint_energy_remainder(&a, &b, &mid).unwrap();
    close(
        midpoint_work + remainder,
        independently_integrated_work,
        2e-14,
    );
}

#[test]
fn g1_represented_midpoint_offsets_retain_quadratic_cubic_and_quartic_work() {
    let eps = f64::EPSILON;
    let linear = scalar(2.0, 0.0, 0.7);
    let a = [0.0, 1.0];
    let b = [0.0, 1.0 + eps];
    let mid = [0.0, f64::midpoint(a[1], b[1])];
    assert_eq!(mid[1], a[1]);
    assert_eq!(
        linear.midpoint_energy_remainder(&a, &b, &mid),
        Some(0.5 * eps * eps)
    );

    let biased = scalar(2.0, 4.0, 0.5);
    let a = [1.0, 0.0];
    let b = [1.0 + eps, 0.0];
    let mid = [f64::midpoint(a[0], b[0]), 0.0];
    assert_eq!(mid[0], a[0]);
    // Independent Taylor expansion at q=1: H''=31, H'''/6=6,
    // H''''/24=1. The omitted leading quadratic term is O(eps^2),
    // whereas the ideal-centered quartic-only expression is O(eps^3).
    let expected = 15.5 * eps * eps + 6.0 * eps.powi(3) + eps.powi(4);
    close(
        biased.midpoint_energy_remainder(&a, &b, &mid).unwrap(),
        expected,
        3e-15,
    );

    let dense = PreloadedModalStorage::new(
        SosModalStorage {
            omegas: vec![1.0, 2.0],
            channels: vec![StressChannel {
                coefficient: 2.0,
                coupling: vec![1.0, 0.5, 0.5, -0.25],
            }],
        },
        &[0.25, -0.5],
    )
    .unwrap();
    let a = [1.0, 1.0, 2.0, -1.0];
    let b = [1.0 + eps, 1.0 + eps, 2.0 + 2.0 * eps, -1.0 + 0.5 * eps];
    let mid = core::array::from_fn::<_, 4, _>(|i| f64::midpoint(a[i], b[i]));
    assert_eq!(mid, a);
    // Exact dyadic expansion along displacement direction [1,2] and
    // momentum direction [1,1/2]: 219/8 eps^2+10 eps^3+2 eps^4.
    let expected = (219.0 / 8.0) * eps * eps + 10.0 * eps.powi(3) + 2.0 * eps.powi(4);
    close(
        dense.midpoint_energy_remainder(&a, &b, &mid).unwrap(),
        expected,
        3e-15,
    );
    close(
        dense.midpoint_energy_remainder(&b, &a, &mid).unwrap(),
        -expected,
        3e-15,
    );

    // The exact midpoint offset also retains a tiny endpoint beside a unit
    // endpoint; subtracting the two rounded half-intervals loses it.
    let remainder = linear
        .midpoint_energy_remainder(&[0.0, 1e-200], &[0.0, 1.0], &[0.0, 0.5])
        .unwrap();
    close(remainder, 0.5e-200, 2e-15);
}

#[test]
fn g1_semidefinite_tangent_does_not_certify_a_nonlinear_energy_minimum() {
    let storage = PreloadedModalStorage::new(
        SosModalStorage {
            omegas: vec![1.0, 1.0],
            channels: vec![StressChannel {
                coefficient: 1.0,
                coupling: vec![1.0, 2.0, 2.0, 0.0],
            }],
        },
        &[1.0, 0.0],
    )
    .unwrap();
    assert_eq!(storage.preload_forces(), &[2.0, 2.0]);
    assert_eq!(storage.local_stiffness(), &[4.0, 6.0, 6.0, 9.0]);
    // The Hessian annihilates [3,-2], but the exact original incremental
    // energy on that direction is 15 e^3 + 225 e^4/4. Negative e arbitrarily
    // close to zero lowers energy: numerical PSD admission is only curvature.
    for eps in [-0.01_f64, -0.005, -0.0025] {
        let energy = storage.hamiltonian(&[3.0 * eps, 0.0, -2.0 * eps, 0.0]);
        let expected = 15.0 * eps.powi(3) + (225.0 / 4.0) * eps.powi(4);
        close(energy, expected, 3e-13);
        assert!(energy < 0.0);
    }
}

#[test]
fn g3_centered_limit_is_the_original_law_and_exactness_is_not_inferred_from_bias() {
    let storage = coupled(&[0.0, 0.0]);
    for state in [[0.2, -0.3, 0.11, 0.04], [-0.1, 0.0, -0.3, 0.23]] {
        close(
            storage.hamiltonian(&state),
            storage.source().hamiltonian(&state),
            2e-15,
        );
        let mut a = [0.0; 4];
        let mut b = [0.0; 4];
        storage.gradient(&state, &mut a);
        storage.source().gradient(&state, &mut b);
        for (actual, expected) in a.into_iter().zip(b) {
            if expected == 0.0 {
                assert_eq!(actual, 0.0);
            } else {
                close(actual, expected, 2e-15);
            }
        }
    }
    assert!(!storage.quadratic_part().unwrap().exact);
    let linear = scalar(2.0, 0.0, 0.7);
    assert!(linear.quadratic_part().unwrap().exact);
    assert_eq!(linear.preload_forces(), &[2.8]);
    assert_eq!(
        linear.midpoint_energy_remainder(&[0.0; 2], &[0.1, -0.2], &[0.05, -0.1]),
        Some(0.0)
    );
}

#[test]
fn g0_admission_refuses_bad_references_asymmetry_overflow_and_negative_curvature() {
    let source = || SosModalStorage {
        omegas: vec![1.0, 1.0],
        channels: vec![StressChannel {
            coefficient: 1.0,
            coupling: vec![1.0, 0.0, 0.0, -1.0],
        }],
    };
    assert!(PreloadedModalStorage::new(source(), &[0.0]).is_err());
    assert!(PreloadedModalStorage::new(source(), &[f64::NAN, 0.0]).is_err());
    assert!(PreloadedModalStorage::new(source(), &[f64::MAX, 0.0]).is_err());
    assert!(matches!(
        PreloadedModalStorage::new(source(), &[0.0, 2.0]),
        Err(fs_nlmodal::NlModalError::Phs(
            fs_phs::PhsError::NotPsd { .. }
        ))
    ));
    let neutral = PreloadedModalStorage::new(source(), &[0.0, 1.0]).unwrap();
    assert_eq!(neutral.local_stiffness()[0], 0.0);
    let mut asymmetric = source();
    asymmetric.channels[0].coupling[1] = 1e-12;
    assert!(asymmetric.validate().is_ok());
    assert!(PreloadedModalStorage::new(asymmetric, &[0.0, 0.0]).is_err());
    let mut out = [17.0; 4];
    assert!(
        neutral
            .hessian_action(&[f64::NAN; 4], &[0.0; 4], &mut out)
            .is_err()
    );
    assert_eq!(out, [17.0; 4]);
    neutral.gradient(&[0.0; 2], &mut out);
    assert!(out.iter().all(|v| v.is_nan()));
}
