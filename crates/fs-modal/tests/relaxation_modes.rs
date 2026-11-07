//! Independent mechanical characteristic equations for causal material loss.

use fs_math::c64::C64;
use fs_modal::relaxation_modes::{
    MaxwellDampingBranch, maxwell_eigenvalues, maxwell_eigenvalues_with_cancel,
};
use fs_modal::{ModalError, quadratic_eigenvalues};

fn branch(stiffness: f64, tau: f64, strain: &[f64]) -> MaxwellDampingBranch {
    MaxwellDampingBranch {
        stiffness,
        relaxation_time_s: tau,
        strain: strain.to_vec(),
    }
}

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual:?} != {expected:?}"
    );
}

fn polynomial_residual(pole: C64, descending: &[f64]) {
    let mut value = C64::ZERO;
    let mut scale = 0.0;
    for &coefficient in descending {
        value = value * pole + C64::from_re(coefficient);
        scale = scale * pole.abs() + coefficient.abs();
    }
    assert!(
        value.abs() <= 2e-11 * scale,
        "pole {pole:?}: {value:?} / {scale}"
    );
}

#[test]
fn material_memory_retains_three_exact_overdamped_poles() {
    // (lambda+1)(lambda+2)(lambda+4), with nonzero passive Maxwell loss.
    let poles = maxwell_eigenvalues(
        &[1.0],
        &[2.0],
        &[8.0 / 5.0],
        1,
        &[branch(12.0 / 5.0, 1.0 / 5.0, &[1.0])],
    )
    .unwrap();
    assert_eq!(poles.len(), 3);
    for (&actual, expected) in poles.iter().zip([-4.0, -2.0, -1.0]) {
        close(actual, C64::from_re(expected), 3e-12);
        polynomial_residual(actual, &[1.0, 7.0, 14.0, 8.0]);
    }
}

#[test]
fn maxwell_loss_does_not_add_equilibrium_stiffness_to_a_free_coordinate() {
    // A Maxwell material relaxes under static strain. Keeping its branch
    // stiffness in K as well would incorrectly erase the exact zero pole.
    let poles =
        maxwell_eigenvalues(&[1.0], &[0.0], &[0.0], 1, &[branch(3.0, 0.25, &[1.0])]).unwrap();
    assert_eq!(poles.len(), 3);
    for (&actual, expected) in poles.iter().zip([-3.0, -1.0, 0.0]) {
        close(actual, C64::from_re(expected), 2e-12);
    }
}

#[test]
fn free_scalar_retains_weak_and_fast_relaxation_poles_relatively() {
    // Both the cancellation-prone 2x2 formula and threshold-based deflation
    // used to fabricate an additional zero. Polynomial residuals cannot
    // detect that error because the rigid-coordinate zero is a genuine root.
    for (stiffness, tau) in [(1e-20, 1.0), (1.0, 1e-20), (1e-300, 1.0), (1.0, 1e-300)] {
        let poles =
            maxwell_eigenvalues(&[1.0], &[0.0], &[0.0], 1, &[branch(stiffness, tau, &[1.0])])
                .unwrap();
        assert_eq!(poles.iter().filter(|s| **s == C64::ZERO).count(), 1);
        for (&actual, expected) in poles[..2].iter().zip([-1.0 / tau, -stiffness * tau]) {
            assert_eq!(actual.im, 0.0);
            assert!(
                (actual.re / expected - 1.0).abs() < 4e-15,
                "{actual:?} != {expected}"
            );
        }
    }
}

#[test]
fn free_scalar_source_products_survive_extreme_units() {
    for exponent in [-800, 800] {
        let mass = 2.0_f64.powi(exponent);
        let stiffness = 2.0_f64.powi(-exponent);
        let poles = maxwell_eigenvalues(
            &[mass],
            &[0.0],
            &[0.0],
            1,
            &[branch(stiffness, 1.0, &[mass])],
        )
        .unwrap();
        // k_b*b²/m = 1, although forming b² directly over/underflows.
        close(poles[0], C64::new(-0.5, -0.75_f64.sqrt()), 3e-15);
        close(poles[1], C64::new(-0.5, 0.75_f64.sqrt()), 3e-15);
        assert_eq!(poles[2], C64::ZERO);
    }
    // Forming m*tau or k_b*tau first would lose the actual source here.
    for value in [1e-300, 1e300] {
        let poles =
            maxwell_eigenvalues(&[value], &[0.0], &[0.0], 1, &[branch(value, value, &[1.0])])
                .unwrap();
        if value < 1.0 {
            assert!((poles[0].re / (-1.0 / value) - 1.0).abs() < 4e-15);
            assert!((poles[1].re / (-value) - 1.0).abs() < 4e-15);
        } else {
            for s in &poles[..2] {
                assert!((s.re / (-0.5 / value) - 1.0).abs() < 4e-15);
                assert!((s.im.abs() - 1.0).abs() < 4e-15);
            }
        }
        assert_eq!(poles[2], C64::ZERO);
    }
}

#[test]
fn free_scalar_preserves_active_viscosity_and_true_zero_multiplicity() {
    let cases = [
        (
            3.0,
            2.0,
            [C64::new(-2.0, -1.0), C64::new(-2.0, 1.0), C64::ZERO],
        ),
        (-4.0, 6.0, [C64::ZERO, C64::ONE, C64::from_re(2.0)]),
        (-0.5, 0.5, [C64::from_re(-0.5), C64::ZERO, C64::ZERO]),
        (-1.0, 1.0, [C64::ZERO; 3]),
    ];
    for (damping, stiffness, expected) in cases {
        let poles = maxwell_eigenvalues(
            &[1.0],
            &[damping],
            &[0.0],
            1,
            &[branch(stiffness, 1.0, &[1.0])],
        )
        .unwrap();
        for (&actual, expected) in poles.iter().zip(expected) {
            close(actual, expected, 3e-15);
        }
    }
}

#[test]
fn unrepresentable_free_scalar_sources_refuse_instead_of_erasing_loss() {
    for material in [
        branch(f64::MAX, 1.0, &[f64::MAX]),
        branch(1.0, f64::from_bits(1), &[1.0]),
        branch(1e-300, 1.0, &[1e-100]),
    ] {
        assert!(matches!(
            maxwell_eigenvalues(&[1.0], &[0.0], &[0.0], 1, &[material]),
            Err(ModalError::QuadraticEig { .. })
        ));
    }
}

#[test]
fn dense_source_zero_counts_refuse_unresolved_nonzero_poles() {
    let coupled = [branch(0.5, 1e-20, &[1.0]), branch(0.5, 1e-20, &[1.0])];
    for (damping, stiffness, branches) in [
        (0.0, 0.0, coupled.to_vec()),
        (1.0, 1e-40, vec![branch(1.0, 1.0, &[1.0])]),
    ] {
        let error =
            maxwell_eigenvalues(&[1.0], &[damping], &[stiffness], 1, &branches).unwrap_err();
        assert!(
            matches!(error, ModalError::QuadraticEig { detail } if detail.contains("exact zero poles"))
        );
    }
    // Two actual rigid poles survive when all material channels are dark.
    let dark = [branch(1.0, 1.0, &[0.0]), branch(1.0, 2.0, &[0.0])];
    let poles = maxwell_eigenvalues(&[1.0], &[0.0], &[0.0], 1, &dark).unwrap();
    assert_eq!(poles.iter().filter(|s| **s == C64::ZERO).count(), 2);
}

#[test]
fn two_relaxation_times_preserve_the_complete_rational_pencil() {
    let (m, c, k) = (2.0, 0.3, 8.0);
    let branches = [branch(5.0, 0.2, &[1.0]), branch(2.0, 0.8, &[-0.7])];
    let poles = maxwell_eigenvalues(&[m], &[c], &[k], 1, &branches).unwrap();
    assert_eq!(poles.len(), 4);
    for &s in &poles {
        let a = s * s.scale(m) + s.scale(c) + C64::from_re(k);
        let d1 = C64::ONE + s.scale(0.2);
        let d2 = C64::ONE + s.scale(0.8);
        let terms = [
            a * d1 * d2,
            s.scale(5.0 * 0.2) * d2,
            s.scale(2.0 * 0.49 * 0.8) * d1,
        ];
        let residual = terms.iter().copied().fold(C64::ZERO, |a, b| a + b);
        assert!(residual.abs() < 2e-11 * terms.iter().map(|x| x.abs()).sum::<f64>());
        assert!(s.re < 0.0, "passive source grew: {s:?}");
    }
    let trace = poles.iter().map(|s| s.re).sum::<f64>();
    assert!((trace - (-c / m - 1.0 / 0.2 - 1.0 / 0.8)).abs() < 2e-12);
}

#[test]
fn shared_strain_coupling_and_gyroscopic_terms_enter_the_full_determinant() {
    let m = [2.0, 0.2, 0.2, 1.3];
    let c = [0.3, 0.4, -0.3, 0.5];
    let k = [8.0, 2.0, 2.0, 5.0];
    let b = [0.7, -0.9];
    let branches = [branch(3.0, 0.4, &b)];
    let poles = maxwell_eigenvalues(&m, &c, &k, 2, &branches).unwrap();
    assert_eq!(poles.len(), 5);
    for s in poles {
        let mut a = [C64::ZERO; 4];
        for i in 0..4 {
            a[i] = (s * s).scale(m[i]) + s.scale(c[i]) + C64::from_re(k[i]);
        }
        // Matrix determinant lemma without inverting the bare pencil:
        // det(A)*(1+tau*s) + k_b*tau*s*b^T adj(A) b = 0.
        let det = a[0] * a[3] - a[1] * a[2];
        let adj =
            a[3].scale(b[0] * b[0]) + a[0].scale(b[1] * b[1]) - (a[1] + a[2]).scale(b[0] * b[1]);
        let left = det * (C64::ONE + s.scale(0.4));
        let right = s.scale(3.0 * 0.4) * adj;
        assert!((left + right).abs() < 3e-10 * (left.abs() + right.abs()));
        assert!(s.re < 0.0);
    }
}

#[test]
fn relaxed_and_unrelaxed_oscillation_limits_follow_physical_memory() {
    for (tau, expected_omega) in [(1e-3, 2.0), (1e3, 3.0)] {
        let poles =
            maxwell_eigenvalues(&[1.0], &[0.0], &[4.0], 1, &[branch(5.0, tau, &[1.0])]).unwrap();
        let oscillation = poles.iter().find(|s| s.im > 1.0).unwrap();
        assert!((oscillation.im - expected_omega).abs() < 2e-5, "{poles:?}");
        assert!(oscillation.re < 0.0);
    }
}

#[test]
fn physical_poles_are_invariant_under_large_coordinate_unit_changes() {
    let m = [2.0, 0.2, 0.2, 1.3];
    let c = [0.3, 0.08, 0.08, 0.5];
    let k = [8.0, 2.0, 2.0, 5.0];
    let b = [0.7, -0.9];
    let reference = maxwell_eigenvalues(&m, &c, &k, 2, &[branch(3.0, 0.4, &b)]).unwrap();
    let factors = [2.0_f64.powi(350), 2.0_f64.powi(-350)];
    let transform = |a: &[f64]| -> Vec<f64> {
        a.iter()
            .enumerate()
            .map(|(i, &v)| v * factors[i / 2] * factors[i % 2])
            .collect()
    };
    let strain = [b[0] * factors[0], b[1] * factors[1]];
    let transformed = maxwell_eigenvalues(
        &transform(&m),
        &transform(&c),
        &transform(&k),
        2,
        &[branch(3.0, 0.4, &strain)],
    )
    .unwrap();
    for (&a, &b) in transformed.iter().zip(&reference) {
        close(a, b, 2e-11);
    }
}

#[test]
fn no_memory_delegates_to_existing_exact_quadratic_paths() {
    let expected = quadratic_eigenvalues(&[1.0], &[1e200], &[1.0], 1).unwrap();
    assert_eq!(
        maxwell_eigenvalues(&[1.0], &[1e200], &[1.0], 1, &[]).unwrap(),
        expected
    );
}

#[test]
fn invalid_materials_shapes_and_cancel_never_return_a_partial_spectrum() {
    for material in [
        branch(-1.0, 0.2, &[1.0]),
        branch(1.0, 0.0, &[1.0]),
        branch(1.0, 0.2, &[f64::NAN]),
        branch(1.0, 0.2, &[1.0, 2.0]),
    ] {
        assert!(matches!(
            maxwell_eigenvalues(&[1.0], &[0.1], &[4.0], 1, &[material]),
            Err(ModalError::DampingInput { .. })
        ));
    }
    assert!(maxwell_eigenvalues(&[1.0], &[0.1], &[4.0], 2, &[]).is_err());
    assert_eq!(
        maxwell_eigenvalues_with_cancel(
            &[1.0],
            &[0.1],
            &[4.0],
            1,
            &[branch(1.0, 0.2, &[1.0])],
            120,
            || true
        ),
        Err(ModalError::Cancelled)
    );
    let mut checks = 0;
    let result = maxwell_eigenvalues_with_cancel(
        &[1.0],
        &[0.1],
        &[4.0],
        1,
        &[branch(1.0, 0.2, &[1.0])],
        120,
        || {
            checks += 1;
            checks > 6
        },
    );
    assert_eq!(result, Err(ModalError::Cancelled));
    assert!(
        maxwell_eigenvalues_with_cancel(
            &[1.0],
            &[0.1],
            &[4.0],
            1,
            &[branch(1.0, 0.2, &[1.0])],
            0,
            || false
        )
        .is_err()
    );
}
