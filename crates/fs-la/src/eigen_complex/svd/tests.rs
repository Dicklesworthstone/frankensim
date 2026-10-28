//! G0/G3: complex reconstruction, unitary vectors, scaling, and typed refusals.
use super::*;

fn inner(a: &[C64], rows: usize, cols: usize, p: usize, q: usize) -> C64 {
    (0..rows).fold(C64::ZERO, |s, i| {
        s + a[i * cols + p].conj() * a[i * cols + q]
    })
}

fn assert_unitary(a: &[C64], rows: usize, cols: usize, tolerance: f64) {
    for p in 0..cols {
        for q in 0..cols {
            let target = if p == q { C64::ONE } else { C64::ZERO };
            assert!(
                (inner(a, rows, cols, p, q) - target).abs() <= tolerance,
                "columns {p},{q}: {:?}",
                inner(a, rows, cols, p, q)
            );
        }
    }
}

fn assert_decomposition(a: &[C64], m: usize, n: usize, d: &ComplexSvd) {
    assert_unitary(&d.u, m, n, 2e-12);
    assert_unitary(&d.v, n, n, 2e-12);
    assert!(d.singular_values.windows(2).all(|p| p[0] >= p[1]));
    let scale = a.iter().map(|z| z.abs()).fold(0.0_f64, f64::max);
    if scale == 0.0 {
        assert!(d.singular_values.iter().all(|&s| s == 0.0));
        return;
    }
    for i in 0..m {
        for j in 0..n {
            let reconstructed = (0..n).fold(C64::ZERO, |sum, k| {
                sum + (d.u[i * n + k] * d.v[j * n + k].conj()).scale(d.singular_values[k] / scale)
            });
            assert!(
                (reconstructed - divide(a[i * n + j], scale)).abs() < 2e-12,
                "entry {i},{j}: {reconstructed:?} vs {:?}",
                divide(a[i * n + j], scale)
            );
        }
    }
}

fn known_rectangular() -> Vec<C64> {
    // Orthonormal columns of a 4x4 Hadamard matrix, with nontrivial row phases.
    let u = [
        C64::new(0.5, 0.0),
        C64::new(0.5, 0.0),
        C64::new(0.5, 0.0),
        C64::new(0.0, 0.5),
        C64::new(0.0, -0.5),
        C64::new(0.0, 0.5),
        C64::new(-0.5, 0.0),
        C64::new(-0.5, 0.0),
        C64::new(0.5, 0.0),
        C64::new(0.0, -0.5),
        C64::new(0.0, 0.5),
        C64::new(0.0, 0.5),
    ];
    let w = C64::new(-0.5, det::sqrt(3.0) / 2.0);
    let v = [
        C64::ONE,
        C64::ONE,
        C64::ONE,
        C64::ONE,
        w,
        w.conj(),
        C64::ONE,
        w.conj(),
        w,
    ]
    .map(|z| divide(z, det::sqrt(3.0)));
    let sigma = [9.0, 3.0, 0.125];
    let mut a = vec![C64::ZERO; 12];
    for i in 0..4 {
        for j in 0..3 {
            for k in 0..3 {
                a[i * 3 + j] = a[i * 3 + j] + (u[i * 3 + k] * v[j * 3 + k].conj()).scale(sigma[k]);
            }
        }
    }
    a
}

#[test]
fn rectangular_complex_reconstruction_and_known_singular_values() {
    let a = known_rectangular();
    let d = svd_complex_with_cancel(&a, 4, 3, 1e-13, 64, || false).unwrap();
    assert_decomposition(&a, 4, 3, &d);
    for (&actual, expected) in d.singular_values.iter().zip([9.0, 3.0, 0.125]) {
        assert!((actual / expected - 1.0).abs() < 1e-13);
    }
    assert!(d.max_correlation <= 1e-13);
    assert_eq!(
        d,
        svd_complex_with_cancel(&a, 4, 3, 1e-13, 64, || false).unwrap()
    );
}

#[test]
fn complex_phase_and_exact_nullspace_receive_unitary_vectors() {
    let a = [
        C64::ONE,
        C64::new(0.0, 1.0),
        C64::ZERO,
        C64::ZERO,
        C64::ZERO,
        C64::ZERO,
        C64::ZERO,
        C64::ZERO,
        C64::ZERO,
        C64::ZERO,
        C64::ZERO,
        C64::ZERO,
    ];
    let d = svd_complex_with_cancel(&a, 4, 3, 1e-13, 64, || false).unwrap();
    assert_decomposition(&a, 4, 3, &d);
    assert!((d.singular_values[0] - det::sqrt(2.0)).abs() < 1e-14);
    assert_eq!(&d.singular_values[1..], &[0.0, 0.0]);
}

#[test]
fn repeated_and_all_zero_singular_values_have_deterministic_completion() {
    let a = [
        C64::ZERO,
        C64::ZERO,
        C64::ZERO,
        C64::ZERO,
        C64::new(0.0, 2.0),
        C64::ZERO,
        C64::ZERO,
        C64::ZERO,
        C64::from_re(2.0),
        C64::ZERO,
        C64::ZERO,
        C64::ZERO,
    ];
    let d = svd_complex_with_cancel(&a, 4, 3, 1e-13, 64, || false).unwrap();
    assert_eq!(d.singular_values, [2.0, 2.0, 0.0]);
    assert_eq!(d.sweeps, 0);
    assert_decomposition(&a, 4, 3, &d);
    let zero = [C64::ZERO; 12];
    let d = svd_complex_with_cancel(&zero, 4, 3, 1e-13, 64, || false).unwrap();
    assert_decomposition(&zero, 4, 3, &d);
    for i in 0..4 {
        for j in 0..3 {
            assert_eq!(d.u[i * 3 + j], if i == j { C64::ONE } else { C64::ZERO });
        }
    }
}

#[test]
fn uniform_extreme_scales_preserve_the_complex_decomposition() {
    let a = known_rectangular();
    for scale in [1e-250, 1e250] {
        let scaled: Vec<_> = a.iter().map(|z| z.scale(scale)).collect();
        let d = svd_complex_with_cancel(&scaled, 4, 3, 1e-13, 64, || false).unwrap();
        assert_decomposition(&scaled, 4, 3, &d);
        for (&actual, expected) in d.singular_values.iter().zip([9.0, 3.0, 0.125]) {
            assert!((actual / scale / expected - 1.0).abs() < 2e-13);
        }
    }
}

#[test]
fn weak_correlated_column_survives_beyond_squared_norm_range() {
    let a = [
        C64::ONE,
        C64::new(0.0, 1e-200),
        C64::ZERO,
        C64::from_re(1e-200),
        C64::ZERO,
        C64::ZERO,
    ];
    let d = svd_complex_with_cancel(&a, 3, 2, 1e-13, 64, || false).unwrap();
    assert_decomposition(&a, 3, 2, &d);
    assert!((d.singular_values[0] - 1.0).abs() < 1e-14);
    assert!((d.singular_values[1] / 1e-200 - 1.0).abs() < 1e-13);
    assert!(d.max_correlation <= 1e-13);
    // Check the weak column relative to its OWN scale, not the largest entry.
    for i in 0..3 {
        let reproduced = (0..2).fold(C64::ZERO, |s, k| {
            s + d.u[i * 2 + k] * (d.v[2 + k].conj().scale(d.singular_values[k] / 1e-200))
        });
        assert!((reproduced - divide(a[i * 2 + 1], 1e-200)).abs() < 1e-13);
    }
}

#[test]
fn bounded_nonconvergence_is_explicit() {
    let a = known_rectangular();
    let failure = svd_complex_with_cancel(&a, 4, 3, 1e-14, 1, || false).unwrap_err();
    assert!(
        matches!(failure,ComplexSvdError::NotConverged{sweeps:1,max_correlation} if max_correlation>1e-14)
    );
}

#[test]
fn cancellation_during_pair_work_does_not_change_input_or_replay() {
    let a: Vec<_> = (0..96)
        .map(|i| {
            C64::new(
                ((i * 17 + 3) % 29) as f64 - 14.0,
                ((i * 11 + 5) % 31) as f64 - 15.0,
            )
        })
        .collect();
    let original = a.clone();
    let mut polls = 0;
    let error = svd_complex_with_cancel(&a, 12, 8, 1e-13, 64, || {
        polls += 1;
        polls == 89
    })
    .unwrap_err();
    assert_eq!(error, ComplexSvdError::Cancelled);
    assert_eq!(polls, 89);
    assert_eq!(a, original);
    let first = svd_complex_with_cancel(&a, 12, 8, 1e-13, 64, || false).unwrap();
    assert_decomposition(&a, 12, 8, &first);
    assert_eq!(
        first,
        svd_complex_with_cancel(&a, 12, 8, 1e-13, 64, || false).unwrap()
    );
}

#[test]
fn admission_and_unrepresentable_ranges_are_typed() {
    for (a, m, n, tol, sweeps) in [
        (vec![C64::ONE], 1, 2, 1e-13, 64),
        (vec![], usize::MAX, 2, 1e-13, 64),
        (vec![], 1025, 0, 1e-13, 64),
        (vec![C64::new(f64::NAN, 0.0)], 1, 1, 1e-13, 64),
        (vec![C64::ONE], 1, 1, 0.0, 64),
        (vec![C64::ONE], 1, 1, 1.0, 64),
        (vec![C64::ONE], 1, 1, 1e-13, 0),
        (
            vec![C64::from_re(f64::MAX), C64::from_re(f64::MIN_POSITIVE)],
            2,
            1,
            1e-13,
            64,
        ),
        (
            vec![C64::from_re(f64::MAX), C64::from_re(f64::MAX)],
            2,
            1,
            1e-13,
            64,
        ),
    ] {
        assert!(matches!(
            svd_complex_with_cancel(&a, m, n, tol, sweeps, || false),
            Err(ComplexSvdError::BadInput { .. })
        ));
    }
    let empty = svd_complex_with_cancel(&[], 0, 0, 1e-13, 64, || false).unwrap();
    assert!(empty.u.is_empty() && empty.v.is_empty() && empty.singular_values.is_empty());
    assert_eq!(
        svd_complex_with_cancel(&[C64::ONE], 1, 1, 1e-13, 64, || true),
        Err(ComplexSvdError::Cancelled)
    );
}
