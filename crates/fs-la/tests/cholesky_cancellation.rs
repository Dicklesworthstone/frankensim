//! G1/G3/G4: the controlled factor preserves the original arithmetic and
//! independently solves a manufactured system without partial publication.

use fs_la::factor::{FactorError, FactorRunError, cholesky, cholesky_cancellable};

fn matrix(n: usize) -> Vec<f64> {
    let mut a = vec![0.0; n * n];
    for i in 0..n {
        a[i * n + i] = 4.0 + (i % 7) as f64 * 0.125;
        for j in i.saturating_sub(2)..i {
            let entry = if i - j == 1 { -0.75 } else { 0.125 };
            a[i * n + j] = entry;
            a[j * n + i] = entry;
        }
    }
    a
}

fn product(a: &[f64], x: &[f64]) -> Vec<f64> {
    a.chunks_exact(x.len())
        .map(|row| row.iter().zip(x).map(|(a, x)| a * x).sum())
        .collect()
}

#[test]
fn g3_controlled_factor_and_solve_retain_original_bits_across_panel_boundaries() {
    for n in [0, 1, 31, 32, 33, 65, 97, 192] {
        let a = matrix(n);
        let original = cholesky(&a, n).unwrap();
        let controlled = cholesky_cancellable(&a, n, || false).unwrap();
        for i in 0..n {
            for j in 0..=i {
                assert_eq!(original.l(i, j).to_bits(), controlled.l(i, j).to_bits());
            }
        }
        let expected: Vec<_> = (0..n)
            .map(|i| ((i * 17 % 23) as f64 - 11.0) / 8.0)
            .collect();
        let rhs = if n == 0 {
            vec![]
        } else {
            product(&a, &expected)
        };
        let mut original_solution = rhs.clone();
        let mut controlled_solution = rhs;
        original.solve(&mut original_solution);
        controlled
            .solve_cancellable(&mut controlled_solution, || false)
            .unwrap();
        for ((a, b), expected) in original_solution
            .iter()
            .zip(&controlled_solution)
            .zip(expected)
        {
            assert_eq!(a.to_bits(), b.to_bits());
            assert!((b - expected).abs() < 3e-14, "n={n}: {b} != {expected}");
        }
    }
}

#[test]
fn g4_cancellation_inside_factorization_and_triangular_solves_is_atomic() {
    let n = 97;
    let a = matrix(n);
    let before = a.clone();
    let mut calls = 0;
    let factor = cholesky_cancellable(&a, n, || {
        calls += 1;
        false
    })
    .unwrap();
    assert!(calls > n);
    for stop in [1, 2, calls / 2, calls] {
        let mut count = 0;
        assert!(matches!(
            cholesky_cancellable(&a, n, || {
                count += 1;
                count == stop
            }),
            Err(FactorRunError::Cancelled)
        ));
        assert_eq!(a, before);
    }
    let rhs: Vec<_> = (0..n).map(|i| i as f64 - 32.0).collect();
    let mut scratch = rhs.clone();
    let mut calls = 0;
    factor
        .solve_cancellable(&mut scratch, || {
            calls += 1;
            false
        })
        .unwrap();
    for stop in [1, 2, n + 2, calls] {
        let mut candidate = rhs.clone();
        let mut count = 0;
        assert_eq!(
            factor.solve_cancellable(&mut candidate, || {
                count += 1;
                count == stop
            }),
            Err(FactorRunError::Cancelled),
        );
        assert!(
            candidate
                .iter()
                .zip(&rhs)
                .all(|(a, b)| a.to_bits() == b.to_bits())
        );
    }
}

#[test]
fn g0_controlled_factor_retains_numerical_refusal_and_empty_cancellation() {
    assert!(matches!(
        cholesky_cancellable(&[1.0, 2.0, 2.0, 1.0], 2, || false),
        Err(FactorRunError::Factor(FactorError::NotSpd { index: 1 }))
    ));
    assert!(matches!(
        cholesky_cancellable(&[], 0, || true),
        Err(FactorRunError::Cancelled)
    ));
}
