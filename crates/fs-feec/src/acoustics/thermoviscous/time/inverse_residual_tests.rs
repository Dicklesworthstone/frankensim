use super::*;

// Dense symmetric positive definite: M_ij = 1 / (1 + |i - j|), a convex decreasing Toeplitz sequence.
fn spd(base: usize) -> Vec<f64> {
    let mut mass = vec![0.0; base * base];
    for i in 0..base {
        for j in 0..base {
            mass[i * base + j] = 1.0 / (1.0 + i.abs_diff(j) as f64);
        }
    }
    mass
}

// The inverse as `prepare_time_*` forms it: one solve per unit column, stored with row stride `n`.
fn inverse_of(mass: &[f64], base: usize, n: usize) -> (Cholesky, Vec<f64>) {
    let factor = cholesky(mass, base).unwrap();
    let mut inverse = vec![0.0; n * n];
    for col in 0..base {
        let mut rhs = vec![0.0; base];
        rhs[col] = 1.0;
        factor.solve(&mut rhs);
        for row in 0..base {
            inverse[row * n + col] = rhs[row];
        }
    }
    (factor, inverse)
}

fn accepted(mass: &[f64], inverse: &[f64], factor: &Cholesky, base: usize, n: usize) -> bool {
    mass_inverse_residual_is_small(mass, inverse, factor, base, n, &mut || false).unwrap()
}

#[test]
fn the_cholesky_inverse_of_a_symmetric_matrix_is_accepted_at_a_padded_stride() {
    let (base, n) = (12, 15);
    let mass = spd(base);
    let (factor, inverse) = inverse_of(&mass, base, n);
    assert!(accepted(&mass, &inverse, &factor, base, n));
}

// Planted negative: an inverse entry that is wrong by 1e-9 (the allowance is 512 eps n times the factor bound, about 1e-11
// at n = 12) is refused, off the diagonal and on it, and so is a nonfinite entry.
#[test]
fn an_inverse_entry_that_is_wrong_by_1e_9_or_not_finite_is_refused() {
    let (base, n) = (12, 15);
    let mass = spd(base);
    let (factor, inverse) = inverse_of(&mass, base, n);
    for (row, col) in [(3, 5), (4, 4)] {
        let mut wrong = inverse.clone();
        wrong[row * n + col] += 1e-9;
        assert!(!accepted(&mass, &wrong, &factor, base, n), "({row}, {col})");
    }
    let mut nan = inverse.clone();
    nan[2 * n + 7] = f64::NAN;
    assert!(!accepted(&mass, &nan, &factor, base, n));
}

// Planted negative: the factorization reads the lower triangle only, so an upper-triangle entry of the mass that is wrong by
// 1e-9 relative, with the lower triangle unchanged, leaves the inverse that of the lower triangle; the residual formed with
// the full matrix must refuse it.
#[test]
fn an_upper_triangle_entry_that_is_wrong_by_1e_9_relative_is_refused() {
    let (base, n) = (12, 15);
    let mass = spd(base);
    let (factor, inverse) = inverse_of(&mass, base, n);
    let mut tampered = mass.clone();
    tampered[2 * base + 7] *= 1.0 + 1e-9;
    assert_eq!(tampered[7 * base + 2], mass[7 * base + 2]);
    assert!(!accepted(&tampered, &inverse, &factor, base, n));
    assert!(accepted(&mass, &inverse, &factor, base, n));
}

#[test]
fn a_cancelled_check_is_an_error_and_not_a_verdict() {
    let (base, n) = (12, 12);
    let mass = spd(base);
    let (factor, inverse) = inverse_of(&mass, base, n);
    assert!(
        mass_inverse_residual_is_small(&mass, &inverse, &factor, base, n, &mut || true).is_err()
    );
}
