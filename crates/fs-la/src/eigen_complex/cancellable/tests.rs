use super::*;
use super::super::{EigFailure, lu_complex};

// Frozen pre-cancellation scalar operation graph; independent of the new API.
fn reference_factor(a: &[C64], n: usize) -> Result<LuComplex, EigFailure> {
    assert_eq!(a.len(), n * n, "a must be n*n = {}", n * n);
    let mut m = a.to_vec();
    let mut perm: Vec<usize> = (0..n).collect();
    for k in 0..n {
        let (mut best, mut best_i) = (m[k * n + k].abs(), k);
        for i in k + 1..n {
            let v = m[i * n + k].abs();
            if v > best {
                best = v;
                best_i = i;
            }
        }
        if best == 0.0 {
            return Err(EigFailure { window_hi: k });
        }
        if best_i != k {
            for c in 0..n {
                m.swap(k * n + c, best_i * n + c);
            }
            perm.swap(k, best_i);
        }
        let piv = m[k * n + k];
        for i in k + 1..n {
            let f = m[i * n + k] / piv;
            m[i * n + k] = f;
            for c in k + 1..n {
                let sub = f * m[k * n + c];
                m[i * n + c] = m[i * n + c] - sub;
            }
        }
    }
    Ok(LuComplex { n, data: m, perm })
}

impl LuComplex {
    /// Solve A·x = b in place.
    fn reference_solve(&self, b: &mut [C64]) {
        let n = self.n;
        assert_eq!(b.len(), n, "b must have length {n}");
        let pb: Vec<C64> = self.perm.iter().map(|&p| b[p]).collect();
        b.copy_from_slice(&pb);
        for i in 0..n {
            let mut v = b[i];
            for (k, &bk) in b.iter().enumerate().take(i) {
                let sub = self.data[i * n + k] * bk;
                v = v - sub;
            }
            b[i] = v;
        }
        for i in (0..n).rev() {
            let mut v = b[i];
            for (k, &bk) in b.iter().enumerate().take(n).skip(i + 1) {
                let sub = self.data[i * n + k] * bk;
                v = v - sub;
            }
            b[i] = v / self.data[i * n + i];
        }
    }
}

fn bits(a: &[C64]) -> Vec<(u64, u64)> {
    a.iter().map(|v| (v.re.to_bits(), v.im.to_bits())).collect()
}
fn matrix(n: usize) -> Vec<C64> {
    // Reverse a diagonally dominant complex matrix's rows to force pivoting.
    (0..n).rev().flat_map(|i| (0..n).map(move |j| {
        if i == j { C64::new(4.0 * n as f64 + 1.0, 0.5) }
        else { C64::new((1 + (i + 2 * j) % 7) as f64 / 16.0, (i as f64 - j as f64) / 32.0) }
    })).collect()
}

#[test]
fn cancellable_factor_and_solve_match_the_original_bits_and_known_solution() {
    for n in [0, 1, 2, 7, 17] {
        let a = matrix(n);
        let original = reference_factor(&a, n).unwrap();
        let mut polls = 0;
        let factor = lu_complex_with_cancel(&a, n, || { polls += 1; false }).unwrap();
        assert_eq!(polls, 2 + n * (n + 1) / 2);
        assert_eq!(factor.perm, original.perm);
        assert_eq!(bits(&factor.data), bits(&original.data));
        let want: Vec<_> = (0..n).map(|i| C64::new(i as f64 - 2.0, 0.25)).collect();
        let rhs: Vec<_> = (0..n).map(|i| (0..n)
            .fold(C64::ZERO, |sum, j| sum + a[i * n + j] * want[j])).collect();
        let mut expected = rhs.clone();
        original.reference_solve(&mut expected);
        let mut actual = rhs.clone();
        factor.solve_with_cancel(&mut actual, || false).unwrap();
        assert_eq!(bits(&actual), bits(&expected));
        for (actual, want) in actual.iter().zip(&want) {
            assert!((*actual - *want).abs() < 1e-12 * (1.0 + want.abs()));
        }
        let mut legacy = rhs;
        lu_complex(&a, n).unwrap().solve(&mut legacy);
        assert_eq!(bits(&legacy), bits(&expected));
    }
    // Equal-magnitude complex pivots must keep the lowest row index.
    let a = [C64::new(0.0, 1.0), C64::from_re(2.0),
        C64::new(0.0, -1.0), C64::from_re(3.0)];
    assert_eq!(lu_complex_with_cancel(&a, 2, || false).unwrap().perm, [0, 1]);
}

#[test]
fn cancellation_at_every_factor_checkpoint_publishes_no_partial_factor() {
    for n in [0, 1, 8] {
        let a = matrix(n);
        let before = bits(&a);
        let mut total = 0;
        lu_complex_with_cancel(&a, n, || { total += 1; false }).unwrap();
        for stop in 1..=total {
            let mut polls = 0;
            let error = lu_complex_with_cancel(&a, n, || { polls += 1; polls == stop }).unwrap_err();
            assert_eq!(error, LuComplexError::Cancelled);
            assert_eq!(polls, stop);
            assert_eq!(bits(&a), before);
        }
    }
}

#[test]
fn cancelled_solves_are_atomic_at_every_checkpoint_and_factors_remain_reusable() {
    for n in [0, 1, 8] {
        let factor = lu_complex(&matrix(n), n).unwrap();
        let before_factor = bits(&factor.data);
        let rhs: Vec<_> = (0..n).map(|i| C64::new(-(i as f64), 0.25)).collect();
        let mut total = 0;
        let mut want = rhs.clone();
        factor.solve_with_cancel(&mut want, || { total += 1; false }).unwrap();
        assert_eq!(total, 2 + 2 * n);
        for stop in 1..=total {
            let mut actual = rhs.clone();
            let mut polls = 0;
            let error = factor.solve_with_cancel(&mut actual, || { polls += 1; polls == stop }).unwrap_err();
            assert_eq!(error, LuComplexError::Cancelled);
            assert_eq!(polls, stop);
            assert_eq!(bits(&actual), bits(&rhs));
            assert_eq!(bits(&factor.data), before_factor);
            factor.solve_with_cancel(&mut actual, || false).unwrap();
            assert_eq!(bits(&actual), bits(&want));
        }
    }
}

#[test]
fn singularity_retains_its_pivot_index_and_is_not_reported_as_cancellation() {
    let a = [2.0, 1.0, 4.0, 2.0].map(C64::from_re);
    assert_eq!(lu_complex_with_cancel(&a, 2, || false).unwrap_err(), LuComplexError::Singular { index: 1 });
    assert_eq!(lu_complex(&a, 2).unwrap_err().window_hi, 1);
    assert_eq!(lu_complex_with_cancel(&a, 2, || true).unwrap_err(), LuComplexError::Cancelled);
}

/// Deterministic complex entries in [-1, 1) from a 64-bit LCG.
fn random(len: usize, seed: u64) -> Vec<C64> {
    let mut s = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    let mut next = || {
        s = s.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
        (s >> 11) as f64 / (1u64 << 52) as f64 - 1.0
    };
    (0..len).map(|_| C64::new(next(), next())).collect()
}

/// A collocation-like exterior Helmholtz matrix: e^{ikr}/(4 pi r) between
/// points on a helix, a dominant complex diagonal, k r up to about 20.
fn bem_like(n: usize) -> Vec<C64> {
    let p: Vec<[f64; 3]> = (0..n)
        .map(|i| {
            let t = i as f64 * 0.37;
            [t.cos() * 0.1, t.sin() * 0.1, i as f64 * 0.004]
        })
        .collect();
    let k = 30.0;
    (0..n * n)
        .map(|e| {
            let (i, j) = (e / n, e % n);
            if i == j {
                return C64::new(0.5, 0.05);
            }
            let r = ((p[i][0] - p[j][0]).powi(2) + (p[i][1] - p[j][1]).powi(2) + (p[i][2] - p[j][2]).powi(2)).sqrt();
            let g = 0.01 / (4.0 * core::f64::consts::PI * r);
            C64::new(g * (k * r).cos(), g * (k * r).sin())
        })
        .collect()
}

#[test]
fn fma_capsule_is_bitwise_the_original_operation_graph() {
    // On an x86-64 FMA host the dispatcher takes the capsule; this compares it
    // with the frozen baseline graph. Elsewhere it re-checks the portable body.
    #[cfg(target_arch = "x86_64")]
    assert_eq!(super::fma::accelerated(), std::arch::is_x86_feature_detected!("fma"));
    for (n, a) in [(64, random(64 * 64, 1)), (200, random(200 * 200, 2)), (150, bem_like(150))] {
        let original = reference_factor(&a, n).unwrap();
        let factor = lu_complex_with_cancel(&a, n, || false).unwrap();
        assert_eq!(factor.perm, original.perm, "n = {n}");
        assert_eq!(bits(&factor.data), bits(&original.data), "n = {n}");
        let rhs = random(n, 3);
        let mut expected = rhs.clone();
        original.reference_solve(&mut expected);
        let mut actual = rhs;
        factor.solve_with_cancel(&mut actual, || false).unwrap();
        assert_eq!(bits(&actual), bits(&expected), "n = {n}");
    }
}
