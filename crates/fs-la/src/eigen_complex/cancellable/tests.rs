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

/// Ok factors match in permutation and data bits; failures in pivot index.
fn assert_same(a: &[C64], n: usize, block: usize, threads: usize) {
    let original = reference_factor(a, n);
    let blocked = lu_complex_blocked_with(a, n, block, threads, || false);
    let at = format!("n = {n}, block = {block}, threads = {threads}");
    match (original, blocked) {
        (Ok(original), Ok(blocked)) => {
            assert_eq!(blocked.perm, original.perm, "{at}");
            assert_eq!(bits(&blocked.data), bits(&original.data), "{at}");
        }
        (Err(original), Err(blocked)) => {
            assert_eq!(blocked, LuComplexError::Singular { index: original.window_hi }, "{at}");
        }
        (original, blocked) => panic!("{at}: original {original:?}, blocked {blocked:?}"),
    }
}

/// Unit-magnitude entries in {1, i, -1, -i}: exact pivot ties at every
/// candidate of the first column, and rows i and i + 4 start identical.
fn ties(n: usize) -> Vec<C64> {
    let units = [C64::ONE, C64::new(0.0, 1.0), C64::from_re(-1.0), C64::new(0.0, -1.0)];
    (0..n * n).map(|e| units[(7 * (e / n) + 3 * (e % n)) % 4]).collect()
}

#[test]
fn blocked_elimination_is_bitwise_the_original_graph_at_every_block_and_thread_count() {
    let small = [(0, matrix(0)), (1, matrix(1)), (7, matrix(7)), (17, matrix(17)),
        (40, ties(40)), (64, random(64 * 64, 4))];
    for (n, a) in &small {
        for block in [1, 2, 3, 16, 64, 1000] {
            for threads in [1, 2, 3, 8] {
                assert_same(a, *n, block, threads);
            }
        }
    }
    for (n, a) in [(150, bem_like(150)), (201, random(201 * 201, 5))] {
        for (block, threads) in [(3, 3), (16, 8), (64, 1), (64, 3), (128, 5)] {
            assert_same(&a, n, block, threads);
        }
        let public = lu_complex_blocked_with_cancel(&a, n, 4, || false).unwrap();
        let original = reference_factor(&a, n).unwrap();
        assert_eq!(bits(&public.data), bits(&original.data), "n = {n}");
        let rhs = random(n, 6);
        let (mut expected, mut actual) = (rhs.clone(), rhs);
        original.reference_solve(&mut expected);
        public.solve_with_cancel(&mut actual, || false).unwrap();
        assert_eq!(bits(&actual), bits(&expected), "n = {n}");
    }
    // Equal-magnitude complex pivots keep the lowest row index.
    let a = [C64::new(0.0, 1.0), C64::from_re(2.0), C64::new(0.0, -1.0), C64::from_re(3.0)];
    assert_eq!(lu_complex_blocked_with(&a, 2, 1, 2, || false).unwrap().perm, [0, 1]);
}

#[test]
fn blocked_singularity_keeps_its_pivot_index_inside_any_panel() {
    // Column 5 is zero, so every path stops at pivot 5: in the second panel
    // for block 4, the first for block 8, and alone for block 1.
    let n = 12;
    let mut a = matrix(n);
    for i in 0..n {
        a[i * n + 5] = C64::ZERO;
    }
    for (block, threads) in [(1, 1), (4, 3), (8, 2), (64, 4)] {
        assert_eq!(
            lu_complex_blocked_with(&a, n, block, threads, || false).unwrap_err(),
            LuComplexError::Singular { index: 5 }
        );
        assert_same(&a, n, block, threads);
    }
    assert_eq!(lu_complex_blocked_with_cancel(&a, n, 2, || true).unwrap_err(), LuComplexError::Cancelled);
}

#[test]
fn blocked_cancellation_at_every_checkpoint_publishes_no_partial_factor() {
    let (n, block) = (10, 4);
    let a = matrix(n);
    let before = bits(&a);
    let mut total = 0;
    lu_complex_blocked_with(&a, n, block, 3, || { total += 1; false }).unwrap();
    // Input and publication, every pivot, and the upper and trailing phases of
    // the two blocks [0, 4) and [4, 8) that have a trailing part.
    assert_eq!(total, 2 + n + 2 * 2);
    for stop in 1..=total {
        let mut polls = 0;
        let error = lu_complex_blocked_with(&a, n, block, 3, || { polls += 1; polls == stop }).unwrap_err();
        assert_eq!(error, LuComplexError::Cancelled);
        assert_eq!(polls, stop);
        assert_eq!(bits(&a), before);
    }
}

/// Timing probe, not a gate. Prints one JSON row per (n, block, threads);
/// every factor is also compared bitwise with the unblocked owner.
#[test]
#[ignore = "timing probe: --release --ignored --nocapture on a quiet many-core host"]
fn blocked_lu_timing_probe() {
    use std::time::Instant;
    for n in [2000, 4096, 6240] {
        let a = bem_like(n);
        let started = Instant::now();
        let serial = lu_complex_with_cancel(&a, n, || false).unwrap();
        let serial_s = started.elapsed().as_secs_f64();
        println!("{{\"probe\":\"lu-blocked\",\"n\":{n},\"path\":\"unblocked\",\"s\":{serial_s:.3}}}");
        for (block, threads) in [(64, 1), (32, 32), (64, 8), (64, 16), (64, 32), (128, 32)] {
            let started = Instant::now();
            let factor = lu_complex_blocked_with(&a, n, block, threads, || false).unwrap();
            let s = started.elapsed().as_secs_f64();
            assert_eq!(factor.perm, serial.perm);
            assert!(bits(&factor.data) == bits(&serial.data), "n = {n}, block = {block}, threads = {threads}");
            println!(
                "{{\"probe\":\"lu-blocked\",\"n\":{n},\"block\":{block},\"threads\":{threads},\"s\":{s:.3},\"speedup\":{:.2},\"bitwise\":true}}",
                serial_s / s
            );
        }
    }
}

#[test]
fn lu_complex_is_the_original_graph_under_any_declared_worker_count() {
    use super::super::with_lu_workers;
    let (n, a) = (150, bem_like(150));
    let original = reference_factor(&a, n).unwrap();
    for workers in [None, Some(1), Some(3), Some(7)] {
        let factor = match workers {
            None => lu_complex(&a, n).unwrap(),
            Some(w) => with_lu_workers(w, || lu_complex(&a, n).unwrap()),
        };
        assert_eq!(factor.perm, original.perm, "{workers:?}");
        assert_eq!(bits(&factor.data), bits(&original.data), "{workers:?}");
    }
    // The declaration is scoped and restored, also after a panic.
    let unwound = std::panic::catch_unwind(|| with_lu_workers(5, || panic!("inside")));
    assert!(unwound.is_err());
    assert_eq!(super::super::LU_WORKERS.with(core::cell::Cell::get), 1);
    assert_eq!(with_lu_workers(4, || with_lu_workers(2, || super::super::LU_WORKERS.with(core::cell::Cell::get))), 2);
    assert_eq!(lu_complex(&[C64::ZERO], 1).unwrap_err().window_hi, 0);
}
