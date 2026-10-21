//! Timing evidence for the complex-LU FMA capsule (WhiteAnchor mail 5404).
//!
//! `lu_complex_fma N...` factors and solves one deterministic diagonally
//! dominant complex matrix per N twice: with the ORIGINAL loop compiled here at
//! baseline codegen (each `C64` product is two `f64::mul_add`, a library `fma()`
//! call on baseline x86-64), and with `fs_la::eigen_complex::lu_complex`, which
//! dispatches to the FMA-compiled capsule when the CPU has FMA. It prints both
//! wall times, their rates over the n^3/3 complex row updates (8 flops each),
//! and whether the two solutions are bitwise identical.
use fs_la::eigen_complex::lu_complex;
use fs_math::c64::C64;
use std::time::Instant;

fn matrix(n: usize) -> Vec<C64> {
    (0..n * n)
        .map(|e| {
            let (i, j) = (e / n, e % n);
            if i == j {
                C64::new(4.0 * n as f64 + 1.0, 0.5)
            } else {
                C64::new(
                    (1 + (i + 2 * j) % 7) as f64 / 16.0,
                    (i as f64 - j as f64) / 32.0,
                )
            }
        })
        .collect()
}

/// The pre-capsule factorization, operation for operation.
fn baseline(a: &[C64], n: usize) -> (Vec<C64>, Vec<usize>) {
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
    (m, perm)
}

fn baseline_solve(data: &[C64], perm: &[usize], n: usize, b: &[C64]) -> Vec<C64> {
    let mut x: Vec<C64> = perm.iter().map(|&p| b[p]).collect();
    for i in 0..n {
        let mut v = x[i];
        for (k, &bk) in x.iter().enumerate().take(i) {
            v = v - data[i * n + k] * bk;
        }
        x[i] = v;
    }
    for i in (0..n).rev() {
        let mut v = x[i];
        for (k, &bk) in x.iter().enumerate().take(n).skip(i + 1) {
            v = v - data[i * n + k] * bk;
        }
        x[i] = v / data[i * n + i];
    }
    x
}

fn main() {
    let sizes: Vec<usize> = std::env::args()
        .skip(1)
        .map(|a| a.parse().expect("N must be a positive integer"))
        .collect();
    assert!(!sizes.is_empty(), "usage: lu_complex_fma N [N...]");
    #[cfg(target_arch = "x86_64")]
    let fma = std::arch::is_x86_feature_detected!("fma");
    #[cfg(not(target_arch = "x86_64"))]
    let fma = false;
    for n in sizes {
        let a = matrix(n);
        let rhs: Vec<C64> = (0..n).map(|i| C64::new(i as f64 - 2.0, 0.25)).collect();
        let t0 = Instant::now();
        let (data, perm) = baseline(&a, n);
        let base_x = baseline_solve(&data, &perm, n, &rhs);
        let t_base = t0.elapsed().as_secs_f64();
        let t1 = Instant::now();
        let lu = lu_complex(&a, n).expect("nonsingular");
        let mut x = rhs.clone();
        lu.solve(&mut x);
        let t_new = t1.elapsed().as_secs_f64();
        let identical = x
            .iter()
            .zip(&base_x)
            .all(|(p, q)| p.re.to_bits() == q.re.to_bits() && p.im.to_bits() == q.im.to_bits());
        let flops = 8.0 * (n as f64).powi(3) / 3.0;
        println!(
            "{{\"n\":{n},\"x86_fma_detected\":{fma},\"baseline_s\":{t_base:.3},\"capsule_s\":{t_new:.3},\"baseline_gflops\":{:.3},\"capsule_gflops\":{:.3},\"speedup\":{:.2},\"solutions_bitwise_identical\":{identical}}}",
            flops / t_base / 1e9,
            flops / t_new / 1e9,
            t_base / t_new
        );
    }
}
