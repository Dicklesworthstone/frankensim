//! Complex eigensolver battery (urvw item 1): analytic spectra
//! (companion/rotation matrices), agreement with the symmetric Jacobi
//! path on embedded real matrices, trace/determinant identities,
//! Hermitian reality, a Chebyshev-companion root fixture, and the
//! cross-ISA golden hash.

use fs_la::eigen_complex::{det_complex, eig};
use fs_math::c64::C64;

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
}

fn cmat(n: usize, mut f: impl FnMut(usize, usize) -> C64) -> Vec<C64> {
    let mut m = vec![C64::ZERO; n * n];
    for i in 0..n {
        for j in 0..n {
            m[i * n + j] = f(i, j);
        }
    }
    m
}

#[test]
fn companion_matrix_gives_roots_of_unity() {
    // Companion of z⁴ − 1: eigenvalues are the 4th roots of unity.
    let n = 4;
    let m = cmat(n, |i, j| {
        if i == 0 && j == n - 1 {
            C64::ONE // −(−1) constant coefficient
        } else if i > 0 && j == i - 1 {
            C64::ONE
        } else {
            C64::ZERO
        }
    });
    let eigs = eig(&m, n).expect("companion converges");
    // SET comparison (canonical ordering is roundoff-fragile when real
    // parts tie at ±1e-16): every expected root has exactly one close
    // computed eigenvalue.
    let want = [
        C64::new(-1.0, 0.0),
        C64::new(0.0, -1.0),
        C64::new(0.0, 1.0),
        C64::new(1.0, 0.0),
    ];
    for w in &want {
        let hits = eigs.iter().filter(|g| (**g - *w).abs() < 1e-12).count();
        assert_eq!(hits, 1, "root {w:?} not found once: {eigs:?}");
    }
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"eigc-companion\",\"verdict\":\"pass\",\"detail\":\"z^4-1 roots of unity to 1e-12\"}}"
    );
}

#[test]
fn rotation_block_gives_conjugate_pair() {
    // Real 2×2 rotation by θ has eigenvalues e^{±iθ}.
    let th = 0.7f64;
    let (c, s) = (fs_math::det::cos(th), fs_math::det::sin(th));
    let m = [
        C64::from_re(c),
        C64::from_re(-s),
        C64::from_re(s),
        C64::from_re(c),
    ];
    let eigs = eig(&m, 2).unwrap();
    assert!((eigs[0] - C64::new(c, -s)).abs() < 1e-14, "{eigs:?}");
    assert!((eigs[1] - C64::new(c, s)).abs() < 1e-14, "{eigs:?}");
}

#[test]
fn embedded_symmetric_matches_jacobi() {
    let n = 16;
    let mut seed = 0x51_u64;
    let mut a = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..=i {
            let v = lcg(&mut seed);
            a[i * n + j] = v;
            a[j * n + i] = v;
        }
    }
    let (jac, _) = fs_la::eigen::jacobi_eigh(&a, n);
    let m: Vec<C64> = a.iter().map(|&v| C64::from_re(v)).collect();
    let mut qr: Vec<f64> = eig(&m, n).unwrap().iter().map(|z| z.re).collect();
    qr.sort_by(f64::total_cmp);
    for (k, (&q, &j)) in qr.iter().zip(&jac).enumerate() {
        assert!(
            (q - j).abs() < 1e-10 * j.abs().max(1.0),
            "eig[{k}]: QR {q} vs Jacobi {j}"
        );
    }
    // Imaginary parts of a real-symmetric spectrum must vanish.
    for z in eig(&m, n).unwrap() {
        assert!(z.im.abs() < 1e-10, "symmetric spectrum must be real: {z:?}");
    }
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"eigc-vs-jacobi\",\"verdict\":\"pass\",\"detail\":\"16x16 symmetric: QR == Jacobi to 1e-10\"}}"
    );
}

#[test]
fn trace_and_determinant_identities() {
    let n = 12;
    let mut seed = 0x7A2_u64;
    let m = cmat(n, |_, _| C64::new(lcg(&mut seed), lcg(&mut seed)));
    let eigs = eig(&m, n).unwrap();
    // Σλ = tr(A).
    let mut tr = C64::ZERO;
    for i in 0..n {
        tr = tr + m[i * n + i];
    }
    let sum = eigs.iter().fold(C64::ZERO, |acc, &z| acc + z);
    assert!((sum - tr).abs() < 1e-10 * tr.abs().max(1.0), "trace: {sum:?} vs {tr:?}");
    // Πλ = det(A) (independent Gaussian-elimination oracle).
    let prod = eigs.iter().fold(C64::ONE, |acc, &z| acc * z);
    let d = det_complex(&m, n);
    assert!(
        (prod - d).abs() < 1e-8 * d.abs().max(1.0),
        "determinant: {prod:?} vs {d:?}"
    );
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"eigc-identities\",\"verdict\":\"pass\",\"detail\":\"trace + det identities on random complex 12x12\"}}"
    );
}

#[test]
fn hermitian_spectrum_is_real() {
    let n = 10;
    let mut seed = 0x4E_u64;
    let mut m = vec![C64::ZERO; n * n];
    for i in 0..n {
        m[i * n + i] = C64::from_re(lcg(&mut seed));
        for j in 0..i {
            let z = C64::new(lcg(&mut seed), lcg(&mut seed));
            m[i * n + j] = z;
            m[j * n + i] = z.conj();
        }
    }
    for z in eig(&m, n).unwrap() {
        assert!(z.im.abs() < 1e-11, "Hermitian eigenvalue not real: {z:?}");
    }
}

#[test]
fn chebyshev_companion_roots() {
    // Colleague-matrix fixture: the degree-6 Chebyshev polynomial T6 has
    // roots cos((2k+1)π/12). Its colleague matrix (Chebyshev companion)
    // is tridiagonal-plus-rank-1; for T_n itself the colleague matrix is
    // exactly the Jacobi matrix with halved corner entries.
    let n = 6;
    let mut m = vec![C64::ZERO; n * n];
    for i in 0..n - 1 {
        let half = if i == 0 { std::f64::consts::FRAC_1_SQRT_2 } else { 0.5 };
        // Standard colleague structure for pure T6: off-diagonals ½ with
        // the √2 correction on the first row/col pair.
        m[i * n + i + 1] = C64::from_re(half);
        m[(i + 1) * n + i] = C64::from_re(half);
    }
    let eigs = eig(&m, n).unwrap();
    let mut want: Vec<f64> = (0..n)
        .map(|k| fs_math::det::cos(std::f64::consts::PI * (2.0 * k as f64 + 1.0) / 12.0))
        .collect();
    want.sort_by(f64::total_cmp);
    for (got, w) in eigs.iter().zip(&want) {
        assert!(got.im.abs() < 1e-12, "T6 roots are real: {got:?}");
        assert!((got.re - w).abs() < 1e-12, "T6 root {} vs {w}", got.re);
    }
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"eigc-colleague\",\"verdict\":\"pass\",\"detail\":\"T6 colleague roots match cos((2k+1)pi/12)\"}}"
    );
}

/// Recorded on aarch64-apple (M4 Pro); must match on x86-64 (trj).
const GOLDEN_HASH: u64 = 0x0; // placeholder: set from first run

#[test]
fn eigc_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let mut seed = 0x601D_u64;
    for trial in 0..6 {
        let n = 6 + trial;
        let m = cmat(n, |_, _| C64::new(lcg(&mut seed), lcg(&mut seed)));
        for z in eig(&m, n).unwrap() {
            feed(z.re);
            feed(z.im);
        }
    }
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"eigc-golden\",\"verdict\":\"info\",\"detail\":\"{acc:#018x}\"}}"
    );
    assert_eq!(
        acc, GOLDEN_HASH,
        "complex eigensolver bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only \
         with semantic justification (golden-evidence policy)"
    );
}
