//! fs-cheb battery: machine-precision recovery with expected degree
//! growth, algebra/calculus identities, plateau robustness, rootfinding,
//! collocation-matrix spectral accuracy, the Dirichlet eigen demo, and
//! the cross-ISA golden hash.

use fs_cheb::{Cheb1, diff_matrix, dirichlet_laplace_eigs, lobatto_points};

#[test]
fn machine_precision_recovery_and_degree_growth() {
    // exp on [-1, 1]: entire function, tiny degree, near-eps accuracy.
    let f = Cheb1::build(&|x: f64| fs_math::det::exp(x), -1.0, 1.0, 1 << 12);
    assert!(
        f.degree() <= 20,
        "exp should resolve at low degree: {}",
        f.degree()
    );
    for k in 0..=50 {
        let x = -1.0 + 2.0 * f64::from(k) / 50.0;
        let err = (f.eval(x) - fs_math::det::exp(x)).abs();
        assert!(err < 1e-13, "exp recovery at {x}: err {err:.2e}");
    }
    // Runge 1/(1+25x²): analytic on [-1,1] but with nearby poles — needs
    // visibly more degree than exp, still converges.
    let runge = Cheb1::build(&|x: f64| 1.0 / (1.0 + 25.0 * x * x), -1.0, 1.0, 1 << 12);
    assert!(
        runge.degree() > f.degree() && runge.degree() <= 300,
        "runge degree {} outside expected band",
        runge.degree()
    );
    for k in 0..=50 {
        let x = -1.0 + 2.0 * f64::from(k) / 50.0;
        let err = (runge.eval(x) - 1.0 / (1.0 + 25.0 * x * x)).abs();
        assert!(err < 1e-12, "runge recovery at {x}: err {err:.2e}");
    }
    // Oscillatory sin(20x) on [0, 3]: degree scales with total phase.
    let osc = Cheb1::build(&|x: f64| fs_math::det::sin(20.0 * x), 0.0, 3.0, 1 << 12);
    assert!(
        osc.degree() > 35 && osc.degree() <= 200,
        "oscillatory degree {} outside expected band",
        osc.degree()
    );
    let err = (osc.eval(1.234_567) - fs_math::det::sin(20.0 * 1.234_567)).abs();
    assert!(err < 1e-11, "oscillatory recovery err {err:.2e}");
    println!(
        "{{\"suite\":\"fs-cheb\",\"case\":\"recovery\",\"verdict\":\"pass\",\"detail\":\"degrees exp={} runge={} sin20x={}\"}}",
        f.degree(),
        runge.degree(),
        osc.degree()
    );
}

#[test]
fn calculus_identities() {
    let f = Cheb1::build(&|x: f64| fs_math::det::exp(x), -1.0, 1.0, 1 << 10);
    // d/dx e^x = e^x.
    let df = f.differentiate();
    for k in 0..=20 {
        let x = -0.95 + 1.9 * f64::from(k) / 20.0;
        assert!(
            (df.eval(x) - fs_math::det::exp(x)).abs() < 1e-11,
            "derivative of exp at {x}"
        );
    }
    // ∫₋₁¹ eˣ = e − 1/e.
    let want = fs_math::det::exp(1.0) - fs_math::det::exp(-1.0);
    assert!(
        (f.integral() - want).abs() < 1e-13,
        "integral {} vs {want}",
        f.integral()
    );
    // Domain scaling: ∫₀³ sin(20x) = (1 − cos 60)/20.
    let osc = Cheb1::build(&|x: f64| fs_math::det::sin(20.0 * x), 0.0, 3.0, 1 << 12);
    let want_osc = (1.0 - fs_math::det::cos(60.0)) / 20.0;
    assert!(
        (osc.integral() - want_osc).abs() < 1e-12,
        "osc integral {} vs {want_osc}",
        osc.integral()
    );
    // Algebra: (f+f)(x) = 2f(x); (f·f)(x) = f(x)² at machine precision.
    let s = f.add(&f);
    let p = f.mul(&f);
    for k in 0..=10 {
        let x = -1.0 + 2.0 * f64::from(k) / 10.0;
        assert!((s.eval(x) - 2.0 * f.eval(x)).abs() < 1e-13);
        assert!((p.eval(x) - f.eval(x) * f.eval(x)).abs() < 1e-12);
    }
    println!(
        "{{\"suite\":\"fs-cheb\",\"case\":\"calculus\",\"verdict\":\"pass\",\"detail\":\"d/dx, integral, algebra identities\"}}"
    );
}

#[test]
fn plateau_detection_survives_noise_floor() {
    // A function evaluated with an artificial ~1e-15 deterministic jitter:
    // the plateau detector must still terminate at a reasonable degree
    // (never chase the noise floor to max_degree).
    let noisy = |x: f64| {
        let bits = x.to_bits() ^ 0x9E37_79B9_7F4A_7C15;
        let jitter = ((bits % 1000) as f64 - 500.0) * 1e-18;
        fs_math::det::exp(x) + jitter
    };
    let f = Cheb1::build(&noisy, -1.0, 1.0, 1 << 12);
    assert!(
        f.degree() <= 64,
        "noise floor chased too far: degree {}",
        f.degree()
    );
    assert!((f.eval(0.3) - fs_math::det::exp(0.3)).abs() < 1e-12);
}

#[test]
fn roots_battery() {
    // sin(20x) on [0, 3]: roots at kπ/20 for k = 0..=19 inside (0, 3)
    // plus the endpoint 0; count sign-change roots k = 1..=19.
    let osc = Cheb1::build(&|x: f64| fs_math::det::sin(20.0 * x), 0.0, 3.0, 1 << 12);
    let roots = osc.roots();
    let interior: Vec<f64> = roots.iter().copied().filter(|&r| r > 1e-9).collect();
    assert_eq!(interior.len(), 19, "sin(20x) interior roots: {interior:?}");
    for (k, &r) in interior.iter().enumerate() {
        let want = std::f64::consts::PI * ((k + 1) as f64) / 20.0;
        assert!((r - want).abs() < 1e-9, "root {k}: {r} vs {want}");
    }
    // Polynomial with known roots: (x−0.5)(x+0.25)x on [−1, 1].
    let poly = Cheb1::build(&|x: f64| (x - 0.5) * (x + 0.25) * x, -1.0, 1.0, 256);
    let r = poly.roots();
    assert_eq!(r.len(), 3, "cubic roots: {r:?}");
    for (got, want) in r.iter().zip(&[-0.25, 0.0, 0.5]) {
        assert!((got - want).abs() < 1e-10, "{got} vs {want}");
    }
    println!(
        "{{\"suite\":\"fs-cheb\",\"case\":\"roots\",\"verdict\":\"pass\",\"detail\":\"19 oscillatory + 3 cubic roots located\"}}"
    );
}

#[test]
fn collocation_matrix_is_spectrally_accurate() {
    // D applied to exp at the Lobatto points ≈ exp (spectral accuracy).
    let n = 24;
    let d = diff_matrix(n);
    let x = lobatto_points(n);
    let m = n + 1;
    let v: Vec<f64> = x.iter().map(|&t| fs_math::det::exp(t)).collect();
    for i in 0..m {
        let mut acc = 0.0f64;
        for j in 0..m {
            acc = d[i * m + j].mul_add(v[j], acc);
        }
        assert!(
            (acc - v[i]).abs() < 1e-10,
            "D·exp at node {i}: {acc} vs {}",
            v[i]
        );
    }
    // Negative-sum trick, stated bitwise-testably: the diagonal is the
    // EXACT negation of the off-diagonal sum computed in construction
    // order (a full-row re-sum in a different order picks up roundoff —
    // that is why the trick fixes accuracy in the first place).
    for i in 0..m {
        let mut off = 0.0f64;
        for j in 0..m {
            if j != i {
                off += d[i * m + j];
            }
        }
        assert_eq!(
            (off + d[i * m + i]).to_bits(),
            0.0f64.to_bits(),
            "row {i}: diagonal must exactly negate the off-diagonal sum"
        );
    }
}

#[test]
fn dirichlet_eigenvalues_match_analytic() {
    // −u″ = λu on [−1, 1], u(±1) = 0: λ_k = (kπ/2)².
    let eigs = dirichlet_laplace_eigs(48, 3);
    let want = [
        std::f64::consts::PI * std::f64::consts::PI / 4.0,
        std::f64::consts::PI * std::f64::consts::PI,
        9.0 * std::f64::consts::PI * std::f64::consts::PI / 4.0,
    ];
    for (k, (&got, &w)) in eigs.iter().zip(&want).enumerate() {
        assert!(
            (got - w).abs() < 1e-6 * w,
            "collocation eig {k}: {got} vs analytic {w}"
        );
    }
    println!(
        "{{\"suite\":\"fs-cheb\",\"case\":\"collocation-eig\",\"verdict\":\"pass\",\"detail\":\"dirichlet laplace eigs {eigs:?} vs (k*pi/2)^2\"}}"
    );
}

/// JUSTIFIED BUMP (bead 27d3, 2026-07-09): fs-fft moved its DCT substrate
/// from radix-2 to mixed radix-4/2 Stockham. The resulting butterfly and
/// twiddle order changes bits while the unchanged DCT oracle, recovery,
/// calculus, root, and eigenvalue tests remain green. Previous radix-2
/// golden: 0xaee4_8002_1eea_9097 (M4 Pro + trj). This radix-4 value is
/// recorded on M4 Pro; the x86-64 row remains armed pending RCH admission.
const GOLDEN_HASH: u64 = 0x22e7_ea21_58c9_e587;

#[test]
fn cheb_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let f = Cheb1::build(&|x: f64| fs_math::det::exp(x), -1.0, 1.0, 1 << 10);
    for &c in f.coeffs() {
        feed(c);
    }
    feed(f.integral());
    let df = f.differentiate();
    feed(df.eval(0.125));
    let osc = Cheb1::build(&|x: f64| fs_math::det::sin(20.0 * x), 0.0, 3.0, 1 << 12);
    for &r in &osc.roots() {
        feed(r);
    }
    for &e in &dirichlet_laplace_eigs(32, 2) {
        feed(e);
    }
    println!(
        "{{\"suite\":\"fs-cheb\",\"case\":\"golden-hash\",\"verdict\":\"info\",\"detail\":\"{acc:#018x}\"}}"
    );
    assert_eq!(
        acc, GOLDEN_HASH,
        "cheb bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}
