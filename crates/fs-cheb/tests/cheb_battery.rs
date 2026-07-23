//! fs-cheb battery: machine-precision recovery with expected degree
//! growth, algebra/calculus identities, plateau robustness, rootfinding,
//! collocation-matrix spectral accuracy, the Dirichlet eigen demo, and
//! the per-ISA deterministic golden hash (the x86-64 row remains armed).

use fs_cheb::cheb2::Cheb2;
use fs_cheb::colleague::{ColleaguePolicy, certified_roots, colleague_roots};
use fs_cheb::{Cheb1, diff_matrix, dirichlet_laplace_eigs, lobatto_points};
use std::panic::{AssertUnwindSafe, catch_unwind};

#[test]
fn extreme_finite_domains_do_not_overflow_the_affine_map() {
    let a = -f64::MAX;
    let b = f64::MAX;
    let constant = Cheb1::from_coeffs(a, b, vec![2.0]);
    for x in [a, 0.0, b] {
        assert_eq!(constant.eval(x), 1.0, "constant failed at {x}");
    }

    let reference_linear = Cheb1::from_coeffs(a, b, vec![0.0, 1.0]);
    assert_eq!(reference_linear.eval(a), -1.0);
    assert_eq!(reference_linear.eval(0.0), 0.0);
    assert_eq!(reference_linear.eval(b), 1.0);
    let derivative = reference_linear.differentiate().eval(0.0);
    assert_eq!(derivative.to_bits(), (1.0 / f64::MAX).to_bits());

    // Coefficient and domain scaling must be combined before either side can
    // overflow: MAX*T1(x/MAX) is the physical identity f(x)=x.
    let physical_identity = Cheb1::from_coeffs(a, b, vec![0.0, f64::MAX]);
    assert!((physical_identity.eval(1.0) - 1.0).abs() <= f64::EPSILON);
    assert_eq!(physical_identity.differentiate().eval(0.0), 1.0);

    // A finite width does not make the traditional doubled numerator safe.
    // Both same-sign and one-sided extreme intervals must take the exceptional
    // center/radius path for large interior coordinates.
    let one_sided = Cheb1::from_coeffs(0.0, f64::MAX, vec![0.0, 1.0]);
    let one_sided_value = one_sided.eval(0.75 * f64::MAX);
    assert!(one_sided_value.is_finite() && (one_sided_value - 0.5).abs() <= f64::EPSILON);
    let same_sign = Cheb1::from_coeffs(f64::MAX / 4.0, f64::MAX, vec![0.0, 1.0]);
    assert!((same_sign.eval(0.875 * f64::MAX) - 2.0 / 3.0).abs() <= 2.0 * f64::EPSILON);

    // The tiny reference-coordinate offset 1/MAX corresponds to physical
    // x=1. Fractional/convex affine formulas absorb it and report a false root
    // at zero; the center/radius path must preserve it end to end.
    let unit_root = Cheb1::from_coeffs(a, b, vec![-2.0 / f64::MAX, 1.0]);
    assert!(unit_root.eval(1.0).abs() <= f64::from_bits(1));
    let colleague = colleague_roots(&unit_root, ColleaguePolicy::default());
    assert_eq!(colleague.len(), 1);
    assert!((colleague[0] - 1.0).abs() <= f64::EPSILON);
    let certified = certified_roots(&unit_root, 1e-12);
    assert!(
        certified
            .iter()
            .any(|root| root.is_certified() && root.interval().contains(1.0)),
        "the physical certificate must contain x=1: {certified:?}"
    );

    let rebuilt = Cheb1::build(
        &|x| {
            assert!(x.is_finite() && a <= x && x <= b);
            1.0
        },
        a,
        b,
        16,
    );
    assert_eq!(rebuilt.eval(a), 1.0);
    assert_eq!(rebuilt.eval(b), 1.0);
    let zero = Cheb1::from_coeffs(a, b, vec![0.0]);
    assert_eq!(zero.integral(), 0.0);

    let min_subnormal = f64::from_bits(1);
    let subnormal_constant = Cheb1::from_coeffs(0.0, min_subnormal, vec![2.0]);
    assert_eq!(
        subnormal_constant.integral().to_bits(),
        min_subnormal.to_bits()
    );

    // The reference-coordinate sum is 4/3*MAX and overflows before a tiny
    // domain radius is applied, although the bounded polynomial and physical
    // integral are both representable.
    let min_normal = f64::MIN_POSITIVE;
    let scaled_integral = Cheb1::from_coeffs(0.0, min_normal, vec![f64::MAX, 0.0, -f64::MAX / 2.0]);
    assert!((scaled_integral.integral() - 8.0 / 3.0).abs() < 2e-15);

    // A normalized fallback must also survive an overflowing prefix followed
    // by cancellation when the final physical integral is representable.
    let cancellation_integral = Cheb1::from_coeffs(
        -0.8,
        0.8,
        vec![
            f64::MAX,
            0.0,
            -f64::MAX / 2.0,
            0.0,
            f64::MAX,
            0.0,
            f64::MAX,
            0.0,
            f64::MAX,
        ],
    );
    let expected = (8.0 / 9.0) * f64::MAX;
    assert!(((cancellation_integral.integral() - expected) / expected).abs() <= 4.0 * f64::EPSILON);

    // Every naive prefix is finite, but (1 + 2^-54) rounds back to one before
    // the later -1 arrives. The exact-real integral is the representable
    // residual 2^-54 and must survive expansion summation.
    let absorbed_residual = fs_math::det::powi(2.0, -54);
    let finite_prefix_cancellation = Cheb1::from_coeffs(
        -1.0,
        1.0,
        vec![1.0, 0.0, -3.0 * fs_math::det::powi(2.0, -55), 0.0, 7.5],
    );
    assert_eq!(
        finite_prefix_cancellation.integral().to_bits(),
        absorbed_residual.to_bits()
    );

    let tiny = Cheb1::from_coeffs(-1.0, 1.0, vec![-2.0 * 0.123e-200, 1e-200]);
    assert!(
        tiny.roots()
            .iter()
            .any(|root| (*root - 0.123).abs() < 1e-12),
        "sign-change detection must not multiply two tiny values"
    );

    let reference_root = 0.123_456_789;
    let extreme_simple = Cheb1::from_coeffs(a, b, vec![-2.0 * reference_root, 1.0]);
    let simple_roots = extreme_simple.roots();
    assert_eq!(simple_roots.len(), 1);
    let expected_physical = reference_root * f64::MAX;
    assert!(
        ((simple_roots[0] - expected_physical) / expected_physical).abs() <= 4.0 * f64::EPSILON
    );

    // A rounded zero of a shifted triple root is not an accurate root
    // certificate. The fallback scanner must refuse its ill-conditioning
    // instead of magnifying a reference error across the extreme domain.
    let r = reference_root;
    let triple = Cheb1::from_coeffs(
        a,
        b,
        vec![
            -3.0 * r - 2.0 * r * r * r,
            0.75 + 3.0 * r * r,
            -1.5 * r,
            0.25,
        ],
    );
    assert!(
        catch_unwind(AssertUnwindSafe(|| triple.roots())).is_err(),
        "multiple/ill-conditioned roots require an error-bearing root API"
    );

    let scaled_colleague = Cheb1::from_coeffs(-1.0, 1.0, vec![f64::MAX, 0.0, f64::MAX]);
    let roots = colleague_roots(&scaled_colleague, ColleaguePolicy::default());
    assert_eq!(roots.len(), 2);
    assert!((roots[0] + 0.5).abs() < 1e-12 && (roots[1] - 0.5).abs() < 1e-12);

    // Trimming is relative to the mathematical series, where the stored c0
    // is halved. Scaling by the stored c0 would incorrectly turn this linear
    // polynomial into a constant and panic rather than return no in-domain
    // candidates.
    let trim_boundary = Cheb1::from_coeffs(-1.0, 1.0, vec![1.0, 7.5e-14]);
    assert!(
        colleague_roots(&trim_boundary, ColleaguePolicy::default()).is_empty(),
        "the retained linear root is far outside the reference domain"
    );

    assert!(
        catch_unwind(AssertUnwindSafe(|| zero.roots())).is_err(),
        "an identically zero polynomial needs a continuum-valued result type"
    );

    let surface = Cheb2::build(&|_, _| 1.0, (a, b, a, b), 1e-12, 2, 16);
    assert_eq!(surface.eval(1.0, -1.0), 1.0);

    // u*v underflows if Cheb2 commits to that pair before applying 1/pivot.
    let small_surface = Cheb2::build(&|_, _| 1e-200, (-1.0, 1.0, -1.0, 1.0), 1e-12, 2, 16);
    assert_eq!(small_surface.rank(), 1);
    assert!((small_surface.eval(0.25, -0.5) / 1e-200 - 1.0).abs() < 1e-12);
}

#[test]
fn algebra_refuses_distinct_domains_and_unrepresentable_coefficients() {
    let unit = Cheb1::from_coeffs(0.0, 1.0, vec![2.0]);
    let adjacent_endpoint =
        Cheb1::from_coeffs(0.0, f64::from_bits(1.0_f64.to_bits() + 1), vec![2.0]);
    assert!(
        catch_unwind(AssertUnwindSafe(|| unit.add(&adjacent_endpoint))).is_err(),
        "algebra must not silently identify adjacent but distinct domains"
    );
    assert!(
        catch_unwind(AssertUnwindSafe(|| unit.mul(&adjacent_endpoint))).is_err(),
        "algebra must not resample an operand on a different domain"
    );

    let maximal = Cheb1::from_coeffs(0.0, 1.0, vec![f64::MAX]);
    assert!(
        catch_unwind(AssertUnwindSafe(|| maximal.add(&maximal))).is_err(),
        "Cheb1 algebra must preserve its finite-coefficient invariant"
    );
}

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

/// LAST ADMITTED VALUE (bead sj31i.55 repair, 2026-07-16): commit 6ee7267
/// deliberately moved reference-coordinate root-refinement and exceptional
/// integral bits (its own golden comment declared this value invalidated and
/// armed); d93fa29 then completed the finite-term-overflow fallback, touching
/// only a previously-panicking path, so no additional bits moved. The
/// unchanged DCT oracle, recovery, calculus, root, and eigenvalue tests stay
/// green. FOUR-QUADRANT VERIFICATION at 0fba65d: aarch64 M4 Pro debug and
/// release, x86-64 ts1 5975WX debug and release — all bit-identical at this
/// value. Never copy a value from a stale local binary.
/// HISTORY: radix-2 0xaee4_8002_1eea_9097 (M4 Pro + trj); radix-4/2
/// 0x22e7_ea21_58c9_e587 (M4 Pro only); radix-8/4/2 pre-root-refinement
/// 0xeea0_4b0a_01de_46cd (M4 Pro + ts1 + ts2, bead obq0).
const GOLDEN_HASH: u64 = 0x5d2f_e305_ce90_06fb;

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
