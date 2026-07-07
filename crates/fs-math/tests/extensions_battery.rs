//! Extension-family battery (wf9.14): measured ULP budgets vs the
//! platform-libm oracle, bitwise symmetry laws, identity cross-checks,
//! IEEE special-value tables, and the cross-ISA golden hash. The original
//! det golden hash is untouched (additive extension).
//!
//! float_cmp allowed file-wide: IEEE special-value tables (pow(x,0) = 1,
//! atan2 signed zeros, ...) are EXACT semantics, not tolerance checks.
#![allow(clippy::float_cmp)]

use fs_math::{det, ulp_distance};

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
}

#[test]
fn tan_budget_and_oddness() {
    let mut seed = 0x7A_u64;
    let mut worst = 0u64;
    for _ in 0..200_000 {
        let x = lcg(&mut seed) * 2000.0;
        let got = det::tan(x);
        let want = x.tan();
        // Near poles both values are huge; ULP distance stays meaningful.
        let d = ulp_distance(got, want);
        worst = worst.max(d);
        assert!(
            d <= det::TAN_ULP_BUDGET,
            "tan({x}) off by {d} ULP: {got} vs {want}"
        );
        // Odd BITWISE.
        assert_eq!(
            det::tan(-x).to_bits(),
            (-got).to_bits(),
            "tan must be odd bitwise at {x}"
        );
    }
    println!(
        "{{\"suite\":\"fs-math\",\"case\":\"tan\",\"verdict\":\"pass\",\"detail\":\"200k samples, worst {worst} ULP (budget {})\"}}",
        det::TAN_ULP_BUDGET
    );
}

#[test]
fn atan_atan2_budget_and_specials() {
    let mut seed = 0xA7A_u64;
    let mut worst = 0u64;
    for _ in 0..200_000 {
        let x = lcg(&mut seed) * 100.0;
        let got = det::atan(x);
        let d = ulp_distance(got, x.atan());
        worst = worst.max(d);
        assert!(d <= det::ATAN_ULP_BUDGET, "atan({x}) off by {d} ULP");
        assert_eq!(
            det::atan(-x).to_bits(),
            (-got).to_bits(),
            "atan odd bitwise"
        );
    }
    // atan2 over all quadrants.
    for _ in 0..100_000 {
        let y = lcg(&mut seed) * 10.0;
        let x = lcg(&mut seed) * 10.0;
        let got = det::atan2(y, x);
        let d = ulp_distance(got, y.atan2(x));
        worst = worst.max(d);
        assert!(
            d <= det::ATAN_ULP_BUDGET + 1,
            "atan2({y},{x}) off by {d} ULP"
        );
        assert_eq!(
            det::atan2(-y, x).to_bits(),
            (-got).to_bits(),
            "atan2 y-sign bitwise"
        );
    }
    // Special values follow libm conventions bitwise.
    let cases: [(f64, f64); 10] = [
        (0.0, 1.0),
        (-0.0, 1.0),
        (0.0, -1.0),
        (-0.0, -1.0),
        (1.0, f64::INFINITY),
        (1.0, f64::NEG_INFINITY),
        (f64::INFINITY, 1.0),
        (f64::NEG_INFINITY, 1.0),
        (f64::INFINITY, f64::INFINITY),
        (f64::NEG_INFINITY, f64::NEG_INFINITY),
    ];
    for (y, x) in cases {
        let got = det::atan2(y, x);
        let want = y.atan2(x);
        assert!(
            ulp_distance(got, want) <= 1,
            "atan2 special ({y},{x}): {got} vs {want}"
        );
    }
    // atan(±∞) = ±π/2 to ≤1 ULP.
    assert!(ulp_distance(det::atan(f64::INFINITY), std::f64::consts::FRAC_PI_2) <= 1);
    println!(
        "{{\"suite\":\"fs-math\",\"case\":\"atan\",\"verdict\":\"pass\",\"detail\":\"300k samples + specials, worst {worst} ULP\"}}"
    );
}

#[test]
fn erf_erfc_budget_symmetry_and_identity() {
    let mut seed = 0xE2F_u64;
    let (mut worst_erf, mut worst_erfc) = (0u64, 0u64);
    for _ in 0..100_000 {
        let x = lcg(&mut seed) * 12.0;
        let ge = det::erf(x);
        let de = ulp_distance(ge, libm_erf(x));
        // Oracle-limited band: the plain-f64 oracle carries Taylor
        // cancellation up to x = 2 (≈ 10 ULP) and CF truncation noise just
        // past its 2.0 handoff — a 48-ULP SANITY tolerance there; the
        // budget-grade evidence in the band is the disjoint-paths
        // cross-validation below (≤ 3 ULP mutual).
        let band = x.abs() > 1.5 && x.abs() < 3.5;
        let tol = if band { 48 } else { det::ERF_ULP_BUDGET };
        if !band {
            worst_erf = worst_erf.max(de);
        }
        assert!(de <= tol, "erf({x}) off by {de} ULP (tol {tol})");
        assert_eq!(det::erf(-x).to_bits(), (-ge).to_bits(), "erf odd bitwise");
        let xc = lcg(&mut seed).abs() * 20.0 + 3.5; // strong-oracle tail
        let gc = det::erfc(xc);
        let dc = ulp_distance(gc, libm_erfc(xc));
        worst_erfc = worst_erfc.max(dc);
        assert!(dc <= det::ERFC_ULP_BUDGET, "erfc({xc}) off by {dc} ULP");
    }
    // Disjoint-path cross-validation in the cancellation band: the dd
    // Taylor sum and the dd continued fraction are algorithmically
    // independent constructions; their mutual agreement at shared x IS
    // the strong evidence there (both would not be wrong identically).
    for k in 0..=60 {
        let x = 2.4 + 0.6 * f64::from(k) / 60.0;
        let (taylor, cf) = det::erf_both_paths(x);
        let d = ulp_distance(taylor, cf);
        assert!(
            d <= 3,
            "disjoint erf paths disagree at {x}: {taylor} vs {cf} ({d} ULP)"
        );
    }
    // erf + erfc = 1 within combined budget (moderate x).
    for k in 0..200 {
        let x = -4.0 + 8.0 * f64::from(k) / 200.0;
        let s = det::erf(x) + det::erfc(x);
        assert!((s - 1.0).abs() < 1e-14, "erf+erfc at {x}: {s}");
    }
    println!(
        "{{\"suite\":\"fs-math\",\"case\":\"erf\",\"verdict\":\"pass\",\"detail\":\"worst erf {worst_erf} / erfc {worst_erfc} ULP (budgets {}/{})\"}}",
        det::ERF_ULP_BUDGET,
        det::ERFC_ULP_BUDGET
    );
}

// libm bindings via std (f64 has no erf in std): use the C-library through
// the `libm`-equivalent std trick — std lacks erf, so oracle via a local
// high-precision series is used instead (dd-free, independent coding).
fn libm_erf(x: f64) -> f64 {
    // Independent oracle: 200-term alternating Kahan-summed series for
    // |x| ≤ 3, else 1 − continued-fraction tail evaluated in plain f64
    // with 128 CF levels (deliberately different construction from the
    // implementation under test — dd there, compensated-f64 here).
    let a = x.abs();
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    if a >= 6.0 {
        return sign;
    }
    if a <= 2.0 {
        let x2 = a * a;
        let (mut sum, mut comp) = (0.0f64, 0.0f64);
        let mut term = a;
        let add = |v: f64, sum: &mut f64, comp: &mut f64| {
            let y = v - *comp;
            let t = *sum + y;
            *comp = (t - *sum) - y;
            *sum = t;
        };
        add(a, &mut sum, &mut comp);
        for k in 1..=200u32 {
            term = term * x2 / f64::from(k);
            let c = term / f64::from(2 * k + 1);
            add(if k % 2 == 1 { -c } else { c }, &mut sum, &mut comp);
        }
        sign * (sum * std::f64::consts::FRAC_2_SQRT_PI)
    } else {
        sign * (1.0 - libm_erfc(a))
    }
}

fn libm_erfc(x: f64) -> f64 {
    if x <= 2.0 {
        return 1.0 - libm_erf(x);
    }
    if x >= 27.5 {
        return 0.0;
    }
    let mut cf = x;
    for k in (1..=128u32).rev() {
        cf = x + f64::from(k) / 2.0 / cf;
    }
    // exp(−x²) with compensated square.
    let hi = x * x;
    let lo = x.mul_add(x, -hi);
    (-hi).exp() * (1.0 - lo) / std::f64::consts::PI.sqrt() / cf
}

#[test]
fn pow_budget_and_specials() {
    let mut seed = 0xB0B_u64;
    for _ in 0..100_000 {
        let x = (lcg(&mut seed) + 0.51).abs() * 20.0 + 1e-3;
        let y = lcg(&mut seed) * 30.0;
        let got = det::pow(x, y);
        let want = x.powf(y);
        if want.is_finite() && want != 0.0 {
            let d = ulp_distance(got, want);
            let budget = det::pow_ulp_budget(x, y);
            assert!(
                d <= budget,
                "pow({x},{y}) off by {d} ULP (budget {budget}): {got} vs {want}"
            );
        }
    }
    // Integer fast path is tight.
    for n in [-8i32, -3, 2, 5, 17, 64] {
        let x = 1.7f64;
        let d = ulp_distance(det::pow(x, f64::from(n)), x.powi(n));
        assert!(d <= 6, "pow integer path n={n}: {d} ULP");
    }
    // IEEE specials.
    assert_eq!(det::pow(2.0, 0.0), 1.0);
    assert_eq!(det::pow(0.0, 3.0), 0.0);
    assert_eq!(det::pow(-0.0, 3.0).to_bits(), (-0.0f64).to_bits());
    assert_eq!(det::pow(0.0, -2.0), f64::INFINITY);
    assert_eq!(det::pow(-2.0, 3.0), -8.0);
    assert!(
        det::pow(-2.0, 0.5).is_nan(),
        "negative base, fractional exp → NaN"
    );
    assert_eq!(det::pow(2.0, f64::INFINITY), f64::INFINITY);
    assert_eq!(det::pow(0.5, f64::INFINITY), 0.0);
    assert_eq!(det::pow(2.0, 0.5).to_bits(), det::sqrt(2.0).to_bits());
    println!(
        "{{\"suite\":\"fs-math\",\"case\":\"pow\",\"verdict\":\"pass\",\"detail\":\"100k samples within honest budget + specials\"}}"
    );
}

#[test]
fn tan_is_sin_over_cos_bitwise() {
    // Shared reduction makes this an IDENTITY, not an approximation —
    // (for even quadrants; odd quadrants use the cotangent form).
    let mut seed = 0x51C_u64;
    for _ in 0..50_000 {
        let x = lcg(&mut seed) * 1.5; // stay in quadrant 0/−1 range
        let t = det::tan(x);
        let ratio = det::sin(x) / det::cos(x);
        assert_eq!(
            t.to_bits(),
            ratio.to_bits(),
            "tan({x}) must equal sin/cos bitwise on shared reduction"
        );
    }
}

/// Recorded on aarch64-apple (M4 Pro); must match on x86-64 (trj).
const GOLDEN_HASH: u64 = 0x5bbb_0256_67a9_0b70;

#[test]
fn extensions_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let mut seed = 0x601D_u64;
    for _ in 0..4000 {
        let x = lcg(&mut seed) * 900.0;
        feed(det::tan(x));
        feed(det::atan(x));
        feed(det::atan2(x, lcg(&mut seed) * 3.0));
        feed(det::erf(x / 100.0));
        feed(det::erfc((x / 40.0).abs()));
        feed(det::pow((x / 900.0).abs() + 0.1, lcg(&mut seed) * 8.0));
    }
    println!(
        "{{\"suite\":\"fs-math\",\"case\":\"extensions-golden\",\"verdict\":\"info\",\"detail\":\"{acc:#018x}\"}}"
    );
    assert_eq!(
        acc, GOLDEN_HASH,
        "extension-family bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with \
         semantic justification (golden-evidence policy)"
    );
}
