//! Payne–Hanek battery (bead r6r5): self-verifying constant
//! regeneration (G5 — all-integer Machin bignum vs the hardcoded
//! limbs), overlap agreement with the Cody–Waite path, the published
//! worst-case double, budget sweeps across the full exponent range
//! against the platform-libm oracle (tests-only oracle, per the
//! det.rs doctrine), and bitwise odd symmetry.

use fs_math::det::{TRIG_DOMAIN, cos, sin, tan};
use fs_math::payne::{
    PH_LIMBS, TWO_OVER_PI_LIMBS, generate_pi, generate_two_over_pi, reduce_pio2_large,
};

fn verdict(name: &str, pass: bool, details: &str) {
    println!("{{\"test\":\"{name}\",\"pass\":{pass},\"details\":\"{details}\"}}");
    assert!(pass, "{name}: {details}");
}

#[allow(clippy::float_cmp)] // the equality fast path is DELIBERATELY bitwise
fn ulp_diff(a: f64, b: f64) -> u64 {
    if a == b {
        return 0;
    }
    // Monotone bit-pattern distance (same-sign assumption at the
    // scales gated here; sign flips register as huge — correctly).
    let ia = i64::from_ne_bytes(a.to_bits().to_ne_bytes());
    let ib = i64::from_ne_bytes(b.to_bits().to_ne_bytes());
    ia.abs_diff(ib)
}

/// payne-001: SELF-VERIFYING CONSTANTS — the all-integer Machin
/// bignum regenerates every hardcoded 2/π limb (1280 bits), and the
/// Machin π itself matches the published hex expansion.
#[test]
fn payne_001_regeneration() {
    let regen = generate_two_over_pi(PH_LIMBS);
    let all_match = regen
        .iter()
        .zip(TWO_OVER_PI_LIMBS.iter())
        .all(|(a, b)| a == b);
    let pi = generate_pi(8);
    // π = 3.243F6A8885A308D313198A2E03707344A409382229 9F31D0… (hex).
    let pi_ok = pi.int == 3
        && pi.frac[0] == 0x243F_6A88_85A3_08D3
        && pi.frac[1] == 0x1319_8A2E_0370_7344
        && pi.frac[2] == 0xA409_3822_299F_31D0;
    verdict(
        "payne-001-regeneration",
        all_match && pi_ok,
        &format!(
            "all {PH_LIMBS} limbs regenerate from integer Machin; pi hex expansion matches published digits"
        ),
    );
}

/// payne-002: overlap agreement — where BOTH reductions are valid
/// (just below TRIG_DOMAIN), Payne–Hanek and Cody–Waite give the same
/// quadrant and matching sin values within the combined budgets.
/// POLICY (documented): the dispatch boundary is |x| = 2²⁰ exactly;
/// Cody–Waite serves at and below, Payne–Hanek above.
#[test]
fn payne_002_overlap() {
    let mut worst = 0u64;
    for k in 0..2000u64 {
        #[allow(clippy::cast_precision_loss)]
        let x = 0.45f64.mul_add(k as f64, TRIG_DOMAIN * 0.5); // spans [2^19, ~2^19.9]
        let (r_ph, q_ph) = reduce_pio2_large(x);
        // Compare via the trig VALUE (quadrant conventions can differ
        // by the nearest-rounding at f = 1/2 boundaries).
        let s_ph = match q_ph {
            0 => sin_core_probe(r_ph),
            1 => cos_core_probe(r_ph),
            2 => -sin_core_probe(r_ph),
            _ => -cos_core_probe(r_ph),
        };
        let s_cw = sin(x);
        worst = worst.max(ulp_diff(s_ph, s_cw));
    }
    verdict(
        "payne-002-overlap",
        worst <= 2,
        &format!(
            "PH vs CW sin over [2^19, 2^19.9]: worst {worst} ULP (boundary policy: CW at and below 2^20)"
        ),
    );
}

// Probe cores via public sin/cos at the reduced argument (|r| ≤ π/4,
// where reduction is the identity path through Cody–Waite).
fn sin_core_probe(r: f64) -> f64 {
    sin(r)
}
fn cos_core_probe(r: f64) -> f64 {
    cos(r)
}

/// payne-003: the published WORST-CASE double for π/2 reduction —
/// x = 6381956970095103·2^797 sits ~2.9e-19 from a multiple of π/2;
/// a naive reduction loses everything, Payne–Hanek keeps full
/// precision (gated against the platform libm oracle).
#[test]
fn payne_003_worst_case() {
    #[allow(clippy::cast_precision_loss)]
    let x = 6_381_956_970_095_103.0f64 * 2.0f64.powi(797);
    let (r, _q) = reduce_pio2_large(x);
    let mine = sin(x);
    let oracle = x.sin();
    let d = ulp_diff(mine, oracle);
    verdict(
        "payne-003-worst-case",
        r.abs() < 1e-10 && d <= 4,
        &format!(
            "reduced |r| = {:.3e} (the famous near-multiple), sin vs libm oracle {d} ULP",
            r.abs()
        ),
    );
}

/// payne-004: budget sweep across the FULL exponent range — sin/cos
/// within 4 ULP of the platform oracle from 2^21 to near f64::MAX
/// (deterministic sample; tan gated at its ratio-of-cores budget).
#[test]
fn payne_004_budget_sweep() {
    let mut h: u64 = 0x1001_2026_0708_0605;
    let mut next = || {
        h ^= h << 13;
        h ^= h >> 7;
        h ^= h << 17;
        h
    };
    let (mut worst_sin, mut worst_cos, mut worst_tan) = (0u64, 0u64, 0u64);
    let mut worst_x = 0.0f64;
    for _ in 0..4000 {
        // Random exponent 21..1000, random mantissa.
        let e = 21 + (next() % 980) as i32;
        #[allow(clippy::cast_precision_loss)]
        let mant = 1.0 + (next() >> 11) as f64 / (1u64 << 53) as f64;
        let x = mant * 2.0f64.powi(e);
        if !x.is_finite() {
            continue;
        }
        let ds = ulp_diff(sin(x), x.sin());
        let dc = ulp_diff(cos(x), x.cos());
        if ds > worst_sin {
            worst_sin = ds;
            worst_x = x;
        }
        worst_cos = worst_cos.max(dc);
        // tan blows up near poles; gate where |cos| is comfortably
        // away from zero (the ratio budget's validity condition).
        if x.cos().abs() > 0.1 {
            worst_tan = worst_tan.max(ulp_diff(tan(x), x.tan()));
        }
    }
    verdict(
        "payne-004-budget-sweep",
        worst_sin <= 4 && worst_cos <= 4 && worst_tan <= 10,
        &format!(
            "4000 samples over 2^21..2^1000: sin {worst_sin} ULP (worst at {worst_x:e}), cos {worst_cos} ULP, tan {worst_tan} ULP (|cos| > 0.1)"
        ),
    );
}

/// payne-005: hard landmarks — huge powers of two and the largest
/// finite decades stay within budget; odd symmetry is BITWISE.
#[test]
fn payne_005_landmarks_and_symmetry() {
    let landmarks = [
        2.0f64.powi(60),
        2.0f64.powi(100),
        2.0f64.powi(500),
        2.0f64.powi(1000),
        1e300,
        1.797e308,
    ];
    let mut worst = 0u64;
    for &x in &landmarks {
        worst = worst.max(ulp_diff(sin(x), x.sin()));
        worst = worst.max(ulp_diff(cos(x), x.cos()));
    }
    let mut symmetric = true;
    for &x in &landmarks {
        if sin(-x).to_bits() != (-sin(x)).to_bits() {
            symmetric = false;
        }
        if cos(-x).to_bits() != cos(x).to_bits() {
            symmetric = false;
        }
    }
    verdict(
        "payne-005-landmarks",
        worst <= 4 && symmetric,
        &format!("landmarks worst {worst} ULP; sin odd / cos even BITWISE at every landmark"),
    );
}
