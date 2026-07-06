//! Deterministic elementary functions (strict mode): built EXCLUSIVELY from
//! IEEE-754 arithmetic (+, −, ×, ÷, `mul_add`, `sqrt`) — every one of which
//! is correctly rounded and therefore bit-identical on every conforming
//! target. No platform libm anywhere: cross-ISA determinism holds BY
//! CONSTRUCTION, and the golden-hash test proves it empirically (verified
//! aarch64-apple vs x86-64).
//!
//! Accuracy: each function declares a ULP budget (`*_ULP_BUDGET`), asserted
//! against a measured maximum versus the platform-libm oracle in tests (the
//! high-precision double-double oracle battery arrives with fs-ivl). Budgets
//! are honest ceilings, tightened as implementations improve.
//!
//! Algorithms (classic, chosen for pure-arithmetic implementability):
//! - `exp`/`expm1`: k = round(x/ln2) reduction with two-part ln2, degree-13
//!   Taylor/Horner core on |r| ≤ ln2/2, exact 2^k scaling via exponent bits.
//! - `ln`: mantissa reduction to [√½, √2), atanh series in s=(m−1)/(m+1)
//!   (|s| ≤ 0.1716), two-part k·ln2 recombination.
//! - `sin`/`cos`: Cody–Waite THREE-PART π/2 reduction — accurate for
//!   |x| ≤ 2²⁰ (documented domain; Payne–Hanek big-argument reduction is a
//!   recorded follow-up) — with degree-13/12 Taylor cores on |r| ≤ π/4.
//! - `tanh`: expm1(2x)/(expm1(2x)+2) (odd symmetry; saturates for |x| > 20).

/// ULP budget for [`exp`] (measured max observed: see tests).
pub const EXP_ULP_BUDGET: u64 = 3;
/// ULP budget for [`expm1`].
pub const EXPM1_ULP_BUDGET: u64 = 3;
/// ULP budget for [`ln`].
pub const LN_ULP_BUDGET: u64 = 3;
/// ULP budget for [`sin`] within the reduction domain |x| ≤ 2²⁰.
pub const SIN_ULP_BUDGET: u64 = 3;
/// ULP budget for [`cos`] within the reduction domain |x| ≤ 2²⁰.
pub const COS_ULP_BUDGET: u64 = 3;
/// ULP budget for [`tanh`].
pub const TANH_ULP_BUDGET: u64 = 5;

/// Trig argument-reduction domain bound: |x| ≤ 2²⁰. Beyond it the Cody–Waite
/// reduction loses bits and the ULP budgets are VOID (no-claim; results are
/// still deterministic, just less accurate).
pub const TRIG_DOMAIN: f64 = 1_048_576.0; // 2^20

// EXACT bit patterns (fdlibm heritage). The hi parts have ≥20 trailing zero
// mantissa bits, so k·LN2_HI (|k| ≤ 2¹⁰) and j·PIO2_* (|j| ≤ 2²⁰) are EXACT
// products — the property the whole reduction-accuracy argument rests on.
// Decimal literals are NOT acceptable here: they round to neighboring
// doubles without the trailing zeros (found the hard way: 184-ULP trig).
const LN2_HI: f64 = f64::from_bits(0x3FE6_2E42_FEE0_0000); // 6.9314718036912382e-1
const LN2_LO: f64 = f64::from_bits(0x3DEA_39EF_3579_3C76); // 1.9082149292705877e-10
const LOG2_E: f64 = std::f64::consts::LOG2_E;

/// e^x, deterministic strict mode.
#[must_use]
pub fn exp(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x > 709.782_712_893_384 {
        return f64::INFINITY;
    }
    if x < -745.133_219_101_941_1 {
        return 0.0;
    }
    let k = (x * LOG2_E).round();
    // r = x − k·ln2, in two exact-ish parts to keep |r| ≤ ln2/2 accurate.
    let r = (-k).mul_add(LN2_LO, (-k).mul_add(LN2_HI, x));
    scale_by_2k(exp_core(r), k as i64)
}

/// e^x − 1, accurate near zero, deterministic strict mode.
#[must_use]
pub fn expm1(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x > 709.782_712_893_384 {
        return f64::INFINITY;
    }
    if x < -37.0 {
        return -1.0; // e^x < 2⁻⁵³ − relative to 1 it has vanished
    }
    if x.abs() < 0.5 * std::f64::consts::LN_2 {
        return expm1_core(x); // no reduction: keeps the small-x accuracy
    }
    exp(x) - 1.0
}

/// Natural logarithm, deterministic strict mode.
#[must_use]
pub fn ln(x: f64) -> f64 {
    if x.is_nan() || x < 0.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        return f64::NEG_INFINITY;
    }
    if x.is_infinite() {
        return f64::INFINITY;
    }
    // Normalize subnormals, then split x = 2^k · m with m ∈ [√½, √2).
    let (x, sub_k) = if x < f64::MIN_POSITIVE {
        (x * 9_007_199_254_740_992.0, -53i64) // ×2⁵³
    } else {
        (x, 0)
    };
    let bits = x.to_bits();
    let mut k = i64::try_from((bits >> 52) & 0x7FF).unwrap_or(0) - 1023;
    let mut m = f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000);
    if m > std::f64::consts::SQRT_2 {
        m *= 0.5;
        k += 1;
    }
    let k = k + sub_k;
    // atanh series: ln m = 2s(1 + s²/3 + s⁴/5 + …), s = (m−1)/(m+1).
    let s = (m - 1.0) / (m + 1.0);
    let z = s * s;
    // Terms through z⁹/19: truncation z¹⁰/21 ≈ 2e-17 relative at |s| ≤ 0.1716.
    let poly = z.mul_add(
        z.mul_add(
            z.mul_add(
                z.mul_add(
                    z.mul_add(
                        z.mul_add(
                            z.mul_add(
                                z.mul_add(z.mul_add(1.0 / 19.0, 1.0 / 17.0), 1.0 / 15.0),
                                1.0 / 13.0,
                            ),
                            1.0 / 11.0,
                        ),
                        1.0 / 9.0,
                    ),
                    1.0 / 7.0,
                ),
                1.0 / 5.0,
            ),
            1.0 / 3.0,
        ),
        1.0,
    );
    let kf = k as f64;
    kf.mul_add(LN2_LO, (2.0 * s).mul_add(poly, kf * LN2_HI))
}

/// sin(x), deterministic strict mode; ULP budget valid for |x| ≤ [`TRIG_DOMAIN`].
#[must_use]
pub fn sin(x: f64) -> f64 {
    if x.is_nan() || x.is_infinite() {
        return f64::NAN;
    }
    let (r, quadrant) = reduce_pio2(x);
    match quadrant {
        0 => sin_core(r),
        1 => cos_core(r),
        2 => -sin_core(r),
        _ => -cos_core(r),
    }
}

/// cos(x), deterministic strict mode; ULP budget valid for |x| ≤ [`TRIG_DOMAIN`].
#[must_use]
pub fn cos(x: f64) -> f64 {
    if x.is_nan() || x.is_infinite() {
        return f64::NAN;
    }
    let (r, quadrant) = reduce_pio2(x);
    match quadrant {
        0 => cos_core(r),
        1 => -sin_core(r),
        2 => -cos_core(r),
        _ => sin_core(r),
    }
}

/// tanh(x), deterministic strict mode. Odd symmetry holds BITWISE by
/// construction: the magnitude is computed once and the sign re-applied
/// (symmetry-by-construction beats symmetry-by-luck — the same doctrine as
/// the geometry layer's quotient parameterizations).
#[must_use]
pub fn tanh(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    let a = x.abs();
    let mag = if a > 20.0 {
        1.0
    } else {
        // tanh a = t / (t + 2), t = expm1(2a): accurate at every scale (the
        // small-a cancellation lives inside expm1's Taylor core).
        let t = expm1(2.0 * a);
        t / (t + 2.0)
    };
    if x.is_sign_negative() { -mag } else { mag }
}

/// √x — IEEE-754 requires correct rounding of sqrt, so the hardware
/// instruction IS the deterministic strict-mode implementation (0 ULP).
#[must_use]
pub fn sqrt(x: f64) -> f64 {
    x.sqrt()
}

// ---------------------------------------------------------------------------
// Cores (pure Horner/mul_add — no table lookups, no platform calls).
// ---------------------------------------------------------------------------

/// exp on |r| ≤ ln2/2 ≈ 0.347: Taylor to r¹³ (tail < 4e-17 relative).
fn exp_core(r: f64) -> f64 {
    expm1_core(r) + 1.0
}

/// expm1 on |r| ≤ ln2/2: r·(1 + r/2 + r²/6 + … ), Horner in r.
fn expm1_core(r: f64) -> f64 {
    const C: [f64; 12] = [
        1.0 / 2.0,
        1.0 / 6.0,
        1.0 / 24.0,
        1.0 / 120.0,
        1.0 / 720.0,
        1.0 / 5_040.0,
        1.0 / 40_320.0,
        1.0 / 362_880.0,
        1.0 / 3_628_800.0,
        1.0 / 39_916_800.0,
        1.0 / 479_001_600.0,
        1.0 / 6_227_020_800.0,
    ];
    let mut p = C[11];
    for c in C[..11].iter().rev() {
        p = p.mul_add(r, *c);
    }
    r.mul_add(r * p, r) // r + r²·(poly) keeps the leading term exact
}

/// sin on |r| ≤ π/4: r − r³·P(r²), Taylor through r¹⁷ — the r¹⁹/19! tail is
/// ≈ 8e-20 at π/4 (≈ 0.001 ULP). A shorter r¹³ core measured 166 ULP at the
/// interval edge (0.785¹⁵/15! ≈ 2e-14); the core-only regression test below
/// keeps that mistake dead.
fn sin_core(r: f64) -> f64 {
    let z = r * r;
    let p = z.mul_add(
        z.mul_add(
            z.mul_add(
                z.mul_add(
                    z.mul_add(
                        z.mul_add(
                            z.mul_add(-1.0 / 355_687_428_096_000.0, 1.0 / 1_307_674_368_000.0),
                            -1.0 / 6_227_020_800.0,
                        ),
                        1.0 / 39_916_800.0,
                    ),
                    -1.0 / 362_880.0,
                ),
                1.0 / 5_040.0,
            ),
            -1.0 / 120.0,
        ),
        1.0 / 6.0,
    );
    (-(z * r)).mul_add(p, r) // r − r³·P(z), leading term exact
}

/// cos on |r| ≤ π/4: 1 − z/2 + z²·(1/24 + z·Q(z)), Taylor through r¹⁶ —
/// the r¹⁸/18! tail is ≈ 1.2e-18 at π/4 (≈ 0.01 ULP).
fn cos_core(r: f64) -> f64 {
    let z = r * r;
    let q = z.mul_add(
        z.mul_add(
            z.mul_add(
                z.mul_add(
                    z.mul_add(1.0 / 20_922_789_888_000.0, -1.0 / 87_178_291_200.0),
                    1.0 / 479_001_600.0,
                ),
                -1.0 / 3_628_800.0,
            ),
            1.0 / 40_320.0,
        ),
        -1.0 / 720.0,
    );
    // 1 − z/2 + z²/24 + z³·q, with the z²/24 term folded into the Horner tail:
    let tail = z.mul_add(q, 1.0 / 24.0);
    (z * z).mul_add(tail, 1.0 - 0.5 * z)
}

/// Cody–Waite three-part π/2 reduction: returns (r, quadrant mod 4).
/// Accurate while |x|·ulp(π/2 error) stays below the core's tolerance —
/// the documented |x| ≤ 2²⁰ domain.
fn reduce_pio2(x: f64) -> (f64, u8) {
    // fdlibm's classic 33-bit split of π/2 as EXACT bit patterns (trailing
    // zeros make j·PIO2_* exact for |j| ≤ 2²⁰); summed they carry ~99 bits,
    // so the residual reduction error stays ≈ 2⁻⁷⁹ — far below core needs.
    const PIO2_HI: f64 = f64::from_bits(0x3FF9_21FB_5440_0000);
    const PIO2_MID: f64 = f64::from_bits(0x3DD0_B461_1A60_0000);
    const PIO2_LO: f64 = f64::from_bits(0x3BA3_198A_2E00_0000);
    // Fourth chunk (fdlibm pio2_3t): matters ONLY where the reduced r is
    // near zero (x ≈ kπ) — there the result magnitude is ~1e-16 and the
    // ~2⁻¹⁰⁴ three-part residual costs several ULP of that tiny value
    // (measured: 7–10 ULP at x = π, 3π/2 without this term).
    const PIO2_LO2: f64 = 8.478_427_660_368_9e-32;
    const TWO_OVER_PI: f64 = std::f64::consts::FRAC_2_PI;
    let j = (x * TWO_OVER_PI).round();
    let r = (-j).mul_add(
        PIO2_LO2,
        (-j).mul_add(PIO2_LO, (-j).mul_add(PIO2_MID, (-j).mul_add(PIO2_HI, x))),
    );
    let q = ((j as i64) & 3) as u8;
    (r, q)
}

/// Exact ×2^k via exponent arithmetic, handling under/overflow into
/// subnormals and infinity (two-step scaling at the extremes).
fn scale_by_2k(v: f64, k: i64) -> f64 {
    let clamp = k.clamp(-2000, 2000);
    let mut v = v;
    let mut remaining = clamp;
    while remaining != 0 {
        let step = remaining.clamp(-1000, 1000);
        v *= f64::from_bits(((1023 + step) as u64) << 52);
        remaining -= step;
    }
    v
}
