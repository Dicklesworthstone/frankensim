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
//! - `sin`/`cos`: Cody–Waite THREE-PART π/2 reduction for |x| ≤ 2²⁰,
//!   Payne–Hanek 1280-bit reduction beyond (self-verifying Machin-generated
//!   limbs, see `payne`) — degree-13/12 Taylor cores on |r| ≤ π/4; budgets
//!   hold for ALL finite arguments.
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

/// Cody–Waite/Payne–Hanek dispatch boundary: |x| ≤ 2²⁰ uses the three-part
/// Cody–Waite reduction; beyond it the `payne` module's 1280-bit reduction
/// takes over (bead r6r5), so the trig ULP budgets now hold for ALL finite
/// arguments (measured 1 ULP across 2²¹..2¹⁰⁰⁰; the large-domain declared
/// ceiling is [`SIN_LARGE_ULP_BUDGET`]).
pub const TRIG_DOMAIN: f64 = 1_048_576.0; // 2^20

/// Declared ULP budget for sin/cos with |x| > [`TRIG_DOMAIN`] (Payne–Hanek
/// path; measured max 1 ULP over the exponent sweep, ceiling kept honest).
pub const SIN_LARGE_ULP_BUDGET: u64 = 4;

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
    let (r, quadrant) = if x.abs() > TRIG_DOMAIN {
        crate::payne::reduce_pio2_large(x)
    } else {
        reduce_pio2(x)
    };
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
    let (r, quadrant) = if x.abs() > TRIG_DOMAIN {
        crate::payne::reduce_pio2_large(x)
    } else {
        reduce_pio2(x)
    };
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

// ---------------------------------------------------------------------------
// Extension family (bead wf9.14): tan, atan, atan2, erf, erfc, pow.
// Additive only — nothing above this line changed (the original golden
// hash is untouched by construction).
// ---------------------------------------------------------------------------

/// Declared ULP budget for [`tan`] (|x| ≤ [`TRIG_DOMAIN`]; the ratio of
/// two ≤3-ULP cores on a SHARED reduced argument).
pub const TAN_ULP_BUDGET: u64 = 8;
/// Declared ULP budget for [`atan`]/[`atan2`].
pub const ATAN_ULP_BUDGET: u64 = 4;
/// Declared ULP budget for [`asin`]/[`acos`]: atan2's budget plus the
/// rounding of the factored complement √((1−x)(1+x)) (≤ 1.5 ULP into
/// the atan argument).
pub const ASIN_ULP_BUDGET: u64 = 6;
/// Declared ULP budget for [`erf`].
pub const ERF_ULP_BUDGET: u64 = 6;
/// Declared ULP budget for [`erfc`].
pub const ERFC_ULP_BUDGET: u64 = 10;
/// [`pow`] budget FORMULA (honest: the y·ln x magnification is real):
/// ≈ 3·(|y·ln x| + 1) + 5 ULP. The dd-ln refinement that removes the
/// magnification is recorded follow-up scope.
#[must_use]
pub fn pow_ulp_budget(x: f64, y: f64) -> u64 {
    let t = (y * ln(x.abs().max(f64::MIN_POSITIVE))).abs();
    (3.0 * (t + 1.0) + 5.0) as u64
}

/// Declared ULP budget for [`hypot`]: the fused `1 + r²` and the final
/// multiply each spend ≤ 1 ULP against libm's (near-correctly-rounded)
/// hypot; measured worst 2 ULP over 200k samples.
pub const HYPOT_ULP_BUDGET: u64 = 2;

/// √(x² + y²) — overflow/underflow-safe and DETERMINISTIC by construction:
/// max-scaling, a fused `1 + r²`, and the correctly-rounded hardware `sqrt`
/// are all IEEE-exact ops, so this is bit-identical on every conforming
/// target — UNLIKE `f64::hypot`, whose platform libm diverges across ISAs
/// (measured: an FNV hash over 250k results differs aarch64 vs x86-64). Route
/// any hypot that feeds solver state or cross-ISA golden bits through here.
/// Symmetric BITWISE (`hypot(x, y) == hypot(y, x)`). Per IEEE-754,
/// `hypot(±∞, y) = +∞` for EVERY `y` (including NaN); any other NaN input
/// yields NaN.
#[must_use]
pub fn hypot(x: f64, y: f64) -> f64 {
    // An infinite operand forces +∞ even against a NaN (IEEE-754), so it must
    // be handled before the magnitude ordering — a NaN would poison the `>=`.
    if x.is_infinite() || y.is_infinite() {
        return f64::INFINITY;
    }
    let (a, b) = (x.abs(), y.abs());
    let (hi, lo) = if a >= b { (a, b) } else { (b, a) };
    if hi == 0.0 {
        return 0.0;
    }
    let r = lo / hi;
    hi * sqrt(r.mul_add(r, 1.0))
}

/// tan(x), deterministic strict mode; budget valid for |x| ≤
/// [`TRIG_DOMAIN`]. Odd BITWISE by construction (shared symmetric
/// reduction + odd/even cores). Pole neighborhoods return the honest
/// ratio (huge but finite until the reduced argument actually hits a
/// core zero).
#[must_use]
pub fn tan(x: f64) -> f64 {
    if x.is_nan() || x.is_infinite() {
        return f64::NAN;
    }
    let (r, quadrant) = if x.abs() > TRIG_DOMAIN {
        crate::payne::reduce_pio2_large(x)
    } else {
        reduce_pio2(x)
    };
    let s = sin_core(r);
    let c = cos_core(r);
    if quadrant & 1 == 0 { s / c } else { -(c / s) }
}

/// Odd Taylor core for atan on |t| ≤ tan(π/12) ≈ 0.2679: 16 terms of
/// t·Σ (−1)ᵏ·z^k/(2k+1); the k = 16 tail is ≈ 0.2679³³/33 ≈ 2.5e-21.
fn atan_core(t: f64) -> f64 {
    let z = t * t;
    let mut p = -1.0 / 31.0;
    // Descending Horner over the exact rational coefficients.
    let coeffs = [
        29.0, 27.0, 25.0, 23.0, 21.0, 19.0, 17.0, 15.0, 13.0, 11.0, 9.0, 7.0, 5.0, 3.0,
    ];
    let mut sign = 1.0;
    for &d in &coeffs {
        p = z.mul_add(p, sign / d);
        sign = -sign;
    }
    let p = z.mul_add(p, 1.0); // leading (k = 0) term: +1

    t * p
}

// π/2 and π/6 as hi/lo pairs (hi = nearest f64; lo = the residual, well
// below hi's ULP): summing (hi − small) + lo keeps the constant's full
// accuracy in the combined result.
const PI_2_HI: f64 = std::f64::consts::FRAC_PI_2;
const PI_2_LO: f64 = 6.123_233_995_736_766e-17;
const PI_6_HI: f64 = std::f64::consts::FRAC_PI_6;
const PI_6_LO: f64 = -1.643_486_876_741_777e-17;
/// tan(π/12) rounded UP (the reduction threshold; overlap is harmless).
const TAN_PI_12: f64 = 0.267_949_192_431_123_05;
const SQRT_3: f64 = 1.732_050_807_568_877_2;
const SQRT_3_LO: f64 = 1.020_508_405_819_248_3e-16;

/// atan(x), deterministic strict mode. Odd BITWISE by construction
/// (sign folded at entry). atan(±∞) = ±π/2.
#[must_use]
pub fn atan(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    let sign = if x.is_sign_negative() { -1.0 } else { 1.0 };
    let a = x.abs();
    if a.is_infinite() {
        return sign * (PI_2_HI + PI_2_LO);
    }
    // Range reduction: a > 1 → π/2 − atan(1/a).
    let (invert, t0) = if a > 1.0 { (true, 1.0 / a) } else { (false, a) };
    // Second reduction to |t| ≤ tan(π/12) via the π/6 rotation:
    // atan(t) = π/6 + atan((√3·t − 1)/(√3 + t)).
    let (rotated, t) = if t0 > TAN_PI_12 {
        let num = SQRT_3.mul_add(t0, -1.0) + SQRT_3_LO * t0;
        (true, num / (SQRT_3 + t0))
    } else {
        (false, t0)
    };
    let base = atan_core(t);
    let mut v = if rotated {
        (PI_6_HI + base) + PI_6_LO
    } else {
        base
    };
    if invert {
        v = (PI_2_HI - v) + PI_2_LO;
    }
    sign * v
}

/// atan2(y, x) with IEEE special-case conventions, deterministic strict
/// mode. Sign symmetries hold bitwise (sign of y folded at entry).
#[must_use]
pub fn atan2(y: f64, x: f64) -> f64 {
    const PI_HI: f64 = std::f64::consts::PI;
    const PI_LO: f64 = 1.224_646_799_147_353_2e-16;
    if x.is_nan() || y.is_nan() {
        return f64::NAN;
    }
    let sign = if y.is_sign_negative() { -1.0 } else { 1.0 };
    let ay = y.abs();
    // Special cases (IEEE 754 / libm conventions).
    if ay == 0.0 {
        return if x.is_sign_negative() {
            sign * (PI_HI + PI_LO)
        } else {
            sign * 0.0
        };
    }
    if x == 0.0 {
        return sign * (PI_2_HI + PI_2_LO);
    }
    if x.is_infinite() {
        return match (x > 0.0, ay.is_infinite()) {
            (true, true) => sign * (PI_2_HI / 2.0 + PI_2_LO / 2.0), // π/4
            (true, false) => sign * 0.0,
            (false, true) => sign * 3.0f64.mul_add(PI_2_HI / 2.0, 3.0 * (PI_2_LO / 2.0)), // 3π/4
            (false, false) => sign * (PI_HI + PI_LO),
        };
    }
    if ay.is_infinite() {
        return sign * (PI_2_HI + PI_2_LO);
    }
    let base = atan(ay / x.abs());
    if x > 0.0 {
        sign * base
    } else {
        sign * ((PI_HI - base) + PI_LO)
    }
}

/// asin(x), deterministic strict mode, via `atan2(x, √((1−x)(1+x)))`.
/// The FACTORED product keeps the complement conditioned at the
/// endpoints (1 − x² cancels catastrophically for |x| → 1; the factors
/// do not). Odd BITWISE (inherited from atan2's sign fold). Domain:
/// |x| ≤ 1; outside → NaN (libm convention). asin(±1) = ±π/2 exactly
/// (nearest f64), through atan2's x = 0 special case.
#[must_use]
pub fn asin(x: f64) -> f64 {
    if x.is_nan() || x.abs() > 1.0 {
        return f64::NAN;
    }
    atan2(x, ((1.0 - x) * (1.0 + x)).sqrt())
}

/// acos(x), deterministic strict mode, via `atan2(√((1−x)(1+x)), x)`.
/// Same factored complement as [`asin`]. Domain: |x| ≤ 1; outside →
/// NaN. acos(1) = +0 and acos(−1) = π exactly (nearest f64), through
/// atan2's y = 0 special case; the reflection acos(−x) = π − acos(x)
/// holds within the declared budget (not bitwise — π rounds).
#[must_use]
pub fn acos(x: f64) -> f64 {
    if x.is_nan() || x.abs() > 1.0 {
        return f64::NAN;
    }
    atan2(((1.0 - x) * (1.0 + x)).sqrt(), x)
}

/// 2/√π and π as double-double values, derived at runtime from exact
/// integer-free dd arithmetic (deterministic; avoids hand-transcribed
/// low limbs).
fn two_over_sqrt_pi_dd() -> crate::dd::Dd {
    use crate::dd::Dd;
    let pi = Dd::from_pair(std::f64::consts::PI, 1.224_646_799_147_353_2e-16);
    Dd::from_f64(2.0) / pi.sqrt()
}

/// erf(x), deterministic strict mode. Odd BITWISE by construction. The
/// Taylor sum runs in DOUBLE-DOUBLE (the reason dd lives at L0): at
/// x = 3 the alternating series cancels ~4 digits, which dd absorbs.
#[must_use]
pub fn erf(x: f64) -> f64 {
    use crate::dd::Dd;
    if x.is_nan() {
        return f64::NAN;
    }
    let sign = if x.is_sign_negative() { -1.0 } else { 1.0 };
    let a = x.abs();
    if a >= 6.0 {
        return sign; // erf(6) − 1 ≈ −2e-17: rounds to ±1
    }
    if a <= 3.0 {
        sign * erf_dd_small(a).to_f64()
    } else {
        // 1 − erfc via the continued fraction (erfc ≤ 2.2e-5 here, so the
        // subtraction from 1 loses nothing).
        sign * (Dd::from_f64(1.0) - erfc_dd_large(a)).to_f64()
    }
}

/// erfc(x), deterministic strict mode; full-precision even in the far
/// tail (exact-dd x², Laplace continued fraction).
#[must_use]
pub fn erfc(x: f64) -> f64 {
    use crate::dd::Dd;
    if x.is_nan() {
        return f64::NAN;
    }
    if x == 0.0 {
        return 1.0;
    }
    if x < 0.0 {
        // erfc(−a) = 2 − erfc(a) = 1 + erf(a): no cancellation.
        let a = -x;
        return if a >= 6.0 {
            2.0
        } else if a <= 3.0 {
            (Dd::from_f64(1.0) + erf_dd_small(a)).to_f64()
        } else {
            (Dd::from_f64(2.0) - erfc_dd_large(a)).to_f64()
        };
    }
    if x <= 3.0 {
        // 1 − erf in dd: the cancellation near x = 3 (erfc ~ 2e-5) is
        // exactly why the dd sum exists.
        (crate::dd::Dd::from_f64(1.0) - erf_dd_small(x)).to_f64()
    } else if x >= 27.5 {
        0.0 // exp(−x²) underflows past every subnormal
    } else {
        erfc_dd_large(x).to_f64()
    }
}

/// BOTH erf paths at the same x (Taylor-dd and CF-dd) — exposed for the
/// conformance battery's path-consistency check: the two constructions
/// are algorithmically disjoint, so their agreement at shared x is
/// genuine cross-validation (the external oracle is weaker than the
/// implementation in the cancellation band).
#[doc(hidden)]
#[must_use]
pub fn erf_both_paths(x: f64) -> (f64, f64) {
    use crate::dd::Dd;
    let taylor = erf_dd_small(x).to_f64();
    let cf = (Dd::from_f64(1.0) - erfc_dd_large(x)).to_f64();
    (taylor, cf)
}

/// dd Taylor sum of erf on [0, 3].
fn erf_dd_small(a: f64) -> crate::dd::Dd {
    use crate::dd::Dd;
    let x = Dd::from_f64(a);
    let x2 = x * x;
    let mut term = x; // x^(2k+1)/k! at k = 0
    let mut sum = x; // k = 0 contribution: x/1
    for k in 1..=90i32 {
        term = term * x2 / Dd::from_f64(f64::from(k));
        let contrib = term / Dd::from_f64(2.0f64.mul_add(f64::from(k), 1.0));
        sum = if k % 2 == 1 {
            sum - contrib
        } else {
            sum + contrib
        };
        if contrib.hi.abs() < 1e-35 * sum.hi.abs() {
            break;
        }
    }
    sum * two_over_sqrt_pi_dd()
}

/// dd Laplace continued fraction for erfc on x > 3:
/// erfc(x) = exp(−x²)/√π · 1/(x + (1/2)/(x + 1/(x + (3/2)/(x + …)))).
fn erfc_dd_large(x: f64) -> crate::dd::Dd {
    use crate::dd::Dd;
    let xd = Dd::from_f64(x);
    let mut cf = xd; // innermost tail ≈ x
    for k in (1..=60i32).rev() {
        cf = xd + Dd::from_f64(f64::from(k) / 2.0) / cf;
    }
    // exp(−x²) with the EXACT dd square (f64 x² alone would cost ~x²/2
    // ULP of relative error — hundreds at x = 25).
    let x2 = xd * xd;
    let e = Dd::from_f64(exp(-x2.hi)) * (Dd::from_f64(1.0) - Dd::from_f64(x2.lo));
    let pi = Dd::from_pair(std::f64::consts::PI, 1.224_646_799_147_353_2e-16);
    e / pi.sqrt() / cf
}

/// pow(x, y), deterministic strict mode, IEEE special cases. Budget is
/// the HONEST formula [`pow_ulp_budget`] (the |y·ln x| magnification is
/// intrinsic to the exp∘ln route; dd-ln refinement recorded). Fast
/// paths: integer y (repeated squaring), y = ±0.5 (exact sqrt).
#[must_use]
// Exact equality against special VALUES (1.0, ±0.5, integers) is the
// IEEE-mandated special-case detection, not tolerance arithmetic.
#[allow(clippy::float_cmp)]
pub fn pow(x: f64, y: f64) -> f64 {
    // IEEE special cases first.
    if y == 0.0 {
        return 1.0;
    }
    if x == 1.0 {
        return 1.0; // pow(1, y) = 1 for every y, including NaN (IEEE-754 §9.2.1)
    }
    if x.is_nan() || y.is_nan() {
        return f64::NAN;
    }
    let y_int = y.fract() == 0.0 && y.abs() < 9.0e15;
    let y_odd = y_int && (y / 2.0).fract() != 0.0;
    if x == 0.0 {
        return match (y > 0.0, y_odd && x.is_sign_negative()) {
            (true, true) => -0.0,
            (true, false) => 0.0,
            (false, true) => f64::NEG_INFINITY,
            (false, false) => f64::INFINITY,
        };
    }
    if x.is_infinite() || y.is_infinite() {
        // Delegate the (finite) magnitude comparison logic.
        let ax = x.abs();
        if y.is_infinite() && ax == 1.0 {
            return 1.0; // pow(±1, ±∞) = 1 (IEEE-754); x = +1 handled above, this is x = -1
        }
        let big = if y.is_infinite() {
            if (ax > 1.0) == (y > 0.0) {
                f64::INFINITY
            } else {
                0.0
            }
        } else if y > 0.0 {
            f64::INFINITY
        } else {
            0.0
        };
        let neg = x.is_sign_negative() && y_odd;
        return if neg { -big } else { big };
    }
    if x < 0.0 {
        if !y_int {
            return f64::NAN; // complex result: refused per IEEE
        }
        let mag = pow(-x, y);
        return if y_odd { -mag } else { mag };
    }
    // Integer fast path: repeated squaring (error ~ log2|y| ULP).
    if y_int && y.abs() <= 512.0 {
        let mut base = if y > 0.0 { x } else { 1.0 / x };
        let mut n = y.abs() as u64;
        let mut acc = 1.0f64;
        while n > 0 {
            if n & 1 == 1 {
                acc *= base;
            }
            base *= base;
            n >>= 1;
        }
        return acc;
    }
    if y == 0.5 {
        return sqrt(x);
    }
    if y == -0.5 {
        return 1.0 / sqrt(x);
    }
    // General path: exp(y·ln x) with the product carried in dd.
    let lx = ln(x);
    let (p_hi, p_lo) = crate::eft::two_prod(y, lx);
    // exp(hi + lo) = exp(hi)·exp(lo) ≈ exp(hi)·(1 + lo): lo ≤ ½ulp(hi),
    // so the first-order correction is exact to ~2⁻¹⁰⁶.
    exp(p_hi) * (1.0 + p_lo)
}

/// xⁿ for integer n, deterministic strict mode.
///
/// `f64::powi` must never feed golden bits: the `llvm.powi` intrinsic
/// has no pinned operation order, so its rounding sequence is
/// optimization-level-dependent (observed 1-ULP debug/release
/// divergence from n = 4 upward — bead 4xnt). This function pins the
/// order in source: LSB-first binary exponentiation (`r *= a` on each
/// set bit, `a *= a` per halving). Negative exponents use one final
/// reciprocal for ordinary results. If the positive intermediate
/// overflows and that reciprocal collapses to zero, the same pinned tree
/// is re-evaluated from the reciprocal base to recover any representable
/// subnormal result such as 2^-1024.
///
/// Semantics: n = 0 returns 1.0 for every x (including NaN and ∞,
/// matching `f64::powi`); n = i32::MIN is handled (|n| carried in
/// unsigned arithmetic). For positive n ≤ 512 this agrees bitwise with
/// [`pow`]'s integer fast path (same order).
///
/// Accuracy: one rounding per executed multiply (popcount(n) +
/// ⌊log₂ n⌋ − 1 of them), NOT correctly rounded; measured ≤ 2 ULP vs
/// platform `powi` for |n| ≤ 64 in the extension battery.
#[must_use]
pub fn powi(x: f64, n: i32) -> f64 {
    fn unsigned(mut a: f64, mut b: u64) -> f64 {
        let mut r = 1.0f64;
        loop {
            if b & 1 == 1 {
                r *= a;
            }
            b /= 2;
            if b == 0 {
                break;
            }
            a *= a;
        }
        r
    }

    let b = i64::from(n).unsigned_abs();
    let positive = unsigned(x, b);
    if n >= 0 {
        return positive;
    }
    let reciprocal = 1.0 / positive;
    if reciprocal == 0.0 && x.is_finite() && x != 0.0 {
        return unsigned(1.0 / x, b);
    }
    reciprocal
}
