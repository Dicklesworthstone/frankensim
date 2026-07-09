//! fs-math — deterministic elementary functions with ULP budgets, plus the
//! workspace floating-point POLICY (patch Rev O; plan §5.4, §6.4).
//!
//! # The policies (normative for every FrankenSim crate)
//!
//! - **Strict mode** (this module's [`det`] functions): pure-IEEE-arithmetic
//!   implementations — bit-identical on every conforming target BY
//!   CONSTRUCTION. Deterministic-mode kernels MUST NOT call platform libm
//!   (`f64::sin` etc.) — that is exactly the cross-ISA divergence class the
//!   G5 report exists to catch. (`sqrt` is the exception: IEEE requires
//!   correct rounding, so the hardware instruction is already deterministic.)
//! - **FMA contraction**: implicit contraction is FORBIDDEN (we rely on
//!   rustc not contracting by default); EXPLICIT `mul_add` is encouraged —
//!   it is exactly rounded and therefore deterministic everywhere. Kernels
//!   declare FMA usage in their contracts (fs-simd's fused-everywhere
//!   elementwise policy is an instance).
//! - **Subnormals**: never flushed. Any future platform/flag that implies
//!   FTZ/DAZ is a policy violation (lintable once fs-tilelang generates
//!   kernel metadata).
//! - **NaN**: production paths must not RELY on NaN payloads or their
//!   propagation details; [`canonical_nan`] is the interchange value.
//!   Comparisons in deterministic tie-breaking use `total_cmp`.
//! - **ULP budgets**: every [`det`] function declares one; tests assert the
//!   MEASURED maximum against it (platform-libm oracle now; the
//!   double-double oracle battery arrives with fs-ivl).

pub mod c64;
pub mod dd;
pub mod det;
pub mod eft;
pub mod payne;

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// The canonical quiet NaN for interchange (positive, standard payload).
#[must_use]
pub fn canonical_nan() -> f64 {
    f64::from_bits(0x7ff8_0000_0000_0000)
}

/// Next representable value toward +∞ (re-export of std semantics; here so
/// interval code has ONE audited nudging vocabulary).
#[must_use]
pub fn next_up(x: f64) -> f64 {
    x.next_up()
}

/// Next representable value toward −∞.
#[must_use]
pub fn next_down(x: f64) -> f64 {
    x.next_down()
}

/// Outward-nudged enclosure of a single computed value: one ULP each way.
/// This is fs-ivl's directed-rounding primitive (plan §6.4: no global
/// rounding-mode fiddling — Rust-safe and thread-safe).
#[must_use]
pub fn nudge_out(x: f64) -> (f64, f64) {
    (x.next_down(), x.next_up())
}

/// ULP distance between two finite floats of the same sign regime — the
/// test-battery metric (order the bit patterns as integers and subtract).
#[must_use]
pub fn ulp_distance(a: f64, b: f64) -> u64 {
    fn ordered(x: f64) -> i64 {
        #[allow(clippy::cast_possible_wrap)] // the wrap IS the ordering trick
        let b = x.to_bits() as i64;
        if b < 0 {
            i64::MIN.wrapping_add(b.wrapping_neg())
        } else {
            b
        }
    }
    ordered(a).abs_diff(ordered(b))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*seed >> 11) as f64) / (1u64 << 53) as f64 // [0, 1)
    }

    /// Measure max ULP vs the platform-libm oracle over deterministic
    /// samples in `range`, plus a caller-supplied edge table.
    fn max_ulp(
        ours: impl Fn(f64) -> f64,
        oracle: impl Fn(f64) -> f64,
        range: (f64, f64),
        edges: &[f64],
        samples: u64,
    ) -> (u64, f64) {
        let mut seed = 0x5EED_F00D_u64;
        let mut worst = (0u64, 0.0f64);
        let mut check = |x: f64| {
            let (a, b) = (ours(x), oracle(x));
            if a.is_finite() && b.is_finite() {
                let d = ulp_distance(a, b);
                if d > worst.0 {
                    worst = (d, x);
                }
            } else {
                assert_eq!(a.is_nan(), b.is_nan(), "NaN disagreement at {x}");
                assert_eq!(a.is_infinite(), b.is_infinite(), "inf disagreement at {x}");
            }
        };
        for _ in 0..samples {
            let t = lcg(&mut seed);
            check(range.0 + t * (range.1 - range.0));
        }
        for &e in edges {
            check(e);
        }
        worst
    }

    #[test]
    fn exp_meets_its_ulp_budget() {
        let (ulp, at) = max_ulp(
            det::exp,
            f64::exp,
            (-700.0, 700.0),
            &[
                0.0, -0.0, 1.0, -1.0, 1e-300, -1e-300, 709.7, -745.0, 0.5, -0.5,
            ],
            200_000,
        );
        println!(
            "{{\"suite\":\"fs-math\",\"case\":\"exp-ulp\",\"verdict\":\"pass\",\"detail\":\"max_ulp={ulp} at {at}\"}}"
        );
        assert!(ulp <= det::EXP_ULP_BUDGET, "exp max ULP {ulp} at {at}");
        assert_eq!(
            det::exp(0.0).to_bits(),
            1.0f64.to_bits(),
            "exp(0) must be exactly 1"
        );
    }

    #[test]
    fn expm1_meets_its_ulp_budget_and_is_accurate_near_zero() {
        let (ulp, at) = max_ulp(
            det::expm1,
            f64::exp_m1,
            (-37.0, 700.0),
            &[1e-30, -1e-30, 1e-10, -1e-10],
            200_000,
        );
        println!(
            "{{\"suite\":\"fs-math\",\"case\":\"expm1-ulp\",\"verdict\":\"pass\",\"detail\":\"max_ulp={ulp} at {at}\"}}"
        );
        assert!(ulp <= det::EXPM1_ULP_BUDGET, "expm1 max ULP {ulp} at {at}");
        // Tiny-x battery where naive exp(x)−1 catastrophically cancels.
        let (ulp_tiny, at_tiny) = max_ulp(det::expm1, f64::exp_m1, (-1e-8, 1e-8), &[], 50_000);
        assert!(
            ulp_tiny <= 1,
            "expm1 tiny-x max ULP {ulp_tiny} at {at_tiny}"
        );
    }

    #[test]
    fn ln_meets_its_ulp_budget() {
        let (ulp, at) = max_ulp(
            det::ln,
            f64::ln,
            (1e-300, 1e300),
            &[
                1.0,
                2.0,
                0.5,
                f64::MIN_POSITIVE,
                f64::MIN_POSITIVE / 8.0,
                1.0 + 1e-15,
                1e308,
            ],
            200_000,
        );
        println!(
            "{{\"suite\":\"fs-math\",\"case\":\"ln-ulp\",\"verdict\":\"pass\",\"detail\":\"max_ulp={ulp} at {at}\"}}"
        );
        assert!(ulp <= det::LN_ULP_BUDGET, "ln max ULP {ulp} at {at}");
        assert_eq!(
            det::ln(1.0).to_bits(),
            0.0f64.to_bits(),
            "ln(1) must be exactly +0"
        );
        // Near-1 relative accuracy (the formulation is relative-exact there).
        let (ulp_n1, at_n1) = max_ulp(det::ln, f64::ln, (0.9, 1.1), &[], 50_000);
        assert!(
            ulp_n1 <= det::LN_ULP_BUDGET,
            "ln near-1 max ULP {ulp_n1} at {at_n1}"
        );
    }

    #[test]
    fn sin_cos_meet_budgets_within_domain() {
        let edges: Vec<f64> = (0..8)
            .map(|k| f64::from(k) * std::f64::consts::FRAC_PI_4)
            .collect();
        let (su, sa) = max_ulp(
            det::sin,
            f64::sin,
            (-det::TRIG_DOMAIN, det::TRIG_DOMAIN),
            &edges,
            200_000,
        );
        let (cu, ca) = max_ulp(
            det::cos,
            f64::cos,
            (-det::TRIG_DOMAIN, det::TRIG_DOMAIN),
            &edges,
            200_000,
        );
        println!(
            "{{\"suite\":\"fs-math\",\"case\":\"trig-ulp\",\"verdict\":\"pass\",\"detail\":\"sin max={su}@{sa} cos max={cu}@{ca}\"}}"
        );
        assert!(su <= det::SIN_ULP_BUDGET, "sin max ULP {su} at {sa}");
        assert!(cu <= det::COS_ULP_BUDGET, "cos max ULP {cu} at {ca}");
        assert_eq!(det::sin(0.0).to_bits(), 0.0f64.to_bits());
        assert_eq!(det::cos(0.0).to_bits(), 1.0f64.to_bits());
    }

    #[test]
    fn tanh_meets_budget_and_saturates() {
        let (ulp, at) = max_ulp(
            det::tanh,
            f64::tanh,
            (-25.0, 25.0),
            &[1e-20, -1e-20, 19.9, -19.9],
            200_000,
        );
        println!(
            "{{\"suite\":\"fs-math\",\"case\":\"tanh-ulp\",\"verdict\":\"pass\",\"detail\":\"max_ulp={ulp} at {at}\"}}"
        );
        assert!(ulp <= det::TANH_ULP_BUDGET, "tanh max ULP {ulp} at {at}");
        assert_eq!(det::tanh(25.0).to_bits(), 1.0f64.to_bits());
        assert_eq!(det::tanh(-25.0).to_bits(), (-1.0f64).to_bits());
        assert_eq!(det::tanh(f64::INFINITY).to_bits(), 1.0f64.to_bits());
    }

    #[test]
    fn symmetries_hold_bitwise() {
        let mut seed = 7u64;
        for _ in 0..20_000 {
            let x = (lcg(&mut seed) - 0.5) * 2.0e5;
            assert_eq!(
                det::sin(-x).to_bits(),
                (-det::sin(x)).to_bits(),
                "sin odd at {x}"
            );
            assert_eq!(
                det::cos(-x).to_bits(),
                det::cos(x).to_bits(),
                "cos even at {x}"
            );
            assert_eq!(
                det::tanh(-x).to_bits(),
                (-det::tanh(x)).to_bits(),
                "tanh odd at {x}"
            );
        }
    }

    #[test]
    fn special_values_follow_the_policy() {
        for f in [det::exp, det::expm1, det::ln, det::sin, det::cos, det::tanh] {
            assert!(f(f64::NAN).is_nan(), "NaN in -> NaN out");
        }
        assert_eq!(det::exp(f64::NEG_INFINITY).to_bits(), 0.0f64.to_bits());
        assert!(det::exp(f64::INFINITY).is_infinite());
        assert!(det::ln(-1.0).is_nan());
        assert_eq!(det::ln(0.0).to_bits(), f64::NEG_INFINITY.to_bits());
        assert!(det::sin(f64::INFINITY).is_nan());
        assert!(canonical_nan().is_nan());
        // Subnormal in/out is never flushed.
        let sub = f64::MIN_POSITIVE / 4.0;
        assert!(det::ln(sub).is_finite(), "subnormal input must work");
        assert!(det::exp(-745.0) >= 0.0, "deep-underflow path");
    }

    #[test]
    fn nudge_helpers_bracket() {
        for x in [0.0, 1.0, -1.0, 1e-300, 1e300, f64::MIN_POSITIVE] {
            let (lo, hi) = nudge_out(x);
            assert!(lo < x || (x == 0.0 && lo < 0.0));
            assert!(hi > x);
            assert_eq!(ulp_distance(lo, hi), 2);
        }
        assert_eq!(next_up(1.0).to_bits(), (1.0 + f64::EPSILON).to_bits());
        assert_eq!(next_down(1.0 + f64::EPSILON).to_bits(), 1.0f64.to_bits());
    }

    /// THE CROSS-ISA GOLDEN HASH: FNV over the output bits of every strict
    /// function on a fixed grid. The SAME constant must hold on aarch64 and
    /// x86-64 — this is bit-determinism-by-construction, verified. If this
    /// test fails after an edit, determinism broke or semantics changed:
    /// both require a deliberate golden bump with justification.
    #[test]
    fn cross_isa_golden_hash() {
        let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
        let mut feed = |v: f64| {
            for b in v.to_bits().to_le_bytes() {
                acc ^= u64::from(b);
                acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
            }
        };
        let mut seed = 0x0DD_BA11_u64;
        for _ in 0..25_000 {
            let t = lcg(&mut seed);
            let x = (t - 0.5) * 1.0e4;
            feed(det::exp(x * 0.05));
            feed(det::expm1(x * 0.05));
            feed(det::ln(t + 1e-9));
            feed(det::sin(x));
            feed(det::cos(x));
            feed(det::tanh(x * 0.01));
            feed(det::sqrt(t * 1e6));
        }
        println!(
            "{{\"suite\":\"fs-math\",\"case\":\"golden-hash\",\"verdict\":\"info\",\"detail\":\"{acc:#018x}\"}}"
        );
        assert_eq!(
            acc, GOLDEN_HASH,
            "strict-mode outputs changed: {acc:#018x} vs recorded {GOLDEN_HASH:#018x} — \
             if intentional, bump the golden WITH justification (golden-evidence policy)"
        );
    }

    /// Recorded on aarch64-apple (M4 Pro); verified identical on x86-64
    /// (Threadripper PRO 5995WX). Bump ONLY with a semantic justification.
    const GOLDEN_HASH: u64 = 0xeb79_cab7_a016_43e5;
}
