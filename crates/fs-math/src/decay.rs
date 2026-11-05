//! Deterministic scalar exponential decay and relaxation.
//!
//! The exposure is the nonnegative integrated hazard `x = integral(rate dt)`;
//! for a constant half-life it is `ln(2) * elapsed / half_life`. These numerical
//! kernels do not prescribe a physical rate, units, or a material's equilibrium.
//! Physical owners admit those quantities before calling this module.
//!
//! The lost fraction is evaluated independently with `expm1`, so a tiny loss is
//! retained even when the surviving fraction rounds to one. Scaled survivors
//! avoid prematurely underflowing an exponential whose product is representable.
//! All operations use the strict elementary-function operation tree.

use crate::{canonical_nan, det};

/// Surviving fraction `exp(-exposure)`, including the infinite-exposure limit.
/// Negative or NaN exposure returns the canonical NaN.
#[must_use]
pub fn remaining_fraction(exposure: f64) -> f64 {
    if exposure >= 0.0 {
        det::exp(-exposure)
    } else {
        canonical_nan()
    }
}

/// Lost fraction `1 - exp(-exposure)`, without cancellation near zero.
/// Negative or NaN exposure returns the canonical NaN.
#[must_use]
pub fn decayed_fraction(exposure: f64) -> f64 {
    if exposure == 0.0 {
        0.0
    } else if exposure > 0.0 {
        -det::expm1(-exposure)
    } else {
        canonical_nan()
    }
}

/// Attenuate a finite signed quantity, preserving representable late tails.
///
/// A cached normal attenuation may be multiplied directly by the physical
/// quantity in a hot loop. When the attenuation would be subnormal, this helper
/// splits the exponent before scaling: `value * exp(-x/2) * exp(-x/2)`.
/// In particular, `1e300 * exp(-800)` is nonzero even though `exp(-800)` is zero.
/// Invalid exposure or a nonfinite value returns the canonical NaN.
#[must_use]
pub fn attenuate(value: f64, exposure: f64) -> f64 {
    if !value.is_finite() || !(exposure >= 0.0) {
        return canonical_nan();
    }
    if exposure == 0.0 || value == 0.0 {
        return value;
    }
    let remaining = remaining_fraction(exposure);
    if remaining >= f64::MIN_POSITIVE {
        value * remaining
    } else {
        let half = det::exp(-0.5 * exposure);
        (value * half) * half
    }
}

/// Relax a finite signed value toward a finite fixed equilibrium.
///
/// This is the exact scalar affine solution in real arithmetic. For short
/// intervals it accumulates the small change onto the initial value; for long
/// intervals it retains the decaying deviation from equilibrium. If forming
/// that deviation overflows, the equivalent convex combination is used.
/// An exactly zero exposure returns the initial value unchanged. Invalid input
/// returns the canonical NaN. No interval enclosure is claimed.
#[must_use]
pub fn relax_toward(initial: f64, equilibrium: f64, exposure: f64) -> f64 {
    if !initial.is_finite() || !equilibrium.is_finite() || !(exposure >= 0.0) {
        return canonical_nan();
    }
    if exposure == 0.0 || initial == equilibrium {
        return initial;
    }
    let difference = equilibrium - initial;
    if difference.is_finite() {
        if exposure <= core::f64::consts::LN_2 {
            initial + difference * decayed_fraction(exposure)
        } else {
            equilibrium + attenuate(-difference, exposure)
        }
    } else {
        attenuate(initial, exposure) + equilibrium * decayed_fraction(exposure)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64, relative: f64) {
        assert!(
            (actual - expected).abs() <= relative * expected.abs(),
            "actual={actual:e}, expected={expected:e}"
        );
    }

    #[test]
    fn g0_half_lives_and_tiny_losses() {
        for halves in [1, 2, 10, 100] {
            close(
                remaining_fraction(f64::from(halves) * core::f64::consts::LN_2),
                2.0_f64.powi(-halves),
                2e-14,
            );
        }
        assert_eq!(remaining_fraction(1e-100), 1.0);
        close(decayed_fraction(1e-100), 1e-100, 1e-15);
        assert_eq!(decayed_fraction(0.0).to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn g0_scaled_survivors_do_not_underflow_prematurely() {
        for exposure in [709.0, 745.0, 800.0, 1000.0, 1400.0] {
            // The oracle exponent remains normal and uses a separate libm path.
            let expected = (1e300_f64.ln() - exposure).exp();
            close(attenuate(1e300, exposure), expected, 2e-13);
            close(attenuate(-1e300, exposure), -expected, 2e-13);
        }
        assert_eq!(attenuate(1e300, f64::INFINITY), 0.0);
        assert_eq!(attenuate(1e300, 1e300), 0.0);
        // A rounded subnormal fraction can also erase a sub-unit donor.
        // 100-digit decimal evaluation rounds this product to the least
        // positive binary64 value, while 0.4 * exp(-744.1) rounds to zero.
        assert_eq!(attenuate(0.4, 744.1), f64::from_bits(1));
    }

    #[test]
    fn g0_small_relaxation_preserves_initial_state_and_resolved_change() {
        assert_eq!(relax_toward(300.0, 1e100, 0.0), 300.0);
        close(relax_toward(300.0, 1e100, 1e-100), 301.0, 1e-15);
        close(
            relax_toward(1e300, 0.0, 800.0),
            attenuate(1e300, 800.0),
            1e-15,
        );
        assert_eq!(relax_toward(300.0, 400.0, core::f64::consts::LN_2), 350.0);
        assert_eq!(relax_toward(300.0, 400.0, f64::INFINITY), 400.0);
        assert_eq!(relax_toward(f64::MAX, f64::MAX, 1.0), f64::MAX);
        assert_eq!(
            relax_toward(f64::MAX, -f64::MAX, core::f64::consts::LN_2),
            0.0
        );
    }

    #[test]
    fn g3_affine_semigroup_and_domains() {
        for (start, equilibrium, first, second) in [
            (300.0, 400.0, 0.01, 0.2),
            (-7.0, 2.0, 2.0, 3.0),
            (1e300, 0.0, 400.0, 400.0),
        ] {
            let split = relax_toward(relax_toward(start, equilibrium, first), equilibrium, second);
            close(
                split,
                relax_toward(start, equilibrium, first + second),
                2e-14,
            );
        }
        for exposure in [-1.0, f64::NEG_INFINITY, f64::NAN] {
            assert!(remaining_fraction(exposure).is_nan());
            assert!(decayed_fraction(exposure).is_nan());
            assert!(attenuate(1.0, exposure).is_nan());
            assert!(relax_toward(1.0, 2.0, exposure).is_nan());
        }
        assert!(attenuate(f64::INFINITY, 1.0).is_nan());
        assert!(relax_toward(1.0, f64::NAN, 0.0).is_nan());
    }
}
