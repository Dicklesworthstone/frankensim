//! Full-cycle norm of a static field plus an elliptical harmonic perturbation.
//!
//! With x=Re(z), y=Im(z), v=(cos(t),sin(t)), the squared norm is
//! c + v^T A v + 2 g^T v, where c=b.b, A=[x y]^T[x y], g=[x y]^T b.
//! In A's principal axes let its largest eigenvalue be a, its gap be d>=0,
//! and the rotated g be h. For ANY delta>0, completing the square gives
//!
//! max ||b+Re(z exp(-i t))||^2 <= c+a+delta+h1^2/delta+h2^2/(d+delta).
//!
//! Minimizing this dual bound is the two-dimensional trust-region secular
//! equation hypot(h1/delta,h2/(d+delta))=1. A bounded search need not reach
//! its root to retain an upper bound. The singular/hard case has an analytic
//! limit. No phase sampling or triangle-bound admission is used.

use fs_math::{c64::C64, dd::Dd, det};

fn dot<const N: usize>(a: &[f64; N], b: &[f64; N]) -> Dd {
    a.iter().zip(b).fold(Dd::ZERO, |sum, (&a, &b)| {
        sum + Dd::from_f64(a) * Dd::from_f64(b)
    })
}

/// Conservative full-cycle bound, in the input vector's units. Invalid inputs
/// or an unrepresentable upper bound return infinity, never NaN or a silent
/// underestimate. Intended for the original three-/six-component cell fields.
///
/// Power-of-two scaling avoids overflow and normalization roundoff; retained
/// Gram/projection products protect cancellation. A binary64 roundoff allowance
/// covers the small eigensystem and dual evaluation, followed by outward rounding.
/// At most 64 scalar bisections run without allocation; the caller checkpoints
/// cancellation at its cell boundary.
pub(super) fn biased_cycle_norm<const N: usize>(mean: [f64; N], amplitude: [C64; N]) -> f64 {
    let mut maximum = 0.0_f64;
    for (&b, z) in mean.iter().zip(&amplitude) {
        if !b.is_finite() || !z.re.is_finite() || !z.im.is_finite() {
            return f64::INFINITY;
        }
        maximum = maximum.max(b.abs()).max(z.re.abs()).max(z.im.abs());
    }
    if maximum == 0.0 {
        return 0.0;
    }
    // Use the leading power of two, or MIN_POSITIVE for subnormal inputs.
    // Normalized components are below two; discarded extreme-ratio underflows
    // are far below the final roundoff allowance relative to the largest input.
    let scale = f64::from_bits((maximum.to_bits() & 0x7ff0_0000_0000_0000).max(1_u64 << 52));
    let b = mean.map(|value| value / scale);
    let x = amplitude.map(|value| value.re / scale);
    let y = amplitude.map(|value| value.im / scale);
    let c = dot(&b, &b).to_f64();
    let aa = dot(&x, &x);
    let ab = dot(&x, &y);
    let bb = dot(&y, &y);
    let gx = dot(&b, &x);
    let gy = dot(&b, &y);
    let difference = (aa - bb).to_f64();
    let twice_cross = 2.0 * ab.to_f64();
    let gap = det::hypot(difference, twice_cross);
    let largest = ((aa + bb + Dd::from_f64(gap)) * Dd::from_f64(0.5)).to_f64();

    // Choose the noncancelling largest-eigenvector formula. Its perpendicular
    // is the other phase axis, even for a rank-one or nearly circular ellipse.
    let (qx, qy) = if gap == 0.0 {
        (1.0, 0.0)
    } else {
        let (u, v) = if difference >= 0.0 {
            (gap + difference, twice_cross)
        } else {
            (twice_cross, gap - difference)
        };
        let length = det::hypot(u, v);
        (u / length, v / length)
    };
    let h1 = (gx * Dd::from_f64(qx) + gy * Dd::from_f64(qy))
        .to_f64()
        .abs();
    let h2 = (gy * Dd::from_f64(qx) - gx * Dd::from_f64(qy))
        .to_f64()
        .abs();
    let magnitude = det::hypot(h1, h2);
    let correction = if gap == 0.0 {
        2.0 * magnitude
    } else if h2 == 0.0 {
        2.0 * h1
    } else if h1 == 0.0 && h2 <= gap {
        // Hard case: the minor-axis optimum leaves room for a major-axis
        // component. The limiting dual delta=0 is finite and attainable.
        h2 * (h2 / gap)
    } else {
        let (mut low, mut high) = (0.0, magnitude);
        for _ in 0..64 {
            let middle = 0.5 * low + 0.5 * high;
            if middle == low || middle == high {
                break;
            }
            if det::hypot(h1 / middle, h2 / (gap + middle)) > 1.0 {
                low = middle;
            } else {
                high = middle;
            }
        }
        // Multiply after division: h_i^2 can underflow even when the quotient
        // is representable. Every term is nonnegative; no dual cancellation.
        high + h1 * (h1 / high) + h2 * (h2 / (gap + high))
    };
    let squared = c + largest + correction;
    // Perturbing the rounded 2x2 eigensystem changes a unit-circle quadratic
    // by at most its coefficient perturbations. Guard those O(eps) terms,
    // not a relative secular-root error (which is ill-conditioned in hard cases).
    let allowance = 64.0
        * (N as f64 + 1.0)
        * f64::EPSILON
        * (c + aa.to_f64() + bb.to_f64() + 2.0 * magnitude + correction);
    let bound = scale * det::sqrt((squared + allowance).next_up());
    if bound.is_finite() {
        bound.next_up()
    } else {
        f64::INFINITY
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const ZERO: C64 = C64::new(0.0, 0.0);

    fn tight(actual: f64, exact: f64) {
        assert!(actual >= exact, "bound {actual:e} below {exact:e}");
        assert!(
            actual <= exact * (1.0 + 5e-13),
            "loose bound {actual:e} for {exact:e}"
        );
    }

    #[test]
    fn analytic_constant_line_circle_and_centered_ellipse() {
        assert_eq!(biased_cycle_norm([0.0; 3], [ZERO; 3]), 0.0);
        tight(biased_cycle_norm([3.0, 4.0, 0.0], [ZERO; 3]), 5.0);
        tight(biased_cycle_norm([2.0], [C64::new(3.0, 4.0)]), 7.0);
        tight(
            biased_cycle_norm([0.0; 3], [C64::new(3.0, 0.0), C64::new(0.0, 4.0), ZERO]),
            4.0,
        );
        tight(
            biased_cycle_norm(
                [0.0, 0.0, 4.0],
                [C64::new(3.0, 0.0), C64::new(0.0, 3.0), ZERO],
            ),
            5.0,
        );
        tight(
            biased_cycle_norm([3.0, 4.0], [C64::new(2.0, 0.0), C64::new(0.0, 2.0)]),
            7.0,
        );
    }

    #[test]
    fn offset_ellipses_admit_cycles_the_triangle_bound_rejects() {
        let ellipse = [C64::new(3.0, 0.0), C64::new(0.0, 1.0)];
        // 9 cos^2(t)+(1+sin(t))^2 peaks at sin(t)=1/8, not an axis.
        let hard = biased_cycle_norm([0.0, 1.0], ellipse);
        tight(hard, det::sqrt(81.0 / 8.0));
        assert!(hard < 3.2); // Triangle bound would be 1+3=4.
        tight(biased_cycle_norm([0.0, 8.0], ellipse), 9.0);
        // A=diag(9,1), lambda=10 and v=(3/5,4/5) satisfy the global
        // positive-definite dual certificate. The maximizing field is (2,8).
        let regular = biased_cycle_norm([0.2, 7.2], ellipse);
        tight(regular, det::sqrt(68.0));
        assert!(regular < 8.25);
    }

    #[test]
    fn phase_rotation_and_nearly_cancelling_projections_preserve_the_bound() {
        let mean = [2.0, -2.0, 1e-12, 0.0, 0.0, 0.0];
        let amplitude = [
            C64::new(3.0, 1.0),
            C64::new(3.0, 1.0),
            C64::new(1e-12, -2e-12),
            ZERO,
            ZERO,
            ZERO,
        ];
        let rotated = amplitude.map(|z| z * C64::new(0.6, 0.8));
        let original = biased_cycle_norm(mean, amplitude);
        let shifted = biased_cycle_norm(mean, rotated);
        // The leading mean is orthogonal to the line: sqrt(8+20).
        tight(original, det::sqrt(28.0));
        assert!((original - shifted).abs() <= 5e-13 * original);
        for i in 0..128 {
            let phase = core::f64::consts::TAU * i as f64 / 128.0;
            let (sine, cosine) = (det::sin(phase), det::cos(phase));
            let norm = mean.iter().zip(amplitude).fold(0.0, |length, (&b, z)| {
                det::hypot(length, b + z.re * cosine + z.im * sine)
            });
            assert!(norm <= original);
        }
    }

    #[test]
    fn scaled_extremes_and_invalid_values_do_not_underestimate() {
        for scale in [1e-308, 1.0, 1e308] {
            let bound = biased_cycle_norm([0.3 * scale], [C64::new(0.4 * scale, 0.0)]);
            let exact = 0.3 * scale + 0.4 * scale;
            assert!(bound.is_finite() && bound >= exact);
            assert!((bound / scale - 0.7).abs() < 1e-12);
        }
        let subnormal = f64::from_bits(1);
        assert!(biased_cycle_norm([subnormal], [ZERO]) >= subnormal);
        assert_eq!(biased_cycle_norm([f64::NAN], [ZERO]), f64::INFINITY);
        assert_eq!(
            biased_cycle_norm([0.0], [C64::new(0.0, f64::INFINITY)]),
            f64::INFINITY
        );
        assert_eq!(
            biased_cycle_norm([f64::MAX], [C64::new(f64::MAX, 0.0)]),
            f64::INFINITY
        );
    }
}
