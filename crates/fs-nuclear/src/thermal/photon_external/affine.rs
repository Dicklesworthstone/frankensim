//! Original affine source-power arithmetic, independent of transport and calorics.

use fs_math::dd::Dd;

// Caller admits finite nonnegative powers and left <= time <= right, left < right.
pub(super) fn interpolate(a: f64, b: f64, left: f64, right: f64, time: f64) -> Dd {
    if time == left || a == b {
        return Dd::from_f64(a);
    }
    if time == right {
        return Dd::from_f64(b);
    }
    let span = Dd::from_f64(right) - Dd::from_f64(left);
    // The nearer endpoint preserves weak positive power at pulse edges.
    if time - left <= right - time {
        Dd::from_f64(a)
            + (Dd::from_f64(b) - Dd::from_f64(a))
                * ((Dd::from_f64(time) - Dd::from_f64(left)) / span)
    } else {
        Dd::from_f64(b)
            + (Dd::from_f64(a) - Dd::from_f64(b))
                * ((Dd::from_f64(right) - Dd::from_f64(time)) / span)
    }
}

// Return integrated joules and observed final scalar rounding; no relative
// certificate is inferred from compensated arithmetic at subnormal boundaries.
pub(super) fn affine_energy(a: Dd, b: Dd, duration: f64) -> Option<(f64, f64)> {
    if !duration.is_finite()
        || duration <= 0.0
        || [a, b].iter().any(|p| {
            !p.hi.is_finite() || !p.lo.is_finite() || p.hi < 0.0 || (p.hi == 0.0 && p.lo < 0.0)
        })
    {
        return None;
    }
    let peak = a.hi.max(b.hi);
    if peak == 0.0 {
        return Some((0.0, 0.0));
    }
    let scale = Dd::from_f64(peak);
    let shape = (a / scale + b / scale) * Dd::from_f64(0.5);
    // Scale before summing endpoints; restore a physical energy before halving
    // subnormal powers. The fallback admits a finite mean when peak*dt overflows.
    let energy = if (peak * duration).is_finite() {
        (scale * Dd::from_f64(duration)) * shape
    } else {
        ((scale * Dd::from_f64(0.5)) * Dd::from_f64(duration)) * (shape * Dd::from_f64(2.0))
    };
    let value = energy.hi;
    if !value.is_finite() || value <= 0.0 || !energy.lo.is_finite() {
        return None;
    }
    let mut error = energy.lo.abs();
    if value < f64::MIN_POSITIVE {
        error = error.next_up();
    }
    Some((value, error))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 4e-15 * expected.abs().max(f64::MIN_POSITIVE),
            "{actual:?} != {expected:?}"
        );
    }

    #[test]
    fn g1_affine_power_matches_an_independent_polynomial_integral() {
        let p = |t| interpolate(2.0, 8.0, 0.0, 2.0, t);
        let (energy, error) = affine_energy(p(0.2), p(1.7), 1.5).unwrap();
        close(
            energy,
            2.0 * (1.7 - 0.2) + 1.5 * (1.7_f64.powi(2) - 0.2_f64.powi(2)),
        );
        assert!(error >= 0.0);
        let negative_slope = affine_energy(p(1.7), p(0.2), 1.5).unwrap().0;
        assert_eq!(energy, negative_slope);
    }

    #[test]
    fn g3_affine_power_subdivision_and_time_translation() {
        let whole = affine_energy(Dd::from_f64(2.0), Dd::from_f64(8.0), 2.0)
            .unwrap()
            .0;
        let mut split = Dd::ZERO;
        for (a, b) in [(0.0, 0.125), (0.125, 0.7), (0.7, 1.4), (1.4, 2.0)] {
            split = split
                + Dd::from_f64(
                    affine_energy(
                        interpolate(2.0, 8.0, 0.0, 2.0, a),
                        interpolate(2.0, 8.0, 0.0, 2.0, b),
                        b - a,
                    )
                    .unwrap()
                    .0,
                );
        }
        close(split.hi, whole);
        close(interpolate(2.0, 8.0, 128.0, 130.0, 128.5).hi, 3.5);
    }

    #[test]
    fn g1_affine_power_retains_subnormal_sources_and_finite_overflow_limits() {
        let tiny = f64::from_bits(1);
        close(
            affine_energy(Dd::from_f64(tiny), Dd::from_f64(tiny), 1e300)
                .unwrap()
                .0,
            tiny * 1e300,
        );
        close(
            affine_energy(Dd::from_f64(tiny), Dd::ZERO, 1e300)
                .unwrap()
                .0,
            (tiny * 1e300) * 0.5,
        );
        close(
            affine_energy(Dd::from_f64(f64::MAX), Dd::ZERO, 1.5)
                .unwrap()
                .0,
            f64::MAX * 0.75,
        );
        assert_eq!(affine_energy(Dd::ZERO, Dd::ZERO, 1.0), Some((0.0, 0.0)));
        assert!(affine_energy(Dd::from_f64(tiny), Dd::ZERO, tiny).is_none());
        assert!(affine_energy(Dd::from_f64(f64::MAX), Dd::from_f64(f64::MAX), 2.0).is_none());
    }

    #[test]
    fn g0_affine_power_preserves_pulse_edges_and_refuses_bad_states() {
        let time = 1.0_f64.next_down();
        let tail = interpolate(1e300, 0.0, 0.0, 1.0, time);
        close(tail.hi, (1.0 - time) * 1e300);
        assert!(tail.hi > 0.0);
        assert_eq!(interpolate(1e300, 1e-300, 0.0, 1.0, 1.0).hi, 1e-300);
        assert!(affine_energy(Dd::from_f64(f64::NAN), Dd::ZERO, 1.0).is_none());
        assert!(affine_energy(Dd::from_f64(-1.0), Dd::ZERO, 1.0).is_none());
        assert!(affine_energy(Dd::ONE, Dd::ONE, 0.0).is_none());
    }
}
