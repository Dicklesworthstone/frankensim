//! Datum-free embedded energy errors and compensated physical heat accounts.

/// Difference between a native endpoint increment and the embedded first-order
/// increment first/stage_weight, divided by their change-based allowance.
/// The embedded increment is never made into an absolute material state.
pub(in crate::conjugate) fn embedded_ratio(
    first: f64,
    endpoint: f64,
    stage_weight: f64,
    absolute: f64,
    relative: f64,
) -> Option<f64> {
    if !first.is_finite()
        || !endpoint.is_finite()
        || !stage_weight.is_finite()
        || !(0.0..1.0).contains(&stage_weight)
        || stage_weight == 0.0
        || !absolute.is_finite()
        || absolute <= 0.0
        || !relative.is_finite()
        || !(0.0..1.0).contains(&relative)
    {
        return None;
    }
    // Multiply both numerator and denominator by the stage weight. This
    // avoids constructing an embedded energy outside original source support
    // and avoids division overflow in an otherwise admissible error ratio.
    let absolute = stage_weight * absolute;
    if absolute <= 0.0 {
        return None;
    }
    let numerator = stage_weight.mul_add(endpoint, -first).abs();
    let movement = (stage_weight * endpoint).abs().max(first.abs());
    let allowance = relative.mul_add(movement, absolute);
    if numerator.is_finite() && allowance.is_finite() {
        return Some(numerator / allowance);
    }
    // Finite physical increments can still overflow their comparison. Scale
    // this exceptional calculation before adding; ordinary small discrepancies
    // above retain the more accurate single fused subtraction.
    let scale = first.abs().max(endpoint.abs()).max(absolute);
    let first = first / scale;
    let endpoint = endpoint / scale;
    let numerator = stage_weight.mul_add(endpoint, -first).abs();
    let allowance = relative.mul_add(
        (stage_weight * endpoint).abs().max(first.abs()),
        absolute / scale,
    );
    Some(if allowance > 0.0 {
        numerator / allowance
    } else if numerator == 0.0 {
        0.0
    } else {
        f64::INFINITY
    })
}

/// Compensate signed physical heat changes, never absolute formation energies.
/// There is no source, integration or acceptance authority in this account.
#[derive(Clone, Copy, Default)]
pub(super) struct Sum {
    value: f64,
    correction: f64,
}

impl Sum {
    pub fn add(&mut self, value: f64) -> Option<()> {
        if !value.is_finite() {
            return None;
        }
        let next = self.value + value;
        if !next.is_finite() {
            return None;
        }
        let lost = if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else {
            (value - next) + self.value
        };
        let correction = self.correction + lost;
        if !correction.is_finite() {
            return None;
        }
        self.value = next;
        self.correction = correction;
        Some(())
    }

    /// Subtract accounts before rounding away either one's compensation.
    pub fn difference(mut self, other: Self) -> Option<f64> {
        self.add(-other.value)?;
        self.add(-other.correction)?;
        self.total()
    }

    pub fn total(self) -> Option<f64> {
        let result = self.value + self.correction;
        result.is_finite().then_some(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const GAMMA: f64 = 1.0 - std::f64::consts::FRAC_1_SQRT_2;

    #[test]
    fn embedded_discrepancy_is_second_order_for_decay() {
        // Original semidiscrete heat mode y'=-y, y0=1.
        let estimate = |dt: f64| {
            let a = GAMMA * dt;
            let y1 = 1.0 / (1.0 + a);
            let y2 = (1.0 - (1.0 - GAMMA) * dt * y1) / (1.0 + a);
            let ratio = embedded_ratio(y1 - 1.0, y2 - 1.0, a / dt, 1.0, 0.0).unwrap();
            let independently = (y2 - (1.0 - dt * y1)).abs();
            assert!((ratio - independently).abs() < 8.0 * f64::EPSILON);
            ratio
        };
        for dt in [0.04, 0.02, 0.01] {
            let quotient = estimate(dt) / estimate(dt * 0.5);
            assert!((3.9..4.1).contains(&quotient), "{quotient}");
        }
    }

    #[test]
    fn original_energy_datums_do_not_relax_the_allowance() {
        // Dyadic states keep this comparison exact even under large datums.
        let reference = embedded_ratio(0.25, 1.0, GAMMA, 0.01, 0.02).unwrap();
        for datum in [0.0, 1e9, -1e12] {
            let first = (datum + 0.25) - datum;
            let endpoint = (datum + 1.0) - datum;
            assert_eq!(
                embedded_ratio(first, endpoint, GAMMA, 0.01, 0.02),
                Some(reference),
            );
        }
        // Temperature can be exactly constant while this original energy
        // estimate still detects a nonzero temporal discrepancy.
        assert!(reference > 1.0);
    }

    #[test]
    fn scale_safe_comparison_does_not_construct_a_large_embedded_state() {
        let actual = embedded_ratio(-f64::MAX, 0.8 * f64::MAX, GAMMA, f64::MAX, 0.9).unwrap();
        let expected = (GAMMA * 0.8 + 1.0) / (GAMMA + 0.9);
        assert!((actual - expected).abs() < 8.0 * f64::EPSILON);
        assert!(actual > 1.0);
        assert_eq!(embedded_ratio(0.0, 0.0, GAMMA, 1e-300, 0.0), Some(0.0),);
        for (a, r) in [(0.0, 0.0), (1.0, -1.0), (1.0, 1.0), (f64::NAN, 0.0)] {
            assert!(embedded_ratio(0.0, 0.0, GAMMA, a, r).is_none());
        }
        assert!(embedded_ratio(0.0, 0.0, GAMMA, f64::from_bits(1), 0.0).is_none());
        assert!(embedded_ratio(f64::INFINITY, 1.0, GAMMA, 1.0, 0.0).is_none());
    }

    #[test]
    fn compensated_heat_retains_small_transfers_between_opposite_channels() {
        let mut heat = Sum::default();
        for value in [1e16, 1.0, -1e16] {
            heat.add(value).unwrap();
        }
        assert_eq!(heat.total(), Some(1.0));
        let saved = heat.total();
        assert!(heat.add(f64::INFINITY).is_none());
        assert_eq!(heat.total(), saved);
        let mut overflow = Sum::default();
        overflow.add(f64::MAX).unwrap();
        assert!(overflow.add(f64::MAX).is_none());
        assert_eq!(overflow.total(), Some(f64::MAX));
    }

    #[test]
    fn signed_net_heat_preserves_closure_between_large_separate_channels() {
        let mut source = Sum::default();
        let mut boundary = Sum::default();
        let mut net = Sum::default();
        for (input, output) in [(1e16, 1e16), (1.0, 0.0)] {
            source.add(input).unwrap();
            boundary.add(output).unwrap();
            net.add(input).unwrap();
            net.add(-output).unwrap();
        }
        // Separately rounded report channels cannot encode this 1 J balance.
        assert_eq!(source.total().unwrap() - boundary.total().unwrap(), 0.0);
        assert_eq!(net.total(), Some(1.0));
        let mut stored = Sum::default();
        stored.add(1.0).unwrap();
        assert_eq!(stored.difference(net), Some(0.0));
        stored.add(0.25).unwrap();
        assert_eq!(stored.difference(net), Some(0.25));

        // Closure also preserves compensation within the two signed accounts.
        let mut large = Sum::default();
        large.add(1e16).unwrap();
        let mut with_small = large;
        with_small.add(1.0).unwrap();
        assert_eq!(with_small.difference(large), Some(1.0));
        assert_eq!(large.difference(with_small), Some(-1.0));
    }
}
