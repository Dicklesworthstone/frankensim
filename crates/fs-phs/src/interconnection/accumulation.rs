//! Two-word contraction arithmetic retained through both physical map factors.
//! No pHS/storage substitute: this is the scalar accumulator used by the native
//! interconnection owner. A refused local scratch sum is discarded by its caller.

#[derive(Clone, Copy, Default)]
pub(crate) struct Sum {
    high: f64,
    low: f64,
}

impl Sum {
    /// A complete checked contraction, retaining both words until its physical
    /// consumer has added the other sources and applied its integration weight.
    pub(crate) fn dot(left: &[f64], right: &[f64]) -> Option<Self> {
        if left.len() != right.len() {
            return None;
        }
        let mut sum = Self::default();
        for (&a, &b) in left.iter().zip(right) {
            sum.product(a, b)?;
        }
        Some(sum)
    }

    fn add(&mut self, value: f64) -> Option<()> {
        let next = self.high + value;
        self.low += if self.high.abs() >= value.abs() {
            (self.high - next) + value
        } else {
            (value - next) + self.high
        };
        self.high = next;
        (self.high.is_finite() && self.low.is_finite()).then_some(())
    }

    /// Retain both the represented product and its original FMA remainder.
    /// Nonzero products that vanish below f64's range refuse, rather than being
    /// silently treated as a missing physical coefficient.
    pub(crate) fn product(&mut self, a: f64, b: f64) -> Option<()> {
        let product = a * b;
        if !product.is_finite() || (product == 0.0 && a != 0.0 && b != 0.0) {
            return None;
        }
        self.add(product)?;
        self.add(a.mul_add(b, -product))
    }

    /// Peak-phasor mean work carries a factor of one half. Halve the larger
    /// scalar before multiplication: a tiny quadrature must not disappear
    /// before a large conjugate restores it, and a finite half-product must
    /// not be refused just because the unhalved product would overflow.
    pub(crate) fn half_product(&mut self, a: f64, b: f64) -> Option<()> {
        let (larger, smaller) = if a.abs() >= b.abs() { (a, b) } else { (b, a) };
        let half = 0.5 * larger;
        if half == 0.0 && larger != 0.0 && smaller != 0.0 { return None; }
        self.product(half, smaller)
    }

    /// Carry BOTH first-stage words through the next physical contraction.
    /// Rounding `source.value()` first could remove an entire coupled mode.
    pub(crate) fn add_scaled(&mut self, source: Self, coefficient: f64) -> Option<()> {
        self.product(source.high, coefficient)?;
        self.product(source.low, coefficient)
    }

    /// Carry both original contraction words through a half-scaled physical
    /// work term without first rounding the contraction or its product.
    pub(crate) fn add_half_scaled(&mut self, source: Self, coefficient: f64) -> Option<()> {
        self.half_product(source.high, coefficient)?;
        self.half_product(source.low, coefficient)
    }

    /// Round only at the final structure entry, whose opposite transpose is
    /// then installed from this same represented value.
    pub(crate) fn value(self) -> Option<f64> {
        let value = self.high + self.low;
        value.is_finite().then_some(value)
    }
}

#[cfg(test)]
mod tests {
    use super::Sum;

    #[test]
    fn unit_coupling_survives_large_first_stage_values_and_either_order() {
        let mut positive = Sum::default();
        positive.product(1.0, 1e16).unwrap();
        positive.product(1.0, 1.0).unwrap();
        let mut negative = Sum::default();
        negative.product(1.0, -1e16).unwrap();
        // Rounding the individual contractions gives cancelling 1e16 terms.
        assert_eq!(positive.value().unwrap() + negative.value().unwrap(), 0.0);
        for order in [[positive, negative], [negative, positive]] {
            let mut coupled = Sum::default();
            for value in order {
                coupled.add_scaled(value, 1.0).unwrap();
            }
            assert_eq!(coupled.value().unwrap(), 1.0);
        }
    }

    #[test]
    fn a_product_rounding_remainder_remains_after_large_cancellation() {
        let delta = 2.0_f64.powi(-27);
        let mut value = Sum::default();
        value.product(1.0 + delta, 1.0 - delta).unwrap();
        value.product(-1.0, 1.0).unwrap();
        assert_eq!(value.value().unwrap(), -2.0_f64.powi(-54));
    }

    #[test]
    fn scaled_intermediate_words_retain_weak_work_conjugate_coefficients() {
        let large = 2.0_f64.powi(53);
        let mut source = Sum::default();
        source.product(large, 1.0).unwrap();
        source.product(1.0, 1.0).unwrap();
        let mut result = Sum::default();
        result.add_scaled(source, 3.0).unwrap();
        result.product(-large, 3.0).unwrap();
        assert_eq!(result.value().unwrap(), 3.0);
        let mut zero = Sum::default();
        zero.add_scaled(source, 0.0).unwrap();
        assert_eq!(zero.value().unwrap(), 0.0);
    }

    #[test]
    fn peak_mean_products_preserve_the_finite_range_and_weak_quadratures() {
        for (a, b, expected) in [
            (f64::MAX, 2.0, f64::MAX),
            (2.0, f64::MAX, f64::MAX),
            (f64::from_bits(1), 2.0_f64.powi(1023), 2.0_f64.powi(-52)),
            (2.0_f64.powi(1023), f64::from_bits(1), 2.0_f64.powi(-52)),
        ] {
            let mut sum = Sum::default();
            sum.half_product(a, b).unwrap();
            assert_eq!(sum.value().unwrap(), expected);
        }
        assert!(Sum::default().half_product(f64::from_bits(1), 1.0).is_none());
        assert!(Sum::default().half_product(f64::INFINITY, 0.0).is_none());
    }

    #[test]
    fn peak_mean_work_keeps_the_unrounded_input_map_words() {
        let large = 2.0_f64.powi(54);
        let mut positive = Sum::default();
        positive.product(large, 1.0).unwrap();
        positive.product(1.0, 1.0).unwrap();
        let mut negative = Sum::default();
        negative.product(-large, 1.0).unwrap();
        let mut work = Sum::default();
        work.add_half_scaled(positive, 2.0).unwrap();
        work.add_half_scaled(negative, 2.0).unwrap();
        assert_eq!(work.value().unwrap(), 1.0);
    }

    #[test]
    fn nonfinite_overflow_and_nonzero_underflow_refuse() {
        for (a, b) in [
            (1e300, 1e300),
            (1e-300, 1e-300),
            (f64::INFINITY, 0.0),
            (f64::NAN, 1.0),
        ] {
            assert!(Sum::default().product(a, b).is_none());
        }
        let mut zero = Sum::default();
        zero.product(0.0, 1e300).unwrap();
        zero.product(1e-300, 0.0).unwrap();
        assert_eq!(zero.value().unwrap(), 0.0);
    }

    #[test]
    fn simultaneous_source_work_retains_small_drive_and_regeneration() {
        let large = 2.0_f64.powi(53);
        for map in [[large, 1.0, -large], [-large, 1.0, large], [1.0, large, -large]] {
            for sign in [-1.0, 1.0] {
                let force = Sum::dot(&map, &[sign; 3]).unwrap();
                assert_eq!(force.value().unwrap(), sign);
                let output = map.map(|entry| 0.625 * entry);
                let power = Sum::dot(&[sign; 3], &output).unwrap();
                let mut work = Sum::default();
                work.add_scaled(power, 0.25).unwrap();
                assert_eq!(work.value().unwrap(), sign * 0.15625);
            }
        }
    }

    #[test]
    fn additional_source_cancels_before_duration_rounding() {
        let large = 2.0_f64.powi(53);
        let mut power = Sum::dot(&[1.0, 1.0], &[large, 1.0]).unwrap();
        power.product(1.0, -large).unwrap();
        let mut work = Sum::default();
        work.add_scaled(power, 0.125).unwrap();
        assert_eq!(work.value().unwrap(), 0.125);
        assert!(Sum::dot(&[1.0], &[]).is_none());
        assert!(Sum::dot(&[0.0], &[f64::NAN]).is_none());
        assert!(Sum::dot(&[f64::MAX], &[2.0]).is_none());
    }
}
