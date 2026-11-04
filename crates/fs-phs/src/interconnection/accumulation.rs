//! Two-word contraction arithmetic retained through both physical map factors.
//! No pHS/storage substitute: this is the scalar accumulator used by the native
//! interconnection owner. A refused local scratch sum is discarded by its caller.

#[derive(Clone, Copy, Default)]
pub(crate) struct Sum {
    high: f64,
    low: f64,
}

impl Sum {
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

    /// Carry BOTH first-stage words through the next physical contraction.
    /// Rounding `source.value()` first could remove an entire coupled mode.
    pub(crate) fn add_scaled(&mut self, source: Self, coefficient: f64) -> Option<()> {
        self.product(source.high, coefficient)?;
        self.product(source.low, coefficient)
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
}
