//! Inclusive zero-power coordinates for the original physical scan study.
use super::{ScanSensitivityError, finite, invalid, positive_coordinate};

/// How selected authored watts map to dimensionless optimizer coordinates.
/// This does not change the physical source law, its derivative or off entries
/// that the caller leaves unselected. All selected powers use the same map.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScanPowerCoordinates {
    /// P = reference_w * exp(x). Existing default; never reaches zero.
    PositiveLog,
    /// P = reference_w * x, x >= 0. The closed half-line includes exact off
    /// states and has dP/dx = reference_w even there. ConstrainedScanStudy owns
    /// the domain inequalities. Negative trials are not cooling or clipped heat.
    NonnegativeLinear,
}

impl ScanPowerCoordinates {
    pub(super) fn initial(self, power: f64, reference: f64) -> Result<f64, ScanSensitivityError> {
        if !power.is_finite() || power < 0.0 || !reference.is_finite() || reference <= 0.0 {
            return Err(invalid("power coordinates require nonnegative finite watts and a positive finite watt scale"));
        }
        match self {
            Self::PositiveLog => {
                if power == 0.0 { return Err(invalid("log-power controls require positive initial watts")); }
                finite(fs_math::det::ln(power) - fs_math::det::ln(reference))
            }
            Self::NonnegativeLinear => {
                let x = finite(power / reference)?;
                if power > 0.0 && x == 0.0 {
                    return Err(invalid("nonnegative power coordinate underflowed; choose a representable watt scale"));
                }
                Ok(if x == 0.0 { 0.0 } else { x })
            }
        }
    }

    pub(super) fn decode(self, x: f64, reference: f64) -> Result<f64, ScanSensitivityError> {
        finite(x)?;
        if !reference.is_finite() || reference <= 0.0 {
            return Err(invalid("power coordinates require a positive finite watt scale"));
        }
        match self {
            Self::PositiveLog => positive_coordinate(x, reference),
            Self::NonnegativeLinear => {
                if x < 0.0 { return Err(invalid("negative heater power is outside the selected design domain")); }
                let power = finite(reference * x)?;
                if x > 0.0 && power == 0.0 {
                    return Err(invalid("nonzero heater power underflowed; choose a representable watt scale"));
                }
                Ok(if power == 0.0 { 0.0 } else { power })
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g1_linear_power_map_includes_off_without_losing_the_activation_derivative() {
        let map = ScanPowerCoordinates::NonnegativeLinear;
        for scale in [1e-200_f64, 0.25, 4.0, 1e200] {
            assert_eq!(map.initial(0.0, scale).unwrap().to_bits(), 0.0_f64.to_bits());
            assert_eq!(map.decode(-0.0, scale).unwrap().to_bits(), 0.0_f64.to_bits());
            let epsilon = 1e-6;
            let derivative = (map.decode(epsilon, scale).unwrap() - map.decode(0.0, scale).unwrap()) / epsilon;
            assert!((derivative / scale - 1.0).abs() < 1e-14);
            for x in [0.25, 1.0, 4.0] {
                let power = map.decode(x, scale).unwrap();
                assert!((map.initial(power, scale).unwrap() / x - 1.0).abs() < 1e-14);
            }
        }
        assert!(ScanPowerCoordinates::PositiveLog.initial(0.0, 1.0).is_err());
        assert!(map.decode(-f64::MIN_POSITIVE, 1.0).is_err());
    }

    #[test]
    fn g0_nonfinite_or_unrepresentable_power_is_never_clipped_to_an_off_state() {
        let map = ScanPowerCoordinates::NonnegativeLinear;
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1.0] {
            assert!(map.decode(bad, 1.0).is_err());
            assert!(map.initial(bad, 1.0).is_err());
        }
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(map.decode(0.0, bad).is_err());
            assert!(map.initial(0.0, bad).is_err());
        }
        assert!(map.decode(1e-200, 1e-200).is_err());
        assert!(map.initial(1e-200, 1e200).is_err());
        assert!(map.decode(1e200, 1e200).is_err());
        assert!(map.initial(1e200, 1e-200).is_err());
    }
}
