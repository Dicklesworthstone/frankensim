//! Positive radiometry and complementary directional proposal weights.

use super::EnvironmentCameraError;
use fs_math::det;

pub(super) fn power_weight(sample: f64, other: f64) -> Result<f64, EnvironmentCameraError> {
    if !sample.is_finite() || sample <= 0.0 || !other.is_finite() || other < 0.0 {
        return Err(EnvironmentCameraError::NumericalRange);
    }
    if other == 0.0 {
        return Ok(1.0);
    }
    let result = if sample >= other {
        let ratio = other / sample;
        1.0 / (1.0 + ratio * ratio)
    } else {
        let ratio = sample / other;
        let square = ratio * ratio;
        square / (1.0 + square)
    };
    if result > 0.0 && result.is_finite() {
        Ok(result)
    } else {
        // A later bright source can rescue a tiny weight. Refuse a positive
        // unrepresentable partition instead of declaring exact extinction.
        Err(EnvironmentCameraError::NumericalRange)
    }
}

pub(super) fn positive_product(
    numerator: &[f64],
    denominator: f64,
) -> Result<f64, EnvironmentCameraError> {
    if !denominator.is_finite()
        || denominator <= 0.0
        || numerator
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err(EnvironmentCameraError::NumericalRange);
    }
    if numerator.contains(&0.0) {
        return Ok(0.0);
    }
    let value = numerator.iter().product::<f64>() / denominator;
    if value.is_finite() && value > 0.0 {
        return Ok(value);
    }
    let exponent =
        numerator.iter().map(|&value| det::ln(value)).sum::<f64>() - det::ln(denominator);
    let value = det::exp(exponent);
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(EnvironmentCameraError::NumericalRange)
    }
}
