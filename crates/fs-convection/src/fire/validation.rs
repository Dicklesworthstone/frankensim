//! Descriptive paired fire-prediction comparisons with explicit uncertainty gaps.
//!
//! Implements the unweighted log-ratio equations in the pinned FDS Validation
//! Guide's "Calculating Model Uncertainty" section. The separated relative
//! standard deviation and bias use that source's delta-method approximation;
//! they are not exact lognormal parameters. No regressogram weighting,
//! confidence interval, normality test, calibration or validation promotion is
//! inferred. Supplying unknown experimental uncertainty as zero is not valid.

use fs_math::det;

/// Exact primary source of the comparison and separation equations.
pub const FIRE_COMPARISON_SOURCE_V1: &str = "https://github.com/firemodels/fds/blob/f3523f3e09b565a5a16fdf36fc40c6ebaeec9789/Manuals/FDS_Validation_Guide/Error_Chapter.tex";

/// One comparable, positive modeled/observed quantity in identical units.
///
/// Match the quantity, location, averaging and event/time definition. For a
/// temperature comparison in the cited method, use a positive rise above the
/// respective ambient, not absolute kelvin or an arbitrary Celsius offset.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FirePredictionPairV1 {
    /// Model prediction in the same units and with the same definition as E.
    pub model: f64,
    /// Experimental observation E, strictly positive and finite.
    pub observed: f64,
}

/// Why model-only scatter and the FDS bias are or are not available.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FireSeparationStatusV1 {
    /// No experimental relative standard uncertainty was supplied.
    UnknownExperimentalUncertainty,
    /// Sample log variance is below twice the supplied experimental variance.
    /// The model scatter cannot be separated under the source's constraint.
    InsufficientScatter,
    /// Sample log variance is at least twice the experimental variance.
    Identified,
}

/// Descriptive statistics; no acceptance band or statistical confidence claim.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FireComparisonV1 {
    /// Number of independent comparisons supplied by the caller.
    pub pair_count: usize,
    /// Arithmetic mean of r=ln(M/E).
    pub mean_log_ratio: f64,
    /// Sample variance sum((r-mean)^2)/(n-1).
    pub sample_log_variance: f64,
    /// exp(mean(r)); this is not the experimental-uncertainty-adjusted FDS bias.
    pub geometric_ratio: f64,
    /// sqrt(sample_log_variance), combining model and experimental scatter.
    pub combined_log_scatter: f64,
    /// Caller-declared combined experimental relative standard uncertainty.
    /// It must include applicable propagated input uncertainty, not just a
    /// sensor's precision, and must be a standard rather than expanded value.
    pub experimental_relative_sigma: Option<f64>,
    /// sqrt(s^2-sigma_E^2), only when separation is identified.
    pub model_relative_sigma: Option<f64>,
    /// exp(mean(r)+(sigma_M^2-sigma_E^2)/2), only when identified.
    pub fds_bias_factor: Option<f64>,
    /// Explicit reason optional separated values are absent or present.
    pub separation: FireSeparationStatusV1,
}

/// An input or finite representation required by the statistics was refused.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FireComparisonErrorV1 {
    /// A sample variance requires at least two pairs.
    TooFewPairs,
    /// Exact rejected value and its zero-based pair index.
    InvalidPair {
        index: usize,
        field: &'static str,
        bits: u64,
    },
    /// Experimental relative standard uncertainty must be finite and >=0.
    InvalidExperimentalUncertainty { bits: u64 },
    /// A requested statistic cannot be represented as a finite admissible f64.
    Arithmetic,
}
impl core::fmt::Display for FireComparisonErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "fire prediction comparison: {self:?}")
    }
}
impl std::error::Error for FireComparisonErrorV1 {}

/// Compute unweighted log-ratio statistics in one pass with no allocations.
///
/// The caller selects comparable independent pairs and a defensible combined
/// experimental standard uncertainty. None leaves separation unavailable.
/// An insufficient sample variance never causes sigma_E to be reduced, sigma_M
/// to be floored, or an adjusted bias to be published. Equality at the declared
/// s^2 >= 2 sigma_E^2 identifiability boundary is admitted explicitly.
pub fn compare_fire_predictions_v1(
    pairs: &[FirePredictionPairV1],
    experimental_relative_sigma: Option<f64>,
) -> Result<FireComparisonV1, FireComparisonErrorV1> {
    if pairs.len() < 2 {
        return Err(FireComparisonErrorV1::TooFewPairs);
    }
    if let Some(sigma) = experimental_relative_sigma {
        if !sigma.is_finite() || sigma < 0.0 {
            return Err(FireComparisonErrorV1::InvalidExperimentalUncertainty {
                bits: sigma.to_bits(),
            });
        }
    }
    let mut mean = 0.0;
    let mut sum_squared_deviation = 0.0;
    for (index, pair) in pairs.iter().enumerate() {
        for (field, value) in [("model", pair.model), ("observed", pair.observed)] {
            if !value.is_finite() || value <= 0.0 {
                return Err(FireComparisonErrorV1::InvalidPair {
                    index,
                    field,
                    bits: value.to_bits(),
                });
            }
        }
        // Separate logarithms avoid overflow/underflow of M/E. Welford's
        // update avoids subtracting nearly equal squared sums.
        let log_ratio = det::ln(pair.model) - det::ln(pair.observed);
        let difference = log_ratio - mean;
        mean += difference / ((index + 1) as f64);
        sum_squared_deviation += difference * (log_ratio - mean);
    }
    let variance = sum_squared_deviation / ((pairs.len() - 1) as f64);
    let geometric_ratio = det::exp(mean);
    if !mean.is_finite()
        || !variance.is_finite()
        || variance < 0.0
        || !geometric_ratio.is_finite()
        || geometric_ratio <= 0.0
    {
        return Err(FireComparisonErrorV1::Arithmetic);
    }
    let mut result = FireComparisonV1 {
        pair_count: pairs.len(),
        mean_log_ratio: mean,
        sample_log_variance: variance,
        geometric_ratio,
        combined_log_scatter: det::sqrt(variance),
        experimental_relative_sigma,
        model_relative_sigma: None,
        fds_bias_factor: None,
        separation: FireSeparationStatusV1::UnknownExperimentalUncertainty,
    };
    if let Some(sigma_e) = experimental_relative_sigma {
        let variance_e = sigma_e * sigma_e;
        if !(2.0 * variance_e).is_finite() || (sigma_e > 0.0 && variance_e == 0.0) {
            return Err(FireComparisonErrorV1::Arithmetic);
        }
        result.separation = FireSeparationStatusV1::InsufficientScatter;
        if variance >= 2.0 * variance_e {
            let variance_m = variance - variance_e;
            let bias = det::exp(mean + 0.5 * (variance_m - variance_e));
            if !bias.is_finite() || bias <= 0.0 {
                return Err(FireComparisonErrorV1::Arithmetic);
            }
            result.model_relative_sigma = Some(det::sqrt(variance_m));
            result.fds_bias_factor = Some(bias);
            result.separation = FireSeparationStatusV1::Identified;
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(value: f64, expected: f64) {
        assert!(
            (value - expected).abs() < 2e-13 * expected.abs().max(1.0),
            "{value} != {expected}"
        );
    }

    fn logarithmic_example() -> [FirePredictionPairV1; 3] {
        let e = core::f64::consts::E;
        [1.0, e, e * e].map(|model| FirePredictionPairV1 {
            model,
            observed: 1.0,
        })
    }

    #[test]
    fn closed_form_log_statistics_and_identified_separation() {
        // Log ratios 0,1,2 give mean1 and unbiased sample variance1 exactly.
        let pairs = logarithmic_example();
        let unknown = compare_fire_predictions_v1(&pairs, None).unwrap();
        assert_eq!(unknown.pair_count, 3);
        close(unknown.mean_log_ratio, 1.0);
        close(unknown.sample_log_variance, 1.0);
        close(unknown.geometric_ratio, core::f64::consts::E);
        close(unknown.combined_log_scatter, 1.0);
        assert_eq!(
            unknown.separation,
            FireSeparationStatusV1::UnknownExperimentalUncertainty
        );
        assert_eq!(unknown.model_relative_sigma, None);
        assert_eq!(unknown.fds_bias_factor, None);
        let known = compare_fire_predictions_v1(&pairs, Some(0.5)).unwrap();
        assert_eq!(known.separation, FireSeparationStatusV1::Identified);
        close(known.model_relative_sigma.unwrap(), 0.75_f64.sqrt());
        close(known.fds_bias_factor.unwrap(), 1.25_f64.exp());
        let rescaled = pairs.map(|p| FirePredictionPairV1 {
            model: p.model * 1000.0,
            observed: p.observed * 1000.0,
        });
        let same = compare_fire_predictions_v1(&rescaled, None).unwrap();
        close(same.mean_log_ratio, unknown.mean_log_ratio);
        close(same.sample_log_variance, unknown.sample_log_variance);
    }

    #[test]
    fn insufficient_scatter_stays_unseparated_and_zero_boundary_is_explicit() {
        let result = compare_fire_predictions_v1(&logarithmic_example(), Some(0.8)).unwrap();
        assert_eq!(
            result.separation,
            FireSeparationStatusV1::InsufficientScatter
        );
        assert_eq!(result.experimental_relative_sigma, Some(0.8));
        assert_eq!(result.model_relative_sigma, None);
        assert_eq!(result.fds_bias_factor, None);
        let identical = [FirePredictionPairV1 {
            model: 2.0,
            observed: 2.0,
        }; 2];
        let zero = compare_fire_predictions_v1(&identical, Some(0.0)).unwrap();
        assert_eq!(zero.separation, FireSeparationStatusV1::Identified);
        assert_eq!(zero.model_relative_sigma, Some(0.0));
        assert_eq!(zero.fds_bias_factor, Some(1.0));
        assert_eq!(
            compare_fire_predictions_v1(&identical, None)
                .unwrap()
                .model_relative_sigma,
            None
        );
    }

    #[test]
    fn invalid_values_refuse_and_extreme_ratios_do_not_overflow() {
        assert_eq!(
            compare_fire_predictions_v1(&[], None),
            Err(FireComparisonErrorV1::TooFewPairs)
        );
        for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let pairs = [FirePredictionPairV1 {
                model: value,
                observed: 1.0,
            }; 2];
            assert!(matches!(
                compare_fire_predictions_v1(&pairs, None),
                Err(FireComparisonErrorV1::InvalidPair { index: 0, .. })
            ));
        }
        for sigma in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                compare_fire_predictions_v1(&logarithmic_example(), Some(sigma)),
                Err(FireComparisonErrorV1::InvalidExperimentalUncertainty { .. })
            ));
        }
        let pairs = [
            FirePredictionPairV1 {
                model: f64::MAX,
                observed: f64::MIN_POSITIVE,
            },
            FirePredictionPairV1 {
                model: f64::MIN_POSITIVE,
                observed: f64::MAX,
            },
        ];
        let extreme = compare_fire_predictions_v1(&pairs, None).unwrap();
        close(extreme.geometric_ratio, 1.0);
        assert!(extreme.combined_log_scatter.is_finite());
        assert!(extreme.combined_log_scatter > 1000.0);
    }
}
