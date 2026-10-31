//! Interpretation of accepted photometric samples and their declared settings.

use crate::tracer::Sampler;

/// Sampling uncertainty of an accepted path estimate, in the units named by
/// the measurement accessor. This does not include numerical, calibration,
/// spectral, geometry/material/medium or omitted-transport error.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SamplingUncertainty {
    /// Fewer than two IID paths cannot supply a sample variance.
    InsufficientIidSamples {
        /// Number of accepted whole paths.
        samples: u32,
    },
    /// Sample-based standard error of the mean, not a confidence interval or
    /// error bound. Zero does not rule out unsampled paths or model error.
    EstimatedIidStandardError {
        /// Standard error in the measurement accessor's stated SI units.
        standard_error: f64,
        /// Number of accepted whole paths.
        samples: u32,
    },
    /// Correlated primary samples do not justify the IID variance formula.
    CorrelatedSamples {
        /// Number of accepted whole paths; increasing it does not create an
        /// IID standard-error claim for this sampler.
        samples: u32,
    },
}

impl SamplingUncertainty {
    pub(super) fn from_estimate(
        sampler: Sampler,
        samples: u32,
        standard_error: Option<f64>,
    ) -> Self {
        match (sampler, standard_error) {
            (Sampler::Iid, Some(standard_error)) => Self::EstimatedIidStandardError {
                standard_error,
                samples,
            },
            (Sampler::Iid, None) => Self::InsufficientIidSamples { samples },
            _ => Self::CorrelatedSamples { samples },
        }
    }
}

/// Declared calculation settings copied from an admitted path job.
///
/// This is neither a complete scene identity nor an accuracy certificate.
/// Replay also requires the original scene, fixtures, medium, measurement
/// geometry and remaining options. Calibration, spectral discretization,
/// geometry/material/medium errors remain unbounded. The separate optional
/// [`super::DiffuseVacuumRemainder`] addresses only a restricted finite-depth
/// model remainder; these settings alone establish no bound.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PathCalculationContext {
    /// Actual primary-sample policy of this job.
    pub sampler: Sampler,
    /// Declared seed used by its counter and scrambled primary streams.
    pub seed: u64,
    /// Declared native path vertex ceiling, not a bound on the omitted tail.
    pub max_depth: u32,
    /// CIE source-data identity and spectral sampling convention.
    pub observer_id: &'static str,
    /// Version of the shared direct photometric model.
    pub photometry_semantics_version: u32,
    /// Version of fixture-path admission, transport and accumulation.
    pub path_semantics_version: u32,
}

impl PathCalculationContext {
    pub(super) const fn new(sampler: Sampler, seed: u64, max_depth: u32) -> Self {
        Self {
            sampler,
            seed,
            max_depth,
            observer_id: crate::photometry::observer::OBSERVER_ID,
            photometry_semantics_version: crate::photometry::PHOTOMETRY_SEMANTICS_VERSION,
            path_semantics_version: super::PHOTOMETRIC_PATH_SEMANTICS_VERSION,
        }
    }
}
