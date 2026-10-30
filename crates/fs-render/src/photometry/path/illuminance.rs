//! Incident illuminance from the same calibrated surface/volume paths as images.

use crate::photometry::area::PhotometricEmitters;
use crate::photometry::source::PhotometricSource;
use crate::photometry::{
    CalculationGrid, IlluminanceProbe, PhotometryError, Receivers, TransportOptions,
};
use crate::tracer::photometric::FixtureTransport;
use crate::tracer::{
    AdaptivePixelAccumulator, DirectStrategy, RenderExecutionConfig, Sampler, Scene, Settings,
};
use fs_exec::Cx;

use super::{PathCalculationContext, PathMedium, PointFixtureSampling, SamplingUncertainty};

/// Fixed sampling policy for surface probes, independent of image dimensions.
#[derive(Debug, Clone, Copy)]
pub struct IlluminanceSampling {
    /// Number of complete incoming paths per probe; must be positive.
    pub samples_per_probe: u32,
    /// Maximum reflecting or scattering vertices after the measurement point.
    /// One includes light reflected once or scattered once. Higher orders are
    /// retained up to this ceiling; the omitted remainder has no error bound.
    pub max_depth: u32,
    /// IID or the existing three-dimensional Owen–Sobol primary sampler.
    /// The latter samples hemisphere direction and wavelength; continuation
    /// uses the established independent counter stream. Full-path Sobol is
    /// not admitted by the calibrated transport entry point.
    pub sampler: Sampler,
    /// Seed retained with the immutable job state for deterministic replay.
    pub seed: u64,
}

impl Default for IlluminanceSampling {
    fn default() -> Self {
        Self {
            samples_per_probe: 256,
            max_depth: 8,
            sampler: Sampler::Iid,
            seed: 0,
        }
    }
}

/// Incident photometric illumination on one side of an infinitesimal probe.
///
/// No receiver BRDF or reflectance factor is applied. This is incident lux,
/// including paths reflected or scattered elsewhere in the scene, rather
/// than outgoing surface luminance. No exposure or display transform applies.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PathIlluminanceEstimate {
    /// Sum of the direct fixture and sampled continuous contributions, in lux.
    pub illuminance_lux: f64,
    /// Unreflected, unscattered point-fixture illumination, evaluated once with
    /// the existing spectral direct transport and geometric visibility.
    pub direct_lux: f64,
    /// Incoming light integrated over continuous directions, in lux. Includes
    /// directly observed finite-area emission, fixture light reflected or
    /// scattered before arrival, and physical LTE volume emission. Point-fixture
    /// delta directions are excluded here and supplied by `direct_lux`, without
    /// double counting. Area NEE and incoming paths form one combined sample.
    pub continuous_lux: f64,
    /// Number of accepted whole incoming paths in this estimate.
    pub samples: u32,
    /// Estimated standard error of the sampled component for IID samples,
    /// available after two paths. Adding the deterministic direct term does
    /// not change it. This is neither a confidence interval nor an error bound;
    /// a zero value does not certify absent model or unsampled transport error.
    pub estimated_standard_error_lux: Option<f64>,
    /// No rigorous arithmetic or spectral-quadrature enclosure is supplied.
    /// Even an estimated sampling standard error of zero leaves this unknown.
    pub numerical: fs_evidence::NumericalCertificate,
    /// Unknown total error, including finite depth, spectral discretization,
    /// geometry, source calibration and numerical arithmetic.
    pub total_error_bound_lux: Option<f64>,
}

#[derive(Clone, Copy)]
struct ProbeMoment {
    direct_lux: Option<f64>,
    continuous: AdaptivePixelAccumulator,
}

impl ProbeMoment {
    const EMPTY: Self = Self {
        direct_lux: None,
        continuous: AdaptivePixelAccumulator::EMPTY,
    };
}

/// Progressive illuminance with surface interreflection and physical media.
///
/// The direct fixture sum is deterministic. Continuous incoming radiance is
/// sampled with density `cos(theta)/pi`, so each path contributes `pi * L`.
/// The original spectral integrator supplies surface and volume interactions;
/// the probe is only a measurement point and consumes no path depth.
///
/// Probes advance in caller order for sample zero, then sample one, and so on.
/// The first direct term and whole-path moment publish together. Cancellation
/// and numerical refusals retain every previously accepted result and resume
/// at the same absolute `(probe, sample)` address with a fresh `Cx`.
/// Scene, fixtures, physical fields and explicit probe arrays remain borrowed
/// and immutable throughout the job. `Scene::camera` is unused and unvalidated.
/// Probes must lie in ambient space, including when finite emitters illuminate
/// them through native glass. An initial interior medium stack is not inferred.
pub struct PathIlluminanceJob<'scene, 'field> {
    transport: FixtureTransport<'scene, 'field>,
    receivers: Receivers<'scene>,
    sampling: IlluminanceSampling,
    moments: Vec<ProbeMoment>,
    completed: u64,
    target: u64,
}

impl<'scene, 'field> PathIlluminanceJob<'scene, 'field> {
    /// Select point-fixture sampling for reflected and scattered paths while
    /// constructing an unstarted job. The cached unscattered direct-lux sum
    /// still evaluates every fixture exactly as before.
    pub fn with_point_sampling(
        mut self,
        policy: PointFixtureSampling,
    ) -> Result<Self, PhotometryError> {
        self.set_point_sampling(policy)?;
        Ok(self)
    }

    /// Change the connection policy only before any path has committed.
    /// Refusals preserve this job and all existing estimates. Source-selection
    /// variance belongs to the whole incoming-path sample and its ordinary
    /// IID standard-error estimate; total physical error remains unknown.
    pub fn set_point_sampling(
        &mut self,
        policy: PointFixtureSampling,
    ) -> Result<(), PhotometryError> {
        if self.completed != 0 {
            return Err(PhotometryError::InvalidOptions);
        }
        self.transport.set_point_sampling(policy)
    }

    /// Retained point-fixture policy; direct probe lux is always a complete sum.
    #[must_use]
    pub fn point_sampling(&self) -> PointFixtureSampling {
        self.transport.point_sampling()
    }

    #[cfg(test)]
    pub(crate) fn disable_acceleration_for_test(&mut self) {
        assert_eq!(self.completed, 0);
        self.transport.disable_acceleration_for_test();
    }

    /// Admit a nonempty set of one-sided receiver points. Supported geometry,
    /// materials, fixtures and media match [`super::PathLuminanceJob`].
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        scene: &'scene Scene,
        sources: &'scene [PhotometricSource],
        medium: PathMedium<'scene, 'field>,
        probes: &'scene [IlluminanceProbe],
        options: TransportOptions,
        sampling: IlluminanceSampling,
        cx: &Cx<'_>,
    ) -> Result<Self, PhotometryError> {
        Self::new_with_emitters(
            scene,
            PhotometricEmitters::points_only(sources),
            medium,
            probes,
            options,
            sampling,
            cx,
        )
    }

    /// Admit probes with calibrated point fixtures and finite mesh emitters.
    /// Finite-source illumination is sampled together with the incoming path;
    /// the deterministic `direct_lux` component remains the point-fixture sum.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_emitters(
        scene: &'scene Scene,
        emitters: PhotometricEmitters<'scene>,
        medium: PathMedium<'scene, 'field>,
        probes: &'scene [IlluminanceProbe],
        options: TransportOptions,
        sampling: IlluminanceSampling,
        cx: &Cx<'_>,
    ) -> Result<Self, PhotometryError> {
        Self::begin(
            scene,
            emitters,
            medium,
            Receivers::Probes(probes),
            probes.len(),
            options,
            sampling,
            cx,
        )
    }

    /// Admit a calculation plane without allocating an intermediate probe
    /// array. Grid values are point samples, not averages over square metres.
    #[allow(clippy::too_many_arguments)]
    pub fn for_grid(
        scene: &'scene Scene,
        sources: &'scene [PhotometricSource],
        medium: PathMedium<'scene, 'field>,
        grid: CalculationGrid,
        options: TransportOptions,
        sampling: IlluminanceSampling,
        cx: &Cx<'_>,
    ) -> Result<Self, PhotometryError> {
        Self::for_grid_with_emitters(
            scene,
            PhotometricEmitters::points_only(sources),
            medium,
            grid,
            options,
            sampling,
            cx,
        )
    }

    /// Admit a calculation plane with calibrated point and finite mesh sources.
    /// As with [`Self::for_grid`], grid values are infinitesimal point samples.
    #[allow(clippy::too_many_arguments)]
    pub fn for_grid_with_emitters(
        scene: &'scene Scene,
        emitters: PhotometricEmitters<'scene>,
        medium: PathMedium<'scene, 'field>,
        grid: CalculationGrid,
        options: TransportOptions,
        sampling: IlluminanceSampling,
        cx: &Cx<'_>,
    ) -> Result<Self, PhotometryError> {
        Self::begin(
            scene,
            emitters,
            medium,
            Receivers::Grid(grid),
            grid.len(),
            options,
            sampling,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn begin(
        scene: &'scene Scene,
        emitters: PhotometricEmitters<'scene>,
        medium: PathMedium<'scene, 'field>,
        receivers: Receivers<'scene>,
        count: usize,
        options: TransportOptions,
        sampling: IlluminanceSampling,
        cx: &Cx<'_>,
    ) -> Result<Self, PhotometryError> {
        cx.checkpoint()?;
        if count == 0 || sampling.samples_per_probe == 0 || sampling.max_depth == 0 {
            return Err(PhotometryError::InvalidOptions);
        }
        if count > options.max_outputs {
            return Err(PhotometryError::OutputBudgetExceeded);
        }
        let probe_count =
            u32::try_from(count).map_err(|_| PhotometryError::OutputBudgetExceeded)?;
        let target = u64::from(probe_count)
            .checked_mul(u64::from(sampling.samples_per_probe))
            .ok_or(PhotometryError::OutputBudgetExceeded)?;
        let settings = Settings {
            width: probe_count,
            height: 1,
            spp: sampling.samples_per_probe,
            max_depth: sampling.max_depth,
            sampler: sampling.sampler,
            strategy: DirectStrategy::Mis,
            seed: sampling.seed,
        };
        let transport = FixtureTransport::try_for_probes_with_emitters(
            scene,
            emitters,
            medium,
            options,
            settings,
            probe_count,
            cx,
        )?;
        let mut moments = Vec::new();
        moments
            .try_reserve_exact(count)
            .map_err(|_| PhotometryError::AllocationRefused)?;
        moments.resize(count, ProbeMoment::EMPTY);
        cx.checkpoint()?;
        Ok(Self {
            transport,
            receivers,
            sampling,
            moments,
            completed: 0,
            target,
        })
    }

    /// Commit at most `max_path_samples` additional complete incoming paths.
    /// The exact direct calculation is cached with the first accepted path of
    /// each probe. A failed first path leaves that probe's estimate absent.
    pub fn advance(
        &mut self,
        cx: &Cx<'_>,
        max_path_samples: usize,
    ) -> Result<usize, PhotometryError> {
        cx.checkpoint()?;
        let start = self.completed;
        let limit = self
            .completed
            .saturating_add(u64::try_from(max_path_samples).unwrap_or(u64::MAX))
            .min(self.target);
        let count = self.moments.len() as u64;
        while self.completed < limit {
            cx.checkpoint()?;
            let index = (self.completed % count) as usize;
            let sample = (self.completed / count) as u32;
            let probe = match self.receivers {
                Receivers::Probes(probes) => probes[index],
                Receivers::Grid(grid) => grid.probe(index)?,
            };
            let mut candidate = self.moments[index];
            let direct_lux = match candidate.direct_lux {
                Some(value) => value,
                None => self.transport.direct_illuminance(cx, &probe)?.lux,
            };
            let xyz = self
                .transport
                .sample_probe(cx, &probe, index as u32, sample)?;
            candidate
                .continuous
                .push(xyz)
                .map_err(|_| PhotometryError::NonFiniteResult)?;
            if !(direct_lux + candidate.continuous.mean_xyz()[1]).is_finite() {
                return Err(PhotometryError::NonFiniteResult);
            }
            candidate.direct_lux = Some(direct_lux);
            cx.checkpoint()?;
            self.moments[index] = candidate;
            self.completed += 1;
        }
        Ok((self.completed - start) as usize)
    }

    /// Advance using the native parked worker crew and bounded private waves.
    ///
    /// Wave sizing, memory accounting and atomic publication match
    /// [`super::PathLuminanceJob::advance_with_execution`]. The exact direct
    /// fixture term is evaluated only at a probe's first absolute sample and
    /// committed with that wave. Incoming samples accumulate in logical order,
    /// so worker count, tile shape and mixed serial/parallel calls preserve lux
    /// and IID standard-error bits. A failed wave leaves the prior prefix intact.
    pub fn advance_with_execution(
        &mut self,
        cx: &Cx<'_>,
        max_path_samples: usize,
        execution: &RenderExecutionConfig,
    ) -> Result<usize, PhotometryError> {
        let transport = &self.transport;
        let receivers = &self.receivers;
        let retained_capacity = self.moments.capacity();
        super::execution::advance(
            cx,
            max_path_samples,
            execution,
            &mut self.moments,
            retained_capacity,
            &mut self.completed,
            self.target,
            |cx, index, sample, moment| {
                let probe = match receivers {
                    Receivers::Probes(probes) => probes[index as usize],
                    Receivers::Grid(grid) => grid.probe(index as usize)?,
                };
                let direct_lux = if sample == 0 && moment.direct_lux.is_none() {
                    Some(transport.direct_illuminance(cx, &probe)?.lux)
                } else {
                    None
                };
                Ok(super::execution::PathSample {
                    xyz: transport.sample_probe(cx, &probe, index, sample)?,
                    direct_lux,
                })
            },
            |moment, sample| {
                let direct_lux = moment
                    .direct_lux
                    .or(sample.direct_lux)
                    .ok_or(PhotometryError::NonFiniteResult)?;
                moment
                    .continuous
                    .push(sample.xyz)
                    .map_err(|_| PhotometryError::NonFiniteResult)?;
                if !(direct_lux + moment.continuous.mean_xyz()[1]).is_finite() {
                    return Err(PhotometryError::NonFiniteResult);
                }
                moment.direct_lux = Some(direct_lux);
                Ok(())
            },
        )
    }

    /// Read an accepted probe estimate. Uncomputed and out-of-range probes
    /// return `None`; they are not presented as zero illumination.
    #[must_use]
    pub fn probe_estimate(&self, index: usize) -> Option<PathIlluminanceEstimate> {
        let moment = *self.moments.get(index)?;
        let direct_lux = moment.direct_lux?;
        let samples = moment.continuous.sample_count();
        let continuous_lux = moment.continuous.mean_xyz()[1];
        let estimated_standard_error_lux = if self.sampling.sampler == Sampler::Iid && samples >= 2
        {
            Some(moment.continuous.dispersion_proxy_xyz()[1])
        } else {
            None
        };
        Some(PathIlluminanceEstimate {
            illuminance_lux: direct_lux + continuous_lux,
            direct_lux,
            continuous_lux,
            samples,
            estimated_standard_error_lux,
            numerical: fs_evidence::NumericalCertificate::no_claim(),
            total_error_bound_lux: None,
        })
    }

    /// Interpret this probe's sampling uncertainty in lux, using exactly the
    /// standard error exposed by [`Self::probe_estimate`]. The deterministic
    /// direct term adds no sampling variance. Uncomputed or out-of-range
    /// probes return `None`; one IID sample and correlated samples are distinct.
    #[must_use]
    pub fn sampling_uncertainty_lux(&self, index: usize) -> Option<SamplingUncertainty> {
        let estimate = self.probe_estimate(index)?;
        Some(SamplingUncertainty::from_estimate(
            self.sampling.sampler,
            estimate.samples,
            estimate.estimated_standard_error_lux,
        ))
    }

    /// Declared observer, model versions and actual retained sampling settings.
    /// This context is not a calibration or a complete scene identity.
    #[must_use]
    pub const fn calculation_context(&self) -> PathCalculationContext {
        PathCalculationContext::new(
            self.sampling.sampler,
            self.sampling.seed,
            self.sampling.max_depth,
        )
    }

    /// Probe-order estimates without another result allocation.
    pub fn estimates(&self) -> impl ExactSizeIterator<Item = Option<PathIlluminanceEstimate>> + '_ {
        (0..self.moments.len()).map(|index| self.probe_estimate(index))
    }

    /// Number of retained receiver estimates.
    #[must_use]
    pub fn probe_count(&self) -> usize {
        self.moments.len()
    }

    /// Fixed sample ceiling for each probe.
    #[must_use]
    pub const fn samples_per_probe(&self) -> u32 {
        self.sampling.samples_per_probe
    }

    /// Number of successfully committed incoming paths across all probes.
    #[must_use]
    pub const fn completed_path_samples(&self) -> u64 {
        self.completed
    }

    /// Whether every probe reached its declared sample ceiling.
    #[must_use]
    pub const fn is_complete(&self) -> bool {
        self.completed == self.target
    }
}
