//! Weighted cross-entropy mean updates from a completed, independent pilot.
//!
//! Fit the unit-covariance Gaussian component to an upper-tail target conditional
//! law: mu_hat = sum_i (p/q)_i 1_elite(i) z_i / sum_i (p/q)_i 1_elite(i).
//! The ranking threshold is a pilot heuristic, not a probability estimate.
//! Subsequent production must use a fresh key and must not include pilot data.
use super::{ImportanceConfig, ImportanceExecution};
use crate::product_sensitivity::Sum;
use crate::UqStatus;

/// Choices fixed before fitting a pilot; no covariance fitting is performed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImportanceFitConfig {
    /// Upper-tail fraction in (0,1]; at least two observations are selected.
    pub elite_fraction: f64,
    /// Convex weight on the fitted mean, in (0,1]; one replaces the old shift.
    pub relaxation: f64,
    /// Numerical Euclidean norm cap on the updated shift, in (0,1e6].
    pub maximum_shift_norm: f64,
}

/// Pilot diagnostics and a new proposal, never a physical probability estimate.
#[derive(Debug, Clone, PartialEq)]
pub struct ImportanceProposalFit {
    pub samples: usize,
    pub selected: usize,
    pub observed_failures: usize,
    /// Rank cutoff, or the original failure ceiling when enough failures exist.
    pub cutoff: f64,
    /// True means the selected event is strictly QoI > the original ceiling.
    /// Otherwise select QoI >= cutoff, including ALL ties at that rank.
    pub targets_failure_event: bool,
    pub effective_samples: f64,
    pub norm_before_limit: f64,
    pub norm_limited: bool,
    /// Defensive fraction is unchanged; only the Gaussian mean is learned.
    pub proposal: ImportanceConfig,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImportanceFitError {
    InvalidConfig,
    IncompletePilot,
    WorkLimit,
    Allocation,
    NumericalRange,
    Cancelled,
}
impl core::fmt::Display for ImportanceFitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(match self {
            Self::InvalidConfig => "invalid importance elite fraction, relaxation or shift limit",
            Self::IncompletePilot => "proposal fitting requires the full successful pilot",
            Self::WorkLimit => "importance pilot replay exceeds 16777216 sample-coordinate pairs",
            Self::Allocation => "importance pilot rank storage allocation refused",
            Self::NumericalRange => "importance pilot fit exceeds finite arithmetic",
            Self::Cancelled => "importance pilot fit cancelled without publishing a proposal",
        })
    }
}
impl std::error::Error for ImportanceFitError {}

impl ImportanceFitConfig {
    pub(super) fn admit(self, samples: usize, dimension: usize) -> Result<(), ImportanceFitError> {
        if !self.elite_fraction.is_finite() || self.elite_fraction <= 0.0
            || self.elite_fraction > 1.0 || !self.relaxation.is_finite()
            || self.relaxation <= 0.0 || self.relaxation > 1.0
            || !self.maximum_shift_norm.is_finite() || self.maximum_shift_norm <= 0.0
            || self.maximum_shift_norm > 1e6
        {
            return Err(ImportanceFitError::InvalidConfig);
        }
        if samples.checked_mul(dimension).is_none_or(|n| n > 16_777_216) {
            return Err(ImportanceFitError::WorkLimit);
        }
        Ok(())
    }
}

impl ImportanceExecution {
    /// Learn a mean shift from ALL of a completed pilot's accepted observations.
    /// This does not change this execution or invoke physical callbacks. Replay
    /// uses its original mixture germs, not whitened physical inputs, retaining
    /// singular/constant Gaussian maps. Prefixes and failed pilots refuse.
    ///
    /// Select at least ceil(N*fraction) upper observations (minimum two), with
    /// all cutoff ties. Once at least that many strict failures exist, fit the
    /// original failure event instead. Target/proposal weights correct the
    /// tilted pilot distribution; normalization is used ONLY for fitting, never
    /// for the final importance probability. Smooth then cap the mean norm.
    /// A single shift can miss disconnected modes or stagnate on flat responses.
    /// Defensive mass is not removed and improved variance is not guaranteed.
    ///
    /// O(N*d + N) bounded work, O(N+d) scratch, two replays of selected germs.
    /// Cancellation brackets rank selection and every replayed sample. No
    /// partial proposal is published. Choose the next experiment independently
    /// of its future outcomes; a new seed alone cannot enforce that obligation.
    pub fn fit_proposal(
        &self,
        config: ImportanceFitConfig,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ImportanceProposalFit, ImportanceFitError> {
        let n = self.values.len();
        let dim = self.plan.parameters.len();
        config.admit(self.plan.budget_max_samples, dim)?;
        if self.status != UqStatus::Complete || n != self.plan.budget_max_samples {
            return Err(ImportanceFitError::IncompletePilot);
        }
        poll(&mut cancelled)?;
        let wanted = ((n as f64 * config.elite_fraction).ceil() as usize).max(2).min(n);
        let mut scores = Vec::new();
        scores.try_reserve_exact(n).map_err(|_| ImportanceFitError::Allocation)?;
        let mut failures = 0;
        for value in &self.values {
            poll(&mut cancelled)?;
            scores.push(value.qoi);
            failures += usize::from(value.failed);
        }
        let targets_failure_event = failures >= wanted;
        let cutoff = if targets_failure_event {
            self.plan.compliance_threshold.expect("admitted threshold")
        } else {
            *scores.select_nth_unstable_by(n - wanted, f64::total_cmp).1
        };
        drop(scores);
        poll(&mut cancelled)?;
        let selected = |i: usize| {
            if targets_failure_event { self.values[i].failed }
            else { self.values[i].qoi >= cutoff }
        };
        let mut count = 0;
        let mut weight_scale = 0.0_f64;
        let mut scales = vec![0.0_f64; dim];
        for (i, value) in self.values.iter().enumerate() {
            poll(&mut cancelled)?;
            if !selected(i) { continue; }
            count += 1;
            weight_scale = weight_scale.max(value.scaled_weight);
            let germs = self.draw_germs(i).map_err(|_| ImportanceFitError::NumericalRange)?;
            for (scale, z) in scales.iter_mut().zip(germs) { *scale = scale.max(z.abs()); }
        }
        if count < 2 || weight_scale <= 0.0 {
            return Err(ImportanceFitError::NumericalRange);
        }
        let mut weights = Sum::default();
        let mut squares = Sum::default();
        let mut numerators: Vec<_> = (0..dim).map(|_| Sum::default()).collect();
        for (i, value) in self.values.iter().enumerate() {
            poll(&mut cancelled)?;
            if !selected(i) { continue; }
            // The common epsilon cancels. Scaling retains at least one weight
            // exactly equal to one without overflowing a sum of p/q weights.
            let w = value.scaled_weight / weight_scale;
            weights.add(w);
            squares.add(w * w);
            let germs = self.draw_germs(i).map_err(|_| ImportanceFitError::NumericalRange)?;
            for j in 0..dim {
                if scales[j] > 0.0 { numerators[j].add(w * (germs[j] / scales[j])); }
            }
        }
        let mut shift = Vec::with_capacity(dim);
        for j in 0..dim {
            // The normalized mean lies in [-1,1]; clamp rounding overshoot only.
            let fitted = (numerators[j].value / weights.value).clamp(-1.0, 1.0) * scales[j];
            let previous = if self.unchanged { 0.0 } else { self.config.germ_shift[j] };
            shift.push((1.0 - config.relaxation) * previous + config.relaxation * fitted);
        }
        let scale = shift.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
        let norm = if scale == 0.0 { 0.0 } else {
            fs_math::det::sqrt(shift.iter().map(|x| (x / scale).powi(2)).sum()) * scale
        };
        if !norm.is_finite() || shift.iter().any(|x| !x.is_finite()) {
            return Err(ImportanceFitError::NumericalRange);
        }
        let limited = norm > config.maximum_shift_norm;
        if limited {
            for x in &mut shift { *x *= config.maximum_shift_norm / norm; }
        }
        poll(&mut cancelled)?;
        Ok(ImportanceProposalFit {
            samples: n, selected: count, observed_failures: failures, cutoff,
            targets_failure_event, effective_samples: weights.value * weights.value / squares.value,
            norm_before_limit: norm, norm_limited: limited,
            proposal: ImportanceConfig {
                germ_shift: shift, defensive_fraction: self.config.defensive_fraction,
            },
        })
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ImportanceFitError> {
    if cancelled() { Err(ImportanceFitError::Cancelled) } else { Ok(()) }
}
