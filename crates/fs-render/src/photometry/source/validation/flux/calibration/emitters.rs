//! Independently bounded emitter gains from controlled lux observations.
//!
//! Row r asserts M_r = B_r + sum_i(d_ri * E_ri * g_i), with known relative
//! outputs d in [0,1] and fixed nonnegative gains g. E comes from the ORIGINAL
//! prepared direct transport, including spectral media and mesh visibility.
//! Beam/spectral shape, placement and medium must stay fixed as outputs change.
//! This is not a DMX response curve, a leakage model or a thermal recalibration.
//!
//! Sparse coefficient rows are traced once. Repeated interval contraction uses
//! positivity and the existing single-gain inverse, never normal equations or
//! an assumed noise distribution. Every update retains every compatible gain
//! vector. Shared input dependencies are relaxed, so nonempty boxes do NOT prove
//! feasibility, uniqueness or physical accuracy. Indistinguishable emitters can
//! remain broad; dark experiments remain unbounded. A finite sweep limit is an
//! explicit termination, not a convergence claim. Total physical error stays
//! unknown. Measurement ranges have the same meaning as in the parent module.

use super::{
    ContentHash, Cx, DomainHasher, EnclosedIlluminance, FluxCalibrationError,
    FluxScaleBound, GridEnclosureBudget, IlluminanceMeasurement, Interval,
    MeasurementCompatibility, PhotometricScene, PhotometryError, Prepared,
    compatible_gain, finite, hash_float, model_range, nonnegative, prepare,
};
use super::super::super::evaluate_with_contributions;

/// Version of controlled-output admission, contraction and result interpretation.
pub const EMITTER_CALIBRATION_SEMANTICS_VERSION: u32 = 1;

/// One meter observation and one known relative output for EACH original source.
///
/// Outputs multiply retained emitted lumens after original dimming/filtering.
/// The job validates length, finiteness and [0,1] range. Exact zero means off
/// only within this fixed-shape multiplicative model; additive leakage is not
/// included. The measurement evidence must justify these operating conditions.
#[derive(Debug, Clone, Copy)]
pub struct EmitterCalibrationObservation<'a> {
    /// Original SI receiver, bounded measured lux and bounded background.
    pub measurement: IlluminanceMeasurement,
    /// Source-order controlled outputs; slices may be shared between rows.
    pub relative_output: &'a [f64],
}

/// Logical storage/admission and total refinement limits. Native geometry,
/// medium and output budgets also remain binding. Allocator metadata is excluded.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EmitterCalibrationBudget {
    /// Maximum original emitters, also bounding each reusable scratch array.
    pub max_sources: usize,
    /// Maximum meter observations and cached row offsets.
    pub max_measurements: usize,
    /// Maximum retained nonzero possible source/observation coefficients.
    pub max_coefficients: usize,
    /// Positive maximum full contraction sweeps after row preparation.
    pub max_sweeps: u32,
}

impl Default for EmitterCalibrationBudget {
    fn default() -> Self {
        Self { max_sources: 1024, max_measurements: 4096,
            max_coefficients: 262_144, max_sweeps: 32 }
    }
}

/// Why inverse iteration stopped. None of these establishes measured accuracy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EmitterCalibrationTermination {
    /// A full sweep made no binary64 endpoint changes; not a feasibility proof.
    Stable,
    /// The declared sweep cap was reached. The retained outer box is still
    /// sound, but further contraction might narrow it or detect a contradiction.
    SweepLimit,
    /// This row proved the current constraints inconsistent. No gains export.
    Inconsistent {
        /// Zero-based row that first proved exclusion; not a minimal conflict set.
        measurement: usize,
    },
}

/// Completed conditional outer gain box, borrowing the immutable job result.
#[derive(Debug, Clone, Copy)]
pub struct EmitterCalibrationResult<'a> {
    scales: Option<&'a [Interval]>,
    /// Distinguishes a stable enclosure, an exhausted refinement cap and exclusion.
    pub termination: EmitterCalibrationTermination,
    /// Full sweeps completed, excluding a partially completed inconsistent sweep.
    pub completed_sweeps: u32,
    /// Binds the ordered emitters, observations, drives, coefficients and outcome.
    pub evidence_identity: ContentHash,
}

impl<'a> EmitterCalibrationResult<'a> {
    /// Source-order necessary gain ranges. Infinite upper endpoints stay explicit;
    /// `None` means inconsistent, not an arbitrary finite prior or a zero vector.
    #[must_use]
    pub const fn scale_enclosures(self) -> Option<&'a [Interval]> { self.scales }

    /// Export one finite outer range to the original forward flux APIs. These
    /// marginal ranges lose joint correlations, conservatively widening predictions.
    /// Sweep-limited ranges remain valid necessary bounds, not converged estimates.
    #[must_use]
    pub fn flux_scale_bound(self, source: usize) -> Option<FluxScaleBound> {
        let scale = *self.scales?.get(source)?;
        let mut hash = DomainHasher::new("org.frankensim.photometry.emitter-calibration.source.v1");
        hash.update(&self.evidence_identity.0);
        hash.update(&(source as u64).to_le_bytes());
        FluxScaleBound::try_new(scale.lo(), scale.hi(), hash.finalize()).ok()
    }
}

/// Holdout prediction at its declared drive levels, without refitting the gains.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EmitterCalibrationCheck {
    /// Native numerical result at the scene's ORIGINAL reference output levels,
    /// not the driven observation. Its unknown total-accuracy claim is unchanged.
    pub reference: EnclosedIlluminance,
    /// Gain- and drive-conditioned direct lux plus declared background.
    pub predicted_measured_lux: Interval,
    /// Original measurement interval, in lux.
    pub measured_lux: Interval,
    /// Overlap is non-exclusion, never validation or statistical coverage.
    pub compatibility: MeasurementCompatibility,
    /// Completed inverse result identity, unchanged by holdout evaluation.
    pub calibration_identity: ContentHash,
    /// Binds this measurement evidence, drive vector and native numerical model.
    pub observation_identity: ContentHash,
}

#[derive(Debug, Clone, Copy)]
struct Coefficient { source: usize, lux: Interval }

/// Progressive independent-emitter calibration. No source or medium is copied.
///
/// Each observation uses the original nominal-plus-enclosure evaluation once.
/// Its bounded sparse coefficient row is cached. Refinement reuses
/// those rows and O(source-count) scratch, with linear work in each row's stored
/// entries. It never subtracts an infinite total from itself. The native nominal
/// evaluation is unchanged. `advance` counts prepared OR contracted rows, so
/// chunking cannot change arithmetic, source addresses or accepted row ordering.
/// Cancellation and allocation/budget refusals preserve every accepted row.
pub struct EmitterCalibrationJob<'a, 'scene> {
    scene: &'a PhotometricScene<'scene>,
    observations: &'a [EmitterCalibrationObservation<'a>],
    prepared: Prepared<'scene>,
    budget: EmitterCalibrationBudget,
    coefficients: Vec<Coefficient>,
    offsets: Vec<usize>,
    pending: Vec<Coefficient>,
    scales: Vec<Interval>,
    candidate: Vec<Interval>,
    suffix: Vec<Interval>,
    next_row: usize,
    sweeps: u32,
    sweep_changed: bool,
    termination: Option<EmitterCalibrationTermination>,
    identity: ContentHash,
}

impl<'a, 'scene> EmitterCalibrationJob<'a, 'scene> {
    /// Admit immutable observations and explicitly identified fixed source states.
    /// Start every gain at [0,+infinity]; no finite prior or preferred solution is
    /// invented. Drive changes require fixed angular/spectral shape and no leakage.
    pub fn new(
        scene: &'a PhotometricScene<'scene>,
        observations: &'a [EmitterCalibrationObservation<'a>],
        scope_evidence: ContentHash, budget: EmitterCalibrationBudget, cx: &Cx<'_>,
    ) -> Result<Self, FluxCalibrationError> {
        cx.checkpoint()?;
        let sources = scene.sources.len();
        if sources == 0 || observations.is_empty() || budget.max_sweeps == 0 {
            return Err(PhotometryError::InvalidOptions.into());
        }
        if sources > budget.max_sources || observations.len() > budget.max_measurements
            || observations.len() > scene.options.max_outputs {
            return Err(FluxCalibrationError::EmitterBudgetExceeded);
        }
        if scope_evidence == ContentHash([0; 32]) { return Err(FluxCalibrationError::MissingEvidence); }
        for row in observations { validate_output(row.relative_output, sources, cx)?; }
        let count = sources.checked_add(1).ok_or(FluxCalibrationError::EmitterBudgetExceeded)?;
        let mut offsets = Vec::new();
        offsets.try_reserve_exact(observations.len().checked_add(1)
            .ok_or(FluxCalibrationError::EmitterBudgetExceeded)?)
            .map_err(|_| PhotometryError::AllocationRefused)?;
        offsets.push(0);
        let mut pending = Vec::new();
        pending.try_reserve_exact(sources).map_err(|_| PhotometryError::AllocationRefused)?;
        let scales = filled(sources, Interval::new(0.0, f64::INFINITY))?;
        let candidate = filled(sources, Interval::point(0.0))?;
        let suffix = filled(count, Interval::point(0.0))?;
        let prepared = prepare(scene, cx)?;
        let mut hash = DomainHasher::new("org.frankensim.photometry.emitter-calibration.v1");
        hash.update(&EMITTER_CALIBRATION_SEMANTICS_VERSION.to_le_bytes());
        hash.update(&scope_evidence.0);
        hash.update(&(sources as u64).to_le_bytes());
        hash.update(&(observations.len() as u64).to_le_bytes());
        cx.checkpoint()?;
        Ok(Self { scene, observations, prepared, budget, coefficients: Vec::new(), offsets,
            pending, scales, candidate, suffix, next_row: 0, sweeps: 0, sweep_changed: false,
            termination: None, identity: hash.finalize() })
    }

    /// Adjust optional medium work only before any coefficient row commits.
    /// Native limits remain binding; unavailable rows are never dropped as outliers.
    pub fn set_grid_enclosure_budget(&mut self, budget: GridEnclosureBudget) -> Result<(), FluxCalibrationError> {
        if self.prepared_measurements() != 0 || budget.max_nodes_per_segment == 0 {
            return Err(PhotometryError::InvalidOptions.into());
        }
        if let Prepared::Available { medium, .. } = &mut self.prepared { medium.set_grid_budget(budget); }
        Ok(())
    }

    /// Perform at most `max_rows` complete preparation or contraction steps.
    /// No partially prepared row or partially contracted gain vector is published.
    pub fn advance(&mut self, cx: &Cx<'_>, max_rows: usize) -> Result<usize, FluxCalibrationError> {
        cx.checkpoint()?;
        let mut done = 0;
        while done < max_rows && self.termination.is_none() {
            if self.prepared_measurements() < self.observations.len() {
                self.prepare_row(cx)?;
            } else {
                let row = self.next_row;
                let entries = &self.coefficients[self.offsets[row]..self.offsets[row + 1]];
                let changed = contract_row(entries, self.observations[row].measurement,
                    &self.scales, &mut self.candidate, &mut self.suffix, cx)?;
                cx.checkpoint()?;
                match changed {
                    None => self.termination = Some(EmitterCalibrationTermination::Inconsistent { measurement: row }),
                    Some(changed) => {
                        for entry in entries { self.scales[entry.source] = self.candidate[entry.source]; }
                        self.sweep_changed |= changed;
                        self.next_row += 1;
                        if self.next_row == self.observations.len() {
                            self.next_row = 0;
                            self.sweeps += 1;
                            if !self.sweep_changed {
                                self.termination = Some(EmitterCalibrationTermination::Stable);
                            } else if self.sweeps == self.budget.max_sweeps {
                                self.termination = Some(EmitterCalibrationTermination::SweepLimit);
                            }
                            self.sweep_changed = false;
                        }
                    }
                }
            }
            done += 1;
        }
        Ok(done)
    }

    fn prepare_row(&mut self, cx: &Cx<'_>) -> Result<(), FluxCalibrationError> {
        let index = self.prepared_measurements();
        let observation = self.observations[index];
        self.pending.clear();
        let pending = &mut self.pending;
        let reference = evaluate_with_contributions(self.scene, cx,
            &observation.measurement.probe(), &self.prepared, &mut |source, value| {
                let lux = driven(value, observation.relative_output[source]);
                if lux.hi() > 0.0 { pending.push(Coefficient { source, lux }); }
            })?;
        let (_, model) = model_range(reference, index)?;
        let required = self.coefficients.len().checked_add(self.pending.len())
            .ok_or(FluxCalibrationError::EmitterBudgetExceeded)?;
        if required > self.budget.max_coefficients { return Err(FluxCalibrationError::EmitterBudgetExceeded); }
        if self.coefficients.capacity() < required {
            let target = required.max(self.coefficients.capacity().saturating_mul(2))
                .min(self.budget.max_coefficients);
            self.coefficients.try_reserve_exact(target - self.coefficients.len())
                .map_err(|_| PhotometryError::AllocationRefused)?;
        }
        let row_identity = observation_identity(observation, model, cx)?;
        let mut hash = DomainHasher::new("org.frankensim.photometry.emitter-calibration.row.v1");
        hash.update(&self.identity.0); hash.update(&row_identity.0);
        hash.update(&(self.pending.len() as u64).to_le_bytes());
        for (i, entry) in self.pending.iter().enumerate() {
            if i % 128 == 0 { cx.checkpoint()?; }
            hash.update(&(entry.source as u64).to_le_bytes());
            hash_float(&mut hash, entry.lux.lo()); hash_float(&mut hash, entry.lux.hi());
        }
        let identity = hash.finalize();
        cx.checkpoint()?;
        self.coefficients.extend_from_slice(&self.pending);
        self.offsets.push(required);
        self.identity = identity;
        Ok(())
    }

    /// Complete cached measurement rows; no retracing is needed during refinement.
    #[must_use]
    pub fn prepared_measurements(&self) -> usize { self.offsets.len() - 1 }
    /// Retained nonzero possible coefficients, excluding reusable scratch capacity.
    #[must_use]
    pub fn coefficient_count(&self) -> usize { self.coefficients.len() }
    /// Explicit stable, sweep-limited or inconsistent terminal outcome is available.
    #[must_use]
    pub fn is_complete(&self) -> bool { self.termination.is_some() }

    /// Export only a completed outer box. Hashing is bounded by admitted sources.
    /// Unbounded components remain infinite and cannot become finite flux bounds.
    #[must_use]
    pub fn result(&self) -> Option<EmitterCalibrationResult<'_>> {
        let termination = self.termination?;
        let mut hash = DomainHasher::new("org.frankensim.photometry.emitter-calibration.result.v1");
        hash.update(&self.identity.0); hash.update(&self.sweeps.to_le_bytes());
        let scales = match termination {
            EmitterCalibrationTermination::Inconsistent { measurement } => {
                hash.update(&[2]); hash.update(&(measurement as u64).to_le_bytes()); None
            }
            _ => {
                hash.update(&[u8::from(termination == EmitterCalibrationTermination::SweepLimit)]);
                for scale in &self.scales { hash_float(&mut hash, scale.lo()); hash_float(&mut hash, scale.hi()); }
                Some(self.scales.as_slice())
            }
        };
        Some(EmitterCalibrationResult { scales, termination,
            completed_sweeps: self.sweeps, evidence_identity: hash.finalize() })
    }

    /// Check additional data at declared outputs without changing or refitting
    /// any gain. Active unbounded gains may make a finite prediction unavailable.
    pub fn check_observation(
        &self, cx: &Cx<'_>, observation: EmitterCalibrationObservation<'_>,
    ) -> Result<EmitterCalibrationCheck, FluxCalibrationError> {
        cx.checkpoint()?;
        validate_output(observation.relative_output, self.scales.len(), cx)?;
        let result = self.result().ok_or(FluxCalibrationError::NoFiniteFit)?;
        let scales = result.scale_enclosures().ok_or(FluxCalibrationError::NoFiniteFit)?;
        let mut predicted = observation.measurement.background_lux();
        let reference = evaluate_with_contributions(self.scene, cx,
            &observation.measurement.probe(), &self.prepared, &mut |source, value| {
                predicted = positive_sum(predicted, positive_product(
                    driven(value, observation.relative_output[source]), scales[source]));
            })?;
        let (_, model) = model_range(reference, self.prepared_measurements())?;
        if !finite(predicted) { return Err(FluxCalibrationError::PredictionRange); }
        let measured = observation.measurement.measured_lux();
        let compatibility = if predicted.intersect(measured).is_some() {
            MeasurementCompatibility::NotExcluded
        } else { MeasurementCompatibility::Contradicted };
        let observation_identity = observation_identity(observation, model, cx)?;
        cx.checkpoint()?;
        Ok(EmitterCalibrationCheck { reference, predicted_measured_lux: predicted,
            measured_lux: measured, compatibility, calibration_identity: result.evidence_identity,
            observation_identity })
    }
}

fn validate_output(values: &[f64], sources: usize, cx: &Cx<'_>) -> Result<(), FluxCalibrationError> {
    if values.len() != sources { return Err(PhotometryError::InvalidOptions.into()); }
    for chunk in values.chunks(128) {
        cx.checkpoint()?;
        if chunk.iter().any(|v| !v.is_finite() || !(0.0..=1.0).contains(v)) {
            return Err(PhotometryError::InvalidOptions.into());
        }
    }
    Ok(())
}

fn observation_identity(row: EmitterCalibrationObservation<'_>, model: ContentHash, cx: &Cx<'_>) -> Result<ContentHash, FluxCalibrationError> {
    let mut hash = DomainHasher::new("org.frankensim.photometry.emitter-observation.v1");
    hash.update(&model.0); hash.update(&row.measurement.evidence_identity().0);
    for range in [row.measurement.measured_lux(), row.measurement.background_lux()] {
        hash_float(&mut hash, range.lo()); hash_float(&mut hash, range.hi());
    }
    hash.update(&(row.relative_output.len() as u64).to_le_bytes());
    for chunk in row.relative_output.chunks(128) {
        cx.checkpoint()?;
        for &value in chunk { hash_float(&mut hash, value); }
    }
    Ok(hash.finalize())
}

fn filled(count: usize, value: Interval) -> Result<Vec<Interval>, FluxCalibrationError> {
    let mut values = Vec::new();
    values.try_reserve_exact(count).map_err(|_| PhotometryError::AllocationRefused)?;
    values.resize(count, value);
    Ok(values)
}

fn driven(value: Interval, output: f64) -> Interval {
    if output == 0.0 || value.hi() == 0.0 { Interval::point(0.0) }
    else if output == 1.0 { value }
    else { nonnegative(value * Interval::point(output)) }
}

// Nonnegative arithmetic needs exact zero handling: an unknown nonnegative
// gain times an identically zero coefficient is ZERO, never 0*infinity/WHOLE.
fn positive_product(a: Interval, b: Interval) -> Interval {
    if a.hi() == 0.0 || b.hi() == 0.0 { return Interval::point(0.0); }
    let lo = if a.lo() == 0.0 || b.lo() == 0.0 { 0.0 }
        else { (Interval::point(a.lo()) * Interval::point(b.lo())).lo().max(0.0) };
    let hi = if a.hi().is_infinite() || b.hi().is_infinite() { f64::INFINITY }
        else { (Interval::point(a.hi()) * Interval::point(b.hi())).hi() };
    Interval::new(lo, hi)
}

fn positive_sum(a: Interval, b: Interval) -> Interval {
    if a.hi() == 0.0 { b } else if b.hi() == 0.0 { a } else { nonnegative(a + b) }
}

// Interval Gauss-Seidel for one positive mixture row. Later terms use a suffix
// of the OLD box; earlier terms use the tightened candidate. Both contain the
// same true vector. Never use `total - own_term`, especially with unbounded
// gains. One sweep is O(nonzero coefficients), not quadratic in source count.
fn contract_row(
    entries: &[Coefficient], measurement: IlluminanceMeasurement, scales: &[Interval],
    candidate: &mut [Interval], suffix: &mut [Interval], cx: &Cx<'_>,
) -> Result<Option<bool>, FluxCalibrationError> {
    suffix[entries.len()] = Interval::point(0.0);
    for (i, entry) in entries.iter().enumerate().rev() {
        if i % 128 == 0 { cx.checkpoint()?; }
        suffix[i] = positive_sum(positive_product(entry.lux, scales[entry.source]), suffix[i + 1]);
    }
    let mut prefix = measurement.background_lux();
    let measured = measurement.measured_lux();
    if positive_sum(prefix, suffix[0]).intersect(measured).is_none() { return Ok(None); }
    let mut changed = false;
    for (i, entry) in entries.iter().enumerate() {
        if i % 128 == 0 { cx.checkpoint()?; }
        let others = positive_sum(prefix, suffix[i + 1]);
        let Some(range) = compatible_gain(entry.lux, measured, others)
            .and_then(|range| range.intersect(scales[entry.source])) else { return Ok(None); };
        changed |= range != scales[entry.source];
        candidate[entry.source] = range;
        prefix = positive_sum(prefix, positive_product(entry.lux, range));
    }
    cx.checkpoint()?;
    Ok(Some(changed))
}

#[cfg(test)]
mod tests;
