//! Infer a common fixture-flux gain from bounded lux observations.
//!
//! For each measurement, the declared model is `M = g E + B`: `E` is the
//! original direct-transport result, `B` is an additive background interval,
//! and ONE nonnegative gain `g` scales every emitter in the scene. This supports
//! a multi-emitter fixture, not independently varying emitter calibrations.
//! Scene geometry, beam shape, spectral shape and the medium remain fixed.
//!
//! Intersect outward necessary gain ranges, without least-squares fitting,
//! Gaussian assumptions, independence assumptions or outlier removal. An empty
//! intersection disproves the common-gain model under the declared intervals.
//! A nonempty intersection only fails to exclude it: interval relaxation loses
//! shared numerical/input dependencies and does not prove physical agreement.
//!
//! Measurement/background ranges must be declared deterministic bounds. Meter
//! standard errors and expanded-uncertainty coverage intervals are not such
//! guarantees (NIST TN 1297, section 6:
//! <https://www.nist.gov/pml/nist-technical-note-1297/nist-tn-1297-6-expanded-uncertainty>).
//! Supplied evidence records must state
//! their meaning, operating conditions and any unbounded systematic effects.
//! Background must include any light omitted by the direct model, or the fit is
//! conditional on its absence. No measured data or calibration authority is
//! supplied here, and total physical error remains unknown.

use super::{FluxScaleBound, GridEnclosureBudget};
use super::super::{
    ContentHash, Cx, DomainHasher, EnclosedIlluminance, EnclosureUnavailable,
    IlluminanceProbe, Interval, PhotometricScene, PhotometryError, Prepared,
    evaluate, finite, hash_float, nonnegative, prepare,
};

/// Independent emitter gains from controlled-output observations.
pub mod emitters;

/// Version of measurement admission, interval inversion and evidence binding.
pub const FLUX_CALIBRATION_SEMANTICS_VERSION: u32 = 1;

/// A refused declaration, unavailable numerical model, or interrupted operation.
#[derive(Debug)]
pub enum FluxCalibrationError {
    /// Lux intervals must be finite and nonnegative.
    InvalidMeasurement,
    /// Independent-emitter count or sparse coefficient storage exceeds its budget.
    EmitterBudgetExceeded,
    /// A nonzero identity must describe the measurement or calibration scope.
    MissingEvidence,
    /// No finite completed fit is available for a prediction check.
    NoFiniteFit,
    /// The conditional prediction cannot be represented by finite endpoints.
    PredictionRange,
    /// A required observation cannot be reduced to an enclosed model value.
    NumericalUnavailable {
        /// Calibration row index; separate checks use the completed-row count.
        measurement: usize,
        /// Original optional-enclosure refusal; no nominal value is certified.
        reason: EnclosureUnavailable,
    },
    /// Original transport, budget or cancellation refusal.
    Transport(PhotometryError),
}

impl core::fmt::Display for FluxCalibrationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "photometric flux calibration refused: {self:?}")
    }
}
impl std::error::Error for FluxCalibrationError {}
impl From<fs_exec::Cancelled> for FluxCalibrationError {
    fn from(_: fs_exec::Cancelled) -> Self { Self::Transport(PhotometryError::Cancelled) }
}
impl From<PhotometryError> for FluxCalibrationError {
    fn from(error: PhotometryError) -> Self { Self::Transport(error) }
}

/// One source-state observation at an original SI probe, with bounded background.
/// No uncertainty is silently inferred from a point estimate or a meter reading.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IlluminanceMeasurement {
    probe: IlluminanceProbe,
    measured: Interval,
    background: Interval,
    evidence: ContentHash,
}

impl IlluminanceMeasurement {
    /// Admit deterministic ranges in lux. The evidence must identify the actual
    /// reading, meter assumptions and background/operating-state declaration.
    pub fn try_new(
        probe: IlluminanceProbe, measured_lux: Interval, background_lux: Interval,
        evidence: ContentHash,
    ) -> Result<Self, FluxCalibrationError> {
        if !finite(measured_lux) || !finite(background_lux)
            || measured_lux.lo() < 0.0 || background_lux.lo() < 0.0 {
            return Err(FluxCalibrationError::InvalidMeasurement);
        }
        if evidence == ContentHash([0; 32]) { return Err(FluxCalibrationError::MissingEvidence); }
        Ok(Self { probe, measured: measured_lux, background: background_lux, evidence })
    }

    /// Original receiver position and normalized receiving normal.
    #[must_use]
    pub const fn probe(self) -> IlluminanceProbe { self.probe }
    /// Declared total incident illuminance range, in lux.
    #[must_use]
    pub const fn measured_lux(self) -> Interval { self.measured }
    /// Declared additive incident background, in lux; not a source gain.
    #[must_use]
    pub const fn background_lux(self) -> Interval { self.background }
    /// Original measurement/assumption record, not authenticated by this API.
    #[must_use]
    pub const fn evidence_identity(self) -> ContentHash { self.evidence }
}

/// Whether the declared observations exclude one common nonnegative gain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FluxCalibrationStatus {
    /// A finite outer gain range remains; this is not proof of a physical fit.
    NotExcluded,
    /// Data leave the gain unbounded; no finite calibration can be exported.
    Unbounded,
    /// This observation contradicts the previous constraints (or the model alone).
    Inconsistent {
        /// Zero-based first conflicting observation in the caller's order.
        measurement: usize,
    },
}

/// A completed fit or an early, conclusive inconsistency. Never a partial fit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FluxCalibrationResult {
    scale: Option<Interval>,
    status: FluxCalibrationStatus,
    processed: usize,
    identity: ContentHash,
}

impl FluxCalibrationResult {
    /// Remaining necessary gain range. `None` means inconsistent, not unknown;
    /// an infinite upper endpoint is explicitly underconstrained.
    #[must_use]
    pub const fn scale_enclosure(self) -> Option<Interval> { self.scale }
    /// Distinguish exclusion, non-exclusion and an unbounded inverse problem.
    #[must_use]
    pub const fn status(self) -> FluxCalibrationStatus { self.status }
    /// Number of complete observations consumed, including a conflicting one.
    #[must_use]
    pub const fn processed_measurements(self) -> usize { self.processed }
    /// Binds scope, ordered observations, numerical models and inferred result.
    #[must_use]
    pub const fn evidence_identity(self) -> ContentHash { self.identity }

    /// Export a finite conditional gain to the existing forward flux API.
    /// Apply the SAME bound to every emitter calibrated together. Reuse at other
    /// fixture settings, spectra or operating states requires new justification.
    /// Empty/unbounded fits never become an arbitrary finite gain or exact one.
    #[must_use]
    pub fn flux_scale_bound(self) -> Option<FluxScaleBound> {
        let scale = self.scale?;
        FluxScaleBound::try_new(scale.lo(), scale.hi(), self.identity).ok()
    }
}

/// Comparison against an additional observation; never a refit or a pass label.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MeasurementCompatibility {
    /// The declared measurement and conditional prediction are disjoint.
    Contradicted,
    /// Overlap alone does not verify fixture accuracy or statistical coverage.
    NotExcluded,
}

/// An original numerical result with a separate conditional measurement check.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FluxCalibrationCheck {
    /// Original unscaled numerical result and unchanged unknown total error.
    pub nominal: EnclosedIlluminance,
    /// Gain-scaled model plus this measurement's additive background, in lux.
    pub predicted_measured_lux: Interval,
    /// Original declared measurement range, in lux.
    pub measured_lux: Interval,
    /// No physical-accuracy or statistical-confidence promotion on overlap.
    pub compatibility: MeasurementCompatibility,
    /// Completed calibration identity used for prediction; unchanged by checks.
    pub calibration_identity: ContentHash,
    /// Additional observation's caller-supplied evidence identity.
    pub measurement_identity: ContentHash,
}

/// Progressive common-gain calibration using the original prepared transport.
///
/// Geometry and spectral metadata are prepared once. Observations stay borrowed;
/// apart from native preparation, retained state and extra work per row are O(1).
/// There is no per-row source-sized array or independent renderer. Work is bounded
/// by `max_outputs` and each `advance` limit. Cancellation/refusal preserves the
/// accepted prefix; a missing numerical enclosure is never skipped as an outlier.
pub struct FluxCalibrationJob<'a, 'scene> {
    scene: &'a PhotometricScene<'scene>,
    measurements: &'a [IlluminanceMeasurement],
    prepared: Prepared<'scene>,
    scale: Option<Interval>,
    conflict: Option<usize>,
    processed: usize,
    identity: ContentHash,
}

impl<'a, 'scene> FluxCalibrationJob<'a, 'scene> {
    /// Scope evidence must declare which fixed-state emitters share the gain.
    /// Empty observations refuse; no finite prior is manufactured to hide an
    /// unidentifiable dark or fully occluded calibration experiment.
    pub fn new(
        scene: &'a PhotometricScene<'scene>, measurements: &'a [IlluminanceMeasurement],
        scope_evidence: ContentHash, cx: &Cx<'_>,
    ) -> Result<Self, FluxCalibrationError> {
        cx.checkpoint()?;
        if measurements.is_empty() { return Err(PhotometryError::InvalidOptions.into()); }
        if measurements.len() > scene.options.max_outputs {
            return Err(PhotometryError::OutputBudgetExceeded.into());
        }
        if scope_evidence == ContentHash([0; 32]) { return Err(FluxCalibrationError::MissingEvidence); }
        let prepared = prepare(scene, cx)?;
        let mut hash = DomainHasher::new("org.frankensim.photometry.flux-calibration.v1");
        hash.update(&FLUX_CALIBRATION_SEMANTICS_VERSION.to_le_bytes());
        hash.update(&scope_evidence.0);
        hash.update(&(measurements.len() as u64).to_le_bytes());
        cx.checkpoint()?;
        Ok(Self { scene, measurements, prepared, scale: Some(Interval::new(0.0, f64::INFINITY)),
            conflict: None, processed: 0, identity: hash.finalize() })
    }

    /// Optional grid work may change only before any observation commits.
    /// The original grid limits remain binding; an unavailable row refuses fit.
    pub fn set_grid_enclosure_budget(&mut self, budget: GridEnclosureBudget) -> Result<(), FluxCalibrationError> {
        if self.processed != 0 || budget.max_nodes_per_segment == 0 {
            return Err(PhotometryError::InvalidOptions.into());
        }
        if let Prepared::Available { medium, .. } = &mut self.prepared { medium.set_grid_budget(budget); }
        Ok(())
    }

    /// Consume complete observations until the limit, completion or the first
    /// contradiction. No later constraint can repair an empty intersection.
    /// The numerical source value, constraint and evidence chain commit together.
    pub fn advance(&mut self, cx: &Cx<'_>, max_measurements: usize) -> Result<usize, FluxCalibrationError> {
        cx.checkpoint()?;
        let start = self.processed;
        let end = start.saturating_add(max_measurements).min(self.measurements.len());
        while self.processed < end && self.conflict.is_none() {
            let index = self.processed;
            let measurement = self.measurements[index];
            let nominal = evaluate(self.scene, cx, &measurement.probe, &self.prepared)?;
            let (model, model_identity) = model_range(nominal, index)?;
            let constraint = compatible_gain(model, measurement.measured, measurement.background);
            let candidate = self.scale.and_then(|old| constraint.and_then(|range| old.intersect(range)));
            let mut hash = DomainHasher::new("org.frankensim.photometry.flux-calibration.row.v1");
            hash.update(&self.identity.0);
            hash.update(&measurement.evidence.0);
            hash.update(&model_identity.0);
            for range in [model, measurement.measured, measurement.background] {
                hash_float(&mut hash, range.lo()); hash_float(&mut hash, range.hi());
            }
            match candidate {
                Some(range) => {
                    hash.update(&[1]); hash_float(&mut hash, range.lo()); hash_float(&mut hash, range.hi());
                }
                None => { hash.update(&[0]); }
            }
            let identity = hash.finalize();
            cx.checkpoint()?;
            self.scale = candidate;
            self.identity = identity;
            self.processed += 1;
            if candidate.is_none() { self.conflict = Some(index); }
        }
        Ok(self.processed - start)
    }

    /// Complete prefix size; no partially evaluated observation is counted.
    #[must_use]
    pub const fn processed_measurements(&self) -> usize { self.processed }
    /// All rows consumed or an early contradiction established.
    #[must_use]
    pub fn is_complete(&self) -> bool { self.conflict.is_some() || self.processed == self.measurements.len() }
    /// A prefix is not exported as a calibration. Complete unbounded and
    /// inconsistent outcomes are explicit rather than omitted or defaulted.
    #[must_use]
    pub fn result(&self) -> Option<FluxCalibrationResult> {
        if !self.is_complete() { return None; }
        let status = match (self.conflict, self.scale) {
            (Some(measurement), _) => FluxCalibrationStatus::Inconsistent { measurement },
            (_, Some(scale)) if finite(scale) => FluxCalibrationStatus::NotExcluded,
            _ => FluxCalibrationStatus::Unbounded,
        };
        Some(FluxCalibrationResult { scale: self.scale, status, processed: self.processed, identity: self.identity })
    }

    /// Check additional data using the SAME borrowed scene and prepared transport.
    /// This never narrows/refits the gain or alters accepted calibration history.
    /// Use independently held-out observations to test generalization; the API
    /// does not assert statistical independence or enforce a data-splitting policy.
    pub fn check_measurement(
        &self, cx: &Cx<'_>, measurement: &IlluminanceMeasurement,
    ) -> Result<FluxCalibrationCheck, FluxCalibrationError> {
        cx.checkpoint()?;
        let fit = self.result().ok_or(FluxCalibrationError::NoFiniteFit)?;
        let scale = fit.flux_scale_bound().ok_or(FluxCalibrationError::NoFiniteFit)?.scale();
        let nominal = evaluate(self.scene, cx, &measurement.probe, &self.prepared)?;
        let (model, _) = model_range(nominal, self.processed)?;
        let scaled = if model.hi() == 0.0 || scale.hi() == 0.0 { Interval::point(0.0) }
            else if scale == Interval::point(1.0) { model } else { nonnegative(model * scale) };
        let predicted = if measurement.background.hi() == 0.0 { scaled }
            else if scaled.hi() == 0.0 { measurement.background }
            else { nonnegative(scaled + measurement.background) };
        if !finite(predicted) { return Err(FluxCalibrationError::PredictionRange); }
        let compatibility = if predicted.intersect(measurement.measured).is_some() {
            MeasurementCompatibility::NotExcluded
        } else { MeasurementCompatibility::Contradicted };
        cx.checkpoint()?;
        Ok(FluxCalibrationCheck { nominal, predicted_measured_lux: predicted,
            measured_lux: measurement.measured, compatibility,
            calibration_identity: fit.identity, measurement_identity: measurement.evidence })
    }
}

fn model_range(value: EnclosedIlluminance, index: usize) -> Result<(Interval, ContentHash), FluxCalibrationError> {
    if let Some(reason) = value.unavailable {
        return Err(FluxCalibrationError::NumericalUnavailable { measurement: index, reason });
    }
    let certificate = value.estimate.numerical;
    if certificate.kind != fs_evidence::NumericalKind::Enclosure
        || value.numerical_model_identity.is_none()
        || !certificate.lo.is_finite() || !certificate.hi.is_finite()
        || certificate.lo < 0.0 || certificate.lo > certificate.hi {
        return Err(FluxCalibrationError::NumericalUnavailable {
            measurement: index, reason: EnclosureUnavailable::ArithmeticRange,
        });
    }
    let identity = value.numerical_model_identity.ok_or(FluxCalibrationError::NumericalUnavailable {
        measurement: index, reason: EnclosureUnavailable::ArithmeticRange,
    })?;
    Ok((Interval::new(certificate.lo, certificate.hi), identity))
}

/// Necessary gain range for M = g E + B, E >= 0, g >= 0. Treat zero support
/// explicitly, never divide by an interval containing zero or replace infinity
/// by a finite prior. Outward basic arithmetic encloses subtraction and division.
fn compatible_gain(model: Interval, measured: Interval, background: Interval) -> Option<Interval> {
    if measured.hi() < background.lo() { return None; }
    if model.hi() == 0.0 {
        return measured.intersect(background).map(|_| Interval::new(0.0, f64::INFINITY));
    }
    let signal_lower = if measured.lo() <= background.hi() { 0.0 }
        else { (Interval::point(measured.lo()) - Interval::point(background.hi())).lo().max(0.0) };
    // Equal binary64 inputs have an exact real difference zero. This preserves
    // a justified zero-gain constraint without promoting rounded underflow.
    let signal_upper = if measured.hi() == background.lo() { 0.0 }
        else { (Interval::point(measured.hi()) - Interval::point(background.lo())).hi() };
    let lower = if signal_lower == 0.0 { 0.0 }
        else { (Interval::point(signal_lower) / Interval::point(model.hi())).lo().max(0.0) };
    let upper = if model.lo() == 0.0 { f64::INFINITY }
        else if signal_upper == 0.0 { 0.0 }
        else { (Interval::point(signal_upper) / Interval::point(model.lo())).hi() };
    (lower <= upper).then(|| Interval::new(lower, upper))
}

#[cfg(test)]
mod tests;
