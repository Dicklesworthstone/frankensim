//! Fixed-frame-lag ensemble smoothing by joint-state conditioning.
//!
//! Retain member-aligned states at explicit times and apply the existing scalar
//! or correlated square-root update to their concatenation. A later reading can
//! update earlier states through sample cross-time covariance, without a dense
//! covariance, new analysis algorithm, or repeated forward integration. This is
//! the direct augmented-ensemble smoother, not an iterative variational solver
//! or the lag-independent fast algorithms of Ravela--McLaughlin (2006).
//!
//! The window bounds a NUMBER of frames, not elapsed seconds. Observations must
//! name retained times exactly: no nearest-time substitution or interpolation.
//! Nonlinear smoothing is approximate and corrected histories need not satisfy
//! the nonlinear dynamics exactly. Only the newest corrected state is forecast.
//! Separate assimilation calls assume independent errors; cross-time correlated
//! readings must enter together as one declared block and must not be reused.

use super::{EnsembleForecast, EnsembleMoments, ForecastReport};
use super::super::{Ensemble, EnsembleAnalysis, EnsembleControl, EnsembleError,
    EnsembleObservation, finite, mean, poll, sum, zeros};
use super::super::correlated::{CorrelatedAnalysis, ObservationBlock};

/// A bounded joint ensemble. Clone forks a numerical checkpoint, including its
/// globally increasing observation cursor, but does not clone the work budget.
#[derive(Debug, Clone, PartialEq)]
pub struct EnsembleSmoother {
    joint: Ensemble,
    times: Vec<f64>,
    dimension: usize,
    max_frames: usize,
}

/// Owned forecast checkpoint. Completed member prefixes use the established
/// EnsembleForecast implementation; partial histories are never observable.
#[derive(Debug, Clone, PartialEq)]
pub struct SmootherForecast {
    forecast: EnsembleForecast,
    times: Vec<f64>,
    dimension: usize,
    max_frames: usize,
    last_observation: Option<u64>,
}

/// Borrowed physical states in observation-row order, including repeated times.
/// It cannot mutate or accidentally expose the concatenated ensemble as a field.
#[derive(Clone, Copy)]
pub struct ObservationStates<'a> {
    values: &'a [f64],
    times: &'a [f64],
    frames: &'a [usize],
    dimension: usize,
}
impl ObservationStates<'_> {
    pub fn len(&self) -> usize { self.frames.len() }
    pub fn is_empty(&self) -> bool { self.frames.is_empty() }
    pub fn time(&self, sample: usize) -> Option<f64> { self.times.get(sample).copied() }
    pub fn state(&self, sample: usize) -> Option<&[f64]> {
        self.frames.get(sample).map(|frame| {
            let start = frame * self.dimension;
            &self.values[start..start+self.dimension]
        })
    }
}

impl EnsembleSmoother {
    /// Consume the initial ensemble without copying its members. The capacity
    /// must be positive; capacity one is filtering with no retained history.
    pub fn new(initial: Ensemble, max_frames: usize) -> Result<Self, EnsembleError> {
        if max_frames == 0 { return Err(EnsembleError::Invalid("smoother needs at least one frame")); }
        initial.dimension.checked_mul(max_frames).and_then(|n| n.checked_mul(initial.count))
            .ok_or(EnsembleError::Invalid("smoother capacity extent overflow"))?;
        let mut times = zeros(1)?; times[0] = initial.time;
        Ok(Self { dimension: initial.dimension, joint: initial, times, max_frames })
    }
    pub fn times(&self) -> &[f64] { &self.times }
    pub fn time(&self) -> f64 { self.joint.time }
    pub fn dimension(&self) -> usize { self.dimension }
    pub fn member_count(&self) -> usize { self.joint.count }
    pub fn last_observation(&self) -> Option<u64> { self.joint.last_observation }
    pub fn member_at(&self, frame: usize, member: usize) -> Option<&[f64]> {
        if frame >= self.times.len() { return None; }
        self.joint.member(member).map(|row| &row[frame*self.dimension..(frame+1)*self.dimension])
    }
    fn frame(&self, time: f64) -> Result<usize, EnsembleError> {
        if !time.is_finite() { return Err(EnsembleError::Invalid("non-finite observation time")); }
        self.times.binary_search_by(|t| t.partial_cmp(&time).expect("finite retained time"))
            .map_err(|_| EnsembleError::Invalid("observation time is not retained in the smoother"))
    }

    /// Marginal sample moments for a retained frame, not confidence bounds.
    /// Uses the same centered/scaled reductions as Ensemble::moments; output
    /// is 2*dimension scalars. No member matrix or covariance is materialized.
    pub fn moments_at<C: FnMut() -> bool>(&self, frame: usize, control: &EnsembleControl,
        cancelled: &mut C) -> Result<EnsembleMoments, EnsembleError>
    {
        poll(cancelled)?;
        if frame >= self.times.len() { return Err(EnsembleError::Invalid("unknown smoother frame")); }
        let required = self.dimension.checked_mul(2).ok_or(EnsembleError::Invalid("moment extent overflow"))?;
        control.admit(0, required)?;
        let mut centers = zeros(self.dimension)?; let mut std = zeros(self.dimension)?;
        let start = frame*self.dimension;
        for i in 0..self.dimension {
            poll(cancelled)?;
            let values = self.joint.values.chunks_exact(self.joint.dimension).map(|x| x[start+i]);
            centers[i] = mean(values.clone(), self.joint.count)?;
            let mut scale = 0.0_f64;
            for value in values.clone() { scale = scale.max(finite(value-centers[i], "smoothing anomaly")?.abs()); }
            if scale != 0.0 {
                let variance = sum(values.map(|value| {
                    let d = (value-centers[i])/scale;
                    (d/(self.joint.count-1) as f64)*d
                }))?;
                std[i] = finite(scale*variance.sqrt(), "smoothing spread")?;
            }
        }
        poll(cancelled)?;
        Ok(EnsembleMoments { mean: centers, std })
    }

    /// Condition every retained frame on a scalar reading acquired at its
    /// declared time, even when that time precedes the latest forecast. IDs
    /// increase globally across forecasts, not in acquisition-time order.
    /// Callers must assign stable ingestion IDs and never renumber reused data.
    pub fn assimilate_scalar<F, C>(&mut self, observation: EnsembleObservation,
        predict: &mut F, control: &mut EnsembleControl, cancelled: &mut C)
        -> Result<EnsembleAnalysis, EnsembleError>
    where F: FnMut(usize, f64, &[f64]) -> Result<f64, String>, C: FnMut() -> bool,
    {
        poll(cancelled)?;
        let frame = self.frame(observation.time)?;
        let start = frame*self.dimension; let end = start+self.dimension;
        let aligned = EnsembleObservation { time: self.time(), ..observation };
        self.joint.assimilate_scalar(aligned, None,
            &mut |member, _, row| predict(member, observation.time, &row[start..end]), control, cancelled)
    }

    /// Condition a complete cross-time block atomically. The block's timestamp
    /// is the current analysis frontier; `sample_times` separately names each
    /// row's acquisition time. Its lower noise factor may couple any of those
    /// rows. Predictions run once per member on the original joint history.
    /// No localization or row-by-row nonlinear relinearization is implied.
    pub fn assimilate_correlated<F, C>(&mut self, block: &ObservationBlock,
        sample_times: &[f64], predict: &mut F, control: &mut EnsembleControl, cancelled: &mut C)
        -> Result<CorrelatedAnalysis, EnsembleError>
    where F: FnMut(usize, ObservationStates<'_>, &mut [f64], &mut dyn FnMut() -> bool) -> Result<(), String>,
        C: FnMut() -> bool,
    {
        poll(cancelled)?;
        if block.time() != self.time() || sample_times.len() != block.ids().len() {
            return Err(EnsembleError::Invalid("smoother block frontier or row count mismatch"));
        }
        let required = self.joint.correlated_workspace(block)?.checked_add(sample_times.len())
            .ok_or(EnsembleError::Invalid("smoother block workspace overflow"))?;
        control.admit(self.member_count(), required)?;
        let mut frames = Vec::new();
        frames.try_reserve_exact(sample_times.len()).map_err(|_| EnsembleError::Allocation)?;
        for time in sample_times { poll(cancelled)?; frames.push(self.frame(*time)?); }
        let dimension = self.dimension;
        self.joint.assimilate_correlated(block, &mut |member, _, values, out, check| {
            predict(member, ObservationStates { values, times: sample_times, frames: &frames, dimension }, out, check)
        }, control, cancelled)
    }

    /// Retain the newest max_frames-1 histories and append a forecast target.
    /// Evicted times are not recoverable. The source remains unchanged, so it
    /// may be inspected or explicitly forked while this checkpoint is running.
    /// Workspace includes two joint member matrices, a callback row and times;
    /// excludes the source smoother, model-owned memory and explicit clones.
    pub fn forecast<C: FnMut() -> bool>(&self, end: f64, control: &EnsembleControl,
        cancelled: &mut C) -> Result<SmootherForecast, EnsembleError>
    {
        poll(cancelled)?;
        if !end.is_finite() || end <= self.time() || !(end-self.time()).is_finite() {
            return Err(EnsembleError::Invalid("smoother forecast target must be finite and later"));
        }
        let keep = self.times.len().min(self.max_frames-1);
        let frames = keep+1;
        let width = frames.checked_mul(self.dimension).ok_or(EnsembleError::Invalid("smoother forecast extent"))?;
        let entries = width.checked_mul(self.member_count()).ok_or(EnsembleError::Invalid("smoother forecast extent"))?;
        let required = forecast_workspace(entries, width, frames)?;
        control.admit(0, required)?;
        let mut source = zeros(entries)?;
        let first = self.times.len()-keep;
        for (member, row) in source.chunks_mut(width).enumerate() {
            poll(cancelled)?;
            let original = self.joint.member(member).expect("bounded member");
            row[..keep*self.dimension].copy_from_slice(&original[first*self.dimension..]);
            // Last source slot holds the current field for the physical model,
            // not an invented history value at the yet-unreached target time.
            row[keep*self.dimension..].copy_from_slice(&original[original.len()-self.dimension..]);
        }
        let candidate = zeros(entries)?;
        let mut times = zeros(frames)?;
        times[..keep].copy_from_slice(&self.times[first..]); times[keep] = end;
        poll(cancelled)?;
        Ok(SmootherForecast { forecast: EnsembleForecast { start: self.time(), end,
            dimension: width, count: self.member_count(), source, candidate, completed: 0 },
            times, dimension: self.dimension, max_frames: self.max_frames,
            last_observation: self.last_observation() })
    }
}
fn forecast_workspace(entries: usize, width: usize, frames: usize) -> Result<usize, EnsembleError> {
    entries.checked_mul(2).and_then(|n| n.checked_add(width)).and_then(|n| n.checked_add(frames))
        .ok_or(EnsembleError::Invalid("smoother forecast workspace overflow"))
}
impl SmootherForecast {
    pub fn completed_members(&self) -> usize { self.forecast.completed_members() }
    pub fn is_complete(&self) -> bool { self.forecast.is_complete() }
    pub fn target_time(&self) -> f64 { self.forecast.target_time() }

    /// The producer sees ONLY one physical field and its start/end times, not
    /// concatenated history. Existing forecast code is reusable unchanged.
    /// Completed prefixes, cancellation latching and call charges are owned by
    /// EnsembleForecast. No observation can enter an incomplete forecast.
    pub fn advance<F, C>(&mut self, max_members: usize, propagate: &mut F,
        control: &mut EnsembleControl, cancelled: &mut C) -> Result<ForecastReport, EnsembleError>
    where F: FnMut(usize, f64, f64, &[f64], &mut [f64], &mut dyn FnMut() -> bool) -> Result<(), String>,
        C: FnMut() -> bool,
    {
        let width = self.forecast.dimension;
        let required = forecast_workspace(self.forecast.source.len(), width, self.times.len())?;
        control.admit(0, required)?;
        let start = width-self.dimension;
        self.forecast.advance(max_members, &mut |member, from, to, row, out, check| {
            out[..start].copy_from_slice(&row[..start]);
            propagate(member, from, to, &row[start..], &mut out[start..], check)
        }, control, cancelled)
    }
    /// Finish without allocation, preserving the global ingestion cursor.
    /// Replacing a newer, independently updated smoother with this result is
    /// an explicit fork rollback, not an automatic merge of observations.
    pub fn finish(self) -> Result<EnsembleSmoother, EnsembleError> {
        let mut joint = self.forecast.finish()?;
        joint.last_observation = self.last_observation;
        Ok(EnsembleSmoother { joint, times: self.times, dimension: self.dimension, max_frames: self.max_frames })
    }
}

#[cfg(test)]
#[path = "smoother_tests.rs"]
mod tests;
