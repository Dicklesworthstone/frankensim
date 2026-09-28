//! Resumable ensemble forecasts and explicit anomaly inflation.
//!
//! Every forecast member starts from the same immutable snapshot time. Finished
//! members are retained across pauses; no mixed-time ensemble is observable.
//! `finish` returns a complete new ensemble instead of mutating the source.
//! Model callbacks must be pure functions of member ID, source, times and fixed
//! model data. Stochastic forecasts must key their own streams by these inputs,
//! not by callback attempt counts, so retries do not change their meaning.

use super::{Ensemble, EnsembleControl, EnsembleError, finite, mean, poll, sum, zeros};

#[path = "smoother.rs"]
pub mod smoother;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ForecastStatus { Complete, MemberLimit, Cancelled }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ForecastReport {
    pub status: ForecastStatus,
    pub completed_this_call: usize,
    pub completed_members: usize,
}

/// Owned, cloneable forecast checkpoint. Only completed member outputs are
/// retained; no model or work budget is cloned. It cannot be used for analysis
/// before every member is at the target time. Clone copies two ensembles.
#[derive(Debug, Clone, PartialEq)]
pub struct EnsembleForecast {
    start: f64,
    end: f64,
    dimension: usize,
    count: usize,
    source: Vec<f64>,
    candidate: Vec<f64>,
    completed: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct EnsembleMoments {
    pub mean: Vec<f64>,
    /// Marginal sample standard deviations, not confidence bounds.
    pub std: Vec<f64>,
}

impl Ensemble {
    /// Allocate a forecast checkpoint. Its workspace ceiling covers the two
    /// owned ensembles plus a single-member transactional callback buffer;
    /// excludes this source object, model memory and explicit checkpoint clones.
    pub fn forecast<C: FnMut() -> bool>(
        &self, end: f64, control: &EnsembleControl, cancelled: &mut C,
    ) -> Result<EnsembleForecast, EnsembleError> {
        poll(cancelled)?;
        if !end.is_finite() || end <= self.time || !(end-self.time).is_finite() {
            return Err(EnsembleError::Invalid("forecast target must be finite and strictly later"));
        }
        let required = self.values.len().checked_mul(2).and_then(|n| n.checked_add(self.dimension))
            .ok_or(EnsembleError::Invalid("forecast workspace overflow"))?;
        control.admit(0, required)?;
        let mut source = zeros(self.values.len())?;
        for (dst, src) in source.chunks_mut(256).zip(self.values.chunks(256)) {
            poll(cancelled)?; dst.copy_from_slice(src);
        }
        let candidate = zeros(self.values.len())?;
        poll(cancelled)?;
        Ok(EnsembleForecast { start: self.time, end, dimension: self.dimension,
            count: self.count, source, candidate, completed: 0 })
    }

    /// Multiply anomaly covariance by a supplied positive factor. This is an
    /// explicit modeling choice, never automatically reapplied on retry or per
    /// observation. It preserves the observation cursor and cannot justify
    /// consuming a measurement again. Factor one is a bitwise no-op.
    pub fn rescale_spread<C: FnMut() -> bool>(
        &mut self, variance_factor: f64, control: &EnsembleControl, cancelled: &mut C,
    ) -> Result<(), EnsembleError> {
        poll(cancelled)?;
        if !variance_factor.is_finite() || variance_factor <= 0.0 {
            return Err(EnsembleError::Invalid("variance multiplier must be finite and positive"));
        }
        if variance_factor == 1.0 { return Ok(()); }
        control.admit(0, self.values.len())?;
        let scale = variance_factor.sqrt();
        let mut candidate = zeros(self.values.len())?;
        for i in 0..self.dimension {
            poll(cancelled)?;
            let center = mean(self.values.chunks_exact(self.dimension).map(|x| x[i]), self.count)?;
            for j in 0..self.count {
                if j % 256 == 0 { poll(cancelled)?; }
                let offset = j*self.dimension+i;
                let deviation = finite(self.values[offset]-center, "inflation anomaly")?;
                candidate[offset] = finite(scale.mul_add(deviation, center), "inflated member")?;
            }
        }
        poll(cancelled)?;
        self.values = candidate;
        Ok(())
    }

    /// Return fieldwise sample moments using O(dimension) output and no dense
    /// covariance. Output scalars are included in the workspace ceiling.
    pub fn moments<C: FnMut() -> bool>(&self, control: &EnsembleControl, cancelled: &mut C)
        -> Result<EnsembleMoments, EnsembleError>
    {
        poll(cancelled)?;
        let required = self.dimension.checked_mul(2).ok_or(EnsembleError::Invalid("moment extent overflow"))?;
        control.admit(0, required)?;
        let mut means = zeros(self.dimension)?; let mut std = zeros(self.dimension)?;
        for i in 0..self.dimension {
            poll(cancelled)?;
            means[i] = mean(self.values.chunks_exact(self.dimension).map(|x| x[i]), self.count)?;
            let mut scale = 0.0_f64;
            for x in self.values.chunks_exact(self.dimension) {
                scale = scale.max(finite(x[i]-means[i], "moment anomaly")?.abs());
            }
            if scale != 0.0 {
                let variance = sum(self.values.chunks_exact(self.dimension).map(|x| {
                    let d = (x[i]-means[i])/scale;
                    (d/(self.count-1) as f64)*d
                }))?;
                std[i] = finite(scale*variance.sqrt(), "sample standard deviation")?;
            }
        }
        poll(cancelled)?;
        Ok(EnsembleMoments { mean: means, std })
    }
}

impl EnsembleForecast {
    pub fn completed_members(&self) -> usize { self.completed }
    pub fn target_time(&self) -> f64 { self.end }
    pub fn is_complete(&self) -> bool { self.completed == self.count }

    /// Complete at most `max_members` more forecasts. The callback overwrites
    /// every output component and receives a latched cancellation check for
    /// long kernels. A failed member is retryable; completed members are not
    /// recalculated, and attempted calls remain charged in `control`.
    pub fn advance<F, C>(
        &mut self, max_members: usize, propagate: &mut F,
        control: &mut EnsembleControl, cancelled: &mut C,
    ) -> Result<ForecastReport, EnsembleError>
    where
        F: FnMut(usize, f64, f64, &[f64], &mut [f64], &mut dyn FnMut() -> bool) -> Result<(), String>,
        C: FnMut() -> bool,
    {
        let initial = self.completed;
        let report = |status, completed| ForecastReport {
            status, completed_this_call: completed-initial, completed_members: completed,
        };
        if cancelled() { return Ok(report(ForecastStatus::Cancelled, self.completed)); }
        if self.is_complete() { return Ok(report(ForecastStatus::Complete, self.completed)); }
        let count = max_members.min(self.count-self.completed);
        let required = self.source.len().checked_mul(2).and_then(|n| n.checked_add(self.dimension))
            .ok_or(EnsembleError::Invalid("forecast workspace overflow"))?;
        control.admit(count, required)?;
        if count == 0 { return Ok(report(ForecastStatus::MemberLimit, self.completed)); }
        let mut out = zeros(self.dimension)?;
        for _ in 0..count {
            if cancelled() { return Ok(report(ForecastStatus::Cancelled, self.completed)); }
            let member = self.completed; let start = member*self.dimension;
            out.fill(f64::NAN);
            control.charge()?;
            let mut stopped = false;
            let mut check = || { stopped |= cancelled(); stopped };
            let result = propagate(member, self.start, self.end,
                &self.source[start..start+self.dimension], &mut out, &mut check);
            if check() { return Ok(report(ForecastStatus::Cancelled, self.completed)); }
            result.map_err(|message| EnsembleError::Model { member, message })?;
            for chunk in out.chunks(256) {
                if cancelled() { return Ok(report(ForecastStatus::Cancelled, self.completed)); }
                if chunk.iter().any(|v| !v.is_finite()) { return Err(EnsembleError::NonFinite("forecast output")); }
            }
            if cancelled() { return Ok(report(ForecastStatus::Cancelled, self.completed)); }
            self.candidate[start..start+self.dimension].copy_from_slice(&out);
            self.completed += 1;
        }
        Ok(report(if self.is_complete() { ForecastStatus::Complete } else { ForecastStatus::MemberLimit }, self.completed))
    }

    /// Consume a COMPLETE checkpoint into a new ensemble. This does no model
    /// work or allocation. The new timestamp permits fresh observation IDs.
    /// Keep a checkpoint instead of calling `finish` on an incomplete forecast.
    pub fn finish(self) -> Result<Ensemble, EnsembleError> {
        if !self.is_complete() { return Err(EnsembleError::Invalid("forecast is incomplete")); }
        Ok(Ensemble { time: self.end, dimension: self.dimension, count: self.count,
            values: self.candidate, last_observation: None })
    }
}

#[cfg(test)]
#[path = "forecast_tests.rs"]
mod tests;
