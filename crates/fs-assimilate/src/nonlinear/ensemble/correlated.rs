//! Same-time correlated observations through the existing square-root update.
//!
//! Supply a lower noise factor L, defining R = L L^T. Whitening is the
//! triangular solve L z = h(x)-y, NOT multiplication by R^-1 and not independent
//! marginal weighting. Append these predicted residuals to a temporary joint
//! ensemble and condition on z=0 with unit noise using the scalar analysis
//! kernel. Remaining predicted residuals are updated along with physical state.
//! No state covariance, inverse, perturbed observations or new eigensolver is
//! used. For linear h this gives the sample-covariance Kalman update; nonlinear
//! h uses one frozen joint ensemble regression, not exact Bayesian inference.
//!
//! Whitening mixes sensor supports. Localization is deliberately not accepted:
//! applying original sensor tapers to mixed observations changes the model.
//! Cross-time correlated errors require an augmented forecast model or smoother,
//! not repeated calls to independent same-time blocks. References: the Kalman
//! conditioning identities and Whitaker--Hamill (2002), MWR 130, 1913--1924.

use super::{Ensemble, EnsembleControl, EnsembleError, EnsembleObservation, finite, poll, zeros};

/// Immutable simultaneous readings and their explicitly supplied noise square
/// root. Rows, columns and prediction outputs use exactly the given ID order.
/// The upper triangle must be zero; positive diagonal gives a nonsingular
/// covariance by construction. This checks a declaration, not noise calibration.
#[derive(Debug, Clone, PartialEq)]
pub struct ObservationBlock {
    time: f64,
    ids: Vec<u64>,
    values: Vec<f64>,
    lower: Vec<f64>,
}
impl ObservationBlock {
    /// Readings must have strictly increasing IDs; no automatic sorting can
    /// silently detach a covariance row from its sensor. `max_components`
    /// bounds owned f64/u64 elements (m*m+2*m), excluding allocator metadata.
    pub fn new(
        time: f64, ids: &[u64], values: &[f64], noise_lower: &[f64],
        max_observations: usize, max_components: usize,
    ) -> Result<Self, EnsembleError> {
        let m = ids.len();
        let matrix = m.checked_mul(m).ok_or(EnsembleError::Invalid("noise factor extent overflow"))?;
        let required = m.checked_mul(2).and_then(|v| matrix.checked_add(v))
            .ok_or(EnsembleError::Invalid("observation block extent overflow"))?;
        if required > max_components {
            return Err(EnsembleError::WorkspaceLimit { required, limit: max_components });
        }
        if m == 0 || m > max_observations || values.len() != m || noise_lower.len() != matrix
            || !time.is_finite() || values.iter().any(|x| !x.is_finite())
        { return Err(EnsembleError::Invalid("finite, bounded, shape-matched observation block required")); }
        if ids.windows(2).any(|pair| pair[0] >= pair[1]) { return Err(EnsembleError::ObservationOrder); }
        for i in 0..m {
            for j in 0..m {
                let value = noise_lower[i*m+j];
                if !value.is_finite() || (j > i && value != 0.0) || (j == i && value <= 0.0) {
                    return Err(EnsembleError::Invalid("noise factor must be finite, lower triangular, positive diagonal"));
                }
            }
        }
        let mut owned_ids = Vec::new();
        owned_ids.try_reserve_exact(m).map_err(|_| EnsembleError::Allocation)?;
        owned_ids.extend_from_slice(ids);
        let mut owned_values = zeros(m)?; owned_values.copy_from_slice(values);
        let mut lower = zeros(matrix)?; lower.copy_from_slice(noise_lower);
        Ok(Self { time, ids: owned_ids, values: owned_values, lower })
    }
    /// Independent errors plus one common reference offset at THIS timestamp:
    /// R_ij = sigma_i^2 * delta_ij + common_sigma^2. The offset is newly drawn
    /// independently for each block, not a bias persistent across forecasts.
    /// All inputs use the same signal unit. No covariance squares are formed.
    ///
    /// Scalar Gaussian conditioning builds its lower factor directly: the
    /// remaining common-source standard deviation is updated after each row.
    /// `max_components` covers construction's two m-by-m buffers and 2*m
    /// retained values/IDs, excluding caller inputs and allocator metadata.
    #[allow(clippy::too_many_arguments)]
    pub fn shared_reference<C: FnMut() -> bool>(
        time: f64, ids: &[u64], values: &[f64], independent_sigma: &[f64], common_sigma: f64,
        max_observations: usize, max_components: usize, cancelled: &mut C,
    ) -> Result<Self, EnsembleError> {
        poll(cancelled)?;
        let m = ids.len();
        let matrix = m.checked_mul(m).ok_or(EnsembleError::Invalid("shared-reference extent overflow"))?;
        let required = matrix.checked_add(m).and_then(|v| v.checked_mul(2))
            .ok_or(EnsembleError::Invalid("shared-reference extent overflow"))?;
        if required > max_components {
            return Err(EnsembleError::WorkspaceLimit { required, limit: max_components });
        }
        if m == 0 || m > max_observations || values.len() != m || independent_sigma.len() != m
            || !time.is_finite() || !common_sigma.is_finite() || common_sigma < 0.0
            || values.iter().any(|v| !v.is_finite())
            || independent_sigma.iter().any(|s| !s.is_finite() || *s <= 0.0)
        { return Err(EnsembleError::Invalid("invalid shared-reference readings or noise scales")); }
        if ids.windows(2).any(|p| p[0] >= p[1]) { return Err(EnsembleError::ObservationOrder); }
        let mut lower = zeros(matrix)?;
        let mut remaining = common_sigma;
        for j in 0..m {
            poll(cancelled)?;
            let independent = independent_sigma[j];
            let diagonal = finite(independent.hypot(remaining), "shared-reference diagonal")?;
            lower[j*m+j] = diagonal;
            let common = finite(remaining*(remaining/diagonal), "shared-reference loading")?;
            for i in j+1..m {
                if i % 256 == 0 { poll(cancelled)?; }
                lower[i*m+j] = common;
            }
            // min * (max/hypot) avoids spurious zero from an underflowed ratio
            // when the two positive scales are many orders of magnitude apart.
            remaining = remaining.min(independent)*(remaining.max(independent)/diagonal);
        }
        poll(cancelled)?;
        let block = Self::new(time, ids, values, &lower, max_observations, max_components)?;
        poll(cancelled)?;
        Ok(block)
    }

    pub fn time(&self) -> f64 { self.time }
    pub fn ids(&self) -> &[u64] { &self.ids }
    pub fn values(&self) -> &[f64] { &self.values }
    pub fn noise_lower(&self) -> &[f64] { &self.lower }

    fn whiten<C: FnMut() -> bool>(&self, values: &mut [f64], cancelled: &mut C)
        -> Result<(), EnsembleError>
    {
        let m = self.ids.len();
        for i in 0..m {
            poll(cancelled)?;
            let diagonal = self.lower[i*m+i];
            // Form dimensionless rows before multiplying. A changed choice of
            // physical units then cannot by itself overflow covariance squares.
            let mut value = finite(values[i] / diagonal, "normalized observation residual")?;
            for j in 0..i {
                if j % 256 == 0 { poll(cancelled)?; }
                let coefficient = finite(self.lower[i*m+j] / diagonal, "normalized noise factor")?;
                value = finite((-coefficient).mul_add(values[j], value), "observation whitening")?;
            }
            values[i] = value;
        }
        Ok(())
    }
}

/// Conditional diagnostics in the declared whitening order, NOT independent
/// physical-sensor innovations. Row i can mix sensors 0..i. These are numerical
/// diagnostics, not a chi-square acceptance test or a coverage certificate.
#[derive(Debug, Clone, PartialEq)]
pub struct CorrelatedAnalysis {
    pub first_observation: u64,
    pub last_observation: u64,
    pub standardized_innovations: Vec<f64>,
    pub square_root_factors: Vec<f64>,
    /// One complete vector prediction per member, not one call per sensor.
    pub model_calls: usize,
}

impl Ensemble {
    /// Conservative simultaneously live scalar envelope, excluding the input
    /// ensemble/block and model-owned work: 2*N*(n+m)+N+2*m. The observation
    /// count also bounds whitening O(N*m*m) and updates O(m*N*(n+m)).
    pub fn correlated_workspace(&self, block: &ObservationBlock) -> Result<usize, EnsembleError> {
        let m = block.ids.len();
        self.dimension.checked_add(m).and_then(|d| d.checked_mul(self.count))
            .and_then(|v| v.checked_mul(2)).and_then(|v| v.checked_add(self.count))
            .and_then(|v| m.checked_mul(2).and_then(|w| v.checked_add(w)))
            .ok_or(EnsembleError::Invalid("correlated workspace extent overflow"))
    }

    /// Assimilate the entire block or publish nothing. `predict` overwrites all
    /// m predictions in the block's ID order. It is called once on each original
    /// forecast member; later rows condition the SAME joint sample regression.
    /// This differs intentionally from scalar nonlinear relinearization.
    ///
    /// Cancellation is latched inside each model call. All attempts stay charged
    /// on error, including a bad member output. The complete block is retryable;
    /// partially applied rows never escape and the observation cursor advances
    /// only after all rows and the final publication checkpoint succeed.
    pub fn assimilate_correlated<F, C>(
        &mut self, block: &ObservationBlock, predict: &mut F,
        control: &mut EnsembleControl, cancelled: &mut C,
    ) -> Result<CorrelatedAnalysis, EnsembleError>
    where
        F: FnMut(usize, f64, &[f64], &mut [f64], &mut dyn FnMut() -> bool) -> Result<(), String>,
        C: FnMut() -> bool,
    {
        poll(cancelled)?;
        if block.time != self.time { return Err(EnsembleError::Invalid("block and ensemble times differ")); }
        if self.last_observation.is_some_and(|last| block.ids[0] <= last) {
            return Err(EnsembleError::ObservationOrder);
        }
        control.admit(self.count, self.correlated_workspace(block)?)?;
        let m = block.ids.len();
        let width = self.dimension + m; // checked in the workspace admission
        let mut joint = Ensemble { time: self.time, dimension: width, count: self.count,
            values: zeros(width*self.count)?, last_observation: self.last_observation };
        for member in 0..self.count {
            poll(cancelled)?;
            let state = self.member(member).expect("bounded member index");
            let row = &mut joint.values[member*width..(member+1)*width];
            row[..self.dimension].copy_from_slice(state);
            let predictions = &mut row[self.dimension..];
            predictions.fill(f64::NAN);
            control.charge()?;
            let mut stopped = false;
            let mut check = || { stopped |= cancelled(); stopped };
            let result = predict(member, self.time, state, predictions, &mut check);
            if check() { return Err(EnsembleError::Cancelled); }
            result.map_err(|message| EnsembleError::Model { member, message })?;
            for (i, value) in predictions.iter_mut().enumerate() {
                if i % 256 == 0 { poll(cancelled)?; }
                finite(*value, "block prediction")?;
                *value = finite(*value-block.values[i], "block prediction residual")?;
            }
            block.whiten(predictions, cancelled)?;
        }
        let mut report = CorrelatedAnalysis {
            first_observation: block.ids[0], last_observation: block.ids[m-1],
            standardized_innovations: zeros(m)?, square_root_factors: zeros(m)?, model_calls: self.count,
        };
        for (i, &id) in block.ids.iter().enumerate() {
            poll(cancelled)?;
            let mut predictions = zeros(self.count)?;
            for (member, value) in predictions.iter_mut().enumerate() {
                if member % 256 == 0 { poll(cancelled)?; }
                *value = joint.values[member*width+self.dimension+i];
            }
            // Reuse the identical scalar kernel. Reading stored predictions is
            // not a model evaluation: do not mint or reset a second call ledger.
            let analysis = joint.analyze_predictions(EnsembleObservation {
                id, time: self.time, value: 0.0, sigma: 1.0,
            }, None, predictions, cancelled)?;
            report.standardized_innovations[i] = analysis.standardized_innovation;
            report.square_root_factors[i] = analysis.square_root_factor;
        }
        let mut candidate = zeros(self.values.len())?;
        for (dst, src) in candidate.chunks_mut(self.dimension).zip(joint.values.chunks(width)) {
            poll(cancelled)?;
            dst.copy_from_slice(&src[..self.dimension]);
        }
        poll(cancelled)?;
        self.values = candidate;
        self.last_observation = Some(report.last_observation);
        Ok(report)
    }
}

#[cfg(test)]
#[path = "correlated_tests.rs"]
mod tests;
