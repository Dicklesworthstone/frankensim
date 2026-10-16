//! Attach real row-addressed readings without changing the planned response.
use super::{ProbeBatchPlan, SharedNoiseStream, copy_string, declaration, reserve};
use crate::groups::{RowId, jacobian::{JacobianBatch, SharedNoiseModel}, stream::StreamUpdate};
use crate::{AssimError, Cx, WorkPlan, WorkProgress, canonicalize_zero, validate_leaf_identity};

/// A real reading identified by acquired dataset row, not by vector position.
/// Response, shared sensitivity, variance and cost come from the sealed plan.
#[derive(Debug, Clone, PartialEq)]
pub struct ProbeReading {
    row: RowId,
    value: f64,
    instrument: String,
}

impl ProbeReading {
    /// Attach an acquired scalar value and explicit instrument identity.
    /// This does not certify the instrument or the declared noise model.
    ///
    /// # Errors
    /// Refuses non-finite readings or an invalid instrument identity.
    pub fn new(row: RowId, value: f64, instrument: impl Into<String>) -> Result<Self, AssimError> {
        let instrument = instrument.into();
        validate_leaf_identity("probe instrument", &instrument)?;
        if !value.is_finite() { return Err(AssimError::NonFiniteObservationValue); }
        Ok(Self { row, value: canonicalize_zero(value), instrument })
    }
    /// The exact row selected by the plan; acquisition order may differ.
    #[must_use]
    pub const fn row(&self) -> &RowId { &self.row }
    /// The actual measured value, never used in sensor selection.
    #[must_use]
    pub const fn value(&self) -> f64 { self.value }
    /// Explicit instrument identity for the accepted-row history.
    #[must_use]
    pub fn instrument(&self) -> &str { &self.instrument }
}

impl ProbeBatchPlan {
    /// Check the live history/model/shape and available row capacity BEFORE
    /// performing an external acquisition. This check reserves no row or cost.
    /// `assimilate_readings` checks again before publishing any new information.
    ///
    /// # Errors
    /// Refuses a changed posterior, another model/physical shape, or a stream
    /// whose remaining row capacity cannot retain the whole selected set.
    pub fn check_current(&self, stream: &SharedNoiseStream) -> Result<(), AssimError> {
        if stream.model_id != self.model_id || stream.physical_dim != self.physical_dim {
            return Err(declaration("acquisition plan differs from the live stream model or shape"));
        }
        if stream.history != self.history {
            return Err(declaration("acquisition plan is stale; re-plan against the current posterior"));
        }
        let count = stream.seen.len() + self.steps.len();
        if count > stream.row_capacity {
            return Err(AssimError::ObservationCountLimit { count, max: stream.row_capacity });
        }
        Ok(())
    }

    /// Assimilate real readings matching the complete planned set, atomically.
    /// Input may arrive in any order: row IDs map values to the exact selected
    /// response and updates run in selection order. Missing, extra, duplicate,
    /// stale or foreign-model data are never silently attached to another row.
    /// The stream's transactional factor update remains the sole numerical
    /// owner; predicted covariance is NOT substituted for a measured update.
    /// Returned update indices refer to `plan.steps()`, not the input order.
    ///
    /// Failure does not consume stream rows or change its posterior/history.
    /// Borrowed readings remain available for retry after cancellation or an
    /// increased compute budget. An already executed nonempty plan is stale;
    /// stream-level exact replay remains available through its ordinary API.
    /// Empty plans accept only empty readings and return a checked no-op.
    ///
    /// Preparation, checked covariance construction and retained assimilation
    /// keep their existing separate compute scopes. Acquisition cost is a
    /// planning limit, not a payment/refund system. External sensing is outside
    /// this function and cannot be undone by numerical rollback.
    ///
    /// # Errors
    /// Propagates admission, row mapping, covariance/factor, budget, allocation
    /// and cancellation failures. No fallible work follows the stream commit.
    pub fn assimilate_readings(
        &self,
        stream: &mut SharedNoiseStream,
        model: &SharedNoiseModel,
        readings: &[ProbeReading],
        cx: &Cx<'_>,
    ) -> Result<StreamUpdate, AssimError> {
        self.check_current(stream)?;
        if model.content_id() != self.model_id {
            return Err(declaration("readings must use the plan's original nuisance model"));
        }
        let m = self.steps.len();
        if readings.len() != m {
            return Err(declaration("readings must cover every selected row exactly once"));
        }
        let (mu, du) = (m as u128, stream.joint.dim() as u128);
        let work = 512 * mu * mu + (mu + 1) * (8 * du + 2_048);
        let mut progress = WorkProgress::new(cx, WorkPlan::checked(0, 0, 0, 0, work, 0, 0)?)?;
        progress.checkpoint("planned acquisition admission")?;
        let mut order = reserve(m)?;
        order.resize(m, None);
        for (index, reading) in readings.iter().enumerate() {
            let mut found = false;
            for (slot, step) in self.steps.iter().enumerate() {
                progress.comparison_bytes("planned acquisition row lookup", 256)?;
                if reading.row == step.candidate.probe.row {
                    if order[slot].is_some() {
                        return Err(AssimError::DuplicateDatasetRow { id: copy_string(reading.row.as_str())? });
                    }
                    order[slot] = Some(index);
                    found = true;
                    break;
                }
            }
            if !found { return Err(declaration("acquired row was not selected by this plan")); }
        }
        let mut records = reserve(m)?;
        for (step, index) in self.steps.iter().zip(order) {
            let reading = &readings[index.ok_or_else(|| declaration("missing selected reading"))?];
            records.push(step.candidate.probe.record(model, reading.value, copy_string(&reading.instrument)?)?);
            progress.scalar("planned acquisition materialization", 4 * du + 512)?;
        }
        progress.checkpoint("planned acquisition prepared")?;
        if records.is_empty() {
            return Ok(StreamUpdate { accepted_indices: Vec::new(), replayed_indices: Vec::new(), diagnostics: None });
        }
        let batch = JacobianBatch::new(model.clone(), records, cx)?;
        progress.checkpoint("planned acquisition ready")?;
        // This call has its own final pre-commit cancellation check. In
        // particular, do not perform fallible diagnostics after it commits.
        stream.assimilate(&batch, cx)
    }
}
