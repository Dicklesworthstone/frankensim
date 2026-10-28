//! Retain a static physical state AND its shared calibration errors across calls.
//!
//! A fresh `JacobianBatch::assimilate` marginalizes its nuisance against an
//! independent prior. Feeding that marginal into another such call would
//! incorrectly reset the shared error. This owner instead retains `[x, b]`,
//! including the learned cross-covariance, and updates with `[H, J]` and only
//! independent instrument noise. Batches may arrive in chunks; exact row replay
//! is a no-op and changed content under an already accepted row ID is refused.
//!
//! Linear physical forecasts preserve the retained shared-error correlations;
//! the nuisance itself remains static. State and nuisance start independent;
//! nuisance has zero prior mean. Callers own units,
//! calibration quality, and independent-noise declarations. State/row tracking
//! is in memory, not a durable replay database or an authenticated certificate.

use fs_blake3::{ContentHash, DomainHasher};

use super::RowId;
use super::jacobian::{JacobianBatch, JacobianRecord, SharedNoiseModel};
use crate::factor::FactorBelief;
use crate::factor::batch::{FactorBatchDiagnostics, assimilate_sequence};
use crate::{
    AssimError, Belief, Cx, MAX_DENSE_OBSERVATIONS, Observation, WorkPlan, WorkProgress,
    canonical_f64_bits, validate_state_dimension,
};

/// Maximum retained row identities. A stream also takes an explicit smaller cap.
pub const MAX_STREAM_ROWS: usize = MAX_DENSE_OBSERVATIONS;

#[derive(Debug, Clone, PartialEq)]
struct SeenRow {
    row: RowId,
    content: ContentHash,
}

/// Joint physical/nuisance posterior with transactional row admission.
#[derive(Debug, Clone, PartialEq)]
pub struct SharedNoiseStream {
    joint: FactorBelief,
    physical_dim: usize,
    model_id: String,
    row_capacity: usize,
    // Sorted by RowId for bounded binary lookup; never evict replay protection.
    seen: Vec<SeenRow>,
    history: ContentHash,
    // None retains the original static-stream semantics until a first forecast.
    time: Option<f64>,
}

/// Indices refer to the supplied batch, in its original declaration order.
#[derive(Debug, Clone, PartialEq)]
pub struct StreamUpdate {
    accepted_indices: Vec<usize>,
    replayed_indices: Vec<usize>,
    diagnostics: Option<FactorBatchDiagnostics>,
}

impl StreamUpdate {
    /// Newly assimilated rows, aligned with the factor diagnostics' receipts.
    #[must_use]
    pub fn accepted_indices(&self) -> &[usize] { &self.accepted_indices }

    /// Exact repeats ignored without a second state update.
    #[must_use]
    pub fn replayed_indices(&self) -> &[usize] { &self.replayed_indices }

    /// Joint-state update diagnostics for NEW rows only. These residuals use
    /// `[H,J]` and independent noise, not marginalized physical-state residuals.
    /// Replay-only calls have no numerical update and therefore no diagnostics.
    #[must_use]
    pub fn diagnostics(&self) -> Option<&FactorBatchDiagnostics> { self.diagnostics.as_ref() }
}

impl SharedNoiseStream {
    /// Start with a checked physical prior and independent zero-mean nuisance.
    /// The augmented dimension must fit the existing factor kernel's envelope.
    ///
    /// # Errors
    /// Refuses a zero/oversized row cap, unsupported joint shape, unresolved
    /// factorization, allocation failure, cancellation, or exhausted work.
    pub fn new(
        prior: &Belief,
        model: &SharedNoiseModel,
        row_capacity: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, AssimError> {
        if row_capacity == 0 || row_capacity > MAX_STREAM_ROWS {
            return Err(AssimError::ObservationCountLimit { count: row_capacity, max: MAX_STREAM_ROWS });
        }
        let n = prior.dim();
        let k = model.parameters().len();
        let d = n.checked_add(k).ok_or(AssimError::WorkPlanOverflow { phase: "shared stream shape" })?;
        validate_state_dimension(d)?;
        let du = d as u128;
        let mut progress = WorkProgress::new(cx, WorkPlan::checked(0, 0, 0, 0, 8 * du * du + 8 * du, 0, 16 * du * du + 1_024)?)?;
        progress.checkpoint("shared stream construction")?;
        let mut hash = DomainHasher::new("org.frankensim.fs-assimilate.shared-stream-initial.v1");
        field(&mut hash, model.content_id().as_bytes(), &mut progress)?;
        field(&mut hash, &(n as u64).to_le_bytes(), &mut progress)?;
        for value in prior.mean() { number(&mut hash, *value, &mut progress)?; }
        for row in prior.covariance() {
            for value in row { number(&mut hash, *value, &mut progress)?; }
        }
        let model_id = copy_string(model.content_id())?;
        // Initial independence is structural: preserve both admitted factors,
        // including exact nuisance null directions, without a larger PSD gate.
        let physical = FactorBelief::from_belief(prior, cx)?;
        let joint = physical.independent_product(&model.factor, cx)?;
        progress.checkpoint("shared stream publication")?;
        Ok(Self { joint, physical_dim: n, model_id, row_capacity, seen: Vec::new(), history: hash.finalize(), time: None })
    }

    /// Joint factor ordering: physical components, then named nuisance parameters.
    /// Do not drop its cross-covariance before processing the next shared batch.
    #[must_use]
    pub const fn joint(&self) -> &FactorBelief { &self.joint }

    /// Physical components occupy the first `physical_dim()` joint entries.
    #[must_use]
    pub const fn physical_dim(&self) -> usize { self.physical_dim }

    /// Immutable nuisance-model content identity admitted at construction.
    #[must_use]
    pub fn model_content_id(&self) -> &str { &self.model_id }

    /// Number of distinct successfully accepted dataset rows.
    #[must_use]
    pub fn accepted_rows(&self) -> usize { self.seen.len() }

    /// Root binding the initial prior/model and accepted ordered row contents.
    /// Chunk boundaries and exact replays do not alter this history root. It is
    /// a content identity, not authentication or a posterior accuracy claim.
    #[must_use]
    pub const fn history_id(&self) -> ContentHash { self.history }

    /// Current forecast frontier, or `None` before the first dynamic step.
    /// New observations passed to `assimilate` belong to this current state;
    /// this API does not backdate or smooth delayed measurements.
    #[must_use]
    pub const fn time(&self) -> Option<f64> { self.time }

    /// Reconstruct a checked physical marginal without discarding the retained
    /// joint factor. Only the physical block is reconstructed and PSD-checked;
    /// all joint factor columns still contribute to its covariance.
    ///
    /// # Errors
    /// Propagates dense reconstruction/PSD, allocation and budget refusals.
    pub fn physical_belief(&self, cx: &Cx<'_>) -> Result<Belief, AssimError> {
        self.joint.prefix_marginal(self.physical_dim, cx)
    }

    /// Assimilate a checked batch under the SAME retained nuisance model.
    /// All conflict/dimension/cap checks and factor updates finish before any
    /// retained state changes. Cancellation or numerical refusal leaves both
    /// the posterior and the replay set intact. No rows are silently evicted.
    ///
    /// Preparation uses a bounded work accountant; factor updates use the
    /// existing cumulative sequence accountant. Batch construction, preparation,
    /// update and optional dense publication retain their separate Cx scopes.
    /// The numerical loop is O(m (n+k)^2); JacobianBatch still performs its own
    /// bounded dense covariance admission before reaching this API.
    ///
    /// # Errors
    /// Refuses another model, changed content for an accepted row, mismatched
    /// physical dimensions, row-cap exhaustion, non-finite arithmetic,
    /// allocation failure, cancellation, or exhausted work/poll budgets.
    pub fn assimilate(&mut self, batch: &JacobianBatch, cx: &Cx<'_>) -> Result<StreamUpdate, AssimError> {
        if batch.model().content_id() != self.model_id {
            return Err(AssimError::SharedSourceDeclaration {
                id: batch.model().identity().to_owned(), reason: "batch differs from retained shared model",
            });
        }
        let m = batch.records().len();
        let d = self.joint.dim();
        // Constructors bound m<=256, d<=256 and retained rows<=4096.
        let (mu, du, su) = (m as u128, d as u128, self.seen.len() as u128);
        let work = mu * (8 * du + 16_384) + 1_024 * (su + mu + 1);
        let hashing = mu * (2_048 + 16 * du) + 256;
        let mut progress = WorkProgress::new(cx, WorkPlan::checked(0, 0, 0, 0, work, 0, hashing)?)?;
        progress.checkpoint("shared stream admission")?;
        let mut accepted_indices = reserve(m)?;
        let mut replayed_indices = reserve(m)?;
        let mut pending = reserve(m)?;
        let mut observations = reserve(m)?;
        let mut history = self.history;
        for (index, record) in batch.records().iter().enumerate() {
            if record.observation().operator().len() != self.physical_dim {
                return Err(AssimError::DimMismatch { state: self.physical_dim, operator: record.observation().operator().len() });
            }
            let content = row_content(record, &self.model_id, &mut progress)?;
            // At most 13 bounded leaf comparisons for the retained row cap.
            progress.comparison_bytes("shared stream row lookup", 13 * 256)?;
            match self.seen.binary_search_by(|seen| seen.row.cmp(record.row())) {
                Ok(position) => {
                    if self.seen[position].content != content {
                        return Err(AssimError::DuplicateDatasetRow { id: record.row().as_str().to_owned() });
                    }
                    replayed_indices.push(index);
                }
                Err(_) => {
                    let count = self.seen.len() + pending.len() + 1;
                    if count > self.row_capacity {
                        return Err(AssimError::ObservationCountLimit { count, max: self.row_capacity });
                    }
                    let mut operator = reserve(d)?;
                    operator.extend_from_slice(record.observation().operator());
                    operator.extend_from_slice(record.jacobian());
                    observations.push(Observation::new(
                        operator, record.observation().value(), record.observation().noise_var(),
                        format!("shared-stream-row:{content}"),
                    )?);
                    pending.push(SeenRow { row: RowId::try_new(copy_string(record.row().as_str())?)?, content });
                    accepted_indices.push(index);
                    let mut hash = DomainHasher::new("org.frankensim.fs-assimilate.shared-stream-history.v1");
                    hash.update(history.as_bytes());
                    hash.update(content.as_bytes());
                    progress.hash_bytes("shared stream history", 64)?;
                    history = hash.finalize();
                }
            }
            progress.scalar("shared stream rows", 4 * du + 1)?;
        }
        if pending.is_empty() {
            progress.checkpoint("shared stream replay")?;
            return Ok(StreamUpdate { accepted_indices, replayed_indices, diagnostics: None });
        }
        // Sort at most one already bounded batch, then merge linearly. Sorting
        // retained history on each call would turn tiny arrivals into long tiles.
        progress.comparison_bytes("shared stream pending order", mu * 32 * 256)?;
        pending.sort_unstable_by(|left, right| left.row.cmp(&right.row));
        let mut next_seen = reserve(self.seen.len() + pending.len())?;
        let mut old = self.seen.iter().peekable();
        for entry in pending {
            while old.peek().is_some_and(|seen| seen.row < entry.row) {
                let seen = old.next().expect("peeked retained row");
                next_seen.push(copy_seen(seen)?);
                progress.comparison_bytes("shared stream retained rows", 512)?;
            }
            next_seen.push(entry);
            progress.comparison_bytes("shared stream retained rows", 512)?;
        }
        for seen in old {
            next_seen.push(copy_seen(seen)?);
            progress.comparison_bytes("shared stream retained rows", 512)?;
        }
        let (joint, diagnostics) = assimilate_sequence(&self.joint, &observations, cx)?.into_parts();
        progress.checkpoint("shared stream commit")?;
        // No fallible work after publication begins.
        self.joint = joint;
        self.seen = next_seen;
        self.history = history;
        Ok(StreamUpdate { accepted_indices, replayed_indices, diagnostics: Some(diagnostics) })
    }
}

fn row_content(record: &JacobianRecord, model: &str, progress: &mut WorkProgress<'_, '_>) -> Result<ContentHash, AssimError> {
    let mut hash = DomainHasher::new("org.frankensim.fs-assimilate.shared-stream-row.v1");
    field(&mut hash, model.as_bytes(), progress)?;
    field(&mut hash, record.row().as_str().as_bytes(), progress)?;
    field(&mut hash, record.observation().instrument().as_bytes(), progress)?;
    field(&mut hash, &(record.observation().operator().len() as u64).to_le_bytes(), progress)?;
    for value in record.observation().operator() { number(&mut hash, *value, progress)?; }
    for value in record.jacobian() { number(&mut hash, *value, progress)?; }
    number(&mut hash, record.observation().value(), progress)?;
    number(&mut hash, record.observation().noise_var(), progress)?;
    Ok(hash.finalize())
}

fn reserve<T>(len: usize) -> Result<Vec<T>, AssimError> {
    let mut values = Vec::new();
    values.try_reserve_exact(len).map_err(|_| AssimError::AllocationRefused {
        stage: "shared stream storage", requested_bytes: len as u128 * core::mem::size_of::<T>() as u128,
    })?;
    Ok(values)
}

fn zeros(len: usize) -> Result<Vec<f64>, AssimError> {
    let mut values = reserve(len)?;
    values.resize(len, 0.0);
    Ok(values)
}

fn copy_string(value: &str) -> Result<String, AssimError> {
    let mut result = String::new();
    result.try_reserve_exact(value.len()).map_err(|_| AssimError::AllocationRefused {
        stage: "shared stream identity", requested_bytes: value.len() as u128,
    })?;
    result.push_str(value);
    Ok(result)
}

fn copy_seen(seen: &SeenRow) -> Result<SeenRow, AssimError> {
    Ok(SeenRow { row: RowId::try_new(copy_string(seen.row.as_str())?)?, content: seen.content })
}

fn field(hash: &mut DomainHasher, bytes: &[u8], progress: &mut WorkProgress<'_, '_>) -> Result<(), AssimError> {
    hash.update(&(bytes.len() as u64).to_le_bytes());
    hash.update(bytes);
    progress.hash_bytes("shared stream identity", bytes.len() as u128 + 8)
}

fn number(hash: &mut DomainHasher, value: f64, progress: &mut WorkProgress<'_, '_>) -> Result<(), AssimError> {
    hash.update(&canonical_f64_bits(value).to_le_bytes());
    progress.hash_bytes("shared stream identity", 8)
}

/// Outcome-free sensor planning under the retained joint posterior.
pub mod design;

mod forecast;
