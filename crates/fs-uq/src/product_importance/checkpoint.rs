//! Proposal-bound recovery of paid physical responses and likelihood weights.
//! A trusted observation stream is retained, not replayed through the model.
use fs_blake3::{ContentHash, hash_domain};
use super::{ImportanceConfig, ImportanceExecution, ImportanceObservation};
use crate::product_checkpoint::execution_identity;
use crate::product_plan::admit_plan;
use crate::{UqCheckpointError, UqPlan, UqStatus};

const MAGIC: &[u8; 8] = b"FSIMP001";
const HEADER: usize = 8 + 32 + 1 + 8;
const FIXED: usize = HEADER + 32;
const RECORD: usize = 16;
const DOMAIN: &str = "org.frankensim.uq.importance.checkpoint.v1";

fn identity(plan: &UqPlan, proposal: &ImportanceConfig, model: ContentHash) -> ContentHash {
    let mut bytes = Vec::with_capacity(48 + 8 * proposal.germ_shift.len());
    bytes.extend_from_slice(execution_identity(plan, model).as_bytes());
    bytes.extend_from_slice(&proposal.defensive_fraction.to_bits().to_le_bytes());
    bytes.extend_from_slice(&(proposal.germ_shift.len() as u64).to_le_bytes());
    for shift in &proposal.germ_shift { bytes.extend_from_slice(&shift.to_bits().to_le_bytes()); }
    // Binds mixture-component addressing, shifted-germ order, log-density
    // arithmetic and the unit-weight target shortcut in addition to the plan.
    hash_domain("org.frankensim.uq.importance.execution.v1", &bytes)
}

impl ImportanceExecution {
    /// Retain every accepted QoI and epsilon-scaled likelihood weight as exact
    /// f64 bits. The strict failure flag and unscaled weight are derived on read.
    /// A complete checkpoint remains terminal; a refused run cannot be saved.
    ///
    /// The model identity must cover all fixed physical inputs, evaluator and
    /// runtime semantics. The complete plan and proposal are also bound. This
    /// checksum detects corruption, NOT authentication or evidence that a model
    /// produced these responses/weights. Restore only trusted checkpoint data.
    pub fn checkpoint(&self, model: ContentHash) -> Result<Vec<u8>, UqCheckpointError> {
        let status = match self.status {
            UqStatus::BudgetTruncated => 0,
            UqStatus::Cancelled => 1,
            UqStatus::Complete => 2,
            UqStatus::Refused => return Err(UqCheckpointError::NotResumable),
        };
        if self.failure.is_some() || self.attempted != self.values.len() {
            return Err(UqCheckpointError::NotResumable);
        }
        let mut bytes = Vec::new();
        bytes.try_reserve_exact(FIXED + RECORD * self.values.len())
            .map_err(|_| UqCheckpointError::InvalidEncoding("importance checkpoint allocation refused"))?;
        bytes.extend_from_slice(MAGIC);
        bytes.extend_from_slice(identity(&self.plan, &self.config, model).as_bytes());
        bytes.push(status);
        bytes.extend_from_slice(&(self.values.len() as u64).to_le_bytes());
        for value in &self.values {
            bytes.extend_from_slice(&value.qoi.to_bits().to_le_bytes());
            bytes.extend_from_slice(&value.scaled_weight.to_bits().to_le_bytes());
        }
        let checksum = hash_domain(DOMAIN, &bytes);
        bytes.extend_from_slice(checksum.as_bytes());
        Ok(bytes)
    }

    /// Restore in O(n) decoding work without any physical callback or sampler
    /// replay. Re-admit the Gaussian law and proposal; bound framing/counts and
    /// verify the checksum/identity before reserving the observation allocation.
    /// Statistics and failure flags are reconstructed, not trusted from a report.
    ///
    /// Weights are retained observations: range checks do not prove that a
    /// trusted producer computed the right likelihood at the claimed ordinal.
    /// No changed proposal, ceiling, law, model, seed or lifetime cap is allowed.
    /// This restores one fixed producer, not an adaptive pilot owner's history.
    pub fn restore(
        plan: &UqPlan, proposal: ImportanceConfig, model: ContentHash, bytes: &[u8],
    ) -> Result<Self, UqCheckpointError> {
        let _ = admit_plan(plan).map_err(UqCheckpointError::InvalidPlan)?;
        if proposal.germ_shift.len() != plan.parameters.len() {
            return Err(UqCheckpointError::InvalidPlan("importance checkpoint shift dimension differs"));
        }
        if bytes.len() < FIXED || bytes.len() > FIXED + RECORD * plan.budget_max_samples {
            return Err(UqCheckpointError::InvalidEncoding("length outside admitted importance budget"));
        }
        if &bytes[..8] != MAGIC {
            return Err(UqCheckpointError::InvalidEncoding("unknown importance checkpoint version"));
        }
        let encoded: [u8; 8] = bytes[41..HEADER].try_into()
            .map_err(|_| UqCheckpointError::InvalidEncoding("missing importance count"))?;
        let count = usize::try_from(u64::from_le_bytes(encoded))
            .map_err(|_| UqCheckpointError::InvalidEncoding("importance count overflows usize"))?;
        if count > plan.budget_max_samples || bytes.len() != FIXED + RECORD * count {
            return Err(UqCheckpointError::InvalidEncoding("importance count or trailing bytes"));
        }
        let status = match bytes[40] {
            0 => UqStatus::BudgetTruncated,
            1 => UqStatus::Cancelled,
            2 => UqStatus::Complete,
            _ => return Err(UqCheckpointError::InvalidEncoding("unknown or refused importance status")),
        };
        if (status == UqStatus::Complete) != (count == plan.budget_max_samples) {
            return Err(UqCheckpointError::InvalidEncoding("importance status disagrees with original budget"));
        }
        let payload_end = bytes.len() - 32;
        if hash_domain(DOMAIN, &bytes[..payload_end]).as_bytes().as_slice() != &bytes[payload_end..] {
            return Err(UqCheckpointError::IntegrityMismatch);
        }
        if identity(plan, &proposal, model).as_bytes().as_slice() != &bytes[8..40] {
            return Err(UqCheckpointError::IdentityMismatch);
        }
        let mut execution = Self::new(plan, proposal).map_err(UqCheckpointError::InvalidPlan)?;
        let epsilon = execution.config.defensive_fraction;
        let threshold = execution.plan.compliance_threshold.expect("admitted ceiling");
        for record in bytes[HEADER..payload_end].chunks_exact(RECORD) {
            let qoi_bits: [u8; 8] = record[..8].try_into()
                .map_err(|_| UqCheckpointError::InvalidEncoding("incomplete importance response"))?;
            let weight_bits: [u8; 8] = record[8..].try_into()
                .map_err(|_| UqCheckpointError::InvalidEncoding("incomplete importance weight"))?;
            let qoi = f64::from_bits(u64::from_le_bytes(qoi_bits));
            let scaled_weight = f64::from_bits(u64::from_le_bytes(weight_bits));
            let likelihood_weight = scaled_weight / epsilon;
            if !qoi.is_finite() || !scaled_weight.is_finite() || scaled_weight <= 0.0
                || scaled_weight > 1.0 || !likelihood_weight.is_finite()
                || likelihood_weight <= 0.0
                || (execution.unchanged && scaled_weight.to_bits() != epsilon.to_bits())
            {
                return Err(UqCheckpointError::InvalidEncoding("invalid importance response or likelihood weight"));
            }
            execution.values.push(ImportanceObservation {
                qoi, failed: qoi > threshold, likelihood_weight, scaled_weight,
            });
        }
        execution.attempted = count;
        execution.status = status;
        Ok(execution)
    }
}
