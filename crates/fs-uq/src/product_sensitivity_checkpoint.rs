//! Durable row-major A/B/hybrid observations for direct-model sensitivity.
//! Only physical evaluations are retained; Philox regenerates the next input.

use fs_blake3::{ContentHash, hash_domain};

use super::SobolExecution;
use crate::product_checkpoint::execution_identity;
use crate::product_plan::admit_plan;
use crate::{UqCheckpointError, UqPlan, UqStatus};

const MAGIC: &[u8; 8] = b"FSSOB001";
const HEADER_LEN: usize = 8 + 32 + 1 + 8;
const FIXED_LEN: usize = HEADER_LEN + 32;
const CHECKSUM_DOMAIN: &str = "org.frankensim.uq.sobol-sensitivity.checkpoint.v1";

fn identity(plan: &UqPlan, model: ContentHash) -> ContentHash {
    // The shared identity binds every plan field and Philox sampler semantics.
    // This additional domain binds A/B ordinals 2*r, 2*r+1 and the declared
    // coordinate order of the d hybrids. Width and row count follow the plan.
    hash_domain(
        "org.frankensim.uq.sobol-sensitivity.execution.v1",
        execution_identity(plan, model).as_bytes(),
    )
}

impl SobolExecution {
    /// Save accepted observations including a partly completed pick-freeze row.
    ///
    /// `model_identity` must cover the evaluator and all fixed inputs/bindings.
    /// The original plan, budget, pairing and IEEE-754 observation bits remain
    /// fixed. The checksum detects corruption, not authenticity or correctness
    /// of the physical observations; only restore checkpoints from a trusted
    /// source. A failed evaluation is terminal and cannot become a valid prefix.
    pub fn checkpoint(&self, model_identity: ContentHash) -> Result<Vec<u8>, UqCheckpointError> {
        let status = match self.status {
            UqStatus::BudgetTruncated => 0,
            UqStatus::Cancelled => 1,
            UqStatus::Complete => 2,
            UqStatus::Refused => return Err(UqCheckpointError::NotResumable),
        };
        if self.failure.is_some() {
            return Err(UqCheckpointError::NotResumable);
        }
        let mut bytes = Vec::with_capacity(FIXED_LEN + self.values.len() * 8);
        bytes.extend_from_slice(MAGIC);
        bytes.extend_from_slice(identity(&self.plan, model_identity).as_bytes());
        bytes.push(status);
        bytes.extend_from_slice(&(self.values.len() as u64).to_le_bytes());
        for value in &self.values {
            bytes.extend_from_slice(&value.to_bits().to_le_bytes());
        }
        let checksum = hash_domain(CHECKSUM_DOMAIN, &bytes);
        bytes.extend_from_slice(checksum.as_bytes());
        Ok(bytes)
    }

    /// Restore the exact next A, B or hybrid input without repeating paid work.
    ///
    /// Re-admits the fixed probability law, validates bounded framing before
    /// allocating observation storage, and recomputes statistics from the
    /// ordered observations. Partial designs produce no sensitivity estimate.
    /// Complete checkpoints remain terminal even when normalization is undefined.
    /// Changed plans/models, other sampling formats, corruption and malformed or
    /// non-finite observations refuse before any evaluator call.
    pub fn restore(
        plan: &UqPlan,
        model_identity: ContentHash,
        bytes: &[u8],
    ) -> Result<Self, UqCheckpointError> {
        let _ = admit_plan(plan).map_err(UqCheckpointError::InvalidPlan)?;
        if bytes.len() < FIXED_LEN || bytes.len() > FIXED_LEN + plan.budget_max_samples * 8 {
            return Err(UqCheckpointError::InvalidEncoding(
                "length outside admitted Sobol budget",
            ));
        }
        if &bytes[..8] != MAGIC {
            return Err(UqCheckpointError::InvalidEncoding(
                "unknown Sobol sensitivity checkpoint version",
            ));
        }
        let encoded: [u8; 8] = bytes[41..HEADER_LEN]
            .try_into()
            .map_err(|_| UqCheckpointError::InvalidEncoding("missing Sobol observation count"))?;
        let count = usize::try_from(u64::from_le_bytes(encoded))
            .map_err(|_| UqCheckpointError::InvalidEncoding("Sobol count overflows usize"))?;
        if count > plan.budget_max_samples || bytes.len() != FIXED_LEN + count * 8 {
            return Err(UqCheckpointError::InvalidEncoding(
                "Sobol observation count or trailing bytes",
            ));
        }
        let status = match bytes[40] {
            0 => UqStatus::BudgetTruncated,
            1 => UqStatus::Cancelled,
            2 => UqStatus::Complete,
            _ => {
                return Err(UqCheckpointError::InvalidEncoding(
                    "unknown or refused Sobol status",
                ));
            }
        };
        if (status == UqStatus::Complete) != (count == plan.budget_max_samples) {
            return Err(UqCheckpointError::InvalidEncoding(
                "Sobol status disagrees with original budget",
            ));
        }
        let payload_end = bytes.len() - 32;
        let checksum = hash_domain(CHECKSUM_DOMAIN, &bytes[..payload_end]);
        if checksum.as_bytes().as_slice() != &bytes[payload_end..] {
            return Err(UqCheckpointError::IntegrityMismatch);
        }
        if identity(plan, model_identity).as_bytes().as_slice() != &bytes[8..40] {
            return Err(UqCheckpointError::IdentityMismatch);
        }
        let mut execution = Self::new(plan).map_err(UqCheckpointError::InvalidPlan)?;
        for chunk in bytes[HEADER_LEN..payload_end].chunks_exact(8) {
            let encoded: [u8; 8] = chunk
                .try_into()
                .map_err(|_| UqCheckpointError::InvalidEncoding("incomplete Sobol observation"))?;
            let value = f64::from_bits(u64::from_le_bytes(encoded));
            if !value.is_finite() {
                return Err(UqCheckpointError::InvalidEncoding(
                    "non-finite Sobol observation",
                ));
            }
            execution.values.push(value);
        }
        execution.attempted = count;
        execution.status = status;
        Ok(execution)
    }
}
