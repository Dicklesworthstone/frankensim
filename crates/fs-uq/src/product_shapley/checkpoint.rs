//! Durable conditional pairs, including paid halves of unfinished coalitions.
use fs_blake3::{ContentHash, hash_domain};

use super::ShapleyExecution;
use crate::product_checkpoint::execution_identity;
use crate::product_plan::admit_plan;
use crate::{UqCheckpointError, UqPlan, UqStatus};

const MAGIC: &[u8; 8] = b"FSSHP001";
const HEADER: usize = 8 + 32 + 1 + 8;
const FIXED: usize = HEADER + 32;
const CHECKSUM: &str = "org.frankensim.uq.shapley.checkpoint.v1";

fn identity(plan: &UqPlan, model: ContentHash) -> ContentHash {
    // Bind original physical order, every permuted conditional law, row germs
    // 2*r/2*r+1, subset-mask/half ordering, and the canonical sampler version.
    hash_domain("org.frankensim.uq.shapley.execution.v1",
        execution_identity(plan, model).as_bytes())
}

impl ShapleyExecution {
    /// Persist all accepted physical responses in their original conditional-pair
    /// order. A half-pair is paid work, not an incomplete observation to discard.
    /// Model identity must cover all fixed inputs, bindings and evaluator semantics.
    /// Checksums detect corruption, not authenticity or physical correctness.
    /// Refused executions cannot publish a reusable selected prefix.
    pub fn checkpoint(&self, model: ContentHash) -> Result<Vec<u8>, UqCheckpointError> {
        let status = match self.status {
            UqStatus::BudgetTruncated => 0,
            UqStatus::Cancelled => 1,
            UqStatus::Complete => 2,
            UqStatus::Refused => return Err(UqCheckpointError::NotResumable),
        };
        if self.failure.is_some() { return Err(UqCheckpointError::NotResumable); }
        let mut bytes = Vec::with_capacity(FIXED + self.values.len() * 8);
        bytes.extend_from_slice(MAGIC);
        bytes.extend_from_slice(identity(&self.plan, model).as_bytes());
        bytes.push(status);
        bytes.extend_from_slice(&(self.values.len() as u64).to_le_bytes());
        for value in &self.values { bytes.extend_from_slice(&value.to_bits().to_le_bytes()); }
        let checksum = hash_domain(CHECKSUM, &bytes);
        bytes.extend_from_slice(checksum.as_bytes());
        Ok(bytes)
    }

    /// Restore without running physics. The next physical call is reconstructed
    /// from the original ordinal, even halfway through a conditional pair.
    /// Re-admit every permuted Gaussian law; do not trust serialized factors,
    /// cached inputs or estimates. Only raw finite response bits are retained.
    /// Length/count/status checks precede observation allocation. The original
    /// one-million-call cap bounds the payload. Completed zero-variance studies
    /// stay completed with undefined effects; they cannot acquire extra samples.
    /// Restore only trusted bytes on the same deterministic runtime profile.
    pub fn restore(plan: &UqPlan, model: ContentHash, bytes: &[u8])
        -> Result<Self, UqCheckpointError>
    {
        let _ = admit_plan(plan).map_err(UqCheckpointError::InvalidPlan)?;
        if bytes.len() < FIXED || bytes.len() > FIXED + plan.budget_max_samples * 8 {
            return Err(UqCheckpointError::InvalidEncoding("length outside Shapley budget"));
        }
        if &bytes[..8] != MAGIC {
            return Err(UqCheckpointError::InvalidEncoding("unknown Shapley checkpoint version"));
        }
        let encoded: [u8; 8] = bytes[41..HEADER].try_into()
            .map_err(|_| UqCheckpointError::InvalidEncoding("missing Shapley count"))?;
        let count = usize::try_from(u64::from_le_bytes(encoded))
            .map_err(|_| UqCheckpointError::InvalidEncoding("Shapley count overflows usize"))?;
        if count > plan.budget_max_samples || bytes.len() != FIXED + count * 8 {
            return Err(UqCheckpointError::InvalidEncoding("Shapley count or trailing bytes"));
        }
        let status = match bytes[40] {
            0 => UqStatus::BudgetTruncated,
            1 => UqStatus::Cancelled,
            2 => UqStatus::Complete,
            _ => return Err(UqCheckpointError::InvalidEncoding("unknown or refused Shapley status")),
        };
        if (status == UqStatus::Complete) != (count == plan.budget_max_samples) {
            return Err(UqCheckpointError::InvalidEncoding("Shapley status disagrees with budget"));
        }
        let end = bytes.len() - 32;
        if hash_domain(CHECKSUM, &bytes[..end]).as_bytes().as_slice() != &bytes[end..] {
            return Err(UqCheckpointError::IntegrityMismatch);
        }
        if identity(plan, model).as_bytes().as_slice() != &bytes[8..40] {
            return Err(UqCheckpointError::IdentityMismatch);
        }
        let mut execution = Self::new(plan).map_err(UqCheckpointError::InvalidPlan)?;
        for chunk in bytes[HEADER..end].chunks_exact(8) {
            let encoded: [u8; 8] = chunk.try_into()
                .map_err(|_| UqCheckpointError::InvalidEncoding("incomplete Shapley value"))?;
            let value = f64::from_bits(u64::from_le_bytes(encoded));
            if !value.is_finite() {
                return Err(UqCheckpointError::InvalidEncoding("non-finite Shapley observation"));
            }
            execution.values.push(value);
        }
        execution.attempted = count;
        execution.status = status;
        // Empty deterministic caches regenerate only input coordinates, never
        // completed physics, under the admitted row/subset/conditional-half map.
        Ok(execution)
    }
}
