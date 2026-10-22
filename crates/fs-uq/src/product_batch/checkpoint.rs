//! Original-ordinal recovery, including successful siblings behind unfinished calls.
use super::UqBatchExecution;
use fs_blake3::{ContentHash, hash_domain};
use crate::product_checkpoint::execution_identity;
use crate::product_plan::admit_plan;
use crate::{UqCheckpointError, UqPlan, UqStatus};

const MAGIC: &[u8; 8] = b"FSUQBT01";
const HEADER: usize = 57;
const FIXED: usize = HEADER + 32;
const DOMAIN: &str = "org.frankensim.uq.batch.checkpoint.v1";
fn invalid(reason: &'static str) -> UqCheckpointError { UqCheckpointError::InvalidEncoding(reason) }
fn identity(plan: &UqPlan, model: ContentHash) -> ContentHash {
    hash_domain("org.frankensim.uq.batch.execution.v1", execution_identity(plan, model).as_bytes())
}

impl UqBatchExecution {
    /// Store raw response bits and pending ordinal slots, not derived statistics.
    /// Model identity must cover the evaluator and all fixed physical inputs.
    /// Worker count is NOT part of the law and may change after restoration.
    /// Checksums detect corruption, not authenticity: restore only trusted data.
    pub fn checkpoint(&self, model: ContentHash) -> Result<Vec<u8>, UqCheckpointError> {
        let report = self.report();
        let status = match report.status {
            UqStatus::BudgetTruncated => 0,
            UqStatus::Cancelled => 1,
            UqStatus::Complete => 2,
            UqStatus::Refused => return Err(UqCheckpointError::NotResumable),
        };
        if self.scalar.failure.is_some() || report.evaluations_attempted != report.evaluations_accepted {
            return Err(UqCheckpointError::NotResumable);
        }
        let mut bytes = Vec::new();
        bytes.try_reserve_exact(FIXED + 8 * self.scalar.values.len() + 9 * self.pending.len())
            .map_err(|_| invalid("batch checkpoint allocation refused"))?;
        bytes.extend_from_slice(MAGIC);
        bytes.extend_from_slice(identity(self.plan(), model).as_bytes());
        bytes.push(status);
        bytes.extend_from_slice(&(self.scalar.values.len() as u64).to_le_bytes());
        bytes.extend_from_slice(&(self.pending.len() as u64).to_le_bytes());
        for &x in &self.scalar.values { bytes.extend_from_slice(&x.to_bits().to_le_bytes()); }
        for x in &self.pending {
            bytes.push(u8::from(x.is_some()));
            bytes.extend_from_slice(&x.unwrap_or(0.0).to_bits().to_le_bytes());
        }
        let checksum = hash_domain(DOMAIN, &bytes);
        bytes.extend_from_slice(checksum.as_bytes());
        Ok(bytes)
    }

    /// Re-admit the original plan; validate framing and identity before reserving
    /// observations. No physics or completed sibling is evaluated again. Final
    /// statistics use the original scalar reducer, not a serialized estimate.
    pub fn restore(plan: &UqPlan, model: ContentHash, bytes: &[u8]) -> Result<Self, UqCheckpointError> {
        let _ = admit_plan(plan).map_err(UqCheckpointError::InvalidPlan)?;
        if bytes.len() < FIXED || bytes.len() > FIXED + 8 * plan.budget_max_samples + 256 {
            return Err(invalid("length outside batch observation envelope"));
        }
        if &bytes[..8] != MAGIC { return Err(invalid("unknown batch checkpoint version")); }
        let read_count = |start: usize| -> Result<usize, UqCheckpointError> {
            let raw: [u8; 8] = bytes[start..start + 8].try_into().map_err(|_| invalid("missing batch count"))?;
            usize::try_from(u64::from_le_bytes(raw)).map_err(|_| invalid("batch count overflows usize"))
        };
        let n = read_count(41)?;
        let pending = read_count(49)?;
        if n > plan.budget_max_samples || pending > 256 || pending > plan.budget_max_samples - n
            || bytes.len() != FIXED + n * 8 + pending * 9
        { return Err(invalid("batch count, window or trailing bytes")); }
        let status = match bytes[40] {
            0 => UqStatus::BudgetTruncated,
            1 => UqStatus::Cancelled,
            2 => UqStatus::Complete,
            _ => return Err(invalid("unknown or refused batch status")),
        };
        if (status == UqStatus::Complete) != (n == plan.budget_max_samples) {
            return Err(invalid("batch status disagrees with original budget"));
        }
        let end = bytes.len() - 32;
        if hash_domain(DOMAIN, &bytes[..end]).as_bytes().as_slice() != &bytes[end..] {
            return Err(UqCheckpointError::IntegrityMismatch);
        }
        if identity(plan, model).as_bytes().as_slice() != &bytes[8..40] {
            return Err(UqCheckpointError::IdentityMismatch);
        }
        let mut execution = Self::new(plan).map_err(UqCheckpointError::InvalidPlan)?;
        let value = |raw: &[u8]| -> Result<f64, UqCheckpointError> {
            let bits: [u8; 8] = raw.try_into().map_err(|_| invalid("incomplete batch value"))?;
            let x = f64::from_bits(u64::from_le_bytes(bits));
            if x.is_finite() { Ok(x) } else { Err(invalid("non-finite batch observation")) }
        };
        for raw in bytes[HEADER..HEADER + 8 * n].chunks_exact(8) { execution.scalar.values.push(value(raw)?); }
        for (i, raw) in bytes[HEADER + 8 * n..end].chunks_exact(9).enumerate() {
            let x = value(&raw[1..])?;
            match raw[0] {
                0 if x.to_bits() == 0 => execution.pending.push(None),
                1 if i != 0 => execution.pending.push(Some(x)),
                _ => return Err(invalid("invalid or noncanonical pending batch slot")),
            }
        }
        execution.scalar.attempted = n + execution.pending.iter().filter(|x| x.is_some()).count();
        execution.scalar.status = status;
        if execution.report().status == UqStatus::Refused {
            return Err(invalid("unrepresentable completed batch statistics"));
        }
        Ok(execution)
    }
}
