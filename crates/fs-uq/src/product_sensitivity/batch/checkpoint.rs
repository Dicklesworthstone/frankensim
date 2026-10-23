//! Original Sobol checkpoint plus retained out-of-order physical responses.
use super::{SobolBatchExecution, SobolExecution};
use crate::product_plan::admit_plan;
use crate::{UqCheckpointError, UqPlan, UqStatus};
use fs_blake3::{ContentHash, hash_domain};

const MAGIC: &[u8; 8] = b"FSSOBB01";
const HEADER: usize = 24;
const INNER_FIXED: usize = 81;
const FIXED: usize = HEADER + INNER_FIXED + 32;
const DOMAIN: &str = "org.frankensim.uq.sobol-batch.checkpoint.v1";
fn invalid(reason: &'static str) -> UqCheckpointError { UqCheckpointError::InvalidEncoding(reason) }

impl SobolBatchExecution {
    /// Persist exact response bits, not derived indices. The embedded original
    /// codec binds every law, budget, pairing and model-identity field. Worker
    /// count is not a model parameter. Checksums detect corruption, not correct
    /// physics or authenticity: restore only trusted files on the same runtime.
    pub fn checkpoint(&self, model: ContentHash) -> Result<Vec<u8>, UqCheckpointError> {
        let accepted = self.scalar.values.len() + self.pending.iter().filter(|x| x.is_some()).count();
        if self.scalar.failure.is_some() || self.scalar.attempted != accepted {
            return Err(UqCheckpointError::NotResumable);
        }
        let inner = self.scalar.checkpoint(model)?;
        let mut bytes = Vec::new();
        bytes.try_reserve_exact(HEADER + inner.len() + 9 * self.pending.len() + 32)
            .map_err(|_| invalid("Sobol batch checkpoint allocation refused"))?;
        bytes.extend_from_slice(MAGIC);
        bytes.extend_from_slice(&(inner.len() as u64).to_le_bytes());
        bytes.extend_from_slice(&(self.pending.len() as u64).to_le_bytes());
        bytes.extend_from_slice(&inner);
        for value in &self.pending {
            bytes.push(u8::from(value.is_some()));
            bytes.extend_from_slice(&value.unwrap_or(0.0).to_bits().to_le_bytes());
        }
        let checksum = hash_domain(DOMAIN, &bytes);
        bytes.extend_from_slice(checksum.as_bytes());
        Ok(bytes)
    }

    /// Validate bounded framing and both codecs before retaining the saved
    /// design. Restoration invokes no physical solver or parameter sampler.
    /// Only missing ordinals are sampled on advance. A partial A/B/hybrid row
    /// remains work, not an estimate; completed undefined-variance reports keep
    /// the original Sobol semantics. Serial and batch formats are distinct.
    pub fn restore(plan: &UqPlan, model: ContentHash, bytes: &[u8]) -> Result<Self, UqCheckpointError> {
        let _ = admit_plan(plan).map_err(UqCheckpointError::InvalidPlan)?;
        if bytes.len() < FIXED || bytes.len() > FIXED + 8 * plan.budget_max_samples + 256 {
            return Err(invalid("length outside Sobol batch envelope"));
        }
        if &bytes[..8] != MAGIC { return Err(invalid("unknown Sobol batch checkpoint version")); }
        let count = |start: usize| -> Result<usize, UqCheckpointError> {
            let raw: [u8; 8] = bytes[start..start + 8].try_into()
                .map_err(|_| invalid("missing Sobol batch count"))?;
            usize::try_from(u64::from_le_bytes(raw)).map_err(|_| invalid("Sobol batch count overflows usize"))
        };
        let inner_len = count(8)?;
        let pending_len = count(16)?;
        if !(INNER_FIXED..=INNER_FIXED + 8 * plan.budget_max_samples).contains(&inner_len)
            || pending_len > 256 || bytes.len() != HEADER + inner_len + 9 * pending_len + 32
        { return Err(invalid("Sobol batch counts or trailing bytes")); }
        let end = bytes.len() - 32;
        if hash_domain(DOMAIN, &bytes[..end]).as_bytes().as_slice() != &bytes[end..] {
            return Err(UqCheckpointError::IntegrityMismatch);
        }
        let inner_end = HEADER + inner_len;
        let mut scalar = SobolExecution::restore(plan, model, &bytes[HEADER..inner_end])?;
        if pending_len > plan.budget_max_samples - scalar.values.len()
            || (scalar.status == UqStatus::Complete && pending_len != 0)
        { return Err(invalid("pending Sobol values exceed original design")); }
        let mut pending = Vec::new();
        pending.try_reserve_exact(pending_len).map_err(|_| invalid("cannot reserve Sobol pending observations"))?;
        for (i, record) in bytes[inner_end..end].chunks_exact(9).enumerate() {
            let raw: [u8; 8] = record[1..].try_into().map_err(|_| invalid("incomplete Sobol pending value"))?;
            let value = f64::from_bits(u64::from_le_bytes(raw));
            match record[0] {
                0 if value.to_bits() == 0 => pending.push(None),
                1 if i != 0 && value.is_finite() => pending.push(Some(value)),
                _ => return Err(invalid("noncanonical or nonfinite Sobol pending observation")),
            }
        }
        scalar.attempted += pending.iter().filter(|x| x.is_some()).count();
        Ok(Self { scalar, pending })
    }
}
