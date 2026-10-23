//! Wrap the canonical QMC prefix codec with retained out-of-order slots.
use super::QmcBatchExecution;
use crate::{QmcConfig, QmcExecution, UqCheckpointError, UqPlan, UqStatus};
use fs_blake3::{ContentHash, hash_domain};

const MAGIC: &[u8; 8] = b"FSQMB001";
const HEADER: usize = 24;
const MIN_SCALAR: usize = 81;
const FIXED: usize = HEADER + MIN_SCALAR + 32;
const DOMAIN: &str = "org.frankensim.uq.qmc-batch.checkpoint.v1";
fn invalid(message: &'static str) -> UqCheckpointError { UqCheckpointError::InvalidEncoding(message) }

impl QmcBatchExecution {
    /// Save exact response bits, the original layout, and each pending slot.
    /// The embedded QMC codec binds the plan, model and sampler semantics.
    /// Worker count may change. Checksums are not authentication or physical
    /// provenance; use trusted checkpoints on the same deterministic runtime.
    pub fn checkpoint(&self, model: ContentHash) -> Result<Vec<u8>, UqCheckpointError> {
        let report = self.report();
        if report.status == UqStatus::Refused
            || report.evaluations_attempted != report.evaluations_accepted
        { return Err(UqCheckpointError::NotResumable); }
        let scalar = self.scalar.checkpoint(model)?;
        let mut bytes = Vec::new();
        bytes.try_reserve_exact(HEADER + scalar.len() + 9 * self.pending.len() + 32)
            .map_err(|_| invalid("QMC batch checkpoint allocation refused"))?;
        bytes.extend_from_slice(MAGIC);
        bytes.extend_from_slice(&(scalar.len() as u64).to_le_bytes());
        bytes.extend_from_slice(&(self.pending.len() as u64).to_le_bytes());
        bytes.extend_from_slice(&scalar);
        for value in &self.pending {
            bytes.push(u8::from(value.is_some()));
            bytes.extend_from_slice(&value.unwrap_or(0.0).to_bits().to_le_bytes());
        }
        let checksum = hash_domain(DOMAIN, &bytes);
        bytes.extend_from_slice(checksum.as_bytes());
        Ok(bytes)
    }

    /// Re-admit the original QMC layout and restore without sampling or physics.
    /// The scalar codec owns its identity/status/value validation; this envelope
    /// validates the bounded pending window before attaching successful siblings.
    pub fn restore(
        plan: &UqPlan, config: QmcConfig, model: ContentHash, bytes: &[u8],
    ) -> Result<Self, UqCheckpointError> {
        if !(2..=1_000_000).contains(&plan.budget_max_samples) {
            return Err(UqCheckpointError::InvalidPlan("QMC batch budget is outside 2..=1000000"));
        }
        if bytes.len() < FIXED || bytes.len() > FIXED + 8 * plan.budget_max_samples + 256 {
            return Err(invalid("QMC batch length exceeds the observation envelope"));
        }
        if &bytes[..8] != MAGIC { return Err(invalid("unknown QMC batch checkpoint version")); }
        let count = |offset: usize| -> Result<usize, UqCheckpointError> {
            let raw: [u8; 8] = bytes[offset..offset + 8].try_into()
                .map_err(|_| invalid("missing QMC batch count"))?;
            usize::try_from(u64::from_le_bytes(raw)).map_err(|_| invalid("QMC batch count overflows usize"))
        };
        let scalar_len = count(8)?;
        let pending_len = count(16)?;
        if scalar_len < MIN_SCALAR || scalar_len > MIN_SCALAR + 8 * plan.budget_max_samples
            || pending_len > 256 || bytes.len() != HEADER + scalar_len + 9 * pending_len + 32
        { return Err(invalid("QMC batch framing or pending window differs")); }
        let end = bytes.len() - 32;
        if hash_domain(DOMAIN, &bytes[..end]).as_bytes().as_slice() != &bytes[end..] {
            return Err(UqCheckpointError::IntegrityMismatch);
        }
        let mut scalar = QmcExecution::restore(plan, config, model, &bytes[HEADER..HEADER + scalar_len])?;
        if pending_len > plan.budget_max_samples - scalar.values.len() {
            return Err(invalid("QMC pending window exceeds the original budget"));
        }
        let mut pending = Vec::new();
        pending.try_reserve_exact(pending_len).map_err(|_| invalid("cannot reserve QMC pending window"))?;
        for (index, raw) in bytes[HEADER + scalar_len..end].chunks_exact(9).enumerate() {
            let bits: [u8; 8] = raw[1..].try_into().map_err(|_| invalid("missing QMC pending value"))?;
            let bits = u64::from_le_bytes(bits);
            let value = f64::from_bits(bits);
            match raw[0] {
                0 if bits == 0 => pending.push(None),
                1 if index != 0 && value.is_finite() => pending.push(Some(value)),
                _ => return Err(invalid("invalid or noncanonical QMC pending slot")),
            }
        }
        scalar.attempted += pending.iter().filter(|value| value.is_some()).count();
        Ok(Self { scalar, pending })
    }
}
