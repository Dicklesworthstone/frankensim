//! Retain the frozen predictor and original observation stream as one checkpoint.
//! The existing raw checkpoint supplies the checksum and all execution decoding;
//! its model identity additionally binds the exact coefficient header. This is
//! integrity, not authentication or proof that coefficients were chosen fairly.

use fs_blake3::{ContentHash, hash_domain};

use super::LinearControlVariate;
use crate::product_checkpoint::execution_identity;
use crate::{UqCheckpointError, UqExecution, UqPlan};

const MAGIC: &[u8; 8] = b"FSUQCV01";
const HEADER_LEN: usize = 16;

impl LinearControlVariate {
    /// Save this pre-sampling control together with an accepted raw MC prefix.
    ///
    /// The original raw checkpoint format and reduction order are unchanged.
    /// Its identity binds the supplied model AND exact coefficient bits, so a
    /// different predictor cannot be attached to the saved observations simply
    /// by replacing the header. Zero-sample checkpoints retain a paid nominal
    /// linearization before the first expensive sample is started.
    ///
    /// `model_identity` must also cover any producer metadata stored by a caller
    /// outside this payload (for example the nominal objective and residual).
    /// Store and restore only trusted checkpoints: checksums cannot establish
    /// physical validity, authenticity or independence from the sample data.
    ///
    /// # Errors
    /// Refuses a different exact plan or an execution the raw owner cannot save.
    pub fn checkpoint(
        &self,
        execution: &UqExecution,
        model_identity: ContentHash,
    ) -> Result<Vec<u8>, UqCheckpointError> {
        if execution_identity(execution.plan(), model_identity)
            != execution_identity(&self.plan, model_identity)
        {
            return Err(UqCheckpointError::IdentityMismatch);
        }
        let mut bytes = Vec::with_capacity(HEADER_LEN + self.gradient.len() * 8);
        bytes.extend_from_slice(MAGIC);
        bytes.extend_from_slice(&(self.gradient.len() as u64).to_le_bytes());
        for value in &self.gradient {
            bytes.extend_from_slice(&value.to_bits().to_le_bytes());
        }
        let identity = bound_model(model_identity, &bytes);
        bytes.extend_from_slice(&execution.checkpoint(identity)?);
        Ok(bytes)
    }

    /// Recover the exact original control and next sample ordinal, without a
    /// nominal solve, model evaluation, refit or additional random draw.
    ///
    /// Validates the current plan and bounded header before decoding observations.
    /// The existing raw owner still decides status, count, checksum and lifetime
    /// budget admission. Means are reconstructed from the declared marginals,
    /// never accepted from the payload. Complete prefixes remain terminal.
    ///
    /// # Errors
    /// Rejects raw-only checkpoints, changed plans/models, substituted coefficient
    /// headers, malformed/nonfinite coefficients and all raw checkpoint refusals.
    pub fn restore(
        plan: &UqPlan,
        model_identity: ContentHash,
        bytes: &[u8],
    ) -> Result<(UqExecution, Self), UqCheckpointError> {
        let fresh = UqExecution::new(plan).map_err(UqCheckpointError::InvalidPlan)?;
        // The admitted plan caps dimensions and samples. Leave at most 1 KiB
        // for the existing raw framing; its decoder checks the EXACT length.
        let prefix_len = HEADER_LEN + plan.parameters.len() * 8;
        if bytes.len() < prefix_len
            || bytes.len() > prefix_len + 1024 + plan.budget_max_samples * 8
        {
            return Err(UqCheckpointError::InvalidEncoding("controlled checkpoint length"));
        }
        if &bytes[..8] != MAGIC {
            return Err(UqCheckpointError::InvalidEncoding("unknown controlled checkpoint version"));
        }
        let count = u64::from_le_bytes(bytes[8..HEADER_LEN].try_into().expect("checked header"));
        if count != plan.parameters.len() as u64 {
            return Err(UqCheckpointError::InvalidEncoding("control coefficient count"));
        }
        let gradient: Vec<f64> = bytes[HEADER_LEN..prefix_len].chunks_exact(8)
            .map(|chunk| f64::from_bits(u64::from_le_bytes(chunk.try_into().expect("eight bytes"))))
            .collect();
        let control = fresh.freeze_linear_control_variate(&gradient)
            .map_err(|_| UqCheckpointError::InvalidEncoding("nonfinite control coefficient"))?;
        let identity = bound_model(model_identity, &bytes[..prefix_len]);
        let execution = UqExecution::restore(plan, identity, &bytes[prefix_len..])?;
        Ok((execution, control))
    }
}

fn bound_model(model: ContentHash, header: &[u8]) -> ContentHash {
    let mut bytes = Vec::with_capacity(32 + header.len());
    bytes.extend_from_slice(&model.0);
    bytes.extend_from_slice(header);
    hash_domain("org.frankensim.uq.frozen-linear-control.v1", &bytes)
}

#[cfg(test)]
mod tests;
