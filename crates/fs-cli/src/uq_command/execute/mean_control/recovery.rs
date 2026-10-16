//! A retained nominal objective/residual bound to the frozen-control checkpoint.
//! No second observation codec: fs-uq owns coefficients, samples and checksum.

use super::{Failure, MeanControl, Result};
use fs_blake3::{ContentHash, hash_domain};
use fs_uq::{LinearControlVariate, UqExecution, UqPlan};

const MAGIC: &[u8; 8] = b"FSUQAD01";
const HEADER_LEN: usize = 24;

impl MeanControl {
    pub(crate) fn checkpoint_bytes(
        &self, execution: &UqExecution, identity: ContentHash,
    ) -> Result<Vec<u8>> {
        let mut bytes = Vec::with_capacity(HEADER_LEN);
        bytes.extend_from_slice(MAGIC);
        bytes.extend_from_slice(&self.nominal_temperature.to_bits().to_le_bytes());
        bytes.extend_from_slice(&self.adjoint_residual.to_bits().to_le_bytes());
        let bound = nominal_identity(identity, &bytes);
        bytes.extend_from_slice(&self.frozen.checkpoint(execution, bound)
            .map_err(|error| failure(error.to_string()))?);
        Ok(bytes)
    }

    /// The original nominal state is read, never recalculated from a new solve.
    /// Altering either diagnostic invalidates the inner plan/model identity.
    pub(crate) fn restore_bytes(
        plan: &UqPlan, identity: ContentHash, bytes: &[u8],
    ) -> Result<(UqExecution, Self)> {
        if bytes.len() < HEADER_LEN || &bytes[..8] != MAGIC {
            return Err(failure("expected a frozen-adjoint checkpoint; raw-only prefixes cannot acquire an adjoint after sampling"));
        }
        let number = |start: usize| f64::from_bits(u64::from_le_bytes(
            bytes[start..start + 8].try_into().expect("checked nominal header")));
        let nominal_temperature = number(8);
        let adjoint_residual = number(16);
        if !nominal_temperature.is_finite()
            || !(adjoint_residual.is_finite() && adjoint_residual >= 0.0)
        {
            return Err(failure("nonfinite nominal objective or invalid adjoint residual"));
        }
        let bound = nominal_identity(identity, &bytes[..HEADER_LEN]);
        let (execution, frozen) = LinearControlVariate::restore(plan, bound, &bytes[HEADER_LEN..])
            .map_err(|error| failure(error.to_string()))?;
        Ok((execution, Self { frozen, nominal_temperature, adjoint_residual }))
    }
}

fn nominal_identity(identity: ContentHash, header: &[u8]) -> ContentHash {
    let mut bytes = Vec::with_capacity(32 + header.len());
    bytes.extend_from_slice(&identity.0);
    bytes.extend_from_slice(header);
    hash_domain("org.frankensim.cooling-uq.nominal-adjoint.v1", &bytes)
}

fn failure(message: impl Into<String>) -> Failure {
    Failure { code: "cooling-network-uq-checkpoint", message: message.into() }
}
