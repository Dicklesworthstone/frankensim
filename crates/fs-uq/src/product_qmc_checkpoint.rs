//! Durable QMC prefixes. The sampler remains random-access: restoring the
//! ordered observations restores the next replicate/point without replay.

use fs_blake3::{ContentHash, hash_domain};

use super::{QmcConfig, QmcExecution};
use crate::product_checkpoint::execution_identity;
use crate::{UqCheckpointError, UqPlan, UqStatus};

const MAGIC: &[u8; 8] = b"FSQMC001";
const HEADER_LEN: usize = 8 + 32 + 1 + 8;
const FIXED_LEN: usize = HEADER_LEN + 32;
const CHECKSUM_DOMAIN: &str = "org.frankensim.uq.qmc.checkpoint.v1";

// Reuse the exact, length-framed plan/model identity, then bind the layout and
// QMC sampler semantics separately. A change to midpoint Sobol addressing,
// scramble keys, or the normal transform requires a new QMC identity version.
fn identity(plan: &UqPlan, config: QmcConfig, model: ContentHash) -> ContentHash {
    let mut bytes = Vec::with_capacity(48);
    bytes.extend_from_slice(&execution_identity(plan, model).0);
    bytes.extend_from_slice(&(config.replicates as u64).to_le_bytes());
    bytes.extend_from_slice(&(config.samples_per_replicate as u64).to_le_bytes());
    hash_domain("org.frankensim.uq.qmc.execution.v1", &bytes)
}

impl QmcExecution {
    /// Persist every accepted point, including an unfinished replicate.
    ///
    /// The caller's model identity must cover the evaluator implementation and
    /// all fixed physical inputs and bindings. IEEE-754 observation bits and
    /// the original layout/budget are retained. The checksum detects corruption;
    /// it does not authenticate a producer or certify the physical observations.
    /// Obtain checkpoints and model identities from a trusted source.
    ///
    /// # Errors
    /// Refused executions are terminal and cannot yield resumable checkpoints.
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
        // Shared admission limits the original budget to one million samples.
        let mut bytes = Vec::with_capacity(FIXED_LEN + self.values.len() * 8);
        bytes.extend_from_slice(MAGIC);
        bytes.extend_from_slice(&identity(&self.plan, self.config, model_identity).0);
        bytes.push(status);
        bytes.extend_from_slice(&(self.values.len() as u64).to_le_bytes());
        for value in &self.values {
            bytes.extend_from_slice(&value.to_bits().to_le_bytes());
        }
        let checksum = hash_domain(CHECKSUM_DOMAIN, &bytes);
        bytes.extend_from_slice(&checksum.0);
        Ok(bytes)
    }

    /// Resume at the exact next Sobol point without re-evaluating paid work.
    ///
    /// Re-admits the probability model and fixed layout, checks framing before
    /// allocating observations, and reconstructs reports from the ordered data.
    /// Incomplete nets remain excluded from replicate estimates. Completed runs
    /// remain terminal; cancelled and truncated runs retain their lifetime cap.
    /// Neither resumption nor repeated checkpoints licenses optional stopping.
    ///
    /// # Errors
    /// Refuses changed plans/layouts/models, MC checkpoints, unknown versions,
    /// malformed counts/statuses, corruption, and non-finite observations.
    pub fn restore(
        plan: &UqPlan,
        config: QmcConfig,
        model_identity: ContentHash,
        bytes: &[u8],
    ) -> Result<Self, UqCheckpointError> {
        let mut execution = Self::new(plan, config).map_err(UqCheckpointError::InvalidPlan)?;
        if bytes.len() < FIXED_LEN || bytes.len() > FIXED_LEN + plan.budget_max_samples * 8 {
            return Err(UqCheckpointError::InvalidEncoding("length outside admitted QMC budget"));
        }
        if &bytes[..8] != MAGIC {
            return Err(UqCheckpointError::InvalidEncoding("unknown QMC checkpoint version"));
        }
        let encoded: [u8; 8] = bytes[41..HEADER_LEN]
            .try_into()
            .map_err(|_| UqCheckpointError::InvalidEncoding("missing QMC observation count"))?;
        let count = usize::try_from(u64::from_le_bytes(encoded))
            .map_err(|_| UqCheckpointError::InvalidEncoding("QMC count overflows usize"))?;
        if count > plan.budget_max_samples || bytes.len() != FIXED_LEN + count * 8 {
            return Err(UqCheckpointError::InvalidEncoding("QMC observation count or trailing bytes"));
        }
        let status = match bytes[40] {
            0 => UqStatus::BudgetTruncated,
            1 => UqStatus::Cancelled,
            2 => UqStatus::Complete,
            _ => return Err(UqCheckpointError::InvalidEncoding("unknown or refused QMC status")),
        };
        if (status == UqStatus::Complete) != (count == plan.budget_max_samples) {
            return Err(UqCheckpointError::InvalidEncoding("QMC status disagrees with original budget"));
        }
        let payload_end = bytes.len() - 32;
        let checksum = hash_domain(CHECKSUM_DOMAIN, &bytes[..payload_end]);
        if checksum.0.as_slice() != &bytes[payload_end..] {
            return Err(UqCheckpointError::IntegrityMismatch);
        }
        if identity(plan, config, model_identity).0.as_slice() != &bytes[8..40] {
            return Err(UqCheckpointError::IdentityMismatch);
        }
        execution.values.reserve(count);
        for chunk in bytes[HEADER_LEN..payload_end].chunks_exact(8) {
            let encoded: [u8; 8] = chunk
                .try_into()
                .map_err(|_| UqCheckpointError::InvalidEncoding("incomplete QMC observation"))?;
            let value = f64::from_bits(u64::from_le_bytes(encoded));
            if !value.is_finite() {
                return Err(UqCheckpointError::InvalidEncoding("non-finite QMC observation"));
            }
            execution.values.push(value);
        }
        execution.attempted = count;
        execution.status = status;
        // QMC uses its own scaled replicate statistics, not MC pointwise
        // variance or confidence-sequence admission. Finite extreme values
        // and partial nets must retain their existing QMC semantics.
        Ok(execution)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{CorrelationModel, ParameterUncertainty, PropagationMethod, UqExecution};

    fn plan() -> UqPlan {
        UqPlan::new("junction", PropagationMethod::QuasiMonteCarlo, 32)
            .with_parameter(ParameterUncertainty::uniform("ambient", 290.0, 310.0, "K"))
            .with_parameter(ParameterUncertainty::gaussian("power", 50.0, 2.0, "W"))
            .with_correlation(CorrelationModel::Independent)
            .with_compliance_threshold(345.0)
    }

    fn layout() -> QmcConfig { QmcConfig { replicates: 4, samples_per_replicate: 8 } }
    fn model() -> ContentHash { hash_domain("test-qmc-model", b"T = ambient + 0.8 * power") }
    fn qoi(x: &[f64]) -> f64 { x[0] + 0.8 * x[1] }

    fn prefix(count: usize) -> QmcExecution {
        let mut execution = QmcExecution::new(&plan(), layout()).unwrap();
        execution.advance(count, || false, |x| Ok::<_, &str>(qoi(x)));
        execution
    }

    fn reseal(bytes: &mut [u8]) {
        let end = bytes.len() - 32;
        let checksum = hash_domain(CHECKSUM_DOMAIN, &bytes[..end]);
        bytes[end..].copy_from_slice(&checksum.0);
    }

    #[test]
    fn every_split_preserves_the_exact_suffix_and_final_report() {
        let full = prefix(32);
        for split in 0..=32 {
            let partial = prefix(split);
            let bytes = partial.checkpoint(model()).unwrap();
            let mut restored = QmcExecution::restore(&plan(), layout(), model(), &bytes).unwrap();
            assert_eq!(restored.report(), partial.report());
            assert_eq!(restored.report().completed_replicates, split / 8);
            let mut calls = 0;
            let result = restored.advance(usize::MAX, || false, |x| {
                calls += 1;
                Ok::<_, &str>(qoi(x))
            });
            assert_eq!(calls, 32 - split);
            assert_eq!(result, full.report());
            assert_eq!(restored.checkpoint(model()), full.checkpoint(model()));
        }
    }

    #[test]
    fn interrupted_point_is_retried_after_durable_restore() {
        let mut execution = prefix(11);
        let mut interrupted = Vec::new();
        execution.advance_interruptible(1, || false, |x| {
            interrupted = x.to_vec();
            Ok::<_, &str>(None)
        });
        let mut restored = QmcExecution::restore(
            &plan(), layout(), model(), &execution.checkpoint(model()).unwrap(),
        ).unwrap();
        assert_eq!(restored.report().status, UqStatus::Cancelled);
        assert_eq!(restored.evaluations_attempted(), 11);
        restored.advance(1, || false, |x| {
            assert_eq!(x, interrupted.as_slice());
            Ok::<_, &str>(qoi(x))
        });
        assert_eq!(restored.report(), prefix(12).report());
    }

    #[test]
    fn model_plan_and_equal_total_but_different_layout_cannot_change() {
        let bytes = prefix(3).checkpoint(model()).unwrap();
        let other_layout = QmcConfig { replicates: 8, samples_per_replicate: 4 };
        assert_eq!(QmcExecution::restore(&plan(), other_layout, model(), &bytes).unwrap_err(), UqCheckpointError::IdentityMismatch);
        assert_eq!(QmcExecution::restore(&plan(), layout(), ContentHash([0; 32]), &bytes).unwrap_err(), UqCheckpointError::IdentityMismatch);
        let mut changes = Vec::new();
        let mut p = plan(); p.seed += 1; changes.push(p);
        let mut p = plan(); p.parameters.swap(0, 1); changes.push(p);
        let mut p = plan(); p.parameters[0].unit.push('x'); changes.push(p);
        let mut p = plan(); p.target_qoi.push('x'); changes.push(p);
        let mut p = plan(); p.compliance_threshold = None; changes.push(p);
        for changed in changes {
            assert_eq!(QmcExecution::restore(&changed, layout(), model(), &bytes).unwrap_err(), UqCheckpointError::IdentityMismatch);
        }
    }

    #[test]
    fn corrupted_truncated_and_forged_observations_refuse() {
        let bytes = prefix(3).checkpoint(model()).unwrap();
        for end in 0..bytes.len() {
            assert!(QmcExecution::restore(&plan(), layout(), model(), &bytes[..end]).is_err());
        }
        let mut extra = bytes.clone(); extra.push(0);
        assert!(QmcExecution::restore(&plan(), layout(), model(), &extra).is_err());
        for index in [0, 8, 40, 41, HEADER_LEN, bytes.len() - 1] {
            let mut bad = bytes.clone(); bad[index] ^= 0x80;
            assert!(QmcExecution::restore(&plan(), layout(), model(), &bad).is_err());
        }
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut bad = bytes.clone();
            bad[HEADER_LEN..HEADER_LEN + 8].copy_from_slice(&value.to_bits().to_le_bytes());
            reseal(&mut bad);
            assert!(QmcExecution::restore(&plan(), layout(), model(), &bad).is_err());
        }
        let mut bad = bytes.clone(); bad[40] = 2; reseal(&mut bad);
        assert!(QmcExecution::restore(&plan(), layout(), model(), &bad).is_err());
        let mut bad = bytes; bad[41..HEADER_LEN].copy_from_slice(&u64::MAX.to_le_bytes()); reseal(&mut bad);
        assert!(QmcExecution::restore(&plan(), layout(), model(), &bad).is_err());
    }

    #[test]
    fn failed_runs_cannot_be_laundered_into_resumable_prefixes() {
        let mut execution = prefix(9);
        execution.advance(1, || false, |_| Err::<f64, _>("solver refused"));
        assert_eq!(execution.checkpoint(model()), Err(UqCheckpointError::NotResumable));
        assert_eq!(execution.evaluations_attempted(), 10);
    }

    #[test]
    fn extreme_finite_values_and_signed_zero_keep_qmc_not_mc_statistics() {
        for value in [-0.0_f64, f64::MAX, -f64::MAX] {
            let mut execution = QmcExecution::new(&plan(), layout()).unwrap();
            execution.advance(32, || false, |_| Ok::<_, &str>(value));
            let bytes = execution.checkpoint(model()).unwrap();
            let mut restored = QmcExecution::restore(&plan(), layout(), model(), &bytes).unwrap();
            assert!(restored.observations().iter().all(|x| x.to_bits() == value.to_bits()));
            assert_eq!(restored.report(), execution.report());
            restored.advance(1, || panic!("complete"), |_| -> Result<f64, &str> { panic!("complete") });
        }
    }

    #[test]
    fn mc_and_qmc_checkpoints_are_not_interchangeable() {
        let mut mc_plan = plan(); mc_plan.method = PropagationMethod::MonteCarlo;
        let mc = UqExecution::new(&mc_plan).unwrap().checkpoint(model()).unwrap();
        assert!(QmcExecution::restore(&plan(), layout(), model(), &mc).is_err());
        assert!(UqExecution::restore(&mc_plan, model(), &prefix(0).checkpoint(model()).unwrap()).is_err());
    }
}
