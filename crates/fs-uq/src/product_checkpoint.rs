//! Durable checkpoints for the product Monte Carlo executor.
//!
//! The envelope binds exact plan fields, original lifetime budget, sampler
//! semantics, caller-supplied model identity, status, and ordered observations.
//! BLAKE3 detects accidental corruption; it does NOT authenticate the producer
//! or prove that an observation came from the named physical model. The caller
//! must obtain checkpoints and model identities from its trusted ledger.

use core::fmt;

use fs_blake3::{ContentHash, hash_domain};

use crate::product_plan::admit_plan;
use crate::{CorrelationModel, PropagationMethod, UncertaintyKind, UqExecution, UqPlan, UqStatus};

const MAGIC: &[u8; 8] = b"FSUQCP01";
const HEADER_LEN: usize = 8 + 32 + 1 + 8;
const FIXED_LEN: usize = HEADER_LEN + 32;
const CHECKSUM_DOMAIN: &str = "org.frankensim.uq.checkpoint.v1";

/// A checkpoint was refused before any model evaluation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum UqCheckpointError {
    /// The expected plan does not satisfy the current executor's admission rules.
    InvalidPlan(&'static str),
    /// Failed executions are terminal and cannot be checkpointed for resumption.
    NotResumable,
    /// The checkpoint's framing, count, status, or observations are invalid.
    InvalidEncoding(&'static str),
    /// The exact plan, original budget, sampler, or model identity differs.
    IdentityMismatch,
    /// The retained payload does not match its checksum.
    IntegrityMismatch,
}

impl fmt::Display for UqCheckpointError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidPlan(reason) => write!(formatter, "UQ checkpoint plan refused: {reason}"),
            Self::NotResumable => formatter.write_str("refused UQ executions cannot resume"),
            Self::InvalidEncoding(reason) => write!(formatter, "invalid UQ checkpoint: {reason}"),
            Self::IdentityMismatch => formatter.write_str("UQ checkpoint plan or model identity differs"),
            Self::IntegrityMismatch => formatter.write_str("UQ checkpoint checksum differs"),
        }
    }
}

impl std::error::Error for UqCheckpointError {}

impl UqExecution {
    /// Encode a versioned checkpoint for storage in a file or ledger artifact.
    ///
    /// `model_identity` must cover the evaluator implementation AND all fixed
    /// inputs (geometry, material cards, parameter bindings, solver tolerances,
    /// and any other state affecting the QoI). It is caller-supplied, not inferred
    /// from the closure. Resumption requires this same identity and exact plan.
    /// IEEE-754 observation bits, including signed zero, are preserved verbatim.
    ///
    /// # Errors
    /// Refused executions cannot resume and return `NotResumable`; retain their
    /// `report()` separately. No failure is turned into a resumable valid prefix.
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
        bytes.extend_from_slice(&execution_identity(&self.plan, model_identity).0);
        bytes.push(status);
        bytes.extend_from_slice(&(self.values.len() as u64).to_le_bytes());
        for value in &self.values {
            bytes.extend_from_slice(&value.to_bits().to_le_bytes());
        }
        let checksum = hash_domain(CHECKSUM_DOMAIN, &bytes);
        bytes.extend_from_slice(&checksum.0);
        Ok(bytes)
    }

    /// Restore the next sample ordinal without re-evaluating the retained prefix.
    ///
    /// The plan is re-admitted and lengths/counts are checked BEFORE allocating
    /// the observation vector. Statistics are recomputed, never trusted from the
    /// payload. Completed checkpoints remain terminal. Cancelled and truncated
    /// checkpoints continue under the original budget with `advance`.
    ///
    /// # Errors
    /// Rejects changed plans/models, unknown versions, truncation/trailing bytes,
    /// count/status inconsistencies, corruption, and non-finite observations.
    pub fn restore(
        plan: &UqPlan,
        model_identity: ContentHash,
        bytes: &[u8],
    ) -> Result<Self, UqCheckpointError> {
        let factor = admit_plan(plan).map_err(UqCheckpointError::InvalidPlan)?;
        if bytes.len() < FIXED_LEN || bytes.len() > FIXED_LEN + plan.budget_max_samples * 8 {
            return Err(UqCheckpointError::InvalidEncoding("length outside admitted budget"));
        }
        if &bytes[..8] != MAGIC {
            return Err(UqCheckpointError::InvalidEncoding("unknown checkpoint version"));
        }
        let count_bytes: [u8; 8] = bytes[41..HEADER_LEN]
            .try_into()
            .map_err(|_| UqCheckpointError::InvalidEncoding("missing observation count"))?;
        let count = usize::try_from(u64::from_le_bytes(count_bytes))
            .map_err(|_| UqCheckpointError::InvalidEncoding("observation count overflows usize"))?;
        if count > plan.budget_max_samples || bytes.len() != FIXED_LEN + count * 8 {
            return Err(UqCheckpointError::InvalidEncoding("observation count or trailing bytes"));
        }
        let status = match bytes[40] {
            0 => UqStatus::BudgetTruncated,
            1 => UqStatus::Cancelled,
            2 => UqStatus::Complete,
            _ => return Err(UqCheckpointError::InvalidEncoding("unknown or refused status")),
        };
        if (status == UqStatus::Complete) != (count == plan.budget_max_samples) {
            return Err(UqCheckpointError::InvalidEncoding("status disagrees with original budget"));
        }
        let payload_end = bytes.len() - 32;
        let checksum = hash_domain(CHECKSUM_DOMAIN, &bytes[..payload_end]);
        if checksum.0.as_slice() != &bytes[payload_end..] {
            return Err(UqCheckpointError::IntegrityMismatch);
        }
        if execution_identity(plan, model_identity).0.as_slice() != &bytes[8..40] {
            return Err(UqCheckpointError::IdentityMismatch);
        }
        let mut values = Vec::with_capacity(count);
        for chunk in bytes[HEADER_LEN..payload_end].chunks_exact(8) {
            let encoded: [u8; 8] = chunk
                .try_into()
                .map_err(|_| UqCheckpointError::InvalidEncoding("incomplete observation"))?;
            let value = f64::from_bits(u64::from_le_bytes(encoded));
            if !value.is_finite() {
                return Err(UqCheckpointError::InvalidEncoding("non-finite observation"));
            }
            values.push(value);
        }
        let execution = Self {
            plan: plan.clone(),
            factor,
            values,
            attempted: count,
            status,
            failure: None,
        };
        if execution.report().status == UqStatus::Refused {
            return Err(UqCheckpointError::InvalidEncoding("unrepresentable completed statistics"));
        }
        Ok(execution)
    }
}

fn text(bytes: &mut Vec<u8>, value: &str) {
    bytes.extend_from_slice(&(value.len() as u64).to_le_bytes());
    bytes.extend_from_slice(value.as_bytes());
}

fn numbers(bytes: &mut Vec<u8>, values: &[f64]) {
    for value in values {
        bytes.extend_from_slice(&value.to_bits().to_le_bytes());
    }
}

fn execution_identity(plan: &UqPlan, model_identity: ContentHash) -> ContentHash {
    // A length-framed, explicitly tagged encoding, not Debug/JSON formatting.
    // This domain also pins the Philox ordinal addressing and normal transform;
    // changing sampler semantics requires a new domain/version and checkpoint.
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&model_identity.0);
    text(&mut bytes, &plan.target_qoi);
    bytes.push(match plan.method {
        PropagationMethod::MonteCarlo => 0,
        PropagationMethod::QuasiMonteCarlo => 1,
        PropagationMethod::PolynomialChaos => 2,
        PropagationMethod::MultilevelMonteCarlo => 3,
        PropagationMethod::EpistemicBounding => 4,
    });
    bytes.extend_from_slice(&(plan.budget_max_samples as u64).to_le_bytes());
    bytes.extend_from_slice(&plan.seed.to_le_bytes());
    match plan.compliance_threshold {
        None => bytes.push(0),
        Some(value) => {
            bytes.push(1);
            numbers(&mut bytes, &[value]);
        }
    }
    bytes.extend_from_slice(&(plan.parameters.len() as u64).to_le_bytes());
    for parameter in &plan.parameters {
        text(&mut bytes, &parameter.name);
        text(&mut bytes, &parameter.unit);
        match parameter.kind {
            UncertaintyKind::AleatoryGaussian { mean, std_dev } => {
                bytes.push(0);
                numbers(&mut bytes, &[mean, std_dev]);
            }
            UncertaintyKind::AleatoryUniform { lo, hi } => {
                bytes.push(1);
                numbers(&mut bytes, &[lo, hi]);
            }
            UncertaintyKind::EpistemicInterval { lo, hi } => {
                bytes.push(2);
                numbers(&mut bytes, &[lo, hi]);
            }
            UncertaintyKind::StatisticalConfidence { estimate, half_width, confidence } => {
                bytes.push(3);
                numbers(&mut bytes, &[estimate, half_width, confidence]);
            }
            UncertaintyKind::Unstated => bytes.push(4),
        }
    }
    let matrix = match &plan.correlation {
        CorrelationModel::Unknown => { bytes.push(0); None }
        CorrelationModel::Independent => { bytes.push(1); None }
        CorrelationModel::Correlated { matrix } => { bytes.push(2); Some(matrix) }
        CorrelationModel::JointGaussian { matrix } => { bytes.push(3); Some(matrix) }
    };
    if let Some(matrix) = matrix {
        bytes.extend_from_slice(&(matrix.len() as u64).to_le_bytes());
        for row in matrix {
            bytes.extend_from_slice(&(row.len() as u64).to_le_bytes());
            numbers(&mut bytes, row);
        }
    }
    hash_domain("org.frankensim.uq.execution.v1", &bytes)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ParameterUncertainty, UqPropagator};

    fn plan() -> UqPlan {
        UqPlan::new("junction", PropagationMethod::MonteCarlo, 23)
            .with_parameter(ParameterUncertainty::uniform("ambient", 290.0, 310.0, "K"))
            .with_parameter(ParameterUncertainty::gaussian("power", 50.0, 2.0, "W"))
            .with_correlation(CorrelationModel::Independent)
            .with_compliance_threshold(345.0)
    }

    fn model() -> ContentHash { hash_domain("test-model", b"T = ambient + 0.8 * power") }

    fn qoi(x: &[f64]) -> f64 { x[0] + 0.8 * x[1] }

    fn checkpoint() -> Vec<u8> {
        let mut execution = UqExecution::new(&plan()).unwrap();
        execution.advance(7, || false, |x| Ok::<_, &str>(qoi(x)));
        execution.checkpoint(model()).unwrap()
    }

    fn reseal(bytes: &mut [u8]) {
        let end = bytes.len() - 32;
        let checksum = hash_domain(CHECKSUM_DOMAIN, &bytes[..end]);
        bytes[end..].copy_from_slice(&checksum.0);
    }

    #[test]
    fn durable_resume_evaluates_only_suffix_and_matches_uninterrupted_run() {
        let plan = plan();
        let mut restored = UqExecution::restore(&plan, model(), &checkpoint()).unwrap();
        assert_eq!(restored.evaluations_attempted(), 7);
        let mut calls = 0;
        let result = restored.advance(usize::MAX, || false, |x| {
            calls += 1;
            Ok::<_, &str>(qoi(x))
        });
        assert_eq!(calls, 16);
        assert_eq!(result, UqPropagator::run(&plan, qoi));
        let mut uninterrupted = UqExecution::new(&plan).unwrap();
        uninterrupted.advance(23, || false, |x| Ok::<_, &str>(qoi(x)));
        assert_eq!(restored.checkpoint(model()), uninterrupted.checkpoint(model()));
    }

    #[test]
    fn every_plan_field_and_model_are_bound_without_ambiguous_framing() {
        let base = plan();
        let bytes = checkpoint();
        let mut mutations = Vec::new();
        let mut p = base.clone(); p.seed += 1; mutations.push(p);
        let mut p = base.clone(); p.budget_max_samples += 1; mutations.push(p);
        let mut p = base.clone(); p.target_qoi.push('x'); mutations.push(p);
        let mut p = base.clone(); p.compliance_threshold = None; mutations.push(p);
        let mut p = base.clone(); p.parameters[0].name.push('x'); mutations.push(p);
        let mut p = base.clone(); p.parameters[0].unit.push('x'); mutations.push(p);
        let mut p = base.clone(); p.parameters.swap(0, 1); mutations.push(p);
        let mut p = base.clone(); p.parameters[0].kind = UncertaintyKind::AleatoryUniform { lo: 289.0, hi: 310.0 }; mutations.push(p);
        for changed in mutations {
            assert_eq!(UqExecution::restore(&changed, model(), &bytes).unwrap_err(), UqCheckpointError::IdentityMismatch);
        }
        assert_eq!(UqExecution::restore(&base, ContentHash([0; 32]), &bytes).unwrap_err(), UqCheckpointError::IdentityMismatch);
        let mut a = base.clone();
        a.parameters[0].name = "ab".into(); a.parameters[0].unit = "c".into();
        let mut b = a.clone();
        b.parameters[0].name = "a".into(); b.parameters[0].unit = "bc".into();
        assert_ne!(execution_identity(&a, model()), execution_identity(&b, model()));
    }

    #[test]
    fn truncated_extended_and_corrupted_payloads_are_rejected() {
        let bytes = checkpoint();
        for end in 0..bytes.len() {
            assert!(UqExecution::restore(&plan(), model(), &bytes[..end]).is_err());
        }
        let mut extended = bytes.clone(); extended.push(0);
        assert!(UqExecution::restore(&plan(), model(), &extended).is_err());
        for index in [0, 8, 40, 41, HEADER_LEN, bytes.len() - 1] {
            let mut corrupted = bytes.clone(); corrupted[index] ^= 0x80;
            assert!(UqExecution::restore(&plan(), model(), &corrupted).is_err());
        }
    }

    #[test]
    fn self_consistent_checksums_do_not_bypass_structural_admission() {
        let bytes = checkpoint();
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut bad = bytes.clone();
            bad[HEADER_LEN..HEADER_LEN + 8].copy_from_slice(&value.to_bits().to_le_bytes());
            reseal(&mut bad);
            assert!(UqExecution::restore(&plan(), model(), &bad).is_err());
        }
        let mut wrong_status = bytes.clone(); wrong_status[40] = 2; reseal(&mut wrong_status);
        assert!(UqExecution::restore(&plan(), model(), &wrong_status).is_err());
        let mut huge_count = bytes; huge_count[41..HEADER_LEN].copy_from_slice(&u64::MAX.to_le_bytes()); reseal(&mut huge_count);
        assert!(UqExecution::restore(&plan(), model(), &huge_count).is_err());
    }

    #[test]
    fn cancelled_empty_and_completed_states_round_trip_without_evaluation() {
        let plan = plan();
        let mut execution = UqExecution::new(&plan).unwrap();
        for _ in 0..2 {
            let bytes = execution.checkpoint(model()).unwrap();
            let restored = UqExecution::restore(&plan, model(), &bytes).unwrap();
            assert_eq!(restored.report(), execution.report());
            execution.advance(1, || true, |x| Ok::<_, &str>(qoi(x)));
        }
        execution.advance(23, || false, |x| Ok::<_, &str>(qoi(x)));
        let bytes = execution.checkpoint(model()).unwrap();
        let mut restored = UqExecution::restore(&plan, model(), &bytes).unwrap();
        assert_eq!(restored.advance(23, || panic!("complete"), |_| -> Result<f64, &str> { panic!("complete") }), execution.report());
    }

    #[test]
    fn signed_zero_observation_bits_survive_round_trip() {
        let mut execution = UqExecution::new(&plan()).unwrap();
        execution.advance(1, || false, |_| Ok::<_, &str>(-0.0));
        let bytes = execution.checkpoint(model()).unwrap();
        let restored = UqExecution::restore(&plan(), model(), &bytes).unwrap();
        assert_eq!(restored.observations()[0].to_bits(), (-0.0_f64).to_bits());
    }

    #[test]
    fn failures_cannot_be_restarted_from_their_successful_prefix() {
        let mut execution = UqExecution::new(&plan()).unwrap();
        execution.advance(7, || false, |x| Ok::<_, &str>(qoi(x)));
        execution.advance(1, || false, |_| Err::<f64, _>("solver failed"));
        assert_eq!(execution.checkpoint(model()), Err(UqCheckpointError::NotResumable));
        assert_eq!(execution.evaluations_attempted(), 8);
    }

    #[test]
    fn dependence_and_exact_float_encodings_participate_in_identity() {
        let mut a = plan();
        a.parameters = vec![
            ParameterUncertainty::gaussian("a", 0.0, 1.0, "1"),
            ParameterUncertainty::gaussian("b", 0.0, 1.0, "1"),
        ];
        a.correlation = CorrelationModel::JointGaussian { matrix: vec![vec![1.0, 0.5], vec![0.5, 1.0]] };
        let bytes = UqExecution::new(&a).unwrap().checkpoint(model()).unwrap();
        let mut b = a.clone();
        b.correlation = CorrelationModel::JointGaussian { matrix: vec![vec![1.0, -0.5], vec![-0.5, 1.0]] };
        assert_eq!(UqExecution::restore(&b, model(), &bytes).unwrap_err(), UqCheckpointError::IdentityMismatch);
        let mut c = a;
        c.parameters[0].kind = UncertaintyKind::AleatoryGaussian { mean: -0.0, std_dev: 1.0 };
        assert_eq!(UqExecution::restore(&c, model(), &bytes).unwrap_err(), UqCheckpointError::IdentityMismatch);
    }
}
