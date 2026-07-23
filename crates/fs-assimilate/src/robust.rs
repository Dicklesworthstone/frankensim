//! Audited robust-observation admission around the scalar Kalman core.
//!
//! The implementation deliberately reuses `assimilate_all`: an exactly
//! diagonal observation covariance takes the existing path unchanged, while a
//! correlated covariance is whitened into ordinary scalar observations. This
//! keeps one owner for Joseph covariance updates and posterior PSD admission.

use fs_blake3::DomainHasher;
use fs_exec::Cx;

use super::{
    AssimError, Belief, Observation, assimilate_all, canonical_f64_bits, canonicalize_zero,
    validate_leaf_identity,
};

const AUDIT_ID_DOMAIN: &str = "org.frankensim.fs-assimilate.robust-observation-audit.v1";
const AUDIT_ID_PREFIX: &str = "robust-observation-audit:v1:";

/// Maximum number of numerically available observations in the dense
/// correlated path.
///
/// Cholesky factorization and whitening are cubic in this dimension. Larger
/// campaigns require a sparse/shared-factor implementation with its own
/// resource contract.
pub const MAX_CORRELATED_OBSERVATIONS: usize = 256;

/// Which side of a threshold contains an unobserved censored value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CensorDirection {
    /// The latent value is at or below the supplied threshold.
    AtOrBelow,
    /// The latent value is at or above the supplied threshold.
    AtOrAbove,
}

/// The pathology attached to a robust observation record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PathologyKind {
    /// A numerically available, ordinary scalar observation.
    Available,
    /// No reading was produced.
    Missing,
    /// Only a one-sided threshold is known.
    Censored,
    /// The instrument clipped a value to its declared range.
    Saturated,
    /// The reading has a declared lag that this v1 core cannot state-augment.
    Delayed,
}

#[derive(Debug, Clone, PartialEq)]
enum ObservationState {
    Available(Observation),
    Missing {
        instrument: String,
    },
    Censored {
        observation: Observation,
        threshold: f64,
        direction: CensorDirection,
    },
    Saturated {
        observation: Observation,
        lower: f64,
        upper: f64,
    },
    Delayed {
        observation: Observation,
        time_constant: f64,
        age: f64,
    },
}

/// One checked record in a robust observation batch.
///
/// Fields and the internal state enum are private so callers cannot bypass
/// identity, bound, or finite-value validation.
#[derive(Debug, Clone, PartialEq)]
pub struct RobustObservation {
    state: ObservationState,
}

impl RobustObservation {
    /// Wrap an ordinary checked scalar observation.
    #[must_use]
    pub fn available(observation: Observation) -> Self {
        Self {
            state: ObservationState::Available(observation),
        }
    }

    /// Declare that an instrument produced no reading.
    ///
    /// Missing records are excluded and reduce the effective observation
    /// degrees of freedom; they do not synthesize a value or variance.
    pub fn missing(instrument: impl Into<String>) -> Result<Self, AssimError> {
        let instrument = instrument.into();
        validate_leaf_identity("instrument", &instrument)?;
        Ok(Self {
            state: ObservationState::Missing { instrument },
        })
    }

    /// Declare a one-sided censored reading.
    ///
    /// The v1 robust core records and refuses this pathology rather than
    /// treating the threshold as an ordinary Gaussian value.
    pub fn censored(
        observation: Observation,
        threshold: f64,
        direction: CensorDirection,
    ) -> Result<Self, AssimError> {
        if !threshold.is_finite() {
            return Err(AssimError::NonFinitePathologyParameter {
                parameter: "censor threshold",
            });
        }
        Ok(Self {
            state: ObservationState::Censored {
                observation,
                threshold: canonicalize_zero(threshold),
                direction,
            },
        })
    }

    /// Declare a range-saturated reading.
    ///
    /// The clipped observation value must equal one declared endpoint. The v1
    /// core audits and refuses the record instead of assimilating the clip.
    pub fn saturated(observation: Observation, lower: f64, upper: f64) -> Result<Self, AssimError> {
        if !lower.is_finite() || !upper.is_finite() {
            return Err(AssimError::NonFinitePathologyParameter {
                parameter: "saturation range",
            });
        }
        if lower > upper {
            return Err(AssimError::InvalidPathologyParameter {
                parameter: "saturation range",
            });
        }
        let value_bits = canonical_f64_bits(observation.value());
        if value_bits != canonical_f64_bits(lower) && value_bits != canonical_f64_bits(upper) {
            return Err(AssimError::InvalidPathologyParameter {
                parameter: "saturated value endpoint",
            });
        }
        Ok(Self {
            state: ObservationState::Saturated {
                observation,
                lower: canonicalize_zero(lower),
                upper: canonicalize_zero(upper),
            },
        })
    }

    /// Declare a delayed reading with a first-order sensor time constant.
    ///
    /// Both values are retained in the audit identity. The current core
    /// refuses publication because it does not own a temporal state history
    /// from which to perform the required lag augmentation.
    pub fn delayed(
        observation: Observation,
        time_constant: f64,
        age: f64,
    ) -> Result<Self, AssimError> {
        if !time_constant.is_finite() {
            return Err(AssimError::NonFinitePathologyParameter {
                parameter: "delay time constant",
            });
        }
        if !age.is_finite() {
            return Err(AssimError::NonFinitePathologyParameter {
                parameter: "delay age",
            });
        }
        if time_constant <= 0.0 {
            return Err(AssimError::InvalidPathologyParameter {
                parameter: "delay time constant",
            });
        }
        if age < 0.0 {
            return Err(AssimError::InvalidPathologyParameter {
                parameter: "delay age",
            });
        }
        Ok(Self {
            state: ObservationState::Delayed {
                observation,
                time_constant,
                age: canonicalize_zero(age),
            },
        })
    }

    /// The calibrated instrument identity.
    #[must_use]
    pub fn instrument(&self) -> &str {
        match &self.state {
            ObservationState::Available(observation)
            | ObservationState::Censored { observation, .. }
            | ObservationState::Saturated { observation, .. }
            | ObservationState::Delayed { observation, .. } => observation.instrument(),
            ObservationState::Missing { instrument } => instrument,
        }
    }

    /// The declared pathology class.
    #[must_use]
    pub fn pathology(&self) -> PathologyKind {
        match self.state {
            ObservationState::Available(_) => PathologyKind::Available,
            ObservationState::Missing { .. } => PathologyKind::Missing,
            ObservationState::Censored { .. } => PathologyKind::Censored,
            ObservationState::Saturated { .. } => PathologyKind::Saturated,
            ObservationState::Delayed { .. } => PathologyKind::Delayed,
        }
    }
}

/// Checked dense robust-observation batch.
///
/// Covariance rows correspond, in declaration order, only to records whose
/// pathology is [`PathologyKind::Available`]. Every diagonal must exactly equal
/// the corresponding [`Observation::noise_var`], leaving one noise authority.
#[derive(Debug, Clone, PartialEq)]
pub struct ObservationBatch {
    records: Vec<RobustObservation>,
    covariance: Vec<Vec<f64>>,
    cholesky: Vec<Vec<f64>>,
    diagonal: bool,
}

impl ObservationBatch {
    /// Construct and numerically gate a bounded dense batch.
    ///
    /// The covariance must be finite, exactly symmetric, and strictly positive
    /// definite under the deterministic binary64 Cholesky computation. A
    /// singular PSD or numerically unresolved matrix is refused without being
    /// described as indefinite.
    pub fn new(
        records: Vec<RobustObservation>,
        covariance: Vec<Vec<f64>>,
        cx: &Cx<'_>,
    ) -> Result<Self, AssimError> {
        if records.is_empty() {
            return Err(AssimError::EmptyObservations);
        }
        if records.len() > super::MAX_DENSE_OBSERVATIONS {
            return Err(AssimError::ObservationCountLimit {
                count: records.len(),
                max: super::MAX_DENSE_OBSERVATIONS,
            });
        }
        cx.checkpoint()
            .map_err(|_| robust_cancelled("batch preflight", 0, 1))?;

        for (index, record) in records.iter().enumerate() {
            validate_leaf_identity("instrument", record.instrument())?;
            if records[..index]
                .iter()
                .any(|prior| prior.instrument() == record.instrument())
            {
                return Err(AssimError::DuplicateObservationInstrument {
                    instrument: record.instrument().to_owned(),
                });
            }
        }

        let available: Vec<&Observation> = records
            .iter()
            .filter_map(|record| match &record.state {
                ObservationState::Available(observation) => Some(observation),
                ObservationState::Missing { .. }
                | ObservationState::Censored { .. }
                | ObservationState::Saturated { .. }
                | ObservationState::Delayed { .. } => None,
            })
            .collect();
        if available.len() > MAX_CORRELATED_OBSERVATIONS {
            return Err(AssimError::ObservationCountLimit {
                count: available.len(),
                max: MAX_CORRELATED_OBSERVATIONS,
            });
        }
        validate_covariance_shape(&covariance, available.len())?;
        let mut diagonal = true;
        for (row_index, row) in covariance.iter().enumerate() {
            for (column_index, value) in row.iter().enumerate() {
                if !value.is_finite() {
                    return Err(AssimError::NonFiniteObservationCovariance {
                        row: row_index,
                        column: column_index,
                    });
                }
                if row_index != column_index && canonical_f64_bits(*value) != 0 {
                    diagonal = false;
                }
            }
            if canonical_f64_bits(row[row_index])
                != canonical_f64_bits(available[row_index].noise_var())
            {
                return Err(AssimError::ObservationCovarianceNoiseMismatch { index: row_index });
            }
        }
        for row in 0..covariance.len() {
            for (column, column_row) in covariance.iter().enumerate().skip(row + 1) {
                if canonical_f64_bits(covariance[row][column])
                    != canonical_f64_bits(column_row[row])
                {
                    return Err(AssimError::NonSymmetricObservationCovariance { row, column });
                }
            }
        }
        if !covariance.is_empty() {
            match Belief::new(vec![0.0; covariance.len()], covariance.clone(), cx) {
                Ok(_) => {}
                Err(AssimError::CovarianceNotPositiveSemidefinite) => {
                    return Err(AssimError::ObservationCovarianceNotPositiveSemidefinite);
                }
                Err(AssimError::CovarianceCertificationUnresolved) => {
                    return Err(AssimError::ObservationCovarianceCertificationUnresolved);
                }
                Err(error) => return Err(error),
            }
        }
        let cholesky = strict_cholesky(&covariance, cx)?;
        Ok(Self {
            records,
            covariance,
            cholesky,
            diagonal,
        })
    }

    /// Read-only record list in covariance-row declaration order.
    #[must_use]
    pub fn records(&self) -> &[RobustObservation] {
        &self.records
    }

    /// Read-only full covariance over available records.
    #[must_use]
    pub fn covariance(&self) -> &[Vec<f64>] {
        &self.covariance
    }

    /// Whether the observation covariance is exactly diagonal.
    #[must_use]
    pub fn is_diagonal(&self) -> bool {
        self.diagonal
    }
}

/// Per-record action retained in the robust audit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ObservationDisposition {
    /// The ordinary observation contributed to the published posterior.
    Accepted,
    /// The record was absent and reduced the effective degrees of freedom.
    ExcludedMissing,
    /// The record was otherwise usable but the batch contained a pathology
    /// that requires fail-closed refusal.
    WithheldByBatchRefusal,
    /// A censored likelihood is not implemented in v1.
    RefusedCensored,
    /// A clipped likelihood is not implemented in v1.
    RefusedSaturated,
    /// Temporal state augmentation is not implemented in v1.
    RefusedDelayed,
}

/// Terminal outcome of robust batch admission.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BatchOutcome {
    /// A posterior was published.
    Updated,
    /// No reading was numerically available.
    RefusedNoUsableObservations,
    /// At least one pathology requires a likelihood/state model not owned by
    /// this v1 core.
    RefusedPathology,
}

/// One immutable audit entry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ObservationAudit {
    entries: Vec<(String, ObservationDisposition)>,
    outcome: BatchOutcome,
    effective_dof: usize,
    receipt_id: String,
}

impl ObservationAudit {
    /// Instrument/disposition pairs in caller declaration order.
    #[must_use]
    pub fn entries(&self) -> &[(String, ObservationDisposition)] {
        &self.entries
    }

    /// Terminal batch outcome.
    #[must_use]
    pub fn outcome(&self) -> BatchOutcome {
        self.outcome
    }

    /// Count of scalar readings that actually contributed to the posterior.
    #[must_use]
    pub fn effective_dof(&self) -> usize {
        self.effective_dof
    }

    /// Domain-separated receipt binding the prior, batch, disposition, and
    /// published posterior (if any).
    #[must_use]
    pub fn receipt_id(&self) -> &str {
        &self.receipt_id
    }
}

/// Audited robust assimilation result.
///
/// A refused outcome carries no posterior, ensuring that callers cannot
/// accidentally consume a partial update while ignoring the audit.
#[derive(Debug, Clone, PartialEq)]
pub struct RobustAssimilation {
    posterior: Option<Belief>,
    audit: ObservationAudit,
}

impl RobustAssimilation {
    /// Published posterior, present only when `audit().outcome() == Updated`.
    #[must_use]
    pub fn posterior(&self) -> Option<&Belief> {
        self.posterior.as_ref()
    }

    /// Complete immutable disposition audit.
    #[must_use]
    pub fn audit(&self) -> &ObservationAudit {
        &self.audit
    }
}

/// Apply a checked robust batch.
///
/// Missing readings are excluded. Censored, saturated, and delayed readings
/// return a fail-closed audited refusal with no posterior. An exactly diagonal
/// covariance delegates the original observations to [`assimilate_all`]
/// unchanged. A correlated covariance is deterministically Cholesky-whitened
/// and then uses that same scalar Joseph path.
pub fn assimilate_observation_batch(
    prior: &Belief,
    batch: &ObservationBatch,
    cx: &Cx<'_>,
) -> Result<RobustAssimilation, AssimError> {
    cx.checkpoint()
        .map_err(|_| robust_cancelled("robust admission", 0, 1))?;
    let has_refused_pathology = batch.records.iter().any(|record| {
        matches!(
            record.state,
            ObservationState::Censored { .. }
                | ObservationState::Saturated { .. }
                | ObservationState::Delayed { .. }
        )
    });
    if has_refused_pathology {
        return build_result(prior, batch, None, BatchOutcome::RefusedPathology, cx);
    }

    let available: Vec<&Observation> = batch
        .records
        .iter()
        .filter_map(|record| match &record.state {
            ObservationState::Available(observation) => Some(observation),
            ObservationState::Missing { .. }
            | ObservationState::Censored { .. }
            | ObservationState::Saturated { .. }
            | ObservationState::Delayed { .. } => None,
        })
        .collect();
    if available.is_empty() {
        return build_result(
            prior,
            batch,
            None,
            BatchOutcome::RefusedNoUsableObservations,
            cx,
        );
    }

    let posterior = if batch.diagonal {
        let observations: Vec<Observation> = available.into_iter().cloned().collect();
        assimilate_all(prior, &observations, cx)?
    } else {
        let whitened = whiten_observations(&available, &batch.cholesky, cx)?;
        assimilate_all(prior, &whitened, cx)?
    };
    let result = build_result(prior, batch, Some(posterior), BatchOutcome::Updated, cx)?;
    cx.checkpoint()
        .map_err(|_| robust_cancelled("robust publication", 1, 2))?;
    Ok(result)
}

fn validate_covariance_shape(
    covariance: &[Vec<f64>],
    observations: usize,
) -> Result<(), AssimError> {
    if covariance.len() != observations {
        return Err(AssimError::ObservationCovarianceDimensionMismatch {
            observations,
            rows: covariance.len(),
        });
    }
    for (row_index, row) in covariance.iter().enumerate() {
        if row.len() != observations {
            return Err(AssimError::ObservationCovarianceRowDimensionMismatch {
                row: row_index,
                expected: observations,
                actual: row.len(),
            });
        }
    }
    Ok(())
}

fn strict_cholesky(matrix: &[Vec<f64>], cx: &Cx<'_>) -> Result<Vec<Vec<f64>>, AssimError> {
    let mut lower: Vec<Vec<f64>> = vec![vec![0.0; matrix.len()]; matrix.len()];
    for row in 0..matrix.len() {
        cx.checkpoint().map_err(|_| {
            robust_cancelled(
                "observation covariance Cholesky",
                row as u128,
                matrix.len() as u128,
            )
        })?;
        for column in 0..=row {
            let mut residual = matrix[row][column];
            for (&row_factor, &column_factor) in
                lower[row][..column].iter().zip(&lower[column][..column])
            {
                residual = (-row_factor).mul_add(column_factor, residual);
            }
            if !residual.is_finite() {
                return Err(AssimError::NonFiniteComputation {
                    stage: "observation covariance Cholesky",
                });
            }
            if row == column {
                if residual <= 0.0 {
                    return Err(AssimError::ObservationCovarianceNotPositiveDefinite {
                        pivot: row,
                    });
                }
                lower[row][column] = residual.sqrt();
            } else {
                lower[row][column] = residual / lower[column][column];
            }
            if !lower[row][column].is_finite() {
                return Err(AssimError::NonFiniteComputation {
                    stage: "observation covariance Cholesky",
                });
            }
        }
    }
    Ok(lower)
}

fn whiten_observations(
    observations: &[&Observation],
    lower: &[Vec<f64>],
    cx: &Cx<'_>,
) -> Result<Vec<Observation>, AssimError> {
    let state_dim = observations[0].operator().len();
    let mut values = vec![0.0; observations.len()];
    let mut operators = vec![vec![0.0; state_dim]; observations.len()];
    for row in 0..observations.len() {
        cx.checkpoint().map_err(|_| {
            robust_cancelled(
                "correlated observation whitening",
                row as u128,
                observations.len() as u128,
            )
        })?;
        let mut value = observations[row].value();
        for prior in 0..row {
            value = (-lower[row][prior]).mul_add(values[prior], value);
        }
        values[row] = value / lower[row][row];
        if !values[row].is_finite() {
            return Err(AssimError::NonFiniteComputation {
                stage: "whitened observation value",
            });
        }
        let (prior_operators, current_and_later) = operators.split_at_mut(row);
        let current_operator = &mut current_and_later[0];
        for (component, output) in current_operator.iter_mut().enumerate() {
            let mut coefficient = observations[row].operator()[component];
            for prior in 0..row {
                coefficient =
                    (-lower[row][prior]).mul_add(prior_operators[prior][component], coefficient);
            }
            *output = coefficient / lower[row][row];
            if !output.is_finite() {
                return Err(AssimError::NonFiniteComputation {
                    stage: "whitened observation operator",
                });
            }
        }
    }

    operators
        .into_iter()
        .zip(values)
        .enumerate()
        .map(|(index, (operator, value))| {
            Observation::new(operator, value, 1.0, format!("robust-whitened-{index:03}"))
        })
        .collect()
}

fn build_result(
    prior: &Belief,
    batch: &ObservationBatch,
    posterior: Option<Belief>,
    outcome: BatchOutcome,
    cx: &Cx<'_>,
) -> Result<RobustAssimilation, AssimError> {
    let entries = batch
        .records
        .iter()
        .map(|record| {
            let disposition = match record.state {
                ObservationState::Available(_) if outcome == BatchOutcome::Updated => {
                    ObservationDisposition::Accepted
                }
                ObservationState::Available(_) => ObservationDisposition::WithheldByBatchRefusal,
                ObservationState::Missing { .. } => ObservationDisposition::ExcludedMissing,
                ObservationState::Censored { .. } => ObservationDisposition::RefusedCensored,
                ObservationState::Saturated { .. } => ObservationDisposition::RefusedSaturated,
                ObservationState::Delayed { .. } => ObservationDisposition::RefusedDelayed,
            };
            (record.instrument().to_owned(), disposition)
        })
        .collect::<Vec<_>>();
    let effective_dof = entries
        .iter()
        .filter(|(_, disposition)| *disposition == ObservationDisposition::Accepted)
        .count();
    let receipt_id = audit_receipt(prior, batch, posterior.as_ref(), outcome, &entries, cx)?;
    Ok(RobustAssimilation {
        posterior,
        audit: ObservationAudit {
            entries,
            outcome,
            effective_dof,
            receipt_id,
        },
    })
}

fn audit_receipt(
    prior: &Belief,
    batch: &ObservationBatch,
    posterior: Option<&Belief>,
    outcome: BatchOutcome,
    entries: &[(String, ObservationDisposition)],
    cx: &Cx<'_>,
) -> Result<String, AssimError> {
    let mut hasher = DomainHasher::new(AUDIT_ID_DOMAIN);
    hash_belief(&mut hasher, prior, cx)?;
    hash_u64(&mut hasher, batch.records.len() as u64);
    for (index, record) in batch.records.iter().enumerate() {
        cx.checkpoint().map_err(|_| {
            robust_cancelled(
                "robust audit records",
                index as u128,
                batch.records.len() as u128,
            )
        })?;
        hash_record(&mut hasher, record);
    }
    hash_matrix(
        &mut hasher,
        &batch.covariance,
        cx,
        "robust audit covariance",
    )?;
    hash_u8(&mut hasher, outcome_tag(outcome));
    for (instrument, disposition) in entries {
        hash_bytes(&mut hasher, instrument.as_bytes());
        hash_u8(&mut hasher, disposition_tag(*disposition));
    }
    match posterior {
        Some(belief) => {
            hash_u8(&mut hasher, 1);
            hash_belief(&mut hasher, belief, cx)?;
        }
        None => hash_u8(&mut hasher, 0),
    }
    cx.checkpoint()
        .map_err(|_| robust_cancelled("robust audit finalize", 1, 2))?;
    Ok(format!("{AUDIT_ID_PREFIX}{}", hasher.finalize()))
}

fn hash_record(hasher: &mut DomainHasher, record: &RobustObservation) {
    hash_bytes(hasher, record.instrument().as_bytes());
    match &record.state {
        ObservationState::Available(observation) => {
            hash_u8(hasher, 0);
            hash_observation(hasher, observation);
        }
        ObservationState::Missing { .. } => hash_u8(hasher, 1),
        ObservationState::Censored {
            observation,
            threshold,
            direction,
        } => {
            hash_u8(hasher, 2);
            hash_observation(hasher, observation);
            hash_f64(hasher, *threshold);
            hash_u8(
                hasher,
                match direction {
                    CensorDirection::AtOrBelow => 0,
                    CensorDirection::AtOrAbove => 1,
                },
            );
        }
        ObservationState::Saturated {
            observation,
            lower,
            upper,
        } => {
            hash_u8(hasher, 3);
            hash_observation(hasher, observation);
            hash_f64(hasher, *lower);
            hash_f64(hasher, *upper);
        }
        ObservationState::Delayed {
            observation,
            time_constant,
            age,
        } => {
            hash_u8(hasher, 4);
            hash_observation(hasher, observation);
            hash_f64(hasher, *time_constant);
            hash_f64(hasher, *age);
        }
    }
}

fn hash_observation(hasher: &mut DomainHasher, observation: &Observation) {
    hash_u64(hasher, observation.operator().len() as u64);
    for coefficient in observation.operator() {
        hash_f64(hasher, *coefficient);
    }
    hash_f64(hasher, observation.value());
    hash_f64(hasher, observation.noise_var());
    hash_bytes(hasher, observation.instrument().as_bytes());
}

fn hash_belief(hasher: &mut DomainHasher, belief: &Belief, cx: &Cx<'_>) -> Result<(), AssimError> {
    hash_u64(hasher, belief.dim() as u64);
    for value in belief.mean() {
        hash_f64(hasher, *value);
    }
    hash_matrix(hasher, belief.covariance(), cx, "robust audit belief")
}

fn hash_matrix(
    hasher: &mut DomainHasher,
    matrix: &[Vec<f64>],
    cx: &Cx<'_>,
    phase: &'static str,
) -> Result<(), AssimError> {
    hash_u64(hasher, matrix.len() as u64);
    for (row_index, row) in matrix.iter().enumerate() {
        cx.checkpoint()
            .map_err(|_| robust_cancelled(phase, row_index as u128, matrix.len() as u128))?;
        hash_u64(hasher, row.len() as u64);
        for value in row {
            hash_f64(hasher, *value);
        }
    }
    Ok(())
}

fn hash_bytes(hasher: &mut DomainHasher, bytes: &[u8]) {
    hash_u64(hasher, bytes.len() as u64);
    hasher.update(bytes);
}

fn hash_f64(hasher: &mut DomainHasher, value: f64) {
    hasher.update(&canonical_f64_bits(value).to_le_bytes());
}

fn hash_u64(hasher: &mut DomainHasher, value: u64) {
    hasher.update(&value.to_le_bytes());
}

fn hash_u8(hasher: &mut DomainHasher, value: u8) {
    hasher.update(&[value]);
}

fn outcome_tag(outcome: BatchOutcome) -> u8 {
    match outcome {
        BatchOutcome::Updated => 0,
        BatchOutcome::RefusedNoUsableObservations => 1,
        BatchOutcome::RefusedPathology => 2,
    }
}

fn disposition_tag(disposition: ObservationDisposition) -> u8 {
    match disposition {
        ObservationDisposition::Accepted => 0,
        ObservationDisposition::ExcludedMissing => 1,
        ObservationDisposition::WithheldByBatchRefusal => 2,
        ObservationDisposition::RefusedCensored => 3,
        ObservationDisposition::RefusedSaturated => 4,
        ObservationDisposition::RefusedDelayed => 5,
    }
}

fn robust_cancelled(phase: &'static str, completed: u128, planned: u128) -> AssimError {
    AssimError::Cancelled {
        phase,
        completed,
        planned,
    }
}
