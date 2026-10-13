//! Shared multivariate calibration/registration errors with signed loadings.
//!
//! The caller supplies a zero-mean nuisance covariance `C` and each reading's
//! nuisance Jacobian `J_i`. The observation covariance is derived as
//! `R = diag(independent_noise) + J C J^T`. In particular, a registration RMS
//! is not a covariance, and parameter cross-correlations are not discarded.
//! The model is exact for additive linear Gaussian nuisances independent of
//! the state and instrument noise; a nonlinear sensor's supplied Jacobian is
//! only its local linearization. Nothing here validates a calibration model.

use super::RowId;
use crate::factor::FactorBelief;
use crate::robust::{
    MAX_CORRELATED_OBSERVATIONS, ObservationBatch, RobustAssimilation, RobustObservation,
    assimilate_observation_batch,
};
use crate::{
    AssimError, Belief, Cx, Observation, WorkPlan, WorkProgress, canonical_f64_bits,
    canonicalize_zero, validate_leaf_identity,
};
use fs_blake3::DomainHasher;

/// Bound on the dense nuisance model, separate from the physical state size.
pub const MAX_SHARED_PARAMETERS: usize = 64;
const MODEL_DOMAIN: &str = "org.frankensim.fs-assimilate.shared-noise-model.v1";
const BATCH_DOMAIN: &str = "org.frankensim.fs-assimilate.jacobian-batch.v1";
const ROW_DOMAIN: &str = "org.frankensim.fs-assimilate.jacobian-row.v1";
const RESULT_DOMAIN: &str = "org.frankensim.fs-assimilate.jacobian-result.v1";

/// One shared latent Gaussian model. Parameter order fixes Jacobian columns.
/// Singular PSD nuisance covariance is allowed; independent instrument noise
/// still makes the resulting observation covariance strictly positive definite
/// in exact arithmetic. Numerical admission can refuse an unresolved matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct SharedNoiseModel {
    identity: String,
    parameters: Vec<String>,
    covariance: Vec<Vec<f64>>,
    factor: FactorBelief,
    content_id: String,
}

impl SharedNoiseModel {
    /// Validate a named, full nuisance covariance; no diagonalization, jitter,
    /// inferred units, or substitution of a scalar registration error occurs.
    ///
    /// # Errors
    /// Refuses invalid/duplicate names, unsupported shapes, non-PSD or
    /// numerically unresolved covariance, cancellation, and budget exhaustion.
    pub fn new(
        identity: impl Into<String>,
        parameters: Vec<String>,
        covariance: Vec<Vec<f64>>,
        cx: &Cx<'_>,
    ) -> Result<Self, AssimError> {
        let identity = identity.into();
        validate_leaf_identity("shared noise model", &identity)?;
        let k = parameters.len();
        if k == 0 || k > MAX_SHARED_PARAMETERS {
            return Err(declaration(&identity, "shared parameter count must be in 1..=64"));
        }
        for (index, parameter) in parameters.iter().enumerate() {
            validate_leaf_identity("shared noise parameter", parameter)?;
            if parameters[..index].contains(parameter) {
                return Err(declaration(parameter, "duplicate shared parameter name"));
            }
        }
        let checked = Belief::new(vec![0.0; k], covariance, cx)?;
        let factor = FactorBelief::from_belief(&checked, cx)?;
        let ku = k as u128;
        // Every admitted identity is bounded to 256 bytes by leaf validation.
        let plan = WorkPlan::checked(0, 0, 0, 0, 0, 0, 8 * ku * ku + 264 * (ku + 1) + 16)?;
        let mut progress = WorkProgress::new(cx, plan)?;
        progress.checkpoint("shared model identity")?;
        let mut hash = DomainHasher::new(MODEL_DOMAIN);
        field(&mut hash, identity.as_bytes(), &mut progress)?;
        field(&mut hash, &(k as u64).to_le_bytes(), &mut progress)?;
        for parameter in &parameters {
            field(&mut hash, parameter.as_bytes(), &mut progress)?;
        }
        for row in checked.covariance() {
            for value in row { number(&mut hash, *value, &mut progress)?; }
        }
        progress.checkpoint("shared model publication")?;
        Ok(Self {
            identity,
            parameters,
            covariance: checked.covariance().to_vec(),
            factor,
            content_id: format!("shared-noise-model:v1:{}", hash.finalize()),
        })
    }

    /// Caller-supplied calibration/registration model identity.
    #[must_use]
    pub fn identity(&self) -> &str { &self.identity }

    /// Named nuisance parameters in Jacobian-column order.
    #[must_use]
    pub fn parameters(&self) -> &[String] { &self.parameters }

    /// The full admitted nuisance covariance, not its diagonal alone.
    #[must_use]
    pub fn covariance(&self) -> &[Vec<f64>] { &self.covariance }

    /// Content identity binds the name, ordered parameters, and covariance.
    #[must_use]
    pub fn content_id(&self) -> &str { &self.content_id }
}

/// An observation with INDEPENDENT instrument noise and a signed nuisance
/// Jacobian bound to a particular model. Do not fold the shared variance into
/// `observation.noise_var()` as well: that would count the same error twice.
#[derive(Debug, Clone, PartialEq)]
pub struct JacobianRecord {
    observation: Observation,
    row: RowId,
    jacobian: Vec<f64>,
    model_id: String,
}

impl JacobianRecord {
    /// Bind a finite signed Jacobian to the model's ordered parameters.
    /// A zero Jacobian is valid: that reading has no shared-error sensitivity.
    ///
    /// # Errors
    /// Refuses a non-finite Jacobian or a mismatched parameter count.
    pub fn new(
        observation: Observation,
        row: RowId,
        model: &SharedNoiseModel,
        mut jacobian: Vec<f64>,
    ) -> Result<Self, AssimError> {
        if jacobian.len() != model.parameters.len() {
            return Err(declaration(row.as_str(), "Jacobian dimension differs from shared model"));
        }
        for value in &mut jacobian {
            if !value.is_finite() {
                return Err(declaration(row.as_str(), "non-finite shared noise Jacobian"));
            }
            *value = canonicalize_zero(*value);
        }
        Ok(Self { observation, row, jacobian, model_id: model.content_id.clone() })
    }

    /// The original observation carrying independent, not total, noise.
    #[must_use]
    pub const fn observation(&self) -> &Observation { &self.observation }

    /// Stable dataset-row identity.
    #[must_use]
    pub const fn row(&self) -> &RowId { &self.row }

    /// Signed sensitivities in the bound model's parameter order.
    #[must_use]
    pub fn jacobian(&self) -> &[f64] { &self.jacobian }
}

/// A row-unique observation batch with covariance derived from the complete
/// nuisance model. Rows retain declaration order. Identity tracking is local
/// to this batch; this API does not deduplicate separate assimilation calls.
#[derive(Debug, Clone, PartialEq)]
pub struct JacobianBatch {
    model: SharedNoiseModel,
    records: Vec<JacobianRecord>,
    numerical_batch: ObservationBatch,
    content_id: String,
}

impl JacobianBatch {
    /// Form `T = J U sqrt(D)` and `R = diag(noise) + T T^T`, retaining signed
    /// cross-covariance without a dense inverse. Construction is O(m k²+m² k),
    /// followed by the existing strict observation-covariance admission.
    ///
    /// # Errors
    /// Refuses duplicate rows, mixed models, mismatched state dimensions,
    /// oversized batches, numerical overflow, budget exhaustion, cancellation,
    /// or unresolved/non-positive-definite observation covariance. No partial
    /// batch is published. Projection has one cumulative accountant; nuisance
    /// validation and robust admission remain separate checked stage scopes.
    pub fn new(
        model: SharedNoiseModel,
        records: Vec<JacobianRecord>,
        cx: &Cx<'_>,
    ) -> Result<Self, AssimError> {
        let m = records.len();
        if m == 0 { return Err(AssimError::EmptyObservations); }
        if m > MAX_CORRELATED_OBSERVATIONS {
            return Err(AssimError::ObservationCountLimit { count: m, max: MAX_CORRELATED_OBSERVATIONS });
        }
        let k = model.parameters.len();
        let n = records[0].observation.operator().len();
        let (mu, ku, nu) = (m as u128, k as u128, n as u128);
        // Shapes already bounded by checked constructors (m,n<=256,k<=64).
        let work = 4 * mu * ku * ku + 4 * mu * mu * ku + 8 * mu * (ku + nu + mu + 1)
            + 256 * mu * mu;
        let hashing = 512 + mu * (2_048 + 8 * (ku + nu));
        let mut progress = WorkProgress::new(cx, WorkPlan::checked(0, mu, 0, 0, work, 0, hashing)?)?;
        progress.checkpoint("shared covariance admission")?;
        for (index, record) in records.iter().enumerate() {
            for earlier in &records[..index] {
                progress.comparison_bytes("shared row uniqueness", 256)?;
                if earlier.row == record.row {
                    return Err(AssimError::DuplicateDatasetRow { id: record.row.as_str().to_owned() });
                }
            }
            if record.model_id != model.content_id {
                return Err(declaration(record.row.as_str(), "Jacobian is bound to a different shared model"));
            }
            if record.observation.operator().len() != n {
                return Err(AssimError::DimMismatch { state: n, operator: record.observation.operator().len() });
            }
            progress.records("shared covariance rows", 1)?;
        }
        let mut projected = zero_matrix(m, k, &mut progress)?;
        for (row, record) in records.iter().enumerate() {
            for column in 0..k {
                let mut value = record.jacobian[column];
                for parameter in 0..column {
                    let u = model.factor.upper_entry(parameter, column)
                        .ok_or_else(|| declaration(model.identity(), "invalid shared factor shape"))?;
                    value = u.mul_add(record.jacobian[parameter], value);
                    progress.scalar("shared covariance projection", 2)?;
                }
                let diagonal = model.factor.diag(column)
                    .ok_or_else(|| declaration(model.identity(), "invalid shared factor diagonal"))?;
                projected[row][column] = value * diagonal.sqrt();
                progress.scalar("shared covariance projection", 2)?;
                finite(projected[row][column])?;
            }
        }
        let mut covariance = zero_matrix(m, m, &mut progress)?;
        for row in 0..m {
            for column in row..m {
                let mut value = if row == column { records[row].observation.noise_var() } else { 0.0 };
                for parameter in 0..k {
                    value = projected[row][parameter].mul_add(projected[column][parameter], value);
                    progress.scalar("shared covariance Gram product", 2)?;
                }
                finite(value)?;
                covariance[row][column] = canonicalize_zero(value);
                covariance[column][row] = canonicalize_zero(value);
            }
        }
        let mut hash = DomainHasher::new(BATCH_DOMAIN);
        field(&mut hash, model.content_id.as_bytes(), &mut progress)?;
        field(&mut hash, &(m as u64).to_le_bytes(), &mut progress)?;
        let mut observations = Vec::with_capacity(m);
        for (index, record) in records.iter().enumerate() {
            let observation = &record.observation;
            field(&mut hash, record.row.as_str().as_bytes(), &mut progress)?;
            field(&mut hash, observation.instrument().as_bytes(), &mut progress)?;
            field(&mut hash, &(n as u64).to_le_bytes(), &mut progress)?;
            for value in observation.operator() { number(&mut hash, *value, &mut progress)?; }
            number(&mut hash, observation.value(), &mut progress)?;
            number(&mut hash, observation.noise_var(), &mut progress)?;
            for value in &record.jacobian { number(&mut hash, *value, &mut progress)?; }
            // Framed hashing avoids ambiguous instrument@row concatenations
            // and accommodates leaf identities already at their length cap.
            let mut row_hash = DomainHasher::new(ROW_DOMAIN);
            field(&mut row_hash, observation.instrument().as_bytes(), &mut progress)?;
            field(&mut row_hash, record.row.as_str().as_bytes(), &mut progress)?;
            observations.push(RobustObservation::available(Observation::new(
                observation.operator().to_vec(), observation.value(), covariance[index][index],
                format!("jacobian-row:{}", row_hash.finalize()),
            )?));
        }
        progress.checkpoint("shared covariance publication")?;
        let numerical_batch = ObservationBatch::new(observations, covariance, cx)?;
        progress.checkpoint("shared covariance publication")?;
        Ok(Self { model, records, numerical_batch, content_id: format!("jacobian-batch:v1:{}", hash.finalize()) })
    }

    /// The nuisance model and its fixed parameter ordering.
    #[must_use]
    pub const fn model(&self) -> &SharedNoiseModel { &self.model }

    /// Original row-unique observations and their model-bound Jacobians.
    #[must_use]
    pub fn records(&self) -> &[JacobianRecord] { &self.records }

    /// Derived total observation covariance, including off-diagonal terms.
    #[must_use]
    pub fn covariance(&self) -> &[Vec<f64>] { self.numerical_batch.covariance() }

    /// Content identity of the model, row identities, readings and Jacobians.
    #[must_use]
    pub fn content_id(&self) -> &str { &self.content_id }

    /// Fuse this model through the checked correlated factor path. The result
    /// identity additionally binds the model-bearing batch to the numerical
    /// audit, which already binds the prior and the published posterior.
    ///
    /// # Errors
    /// Propagates dimensional, numerical, budget, and cancellation refusals.
    pub fn assimilate(&self, prior: &Belief, cx: &Cx<'_>) -> Result<JacobianAssimilation, AssimError> {
        let robust = assimilate_observation_batch(prior, &self.numerical_batch, cx)?;
        let mut progress = WorkProgress::new(cx, WorkPlan::checked(0, 0, 0, 0, 0, 0, 512)?)?;
        progress.checkpoint("shared posterior identity")?;
        let mut hash = DomainHasher::new(RESULT_DOMAIN);
        field(&mut hash, self.content_id.as_bytes(), &mut progress)?;
        field(&mut hash, robust.audit().receipt_id().as_bytes(), &mut progress)?;
        progress.checkpoint("shared posterior publication")?;
        Ok(JacobianAssimilation { robust, content_id: format!("jacobian-result:v1:{}", hash.finalize()) })
    }
}

/// A posterior with both numerical and nuisance-model identity retained.
#[derive(Debug, Clone, PartialEq)]
pub struct JacobianAssimilation {
    robust: RobustAssimilation,
    content_id: String,
}

impl JacobianAssimilation {
    /// Published posterior and numerical audit (including factor diagnostics).
    #[must_use]
    pub const fn robust(&self) -> &RobustAssimilation { &self.robust }

    /// Identity binding the nuisance model and row Jacobians to the update.
    #[must_use]
    pub fn content_id(&self) -> &str { &self.content_id }
}

fn declaration(id: &str, reason: &'static str) -> AssimError {
    AssimError::SharedSourceDeclaration { id: id.to_owned(), reason }
}

fn finite(value: f64) -> Result<(), AssimError> {
    if value.is_finite() { Ok(()) }
    else { Err(AssimError::NonFiniteComputation { stage: "shared covariance projection" }) }
}

fn zero_matrix(rows: usize, columns: usize, progress: &mut WorkProgress<'_, '_>) -> Result<Vec<Vec<f64>>, AssimError> {
    let refused = || AssimError::AllocationRefused {
        stage: "shared covariance matrix", requested_bytes: rows as u128 * (core::mem::size_of::<Vec<f64>>() as u128 + 8 * columns as u128),
    };
    let mut result = Vec::new();
    result.try_reserve_exact(rows).map_err(|_| refused())?;
    for _ in 0..rows {
        let mut row = Vec::new();
        row.try_reserve_exact(columns).map_err(|_| refused())?;
        row.resize(columns, 0.0);
        progress.scalar("shared covariance allocation", columns as u128)?;
        result.push(row);
    }
    Ok(result)
}

fn field(hash: &mut DomainHasher, bytes: &[u8], progress: &mut WorkProgress<'_, '_>) -> Result<(), AssimError> {
    hash.update(&(bytes.len() as u64).to_le_bytes());
    hash.update(bytes);
    progress.hash_bytes("shared noise identity", bytes.len() as u128 + 8)
}

fn number(hash: &mut DomainHasher, value: f64, progress: &mut WorkProgress<'_, '_>) -> Result<(), AssimError> {
    hash.update(&canonical_f64_bits(value).to_le_bytes());
    progress.hash_bytes("shared noise identity", 8)
}
