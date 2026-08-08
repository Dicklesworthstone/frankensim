//! Dimensional adoption for the assimilation lane (bead sj31i.7.3):
//! mechanical binding of the shared inference dimensional core
//! (`fs_qty::inference`) to beliefs, observations, H/R algebra,
//! innovations, gains, residuals, misfit, and receipts.
//!
//! The raw `Belief`/`Observation` constructors stay untouched (candidate
//! identity v4 pins their exact behavior); this module is the dimensionally
//! typed front door. A [`DimensionedBelief`] binds a [`StateSchema`] to a
//! validated belief; a [`DimensionedObservation`] is admitted only when its
//! operator and noise variance satisfy the schema algebra — affine
//! absolute-temperature state slots are refused through linear operators,
//! covariance/information dimensions are derived mechanically, and the
//! weighted misfit is verified dimensionless. The assimilation receipt
//! binds the schema identity, the numeric policy version, and the exact
//! posterior bits.

use fs_qty::Dims;
use fs_qty::inference::{
    CovarianceSchema, InferenceError, ObservationSchema, OperatorSchema, StateSchema,
};
use fs_qty::semantic::QuantitySpec;

use crate::{AssimError, Belief, Observation, PSD_ADMISSION_POLICY_VERSION, assimilate};

/// Wire prefix of the dimensioned assimilation receipt identity.
pub const DIMENSIONED_RECEIPT_PREFIX: &str = "dimensioned-assimilation:v1:";

const RECEIPT_DOMAIN: &str = "org.frankensim.fs-assimilate.dimensioned-assimilation.v1";

fn dimensional(error: InferenceError) -> AssimError {
    AssimError::Dimensional(error)
}

/// A belief with a mechanically bound state schema: mean and covariance
/// shapes are checked against the schema, and every downstream dimensional
/// question (covariance entry, information entry, gain, residual) is
/// answered by the schema algebra rather than by convention.
#[derive(Debug, Clone, PartialEq)]
pub struct DimensionedBelief {
    state_schema: StateSchema,
    belief: Belief,
}

impl DimensionedBelief {
    /// Bind a validated belief to a state schema. The schema declares the
    /// dimension of every state slot; the belief must match it in shape.
    ///
    /// # Errors
    /// Returns [`AssimError::DimMismatch`] when the belief dimension
    /// differs from the schema slot count.
    pub fn try_new(state_schema: StateSchema, belief: Belief) -> Result<Self, AssimError> {
        if belief.dim() != state_schema.len() {
            return Err(AssimError::DimMismatch {
                state: state_schema.len(),
                operator: belief.dim(),
            });
        }
        Ok(Self {
            state_schema,
            belief,
        })
    }

    /// The bound state schema.
    #[must_use]
    pub const fn state_schema(&self) -> &StateSchema {
        &self.state_schema
    }

    /// The wrapped belief.
    #[must_use]
    pub const fn belief(&self) -> &Belief {
        &self.belief
    }

    /// Mechanically derived dimensions of covariance entry `(i, j)`:
    /// `dims(i) * dims(j)`.
    ///
    /// # Errors
    /// Returns the typed dimensional error for out-of-range slots or
    /// exponent overflow.
    pub fn covariance_entry_dims(&self, row: usize, column: usize) -> Result<Dims, AssimError> {
        CovarianceSchema::over(self.state_schema.clone())
            .entry_dims(row, column)
            .map_err(dimensional)
    }

    /// Mechanically derived dimensions of information (inverse covariance)
    /// entry `(i, j)`.
    ///
    /// # Errors
    /// See [`Self::covariance_entry_dims`].
    pub fn information_entry_dims(&self, row: usize, column: usize) -> Result<Dims, AssimError> {
        CovarianceSchema::over(self.state_schema.clone())
            .information_entry_dims(row, column)
            .map_err(dimensional)
    }

    /// State slots carrying affine absolute temperature: these hold affine
    /// points, while increments, gains, and covariances live in the
    /// difference algebra. An observation operator cannot touch them until
    /// the slot is converted with
    /// [`fs_qty::inference::SlotSchema::as_difference`].
    #[must_use]
    pub fn affine_absolute_slots(&self) -> Vec<usize> {
        self.state_schema
            .slots()
            .iter()
            .enumerate()
            .filter_map(|(slot, schema)| schema.is_affine_absolute_temperature().then_some(slot))
            .collect()
    }
}

/// A dimensionally admitted scalar observation: the reading schema, the
/// operator schema over a state schema, and the checked raw observation.
/// Construction refuses dimensionally inconsistent or affine-trapped
/// declarations; there is no untyped path into a dimensioned update.
#[derive(Debug, Clone, PartialEq)]
pub struct DimensionedObservation {
    observation: Observation,
    observation_schema: ObservationSchema,
    operator_schema: OperatorSchema,
}

impl DimensionedObservation {
    /// Admit one scalar observation against a state schema. The operator
    /// row maps state to reading: column `i` carries
    /// `reading_dims - state_dims(i)`, the noise variance carries the
    /// squared reading dimensions, and an affine absolute-temperature state
    /// slot refuses admission outright.
    ///
    /// # Errors
    /// Returns [`AssimError`] for operator-length mismatch, a malformed raw
    /// observation, or the typed dimensional refusals of the schema
    /// algebra.
    pub fn try_new(
        state_schema: &StateSchema,
        output: QuantitySpec,
        operator: Vec<f64>,
        value: f64,
        noise_variance: f64,
        instrument: impl Into<String>,
    ) -> Result<Self, AssimError> {
        if operator.len() != state_schema.len() {
            return Err(AssimError::DimMismatch {
                state: state_schema.len(),
                operator: operator.len(),
            });
        }
        let observation_schema = ObservationSchema::new(output);
        let operator_schema = OperatorSchema::try_new(observation_schema, state_schema.clone())
            .map_err(dimensional)?;
        let observation = Observation::new(operator, value, noise_variance, instrument)?;
        Ok(Self {
            observation,
            observation_schema,
            operator_schema,
        })
    }

    /// The checked raw observation.
    #[must_use]
    pub const fn observation(&self) -> &Observation {
        &self.observation
    }

    /// The reading schema.
    #[must_use]
    pub const fn observation_schema(&self) -> ObservationSchema {
        self.observation_schema
    }

    /// Reading dimensions.
    #[must_use]
    pub fn reading_dims(&self) -> Dims {
        self.observation_schema.value_dims()
    }

    /// Noise-variance dimensions (the squared reading dimensions).
    ///
    /// # Errors
    /// Returns the typed dimensional error on exponent overflow.
    pub fn noise_variance_dims(&self) -> Result<Dims, AssimError> {
        self.observation_schema
            .noise_variance_dims()
            .map_err(dimensional)
    }

    /// Residual/innovation dimensions (the reading dimensions).
    #[must_use]
    pub fn residual_dims(&self) -> Dims {
        self.observation_schema.residual_dims()
    }

    /// Mechanically derived dimensions of operator column `i`.
    ///
    /// # Errors
    /// Returns the typed dimensional error for out-of-range slots or
    /// exponent overflow.
    pub fn operator_column_dims(&self, slot: usize) -> Result<Dims, AssimError> {
        self.operator_schema.column_dims(slot).map_err(dimensional)
    }

    /// Mechanically derived Kalman-gain dimensions for state slot `i`:
    /// `state_dims(i) - reading_dims`.
    ///
    /// # Errors
    /// Returns the typed dimensional error for out-of-range slots or
    /// exponent overflow.
    pub fn gain_dims(&self, slot: usize) -> Result<Dims, AssimError> {
        let state_dims = self
            .operator_schema
            .state()
            .slot(slot)
            .map_err(dimensional)?
            .dims();
        state_dims
            .checked_minus(self.reading_dims())
            .ok_or(AssimError::Dimensional(InferenceError::DimensionOverflow {
                stage: "gain",
            }))
    }
}

/// The dimensioned assimilation product: the posterior belief with its
/// state schema (unchanged by a linear update), the verified-dimensionless
/// weighted misfits, and a receipt identity binding the schema, the numeric
/// policy, and the exact posterior bits.
#[derive(Debug, Clone, PartialEq)]
pub struct DimensionedAssimilation {
    posterior: DimensionedBelief,
    innovation_dims: Dims,
    weighted_misfit_before: f64,
    weighted_misfit_after: f64,
    receipt_identity: String,
}

impl DimensionedAssimilation {
    /// The posterior belief with its bound schema.
    #[must_use]
    pub const fn posterior(&self) -> &DimensionedBelief {
        &self.posterior
    }

    /// Innovation/residual dimensions of the fused observation.
    #[must_use]
    pub const fn innovation_dims(&self) -> Dims {
        self.innovation_dims
    }

    /// Weighted scalar misfit before the update (verified dimensionless).
    #[must_use]
    pub const fn weighted_misfit_before(&self) -> f64 {
        self.weighted_misfit_before
    }

    /// Weighted scalar misfit after the update (verified dimensionless).
    #[must_use]
    pub const fn weighted_misfit_after(&self) -> f64 {
        self.weighted_misfit_after
    }

    /// The dimensioned receipt identity
    /// `dimensioned-assimilation:v1:<64 lowercase hex>`.
    #[must_use]
    pub fn receipt_identity(&self) -> &str {
        &self.receipt_identity
    }
}

fn receipt_identity(
    prior: &DimensionedBelief,
    obs: &DimensionedObservation,
    posterior: &Belief,
    mode: &'static str,
) -> String {
    let mut hasher = fs_blake3::DomainHasher::new(RECEIPT_DOMAIN);
    hasher.update(prior.state_schema().identity().as_bytes());
    hasher.update(&(obs.observation().operator().len() as u64).to_le_bytes());
    for entry in obs.observation().operator() {
        hasher.update(&entry.to_le_bytes());
    }
    hasher.update(&obs.observation().value().to_le_bytes());
    hasher.update(&obs.observation().noise_var().to_le_bytes());
    hasher.update(obs.observation().instrument().as_bytes());
    hasher.update(&PSD_ADMISSION_POLICY_VERSION.to_le_bytes());
    hasher.update(mode.as_bytes());
    for entry in posterior.mean() {
        hasher.update(&entry.to_le_bytes());
    }
    for row in posterior.covariance() {
        for entry in row {
            hasher.update(&entry.to_le_bytes());
        }
    }
    format!("{DIMENSIONED_RECEIPT_PREFIX}{}", hasher.finalize())
}

/// Assimilate one dimensionally admitted observation into a dimensioned
/// belief through the production scalar updater, returning the posterior
/// with its schema, the verified-dimensionless weighted misfits, and the
/// dimensionally bound receipt identity.
///
/// The weighted misfit law is mechanical: residual squared has the squared
/// reading dimensions, the noise variance has the squared reading
/// dimensions, and their quotient is verified dimensionless before the
/// update runs — the check is schema algebra, not a prose claim. The state
/// schema is invariant under the linear update, so the posterior shares the
/// prior's schema.
///
/// # Errors
/// Returns [`AssimError`] from the production update (malformed state,
/// degenerate innovation, budget, or cancellation) and from the
/// dimensional-misfit check.
pub fn assimilate_dimensioned(
    prior: &DimensionedBelief,
    obs: &DimensionedObservation,
    cx: &crate::Cx<'_>,
) -> Result<DimensionedAssimilation, AssimError> {
    // Mechanical misfit dimension check: residual^2 / noise_var must be
    // dimensionless before any arithmetic runs.
    let residual_squared = obs
        .residual_dims()
        .checked_times(2)
        .ok_or(AssimError::Dimensional(InferenceError::DimensionOverflow {
            stage: "misfit residual square",
        }))?;
    let noise_dims = obs.noise_variance_dims()?;
    if residual_squared != noise_dims {
        return Err(AssimError::Dimensional(InferenceError::DimensionMismatch {
            role: "weighted misfit",
            expected: residual_squared.unit_string(),
            actual: noise_dims.unit_string(),
        }));
    }
    let innovation_dims = obs.residual_dims();

    let misfit_before =
        crate::misfit(prior.belief(), core::slice::from_ref(obs.observation()), cx)?;
    let posterior = assimilate(prior.belief(), obs.observation(), cx)?;
    let misfit_after = crate::misfit(&posterior, core::slice::from_ref(obs.observation()), cx)?;
    let identity = receipt_identity(prior, obs, &posterior, cx.mode().name());
    Ok(DimensionedAssimilation {
        posterior: DimensionedBelief {
            state_schema: prior.state_schema.clone(),
            belief: posterior,
        },
        innovation_dims,
        weighted_misfit_before: misfit_before,
        weighted_misfit_after: misfit_after,
        receipt_identity: identity,
    })
}
