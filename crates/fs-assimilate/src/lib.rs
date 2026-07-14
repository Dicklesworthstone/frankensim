//! fs-assimilate — validation as a living belief (plan addendum, Proposal 11).
//! Layer: L4.
//!
//! Strain-gauge and wind-tunnel data update the MODEL-FORM POSTERIOR that
//! Proposal 3 tracks per regime, so "validated" stops being a one-time stamp
//! and becomes a living belief state. A sensor readout is a TRACE of the field
//! onto the sensor's support — an observation operator expressed in the same
//! restriction-map algebra as the sheaf.
//!
//! This crate is the linear-Gaussian core of that assimilation: a [`Belief`]
//! (Gaussian state) is updated by [`Observation`]s (a restriction-map row + a
//! reading + its instrument noise) via sequential Kalman fusion. Two honest
//! properties:
//! - POINT SENSORS ([`point_sensor`]) are the REGISTRATION-FREE path (the R8
//!   fallback): their observation operator picks a state component directly, so
//!   they work even where full-field scan integration is premature. Scan
//!   observations ([`scan_observation`]) carry the registration variance too.
//! - Assimilation produces an **estimated candidate** tied to a proposed regime
//!   ([`assimilate_colored`]). Experimental validation is a separate admission
//!   act requiring calibrated data and an external authenticated authority.
//!
//! Deterministic and cooperatively cancellable through [`fs_exec::Cx`].

use core::fmt;
use core::mem::size_of;

pub use fs_evidence::{Color, ValidityDomain};
use fs_exec::Cx;

const CANDIDATE_ID_DOMAIN: &str = "org.frankensim.fs-assimilate.candidate.v3";
const CANDIDATE_ID_PREFIX: &str = "assimilation-candidate:v3:";
const SCALAR_POLL_STRIDE: u128 = 256;
const RECORD_POLL_STRIDE: u128 = 16;
const CANONICAL_COMPARE_BYTE_POLL_STRIDE: u128 = 1_024;
const HASH_BYTE_POLL_STRIDE: usize = 1_024;
const POLL_POLICY_ID: &str = "fixed-stride:v2";
/// Maximum state dimension admitted by the synchronous dense v0 core.
///
/// The Joseph update is `O(n^3)` and owns several `n x n` work matrices. Larger
/// states belong on a sparse or matrix-free, cancellable assimilation path.
pub const MAX_DENSE_STATE_DIM: usize = 256;
/// Maximum observations admitted by one synchronous dense aggregate call.
///
/// This also bounds canonical-order sorting and candidate-identity materialization
/// for low-dimensional campaigns. High-rate streams belong in a cancellable,
/// incremental assimilation session rather than one monolithic call.
pub const MAX_DENSE_OBSERVATIONS: usize = 4_096;
/// Maximum `observation_count * state_dimension^3` work proxy admitted by one
/// dense aggregate update.
///
/// The Joseph covariance update is cubic in the state dimension. A count cap by
/// itself would still admit tens of billions of dense operations at the largest
/// state, so aggregate admission must bound the multiplicative workload too.
pub const MAX_DENSE_UPDATE_CUBIC_WORK: u128 = 4 * 256_u128 * 256_u128 * 256_u128;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct WorkPlan {
    validation_psd: u128,
    record_materialization: u128,
    canonical_ordering: u128,
    misfit_passes: u128,
    joseph_update: u128,
    psd_revalidation: u128,
    hashing: u128,
    total: u128,
}

impl WorkPlan {
    fn checked(
        validation_psd: u128,
        record_materialization: u128,
        canonical_ordering: u128,
        misfit_passes: u128,
        joseph_update: u128,
        psd_revalidation: u128,
        hashing: u128,
    ) -> Result<Self, AssimError> {
        let total = [
            validation_psd,
            record_materialization,
            canonical_ordering,
            misfit_passes,
            joseph_update,
            psd_revalidation,
            hashing,
        ]
        .into_iter()
        .try_fold(0_u128, u128::checked_add)
        .ok_or(AssimError::WorkPlanOverflow { phase: "total" })?;
        Ok(Self {
            validation_psd,
            record_materialization,
            canonical_ordering,
            misfit_passes,
            joseph_update,
            psd_revalidation,
            hashing,
            total,
        })
    }
}

struct WorkProgress<'a, 's> {
    cx: &'a Cx<'s>,
    completed: u128,
    planned: u128,
    scalar_since_poll: u128,
    records_since_poll: u128,
    comparison_bytes_since_poll: u128,
    hash_bytes_since_poll: u128,
    initial_poll_quota: u32,
    polls_remaining: u32,
    shared_polls_remaining: Option<&'a mut u32>,
}

impl<'a, 's> WorkProgress<'a, 's> {
    fn new(cx: &'a Cx<'s>, plan: WorkPlan) -> Self {
        Self {
            cx,
            completed: 0,
            planned: plan.total,
            scalar_since_poll: 0,
            records_since_poll: 0,
            comparison_bytes_since_poll: 0,
            hash_bytes_since_poll: 0,
            initial_poll_quota: cx.budget().poll_quota,
            polls_remaining: cx.budget().poll_quota,
            shared_polls_remaining: None,
        }
    }

    fn new_with_shared_poll_quota(
        cx: &'a Cx<'s>,
        plan: WorkPlan,
        polls_remaining: &'a mut u32,
    ) -> Self {
        Self {
            cx,
            completed: 0,
            planned: plan.total,
            scalar_since_poll: 0,
            records_since_poll: 0,
            comparison_bytes_since_poll: 0,
            hash_bytes_since_poll: 0,
            initial_poll_quota: *polls_remaining,
            polls_remaining: *polls_remaining,
            shared_polls_remaining: Some(polls_remaining),
        }
    }

    fn checkpoint(&mut self, phase: &'static str) -> Result<(), AssimError> {
        if self.polls_remaining == 0 {
            return Err(self.cancelled(phase));
        }
        if self.polls_remaining != u32::MAX {
            self.polls_remaining -= 1;
            if let Some(shared) = self.shared_polls_remaining.as_deref_mut() {
                *shared = self.polls_remaining;
            }
        }
        self.cx.checkpoint().map_err(|_| self.cancelled(phase))
    }

    fn scalar(&mut self, phase: &'static str, units: u128) -> Result<(), AssimError> {
        self.advance(phase, units)?;
        self.scalar_since_poll = self
            .scalar_since_poll
            .checked_add(units)
            .ok_or(AssimError::WorkPlanOverflow { phase })?;
        while self.scalar_since_poll >= SCALAR_POLL_STRIDE {
            self.scalar_since_poll -= SCALAR_POLL_STRIDE;
            self.checkpoint(phase)?;
        }
        Ok(())
    }

    fn records(&mut self, phase: &'static str, units: u128) -> Result<(), AssimError> {
        self.advance(phase, units)?;
        self.records_since_poll = self
            .records_since_poll
            .checked_add(units)
            .ok_or(AssimError::WorkPlanOverflow { phase })?;
        while self.records_since_poll >= RECORD_POLL_STRIDE {
            self.records_since_poll -= RECORD_POLL_STRIDE;
            self.checkpoint(phase)?;
        }
        Ok(())
    }

    fn hash_bytes(&mut self, phase: &'static str, units: u128) -> Result<(), AssimError> {
        self.advance(phase, units)?;
        self.hash_bytes_since_poll = self
            .hash_bytes_since_poll
            .checked_add(units)
            .ok_or(AssimError::WorkPlanOverflow { phase })?;
        while self.hash_bytes_since_poll >= HASH_BYTE_POLL_STRIDE as u128 {
            self.hash_bytes_since_poll -= HASH_BYTE_POLL_STRIDE as u128;
            self.checkpoint(phase)?;
        }
        Ok(())
    }

    fn comparison_bytes(&mut self, phase: &'static str, units: u128) -> Result<(), AssimError> {
        self.advance(phase, units)?;
        self.comparison_bytes_since_poll = self
            .comparison_bytes_since_poll
            .checked_add(units)
            .ok_or(AssimError::WorkPlanOverflow { phase })?;
        while self.comparison_bytes_since_poll >= CANONICAL_COMPARE_BYTE_POLL_STRIDE {
            self.comparison_bytes_since_poll -= CANONICAL_COMPARE_BYTE_POLL_STRIDE;
            self.checkpoint(phase)?;
        }
        Ok(())
    }

    fn advance(&mut self, phase: &'static str, units: u128) -> Result<(), AssimError> {
        let attempted = self
            .completed
            .checked_add(units)
            .ok_or(AssimError::WorkPlanOverflow { phase })?;
        if attempted > self.planned {
            return Err(AssimError::WorkPlanExceeded {
                phase,
                attempted,
                planned: self.planned,
            });
        }
        self.completed = attempted;
        Ok(())
    }

    fn cancelled(&self, phase: &'static str) -> AssimError {
        AssimError::Cancelled {
            phase,
            completed: self.completed,
            planned: self.planned,
        }
    }
}

/// A Gaussian belief over an `n`-dimensional state.
///
/// Construction is checked so the mean is finite and the covariance is finite,
/// square, symmetric, and positive semidefinite. Fields stay private so a
/// checked belief cannot later be made ragged or non-finite.
#[derive(Debug, Clone, PartialEq)]
pub struct Belief {
    mean: Vec<f64>,
    cov: Vec<Vec<f64>>,
}

impl Belief {
    /// Construct a checked belief from a mean and full covariance matrix.
    ///
    /// # Errors
    /// Returns [`AssimError`] when the state is empty or any covariance
    /// invariant is violated.
    pub fn new(
        mut mean: Vec<f64>,
        mut cov: Vec<Vec<f64>>,
        cx: &Cx<'_>,
    ) -> Result<Self, AssimError> {
        preflight_belief_shape(&mean, &cov)?;
        let plan = belief_validation_work_plan(mean.len())?;
        let mut progress = WorkProgress::new(cx, plan);
        progress.checkpoint("initial")?;
        validate_belief_parts(&mean, &cov, "belief-validation", &mut progress)?;
        canonicalize_belief_zeros(&mut mean, &mut cov, &mut progress)?;
        progress.checkpoint("finalize")?;
        Ok(Self { mean, cov })
    }

    /// Construct a checked 1-D belief `N(mean, var)`.
    ///
    /// # Errors
    /// Returns [`AssimError`] when `mean` is non-finite or `var` is non-finite
    /// or negative.
    pub fn scalar(mean: f64, var: f64) -> Result<Self, AssimError> {
        if !mean.is_finite() {
            return Err(AssimError::NonFiniteMean { index: 0 });
        }
        if !var.is_finite() {
            return Err(AssimError::NonFiniteCovariance { row: 0, column: 0 });
        }
        if var < 0.0 {
            return Err(AssimError::NegativeVariance { index: 0 });
        }
        Ok(Self {
            mean: vec![canonicalize_zero(mean)],
            cov: vec![vec![canonicalize_zero(var)]],
        })
    }

    /// Construct a checked independent (diagonal-covariance) belief.
    ///
    /// # Errors
    /// Returns [`AssimError`] when the vectors have different lengths, are
    /// empty, contain non-finite values, or contain a negative variance.
    pub fn diagonal(mut means: Vec<f64>, vars: &[f64], cx: &Cx<'_>) -> Result<Self, AssimError> {
        if means.len() != vars.len() {
            return Err(AssimError::DiagonalDimensionMismatch {
                means: means.len(),
                variances: vars.len(),
            });
        }
        if means.is_empty() {
            return Err(AssimError::EmptyBelief);
        }
        validate_state_dimension(means.len())?;
        let _matrix_entries = checked_square_usize(means.len(), "diagonal covariance")?;
        let plan = belief_validation_work_plan(means.len())?;
        let mut progress = WorkProgress::new(cx, plan);
        progress.checkpoint("initial")?;
        for (index, value) in means.iter().enumerate() {
            progress.scalar("belief-validation", 1)?;
            if !value.is_finite() {
                return Err(AssimError::NonFiniteMean { index });
            }
        }
        for (index, variance) in vars.iter().enumerate() {
            progress.scalar("belief-validation", 1)?;
            if !variance.is_finite() {
                return Err(AssimError::NonFiniteCovariance {
                    row: index,
                    column: index,
                });
            }
            if *variance < 0.0 {
                return Err(AssimError::NegativeVariance { index });
            }
        }
        let mut cov = zero_matrix(means.len(), "diagonal-materialization", &mut progress)?;
        for (i, &variance) in vars.iter().enumerate() {
            cov[i][i] = variance;
            progress.scalar("diagonal-materialization", 1)?;
        }
        validate_belief_parts(&means, &cov, "belief-validation", &mut progress)?;
        canonicalize_belief_zeros(&mut means, &mut cov, &mut progress)?;
        progress.checkpoint("finalize")?;
        Ok(Self { mean: means, cov })
    }

    /// Recheck every structural and numerical belief invariant.
    ///
    /// # Errors
    /// Returns the first violated invariant.
    pub fn validate(&self, cx: &Cx<'_>) -> Result<(), AssimError> {
        preflight_belief_shape(&self.mean, &self.cov)?;
        let plan = belief_validation_work_plan(self.dim())?;
        let mut progress = WorkProgress::new(cx, plan);
        progress.checkpoint("initial")?;
        validate_belief_parts(&self.mean, &self.cov, "belief-validation", &mut progress)?;
        progress.checkpoint("finalize")
    }

    fn from_covariance_preserving_update(
        mut mean: Vec<f64>,
        mut cov: Vec<Vec<f64>>,
        progress: &mut WorkProgress<'_, '_>,
    ) -> Result<Self, AssimError> {
        // Floating-point evaluation does not inherit the exact-arithmetic PSD
        // closure law automatically. Route every computed posterior through
        // the same fail-closed boundary as an externally supplied belief.
        validate_belief_parts(&mean, &cov, "posterior-psd", progress)?;
        canonicalize_belief_zeros(&mut mean, &mut cov, progress)?;
        Ok(Self { mean, cov })
    }

    /// The state dimension.
    #[must_use]
    pub fn dim(&self) -> usize {
        self.mean.len()
    }

    /// Read-only view of the state mean.
    #[must_use]
    pub fn mean(&self) -> &[f64] {
        &self.mean
    }

    /// Read-only view of the covariance matrix.
    #[must_use]
    pub fn covariance(&self) -> &[Vec<f64>] {
        &self.cov
    }

    /// The mean of state component `component`.
    ///
    /// # Errors
    /// Returns [`AssimError::ComponentOutOfRange`] for an invalid component.
    pub fn component_mean(&self, component: usize) -> Result<f64, AssimError> {
        self.mean
            .get(component)
            .copied()
            .ok_or(AssimError::ComponentOutOfRange {
                component,
                dim: self.dim(),
            })
    }

    /// The variance of state component `component`.
    ///
    /// # Errors
    /// Returns [`AssimError::ComponentOutOfRange`] for an invalid component.
    pub fn variance(&self, component: usize) -> Result<f64, AssimError> {
        self.cov
            .get(component)
            .and_then(|row| row.get(component))
            .copied()
            .ok_or(AssimError::ComponentOutOfRange {
                component,
                dim: self.dim(),
            })
    }
}

/// One scalar observation: `value = operator · state + noise`, where `operator`
/// is the restriction-map row (the sensor's trace) and `noise_var` is the
/// instrument (+ registration) variance.
///
/// Construction is checked and fields stay private, preventing a valid
/// observation from being mutated into an empty, non-finite, or unanchored one.
#[derive(Debug, Clone, PartialEq)]
pub struct Observation {
    operator: Vec<f64>,
    value: f64,
    noise_var: f64,
    instrument: String,
}

impl Observation {
    /// Construct a checked scalar observation.
    ///
    /// # Errors
    /// Returns [`AssimError`] for an empty, oversized, zero, or non-finite
    /// operator; a non-finite reading; non-positive noise; or an unusable
    /// instrument identity.
    pub fn new(
        operator: Vec<f64>,
        value: f64,
        noise_var: f64,
        instrument: impl Into<String>,
    ) -> Result<Self, AssimError> {
        let observation = Self {
            operator: operator.into_iter().map(canonicalize_zero).collect(),
            value: canonicalize_zero(value),
            noise_var: canonicalize_zero(noise_var),
            instrument: instrument.into(),
        };
        observation.validate()?;
        Ok(observation)
    }

    /// Recheck every observation invariant except equality with a particular
    /// belief dimension.
    ///
    /// # Errors
    /// Returns the first violated invariant.
    pub fn validate(&self) -> Result<(), AssimError> {
        if self.operator.is_empty() {
            return Err(AssimError::EmptyObservationOperator);
        }
        validate_state_dimension(self.operator.len())?;
        for (index, coefficient) in self.operator.iter().enumerate() {
            if !coefficient.is_finite() {
                return Err(AssimError::NonFiniteObservationOperator { index });
            }
        }
        if self
            .operator
            .iter()
            .all(|coefficient| canonical_f64_bits(*coefficient) == 0)
        {
            return Err(AssimError::ZeroObservationOperator);
        }
        if !self.value.is_finite() {
            return Err(AssimError::NonFiniteObservationValue);
        }
        validate_noise(self.noise_var)?;
        validate_leaf_identity("instrument", &self.instrument)
    }

    /// Read-only view of the observation operator.
    #[must_use]
    pub fn operator(&self) -> &[f64] {
        &self.operator
    }

    /// The observed scalar value.
    #[must_use]
    pub fn value(&self) -> f64 {
        self.value
    }

    /// The total observation noise variance.
    #[must_use]
    pub fn noise_var(&self) -> f64 {
        self.noise_var
    }

    /// The calibrated instrument identity.
    #[must_use]
    pub fn instrument(&self) -> &str {
        &self.instrument
    }
}

/// A registration-free point-sensor observation of state component `component`
/// (a strain gauge / thermocouple): its operator is the unit row `e_component`.
///
/// # Errors
/// Returns [`AssimError`] for a zero or oversized dimension, an out-of-range
/// component, or any malformed reading, noise, or instrument identity.
pub fn point_sensor(
    component: usize,
    dim: usize,
    value: f64,
    instrument_noise: f64,
    instrument: impl Into<String>,
) -> Result<Observation, AssimError> {
    if dim == 0 {
        return Err(AssimError::EmptyStateDimension);
    }
    validate_state_dimension(dim)?;
    if component >= dim {
        return Err(AssimError::ComponentOutOfRange { component, dim });
    }
    let mut operator = vec![0.0; dim];
    operator[component] = 1.0;
    Observation::new(operator, value, instrument_noise, instrument)
}

/// A full-field scan observation whose noise carries registration variance on
/// top of the strictly positive instrument variance (R8).
///
/// # Errors
/// Returns [`AssimError`] for malformed observation data, non-positive
/// instrument noise, negative registration variance, or an overflowing total.
pub fn scan_observation(
    operator: Vec<f64>,
    value: f64,
    instrument_noise: f64,
    registration_var: f64,
    instrument: impl Into<String>,
) -> Result<Observation, AssimError> {
    validate_noise(instrument_noise)?;
    if !registration_var.is_finite() {
        return Err(AssimError::NonFiniteRegistrationVariance);
    }
    if registration_var < 0.0 {
        return Err(AssimError::NegativeRegistrationVariance);
    }
    let noise_var = instrument_noise + registration_var;
    if !noise_var.is_finite() {
        return Err(AssimError::NonFiniteNoise);
    }
    Observation::new(operator, value, noise_var, instrument)
}

/// A structured assimilation failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AssimError {
    /// A belief cannot have zero state dimensions.
    EmptyBelief,
    /// A sensor cannot declare a zero state dimension.
    EmptyStateDimension,
    /// The dense synchronous v0 core refuses a state beyond its declared
    /// memory/compute envelope.
    StateDimensionLimit {
        /// Requested dimension.
        dim: usize,
        /// Maximum admitted dense dimension.
        max: usize,
    },
    /// An aggregate supplied too many observations for synchronous sorting,
    /// hashing, and iteration.
    ObservationCountLimit {
        /// Requested observation count.
        count: usize,
        /// Maximum admitted count.
        max: usize,
    },
    /// The observation-count/state-size product exceeds the bounded dense
    /// Joseph-update work envelope.
    AssimilationWorkLimit {
        /// Requested `observation_count * state_dimension^3` proxy units.
        requested: u128,
        /// Maximum admitted proxy units.
        max: u128,
    },
    /// Checked arithmetic could not represent the complete declared work
    /// shape before allocation.
    WorkPlanOverflow {
        /// Stable phase whose work expression overflowed.
        phase: &'static str,
    },
    /// Executed work exceeded the preflighted shape. This is an internal
    /// fail-closed accounting defect, never a partial-success condition.
    WorkPlanExceeded {
        /// Stable phase that attempted unplanned work.
        phase: &'static str,
        /// Work total that would have been recorded.
        attempted: u128,
        /// Complete preflighted work plan.
        planned: u128,
    },
    /// A compositional caller supplied more remaining polls than the ambient
    /// execution context admitted.
    PollQuotaExceedsAmbient {
        /// Caller-supplied remaining poll slice.
        requested: u32,
        /// Maximum slice admitted by the ambient execution context.
        ambient: u32,
    },
    /// The operation observed cancellation or exhausted its poll quota at a
    /// deterministic checkpoint. No partial belief or candidate is published.
    Cancelled {
        /// Stable phase containing the checkpoint.
        phase: &'static str,
        /// Declared work units completed before the checkpoint.
        completed: u128,
        /// Complete preflighted work plan.
        planned: u128,
    },
    /// The diagonal constructor received a different count of means and
    /// variances.
    DiagonalDimensionMismatch {
        /// Number of means.
        means: usize,
        /// Number of variances.
        variances: usize,
    },
    /// The covariance row count differs from the mean dimension.
    CovarianceDimensionMismatch {
        /// State dimension from the mean.
        state: usize,
        /// Covariance row count.
        rows: usize,
    },
    /// A covariance row is ragged.
    CovarianceRowDimensionMismatch {
        /// Offending row.
        row: usize,
        /// Required column count.
        expected: usize,
        /// Actual column count.
        actual: usize,
    },
    /// A mean component is NaN or infinite.
    NonFiniteMean {
        /// Offending component.
        index: usize,
    },
    /// A covariance entry is NaN or infinite.
    NonFiniteCovariance {
        /// Offending row.
        row: usize,
        /// Offending column.
        column: usize,
    },
    /// A diagonal covariance entry is negative.
    NegativeVariance {
        /// Offending component.
        index: usize,
    },
    /// A covariance pair is not exactly symmetric.
    NonSymmetricCovariance {
        /// Row of the upper-triangular entry.
        row: usize,
        /// Column of the upper-triangular entry.
        column: usize,
    },
    /// The symmetric covariance is not positive semidefinite.
    CovarianceNotPositiveSemidefinite,
    /// An observation operator has no coefficients.
    EmptyObservationOperator,
    /// An observation operator contains no state sensitivity.
    ZeroObservationOperator,
    /// An observation-operator coefficient is NaN or infinite.
    NonFiniteObservationOperator {
        /// Offending coefficient.
        index: usize,
    },
    /// The observed scalar value is NaN or infinite.
    NonFiniteObservationValue,
    /// An observation operator's length differs from the state dimension.
    DimMismatch {
        /// State dimension.
        state: usize,
        /// Operator length.
        operator: usize,
    },
    /// A requested state component is outside the declared dimension.
    ComponentOutOfRange {
        /// Requested component.
        component: usize,
        /// State dimension.
        dim: usize,
    },
    /// Observation noise is zero or negative.
    NonPositiveNoise,
    /// Observation noise is NaN or infinite, including overflow while combining
    /// instrument and registration variances.
    NonFiniteNoise,
    /// Registration variance is negative.
    NegativeRegistrationVariance,
    /// Registration variance is NaN or infinite.
    NonFiniteRegistrationVariance,
    /// An instrument identity is blank.
    EmptyInstrument,
    /// A regime-axis identity is blank.
    EmptyRegime,
    /// A machine-readable identity violates the shared evidence grammar.
    InvalidIdentity {
        /// Identity role (`instrument` or `regime_param`).
        field: &'static str,
        /// Stable rejection reason from `fs-evidence`.
        reason: &'static str,
    },
    /// An aggregate operation requires at least one observation.
    EmptyObservations,
    /// A regime bound is NaN or infinite.
    NonFiniteRegimeBounds,
    /// The regime lower bound exceeds its upper bound.
    InvertedRegimeBounds,
    /// The innovation covariance was non-positive (degenerate).
    SingularInnovation,
    /// Finite inputs overflowed or otherwise produced a non-finite intermediate.
    NonFiniteComputation {
        /// Stable computation stage.
        stage: &'static str,
    },
}

impl fmt::Display for AssimError {
    #[allow(clippy::too_many_lines)] // exhaustive structured diagnostics stay co-located
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyBelief => write!(f, "belief state dimension must be non-zero"),
            Self::EmptyStateDimension => write!(f, "sensor state dimension must be non-zero"),
            Self::StateDimensionLimit { dim, max } => {
                write!(f, "dense assimilation dimension {dim} exceeds limit {max}")
            }
            Self::ObservationCountLimit { count, max } => write!(
                f,
                "dense assimilation observation count {count} exceeds limit {max}"
            ),
            Self::AssimilationWorkLimit { requested, max } => write!(
                f,
                "dense assimilation work proxy {requested} exceeds limit {max}"
            ),
            Self::WorkPlanOverflow { phase } => {
                write!(
                    f,
                    "assimilation work-plan arithmetic overflowed during {phase}"
                )
            }
            Self::WorkPlanExceeded {
                phase,
                attempted,
                planned,
            } => write!(
                f,
                "assimilation work accounting exceeded its preflight during {phase}: \
                 attempted {attempted} units with {planned} planned"
            ),
            Self::PollQuotaExceedsAmbient { requested, ambient } => write!(
                f,
                "shared assimilation poll quota {requested} exceeds ambient quota {ambient}"
            ),
            Self::Cancelled {
                phase,
                completed,
                planned,
            } => write!(
                f,
                "assimilation cancelled during {phase} after {completed} of {planned} declared work units"
            ),
            Self::DiagonalDimensionMismatch { means, variances } => write!(
                f,
                "diagonal belief has {means} means but {variances} variances"
            ),
            Self::CovarianceDimensionMismatch { state, rows } => write!(
                f,
                "belief dimension is {state} but covariance has {rows} rows"
            ),
            Self::CovarianceRowDimensionMismatch {
                row,
                expected,
                actual,
            } => write!(
                f,
                "covariance row {row} has {actual} columns; expected {expected}"
            ),
            Self::NonFiniteMean { index } => {
                write!(f, "belief mean component {index} is non-finite")
            }
            Self::NonFiniteCovariance { row, column } => {
                write!(f, "covariance entry ({row}, {column}) is non-finite")
            }
            Self::NegativeVariance { index } => {
                write!(f, "covariance diagonal {index} is negative")
            }
            Self::NonSymmetricCovariance { row, column } => write!(
                f,
                "covariance entries ({row}, {column}) and ({column}, {row}) differ"
            ),
            Self::CovarianceNotPositiveSemidefinite => {
                write!(f, "covariance is not positive semidefinite")
            }
            Self::EmptyObservationOperator => write!(f, "observation operator must not be empty"),
            Self::ZeroObservationOperator => {
                write!(
                    f,
                    "observation operator must contain a non-zero coefficient"
                )
            }
            Self::NonFiniteObservationOperator { index } => {
                write!(f, "observation operator coefficient {index} is non-finite")
            }
            Self::NonFiniteObservationValue => write!(f, "observation value is non-finite"),
            Self::DimMismatch { state, operator } => write!(
                f,
                "state dimension is {state} but observation operator length is {operator}"
            ),
            Self::ComponentOutOfRange { component, dim } => {
                write!(f, "component {component} is outside state dimension {dim}")
            }
            Self::NonPositiveNoise => write!(f, "observation noise must be strictly positive"),
            Self::NonFiniteNoise => write!(f, "observation noise is non-finite"),
            Self::NegativeRegistrationVariance => {
                write!(f, "registration variance must be non-negative")
            }
            Self::NonFiniteRegistrationVariance => {
                write!(f, "registration variance is non-finite")
            }
            Self::EmptyInstrument => write!(f, "instrument identity must not be blank"),
            Self::EmptyRegime => write!(f, "regime axis identity must not be blank"),
            Self::InvalidIdentity { field, reason } => {
                write!(f, "invalid {field} identity: {reason}")
            }
            Self::EmptyObservations => write!(f, "at least one observation is required"),
            Self::NonFiniteRegimeBounds => write!(f, "regime bounds must be finite"),
            Self::InvertedRegimeBounds => {
                write!(f, "regime lower bound must not exceed its upper bound")
            }
            Self::SingularInnovation => {
                write!(f, "innovation covariance is non-positive")
            }
            Self::NonFiniteComputation { stage } => {
                write!(f, "assimilation produced a non-finite value during {stage}")
            }
        }
    }
}

impl std::error::Error for AssimError {}

/// The model-data misfit `Σⱼ (hⱼ·mean − yⱼ)² / rⱼ` — the weighted squared
/// residual assimilation seeks to reduce.
///
/// # Errors
/// Returns [`AssimError`] for an empty observation set, malformed input, a
/// dimension mismatch, or a non-finite computed term or sum.
pub fn misfit(
    belief: &Belief,
    observations: &[Observation],
    cx: &Cx<'_>,
) -> Result<f64, AssimError> {
    let plan = assimilation_work_plan(
        belief,
        observations,
        1,
        false,
        false,
        None,
        cx.mode().name().len(),
    )?;
    let mut progress = WorkProgress::new(cx, plan);
    progress.checkpoint("initial")?;
    let observations = validated_canonical_observations(observations, belief.dim(), &mut progress)?;
    let total = misfit_canonical(belief, &observations, "misfit", &mut progress)?;
    progress.checkpoint("finalize")?;
    Ok(total)
}

fn misfit_canonical(
    belief: &Belief,
    observations: &CanonicalObservations<'_>,
    phase: &'static str,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<f64, AssimError> {
    let mut total = 0.0;
    for record_index in observations.ordered_indices() {
        let observation = observations.records[record_index].observation;
        let predicted = checked_dot(
            &observation.operator,
            &belief.mean,
            "misfit prediction",
            phase,
            progress,
        )?;
        let residual = predicted - observation.value;
        progress.scalar(phase, 1)?;
        if !residual.is_finite() {
            return Err(AssimError::NonFiniteComputation {
                stage: "misfit residual",
            });
        }
        let term = residual * residual / observation.noise_var;
        progress.scalar(phase, 2)?;
        if !term.is_finite() {
            return Err(AssimError::NonFiniteComputation {
                stage: "misfit term",
            });
        }
        total += term;
        progress.scalar(phase, 1)?;
        if !total.is_finite() {
            return Err(AssimError::NonFiniteComputation {
                stage: "misfit sum",
            });
        }
    }
    Ok(total)
}

/// Fuse one observation into the belief by the scalar Kalman update. For a
/// valid covariance, every posterior component variance is at most its prior
/// value (information only increases).
///
/// # Errors
/// Returns [`AssimError`] for malformed input, a dimension mismatch, a
/// degenerate innovation, or a non-finite computed intermediate.
pub fn assimilate(prior: &Belief, obs: &Observation, cx: &Cx<'_>) -> Result<Belief, AssimError> {
    let observations = core::slice::from_ref(obs);
    let plan = assimilation_work_plan(
        prior,
        observations,
        0,
        true,
        false,
        None,
        cx.mode().name().len(),
    )?;
    let mut progress = WorkProgress::new(cx, plan);
    progress.checkpoint("initial")?;
    let observations = validated_canonical_observations(observations, prior.dim(), &mut progress)?;
    let belief = assimilate_canonical(prior, &observations, &mut progress)?;
    progress.checkpoint("finalize")?;
    Ok(belief)
}

/// Fuse all observations in their canonical content order. The mathematical
/// linear-Gaussian posterior is order-independent; canonical evaluation also
/// makes the floating-point result bit-stable across input permutations.
///
/// # Errors
/// Returns [`AssimError`] for an empty observation set or any error described by
/// [`assimilate`].
pub fn assimilate_all(
    prior: &Belief,
    observations: &[Observation],
    cx: &Cx<'_>,
) -> Result<Belief, AssimError> {
    let plan = assimilation_work_plan(
        prior,
        observations,
        0,
        true,
        false,
        None,
        cx.mode().name().len(),
    )?;
    let mut progress = WorkProgress::new(cx, plan);
    progress.checkpoint("initial")?;
    let observations = validated_canonical_observations(observations, prior.dim(), &mut progress)?;
    let belief = assimilate_canonical(prior, &observations, &mut progress)?;
    progress.checkpoint("finalize")?;
    Ok(belief)
}

fn assimilate_canonical(
    prior: &Belief,
    observations: &CanonicalObservations<'_>,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<Belief, AssimError> {
    validate_assimilation_work(prior.dim(), observations.len())?;
    let mut belief = clone_belief_for_update(prior, progress)?;
    for record_index in observations.ordered_indices() {
        let observation = observations.records[record_index].observation;
        belief = assimilate_checked(&belief, observation, progress)?;
    }
    Ok(belief)
}

/// An estimated, regime-tagged assimilated-posterior candidate.
///
/// The fields are read-only so this crate cannot accidentally expose a mutable
/// route from its honest estimated output to a stronger evidence color.
#[derive(Debug, Clone, PartialEq)]
pub struct AssimilatedPosterior {
    belief: Belief,
    color: Color,
    regime: ValidityDomain,
    misfit_before: f64,
    misfit_after: f64,
}

impl AssimilatedPosterior {
    /// The updated belief.
    #[must_use]
    pub fn belief(&self) -> &Belief {
        &self.belief
    }

    /// The honest estimated color of this candidate.
    #[must_use]
    pub fn color(&self) -> &Color {
        &self.color
    }

    /// The proposed regime for later experimental validation.
    #[must_use]
    pub fn regime(&self) -> &ValidityDomain {
        &self.regime
    }

    /// Model-data misfit before assimilation.
    #[must_use]
    pub fn misfit_before(&self) -> f64 {
        self.misfit_before
    }

    /// Model-data misfit after assimilation.
    #[must_use]
    pub fn misfit_after(&self) -> f64 {
        self.misfit_after
    }
}

/// Assimilate observations and return an instrument-bound **estimated**
/// candidate for a named regime — Proposal 3's living-belief update.
///
/// The candidate identity is a bounded domain-separated BLAKE3 digest over the
/// complete prior, the observation multiset (canonicalized independent of
/// input ordering), the proposed regime, the execution mode and logical stream
/// identity, every ambient budget field, the effective poll-quota slice, the
/// complete work plan, and the fixed poll policy. This function does not claim
/// that seeing data is itself
/// validation. Promotion to
/// [`Color::Validated`] belongs at an external admission boundary that
/// authenticates calibrated dataset provenance and validation authority.
///
/// # Errors
/// Returns [`AssimError`] for an invalid regime, an empty observation set, or
/// any malformed/non-finite assimilation input or result.
pub fn assimilate_colored(
    prior: &Belief,
    observations: &[Observation],
    regime_param: &str,
    regime_lo: f64,
    regime_hi: f64,
    cx: &Cx<'_>,
) -> Result<AssimilatedPosterior, AssimError> {
    let plan = assimilation_work_plan(
        prior,
        observations,
        2,
        true,
        true,
        Some((regime_param, regime_lo, regime_hi)),
        cx.mode().name().len(),
    )?;
    let mut progress = WorkProgress::new(cx, plan);
    assimilate_colored_planned(
        prior,
        observations,
        regime_param,
        regime_lo,
        regime_hi,
        cx,
        plan,
        &mut progress,
    )
}

/// Assimilate an estimated candidate while consuming a caller-owned remaining
/// poll quota in place.
///
/// This is the compositional form for a parent workflow that already charged
/// checkpoints against the same ambient [`Cx`] budget. The effective slice is
/// bound into the candidate identity in addition to the ambient budget, and
/// values above the ambient quota fail closed.
///
/// This low-level seam cannot authenticate the provenance of a raw counter. A
/// parent claiming one invocation-global quota must encapsulate the counter and
/// pass the same monotonically decreasing value to every nested call; replacing
/// or increasing it between calls starts a caller-authored slice outside that
/// no-reissue claim.
///
/// # Errors
/// Returns [`AssimError::PollQuotaExceedsAmbient`] when the supplied remaining
/// slice exceeds the ambient context's admitted quota, or the same structured
/// failures as [`assimilate_colored`] after admission.
pub fn assimilate_colored_with_shared_poll_quota(
    prior: &Belief,
    observations: &[Observation],
    regime_param: &str,
    regime_lo: f64,
    regime_hi: f64,
    cx: &Cx<'_>,
    polls_remaining: &mut u32,
) -> Result<AssimilatedPosterior, AssimError> {
    if *polls_remaining > cx.budget().poll_quota {
        return Err(AssimError::PollQuotaExceedsAmbient {
            requested: *polls_remaining,
            ambient: cx.budget().poll_quota,
        });
    }
    let plan = assimilation_work_plan(
        prior,
        observations,
        2,
        true,
        true,
        Some((regime_param, regime_lo, regime_hi)),
        cx.mode().name().len(),
    )?;
    let mut progress = WorkProgress::new_with_shared_poll_quota(cx, plan, polls_remaining);
    assimilate_colored_planned(
        prior,
        observations,
        regime_param,
        regime_lo,
        regime_hi,
        cx,
        plan,
        &mut progress,
    )
}

#[allow(clippy::too_many_arguments)]
fn assimilate_colored_planned(
    prior: &Belief,
    observations: &[Observation],
    regime_param: &str,
    regime_lo: f64,
    regime_hi: f64,
    cx: &Cx<'_>,
    plan: WorkPlan,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<AssimilatedPosterior, AssimError> {
    progress.checkpoint("initial")?;
    validate_regime(regime_param, regime_lo, regime_hi)?;
    let regime_lo = canonicalize_zero(regime_lo);
    let regime_hi = canonicalize_zero(regime_hi);
    let observations = validated_canonical_observations(observations, prior.dim(), progress)?;
    let misfit_before = misfit_canonical(prior, &observations, "misfit-before", progress)?;
    let belief = assimilate_canonical(prior, &observations, progress)?;
    let misfit_after = misfit_canonical(&belief, &observations, "misfit-after", progress)?;
    let estimator = candidate_identity(
        prior,
        &observations,
        regime_param,
        regime_lo,
        regime_hi,
        cx,
        plan,
        progress,
    )?;
    debug_assert!(fs_evidence::color_leaf_identity_reason(&estimator).is_none());

    let regime = ValidityDomain::unconstrained().with(regime_param, regime_lo, regime_hi);
    progress.checkpoint("finalize")?;

    Ok(AssimilatedPosterior {
        belief,
        color: Color::Estimated {
            estimator,
            dispersion: f64::INFINITY,
        },
        regime,
        misfit_before,
        misfit_after,
    })
}

fn checked_work_mul(left: u128, right: u128, phase: &'static str) -> Result<u128, AssimError> {
    left.checked_mul(right)
        .ok_or(AssimError::WorkPlanOverflow { phase })
}

fn checked_work_add(left: u128, right: u128, phase: &'static str) -> Result<u128, AssimError> {
    left.checked_add(right)
        .ok_or(AssimError::WorkPlanOverflow { phase })
}

fn checked_square_usize(value: usize, phase: &'static str) -> Result<usize, AssimError> {
    value
        .checked_mul(value)
        .and_then(|entries| entries.checked_mul(size_of::<f64>()))
        .ok_or(AssimError::WorkPlanOverflow { phase })
}

fn zero_matrix(
    dimension: usize,
    phase: &'static str,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<Vec<Vec<f64>>, AssimError> {
    let mut matrix = Vec::with_capacity(dimension);
    for _ in 0..dimension {
        let mut row = Vec::with_capacity(dimension);
        for _ in 0..dimension {
            row.push(0.0);
            progress.scalar(phase, 1)?;
        }
        matrix.push(row);
    }
    Ok(matrix)
}

fn clone_belief_for_update(
    prior: &Belief,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<Belief, AssimError> {
    let mut mean = Vec::with_capacity(prior.mean.len());
    for value in &prior.mean {
        mean.push(*value);
        progress.scalar("posterior-clone", 1)?;
    }
    let mut cov = Vec::with_capacity(prior.cov.len());
    for source_row in &prior.cov {
        let mut row = Vec::with_capacity(source_row.len());
        for value in source_row {
            row.push(*value);
            progress.scalar("posterior-clone", 1)?;
        }
        cov.push(row);
    }
    Ok(Belief { mean, cov })
}

fn preflight_belief_shape(mean: &[f64], cov: &[Vec<f64>]) -> Result<(), AssimError> {
    let n = mean.len();
    if n == 0 {
        return Err(AssimError::EmptyBelief);
    }
    validate_state_dimension(n)?;
    if cov.len() != n {
        return Err(AssimError::CovarianceDimensionMismatch {
            state: n,
            rows: cov.len(),
        });
    }
    for (row, values) in cov.iter().enumerate() {
        if values.len() != n {
            return Err(AssimError::CovarianceRowDimensionMismatch {
                row,
                expected: n,
                actual: values.len(),
            });
        }
    }
    let _dense_bytes = checked_square_usize(n, "belief covariance")?;
    Ok(())
}

fn belief_validation_work_plan(n: usize) -> Result<WorkPlan, AssimError> {
    let n = n as u128;
    let n2 = checked_work_mul(n, n, "belief validation")?;
    let n3 = checked_work_mul(n2, n, "belief PSD")?;
    let validation_psd = checked_work_add(
        checked_work_add(
            checked_work_mul(6, n, "belief validation")?,
            checked_work_mul(8, n2, "belief validation")?,
            "belief validation",
        )?,
        n3,
        "belief PSD",
    )?;
    WorkPlan::checked(validation_psd, 0, 0, 0, 0, 0, 0)
}

fn assimilation_work_plan(
    prior: &Belief,
    observations: &[Observation],
    misfit_passes: u128,
    include_update: bool,
    include_hash: bool,
    regime: Option<(&str, f64, f64)>,
    mode_name_len: usize,
) -> Result<WorkPlan, AssimError> {
    validate_assimilation_work(prior.dim(), observations.len())?;
    let _dense_matrix_bytes = checked_square_usize(prior.dim(), "Joseph matrices")?;
    let _order_bytes = observations
        .len()
        .checked_mul(size_of::<usize>())
        .and_then(|bytes| bytes.checked_mul(2))
        .ok_or(AssimError::WorkPlanOverflow {
            phase: "canonical ordering",
        })?;
    let n = prior.dim() as u128;
    let count = observations.len() as u128;
    let n2 = checked_work_mul(n, n, "assimilation preflight")?;
    let n3 = checked_work_mul(n2, n, "assimilation preflight")?;

    let validation_psd = checked_work_mul(
        count,
        checked_work_add(n, 4, "observation validation")?,
        "observation validation",
    )?;
    let (record_materialization, maximum_record_bytes) =
        observations
            .iter()
            .try_fold((0_u128, 0_u128), |(sum, maximum), observation| {
                let bytes = canonical_observation_size(observation)? as u128;
                Ok::<_, AssimError>((
                    checked_work_add(sum, bytes, "canonical record materialization")?,
                    maximum.max(bytes),
                ))
            })?;
    let merge_levels = if observations.len() <= 1 {
        0
    } else {
        u128::from(usize::BITS - (observations.len() - 1).leading_zeros())
    };
    let merge_slots = checked_work_mul(count, merge_levels, "canonical ordering")?;
    // Each occupied merge slot performs at most one lexicographic comparison,
    // and each comparison scans no more than the largest admitted record.
    let comparison_byte_budget =
        checked_work_mul(merge_slots, maximum_record_bytes, "canonical comparison")?;
    let canonical_ordering = checked_work_add(
        checked_work_add(
            checked_work_mul(2, count, "canonical ordering")?,
            merge_slots,
            "canonical ordering",
        )?,
        comparison_byte_budget,
        "canonical ordering",
    )?;
    let one_misfit = checked_work_mul(
        count,
        checked_work_add(checked_work_mul(2, n, "misfit")?, 4, "misfit")?,
        "misfit",
    )?;
    let misfit_work = checked_work_mul(one_misfit, misfit_passes, "misfit passes")?;
    let (joseph_update, psd_revalidation) = if include_update {
        // The posterior is cloned once before the sequential updates. Each
        // observation then performs dense products plus one computed upper
        // triangle in the Joseph covariance. Zero-initialization is explicit
        // work and the triangular term includes both its n-term dot product
        // and the final noise update.
        let clone_work = checked_work_add(n, n2, "posterior clone")?;
        let upper_triangle = checked_work_mul(
            n,
            checked_work_add(n, 1, "Joseph upper triangle")?,
            "Joseph upper triangle",
        )? / 2;
        let joseph_per_observation = checked_work_add(
            checked_work_add(
                checked_work_add(
                    n3,
                    checked_work_mul(6, n2, "Joseph update")?,
                    "Joseph update",
                )?,
                checked_work_mul(8, n, "Joseph update")?,
                "Joseph update",
            )?,
            checked_work_add(
                2,
                checked_work_mul(
                    upper_triangle,
                    checked_work_add(n, 1, "Joseph upper triangle")?,
                    "Joseph upper triangle",
                )?,
                "Joseph update",
            )?,
            "Joseph update",
        )?;
        let joseph_update = checked_work_add(
            clone_work,
            checked_work_mul(count, joseph_per_observation, "Joseph update")?,
            "Joseph update",
        )?;
        // 6n²+n³ bounds the full structural/PSD/canonicalization pass for
        // n>=2. The two-unit scalar correction makes the same bound sound at
        // n=1, where fixed per-belief work otherwise dominates the cubic term.
        let psd_per_observation = checked_work_add(
            checked_work_add(
                checked_work_mul(6, n2, "PSD revalidation")?,
                n3,
                "PSD revalidation",
            )?,
            2,
            "PSD revalidation",
        )?;
        (
            joseph_update,
            checked_work_mul(count, psd_per_observation, "PSD revalidation")?,
        )
    } else {
        (0, 0)
    };
    let hashing = if include_hash {
        let (regime_param, _, _) = regime.ok_or(AssimError::WorkPlanOverflow {
            phase: "candidate regime",
        })?;
        candidate_identity_work_size(prior, observations, regime_param, mode_name_len)?
    } else {
        0
    };

    WorkPlan::checked(
        validation_psd,
        record_materialization,
        canonical_ordering,
        misfit_work,
        joseph_update,
        psd_revalidation,
        hashing,
    )
}

fn validate_belief_parts(
    mean: &[f64],
    cov: &[Vec<f64>],
    phase: &'static str,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<(), AssimError> {
    validate_belief_structure(mean, cov, phase, progress)?;
    if !covariance_is_positive_semidefinite(cov, phase, progress)? {
        return Err(AssimError::CovarianceNotPositiveSemidefinite);
    }
    Ok(())
}

fn validate_belief_structure(
    mean: &[f64],
    cov: &[Vec<f64>],
    phase: &'static str,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<(), AssimError> {
    let n = mean.len();
    if n == 0 {
        return Err(AssimError::EmptyBelief);
    }
    validate_state_dimension(n)?;
    if cov.len() != n {
        return Err(AssimError::CovarianceDimensionMismatch {
            state: n,
            rows: cov.len(),
        });
    }
    for (index, value) in mean.iter().enumerate() {
        progress.scalar(phase, 1)?;
        if !value.is_finite() {
            return Err(AssimError::NonFiniteMean { index });
        }
    }
    for (row_index, row) in cov.iter().enumerate() {
        if row.len() != n {
            return Err(AssimError::CovarianceRowDimensionMismatch {
                row: row_index,
                expected: n,
                actual: row.len(),
            });
        }
        for (column_index, value) in row.iter().enumerate() {
            progress.scalar(phase, 1)?;
            if !value.is_finite() {
                return Err(AssimError::NonFiniteCovariance {
                    row: row_index,
                    column: column_index,
                });
            }
        }
        progress.scalar(phase, 1)?;
        if row[row_index] < 0.0 {
            return Err(AssimError::NegativeVariance { index: row_index });
        }
    }
    for (row_index, row) in cov.iter().enumerate() {
        for (column_index, column) in cov.iter().enumerate().skip(row_index + 1) {
            progress.scalar(phase, 1)?;
            if canonical_f64_bits(row[column_index]) != canonical_f64_bits(column[row_index]) {
                return Err(AssimError::NonSymmetricCovariance {
                    row: row_index,
                    column: column_index,
                });
            }
        }
    }
    Ok(())
}

fn validate_state_dimension(dim: usize) -> Result<(), AssimError> {
    if dim > MAX_DENSE_STATE_DIM {
        Err(AssimError::StateDimensionLimit {
            dim,
            max: MAX_DENSE_STATE_DIM,
        })
    } else {
        Ok(())
    }
}

fn validate_observation_count(count: usize) -> Result<(), AssimError> {
    if count == 0 {
        Err(AssimError::EmptyObservations)
    } else if count > MAX_DENSE_OBSERVATIONS {
        Err(AssimError::ObservationCountLimit {
            count,
            max: MAX_DENSE_OBSERVATIONS,
        })
    } else {
        Ok(())
    }
}

fn validate_assimilation_work(dim: usize, observation_count: usize) -> Result<(), AssimError> {
    validate_observation_count(observation_count)?;
    let dim = dim as u128;
    let requested = dim
        .checked_mul(dim)
        .and_then(|value| value.checked_mul(dim))
        .and_then(|value| value.checked_mul(observation_count as u128))
        .ok_or(AssimError::WorkPlanOverflow {
            phase: "dense assimilation envelope",
        })?;
    if requested > MAX_DENSE_UPDATE_CUBIC_WORK {
        Err(AssimError::AssimilationWorkLimit {
            requested,
            max: MAX_DENSE_UPDATE_CUBIC_WORK,
        })
    } else {
        Ok(())
    }
}

fn covariance_is_positive_semidefinite(
    cov: &[Vec<f64>],
    phase: &'static str,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<bool, AssimError> {
    // Scaling to a correlation matrix makes the Schur-complement test
    // dimensionless. Without this step, one enormous variance can hide an
    // invalid correlation involving a much smaller variance. This boundary is
    // deliberately fail-closed: unlike a solver convergence test, no negative
    // pivot is a harmless tolerance event. Ambiguous roundoff is rejected
    // rather than silently relabelled as zero curvature.
    let mut active = Vec::with_capacity(cov.len());
    for (index, row) in cov.iter().enumerate() {
        progress.scalar(phase, 1)?;
        if canonical_f64_bits(row[index]) == 0 {
            // A PSD matrix with zero variance has an exactly zero row/column.
            for entry in row {
                progress.scalar(phase, 1)?;
                if canonical_f64_bits(*entry) != 0 {
                    return Ok(false);
                }
            }
        } else {
            active.push(index);
        }
    }
    let n = active.len();
    if n == 0 {
        return Ok(true);
    }

    let mut scaled = zero_matrix(n, phase, progress)?;
    for (scaled_row, &source_row) in active.iter().enumerate() {
        scaled[scaled_row][scaled_row] = 1.0;
        progress.scalar(phase, 1)?;
        for (scaled_column, &source_column) in active.iter().enumerate().skip(scaled_row + 1) {
            progress.scalar(phase, 1)?;
            if !square_is_at_most_product(
                cov[source_row][source_column].abs(),
                cov[source_row][source_row],
                cov[source_column][source_column],
            ) {
                // This exact binary-rational comparison enforces every 2x2
                // principal minor before square roots and divisions can round
                // an invalid correlation back onto the unit boundary.
                return Ok(false);
            }
            let row_scale = cov[source_row][source_row].sqrt();
            let column_scale = cov[source_column][source_column].sqrt();
            let (first_divisor, second_divisor) = if row_scale >= column_scale {
                (row_scale, column_scale)
            } else {
                (column_scale, row_scale)
            };
            let correlation = cov[source_row][source_column] / first_divisor / second_divisor;
            if !correlation.is_finite() || correlation.abs() > 1.0 {
                return Ok(false);
            }
            scaled[scaled_row][scaled_column] = correlation;
            scaled[scaled_column][scaled_row] = correlation;
        }
    }

    // Symmetric diagonal pivoting avoids dividing by a small Schur pivot when
    // a better-conditioned one remains. The transformation is a sequence of
    // congruences, so a negative diagonal in any Schur complement is direct
    // evidence of negative curvature.
    for pivot_index in 0..n {
        let mut selected = pivot_index;
        for candidate in (pivot_index + 1)..n {
            progress.scalar(phase, 1)?;
            if scaled[candidate][candidate] > scaled[selected][selected] {
                selected = candidate;
            }
        }
        if selected != pivot_index {
            scaled.swap(selected, pivot_index);
            for row in &mut scaled {
                row.swap(selected, pivot_index);
            }
        }

        let pivot = scaled[pivot_index][pivot_index];
        if !pivot.is_finite() || pivot < 0.0 {
            return Ok(false);
        }
        if canonical_f64_bits(pivot) == 0 {
            // A PSD matrix with a zero diagonal has an exactly zero row and
            // column. Accept exact singular structure, but never manufacture
            // it by tolerance-clamping a negative pivot.
            for entry in &scaled[pivot_index][(pivot_index + 1)..] {
                progress.scalar(phase, 1)?;
                if canonical_f64_bits(*entry) != 0 {
                    return Ok(false);
                }
            }
            continue;
        }

        let mut pivot_column = Vec::with_capacity(n);
        for row in &scaled {
            pivot_column.push(row[pivot_index]);
            progress.scalar(phase, 1)?;
        }
        for row in (pivot_index + 1)..n {
            let multiplier = pivot_column[row] / pivot;
            progress.scalar(phase, 1)?;
            if !multiplier.is_finite() {
                return Ok(false);
            }
            for (column, column_pivot) in pivot_column.iter().enumerate().skip(row) {
                let updated = (-multiplier).mul_add(*column_pivot, scaled[column][row]);
                progress.scalar(phase, 1)?;
                if !updated.is_finite() || (row == column && updated < 0.0) {
                    return Ok(false);
                }
                scaled[column][row] = updated;
                scaled[row][column] = updated;
            }
        }
    }
    Ok(true)
}

fn square_is_at_most_product(value: f64, left: f64, right: f64) -> bool {
    let square = binary_product(value, value);
    let diagonal_product = binary_product(left, right);
    compare_binary_products(square, diagonal_product) != core::cmp::Ordering::Greater
}

fn binary_product(left: f64, right: f64) -> (u128, i32) {
    let (left_significand, left_exponent) = binary_significand_and_exponent(left);
    let (right_significand, right_exponent) = binary_significand_and_exponent(right);
    (
        u128::from(left_significand) * u128::from(right_significand),
        left_exponent + right_exponent,
    )
}

fn binary_significand_and_exponent(value: f64) -> (u64, i32) {
    debug_assert!(value.is_finite() && value >= 0.0);
    let bits = value.to_bits();
    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
    let exponent_bits = ((bits >> 52) & 0x7ff) as i32; // Masked to eleven bits.
    let fraction = bits & ((1_u64 << 52) - 1);
    if exponent_bits == 0 {
        (fraction, -1074)
    } else {
        ((1_u64 << 52) | fraction, exponent_bits - 1023 - 52)
    }
}

fn compare_binary_products(left: (u128, i32), right: (u128, i32)) -> core::cmp::Ordering {
    let ((left_significand, left_exponent), (right_significand, right_exponent)) = (left, right);
    match (left_significand == 0, right_significand == 0) {
        (true, true) => return core::cmp::Ordering::Equal,
        (true, false) => return core::cmp::Ordering::Less,
        (false, true) => return core::cmp::Ordering::Greater,
        (false, false) => {}
    }

    let left_top_bit = i64::from(left_significand.ilog2()) + i64::from(left_exponent);
    let right_top_bit = i64::from(right_significand.ilog2()) + i64::from(right_exponent);
    match left_top_bit.cmp(&right_top_bit) {
        core::cmp::Ordering::Equal => {
            if left_exponent >= right_exponent {
                (left_significand << (left_exponent - right_exponent).unsigned_abs())
                    .cmp(&right_significand)
            } else {
                left_significand
                    .cmp(&(right_significand << (right_exponent - left_exponent).unsigned_abs()))
            }
        }
        ordering => ordering,
    }
}

fn validate_noise(noise_var: f64) -> Result<(), AssimError> {
    if !noise_var.is_finite() {
        Err(AssimError::NonFiniteNoise)
    } else if noise_var <= 0.0 {
        Err(AssimError::NonPositiveNoise)
    } else {
        Ok(())
    }
}

fn validate_leaf_identity(field: &'static str, identity: &str) -> Result<(), AssimError> {
    if identity.trim().is_empty() {
        return match field {
            "instrument" => Err(AssimError::EmptyInstrument),
            "regime_param" => Err(AssimError::EmptyRegime),
            _ => Err(AssimError::InvalidIdentity {
                field,
                reason: "blank",
            }),
        };
    }
    if let Some(reason) = fs_evidence::color_leaf_identity_reason(identity) {
        return Err(AssimError::InvalidIdentity { field, reason });
    }
    Ok(())
}

fn validate_regime(regime_param: &str, lo: f64, hi: f64) -> Result<(), AssimError> {
    validate_leaf_identity("regime_param", regime_param)?;
    if !lo.is_finite() || !hi.is_finite() {
        return Err(AssimError::NonFiniteRegimeBounds);
    }
    if lo > hi {
        return Err(AssimError::InvertedRegimeBounds);
    }
    Ok(())
}

struct CanonicalObservation<'a> {
    bytes: Vec<u8>,
    observation: &'a Observation,
}

struct CanonicalObservations<'a> {
    records: Vec<CanonicalObservation<'a>>,
    order: Vec<usize>,
}

impl CanonicalObservations<'_> {
    fn len(&self) -> usize {
        self.order.len()
    }

    fn ordered_indices(&self) -> impl Iterator<Item = usize> + '_ {
        self.order.iter().copied()
    }
}

fn validate_observation_for_dim(
    observation: &Observation,
    state_dim: usize,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<(), AssimError> {
    if observation.operator.is_empty() {
        return Err(AssimError::EmptyObservationOperator);
    }
    validate_state_dimension(observation.operator.len())?;
    let mut any_nonzero = false;
    for (index, coefficient) in observation.operator.iter().enumerate() {
        progress.scalar("observation-validation", 1)?;
        if !coefficient.is_finite() {
            return Err(AssimError::NonFiniteObservationOperator { index });
        }
        any_nonzero |= canonical_f64_bits(*coefficient) != 0;
    }
    if !any_nonzero {
        return Err(AssimError::ZeroObservationOperator);
    }
    if !observation.value.is_finite() {
        return Err(AssimError::NonFiniteObservationValue);
    }
    validate_noise(observation.noise_var)?;
    validate_leaf_identity("instrument", &observation.instrument)?;
    if observation.operator.len() != state_dim {
        return Err(AssimError::DimMismatch {
            state: state_dim,
            operator: observation.operator.len(),
        });
    }
    progress.scalar("observation-validation", 4)?;
    Ok(())
}

fn atom_encoded_size(label_len: usize, value_len: usize) -> Result<usize, AssimError> {
    size_of::<u128>()
        .checked_mul(2)
        .and_then(|framing| framing.checked_add(label_len))
        .and_then(|with_label| with_label.checked_add(value_len))
        .ok_or(AssimError::WorkPlanOverflow {
            phase: "canonical atom",
        })
}

fn canonical_observation_size(observation: &Observation) -> Result<usize, AssimError> {
    let mut size = atom_encoded_size(b"operator-length".len(), size_of::<u128>())?;
    let coefficient_atom = atom_encoded_size(b"operator-coefficient".len(), size_of::<u64>())?;
    let value_atom = atom_encoded_size(b"value".len(), size_of::<u64>())?;
    let noise_atom = atom_encoded_size(b"noise-variance".len(), size_of::<u64>())?;
    let instrument_atom = atom_encoded_size(b"instrument".len(), observation.instrument.len())?;
    size = size
        .checked_add(
            coefficient_atom
                .checked_mul(observation.operator.len())
                .ok_or(AssimError::WorkPlanOverflow {
                    phase: "canonical operator",
                })?,
        )
        .and_then(|value| value.checked_add(value_atom))
        .and_then(|value| value.checked_add(noise_atom))
        .and_then(|value| value.checked_add(instrument_atom))
        .ok_or(AssimError::WorkPlanOverflow {
            phase: "canonical observation",
        })?;
    Ok(size)
}

fn push_record_atom(
    buffer: &mut Vec<u8>,
    label: &[u8],
    value: &[u8],
    progress: &mut WorkProgress<'_, '_>,
) -> Result<(), AssimError> {
    append_record_bytes(buffer, &usize_bytes(label.len()), progress)?;
    append_record_bytes(buffer, label, progress)?;
    append_record_bytes(buffer, &usize_bytes(value.len()), progress)?;
    append_record_bytes(buffer, value, progress)
}

fn append_record_bytes(
    buffer: &mut Vec<u8>,
    bytes: &[u8],
    progress: &mut WorkProgress<'_, '_>,
) -> Result<(), AssimError> {
    for chunk in bytes.chunks(RECORD_POLL_STRIDE as usize) {
        buffer.extend_from_slice(chunk);
        progress.records("canonical-materialization", chunk.len() as u128)?;
    }
    Ok(())
}

fn canonical_observation_bytes(
    observation: &Observation,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<Vec<u8>, AssimError> {
    let mut record = Vec::with_capacity(canonical_observation_size(observation)?);
    push_record_atom(
        &mut record,
        b"operator-length",
        &usize_bytes(observation.operator.len()),
        progress,
    )?;
    for coefficient in &observation.operator {
        push_record_atom(
            &mut record,
            b"operator-coefficient",
            &canonical_f64_bits(*coefficient).to_le_bytes(),
            progress,
        )?;
    }
    push_record_atom(
        &mut record,
        b"value",
        &canonical_f64_bits(observation.value).to_le_bytes(),
        progress,
    )?;
    push_record_atom(
        &mut record,
        b"noise-variance",
        &canonical_f64_bits(observation.noise_var).to_le_bytes(),
        progress,
    )?;
    push_record_atom(
        &mut record,
        b"instrument",
        observation.instrument.as_bytes(),
        progress,
    )?;
    debug_assert_eq!(canonical_observation_size(observation), Ok(record.len()));
    Ok(record)
}

fn compare_canonical_records(
    left: &[u8],
    right: &[u8],
    progress: &mut WorkProgress<'_, '_>,
) -> Result<core::cmp::Ordering, AssimError> {
    let common_bytes = left.len().min(right.len());
    let stride = CANONICAL_COMPARE_BYTE_POLL_STRIDE as usize;
    let mut offset = 0;
    while offset < common_bytes {
        // Honor the accumulated boundary across record comparisons, rather
        // than allowing two short tails to create a larger unpolled window.
        let bytes_until_poll = stride - progress.comparison_bytes_since_poll as usize;
        let end = offset.saturating_add(bytes_until_poll).min(common_bytes);
        let ordering = left[offset..end].cmp(&right[offset..end]);
        progress.comparison_bytes("canonical-compare", (end - offset) as u128)?;
        if ordering != core::cmp::Ordering::Equal {
            return Ok(ordering);
        }
        offset = end;
    }
    Ok(left.len().cmp(&right.len()))
}

fn validated_canonical_observations<'a>(
    observations: &'a [Observation],
    state_dim: usize,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<CanonicalObservations<'a>, AssimError> {
    validate_observation_count(observations.len())?;
    for observation in observations {
        validate_observation_for_dim(observation, state_dim, progress)?;
    }

    let mut records = Vec::with_capacity(observations.len());
    for observation in observations {
        records.push(CanonicalObservation {
            bytes: canonical_observation_bytes(observation, progress)?,
            observation,
        });
    }

    let mut order = Vec::with_capacity(records.len());
    for index in 0..records.len() {
        order.push(index);
        progress.records("canonical-ordering", 1)?;
    }
    let mut scratch = Vec::with_capacity(records.len());
    for _ in 0..records.len() {
        scratch.push(0);
        progress.records("canonical-ordering", 1)?;
    }
    let mut width = 1_usize;
    while width < order.len() {
        let mut start = 0_usize;
        while start < order.len() {
            let middle = start.saturating_add(width).min(order.len());
            let end = middle.saturating_add(width).min(order.len());
            let (mut left, mut right, mut output) = (start, middle, start);
            while left < middle && right < end {
                let left_index = order[left];
                let right_index = order[right];
                if compare_canonical_records(
                    &records[left_index].bytes,
                    &records[right_index].bytes,
                    progress,
                )? != core::cmp::Ordering::Greater
                {
                    scratch[output] = left_index;
                    left += 1;
                } else {
                    scratch[output] = right_index;
                    right += 1;
                }
                output += 1;
                progress.records("canonical-ordering", 1)?;
            }
            while left < middle {
                scratch[output] = order[left];
                left += 1;
                output += 1;
                progress.records("canonical-ordering", 1)?;
            }
            while right < end {
                scratch[output] = order[right];
                right += 1;
                output += 1;
                progress.records("canonical-ordering", 1)?;
            }
            start = end;
        }
        core::mem::swap(&mut order, &mut scratch);
        width = width.checked_mul(2).ok_or(AssimError::WorkPlanOverflow {
            phase: "canonical ordering",
        })?;
    }
    Ok(CanonicalObservations { records, order })
}

fn assimilate_checked(
    prior: &Belief,
    obs: &Observation,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<Belief, AssimError> {
    let n = prior.dim();
    let h = &obs.operator;

    let mut ph = Vec::with_capacity(n);
    for row in &prior.cov {
        ph.push(checked_dot(
            row,
            h,
            "covariance-times-operator",
            "joseph-update",
            progress,
        )?);
    }
    let innovation_variance =
        checked_dot(h, &ph, "innovation variance", "joseph-update", progress)? + obs.noise_var;
    progress.scalar("joseph-update", 1)?;
    if !innovation_variance.is_finite() {
        return Err(AssimError::NonFiniteComputation {
            stage: "innovation variance",
        });
    }
    if innovation_variance <= 0.0 {
        return Err(AssimError::SingularInnovation);
    }
    let mut gain = Vec::with_capacity(n);
    for entry in &ph {
        let gain_entry = entry / innovation_variance;
        progress.scalar("joseph-update", 1)?;
        if !gain_entry.is_finite() {
            return Err(AssimError::NonFiniteComputation {
                stage: "Kalman gain",
            });
        }
        gain.push(gain_entry);
    }

    let predicted = checked_dot(
        h,
        &prior.mean,
        "observation prediction",
        "joseph-update",
        progress,
    )?;
    let innovation = obs.value - predicted;
    progress.scalar("joseph-update", 1)?;
    if !innovation.is_finite() {
        return Err(AssimError::NonFiniteComputation {
            stage: "observation innovation",
        });
    }

    let mut mean = Vec::with_capacity(n);
    for (prior_mean, gain_entry) in prior.mean.iter().zip(&gain) {
        let updated = prior_mean + gain_entry * innovation;
        progress.scalar("joseph-update", 2)?;
        if !updated.is_finite() {
            return Err(AssimError::NonFiniteComputation {
                stage: "posterior mean",
            });
        }
        mean.push(updated);
    }

    let cov = joseph_covariance(prior, h, obs.noise_var, &gain, progress)?;
    Belief::from_covariance_preserving_update(mean, cov, progress)
}

fn joseph_covariance(
    prior: &Belief,
    observation_operator: &[f64],
    noise_variance: f64,
    gain: &[f64],
    progress: &mut WorkProgress<'_, '_>,
) -> Result<Vec<Vec<f64>>, AssimError> {
    let n = prior.dim();
    // Joseph form, P' = (I-KH)P(I-KH)^T + KRK^T, retains both PSD terms
    // instead of relying on a cancellation-prone rank-one subtraction. The
    // final matrix is mirrored from one computed triangle for exact symmetry
    // and then passes through the full public Belief validator.
    let mut transform = zero_matrix(n, "joseph-update", progress)?;
    for (row, transform_row) in transform.iter_mut().enumerate() {
        for (column, entry) in transform_row.iter_mut().enumerate() {
            let identity = if row == column { 1.0 } else { 0.0 };
            *entry = (-gain[row]).mul_add(observation_operator[column], identity);
            progress.scalar("joseph-update", 1)?;
            if !entry.is_finite() {
                return Err(AssimError::NonFiniteComputation {
                    stage: "Joseph transform",
                });
            }
        }
    }

    let mut transformed_prior = zero_matrix(n, "joseph-update", progress)?;
    for (row, transformed_row) in transformed_prior.iter_mut().enumerate() {
        for (column, transformed_entry) in transformed_row.iter_mut().enumerate() {
            let mut entry = 0.0;
            for (transform_entry, prior_row) in transform[row].iter().zip(&prior.cov) {
                entry = transform_entry.mul_add(prior_row[column], entry);
                progress.scalar("joseph-update", 1)?;
                if !entry.is_finite() {
                    return Err(AssimError::NonFiniteComputation {
                        stage: "Joseph left product",
                    });
                }
            }
            *transformed_entry = entry;
        }
    }

    let noise_scale = noise_variance.sqrt();
    let mut noise_factor = Vec::with_capacity(n);
    for gain_entry in gain {
        let factor = gain_entry * noise_scale;
        progress.scalar("joseph-update", 1)?;
        if !factor.is_finite() {
            return Err(AssimError::NonFiniteComputation {
                stage: "Joseph noise factor",
            });
        }
        noise_factor.push(factor);
    }

    let mut cov = zero_matrix(n, "joseph-update", progress)?;
    for row in 0..n {
        for column in row..n {
            let propagated = checked_dot_fma(
                &transformed_prior[row],
                &transform[column],
                "Joseph propagated covariance",
                "joseph-update",
                progress,
            )?;
            let updated = noise_factor[row].mul_add(noise_factor[column], propagated);
            progress.scalar("joseph-update", 1)?;
            if !updated.is_finite() {
                return Err(AssimError::NonFiniteComputation {
                    stage: "posterior covariance",
                });
            }
            cov[row][column] = updated;
            cov[column][row] = updated;
        }
    }
    Ok(cov)
}

fn checked_dot_fma(
    a: &[f64],
    b: &[f64],
    stage: &'static str,
    phase: &'static str,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<f64, AssimError> {
    debug_assert_eq!(a.len(), b.len());
    let mut total = 0.0;
    for (left, right) in a.iter().zip(b) {
        total = left.mul_add(*right, total);
        progress.scalar(phase, 1)?;
        if !total.is_finite() {
            return Err(AssimError::NonFiniteComputation { stage });
        }
    }
    Ok(total)
}

fn checked_dot(
    a: &[f64],
    b: &[f64],
    stage: &'static str,
    phase: &'static str,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<f64, AssimError> {
    debug_assert_eq!(a.len(), b.len());
    let mut total = 0.0;
    for (left, right) in a.iter().zip(b) {
        let product = left * right;
        progress.scalar(phase, 1)?;
        if !product.is_finite() {
            return Err(AssimError::NonFiniteComputation { stage });
        }
        total += product;
        progress.scalar(phase, 1)?;
        if !total.is_finite() {
            return Err(AssimError::NonFiniteComputation { stage });
        }
    }
    Ok(total)
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn candidate_identity(
    prior: &Belief,
    observations: &CanonicalObservations<'_>,
    regime_param: &str,
    regime_lo: f64,
    regime_hi: f64,
    cx: &Cx<'_>,
    plan: WorkPlan,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<String, AssimError> {
    let mut hasher = fs_blake3::Blake3::new();
    hash_atom(
        &mut hasher,
        b"identity-domain",
        CANDIDATE_ID_DOMAIN.as_bytes(),
        progress,
    )?;
    hash_atom(
        &mut hasher,
        b"state-dimension",
        &usize_bytes(prior.dim()),
        progress,
    )?;
    for value in &prior.mean {
        hash_atom(
            &mut hasher,
            b"prior-mean",
            &canonical_f64_bits(*value).to_le_bytes(),
            progress,
        )?;
    }
    for row in &prior.cov {
        for value in row {
            hash_atom(
                &mut hasher,
                b"prior-covariance",
                &canonical_f64_bits(*value).to_le_bytes(),
                progress,
            )?;
        }
    }

    for record_index in observations.ordered_indices() {
        hash_atom(
            &mut hasher,
            b"observation",
            &observations.records[record_index].bytes,
            progress,
        )?;
    }
    hash_atom(
        &mut hasher,
        b"regime-axis",
        regime_param.as_bytes(),
        progress,
    )?;
    hash_atom(
        &mut hasher,
        b"regime-lo",
        &canonical_f64_bits(regime_lo).to_le_bytes(),
        progress,
    )?;
    hash_atom(
        &mut hasher,
        b"regime-hi",
        &canonical_f64_bits(regime_hi).to_le_bytes(),
        progress,
    )?;

    let budget = cx.budget();
    hash_atom(
        &mut hasher,
        b"cx-mode",
        cx.mode().name().as_bytes(),
        progress,
    )?;
    let stream = cx.stream_key();
    for (label, value) in [
        (b"cx-stream-seed".as_slice(), stream.seed),
        (b"cx-stream-kernel".as_slice(), stream.kernel_id),
        (b"cx-stream-tile".as_slice(), stream.tile),
        (b"cx-stream-iteration".as_slice(), stream.iteration),
    ] {
        hash_atom(&mut hasher, label, &value.to_le_bytes(), progress)?;
    }
    hash_atom(
        &mut hasher,
        b"budget-deadline-present",
        &[u8::from(budget.deadline.is_some())],
        progress,
    )?;
    hash_atom(
        &mut hasher,
        b"budget-deadline-ns",
        &match budget.deadline {
            Some(deadline) => deadline.as_nanos(),
            None => 0,
        }
        .to_le_bytes(),
        progress,
    )?;
    hash_atom(
        &mut hasher,
        b"budget-poll-quota",
        &budget.poll_quota.to_le_bytes(),
        progress,
    )?;
    let effective_poll_quota = progress.initial_poll_quota;
    hash_atom(
        &mut hasher,
        b"effective-poll-quota",
        &effective_poll_quota.to_le_bytes(),
        progress,
    )?;
    hash_atom(
        &mut hasher,
        b"budget-cost-present",
        &[u8::from(budget.cost_quota.is_some())],
        progress,
    )?;
    hash_atom(
        &mut hasher,
        b"budget-cost-quota",
        &budget.cost_quota.unwrap_or(0).to_le_bytes(),
        progress,
    )?;
    hash_atom(
        &mut hasher,
        b"budget-priority",
        &[budget.priority],
        progress,
    )?;

    for (label, units) in [
        (b"plan-validation-psd".as_slice(), plan.validation_psd),
        (
            b"plan-record-materialization".as_slice(),
            plan.record_materialization,
        ),
        (
            b"plan-canonical-ordering".as_slice(),
            plan.canonical_ordering,
        ),
        (b"plan-misfit-passes".as_slice(), plan.misfit_passes),
        (b"plan-joseph-update".as_slice(), plan.joseph_update),
        (b"plan-psd-revalidation".as_slice(), plan.psd_revalidation),
        (b"plan-hashing".as_slice(), plan.hashing),
        (b"plan-total".as_slice(), plan.total),
    ] {
        hash_atom(&mut hasher, label, &units.to_le_bytes(), progress)?;
    }
    hash_atom(
        &mut hasher,
        b"poll-policy",
        POLL_POLICY_ID.as_bytes(),
        progress,
    )?;
    hash_atom(
        &mut hasher,
        b"scalar-poll-stride",
        &SCALAR_POLL_STRIDE.to_le_bytes(),
        progress,
    )?;
    hash_atom(
        &mut hasher,
        b"record-poll-stride",
        &RECORD_POLL_STRIDE.to_le_bytes(),
        progress,
    )?;
    hash_atom(
        &mut hasher,
        b"canonical-compare-byte-poll-stride",
        &CANONICAL_COMPARE_BYTE_POLL_STRIDE.to_le_bytes(),
        progress,
    )?;
    hash_atom(
        &mut hasher,
        b"hash-byte-poll-stride",
        &(HASH_BYTE_POLL_STRIDE as u128).to_le_bytes(),
        progress,
    )?;

    Ok(format!("{CANDIDATE_ID_PREFIX}{}", hasher.finalize()))
}

fn hash_update(
    hasher: &mut fs_blake3::Blake3,
    bytes: &[u8],
    progress: &mut WorkProgress<'_, '_>,
) -> Result<(), AssimError> {
    for chunk in bytes.chunks(HASH_BYTE_POLL_STRIDE) {
        hasher.update(chunk);
        progress.hash_bytes("candidate-hash", chunk.len() as u128)?;
    }
    Ok(())
}

fn hash_atom(
    hasher: &mut fs_blake3::Blake3,
    label: &[u8],
    value: &[u8],
    progress: &mut WorkProgress<'_, '_>,
) -> Result<(), AssimError> {
    hash_update(hasher, &usize_bytes(label.len()), progress)?;
    hash_update(hasher, label, progress)?;
    hash_update(hasher, &usize_bytes(value.len()), progress)?;
    hash_update(hasher, value, progress)
}

fn candidate_identity_work_size(
    prior: &Belief,
    observations: &[Observation],
    regime_param: &str,
    mode_name_len: usize,
) -> Result<u128, AssimError> {
    let mut total = 0_usize;
    add_identity_atoms(&mut total, b"identity-domain", CANDIDATE_ID_DOMAIN.len(), 1)?;
    add_identity_atoms(&mut total, b"state-dimension", size_of::<u128>(), 1)?;
    add_identity_atoms(
        &mut total,
        b"prior-mean",
        size_of::<u64>(),
        prior.mean.len(),
    )?;
    let covariance_entries =
        prior
            .dim()
            .checked_mul(prior.dim())
            .ok_or(AssimError::WorkPlanOverflow {
                phase: "candidate prior covariance",
            })?;
    add_identity_atoms(
        &mut total,
        b"prior-covariance",
        size_of::<u64>(),
        covariance_entries,
    )?;
    for observation in observations {
        add_identity_atoms(
            &mut total,
            b"observation",
            canonical_observation_size(observation)?,
            1,
        )?;
    }
    add_identity_atoms(&mut total, b"regime-axis", regime_param.len(), 1)?;
    add_identity_atoms(&mut total, b"regime-lo", size_of::<u64>(), 1)?;
    add_identity_atoms(&mut total, b"regime-hi", size_of::<u64>(), 1)?;
    add_identity_atoms(&mut total, b"cx-mode", mode_name_len, 1)?;
    for label in [
        b"cx-stream-seed".as_slice(),
        b"cx-stream-kernel".as_slice(),
        b"cx-stream-tile".as_slice(),
        b"cx-stream-iteration".as_slice(),
    ] {
        add_identity_atoms(&mut total, label, size_of::<u64>(), 1)?;
    }
    add_identity_atoms(&mut total, b"budget-deadline-present", size_of::<u8>(), 1)?;
    add_identity_atoms(&mut total, b"budget-deadline-ns", size_of::<u64>(), 1)?;
    add_identity_atoms(&mut total, b"budget-poll-quota", size_of::<u32>(), 1)?;
    add_identity_atoms(&mut total, b"effective-poll-quota", size_of::<u32>(), 1)?;
    add_identity_atoms(&mut total, b"budget-cost-present", size_of::<u8>(), 1)?;
    add_identity_atoms(&mut total, b"budget-cost-quota", size_of::<u64>(), 1)?;
    add_identity_atoms(&mut total, b"budget-priority", size_of::<u8>(), 1)?;
    for label in [
        b"plan-validation-psd".as_slice(),
        b"plan-record-materialization".as_slice(),
        b"plan-canonical-ordering".as_slice(),
        b"plan-misfit-passes".as_slice(),
        b"plan-joseph-update".as_slice(),
        b"plan-psd-revalidation".as_slice(),
        b"plan-hashing".as_slice(),
        b"plan-total".as_slice(),
    ] {
        add_identity_atoms(&mut total, label, size_of::<u128>(), 1)?;
    }
    add_identity_atoms(&mut total, b"poll-policy", POLL_POLICY_ID.len(), 1)?;
    add_identity_atoms(&mut total, b"scalar-poll-stride", size_of::<u128>(), 1)?;
    add_identity_atoms(&mut total, b"record-poll-stride", size_of::<u128>(), 1)?;
    add_identity_atoms(
        &mut total,
        b"canonical-compare-byte-poll-stride",
        size_of::<u128>(),
        1,
    )?;
    add_identity_atoms(&mut total, b"hash-byte-poll-stride", size_of::<u128>(), 1)?;
    Ok(total as u128)
}

fn add_identity_atoms(
    total: &mut usize,
    label: &[u8],
    value_len: usize,
    count: usize,
) -> Result<(), AssimError> {
    let repeated = atom_encoded_size(label.len(), value_len)?
        .checked_mul(count)
        .ok_or(AssimError::WorkPlanOverflow {
            phase: "candidate identity",
        })?;
    *total = total
        .checked_add(repeated)
        .ok_or(AssimError::WorkPlanOverflow {
            phase: "candidate identity",
        })?;
    Ok(())
}

fn usize_bytes(value: usize) -> [u8; 16] {
    (value as u128).to_le_bytes()
}

fn canonical_f64_bits(value: f64) -> u64 {
    const SIGN_BIT: u64 = 1_u64 << 63;
    match value.to_bits() {
        SIGN_BIT => 0,
        bits => bits,
    }
}

fn canonicalize_zero(value: f64) -> f64 {
    f64::from_bits(canonical_f64_bits(value))
}

fn canonicalize_belief_zeros(
    mean: &mut [f64],
    cov: &mut [Vec<f64>],
    progress: &mut WorkProgress<'_, '_>,
) -> Result<(), AssimError> {
    for value in mean {
        *value = canonicalize_zero(*value);
        progress.scalar("belief-canonicalization", 1)?;
    }
    for value in cov.iter_mut().flatten() {
        *value = canonicalize_zero(*value);
        progress.scalar("belief-canonicalization", 1)?;
    }
    Ok(())
}
