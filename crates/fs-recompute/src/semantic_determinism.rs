//! Semantic computation keys, tolerance roles, determinism classes, and legacy
//! migration (bead `frankensim-mkfvu.1`).
//!
//! # Core Contract and Architecture
//!
//! A computation's identity and its determinism contract depend on whether
//! `required_tolerance` legitimately alters the executed algorithm:
//!
//! 1. **`ToleranceRole::StoppingCriterion` / `ToleranceRole::InputParameter`**:
//!    The required tolerance directly determines iteration count, grid sizing,
//!    or polynomial order. Different tolerances produce different computations
//!    and different artifacts by design. The tolerance is bound into the
//!    [`ComputationKey`].
//! 2. **`ToleranceRole::QueryThreshold`**:
//!    The computation is fixed and pure (e.g. evaluating a closed-form expression
//!    or fixed-order quadrature); tolerance is only a post-hoc query criterion
//!    for memoized cache hit testing. The tolerance is NOT part of the
//!    [`ComputationKey`]; [`OutputObservation::achieved_error`] satisfies or
//!    fails the query slack.
//! 3. **`DeterminismClass`**:
//!    Explicitly models whether bit-identical results are guaranteed across runs
//!    ([`DeterminismClass::ExactDeterministic`] and
//!    [`DeterminismClass::ToleranceDependentDeterministic`]) or whether the op
//!    operates under relaxed heuristics ([`DeterminismClass::Nondeterministic`]).
//!    Nondeterministic operations cannot mint verified or certified authority.

use core::fmt;
use std::collections::BTreeMap;

use fs_ledger::{ContentHash, hash_bytes};

use crate::PinReason;

/// Exact domain string for semantic computation keys.
pub const COMPUTATION_KEY_DOMAIN: &str = "org.frankensim.fs-recompute.computation-key.v1";

/// Exact domain string for output observations.
pub const OUTPUT_OBSERVATION_DOMAIN: &str = "org.frankensim.fs-recompute.output-observation.v1";

/// Semantic version of the computation key encoding.
pub const COMPUTATION_KEY_VERSION: u32 = 1;

/// Determinism guarantee class for an operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DeterminismClass {
    /// Bit-identical output across runs, thread counts, and completions on the same ISA.
    ExactDeterministic,
    /// Bit-identical output for fixed `(inputs, params, policy, seed, tolerance)`.
    ToleranceDependentDeterministic,
    /// Relaxed or stochastic execution mode with no bit-level reproducibility contract.
    Nondeterministic,
}

impl DeterminismClass {
    /// Whether this class provides certified deterministic reproducibility.
    #[must_use]
    pub const fn is_deterministic(self) -> bool {
        match self {
            Self::ExactDeterministic | Self::ToleranceDependentDeterministic => true,
            Self::Nondeterministic => false,
        }
    }

    /// Machine-readable label.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::ExactDeterministic => "exact-deterministic",
            Self::ToleranceDependentDeterministic => "tolerance-dependent-deterministic",
            Self::Nondeterministic => "nondeterministic",
        }
    }
}

impl fmt::Display for DeterminismClass {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

/// The role that tolerance plays in a given operation family.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ToleranceRole {
    /// Tolerance is an explicit algorithmic input (e.g. adaptive mesh edge length).
    InputParameter,
    /// Tolerance is an iterative stopping criterion (e.g. Krylov solver relative residual).
    StoppingCriterion,
    /// Tolerance is a post-hoc query threshold for cache hits, not an input to the kernel.
    QueryThreshold,
    /// Operation does not use or accept tolerances.
    None,
}

impl ToleranceRole {
    /// Whether the tolerance alters the kernel's execution and must be part of the key.
    #[must_use]
    pub const fn affects_computation(self) -> bool {
        match self {
            Self::InputParameter | Self::StoppingCriterion => true,
            Self::QueryThreshold | Self::None => false,
        }
    }

    /// Machine-readable label.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::InputParameter => "input-parameter",
            Self::StoppingCriterion => "stopping-criterion",
            Self::QueryThreshold => "query-threshold",
            Self::None => "none",
        }
    }
}

/// Supported operation family catalog.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum OperationFamily {
    /// Exact geometric reductions and discrete topology evaluations.
    DiscreteGeometry,
    /// Adaptive geometric representation conversion and tessellation.
    AdaptiveConversion,
    /// Linear Krylov and sparse system solves.
    LinearSolve,
    /// Nonlinear root searches and optimization linesearches.
    NonlinearSolve,
    /// Deterministic pseudo-random simulation (counter-based Philox RNG).
    DeterministicStochastic,
    /// Unmonitored or exploratory fast heuristics.
    FastHeuristic,
}

impl OperationFamily {
    /// Default determinism class for this family.
    #[must_use]
    pub const fn default_determinism_class(self) -> DeterminismClass {
        match self {
            Self::DiscreteGeometry => DeterminismClass::ExactDeterministic,
            Self::AdaptiveConversion | Self::LinearSolve | Self::NonlinearSolve => {
                DeterminismClass::ToleranceDependentDeterministic
            }
            Self::DeterministicStochastic => DeterminismClass::ExactDeterministic,
            Self::FastHeuristic => DeterminismClass::Nondeterministic,
        }
    }

    /// Default tolerance role for this family.
    #[must_use]
    pub const fn default_tolerance_role(self) -> ToleranceRole {
        match self {
            Self::DiscreteGeometry => ToleranceRole::None,
            Self::AdaptiveConversion => ToleranceRole::InputParameter,
            Self::LinearSolve | Self::NonlinearSolve => ToleranceRole::StoppingCriterion,
            Self::DeterministicStochastic => ToleranceRole::QueryThreshold,
            Self::FastHeuristic => ToleranceRole::None,
        }
    }
}

/// Result of evaluating an operation against the determinism contract.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DeterminismDisposition {
    /// Execution satisfied the determinism contract (bit-identical repeat or valid first run).
    Satisfied,
    /// Bit mismatch detected between identical computation keys.
    Violation {
        /// Content hash recorded previously.
        expected_artifact: ContentHash,
        /// Content hash produced by current run.
        actual_artifact: ContentHash,
        /// Actionable diagnosis of probable cause.
        diagnosis: &'static str,
    },
    /// Policy mismatch (e.g. nondeterministic run attempted on deterministic contract).
    PolicyMismatch {
        /// Expected determinism class.
        expected: DeterminismClass,
        /// Actual execution policy mode.
        actual: String,
    },
    /// Nondeterministic execution without authority claims.
    UncheckedNondeterministic,
}

/// Execution policy controlling algorithm behavior and tolerance usage.
#[derive(Debug, Clone, PartialEq)]
pub struct ExecutionPolicy {
    /// Determinism guarantee class.
    pub determinism_class: DeterminismClass,
    /// Role of required tolerance.
    pub tolerance_role: ToleranceRole,
    /// Required tolerance value, if applicable. Must be finite and positive if set.
    pub required_tolerance: Option<f64>,
    /// RNG seed for counter-based streams.
    pub rng_seed: u64,
    /// Code version hash of the executing kernel.
    pub code_version_hash: ContentHash,
    /// Max iterations or budget ceiling, if specified.
    pub max_iterations: Option<u64>,
}

impl ExecutionPolicy {
    /// Construct a strictly validated execution policy.
    ///
    /// # Errors
    /// Refuses non-finite or non-positive tolerances when required by the tolerance role.
    pub fn try_new(
        determinism_class: DeterminismClass,
        tolerance_role: ToleranceRole,
        required_tolerance: Option<f64>,
        rng_seed: u64,
        code_version_hash: ContentHash,
        max_iterations: Option<u64>,
    ) -> Result<Self, SemanticKeyError> {
        if let Some(tol) = required_tolerance {
            if !tol.is_finite() || tol <= 0.0 {
                return Err(SemanticKeyError::InvalidTolerance {
                    value: tol,
                    reason: "required tolerance must be finite and strictly positive",
                });
            }
        } else if tolerance_role.affects_computation() {
            return Err(SemanticKeyError::MissingTolerance {
                role: tolerance_role,
            });
        }

        Ok(Self {
            determinism_class,
            tolerance_role,
            required_tolerance,
            rng_seed,
            code_version_hash,
            max_iterations,
        })
    }

    /// Construct a default exact deterministic policy with no tolerance dependency.
    #[must_use]
    pub fn exact_deterministic(code_version_hash: ContentHash, rng_seed: u64) -> Self {
        Self {
            determinism_class: DeterminismClass::ExactDeterministic,
            tolerance_role: ToleranceRole::None,
            required_tolerance: None,
            rng_seed,
            code_version_hash,
            max_iterations: None,
        }
    }
}

/// Canonical semantic key identifying a unit of computation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ComputationKey {
    /// Stable operation name/verb.
    pub op_id: String,
    /// Ordered content hashes of all upstream inputs.
    pub input_hashes: Vec<ContentHash>,
    /// Canonical parameter map (sorted lexicographically by key).
    pub params: BTreeMap<String, String>,
    /// Execution policy governing this computation.
    pub policy_determinism_class: DeterminismClass,
    /// Tolerance role for this computation.
    pub policy_tolerance_role: ToleranceRole,
    /// Required tolerance bits if the tolerance affects computation; 0 if not.
    pub effective_tolerance_bits: u64,
    /// RNG seed.
    pub rng_seed: u64,
    /// Code version hash.
    pub code_version_hash: ContentHash,
    /// Max iterations if specified.
    pub max_iterations: Option<u64>,
}

impl ComputationKey {
    /// Construct and validate a computation key.
    ///
    /// # Errors
    /// Refuses empty op_id, invalid parameters, or malformed policy values.
    pub fn try_new(
        op_id: impl Into<String>,
        input_hashes: Vec<ContentHash>,
        params: BTreeMap<String, String>,
        policy: &ExecutionPolicy,
    ) -> Result<Self, SemanticKeyError> {
        let op_id = op_id.into();
        if op_id.trim().is_empty() {
            return Err(SemanticKeyError::EmptyOpId);
        }

        let effective_tolerance_bits = if policy.tolerance_role.affects_computation() {
            let tol = policy
                .required_tolerance
                .ok_or(SemanticKeyError::MissingTolerance {
                    role: policy.tolerance_role,
                })?;
            tol.to_bits()
        } else {
            0
        };

        Ok(Self {
            op_id,
            input_hashes,
            params,
            policy_determinism_class: policy.determinism_class,
            policy_tolerance_role: policy.tolerance_role,
            effective_tolerance_bits,
            rng_seed: policy.rng_seed,
            code_version_hash: policy.code_version_hash,
            max_iterations: policy.max_iterations,
        })
    }

    /// Compute the canonical content hash over this computation key.
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        let mut bytes = Vec::new();
        push_string(&mut bytes, COMPUTATION_KEY_DOMAIN);
        push_u64(&mut bytes, u64::from(COMPUTATION_KEY_VERSION));
        push_string(&mut bytes, &self.op_id);

        push_u64(&mut bytes, self.input_hashes.len() as u64);
        for hash in &self.input_hashes {
            push_bytes(&mut bytes, hash.as_bytes());
        }

        push_u64(&mut bytes, self.params.len() as u64);
        for (k, v) in &self.params {
            push_string(&mut bytes, k);
            push_string(&mut bytes, v);
        }

        push_string(&mut bytes, self.policy_determinism_class.as_str());
        push_string(&mut bytes, self.policy_tolerance_role.as_str());
        push_u64(&mut bytes, self.effective_tolerance_bits);
        push_u64(&mut bytes, self.rng_seed);
        push_bytes(&mut bytes, self.code_version_hash.as_bytes());
        push_u64(&mut bytes, self.max_iterations.unwrap_or(0));

        hash_bytes(&bytes)
    }
}

/// Observed output metrics and artifact address from a completed computation.
#[derive(Debug, Clone, PartialEq)]
pub struct OutputObservation {
    /// Content-addressed identity of produced artifact.
    pub artifact_hash: ContentHash,
    /// Measured error or residual achieved, if applicable.
    pub achieved_error: Option<f64>,
    /// Monotonic wall time consumed in seconds.
    pub wall_time_s: Option<f64>,
    /// Peak memory footprint in bytes.
    pub peak_memory_bytes: Option<u64>,
}

impl OutputObservation {
    /// Construct a new output observation.
    ///
    /// # Errors
    /// Refuses non-finite achieved error or negative timing/memory.
    pub fn try_new(
        artifact_hash: ContentHash,
        achieved_error: Option<f64>,
        wall_time_s: Option<f64>,
        peak_memory_bytes: Option<u64>,
    ) -> Result<Self, SemanticKeyError> {
        if let Some(err) = achieved_error {
            if !err.is_finite() || err < 0.0 {
                return Err(SemanticKeyError::InvalidAchievedError {
                    value: err,
                    reason: "achieved error must be finite and non-negative",
                });
            }
        }
        if let Some(wall) = wall_time_s {
            if !wall.is_finite() || wall < 0.0 {
                return Err(SemanticKeyError::InvalidTiming { value: wall });
            }
        }
        Ok(Self {
            artifact_hash,
            achieved_error,
            wall_time_s,
            peak_memory_bytes,
        })
    }
}

/// Structured errors for semantic computation keys and determinism policies.
#[derive(Debug, Clone, PartialEq)]
pub enum SemanticKeyError {
    /// Empty operation identifier.
    EmptyOpId,
    /// Invalid tolerance specification.
    InvalidTolerance {
        /// Provided value.
        value: f64,
        /// Reason for refusal.
        reason: &'static str,
    },
    /// Required tolerance is missing for an operation whose tolerance role affects computation.
    MissingTolerance {
        /// Declared role.
        role: ToleranceRole,
    },
    /// Invalid achieved error measurement.
    InvalidAchievedError {
        /// Provided value.
        value: f64,
        /// Reason for refusal.
        reason: &'static str,
    },
    /// Invalid execution timing.
    InvalidTiming {
        /// Non-finite or negative timing.
        value: f64,
    },
}

impl fmt::Display for SemanticKeyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyOpId => write!(f, "operation identifier must not be empty"),
            Self::InvalidTolerance { value, reason } => {
                write!(f, "invalid required tolerance {value}: {reason}")
            }
            Self::MissingTolerance { role } => {
                write!(
                    f,
                    "tolerance role `{}` requires an explicit tolerance value",
                    role.as_str()
                )
            }
            Self::InvalidAchievedError { value, reason } => {
                write!(f, "invalid achieved error {value}: {reason}")
            }
            Self::InvalidTiming { value } => {
                write!(
                    f,
                    "invalid wall timing {value}: must be finite and non-negative"
                )
            }
        }
    }
}

impl std::error::Error for SemanticKeyError {}

fn push_u64(bytes: &mut Vec<u8>, v: u64) {
    bytes.extend_from_slice(&v.to_le_bytes());
}

fn push_bytes(bytes: &mut Vec<u8>, slice: &[u8]) {
    push_u64(bytes, slice.len() as u64);
    bytes.extend_from_slice(slice);
}

fn push_string(bytes: &mut Vec<u8>, s: &str) {
    push_bytes(bytes, s.as_bytes());
}

/// A stored semantic computation node.
#[derive(Debug, Clone, PartialEq)]
pub struct StoredSemanticComputation {
    /// Canonical computation key.
    pub key: ComputationKey,
    /// Output observation recorded from execution.
    pub observation: OutputObservation,
    /// Pinned reasons (empty = evictable).
    pub pins: Vec<PinReason>,
    /// Insertion order sequence.
    pub seq: u64,
}

/// Migrate a legacy 7-field [`crate::NodeRecord`] to a versioned [`ComputationKey`]
/// and [`OutputObservation`].
///
/// Conservatively assumes [`ToleranceRole::StoppingCriterion`] when `required_tolerance > 0`,
/// preventing false-positive determinism trip-wire trips across legitimate discretization levels.
#[must_use]
pub fn migrate_legacy_node_record(
    record: &crate::NodeRecord,
    artifact_hash: ContentHash,
) -> (ComputationKey, OutputObservation) {
    let mut params = BTreeMap::new();
    for (k, v) in &record.params {
        let val_str = match v {
            crate::ParamValue::F64(bits) => format!("f64:{bits:016X}"),
            crate::ParamValue::Int(i) => format!("int:{i}"),
            crate::ParamValue::Str(s) => format!("str:{s}"),
        };
        params.insert(k.clone(), val_str);
    }

    let (tolerance_role, required_tolerance) = if record.required_tolerance > 0.0 {
        (
            ToleranceRole::StoppingCriterion,
            Some(record.required_tolerance),
        )
    } else {
        (ToleranceRole::None, None)
    };

    let policy = ExecutionPolicy {
        determinism_class: DeterminismClass::ToleranceDependentDeterministic,
        tolerance_role,
        required_tolerance,
        rng_seed: record.rng_seed,
        code_version_hash: record.code_version_hash,
        max_iterations: None,
    };

    let comp_key =
        ComputationKey::try_new(&record.op_id, record.input_hashes.clone(), params, &policy)
            .expect("valid legacy record migration");

    let obs = OutputObservation {
        artifact_hash,
        achieved_error: Some(record.achieved_error),
        wall_time_s: None,
        peak_memory_bytes: None,
    };

    (comp_key, obs)
}
