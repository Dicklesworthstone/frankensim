//! fs-diffreal-e2e — the differentiation & reality end-to-end suite (plan
//! addendum, Proposal 11 / Layer-3 conformance). Layer: L6.
//!
//! A runnable battery for selected Layer-3 integration fixtures: a shared
//! transpose engine over a local scalar path, a synthetic as-built loop, and
//! tolerance allocation. It records whether those fixed fixtures fail safe.
//! Four stages emit structured log events (returned as data, never printed):
//!
//! 1. **Differentiation** — the production `fs-adjoint` tape/VJP path agrees
//!    with independent dual-number and two-step finite-difference oracles, and
//!    a path with a MISSING VJP (a forced remesh) raises a structured error
//!    that BLOCKS the gradient — never a silent zero.
//! 2. **As-built loop** — register a scanned fixture (error carried forward),
//!    compute an estimated as-built δ carrying calibration provenance,
//!    LOCALIZE a seeded defect, and run registration-free point-sensor
//!    assimilation that reduces the model-data misfit ([`fs_asbuilt`],
//!    [`fs_assimilate`]).
//! 3. **Tolerance allocation** — a GD&T report consumes sealed sensitivities,
//!    tightens the high-sensitivity feature, loosens the low one, and reports
//!    only whether caller-supplied band samples agree with the linearization.
//!    It deliberately makes no probability claim ([`fs_toleralloc`]).
//! 4. **(Gated) spacetime** — the temporal-complex capability exists in
//!    `fs-time`, but its coupled end-to-end fixture is not integrated and
//!    activated in this battery; it is reported as gated, not silently passed.
//!
//! [`run_battery`] runs all four under an explicit [`Cx`] and returns a
//! structured [`DiffRealReport`] only after every cancellation-aware stage has
//! finalized.

use fs_ad::dual::{Dual64, gradient as dual_gradient};
use fs_adjoint::transpose::{Tape, TransposeError, Vjp, VjpRegistry, fd_falsifier};
mod invocation;

use fs_asbuilt::{Fiducial, Point2, as_built_diff_budgeted, register_budgeted};
use fs_assimilate::{AssimError, Belief, Observation, assimilate_colored_budgeted, point_sensor};
use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::Color;
use fs_exec::{Budget, Cx, ExecMode, StreamKey};
use fs_toleralloc::{
    Action, ColorRank, Feature, allocate, gdt_report, robustness_check, variance_budget,
};
use std::collections::BTreeSet;
use std::sync::Arc;

/// Stable name of the differentiation stage.
pub const DIFFERENTIATION_STAGE: &str = "differentiation";
/// Stable name of the as-built/assimilation stage.
pub const AS_BUILT_STAGE: &str = "as-built-loop";
/// Stable name of the tolerance-allocation stage.
pub const TOLERANCE_STAGE: &str = "tolerance-allocation";
/// Stable name of the spacetime-integration stage.
pub const SPACETIME_STAGE: &str = "spacetime-gated";

/// Versioned fixture identity expected for the differentiation stage.
pub const DIFFERENTIATION_EVIDENCE_IDENTITY: &str = "fs-diffreal-e2e/differentiation-fixture/v2";
/// Versioned fixture identity expected for the as-built/assimilation stage.
pub const AS_BUILT_EVIDENCE_IDENTITY: &str = "fs-diffreal-e2e/as-built-fixture/v2";
/// Versioned fixture identity expected for the tolerance-allocation stage.
pub const TOLERANCE_EVIDENCE_IDENTITY: &str = "fs-diffreal-e2e/tolerance-allocation-fixture/v3";
/// Versioned fixture identity expected for the spacetime-integration stage.
pub const SPACETIME_EVIDENCE_IDENTITY: &str = "fs-diffreal-e2e/spacetime-integration-gate/v1";

/// Version of the production-path sensitivity sealing policy.
pub const SENSITIVITY_POLICY_VERSION: &str = "fs-diffreal-e2e/sensitivity-policy/v1";
/// Version of the fixed local affine/square/identity VJP registry semantics.
pub const DIFFERENTIATION_REGISTRY_POLICY: &str = "fs-diffreal-e2e/production-vjp-registry/v1";
/// Version of the canonical stage-receipt schema and verification policy.
pub const STAGE_RECEIPT_POLICY_VERSION: u32 = 2;
/// Version of the ordered report-root schema and verification policy.
pub const REPORT_RECEIPT_POLICY_VERSION: u32 = 2;
/// Versioned battery readiness policy bound into every report root.
pub const REPORT_PROMOTION_POLICY: &str = "fs-diffreal-e2e/promotion-policy/v3";

/// The production differentiation fixture's operator path.
pub const PRODUCTION_DIFFERENTIATION_PATH: [&str; 3] = ["sdf", "spline", "solve"];

const SENSITIVITY_IDENTITY_DOMAIN: &str = "frankensim.fs-diffreal-e2e.sensitivity.v1";
const FIXTURE_INPUT_IDENTITY_DOMAIN: &str = "frankensim.fs-diffreal-e2e.fixture-inputs.v1";
const STAGE_RECEIPT_IDENTITY_DOMAIN: &str = "frankensim.fs-diffreal-e2e.stage-receipt.v2";
const STAGE_RESULT_IDENTITY_DOMAIN: &str = "frankensim.fs-diffreal-e2e.stage-results.v1";
const REPORT_RECEIPT_IDENTITY_DOMAIN: &str = "frankensim.fs-diffreal-e2e.report-receipt.v2";
const PROMOTION_POLICY_FINGERPRINT_DOMAIN: &str =
    "frankensim.fs-diffreal-e2e.promotion-policy-fingerprint.v1";
const PROMOTION_VERIFICATION_SUBJECT_DOMAIN: &str =
    "frankensim.fs-diffreal-e2e.promotion-verification-subject.v1";
const PROMOTION_VERIFICATION_PURPOSE: &str = "frankensim.diffreal-promotion.v1";
const ASSIMILATION_POLL_POLICY: &str = "fixed-stride:v3";
const MAX_DIFFERENTIATION_OPS: usize = 16;
const MAX_OP_NAME_BYTES: usize = 64;
const DIFFERENTIATION_WORK_UNITS: u64 = 12;
const DIFFERENTIATION_STAGE_WORK_UNITS: u64 = 24;
const TOLERANCE_WORK_UNITS: u64 = 32;
const SPACETIME_WORK_UNITS: u64 = 1;

const DIFFERENTIATION_FIXTURE_INPUT: f64 = 1.5;
const MISSING_VJP_FIXTURE_PATH: [&str; 3] = ["sdf", "remesh", "solve"];
const AS_BUILT_DESIGN_POINTS: [(f64, f64); 3] = [(0.0, 0.0), (2.0, 0.0), (0.0, 2.0)];
const AS_BUILT_ROTATION_RADIANS: f64 = 0.3;
const AS_BUILT_TRANSLATION: (f64, f64) = (4.0, 1.0);
const AS_BUILT_DEFECT_INDEX: usize = 1;
const AS_BUILT_DEFECT_X: f64 = 0.3;
const AS_BUILT_DESIGN_TOLERANCE: f64 = 0.5;
const AS_BUILT_MEASUREMENT_NOISE: f64 = 0.02;
const AS_BUILT_CALIBRATION_CANDIDATE: &str = "cmm-cal-2026";
const AS_BUILT_PRIOR_MEAN: [f64; 2] = [20.0, 20.0];
const AS_BUILT_PRIOR_DIAGONAL_COVARIANCE: [f64; 2] = [9.0, 9.0];
const AS_BUILT_OBSERVATIONS: [(usize, f64, f64, &str); 2] = [
    (0, 24.0, 0.25, "thermocouple-1"),
    (1, 18.5, 0.25, "thermocouple-2"),
];
const AS_BUILT_ASSIMILATION_PARAMETER: &str = "Re";
const AS_BUILT_ASSIMILATION_BOUNDS: (f64, f64) = (1.0e5, 3.0e5);
const TOLERANCE_SENSITIVITY_INPUTS: [f64; 2] = [1.0, -0.475];
const TOLERANCE_EXTREME_QOIS: [f64; 3] = [0.9, -0.8, 0.5];
const TOLERANCE_PERFORMANCE_TOLERANCE: f64 = 1.0;
const TOLERANCE_TARGET_PROBABILITY: f64 = 0.99;
const TOLERANCE_SIGMA_MULTIPLIER: f64 = 3.0;
const TOLERANCE_NOMINAL_QOI: f64 = 0.0;
const TOLERANCE_ROBUSTNESS_MARGIN: f64 = 0.2;
const TOLERANCE_FEATURE_NAMES: [&str; 2] = ["critical", "slack"];
const TOLERANCE_FEATURE_COST_COEFFICIENT: f64 = 1.0;
const TOLERANCE_FEATURE_BASELINE: f64 = 0.5;

/// Typed refusal from the production differentiation and independent-oracle
/// path. Floating-point payloads are retained as exact bits so errors compare
/// and replay deterministically even for NaN inputs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DifferentiationError {
    /// Cancellation was observed at a bounded stage boundary.
    Cancelled,
    /// The ambient cost quota cannot admit the fixed work plan.
    WorkBudgetExceeded {
        /// Work units required by the operation.
        required: u64,
        /// Work units made available by the ambient context.
        available: u64,
    },
    /// An empty path has no production operator semantics to verify.
    EmptyPath,
    /// The independent fixture oracles are defined only for the canonical
    /// affine→square→identity path.
    OraclePathMismatch {
        /// Expected operator count.
        expected_len: usize,
        /// Supplied operator count.
        observed_len: usize,
        /// First differing position, or the shared prefix length when only the
        /// lengths differ.
        first_mismatch: usize,
    },
    /// A caller supplied more operator nodes than the bounded tape admits.
    PathTooLong {
        /// Maximum admitted operator count.
        limit: usize,
        /// Supplied operator count.
        observed: usize,
    },
    /// One operator name exceeded the bounded registry/token limit.
    OpNameTooLong {
        /// Position in the operator path or registry insertion.
        index: usize,
        /// Maximum admitted bytes.
        limit: usize,
        /// Supplied bytes.
        observed: usize,
    },
    /// An empty operator name cannot identify a VJP or forward semantic.
    EmptyOpName {
        /// Position in the operator path or registry insertion.
        index: usize,
    },
    /// Coverage is checked before numeric evaluation, so a structural seam
    /// cannot be hidden behind a NaN or overflow refusal.
    MissingVjp {
        /// First uncovered operator in forward path order.
        op: String,
    },
    /// The local production fixture has no forward semantic for this name.
    UnsupportedOperator {
        /// Registered but unsupported forward operator.
        op: String,
    },
    /// The input is NaN or infinite.
    NonFiniteInput {
        /// Exact rejected bits.
        bits: u64,
    },
    /// A finite input produced a non-finite primal at one operator.
    NonFinitePrimal {
        /// Operator that produced the value.
        op: String,
        /// Exact rejected bits.
        bits: u64,
    },
    /// The shared production transpose refused a declared seam.
    Transpose(TransposeError),
    /// The shared transpose returned no cotangent for the input leaf.
    MissingLeafGradient,
    /// The scalar fixture received a non-scalar leaf cotangent.
    InvalidGradientShape {
        /// Returned cotangent width.
        observed: usize,
    },
    /// The production VJP sweep returned NaN or infinity.
    NonFiniteGradient {
        /// Exact rejected bits.
        bits: u64,
    },
    /// Production, dual, and two-step finite-difference evidence disagreed.
    OracleDisagreement {
        /// Production reverse-mode gradient bits.
        production_bits: u64,
        /// Independent dual gradient bits.
        dual_bits: u64,
        /// Fine central-difference gradient bits.
        fd_fine_bits: u64,
        /// Conditioning-aware accepted-difference bits.
        tolerance_bits: u64,
    },
    /// A unit scale must be finite and strictly positive.
    InvalidInputScale {
        /// Exact rejected bits.
        bits: u64,
    },
    /// Applying an admitted unit scale overflowed the sensitivity.
    NonFiniteRescaledGradient {
        /// Exact unrepresentable result bits.
        bits: u64,
    },
    /// A sealed sensitivity no longer matches its fixed-schema identity.
    SensitivityIntegrityMismatch {
        /// Identity carried by the rejected receipt.
        identity: ContentHash,
    },
}

impl DifferentiationError {
    fn is_runtime_refusal(&self) -> bool {
        matches!(self, Self::Cancelled | Self::WorkBudgetExceeded { .. })
    }
}

impl core::fmt::Display for DifferentiationError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Cancelled => {
                formatter.write_str("differentiation cancelled at a bounded checkpoint")
            }
            Self::WorkBudgetExceeded {
                required,
                available,
            } => write!(
                formatter,
                "differentiation requires {required} work units but the context admits {available}"
            ),
            Self::EmptyPath => {
                formatter.write_str("differentiation path contains no production operators")
            }
            Self::OraclePathMismatch {
                expected_len,
                observed_len,
                first_mismatch,
            } => write!(
                formatter,
                "independent oracles require the canonical {expected_len}-operator path; observed {observed_len} operators with the first mismatch at index {first_mismatch}"
            ),
            Self::PathTooLong { limit, observed } => write!(
                formatter,
                "differentiation path has {observed} operators; the admitted limit is {limit}"
            ),
            Self::OpNameTooLong {
                index,
                limit,
                observed,
            } => write!(
                formatter,
                "operator {index} has {observed} name bytes; the admitted limit is {limit}"
            ),
            Self::EmptyOpName { index } => {
                write!(formatter, "operator {index} has an empty registry name")
            }
            Self::MissingVjp { op } => write!(
                formatter,
                "op '{op}' has no registered VJP: the gradient is BLOCKED (never silent-zero)"
            ),
            Self::UnsupportedOperator { op } => write!(
                formatter,
                "op '{op}' has a VJP registration but no admitted forward semantic in this fixture"
            ),
            Self::NonFiniteInput { bits } => {
                write!(
                    formatter,
                    "differentiation input is non-finite (bits=0x{bits:016x})"
                )
            }
            Self::NonFinitePrimal { op, bits } => write!(
                formatter,
                "op '{op}' produced a non-finite primal (bits=0x{bits:016x})"
            ),
            Self::Transpose(error) => write!(formatter, "production transpose refused: {error}"),
            Self::MissingLeafGradient => {
                formatter.write_str("production transpose returned no input-leaf gradient")
            }
            Self::InvalidGradientShape { observed } => write!(
                formatter,
                "production transpose returned {observed} input cotangents; expected exactly one"
            ),
            Self::NonFiniteGradient { bits } => write!(
                formatter,
                "production transpose returned a non-finite gradient (bits=0x{bits:016x})"
            ),
            Self::OracleDisagreement {
                production_bits,
                dual_bits,
                fd_fine_bits,
                tolerance_bits,
            } => write!(
                formatter,
                "production/dual/FD gradients disagree: production=0x{production_bits:016x}, dual=0x{dual_bits:016x}, fd=0x{fd_fine_bits:016x}, tolerance=0x{tolerance_bits:016x}"
            ),
            Self::InvalidInputScale { bits } => write!(
                formatter,
                "input-unit scale must be finite and positive (bits=0x{bits:016x})"
            ),
            Self::NonFiniteRescaledGradient { bits } => write!(
                formatter,
                "input-unit rescaling produced a non-finite gradient (bits=0x{bits:016x})"
            ),
            Self::SensitivityIntegrityMismatch { identity } => write!(
                formatter,
                "sealed sensitivity no longer matches its fixed-schema identity {identity}"
            ),
        }
    }
}

impl std::error::Error for DifferentiationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Transpose(error) => Some(error),
            _ => None,
        }
    }
}

/// Safe local façade over the shared production VJP registry. The companion
/// name set makes missing-VJP admission possible before any numeric work.
#[derive(Debug, Default)]
pub struct DifferentiationRegistry {
    inner: VjpRegistry,
    registered: BTreeSet<String>,
}

impl DifferentiationRegistry {
    /// Construct an empty registry.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Register a differentiable scalar operator.
    ///
    /// # Errors
    /// Refuses empty or overlong names before cloning them.
    pub fn register<V>(&mut self, op: &str, vjp: V) -> Result<(), DifferentiationError>
    where
        V: Vjp + 'static,
    {
        admit_op_name(0, op)?;
        self.inner.register(op, Arc::new(vjp));
        self.registered.insert(op.to_string());
        Ok(())
    }

    /// Declare an operator non-differentiable with an explicit consequence.
    ///
    /// # Errors
    /// Refuses empty or overlong names before cloning them.
    pub fn declare_non_differentiable(
        &mut self,
        op: &str,
        reason: &str,
        consequence: &str,
    ) -> Result<(), DifferentiationError> {
        admit_op_name(0, op)?;
        self.inner
            .declare_non_differentiable(op, reason, consequence);
        self.registered.insert(op.to_string());
        Ok(())
    }
}

/// Exact deterministic result from the production tape/VJP sweep.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PathDerivative {
    value_bits: u64,
    gradient_bits: u64,
}

impl PathDerivative {
    /// Production primal value.
    #[must_use]
    pub fn value(self) -> f64 {
        f64::from_bits(self.value_bits)
    }

    /// Production reverse-mode gradient.
    #[must_use]
    pub fn gradient(self) -> f64 {
        f64::from_bits(self.gradient_bits)
    }
}

/// Opaque sensitivity evidence minted only after the production reverse sweep
/// agrees with independent dual and two-step finite-difference oracles.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SealedSensitivity {
    ops: Vec<String>,
    input_bits: u64,
    value_bits: u64,
    production_gradient_bits: u64,
    dual_gradient_bits: u64,
    fd_coarse_bits: u64,
    fd_fine_bits: u64,
    fd_tolerance_bits: u64,
    identity: ContentHash,
}

impl SealedSensitivity {
    /// Production primal value.
    #[must_use]
    pub fn value(&self) -> f64 {
        f64::from_bits(self.value_bits)
    }

    /// Independently checked production gradient.
    #[must_use]
    pub fn gradient(&self) -> f64 {
        f64::from_bits(self.production_gradient_bits)
    }

    /// Independent dual-number gradient.
    #[must_use]
    pub fn dual_gradient(&self) -> f64 {
        f64::from_bits(self.dual_gradient_bits)
    }

    /// Coarse central-difference result.
    #[must_use]
    pub fn fd_coarse(&self) -> f64 {
        f64::from_bits(self.fd_coarse_bits)
    }

    /// Fine central-difference result.
    #[must_use]
    pub fn fd_fine(&self) -> f64 {
        f64::from_bits(self.fd_fine_bits)
    }

    /// Conditioning-aware FD acceptance tolerance.
    #[must_use]
    pub fn fd_tolerance(&self) -> f64 {
        f64::from_bits(self.fd_tolerance_bits)
    }

    /// Content identity binding path, input, all oracle values, and policy.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Recompute the fixed-schema content identity.
    #[must_use]
    pub fn verifies_integrity(&self) -> bool {
        self.identity
            == sensitivity_identity(
                &self.ops,
                [
                    self.input_bits,
                    self.value_bits,
                    self.production_gradient_bits,
                    self.dual_gradient_bits,
                    self.fd_coarse_bits,
                    self.fd_fine_bits,
                    self.fd_tolerance_bits,
                ],
            )
    }

    /// Express the gradient in caller units when one caller input unit equals
    /// `canonical_per_input_unit` canonical units.
    ///
    /// # Errors
    /// Refuses a non-finite/non-positive scale or an unrepresentable result.
    pub fn gradient_in_input_units(
        &self,
        canonical_per_input_unit: f64,
    ) -> Result<f64, DifferentiationError> {
        if !canonical_per_input_unit.is_finite() || canonical_per_input_unit <= 0.0 {
            return Err(DifferentiationError::InvalidInputScale {
                bits: canonical_per_input_unit.to_bits(),
            });
        }
        let gradient = self.gradient() * canonical_per_input_unit;
        if !gradient.is_finite() {
            return Err(DifferentiationError::NonFiniteRescaledGradient {
                bits: gradient.to_bits(),
            });
        }
        Ok(gradient)
    }
}

/// Typed event emitted by a battery stage. Crate-authored fixed-fixture events
/// are deterministic and floating-point fields use exact bits. Public callers
/// may construct diagnostic gate/refusal details of arbitrary length, but such
/// data carries no battery-report authority.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StageEvent {
    /// A production gradient passed both independent oracles.
    GradientVerified {
        /// Content identity of the sealed sensitivity.
        receipt: ContentHash,
        /// Input bits.
        input_bits: u64,
        /// Production primal bits.
        value_bits: u64,
        /// Production reverse-mode gradient bits.
        production_bits: u64,
        /// Independent dual gradient bits.
        dual_bits: u64,
        /// Coarse central-difference bits.
        fd_coarse_bits: u64,
        /// Fine central-difference bits.
        fd_fine_bits: u64,
        /// FD acceptance tolerance bits.
        tolerance_bits: u64,
    },
    /// A production or independent-oracle check refused the gradient.
    DifferentiationRejected {
        /// Typed cause.
        error: DifferentiationError,
    },
    /// The deliberate missing-VJP falsifier disposition.
    MissingVjpProbe {
        /// Missing operator name.
        op: String,
        /// Whether the production path blocked it.
        blocked: bool,
    },
    /// Registration assertion.
    Registration {
        /// RMS residual bits.
        residual_bits: u64,
        /// Fixed-fixture bound result.
        within_tolerance: bool,
    },
    /// As-built defect localization without color upgrade.
    AsBuiltDelta {
        /// Maximum deviation bits.
        max_deviation_bits: u64,
        /// Deterministic maximum index.
        defect_index: Option<usize>,
        /// Whether the candidate remained Estimated.
        estimated: bool,
    },
    /// Before/after assimilation misfit.
    Assimilation {
        /// Before-misfit bits.
        before_bits: u64,
        /// After-misfit bits.
        after_bits: u64,
        /// Whether the checked misfit decreased.
        reduced: bool,
    },
    /// Direction chosen for the two tolerance features.
    ToleranceActions {
        /// Critical-feature action.
        critical: Option<Action>,
        /// Slack-feature action.
        slack: Option<Action>,
    },
    /// GD&T loosening justification disposition.
    GdtJustification {
        /// Number of loosened features.
        loosened: usize,
        /// Whether every loosening used sealed Verified sensitivity.
        all_verified: bool,
    },
    /// Sampled linearization check; deliberately not a probability statement.
    SampledLinearization {
        /// Number of caller-provided band samples.
        samples: usize,
        /// Whether those samples stayed inside the linearized bound.
        confirmed: bool,
        /// Linearized standard-deviation bits.
        linearized_std_bits: u64,
        /// Always false in crate-authored tolerance events: samples do not
        /// prove probability. Caller-authored diagnostics have no authority.
        probability_claimed: bool,
    },
    /// Structured deliberate capability gate.
    Gate {
        /// Stable code.
        code: &'static str,
        /// Deterministic detail.
        detail: String,
    },
    /// Structured inability to evaluate a scientific assertion.
    Refusal {
        /// Stable code.
        code: &'static str,
        /// Deterministic detail.
        detail: String,
    },
}

impl StageEvent {
    /// Stable event-kind code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::GradientVerified { .. } => "gradient-verified",
            Self::DifferentiationRejected { .. } => "differentiation-rejected",
            Self::MissingVjpProbe { .. } => "missing-vjp-probe",
            Self::Registration { .. } => "registration",
            Self::AsBuiltDelta { .. } => "as-built-delta",
            Self::Assimilation { .. } => "assimilation",
            Self::ToleranceActions { .. } => "tolerance-actions",
            Self::GdtJustification { .. } => "gdt-justification",
            Self::SampledLinearization { .. } => "sampled-linearization",
            Self::Gate { .. } => "gate",
            Self::Refusal { .. } => "refusal",
        }
    }

    fn is_well_formed(&self) -> bool {
        match self {
            Self::MissingVjpProbe { op, .. } => !op.trim().is_empty(),
            Self::Gate { code, detail } | Self::Refusal { code, detail } => {
                !code.trim().is_empty() && !detail.trim().is_empty()
            }
            _ => true,
        }
    }
}

impl core::fmt::Display for StageEvent {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::GradientVerified {
                receipt,
                input_bits,
                value_bits,
                production_bits,
                dual_bits,
                fd_coarse_bits,
                fd_fine_bits,
                tolerance_bits,
            } => write!(
                formatter,
                "event={} receipt={} input=0x{input_bits:016x} value=0x{value_bits:016x} production=0x{production_bits:016x} dual=0x{dual_bits:016x} fd_coarse=0x{fd_coarse_bits:016x} fd_fine=0x{fd_fine_bits:016x} tolerance=0x{tolerance_bits:016x}",
                self.code(),
                receipt
            ),
            Self::DifferentiationRejected { error } => {
                write!(formatter, "event={} error={error}", self.code())
            }
            Self::MissingVjpProbe { op, blocked } => {
                write!(formatter, "event={} op={op} blocked={blocked}", self.code())
            }
            Self::Registration {
                residual_bits,
                within_tolerance,
            } => write!(
                formatter,
                "event={} residual=0x{residual_bits:016x} within_tolerance={within_tolerance}",
                self.code()
            ),
            Self::AsBuiltDelta {
                max_deviation_bits,
                defect_index,
                estimated,
            } => write!(
                formatter,
                "event={} max_deviation=0x{max_deviation_bits:016x} defect_index={defect_index:?} estimated={estimated}",
                self.code()
            ),
            Self::Assimilation {
                before_bits,
                after_bits,
                reduced,
            } => write!(
                formatter,
                "event={} before=0x{before_bits:016x} after=0x{after_bits:016x} reduced={reduced}",
                self.code()
            ),
            Self::ToleranceActions { critical, slack } => write!(
                formatter,
                "event={} critical={critical:?} slack={slack:?}",
                self.code()
            ),
            Self::GdtJustification {
                loosened,
                all_verified,
            } => write!(
                formatter,
                "event={} loosened={loosened} all_verified={all_verified}",
                self.code()
            ),
            Self::SampledLinearization {
                samples,
                confirmed,
                linearized_std_bits,
                probability_claimed,
            } => write!(
                formatter,
                "event={} samples={samples} confirmed={confirmed} linearized_std=0x{linearized_std_bits:016x} probability_claimed={probability_claimed}",
                self.code()
            ),
            Self::Gate { code, detail } | Self::Refusal { code, detail } => {
                write!(
                    formatter,
                    "event={} code={code} detail={detail}",
                    self.code()
                )
            }
        }
    }
}

/// Whether a stage is load-bearing for this battery's promotion decision.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StageRequirement {
    /// The report is incomplete until the stage has actually run.
    Required,
    /// The stage is diagnostic and does not block the required-stage decision.
    Optional,
}

impl core::fmt::Display for StageRequirement {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Required => formatter.write_str("required"),
            Self::Optional => formatter.write_str("optional"),
        }
    }
}

/// Stable machine code plus deterministic human-readable detail for a stage
/// that did not pass.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StageReason {
    /// Stable reason code for ledgers and programmatic diagnostics.
    pub code: &'static str,
    /// Human-readable detail. This is diagnostic data, never printed here.
    pub detail: String,
}

impl StageReason {
    /// Construct a structured reason.
    #[must_use]
    pub fn new(code: &'static str, detail: impl Into<String>) -> Self {
        Self {
            code,
            detail: detail.into(),
        }
    }

    fn is_well_formed(&self) -> bool {
        !self.code.trim().is_empty() && !self.detail.trim().is_empty()
    }
}

impl core::fmt::Display for StageReason {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(formatter, "[{}]: {}", self.code, self.detail)
    }
}

/// Scientific disposition of one stage.
///
/// `Failed` means the stage ran and an assertion was false. `Gated` and
/// `Refused` mean the assertion was not validly evaluated, so neither can
/// satisfy report completeness.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StageStatus {
    /// Every load-bearing assertion ran and passed.
    Passed,
    /// The stage ran, but at least one load-bearing assertion was false.
    Failed(StageReason),
    /// The capability or integration is deliberately unavailable.
    Gated(StageReason),
    /// The stage declined to evaluate because an admissibility condition,
    /// budget, or cancellation condition prevented a trustworthy result.
    Refused(StageReason),
}

impl StageStatus {
    /// Stable lowercase status code for deterministic records.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::Passed => "passed",
            Self::Failed(_) => "failed",
            Self::Gated(_) => "gated",
            Self::Refused(_) => "refused",
        }
    }

    /// Did the stage actually run to a scientific pass/fail decision?
    #[must_use]
    pub const fn is_evaluated(&self) -> bool {
        matches!(self, Self::Passed | Self::Failed(_))
    }

    /// Did the stage actually run and pass?
    #[must_use]
    pub const fn is_passed(&self) -> bool {
        matches!(self, Self::Passed)
    }

    /// Structured reason for every non-passing disposition.
    #[must_use]
    pub const fn reason(&self) -> Option<&StageReason> {
        match self {
            Self::Passed => None,
            Self::Failed(reason) | Self::Gated(reason) | Self::Refused(reason) => Some(reason),
        }
    }

    fn is_well_formed(&self) -> bool {
        self.reason().is_none_or(StageReason::is_well_formed)
    }
}

impl core::fmt::Display for StageStatus {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        formatter.write_str(self.code())?;
        if let Some(reason) = self.reason() {
            write!(formatter, "{reason}")?;
        }
        Ok(())
    }
}

/// One stage's structured diagnostic result.
///
/// A `StageLog` is freely constructible DATA. By itself it carries no
/// promotion authority and cannot be inserted into an opaque
/// [`DiffRealReport`] by downstream callers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StageLog {
    /// The stage name.
    pub stage: &'static str,
    /// Whether this stage participates in the required-stage decision.
    pub requirement: StageRequirement,
    /// Typed scientific disposition; unavailable work is never a pass.
    pub status: StageStatus,
    /// Versioned identity of the fixture/schema whose result this log records.
    /// This is a diagnostic identity binding, not a content hash, proof
    /// certificate, independent verification receipt, or authorization.
    pub evidence_identity: &'static str,
    /// Typed deterministic log events.
    pub events: Vec<StageEvent>,
}

impl StageLog {
    /// Construct one plain diagnostic stage record.
    ///
    /// Construction does not confer authority or add the record to a
    /// [`DiffRealReport`].
    #[must_use]
    pub fn new(
        stage: &'static str,
        requirement: StageRequirement,
        status: StageStatus,
        evidence_identity: &'static str,
        events: Vec<StageEvent>,
    ) -> Self {
        Self {
            stage,
            requirement,
            status,
            evidence_identity,
            events,
        }
    }

    /// Did this stage actually run and pass?
    #[must_use]
    pub const fn passed(&self) -> bool {
        self.status.is_passed()
    }

    fn is_well_formed(&self) -> bool {
        !self.stage.trim().is_empty()
            && !self.evidence_identity.trim().is_empty()
            && !self.events.is_empty()
            && self.events.iter().all(StageEvent::is_well_formed)
            && self.status.is_well_formed()
    }
}

impl core::fmt::Display for StageLog {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            formatter,
            "stage={} requirement={} status={} evidence_identity={}",
            self.stage, self.requirement, self.status, self.evidence_identity
        )
    }
}

fn push_identity_field(output: &mut Vec<u8>, label: &str, value: &[u8]) {
    output.extend_from_slice(&(label.len() as u64).to_le_bytes());
    output.extend_from_slice(label.as_bytes());
    output.extend_from_slice(&(value.len() as u64).to_le_bytes());
    output.extend_from_slice(value);
}

fn push_identity_str(output: &mut Vec<u8>, label: &str, value: &str) {
    push_identity_field(output, label, value.as_bytes());
}

fn push_identity_u64(output: &mut Vec<u8>, label: &str, value: u64) {
    push_identity_field(output, label, &value.to_le_bytes());
}

fn push_identity_f64(output: &mut Vec<u8>, label: &str, value: f64) {
    push_identity_u64(output, label, value.to_bits());
}

fn identity_usize(value: usize) -> u64 {
    u64::try_from(value).expect("a Rust allocation length and index fit u64")
}

fn push_identity_usize(output: &mut Vec<u8>, label: &str, value: usize) {
    push_identity_u64(output, label, identity_usize(value));
}

fn push_optional_action(output: &mut Vec<u8>, label: &str, action: Option<Action>) {
    match action {
        Some(action) => {
            push_identity_field(output, &format!("{label}.present"), &[1]);
            push_action(output, &format!("{label}.value"), action);
        }
        None => push_identity_field(output, &format!("{label}.present"), &[0]),
    }
}

fn encode_transpose_error(output: &mut Vec<u8>, label: &str, error: &TransposeError) {
    match error {
        TransposeError::MissingVjp { op } => {
            push_identity_field(output, &format!("{label}.tag"), &[0]);
            push_identity_str(output, &format!("{label}.op"), op);
        }
        TransposeError::NonDifferentiableInPath {
            op,
            reason,
            color_consequence,
        } => {
            push_identity_field(output, &format!("{label}.tag"), &[1]);
            push_identity_str(output, &format!("{label}.op"), op);
            push_identity_str(output, &format!("{label}.reason"), reason);
            push_identity_str(
                output,
                &format!("{label}.color-consequence"),
                color_consequence,
            );
        }
    }
}

#[allow(clippy::too_many_lines)]
fn encode_differentiation_error(output: &mut Vec<u8>, label: &str, error: &DifferentiationError) {
    match error {
        DifferentiationError::Cancelled => {
            push_identity_field(output, &format!("{label}.tag"), &[0]);
        }
        DifferentiationError::WorkBudgetExceeded {
            required,
            available,
        } => {
            push_identity_field(output, &format!("{label}.tag"), &[1]);
            push_identity_u64(output, &format!("{label}.required"), *required);
            push_identity_u64(output, &format!("{label}.available"), *available);
        }
        DifferentiationError::EmptyPath => {
            push_identity_field(output, &format!("{label}.tag"), &[2]);
        }
        DifferentiationError::OraclePathMismatch {
            expected_len,
            observed_len,
            first_mismatch,
        } => {
            push_identity_field(output, &format!("{label}.tag"), &[3]);
            push_identity_usize(output, &format!("{label}.expected-len"), *expected_len);
            push_identity_usize(output, &format!("{label}.observed-len"), *observed_len);
            push_identity_usize(output, &format!("{label}.first-mismatch"), *first_mismatch);
        }
        DifferentiationError::PathTooLong { limit, observed } => {
            push_identity_field(output, &format!("{label}.tag"), &[4]);
            push_identity_usize(output, &format!("{label}.limit"), *limit);
            push_identity_usize(output, &format!("{label}.observed"), *observed);
        }
        DifferentiationError::OpNameTooLong {
            index,
            limit,
            observed,
        } => {
            push_identity_field(output, &format!("{label}.tag"), &[5]);
            push_identity_usize(output, &format!("{label}.index"), *index);
            push_identity_usize(output, &format!("{label}.limit"), *limit);
            push_identity_usize(output, &format!("{label}.observed"), *observed);
        }
        DifferentiationError::EmptyOpName { index } => {
            push_identity_field(output, &format!("{label}.tag"), &[6]);
            push_identity_usize(output, &format!("{label}.index"), *index);
        }
        DifferentiationError::MissingVjp { op } => {
            push_identity_field(output, &format!("{label}.tag"), &[7]);
            push_identity_str(output, &format!("{label}.op"), op);
        }
        DifferentiationError::UnsupportedOperator { op } => {
            push_identity_field(output, &format!("{label}.tag"), &[8]);
            push_identity_str(output, &format!("{label}.op"), op);
        }
        DifferentiationError::NonFiniteInput { bits } => {
            push_identity_field(output, &format!("{label}.tag"), &[9]);
            push_identity_u64(output, &format!("{label}.bits"), *bits);
        }
        DifferentiationError::NonFinitePrimal { op, bits } => {
            push_identity_field(output, &format!("{label}.tag"), &[10]);
            push_identity_str(output, &format!("{label}.op"), op);
            push_identity_u64(output, &format!("{label}.bits"), *bits);
        }
        DifferentiationError::Transpose(error) => {
            push_identity_field(output, &format!("{label}.tag"), &[11]);
            encode_transpose_error(output, &format!("{label}.transpose"), error);
        }
        DifferentiationError::MissingLeafGradient => {
            push_identity_field(output, &format!("{label}.tag"), &[12]);
        }
        DifferentiationError::InvalidGradientShape { observed } => {
            push_identity_field(output, &format!("{label}.tag"), &[13]);
            push_identity_usize(output, &format!("{label}.observed"), *observed);
        }
        DifferentiationError::NonFiniteGradient { bits } => {
            push_identity_field(output, &format!("{label}.tag"), &[14]);
            push_identity_u64(output, &format!("{label}.bits"), *bits);
        }
        DifferentiationError::OracleDisagreement {
            production_bits,
            dual_bits,
            fd_fine_bits,
            tolerance_bits,
        } => {
            push_identity_field(output, &format!("{label}.tag"), &[15]);
            push_identity_u64(
                output,
                &format!("{label}.production-bits"),
                *production_bits,
            );
            push_identity_u64(output, &format!("{label}.dual-bits"), *dual_bits);
            push_identity_u64(output, &format!("{label}.fd-fine-bits"), *fd_fine_bits);
            push_identity_u64(output, &format!("{label}.tolerance-bits"), *tolerance_bits);
        }
        DifferentiationError::InvalidInputScale { bits } => {
            push_identity_field(output, &format!("{label}.tag"), &[16]);
            push_identity_u64(output, &format!("{label}.bits"), *bits);
        }
        DifferentiationError::NonFiniteRescaledGradient { bits } => {
            push_identity_field(output, &format!("{label}.tag"), &[17]);
            push_identity_u64(output, &format!("{label}.bits"), *bits);
        }
        DifferentiationError::SensitivityIntegrityMismatch { identity } => {
            push_identity_field(output, &format!("{label}.tag"), &[18]);
            push_identity_field(output, &format!("{label}.identity"), identity.as_bytes());
        }
    }
}

#[allow(clippy::too_many_lines)]
fn encode_stage_event(output: &mut Vec<u8>, label: &str, event: &StageEvent) {
    match event {
        StageEvent::GradientVerified {
            receipt,
            input_bits,
            value_bits,
            production_bits,
            dual_bits,
            fd_coarse_bits,
            fd_fine_bits,
            tolerance_bits,
        } => {
            push_identity_field(output, &format!("{label}.tag"), &[0]);
            push_identity_field(output, &format!("{label}.receipt"), receipt.as_bytes());
            for (field, value) in [
                ("input-bits", *input_bits),
                ("value-bits", *value_bits),
                ("production-bits", *production_bits),
                ("dual-bits", *dual_bits),
                ("fd-coarse-bits", *fd_coarse_bits),
                ("fd-fine-bits", *fd_fine_bits),
                ("tolerance-bits", *tolerance_bits),
            ] {
                push_identity_u64(output, &format!("{label}.{field}"), value);
            }
        }
        StageEvent::DifferentiationRejected { error } => {
            push_identity_field(output, &format!("{label}.tag"), &[1]);
            encode_differentiation_error(output, &format!("{label}.error"), error);
        }
        StageEvent::MissingVjpProbe { op, blocked } => {
            push_identity_field(output, &format!("{label}.tag"), &[2]);
            push_identity_str(output, &format!("{label}.op"), op);
            push_identity_field(output, &format!("{label}.blocked"), &[u8::from(*blocked)]);
        }
        StageEvent::Registration {
            residual_bits,
            within_tolerance,
        } => {
            push_identity_field(output, &format!("{label}.tag"), &[3]);
            push_identity_u64(output, &format!("{label}.residual-bits"), *residual_bits);
            push_identity_field(
                output,
                &format!("{label}.within-tolerance"),
                &[u8::from(*within_tolerance)],
            );
        }
        StageEvent::AsBuiltDelta {
            max_deviation_bits,
            defect_index,
            estimated,
        } => {
            push_identity_field(output, &format!("{label}.tag"), &[4]);
            push_identity_u64(
                output,
                &format!("{label}.max-deviation-bits"),
                *max_deviation_bits,
            );
            match defect_index {
                Some(index) => {
                    push_identity_field(output, &format!("{label}.defect-index-present"), &[1]);
                    push_identity_usize(output, &format!("{label}.defect-index"), *index);
                }
                None => push_identity_field(output, &format!("{label}.defect-index-present"), &[0]),
            }
            push_identity_field(
                output,
                &format!("{label}.estimated"),
                &[u8::from(*estimated)],
            );
        }
        StageEvent::Assimilation {
            before_bits,
            after_bits,
            reduced,
        } => {
            push_identity_field(output, &format!("{label}.tag"), &[5]);
            push_identity_u64(output, &format!("{label}.before-bits"), *before_bits);
            push_identity_u64(output, &format!("{label}.after-bits"), *after_bits);
            push_identity_field(output, &format!("{label}.reduced"), &[u8::from(*reduced)]);
        }
        StageEvent::ToleranceActions { critical, slack } => {
            push_identity_field(output, &format!("{label}.tag"), &[6]);
            push_optional_action(output, &format!("{label}.critical"), *critical);
            push_optional_action(output, &format!("{label}.slack"), *slack);
        }
        StageEvent::GdtJustification {
            loosened,
            all_verified,
        } => {
            push_identity_field(output, &format!("{label}.tag"), &[7]);
            push_identity_usize(output, &format!("{label}.loosened"), *loosened);
            push_identity_field(
                output,
                &format!("{label}.all-verified"),
                &[u8::from(*all_verified)],
            );
        }
        StageEvent::SampledLinearization {
            samples,
            confirmed,
            linearized_std_bits,
            probability_claimed,
        } => {
            push_identity_field(output, &format!("{label}.tag"), &[8]);
            push_identity_usize(output, &format!("{label}.samples"), *samples);
            push_identity_field(
                output,
                &format!("{label}.confirmed"),
                &[u8::from(*confirmed)],
            );
            push_identity_u64(
                output,
                &format!("{label}.linearized-std-bits"),
                *linearized_std_bits,
            );
            push_identity_field(
                output,
                &format!("{label}.probability-claimed"),
                &[u8::from(*probability_claimed)],
            );
        }
        StageEvent::Gate { code, detail } => {
            push_identity_field(output, &format!("{label}.tag"), &[9]);
            push_identity_str(output, &format!("{label}.code"), code);
            push_identity_str(output, &format!("{label}.detail"), detail);
        }
        StageEvent::Refusal { code, detail } => {
            push_identity_field(output, &format!("{label}.tag"), &[10]);
            push_identity_str(output, &format!("{label}.code"), code);
            push_identity_str(output, &format!("{label}.detail"), detail);
        }
    }
}

/// Exact execution provenance bound into every stage receipt and the report
/// root. Construction is crate-controlled through [`Cx`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DiffRealExecutionIdentity {
    stream_key: StreamKey,
    budget: Budget,
    mode: ExecMode,
    operation_memory_limit_bytes: Option<u64>,
}

impl DiffRealExecutionIdentity {
    fn from_cx(cx: &Cx<'_>) -> Self {
        Self {
            stream_key: cx.stream_key(),
            budget: cx.budget(),
            mode: cx.mode(),
            operation_memory_limit_bytes: cx.lease().and_then(|lease| lease.limit_bytes()),
        }
    }

    /// Logical stream identity used by the battery run.
    #[must_use]
    pub const fn stream_key(&self) -> StreamKey {
        self.stream_key
    }

    /// Ambient budget captured before the battery ran.
    #[must_use]
    pub const fn budget(&self) -> Budget {
        self.budget
    }

    /// Execution mode captured before the battery ran.
    #[must_use]
    pub const fn mode(&self) -> ExecMode {
        self.mode
    }

    /// Configured operation-memory limit in force for the source context.
    /// `None` identifies a manually constructed or unbounded context.
    #[must_use]
    pub const fn operation_memory_limit_bytes(&self) -> Option<u64> {
        self.operation_memory_limit_bytes
    }

    fn matches_cx(&self, cx: &Cx<'_>) -> bool {
        *self == Self::from_cx(cx)
    }
}

fn encode_execution_identity(output: &mut Vec<u8>, execution: &DiffRealExecutionIdentity) {
    push_identity_field(
        output,
        "execution.mode-tag",
        &[match execution.mode {
            ExecMode::Deterministic => 0,
            ExecMode::Fast => 1,
        }],
    );
    push_identity_u64(output, "execution.stream.seed", execution.stream_key.seed);
    push_identity_u64(
        output,
        "execution.stream.kernel-id",
        execution.stream_key.kernel_id,
    );
    push_identity_u64(output, "execution.stream.tile", execution.stream_key.tile);
    push_identity_u64(
        output,
        "execution.stream.iteration",
        execution.stream_key.iteration,
    );
    match execution.budget.deadline {
        Some(deadline) => {
            push_identity_field(output, "execution.budget.deadline-present", &[1]);
            push_identity_u64(
                output,
                "execution.budget.deadline-nanos",
                deadline.as_nanos(),
            );
        }
        None => push_identity_field(output, "execution.budget.deadline-present", &[0]),
    }
    push_identity_field(
        output,
        "execution.budget.poll-quota",
        &execution.budget.poll_quota.to_le_bytes(),
    );
    match execution.budget.cost_quota {
        Some(cost) => {
            push_identity_field(output, "execution.budget.cost-quota-present", &[1]);
            push_identity_u64(output, "execution.budget.cost-quota", cost);
        }
        None => push_identity_field(output, "execution.budget.cost-quota-present", &[0]),
    }
    push_identity_field(
        output,
        "execution.budget.priority",
        &[execution.budget.priority],
    );
    match execution.operation_memory_limit_bytes {
        Some(limit) => {
            push_identity_field(output, "execution.memory-limit-present", &[1]);
            push_identity_u64(output, "execution.memory-limit-bytes", limit);
        }
        None => push_identity_field(output, "execution.memory-limit-present", &[0]),
    }
}

fn push_fixture_path(output: &mut Vec<u8>, label: &str, path: &[&str]) {
    push_identity_u64(output, &format!("{label}.length"), path.len() as u64);
    for (index, op) in path.iter().enumerate() {
        push_identity_str(output, &format!("{label}.{index}"), op);
    }
}

#[allow(clippy::too_many_lines)]
fn fixture_input_identity(stage: &StageLog) -> ContentHash {
    let mut canonical = Vec::new();
    push_identity_u64(
        &mut canonical,
        "stage-receipt-policy-version",
        u64::from(STAGE_RECEIPT_POLICY_VERSION),
    );
    push_identity_str(&mut canonical, "stage", stage.stage);
    push_identity_str(&mut canonical, "evidence-identity", stage.evidence_identity);

    match stage.stage {
        DIFFERENTIATION_STAGE => {
            push_fixture_path(
                &mut canonical,
                "production-path",
                &PRODUCTION_DIFFERENTIATION_PATH,
            );
            push_identity_f64(
                &mut canonical,
                "production-input",
                DIFFERENTIATION_FIXTURE_INPUT,
            );
            push_fixture_path(
                &mut canonical,
                "missing-vjp-path",
                &MISSING_VJP_FIXTURE_PATH,
            );
            push_identity_f64(
                &mut canonical,
                "missing-vjp-input",
                DIFFERENTIATION_FIXTURE_INPUT,
            );
            push_identity_str(
                &mut canonical,
                "registry-policy",
                DIFFERENTIATION_REGISTRY_POLICY,
            );
            push_identity_str(
                &mut canonical,
                "sensitivity-policy",
                SENSITIVITY_POLICY_VERSION,
            );
            push_identity_u64(
                &mut canonical,
                "admitted-work-units",
                DIFFERENTIATION_STAGE_WORK_UNITS,
            );
        }
        AS_BUILT_STAGE => {
            for (index, (x, y)) in AS_BUILT_DESIGN_POINTS.iter().copied().enumerate() {
                push_identity_f64(&mut canonical, &format!("design.{index}.x"), x);
                push_identity_f64(&mut canonical, &format!("design.{index}.y"), y);
            }
            push_identity_f64(
                &mut canonical,
                "registration.rotation-radians",
                AS_BUILT_ROTATION_RADIANS,
            );
            push_identity_f64(
                &mut canonical,
                "registration.translation-x",
                AS_BUILT_TRANSLATION.0,
            );
            push_identity_f64(
                &mut canonical,
                "registration.translation-y",
                AS_BUILT_TRANSLATION.1,
            );
            push_identity_u64(&mut canonical, "defect.index", AS_BUILT_DEFECT_INDEX as u64);
            push_identity_f64(&mut canonical, "defect.delta-x", AS_BUILT_DEFECT_X);
            push_identity_f64(
                &mut canonical,
                "design-tolerance",
                AS_BUILT_DESIGN_TOLERANCE,
            );
            push_identity_f64(
                &mut canonical,
                "measurement-noise",
                AS_BUILT_MEASUREMENT_NOISE,
            );
            push_identity_str(
                &mut canonical,
                "calibration-candidate",
                AS_BUILT_CALIBRATION_CANDIDATE,
            );
            for (index, value) in AS_BUILT_PRIOR_MEAN.iter().copied().enumerate() {
                push_identity_f64(&mut canonical, &format!("prior.mean.{index}"), value);
            }
            for (index, value) in AS_BUILT_PRIOR_DIAGONAL_COVARIANCE
                .iter()
                .copied()
                .enumerate()
            {
                push_identity_f64(
                    &mut canonical,
                    &format!("prior.diagonal-covariance.{index}"),
                    value,
                );
            }
            for (index, (state_index, value, variance, label)) in
                AS_BUILT_OBSERVATIONS.iter().copied().enumerate()
            {
                push_identity_u64(
                    &mut canonical,
                    &format!("observation.{index}.state-dimension"),
                    AS_BUILT_PRIOR_MEAN.len() as u64,
                );
                push_identity_u64(
                    &mut canonical,
                    &format!("observation.{index}.state-index"),
                    state_index as u64,
                );
                push_identity_f64(&mut canonical, &format!("observation.{index}.value"), value);
                push_identity_f64(
                    &mut canonical,
                    &format!("observation.{index}.variance"),
                    variance,
                );
                push_identity_str(&mut canonical, &format!("observation.{index}.label"), label);
            }
            push_identity_str(
                &mut canonical,
                "assimilation.parameter",
                AS_BUILT_ASSIMILATION_PARAMETER,
            );
            push_identity_f64(
                &mut canonical,
                "assimilation.lower",
                AS_BUILT_ASSIMILATION_BOUNDS.0,
            );
            push_identity_f64(
                &mut canonical,
                "assimilation.upper",
                AS_BUILT_ASSIMILATION_BOUNDS.1,
            );
            push_identity_u64(
                &mut canonical,
                "as-built-work-plan-version",
                u64::from(fs_asbuilt::AS_BUILT_WORK_PLAN_VERSION),
            );
            push_identity_u64(
                &mut canonical,
                "as-built-poll-policy-version",
                u64::from(fs_asbuilt::AS_BUILT_POLL_POLICY_VERSION),
            );
            push_identity_u64(
                &mut canonical,
                "assimilation-psd-policy-version",
                u64::from(fs_assimilate::PSD_ADMISSION_POLICY_VERSION),
            );
            push_identity_str(
                &mut canonical,
                "assimilation-poll-policy",
                ASSIMILATION_POLL_POLICY,
            );
            push_identity_u64(
                &mut canonical,
                "color-algebra-version",
                u64::from(fs_evidence::COLOR_ALGEBRA_VERSION),
            );
            push_identity_str(
                &mut canonical,
                "invocation-budget-policy",
                invocation::AS_BUILT_INVOCATION_POLICY,
            );
            push_identity_u64(
                &mut canonical,
                "invocation-receipt-version",
                u64::from(fs_exec::INVOCATION_RECEIPT_VERSION),
            );
        }
        TOLERANCE_STAGE => {
            push_fixture_path(
                &mut canonical,
                "production-path",
                &PRODUCTION_DIFFERENTIATION_PATH,
            );
            for (index, input) in TOLERANCE_SENSITIVITY_INPUTS.iter().copied().enumerate() {
                push_identity_f64(&mut canonical, &format!("sensitivity-input.{index}"), input);
            }
            for (index, name) in TOLERANCE_FEATURE_NAMES.iter().copied().enumerate() {
                push_identity_str(&mut canonical, &format!("feature.{index}.name"), name);
                push_identity_f64(
                    &mut canonical,
                    &format!("feature.{index}.cost-coefficient"),
                    TOLERANCE_FEATURE_COST_COEFFICIENT,
                );
                push_identity_f64(
                    &mut canonical,
                    &format!("feature.{index}.baseline-tolerance"),
                    TOLERANCE_FEATURE_BASELINE,
                );
            }
            for (index, sample) in TOLERANCE_EXTREME_QOIS.iter().copied().enumerate() {
                push_identity_f64(&mut canonical, &format!("extreme-qoi.{index}"), sample);
            }
            for (label, value) in [
                (
                    "variance.performance-tolerance",
                    TOLERANCE_PERFORMANCE_TOLERANCE,
                ),
                ("variance.target-probability", TOLERANCE_TARGET_PROBABILITY),
                ("allocation.sigma-multiplier", TOLERANCE_SIGMA_MULTIPLIER),
                ("robustness.nominal-qoi", TOLERANCE_NOMINAL_QOI),
                ("robustness.sigma-multiplier", TOLERANCE_SIGMA_MULTIPLIER),
                ("robustness.relative-margin", TOLERANCE_ROBUSTNESS_MARGIN),
            ] {
                push_identity_f64(&mut canonical, label, value);
            }
            push_identity_str(
                &mut canonical,
                "sensitivity-policy",
                SENSITIVITY_POLICY_VERSION,
            );
            push_identity_u64(
                &mut canonical,
                "color-algebra-version",
                u64::from(fs_evidence::COLOR_ALGEBRA_VERSION),
            );
            push_identity_u64(&mut canonical, "admitted-work-units", TOLERANCE_WORK_UNITS);
        }
        SPACETIME_STAGE => {
            push_identity_str(
                &mut canonical,
                "gate-code",
                "diffreal.spacetime.integration-not-activated",
            );
            push_identity_str(
                &mut canonical,
                "declared-capability",
                "fs-time/temporal-complex",
            );
            push_identity_str(
                &mut canonical,
                "dependency-bead",
                "frankensim-epic-coupling-bk0o.7",
            );
            push_identity_u64(&mut canonical, "admitted-work-units", SPACETIME_WORK_UNITS);
        }
        _ => {
            push_identity_str(
                &mut canonical,
                "optional-stage-input-boundary",
                "diagnostic-only; no required-stage promotion authority",
            );
        }
    }

    hash_domain(FIXTURE_INPUT_IDENTITY_DOMAIN, &canonical)
}

fn push_validity_domain(output: &mut Vec<u8>, label: &str, regime: &fs_evidence::ValidityDomain) {
    push_identity_u64(
        output,
        &format!("{label}.axis-count"),
        regime.bounds().len() as u64,
    );
    for (index, (axis, (lo, hi))) in regime.bounds().iter().enumerate() {
        push_identity_str(output, &format!("{label}.{index}.axis"), axis);
        push_identity_f64(output, &format!("{label}.{index}.lo"), *lo);
        push_identity_f64(output, &format!("{label}.{index}.hi"), *hi);
    }
}

fn push_color(output: &mut Vec<u8>, label: &str, color: &Color) {
    match color {
        Color::Verified { lo, hi } => {
            push_identity_field(output, &format!("{label}.rank-tag"), &[0]);
            push_identity_f64(output, &format!("{label}.lo"), *lo);
            push_identity_f64(output, &format!("{label}.hi"), *hi);
        }
        Color::Validated { regime, dataset } => {
            push_identity_field(output, &format!("{label}.rank-tag"), &[1]);
            push_identity_str(output, &format!("{label}.dataset"), dataset);
            push_validity_domain(output, &format!("{label}.regime"), regime);
        }
        Color::Estimated {
            estimator,
            dispersion,
        } => {
            push_identity_field(output, &format!("{label}.rank-tag"), &[2]);
            push_identity_str(output, &format!("{label}.estimator"), estimator);
            push_identity_f64(output, &format!("{label}.dispersion"), *dispersion);
        }
    }
}

fn push_color_rank(output: &mut Vec<u8>, label: &str, color: ColorRank) {
    push_identity_field(
        output,
        label,
        &[match color {
            ColorRank::Estimated => 0,
            ColorRank::Validated => 1,
            ColorRank::Verified => 2,
        }],
    );
}

fn push_action(output: &mut Vec<u8>, label: &str, action: Action) {
    push_identity_field(
        output,
        label,
        &[match action {
            Action::Tighten => 0,
            Action::Loosen => 1,
            Action::Unchanged => 2,
        }],
    );
}

fn push_sensitivity(output: &mut Vec<u8>, label: &str, sensitivity: &SealedSensitivity) {
    push_identity_u64(
        output,
        &format!("{label}.operator-count"),
        sensitivity.ops.len() as u64,
    );
    for (index, op) in sensitivity.ops.iter().enumerate() {
        push_identity_str(output, &format!("{label}.operator.{index}"), op);
    }
    for (field, value) in [
        ("input", sensitivity.input_bits),
        ("value", sensitivity.value_bits),
        ("production-gradient", sensitivity.production_gradient_bits),
        ("dual-gradient", sensitivity.dual_gradient_bits),
        ("fd-coarse", sensitivity.fd_coarse_bits),
        ("fd-fine", sensitivity.fd_fine_bits),
        ("fd-tolerance", sensitivity.fd_tolerance_bits),
    ] {
        push_identity_u64(output, &format!("{label}.{field}-bits"), value);
    }
    push_identity_field(
        output,
        &format!("{label}.identity"),
        sensitivity.identity.as_bytes(),
    );
}

fn finish_stage_result(stage: &str, canonical: Vec<u8>) -> ContentHash {
    let mut framed = Vec::new();
    push_identity_u64(
        &mut framed,
        "stage-receipt-policy-version",
        u64::from(STAGE_RECEIPT_POLICY_VERSION),
    );
    push_identity_str(&mut framed, "stage", stage);
    push_identity_field(&mut framed, "result-payload", &canonical);
    hash_domain(STAGE_RESULT_IDENTITY_DOMAIN, &framed)
}

fn diagnostic_result_identity(stage: &StageLog) -> ContentHash {
    let mut canonical = Vec::new();
    let status_tag = match &stage.status {
        StageStatus::Passed => 0,
        StageStatus::Failed(_) => 1,
        StageStatus::Gated(_) => 2,
        StageStatus::Refused(_) => 3,
    };
    push_identity_field(&mut canonical, "status-tag", &[status_tag]);
    if let Some(reason) = stage.status.reason() {
        push_identity_str(&mut canonical, "reason.code", reason.code);
        push_identity_str(&mut canonical, "reason.detail", &reason.detail);
    }
    push_identity_u64(&mut canonical, "event-count", stage.events.len() as u64);
    for (index, event) in stage.events.iter().enumerate() {
        encode_stage_event(&mut canonical, &format!("event.{index}"), event);
    }
    finish_stage_result(stage.stage, canonical)
}

fn differentiation_result_identity(
    sensitivity: &SealedSensitivity,
    missing_probe: &Result<PathDerivative, DifferentiationError>,
) -> ContentHash {
    let mut canonical = Vec::new();
    push_sensitivity(&mut canonical, "sealed-sensitivity", sensitivity);
    match missing_probe {
        Ok(derivative) => {
            push_identity_str(
                &mut canonical,
                "missing-probe.disposition",
                "unexpected-pass",
            );
            push_identity_u64(
                &mut canonical,
                "missing-probe.value-bits",
                derivative.value_bits,
            );
            push_identity_u64(
                &mut canonical,
                "missing-probe.gradient-bits",
                derivative.gradient_bits,
            );
        }
        Err(error) => {
            push_identity_str(&mut canonical, "missing-probe.disposition", "refused");
            encode_differentiation_error(&mut canonical, "missing-probe.typed-error", error);
        }
    }
    finish_stage_result(DIFFERENTIATION_STAGE, canonical)
}

fn as_built_result_identity(
    registration: &fs_asbuilt::Registration,
    difference: &fs_asbuilt::AsBuiltDiff,
    posterior: &fs_assimilate::AssimilatedPosterior,
) -> ContentHash {
    let mut canonical = Vec::new();
    for (label, value) in [
        ("registration.rotation-radians", registration.rotation_rad()),
        ("registration.translation-x", registration.tx()),
        ("registration.translation-y", registration.ty()),
        ("registration.residual-rms", registration.residual_rms()),
    ] {
        push_identity_f64(&mut canonical, label, value);
    }
    push_identity_u64(
        &mut canonical,
        "difference.deviation-count",
        difference.deviations().len() as u64,
    );
    for (index, deviation) in difference.deviations().iter().copied().enumerate() {
        push_identity_f64(
            &mut canonical,
            &format!("difference.deviation.{index}"),
            deviation,
        );
    }
    push_identity_f64(
        &mut canonical,
        "difference.max-deviation",
        difference.max_deviation(),
    );
    push_identity_u64(
        &mut canonical,
        "difference.max-deviation-index",
        difference.max_deviation_index() as u64,
    );
    push_identity_field(
        &mut canonical,
        "difference.within-tolerance",
        &[u8::from(difference.within_tolerance())],
    );
    push_identity_field(
        &mut canonical,
        "difference.above-noise-floor",
        &[u8::from(difference.above_noise_floor())],
    );
    push_validity_domain(
        &mut canonical,
        "difference.proposed-regime",
        difference.proposed_regime(),
    );
    push_color(&mut canonical, "difference.color", difference.color());

    let belief = posterior.belief();
    push_identity_u64(&mut canonical, "posterior.dimension", belief.dim() as u64);
    for (index, value) in belief.mean().iter().copied().enumerate() {
        push_identity_f64(&mut canonical, &format!("posterior.mean.{index}"), value);
    }
    for (row, values) in belief.covariance().iter().enumerate() {
        for (column, value) in values.iter().copied().enumerate() {
            push_identity_f64(
                &mut canonical,
                &format!("posterior.covariance.{row}.{column}"),
                value,
            );
        }
    }
    push_color(&mut canonical, "posterior.color", posterior.color());
    push_validity_domain(&mut canonical, "posterior.regime", posterior.regime());
    push_identity_f64(
        &mut canonical,
        "posterior.misfit-before",
        posterior.misfit_before(),
    );
    push_identity_f64(
        &mut canonical,
        "posterior.misfit-after",
        posterior.misfit_after(),
    );
    finish_stage_result(AS_BUILT_STAGE, canonical)
}

fn tolerance_result_identity(
    critical: &SealedSensitivity,
    slack: &SealedSensitivity,
    allocation: &fs_toleralloc::Allocation,
    report: &[fs_toleralloc::Suggestion],
    verdict: &fs_toleralloc::RobustnessVerdict,
) -> ContentHash {
    let mut canonical = Vec::new();
    push_sensitivity(&mut canonical, "critical", critical);
    push_sensitivity(&mut canonical, "slack", slack);
    push_identity_u64(
        &mut canonical,
        "allocation.item-count",
        allocation.items.len() as u64,
    );
    for (index, item) in allocation.items.iter().enumerate() {
        push_identity_str(
            &mut canonical,
            &format!("allocation.item.{index}.name"),
            &item.name,
        );
        push_identity_f64(
            &mut canonical,
            &format!("allocation.item.{index}.tolerance"),
            item.tolerance,
        );
        push_identity_f64(
            &mut canonical,
            &format!("allocation.item.{index}.sensitivity"),
            item.sensitivity,
        );
        push_color_rank(
            &mut canonical,
            &format!("allocation.item.{index}.color"),
            item.sensitivity_color,
        );
        push_action(
            &mut canonical,
            &format!("allocation.item.{index}.action"),
            item.action,
        );
    }
    push_identity_f64(
        &mut canonical,
        "allocation.total-cost",
        allocation.total_cost,
    );
    push_identity_f64(
        &mut canonical,
        "allocation.achieved-variance",
        allocation.achieved_variance,
    );
    push_identity_u64(&mut canonical, "gdt-row-count", report.len() as u64);
    for (index, row) in report.iter().enumerate() {
        push_identity_str(&mut canonical, &format!("gdt.{index}.name"), &row.name);
        push_identity_f64(
            &mut canonical,
            &format!("gdt.{index}.tolerance"),
            row.tolerance,
        );
        push_action(&mut canonical, &format!("gdt.{index}.action"), row.action);
        push_identity_f64(
            &mut canonical,
            &format!("gdt.{index}.certified-sensitivity"),
            row.certified_sensitivity,
        );
        push_color_rank(&mut canonical, &format!("gdt.{index}.color"), row.color);
    }
    push_identity_f64(
        &mut canonical,
        "robustness.linearized-std",
        verdict.linearized_std,
    );
    push_identity_f64(
        &mut canonical,
        "robustness.sampled-max-deviation",
        verdict.sampled_max_deviation,
    );
    push_identity_field(
        &mut canonical,
        "robustness.confirmed",
        &[u8::from(verdict.confirmed)],
    );
    finish_stage_result(TOLERANCE_STAGE, canonical)
}

#[derive(Debug, Clone)]
struct StageExecution {
    log: StageLog,
    inputs: ContentHash,
    result: ContentHash,
    invocation: Option<fs_exec::InvocationReceipt>,
}

impl StageExecution {
    fn diagnostic(log: StageLog) -> Self {
        let inputs = fixture_input_identity(&log);
        let result = diagnostic_result_identity(&log);
        Self {
            log,
            inputs,
            result,
            invocation: None,
        }
    }
}

fn stage_receipt_identity(
    version: u32,
    stage: &StageLog,
    execution: &DiffRealExecutionIdentity,
    fixture_inputs: ContentHash,
    result: ContentHash,
    invocation: Option<&fs_exec::InvocationReceipt>,
) -> ContentHash {
    let mut canonical = Vec::new();
    push_identity_u64(&mut canonical, "receipt-version", u64::from(version));
    push_identity_str(&mut canonical, "stage", stage.stage);
    push_identity_field(
        &mut canonical,
        "requirement-tag",
        &[match stage.requirement {
            StageRequirement::Required => 0,
            StageRequirement::Optional => 1,
        }],
    );
    push_identity_field(
        &mut canonical,
        "status-tag",
        &[match &stage.status {
            StageStatus::Passed => 0,
            StageStatus::Failed(_) => 1,
            StageStatus::Gated(_) => 2,
            StageStatus::Refused(_) => 3,
        }],
    );
    if let Some(reason) = stage.status.reason() {
        push_identity_str(&mut canonical, "status.reason.code", reason.code);
        push_identity_str(&mut canonical, "status.reason.detail", &reason.detail);
    } else {
        push_identity_str(&mut canonical, "status.reason", "none");
    }
    push_identity_str(&mut canonical, "evidence-identity", stage.evidence_identity);
    push_identity_field(
        &mut canonical,
        "fixture-input-root",
        fixture_inputs.as_bytes(),
    );
    push_identity_field(&mut canonical, "stage-result-root", result.as_bytes());
    push_identity_field(
        &mut canonical,
        "invocation-receipt-present",
        &[u8::from(invocation.is_some())],
    );
    if let Some(invocation) = invocation {
        push_identity_u64(
            &mut canonical,
            "invocation-receipt-version",
            u64::from(invocation.version()),
        );
        push_identity_field(
            &mut canonical,
            "invocation-id",
            invocation.invocation_id().as_bytes(),
        );
        push_identity_field(
            &mut canonical,
            "invocation-receipt-root",
            invocation.root().as_bytes(),
        );
        push_identity_field(
            &mut canonical,
            "invocation-disposition",
            &[match invocation.disposition() {
                fs_exec::InvocationDisposition::Completed => 0,
                fs_exec::InvocationDisposition::Cancelled => 1,
                fs_exec::InvocationDisposition::Refused => 2,
            }],
        );
    }
    push_identity_u64(&mut canonical, "event-count", stage.events.len() as u64);
    for (index, event) in stage.events.iter().enumerate() {
        encode_stage_event(&mut canonical, &format!("event.{index}"), event);
    }
    encode_execution_identity(&mut canonical, execution);
    push_identity_str(&mut canonical, "promotion-policy", REPORT_PROMOTION_POLICY);
    hash_domain(STAGE_RECEIPT_IDENTITY_DOMAIN, &canonical)
}

/// Opaque, content-addressed receipt for one crate-authored stage execution.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StageReceipt {
    version: u32,
    stage: &'static str,
    fixture_inputs: ContentHash,
    result: ContentHash,
    invocation: Option<fs_exec::InvocationReceipt>,
    root: ContentHash,
}

impl StageReceipt {
    fn seal(stage: &StageExecution, execution: &DiffRealExecutionIdentity) -> Self {
        let fixture_inputs = stage.inputs;
        let invocation = stage.invocation.clone();
        let root = stage_receipt_identity(
            STAGE_RECEIPT_POLICY_VERSION,
            &stage.log,
            execution,
            fixture_inputs,
            stage.result,
            invocation.as_ref(),
        );
        Self {
            version: STAGE_RECEIPT_POLICY_VERSION,
            stage: stage.log.stage,
            fixture_inputs,
            result: stage.result,
            invocation,
            root,
        }
    }

    /// Receipt schema/policy version.
    #[must_use]
    pub const fn version(&self) -> u32 {
        self.version
    }

    /// Stable stage name authenticated by this receipt.
    #[must_use]
    pub const fn stage(&self) -> &'static str {
        self.stage
    }

    /// Content address of the exact fixed-fixture inputs and policy versions.
    #[must_use]
    pub const fn fixture_inputs(&self) -> ContentHash {
        self.fixture_inputs
    }

    /// Content address of the complete stage result retained by the fixed
    /// battery implementation, including values summarized out of diagnostics.
    #[must_use]
    pub const fn result(&self) -> ContentHash {
        self.result
    }

    /// Full affine invocation receipt for stages that execute a composed
    /// scientific transaction.
    #[must_use]
    pub const fn invocation(&self) -> Option<&fs_exec::InvocationReceipt> {
        self.invocation.as_ref()
    }

    /// Content address of the inputs, result log, execution identity, budgets,
    /// and policy versions for this stage.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }

    fn verifies(&self, stage: &StageLog, execution: &DiffRealExecutionIdentity) -> bool {
        let invocation_valid = match (&self.invocation, stage.stage) {
            (Some(invocation), AS_BUILT_STAGE) => {
                invocation::AsBuiltInvocationPlan::preflight(execution)
                    .is_ok_and(|plan| plan.verifies_receipt(invocation, execution))
            }
            (None, AS_BUILT_STAGE) | (Some(_), _) => false,
            (None, _) => true,
        };
        self.version == STAGE_RECEIPT_POLICY_VERSION
            && self.stage == stage.stage
            && invocation_valid
            && self.root
                == stage_receipt_identity(
                    self.version,
                    stage,
                    execution,
                    self.fixture_inputs,
                    self.result,
                    self.invocation.as_ref(),
                )
    }
}

fn report_receipt_identity(
    version: u32,
    execution: &DiffRealExecutionIdentity,
    receipts: &[StageReceipt],
) -> ContentHash {
    let mut canonical = Vec::new();
    push_identity_u64(&mut canonical, "receipt-version", u64::from(version));
    push_identity_str(&mut canonical, "promotion-policy", REPORT_PROMOTION_POLICY);
    encode_execution_identity(&mut canonical, execution);
    push_identity_u64(&mut canonical, "stage-count", receipts.len() as u64);
    for (index, receipt) in receipts.iter().enumerate() {
        push_identity_str(
            &mut canonical,
            &format!("stage.{index}.name"),
            receipt.stage,
        );
        push_identity_field(
            &mut canonical,
            &format!("stage.{index}.receipt-root"),
            receipt.root.as_bytes(),
        );
    }
    hash_domain(REPORT_RECEIPT_IDENTITY_DOMAIN, &canonical)
}

/// Fingerprint of every local rule whose agreement is required before an
/// external authority may authenticate a report root.
#[must_use]
pub fn promotion_policy_fingerprint() -> ContentHash {
    let mut canonical = Vec::new();
    push_identity_str(&mut canonical, "promotion-policy", REPORT_PROMOTION_POLICY);
    push_identity_u64(
        &mut canonical,
        "stage-receipt-version",
        u64::from(STAGE_RECEIPT_POLICY_VERSION),
    );
    push_identity_u64(
        &mut canonical,
        "report-receipt-version",
        u64::from(REPORT_RECEIPT_POLICY_VERSION),
    );
    push_identity_str(
        &mut canonical,
        "sensitivity-policy",
        SENSITIVITY_POLICY_VERSION,
    );
    push_identity_str(
        &mut canonical,
        "differentiation-registry-policy",
        DIFFERENTIATION_REGISTRY_POLICY,
    );
    push_identity_u64(
        &mut canonical,
        "as-built-work-plan-version",
        u64::from(fs_asbuilt::AS_BUILT_WORK_PLAN_VERSION),
    );
    push_identity_u64(
        &mut canonical,
        "as-built-poll-policy-version",
        u64::from(fs_asbuilt::AS_BUILT_POLL_POLICY_VERSION),
    );
    push_identity_str(
        &mut canonical,
        "as-built-invocation-policy",
        invocation::AS_BUILT_INVOCATION_POLICY,
    );
    push_identity_u64(
        &mut canonical,
        "invocation-receipt-version",
        u64::from(fs_exec::INVOCATION_RECEIPT_VERSION),
    );
    push_identity_u64(
        &mut canonical,
        "assimilation-psd-policy-version",
        u64::from(fs_assimilate::PSD_ADMISSION_POLICY_VERSION),
    );
    push_identity_str(
        &mut canonical,
        "assimilation-poll-policy",
        ASSIMILATION_POLL_POLICY,
    );
    push_identity_u64(
        &mut canonical,
        "color-algebra-version",
        u64::from(fs_evidence::COLOR_ALGEBRA_VERSION),
    );
    for (index, required) in REQUIRED_STAGES.iter().enumerate() {
        push_identity_str(
            &mut canonical,
            &format!("required-stage.{index}.name"),
            required.name,
        );
        push_identity_str(
            &mut canonical,
            &format!("required-stage.{index}.evidence-identity"),
            required.evidence_identity,
        );
    }
    hash_domain(PROMOTION_POLICY_FINGERPRINT_DOMAIN, &canonical)
}

fn promotion_verification_subject(
    report_root: ContentHash,
    execution: &DiffRealExecutionIdentity,
    report_schema_version: u32,
    stage_schema_version: u32,
    policy_fingerprint: ContentHash,
) -> ContentHash {
    let mut canonical = Vec::new();
    push_identity_str(&mut canonical, "purpose", PROMOTION_VERIFICATION_PURPOSE);
    push_identity_field(&mut canonical, "report-root", report_root.as_bytes());
    encode_execution_identity(&mut canonical, execution);
    push_identity_u64(
        &mut canonical,
        "report-schema-version",
        u64::from(report_schema_version),
    );
    push_identity_u64(
        &mut canonical,
        "stage-schema-version",
        u64::from(stage_schema_version),
    );
    push_identity_field(
        &mut canonical,
        "policy-fingerprint",
        policy_fingerprint.as_bytes(),
    );
    hash_domain(PROMOTION_VERIFICATION_SUBJECT_DOMAIN, &canonical)
}

/// Detached report-root attestation supplied by an external authority.
/// Construction alone grants no authority; a verifier must accept it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PromotionAttestation {
    key_id: String,
    signature: Vec<u8>,
    policy_fingerprint: ContentHash,
}

impl PromotionAttestation {
    /// Assemble detached attestation data for an injected authority verifier.
    #[must_use]
    pub fn new(
        key_id: impl Into<String>,
        signature: impl Into<Vec<u8>>,
        policy_fingerprint: ContentHash,
    ) -> Self {
        Self {
            key_id: key_id.into(),
            signature: signature.into(),
            policy_fingerprint,
        }
    }

    /// Authority key identity interpreted by the injected verifier.
    #[must_use]
    pub fn key_id(&self) -> &str {
        &self.key_id
    }

    /// Detached authentication bytes interpreted by the injected verifier.
    #[must_use]
    pub fn signature(&self) -> &[u8] {
        &self.signature
    }

    /// Authority policy fingerprint claimed by the attestation.
    #[must_use]
    pub const fn policy_fingerprint(&self) -> ContentHash {
        self.policy_fingerprint
    }
}

/// Immutable request presented exactly once to a promotion authority.
#[derive(Debug, Clone, Copy)]
pub struct PromotionVerificationRequest<'a> {
    purpose: &'static str,
    report_root: ContentHash,
    execution: &'a DiffRealExecutionIdentity,
    report_schema_version: u32,
    stage_schema_version: u32,
    policy_fingerprint: ContentHash,
    subject: ContentHash,
}

impl PromotionVerificationRequest<'_> {
    /// Domain-separated authority purpose.
    #[must_use]
    pub const fn purpose(&self) -> &'static str {
        self.purpose
    }

    /// Ordered report root presented for authentication.
    #[must_use]
    pub const fn report_root(&self) -> ContentHash {
        self.report_root
    }

    /// Exact execution identity bound by the report.
    #[must_use]
    pub const fn execution(&self) -> &DiffRealExecutionIdentity {
        self.execution
    }

    /// Report receipt schema version.
    #[must_use]
    pub const fn report_schema_version(&self) -> u32 {
        self.report_schema_version
    }

    /// Stage receipt schema version.
    #[must_use]
    pub const fn stage_schema_version(&self) -> u32 {
        self.stage_schema_version
    }

    /// Exact local policy fingerprint the verifier must echo in its decision.
    #[must_use]
    pub const fn policy_fingerprint(&self) -> ContentHash {
        self.policy_fingerprint
    }

    /// Domain-separated subject that an authority must authenticate. It binds
    /// the purpose, ordered report root, full execution identity, both receipt
    /// versions, and the exact policy fingerprint.
    #[must_use]
    pub const fn subject(&self) -> ContentHash {
        self.subject
    }
}

/// Atomic authority verdict. Non-authorizing causes remain distinguishable.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PromotionVerdict {
    /// The detached attestation authenticates this exact request.
    Authorized,
    /// Authentication bytes did not verify.
    WrongSignature,
    /// The key identity is not known to the authority.
    UnknownKey,
    /// The key is known but no longer authorized.
    RevokedKey,
    /// The authority evaluated a different policy fingerprint.
    PolicyMismatch,
}

/// One indivisible verifier decision and the policy under which it was made.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PromotionVerificationDecision {
    verdict: PromotionVerdict,
    policy_fingerprint: ContentHash,
}

impl PromotionVerificationDecision {
    /// Construct an atomic verifier decision.
    #[must_use]
    pub const fn new(verdict: PromotionVerdict, policy_fingerprint: ContentHash) -> Self {
        Self {
            verdict,
            policy_fingerprint,
        }
    }

    /// Authority verdict.
    #[must_use]
    pub const fn verdict(self) -> PromotionVerdict {
        self.verdict
    }

    /// Policy fingerprint used by the authority.
    #[must_use]
    pub const fn policy_fingerprint(self) -> ContentHash {
        self.policy_fingerprint
    }
}

/// Injected capability that authenticates one ordered DiffReal report root.
pub trait PromotionReceiptVerifier {
    /// Verify the detached attestation against [`PromotionVerificationRequest::subject`],
    /// never against the raw report root alone.
    fn verify(
        &self,
        request: &PromotionVerificationRequest<'_>,
        attestation: &PromotionAttestation,
    ) -> PromotionVerificationDecision;
}

/// Default deny-all authority: local content hashes never self-promote.
#[derive(Debug, Default, Clone, Copy)]
pub struct NoPromotionReceiptVerifier;

impl PromotionReceiptVerifier for NoPromotionReceiptVerifier {
    fn verify(
        &self,
        request: &PromotionVerificationRequest<'_>,
        _attestation: &PromotionAttestation,
    ) -> PromotionVerificationDecision {
        PromotionVerificationDecision::new(PromotionVerdict::UnknownKey, request.policy_fingerprint)
    }
}

/// Fail-closed report authentication error.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PromotionReceiptError {
    /// The report uses an unknown ordered-root schema.
    UnknownReportVersion {
        /// Version carried by the rejected report.
        observed: u32,
    },
    /// One stage uses an unknown receipt schema.
    UnknownStageVersion {
        /// Stable stage name.
        stage: &'static str,
        /// Version carried by the rejected stage receipt.
        observed: u32,
    },
    /// A stage/result/root mutation or ordered-root mismatch was detected.
    IntegrityMismatch,
    /// The expected replay Cx does not match the receipt-bound provenance.
    ReplayContextMismatch,
    /// Required stage structure or semantic result evidence is malformed.
    InvalidStageSemantics,
    /// The attestation claims a different authority policy.
    AttestationPolicyMismatch,
    /// The verifier decision was made under a different policy.
    DecisionPolicyMismatch,
    /// The verifier explicitly refused authority.
    Unauthorized {
        /// Atomic refusal returned by the verifier.
        verdict: PromotionVerdict,
    },
}

impl core::fmt::Display for PromotionReceiptError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::UnknownReportVersion { observed } => {
                write!(
                    formatter,
                    "unknown DiffReal report receipt version {observed}"
                )
            }
            Self::UnknownStageVersion { stage, observed } => write!(
                formatter,
                "stage '{stage}' uses unknown DiffReal receipt version {observed}"
            ),
            Self::IntegrityMismatch => {
                formatter.write_str("DiffReal stage or ordered report receipt integrity failed")
            }
            Self::ReplayContextMismatch => formatter.write_str(
                "DiffReal report was replayed under a different execution identity or budget",
            ),
            Self::InvalidStageSemantics => formatter
                .write_str("DiffReal required-stage structure or semantic transcript is invalid"),
            Self::AttestationPolicyMismatch => formatter
                .write_str("DiffReal promotion attestation names a different policy fingerprint"),
            Self::DecisionPolicyMismatch => formatter.write_str(
                "DiffReal promotion verifier decided under a different policy fingerprint",
            ),
            Self::Unauthorized { verdict } => {
                write!(
                    formatter,
                    "DiffReal promotion authority refused: {verdict:?}"
                )
            }
        }
    }
}

impl std::error::Error for PromotionReceiptError {}

/// The full crate-authored Layer-3 battery report.
///
/// Construction is intentionally private: downstream callers may inspect the
/// stage diagnostics, but cannot assemble caller-supplied rows into a report
/// whose battery-local readiness predicates pass.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DiffRealReport {
    /// Stage logs. The four required stages have a fixed relative order;
    /// additional stages must be explicitly optional.
    stages: Vec<StageLog>,
    receipts: Vec<StageReceipt>,
    execution: DiffRealExecutionIdentity,
    receipt_version: u32,
    receipt_root: ContentHash,
}

/// A report whose ordered root was accepted by an injected authority under
/// the exact local policy fingerprint. This wrapper is the sole promotion
/// authority surface.
#[derive(Debug)]
pub struct AuthenticatedDiffRealReport<'a> {
    report: &'a DiffRealReport,
    attestation: PromotionAttestation,
    decision: PromotionVerificationDecision,
}

impl DiffRealReport {
    fn seal(stages: Vec<StageExecution>, execution: DiffRealExecutionIdentity) -> Self {
        let receipts: Vec<_> = stages
            .iter()
            .map(|stage| StageReceipt::seal(stage, &execution))
            .collect();
        let receipt_root =
            report_receipt_identity(REPORT_RECEIPT_POLICY_VERSION, &execution, &receipts);
        Self {
            stages: stages.into_iter().map(|stage| stage.log).collect(),
            receipts,
            execution,
            receipt_version: REPORT_RECEIPT_POLICY_VERSION,
            receipt_root,
        }
    }

    /// Ordered read-only stage diagnostics produced by this battery run.
    #[must_use]
    pub fn stages(&self) -> &[StageLog] {
        &self.stages
    }

    /// Ordered stage receipts bound into the report root.
    #[must_use]
    pub fn receipts(&self) -> &[StageReceipt] {
        &self.receipts
    }

    /// Receipt for a named stage.
    #[must_use]
    pub fn receipt(&self, name: &str) -> Option<&StageReceipt> {
        self.receipts.iter().find(|receipt| receipt.stage == name)
    }

    /// Execution provenance authenticated by every receipt.
    #[must_use]
    pub const fn execution_identity(&self) -> DiffRealExecutionIdentity {
        self.execution
    }

    /// Ordered report content root.
    #[must_use]
    pub const fn receipt_root(&self) -> ContentHash {
        self.receipt_root
    }

    /// Domain-separated subject for an external promotion attestation. This is
    /// not authority by itself; it is the exact byte-identity a verifier must
    /// authenticate.
    #[must_use]
    pub fn promotion_verification_subject(&self) -> ContentHash {
        promotion_verification_subject(
            self.receipt_root,
            &self.execution,
            self.receipt_version,
            STAGE_RECEIPT_POLICY_VERSION,
            promotion_policy_fingerprint(),
        )
    }

    /// Recompute every stage receipt and the ordered report root.
    #[must_use]
    pub fn verifies_integrity(&self) -> bool {
        self.receipt_version == REPORT_RECEIPT_POLICY_VERSION
            && self.receipts.len() == self.stages.len()
            && self
                .receipts
                .iter()
                .zip(&self.stages)
                .all(|(receipt, stage)| receipt.verifies(stage, &self.execution))
            && self.receipt_root
                == report_receipt_identity(self.receipt_version, &self.execution, &self.receipts)
    }

    /// Verify receipt integrity and reject replay under a different Cx stream,
    /// mode, deadline, poll quota, cost quota, or priority.
    #[must_use]
    pub fn verifies_for(&self, cx: &Cx<'_>) -> bool {
        self.execution.matches_cx(cx) && self.verifies_integrity()
    }

    /// Authenticate this exact ordered report root under an injected external
    /// authority. Local integrity, replay context, versions, and semantic
    /// structure are checked before the verifier is invoked.
    ///
    /// # Errors
    /// Returns a precise fail-closed cause for unknown schemas, mutation,
    /// replay mismatch, semantic contradiction, policy drift, or refusal.
    pub fn authenticate<'a>(
        &'a self,
        expected_cx: &Cx<'_>,
        attestation: PromotionAttestation,
        verifier: &dyn PromotionReceiptVerifier,
    ) -> Result<AuthenticatedDiffRealReport<'a>, PromotionReceiptError> {
        if self.receipt_version != REPORT_RECEIPT_POLICY_VERSION {
            return Err(PromotionReceiptError::UnknownReportVersion {
                observed: self.receipt_version,
            });
        }
        if let Some(receipt) = self
            .receipts
            .iter()
            .find(|receipt| receipt.version != STAGE_RECEIPT_POLICY_VERSION)
        {
            return Err(PromotionReceiptError::UnknownStageVersion {
                stage: receipt.stage,
                observed: receipt.version,
            });
        }
        if !self.verifies_integrity() {
            return Err(PromotionReceiptError::IntegrityMismatch);
        }
        if !self.execution.matches_cx(expected_cx) {
            return Err(PromotionReceiptError::ReplayContextMismatch);
        }
        if !self.required_schema_is_valid() || !self.stages.iter().all(stage_semantics_are_valid) {
            return Err(PromotionReceiptError::InvalidStageSemantics);
        }

        let policy_fingerprint = promotion_policy_fingerprint();
        if attestation.policy_fingerprint != policy_fingerprint {
            return Err(PromotionReceiptError::AttestationPolicyMismatch);
        }
        let subject = promotion_verification_subject(
            self.receipt_root,
            &self.execution,
            self.receipt_version,
            STAGE_RECEIPT_POLICY_VERSION,
            policy_fingerprint,
        );
        let request = PromotionVerificationRequest {
            purpose: PROMOTION_VERIFICATION_PURPOSE,
            report_root: self.receipt_root,
            execution: &self.execution,
            report_schema_version: self.receipt_version,
            stage_schema_version: STAGE_RECEIPT_POLICY_VERSION,
            policy_fingerprint,
            subject,
        };
        let decision = verifier.verify(&request, &attestation);
        if decision.policy_fingerprint != policy_fingerprint {
            return Err(PromotionReceiptError::DecisionPolicyMismatch);
        }
        if decision.verdict != PromotionVerdict::Authorized {
            return Err(PromotionReceiptError::Unauthorized {
                verdict: decision.verdict,
            });
        }
        Ok(AuthenticatedDiffRealReport {
            report: self,
            attestation,
            decision,
        })
    }

    /// Is every required stage present exactly once with the expected evidence
    /// identity and an evaluated (`Passed` or `Failed`) result?
    #[must_use]
    pub fn complete(&self) -> bool {
        self.verifies_integrity()
            && self.required_schema_is_valid()
            && REQUIRED_STAGES.iter().all(|required| {
                self.stage(required.name)
                    .is_some_and(|stage| stage.status.is_evaluated())
            })
    }

    /// Did every required stage actually run and pass?
    ///
    /// Missing, duplicated, gated, refused, identity-mismatched, or malformed
    /// required records all return `false`.
    #[must_use]
    pub fn all_required_passed(&self) -> bool {
        self.verifies_integrity()
            && self.required_schema_is_valid()
            && REQUIRED_STAGES
                .iter()
                .all(|required| self.stage(required.name).is_some_and(StageLog::passed))
    }

    /// Is the untrusted report structurally and semantically ready for an
    /// authority decision? This is not promotion authority.
    #[must_use]
    pub fn structurally_ready(&self) -> bool {
        self.complete() && self.all_required_passed()
    }

    /// A named stage.
    #[must_use]
    pub fn stage(&self, name: &str) -> Option<&StageLog> {
        self.stages.iter().find(|s| s.stage == name)
    }

    fn required_schema_is_valid(&self) -> bool {
        if self.stages.iter().any(|stage| !stage.is_well_formed())
            || self.stages.iter().enumerate().any(|(index, stage)| {
                self.stages[index + 1..]
                    .iter()
                    .any(|other| stage.stage == other.stage)
            })
        {
            return false;
        }

        if self
            .stages
            .iter()
            .filter(|stage| stage.requirement == StageRequirement::Required)
            .count()
            != REQUIRED_STAGES.len()
        {
            return false;
        }

        self.stages
            .iter()
            .filter(|stage| stage.requirement == StageRequirement::Required)
            .zip(REQUIRED_STAGES.iter())
            .all(|(stage, required)| {
                stage.stage == required.name
                    && stage.evidence_identity == required.evidence_identity
            })
    }
}

impl AuthenticatedDiffRealReport<'_> {
    /// The authenticated report root.
    #[must_use]
    pub const fn receipt_root(&self) -> ContentHash {
        self.report.receipt_root
    }

    /// Detached authority attestation retained with this decision.
    #[must_use]
    pub const fn attestation(&self) -> &PromotionAttestation {
        &self.attestation
    }

    /// Atomic authority decision retained with this report.
    #[must_use]
    pub const fn decision(&self) -> PromotionVerificationDecision {
        self.decision
    }

    /// Did an authenticated report also satisfy every required scientific
    /// stage? Only this opaque wrapper exposes promotion readiness.
    #[must_use]
    pub fn promotion_ready(&self) -> bool {
        self.report.structurally_ready()
    }

    /// Read-only diagnostic report.
    #[must_use]
    pub const fn report(&self) -> &DiffRealReport {
        self.report
    }
}

#[derive(Clone, Copy)]
struct RequiredStage {
    name: &'static str,
    evidence_identity: &'static str,
}

const REQUIRED_STAGES: [RequiredStage; 4] = [
    RequiredStage {
        name: DIFFERENTIATION_STAGE,
        evidence_identity: DIFFERENTIATION_EVIDENCE_IDENTITY,
    },
    RequiredStage {
        name: AS_BUILT_STAGE,
        evidence_identity: AS_BUILT_EVIDENCE_IDENTITY,
    },
    RequiredStage {
        name: TOLERANCE_STAGE,
        evidence_identity: TOLERANCE_EVIDENCE_IDENTITY,
    },
    RequiredStage {
        name: SPACETIME_STAGE,
        evidence_identity: SPACETIME_EVIDENCE_IDENTITY,
    },
];

fn gradient_event_is_valid(event: &StageEvent, expected_input: f64) -> bool {
    let StageEvent::GradientVerified {
        receipt,
        input_bits,
        value_bits,
        production_bits,
        dual_bits,
        fd_coarse_bits,
        fd_fine_bits,
        tolerance_bits,
    } = event
    else {
        return false;
    };
    let bits = [
        *input_bits,
        *value_bits,
        *production_bits,
        *dual_bits,
        *fd_coarse_bits,
        *fd_fine_bits,
        *tolerance_bits,
    ];
    if *input_bits != expected_input.to_bits()
        || bits.iter().any(|bits| !f64::from_bits(*bits).is_finite())
        || f64::from_bits(*tolerance_bits) < 0.0
    {
        return false;
    }
    let ops: Vec<String> = PRODUCTION_DIFFERENTIATION_PATH
        .iter()
        .map(|op| (*op).to_string())
        .collect();
    let tolerance = f64::from_bits(*tolerance_bits);
    let production = f64::from_bits(*production_bits);
    let dual = f64::from_bits(*dual_bits);
    let fd_coarse = f64::from_bits(*fd_coarse_bits);
    let fd_fine = f64::from_bits(*fd_fine_bits);
    *receipt == sensitivity_identity(&ops, bits)
        && (production - dual).abs() <= tolerance
        && (production - fd_fine).abs() <= tolerance
        && 3.0 * (fd_coarse - fd_fine).abs() <= tolerance
}

fn differentiation_passed_semantics(stage: &StageLog) -> bool {
    stage.events.len() == 2
        && gradient_event_is_valid(&stage.events[0], DIFFERENTIATION_FIXTURE_INPUT)
        && matches!(
            &stage.events[1],
            StageEvent::MissingVjpProbe { op, blocked: true } if op == "remesh"
        )
}

fn as_built_passed_semantics(stage: &StageLog) -> bool {
    if stage.events.len() != 3 {
        return false;
    }
    matches!(
        &stage.events[0],
        StageEvent::Registration {
            residual_bits,
            within_tolerance: true,
        } if f64::from_bits(*residual_bits).is_finite()
    ) && matches!(
        &stage.events[1],
        StageEvent::AsBuiltDelta {
            max_deviation_bits,
            defect_index: Some(AS_BUILT_DEFECT_INDEX),
            estimated: true,
        } if f64::from_bits(*max_deviation_bits).is_finite()
    ) && matches!(
        &stage.events[2],
        StageEvent::Assimilation {
            before_bits,
            after_bits,
            reduced: true,
        } if f64::from_bits(*before_bits).is_finite()
            && f64::from_bits(*after_bits).is_finite()
            && f64::from_bits(*after_bits) <= f64::from_bits(*before_bits)
    )
}

fn tolerance_passed_semantics(stage: &StageLog) -> bool {
    if stage.events.len() != 5 {
        return false;
    }
    gradient_event_is_valid(&stage.events[0], TOLERANCE_SENSITIVITY_INPUTS[0])
        && gradient_event_is_valid(&stage.events[1], TOLERANCE_SENSITIVITY_INPUTS[1])
        && matches!(
            &stage.events[2],
            StageEvent::ToleranceActions {
                critical: Some(Action::Tighten),
                slack: Some(Action::Loosen),
            }
        )
        && matches!(
            &stage.events[3],
            StageEvent::GdtJustification {
                loosened,
                all_verified: true,
            } if *loosened > 0
        )
        && matches!(
            &stage.events[4],
            StageEvent::SampledLinearization {
                samples,
                confirmed: true,
                linearized_std_bits,
                probability_claimed: false,
            } if *samples == TOLERANCE_EXTREME_QOIS.len()
                && f64::from_bits(*linearized_std_bits).is_finite()
                && f64::from_bits(*linearized_std_bits) >= 0.0
        )
}

fn passed_stage_semantics(stage: &StageLog) -> bool {
    match stage.stage {
        DIFFERENTIATION_STAGE => differentiation_passed_semantics(stage),
        AS_BUILT_STAGE => as_built_passed_semantics(stage),
        TOLERANCE_STAGE => tolerance_passed_semantics(stage),
        // There is no activated spacetime success transcript in receipt v1.
        SPACETIME_STAGE => false,
        _ => stage.requirement == StageRequirement::Optional,
    }
}

fn stage_semantics_are_valid(stage: &StageLog) -> bool {
    if !stage.is_well_formed() {
        return false;
    }
    match &stage.status {
        StageStatus::Passed => passed_stage_semantics(stage),
        StageStatus::Failed(_) => !passed_stage_semantics(stage),
        StageStatus::Gated(reason) => {
            stage.stage == SPACETIME_STAGE
                && reason.code == "diffreal.spacetime.integration-not-activated"
                && stage.events.iter().any(|event| {
                    matches!(
                        event,
                        StageEvent::Gate { code, .. }
                            if *code == "diffreal.spacetime.integration-not-activated"
                    )
                })
        }
        StageStatus::Refused(reason) => stage
            .events
            .iter()
            .any(|event| matches!(event, StageEvent::Refusal { code, .. } if *code == reason.code)),
    }
}

/// Typed refusal that prevents publication of a partial battery report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DiffRealError {
    /// Invocation preflight, deadline, cancellation, or affine accounting
    /// refused the composed scientific transaction.
    Invocation(fs_exec::InvocationError),
    /// A drained transaction produced a non-acceptance terminal receipt.
    InvocationDidNotComplete {
        /// Derived terminal disposition.
        disposition: fs_exec::InvocationDisposition,
        /// Canonical terminal accounting root.
        receipt_root: ContentHash,
        /// Full terminal accounting evidence retained after drain.
        receipt: Box<fs_exec::InvocationReceipt>,
        /// Scientific/runtime cause observed before drain, when one exists.
        cause: Option<Box<DiffRealError>>,
    },
    /// Differentiation admission or runtime work could not produce a stage
    /// disposition.
    Differentiation(DifferentiationError),
    /// Registration, registered geometry, or as-built comparison failed.
    AsBuilt(fs_asbuilt::RegError),
    /// Belief construction, observation declaration, or assimilation failed.
    Assimilation(AssimError),
    /// A stage observed cancellation at a bounded boundary.
    Cancelled {
        /// Stable stage name.
        stage: &'static str,
    },
    /// A fixed stage plan exceeded the ambient cost quota before work began.
    WorkBudgetExceeded {
        /// Stable stage name.
        stage: &'static str,
        /// Fixed required work units.
        required: u64,
        /// Ambient available work units.
        available: u64,
    },
}

impl core::fmt::Display for DiffRealError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invocation(error) => write!(formatter, "invocation refused: {error}"),
            Self::InvocationDidNotComplete {
                disposition,
                receipt_root,
                cause,
                ..
            } => {
                write!(
                    formatter,
                    "invocation drained with {disposition:?} disposition (receipt {receipt_root})"
                )?;
                if let Some(cause) = cause {
                    write!(formatter, ": {cause}")?;
                }
                Ok(())
            }
            Self::Differentiation(error) => {
                write!(formatter, "differentiation stage failed: {error}")
            }
            Self::AsBuilt(error) => write!(formatter, "as-built stage failed: {error}"),
            Self::Assimilation(error) => write!(formatter, "assimilation stage failed: {error}"),
            Self::Cancelled { stage } => {
                write!(
                    formatter,
                    "stage '{stage}' cancelled at a bounded checkpoint"
                )
            }
            Self::WorkBudgetExceeded {
                stage,
                required,
                available,
            } => write!(
                formatter,
                "stage '{stage}' requires {required} work units but the context admits {available}"
            ),
        }
    }
}

impl std::error::Error for DiffRealError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Invocation(error) => Some(error),
            Self::Differentiation(error) => Some(error),
            Self::AsBuilt(error) => Some(error),
            Self::Assimilation(error) => Some(error),
            Self::InvocationDidNotComplete { cause, .. } => cause
                .as_deref()
                .map(|cause| cause as &(dyn std::error::Error + 'static)),
            Self::Cancelled { .. } | Self::WorkBudgetExceeded { .. } => None,
        }
    }
}

impl DiffRealError {
    /// Terminal invocation evidence for runtime failures that reached drain.
    #[must_use]
    pub fn terminal_invocation_receipt(&self) -> Option<&fs_exec::InvocationReceipt> {
        match self {
            Self::InvocationDidNotComplete { receipt, .. } => Some(receipt),
            _ => None,
        }
    }
}

impl From<fs_asbuilt::RegError> for DiffRealError {
    fn from(error: fs_asbuilt::RegError) -> Self {
        Self::AsBuilt(error)
    }
}

impl From<fs_exec::InvocationError> for DiffRealError {
    fn from(error: fs_exec::InvocationError) -> Self {
        Self::Invocation(error)
    }
}

impl From<DifferentiationError> for DiffRealError {
    fn from(error: DifferentiationError) -> Self {
        Self::Differentiation(error)
    }
}

impl From<AssimError> for DiffRealError {
    fn from(error: AssimError) -> Self {
        Self::Assimilation(error)
    }
}

/// Run the full Layer-3 battery.
///
/// # Errors
/// Propagates structured cancellation, ambient-budget refusal, or a lower-layer
/// as-built/assimilation error. No partial battery report is published.
pub fn run_battery(cx: &Cx<'_>) -> Result<DiffRealReport, DiffRealError> {
    let clock = fs_exec::WallClock::new();
    run_battery_with_clock(cx, &clock)
}

/// Run the full battery with an injected monotonic clock for deterministic
/// deadline-edge and replay verification.
///
/// # Errors
/// As [`run_battery`].
pub fn run_battery_with_clock(
    cx: &Cx<'_>,
    clock: &dyn fs_exec::TimeSource,
) -> Result<DiffRealReport, DiffRealError> {
    let execution = DiffRealExecutionIdentity::from_cx(cx);
    Ok(DiffRealReport::seal(
        vec![
            execute_differentiation(cx)?,
            execute_as_built_loop(cx, clock)?,
            execute_tolerance_allocation(cx)?,
            execute_spacetime_gated(cx)?,
        ],
        execution,
    ))
}

// -- Stage 1: differentiation ----------------------------------------------

/// The fixture composite `f(x) = (2x + 1)²` used only by the independent
/// finite-difference oracle.
fn composite(x: f64) -> f64 {
    let h = 2.0 * x + 1.0;
    h * h
}

fn admit_op_name(index: usize, op: &str) -> Result<(), DifferentiationError> {
    if op.is_empty() {
        return Err(DifferentiationError::EmptyOpName { index });
    }
    if op.len() > MAX_OP_NAME_BYTES {
        return Err(DifferentiationError::OpNameTooLong {
            index,
            limit: MAX_OP_NAME_BYTES,
            observed: op.len(),
        });
    }
    Ok(())
}

fn admit_path(ops: &[&str]) -> Result<(), DifferentiationError> {
    if ops.is_empty() {
        return Err(DifferentiationError::EmptyPath);
    }
    if ops.len() > MAX_DIFFERENTIATION_OPS {
        return Err(DifferentiationError::PathTooLong {
            limit: MAX_DIFFERENTIATION_OPS,
            observed: ops.len(),
        });
    }
    for (index, op) in ops.iter().enumerate() {
        admit_op_name(index, op)?;
    }
    Ok(())
}

fn differentiation_checkpoint(cx: &Cx<'_>) -> Result<(), DifferentiationError> {
    cx.checkpoint().map_err(|_| DifferentiationError::Cancelled)
}

fn admit_differentiation_work(cx: &Cx<'_>, required: u64) -> Result<(), DifferentiationError> {
    differentiation_checkpoint(cx)?;
    if let Some(available) = cx.budget().cost_quota
        && available < required
    {
        return Err(DifferentiationError::WorkBudgetExceeded {
            required,
            available,
        });
    }
    Ok(())
}

fn admit_stage_work(cx: &Cx<'_>, stage: &'static str, required: u64) -> Result<(), DiffRealError> {
    cx.checkpoint()
        .map_err(|_| DiffRealError::Cancelled { stage })?;
    if let Some(available) = cx.budget().cost_quota
        && available < required
    {
        return Err(DiffRealError::WorkBudgetExceeded {
            stage,
            required,
            available,
        });
    }
    Ok(())
}

fn stage_checkpoint(cx: &Cx<'_>, stage: &'static str) -> Result<(), DiffRealError> {
    cx.checkpoint()
        .map_err(|_| DiffRealError::Cancelled { stage })
}

fn stage_runtime_refusal(stage: &'static str, error: DifferentiationError) -> DiffRealError {
    match error {
        DifferentiationError::Cancelled => DiffRealError::Cancelled { stage },
        DifferentiationError::WorkBudgetExceeded {
            required,
            available,
        } => DiffRealError::WorkBudgetExceeded {
            stage,
            required,
            available,
        },
        other => DiffRealError::Differentiation(other),
    }
}

fn scalar_primal(primal_inputs: &[&[f64]]) -> f64 {
    primal_inputs
        .first()
        .and_then(|values| values.first())
        .copied()
        .unwrap_or(f64::NAN)
}

fn scalar_cotangent(out_cotangent: &[f64]) -> f64 {
    out_cotangent.first().copied().unwrap_or(f64::NAN)
}

#[derive(Debug)]
struct AffineVjp;

impl Vjp for AffineVjp {
    fn vjp(&self, _primal_inputs: &[&[f64]], out_cotangent: &[f64]) -> Vec<Vec<f64>> {
        vec![vec![2.0 * scalar_cotangent(out_cotangent)]]
    }
}

#[derive(Debug)]
struct SquareVjp;

impl Vjp for SquareVjp {
    fn vjp(&self, primal_inputs: &[&[f64]], out_cotangent: &[f64]) -> Vec<Vec<f64>> {
        vec![vec![
            (2.0 * scalar_primal(primal_inputs)) * scalar_cotangent(out_cotangent),
        ]]
    }
}

#[derive(Debug)]
struct IdentityVjp;

impl Vjp for IdentityVjp {
    fn vjp(&self, _primal_inputs: &[&[f64]], out_cotangent: &[f64]) -> Vec<Vec<f64>> {
        vec![vec![scalar_cotangent(out_cotangent)]]
    }
}

/// Construct the production fixture's shared tape/VJP registry.
///
/// # Errors
/// Returns a structured name-admission error if a future fixed operator name
/// violates the public registry bounds.
pub fn production_vjp_registry() -> Result<DifferentiationRegistry, DifferentiationError> {
    let mut registry = DifferentiationRegistry::new();
    registry.register("sdf", AffineVjp)?;
    registry.register("spline", SquareVjp)?;
    registry.register("solve", IdentityVjp)?;
    Ok(registry)
}

fn apply_fixture_operator(op: &str, input: f64) -> Result<f64, DifferentiationError> {
    match op {
        "sdf" => Ok(2.0 * input + 1.0),
        "spline" => Ok(input * input),
        "solve" => Ok(input),
        _ => Err(DifferentiationError::UnsupportedOperator { op: op.to_string() }),
    }
}

/// Differentiate an admitted scalar operator path through the shared
/// `fs-adjoint` tape and VJP transpose. Missing registrations are checked in
/// forward order before input arithmetic, so they take precedence over a
/// hostile non-finite value.
///
/// # Errors
/// Returns a typed admission, cancellation, forward, transpose, or
/// representability refusal. It never substitutes a silent zero.
pub fn differentiate_path(
    ops: &[&str],
    registry: &DifferentiationRegistry,
    x: f64,
    cx: &Cx<'_>,
) -> Result<PathDerivative, DifferentiationError> {
    differentiation_checkpoint(cx)?;
    admit_path(ops)?;
    for op in ops {
        if !registry.registered.contains(*op) {
            return Err(DifferentiationError::MissingVjp {
                op: (*op).to_string(),
            });
        }
    }
    if !x.is_finite() {
        return Err(DifferentiationError::NonFiniteInput { bits: x.to_bits() });
    }
    let path_units = u64::try_from(ops.len())
        .unwrap_or(u64::MAX)
        .saturating_add(2);
    admit_differentiation_work(cx, path_units)?;

    let mut tape = Tape::new();
    let leaf = tape.leaf(vec![x]);
    let mut current = leaf;
    let mut value = x;
    for op in ops {
        differentiation_checkpoint(cx)?;
        value = apply_fixture_operator(op, value)?;
        if !value.is_finite() {
            return Err(DifferentiationError::NonFinitePrimal {
                op: (*op).to_string(),
                bits: value.to_bits(),
            });
        }
        current = tape.apply(op, &[current], vec![value]);
    }

    differentiation_checkpoint(cx)?;
    let gradients = tape
        .transpose(&registry.inner, current, &[1.0])
        .map_err(DifferentiationError::Transpose)?;
    let input_gradient = gradients
        .get(&leaf)
        .ok_or(DifferentiationError::MissingLeafGradient)?;
    if input_gradient.len() != 1 {
        return Err(DifferentiationError::InvalidGradientShape {
            observed: input_gradient.len(),
        });
    }
    let gradient = input_gradient[0];
    if !gradient.is_finite() {
        return Err(DifferentiationError::NonFiniteGradient {
            bits: gradient.to_bits(),
        });
    }
    differentiation_checkpoint(cx)?;
    Ok(PathDerivative {
        value_bits: value.to_bits(),
        gradient_bits: gradient.to_bits(),
    })
}

fn dual_fixture(x: f64) -> (f64, f64) {
    let (value, [gradient]) = dual_gradient([x], |[x]| {
        let affine = x * Dual64::constant(2.0) + Dual64::constant(1.0);
        affine * affine
    });
    (value, gradient)
}

fn sensitivity_identity(ops: &[String], fields: [u64; 7]) -> ContentHash {
    let policy = SENSITIVITY_POLICY_VERSION.as_bytes();
    let op_bytes: usize = ops.iter().map(String::len).sum();
    let mut preimage = Vec::with_capacity(
        8 + policy.len() + 8 + ops.len().saturating_mul(8) + op_bytes + fields.len() * 8,
    );
    preimage.extend_from_slice(
        &u64::try_from(policy.len())
            .unwrap_or(u64::MAX)
            .to_le_bytes(),
    );
    preimage.extend_from_slice(policy);
    preimage.extend_from_slice(&u64::try_from(ops.len()).unwrap_or(u64::MAX).to_le_bytes());
    for op in ops {
        preimage.extend_from_slice(&u64::try_from(op.len()).unwrap_or(u64::MAX).to_le_bytes());
        preimage.extend_from_slice(op.as_bytes());
    }
    for field in fields {
        preimage.extend_from_slice(&field.to_le_bytes());
    }
    hash_domain(SENSITIVITY_IDENTITY_DOMAIN, &preimage)
}

/// Run the production reverse sweep for [`PRODUCTION_DIFFERENTIATION_PATH`] and
/// independently seal its result with a dual-number oracle plus a two-step
/// conditioning-aware FD falsifier. Other paths are rejected because these
/// fixture oracles do not define their semantics.
///
/// # Errors
/// Returns a typed path, runtime, representability, or oracle disagreement.
pub fn verify_sensitivity(
    ops: &[&str],
    registry: &DifferentiationRegistry,
    x: f64,
    cx: &Cx<'_>,
) -> Result<SealedSensitivity, DifferentiationError> {
    differentiation_checkpoint(cx)?;
    admit_path(ops)?;
    let first_mismatch = ops
        .iter()
        .zip(PRODUCTION_DIFFERENTIATION_PATH.iter())
        .position(|(observed, expected)| observed != expected)
        .unwrap_or_else(|| ops.len().min(PRODUCTION_DIFFERENTIATION_PATH.len()));
    if ops != PRODUCTION_DIFFERENTIATION_PATH.as_slice() {
        return Err(DifferentiationError::OraclePathMismatch {
            expected_len: PRODUCTION_DIFFERENTIATION_PATH.len(),
            observed_len: ops.len(),
            first_mismatch,
        });
    }
    admit_differentiation_work(cx, DIFFERENTIATION_WORK_UNITS)?;
    let production = differentiate_path(ops, registry, x, cx)?;
    differentiation_checkpoint(cx)?;

    let (dual_value, dual) = dual_fixture(x);
    // FD step choice (7namd repair): the production composite is a
    // QUADRATIC, so central differences carry ZERO truncation error at
    // any step — a larger step strictly shrinks the rounding noise of
    // the difference quotient. At h = 1e-4 the noise near the
    // composite's root (x ≈ -0.5, where 2x+1 cancels catastrophically)
    // reaches ~2e-14 and crosses the 64ε acceptance band; h = 1e-2
    // keeps the same exact analytic derivative with ~100× less noise.
    let h = 1.0e-2 * x.abs().max(1.0);
    // 1024ε: after the libm-doctrine routing (2b2cbc4) the composite's
    // det-routed transcendentals shift at the last ULP, and the central
    // difference quotient amplifies that to ~1e-13 near the cancellation
    // root — above the old 64ε band. This band is a ROUNDING-NOISE
    // budget, not a correctness tolerance: the fixture is quadratic
    // (zero truncation), so a genuinely wrong gradient still misses by
    // O(h) ≈ 1e-2, five orders above the band; falsifier power is
    // unaffected.
    let fd = fd_falsifier(
        &|point| composite(point[0]),
        &[x],
        &[1.0],
        production.gradient(),
        h,
        1024.0 * f64::EPSILON,
    );
    let value_tolerance =
        64.0 * f64::EPSILON * production.value().abs().max(dual_value.abs()).max(1.0);
    let gradient_tolerance =
        64.0 * f64::EPSILON * production.gradient().abs().max(dual.abs()).max(1.0);
    let tolerance = fd.tolerance.max(gradient_tolerance);
    let oracles_finite = dual_value.is_finite()
        && dual.is_finite()
        && fd.fd_coarse.is_finite()
        && fd.fd_fine.is_finite()
        && tolerance.is_finite();
    let values_agree = (production.value() - dual_value).abs() <= value_tolerance;
    let dual_agrees = (production.gradient() - dual).abs() <= gradient_tolerance;
    if !oracles_finite || !values_agree || !dual_agrees || !fd.consistent {
        return Err(DifferentiationError::OracleDisagreement {
            production_bits: production.gradient().to_bits(),
            dual_bits: dual.to_bits(),
            fd_fine_bits: fd.fd_fine.to_bits(),
            tolerance_bits: tolerance.to_bits(),
        });
    }
    differentiation_checkpoint(cx)?;

    let ops: Vec<String> = ops.iter().map(|op| (*op).to_string()).collect();
    let fields = [
        x.to_bits(),
        production.value().to_bits(),
        production.gradient().to_bits(),
        dual.to_bits(),
        fd.fd_coarse.to_bits(),
        fd.fd_fine.to_bits(),
        tolerance.to_bits(),
    ];
    let identity = sensitivity_identity(&ops, fields);
    Ok(SealedSensitivity {
        ops,
        input_bits: fields[0],
        value_bits: fields[1],
        production_gradient_bits: fields[2],
        dual_gradient_bits: fields[3],
        fd_coarse_bits: fields[4],
        fd_fine_bits: fields[5],
        fd_tolerance_bits: fields[6],
        identity,
    })
}

fn sensitivity_event(receipt: &SealedSensitivity) -> StageEvent {
    StageEvent::GradientVerified {
        receipt: receipt.identity(),
        input_bits: receipt.input_bits,
        value_bits: receipt.value_bits,
        production_bits: receipt.production_gradient_bits,
        dual_bits: receipt.dual_gradient_bits,
        fd_coarse_bits: receipt.fd_coarse_bits,
        fd_fine_bits: receipt.fd_fine_bits,
        tolerance_bits: receipt.fd_tolerance_bits,
    }
}

/// Stage 1 with an injectable registry for independent kill tests.
///
/// # Errors
/// Cancellation and ambient work-budget refusals suppress the partial report.
fn execute_differentiation_with_registry(
    cx: &Cx<'_>,
    registry: &DifferentiationRegistry,
) -> Result<StageExecution, DiffRealError> {
    admit_stage_work(cx, DIFFERENTIATION_STAGE, DIFFERENTIATION_STAGE_WORK_UNITS)?;
    let sensitivity = match verify_sensitivity(
        &PRODUCTION_DIFFERENTIATION_PATH,
        registry,
        DIFFERENTIATION_FIXTURE_INPUT,
        cx,
    ) {
        Ok(sensitivity) => sensitivity,
        Err(error) if error.is_runtime_refusal() => {
            return Err(stage_runtime_refusal(DIFFERENTIATION_STAGE, error));
        }
        Err(error) => {
            return Ok(StageExecution::diagnostic(StageLog::new(
                DIFFERENTIATION_STAGE,
                StageRequirement::Required,
                StageStatus::Failed(StageReason::new(
                    "diffreal.differentiation.production-rejected",
                    "production differentiation or an independent oracle was rejected; inspect the typed event",
                )),
                DIFFERENTIATION_EVIDENCE_IDENTITY,
                vec![StageEvent::DifferentiationRejected { error }],
            )));
        }
    };

    let mut events = vec![sensitivity_event(&sensitivity)];
    let missing_probe = differentiate_path(
        &MISSING_VJP_FIXTURE_PATH,
        registry,
        DIFFERENTIATION_FIXTURE_INPUT,
        cx,
    );
    let blocked = matches!(
        &missing_probe,
        Err(DifferentiationError::MissingVjp { op }) if op == "remesh"
    );
    events.push(StageEvent::MissingVjpProbe {
        op: "remesh".to_string(),
        blocked,
    });
    if let Err(error) = &missing_probe {
        if error.is_runtime_refusal() {
            return Err(stage_runtime_refusal(DIFFERENTIATION_STAGE, error.clone()));
        }
        if !matches!(error, DifferentiationError::MissingVjp { .. }) {
            events.push(StageEvent::DifferentiationRejected {
                error: error.clone(),
            });
        }
    }

    let status = if sensitivity.verifies_integrity() && blocked {
        StageStatus::Passed
    } else {
        StageStatus::Failed(StageReason::new(
            "diffreal.differentiation.assertion-failed",
            "production gradient sealing or the missing-VJP kill assertion failed; inspect typed events",
        ))
    };
    let log = StageLog::new(
        DIFFERENTIATION_STAGE,
        StageRequirement::Required,
        status,
        DIFFERENTIATION_EVIDENCE_IDENTITY,
        events,
    );
    Ok(StageExecution {
        inputs: fixture_input_identity(&log),
        result: differentiation_result_identity(&sensitivity, &missing_probe),
        log,
        invocation: None,
    })
}

/// Stage 1 with an injectable registry for independent kill tests.
///
/// This diagnostic seam never mints a stage receipt; only the fixed
/// [`run_battery`] path can add a result to a sealed report.
///
/// # Errors
/// Cancellation and ambient work-budget refusals suppress the partial report.
pub fn stage_differentiation_with_registry(
    cx: &Cx<'_>,
    registry: &DifferentiationRegistry,
) -> Result<StageLog, DiffRealError> {
    execute_differentiation_with_registry(cx, registry).map(|stage| stage.log)
}

/// Stage 1: production tape/VJP gradient, independent dual/FD sealing, and a
/// missing-VJP kill probe.
///
/// # Errors
/// Cancellation and ambient work-budget refusals suppress the partial report.
pub fn stage_differentiation(cx: &Cx<'_>) -> Result<StageLog, DiffRealError> {
    let registry = production_vjp_registry()?;
    execute_differentiation_with_registry(cx, &registry).map(|stage| stage.log)
}

fn execute_differentiation(cx: &Cx<'_>) -> Result<StageExecution, DiffRealError> {
    let registry = production_vjp_registry()?;
    execute_differentiation_with_registry(cx, &registry)
}

// -- Stage 2: as-built loop -------------------------------------------------

/// Stage 2: register a scan, estimate as-built δ, localize a defect, assimilate.
///
/// # Errors
/// Propagates fixed-work admission, cancellation, or a structured lower-layer
/// refusal and publishes no partial stage log.
fn run_as_built_phase<T>(
    parent: &mut fs_exec::ChildBudget<'_, '_>,
    phase: &'static str,
    grant: fs_exec::InvocationResources,
    refusal_reason: &'static str,
    operation: impl FnOnce(&mut fs_exec::ChildBudget<'_, '_>) -> Result<T, DiffRealError>,
) -> Result<T, DiffRealError> {
    let mut child = parent.split_child(phase, grant)?;
    let outcome = operation(&mut child);
    if outcome.is_err() {
        child.refuse(phase, invocation::domain_refusal(phase, refusal_reason));
    }
    let disposition = child.finish()?;
    match outcome {
        Ok(value) if disposition == fs_exec::InvocationDisposition::Completed => Ok(value),
        Ok(_) => Err(DiffRealError::Cancelled {
            stage: AS_BUILT_STAGE,
        }),
        Err(error) => Err(error),
    }
}

fn as_built_observations() -> Result<Vec<Observation>, DiffRealError> {
    AS_BUILT_OBSERVATIONS
        .iter()
        .map(|&(state_index, value, variance, label)| {
            point_sensor(
                state_index,
                AS_BUILT_PRIOR_MEAN.len(),
                value,
                variance,
                label,
            )
            .map_err(DiffRealError::from)
        })
        .collect()
}

#[allow(clippy::too_many_lines)]
fn execute_as_built_science(
    cx: &Cx<'_>,
    transaction: &mut fs_exec::ChildBudget<'_, '_>,
    plan: &invocation::AsBuiltInvocationPlan,
) -> Result<(StageLog, ContentHash, ContentHash), DiffRealError> {
    transaction.charge_work(plan.setup.work())?;
    transaction.charge_cost(plan.setup.cost())?;
    transaction.charge_evaluations(plan.setup.evaluations())?;
    transaction.poll("as-built.setup.begin")?;
    let mut setup_memory =
        transaction.reserve_memory("as-built-fixed-inputs", plan.setup.memory())?;
    let mut events = Vec::new();
    let mut assertions_passed = true;

    // a scanned fixture: design datums transformed by a known rigid motion.
    let design = [
        Point2::new(AS_BUILT_DESIGN_POINTS[0].0, AS_BUILT_DESIGN_POINTS[0].1)?,
        Point2::new(AS_BUILT_DESIGN_POINTS[1].0, AS_BUILT_DESIGN_POINTS[1].1)?,
        Point2::new(AS_BUILT_DESIGN_POINTS[2].0, AS_BUILT_DESIGN_POINTS[2].1)?,
    ];
    let theta = AS_BUILT_ROTATION_RADIANS;
    let (tx, ty) = AS_BUILT_TRANSLATION;
    let xf = |p: Point2| {
        let (s, c) = theta.sin_cos();
        Point2::new(c * p.x() - s * p.y() + tx, s * p.x() + c * p.y() + ty)
    };
    let fids: Vec<Fiducial> = design
        .iter()
        .map(|&datum| Ok(Fiducial::new(datum, xf(datum)?)))
        .collect::<Result<_, fs_asbuilt::RegError>>()?;
    setup_memory.budget().poll("as-built.setup.complete")?;
    let reg = run_as_built_phase(
        setup_memory.budget(),
        "as-built.registration",
        plan.registration,
        "registration refused its fixed fixture",
        |budget| register_budgeted(&fids, budget).map_err(DiffRealError::from),
    )?;
    let reg_ok = reg.residual_rms() < 1e-9;
    events.push(StageEvent::Registration {
        residual_bits: reg.residual_rms().to_bits(),
        within_tolerance: reg_ok,
    });
    assertions_passed &= reg_ok;

    // as-built δ with a SEEDED DEFECT on the middle point.
    let design_pts = vec![design[0], design[1], design[2]];
    let mut scanned: Vec<Point2> = design_pts
        .iter()
        .map(|&point| reg.apply(point))
        .collect::<Result<_, _>>()?;
    scanned[AS_BUILT_DEFECT_INDEX] = Point2::new(
        scanned[AS_BUILT_DEFECT_INDEX].x() + AS_BUILT_DEFECT_X,
        scanned[AS_BUILT_DEFECT_INDEX].y(),
    )?;
    let diff = run_as_built_phase(
        setup_memory.budget(),
        "as-built.comparison",
        plan.difference,
        "as-built comparison refused its fixed fixture",
        |budget| {
            as_built_diff_budgeted(
                &reg,
                &design_pts,
                &scanned,
                AS_BUILT_DESIGN_TOLERANCE,
                AS_BUILT_MEASUREMENT_NOISE,
                AS_BUILT_CALIBRATION_CANDIDATE,
                cx,
                budget,
            )
            .map_err(DiffRealError::from)
        },
    )?;
    // The comparison retains the deterministic argmax; no second traversal.
    let defect_idx = Some(diff.max_deviation_index());
    let localized = defect_idx == Some(AS_BUILT_DEFECT_INDEX)
        && (diff.max_deviation() - AS_BUILT_DEFECT_X).abs() < 1e-9;
    let estimated = matches!(diff.color(), Color::Estimated { .. });
    events.push(StageEvent::AsBuiltDelta {
        max_deviation_bits: diff.max_deviation().to_bits(),
        defect_index: defect_idx,
        estimated,
    });
    assertions_passed &= localized && estimated;

    // registration-free point-sensor 4D-Var: misfit reduction.
    let prior = run_as_built_phase(
        setup_memory.budget(),
        "as-built.prior",
        plan.belief,
        "prior construction refused its fixed fixture",
        |budget| {
            Belief::diagonal_budgeted(
                AS_BUILT_PRIOR_MEAN.to_vec(),
                &AS_BUILT_PRIOR_DIAGONAL_COVARIANCE,
                cx,
                budget,
            )
            .map_err(DiffRealError::from)
        },
    )?;
    let obs = as_built_observations()?;
    let assimilated = run_as_built_phase(
        setup_memory.budget(),
        "as-built.assimilation",
        plan.assimilation,
        "colored assimilation refused its fixed fixture",
        |budget| {
            let actual = fs_assimilate::colored_assimilation_invocation_resources(
                &prior,
                &obs,
                AS_BUILT_ASSIMILATION_PARAMETER,
                AS_BUILT_ASSIMILATION_BOUNDS.0,
                AS_BUILT_ASSIMILATION_BOUNDS.1,
                cx,
            )?;
            if actual != plan.assimilation_shape {
                return Err(DiffRealError::Invocation(
                    fs_exec::InvocationError::ExplicitRefusal {
                        phase: "as-built.assimilation.preflight",
                        reason: invocation::domain_refusal(
                            "as-built.assimilation.preflight",
                            "fixed exact consumables drifted",
                        ),
                    },
                ));
            }
            assimilate_colored_budgeted(
                &prior,
                &obs,
                AS_BUILT_ASSIMILATION_PARAMETER,
                AS_BUILT_ASSIMILATION_BOUNDS.0,
                AS_BUILT_ASSIMILATION_BOUNDS.1,
                cx,
                budget,
            )
            .map_err(DiffRealError::from)
        },
    )?;
    let reduced = assimilated.misfit_after() < assimilated.misfit_before();
    events.push(StageEvent::Assimilation {
        before_bits: assimilated.misfit_before().to_bits(),
        after_bits: assimilated.misfit_after().to_bits(),
        reduced,
    });
    assertions_passed &= reduced;

    let status = if assertions_passed {
        StageStatus::Passed
    } else {
        StageStatus::Failed(StageReason::new(
            "diffreal.as-built.assertion-failed",
            "registration, defect-localization, evidence-color, or assimilation assertion failed; inspect events",
        ))
    };
    run_as_built_phase(
        setup_memory.budget(),
        "as-built.publication",
        plan.publication,
        "stage log, result identity, or retained output publication failed",
        |budget| {
            budget.charge_work(plan.publication.work())?;
            budget.charge_cost(plan.publication.cost())?;
            budget.charge_evaluations(plan.publication.evaluations())?;
            budget.poll("as-built.publication.begin")?;
            let mut publication_memory =
                budget.reserve_memory("as-built-stage-publication", plan.publication.memory())?;
            let log = StageLog::new(
                AS_BUILT_STAGE,
                StageRequirement::Required,
                status,
                AS_BUILT_EVIDENCE_IDENTITY,
                events,
            );
            let inputs = fixture_input_identity(&log);
            let result = as_built_result_identity(&reg, &diff, &assimilated);
            publication_memory
                .budget()
                .poll("as-built.publication.finalize")?;
            publication_memory
                .budget()
                .publish_output(plan.publication.output())?;
            drop(publication_memory);
            Ok((log, inputs, result))
        },
    )
}

fn execute_as_built_loop(
    cx: &Cx<'_>,
    clock: &dyn fs_exec::TimeSource,
) -> Result<StageExecution, DiffRealError> {
    let execution = DiffRealExecutionIdentity::from_cx(cx);
    let plan = invocation::AsBuiltInvocationPlan::preflight(&execution)?;
    let admission = fs_exec::InvocationAdmitter::new().admit(
        invocation::invocation_id(&execution),
        plan.limits(&execution),
        plan.total,
    )?;
    let mut root = admission.begin(cx, clock)?;
    let mut transaction = root.split_child("as-built.transaction", plan.total)?;
    let outcome = execute_as_built_science(cx, &mut transaction, &plan);
    if outcome.is_err() {
        transaction.refuse(
            AS_BUILT_STAGE,
            invocation::domain_refusal(AS_BUILT_STAGE, "scientific transaction failed"),
        );
    }
    let transaction_disposition = transaction.finish()?;
    let invocation = root.finish()?;
    let (log, inputs, result) = match outcome {
        Ok(stage)
            if transaction_disposition == fs_exec::InvocationDisposition::Completed
                && plan.verifies_receipt(&invocation, &execution) =>
        {
            stage
        }
        Ok(_) => {
            return Err(DiffRealError::InvocationDidNotComplete {
                disposition: invocation.disposition(),
                receipt_root: invocation.root(),
                receipt: Box::new(invocation),
                cause: None,
            });
        }
        Err(error) => {
            return Err(DiffRealError::InvocationDidNotComplete {
                disposition: invocation.disposition(),
                receipt_root: invocation.root(),
                receipt: Box::new(invocation),
                cause: Some(Box::new(error)),
            });
        }
    };
    Ok(StageExecution {
        log,
        inputs,
        result,
        invocation: Some(invocation),
    })
}

/// Stage 2: register a scan, estimate as-built δ, localize a defect, assimilate.
///
/// # Errors
/// Propagates fixed-work admission, cancellation, or a structured lower-layer
/// refusal and publishes no partial stage log.
pub fn stage_as_built_loop(cx: &Cx<'_>) -> Result<StageLog, DiffRealError> {
    let clock = fs_exec::WallClock::new();
    stage_as_built_loop_with_clock(cx, &clock)
}

/// Run only the as-built stage with an injected monotonic clock.
///
/// # Errors
/// As [`stage_as_built_loop`].
pub fn stage_as_built_loop_with_clock(
    cx: &Cx<'_>,
    clock: &dyn fs_exec::TimeSource,
) -> Result<StageLog, DiffRealError> {
    execute_as_built_loop(cx, clock).map(|stage| stage.log)
}

// -- Stage 3: tolerance allocation ------------------------------------------

fn tolerance_refusal(code: &'static str, detail: String) -> StageLog {
    StageLog::new(
        TOLERANCE_STAGE,
        StageRequirement::Required,
        StageStatus::Refused(StageReason::new(code, detail.clone())),
        TOLERANCE_EVIDENCE_IDENTITY,
        vec![StageEvent::Refusal { code, detail }],
    )
}

fn tolerance_derivative_failure(error: DifferentiationError) -> StageLog {
    StageLog::new(
        TOLERANCE_STAGE,
        StageRequirement::Required,
        StageStatus::Failed(StageReason::new(
            "diffreal.tolerance.sensitivity-rejected",
            "a sealed sensitivity was rejected; inspect the typed event",
        )),
        TOLERANCE_EVIDENCE_IDENTITY,
        vec![StageEvent::DifferentiationRejected { error }],
    )
}

/// Stage 3 with caller-supplied QoI samples at selected tolerance-band
/// extremes. The samples test only the local linearization; they do not
/// establish a probability or complete-corner claim.
///
/// # Errors
/// Cancellation and ambient work-budget refusals suppress the partial report.
fn execute_tolerance_allocation_with_samples(
    cx: &Cx<'_>,
    extreme_qois: &[f64],
) -> Result<StageExecution, DiffRealError> {
    admit_stage_work(cx, TOLERANCE_STAGE, TOLERANCE_WORK_UNITS)?;
    let mut events = Vec::new();
    let mut assertions_passed = true;

    let registry = production_vjp_registry()?;
    let critical = match verify_sensitivity(
        &PRODUCTION_DIFFERENTIATION_PATH,
        &registry,
        TOLERANCE_SENSITIVITY_INPUTS[0],
        cx,
    ) {
        Ok(receipt) => receipt,
        Err(error) if error.is_runtime_refusal() => {
            return Err(stage_runtime_refusal(TOLERANCE_STAGE, error));
        }
        Err(error) => {
            return Ok(StageExecution::diagnostic(tolerance_derivative_failure(
                error,
            )));
        }
    };
    let slack = match verify_sensitivity(
        &PRODUCTION_DIFFERENTIATION_PATH,
        &registry,
        TOLERANCE_SENSITIVITY_INPUTS[1],
        cx,
    ) {
        Ok(receipt) => receipt,
        Err(error) if error.is_runtime_refusal() => {
            return Err(stage_runtime_refusal(TOLERANCE_STAGE, error));
        }
        Err(error) => {
            return Ok(StageExecution::diagnostic(tolerance_derivative_failure(
                error,
            )));
        }
    };
    stage_checkpoint(cx, TOLERANCE_STAGE)?;
    events.push(sensitivity_event(&critical));
    events.push(sensitivity_event(&slack));

    if !critical.verifies_integrity() || !slack.verifies_integrity() {
        let identity = if !critical.verifies_integrity() {
            critical.identity()
        } else {
            slack.identity()
        };
        return Ok(StageExecution::diagnostic(tolerance_derivative_failure(
            DifferentiationError::SensitivityIntegrityMismatch { identity },
        )));
    }

    let feat = |name: &str, receipt: &SealedSensitivity| Feature {
        name: name.into(),
        sensitivity: receipt.gradient().abs(),
        sensitivity_color: ColorRank::Verified,
        cost_coeff: TOLERANCE_FEATURE_COST_COEFFICIENT,
        baseline_tolerance: TOLERANCE_FEATURE_BASELINE,
    };
    let budget = match variance_budget(
        TOLERANCE_PERFORMANCE_TOLERANCE,
        TOLERANCE_TARGET_PROBABILITY,
    ) {
        Ok(budget) => budget,
        Err(_error) => {
            let code = "diffreal.tolerance.invalid-budget-fixture";
            let detail =
                "the fixed tolerance-budget fixture was refused by fs-toleralloc".to_string();
            return Ok(StageExecution::diagnostic(tolerance_refusal(code, detail)));
        }
    };
    stage_checkpoint(cx, TOLERANCE_STAGE)?;
    let alloc = match allocate(
        &[
            feat(TOLERANCE_FEATURE_NAMES[0], &critical),
            feat(TOLERANCE_FEATURE_NAMES[1], &slack),
        ],
        budget,
        TOLERANCE_SIGMA_MULTIPLIER,
    ) {
        Ok(allocation) => allocation,
        Err(_error) => {
            let code = "diffreal.tolerance.allocation-refused";
            let detail =
                "the fixed tolerance-allocation fixture was refused by fs-toleralloc".to_string();
            return Ok(StageExecution::diagnostic(tolerance_refusal(code, detail)));
        }
    };
    // tighten where sensitivity is large, loosen where small.
    let critical_action = alloc
        .items
        .iter()
        .find(|item| item.name == TOLERANCE_FEATURE_NAMES[0])
        .map(|item| item.action);
    let slack_action = alloc
        .items
        .iter()
        .find(|item| item.name == TOLERANCE_FEATURE_NAMES[1])
        .map(|item| item.action);
    let tighten_high = critical_action == Some(Action::Tighten);
    let loosen_low = slack_action == Some(Action::Loosen);
    events.push(StageEvent::ToleranceActions {
        critical: critical_action,
        slack: slack_action,
    });
    assertions_passed &= tighten_high && loosen_low;

    // the GD&T report attaches a certified sensitivity to every loosened tol.
    stage_checkpoint(cx, TOLERANCE_STAGE)?;
    let report = match gdt_report(&alloc) {
        Ok(report) => report,
        Err(_error) => {
            let code = "diffreal.tolerance.report-refused";
            let detail = "the fixed GD&T report fixture was refused by fs-toleralloc".to_string();
            return Ok(StageExecution::diagnostic(tolerance_refusal(code, detail)));
        }
    };
    let loosened = report
        .iter()
        .filter(|suggestion| suggestion.action == Action::Loosen)
        .count();
    let justified = loosened > 0
        && report
            .iter()
            .filter(|suggestion| suggestion.action == Action::Loosen)
            .all(|suggestion| {
                suggestion.certified_sensitivity > 0.0 && suggestion.color == ColorRank::Verified
            });
    events.push(StageEvent::GdtJustification {
        loosened,
        all_verified: justified,
    });
    assertions_passed &= justified;

    // This checks only the supplied sample set. It is neither a complete corner
    // enumeration nor a probabilistic conformance certificate.
    stage_checkpoint(cx, TOLERANCE_STAGE)?;
    let verdict = match robustness_check(
        &alloc,
        extreme_qois,
        TOLERANCE_NOMINAL_QOI,
        TOLERANCE_SIGMA_MULTIPLIER,
        TOLERANCE_ROBUSTNESS_MARGIN,
    ) {
        Ok(verdict) => verdict,
        Err(_error) => {
            let code = "diffreal.tolerance.robustness-refused";
            let detail =
                "the fixed tolerance-robustness fixture was refused by fs-toleralloc".to_string();
            return Ok(StageExecution::diagnostic(tolerance_refusal(code, detail)));
        }
    };
    events.push(StageEvent::SampledLinearization {
        samples: extreme_qois.len(),
        confirmed: verdict.confirmed,
        linearized_std_bits: verdict.linearized_std.to_bits(),
        probability_claimed: false,
    });
    assertions_passed &= verdict.confirmed;

    let status = if assertions_passed {
        StageStatus::Passed
    } else {
        StageStatus::Failed(StageReason::new(
            "diffreal.tolerance.assertion-failed",
            "allocation direction, sensitivity justification, or sampled-extremes assertion failed; inspect events",
        ))
    };
    let log = StageLog::new(
        TOLERANCE_STAGE,
        StageRequirement::Required,
        status,
        TOLERANCE_EVIDENCE_IDENTITY,
        events,
    );
    Ok(StageExecution {
        inputs: fixture_input_identity(&log),
        result: tolerance_result_identity(&critical, &slack, &alloc, &report, &verdict),
        log,
        invocation: None,
    })
}

/// Stage 3 with caller-supplied QoI samples at selected tolerance-band
/// extremes. This is a diagnostic falsifier seam and cannot mint a receipt.
///
/// # Errors
/// Cancellation and ambient work-budget refusals suppress the partial report.
pub fn stage_tolerance_allocation_with_samples(
    cx: &Cx<'_>,
    extreme_qois: &[f64],
) -> Result<StageLog, DiffRealError> {
    execute_tolerance_allocation_with_samples(cx, extreme_qois).map(|stage| stage.log)
}

/// Stage 3: adjoint-driven GD&T using sealed local-fixture sensitivities and a
/// fixed, explicitly sampled linearization check.
///
/// # Errors
/// Cancellation and ambient work-budget refusals suppress the partial report.
pub fn stage_tolerance_allocation(cx: &Cx<'_>) -> Result<StageLog, DiffRealError> {
    execute_tolerance_allocation_with_samples(cx, &TOLERANCE_EXTREME_QOIS).map(|stage| stage.log)
}

fn execute_tolerance_allocation(cx: &Cx<'_>) -> Result<StageExecution, DiffRealError> {
    execute_tolerance_allocation_with_samples(cx, &TOLERANCE_EXTREME_QOIS)
}

// -- Stage 4: gated spacetime -----------------------------------------------

/// Stage 4: the spacetime-complex capability is not integrated and activated
/// in this battery (honestly gated, never silently passed).
///
/// # Errors
/// Cancellation or an ambient work budget below the fixed gate-recording cost
/// suppresses the partial report.
fn execute_spacetime_gated(cx: &Cx<'_>) -> Result<StageExecution, DiffRealError> {
    admit_stage_work(cx, SPACETIME_STAGE, SPACETIME_WORK_UNITS)?;
    Ok(StageExecution::diagnostic(StageLog::new(
        SPACETIME_STAGE,
        StageRequirement::Required,
        StageStatus::Gated(StageReason::new(
            "diffreal.spacetime.integration-not-activated",
            "fs-time temporal-complex support exists, but the coupled end-to-end fixture is not integrated and activated in this battery",
        )),
        SPACETIME_EVIDENCE_IDENTITY,
        vec![StageEvent::Gate {
            code: "diffreal.spacetime.integration-not-activated",
            detail: "temporal-complex dependency frankensim-epic-coupling-bk0o.7 is shipped, but this battery has no activated coupled spacetime fixture; stage not asserted"
                .to_string(),
        }],
    )))
}

/// Stage 4: record the deliberately unavailable spacetime integration gate.
///
/// # Errors
/// Cancellation or an ambient work budget below the fixed gate-recording cost
/// suppresses the partial report.
pub fn stage_spacetime_gated(cx: &Cx<'_>) -> Result<StageLog, DiffRealError> {
    execute_spacetime_gated(cx).map(|stage| stage.log)
}

#[cfg(test)]
mod report_policy_tests {
    use super::*;

    fn test_execution() -> DiffRealExecutionIdentity {
        DiffRealExecutionIdentity {
            stream_key: StreamKey {
                seed: 7,
                kernel_id: 11,
                tile: 13,
                iteration: 17,
            },
            budget: Budget::INFINITE,
            mode: ExecMode::Deterministic,
            operation_memory_limit_bytes: None,
        }
    }

    /// Run `f` under a context whose [`DiffRealExecutionIdentity`]
    /// matches [`test_execution`], so receipts sealed against the test
    /// identity verify against invocations executed here.
    fn with_test_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = fs_exec::CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 7,
                    kernel_id: 11,
                    tile: 13,
                    iteration: 17,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    fn diagnostic_passed(stage: &'static str, identity: &'static str) -> StageLog {
        StageLog::new(
            stage,
            StageRequirement::Required,
            StageStatus::Passed,
            identity,
            vec![StageEvent::Gate {
                code: "test.fixture.executed",
                detail: format!("{stage} fixture executed"),
            }],
        )
    }

    /// Seal a policy-fixture report. The as-built stage runs its REAL
    /// composed invocation (receipts for that stage refuse `None`
    /// invocations since the affine-invocation hardening), while the
    /// other stages stay cheap diagnostics.
    fn sealed_report(spacetime_status: StageStatus, extra: Option<StageLog>) -> DiffRealReport {
        with_test_cx(|cx| {
            let clock = fs_exec::VirtualClock::new();
            let as_built = execute_as_built_loop(cx, &clock)
                .expect("the fixed as-built loop completes under an infinite budget");
            let mut stages = vec![
                StageExecution::diagnostic(diagnostic_passed(
                    DIFFERENTIATION_STAGE,
                    DIFFERENTIATION_EVIDENCE_IDENTITY,
                )),
                as_built,
                StageExecution::diagnostic(diagnostic_passed(
                    TOLERANCE_STAGE,
                    TOLERANCE_EVIDENCE_IDENTITY,
                )),
                StageExecution::diagnostic(StageLog::new(
                    SPACETIME_STAGE,
                    StageRequirement::Required,
                    spacetime_status,
                    SPACETIME_EVIDENCE_IDENTITY,
                    vec![StageEvent::Gate {
                        code: "test.spacetime.disposition",
                        detail: "spacetime fixture disposition recorded".to_string(),
                    }],
                )),
            ];
            if let Some(log) = extra {
                stages.push(StageExecution::diagnostic(log));
            }
            DiffRealReport::seal(stages, DiffRealExecutionIdentity::from_cx(cx))
        })
    }

    fn required_status_report(spacetime_status: StageStatus) -> DiffRealReport {
        sealed_report(spacetime_status, None)
    }

    #[test]
    fn all_passed_required_stages_are_structurally_ready_but_not_authenticated() {
        let report = required_status_report(StageStatus::Passed);
        assert!(report.complete());
        assert!(report.all_required_passed());
        assert!(report.structurally_ready());
    }

    #[test]
    fn a_failed_required_stage_is_complete_but_not_promotion_ready() {
        let report = required_status_report(StageStatus::Failed(StageReason::new(
            "test.spacetime.assertion-failed",
            "the spacetime fixture ran and violated its asserted bound",
        )));
        assert!(
            report.complete(),
            "failed is an evaluated scientific outcome"
        );
        assert!(!report.all_required_passed());
        assert!(!report.structurally_ready());
    }

    #[test]
    fn a_gated_required_stage_is_neither_complete_nor_promotion_ready() {
        let report = required_status_report(StageStatus::Gated(StageReason::new(
            "test.spacetime.gated",
            "the required capability is unavailable",
        )));
        assert!(!report.complete());
        assert!(!report.all_required_passed());
        assert!(!report.structurally_ready());
    }

    #[test]
    fn a_refused_required_stage_is_neither_complete_nor_promotion_ready() {
        let report = required_status_report(StageStatus::Refused(StageReason::new(
            "test.spacetime.refused",
            "the stage exhausted its admitted budget before evaluation",
        )));
        assert!(!report.complete());
        assert!(!report.all_required_passed());
        assert!(!report.structurally_ready());
    }

    #[test]
    fn an_explicit_optional_gate_does_not_block_required_stage_promotion() {
        let report = sealed_report(
            StageStatus::Passed,
            Some(StageLog::new(
                "diagnostic-only",
                StageRequirement::Optional,
                StageStatus::Gated(StageReason::new(
                    "test.optional.gated",
                    "the optional diagnostic backend is unavailable",
                )),
                "fs-diffreal-e2e/optional-diagnostic/v1",
                vec![StageEvent::Gate {
                    code: "test.optional.gated",
                    detail: "optional diagnostic gate retained".to_string(),
                }],
            )),
        );
        assert!(report.complete());
        assert!(report.all_required_passed());
        assert!(report.structurally_ready());
    }

    #[test]
    fn receipt_mutation_reorder_omission_and_unknown_versions_fail_closed() {
        let report = required_status_report(StageStatus::Passed);
        assert!(report.verifies_integrity());

        let mut result_mutation = report.clone();
        result_mutation.receipts[0].result.0[0] ^= 1;
        assert!(!result_mutation.verifies_integrity());

        let mut input_mutation = report.clone();
        input_mutation.receipts[0].fixture_inputs.0[0] ^= 1;
        assert!(!input_mutation.verifies_integrity());

        let mut stage_root_mutation = report.clone();
        stage_root_mutation.receipts[0].root.0[0] ^= 1;
        assert!(!stage_root_mutation.verifies_integrity());

        let mut report_root_mutation = report.clone();
        report_root_mutation.receipt_root.0[0] ^= 1;
        assert!(!report_root_mutation.verifies_integrity());

        let mut reordered = report.clone();
        reordered.receipts.swap(0, 1);
        assert!(!reordered.verifies_integrity());

        let mut omitted = report.clone();
        omitted.receipts.pop();
        assert!(!omitted.verifies_integrity());

        let mut unknown_stage = report.clone();
        unknown_stage.receipts[0].version = STAGE_RECEIPT_POLICY_VERSION + 1;
        assert!(!unknown_stage.verifies_integrity());

        let mut unknown_report = report;
        unknown_report.receipt_version = REPORT_RECEIPT_POLICY_VERSION + 1;
        assert!(!unknown_report.verifies_integrity());
    }

    #[test]
    fn canonical_stage_and_execution_fields_are_receipt_semantic() {
        let report = required_status_report(StageStatus::Passed);

        let mut status = report.clone();
        status.stages[0].status = StageStatus::Failed(StageReason::new(
            "test.mutated",
            "status changed after sealing",
        ));
        assert!(!status.verifies_integrity());

        let mut reason = required_status_report(StageStatus::Failed(StageReason::new(
            "test.failed",
            "original reason",
        )));
        let StageStatus::Failed(reason_payload) = &mut reason.stages[3].status else {
            panic!("synthetic failure has a reason")
        };
        reason_payload.detail.push_str(" mutated");
        assert!(!reason.verifies_integrity());

        let mut evidence_identity = report.clone();
        evidence_identity.stages[0].evidence_identity = "mutated/v1";
        assert!(!evidence_identity.verifies_integrity());

        let mut event_payload = report.clone();
        let StageEvent::Gate { detail, .. } = &mut event_payload.stages[0].events[0] else {
            panic!("synthetic stage uses one gate-shaped diagnostic")
        };
        detail.push_str(" mutated");
        assert!(!event_payload.verifies_integrity());

        let mut event_order = report.clone();
        event_order.stages[0].events.push(StageEvent::Gate {
            code: "test.second",
            detail: "second event".to_string(),
        });
        event_order.stages[0].events.swap(0, 1);
        assert!(!event_order.verifies_integrity());

        let mut execution = report.clone();
        execution.execution.stream_key.seed ^= 1;
        assert!(!execution.verifies_integrity());

        let mut budget = report;
        budget.execution.budget = budget.execution.budget.with_priority(9);
        assert!(!budget.verifies_integrity());
    }

    #[test]
    fn authentication_rejects_unknown_versions_before_authority_dispatch() {
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        let gate = fs_exec::CancelGate::new();
        pool.scope(|arena| {
            let execution = test_execution();
            let cx = Cx::new(
                &gate,
                arena,
                execution.stream_key,
                execution.budget,
                execution.mode,
            );
            let report = run_battery(&cx).expect("fixed battery produces a report");
            let attestation =
                PromotionAttestation::new("unused", Vec::new(), promotion_policy_fingerprint());

            let mut unknown_report = report.clone();
            unknown_report.receipt_version = REPORT_RECEIPT_POLICY_VERSION + 1;
            assert!(matches!(
                unknown_report.authenticate(&cx, attestation.clone(), &NoPromotionReceiptVerifier,),
                Err(PromotionReceiptError::UnknownReportVersion { .. })
            ));

            let mut unknown_stage = report;
            unknown_stage.receipts[0].version = STAGE_RECEIPT_POLICY_VERSION + 1;
            assert!(matches!(
                unknown_stage.authenticate(&cx, attestation, &NoPromotionReceiptVerifier),
                Err(PromotionReceiptError::UnknownStageVersion { .. })
            ));
        });
    }

    #[derive(Debug)]
    struct MustNotDispatch;

    impl PromotionReceiptVerifier for MustNotDispatch {
        fn verify(
            &self,
            _request: &PromotionVerificationRequest<'_>,
            _attestation: &PromotionAttestation,
        ) -> PromotionVerificationDecision {
            panic!("semantic validation must precede authority dispatch")
        }
    }

    #[test]
    fn contradictory_but_rehashed_stage_fails_before_authority_dispatch() {
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        let gate = fs_exec::CancelGate::new();
        pool.scope(|arena| {
            let execution = test_execution();
            let cx = Cx::new(
                &gate,
                arena,
                execution.stream_key,
                execution.budget,
                execution.mode,
            );
            let mut report = run_battery(&cx).expect("fixed battery produces a report");
            let StageEvent::MissingVjpProbe { blocked, .. } = &mut report.stages[0].events[1]
            else {
                panic!("differentiation transcript retains the kill probe")
            };
            *blocked = false;
            let receipt = &mut report.receipts[0];
            receipt.root = stage_receipt_identity(
                receipt.version,
                &report.stages[0],
                &report.execution,
                receipt.fixture_inputs,
                receipt.result,
                receipt.invocation.as_ref(),
            );
            report.receipt_root = report_receipt_identity(
                report.receipt_version,
                &report.execution,
                &report.receipts,
            );
            assert!(
                report.verifies_integrity(),
                "the forged root is self-consistent"
            );
            let attestation =
                PromotionAttestation::new("unused", Vec::new(), promotion_policy_fingerprint());
            assert!(matches!(
                report.authenticate(&cx, attestation, &MustNotDispatch),
                Err(PromotionReceiptError::InvalidStageSemantics)
            ));
        });
    }

    #[test]
    fn sealed_sensitivity_integrity_rejects_field_tampering() {
        let ops: Vec<String> = PRODUCTION_DIFFERENTIATION_PATH
            .iter()
            .map(|op| (*op).to_string())
            .collect();
        let fields = [
            1.5_f64.to_bits(),
            16.0_f64.to_bits(),
            16.0_f64.to_bits(),
            16.0_f64.to_bits(),
            16.0_f64.to_bits(),
            16.0_f64.to_bits(),
            1.0e-12_f64.to_bits(),
        ];
        let identity = sensitivity_identity(&ops, fields);
        let mut receipt = SealedSensitivity {
            ops,
            input_bits: fields[0],
            value_bits: fields[1],
            production_gradient_bits: fields[2],
            dual_gradient_bits: fields[3],
            fd_coarse_bits: fields[4],
            fd_fine_bits: fields[5],
            fd_tolerance_bits: fields[6],
            identity,
        };
        assert!(receipt.verifies_integrity());
        receipt.production_gradient_bits ^= 1;
        assert!(!receipt.verifies_integrity());
    }

    #[test]
    fn malformed_or_schema_incomplete_reports_fail_closed() {
        let all_passed = required_status_report(StageStatus::Passed);

        let mut missing = all_passed.clone();
        missing
            .stages
            .retain(|stage| stage.stage != SPACETIME_STAGE);
        assert!(!missing.complete());
        assert!(!missing.all_required_passed());

        let mut duplicate = all_passed.clone();
        duplicate.stages.push(duplicate.stages[0].clone());
        assert!(!duplicate.complete());
        assert!(!duplicate.all_required_passed());

        let mut mismatched_identity = all_passed.clone();
        mismatched_identity.stages[0].evidence_identity = "wrong-fixture/v1";
        assert!(!mismatched_identity.complete());
        assert!(!mismatched_identity.all_required_passed());

        let mut reordered = all_passed.clone();
        reordered.stages.swap(0, 1);
        assert!(!reordered.complete());
        assert!(!reordered.all_required_passed());

        let mut blank_reason = all_passed.clone();
        blank_reason.stages[3].status = StageStatus::Failed(StageReason::new("", ""));
        assert!(!blank_reason.complete());
        assert!(!blank_reason.all_required_passed());

        let mut empty_log = all_passed;
        empty_log.stages[0].events.clear();
        assert!(!empty_log.complete());
        assert!(!empty_log.all_required_passed());
    }
}
