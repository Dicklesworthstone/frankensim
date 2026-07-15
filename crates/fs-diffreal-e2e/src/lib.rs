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
use fs_asbuilt::{Fiducial, Point2, as_built_diff, register};
use fs_assimilate::{AssimError, Belief, assimilate_colored, misfit, point_sensor};
use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::Color;
use fs_exec::Cx;
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
pub const AS_BUILT_EVIDENCE_IDENTITY: &str = "fs-diffreal-e2e/as-built-fixture/v1";
/// Versioned fixture identity expected for the tolerance-allocation stage.
pub const TOLERANCE_EVIDENCE_IDENTITY: &str = "fs-diffreal-e2e/tolerance-allocation-fixture/v3";
/// Versioned fixture identity expected for the spacetime-integration stage.
pub const SPACETIME_EVIDENCE_IDENTITY: &str = "fs-diffreal-e2e/spacetime-integration-gate/v1";

/// Version of the production-path sensitivity sealing policy.
pub const SENSITIVITY_POLICY_VERSION: &str = "fs-diffreal-e2e/sensitivity-policy/v1";

/// The production differentiation fixture's operator path.
pub const PRODUCTION_DIFFERENTIATION_PATH: [&str; 3] = ["sdf", "spline", "solve"];

const SENSITIVITY_IDENTITY_DOMAIN: &str = "frankensim.fs-diffreal-e2e.sensitivity.v1";
const MAX_DIFFERENTIATION_OPS: usize = 16;
const MAX_OP_NAME_BYTES: usize = 64;
const DIFFERENTIATION_WORK_UNITS: u64 = 12;
const DIFFERENTIATION_STAGE_WORK_UNITS: u64 = 24;
const AS_BUILT_WORK_UNITS: u64 = 64;
const TOLERANCE_WORK_UNITS: u64 = 32;
const SPACETIME_WORK_UNITS: u64 = 1;

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
}

impl DiffRealReport {
    /// Ordered read-only stage diagnostics produced by this battery run.
    #[must_use]
    pub fn stages(&self) -> &[StageLog] {
        &self.stages
    }

    /// Is every required stage present exactly once with the expected evidence
    /// identity and an evaluated (`Passed` or `Failed`) result?
    #[must_use]
    pub fn complete(&self) -> bool {
        self.required_schema_is_valid()
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
        self.required_schema_is_valid()
            && REQUIRED_STAGES
                .iter()
                .all(|required| self.stage(required.name).is_some_and(StageLog::passed))
    }

    /// Did this crate-authored fixed battery complete every required fixture
    /// and pass every required assertion?
    ///
    /// This is battery-local readiness only. It is not scientific release
    /// admission, external validation, or an authenticated promotion receipt.
    #[must_use]
    pub fn promotion_ready(&self) -> bool {
        self.complete() && self.all_required_passed()
    }

    /// Fail-closed compatibility alias for [`Self::promotion_ready`].
    #[must_use]
    pub fn passed(&self) -> bool {
        self.promotion_ready()
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

/// Typed refusal that prevents publication of a partial battery report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DiffRealError {
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
            Self::Differentiation(error) => Some(error),
            Self::AsBuilt(error) => Some(error),
            Self::Assimilation(error) => Some(error),
            Self::Cancelled { .. } | Self::WorkBudgetExceeded { .. } => None,
        }
    }
}

impl From<fs_asbuilt::RegError> for DiffRealError {
    fn from(error: fs_asbuilt::RegError) -> Self {
        Self::AsBuilt(error)
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
    Ok(DiffRealReport {
        stages: vec![
            stage_differentiation(cx)?,
            stage_as_built_loop(cx)?,
            stage_tolerance_allocation(cx)?,
            stage_spacetime_gated(cx)?,
        ],
    })
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
    let h = 1.0e-4 * x.abs().max(1.0);
    let fd = fd_falsifier(
        &|point| composite(point[0]),
        &[x],
        &[1.0],
        production.gradient(),
        h,
        64.0 * f64::EPSILON,
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
pub fn stage_differentiation_with_registry(
    cx: &Cx<'_>,
    registry: &DifferentiationRegistry,
) -> Result<StageLog, DiffRealError> {
    admit_stage_work(cx, DIFFERENTIATION_STAGE, DIFFERENTIATION_STAGE_WORK_UNITS)?;
    let sensitivity = match verify_sensitivity(&PRODUCTION_DIFFERENTIATION_PATH, registry, 1.5, cx)
    {
        Ok(sensitivity) => sensitivity,
        Err(error) if error.is_runtime_refusal() => {
            return Err(stage_runtime_refusal(DIFFERENTIATION_STAGE, error));
        }
        Err(error) => {
            return Ok(StageLog::new(
                DIFFERENTIATION_STAGE,
                StageRequirement::Required,
                StageStatus::Failed(StageReason::new(
                    "diffreal.differentiation.production-rejected",
                    error.to_string(),
                )),
                DIFFERENTIATION_EVIDENCE_IDENTITY,
                vec![StageEvent::DifferentiationRejected { error }],
            ));
        }
    };

    let mut events = vec![sensitivity_event(&sensitivity)];
    let remesh = ["sdf", "remesh", "solve"];
    let missing_probe = differentiate_path(&remesh, registry, 1.5, cx);
    let blocked = matches!(
        &missing_probe,
        Err(DifferentiationError::MissingVjp { op }) if op == "remesh"
    );
    events.push(StageEvent::MissingVjpProbe {
        op: "remesh".to_string(),
        blocked,
    });
    if let Err(error) = missing_probe {
        if error.is_runtime_refusal() {
            return Err(stage_runtime_refusal(DIFFERENTIATION_STAGE, error));
        }
        if !matches!(error, DifferentiationError::MissingVjp { .. }) {
            events.push(StageEvent::DifferentiationRejected { error });
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
    Ok(StageLog::new(
        DIFFERENTIATION_STAGE,
        StageRequirement::Required,
        status,
        DIFFERENTIATION_EVIDENCE_IDENTITY,
        events,
    ))
}

/// Stage 1: production tape/VJP gradient, independent dual/FD sealing, and a
/// missing-VJP kill probe.
///
/// # Errors
/// Cancellation and ambient work-budget refusals suppress the partial report.
pub fn stage_differentiation(cx: &Cx<'_>) -> Result<StageLog, DiffRealError> {
    let registry = production_vjp_registry()?;
    stage_differentiation_with_registry(cx, &registry)
}

// -- Stage 2: as-built loop -------------------------------------------------

/// Stage 2: register a scan, estimate as-built δ, localize a defect, assimilate.
///
/// # Errors
/// Propagates fixed-work admission, cancellation, or a structured lower-layer
/// refusal and publishes no partial stage log.
pub fn stage_as_built_loop(cx: &Cx<'_>) -> Result<StageLog, DiffRealError> {
    admit_stage_work(cx, AS_BUILT_STAGE, AS_BUILT_WORK_UNITS)?;
    let mut events = Vec::new();
    let mut assertions_passed = true;

    // a scanned fixture: design datums transformed by a known rigid motion.
    let design = [
        Point2::new(0.0, 0.0)?,
        Point2::new(2.0, 0.0)?,
        Point2::new(0.0, 2.0)?,
    ];
    let (theta, tx, ty) = (0.3_f64, 4.0, 1.0);
    let xf = |p: Point2| {
        let (s, c) = theta.sin_cos();
        Point2::new(c * p.x() - s * p.y() + tx, s * p.x() + c * p.y() + ty)
    };
    let fids: Vec<Fiducial> = design
        .iter()
        .map(|&datum| Ok(Fiducial::new(datum, xf(datum)?)))
        .collect::<Result<_, fs_asbuilt::RegError>>()?;
    let reg = register(&fids, cx)?;
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
    scanned[1] = Point2::new(scanned[1].x() + 0.3, scanned[1].y())?;
    let diff = as_built_diff(&reg, &design_pts, &scanned, 0.5, 0.02, "cmm-cal-2026", cx)?;
    // localize the defect: the argmax deviation is the seeded point (index 1).
    let defect_idx = diff
        .deviations()
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(b.1))
        .map(|(i, _)| i);
    let localized = defect_idx == Some(1) && (diff.max_deviation() - 0.3).abs() < 1e-9;
    let estimated = matches!(diff.color(), Color::Estimated { .. });
    events.push(StageEvent::AsBuiltDelta {
        max_deviation_bits: diff.max_deviation().to_bits(),
        defect_index: defect_idx,
        estimated,
    });
    assertions_passed &= localized && estimated;

    // registration-free point-sensor 4D-Var: misfit reduction.
    let prior = Belief::diagonal(vec![20.0, 20.0], &[9.0, 9.0], cx)?;
    let obs = vec![
        point_sensor(0, 2, 24.0, 0.25, "thermocouple-1")?,
        point_sensor(1, 2, 18.5, 0.25, "thermocouple-2")?,
    ];
    let assimilated = assimilate_colored(&prior, &obs, "Re", 1e5, 3e5, cx)?;
    let misfit_reduced = assimilated.misfit_after() < assimilated.misfit_before();
    let checked_after = misfit(assimilated.belief(), &obs, cx)?;
    let checked_before = misfit(&prior, &obs, cx)?;
    let reduced = misfit_reduced && checked_after <= checked_before;
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
    Ok(StageLog::new(
        AS_BUILT_STAGE,
        StageRequirement::Required,
        status,
        AS_BUILT_EVIDENCE_IDENTITY,
        events,
    ))
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
            error.to_string(),
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
pub fn stage_tolerance_allocation_with_samples(
    cx: &Cx<'_>,
    extreme_qois: &[f64],
) -> Result<StageLog, DiffRealError> {
    admit_stage_work(cx, TOLERANCE_STAGE, TOLERANCE_WORK_UNITS)?;
    let mut events = Vec::new();
    let mut assertions_passed = true;

    let registry = production_vjp_registry()?;
    let critical = match verify_sensitivity(&PRODUCTION_DIFFERENTIATION_PATH, &registry, 1.0, cx) {
        Ok(receipt) => receipt,
        Err(error) if error.is_runtime_refusal() => {
            return Err(stage_runtime_refusal(TOLERANCE_STAGE, error));
        }
        Err(error) => return Ok(tolerance_derivative_failure(error)),
    };
    let slack = match verify_sensitivity(&PRODUCTION_DIFFERENTIATION_PATH, &registry, -0.475, cx) {
        Ok(receipt) => receipt,
        Err(error) if error.is_runtime_refusal() => {
            return Err(stage_runtime_refusal(TOLERANCE_STAGE, error));
        }
        Err(error) => return Ok(tolerance_derivative_failure(error)),
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
        return Ok(tolerance_derivative_failure(
            DifferentiationError::SensitivityIntegrityMismatch { identity },
        ));
    }

    let feat = |name: &str, receipt: &SealedSensitivity| Feature {
        name: name.into(),
        sensitivity: receipt.gradient().abs(),
        sensitivity_color: ColorRank::Verified,
        cost_coeff: 1.0,
        baseline_tolerance: 0.5,
    };
    let budget = match variance_budget(1.0, 0.99) {
        Ok(budget) => budget,
        Err(error) => {
            let code = "diffreal.tolerance.invalid-budget-fixture";
            let detail = format!("the fixed tolerance-budget fixture was refused: {error:?}");
            return Ok(tolerance_refusal(code, detail));
        }
    };
    stage_checkpoint(cx, TOLERANCE_STAGE)?;
    let alloc = match allocate(
        &[feat("critical", &critical), feat("slack", &slack)],
        budget,
        3.0,
    ) {
        Ok(allocation) => allocation,
        Err(error) => {
            let code = "diffreal.tolerance.allocation-refused";
            let detail = format!("the fixed tolerance-allocation fixture was refused: {error:?}");
            return Ok(tolerance_refusal(code, detail));
        }
    };
    // tighten where sensitivity is large, loosen where small.
    let critical_action = alloc
        .items
        .iter()
        .find(|item| item.name == "critical")
        .map(|item| item.action);
    let slack_action = alloc
        .items
        .iter()
        .find(|item| item.name == "slack")
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
        Err(error) => {
            let code = "diffreal.tolerance.report-refused";
            let detail = format!("the fixed GD&T report fixture was refused: {error:?}");
            return Ok(tolerance_refusal(code, detail));
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
    let verdict = match robustness_check(&alloc, extreme_qois, 0.0, 3.0, 0.2) {
        Ok(verdict) => verdict,
        Err(error) => {
            let code = "diffreal.tolerance.robustness-refused";
            let detail = format!("the fixed tolerance-robustness fixture was refused: {error:?}");
            return Ok(tolerance_refusal(code, detail));
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
    Ok(StageLog::new(
        TOLERANCE_STAGE,
        StageRequirement::Required,
        status,
        TOLERANCE_EVIDENCE_IDENTITY,
        events,
    ))
}

/// Stage 3: adjoint-driven GD&T using sealed local-fixture sensitivities and a
/// fixed, explicitly sampled linearization check.
///
/// # Errors
/// Cancellation and ambient work-budget refusals suppress the partial report.
pub fn stage_tolerance_allocation(cx: &Cx<'_>) -> Result<StageLog, DiffRealError> {
    stage_tolerance_allocation_with_samples(cx, &[0.9, -0.8, 0.5])
}

// -- Stage 4: gated spacetime -----------------------------------------------

/// Stage 4: the spacetime-complex capability is not integrated and activated
/// in this battery (honestly gated, never silently passed).
///
/// # Errors
/// Cancellation or an ambient work budget below the fixed gate-recording cost
/// suppresses the partial report.
pub fn stage_spacetime_gated(cx: &Cx<'_>) -> Result<StageLog, DiffRealError> {
    admit_stage_work(cx, SPACETIME_STAGE, SPACETIME_WORK_UNITS)?;
    Ok(StageLog::new(
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
    ))
}

#[cfg(test)]
mod report_policy_tests {
    use super::*;

    fn required_status_report(spacetime_status: StageStatus) -> DiffRealReport {
        let passed = |stage, identity| {
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
        };
        DiffRealReport {
            stages: vec![
                passed(DIFFERENTIATION_STAGE, DIFFERENTIATION_EVIDENCE_IDENTITY),
                passed(AS_BUILT_STAGE, AS_BUILT_EVIDENCE_IDENTITY),
                passed(TOLERANCE_STAGE, TOLERANCE_EVIDENCE_IDENTITY),
                StageLog::new(
                    SPACETIME_STAGE,
                    StageRequirement::Required,
                    spacetime_status,
                    SPACETIME_EVIDENCE_IDENTITY,
                    vec![StageEvent::Gate {
                        code: "test.spacetime.disposition",
                        detail: "spacetime fixture disposition recorded".to_string(),
                    }],
                ),
            ],
        }
    }

    #[test]
    fn all_passed_required_stages_are_complete_and_promotion_ready() {
        let report = required_status_report(StageStatus::Passed);
        assert!(report.complete());
        assert!(report.all_required_passed());
        assert!(report.promotion_ready());
        assert!(report.passed());
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
        assert!(!report.promotion_ready());
        assert!(!report.passed());
    }

    #[test]
    fn a_gated_required_stage_is_neither_complete_nor_promotion_ready() {
        let report = required_status_report(StageStatus::Gated(StageReason::new(
            "test.spacetime.gated",
            "the required capability is unavailable",
        )));
        assert!(!report.complete());
        assert!(!report.all_required_passed());
        assert!(!report.promotion_ready());
        assert!(!report.passed());
    }

    #[test]
    fn a_refused_required_stage_is_neither_complete_nor_promotion_ready() {
        let report = required_status_report(StageStatus::Refused(StageReason::new(
            "test.spacetime.refused",
            "the stage exhausted its admitted budget before evaluation",
        )));
        assert!(!report.complete());
        assert!(!report.all_required_passed());
        assert!(!report.promotion_ready());
        assert!(!report.passed());
    }

    #[test]
    fn an_explicit_optional_gate_does_not_block_required_stage_promotion() {
        let mut report = required_status_report(StageStatus::Passed);
        report.stages.push(StageLog::new(
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
        ));
        assert!(report.complete());
        assert!(report.all_required_passed());
        assert!(report.promotion_ready());
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
