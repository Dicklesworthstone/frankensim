//! fs-viz — scientific visualization primitives. Layer: L5.
//!
//! A 10⁸-cell field is illegible as raw numbers; the job of visualization is to
//! compress it into a picture an agent can read in one glance. This v0 is the
//! verifiable core of that job — the pieces with ANALYTIC ground truth, so a
//! rendered topology is a certificate rather than a guess:
//!
//! - [`streamline_with_cx`] — admitted, cancellable, atomically published RK4
//!   integration of a flow field; on a rotation field the streamline is a
//!   circle (radius conserved), on a saddle field it is a hyperbola (`xy`
//!   conserved). [`streamline`] is the no-authority compatibility wrapper;
//! - [`classify_hessian`] — the Morse critical-point type (min / max / saddle)
//!   and index from a scalar field's Hessian — the atom of a Morse–Smale
//!   complex;
//! - [`Grid2::isocontour_crossings`] — bounded, fail-closed isocontour edge
//!   intersections of an admitted scalar grid, which on a circle SDF all lie
//!   on the circle.
//! - [`Grid3::isosurface`] — bounded deterministic marching tetrahedra from an
//!   owned x-fastest scalar grid to a renderer-ready indexed triangle mesh.
//! - [`ScalarField3`] — a bounded versioned artifact codec with explicit
//!   node/cell centering, quantity, and units for ledger composition.
//!
//! Deterministic for deterministic callbacks. Scoped streamline and contour
//! extraction consume the lower-layer [`fs_exec::Cx`] contract and identify
//! published artifacts with the workspace's [`fs_blake3`] owner.

mod isosurface;
mod scalar_field;

pub use isosurface::{Grid3, Grid3Error, IsoMesh3, IsoSurfaceError, Vec3};
pub use scalar_field::{
    SCALAR_FIELD3_ARTIFACT_KIND, SCALAR_FIELD3_SCHEMA_VERSION, ScalarField3, ScalarField3Error,
    ScalarFieldSemantics, ScalarLayout3,
};

use core::mem::size_of;
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::{AdmittedBudget, BudgetRefusal, Cx};

/// Domain for the canonical identity of a published Grid2 crossing artifact.
pub const ISO_CONTOUR_ARTIFACT_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-viz.isocontour-crossings.v1";

/// Version of the canonical Grid2 crossing-artifact preimage.
pub const ISO_CONTOUR_ARTIFACT_IDENTITY_VERSION: u32 = 1;

/// A 2-D point / vector.
pub type Vec2 = [f64; 2];

/// Domain for canonical identities of published fixed-step streamline artifacts.
pub const STREAMLINE_ARTIFACT_IDENTITY_DOMAIN: &str = "org.frankensim.fs-viz.streamline-rk4.v1";

/// Version of the canonical streamline request/result preimage.
pub const STREAMLINE_ARTIFACT_IDENTITY_VERSION: u32 = 1;

/// Fixed integration method admitted by the v1 streamline service.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StreamlineMethod {
    /// Classical explicit four-stage, fourth-order Runge-Kutta.
    Rk4,
}

impl StreamlineMethod {
    const fn identity_tag(self) -> u8 {
        match self {
            Self::Rk4 => 1,
        }
    }
}

/// Explicit unit interpretation for a streamline request.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StreamlineUnits {
    /// Coordinates, vector values, and integration time are dimensionless.
    Dimensionless,
}

impl StreamlineUnits {
    const fn identity_tag(self) -> u8 {
        match self {
            Self::Dimensionless => 1,
        }
    }
}

/// Axis-aligned closed domain admitted for an optional boundary policy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StreamlineDomain2 {
    /// Inclusive lower coordinate bounds.
    pub lower: Vec2,
    /// Inclusive upper coordinate bounds.
    pub upper: Vec2,
}

/// Action when a finite RK4 candidate leaves the admitted closed domain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StreamlineBoundaryPolicy {
    /// Refuse the entire operation; no partial trajectory publishes.
    RefuseExit,
    /// Publish the prefix ending at the last point inside the domain.
    StopBeforeExit,
}

impl StreamlineBoundaryPolicy {
    const fn identity_tag(self) -> u8 {
        match self {
            Self::RefuseExit => 1,
            Self::StopBeforeExit => 2,
        }
    }
}

/// Action when one RK4 step rounds back to the current point bit-for-bit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StreamlineStagnationPolicy {
    /// Retain the repeated point and continue through the requested step count.
    RetainRepeatedPoints,
    /// Publish the unique prefix without appending the repeated candidate.
    StopBeforeRepeat,
}

impl StreamlineStagnationPolicy {
    const fn identity_tag(self) -> u8 {
        match self {
            Self::RetainRepeatedPoints => 1,
            Self::StopBeforeRepeat => 2,
        }
    }
}

/// Fully explicit semantic and polling request for one fixed-step streamline.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StreamlineRequest {
    /// Initial point, always retained on successful zero-step work.
    pub seed: Vec2,
    /// Finite, nonzero signed step. Negative values request reverse time.
    pub dt: f64,
    /// Number of complete RK4 steps requested.
    pub steps: usize,
    /// Versioned numerical method.
    pub method: StreamlineMethod,
    /// Requested method/receipt schema version.
    pub method_version: u32,
    /// Explicit unit interpretation.
    pub units: StreamlineUnits,
    /// Optional closed integration domain.
    pub domain: Option<StreamlineDomain2>,
    /// Policy for a candidate outside `domain`.
    pub boundary_policy: StreamlineBoundaryPolicy,
    /// Policy for a bit-identical repeated point.
    pub stagnation_policy: StreamlineStagnationPolicy,
    /// Maximum complete RK4 steps or identity points between checkpoints.
    pub items_per_poll: usize,
}

impl StreamlineRequest {
    /// Construct the dimensionless RK4 request used by the compatibility path.
    #[must_use]
    pub const fn dimensionless_rk4(
        seed: Vec2,
        dt: f64,
        steps: usize,
        items_per_poll: usize,
    ) -> Self {
        Self {
            seed,
            dt,
            steps,
            method: StreamlineMethod::Rk4,
            method_version: STREAMLINE_ARTIFACT_IDENTITY_VERSION,
            units: StreamlineUnits::Dimensionless,
            domain: None,
            boundary_policy: StreamlineBoundaryPolicy::RefuseExit,
            stagnation_policy: StreamlineStagnationPolicy::RetainRepeatedPoints,
            items_per_poll,
        }
    }
}

/// Separately admitted resource in one streamline transaction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StreamlineResource {
    /// Requested RK4 steps.
    Steps,
    /// Vector-field callback attempts.
    FieldEvaluations,
    /// Seed plus staged step results.
    OutputPoints,
    /// Logical point-payload bytes.
    OutputBytes,
    /// Fixed RK4 and identity scratch bytes.
    ScratchBytes,
    /// Terminal diagnostic records.
    DiagnosticRecords,
    /// Fixed report/error bytes.
    DiagnosticBytes,
    /// Simultaneously live request, output, scratch, and diagnostic bytes.
    LiveBytes,
    /// Canonical artifact-identity bytes hashed.
    IdentityBytes,
    /// Bounded cancellation checkpoints.
    Polls,
    /// Deterministic scalar work charged to the ambient execution budget.
    WorkUnits,
}

impl core::fmt::Display for StreamlineResource {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        formatter.write_str(match self {
            Self::Steps => "steps",
            Self::FieldEvaluations => "field evaluations",
            Self::OutputPoints => "output points",
            Self::OutputBytes => "output bytes",
            Self::ScratchBytes => "scratch bytes",
            Self::DiagnosticRecords => "diagnostic records",
            Self::DiagnosticBytes => "diagnostic bytes",
            Self::LiveBytes => "simultaneously live bytes",
            Self::IdentityBytes => "identity bytes",
            Self::Polls => "polls",
            Self::WorkUnits => "work units",
        })
    }
}

/// Complete caller envelope for one streamline transaction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StreamlineBudget {
    /// Maximum admitted RK4 steps.
    pub step_limit: usize,
    /// Maximum vector-field callback attempts.
    pub field_evaluation_limit: usize,
    /// Maximum seed-plus-step output count.
    pub output_point_limit: usize,
    /// Maximum logical point-payload bytes.
    pub output_byte_limit: usize,
    /// Maximum fixed temporary bytes.
    pub scratch_byte_limit: usize,
    /// Maximum terminal diagnostic records.
    pub diagnostic_record_limit: usize,
    /// Maximum fixed report/error bytes.
    pub diagnostic_byte_limit: usize,
    /// Maximum simultaneously live accounted bytes.
    pub live_byte_limit: usize,
    /// Maximum canonical identity bytes hashed.
    pub identity_byte_limit: usize,
    /// Maximum checkpoint attempts.
    pub poll_limit: usize,
    /// Maximum deterministic work units.
    pub work_unit_limit: u64,
}

/// Checked conservative plan for a complete streamline request.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StreamlinePlan {
    /// Requested RK4 steps.
    pub steps: usize,
    /// Four callback attempts per complete RK4 step.
    pub field_evaluations: usize,
    /// Seed plus one possible result per step.
    pub output_points: usize,
    /// Bytes retained by the request and operation envelope.
    pub input_bytes: usize,
    /// Logical output payload bytes requested from the allocator.
    pub output_bytes: usize,
    /// Fixed RK4/identity scratch bytes.
    pub scratch_bytes: usize,
    /// Number of required terminal records.
    pub diagnostic_records: usize,
    /// Fixed report/error bytes.
    pub diagnostic_bytes: usize,
    /// Maximum canonical identity bytes hashed.
    pub identity_bytes: usize,
    /// Requested simultaneously live bytes before allocator slack.
    pub live_bytes: usize,
    /// Maximum checkpoint attempts for the full trajectory.
    pub polls: usize,
    /// Conservative deterministic work charged to the ambient Cx budget.
    pub work_units: u64,
    /// Admitted maximum steps/identity points between checkpoints.
    pub items_per_poll: usize,
}

/// RK4 stage associated with a typed callback or arithmetic refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StreamlineStage {
    /// First field evaluation at the current point.
    K1,
    /// Arithmetic constructing the second-stage point.
    K2Input,
    /// Second field evaluation.
    K2,
    /// Arithmetic constructing the third-stage point.
    K3Input,
    /// Third field evaluation.
    K3,
    /// Arithmetic constructing the fourth-stage point.
    K4Input,
    /// Fourth field evaluation.
    K4,
    /// Weighted RK4 state update.
    Advance,
}

/// Why a successfully published trajectory stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StreamlineTermination {
    /// All requested steps completed.
    StepsComplete,
    /// The zero-based attempted step would have left the admitted domain.
    DomainExit { attempted_step: usize },
    /// The zero-based attempted step rounded to the current point.
    Stagnated { attempted_step: usize },
}

impl StreamlineTermination {
    const fn identity_tag(self) -> u8 {
        match self {
            Self::StepsComplete => 1,
            Self::DomainExit { .. } => 2,
            Self::Stagnated { .. } => 3,
        }
    }

    const fn identity_step(self, requested_steps: usize) -> usize {
        match self {
            Self::StepsComplete => requested_steps,
            Self::DomainExit { attempted_step } | Self::Stagnated { attempted_step } => {
                attempted_step
            }
        }
    }
}

/// Terminal disposition of a scoped streamline transaction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StreamlineDisposition {
    /// Every requested step completed and published.
    Completed,
    /// A declared boundary/stagnation policy published a valid prefix.
    Terminated,
    /// Admission, arithmetic, callback, allocation, or a budget refused.
    Refused,
    /// Caller-owned cancellation was observed at a checkpoint.
    Cancelled,
}

/// Fixed structured evidence for one streamline attempt.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StreamlineReport {
    /// Exact semantic request.
    pub request: StreamlineRequest,
    /// Caller envelope presented for admission.
    pub operation_budget: Option<StreamlineBudget>,
    /// Checked plan, absent when request validation or plan arithmetic failed.
    pub plan: Option<StreamlinePlan>,
    /// Successfully appended RK4 steps.
    pub completed_steps: usize,
    /// Field callback attempts, including a panicking or non-finite callback.
    pub field_evaluations: usize,
    /// Points staged in the private output, including the seed.
    pub output_points: usize,
    /// Checkpoint attempts.
    pub polls: usize,
    /// Actual deterministic work charged or awaiting its terminal charge.
    pub work_units: u64,
    /// Allocator-reported point capacity in payload bytes.
    pub reserved_output_bytes: usize,
    /// Highest accounted simultaneously live byte count.
    pub peak_live_bytes: usize,
    /// Canonical identity bytes absorbed by the streaming hasher.
    pub identity_bytes_hashed: usize,
    /// Successful termination reason.
    pub termination: Option<StreamlineTermination>,
    /// Embedded local error estimate. Fixed-step RK4 v1 has none.
    pub error_estimate: Option<f64>,
    /// Number of retained terminal diagnostic records.
    pub diagnostic_records: usize,
    /// Whether a terminal success/refusal/cancellation record was finalized.
    pub terminal: bool,
    /// Whether the private trajectory passed the final checkpoint and escaped.
    pub published: bool,
    /// Terminal disposition.
    pub disposition: StreamlineDisposition,
    /// Domain-separated identity of the published trajectory.
    pub artifact_identity: Option<ContentHash>,
}

impl StreamlineReport {
    fn empty(request: StreamlineRequest, budget: StreamlineBudget) -> Self {
        Self {
            request,
            operation_budget: Some(budget),
            plan: None,
            completed_steps: 0,
            field_evaluations: 0,
            output_points: 0,
            polls: 0,
            work_units: 0,
            reserved_output_bytes: 0,
            peak_live_bytes: 0,
            identity_bytes_hashed: 0,
            termination: None,
            error_estimate: None,
            diagnostic_records: 0,
            terminal: false,
            published: false,
            disposition: StreamlineDisposition::Refused,
            artifact_identity: None,
        }
    }
}

/// Atomically published streamline points and their exact terminal report.
#[derive(Debug, Clone, PartialEq)]
pub struct StreamlineOutput {
    points: Vec<Vec2>,
    report: StreamlineReport,
}

impl StreamlineOutput {
    /// Published points in integration order, including the seed.
    #[must_use]
    pub fn points(&self) -> &[Vec2] {
        &self.points
    }

    /// Final resource, termination, cancellation, and identity report.
    #[must_use]
    pub const fn report(&self) -> &StreamlineReport {
        &self.report
    }

    /// Consume the wrapper into its point vector and report.
    #[must_use]
    pub fn into_parts(self) -> (Vec<Vec2>, StreamlineReport) {
        (self.points, self.report)
    }
}

/// Scoped streamline failure plus the terminal no-publication report.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StreamlineRunError {
    /// Typed root refusal.
    pub error: StreamlineError,
    /// Terminal counters and checked plan, when constructible.
    pub report: StreamlineReport,
}

impl core::fmt::Display for StreamlineRunError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            formatter,
            "{} (steps={}, calls={}, points={}, polls={}, published={})",
            self.error,
            self.report.completed_steps,
            self.report.field_evaluations,
            self.report.output_points,
            self.report.polls,
            self.report.published
        )
    }
}

impl std::error::Error for StreamlineRunError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.error)
    }
}

/// Typed refusal from fixed-step streamline admission or execution.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum StreamlineError {
    /// A seed component is non-finite.
    NonFiniteSeed {
        /// Cartesian component, `0..2`.
        component: usize,
        /// Rejected value.
        value: f64,
    },
    /// `dt` must be finite and nonzero; negative finite values are supported.
    InvalidStepSize {
        /// Rejected signed step.
        dt: f64,
    },
    /// The request names an unsupported method/receipt schema version.
    UnsupportedMethodVersion {
        /// Rejected version.
        version: u32,
    },
    /// One optional domain axis is non-finite, empty, or reversed.
    InvalidDomain {
        /// Cartesian axis, `0..2`.
        axis: usize,
        /// Rejected inclusive lower bound.
        lower: f64,
        /// Rejected inclusive upper bound.
        upper: f64,
    },
    /// The seed lies outside an otherwise valid closed domain.
    SeedOutsideDomain {
        /// Cartesian axis, `0..2`.
        axis: usize,
        /// Rejected seed coordinate.
        seed: f64,
        /// Inclusive lower bound.
        lower: f64,
        /// Inclusive upper bound.
        upper: f64,
    },
    /// Polling requires a positive step/identity-point stride.
    InvalidPollStride {
        /// Rejected stride.
        items_per_poll: usize,
    },
    /// Checked construction of the complete plan overflowed.
    PlanOverflow {
        /// Resource whose derivation was unrepresentable.
        resource: StreamlineResource,
    },
    /// One checked requirement exceeds its caller limit.
    OperationBudgetExceeded {
        /// Refused resource.
        resource: StreamlineResource,
        /// Checked requirement.
        required: u128,
        /// Caller-provided limit.
        limit: u128,
    },
    /// The ambient Cx cancellation/deadline/poll/cost contract refused.
    ExecutionBudgetRefused {
        /// Exact shared-accountant refusal.
        refusal: BudgetRefusal,
    },
    /// Storage for the complete private point vector could not be reserved.
    AllocationFailed {
        /// Required point capacity.
        required: usize,
    },
    /// The user callback unwound at a named RK4 stage; the panic was contained.
    CallbackPanicked {
        /// Zero-based attempted RK4 step.
        step: usize,
        /// Field-evaluation stage.
        stage: StreamlineStage,
    },
    /// A field callback returned a non-finite component.
    NonFiniteFieldValue {
        /// Zero-based attempted RK4 step.
        step: usize,
        /// Field-evaluation stage.
        stage: StreamlineStage,
        /// Cartesian component, `0..2`.
        component: usize,
        /// Rejected value.
        value: f64,
    },
    /// Finite inputs overflowed or otherwise produced non-finite RK4 geometry.
    NonFiniteIntermediate {
        /// Zero-based attempted RK4 step.
        step: usize,
        /// Arithmetic stage.
        stage: StreamlineStage,
        /// Cartesian component, `0..2`.
        component: usize,
        /// Rejected value.
        value: f64,
    },
    /// A candidate left the closed domain under `RefuseExit`.
    DomainExit {
        /// Zero-based attempted RK4 step.
        step: usize,
        /// Exact candidate coordinate bits.
        point_bits: [u64; 2],
    },
}

impl core::fmt::Display for StreamlineError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NonFiniteSeed { component, value } => {
                write!(
                    formatter,
                    "streamline seed component {component} is non-finite ({value})"
                )
            }
            Self::InvalidStepSize { dt } => write!(
                formatter,
                "streamline step size must be finite and nonzero (got {dt})"
            ),
            Self::UnsupportedMethodVersion { version } => write!(
                formatter,
                "streamline method/receipt version {version} is unsupported"
            ),
            Self::InvalidDomain { axis, lower, upper } => write!(
                formatter,
                "streamline domain axis {axis} must be finite and increasing (got {lower}..{upper})"
            ),
            Self::SeedOutsideDomain {
                axis,
                seed,
                lower,
                upper,
            } => write!(
                formatter,
                "streamline seed coordinate {seed} lies outside domain axis {axis} bounds {lower}..{upper}"
            ),
            Self::InvalidPollStride { items_per_poll } => write!(
                formatter,
                "streamline items-per-poll must be positive (got {items_per_poll})"
            ),
            Self::PlanOverflow { resource } => {
                write!(formatter, "streamline {resource} plan overflowed")
            }
            Self::OperationBudgetExceeded {
                resource,
                required,
                limit,
            } => write!(
                formatter,
                "streamline requires {required} {resource}, exceeding the explicit limit {limit}"
            ),
            Self::ExecutionBudgetRefused { refusal } => {
                write!(formatter, "streamline execution budget refused: {refusal}")
            }
            Self::AllocationFailed { required } => write!(
                formatter,
                "streamline could not reserve storage for {required} points"
            ),
            Self::CallbackPanicked { step, stage } => write!(
                formatter,
                "streamline callback panicked at step {step}, stage {stage:?}"
            ),
            Self::NonFiniteFieldValue {
                step,
                stage,
                component,
                value,
            } => write!(
                formatter,
                "streamline callback returned non-finite component {component} ({value}) at step {step}, stage {stage:?}"
            ),
            Self::NonFiniteIntermediate {
                step,
                stage,
                component,
                value,
            } => write!(
                formatter,
                "streamline arithmetic produced non-finite component {component} ({value}) at step {step}, stage {stage:?}"
            ),
            Self::DomainExit { step, point_bits } => write!(
                formatter,
                "streamline candidate at step {step} left the admitted domain (point bits {point_bits:?})"
            ),
        }
    }
}

impl std::error::Error for StreamlineError {}

/// Derive the exact conservative resource envelope for `request`.
///
/// # Errors
/// Returns a typed request-validation or checked-arithmetic refusal.
pub fn required_streamline_budget(
    request: StreamlineRequest,
) -> Result<StreamlineBudget, StreamlineError> {
    let plan = streamline_requirements(request)?;
    Ok(StreamlineBudget {
        step_limit: plan.steps,
        field_evaluation_limit: plan.field_evaluations,
        output_point_limit: plan.output_points,
        output_byte_limit: plan.output_bytes,
        scratch_byte_limit: plan.scratch_bytes,
        diagnostic_record_limit: plan.diagnostic_records,
        diagnostic_byte_limit: plan.diagnostic_bytes,
        live_byte_limit: plan.live_bytes,
        identity_byte_limit: plan.identity_bytes,
        poll_limit: plan.polls,
        work_unit_limit: plan.work_units,
    })
}

/// Build and admit the checked plan for an explicit streamline envelope.
///
/// # Errors
/// Returns a typed request, arithmetic, or per-resource budget refusal.
pub fn streamline_plan(
    request: StreamlineRequest,
    budget: StreamlineBudget,
) -> Result<StreamlinePlan, StreamlineError> {
    let plan = streamline_requirements(request)?;
    admit_streamline_plan(plan, budget)?;
    Ok(plan)
}

/// Run a fixed-step RK4 streamline under a caller-owned execution context.
///
/// The callback is invoked only after complete plan admission, the first Cx
/// checkpoint, and private output reservation. A callback unwind is contained
/// as a typed refusal; callback side effects remain the caller's responsibility.
///
/// # Errors
/// Returns [`StreamlineRunError`] with structured terminal evidence. No error,
/// cancellation, or unwind publishes a partial trajectory or identity.
pub fn streamline_with_cx(
    cx: &Cx<'_>,
    mut field: impl FnMut(Vec2) -> Vec2,
    request: StreamlineRequest,
    budget: StreamlineBudget,
) -> Result<StreamlineOutput, StreamlineRunError> {
    run_streamline_with(
        Some(cx),
        &mut field,
        request,
        budget,
        |_| {},
        |points, required| points.try_reserve_exact(required).map_err(|_| ()),
    )
}

/// No-authority compatibility wrapper for fixed-step RK4.
///
/// Valid finite requests preserve the original ordered seed-plus-step result.
/// Any invalid request, checked overflow, allocation refusal, non-finite field
/// result/arithmetic, domain refusal, or contained callback unwind returns an
/// empty vector. It has no caller-owned cancellation, receipt, or error channel;
/// authoritative callers must use [`streamline_with_cx`]. Negative finite `dt`
/// is supported, zero `dt` is rejected, and zero steps returns the finite seed.
#[must_use]
pub fn streamline(
    mut field: impl FnMut(Vec2) -> Vec2,
    seed: Vec2,
    dt: f64,
    steps: usize,
) -> Vec<Vec2> {
    let request = StreamlineRequest::dimensionless_rk4(seed, dt, steps, usize::MAX);
    let Ok(budget) = required_streamline_budget(request) else {
        return Vec::new();
    };
    run_streamline_with(
        None,
        &mut field,
        request,
        budget,
        |_| {},
        |points, required| points.try_reserve_exact(required).map_err(|_| ()),
    )
    .map(|output| output.points)
    .unwrap_or_default()
}

const STREAMLINE_IDENTITY_FIXED_PAYLOAD_BYTES: usize =
    size_of::<u32>() + 7 * size_of::<u8>() + 13 * size_of::<u64>();
const STREAMLINE_POINT_IDENTITY_BYTES: usize = 2 * size_of::<u64>();
const STREAMLINE_POINT_IDENTITY_WORK_UNITS: u64 = 16;
const STREAMLINE_RK4_FIELD_CALLS_PER_STEP: usize = 4;
const STREAMLINE_RK4_ARITHMETIC_WORK_PER_STEP: usize = 4;

fn checked_streamline_mul(
    left: usize,
    right: usize,
    resource: StreamlineResource,
) -> Result<usize, StreamlineError> {
    left.checked_mul(right)
        .ok_or(StreamlineError::PlanOverflow { resource })
}

fn checked_streamline_add(
    left: usize,
    right: usize,
    resource: StreamlineResource,
) -> Result<usize, StreamlineError> {
    left.checked_add(right)
        .ok_or(StreamlineError::PlanOverflow { resource })
}

fn checked_streamline_sum(
    values: impl IntoIterator<Item = usize>,
    resource: StreamlineResource,
) -> Result<usize, StreamlineError> {
    values.into_iter().try_fold(0usize, |sum, value| {
        checked_streamline_add(sum, value, resource)
    })
}

fn streamline_chunk_count(items: usize, items_per_poll: usize) -> Result<usize, StreamlineError> {
    if items == 0 {
        return Ok(0);
    }
    let preceding = items.checked_sub(1).ok_or(StreamlineError::PlanOverflow {
        resource: StreamlineResource::Polls,
    })?;
    checked_streamline_add(preceding / items_per_poll, 1, StreamlineResource::Polls)
}

fn validate_streamline_request(request: StreamlineRequest) -> Result<(), StreamlineError> {
    for (component, value) in request.seed.into_iter().enumerate() {
        if !value.is_finite() {
            return Err(StreamlineError::NonFiniteSeed { component, value });
        }
    }
    if !request.dt.is_finite() || request.dt == 0.0 {
        return Err(StreamlineError::InvalidStepSize { dt: request.dt });
    }
    if request.method_version != STREAMLINE_ARTIFACT_IDENTITY_VERSION {
        return Err(StreamlineError::UnsupportedMethodVersion {
            version: request.method_version,
        });
    }
    if request.items_per_poll == 0 {
        return Err(StreamlineError::InvalidPollStride {
            items_per_poll: request.items_per_poll,
        });
    }
    if let Some(domain) = request.domain {
        for axis in 0..2 {
            let extent = domain.upper[axis] - domain.lower[axis];
            if !(domain.lower[axis].is_finite()
                && domain.upper[axis].is_finite()
                && domain.upper[axis] > domain.lower[axis]
                && extent.is_finite())
            {
                return Err(StreamlineError::InvalidDomain {
                    axis,
                    lower: domain.lower[axis],
                    upper: domain.upper[axis],
                });
            }
            if request.seed[axis] < domain.lower[axis] || request.seed[axis] > domain.upper[axis] {
                return Err(StreamlineError::SeedOutsideDomain {
                    axis,
                    seed: request.seed[axis],
                    lower: domain.lower[axis],
                    upper: domain.upper[axis],
                });
            }
        }
    }
    Ok(())
}

fn streamline_requirements(request: StreamlineRequest) -> Result<StreamlinePlan, StreamlineError> {
    validate_streamline_request(request)?;
    u64::try_from(request.steps).map_err(|_| StreamlineError::PlanOverflow {
        resource: StreamlineResource::Steps,
    })?;
    let output_points = checked_streamline_add(request.steps, 1, StreamlineResource::OutputPoints)?;
    let field_evaluations = checked_streamline_mul(
        request.steps,
        STREAMLINE_RK4_FIELD_CALLS_PER_STEP,
        StreamlineResource::FieldEvaluations,
    )?;
    u64::try_from(field_evaluations).map_err(|_| StreamlineError::PlanOverflow {
        resource: StreamlineResource::FieldEvaluations,
    })?;
    u64::try_from(output_points).map_err(|_| StreamlineError::PlanOverflow {
        resource: StreamlineResource::OutputPoints,
    })?;
    let input_bytes = checked_streamline_add(
        size_of::<StreamlineRequest>(),
        size_of::<StreamlineBudget>(),
        StreamlineResource::LiveBytes,
    )?;
    let output_bytes = checked_streamline_mul(
        output_points,
        size_of::<Vec2>(),
        StreamlineResource::OutputBytes,
    )?;
    let rk4_scratch_bytes = size_of::<[Vec2; 6]>();
    let identity_scratch_bytes = checked_streamline_add(
        size_of::<DomainHasher>(),
        size_of::<u64>(),
        StreamlineResource::ScratchBytes,
    )?;
    let scratch_bytes = rk4_scratch_bytes.max(identity_scratch_bytes);
    let diagnostic_records = 1;
    let diagnostic_bytes = size_of::<StreamlineOutput>().max(size_of::<StreamlineRunError>());
    let point_identity_bytes = checked_streamline_mul(
        output_points,
        STREAMLINE_POINT_IDENTITY_BYTES,
        StreamlineResource::IdentityBytes,
    )?;
    let identity_bytes = checked_streamline_sum(
        [
            STREAMLINE_ARTIFACT_IDENTITY_DOMAIN.len(),
            STREAMLINE_IDENTITY_FIXED_PAYLOAD_BYTES,
            point_identity_bytes,
        ],
        StreamlineResource::IdentityBytes,
    )?;
    let live_bytes = checked_streamline_sum(
        [input_bytes, output_bytes, scratch_bytes, diagnostic_bytes],
        StreamlineResource::LiveBytes,
    )?;
    let step_chunks = streamline_chunk_count(request.steps, request.items_per_poll)?;
    let identity_chunks = streamline_chunk_count(output_points, request.items_per_poll)?;
    let polls =
        checked_streamline_sum([2, step_chunks, identity_chunks], StreamlineResource::Polls)?;
    let arithmetic_work = checked_streamline_mul(
        request.steps,
        STREAMLINE_RK4_ARITHMETIC_WORK_PER_STEP,
        StreamlineResource::WorkUnits,
    )?;
    let work_units_usize = checked_streamline_sum(
        [
            field_evaluations,
            arithmetic_work,
            request.steps,
            1,
            identity_bytes,
            1,
        ],
        StreamlineResource::WorkUnits,
    )?;
    let work_units =
        u64::try_from(work_units_usize).map_err(|_| StreamlineError::PlanOverflow {
            resource: StreamlineResource::WorkUnits,
        })?;

    Ok(StreamlinePlan {
        steps: request.steps,
        field_evaluations,
        output_points,
        input_bytes,
        output_bytes,
        scratch_bytes,
        diagnostic_records,
        diagnostic_bytes,
        identity_bytes,
        live_bytes,
        polls,
        work_units,
        items_per_poll: request.items_per_poll,
    })
}

fn admit_streamline_plan(
    plan: StreamlinePlan,
    budget: StreamlineBudget,
) -> Result<(), StreamlineError> {
    let requirements = [
        (
            StreamlineResource::Steps,
            usize_to_u128(plan.steps),
            usize_to_u128(budget.step_limit),
        ),
        (
            StreamlineResource::FieldEvaluations,
            usize_to_u128(plan.field_evaluations),
            usize_to_u128(budget.field_evaluation_limit),
        ),
        (
            StreamlineResource::OutputPoints,
            usize_to_u128(plan.output_points),
            usize_to_u128(budget.output_point_limit),
        ),
        (
            StreamlineResource::OutputBytes,
            usize_to_u128(plan.output_bytes),
            usize_to_u128(budget.output_byte_limit),
        ),
        (
            StreamlineResource::ScratchBytes,
            usize_to_u128(plan.scratch_bytes),
            usize_to_u128(budget.scratch_byte_limit),
        ),
        (
            StreamlineResource::DiagnosticRecords,
            usize_to_u128(plan.diagnostic_records),
            usize_to_u128(budget.diagnostic_record_limit),
        ),
        (
            StreamlineResource::DiagnosticBytes,
            usize_to_u128(plan.diagnostic_bytes),
            usize_to_u128(budget.diagnostic_byte_limit),
        ),
        (
            StreamlineResource::LiveBytes,
            usize_to_u128(plan.live_bytes),
            usize_to_u128(budget.live_byte_limit),
        ),
        (
            StreamlineResource::IdentityBytes,
            usize_to_u128(plan.identity_bytes),
            usize_to_u128(budget.identity_byte_limit),
        ),
        (
            StreamlineResource::Polls,
            usize_to_u128(plan.polls),
            usize_to_u128(budget.poll_limit),
        ),
        (
            StreamlineResource::WorkUnits,
            u128::from(plan.work_units),
            u128::from(budget.work_unit_limit),
        ),
    ];
    for (resource, required, limit) in requirements {
        if required > limit {
            return Err(StreamlineError::OperationBudgetExceeded {
                resource,
                required,
                limit,
            });
        }
    }
    Ok(())
}

fn terminal_streamline_error(
    error: StreamlineError,
    mut report: StreamlineReport,
) -> StreamlineRunError {
    report.diagnostic_records = 1;
    report.terminal = true;
    report.published = false;
    report.artifact_identity = None;
    report.disposition = if matches!(
        error,
        StreamlineError::ExecutionBudgetRefused {
            refusal: BudgetRefusal::Cancelled { .. }
        }
    ) {
        StreamlineDisposition::Cancelled
    } else {
        StreamlineDisposition::Refused
    };
    StreamlineRunError { error, report }
}

fn streamline_checkpoint<'clock, P>(
    cx: Option<&Cx<'clock>>,
    admitted: &mut Option<AdmittedBudget<'clock>>,
    report: &mut StreamlineReport,
    phase: &'static str,
    before_checkpoint: &mut P,
) -> Result<(), StreamlineError>
where
    P: FnMut(&StreamlineReport),
{
    let Some(cx) = cx else {
        return Ok(());
    };
    before_checkpoint(report);
    report.polls += 1;
    let Some(admitted) = admitted.as_mut() else {
        return Err(StreamlineError::PlanOverflow {
            resource: StreamlineResource::WorkUnits,
        });
    };
    admitted
        .checkpoint(phase, cx)
        .map_err(|refusal| StreamlineError::ExecutionBudgetRefused { refusal })
}

fn streamline_charge(
    admitted: &mut Option<AdmittedBudget<'_>>,
    pending_work: &mut u64,
    phase: &'static str,
) -> Result<(), StreamlineError> {
    if *pending_work == 0 {
        return Ok(());
    }
    if let Some(admitted) = admitted {
        admitted
            .charge_cost(phase, *pending_work)
            .map_err(|refusal| StreamlineError::ExecutionBudgetRefused { refusal })?;
    }
    *pending_work = 0;
    Ok(())
}

fn checked_streamline_intermediate(
    point: Vec2,
    step: usize,
    stage: StreamlineStage,
    report: &mut StreamlineReport,
    pending_work: &mut u64,
) -> Result<Vec2, StreamlineError> {
    report.work_units += 1;
    *pending_work += 1;
    for (component, value) in point.into_iter().enumerate() {
        if !value.is_finite() {
            return Err(StreamlineError::NonFiniteIntermediate {
                step,
                stage,
                component,
                value,
            });
        }
    }
    Ok(point)
}

fn evaluate_streamline_field<F>(
    field: &mut F,
    point: Vec2,
    step: usize,
    stage: StreamlineStage,
    report: &mut StreamlineReport,
    pending_work: &mut u64,
) -> Result<Vec2, StreamlineError>
where
    F: FnMut(Vec2) -> Vec2,
{
    report.field_evaluations += 1;
    report.work_units += 1;
    *pending_work += 1;
    let value = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| field(point)))
        .map_err(|_| StreamlineError::CallbackPanicked { step, stage })?;
    for (component, value_component) in value.into_iter().enumerate() {
        if !value_component.is_finite() {
            return Err(StreamlineError::NonFiniteFieldValue {
                step,
                stage,
                component,
                value: value_component,
            });
        }
    }
    Ok(value)
}

fn rk4_streamline_step<F>(
    field: &mut F,
    point: Vec2,
    dt: f64,
    step: usize,
    report: &mut StreamlineReport,
    pending_work: &mut u64,
) -> Result<Vec2, StreamlineError>
where
    F: FnMut(Vec2) -> Vec2,
{
    let k1 = evaluate_streamline_field(
        field,
        point,
        step,
        StreamlineStage::K1,
        report,
        pending_work,
    )?;
    let k2_point = checked_streamline_intermediate(
        [point[0] + 0.5 * dt * k1[0], point[1] + 0.5 * dt * k1[1]],
        step,
        StreamlineStage::K2Input,
        report,
        pending_work,
    )?;
    let k2 = evaluate_streamline_field(
        field,
        k2_point,
        step,
        StreamlineStage::K2,
        report,
        pending_work,
    )?;
    let k3_point = checked_streamline_intermediate(
        [point[0] + 0.5 * dt * k2[0], point[1] + 0.5 * dt * k2[1]],
        step,
        StreamlineStage::K3Input,
        report,
        pending_work,
    )?;
    let k3 = evaluate_streamline_field(
        field,
        k3_point,
        step,
        StreamlineStage::K3,
        report,
        pending_work,
    )?;
    let k4_point = checked_streamline_intermediate(
        [point[0] + dt * k3[0], point[1] + dt * k3[1]],
        step,
        StreamlineStage::K4Input,
        report,
        pending_work,
    )?;
    let k4 = evaluate_streamline_field(
        field,
        k4_point,
        step,
        StreamlineStage::K4,
        report,
        pending_work,
    )?;
    checked_streamline_intermediate(
        [
            point[0] + dt / 6.0 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]),
            point[1] + dt / 6.0 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]),
        ],
        step,
        StreamlineStage::Advance,
        report,
        pending_work,
    )
}

fn point_outside_streamline_domain(point: Vec2, domain: StreamlineDomain2) -> bool {
    (0..2).any(|axis| point[axis] < domain.lower[axis] || point[axis] > domain.upper[axis])
}

fn same_streamline_point(left: Vec2, right: Vec2) -> bool {
    left.map(f64::to_bits) == right.map(f64::to_bits)
}

fn streamline_identity_usize(
    value: usize,
    resource: StreamlineResource,
) -> Result<u64, StreamlineError> {
    u64::try_from(value).map_err(|_| StreamlineError::PlanOverflow { resource })
}

#[allow(clippy::too_many_lines)] // One transaction keeps callbacks, charging, and publication ordered.
fn run_streamline_with<'clock, F, P, R>(
    cx: Option<&Cx<'clock>>,
    field: &mut F,
    request: StreamlineRequest,
    budget: StreamlineBudget,
    mut before_checkpoint: P,
    mut reserve_output: R,
) -> Result<StreamlineOutput, StreamlineRunError>
where
    F: FnMut(Vec2) -> Vec2,
    P: FnMut(&StreamlineReport),
    R: FnMut(&mut Vec<Vec2>, usize) -> Result<(), ()>,
{
    let mut report = StreamlineReport::empty(request, budget);
    let plan = match streamline_requirements(request) {
        Ok(plan) => plan,
        Err(error) => return Err(terminal_streamline_error(error, report)),
    };
    report.plan = Some(plan);
    report.peak_live_bytes = plan.input_bytes + plan.diagnostic_bytes;
    if let Err(error) = admit_streamline_plan(plan, budget) {
        return Err(terminal_streamline_error(error, report));
    }

    let mut admitted = if let Some(cx) = cx {
        match AdmittedBudget::admit_ambient(cx, plan.work_units) {
            Ok(admitted) => Some(admitted),
            Err(refusal) => {
                return Err(terminal_streamline_error(
                    StreamlineError::ExecutionBudgetRefused { refusal },
                    report,
                ));
            }
        }
    } else {
        None
    };
    if let Err(error) = streamline_checkpoint(
        cx,
        &mut admitted,
        &mut report,
        "fs-viz.streamline.admission",
        &mut before_checkpoint,
    ) {
        return Err(terminal_streamline_error(error, report));
    }

    let mut points = Vec::new();
    if reserve_output(&mut points, plan.output_points).is_err() {
        return Err(terminal_streamline_error(
            StreamlineError::AllocationFailed {
                required: plan.output_points,
            },
            report,
        ));
    }
    let reserved_output_bytes = match points.capacity().checked_mul(size_of::<Vec2>()) {
        Some(bytes) => bytes,
        None => {
            return Err(terminal_streamline_error(
                StreamlineError::PlanOverflow {
                    resource: StreamlineResource::OutputBytes,
                },
                report,
            ));
        }
    };
    report.reserved_output_bytes = reserved_output_bytes;
    if points.capacity() < plan.output_points {
        return Err(terminal_streamline_error(
            StreamlineError::AllocationFailed {
                required: plan.output_points,
            },
            report,
        ));
    }
    if reserved_output_bytes > budget.output_byte_limit {
        return Err(terminal_streamline_error(
            StreamlineError::OperationBudgetExceeded {
                resource: StreamlineResource::OutputBytes,
                required: usize_to_u128(reserved_output_bytes),
                limit: usize_to_u128(budget.output_byte_limit),
            },
            report,
        ));
    }
    let non_output_live = plan.live_bytes - plan.output_bytes;
    let actual_live_bytes = match non_output_live.checked_add(reserved_output_bytes) {
        Some(bytes) => bytes,
        None => {
            return Err(terminal_streamline_error(
                StreamlineError::PlanOverflow {
                    resource: StreamlineResource::LiveBytes,
                },
                report,
            ));
        }
    };
    report.peak_live_bytes = actual_live_bytes;
    if actual_live_bytes > budget.live_byte_limit {
        return Err(terminal_streamline_error(
            StreamlineError::OperationBudgetExceeded {
                resource: StreamlineResource::LiveBytes,
                required: usize_to_u128(actual_live_bytes),
                limit: usize_to_u128(budget.live_byte_limit),
            },
            report,
        ));
    }

    let mut pending_work = 1u64;
    report.work_units = 1;
    points.push(request.seed);
    report.output_points = 1;
    let mut steps_since_poll = plan.items_per_poll;
    for step in 0..request.steps {
        if steps_since_poll == plan.items_per_poll {
            if let Err(error) = streamline_charge(
                &mut admitted,
                &mut pending_work,
                "fs-viz.streamline.step-chunk",
            ) {
                return Err(terminal_streamline_error(error, report));
            }
            if let Err(error) = streamline_checkpoint(
                cx,
                &mut admitted,
                &mut report,
                "fs-viz.streamline.step-chunk",
                &mut before_checkpoint,
            ) {
                return Err(terminal_streamline_error(error, report));
            }
            steps_since_poll = 0;
        }

        let point = points[points.len() - 1];
        let candidate = match rk4_streamline_step(
            field,
            point,
            request.dt,
            step,
            &mut report,
            &mut pending_work,
        ) {
            Ok(candidate) => candidate,
            Err(error) => {
                if let Err(budget_error) = streamline_charge(
                    &mut admitted,
                    &mut pending_work,
                    "fs-viz.streamline.step-refusal",
                ) {
                    return Err(terminal_streamline_error(budget_error, report));
                }
                return Err(terminal_streamline_error(error, report));
            }
        };

        if let Some(domain) = request.domain
            && point_outside_streamline_domain(candidate, domain)
        {
            match request.boundary_policy {
                StreamlineBoundaryPolicy::RefuseExit => {
                    if let Err(budget_error) = streamline_charge(
                        &mut admitted,
                        &mut pending_work,
                        "fs-viz.streamline.domain-refusal",
                    ) {
                        return Err(terminal_streamline_error(budget_error, report));
                    }
                    return Err(terminal_streamline_error(
                        StreamlineError::DomainExit {
                            step,
                            point_bits: candidate.map(f64::to_bits),
                        },
                        report,
                    ));
                }
                StreamlineBoundaryPolicy::StopBeforeExit => {
                    report.termination = Some(StreamlineTermination::DomainExit {
                        attempted_step: step,
                    });
                    break;
                }
            }
        }

        if same_streamline_point(candidate, point)
            && request.stagnation_policy == StreamlineStagnationPolicy::StopBeforeRepeat
        {
            report.termination = Some(StreamlineTermination::Stagnated {
                attempted_step: step,
            });
            break;
        }

        points.push(candidate);
        report.completed_steps += 1;
        report.output_points = points.len();
        report.work_units += 1;
        pending_work += 1;
        steps_since_poll += 1;
    }
    if report.termination.is_none() {
        report.termination = Some(StreamlineTermination::StepsComplete);
    }

    if let Err(error) = streamline_charge(
        &mut admitted,
        &mut pending_work,
        "fs-viz.streamline.step-finalize",
    ) {
        return Err(terminal_streamline_error(error, report));
    }
    if let Err(error) = streamline_checkpoint(
        cx,
        &mut admitted,
        &mut report,
        "fs-viz.streamline.identity",
        &mut before_checkpoint,
    ) {
        return Err(terminal_streamline_error(error, report));
    }

    let requested_steps = match streamline_identity_usize(request.steps, StreamlineResource::Steps)
    {
        Ok(value) => value,
        Err(error) => return Err(terminal_streamline_error(error, report)),
    };
    let completed_steps =
        match streamline_identity_usize(report.completed_steps, StreamlineResource::Steps) {
            Ok(value) => value,
            Err(error) => return Err(terminal_streamline_error(error, report)),
        };
    let field_evaluations = match streamline_identity_usize(
        report.field_evaluations,
        StreamlineResource::FieldEvaluations,
    ) {
        Ok(value) => value,
        Err(error) => return Err(terminal_streamline_error(error, report)),
    };
    let point_count =
        match streamline_identity_usize(points.len(), StreamlineResource::OutputPoints) {
            Ok(value) => value,
            Err(error) => return Err(terminal_streamline_error(error, report)),
        };
    let termination = match report.termination {
        Some(termination) => termination,
        None => {
            return Err(terminal_streamline_error(
                StreamlineError::PlanOverflow {
                    resource: StreamlineResource::IdentityBytes,
                },
                report,
            ));
        }
    };
    let termination_step = match streamline_identity_usize(
        termination.identity_step(request.steps),
        StreamlineResource::Steps,
    ) {
        Ok(value) => value,
        Err(error) => return Err(terminal_streamline_error(error, report)),
    };
    let (domain_present, domain) = match request.domain {
        Some(domain) => (1u8, domain),
        None => (
            0u8,
            StreamlineDomain2 {
                lower: [0.0; 2],
                upper: [0.0; 2],
            },
        ),
    };
    let mut hasher = DomainHasher::new(STREAMLINE_ARTIFACT_IDENTITY_DOMAIN);
    report.identity_bytes_hashed = STREAMLINE_ARTIFACT_IDENTITY_DOMAIN.len();
    hasher.update(&STREAMLINE_ARTIFACT_IDENTITY_VERSION.to_le_bytes());
    hasher.update(&[request.method.identity_tag()]);
    hasher.update(&[request.units.identity_tag()]);
    hasher.update(&[request.boundary_policy.identity_tag()]);
    hasher.update(&[request.stagnation_policy.identity_tag()]);
    for coordinate in request.seed {
        hasher.update(&coordinate.to_bits().to_le_bytes());
    }
    hasher.update(&request.dt.to_bits().to_le_bytes());
    hasher.update(&requested_steps.to_le_bytes());
    hasher.update(&[domain_present]);
    for coordinate in domain.lower.into_iter().chain(domain.upper) {
        hasher.update(&coordinate.to_bits().to_le_bytes());
    }
    hasher.update(&[termination.identity_tag()]);
    hasher.update(&termination_step.to_le_bytes());
    hasher.update(&completed_steps.to_le_bytes());
    hasher.update(&field_evaluations.to_le_bytes());
    hasher.update(&point_count.to_le_bytes());
    hasher.update(&[0]);
    hasher.update(&0u64.to_le_bytes());
    report.identity_bytes_hashed += STREAMLINE_IDENTITY_FIXED_PAYLOAD_BYTES;
    let identity_header_work = match u64::try_from(report.identity_bytes_hashed) {
        Ok(value) => value,
        Err(_) => {
            return Err(terminal_streamline_error(
                StreamlineError::PlanOverflow {
                    resource: StreamlineResource::WorkUnits,
                },
                report,
            ));
        }
    };
    report.work_units += identity_header_work;
    pending_work += identity_header_work;

    for (index, point) in points.iter().enumerate() {
        if index > 0 && index % plan.items_per_poll == 0 {
            if let Err(error) = streamline_charge(
                &mut admitted,
                &mut pending_work,
                "fs-viz.streamline.identity-chunk",
            ) {
                return Err(terminal_streamline_error(error, report));
            }
            if let Err(error) = streamline_checkpoint(
                cx,
                &mut admitted,
                &mut report,
                "fs-viz.streamline.identity-chunk",
                &mut before_checkpoint,
            ) {
                return Err(terminal_streamline_error(error, report));
            }
        }
        hasher.update(&point[0].to_bits().to_le_bytes());
        hasher.update(&point[1].to_bits().to_le_bytes());
        report.identity_bytes_hashed += STREAMLINE_POINT_IDENTITY_BYTES;
        report.work_units += STREAMLINE_POINT_IDENTITY_WORK_UNITS;
        pending_work += STREAMLINE_POINT_IDENTITY_WORK_UNITS;
    }
    let identity = hasher.finalize();
    report.work_units += 1;
    pending_work += 1;
    if let Err(error) = streamline_charge(
        &mut admitted,
        &mut pending_work,
        "fs-viz.streamline.identity-finalize",
    ) {
        return Err(terminal_streamline_error(error, report));
    }
    if let Err(error) = streamline_checkpoint(
        cx,
        &mut admitted,
        &mut report,
        "fs-viz.streamline.publish",
        &mut before_checkpoint,
    ) {
        return Err(terminal_streamline_error(error, report));
    }

    report.terminal = true;
    report.published = true;
    report.diagnostic_records = 1;
    report.disposition = if termination == StreamlineTermination::StepsComplete {
        StreamlineDisposition::Completed
    } else {
        StreamlineDisposition::Terminated
    };
    report.artifact_identity = Some(identity);
    Ok(StreamlineOutput { points, report })
}

/// The Morse type of a critical point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CriticalKind {
    /// A local minimum (Morse index 0).
    Minimum,
    /// A saddle (Morse index 1 in 2-D).
    Saddle,
    /// A local maximum (Morse index 2 in 2-D).
    Maximum,
    /// Degenerate or unclassifiable (a zero eigenvalue or invalid numeric
    /// input) — not a Morse critical point.
    Degenerate,
}

/// A classified critical point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CriticalPoint {
    /// The Morse type.
    pub kind: CriticalKind,
    /// The Morse index (number of negative Hessian eigenvalues).
    pub morse_index: usize,
}

/// Classify a critical point from its (symmetric) `2×2` Hessian: the Morse index
/// is the number of negative eigenvalues; a near-zero eigenvalue is degenerate.
/// Non-finite entries and non-finite or negative tolerances fail closed as
/// degenerate with no claimed negative eigenvalues.
#[must_use]
pub fn classify_hessian(hessian: [[f64; 2]; 2], tol: f64) -> CriticalPoint {
    if !tol.is_finite() || tol < 0.0 || hessian.iter().flatten().any(|value| !value.is_finite()) {
        return CriticalPoint {
            kind: CriticalKind::Degenerate,
            morse_index: 0,
        };
    }
    let (a, b, c) = (hessian[0][0], hessian[0][1], hessian[1][1]);
    let mean = f64::midpoint(a, c);
    let half_diff = f64::midpoint(a, -c);
    let scale = half_diff.abs().max(b.abs());
    let disc = if scale == 0.0 {
        0.0
    } else {
        let x = half_diff / scale;
        let y = b / scale;
        scale * (x * x + y * y).sqrt()
    };
    let (l1, l2) = (mean - disc, mean + disc);
    if l1.abs() <= tol || l2.abs() <= tol {
        return CriticalPoint {
            kind: CriticalKind::Degenerate,
            morse_index: usize::from(l1 < -tol) + usize::from(l2 < -tol),
        };
    }
    let morse_index = usize::from(l1 < 0.0) + usize::from(l2 < 0.0);
    let kind = match morse_index {
        0 => CriticalKind::Minimum,
        1 => CriticalKind::Saddle,
        _ => CriticalKind::Maximum,
    };
    CriticalPoint { kind, morse_index }
}

/// A scalar field sampled on a regular 2-D grid.
#[derive(Debug, Clone)]
pub struct Grid2 {
    nx: usize,
    ny: usize,
    lo: Vec2,
    hi: Vec2,
    values: Vec<f64>,
}

/// Admission failures for a regular 2-D scalar grid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Grid2Error {
    /// Both dimensions must contain at least two sample nodes.
    InvalidDimensions {
        /// Rejected `[nx, ny]` dimensions.
        dimensions: [usize; 2],
    },
    /// The dimension product does not fit in `usize`.
    NodeCountOverflow {
        /// Rejected `[nx, ny]` dimensions.
        dimensions: [usize; 2],
    },
    /// The explicit sampling budget is smaller than the grid.
    NodeBudgetExceeded {
        /// Required sample count.
        required: usize,
        /// Caller-provided limit.
        limit: usize,
    },
    /// A world-space interval is non-finite, non-increasing, or has a
    /// non-finite extent.
    InvalidBounds {
        /// Cartesian axis, `0..2`.
        axis: usize,
        /// Rejected lower bound.
        lower: f64,
        /// Rejected upper bound.
        upper: f64,
    },
    /// Distinct logical nodes collapse onto the same or a reversed floating
    /// coordinate at the requested resolution.
    UnrepresentableCoordinates {
        /// Cartesian axis, `0..2`.
        axis: usize,
        /// Earlier logical node index.
        first_index: usize,
        /// Earlier generated coordinate.
        first: f64,
        /// Later logical node index.
        second_index: usize,
        /// Later generated coordinate.
        second: f64,
    },
    /// A sampled scalar is non-finite.
    NonFiniteValue {
        /// Linear x-fastest node index.
        index: usize,
        /// Rejected value.
        value: f64,
    },
    /// The requested sample allocation could not be reserved.
    AllocationFailed {
        /// Requested node count.
        nodes: usize,
    },
}

impl core::fmt::Display for Grid2Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidDimensions { dimensions } => write!(
                f,
                "Grid2 dimensions must each be at least two (got {dimensions:?})"
            ),
            Self::NodeCountOverflow { dimensions } => write!(
                f,
                "Grid2 node count overflows for dimensions {dimensions:?}"
            ),
            Self::NodeBudgetExceeded { required, limit } => write!(
                f,
                "Grid2 requires {required} nodes, exceeding the explicit limit {limit}"
            ),
            Self::InvalidBounds { axis, lower, upper } => write!(
                f,
                "Grid2 axis {axis} bounds must be finite and increasing with finite extent (got {lower}..{upper})"
            ),
            Self::UnrepresentableCoordinates {
                axis,
                first_index,
                first,
                second_index,
                second,
            } => write!(
                f,
                "Grid2 axis {axis} nodes {first_index} and {second_index} are not strictly ordered ({first}, {second})"
            ),
            Self::NonFiniteValue { index, value } => {
                write!(f, "Grid2 value {index} is non-finite ({value})")
            }
            Self::AllocationFailed { nodes } => {
                write!(f, "Grid2 could not reserve storage for {nodes} nodes")
            }
        }
    }
}

impl std::error::Error for Grid2Error {}

/// A separately admitted resource in one Grid2 contour operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IsoContourResource {
    /// Logical grid cells.
    Cells,
    /// Horizontal plus vertical grid-edge visits.
    EdgeVisits,
    /// Constant-time exact-node ownership decisions.
    ExactOwnershipChecks,
    /// Strict-crossing interpolation attempts.
    Interpolations,
    /// Caller-visible crossing points.
    OutputCrossings,
    /// Logical crossing-payload bytes.
    OutputBytes,
    /// Fixed temporary state used while extracting and identifying output.
    ScratchBytes,
    /// Terminal diagnostic records.
    DiagnosticRecords,
    /// Fixed terminal report/error storage.
    DiagnosticBytes,
    /// Input, output, scratch, and diagnostic bytes simultaneously live.
    LiveBytes,
    /// Canonical artifact-identity bytes hashed.
    IdentityBytes,
    /// Bounded cancellation checkpoints.
    Polls,
    /// Deterministic scalar work charged to the ambient execution budget.
    WorkUnits,
}

impl core::fmt::Display for IsoContourResource {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        formatter.write_str(match self {
            Self::Cells => "cells",
            Self::EdgeVisits => "edge visits",
            Self::ExactOwnershipChecks => "exact ownership checks",
            Self::Interpolations => "interpolations",
            Self::OutputCrossings => "output crossings",
            Self::OutputBytes => "output bytes",
            Self::ScratchBytes => "scratch bytes",
            Self::DiagnosticRecords => "diagnostic records",
            Self::DiagnosticBytes => "diagnostic bytes",
            Self::LiveBytes => "simultaneously live bytes",
            Self::IdentityBytes => "identity bytes",
            Self::Polls => "polls",
            Self::WorkUnits => "work units",
        })
    }
}

/// Complete caller envelope for one cancellable Grid2 contour operation.
///
/// Count limits cover the conservative full-grid plan rather than only the
/// path eventually selected by sample values. This permits refusal before
/// output allocation or edge work.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IsoContourBudget {
    /// Maximum admitted logical cells.
    pub cell_limit: usize,
    /// Maximum admitted horizontal-plus-vertical edge visits.
    pub edge_visit_limit: usize,
    /// Maximum admitted exact-node ownership decisions.
    pub exact_ownership_limit: usize,
    /// Maximum admitted strict-crossing interpolation attempts.
    pub interpolation_limit: usize,
    /// Maximum caller-visible crossing points.
    pub crossing_limit: usize,
    /// Maximum logical bytes reserved for crossing payloads.
    pub output_byte_limit: usize,
    /// Maximum fixed scratch bytes.
    pub scratch_byte_limit: usize,
    /// Maximum terminal diagnostic records.
    pub diagnostic_record_limit: usize,
    /// Maximum fixed report/error bytes.
    pub diagnostic_byte_limit: usize,
    /// Maximum simultaneously live input, output, scratch, and diagnostic bytes.
    pub live_byte_limit: usize,
    /// Maximum canonical identity bytes hashed.
    pub identity_byte_limit: usize,
    /// Maximum cancellation checkpoint attempts.
    pub poll_limit: usize,
    /// Maximum deterministic work units.
    pub work_unit_limit: u64,
    /// Maximum edge or identity-point items between checkpoints.
    pub items_per_poll: usize,
}

/// Checked worst-case plan admitted before one Grid2 contour operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IsoContourPlan {
    /// Admitted `[nx, ny]` dimensions.
    pub dimensions: [usize; 2],
    /// Borrowed scalar nodes retained for the operation.
    pub nodes: usize,
    /// Logical cells in the grid.
    pub cells: usize,
    /// Horizontal plus vertical edges visited exactly once.
    pub edge_visits: usize,
    /// Conservative maximum exact-node ownership checks.
    pub exact_ownership_checks: usize,
    /// Conservative maximum strict-crossing interpolations.
    pub interpolations: usize,
    /// Caller-visible output capacity.
    pub max_crossings: usize,
    /// Bytes retained by the borrowed Grid2 value and payload.
    pub input_bytes: usize,
    /// Logical output payload bytes requested from the allocator.
    pub output_bytes: usize,
    /// Fixed scratch bytes used by the streaming identity hasher.
    pub scratch_bytes: usize,
    /// Number of terminal diagnostic records required.
    pub diagnostic_records: usize,
    /// Fixed terminal report/error bytes.
    pub diagnostic_bytes: usize,
    /// Maximum canonical identity bytes hashed.
    pub identity_bytes: usize,
    /// Requested simultaneously live bytes before allocator slack.
    pub live_bytes: usize,
    /// Maximum checkpoint attempts for the full-output path.
    pub polls: usize,
    /// Conservative scalar work charged to the ambient Cx budget.
    pub work_units: u64,
    /// Admitted maximum items between checkpoints.
    pub items_per_poll: usize,
}

/// Terminal state of a scoped contour operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IsoContourDisposition {
    /// The complete private result passed the final checkpoint and published.
    Completed,
    /// Admission, numerics, allocation, or an execution budget refused.
    Refused,
    /// Caller-owned cancellation was observed at a bounded checkpoint.
    Cancelled,
}

/// Fixed, structured telemetry for one contour attempt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IsoContourReport {
    /// Caller envelope presented for admission; absent only in internal
    /// edge-level test accounting that did not start an operation.
    pub operation_budget: Option<IsoContourBudget>,
    /// Checked plan, absent only when plan construction itself overflowed or
    /// the level/poll stride was invalid.
    pub plan: Option<IsoContourPlan>,
    /// Cells whose lower-left traversal point was reached.
    pub cell_visits: usize,
    /// Row-major grid nodes visited by the traversal.
    pub node_visits: usize,
    /// Edges examined.
    pub edge_visits: usize,
    /// Exact-node ownership decisions performed.
    pub exact_ownership_checks: usize,
    /// Exact-node candidates discarded because another incident edge owns them.
    pub unowned_exact_candidates: usize,
    /// Strict-crossing interpolations attempted.
    pub interpolations: usize,
    /// Exact-level points staged by their canonical owner edges.
    pub exact_crossings: usize,
    /// Strictly interpolated points staged in the private output.
    pub interpolated_crossings: usize,
    /// Points staged in the private output vector.
    pub crossings: usize,
    /// Cancellation checkpoint attempts.
    pub polls: usize,
    /// Actual deterministic work charged or awaiting the terminal charge.
    pub work_units: u64,
    /// Allocator-reported output capacity in payload bytes.
    pub reserved_output_bytes: usize,
    /// Highest accounted simultaneously live byte count.
    pub peak_live_bytes: usize,
    /// Canonical identity bytes absorbed by the streaming hasher.
    pub identity_bytes_hashed: usize,
    /// Number of retained terminal diagnostic records.
    pub diagnostic_records: usize,
    /// Whether a terminal success/refusal/cancellation record was finalized.
    pub terminal: bool,
    /// Whether the private output passed the final checkpoint and escaped.
    pub published: bool,
    /// Terminal disposition.
    pub disposition: IsoContourDisposition,
    /// Domain-separated identity of the published output; absent on failure.
    pub artifact_identity: Option<ContentHash>,
}

impl IsoContourReport {
    fn empty() -> Self {
        Self {
            operation_budget: None,
            plan: None,
            cell_visits: 0,
            node_visits: 0,
            edge_visits: 0,
            exact_ownership_checks: 0,
            unowned_exact_candidates: 0,
            interpolations: 0,
            exact_crossings: 0,
            interpolated_crossings: 0,
            crossings: 0,
            polls: 0,
            work_units: 0,
            reserved_output_bytes: 0,
            peak_live_bytes: 0,
            identity_bytes_hashed: 0,
            diagnostic_records: 0,
            terminal: false,
            published: false,
            disposition: IsoContourDisposition::Refused,
            artifact_identity: None,
        }
    }
}

/// Atomically published Grid2 crossings and the exact report for their run.
#[derive(Debug, Clone, PartialEq)]
pub struct IsoContourOutput {
    crossings: Vec<Vec2>,
    report: IsoContourReport,
}

impl IsoContourOutput {
    /// Published points in deterministic traversal order.
    #[must_use]
    pub fn crossings(&self) -> &[Vec2] {
        &self.crossings
    }

    /// Final resource, cancellation, and identity report.
    #[must_use]
    pub const fn report(&self) -> &IsoContourReport {
        &self.report
    }

    /// Consume the wrapper into its points and report.
    #[must_use]
    pub fn into_parts(self) -> (Vec<Vec2>, IsoContourReport) {
        (self.crossings, self.report)
    }
}

/// Scoped contour failure plus the terminal no-publication report.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IsoContourRunError {
    /// Typed root refusal.
    pub error: IsoContourError,
    /// Terminal counters and admitted plan, when one was constructible.
    pub report: IsoContourReport,
}

impl core::fmt::Display for IsoContourRunError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            formatter,
            "{} (edges={}, crossings={}, polls={}, published={})",
            self.error,
            self.report.edge_visits,
            self.report.crossings,
            self.report.polls,
            self.report.published
        )
    }
}

impl std::error::Error for IsoContourRunError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.error)
    }
}

/// Failures from bounded 2-D isocontour edge extraction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IsoContourError {
    /// The requested isovalue is non-finite.
    NonFiniteIso {
        /// Rejected isovalue.
        iso: f64,
    },
    /// At least one output crossing must be admitted by the explicit budget.
    ZeroCrossingLimit,
    /// A scoped operation must declare a positive item stride between polls.
    InvalidPollStride {
        /// Rejected item stride.
        items_per_poll: usize,
    },
    /// Checked construction of the complete operation plan overflowed.
    PlanOverflow {
        /// Resource whose derivation was not representable.
        resource: IsoContourResource,
    },
    /// A checked worst-case resource requirement exceeds its caller limit.
    OperationBudgetExceeded {
        /// Refused resource.
        resource: IsoContourResource,
        /// Checked requirement.
        required: u128,
        /// Caller-provided limit.
        limit: u128,
    },
    /// The ambient Cx cancellation/deadline/poll/cost contract refused.
    ExecutionBudgetRefused {
        /// Exact refusal returned by the shared fs-exec accountant.
        refusal: BudgetRefusal,
    },
    /// Extraction found more crossings than the caller admitted.
    CrossingBudgetExceeded {
        /// Caller-provided maximum crossing count.
        limit: usize,
    },
    /// A whole grid edge lies on the requested level, so its intersection is a
    /// segment that cannot be represented by a point-only result.
    CoincidentLevelEdge {
        /// First endpoint as `[i, j]`.
        first: [usize; 2],
        /// Second endpoint as `[i, j]`.
        second: [usize; 2],
    },
    /// Binary64 interpolation of a strict real edge crossing did not yield a
    /// point that remains strictly interior on every varying coordinate.
    UnrepresentableIntersection {
        /// First endpoint as `[i, j]` in deterministic traversal order.
        first: [usize; 2],
        /// Second endpoint as `[i, j]` in deterministic traversal order.
        second: [usize; 2],
        /// Exact coordinate bits of the first endpoint.
        first_point_bits: [u64; 2],
        /// Exact coordinate bits of the second endpoint.
        second_point_bits: [u64; 2],
        /// Exact sampled-value bits at the first endpoint.
        first_value_bits: u64,
        /// Exact sampled-value bits at the second endpoint.
        second_value_bits: u64,
        /// Exact requested isovalue bits.
        iso_bits: u64,
        /// Scaled distance from the first endpoint value to the isovalue.
        first_distance_bits: u64,
        /// Scaled distance from the second endpoint value to the isovalue.
        second_distance_bits: u64,
        /// Computed interpolation-parameter bits.
        interpolation_bits: u64,
        /// Computed point bits that collapsed onto or outside the edge.
        point_bits: [u64; 2],
        /// First coordinate that was not representably interior/on-edge.
        collapsed_axis: usize,
    },
    /// Storage for the next crossing could not be reserved.
    AllocationFailed {
        /// Number of crossings, including the one that could not be reserved.
        required: usize,
    },
    /// Interpolation produced a non-finite point or invalid interpolation
    /// fraction from otherwise finite inputs.
    NonFiniteGeometry,
}

impl core::fmt::Display for IsoContourError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NonFiniteIso { iso } => {
                write!(f, "isocontour level must be finite (got {iso})")
            }
            Self::ZeroCrossingLimit => {
                write!(f, "isocontour crossing limit must be positive")
            }
            Self::InvalidPollStride { items_per_poll } => write!(
                f,
                "isocontour items-per-poll must be positive (got {items_per_poll})"
            ),
            Self::PlanOverflow { resource } => {
                write!(f, "isocontour {resource} plan overflowed")
            }
            Self::OperationBudgetExceeded {
                resource,
                required,
                limit,
            } => write!(
                f,
                "isocontour requires {required} {resource}, exceeding the explicit limit {limit}"
            ),
            Self::ExecutionBudgetRefused { refusal } => {
                write!(f, "isocontour execution budget refused: {refusal}")
            }
            Self::CrossingBudgetExceeded { limit } => {
                write!(f, "isocontour exceeded the explicit {limit}-crossing limit")
            }
            Self::CoincidentLevelEdge { first, second } => write!(
                f,
                "isocontour edge {first:?}..{second:?} lies wholly on the requested level"
            ),
            Self::UnrepresentableIntersection {
                first,
                second,
                first_point_bits,
                second_point_bits,
                first_value_bits,
                second_value_bits,
                iso_bits,
                first_distance_bits,
                second_distance_bits,
                interpolation_bits,
                point_bits,
                collapsed_axis,
            } => write!(
                f,
                "isocontour edge {first:?}..{second:?} produced no admitted representably interior binary64 crossing on axis {collapsed_axis} (endpoint bits {first_point_bits:?}..{second_point_bits:?}, value bits {first_value_bits:016x}..{second_value_bits:016x}, iso {iso_bits:016x}, scaled distances {first_distance_bits:016x}/{second_distance_bits:016x}, t {interpolation_bits:016x}, point bits {point_bits:?})"
            ),
            Self::AllocationFailed { required } => write!(
                f,
                "isocontour could not reserve storage for {required} crossings"
            ),
            Self::NonFiniteGeometry => {
                write!(f, "isocontour interpolation produced non-finite geometry")
            }
        }
    }
}

impl std::error::Error for IsoContourError {}

impl Grid2 {
    /// Sample `f` on an admitted `nx × ny` regular grid spanning `[lo, hi]`.
    ///
    /// No sample is evaluated until dimensions, the explicit node budget, and
    /// world bounds have been admitted, sample storage has been reserved, and
    /// every generated axis coordinate is strictly increasing. Sampling order
    /// is row-major with x fastest.
    ///
    /// # Errors
    /// [`Grid2Error`] for invalid dimensions/bounds, count overflow,
    /// coordinate collapse, budget or allocation refusal, or the first
    /// non-finite sample.
    pub fn from_fn(
        nx: usize,
        ny: usize,
        lo: Vec2,
        hi: Vec2,
        node_limit: usize,
        mut f: impl FnMut(Vec2) -> f64,
    ) -> Result<Grid2, Grid2Error> {
        let node_count = validate_grid2_layout(nx, ny, lo, hi, node_limit)?;
        let mut values = Vec::new();
        values
            .try_reserve_exact(node_count)
            .map_err(|_| Grid2Error::AllocationFailed { nodes: node_count })?;
        validate_grid2_coordinates(nx, ny, lo, hi)?;
        for j in 0..ny {
            for i in 0..nx {
                let p = grid_point(nx, ny, lo, hi, i, j);
                let value = f(p);
                if !value.is_finite() {
                    return Err(Grid2Error::NonFiniteValue {
                        index: values.len(),
                        value,
                    });
                }
                values.push(value);
            }
        }
        Ok(Grid2 {
            nx,
            ny,
            lo,
            hi,
            values,
        })
    }

    /// The scalar value at grid node `(i, j)`.
    ///
    /// # Panics
    /// Panics when either index is outside the admitted grid.
    #[must_use]
    pub fn at(&self, i: usize, j: usize) -> f64 {
        assert!(i < self.nx && j < self.ny, "Grid2 node index out of bounds");
        self.values[j * self.nx + i]
    }

    /// The world coordinate of grid node `(i, j)`.
    ///
    /// # Panics
    /// Callers must keep both indices inside the admitted grid. Out-of-range
    /// coordinates are not part of the public grid domain.
    #[must_use]
    pub fn point(&self, i: usize, j: usize) -> Vec2 {
        assert!(i < self.nx && j < self.ny, "Grid2 node index out of bounds");
        grid_point(self.nx, self.ny, self.lo, self.hi, i, j)
    }

    /// Derive the exact conservative resource envelope required for a caller's
    /// output cap and cancellation stride.
    ///
    /// The derived envelope admits the complete grid independently of sampled
    /// signs: every edge could require either exact ownership or interpolation,
    /// and every output slot could participate in artifact identity.
    ///
    /// # Errors
    /// [`IsoContourError::ZeroCrossingLimit`],
    /// [`IsoContourError::InvalidPollStride`], or checked plan overflow.
    pub fn required_isocontour_budget(
        &self,
        crossing_limit: usize,
        items_per_poll: usize,
    ) -> Result<IsoContourBudget, IsoContourError> {
        let plan = contour_requirements(self, crossing_limit, items_per_poll)?;
        Ok(IsoContourBudget {
            cell_limit: plan.cells,
            edge_visit_limit: plan.edge_visits,
            exact_ownership_limit: plan.exact_ownership_checks,
            interpolation_limit: plan.interpolations,
            crossing_limit: plan.max_crossings,
            output_byte_limit: plan.output_bytes,
            scratch_byte_limit: plan.scratch_bytes,
            diagnostic_record_limit: plan.diagnostic_records,
            diagnostic_byte_limit: plan.diagnostic_bytes,
            live_byte_limit: plan.live_bytes,
            identity_byte_limit: plan.identity_bytes,
            poll_limit: plan.polls,
            work_unit_limit: plan.work_units,
            items_per_poll: plan.items_per_poll,
        })
    }

    /// Build and admit the checked plan for an explicit contour envelope.
    ///
    /// # Errors
    /// Returns a typed plan or resource-budget refusal before allocation or
    /// edge traversal.
    pub fn isocontour_plan(
        &self,
        budget: IsoContourBudget,
    ) -> Result<IsoContourPlan, IsoContourError> {
        let plan = contour_requirements(self, budget.crossing_limit, budget.items_per_poll)?;
        admit_contour_plan(plan, budget)?;
        Ok(plan)
    }

    /// Compatibility entry point for bounded isocontour edge crossings.
    ///
    /// Strict sign changes yield one scaled-interpolation point. An exact-level
    /// endpoint is emitted once at its original bits even when several incident
    /// edges meet there. A wholly coincident edge is refused because its
    /// intersection is a segment, not a unique point. Traversal and first-node
    /// ownership are row-major, with the positive-x edge before positive-y.
    ///
    /// # Errors
    /// [`IsoContourError`] for a non-finite level, zero or exceeded output
    /// budget, a coincident level edge, an unrepresentable strict intersection,
    /// allocation refusal, or non-finite interpolated geometry. Every failure
    /// returns no partial crossing vector. This compatibility entry point
    /// derives a complete exact envelope but has no caller-owned cancellation;
    /// production scoped work should use [`Grid2::isocontour_crossings_with_cx`].
    pub fn isocontour_crossings(
        &self,
        iso: f64,
        crossing_limit: usize,
    ) -> Result<Vec<Vec2>, IsoContourError> {
        if !iso.is_finite() {
            return Err(IsoContourError::NonFiniteIso { iso });
        }
        let budget = self.required_isocontour_budget(crossing_limit, usize::MAX)?;
        run_isocontour_with(
            self,
            None,
            iso,
            budget,
            |_| {},
            |output, required| output.try_reserve_exact(required).map_err(|_| ()),
        )
        .map(|output| output.crossings)
        .map_err(|error| error.error)
    }

    /// Run a fully planned contour extraction under a caller-owned execution
    /// context. Output remains private until the final checkpoint and identity
    /// seal both succeed.
    ///
    /// # Errors
    /// Returns [`IsoContourRunError`] with the checked plan and exact terminal
    /// counters. Cancellation, deadline, poll-quota, and cost-quota refusals
    /// never publish a partial vector or artifact identity.
    pub fn isocontour_crossings_with_cx(
        &self,
        cx: &Cx<'_>,
        iso: f64,
        budget: IsoContourBudget,
    ) -> Result<IsoContourOutput, IsoContourRunError> {
        run_isocontour_with(
            self,
            Some(cx),
            iso,
            budget,
            |_| {},
            |output, required| output.try_reserve_exact(required).map_err(|_| ()),
        )
    }
}

const ISO_CONTOUR_IDENTITY_FIXED_PAYLOAD_BYTES: usize =
    size_of::<u32>() + 2 * size_of::<u64>() + 4 * size_of::<u64>() + 2 * size_of::<u64>();
const ISO_CONTOUR_POINT_IDENTITY_BYTES: usize = 2 * size_of::<u64>();
// Two canonical u64 coordinate-bit patterns are charged byte-for-byte.
const ISO_CONTOUR_POINT_IDENTITY_WORK_UNITS: u64 = 16;

fn usize_to_u128(value: usize) -> u128 {
    u128::try_from(value).unwrap_or(u128::MAX)
}

fn checked_plan_mul(
    left: usize,
    right: usize,
    resource: IsoContourResource,
) -> Result<usize, IsoContourError> {
    left.checked_mul(right)
        .ok_or(IsoContourError::PlanOverflow { resource })
}

fn checked_plan_add(
    left: usize,
    right: usize,
    resource: IsoContourResource,
) -> Result<usize, IsoContourError> {
    left.checked_add(right)
        .ok_or(IsoContourError::PlanOverflow { resource })
}

fn checked_plan_sum(
    values: impl IntoIterator<Item = usize>,
    resource: IsoContourResource,
) -> Result<usize, IsoContourError> {
    values
        .into_iter()
        .try_fold(0usize, |sum, value| checked_plan_add(sum, value, resource))
}

fn chunk_count(items: usize, items_per_poll: usize) -> Result<usize, IsoContourError> {
    let preceding = items.checked_sub(1).ok_or(IsoContourError::PlanOverflow {
        resource: IsoContourResource::Polls,
    })?;
    checked_plan_add(preceding / items_per_poll, 1, IsoContourResource::Polls)
}

#[allow(clippy::too_many_lines)] // One checked derivation keeps every coupled resource auditable.
fn contour_requirements(
    grid: &Grid2,
    crossing_limit: usize,
    items_per_poll: usize,
) -> Result<IsoContourPlan, IsoContourError> {
    if crossing_limit == 0 {
        return Err(IsoContourError::ZeroCrossingLimit);
    }
    if items_per_poll == 0 {
        return Err(IsoContourError::InvalidPollStride { items_per_poll });
    }
    u64::try_from(grid.nx).map_err(|_| IsoContourError::PlanOverflow {
        resource: IsoContourResource::IdentityBytes,
    })?;
    u64::try_from(grid.ny).map_err(|_| IsoContourError::PlanOverflow {
        resource: IsoContourResource::IdentityBytes,
    })?;
    u64::try_from(crossing_limit).map_err(|_| IsoContourError::PlanOverflow {
        resource: IsoContourResource::OutputCrossings,
    })?;

    let nx_edges = grid.nx - 1;
    let ny_edges = grid.ny - 1;
    let nodes = checked_plan_mul(grid.nx, grid.ny, IsoContourResource::LiveBytes)?;
    let cells = checked_plan_mul(nx_edges, ny_edges, IsoContourResource::Cells)?;
    let horizontal = checked_plan_mul(nx_edges, grid.ny, IsoContourResource::EdgeVisits)?;
    let vertical = checked_plan_mul(grid.nx, ny_edges, IsoContourResource::EdgeVisits)?;
    let edge_visits = checked_plan_add(horizontal, vertical, IsoContourResource::EdgeVisits)?;
    let node_payload_bytes = checked_plan_mul(
        grid.values.capacity(),
        size_of::<f64>(),
        IsoContourResource::LiveBytes,
    )?;
    let input_bytes = checked_plan_add(
        size_of::<Grid2>(),
        node_payload_bytes,
        IsoContourResource::LiveBytes,
    )?;
    let output_bytes = checked_plan_mul(
        crossing_limit,
        size_of::<Vec2>(),
        IsoContourResource::OutputBytes,
    )?;
    let scratch_bytes = checked_plan_add(
        size_of::<DomainHasher>(),
        size_of::<u64>(),
        IsoContourResource::ScratchBytes,
    )?;
    let diagnostic_records = 1;
    let diagnostic_bytes = size_of::<IsoContourOutput>().max(size_of::<IsoContourRunError>());
    let point_identity_bytes = checked_plan_mul(
        crossing_limit,
        ISO_CONTOUR_POINT_IDENTITY_BYTES,
        IsoContourResource::IdentityBytes,
    )?;
    let identity_bytes = checked_plan_sum(
        [
            ISO_CONTOUR_ARTIFACT_IDENTITY_DOMAIN.len(),
            ISO_CONTOUR_IDENTITY_FIXED_PAYLOAD_BYTES,
            point_identity_bytes,
        ],
        IsoContourResource::IdentityBytes,
    )?;
    let live_bytes = checked_plan_sum(
        [input_bytes, output_bytes, scratch_bytes, diagnostic_bytes],
        IsoContourResource::LiveBytes,
    )?;
    let edge_chunks = chunk_count(edge_visits, items_per_poll)?;
    let identity_chunks = chunk_count(crossing_limit, items_per_poll)?;
    let polls = checked_plan_sum([2, edge_chunks, identity_chunks], IsoContourResource::Polls)?;
    let work_units_usize = checked_plan_sum(
        [
            nodes,
            cells,
            edge_visits,
            edge_visits,
            edge_visits,
            crossing_limit,
            identity_bytes,
            1,
        ],
        IsoContourResource::WorkUnits,
    )?;
    let work_units =
        u64::try_from(work_units_usize).map_err(|_| IsoContourError::PlanOverflow {
            resource: IsoContourResource::WorkUnits,
        })?;

    Ok(IsoContourPlan {
        dimensions: [grid.nx, grid.ny],
        nodes,
        cells,
        edge_visits,
        exact_ownership_checks: edge_visits,
        interpolations: edge_visits,
        max_crossings: crossing_limit,
        input_bytes,
        output_bytes,
        scratch_bytes,
        diagnostic_records,
        diagnostic_bytes,
        identity_bytes,
        live_bytes,
        polls,
        work_units,
        items_per_poll,
    })
}

fn admit_contour_plan(
    plan: IsoContourPlan,
    budget: IsoContourBudget,
) -> Result<(), IsoContourError> {
    let requirements = [
        (
            IsoContourResource::Cells,
            usize_to_u128(plan.cells),
            usize_to_u128(budget.cell_limit),
        ),
        (
            IsoContourResource::EdgeVisits,
            usize_to_u128(plan.edge_visits),
            usize_to_u128(budget.edge_visit_limit),
        ),
        (
            IsoContourResource::ExactOwnershipChecks,
            usize_to_u128(plan.exact_ownership_checks),
            usize_to_u128(budget.exact_ownership_limit),
        ),
        (
            IsoContourResource::Interpolations,
            usize_to_u128(plan.interpolations),
            usize_to_u128(budget.interpolation_limit),
        ),
        (
            IsoContourResource::OutputBytes,
            usize_to_u128(plan.output_bytes),
            usize_to_u128(budget.output_byte_limit),
        ),
        (
            IsoContourResource::ScratchBytes,
            usize_to_u128(plan.scratch_bytes),
            usize_to_u128(budget.scratch_byte_limit),
        ),
        (
            IsoContourResource::DiagnosticRecords,
            usize_to_u128(plan.diagnostic_records),
            usize_to_u128(budget.diagnostic_record_limit),
        ),
        (
            IsoContourResource::DiagnosticBytes,
            usize_to_u128(plan.diagnostic_bytes),
            usize_to_u128(budget.diagnostic_byte_limit),
        ),
        (
            IsoContourResource::LiveBytes,
            usize_to_u128(plan.live_bytes),
            usize_to_u128(budget.live_byte_limit),
        ),
        (
            IsoContourResource::IdentityBytes,
            usize_to_u128(plan.identity_bytes),
            usize_to_u128(budget.identity_byte_limit),
        ),
        (
            IsoContourResource::Polls,
            usize_to_u128(plan.polls),
            usize_to_u128(budget.poll_limit),
        ),
        (
            IsoContourResource::WorkUnits,
            u128::from(plan.work_units),
            u128::from(budget.work_unit_limit),
        ),
    ];
    for (resource, required, limit) in requirements {
        if required > limit {
            return Err(IsoContourError::OperationBudgetExceeded {
                resource,
                required,
                limit,
            });
        }
    }
    Ok(())
}

fn terminal_contour_error(
    error: IsoContourError,
    mut report: IsoContourReport,
) -> IsoContourRunError {
    report.diagnostic_records = 1;
    report.terminal = true;
    report.published = false;
    report.artifact_identity = None;
    report.disposition = if matches!(
        error,
        IsoContourError::ExecutionBudgetRefused {
            refusal: BudgetRefusal::Cancelled { .. }
        }
    ) {
        IsoContourDisposition::Cancelled
    } else {
        IsoContourDisposition::Refused
    };
    IsoContourRunError { error, report }
}

fn contour_checkpoint<'clock, P>(
    cx: Option<&Cx<'clock>>,
    admitted: &mut Option<AdmittedBudget<'clock>>,
    report: &mut IsoContourReport,
    phase: &'static str,
    before_checkpoint: &mut P,
) -> Result<(), IsoContourError>
where
    P: FnMut(&IsoContourReport),
{
    let Some(cx) = cx else {
        return Ok(());
    };
    before_checkpoint(report);
    report.polls += 1;
    let Some(admitted) = admitted.as_mut() else {
        return Err(IsoContourError::PlanOverflow {
            resource: IsoContourResource::WorkUnits,
        });
    };
    admitted
        .checkpoint(phase, cx)
        .map_err(|refusal| IsoContourError::ExecutionBudgetRefused { refusal })
}

fn contour_charge(
    admitted: &mut Option<AdmittedBudget<'_>>,
    pending_work: &mut u64,
    phase: &'static str,
) -> Result<(), IsoContourError> {
    if *pending_work == 0 {
        return Ok(());
    }
    if let Some(admitted) = admitted {
        admitted
            .charge_cost(phase, *pending_work)
            .map_err(|refusal| IsoContourError::ExecutionBudgetRefused { refusal })?;
    }
    *pending_work = 0;
    Ok(())
}

#[allow(clippy::too_many_lines)] // One atomic transaction keeps charging and publication ordered.
fn run_isocontour_with<'clock, P, R>(
    grid: &Grid2,
    cx: Option<&Cx<'clock>>,
    iso: f64,
    budget: IsoContourBudget,
    mut before_checkpoint: P,
    mut reserve_output: R,
) -> Result<IsoContourOutput, IsoContourRunError>
where
    P: FnMut(&IsoContourReport),
    R: FnMut(&mut Vec<Vec2>, usize) -> Result<(), ()>,
{
    let mut report = IsoContourReport::empty();
    report.operation_budget = Some(budget);
    if !iso.is_finite() {
        return Err(terminal_contour_error(
            IsoContourError::NonFiniteIso { iso },
            report,
        ));
    }
    let plan = match contour_requirements(grid, budget.crossing_limit, budget.items_per_poll) {
        Ok(plan) => plan,
        Err(error) => return Err(terminal_contour_error(error, report)),
    };
    report.plan = Some(plan);
    report.peak_live_bytes = plan.input_bytes + plan.diagnostic_bytes;
    if let Err(error) = admit_contour_plan(plan, budget) {
        return Err(terminal_contour_error(error, report));
    }

    let mut admitted = if let Some(cx) = cx {
        match AdmittedBudget::admit_ambient(cx, plan.work_units) {
            Ok(admitted) => Some(admitted),
            Err(refusal) => {
                return Err(terminal_contour_error(
                    IsoContourError::ExecutionBudgetRefused { refusal },
                    report,
                ));
            }
        }
    } else {
        None
    };
    if let Err(error) = contour_checkpoint(
        cx,
        &mut admitted,
        &mut report,
        "fs-viz.isocontour.admission",
        &mut before_checkpoint,
    ) {
        return Err(terminal_contour_error(error, report));
    }

    let mut crossings = Vec::new();
    if reserve_output(&mut crossings, plan.max_crossings).is_err() {
        return Err(terminal_contour_error(
            IsoContourError::AllocationFailed {
                required: plan.max_crossings,
            },
            report,
        ));
    }
    let reserved_output_bytes = match crossings.capacity().checked_mul(size_of::<Vec2>()) {
        Some(bytes) => bytes,
        None => {
            return Err(terminal_contour_error(
                IsoContourError::PlanOverflow {
                    resource: IsoContourResource::OutputBytes,
                },
                report,
            ));
        }
    };
    report.reserved_output_bytes = reserved_output_bytes;
    if crossings.capacity() < plan.max_crossings {
        return Err(terminal_contour_error(
            IsoContourError::AllocationFailed {
                required: plan.max_crossings,
            },
            report,
        ));
    }
    if reserved_output_bytes > budget.output_byte_limit {
        return Err(terminal_contour_error(
            IsoContourError::OperationBudgetExceeded {
                resource: IsoContourResource::OutputBytes,
                required: usize_to_u128(reserved_output_bytes),
                limit: usize_to_u128(budget.output_byte_limit),
            },
            report,
        ));
    }
    let non_output_live = plan.live_bytes - plan.output_bytes;
    let actual_live_bytes = match non_output_live.checked_add(reserved_output_bytes) {
        Some(bytes) => bytes,
        None => {
            return Err(terminal_contour_error(
                IsoContourError::PlanOverflow {
                    resource: IsoContourResource::LiveBytes,
                },
                report,
            ));
        }
    };
    report.peak_live_bytes = actual_live_bytes;
    if actual_live_bytes > budget.live_byte_limit {
        return Err(terminal_contour_error(
            IsoContourError::OperationBudgetExceeded {
                resource: IsoContourResource::LiveBytes,
                required: usize_to_u128(actual_live_bytes),
                limit: usize_to_u128(budget.live_byte_limit),
            },
            report,
        ));
    }

    let mut pending_work = 0u64;
    let mut edges_since_poll = plan.items_per_poll;
    for j in 0..grid.ny {
        for i in 0..grid.nx {
            let node_has_edge = i + 1 < grid.nx || j + 1 < grid.ny;
            if node_has_edge && edges_since_poll == plan.items_per_poll {
                if let Err(error) = contour_charge(
                    &mut admitted,
                    &mut pending_work,
                    "fs-viz.isocontour.edge-chunk",
                ) {
                    return Err(terminal_contour_error(error, report));
                }
                if let Err(error) = contour_checkpoint(
                    cx,
                    &mut admitted,
                    &mut report,
                    "fs-viz.isocontour.edge-chunk",
                    &mut before_checkpoint,
                ) {
                    return Err(terminal_contour_error(error, report));
                }
                edges_since_poll = 0;
            }

            report.node_visits += 1;
            report.work_units += 1;
            pending_work += 1;
            if i + 1 < grid.nx && j + 1 < grid.ny {
                report.cell_visits += 1;
                report.work_units += 1;
                pending_work += 1;
            }

            let (point, value) = (grid.point(i, j), grid.at(i, j));
            if i + 1 < grid.nx {
                report.edge_visits += 1;
                report.work_units += 1;
                pending_work += 1;
                let crossing = match edge_crossing(
                    [i, j],
                    point,
                    value,
                    [i + 1, j],
                    grid.point(i + 1, j),
                    grid.at(i + 1, j),
                    iso,
                    &mut report,
                    &mut pending_work,
                ) {
                    Ok(crossing) => crossing,
                    Err(error) => {
                        if let Err(budget_error) = contour_charge(
                            &mut admitted,
                            &mut pending_work,
                            "fs-viz.isocontour.edge-refusal",
                        ) {
                            return Err(terminal_contour_error(budget_error, report));
                        }
                        return Err(terminal_contour_error(error, report));
                    }
                };
                if let Some(crossing) = crossing {
                    if let Err(error) =
                        push_crossing(&mut crossings, crossing, plan.max_crossings, &mut report)
                    {
                        if let Err(budget_error) = contour_charge(
                            &mut admitted,
                            &mut pending_work,
                            "fs-viz.isocontour.output-refusal",
                        ) {
                            return Err(terminal_contour_error(budget_error, report));
                        }
                        return Err(terminal_contour_error(error, report));
                    }
                    report.work_units += 1;
                    pending_work += 1;
                }
                edges_since_poll += 1;
            }
            if j + 1 < grid.ny {
                if edges_since_poll == plan.items_per_poll {
                    if let Err(error) = contour_charge(
                        &mut admitted,
                        &mut pending_work,
                        "fs-viz.isocontour.edge-chunk",
                    ) {
                        return Err(terminal_contour_error(error, report));
                    }
                    if let Err(error) = contour_checkpoint(
                        cx,
                        &mut admitted,
                        &mut report,
                        "fs-viz.isocontour.edge-chunk",
                        &mut before_checkpoint,
                    ) {
                        return Err(terminal_contour_error(error, report));
                    }
                    edges_since_poll = 0;
                }
                report.edge_visits += 1;
                report.work_units += 1;
                pending_work += 1;
                let crossing = match edge_crossing(
                    [i, j],
                    point,
                    value,
                    [i, j + 1],
                    grid.point(i, j + 1),
                    grid.at(i, j + 1),
                    iso,
                    &mut report,
                    &mut pending_work,
                ) {
                    Ok(crossing) => crossing,
                    Err(error) => {
                        if let Err(budget_error) = contour_charge(
                            &mut admitted,
                            &mut pending_work,
                            "fs-viz.isocontour.edge-refusal",
                        ) {
                            return Err(terminal_contour_error(budget_error, report));
                        }
                        return Err(terminal_contour_error(error, report));
                    }
                };
                if let Some(crossing) = crossing {
                    if let Err(error) =
                        push_crossing(&mut crossings, crossing, plan.max_crossings, &mut report)
                    {
                        if let Err(budget_error) = contour_charge(
                            &mut admitted,
                            &mut pending_work,
                            "fs-viz.isocontour.output-refusal",
                        ) {
                            return Err(terminal_contour_error(budget_error, report));
                        }
                        return Err(terminal_contour_error(error, report));
                    }
                    report.work_units += 1;
                    pending_work += 1;
                }
                edges_since_poll += 1;
            }
        }
    }

    if let Err(error) = contour_charge(
        &mut admitted,
        &mut pending_work,
        "fs-viz.isocontour.edge-finalize",
    ) {
        return Err(terminal_contour_error(error, report));
    }
    if let Err(error) = contour_checkpoint(
        cx,
        &mut admitted,
        &mut report,
        "fs-viz.isocontour.identity",
        &mut before_checkpoint,
    ) {
        return Err(terminal_contour_error(error, report));
    }

    let nx = match u64::try_from(grid.nx) {
        Ok(value) => value,
        Err(_) => {
            return Err(terminal_contour_error(
                IsoContourError::PlanOverflow {
                    resource: IsoContourResource::IdentityBytes,
                },
                report,
            ));
        }
    };
    let ny = match u64::try_from(grid.ny) {
        Ok(value) => value,
        Err(_) => {
            return Err(terminal_contour_error(
                IsoContourError::PlanOverflow {
                    resource: IsoContourResource::IdentityBytes,
                },
                report,
            ));
        }
    };
    let crossing_count = match u64::try_from(crossings.len()) {
        Ok(value) => value,
        Err(_) => {
            return Err(terminal_contour_error(
                IsoContourError::PlanOverflow {
                    resource: IsoContourResource::IdentityBytes,
                },
                report,
            ));
        }
    };
    let mut hasher = DomainHasher::new(ISO_CONTOUR_ARTIFACT_IDENTITY_DOMAIN);
    report.identity_bytes_hashed = ISO_CONTOUR_ARTIFACT_IDENTITY_DOMAIN.len();
    hasher.update(&ISO_CONTOUR_ARTIFACT_IDENTITY_VERSION.to_le_bytes());
    hasher.update(&nx.to_le_bytes());
    hasher.update(&ny.to_le_bytes());
    for coordinate in grid.lo.into_iter().chain(grid.hi) {
        hasher.update(&coordinate.to_bits().to_le_bytes());
    }
    hasher.update(&iso.to_bits().to_le_bytes());
    hasher.update(&crossing_count.to_le_bytes());
    report.identity_bytes_hashed += ISO_CONTOUR_IDENTITY_FIXED_PAYLOAD_BYTES;
    let identity_header_work = match u64::try_from(report.identity_bytes_hashed) {
        Ok(value) => value,
        Err(_) => {
            return Err(terminal_contour_error(
                IsoContourError::PlanOverflow {
                    resource: IsoContourResource::WorkUnits,
                },
                report,
            ));
        }
    };
    report.work_units += identity_header_work;
    pending_work += identity_header_work;

    for (index, point) in crossings.iter().enumerate() {
        if index > 0 && index % plan.items_per_poll == 0 {
            if let Err(error) = contour_charge(
                &mut admitted,
                &mut pending_work,
                "fs-viz.isocontour.identity-chunk",
            ) {
                return Err(terminal_contour_error(error, report));
            }
            if let Err(error) = contour_checkpoint(
                cx,
                &mut admitted,
                &mut report,
                "fs-viz.isocontour.identity-chunk",
                &mut before_checkpoint,
            ) {
                return Err(terminal_contour_error(error, report));
            }
        }
        hasher.update(&point[0].to_bits().to_le_bytes());
        hasher.update(&point[1].to_bits().to_le_bytes());
        report.identity_bytes_hashed += ISO_CONTOUR_POINT_IDENTITY_BYTES;
        report.work_units += ISO_CONTOUR_POINT_IDENTITY_WORK_UNITS;
        pending_work += ISO_CONTOUR_POINT_IDENTITY_WORK_UNITS;
    }
    let identity = hasher.finalize();
    report.work_units += 1;
    pending_work += 1;
    if let Err(error) = contour_charge(
        &mut admitted,
        &mut pending_work,
        "fs-viz.isocontour.identity-finalize",
    ) {
        return Err(terminal_contour_error(error, report));
    }
    if let Err(error) = contour_checkpoint(
        cx,
        &mut admitted,
        &mut report,
        "fs-viz.isocontour.publish",
        &mut before_checkpoint,
    ) {
        return Err(terminal_contour_error(error, report));
    }

    report.terminal = true;
    report.published = true;
    report.diagnostic_records = 1;
    report.disposition = IsoContourDisposition::Completed;
    report.artifact_identity = Some(identity);
    Ok(IsoContourOutput { crossings, report })
}

fn validate_grid2_layout(
    nx: usize,
    ny: usize,
    lo: Vec2,
    hi: Vec2,
    node_limit: usize,
) -> Result<usize, Grid2Error> {
    let dimensions = [nx, ny];
    if dimensions.into_iter().any(|dimension| dimension < 2) {
        return Err(Grid2Error::InvalidDimensions { dimensions });
    }
    let node_count = nx
        .checked_mul(ny)
        .ok_or(Grid2Error::NodeCountOverflow { dimensions })?;
    if node_count > node_limit {
        return Err(Grid2Error::NodeBudgetExceeded {
            required: node_count,
            limit: node_limit,
        });
    }
    for axis in 0..2 {
        let extent = hi[axis] - lo[axis];
        if !(lo[axis].is_finite()
            && hi[axis].is_finite()
            && hi[axis] > lo[axis]
            && extent.is_finite())
        {
            return Err(Grid2Error::InvalidBounds {
                axis,
                lower: lo[axis],
                upper: hi[axis],
            });
        }
    }
    Ok(node_count)
}

fn validate_grid2_coordinates(nx: usize, ny: usize, lo: Vec2, hi: Vec2) -> Result<(), Grid2Error> {
    for (axis, nodes) in [nx, ny].into_iter().enumerate() {
        let mut previous = grid_axis_point(nodes, lo[axis], hi[axis], 0);
        for index in 1..nodes {
            let current = grid_axis_point(nodes, lo[axis], hi[axis], index);
            if !current.is_finite() || current <= previous {
                return Err(Grid2Error::UnrepresentableCoordinates {
                    axis,
                    first_index: index - 1,
                    first: previous,
                    second_index: index,
                    second: current,
                });
            }
            previous = current;
        }
    }
    Ok(())
}

#[derive(Debug, Clone, Copy)]
enum EdgeCrossing2 {
    Exact(Vec2),
    Interpolated(Vec2),
}

#[allow(clippy::too_many_arguments)] // Endpoint evidence stays explicit for typed refusals.
fn edge_crossing(
    a_index: [usize; 2],
    a: Vec2,
    va: f64,
    b_index: [usize; 2],
    b: Vec2,
    vb: f64,
    iso: f64,
    report: &mut IsoContourReport,
    pending_work: &mut u64,
) -> Result<Option<EdgeCrossing2>, IsoContourError> {
    let a_exact = va == iso;
    let b_exact = vb == iso;
    if a_exact && b_exact {
        return Err(IsoContourError::CoincidentLevelEdge {
            first: a_index,
            second: b_index,
        });
    }
    if a_exact {
        report.exact_ownership_checks += 1;
        report.work_units += 1;
        *pending_work += 1;
        if edge_owns_exact_node(a_index, b_index, a_index) {
            return Ok(Some(EdgeCrossing2::Exact(a)));
        }
        report.unowned_exact_candidates += 1;
        return Ok(None);
    }
    if b_exact {
        report.exact_ownership_checks += 1;
        report.work_units += 1;
        *pending_work += 1;
        if edge_owns_exact_node(a_index, b_index, b_index) {
            return Ok(Some(EdgeCrossing2::Exact(b)));
        }
        report.unowned_exact_candidates += 1;
        return Ok(None);
    }
    if (va < iso) == (vb < iso) {
        return Ok(None);
    }
    report.interpolations += 1;
    report.work_units += 1;
    *pending_work += 1;
    let scale = va.abs().max(vb.abs()).max(iso.abs());
    let scaled_iso = iso / scale;
    let a_distance = (va / scale - scaled_iso).abs();
    let b_distance = (vb / scale - scaled_iso).abs();
    let t = a_distance / (a_distance + b_distance);
    if !t.is_finite() {
        return Err(IsoContourError::NonFiniteGeometry);
    }
    let point = [
        (b[0] - a[0]).mul_add(t, a[0]),
        (b[1] - a[1]).mul_add(t, a[1]),
    ];
    if point.into_iter().any(|coordinate| !coordinate.is_finite()) {
        return Err(IsoContourError::NonFiniteGeometry);
    }
    if let Some(collapsed_axis) = first_unrepresentable_intersection_axis(a, b, point) {
        return Err(IsoContourError::UnrepresentableIntersection {
            first: a_index,
            second: b_index,
            first_point_bits: a.map(f64::to_bits),
            second_point_bits: b.map(f64::to_bits),
            first_value_bits: va.to_bits(),
            second_value_bits: vb.to_bits(),
            iso_bits: iso.to_bits(),
            first_distance_bits: a_distance.to_bits(),
            second_distance_bits: b_distance.to_bits(),
            interpolation_bits: t.to_bits(),
            point_bits: point.map(f64::to_bits),
            collapsed_axis,
        });
    }
    Ok(Some(EdgeCrossing2::Interpolated(point)))
}

/// Return whether this positive-axis edge is the deterministic owner of an
/// exact-level endpoint.
///
/// Edge traversal is row-major by its first endpoint, with positive x before
/// positive y. Consequently the first incident edge of `(i, j)` is its edge
/// from the row below when `j > 0`, otherwise its edge from the left when
/// `i > 0`, and otherwise the origin's positive-x edge. Selecting that edge
/// directly makes ownership constant-time and needs no global deduplication
/// search or marker allocation.
fn edge_owns_exact_node(first: [usize; 2], second: [usize; 2], exact: [usize; 2]) -> bool {
    let owner = if exact[1] > 0 {
        ([exact[0], exact[1] - 1], exact)
    } else if exact[0] > 0 {
        ([exact[0] - 1, 0], exact)
    } else {
        (exact, [1, 0])
    };
    (first, second) == owner
}

fn first_unrepresentable_intersection_axis(a: Vec2, b: Vec2, point: Vec2) -> Option<usize> {
    for axis in 0..2 {
        match a[axis].total_cmp(&b[axis]) {
            core::cmp::Ordering::Equal => {
                if point[axis].to_bits() != a[axis].to_bits() {
                    return Some(axis);
                }
            }
            core::cmp::Ordering::Less => {
                if !(a[axis].total_cmp(&point[axis]).is_lt()
                    && point[axis].total_cmp(&b[axis]).is_lt())
                {
                    return Some(axis);
                }
            }
            core::cmp::Ordering::Greater => {
                if !(b[axis].total_cmp(&point[axis]).is_lt()
                    && point[axis].total_cmp(&a[axis]).is_lt())
                {
                    return Some(axis);
                }
            }
        }
    }
    None
}

fn push_crossing(
    crossings: &mut Vec<Vec2>,
    crossing: EdgeCrossing2,
    crossing_limit: usize,
    report: &mut IsoContourReport,
) -> Result<(), IsoContourError> {
    let (point, exact_increment, interpolated_increment) = match crossing {
        EdgeCrossing2::Exact(point) => (point, 1, 0),
        EdgeCrossing2::Interpolated(point) => (point, 0, 1),
    };
    if crossings.len() == crossing_limit {
        return Err(IsoContourError::CrossingBudgetExceeded {
            limit: crossing_limit,
        });
    }
    let next_exact_crossings = report.exact_crossings.checked_add(exact_increment).ok_or(
        IsoContourError::PlanOverflow {
            resource: IsoContourResource::OutputCrossings,
        },
    )?;
    let next_interpolated_crossings = report
        .interpolated_crossings
        .checked_add(interpolated_increment)
        .ok_or(IsoContourError::PlanOverflow {
            resource: IsoContourResource::OutputCrossings,
        })?;
    let required = crossings
        .len()
        .checked_add(1)
        .ok_or(IsoContourError::AllocationFailed {
            required: usize::MAX,
        })?;
    debug_assert!(crossings.capacity() >= required);
    crossings.push(point);
    report.exact_crossings = next_exact_crossings;
    report.interpolated_crossings = next_interpolated_crossings;
    report.crossings = crossings.len();
    debug_assert_eq!(
        report.crossings,
        report.exact_crossings + report.interpolated_crossings
    );
    Ok(())
}

fn grid_point(nx: usize, ny: usize, lo: Vec2, hi: Vec2, i: usize, j: usize) -> Vec2 {
    [
        grid_axis_point(nx, lo[0], hi[0], i),
        grid_axis_point(ny, lo[1], hi[1], j),
    ]
}

fn grid_axis_point(nodes: usize, lower: f64, upper: f64, index: usize) -> f64 {
    if index == 0 {
        lower
    } else if index + 1 == nodes {
        upper
    } else {
        let t = index as f64 / (nodes - 1) as f64;
        (upper - lower).mul_add(t, lower)
    }
}

#[cfg(test)]
mod viz_fault_tests {
    use super::*;
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    use std::panic::{AssertUnwindSafe, catch_unwind};

    fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                gate,
                arena,
                StreamKey {
                    seed: 0x1_50_C0,
                    kernel_id: 5,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    fn affine_grid() -> Grid2 {
        Grid2::from_fn(9, 9, [-1.0; 2], [1.0; 2], 81, |point| {
            point[0] + 0.25 * point[1] - 0.03
        })
        .expect("finite affine fixture")
    }

    fn streamline_fixture() -> (StreamlineRequest, StreamlineBudget) {
        let request = StreamlineRequest::dimensionless_rk4([0.0; 2], 0.125, 8, 1);
        let budget = required_streamline_budget(request).expect("checked streamline fixture");
        (request, budget)
    }

    #[test]
    fn g4_streamline_injected_allocation_refusal_precedes_callback_work() {
        let (request, budget) = streamline_fixture();
        let gate = CancelGate::new();
        let calls = std::cell::Cell::new(0usize);
        let refusal = with_cx(&gate, |cx| {
            let mut field = |_| {
                calls.set(calls.get() + 1);
                [1.0, 0.0]
            };
            run_streamline_with(
                Some(cx),
                &mut field,
                request,
                budget,
                |_| {},
                |_, _| Err(()),
            )
            .expect_err("injected point reservation refusal")
        });
        assert_eq!(
            refusal.error,
            StreamlineError::AllocationFailed {
                required: request.steps + 1,
            }
        );
        assert_eq!(calls.get(), 0);
        assert_eq!(refusal.report.field_evaluations, 0);
        assert_eq!(refusal.report.output_points, 0);
        assert_eq!(refusal.report.reserved_output_bytes, 0);
        assert!(refusal.report.terminal && !refusal.report.published);
    }

    #[test]
    fn g4_streamline_injected_midstage_cancellation_is_retryable() {
        let (request, budget) = streamline_fixture();
        let gate = CancelGate::new();
        let mut requested = false;
        let cancellation = with_cx(&gate, |cx| {
            let mut field = |_| [1.0, -0.5];
            run_streamline_with(
                Some(cx),
                &mut field,
                request,
                budget,
                |report| {
                    if !requested && report.completed_steps >= 2 {
                        requested = true;
                        gate.request();
                    }
                },
                |points, required| points.try_reserve_exact(required).map_err(|_| ()),
            )
            .expect_err("injected midstage cancellation")
        });
        assert!(requested);
        assert!(matches!(
            cancellation.error,
            StreamlineError::ExecutionBudgetRefused {
                refusal: BudgetRefusal::Cancelled { .. }
            }
        ));
        assert_eq!(cancellation.report.completed_steps, 2);
        assert_eq!(
            cancellation.report.disposition,
            StreamlineDisposition::Cancelled
        );
        assert!(!cancellation.report.published);
        assert!(cancellation.report.artifact_identity.is_none());

        let retry_gate = CancelGate::new();
        let retry = with_cx(&retry_gate, |cx| {
            streamline_with_cx(cx, |_| [1.0, -0.5], request, budget)
                .expect("retry after cancellation")
        });
        let direct_gate = CancelGate::new();
        let direct = with_cx(&direct_gate, |cx| {
            streamline_with_cx(cx, |_| [1.0, -0.5], request, budget)
                .expect("direct after cancellation")
        });
        assert_eq!(retry, direct);
    }

    #[test]
    fn g4_injected_allocation_refusal_precedes_edge_work() {
        let grid = affine_grid();
        let budget = grid
            .required_isocontour_budget(32, 4)
            .expect("checked fixture budget");
        let gate = CancelGate::new();
        let refusal = with_cx(&gate, |cx| {
            run_isocontour_with(&grid, Some(cx), 0.0, budget, |_| {}, |_, _| Err(()))
                .expect_err("injected reserve refusal")
        });
        assert_eq!(
            refusal.error,
            IsoContourError::AllocationFailed { required: 32 }
        );
        assert_eq!(refusal.report.edge_visits, 0);
        assert_eq!(refusal.report.crossings, 0);
        assert_eq!(refusal.report.reserved_output_bytes, 0);
        assert!(refusal.report.terminal && !refusal.report.published);
    }

    #[test]
    fn g4_injected_midstage_cancellation_is_atomic_and_retryable() {
        let grid = affine_grid();
        let budget = grid
            .required_isocontour_budget(32, 2)
            .expect("checked fixture budget");
        let gate = CancelGate::new();
        let mut requested = false;
        let cancellation = with_cx(&gate, |cx| {
            run_isocontour_with(
                &grid,
                Some(cx),
                0.0,
                budget,
                |report| {
                    if !requested && report.edge_visits >= 2 {
                        requested = true;
                        gate.request();
                    }
                },
                |output, required| output.try_reserve_exact(required).map_err(|_| ()),
            )
            .expect_err("injected stage cancellation")
        });
        assert!(requested);
        assert!(matches!(
            cancellation.error,
            IsoContourError::ExecutionBudgetRefused {
                refusal: BudgetRefusal::Cancelled { .. }
            }
        ));
        assert_eq!(cancellation.report.edge_visits, 2);
        assert_eq!(
            cancellation.report.disposition,
            IsoContourDisposition::Cancelled
        );
        assert!(!cancellation.report.published);
        assert!(cancellation.report.artifact_identity.is_none());

        let retry_gate = CancelGate::new();
        let retry = with_cx(&retry_gate, |cx| {
            grid.isocontour_crossings_with_cx(cx, 0.0, budget)
                .expect("retry succeeds")
        });
        let direct_gate = CancelGate::new();
        let direct = with_cx(&direct_gate, |cx| {
            grid.isocontour_crossings_with_cx(cx, 0.0, budget)
                .expect("direct succeeds")
        });
        assert_eq!(retry, direct);
    }

    #[test]
    fn g4_checkpoint_panic_drops_private_output_and_retry_is_identical() {
        let grid = affine_grid();
        let budget = grid
            .required_isocontour_budget(32, 2)
            .expect("checked fixture budget");
        let panic_gate = CancelGate::new();
        let panicked = catch_unwind(AssertUnwindSafe(|| {
            with_cx(&panic_gate, |cx| {
                let _ = run_isocontour_with(
                    &grid,
                    Some(cx),
                    0.0,
                    budget,
                    |report| {
                        assert!(report.edge_visits < 2, "injected checkpoint panic");
                    },
                    |output, required| output.try_reserve_exact(required).map_err(|_| ()),
                );
            });
        }));
        assert!(panicked.is_err());

        let retry_gate = CancelGate::new();
        let retry = with_cx(&retry_gate, |cx| {
            grid.isocontour_crossings_with_cx(cx, 0.0, budget)
                .expect("retry after unwind")
        });
        let direct_gate = CancelGate::new();
        let direct = with_cx(&direct_gate, |cx| {
            grid.isocontour_crossings_with_cx(cx, 0.0, budget)
                .expect("direct after unwind")
        });
        assert_eq!(retry, direct);
    }

    #[test]
    fn g3_interpolation_is_endpoint_reversal_and_translation_equivariant() {
        fn interpolate(a: Vec2, va: f64, b: Vec2, vb: f64) -> Vec2 {
            let mut report = IsoContourReport::empty();
            let mut pending = 0;
            match edge_crossing([0, 0], a, va, [1, 0], b, vb, 0.0, &mut report, &mut pending)
                .expect("representable strict crossing")
                .expect("opposite signs cross")
            {
                EdgeCrossing2::Exact(_) => panic!("strict crossing was relabeled exact"),
                EdgeCrossing2::Interpolated(point) => point,
            }
        }

        let forward = interpolate([0.0, 5.0], -1.0, [2.0, 5.0], 1.0);
        let reversed = interpolate([2.0, 5.0], 1.0, [0.0, 5.0], -1.0);
        let translated = interpolate([8.0, -3.0], -1.0, [10.0, -3.0], 1.0);
        assert_eq!(
            forward.map(f64::to_bits),
            [1.0_f64.to_bits(), 5.0_f64.to_bits()]
        );
        assert_eq!(reversed.map(f64::to_bits), forward.map(f64::to_bits));
        assert_eq!(
            translated.map(f64::to_bits),
            [9.0_f64.to_bits(), (-3.0_f64).to_bits()]
        );
    }
}
