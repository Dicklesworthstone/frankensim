//! fs-neuroshape-e2e — NeuroShapeCert: certified facts about a neural implicit.
//! Layer: L5 (LUMEN).
//!
//! # The campaign
//!
//! A learned neural SDF renders a shape, but gives no guarantees: how far can a
//! sphere-tracing ray step without tunneling through a thin feature, and which
//! topology facts are actually certified? This campaign proves a safe step and
//! the existence of at least one enclosed negative component. It deliberately
//! makes no exact component-count claim, composing crates never designed to meet:
//!
//! - **The field** ([`fs_rep_neural`]): a small `tanh`-MLP SDF whose
//!   spectral-normalized effective form is `≈ 2.12·Σ tanh(3(±coord − 0.7)) + 6.5`
//!   — provably negative near the origin, provably positive on a surrounding ring.
//! - **A certified Lipschitz constant** — `L = Π Uᵢ` (product of outward
//!   spectral-norm upper bounds). A degenerate IBP enclosure at the origin
//!   supplies a certified lower sign margin; dividing it downward by `L`
//!   yields a sphere-trace step that CANNOT tunnel through the surface.
//! - **A topology certificate by interval arithmetic** — the network's sound
//!   Interval Bound Propagation (`eval_interval`) proves a central box is
//!   strictly inside (`hi < 0`) and that the FOUR edge strips of a bounding box
//!   are strictly outside (`lo > 0`). Those strips tile the box boundary into a
//!   CLOSED frame (corners overlap), so the component meeting the negative
//!   central box cannot cross it: at least one component is proven to exist and
//!   be ENCLOSED — a proof, not a mesh. (Discrete ring boxes would leave angular
//!   gaps and prove no enclosure theorem.)
//! - **Typed topology evidence**: the negative central box and positive closed
//!   frame construct a [`CertifiedEnclosedComponentExists`]. Its public
//!   [`ComponentCountEvidence`] reports only the global lower bound `>= 1`;
//!   disconnected interior wells or negative exterior regions remain possible.
//! - **A curvature cross-check** ([`fs_viz`]): the origin has a positive-definite
//!   finite-difference Hessian. Without a certified zero gradient this is not a
//!   critical-point or minimum theorem, and never a component-count proof.
//!   `isocontour_crossings` separately localizes the sampled zero set.
//! Same-build bit-deterministic replay; retained cross-ISA evidence is required
//! before attaching a portable G5 receipt.

use fs_rep_neural::{
    Layer, MLP_ACTIVATION_SEMANTICS, MLP_ACTIVATION_SEMANTICS_VERSION, MLP_ACTIVATION_ULP_BUDGET,
    MlpSdf, NeuralFieldIdentity, SAFE_STEP_POLICY, SAFE_STEP_POLICY_VERSION, SafeStepDerivation,
    derive_safe_step,
};
use fs_exec::BudgetRefusal;
use fs_viz::{CriticalKind, Grid2, Grid2Error, IsoContourError, Vec2, classify_hessian};
use std::fmt;

/// Version of the public component-evidence semantics carried by
/// [`NeuroShapeReport`].
///
/// Version 1 means that enclosed-component evidence carries only a global
/// lower bound, while an exact component count remains unavailable. Adapters
/// serializing these fields must carry this value so consumers can refuse
/// layouts whose topology semantics they do not understand.
pub const NEUROSHAPE_COMPONENT_EVIDENCE_SCHEMA_VERSION: u32 = 1;

/// Stable schema version of [`NeuroShapeReport::surface_localization`] and of
/// every wire code derived from it (bead frankensim-o33vo).
///
/// Version 1 pins [`SurfaceLocalizationStatus`] codes `1..=8`, the
/// [`LocalizationStage`] codes `1..=2` (`0` meaning "no single stage"), and
/// the [`LocalizationDiagnostic`] codes `1..=18`. Serializers must carry this
/// value and consumers must refuse an unrecognized code instead of
/// reinterpreting it.
pub const NEUROSHAPE_LOCALIZATION_SCHEMA_VERSION: u32 = 1;

/// Sentinel stored in [`StageDetail`] numeric slots that carry no meaning for
/// the active [`LocalizationDiagnostic`].
pub const LOCALIZATION_DETAIL_UNDEFINED: u64 = u64::MAX;
/// `u32` companion of [`LOCALIZATION_DETAIL_UNDEFINED`].
pub const LOCALIZATION_DETAIL_UNDEFINED_U32: u32 = u32::MAX;

/// Exact fs-viz stage that produced a sampled zero-set localization outcome.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LocalizationStage {
    /// `fs_viz::Grid2::from_fn` admitted sampling.
    GridConstruction,
    /// `Grid2::isocontour_crossings` bounded extraction.
    IsoContourExtraction,
}

impl LocalizationStage {
    /// Stable wire code under [`NEUROSHAPE_LOCALIZATION_SCHEMA_VERSION`].
    #[must_use]
    pub const fn code(self) -> u32 {
        match self {
            Self::GridConstruction => 1,
            Self::IsoContourExtraction => 2,
        }
    }
}

/// Stable coarse outcome class of sampled zero-set localization.
///
/// The discriminants ARE the wire codes under
/// [`NEUROSHAPE_LOCALIZATION_SCHEMA_VERSION`]; `0` stays reserved as the
/// no-claim convention shared with the other NeuroShape status slots.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u32)]
pub enum SurfaceLocalizationStatus {
    /// A valid grid produced strict level crossings.
    Localized = 1,
    /// A valid grid contains no strict level crossing.
    ValidEmpty = 2,
    /// Malformed caller input was refused before any extraction work.
    InvalidInput = 3,
    /// Finite input produced values or geometry without a representable
    /// point result (non-finite samples, coincident level edges, collapsed
    /// coordinates, unrepresentable intersections).
    Unrepresentable = 4,
    /// A checked requirement exceeded its explicit caller envelope.
    ResourceRefused = 5,
    /// The ambient fs-exec cancellation/deadline/poll/cost contract refused.
    Cancelled = 6,
    /// Output or input storage could not be reserved.
    AllocationRefused = 7,
    /// Checked arithmetic inside the kernel failed; never expected from an
    /// admitted plan, always reported rather than inferred.
    InternalFault = 8,
}

impl SurfaceLocalizationStatus {
    /// Stable wire code under [`NEUROSHAPE_LOCALIZATION_SCHEMA_VERSION`].
    #[must_use]
    pub const fn code(self) -> u32 {
        self as u32
    }
}

/// Stable diagnostic naming the exact underlying fs-viz refusal variant.
///
/// Discriminants are the wire codes; they are frozen by
/// [`NEUROSHAPE_LOCALIZATION_SCHEMA_VERSION`] so downstream consumers can
/// branch on them without display-string parsing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u32)]
pub enum LocalizationDiagnostic {
    /// `Grid2Error::InvalidDimensions`.
    GridInvalidDimensions = 1,
    /// `Grid2Error::NodeCountOverflow`.
    GridNodeCountOverflow = 2,
    /// `Grid2Error::NodeBudgetExceeded`.
    GridNodeBudgetExceeded = 3,
    /// `Grid2Error::InvalidBounds`.
    GridInvalidBounds = 4,
    /// `Grid2Error::UnrepresentableCoordinates`.
    GridUnrepresentableCoordinates = 5,
    /// `Grid2Error::NonFiniteValue`.
    GridNonFiniteValue = 6,
    /// `Grid2Error::AllocationFailed`.
    GridAllocationFailed = 7,
    /// `IsoContourError::NonFiniteIso`.
    IsoNonFiniteLevel = 8,
    /// `IsoContourError::ZeroCrossingLimit`.
    IsoZeroCrossingLimit = 9,
    /// `IsoContourError::InvalidPollStride`.
    IsoInvalidPollStride = 10,
    /// `IsoContourError::PlanOverflow`.
    IsoPlanOverflow = 11,
    /// `IsoContourError::OperationBudgetExceeded`.
    IsoOperationBudgetExceeded = 12,
    /// `IsoContourError::ExecutionBudgetRefused`.
    IsoExecutionBudgetRefused = 13,
    /// `IsoContourError::CrossingBudgetExceeded`.
    IsoCrossingBudgetExceeded = 14,
    /// `IsoContourError::CoincidentLevelEdge`.
    IsoCoincidentLevelEdge = 15,
    /// `IsoContourError::UnrepresentableIntersection`.
    IsoUnrepresentableIntersection = 16,
    /// `IsoContourError::AllocationFailed`.
    IsoAllocationFailed = 17,
    /// `IsoContourError::NonFiniteGeometry`.
    IsoNonFiniteGeometry = 18,
}

impl LocalizationDiagnostic {
    /// Stable wire code under [`NEUROSHAPE_LOCALIZATION_SCHEMA_VERSION`].
    #[must_use]
    pub const fn code(self) -> u32 {
        self as u32
    }
}

/// Bounded exact context retained by a refusal outcome.
///
/// Slots that carry no meaning for the active [`LocalizationDiagnostic`]
/// hold the [`LOCALIZATION_DETAIL_UNDEFINED`] sentinels. Edge endpoints are
/// packed losslessly as `(i << 32) | j`; `required`/`limit` saturate at
/// [`u64::MAX`] for `u128` requirements beyond that range. No field ever
/// carries attacker-sized text.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StageDetail {
    /// Stage whose kernel produced the refusal.
    pub stage: LocalizationStage,
    /// Exact underlying refusal variant.
    pub diagnostic: LocalizationDiagnostic,
    /// Offending Cartesian axis or interpolation collapse axis.
    pub axis: u32,
    /// First offender node index, or packed first edge endpoint `[i, j]`.
    pub first_index: u64,
    /// Second offender node index, or packed second edge endpoint `[i, j]`.
    pub second_index: u64,
    /// Exact binary64 bits of the first offending scalar.
    pub scalar_bits: u64,
    /// Exact binary64 bits of the second offending scalar.
    pub second_bits: u64,
    /// Checked requirement, when the refusal compares one against a limit.
    pub required: u64,
    /// Caller-provided limit, when one exists.
    pub limit: u64,
    /// Diagnostic-specific auxiliary code (for example the
    /// `fs_viz::IsoContourResource` ordinal `1..=13` for plan overflow).
    pub aux: u32,
}

impl StageDetail {
    /// A detail whose numeric slots are all undefined except `stage` and
    /// `diagnostic`.
    #[must_use]
    pub const fn new(
        stage: LocalizationStage,
        diagnostic: LocalizationDiagnostic,
    ) -> Self {
        Self {
            stage,
            diagnostic,
            axis: LOCALIZATION_DETAIL_UNDEFINED_U32,
            first_index: LOCALIZATION_DETAIL_UNDEFINED,
            second_index: LOCALIZATION_DETAIL_UNDEFINED,
            scalar_bits: LOCALIZATION_DETAIL_UNDEFINED,
            second_bits: LOCALIZATION_DETAIL_UNDEFINED,
            required: LOCALIZATION_DETAIL_UNDEFINED,
            limit: LOCALIZATION_DETAIL_UNDEFINED,
            aux: LOCALIZATION_DETAIL_UNDEFINED_U32,
        }
    }
}

/// Stable kind of an ambient fs-exec budget/cancellation refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u32)]
pub enum CancellationKind {
    /// `BudgetRefusal::Cancelled`: cancellation observed at a checkpoint.
    CancelledAtCheckpoint = 1,
    /// `BudgetRefusal::DeadlineExpiredAtAdmission`.
    DeadlineExpiredAtAdmission = 2,
    /// `BudgetRefusal::DeadlineWithoutClock`.
    DeadlineWithoutClock = 3,
    /// `BudgetRefusal::CostPlanExceedsQuota`.
    CostPlanExceedsQuota = 4,
    /// `BudgetRefusal::DeadlineExpired`: deadline passed mid-run.
    DeadlineExpiredMidRun = 5,
    /// `BudgetRefusal::PollsExhausted`.
    PollsExhausted = 6,
    /// `BudgetRefusal::CostExhausted`.
    CostQuotaExhausted = 7,
}

/// Ambient execution-contract refusal retained with its exact stage and
/// finalization phase. `phase` is a compile-time `'static` checkpoint name
/// from fs-exec, never attacker-sized text.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CancellationDetail {
    /// Stage whose kernel observed the refusal.
    pub stage: LocalizationStage,
    /// Stable refusal kind.
    pub kind: CancellationKind,
    /// Stable checkpoint phase; empty when refused at admission time.
    pub phase: &'static str,
    /// Deadline nanoseconds for deadline kinds; otherwise undefined.
    pub deadline_ns: u64,
    /// Observed clock nanoseconds for mid-run/admission deadline refusals;
    /// otherwise undefined.
    pub observed_ns: u64,
    /// Requested cost units for cost kinds, admitted poll quota for
    /// exhausted polls; otherwise undefined.
    pub quota_context_a: u64,
    /// Remaining units before an exhausted cost charge; otherwise undefined.
    pub quota_context_b: u64,
    /// Admitted cost quota for cost kinds; otherwise undefined.
    pub quota_context_c: u64,
}

/// Typed outcome of the report's sampled zero-set localization
/// (`fs-viz` Grid2 construction plus isocontour extraction).
///
/// This is the authoritative localization record. The legacy
/// `surface_crossings` / `max_crossing_radius` /
/// `nearest_surface_radius` fields remain only as derived, non-authoritative
/// compatibility views of it: `NaN` sentinels never carry status, and
/// [`SurfaceLocalization::ValidEmpty`] is distinct from every refusal.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SurfaceLocalization {
    /// A valid grid produced `crossings` strict level crossings with the
    /// stated extreme radii.
    Localized {
        /// Number of strict crossings found.
        crossings: usize,
        /// Largest crossing radius.
        max_radius: f64,
        /// Smallest crossing radius.
        nearest_radius: f64,
    },
    /// A valid grid contains no strict level crossing anywhere.
    ValidEmpty,
    /// Malformed caller input refused before any extraction work.
    InvalidInput(StageDetail),
    /// Finite input without a representable point result.
    Unrepresentable(StageDetail),
    /// A checked requirement exceeded its explicit caller envelope.
    ResourceRefused(StageDetail),
    /// Ambient cancellation/deadline/poll/cost refusal.
    Cancelled(CancellationDetail),
    /// Storage could not be reserved.
    AllocationRefused(StageDetail),
    /// Checked kernel arithmetic failed; always reported, never inferred.
    InternalFault(StageDetail),
}

impl SurfaceLocalization {
    /// Coarse stable status class.
    #[must_use]
    pub const fn status(&self) -> SurfaceLocalizationStatus {
        match self {
            Self::Localized { .. } => SurfaceLocalizationStatus::Localized,
            Self::ValidEmpty => SurfaceLocalizationStatus::ValidEmpty,
            Self::InvalidInput(_) => SurfaceLocalizationStatus::InvalidInput,
            Self::Unrepresentable(_) => SurfaceLocalizationStatus::Unrepresentable,
            Self::ResourceRefused(_) => SurfaceLocalizationStatus::ResourceRefused,
            Self::Cancelled(_) => SurfaceLocalizationStatus::Cancelled,
            Self::AllocationRefused(_) => SurfaceLocalizationStatus::AllocationRefused,
            Self::InternalFault(_) => SurfaceLocalizationStatus::InternalFault,
        }
    }

    /// Exact producing stage; `None` for success outcomes.
    #[must_use]
    pub const fn stage(&self) -> Option<LocalizationStage> {
        match self {
            Self::Localized { .. } | Self::ValidEmpty => None,
            Self::InvalidInput(detail)
            | Self::Unrepresentable(detail)
            | Self::ResourceRefused(detail)
            | Self::AllocationRefused(detail)
            | Self::InternalFault(detail) => Some(detail.stage),
            Self::Cancelled(detail) => Some(detail.stage),
        }
    }
}

/// Public campaign input named by a structured admission refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CampaignParameter {
    /// Outer half-width used by the frame and visualization grid.
    RingRadius,
    /// Half-width of the central interval-negative square.
    InnerHalfWidth,
}

impl fmt::Display for CampaignParameter {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::RingRadius => "ring_r",
            Self::InnerHalfWidth => "inner",
        })
    }
}

/// Structured campaign admission refusal for untrusted callers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CampaignError {
    /// NeuroShapeCert is currently a two-dimensional theorem/campaign.
    InputDimension {
        /// Input dimension the campaign requires (always `2` in this tranche).
        expected: usize,
        /// Input dimension the supplied network actually declares.
        actual: usize,
    },
    /// A geometric parameter was NaN or infinite.
    NonFiniteParameter(CampaignParameter),
    /// `ring_r` was non-positive or `inner` was negative.
    OutOfRangeParameter(CampaignParameter),
}

impl fmt::Display for CampaignError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InputDimension { expected, actual } => write!(
                formatter,
                "NeuroShape input dimension mismatch: expected {expected}, got {actual}"
            ),
            Self::NonFiniteParameter(parameter) => {
                write!(formatter, "NeuroShape parameter {parameter} must be finite")
            }
            Self::OutOfRangeParameter(CampaignParameter::RingRadius) => {
                formatter.write_str("NeuroShape parameter ring_r must be positive")
            }
            Self::OutOfRangeParameter(CampaignParameter::InnerHalfWidth) => {
                formatter.write_str("NeuroShape parameter inner must be non-negative")
            }
        }
    }
}

impl std::error::Error for CampaignError {}

/// The blob SDF network. `MlpSdf::new` spectral-normalizes every layer to
/// exactly `bound`, so with `bound = √18` the effective hidden slope is
/// `√18/√2 = 3` (a wall at `|coord| = 0.7`) and the effective output weight is
/// `√18/2 ≈ 2.12`. With the biases below the effective field is
/// `f ≈ 2.12·Σ tanh(3(±coord − 0.7)) + 6.5`: negative near the origin, positive
/// on a surrounding ring. `L = bound² = 18`.
#[must_use]
pub fn blob_sdf_net() -> MlpSdf {
    // Hidden layer: one tanh wall per ±axis direction (bias −2.1 ⇒ wall at 0.7).
    let l1 = Layer::new(
        vec![
            vec![3.0, 0.0],
            vec![-3.0, 0.0],
            vec![0.0, 3.0],
            vec![0.0, -3.0],
        ],
        vec![-2.1, -2.1, -2.1, -2.1],
    );
    // Linear output: sum the walls, lift by +6.5 (bias is not normalized).
    let l2 = Layer::new(vec![vec![1.0, 1.0, 1.0, 1.0]], vec![6.5]);
    MlpSdf::new(vec![l1, l2], (18.0_f64).sqrt())
}

fn is_finite_ordered_interval((lo, hi): (f64, f64)) -> bool {
    lo.is_finite() && hi.is_finite() && lo <= hi
}

/// A constructor-sealed, campaign-local witness that at least one connected
/// component of `{f < 0}` exists and is enclosed by the certified-positive
/// boundary frame.
///
/// `MlpSdf` is a continuous composition of affine maps and `tanh`: the connected
/// negative central square therefore lies in one negative component, and any
/// path from that component to the exterior must cross the positive frame.
///
/// The private fields are important: callers can inspect or clone a witness
/// produced by [`run_campaign`], but cannot manufacture one through safe public
/// constructors from booleans or a sampled contour. This value has no field,
/// source, unit, budget, or receipt identity and therefore is not portable
/// authority. It proves neither that the full negative set is bounded nor that
/// its global component count is exactly one.
#[derive(Debug, Clone, PartialEq)]
pub struct CertifiedEnclosedComponentExists {
    central_box_half_width: f64,
    central_box_interval: (f64, f64),
    boundary_frame_outer_half_width: f64,
    boundary_frame_inner_half_width: f64,
    boundary_strip_intervals: [(f64, f64); 4],
}

impl CertifiedEnclosedComponentExists {
    fn from_interval_frame(
        central_box_half_width: f64,
        central_box_interval: (f64, f64),
        boundary_frame_outer_half_width: f64,
        boundary_frame_width: f64,
        boundary_strip_intervals: [(f64, f64); 4],
    ) -> Option<Self> {
        let boundary_frame_inner_half_width =
            boundary_frame_outer_half_width - boundary_frame_width;
        if !central_box_half_width.is_finite()
            || central_box_half_width < 0.0
            || !boundary_frame_outer_half_width.is_finite()
            || !boundary_frame_width.is_finite()
            || boundary_frame_width <= 0.0
            || !boundary_frame_inner_half_width.is_finite()
            || boundary_frame_inner_half_width <= central_box_half_width
            || !is_finite_ordered_interval(central_box_interval)
            || central_box_interval.1 >= 0.0
            || boundary_strip_intervals
                .iter()
                .any(|&interval| !is_finite_ordered_interval(interval) || interval.0 <= 0.0)
        {
            return None;
        }

        Some(Self {
            central_box_half_width,
            central_box_interval,
            boundary_frame_outer_half_width,
            boundary_frame_inner_half_width,
            boundary_strip_intervals,
        })
    }

    /// Half-width of the central square certified strictly negative.
    #[must_use]
    pub const fn central_box_half_width(&self) -> f64 {
        self.central_box_half_width
    }

    /// Sound IBP enclosure over the central square.
    #[must_use]
    pub const fn central_box_interval(&self) -> (f64, f64) {
        self.central_box_interval
    }

    /// Outer half-width of the square boundary frame.
    #[must_use]
    pub const fn boundary_frame_outer_half_width(&self) -> f64 {
        self.boundary_frame_outer_half_width
    }

    /// Inner half-width of the square boundary frame.
    #[must_use]
    pub const fn boundary_frame_inner_half_width(&self) -> f64 {
        self.boundary_frame_inner_half_width
    }

    /// Sound IBP enclosures for the top, bottom, left, and right frame strips.
    #[must_use]
    pub const fn boundary_strip_intervals(&self) -> &[(f64, f64); 4] {
        &self.boundary_strip_intervals
    }
}

/// What the campaign can state about the global number of negative components.
///
/// This enum is non-exhaustive so a future global topology certificate can add
/// an exact-count state without turning today's lower-bound witness into one.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ComponentCountEvidence {
    /// No positive global component-count statement is certified.
    Unknown,
    /// The closed interval frame certifies that at least one enclosed component
    /// exists. The upper bound remains unknown.
    LowerBound(CertifiedEnclosedComponentExists),
}

impl ComponentCountEvidence {
    /// Certified lower bound on the global component count.
    #[must_use]
    pub const fn lower_bound(&self) -> usize {
        match self {
            Self::Unknown => 0,
            Self::LowerBound(_) => 1,
        }
    }

    /// Certified exact global component count, when available.
    ///
    /// Phase 0 exposes no exact-count certificate, so this is always `None`.
    #[must_use]
    pub const fn exact_count(&self) -> Option<usize> {
        match self {
            Self::Unknown | Self::LowerBound(_) => None,
        }
    }
}

/// Stable ordinal for `fs_viz::IsoContourResource`, frozen under
/// [`NEUROSHAPE_LOCALIZATION_SCHEMA_VERSION`] and used as the auxiliary code
/// for [`LocalizationDiagnostic::IsoPlanOverflow`].
#[must_use]
pub const fn iso_contour_resource_code(resource: fs_viz::IsoContourResource) -> u32 {
    use fs_viz::IsoContourResource as R;
    match resource {
        R::Cells => 1,
        R::EdgeVisits => 2,
        R::ExactOwnershipChecks => 3,
        R::Interpolations => 4,
        R::OutputCrossings => 5,
        R::OutputBytes => 6,
        R::ScratchBytes => 7,
        R::DiagnosticRecords => 8,
        R::DiagnosticBytes => 9,
        R::LiveBytes => 10,
        R::IdentityBytes => 11,
        R::Polls => 12,
        R::WorkUnits => 13,
    }
}

/// Lossless packing of a grid edge endpoint `[i, j]` into one `u64`.
const fn pack_edge(endpoint: [usize; 2]) -> u64 {
    ((endpoint[0] as u64) << 32) | endpoint[1] as u64
}
/// Saturating narrowing used for `u128` requirements/limits.
const fn saturate(value: u128) -> u64 {
    if value > u64::MAX as u128 {
        u64::MAX
    } else {
        value as u64
    }
}

impl From<&BudgetRefusal> for CancellationDetail {
    fn from(refusal: &BudgetRefusal) -> Self {
        let undefined = LOCALIZATION_DETAIL_UNDEFINED;
        match *refusal {
            BudgetRefusal::Cancelled { phase } => Self {
                stage: LocalizationStage::IsoContourExtraction,
                kind: CancellationKind::CancelledAtCheckpoint,
                phase,
                deadline_ns: undefined,
                observed_ns: undefined,
                quota_context_a: undefined,
                quota_context_b: undefined,
                quota_context_c: undefined,
            },
            BudgetRefusal::DeadlineExpiredAtAdmission {
                deadline_ns,
                observed_ns,
            } => Self {
                stage: LocalizationStage::IsoContourExtraction,
                kind: CancellationKind::DeadlineExpiredAtAdmission,
                phase: "",
                deadline_ns,
                observed_ns,
                quota_context_a: undefined,
                quota_context_b: undefined,
                quota_context_c: undefined,
            },
            BudgetRefusal::DeadlineWithoutClock { deadline_ns } => Self {
                stage: LocalizationStage::IsoContourExtraction,
                kind: CancellationKind::DeadlineWithoutClock,
                phase: "",
                deadline_ns,
                observed_ns: undefined,
                quota_context_a: undefined,
                quota_context_b: undefined,
                quota_context_c: undefined,
            },
            BudgetRefusal::CostPlanExceedsQuota { planned, quota } => Self {
                stage: LocalizationStage::IsoContourExtraction,
                kind: CancellationKind::CostPlanExceedsQuota,
                phase: "",
                deadline_ns: undefined,
                observed_ns: undefined,
                quota_context_a: planned,
                quota_context_b: undefined,
                quota_context_c: quota,
            },
            BudgetRefusal::DeadlineExpired {
                phase,
                deadline_ns,
                observed_ns,
            } => Self {
                stage: LocalizationStage::IsoContourExtraction,
                kind: CancellationKind::DeadlineExpiredMidRun,
                phase,
                deadline_ns,
                observed_ns,
                quota_context_a: undefined,
                quota_context_b: undefined,
                quota_context_c: undefined,
            },
            BudgetRefusal::PollsExhausted { phase, quota } => Self {
                stage: LocalizationStage::IsoContourExtraction,
                kind: CancellationKind::PollsExhausted,
                phase,
                deadline_ns: undefined,
                observed_ns: undefined,
                quota_context_a: u64::from(quota),
                quota_context_b: undefined,
                quota_context_c: u64::from(quota),
            },
            BudgetRefusal::CostExhausted {
                phase,
                requested,
                remaining,
                quota,
            } => Self {
                stage: LocalizationStage::IsoContourExtraction,
                kind: CancellationKind::CostQuotaExhausted,
                phase,
                deadline_ns: undefined,
                observed_ns: undefined,
                quota_context_a: requested,
                quota_context_b: remaining,
                quota_context_c: quota,
            },
        }
    }
}

impl From<Grid2Error> for SurfaceLocalization {
    fn from(error: Grid2Error) -> Self {
        let stage = LocalizationStage::GridConstruction;
        match error {
            Grid2Error::InvalidDimensions { dimensions } => {
                let mut detail = StageDetail::new(stage, LocalizationDiagnostic::GridInvalidDimensions);
                detail.first_index = dimensions[0] as u64;
                detail.second_index = dimensions[1] as u64;
                Self::InvalidInput(detail)
            }
            // A dimension product that overflows `usize` is an unadmittable
            // caller-requested scale, not a kernel fault.
            Grid2Error::NodeCountOverflow { dimensions } => {
                let mut detail = StageDetail::new(stage, LocalizationDiagnostic::GridNodeCountOverflow);
                detail.first_index = dimensions[0] as u64;
                detail.second_index = dimensions[1] as u64;
                Self::InvalidInput(detail)
            }
            Grid2Error::NodeBudgetExceeded { required, limit } => {
                let mut detail = StageDetail::new(stage, LocalizationDiagnostic::GridNodeBudgetExceeded);
                detail.required = required as u64;
                detail.limit = limit as u64;
                Self::ResourceRefused(detail)
            }
            Grid2Error::InvalidBounds { axis, lower, upper } => {
                let mut detail = StageDetail::new(stage, LocalizationDiagnostic::GridInvalidBounds);
                detail.axis = axis as u32;
                detail.scalar_bits = lower.to_bits();
                detail.second_bits = upper.to_bits();
                Self::InvalidInput(detail)
            }
            Grid2Error::UnrepresentableCoordinates {
                axis,
                first_index,
                first,
                second_index,
                second,
            } => {
                let mut detail =
                    StageDetail::new(stage, LocalizationDiagnostic::GridUnrepresentableCoordinates);
                detail.axis = axis as u32;
                detail.first_index = first_index as u64;
                detail.second_index = second_index as u64;
                detail.scalar_bits = first.to_bits();
                detail.second_bits = second.to_bits();
                Self::Unrepresentable(detail)
            }
            Grid2Error::NonFiniteValue { index, value } => {
                let mut detail = StageDetail::new(stage, LocalizationDiagnostic::GridNonFiniteValue);
                detail.first_index = index as u64;
                detail.scalar_bits = value.to_bits();
                Self::Unrepresentable(detail)
            }
            Grid2Error::AllocationFailed { nodes } => {
                let mut detail = StageDetail::new(stage, LocalizationDiagnostic::GridAllocationFailed);
                detail.required = nodes as u64;
                Self::AllocationRefused(detail)
            }
        }
    }
}

impl From<IsoContourError> for SurfaceLocalization {
    fn from(error: IsoContourError) -> Self {
        let stage = LocalizationStage::IsoContourExtraction;
        match error {
            IsoContourError::NonFiniteIso { iso } => {
                let mut detail = StageDetail::new(stage, LocalizationDiagnostic::IsoNonFiniteLevel);
                detail.scalar_bits = iso.to_bits();
                Self::InvalidInput(detail)
            }
            IsoContourError::ZeroCrossingLimit => {
                let mut detail = StageDetail::new(stage, LocalizationDiagnostic::IsoZeroCrossingLimit);
                detail.required = 1;
                detail.limit = 0;
                Self::InvalidInput(detail)
            }
            IsoContourError::InvalidPollStride { items_per_poll } => {
                let mut detail = StageDetail::new(stage, LocalizationDiagnostic::IsoInvalidPollStride);
                detail.required = items_per_poll as u64;
                Self::InvalidInput(detail)
            }
            IsoContourError::PlanOverflow { resource } => {
                let mut detail = StageDetail::new(stage, LocalizationDiagnostic::IsoPlanOverflow);
                detail.aux = iso_contour_resource_code(resource);
                Self::InternalFault(detail)
            }
            IsoContourError::OperationBudgetExceeded {
                resource,
                required,
                limit,
            } => {
                let mut detail =
                    StageDetail::new(stage, LocalizationDiagnostic::IsoOperationBudgetExceeded);
                detail.aux = iso_contour_resource_code(resource);
                detail.required = saturate(required);
                detail.limit = saturate(limit);
                Self::ResourceRefused(detail)
            }
            IsoContourError::ExecutionBudgetRefused { refusal } => {
                let mut cancellation = CancellationDetail::from(&refusal);
                cancellation.stage = stage;
                Self::Cancelled(cancellation)
            }
            IsoContourError::CrossingBudgetExceeded { limit } => {
                let mut detail =
                    StageDetail::new(stage, LocalizationDiagnostic::IsoCrossingBudgetExceeded);
                detail.limit = limit as u64;
                Self::ResourceRefused(detail)
            }
            IsoContourError::CoincidentLevelEdge { first, second } => {
                let mut detail = StageDetail::new(stage, LocalizationDiagnostic::IsoCoincidentLevelEdge);
                detail.first_index = pack_edge(first);
                detail.second_index = pack_edge(second);
                Self::Unrepresentable(detail)
            }
            IsoContourError::UnrepresentableIntersection {
                first,
                second,
                iso_bits,
                interpolation_bits,
                collapsed_axis,
                ..
            } => {
                // The full endpoint/value bit dump stays in the fs-viz error at
                // the source stage; this boundary retains offender identity,
                // level bits, interpolation bits, and the collapse axis.
                let mut detail =
                    StageDetail::new(stage, LocalizationDiagnostic::IsoUnrepresentableIntersection);
                detail.first_index = pack_edge(first);
                detail.second_index = pack_edge(second);
                detail.scalar_bits = iso_bits;
                detail.second_bits = interpolation_bits;
                detail.axis = collapsed_axis as u32;
                Self::Unrepresentable(detail)
            }
            IsoContourError::AllocationFailed { required } => {
                let mut detail = StageDetail::new(stage, LocalizationDiagnostic::IsoAllocationFailed);
                detail.required = required as u64;
                Self::AllocationRefused(detail)
            }
            IsoContourError::NonFiniteGeometry => Self::Unrepresentable(StageDetail::new(
                stage,
                LocalizationDiagnostic::IsoNonFiniteGeometry,
            )),
        }
    }
}

/// The campaign report.
#[derive(Debug, Clone)]
pub struct NeuroShapeReport {
    /// Content identity of the normalized MLP and certificate arithmetic.
    pub field_identity: NeuralFieldIdentity,
    /// Governed hidden-activation semantic version.
    pub activation_semantics_version: u32,
    /// Governed hidden-activation semantic name.
    pub activation_semantics: &'static str,
    /// ULP budget used by the interval activation enclosure.
    pub activation_ulp_budget: u64,
    /// Safe-step derivation semantic version.
    pub safe_step_policy_version: u32,
    /// Safe-step derivation semantic name.
    pub safe_step_policy: &'static str,
    /// The certified global Lipschitz constant `L`.
    pub lipschitz: f64,
    /// Nominal field value at the origin, retained for visualization only.
    pub origin_value: f64,
    /// Interval-derived sign margin and downward-rounded no-tunnel step.
    pub safe_step: SafeStepDerivation,
    /// IBP enclosure of `f` over the central box.
    pub inside_interval: (f64, f64),
    /// Is the central box certified strictly inside (`hi < 0`)?
    pub certified_inside: bool,
    /// How many of the box-boundary strips are certified strictly outside.
    pub boundary_certified: usize,
    /// Total boundary strips (4 — a CLOSED frame around the box).
    pub boundary_segments: usize,
    /// Is every strip in the closed boundary frame certified strictly positive?
    /// This is a local frame fact, not a claim that the full negative set is
    /// bounded.
    pub boundary_frame_certified: bool,
    /// Typed component-count evidence. A certified frame yields only a lower
    /// bound of one and never an exact count in this tranche.
    pub component_count_evidence: ComponentCountEvidence,
    /// Is the origin's finite-difference Hessian positive definite under the
    /// classifier tolerance? This is curvature corroboration only: without a
    /// certified zero gradient it does not establish a critical point or local
    /// minimum.
    pub origin_hessian_positive_definite: bool,
    /// Typed outcome of visualization zero-set localization. This is the
    /// authoritative localization record carried by the report and every
    /// native/WASM serialization of it.
    pub surface_localization: SurfaceLocalization,
    /// DERIVED non-authoritative view of [`NeuroShapeReport::surface_localization`]:
    /// crossing count for [`SurfaceLocalizationStatus::Localized`], otherwise
    /// `0`. Zero never distinguishes valid-empty from a refusal.
    pub surface_crossings: usize,
    /// DERIVED view: largest crossing radius when localized, `0` for a valid
    /// empty result, or `NaN` for every refusal.
    pub max_crossing_radius: f64,
    /// DERIVED view: smallest crossing radius (nearest surface point) when
    /// localized, `+inf` for a valid empty result, or `NaN` for every
    /// refusal. The safe step radius must under-estimate it when finite.
    pub nearest_surface_radius: f64,
}

fn radius(p: Vec2) -> f64 {
    p[0].hypot(p[1])
}

/// Run the NeuroShapeCert campaign on `net` with a bounding box of half-width
/// `ring_r` (its four edge strips form the closed barrier) and a central
/// certified-inside box of half-width `inner`.
///
/// # Panics
///
/// Panics on invalid campaign inputs. Use [`try_run_campaign`] at untrusted
/// boundaries that need a structured non-trapping refusal.
#[must_use]
pub fn run_campaign(net: &MlpSdf, ring_r: f64, inner: f64) -> NeuroShapeReport {
    try_run_campaign(net, ring_r, inner).unwrap_or_else(|error| panic!("{error}"))
}

/// Fallible NeuroShape campaign admission and execution.
pub fn try_run_campaign(
    net: &MlpSdf,
    ring_r: f64,
    inner: f64,
) -> Result<NeuroShapeReport, CampaignError> {
    if net.input_dim() != 2 {
        return Err(CampaignError::InputDimension {
            expected: 2,
            actual: net.input_dim(),
        });
    }
    for (parameter, value) in [
        (CampaignParameter::RingRadius, ring_r),
        (CampaignParameter::InnerHalfWidth, inner),
    ] {
        if !value.is_finite() {
            return Err(CampaignError::NonFiniteParameter(parameter));
        }
    }
    if ring_r <= 0.0 {
        return Err(CampaignError::OutOfRangeParameter(
            CampaignParameter::RingRadius,
        ));
    }
    if inner < 0.0 {
        return Err(CampaignError::OutOfRangeParameter(
            CampaignParameter::InnerHalfWidth,
        ));
    }

    let lipschitz = net.lipschitz();
    let origin_value = net.eval(&[0.0, 0.0]);
    let origin_interval = net.eval_interval(&[0.0, 0.0], &[0.0, 0.0]);
    let safe_step = derive_safe_step(origin_interval, lipschitz);

    // Interval topology certificate.
    let inside_interval = net.eval_interval(&[-inner, -inner], &[inner, inner]);
    let certified_inside = is_finite_ordered_interval(inside_interval) && inside_interval.1 < 0.0;

    // A CLOSED barrier: the four edge strips of the box [−R, R]² tile the whole
    // boundary frame (corners overlap), so certifying every strip strictly
    // outside (lo > 0) RIGOROUSLY traps the negative component meeting the
    // central box. It does not exclude other interior or exterior components.
    // Eight discrete boxes would leave angular gaps and prove no enclosure.
    let r = ring_r;
    let w = 0.4;
    let strips = [
        ([-r, r - w], [r, r]),   // top
        ([-r, -r], [r, -r + w]), // bottom
        ([-r, -r], [-r + w, r]), // left
        ([r - w, -r], [r, r]),   // right
    ];
    let boundary_segments = strips.len();
    let boundary_strip_intervals = strips.map(|(lo_pt, hi_pt)| net.eval_interval(&lo_pt, &hi_pt));
    let boundary_certified = boundary_strip_intervals
        .iter()
        .filter(|&&interval| is_finite_ordered_interval(interval) && interval.0 > 0.0)
        .count();
    let boundary_frame_certified = boundary_certified == boundary_segments;
    let component_count_evidence = CertifiedEnclosedComponentExists::from_interval_frame(
        inner,
        inside_interval,
        ring_r,
        w,
        boundary_strip_intervals,
    )
    .map_or(
        ComponentCountEvidence::Unknown,
        ComponentCountEvidence::LowerBound,
    );

    // Curvature cross-check at the origin (Hessian by finite difference). This
    // does not establish criticality because the gradient is not certified zero.
    let h = 1e-3;
    let f00 = origin_value;
    let fxx = (net.eval(&[h, 0.0]) - 2.0 * f00 + net.eval(&[-h, 0.0])) / (h * h);
    let fyy = (net.eval(&[0.0, h]) - 2.0 * f00 + net.eval(&[0.0, -h])) / (h * h);
    let fxy = (net.eval(&[h, h]) - net.eval(&[h, -h]) - net.eval(&[-h, h]) + net.eval(&[-h, -h]))
        / (4.0 * h * h);
    let crit = classify_hessian([[fxx, fxy], [fxy, fyy]], 1e-6);
    let origin_hessian_positive_definite = crit.kind == CriticalKind::Minimum;

    // Localize the zero set on a visualization grid. The typed outcome below
    // is authoritative; the legacy crossing/radius fields are derived views
    // of it (bead frankensim-o33vo).
    const GRID_N: usize = 81;
    const CROSSING_LIMIT: usize = 2 * GRID_N * (GRID_N - 1);
    let surface_localization = match Grid2::from_fn(
        GRID_N,
        GRID_N,
        [-ring_r - 0.5, -ring_r - 0.5],
        [ring_r + 0.5, ring_r + 0.5],
        GRID_N * GRID_N,
        |p| net.eval(&[p[0], p[1]]),
    ) {
        Ok(grid) => match grid.isocontour_crossings(0.0, CROSSING_LIMIT) {
            Ok(crossings) if crossings.is_empty() => SurfaceLocalization::ValidEmpty,
            Ok(crossings) => {
                let mut max_radius = 0.0_f64;
                let mut nearest_radius = f64::INFINITY;
                for point in &crossings {
                    let r = radius(*point);
                    max_radius = max_radius.max(r);
                    nearest_radius = nearest_radius.min(r);
                }
                SurfaceLocalization::Localized {
                    crossings: crossings.len(),
                    max_radius,
                    nearest_radius,
                }
            }
            Err(error) => SurfaceLocalization::from(error),
        },
        Err(error) => SurfaceLocalization::from(error),
    };
    let (surface_crossings, max_crossing_radius, nearest_surface_radius) =
        match &surface_localization {
            SurfaceLocalization::Localized {
                crossings,
                max_radius,
                nearest_radius,
            } => (*crossings, *max_radius, *nearest_radius),
            // Valid empty keeps its historical `0` / `+inf` sentinels; every
            // refusal reports `NaN` radii so NaN never carries status alone.
            SurfaceLocalization::ValidEmpty => (0, 0.0, f64::INFINITY),
            _ => (0, f64::NAN, f64::NAN),
        };

    Ok(NeuroShapeReport {
        field_identity: net.identity(),
        activation_semantics_version: MLP_ACTIVATION_SEMANTICS_VERSION,
        activation_semantics: MLP_ACTIVATION_SEMANTICS,
        activation_ulp_budget: MLP_ACTIVATION_ULP_BUDGET,
        safe_step_policy_version: SAFE_STEP_POLICY_VERSION,
        safe_step_policy: SAFE_STEP_POLICY,
        lipschitz,
        origin_value,
        safe_step,
        inside_interval,
        certified_inside,
        boundary_certified,
        boundary_segments,
        boundary_frame_certified,
        component_count_evidence,
        origin_hessian_positive_definite,
        surface_crossings,
        surface_localization,
        max_crossing_radius,
        nearest_surface_radius,
    })
}
