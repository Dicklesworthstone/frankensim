//! fs-oed-e2e — SensorForge: optimal experimental design that knows when to
//! stop. Layer: L4 (ASCENT).
//!
//! # The campaign
//!
//! You must pick the best of several designs, but their performances are only
//! estimated; you can spend sensors to sharpen them. Which do you measure, and
//! when have you measured enough? This answers both with evidence, composing
//! crates never designed to meet:
//!
//! - **Kalman fusion** ([`fs_assimilate`]): each candidate is a Gaussian belief;
//!   a sensor reading is fused with the exact scalar Kalman update, shrinking that
//!   candidate's posterior variance.
//! - **Value of information** ([`fs_voi`]): at each step the Expected Value of
//!   Perfect Information scores the decision's ambiguity; the campaign's
//!   cancellation-aware action-value reduction places the next sensor on the
//!   candidate whose measurement most sharpens the DECISION (not the
//!   most-uncertain candidate), and says STOP the instant EVPI falls below
//!   threshold — the design choice is already robust.
//! - **Budget allocation** ([`fs_toleralloc`]): the measurement-precision budget
//!   is then distributed cost-optimally across candidates by sensitivity.
//! - **Honest colors** ([`fs_evidence`]): posterior variance and EVPI remain
//!   `Estimated`; their bounded identities commit to every campaign input and
//!   every instrument-bound assimilation candidate.
//!
//! Deterministic (sensor readings hit each candidate's true value; the Kalman
//! variance update is observation-independent). Public inference values carry
//! `fs-qty` dimensions and optional semantic kinds; lower scalar kernels see
//! coherent-SI `f64` only after admission.

use core::fmt;
use std::collections::{BTreeMap, BTreeSet};

use fs_assimilate::{AssimError, Belief, assimilate_colored_with_shared_poll_quota, point_sensor};
use fs_evidence::{Color, ColorRank, color_leaf_identity_reason};
use fs_exec::Cx;
use fs_qty::semantic::{QuantityKind, SemanticQty, SemanticType, ValueForm};
use fs_qty::{Dims, QUANTITY_SPEC_ENCODED_LEN, QtyAny, QuantitySpec};
use fs_toleralloc::{Feature, allocate};
use fs_voi::{
    Action, ActionKind, ActionValue, DesignEstimate, Recommendation, Uncertainty,
    expected_opportunity_loss_by,
};

/// Maximum accepted candidate-name length.
pub const MAX_CANDIDATE_NAME_BYTES: usize = 128;
/// Maximum number of candidates in one synchronous campaign.
pub const MAX_CAMPAIGN_CANDIDATES: usize = 256;
/// Maximum number of sensor placements in one synchronous campaign.
pub const MAX_CAMPAIGN_SENSORS: usize = 4_096;
/// Maximum admitted action-design work units. One sensor-action score evaluates
/// the decision model at every retained normal-quadrature point.
pub const MAX_CAMPAIGN_EVALUATIONS: usize = 10_500_000;

/// Semantic version of the sealed SensorForge report estimator identities.
// v8 (bead sj31i.7): the identity preimage additionally binds the exact
// objective dimension/semantic schema. Numerically identical campaigns under
// incompatible quantity kinds therefore cannot share evidence identities.
// v7 (bead sj31i.62): the identity preimage additionally binds the
// admitted byte plan and the byte-accounting policy version, so two
// campaigns that agree on science but were admitted under different
// byte envelopes carry distinct identities.
// v6 (bead sj31i.62): the campaign canonicalizes candidate order at
// admission, so the identity preimage binds the CANONICAL declaration
// sequence — permuted caller menus now collapse to one identity where
// v5 deliberately kept declaration order identity-semantic.
pub const OED_REPORT_IDENTITY_VERSION: u64 = 8;

const REPORT_ID_DOMAIN: &str = "org.frankensim.fs-oed-e2e.report.v6";
// v4 (bead sj31i.5): robustness-bearing EVPI evaluations and action
// ranking use the same full multi-alternative opportunity loss (one work
// unit still means one record's participation in one bounded evaluation;
// the quadrature depth is a fixed constant).
const CAMPAIGN_PLANNING_POLICY_VERSION: u64 = 4;
const CAMPAIGN_POLL_POLICY_VERSION: u64 = 2;
const CAMPAIGN_RECORD_POLL_STRIDE: usize = 256;
const CAMPAIGN_ACTION_POLL_STRIDE: usize = 1;
// Byte-accounting policy (bead sj31i.62): every resource-driving
// operation is charged in BYTE UNITS against a worst-case plan
// preflighted at admission with checked arithmetic. A byte unit is a
// deterministic accounting bound on bytes visited, compared, hashed,
// or retained at one bounded seam — charges are formula-based upper
// bounds evaluated on the actual shape, never measured allocator
// traffic, so ledgers replay bit-identically across hosts.
// Byte policy v2 (bead sj31i.5): full-EOL evaluations charge one menu
// sweep per quadrature node (FULL_EOL_SCAN_SWEEPS); the pre-action
// baseline scan is charged at the action-value seam.
// v3 binds and charges the 12-byte objective schema at admission and in each
// report preimage, and carries that schema token into instrument identities.
const CAMPAIGN_BYTE_POLICY_VERSION: u64 = 3;
// One retained candidate/estimate/posterior record's fixed scalar
// payload (means, variances, uncertainty components), name excluded.
const RECORD_SCALAR_BYTE_UNITS: u128 = 32;
// One EVPI/summary scan reads a mean and a total deviation per record.
const SCAN_READ_BYTE_UNITS: u128 = 16;
// One FULL expected-opportunity-loss evaluation (bead sj31i.5) sweeps
// the menu once per quadrature node plus the best-scan and window
// passes.
const FULL_EOL_SCAN_SWEEPS: u128 = fs_voi::EOL_QUADRATURE_PANELS as u128 + 3;
// `measure-` / `sensor-` identity prefixes added to candidate names.
const ACTION_PREFIX_BYTE_UNITS: u128 = 8;
const OBJECTIVE_SCHEMA_BYTES: u128 = QUANTITY_SPEC_ENCODED_LEN as u128;
const OBJECTIVE_SCHEMA_HEX_BYTES: u128 = OBJECTIVE_SCHEMA_BYTES * 2;
const SENSOR_IDENTITY_FIXED_BYTES: u128 = 7 + 2 + OBJECTIVE_SCHEMA_HEX_BYTES;

// Nine-point Gauss-Hermite rule, transformed and normalized for expectations
// under N(0, 1). The policy is deterministic and substantially more faithful
// than evaluating only at the unchanged posterior mean. It remains an
// Estimated decision model: no quadrature-remainder certificate is claimed.
const NORMAL_EXPECTATION_RULE: [(f64, f64); 9] = [
    (-4.512_745_863_399_783_5, 2.234_584_400_774_658_3e-5),
    (-3.205_429_002_856_470_3, 0.002_789_141_321_231_769),
    (-2.076_847_978_677_83, 0.049_916_406_765_217_88),
    (-1.023_255_663_789_132_6, 0.244_097_502_894_939_45),
    (0.0, 0.406_349_206_349_206_35),
    (1.023_255_663_789_132_6, 0.244_097_502_894_939_45),
    (2.076_847_978_677_83, 0.049_916_406_765_217_88),
    (3.205_429_002_856_470_3, 0.002_789_141_321_231_769),
    (4.512_745_863_399_783_5, 2.234_584_400_774_658_3e-5),
];
const ACTION_EVALUATION_FACTOR: usize = NORMAL_EXPECTATION_RULE.len() + 2;

fn canonicalize_zero(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

/// Exact quantity schema for a scalar OED objective.
///
/// A dimension-only schema is an explicit no-kind claim, not a wildcard.
/// Semantic schemas retain the exact quantity kind and scalar form, so
/// dimensionally equal pressure/stress or absolute/delta-temperature values
/// remain incompatible at campaign admission.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ObjectiveSpec {
    dims: Dims,
    semantic_type: Option<SemanticType>,
}

impl ObjectiveSpec {
    /// Construct an explicit dimension-only objective schema.
    #[must_use]
    pub const fn dimensional(dims: Dims) -> Self {
        Self {
            dims,
            semantic_type: None,
        }
    }

    /// Construct the exact schema retained by an already-validated semantic
    /// quantity. Keep this private: `SemanticType` is also a diagnostic
    /// descriptor, while public construction must pass through
    /// `SemanticQty`'s kind/form/value validation.
    const fn from_validated_semantic(semantic_type: SemanticType) -> Self {
        Self {
            dims: semantic_type.expected_dims(),
            semantic_type: Some(semantic_type),
        }
    }

    /// Six-base coherent-SI dimension vector.
    #[must_use]
    pub const fn dims(self) -> Dims {
        self.dims
    }

    /// Exact semantic kind/form, or `None` for an explicit dimension-only
    /// declaration.
    #[must_use]
    pub const fn semantic_type(self) -> Option<SemanticType> {
        self.semantic_type
    }

    /// Reusable fs-qty schema retained by this OED-specific wrapper.
    #[must_use]
    pub const fn quantity_spec(self) -> QuantitySpec {
        match self.semantic_type {
            Some(semantic_type) => QuantitySpec::semantic(semantic_type),
            None => QuantitySpec::dimensional(self.dims),
        }
    }

    /// Dimension vector required by prior, sensor-noise, and posterior
    /// variances (`Q²`).
    pub fn variance_dims(self) -> Result<Dims, CandidateError> {
        QtyAny::new(1.0, self.dims)
            .powi(2)
            .map(|quantity| quantity.dims)
            .map_err(|_| CandidateError::DimensionOverflow {
                field: "objective variance",
                dims: self.dims,
                exponent: 2,
            })
    }

    /// Schema of a decision difference such as EVPI or its stop threshold.
    ///
    /// Absolute and delta temperatures map to delta temperature. Other
    /// semantic kinds intentionally become dimension-only: `fs-qty` does not
    /// currently define a general decision-loss/difference kind, and tagging
    /// an expected loss as pressure RMS, a composition fraction, or another
    /// measured-value form would be false semantics. The original objective
    /// schema remains separately retained and identity-bound.
    #[must_use]
    pub fn decision_spec(self) -> Self {
        let semantic_type = self.semantic_type.and_then(|semantic_type| {
            matches!(
                semantic_type.kind(),
                QuantityKind::AbsoluteTemperature | QuantityKind::TemperatureDifference
            )
            .then(|| SemanticType::new(QuantityKind::TemperatureDifference, ValueForm::Static))
        });
        Self {
            dims: self.dims,
            semantic_type,
        }
    }

    /// Admit a finite EVPI/threshold value under this objective's derived
    /// decision-difference schema.
    pub fn decision_value(self, value: f64) -> Result<ObjectiveValue, CandidateError> {
        if !value.is_finite() {
            return Err(CandidateError::InvalidNumber {
                field: "decision value",
                requirement: "finite",
            });
        }
        Ok(ObjectiveValue::from_raw(value, self.decision_spec()))
    }

    fn canonical_bytes(self) -> [u8; QUANTITY_SPEC_ENCODED_LEN] {
        self.quantity_spec().canonical_bytes()
    }

    fn identity_hex(self) -> String {
        const HEX: &[u8; 16] = b"0123456789abcdef";
        let bytes = self.canonical_bytes();
        let mut encoded = String::with_capacity(bytes.len() * 2);
        for byte in bytes {
            encoded.push(HEX[usize::from(byte >> 4)] as char);
            encoded.push(HEX[usize::from(byte & 0x0f)] as char);
        }
        encoded
    }
}

/// A finite coherent-SI objective scalar carrying its exact inference schema.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ObjectiveValue {
    quantity: QtyAny,
    spec: ObjectiveSpec,
}

impl ObjectiveValue {
    /// Admit a finite scalar under an explicit dimension-only schema.
    pub fn dimensional(quantity: QtyAny) -> Result<Self, CandidateError> {
        if !quantity.value.is_finite() {
            return Err(CandidateError::InvalidNumber {
                field: "objective value",
                requirement: "finite",
            });
        }
        Ok(Self {
            quantity: QtyAny::new(canonicalize_zero(quantity.value), quantity.dims),
            spec: ObjectiveSpec::dimensional(quantity.dims),
        })
    }

    /// Admit an already-validated semantic quantity without erasing its kind
    /// or scalar form.
    #[must_use]
    pub fn semantic(quantity: SemanticQty) -> Self {
        Self {
            quantity: QtyAny::new(
                canonicalize_zero(quantity.value()),
                quantity.quantity().dims,
            ),
            spec: ObjectiveSpec::from_validated_semantic(quantity.semantic_type()),
        }
    }

    /// Admit a finite dimensionless scalar under an explicit no-kind schema.
    pub fn dimensionless(value: f64) -> Result<Self, CandidateError> {
        Self::dimensional(QtyAny::dimensionless(value))
    }

    /// Coherent-SI value and dimension vector.
    #[must_use]
    pub const fn quantity(self) -> QtyAny {
        self.quantity
    }

    /// Raw coherent-SI scalar. The schema remains available through
    /// [`ObjectiveValue::spec`].
    #[must_use]
    pub const fn value(self) -> f64 {
        self.quantity.value
    }

    /// Raw coherent-SI bit pattern for deterministic evidence comparisons.
    #[must_use]
    pub const fn to_bits(self) -> u64 {
        self.quantity.value.to_bits()
    }

    /// Exact dimension/semantic schema.
    #[must_use]
    pub const fn spec(self) -> ObjectiveSpec {
        self.spec
    }

    fn from_raw(value: f64, spec: ObjectiveSpec) -> Self {
        debug_assert!(value.is_finite());
        debug_assert!(spec.semantic_type.is_none_or(|semantic_type| {
            SemanticQty::new(QtyAny::new(value, spec.dims), semantic_type).is_ok()
        }));
        Self {
            quantity: QtyAny::new(canonicalize_zero(value), spec.dims),
            spec,
        }
    }
}

/// A rejected candidate declaration.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CandidateError {
    /// The candidate name cannot serve as a bounded provenance identity.
    InvalidName {
        /// Structural rejection reason.
        reason: &'static str,
    },
    /// A numeric field violates its declared domain.
    InvalidNumber {
        /// Offending field.
        field: &'static str,
        /// Required domain.
        requirement: &'static str,
    },
    /// Two objective scalars carry incompatible dimension/semantic schemas.
    ObjectiveSchemaMismatch {
        /// Offending field.
        field: &'static str,
        /// Schema supplied by the caller.
        actual: ObjectiveSpec,
        /// Schema required by the declaration.
        expected: ObjectiveSpec,
    },
    /// A variance/noise/cost field carries incompatible dimensions.
    DimensionMismatch {
        /// Offending field.
        field: &'static str,
        /// Dimension vector supplied by the caller.
        actual: Dims,
        /// Dimension vector required by the declaration.
        expected: Dims,
    },
    /// Deriving a required power of the objective dimension overflowed the
    /// admitted exponent domain.
    DimensionOverflow {
        /// Derived field.
        field: &'static str,
        /// Source dimensions.
        dims: Dims,
        /// Requested integer exponent.
        exponent: i8,
    },
}

impl fmt::Display for CandidateError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidName { reason } => {
                write!(f, "candidate name is not an admissible identity: {reason}")
            }
            Self::InvalidNumber { field, requirement } => {
                write!(f, "candidate `{field}` must be {requirement}")
            }
            Self::ObjectiveSchemaMismatch {
                field,
                actual,
                expected,
            } => write!(
                f,
                "candidate `{field}` objective schema {actual:?} does not match {expected:?}"
            ),
            Self::DimensionMismatch {
                field,
                actual,
                expected,
            } => write!(
                f,
                "candidate `{field}` dimensions {actual:?} do not match required {expected:?}"
            ),
            Self::DimensionOverflow {
                field,
                dims,
                exponent,
            } => write!(
                f,
                "candidate `{field}` dimensions {dims:?}^{exponent} overflow the supported exponent domain"
            ),
        }
    }
}

impl std::error::Error for CandidateError {}

/// A candidate design under measurement.
#[derive(Debug, Clone, PartialEq)]
pub struct Candidate {
    name: String,
    objective_spec: ObjectiveSpec,
    truth: f64,
    prior_mean: f64,
    prior_var: f64,
    sensor_noise_variance: f64,
    sensor_cost: f64,
}

impl Candidate {
    /// Construct a checked candidate.
    ///
    /// # Errors
    /// Returns [`CandidateError`] for an unusable name, incompatible objective
    /// schemas, non-`Q²` variance/noise, a non-dimensionless cost, or an invalid
    /// numeric domain.
    pub fn new(
        name: impl Into<String>,
        truth: ObjectiveValue,
        prior_mean: ObjectiveValue,
        prior_variance: QtyAny,
        sensor_noise_variance: QtyAny,
        sensor_cost: QtyAny,
    ) -> Result<Self, CandidateError> {
        let name = name.into();
        let name_reason = if name.len() > MAX_CANDIDATE_NAME_BYTES {
            Some("too-long")
        } else {
            color_leaf_identity_reason(&name)
        };
        if let Some(reason) = name_reason {
            return Err(CandidateError::InvalidName { reason });
        }
        let objective_spec = truth.spec;
        if prior_mean.spec != objective_spec {
            return Err(CandidateError::ObjectiveSchemaMismatch {
                field: "prior_mean",
                actual: prior_mean.spec,
                expected: objective_spec,
            });
        }
        let variance_dims = objective_spec.variance_dims()?;
        for (field, quantity) in [
            ("prior_variance", prior_variance),
            ("sensor_noise_variance", sensor_noise_variance),
        ] {
            if quantity.dims != variance_dims {
                return Err(CandidateError::DimensionMismatch {
                    field,
                    actual: quantity.dims,
                    expected: variance_dims,
                });
            }
        }
        if !prior_variance.value.is_finite() || prior_variance.value < 0.0 {
            return Err(CandidateError::InvalidNumber {
                field: "prior_variance",
                requirement: "finite and non-negative",
            });
        }
        if !sensor_noise_variance.value.is_finite() || sensor_noise_variance.value <= 0.0 {
            return Err(CandidateError::InvalidNumber {
                field: "sensor_noise_variance",
                requirement: "finite and positive",
            });
        }
        if sensor_cost.dims != Dims::NONE {
            return Err(CandidateError::DimensionMismatch {
                field: "sensor_cost",
                actual: sensor_cost.dims,
                expected: Dims::NONE,
            });
        }
        if !sensor_cost.value.is_finite() || sensor_cost.value <= 0.0 {
            return Err(CandidateError::InvalidNumber {
                field: "sensor_cost",
                requirement: "finite, positive, and dimensionless",
            });
        }
        Ok(Self {
            name,
            objective_spec,
            truth: truth.quantity.value,
            prior_mean: prior_mean.quantity.value,
            prior_var: canonicalize_zero(prior_variance.value),
            sensor_noise_variance: canonicalize_zero(sensor_noise_variance.value),
            sensor_cost: canonicalize_zero(sensor_cost.value),
        })
    }

    /// Candidate identity.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Sensor reading used by this deterministic worked campaign.
    #[must_use]
    pub fn truth(&self) -> ObjectiveValue {
        ObjectiveValue::from_raw(self.truth, self.objective_spec)
    }

    /// Prior objective mean.
    #[must_use]
    pub fn prior_mean(&self) -> ObjectiveValue {
        ObjectiveValue::from_raw(self.prior_mean, self.objective_spec)
    }

    /// Prior objective variance.
    #[must_use]
    pub fn prior_variance(&self) -> QtyAny {
        QtyAny::new(
            self.prior_var,
            self.objective_spec
                .variance_dims()
                .expect("candidate construction sealed objective variance dimensions"),
        )
    }

    /// Sensor noise variance.
    #[must_use]
    pub fn sensor_noise_variance(&self) -> QtyAny {
        QtyAny::new(
            self.sensor_noise_variance,
            self.objective_spec
                .variance_dims()
                .expect("candidate construction sealed objective variance dimensions"),
        )
    }

    /// Cost of one measurement.
    #[must_use]
    pub const fn sensor_cost(&self) -> QtyAny {
        QtyAny::dimensionless(self.sensor_cost)
    }

    /// Exact objective schema carried by means, variances, EVPI, and evidence.
    #[must_use]
    pub const fn objective_spec(&self) -> ObjectiveSpec {
        self.objective_spec
    }
}

/// A rejected campaign or failed campaign computation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OedError {
    /// The ambient admitted budget refused admission or further work
    /// (deadline or cost); the typed refusal is retained verbatim and
    /// no partial report is published (bead sj31i.6).
    BudgetRefused(fs_exec::BudgetRefusal),
    /// At least one candidate is required.
    NoCandidates,
    /// The synchronous candidate cap was exceeded.
    TooManyCandidates {
        /// Requested count.
        count: usize,
        /// Accepted maximum.
        max: usize,
    },
    /// The synchronous placement cap was exceeded.
    TooManySensors {
        /// Requested count.
        count: usize,
        /// Accepted maximum.
        max: usize,
    },
    /// The requested planning work exceeds the synchronous campaign budget.
    WorkBudgetExceeded {
        /// Candidate count.
        candidates: usize,
        /// Requested placement cap.
        max_sensors: usize,
        /// Requested action-design evaluations.
        evaluations: usize,
        /// Accepted maximum product.
        max_evaluations: usize,
    },
    /// Cancellation or poll-quota exhaustion was observed at a deterministic
    /// campaign boundary.
    Cancelled {
        /// Phase whose boundary observed the request.
        phase: &'static str,
        /// Placements committed before the request was observed.
        completed_placements: usize,
        /// Logical work units completed before the request was observed.
        completed_work_units: u128,
        /// Admitted worst-case work bound for the requested campaign cap.
        admitted_work_units: u128,
    },
    /// A lower-layer assimilation observed cancellation after this campaign
    /// had completed the selected sensor observation but before committing a
    /// posterior or placement.
    AssimilationCancelled {
        /// Candidate whose posterior update was cancelled.
        candidate: String,
        /// Placements committed before the lower-layer request was observed.
        completed_placements: usize,
        /// Campaign work completed before entering the cancelled update.
        completed_work_units: u128,
        /// Admitted worst-case campaign work bound.
        admitted_work_units: u128,
        /// Structured lower-layer phase and progress evidence.
        source: Box<AssimError>,
    },
    /// Executed logical work did not match the exact realized shape or exceeded
    /// the admitted worst-case bound.
    WorkPlanMismatch {
        /// Work credited by the execution ledger.
        completed_work_units: u128,
        /// Exact work implied by the realized early-stop path.
        realized_work_units: u128,
        /// Worst-case work admitted before scientific execution.
        admitted_work_units: u128,
    },
    /// The worst-case byte plan overflowed checked admission arithmetic.
    BytePlanOverflow {
        /// Candidate count.
        candidates: usize,
        /// Requested placement cap.
        max_sensors: usize,
    },
    /// A byte charge at a bounded seam would exceed the admitted
    /// worst-case byte plan. This is a fail-closed accounting defect
    /// surface: admission upper-bounds every charging seam, so an
    /// exceeded plan means a seam charges more than admission declared.
    ByteBudgetExceeded {
        /// Seam that attempted the charge.
        at: &'static str,
        /// Ledger total the charge would have reached.
        charged_byte_units: u128,
        /// Worst-case byte units admitted before scientific execution.
        admitted_byte_units: u128,
    },
    /// A retained-output allocation was refused by the allocator.
    OutputAllocationRefused {
        /// Output whose reservation failed.
        what: &'static str,
    },
    /// The EVPI stop threshold must be finite and non-negative.
    InvalidThreshold,
    /// Candidate objective schemas must agree before any scientific work.
    ObjectiveSchemaMismatch {
        /// Candidate carrying the incompatible schema.
        candidate: String,
        /// Schema supplied by that candidate.
        actual: ObjectiveSpec,
        /// Schema established by the canonical campaign declaration.
        expected: ObjectiveSpec,
    },
    /// The stop threshold must carry the objective's decision-difference
    /// schema (absolute temperature objectives require delta-temperature
    /// thresholds).
    ThresholdSchemaMismatch {
        /// Threshold schema supplied by the caller.
        actual: ObjectiveSpec,
        /// Required decision-difference schema.
        expected: ObjectiveSpec,
    },
    /// Candidate identities must be unique because actions address them by name.
    DuplicateCandidate {
        /// Repeated identity.
        name: String,
    },
    /// A checked scalar belief unexpectedly rejected an internal access.
    BeliefInvariant(AssimError),
    /// An observation or posterior update failed.
    Assimilation {
        /// Candidate being measured.
        candidate: String,
        /// Structured lower-layer failure.
        source: AssimError,
    },
    /// The bounded VoI reduction returned an action outside its own menu.
    UnknownRecommendation {
        /// Returned action identity.
        action: String,
    },
    /// A deterministic derived quantity overflowed or became NaN.
    NonFiniteComputation {
        /// Quantity whose contract failed.
        quantity: &'static str,
    },
    /// The tolerance allocator rejected checked positive-sensitivity inputs.
    AllocationFailed,
    /// The allocator omitted a checked positive-sensitivity candidate.
    MissingAllocation {
        /// Missing candidate identity.
        candidate: String,
    },
    /// A design menu presented to the canonical-order constructor was
    /// not in strict canonical (name-ascending, duplicate-free) order.
    CanonicalOrderViolated {
        /// First position whose entry breaks the order.
        position: usize,
    },
    /// A mean-override view was constructed with an out-of-range index
    /// or a non-finite override payload.
    OverrideInvalid {
        /// What is wrong.
        what: &'static str,
    },
}

impl fmt::Display for OedError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::BudgetRefused(refusal) => {
                write!(f, "campaign budget refused: {refusal}")
            }
            Self::NoCandidates => write!(f, "SensorForge needs at least one candidate"),
            Self::TooManyCandidates { count, max } => {
                write!(f, "candidate count {count} exceeds synchronous cap {max}")
            }
            Self::TooManySensors { count, max } => {
                write!(f, "sensor cap {count} exceeds synchronous cap {max}")
            }
            Self::WorkBudgetExceeded {
                candidates,
                max_sensors,
                evaluations,
                max_evaluations,
            } => write!(
                f,
                "campaign work {candidates}^2 x {max_sensors} x \
                 {ACTION_EVALUATION_FACTOR} = {evaluations} exceeds \
                {max_evaluations} action-design evaluations"
            ),
            Self::Cancelled {
                phase,
                completed_placements,
                completed_work_units,
                admitted_work_units,
            } => write!(
                f,
                "campaign cancelled or poll budget exhausted during {phase} after \
                 {completed_placements} placements and \
                {completed_work_units}/{admitted_work_units} admitted logical work units"
            ),
            Self::AssimilationCancelled {
                candidate,
                completed_placements,
                completed_work_units,
                admitted_work_units,
                source,
            } => write!(
                f,
                "assimilation for candidate `{candidate}` cancelled after \
                 {completed_placements} committed placements and \
                 {completed_work_units}/{admitted_work_units} admitted campaign work units: {source}"
            ),
            Self::WorkPlanMismatch {
                completed_work_units,
                realized_work_units,
                admitted_work_units,
            } => write!(
                f,
                "campaign work ledger mismatch: completed {completed_work_units}, \
                 realized {realized_work_units}, admitted {admitted_work_units}"
            ),
            Self::BytePlanOverflow {
                candidates,
                max_sensors,
            } => write!(
                f,
                "worst-case byte plan for {candidates} candidate(s) and \
                 {max_sensors} placement(s) overflowed checked admission arithmetic"
            ),
            Self::ByteBudgetExceeded {
                at,
                charged_byte_units,
                admitted_byte_units,
            } => write!(
                f,
                "byte charge at `{at}` would reach {charged_byte_units} of \
                 {admitted_byte_units} admitted byte units — a seam charges more \
                 than admission declared"
            ),
            Self::OutputAllocationRefused { what } => {
                write!(f, "retained-output allocation for {what} was refused")
            }
            Self::InvalidThreshold => {
                write!(f, "EVPI threshold must be finite and non-negative")
            }
            Self::ObjectiveSchemaMismatch {
                candidate,
                actual,
                expected,
            } => write!(
                f,
                "candidate `{candidate}` objective schema {actual:?} does not match campaign schema {expected:?}"
            ),
            Self::ThresholdSchemaMismatch { actual, expected } => write!(
                f,
                "EVPI threshold schema {actual:?} does not match decision schema {expected:?}"
            ),
            Self::DuplicateCandidate { name } => {
                write!(f, "candidate identity `{name}` is duplicated")
            }
            Self::BeliefInvariant(source) => write!(f, "scalar belief invariant failed: {source}"),
            Self::Assimilation { candidate, source } => {
                write!(
                    f,
                    "assimilation failed for candidate `{candidate}`: {source}"
                )
            }
            Self::UnknownRecommendation { action } => {
                write!(f, "VoI reduction returned unknown action `{action}`")
            }
            Self::NonFiniteComputation { quantity } => {
                write!(f, "campaign produced non-finite `{quantity}`")
            }
            Self::AllocationFailed => write!(f, "precision allocation failed"),
            Self::MissingAllocation { candidate } => {
                write!(f, "precision allocation omitted candidate `{candidate}`")
            }
            Self::CanonicalOrderViolated { position } => write!(
                f,
                "design menu is not in canonical name order at position {position}; \
                 canonicalize once at campaign admission"
            ),
            Self::OverrideInvalid { what } => {
                write!(f, "mean-override view rejected: {what}")
            }
        }
    }
}

impl std::error::Error for OedError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::BeliefInvariant(source) | Self::Assimilation { source, .. } => Some(source),
            Self::BudgetRefused(refusal) => Some(refusal),
            Self::AssimilationCancelled { source, .. } => Some(source.as_ref()),
            _ => None,
        }
    }
}

/// Final scalar posterior for one candidate.
#[derive(Debug, Clone, PartialEq)]
pub struct PosteriorSummary {
    /// Candidate identity.
    name: String,
    /// Exact schema of the posterior mean.
    objective_spec: ObjectiveSpec,
    /// Posterior mean.
    mean: f64,
    /// Posterior variance.
    variance: f64,
}

impl PosteriorSummary {
    /// Candidate identity.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Posterior mean with its exact objective schema.
    #[must_use]
    pub fn mean(&self) -> ObjectiveValue {
        ObjectiveValue::from_raw(self.mean, self.objective_spec)
    }

    /// Posterior variance with dimensions `Q²`.
    #[must_use]
    pub fn variance(&self) -> QtyAny {
        QtyAny::new(
            self.variance,
            self.objective_spec
                .variance_dims()
                .expect("posterior construction sealed variance dimensions"),
        )
    }
}

/// The campaign report.
///
/// Reports are sealed outputs of [`run_campaign`]. Read-only accessors expose
/// the complete result without permitting callers to replace science fields or
/// evidence identities independently.
#[derive(Debug, Clone, PartialEq)]
pub struct OedReport {
    /// Exact schema of candidate truth/prior/posterior means.
    objective_spec: ObjectiveSpec,
    /// Candidate names in the order sensors were placed.
    placements: Vec<String>,
    /// Number of sensors placed.
    sensors_placed: usize,
    /// Total prior variance across candidates.
    prior_total_variance: f64,
    /// Total posterior variance across candidates.
    posterior_total_variance: f64,
    /// Fractional variance reduction.
    variance_reduction: f64,
    /// EVPI before any sensor.
    initial_evpi: f64,
    /// EVPI after the campaign stopped.
    final_evpi: f64,
    /// Did the decision become robust (planner chose to STOP)?
    decision_robust: bool,
    /// The finally-chosen (lowest-cost posterior) design.
    chosen_design: String,
    /// The cost-optimal tolerance allocation `(name, tolerance)`.
    /// A zero-sensitivity candidate receives `+infinity`, the exact unconstrained
    /// optimum under the first-order allocation model.
    allocation: Vec<(String, f64)>,
    /// EVPI before sensing and after every completed placement.
    evpi_trace: Vec<f64>,
    /// Final scalar posterior in candidate order.
    posteriors: Vec<PosteriorSummary>,
    /// Instrument-bound estimated candidate emitted by each assimilation.
    assimilation_colors: Vec<Color>,
    /// The posterior-variance color (`Estimated` until independently certified).
    variance_color: Color,
    /// The EVPI color (`Estimated` — decision-theoretic).
    evpi_color: Color,
    /// Worst-case byte units admitted before scientific work started.
    admitted_byte_units: u128,
    /// Deterministic byte units charged across every bounded seam. This
    /// upper-bounds peak transient traffic; it is an accounting bound,
    /// not measured allocator bytes.
    consumed_byte_units: u128,
    /// The charged subset that remains live in this published report.
    retained_byte_units: u128,
}

impl OedReport {
    /// Exact schema retained by objective means and derived decision values.
    #[must_use]
    pub const fn objective_spec(&self) -> ObjectiveSpec {
        self.objective_spec
    }

    /// Schema of EVPI and its stop threshold. Absolute-temperature objectives
    /// expose delta-temperature decision values.
    #[must_use]
    pub fn decision_spec(&self) -> ObjectiveSpec {
        self.objective_spec.decision_spec()
    }

    /// Candidate names in placement order.
    #[must_use]
    pub fn placements(&self) -> &[String] {
        &self.placements
    }

    /// Number of completed sensor placements.
    #[must_use]
    pub const fn sensors_placed(&self) -> usize {
        self.sensors_placed
    }

    /// Total variance before sensing.
    #[must_use]
    pub fn prior_total_variance(&self) -> QtyAny {
        QtyAny::new(
            self.prior_total_variance,
            self.objective_spec
                .variance_dims()
                .expect("report construction sealed variance dimensions"),
        )
    }

    /// Total variance after sensing.
    #[must_use]
    pub fn posterior_total_variance(&self) -> QtyAny {
        QtyAny::new(
            self.posterior_total_variance,
            self.objective_spec
                .variance_dims()
                .expect("report construction sealed variance dimensions"),
        )
    }

    /// Fractional reduction in total variance.
    #[must_use]
    pub const fn variance_reduction(&self) -> f64 {
        self.variance_reduction
    }

    /// EVPI before the first placement.
    #[must_use]
    pub fn initial_evpi(&self) -> ObjectiveValue {
        ObjectiveValue::from_raw(self.initial_evpi, self.objective_spec.decision_spec())
    }

    /// EVPI when the campaign stopped.
    #[must_use]
    pub fn final_evpi(&self) -> ObjectiveValue {
        ObjectiveValue::from_raw(self.final_evpi, self.objective_spec.decision_spec())
    }

    /// Whether the modeled EVPI met the requested stop threshold.
    #[must_use]
    pub const fn decision_robust(&self) -> bool {
        self.decision_robust
    }

    /// Finally chosen design identity.
    #[must_use]
    pub fn chosen_design(&self) -> &str {
        &self.chosen_design
    }

    /// Cost-optimal tolerance allocation in candidate order.
    #[must_use]
    pub fn allocation(&self) -> &[(String, f64)] {
        &self.allocation
    }

    /// EVPI before sensing and after each completed placement.
    #[must_use]
    pub fn evpi_trace(&self) -> impl ExactSizeIterator<Item = ObjectiveValue> + '_ {
        let spec = self.objective_spec.decision_spec();
        self.evpi_trace
            .iter()
            .copied()
            .map(move |value| ObjectiveValue::from_raw(value, spec))
    }

    /// Final scalar posterior summaries in candidate order.
    #[must_use]
    pub fn posteriors(&self) -> &[PosteriorSummary] {
        &self.posteriors
    }

    /// Instrument-bound colors emitted by completed assimilations.
    #[must_use]
    pub fn assimilation_colors(&self) -> &[Color] {
        &self.assimilation_colors
    }

    /// Sealed posterior-variance evidence color.
    #[must_use]
    pub const fn variance_color(&self) -> &Color {
        &self.variance_color
    }

    /// Sealed EVPI evidence color.
    #[must_use]
    pub const fn evpi_color(&self) -> &Color {
        &self.evpi_color
    }

    /// Worst-case byte units admitted before scientific work started.
    #[must_use]
    pub const fn admitted_byte_units(&self) -> u128 {
        self.admitted_byte_units
    }

    /// Deterministic byte units charged across every bounded seam
    /// (an upper bound on peak transient traffic, not measured
    /// allocator bytes).
    #[must_use]
    pub const fn consumed_byte_units(&self) -> u128 {
        self.consumed_byte_units
    }

    /// The charged subset that remains live in this published report.
    #[must_use]
    pub const fn retained_byte_units(&self) -> u128 {
        self.retained_byte_units
    }
}

/// The preflighted worst-case logical work bound. A unit is one bounded
/// candidate/color record visit, one scalar assimilation transaction, or one
/// retained hash; it is a deterministic scheduling/accounting unit, not an
/// instruction count. Early STOP paths realize fewer units and are checked
/// separately rather than padded with phantom work.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CampaignWorkPlan {
    candidates: usize,
    max_sensors: usize,
    action_design_evaluations: usize,
    setup_work_units: u128,
    per_placement_work_units: u128,
    maximum_finalization_work_units: u128,
    admitted_work_units: u128,
}

/// Exact logical work implied by one realized early-stop path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CampaignRealizedWorkPlan {
    completed_placements: u128,
    action_rounds: u128,
    positive_prior_candidates: u128,
    setup_work_units: u128,
    placement_work_units: u128,
    incomplete_action_work_units: u128,
    finalization_work_units: u128,
    realized_work_units: u128,
}

impl CampaignRealizedWorkPlan {
    const fn identity_fields(self) -> [u128; 8] {
        [
            self.completed_placements,
            self.action_rounds,
            self.positive_prior_candidates,
            self.setup_work_units,
            self.placement_work_units,
            self.incomplete_action_work_units,
            self.finalization_work_units,
            self.realized_work_units,
        ]
    }
}

impl CampaignWorkPlan {
    fn checked(candidates: usize, max_sensors: usize) -> Result<Self, OedError> {
        let action_design_pairs =
            candidates
                .checked_mul(candidates)
                .ok_or(OedError::WorkBudgetExceeded {
                    candidates,
                    max_sensors,
                    evaluations: usize::MAX,
                    max_evaluations: MAX_CAMPAIGN_EVALUATIONS,
                })?;
        let action_design_evaluations = action_design_pairs
            .checked_mul(max_sensors)
            .and_then(|work| work.checked_mul(ACTION_EVALUATION_FACTOR))
            .ok_or(OedError::WorkBudgetExceeded {
                candidates,
                max_sensors,
                evaluations: usize::MAX,
                max_evaluations: MAX_CAMPAIGN_EVALUATIONS,
            })?;
        if action_design_evaluations > MAX_CAMPAIGN_EVALUATIONS {
            return Err(OedError::WorkBudgetExceeded {
                candidates,
                max_sensors,
                evaluations: action_design_evaluations,
                max_evaluations: MAX_CAMPAIGN_EVALUATIONS,
            });
        }

        let n = candidates as u128;
        let placements = max_sensors as u128;
        // Setup: validation, belief construction, prior variance, and initial
        // estimates/EVPI (five candidate scans). Sensor actions are rebuilt
        // after every posterior update because their effect depends on P and R.
        let setup_work_units = n.checked_mul(5).ok_or(OedError::WorkBudgetExceeded {
            candidates,
            max_sensors,
            evaluations: usize::MAX,
            max_evaluations: MAX_CAMPAIGN_EVALUATIONS,
        })?;
        // Each placement: one action construction scan; for every action, a
        // target lookup, posterior-template construction, and one EVPI scan per
        // normal-expectation node; then the action-menu record, chosen-action
        // lookup, sensor+assimilation transaction, and refreshed estimates/EVPI.
        let per_placement_work_units = n
            .checked_mul(n)
            .and_then(|work| work.checked_mul(ACTION_EVALUATION_FACTOR as u128))
            .and_then(|work| work.checked_add(n.checked_mul(5)?))
            .and_then(|work| work.checked_add(2))
            .ok_or(OedError::WorkBudgetExceeded {
                candidates,
                max_sensors,
                evaluations: usize::MAX,
                max_evaluations: MAX_CAMPAIGN_EVALUATIONS,
            })?;
        // Worst-case final science consumes twelve candidate scans. Each full
        // report identity reserves three candidate-sized and three
        // max-placement-sized sequences, the trace's initial value, and one
        // bounded hash. The realized plan later substitutes the actual positive
        // priors and placements. The last unit is publication.
        let maximum_finalization_work_units = n
            .checked_mul(18)
            .and_then(|work| work.checked_add(placements.checked_mul(6)?))
            .and_then(|work| work.checked_add(5))
            .ok_or(OedError::WorkBudgetExceeded {
                candidates,
                max_sensors,
                evaluations: usize::MAX,
                max_evaluations: MAX_CAMPAIGN_EVALUATIONS,
            })?;
        let admitted_work_units = per_placement_work_units
            .checked_mul(placements)
            .and_then(|work| work.checked_add(setup_work_units))
            .and_then(|work| work.checked_add(maximum_finalization_work_units))
            .ok_or(OedError::WorkBudgetExceeded {
                candidates,
                max_sensors,
                evaluations: usize::MAX,
                max_evaluations: MAX_CAMPAIGN_EVALUATIONS,
            })?;

        Ok(Self {
            candidates,
            max_sensors,
            action_design_evaluations,
            setup_work_units,
            per_placement_work_units,
            maximum_finalization_work_units,
            admitted_work_units,
        })
    }

    fn realized(
        self,
        completed_placements: usize,
        action_rounds: usize,
        positive_prior_candidates: usize,
    ) -> Option<CampaignRealizedWorkPlan> {
        if completed_placements > self.max_sensors
            || action_rounds < completed_placements
            || action_rounds.checked_sub(completed_placements)? > 1
            || positive_prior_candidates > self.candidates
        {
            return None;
        }
        let n = self.candidates as u128;
        let completed_placements = completed_placements as u128;
        let action_rounds = action_rounds as u128;
        let positive_prior_candidates = positive_prior_candidates as u128;
        let placement_work_units = self
            .per_placement_work_units
            .checked_mul(completed_placements)?;
        let incomplete_rounds = action_rounds.checked_sub(completed_placements)?;
        let incomplete_action_work_units = n
            .checked_mul(n)?
            .checked_mul(ACTION_EVALUATION_FACTOR as u128)?
            .checked_add(n.checked_mul(2)?)?
            .checked_mul(incomplete_rounds)?;
        // Final science is 7n + 5m, where m is the number of positive-prior
        // candidates accepted by fs-toleralloc. Two report identities add
        // 6n + 6s + 4, and publication adds one.
        let finalization_work_units = n
            .checked_mul(13)?
            .checked_add(positive_prior_candidates.checked_mul(5)?)?
            .checked_add(completed_placements.checked_mul(6)?)?
            .checked_add(5)?;
        let realized_work_units = self
            .setup_work_units
            .checked_add(placement_work_units)?
            .checked_add(incomplete_action_work_units)?
            .checked_add(finalization_work_units)?;
        if realized_work_units > self.admitted_work_units {
            return None;
        }
        Some(CampaignRealizedWorkPlan {
            completed_placements,
            action_rounds,
            positive_prior_candidates,
            setup_work_units: self.setup_work_units,
            placement_work_units,
            incomplete_action_work_units,
            finalization_work_units,
            realized_work_units,
        })
    }
}

/// The preflighted worst-case byte bound (bead sj31i.62). A byte unit
/// is a deterministic upper bound on bytes visited, compared, hashed,
/// or retained at one bounded seam. Seams charge formula-based bounds
/// evaluated on the ACTUAL shape; admission evaluates the same
/// formulas at the worst-case shape with checked arithmetic, so every
/// seam charge is covered before scientific work starts and the
/// in-flight ledger can exceed the plan only through an accounting
/// defect — which refuses typed instead of publishing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CampaignBytePlan {
    admission_byte_units: u128,
    setup_byte_units: u128,
    per_placement_byte_units: u128,
    maximum_finalization_byte_units: u128,
    admitted_byte_units: u128,
}

/// One canonical-menu rebuild: fresh estimate records (scalar payloads
/// plus cloned names) and the O(n) strict-order window verification
/// (two name reads per adjacent pair).
fn menu_build_byte_units(n: u128, sum_name: u128, max_name: u128) -> Option<u128> {
    n.checked_mul(RECORD_SCALAR_BYTE_UNITS)?
        .checked_add(sum_name)?
        .checked_add(n.saturating_sub(1).checked_mul(max_name.checked_mul(2)?)?)
}

/// The actual campaign shape's byte coefficients, computed once after
/// admission so every seam charges the same deterministic formulas the
/// admitted plan bounded. Saturation cannot launder an overflow: the
/// checked plan already refused any shape whose formulas overflow, and
/// a saturated charge exceeding the plan refuses typed at the seam.
#[derive(Debug, Clone, Copy)]
struct ByteShape {
    n: u128,
    sum_name: u128,
    max_name: u128,
}

impl ByteShape {
    fn of(candidates: &[Candidate]) -> Self {
        Self {
            n: candidates.len() as u128,
            sum_name: candidates
                .iter()
                .map(|candidate| candidate.name.len() as u128)
                .sum(),
            max_name: candidates
                .iter()
                .map(|candidate| candidate.name.len() as u128)
                .max()
                .unwrap_or(0),
        }
    }

    /// One record's fixed scalar payload per candidate.
    fn records(self) -> u128 {
        self.n.saturating_mul(RECORD_SCALAR_BYTE_UNITS)
    }

    /// One scalar read per candidate.
    fn scalar_reads(self) -> u128 {
        self.n.saturating_mul(8)
    }

    /// One EVPI/summary scan over the whole menu.
    fn scan(self) -> u128 {
        self.n.saturating_mul(SCAN_READ_BYTE_UNITS)
    }

    /// One FULL expected-opportunity-loss evaluation over the menu.
    fn full_scan(self) -> u128 {
        self.scan().saturating_mul(FULL_EOL_SCAN_SWEEPS)
    }

    /// One canonical-menu rebuild including window verification.
    fn menu_build(self) -> u128 {
        menu_build_byte_units(self.n, self.sum_name, self.max_name).unwrap_or(u128::MAX)
    }

    /// One sensor-action construction sweep (records plus prefixed
    /// action names).
    fn action_construction(self) -> u128 {
        self.n
            .saturating_mul(RECORD_SCALAR_BYTE_UNITS.saturating_add(ACTION_PREFIX_BYTE_UNITS))
            .saturating_add(self.sum_name)
    }

    /// One action's quadrature-view evaluation: the view construction
    /// reads three scalars, then every expectation node runs one FULL
    /// opportunity-loss evaluation over the menu (bead sj31i.5).
    fn action_evaluation(self) -> u128 {
        (ACTION_EVALUATION_FACTOR as u128)
            .saturating_mul(self.full_scan())
            .saturating_add(24)
    }

    /// The chosen-action linear lookup bound.
    fn chosen_lookup(self) -> u128 {
        self.n
            .saturating_mul(self.max_name.saturating_add(ACTION_PREFIX_BYTE_UNITS))
    }
}

impl CampaignBytePlan {
    /// Never-undercounting comparison count for the one-time admission
    /// sort: `n * bit_width(n) + n` dominates merge sort's worst case.
    fn sort_comparison_bound(n: u128) -> Option<u128> {
        let width = u128::from(u128::BITS - n.max(1).leading_zeros());
        n.checked_mul(width)?.checked_add(n)
    }

    /// Worst-case one-identity preimage length: a fixed envelope for
    /// version/mode/stream/budget/plan/policy/rule fields plus
    /// length-prefixed candidate, placement, allocation, trace,
    /// posterior, and color sequences at their maximum shapes. Every
    /// term dominates the exact builder's contribution.
    fn identity_preimage_bound(
        n: u128,
        placements: u128,
        sum_name: u128,
        max_name: u128,
    ) -> Option<u128> {
        let candidate_rows = n.checked_mul(48)?.checked_add(sum_name)?;
        let placement_rows = placements.checked_mul(max_name.checked_add(8)?)?;
        let allocation_rows = n.checked_mul(16)?.checked_add(sum_name)?;
        let trace_rows = placements.checked_add(2)?.checked_mul(8)?;
        let posterior_rows = n.checked_mul(24)?.checked_add(sum_name)?;
        let color_rows = placements.checked_mul(256)?;
        1024u128
            .checked_add(OBJECTIVE_SCHEMA_BYTES)?
            .checked_add(candidate_rows)?
            .checked_add(placement_rows)?
            .checked_add(allocation_rows)?
            .checked_add(trace_rows)?
            .checked_add(posterior_rows)?
            .checked_add(color_rows)?
            .checked_add(max_name.checked_add(8)?)
    }

    fn checked(
        candidates: usize,
        max_sensors: usize,
        sum_name: usize,
        max_name: usize,
    ) -> Result<Self, OedError> {
        let overflow = OedError::BytePlanOverflow {
            candidates,
            max_sensors,
        };
        let n = candidates as u128;
        let placements = max_sensors as u128;
        let sum_name = sum_name as u128;
        let max_name = max_name as u128;
        let evaluation_factor = ACTION_EVALUATION_FACTOR as u128;

        // Admission: the duplicate-identity scan reads every name once;
        // the one-time canonical sort compares two names per comparison;
        // and every exact schema comparison reads both 12-byte operands.
        let admission_byte_units = (|| {
            let schema_comparisons = n
                .checked_add(1)?
                .checked_mul(2)?
                .checked_mul(OBJECTIVE_SCHEMA_BYTES)?;
            Self::sort_comparison_bound(n)?
                .checked_mul(max_name.checked_mul(2)?)?
                .checked_add(sum_name)?
                .checked_add(schema_comparisons)
        })()
        .ok_or_else(|| overflow.clone())?;

        // Setup mirrors the five admitted setup scans: belief records,
        // prior-variance reads, the initial menu build, and the initial
        // EVPI scan.
        let setup_byte_units = (|| {
            n.checked_mul(RECORD_SCALAR_BYTE_UNITS)?
                .checked_add(n.checked_mul(8)?)?
                .checked_add(menu_build_byte_units(n, sum_name, max_name)?)?
                .checked_add(
                    n.checked_mul(SCAN_READ_BYTE_UNITS)?
                        .checked_mul(FULL_EOL_SCAN_SWEEPS)?,
                )
        })()
        .ok_or_else(|| overflow.clone())?;

        // One full placement round: action construction (records plus
        // prefixed action names), the quadrature-view evaluation of
        // every action (each evaluation scans the whole menu; view
        // construction reads three scalars), the chosen-action lookup,
        // the observation record, the committed posterior/color/
        // placement/trace retention, and the menu + EVPI refreshes.
        let per_placement_byte_units = (|| {
            let action_construction = n
                .checked_mul(RECORD_SCALAR_BYTE_UNITS.checked_add(ACTION_PREFIX_BYTE_UNITS)?)?
                .checked_add(sum_name)?;
            let action_evaluation = n.checked_mul(
                evaluation_factor
                    .checked_mul(
                        n.checked_mul(SCAN_READ_BYTE_UNITS)?
                            .checked_mul(FULL_EOL_SCAN_SWEEPS)?,
                    )?
                    .checked_add(24)?,
            )?;
            let chosen_lookup = n.checked_mul(max_name.checked_add(ACTION_PREFIX_BYTE_UNITS)?)?;
            let observation = SENSOR_IDENTITY_FIXED_BYTES
                .checked_add(max_name)?
                .checked_add(RECORD_SCALAR_BYTE_UNITS)?;
            let commit_retention = RECORD_SCALAR_BYTE_UNITS
                .checked_add(32)?
                .checked_add(max_name)?
                .checked_add(8)?;
            action_construction
                .checked_add(action_evaluation)?
                .checked_add(chosen_lookup)?
                .checked_add(observation)?
                .checked_add(commit_retention)?
                .checked_add(menu_build_byte_units(n, sum_name, max_name)?)?
                .checked_add(
                    n.checked_mul(SCAN_READ_BYTE_UNITS)?
                        .checked_mul(FULL_EOL_SCAN_SWEEPS)?,
                )
        })()
        .ok_or_else(|| overflow.clone())?;

        // Finalization: final menu/EVPI/variance/chosen-design scans,
        // the retained allocation and posterior summaries, and two
        // report identities (preimage hash bytes plus the retained
        // identity strings, bounded at 128 each).
        let maximum_finalization_byte_units = (|| {
            let scans = menu_build_byte_units(n, sum_name, max_name)?
                .checked_add(
                    n.checked_mul(SCAN_READ_BYTE_UNITS)?
                        .checked_mul(FULL_EOL_SCAN_SWEEPS)?,
                )?
                .checked_add(n.checked_mul(8)?)?
                .checked_add(n.checked_mul(max_name.checked_add(8)?)?)?
                .checked_add(max_name)?;
            let allocation = n
                .checked_mul(SCAN_READ_BYTE_UNITS)?
                .checked_add(n.checked_mul(max_name.checked_add(8)?)?)?;
            let posteriors = n.checked_mul(
                max_name
                    .checked_add(16)?
                    .checked_add(OBJECTIVE_SCHEMA_BYTES)?,
            )?;
            let identities = Self::identity_preimage_bound(n, placements, sum_name, max_name)?
                .checked_add(64)?
                .checked_add(128)?
                .checked_mul(2)?;
            scans
                .checked_add(allocation)?
                .checked_add(posteriors)?
                .checked_add(OBJECTIVE_SCHEMA_BYTES)?
                .checked_add(identities)
        })()
        .ok_or_else(|| overflow.clone())?;

        let admitted_byte_units = per_placement_byte_units
            .checked_mul(placements)
            .and_then(|total| total.checked_add(admission_byte_units))
            .and_then(|total| total.checked_add(setup_byte_units))
            .and_then(|total| total.checked_add(maximum_finalization_byte_units))
            .ok_or(overflow)?;

        Ok(Self {
            admission_byte_units,
            setup_byte_units,
            per_placement_byte_units,
            maximum_finalization_byte_units,
            admitted_byte_units,
        })
    }

    const fn identity_fields(self) -> [u128; 5] {
        [
            self.admission_byte_units,
            self.setup_byte_units,
            self.per_placement_byte_units,
            self.maximum_finalization_byte_units,
            self.admitted_byte_units,
        ]
    }
}

#[derive(Debug)]
struct CampaignProgress<'s> {
    completed_placements: usize,
    completed_work_units: u128,
    // Private invocation-global ledger: only `checkpoint` and the nested
    // assimilation transaction can decrease this value, and no caller can
    // replace it between campaign phases.
    polls_remaining: u32,
    // Byte ledger (bead sj31i.62): deterministic seam charges against
    // the admitted worst-case byte plan, with the retained subset
    // tracked separately. Charges refuse typed instead of exceeding
    // the plan.
    byte_plan: CampaignBytePlan,
    charged_byte_units: u128,
    retained_byte_units: u128,
    // Ambient deadline/cost authority (bead sj31i.6); polls stay in the
    // raw ledger above because nested assimilation shares it.
    ambient: fs_exec::AdmittedBudget<'s>,
}

impl<'s> CampaignProgress<'s> {
    fn admit(
        cx: &Cx<'s>,
        completed_work_units: u128,
        planned_cost: u64,
        byte_plan: CampaignBytePlan,
    ) -> Result<Self, OedError> {
        let ambient = fs_exec::AdmittedBudget::admit_ambient(cx, planned_cost)
            .map_err(OedError::BudgetRefused)?;
        let mut progress = Self {
            completed_placements: 0,
            completed_work_units,
            polls_remaining: cx.budget().poll_quota,
            byte_plan,
            charged_byte_units: 0,
            retained_byte_units: 0,
            ambient,
        };
        // The admission scan and one-time canonical sort already ran;
        // their declared bound is the ledger's opening charge.
        progress.charge_bytes("campaign admission", byte_plan.admission_byte_units)?;
        Ok(progress)
    }

    /// Charge one bounded seam's deterministic byte bound. Exceeding
    /// the admitted plan is an accounting defect and refuses typed —
    /// no partial report can be published past a refused charge.
    fn charge_bytes(&mut self, at: &'static str, bytes: u128) -> Result<(), OedError> {
        let charged =
            self.charged_byte_units
                .checked_add(bytes)
                .ok_or(OedError::ByteBudgetExceeded {
                    at,
                    charged_byte_units: u128::MAX,
                    admitted_byte_units: self.byte_plan.admitted_byte_units,
                })?;
        if charged > self.byte_plan.admitted_byte_units {
            return Err(OedError::ByteBudgetExceeded {
                at,
                charged_byte_units: charged,
                admitted_byte_units: self.byte_plan.admitted_byte_units,
            });
        }
        self.charged_byte_units = charged;
        Ok(())
    }

    /// Charge bytes that stay live in the published report.
    fn retain_bytes(&mut self, at: &'static str, bytes: u128) -> Result<(), OedError> {
        self.charge_bytes(at, bytes)?;
        self.retained_byte_units = self
            .retained_byte_units
            .checked_add(bytes)
            .expect("retained bytes are a subset of the checked charge ledger");
        Ok(())
    }

    fn advance(&mut self, units: u128) {
        self.completed_work_units = self
            .completed_work_units
            .checked_add(units)
            .expect("admitted campaign progress cannot exceed u128");
        // Cost accrues with admitted work; exhaustion beyond the admitted
        // plan is impossible because admission bounded the plan, so the
        // charge only feeds retained consumption evidence.
        let _ = self
            .ambient
            .charge_cost("campaign-work", u64::try_from(units).unwrap_or(u64::MAX));
    }

    fn checkpoint(
        &mut self,
        cx: &Cx<'_>,
        plan: CampaignWorkPlan,
        phase: &'static str,
    ) -> Result<(), OedError> {
        // Deadline and cancellation first (sj31i.6); the poll ledger below
        // stays a raw counter because nested assimilation shares it.
        self.ambient
            .observe_deadline(phase, cx)
            .map_err(|refusal| match refusal {
                fs_exec::BudgetRefusal::Cancelled { .. } => self.cancelled(plan, phase),
                other => OedError::BudgetRefused(other),
            })?;
        if self.polls_remaining == 0 {
            return Err(self.cancelled(plan, phase));
        }
        if self.polls_remaining != u32::MAX {
            self.polls_remaining -= 1;
        }
        cx.checkpoint().map_err(|_| self.cancelled(plan, phase))
    }

    fn finish(
        &self,
        plan: CampaignWorkPlan,
        realized: CampaignRealizedWorkPlan,
    ) -> Result<(), OedError> {
        if self.completed_work_units == realized.realized_work_units
            && realized.realized_work_units <= plan.admitted_work_units
        {
            Ok(())
        } else {
            Err(OedError::WorkPlanMismatch {
                completed_work_units: self.completed_work_units,
                realized_work_units: realized.realized_work_units,
                admitted_work_units: plan.admitted_work_units,
            })
        }
    }

    fn cancelled(&self, plan: CampaignWorkPlan, phase: &'static str) -> OedError {
        OedError::Cancelled {
            phase,
            completed_placements: self.completed_placements,
            completed_work_units: self.completed_work_units,
            admitted_work_units: plan.admitted_work_units,
        }
    }
}

fn to_estimates(
    candidates: &[Candidate],
    beliefs: &[Belief],
) -> Result<Vec<DesignEstimate>, OedError> {
    if candidates.len() != beliefs.len() {
        return Err(OedError::NonFiniteComputation {
            quantity: "candidate/belief cardinality",
        });
    }
    candidates
        .iter()
        .zip(beliefs)
        .map(|(c, b)| {
            let mean = b.component_mean(0).map_err(OedError::BeliefInvariant)?;
            let variance = b.variance(0).map_err(OedError::BeliefInvariant)?;
            Ok(DesignEstimate::new(
                c.name.clone(),
                mean,
                Uncertainty {
                    numerical: 0.0,
                    statistical: variance.sqrt(),
                    model: 0.0,
                },
            ))
        })
        .collect()
}

fn total_variance(beliefs: &[Belief]) -> Result<f64, OedError> {
    beliefs.iter().try_fold(0.0, |total, belief| {
        let variance = belief.variance(0).map_err(OedError::BeliefInvariant)?;
        let next = total + variance;
        if next.is_finite() {
            Ok(next)
        } else {
            Err(OedError::NonFiniteComputation {
                quantity: "total variance",
            })
        }
    })
}

/// The campaign's design menu in CANONICAL identity order (bead
/// sj31i.62). Identity and order are validated ONCE at construction —
/// strictly ascending unique names — and are immutable thereafter:
/// values refresh only through [`CanonicalDesignMenu::from_canonical`]
/// on estimates that are already in canonical order (the campaign
/// canonicalizes its candidates once at admission, so every derived
/// estimate vector inherits the order for free). Evaluation over the
/// menu neither clones nor sorts; canonical order supplies the
/// equal-mean tie-break the old clone-and-sort imposed per call.
/// Robustness-bearing values (initial/final EVPI, the STOP gate, the
/// trace) and action ranking come from the same FULL multi-alternative
/// [`fs_voi::expected_opportunity_loss_by`] algebra (bead sj31i.5).
#[derive(Debug)]
struct CanonicalDesignMenu {
    estimates: Vec<DesignEstimate>,
}

impl CanonicalDesignMenu {
    /// Wrap estimates that are ALREADY in canonical order, verifying
    /// the representation invariant in one O(n) window scan (no sort,
    /// no allocation).
    fn from_canonical(estimates: Vec<DesignEstimate>) -> Result<Self, OedError> {
        if let Some(position) = estimates
            .windows(2)
            .position(|pair| pair[0].name >= pair[1].name)
        {
            return Err(OedError::CanonicalOrderViolated {
                position: position + 1,
            });
        }
        Ok(Self { estimates })
    }

    fn estimates(&self) -> &[DesignEstimate] {
        &self.estimates
    }

    fn len(&self) -> usize {
        self.estimates.len()
    }

    /// O(log n) identity lookup — canonical order makes the old linear
    /// scan unnecessary.
    fn index_of(&self, name: &str) -> Option<usize> {
        self.estimates
            .binary_search_by(|estimate| estimate.name.as_str().cmp(name))
            .ok()
    }

    /// Allocation-free, sort-free FULL multi-alternative expected
    /// opportunity loss (bead sj31i.5) with the campaign's finiteness
    /// and sign contract — the only value allowed to feed robustness
    /// decisions, traces, and report science.
    fn full_opportunity_loss_checked(&self) -> Result<f64, OedError> {
        let value = expected_opportunity_loss_by(
            self.estimates.len(),
            &|idx| self.estimates[idx].mean,
            &|idx| self.estimates[idx].uncertainty.total_std(),
        );
        if value.is_finite() && value >= 0.0 {
            Ok(canonicalize_zero(value))
        } else {
            Err(OedError::NonFiniteComputation { quantity: "EVPI" })
        }
    }
}

/// A typed, NON-OWNING one-index substitution over an immutable
/// [`CanonicalDesignMenu`]: the predictive-quadrature EVPI evaluation
/// sees the target's overridden mean and posterior statistical
/// uncertainty without cloning, mutating, or restoring shared scratch
/// — stale restoration and cancellation-corrupted menu state are
/// unrepresentable because nothing is ever written. The borrow ties
/// the view's lifetime to the menu, so it cannot escape the call.
struct MeanOverrideView<'menu> {
    menu: &'menu CanonicalDesignMenu,
    index: usize,
    mean: f64,
    total_std: f64,
}

impl<'menu> MeanOverrideView<'menu> {
    /// Validate the selected index and finite override payload, and
    /// precompute the target's overridden total uncertainty (its
    /// numerical and model components are read from the immutable
    /// menu; only the statistical component is substituted).
    fn new(
        menu: &'menu CanonicalDesignMenu,
        index: usize,
        mean: f64,
        statistical_std: f64,
    ) -> Result<Self, OedError> {
        let Some(target) = menu.estimates.get(index) else {
            return Err(OedError::OverrideInvalid {
                what: "override index is outside the canonical menu",
            });
        };
        if !mean.is_finite() {
            return Err(OedError::OverrideInvalid {
                what: "overridden mean must be finite",
            });
        }
        if !(statistical_std.is_finite() && statistical_std >= 0.0) {
            return Err(OedError::OverrideInvalid {
                what: "overridden statistical uncertainty must be finite and nonnegative",
            });
        }
        let overridden = Uncertainty {
            numerical: target.uncertainty.numerical,
            statistical: statistical_std,
            model: target.uncertainty.model,
        };
        Ok(Self {
            menu,
            index,
            mean,
            total_std: overridden.total_std(),
        })
    }

    /// FULL expected opportunity loss over the menu with this view's
    /// substitution — same evaluator, same tie-break, zero allocation
    /// (bead sj31i.5: action valuation and the robustness gate share
    /// ONE decision algebra, so an action's value is exactly the full
    /// loss it removes).
    fn opportunity_loss_checked(&self) -> Result<f64, OedError> {
        let value = expected_opportunity_loss_by(
            self.menu.len(),
            &|idx| {
                if idx == self.index {
                    self.mean
                } else {
                    self.menu.estimates[idx].mean
                }
            },
            &|idx| {
                if idx == self.index {
                    self.total_std
                } else {
                    self.menu.estimates[idx].uncertainty.total_std()
                }
            },
        );
        if value.is_finite() && value >= 0.0 {
            Ok(canonicalize_zero(value))
        } else {
            Err(OedError::NonFiniteComputation { quantity: "EVPI" })
        }
    }
}

/// Stable scalar Kalman variance update `P' = P R / (P + R)` without the
/// overflowing intermediate `P * R`. Candidate construction and Belief enforce
/// the input domains; this check protects the independently callable planner
/// path from a derived floating-point failure.
fn predicted_posterior_variance(prior: f64, noise: f64) -> Result<f64, OedError> {
    if prior == 0.0 {
        return Ok(0.0);
    }
    if !prior.is_finite() || prior < 0.0 || !noise.is_finite() || noise <= 0.0 {
        return Err(OedError::NonFiniteComputation {
            quantity: "sensor posterior variance inputs",
        });
    }
    // Divide the smaller operand only by a value at least as large. This is
    // algebraically identical to `P R / (P + R)`, avoids the overflowing
    // product, and—unlike scaling both inputs by `max(P, R)`—does not erase a
    // representable subnormal posterior when P and R span the full exponent
    // range.
    let posterior = if prior <= noise {
        prior / (1.0 + prior / noise)
    } else {
        noise / (1.0 + noise / prior)
    };
    if !posterior.is_finite() || posterior < 0.0 || posterior > prior {
        return Err(OedError::NonFiniteComputation {
            quantity: "predicted sensor posterior variance",
        });
    }
    Ok(canonicalize_zero(posterior))
}

fn sensor_actions(candidates: &[Candidate], beliefs: &[Belief]) -> Result<Vec<Action>, OedError> {
    if candidates.len() != beliefs.len() {
        return Err(OedError::NonFiniteComputation {
            quantity: "candidate/belief cardinality",
        });
    }
    candidates
        .iter()
        .zip(beliefs)
        .map(|(candidate, belief)| {
            let prior = belief.variance(0).map_err(OedError::BeliefInvariant)?;
            let posterior = predicted_posterior_variance(prior, candidate.sensor_noise_variance)?;
            let reduction = if prior == 0.0 {
                0.0
            } else {
                let std_ratio = (posterior / prior).clamp(0.0, 1.0).sqrt();
                canonicalize_zero((1.0 - std_ratio).clamp(0.0, 1.0))
            };
            Ok(Action {
                name: format!("measure-{}", candidate.name),
                kind: ActionKind::Sample,
                target_design: candidate.name.clone(),
                reduction,
                cost: candidate.sensor_cost,
            })
        })
        .collect()
}

/// Outcome-integrated value of a scalar Gaussian sensor action. The posterior
/// variance is the exact declared Kalman-model update. Posterior-mean movement
/// is integrated under its pre-posterior Gaussian distribution with the fixed
/// rule above, so a noisy sensor cannot inherit a fictitious universal effect.
fn expected_sensor_action_value(
    menu: &CanonicalDesignMenu,
    action: &Action,
    before: f64,
) -> Result<ActionValue, OedError> {
    let target =
        menu.index_of(&action.target_design)
            .ok_or_else(|| OedError::UnknownRecommendation {
                action: action.name.clone(),
            })?;
    let prior_mean = menu.estimates()[target].mean;
    if action.kind != ActionKind::Sample {
        return Err(OedError::NonFiniteComputation {
            quantity: "non-sensor action in sensor planner",
        });
    }
    let prior_std = menu.estimates()[target].uncertainty.total_std();
    let posterior_statistical =
        menu.estimates()[target].uncertainty.statistical * (1.0 - action.reduction).clamp(0.0, 1.0);
    // The overridden target's total uncertainty is fixed across the
    // quadrature; only its mean varies per node. Precompute it through
    // one throwaway view so the value is IDENTICAL to what each node's
    // view uses.
    let posterior_std =
        MeanOverrideView::new(menu, target, prior_mean, posterior_statistical)?.total_std;
    let mean_shift_variance = (prior_std * prior_std - posterior_std * posterior_std).max(0.0);
    let mean_shift_std = mean_shift_variance.sqrt();

    let mut expected_remaining_evpi = 0.0;
    for (normal_node, probability_weight) in NORMAL_EXPECTATION_RULE {
        let posterior_mean = prior_mean + normal_node * mean_shift_std;
        if !posterior_mean.is_finite() {
            return Err(OedError::NonFiniteComputation {
                quantity: "predictive posterior mean",
            });
        }
        // One-index substitution over the immutable menu: no clone, no
        // sort, no scratch mutation to restore on any failure path.
        let view = MeanOverrideView::new(
            menu,
            target,
            canonicalize_zero(posterior_mean),
            posterior_statistical,
        )?;
        expected_remaining_evpi += probability_weight * view.opportunity_loss_checked()?;
    }
    // Preserve the exact identity map after executing the declared fixed-shape
    // quadrature work. Summing nine identical weighted EVPI values can land a
    // single ulp below `before`; that rounding artifact is not sensor value.
    if action.reduction == 0.0 {
        expected_remaining_evpi = before;
    }
    if !expected_remaining_evpi.is_finite() || expected_remaining_evpi < 0.0 {
        return Err(OedError::NonFiniteComputation {
            quantity: "expected posterior EVPI",
        });
    }
    let value = canonicalize_zero((before - expected_remaining_evpi).max(0.0));
    let value_per_cost = if action.cost.is_finite() && action.cost > 0.0 {
        value / action.cost
    } else {
        0.0
    };
    if !value.is_finite() || !value_per_cost.is_finite() {
        return Err(OedError::NonFiniteComputation {
            quantity: "sensor action value",
        });
    }
    Ok(ActionValue {
        action: action.name.clone(),
        value,
        cost: action.cost,
        value_per_cost,
    })
}

fn precision_allocation(candidates: &[Candidate]) -> Result<(Vec<(String, f64)>, usize), OedError> {
    let features: Vec<Feature> = candidates
        .iter()
        .filter(|candidate| candidate.prior_var > 0.0)
        .map(|candidate| Feature {
            name: candidate.name.clone(),
            sensitivity: candidate.prior_var.sqrt(),
            sensitivity_color: ColorRank::Estimated,
            cost_coeff: candidate.sensor_cost,
            baseline_tolerance: 0.1,
        })
        .collect();
    let positive_prior_candidates = features.len();
    let allocated: BTreeMap<String, f64> = if features.is_empty() {
        BTreeMap::new()
    } else {
        allocate(&features, 0.02, 3.0)
            .map_err(|_| OedError::AllocationFailed)?
            .items
            .into_iter()
            .map(|item| (item.name, item.tolerance))
            .collect()
    };

    let allocation = candidates
        .iter()
        .map(|candidate| {
            if candidate.prior_var == 0.0 {
                Ok((candidate.name.clone(), f64::INFINITY))
            } else {
                let tolerance = allocated.get(&candidate.name).copied().ok_or_else(|| {
                    OedError::MissingAllocation {
                        candidate: candidate.name.clone(),
                    }
                })?;
                if !tolerance.is_finite() || tolerance <= 0.0 {
                    return Err(OedError::NonFiniteComputation {
                        quantity: "allocated tolerance",
                    });
                }
                Ok((candidate.name.clone(), tolerance))
            }
        })
        .collect::<Result<Vec<_>, OedError>>()?;
    Ok((allocation, positive_prior_candidates))
}

fn push_bytes(output: &mut Vec<u8>, value: &[u8]) {
    output.extend_from_slice(&(value.len() as u64).to_le_bytes());
    output.extend_from_slice(value);
}

fn push_str(output: &mut Vec<u8>, value: &str) {
    push_bytes(output, value.as_bytes());
}

#[derive(Debug, Clone, Copy)]
struct ReportIdentityOutputs<'a> {
    placements: &'a [String],
    sensors_placed: usize,
    prior_total_variance: f64,
    posterior_total_variance: f64,
    variance_reduction: f64,
    initial_evpi: f64,
    final_evpi: f64,
    decision_robust: bool,
    chosen_design: &'a str,
    allocation: &'a [(String, f64)],
    evpi_trace: &'a [f64],
    posteriors: &'a [PosteriorSummary],
    assimilation_colors: &'a [Color],
    variance_color_dispersion: f64,
    evpi_color_dispersion: f64,
}

#[derive(Debug, Clone, Copy)]
struct ReportIdentitySource<'a> {
    candidates: &'a [Candidate],
    objective_spec: ObjectiveSpec,
    threshold: f64,
    max_sensors: usize,
    outputs: ReportIdentityOutputs<'a>,
    plan: CampaignWorkPlan,
    byte_plan: CampaignBytePlan,
    realized: CampaignRealizedWorkPlan,
}

fn push_estimated_color_descriptor(output: &mut Vec<u8>, quantity: &str, dispersion: f64) {
    push_str(output, quantity);
    output.extend_from_slice(&fs_evidence::COLOR_ALGEBRA_VERSION.to_le_bytes());
    push_str(output, "Estimated");
    output.extend_from_slice(&dispersion.to_bits().to_le_bytes());
}

fn report_identity(
    quantity: &str,
    source: &ReportIdentitySource<'_>,
    progress: &mut CampaignProgress,
    cx: &Cx<'_>,
) -> Result<String, OedError> {
    report_identity_with_rule(quantity, source, &NORMAL_EXPECTATION_RULE, progress, cx)
}

#[allow(clippy::too_many_lines)] // One canonical manifest keeps field order auditable.
fn report_identity_with_rule(
    quantity: &str,
    source: &ReportIdentitySource<'_>,
    expectation_rule: &[(f64, f64)],
    progress: &mut CampaignProgress,
    cx: &Cx<'_>,
) -> Result<String, OedError> {
    let candidates = source.candidates;
    let objective_spec = source.objective_spec;
    let threshold = source.threshold;
    let max_sensors = source.max_sensors;
    let outputs = source.outputs;
    let plan = source.plan;
    let realized = source.realized;
    let mut canonical = Vec::new();
    canonical.extend_from_slice(&OED_REPORT_IDENTITY_VERSION.to_le_bytes());
    canonical.extend_from_slice(&objective_spec.canonical_bytes());
    push_str(&mut canonical, cx.mode().name());
    let stream = cx.stream_key();
    for value in [stream.seed, stream.kernel_id, stream.tile, stream.iteration] {
        canonical.extend_from_slice(&value.to_le_bytes());
    }
    let budget = cx.budget();
    match budget.deadline {
        Some(deadline) => {
            canonical.push(1);
            canonical.extend_from_slice(&deadline.as_nanos().to_le_bytes());
        }
        None => canonical.push(0),
    }
    canonical.extend_from_slice(&budget.poll_quota.to_le_bytes());
    match budget.cost_quota {
        Some(cost_quota) => {
            canonical.push(1);
            canonical.extend_from_slice(&cost_quota.to_le_bytes());
        }
        None => canonical.push(0),
    }
    canonical.push(budget.priority);
    for value in [
        plan.candidates as u128,
        plan.max_sensors as u128,
        plan.action_design_evaluations as u128,
        plan.setup_work_units,
        plan.per_placement_work_units,
        plan.maximum_finalization_work_units,
        plan.admitted_work_units,
    ] {
        canonical.extend_from_slice(&value.to_le_bytes());
    }
    for value in realized.identity_fields() {
        canonical.extend_from_slice(&value.to_le_bytes());
    }
    // v7: the admitted byte envelope and its policy are identity-
    // semantic — agreement on science under a different byte plan is a
    // different admitted artifact.
    for value in source.byte_plan.identity_fields() {
        canonical.extend_from_slice(&value.to_le_bytes());
    }
    canonical.extend_from_slice(&CAMPAIGN_BYTE_POLICY_VERSION.to_le_bytes());
    canonical.extend_from_slice(&CAMPAIGN_PLANNING_POLICY_VERSION.to_le_bytes());
    canonical.extend_from_slice(&CAMPAIGN_POLL_POLICY_VERSION.to_le_bytes());
    canonical.extend_from_slice(&fs_voi::EVPI_SEMANTICS_VERSION.to_le_bytes());
    canonical.extend_from_slice(&(expectation_rule.len() as u64).to_le_bytes());
    for &(node, weight) in expectation_rule {
        canonical.extend_from_slice(&node.to_bits().to_le_bytes());
        canonical.extend_from_slice(&weight.to_bits().to_le_bytes());
    }
    canonical.extend_from_slice(&(CAMPAIGN_RECORD_POLL_STRIDE as u64).to_le_bytes());
    canonical.extend_from_slice(&(CAMPAIGN_ACTION_POLL_STRIDE as u64).to_le_bytes());
    push_str(&mut canonical, quantity);
    canonical.extend_from_slice(&(candidates.len() as u64).to_le_bytes());
    for (index, candidate) in candidates.iter().enumerate() {
        if index.is_multiple_of(CAMPAIGN_RECORD_POLL_STRIDE) {
            progress.checkpoint(cx, plan, "report identity candidates")?;
        }
        push_str(&mut canonical, &candidate.name);
        for value in [
            candidate.truth,
            candidate.prior_mean,
            candidate.prior_var,
            candidate.sensor_noise_variance,
            candidate.sensor_cost,
        ] {
            canonical.extend_from_slice(&value.to_bits().to_le_bytes());
        }
        progress.advance(1);
    }
    canonical.extend_from_slice(&threshold.to_bits().to_le_bytes());
    canonical.extend_from_slice(&(max_sensors as u64).to_le_bytes());
    canonical.extend_from_slice(&(outputs.placements.len() as u64).to_le_bytes());
    for (index, placement) in outputs.placements.iter().enumerate() {
        if index.is_multiple_of(CAMPAIGN_RECORD_POLL_STRIDE) {
            progress.checkpoint(cx, plan, "report identity placements")?;
        }
        push_str(&mut canonical, placement);
        progress.advance(1);
    }
    canonical.extend_from_slice(&(outputs.sensors_placed as u64).to_le_bytes());
    for value in [
        outputs.prior_total_variance,
        outputs.posterior_total_variance,
        outputs.variance_reduction,
        outputs.initial_evpi,
        outputs.final_evpi,
    ] {
        canonical.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    canonical.push(u8::from(outputs.decision_robust));
    push_str(&mut canonical, outputs.chosen_design);
    canonical.extend_from_slice(&(outputs.allocation.len() as u64).to_le_bytes());
    for (index, (name, tolerance)) in outputs.allocation.iter().enumerate() {
        if index.is_multiple_of(CAMPAIGN_RECORD_POLL_STRIDE) {
            progress.checkpoint(cx, plan, "report identity allocation")?;
        }
        push_str(&mut canonical, name);
        canonical.extend_from_slice(&tolerance.to_bits().to_le_bytes());
        progress.advance(1);
    }
    canonical.extend_from_slice(&(outputs.evpi_trace.len() as u64).to_le_bytes());
    for (index, value) in outputs.evpi_trace.iter().enumerate() {
        if index.is_multiple_of(CAMPAIGN_RECORD_POLL_STRIDE) {
            progress.checkpoint(cx, plan, "report identity EVPI trace")?;
        }
        canonical.extend_from_slice(&value.to_bits().to_le_bytes());
        progress.advance(1);
    }
    canonical.extend_from_slice(&(outputs.posteriors.len() as u64).to_le_bytes());
    for (index, posterior) in outputs.posteriors.iter().enumerate() {
        if index.is_multiple_of(CAMPAIGN_RECORD_POLL_STRIDE) {
            progress.checkpoint(cx, plan, "report identity posteriors")?;
        }
        push_str(&mut canonical, &posterior.name);
        canonical.extend_from_slice(&posterior.mean.to_bits().to_le_bytes());
        canonical.extend_from_slice(&posterior.variance.to_bits().to_le_bytes());
        progress.advance(1);
    }
    canonical.extend_from_slice(&(outputs.assimilation_colors.len() as u64).to_le_bytes());
    for (index, color) in outputs.assimilation_colors.iter().enumerate() {
        if index.is_multiple_of(CAMPAIGN_RECORD_POLL_STRIDE) {
            progress.checkpoint(cx, plan, "report identity assimilation colors")?;
        }
        push_bytes(&mut canonical, &color.canonical_bytes());
        progress.advance(1);
    }
    // The estimator strings are the hashes being derived and therefore cannot
    // appear in their own preimage. Bind the stable color algebra, variant, and
    // both dispersions; the sealed report is constructed only after both
    // estimator strings have been derived from this complete source.
    push_estimated_color_descriptor(
        &mut canonical,
        "posterior-variance",
        outputs.variance_color_dispersion,
    );
    push_estimated_color_descriptor(&mut canonical, "evpi", outputs.evpi_color_dispersion);
    progress.checkpoint(cx, plan, "report identity hash")?;
    // Hash bytes are charged at their exact deterministic preimage
    // length plus the digest emission (bead sj31i.62).
    progress.charge_bytes(
        "report identity hash",
        (canonical.len() as u128).saturating_add(64),
    )?;
    let identity = format!(
        "sensorforge-{quantity}:v{OED_REPORT_IDENTITY_VERSION}:{}",
        fs_blake3::hash_domain(REPORT_ID_DOMAIN, &canonical)
    );
    progress.retain_bytes("report identity hash", identity.len() as u128)?;
    progress.advance(1);
    debug_assert!(color_leaf_identity_reason(&identity).is_none());
    Ok(identity)
}

fn validate_campaign(
    candidates: &[Candidate],
    threshold: ObjectiveValue,
    max_sensors: usize,
) -> Result<(f64, CampaignWorkPlan, CampaignBytePlan, Vec<Candidate>), OedError> {
    if candidates.is_empty() {
        return Err(OedError::NoCandidates);
    }
    if candidates.len() > MAX_CAMPAIGN_CANDIDATES {
        return Err(OedError::TooManyCandidates {
            count: candidates.len(),
            max: MAX_CAMPAIGN_CANDIDATES,
        });
    }
    if max_sensors > MAX_CAMPAIGN_SENSORS {
        return Err(OedError::TooManySensors {
            count: max_sensors,
            max: MAX_CAMPAIGN_SENSORS,
        });
    }
    let plan = CampaignWorkPlan::checked(candidates.len(), max_sensors)?;
    if threshold.quantity.value < 0.0 {
        return Err(OedError::InvalidThreshold);
    }
    let objective_spec = candidates[0].objective_spec;
    let decision_spec = objective_spec.decision_spec();
    if threshold.spec != decision_spec {
        return Err(OedError::ThresholdSchemaMismatch {
            actual: threshold.spec,
            expected: decision_spec,
        });
    }
    let mut names = BTreeSet::new();
    let mut sum_name = 0usize;
    let mut max_name = 0usize;
    for candidate in candidates {
        if candidate.objective_spec != objective_spec {
            return Err(OedError::ObjectiveSchemaMismatch {
                candidate: candidate.name.clone(),
                actual: candidate.objective_spec,
                expected: objective_spec,
            });
        }
        if !names.insert(candidate.name.as_str()) {
            return Err(OedError::DuplicateCandidate {
                name: candidate.name.clone(),
            });
        }
        sum_name = sum_name.saturating_add(candidate.name.len());
        max_name = max_name.max(candidate.name.len());
    }
    // The byte envelope is preflighted from the ACTUAL validated names
    // (bead sj31i.62): admission evaluates the same charge formulas
    // every later seam uses, at the worst-case shape.
    let byte_plan = CampaignBytePlan::checked(candidates.len(), max_sensors, sum_name, max_name)?;
    // Canonicalize unique candidate identity and order EXACTLY ONCE at
    // admission (bead sj31i.62); every derived belief/estimate/action
    // sequence inherits this order, so no later phase re-sorts. The
    // sort shares the validation scan's accounting unit — a unit is a
    // bounded record visit, not an instruction count, and the per-call
    // clone-and-sort work this replaces was never separately charged.
    // Its byte bound is the ledger's opening charge.
    let mut canonical = Vec::new();
    canonical.try_reserve_exact(candidates.len()).map_err(|_| {
        OedError::OutputAllocationRefused {
            what: "canonical candidate menu",
        }
    })?;
    canonical.extend_from_slice(candidates);
    canonical.sort_by(|left, right| left.name.cmp(&right.name));
    Ok((threshold.quantity.value, plan, byte_plan, canonical))
}

struct CampaignState {
    beliefs: Vec<Belief>,
    placements: Vec<String>,
    assimilation_colors: Vec<Color>,
    evpi_trace: Vec<f64>,
    decision_robust: bool,
    action_rounds: usize,
}

fn recommend_with_cancellation(
    menu: &CanonicalDesignMenu,
    actions: &[Action],
    current_evpi: f64,
    threshold: f64,
    plan: CampaignWorkPlan,
    shape: ByteShape,
    progress: &mut CampaignProgress,
    cx: &Cx<'_>,
) -> Result<Recommendation, OedError> {
    if current_evpi <= threshold {
        return Ok(Recommendation::Stop {
            reason: format!("decision robust: EVPI {current_evpi:.3e} <= {threshold:.3e}"),
        });
    }

    let mut best = None;
    for action in actions {
        progress.checkpoint(cx, plan, "action-value tile")?;
        let value = expected_sensor_action_value(menu, action, current_evpi)?;
        progress.advance((menu.len() as u128) * (ACTION_EVALUATION_FACTOR as u128) + 1);
        progress.charge_bytes("action-value tile", shape.action_evaluation())?;
        if value.value <= 0.0 || value.value_per_cost <= 0.0 {
            continue;
        }
        let replace = best.as_ref().is_none_or(|current: &ActionValue| {
            match value.value_per_cost.total_cmp(&current.value_per_cost) {
                std::cmp::Ordering::Greater => true,
                std::cmp::Ordering::Equal => value.action < current.action,
                std::cmp::Ordering::Less => false,
            }
        });
        if replace {
            best = Some(value);
        }
    }
    progress.checkpoint(cx, plan, "action-value drain")?;

    Ok(match best {
        Some(value) => Recommendation::Act {
            action: value.action,
            value_per_cost: value.value_per_cost,
        },
        None => Recommendation::Stop {
            reason: "no action changes the decision".to_string(),
        },
    })
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn execute_placements(
    candidates: &[Candidate],
    threshold: f64,
    max_sensors: usize,
    mut beliefs: Vec<Belief>,
    mut menu: CanonicalDesignMenu,
    initial_evpi: f64,
    plan: CampaignWorkPlan,
    shape: ByteShape,
    progress: &mut CampaignProgress,
    cx: &Cx<'_>,
) -> Result<CampaignState, OedError> {
    // Retained campaign outputs reserve their admitted capacity
    // fallibly before any placement runs (bead sj31i.62).
    let mut placements: Vec<String> = Vec::new();
    placements
        .try_reserve_exact(max_sensors)
        .map_err(|_| OedError::OutputAllocationRefused {
            what: "sensor placements",
        })?;
    let mut assimilation_colors: Vec<Color> = Vec::new();
    assimilation_colors
        .try_reserve_exact(max_sensors)
        .map_err(|_| OedError::OutputAllocationRefused {
            what: "assimilation colors",
        })?;
    let mut evpi_trace: Vec<f64> = Vec::new();
    evpi_trace
        .try_reserve_exact(max_sensors.saturating_add(1))
        .map_err(|_| OedError::OutputAllocationRefused { what: "EVPI trace" })?;
    evpi_trace.push(initial_evpi);
    let mut decision_robust = false;
    let mut action_rounds = 0;
    let mut current_evpi = initial_evpi;

    loop {
        if current_evpi <= threshold {
            decision_robust = true;
            break;
        }
        if placements.len() >= max_sensors {
            break;
        }
        progress.checkpoint(cx, plan, "action construction")?;
        let actions = sensor_actions(candidates, &beliefs)?;
        progress.advance(candidates.len() as u128);
        progress.charge_bytes("action construction", shape.action_construction())?;
        progress.checkpoint(cx, plan, "action construction drain")?;
        let recommendation = recommend_with_cancellation(
            &menu,
            &actions,
            current_evpi,
            threshold,
            plan,
            shape,
            progress,
            cx,
        )?;
        action_rounds += 1;
        let Recommendation::Act { action, .. } = recommendation else {
            break;
        };
        progress.checkpoint(cx, plan, "chosen-action lookup")?;
        let idx = actions
            .iter()
            .position(|candidate| candidate.name == action)
            .ok_or(OedError::UnknownRecommendation { action })?;
        progress.advance(candidates.len() as u128);
        progress.charge_bytes("chosen-action lookup", shape.chosen_lookup())?;
        let observation = point_sensor(
            0,
            1,
            candidates[idx].truth,
            candidates[idx].sensor_noise_variance,
            format!(
                "sensor-{}-q{}",
                candidates[idx].name,
                candidates[idx].objective_spec.identity_hex()
            ),
        )
        .map_err(|source| OedError::Assimilation {
            candidate: candidates[idx].name.clone(),
            source,
        })?;
        progress.advance(1);
        progress.charge_bytes(
            "sensor observation",
            SENSOR_IDENTITY_FIXED_BYTES
                .saturating_add(candidates[idx].name.len() as u128)
                .saturating_add(RECORD_SCALAR_BYTE_UNITS),
        )?;
        let next_count = placements.len() + 1;
        let posterior = assimilate_colored_with_shared_poll_quota(
            &beliefs[idx],
            std::slice::from_ref(&observation),
            "sensor_count",
            0.0,
            next_count as f64,
            cx,
            &mut progress.polls_remaining,
        )
        .map_err(|source| {
            if matches!(source, AssimError::Cancelled { .. }) {
                OedError::AssimilationCancelled {
                    candidate: candidates[idx].name.clone(),
                    completed_placements: progress.completed_placements,
                    completed_work_units: progress.completed_work_units,
                    admitted_work_units: plan.admitted_work_units,
                    source: Box::new(source),
                }
            } else {
                OedError::Assimilation {
                    candidate: candidates[idx].name.clone(),
                    source,
                }
            }
        })?;
        progress.advance(1);
        // Request -> drain -> finalize: do not publish the scratch posterior
        // into campaign state until the lower-layer transaction has drained
        // and this deterministic commit boundary is still live.
        progress.checkpoint(cx, plan, "placement commit")?;
        progress.charge_bytes("placement commit", RECORD_SCALAR_BYTE_UNITS)?;
        progress.retain_bytes(
            "placement commit",
            32u128
                .saturating_add(candidates[idx].name.len() as u128)
                .saturating_add(8),
        )?;
        beliefs[idx] = posterior.belief().clone();
        assimilation_colors.push(posterior.color().clone());
        placements.push(candidates[idx].name.clone());
        progress.completed_placements = placements.len();

        progress.checkpoint(cx, plan, "posterior estimate refresh")?;
        menu = CanonicalDesignMenu::from_canonical(to_estimates(candidates, &beliefs)?)?;
        progress.advance(candidates.len() as u128);
        progress.charge_bytes("posterior estimate refresh", shape.menu_build())?;
        progress.checkpoint(cx, plan, "posterior EVPI refresh")?;
        current_evpi = menu.full_opportunity_loss_checked()?;
        progress.advance(candidates.len() as u128);
        progress.charge_bytes("posterior EVPI refresh", shape.full_scan())?;
        evpi_trace.push(current_evpi);
    }

    Ok(CampaignState {
        beliefs,
        placements,
        assimilation_colors,
        evpi_trace,
        decision_robust,
        action_rounds,
    })
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn finish_report(
    candidates: &[Candidate],
    threshold: f64,
    max_sensors: usize,
    prior_total_variance: f64,
    initial_evpi: f64,
    state: CampaignState,
    plan: CampaignWorkPlan,
    shape: ByteShape,
    progress: &mut CampaignProgress,
    cx: &Cx<'_>,
) -> Result<OedReport, OedError> {
    progress.checkpoint(cx, plan, "final estimate summary")?;
    let menu = CanonicalDesignMenu::from_canonical(to_estimates(candidates, &state.beliefs)?)?;
    progress.advance(candidates.len() as u128);
    progress.charge_bytes("final estimate summary", shape.menu_build())?;
    progress.checkpoint(cx, plan, "final EVPI")?;
    let final_evpi = menu.full_opportunity_loss_checked()?;
    progress.advance(candidates.len() as u128);
    progress.charge_bytes("final EVPI", shape.full_scan())?;
    progress.checkpoint(cx, plan, "final variance")?;
    let posterior_total_variance = total_variance(&state.beliefs)?;
    progress.advance(candidates.len() as u128);
    progress.charge_bytes("final variance", shape.scalar_reads())?;
    progress.checkpoint(cx, plan, "chosen-design reduction")?;
    let chosen_design = menu
        .estimates()
        .iter()
        .min_by(|a, b| a.mean.total_cmp(&b.mean).then_with(|| a.name.cmp(&b.name)))
        .map(|design| design.name.clone())
        .ok_or(OedError::NonFiniteComputation {
            quantity: "chosen design",
        })?;
    progress.advance(candidates.len() as u128);
    progress.charge_bytes(
        "chosen-design reduction",
        shape.n.saturating_mul(shape.max_name.saturating_add(8)),
    )?;
    progress.retain_bytes("chosen-design reduction", chosen_design.len() as u128)?;
    progress.checkpoint(cx, plan, "precision allocation")?;
    let (allocation, positive_prior_candidates) = precision_allocation(candidates)?;
    progress.advance((candidates.len() as u128) * 2 + (positive_prior_candidates as u128) * 5);
    progress.charge_bytes("precision allocation", shape.scan())?;
    let allocation_retained: u128 = allocation
        .iter()
        .map(|(name, _)| (name.len() as u128).saturating_add(8))
        .sum();
    progress.retain_bytes("precision allocation", allocation_retained)?;
    progress.checkpoint(cx, plan, "posterior summaries")?;
    let mut posteriors: Vec<PosteriorSummary> = Vec::new();
    posteriors
        .try_reserve_exact(candidates.len())
        .map_err(|_| OedError::OutputAllocationRefused {
            what: "posterior summaries",
        })?;
    for (candidate, belief) in candidates.iter().zip(&state.beliefs) {
        posteriors.push(PosteriorSummary {
            name: candidate.name.clone(),
            objective_spec: candidate.objective_spec,
            mean: belief
                .component_mean(0)
                .map_err(OedError::BeliefInvariant)?,
            variance: belief.variance(0).map_err(OedError::BeliefInvariant)?,
        });
    }
    progress.advance(candidates.len() as u128);
    progress.retain_bytes(
        "posterior summaries",
        shape.sum_name.saturating_add(
            shape
                .n
                .saturating_mul(16u128.saturating_add(OBJECTIVE_SCHEMA_BYTES)),
        ),
    )?;
    let variance_reduction = if prior_total_variance == 0.0 {
        0.0
    } else {
        let reduction = (prior_total_variance - posterior_total_variance) / prior_total_variance;
        if !reduction.is_finite() {
            return Err(OedError::NonFiniteComputation {
                quantity: "variance reduction",
            });
        }
        canonicalize_zero(reduction)
    };
    let sensors_placed = state.placements.len();
    let realized = plan
        .realized(
            sensors_placed,
            state.action_rounds,
            positive_prior_candidates,
        )
        .ok_or(OedError::WorkPlanMismatch {
            completed_work_units: progress.completed_work_units,
            realized_work_units: u128::MAX,
            admitted_work_units: plan.admitted_work_units,
        })?;
    let variance_color_dispersion = f64::INFINITY;
    let evpi_color_dispersion = final_evpi;
    let (variance_identity, evpi_identity) = {
        let source = ReportIdentitySource {
            candidates,
            objective_spec: candidates[0].objective_spec,
            threshold,
            max_sensors,
            outputs: ReportIdentityOutputs {
                placements: &state.placements,
                sensors_placed,
                prior_total_variance,
                posterior_total_variance,
                variance_reduction,
                initial_evpi,
                final_evpi,
                decision_robust: state.decision_robust,
                chosen_design: &chosen_design,
                allocation: &allocation,
                evpi_trace: &state.evpi_trace,
                posteriors: &posteriors,
                assimilation_colors: &state.assimilation_colors,
                variance_color_dispersion,
                evpi_color_dispersion,
            },
            plan,
            byte_plan: progress.byte_plan,
            realized,
        };
        let variance_identity = report_identity("posterior-variance", &source, progress, cx)?;
        let evpi_identity = report_identity("evpi", &source, progress, cx)?;
        (variance_identity, evpi_identity)
    };

    progress.checkpoint(cx, plan, "report publication")?;
    progress.retain_bytes("objective schema", OBJECTIVE_SCHEMA_BYTES)?;
    progress.advance(1);
    progress.finish(plan, realized)?;

    Ok(OedReport {
        objective_spec: candidates[0].objective_spec,
        admitted_byte_units: progress.byte_plan.admitted_byte_units,
        consumed_byte_units: progress.charged_byte_units,
        retained_byte_units: progress.retained_byte_units,
        sensors_placed,
        placements: state.placements,
        prior_total_variance,
        posterior_total_variance,
        variance_reduction,
        initial_evpi,
        final_evpi,
        decision_robust: state.decision_robust,
        chosen_design,
        allocation,
        evpi_trace: state.evpi_trace,
        posteriors,
        assimilation_colors: state.assimilation_colors,
        variance_color: Color::Estimated {
            estimator: variance_identity,
            dispersion: variance_color_dispersion,
        },
        evpi_color: Color::Estimated {
            estimator: evpi_identity,
            dispersion: evpi_color_dispersion,
        },
    })
}

/// Run the SensorForge campaign under an explicit execution context; stop when
/// EVPI <= `threshold` or after `max_sensors` placements. The threshold must
/// carry the objective's decision-difference schema.
///
/// The complete worst-case work bound is checked before scientific work starts,
/// and the exact realized early-stop shape is checked before publication. The
/// initial STOP condition is evaluated even when `max_sensors == 0`, and
/// cancellation is polled at deterministic action/record boundaries.
///
/// # Errors
/// Returns [`OedError`] for invalid campaign bounds, duplicate candidate names,
/// observed cancellation, a lower-layer assimilation/allocation failure, or a
/// non-finite derived value. A cancellation never returns a partial report.
pub fn run_campaign(
    candidates: &[Candidate],
    threshold: ObjectiveValue,
    max_sensors: usize,
    cx: &Cx<'_>,
) -> Result<OedReport, OedError> {
    let (threshold, plan, byte_plan, candidates) =
        validate_campaign(candidates, threshold, max_sensors)?;
    let candidates = candidates.as_slice();
    let shape = ByteShape::of(candidates);
    let mut progress = CampaignProgress::admit(
        cx,
        candidates.len() as u128,
        u64::try_from(plan.admitted_work_units).unwrap_or(u64::MAX),
        byte_plan,
    )?;
    progress.checkpoint(cx, plan, "campaign admission")?;
    let beliefs: Vec<Belief> = candidates
        .iter()
        .map(|c| Belief::scalar(c.prior_mean, c.prior_var))
        .collect::<Result<Vec<_>, _>>()
        .map_err(OedError::BeliefInvariant)?;
    progress.advance(candidates.len() as u128);
    progress.charge_bytes("prior beliefs", shape.records())?;
    progress.checkpoint(cx, plan, "prior variance")?;
    let prior_total_variance = total_variance(&beliefs)?;
    progress.advance(candidates.len() as u128);
    progress.charge_bytes("prior variance", shape.scalar_reads())?;
    progress.checkpoint(cx, plan, "initial estimates")?;
    let menu = CanonicalDesignMenu::from_canonical(to_estimates(candidates, &beliefs)?)?;
    progress.advance(candidates.len() as u128);
    progress.charge_bytes("initial estimates", shape.menu_build())?;
    progress.checkpoint(cx, plan, "initial EVPI")?;
    let initial_evpi = menu.full_opportunity_loss_checked()?;
    progress.advance(candidates.len() as u128);
    progress.charge_bytes("initial EVPI", shape.full_scan())?;
    let state = execute_placements(
        candidates,
        threshold,
        max_sensors,
        beliefs,
        menu,
        initial_evpi,
        plan,
        shape,
        &mut progress,
        cx,
    )?;
    finish_report(
        candidates,
        threshold,
        max_sensors,
        prior_total_variance,
        initial_evpi,
        state,
        plan,
        shape,
        &mut progress,
        cx,
    )
}

/// The worked scenario: four designs with uncertain COST (lower is better). The
/// two cheapest (A, B) are close and uncertain — the decision hinges on
/// measuring THEM, not the clearly-costlier C or D.
pub fn demo_candidates() -> Result<Vec<Candidate>, CandidateError> {
    [
        ("A", 0.60, 0.60, 0.10, 0.01, 1.0),
        ("B", 0.65, 0.65, 0.12, 0.01, 1.0),
        ("C", 0.85, 0.85, 0.06, 0.01, 1.0),
        ("D", 1.10, 1.10, 0.04, 0.01, 1.0),
    ]
    .into_iter()
    .map(
        |(name, truth, prior_mean, prior_var, sensor_noise, sensor_cost)| {
            Candidate::new(
                name,
                ObjectiveValue::dimensionless(truth)?,
                ObjectiveValue::dimensionless(prior_mean)?,
                QtyAny::dimensionless(prior_var),
                QtyAny::dimensionless(sensor_noise),
                QtyAny::dimensionless(sensor_cost),
            )
        },
    )
    .collect()
}

#[cfg(test)]
mod tests {
    use super::{
        CAMPAIGN_PLANNING_POLICY_VERSION, CAMPAIGN_POLL_POLICY_VERSION, CampaignBytePlan,
        CampaignProgress, CampaignRealizedWorkPlan, CampaignWorkPlan, Candidate,
        NORMAL_EXPECTATION_RULE, OED_REPORT_IDENTITY_VERSION, ObjectiveSpec, PosteriorSummary,
        REPORT_ID_DOMAIN, ReportIdentityOutputs, ReportIdentitySource, canonicalize_zero,
        demo_candidates, predicted_posterior_variance, report_identity_with_rule,
    };
    use fs_evidence::Color;
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

    #[derive(Clone)]
    struct IdentityFixture {
        candidates: Vec<Candidate>,
        threshold: f64,
        max_sensors: usize,
        placements: Vec<String>,
        sensors_placed: usize,
        prior_total_variance: f64,
        posterior_total_variance: f64,
        variance_reduction: f64,
        initial_evpi: f64,
        final_evpi: f64,
        decision_robust: bool,
        chosen_design: String,
        allocation: Vec<(String, f64)>,
        evpi_trace: Vec<f64>,
        posteriors: Vec<PosteriorSummary>,
        assimilation_colors: Vec<Color>,
        variance_color_dispersion: f64,
        evpi_color_dispersion: f64,
        realized: CampaignRealizedWorkPlan,
    }

    impl IdentityFixture {
        fn new() -> Self {
            let candidates = demo_candidates().expect("demo candidates");
            let allocation = candidates
                .iter()
                .enumerate()
                .map(|(index, candidate)| (candidate.name().to_string(), 0.1 + index as f64 * 0.01))
                .collect();
            let posteriors = candidates
                .iter()
                .map(|candidate| PosteriorSummary {
                    name: candidate.name().to_string(),
                    objective_spec: candidate.objective_spec(),
                    mean: candidate.prior_mean().quantity().value,
                    variance: candidate.prior_variance().value,
                })
                .collect();
            let plan = CampaignWorkPlan::checked(candidates.len(), 1).expect("fixture work plan");
            let realized = plan
                .realized(1, 1, candidates.len())
                .expect("fixture realized work plan");
            Self {
                candidates,
                threshold: 0.1,
                max_sensors: 1,
                placements: vec!["A".to_string()],
                sensors_placed: 1,
                prior_total_variance: 0.32,
                posterior_total_variance: 0.20,
                variance_reduction: 0.375,
                initial_evpi: 0.4,
                final_evpi: 0.2,
                decision_robust: false,
                chosen_design: "A".to_string(),
                allocation,
                evpi_trace: vec![0.4, 0.2],
                posteriors,
                assimilation_colors: vec![Color::Estimated {
                    estimator: "sensor-A-v1".to_string(),
                    dispersion: 0.01,
                }],
                variance_color_dispersion: f64::INFINITY,
                evpi_color_dispersion: 0.2,
                realized,
            }
        }

        fn source(&self) -> ReportIdentitySource<'_> {
            ReportIdentitySource {
                candidates: &self.candidates,
                objective_spec: self.candidates[0].objective_spec,
                threshold: self.threshold,
                max_sensors: self.max_sensors,
                outputs: ReportIdentityOutputs {
                    placements: &self.placements,
                    sensors_placed: self.sensors_placed,
                    prior_total_variance: self.prior_total_variance,
                    posterior_total_variance: self.posterior_total_variance,
                    variance_reduction: self.variance_reduction,
                    initial_evpi: self.initial_evpi,
                    final_evpi: self.final_evpi,
                    decision_robust: self.decision_robust,
                    chosen_design: &self.chosen_design,
                    allocation: &self.allocation,
                    evpi_trace: &self.evpi_trace,
                    posteriors: &self.posteriors,
                    assimilation_colors: &self.assimilation_colors,
                    variance_color_dispersion: self.variance_color_dispersion,
                    evpi_color_dispersion: self.evpi_color_dispersion,
                },
                plan: CampaignWorkPlan::checked(self.candidates.len(), self.max_sensors)
                    .expect("fixture work plan"),
                byte_plan: self.byte_plan(),
                realized: self.realized,
            }
        }

        fn byte_plan(&self) -> CampaignBytePlan {
            let sum_name: usize = self
                .candidates
                .iter()
                .map(|candidate| candidate.name.len())
                .sum();
            let max_name = self
                .candidates
                .iter()
                .map(|candidate| candidate.name.len())
                .max()
                .unwrap_or(0);
            CampaignBytePlan::checked(self.candidates.len(), self.max_sensors, sum_name, max_name)
                .expect("fixture byte plan")
        }
    }

    fn with_test_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 1,
                    kernel_id: 2,
                    tile: 3,
                    iteration: 4,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    fn report_identities(fixture: &IdentityFixture, rule: &[(f64, f64)]) -> (String, String) {
        with_test_cx(|cx| {
            let source = fixture.source();
            let mut progress = CampaignProgress::admit(cx, 0, 0, fixture.byte_plan())
                .expect("identity fixture budget admits");
            let variance =
                report_identity_with_rule("posterior-variance", &source, rule, &mut progress, cx)
                    .expect("variance identity");
            let evpi = report_identity_with_rule("evpi", &source, rule, &mut progress, cx)
                .expect("EVPI identity");
            (variance, evpi)
        })
    }

    #[test]
    fn scalar_variance_update_preserves_representable_extreme_posteriors() {
        for (prior, noise) in [
            (1.0e300, 1.0e-300),
            (1.0e-300, 1.0e300),
            (f64::MAX, f64::from_bits(1)),
        ] {
            let posterior = predicted_posterior_variance(prior, noise)
                .expect("positive finite scalar variances have a posterior");
            assert!(posterior > 0.0, "a representable posterior was erased");
            assert!(posterior <= prior.min(noise));
        }
    }

    #[test]
    fn normal_expectation_rule_is_positive_normalized_and_symmetric() {
        let weight_sum: f64 = NORMAL_EXPECTATION_RULE
            .iter()
            .map(|(_, weight)| weight)
            .sum();
        assert!((weight_sum - 1.0).abs() <= 8.0 * f64::EPSILON);
        for (left, right) in NORMAL_EXPECTATION_RULE
            .iter()
            .zip(NORMAL_EXPECTATION_RULE.iter().rev())
        {
            assert!(left.1 > 0.0);
            assert_eq!(left.0.to_bits(), canonicalize_zero(-right.0).to_bits());
            assert_eq!(left.1.to_bits(), right.1.to_bits());
        }
    }

    #[test]
    fn report_identity_versions_and_final_work_shape_are_locked() {
        // v8 (bead sj31i.7): deliberate bump — the preimage now binds
        // the exact objective dimension and semantic schema.
        assert_eq!(OED_REPORT_IDENTITY_VERSION, 8);
        assert_eq!(REPORT_ID_DOMAIN, "org.frankensim.fs-oed-e2e.report.v6");
        // v4 (bead sj31i.5): full-EOL robustness evaluations.
        assert_eq!(CAMPAIGN_PLANNING_POLICY_VERSION, 4);
        assert_eq!(CAMPAIGN_POLL_POLICY_VERSION, 2);
        assert_eq!(super::CAMPAIGN_BYTE_POLICY_VERSION, 3);
        // v2 (bead sj31i.5): full multi-alternative opportunity loss.
        assert_eq!(fs_voi::EVPI_SEMANTICS_VERSION, 2);

        let plan = CampaignWorkPlan::checked(4, 12).expect("admitted work plan");
        assert_eq!(plan.maximum_finalization_work_units, 18 * 4 + 6 * 12 + 5);
        assert_eq!(plan.admitted_work_units, 2_545);
        assert_eq!(
            plan.realized(0, 0, 4)
                .expect("immediate STOP shape")
                .realized_work_units,
            97
        );
        assert_eq!(
            plan.realized(0, 1, 4)
                .expect("one completed zero-value action round")
                .realized_work_units,
            281
        );
        assert_eq!(
            plan.realized(12, 12, 4)
                .expect("full placement shape")
                .realized_work_units,
            plan.admitted_work_units
        );

        let identities = report_identities(&IdentityFixture::new(), &NORMAL_EXPECTATION_RULE);
        assert!(
            identities
                .0
                .starts_with("sensorforge-posterior-variance:v8:")
        );
        assert!(identities.1.starts_with("sensorforge-evpi:v8:"));
    }

    #[test]
    fn objective_schema_encoding_binds_all_six_axes_and_semantic_parameters() {
        use fs_qty::Dims;
        use fs_qty::semantic::{
            AngleDomain, CompositionBasis, QuantityKind, QuantitySpec, SemanticType, StrainBasis,
            StrainComponent, ValueForm,
        };

        let base = ObjectiveSpec::dimensional(Dims([1, -2, 3, -4, 5, -6]));
        assert_eq!(
            base.canonical_bytes(),
            [1, 1, 254, 3, 252, 5, 250, 0, 0, 0, 0, 0]
        );
        assert_eq!(
            QuantitySpec::from_canonical_bytes(&base.canonical_bytes()),
            Ok(base.quantity_spec())
        );
        for axis in 0..6 {
            let mut exponents = base.dims().0;
            exponents[axis] += 1;
            assert_ne!(
                base.canonical_bytes(),
                ObjectiveSpec::dimensional(Dims(exponents)).canonical_bytes()
            );
        }

        let semantic_types = [
            SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static),
            SemanticType::new(QuantityKind::TemperatureDifference, ValueForm::Static),
            SemanticType::new(
                QuantityKind::Angle(AngleDomain::Mechanical),
                ValueForm::Static,
            ),
            SemanticType::new(
                QuantityKind::Angle(AngleDomain::Electrical),
                ValueForm::Static,
            ),
            SemanticType::new(
                QuantityKind::AngularVelocity(AngleDomain::Mechanical),
                ValueForm::Static,
            ),
            SemanticType::new(
                QuantityKind::AngularVelocity(AngleDomain::Electrical),
                ValueForm::Static,
            ),
            SemanticType::new(QuantityKind::Torque, ValueForm::Static),
            SemanticType::new(QuantityKind::Energy, ValueForm::Static),
            SemanticType::new(QuantityKind::Pressure, ValueForm::Static),
            SemanticType::new(QuantityKind::Stress, ValueForm::Static),
            SemanticType::new(
                QuantityKind::Strain {
                    basis: StrainBasis::Tensor,
                    component: StrainComponent::Normal,
                },
                ValueForm::Static,
            ),
            SemanticType::new(
                QuantityKind::Strain {
                    basis: StrainBasis::Tensor,
                    component: StrainComponent::Shear,
                },
                ValueForm::Static,
            ),
            SemanticType::new(
                QuantityKind::Strain {
                    basis: StrainBasis::Engineering,
                    component: StrainComponent::Normal,
                },
                ValueForm::Static,
            ),
            SemanticType::new(
                QuantityKind::Strain {
                    basis: StrainBasis::Engineering,
                    component: StrainComponent::Shear,
                },
                ValueForm::Static,
            ),
            SemanticType::new(
                QuantityKind::Composition(CompositionBasis::MassFraction),
                ValueForm::Static,
            ),
            SemanticType::new(
                QuantityKind::Composition(CompositionBasis::MoleFraction),
                ValueForm::Static,
            ),
            SemanticType::new(
                QuantityKind::Composition(CompositionBasis::VolumeFraction),
                ValueForm::Static,
            ),
            SemanticType::new(QuantityKind::Mass, ValueForm::Static),
            SemanticType::new(QuantityKind::Amount, ValueForm::Static),
            SemanticType::new(QuantityKind::MolarMass, ValueForm::Static),
            SemanticType::new(QuantityKind::MassConcentration, ValueForm::Static),
            SemanticType::new(QuantityKind::AmountConcentration, ValueForm::Static),
            SemanticType::new(QuantityKind::Entropy, ValueForm::Static),
            SemanticType::new(QuantityKind::HeatCapacity, ValueForm::Static),
            SemanticType::new(QuantityKind::AcousticPressure, ValueForm::Static),
            SemanticType::new(QuantityKind::AcousticPower, ValueForm::Static),
            SemanticType::new(QuantityKind::Pressure, ValueForm::Instantaneous),
            SemanticType::new(QuantityKind::Pressure, ValueForm::Peak),
            SemanticType::new(QuantityKind::Pressure, ValueForm::Rms),
        ];
        let encodings: std::collections::BTreeSet<_> = semantic_types
            .into_iter()
            .map(|semantic_type| {
                ObjectiveSpec::from_validated_semantic(semantic_type).canonical_bytes()
            })
            .collect();
        assert_eq!(encodings.len(), semantic_types.len());
        assert!(!encodings.contains(&base.canonical_bytes()));
    }

    #[test]
    fn equal_bits_under_pressure_and_stress_schemas_move_identity() {
        let baseline_fixture = IdentityFixture::new();
        let baseline = report_identities(&baseline_fixture, &NORMAL_EXPECTATION_RULE);
        let pressure = ObjectiveSpec::from_validated_semantic(fs_qty::semantic::SemanticType::new(
            fs_qty::semantic::QuantityKind::Pressure,
            fs_qty::semantic::ValueForm::Static,
        ));
        let stress = ObjectiveSpec::from_validated_semantic(fs_qty::semantic::SemanticType::new(
            fs_qty::semantic::QuantityKind::Stress,
            fs_qty::semantic::ValueForm::Static,
        ));
        assert_eq!(pressure.dims(), stress.dims());
        assert_ne!(pressure, stress);

        let mut pressure_fixture = baseline_fixture.clone();
        for candidate in &mut pressure_fixture.candidates {
            candidate.objective_spec = pressure;
        }
        let mut stress_fixture = baseline_fixture;
        for candidate in &mut stress_fixture.candidates {
            candidate.objective_spec = stress;
        }
        let pressure_id = report_identities(&pressure_fixture, &NORMAL_EXPECTATION_RULE);
        let stress_id = report_identities(&stress_fixture, &NORMAL_EXPECTATION_RULE);
        assert_ne!(baseline, pressure_id);
        assert_ne!(pressure_id, stress_id);
    }

    #[test]
    fn exact_normal_expectation_rule_bits_are_identity_semantic() {
        let fixture = IdentityFixture::new();
        let baseline = report_identities(&fixture, &NORMAL_EXPECTATION_RULE);

        let mut changed_node = NORMAL_EXPECTATION_RULE;
        changed_node[0].0 = f64::from_bits(changed_node[0].0.to_bits() ^ 1);
        let node_identity = report_identities(&fixture, &changed_node);
        assert_ne!(baseline.0, node_identity.0);
        assert_ne!(baseline.1, node_identity.1);

        let mut changed_weight = NORMAL_EXPECTATION_RULE;
        changed_weight[0].1 = f64::from_bits(changed_weight[0].1.to_bits() ^ 1);
        let weight_identity = report_identities(&fixture, &changed_weight);
        assert_ne!(baseline.0, weight_identity.0);
        assert_ne!(baseline.1, weight_identity.1);
    }

    #[test]
    fn every_sealed_report_output_moves_both_identities() {
        let baseline_fixture = IdentityFixture::new();
        let baseline = report_identities(&baseline_fixture, &NORMAL_EXPECTATION_RULE);
        let mut mutations = Vec::new();

        let mut changed = baseline_fixture.clone();
        changed.placements[0] = "B".to_string();
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.sensors_placed = 2;
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.prior_total_variance = 0.33;
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.posterior_total_variance = 0.21;
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.variance_reduction = 0.376;
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.initial_evpi = 0.41;
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.final_evpi = 0.21;
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.decision_robust = true;
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.chosen_design = "B".to_string();
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.allocation[0].0 = "B".to_string();
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.allocation[0].1 = 0.11;
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.evpi_trace[0] = 0.41;
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.posteriors[0].name = "B".to_string();
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.posteriors[0].mean = 0.61;
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.posteriors[0].variance = 0.11;
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.assimilation_colors[0] = Color::Estimated {
            estimator: "sensor-B-v1".to_string(),
            dispersion: 0.01,
        };
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.assimilation_colors[0] = Color::Estimated {
            estimator: "sensor-A-v1".to_string(),
            dispersion: 0.02,
        };
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.variance_color_dispersion = 1.0e300;
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.evpi_color_dispersion = 0.21;
        mutations.push(changed);

        for (index, mutation) in mutations.iter().enumerate() {
            let identity = report_identities(mutation, &NORMAL_EXPECTATION_RULE);
            assert_ne!(
                baseline.0, identity.0,
                "variance identity ignored output mutation {index}"
            );
            assert_ne!(
                baseline.1, identity.1,
                "EVPI identity ignored output mutation {index}"
            );
        }

        for index in 0..8 {
            let mut changed = baseline_fixture.clone();
            match index {
                0 => changed.realized.completed_placements += 1,
                1 => changed.realized.action_rounds += 1,
                2 => changed.realized.positive_prior_candidates += 1,
                3 => changed.realized.setup_work_units += 1,
                4 => changed.realized.placement_work_units += 1,
                5 => changed.realized.incomplete_action_work_units += 1,
                6 => changed.realized.finalization_work_units += 1,
                7 => changed.realized.realized_work_units += 1,
                _ => unreachable!("retained realized-work field count"),
            }
            let identity = report_identities(&changed, &NORMAL_EXPECTATION_RULE);
            assert_ne!(
                baseline.0, identity.0,
                "variance identity ignored realized-work field {index}"
            );
            assert_ne!(
                baseline.1, identity.1,
                "EVPI identity ignored realized-work field {index}"
            );
        }
    }

    /// The retired clone-and-sort path, retained as the INDEPENDENT
    /// ORACLE for the canonical-menu evaluator (bead sj31i.62 G3,
    /// upgraded with sj31i.5 to the full multi-alternative loss): same
    /// fs_voi evaluator, same name sort, same finiteness/sign contract.
    fn oracle_checked_evpi(estimates: &[super::DesignEstimate]) -> Result<f64, super::OedError> {
        let mut canonical = estimates.to_vec();
        canonical.sort_by(|left, right| left.name.cmp(&right.name));
        let value = fs_voi::expected_opportunity_loss(&canonical);
        if value.is_finite() && value >= 0.0 {
            Ok(super::canonicalize_zero(value))
        } else {
            Err(super::OedError::NonFiniteComputation { quantity: "EVPI" })
        }
    }

    fn fixture_estimates() -> Vec<super::DesignEstimate> {
        // Includes an exact equal-mean tie (A/B) so the canonical
        // tie-break is exercised, plus a non-finite mean the scan must
        // skip exactly as the oracle does.
        [
            ("alpha", 0.60, 0.10),
            ("beta", 0.60, 0.12),
            ("gamma", 0.85, 0.06),
            ("delta", f64::NAN, 0.04),
            ("epsilon", 1.10, 0.02),
        ]
        .into_iter()
        .map(|(name, mean, statistical)| {
            super::DesignEstimate::new(
                name,
                mean,
                fs_voi::Uncertainty {
                    numerical: 0.0,
                    statistical,
                    model: 0.0,
                },
            )
        })
        .collect()
    }

    fn canonical_sorted(mut estimates: Vec<super::DesignEstimate>) -> Vec<super::DesignEstimate> {
        estimates.sort_by(|left, right| left.name.cmp(&right.name));
        estimates
    }

    /// sj31i.62 G0: canonical-order admission refuses unsorted and
    /// duplicate menus with the breaking position named; empty and
    /// singleton menus admit with EVPI exactly 0.
    #[test]
    fn canonical_menu_admission_g0() {
        let unsorted = fixture_estimates();
        let refusal = super::CanonicalDesignMenu::from_canonical(unsorted)
            .expect_err("declaration order is not canonical");
        // Declaration order is alpha, beta, gamma, delta, epsilon: the
        // first descending pair is (gamma, delta), and `position` names
        // the offending LATER element — delta at index 3.
        assert!(matches!(
            refusal,
            super::OedError::CanonicalOrderViolated { position: 3 }
        ));
        let mut duplicated = canonical_sorted(fixture_estimates());
        let duplicate_name = duplicated[0].name.clone();
        duplicated[1].name = duplicate_name;
        duplicated.sort_by(|left, right| left.name.cmp(&right.name));
        assert!(matches!(
            super::CanonicalDesignMenu::from_canonical(duplicated),
            Err(super::OedError::CanonicalOrderViolated { .. })
        ));
        let empty = super::CanonicalDesignMenu::from_canonical(Vec::new()).expect("empty menu");
        assert_eq!(
            empty
                .full_opportunity_loss_checked()
                .expect("empty full EOL"),
            0.0
        );
        let singleton =
            super::CanonicalDesignMenu::from_canonical(vec![fixture_estimates().remove(0)])
                .expect("singleton menu");
        assert_eq!(
            singleton
                .full_opportunity_loss_checked()
                .expect("singleton full EOL"),
            0.0
        );
    }

    /// sj31i.62 G3: the no-sort canonical evaluator is BITWISE equal to
    /// the retired clone-and-sort oracle, on canonical fixtures and
    /// under input permutations (which the oracle absorbs by sorting
    /// and the menu absorbs by canonical admission).
    #[test]
    fn canonical_menu_matches_oracle_bitwise() {
        let base = fixture_estimates();
        let oracle = oracle_checked_evpi(&base).expect("oracle EVPI");
        // Deterministic permutations: rotations of the declaration order.
        for rotation in 0..base.len() {
            let mut permuted = base.clone();
            permuted.rotate_left(rotation);
            assert_eq!(
                oracle_checked_evpi(&permuted)
                    .expect("oracle is order-independent")
                    .to_bits(),
                oracle.to_bits(),
            );
            let menu = super::CanonicalDesignMenu::from_canonical(canonical_sorted(permuted))
                .expect("canonical menu");
            assert_eq!(
                menu.full_opportunity_loss_checked()
                    .expect("menu EVPI")
                    .to_bits(),
                oracle.to_bits(),
                "no-sort evaluator must match the oracle bitwise (rotation {rotation})"
            );
        }
    }

    /// sj31i.62 G0+G3: the one-index override view validates its index
    /// and payload, cannot mutate the menu, and its EVPI is bitwise
    /// equal to the oracle run on an explicitly rebuilt overridden menu.
    #[test]
    fn override_view_matches_rebuilt_menu_bitwise() {
        let canonical = canonical_sorted(fixture_estimates());
        let menu = super::CanonicalDesignMenu::from_canonical(canonical.clone()).expect("menu");
        assert!(matches!(
            super::MeanOverrideView::new(&menu, canonical.len(), 0.5, 0.1),
            Err(super::OedError::OverrideInvalid { .. })
        ));
        assert!(matches!(
            super::MeanOverrideView::new(&menu, 0, f64::NAN, 0.1),
            Err(super::OedError::OverrideInvalid { .. })
        ));
        assert!(matches!(
            super::MeanOverrideView::new(&menu, 0, 0.5, f64::NAN),
            Err(super::OedError::OverrideInvalid { .. })
        ));
        assert!(matches!(
            super::MeanOverrideView::new(&menu, 0, 0.5, -0.1),
            Err(super::OedError::OverrideInvalid { .. })
        ));
        for index in 0..canonical.len() {
            for (mean, statistical) in [(0.55, 0.05), (0.60, 0.0), (2.0, 0.3)] {
                let view = super::MeanOverrideView::new(&menu, index, mean, statistical)
                    .expect("valid override view");
                let mut rebuilt = canonical.clone();
                rebuilt[index].mean = mean;
                rebuilt[index].uncertainty.statistical = statistical;
                assert_eq!(
                    view.opportunity_loss_checked()
                        .expect("view EVPI")
                        .to_bits(),
                    oracle_checked_evpi(&rebuilt)
                        .expect("oracle EVPI")
                        .to_bits(),
                    "override view must equal an independently rebuilt menu \
                     (index {index}, mean {mean}, stat {statistical})"
                );
            }
        }
        // The menu is observably unchanged after every view: identity
        // order and values are exactly the admitted ones. The fixture
        // deliberately carries a NaN mean (delta), so slice PartialEq is
        // unusable (NaN != NaN); the bitwise intent needs exact bit
        // patterns per field.
        assert_eq!(menu.estimates().len(), canonical.len());
        for (kept, admitted) in menu.estimates().iter().zip(&canonical) {
            assert_eq!(kept.name, admitted.name);
            assert_eq!(kept.mean.to_bits(), admitted.mean.to_bits());
            assert_eq!(
                kept.uncertainty.numerical.to_bits(),
                admitted.uncertainty.numerical.to_bits()
            );
            assert_eq!(
                kept.uncertainty.statistical.to_bits(),
                admitted.uncertainty.statistical.to_bits()
            );
            assert_eq!(
                kept.uncertainty.model.to_bits(),
                admitted.uncertainty.model.to_bits()
            );
        }
    }

    #[test]
    fn byte_plan_admission_arithmetic_is_checked() {
        // A shape whose charge formulas overflow u128 refuses typed at
        // admission instead of wrapping into a fictitious envelope.
        assert!(matches!(
            CampaignBytePlan::checked(usize::MAX, usize::MAX, usize::MAX, usize::MAX),
            Err(super::OedError::BytePlanOverflow {
                candidates: usize::MAX,
                max_sensors: usize::MAX,
            })
        ));
        // Real campaign bounds always admit: the worst legal shape.
        let plan = CampaignBytePlan::checked(
            super::MAX_CAMPAIGN_CANDIDATES,
            super::MAX_CAMPAIGN_SENSORS,
            super::MAX_CAMPAIGN_CANDIDATES * super::MAX_CANDIDATE_NAME_BYTES,
            super::MAX_CANDIDATE_NAME_BYTES,
        )
        .expect("the maximum legal campaign shape has a checked byte envelope");
        assert!(plan.admitted_byte_units > 0);
    }

    #[test]
    fn byte_charges_refuse_typed_beyond_the_admitted_plan() {
        with_test_cx(|cx| {
            let plan = CampaignBytePlan::checked(2, 1, 2, 1).expect("tiny byte plan");
            let mut progress =
                CampaignProgress::admit(cx, 0, 0, plan).expect("tiny fixture admits");
            let opening = progress.charged_byte_units;
            assert_eq!(opening, plan.admission_byte_units);
            // A charge past the admitted envelope refuses typed and
            // leaves the ledger untouched — no partial accounting.
            let refused = progress.charge_bytes("adversarial seam", u128::MAX);
            assert!(matches!(
                refused,
                Err(super::OedError::ByteBudgetExceeded {
                    at: "adversarial seam",
                    ..
                })
            ));
            assert_eq!(progress.charged_byte_units, opening);
            assert_eq!(progress.retained_byte_units, 0);
            // Retained charges land in both ledgers.
            progress
                .retain_bytes("retained seam", 8)
                .expect("a covered retained charge is accepted");
            assert_eq!(progress.charged_byte_units, opening + 8);
            assert_eq!(progress.retained_byte_units, 8);
        });
    }
}
