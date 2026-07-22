//! Versioned eight-term engineering uncertainty budgets.
//!
//! The legacy [`crate::UncertaintyBreakdown`] is intentionally retained as a
//! decision-facing three-slice view. This module records the richer source
//! accounting needed to answer which engineering evidence should be bought
//! next without forcing every term into an additive Gaussian model.
//!
//! Correlation is honored only when an explicit, validated covariance block
//! names every participating term. Everything else composes by conservative
//! linear addition. An unknown term never becomes zero: it changes the total
//! into [`BudgetTotal::Unknown`] while retaining the known contribution.

use core::fmt;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

use fs_blake3::{ContentHash, hash_domain};

use crate::{UncertaintyBreakdown, balance::BoundedId};

/// Current canonical transport version.
pub const ENGINEERING_UNCERTAINTY_SCHEMA_VERSION: u32 = 1;
/// Number of mandatory engineering uncertainty sources.
pub const ENGINEERING_UNCERTAINTY_TERM_COUNT: usize = 8;
/// Maximum retained explanation or justification bytes.
pub const MAX_UNCERTAINTY_TEXT_BYTES: usize = 1024;
/// Maximum accepted canonical transport size.
pub const MAX_UNCERTAINTY_CANONICAL_BYTES: usize = 1024 * 1024;

const MAGIC: &[u8; 4] = b"FSEU";
const IDENTITY_DOMAIN: &str = "org.frankensim.fs-evidence.engineering-uncertainty.v1";
const BLOCK_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-evidence.engineering-uncertainty.covariance-block.v1";
const COMPOSITION_DOMAIN: &str =
    "org.frankensim.fs-evidence.engineering-uncertainty.composition.v1";
const DOMINANCE_ORDER: [EngineeringUncertaintyKind; ENGINEERING_UNCERTAINTY_TERM_COUNT] = [
    EngineeringUncertaintyKind::ModelForm,
    EngineeringUncertaintyKind::Measurement,
    EngineeringUncertaintyKind::BoundaryConditions,
    EngineeringUncertaintyKind::Parameters,
    EngineeringUncertaintyKind::Geometry,
    EngineeringUncertaintyKind::Discretization,
    EngineeringUncertaintyKind::SolverAlgebraic,
    EngineeringUncertaintyKind::Roundoff,
];

/// Stable named rule violated by budget admission.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum UncertaintyRule {
    /// A bounded identifier was malformed.
    IdBounds,
    /// A required explanation or justification was empty or oversized.
    TextBounds,
    /// A numeric value was non-finite, negative, or otherwise inadmissible.
    NumericDomain,
    /// The eight-term source set was incomplete or duplicated.
    TermSet,
    /// A covariance block did not match the terms that reference it.
    CovarianceMembership,
    /// A covariance matrix was malformed, asymmetric, or not positive semidefinite.
    CovarianceMatrix,
    /// Two budgets described different QoIs or units.
    IncompatibleBudget,
    /// A numerical-only update named a non-numerical source or omitted one of
    /// the three numerical sources.
    NumericalUpdate,
    /// A requirement, plausibility bound, or flip-analysis input was invalid.
    RequirementAssessment,
    /// A collection or transport exceeded its declared bound.
    CollectionBudget,
    /// The transport schema version is unknown.
    SchemaVersion,
}

impl UncertaintyRule {
    /// Stable machine-readable rule slug.
    #[must_use]
    pub const fn slug(self) -> &'static str {
        match self {
            Self::IdBounds => "uncertainty-id-bounds",
            Self::TextBounds => "uncertainty-text-bounds",
            Self::NumericDomain => "uncertainty-numeric-domain",
            Self::TermSet => "uncertainty-term-set",
            Self::CovarianceMembership => "uncertainty-covariance-membership",
            Self::CovarianceMatrix => "uncertainty-covariance-matrix",
            Self::IncompatibleBudget => "uncertainty-incompatible-budget",
            Self::NumericalUpdate => "uncertainty-numerical-update",
            Self::RequirementAssessment => "uncertainty-requirement-assessment",
            Self::CollectionBudget => "uncertainty-collection-budget",
            Self::SchemaVersion => "uncertainty-schema-version",
        }
    }
}

/// Typed admission refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UncertaintyError {
    rule: UncertaintyRule,
    detail: String,
}

impl UncertaintyError {
    fn new(rule: UncertaintyRule, detail: impl Into<String>) -> Self {
        Self {
            rule,
            detail: detail.into(),
        }
    }

    /// Violated rule.
    #[must_use]
    pub const fn rule(&self) -> UncertaintyRule {
        self.rule
    }

    /// Human-readable diagnostic detail.
    #[must_use]
    pub fn detail(&self) -> &str {
        &self.detail
    }
}

impl fmt::Display for UncertaintyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.rule.slug(), self.detail)
    }
}

impl std::error::Error for UncertaintyError {}

fn refuse<T>(rule: UncertaintyRule, detail: impl Into<String>) -> Result<T, UncertaintyError> {
    Err(UncertaintyError::new(rule, detail))
}

/// Canonical transport refusal with a byte offset.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UncertaintyCodecError {
    offset: usize,
    detail: String,
}

impl UncertaintyCodecError {
    fn at(offset: usize, detail: impl Into<String>) -> Self {
        Self {
            offset,
            detail: detail.into(),
        }
    }

    /// Byte offset at which decoding refused.
    #[must_use]
    pub const fn offset(&self) -> usize {
        self.offset
    }

    /// Human-readable diagnostic detail.
    #[must_use]
    pub fn detail(&self) -> &str {
        &self.detail
    }
}

impl fmt::Display for UncertaintyCodecError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "uncertainty-codec@{}: {}", self.offset, self.detail)
    }
}

impl std::error::Error for UncertaintyCodecError {}

/// The mandatory eight engineering uncertainty sources.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum EngineeringUncertaintyKind {
    /// Floating-point and arithmetic rounding.
    Roundoff,
    /// Algebraic and iterative/nonlinear solver termination.
    SolverAlgebraic,
    /// Spatial and temporal discretization.
    Discretization,
    /// Geometry representation, import, registration, and as-built state.
    Geometry,
    /// Material, manufacturing, and other parameter uncertainty.
    Parameters,
    /// Boundary and operating-condition uncertainty.
    BoundaryConditions,
    /// Model discrepancy and closure inadequacy.
    ModelForm,
    /// Sensor, preprocessing, and comparison-data uncertainty.
    Measurement,
}

impl EngineeringUncertaintyKind {
    /// Canonical wire and presentation order.
    pub const ALL: [Self; ENGINEERING_UNCERTAINTY_TERM_COUNT] = [
        Self::Roundoff,
        Self::SolverAlgebraic,
        Self::Discretization,
        Self::Geometry,
        Self::Parameters,
        Self::BoundaryConditions,
        Self::ModelForm,
        Self::Measurement,
    ];

    /// Stable source name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Roundoff => "roundoff",
            Self::SolverAlgebraic => "solver-algebraic",
            Self::Discretization => "discretization",
            Self::Geometry => "geometry",
            Self::Parameters => "parameters",
            Self::BoundaryConditions => "boundary-conditions",
            Self::ModelForm => "model-form",
            Self::Measurement => "measurement",
        }
    }

    const fn code(self) -> u8 {
        match self {
            Self::Roundoff => 1,
            Self::SolverAlgebraic => 2,
            Self::Discretization => 3,
            Self::Geometry => 4,
            Self::Parameters => 5,
            Self::BoundaryConditions => 6,
            Self::ModelForm => 7,
            Self::Measurement => 8,
        }
    }

    fn from_code(code: u8) -> Option<Self> {
        Self::ALL.into_iter().find(|kind| kind.code() == code)
    }

    const fn legacy_slice(self) -> LegacySlice {
        match self {
            Self::Roundoff | Self::SolverAlgebraic | Self::Discretization | Self::Geometry => {
                LegacySlice::Numerical
            }
            Self::Parameters | Self::BoundaryConditions | Self::Measurement => {
                LegacySlice::Statistical
            }
            Self::ModelForm => LegacySlice::Model,
        }
    }
}

/// A content-addressed artifact with a bounded semantic role.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct UncertaintyArtifactRef {
    role: BoundedId,
    digest: ContentHash,
}

impl UncertaintyArtifactRef {
    /// Admit an artifact reference.
    pub fn new(role: &str, digest: ContentHash) -> Result<Self, UncertaintyError> {
        let role = BoundedId::new(role)
            .map_err(|error| UncertaintyError::new(UncertaintyRule::IdBounds, error.to_string()))?;
        Ok(Self { role, digest })
    }

    /// Semantic role of the referenced artifact.
    #[must_use]
    pub fn role(&self) -> &str {
        self.role.as_str()
    }

    /// Content identity of the referenced artifact.
    #[must_use]
    pub const fn digest(&self) -> ContentHash {
        self.digest
    }
}

/// Retained distribution summary plus a replayable artifact reference.
#[derive(Debug, Clone, PartialEq)]
pub struct DistributionTerm {
    /// Distribution mean in the budget unit.
    pub mean: f64,
    /// Standard deviation in the budget unit.
    pub standard_deviation: f64,
    /// Conservative half-width used by budget aggregation.
    pub conservative_half_width: f64,
    /// Coverage/confidence level for the stated half-width.
    pub level: f64,
    /// Replayable samples/fit/posterior artifact.
    pub replay: UncertaintyArtifactRef,
}

impl DistributionTerm {
    fn validate(&self) -> Result<(), UncertaintyError> {
        if !self.mean.is_finite()
            || !non_negative_finite(self.standard_deviation)
            || !non_negative_finite(self.conservative_half_width)
            || !self.level.is_finite()
            || self.level <= 0.0
            || self.level >= 1.0
        {
            return refuse(
                UncertaintyRule::NumericDomain,
                "distribution mean, deviation, half-width, or level is outside its domain",
            );
        }
        Ok(())
    }
}

/// Retained ensemble envelope plus its replayable member artifact.
#[derive(Debug, Clone, PartialEq)]
pub struct EnsembleTerm {
    /// Number of retained members.
    pub member_count: u32,
    /// Conservative half-width of the ensemble envelope.
    pub conservative_half_width: f64,
    /// Replayable member/output artifact.
    pub replay: UncertaintyArtifactRef,
}

impl EnsembleTerm {
    fn validate(&self) -> Result<(), UncertaintyError> {
        if self.member_count == 0 || !non_negative_finite(self.conservative_half_width) {
            return refuse(
                UncertaintyRule::NumericDomain,
                "ensemble needs a positive member count and finite non-negative half-width",
            );
        }
        Ok(())
    }
}

/// Explicit covariance block spanning two or more engineering sources.
#[derive(Debug, Clone, PartialEq)]
pub struct CovarianceBlock {
    id: BoundedId,
    covariance_artifact: UncertaintyArtifactRef,
    members: Vec<EngineeringUncertaintyKind>,
    covariance: Vec<f64>,
}

impl CovarianceBlock {
    /// Admit a covariance matrix in row-major member order.
    pub fn try_new(
        id: &str,
        covariance_artifact: UncertaintyArtifactRef,
        members: Vec<EngineeringUncertaintyKind>,
        covariance: Vec<f64>,
    ) -> Result<Self, UncertaintyError> {
        let id = BoundedId::new(id)
            .map_err(|error| UncertaintyError::new(UncertaintyRule::IdBounds, error.to_string()))?;
        if !(2..=ENGINEERING_UNCERTAINTY_TERM_COUNT).contains(&members.len()) {
            return refuse(
                UncertaintyRule::CovarianceMembership,
                "covariance block must span 2..=8 terms",
            );
        }
        if !members.windows(2).all(|pair| pair[0] < pair[1]) {
            return refuse(
                UncertaintyRule::CovarianceMembership,
                "covariance members must be unique and in canonical source order",
            );
        }
        let expected = members.len().checked_mul(members.len()).ok_or_else(|| {
            UncertaintyError::new(
                UncertaintyRule::CollectionBudget,
                "covariance dimensions overflow",
            )
        })?;
        if covariance.len() != expected {
            return refuse(
                UncertaintyRule::CovarianceMatrix,
                format!(
                    "covariance has {} entries; expected {expected}",
                    covariance.len()
                ),
            );
        }
        validate_covariance(&covariance, members.len())?;
        Ok(Self {
            id,
            covariance_artifact,
            members,
            covariance,
        })
    }

    /// Stable block label.
    #[must_use]
    pub fn id(&self) -> &str {
        self.id.as_str()
    }

    /// Artifact holding the covariance authority.
    #[must_use]
    pub const fn covariance_artifact(&self) -> &UncertaintyArtifactRef {
        &self.covariance_artifact
    }

    /// Canonically ordered participating sources.
    #[must_use]
    pub fn members(&self) -> &[EngineeringUncertaintyKind] {
        &self.members
    }

    /// Row-major covariance values.
    #[must_use]
    pub fn covariance(&self) -> &[f64] {
        &self.covariance
    }

    /// Content identity of the exact block declaration.
    #[must_use]
    pub fn content_id(&self) -> ContentHash {
        let mut bytes = Vec::new();
        encode_block(&mut bytes, self);
        hash_domain(BLOCK_IDENTITY_DOMAIN, &bytes)
    }

    fn marginal_half_width(&self, kind: EngineeringUncertaintyKind) -> Option<f64> {
        let index = self.members.iter().position(|member| *member == kind)?;
        Some(next_up(
            self.covariance[index * self.members.len() + index].sqrt(),
        ))
    }

    fn combined_half_width(&self) -> f64 {
        let mut variance = 0.0;
        for value in &self.covariance {
            variance = add_up(variance, *value);
        }
        next_up(variance.max(0.0).sqrt())
    }
}

/// Representation of one uncertainty contribution.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum TermValue {
    /// Certified non-negative interval for a half-width.
    IntervalBound {
        /// Known lower bound on the half-width.
        lower: f64,
        /// Conservative upper bound used for aggregation.
        upper: f64,
    },
    /// Distribution summary with a replayable authority.
    Distribution(DistributionTerm),
    /// Ensemble envelope with a replayable authority.
    Ensemble(EnsembleTerm),
    /// Explicit covariance block shared by every participating term.
    CorrelatedBlock(CovarianceBlock),
    /// No finite upper bound is currently justified.
    Unknown {
        /// Named evidence gap; never empty.
        reason: String,
    },
    /// Explicit zero contribution with a non-empty justification.
    Negligible {
        /// Why zero is admissible for this source.
        justification: String,
    },
}

impl TermValue {
    /// Construct an interval-bound term.
    pub fn interval(lower: f64, upper: f64) -> Result<Self, UncertaintyError> {
        let value = Self::IntervalBound { lower, upper };
        value.validate()?;
        Ok(value)
    }

    /// Construct an explicitly unknown term.
    pub fn unknown(reason: impl Into<String>) -> Result<Self, UncertaintyError> {
        let value = Self::Unknown {
            reason: reason.into(),
        };
        value.validate()?;
        Ok(value)
    }

    /// Construct an explicitly negligible term.
    pub fn negligible(justification: impl Into<String>) -> Result<Self, UncertaintyError> {
        let value = Self::Negligible {
            justification: justification.into(),
        };
        value.validate()?;
        Ok(value)
    }

    fn validate(&self) -> Result<(), UncertaintyError> {
        match self {
            Self::IntervalBound { lower, upper } => {
                if !non_negative_finite(*lower) || !non_negative_finite(*upper) || lower > upper {
                    return refuse(
                        UncertaintyRule::NumericDomain,
                        "interval half-width bounds must be finite, non-negative, and ordered",
                    );
                }
            }
            Self::Distribution(summary) => summary.validate()?,
            Self::Ensemble(summary) => summary.validate()?,
            Self::CorrelatedBlock(block) => {
                validate_covariance(&block.covariance, block.members.len())?;
            }
            Self::Unknown { reason } => admit_text("unknown reason", reason)?,
            Self::Negligible { justification } => {
                admit_text("negligible justification", justification)?;
            }
        }
        Ok(())
    }

    fn marginal_half_width(&self, kind: EngineeringUncertaintyKind) -> Option<f64> {
        match self {
            Self::IntervalBound { upper, .. } => Some(*upper),
            Self::Distribution(summary) => Some(summary.conservative_half_width),
            Self::Ensemble(summary) => Some(summary.conservative_half_width),
            Self::CorrelatedBlock(block) => block.marginal_half_width(kind),
            Self::Unknown { .. } => None,
            Self::Negligible { .. } => Some(0.0),
        }
    }
}

/// One mandatory source term and its provenance authority.
#[derive(Debug, Clone, PartialEq)]
pub struct EngineeringUncertaintyTerm {
    /// Engineering source category.
    kind: EngineeringUncertaintyKind,
    /// Rich uncertainty representation.
    value: TermValue,
    /// Evidence/provenance artifact supporting this declaration.
    provenance: UncertaintyArtifactRef,
}

impl EngineeringUncertaintyTerm {
    /// Admit one term.
    pub fn try_new(
        kind: EngineeringUncertaintyKind,
        value: TermValue,
        provenance: UncertaintyArtifactRef,
    ) -> Result<Self, UncertaintyError> {
        value.validate()?;
        if let TermValue::CorrelatedBlock(block) = &value
            && !block.members.contains(&kind)
        {
            return refuse(
                UncertaintyRule::CovarianceMembership,
                format!("block {} does not include {}", block.id(), kind.name()),
            );
        }
        Ok(Self {
            kind,
            value,
            provenance,
        })
    }

    /// Engineering source category fixed at admission.
    #[must_use]
    pub const fn kind(&self) -> EngineeringUncertaintyKind {
        self.kind
    }

    /// Admitted representation. No mutable accessor is exposed: changing a
    /// value must re-enter [`Self::try_new`].
    #[must_use]
    pub const fn value(&self) -> &TermValue {
        &self.value
    }

    /// Provenance authority fixed together with the source and value.
    #[must_use]
    pub const fn provenance(&self) -> &UncertaintyArtifactRef {
        &self.provenance
    }
}

/// Sealed update containing exactly the numerical certificate sources.
///
/// Model-form, measurement, geometry, parameter, and boundary-condition
/// evidence cannot be represented by this type. Applying one therefore
/// changes only roundoff, solver/algebraic, and discretization terms; callers
/// that acquired better numerical evidence cannot accidentally rewrite a
/// model or experimental authority term in the same operation.
#[derive(Debug, Clone, PartialEq)]
pub struct NumericalUncertaintyUpdate {
    terms: [EngineeringUncertaintyTerm; 3],
}

impl NumericalUncertaintyUpdate {
    /// Admit exactly one roundoff, solver/algebraic, and discretization term.
    pub fn try_new(terms: Vec<EngineeringUncertaintyTerm>) -> Result<Self, UncertaintyError> {
        const KINDS: [EngineeringUncertaintyKind; 3] = [
            EngineeringUncertaintyKind::Roundoff,
            EngineeringUncertaintyKind::SolverAlgebraic,
            EngineeringUncertaintyKind::Discretization,
        ];
        if terms.len() != KINDS.len() {
            return refuse(
                UncertaintyRule::NumericalUpdate,
                format!("numerical update has {} terms; expected 3", terms.len()),
            );
        }
        let mut by_kind = BTreeMap::new();
        for term in terms {
            if !KINDS.contains(&term.kind) {
                return refuse(
                    UncertaintyRule::NumericalUpdate,
                    format!("{} is not a numerical update source", term.kind.name()),
                );
            }
            let kind = term.kind;
            if by_kind.insert(kind, term).is_some() {
                return refuse(
                    UncertaintyRule::NumericalUpdate,
                    format!("duplicate {} numerical update term", kind.name()),
                );
            }
        }
        let terms = KINDS.map(|kind| {
            by_kind.remove(&kind).ok_or_else(|| {
                UncertaintyError::new(
                    UncertaintyRule::NumericalUpdate,
                    format!("missing {} numerical update term", kind.name()),
                )
            })
        });
        let terms = terms.into_iter().collect::<Result<Vec<_>, _>>()?;
        let terms = terms.try_into().map_err(|_| {
            UncertaintyError::new(
                UncertaintyRule::NumericalUpdate,
                "numerical update term count changed",
            )
        })?;
        Ok(Self { terms })
    }

    /// Canonically ordered numerical terms.
    #[must_use]
    pub const fn terms(&self) -> &[EngineeringUncertaintyTerm; 3] {
        &self.terms
    }
}

/// Derived total without laundering unknown terms into zero.
#[derive(Debug, Clone, PartialEq)]
pub enum BudgetTotal {
    /// Every source has a finite conservative contribution.
    Bounded {
        /// Conservative half-width in the budget unit.
        conservative_half_width: f64,
    },
    /// At least one source has no justified finite upper bound.
    Unknown {
        /// Conservative contribution from the fully specified sources.
        known_conservative_half_width: f64,
        /// Canonically ordered unresolved sources.
        unknown_terms: Vec<EngineeringUncertaintyKind>,
    },
    /// Fully specified finite inputs overflowed the representable aggregate.
    Unbounded {
        /// Stable diagnostic; no finite total may be inferred.
        reason: &'static str,
    },
}

/// Generalized dominant-source result.
#[derive(Debug, Clone, PartialEq)]
pub enum DominantEngineeringTerm {
    /// One finite source has the largest marginal contribution.
    Known {
        /// Dominant source.
        kind: EngineeringUncertaintyKind,
        /// Conservative marginal half-width.
        conservative_half_width: f64,
    },
    /// A finite dominance result is impossible while sources are unknown.
    Unknown {
        /// Canonically ordered unresolved sources.
        terms: Vec<EngineeringUncertaintyKind>,
    },
}

/// Direction of a scalar engineering requirement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RequirementRelation {
    /// The quantity must remain below the stated limit.
    AtMost,
    /// The quantity must remain above the stated limit.
    AtLeast,
}

/// One sourced scalar requirement evaluated against an uncertainty budget.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ScalarRequirement {
    id: BoundedId,
    qoi: BoundedId,
    unit: BoundedId,
    relation: RequirementRelation,
    limit_bits: u64,
    provenance: UncertaintyArtifactRef,
}

impl ScalarRequirement {
    /// Admit a finite scalar requirement with an explicit source artifact.
    pub fn try_new(
        id: &str,
        qoi: &str,
        unit: &str,
        relation: RequirementRelation,
        limit: f64,
        provenance: UncertaintyArtifactRef,
    ) -> Result<Self, UncertaintyError> {
        if !limit.is_finite() {
            return refuse(
                UncertaintyRule::RequirementAssessment,
                "requirement limit must be finite",
            );
        }
        let admit_id = |label: &'static str, value: &str| {
            BoundedId::new(value).map_err(|error| {
                UncertaintyError::new(
                    UncertaintyRule::RequirementAssessment,
                    format!("invalid requirement {label}: {error}"),
                )
            })
        };
        Ok(Self {
            id: admit_id("id", id)?,
            qoi: admit_id("qoi", qoi)?,
            unit: admit_id("unit", unit)?,
            relation,
            limit_bits: limit.to_bits(),
            provenance,
        })
    }

    /// Stable requirement identity.
    #[must_use]
    pub fn id(&self) -> &str {
        self.id.as_str()
    }

    /// Quantity of interest governed by this requirement.
    #[must_use]
    pub fn qoi(&self) -> &str {
        self.qoi.as_str()
    }

    /// Unit shared with the governed uncertainty budget.
    #[must_use]
    pub fn unit(&self) -> &str {
        self.unit.as_str()
    }

    /// Requirement direction.
    #[must_use]
    pub const fn relation(&self) -> RequirementRelation {
        self.relation
    }

    /// Finite requirement limit.
    #[must_use]
    pub fn limit(&self) -> f64 {
        f64::from_bits(self.limit_bits)
    }

    /// Source artifact for the exact requirement.
    #[must_use]
    pub const fn provenance(&self) -> &UncertaintyArtifactRef {
        &self.provenance
    }
}

/// A sourced symmetric plausibility bound for an otherwise unknown term.
///
/// The budget term remains [`TermValue::Unknown`]: this record says only that
/// a decision audit may bound its absolute effect by `maximum_abs_effect`.
/// It does not promote the term to an interval certificate.
#[derive(Debug, Clone, PartialEq)]
pub struct UnknownPlausibilityBound {
    kind: EngineeringUncertaintyKind,
    maximum_abs_effect: f64,
    provenance: UncertaintyArtifactRef,
}

impl UnknownPlausibilityBound {
    /// Admit a finite non-negative plausibility half-width and its authority.
    pub fn try_new(
        kind: EngineeringUncertaintyKind,
        maximum_abs_effect: f64,
        provenance: UncertaintyArtifactRef,
    ) -> Result<Self, UncertaintyError> {
        if !non_negative_finite(maximum_abs_effect) {
            return refuse(
                UncertaintyRule::RequirementAssessment,
                format!(
                    "{} plausibility bound must be finite, non-negative, and not negative zero",
                    kind.name()
                ),
            );
        }
        Ok(Self {
            kind,
            maximum_abs_effect,
            provenance,
        })
    }

    /// Unknown source constrained by this declaration.
    #[must_use]
    pub const fn kind(&self) -> EngineeringUncertaintyKind {
        self.kind
    }

    /// Maximum absolute effect in the budget unit.
    #[must_use]
    pub const fn maximum_abs_effect(&self) -> f64 {
        self.maximum_abs_effect
    }

    /// Authority for the plausibility declaration.
    #[must_use]
    pub const fn provenance(&self) -> &UncertaintyArtifactRef {
        &self.provenance
    }
}

/// Whether one flipping source had a sourced finite plausibility bound.
#[derive(Debug, Clone, PartialEq)]
pub enum FlipBound {
    /// No finite plausibility authority was supplied.
    Unbounded,
    /// A sourced symmetric absolute-effect ceiling was supplied.
    Bounded(UnknownPlausibilityBound),
}

/// One unknown source that can participate in changing the verdict.
#[derive(Debug, Clone, PartialEq)]
pub struct FlippingUnknown {
    kind: EngineeringUncertaintyKind,
    reason: String,
    required_magnitude: f64,
    bound: FlipBound,
    suggested_action: crate::action::ActionKind,
}

impl FlippingUnknown {
    /// Uncertainty source that can change the verdict.
    #[must_use]
    pub const fn kind(&self) -> EngineeringUncertaintyKind {
        self.kind
    }

    /// Original named evidence gap from the budget.
    #[must_use]
    pub fn reason(&self) -> &str {
        &self.reason
    }

    /// Smallest additional absolute effect needed from this source when every
    /// other bounded unknown is allowed to move adversarially.
    #[must_use]
    pub const fn required_magnitude(&self) -> f64 {
        self.required_magnitude
    }

    /// Finite plausibility authority, or an explicit unbounded state.
    #[must_use]
    pub const fn bound(&self) -> &FlipBound {
        &self.bound
    }

    /// Default evidence-action class for resolving this source. Cost-aware
    /// ranking remains an L4 `fs-voi` responsibility.
    #[must_use]
    pub const fn suggested_action(&self) -> crate::action::ActionKind {
        self.suggested_action
    }
}

/// Tri-state requirement result. An indeterminate result is a first-class
/// outcome and never aliases compliance or non-compliance.
#[derive(Debug, Clone, PartialEq)]
pub enum ComplianceVerdict {
    /// The complete admitted effect envelope remains strictly inside the
    /// requirement, with a positive residual margin.
    Compliant {
        /// Conservative distance from the worst admitted case to the limit.
        margin: f64,
        /// Identity of the exact eight-term budget.
        budget: ContentHash,
        /// Exact sourced requirement evaluated.
        requirement: ScalarRequirement,
    },
    /// The complete admitted effect envelope remains strictly outside the
    /// requirement, with a positive residual shortfall.
    NonCompliant {
        /// Conservative distance from the best admitted case to the limit.
        shortfall: f64,
        /// Identity of the exact eight-term budget.
        budget: ContentHash,
        /// Exact sourced requirement evaluated.
        requirement: ScalarRequirement,
    },
    /// The known band touches/straddles the limit or one or more unknowns can
    /// change the baseline verdict.
    Indeterminate {
        /// Outward-rounded lower endpoint from fully specified terms only.
        known_lower: f64,
        /// Outward-rounded upper endpoint from fully specified terms only.
        known_upper: f64,
        /// Unknown sources that can participate in a verdict flip.
        flipping_unknowns: Vec<FlippingUnknown>,
        /// Identity of the exact eight-term budget.
        budget: ContentHash,
        /// Exact sourced requirement evaluated.
        requirement: ScalarRequirement,
    },
}

impl ComplianceVerdict {
    /// Exact sourced requirement consumed by this verdict.
    #[must_use]
    pub const fn requirement(&self) -> &ScalarRequirement {
        match self {
            Self::Compliant { requirement, .. }
            | Self::NonCompliant { requirement, .. }
            | Self::Indeterminate { requirement, .. } => requirement,
        }
    }

    /// Identity of the exact uncertainty budget consumed by this verdict.
    #[must_use]
    pub const fn budget(&self) -> ContentHash {
        match self {
            Self::Compliant { budget, .. }
            | Self::NonCompliant { budget, .. }
            | Self::Indeterminate { budget, .. } => *budget,
        }
    }

    /// Unknown sources whose admitted effects can change the verdict.
    /// Binary verdicts have no flipping unknowns.
    #[must_use]
    pub fn flipping_unknowns(&self) -> &[FlippingUnknown] {
        match self {
            Self::Indeterminate {
                flipping_unknowns, ..
            } => flipping_unknowns,
            Self::Compliant { .. } | Self::NonCompliant { .. } => &[],
        }
    }

    /// Deterministic reviewer-facing audit trail. This rendering is an
    /// explanation, not a separate scientific authority.
    #[must_use]
    pub fn render_report(&self) -> String {
        let requirement = self.requirement();
        let relation = match requirement.relation() {
            RequirementRelation::AtMost => "at-most",
            RequirementRelation::AtLeast => "at-least",
        };
        let mut output = format!(
            "requirement={} qoi={} unit={} relation={} limit={} requirement-provenance={}@{} budget={}\n",
            requirement.id(),
            requirement.qoi(),
            requirement.unit(),
            relation,
            requirement.limit(),
            requirement.provenance().role(),
            requirement.provenance().digest(),
            self.budget()
        );
        match self {
            Self::Compliant { margin, .. } => {
                let _ = writeln!(output, "verdict=compliant residual-margin={margin}");
            }
            Self::NonCompliant { shortfall, .. } => {
                let _ = writeln!(
                    output,
                    "verdict=non-compliant residual-shortfall={shortfall}"
                );
            }
            Self::Indeterminate {
                known_lower,
                known_upper,
                flipping_unknowns,
                ..
            } => {
                let _ = writeln!(
                    output,
                    "verdict=indeterminate known-band=[{known_lower},{known_upper}] flipping-unknowns={}",
                    flipping_unknowns.len()
                );
                for unknown in flipping_unknowns {
                    let bound = match unknown.bound() {
                        FlipBound::Unbounded => "unbounded".to_string(),
                        FlipBound::Bounded(bound) => format!(
                            "bounded:{}:{}@{}",
                            bound.maximum_abs_effect(),
                            bound.provenance().role(),
                            bound.provenance().digest()
                        ),
                    };
                    let _ = writeln!(
                        output,
                        "- unknown={} reason={} required-to-flip={} plausibility={} suggested-action={}",
                        unknown.kind().name(),
                        unknown.reason(),
                        unknown.required_magnitude(),
                        bound,
                        action_kind_name(unknown.suggested_action())
                    );
                }
            }
        }
        output
    }
}

/// Audit-preserving projection into the legacy three-slice API.
#[derive(Debug, Clone, PartialEq)]
pub struct LegacyUncertaintyProjection {
    original_budget: ContentHash,
    reference_magnitude: f64,
    breakdown: UncertaintyBreakdown,
    numerical_sources: Vec<EngineeringUncertaintyKind>,
    statistical_sources: Vec<EngineeringUncertaintyKind>,
    model_sources: Vec<EngineeringUncertaintyKind>,
}

impl LegacyUncertaintyProjection {
    /// Identity of the exact eight-term source budget.
    #[must_use]
    pub const fn original_budget(&self) -> ContentHash {
        self.original_budget
    }

    /// Finite nonzero QoI magnitude used to convert absolute half-widths to
    /// the relative currency of [`UncertaintyBreakdown`].
    #[must_use]
    pub const fn reference_magnitude(&self) -> f64 {
        self.reference_magnitude
    }

    /// Existing three-slice view. Unknowns map to infinity, never zero.
    #[must_use]
    pub const fn breakdown(&self) -> UncertaintyBreakdown {
        self.breakdown
    }

    /// Sources accounted into the numerical slice.
    #[must_use]
    pub fn numerical_sources(&self) -> &[EngineeringUncertaintyKind] {
        &self.numerical_sources
    }

    /// Sources accounted into the statistical slice.
    #[must_use]
    pub fn statistical_sources(&self) -> &[EngineeringUncertaintyKind] {
        &self.statistical_sources
    }

    /// Sources accounted into the model slice.
    #[must_use]
    pub fn model_sources(&self) -> &[EngineeringUncertaintyKind] {
        &self.model_sources
    }
}

/// Complete, versioned eight-term uncertainty budget for one QoI and unit.
#[derive(Debug, Clone, PartialEq)]
pub struct EngineeringUncertaintyBudget {
    schema_version: u32,
    qoi: BoundedId,
    unit: BoundedId,
    terms: [EngineeringUncertaintyTerm; ENGINEERING_UNCERTAINTY_TERM_COUNT],
}

impl EngineeringUncertaintyBudget {
    /// Admit exactly one term from each mandatory source category.
    pub fn try_new(
        qoi: &str,
        unit: &str,
        terms: Vec<EngineeringUncertaintyTerm>,
    ) -> Result<Self, UncertaintyError> {
        if terms.len() != ENGINEERING_UNCERTAINTY_TERM_COUNT {
            return refuse(
                UncertaintyRule::TermSet,
                format!("budget has {} terms; expected 8", terms.len()),
            );
        }
        let qoi = BoundedId::new(qoi)
            .map_err(|error| UncertaintyError::new(UncertaintyRule::IdBounds, error.to_string()))?;
        let unit = BoundedId::new(unit)
            .map_err(|error| UncertaintyError::new(UncertaintyRule::IdBounds, error.to_string()))?;
        let mut by_kind = BTreeMap::new();
        for term in terms {
            term.value.validate()?;
            let kind = term.kind;
            if by_kind.insert(kind, term).is_some() {
                return refuse(
                    UncertaintyRule::TermSet,
                    format!("duplicate {} term", kind.name()),
                );
            }
        }
        let ordered = EngineeringUncertaintyKind::ALL.map(|kind| {
            by_kind.remove(&kind).ok_or_else(|| {
                UncertaintyError::new(
                    UncertaintyRule::TermSet,
                    format!("missing {} term", kind.name()),
                )
            })
        });
        let terms = ordered.into_iter().collect::<Result<Vec<_>, _>>()?;
        let terms: [EngineeringUncertaintyTerm; ENGINEERING_UNCERTAINTY_TERM_COUNT] = terms
            .try_into()
            .map_err(|_| UncertaintyError::new(UncertaintyRule::TermSet, "term count changed"))?;
        validate_block_membership(&terms)?;
        Ok(Self {
            schema_version: ENGINEERING_UNCERTAINTY_SCHEMA_VERSION,
            qoi,
            unit,
            terms,
        })
    }

    /// Canonical schema version.
    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Quantity-of-interest identity.
    #[must_use]
    pub fn qoi(&self) -> &str {
        self.qoi.as_str()
    }

    /// Unit identity shared by every numeric term.
    #[must_use]
    pub fn unit(&self) -> &str {
        self.unit.as_str()
    }

    /// Terms in canonical source order.
    #[must_use]
    pub fn terms(&self) -> &[EngineeringUncertaintyTerm; ENGINEERING_UNCERTAINTY_TERM_COUNT] {
        &self.terms
    }

    /// Look up one mandatory source term.
    #[must_use]
    pub fn term(&self, kind: EngineeringUncertaintyKind) -> &EngineeringUncertaintyTerm {
        &self.terms[usize::from(kind.code() - 1)]
    }

    /// Replace only roundoff, solver/algebraic, and discretization evidence.
    /// Every other term is cloned bit-for-bit from this budget.
    pub fn apply_numerical_update(
        &self,
        update: &NumericalUncertaintyUpdate,
    ) -> Result<Self, UncertaintyError> {
        let mut terms = self.terms.clone();
        for term in update.terms() {
            terms[usize::from(term.kind.code() - 1)] = term.clone();
        }
        Self::try_new(self.qoi(), self.unit(), Vec::from(terms))
    }

    /// Conservative total honoring only explicit covariance blocks.
    #[must_use]
    pub fn total(&self) -> BudgetTotal {
        let mut known = 0.0;
        let mut unknown_terms = Vec::new();
        let mut counted_blocks = BTreeSet::new();
        for term in &self.terms {
            match &term.value {
                TermValue::Unknown { .. } => unknown_terms.push(term.kind),
                TermValue::CorrelatedBlock(block) => {
                    if counted_blocks.insert(block.content_id()) {
                        known = add_up(known, block.combined_half_width());
                    }
                }
                value => {
                    let Some(half_width) = value.marginal_half_width(term.kind) else {
                        unknown_terms.push(term.kind);
                        continue;
                    };
                    known = add_up(known, half_width);
                }
            }
        }
        if !known.is_finite() {
            BudgetTotal::Unbounded {
                reason: "finite term aggregation overflowed",
            }
        } else if unknown_terms.is_empty() {
            BudgetTotal::Bounded {
                conservative_half_width: known,
            }
        } else {
            BudgetTotal::Unknown {
                known_conservative_half_width: known,
                unknown_terms,
            }
        }
    }

    /// Evaluate a sourced scalar requirement without treating an unknown term
    /// as zero.
    ///
    /// Fully specified terms first form an outward-rounded band around
    /// `nominal`. Unknown terms are unbounded unless the caller supplies a
    /// matching [`UnknownPlausibilityBound`]. Finite plausibility bounds are
    /// summed conservatively because independence is not assumed. Touching the
    /// requirement limit is always [`ComplianceVerdict::Indeterminate`]; only
    /// a strictly positive residual margin or shortfall produces a binary
    /// verdict.
    pub fn assess_requirement(
        &self,
        nominal: f64,
        requirement: &ScalarRequirement,
        plausibility_bounds: &[UnknownPlausibilityBound],
    ) -> Result<ComplianceVerdict, UncertaintyError> {
        let bounds_by_kind =
            validate_requirement_inputs(self, nominal, requirement, plausibility_bounds)?;
        let (known_lower, known_upper) = known_requirement_band(self, nominal);
        let baseline = requirement_baseline(requirement, known_lower, known_upper);
        let bounded_sum = bounds_by_kind
            .values()
            .fold(0.0, |sum, bound| add_up(sum, bound.maximum_abs_effect));
        let has_unbounded = self.terms.iter().any(|term| {
            matches!(term.value, TermValue::Unknown { .. })
                && !bounds_by_kind.contains_key(&term.kind)
        });

        if let Some((compliant, distance)) = baseline
            && !has_unbounded
            && bounded_sum < distance
        {
            let residual = sub_down(distance, bounded_sum);
            return Ok(if compliant {
                ComplianceVerdict::Compliant {
                    margin: residual,
                    budget: self.content_id(),
                    requirement: requirement.clone(),
                }
            } else {
                ComplianceVerdict::NonCompliant {
                    shortfall: residual,
                    budget: self.content_id(),
                    requirement: requirement.clone(),
                }
            });
        }

        let baseline_distance = baseline.map_or(0.0, |(_, distance)| distance);
        Ok(ComplianceVerdict::Indeterminate {
            known_lower,
            known_upper,
            flipping_unknowns: requirement_flipping_unknowns(
                self,
                &bounds_by_kind,
                baseline_distance,
            ),
            budget: self.content_id(),
            requirement: requirement.clone(),
        })
    }

    /// Deterministic dominant source. Unknown sources refuse finite ranking;
    /// otherwise ties prefer the harder-to-shrink source.
    #[must_use]
    pub fn dominant(&self) -> DominantEngineeringTerm {
        let unknown = self
            .terms
            .iter()
            .filter_map(|term| matches!(term.value, TermValue::Unknown { .. }).then_some(term.kind))
            .collect::<Vec<_>>();
        if !unknown.is_empty() {
            return DominantEngineeringTerm::Unknown { terms: unknown };
        }
        let first = DOMINANCE_ORDER[0];
        let mut best = (
            first,
            self.term(first)
                .value
                .marginal_half_width(first)
                .unwrap_or(f64::INFINITY),
        );
        for kind in &DOMINANCE_ORDER[1..] {
            let half_width = self
                .term(*kind)
                .value
                .marginal_half_width(*kind)
                .unwrap_or(f64::INFINITY);
            if half_width > best.1 {
                best = (*kind, half_width);
            }
        }
        DominantEngineeringTerm::Known {
            kind: best.0,
            conservative_half_width: best.1,
        }
    }

    /// Project every source into the legacy three-slice API without dropping
    /// accounting identity. Marginals are linearly summed inside each slice;
    /// an unknown source maps that slice to infinity.
    pub fn project_legacy(
        &self,
        reference_magnitude: f64,
    ) -> Result<LegacyUncertaintyProjection, UncertaintyError> {
        if !reference_magnitude.is_finite() || reference_magnitude == 0.0 {
            return refuse(
                UncertaintyRule::NumericDomain,
                "legacy projection needs a finite nonzero QoI reference magnitude",
            );
        }
        let scale = reference_magnitude.abs();
        let mut numerical_rel = 0.0;
        let mut statistical_rel = 0.0;
        let mut model_rel = 0.0;
        let mut numerical_sources = Vec::new();
        let mut statistical_sources = Vec::new();
        let mut model_sources = Vec::new();
        for term in &self.terms {
            let absolute = term
                .value
                .marginal_half_width(term.kind)
                .unwrap_or(f64::INFINITY);
            let contribution = divide_up(absolute, scale);
            match term.kind.legacy_slice() {
                LegacySlice::Numerical => {
                    numerical_sources.push(term.kind);
                    numerical_rel = add_up(numerical_rel, contribution);
                }
                LegacySlice::Statistical => {
                    statistical_sources.push(term.kind);
                    statistical_rel = add_up(statistical_rel, contribution);
                }
                LegacySlice::Model => {
                    model_sources.push(term.kind);
                    model_rel = add_up(model_rel, contribution);
                }
            }
        }
        Ok(LegacyUncertaintyProjection {
            original_budget: self.content_id(),
            reference_magnitude,
            breakdown: UncertaintyBreakdown {
                numerical_rel,
                statistical_rel,
                model_rel,
            },
            numerical_sources,
            statistical_sources,
            model_sources,
        })
    }

    /// Conservatively compose two budgets for the same QoI/unit. Any mixed
    /// rich representation degrades to an interval upper bound; unknown is
    /// absorbing and negligible is the additive identity.
    pub fn compose(&self, other: &Self) -> Result<Self, UncertaintyError> {
        if self.qoi != other.qoi || self.unit != other.unit {
            return refuse(
                UncertaintyRule::IncompatibleBudget,
                format!(
                    "cannot compose {}/{} with {}/{}",
                    self.qoi(),
                    self.unit(),
                    other.qoi(),
                    other.unit()
                ),
            );
        }
        let left_budget_id = self.content_id();
        let right_budget_id = other.content_id();
        let mut terms = Vec::with_capacity(ENGINEERING_UNCERTAINTY_TERM_COUNT);
        for kind in EngineeringUncertaintyKind::ALL {
            let left = self.term(kind);
            let right = other.term(kind);
            let value = compose_values(kind, &left.value, &right.value)?;
            let mut provenance_bytes = Vec::with_capacity(65);
            provenance_bytes.push(kind.code());
            provenance_bytes.extend_from_slice(left_budget_id.as_bytes());
            provenance_bytes.extend_from_slice(right_budget_id.as_bytes());
            let provenance = UncertaintyArtifactRef::new(
                "uncertainty:composed",
                hash_domain(COMPOSITION_DOMAIN, &provenance_bytes),
            )?;
            terms.push(EngineeringUncertaintyTerm::try_new(
                kind, value, provenance,
            )?);
        }
        Self::try_new(self.qoi(), self.unit(), terms)
    }

    /// Canonical versioned transport bytes.
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(MAGIC);
        out.extend_from_slice(&self.schema_version.to_be_bytes());
        push_str(&mut out, self.qoi.as_str());
        push_str(&mut out, self.unit.as_str());
        out.push(ENGINEERING_UNCERTAINTY_TERM_COUNT as u8);
        for term in &self.terms {
            out.push(term.kind.code());
            encode_artifact(&mut out, &term.provenance);
            encode_value(&mut out, &term.value);
        }
        out
    }

    /// Domain-separated identity over exact canonical bytes.
    #[must_use]
    pub fn content_id(&self) -> ContentHash {
        hash_domain(IDENTITY_DOMAIN, &self.canonical_bytes())
    }

    /// Deterministic human-readable report with one provenance link per term.
    #[must_use]
    pub fn render_report(&self) -> String {
        let mut output = String::new();
        let _ = writeln!(
            output,
            "engineering-uncertainty-v{} qoi={} unit={} identity={}",
            self.schema_version,
            self.qoi(),
            self.unit(),
            self.content_id()
        );
        for term in &self.terms {
            let representation = match &term.value {
                TermValue::IntervalBound { lower, upper } => {
                    format!("interval-half-width=[{lower},{upper}]")
                }
                TermValue::Distribution(summary) => format!(
                    "distribution-half-width={} level={} replay={}@{}",
                    summary.conservative_half_width,
                    summary.level,
                    summary.replay.role(),
                    summary.replay.digest()
                ),
                TermValue::Ensemble(summary) => format!(
                    "ensemble-half-width={} members={} replay={}@{}",
                    summary.conservative_half_width,
                    summary.member_count,
                    summary.replay.role(),
                    summary.replay.digest()
                ),
                TermValue::CorrelatedBlock(block) => format!(
                    "correlated-block={} covariance={}@{}",
                    block.id(),
                    block.covariance_artifact().role(),
                    block.covariance_artifact().digest()
                ),
                TermValue::Unknown { reason } => format!("unknown={reason}"),
                TermValue::Negligible { justification } => {
                    format!("negligible={justification}")
                }
            };
            let _ = writeln!(
                output,
                "- {}: {}; provenance={}@{}",
                term.kind.name(),
                representation,
                term.provenance.role(),
                term.provenance.digest()
            );
        }
        match self.total() {
            BudgetTotal::Bounded {
                conservative_half_width,
            } => {
                let _ = writeln!(output, "total=bounded:{conservative_half_width}");
            }
            BudgetTotal::Unknown {
                known_conservative_half_width,
                unknown_terms,
            } => {
                let names = unknown_terms
                    .iter()
                    .map(|kind| kind.name())
                    .collect::<Vec<_>>()
                    .join(",");
                let _ = writeln!(
                    output,
                    "total=unknown:known-half-width={known_conservative_half_width}:terms={names}"
                );
            }
            BudgetTotal::Unbounded { reason } => {
                let _ = writeln!(output, "total=unbounded:{reason}");
            }
        }
        output
    }

    /// Decode, revalidate, and require a canonical byte-for-byte fixed point.
    pub fn decode(bytes: &[u8]) -> Result<Self, UncertaintyCodecError> {
        if bytes.len() > MAX_UNCERTAINTY_CANONICAL_BYTES {
            return Err(UncertaintyCodecError::at(
                0,
                "transport exceeds the 1 MiB budget",
            ));
        }
        let mut reader = Reader { bytes, pos: 0 };
        if reader.take(4)? != MAGIC {
            return Err(UncertaintyCodecError::at(0, "bad magic"));
        }
        let schema_version = reader.u32()?;
        if schema_version != ENGINEERING_UNCERTAINTY_SCHEMA_VERSION {
            return Err(UncertaintyCodecError::at(
                reader.pos,
                format!("unknown schema version {schema_version}"),
            ));
        }
        let qoi = reader.string()?;
        let unit = reader.string()?;
        let count = usize::from(reader.u8()?);
        if count != ENGINEERING_UNCERTAINTY_TERM_COUNT {
            return Err(UncertaintyCodecError::at(
                reader.pos,
                format!("term count {count} is not 8"),
            ));
        }
        let mut terms = Vec::with_capacity(count);
        for _ in 0..count {
            let kind = EngineeringUncertaintyKind::from_code(reader.u8()?)
                .ok_or_else(|| UncertaintyCodecError::at(reader.pos, "bad source tag"))?;
            let provenance = decode_artifact(&mut reader)?;
            let value = decode_value(&mut reader)?;
            let term = EngineeringUncertaintyTerm::try_new(kind, value, provenance)
                .map_err(|error| UncertaintyCodecError::at(reader.pos, error.to_string()))?;
            terms.push(term);
        }
        if reader.pos != bytes.len() {
            return Err(UncertaintyCodecError::at(reader.pos, "trailing bytes"));
        }
        let budget = Self::try_new(&qoi, &unit, terms)
            .map_err(|error| UncertaintyCodecError::at(reader.pos, error.to_string()))?;
        if budget.canonical_bytes() != bytes {
            return Err(UncertaintyCodecError::at(
                reader.pos,
                "transport is not canonical",
            ));
        }
        Ok(budget)
    }
}

#[derive(Debug, Clone, Copy)]
enum LegacySlice {
    Numerical,
    Statistical,
    Model,
}

fn compose_values(
    kind: EngineeringUncertaintyKind,
    left: &TermValue,
    right: &TermValue,
) -> Result<TermValue, UncertaintyError> {
    match (left, right) {
        (TermValue::Unknown { reason: left }, TermValue::Unknown { reason: right }) => {
            TermValue::unknown(combine_text("left", left, "right", right)?)
        }
        (TermValue::Unknown { reason }, _) => {
            TermValue::unknown(combine_text("left", reason, "right", "specified")?)
        }
        (_, TermValue::Unknown { reason }) => {
            TermValue::unknown(combine_text("left", "specified", "right", reason)?)
        }
        (
            TermValue::Negligible {
                justification: left,
            },
            TermValue::Negligible {
                justification: right,
            },
        ) => TermValue::negligible(combine_text("left", left, "right", right)?),
        (TermValue::Negligible { .. }, value) | (value, TermValue::Negligible { .. }) => {
            Ok(value.clone())
        }
        (left, right) => {
            let left = left.marginal_half_width(kind).ok_or_else(|| {
                UncertaintyError::new(UncertaintyRule::NumericDomain, "missing left marginal")
            })?;
            let right = right.marginal_half_width(kind).ok_or_else(|| {
                UncertaintyError::new(UncertaintyRule::NumericDomain, "missing right marginal")
            })?;
            let upper = add_up(left, right);
            if upper.is_finite() {
                TermValue::interval(0.0, upper)
            } else {
                TermValue::unknown("finite term composition overflowed")
            }
        }
    }
}

fn combine_text(
    left_label: &str,
    left: &str,
    right_label: &str,
    right: &str,
) -> Result<String, UncertaintyError> {
    let combined = format!("{left_label}: {left}; {right_label}: {right}");
    admit_text("composed explanation", &combined)?;
    Ok(combined)
}

fn validate_block_membership(
    terms: &[EngineeringUncertaintyTerm; ENGINEERING_UNCERTAINTY_TERM_COUNT],
) -> Result<(), UncertaintyError> {
    let mut seen: BTreeMap<ContentHash, (&CovarianceBlock, BTreeSet<EngineeringUncertaintyKind>)> =
        BTreeMap::new();
    for term in terms {
        let TermValue::CorrelatedBlock(block) = &term.value else {
            continue;
        };
        if !block.members.contains(&term.kind) {
            return refuse(
                UncertaintyRule::CovarianceMembership,
                format!("block {} omits {}", block.id(), term.kind.name()),
            );
        }
        let id = block.content_id();
        let entry = seen.entry(id).or_insert_with(|| (block, BTreeSet::new()));
        entry.1.insert(term.kind);
    }
    for (block, actual) in seen.values() {
        let expected = block.members.iter().copied().collect::<BTreeSet<_>>();
        if *actual != expected {
            return refuse(
                UncertaintyRule::CovarianceMembership,
                format!(
                    "block {} declares {:?} but is referenced by {:?}",
                    block.id(),
                    expected,
                    actual
                ),
            );
        }
    }
    Ok(())
}

fn validate_covariance(values: &[f64], n: usize) -> Result<(), UncertaintyError> {
    if !(2..=ENGINEERING_UNCERTAINTY_TERM_COUNT).contains(&n) || values.len() != n * n {
        return refuse(
            UncertaintyRule::CovarianceMatrix,
            "covariance dimensions are outside the admitted shape",
        );
    }
    if values.iter().any(|value| !value.is_finite()) {
        return refuse(
            UncertaintyRule::CovarianceMatrix,
            "covariance contains a non-finite value",
        );
    }
    for row in 0..n {
        if values[row * n + row] < 0.0 {
            return refuse(
                UncertaintyRule::CovarianceMatrix,
                "covariance diagonal is negative",
            );
        }
        for column in row + 1..n {
            if values[row * n + column].to_bits() != values[column * n + row].to_bits() {
                return refuse(
                    UncertaintyRule::CovarianceMatrix,
                    "covariance is not bit-symmetric",
                );
            }
        }
    }

    // Pivot-free LDL^T semidefinite test. A zero pivot is admitted only when
    // the corresponding residual column is zero; otherwise the matrix has a
    // negative direction. The tolerance is scale-relative and used only for
    // rejecting roundoff-sized decomposition residue, never for changing the
    // retained covariance bytes.
    let scale = values
        .iter()
        .fold(1.0_f64, |acc, value| acc.max(value.abs()));
    let tolerance = scale * f64::EPSILON * 64.0 * n as f64;
    let mut lower = vec![0.0; n * n];
    let mut diagonal = vec![0.0; n];
    for i in 0..n {
        let mut pivot = values[i * n + i];
        for k in 0..i {
            pivot -= lower[i * n + k] * lower[i * n + k] * diagonal[k];
        }
        if pivot < -tolerance {
            return refuse(
                UncertaintyRule::CovarianceMatrix,
                "covariance is not positive semidefinite",
            );
        }
        diagonal[i] = if pivot.abs() <= tolerance { 0.0 } else { pivot };
        lower[i * n + i] = 1.0;
        for j in i + 1..n {
            let mut residual = values[j * n + i];
            for k in 0..i {
                residual -= lower[j * n + k] * lower[i * n + k] * diagonal[k];
            }
            if diagonal[i] == 0.0 {
                if residual.abs() > tolerance {
                    return refuse(
                        UncertaintyRule::CovarianceMatrix,
                        "covariance has a nonzero residual under a zero pivot",
                    );
                }
                lower[j * n + i] = 0.0;
            } else {
                lower[j * n + i] = residual / diagonal[i];
            }
        }
    }
    Ok(())
}

fn non_negative_finite(value: f64) -> bool {
    value.is_finite() && value >= 0.0 && !is_negative_zero(value)
}

fn is_negative_zero(value: f64) -> bool {
    value.to_bits() == (-0.0_f64).to_bits()
}

fn admit_text(label: &str, value: &str) -> Result<(), UncertaintyError> {
    if value.trim().is_empty()
        || value.len() > MAX_UNCERTAINTY_TEXT_BYTES
        || value.chars().any(char::is_control)
    {
        return refuse(
            UncertaintyRule::TextBounds,
            format!(
                "{label} must contain 1..={MAX_UNCERTAINTY_TEXT_BYTES} non-control UTF-8 bytes"
            ),
        );
    }
    Ok(())
}

fn next_up(value: f64) -> f64 {
    if value.is_nan() || value == f64::INFINITY {
        return value;
    }
    if value == -0.0 {
        return f64::from_bits(1);
    }
    let bits = value.to_bits();
    f64::from_bits(if value >= 0.0 { bits + 1 } else { bits - 1 })
}

fn next_down(value: f64) -> f64 {
    -next_up(-value)
}

fn sub_down(left: f64, right: f64) -> f64 {
    if right == 0.0 {
        left
    } else {
        next_down(left - right)
    }
}

fn validate_requirement_inputs(
    budget: &EngineeringUncertaintyBudget,
    nominal: f64,
    requirement: &ScalarRequirement,
    plausibility_bounds: &[UnknownPlausibilityBound],
) -> Result<BTreeMap<EngineeringUncertaintyKind, UnknownPlausibilityBound>, UncertaintyError> {
    if !nominal.is_finite() {
        return refuse(
            UncertaintyRule::RequirementAssessment,
            "requirement assessment nominal must be finite",
        );
    }
    if budget.qoi() != requirement.qoi() || budget.unit() != requirement.unit() {
        return refuse(
            UncertaintyRule::RequirementAssessment,
            format!(
                "requirement {} governs {} [{}], but budget is {} [{}]",
                requirement.id(),
                requirement.qoi(),
                requirement.unit(),
                budget.qoi(),
                budget.unit()
            ),
        );
    }

    let mut bounds_by_kind = BTreeMap::new();
    for bound in plausibility_bounds {
        if !matches!(budget.term(bound.kind).value, TermValue::Unknown { .. }) {
            return refuse(
                UncertaintyRule::RequirementAssessment,
                format!(
                    "{} has a plausibility bound but its budget term is not Unknown",
                    bound.kind.name()
                ),
            );
        }
        if bounds_by_kind.insert(bound.kind, bound.clone()).is_some() {
            return refuse(
                UncertaintyRule::RequirementAssessment,
                format!("duplicate {} plausibility bound", bound.kind.name()),
            );
        }
    }
    Ok(bounds_by_kind)
}

fn known_requirement_band(budget: &EngineeringUncertaintyBudget, nominal: f64) -> (f64, f64) {
    let half_width = match budget.total() {
        BudgetTotal::Bounded {
            conservative_half_width,
        }
        | BudgetTotal::Unknown {
            known_conservative_half_width: conservative_half_width,
            ..
        } => conservative_half_width,
        BudgetTotal::Unbounded { .. } => f64::INFINITY,
    };
    if half_width.is_finite() {
        (sub_down(nominal, half_width), add_up(nominal, half_width))
    } else {
        (f64::NEG_INFINITY, f64::INFINITY)
    }
}

fn requirement_baseline(
    requirement: &ScalarRequirement,
    known_lower: f64,
    known_upper: f64,
) -> Option<(bool, f64)> {
    let limit = requirement.limit();
    match requirement.relation() {
        RequirementRelation::AtMost if known_upper < limit => {
            Some((true, sub_down(limit, known_upper)))
        }
        RequirementRelation::AtMost if known_lower > limit => {
            Some((false, sub_down(known_lower, limit)))
        }
        RequirementRelation::AtLeast if known_lower > limit => {
            Some((true, sub_down(known_lower, limit)))
        }
        RequirementRelation::AtLeast if known_upper < limit => {
            Some((false, sub_down(limit, known_upper)))
        }
        RequirementRelation::AtMost | RequirementRelation::AtLeast => None,
    }
}

fn requirement_flipping_unknowns(
    budget: &EngineeringUncertaintyBudget,
    bounds_by_kind: &BTreeMap<EngineeringUncertaintyKind, UnknownPlausibilityBound>,
    baseline_distance: f64,
) -> Vec<FlippingUnknown> {
    budget
        .terms
        .iter()
        .filter_map(|term| {
            let TermValue::Unknown { reason } = &term.value else {
                return None;
            };
            let bound = bounds_by_kind.get(&term.kind).cloned();
            let other_bounded = bounds_by_kind
                .iter()
                .filter(|(other, _)| **other != term.kind)
                .fold(0.0, |sum, (_, other)| add_up(sum, other.maximum_abs_effect));
            let required_magnitude = if baseline_distance > other_bounded {
                sub_down(baseline_distance, other_bounded)
            } else {
                0.0
            };
            bound
                .as_ref()
                .is_none_or(|bound| bound.maximum_abs_effect >= required_magnitude)
                .then(|| FlippingUnknown {
                    kind: term.kind,
                    reason: reason.clone(),
                    required_magnitude,
                    bound: bound.map_or(FlipBound::Unbounded, FlipBound::Bounded),
                    suggested_action: suggested_resolution_action(term.kind),
                })
        })
        .collect()
}

const fn suggested_resolution_action(
    kind: EngineeringUncertaintyKind,
) -> crate::action::ActionKind {
    match kind {
        EngineeringUncertaintyKind::Roundoff | EngineeringUncertaintyKind::SolverAlgebraic => {
            crate::action::ActionKind::SolverTolerance
        }
        EngineeringUncertaintyKind::Discretization => crate::action::ActionKind::MeshRefinement,
        EngineeringUncertaintyKind::Geometry => crate::action::ActionKind::RepresentationEscalation,
        EngineeringUncertaintyKind::Parameters => crate::action::ActionKind::MaterialCouponTest,
        EngineeringUncertaintyKind::BoundaryConditions
        | EngineeringUncertaintyKind::Measurement => crate::action::ActionKind::SensorCampaign,
        EngineeringUncertaintyKind::ModelForm => crate::action::ActionKind::Falsification,
    }
}

const fn action_kind_name(kind: crate::action::ActionKind) -> &'static str {
    match kind {
        crate::action::ActionKind::SolverTolerance => "solver-tolerance",
        crate::action::ActionKind::MeshRefinement => "mesh-refinement",
        crate::action::ActionKind::TimeRefinement => "time-refinement",
        crate::action::ActionKind::RepresentationEscalation => "representation-escalation",
        crate::action::ActionKind::UqSamples => "uq-samples",
        crate::action::ActionKind::MaterialCouponTest => "material-coupon-test",
        crate::action::ActionKind::SensorCampaign => "sensor-campaign",
        crate::action::ActionKind::Falsification => "falsification",
        crate::action::ActionKind::StandardsObligation => "standards-obligation",
        crate::action::ActionKind::Refusal => "refusal",
    }
}

fn add_up(left: f64, right: f64) -> f64 {
    if left == f64::INFINITY || right == f64::INFINITY {
        f64::INFINITY
    } else if left == 0.0 {
        right
    } else if right == 0.0 {
        left
    } else {
        next_up(left + right)
    }
}

fn divide_up(value: f64, positive_scale: f64) -> f64 {
    if value == 0.0 || value == f64::INFINITY {
        value
    } else {
        next_up(value / positive_scale)
    }
}

fn push_str(out: &mut Vec<u8>, value: &str) {
    out.extend_from_slice(&(value.len() as u32).to_be_bytes());
    out.extend_from_slice(value.as_bytes());
}

fn encode_artifact(out: &mut Vec<u8>, artifact: &UncertaintyArtifactRef) {
    push_str(out, artifact.role.as_str());
    out.extend_from_slice(artifact.digest.as_bytes());
}

fn encode_block(out: &mut Vec<u8>, block: &CovarianceBlock) {
    push_str(out, block.id.as_str());
    encode_artifact(out, &block.covariance_artifact);
    out.push(block.members.len() as u8);
    for member in &block.members {
        out.push(member.code());
    }
    out.push(block.members.len() as u8);
    for value in &block.covariance {
        out.extend_from_slice(&value.to_bits().to_be_bytes());
    }
}

fn encode_value(out: &mut Vec<u8>, value: &TermValue) {
    match value {
        TermValue::IntervalBound { lower, upper } => {
            out.push(1);
            out.extend_from_slice(&lower.to_bits().to_be_bytes());
            out.extend_from_slice(&upper.to_bits().to_be_bytes());
        }
        TermValue::Distribution(summary) => {
            out.push(2);
            for value in [
                summary.mean,
                summary.standard_deviation,
                summary.conservative_half_width,
                summary.level,
            ] {
                out.extend_from_slice(&value.to_bits().to_be_bytes());
            }
            encode_artifact(out, &summary.replay);
        }
        TermValue::Ensemble(summary) => {
            out.push(3);
            out.extend_from_slice(&summary.member_count.to_be_bytes());
            out.extend_from_slice(&summary.conservative_half_width.to_bits().to_be_bytes());
            encode_artifact(out, &summary.replay);
        }
        TermValue::CorrelatedBlock(block) => {
            out.push(4);
            encode_block(out, block);
        }
        TermValue::Unknown { reason } => {
            out.push(5);
            push_str(out, reason);
        }
        TermValue::Negligible { justification } => {
            out.push(6);
            push_str(out, justification);
        }
    }
}

fn decode_artifact(
    reader: &mut Reader<'_>,
) -> Result<UncertaintyArtifactRef, UncertaintyCodecError> {
    let role = reader.string()?;
    let digest = ContentHash::from_slice(reader.take(32)?)
        .ok_or_else(|| UncertaintyCodecError::at(reader.pos, "bad artifact digest"))?;
    UncertaintyArtifactRef::new(&role, digest)
        .map_err(|error| UncertaintyCodecError::at(reader.pos, error.to_string()))
}

fn decode_block(reader: &mut Reader<'_>) -> Result<CovarianceBlock, UncertaintyCodecError> {
    let id = reader.string()?;
    let artifact = decode_artifact(reader)?;
    let count = usize::from(reader.u8()?);
    if !(2..=ENGINEERING_UNCERTAINTY_TERM_COUNT).contains(&count) {
        return Err(UncertaintyCodecError::at(
            reader.pos,
            "covariance member count outside 2..=8",
        ));
    }
    let mut members = Vec::with_capacity(count);
    for _ in 0..count {
        members.push(
            EngineeringUncertaintyKind::from_code(reader.u8()?)
                .ok_or_else(|| UncertaintyCodecError::at(reader.pos, "bad covariance member"))?,
        );
    }
    let matrix_count = usize::from(reader.u8()?);
    if matrix_count != count {
        return Err(UncertaintyCodecError::at(
            reader.pos,
            "covariance matrix dimension disagrees with member count",
        ));
    }
    let entries = count
        .checked_mul(count)
        .ok_or_else(|| UncertaintyCodecError::at(reader.pos, "covariance size overflow"))?;
    let mut covariance = Vec::with_capacity(entries);
    for _ in 0..entries {
        covariance.push(reader.f64()?);
    }
    CovarianceBlock::try_new(&id, artifact, members, covariance)
        .map_err(|error| UncertaintyCodecError::at(reader.pos, error.to_string()))
}

fn decode_value(reader: &mut Reader<'_>) -> Result<TermValue, UncertaintyCodecError> {
    let value = match reader.u8()? {
        1 => TermValue::IntervalBound {
            lower: reader.f64()?,
            upper: reader.f64()?,
        },
        2 => TermValue::Distribution(DistributionTerm {
            mean: reader.f64()?,
            standard_deviation: reader.f64()?,
            conservative_half_width: reader.f64()?,
            level: reader.f64()?,
            replay: decode_artifact(reader)?,
        }),
        3 => TermValue::Ensemble(EnsembleTerm {
            member_count: reader.u32()?,
            conservative_half_width: reader.f64()?,
            replay: decode_artifact(reader)?,
        }),
        4 => TermValue::CorrelatedBlock(decode_block(reader)?),
        5 => TermValue::Unknown {
            reason: reader.string()?,
        },
        6 => TermValue::Negligible {
            justification: reader.string()?,
        },
        tag => {
            return Err(UncertaintyCodecError::at(
                reader.pos,
                format!("unknown term-value tag {tag}"),
            ));
        }
    };
    value
        .validate()
        .map_err(|error| UncertaintyCodecError::at(reader.pos, error.to_string()))?;
    Ok(value)
}

struct Reader<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    fn take(&mut self, count: usize) -> Result<&'a [u8], UncertaintyCodecError> {
        let end = self
            .pos
            .checked_add(count)
            .ok_or_else(|| UncertaintyCodecError::at(self.pos, "offset overflow"))?;
        let bytes = self
            .bytes
            .get(self.pos..end)
            .ok_or_else(|| UncertaintyCodecError::at(self.pos, "truncated transport"))?;
        self.pos = end;
        Ok(bytes)
    }

    fn u8(&mut self) -> Result<u8, UncertaintyCodecError> {
        Ok(self.take(1)?[0])
    }

    fn u32(&mut self) -> Result<u32, UncertaintyCodecError> {
        let bytes: [u8; 4] = self
            .take(4)?
            .try_into()
            .map_err(|_| UncertaintyCodecError::at(self.pos, "bad u32"))?;
        Ok(u32::from_be_bytes(bytes))
    }

    fn f64(&mut self) -> Result<f64, UncertaintyCodecError> {
        let bytes: [u8; 8] = self
            .take(8)?
            .try_into()
            .map_err(|_| UncertaintyCodecError::at(self.pos, "bad f64"))?;
        Ok(f64::from_bits(u64::from_be_bytes(bytes)))
    }

    fn string(&mut self) -> Result<String, UncertaintyCodecError> {
        let length = usize::try_from(self.u32()?)
            .map_err(|_| UncertaintyCodecError::at(self.pos, "string length overflow"))?;
        if length > MAX_UNCERTAINTY_TEXT_BYTES {
            return Err(UncertaintyCodecError::at(
                self.pos,
                "string exceeds the 1024-byte field budget",
            ));
        }
        let start = self.pos;
        let bytes = self.take(length)?;
        std::str::from_utf8(bytes)
            .map(str::to_owned)
            .map_err(|_| UncertaintyCodecError::at(start, "string is not UTF-8"))
    }
}
