//! L6 decision projection over already-admitted evidence artifacts.
//!
//! This module deliberately does not estimate a quantity, interpret a safety
//! factor, recompute compliance, price an action, or verify scientific
//! evidence. It checks that lower-layer artifacts agree and then gives that
//! exact decision view a deterministic replay identity.

use core::{fmt, marker::PhantomData};
use std::fmt::Write as _;

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::{
    action::ActionKind,
    uncertainty::{
        BudgetAttribution, ComplianceVerdict, DecisionAttribution, EngineeringUncertaintyBudget,
        FlippingUnknown, RequirementRelation, ScalarRequirement, UncertaintyArtifactRef,
        UncertaintyAttribution,
    },
    vv::{ArtifactKind, ArtifactRef},
};
use fs_package::VerifiedPackage;
use fs_voi::{RecommendedEvidence, UnknownResolutionRecommendation};

/// Semantic identity version of an L6 decision assessment.
pub const DECISION_ASSESSMENT_IDENTITY_VERSION: u32 = 2;
/// Domain-separated identity namespace of an L6 decision assessment.
pub const DECISION_ASSESSMENT_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-session.decision-assessment.v2";

const MAX_DECISION_ID_BYTES: usize = 128;
const MAX_AUTHORITY_FIELD_BYTES: usize = 512;

/// Closed requirement-authority families understood by the L6 decision schema.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RequirementAuthorityKind {
    /// Published standard or code clause.
    Standard,
    /// Manufacturer datasheet or product specification.
    Datasheet,
    /// Versioned internal engineering policy.
    InternalPolicy,
    /// Explicit user declaration without external-document authority.
    UserDeclaration,
}

impl RequirementAuthorityKind {
    const fn wire_tag(self) -> u8 {
        match self {
            Self::Standard => 1,
            Self::Datasheet => 2,
            Self::InternalPolicy => 3,
            Self::UserDeclaration => 4,
        }
    }

    /// Stable reviewer-facing authority-family slug.
    #[must_use]
    pub const fn slug(self) -> &'static str {
        match self {
            Self::Standard => "standard",
            Self::Datasheet => "datasheet",
            Self::InternalPolicy => "internal-policy",
            Self::UserDeclaration => "user-declaration",
        }
    }
}

/// Human-auditable lineage retained alongside one content-bound authority.
///
/// The fields identify what the caller declared; they do not authenticate the
/// document or establish that its clause applies.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RequirementAuthority {
    kind: RequirementAuthorityKind,
    document: String,
    version: String,
    locator: String,
}

impl RequirementAuthority {
    /// Admit trim-canonical, control-free lineage with bounded retained bytes.
    pub fn try_new(
        kind: RequirementAuthorityKind,
        document: impl Into<String>,
        version: impl Into<String>,
        locator: impl Into<String>,
    ) -> Result<Self, DecisionAssessmentError> {
        Ok(Self {
            kind,
            document: admit_authority_field("requirement.source.document", document.into())?,
            version: admit_authority_field("requirement.source.version", version.into())?,
            locator: admit_authority_field("requirement.source.locator", locator.into())?,
        })
    }

    /// Authority family.
    #[must_use]
    pub const fn kind(&self) -> RequirementAuthorityKind {
        self.kind
    }

    /// Stable source document or declaration identity.
    #[must_use]
    pub fn document(&self) -> &str {
        &self.document
    }

    /// Exact edition, revision, or semantic version.
    #[must_use]
    pub fn version(&self) -> &str {
        &self.version
    }

    /// Clause, table, section, or declaration locator.
    #[must_use]
    pub fn locator(&self) -> &str {
        &self.locator
    }
}

/// Typed, content-bound reference to the quantity evidence consumed by a
/// decision assessment.
///
/// `Q` distinguishes quantities at compile time. The explicit QoI, unit, and
/// schema strings remain in the runtime identity so type erasure cannot erase
/// their engineering meaning.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EvidenceRef<Q> {
    qoi: String,
    unit: String,
    schema: String,
    artifact: ContentHash,
    marker: PhantomData<fn() -> Q>,
}

impl<Q> EvidenceRef<Q> {
    /// Admit a nonzero content identity with bounded canonical semantics.
    pub fn try_new(
        qoi: impl Into<String>,
        unit: impl Into<String>,
        schema: impl Into<String>,
        artifact: ContentHash,
    ) -> Result<Self, DecisionAssessmentError> {
        let qoi = admit_id("quantity.qoi", qoi.into())?;
        let unit = admit_id("quantity.unit", unit.into())?;
        let schema = admit_id("quantity.schema", schema.into())?;
        require_nonzero_hash("quantity.artifact", artifact)?;
        Ok(Self {
            qoi,
            unit,
            schema,
            artifact,
            marker: PhantomData,
        })
    }

    /// Quantity-of-interest identity.
    #[must_use]
    pub fn qoi(&self) -> &str {
        &self.qoi
    }

    /// Exact engineering unit.
    #[must_use]
    pub fn unit(&self) -> &str {
        &self.unit
    }

    /// Semantic schema of the referenced quantity artifact.
    #[must_use]
    pub fn schema(&self) -> &str {
        &self.schema
    }

    /// Ledger/content-store key of the exact quantity evidence.
    #[must_use]
    pub const fn artifact(&self) -> ContentHash {
        self.artifact
    }
}

/// Safety factor already applied by the sourced requirement authority.
///
/// L6 records the value and policy artifact but intentionally does not invent
/// a universal multiplication/division rule. The nested scalar requirement's
/// limit is the effective limit that the compliance evaluator consumed.
#[derive(Debug, Clone, PartialEq)]
pub struct AppliedSafetyFactor {
    value: f64,
    policy: UncertaintyArtifactRef,
}

impl AppliedSafetyFactor {
    /// Admit a finite factor of at least one and a nonzero policy artifact.
    pub fn try_new(
        value: f64,
        policy: UncertaintyArtifactRef,
    ) -> Result<Self, DecisionAssessmentError> {
        if !value.is_finite() || value < 1.0 {
            return Err(DecisionAssessmentError::InvalidSafetyFactor);
        }
        require_nonzero_hash("requirement.safety-factor-policy", policy.digest())?;
        Ok(Self { value, policy })
    }

    /// Declared factor already reflected in the effective requirement limit.
    #[must_use]
    pub const fn value(&self) -> f64 {
        self.value
    }

    /// Authority defining how the factor was applied.
    #[must_use]
    pub const fn policy(&self) -> &UncertaintyArtifactRef {
        &self.policy
    }
}

/// Effective sourced scalar requirement plus its safety-factor authority.
#[derive(Debug, Clone, PartialEq)]
pub struct DecisionRequirement {
    scalar: ScalarRequirement,
    source: RequirementAuthority,
    safety_factor: AppliedSafetyFactor,
    safety_factor_source: RequirementAuthority,
}

impl DecisionRequirement {
    /// Pair an already-effective scalar requirement with complete source and
    /// factor-policy lineage.
    pub fn try_new(
        scalar: ScalarRequirement,
        source: RequirementAuthority,
        safety_factor: AppliedSafetyFactor,
        safety_factor_source: RequirementAuthority,
    ) -> Result<Self, DecisionAssessmentError> {
        require_nonzero_hash("requirement.provenance", scalar.provenance().digest())?;
        Ok(Self {
            scalar,
            source,
            safety_factor,
            safety_factor_source,
        })
    }

    /// Exact requirement consumed by lower-layer compliance evaluation.
    #[must_use]
    pub const fn scalar(&self) -> &ScalarRequirement {
        &self.scalar
    }

    /// Exact declared source lineage for the effective requirement.
    #[must_use]
    pub const fn source(&self) -> &RequirementAuthority {
        &self.source
    }

    /// Declared factor and policy authority.
    #[must_use]
    pub const fn safety_factor(&self) -> &AppliedSafetyFactor {
        &self.safety_factor
    }

    /// Exact declared source lineage for the factor application policy.
    #[must_use]
    pub const fn safety_factor_source(&self) -> &RequirementAuthority {
        &self.safety_factor_source
    }
}

/// Pure L6 projection of one engineering decision and its replay authorities.
#[derive(Debug, Clone, PartialEq)]
pub struct DecisionAssessment<Q> {
    quantity: EvidenceRef<Q>,
    requirement: DecisionRequirement,
    context: ArtifactRef,
    compliance: ComplianceVerdict,
    budget: EngineeringUncertaintyBudget,
    attribution: UncertaintyAttribution,
    actions: Vec<UnknownResolutionRecommendation>,
    replay_package: ContentHash,
    content_hash: ContentHash,
}

impl<Q> DecisionAssessment<Q> {
    /// Assemble a decision projection from already-created lower-layer
    /// artifacts. This is offline and side-effect free; no authority is
    /// synthesized when an input is absent or inconsistent.
    #[allow(clippy::too_many_arguments)] // Every decision authority remains explicit.
    pub fn try_assemble(
        quantity: EvidenceRef<Q>,
        requirement: DecisionRequirement,
        context: ArtifactRef,
        compliance: ComplianceVerdict,
        budget: EngineeringUncertaintyBudget,
        attribution: UncertaintyAttribution,
        actions: Vec<UnknownResolutionRecommendation>,
        replay_package: &VerifiedPackage,
    ) -> Result<Self, DecisionAssessmentError> {
        validate_bindings(
            &quantity,
            &requirement,
            &context,
            &compliance,
            &budget,
            &attribution,
            &actions,
            replay_package,
        )?;
        let replay_package = replay_package.report().merkle_root();
        let mut assessment = Self {
            quantity,
            requirement,
            context,
            compliance,
            budget,
            attribution,
            actions,
            replay_package,
            content_hash: ContentHash([0; 32]),
        };
        assessment.content_hash = assessment.recompute_content_hash();
        Ok(assessment)
    }

    /// Typed quantity evidence reference.
    #[must_use]
    pub const fn quantity(&self) -> &EvidenceRef<Q> {
        &self.quantity
    }

    /// Effective requirement and safety-factor authority.
    #[must_use]
    pub const fn requirement(&self) -> &DecisionRequirement {
        &self.requirement
    }

    /// Exact context-of-use artifact reference.
    #[must_use]
    pub const fn context(&self) -> &ArtifactRef {
        &self.context
    }

    /// Lower-layer tri-state compliance result.
    #[must_use]
    pub const fn compliance(&self) -> &ComplianceVerdict {
        &self.compliance
    }

    /// Exact eight-term uncertainty budget.
    #[must_use]
    pub const fn budget(&self) -> &EngineeringUncertaintyBudget {
        &self.budget
    }

    /// Paired budget-magnitude and decision-influence attribution.
    #[must_use]
    pub const fn attribution(&self) -> &UncertaintyAttribution {
        &self.attribution
    }

    /// Unknown evidence gaps and the minimum adverse effects that can change
    /// the current verdict.
    #[must_use]
    pub fn flip_conditions(&self) -> &[FlippingUnknown] {
        self.compliance.flipping_unknowns()
    }

    /// Largest finite source/group in the conservative budget-magnitude view.
    /// Explicit unknown budget links remain separately available through
    /// [`UncertaintyAttribution::unknown_budget`].
    #[must_use]
    pub fn largest_known_budget_link(&self) -> Option<&BudgetAttribution> {
        self.attribution.known_budget_ranked().first()
    }

    /// Source/group with the strongest one-group-at-a-time requirement
    /// influence.
    #[must_use]
    pub fn strongest_decision_link(&self) -> Option<&DecisionAttribution> {
        self.attribution.decision_ranked().first()
    }

    /// Cost-aware or explicitly unpriced actions for verdict-flipping unknowns.
    #[must_use]
    pub fn actions(&self) -> &[UnknownResolutionRecommendation] {
        &self.actions
    }

    /// Root of the structurally and policy-bound package needed for replay.
    #[must_use]
    pub const fn replay_package(&self) -> ContentHash {
        self.replay_package
    }

    /// Versioned, domain-separated identity of the complete projection.
    #[must_use]
    pub const fn content_hash(&self) -> ContentHash {
        self.content_hash
    }

    /// Recompute and compare the retained projection identity.
    #[must_use]
    pub fn validate_content_hash(&self) -> bool {
        self.content_hash == self.recompute_content_hash()
    }

    /// Deterministic reviewer-facing explanation. This is a projection of the
    /// retained artifacts, not a new compliance or scientific claim.
    #[must_use]
    pub fn render_explain(&self) -> String {
        let scalar = self.requirement.scalar();
        let mut output = format!(
            "decision-assessment-v{} identity={}\nquantity={} unit={} schema={} artifact={}\nrequirement={} effective-limit={} requirement-source-kind={} requirement-document={} requirement-version={} requirement-locator={} requirement-artifact={}@{}\nsafety-factor={} safety-factor-source-kind={} safety-factor-document={} safety-factor-version={} safety-factor-locator={} safety-factor-artifact={}@{}\ncontext={} hash={}\nreplay-package={}\n",
            DECISION_ASSESSMENT_IDENTITY_VERSION,
            self.content_hash,
            self.quantity.qoi(),
            self.quantity.unit(),
            self.quantity.schema(),
            self.quantity.artifact(),
            scalar.id(),
            scalar.limit(),
            self.requirement.source().kind().slug(),
            self.requirement.source().document(),
            self.requirement.source().version(),
            self.requirement.source().locator(),
            scalar.provenance().role(),
            scalar.provenance().digest(),
            self.requirement.safety_factor().value(),
            self.requirement.safety_factor_source().kind().slug(),
            self.requirement.safety_factor_source().document(),
            self.requirement.safety_factor_source().version(),
            self.requirement.safety_factor_source().locator(),
            self.requirement.safety_factor().policy().role(),
            self.requirement.safety_factor().policy().digest(),
            self.context.id(),
            self.context.hash(),
            self.replay_package,
        );
        output.push_str(&self.compliance.render_report());
        output.push_str(&self.attribution.render_report());
        if self.actions.is_empty() {
            output.push_str("next-actions=none-required-by-current-verdict\n");
        } else {
            output.push_str("next-actions:\n");
            for action in &self.actions {
                let _ = write!(
                    output,
                    "- unknown={} required-to-flip={} reason={} ",
                    action.unknown.name(),
                    action.required_magnitude,
                    action.reason
                );
                match &action.recommended_evidence {
                    RecommendedEvidence::Priced {
                        action,
                        action_kind,
                        decision_value,
                        cost,
                        value_per_cost,
                    } => {
                        let _ = writeln!(
                            output,
                            "action={} kind={} decision-value={} cost={} value-per-cost={}",
                            action,
                            action_kind_slug(*action_kind).unwrap_or("unsupported"),
                            decision_value,
                            cost,
                            value_per_cost
                        );
                    }
                    RecommendedEvidence::Unpriced { suggested_action } => {
                        let _ = writeln!(
                            output,
                            "unpriced-kind={}",
                            action_kind_slug(*suggested_action).unwrap_or("unsupported")
                        );
                    }
                }
            }
        }
        output
    }

    fn recompute_content_hash(&self) -> ContentHash {
        let mut encoder = Encoder::new();
        encoder.u32(DECISION_ASSESSMENT_IDENTITY_VERSION);
        encoder.string(self.quantity.qoi());
        encoder.string(self.quantity.unit());
        encoder.string(self.quantity.schema());
        encoder.hash(self.quantity.artifact());
        encode_requirement(&mut encoder, &self.requirement);
        encoder.u8(self.context.kind().canonical_wire_tag());
        encoder.string(self.context.id().as_str());
        encoder.hash(self.context.hash());
        encoder.hash(self.budget.content_id());
        encoder.string(&self.compliance.render_report());
        encoder.string(&self.attribution.render_report());
        encoder.u64(self.actions.len() as u64);
        for action in &self.actions {
            encode_action(&mut encoder, action);
        }
        encoder.hash(self.replay_package);
        hash_domain(DECISION_ASSESSMENT_IDENTITY_DOMAIN, &encoder.finish())
    }
}

#[allow(clippy::too_many_arguments)]
fn validate_bindings<Q>(
    quantity: &EvidenceRef<Q>,
    requirement: &DecisionRequirement,
    context: &ArtifactRef,
    compliance: &ComplianceVerdict,
    budget: &EngineeringUncertaintyBudget,
    attribution: &UncertaintyAttribution,
    actions: &[UnknownResolutionRecommendation],
    replay_package: &VerifiedPackage,
) -> Result<(), DecisionAssessmentError> {
    if context.kind() != ArtifactKind::ContextOfUse {
        return Err(DecisionAssessmentError::WrongContextKind);
    }
    require_nonzero_hash("context.hash", context.hash())?;
    if quantity.qoi() != requirement.scalar().qoi()
        || quantity.qoi() != budget.qoi()
        || quantity.unit() != requirement.scalar().unit()
        || quantity.unit() != budget.unit()
    {
        return Err(DecisionAssessmentError::QuantityMismatch);
    }
    if compliance.requirement() != requirement.scalar() {
        return Err(DecisionAssessmentError::RequirementMismatch);
    }
    if compliance.budget() != budget.content_id() {
        return Err(DecisionAssessmentError::BudgetMismatch);
    }
    if attribution.baseline() != compliance {
        return Err(DecisionAssessmentError::AttributionMismatch);
    }
    validate_actions(compliance, actions)?;
    if !replay_package.validate_binding() {
        return Err(DecisionAssessmentError::ReplayPackageInvalid);
    }
    require_nonzero_hash("replay-package.root", replay_package.report().merkle_root())?;
    Ok(())
}

fn validate_actions(
    compliance: &ComplianceVerdict,
    actions: &[UnknownResolutionRecommendation],
) -> Result<(), DecisionAssessmentError> {
    let unknowns = compliance.flipping_unknowns();
    if unknowns.len() != actions.len() {
        return Err(DecisionAssessmentError::ActionMismatch);
    }
    for (unknown, action) in unknowns.iter().zip(actions) {
        if unknown.kind() != action.unknown
            || unknown.reason() != action.reason
            || unknown.required_magnitude().to_bits() != action.required_magnitude.to_bits()
        {
            return Err(DecisionAssessmentError::ActionMismatch);
        }
        match &action.recommended_evidence {
            RecommendedEvidence::Priced {
                action,
                action_kind,
                decision_value,
                cost,
                value_per_cost,
            } => {
                validate_id("action.id", action)?;
                action_kind_slug(*action_kind)?;
                if !decision_value.is_finite()
                    || *decision_value <= 0.0
                    || !cost.is_finite()
                    || *cost < 0.0
                    || value_per_cost.is_nan()
                    || *value_per_cost <= 0.0
                {
                    return Err(DecisionAssessmentError::InvalidActionValue);
                }
                // Keep the stored ranking faithful to fs-voi's action-value
                // semantics: positive costs divide normally, while an
                // admissible zero cost has unbounded value per cost.
                let canonical_value_per_cost = if *cost > 0.0 {
                    *decision_value / *cost
                } else {
                    f64::INFINITY
                };
                if value_per_cost.to_bits() != canonical_value_per_cost.to_bits() {
                    return Err(DecisionAssessmentError::InvalidActionValue);
                }
            }
            RecommendedEvidence::Unpriced { suggested_action } => {
                action_kind_slug(*suggested_action)?;
                if *suggested_action != unknown.suggested_action() {
                    return Err(DecisionAssessmentError::ActionMismatch);
                }
            }
        }
    }
    Ok(())
}

fn encode_requirement(encoder: &mut Encoder, requirement: &DecisionRequirement) {
    let scalar = requirement.scalar();
    encoder.string(scalar.id());
    encoder.string(scalar.qoi());
    encoder.string(scalar.unit());
    encoder.u8(match scalar.relation() {
        RequirementRelation::AtMost => 1,
        RequirementRelation::AtLeast => 2,
    });
    encoder.u64(scalar.limit().to_bits());
    encoder.string(scalar.provenance().role());
    encoder.hash(scalar.provenance().digest());
    encode_authority(encoder, requirement.source());
    encoder.u64(requirement.safety_factor().value().to_bits());
    encoder.string(requirement.safety_factor().policy().role());
    encoder.hash(requirement.safety_factor().policy().digest());
    encode_authority(encoder, requirement.safety_factor_source());
}

fn encode_authority(encoder: &mut Encoder, authority: &RequirementAuthority) {
    encoder.u8(authority.kind().wire_tag());
    encoder.string(authority.document());
    encoder.string(authority.version());
    encoder.string(authority.locator());
}

fn encode_action(encoder: &mut Encoder, action: &UnknownResolutionRecommendation) {
    encoder.string(action.unknown.name());
    encoder.string(&action.reason);
    encoder.u64(action.required_magnitude.to_bits());
    match &action.recommended_evidence {
        RecommendedEvidence::Priced {
            action,
            action_kind,
            decision_value,
            cost,
            value_per_cost,
        } => {
            encoder.u8(1);
            encoder.string(action);
            encoder.string(action_kind_slug(*action_kind).expect("validated action kind"));
            encoder.u64(decision_value.to_bits());
            encoder.u64(cost.to_bits());
            encoder.u64(value_per_cost.to_bits());
        }
        RecommendedEvidence::Unpriced { suggested_action } => {
            encoder.u8(2);
            encoder.string(action_kind_slug(*suggested_action).expect("validated action kind"));
        }
    }
}

fn action_kind_slug(kind: ActionKind) -> Result<&'static str, DecisionAssessmentError> {
    match kind {
        ActionKind::SolverTolerance => Ok("solver-tolerance"),
        ActionKind::MeshRefinement => Ok("mesh-refinement"),
        ActionKind::TimeRefinement => Ok("time-refinement"),
        ActionKind::RepresentationEscalation => Ok("representation-escalation"),
        ActionKind::UqSamples => Ok("uq-samples"),
        ActionKind::MaterialCouponTest => Ok("material-coupon-test"),
        ActionKind::SensorCampaign => Ok("sensor-campaign"),
        ActionKind::Falsification => Ok("falsification"),
        ActionKind::StandardsObligation => Ok("standards-obligation"),
        ActionKind::Refusal => Ok("refusal"),
        _ => Err(DecisionAssessmentError::UnsupportedActionKind),
    }
}

fn admit_id(field: &'static str, value: String) -> Result<String, DecisionAssessmentError> {
    validate_id(field, &value)?;
    Ok(value)
}

fn admit_authority_field(
    field: &'static str,
    value: String,
) -> Result<String, DecisionAssessmentError> {
    if value.is_empty()
        || value.len() > MAX_AUTHORITY_FIELD_BYTES
        || value.trim() != value
        || value.chars().any(char::is_control)
    {
        return Err(DecisionAssessmentError::InvalidField { field });
    }
    Ok(value)
}

fn validate_id(field: &'static str, value: &str) -> Result<(), DecisionAssessmentError> {
    if value.is_empty()
        || value.len() > MAX_DECISION_ID_BYTES
        || !value.bytes().all(|byte| byte.is_ascii_graphic())
    {
        return Err(DecisionAssessmentError::InvalidField { field });
    }
    Ok(())
}

fn require_nonzero_hash(
    field: &'static str,
    hash: ContentHash,
) -> Result<(), DecisionAssessmentError> {
    if hash.as_bytes() == &[0; 32] {
        return Err(DecisionAssessmentError::MissingArtifact { field });
    }
    Ok(())
}

struct Encoder(Vec<u8>);

impl Encoder {
    fn new() -> Self {
        Self(Vec::new())
    }

    fn u8(&mut self, value: u8) {
        self.0.push(value);
    }

    fn u32(&mut self, value: u32) {
        self.0.extend_from_slice(&value.to_le_bytes());
    }

    fn u64(&mut self, value: u64) {
        self.0.extend_from_slice(&value.to_le_bytes());
    }

    fn string(&mut self, value: &str) {
        self.u64(value.len() as u64);
        self.0.extend_from_slice(value.as_bytes());
    }

    fn hash(&mut self, value: ContentHash) {
        self.0.extend_from_slice(value.as_bytes());
    }

    fn finish(self) -> Vec<u8> {
        self.0
    }
}

/// Typed refusal at the decision-projection boundary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DecisionAssessmentError {
    /// A bounded identity field was blank, non-ASCII-graphic, or oversized.
    InvalidField {
        /// Stable field name; rejected input bytes are not retained.
        field: &'static str,
    },
    /// A required content-addressed artifact had the all-zero sentinel.
    MissingArtifact {
        /// Stable role of the absent content identity.
        field: &'static str,
    },
    /// A safety factor was non-finite or below one.
    InvalidSafetyFactor,
    /// The supplied V&V reference was not a context-of-use artifact.
    WrongContextKind,
    /// Quantity, requirement, and budget QoI/unit semantics differ.
    QuantityMismatch,
    /// The compliance verdict was derived from another requirement.
    RequirementMismatch,
    /// The compliance verdict was derived from another uncertainty budget.
    BudgetMismatch,
    /// Attribution was derived from another compliance replay.
    AttributionMismatch,
    /// Actions do not correspond one-for-one with verdict-flipping unknowns.
    ActionMismatch,
    /// A priced action carried an invalid decision value, cost, or ratio.
    InvalidActionValue,
    /// A newer action taxonomy variant has no identity tag in this schema.
    UnsupportedActionKind,
    /// The replay package/report pair no longer validates its binding.
    ReplayPackageInvalid,
}

impl fmt::Display for DecisionAssessmentError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidField { field } => write!(f, "invalid bounded decision field {field}"),
            Self::MissingArtifact { field } => write!(f, "required artifact {field} is absent"),
            Self::InvalidSafetyFactor => {
                f.write_str("safety factor must be finite and at least one")
            }
            Self::WrongContextKind => f.write_str("context reference is not context-of-use"),
            Self::QuantityMismatch => {
                f.write_str("quantity, requirement, and budget QoI/unit semantics differ")
            }
            Self::RequirementMismatch => {
                f.write_str("compliance verdict names a different requirement")
            }
            Self::BudgetMismatch => {
                f.write_str("compliance verdict names a different uncertainty budget")
            }
            Self::AttributionMismatch => {
                f.write_str("attribution baseline differs from the compliance verdict")
            }
            Self::ActionMismatch => {
                f.write_str("next actions do not match the verdict-flipping unknowns")
            }
            Self::InvalidActionValue => {
                f.write_str("priced next action is outside fs-voi's eligible value/cost domain")
            }
            Self::UnsupportedActionKind => {
                f.write_str("action taxonomy variant is unsupported by identity schema v1")
            }
            Self::ReplayPackageInvalid => {
                f.write_str("replay package no longer validates its package/report binding")
            }
        }
    }
}

impl std::error::Error for DecisionAssessmentError {}
