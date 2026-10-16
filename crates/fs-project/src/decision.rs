//! Pure project-to-decision projection for requirement-bearing QoIs.
//!
//! This module binds admitted `.fsim` intent to `fs-session`'s generic L6
//! decision inputs. It does not solve physics, recompute compliance, resolve a
//! source document, or upgrade evidence authority.

use core::fmt;
use std::collections::BTreeSet;

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::{
    uncertainty::{
        ComplianceVerdict, EngineeringUncertaintyBudget, RequirementRelation, ScalarRequirement,
        UncertaintyArtifactRef, UncertaintyAttribution,
    },
    vv::{ArtifactId, ArtifactKind, ArtifactRef},
};
use fs_package::VerifiedPackage;
use fs_scenario::Violation;
use fs_session::{
    AppliedSafetyFactor, DecisionAssessment, DecisionAssessmentError, DecisionRequirement,
    EvidenceRef, RequirementAuthority, RequirementAuthorityKind,
};
use fs_voi::UnknownResolutionRecommendation;

use crate::{
    FSIM_VERSION,
    spec::{
        ConsequenceClass, DecisionGate, Metadata, ProjectSpec, RequirementDirection,
        RequirementSeverity, RequirementSource, RequirementSourceKind, ThermalLimit,
    },
};

/// Identity domain for one exact project requirement declaration.
pub const PROJECT_REQUIREMENT_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-project.decision-requirement.v1";
/// Identity domain for one exact safety-factor policy declaration.
pub const PROJECT_SAFETY_FACTOR_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-project.safety-factor-policy.v1";
/// Identity domain for one exact project context-of-use declaration.
pub const PROJECT_DECISION_CONTEXT_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-project.decision-context.v1";

/// Content-bound project context consumed by one decision assessment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProjectDecisionContext {
    project_name: String,
    created: String,
    context_of_use: String,
    intended_decision: String,
    decision_gate: DecisionGate,
    consequence: ConsequenceClass,
    artifact: ArtifactRef,
}

impl ProjectDecisionContext {
    /// Project name retained in the context identity.
    #[must_use]
    pub fn project_name(&self) -> &str {
        &self.project_name
    }

    /// Project creation date retained in the context identity.
    #[must_use]
    pub fn created(&self) -> &str {
        &self.created
    }

    /// Declared context in which the result may be relied on.
    #[must_use]
    pub fn context_of_use(&self) -> &str {
        &self.context_of_use
    }

    /// Engineering decision the result is intended to inform.
    #[must_use]
    pub fn intended_decision(&self) -> &str {
        &self.intended_decision
    }

    /// Declared decision gate.
    #[must_use]
    pub const fn decision_gate(&self) -> DecisionGate {
        self.decision_gate
    }

    /// Consequence framing for an unsupported decision.
    #[must_use]
    pub const fn consequence(&self) -> ConsequenceClass {
        self.consequence
    }

    /// Whether this context may retain an explicit indeterminate assessment.
    #[must_use]
    pub const fn permits_indeterminate(&self) -> bool {
        matches!(self.decision_gate, DecisionGate::ScopingEstimate)
            && !matches!(self.consequence, ConsequenceClass::SafetyCritical)
    }

    /// Exact content-bound `ContextOfUse` reference passed to `fs-session`.
    #[must_use]
    pub const fn artifact(&self) -> &ArtifactRef {
        &self.artifact
    }
}

/// Project-owned authority inputs for one requirement-bearing QoI.
#[derive(Debug, Clone, PartialEq)]
pub struct ProjectDecisionAuthority {
    requirement: DecisionRequirement,
    context: ProjectDecisionContext,
}

impl ProjectDecisionAuthority {
    /// Bind one already-admitted project metadata/requirement pair.
    ///
    /// This checks the decision-specific text, dimensions, numeric domains,
    /// and source lineages. [`project_decision_authorities`] remains the
    /// stronger entry point when a complete [`ProjectSpec`] is available.
    pub fn try_from_project_parts(
        metadata: &Metadata,
        limit: &ThermalLimit,
    ) -> Result<Self, ProjectDecisionError> {
        validate_decision_parts(metadata, limit)?;
        Ok(Self {
            requirement: decision_requirement(limit)?,
            context: project_context(metadata)?,
        })
    }

    /// Exact requirement and both source lineages.
    #[must_use]
    pub const fn requirement(&self) -> &DecisionRequirement {
        &self.requirement
    }

    /// Exact project context and gate.
    #[must_use]
    pub const fn context(&self) -> &ProjectDecisionContext {
        &self.context
    }

    /// Assemble already-created lower-layer artifacts, then apply the project
    /// context gate to the otherwise valid tri-state assessment.
    #[allow(clippy::too_many_arguments)]
    pub fn try_assemble<Q>(
        &self,
        quantity: EvidenceRef<Q>,
        compliance: ComplianceVerdict,
        budget: EngineeringUncertaintyBudget,
        attribution: UncertaintyAttribution,
        actions: Vec<UnknownResolutionRecommendation>,
        replay_package: &VerifiedPackage,
    ) -> Result<DecisionAssessment<Q>, ProjectDecisionError> {
        let assessment = DecisionAssessment::try_assemble(
            quantity,
            self.requirement.clone(),
            self.context.artifact.clone(),
            compliance,
            budget,
            attribution,
            actions,
            replay_package,
        )?;
        if matches!(
            assessment.compliance(),
            ComplianceVerdict::Indeterminate { .. }
        ) && !self.context.permits_indeterminate()
        {
            return Err(ProjectDecisionError::IndeterminateRefused {
                qoi: self.requirement.scalar().qoi().to_string(),
                decision_gate: self.context.decision_gate,
                consequence: self.context.consequence,
            });
        }
        Ok(assessment)
    }
}

/// Build deterministic authorities for every requirement-bearing project QoI.
///
/// The project must be fully admissible, and one QoI may carry at most one
/// scalar requirement in this thermal vertical.
pub fn project_decision_authorities(
    project: &ProjectSpec,
) -> Result<Vec<ProjectDecisionAuthority>, ProjectDecisionError> {
    let violations = project.validate();
    if !violations.is_empty() {
        return Err(ProjectDecisionError::InvalidProject { violations });
    }
    let metadata = project
        .metadata
        .as_ref()
        .ok_or(ProjectDecisionError::ProjectionInput {
            field: "metadata",
            detail: "validated project omitted metadata".to_string(),
        })?;
    let requirements =
        project
            .requirements
            .as_ref()
            .ok_or(ProjectDecisionError::ProjectionInput {
                field: "requirements",
                detail: "validated project omitted requirements".to_string(),
            })?;

    let mut seen = BTreeSet::new();
    for limit in requirements {
        if !seen.insert(limit.qoi.as_str()) {
            return Err(ProjectDecisionError::DuplicateRequirementQoi {
                qoi: limit.qoi.clone(),
            });
        }
    }

    let context = project_context(metadata)?;
    let mut authorities = requirements
        .iter()
        .map(|limit| {
            validate_decision_parts(metadata, limit)?;
            Ok(ProjectDecisionAuthority {
                requirement: decision_requirement(limit)?,
                context: context.clone(),
            })
        })
        .collect::<Result<Vec<_>, ProjectDecisionError>>()?;
    authorities.sort_by(|left, right| {
        left.requirement
            .scalar()
            .qoi()
            .cmp(right.requirement.scalar().qoi())
    });
    Ok(authorities)
}

/// Select exactly one project decision authority by QoI identity.
pub fn project_decision_authority(
    project: &ProjectSpec,
    qoi: &str,
) -> Result<ProjectDecisionAuthority, ProjectDecisionError> {
    project_decision_authorities(project)?
        .into_iter()
        .find(|authority| authority.requirement.scalar().qoi() == qoi)
        .ok_or_else(|| ProjectDecisionError::RequirementMissing {
            qoi: qoi.to_string(),
        })
}

fn validate_decision_parts(
    metadata: &Metadata,
    limit: &ThermalLimit,
) -> Result<(), ProjectDecisionError> {
    const MAX_PROJECT_DECISION_TEXT_BYTES: usize = 4096;
    for (field, value) in [
        ("metadata.name", metadata.name.as_str()),
        ("metadata.created", metadata.created.as_str()),
        ("metadata.context-of-use", metadata.context_of_use.as_str()),
        (
            "metadata.intended-decision",
            metadata.intended_decision.as_str(),
        ),
        ("requirement.qoi", limit.qoi.as_str()),
        ("requirement.class", limit.class.as_str()),
        ("requirement.region", limit.region.as_str()),
    ] {
        if value.is_empty()
            || value.len() > MAX_PROJECT_DECISION_TEXT_BYTES
            || value.trim() != value
            || value.chars().any(char::is_control)
        {
            return Err(ProjectDecisionError::ProjectionInput {
                field,
                detail: "value must be nonempty, bounded, trim-canonical, and control-free"
                    .to_string(),
            });
        }
    }
    if limit.limit.dims != crate::spec::dims::TEMPERATURE
        || limit.margin.dims != crate::spec::dims::TEMPERATURE
        || !limit.limit.value.is_finite()
        || !limit.margin.value.is_finite()
        || limit.margin.value < 0.0
    {
        return Err(ProjectDecisionError::ProjectionInput {
            field: "requirement.quantity",
            detail: "thermal limit and non-negative margin must be finite kelvin quantities"
                .to_string(),
        });
    }
    if !limit.safety_factor.factor.is_finite() || limit.safety_factor.factor < 1.0 {
        return Err(ProjectDecisionError::ProjectionInput {
            field: "requirement.safety-factor",
            detail: "factor must be finite and at least one".to_string(),
        });
    }
    Ok(())
}

fn decision_requirement(limit: &ThermalLimit) -> Result<DecisionRequirement, ProjectDecisionError> {
    let source = requirement_authority(&limit.source)?;
    let safety_factor_source = requirement_authority(&limit.safety_factor.source)?;
    let requirement_digest = hash_requirement(limit);
    let safety_factor_digest = hash_safety_factor(limit);
    let requirement_artifact =
        UncertaintyArtifactRef::new("project-requirement-source", requirement_digest).map_err(
            |error| ProjectDecisionError::ProjectionInput {
                field: "requirement.source",
                detail: error.to_string(),
            },
        )?;
    let policy_artifact =
        UncertaintyArtifactRef::new("project-safety-factor-policy", safety_factor_digest).map_err(
            |error| ProjectDecisionError::ProjectionInput {
                field: "requirement.safety-factor",
                detail: error.to_string(),
            },
        )?;
    let scalar_id = format!("project-requirement:{requirement_digest}");
    let relation = match limit.direction {
        RequirementDirection::AtMost => RequirementRelation::AtMost,
        RequirementDirection::AtLeast => RequirementRelation::AtLeast,
    };
    let scalar = ScalarRequirement::try_new(
        &scalar_id,
        &limit.qoi,
        "kelvin",
        relation,
        limit.limit.value,
        requirement_artifact,
    )
    .map_err(|error| ProjectDecisionError::ProjectionInput {
        field: "requirement.scalar",
        detail: error.to_string(),
    })?;
    let safety_factor = AppliedSafetyFactor::try_new(limit.safety_factor.factor, policy_artifact)?;
    DecisionRequirement::try_new(scalar, source, safety_factor, safety_factor_source)
        .map_err(ProjectDecisionError::from)
}

fn requirement_authority(
    source: &RequirementSource,
) -> Result<RequirementAuthority, ProjectDecisionError> {
    let kind = match source.kind {
        RequirementSourceKind::Standard => RequirementAuthorityKind::Standard,
        RequirementSourceKind::Datasheet => RequirementAuthorityKind::Datasheet,
        RequirementSourceKind::InternalPolicy => RequirementAuthorityKind::InternalPolicy,
        RequirementSourceKind::UserDeclaration => RequirementAuthorityKind::UserDeclaration,
    };
    RequirementAuthority::try_new(
        kind,
        source.document.clone(),
        source.version.clone(),
        source.locator.clone(),
    )
    .map_err(ProjectDecisionError::from)
}

fn project_context(metadata: &Metadata) -> Result<ProjectDecisionContext, ProjectDecisionError> {
    let hash = hash_context(metadata);
    let id = ArtifactId::try_new(&format!("project-context:{hash}")).map_err(|error| {
        ProjectDecisionError::ProjectionInput {
            field: "metadata.context-of-use",
            detail: error.to_string(),
        }
    })?;
    Ok(ProjectDecisionContext {
        project_name: metadata.name.clone(),
        created: metadata.created.clone(),
        context_of_use: metadata.context_of_use.clone(),
        intended_decision: metadata.intended_decision.clone(),
        decision_gate: metadata.decision_gate,
        consequence: metadata.consequence,
        artifact: ArtifactRef::new(ArtifactKind::ContextOfUse, id, hash),
    })
}

fn hash_requirement(limit: &ThermalLimit) -> ContentHash {
    let mut encoder = Encoder::new();
    encoder.u32(FSIM_VERSION);
    encoder.string(&limit.qoi);
    encoder.string(&limit.class);
    encoder.string(&limit.region);
    encoder.string(limit.direction.slug());
    encoder.qty(limit.limit);
    encoder.qty(limit.margin);
    encoder.source(&limit.source);
    encoder.u64(limit.safety_factor.factor.to_bits());
    encoder.source(&limit.safety_factor.source);
    encoder.string(limit.severity.slug());
    hash_domain(PROJECT_REQUIREMENT_IDENTITY_DOMAIN, &encoder.finish())
}

fn hash_safety_factor(limit: &ThermalLimit) -> ContentHash {
    let mut encoder = Encoder::new();
    encoder.u32(FSIM_VERSION);
    encoder.u64(limit.safety_factor.factor.to_bits());
    encoder.source(&limit.safety_factor.source);
    hash_domain(PROJECT_SAFETY_FACTOR_IDENTITY_DOMAIN, &encoder.finish())
}

fn hash_context(metadata: &Metadata) -> ContentHash {
    let mut encoder = Encoder::new();
    encoder.u32(FSIM_VERSION);
    encoder.string(&metadata.name);
    encoder.string(&metadata.created);
    encoder.string(&metadata.context_of_use);
    encoder.string(&metadata.intended_decision);
    encoder.string(metadata.decision_gate.slug());
    encoder.string(metadata.consequence.slug());
    hash_domain(PROJECT_DECISION_CONTEXT_IDENTITY_DOMAIN, &encoder.finish())
}

struct Encoder(Vec<u8>);

impl Encoder {
    fn new() -> Self {
        Self(Vec::new())
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

    fn qty(&mut self, value: fs_qty::QtyAny) {
        self.u64(value.value.to_bits());
        for exponent in value.dims.0 {
            self.0.extend_from_slice(&exponent.to_le_bytes());
        }
    }

    fn source(&mut self, source: &RequirementSource) {
        self.string(source.kind.slug());
        self.string(&source.document);
        self.string(&source.version);
        self.string(&source.locator);
    }

    fn finish(self) -> Vec<u8> {
        self.0
    }
}

impl RequirementDirection {
    const fn slug(self) -> &'static str {
        match self {
            Self::AtMost => "at-most",
            Self::AtLeast => "at-least",
        }
    }
}

impl RequirementSeverity {
    const fn slug(self) -> &'static str {
        match self {
            Self::ReliabilityDerating => "reliability-derating",
            Self::DamageLimit => "damage-limit",
            Self::SafetyCritical => "safety-critical",
        }
    }
}

/// Typed refusal while binding project intent to decision artifacts.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProjectDecisionError {
    /// The project failed its ordinary schema/dimension/reference admission.
    InvalidProject {
        /// Complete project validation findings.
        violations: Vec<Violation>,
    },
    /// No project requirement carries the requested QoI.
    RequirementMissing {
        /// Requested quantity identity.
        qoi: String,
    },
    /// More than one project requirement claims the same QoI.
    DuplicateRequirementQoi {
        /// Ambiguous quantity identity.
        qoi: String,
    },
    /// A project field could not enter the stricter decision identity domain.
    ProjectionInput {
        /// Stable failing field.
        field: &'static str,
        /// Lower-layer refusal detail.
        detail: String,
    },
    /// The project context forbids using an otherwise valid indeterminate result.
    IndeterminateRefused {
        /// Governed quantity identity.
        qoi: String,
        /// Gate that requires a determinate result.
        decision_gate: DecisionGate,
        /// Consequence framing participating in the refusal.
        consequence: ConsequenceClass,
    },
    /// The generic decision projection rejected inconsistent lower artifacts.
    Assessment(DecisionAssessmentError),
}

impl From<DecisionAssessmentError> for ProjectDecisionError {
    fn from(error: DecisionAssessmentError) -> Self {
        Self::Assessment(error)
    }
}

impl fmt::Display for ProjectDecisionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidProject { violations } => write!(
                f,
                "project is not admissible for decision assembly ({} finding(s), first={})",
                violations.len(),
                violations
                    .first()
                    .map_or("none", |violation| violation.code)
            ),
            Self::RequirementMissing { qoi } => {
                write!(f, "project has no requirement for QoI {qoi:?}")
            }
            Self::DuplicateRequirementQoi { qoi } => {
                write!(f, "project has multiple requirements for QoI {qoi:?}")
            }
            Self::ProjectionInput { field, detail } => {
                write!(f, "project decision field {field} was refused: {detail}")
            }
            Self::IndeterminateRefused {
                qoi,
                decision_gate,
                consequence,
            } => write!(
                f,
                "indeterminate QoI {qoi:?} is unusable for gate {} with consequence {}",
                decision_gate.slug(),
                consequence.slug()
            ),
            Self::Assessment(error) => write!(f, "decision assessment refused: {error}"),
        }
    }
}

impl std::error::Error for ProjectDecisionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Assessment(error) => Some(error),
            _ => None,
        }
    }
}
