//! Adversarial validation cases and honesty-first assessment.
//!
//! A registered case is a challenge specification, not evidence that the
//! challenge was executed. Retained corpus bindings and planned evidence are
//! distinct types, and the scorecard renders every unexecuted cell as
//! `NO-DATA`. A prediction outside its declared envelope is a false acceptance;
//! a refusal or demotion passes only when it names the case's expected
//! dominant uncertainty.

use core::fmt;
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::sync::LazyLock;

use fs_blake3::{ContentHash, hash_domain};

use crate::corpus::{EvidenceLevel, corpus};

/// Canonical schema for adversarial case, assessment, and scorecard identity.
pub const ADVERSARIAL_CASE_SCHEMA_VERSION: u32 = 1;
/// Maximum number of cases accepted in one registry.
pub const MAX_ADVERSARIAL_CASES: usize = 64;
/// Maximum UTF-8 bytes accepted in a case text field.
pub const MAX_ADVERSARIAL_TEXT_BYTES: usize = 1_024;

const REGISTRY_IDENTITY_DOMAIN: &str = "org.frankensim.fs-vvreg.adversarial-registry.v1";
const CASE_IDENTITY_DOMAIN: &str = "org.frankensim.fs-vvreg.adversarial-case.v1";
const ASSESSMENT_IDENTITY_DOMAIN: &str = "org.frankensim.fs-vvreg.adversarial-assessment.v1";

/// Source family intended to challenge the attacked assumption.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum AdversarialEvidenceBasis {
    /// Closed-form or manufactured reference.
    Analytic,
    /// Separately implemented solver or code.
    CrossCode,
    /// Controlled published experiment.
    ControlledExperiment,
    /// Instrumented in-house rig or operational fixture.
    InstrumentedRig,
}

impl AdversarialEvidenceBasis {
    /// Stable scorecard spelling.
    #[must_use]
    pub const fn slug(self) -> &'static str {
        match self {
            Self::Analytic => "analytic",
            Self::CrossCode => "cross-code",
            Self::ControlledExperiment => "controlled-experiment",
            Self::InstrumentedRig => "instrumented-rig",
        }
    }

    const fn tag(self) -> u8 {
        match self {
            Self::Analytic => 1,
            Self::CrossCode => 2,
            Self::ControlledExperiment => 3,
            Self::InstrumentedRig => 4,
        }
    }
}

/// Whether the challenge evidence is actually retained now.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum AdversarialEvidence {
    /// Exact dataset id in the authoritative [`crate::corpus::corpus`].
    Retained {
        /// Corpus dataset identity.
        dataset_id: &'static str,
    },
    /// Required evidence has not yet been retained.
    Planned {
        /// Bead governing acquisition or import.
        tracking_bead: &'static str,
        /// Stable no-data reason.
        reason: &'static str,
    },
}

impl AdversarialEvidence {
    /// Exact retained dataset id, if one exists.
    #[must_use]
    pub const fn dataset_id(self) -> Option<&'static str> {
        match self {
            Self::Retained { dataset_id } => Some(dataset_id),
            Self::Planned { .. } => None,
        }
    }

    /// Whether this challenge has retained corpus evidence.
    #[must_use]
    pub const fn is_retained(self) -> bool {
        matches!(self, Self::Retained { .. })
    }

    const fn tag(self) -> u8 {
        match self {
            Self::Retained { .. } => 1,
            Self::Planned { .. } => 2,
        }
    }
}

/// Assumption that a case is constructed to attack.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum AttackedAssumption {
    /// Correlation-rung flow remains attached and one-dimensional.
    AttachedFlow,
    /// Thermal contact resistance is zero or known exactly.
    KnownContactResistance,
    /// Convection dominates surface radiation.
    ConvectionDominates,
    /// A fan curve has one stable operating-point intersection.
    StableFanOperatingPoint,
    /// Vent leakage and enclosure boundary conditions are known.
    KnownVentLeakage,
    /// One nominal material-property value represents manufactured lots.
    FixedMaterialProperties,
    /// A lumped temperature is valid outside its small-Biot domain.
    LumpedTemperature,
    /// Forced-convection closure remains valid in buoyancy-driven reversal.
    ForcedConvectionClosure,
}

impl AttackedAssumption {
    /// Stable scorecard spelling.
    #[must_use]
    pub const fn slug(self) -> &'static str {
        match self {
            Self::AttachedFlow => "attached-flow",
            Self::KnownContactResistance => "known-contact-resistance",
            Self::ConvectionDominates => "convection-dominates",
            Self::StableFanOperatingPoint => "stable-fan-operating-point",
            Self::KnownVentLeakage => "known-vent-leakage",
            Self::FixedMaterialProperties => "fixed-material-properties",
            Self::LumpedTemperature => "lumped-temperature",
            Self::ForcedConvectionClosure => "forced-convection-closure",
        }
    }

    const fn tag(self) -> u8 {
        self as u8 + 1
    }
}

/// Dominant uncertainty that an honest refusal or demotion must name.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum DominantUncertainty {
    /// Separated flow and unresolved recirculation topology.
    FlowTopology,
    /// Interface-card value, state, or placement.
    ContactResistance,
    /// Radiation model, emissivity, or view factor.
    RadiationModel,
    /// Fan curve and operating-point stability.
    FanOperatingPoint,
    /// Leakage area or boundary-condition state.
    BoundaryCondition,
    /// Material state, lot, or property distribution.
    MaterialProperty,
    /// Spatial temperature gradients excluded by a lumped model.
    SpatialTemperature,
    /// Buoyancy/forced-flow regime transition.
    MixedConvectionRegime,
}

impl DominantUncertainty {
    /// Stable log spelling.
    #[must_use]
    pub const fn slug(self) -> &'static str {
        match self {
            Self::FlowTopology => "flow-topology",
            Self::ContactResistance => "contact-resistance",
            Self::RadiationModel => "radiation-model",
            Self::FanOperatingPoint => "fan-operating-point",
            Self::BoundaryCondition => "boundary-condition",
            Self::MaterialProperty => "material-property",
            Self::SpatialTemperature => "spatial-temperature",
            Self::MixedConvectionRegime => "mixed-convection-regime",
        }
    }

    const fn tag(self) -> u8 {
        self as u8 + 1
    }
}

/// One immutable adversarial challenge declaration.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct AdversarialCase {
    /// Stable lowercase slug.
    pub id: &'static str,
    /// Human-readable title.
    pub title: &'static str,
    /// Exact regime label used by reports.
    pub regime: &'static str,
    /// Assumption intentionally attacked.
    pub attacked_assumption: AttackedAssumption,
    /// Uncertainty an honest refusal/demotion must attribute.
    pub expected_dominant_uncertainty: DominantUncertainty,
    /// Intended evidence family.
    pub evidence_basis: AdversarialEvidenceBasis,
    /// Retained binding or explicit planned state.
    pub evidence: AdversarialEvidence,
    /// Public limitation statement for the scorecard.
    pub regime_limitation: &'static str,
}

/// Validated, deterministic adversarial challenge registry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdversarialRegistry {
    cases: Vec<AdversarialCase>,
    case_identities: Vec<ContentHash>,
    identity: ContentHash,
}

impl AdversarialRegistry {
    /// Validate, sort, and identity-bind challenge declarations.
    pub fn build(mut cases: Vec<AdversarialCase>) -> Result<Self, AdversarialRegistryError> {
        if cases.len() > MAX_ADVERSARIAL_CASES {
            return Err(AdversarialRegistryError::ResourceLimit {
                limit: MAX_ADVERSARIAL_CASES,
                observed: cases.len(),
            });
        }
        for case in &cases {
            validate_case(case)?;
        }
        cases.sort_by_key(|case| case.id);
        for pair in cases.windows(2) {
            if pair[0].id == pair[1].id {
                return Err(AdversarialRegistryError::DuplicateCaseId {
                    id: pair[0].id.to_string(),
                });
            }
        }
        let case_identities = cases.iter().map(case_identity).collect::<Vec<_>>();
        let identity = registry_identity(&case_identities);
        Ok(Self {
            cases,
            case_identities,
            identity,
        })
    }

    /// Canonically ordered cases.
    #[must_use]
    pub fn cases(&self) -> &[AdversarialCase] {
        &self.cases
    }

    /// Registry identity covering every complete case declaration.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Identity of one named case.
    #[must_use]
    pub fn case_identity(&self, id: &str) -> Option<ContentHash> {
        self.cases
            .binary_search_by_key(&id, |case| case.id)
            .ok()
            .map(|index| self.case_identities[index])
    }

    /// Evaluate one prediction/refusal/demotion against the registered honesty
    /// rule. A result is a diagnostic receipt, never evidence authority.
    pub fn assess(
        &self,
        case_id: &str,
        outcome: AdversarialOutcome,
    ) -> Result<AdversarialAssessment, AdversarialAssessmentRefusal> {
        let index = self
            .cases
            .binary_search_by_key(&case_id, |case| case.id)
            .map_err(|_| AdversarialAssessmentRefusal::UnknownCase {
                id: case_id.to_string(),
            })?;
        validate_outcome(outcome)?;
        let case = self.cases[index];
        let (verdict, false_acceptance) = match outcome {
            AdversarialOutcome::Prediction {
                absolute_error,
                allowed_error,
                ..
            } => {
                if absolute_error <= allowed_error {
                    (HonestyVerdict::Pass, false)
                } else {
                    (HonestyVerdict::Fail, true)
                }
            }
            AdversarialOutcome::Refused { dominant } | AdversarialOutcome::Demoted { dominant } => {
                if dominant == case.expected_dominant_uncertainty {
                    (HonestyVerdict::Pass, false)
                } else {
                    (HonestyVerdict::Fail, false)
                }
            }
            AdversarialOutcome::NoData => (HonestyVerdict::NoData, false),
        };
        let identity = assessment_identity(
            self.identity,
            self.case_identities[index],
            outcome,
            verdict,
            false_acceptance,
        );
        Ok(AdversarialAssessment {
            registry_identity: self.identity,
            case_identity: self.case_identities[index],
            case_id: case.id,
            outcome,
            verdict,
            false_acceptance,
            identity,
        })
    }

    /// Render the future public scorecard's deterministic adversarial
    /// regime-limitation section.
    ///
    /// Missing assessments and explicitly planned evidence remain `NO-DATA`;
    /// neither is coerced to a zero error or a passing result.
    pub fn render_regime_limitations(
        &self,
        assessments: &[AdversarialAssessment],
    ) -> Result<String, AdversarialScorecardError> {
        let mut by_case = BTreeMap::new();
        for assessment in assessments {
            if assessment.registry_identity != self.identity {
                return Err(AdversarialScorecardError::ForeignAssessment {
                    case_id: assessment.case_id.to_string(),
                });
            }
            if by_case.insert(assessment.case_id, assessment).is_some() {
                return Err(AdversarialScorecardError::DuplicateAssessment {
                    case_id: assessment.case_id.to_string(),
                });
            }
        }

        let false_acceptances = assessments
            .iter()
            .filter(|assessment| assessment.false_acceptance)
            .count();
        let mut report = format!(
            "## Adversarial regime limitations\n\nschema: {ADVERSARIAL_CASE_SCHEMA_VERSION}\nregistry: {}\nfalse_acceptance_count: {false_acceptances}\n\n",
            hex(self.identity)
        );
        report.push_str(
            "| case | regime | attacked assumption | evidence | assessment | dominant uncertainty | limitation |\n",
        );
        report.push_str("| --- | --- | --- | --- | --- | --- | --- |\n");
        for case in &self.cases {
            let assessment = by_case.get(case.id).copied();
            let evidence = match case.evidence {
                AdversarialEvidence::Retained { dataset_id } => {
                    format!("retained:{dataset_id}")
                }
                AdversarialEvidence::Planned { tracking_bead, .. } => {
                    format!("NO-DATA:{tracking_bead}")
                }
            };
            let status = assessment.map_or("NO-DATA", |value| value.verdict.slug());
            let dominant = assessment
                .and_then(|value| value.outcome.dominant())
                .map_or("NO-DATA", DominantUncertainty::slug);
            let _ = writeln!(
                report,
                "| {} | {} | {} | {} | {} | {} | {} |",
                case.id,
                case.regime,
                case.attacked_assumption.slug(),
                evidence,
                status,
                dominant,
                case.regime_limitation
            );
        }
        Ok(report)
    }
}

/// System response to one adversarial challenge.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AdversarialOutcome {
    /// The system accepted and reported a prediction.
    Prediction {
        /// Absolute error against the challenge reference.
        absolute_error: f64,
        /// Inclusive challenge envelope.
        allowed_error: f64,
        /// Reported dominant uncertainty.
        dominant: DominantUncertainty,
    },
    /// The system refused to predict.
    Refused {
        /// Reported dominant uncertainty.
        dominant: DominantUncertainty,
    },
    /// The system demoted its claim.
    Demoted {
        /// Reported dominant uncertainty.
        dominant: DominantUncertainty,
    },
    /// No challenge result exists.
    NoData,
}

impl AdversarialOutcome {
    /// Reported dominant uncertainty, if the outcome has one.
    #[must_use]
    pub const fn dominant(self) -> Option<DominantUncertainty> {
        match self {
            Self::Prediction { dominant, .. }
            | Self::Refused { dominant }
            | Self::Demoted { dominant } => Some(dominant),
            Self::NoData => None,
        }
    }

    const fn tag(self) -> u8 {
        match self {
            Self::Prediction { .. } => 1,
            Self::Refused { .. } => 2,
            Self::Demoted { .. } => 3,
            Self::NoData => 4,
        }
    }
}

/// Honesty result, kept separate from solver accuracy or evidence color.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HonestyVerdict {
    /// In-envelope prediction or correctly attributed refusal/demotion.
    Pass,
    /// Out-of-envelope accepted prediction or misattributed refusal/demotion.
    Fail,
    /// Challenge was not executed.
    NoData,
}

impl HonestyVerdict {
    /// Stable scorecard spelling.
    #[must_use]
    pub const fn slug(self) -> &'static str {
        match self {
            Self::Pass => "PASS",
            Self::Fail => "FAIL",
            Self::NoData => "NO-DATA",
        }
    }
}

/// Sealed diagnostic record returned by [`AdversarialRegistry::assess`].
#[derive(Debug, Clone, PartialEq)]
pub struct AdversarialAssessment {
    registry_identity: ContentHash,
    case_identity: ContentHash,
    case_id: &'static str,
    outcome: AdversarialOutcome,
    verdict: HonestyVerdict,
    false_acceptance: bool,
    identity: ContentHash,
}

impl AdversarialAssessment {
    /// Exact case id.
    #[must_use]
    pub const fn case_id(&self) -> &'static str {
        self.case_id
    }

    /// Supplied system outcome.
    #[must_use]
    pub const fn outcome(&self) -> AdversarialOutcome {
        self.outcome
    }

    /// Honesty verdict.
    #[must_use]
    pub const fn verdict(&self) -> HonestyVerdict {
        self.verdict
    }

    /// Whether an accepted prediction was outside its challenge envelope.
    #[must_use]
    pub const fn is_false_acceptance(&self) -> bool {
        self.false_acceptance
    }

    /// Canonical assessment identity.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Deterministic bounded audit row.
    #[must_use]
    pub fn render_log(&self) -> String {
        format!(
            "schema={} case={} case_identity={} outcome={} dominant={} verdict={} false_acceptance={} assessment={}",
            ADVERSARIAL_CASE_SCHEMA_VERSION,
            self.case_id,
            hex(self.case_identity),
            outcome_slug(self.outcome),
            self.outcome
                .dominant()
                .map_or("none", DominantUncertainty::slug),
            self.verdict.slug(),
            self.false_acceptance,
            hex(self.identity)
        )
    }
}

/// Registry-construction refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AdversarialRegistryError {
    /// A bounded field is invalid.
    InvalidField {
        /// Stable field name.
        field: &'static str,
        /// Stable reason.
        reason: &'static str,
    },
    /// Too many cases were supplied.
    ResourceLimit {
        /// Maximum accepted count.
        limit: usize,
        /// Supplied count.
        observed: usize,
    },
    /// Two cases use the same id.
    DuplicateCaseId {
        /// Duplicate id.
        id: String,
    },
    /// A retained corpus id does not exist.
    UnknownDataset {
        /// Case id.
        case_id: String,
        /// Missing dataset id.
        dataset_id: String,
    },
    /// Retained dataset coordinates do not match the declared basis.
    EvidenceBasisMismatch {
        /// Case id.
        case_id: String,
        /// Declared basis.
        basis: AdversarialEvidenceBasis,
        /// Dataset's actual legacy coordinate tag.
        level: EvidenceLevel,
    },
}

impl fmt::Display for AdversarialRegistryError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidField { field, reason } => {
                write!(formatter, "adversarial case field `{field}` {reason}")
            }
            Self::ResourceLimit { limit, observed } => write!(
                formatter,
                "adversarial case count {observed} exceeds limit {limit}"
            ),
            Self::DuplicateCaseId { id } => write!(formatter, "duplicate adversarial case `{id}`"),
            Self::UnknownDataset {
                case_id,
                dataset_id,
            } => write!(
                formatter,
                "adversarial case `{case_id}` names unknown dataset `{dataset_id}`"
            ),
            Self::EvidenceBasisMismatch {
                case_id,
                basis,
                level,
            } => write!(
                formatter,
                "adversarial case `{case_id}` declares {} but retained dataset has legacy coordinate {level:?}",
                basis.slug()
            ),
        }
    }
}

impl std::error::Error for AdversarialRegistryError {}

/// Assessment-input refusal. Failed honesty assessments are successful
/// evaluations and therefore are not represented by this type.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AdversarialAssessmentRefusal {
    /// Case id was not registered.
    UnknownCase {
        /// Supplied id.
        id: String,
    },
    /// Prediction arithmetic was invalid.
    InvalidPrediction {
        /// Stable reason.
        reason: &'static str,
    },
}

impl fmt::Display for AdversarialAssessmentRefusal {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownCase { id } => write!(formatter, "unknown adversarial case `{id}`"),
            Self::InvalidPrediction { reason } => {
                write!(formatter, "invalid adversarial prediction: {reason}")
            }
        }
    }
}

impl std::error::Error for AdversarialAssessmentRefusal {}

/// Scorecard assembly refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AdversarialScorecardError {
    /// An assessment belongs to another registry identity.
    ForeignAssessment {
        /// Assessment case id.
        case_id: String,
    },
    /// More than one assessment was supplied for a case.
    DuplicateAssessment {
        /// Duplicate case id.
        case_id: String,
    },
}

impl fmt::Display for AdversarialScorecardError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ForeignAssessment { case_id } => {
                write!(
                    formatter,
                    "assessment for `{case_id}` belongs to another registry"
                )
            }
            Self::DuplicateAssessment { case_id } => {
                write!(formatter, "duplicate assessment for `{case_id}`")
            }
        }
    }
}

impl std::error::Error for AdversarialScorecardError {}

/// Authoritative built-in adversarial thermal challenge registry.
#[must_use]
pub fn adversarial_registry() -> &'static AdversarialRegistry {
    static REGISTRY: LazyLock<AdversarialRegistry> = LazyLock::new(|| {
        AdversarialRegistry::build(SEEDED_CASES.to_vec())
            .expect("built-in adversarial cases satisfy their contract")
    });
    &REGISTRY
}

static SEEDED_CASES: [AdversarialCase; 8] = [
    AdversarialCase {
        id: "recirculation-behind-strip-fins",
        title: "Separated forced flow behind inline strip fins",
        regime: "forced-air-strip-fin-re-810-3800",
        attacked_assumption: AttackedAssumption::AttachedFlow,
        expected_dominant_uncertainty: DominantUncertainty::FlowTopology,
        evidence_basis: AdversarialEvidenceBasis::ControlledExperiment,
        evidence: AdversarialEvidence::Retained {
            dataset_id: "pires-fonseca-2024-flat-strip-fins",
        },
        regime_limitation: "Attached-flow correlations are not validated inside separated fin wakes; predict within the retained experimental envelope or demote for flow-topology uncertainty.",
    },
    AdversarialCase {
        id: "contact-dominated-two-layer-stack",
        title: "Two-layer stack dominated by interface resistance",
        regime: "series-contact-resistance-dominant",
        attacked_assumption: AttackedAssumption::KnownContactResistance,
        expected_dominant_uncertainty: DominantUncertainty::ContactResistance,
        evidence_basis: AdversarialEvidenceBasis::Analytic,
        evidence: AdversarialEvidence::Retained {
            dataset_id: "thermal-a-contact-series",
        },
        regime_limitation: "Perfect-contact or silently defaulted interface laws are outside scope; the interface card and its uncertainty must dominate the reported budget when appropriate.",
    },
    AdversarialCase {
        id: "radiation-dominated-low-flow-enclosure",
        title: "Low-flow enclosure with radiative exchange dominance",
        regime: "low-flow-high-emissivity-enclosure",
        attacked_assumption: AttackedAssumption::ConvectionDominates,
        expected_dominant_uncertainty: DominantUncertainty::RadiationModel,
        evidence_basis: AdversarialEvidenceBasis::Analytic,
        evidence: AdversarialEvidence::Retained {
            dataset_id: "thermal-a-parallel-plate-view-factor",
        },
        regime_limitation: "Convection-only predictions are outside scope when radiative exchange is material; emissivity, view-factor, and nonlinear-radiation discrepancy must remain visible.",
    },
    AdversarialCase {
        id: "fan-stall-multiple-operating-points",
        title: "Fan-stall curve with multiple or unstable intersections",
        regime: "fan-stall-negative-slope",
        attacked_assumption: AttackedAssumption::StableFanOperatingPoint,
        expected_dominant_uncertainty: DominantUncertainty::FanOperatingPoint,
        evidence_basis: AdversarialEvidenceBasis::InstrumentedRig,
        evidence: AdversarialEvidence::Planned {
            tracking_bead: "frankensim-extreal-program-f85xj.4.5",
            reason: "the instrumented fan-stall sweep and raw pressure-flow histories are not retained",
        },
        regime_limitation: "A single steady fan-curve intersection is not admissible through stall or hysteresis; unresolved operating-point multiplicity requires refusal or demotion.",
    },
    AdversarialCase {
        id: "uncertain-blockable-vent-leakage",
        title: "Blockable enclosure vents with uncertain leakage",
        regime: "sealed-to-leaky-enclosure-transition",
        attacked_assumption: AttackedAssumption::KnownVentLeakage,
        expected_dominant_uncertainty: DominantUncertainty::BoundaryCondition,
        evidence_basis: AdversarialEvidenceBasis::InstrumentedRig,
        evidence: AdversarialEvidence::Planned {
            tracking_bead: "frankensim-extreal-program-f85xj.4.5",
            reason: "as-built leakage metrology and retained vent-state sweeps do not yet exist",
        },
        regime_limitation: "Nominal vent geometry cannot stand in for as-built leakage; unknown leakage area must remain a boundary-condition uncertainty.",
    },
    AdversarialCase {
        id: "material-lot-property-variability",
        title: "Manufactured lots with strong thermal-property variability",
        regime: "multi-lot-as-manufactured-material-state",
        attacked_assumption: AttackedAssumption::FixedMaterialProperties,
        expected_dominant_uncertainty: DominantUncertainty::MaterialProperty,
        evidence_basis: AdversarialEvidenceBasis::InstrumentedRig,
        evidence: AdversarialEvidence::Planned {
            tracking_bead: "frankensim-extreal-program-f85xj.4.5",
            reason: "lot-resolved coupons, calibration records, and raw property measurements are not retained",
        },
        regime_limitation: "Nominal handbook values do not establish an as-manufactured lot; unresolved lot variability must remain in the material-property budget.",
    },
    AdversarialCase {
        id: "biot-extremes-lumped-breakdown",
        title: "Low- and high-Biot extremes around lumped-model validity",
        regime: "biot-number-validity-boundary",
        attacked_assumption: AttackedAssumption::LumpedTemperature,
        expected_dominant_uncertainty: DominantUncertainty::SpatialTemperature,
        evidence_basis: AdversarialEvidenceBasis::Analytic,
        evidence: AdversarialEvidence::Retained {
            dataset_id: "thermal-a-lumped-transient",
        },
        regime_limitation: "Lumped-capacitance predictions are restricted to their declared small-Biot domain; high-Biot cases require spatial resolution or an explicit demotion.",
    },
    AdversarialCase {
        id: "natural-convection-cavity-reversal",
        title: "Natural-convection cavity with buoyant flow reversal",
        regime: "buoyancy-dominated-enclosure-cavity",
        attacked_assumption: AttackedAssumption::ForcedConvectionClosure,
        expected_dominant_uncertainty: DominantUncertainty::MixedConvectionRegime,
        evidence_basis: AdversarialEvidenceBasis::CrossCode,
        evidence: AdversarialEvidence::Planned {
            tracking_bead: "frankensim-extreal-program-f85xj.4.3",
            reason: "the independently generated cavity deck, field, probes, and extraction receipt are not retained",
        },
        regime_limitation: "Forced-convection correlations do not cover buoyancy-driven reversal; use an independently retained cavity reference or refuse the regime.",
    },
];

fn validate_case(case: &AdversarialCase) -> Result<(), AdversarialRegistryError> {
    validate_slug(case.id)?;
    validate_text("title", case.title)?;
    validate_text("regime", case.regime)?;
    validate_text("regime_limitation", case.regime_limitation)?;
    match case.evidence {
        AdversarialEvidence::Retained { dataset_id } => {
            validate_slug(dataset_id)?;
            let dataset = corpus().dataset(dataset_id).ok_or_else(|| {
                AdversarialRegistryError::UnknownDataset {
                    case_id: case.id.to_string(),
                    dataset_id: dataset_id.to_string(),
                }
            })?;
            if !basis_matches_level(case.evidence_basis, dataset.evidence_level()) {
                return Err(AdversarialRegistryError::EvidenceBasisMismatch {
                    case_id: case.id.to_string(),
                    basis: case.evidence_basis,
                    level: dataset.evidence_level(),
                });
            }
        }
        AdversarialEvidence::Planned {
            tracking_bead,
            reason,
        } => {
            validate_text("tracking_bead", tracking_bead)?;
            validate_text("planned_reason", reason)?;
        }
    }
    Ok(())
}

fn validate_slug(value: &str) -> Result<(), AdversarialRegistryError> {
    if value.is_empty()
        || value.len() > MAX_ADVERSARIAL_TEXT_BYTES
        || !value
            .bytes()
            .all(|byte| byte.is_ascii_lowercase() || byte.is_ascii_digit() || byte == b'-')
    {
        return Err(AdversarialRegistryError::InvalidField {
            field: "id",
            reason: "must be a bounded lowercase slug",
        });
    }
    Ok(())
}

fn validate_text(field: &'static str, value: &str) -> Result<(), AdversarialRegistryError> {
    if value.trim().is_empty() {
        return Err(AdversarialRegistryError::InvalidField {
            field,
            reason: "is blank",
        });
    }
    if value.len() > MAX_ADVERSARIAL_TEXT_BYTES {
        return Err(AdversarialRegistryError::InvalidField {
            field,
            reason: "exceeds the byte limit",
        });
    }
    if value.chars().any(char::is_control) {
        return Err(AdversarialRegistryError::InvalidField {
            field,
            reason: "contains a control character",
        });
    }
    Ok(())
}

const fn basis_matches_level(basis: AdversarialEvidenceBasis, level: EvidenceLevel) -> bool {
    matches!(
        (basis, level),
        (AdversarialEvidenceBasis::Analytic, EvidenceLevel::Analytic)
            | (
                AdversarialEvidenceBasis::CrossCode,
                EvidenceLevel::CrossCode
            )
            | (
                AdversarialEvidenceBasis::ControlledExperiment,
                EvidenceLevel::PublishedExperiment | EvidenceLevel::Blind
            )
            | (
                AdversarialEvidenceBasis::InstrumentedRig,
                EvidenceLevel::PublishedExperiment | EvidenceLevel::Blind | EvidenceLevel::Field
            )
    )
}

fn validate_outcome(outcome: AdversarialOutcome) -> Result<(), AdversarialAssessmentRefusal> {
    if let AdversarialOutcome::Prediction {
        absolute_error,
        allowed_error,
        ..
    } = outcome
    {
        if !absolute_error.is_finite() || !allowed_error.is_finite() {
            return Err(AdversarialAssessmentRefusal::InvalidPrediction {
                reason: "error and envelope must be finite",
            });
        }
        if absolute_error < 0.0 || allowed_error < 0.0 {
            return Err(AdversarialAssessmentRefusal::InvalidPrediction {
                reason: "error and envelope must be nonnegative",
            });
        }
    }
    Ok(())
}

fn case_identity(case: &AdversarialCase) -> ContentHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&ADVERSARIAL_CASE_SCHEMA_VERSION.to_le_bytes());
    push_text(&mut bytes, case.id);
    push_text(&mut bytes, case.title);
    push_text(&mut bytes, case.regime);
    bytes.push(case.attacked_assumption.tag());
    bytes.push(case.expected_dominant_uncertainty.tag());
    bytes.push(case.evidence_basis.tag());
    bytes.push(case.evidence.tag());
    match case.evidence {
        AdversarialEvidence::Retained { dataset_id } => push_text(&mut bytes, dataset_id),
        AdversarialEvidence::Planned {
            tracking_bead,
            reason,
        } => {
            push_text(&mut bytes, tracking_bead);
            push_text(&mut bytes, reason);
        }
    }
    push_text(&mut bytes, case.regime_limitation);
    hash_domain(CASE_IDENTITY_DOMAIN, &bytes)
}

fn registry_identity(case_identities: &[ContentHash]) -> ContentHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&ADVERSARIAL_CASE_SCHEMA_VERSION.to_le_bytes());
    push_len(&mut bytes, case_identities.len());
    for identity in case_identities {
        bytes.extend_from_slice(&identity.0);
    }
    hash_domain(REGISTRY_IDENTITY_DOMAIN, &bytes)
}

fn assessment_identity(
    registry_identity: ContentHash,
    case_identity: ContentHash,
    outcome: AdversarialOutcome,
    verdict: HonestyVerdict,
    false_acceptance: bool,
) -> ContentHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&ADVERSARIAL_CASE_SCHEMA_VERSION.to_le_bytes());
    bytes.extend_from_slice(&registry_identity.0);
    bytes.extend_from_slice(&case_identity.0);
    bytes.push(outcome.tag());
    if let AdversarialOutcome::Prediction {
        absolute_error,
        allowed_error,
        ..
    } = outcome
    {
        bytes.extend_from_slice(&absolute_error.to_bits().to_le_bytes());
        bytes.extend_from_slice(&allowed_error.to_bits().to_le_bytes());
    }
    if let Some(dominant) = outcome.dominant() {
        bytes.push(dominant.tag());
    }
    bytes.push(match verdict {
        HonestyVerdict::Pass => 1,
        HonestyVerdict::Fail => 2,
        HonestyVerdict::NoData => 3,
    });
    bytes.push(u8::from(false_acceptance));
    hash_domain(ASSESSMENT_IDENTITY_DOMAIN, &bytes)
}

fn outcome_slug(outcome: AdversarialOutcome) -> &'static str {
    match outcome {
        AdversarialOutcome::Prediction { .. } => "prediction",
        AdversarialOutcome::Refused { .. } => "refused",
        AdversarialOutcome::Demoted { .. } => "demoted",
        AdversarialOutcome::NoData => "no-data",
    }
}

fn push_len(bytes: &mut Vec<u8>, len: usize) {
    bytes.extend_from_slice(&u64::try_from(len).unwrap_or(u64::MAX).to_le_bytes());
}

fn push_text(bytes: &mut Vec<u8>, value: &str) {
    push_len(bytes, value.len());
    bytes.extend_from_slice(value.as_bytes());
}

fn hex(hash: ContentHash) -> String {
    const DIGITS: &[u8; 16] = b"0123456789abcdef";
    let mut output = String::with_capacity(64);
    for byte in hash.0 {
        output.push(char::from(DIGITS[usize::from(byte >> 4)]));
        output.push(char::from(DIGITS[usize::from(byte & 0x0f)]));
    }
    output
}
