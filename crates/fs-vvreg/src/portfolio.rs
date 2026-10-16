//! Claim-scoped external-evidence portfolios.
//!
//! The axes in this module are coordinates, never rungs in a quality ladder.
//! In particular, field monitoring cannot substitute for a controlled
//! experiment, and repeating evidence on one axis cannot manufacture a
//! missing axis.

use core::fmt;

use fs_blake3::{ContentHash, hash_domain};

use crate::corpus::EvidenceLevel;

/// Canonical schema for portfolio and admission identities.
pub const EVIDENCE_PORTFOLIO_SCHEMA_VERSION: u32 = 1;
/// Maximum observations admitted into one bounded portfolio.
pub const MAX_PORTFOLIO_OBSERVATIONS: usize = 4_096;
/// Maximum UTF-8 bytes in a QoI or regime identifier.
pub const MAX_PORTFOLIO_ID_BYTES: usize = 256;

const PORTFOLIO_IDENTITY_DOMAIN: &str = "org.frankensim.fs-vvreg.evidence-portfolio.v1";
const ADMISSION_IDENTITY_DOMAIN: &str = "org.frankensim.fs-vvreg.evidence-portfolio-admission.v1";

/// Independent external-evidence coordinates.
///
/// Deriving `Ord` gives the canonical wire/report order. It does not define an
/// epistemic ranking.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EvidenceAxis {
    /// Analytic, manufactured, or independently checked numerical behavior.
    NumericalVerification,
    /// Agreement with a separately implemented code or solver.
    CrossCodeAgreement,
    /// Comparison with a controlled physical experiment.
    ControlledExperimentalValidation,
    /// A prediction evaluated after a preregistered blind split was frozen.
    BlindPredictiveValidation,
    /// Monitoring under field or operational conditions.
    FieldMonitoring,
    /// Evidence that a claim transfers across declared regimes.
    TransferabilityAcrossRegimes,
    /// Reproduction by an independent team or implementation lineage.
    IndependentReproduction,
}

impl EvidenceAxis {
    /// Canonical axis order for scorecards and identities.
    pub const ALL: [Self; 7] = [
        Self::NumericalVerification,
        Self::CrossCodeAgreement,
        Self::ControlledExperimentalValidation,
        Self::BlindPredictiveValidation,
        Self::FieldMonitoring,
        Self::TransferabilityAcrossRegimes,
        Self::IndependentReproduction,
    ];

    /// Stable machine-readable axis name.
    #[must_use]
    pub const fn slug(self) -> &'static str {
        match self {
            Self::NumericalVerification => "numerical-verification",
            Self::CrossCodeAgreement => "cross-code-agreement",
            Self::ControlledExperimentalValidation => "controlled-experimental-validation",
            Self::BlindPredictiveValidation => "blind-predictive-validation",
            Self::FieldMonitoring => "field-monitoring",
            Self::TransferabilityAcrossRegimes => "transferability-across-regimes",
            Self::IndependentReproduction => "independent-reproduction",
        }
    }

    /// Stable zero-based position in [`Self::ALL`].
    #[must_use]
    pub const fn index(self) -> usize {
        match self {
            Self::NumericalVerification => 0,
            Self::CrossCodeAgreement => 1,
            Self::ControlledExperimentalValidation => 2,
            Self::BlindPredictiveValidation => 3,
            Self::FieldMonitoring => 4,
            Self::TransferabilityAcrossRegimes => 5,
            Self::IndependentReproduction => 6,
        }
    }

    const fn tag(self) -> u8 {
        match self {
            Self::NumericalVerification => 1,
            Self::CrossCodeAgreement => 2,
            Self::ControlledExperimentalValidation => 3,
            Self::BlindPredictiveValidation => 4,
            Self::FieldMonitoring => 5,
            Self::TransferabilityAcrossRegimes => 6,
            Self::IndependentReproduction => 7,
        }
    }
}

/// Legacy A-E corpus tags interpreted as one or more portfolio coordinates.
///
/// Level D is both a controlled experiment and a blind prediction. Level E is
/// field monitoring only; it is deliberately not an experimental-validation
/// coordinate.
#[must_use]
pub const fn axes_for_level(level: EvidenceLevel) -> &'static [EvidenceAxis] {
    match level {
        EvidenceLevel::Analytic => &[EvidenceAxis::NumericalVerification],
        EvidenceLevel::CrossCode => &[EvidenceAxis::CrossCodeAgreement],
        EvidenceLevel::PublishedExperiment => &[EvidenceAxis::ControlledExperimentalValidation],
        EvidenceLevel::Blind => &[
            EvidenceAxis::ControlledExperimentalValidation,
            EvidenceAxis::BlindPredictiveValidation,
        ],
        EvidenceLevel::Field => &[EvidenceAxis::FieldMonitoring],
    }
}

/// Claim classes with explicit, non-substitutable axis requirements.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PortfolioClaimClass {
    /// A numerical implementation behaves as declared.
    NumericallyVerified,
    /// A numerical result agrees with a distinct implementation.
    CrossCodeConsistent,
    /// A physical prediction is validated in the named regime.
    ValidatedPrediction,
    /// A physical prediction is validated and evaluated blind.
    BlindValidatedPrediction,
    /// A physically validated prediction is also supported in field use.
    FieldSupportedPrediction,
    /// A physically validated prediction transfers across named regimes.
    TransferablePrediction,
    /// A physically validated prediction was independently reproduced.
    IndependentlyReproducedPrediction,
}

impl PortfolioClaimClass {
    /// Stable machine-readable claim-class name.
    #[must_use]
    pub const fn slug(self) -> &'static str {
        match self {
            Self::NumericallyVerified => "numerically-verified",
            Self::CrossCodeConsistent => "cross-code-consistent",
            Self::ValidatedPrediction => "validated-prediction",
            Self::BlindValidatedPrediction => "blind-validated-prediction",
            Self::FieldSupportedPrediction => "field-supported-prediction",
            Self::TransferablePrediction => "transferable-prediction",
            Self::IndependentlyReproducedPrediction => "independently-reproduced-prediction",
        }
    }

    /// Exact axes required by this claim class.
    #[must_use]
    pub const fn required_axes(self) -> &'static [EvidenceAxis] {
        match self {
            Self::NumericallyVerified => &[EvidenceAxis::NumericalVerification],
            Self::CrossCodeConsistent => &[
                EvidenceAxis::NumericalVerification,
                EvidenceAxis::CrossCodeAgreement,
            ],
            Self::ValidatedPrediction => &[EvidenceAxis::ControlledExperimentalValidation],
            Self::BlindValidatedPrediction => &[
                EvidenceAxis::ControlledExperimentalValidation,
                EvidenceAxis::BlindPredictiveValidation,
            ],
            Self::FieldSupportedPrediction => &[
                EvidenceAxis::ControlledExperimentalValidation,
                EvidenceAxis::FieldMonitoring,
            ],
            Self::TransferablePrediction => &[
                EvidenceAxis::ControlledExperimentalValidation,
                EvidenceAxis::TransferabilityAcrossRegimes,
            ],
            Self::IndependentlyReproducedPrediction => &[
                EvidenceAxis::ControlledExperimentalValidation,
                EvidenceAxis::IndependentReproduction,
            ],
        }
    }

    const fn tag(self) -> u8 {
        match self {
            Self::NumericallyVerified => 1,
            Self::CrossCodeConsistent => 2,
            Self::ValidatedPrediction => 3,
            Self::BlindValidatedPrediction => 4,
            Self::FieldSupportedPrediction => 5,
            Self::TransferablePrediction => 6,
            Self::IndependentlyReproducedPrediction => 7,
        }
    }
}

/// One exact claim-scoped coordinate supplied to a portfolio.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct PortfolioObservation {
    axis: EvidenceAxis,
    qoi: String,
    regime: String,
    source: ContentHash,
    independence_group: ContentHash,
}

impl PortfolioObservation {
    /// Validate and bind one external-evidence coordinate.
    pub fn try_new(
        axis: EvidenceAxis,
        qoi: impl Into<String>,
        regime: impl Into<String>,
        source: ContentHash,
        independence_group: ContentHash,
    ) -> Result<Self, PortfolioRefusal> {
        let qoi = qoi.into();
        let regime = regime.into();
        validate_id("qoi", &qoi)?;
        validate_id("regime", &regime)?;
        if source.0 == [0; 32] {
            return Err(PortfolioRefusal::ZeroHash { field: "source" });
        }
        if independence_group.0 == [0; 32] {
            return Err(PortfolioRefusal::ZeroHash {
                field: "independence_group",
            });
        }
        Ok(Self {
            axis,
            qoi,
            regime,
            source,
            independence_group,
        })
    }

    /// Evidence coordinate.
    #[must_use]
    pub const fn axis(&self) -> EvidenceAxis {
        self.axis
    }

    /// Exact QoI identifier.
    #[must_use]
    pub fn qoi(&self) -> &str {
        &self.qoi
    }

    /// Exact regime identifier.
    #[must_use]
    pub fn regime(&self) -> &str {
        &self.regime
    }

    /// Exact source artifact identity.
    #[must_use]
    pub const fn source(&self) -> ContentHash {
        self.source
    }

    /// Declared independence group used by the reproduction guard.
    #[must_use]
    pub const fn independence_group(&self) -> ContentHash {
        self.independence_group
    }
}

/// Bounded, canonical collection of claim-scoped evidence coordinates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EvidencePortfolio {
    observations: Vec<PortfolioObservation>,
    identity: ContentHash,
}

impl EvidencePortfolio {
    /// Validate, canonicalize, and deduplicate exact replay rows.
    pub fn try_new(mut observations: Vec<PortfolioObservation>) -> Result<Self, PortfolioRefusal> {
        if observations.len() > MAX_PORTFOLIO_OBSERVATIONS {
            return Err(PortfolioRefusal::ResourceLimit {
                resource: "observations",
                limit: MAX_PORTFOLIO_OBSERVATIONS,
                observed: observations.len(),
            });
        }
        observations.sort();
        observations.dedup();
        let identity = hash_domain(
            PORTFOLIO_IDENTITY_DOMAIN,
            &encode_observations(&observations),
        );
        Ok(Self {
            observations,
            identity,
        })
    }

    /// Canonical, exact observations.
    #[must_use]
    pub fn observations(&self) -> &[PortfolioObservation] {
        &self.observations
    }

    /// Canonical portfolio identity.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Admit a claim only when every exact axis requirement is present for
    /// the same QoI and regime.
    pub fn admit(
        &self,
        claim_class: PortfolioClaimClass,
        qoi: &str,
        regime: &str,
    ) -> Result<PortfolioAdmission, PortfolioRefusal> {
        validate_id("qoi", qoi)?;
        validate_id("regime", regime)?;

        let mut support = Vec::with_capacity(claim_class.required_axes().len());
        for &axis in claim_class.required_axes() {
            let observation = self
                .observations
                .iter()
                .find(|observation| {
                    observation.axis == axis
                        && observation.qoi == qoi
                        && observation.regime == regime
                })
                .ok_or_else(|| PortfolioRefusal::MissingAxis {
                    claim_class,
                    axis,
                    qoi: qoi.to_string(),
                    regime: regime.to_string(),
                })?;
            support.push(observation.clone());
        }

        if claim_class == PortfolioClaimClass::IndependentlyReproducedPrediction {
            let mut experiment_candidates = self.observations.iter().filter(|observation| {
                observation.axis == EvidenceAxis::ControlledExperimentalValidation
                    && observation.qoi == qoi
                    && observation.regime == regime
            });
            let pair = experiment_candidates.find_map(|experiment| {
                self.observations
                    .iter()
                    .find(|observation| {
                        observation.axis == EvidenceAxis::IndependentReproduction
                            && observation.qoi == qoi
                            && observation.regime == regime
                            && observation.source != experiment.source
                            && observation.independence_group != experiment.independence_group
                    })
                    .map(|reproduction| (experiment, reproduction))
            });
            let experiment_group = support
                .iter()
                .find(|observation| {
                    observation.axis == EvidenceAxis::ControlledExperimentalValidation
                })
                .map(|observation| observation.independence_group)
                .ok_or_else(|| PortfolioRefusal::MissingAxis {
                    claim_class,
                    axis: EvidenceAxis::ControlledExperimentalValidation,
                    qoi: qoi.to_string(),
                    regime: regime.to_string(),
                })?;
            let (experiment, reproduction) =
                pair.ok_or_else(|| PortfolioRefusal::IndependenceNotEstablished {
                    qoi: qoi.to_string(),
                    regime: regime.to_string(),
                    experiment_group,
                })?;
            if let Some(slot) = support.iter_mut().find(|observation| {
                observation.axis == EvidenceAxis::ControlledExperimentalValidation
            }) {
                *slot = experiment.clone();
            }
            if let Some(slot) = support
                .iter_mut()
                .find(|observation| observation.axis == EvidenceAxis::IndependentReproduction)
            {
                *slot = reproduction.clone();
            }
        }

        let identity = admission_identity(claim_class, qoi, regime, &support);
        Ok(PortfolioAdmission {
            claim_class,
            qoi: qoi.to_string(),
            regime: regime.to_string(),
            support,
            identity,
        })
    }
}

/// Opaque proof that the structural per-axis admission rule ran.
///
/// This is not source authentication and does not mint an evidence color.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PortfolioAdmission {
    claim_class: PortfolioClaimClass,
    qoi: String,
    regime: String,
    support: Vec<PortfolioObservation>,
    identity: ContentHash,
}

impl PortfolioAdmission {
    /// Admitted claim class.
    #[must_use]
    pub const fn claim_class(&self) -> PortfolioClaimClass {
        self.claim_class
    }

    /// Exact QoI.
    #[must_use]
    pub fn qoi(&self) -> &str {
        &self.qoi
    }

    /// Exact regime.
    #[must_use]
    pub fn regime(&self) -> &str {
        &self.regime
    }

    /// One selected observation for every required axis.
    #[must_use]
    pub fn support(&self) -> &[PortfolioObservation] {
        &self.support
    }

    /// Canonical admission identity.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Deterministic bounded log row for reviewers.
    #[must_use]
    pub fn render_log(&self) -> String {
        let axes = self
            .support
            .iter()
            .map(|observation| observation.axis.slug())
            .collect::<Vec<_>>()
            .join(",");
        format!(
            "schema={} claim_class={} qoi={} regime={} axes={} admission={}",
            EVIDENCE_PORTFOLIO_SCHEMA_VERSION,
            self.claim_class.slug(),
            self.qoi,
            self.regime,
            axes,
            hex(self.identity)
        )
    }
}

/// Typed portfolio-construction or admission refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PortfolioRefusal {
    /// A QoI or regime identifier is empty, too large, or control-bearing.
    InvalidIdentifier {
        /// Stable field name.
        field: &'static str,
        /// Stable refusal reason.
        reason: &'static str,
    },
    /// A content identity used the all-zero absence sentinel.
    ZeroHash {
        /// Stable field name.
        field: &'static str,
    },
    /// A bounded collection exceeded its hard cap.
    ResourceLimit {
        /// Capped resource.
        resource: &'static str,
        /// Maximum admitted value.
        limit: usize,
        /// Supplied value.
        observed: usize,
    },
    /// The named claim is missing one exact, non-substitutable axis.
    MissingAxis {
        /// Claim whose rule was evaluated.
        claim_class: PortfolioClaimClass,
        /// Required axis that had no same-QoI/same-regime observation.
        axis: EvidenceAxis,
        /// Exact QoI.
        qoi: String,
        /// Exact regime.
        regime: String,
    },
    /// No reproduction had both a distinct source and independence group.
    IndependenceNotEstablished {
        /// Exact QoI.
        qoi: String,
        /// Exact regime.
        regime: String,
        /// Reference experiment group retained for diagnosis.
        experiment_group: ContentHash,
    },
}

impl fmt::Display for PortfolioRefusal {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidIdentifier { field, reason } => {
                write!(formatter, "portfolio {field} identifier {reason}")
            }
            Self::ZeroHash { field } => {
                write!(formatter, "portfolio {field} cannot use the zero hash")
            }
            Self::ResourceLimit {
                resource,
                limit,
                observed,
            } => write!(
                formatter,
                "portfolio {resource} count {observed} exceeds limit {limit}"
            ),
            Self::MissingAxis {
                claim_class,
                axis,
                qoi,
                regime,
            } => write!(
                formatter,
                "portfolio claim {} for qoi={qoi} regime={regime} is missing axis {}",
                claim_class.slug(),
                axis.slug()
            ),
            Self::IndependenceNotEstablished { qoi, regime, .. } => write!(
                formatter,
                "portfolio qoi={qoi} regime={regime} has no reproduction with both a distinct source and independence group"
            ),
        }
    }
}

impl std::error::Error for PortfolioRefusal {}

fn validate_id(field: &'static str, value: &str) -> Result<(), PortfolioRefusal> {
    if value.trim().is_empty() {
        return Err(PortfolioRefusal::InvalidIdentifier {
            field,
            reason: "is blank",
        });
    }
    if value.len() > MAX_PORTFOLIO_ID_BYTES {
        return Err(PortfolioRefusal::InvalidIdentifier {
            field,
            reason: "exceeds the byte limit",
        });
    }
    if value.chars().any(char::is_control) {
        return Err(PortfolioRefusal::InvalidIdentifier {
            field,
            reason: "contains a control character",
        });
    }
    Ok(())
}

fn encode_observations(observations: &[PortfolioObservation]) -> Vec<u8> {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&EVIDENCE_PORTFOLIO_SCHEMA_VERSION.to_le_bytes());
    push_len(&mut bytes, observations.len());
    for observation in observations {
        bytes.push(observation.axis.tag());
        push_text(&mut bytes, &observation.qoi);
        push_text(&mut bytes, &observation.regime);
        bytes.extend_from_slice(&observation.source.0);
        bytes.extend_from_slice(&observation.independence_group.0);
    }
    bytes
}

fn admission_identity(
    claim_class: PortfolioClaimClass,
    qoi: &str,
    regime: &str,
    support: &[PortfolioObservation],
) -> ContentHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&EVIDENCE_PORTFOLIO_SCHEMA_VERSION.to_le_bytes());
    bytes.push(claim_class.tag());
    push_text(&mut bytes, qoi);
    push_text(&mut bytes, regime);
    bytes.extend_from_slice(&encode_observations(support));
    hash_domain(ADMISSION_IDENTITY_DOMAIN, &bytes)
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
    let mut out = String::with_capacity(64);
    for byte in hash.0 {
        out.push(char::from(DIGITS[usize::from(byte >> 4)]));
        out.push(char::from(DIGITS[usize::from(byte & 0x0f)]));
    }
    out
}
