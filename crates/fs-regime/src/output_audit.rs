//! Deterministic final-envelope audits for user-facing QoI claims.
//!
//! A validated model card is regional evidence. This module checks every
//! operating point against every card consumed by a QoI, partitions partial
//! sweeps exactly, and emits a typed receipt. Any exit forces the effective
//! color to an unbounded estimate; an explicit override is retained only as
//! an acknowledgement and can never restore claim strength.

use crate::cards::axis_distance_to_validity;
use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::uncertainty::{
    EngineeringUncertaintyBudget, EngineeringUncertaintyKind, EngineeringUncertaintyTerm,
    TermValue, UncertaintyArtifactRef, UncertaintyError,
};
use fs_evidence::{Color, ModelCard, ProvenanceHash, validate_color_payload};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

const OUTPUT_CLAIM_RECEIPT_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-regime.output-claim-receipt.v1";

/// One named operating point in dimensionless-group space.
#[derive(Debug, Clone, PartialEq)]
pub struct OperatingPoint {
    /// Stable point identity used in receipts and exact sweep partitions.
    pub id: String,
    /// Dimensionless group values at this point.
    pub groups: BTreeMap<String, f64>,
}

/// Minimal card projection required by the product-boundary regime audit.
///
/// This snapshot deliberately carries only identity and validity. It lets a
/// higher layer audit cards owned by other domains without inventing an
/// [`fs_evidence::Ambition`], discrepancy band, calibration artifact, or
/// failure catalog merely to satisfy the richer [`ModelCard`] schema.
#[derive(Debug, Clone, PartialEq)]
pub struct RegimeAuditCard {
    /// Stable card identity referenced by [`QoiClaim::model_cards`].
    pub name: String,
    /// Semantic version or other owner-defined version identity.
    pub version: String,
    /// Exact validity box supplied by the owning card type.
    pub validity: fs_evidence::ValidityDomain,
}

impl RegimeAuditCard {
    /// Project one owner-defined card identity and validity box into the audit.
    ///
    /// Structural identity validation remains centralized in
    /// [`audit_product_output_with_cards`], so malformed inputs return the
    /// same typed errors as evidence-model cards.
    #[must_use]
    pub fn new(
        name: impl Into<String>,
        version: impl Into<String>,
        validity: fs_evidence::ValidityDomain,
    ) -> Self {
        Self {
            name: name.into(),
            version: version.into(),
            validity,
        }
    }
}

impl From<&ModelCard> for RegimeAuditCard {
    fn from(card: &ModelCard) -> Self {
        Self::new(&card.name, &card.version, card.validity.clone())
    }
}

/// One user-facing quantity and the model cards consumed to produce it.
#[derive(Debug, Clone, PartialEq)]
pub struct QoiClaim {
    /// Stable quantity identity.
    pub qoi: String,
    /// Color before the final operating-envelope audit.
    pub color: Color,
    /// Names of every consumed model card.
    pub model_cards: Vec<String>,
    /// Explicit authorization to proceed despite a demotion, if supplied.
    pub override_acknowledgement: Option<OverrideAcknowledgement>,
}

/// An explicit decision to proceed with a demoted output.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OverrideAcknowledgement {
    /// Identity of the actor accepting the limitation.
    pub actor: String,
    /// Human-readable reason retained in the receipt.
    pub reason: String,
}

/// Exact coverage class for the audited operating points.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum EnvelopeCoverage {
    /// Every point satisfies every consumed model card.
    FullyInDomain,
    /// At least one point is in-domain and at least one is out-of-domain.
    Partial,
    /// No point satisfies every consumed model card.
    FullyOutOfDomain,
}

impl EnvelopeCoverage {
    fn name(self) -> &'static str {
        match self {
            Self::FullyInDomain => "fully-in-domain",
            Self::Partial => "partial",
            Self::FullyOutOfDomain => "fully-out-of-domain",
        }
    }
}

/// Why one card axis rejected one operating point.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum AxisViolationKind {
    /// The operating point omitted a constrained axis.
    Missing,
    /// The operating point supplied NaN or infinity.
    NonFinite,
    /// The value is below the inclusive lower bound.
    Below,
    /// The value is above the inclusive upper bound.
    Above,
    /// The model card itself contains unusable bounds.
    InvalidBounds,
}

impl AxisViolationKind {
    fn name(self) -> &'static str {
        match self {
            Self::Missing => "missing",
            Self::NonFinite => "non-finite",
            Self::Below => "below",
            Self::Above => "above",
            Self::InvalidBounds => "invalid-bounds",
        }
    }
}

/// One named, distance-scored model-card violation.
#[derive(Debug, Clone, PartialEq)]
pub struct RegimeViolation {
    /// Operating point that violated the card.
    pub point: String,
    /// Model card name.
    pub card: String,
    /// Model card semantic version.
    pub card_version: String,
    /// Constrained dimensionless axis.
    pub axis: String,
    /// Observed value, or `None` when the point omitted this axis.
    pub observed: Option<f64>,
    /// Inclusive lower validity bound.
    pub lo: f64,
    /// Inclusive upper validity bound.
    pub hi: f64,
    /// Stable reason class.
    pub kind: AxisViolationKind,
    /// Per-axis log-space distance from [`axis_distance_to_validity`].
    pub distance: f64,
}

/// Stable identity of one model card consumed by a QoI.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct ConsumedModelCard {
    /// Model card name.
    pub name: String,
    /// Model card semantic version.
    pub version: String,
}

/// Final audit receipt for one QoI.
#[derive(Debug, Clone, PartialEq)]
pub struct OutputClaimReceipt {
    /// Stable quantity identity.
    pub qoi: String,
    /// Color entering the final audit.
    pub original_color: Color,
    /// Color allowed at the product boundary.
    pub effective_color: Color,
    /// Incoming color retained for the exact in-domain partition, when nonempty.
    pub in_domain_color: Option<Color>,
    /// Demoted color assigned to the exact out-of-domain partition, when nonempty.
    pub out_of_domain_color: Option<Color>,
    /// Exact envelope coverage class.
    pub coverage: EnvelopeCoverage,
    /// Points satisfying every consumed card, sorted by identity.
    pub in_domain_points: Vec<String>,
    /// Points violating at least one consumed card, sorted by identity.
    pub out_of_domain_points: Vec<String>,
    /// Consumed card identities, sorted by name and deduplicated by registry law.
    pub model_cards: Vec<ConsumedModelCard>,
    /// Named violations, canonically sorted.
    pub violations: Vec<RegimeViolation>,
    /// Optional acknowledgement; never affects `effective_color`.
    pub override_acknowledgement: Option<OverrideAcknowledgement>,
}

impl OutputClaimReceipt {
    /// Whether the final audit weakened the output to an unbounded estimate.
    #[must_use]
    pub fn demoted(&self) -> bool {
        !self.out_of_domain_points.is_empty()
    }

    /// Domain-separated identity of the exact canonical receipt.
    #[must_use]
    pub fn content_id(&self) -> ContentHash {
        hash_domain(
            OUTPUT_CLAIM_RECEIPT_IDENTITY_DOMAIN,
            self.to_canonical_json().as_bytes(),
        )
    }

    /// Canonical JSON projection for ledger/report handoff.
    #[must_use]
    pub fn to_canonical_json(&self) -> String {
        let mut out = String::with_capacity(512);
        let _ = write!(
            out,
            "{{\"qoi\":{},\"coverage\":{},\"original_color\":{},\"effective_color\":{},\"in_domain_color\":{},\"out_of_domain_color\":{},",
            json_string(&self.qoi),
            json_string(self.coverage.name()),
            color_json(&self.original_color),
            color_json(&self.effective_color),
            optional_color_json(self.in_domain_color.as_ref()),
            optional_color_json(self.out_of_domain_color.as_ref()),
        );
        out.push_str("\"model_cards\":[");
        for (index, card) in self.model_cards.iter().enumerate() {
            if index > 0 {
                out.push(',');
            }
            let _ = write!(
                out,
                "{{\"name\":{},\"version\":{}}}",
                json_string(&card.name),
                json_string(&card.version)
            );
        }
        out.push(']');
        out.push(',');
        push_string_array(&mut out, "in_domain_points", &self.in_domain_points);
        out.push(',');
        push_string_array(&mut out, "out_of_domain_points", &self.out_of_domain_points);
        out.push_str(",\"violations\":[");
        for (index, violation) in self.violations.iter().enumerate() {
            if index > 0 {
                out.push(',');
            }
            let _ = write!(
                out,
                "{{\"point\":{},\"card\":{},\"card_version\":{},\"axis\":{},\"observed\":{},\"bounds\":[{},{}],\"kind\":{},\"distance\":{}}}",
                json_string(&violation.point),
                json_string(&violation.card),
                json_string(&violation.card_version),
                json_string(&violation.axis),
                violation
                    .observed
                    .map_or_else(|| "null".to_string(), json_f64),
                json_f64(violation.lo),
                json_f64(violation.hi),
                json_string(violation.kind.name()),
                json_f64(violation.distance),
            );
        }
        out.push_str("],\"override_acknowledgement\":");
        if let Some(acknowledgement) = &self.override_acknowledgement {
            let _ = write!(
                out,
                "{{\"actor\":{},\"reason\":{}}}",
                json_string(&acknowledgement.actor),
                json_string(&acknowledgement.reason),
            );
        } else {
            out.push_str("null");
        }
        out.push('}');
        out
    }

    /// Deterministic report-ready no-claim entry.
    #[must_use]
    pub fn no_claim_markdown(&self) -> Option<String> {
        if !self.demoted() {
            return None;
        }
        let mut out = format!(
            "- `{}`: **estimated / no dispersion claim**; coverage `{}`; {} of {} operating points outside all-card intersection.",
            self.qoi,
            self.coverage.name(),
            self.out_of_domain_points.len(),
            self.in_domain_points.len() + self.out_of_domain_points.len(),
        );
        for violation in &self.violations {
            let observed = violation
                .observed
                .map_or_else(|| "missing".to_string(), |value| format!("{value:.6e}"));
            let _ = write!(
                out,
                " `{}` / `{}` / `{}`: {} vs [{:.6e}, {:.6e}], distance {:.6e}.",
                violation.point,
                violation.card,
                violation.axis,
                observed,
                violation.lo,
                violation.hi,
                violation.distance,
            );
        }
        if let Some(acknowledgement) = &self.override_acknowledgement {
            let _ = write!(
                out,
                " Override acknowledged by `{}`: {}; acknowledgement does not restore color.",
                acknowledgement.actor, acknowledgement.reason,
            );
        }
        Some(out)
    }
}

/// Deterministic collection of final QoI audit receipts.
#[derive(Debug, Clone, PartialEq)]
pub struct ProductOutputAudit {
    /// One receipt per QoI, sorted by QoI identity.
    pub receipts: Vec<OutputClaimReceipt>,
    /// Deterministic provenance of the canonical receipt collection.
    pub provenance: ProvenanceHash,
}

/// Structural input error for a final-output audit.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OutputAuditError {
    /// This audit's dimensionless numeric operating points cannot carry typed axes.
    UnsupportedTypedAxes {
        /// Model card requiring a typed operating-point adapter.
        card: String,
    },
    /// At least one operating point is required.
    NoOperatingPoints,
    /// An identity is empty or contains surrounding whitespace.
    InvalidIdentity {
        /// Input field.
        field: &'static str,
        /// Offending value.
        value: String,
    },
    /// A supposedly unique identity was repeated.
    DuplicateIdentity {
        /// Input field.
        field: &'static str,
        /// Repeated value.
        value: String,
    },
    /// A QoI did not declare any model cards.
    NoModelCards {
        /// Offending QoI.
        qoi: String,
    },
    /// A QoI references a card absent from the supplied registry.
    UnknownModelCard {
        /// Offending QoI.
        qoi: String,
        /// Missing card name.
        card: String,
    },
    /// A QoI entered the product boundary with a malformed color payload.
    InvalidColor {
        /// Offending QoI.
        qoi: String,
        /// Structural color diagnosis.
        reason: String,
    },
}

impl core::fmt::Display for OutputAuditError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NoOperatingPoints => {
                write!(f, "output audit requires at least one operating point")
            }
            Self::UnsupportedTypedAxes { card } => write!(
                f,
                "model card {card:?} requires typed operating points; this audit profile cannot represent them"
            ),
            Self::InvalidIdentity { field, value } => {
                write!(f, "invalid {field} identity {value:?}")
            }
            Self::DuplicateIdentity { field, value } => {
                write!(f, "duplicate {field} identity {value:?}")
            }
            Self::NoModelCards { qoi } => write!(f, "QoI {qoi:?} consumed no model cards"),
            Self::UnknownModelCard { qoi, card } => {
                write!(f, "QoI {qoi:?} references unknown model card {card:?}")
            }
            Self::InvalidColor { qoi, reason } => {
                write!(f, "QoI {qoi:?} has invalid color: {reason}")
            }
        }
    }
}

impl std::error::Error for OutputAuditError {}

/// Refusal to apply a final-envelope receipt to an engineering budget.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OutputAuditBudgetError {
    /// The receipt and budget describe different quantities.
    QoiMismatch {
        /// Quantity named by the receipt.
        receipt_qoi: String,
        /// Quantity named by the budget.
        budget_qoi: String,
    },
    /// The derived eight-term budget failed its normal admission rules.
    Uncertainty(UncertaintyError),
}

impl core::fmt::Display for OutputAuditBudgetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::QoiMismatch {
                receipt_qoi,
                budget_qoi,
            } => write!(
                f,
                "cannot apply regime receipt for {receipt_qoi:?} to budget for {budget_qoi:?}"
            ),
            Self::Uncertainty(error) => write!(f, "regime-derived budget refused: {error}"),
        }
    }
}

impl std::error::Error for OutputAuditBudgetError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::QoiMismatch { .. } => None,
            Self::Uncertainty(error) => Some(error),
        }
    }
}

impl From<UncertaintyError> for OutputAuditBudgetError {
    fn from(error: UncertaintyError) -> Self {
        Self::Uncertainty(error)
    }
}

/// Apply one final-envelope receipt to the matching eight-term budget.
///
/// A fully in-domain receipt returns an exact clone. A demoted receipt makes
/// the model-form term explicitly unknown and provenance-links that declaration
/// to the exact receipt. If model form participated in a covariance block, all
/// members of that block become unknown because retaining their joint finite
/// representation after removing one member would be unsound.
///
/// # Errors
/// Refuses a receipt/budget QoI mismatch or a derived term that exceeds the
/// engineering-budget admission bounds.
pub fn apply_output_audit_to_budget(
    receipt: &OutputClaimReceipt,
    budget: &EngineeringUncertaintyBudget,
) -> Result<EngineeringUncertaintyBudget, OutputAuditBudgetError> {
    if receipt.qoi != budget.qoi() {
        return Err(OutputAuditBudgetError::QoiMismatch {
            receipt_qoi: receipt.qoi.clone(),
            budget_qoi: budget.qoi().to_string(),
        });
    }
    if !receipt.demoted() {
        return Ok(budget.clone());
    }

    let audit_id = receipt.content_id();
    let audit_provenance = UncertaintyArtifactRef::new("regime-output-audit", audit_id)?;
    let model_term = budget.term(EngineeringUncertaintyKind::ModelForm);
    let coupled_kinds = match model_term.value() {
        TermValue::CorrelatedBlock(block) => Some((block.id(), block.members())),
        _ => None,
    };
    let model_reason = model_form_reason(receipt, model_term);
    let mut terms = Vec::with_capacity(EngineeringUncertaintyKind::ALL.len());
    for term in budget.terms() {
        let value = if term.kind() == EngineeringUncertaintyKind::ModelForm {
            TermValue::unknown(&model_reason)?
        } else if let Some((block_id, members)) = coupled_kinds
            && members.contains(&term.kind())
        {
            TermValue::unknown(format!(
                "fs-regime audit {audit_id} invalidated covariance block {block_id:?}: its model-form member exited the validated domain"
            ))?
        } else {
            terms.push(term.clone());
            continue;
        };
        terms.push(EngineeringUncertaintyTerm::try_new(
            term.kind(),
            value,
            audit_provenance.clone(),
        )?);
    }
    EngineeringUncertaintyBudget::try_new(budget.qoi(), budget.unit(), terms).map_err(Into::into)
}

fn model_form_reason(receipt: &OutputClaimReceipt, prior: &EngineeringUncertaintyTerm) -> String {
    let audit_id = receipt.content_id();
    let mut reason = format!(
        "fs-regime audit {audit_id}: qoi={:?}, coverage={}, violations=[",
        receipt.qoi,
        receipt.coverage.name(),
    );
    for (index, violation) in receipt.violations.iter().enumerate() {
        if index > 0 {
            reason.push_str(", ");
        }
        let observed = violation
            .observed
            .map_or_else(|| "missing".to_string(), |value| format!("{value:.17e}"));
        let _ = write!(
            reason,
            "point={:?}/card={:?}@{:?}/axis={:?}/kind={}/observed={}/bounds=[{:.17e},{:.17e}]/distance={:.17e}",
            violation.point,
            violation.card,
            violation.card_version,
            violation.axis,
            violation.kind.name(),
            observed,
            violation.lo,
            violation.hi,
            violation.distance,
        );
    }
    reason.push_str("]; prior-model-form=");
    match prior.value() {
        TermValue::IntervalBound { lower, upper } => {
            let _ = write!(reason, "interval[{lower:.17e},{upper:.17e}]");
        }
        TermValue::Distribution(_) => reason.push_str("distribution"),
        TermValue::Ensemble(_) => reason.push_str("ensemble"),
        TermValue::CorrelatedBlock(block) => {
            let _ = write!(reason, "covariance-block({:?})", block.id());
        }
        TermValue::Unknown {
            reason: prior_reason,
        } => {
            let _ = write!(reason, "unknown({prior_reason:?})");
        }
        TermValue::Negligible { justification } => {
            let _ = write!(reason, "negligible({justification:?})");
        }
        _ => reason.push_str("future-representation"),
    }
    let _ = write!(
        reason,
        "; prior-authority={:?}@{}",
        prior.provenance().role(),
        prior.provenance().digest(),
    );
    if receipt.override_acknowledgement.is_some() {
        reason.push_str("; override-acknowledged-without-color-restoration");
    }
    reason
}

/// Audit every QoI/card/operating-point combination at the product boundary.
///
/// # Errors
/// Refuses missing points, malformed or duplicate identities, cardless QoIs,
/// duplicate registry names, and references to absent model cards.
pub fn audit_product_output(
    registry: &[ModelCard],
    operating_points: &[OperatingPoint],
    claims: &[QoiClaim],
) -> Result<ProductOutputAudit, OutputAuditError> {
    let audit_cards = registry
        .iter()
        .map(RegimeAuditCard::from)
        .collect::<Vec<_>>();
    audit_product_output_with_cards(&audit_cards, operating_points, claims)
}

/// Audit every QoI against owner-neutral card identity/validity projections.
///
/// Use this entry point for card classes whose authoritative schema is not
/// [`ModelCard`]. The projection must preserve the owning card's exact identity,
/// version, and validity box; this API does not assert or synthesize any other
/// evidence-card fields.
///
/// # Errors
/// Refuses missing points, malformed or duplicate identities, cardless QoIs,
/// duplicate registry names, and references to absent model cards.
pub fn audit_product_output_with_cards(
    registry: &[RegimeAuditCard],
    operating_points: &[OperatingPoint],
    claims: &[QoiClaim],
) -> Result<ProductOutputAudit, OutputAuditError> {
    if operating_points.is_empty() {
        return Err(OutputAuditError::NoOperatingPoints);
    }
    let cards = unique_cards(registry)?;
    let points = unique_points(operating_points)?;
    let claims = unique_claims(claims)?;

    let mut receipts = Vec::with_capacity(claims.len());
    for claim in claims {
        let model_cards = canonical_card_names(claim)?;
        let selected_cards = model_cards
            .iter()
            .map(|name| {
                cards
                    .get(name)
                    .copied()
                    .ok_or_else(|| OutputAuditError::UnknownModelCard {
                        qoi: claim.qoi.clone(),
                        card: name.clone(),
                    })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let consumed_cards = selected_cards
            .iter()
            .map(|card| ConsumedModelCard {
                name: card.name.clone(),
                version: card.version.clone(),
            })
            .collect::<Vec<_>>();
        receipts.push(audit_claim(
            claim,
            &points,
            &consumed_cards,
            &selected_cards,
        ));
    }
    receipts.sort_by(|left, right| left.qoi.cmp(&right.qoi));
    let canonical = receipts
        .iter()
        .map(OutputClaimReceipt::to_canonical_json)
        .collect::<Vec<_>>()
        .join("\n");
    Ok(ProductOutputAudit {
        receipts,
        provenance: ProvenanceHash::of_bytes(canonical.as_bytes()),
    })
}

fn unique_cards(
    registry: &[RegimeAuditCard],
) -> Result<BTreeMap<String, &RegimeAuditCard>, OutputAuditError> {
    let mut cards = BTreeMap::new();
    for card in registry {
        if card.validity.has_typed_axes() {
            return Err(OutputAuditError::UnsupportedTypedAxes {
                card: card.name.clone(),
            });
        }
        validate_identity("model-card", &card.name)?;
        validate_identity("model-card-version", &card.version)?;
        for axis in card.validity.bounds().keys() {
            validate_identity("model-card-axis", axis)?;
        }
        if cards.insert(card.name.clone(), card).is_some() {
            return Err(OutputAuditError::DuplicateIdentity {
                field: "model-card",
                value: card.name.clone(),
            });
        }
    }
    Ok(cards)
}

fn unique_points(
    operating_points: &[OperatingPoint],
) -> Result<Vec<&OperatingPoint>, OutputAuditError> {
    let mut by_id = BTreeMap::new();
    for point in operating_points {
        validate_identity("operating-point", &point.id)?;
        for axis in point.groups.keys() {
            validate_identity("operating-point-axis", axis)?;
        }
        if by_id.insert(point.id.clone(), point).is_some() {
            return Err(OutputAuditError::DuplicateIdentity {
                field: "operating-point",
                value: point.id.clone(),
            });
        }
    }
    Ok(by_id.into_values().collect())
}

fn unique_claims(claims: &[QoiClaim]) -> Result<Vec<&QoiClaim>, OutputAuditError> {
    let mut by_qoi = BTreeMap::new();
    for claim in claims {
        validate_identity("qoi", &claim.qoi)?;
        validate_color_payload(&claim.color).map_err(|error| OutputAuditError::InvalidColor {
            qoi: claim.qoi.clone(),
            reason: error.to_string(),
        })?;
        if by_qoi.insert(claim.qoi.clone(), claim).is_some() {
            return Err(OutputAuditError::DuplicateIdentity {
                field: "qoi",
                value: claim.qoi.clone(),
            });
        }
        if let Some(acknowledgement) = &claim.override_acknowledgement {
            validate_identity("override-actor", &acknowledgement.actor)?;
            validate_identity("override-reason", &acknowledgement.reason)?;
        }
    }
    Ok(by_qoi.into_values().collect())
}

fn canonical_card_names(claim: &QoiClaim) -> Result<Vec<String>, OutputAuditError> {
    if claim.model_cards.is_empty() {
        return Err(OutputAuditError::NoModelCards {
            qoi: claim.qoi.clone(),
        });
    }
    let mut names = BTreeSet::new();
    for name in &claim.model_cards {
        validate_identity("model-card-reference", name)?;
        names.insert(name.clone());
    }
    Ok(names.into_iter().collect())
}

fn validate_identity(field: &'static str, value: &str) -> Result<(), OutputAuditError> {
    if value.is_empty() || value.trim() != value {
        Err(OutputAuditError::InvalidIdentity {
            field,
            value: value.to_string(),
        })
    } else {
        Ok(())
    }
}

fn audit_claim(
    claim: &QoiClaim,
    points: &[&OperatingPoint],
    model_cards: &[ConsumedModelCard],
    selected_cards: &[&RegimeAuditCard],
) -> OutputClaimReceipt {
    let mut in_domain_points = Vec::new();
    let mut out_of_domain_points = Vec::new();
    let mut violations = Vec::new();
    for point in points {
        let before = violations.len();
        for card in selected_cards {
            collect_violations(point, card, &mut violations);
        }
        if before == violations.len() {
            in_domain_points.push(point.id.clone());
        } else {
            out_of_domain_points.push(point.id.clone());
        }
    }
    let coverage = if out_of_domain_points.is_empty() {
        EnvelopeCoverage::FullyInDomain
    } else if in_domain_points.is_empty() {
        EnvelopeCoverage::FullyOutOfDomain
    } else {
        EnvelopeCoverage::Partial
    };
    violations.sort_by(|left, right| {
        left.point
            .cmp(&right.point)
            .then(left.card.cmp(&right.card))
            .then(left.axis.cmp(&right.axis))
    });
    let out_of_domain_color = if out_of_domain_points.is_empty() {
        None
    } else {
        Some(demoted_color(claim, model_cards, &violations))
    };
    let effective_color = out_of_domain_color
        .clone()
        .unwrap_or_else(|| claim.color.clone());
    let in_domain_color = (!in_domain_points.is_empty()).then(|| claim.color.clone());
    OutputClaimReceipt {
        qoi: claim.qoi.clone(),
        original_color: claim.color.clone(),
        effective_color,
        in_domain_color,
        out_of_domain_color,
        coverage,
        in_domain_points,
        out_of_domain_points,
        model_cards: model_cards.to_vec(),
        violations,
        override_acknowledgement: claim.override_acknowledgement.clone(),
    }
}

fn collect_violations(
    point: &OperatingPoint,
    card: &RegimeAuditCard,
    violations: &mut Vec<RegimeViolation>,
) {
    for (axis, &(lo, hi)) in card.validity.bounds() {
        let observed = point.groups.get(axis).copied();
        let kind = if !lo.is_finite() || !hi.is_finite() || lo > hi {
            Some(AxisViolationKind::InvalidBounds)
        } else {
            match observed {
                None => Some(AxisViolationKind::Missing),
                Some(value) if !value.is_finite() => Some(AxisViolationKind::NonFinite),
                Some(value) if value < lo => Some(AxisViolationKind::Below),
                Some(value) if value > hi => Some(AxisViolationKind::Above),
                Some(_) => None,
            }
        };
        if let Some(kind) = kind {
            violations.push(RegimeViolation {
                point: point.id.clone(),
                card: card.name.clone(),
                card_version: card.version.clone(),
                axis: axis.clone(),
                observed,
                lo,
                hi,
                kind,
                distance: axis_distance_to_validity(observed, lo, hi),
            });
        }
    }
}

fn demoted_color(
    claim: &QoiClaim,
    model_cards: &[ConsumedModelCard],
    violations: &[RegimeViolation],
) -> Color {
    let mut identity = Vec::new();
    push_identity_field(&mut identity, b"regime-output-audit-v1");
    push_identity_field(&mut identity, claim.qoi.as_bytes());
    push_identity_field(&mut identity, &claim.color.canonical_bytes());
    for card in model_cards {
        push_identity_field(&mut identity, card.name.as_bytes());
        push_identity_field(&mut identity, card.version.as_bytes());
    }
    for violation in violations {
        push_identity_field(&mut identity, violation.point.as_bytes());
        push_identity_field(&mut identity, violation.card.as_bytes());
        push_identity_field(&mut identity, violation.card_version.as_bytes());
        push_identity_field(&mut identity, violation.axis.as_bytes());
        match violation.observed {
            Some(value) => {
                push_identity_field(&mut identity, b"present");
                push_identity_field(&mut identity, &value.to_bits().to_le_bytes());
            }
            None => push_identity_field(&mut identity, b"missing"),
        }
        push_identity_field(&mut identity, &violation.lo.to_bits().to_le_bytes());
        push_identity_field(&mut identity, &violation.hi.to_bits().to_le_bytes());
        push_identity_field(&mut identity, violation.kind.name().as_bytes());
        push_identity_field(&mut identity, &violation.distance.to_bits().to_le_bytes());
    }
    Color::Estimated {
        estimator: format!(
            "regime-audit:{:016x}",
            ProvenanceHash::of_bytes(&identity).0
        ),
        dispersion: f64::INFINITY,
    }
}

fn push_identity_field(identity: &mut Vec<u8>, bytes: &[u8]) {
    identity.extend_from_slice(bytes.len().to_string().as_bytes());
    identity.push(b':');
    identity.extend_from_slice(bytes);
}

fn push_string_array(out: &mut String, key: &str, values: &[String]) {
    let _ = write!(out, "{}:[", json_string(key));
    for (index, value) in values.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        out.push_str(&json_string(value));
    }
    out.push(']');
}

fn json_f64(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.17e}")
    } else {
        json_string(&format!("non-finite:{value}"))
    }
}

fn color_payload_json(color: &Color) -> String {
    if matches!(color, Color::Validated { regime, .. } if regime.has_typed_axes()) {
        return color.payload_json();
    }
    match color {
        Color::Verified { lo, hi } => {
            format!("{{\"interval\":[{},{}]}}", json_f64(*lo), json_f64(*hi))
        }
        Color::Validated { regime, dataset } => {
            let mut axes = String::new();
            for (index, (axis, (lo, hi))) in regime.bounds().iter().enumerate() {
                if index > 0 {
                    axes.push(',');
                }
                let _ = write!(
                    axes,
                    "{}:[{},{}]",
                    json_string(axis),
                    json_f64(*lo),
                    json_f64(*hi)
                );
            }
            format!(
                "{{\"dataset\":{},\"regime\":{{{axes}}}}}",
                json_string(dataset)
            )
        }
        Color::Estimated {
            estimator,
            dispersion,
        } => format!(
            "{{\"estimator\":{},\"dispersion\":{}}}",
            json_string(estimator),
            json_f64(*dispersion)
        ),
    }
}

fn color_json(color: &Color) -> String {
    format!(
        "{{\"name\":{},\"payload\":{}}}",
        json_string(color.name()),
        color_payload_json(color)
    )
}

fn optional_color_json(color: Option<&Color>) -> String {
    color.map_or_else(|| "null".to_string(), color_json)
}

fn json_string(value: &str) -> String {
    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            character if u32::from(character) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", u32::from(character));
            }
            character => out.push(character),
        }
    }
    out.push('"');
    out
}
