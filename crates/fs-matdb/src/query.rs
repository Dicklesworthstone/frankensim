//! The query path (bead 5hmy, PR-4 of 5): every material answer is
//! `Evidence<PropertySample>` PLUS a [`PropertyUsageReceipt`] — never a
//! bare number.
//!
//! Discipline, in order:
//! - the query point is validated (finite, named axes);
//! - only claims whose [`fs_evidence::ValidityDomain`] CONTAINS the
//!   point are candidates — evaluation outside validity is a typed
//!   refusal, never a silent extrapolation;
//! - selection among candidates is an EXPLICIT policy; conflicting
//!   claims are never averaged into an invented canonical value —
//!   ambiguity refuses and names the candidates;
//! - the evidence slices map honestly: the datum's band is the STATED
//!   uncertainty (statistical slice); `Unstated` uncertainty maps to an
//!   explicit numerical no-claim, not a manufactured certificate; the
//!   claim's validity and in-domain fact live in the model slice; and
//!   the receipt records what was considered, selected, and decided.

use std::collections::BTreeMap;

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::{
    Evidence, ModelEvidence, NumericalCertificate, ProvenanceHash, SensitivitySummary,
    StatisticalCertificate,
};
use fs_qty::Dims;

use crate::{ClaimId, ClaimSet, InterpolationPolicy, MatDbError, PropertyValue, UncertaintyModel};

/// Hash domain for property-usage-receipt canonical identity.
const RECEIPT_HASH_DOMAIN: &str = "org.frankensim.fs-matdb.property-usage-receipt.v1";

/// The query evaluator's semantic version (recorded in every receipt;
/// bumped when selection or evaluation semantics change).
pub const MATDB_EVALUATOR_VERSION: u32 = 1;

/// A named-axis query point ("T" → 293.15, "normal_pressure" → 2.0e5).
/// Axis names match [`fs_evidence::ValidityDomain`] axis names.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct QueryPoint {
    axes: BTreeMap<String, f64>,
}

impl QueryPoint {
    /// An empty point (only unconstrained claims can answer it).
    #[must_use]
    pub fn new() -> QueryPoint {
        QueryPoint::default()
    }

    /// Set one named axis.
    ///
    /// # Errors
    /// [`MatDbError::NonFiniteQueryPoint`] for a non-finite coordinate.
    pub fn with(mut self, axis: impl Into<String>, value: f64) -> Result<QueryPoint, MatDbError> {
        let axis = axis.into();
        if !value.is_finite() {
            return Err(MatDbError::NonFiniteQueryPoint {
                axis,
                bits: value.to_bits(),
            });
        }
        self.axes.insert(axis, value);
        Ok(self)
    }

    /// The named coordinates.
    #[must_use]
    pub fn axes(&self) -> &BTreeMap<String, f64> {
        &self.axes
    }
}

/// How the answer chooses among in-domain candidate claims. Fusion is
/// explicit: no policy invents a canonical value from disagreeing
/// claims — ambiguity is a typed refusal naming the candidates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SelectionPolicy {
    /// Exactly one in-domain claim may exist; two or more refuse.
    SingleClaimOnly,
    /// Observation-backed claims outrank citation-only claims; the
    /// surviving set must still be a singleton.
    PreferObservationBacked,
}

impl SelectionPolicy {
    /// Stable receipt tag.
    #[must_use]
    pub fn tag(self) -> &'static str {
        match self {
            SelectionPolicy::SingleClaimOnly => "single-claim-only",
            SelectionPolicy::PreferObservationBacked => "prefer-observation-backed",
        }
    }
}

/// How the selected claim was evaluated at the point (successful paths
/// only — extrapolation never succeeds, it refuses).
#[derive(Debug, Clone, PartialEq)]
pub enum EvaluationDecision {
    /// A scalar claim valid across its whole validity box.
    ConstantWithinValidity,
    /// A single tabulated scalar (no abscissa involved).
    ExactScalar,
    /// A tabulated value hit exactly (bit-equal abscissa).
    ExactTabulated {
        /// The matched abscissa.
        at: f64,
    },
    /// Piecewise-linear interpolation strictly inside the knot span.
    LinearInside {
        /// Left bracketing knot abscissa.
        x_lo: f64,
        /// Right bracketing knot abscissa.
        x_hi: f64,
    },
}

impl EvaluationDecision {
    fn receipt_bytes(&self) -> Vec<u8> {
        match self {
            EvaluationDecision::ConstantWithinValidity => b"constant-within-validity".to_vec(),
            EvaluationDecision::ExactScalar => b"exact-scalar".to_vec(),
            EvaluationDecision::ExactTabulated { at } => {
                let mut out = b"exact-tabulated".to_vec();
                out.extend_from_slice(&at.to_bits().to_le_bytes());
                out
            }
            EvaluationDecision::LinearInside { x_lo, x_hi } => {
                let mut out = b"linear-inside".to_vec();
                out.extend_from_slice(&x_lo.to_bits().to_le_bytes());
                out.extend_from_slice(&x_hi.to_bits().to_le_bytes());
                out
            }
        }
    }
}

/// The evaluated sample a query returns (inside `Evidence<_>`).
#[derive(Debug, Clone, PartialEq)]
pub struct PropertySample {
    /// The evaluated SI value.
    pub value: f64,
    /// The value's dimensions.
    pub dims: Dims,
    /// The stated uncertainty model it inherits from its claim.
    pub uncertainty: UncertaintyModel,
}

/// The receipt every answer carries: what was asked, what was
/// considered, what was selected under which policy, how it was
/// evaluated, and which sources are load-bearing.
#[derive(Debug, Clone, PartialEq)]
pub struct PropertyUsageReceipt {
    /// The property name queried.
    pub property: String,
    /// The query point (named axes, canonical order).
    pub query_point: Vec<(String, f64)>,
    /// Every claim whose key matched the name, in insertion order —
    /// including out-of-domain claims (the receipt shows what was NOT
    /// eligible, so a narrow answer cannot masquerade as consensus).
    pub considered: Vec<ClaimId>,
    /// The in-domain candidates after validity filtering.
    pub in_domain: Vec<ClaimId>,
    /// The selected claim.
    pub selected: ClaimId,
    /// The selection policy's stable tag.
    pub policy: &'static str,
    /// How the value was produced.
    pub decision: EvaluationDecision,
    /// Whether the selected claim is observation-backed (specimen and
    /// process context exist). Citation-only answers can never be
    /// Validated-class downstream.
    pub observation_backed: bool,
    /// The evaluator's semantic version.
    pub evaluator_version: u32,
    /// Content hashes of the selected claim and its observations.
    pub source_hashes: Vec<ContentHash>,
}

impl PropertyUsageReceipt {
    /// Canonical receipt identity over every field.
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        let mut payload = Vec::new();
        let mut push = |part: &[u8]| {
            payload.extend_from_slice(&u64::try_from(part.len()).unwrap_or(u64::MAX).to_le_bytes());
            payload.extend_from_slice(part);
        };
        push(self.property.as_bytes());
        for (axis, value) in &self.query_point {
            push(axis.as_bytes());
            push(&value.to_bits().to_le_bytes());
        }
        for id in &self.considered {
            push(&id.0.0);
        }
        for id in &self.in_domain {
            push(&id.0.0);
        }
        push(&self.selected.0.0);
        push(self.policy.as_bytes());
        push(&self.decision.receipt_bytes());
        push(&[u8::from(self.observation_backed)]);
        push(&self.evaluator_version.to_le_bytes());
        for hash in &self.source_hashes {
            push(&hash.0);
        }
        hash_domain(RECEIPT_HASH_DOMAIN, &payload)
    }
}

/// A complete material answer: the evidence-carried sample plus its
/// usage receipt.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialAnswer {
    /// The sample with its honest evidence slices.
    pub evidence: Evidence<PropertySample>,
    /// The usage receipt.
    pub receipt: PropertyUsageReceipt,
}

impl ClaimSet {
    /// Answer a property query at a point under an explicit selection
    /// policy.
    ///
    /// # Errors
    /// [`MatDbError::UnknownProperty`] when no claim carries the name;
    /// [`MatDbError::NoClaimInDomain`] when the point is outside every
    /// claim's validity (THE extrapolation refusal);
    /// [`MatDbError::AmbiguousSelection`] when the policy cannot narrow
    /// the candidates to one claim (fusion must be explicit);
    /// [`MatDbError::MissingQueryAxis`] / [`MatDbError::OutsideKnotSpan`]
    /// / [`MatDbError::UnsupportedEvaluation`] from curve evaluation.
    pub fn query(
        &self,
        property: &str,
        point: &QueryPoint,
        policy: SelectionPolicy,
    ) -> Result<MaterialAnswer, MatDbError> {
        let considered_pairs = self.claims_for(property);
        if considered_pairs.is_empty() {
            return Err(MatDbError::UnknownProperty {
                property: property.to_string(),
            });
        }
        let considered: Vec<ClaimId> = considered_pairs.iter().map(|(id, _)| *id).collect();
        let in_domain_pairs: Vec<_> = considered_pairs
            .iter()
            .filter(|(_, claim)| claim.validity.contains(point.axes()))
            .collect();
        if in_domain_pairs.is_empty() {
            return Err(MatDbError::NoClaimInDomain {
                property: property.to_string(),
                considered: considered.len(),
            });
        }
        let in_domain: Vec<ClaimId> = in_domain_pairs.iter().map(|(id, _)| *id).collect();
        let selected_pairs: Vec<_> = match policy {
            SelectionPolicy::SingleClaimOnly => in_domain_pairs.clone(),
            SelectionPolicy::PreferObservationBacked => {
                let backed: Vec<_> = in_domain_pairs
                    .iter()
                    .filter(|(_, claim)| !claim.observations.is_empty())
                    .copied()
                    .collect();
                if backed.is_empty() {
                    in_domain_pairs.clone()
                } else {
                    backed
                }
            }
        };
        if selected_pairs.len() != 1 {
            return Err(MatDbError::AmbiguousSelection {
                property: property.to_string(),
                candidates: selected_pairs.iter().map(|(id, _)| *id).collect(),
            });
        }
        let (selected_id, claim) = selected_pairs[0];
        let (value, decision) = evaluate(&claim.value, claim.interpolation, point)?;

        let mut source_hashes = vec![selected_id.0];
        for observation in &claim.observations {
            source_hashes.push(observation.0);
        }
        let receipt = PropertyUsageReceipt {
            property: property.to_string(),
            query_point: point
                .axes()
                .iter()
                .map(|(axis, &v)| (axis.clone(), v))
                .collect(),
            considered,
            in_domain,
            selected: *selected_id,
            policy: policy.tag(),
            decision,
            observation_backed: !claim.observations.is_empty(),
            evaluator_version: MATDB_EVALUATOR_VERSION,
            source_hashes,
        };

        // Honest slice mapping. Numerical: the stated band as an
        // ESTIMATE around the value (never Exact/Enclosure — a datum is
        // not interval-certified numerics), or an explicit no-claim for
        // Unstated uncertainty. Statistical: the stated half-width.
        // Model: the claim's validity with the verified in-domain fact.
        let (numerical, statistical) = match claim.uncertainty {
            UncertaintyModel::Unstated => (
                NumericalCertificate::no_claim(),
                StatisticalCertificate::None,
            ),
            UncertaintyModel::HalfWidth {
                half_width,
                confidence,
            } => (
                NumericalCertificate::estimate(value - half_width, value + half_width),
                StatisticalCertificate::HalfWidth {
                    half_width,
                    confidence,
                },
            ),
            UncertaintyModel::RelativeHalfWidth {
                fraction,
                confidence,
            } => {
                let half_width = fraction * value.abs();
                (
                    NumericalCertificate::estimate(value - half_width, value + half_width),
                    StatisticalCertificate::HalfWidth {
                        half_width,
                        confidence,
                    },
                )
            }
        };
        let model = ModelEvidence {
            cards: vec![format!("fs-matdb:{property}")],
            assumptions: vec![format!(
                "claim provenance: {} ({})",
                claim.provenance.source, claim.provenance.license
            )],
            validity: claim.validity.clone(),
            discrepancy_rel: 0.0,
            in_domain: true,
        };
        let provenance = ProvenanceHash(u64::from_le_bytes(
            receipt.content_hash().0[..8]
                .try_into()
                .expect("hash has at least 8 bytes"),
        ));
        let evidence = Evidence {
            value: PropertySample {
                value,
                dims: claim.value.dims(),
                uncertainty: claim.uncertainty.clone(),
            },
            qoi: value,
            numerical,
            statistical,
            model,
            sensitivity: SensitivitySummary::default(),
            provenance,
            adjoint_ref: None,
        };
        Ok(MaterialAnswer { evidence, receipt })
    }
}

/// Evaluate a claim payload at a point under its interpolation policy.
fn evaluate(
    value: &PropertyValue,
    policy: InterpolationPolicy,
    point: &QueryPoint,
) -> Result<(f64, EvaluationDecision), MatDbError> {
    match (value, policy) {
        (PropertyValue::Scalar { value, .. }, InterpolationPolicy::ConstantWithinValidity) => {
            Ok((*value, EvaluationDecision::ConstantWithinValidity))
        }
        (PropertyValue::Scalar { value, .. }, InterpolationPolicy::TabulatedOnly) => {
            Ok((*value, EvaluationDecision::ExactScalar))
        }
        (PropertyValue::Scalar { .. }, InterpolationPolicy::LinearInside) => {
            Err(MatDbError::UnsupportedEvaluation {
                reason: "a scalar claim has no knot span to interpolate inside",
            })
        }
        (
            PropertyValue::Curve {
                abscissa, knots, ..
            },
            InterpolationPolicy::LinearInside,
        ) => {
            let x = *point
                .axes()
                .get(abscissa)
                .ok_or_else(|| MatDbError::MissingQueryAxis {
                    axis: abscissa.clone(),
                })?;
            let first = knots[0].0;
            let last = knots[knots.len() - 1].0;
            if x < first || x > last {
                return Err(MatDbError::OutsideKnotSpan {
                    axis: abscissa.clone(),
                    requested: x,
                    lo: first,
                    hi: last,
                });
            }
            if let Some(&(kx, ky)) = knots.iter().find(|&&(kx, _)| kx.to_bits() == x.to_bits()) {
                return Ok((ky, EvaluationDecision::ExactTabulated { at: kx }));
            }
            let window = knots
                .windows(2)
                .find(|w| w[0].0 <= x && x <= w[1].0)
                .expect("span containment guarantees a bracketing window");
            let (x0, y0) = window[0];
            let (x1, y1) = window[1];
            let t = (x - x0) / (x1 - x0);
            let y = y0 + t * (y1 - y0);
            Ok((y, EvaluationDecision::LinearInside { x_lo: x0, x_hi: x1 }))
        }
        (
            PropertyValue::Curve {
                abscissa, knots, ..
            },
            InterpolationPolicy::TabulatedOnly,
        ) => {
            let x = *point
                .axes()
                .get(abscissa)
                .ok_or_else(|| MatDbError::MissingQueryAxis {
                    axis: abscissa.clone(),
                })?;
            knots
                .iter()
                .find(|&&(kx, _)| kx.to_bits() == x.to_bits())
                .map(|&(kx, ky)| (ky, EvaluationDecision::ExactTabulated { at: kx }))
                .ok_or(MatDbError::UnsupportedEvaluation {
                    reason: "tabulated-only claim has no knot at the requested abscissa",
                })
        }
        (PropertyValue::Curve { .. }, InterpolationPolicy::ConstantWithinValidity) => {
            Err(MatDbError::UnsupportedEvaluation {
                reason: "a curve claim cannot be evaluated as a validity-wide constant",
            })
        }
    }
}
