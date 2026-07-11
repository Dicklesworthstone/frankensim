//! Discrepancy models and model bracketing (patch Rev B mechanisms 2–3).
//!
//! A [`DiscrepancyModel`] is fit from paired two-fidelity evaluations (the
//! ledger's accumulating corpus) and answers "how wrong is the cheap model
//! HERE" — refusing with a teaching [`OutOfDomain`] when asked outside the
//! region it has data for or with a parameter key set that differs from the
//! exact training schema (the surrogate out-of-distribution guard). v1 is
//! deliberately statistics-light: an observed parameter box plus
//! mean/max relative discrepancy — honest bookkeeping, not learning
//! (learned discrepancy models arrive with FrankenTorch — CONTRACT
//! no-claims).
//!
//! A [`ModelBracket`] handles weakly-understood physics (the vessel
//! flagship's contact-line mitigation): run EVERY plausible model, report
//! the QoI spread as an enclosure plus a model-form band — sensitivity to
//! the modeling choice, not pretended certainty.

use crate::{
    Evidence, ModelEvidence, NumericalCertificate, ProvenanceHash, SensitivitySummary,
    StatisticalCertificate, ValidityDomain, color_identity_reason, color_leaf_identity_reason,
};
use core::fmt;
use std::collections::BTreeMap;

const MAX_TRAINING_PAIRS: usize = 65_536;
const MAX_TRAINING_PARAMETERS: usize = 1_024;
const MAX_TRAINING_COORDINATES: usize = 1_048_576;
const MIN_BRACKET_MEMBERS: usize = 2;
const MAX_BRACKET_MEMBERS: usize = 1_024;

/// One paired two-fidelity evaluation at a parameter point.
#[derive(Debug, Clone, PartialEq)]
pub struct FidelityPair {
    /// Where in parameter space the pair was evaluated.
    pub params: BTreeMap<String, f64>,
    /// Low-fidelity QoI.
    pub lo_fi: f64,
    /// High-fidelity QoI (the reference).
    pub hi_fi: f64,
}

/// The queried band at an in-domain point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DiscrepancyBand {
    /// Mean relative discrepancy across training pairs.
    pub mean_rel: f64,
    /// Worst observed relative discrepancy (the conservative band).
    pub max_rel: f64,
}

/// Structured discrepancy-model or model-bracket refusal.
///
/// The historical name is retained for API compatibility; bracket construction
/// uses the same error so every public discrepancy-evidence path fails with a
/// nameable, deterministic reason.
#[derive(Debug, Clone, PartialEq)]
pub enum FitError {
    /// No paired observations were supplied.
    EmptyTrainingSet,
    /// The paired corpus exceeds the bounded v1 bookkeeping budget.
    TooManyTrainingPairs {
        /// Supplied pair count.
        count: usize,
        /// Maximum admitted pair count.
        maximum: usize,
    },
    /// A pair did not declare any parameter coordinates.
    EmptyParameterSchema,
    /// The parameter schema exceeds the bounded validity-domain budget.
    TooManyTrainingParameters {
        /// Supplied parameter count.
        count: usize,
        /// Maximum admitted parameter count.
        maximum: usize,
    },
    /// Pair count times parameter count exceeds the bounded fit-work budget.
    TooManyTrainingCoordinates {
        /// Supplied pair count.
        pairs: usize,
        /// Parameters per pair.
        parameters: usize,
        /// Maximum admitted coordinate count.
        maximum: usize,
    },
    /// A pair's parameter names differ from the first pair's schema.
    InconsistentParameterSchema {
        /// Zero-based pair index.
        pair_index: usize,
    },
    /// A parameter name cannot become a validity-domain identity.
    InvalidParameterIdentity {
        /// Zero-based pair index.
        pair_index: usize,
        /// Shared identity-grammar rejection reason.
        reason: &'static str,
    },
    /// A training pair contains a NaN or infinite QoI.
    NonFiniteTrainingQoi {
        /// Zero-based pair index.
        pair_index: usize,
    },
    /// A training coordinate is NaN or infinite.
    NonFiniteTrainingParameter {
        /// Zero-based pair index.
        pair_index: usize,
        /// Bounded, already-validated parameter identity.
        param: String,
    },
    /// `evidence_at` was given an unusable model-card identity.
    InvalidCardIdentity {
        /// Shared leaf-identity rejection reason.
        reason: &'static str,
    },
    /// An otherwise valid query is outside the observed parameter box.
    QueryOutOfDomain(OutOfDomain),
    /// A bracket cannot measure model-form spread with fewer than two models.
    TooFewBracketMembers {
        /// Number of admitted unique members.
        count: usize,
        /// Required minimum.
        minimum: usize,
    },
    /// A bracket exceeded its bounded member budget.
    TooManyBracketMembers {
        /// Maximum admitted member count.
        maximum: usize,
    },
    /// A member name cannot become a model-card identity.
    InvalidBracketMemberIdentity {
        /// Shared leaf-identity rejection reason.
        reason: &'static str,
    },
    /// Two member QoIs claimed the same model identity.
    DuplicateBracketMember {
        /// The duplicate, bounded identity.
        name: String,
    },
    /// A member QoI was NaN or infinite.
    NonFiniteBracketQoi {
        /// The member's bounded identity.
        name: String,
    },
}

impl fmt::Display for FitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FitError::EmptyTrainingSet => {
                write!(f, "discrepancy fit refused: the training set is empty")
            }
            FitError::TooManyTrainingPairs { count, maximum } => write!(
                f,
                "discrepancy fit refused: {count} training pairs exceed the bounded maximum {maximum}"
            ),
            FitError::EmptyParameterSchema => write!(
                f,
                "discrepancy fit refused: training pairs must share a non-empty parameter schema"
            ),
            FitError::TooManyTrainingParameters { count, maximum } => write!(
                f,
                "discrepancy fit refused: {count} parameters exceed the bounded maximum {maximum}"
            ),
            FitError::TooManyTrainingCoordinates {
                pairs,
                parameters,
                maximum,
            } => write!(
                f,
                "discrepancy fit refused: {pairs} pairs x {parameters} parameters exceed the bounded {maximum}-coordinate work budget"
            ),
            FitError::InconsistentParameterSchema { pair_index } => write!(
                f,
                "discrepancy fit refused: training pair {pair_index} does not have exactly the first pair's parameter schema"
            ),
            FitError::InvalidParameterIdentity { pair_index, reason } => write!(
                f,
                "discrepancy fit refused: training pair {pair_index} has an invalid parameter identity ({reason})"
            ),
            FitError::NonFiniteTrainingQoi { pair_index } => write!(
                f,
                "discrepancy fit refused: training pair {pair_index} has a non-finite QoI"
            ),
            FitError::NonFiniteTrainingParameter { pair_index, param } => write!(
                f,
                "discrepancy fit refused: training pair {pair_index} parameter `{param}` is non-finite"
            ),
            FitError::InvalidCardIdentity { reason } => write!(
                f,
                "discrepancy evidence refused: the model-card identity is invalid ({reason})"
            ),
            FitError::QueryOutOfDomain(error) => error.fmt(f),
            FitError::TooFewBracketMembers { count, minimum } => write!(
                f,
                "model bracket refused: {count} unique model member(s) cannot measure model-choice spread; supply at least {minimum}"
            ),
            FitError::TooManyBracketMembers { maximum } => write!(
                f,
                "model bracket refused: member count exceeds the bounded maximum {maximum}"
            ),
            FitError::InvalidBracketMemberIdentity { reason } => write!(
                f,
                "model bracket refused: a member identity is invalid ({reason})"
            ),
            FitError::DuplicateBracketMember { name } => write!(
                f,
                "model bracket refused: duplicate model identity `{name}` is ambiguous"
            ),
            FitError::NonFiniteBracketQoi { name } => write!(
                f,
                "model bracket refused: member `{name}` has a non-finite QoI"
            ),
        }
    }
}

impl core::error::Error for FitError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            FitError::QueryOutOfDomain(error) => Some(error),
            _ => None,
        }
    }
}

/// The out-of-distribution refusal, naming the violated parameter (the
/// diagnosis an agent needs to decide between gathering data and
/// escalating fidelity).
#[derive(Debug, Clone, PartialEq)]
pub struct OutOfDomain {
    /// The parameter outside the trained box, missing from the query, or absent
    /// from the training schema.
    pub param: String,
    /// The queried value (`None` = the query omitted the parameter).
    pub value: Option<f64>,
    /// The trained box for that parameter. `None` means the query supplied an
    /// unexpected parameter that was never part of the training schema.
    pub trained: Option<(f64, f64)>,
}

impl fmt::Display for OutOfDomain {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match (self.value, self.trained) {
            (Some(v), Some(trained)) => write!(
                f,
                "discrepancy query out of domain: `{}` = {v} lies outside the trained box \
                 [{}, {}] — the band would be extrapolation; gather pairs there or escalate \
                 fidelity",
                self.param, trained.0, trained.1
            ),
            (None, Some(trained)) => write!(
                f,
                "discrepancy query out of domain: `{}` was not supplied but the model was \
                 trained on it (trained box [{}, {}])",
                self.param, trained.0, trained.1
            ),
            (Some(v), None) => write!(
                f,
                "discrepancy query out of domain: unexpected parameter `{}` = {v} was not \
                 present in the exact training schema; remove it or refit the discrepancy model \
                 with that dimension",
                self.param
            ),
            (None, None) => write!(
                f,
                "discrepancy query out of domain: unexpected parameter `{}` was not present in \
                 the exact training schema",
                self.param
            ),
        }
    }
}

impl core::error::Error for OutOfDomain {}

/// A two-fidelity discrepancy model: observed box + mean/max relative
/// discrepancy (v1 bookkeeping — see module docs).
#[derive(Debug, Clone, PartialEq)]
pub struct DiscrepancyModel {
    observed: ValidityDomain,
    trained_bounds: BTreeMap<String, (f64, f64)>,
    band: DiscrepancyBand,
    pairs: usize,
}

fn validate_training_shape(pair_count: usize, parameter_count: usize) -> Result<(), FitError> {
    if pair_count == 0 {
        return Err(FitError::EmptyTrainingSet);
    }
    if pair_count > MAX_TRAINING_PAIRS {
        return Err(FitError::TooManyTrainingPairs {
            count: pair_count,
            maximum: MAX_TRAINING_PAIRS,
        });
    }
    if parameter_count == 0 {
        return Err(FitError::EmptyParameterSchema);
    }
    if parameter_count > MAX_TRAINING_PARAMETERS {
        return Err(FitError::TooManyTrainingParameters {
            count: parameter_count,
            maximum: MAX_TRAINING_PARAMETERS,
        });
    }
    if pair_count
        .checked_mul(parameter_count)
        .is_none_or(|coordinates| coordinates > MAX_TRAINING_COORDINATES)
    {
        return Err(FitError::TooManyTrainingCoordinates {
            pairs: pair_count,
            parameters: parameter_count,
            maximum: MAX_TRAINING_COORDINATES,
        });
    }
    Ok(())
}

impl DiscrepancyModel {
    /// Fit from paired evaluations. The observed box is the per-parameter
    /// min/max over training points; the band is mean/max of
    /// `|hi - lo| / max(|hi|, tiny)`.
    ///
    /// # Errors
    /// [`FitError`] when the corpus is empty, oversized, non-finite, uses an
    /// invalid parameter identity, or does not have one exact shared non-empty
    /// parameter schema.
    pub fn fit(pairs: &[FidelityPair]) -> Result<Self, FitError> {
        let Some(first_pair) = pairs.first() else {
            return Err(FitError::EmptyTrainingSet);
        };
        let parameter_count = first_pair.params.len();
        validate_training_shape(pairs.len(), parameter_count)?;

        let mut bounds: BTreeMap<String, (f64, f64)> = BTreeMap::new();
        for (param, &value) in &first_pair.params {
            if let Some(reason) = color_identity_reason(param) {
                return Err(FitError::InvalidParameterIdentity {
                    pair_index: 0,
                    reason,
                });
            }
            bounds.insert(param.clone(), (value, value));
        }

        let mut mean_rel = 0.0_f64;
        let mut max_rel = 0.0_f64;
        for (pair_index, pair) in pairs.iter().enumerate() {
            if pair.params.len() != bounds.len() || !pair.params.keys().eq(bounds.keys()) {
                return Err(FitError::InconsistentParameterSchema { pair_index });
            }
            if !pair.lo_fi.is_finite() || !pair.hi_fi.is_finite() {
                return Err(FitError::NonFiniteTrainingQoi { pair_index });
            }
            for (param, &value) in &pair.params {
                if !value.is_finite() {
                    return Err(FitError::NonFiniteTrainingParameter {
                        pair_index,
                        param: param.clone(),
                    });
                }
                let Some((lo, hi)) = bounds.get_mut(param) else {
                    return Err(FitError::InconsistentParameterSchema { pair_index });
                };
                *lo = lo.min(value);
                *hi = hi.max(value);
            }

            let rel = (pair.hi_fi - pair.lo_fi).abs() / pair.hi_fi.abs().max(f64::MIN_POSITIVE);
            max_rel = max_rel.max(rel);
            if rel.is_infinite() {
                mean_rel = f64::INFINITY;
            } else if mean_rel.is_finite() {
                let count = (pair_index + 1) as f64;
                mean_rel += (rel - mean_rel) / count;
                mean_rel = mean_rel.min(max_rel);
            }
        }
        let mut observed = ValidityDomain::unconstrained();
        for (k, &(lo, hi)) in &bounds {
            observed = observed.with(k.clone(), lo, hi);
        }
        Ok(DiscrepancyModel {
            observed,
            trained_bounds: bounds,
            band: DiscrepancyBand { mean_rel, max_rel },
            pairs: pairs.len(),
        })
    }

    /// Number of training pairs.
    #[must_use]
    pub fn pairs(&self) -> usize {
        self.pairs
    }

    /// The observed (trained) parameter box.
    #[must_use]
    pub fn trained_domain(&self) -> &ValidityDomain {
        &self.observed
    }

    /// Query the band at `point`, refusing out-of-distribution use.
    ///
    /// # Errors
    /// [`OutOfDomain`] naming the first missing, unexpected, non-finite, or
    /// out-of-range parameter (BTreeMap order — deterministic diagnosis). The
    /// query key set must equal the training schema exactly; silently ignoring
    /// an extra physical dimension would make the in-domain claim unsound.
    pub fn query(&self, point: &BTreeMap<String, f64>) -> Result<DiscrepancyBand, OutOfDomain> {
        for (param, &(lo, hi)) in &self.trained_bounds {
            match point.get(param) {
                None => {
                    return Err(OutOfDomain {
                        param: param.clone(),
                        value: None,
                        trained: Some((lo, hi)),
                    });
                }
                Some(&v) if !v.is_finite() || v < lo || v > hi => {
                    return Err(OutOfDomain {
                        param: param.clone(),
                        value: Some(v),
                        trained: Some((lo, hi)),
                    });
                }
                Some(_) => {}
            }
        }
        if let Some((param, &value)) = point
            .iter()
            .find(|(param, _)| !self.trained_bounds.contains_key(param.as_str()))
        {
            return Err(OutOfDomain {
                param: param.clone(),
                value: Some(value),
                trained: None,
            });
        }
        Ok(self.band)
    }

    /// Model evidence for an in-domain use of the LOW-fidelity model,
    /// carrying the conservative (max) band.
    ///
    /// # Errors
    /// [`FitError::InvalidCardIdentity`] for an unusable evidence identity, or
    /// [`FitError::QueryOutOfDomain`] as in [`DiscrepancyModel::query`].
    pub fn evidence_at(
        &self,
        card_name: &str,
        point: &BTreeMap<String, f64>,
    ) -> Result<ModelEvidence, FitError> {
        if let Some(reason) = color_leaf_identity_reason(card_name) {
            return Err(FitError::InvalidCardIdentity { reason });
        }
        let band = self.query(point).map_err(FitError::QueryOutOfDomain)?;
        Ok(ModelEvidence {
            cards: vec![card_name.to_string()],
            assumptions: vec![format!(
                "lo-fi accuracy from a {}-pair two-fidelity discrepancy model",
                self.pairs
            )],
            validity: self.observed.clone(),
            discrepancy_rel: band.max_rel,
            in_domain: true,
        })
    }
}

/// Model bracketing: N plausible models of weakly-understood physics, one
/// QoI each; the evidence is the SPREAD, not a pretended point value.
#[derive(Debug, Clone, PartialEq)]
pub struct ModelBracket {
    members: Vec<(String, f64)>,
}

impl ModelBracket {
    /// Start a bracket.
    #[must_use]
    pub fn new() -> Self {
        ModelBracket {
            members: Vec::new(),
        }
    }

    fn push_member(&mut self, name: &str, qoi: f64) -> Result<(), FitError> {
        if self.members.len() >= MAX_BRACKET_MEMBERS {
            return Err(FitError::TooManyBracketMembers {
                maximum: MAX_BRACKET_MEMBERS,
            });
        }
        if let Some(reason) = color_leaf_identity_reason(name) {
            return Err(FitError::InvalidBracketMemberIdentity { reason });
        }
        if !qoi.is_finite() {
            return Err(FitError::NonFiniteBracketQoi {
                name: name.to_string(),
            });
        }
        if self.members.iter().any(|(member, _)| member == name) {
            return Err(FitError::DuplicateBracketMember {
                name: name.to_string(),
            });
        }
        self.members.push((name.to_string(), qoi));
        Ok(())
    }

    /// Add a member model's QoI and report admission failure immediately.
    ///
    /// # Errors
    /// [`FitError`] when the member identity or QoI is invalid, duplicated, or
    /// exceeds the bounded bracket-member budget.
    pub fn try_with_member(mut self, name: impl AsRef<str>, qoi: f64) -> Result<Self, FitError> {
        self.push_member(name.as_ref(), qoi)?;
        Ok(self)
    }

    /// Collapse the bracket into evidence: the numerical slice encloses
    /// every member's QoI (outward-rounded); the model slice records the
    /// bracketing and carries the relative spread as its band; the
    /// representative value is the member MIDRANGE (deterministic). At
    /// least two uniquely named finite members are required.
    ///
    /// # Errors
    /// [`FitError::TooFewBracketMembers`] when fewer than two valid models were
    /// supplied.
    pub fn evidence(&self, provenance: ProvenanceHash) -> Result<Evidence<f64>, FitError> {
        if self.members.len() < MIN_BRACKET_MEMBERS {
            return Err(FitError::TooFewBracketMembers {
                count: self.members.len(),
                minimum: MIN_BRACKET_MEMBERS,
            });
        }
        let mut qois = self.members.iter().map(|(_, qoi)| *qoi);
        let Some(first) = qois.next() else {
            return Err(FitError::TooFewBracketMembers {
                count: 0,
                minimum: MIN_BRACKET_MEMBERS,
            });
        };
        let (lo, hi) = qois.fold((first, first), |(lo, hi), qoi| (lo.min(qoi), hi.max(qoi)));
        let mid = f64::midpoint(lo, hi);
        let spread_rel = (hi - lo) / mid.abs().max(f64::MIN_POSITIVE);
        let mut names: Vec<String> = self.members.iter().map(|(n, _)| n.clone()).collect();
        names.sort_unstable();
        let mut sensitivity = SensitivitySummary::default();
        sensitivity
            .d_qoi
            .insert("model-choice(bracket-spread)".to_string(), hi - lo);
        let enclosure_lo = lo.next_down();
        let enclosure_hi = hi.next_up();
        Ok(Evidence {
            value: mid,
            qoi: mid,
            numerical: NumericalCertificate::enclosure(
                if enclosure_lo.is_finite() {
                    enclosure_lo
                } else {
                    lo
                },
                if enclosure_hi.is_finite() {
                    enclosure_hi
                } else {
                    hi
                },
            ),
            statistical: StatisticalCertificate::None,
            model: ModelEvidence {
                assumptions: vec![format!("model-bracketed over: {}", names.join(", "))],
                cards: names,
                validity: ValidityDomain::unconstrained(),
                discrepancy_rel: spread_rel,
                in_domain: true,
            },
            sensitivity,
            provenance,
            adjoint_ref: None,
        })
    }
}

impl Default for ModelBracket {
    fn default() -> Self {
        ModelBracket::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pt(pairs: &[(&str, f64)]) -> BTreeMap<String, f64> {
        pairs.iter().map(|&(k, v)| (k.to_string(), v)).collect()
    }

    fn bracket(members: &[(&str, f64)]) -> ModelBracket {
        members
            .iter()
            .try_fold(ModelBracket::new(), |bracket, &(name, qoi)| {
                bracket.try_with_member(name, qoi)
            })
            .expect("valid bracket fixture")
    }

    #[test]
    fn fit_refuses_empty_and_non_finite_training_sets() {
        let err = DiscrepancyModel::fit(&[]).expect_err("empty");
        assert!(matches!(&err, FitError::EmptyTrainingSet));
        assert!(err.to_string().contains("empty"), "{err}");
        let err = DiscrepancyModel::fit(&[FidelityPair {
            params: pt(&[("Re", 1e4)]),
            lo_fi: f64::NAN,
            hi_fi: 1.0,
        }])
        .expect_err("nan");
        assert!(err.to_string().contains("non-finite"), "{err}");
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let err = DiscrepancyModel::fit(&[FidelityPair {
                params: pt(&[("Re", value)]),
                lo_fi: 1.0,
                hi_fi: 1.0,
            }])
            .expect_err("non-finite parameter");
            assert!(err.to_string().contains("parameter"), "{err}");
        }
    }

    #[test]
    fn fit_requires_one_exact_bounded_parameter_schema() {
        let empty_schema = FidelityPair {
            params: BTreeMap::new(),
            lo_fi: 1.0,
            hi_fi: 1.0,
        };
        assert!(matches!(
            DiscrepancyModel::fit(&[empty_schema]),
            Err(FitError::EmptyParameterSchema)
        ));

        let inconsistent = [
            FidelityPair {
                params: pt(&[("Re", 1.0), ("Ma", 0.1)]),
                lo_fi: 1.0,
                hi_fi: 1.0,
            },
            FidelityPair {
                params: pt(&[("Re", 2.0)]),
                lo_fi: 1.0,
                hi_fi: 1.0,
            },
        ];
        assert!(matches!(
            DiscrepancyModel::fit(&inconsistent),
            Err(FitError::InconsistentParameterSchema { pair_index: 1 })
        ));

        for param in ["", " Re", "pending", "control\naxis"] {
            assert!(matches!(
                DiscrepancyModel::fit(&[FidelityPair {
                    params: pt(&[(param, 1.0)]),
                    lo_fi: 1.0,
                    hi_fi: 1.0,
                }]),
                Err(FitError::InvalidParameterIdentity { pair_index: 0, .. })
            ));
        }
        let oversized = "x".repeat(crate::MAX_COLOR_IDENTITY_BYTES + 1);
        assert!(matches!(
            DiscrepancyModel::fit(&[FidelityPair {
                params: pt(&[(oversized.as_str(), 1.0)]),
                lo_fi: 1.0,
                hi_fi: 1.0,
            }]),
            Err(FitError::InvalidParameterIdentity {
                pair_index: 0,
                reason: "too-long"
            })
        ));

        assert!(matches!(
            validate_training_shape(MAX_TRAINING_PAIRS + 1, 1),
            Err(FitError::TooManyTrainingPairs { .. })
        ));
        assert!(matches!(
            validate_training_shape(1, MAX_TRAINING_PARAMETERS + 1),
            Err(FitError::TooManyTrainingParameters { .. })
        ));
        assert!(matches!(
            validate_training_shape(
                MAX_TRAINING_PAIRS,
                MAX_TRAINING_COORDINATES / MAX_TRAINING_PAIRS + 1,
            ),
            Err(FitError::TooManyTrainingCoordinates { .. })
        ));
    }

    #[test]
    fn fit_mean_is_bounded_by_max_even_when_naive_sum_would_overflow() {
        let pairs = [
            FidelityPair {
                params: pt(&[("x", 0.0)]),
                lo_fi: -f64::MAX,
                hi_fi: 1.0,
            },
            FidelityPair {
                params: pt(&[("x", 1.0)]),
                lo_fi: -f64::MAX,
                hi_fi: 1.0,
            },
        ];
        let model = DiscrepancyModel::fit(&pairs).expect("finite extreme fit");
        let band = model.query(&pt(&[("x", 0.5)])).expect("in domain");
        assert!(band.mean_rel.is_finite());
        assert!(band.mean_rel <= band.max_rel);

        let unbounded = DiscrepancyModel::fit(&[
            FidelityPair {
                params: pt(&[("x", 0.0)]),
                lo_fi: -f64::MAX,
                hi_fi: f64::MAX,
            },
            FidelityPair {
                params: pt(&[("x", 1.0)]),
                lo_fi: 1.0,
                hi_fi: 1.0,
            },
        ])
        .expect("derived overflow is an honest unbounded band");
        let band = unbounded.query(&pt(&[("x", 0.5)])).expect("in domain");
        assert!(band.mean_rel.is_infinite() && band.max_rel.is_infinite());
    }

    #[test]
    fn in_domain_queries_report_the_band_and_out_of_domain_refuses() {
        let pairs: Vec<FidelityPair> = (0..10)
            .map(|i| {
                let re = 1e4 + f64::from(i) * 1e4;
                FidelityPair {
                    params: pt(&[("Re", re)]),
                    lo_fi: 1.0 + 0.05 * f64::from(i % 3),
                    hi_fi: 1.0,
                }
            })
            .collect();
        let model = DiscrepancyModel::fit(&pairs).expect("fit");
        let band = model.query(&pt(&[("Re", 5e4)])).expect("in domain");
        assert!(band.max_rel >= band.mean_rel && band.max_rel <= 0.2);
        let err = model.query(&pt(&[("Re", 1e6)])).expect_err("extrapolation");
        assert_eq!(err.param, "Re");
        assert!(err.to_string().contains("extrapolation"), "{err}");
        let err = model.query(&pt(&[("Ma", 0.1)])).expect_err("missing param");
        assert!(err.to_string().contains("not supplied"), "{err}");
        let err = model
            .query(&pt(&[("Re", 5e4), ("Mach", 0.1)]))
            .expect_err("an untrained query dimension must not be ignored");
        assert_eq!(err.param, "Mach");
        assert_eq!(err.trained, None);
        assert!(err.to_string().contains("exact training schema"), "{err}");
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let err = model
                .query(&pt(&[("Re", value)]))
                .expect_err("non-finite query");
            assert_eq!(err.param, "Re");
            assert_eq!(err.value.map(f64::to_bits), Some(value.to_bits()));
        }
        assert!(matches!(
            model.evidence_at("pending", &pt(&[("Re", 5e4)])),
            Err(FitError::InvalidCardIdentity {
                reason: "placeholder"
            })
        ));
        assert!(matches!(
            model.evidence_at("panel-vs-les", &pt(&[("Re", 1e6)])),
            Err(FitError::QueryOutOfDomain(OutOfDomain { ref param, .. }))
                if param == "Re"
        ));
    }

    #[test]
    fn brackets_enclose_every_member_and_report_the_spread() {
        let bracket = bracket(&[
            ("contact-angle-60", 0.90),
            ("contact-angle-90", 1.00),
            ("contact-angle-120", 1.16),
        ]);
        let ev = bracket
            .evidence(ProvenanceHash::of_bytes(b"vessel-lip"))
            .expect("nonempty bracket");
        assert!(ev.numerical.lo <= 0.90 && ev.numerical.hi >= 1.16);
        assert!(
            ev.model.discrepancy_rel > 0.2,
            "{}",
            ev.model.discrepancy_rel
        );
        assert!(ev.model.assumptions[0].contains("model-bracketed"));
        assert!(matches!(
            crate::color_of(&ev.numerical, &ev.model),
            crate::Color::Estimated { dispersion, .. }
                if dispersion.to_bits() == ev.model.discrepancy_rel.to_bits()
        ));
        assert!(matches!(
            ModelBracket::new().evidence(ProvenanceHash(0)),
            Err(FitError::TooFewBracketMembers {
                count: 0,
                minimum: MIN_BRACKET_MEMBERS
            })
        ));
        let single = ModelBracket::new()
            .try_with_member("only-model", 1.0)
            .expect("valid member");
        assert!(matches!(
            single.evidence(ProvenanceHash(0)),
            Err(FitError::TooFewBracketMembers { count: 1, .. })
        ));
    }

    #[test]
    fn bracket_refuses_non_finite_duplicate_and_invalid_members() {
        for qoi in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                ModelBracket::new().try_with_member("bad-model", qoi),
                Err(FitError::NonFiniteBracketQoi { ref name }) if name == "bad-model"
            ));
        }
        let one_member = ModelBracket::new()
            .try_with_member("same-model", 1.0)
            .expect("first identity");
        assert!(matches!(
            one_member.try_with_member("same-model", 2.0),
            Err(FitError::DuplicateBracketMember { ref name }) if name == "same-model"
        ));
        for name in ["", " pending", "pending", "derived:v2:forged"] {
            assert!(matches!(
                ModelBracket::new().try_with_member(name, 1.0),
                Err(FitError::InvalidBracketMemberIdentity { .. })
            ));
        }
        let two_members = bracket(&[("a", 1.0), ("b", 2.0)]);
        assert!(matches!(
            two_members.try_with_member("b", 3.0),
            Err(FitError::DuplicateBracketMember { .. })
        ));
        let full = ModelBracket {
            members: (0..MAX_BRACKET_MEMBERS)
                .map(|index| (format!("model-{index}"), index as f64))
                .collect(),
        };
        assert!(matches!(
            full.try_with_member("one-too-many", 0.0),
            Err(FitError::TooManyBracketMembers {
                maximum: MAX_BRACKET_MEMBERS
            })
        ));
    }

    #[test]
    fn finite_extreme_bracket_keeps_finite_numerical_evidence() {
        let evidence = ModelBracket::new()
            .try_with_member("negative-extreme", -f64::MAX)
            .expect("first member")
            .try_with_member("positive-extreme", f64::MAX)
            .expect("second member")
            .evidence(ProvenanceHash(0))
            .expect("finite bracket");
        assert!(evidence.numerical.lo.is_finite());
        assert!(evidence.numerical.hi.is_finite());
        assert!(evidence.model.discrepancy_rel.is_infinite());
        assert!(evidence.clone().certified().is_ok());
        assert!(evidence.breakdown().model_rel.is_infinite());
    }
}
