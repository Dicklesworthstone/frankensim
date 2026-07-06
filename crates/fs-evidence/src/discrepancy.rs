//! Discrepancy models and model bracketing (patch Rev B mechanisms 2–3).
//!
//! A [`DiscrepancyModel`] is fit from paired two-fidelity evaluations (the
//! ledger's accumulating corpus) and answers "how wrong is the cheap model
//! HERE" — refusing with a teaching [`OutOfDomain`] when asked outside the
//! region it has data for (the surrogate out-of-distribution guard). v1 is
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
    StatisticalCertificate, ValidityDomain,
};
use core::fmt;
use std::collections::BTreeMap;

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

/// Structured fit failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FitError {
    /// What was wrong with the training set.
    pub detail: &'static str,
}

impl fmt::Display for FitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "discrepancy fit refused: {}; supply at least one (lo-fi, hi-fi) pair with finite \
             values",
            self.detail
        )
    }
}

impl core::error::Error for FitError {}

/// The out-of-distribution refusal, naming the violated parameter (the
/// diagnosis an agent needs to decide between gathering data and
/// escalating fidelity).
#[derive(Debug, Clone, PartialEq)]
pub struct OutOfDomain {
    /// The parameter outside the trained box (or missing from the query).
    pub param: String,
    /// The queried value (`None` = the query omitted the parameter).
    pub value: Option<f64>,
    /// The trained box for that parameter.
    pub trained: (f64, f64),
}

impl fmt::Display for OutOfDomain {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.value {
            Some(v) => write!(
                f,
                "discrepancy query out of domain: `{}` = {v} lies outside the trained box \
                 [{}, {}] — the band would be extrapolation; gather pairs there or escalate \
                 fidelity",
                self.param, self.trained.0, self.trained.1
            ),
            None => write!(
                f,
                "discrepancy query out of domain: `{}` was not supplied but the model was \
                 trained on it (trained box [{}, {}])",
                self.param, self.trained.0, self.trained.1
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

impl DiscrepancyModel {
    /// Fit from paired evaluations. The observed box is the per-parameter
    /// min/max over training points; the band is mean/max of
    /// `|hi - lo| / max(|hi|, tiny)`.
    ///
    /// # Errors
    /// [`FitError`] on an empty or non-finite training set.
    pub fn fit(pairs: &[FidelityPair]) -> Result<Self, FitError> {
        if pairs.is_empty() {
            return Err(FitError {
                detail: "the training set is empty",
            });
        }
        let mut bounds: BTreeMap<String, (f64, f64)> = BTreeMap::new();
        let mut rels = Vec::with_capacity(pairs.len());
        for p in pairs {
            if !p.lo_fi.is_finite() || !p.hi_fi.is_finite() {
                return Err(FitError {
                    detail: "a training pair has a non-finite QoI",
                });
            }
            for (k, &v) in &p.params {
                bounds
                    .entry(k.clone())
                    .and_modify(|(lo, hi)| {
                        *lo = lo.min(v);
                        *hi = hi.max(v);
                    })
                    .or_insert((v, v));
            }
            rels.push((p.hi_fi - p.lo_fi).abs() / p.hi_fi.abs().max(f64::MIN_POSITIVE));
        }
        let mean_rel = rels.iter().sum::<f64>() / rels.len() as f64;
        let max_rel = rels.iter().copied().fold(0.0f64, f64::max);
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
    /// [`OutOfDomain`] naming the first violated parameter (BTreeMap
    /// order — deterministic diagnosis).
    pub fn query(&self, point: &BTreeMap<String, f64>) -> Result<DiscrepancyBand, OutOfDomain> {
        for (param, &(lo, hi)) in &self.trained_bounds {
            match point.get(param) {
                None => {
                    return Err(OutOfDomain {
                        param: param.clone(),
                        value: None,
                        trained: (lo, hi),
                    });
                }
                Some(&v) if v < lo || v > hi => {
                    return Err(OutOfDomain {
                        param: param.clone(),
                        value: Some(v),
                        trained: (lo, hi),
                    });
                }
                Some(_) => {}
            }
        }
        Ok(self.band)
    }

    /// Model evidence for an in-domain use of the LOW-fidelity model,
    /// carrying the conservative (max) band.
    ///
    /// # Errors
    /// [`OutOfDomain`] as in [`DiscrepancyModel::query`].
    pub fn evidence_at(
        &self,
        card_name: &str,
        point: &BTreeMap<String, f64>,
    ) -> Result<ModelEvidence, OutOfDomain> {
        let band = self.query(point)?;
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

    /// Add a member model's QoI.
    #[must_use]
    pub fn with_member(mut self, name: impl Into<String>, qoi: f64) -> Self {
        self.members.push((name.into(), qoi));
        self
    }

    /// Collapse the bracket into evidence: the numerical slice encloses
    /// every member's QoI (outward-rounded); the model slice records the
    /// bracketing and carries the relative spread as its band; the
    /// representative value is the member MIDRANGE (deterministic). At
    /// least one member is required — the teaching `None` otherwise.
    #[must_use]
    pub fn evidence(&self, provenance: ProvenanceHash) -> Option<Evidence<f64>> {
        if self.members.is_empty() {
            return None;
        }
        let lo = self
            .members
            .iter()
            .map(|&(_, q)| q)
            .fold(f64::INFINITY, f64::min);
        let hi = self
            .members
            .iter()
            .map(|&(_, q)| q)
            .fold(f64::NEG_INFINITY, f64::max);
        let mid = 0.5 * (lo + hi);
        let spread_rel = (hi - lo) / mid.abs().max(f64::MIN_POSITIVE);
        let mut names: Vec<String> = self.members.iter().map(|(n, _)| n.clone()).collect();
        names.sort_unstable();
        let mut sensitivity = SensitivitySummary::default();
        sensitivity
            .d_qoi
            .insert("model-choice(bracket-spread)".to_string(), hi - lo);
        Some(Evidence {
            value: mid,
            qoi: mid,
            numerical: NumericalCertificate::enclosure(lo.next_down(), hi.next_up()),
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

    #[test]
    fn fit_refuses_empty_and_non_finite_training_sets() {
        let err = DiscrepancyModel::fit(&[]).expect_err("empty");
        assert!(err.to_string().contains("empty"), "{err}");
        let err = DiscrepancyModel::fit(&[FidelityPair {
            params: pt(&[("Re", 1e4)]),
            lo_fi: f64::NAN,
            hi_fi: 1.0,
        }])
        .expect_err("nan");
        assert!(err.to_string().contains("non-finite"), "{err}");
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
    }

    #[test]
    fn brackets_enclose_every_member_and_report_the_spread() {
        let bracket = ModelBracket::new()
            .with_member("contact-angle-60", 0.90)
            .with_member("contact-angle-90", 1.00)
            .with_member("contact-angle-120", 1.16);
        let ev = bracket
            .evidence(ProvenanceHash::of_bytes(b"vessel-lip"))
            .expect("nonempty bracket");
        assert!(ev.numerical.lo <= 0.90 && ev.numerical.hi >= 1.16);
        assert!(ev.model.discrepancy_rel > 0.2, "{}", ev.model.discrepancy_rel);
        assert!(ev.model.assumptions[0].contains("model-bracketed"));
        assert!(ModelBracket::new().evidence(ProvenanceHash(0)).is_none());
    }
}
