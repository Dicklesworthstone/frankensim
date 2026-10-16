//! DWR screening of the declared weighted compliance objective, after independent
//! equilibrium solves. These public values are numerical data, not certificates.
//! Stress remains an unweighted, sample-scoped constraint of the original owner.

use fs_cutfem::CutFemError;
use crate::goal_resolution::{DwrGainDecision, DwrGatePolicy, DwrMeasurement, compare_dwr_gains};
use crate::robust_descent::{MultiLoadProjectedProgress, MultiLoadProjectedStage, MultiLoadProjectedState};
use crate::{RobustLoadCase, WeightedComplianceDwrAssessment, WeightedComplianceDwrStage};

fn invalid(what: impl Into<String>) -> CutFemError {
    CutFemError::InvalidElasticityInput { what: what.into() }
}

/// Numerical sensitivity allowed relative to the proposed aggregate gain.
/// Extra DWR solves consume the optimizer's EXISTING lifetime `max_solves`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MultiLoadDwrPolicy {
    /// Strictly between zero and one. Passing does not certify continuum error.
    pub max_error_fraction: f64,
}
impl MultiLoadDwrPolicy {
    /// Admit the canonical two-level estimator's supported geometry/load family.
    pub fn validate(self, level: u32, cases: usize) -> Result<(), CutFemError> {
        if !(1..=5).contains(&level) || !(1..=16).contains(&cases)
            || !self.max_error_fraction.is_finite()
            || self.max_error_fraction <= 0.0 || self.max_error_fraction >= 1.0 {
            return Err(invalid("weighted DWR checks require level 1..=5, 1..=16 cases and an error fraction in (0,1)"));
        }
        Ok(())
    }
}

/// Every load's actual coarse/enriched solve on the same unchanged geometry.
/// Keeping the load declaration beside each row prevents reorder/weight aliases.
#[derive(Debug, Clone, PartialEq)]
pub struct MultiLoadDwrMeasurement {
    /// Declaration-order data, including zero-weight operating conditions.
    pub cases: Vec<(RobustLoadCase, DwrMeasurement)>,
}
impl MultiLoadDwrMeasurement {
    /// Reconstruct the unnormalized weighted objective from complete finite data.
    /// Caller-constructed values do not acquire solver/provenance authority here.
    pub fn aggregate(&self, loads: &[RobustLoadCase]) -> Result<DwrMeasurement, CutFemError> {
        if !(1..=16).contains(&loads.len()) || self.cases.len() != loads.len()
            || !loads.iter().any(|load| load.weight() > 0.0) {
            return Err(invalid("weighted DWR requires a complete load family and a positive objective weight"));
        }
        let first = self.cases[0].1;
        let mut total = DwrMeasurement { coarse_compliance: 0.0, enriched_compliance: 0.0,
            eta_signed: 0.0, coarse_relative_residual: 0.0, enriched_relative_residual: 0.0, ..first };
        for ((declared, measured), expected) in self.cases.iter().zip(loads) {
            measured.validate()?;
            if declared != expected || measured.level != first.level
                || measured.snapshot != first.snapshot || measured.area.to_bits() != first.area.to_bits() {
                return Err(invalid("weighted DWR mixes loads, levels or geometry identities"));
            }
            let weight = declared.weight();
            total.coarse_compliance += weight * measured.coarse_compliance;
            total.enriched_compliance += weight * measured.enriched_compliance;
            total.eta_signed += weight * measured.eta_signed;
            total.coarse_relative_residual = total.coarse_relative_residual.max(measured.coarse_relative_residual);
            total.enriched_relative_residual = total.enriched_relative_residual.max(measured.enriched_relative_residual);
        }
        total.validate()?;
        Ok(total)
    }

    /// Bind EVERY coarse response, not merely its aggregate, to existing mechanics.
    pub fn matches_state(&self, loads: &[RobustLoadCase], state: &MultiLoadProjectedState,
        level: u32) -> Result<(), CutFemError>
    {
        let total = self.aggregate(loads)?;
        if total.level != level || total.snapshot != state.snapshot
            || total.area.to_bits() != state.volume.to_bits()
            || total.coarse_compliance.to_bits() != state.objective.to_bits()
            || state.active_case.is_some() || state.case_compliances.len() != loads.len()
            || self.cases.iter().zip(&state.case_compliances)
                .any(|((_, m), c)| m.coarse_compliance.to_bits() != c.to_bits()) {
            return Err(invalid("weighted DWR coarse solves do not reproduce the exact accepted/projected load mechanics"));
        }
        Ok(())
    }

    pub(crate) fn from_assessment(a: WeightedComplianceDwrAssessment, loads: &[RobustLoadCase],
        state: &MultiLoadProjectedState, level: u32) -> Result<Self, CutFemError>
    {
        let measured = Self { cases: a.cases.into_iter().map(|case| (case.load, DwrMeasurement {
            snapshot: a.snapshot, level: a.level, area: a.volume,
            coarse_compliance: case.estimate.j_primal, enriched_compliance: case.estimate.j_enriched,
            eta_signed: case.estimate.eta_signed,
            coarse_relative_residual: case.estimate.primal_relative_residual,
            enriched_relative_residual: case.estimate.enriched_relative_residual,
        })).collect() };
        measured.matches_state(loads, state, level)?;
        let total = measured.aggregate(loads)?;
        if total.coarse_compliance.to_bits() != a.coarse_compliance.to_bits()
            || total.enriched_compliance.to_bits() != a.enriched_compliance.to_bits()
            || total.eta_signed.to_bits() != a.eta_signed.to_bits() {
            return Err(invalid("weighted DWR aggregate disagrees with its independent cases"));
        }
        Ok(measured)
    }
}

/// Arithmetic of an Estimated acceptance decision. Individual objectives may
/// trade off; opposite per-case error changes cannot hide behind an aggregate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MultiLoadDwrDecision {
    pub aggregate: DwrGainDecision,
    /// Sum_i w_i * |eta_i(baseline) - eta_i(candidate)|.
    pub weighted_correction_variation: f64,
    /// Sum_i w_i * |fine_gain_i - coarse_gain_i|.
    pub weighted_gain_variation: f64,
    pub refusal: Option<&'static str>,
}

/// Apply weights AFTER solving each case and BEFORE combining numerical changes.
/// This is conservative screening, not a proof of continuum improvement; shared
/// bias can survive. Absolute cell-indicator mass is never an interval radius.
pub fn compare_weighted_dwr_gains(baseline: &MultiLoadDwrMeasurement,
    candidate: &MultiLoadDwrMeasurement, loads: &[RobustLoadCase],
    policy: MultiLoadDwrPolicy, decrease: f64) -> Result<MultiLoadDwrDecision, CutFemError>
{
    let old = baseline.aggregate(loads)?;
    let new = candidate.aggregate(loads)?;
    policy.validate(old.level, loads.len())?;
    // Reuse the single scalar arithmetic owner. Four is its representation
    // admission minimum, NOT a work grant; real work uses the lifetime budget.
    let aggregate = compare_dwr_gains(old, new, DwrGatePolicy {
        max_error_fraction: policy.max_error_fraction, max_solves_per_attempt: 4,
    }, decrease)?;
    let mut correction = 0.0;
    let mut gain_change = 0.0;
    for ((load, a), (_, b)) in baseline.cases.iter().zip(&candidate.cases) {
        // Still validate zero-weight rows above, but do not evaluate 0 * overflow.
        let w = load.weight();
        if w == 0.0 { continue; }
        correction += w * (a.eta_signed - b.eta_signed).abs();
        gain_change += w * ((a.enriched_compliance - b.enriched_compliance)
            - (a.coarse_compliance - b.coarse_compliance)).abs();
    }
    if !correction.is_finite() || !gain_change.is_finite() {
        return Err(invalid("weighted DWR gain variation overflow"));
    }
    let refusal = aggregate.refusal.or_else(||
        (correction.max(gain_change) > aggregate.allowance)
            .then_some("independent-load-discretization-dominates-gain"));
    Ok(MultiLoadDwrDecision { aggregate, weighted_correction_variation: correction,
        weighted_gain_variation: gain_change, refusal })
}

/// Existing optimizer boundaries plus whole-case DWR cancellation points.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MultiLoadDwrStage {
    Prepare,
    Baseline(WeightedComplianceDwrStage),
    Candidate { index: usize, stage: WeightedComplianceDwrStage },
    Optimizer(MultiLoadProjectedStage),
}
#[derive(Debug, Clone)]
pub struct MultiLoadDwrCandidate {
    pub index: usize,
    pub measurement: Option<MultiLoadDwrMeasurement>,
    pub decision: Option<MultiLoadDwrDecision>,
    pub refusal: Option<String>,
}
#[derive(Debug, Clone)]
pub enum MultiLoadDwrProgress {
    IterationLimit,
    SolveBudget,
    Searched {
        baseline: MultiLoadDwrMeasurement,
        candidates: Vec<MultiLoadDwrCandidate>,
        progress: MultiLoadProjectedProgress,
    },
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DesignBoxEdge;
    fn load(weight: f64, vertical: f64) -> RobustLoadCase {
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, vertical], weight).unwrap()
    }
    fn measured(loads: &[RobustLoadCase], coarse: [f64; 2], fine: [f64; 2], eta: [f64; 2]) -> MultiLoadDwrMeasurement {
        MultiLoadDwrMeasurement { cases: loads.iter().enumerate().map(|(i, &load)| (load,
            DwrMeasurement { snapshot: 7, level: 3, area: 0.6, coarse_compliance: coarse[i],
                enriched_compliance: fine[i], eta_signed: eta[i],
                coarse_relative_residual: 1e-14, enriched_relative_residual: 1e-14 })).collect() }
    }
    fn policy() -> MultiLoadDwrPolicy { MultiLoadDwrPolicy { max_error_fraction: 0.5 } }

    #[test]
    fn weighted_dwr_opposite_error_changes_cannot_cancel_into_an_acceptance() {
        let loads = [load(1.0, -1.0), load(1.0, 1.0)];
        let old = measured(&loads, [10.0; 2], [10.0; 2], [0.0; 2]);
        let new = measured(&loads, [9.0; 2], [10.0, 8.0], [1.0, -1.0]);
        let decision = compare_weighted_dwr_gains(&old, &new, &loads, policy(), 0.0).unwrap();
        assert!(decision.aggregate.refusal.is_none());
        assert_eq!(decision.aggregate.signed_correction_change, 0.0);
        assert_eq!(decision.weighted_correction_variation, 2.0);
        assert_eq!(decision.weighted_gain_variation, 2.0);
        assert_eq!(decision.refusal, Some("independent-load-discretization-dominates-gain"));
        let consistent = measured(&loads, [9.0; 2], [9.0; 2], [0.0; 2]);
        assert!(compare_weighted_dwr_gains(&old, &consistent, &loads, policy(), 0.0).unwrap().refusal.is_none());
    }

    #[test]
    fn weighted_dwr_allows_objective_tradeoffs_without_normalizing_weights() {
        let loads = [load(2.0, -1.0), load(0.5, 1.0)];
        let old = measured(&loads, [10.0; 2], [11.0; 2], [1.0; 2]);
        let new = measured(&loads, [9.0, 11.0], [10.0, 12.0], [1.0; 2]);
        let decision = compare_weighted_dwr_gains(&old, &new, &loads, policy(), 0.0).unwrap();
        assert_eq!(old.aggregate(&loads).unwrap().coarse_compliance, 25.0);
        assert_eq!(decision.aggregate.coarse_gain, 1.5);
        assert!(decision.refusal.is_none());
        let reversed = measured(&loads, [9.0, 11.0], [12.0, 13.0], [1.0; 2]);
        assert_eq!(compare_weighted_dwr_gains(&old, &reversed, &loads, policy(), 0.0).unwrap().refusal,
            Some("enriched-decrease-not-met"));
    }

    #[test]
    fn weighted_dwr_binds_each_case_and_validates_zero_weight_evidence() {
        let loads = [load(1.0, -1.0), load(0.0, 4.0)];
        let old = measured(&loads, [10.0, 160.0], [11.0, 176.0], [1.0, 16.0]);
        let state = MultiLoadProjectedState { objective: 10.0, case_compliances: vec![10.0, 160.0],
            active_case: None, volume: 0.6, snapshot: 7 };
        assert!(old.matches_state(&loads, &state, 3).is_ok());
        let mut bad = old.clone(); bad.cases[1].1.coarse_compliance += 1.0;
        assert!(bad.matches_state(&loads, &state, 3).is_err());
        for mutate in 0..4 {
            let mut bad = old.clone();
            match mutate {
                0 => bad.cases.swap(0, 1),
                1 => { bad.cases.pop(); },
                2 => bad.cases[1].1.enriched_relative_residual = f64::NAN,
                _ => bad.cases[1].1.snapshot += 1,
            }
            assert!(bad.aggregate(&loads).is_err());
        }
    }
}
