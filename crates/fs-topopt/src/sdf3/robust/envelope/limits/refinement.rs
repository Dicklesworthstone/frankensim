//! Finer-grid checks of every authored eroded compliance requirement.
//! The probe inherits physical stiffness, not refiltered fine raw controls.
use std::{collections::BTreeMap, ops::ControlFlow};
use fs_cutfem::elastic3::{
    adaptive::{AdaptiveElasticity3, enrichment::AdaptiveTransfer3},
    surface::ReferenceLoad3,
};
use fs_cutfem::octree3::Octant3;
use fs_dwr::elasticity3::{GoalError3, GoalMarking3, dorfler3};
use super::ComplianceCaps3;
use crate::sdf3::{AdaptiveSdf3Elasticity, CutDensityStudy3};
use crate::sdf3::robust::{RobustOptions3, goal::RobustComplianceRefinement3};
use crate::sdf3_goal::{GoalReferenceLoad3, GoalRefinementError3, GoalRefinementOptions3};
use crate::{EvaluationStop, SimpParams, SolveControl, SolveWork};

/// Physical compliance values and their cap in original load order.
/// Normalized differences measure two-grid change, NOT remaining error.
#[derive(Debug, Clone, PartialEq)]
pub struct ComplianceCapCheck3 {
    pub cap: f64,
    pub coarse: f64,
    pub fine: f64,
    pub normalized_change: f64,
}
impl ComplianceCapCheck3 {
    #[must_use]
    pub fn satisfied(&self) -> bool {
        self.cap.is_finite() && self.cap > 0.0
            && [self.coarse, self.fine].iter().all(|c| c.is_finite() && *c >= 0.0 && *c <= self.cap)
    }
}

#[derive(Debug, Clone)]
pub struct ComplianceCapRefinement3 {
    pub caps: ComplianceCaps3,
    /// Complete eroded solve and two-grid decomposition for EACH load.
    /// Its uniform-mean diagnostic is not the cap-admission criterion.
    pub assessment: RobustComplianceRefinement3,
    pub checks: Vec<ComplianceCapCheck3>,
    pub coarse_envelope: f64,
    pub fine_envelope: f64,
    /// Sum of per-case absolute local DWR masses divided by that case's cap.
    /// Different physical scales and opposing signed corrections cannot hide
    /// a requirement behind a weighted-average compliance.
    pub marking_mass: BTreeMap<Octant3, f64>,
    pub work: SolveWork,
}
impl ComplianceCapRefinement3 {
    /// Both grids meet EVERY physical cap and the smooth constraint, with all
    /// cap-normalized changes within the separate requested tolerance.
    /// No continuum, stationarity, stress, or dilated-volume error claim.
    #[must_use]
    pub fn stable_and_feasible(&self, smooth_tolerance: f64, change_tolerance: f64) -> bool {
        smooth_tolerance.is_finite() && (0.0..1.0).contains(&smooth_tolerance)
            && change_tolerance.is_finite() && change_tolerance >= 0.0
            && self.checks.len() == self.caps.caps().len()
            && self.coarse_envelope.is_finite() && self.fine_envelope.is_finite()
            && self.coarse_envelope - 1.0 <= smooth_tolerance
            && self.fine_envelope - 1.0 <= smooth_tolerance
            && self.checks.iter().zip(self.caps.caps()).all(|(c, cap)|
                c.cap == *cap && c.satisfied() && c.normalized_change.is_finite()
                && c.normalized_change >= 0.0 && c.normalized_change <= change_tolerance)
    }
    pub fn mark(&self, fraction: f64, max_marks: usize, checkpoint: impl FnMut() -> ControlFlow<()>)
        -> Result<GoalMarking3, GoalError3> {
        dorfler3(&self.marking_mass, fraction, max_marks, checkpoint)
    }
}

fn invalid(message: &'static str) -> GoalRefinementError3 { GoalRefinementError3::Invalid(message) }

impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Re-solve the current design on its explicit eroded model and on a
    /// genuinely enriched background. Every reference law is a constraint:
    /// there are no zero weights, inferred caps, or dropped load cases.
    ///
    /// Fine PHYSICAL stiffness is parent-inherited. Fine raw filtering and
    /// optimization are a DIFFERENT experiment. Geometry, base material,
    /// supports, integration policy and pure density-independent reference
    /// laws must describe the same problem. The transfer checks compatibility,
    /// but does not certify SDF identity. Homogeneous supports only.
    ///
    /// All physics and the extra transfer share the supplied control. Success
    /// and returned errors leave incoming nominal parameters/scales unchanged.
    /// Failed enrichment returns no partially checked load family.
    #[allow(clippy::too_many_arguments)]
    pub fn estimate_compliance_cap_enrichment(
        &mut self, enriched: AdaptiveElasticity3, laws: &[ReferenceLoad3<'_>],
        rho: &[f64], params: SimpParams, robustness: RobustOptions3, caps: &ComplianceCaps3,
        options: GoalRefinementOptions3, control: &mut SolveControl<'_>,
    ) -> Result<ComplianceCapRefinement3, GoalRefinementError3> {
        self.require_material_body_recipe("sdf3-cap-enrichment-requires-material-acceleration")?;
        if laws.len() != caps.caps().len() || laws.len() > options.max_load_cases {
            return Err(invalid("one reference load per compliance cap is required"));
        }
        robustness.assert_valid(params);
        if rho.len() != self.cells() || rho.iter().any(|r| !r.is_finite() || !(0.0..=1.0).contains(r)) {
            return Err(invalid("invalid compliance-cap raw design"));
        }
        control.checkpoint("sdf3-cap-enrichment-start")?;
        let mut stop = None;
        let transfer = AdaptiveTransfer3::new(self.operator().adaptive(), &enriched, options.max_transfer_terms, || {
            match control.checkpoint("sdf3-cap-enrichment-transfer") {
                Ok(()) => ControlFlow::Continue(()),
                Err(error) => { stop = Some(error); ControlFlow::Break(()) }
            }
        });
        if let Some(error) = stop { return Err(error.into()); }
        let transfer = transfer?;
        if !enriched.leaves().iter().zip(transfer.parents()).any(|(leaf, &parent)|
            leaf.level() > self.operator().adaptive().leaves()[parent].level()) {
            return Err(invalid("compliance-cap assessment requires actual active-space enrichment"));
        }
        drop(transfer);
        let loads: Vec<_> = laws.iter().map(|&load| GoalReferenceLoad3 {
            load, weight: 1.0 / laws.len() as f64,
        }).collect();
        let assessment = self.estimate_robust_reference_compliance_enrichment(enriched, &loads, rho,
            params, robustness, options, control)?;
        let coarse: Vec<_> = assessment.estimate.cases.iter().map(|c| c.coarse_value).collect();
        let fine: Vec<_> = assessment.estimate.cases.iter().map(|c| c.fine_value).collect();
        let coarse_envelope = caps.envelope().aggregate(&coarse)?.value;
        let fine_envelope = caps.envelope().aggregate(&fine)?.value;
        let mut checks = Vec::with_capacity(laws.len());
        let mut marking_mass: BTreeMap<Octant3, f64> = BTreeMap::new();
        for ((case, &cap), (&coarse, &fine)) in assessment.estimate.cases.iter().zip(caps.caps())
            .zip(coarse.iter().zip(&fine)) {
            control.checkpoint("sdf3-cap-enrichment-check")?;
            let normalized_change = (fine - coarse).abs() / cap;
            if !normalized_change.is_finite() { return Err(invalid("unrepresentable compliance-cap change")); }
            checks.push(ComplianceCapCheck3 { cap, coarse, fine, normalized_change });
            for (&cell, local) in &case.cells {
                control.checkpoint("sdf3-cap-enrichment-marking")?;
                let scaled = local.marking_mass / cap;
                let value = marking_mass.entry(cell).or_insert(0.0);
                *value += scaled;
                if !scaled.is_finite() || scaled < 0.0 || !value.is_finite() {
                    return Err(EvaluationStop::Breakdown { stage: "sdf3-cap-enrichment-marking" }.into());
                }
            }
        }
        control.checkpoint("sdf3-cap-enrichment-publish")?;
        Ok(ComplianceCapRefinement3 { caps: caps.clone(), assessment, checks,
            coarse_envelope, fine_envelope, marking_mass, work: control.work() })
    }
}

#[cfg(test)]
#[path = "refinement_tests.rs"]
mod tests;
