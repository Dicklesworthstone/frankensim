//! Two-grid verification of the same binary-eroded physical material.
pub mod scan;
use super::*;
use crate::sdf3::robust::envelope::limits::{ComplianceCaps3, binary::{
    BinaryComplianceOptions3, BinaryComplianceEvaluation3, BinaryComplianceError3,
}};

#[derive(Debug, Clone)]
pub struct BinaryLoadSpaceCapRefinement3 {
    pub cap: LoadSpaceCap3,
    pub rho: Vec<f64>,
    /// Coarse masks and freshly solved binary-eroded basis fields. The embedded
    /// finite-case envelope is diagnostic, not the full load-set requirement.
    pub physical: BinaryComplianceEvaluation3,
    pub coarse: LoadSpaceGridCheck3,
    pub fine: LoadSpaceGridCheck3,
    pub normalized_gram_change: f64,
    pub basis_assessment: ComplianceRefinement3,
    pub work: SolveWork,
}
impl BinaryLoadSpaceCapRefinement3 {
    /// Exactly the same physical and full-matrix change checks as the gray
    /// two-grid assessment. No volume, topology or zero-void limit is checked.
    #[must_use]
    pub fn stable_and_feasible(&self, smooth_tolerance: f64, change_tolerance: f64) -> bool {
        pair_admits(&self.cap, &self.coarse, &self.fine, self.normalized_gram_change,
            smooth_tolerance, change_tolerance)
    }
}
fn binary_error(error: BinaryComplianceError3) -> GoalRefinementError3 {
    match error {
        BinaryComplianceError3::Invalid(message) => invalid(message),
        BinaryComplianceError3::Evaluation(stop) => stop.into(),
    }
}

impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Freshly threshold and solve the COARSE binary material, then inherit its
    /// physical solid/void stiffness on the enriched operator. Fine raw controls
    /// are never filtered and the fine field is never rethresholded. All cross
    /// compliances and worst directions are recomputed independently per grid.
    ///
    /// Reuses binary elasticity and the shared reference-load/DWR solve owner.
    /// Params explicitly identify the gray source of the thresholded masks; the
    /// binary void stiffness is the independently declared endpoint option.
    /// Protected regions, fixed supports and loads in voids retain their original
    /// semantics. This does not demonstrate a zero-void structural limit.
    ///
    /// Same geometry/load-law obligations as estimate_load_space_cap_enrichment.
    /// All work uses control and every returned outcome restores incoming study
    /// state. No cached caller-supplied fields are trusted or partly published.
    #[allow(clippy::too_many_arguments)]
    pub fn estimate_binary_load_space_cap_enrichment(&mut self, enriched: AdaptiveElasticity3,
        laws: &[ReferenceLoad3<'_>], rho: &[f64], params: SimpParams, cap: &LoadSpaceCap3,
        binary: BinaryComplianceOptions3, options: GoalRefinementOptions3,
        control: &mut SolveControl<'_>) -> Result<BinaryLoadSpaceCapRefinement3, GoalRefinementError3> {
        // Refuse malformed endpoint controls before invoking any user load law.
        if !binary.threshold.is_finite() || binary.threshold <= 0.0 || binary.threshold > 1.0
            || !binary.void_stiffness.is_finite() || binary.void_stiffness <= 0.0 || binary.void_stiffness >= 1.0
            || laws.len() > binary.max_load_cases {
            return Err(invalid("invalid binary load-space enrichment endpoint policy"));
        }
        admit_enrichment(self, &enriched, laws, rho, params, binary.robustness, cap, options, control)?;
        let caps = ComplianceCaps3::new(vec![cap.cap(); laws.len()], cap.policy().temperature(), laws.len())
            .map_err(invalid)?;
        let incoming_params = self.params;
        let incoming_scales = self.operator().scales().to_vec();
        let result = (|| {
            let mut forces = Vec::with_capacity(laws.len());
            for &law in laws {
                forces.push(self.operator().adaptive().reference_load(law,
                    || poll(control, "sdf3-binary-load-space-grid-load"))?);
            }
            let basis: Vec<_> = forces.iter().map(Vec::as_slice).collect();
            self.params = params;
            let physical = self.evaluate_binary_compliance_caps(rho, &basis, &caps, binary, control)
                .map_err(binary_error)?;
            self.operator.set_scales(&physical.eroded_stiffness_scales())?;
            let (coarse, fine, normalized_gram_change, basis_assessment) =
                self.assess_load_space_grids(enriched, laws, &forces, &physical.displacements,
                    cap, options, control)?;
            control.checkpoint("sdf3-binary-load-space-grid-publish")?;
            Ok(BinaryLoadSpaceCapRefinement3 { cap: cap.clone(), rho: rho.to_vec(), physical,
                coarse, fine, normalized_gram_change, basis_assessment, work: control.work() })
        })();
        self.params = incoming_params;
        self.operator.set_scales(&incoming_scales).expect("incoming scales remain admitted");
        result
    }
}

#[cfg(test)]
#[path = "binary_tests.rs"]
mod tests;
