//! Check an entire simultaneous-load requirement on an enriched physical grid.
//! Both grids solve the same authored load laws and inherited eroded material.
pub mod binary;
use std::ops::ControlFlow;
use fs_cutfem::elastic3::{adaptive::{AdaptiveElasticity3, enrichment::AdaptiveTransfer3}, surface::ReferenceLoad3};
use super::LoadSpaceCap3;
use super::super::super::{LoadSpaceSpectrum3, spectrum, SPECTRAL_TOLERANCE};
use crate::pipeline::LoadCase;
use crate::sdf3::{AdaptiveSdf3Elasticity, CutDensityStudy3, CutDensityEvaluation3, Sdf3Elasticity, failure};
use crate::sdf3::robust::RobustOptions3;
use crate::sdf3_goal::{ComplianceRefinement3, GoalReferenceLoad3, GoalRefinementError3, GoalRefinementOptions3};
use crate::{EvaluationStop, SimpParams, SolveControl, SolveWork};

/// Worst combination on ONE grid. Force and displacement use that grid's node
/// order. A repeated leading eigenvalue need not have a unique maximizing q.
#[derive(Debug, Clone)]
pub struct LoadSpaceGridCheck3 {
    pub spectrum: LoadSpaceSpectrum3,
    pub worst_coefficients: Vec<f64>,
    pub worst_force: Vec<f64>,
    pub worst_displacement: Vec<f64>,
    pub worst_compliance: f64,
    pub witness_relative_defect: f64,
}

#[derive(Debug, Clone)]
pub struct LoadSpaceCapRefinement3 {
    pub cap: LoadSpaceCap3,
    /// Nominal model and raw controls that generated the coarse eroded field.
    pub params: SimpParams,
    pub robustness: RobustOptions3,
    pub rho: Vec<f64>,
    pub eroded: CutDensityEvaluation3,
    pub coarse: LoadSpaceGridCheck3,
    pub fine: LoadSpaceGridCheck3,
    /// Frobenius norm of G_fine-G_coarse, with BOTH Grams normalized by cap.
    /// In exact arithmetic this bounds |q^T (G_fine-G_coarse) q| for ||q||<=1.
    /// It detects changing adverse directions even if the maxima are equal.
    /// This numerical observation is NOT an interval or continuum error bound.
    pub normalized_gram_change: f64,
    /// Existing basis-wise DWR decompositions. Their uniform-mean aggregate
    /// and marking masses are diagnostics, NOT a full-load-set certificate.
    pub basis_assessment: ComplianceRefinement3,
    pub work: SolveWork,
}
impl LoadSpaceCapRefinement3 {
    /// Both freshly computed load-set maxima must meet the physical cap; the
    /// smooth tolerance never waives physical excess. The independent change
    /// tolerance applies to the WHOLE Gram difference, not only its diagonal
    /// or the difference between two maxima. No dilated-volume error is tested.
    #[must_use]
    pub fn stable_and_feasible(&self, smooth_tolerance: f64, change_tolerance: f64) -> bool {
        pair_admits(&self.cap, &self.coarse, &self.fine, self.normalized_gram_change,
            smooth_tolerance, change_tolerance)
    }
}

fn pair_admits(cap: &LoadSpaceCap3, coarse: &LoadSpaceGridCheck3, fine: &LoadSpaceGridCheck3,
    change: f64, smooth_tolerance: f64, change_tolerance: f64) -> bool {
    change_tolerance.is_finite() && change_tolerance >= 0.0
        && change.is_finite() && change >= 0.0 && change <= change_tolerance
        && [coarse, fine].iter().all(|grid| cap.admits_spectrum(
            cap.policy(), &grid.spectrum, grid.worst_compliance, smooth_tolerance))
}

fn poll(control: &mut SolveControl<'_>, stage: &'static str) -> ControlFlow<()> {
    if control.checkpoint(stage).is_ok() { ControlFlow::Continue(()) } else { ControlFlow::Break(()) }
}
fn invalid(message: &'static str) -> GoalRefinementError3 { GoalRefinementError3::Invalid(message) }

// Finite, scale-safe numerical Frobenius norm; do not discard off-diagonals.
fn gram_change(a: &[f64], b: &[f64]) -> Result<f64, EvaluationStop> {
    if a.len() != b.len() || a.is_empty() { return Err(failure("sdf3-load-space-grid-shape")); }
    let delta: Vec<_> = a.iter().zip(b).map(|(a,b)| a-b).collect();
    if delta.iter().any(|x| !x.is_finite()) { return Err(failure("sdf3-load-space-grid-change")); }
    let scale = delta.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
    let value = if scale == 0.0 { 0.0 } else {
        scale * delta.iter().map(|x| (x/scale)*(x/scale)).sum::<f64>().sqrt()
    };
    if !value.is_finite() { return Err(failure("sdf3-load-space-grid-change")); }
    Ok(value)
}

fn grid_check(forces: &[Vec<f64>], states: &[Vec<f64>], fixed: &[bool],
    cap: &LoadSpaceCap3, control: &mut SolveControl<'_>) -> Result<LoadSpaceGridCheck3, EvaluationStop> {
    let n = cap.policy().dimension();
    let dofs = fixed.len().checked_mul(3).ok_or_else(|| failure("sdf3-load-space-grid-shape"))?;
    if forces.len() != n || states.len() != n || forces.iter().chain(states).any(|v|
        v.len() != dofs || v.iter().any(|x| !x.is_finite())) {
        return Err(failure("sdf3-load-space-grid-shape"));
    }
    let mut gram = vec![0.0; n*n];
    for i in 0..n { for j in 0..n {
        let mut work = 0.0;
        for d in 0..dofs {
            if d % 256 == 0 { control.checkpoint("sdf3-load-space-grid-cross-work")?; }
            if !fixed[d/3] { work += forces[i][d]*states[j][d]; }
        }
        gram[i*n+j] = work/cap.cap();
    }}
    control.checkpoint("sdf3-load-space-grid-spectrum")?;
    let spectrum = spectrum(gram, cap.policy())?;
    control.checkpoint("sdf3-load-space-grid-spectrum")?;
    let worst_coefficients: Vec<_> = (0..n).map(|i| spectrum.eigenvectors[i*n+n-1]).collect();
    let mut worst_force = vec![0.0; dofs];
    let mut worst_displacement = vec![0.0; dofs];
    let mut worst_compliance = 0.0;
    for d in 0..dofs {
        if d % 256 == 0 { control.checkpoint("sdf3-load-space-grid-witness")?; }
        worst_force[d] = (0..n).map(|i| worst_coefficients[i]*forces[i][d]).sum();
        worst_displacement[d] = (0..n).map(|i| worst_coefficients[i]*states[i][d]).sum();
        if !fixed[d/3] { worst_compliance += worst_force[d]*worst_displacement[d]; }
    }
    let normalized_work = worst_compliance/cap.cap();
    let scale = normalized_work.abs().max(spectrum.maximum.abs());
    let witness_relative_defect = if scale == 0.0 { 0.0 }
        else { (normalized_work/scale-spectrum.maximum/scale).abs() };
    if !worst_compliance.is_finite() || worst_compliance < 0.0 || !normalized_work.is_finite()
        || worst_force.iter().chain(&worst_displacement).any(|v| !v.is_finite())
        || !witness_relative_defect.is_finite() || witness_relative_defect > SPECTRAL_TOLERANCE {
        return Err(failure("sdf3-load-space-grid-witness"));
    }
    Ok(LoadSpaceGridCheck3 { spectrum, worst_coefficients, worst_force, worst_displacement,
        worst_compliance, witness_relative_defect })
}

#[allow(clippy::too_many_arguments)]
fn admit_enrichment<O: AdaptiveSdf3Elasticity>(study: &CutDensityStudy3<O>, enriched: &AdaptiveElasticity3,
    laws: &[ReferenceLoad3<'_>], rho: &[f64], params: SimpParams, robustness: RobustOptions3,
    cap: &LoadSpaceCap3, options: GoalRefinementOptions3, control: &mut SolveControl<'_>)
    -> Result<(), GoalRefinementError3> {
    study.require_homogeneous_load_space().map_err(invalid)?;
    robustness.assert_valid(params);
    if laws.len() != cap.policy().dimension() || laws.len() > options.max_load_cases
        || options.max_transfer_terms == 0
        || ![options.numerical.residual_tolerance, options.numerical.identity_tolerance]
            .iter().all(|x| x.is_finite() && *x > 0.0 && *x < 1.0)
        || rho.len() != study.cells() || rho.iter().any(|r| !r.is_finite() || !(0.0..=1.0).contains(r)) {
        return Err(invalid("invalid load-space refinement model, basis or resource policy"));
    }
    control.checkpoint("sdf3-load-space-grid-start")?;
    let transfer = AdaptiveTransfer3::new(study.operator().adaptive(), enriched,
        options.max_transfer_terms, || poll(control, "sdf3-load-space-grid-transfer"))?;
    if !enriched.leaves().iter().zip(transfer.parents()).any(|(leaf, &parent)|
        leaf.level() > study.operator().adaptive().leaves()[parent].level()) {
        return Err(invalid("load-space assessment requires actual active-space enrichment"));
    }
    drop(transfer);
    Ok(())
}

impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    #[allow(clippy::too_many_arguments)]
    fn assess_load_space_grids(&self, enriched: AdaptiveElasticity3, laws: &[ReferenceLoad3<'_>],
        forces: &[Vec<f64>], states: &[Vec<f64>], cap: &LoadSpaceCap3,
        options: GoalRefinementOptions3, control: &mut SolveControl<'_>)
        -> Result<(LoadSpaceGridCheck3, LoadSpaceGridCheck3, f64, ComplianceRefinement3), GoalRefinementError3> {
        let coarse = grid_check(forces, states, self.operator().fixed(), cap, control)?;
        let loads: Vec<_> = laws.iter().map(|&load| GoalReferenceLoad3 {
            load, weight: 1.0/laws.len() as f64 }).collect();
        let fine_fixed = enriched.fixed().to_vec();
        let (basis_assessment, fields) = self.estimate_reference_compliance_with_fields(
            enriched, &loads, states, options, control)?;
        let fine = grid_check(&fields.forces, &fields.displacements, &fine_fixed, cap, control)?;
        let change = gram_change(&fine.spectrum.gram, &coarse.spectrum.gram)?;
        Ok((coarse, fine, change, basis_assessment))
    }

    /// Re-solve the authored basis on both grids and check the full combined-load
    /// cap. The finer grid inherits the COARSE ERODED PHYSICAL stiffness; it does
    /// not inherit/refilter raw controls or replace the current design.
    ///
    /// Reuse the existing reference-load DWR solver and its fixed preconditioner,
    /// retaining its actual fine fields to recover every cross compliance. No
    /// pairwise polarization solves or sampling of directions is required. Fine
    /// and coarse worst directions are recomputed independently. Basis DWR
    /// results stay diagnostic; no combined-load localization claim is made.
    ///
    /// Pure reference laws, SDF, bounds, base material, integration policy and
    /// homogeneous supports must describe the same physical problem. Transfer
    /// checks compatibility and actual active enrichment, not SDF identity.
    /// All solves share control; success and returned errors restore incoming
    /// nominal parameters/scales and publish no partial load-set assessment.
    #[allow(clippy::too_many_arguments)]
    pub fn estimate_load_space_cap_enrichment(&mut self, enriched: AdaptiveElasticity3,
        laws: &[ReferenceLoad3<'_>], rho: &[f64], params: SimpParams, robustness: RobustOptions3,
        cap: &LoadSpaceCap3, options: GoalRefinementOptions3, control: &mut SolveControl<'_>)
        -> Result<LoadSpaceCapRefinement3, GoalRefinementError3> {
        admit_enrichment(self, &enriched, laws, rho, params, robustness, cap, options, control)?;
        let incoming_params = self.params;
        let incoming_scales = self.operator().scales().to_vec();
        let result = (|| {
            let mut forces = Vec::with_capacity(laws.len());
            for &law in laws {
                forces.push(self.operator().adaptive().reference_load(law,
                    || poll(control, "sdf3-load-space-grid-load"))?);
            }
            let nodal: Vec<_> = forces.iter().map(|force| LoadCase {
                force, weight: 1.0/laws.len() as f64 }).collect();
            self.params = params;
            self.params.eta += robustness.threshold_offset;
            // Existing evaluator installs this same eroded design for the DWR
            // owner to inherit, including protected regions and actual cut mass.
            let eroded = self.evaluate(rho, &nodal, control)?;
            let (coarse, fine, normalized_gram_change, basis_assessment) =
                self.assess_load_space_grids(enriched, laws, &forces,
                    &eroded.objective.displacements, cap, options, control)?;
            control.checkpoint("sdf3-load-space-grid-publish")?;
            Ok(LoadSpaceCapRefinement3 { cap: cap.clone(), params, robustness, rho: rho.to_vec(),
                eroded, coarse, fine, normalized_gram_change, basis_assessment, work: control.work() })
        })();
        self.params = incoming_params;
        self.operator.set_scales(&incoming_scales).expect("incoming scales remain admitted");
        result
    }
}

#[cfg(test)]
#[path = "refinement_tests.rs"]
mod tests;
