//! Two-grid error in the actual squared response-fitting objective, including
//! imposed motion. The secant weight w*(r_f+r_c-2*t)/(2*scale^2) makes the
//! response difference identity EXACT for quadratic misfits. A zero coarse
//! misfit therefore does not suppress a nonzero enriched discrepancy.
//! This is a numerical two-grid comparison, not a continuum error bound.
use super::*;
use std::collections::BTreeMap;
use fs_cutfem::elastic3::adaptive::dirichlet::PrescribedMotion3;
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::elastic3::surface::{ReferenceLoad3, SurfaceForce3};
use fs_cutfem::octree3::Octant3;
use fs_dwr::elasticity3::{estimate_motion_goal3, dorfler3, GoalError3, GoalEstimate3, GoalFields3, GoalMarking3, GoalOptions3};
use crate::sdf3_goal::GoalRefinementError3;

/// The same linear observation law on either grid, never a transferred nodal q.
#[derive(Clone, Copy)]
pub struct ReferenceResponseTarget3<'a> {
    pub observation: ReferenceLoad3<'a>,
    pub target: f64,
    pub scale: f64,
    pub weight: f64,
}
/// Reference laws corresponding to one independently solved experiment.
#[derive(Clone, Copy)]
pub struct ReferenceResponseCase3<'a> {
    pub load: ReferenceLoad3<'a>,
    pub prescribed: Option<&'a PrescribedMotion3<'a>>,
    pub targets: &'a [ReferenceResponseTarget3<'a>],
}
#[derive(Debug, Clone, Copy)]
pub struct ResponseRefinementOptions3 {
    /// Same response weights and volume coefficient as the optimization model.
    pub response: ResponseOptions3,
    pub numerical: GoalOptions3,
    pub max_transfer_terms: usize,
}
impl Default for ResponseRefinementOptions3 {
    fn default() -> Self {
        Self { response: ResponseOptions3::default(), numerical: GoalOptions3::default(), max_transfer_terms: 2_000_000 }
    }
}
#[derive(Debug, Clone)]
pub struct ResponseCaseRefinement3 {
    pub coarse_responses: Vec<f64>,
    pub fine_responses: Vec<f64>,
    /// Coefficients of a common linear observation law, chosen AFTER both
    /// primal solves. These are NOT the optimizer's coarse-point derivatives.
    pub secant_weights: Vec<f64>,
    /// Linearized two-grid evidence whose correction equals this case's
    /// quadratic misfit difference. The goal VALUES are linear secant values,
    /// not the nonnegative squared misfits reported at family level.
    pub goal: GoalEstimate3,
    /// Enriched field at inherited physical stiffness, not a reoptimized state.
    pub fine_displacement: Vec<f64>,
}
#[derive(Debug, Clone)]
pub struct ResponseRefinement3 {
    pub coarse_objective: f64,
    pub fine_objective: f64,
    pub coarse_volume_fraction: f64,
    pub fine_volume_fraction: f64,
    /// Change of numerical measure for the same inherited projected density.
    pub volume_correction: f64,
    pub cases: Vec<ResponseCaseRefinement3>,
    pub marking_mass: BTreeMap<Octant3, f64>,
    pub identity_relative_defect: f64,
    pub work: SolveWork,
}
impl ResponseRefinement3 {
    pub fn correction(&self) -> f64 {
        self.volume_correction + self.cases.iter().map(|c|c.goal.correction()).sum::<f64>()
    }
    /// Marks the ORIGINAL grid from summed absolute case-local contributions.
    /// Volume quadrature discrepancy is reported separately, not fabricated as
    /// an adjoint indicator. Empty marks do not certify accuracy.
    pub fn mark(&self, theta: f64, max_marks: usize, checkpoint: impl FnMut()->ControlFlow<()>)
        -> Result<GoalMarking3, GoalError3> {
        dorfler3(&self.marking_mass,theta,max_marks,checkpoint)
    }
}
fn checked(x: f64) -> Result<f64, GoalRefinementError3> {
    if x.is_finite() { Ok(x) } else { Err(GoalRefinementError3::Invalid("response enrichment overflow")) }
}
fn matches_value(a: f64, b: f64, tolerance: f64) -> bool {
    if !a.is_finite() || !b.is_finite() { return false; }
    let scale=a.abs().max(b.abs());
    scale==0.0 || (a/scale-b/scale).abs()<=tolerance
}
fn weighted_body(targets: &[ReferenceResponseTarget3<'_>], weights: &[f64], p: [f64;3]) -> [f64;3] {
    let mut sum=[0.0;3];
    for (t,&w) in targets.iter().zip(weights) {
        if w==0.0 { continue; }
        if let Some(f)=t.observation.body { let v=f(p); for c in 0..3 {sum[c]+=w*v[c];} }
    }
    sum
}
fn weighted_surface(targets: &[ReferenceResponseTarget3<'_>], weights: &[f64], p: [f64;3], n: [f64;3]) -> [f64;3] {
    let mut sum=[0.0;3];
    for (t,&w) in targets.iter().zip(weights) {
        if w==0.0 {continue;}
        let v=match t.observation.surface {
            Some(SurfaceForce3::Traction(f))=>f(p,n),
            Some(SurfaceForce3::Pressure(f))=>{let pressure=f(p);n.map(|v|-pressure*v)},
            None=>continue,
        };
        for c in 0..3 {sum[c]+=w*v[c];}
    }
    sum
}

impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Compare an explicitly retained response evaluation with an enriched
    /// physical model. Both preparation policies are supplied by their existing
    /// backends (identity/Jacobi/two-level/multilevel); each is prepared once.
    /// All solves and setup consume the SAME caller control. No solver is added.
    ///
    /// Reconstruct the retained raw design and check the entire coarse family
    /// before preparing or solving enriched physics. Each observation, external
    /// force and prescribed motion must be the same pure physical law that made
    /// `accepted`. Cached aggregate adjoints are NOT reused: a secant observation
    /// needs its own coarse/fine adjoints, including at an exact coarse fit.
    ///
    /// The enriched operator is consumed and receives parent-inherited stiffness.
    /// No raw-density refiltering occurs on it. Incoming coarse scales are restored
    /// on success AND returned errors, so a failed estimate cannot replace an
    /// accepted optimizer state. Panicking callbacks are not caught. The caller
    /// still owns equality of the implicit geometry and selected support patches.
    pub fn estimate_response_enrichment<F: AdaptiveSdf3Elasticity>(&mut self, mut enriched: F,
        accepted: &ResponseEvaluation3, cases: &[ReferenceResponseCase3<'_>],
        options: ResponseRefinementOptions3, control: &mut SolveControl<'_>)
        -> Result<ResponseRefinement3, GoalRefinementError3> {
        control.checkpoint("response-goal-start")?;
        if accepted.reaction_responses.iter().any(|r| !r.is_empty()) {
            return Err(GoalRefinementError3::Invalid("reaction targets require a reaction-aware refinement goal"));
        }
        let mut count=0usize;
        if cases.is_empty() || cases.len()>options.response.max_cases
            || cases.len()!=accepted.displacements.len() || cases.len()!=accepted.responses.len()
            || !options.response.volume_weight.is_finite() || options.response.volume_weight<0.0
            || ![options.numerical.residual_tolerance,options.numerical.identity_tolerance]
                .iter().all(|v|v.is_finite() && *v>0.0 && *v<1.0)
            || accepted.rho.len()!=self.cells() || accepted.rho.iter().any(|r|!r.is_finite() || !(0.0..=1.0).contains(r)) {
            return Err(GoalRefinementError3::Invalid("invalid response enrichment family/design/options"));
        }
        let mut positive=options.response.volume_weight>0.0;
        for (case,responses) in cases.iter().zip(&accepted.responses) {
            count=count.checked_add(case.targets.len()).filter(|n|*n<=options.response.max_observations)
                .ok_or(GoalRefinementError3::Invalid("response enrichment observation budget"))?;
            if case.targets.is_empty() || case.targets.len()!=responses.len() {
                return Err(GoalRefinementError3::Invalid("response enrichment target count"));
            }
            for t in case.targets {
                if !t.target.is_finite() || !t.scale.is_finite() || t.scale<=0.0 || !t.weight.is_finite() || t.weight<0.0 {
                    return Err(GoalRefinementError3::Invalid("invalid response enrichment target"));
                }
                positive |= t.weight>0.0;
            }
        }
        if !positive {return Err(GoalRefinementError3::Invalid("zero response enrichment objective"));}
        let design=self.design(&accepted.rho,control)?;
        if design.scales!=accepted.scales || design.projected!=accepted.projected_rho {
            return Err(GoalRefinementError3::Invalid("retained response design does not match current pipeline"));
        }
        let previous=self.operator.scales().to_vec();
        self.operator.set_scales(&design.scales)?;
        let result=(|| {
            let coarse=self.operator.adaptive();
            // Validate all accepted fields and target laws before numerical setup.
            for ((case,u),responses) in cases.iter().zip(&accepted.displacements).zip(&accepted.responses) {
                let rhs=coarse.reference_load_with_motion(case.load,case.prescribed,||poll(control))?;
                let residual=coarse.field_residual(u,&rhs,||poll(control))?;
                if residual>options.numerical.residual_tolerance {
                    return Err(GoalError3::FieldResidual{field:"coarse-response-primal",value:residual}.into());
                }
                for (t,&expected) in case.targets.iter().zip(responses) {
                    let q=coarse.reference_load(t.observation,||poll(control))?;
                    if !matches_value(dot(&q,u),expected,options.numerical.identity_tolerance) {
                        return Err(GoalRefinementError3::Invalid("coarse observation law/response mismatch"));
                    }
                }
            }
            let inherited=AdaptiveTransfer3::new(coarse,enriched.adaptive(),options.max_transfer_terms,||poll(control))?.inherited_scales();
            enriched.set_scales(&inherited)?;
            let transfer=AdaptiveTransfer3::new(coarse,enriched.adaptive(),options.max_transfer_terms,||poll(control))?;
            if !transfer.fine().leaves().iter().zip(transfer.parents())
                .any(|(f,&c)| f.level()>coarse.leaves()[c].level()) {
                return Err(GoalRefinementError3::Invalid("response goal requires active-space enrichment"));
            }
            let fine_volumes=enriched.volumes();let total=checked(fine_volumes.iter().sum())?;
            if total<=0.0 {return Err(GoalRefinementError3::Invalid("empty enriched volume"));}
            let fine_volume=checked(fine_volumes.iter().zip(transfer.parents())
                .map(|(v,&parent)|(v/total)*design.projected[parent]).sum())?;
            let alpha=options.response.volume_weight;
            let mut report=ResponseRefinement3 {
                coarse_objective:checked(alpha*design.volume)?,fine_objective:checked(alpha*fine_volume)?,
                coarse_volume_fraction:design.volume,fine_volume_fraction:fine_volume,
                volume_correction:checked(alpha*(fine_volume-design.volume))?,cases:Vec::with_capacity(cases.len()),
                marking_mass:BTreeMap::new(),identity_relative_defect:0.0,work:control.work(),
            };
            let coarse_prepared=self.operator.prepare_elasticity(control)?;
            let fine_prepared=enriched.prepare_elasticity(control)?;
            for (case,uc) in cases.iter().zip(&accepted.displacements) {
                control.checkpoint("response-goal-case")?;
                let rhs=transfer.fine().reference_load_with_motion(case.load,case.prescribed,||poll(control))?;
                let uf=checked_solve_preconditioned(&enriched,&fine_prepared,&rhs,1e-11,"response-goal-primal",control)?;
                let mut vc=Vec::with_capacity(case.targets.len());let mut vf=Vec::with_capacity(case.targets.len());
                let mut weights=Vec::with_capacity(case.targets.len());
                for target in case.targets {
                    let qc=coarse.reference_load(target.observation,||poll(control))?;
                    let qf=transfer.fine().reference_load(target.observation,||poll(control))?;
                    let c=checked(dot(&qc,uc))?;let f=checked(dot(&qf,&uf))?;
                    vc.push(c);vf.push(f);
                    if target.weight==0.0 {weights.push(0.0);continue;}
                    let ec=checked((c-target.target)/target.scale)?;let ef=checked((f-target.target)/target.scale)?;
                    report.coarse_objective=checked(report.coarse_objective+0.5*target.weight*ec*ec)?;
                    report.fine_objective=checked(report.fine_objective+0.5*target.weight*ef*ef)?;
                    weights.push(checked(target.weight*((0.5*ec+0.5*ef)/target.scale))?);
                }
                let body=|p|weighted_body(case.targets,&weights,p);
                let surface=|p,n|weighted_surface(case.targets,&weights,p,n);
                let goal=ReferenceLoad3 {
                    body:if case.targets.iter().any(|t|t.observation.body.is_some()) {Some(&body)} else {None},
                    surface:if case.targets.iter().any(|t|t.observation.surface.is_some()) {Some(SurfaceForce3::Traction(&surface))} else {None},
                };
                let qc=coarse.reference_load(goal,||poll(control))?;
                let qf=transfer.fine().reference_load(goal,||poll(control))?;
                let zc=checked_solve_preconditioned(&self.operator,&coarse_prepared,&qc,1e-11,"response-goal-coarse-adjoint",control)?;
                let zf=checked_solve_preconditioned(&enriched,&fine_prepared,&qf,1e-11,"response-goal-fine-adjoint",control)?;
                let evidence=estimate_motion_goal3(&transfer,case.load,case.prescribed,goal,
                    GoalFields3{coarse_primal:uc,fine_primal:&uf,coarse_adjoint:&zc,fine_adjoint:&zf},
                    options.numerical,||poll(control))?;
                for (&cell,local) in &evidence.cells {
                    let value=report.marking_mass.entry(cell).or_insert(0.0);
                    *value=checked(*value+local.marking_mass)?;
                }
                report.cases.push(ResponseCaseRefinement3{coarse_responses:vc,fine_responses:vf,
                    secant_weights:weights,goal:evidence,fine_displacement:uf});
            }
            let correction=checked(report.correction())?;
            let scale=report.coarse_objective.abs().max(report.fine_objective.abs()).max(correction.abs());
            if scale>0.0 {report.identity_relative_defect=((report.fine_objective/scale-report.coarse_objective/scale)-correction/scale).abs();}
            if !report.identity_relative_defect.is_finite() || report.identity_relative_defect>options.numerical.identity_tolerance {
                return Err(GoalError3::Identity{relative_defect:report.identity_relative_defect}.into());
            }
            control.checkpoint("response-goal-publish")?;report.work=control.work();Ok(report)
        })();
        self.operator.set_scales(&previous).expect("incoming scales were admitted");
        result
    }
}
