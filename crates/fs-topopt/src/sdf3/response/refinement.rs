//! Two-grid error in squared displacement AND embedded-reaction response loss.
//! Secant weights retain a nonzero enriched discrepancy at an exact coarse fit.
//! Reaction goals are affine in displacement; their boundary-offset transfer
//! is kept explicitly, not silently discarded or relabeled as residual error.
//! This is a numerical two-grid comparison, not a continuum error bound.
use super::*;
use std::collections::BTreeMap;
use fs_cutfem::elastic3::adaptive::dirichlet::PrescribedMotion3;
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::elastic3::surface::{ReferenceLoad3, SurfaceForce3};
use fs_cutfem::octree3::Octant3;
use fs_dwr::elasticity3::{estimate_motion_goal3, dorfler3, GoalError3, GoalEstimate3, GoalFields3, GoalMarking3, GoalOptions3};
use fs_dwr::elasticity3::reaction::{AffineMotionGoal3, assemble_affine_motion_goal3, estimate_affine_motion_goal3};
use fs_dwr::elasticity3::response::EquilibriumLoad3;
use crate::sdf3_goal::GoalRefinementError3;

/// The same observation law on either grid, never a transferred nodal q.
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
    /// Chosen after both primal solves, not coarse fitting derivatives.
    pub secant_weights: Vec<f64>,
    pub coarse_reactions: Vec<f64>,
    pub fine_reactions: Vec<f64>,
    pub reaction_secant_weights: Vec<f64>,
    /// Derivative-work decomposition; its goal values are NOT squared loss.
    pub goal: GoalEstimate3,
    /// Fine minus coarse constant term of the secant-weighted reaction goal.
    /// Required in the objective identity; NOT an adjoint cell indicator.
    pub reaction_offset_correction: f64,
    pub fine_displacement: Vec<f64>,
}
#[derive(Debug, Clone)]
pub struct ResponseRefinement3 {
    pub coarse_objective: f64,
    pub fine_objective: f64,
    pub coarse_volume_fraction: f64,
    pub fine_volume_fraction: f64,
    /// Numerical measure change for inherited projected density, not refiltering.
    pub volume_correction: f64,
    pub cases: Vec<ResponseCaseRefinement3>,
    pub marking_mass: BTreeMap<Octant3, f64>,
    /// Scaled by objective and, for affine reactions, the actual cancelling
    /// identity terms. This is not relative accuracy of a near-zero loss.
    pub identity_relative_defect: f64,
    pub work: SolveWork,
}
impl ResponseRefinement3 {
    pub fn correction(&self) -> f64 {
        self.volume_correction + self.cases.iter().map(|c|c.goal.correction()+c.reaction_offset_correction).sum::<f64>()
    }
    /// Mark ORIGINAL-grid cells using summed absolute local residuals. Neither
    /// measure nor reaction-offset transfer is fabricated into an indicator.
    /// Empty marks do not certify accuracy, including when transfer is nonzero.
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
fn weighted_reaction(targets: &[ReactionTarget3<'_>], weights: &[f64], p: [f64;3], n: [f64;3]) -> [f64;3] {
    let mut sum=[0.0;3];
    for (t,&w) in targets.iter().zip(weights) {
        if w==0.0 {continue;}
        let h=(t.mode)(p,n);for c in 0..3 {sum[c]+=w*h[c];}
    }
    sum
}

impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Original displacement-only API. Reaction-bearing evidence is refused,
    /// never silently stripped. Both entry points share the same implementation.
    pub fn estimate_response_enrichment<F: AdaptiveSdf3Elasticity>(&mut self, enriched: F,
        accepted: &ResponseEvaluation3, cases: &[ReferenceResponseCase3<'_>],
        options: ResponseRefinementOptions3, control: &mut SolveControl<'_>)
        -> Result<ResponseRefinement3, GoalRefinementError3> {
        if accepted.reaction_responses.iter().any(|r|!r.is_empty()) {
            return Err(GoalRefinementError3::Invalid("reaction-bearing objectives require reaction-aware refinement"));
        }
        self.estimate_response_family(enriched,accepted,cases,None,options,control)
    }

    /// Estimate the COMPLETE mixed fitting objective, including reaction-only
    /// cases. Supply one reaction-target slice per experiment in its original
    /// order. Modes, external laws, imposed motion, targets, scales and weights
    /// must match the retained evaluation. Reactions are reintegrated, never
    /// transferred nodally or differentiated through the iterative solver.
    ///
    /// Each backend prepares once for its family. Fine stiffness is inherited,
    /// not obtained by refiltering raw density. The entire coarse family is
    /// revalidated before preparation. All work uses the shared control; incoming
    /// coarse scales are restored on success AND returned errors. No continuum,
    /// actuator-work, follower-load, or shape-derivative claim.
    #[allow(clippy::too_many_arguments)]
    pub fn estimate_response_enrichment_with_reactions<F: AdaptiveSdf3Elasticity>(&mut self, enriched: F,
        accepted: &ResponseEvaluation3, cases: &[ReferenceResponseCase3<'_>], reactions: &[&[ReactionTarget3<'_>]],
        options: ResponseRefinementOptions3, control: &mut SolveControl<'_>)
        -> Result<ResponseRefinement3, GoalRefinementError3> {
        self.estimate_response_family(enriched,accepted,cases,Some(reactions),options,control)
    }

    #[allow(clippy::too_many_arguments)]
    fn estimate_response_family<F: AdaptiveSdf3Elasticity>(&mut self, mut enriched: F,
        accepted: &ResponseEvaluation3, cases: &[ReferenceResponseCase3<'_>], reactions: Option<&[&[ReactionTarget3<'_>]]>,
        options: ResponseRefinementOptions3, control: &mut SolveControl<'_>)
        -> Result<ResponseRefinement3, GoalRefinementError3> {
        control.checkpoint("response-goal-start")?;
        let mut count=0usize;
        if cases.is_empty() || cases.len()>options.response.max_cases
            || cases.len()!=accepted.displacements.len() || cases.len()!=accepted.responses.len()
            || cases.len()!=accepted.reaction_responses.len() || reactions.is_some_and(|r|r.len()!=cases.len())
            || !options.response.volume_weight.is_finite() || options.response.volume_weight<0.0
            || ![options.numerical.residual_tolerance,options.numerical.identity_tolerance]
                .iter().all(|v|v.is_finite() && *v>0.0 && *v<1.0)
            || accepted.rho.len()!=self.cells() || accepted.rho.iter().any(|r|!r.is_finite() || !(0.0..=1.0).contains(r)) {
            return Err(GoalRefinementError3::Invalid("invalid response enrichment family/design/options"));
        }
        let mut positive=options.response.volume_weight>0.0;
        for (index,case) in cases.iter().enumerate() {
            let rt=reactions.map_or(&[][..],|r|r[index]);
            count=count.checked_add(case.targets.len()).and_then(|n|n.checked_add(rt.len()))
                .filter(|n|*n<=options.response.max_observations)
                .ok_or(GoalRefinementError3::Invalid("response enrichment observation budget"))?;
            if (case.targets.is_empty()&&rt.is_empty()) || case.targets.len()!=accepted.responses[index].len()
                || rt.len()!=accepted.reaction_responses[index].len() {
                return Err(GoalRefinementError3::Invalid("response enrichment target count"));
            }
            for (target,scale,weight) in case.targets.iter().map(|t|(t.target,t.scale,t.weight))
                .chain(rt.iter().map(|t|(t.target,t.scale,t.weight))) {
                if !target.is_finite() || !scale.is_finite() || scale<=0.0 || !weight.is_finite() || weight<0.0 {
                    return Err(GoalRefinementError3::Invalid("invalid response enrichment target"));
                }
                positive |= weight>0.0;
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
            // Validate the COMPLETE family before preparation or fine solves.
            let mut retained_loss=0.0;
            for (index,case) in cases.iter().enumerate() {
                let u=&accepted.displacements[index];
                let rhs=coarse.reference_load_with_motion(case.load,case.prescribed,||poll(control))?;
                let residual=coarse.field_residual(u,&rhs,||poll(control))?;
                if residual>options.numerical.residual_tolerance {
                    return Err(GoalError3::FieldResidual{field:"coarse-response-primal",value:residual}.into());
                }
                for (t,&expected) in case.targets.iter().zip(&accepted.responses[index]) {
                    let q=coarse.reference_load(t.observation,||poll(control))?;
                    if !matches_value(dot(&q,u),expected,options.numerical.identity_tolerance) {
                        return Err(GoalRefinementError3::Invalid("coarse observation law/response mismatch"));
                    }
                    if t.weight>0.0 {let e=checked((expected-t.target)/t.scale)?;retained_loss=checked(retained_loss+0.5*t.weight*e*e)?;}
                }
                for (t,&expected) in reactions.map_or(&[][..],|r|r[index]).iter().zip(&accepted.reaction_responses[index]) {
                    let actual=coarse.embedded_reaction(u,case.prescribed,t.mode,||poll(control))?.value;
                    if !matches_value(actual,expected,options.numerical.identity_tolerance) {
                        return Err(GoalRefinementError3::Invalid("coarse reaction law/response mismatch"));
                    }
                    if t.weight>0.0 {let e=checked((expected-t.target)/t.scale)?;retained_loss=checked(retained_loss+0.5*t.weight*e*e)?;}
                }
            }
            // Mixed goals must not accept target/weight changes under cached loss.
            if reactions.is_some() && !matches_value(checked(retained_loss+options.response.volume_weight*design.volume)?,
                accepted.objective,options.numerical.identity_tolerance) {
                return Err(GoalRefinementError3::Invalid("retained objective differs from supplied response targets"));
            }
            let inherited=AdaptiveTransfer3::new(coarse,enriched.adaptive(),options.max_transfer_terms,||poll(control))?.inherited_scales();
            enriched.set_scales(&inherited)?;
            let transfer=AdaptiveTransfer3::new(coarse,enriched.adaptive(),options.max_transfer_terms,||poll(control))?;
            if !transfer.fine().leaves().iter().zip(transfer.parents()).any(|(f,&c)| f.level()>coarse.leaves()[c].level()) {
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
            for (index,case) in cases.iter().enumerate() {
                control.checkpoint("response-goal-case")?;
                let uc=&accepted.displacements[index];let rt=reactions.map_or(&[][..],|r|r[index]);
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
                let mut rc=Vec::with_capacity(rt.len());let mut rf=Vec::with_capacity(rt.len());let mut rw=Vec::with_capacity(rt.len());
                for target in rt {
                    let c=coarse.embedded_reaction(uc,case.prescribed,target.mode,||poll(control))?.value;
                    let f=transfer.fine().embedded_reaction(&uf,case.prescribed,target.mode,||poll(control))?.value;
                    rc.push(c);rf.push(f);
                    if target.weight==0.0 {rw.push(0.0);continue;}
                    let ec=checked((c-target.target)/target.scale)?;let ef=checked((f-target.target)/target.scale)?;
                    report.coarse_objective=checked(report.coarse_objective+0.5*target.weight*ec*ec)?;
                    report.fine_objective=checked(report.fine_objective+0.5*target.weight*ef*ef)?;
                    rw.push(checked(target.weight*((0.5*ec+0.5*ef)/target.scale))?);
                }
                let body=|p|weighted_body(case.targets,&weights,p);
                let surface=|p,n|weighted_surface(case.targets,&weights,p,n);
                let goal=ReferenceLoad3 {
                    body:if case.targets.iter().any(|t|t.observation.body.is_some()) {Some(&body)} else {None},
                    surface:if case.targets.iter().any(|t|t.observation.surface.is_some()) {Some(SurfaceForce3::Traction(&surface))} else {None},
                };
                let mode=|p,n|weighted_reaction(rt,&rw,p,n);
                let affine=AffineMotionGoal3{displacement:goal,reaction_mode:Some(&mode)};
                let (qc,qf)=if rt.is_empty() {
                    (coarse.reference_load(goal,||poll(control))?,transfer.fine().reference_load(goal,||poll(control))?)
                } else {
                    (assemble_affine_motion_goal3(coarse,affine,case.prescribed,||poll(control))?.gradient,
                     assemble_affine_motion_goal3(transfer.fine(),affine,case.prescribed,||poll(control))?.gradient)
                };
                let zc=checked_solve_preconditioned(&self.operator,&coarse_prepared,&qc,1e-11,"response-goal-coarse-adjoint",control)?;
                let zf=checked_solve_preconditioned(&enriched,&fine_prepared,&qf,1e-11,"response-goal-fine-adjoint",control)?;
                let fields=GoalFields3{coarse_primal:uc,fine_primal:&uf,coarse_adjoint:&zc,fine_adjoint:&zf};
                let (evidence,offset)=if rt.is_empty() {
                    (estimate_motion_goal3(&transfer,case.load,case.prescribed,goal,fields,options.numerical,||poll(control))?,0.0)
                } else {
                    let e=estimate_affine_motion_goal3(&transfer,EquilibriumLoad3{external:case.load,prescribed:case.prescribed},
                        affine,fields,options.numerical,||poll(control))?;
                    (e.linearized,e.offset_transfer)
                };
                for (&cell,local) in &evidence.cells {
                    let value=report.marking_mass.entry(cell).or_insert(0.0);*value=checked(*value+local.marking_mass)?;
                }
                report.cases.push(ResponseCaseRefinement3{coarse_responses:vc,fine_responses:vf,secant_weights:weights,
                    coarse_reactions:rc,fine_reactions:rf,reaction_secant_weights:rw,
                    goal:evidence,reaction_offset_correction:offset,fine_displacement:uf});
            }
            let correction=checked(report.correction())?;
            let mut scale=report.coarse_objective.abs().max(report.fine_objective.abs()).max(correction.abs());
            for case in &report.cases {
                if !case.coarse_reactions.is_empty() {
                    for term in [case.goal.dwr,case.goal.coarse_space,case.goal.goal_transfer,case.goal.algebraic,case.reaction_offset_correction] {
                        scale=scale.max(term.abs());
                    }
                }
            }
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
