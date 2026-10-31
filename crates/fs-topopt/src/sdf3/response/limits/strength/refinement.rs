//! Enriched-grid assessment of the retained joint response/stress design.
//! Reuse response-limit DWR primals; sample stress from those same fine fields.
use super::*;
use super::super::refinement::{ReferenceResponseLimitCase3, ResponseLimitRefinement3,
    ResponseLimitRefinementOptions3};
use crate::sdf3::Sdf3Elasticity;
use crate::sdf3::stress::{StressError3, StressMeasure3, StressOptions3, von_mises};
use crate::sdf3_goal::GoalRefinementError3;
use fs_dwr::elasticity3::GoalOptions3;
use std::ops::ControlFlow;

#[derive(Debug, Clone, Copy)]
pub struct ResponseStressRefinementOptions3 {
    pub joint: ResponseStressOptions3,
    pub numerical: GoalOptions3,
    pub max_transfer_terms: usize,
}
impl Default for ResponseStressRefinementOptions3 {
    fn default() -> Self {
        Self { joint: ResponseStressOptions3::default(), numerical: GoalOptions3::default(),
            max_transfer_terms: 2_000_000 }
    }
}
#[derive(Debug, Clone, PartialEq)]
pub enum ResponseStressRefinementError3 {
    Invalid(&'static str),
    Response(GoalRefinementError3),
    Stress(StressError3),
    Evaluation(EvaluationStop),
}
impl From<GoalRefinementError3> for ResponseStressRefinementError3 {
    fn from(e: GoalRefinementError3) -> Self { Self::Response(e) }
}
impl From<StressError3> for ResponseStressRefinementError3 {
    fn from(e: StressError3) -> Self {
        match e { StressError3::Evaluation(e) => Self::Evaluation(e), other => Self::Stress(other) }
    }
}
impl From<EvaluationStop> for ResponseStressRefinementError3 {
    fn from(e: EvaluationStop) -> Self { Self::Evaluation(e) }
}
impl std::fmt::Display for ResponseStressRefinementError3 {
    fn fmt(&self, f:&mut std::fmt::Formatter<'_>)->std::fmt::Result {
        write!(f,"joint response/stress refinement refused: {self:?}")
    }
}
impl std::error::Error for ResponseStressRefinementError3 {}

#[derive(Debug, Clone)]
pub struct StressGridSummary3 {
    pub aggregate: f64,
    pub sampled_relaxed_max: f64,
    pub sampled_physical_max: f64,
    pub point_count: usize,
}
#[derive(Debug, Clone)]
pub struct ResponseStressRefinement3 {
    pub response: ResponseLimitRefinement3,
    pub coarse_stress: StressGridSummary3,
    pub fine_stress: StressGridSummary3,
    pub coarse_constraint: f64,
    pub fine_constraint: f64,
    pub coarse_stress_residual: f64,
    pub fine_stress_residual: f64,
    pub coarse_worst_residual: f64,
    pub fine_worst_residual: f64,
    pub work: SolveWork,
}
impl ResponseStressRefinement3 {
    #[must_use]
    pub fn fine_feasible(&self, tolerance:f64)->bool {
        if !tolerance.is_finite() || tolerance<0.0 || !self.response.fine_feasible(tolerance) {
            return false;
        }
        let relaxed=self.fine_stress.sampled_relaxed_max/self.response_stress_limit()-1.0;
        let physical=self.fine_stress.sampled_physical_max/self.response_stress_limit()-1.0;
        [self.fine_constraint,self.fine_stress_residual,relaxed,physical,self.fine_worst_residual]
            .iter().all(|v|v.is_finite()&&*v<=tolerance)
    }
    fn response_stress_limit(&self)->f64 {
        self.fine_stress.aggregate/(1.0+self.fine_stress_residual)
    }
    #[must_use]
    pub fn max_normalized_change(&self)->f64 {
        let cap=self.response_stress_limit();
        self.response.max_normalized_change()
            .max((self.fine_stress_residual-self.coarse_stress_residual).abs())
            .max(((self.fine_stress.sampled_relaxed_max-self.coarse_stress.sampled_relaxed_max)/cap).abs())
            .max(((self.fine_stress.sampled_physical_max-self.coarse_stress.sampled_physical_max)/cap).abs())
    }
}

#[derive(Default)]
struct StressNorm { maximum:f64, sum:f64 }
impl StressNorm {
    fn add(&mut self,value:f64,p:f64)->Result<(),ResponseStressRefinementError3> {
        if !value.is_finite()||value<0.0 {return Err(ResponseStressRefinementError3::Invalid("invalid enriched stress sample"));}
        if value>self.maximum {
            self.sum=if self.maximum==0.0 {1.0}
                else {self.sum*fs_math::det::pow(self.maximum/value,p)+1.0};
            self.maximum=value;
        } else if value>0.0 {
            self.sum+=fs_math::det::pow(value/self.maximum,p);
        }
        if !self.sum.is_finite(){return Err(ResponseStressRefinementError3::Invalid("enriched stress norm overflow"));}
        Ok(())
    }
    fn value(&self,p:f64)->Result<f64,ResponseStressRefinementError3> {
        let v=if self.maximum==0.0 {0.0}else{self.maximum*fs_math::det::pow(self.sum,1.0/p)};
        if !v.is_finite()||v<self.maximum {return Err(ResponseStressRefinementError3::Invalid("invalid enriched stress aggregate"));}
        Ok(v)
    }
}
fn poll(control:&mut SolveControl<'_>)->ControlFlow<()> {
    if control.checkpoint("response-stress-refinement-physics").is_ok(){ControlFlow::Continue(())}
    else{ControlFlow::Break(())}
}
fn near(a:f64,b:f64,scale:f64,tol:f64)->bool {
    a.is_finite()&&b.is_finite()&&((a-b)/scale.max(1e-300)).abs()<=tol
}
fn stress_summary<'a,O,I>(operator:&O,projected:&[f64],fields:I,params:crate::SimpParams,
    options:StressOptions3,control:&mut SolveControl<'_>)
    ->Result<StressGridSummary3,ResponseStressRefinementError3>
where O:Sdf3Elasticity,I:IntoIterator<Item=&'a [f64]> {
    let fields:Vec<&[f64]>=fields.into_iter().collect();
    if projected.len()!=operator.cells() || projected.iter().any(|r|!r.is_finite()||!(0.0..=1.0).contains(r))
        || fields.is_empty() || fields.len()>options.max_cases || options.max_points==0
        || !options.relaxation_power.is_finite() || options.relaxation_power<1.0
        || options.relaxation_power>=params.penal || !options.aggregation_power.is_finite()
        || !(2.0..=64.0).contains(&options.aggregation_power)
        || fields.iter().any(|u|u.len()!=operator.n()||u.iter().any(|v|!v.is_finite())) {
        return Err(ResponseStressRefinementError3::Invalid("invalid enriched stress family or policy"));
    }
    let mut norm=StressNorm::default();let mut point_count=0usize;
    let mut relaxed_max=0.0f64;let mut physical_max=0.0f64;
    for u in fields {
        control.checkpoint("response-stress-refinement-stress")?;
        let points=operator.bulk_stress(u,options.max_points-point_count,||poll(control))
            .map_err(StressError3::from)?;
        point_count=point_count.checked_add(points.len())
            .filter(|n|*n<=options.max_points).ok_or(StressError3::PointBudget)?;
        for (i,point) in points.iter().enumerate() {
            if i%256==0 {control.checkpoint("response-stress-refinement-stress")?;}
            let vm=von_mises(&point.reference_stress).0;
            let r=projected[point.cell];
            let relaxed=if r==0.0 {0.0}else{fs_math::det::pow(r,options.relaxation_power)*vm};
            let physical=von_mises(&point.stress).0;
            norm.add(relaxed,options.aggregation_power)?;norm.add(physical,options.aggregation_power)?;
            relaxed_max=relaxed_max.max(relaxed);physical_max=physical_max.max(physical);
        }
    }
    if point_count==0 {return Err(ResponseStressRefinementError3::Invalid("empty enriched stress quadrature"));}
    Ok(StressGridSummary3{aggregate:norm.value(options.aggregation_power)?,
        sampled_relaxed_max:relaxed_max,sampled_physical_max:physical_max,point_count})
}

impl<O:AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Recheck every response inequality and sampled stress on an enriched grid
    /// at the retained physical material. Fine stress reuses the response DWR
    /// primal fields; it starts no additional equilibrium or adjoint solves.
    pub fn estimate_response_stress_enrichment<F:AdaptiveSdf3Elasticity>(&mut self,enriched:F,
        accepted:&ResponseStressEvaluation3,cases:&[ReferenceResponseLimitCase3<'_>],
        options:ResponseStressRefinementOptions3,control:&mut SolveControl<'_>)
        ->Result<ResponseStressRefinement3,ResponseStressRefinementError3> {
        control.checkpoint("response-stress-refinement-start")?;
        if !options.joint.stress_limit.is_finite()||options.joint.stress_limit<=0.0
            || accepted.stress_limit!=options.joint.stress_limit
            || accepted.stress.measure!=StressMeasure3::SampledPeakBound
            || accepted.stress.relaxation_power!=options.joint.stress.relaxation_power
            || accepted.stress.aggregation_power!=options.joint.stress.aggregation_power
            || accepted.response.rho!=accepted.stress.rho
            || accepted.response.projected_rho!=accepted.stress.projected_rho
            || accepted.response.displacements!=accepted.stress.displacements {
            return Err(ResponseStressRefinementError3::Invalid("retained joint evidence or stress policy mismatch"));
        }
        let response_options=ResponseLimitRefinementOptions3{limits:options.joint.response,
            numerical:options.numerical,max_transfer_terms:options.max_transfer_terms,
            ..Default::default()};
        let (response,enriched)=self.estimate_response_limit_enrichment_retaining(enriched,
            &accepted.response,cases,response_options,control)?;
        if response.fine_projected_rho.len()!=enriched.cells() {
            return Err(ResponseStressRefinementError3::Invalid("enriched projected-density map mismatch"));
        }
        let previous=self.operator.scales().to_vec();
        self.operator.set_scales(&accepted.stress.scales)
            .map_err(|e|ResponseStressError3::Response(ResponseError3::Physics(e)))
            .map_err(|_|ResponseStressRefinementError3::Invalid("failed to install retained stress material"))?;
        let coarse_result=stress_summary(&self.operator,&accepted.response.projected_rho,
            accepted.response.displacements.iter().map(Vec::as_slice),self.params,
            options.joint.stress,control);
        self.operator.set_scales(&previous).expect("incoming refinement material remains valid");
        let coarse_stress=coarse_result?;
        let tol=options.numerical.identity_tolerance;
        if !near(coarse_stress.aggregate,accepted.stress.aggregate,accepted.stress.aggregate.abs().max(1.0),tol)
            || !near(coarse_stress.sampled_relaxed_max,accepted.stress.sampled_relaxed_max,
                accepted.stress.aggregate.abs().max(1.0),tol)
            || !near(coarse_stress.sampled_physical_max,accepted.stress.sampled_physical_max,
                accepted.stress.aggregate.abs().max(1.0),tol) {
            return Err(ResponseStressRefinementError3::Invalid("retained coarse stress evidence changed"));
        }
        let fine_stress=stress_summary(&enriched,&response.fine_projected_rho,
            response.cases.iter().map(|c|c.fine_displacement.as_slice()),self.params,
            options.joint.stress,control)?;
        let cap=options.joint.stress_limit;
        let coarse_stress_residual=coarse_stress.aggregate/cap-1.0;
        let fine_stress_residual=fine_stress.aggregate/cap-1.0;
        let coarse_constraint=aggregate(&[response.coarse_constraint,coarse_stress_residual],
            options.joint.response.sharpness).map_err(ResponseStressError3::from)
            .map_err(|_|ResponseStressRefinementError3::Invalid("invalid coarse joint aggregate"))?.0;
        let fine_constraint=aggregate(&[response.fine_constraint,fine_stress_residual],
            options.joint.response.sharpness).map_err(ResponseStressError3::from)
            .map_err(|_|ResponseStressRefinementError3::Invalid("invalid fine joint aggregate"))?.0;
        let coarse_worst_residual=response.coarse_worst_residual
            .max(coarse_stress.sampled_relaxed_max/cap-1.0)
            .max(coarse_stress.sampled_physical_max/cap-1.0);
        let fine_worst_residual=response.fine_worst_residual
            .max(fine_stress.sampled_relaxed_max/cap-1.0)
            .max(fine_stress.sampled_physical_max/cap-1.0);
        if !near(coarse_constraint,accepted.constraint,1.0,tol)
            || !near(coarse_stress_residual,accepted.stress_residual,1.0,tol)
            || !near(coarse_worst_residual,accepted.worst_residual,1.0,tol) {
            return Err(ResponseStressRefinementError3::Invalid("retained coarse joint evidence changed"));
        }
        control.checkpoint("response-stress-refinement-publish")?;
        Ok(ResponseStressRefinement3{response,coarse_stress,fine_stress,coarse_constraint,
            fine_constraint,coarse_stress_residual,fine_stress_residual,coarse_worst_residual,
            fine_worst_residual,work:control.work()})
    }
}
