//! Filtered SIMP optimization on real 3-D raw-implicit CutFEM operators.
//!
//! Uniform Cartesian and Q1-conforming octree backgrounds share the SAME
//! filter/projection/SIMP chain, independent-load evaluator and OC driver.
//! Geometry and quadrature stay fixed inside each optimization call. Refined
//! studies reassemble geometry between calls, inherit only raw densities, and
//! must restore volume feasibility and re-solve their new baseline.
//! `continuation` advances SIMP/projection models on the SAME retained geometry,
//! restores stage feasibility and can gate each baseline with numerical gradients.
//! `adaptive_continuation` joins those stages with actual compliance-DWR
//! background refinement and atomic, freshly solved raw-design transfers.
//!
//! Numerical cut volumes, discrete sensitivities and energy-based marking are
//! not continuum certificates or proofs of optimality. The adaptive driver uses
//! the separate `sdf3_goal` residual estimator, not energy-based marking. The
//! graph Helmholtz filter is not the tetrahedral FEEC filter.

mod operator;
mod refine;
pub mod continuation;
pub mod adaptive_continuation;
pub mod response;
pub mod stress;
pub mod design;
pub use operator::{AdaptiveSdf3Elasticity, Sdf3Elasticity};
pub use refine::inherit_raw_densities3;
use fs_cutfem::elastic3::CutElasticity3;
use fs_solver::op::{CsrOp, LinearOp};
use fs_sparse::precond::{IdentityPrecond, Precond};
use crate::control::{EvaluationStop, SolveControl};
use crate::filter::{heaviside, heaviside_derivative};
use crate::multi_load::{MultiLoadOcIteration, MultiLoadOcOptions, MultiLoadOcReport, MultiLoadOcTermination};
use crate::pipeline::{LoadCase, MultiLoadCompliance, SimpParams};

fn failure(stage: &'static str) -> EvaluationStop { EvaluationStop::Breakdown { stage } }

fn assert_loads<O: Sdf3Elasticity>(study: &CutDensityStudy3<O>, loads: &[LoadCase<'_>]) {
    assert!(!loads.is_empty() && loads.iter().any(|l| l.weight > 0.0), "at least one positive load weight required");
    for load in loads {
        assert!(load.weight.is_finite() && load.weight >= 0.0, "invalid load weight");
        assert_eq!(load.force.len(), study.operator.n(), "load shape mismatch");
        assert!(load.force.iter().all(|f| f.is_finite()), "nonfinite load");
    }
}

fn assert_oc_inputs<O: Sdf3Elasticity>(study: &CutDensityStudy3<O>, loads: &[LoadCase<'_>],
    rho0: &[f64], options: MultiLoadOcOptions) {
    assert_loads(study, loads);
    assert_eq!(rho0.len(),study.cells(),"initial density shape mismatch");
    assert!(rho0.iter().all(|r|r.is_finite()&&(1e-3..=1.0).contains(r)),"OC densities must lie in [0.001,1]");
    assert!(options.volume_fraction.is_finite()&&options.volume_fraction>0.0&&options.volume_fraction<=1.0,"invalid volume cap");
    assert!(options.move_limit.is_finite()&&options.move_limit>0.0&&options.move_limit<=1.0,"invalid move limit");
    assert!(options.volume_tolerance.is_finite()&&options.volume_tolerance>=0.0&&options.volume_tolerance<options.volume_fraction,"invalid volume tolerance");
    assert!(options.change_tolerance.is_finite()&&options.change_tolerance>=0.0&&options.max_backtracks<=64,"invalid stopping/backtracking policy");
}

/// Geometry-bound density pipeline. The default backend preserves existing
/// Cartesian callers; an adaptive backend uses independent master-node DOFs.
pub struct CutDensityStudy3<O: Sdf3Elasticity = CutElasticity3> {
    operator: O,
    filter: CsrOp,
    mass: Vec<f64>,
    params: SimpParams,
    filter_radius: f64,
}

struct Design3 { projected: Vec<f64>, slope: Vec<f64>, scales: Vec<f64>, volume: f64 }

/// Complete equilibrium and exact discrete derivatives of one design.
#[derive(Debug, Clone)]
pub struct CutDensityEvaluation3 {
    /// Weighted independently solved compliance, per-load fields and gradient.
    pub objective: MultiLoadCompliance,
    /// Projected densities at exactly the solved design.
    pub projected_rho: Vec<f64>,
    /// Numerical material-volume fraction of the raw-SDF design domain.
    pub volume_fraction: f64,
    /// Derivative of projected volume through the same filter/projection map.
    pub volume_gradient: Vec<f64>,
}

impl<O: Sdf3Elasticity> CutDensityStudy3<O> {
    /// Bind the cut operator and build `A=M+r^2 L` on face-adjacent active cells.
    /// `M` contains normalized positive cut volumes. Each symmetric graph edge
    /// has weight harmonic_mean(M_i,M_j)/normal_center_separation^2. Each
    /// coarse/fine face patch contributes one edge, with no diagonal neighbors.
    /// Thus F=A^-1 M preserves constants and M-weighted volume in exact
    /// arithmetic, and its pullback is M A^-1, NOT A^-1 M.
    /// `radius` has the same length unit as the supplied operator coordinates.
    pub fn new(operator: O, radius: f64, params: SimpParams) -> Self {
        params.assert_valid();
        assert!(radius.is_finite() && radius >= 0.0 && (radius*radius).is_finite(), "invalid filter radius");
        let volumes = operator.volumes();
        let total: f64 = volumes.iter().sum();
        assert!(total.is_finite() && total > 0.0, "invalid cut-domain volume");
        let mass: Vec<f64> = volumes.iter().map(|v| v/total).collect();
        assert!(mass.iter().all(|v| v.is_finite() && *v > 0.0), "unrepresentable normalized cut volumes");
        let mut coo = fs_sparse::Coo::new(mass.len(), mass.len());
        for (i,&m) in mass.iter().enumerate() { coo.push(i,i,m); }
        for (i,j,distance) in operator.filter_edges() {
            let ratio = radius/distance;
            let low = mass[i].min(mass[j]); let high = mass[i].max(mass[j]);
            let weight = (2.0*low/(1.0+low/high))*ratio*ratio;
            assert!(distance.is_finite() && distance > 0.0 && weight.is_finite(), "unrepresentable graph filter edge");
            coo.push(i,i,weight); coo.push(j,j,weight);
            coo.push(i,j,-weight); coo.push(j,i,-weight);
        }
        Self { operator, filter: CsrOp::symmetric(coo.assemble()), mass, params, filter_radius: radius }
    }

    /// Read-only access to the geometry, accepted scales, and field ordering.
    #[must_use]
    pub const fn operator(&self) -> &O { &self.operator }
    /// Current material/projection model, matching the last accepted stage.
    #[must_use]
    pub const fn params(&self) -> SimpParams { self.params }
    /// Raw cell count, also the design dimension.
    #[must_use]
    pub fn cells(&self) -> usize { self.mass.len() }

    fn filter_apply(&self, rho: &[f64], control: &mut SolveControl<'_>) -> Result<Vec<f64>,EvaluationStop> {
        let rhs: Vec<f64> = rho.iter().zip(&self.mass).map(|(r,m)|r*m).collect();
        checked_solve(&self.filter,&rhs,1e-12,"sdf3-filter",control)
    }
    fn pullback(&self, local: &[f64], control: &mut SolveControl<'_>) -> Result<Vec<f64>,EvaluationStop> {
        let solved = checked_solve(&self.filter,local,1e-12,"sdf3-filter-transpose",control)?;
        Ok(solved.iter().zip(&self.mass).map(|(v,m)|v*m).collect())
    }
    fn design(&self,rho:&[f64],control:&mut SolveControl<'_>) -> Result<Design3,EvaluationStop> {
        control.checkpoint("sdf3-design")?;
        assert_eq!(rho.len(),self.cells(),"one raw density per active cut cell required");
        assert!(rho.iter().all(|r|r.is_finite()&&(0.0..=1.0).contains(r)),"invalid raw density");
        let filtered = self.filter_apply(rho,control)?;
        let p = self.params;
        let mut projected = Vec::with_capacity(rho.len());
        let mut slope = Vec::with_capacity(rho.len());
        let mut scales = Vec::with_capacity(rho.len());
        for &r in &filtered {
            let h = heaviside(r,p.beta,p.eta);
            if !h.is_finite() { return Err(failure("sdf3-projection")); }
            let physical = h.clamp(0.0,1.0);
            let derivative = if h < 0.0 || h > 1.0 { 0.0 } else { heaviside_derivative(r,p.beta,p.eta) };
            // Unlike a hidden positive rho floor, this value and derivative
            // are of the SAME map, including the exact zero-density limit.
            let power = if physical == 0.0 { 0.0 } else { fs_math::det::pow(physical,p.penal) };
            let scale = p.e_min+(1.0-p.e_min)*power;
            if !derivative.is_finite() || !scale.is_finite() { return Err(failure("sdf3-simp")); }
            projected.push(physical); slope.push(derivative); scales.push(scale);
        }
        let volume = projected.iter().zip(&self.mass).map(|(r,m)|r*m).sum();
        control.checkpoint("sdf3-design-publish")?;
        Ok(Design3 { projected, slope, scales, volume })
    }

    /// Volume and its exact discrete raw-density gradient, with no operator mutation.
    pub fn volume_and_gradient(&self,rho:&[f64],control:&mut SolveControl<'_>) -> Result<(f64,Vec<f64>),EvaluationStop> {
        let d = self.design(rho,control)?;
        let local: Vec<f64> = d.slope.iter().zip(&self.mass).map(|(s,m)|s*m).collect();
        Ok((d.volume,self.pullback(&local,control)?))
    }

    /// Evaluate real independent equilibria and the full filter/projection/SIMP
    /// chain. All filter, load and pullback solves share the existing SolveControl.
    /// Invalid inputs panic before mutation. Any returned stop restores scales;
    /// no partial load family or gradient is usable as an evaluation.
    /// An explicitly preconditioned backend prepares once per evaluated design,
    /// not once per load. Its setup work is retained even on a rejected trial.
    /// Bare backends preserve the original identity-preconditioned arithmetic.
    pub fn evaluate(&mut self,rho:&[f64],loads:&[LoadCase<'_>],control:&mut SolveControl<'_>)
        -> Result<CutDensityEvaluation3,EvaluationStop> {
        control.checkpoint("sdf3-evaluation")?;
        assert_loads(self, loads);
        let design = self.design(rho,control)?;
        let previous = self.operator.scales().to_vec();
        self.operator.set_scales(&design.scales).map_err(|_|failure("sdf3-scales"))?;
        let result = (|| {
            // Prepared once at this exact density, then shared by all loads.
            // Its immutable borrow ends before any rollback or next trial.
            let prepared = self.operator.prepare_elasticity(control)?;
            let mut compliance = 0.0;
            let mut local = vec![0.0;self.cells()];
            let mut displacements = Vec::with_capacity(loads.len());
            let mut case_compliances = Vec::with_capacity(loads.len());
            for load in loads {
                control.checkpoint("sdf3-elasticity")?;
                let rhs: Vec<f64> = load.force.iter().enumerate().map(|(i,&f)|if self.operator.fixed()[i/3] {0.0}else{f}).collect();
                let u = checked_solve_preconditioned(&self.operator,&prepared,&rhs,1e-11,"sdf3-elasticity",control)?;
                let c: f64 = rhs.iter().zip(&u).map(|(f,u)|f*u).sum();
                compliance += load.weight*c;
                let energies = self.operator.scale_quadratic_forms(&u).map_err(|_|failure("sdf3-energy"))?;
                for (value,energy) in local.iter_mut().zip(energies) { *value -= load.weight*energy; }
                case_compliances.push(c); displacements.push(u);
            }
            let p = self.params;
            for ((value,&r),&slope) in local.iter_mut().zip(&design.projected).zip(&design.slope) {
                let power = if p.penal == 1.0 {1.0} else if r == 0.0 {0.0} else {fs_math::det::pow(r,p.penal-1.0)};
                *value *= (1.0-p.e_min)*p.penal*power*slope;
            }
            if !compliance.is_finite() || !local.iter().all(|v|v.is_finite()) {return Err(failure("sdf3-compliance"));}
            let gradient = self.pullback(&local,control)?;
            let local_volume: Vec<f64> = self.mass.iter().zip(&design.slope).map(|(m,s)|m*s).collect();
            let volume_gradient = self.pullback(&local_volume,control)?;
            control.checkpoint("sdf3-evaluation-publish")?;
            Ok(CutDensityEvaluation3 {
                objective: MultiLoadCompliance { compliance, case_compliances, displacements, gradient },
                projected_rho: design.projected, volume_fraction: design.volume, volume_gradient,
            })
        })();
        if result.is_err() { self.operator.set_scales(&previous).expect("previous scales valid"); }
        result
    }
}

// Reuse the existing cumulative CG control, but do not admit a recurrence-only
// result as physics: verify b-Ax explicitly. Unlike native elastic3's correction
// solve, this adapter refuses a failed residual gate instead of restarting.
fn checked_solve(op:&impl LinearOp,rhs:&[f64],tol:f64,stage:&'static str,control:&mut SolveControl<'_>)
    -> Result<Vec<f64>,EvaluationStop> {
    checked_solve_preconditioned(op,&IdentityPrecond,rhs,tol,stage,control)
}
fn checked_solve_preconditioned(op:&impl LinearOp,preconditioner:&impl Precond,rhs:&[f64],tol:f64,
    stage:&'static str,control:&mut SolveControl<'_>) -> Result<Vec<f64>,EvaluationStop> {
    let x = control.solve_preconditioned(op,preconditioner,rhs,0.1*tol,50_000,stage)?;
    control.checkpoint("sdf3-true-residual")?;
    let mut applied = vec![0.0;rhs.len()]; op.apply(&x,&mut applied);
    let scale = rhs.iter().map(|v|v.abs()).fold(0.0_f64,f64::max);
    if scale > 0.0 {
        let mut nr = 0.0; let mut nb = 0.0;
        for (&b,a) in rhs.iter().zip(applied) {
            let r = (b-a)/scale; nr += r*r; nb += (b/scale)*(b/scale);
        }
        if !nr.is_finite() || !(nr <= tol*tol*nb) { return Err(failure("sdf3-true-residual")); }
    } else if applied.iter().any(|a|*a != 0.0) { return Err(failure("sdf3-zero-residual")); }
    control.checkpoint("sdf3-true-residual-publish")?;
    Ok(x)
}

fn trial(rho:&[f64],ratios:&[f64],lambda:f64,step:f64)->Vec<f64> {
    rho.iter().zip(ratios).map(|(&r,&q)| {
        let low=(r-step).max(1e-3);let high=(r+step).min(1.0);let exponent=0.5*(q-lambda);
        if exponent<=fs_math::det::ln(low/r) {low}
        else if exponent>=fs_math::det::ln(high/r) {high}
        else {(r*fs_math::det::exp(exponent)).clamp(low,high)}
    }).collect()
}

/// Volume-feasible, monotone-compliance OC updates on the fixed raw-SDF domain.
/// The initial projected design must be feasible. The result reuses the existing
/// report/termination types; DesignChange is not stationarity. Only complete,
/// independently solved accepted designs replace the current operator and fields.
/// Bisection and rejected solves count against the shared control. No reintegration
/// or remeshing is performed. Parameters are fixed throughout this call.
pub fn controlled_sdf3_optimality_criteria<O:Sdf3Elasticity>(study:&mut CutDensityStudy3<O>,loads:&[LoadCase<'_>],rho0:&[f64],
    options:MultiLoadOcOptions,control:&mut SolveControl<'_>)->MultiLoadOcReport {
    assert_oc_inputs(study, loads, rho0, options);
    let mut report=MultiLoadOcReport {rho:rho0.to_vec(),projected_rho:Vec::new(),displacements:Vec::new(),history:Vec::new(),
        termination:MultiLoadOcTermination::IterationBudget,evaluation_stop:None,work:control.work()};
    let mut accepted_scales=study.operator.scales().to_vec();
    let cap=options.volume_fraction+options.volume_tolerance;
    let outcome=(||->Result<(),EvaluationStop> {
        assert!(study.design(rho0,control)?.volume<=cap,"starting projected density exceeds volume cap");
        let mut current=study.evaluate(rho0,loads,control)?;
        accepted_scales=study.operator.scales().to_vec();
        retain(&mut report,&current,0,0.0);
        for iteration in 1..=options.max_iterations {
            control.checkpoint("sdf3-optimizer")?;
            if !current.volume_gradient.iter().all(|v|v.is_finite()&&*v>0.0)
                || !current.objective.gradient.iter().all(|g|g.is_finite()&&*g<=0.0) {return Err(failure("sdf3-oc-gradient"));}
            let ratios:Vec<f64>=current.objective.gradient.iter().zip(&current.volume_gradient).map(|(&g,&v)|
                if g<0.0 {fs_math::det::ln(-g)-fs_math::det::ln(v)} else {f64::NEG_INFINITY}).collect();
            let low=ratios.iter().copied().filter(|x|x.is_finite()).fold(f64::INFINITY,f64::min)-32.0;
            let high=ratios.iter().copied().filter(|x|x.is_finite()).fold(f64::NEG_INFINITY,f64::max)+32.0;
            if !low.is_finite()||!high.is_finite() {report.termination=MultiLoadOcTermination::NoAcceptableStep;break;}
            let mut accepted=None;let mut step=options.move_limit;
            for _ in 0..=options.max_backtracks {
                control.checkpoint("sdf3-backtrack")?;
                let mut candidate=trial(&report.rho,&ratios,high,step);
                if study.design(&candidate,control)?.volume<=cap {
                    let (mut a,mut b)=(low,high);
                    for _ in 0..64 {
                        control.checkpoint("sdf3-multiplier")?;
                        let mid=f64::midpoint(a,b);
                        if mid<=a||mid>=b {break;}
                        let probe=trial(&report.rho,&ratios,mid,step);
                        if study.design(&probe,control)?.volume>cap {a=mid;} else {b=mid;candidate=probe;}
                    }
                    let evaluated=study.evaluate(&candidate,loads,control)?;
                    if evaluated.volume_fraction<=cap && evaluated.objective.compliance<=current.objective.compliance {
                        accepted=Some((candidate,evaluated));break;
                    }
                }
                step*=0.5;
            }
            let Some((rho,next))=accepted else {report.termination=MultiLoadOcTermination::NoAcceptableStep;break;};
            let change=rho.iter().zip(&report.rho).map(|(a,b)|(a-b).abs()).fold(0.0_f64,f64::max);
            // Poll before publication: a stopped candidate is not an accepted design.
            control.checkpoint("sdf3-accept")?;
            report.rho=rho;current=next;accepted_scales=study.operator.scales().to_vec();
            retain(&mut report,&current,iteration,change);
            if change<=options.change_tolerance {report.termination=MultiLoadOcTermination::DesignChange;break;}
        }
        Ok(())
    })();
    study.operator.set_scales(&accepted_scales).expect("accepted scales remain valid");
    if let Err(stop)=outcome {
        report.termination=match &stop {
            EvaluationStop::Cancelled=>MultiLoadOcTermination::Cancelled,
            EvaluationStop::LinearBudget{..}|EvaluationStop::TotalBudget{..}=>MultiLoadOcTermination::LinearBudget,
            EvaluationStop::Breakdown{..}=>MultiLoadOcTermination::NumericalFailure,
        };
        report.evaluation_stop=Some(stop);
    }
    report.work=control.work();report
}
fn retain(report:&mut MultiLoadOcReport,current:&CutDensityEvaluation3,iteration:usize,change:f64) {
    report.projected_rho.clone_from(&current.projected_rho);
    report.displacements.clone_from(&current.objective.displacements);
    report.history.push(MultiLoadOcIteration {iteration,compliance:current.objective.compliance,
        case_compliances:current.objective.case_compliances.clone(),volume_fraction:current.volume_fraction,max_change:change});
}
