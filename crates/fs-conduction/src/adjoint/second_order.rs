//! Matrix-free design Hessians for the production linear P1 thermal model.
//!
//! For A(rho) T = b(rho), J = w^T T_free, and fixed w/boundary data:
//! A dT = -(sum_e v_e K_e) T_full |_free;
//! A^T lambda = w; A^T dlambda = -(sum_e v_e K_e) lambda_full |_free;
//! (Hv)_e = -dlambda^T K_e T_full - lambda^T K_e dT_full.
//! The Dirichlet lift is differentiated exactly: dT and both adjoints vanish
//! on prescribed vertices, but the primal in K_e T uses the FULL temperature.
//! No differentiated Krylov iterations, finite differences, dense parameter
//! Hessian, shape derivative, k(T) derivative or continuum certificate.

use super::{ConductivityDesign, DesignSolution};
use crate::{ConductionError, ASSEMBLY_TILE};
use fs_exec::Cx;
use fs_solver::{CgState, CsrOp, norm2};

/// Recomputed relative residual and work for one linear sensitivity solve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SensitivitySolve {
    pub relative_residual: f64,
    pub iterations: usize,
}

/// A complete direction, including its physical temperature response.
#[derive(Debug, Clone, PartialEq)]
pub struct ConductivityHessianProduct {
    /// One total second derivative per element multiplier.
    pub values: Vec<f64>,
    /// Full nodal dT, exactly zero on the prescribed boundary.
    pub temperature_direction: Vec<f64>,
    pub tangent: SensitivitySolve,
    pub incremental_adjoint: SensitivitySolve,
}

/// Immutable primal/adjoint pair and the SAME assembled operator at one design.
/// Retain this object across Krylov products: each nonzero direction requires
/// two linear solves, not another primal, adjoint or sparse assembly. Storage
/// is sparse in spatial dofs and linear in element parameters; no parameter
/// Jacobian or Hessian is assembled. The caller owns cumulative solve budgets.
/// Each solve is limited by the design's LinearConfig and polls Cx between
/// Krylov iterations. Assembly, sparse matvec and allocation phases retain
/// their production granularity; no whole-process byte/latency bound is implied.
pub struct ConductivityGoal<'d, 'm> {
    design: &'d ConductivityDesign<'m>,
    matrix: fs_sparse::Csr,
    solution: DesignSolution,
    value: f64,
    gradient: Vec<f64>,
    lambda: Vec<f64>,
    dual: SensitivitySolve,
}
impl core::fmt::Debug for ConductivityGoal<'_, '_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("ConductivityGoal").field("value",&self.value)
            .field("parameters",&self.gradient.len()).field("dual",&self.dual)
            .finish_non_exhaustive()
    }
}
fn poll(cx: &Cx<'_>, at: usize) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled { stage: "conductivity-hessian",at })
}
fn check(field: &'static str, values: &[f64], expected: usize) -> Result<(), ConductionError> {
    if values.len()!=expected {
        return Err(ConductionError::FieldLength { field,expected,found: values.len() });
    }
    if let Some(value) = values.iter().find(|v| !v.is_finite()) {
        return Err(ConductionError::NonFinite { field,bits: value.to_bits() });
    }
    Ok(())
}
fn finite(field: &'static str, value: f64) -> Result<f64, ConductionError> {
    if value.is_finite() { Ok(value) }
    else { Err(ConductionError::NonFinite { field,bits: value.to_bits() }) }
}

/// Normalize RHS units before Krylov norms, measure the TRUE residual, and
/// never accept NaN, a stalled recurrence or a merely recursive residual.
fn solve(cx: &Cx<'_>, matrix: &fs_sparse::Csr, rhs: &[f64],
    config: crate::LinearConfig,
) -> Result<(Vec<f64>,SensitivitySolve),ConductionError> {
    poll(cx,0)?;
    check("sensitivity right hand side",rhs,rhs.len())?;
    let scale = rhs.iter().map(|x| x.abs()).fold(0.0_f64,f64::max);
    if scale==0.0 {
        return Ok((vec![0.0;rhs.len()],SensitivitySolve { relative_residual: 0.0,iterations: 0 }));
    }
    let normalized: Vec<_> = rhs.iter().map(|x| x/scale).collect();
    let op = CsrOp::symmetric(matrix.clone());
    let pre = crate::solve::spd_preconditioner(matrix);
    let mut state = CgState::new(&op,&pre,&normalized);
    while state.rel_residual()>=config.tolerance && state.iters<config.max_iterations {
        poll(cx,state.iters)?;
        let before = state.iters;
        state.run(&op,&pre,config.tolerance,1);
        if state.iters==before { break; }
    }
    poll(cx,state.iters)?;
    check("normalized sensitivity solution",&state.x,rhs.len())?;
    let mut applied = vec![0.0;rhs.len()]; matrix.spmv(&state.x,&mut applied);
    let residual: Vec<_> = normalized.iter().zip(&applied).map(|(b,a)| b-a).collect();
    let relative = norm2(&residual)/norm2(&normalized);
    if !relative.is_finite() || relative>=config.tolerance {
        return Err(ConductionError::LinearSolveFailed { iteration: 0,
            krylov_iterations: state.iters,true_relative_residual: relative,tolerance: config.tolerance });
    }
    let values: Vec<_> = state.x.iter().map(|x| x*scale).collect();
    check("sensitivity solution",&values,rhs.len())?;
    poll(cx,state.iters)?;
    Ok((values,SensitivitySolve { relative_residual: relative,iterations: state.iters }))
}

fn homogeneous(design: &ConductivityDesign<'_>, free: &[f64]) -> Vec<f64> {
    let mut full = vec![0.0;design.problem.mesh.vertex_count()];
    for (&v,&x) in design.dofs.free().iter().zip(free) { full[v]=x; }
    full
}
fn element_action(design: &ConductivityDesign<'_>, e: usize, full: &[f64]) -> [f64;4] {
    let mesh = design.problem.mesh;
    let tet = mesh.complex().tets[e];
    let ke = crate::assemble::element_stiffness(mesh,e,&design.element_tensors[e]);
    let mut values = [0.0;4];
    for a in 0..4 { for b in 0..4 {
        values[a] = ke[a][b].mul_add(full[tet[b] as usize],values[a]);
    }}
    values
}
fn pullback(cx: &Cx<'_>, design: &ConductivityDesign<'_>, left: &[f64], full: &[f64])
    -> Result<Vec<f64>,ConductionError>
{
    let mut values = vec![0.0;design.parameter_count()];
    for (e,value) in values.iter_mut().enumerate() {
        if e%ASSEMBLY_TILE==0 { poll(cx,e)?; }
        let tet = design.problem.mesh.complex().tets[e];
        let applied = element_action(design,e,full);
        for a in 0..4 { if let Some(i)=design.dofs.slot_of(tet[a] as usize) {
            *value = left[i].mul_add(applied[a],*value);
        }}
        finite("conductivity derivative contraction",*value)?;
    }
    Ok(values)
}
fn directional_rhs(cx: &Cx<'_>, design: &ConductivityDesign<'_>, direction: &[f64], full: &[f64])
    -> Result<Vec<f64>,ConductionError>
{
    let mut rhs = vec![0.0;design.dofs.n()];
    for (e,&v) in direction.iter().enumerate() {
        if e%ASSEMBLY_TILE==0 { poll(cx,e)?; }
        if v==0.0 { continue; }
        let tet = design.problem.mesh.complex().tets[e];
        let applied = element_action(design,e,full);
        for a in 0..4 { if let Some(i)=design.dofs.slot_of(tet[a] as usize) {
            rhs[i] = finite("conductivity directional residual",(-v).mul_add(applied[a],rhs[i]))?;
        }}
    }
    Ok(rhs)
}

impl<'m> ConductivityDesign<'m> {
    /// Prepare J=sum(weights*T_free), its gradient, and repeatable exact
    /// implicit Hessian actions. Uses this design's admitted fixed anisotropic
    /// element tensors, fixed source and Robin/Dirichlet data. Weights must be
    /// finite and use free-dof order; they may be signed or identically zero.
    /// Performs one assembly and two residual-checked linear solves. Provider
    /// or optimizer code must charge those solves before calling this method.
    pub fn linearize_goal<'d>(&'d self, cx: &Cx<'_>, rho: &[f64], weights: &[f64])
        -> Result<ConductivityGoal<'d,'m>,ConductionError>
    {
        poll(cx,0)?;
        if !self.linear.tolerance.is_finite() || self.linear.tolerance<=0.0
            || self.linear.tolerance>=1.0 || self.linear.max_iterations==0 {
            return Err(ConductionError::Config { parameter: "second-order linear policy",
                what: "require tolerance in (0,1) and positive iteration allowance".to_string() });
        }
        self.check_rho(rho)?;
        check("goal weights",weights,self.dofs.n())?;
        let (matrix,rhs) = self.system(cx,rho)?;
        let (free,primal) = solve(cx,&matrix,&rhs,self.linear)?;
        let temperature = self.dofs.scatter(&free);
        self.check_temperature_support(&temperature)?;
        let value = finite("conductivity goal",weights.iter().zip(&free).map(|(w,t)| w*t).sum())?;
        let (lambda,dual) = solve(cx,&matrix,weights,self.linear)?;
        let mut gradient = pullback(cx,self,&lambda,&temperature)?;
        for value in &mut gradient { *value = -*value; }
        poll(cx,0)?;
        Ok(ConductivityGoal { design: self,matrix,value,gradient,lambda,dual,
            solution: DesignSolution { free_temperature: free,temperature,
                primal_relative_residual: primal.relative_residual,primal_iterations: primal.iterations } })
    }
}
impl ConductivityGoal<'_, '_> {
    pub fn value(&self) -> f64 { self.value }
    pub fn gradient(&self) -> &[f64] { &self.gradient }
    pub fn solution(&self) -> &DesignSolution { &self.solution }
    pub fn dual_solve(&self) -> SensitivitySolve { self.dual }

    /// Apply the total reduced Hessian, not a Gauss--Newton approximation.
    /// No primal solve or sparse assembly is repeated. The direction is in
    /// rho coordinates, may have either sign, and does not change the point.
    /// A failure or cancellation publishes no partial product and leaves this
    /// linearization usable with a new Cx. Nonzero directions cost two solves;
    /// the identically-zero direction costs none.
    pub fn hessian_vector_product(&self, cx: &Cx<'_>, direction: &[f64])
        -> Result<ConductivityHessianProduct,ConductionError>
    {
        poll(cx,0)?;
        check("conductivity Hessian direction",direction,self.design.parameter_count())?;
        let zero = SensitivitySolve { relative_residual: 0.0,iterations: 0 };
        if direction.iter().all(|&v| v==0.0) {
            return Ok(ConductivityHessianProduct { values: vec![0.0;direction.len()],
                temperature_direction: vec![0.0;self.solution.temperature.len()],
                tangent: zero,incremental_adjoint: zero });
        }
        let rhs = directional_rhs(cx,self.design,direction,&self.solution.temperature)?;
        let (dt,tangent) = solve(cx,&self.matrix,&rhs,self.design.linear)?;
        let full_dt = homogeneous(self.design,&dt);
        let full_lambda = homogeneous(self.design,&self.lambda);
        let dual_rhs = directional_rhs(cx,self.design,direction,&full_lambda)?;
        let (dlambda,incremental_adjoint) = solve(cx,&self.matrix,&dual_rhs,self.design.linear)?;
        let mut values = pullback(cx,self.design,&dlambda,&self.solution.temperature)?;
        let other = pullback(cx,self.design,&self.lambda,&full_dt)?;
        for (e,(value,term)) in values.iter_mut().zip(other).enumerate() {
            if e%ASSEMBLY_TILE==0 { poll(cx,e)?; }
            *value = finite("conductivity Hessian product",-*value-term)?;
        }
        poll(cx,0)?;
        Ok(ConductivityHessianProduct { values,temperature_direction: full_dt,tangent,incremental_adjoint })
    }
}
