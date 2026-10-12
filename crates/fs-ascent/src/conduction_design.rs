//! Native linear thermal design for the live optimization IR and Newton study.
//! `rho = exp(x)` keeps every element conductivity multiplier positive. The
//! exact coordinate Hessian includes BOTH `rho * H_rho(rho*v)` and
//! `rho * gradient_rho * v`; omitting the latter is not Newton in log space.

use fs_conduction::adjoint::{ConductivityDesign, second_order::ConductivityGoal};
use fs_conduction::{ConductionError, TEMPERATURE_DIMS};
use fs_exec::Cx;
use fs_opt::reverse::physics::{PhysicsError, PhysicsHessian, PhysicsModel, PhysicsSample, PhysicsSignature};
use std::cell::{Cell, RefCell};

struct Cached<'d, 'm> {
    point_bits: Vec<u64>,
    rho: Vec<f64>,
    goal: ConductivityGoal<'d, 'm>,
}

/// A scalar weighted free-node temperature, in kelvin, over log-conductivity
/// coordinates. Weights are dimensionless and fixed for the model's lifetime.
/// Bind to one `PdeResidual` with the same study name and `Rn(element_count)`.
/// A real `Cx` is required, including for primal evaluation.
///
/// One point-keyed linearization is retained. Cache misses reserve two linear
/// solve slots (primal + adjoint); nonzero Hessian actions reserve two more
/// (tangent + incremental adjoint). These are conservative reservations, NOT
/// measured successful solves or Krylov iteration counts. Failed preparations,
/// cancellations after admission, zero RHS solves and failed actions receive
/// no refund. Cached values/gradients and zero Hessian directions reserve no
/// additional sensitivity slots. Cloning the optimizer cannot reset this
/// provider-owned cumulative budget. Increase the ceiling explicitly to resume.
///
/// Rejected trials can evict a cached point, but cannot cause a Hessian at the
/// wrong point: older accepted tapes re-prepare using their exact point bits.
/// The previous cache remains intact if the replacement preparation fails.
/// No global memory admission, process checkpoint, continuum error bound,
/// nonlinear k(T), shape derivative or material-card selection is implied.
pub struct ConductivityNewtonModel<'d, 'm> {
    design: &'d ConductivityDesign<'m>,
    study: String,
    weights: Vec<f64>,
    cache: RefCell<Option<Cached<'d, 'm>>>,
    reserved_solves: Cell<usize>,
    solve_limit: Cell<usize>,
    last_native_error: RefCell<Option<ConductionError>>,
}
impl core::fmt::Debug for ConductivityNewtonModel<'_, '_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("ConductivityNewtonModel").field("study",&self.study)
            .field("parameters",&self.design.parameter_count())
            .field("reserved_solves",&self.reserved_solves.get())
            .field("solve_limit",&self.solve_limit.get()).finish_non_exhaustive()
    }
}
fn poll(cx: &Cx<'_>) -> Result<(),PhysicsError> {
    cx.checkpoint().map_err(|_| PhysicsError::Cancelled)
}
fn finite(quantity: &'static str, value: f64) -> Result<f64,PhysicsError> {
    if value.is_finite() { Ok(value) }
    else { Err(PhysicsError::NonFinite { quantity,bits: value.to_bits() }) }
}
impl<'d, 'm> ConductivityNewtonModel<'d, 'm> {
    /// Construct without assembly or a solve. A zero allowance is valid and
    /// refuses the first uncached numerical call, rather than implying success.
    pub fn new(design: &'d ConductivityDesign<'m>, study: impl Into<String>,
        weights: Vec<f64>, solve_limit: usize,
    ) -> Result<Self,PhysicsError> {
        let study=study.into();
        if study.is_empty() { return Err(PhysicsError::InvalidInput("empty thermal study name")); }
        if weights.len()!=design.dofs().n() {
            return Err(PhysicsError::InvalidInput("thermal goal weights must use free-dof order"));
        }
        for &w in &weights { finite("thermal goal weight",w)?; }
        Ok(Self { design,study,weights,cache: RefCell::new(None),
            reserved_solves: Cell::new(0),solve_limit: Cell::new(solve_limit),
            last_native_error: RefCell::new(None) })
    }
    pub fn reserved_solves(&self) -> usize { self.reserved_solves.get() }
    pub fn solve_limit(&self) -> usize { self.solve_limit.get() }
    /// Change the absolute ceiling, never the spent counter or model meaning.
    pub fn set_solve_limit(&self, limit: usize) { self.solve_limit.set(limit); }
    /// Preserve the complete native refusal beyond the shared physics protocol.
    pub fn last_native_error(&self) -> Option<ConductionError> { self.last_native_error.borrow().clone() }

    fn reserve(&self, slots: usize) -> Result<(),PhysicsError> {
        let used=self.reserved_solves.get();let limit=self.solve_limit.get();
        if used>limit || slots>limit-used {
            return Err(PhysicsError::Budget { resource: "thermal linear solve reservations",used,limit });
        }
        self.reserved_solves.set(used+slots);
        Ok(())
    }
    fn native_error(&self, error: ConductionError) -> PhysicsError {
        let result=match &error {
            ConductionError::Cancelled { .. } => PhysicsError::Cancelled,
            ConductionError::NonFinite { field,bits } => PhysicsError::NonFinite { quantity: *field,bits: *bits },
            ConductionError::LinearSolveFailed { krylov_iterations,true_relative_residual,tolerance,.. } =>
                PhysicsError::NotConverged { phase: "thermal primal/adjoint/sensitivity",
                    iterations: *krylov_iterations,residual_bits: true_relative_residual.to_bits(),
                    tolerance_bits: tolerance.to_bits() },
            ConductionError::Config { .. } | ConductionError::FieldLength { .. } =>
                PhysicsError::InvalidInput("native thermal configuration; see last_native_error"),
            _ => PhysicsError::Domain("native thermal model refusal; see last_native_error"),
        };
        *self.last_native_error.borrow_mut()=Some(error);
        result
    }
    fn at_point<R>(&self, point: &[f64], cx: &Cx<'_>,
        action: impl FnOnce(&Cached<'d,'m>) -> Result<R,PhysicsError>,
    ) -> Result<R,PhysicsError> {
        poll(cx)?;
        if point.len()!=self.design.parameter_count() {
            return Err(PhysicsError::InvalidInput("log-conductivity point length mismatch"));
        }
        let mut rho=Vec::with_capacity(point.len());
        for (i,&x) in point.iter().enumerate() {
            if i%256==0 { poll(cx)?; }
            finite("log-conductivity coordinate",x)?;
            let r=finite("conductivity multiplier",fs_math::det::exp(x))?;
            if r<=0.0 { return Err(PhysicsError::Domain("conductivity multiplier underflowed to zero")); }
            rho.push(r);
        }
        let mut cache=self.cache.borrow_mut();
        let matches=cache.as_ref().is_some_and(|c| c.point_bits.iter().zip(point)
            .all(|(&bits,x)| bits==x.to_bits()));
        if !matches {
            self.reserve(2)?;
            let goal=self.design.linearize_goal(cx,&rho,&self.weights).map_err(|e| self.native_error(e))?;
            poll(cx)?;
            *cache=Some(Cached { point_bits: point.iter().map(|x| x.to_bits()).collect(),rho,goal });
        }
        let result=action(cache.as_ref().expect("successful point preparation"))?;
        poll(cx)?;
        Ok(result)
    }
}
impl PhysicsModel for ConductivityNewtonModel<'_, '_> {
    fn signature(&self) -> PhysicsSignature<'_> {
        PhysicsSignature { study: &self.study,parameter_count: self.design.parameter_count(),
            parameter_dims: TEMPERATURE_DIMS.checked_times(0).expect("zero dimension exponents"),
            residual_dims: TEMPERATURE_DIMS }
    }
    fn value_gradient(&self, point: &[f64], cx: Option<&Cx<'_>>) -> Result<PhysicsSample,PhysicsError> {
        let cx=cx.ok_or(PhysicsError::InvalidInput("native thermal design requires a Cx"))?;
        self.at_point(point,cx,|c| {
            let mut gradient=Vec::with_capacity(c.rho.len());
            for (i,(&r,&g)) in c.rho.iter().zip(c.goal.gradient()).enumerate() {
                if i%256==0 { poll(cx)?; }
                gradient.push(finite("log-conductivity gradient",r*g)?);
            }
            Ok(PhysicsSample { value: c.goal.value(),gradient })
        })
    }
    fn second_order(&self) -> Option<&dyn PhysicsHessian> { Some(self) }
}
impl PhysicsHessian for ConductivityNewtonModel<'_, '_> {
    fn hessian_vector_product(&self, point: &[f64], direction: &[f64], cx: Option<&Cx<'_>>)
        -> Result<Vec<f64>,PhysicsError>
    {
        let cx=cx.ok_or(PhysicsError::InvalidInput("native thermal design requires a Cx"))?;
        poll(cx)?;
        if direction.len()!=self.design.parameter_count() {
            return Err(PhysicsError::InvalidInput("log-conductivity Hessian direction length mismatch"));
        }
        for (i,&v) in direction.iter().enumerate() { if i%256==0 { poll(cx)?; } finite("Hessian direction",v)?; }
        self.at_point(point,cx,|c| {
            if direction.iter().all(|&v| v==0.0) { return Ok(vec![0.0;direction.len()]); }
            let mut physical=Vec::with_capacity(direction.len());
            for (i,(&r,&v)) in c.rho.iter().zip(direction).enumerate() {
                if i%256==0 { poll(cx)?; }
                physical.push(finite("conductivity tangent",r*v)?);
            }
            self.reserve(2)?;
            let hv=c.goal.hessian_vector_product(cx,&physical).map_err(|e| self.native_error(e))?;
            let mut product=Vec::with_capacity(direction.len());
            for i in 0..direction.len() {
                if i%256==0 { poll(cx)?; }
                product.push(finite("log-conductivity Hessian action",
                    c.rho[i]*(hv.values[i]+c.goal.gradient()[i]*direction[i]))?);
            }
            Ok(product)
        })
    }
}
