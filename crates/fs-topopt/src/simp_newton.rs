//! Joint live-IR compliance and physical-volume derivatives of raw SIMP densities.
//!
//! Both views borrow one owner: one exact-point linearization, immutable loads,
//! and one lifetime work budget. The volume is the cell-volume-weighted projected
//! density, never penalized stiffness. The native filter and elasticity solvers
//! remain the only derivative authority. This module does not optimize by itself.

use std::cell::{Cell, RefCell};
use std::ops::ControlFlow;

use fs_exec::Cx;
use fs_opt::reverse::physics::{PhysicsError, PhysicsHessian, PhysicsModel, PhysicsSample, PhysicsSignature};
use fs_qty::Dims;

use crate::compliance_newton::{COMPLIANCE_DIMS, ComplianceNewtonBudget, ComplianceNewtonWork};
use crate::elasticity::second_order::{ComplianceError, ComplianceSolveOptions};
use crate::pipeline::LoadCase;
use crate::simp_second_order::SimpComplianceGoal;
use crate::{DensityElasticity, DesignPipeline, EvaluationStop, SolveBudget, SolveControl};

struct Cached<'a> {
    bits: Vec<u64>,
    goal: SimpComplianceGoal<'a>,
    volume: PhysicsSample,
}

/// Joint raw-density physics for compliance objectives / volume constraints, or
/// minimum-volume objectives / compliance constraints. Create the scalar views
/// with [`Self::compliance`] and [`Self::volume`], then bind both into the IR.
///
/// Preparation reserves one forward filter, each independent load, and two
/// transpose filters. Compliance products reserve two filters and each active
/// load tangent; nonlinear volume products reserve only two filters. Reservations
/// are conservative ceilings, not actual solve counts. All Krylov work, including
/// failures, is charged once to this owner rather than to either scalar view.
/// Cached reads and cached zero directions do not spend additional work.
///
/// A failed preparation never replaces the previous cache. A successful trial
/// may evict it; evaluating an older accepted tape re-prepares that exact point.
/// Loads, parameters, geometry and weights cannot change during the model's life.
/// Raw densities are in [0,1]; the native second-order goal additionally requires
/// projected densities strictly inside (1e-12,1). No endpoint Hessian, global byte
/// admission, process checkpoint, or continuum certificate is claimed.
pub struct SimpNewtonModel<'a> {
    pipeline: &'a DesignPipeline,
    elasticity: &'a DensityElasticity,
    loads: Vec<(Vec<f64>, f64)>,
    volumes: Vec<f64>,
    options: ComplianceSolveOptions,
    budget: Cell<ComplianceNewtonBudget>,
    work: Cell<ComplianceNewtonWork>,
    cache: RefCell<Option<Cached<'a>>>,
    last_native_error: RefCell<Option<ComplianceError>>,
}

impl std::fmt::Debug for SimpNewtonModel<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SimpNewtonModel")
            .field("parameters", &self.elasticity.cells())
            .field("loads", &self.loads.len())
            .field("work", &self.work.get())
            .field("budget", &self.budget.get())
            .finish_non_exhaustive()
    }
}

#[derive(Debug, Clone, Copy)]
enum Output { Compliance, Volume }

/// A named scalar view. Cloning a view shares, and cannot reset, its owner's
/// cache or counters. Distinct IR study names are supplied explicitly by callers.
#[derive(Debug, Clone)]
pub struct SimpNewtonGoal<'m, 'a> {
    model: &'m SimpNewtonModel<'a>,
    study: String,
    output: Output,
}

fn poll(cx: Option<&Cx<'_>>) -> Result<(), PhysicsError> {
    if let Some(cx) = cx { cx.checkpoint().map_err(|_| PhysicsError::Cancelled)?; }
    Ok(())
}
fn finite(value: f64, quantity: &'static str) -> Result<f64, PhysicsError> {
    if value.is_finite() { Ok(value) }
    else { Err(PhysicsError::NonFinite { quantity, bits: value.to_bits() }) }
}

impl<'a> SimpNewtonModel<'a> {
    /// Admit immutable modeling data without performing a filter or PDE solve.
    /// `volumes` must be the physical design-cell volumes, in cell order.
    pub fn new(pipeline: &'a DesignPipeline, elasticity: &'a DensityElasticity,
        loads: &[LoadCase<'_>], volumes: &[f64], options: ComplianceSolveOptions,
        budget: ComplianceNewtonBudget,
    ) -> Result<Self, PhysicsError> {
        let n = elasticity.cells();
        if n == 0 || pipeline.filter.cells() != n || volumes.len() != n {
            return Err(PhysicsError::InvalidInput("SIMP mesh, filter and volume shapes must agree"));
        }
        let p = pipeline.params;
        if !(p.e_min.is_finite() && p.e_min > 0.0 && p.e_min < 1.0
            && p.penal.is_finite() && p.penal >= 1.0 && p.beta.is_finite() && p.beta >= 0.0
            && p.eta.is_finite() && (0.0..=1.0).contains(&p.eta)) {
            return Err(PhysicsError::InvalidInput("invalid SIMP/projection parameters"));
        }
        if !(options.tolerance.is_finite() && options.tolerance > 0.0 && options.tolerance < 1.0) {
            return Err(PhysicsError::InvalidInput("relative tolerance must lie in (0, 1)"));
        }
        let mut total = 0.0;
        for &v in volumes {
            if !(v.is_finite() && v > 0.0) {
                return Err(PhysicsError::InvalidInput("physical cell volumes must be finite and positive"));
            }
            total = finite(total + v, "total physical volume")?;
        }
        if loads.is_empty() || !loads.iter().any(|l| l.weight > 0.0) {
            return Err(PhysicsError::InvalidInput("at least one positive-weight load required"));
        }
        for load in loads {
            if load.force.len() != elasticity.n() || !(load.weight.is_finite() && load.weight >= 0.0) {
                return Err(PhysicsError::InvalidInput("invalid structural load shape or weight"));
            }
            for &force in load.force { finite(force, "structural force")?; }
        }
        loads.len().checked_add(3).ok_or(PhysicsError::InvalidInput("too many loads"))?;
        Ok(Self { pipeline, elasticity, loads: loads.iter().map(|l| (l.force.to_vec(), l.weight)).collect(),
            volumes: volumes.to_vec(), options, budget: Cell::new(budget),
            work: Cell::new(ComplianceNewtonWork::default()), cache: RefCell::new(None),
            last_native_error: RefCell::new(None) })
    }

    pub fn compliance(&self, study: impl Into<String>) -> Result<SimpNewtonGoal<'_, 'a>, PhysicsError> {
        self.view(study.into(), Output::Compliance)
    }
    pub fn volume(&self, study: impl Into<String>) -> Result<SimpNewtonGoal<'_, 'a>, PhysicsError> {
        self.view(study.into(), Output::Volume)
    }
    fn view(&self, study: String, output: Output) -> Result<SimpNewtonGoal<'_, 'a>, PhysicsError> {
        if study.is_empty() { return Err(PhysicsError::InvalidInput("nonempty study name required")); }
        Ok(SimpNewtonGoal { model: self, study, output })
    }
    pub fn work(&self) -> ComplianceNewtonWork { self.work.get() }
    pub fn budget(&self) -> ComplianceNewtonBudget { self.budget.get() }
    pub fn set_budget(&self, budget: ComplianceNewtonBudget) { self.budget.set(budget); }
    /// A later success does not erase the most recent native refusal.
    pub fn last_native_error(&self) -> Option<ComplianceError> { self.last_native_error.borrow().clone() }

    fn native_error(&self, error: ComplianceError) -> PhysicsError {
        let result = match &error {
            ComplianceError::InvalidInput(message) => PhysicsError::InvalidInput(*message),
            ComplianceError::Stopped(EvaluationStop::Cancelled) => PhysicsError::Cancelled,
            ComplianceError::Stopped(EvaluationStop::TotalBudget { .. }) => PhysicsError::Budget {
                resource: "SIMP Krylov iterations", used: self.work.get().linear_iterations,
                limit: self.budget.get().linear_iterations,
            },
            ComplianceError::Stopped(EvaluationStop::LinearBudget { stage, iterations, residual_estimate }) =>
                PhysicsError::NotConverged { phase: *stage, iterations: *iterations,
                    residual_bits: residual_estimate.to_bits(), tolerance_bits: self.options.tolerance.to_bits() },
            _ => PhysicsError::Domain("native SIMP refusal; see last_native_error"),
        };
        *self.last_native_error.borrow_mut() = Some(error);
        result
    }

    fn controlled<R>(&self, slots: usize, cx: Option<&Cx<'_>>,
        action: impl FnOnce(&mut SolveControl<'_>) -> Result<R, ComplianceError>,
    ) -> Result<R, PhysicsError> {
        poll(cx)?;
        let mut work = self.work.get(); let budget = self.budget.get();
        if work.solve_reservations > budget.solve_reservations
            || slots > budget.solve_reservations - work.solve_reservations {
            return Err(PhysicsError::Budget { resource: "SIMP solve reservations",
                used: work.solve_reservations, limit: budget.solve_reservations });
        }
        if work.linear_iterations > budget.linear_iterations {
            return Err(PhysicsError::Budget { resource: "SIMP Krylov iterations",
                used: work.linear_iterations, limit: budget.linear_iterations });
        }
        work.solve_reservations += slots;
        self.work.set(work);
        let mut callback = |_| match poll(cx) {
            Ok(()) => ControlFlow::Continue(()), Err(_) => ControlFlow::Break(()),
        };
        let mut control = SolveControl::new(SolveBudget {
            per_solve_iterations: self.options.max_iterations,
            total_iterations: budget.linear_iterations - work.linear_iterations,
        }, &mut callback);
        let result = action(&mut control);
        let spent = control.work();
        debug_assert!(spent.linear_solves <= slots);
        work.started_solves += spent.linear_solves;
        work.linear_iterations += spent.linear_iterations;
        self.work.set(work);
        let result = result.map_err(|e| self.native_error(e))?;
        poll(cx)?;
        Ok(result)
    }

    fn at_point<R>(&self, point: &[f64], cx: Option<&Cx<'_>>,
        action: impl FnOnce(&Cached<'a>) -> Result<R, PhysicsError>,
    ) -> Result<R, PhysicsError> {
        poll(cx)?;
        if point.len() != self.elasticity.cells() {
            return Err(PhysicsError::InvalidInput("raw-density point shape mismatch"));
        }
        for (e, &rho) in point.iter().enumerate() {
            if e % 256 == 0 { poll(cx)?; }
            finite(rho, "raw density")?;
            if !(0.0..=1.0).contains(&rho) {
                return Err(PhysicsError::Domain("raw densities must lie in [0,1]"));
            }
        }
        let mut cache = self.cache.borrow_mut();
        if !cache.as_ref().is_some_and(|c| c.bits.iter().zip(point).all(|(&b, x)| b == x.to_bits())) {
            let loads: Vec<_> = self.loads.iter().map(|(force, weight)| LoadCase { force, weight: *weight }).collect();
            let prepared = self.controlled(loads.len() + 3, cx, |control| {
                let goal = SimpComplianceGoal::new(self.pipeline, self.elasticity, point, &loads, self.options, control)?;
                let (value, gradient) = goal.volume_value_gradient(&self.volumes, control)?;
                Ok(Cached { bits: point.iter().map(|x| x.to_bits()).collect(), goal,
                    volume: PhysicsSample { value, gradient } })
            })?;
            poll(cx)?;
            *cache = Some(prepared);
        }
        let result = action(cache.as_ref().expect("successful SIMP preparation"))?;
        poll(cx)?;
        Ok(result)
    }
}

impl PhysicsModel for SimpNewtonGoal<'_, '_> {
    fn signature(&self) -> PhysicsSignature<'_> {
        PhysicsSignature { study: &self.study, parameter_count: self.model.elasticity.cells(),
            parameter_dims: Dims::NONE, residual_dims: match self.output {
                Output::Compliance => COMPLIANCE_DIMS, Output::Volume => Dims::NONE,
            } }
    }
    fn value_gradient(&self, point: &[f64], cx: Option<&Cx<'_>>) -> Result<PhysicsSample, PhysicsError> {
        self.model.at_point(point, cx, |c| Ok(match self.output {
            Output::Compliance => PhysicsSample { value: c.goal.value(), gradient: c.goal.gradient().to_vec() },
            Output::Volume => c.volume.clone(),
        }))
    }
    fn second_order(&self) -> Option<&dyn PhysicsHessian> { Some(self) }
}
impl PhysicsHessian for SimpNewtonGoal<'_, '_> {
    fn hessian_vector_product(&self, point: &[f64], direction: &[f64], cx: Option<&Cx<'_>>)
        -> Result<Vec<f64>, PhysicsError>
    {
        poll(cx)?;
        if direction.len() != self.model.elasticity.cells() {
            return Err(PhysicsError::InvalidInput("raw-density Hessian direction shape mismatch"));
        }
        for (e, &v) in direction.iter().enumerate() {
            if e % 256 == 0 { poll(cx)?; } finite(v, "raw-density direction")?;
        }
        self.model.at_point(point, cx, |c| {
            if direction.iter().all(|&v| v == 0.0) { return Ok(vec![0.0; direction.len()]); }
            match self.output {
                Output::Compliance => {
                    let slots = 2 + self.model.loads.iter().filter(|(_, w)| *w > 0.0).count();
                    self.model.controlled(slots, cx, |control| c.goal.hessian_vector_product(direction, control))
                }
                Output::Volume => {
                    if self.model.pipeline.params.beta == 0.0 { return Ok(vec![0.0; direction.len()]); }
                    self.model.controlled(2, cx, |control|
                        c.goal.volume_hessian_vector_product(&self.model.volumes, direction, control))
                }
            }
        })
    }
}
