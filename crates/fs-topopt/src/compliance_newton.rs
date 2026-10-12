//! Native, live-IR structural design in log-stiffness coordinates.
//!
//! The base tetrahedral operator supplies fixed geometry, homogeneous supports
//! and a reference elastic material. For `m = exp(x)`, the coordinate Hessian
//! is `diag(m) H_m diag(m) + diag(m * g_m)`, NOT just the first term.
//! SI geometry, force and reference modulus give compliance in joules.
//! This is stiffness sizing, not a SIMP/filter/shape Hessian or a stress model.

use std::cell::{Cell, RefCell};
use std::ops::ControlFlow;
use fs_exec::Cx;
use fs_opt::reverse::physics::{PhysicsError, PhysicsHessian, PhysicsModel, PhysicsSample, PhysicsSignature};
use fs_qty::Dims;
use crate::{DensityElasticity, EvaluationStop, SolveBudget, SolveControl};
use crate::elasticity::second_order::{ComplianceError, ComplianceGoal, ComplianceSolveOptions};
use crate::pipeline::LoadCase;

/// Work-conjugate force times displacement, not displacement or squared force.
pub const COMPLIANCE_DIMS: Dims = Dims([2, 1, -2, 0, 0, 0]);

/// Absolute lifetime ceilings. Raising a ceiling never refunds admitted work.
#[derive(Debug, Clone, Copy)]
pub struct ComplianceNewtonBudget {
    /// Conservative reservation: every prepared load, or every positive-weight
    /// tangent, costs one slot. Failed or cancelled reservations are not refunded.
    pub solve_reservations: usize,
    /// Actual cumulative CG iterations, including incomplete evaluations.
    pub linear_iterations: usize,
}

/// Reservation counts and actual started solves are intentionally distinct.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ComplianceNewtonWork {
    pub solve_reservations: usize,
    pub started_solves: usize,
    pub linear_iterations: usize,
}

struct Cached<'a> {
    bits: Vec<u64>,
    moduli: Vec<f64>,
    goal: ComplianceGoal<'a>,
}

/// A real native compliance provider for `PdeResidual` and constrained Newton.
///
/// The model owns a single exact-point cache and cumulative physics work. A
/// rejected trial may evict it; an older accepted tape will re-prepare its own
/// exact point, never use a trial Hessian. Failed preparations leave the previous
/// cache intact. Cached values/gradients and cached zero directions need no
/// additional solves. Cloning an optimizer cannot reset these provider counters.
/// `None` explicitly disables Cx cancellation; work ceilings still apply.
///
/// This adapter does not impose a global byte cap or serialize a checkpoint.
/// Results are numerical fixed-mesh responses, not continuum certificates.
pub struct ComplianceNewtonModel<'a> {
    elasticity: &'a DensityElasticity,
    study: String,
    loads: Vec<(Vec<f64>, f64)>,
    options: ComplianceSolveOptions,
    budget: Cell<ComplianceNewtonBudget>,
    work: Cell<ComplianceNewtonWork>,
    cache: RefCell<Option<Cached<'a>>>,
    last_native_error: RefCell<Option<ComplianceError>>,
}

impl std::fmt::Debug for ComplianceNewtonModel<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ComplianceNewtonModel").field("study", &self.study)
            .field("parameters", &self.elasticity.cells()).field("loads", &self.loads.len())
            .field("work", &self.work.get()).field("budget", &self.budget.get())
            .finish_non_exhaustive()
    }
}
fn poll(cx: Option<&Cx<'_>>) -> Result<(), PhysicsError> {
    if let Some(cx) = cx { cx.checkpoint().map_err(|_| PhysicsError::Cancelled)?; }
    Ok(())
}
fn finite(value: f64, quantity: &'static str) -> Result<f64, PhysicsError> {
    if value.is_finite() { Ok(value) }
    else { Err(PhysicsError::NonFinite { quantity, bits: value.to_bits() }) }
}

impl<'a> ComplianceNewtonModel<'a> {
    /// Validate and retain immutable load data without assembling or solving.
    /// The initial public `elasticity.moduli` does not define this parameter
    /// mapping: each `exp(x_e)` multiplies the original element block directly.
    pub fn new(elasticity: &'a DensityElasticity, study: impl Into<String>,
        loads: &[LoadCase<'_>], options: ComplianceSolveOptions, budget: ComplianceNewtonBudget,
    ) -> Result<Self, PhysicsError> {
        let study = study.into();
        if study.is_empty() || elasticity.cells() == 0 {
            return Err(PhysicsError::InvalidInput("nonempty study and elasticity mesh required"));
        }
        if !(options.tolerance.is_finite() && options.tolerance > 0.0 && options.tolerance < 1.0) {
            return Err(PhysicsError::InvalidInput("relative tolerance must lie in (0, 1)"));
        }
        if loads.is_empty() || !loads.iter().any(|load| load.weight > 0.0) {
            return Err(PhysicsError::InvalidInput("at least one positive-weight structural load required"));
        }
        for load in loads {
            if load.force.len() != elasticity.n() || !(load.weight.is_finite() && load.weight >= 0.0) {
                return Err(PhysicsError::InvalidInput("structural load shape or weight is invalid"));
            }
            for &f in load.force { finite(f, "structural force")?; }
        }
        Ok(Self { elasticity, study, loads: loads.iter().map(|l| (l.force.to_vec(), l.weight)).collect(),
            options, budget: Cell::new(budget), work: Cell::new(ComplianceNewtonWork::default()),
            cache: RefCell::new(None), last_native_error: RefCell::new(None) })
    }

    pub fn work(&self) -> ComplianceNewtonWork { self.work.get() }
    pub fn budget(&self) -> ComplianceNewtonBudget { self.budget.get() }
    pub fn set_budget(&self, budget: ComplianceNewtonBudget) { self.budget.set(budget); }
    /// Most recent native refusal; a later success does not erase it.
    pub fn last_native_error(&self) -> Option<ComplianceError> { self.last_native_error.borrow().clone() }

    fn native_error(&self, error: ComplianceError) -> PhysicsError {
        let result = match &error {
            ComplianceError::InvalidInput(message) => PhysicsError::InvalidInput(*message),
            ComplianceError::Stopped(EvaluationStop::Cancelled) => PhysicsError::Cancelled,
            ComplianceError::Stopped(EvaluationStop::TotalBudget { .. }) => PhysicsError::Budget {
                resource: "structural Krylov iterations", used: self.work.get().linear_iterations,
                limit: self.budget.get().linear_iterations,
            },
            ComplianceError::Stopped(EvaluationStop::LinearBudget { stage, iterations, residual_estimate }) =>
                PhysicsError::NotConverged { phase: *stage, iterations: *iterations,
                    residual_bits: residual_estimate.to_bits(), tolerance_bits: self.options.tolerance.to_bits() },
            // Preserve the precise native error rather than inventing a missing
            // iteration count or non-finite bit pattern in the common protocol.
            _ => PhysicsError::Domain("native structural refusal; see last_native_error"),
        };
        *self.last_native_error.borrow_mut() = Some(error);
        result
    }

    fn controlled<R>(&self, slots: usize, cx: Option<&Cx<'_>>,
        action: impl FnOnce(&mut SolveControl<'_>) -> Result<R, ComplianceError>,
    ) -> Result<R, PhysicsError> {
        poll(cx)?;
        let mut work = self.work.get();
        let budget = self.budget.get();
        if work.solve_reservations > budget.solve_reservations
            || slots > budget.solve_reservations - work.solve_reservations {
            return Err(PhysicsError::Budget { resource: "structural solve reservations",
                used: work.solve_reservations, limit: budget.solve_reservations });
        }
        if work.linear_iterations > budget.linear_iterations {
            return Err(PhysicsError::Budget { resource: "structural Krylov iterations",
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
        // Each started solve has already been reserved; each completed iteration
        // was admitted by the remaining absolute ceiling. Neither sum overflows.
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
            return Err(PhysicsError::InvalidInput("log-stiffness point shape mismatch"));
        }
        let mut moduli = Vec::with_capacity(point.len());
        for (e, &x) in point.iter().enumerate() {
            if e % 256 == 0 { poll(cx)?; }
            finite(x, "log-stiffness coordinate")?;
            let m = finite(fs_math::det::exp(x), "stiffness multiplier")?;
            if m <= 0.0 { return Err(PhysicsError::Domain("stiffness multiplier underflowed to zero")); }
            moduli.push(m);
        }
        let mut cache = self.cache.borrow_mut();
        if !cache.as_ref().is_some_and(|c| c.bits.iter().zip(point).all(|(&bits, x)| bits == x.to_bits())) {
            let loads: Vec<LoadCase<'_>> = self.loads.iter().map(|(force, weight)| LoadCase { force, weight: *weight }).collect();
            let goal = self.controlled(loads.len(), cx, |control|
                ComplianceGoal::new(self.elasticity, &moduli, &loads, self.options, control))?;
            poll(cx)?;
            *cache = Some(Cached { bits: point.iter().map(|x| x.to_bits()).collect(), moduli, goal });
        }
        let result = action(cache.as_ref().expect("successful structural preparation"))?;
        poll(cx)?;
        Ok(result)
    }
}

impl PhysicsModel for ComplianceNewtonModel<'_> {
    fn signature(&self) -> PhysicsSignature<'_> {
        PhysicsSignature { study: &self.study, parameter_count: self.elasticity.cells(),
            parameter_dims: Dims::NONE, residual_dims: COMPLIANCE_DIMS }
    }
    fn value_gradient(&self, point: &[f64], cx: Option<&Cx<'_>>) -> Result<PhysicsSample, PhysicsError> {
        self.at_point(point, cx, |c| {
            let mut gradient = Vec::with_capacity(c.moduli.len());
            for (e, (&m, &g)) in c.moduli.iter().zip(c.goal.gradient()).enumerate() {
                if e % 256 == 0 { poll(cx)?; }
                gradient.push(finite(m * g, "log-stiffness gradient")?);
            }
            Ok(PhysicsSample { value: c.goal.value(), gradient })
        })
    }
    fn second_order(&self) -> Option<&dyn PhysicsHessian> { Some(self) }
}
impl PhysicsHessian for ComplianceNewtonModel<'_> {
    fn hessian_vector_product(&self, point: &[f64], direction: &[f64], cx: Option<&Cx<'_>>)
        -> Result<Vec<f64>, PhysicsError>
    {
        poll(cx)?;
        if direction.len() != self.elasticity.cells() {
            return Err(PhysicsError::InvalidInput("log-stiffness Hessian direction shape mismatch"));
        }
        for (e, &v) in direction.iter().enumerate() {
            if e % 256 == 0 { poll(cx)?; } finite(v, "Hessian direction")?;
        }
        self.at_point(point, cx, |c| {
            if direction.iter().all(|&v| v == 0.0) { return Ok(vec![0.0; direction.len()]); }
            let mut tangent = Vec::with_capacity(direction.len());
            for (e, (&m, &v)) in c.moduli.iter().zip(direction).enumerate() {
                if e % 256 == 0 { poll(cx)?; }
                tangent.push(finite(m * v, "structural modulus tangent")?);
            }
            let slots = self.loads.iter().filter(|(_, weight)| *weight > 0.0).count();
            let physical = self.controlled(slots, cx, |control| c.goal.hessian_vector_product(&tangent, control))?;
            let mut product = Vec::with_capacity(direction.len());
            for e in 0..direction.len() {
                if e % 256 == 0 { poll(cx)?; }
                product.push(finite(c.moduli[e] * (physical[e] + c.goal.gradient()[e] * direction[e]),
                    "log-stiffness Hessian action")?);
            }
            Ok(product)
        })
    }
}
