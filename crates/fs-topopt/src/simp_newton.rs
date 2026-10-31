//! Joint live-IR compliance, physical-volume and material-mass derivatives.
//!
//! Both views borrow one owner: one exact-point linearization, immutable loads,
//! and one lifetime work budget. Ordinary SIMP shares a projection; robust SIMP
//! uses eroded compliance and dilated volume on one common filtered field.
//! Volume is always cell-volume-weighted projected density, never stiffness.

use std::cell::{Cell, RefCell};
use std::ops::ControlFlow;
use std::rc::Rc;

use fs_exec::Cx;
use fs_opt::reverse::physics::{PhysicsError, PhysicsHessian, PhysicsModel, PhysicsSample, PhysicsSignature};
use fs_qty::Dims;

use crate::compliance_newton::{COMPLIANCE_DIMS, ComplianceNewtonBudget, ComplianceNewtonWork};
use crate::elasticity::second_order::{ComplianceError, ComplianceSolveOptions};
use crate::filter::DensityFilter;
use crate::pipeline::{LoadCase, SimpParams};
use crate::simp_second_order::{SimpComplianceGoal, robust_parameters};
use crate::{DensityElasticity, DesignPipeline, EvaluationStop, RobustPipeline, SolveBudget, SolveControl};

/// SI kilograms for a physical material inventory objective or constraint.
pub const MASS_DIMS: Dims = Dims([0, 1, 0, 0, 0, 0]);

/// Shared lifetime admission for a family of immutable SIMP continuation stages.
/// Clones share counters and absolute ceilings, not numerical caches or model
/// meaning. Construct all stages with the SAME handle to avoid resetting work
/// when projection sharpness or penalization changes. This handle is deliberately
/// single-threaded, like the provider's RefCell cache; it is not a global byte cap.
#[derive(Debug, Clone)]
pub struct SimpNewtonResources(Rc<Resources>);

#[derive(Debug)]
struct Resources {
    budget: Cell<ComplianceNewtonBudget>,
    work: Cell<ComplianceNewtonWork>,
}
impl SimpNewtonResources {
    /// Start an empty work account. Zero ceilings are valid admission policies.
    pub fn new(budget: ComplianceNewtonBudget) -> Self {
        Self(Rc::new(Resources { budget: Cell::new(budget),
            work: Cell::new(ComplianceNewtonWork::default()) }))
    }
    pub fn work(&self) -> ComplianceNewtonWork { self.0.work.get() }
    pub fn budget(&self) -> ComplianceNewtonBudget { self.0.budget.get() }
    /// Change absolute allowances for ALL attached stages, never spent work.
    pub fn set_budget(&self, budget: ComplianceNewtonBudget) { self.0.budget.set(budget); }
}

struct Cached<'a> {
    bits: Vec<u64>,
    goal: SimpComplianceGoal<'a>,
    volume: PhysicsSample,
    mass: Option<PhysicsSample>,
}

/// Joint raw-density physics for compliance objectives / volume constraints, or
/// minimum-volume objectives / compliance constraints. Create the scalar views
/// with [`Self::compliance`], [`Self::volume`] and [`Self::mass`], then bind
/// them into the IR. Material mass uses the stiffness owner's source-bound
/// densities and geometric cell volumes, without normalizing kilograms away.
/// [`Self::new_robust`] selects eroded stiffness and dilated volume AND mass.
///
/// Preparation reserves one forward filter, each independent load, and two
/// transpose filters, plus one transpose for material-bound mass. Compliance products reserve two filters and each active
/// load tangent; nonlinear volume products reserve only two filters. Reservations
/// are conservative ceilings, not actual solve counts. All Krylov work, including
/// failures, is charged once to this owner rather than to either scalar view.
/// Cached reads and cached zero directions do not spend additional work.
///
/// A failed preparation never replaces the previous cache. A successful trial
/// may evict it; evaluating an older accepted tape re-prepares that exact point.
/// Loads, parameters, geometry and weights cannot change during the model's life.
/// Raw densities are in [0,1]; the native second-order goal additionally requires
/// stiffness projected densities strictly inside (1e-12,1). No endpoint Hessian,
/// certified feature size, process checkpoint, or continuum certificate is claimed.
pub struct SimpNewtonModel<'a> {
    filter: &'a DensityFilter,
    params: SimpParams,
    volume_eta: f64,
    elasticity: &'a DensityElasticity,
    loads: Vec<(Vec<f64>, f64)>,
    volumes: Vec<f64>,
    options: ComplianceSolveOptions,
    resources: SimpNewtonResources,
    cache: RefCell<Option<Cached<'a>>>,
    last_native_error: RefCell<Option<ComplianceError>>,
}

impl std::fmt::Debug for SimpNewtonModel<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SimpNewtonModel")
            .field("parameters", &self.elasticity.cells())
            .field("loads", &self.loads.len())
            .field("compliance_threshold", &self.params.eta)
            .field("volume_threshold", &self.volume_eta)
            .field("work", &self.resources.0.work.get())
            .field("budget", &self.resources.0.budget.get())
            .finish_non_exhaustive()
    }
}

#[derive(Debug, Clone, Copy)]
enum Output {
    Compliance,
    Volume,
    Mass,
    LoadCase(usize),
}

/// A named scalar view. Cloning shares, and cannot reset, its owner's state.
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
    /// Ordinary SIMP: admit immutable modeling data without numerical solves.
    /// `volumes` must be the physical design-cell volumes, in cell order.
    pub fn new(pipeline: &'a DesignPipeline, elasticity: &'a DensityElasticity,
        loads: &[LoadCase<'_>], volumes: &[f64], options: ComplianceSolveOptions,
        budget: ComplianceNewtonBudget,
    ) -> Result<Self, PhysicsError> {
        Self::from_parts(&pipeline.filter, pipeline.params, pipeline.params.eta,
            elasticity, loads, volumes, options, budget)
    }

    /// Robust three-field design using `eta +/- eta_offset` from the existing
    /// RobustPipeline. Compliance uses erosion; the volume view uses dilation.
    /// The nominal volume is not the constraint and no adaptive target is hidden
    /// here: the caller declares its dilated-volume allowance explicitly in IR.
    ///
    /// No extra Helmholtz or elasticity solves are needed compared with ordinary
    /// joint SIMP preparation. The projections differ only in their pointwise
    /// values/derivatives; all caches, work ceilings and continuation semantics
    /// remain those of this owner. Zero offset reduces to ordinary SIMP.
    pub fn new_robust(pipeline: &'a RobustPipeline, elasticity: &'a DensityElasticity,
        loads: &[LoadCase<'_>], volumes: &[f64], options: ComplianceSolveOptions,
        budget: ComplianceNewtonBudget,
    ) -> Result<Self, PhysicsError> {
        let (params, volume_eta) = robust_parameters(pipeline.params, pipeline.eta_offset)
            .map_err(|error| match error {
                ComplianceError::InvalidInput(message) => PhysicsError::InvalidInput(message),
                _ => PhysicsError::InvalidInput("invalid robust thresholds"),
            })?;
        Self::from_parts(&pipeline.filter, params, volume_eta, elasticity, loads, volumes, options, budget)
    }

    /// Ordinary SIMP stage with a caller-owned shared lifetime work account.
    /// Admission does no numerical work; the stage retains a separate cache.
    pub fn new_with_resources(pipeline: &'a DesignPipeline, elasticity: &'a DensityElasticity,
        loads: &[LoadCase<'_>], volumes: &[f64], options: ComplianceSolveOptions,
        resources: SimpNewtonResources,
    ) -> Result<Self, PhysicsError> {
        let mut model = Self::new(pipeline, elasticity, loads, volumes, options, resources.budget())?;
        model.resources = resources;
        Ok(model)
    }

    /// Robust stage with the SAME work account as other continuation stages.
    /// The threshold and SIMP parameters remain immutable per stage; sharing
    /// admission must never mean reusing another stage's primal or derivatives.
    pub fn new_robust_with_resources(pipeline: &'a RobustPipeline, elasticity: &'a DensityElasticity,
        loads: &[LoadCase<'_>], volumes: &[f64], options: ComplianceSolveOptions,
        resources: SimpNewtonResources,
    ) -> Result<Self, PhysicsError> {
        let mut model = Self::new_robust(pipeline, elasticity, loads, volumes, options, resources.budget())?;
        model.resources = resources;
        Ok(model)
    }

    /// Borrow this stage's account for construction of further immutable stages.
    /// Counters reported by every attached model are family-wide, not per-stage.
    pub fn resources(&self) -> SimpNewtonResources { self.resources.clone() }

    #[allow(clippy::too_many_arguments)]
    fn from_parts(
        filter: &'a DensityFilter,
        p: SimpParams,
        volume_eta: f64,
        elasticity: &'a DensityElasticity,
        loads: &[LoadCase<'_>],
        volumes: &[f64],
        options: ComplianceSolveOptions,
        budget: ComplianceNewtonBudget,
    ) -> Result<Self, PhysicsError> {
        let n = elasticity.cells();
        if n == 0 || filter.cells() != n || volumes.len() != n {
            return Err(PhysicsError::InvalidInput("SIMP mesh, filter and volume shapes must agree"));
        }
        if !(p.e_min.is_finite() && p.e_min > 0.0 && p.e_min < 1.0
            && p.penal.is_finite() && p.penal >= 1.0 && p.beta.is_finite() && p.beta >= 0.0
            && p.eta.is_finite() && (0.0..=1.0).contains(&p.eta)
            && volume_eta.is_finite() && (0.0..=1.0).contains(&volume_eta)) {
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
        if loads.is_empty() {
            return Err(PhysicsError::InvalidInput("at least one structural load case required"));
        }
        for load in loads {
            if load.force.len() != elasticity.n() || !(load.weight.is_finite() && load.weight >= 0.0) {
                return Err(PhysicsError::InvalidInput("invalid structural load shape or weight"));
            }
            for &force in load.force { finite(force, "structural force")?; }
        }
        loads
            .len()
            .checked_add(3 + usize::from(elasticity.material_field().is_some()))
            .ok_or(PhysicsError::InvalidInput("too many loads"))?;
        Ok(Self { filter, params: p, volume_eta, elasticity,
            loads: loads.iter().map(|l| (l.force.to_vec(), l.weight)).collect(),
            volumes: volumes.to_vec(), options, resources: SimpNewtonResources::new(budget),
            cache: RefCell::new(None),
            last_native_error: RefCell::new(None) })
    }

    pub fn compliance(&self, study: impl Into<String>) -> Result<SimpNewtonGoal<'_, 'a>, PhysicsError> {
        self.view(study.into(), Output::Compliance)
    }
    pub fn volume(&self, study: impl Into<String>) -> Result<SimpNewtonGoal<'_, 'a>, PhysicsError> {
        self.view(study.into(), Output::Volume)
    }
    /// Kilograms of projected material, for minimum-mass design or a physical
    /// mass limit. Shares the same primal state and lifetime work account as
    /// compliance and volume; requires an assigned material field.
    pub fn mass(&self, study: impl Into<String>) -> Result<SimpNewtonGoal<'_, 'a>, PhysicsError> {
        if self.elasticity.reference_cell_masses_kg().is_none() {
            return Err(PhysicsError::InvalidInput(
                "physical mass requires an assigned material field",
            ));
        }
        self.view(study.into(), Output::Mass)
    }
    /// Expose one independently solved structural load case as its own scalar
    /// compliance physics view. The original aggregate weight is not applied.
    /// All views share this model's primal cache and lifetime work account.
    pub fn load_case_compliance(&self, index: usize, study: impl Into<String>)
        -> Result<SimpNewtonGoal<'_, 'a>, PhysicsError>
    {
        if index >= self.loads.len() {
            return Err(PhysicsError::InvalidInput("load-case index out of range"));
        }
        self.view(study.into(), Output::LoadCase(index))
    }
    fn view(&self, study: String, output: Output) -> Result<SimpNewtonGoal<'_, 'a>, PhysicsError> {
        if study.is_empty() { return Err(PhysicsError::InvalidInput("nonempty study name required")); }
        Ok(SimpNewtonGoal { model: self, study, output })
    }
    /// Actual stiffness and physical-inventory (volume and mass) thresholds.
    pub fn projection_thresholds(&self) -> (f64, f64) { (self.params.eta, self.volume_eta) }
    pub fn work(&self) -> ComplianceNewtonWork { self.resources.0.work.get() }
    pub fn budget(&self) -> ComplianceNewtonBudget { self.resources.0.budget.get() }
    pub fn set_budget(&self, budget: ComplianceNewtonBudget) { self.resources.0.budget.set(budget); }
    /// A later success does not erase the most recent native refusal.
    pub fn last_native_error(&self) -> Option<ComplianceError> { self.last_native_error.borrow().clone() }

    fn native_error(&self, error: ComplianceError) -> PhysicsError {
        let result = match &error {
            ComplianceError::InvalidInput(message) => PhysicsError::InvalidInput(*message),
            ComplianceError::Stopped(EvaluationStop::Cancelled) => PhysicsError::Cancelled,
            ComplianceError::Stopped(EvaluationStop::TotalBudget { .. }) => PhysicsError::Budget {
                resource: "SIMP Krylov iterations", used: self.resources.0.work.get().linear_iterations,
                limit: self.resources.0.budget.get().linear_iterations,
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
        let mut work = self.resources.0.work.get(); let budget = self.resources.0.budget.get();
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
        self.resources.0.work.set(work);
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
        self.resources.0.work.set(work);
        let result = result.map_err(|e| self.native_error(e))?;
        poll(cx)?;
        Ok(result)
    }

    fn at_point<R>(
        &self,
        point: &[f64],
        cx: Option<&Cx<'_>>,
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
            let has_mass = self.elasticity.reference_cell_masses_kg().is_some();
            let prepared = self.controlled(loads.len() + 3 + usize::from(has_mass), cx, |control| {
                let goal = SimpComplianceGoal::prepare(self.filter, self.params, self.volume_eta,
                    self.elasticity, point, &loads, self.options, control)?;
                let (value, gradient) = goal.volume_value_gradient(&self.volumes, control)?;
                let mass = if has_mass {
                    let (value, gradient) = goal.mass_value_gradient(control)?;
                    Some(PhysicsSample { value, gradient })
                } else { None };
                Ok(Cached { bits: point.iter().map(|x| x.to_bits()).collect(), goal,
                    volume: PhysicsSample { value, gradient }, mass })
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
                Output::Compliance | Output::LoadCase(_) => COMPLIANCE_DIMS, Output::Volume => Dims::NONE,
                Output::Mass => MASS_DIMS,
            } }
    }
    fn value_gradient(
        &self,
        point: &[f64],
        cx: Option<&Cx<'_>>,
    ) -> Result<PhysicsSample, PhysicsError> {
        self.model.at_point(point, cx, |c| match self.output {
            Output::Compliance => Ok(PhysicsSample { value: c.goal.value(), gradient: c.goal.gradient().to_vec() }),
            Output::Volume => Ok(c.volume.clone()),
            Output::Mass => Ok(c.mass.as_ref().expect("admitted material mass view").clone()),
            Output::LoadCase(index) => self.model.controlled(1, cx, |control| {
                let (value, gradient) = c.goal.load_case_value_gradient(index, control)?;
                Ok(PhysicsSample { value, gradient })
            }),
        })
    }
    fn second_order(&self) -> Option<&dyn PhysicsHessian> { Some(self) }
}
impl PhysicsHessian for SimpNewtonGoal<'_, '_> {
    fn hessian_vector_product(
        &self,
        point: &[f64],
        direction: &[f64],
        cx: Option<&Cx<'_>>,
    ) -> Result<Vec<f64>, PhysicsError> {
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
                    if self.model.params.beta == 0.0 { return Ok(vec![0.0; direction.len()]); }
                    self.model.controlled(2, cx, |control|
                        c.goal.volume_hessian_vector_product(&self.model.volumes, direction, control))
                }
                Output::Mass => {
                    if self.model.params.beta == 0.0 { return Ok(vec![0.0; direction.len()]); }
                    self.model.controlled(2, cx, |control|
                        c.goal.mass_hessian_vector_product(direction, control))
                }
                Output::LoadCase(index) => self.model.controlled(3, cx, |control|
                    c.goal.load_case_hessian_vector_product(index, direction, control)),
            }
        })
    }
}
