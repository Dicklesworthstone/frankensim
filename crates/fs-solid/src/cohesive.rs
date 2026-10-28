//! Quasistatic cohesive interfaces coupled to a fixed reduced elastic operator.
//!
//! Solves `K u + sum B_i^T A_i t_i(B_i u + g_i; history_i) = f` through
//! the existing Newton/FGMRES owner. `u` and prescribed jumps `g` are metres,
//! `f` is newtons, `A` is original bonded area, and `B` is dimensionless.
//! A sparse jump map may represent a paired-node interface or a supplied trace
//! interpolation. Geometry, fixed local frames, quadrature/area, eliminated
//! constraints and the elastic operator are the caller's responsibility.
//!
//! All trial tractions use the SAME last accepted damage. Only a converged,
//! independently residual-checked endpoint publishes new displacement/history.
//! The material's selected piecewise tangent is used at switches. Softening
//! can be indefinite; FGMRES is intentional, not an SPD/CG approximation.
//!
//! This is equilibrium continuation, not a dynamic snap-back solver, stability
//! proof or unique-root selection. Subdivide load paths to resolve events and
//! reversals. Kinematics and bond topology stay fixed; no thermal contact,
//! healing, friction, automatic bond creation, or calibrated print prediction
//! is inferred. Reported energy is COHESIVE-ONLY straight-path work/storage/
//! fracture plus separately reported physical bond viscosity, not a complete
//! structural, actuator, thermal or entropy balance.

use std::cell::{Cell, RefCell};
use std::collections::BTreeSet;

use fs_exec::Cx;
use fs_material::cohesive::{BilinearCohesive, CohesiveEnergy, CohesiveError, CohesivePoint, CohesiveRateStep};
use fs_material::cohesive::source::ResolvedBilinearCohesive;
use fs_material::state_point::MaterialStatePointError;
use fs_matdb::MaterialStateId;
use fs_solver::{LinearOp, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, NonlinearProblem, SolverRunProgress};

/// Geometry-backed matching P1 interfaces and displacement transfer.
pub mod p1;

/// Thermal shrinkage and Maxwell history driving a restrained cohesive joint.
pub mod thermal_joint;

/// One sparse column of a local displacement-jump map.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct JumpTerm {
    pub dof: usize,
    /// Normal, shear-1 and shear-2 jump per metre of this displacement DOF.
    pub local: [f64; 3],
}

/// One originally bonded quadrature patch in a caller-supplied fixed local frame.
#[derive(Debug, Clone)]
pub struct CohesivePatch<'a> {
    id: u64,
    area_m2: f64,
    law: &'a BilinearCohesive,
    source: Option<&'a ResolvedBilinearCohesive>,
    terms: Vec<JumpTerm>,
}
impl<'a> CohesivePatch<'a> {
    /// Admit at most 64 distinct sparse displacement columns. Duplicate DOFs
    /// refuse rather than silently changing the caller's interpolation. An empty
    /// map is a fully prescribed interface: its history/energy still evolves.
    pub fn new(id: u64, area_m2: f64, law: &'a BilinearCohesive,
        mut terms: Vec<JumpTerm>) -> Result<Self, Error> {
        positive(area_m2)?;
        if terms.len() > 64 { return Err(Error::Invalid("patch admits at most 64 jump columns")); }
        for term in &terms {
            for value in term.local { finite(value)?; }
            if term.local == [0.0; 3] { return Err(Error::Invalid("zero jump column")); }
        }
        terms.sort_by_key(|term| term.dof);
        if terms.windows(2).any(|p| p[0].dof == p[1].dof) { return Err(Error::Invalid("duplicate patch DOF")); }
        Ok(Self { id, area_m2, law, source: None, terms })
    }
    /// Bind ordered interface data to the same work-conjugate mechanics owner.
    /// Actual receiving material states must match the source surface roles.
    /// Source conditions remain fixed for every damage history of this patch.
    pub fn from_interface(id: u64, area_m2: f64, source: &'a ResolvedBilinearCohesive,
        surface_materials: [&MaterialStateId; 2], terms: Vec<JumpTerm>) -> Result<Self, Error> {
        source.require_surface_materials(surface_materials[0], surface_materials[1])
            .map_err(|error| Error::Source { patch: id, error })?;
        let mut patch = Self::new(id, area_m2, source.law(), terms)?;
        patch.source = Some(source);
        Ok(patch)
    }
    /// Original source bundle when this patch was constructed from interface data.
    #[must_use]
    pub fn source(&self) -> Option<&'a ResolvedBilinearCohesive> { self.source }
    #[must_use]
    pub fn id(&self) -> u64 { self.id }
    #[must_use]
    pub fn area_m2(&self) -> f64 { self.area_m2 }

    fn jump(&self, u: &[f64], offset: [f64; 3]) -> Result<[f64; 3], Error> {
        let mut jump = offset;
        for term in &self.terms {
            for (value, weight) in jump.iter_mut().zip(term.local) { *value = finite(*value + weight*u[term.dof])?; }
        }
        Ok(jump)
    }
    fn scatter(&self, traction: [f64; 3], output: &mut [f64]) -> Result<(), Error> {
        for term in &self.terms {
            let force = self.area_m2*dot(term.local, traction);
            output[term.dof] = finite(output[term.dof]+force)?;
        }
        Ok(())
    }
}

/// Allocation/work admission. Does not account for storage internal to K.
#[derive(Debug, Clone, Copy)]
pub struct Budget {
    pub max_dofs: usize,
    pub max_patches: usize,
    pub max_jump_terms: usize,
    /// Maximum `2*(restart+1)*dofs + (restart+1)*restart` basis/Hessenberg
    /// scalar entries. Other solver vectors/history are separately bounded by
    /// admitted dimensions and iteration counts; this is not a total byte cap.
    pub max_krylov_basis_entries: usize,
}

/// Fixed geometry/operator binding. K must stay pure, linear and unchanged for
/// every state/solve derived from this object, including through interior mutation.
pub struct CohesiveSystem<'a> {
    elastic: &'a dyn LinearOp,
    patches: Vec<CohesivePatch<'a>>,
    n: usize,
    budget: Budget,
}
impl<'a> CohesiveSystem<'a> {
    /// Preserve input patch order for prescribed-jump arrays and observations.
    /// Callers must remove rigid modes or otherwise supply solvable restraints.
    pub fn new(elastic: &'a dyn LinearOp, patches: Vec<CohesivePatch<'a>>, budget: Budget)
        -> Result<Self, Error> {
        let n = elastic.n();
        if n == 0 || n > budget.max_dofs || patches.is_empty() || patches.len() > budget.max_patches {
            return Err(Error::Budget("degrees of freedom or interface patches"));
        }
        let mut terms = 0usize;
        let mut ids = BTreeSet::new();
        for patch in &patches {
            if !ids.insert(patch.id) { return Err(Error::Invalid("duplicate patch identity")); }
            terms = terms.checked_add(patch.terms.len()).ok_or(Error::Budget("jump map overflow"))?;
            if terms > budget.max_jump_terms { return Err(Error::Budget("jump columns")); }
            if patch.terms.iter().any(|term| term.dof >= n) { return Err(Error::Invalid("jump DOF outside elastic operator")); }
        }
        Ok(Self { elastic, patches, n, budget })
    }
    #[must_use]
    pub fn patches(&self) -> &[CohesivePatch<'a>] { &self.patches }
    /// Explicitly initialize bonded zero-jump history. This does not infer a
    /// bond from proximity, temperature or a thermal conductance coefficient.
    #[must_use]
    pub fn bonded_state(&self) -> CohesiveState<'_, 'a> {
        CohesiveState { system: self, displacements: vec![0.0; self.n],
            points: self.patches.iter().map(|patch| patch.law.bonded_point()).collect() }
    }
}

/// Immutable accepted displacement and damage, tied to the same system.
#[derive(Clone)]
pub struct CohesiveState<'m, 'a> {
    system: &'m CohesiveSystem<'a>,
    displacements: Vec<f64>,
    points: Vec<CohesivePoint<'a>>,
}

/// Explicit numerical and applicability controls; no constitutive defaults.
#[derive(Debug, Clone, Copy)]
pub struct Config {
    pub newton: NewtonKrylovConfig,
    pub max_newton_iterations: usize,
    /// Counts complete attempted residual/Jacobian actions, including the
    /// initial residual and independent final residual check.
    pub max_operator_evaluations: usize,
    /// Independent absolute infinity-norm equilibrium gate, N.
    pub force_tolerance_n: f64,
    pub local_energy_tolerance_j_m2: f64,
    pub total_cohesive_energy_tolerance_j: f64,
    /// Caller-declared reduced-coordinate range; not an inferred strain bound.
    pub maximum_displacement_m: f64,
}

#[derive(Debug, Clone)]
pub enum Error {
    Invalid(&'static str),
    Budget(&'static str),
    Material { patch: u64, error: CohesiveError },
    Source { patch: u64, error: MaterialStatePointError },
    Cancelled,
    Solver(String),
    NotConverged(NewtonReport),
    ForceResidual(f64),
    CohesiveEnergyDefect(f64),
}
impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result { write!(f, "cohesive equilibrium refused: {self:?}") }
}
impl std::error::Error for Error {}

/// Successful endpoint. Publish `state` only after this complete result exists.
pub struct EquilibriumStep<'m, 'a> {
    pub state: CohesiveState<'m, 'a>,
    pub nonlinear: NewtonReport,
    pub maximum_force_residual_n: f64,
    pub operator_evaluations: usize,
    /// Area-integrated material energy, in joules (not the full solid balance).
    pub cohesive_energy_j: CohesiveEnergy,
    /// Physical bond-viscosity loss, J, distinct from intrinsic fracture and
    /// solver damping. Total cohesive path work is energy.work plus this loss.
    pub viscous_dissipated_j: f64,
}

impl<'m, 'a> CohesiveState<'m, 'a> {
    #[must_use]
    pub fn displacements_m(&self) -> &[f64] { &self.displacements }
    /// Same ordering as system.patches(); inspect jump, traction, damage and
    /// failed status before an explicit downstream contact/thermal update.
    #[must_use]
    pub fn points(&self) -> &[CohesivePoint<'a>] { &self.points }

    pub fn solve(&self, force_n: &[f64], prescribed_jumps_m: &[[f64; 3]], config: Config, cx: &Cx<'_>)
        -> Result<EquilibriumStep<'m, 'a>, Error> {
        self.solve_controlled(force_n, prescribed_jumps_m, config, &mut || cx.checkpoint().is_err())
    }

    /// Predicate form of solve with identical physics and cancellation points.
    /// Polls at least every 256 vector entries, every patch, around each elastic
    /// action and inside the native Newton/FGMRES driver. K itself must have a
    /// bounded apply latency; this wrapper cannot preempt an opaque callback.
    /// Cancellation, out-of-domain trials and exhausted work all refuse the
    /// whole load increment, leaving accepted physical history unchanged.
    pub fn solve_controlled<C: FnMut() -> bool>(&self, force_n: &[f64], prescribed_jumps_m: &[[f64; 3]],
        config: Config, cancelled: &mut C) -> Result<EquilibriumStep<'m, 'a>, Error> {
        self.solve_interval(force_n, prescribed_jumps_m, None, config, cancelled)
    }

    /// Native implicit equilibrium with physical bond rates over an explicit dt.
    /// Traction and tangent use the same accepted history and endpoint interval.
    /// The jump path is straight at constant rate; subdivide actual reversals.
    pub fn solve_over_time(&self, force_n: &[f64], prescribed_jumps_m: &[[f64; 3]],
        dt_s: f64, config: Config, cx: &Cx<'_>) -> Result<EquilibriumStep<'m, 'a>, Error> {
        self.solve_over_time_controlled(force_n, prescribed_jumps_m, dt_s, config,
            &mut || cx.checkpoint().is_err())
    }

    /// Predicate form with the same atomicity and native cancellation points.
    pub fn solve_over_time_controlled<C: FnMut() -> bool>(&self, force_n: &[f64],
        prescribed_jumps_m: &[[f64; 3]], dt_s: f64, config: Config,
        cancelled: &mut C) -> Result<EquilibriumStep<'m, 'a>, Error> {
        positive(dt_s)?;
        self.solve_interval(force_n, prescribed_jumps_m, Some(dt_s), config, cancelled)
    }

    fn solve_interval<C: FnMut() -> bool>(&self, force_n: &[f64],
        prescribed_jumps_m: &[[f64; 3]], dt_s: Option<f64>, config: Config,
        cancelled: &mut C) -> Result<EquilibriumStep<'m, 'a>, Error> {
        if cancelled() { return Err(Error::Cancelled); }
        let system = self.system;
        config.newton.validate().map_err(|e| Error::Solver(e.to_string()))?;
        for value in [config.force_tolerance_n, config.local_energy_tolerance_j_m2,
            config.total_cohesive_energy_tolerance_j, config.maximum_displacement_m] { positive(value)?; }
        if config.max_newton_iterations == 0 || config.max_operator_evaluations == 0 {
            return Err(Error::Budget("nonzero solver work allowances required"));
        }
        let r = config.newton.linear_restart;
        let entries = r.checked_add(1).and_then(|rp| system.n.checked_mul(2)
            .and_then(|n2| n2.checked_add(r)).and_then(|width| rp.checked_mul(width)))
            .ok_or(Error::Budget("Krylov basis overflow"))?;
        if entries > system.budget.max_krylov_basis_entries { return Err(Error::Budget("Krylov basis")); }
        if force_n.len() != system.n || prescribed_jumps_m.len() != system.patches.len()
            || system.elastic.n() != system.n { return Err(Error::Invalid("operator/load/patch dimensions differ")); }
        for chunk in force_n.chunks(256) {
            if cancelled() { return Err(Error::Cancelled); }
            for &value in chunk { finite(value)?; }
        }
        for offset in prescribed_jumps_m {
            if cancelled() { return Err(Error::Cancelled); }
            for &value in offset { finite(value)?; }
        }
        let problem = Equilibrium { state: self, force: force_n, offsets: prescribed_jumps_m, dt_s, config,
            cancelled: RefCell::new(cancelled), fault: RefCell::new(None), evaluations: Cell::new(0) };
        let solver = NewtonKrylovState::new(&problem, self.displacements.clone(), config.newton);
        problem.poll()?;
        let mut solver = solver.map_err(|e| Error::Solver(e.to_string()))?;
        let run = solver.run_controlled(&problem, config.max_newton_iterations, &mut || problem.poll().is_err());
        problem.poll()?;
        if run.progress != SolverRunProgress::Complete { return Err(Error::Solver(format!("{:?}", run.progress))); }
        if !run.report.converged { return Err(Error::NotConverged(run.report)); }
        let mut residual = vec![0.0; system.n];
        problem.action(&solver.x, None, &mut residual)?;
        let maximum_force_residual_n = residual.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
        if maximum_force_residual_n > config.force_tolerance_n { return Err(Error::ForceResidual(maximum_force_residual_n)); }
        let mut points = Vec::with_capacity(system.patches.len());
        let mut energy = CohesiveEnergy { work: 0.0, stored_change: 0.0, dissipated: 0.0, defect: 0.0 };
        let mut viscous_dissipated_j = 0.0;
        for (index, patch) in system.patches.iter().enumerate() {
            problem.poll()?;
            let step = problem.patch_trial(index, &solver.x)?;
            energy.work = finite(energy.work+patch.area_m2*step.equilibrium.energy.work)?;
            energy.stored_change = finite(energy.stored_change+patch.area_m2*step.equilibrium.energy.stored_change)?;
            energy.dissipated = finite(energy.dissipated+patch.area_m2*step.equilibrium.energy.dissipated)?;
            viscous_dissipated_j = finite(viscous_dissipated_j+patch.area_m2*step.viscous_dissipated_j_m2)?;
            points.push(step.equilibrium.point);
        }
        finite(energy.work+viscous_dissipated_j)?;
        energy.defect = finite(energy.work-energy.stored_change-energy.dissipated)?;
        if energy.defect.abs() > config.total_cohesive_energy_tolerance_j { return Err(Error::CohesiveEnergyDefect(energy.defect)); }
        problem.poll()?;
        Ok(EquilibriumStep {
            state: Self { system, displacements: solver.x, points }, nonlinear: run.report,
            maximum_force_residual_n, operator_evaluations: problem.evaluations.get(), cohesive_energy_j: energy, viscous_dissipated_j,
        })
    }
}

struct Equilibrium<'p, 'm, 'a> {
    state: &'p CohesiveState<'m, 'a>,
    force: &'p [f64],
    offsets: &'p [[f64; 3]],
    dt_s: Option<f64>,
    config: Config,
    cancelled: RefCell<&'p mut dyn FnMut() -> bool>,
    fault: RefCell<Option<Error>>,
    evaluations: Cell<usize>,
}
impl<'p, 'm, 'a> Equilibrium<'p, 'm, 'a> {
    fn fail(&self, error: Error) {
        if self.fault.borrow().is_none() { *self.fault.borrow_mut() = Some(error); }
    }
    fn poll(&self) -> Result<(), Error> {
        if let Some(error) = self.fault.borrow().as_ref() { return Err(error.clone()); }
        if (self.cancelled.borrow_mut())() {
            self.fail(Error::Cancelled);
            return Err(Error::Cancelled);
        }
        Ok(())
    }
    fn patch_trial(&self, index: usize, x: &[f64]) -> Result<CohesiveRateStep<'a>, Error> {
        let patch = &self.state.system.patches[index];
        let jump = patch.jump(x, self.offsets[index])?;
        match self.dt_s {
            Some(dt) => self.state.points[index].trial_over_time(jump, dt, self.config.local_energy_tolerance_j_m2),
            None => self.state.points[index].trial(jump, self.config.local_energy_tolerance_j_m2).map(Into::into),
        }.map_err(|error| Error::Material { patch: patch.id, error })
    }
    fn action(&self, x: &[f64], direction: Option<&[f64]>, output: &mut [f64]) -> Result<(), Error> {
        self.poll()?;
        if self.evaluations.get() >= self.config.max_operator_evaluations { return Err(Error::Budget("residual/Jacobian evaluations")); }
        self.evaluations.set(self.evaluations.get()+1);
        let system = self.state.system;
        if x.len() != system.n || output.len() != system.n || direction.is_some_and(|v| v.len() != system.n) {
            return Err(Error::Invalid("callback dimensions"));
        }
        for chunk in x.chunks(256) {
            self.poll()?;
            for &value in chunk {
                finite(value)?;
                if value.abs() > self.config.maximum_displacement_m { return Err(Error::Invalid("displacement applicability limit")); }
            }
        }
        let input = direction.unwrap_or(x);
        for chunk in input.chunks(256) { self.poll()?; for &v in chunk { finite(v)?; } }
        output.fill(f64::NAN);
        system.elastic.apply(input, output);
        self.poll()?;
        for chunk in output.chunks(256) { self.poll()?; for &v in chunk { finite(v)?; } }
        for (index, patch) in system.patches.iter().enumerate() {
            self.poll()?;
            let step = self.patch_trial(index, x)?;
            let traction = match direction {
                None => step.traction_pa,
                Some(v) => {
                    let jump_direction = patch.jump(v, [0.0; 3])?;
                    step.branch_tangent_pa_m.map(|row| dot(row, jump_direction))
                }
            };
            patch.scatter(traction, output)?;
        }
        if direction.is_none() {
            for (block, chunk) in output.chunks_mut(256).enumerate() {
                self.poll()?;
                for (i, value) in chunk.iter_mut().enumerate() { *value = finite(*value-self.force[block*256+i])?; }
            }
        }
        self.poll()
    }
}
impl NonlinearProblem for Equilibrium<'_, '_, '_> {
    fn dimension(&self) -> usize { self.state.system.n }
    fn residual(&self, x: &[f64], residual: &mut [f64]) {
        if let Err(error) = self.action(x, None, residual) { self.fail(error); residual.fill(f64::NAN); }
    }
    fn jacobian_apply(&self, x: &[f64], direction: &[f64], output: &mut [f64]) {
        if let Err(error) = self.action(x, Some(direction), output) { self.fail(error); output.fill(f64::NAN); }
    }
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 { a[0]*b[0]+a[1]*b[1]+a[2]*b[2] }
fn finite(value: f64) -> Result<f64, Error> {
    if value.is_finite() { Ok(value) } else { Err(Error::Invalid("nonfinite data/arithmetic")) }
}
fn positive(value: f64) -> Result<(), Error> {
    if value.is_finite() && value > 0.0 { Ok(()) } else { Err(Error::Invalid("positive finite data required")) }
}

#[cfg(test)]
mod tests;
