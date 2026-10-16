//! Boundary elimination and fallible physics adapter for the original solver.
use super::{MembraneError, MembranePatch, finite_vectors, poll};
use fs_solver::nonlinear::{NewtonKrylovConfig, NewtonKrylovState, NewtonReport, NonlinearProblem, SolverRunProgress};
use std::cell::{Cell, RefCell};

/// Explicit nonlinear accuracy/work controls; no stiffness regularization.
#[derive(Clone, Copy, Debug)]
pub struct MembraneSolveSettings {
    /// Original nonlinear solver configuration, in displacement/force SI units.
    pub newton: NewtonKrylovConfig,
    /// Maximum outer Newton attempts in this solve.
    pub max_iterations: usize,
    /// Maximum combined residual and tangent evaluations, including rejected
    /// trials, the initial state and the independent final residual evaluation.
    pub max_evaluations: usize,
    /// Allocation cap on the requested FGMRES restart columns.
    pub max_krylov_restart: usize,
}

/// A freshly re-evaluated equilibrium, not a stability or mesh-convergence claim.
#[derive(Clone, Debug)]
pub struct MembraneSolution {
    /// Three Cartesian displacements per original reference node, in metres.
    pub displacement_m: Vec<[f64; 3]>,
    /// Internal gradient minus dead load. At prescribed DOFs this is the
    /// applied support reaction; free entries meet the declared residual gate.
    pub residual_n: Vec<[f64; 3]>,
    /// Integrated reference-area stored energy, in joules.
    pub strain_energy_j: f64,
    /// Stored energy minus dead-load work, in joules.
    pub potential_energy_j: f64,
    /// Original Newton report, with the actual convergence result/history.
    pub newton: NewtonReport,
    /// Actual residual-plus-tangent evaluations consumed by this solve.
    pub evaluations: usize,
}

struct Adapter<'a, C> {
    patch: &'a MembranePatch,
    loads: &'a [[f64; 3]],
    fixed: Vec<Option<f64>>,
    free: Vec<usize>,
    cancelled: RefCell<&'a mut C>,
    stopped: RefCell<Option<MembraneError>>,
    evaluations: Cell<usize>,
    maximum: usize,
}
impl<C: FnMut() -> bool> Adapter<'_, C> {
    fn cancel(&self) -> bool {
        if self.stopped.borrow().is_some() { return true; }
        if (self.cancelled.borrow_mut())() {
            *self.stopped.borrow_mut() = Some(MembraneError::Cancelled);
            return true;
        }
        false
    }
    fn begin(&self) -> bool {
        if self.cancel() { return false; }
        if self.evaluations.get() >= self.maximum {
            *self.stopped.borrow_mut() = Some(MembraneError::Budget("residual/tangent evaluations"));
            return false;
        }
        self.evaluations.set(self.evaluations.get()+1);
        true
    }
    fn expand(&self, x: &[f64], direction: bool) -> Vec<[f64; 3]> {
        let mut u = vec![[0.0;3];self.patch.mesh.nodes.len()];
        if !direction {
            for (dof,value) in self.fixed.iter().enumerate() {
                if let Some(value) = value { u[dof/3][dof%3] = *value; }
            }
        }
        for (&dof,&value) in self.free.iter().zip(x) { u[dof/3][dof%3] = value; }
        u
    }
}
impl<C: FnMut() -> bool> NonlinearProblem for Adapter<'_, C> {
    fn dimension(&self) -> usize { self.free.len() }
    fn residual(&self, x: &[f64], output: &mut [f64]) {
        // The original driver rejects nonfinite trial residuals and backtracks.
        // No state/source refusal is replaced with a fabricated physical force.
        output.fill(f64::NAN);
        if x.len()!=self.free.len() || output.len()!=self.free.len() || !self.begin() { return; }
        let u = self.expand(x,false);
        if let Ok(value) = self.patch.evaluate(&u,||self.cancel()) {
            for (out,&dof) in output.iter_mut().zip(&self.free) {
                *out = value.internal_gradient_n[dof/3][dof%3]-self.loads[dof/3][dof%3];
            }
        }
    }
    fn jacobian_apply(&self, x: &[f64], direction: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if x.len()!=self.free.len() || direction.len()!=self.free.len()
            || output.len()!=self.free.len() || !self.begin() { return; }
        let u = self.expand(x,false);
        let v = self.expand(direction,true);
        if let Ok(value) = self.patch.tangent_apply(&u,&v,||self.cancel()) {
            for (out,&dof) in output.iter_mut().zip(&self.free) { *out = value[dof/3][dof%3]; }
        }
    }
}

impl MembranePatch {
    /// Solve static force balance through the original matrix-free Newton--Krylov
    /// owner. Prescribed entries are `(3*node+axis, displacement_m)` and are
    /// eliminated, not penalized. Their values override the initial guess.
    /// Loads are Cartesian dead forces [N], NOT follower pressure. The caller
    /// supplies any prestretch through boundary displacement, never through an
    /// invented membrane stiffness. Initial guesses are checked before solving.
    ///
    /// Cancellation, exhausted work, or failure returns no partially converged
    /// solution. This convenience call does not expose a resumable checkpoint;
    /// it does reuse the shared driver's globalization and staged iteration.
    pub fn solve(
        &self, dead_load_n: &[[f64;3]], prescribed: &[(usize,f64)],
        initial_displacement_m: &[[f64;3]], settings: MembraneSolveSettings,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MembraneSolution, MembraneError> {
        poll(&mut cancelled)?;
        settings.newton.validate().map_err(MembraneError::Newton)?;
        let n = self.mesh.nodes.len();
        if !finite_vectors(dead_load_n,n) || !finite_vectors(initial_displacement_m,n)
            || settings.max_iterations == 0 || settings.max_evaluations == 0 {
            return Err(MembraneError::Invalid("finite per-node loads/initial state and positive solve budgets required"));
        }
        if settings.newton.linear_restart > settings.max_krylov_restart {
            return Err(MembraneError::Budget("FGMRES restart allocation"));
        }
        let mut fixed = vec![None;3*n];
        for &(dof,value) in prescribed {
            poll(&mut cancelled)?;
            if dof >= fixed.len() || !value.is_finite() || fixed[dof].is_some() {
                return Err(MembraneError::Invalid("unique in-range finite prescribed displacements required"));
            }
            fixed[dof] = Some(value);
        }
        let free: Vec<_> = (0..3*n).filter(|&dof|fixed[dof].is_none()).collect();
        let x: Vec<_> = free.iter().map(|&dof|initial_displacement_m[dof/3][dof%3]).collect();
        let adapter = Adapter { patch:self, loads:dead_load_n, fixed, free,
            cancelled:RefCell::new(&mut cancelled), stopped:RefCell::new(None),
            evaluations:Cell::new(0), maximum:settings.max_evaluations };
        // Attribute invalid initial geometry to its physical element, rather
        // than returning only the generic solver's nonfinite-residual error.
        if !adapter.begin() { return Err(adapter.stopped.borrow_mut().take().expect("stopped adapter")); }
        self.evaluate(&adapter.expand(&x,false),||adapter.cancel())?;
        let initial = NewtonKrylovState::new(&adapter,x,settings.newton);
        if let Some(error) = adapter.stopped.borrow_mut().take() { return Err(error); }
        let mut state = initial.map_err(MembraneError::Newton)?;
        let initial_norm = state.residual_norm();
        let result = state.run_controlled(&adapter,settings.max_iterations,&mut ||adapter.cancel());
        if let Some(error) = adapter.stopped.borrow_mut().take() { return Err(error); }
        if result.progress != SolverRunProgress::Complete {
            return Err(MembraneError::Invalid("interrupted shared nonlinear callback; no equilibrium published"));
        }
        if !result.report.converged { return Err(MembraneError::Unconverged(result.report)); }
        if !adapter.begin() { return Err(adapter.stopped.borrow_mut().take().expect("stopped adapter")); }
        let displacement_m = adapter.expand(&state.x,false);
        let value = self.evaluate(&displacement_m,||adapter.cancel())?;
        let residual_n: Vec<[f64;3]> = value.internal_gradient_n.iter().zip(dead_load_n)
            .map(|(g,f)|std::array::from_fn(|i|g[i]-f[i])).collect();
        let free_norm = adapter.free.iter().fold(0.0_f64,|r,&dof|r.hypot(residual_n[dof/3][dof%3]));
        let target = settings.newton.absolute_tolerance.max(settings.newton.relative_tolerance*initial_norm);
        let work: f64 = dead_load_n.iter().zip(&displacement_m)
            .map(|(f,u)|f.iter().zip(u).map(|(f,u)|f*u).sum::<f64>()).sum();
        let potential_energy_j = value.strain_energy_j-work;
        if !finite_vectors(&residual_n,n) || !free_norm.is_finite() || free_norm > target
            || !potential_energy_j.is_finite() {
            return Err(MembraneError::Invalid("fresh physical residual or energy did not satisfy the solve contract"));
        }
        if adapter.cancel() { return Err(adapter.stopped.borrow_mut().take().expect("stopped adapter")); }
        Ok(MembraneSolution { displacement_m, residual_n, strain_energy_j:value.strain_energy_j,
            potential_energy_j, newton:result.report, evaluations:adapter.evaluations.get() })
    }
}