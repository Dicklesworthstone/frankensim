//! Fallible physical adapter for the existing Newton--Krylov owner.
use super::*;
use fs_solid::membrane::MembraneSolveSettings;
use fs_solver::{NewtonKrylovState, NonlinearProblem};
use std::cell::{Cell, RefCell};

/// Original material and membrane solver controls; neither tolerance is inferred.
#[derive(Clone, Copy, Debug)]
pub struct GasMembraneSolve {
    /// First-parcel entropy inversion, including its initial temperature.
    /// A sealed exterior retains its own independent source controls.
    pub gas: GasEntropySolve,
    /// Original force/displacement Newton controls and work/allocation bounds.
    pub membrane: MembraneSolveSettings,
}

/// Independently re-evaluated coupled equilibrium and original solver report.
#[derive(Clone, Debug)]
pub struct GasMembraneSolution {
    /// Fresh geometry, gas state, full force residual and conservative potential.
    pub equilibrium: GasMembraneEvaluation,
    /// Original nonlinear convergence history.
    pub newton: NewtonReport,
    /// Residual/tangent callbacks plus initial and fresh final physical checks.
    pub evaluations: usize,
    /// Coupled geometry proposals, each reconstructing one or two source parcels.
    /// Repeated Krylov actions reuse both gases at the unchanged trial shape.
    pub gas_trials: usize,
}

struct Adapter<'a, 'p, C> {
    model: &'a GasMembrane<'p>,
    loads: &'a [[f64; 3]],
    controls: GasEntropySolve,
    cancelled: RefCell<C>,
    stopped: RefCell<Option<GasMembraneError>>,
    cache: RefCell<Option<Trial>>,
    evaluations: Cell<usize>,
    gas_trials: Cell<usize>,
    maximum: usize,
}
impl<C: FnMut() -> bool> Adapter<'_, '_, C> {
    fn cancel(&self) -> bool {
        if self.stopped.borrow().is_some() { return true; }
        if (self.cancelled.borrow_mut())() {
            *self.stopped.borrow_mut() = Some(GasMembraneError::Cancelled);
            return true;
        }
        false
    }
    fn begin(&self) -> bool {
        if self.cancel() { return false; }
        if self.evaluations.get() >= self.maximum {
            *self.stopped.borrow_mut() = Some(GasMembraneError::EvaluationBudget);
            return false;
        }
        self.evaluations.set(self.evaluations.get() + 1);
        true
    }
    fn fresh(&self, displacement: &[[f64; 3]]) -> Result<Trial, GasMembraneError> {
        self.gas_trials.set(self.gas_trials.get() + 1);
        self.model.trial(displacement, self.loads, self.controls, &mut || self.cancel())
    }
    fn prepare(&self, x: &[f64]) -> bool {
        if !self.begin() || x.len() != self.model.free_dofs.len() || x.iter().any(|v| !v.is_finite()) {
            return false;
        }
        let matches = self.cache.borrow().as_ref().is_some_and(|trial| {
            self.model.free_dofs.iter().zip(x).all(|(&i, &value)| {
                trial.value.displacement_m[i / 3][i % 3] == value
            })
        });
        if !matches {
            match self.fresh(&self.model.expand(x, false)) {
                Ok(trial) => *self.cache.borrow_mut() = Some(trial),
                // The native line search rejects a nonfinite residual and can
                // backtrack into the admitted geometry/material domain. Source
                // refusals never become fabricated forces or accepted states.
                Err(GasMembraneError::Cancelled) => {
                    *self.stopped.borrow_mut() = Some(GasMembraneError::Cancelled);
                    return false;
                }
                Err(_) => return false,
            }
        }
        true
    }
}
impl<C: FnMut() -> bool> NonlinearProblem for Adapter<'_, '_, C> {
    fn dimension(&self) -> usize { self.model.free_dofs.len() }
    fn residual(&self, x: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if output.len() != self.dimension() || !self.prepare(x) { return; }
        let cache = self.cache.borrow();
        let value = &cache.as_ref().expect("prepared trial").value;
        for (out, &i) in output.iter_mut().zip(&self.model.free_dofs) {
            *out = value.residual_n[i / 3][i % 3];
        }
    }
    fn jacobian_apply(&self, x: &[f64], direction: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if output.len() != self.dimension() || !self.prepare(x) { return; }
        let cache = self.cache.borrow();
        if let Ok(action) = self.model.tangent_at(
            cache.as_ref().expect("prepared trial"), direction, &mut || self.cancel(),
        ) {
            output.copy_from_slice(&action);
        }
    }
}

impl GasMembrane<'_> {
    /// Solve the original material force balance and original gas isentrope
    /// together. This immutable proposal has no partial publication, pressure
    /// lag, separate gas iteration or replacement membrane constitutive law.
    pub fn solve(
        &self,
        dead_load_n: &[[f64; 3]],
        initial_displacement_m: &[[f64; 3]],
        settings: GasMembraneSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasMembraneSolution, GasMembraneError> {
        poll(&mut cancelled)?;
        settings.membrane.newton.validate().map_err(GasMembraneError::Newton)?;
        if settings.membrane.max_iterations == 0 || settings.membrane.max_evaluations == 0
            || settings.membrane.newton.linear_restart > settings.membrane.max_krylov_restart
        {
            return Err(GasMembraneError::Invalid("positive solve budgets and admitted Krylov restart required"));
        }
        self.admit_shape(initial_displacement_m)?;
        let adapter = Adapter {
            model: self, loads: dead_load_n, controls: settings.gas,
            cancelled: RefCell::new(cancelled), stopped: RefCell::new(None),
            cache: RefCell::new(None), evaluations: Cell::new(0), gas_trials: Cell::new(0),
            maximum: settings.membrane.max_evaluations,
        };
        if !adapter.begin() { return Err(adapter.stopped.borrow_mut().take().expect("stopped adapter")); }
        let initial_trial = adapter.fresh(initial_displacement_m)?;
        *adapter.cache.borrow_mut() = Some(initial_trial);
        let x = self.free_dofs.iter().map(|&i| initial_displacement_m[i / 3][i % 3]).collect();
        let initial = NewtonKrylovState::new(&adapter, x, settings.membrane.newton);
        if let Some(error) = adapter.stopped.borrow_mut().take() { return Err(error); }
        let mut state = initial.map_err(GasMembraneError::Newton)?;
        let initial_norm = state.residual_norm();
        let result = state.run_controlled(&adapter, settings.membrane.max_iterations, &mut || adapter.cancel());
        if let Some(error) = adapter.stopped.borrow_mut().take() { return Err(error); }
        if result.progress != SolverRunProgress::Complete {
            return Err(GasMembraneError::Interrupted(result.progress));
        }
        if !result.report.converged { return Err(GasMembraneError::Unconverged(result.report)); }
        if !adapter.begin() { return Err(adapter.stopped.borrow_mut().take().expect("stopped adapter")); }
        let final_trial = adapter.fresh(&self.expand(&state.x, false))?;
        let residual = self.free_dofs.iter().fold(0.0_f64, |norm, &i| {
            norm.hypot(final_trial.value.residual_n[i / 3][i % 3])
        });
        let target = settings.membrane.newton.absolute_tolerance
            .max(settings.membrane.newton.relative_tolerance * initial_norm);
        if !residual.is_finite() || residual > target {
            return Err(GasMembraneError::Invalid("fresh coupled physical residual missed the original Newton tolerance"));
        }
        if adapter.cancel() { return Err(adapter.stopped.borrow_mut().take().expect("stopped adapter")); }
        Ok(GasMembraneSolution {
            equilibrium: final_trial.value, newton: result.report,
            evaluations: adapter.evaluations.get(), gas_trials: adapter.gas_trials.get(),
        })
    }
}
