//! Resolve physical operating points with the existing Newton--Krylov owner,
//! then admit the original complete resonance image. No second nonlinear solver,
//! relaxation-to-rest simulation, or accepted-state mutation is introduced.
use fs_phs::harmonic::equilibrium::{EquilibriumLinearization, LinearizationError};
use fs_phs::harmonic::equilibrium::stationary::{StationaryProblem, StationaryWorkspace};
use fs_solver::{NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport,
    NonlinearProblem, SolverRunProgress};
use std::cell::{Cell, RefCell};

/// Numerical-owner controls and hard call/storage ceilings. Source callbacks
/// have their own applicability and cost; these counts are not a wall-time cap.
#[derive(Clone, Copy, Debug)]
pub struct EquilibriumSolveOptions {
    /// Uses scaled original balances. Relative tolerance must be zero: a large
    /// initial imbalance cannot enlarge the final physical acceptance target.
    pub newton: NewtonKrylovConfig,
    pub relative_balance_tolerance: f64,
    pub max_iterations: usize,
    /// Includes the initial residual and rejected line-search trials.
    pub max_residual_calls: usize,
    /// Each action invokes at most one original storage Hessian and loss tangent.
    /// Final dense linearization has its separate source action budget.
    pub max_jacobian_calls: usize,
    /// Conservative live solver + physical scratch payload, not process RSS.
    pub max_workspace_bytes: usize,
}

#[derive(Debug)]
pub enum EquilibriumSolveError {
    Invalid(&'static str),
    Physical(LinearizationError),
    Newton(NewtonError),
    Budget { residual_calls: usize, jacobian_calls: usize },
    Stalled { report: NewtonReport, last_source_refusal: Option<LinearizationError> },
    Interrupted(SolverRunProgress),
    Cancelled,
}
impl core::fmt::Display for EquilibriumSolveError {
    fn fmt(&self,f:&mut core::fmt::Formatter<'_>)->core::fmt::Result {
        write!(f,"stationary resonance preparation refused: {self:?}")
    }
}
impl std::error::Error for EquilibriumSolveError {}

/// Validated original operating point and its complete linear resonance model.
/// No candidate state is installed into its source; callers retain the original
/// time/history owner and can inspect the stationary branch before any restart.
pub struct SolvedEquilibrium {
    pub linearization: EquilibriumLinearization,
    pub solve: NewtonReport,
    pub residual_calls: usize,
    pub jacobian_calls: usize,
}

struct Adapter<'p,'s,'c> {
    problem: &'p StationaryProblem<'s>,
    workspace: RefCell<StationaryWorkspace<'p,'s>>,
    cancel: RefCell<&'c mut dyn FnMut()->bool>,
    cancelled: Cell<bool>,
    exhausted: Cell<bool>,
    residuals: Cell<usize>,
    jacobians: Cell<usize>,
    options: EquilibriumSolveOptions,
    refusal: RefCell<Option<LinearizationError>>,
}
impl Adapter<'_, '_, '_> {
    fn poll(&self)->bool {
        if !self.cancelled.get() && (self.cancel.borrow_mut())() {self.cancelled.set(true);}
        self.cancelled.get() || self.exhausted.get()
    }
    fn charge(&self,counter:&Cell<usize>,limit:usize)->bool {
        if self.poll() {return false;}
        if counter.get()>=limit {self.exhausted.set(true);return false;}
        counter.set(counter.get()+1);true
    }
    fn interrupted(&self)->Option<EquilibriumSolveError> {
        if self.exhausted.get() {Some(EquilibriumSolveError::Budget {
            residual_calls:self.residuals.get(),jacobian_calls:self.jacobians.get()})}
        else if self.cancelled.get() {Some(EquilibriumSolveError::Cancelled)} else {None}
    }
}
impl NonlinearProblem for Adapter<'_, '_, '_> {
    fn dimension(&self)->usize {self.problem.dimension()}
    fn residual(&self,z:&[f64],out:&mut[f64]) {
        out.fill(f64::NAN);
        if !self.charge(&self.residuals,self.options.max_residual_calls) {return;}
        let result=self.problem.residual_into(z,&mut self.workspace.borrow_mut(),out,||self.poll());
        if let Err(error)=result {*self.refusal.borrow_mut()=Some(error);}
        // A domain refusal remains NaN, so the original line search may shorten
        // the trial. Call exhaustion and cancellation are latched separately.
    }
    fn jacobian_apply(&self,z:&[f64],d:&[f64],out:&mut[f64]) {
        out.fill(f64::NAN);
        if !self.charge(&self.jacobians,self.options.max_jacobian_calls) {return;}
        let result=self.problem.jacobian_action_into(z,d,&mut self.workspace.borrow_mut(),out,||self.poll());
        if let Err(error)=result {*self.refusal.borrow_mut()=Some(error);}
    }
}

/// Solve a held-input stationary balance in explicit scaled coordinates, using
/// fs-solver's unchanged globalization, forcing, Krylov and cancellation path.
/// All ordinary ODE/descriptor/material coordinates remain present. The final
/// original-equation and Hessian admissions are independent of solver telemetry.
/// A converged root may be unstable or one of several branches; no root selection
/// or uniqueness certificate is inferred. Domain callbacks must be pure.
pub fn solve_equilibrium(
    problem:&StationaryProblem<'_>,initial_state:&[f64],options:EquilibriumSolveOptions,
    mut cancelled:impl FnMut()->bool,
)->Result<SolvedEquilibrium,EquilibriumSolveError> {
    if cancelled() {return Err(EquilibriumSolveError::Cancelled);}
    options.newton.validate().map_err(EquilibriumSolveError::Newton)?;
    let n=problem.dimension();let k=options.newton.linear_restart;
    if options.newton.relative_tolerance!=0.0
        || !options.relative_balance_tolerance.is_finite()
        || options.relative_balance_tolerance<=0.0 || options.relative_balance_tolerance>=1.0
        || options.newton.absolute_tolerance>options.relative_balance_tolerance
        || options.max_iterations==0 || options.max_iterations>10_000
        || options.max_residual_calls==0 || options.max_jacobian_calls==0 || k>1024 {
        return Err(EquilibriumSolveError::Invalid("invalid stationary accuracy or bounded solver controls"));
    }
    let (nn,kk,it)=(n as u128,k as u128,options.max_iterations as u128);
    let bytes=128*nn*(kk+16)+32*(kk+4)*(kk+4)+1024*(it+1);
    if bytes>options.max_workspace_bytes as u128 {
        return Err(EquilibriumSolveError::Invalid("stationary Krylov/workspace payload exceeds capacity"));
    }
    let z=problem.scaled_state(initial_state).map_err(EquilibriumSolveError::Physical)?;
    let adapter=Adapter {problem,workspace:RefCell::new(problem.workspace()),
        cancel:RefCell::new(&mut cancelled),cancelled:Cell::new(false),exhausted:Cell::new(false),
        residuals:Cell::new(0),jacobians:Cell::new(0),options,refusal:RefCell::new(None)};
    let construction=NewtonKrylovState::new(&adapter,z,options.newton);
    if let Some(error)=adapter.interrupted() {return Err(error);}
    let mut solver=construction.map_err(|e| match adapter.refusal.borrow_mut().take() {
        Some(source)=>EquilibriumSolveError::Physical(source),None=>EquilibriumSolveError::Newton(e),
    })?;
    let outcome=solver.run_controlled(&adapter,options.max_iterations,&mut ||adapter.poll());
    if let Some(error)=adapter.interrupted() {return Err(error);}
    if outcome.progress!=SolverRunProgress::Complete {
        return Err(EquilibriumSolveError::Interrupted(outcome.progress));
    }
    if !outcome.report.converged {
        return Err(EquilibriumSolveError::Stalled {report:outcome.report,
            last_source_refusal:adapter.refusal.borrow_mut().take()});
    }
    let image=problem.linearize(&solver.x,options.relative_balance_tolerance,||adapter.poll());
    if let Some(error)=adapter.interrupted() {return Err(error);}
    let linearization=image.map_err(EquilibriumSolveError::Physical)?;
    if adapter.poll() {return Err(EquilibriumSolveError::Cancelled);}
    Ok(SolvedEquilibrium {linearization,solve:outcome.report,
        residual_calls:adapter.residuals.get(),jacobian_calls:adapter.jacobians.get()})
}
