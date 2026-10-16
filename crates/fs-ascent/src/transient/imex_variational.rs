//! Weak-constraint reconstruction through the production operator IMEX adjoint.
//!
//! Autonomous u' = L(p)u + N(u,p) uses the existing ARS(2,2,2) stage solves and
//! checkpointed transpose sweep. The variational objective and L-BFGS engine
//! are shared with RK45 reconstruction. This module never differentiates Krylov
//! iterations, approximates a transpose, or implements another time integrator.

use crate::transient::variational::WindowError;
use crate::transient::variational::intervals::{IntervalScheme, IntervalTape};
use fs_solver::FlexiblePreconditioner;
use fs_time::adaptive::adjoint::trajectory::TrajectoryGradient;
use fs_time::stiff::{ImexSolveConfig, ImexSolveError, OperatorImex2};
use fs_time::stiff::adjoint::{ImexAdjointError, ImexVjp};
use fs_time::stiff::adjoint::trajectory::{
    ImexRecordingConfig, ImexRecordingStatus, ImexReplayBudget, ImexTrajectoryError, RecordedImex2,
};

/// Fixed method and per-interval resource allowances. Numerical accuracy is
/// controlled by the chosen step and primal/adjoint residuals, not by a new
/// adaptive error estimator. Model-error sigma still describes endpoint-state
/// uncertainty; neither the step nor the number of substeps rescales it.
#[derive(Debug, Clone, Copy)]
pub struct ImexWindowConfig {
    pub step: f64,
    pub solve: ImexSolveConfig,
    pub max_workspace_components: usize,
    pub max_forward_steps: usize,
    pub max_records: usize,
    pub replay: ImexReplayBudget,
    /// Sum of all requested interval step counts; bounds clock construction.
    pub max_clock_steps: usize,
    /// Owned time/count entries (2*intervals+1), excluding producer memory.
    pub max_grid_components: usize,
}

/// Immutable fixed-step clock, with separate primal and adjoint preconditioners.
/// Use `times()` to construct the WeakConstraintWindow. Times are obtained by
/// the SAME repeated binary64 addition as the production method, not a rounded
/// `(end-start)/h` step count. Observations must target these exact endpoints.
/// An autonomous model and both preconditioners must stay unchanged on replay.
pub struct ImexWindowPolicy<'a, P, Q> {
    times: Vec<f64>,
    counts: Vec<usize>,
    config: ImexWindowConfig,
    primal: &'a P,
    adjoint: &'a Q,
}
impl<P, Q> Clone for ImexWindowPolicy<'_, P, Q> {
    fn clone(&self) -> Self {
        Self { times:self.times.clone(), counts:self.counts.clone(), config:self.config,
            primal:self.primal, adjoint:self.adjoint }
    }
}
impl<'a, P, Q> ImexWindowPolicy<'a, P, Q> {
    pub fn new(start: f64, counts: &[usize], config: ImexWindowConfig,
        primal: &'a P, adjoint: &'a Q, cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, WindowError> {
        poll(cancelled)?;
        let required=counts.len().checked_mul(2).and_then(|n|n.checked_add(1))
            .ok_or(WindowError::Invalid("IMEX grid extent overflow"))?;
        if required>config.max_grid_components {
            return Err(WindowError::WorkspaceLimit{required,limit:config.max_grid_components});
        }
        if counts.is_empty() || !start.is_finite() || !config.step.is_finite() || config.step<=0.0
            || !config.solve.tolerance.is_finite() || config.solve.tolerance<=0.0 || config.solve.tolerance>=1.0
            || config.solve.restart==0 || config.solve.max_cycles==0
        { return Err(WindowError::Invalid("finite IMEX clock, positive step and bounded solver policy required")); }
        let mut total=0usize;
        for &count in counts {
            poll(cancelled)?;
            if count==0 { return Err(WindowError::Invalid("each IMEX interval needs at least one step")); }
            total=total.checked_add(count).ok_or(WindowError::Invalid("IMEX clock-work overflow"))?;
            if total>config.max_clock_steps { return Err(WindowError::Invalid("IMEX clock-work limit")); }
        }
        let mut times=zeros(counts.len()+1)?;times[0]=start;
        let mut time=start;
        for (i,&count) in counts.iter().enumerate() {
            for j in 0..count {
                if j%256==0 {poll(cancelled)?;}
                let next=time+config.step;
                if !next.is_finite() || next<=time {
                    return Err(WindowError::Invalid("IMEX step cannot advance the declared clock"));
                }
                time=next;
            }
            if !(time-times[i]).is_finite() {return Err(WindowError::Invalid("IMEX interval duration overflow"));}
            times[i+1]=time;
        }
        let mut owned=Vec::new();owned.try_reserve_exact(counts.len()).map_err(|_|WindowError::Allocation)?;
        owned.extend_from_slice(counts);poll(cancelled)?;
        Ok(Self{times,counts:owned,config,primal,adjoint})
    }
    pub fn times(&self)->&[f64] {&self.times}
    pub fn step_counts(&self)->&[usize] {&self.counts}
}

fn poll(cancelled:&mut dyn FnMut()->bool)->Result<(),WindowError> {
    if cancelled() {Err(WindowError::Cancelled)} else {Ok(())}
}
fn zeros(n:usize)->Result<Vec<f64>,WindowError> {
    let mut values=Vec::new();values.try_reserve_exact(n).map_err(|_|WindowError::Allocation)?;
    values.resize(n,0.0);Ok(values)
}

fn refusal(interval:usize,phase:&'static str,error:ImexTrajectoryError)->WindowError {
    match error {
        ImexTrajectoryError::Step(ImexAdjointError::Step(ImexSolveError::Cancelled))=>WindowError::Cancelled,
        other=>WindowError::Integrator{interval,phase,diagnostic:other.to_string()},
    }
}

/// Complete native recording; retained gradients refer to its exact endpoint.
/// Original IMEX numerical refusals retain their full diagnostic text and phase
/// in WindowError::Integrator. Cancellation remains separately attributable.
pub struct ImexWindowTape<'a,M,P,Q> {
    tape:RecordedImex2<'a,M,P>,
    adjoint:&'a Q,
    budget:ImexReplayBudget,
    interval:usize,
}
impl<M:ImexVjp,P:FlexiblePreconditioner,Q:FlexiblePreconditioner> IntervalScheme<M> for ImexWindowPolicy<'_,P,Q> {
    type Tape<'a> = ImexWindowTape<'a,M,P,Q> where Self:'a,M:'a;
    fn dimension(&self,model:&M)->usize {model.n()}
    fn parameter_count(&self,model:&M)->usize {model.parameter_count()}
    fn validate(&self,times:&[f64])->Result<(),WindowError> {
        if times!=self.times.as_slice() {return Err(WindowError::Invalid("window times must equal the declared IMEX clock"));}
        Ok(())
    }
    fn record<'a>(&'a self,model:&'a M,interval:usize,start:f64,end:f64,initial:&[f64],
        cancelled:&mut dyn FnMut()->bool)->Result<Self::Tape<'a>,WindowError>
    {
        poll(cancelled)?;
        if interval>=self.counts.len() || start!=self.times[interval] || end!=self.times[interval+1]
            || model.n()==0 || initial.len()!=model.n() || initial.iter().any(|v|!v.is_finite())
        {return Err(WindowError::Invalid("IMEX interval clock or model-state shape mismatch"));}
        let method=OperatorImex2::new(model.n(),self.config.step,self.config.solve);
        let mut tape=RecordedImex2::new(method,model,self.primal,start,initial,ImexRecordingConfig{
            steps:self.counts[interval],max_workspace_components:self.config.max_workspace_components,
        }).map_err(|e|refusal(interval,"IMEX recording",e))?;
        let report=tape.advance(self.config.max_forward_steps,self.config.max_records,&mut||cancelled())
            .map_err(|e|refusal(interval,"IMEX forward",e))?;
        if report.status==ImexRecordingStatus::Cancelled {return Err(WindowError::Cancelled);}
        if report.status!=ImexRecordingStatus::ReachedEnd {
            return Err(WindowError::Integrator{interval,phase:"IMEX forward",diagnostic:format!("incomplete recording: {:?}",report.status)});
        }
        if tape.time()!=end {return Err(WindowError::IntervalOutput{interval,what:"IMEX endpoint clock"});}
        poll(cancelled)?;
        Ok(ImexWindowTape{tape,adjoint:self.adjoint,budget:self.config.replay,interval})
    }
}
impl<M:ImexVjp,P:FlexiblePreconditioner,Q:FlexiblePreconditioner> IntervalTape for ImexWindowTape<'_,M,P,Q> {
    fn endpoint(&self)->&[f64] {self.tape.state()}
    fn end_time(&self)->f64 {self.tape.time()}
    fn accepted_steps(&self)->usize {self.tape.accepted_steps()}
    fn pullback(&self,seed:&[f64],direct_parameters:&[f64],cancelled:&mut dyn FnMut()->bool)
        ->Result<TrajectoryGradient,WindowError>
    {
        let gradient=self.tape.pullback(seed,direct_parameters,self.adjoint,self.budget,&mut||cancelled())
            .map_err(|e|refusal(self.interval,"IMEX adjoint",e))?;
        Ok(TrajectoryGradient{initial:gradient.initial,parameters:gradient.parameters,
            replayed_steps:gradient.replayed_steps,peak_checkpoints:gradient.peak_checkpoints})
    }
}

#[cfg(test)]
mod tests;
