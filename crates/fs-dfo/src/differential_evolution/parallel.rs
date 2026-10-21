//! Concurrent physical evaluations, with original DE sampling and selection.
//!
//! Each wave stops at a generation boundary. Coordinates use the same Philox
//! member/generation keys as the serial path. Scores are consumed in member
//! order, not completion order. With no target stop, equal callback allowances
//! give the same populations for deterministic objectives, regardless of wave
//! width or worker count. Target stopping occurs BETWEEN waves: already launched
//! work is charged and retained, not retroactively cancelled or refunded.
//!
//! The existing TilePool bridge propagates the caller's cancellation authority
//! and drains its workers before returning. Successful callback values are cached
//! even when a sibling fails or cancellation occurs. The serial API can resume
//! that same state, retrying only missing/failed values. No worker mutates parents.
//! Heap owned by the physics callback is not bounded by the executor lease; the
//! callback must retain its own solver/memory admission. No speedup is promised.

use super::{DifferentialEvolutionError, DifferentialEvolutionReport,
    DifferentialEvolutionState, DifferentialEvolutionStop, Pending, checkpoint,
    finite, invalid, reserved};
use core::ops::ControlFlow;
use fs_alloc::OperationMemoryLease;
use fs_exec::{Cancelled, Cx, RunId, TileKernel, TilePlan, TilePool};
use std::collections::VecDeque;
use std::fmt;
use std::sync::{Mutex, atomic::{AtomicUsize, Ordering}};

/// Parallel execution never replaces an original simulation failure with fitness.
#[derive(Debug)]
pub enum DifferentialEvolutionParallelError<E> {
    Solver(DifferentialEvolutionError),
    /// Lowest failed member in this wave, independent of worker completion order.
    Evaluation { member: usize, source: E },
    /// Drained executor refusal or worker panic; completed scores remain cached.
    Executor(String),
}
impl<E> From<DifferentialEvolutionError> for DifferentialEvolutionParallelError<E> {
    fn from(error: DifferentialEvolutionError) -> Self { Self::Solver(error) }
}
impl<E: fmt::Display> fmt::Display for DifferentialEvolutionParallelError<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Solver(error) => write!(f, "{error}"),
            Self::Evaluation { member, source } => write!(f, "DE member {member}: {source}"),
            Self::Executor(detail) => write!(f, "DE population executor: {detail}"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for DifferentialEvolutionParallelError<E> {}

impl DifferentialEvolutionState {
    /// Evaluate bounded waves concurrently on a caller-owned production pool.
    ///
    /// `max_wave` limits simultaneous candidate slots, not the number of worker
    /// threads. Both the per-call and lifetime limits count ACTUAL callback
    /// starts, including failed/panicking attempts. Pre-cancelled/unlaunched tiles
    /// spend no callback slot. Complete parents are still published atomically.
    ///
    /// The callback receives its worker's Cx, never another thread's arena. It
    /// must be deterministic from the decision and immutable physical inputs;
    /// do not use worker/wave identity as an additional random model parameter.
    /// The same immutable objective must be used after resuming or switching to
    /// serial execution. Calls on one state are exclusive via `&mut self`.
    ///
    /// # Errors
    /// Invalid wave width, cancellation, original callback error, nonfinite
    /// score, allocation or executor failure. All workers drain before return.
    /// Successful sibling scores survive every error; failed members remain for
    /// an explicitly charged retry. If several callbacks fail, the first member's
    /// error is returned; later errors are dropped but their attempts stay spent.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_parallel<E: Send, F>(
        &mut self, cx: &Cx<'_>, pool: &TilePool, lease: &OperationMemoryLease,
        additional_evaluations: usize, max_wave: usize, evaluate: &F,
    ) -> Result<DifferentialEvolutionReport, DifferentialEvolutionParallelError<E>>
    where F: Fn(&Cx<'_>, &[f64]) -> Result<f64, E> + Sync,
    {
        let mut poll = || cx.checkpoint().is_ok();
        checkpoint(&mut poll)?;
        if max_wave == 0 { return Err(invalid("DE parallel wave width must be positive").into()); }
        let mut spent = 0;
        loop {
            // This is the same ordered publication path used by advance().
            // No new callback can run while an earlier cached result is pending.
            loop {
                checkpoint(&mut poll)?;
                if self.pending.is_none() { self.pending = self.queued.pop_front(); }
                if !self.pending.as_ref().is_some_and(|p| p.value.is_some()) { break; }
                self.consume(&mut poll)?;
            }
            let stop = if self.best.as_ref().is_some_and(|p|
                self.config.target.is_some_and(|target| p.value <= target)) {
                Some(DifferentialEvolutionStop::TargetReached)
            } else if self.evaluations >= self.config.max_evaluations {
                Some(DifferentialEvolutionStop::EvaluationLimit)
            } else if spent >= additional_evaluations {
                Some(DifferentialEvolutionStop::Paused)
            } else { None };
            if let Some(stop) = stop { return Ok(self.report(stop, &mut poll)?); }
            let allowance = max_wave.min(additional_evaluations-spent)
                .min(self.config.max_evaluations-self.evaluations);
            self.prepare_wave(allowance, &mut poll)?;
            let member = if self.population.len() < self.config.population_size {
                self.population.len()
            } else { self.next.len() };
            // References are scoped to the drained run. Only scalar outcomes and
            // local queue indices escape, so no parent/candidate borrow is mutated.
            let (execution, started, slots) = {
                let mut jobs = reserved(allowance.min(1+self.queued.len()))?;
                for (offset, pending) in self.pending.iter().chain(&self.queued).enumerate() {
                    if offset % 256 == 0 { checkpoint(&mut poll)?; }
                    if pending.value.is_none() {
                        jobs.push((offset, pending.coordinates.as_slice()));
                        if jobs.len() == allowance { break; }
                    }
                }
                if jobs.is_empty() { return Err(invalid("DE parallel wave has no unresolved member").into()); }
                let mut slots = reserved(jobs.len())?;
                for &(offset, _) in &jobs { slots.push((offset, Mutex::new(None))); }
                let kernel = Wave { jobs, slots, evaluate, started: AtomicUsize::new(0) };
                let (execution, _) = pool.run_declared_leased_with_cx(
                    cx, &kernel, RunId(cx.stream_key().kernel_id), lease);
                (execution, kernel.started.load(Ordering::Relaxed), kernel.slots)
            };
            // starts <= jobs <= allowance, so both sums fit their admitted caps.
            self.evaluations += started;
            spent += started;
            let mut failure = None;
            for (offset, slot) in slots {
                let outcome = match slot.into_inner() {
                    Ok(outcome) => outcome,
                    Err(poisoned) => {
                        if failure.is_none() {
                            failure = Some(DifferentialEvolutionParallelError::Executor(
                                "poisoned callback outcome slot".into()));
                        }
                        poisoned.into_inner()
                    }
                };
                match outcome {
                    Some(Ok(value)) => match finite(value) {
                        Ok(value) => {
                            let pending = if offset == 0 {
                                self.pending.as_mut().expect("retained first DE member")
                            } else { &mut self.queued[offset-1] };
                            pending.value = Some(value);
                        }
                        Err(error) if failure.is_none() => failure = Some(error.into()),
                        Err(_) => {},
                    },
                    Some(Err(source)) if failure.is_none() => failure = Some(
                        DifferentialEvolutionParallelError::Evaluation { member: member+offset, source }),
                    _ => {},
                }
            }
            // Cache all completed finite values BEFORE cancellation/error escape.
            checkpoint(&mut poll)?;
            execution.map_err(|error| DifferentialEvolutionParallelError::Executor(format!("{error:?}")))?;
            if let Some(failure) = failure { return Err(failure); }
            if started == 0 { return Err(invalid("DE executor completed without evaluating its unresolved wave").into()); }
        }
    }

    fn prepare_wave(&mut self, count: usize, poll: &mut impl FnMut() -> bool)
        -> Result<(), DifferentialEvolutionError>
    {
        // A retained wave may contain gaps from failed or unlaunched callbacks.
        // Finish it before creating any later members or changing its width.
        if !self.queued.is_empty() { return Ok(()); }
        let member = if self.population.len() < self.config.population_size {
            self.population.len()
        } else { self.next.len() };
        let count = count.min(self.config.population_size-member);
        let mut queued = VecDeque::new();
        queued.try_reserve_exact(count).map_err(|_| DifferentialEvolutionError::Allocation)?;
        for offset in usize::from(self.pending.is_some())..count {
            checkpoint(poll)?;
            queued.push_back(Pending { coordinates: self.candidate_at(member+offset, poll)?, value: None });
        }
        checkpoint(poll)?;
        if self.pending.is_none() { self.pending = queued.pop_front(); }
        self.queued = queued;
        Ok(())
    }
}

struct Wave<'a, F, E> {
    jobs: Vec<(usize, &'a [f64])>,
    slots: Vec<(usize, Mutex<Option<Result<f64, E>>>)>,
    evaluate: &'a F,
    started: AtomicUsize,
}
impl<E: Send, F> TileKernel for Wave<'_, F, E>
where F: Fn(&Cx<'_>, &[f64]) -> Result<f64, E> + Sync,
{
    type Out = u64;
    fn tiles(&self) -> TilePlan {
        TilePlan::new("fs-dfo/differential-evolution-wave", self.jobs.len() as u64)
    }
    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, u64> {
        if cx.checkpoint().is_err() { return ControlFlow::Break(Cancelled); }
        let i = tile as usize;
        self.started.fetch_add(1, Ordering::Relaxed);
        let result = (self.evaluate)(cx, self.jobs[i].1);
        // No user code runs under the lock. Retain completed work even if the
        // callback requested cancellation immediately before returning.
        match self.slots[i].1.lock() {
            Ok(mut slot) => *slot = Some(result),
            Err(_) => return ControlFlow::Break(Cancelled),
        }
        ControlFlow::Continue(1)
    }
}

#[cfg(test)]
mod tests;
