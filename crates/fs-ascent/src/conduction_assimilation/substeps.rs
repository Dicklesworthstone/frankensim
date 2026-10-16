//! Conduction-specific schedule binding over the shared substep replay engine.
//!
//! Refinement changes numerical time resolution, NOT the observation timetable,
//! knot-state controls, or endpoint model-error penalties. The existing
//! CheckpointedIntervals owner records, hashes and reverses composed maps.
//! This adapter only maps its global fine-step addresses back to the original
//! observation interval and local substep for source/parameter callbacks.

use super::{ConductionInterval, ConductionSubstep, ConductionWindowModel,
    ConductionWindowPolicy, finite, poll};
use crate::transient::variational::{WindowError, intervals::{IntervalScheme, IntervalTape}};
use crate::transient::variational::intervals::substeps::{
    CheckpointedIntervals, SubstepBudget, SubstepGrid, SubstepTape,
};
use fs_time::adaptive::adjoint::trajectory::TrajectoryGradient;

#[derive(Debug, Clone, Copy)]
pub struct SubstepLimits {
    /// Sum of subdivisions across the whole observation timetable. Each is
    /// positive; fixed-clock storage is linear in this declared bound.
    pub max_steps: usize,
    /// Conservative retained entries per interval: two free fields plus four
    /// u64 hash words and one step count per substep. Policy clock storage, Vec
    /// descriptors and the one live sparse PDE linearization are separate.
    pub max_record_components: usize,
    /// Parked replay states including the borrowed initial state. A multistep
    /// interval needs floor(log2(steps))+1; a single step retains its original
    /// linearization and does not replay. Live cotangents are separate.
    pub checkpoints: usize,
    /// Maximum primal step reconstructions in EACH interval reverse sweep.
    /// Every reconstruction has the physical policy's own solve/work limits.
    pub replayed_steps: usize,
}

/// Refines a base conduction policy without reinterpreting its source indices.
/// All mesh, boundary, physical Cx and covariance semantics remain unchanged.
/// Clone copies numerical policy, not the variational work allowance.
#[derive(Clone)]
pub struct ConductionSubsteps<'a, 'cx> {
    composed: CheckpointedIntervals<FineConduction<'a, 'cx>>,
}
#[derive(Clone)]
struct FineConduction<'a, 'cx> {
    base: ConductionWindowPolicy<'a, 'cx>,
    grid: SubstepGrid,
}
impl<'a, 'cx> ConductionSubsteps<'a, 'cx> {
    pub fn new(base: ConductionWindowPolicy<'a, 'cx>, subdivisions: &[usize],
        limits: SubstepLimits, cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, WindowError> {
        poll(cancelled)?;
        if subdivisions.len() != base.times().len()-1 {
            return Err(WindowError::Invalid("one subdivision count per observation interval required"));
        }
        let mut total = 0usize;
        for &count in subdivisions {
            poll(cancelled)?;
            if count == 0 { return Err(WindowError::Invalid("subdivision counts must be positive")); }
            total = total.checked_add(count).ok_or(WindowError::Invalid("substep count overflow"))?;
            if total > limits.max_steps { return Err(WindowError::Invalid("conduction substep limit")); }
            let required = base.dimension().checked_mul(2)
                .and_then(|n| count.checked_mul(5).and_then(|r| n.checked_add(r)))
                .ok_or(WindowError::Invalid("substep record extent overflow"))?;
            if required > limits.max_record_components {
                return Err(WindowError::WorkspaceLimit { required, limit: limits.max_record_components });
            }
            let slots = if count == 1 { 0 } else { (usize::BITS-count.leading_zeros()) as usize };
            if limits.checkpoints < slots {
                return Err(WindowError::Invalid("insufficient conduction replay checkpoints"));
            }
        }
        let mut fine = Vec::new(); let mut indices = Vec::new();
        fine.try_reserve_exact(total.checked_add(1).ok_or(WindowError::Invalid("substep grid overflow"))?)
            .map_err(|_| WindowError::Allocation)?;
        indices.try_reserve_exact(base.times().len()).map_err(|_| WindowError::Allocation)?;
        fine.push(base.times()[0]); indices.push(0);
        for (k, &count) in subdivisions.iter().enumerate() {
            let (start, end) = (base.times()[k], base.times()[k+1]);
            for j in 1..=count {
                if j % 256 == 0 { poll(cancelled)?; }
                fine.push(if j == count { end } else { (end-start).mul_add(j as f64/count as f64, start) });
            }
            indices.push(fine.len()-1);
        }
        // The existing grid owner admits ordered representable endpoints. Use
        // its fallible constructor for both copies instead of infallible clone.
        let grid = SubstepGrid::new(&fine, &indices, limits.max_steps)?;
        let schedule_grid = SubstepGrid::new(&fine, &indices, limits.max_steps)?;
        let budget = SubstepBudget { max_state_components: base.dimension(),
            max_parameters: base.config.max_parameters, checkpoints: limits.checkpoints,
            replayed_substeps: limits.replayed_steps };
        let inner = FineConduction { base, grid: schedule_grid };
        poll(cancelled)?;
        Ok(Self { composed: CheckpointedIntervals::new(inner, grid, budget) })
    }
    pub fn base(&self) -> &ConductionWindowPolicy<'a, 'cx> { &self.composed.inner().base }
    pub fn times(&self) -> &[f64] { self.composed.grid().knot_times() }
    /// Exact internal clock, including the original observation endpoints.
    pub fn substep_times(&self, interval: usize) -> Option<&[f64]> {
        let grid = self.composed.grid(); let indices = grid.knot_indices();
        let end = interval.checked_add(1)?;
        if end >= indices.len() { return None; }
        Some(&grid.fine_times()[indices[interval]..=indices[end]])
    }
}
impl<'cx, M: ConductionWindowModel> IntervalScheme<M> for FineConduction<'_, 'cx> {
    type Tape<'a> = ConductionInterval<'a, 'cx, M> where Self: 'a, M: 'a;
    fn dimension(&self, _model: &M) -> usize { self.base.dimension() }
    fn parameter_count(&self, model: &M) -> usize { model.parameter_count() }
    fn validate(&self, fine: &[f64]) -> Result<(), WindowError> {
        if fine != self.grid.fine_times() { return Err(WindowError::Invalid("conduction fine clock changed")); }
        Ok(())
    }
    fn record<'a>(&'a self, model: &'a M, index: usize, start: f64, end: f64,
        initial: &[f64], check: &mut dyn FnMut() -> bool,
    ) -> Result<Self::Tape<'a>, WindowError> {
        poll(check)?;
        let fine = self.grid.fine_times();
        if index >= fine.len()-1 || start != fine[index] || end != fine[index+1] {
            return Err(WindowError::Invalid("conduction substep address mismatch"));
        }
        let knots = self.grid.knot_indices();
        let interval = knots.partition_point(|&i| i <= index)-1;
        let time = ConductionSubstep { interval, index: index-knots[interval], start, end };
        self.base.record_step(model, interval, start, end, initial, Some(time), check)
    }
}

/// Complete inner map; the one-step branch preserves the original tape,
/// arithmetic and zero-replay behavior. Multistep sweeps use the shared owner.
pub struct ConductionSubstepTape<'a, 'cx, M> { inner: TapeKind<'a, 'cx, M> }
enum TapeKind<'a, 'cx, M> {
    Single(ConductionInterval<'a, 'cx, M>),
    Multiple(SubstepTape<'a, M, FineConduction<'a, 'cx>>),
}
impl<'cx, M: ConductionWindowModel> IntervalScheme<M> for ConductionSubsteps<'_, 'cx> {
    type Tape<'a> = ConductionSubstepTape<'a, 'cx, M> where Self: 'a, M: 'a;
    fn dimension(&self, _model: &M) -> usize { self.base().dimension() }
    fn parameter_count(&self, model: &M) -> usize { model.parameter_count() }
    fn validate(&self, times: &[f64]) -> Result<(), WindowError> {
        if times != self.times() { return Err(WindowError::Invalid("window must use the original observation clock")); }
        Ok(())
    }
    fn record<'a>(&'a self, model: &'a M, interval: usize, start: f64, end: f64,
        initial: &[f64], check: &mut dyn FnMut() -> bool,
    ) -> Result<Self::Tape<'a>, WindowError> {
        poll(check)?;
        let grid = self.substep_times(interval).ok_or(WindowError::Invalid("unknown conduction interval"))?;
        if start != grid[0] || end != grid[grid.len()-1] || initial.len() != self.base().dimension()
            || !finite(initial) || model.parameter_count() > self.base().config.max_parameters
        { return Err(WindowError::Invalid("substep interval, state or parameter mismatch")); }
        let inner = if grid.len() == 2 {
            let index = self.composed.grid().knot_indices()[interval];
            TapeKind::Single(self.composed.inner().record(model, index, start, end, initial, check)?)
        } else {
            TapeKind::Multiple(self.composed.record(model, interval, start, end, initial, check)?)
        };
        poll(check)?;
        Ok(ConductionSubstepTape { inner })
    }
}
impl<M: ConductionWindowModel> IntervalTape for ConductionSubstepTape<'_, '_, M> {
    fn endpoint(&self) -> &[f64] {
        match &self.inner { TapeKind::Single(t) => t.endpoint(), TapeKind::Multiple(t) => t.endpoint() }
    }
    fn end_time(&self) -> f64 {
        match &self.inner { TapeKind::Single(t) => t.end_time(), TapeKind::Multiple(t) => t.end_time() }
    }
    fn accepted_steps(&self) -> usize {
        match &self.inner { TapeKind::Single(t) => t.accepted_steps(), TapeKind::Multiple(t) => t.accepted_steps() }
    }
    fn pullback(&self, seed: &[f64], direct: &[f64], check: &mut dyn FnMut() -> bool)
        -> Result<TrajectoryGradient, WindowError> {
        match &self.inner {
            TapeKind::Single(t) => t.pullback(seed, direct, check),
            TapeKind::Multiple(t) => t.pullback(seed, direct, check),
        }
    }
}

#[cfg(test)]
mod tests;
