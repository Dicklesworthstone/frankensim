//! Deterministic, transactional internal refinement of a physical output tick.
//! Same prepared solver and material acceptance; no reset, clipping or skipped time.
use super::{PreparedImpactSystem, ImpactError, ImpactFrame, ImpactSource, CancelGate};
use super::super::{invalid, ImpactConfig};
use fs_material::fiber::WoolFeltState;
use std::ops::Deref;

/// Explicit work bound, not an accuracy estimator or hard-real-time guarantee.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ImpactSubstepConfig {
    /// Smallest trial duration is the output period divided by 2^max_depth.
    /// At most 1024 accepted leaves per output tick (depth 0..=10).
    pub max_depth: u8,
    /// Total nonlinear solve attempts, including refused parents (1..=2047).
    /// A smaller limit may refuse even before maximum depth is reached.
    pub max_attempts: usize,
}
impl ImpactSubstepConfig {
    fn validate(self, dt: f64) -> Result<(), ImpactError> {
        if self.max_depth > 10 || self.max_attempts == 0 || self.max_attempts > 2047
            || !dt.is_finite() || dt <= 0.0 || dt / f64::from(1u32 << self.max_depth) == 0.0 {
            return Err(invalid("impact substeps require depth 0..10, 1..2047 attempts and a representable positive leaf period"));
        }
        Ok(())
    }
}

/// Work used by the last completely accepted output tick; refusals leave it alone.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ImpactSubstepReport {
    pub attempted_solves: usize,
    pub accepted_substeps: usize,
    pub deepest_level: u8,
}

/// Nonlinear mechanics with internal dyadic refinement and an unchanged output clock.
///
/// The held force is applied to every leaf for its actual duration. Felt history
/// advances at each accepted leaf; a failed *output tick* restores all such
/// changes. No observer, force program or audio filter sees a partial tick.
/// Storage and solver callbacks must avoid allocation for allocation-free steps.
/// This controls solve failure, NOT temporal discretization or aliasing error.
pub struct SubsteppedImpactSystem {
    prepared: PreparedImpactSystem,
    config: ImpactSubstepConfig,
    rollback_state: Vec<f64>,
    rollback_history: Vec<WoolFeltState>,
    report: ImpactSubstepReport,
}
impl PreparedImpactSystem {
    /// Enable bounded internal refinement without changing motion, history or time.
    pub fn with_substeps(self, config: ImpactSubstepConfig) -> Result<SubsteppedImpactSystem, ImpactError> {
        config.validate(self.inner.config.dt_s)?;
        let divisor = f64::from(1u32 << config.max_depth);
        if self.inner.config.energy_absolute_tolerance_j / divisor == 0.0
            || self.inner.config.energy_relative_tolerance / divisor == 0.0 {
            return Err(invalid("substep energy allowances are not representable"));
        }
        let rollback_state = self.inner.x.clone();
        let rollback_history = self.inner.histories.borrow().clone();
        Ok(SubsteppedImpactSystem { prepared: self, config, rollback_state, rollback_history,
            report: ImpactSubstepReport::default() })
    }
}
impl Deref for SubsteppedImpactSystem {
    type Target = PreparedImpactSystem;
    fn deref(&self) -> &Self::Target { &self.prepared }
}
impl SubsteppedImpactSystem {
    pub fn last_substeps(&self) -> ImpactSubstepReport { self.report }
    /// Disable internal refinement without resetting the accepted instrument.
    pub fn into_prepared(self) -> PreparedImpactSystem { self.prepared }
    pub fn set_analytic_newton(&mut self, enabled: bool) { self.prepared.set_analytic_newton(enabled); }
    pub fn set_iteration_limit(&mut self, limit: usize) -> Result<(), ImpactError> {
        self.prepared.set_iteration_limit(limit)
    }

    /// Advance exactly one nominal sample, or publish nothing. Only numerical
    /// Newton/energy refusals cause bisection; invalid physics is never repaired.
    /// All attempts share the original held input. Dyadic integer addresses
    /// prevent gaps, overlaps and accumulation of floating-point clock drift.
    pub fn step(&mut self, external: &[f64], gate: &CancelGate) -> Result<ImpactFrame, ImpactError> {
        if gate.is_requested() { return Err(ImpactError::Cancelled); }
        let inner = &self.prepared.inner;
        if inner.sample >= inner.config.max_steps { return Err(ImpactError::Budget); }
        if external.len() != inner.modes || external.iter().any(|v|
            !v.is_finite() || v.abs() > inner.config.maximum_generalized_force) {
            return Err(invalid("external generalized force shape or ceiling failed"));
        }
        let sample = inner.sample;
        let nominal = inner.config;
        let before = inner.stored_energy_j();
        self.rollback_state.copy_from_slice(&inner.x);
        self.rollback_history.clone_from_slice(&inner.histories.borrow());
        let result = self.advance_tick(external, gate, sample, nominal, before);
        // The inner sample counter counts OUTPUT ticks, not solver trials.
        // Its temporary leaf frame times are never returned to a consumer.
        self.prepared.inner.config = nominal;
        match result {
            Ok((frame, report)) => {
                self.prepared.inner.sample = sample + 1;
                self.report = report;
                Ok(frame)
            }
            Err(error) => {
                self.prepared.inner.x.copy_from_slice(&self.rollback_state);
                self.prepared.inner.histories.borrow_mut().clone_from_slice(&self.rollback_history);
                self.prepared.inner.sample = sample;
                Err(error)
            }
        }
    }

    fn advance_tick(&mut self, external: &[f64], gate: &CancelGate, sample: u64,
        nominal: ImpactConfig, before: f64,
    ) -> Result<(ImpactFrame, ImpactSubstepReport), ImpactError> {
        let total = 1u32 << self.config.max_depth;
        let (mut at, mut width) = (0u32, total);
        let mut report = ImpactSubstepReport::default();
        let mut sums = [0.0;3];
        let mut corrections = [0.0;3];
        let mut residual = 0.0_f64;
        let mut after = before;
        while at < total {
            if gate.is_requested() { return Err(ImpactError::Cancelled); }
            if report.attempted_solves >= self.config.max_attempts {
                return Err(ImpactError::SubstepBudget { attempted_solves: report.attempted_solves,
                    accepted_substeps: report.accepted_substeps, deepest_level: report.deepest_level });
            }
            let level = (total / width).trailing_zeros() as u8;
            report.deepest_level = report.deepest_level.max(level);
            report.attempted_solves += 1;
            let fraction = f64::from(width) / f64::from(total);
            self.prepared.inner.config = nominal;
            self.prepared.inner.config.dt_s = nominal.dt_s * fraction;
            // Subdivision must not multiply the parent tick's energy allowance.
            self.prepared.inner.config.energy_absolute_tolerance_j *= fraction;
            self.prepared.inner.config.energy_relative_tolerance *= fraction;
            self.prepared.inner.sample = sample;
            match self.prepared.step(external, gate) {
                Ok(frame) => {
                    report.accepted_substeps += 1;
                    after = frame.stored_energy_j;
                    residual = residual.max(frame.solver_residual);
                    // Compensated sums keep long contact/felt subdivisions from
                    // spending the output tick's energy tolerance on bookkeeping.
                    for (i, value) in [frame.dissipated_energy_j, frame.felt_crush_loss_j,
                        frame.supplied_work_j].into_iter().enumerate() {
                        let y = value - corrections[i];
                        let t = sums[i] + y;
                        corrections[i] = (t - sums[i]) - y;
                        sums[i] = t;
                    }
                    at += width;
                    // Return to a coarser sibling only at its exact boundary.
                    while width < total && at % (2 * width) == 0 && 2 * width <= total - at {
                        width *= 2;
                    }
                }
                Err(error) => {
                    if !matches!(&error, ImpactError::PreparedSolve(fs_phs::PhsError::NewtonStalled { .. })
                        | ImpactError::Energy { .. }) {
                        return Err(error);
                    }
                    if width == 1 {
                        return Err(error); // Preserve the actual terminal numerical refusal.
                    }
                    width /= 2;
                }
            }
        }
        let balance = after - before + sums[0] - sums[2];
        let config = nominal;
        let tolerance = config.energy_absolute_tolerance_j + config.energy_relative_tolerance
            * (before.abs() + after.abs() + sums[0].abs() + sums[2].abs());
        if !balance.is_finite() || balance.abs() > tolerance {
            return Err(ImpactError::Energy { residual_j: balance, tolerance_j: tolerance });
        }
        if gate.is_requested() { return Err(ImpactError::Cancelled); }
        Ok((ImpactFrame { sample: sample + 1, time_s: (sample + 1) as f64 * nominal.dt_s,
            stored_energy_j: after, dissipated_energy_j: sums[0], felt_crush_loss_j: sums[1],
            supplied_work_j: sums[2], balance_residual_j: balance, solver_residual: residual }, report))
    }
}
impl ImpactSource for SubsteppedImpactSystem {
    type Frame = ImpactFrame;
    fn state(&self) -> &[f64] { self.prepared.state() }
    fn mode_count(&self) -> usize { self.prepared.mode_count() }
    fn samples_rendered(&self) -> u64 { self.prepared.samples_rendered() }
    fn sample_period_s(&self) -> f64 { self.prepared.sample_period_s() }
    fn remaining_steps(&self) -> u64 { self.prepared.remaining_steps() }
    fn maximum_generalized_force(&self) -> f64 { self.prepared.maximum_generalized_force() }
    fn advance(&mut self, forces: &[f64], gate: &CancelGate) -> Result<Self::Frame, ImpactError> {
        self.step(forces, gate)
    }
}

#[cfg(test)]
mod tests;
