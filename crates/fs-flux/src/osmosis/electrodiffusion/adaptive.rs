//! Conservative step doubling with whole-window rollback and physical budgets.
//!
//! Snapshots contain solution inventories and diagnostics only. Original
//! mechanics, membrane, fixed charge and trapped particles always remain with
//! the pair. Neither a coarse attempt nor a failed window changes that physics.

use super::*;

#[derive(Clone, Copy, Debug)]
pub struct IonicAdvanceControls {
    pub maximum_step_s: f64,
    pub minimum_step_s: f64,
    pub relative_tolerance: f64,
    pub absolute_amount_tolerance_mol: f64,
    pub absolute_water_tolerance_m3: f64,
    pub maximum_steps: usize,
    pub maximum_rejections: usize,
    /// Physical whole-window budgets and algebraic limits. The potential-work
    /// counter is shared by all coarse, fine and refused implicit calls.
    pub step: IonicStepControls,
}
impl Default for IonicAdvanceControls {
    fn default() -> Self {
        Self {
            maximum_step_s: 1.0,
            minimum_step_s: 1e-12,
            relative_tolerance: 1e-4,
            absolute_amount_tolerance_mol: 1e-14,
            absolute_water_tolerance_m3: 1e-15,
            maximum_steps: 10_000,
            maximum_rejections: 128,
            step: IonicStepControls::default(),
        }
    }
}
impl IonicAdvanceControls {
    fn validate(self) -> Result<()> {
        self.step.validate()?;
        if !positive(self.maximum_step_s)
            || !positive(self.minimum_step_s)
            || self.minimum_step_s > self.maximum_step_s
            || !positive(self.relative_tolerance)
            || self.relative_tolerance >= 1.0
            || !positive(self.absolute_amount_tolerance_mol)
            || !positive(self.absolute_water_tolerance_m3)
            || self.maximum_steps == 0
        {
            return Err(invalid(
                "finite ionic advancement needs positive temporal accuracy and bounded steps",
            ));
        }
        Ok(())
    }
    fn inner(self, fraction: f64) -> Result<IonicStepControls> {
        let mut step = self.step;
        step.absolute_water_tolerance_m3 = step
            .absolute_water_tolerance_m3
            .min(0.01 * self.absolute_water_tolerance_m3 * fraction);
        step.absolute_amount_tolerance_mol = step
            .absolute_amount_tolerance_mol
            .min(0.01 * self.absolute_amount_tolerance_mol * fraction);
        step.relative_transfer_tolerance = step
            .relative_transfer_tolerance
            .min(0.01 * self.relative_tolerance);
        step.validate()?;
        Ok(step)
    }
}

#[derive(Clone, Debug)]
pub struct IonicAdvanceReport {
    pub interval_s: [f64; 2],
    pub accepted_steps: usize,
    pub rejected_steps: usize,
    pub step_evaluations: usize,
    pub water_evaluations: usize,
    pub salt_evaluations: usize,
    pub potential_evaluations: usize,
    pub maximum_scaled_error: f64,
    pub minimum_accepted_step_s: f64,
    pub maximum_accepted_step_s: f64,
    pub solvent_transfer_m3: f64,
    pub solute_transfer_mol: IonicSpeciesValues,
    /// Exact original-to-final mechanical primitive, not a sum of trial work.
    pub pressure_work_into_solution_j: f64,
    pub total_potential_change_j: f64,
    pub water_balance_residual_m3: f64,
    pub solute_balance_residual_mol: IonicSpeciesValues,
    pub endpoint: IonicStepReport,
}

struct State {
    solutions: [AqueousSolution; 2],
    transport: ElectrodiffusionReport,
    time_s: f64,
}
impl State {
    fn new(pair: &IonicOsmoticPair) -> Self {
        Self {
            solutions: pair.compartments.each_ref().map(|c| c.solution().clone()),
            transport: pair.endpoint,
            time_s: pair.time_s,
        }
    }
    fn capture(&mut self, pair: &IonicOsmoticPair) -> Result<()> {
        for side in 0..2 {
            let from = pair.compartments[side].solution();
            self.solutions[side]
                .replace_inventories(from.solvent_volume_m3(), from.solute_moles())?;
        }
        self.transport = pair.endpoint;
        self.time_s = pair.time_s;
        Ok(())
    }
    fn restore(&self, pair: &mut IonicOsmoticPair) -> Result<()> {
        for side in 0..2 {
            let from = &self.solutions[side];
            pair.compartments[side]
                .solution
                .replace_inventories(from.solvent_volume_m3(), from.solute_moles())?;
        }
        pair.endpoint = self.transport;
        pair.time_s = self.time_s;
        Ok(())
    }
    fn swap_with(&mut self, pair: &mut IonicOsmoticPair) {
        for side in 0..2 {
            core::mem::swap(
                &mut self.solutions[side],
                &mut pair.compartments[side].solution,
            );
        }
        core::mem::swap(&mut self.transport, &mut pair.endpoint);
        core::mem::swap(&mut self.time_s, &mut pair.time_s);
    }
}

/// Prepared snapshots. Reuse requires the same original solution model; other
/// pair physics is never stored here and cannot be substituted by restoration.
pub struct IonicAdvanceWorkspace {
    original: State,
    start: State,
    coarse: State,
    bytes: usize,
}
impl IonicAdvanceWorkspace {
    pub const fn workspace_bytes(&self) -> usize {
        self.bytes
    }
    fn compatible(&self, pair: &IonicOsmoticPair) -> bool {
        self.original
            .solutions
            .iter()
            .zip(&pair.compartments)
            .all(|(a, cell)| {
                let b = cell.solution();
                a.solute_moles().len() == b.solute_moles().len()
                    && a.temperature_k().to_bits() == b.temperature_k().to_bits()
                    && a.solvent_molar_volume_m3_mol().to_bits()
                        == b.solvent_molar_volume_m3_mol().to_bits()
                    && a.nonideal_model() == b.nonideal_model()
            })
    }
}

impl IonicOsmoticPair {
    pub fn adaptive_workspace(
        &self,
        maximum_workspace_bytes: usize,
    ) -> Result<IonicAdvanceWorkspace> {
        let bytes = core::mem::size_of::<IonicAdvanceWorkspace>()
            + 3 * self
                .compartments
                .iter()
                .map(|c| c.solution().heap_bytes())
                .sum::<usize>();
        if bytes > maximum_workspace_bytes {
            return Err(OsmoticNetworkError::Workspace {
                required_bytes: bytes,
                budget_bytes: maximum_workspace_bytes,
            }
            .into());
        }
        Ok(IonicAdvanceWorkspace {
            original: State::new(self),
            start: State::new(self),
            coarse: State::new(self),
            bytes,
        })
    }

    /// Advance to the representable target with local error control on both
    /// water volumes and every ion amount. Only two actual half steps are
    /// accepted. Any refusal restores all original inventories, diagnostics and
    /// time, including failures of the final complete-window physical checks.
    /// Successful prepared advancement allocates no heap storage.
    pub fn advance(
        &mut self,
        duration_s: f64,
        controls: IonicAdvanceControls,
        workspace: &mut IonicAdvanceWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IonicAdvanceReport> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if !positive(duration_s) || !workspace.compatible(self) {
            return Err(invalid(
                "positive ionic duration and compatible original solution workspace required",
            ));
        }
        let target = finite(self.time_s + duration_s)?;
        if target <= self.time_s {
            return Err(OsmoticNetworkError::StepUnderflow.into());
        }
        workspace.original.capture(self)?;
        let result = self.advance_inner(target, controls, workspace, &mut cancelled);
        if result.is_err() {
            workspace.original.swap_with(self);
        }
        result
    }

    fn advance_inner(
        &mut self,
        target: f64,
        controls: IonicAdvanceControls,
        workspace: &mut IonicAdvanceWorkspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<IonicAdvanceReport> {
        let start_time = self.time_s;
        let duration = target - start_time;
        let mut proposed = controls.maximum_step_s.min(duration);
        let (mut accepted, mut rejected, mut evaluations) = (0, 0, 0);
        let (mut maximum_error, mut minimum_step, mut maximum_step) =
            (0.0_f64, f64::INFINITY, 0.0_f64);
        let mut work = Work::default();
        let mut endpoint = None;
        while self.time_s < target {
            poll(cancelled)?;
            if accepted >= controls.maximum_steps {
                return Err(OsmoticNetworkError::WorkBudget.into());
            }
            workspace.start.capture(self)?;
            let planned = if proposed >= target - self.time_s {
                target
            } else {
                self.time_s + proposed
            };
            let end = if planned < target && planned.next_up() >= target {
                target
            } else {
                planned
            };
            let dt = end - self.time_s;
            let midpoint = self.time_s + 0.5 * dt;
            if midpoint <= self.time_s || midpoint >= end || end > target {
                return Err(OsmoticNetworkError::StepUnderflow.into());
            }
            let (first_dt, second_dt) = (midpoint - self.time_s, end - midpoint);
            let (u, v) = (first_dt / dt, second_dt / dt);
            let factor = (u * u + v * v) / (2.0 * u * v);
            let trial = (|| -> Result<_> {
                bump(&mut evaluations)?;
                self.step_counted(dt, controls.inner(dt / duration)?, &mut work, cancelled)?;
                workspace.coarse.capture(self)?;
                workspace.start.restore(self)?;
                bump(&mut evaluations)?;
                self.step_counted(
                    first_dt,
                    controls.inner(first_dt / duration)?,
                    &mut work,
                    cancelled,
                )?;
                if self.time_s != midpoint {
                    return Err(OsmoticNetworkError::StepUnderflow.into());
                }
                bump(&mut evaluations)?;
                let second = self.step_counted(
                    second_dt,
                    controls.inner(second_dt / duration)?,
                    &mut work,
                    cancelled,
                )?;
                if self.time_s != end {
                    return Err(OsmoticNetworkError::StepUnderflow.into());
                }
                let error =
                    temporal_error(self, &workspace.start, &workspace.coarse, controls)? * factor;
                Ok((error, second))
            })();
            let (error, second) = match trial {
                Ok(value) => value,
                Err(error) if retryable(&error) => {
                    workspace.start.restore(self)?;
                    if rejected >= controls.maximum_rejections || 0.5 * dt < controls.minimum_step_s
                    {
                        return Err(error);
                    }
                    rejected += 1;
                    proposed = 0.5 * dt;
                    continue;
                }
                Err(error) => return Err(error),
            };
            if error > 1.0 {
                workspace.start.restore(self)?;
                if 0.5 * dt < controls.minimum_step_s {
                    return Err(IonicOsmoticError::TemporalAccuracy {
                        step_s: dt,
                        scaled_error: error,
                    });
                }
                if rejected >= controls.maximum_rejections {
                    return Err(OsmoticNetworkError::WorkBudget.into());
                }
                rejected += 1;
                proposed = 0.5 * dt;
                continue;
            }
            accepted += 1;
            maximum_error = maximum_error.max(error);
            minimum_step = minimum_step.min(dt);
            maximum_step = maximum_step.max(dt);
            endpoint = Some(second);
            proposed = if error < 0.125 {
                (2.0 * dt).min(controls.maximum_step_s)
            } else {
                dt
            };
        }
        let before = workspace.original.solutions.each_ref();
        let after = self.compartments.each_ref().map(|c| c.solution());
        let (water, solute) = self.balance(before, after, controls.step)?;
        let (pressure_work, potential_change) = self.potential_change(before, after)?;
        if potential_change > controls.step.free_energy_tolerance_j {
            return Err(OsmoticNetworkError::CoupledFreeEnergy {
                change_j: potential_change,
                maximum_j: controls.step.free_energy_tolerance_j,
            }
            .into());
        }
        let report = IonicAdvanceReport {
            interval_s: [start_time, target],
            accepted_steps: accepted,
            rejected_steps: rejected,
            step_evaluations: evaluations,
            water_evaluations: work.water,
            salt_evaluations: work.salt,
            potential_evaluations: work.potential,
            maximum_scaled_error: maximum_error,
            minimum_accepted_step_s: minimum_step,
            maximum_accepted_step_s: maximum_step,
            solvent_transfer_m3: before[0].solvent_volume_m3() - after[0].solvent_volume_m3(),
            solute_transfer_mol: IonicSpeciesValues::from_fn(self.species_count(), |i| {
                before[0].solute_moles()[i] - after[0].solute_moles()[i]
            }),
            pressure_work_into_solution_j: pressure_work,
            total_potential_change_j: potential_change,
            water_balance_residual_m3: water,
            solute_balance_residual_mol: solute,
            endpoint: endpoint.ok_or(OsmoticNetworkError::StepUnderflow)?,
        };
        poll(cancelled)?;
        Ok(report)
    }
}

fn temporal_error(
    pair: &IonicOsmoticPair,
    start: &State,
    coarse: &State,
    controls: IonicAdvanceControls,
) -> Result<f64> {
    let mut error = 0.0_f64;
    for side in 0..2 {
        let (a, b, c) = (
            &start.solutions[side],
            &coarse.solutions[side],
            pair.compartments[side].solution(),
        );
        let mut sample = |a: f64, b: f64, c: f64, absolute| -> Result<()> {
            let allowance =
                finite(absolute + controls.relative_tolerance * a.abs().max(b.abs()).max(c.abs()))?;
            error = error.max(finite((b - c).abs() / allowance)?);
            Ok(())
        };
        sample(
            a.solvent_volume_m3(),
            b.solvent_volume_m3(),
            c.solvent_volume_m3(),
            controls.absolute_water_tolerance_m3,
        )?;
        for i in 0..pair.species_count() {
            sample(
                a.solute_moles()[i],
                b.solute_moles()[i],
                c.solute_moles()[i],
                controls.absolute_amount_tolerance_mol,
            )?;
        }
    }
    Ok(error)
}
fn retryable(error: &IonicOsmoticError) -> bool {
    matches!(
        error,
        IonicOsmoticError::Network(
            OsmoticNetworkError::Unbracketed
                | OsmoticNetworkError::NonlinearConvergence
                | OsmoticNetworkError::CoupledStepFreeEnergy { .. }
        ) | IonicOsmoticError::Membrane(ElectrodiffusionError::NonConvergence)
    )
}
