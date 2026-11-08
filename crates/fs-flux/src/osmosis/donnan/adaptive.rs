//! Prepared temporal error control for finite Donnan water/ion relaxation.
//!
//! One backward-Euler step is compared with two half steps from the same
//! inventories. Only the conservative fine state is retained: extrapolation
//! would discard positivity and the original fast-ion equilibrium closure.
//! Every mobile ion amount and both solvent volumes enter the local error
//! norm. This estimate is neither a global-error nor a model-validity bound.

use super::super::compensated_add;
use super::*;

/// Temporal accuracy and bounded work, separate from each algebraic solve.
#[derive(Clone, Copy, Debug)]
pub struct DonnanAdvanceControls {
    /// Largest proposed interval; one terminal clock ULP may be folded into it.
    pub maximum_step_s: f64,
    /// A final remainder may be smaller, but must admit two distinct halves.
    pub minimum_step_s: f64,
    /// Relative local-error allowance for every volume and mobile ion amount.
    pub relative_tolerance: f64,
    /// Absolute allowance for EACH mobile ion on EACH side, including zeros.
    pub absolute_amount_tolerance_mol: f64,
    /// Absolute local-error allowance for each solvent volume.
    pub absolute_water_tolerance_m3: f64,
    pub maximum_steps: usize,
    pub maximum_rejections: usize,
    /// Algebraic accuracy/work caps. Advancement only tightens these: absolute
    /// water and charge allowances and relative charge accuracy are distributed
    /// over the window at 1% of corresponding temporal tolerances. Charge
    /// accuracy retains the original charge-scale definition; it is not a
    /// separate certificate of each ion's relative partition accuracy.
    pub step: DonnanStepControls,
}
impl Default for DonnanAdvanceControls {
    fn default() -> Self {
        Self {
            maximum_step_s: 1.0,
            minimum_step_s: 1e-12,
            relative_tolerance: 1e-4,
            absolute_amount_tolerance_mol: 1e-14,
            absolute_water_tolerance_m3: 1e-15,
            maximum_steps: 10_000,
            maximum_rejections: 128,
            step: DonnanStepControls::default(),
        }
    }
}

/// Complete accepted window with original-to-final conservation residuals.
#[derive(Clone, Debug)]
pub struct DonnanAdvanceReport {
    pub interval_s: [f64; 2],
    pub accepted_steps: usize,
    pub rejected_steps: usize,
    /// All implicit calls, including coarse and refused attempts.
    pub step_evaluations: usize,
    /// All attempted water-residual evaluations, including refused attempts.
    pub root_evaluations: usize,
    /// Evaluations from completed ion equilibria, as in [`DonnanStepReport`].
    pub ion_root_evaluations: usize,
    pub maximum_scaled_error: f64,
    pub minimum_accepted_step_s: f64,
    pub maximum_accepted_step_s: f64,
    pub solvent_transfer_m3: f64,
    /// Compensated sum over accepted half steps only; not a caloric balance.
    pub pressure_work_into_solution_j: f64,
    pub water_balance_residual_m3: f64,
    pub solute_balance_residual_mol: [f64; MAX_DONNAN_SPECIES],
    pub species_count: usize,
    /// Last accepted HALF step, including its actual equilibrium and residuals.
    pub endpoint: DonnanStepReport,
}

struct State {
    compartments: [OsmoticCompartment; 2],
    equilibrium: DonnanReport,
    time_s: f64,
}
impl State {
    fn new(pair: &DonnanOsmoticPair) -> Self {
        Self {
            compartments: pair.compartments.clone(),
            equilibrium: pair.equilibrium,
            time_s: pair.time_s,
        }
    }
    fn capture(&mut self, pair: &DonnanOsmoticPair) -> Result<()> {
        copy_cells(&pair.compartments, &mut self.compartments)?;
        self.equilibrium = pair.equilibrium;
        self.time_s = pair.time_s;
        Ok(())
    }
    fn restore(&self, pair: &mut DonnanOsmoticPair) -> Result<()> {
        copy_cells(&self.compartments, &mut pair.compartments)?;
        pair.equilibrium = self.equilibrium;
        pair.time_s = self.time_s;
        Ok(())
    }
    // Infallible rollback after even a late cancellation or partially restored
    // trial; private nonlinear scratch is overwritten by the next step.
    fn swap_with(&mut self, pair: &mut DonnanOsmoticPair) {
        core::mem::swap(&mut self.compartments, &mut pair.compartments);
        core::mem::swap(&mut self.equilibrium, &mut pair.equilibrium);
        core::mem::swap(&mut self.time_s, &mut pair.time_s);
    }
}
fn copy_cells(from: &[OsmoticCompartment; 2], to: &mut [OsmoticCompartment; 2]) -> Result<()> {
    for (source, target) in from.iter().zip(to) {
        target.solution.replace_inventories(
            source.solution().solvent_volume_m3(),
            source.solution().solute_moles(),
        )?;
        target.pressure = source.pressure;
        target.minimum_volume_m3 = source.minimum_volume_m3;
    }
    Ok(())
}

/// Three prepared snapshots; single-step callers pay no additional storage.
pub struct DonnanAdvanceWorkspace {
    original: State,
    start: State,
    coarse: State,
    species_count: usize,
    bytes: usize,
}
impl DonnanAdvanceWorkspace {
    /// Additional owned storage, including the workspace itself.
    pub const fn workspace_bytes(&self) -> usize {
        self.bytes
    }
    fn compatible(&self, pair: &DonnanOsmoticPair) -> bool {
        self.species_count == pair.mobile_totals_mol.len()
            && self
                .original
                .compartments
                .iter()
                .zip(&pair.compartments)
                .all(|(a, b)| {
                    let a = a.solution();
                    let b = b.solution();
                    a.temperature_k().to_bits() == b.temperature_k().to_bits()
                        && a.solvent_molar_volume_m3_mol().to_bits()
                            == b.solvent_molar_volume_m3_mol().to_bits()
                })
    }
}

impl DonnanAdvanceControls {
    fn validate(&self) -> Result<()> {
        let ions = self.step.ion_equilibrium;
        if !positive(self.maximum_step_s)
            || !positive(self.minimum_step_s)
            || self.minimum_step_s > self.maximum_step_s
            || !positive(self.relative_tolerance)
            || self.relative_tolerance >= 1.0
            || !positive(self.absolute_amount_tolerance_mol)
            || !positive(self.absolute_water_tolerance_m3)
            || self.maximum_steps == 0
            || !positive(self.step.absolute_water_tolerance_m3)
            || !positive(self.step.relative_water_tolerance)
            || self.step.relative_water_tolerance >= 1.0
            || !(2..=4096).contains(&self.step.maximum_root_evaluations)
            || !ions.absolute_charge_tolerance_mol.is_finite()
            || ions.absolute_charge_tolerance_mol < 0.0
            || !positive(ions.relative_charge_tolerance)
            || ions.relative_charge_tolerance >= 1.0
            || !(1..=4096).contains(&ions.maximum_iterations)
        {
            return Err(invalid("invalid adaptive Donnan controls"));
        }
        Ok(())
    }
    fn inner(&self, fraction: f64) -> Result<DonnanStepControls> {
        let mut inner = self.step;
        inner.absolute_water_tolerance_m3 = inner
            .absolute_water_tolerance_m3
            .min(0.01 * self.absolute_water_tolerance_m3 * fraction);
        inner.relative_water_tolerance = inner
            .relative_water_tolerance
            .min(0.01 * self.relative_tolerance);
        inner.ion_equilibrium.absolute_charge_tolerance_mol = inner
            .ion_equilibrium
            .absolute_charge_tolerance_mol
            .min(0.01 * self.absolute_amount_tolerance_mol * fraction);
        inner.ion_equilibrium.relative_charge_tolerance = inner
            .ion_equilibrium
            .relative_charge_tolerance
            .min(0.01 * self.relative_tolerance * fraction);
        if !positive(inner.absolute_water_tolerance_m3)
            || !positive(inner.relative_water_tolerance)
            || !positive(inner.ion_equilibrium.relative_charge_tolerance)
        {
            return Err(invalid(
                "adaptive Donnan algebraic tolerances are not representable",
            ));
        }
        Ok(inner)
    }
}

fn scaled_error(start: f64, coarse: f64, fine: f64, absolute: f64, relative: f64) -> Result<f64> {
    let allowance = finite(
        absolute + relative * start.abs().max(coarse.abs()).max(fine.abs()),
        "nonfinite Donnan temporal allowance",
    )?;
    finite(
        (coarse - fine).abs() / allowance,
        "nonfinite Donnan temporal error",
    )
}
fn temporal_error(
    pair: &DonnanOsmoticPair,
    start: &State,
    coarse: &State,
    controls: DonnanAdvanceControls,
) -> Result<f64> {
    let mut error: f64 = 0.0;
    for side in 0..2 {
        let (a, b, c) = (
            start.compartments[side].solution(),
            coarse.compartments[side].solution(),
            pair.compartments[side].solution(),
        );
        error = error.max(scaled_error(
            a.solvent_volume_m3(),
            b.solvent_volume_m3(),
            c.solvent_volume_m3(),
            controls.absolute_water_tolerance_m3,
            controls.relative_tolerance,
        )?);
        for i in 0..pair.mobile_totals_mol.len() {
            error = error.max(scaled_error(
                a.solute_moles()[i],
                b.solute_moles()[i],
                c.solute_moles()[i],
                controls.absolute_amount_tolerance_mol,
                controls.relative_tolerance,
            )?);
        }
    }
    Ok(error)
}
fn retryable(error: &DonnanWaterError) -> bool {
    matches!(
        error,
        DonnanWaterError::Network(
            OsmoticNetworkError::Unbracketed | OsmoticNetworkError::NonlinearConvergence
        ) | DonnanWaterError::Ions(DonnanError::NonConvergence)
            | DonnanWaterError::InventoryBalance { .. }
    )
}

impl DonnanOsmoticPair {
    /// Prepare optional snapshots under a separate byte budget. The workspace
    /// may be reused with the same state shape and ideal solvent model; physical
    /// coefficients, geometry, ions and fixed charges remain owned by the pair.
    pub fn adaptive_workspace(
        &self,
        maximum_workspace_bytes: usize,
    ) -> Result<DonnanAdvanceWorkspace> {
        let bytes = core::mem::size_of::<DonnanAdvanceWorkspace>()
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
        Ok(DonnanAdvanceWorkspace {
            original: State::new(self),
            start: State::new(self),
            coarse: State::new(self),
            species_count: self.mobile_totals_mol.len(),
            bytes,
        })
    }

    /// Advance exactly to the representable `time_s() + duration_s`. The local
    /// estimator compares the full step with two actual representable half
    /// intervals, accounting for unequal halves near clock resolution. Any
    /// failure restores all original inventories, equilibrium diagnostics and
    /// time. Prepared successful advancement performs no heap allocation.
    pub fn advance(
        &mut self,
        duration_s: f64,
        controls: DonnanAdvanceControls,
        workspace: &mut DonnanAdvanceWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DonnanAdvanceReport> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if !positive(duration_s) || !workspace.compatible(self) {
            return Err(invalid(
                "invalid adaptive Donnan duration or incompatible workspace",
            ));
        }
        let target = finite(
            self.time_s + duration_s,
            "Donnan target time is not representable",
        )?;
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
        controls: DonnanAdvanceControls,
        workspace: &mut DonnanAdvanceWorkspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DonnanAdvanceReport> {
        let start_time = self.time_s;
        let duration = target - start_time;
        let mut step = controls.maximum_step_s.min(duration);
        let (mut accepted, mut rejected, mut evaluations, mut roots, mut ion_roots) =
            (0, 0, 0usize, 0, 0);
        let (mut maximum_error, mut minimum_step, mut maximum_step): (f64, f64, f64) =
            (0.0, f64::INFINITY, 0.0);
        let (mut pressure_work, mut work_correction) = (0.0, 0.0);
        let mut endpoint = None;
        while self.time_s < target {
            poll(cancelled)?;
            if accepted >= controls.maximum_steps {
                return Err(OsmoticNetworkError::WorkBudget.into());
            }
            workspace.start.capture(self)?;
            let planned = if step >= target - self.time_s {
                target
            } else {
                self.time_s + step
            };
            // Solve a one-ULP terminal remainder as part of this physical step;
            // never merely move the clock to skip unintegrated transport.
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
            let coarse_controls = controls.inner(dt / duration)?;
            let first_controls = controls.inner(first_dt / duration)?;
            let second_controls = controls.inner(second_dt / duration)?;
            let (u, v) = (first_dt / dt, second_dt / dt);
            let factor = (u * u + v * v) / (2.0 * u * v);
            let trial = (|| -> Result<_> {
                evaluations = evaluations
                    .checked_add(1)
                    .ok_or(OsmoticNetworkError::WorkBudget)?;
                self.step_counted(
                    dt,
                    coarse_controls,
                    &mut roots,
                    &mut ion_roots,
                    &mut *cancelled,
                )?;
                workspace.coarse.capture(self)?;
                workspace.start.restore(self)?;
                evaluations = evaluations
                    .checked_add(1)
                    .ok_or(OsmoticNetworkError::WorkBudget)?;
                let first = self.step_counted(
                    first_dt,
                    first_controls,
                    &mut roots,
                    &mut ion_roots,
                    &mut *cancelled,
                )?;
                if self.time_s != midpoint {
                    return Err(OsmoticNetworkError::StepUnderflow.into());
                }
                evaluations = evaluations
                    .checked_add(1)
                    .ok_or(OsmoticNetworkError::WorkBudget)?;
                let second = self.step_counted(
                    second_dt,
                    second_controls,
                    &mut roots,
                    &mut ion_roots,
                    &mut *cancelled,
                )?;
                if self.time_s != end {
                    return Err(OsmoticNetworkError::StepUnderflow.into());
                }
                let error =
                    temporal_error(self, &workspace.start, &workspace.coarse, controls)? * factor;
                Ok((error, first.pressure_work_into_solution_j, second))
            })();
            let (error, first_work, second) = match trial {
                Ok(value) => value,
                Err(error) if retryable(&error) => {
                    workspace.start.restore(self)?;
                    if rejected >= controls.maximum_rejections || 0.5 * dt < controls.minimum_step_s
                    {
                        return Err(error);
                    }
                    rejected += 1;
                    step = 0.5 * dt;
                    continue;
                }
                Err(error) => return Err(error),
            };
            if error > 1.0 {
                workspace.start.restore(self)?;
                if 0.5 * dt < controls.minimum_step_s {
                    return Err(DonnanWaterError::TemporalAccuracy {
                        step_s: dt,
                        scaled_error: error,
                    });
                }
                if rejected >= controls.maximum_rejections {
                    return Err(OsmoticNetworkError::WorkBudget.into());
                }
                rejected += 1;
                step = 0.5 * dt;
                continue;
            }
            accepted += 1;
            maximum_error = maximum_error.max(error);
            minimum_step = minimum_step.min(dt);
            maximum_step = maximum_step.max(dt);
            compensated_add(&mut pressure_work, &mut work_correction, first_work);
            compensated_add(
                &mut pressure_work,
                &mut work_correction,
                second.pressure_work_into_solution_j,
            );
            endpoint = Some(second);
            step = if error < 0.125 {
                (2.0 * dt).min(controls.maximum_step_s)
            } else {
                dt
            };
        }
        let original = &workspace.original;
        let delta = [0, 1].map(|side| {
            self.compartments[side].solution().solvent_volume_m3()
                - original.compartments[side].solution().solvent_volume_m3()
        });
        let water_residual = finite(
            delta[0] + delta[1],
            "nonfinite adaptive Donnan water balance",
        )?;
        let total_water = original
            .compartments
            .iter()
            .map(|c| c.solution().solvent_volume_m3())
            .sum::<f64>();
        if water_residual.abs()
            > controls.step.absolute_water_tolerance_m3 + 64.0 * f64::EPSILON * total_water
        {
            return Err(invalid("adaptive Donnan water inventory balance failed"));
        }
        let mut balances = [0.0; MAX_DONNAN_SPECIES];
        for (i, residual) in balances
            .iter_mut()
            .enumerate()
            .take(self.mobile_totals_mol.len())
        {
            *residual = (self.compartments[0].solution().solute_moles()[i]
                - original.compartments[0].solution().solute_moles()[i])
                + (self.compartments[1].solution().solute_moles()[i]
                    - original.compartments[1].solution().solute_moles()[i]);
            if !residual.is_finite()
                || residual.abs() > 64.0 * f64::EPSILON * self.mobile_totals_mol[i]
            {
                return Err(DonnanWaterError::InventoryBalance {
                    species: i,
                    residual_mol: *residual,
                });
            }
        }
        let report = DonnanAdvanceReport {
            interval_s: [start_time, target],
            accepted_steps: accepted,
            rejected_steps: rejected,
            step_evaluations: evaluations,
            root_evaluations: roots,
            ion_root_evaluations: ion_roots,
            maximum_scaled_error: maximum_error,
            minimum_accepted_step_s: minimum_step,
            maximum_accepted_step_s: maximum_step,
            solvent_transfer_m3: -delta[0],
            pressure_work_into_solution_j: finite(pressure_work, "nonfinite adaptive Donnan work")?,
            water_balance_residual_m3: water_residual,
            solute_balance_residual_mol: balances,
            species_count: self.mobile_totals_mol.len(),
            endpoint: endpoint.ok_or(OsmoticNetworkError::StepUnderflow)?,
        };
        poll(cancelled)?;
        Ok(report)
    }
}
