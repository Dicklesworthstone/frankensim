//! Reversible native intervals with exact local wall heat exchange.

use super::{EnergyLedger, NsfError, NsfGrid3, NsfThermalWallFace, body_force};
use crate::lattice3::Lattice3;

/// Heat on one authored impermeable face over original lattice steps.
/// Positive heat enters the fluid for either orientation. This is the same
/// sum of old/new molecular flux pairs used in its conservative energy equation;
/// moving-wall stress work and unlisted global walls remain separate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfThermalWallExchange {
    pub cell: [usize; 3],
    pub axis: usize,
    pub high: bool,
    pub heat_into_fluid_j: f64,
}

/// An unpublished physical fluid interval, including its exact per-face heat.
/// Dropping it preserves every accepted field, material cache, boundary
/// control, force/energy history, clock and cumulative account.
pub struct NsfThermalWallTrial<'a, L: Lattice3<27>> {
    owner: &'a mut NsfGrid3<L>,
    candidate: NsfGrid3<L>,
    ledger: EnergyLedger,
    exchange: Vec<NsfThermalWallExchange>,
    previous_energy_j: f64,
}

// Exact fluxes retained from the accepted conservative corrector. A multi-step
// interval uses next_flux as the next native step's previous_flux, without
// refreshing history or querying material at its corrected final temperature.
pub(super) struct ThermalWallStepRecord<'a> {
    pub(super) previous_flux: &'a [f64],
    pub(super) next_flux: &'a mut Vec<f64>,
    pub(super) exchanges: &'a mut Vec<NsfThermalWallExchange>,
}

impl<L: Lattice3<27>> NsfThermalWallTrial<'_, L> {
    pub const fn grid(&self) -> &NsfGrid3<L> {
        &self.candidate
    }

    /// Summed original step accounts, with final storage in joules. Multi-step
    /// residuals close this whole interval; one step retains its native ledger.
    pub const fn energy_ledger(&self) -> &EnergyLedger {
        &self.ledger
    }

    /// Summed native trapezoids in the canonical thermal-wall table order.
    /// These cover only explicitly listed faces; their sum need not equal
    /// aggregate wall heat when unlisted global walls also exchange heat.
    pub fn wall_exchange(&self) -> &[NsfThermalWallExchange] {
        &self.exchange
    }

    pub const fn previous_energy_j(&self) -> f64 {
        self.previous_energy_j
    }

    /// Publish this already admitted candidate. No fallible work follows.
    pub fn accept(self) {
        *self.owner = self.candidate;
    }
}

impl<L: Lattice3<27>> NsfGrid3<L> {
    /// Preview one original lattice step under proposed local thermal walls.
    /// The controls are held over this step; the recorded setter (`record = true`
    /// in `set_thermal_walls_impl`) refreshes its starting FV history even when
    /// the proposed table is unchanged, which the public setter does not (it
    /// returns early then), and its
    /// actual old/new face fluxes supply the heat
    /// account. No post-step material query reconstructs the exchanged joules.
    ///
    /// Equivalent to [`Self::preview_thermal_wall_interval`] with one step and
    /// `maximum_cell_steps = maximum_cells`. It retains the original one-step
    /// ledger and numerical path, including the explicit starting-history rebase.
    pub fn preview_thermal_wall_step(
        &mut self,
        walls: Vec<NsfThermalWallFace>,
        maximum_cells: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<NsfThermalWallTrial<'_, L>, NsfError> {
        self.preview_thermal_wall_interval(walls, 1, maximum_cells, maximum_cells, cancelled)
    }

    /// Preview a positive, bounded number of original lattice steps under one
    /// held local wall table. Copy physical state and refresh its starting wall
    /// history once. Every later step consumes the previous accepted corrector's
    /// exact sparse face fluxes; no intermediate material requery or history
    /// reset changes the original trajectory or its consumed heat quadrature.
    ///
    /// Closed domains, existing held openings, molecular/local/dynamic SGS,
    /// gravity and held body forces retain their original owners. Manufactured
    /// `SourceFn` callbacks are refused because their external mutable state
    /// cannot be rolled back. The original source material, Mach, caloric,
    /// stability and endpoint convergence admissions apply to every step.
    ///
    /// Periodic/channel dynamic filters share only their already-budgeted,
    /// mutex-protected scratch with this exclusive trial. Every fit prepares
    /// the current candidate fields before consuming that scratch. Lagrangian
    /// filters instead copy accepted momentum/heat history into an independent
    /// workspace, bounded by the copied cell count and the original per-owner
    /// workspace admission; they add 264 nominal bytes per candidate cell.
    /// Candidate native steps advance only that candidate's history. Accepted
    /// transport coefficients, controls and reports remain separate. No second
    /// fluid solver is introduced. A dropped/refused
    /// trial may leave scratch overwritten, as an ordinary failed native step
    /// can, but cannot publish any candidate physical/dynamic coefficient.
    /// The original wall update refits momentum/heat for the proposed boundary;
    /// modeled heat remains on interior faces and molecular conduction alone
    /// owns the wall exchange.
    ///
    /// `maximum_cells` bounds the full copied grid, including solid cells.
    /// `maximum_cell_steps` separately bounds `cell_count() * steps`, checked
    /// before copying. Candidate storage is O(cells); work is O(cells * steps)
    /// times the existing bounded corrector passes. Cancellation is polled
    /// before/after the copy and wall update, between original steps, at native
    /// work boundaries and while accumulating faces. A late refusal discards
    /// the entire interval, including a successfully computed prefix.
    /// Callers coupling a finite wall retain their own time-resolution and
    /// joint energy acceptance; this is not a separate fluid integrator.
    pub fn preview_thermal_wall_interval(
        &mut self,
        walls: Vec<NsfThermalWallFace>,
        steps: usize,
        maximum_cells: usize,
        maximum_cell_steps: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NsfThermalWallTrial<'_, L>, NsfError> {
        body_force::checkpoint(&mut cancelled)?;
        if self.n > maximum_cells {
            return Err(NsfError::Cap {
                requested: self.n,
                maximum: maximum_cells,
            });
        }
        if steps == 0 {
            return Err(NsfError::Input(
                "thermal wall interval needs positive native steps",
            ));
        }
        let cell_steps = self.n.checked_mul(steps).ok_or(NsfError::Input(
            "thermal wall interval cell-step work overflow",
        ))?;
        if cell_steps > maximum_cell_steps {
            return Err(NsfError::Cap {
                requested: cell_steps,
                maximum: maximum_cell_steps,
            });
        }
        if self.sources.is_some() {
            return Err(NsfError::Input(
                "thermal wall trials cannot retain manufactured source callbacks",
            ));
        }
        if self.steps.checked_add(steps).is_none() {
            return Err(NsfError::Input("thermal wall trial clock overflow"));
        }
        let previous_energy_j = self.total_energy_j();
        if !previous_energy_j.is_finite() {
            return Err(NsfError::Input(
                "thermal wall trial starts with nonfinite energy",
            ));
        }
        // Retain every accepted field with the existing exact native copy.
        // Characteristic topology cannot carry a dynamic closure, so restore
        // its original controls and shared scratch explicitly for wall trials.
        let mut candidate = self.characteristic_copy();
        body_force::checkpoint(&mut cancelled)?;
        candidate.dynamic_sgs = self
            .dynamic_sgs
            .as_ref()
            .map(|dynamic| dynamic.shared_for_trial())
            .transpose()?;
        body_force::checkpoint(&mut cancelled)?;
        let mut before = candidate.set_thermal_walls_recorded(walls, &mut cancelled)?;
        let mut after = Vec::new();
        let mut step_exchange = Vec::new();
        let mut exchange = Vec::new();
        let mut ledger = EnergyLedger::default();
        for step_index in 0..steps {
            body_force::checkpoint(&mut cancelled)?;
            let step = candidate.step_with_wall_exchange(
                None,
                &mut cancelled,
                Some(ThermalWallStepRecord {
                    previous_flux: &before,
                    next_flux: &mut after,
                    exchanges: &mut step_exchange,
                }),
            )?;
            before = core::mem::take(&mut after);
            if step_index == 0 {
                // Preserve the original one-step return exactly, including
                // its residual normalization and signed zero arithmetic.
                ledger = step;
                exchange = core::mem::take(&mut step_exchange);
            } else {
                accumulate_ledger(&mut ledger, step);
                if exchange.len() != step_exchange.len() {
                    return Err(NsfError::Input("thermal wall interval face count changed"));
                }
                for (total, current) in exchange.iter_mut().zip(&step_exchange) {
                    body_force::checkpoint(&mut cancelled)?;
                    if (total.cell, total.axis, total.high)
                        != (current.cell, current.axis, current.high)
                    {
                        return Err(NsfError::Input(
                            "thermal wall interval face identity changed",
                        ));
                    }
                    total.heat_into_fluid_j += current.heat_into_fluid_j;
                    if !total.heat_into_fluid_j.is_finite() {
                        return Err(NsfError::Input(
                            "thermal wall interval face heat overflowed",
                        ));
                    }
                }
            }
            if !finite_ledger(ledger) {
                return Err(NsfError::Input(
                    "thermal wall interval energy account overflowed",
                ));
            }
        }
        let energy = candidate.total_energy_j();
        if steps > 1 {
            ledger.residual_rel = (energy
                - previous_energy_j
                - ledger.wall_heat_j
                - ledger.wall_work_j
                - ledger.open_boundary_energy_j
                - ledger.open_energy_j
                - ledger.source_energy_j
                - ledger.gravity_work_j
                - ledger.body_force_work_j)
                .abs()
                / previous_energy_j.abs().max(energy.abs()).max(1e-300);
        }
        if !energy.is_finite() || !ledger.residual_rel.is_finite() {
            return Err(NsfError::Input(
                "thermal wall trial produced nonfinite energy",
            ));
        }
        body_force::checkpoint(&mut cancelled)?;
        Ok(NsfThermalWallTrial {
            owner: self,
            candidate,
            ledger,
            exchange,
            previous_energy_j,
        })
    }
}

fn accumulate_ledger(total: &mut EnergyLedger, step: EnergyLedger) {
    total.internal_j = step.internal_j;
    total.kinetic_j = step.kinetic_j;
    total.viscous_heating_j += step.viscous_heating_j;
    total.unresolved_internal_exchange_j += step.unresolved_internal_exchange_j;
    total.pressure_work_j += step.pressure_work_j;
    total.wall_heat_j += step.wall_heat_j;
    total.wall_work_j += step.wall_work_j;
    total.open_boundary_energy_j += step.open_boundary_energy_j;
    total.open_energy_j += step.open_energy_j;
    total.source_energy_j += step.source_energy_j;
    total.gravity_work_j += step.gravity_work_j;
    total.body_force_work_j += step.body_force_work_j;
}

fn finite_ledger(ledger: EnergyLedger) -> bool {
    [
        ledger.internal_j,
        ledger.kinetic_j,
        ledger.viscous_heating_j,
        ledger.unresolved_internal_exchange_j,
        ledger.pressure_work_j,
        ledger.wall_heat_j,
        ledger.wall_work_j,
        ledger.open_boundary_energy_j,
        ledger.open_energy_j,
        ledger.source_energy_j,
        ledger.gravity_work_j,
        ledger.body_force_work_j,
        ledger.residual_rel,
    ]
    .into_iter()
    .all(f64::is_finite)
}

#[cfg(test)]
mod tests;
