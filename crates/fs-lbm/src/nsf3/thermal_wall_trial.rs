//! One reversible native step with exact local wall heat exchange.

use super::{EnergyLedger, NsfError, NsfGrid3, NsfThermalWallFace, body_force};
use crate::lattice3::Lattice3;

/// Heat on one authored impermeable face over the original lattice step.
/// Positive heat enters the fluid for either orientation. This is the same
/// old/new molecular flux pair used in its conservative energy equation;
/// moving-wall stress work and unlisted global walls remain separate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfThermalWallExchange {
    pub cell: [usize; 3],
    pub axis: usize,
    pub high: bool,
    pub heat_into_fluid_j: f64,
}

/// An unpublished physical fluid step, including its exact per-face heat.
/// Dropping it preserves every accepted field, material cache, boundary
/// control, force/energy history, clock and cumulative account.
pub struct NsfThermalWallTrial<'a, L: Lattice3<27>> {
    owner: &'a mut NsfGrid3<L>,
    candidate: NsfGrid3<L>,
    ledger: EnergyLedger,
    exchange: Vec<NsfThermalWallExchange>,
    previous_energy_j: f64,
}

impl<L: Lattice3<27>> NsfThermalWallTrial<'_, L> {
    pub const fn grid(&self) -> &NsfGrid3<L> {
        &self.candidate
    }

    /// The original complete one-step account, with final storage in joules.
    pub const fn energy_ledger(&self) -> &EnergyLedger {
        &self.ledger
    }

    /// Exact native trapezoids in the canonical thermal-wall table order.
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
    /// The controls are held over this step; the original setter refreshes its
    /// starting FV history even when the proposed table is unchanged, and its
    /// actual old/new face fluxes supply the heat
    /// account. No post-step material query reconstructs the exchanged joules.
    ///
    /// This supports closed domains, existing held openings, molecular, local
    /// and dynamic SGS transport, gravity and held body forces. The dynamic
    /// closures retain their original controls, filter geometry and accepted
    /// transport; only their deterministic scratch is rebuilt. The native wall
    /// update refits momentum/heat for the proposed thermal boundary before the
    /// step. Modeled heat remains on interior faces; molecular conduction owns
    /// the wall exchange. Manufactured `SourceFn` callbacks remain refused
    /// because their external state cannot participate in this transaction.
    /// The native source material, Mach, caloric, stability and endpoint
    /// convergence admissions still apply.
    ///
    /// `maximum_cells` bounds the full copied grid, including solid cells.
    /// Candidate storage and native step work are O(cells), plus the existing
    /// bounded corrector passes. A dynamic candidate additionally owns one
    /// filter workspace, bounded by the original retained workspace bytes.
    /// Cancellation is polled before and after the
    /// copy and wall update, at native work boundaries, and before return.
    /// Callers coupling a finite wall retain their own time-resolution and
    /// joint energy acceptance; this is not a separate fluid integrator.
    pub fn preview_thermal_wall_step(
        &mut self,
        walls: Vec<NsfThermalWallFace>,
        maximum_cells: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NsfThermalWallTrial<'_, L>, NsfError> {
        body_force::checkpoint(&mut cancelled)?;
        if self.n > maximum_cells {
            return Err(NsfError::Cap {
                requested: self.n,
                maximum: maximum_cells,
            });
        }
        if self.sources.is_some() {
            return Err(NsfError::Input(
                "thermal wall trials require fluid state without SourceFn callbacks",
            ));
        }
        if self.steps.checked_add(1).is_none() {
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
        // that owner's controls and independent scratch explicitly for walls.
        let mut candidate = self.characteristic_copy();
        body_force::checkpoint(&mut cancelled)?;
        candidate.dynamic_sgs = self
            .dynamic_sgs
            .as_ref()
            .map(|dynamic| dynamic.fresh_for_trial(self.spec.dims))
            .transpose()?;
        body_force::checkpoint(&mut cancelled)?;
        let before = candidate.set_thermal_walls_recorded(walls, &mut cancelled)?;
        let mut exchange = Vec::new();
        let ledger = candidate.step_with_wall_exchange(
            None,
            &mut cancelled,
            Some((&before, &mut exchange)),
        )?;
        if !candidate.total_energy_j().is_finite() || !ledger.residual_rel.is_finite() {
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

#[cfg(test)]
mod tests;
