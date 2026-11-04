//! Instantaneous loading of the original dense liquid mass inventory.

use super::FreeSurface;
use crate::{
    Cell, E, W,
    impulse::{ImpulseControls, ImpulseError, ImpulseReport, ImpulseSum, poll},
};
use fs_exec::CancelGate;

impl FreeSurface {
    /// Apply one signed impulse per original cell, in lattice momentum units
    /// per whole geometric cell (per unit depth in this two-dimensional owner).
    ///
    /// A fluid cell has mass `sum(f)`; an interface cell uses its original
    /// tracked liquid mass. Thus `delta_u = impulse / liquid_mass`, including
    /// partial fill. Inputs are NOT automatically multiplied by fill: a source
    /// specified per occupied liquid volume must supply that factor itself.
    /// Gas and wall entries must be zero. Unlocated conversion carry is retained
    /// and cannot receive a cell-addressed impulse.
    ///
    /// The shared hydrodynamic kernel retains the Guo force frame, intrinsic
    /// density and viscous stress. Canonical density closes exactly before
    /// readback, so mass, fill, cell classification, capillarity and wetting
    /// history are unchanged. No time, streaming, conversion, damping or heat
    /// is added. This is a resolved low-Mach momentum input; pressure transport
    /// and interface motion remain with subsequent native evolution.
    ///
    /// Reuses the existing collision scratch. The whole proposal is checked
    /// for speed, actual liquid momentum, signed work and cancellation before
    /// publication. Invalid or unresolved partial mass refuses without clamping.
    /// Errors may change scratch, but never accepted liquid state.
    pub fn apply_impulse_density(
        &mut self,
        impulses: &[[f64; 2]],
        controls: ImpulseControls,
        gate: &CancelGate,
    ) -> Result<ImpulseReport<2>, ImpulseError> {
        poll(gate)?;
        let count = self
            .grid
            .nx
            .checked_mul(self.grid.ny)
            .ok_or(ImpulseError::Invalid("invalid liquid grid size"))?;
        controls.admit(count)?;
        if self.grid.f.len() != count
            || self.grid.flags.len() != count
            || self.grid.fext.len() != count
            || self.mass.len() != count
            || impulses.len() != count
        {
            return Err(ImpulseError::Invalid(
                "impulse and liquid fields must cover the original grid",
            ));
        }
        self.scratch.clear();
        self.scratch
            .try_reserve(count)
            .map_err(|_| ImpulseError::Invalid("liquid impulse scratch allocation refused"))?;
        let directions = E.map(|(x, y)| [f64::from(x), f64::from(y)]);
        let mut sum = ImpulseSum::new();
        for (index, (&original, &impulse)) in self.grid.f.iter().zip(impulses).enumerate() {
            if index % 64 == 0 {
                poll(gate)?;
            }
            let mut next = original;
            if impulse != [0.0; 2] {
                let mass = match self.grid.flags[index] {
                    Cell::Fluid => original.iter().sum(),
                    Cell::Interface => self.mass[index],
                    Cell::Gas | Cell::Wall => {
                        return Err(ImpulseError::Invalid(
                            "nonzero impulse requires original liquid mass",
                        ));
                    }
                };
                sum.kick_with_mass(
                    &mut next,
                    &directions,
                    &W,
                    self.grid.g,
                    self.grid.fext[index],
                    impulse,
                    mass,
                    controls,
                    index,
                )?;
            }
            self.scratch.push(next);
        }
        let report = sum.finish(controls)?;
        poll(gate)?;
        if report.affected_cells != 0 {
            core::mem::swap(&mut self.grid.f, &mut self.scratch);
        }
        Ok(report)
    }
}
