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
        self.apply_impulse_density_with_receipt(impulses, controls, gate, Ok)
    }

    /// Admit a dependent receipt before publishing any original populations.
    /// Physical wrappers use this for fallible unit conversion without cloning
    /// the liquid owner or permitting a failure after accepted state changes.
    pub(crate) fn apply_impulse_density_with_receipt<R>(
        &mut self,
        impulses: &[[f64; 2]],
        controls: ImpulseControls,
        gate: &CancelGate,
        admit: impl FnOnce(ImpulseReport<2>) -> Result<R, ImpulseError>,
    ) -> Result<R, ImpulseError> {
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
        let changed = report.affected_cells != 0;
        let report = admit(report)?;
        poll(gate)?;
        if changed {
            core::mem::swap(&mut self.grid.f, &mut self.scratch);
        }
        Ok(report)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ContactModel, Grid};

    #[test]
    fn dependent_receipt_refusal_does_not_publish_native_populations() {
        let mut liquid = FreeSurface::new(Grid::uniform(2, 2, 0.9), 0.0, ContactModel::Neutral);
        let before = liquid.clone();
        let controls = ImpulseControls {
            max_cells: 4,
            maximum_speed: 0.1,
            maximum_momentum_defect: 1e-12,
            maximum_work_defect: 1e-12,
        };
        let gate = CancelGate::new_clock_free();
        let result = liquid.apply_impulse_density_with_receipt(
            &[[0.001, 0.0]; 4],
            controls,
            &gate,
            |report| {
                assert_eq!(report.affected_cells, 4);
                assert!(report.work > 0.0);
                Err::<(), _>(ImpulseError::Invalid("dependent SI receipt refused"))
            },
        );
        assert_eq!(result, Err(ImpulseError::Invalid("dependent SI receipt refused")));
        assert_eq!(liquid.grid.f, before.grid.f);
        assert_eq!(liquid.mass, before.mass);
        assert_eq!(liquid.grid.flags, before.grid.flags);
        assert_eq!(liquid.carry, before.carry);
        assert!(liquid.apply_impulse_density(&[[0.001, 0.0]; 4], controls, &gate).is_ok());
    }
}
