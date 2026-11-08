//! Instantaneous body impulse on fully occupied native D2Q9 fluid cells.

use super::{Cell, Grid};
use crate::{
    E, Q, W,
    impulse::{
        ImpulseControls, ImpulseError, ImpulseReport, ImpulseSum, WrenchImpulse2,
        WrenchImpulseControls2, WrenchImpulseReport2, poll,
        wrench::{WrenchCell2, WrenchPlan2, planar_momentum_change},
    },
};
use fs_exec::CancelGate;

impl Grid {
    /// Apply a prescribed impulse density in lattice momentum units, once.
    ///
    /// `impulses` covers the original x-fastest grid. Nonzero entries are
    /// admitted only at `Cell::Fluid`: wall/gas cells have no fluid inertia,
    /// and interface mass belongs to the free-surface owner. This method must
    /// not be used as a substitute for that owner's partial-volume accounting.
    /// Density and physical nonequilibrium stress survive; physical velocity
    /// changes by `impulse/rho`. Gravity, external forcing, relaxation times,
    /// flags and boundary data are retained. No time or collision is advanced.
    /// The returned work is signed kinetic work, never an inferred heat source.
    ///
    /// Reuses the caller's collision scratch. Errors/cancellation may change
    /// scratch but never accepted grid state. Zero impulse preserves all bits.
    pub fn apply_impulse_density(
        &mut self,
        impulses: &[[f64; 2]],
        controls: ImpulseControls,
        scratch: &mut Vec<[f64; Q]>,
        gate: &CancelGate,
    ) -> Result<ImpulseReport<2>, ImpulseError> {
        poll(gate)?;
        if impulses.len() != self.f.len() {
            return Err(ImpulseError::Invalid(
                "impulse must cover the original grid",
            ));
        }
        self.apply_impulse_with_receipt(
            controls,
            scratch,
            gate,
            |index| impulses[index],
            |report, _, _| Ok(report),
        )
    }

    /// Deposit a total planar force-and-couple impulse on selected fluid cells.
    ///
    /// Strictly increasing original x-fastest indices select actual population
    /// mass at cell centers `(x+1/2,y+1/2)`. Translation plus rigid rotation is
    /// the least mass-weighted velocity increment realizing the prescribed
    /// linear and angular impulse. One cell can receive a linear impulse, but
    /// cannot carry an independent couple. Nonperiodic topology is required;
    /// this method does not invent an unwrapped angular coordinate chart.
    ///
    /// Density, physical viscous stress and the original Guo force frame are
    /// retained. Actual linear/angular momentum and signed work are checked
    /// before publication. No time, collision, heat or implicit body reaction
    /// is added. Reuses the caller's population scratch and one compact cell
    /// stencil, without allocating a full-grid impulse array.
    pub fn apply_wrench_impulse(
        &mut self,
        indices: &[usize],
        wrench: WrenchImpulse2,
        controls: WrenchImpulseControls2,
        scratch: &mut Vec<[f64; Q]>,
        gate: &CancelGate,
    ) -> Result<WrenchImpulseReport2, ImpulseError> {
        poll(gate)?;
        let count = self
            .nx
            .checked_mul(self.ny)
            .ok_or(ImpulseError::Invalid("invalid grid size"))?;
        controls.impulse.admit(count)?;
        if self.periodic_x
            || self.periodic_y
            || indices.is_empty()
            || indices.len() > count
            || self.f.len() != count
            || self.flags.len() != count
            || self.fext.len() != count
        {
            return Err(ImpulseError::Invalid(
                "wrench requires complete nonperiodic fields and a nonempty bounded fluid stencil",
            ));
        }
        let mut cells = Vec::new();
        cells
            .try_reserve_exact(indices.len())
            .map_err(|_| ImpulseError::Invalid("wrench stencil allocation refused"))?;
        for (i, &index) in indices.iter().enumerate() {
            if i % 64 == 0 {
                poll(gate)?;
            }
            if index >= count || self.flags[index] != Cell::Fluid {
                return Err(ImpulseError::Invalid(
                    "wrench stencil contains a nonfluid cell",
                ));
            }
            cells.push(WrenchCell2 {
                index,
                position: [
                    (index % self.nx) as f64 + 0.5,
                    (index / self.nx) as f64 + 0.5,
                ],
                mass: self.f[index].iter().sum(),
            });
        }
        let plan = WrenchPlan2::prepare(&cells, wrench, controls, gate)?;
        let mut next_cell = 0;
        self.apply_impulse_with_receipt(
            controls.impulse,
            scratch,
            gate,
            |index| {
                if next_cell < cells.len() && cells[next_cell].index == index {
                    let impulse = plan.impulse(cells[next_cell]);
                    next_cell += 1;
                    impulse
                } else {
                    [0.0; 2]
                }
            },
            |report, before, after| {
                plan.finish(&cells, report, controls, gate, |index| {
                    planar_momentum_change(
                        &before.f[index],
                        &after[index],
                        before.g,
                        before.fext[index],
                        None,
                    )
                })
            },
        )
    }

    fn apply_impulse_with_receipt<R>(
        &mut self,
        controls: ImpulseControls,
        scratch: &mut Vec<[f64; Q]>,
        gate: &CancelGate,
        mut impulse_at: impl FnMut(usize) -> [f64; 2],
        admit: impl FnOnce(ImpulseReport<2>, &Self, &[[f64; Q]]) -> Result<R, ImpulseError>,
    ) -> Result<R, ImpulseError> {
        poll(gate)?;
        let count = self
            .nx
            .checked_mul(self.ny)
            .ok_or(ImpulseError::Invalid("invalid grid size"))?;
        controls.admit(count)?;
        if self.f.len() != count || self.flags.len() != count || self.fext.len() != count {
            return Err(ImpulseError::Invalid(
                "impulse and original fields must cover the grid",
            ));
        }
        scratch.clear();
        scratch
            .try_reserve(count)
            .map_err(|_| ImpulseError::Invalid("impulse scratch allocation refused"))?;
        let directions = E.map(|(x, y)| [f64::from(x), f64::from(y)]);
        let mut sum = ImpulseSum::new();
        for (index, &original) in self.f.iter().enumerate() {
            if index % 64 == 0 {
                poll(gate)?;
            }
            let impulse = impulse_at(index);
            if self.flags[index] != Cell::Fluid && impulse != [0.0; 2] {
                return Err(ImpulseError::Invalid(
                    "nonzero impulse requires a fully occupied fluid cell",
                ));
            }
            let mut next = original;
            sum.kick(
                &mut next,
                &directions,
                &W,
                self.g,
                self.fext[index],
                impulse,
                controls,
                index,
            )?;
            scratch.push(next);
        }
        let report = sum.finish(controls)?;
        let changed = report.affected_cells != 0;
        let report = admit(report, self, scratch)?;
        poll(gate)?;
        if changed {
            core::mem::swap(&mut self.f, scratch);
        }
        Ok(report)
    }
}
