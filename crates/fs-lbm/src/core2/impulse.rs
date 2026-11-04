//! Instantaneous body impulse on fully occupied native D2Q9 fluid cells.

use super::{Cell, Grid};
use crate::{
    E, Q, W,
    impulse::{ImpulseControls, ImpulseError, ImpulseReport, ImpulseSum, poll},
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
        let count = self
            .nx
            .checked_mul(self.ny)
            .ok_or(ImpulseError::Invalid("invalid grid size"))?;
        controls.admit(count)?;
        if self.f.len() != count
            || self.flags.len() != count
            || self.fext.len() != count
            || impulses.len() != count
        {
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
        for (index, (&original, &impulse)) in self.f.iter().zip(impulses).enumerate() {
            if index % 64 == 0 {
                poll(gate)?;
            }
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
        poll(gate)?;
        if report.affected_cells != 0 {
            core::mem::swap(&mut self.f, scratch);
        }
        Ok(report)
    }
}
