//! Compact, borrowed impulse staging for reciprocal liquid/structure sources.
//!
//! The borrow prevents fluid mutation between sample, population readback and
//! publication. Dropping a stage or candidate leaves the original owner intact.

use super::super::{Cell3, FreeSurface3};
use crate::d3q19::{E3, W3};
use crate::impulse::{ImpulseControls, ImpulseError, ImpulseReport, ImpulseSum, poll};
use fs_exec::CancelGate;

/// Original located liquid state, in lattice units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceImpulseCell3 {
    /// X-fastest index in the complete fixed domain.
    pub index: usize,
    /// Cell centre measured from the original domain corner.
    pub position: [f64; 3],
    /// Actual retained liquid mass, not full-cell population density.
    pub mass: f64,
    /// Original physical velocity, including the held Guo half-force.
    pub velocity: [f64; 3],
}

/// Actual population readback at the same original position and mass.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceImpulseReadback3 {
    pub before: FreeSurfaceImpulseCell3,
    pub velocity_after: [f64; 3],
    pub momentum_change: [f64; 3],
}

struct Cell {
    sample: FreeSurfaceImpulseCell3,
    slot: usize,
    lane: usize,
    populations: [f64; 19],
    force: [f64; 3],
    velocity_after: [f64; 3],
}

/// A borrowed original wet stencil; no population has changed.
#[must_use]
pub struct FreeSurfaceImpulseStencil3<'a> {
    owner: &'a mut FreeSurface3,
    cells: Vec<Cell>,
    controls: ImpulseControls,
}

/// Complete compact candidate, ready for a coupled physical acceptance gate.
/// The original owner is borrowed and unchanged until commit succeeds.
#[must_use]
pub struct FreeSurfaceImpulseCandidate3<'a> {
    owner: &'a mut FreeSurface3,
    cells: Vec<Cell>,
    report: ImpulseReport<3>,
}

fn velocity(populations: &[f64; 19], force: [f64; 3]) -> [f64; 3] {
    let density = populations.iter().sum::<f64>();
    core::array::from_fn(|axis| {
        let raw = populations
            .iter()
            .zip(E3)
            .map(|(f, e)| f * [f64::from(e.0), f64::from(e.1), f64::from(e.2)][axis])
            .sum::<f64>();
        (raw + 0.5 * force[axis]) / density
    })
}

impl FreeSurface3 {
    /// Borrow a compact, strictly increasing stencil of original wet cells.
    ///
    /// Storage and work scale with selected cells, with native tile lookups;
    /// no dense impulse field or whole-grid clone is made. The original complete
    /// domain must fit controls.max_cells. Cells must carry located positive
    /// liquid inventory; gas, solids, piston cells and unresolved interface mass
    /// refuse. Unlocated carry has no invented position or inertia.
    ///
    /// The returned samples use the exact retained mass and physical Guo-frame
    /// velocity. Stage impulses from those samples, inspect their independent
    /// readback, then publish only after the other participant is also ready.
    pub fn prepare_impulse_stencil(
        &mut self,
        indices: &[usize],
        controls: ImpulseControls,
        gate: &CancelGate,
    ) -> Result<FreeSurfaceImpulseStencil3<'_>, ImpulseError> {
        poll(gate)?;
        let [nx, ny, nz] = self.domain_cells();
        let count =
            nx.checked_mul(ny)
                .and_then(|n| n.checked_mul(nz))
                .ok_or(ImpulseError::Invalid(
                    "liquid impulse domain size overflows",
                ))?;
        controls.admit(count)?;
        if indices.is_empty() || indices.len() > count {
            return Err(ImpulseError::Invalid(
                "liquid impulse stencil must be nonempty and bounded",
            ));
        }
        let mut cells = Vec::new();
        cells
            .try_reserve_exact(indices.len())
            .map_err(|_| ImpulseError::Invalid("liquid impulse stencil allocation refused"))?;
        for (i, &index) in indices.iter().enumerate() {
            if i % 64 == 0 {
                poll(gate)?;
            }
            if index >= count || (i != 0 && indices[i - 1] >= index) {
                return Err(ImpulseError::Invalid(
                    "liquid impulse indices must be in-domain and strictly increasing",
                ));
            }
            let [x, y, z] = [index % nx, (index / nx) % ny, index / nx / ny];
            let position = [x as i64, y as i64, z as i64];
            if self.planar_piston_contains(position) {
                return Err(ImpulseError::Invalid(
                    "liquid impulse stencil intersects a piston",
                ));
            }
            let (slot, lane) = self
                .grid
                .resolve_source(position[0], position[1], position[2])
                .ok_or(ImpulseError::Invalid(
                    "liquid impulse stencil is not located in liquid",
                ))?;
            let key = self.grid.active_keys()[slot];
            let tile = &self.tiles[&key];
            let populations = self.grid.populations(slot, lane);
            let density = populations.iter().sum::<f64>();
            let mass = match tile.cells[lane] {
                Cell3::Fluid => density,
                Cell3::Interface => tile.mass[lane],
                Cell3::Gas => {
                    return Err(ImpulseError::Invalid("liquid impulse stencil contains gas"));
                }
            };
            let force = self.force_at_lane(key, lane);
            let velocity = velocity(&populations, force);
            let speed = velocity[0].hypot(velocity[1]).hypot(velocity[2]);
            if !mass.is_finite()
                || mass <= 0.0
                || !density.is_finite()
                || mass > density
                || populations.iter().any(|f| !f.is_finite() || *f <= 0.0)
                || force.iter().any(|f| !f.is_finite())
                || !speed.is_finite()
                || speed > controls.maximum_speed
            {
                return Err(ImpulseError::Cell { index });
            }
            cells.push(Cell {
                sample: FreeSurfaceImpulseCell3 {
                    index,
                    position: [x as f64 + 0.5, y as f64 + 0.5, z as f64 + 0.5],
                    mass,
                    velocity,
                },
                slot,
                lane,
                populations,
                force,
                velocity_after: velocity,
            });
        }
        poll(gate)?;
        Ok(FreeSurfaceImpulseStencil3 {
            owner: self,
            cells,
            controls,
        })
    }
}

impl<'a> FreeSurfaceImpulseStencil3<'a> {
    /// Original samples in the caller's canonical increasing index order.
    pub fn cells(&self) -> impl ExactSizeIterator<Item = FreeSurfaceImpulseCell3> + '_ {
        self.cells.iter().map(|cell| cell.sample)
    }

    /// Stage one total liquid impulse per selected cell, in the same order.
    ///
    /// The original population kernel uses retained liquid mass, preserves
    /// density and physical nonequilibrium stress, and retains the held force.
    /// No collision, transport, fill/classification, geometry or clock advances.
    /// An error drops only staged storage. A zero field retains original bits.
    pub fn with_impulses(
        mut self,
        impulses: &[[f64; 3]],
        gate: &CancelGate,
    ) -> Result<FreeSurfaceImpulseCandidate3<'a>, ImpulseError> {
        poll(gate)?;
        if impulses.len() != self.cells.len() {
            return Err(ImpulseError::Invalid(
                "compact liquid impulse field has the wrong length",
            ));
        }
        let directions = E3.map(|(x, y, z)| [f64::from(x), f64::from(y), f64::from(z)]);
        let mut sum = ImpulseSum::new();
        for (i, (cell, &impulse)) in self.cells.iter_mut().zip(impulses).enumerate() {
            if i % 64 == 0 {
                poll(gate)?;
            }
            sum.kick_with_mass(
                &mut cell.populations,
                &directions,
                &W3,
                [0.0; 3],
                cell.force,
                impulse,
                cell.sample.mass,
                self.controls,
                cell.sample.index,
            )?;
            cell.velocity_after = velocity(&cell.populations, cell.force);
        }
        let report = sum.finish(self.controls)?;
        poll(gate)?;
        Ok(FreeSurfaceImpulseCandidate3 {
            owner: self.owner,
            cells: self.cells,
            report,
        })
    }
}

impl FreeSurfaceImpulseCandidate3<'_> {
    /// Original population-level linear momentum and signed-work receipt.
    pub const fn report(&self) -> &ImpulseReport<3> {
        &self.report
    }

    /// Actual original/staged physical velocities, using the SAME retained mass.
    /// Coupling owners can form torque about their own physical origin without
    /// inventing a body reaction or confusing empty interface volume with mass.
    pub fn readback(&self) -> impl ExactSizeIterator<Item = FreeSurfaceImpulseReadback3> + '_ {
        self.cells.iter().map(|cell| FreeSurfaceImpulseReadback3 {
            before: cell.sample,
            velocity_after: cell.velocity_after,
            momentum_change: core::array::from_fn(|axis| {
                cell.sample.mass * (cell.velocity_after[axis] - cell.sample.velocity[axis])
            }),
        })
    }

    /// Publish only the selected populations after the final cancellation poll.
    /// The whole candidate is already validated; no fallible work follows the
    /// first write. The composing owner publishes its prepared reaction too.
    pub fn commit(self, gate: &CancelGate) -> Result<ImpulseReport<3>, ImpulseError> {
        poll(gate)?;
        if self.report.affected_cells != 0 {
            for cell in self.cells {
                self.owner
                    .grid
                    .set_populations(cell.slot, cell.lane, cell.populations);
            }
        }
        Ok(self.report)
    }
}

#[cfg(test)]
mod tests;
