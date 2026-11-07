//! Compact force-and-couple loading of the original located liquid inventory.

use super::super::{Cell3, FreeSurface3};
use crate::d3q19::{E3, W3};
use crate::impulse::{
    ImpulseError, ImpulseSum, WrenchImpulse3, WrenchImpulseControls3, WrenchImpulseReport3, poll,
    wrench::{WrenchCell3, WrenchPlan3},
};
use fs_exec::CancelGate;

impl FreeSurface3 {
    /// Deposit one total linear/angular impulse on selected original wet cells.
    ///
    /// `indices` are strictly increasing x-fastest indices in the complete
    /// fixed domain. Positions are cell centers `(x+1/2,y+1/2,z+1/2)` in lattice
    /// lengths. The selected liquid's actual retained mass determines its
    /// center and inertia: full `Fluid` cells use population density, while
    /// `Interface` cells use tracked liquid mass, never the gas-filled volume.
    /// Unlocated redistribution carry supplies neither a position nor inertia.
    ///
    /// The shared rank-aware wrench solve chooses the least mass-weighted
    /// squared velocity change. A collinear cloud carries transverse couples;
    /// an unresolved axial couple refuses. Each population kick preserves
    /// intrinsic density and physical nonequilibrium stress in its held Guo
    /// force frame. Actual liquid momentum, angular momentum and signed work
    /// are admitted before any population is published.
    ///
    /// Uses two compact stencil vectors and bounded native tile lookups, with
    /// no whole-domain scan or impulse field. Only selected populations change;
    /// fill, mass, carry, classification, geometry, force and time stay fixed.
    /// This is neither a collision/contact law nor a pressure projection. A
    /// coupled solid must receive its opposite wrench from the coupling owner.
    /// Refusal or cancellation preserves the complete accepted liquid state.
    pub fn apply_wrench_impulse(
        &mut self,
        indices: &[usize],
        wrench: WrenchImpulse3,
        controls: WrenchImpulseControls3,
        gate: &CancelGate,
    ) -> Result<WrenchImpulseReport3, ImpulseError> {
        self.apply_wrench_impulse_with_receipt(indices, wrench, controls, gate, Ok)
    }

    /// Physical adapters admit their converted receipt inside this transaction.
    pub(crate) fn apply_wrench_impulse_with_receipt<R>(
        &mut self,
        indices: &[usize],
        wrench: WrenchImpulse3,
        controls: WrenchImpulseControls3,
        gate: &CancelGate,
        admit: impl FnOnce(WrenchImpulseReport3) -> Result<R, ImpulseError>,
    ) -> Result<R, ImpulseError> {
        poll(gate)?;
        let [nx, ny, nz] = self.domain_cells();
        let count = nx
            .checked_mul(ny)
            .and_then(|n| n.checked_mul(nz))
            .ok_or(ImpulseError::Invalid("invalid free-surface domain size"))?;
        controls.impulse.admit(count)?;
        if indices.is_empty() || indices.len() > count {
            return Err(ImpulseError::Invalid(
                "liquid wrench needs a nonempty bounded wet-cell stencil",
            ));
        }
        let mut cells = Vec::new();
        let mut replacements = Vec::new();
        cells
            .try_reserve_exact(indices.len())
            .map_err(|_| ImpulseError::Invalid("liquid wrench stencil allocation refused"))?;
        replacements
            .try_reserve_exact(indices.len())
            .map_err(|_| ImpulseError::Invalid("liquid wrench population allocation refused"))?;
        for (i, &index) in indices.iter().enumerate() {
            if i % 64 == 0 {
                poll(gate)?;
            }
            if index >= count || (i != 0 && indices[i - 1] >= index) {
                return Err(ImpulseError::Invalid(
                    "liquid wrench indices must be in-domain and strictly increasing",
                ));
            }
            let [x, y, z] = [index % nx, (index / nx) % ny, index / nx / ny];
            let position = [x as i64, y as i64, z as i64];
            if self.planar_piston_contains(position) {
                return Err(ImpulseError::Invalid(
                    "liquid wrench stencil intersects a piston",
                ));
            }
            let (slot, lane) = self
                .grid
                .resolve_source(position[0], position[1], position[2])
                .ok_or(ImpulseError::Invalid(
                    "liquid wrench stencil is not located in liquid",
                ))?;
            let key = self.grid.active_keys()[slot];
            let tile = &self.tiles[&key];
            let populations = self.grid.populations(slot, lane);
            let density = populations.iter().sum::<f64>();
            let mass = match tile.cells[lane] {
                Cell3::Fluid => density,
                Cell3::Interface => tile.mass[lane],
                Cell3::Gas => {
                    return Err(ImpulseError::Invalid("liquid wrench stencil contains gas"));
                }
            };
            if !mass.is_finite()
                || mass <= 0.0
                || !density.is_finite()
                || mass > density
                || populations.iter().any(|f| !f.is_finite() || *f <= 0.0)
            {
                return Err(ImpulseError::Cell { index });
            }
            cells.push(WrenchCell3 {
                index,
                position: [x as f64 + 0.5, y as f64 + 0.5, z as f64 + 0.5],
                mass,
            });
            replacements.push((slot, lane, populations));
        }
        let plan = WrenchPlan3::prepare(&cells, wrench, controls, gate)?;
        let directions = E3.map(|(x, y, z)| [f64::from(x), f64::from(y), f64::from(z)]);
        let mut sum = ImpulseSum::new();
        for (i, (&cell, (slot, lane, populations))) in
            cells.iter().zip(&mut replacements).enumerate()
        {
            if i % 64 == 0 {
                poll(gate)?;
            }
            let key = self.grid.active_keys()[*slot];
            sum.kick_with_mass(
                populations,
                &directions,
                &W3,
                [0.0; 3],
                self.force_at_lane(key, *lane),
                plan.impulse(cell),
                cell.mass,
                controls.impulse,
                cell.index,
            )?;
        }
        let report = sum.finish(controls.impulse)?;
        let changed = report.affected_cells != 0;
        let mut cursor = 0;
        let report = plan.finish(&cells, report, controls, gate, |index| {
            let cell = cells[cursor];
            debug_assert_eq!(cell.index, index);
            let (slot, lane, after) = &replacements[cursor];
            cursor += 1;
            let before = self.grid.populations(*slot, *lane);
            let key = self.grid.active_keys()[*slot];
            let force = self.force_at_lane(key, *lane);
            // Read the same physical velocity as the population owner, then
            // integrate with RETAINED mass. Population momentum would give
            // fictitious inertia to the empty fraction of an interface cell.
            let velocity = |f: &[f64; 19]| {
                let density = f.iter().sum::<f64>();
                core::array::from_fn::<_, 3, _>(|a| {
                    (f.iter().zip(directions).map(|(f, e)| f * e[a]).sum::<f64>() + 0.5 * force[a])
                        / density
                })
            };
            let old = velocity(&before);
            let new = velocity(after);
            core::array::from_fn(|a| cell.mass * (new[a] - old[a]))
        })?;
        let report = admit(report)?;
        poll(gate)?;
        if changed {
            for (slot, lane, populations) in replacements {
                self.grid.set_populations(slot, lane, populations);
            }
        }
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
