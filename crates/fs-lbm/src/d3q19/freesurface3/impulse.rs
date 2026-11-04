//! Instantaneous loading of the original located liquid mass on sparse tiles.

use super::{Cell3, FreeSurface3, TILE_CELLS};
use crate::d3q19::{E3, W3};
use crate::impulse::{ImpulseControls, ImpulseError, ImpulseReport, ImpulseSum, poll};
use fs_exec::CancelGate;

impl FreeSurface3 {
    /// Apply one instantaneous liquid impulse per original geometric cell.
    ///
    /// The field covers the complete fixed domain in x-fastest order, including
    /// inactive cells. Entries are lattice momentum per whole geometric cell;
    /// callers supplying an impulse per liquid volume must include their fill
    /// weighting. Nonzero entries require an original, located wet cell. Gas,
    /// inactive/solid tiles and cells occupied by a piston must receive zero.
    ///
    /// Native liquid inertia is the original population density for `Fluid`
    /// and the explicitly tracked mass for `Interface`: `delta_u = J / mass`.
    /// An unresolved interface inventory outside `0 < mass <= density` refuses
    /// rather than clamping its fill or assigning fictitious fluid inertia.
    /// The shared population kernel retains canonical density exactly, preserves
    /// physical nonequilibrium stress to roundoff, and keeps the held Guo force.
    /// The receipt reads back actual liquid momentum and signed kinetic work.
    ///
    /// This advances no time or collision/streaming step and changes no tracked
    /// mass, fill, carry, classification, topology, force field or geometry epoch.
    /// Later native evolution owns pressure propagation and boundary reactions;
    /// this operation does not apply a paired impulse to a solid or assign heat.
    ///
    /// Only affected populations are staged, in canonical Morton/lane order.
    /// Allocation is bounded by the original cell cap and nonzero input count.
    /// Any refusal or observed cancellation leaves the entire accepted owner
    /// intact. An all-zero field allocates no replacement storage and preserves
    /// all population bits.
    pub fn apply_impulse_density(
        &mut self,
        impulses: &[[f64; 3]],
        controls: ImpulseControls,
        gate: &CancelGate,
    ) -> Result<ImpulseReport<3>, ImpulseError> {
        self.apply_impulse_density_with_receipt(impulses, controls, gate, Ok)
    }

    /// Admit a dependent receipt before publishing any original populations.
    /// Physical wrappers use this for fallible unit conversion without cloning
    /// the liquid owner or permitting a failure after accepted state changes.
    pub(crate) fn apply_impulse_density_with_receipt<R>(
        &mut self,
        impulses: &[[f64; 3]],
        controls: ImpulseControls,
        gate: &CancelGate,
        admit: impl FnOnce(ImpulseReport<3>) -> Result<R, ImpulseError>,
    ) -> Result<R, ImpulseError> {
        poll(gate)?;
        let [nx, ny, nz] = self.domain_cells();
        let count = nx
            .checked_mul(ny)
            .and_then(|n| n.checked_mul(nz))
            .ok_or(ImpulseError::Invalid("invalid free-surface domain size"))?;
        controls.admit(count)?;
        if impulses.len() != count {
            return Err(ImpulseError::Invalid(
                "impulse field must cover the original free-surface domain",
            ));
        }

        // Count the COMPLETE input, not just active tiles. A requested impulse
        // on an inactive or dry cell must never silently disappear from a
        // sparse traversal. Every admissible nonzero entry is staged once below.
        let mut requested = 0;
        for (index, impulse) in impulses.iter().enumerate() {
            if index % TILE_CELLS == 0 {
                poll(gate)?;
            }
            if !impulse.iter().all(|j| j.is_finite()) {
                return Err(ImpulseError::Cell { index });
            }
            requested += usize::from(*impulse != [0.0; 3]);
        }
        let mut sum = ImpulseSum::new();
        if requested == 0 {
            let report = admit(sum.finish(controls)?)?;
            poll(gate)?;
            return Ok(report);
        }
        let mut replacements = Vec::new();
        replacements
            .try_reserve_exact(requested)
            .map_err(|_| ImpulseError::Invalid("free-surface impulse staging allocation failed"))?;
        let directions = E3.map(|(x, y, z)| [f64::from(x), f64::from(y), f64::from(z)]);

        for (&key, tile) in &self.tiles {
            poll(gate)?;
            let slot = self.grid.slot_of(key).expect("active free-surface tile");
            for lane in 0..TILE_CELLS {
                let (x, y, z) = Self::global_of(key, lane);
                let index = (z as usize * ny + y as usize) * nx + x as usize;
                let impulse = impulses[index];
                if impulse == [0.0; 3] || self.planar_piston_contains([x, y, z]) {
                    continue;
                }
                let (mut populations, mass) = match tile.cells[lane] {
                    Cell3::Gas => continue,
                    Cell3::Fluid => {
                        let populations = self.grid.populations(slot, lane);
                        let mass = populations.iter().sum::<f64>();
                        (populations, mass)
                    }
                    Cell3::Interface => (self.grid.populations(slot, lane), tile.mass[lane]),
                };
                sum.kick_with_mass(
                    &mut populations,
                    &directions,
                    &W3,
                    [0.0; 3],
                    self.force_at_lane(key, lane),
                    impulse,
                    mass,
                    controls,
                    index,
                )?;
                replacements.push((slot, lane, populations));
            }
        }
        if replacements.len() != requested {
            return Err(ImpulseError::Invalid(
                "nonzero impulse requires an original located liquid cell",
            ));
        }
        let report = admit(sum.finish(controls)?)?;
        poll(gate)?;
        for (slot, lane, populations) in replacements {
            self.grid.set_populations(slot, lane, populations);
        }
        Ok(report)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::d3q19::freesurface3::ContactModel3;

    #[test]
    fn cancellation_after_dependent_receipt_keeps_all_sparse_populations() {
        let mut liquid = FreeSurface3::new(
            4, 4, 4, 0.9, [0.0; 3], 0.0, ContactModel3::Neutral, |_, _, _| true,
        ).unwrap();
        let before = liquid.grid().state_bits();
        let mass = liquid.ledger_mass();
        let controls = ImpulseControls {
            max_cells: 64,
            maximum_speed: 0.1,
            maximum_momentum_defect: 1e-12,
            maximum_work_defect: 1e-12,
        };
        let gate = CancelGate::new_clock_free();
        let result = liquid.apply_impulse_density_with_receipt(
            &[[0.001, 0.0, 0.0]; 64],
            controls,
            &gate,
            |report| {
                assert_eq!(report.affected_cells, 64);
                gate.request();
                Ok(report)
            },
        );
        assert_eq!(result, Err(ImpulseError::Cancelled));
        assert_eq!(liquid.grid().state_bits(), before);
        assert_eq!(liquid.ledger_mass(), mass);
        assert_eq!(liquid.grid().steps(), 0);
    }
}
