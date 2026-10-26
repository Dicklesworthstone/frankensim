//! Held local forces in the original Guo population frame.

use super::*;
use crate::d3q19::{W3, guo_projection3};

pub(super) fn source(q: usize, velocity: [f64; 3], force: [f64; 3]) -> f64 {
    let e = [f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2)];
    let eu = e.iter().zip(velocity).map(|(e, u)| e * u).sum();
    W3[q] * guo_projection3(e, velocity, eu, force)
}

pub(super) fn reframe(
    mut populations: [f64; Q3],
    physical_velocity: [f64; 3],
    previous: [f64; 3],
    next: [f64; 3],
) -> Result<[f64; Q3], FreeSurfaceError3> {
    if previous == next {
        return Ok(populations);
    }
    let rho = populations.iter().sum::<f64>();
    let delta = core::array::from_fn(|d| previous[d] - next[d]);
    for (q, value) in populations.iter_mut().enumerate() {
        *value += 0.5 * source(q, physical_velocity, delta);
    }
    // The analytic source has zero mass. Keep its floating-point reduction
    // defect in the rest population, never in VOF mass or a heat source.
    populations[0] += rho - populations.iter().sum::<f64>();
    if populations.iter().any(|v| !v.is_finite()) {
        return Err(FreeSurfaceError3::Transport(
            "unrepresentable force-frame change",
        ));
    }
    Ok(populations)
}

impl FreeSurface3 {
    pub(crate) fn reference_force(&self) -> [f64; 3] {
        self.grid.collision().1
    }

    pub(super) fn force_at_lane(&self, key: u64, lane: usize) -> [f64; 3] {
        self.local_force
            .as_ref()
            .map_or_else(|| self.grid.collision().1, |forces| forces[&(key, lane)])
    }

    /// Actual force density in lattice units at an accepted wet cell. Its
    /// half-source correction is already included by `sample_at`; adding it
    /// again would count the same force twice. Gas/wall/outside returns None.
    #[must_use]
    pub fn force_density_at(&self, x: usize, y: usize, z: usize) -> Option<[f64; 3]> {
        let (slot, lane) = self.grid.resolve_source(
            i64::try_from(x).ok()?,
            i64::try_from(y).ok()?,
            i64::try_from(z).ok()?,
        )?;
        let key = self.grid.active_keys()[slot];
        (self.tiles.get(&key)?.cells[lane] != Cell3::Gas).then(|| self.force_at_lane(key, lane))
    }

    /// Change only the source frame, preserving physical velocity and the
    /// force-corrected second moment. The caller owns the enclosing candidate.
    /// No elapsed time, impulse, heat, mass route or mechanical work is added.
    pub(crate) fn rebind_force_field(
        &mut self,
        forces: BTreeMap<(u64, usize), [f64; 3]>,
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        check_cancel(gate)?;
        let wet = self
            .tiles
            .values()
            .map(|tile| {
                tile.cells
                    .iter()
                    .filter(|&&cell| cell != Cell3::Gas)
                    .count()
            })
            .sum::<usize>();
        if forces.len() != wet {
            return Err(FreeSurfaceError3::Transport(
                "local force must cover each wet cell once",
            ));
        }
        let uniform = self.grid.collision().1;
        let mut replacements = Vec::new();
        for (&key, tile) in &self.tiles {
            check_cancel(gate)?;
            let slot = self.grid.slot_of(key).expect("active");
            for lane in 0..TILE_CELLS {
                if tile.cells[lane] == Cell3::Gas {
                    continue;
                }
                let next = *forces
                    .get(&(key, lane))
                    .ok_or(FreeSurfaceError3::Transport(
                        "missing local force at a wet cell",
                    ))?;
                if next.iter().any(|f| !f.is_finite()) {
                    return Err(FreeSurfaceError3::Transport("nonfinite local force"));
                }
                let previous = self.force_at_lane(key, lane);
                // Preserve the original population bits when the force does
                // not change, including the constant-density gravity limit.
                if previous == next {
                    continue;
                }
                let populations = self.grid.populations(slot, lane);
                let (rho, raw) = self.grid.cell_macros(slot, lane);
                if !rho.is_finite() || rho <= 0.0 {
                    return Err(FreeSurfaceError3::Transport(
                        "invalid density in force-frame change",
                    ));
                }
                let velocity = core::array::from_fn(|d| raw[d] + 0.5 * previous[d] / rho);
                let populations = reframe(populations, velocity, previous, next)?;
                replacements.push((slot, lane, populations));
            }
        }
        check_cancel(gate)?;
        for (slot, lane, populations) in replacements {
            self.grid.set_populations(slot, lane, populations);
        }
        self.local_force = if forces.values().all(|&force| force == uniform) {
            None
        } else {
            Some(forces)
        };
        Ok(())
    }
}
