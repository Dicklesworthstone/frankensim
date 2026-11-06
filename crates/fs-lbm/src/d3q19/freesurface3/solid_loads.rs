//! Actual pressure/viscous reaction on closed, fully submerged stationary solids.
//! The original stream supplies every reflected population; no second collision
//! or reconstructed stress field supplies the load.

use super::*;
use crate::MACH_LIMIT;
use crate::lattice3::wall_loads::StationaryWallMomentum3;

/// One stationary immersed object, possibly a disconnected rigid assembly.
/// Tile coordinates refer to the original permanently solid tile geometry.
#[derive(Clone, Copy, Debug)]
pub struct SubmergedSolid3<'a> {
    pub solid_tiles: &'a [(u32, u32, u32)],
    /// Moment origin in lattice lengths; cell centers are integer + 1/2.
    pub moment_origin: [f64; 3],
}

/// Explicit work/storage and conservation bounds for one load-producing step.
#[derive(Clone, Copy, Debug)]
pub struct SubmergedSolidLoadControls3 {
    /// Sum of selected solid tiles, including all objects.
    pub maximum_solid_tiles: usize,
    /// Total reflected links on their complete exposed surfaces.
    pub maximum_boundary_links: usize,
    pub maximum_relative_mass_defect: f64,
}

/// Actual stationary-body reactions during the same native free-surface step.
#[derive(Clone, Debug, PartialEq)]
pub struct SubmergedSolidLoads3 {
    /// Input-object order. The impulses are ON each solid; their negatives
    /// are the reciprocal liquid impulses. Work is zero at a stationary wall.
    pub bodies: Vec<StationaryWallMomentum3>,
    /// The fixed geometry epoch used by both populations and reaction arms.
    pub geometry_epoch: u64,
    pub mass_before: f64,
    pub mass_after: f64,
    pub relative_mass_defect: f64,
    /// Maximum physical lattice speed at either admitted wet endpoint.
    pub maximum_speed: f64,
}

fn invalid(reason: &'static str) -> FreeSurfaceError3 {
    FreeSurfaceError3::Transport(reason)
}

#[derive(Clone, Default)]
struct LoadSum {
    sum: [f64; 6],
    correction: [f64; 6],
    links: usize,
}

pub(super) struct SolidLoadStep3<'a> {
    owners: BTreeMap<u64, usize>,
    bodies: &'a [SubmergedSolid3<'a>],
    totals: Vec<LoadSum>,
    dims: [usize; 3],
    maximum_links: usize,
}

impl<'a> SolidLoadStep3<'a> {
    fn new(
        grid: &FreeSurface3,
        bodies: &'a [SubmergedSolid3<'a>],
        controls: SubmergedSolidLoadControls3,
        gate: &CancelGate,
    ) -> Result<Self, FreeSurfaceError3> {
        if bodies.is_empty()
            || bodies.len() > controls.maximum_solid_tiles
            || controls.maximum_boundary_links == 0
            || !controls.maximum_relative_mass_defect.is_finite()
            || controls.maximum_relative_mass_defect < 0.0
        {
            return Err(invalid(
                "submerged-solid load bounds or selection are invalid",
            ));
        }
        let mut count = 0usize;
        for body in bodies {
            check_cancel(gate)?;
            count = count
                .checked_add(body.solid_tiles.len())
                .ok_or_else(|| invalid("submerged-solid tile count overflow"))?;
            if count > controls.maximum_solid_tiles
                || body.solid_tiles.is_empty()
                || body.moment_origin.iter().any(|x| !x.is_finite())
            {
                return Err(invalid(
                    "submerged-solid tile budget, empty object or moment origin",
                ));
            }
        }
        let dims = grid.domain_cells();
        let mut owners = BTreeMap::new();
        for (body, selection) in bodies.iter().enumerate() {
            for &(x, y, z) in selection.solid_tiles {
                check_cancel(gate)?;
                let tile = [x as usize, y as usize, z as usize];
                let key = morton3(x, y, z);
                if (0..3).any(|a| tile[a] >= dims[a] / TILE)
                    || !grid.solid.contains(&key)
                    || owners.insert(key, body).is_some()
                {
                    return Err(invalid(
                        "submerged-solid tiles must be distinct original in-domain solids",
                    ));
                }
            }
        }
        Ok(Self {
            owners,
            bodies,
            totals: vec![LoadSum::default(); bodies.len()],
            dims,
            maximum_links: controls.maximum_boundary_links,
        })
    }

    /// Check every exposed D3Q19 link, including diagonal edge links. Adjacent
    /// tiles of the SAME object are internal. A wall, another object, gas or an
    /// interface at its boundary cannot supply a complete wetted pressure load.
    fn admit_surface(
        &self,
        grid: &FreeSurface3,
        gate: &CancelGate,
    ) -> Result<Vec<usize>, FreeSurfaceError3> {
        let mut counts = vec![0usize; self.bodies.len()];
        let mut total = 0usize;
        for (&key, &body) in &self.owners {
            check_cancel(gate)?;
            if !grid.solid.contains(&key) {
                return Err(invalid("submerged-solid geometry changed"));
            }
            for lane in 0..TILE_CELLS {
                let (x, y, z) = FreeSurface3::global_of(key, lane);
                for q in 1..Q3 {
                    let neighbor = [
                        x + i64::from(E3[q].0),
                        y + i64::from(E3[q].1),
                        z + i64::from(E3[q].2),
                    ];
                    if (0..3).any(|a| neighbor[a] < 0 || neighbor[a] >= self.dims[a] as i64) {
                        return Err(invalid(
                            "submerged solid must be separated from every domain wall",
                        ));
                    }
                    let neighbor_key = morton3(
                        (neighbor[0] as usize / TILE) as u32,
                        (neighbor[1] as usize / TILE) as u32,
                        (neighbor[2] as usize / TILE) as u32,
                    );
                    if self.owners.get(&neighbor_key) == Some(&body) {
                        continue;
                    }
                    if !matches!(
                        grid.neighbor_state(key, lane, q),
                        Neighbor::Active(_, _, Cell3::Fluid)
                    ) {
                        return Err(invalid(
                            "submerged-solid boundary must remain completely flooded by full liquid",
                        ));
                    }
                    if total == self.maximum_links {
                        return Err(invalid("submerged-solid boundary link budget exceeded"));
                    }
                    total += 1;
                    counts[body] += 1;
                }
            }
        }
        Ok(counts)
    }

    pub(super) fn record(
        &mut self,
        fluid: [i64; 3],
        source: [i64; 3],
        q: usize,
        incoming: f64,
        outgoing: f64,
    ) -> Result<(), FreeSurfaceError3> {
        if (0..3).any(|a| source[a] < 0 || source[a] >= self.dims[a] as i64) {
            return Ok(());
        }
        let key = morton3(
            (source[0] as usize / TILE) as u32,
            (source[1] as usize / TILE) as u32,
            (source[2] as usize / TILE) as u32,
        );
        let Some(&body) = self.owners.get(&key) else {
            return Ok(());
        };
        if !incoming.is_finite() || incoming <= 0.0 || !outgoing.is_finite() || outgoing <= 0.0 {
            return Err(invalid(
                "submerged-solid exchange requires positive finite reflected populations",
            ));
        }
        let e = [f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2)];
        let impulse = e.map(|c| -(incoming + outgoing) * c);
        let r: [f64; 3] = core::array::from_fn(|a| {
            fluid[a] as f64 + 0.5 - 0.5 * e[a] - self.bodies[body].moment_origin[a]
        });
        let torque = [
            r[1] * impulse[2] - r[2] * impulse[1],
            r[2] * impulse[0] - r[0] * impulse[2],
            r[0] * impulse[1] - r[1] * impulse[0],
        ];
        let total = &mut self.totals[body];
        // Canonical compensated sums retain small buoyancy/drag after large
        // equilibrium pressure contributions cancel around the closed object.
        for (a, v) in impulse.into_iter().chain(torque).enumerate() {
            let next = total.sum[a] + v;
            total.correction[a] += if total.sum[a].abs() >= v.abs() {
                (total.sum[a] - next) + v
            } else {
                (v - next) + total.sum[a]
            };
            total.sum[a] = next;
            if !next.is_finite() || !total.correction[a].is_finite() {
                return Err(invalid("submerged-solid reaction is unrepresentable"));
            }
        }
        total.links += 1;
        Ok(())
    }

    fn finish(self, expected: &[usize]) -> Result<Vec<StationaryWallMomentum3>, FreeSurfaceError3> {
        self.totals
            .into_iter()
            .zip(expected)
            .map(|(total, &expected)| {
                let values: [f64; 6] = core::array::from_fn(|a| total.sum[a] + total.correction[a]);
                if total.links != expected || values.iter().any(|v| !v.is_finite()) {
                    return Err(invalid(
                        "submerged-solid reflected-link coverage or total is invalid",
                    ));
                }
                Ok(StationaryWallMomentum3 {
                    wall_impulse: [values[0], values[1], values[2]],
                    wall_angular_impulse: [values[3], values[4], values[5]],
                    measured_links: total.links,
                })
            })
            .collect()
    }
}

impl FreeSurface3 {
    fn solid_load_endpoint_speed(&self, gate: &CancelGate) -> Result<f64, FreeSurfaceError3> {
        let mut maximum = 0.0_f64;
        for (&key, tile) in &self.tiles {
            check_cancel(gate)?;
            let slot = self.grid.slot_of(key).expect("active wet tile");
            for (lane, &cell) in tile.cells.iter().enumerate() {
                if cell == Cell3::Gas {
                    continue;
                }
                if self
                    .grid
                    .populations(slot, lane)
                    .iter()
                    .any(|p| !p.is_finite() || *p <= 0.0)
                {
                    return Err(invalid(
                        "submerged-solid step needs positive finite wet populations",
                    ));
                }
                let (x, y, z) = Self::global_of(key, lane);
                let sample = self
                    .sample_at(x as usize, y as usize, z as usize)
                    .ok_or_else(|| invalid("submerged-solid step lost a wet endpoint"))?;
                let v = sample.velocity[0]
                    .hypot(sample.velocity[1])
                    .hypot(sample.velocity[2]);
                if !sample.density.is_finite()
                    || sample.density <= 0.0
                    || !v.is_finite()
                    || v / CS2_3.sqrt() >= MACH_LIMIT
                {
                    return Err(invalid(
                        "submerged-solid step exceeds the physical low-Mach domain",
                    ));
                }
                maximum = maximum.max(v);
            }
        }
        Ok(maximum)
    }

    /// Advance once and return actual buoyancy, pressure and viscous impulses
    /// on every selected immersed object, about its supplied moment origin.
    ///
    /// Uses the original collision, gas reconstruction, gravity, capillarity
    /// and VOF conversion. Each solid reflection is observed exactly once;
    /// no auxiliary collision, guessed stress or added body force is used.
    /// All exposed links of every object must remain surrounded by full liquid.
    /// Uniform pressure has zero resultant and torque on these closed surfaces,
    /// so total loads require no absolute pressure datum. Raw population
    /// impulses retain the native numerical pressure for momentum accounting.
    ///
    /// Geometry and objects remain fixed. Their work is zero. These loads can
    /// drive a coupled mechanics owner, but this method does not supply body
    /// inertia, position updates, partial wetting or moving-geometry remap.
    /// Reaction sampling uses O(selected tiles + objects) auxiliary storage,
    /// with no population buffer beyond the single native candidate.
    ///
    /// # Errors
    /// Refuses incomplete wet surfaces, duplicate/invalid tiles, the explicit
    /// tile/link/mass bounds, nonfinite loads, native failure or cancellation.
    /// All refusals retain the original populations, force, mass and geometry.
    pub fn step_with_submerged_solid_loads(
        &mut self,
        bodies: &[SubmergedSolid3<'_>],
        controls: SubmergedSolidLoadControls3,
        gate: &CancelGate,
    ) -> Result<SubmergedSolidLoads3, FreeSurfaceError3> {
        let (candidate, report) = self.submerged_solid_load_candidate(bodies, controls, gate)?;
        check_cancel(gate)?;
        *self = candidate;
        Ok(report)
    }

    /// Unpublished version of [`Self::step_with_submerged_solid_loads`] for
    /// coupled-state acceptance. Dropping the candidate changes no fluid state.
    pub fn submerged_solid_load_candidate(
        &self,
        bodies: &[SubmergedSolid3<'_>],
        controls: SubmergedSolidLoadControls3,
        gate: &CancelGate,
    ) -> Result<(Self, SubmergedSolidLoads3), FreeSurfaceError3> {
        check_cancel(gate)?;
        let mut exchange = SolidLoadStep3::new(self, bodies, controls, gate)?;
        let expected = exchange.admit_surface(self, gate)?;
        let before_speed = self.solid_load_endpoint_speed(gate)?;
        let mass_before = self.ledger_mass();
        if !mass_before.is_finite() || mass_before <= 0.0 {
            return Err(invalid(
                "submerged-solid step requires positive finite liquid mass",
            ));
        }
        let mut candidate = self.clone();
        candidate.step_inner_with_boundaries(gate, None, None, None, None, Some(&mut exchange))?;
        exchange.admit_surface(&candidate, gate)?;
        let maximum_speed = before_speed.max(candidate.solid_load_endpoint_speed(gate)?);
        let mass_after = candidate.ledger_mass();
        let relative_mass_defect = ((mass_after - mass_before) / mass_before).abs();
        if !mass_after.is_finite()
            || mass_after <= 0.0
            || !relative_mass_defect.is_finite()
            || relative_mass_defect > controls.maximum_relative_mass_defect
        {
            return Err(invalid(
                "submerged-solid step exceeds its liquid mass-defect bound",
            ));
        }
        let report = SubmergedSolidLoads3 {
            bodies: exchange.finish(&expected)?,
            geometry_epoch: self.geometry_epoch,
            mass_before,
            mass_after,
            relative_mass_defect,
            maximum_speed,
        };
        check_cancel(gate)?;
        Ok((candidate, report))
    }
}

#[cfg(test)]
mod tests;
