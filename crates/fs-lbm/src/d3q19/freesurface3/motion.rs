//! Tangential motion of one completely flooded planar domain wall.
//!
//! The halfway wall has zero normal velocity and unchanged geometry. This
//! deliberately does not model cover/vacate events, a piston, a moving contact
//! line, or a partially wetted wall face. Original VOF collision, transport,
//! capillarity and conversion remain the only fluid evolution operators.

use super::*;
use crate::MACH_LIMIT;
use crate::d3q19::W3;

/// One accepted exact one-cell geometry translation of a complete planar piston.
///
/// This is a conservative geometry/remap receipt, not a pressure solve. The
/// liquid's native populations and interface masses are translated unchanged;
/// therefore linear momentum is retained and the reported boundary impulse/work
/// are exactly zero. Angular momentum about the fixed grid origin changes only
/// by the geometric translation crossed with retained linear momentum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarPistonRemapReport3 {
    /// Domain side from which the piston occupies complete cell layers.
    pub wall: PlanarWall3,
    /// True when the piston advanced one cell into the domain; false on retreat.
    pub advanced: bool,
    /// Covered whole-cell layers before and after this publication.
    pub depth_before_cells: usize,
    pub depth_after_cells: usize,
    /// Moving-geometry epochs before and after this publication.
    pub geometry_epoch_before: u64,
    pub geometry_epoch_after: u64,
    /// Integer translation applied to every wet cell.
    pub liquid_translation_cells: [i32; 3],
    /// Number of fluid/interface cells translated.
    pub wet_cells: usize,
    /// Native liquid ledger around the geometry publication.
    pub mass_before: f64,
    pub mass_after: f64,
    pub relative_mass_defect: f64,
    /// Physical population momentum in lattice units before/after.
    pub momentum_before: [f64; 3],
    pub momentum_after: [f64; 3],
    /// Pure remap has no hydrodynamic impulse; a later dynamics owner must
    /// supply pressure/traction work for a finite-speed piston.
    pub fluid_impulse: [f64; 3],
    pub piston_impulse: [f64; 3],
    pub fluid_work: f64,
    pub piston_work: f64,
}

#[derive(Clone, Copy)]
struct PistonWetCell3 {
    destination: [usize; 3],
    cell: Cell3,
    mass: f64,
    populations: [f64; Q3],
    local_sigma: Option<f64>,
    local_force: Option<[f64; 3]>,
}

/// Exchange with one prescribed wall during one native lattice timestep.
///
/// The raw impulse contains the numerical density-one pressure. Subtract the
/// separately reported reference-pressure impulse before interpreting a
/// pressure-perturbation load. The tangential work is independent of this datum.
/// This is boundary mechanical work, not viscous heat or a kinetic-energy
/// conservation certificate for the isothermal VOF discretization.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarWallMotionReport3 {
    /// Selected stationary plane; its material moves tangentially.
    pub wall: PlanarWall3,
    /// Prescribed wall velocity in the existing lattice inertial frame.
    pub velocity: [f64; 3],
    /// Actual incoming-minus-outgoing population momentum delivered to liquid.
    pub fluid_impulse: [f64; 3],
    /// Exactly opposite impulse on the prescribed moving substrate.
    pub wall_impulse: [f64; 3],
    /// Density-one pressure impulse into liquid, cs² times wetted area/normal.
    pub reference_pressure_fluid_impulse: [f64; 3],
    /// Wall work delivered to the liquid, velocity dot fluid impulse.
    pub fluid_work: f64,
    /// Equal and opposite work received by the wall's external actuator.
    pub wall_work: f64,
    /// Complete flooded plane area in square lattice spacings.
    pub wetted_area: f64,
    /// Native liquid mass before this accepted timestep.
    pub mass_before: f64,
    /// Native liquid mass after this accepted timestep, including native carry.
    pub mass_after: f64,
    /// Absolute relative change in native liquid mass.
    pub relative_mass_defect: f64,
    /// Maximum physical wet-cell speed at either accepted endpoint.
    pub maximum_speed: f64,
}

fn refused(coordinates: Option<[usize; 3]>, reason: &'static str) -> FreeSurfaceError3 {
    FreeSurfaceError3::InvalidWallMotion {
        coordinates,
        reason,
    }
}

fn speed(velocity: [f64; 3]) -> f64 {
    velocity[0].hypot(velocity[1]).hypot(velocity[2])
}

pub(super) struct PlanarWallStep3 {
    wall: PlanarWall3,
    extent: i64,
    velocity: [f64; 3],
    impulse: [f64; 3],
}

impl PlanarWallStep3 {
    /// Only a source outside the selected domain plane receives its motion.
    /// At domain edges this plane owns every crossing link, retaining both
    /// tangentially opposed diagonals and therefore zero added mass.
    pub(super) fn reflect(
        &mut self,
        source: [i64; 3],
        q: usize,
        density: f64,
        outgoing: f64,
    ) -> Result<f64, FreeSurfaceError3> {
        let outside = if self.wall.upper() {
            source[self.wall.axis()] >= self.extent
        } else {
            source[self.wall.axis()] < 0
        };
        if !outside {
            return Ok(outgoing);
        }
        let direction = [E3[q].0 as f64, E3[q].1 as f64, E3[q].2 as f64];
        let projected = (0..3).map(|d| direction[d] * self.velocity[d]).sum::<f64>();
        let correction = 2.0 * W3[q] * density * projected / CS2_3;
        // Preserve the exact stationary reflection when the correction is zero.
        let incoming = if correction == 0.0 {
            outgoing
        } else {
            outgoing + correction
        };
        if !outgoing.is_finite() || outgoing <= 0.0 || !incoming.is_finite() || incoming <= 0.0 {
            return Err(refused(
                None,
                "moving-wall reflection requires positive finite populations",
            ));
        }
        for (axis, component) in self.impulse.iter_mut().enumerate() {
            *component += (incoming + outgoing) * direction[axis];
            if !component.is_finite() {
                return Err(refused(
                    None,
                    "moving-wall momentum exchange is unrepresentable",
                ));
            }
        }
        Ok(incoming)
    }
}

impl FreeSurface3 {
    fn liquid_linear_momentum(&self, gate: &CancelGate) -> Result<[f64; 3], FreeSurfaceError3> {
        let mut total = [0.0; 3];
        for (&key, tile) in &self.tiles {
            check_cancel(gate)?;
            let slot = self.grid.slot_of(key).expect("active free-surface tile");
            for (lane, &cell) in tile.cells.iter().enumerate() {
                if cell == Cell3::Gas {
                    continue;
                }
                let p = Self::global_of(key, lane);
                let sample = self
                    .sample_at(p.0 as usize, p.1 as usize, p.2 as usize)
                    .ok_or_else(|| refused(None, "wet cell disappeared from piston momentum audit"))?;
                let mass = if cell == Cell3::Fluid {
                    sample.density
                } else {
                    tile.mass[lane]
                };
                for axis in 0..3 {
                    total[axis] += mass * sample.velocity[axis];
                }
            }
        }
        if total.into_iter().all(f64::is_finite) {
            Ok(total)
        } else {
            Err(refused(None, "piston momentum audit overflowed"))
        }
    }

    /// Move a complete planar piston by exactly one Cartesian cell and translate
    /// every wet VOF cell by the same integer offset. This is the conservative
    /// geometry-epoch primitive for the plug-translation limit of MR69.
    ///
    /// Advancing occupies one additional whole-cell slab at wall; retreating
    /// vacates one. Actual native populations, tracked interface mass, local
    /// force frame and local surface coefficient travel with their liquid cell.
    /// No equilibrium reset is used for newly uncovered cells. A destination
    /// outside the domain or inside a permanent solid tile refuses atomically.
    ///
    /// The remap itself changes no population momentum and performs no physical
    /// pressure work. A finite-speed piston needs a separate dynamics/traction
    /// closure; reporting zero impulse/work here prevents geometric relocation
    /// from being mislabeled as physical force.
    pub fn remap_planar_piston_one_cell(
        &mut self,
        wall: PlanarWall3,
        advanced: bool,
        maximum_relative_mass_defect: f64,
        gate: &CancelGate,
    ) -> Result<PlanarPistonRemapReport3, FreeSurfaceError3> {
        check_cancel(gate)?;
        if !maximum_relative_mass_defect.is_finite() || maximum_relative_mass_defect < 0.0 {
            return Err(refused(None, "piston remap mass tolerance must be finite and nonnegative"));
        }
        if self.planar_contact.is_some() || self.planar_contact_field.is_some() {
            return Err(refused(None, "moving piston does not yet admit a wall contact-angle field"));
        }
        let (bound_wall, depth_before) = self
            .planar_piston
            .map_or((wall, 0usize), |p| (p.wall, p.depth_cells));
        if bound_wall != wall {
            return Err(refused(None, "one free-surface state supports one planar piston wall"));
        }
        if !advanced && depth_before == 0 {
            return Err(refused(None, "planar piston cannot retreat beyond the original domain face"));
        }
        let dims = self.grid.tile_dims();
        let extents = [dims.0 * TILE, dims.1 * TILE, dims.2 * TILE];
        let axis = wall.axis();
        let depth_after = if advanced {
            depth_before.checked_add(1).ok_or_else(|| refused(None, "piston depth overflowed"))?
        } else {
            depth_before - 1
        };
        if depth_after >= extents[axis] {
            return Err(refused(None, "piston would occupy the complete fluid domain"));
        }
        let sign = if wall.upper() { -1i32 } else { 1i32 };
        let shift = if advanced { sign } else { -sign };
        let mut translation = [0i32; 3];
        translation[axis] = shift;
        let mass_before = self.ledger_mass();
        if !mass_before.is_finite() || mass_before <= 0.0 {
            return Err(refused(None, "piston remap requires positive finite liquid mass"));
        }
        let momentum_before = self.liquid_linear_momentum(gate)?;
        let mut wet = Vec::new();
        wet.try_reserve(self.tiles.len().saturating_mul(TILE_CELLS))
            .map_err(|_| refused(None, "piston remap snapshot allocation refused"))?;
        for (&key, tile) in &self.tiles {
            check_cancel(gate)?;
            let slot = self.grid.slot_of(key).expect("active free-surface tile");
            for (lane, &cell) in tile.cells.iter().enumerate() {
                if cell == Cell3::Gas {
                    continue;
                }
                let g = Self::global_of(key, lane);
                let source = [g.0 as i64, g.1 as i64, g.2 as i64];
                let mut destination = source;
                destination[axis] += i64::from(shift);
                if (0..3).any(|d| destination[d] < 0 || destination[d] >= extents[d] as i64) {
                    return Err(refused(
                        Some(source.map(|v| v as usize)),
                        "piston translation would move liquid outside the domain",
                    ));
                }
                let destination = destination.map(|v| v as usize);
                let tile_coordinates = destination.map(|v| (v / TILE) as u32);
                if self.solid.contains(&morton3(
                    tile_coordinates[0],
                    tile_coordinates[1],
                    tile_coordinates[2],
                )) {
                    return Err(refused(
                        Some(destination),
                        "piston translation would move liquid into a permanent solid tile",
                    ));
                }
                wet.push(PistonWetCell3 {
                    destination,
                    cell,
                    mass: tile.mass[lane],
                    populations: self.grid.populations(slot, lane),
                    local_sigma: self.local_sigma.as_ref().and_then(|m| m.get(&(key, lane)).copied()),
                    local_force: self.local_force.as_ref().and_then(|m| m.get(&(key, lane)).copied()),
                });
            }
        }
        if wet.is_empty() {
            return Err(refused(None, "piston remap has no wet cells to translate"));
        }
        let geometry_epoch_after = self
            .geometry_epoch
            .checked_add(1)
            .ok_or_else(|| refused(None, "moving-geometry epoch overflowed"))?;
        let mut candidate = self.clone();
        candidate.planar_piston = (depth_after > 0).then_some(PlanarPiston3 {
            wall,
            depth_cells: depth_after,
        });
        candidate.geometry_epoch = geometry_epoch_after;

        for cell in &wet {
            check_cancel(gate)?;
            candidate.ensure_active_for_advance(
                cell.destination[0] as i64,
                cell.destination[1] as i64,
                cell.destination[2] as i64,
            )?;
            if candidate.planar_piston_contains(cell.destination.map(|v| v as i64)) {
                return Err(refused(Some(cell.destination), "translated liquid remains inside piston geometry"));
            }
            if candidate
                .grid
                .resolve_source(
                    cell.destination[0] as i64,
                    cell.destination[1] as i64,
                    cell.destination[2] as i64,
                )
                .is_none()
            {
                return Err(refused(Some(cell.destination), "piston destination could not be activated"));
            }
        }

        let old_keys: Vec<u64> = candidate.tiles.keys().copied().collect();
        for key in old_keys {
            check_cancel(gate)?;
            let tile = candidate.tiles.get_mut(&key).expect("candidate tile");
            for lane in 0..TILE_CELLS {
                if tile.cells[lane] != Cell3::Gas {
                    tile.cells[lane] = Cell3::Gas;
                    tile.mass[lane] = 0.0;
                }
            }
        }
        if let Some(field) = &mut candidate.local_sigma {
            field.clear();
        }
        if let Some(field) = &mut candidate.local_force {
            field.clear();
        }
        for cell in &wet {
            check_cancel(gate)?;
            let (slot, lane) = candidate
                .grid
                .resolve_source(
                    cell.destination[0] as i64,
                    cell.destination[1] as i64,
                    cell.destination[2] as i64,
                )
                .expect("piston destination admitted above");
            let key = candidate.grid.active_keys()[slot];
            let target = candidate.tiles.get_mut(&key).expect("active destination tile");
            if target.cells[lane] != Cell3::Gas {
                return Err(refused(Some(cell.destination), "piston remap destinations alias"));
            }
            target.cells[lane] = cell.cell;
            target.mass[lane] = cell.mass;
            candidate.grid.set_populations(slot, lane, cell.populations);
            if let (Some(field), Some(value)) = (&mut candidate.local_sigma, cell.local_sigma) {
                field.insert((key, lane), value);
            }
            if let (Some(field), Some(value)) = (&mut candidate.local_force, cell.local_force) {
                field.insert((key, lane), value);
            }
        }
        candidate.fill_smooth.clear();
        candidate.assert_closure();
        let mass_after = candidate.ledger_mass();
        let relative_mass_defect = ((mass_after - mass_before) / mass_before).abs();
        if !mass_after.is_finite()
            || !relative_mass_defect.is_finite()
            || relative_mass_defect > maximum_relative_mass_defect
        {
            return Err(refused(None, "piston remap exceeds its liquid mass-conservation tolerance"));
        }
        let momentum_after = candidate.liquid_linear_momentum(gate)?;
        let momentum_scale = momentum_before
            .into_iter()
            .chain(momentum_after)
            .map(f64::abs)
            .fold(1.0, f64::max);
        if (0..3).any(|d| {
            (momentum_after[d] - momentum_before[d]).abs()
                > 256.0 * f64::EPSILON * momentum_scale
        }) {
            return Err(refused(None, "piston geometry remap changed physical liquid momentum"));
        }
        check_cancel(gate)?;
        let report = PlanarPistonRemapReport3 {
            wall,
            advanced,
            depth_before_cells: depth_before,
            depth_after_cells: depth_after,
            geometry_epoch_before: self.geometry_epoch,
            geometry_epoch_after,
            liquid_translation_cells: translation,
            wet_cells: wet.len(),
            mass_before,
            mass_after,
            relative_mass_defect,
            momentum_before,
            momentum_after,
            fluid_impulse: [0.0; 3],
            piston_impulse: [0.0; 3],
            fluid_work: 0.0,
            piston_work: 0.0,
        };
        *self = candidate;
        Ok(report)
    }

    fn admit_planar_wall_state(
        &self,
        wall: PlanarWall3,
        velocity: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(f64, f64), FreeSurfaceError3> {
        check_cancel(gate)?;
        let axis = wall.axis();
        let dims = self.grid.tile_dims();
        let extents = [dims.0 * TILE, dims.1 * TILE, dims.2 * TILE];
        let face = if wall.upper() { extents[axis] - 1 } else { 0 };
        let tangents = [(axis + 1) % 3, (axis + 2) % 3];
        let area = extents[tangents[0]]
            .checked_mul(extents[tangents[1]])
            .ok_or_else(|| refused(None, "moving-wall plane area overflows"))?;
        for a in 0..extents[tangents[0]] {
            check_cancel(gate)?;
            for b in 0..extents[tangents[1]] {
                if b % TILE == 0 {
                    check_cancel(gate)?;
                }
                let mut p = [0usize; 3];
                p[axis] = face;
                p[tangents[0]] = a;
                p[tangents[1]] = b;
                let Some(sample) = self.sample_at(p[0], p[1], p[2]) else {
                    return Err(refused(
                        Some(p),
                        "moving plane must be completely flooded; gas, inactive or solid wall cells are unsupported",
                    ));
                };
                if sample.cell != Cell3::Fluid || sample.fill_fraction != 1.0 {
                    return Err(refused(
                        Some(p),
                        "moving plane requires full liquid cells; partial wetting and moving contact lines are unsupported",
                    ));
                }
                let relative = core::array::from_fn(|d| sample.velocity[d] - velocity[d]);
                if !speed(relative).is_finite() || speed(relative) / CS2_3.sqrt() >= MACH_LIMIT {
                    return Err(refused(
                        Some(p),
                        "liquid/wall relative velocity exceeds the lattice Mach domain",
                    ));
                }
            }
        }
        // The optional moving-wall entry admits the complete physical endpoint,
        // not just an average speed that could conceal a local high-Mach cell.
        let mut maximum_speed = 0.0_f64;
        for (&key, tile) in &self.tiles {
            check_cancel(gate)?;
            let slot = self.grid.slot_of(key).expect("active free-surface tile");
            for (lane, &cell) in tile.cells.iter().enumerate() {
                if cell == Cell3::Gas {
                    continue;
                }
                let p = Self::global_of(key, lane);
                let p = [p.0 as usize, p.1 as usize, p.2 as usize];
                let f = self.grid.populations(slot, lane);
                if f.iter().any(|v| !v.is_finite() || *v <= 0.0) {
                    return Err(refused(
                        Some(p),
                        "moving-wall endpoint requires positive finite wet populations",
                    ));
                }
                let sample = self.sample_at(p[0], p[1], p[2]).expect("wet lane");
                let actual = speed(sample.velocity);
                if !sample.density.is_finite()
                    || sample.density <= 0.0
                    || !actual.is_finite()
                    || actual / CS2_3.sqrt() >= MACH_LIMIT
                {
                    return Err(refused(
                        Some(p),
                        "liquid velocity exceeds the lattice Mach domain",
                    ));
                }
                maximum_speed = maximum_speed.max(actual);
            }
        }
        Ok((area as f64, maximum_speed))
    }

    /// One original VOF timestep with prescribed tangential planar-wall motion.
    ///
    /// The entire selected domain plane must remain flooded with full liquid
    /// cells. This supports, for example, a moving bottom below a resolved free
    /// surface. A partial wall face, gas contact, normal wall velocity, evolving
    /// coverage or moving contact line refuses rather than inferring wetted
    /// face area from volume fraction. Velocity is constant over this timestep
    /// in the SAME inertial frame as the Guo-corrected liquid velocity.
    ///
    /// The reflected populations, reaction impulse and actuator work use the
    /// same link exchange. The original collision, mass routes, conversion and
    /// capillary operator advance once. Changed wall velocities are admitted
    /// afresh on every call. Zero velocity preserves stationary population bits.
    ///
    /// All admission, streaming, conversion, endpoint, conservation and
    /// cancellation checks occur on one private candidate before publication.
    /// Geometry does not change; this is not the cover/vacate part of MR69.
    ///
    /// # Errors
    /// Refuses unsupported geometry, excessive absolute/relative Mach,
    /// nonpositive populations, a mass defect above the explicit nonnegative
    /// bound, an original native failure or observed cancellation.
    pub fn step_with_planar_wall_motion(
        &mut self,
        wall: PlanarWall3,
        velocity: [f64; 3],
        maximum_relative_mass_defect: f64,
        gate: &CancelGate,
    ) -> Result<PlanarWallMotionReport3, FreeSurfaceError3> {
        let (candidate, report) =
            self.planar_wall_motion_candidate(wall, velocity, maximum_relative_mass_defect, gate)?;
        check_cancel(gate)?;
        *self = candidate;
        Ok(report)
    }

    // The physical-units wrapper finishes its own representability/speed
    // admission before installing this SAME candidate; no second clone/step.
    pub(crate) fn planar_wall_motion_candidate(
        &self,
        wall: PlanarWall3,
        velocity: [f64; 3],
        maximum_relative_mass_defect: f64,
        gate: &CancelGate,
    ) -> Result<(Self, PlanarWallMotionReport3), FreeSurfaceError3> {
        check_cancel(gate)?;
        if velocity[wall.axis()] != 0.0 {
            return Err(refused(
                None,
                "normal wall motion requires a conservative cut-volume/remap owner",
            ));
        }
        let wall_speed = speed(velocity);
        if !wall_speed.is_finite() || wall_speed / CS2_3.sqrt() >= MACH_LIMIT {
            return Err(refused(
                None,
                "wall velocity exceeds the lattice Mach domain",
            ));
        }
        if !maximum_relative_mass_defect.is_finite() || maximum_relative_mass_defect < 0.0 {
            return Err(refused(
                None,
                "moving-wall relative mass tolerance must be finite and nonnegative",
            ));
        }
        let (wetted_area, before_speed) = self.admit_planar_wall_state(wall, velocity, gate)?;
        let mass_before = self.ledger_mass();
        if !mass_before.is_finite() || mass_before <= 0.0 {
            return Err(refused(
                None,
                "moving-wall liquid mass must be finite and positive",
            ));
        }
        let dims = self.grid.tile_dims();
        let extent = [dims.0 * TILE, dims.1 * TILE, dims.2 * TILE][wall.axis()] as i64;
        let mut exchange = PlanarWallStep3 {
            wall,
            extent,
            velocity,
            impulse: [0.0; 3],
        };
        let mut candidate = self.clone();
        candidate.step_inner_with_wall(gate, None, None, None, Some(&mut exchange))?;
        let (_, after_speed) = candidate.admit_planar_wall_state(wall, velocity, gate)?;
        let mass_after = candidate.ledger_mass();
        let relative_mass_defect = ((mass_after - mass_before) / mass_before).abs();
        if !mass_after.is_finite()
            || !relative_mass_defect.is_finite()
            || relative_mass_defect > maximum_relative_mass_defect
        {
            return Err(refused(
                None,
                "moving-wall timestep exceeds its liquid mass-conservation tolerance",
            ));
        }
        let mut reference_pressure_fluid_impulse = [0.0; 3];
        reference_pressure_fluid_impulse[wall.axis()] = if wall.upper() {
            -CS2_3 * wetted_area
        } else {
            CS2_3 * wetted_area
        };
        let fluid_work = (0..3)
            .map(|d| velocity[d] * exchange.impulse[d])
            .sum::<f64>();
        if !fluid_work.is_finite() {
            return Err(refused(
                None,
                "moving-wall actuator work is unrepresentable",
            ));
        }
        let report = PlanarWallMotionReport3 {
            wall,
            velocity,
            fluid_impulse: exchange.impulse,
            wall_impulse: exchange.impulse.map(|v| -v),
            reference_pressure_fluid_impulse,
            fluid_work,
            wall_work: -fluid_work,
            wetted_area,
            mass_before,
            mass_after,
            relative_mass_defect,
            maximum_speed: before_speed.max(after_speed),
        };
        check_cancel(gate)?;
        Ok((candidate, report))
    }
}

#[cfg(test)]
mod tests;
