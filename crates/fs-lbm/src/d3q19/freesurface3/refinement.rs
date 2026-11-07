//! Uniform conservative resolution transfer of the existing sparse VOF owner.

use super::*;
use crate::d3q19::sparse::MORTON_COORD_BITS;
use crate::d3q19::{CollisionModel3, force_corrected_stress3};
use fs_exec::reduce::Compensated;

/// Explicit target size, physical-state and remap-error limits.
/// All dimensional values use the ORIGINAL coarse lattice units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DyadicRefinementControls3 {
    pub max_cells: usize,
    pub max_active_tiles: usize,
    pub maximum_speed: f64,
    pub maximum_mass_defect: f64,
    pub maximum_first_moment_defect: f64,
    pub maximum_momentum_defect: f64,
    pub maximum_angular_momentum_defect: f64,
    pub maximum_kinetic_energy_defect: f64,
    /// Maximum cellwise L1 error of the six physical BGK viscous stress modes.
    pub maximum_viscous_stress_defect: f64,
}

/// Actual located liquid inventory in ORIGINAL coarse lattice units.
/// First and angular moments are about the fixed domain's lower corner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DyadicRefinementInventory3 {
    pub mass: f64,
    pub first_mass_moment: [f64; 3],
    pub momentum: [f64; 3],
    pub angular_momentum: [f64; 3],
    pub kinetic_energy: f64,
}

/// Readback of the same physical inventory before/after resolution transfer.
#[derive(Clone, Copy, Debug, PartialEq)]
#[must_use]
pub struct DyadicRefinementReport3 {
    pub source_cells: [usize; 3],
    pub target_cells: [usize; 3],
    pub source_active_tiles: usize,
    pub target_active_tiles: usize,
    pub before: DyadicRefinementInventory3,
    pub after: DyadicRefinementInventory3,
    /// Sum of local absolute defects plus the absolute global defect.
    pub absolute_mass_defect: f64,
    pub absolute_first_moment_defect: f64,
    pub absolute_momentum_defect: f64,
    pub absolute_angular_momentum_defect: f64,
    /// Signed numerical change, evaluated without subtracting large energies.
    /// This is remap error, never physical work or heat.
    pub kinetic_energy_change: f64,
    pub absolute_kinetic_energy_defect: f64,
    pub maximum_viscous_stress_defect: f64,
    pub source_geometry_epoch: u64,
    pub target_geometry_epoch: u64,
}

fn invalid(reason: &'static str) -> FreeSurfaceError3 {
    FreeSurfaceError3::Transport(reason)
}

fn scale(value: f64, factor: f64) -> Result<f64, FreeSurfaceError3> {
    let result = value * factor;
    if !result.is_finite() || (value != 0.0 && result == 0.0) {
        return Err(invalid("dyadic physical scaling is unrepresentable"));
    }
    Ok(result)
}

fn vector(value: [f64; 3], factor: f64) -> Result<[f64; 3], FreeSurfaceError3> {
    Ok([
        scale(value[0], factor)?,
        scale(value[1], factor)?,
        scale(value[2], factor)?,
    ])
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1].mul_add(b[2], -a[2] * b[1]),
        a[2].mul_add(b[0], -a[0] * b[2]),
        a[0].mul_add(b[1], -a[1] * b[0]),
    ]
}

fn l1(values: [f64; 3]) -> f64 {
    values.into_iter().map(f64::abs).sum()
}

fn acceleration(
    field: BodyAccelerationField3,
) -> Result<BodyAccelerationField3, FreeSurfaceError3> {
    Ok(match field {
        BodyAccelerationField3::Uniform(value) => {
            BodyAccelerationField3::Uniform(vector(value, 0.5)?)
        }
        BodyAccelerationField3::Affine {
            origin,
            acceleration,
            gradient,
        } => BodyAccelerationField3::Affine {
            origin: vector(origin, 2.0)?,
            acceleration: vector(acceleration, 0.5)?,
            gradient: [
                vector(gradient[0], 0.25)?,
                vector(gradient[1], 0.25)?,
                vector(gradient[2], 0.25)?,
            ],
        },
    })
}

/// Canonical density and physical velocity, in the explicitly retained force frame.
fn state(
    populations: &[f64; Q3],
    force: [f64; 3],
    speed: f64,
) -> Result<(f64, [f64; 3]), FreeSurfaceError3> {
    let rho = populations.iter().sum::<f64>();
    let velocity: [f64; 3] = core::array::from_fn(|a| {
        let momentum = populations
            .iter()
            .zip(E3)
            .map(|(f, e)| f * f64::from([e.0, e.1, e.2][a]))
            .sum::<f64>();
        (momentum + 0.5 * force[a]) / rho
    });
    if !rho.is_finite()
        || rho <= 0.0
        || populations.iter().any(|f| !f.is_finite() || *f <= 0.0)
        || force.into_iter().chain(velocity).any(|x| !x.is_finite())
        || velocity[0].hypot(velocity[1]).hypot(velocity[2]) > speed
    {
        return Err(invalid(
            "dyadic transfer needs positive populations and admitted physical speed",
        ));
    }
    Ok((rho, velocity))
}

impl FreeSurface3 {
    /// Refine the SAME bounded physical domain by two in every spatial axis.
    ///
    /// Convective scaling halves both dx and dt. Native velocity and pressure
    /// units therefore stay fixed, BGK tau becomes 2*tau-1/2, force density
    /// halves, and surface tension doubles. The physical Guo population frame
    /// is rescaled by tau_f/(2*tau_c), then rebound to the fine held force.
    /// This preserves the physical viscous stress (the nonequilibrium second
    /// moment multiplied by 1-1/(2*tau)), rather than copying raw populations.
    ///
    /// Each original voxel-average liquid state supplies eight symmetric
    /// children, with the SAME fill and intrinsic density. Their physical
    /// masses are one eighth of the parent. This preserves the P0 inventory's
    /// first mass moment, linear/angular momentum and kinetic storage; it does
    /// not reconstruct a hidden subcell interface or higher velocity modes.
    /// The subsequent finer discretization need not follow the same trajectory.
    ///
    /// Held cellwise forces are prolonged piecewise constantly in physical
    /// space. Uniform/affine acceleration retains its original ownership and
    /// is sampled at the actual child centers, with a zero-time Guo reframe.
    /// Static solid tiles, uniform contact angle and integer piston geometry
    /// retain their physical location. Source-local capillarity/bulk relaxation,
    /// spatial contact-angle history, non-BGK collision and unlocated carry
    /// require their own transfer owners and refuse here.
    ///
    /// Constructs one sparse target with eight times the active tiles, without
    /// scanning the full inactive domain. Source state is immutable on every
    /// outcome. Native balances, positivity, Mach/capillary limits and bounded
    /// cancellation all precede returning the candidate. Geometry epoch advances
    /// because cell identities changed; no time, work, heat or conversion event
    /// advances. Cumulative conversion counts remain; the old mesh's last
    /// conversion-only receipt stays on the source, not the new generation.
    pub fn refine_dyadic(
        &self,
        controls: DyadicRefinementControls3,
        gate: &CancelGate,
    ) -> Result<(Self, DyadicRefinementReport3), FreeSurfaceError3> {
        check_cancel(gate)?;
        let budgets = [
            controls.maximum_mass_defect,
            controls.maximum_first_moment_defect,
            controls.maximum_momentum_defect,
            controls.maximum_angular_momentum_defect,
            controls.maximum_kinetic_energy_defect,
            controls.maximum_viscous_stress_defect,
        ];
        if !controls.maximum_speed.is_finite()
            || controls.maximum_speed <= 0.0
            || controls.maximum_speed > crate::MACH_LIMIT * CS2_3.sqrt()
            || budgets.into_iter().any(|v| !v.is_finite() || v < 0.0)
        {
            return Err(invalid("invalid dyadic speed or physical-error budget"));
        }
        if self.carry != 0.0
            || self.local_sigma.is_some()
            || self.local_bulk_relaxation.is_some()
            || self.planar_contact_field.is_some()
        {
            return Err(invalid(
                "dyadic refinement cannot transfer unlocated mass or source-local interface/bulk history",
            ));
        }
        let CollisionModel3::Bgk { tau } = self.collision_model() else {
            return Err(invalid(
                "dyadic refinement requires the original uniform BGK collision",
            ));
        };
        let fine_tau = 2.0_f64.mul_add(tau, -0.5);
        if !tau.is_finite() || tau <= 0.5 || !fine_tau.is_finite() || fine_tau <= 0.5 {
            return Err(invalid("dyadic BGK relaxation is unrepresentable"));
        }
        let source_cells = self.domain_cells();
        let mut target_cells = [0; 3];
        for a in 0..3 {
            target_cells[a] = source_cells[a]
                .checked_mul(2)
                .ok_or_else(|| invalid("dyadic domain size overflows"))?;
            if target_cells[a] / TILE > (1usize << MORTON_COORD_BITS) {
                return Err(invalid("dyadic domain exceeds native Morton coordinates"));
            }
        }
        let count = target_cells
            .into_iter()
            .try_fold(1usize, |n, m| n.checked_mul(m))
            .ok_or_else(|| invalid("dyadic target cell count overflows"))?;
        let target_tiles = self
            .tiles
            .len()
            .checked_mul(8)
            .ok_or_else(|| invalid("dyadic active tile count overflows"))?;
        if count == 0
            || count > controls.max_cells
            || target_tiles > controls.max_active_tiles
            || target_tiles > isize::MAX as usize / super::super::sparse::state_bytes_per_tile()
        {
            return Err(invalid(
                "dyadic target exceeds its cell or active-tile budget",
            ));
        }
        // The original binary constructor retains authored solid keys even
        // outside its physical domain. Validate those keys before allocating
        // children: doubling an otherwise valid far Morton key can overflow
        // the encoder's range and panic instead of producing a candidate.
        let source_tiles = source_cells.map(|n| n / TILE);
        let parent_morton_limit = (1u32 << MORTON_COORD_BITS) / 2;
        for &key in &self.solid {
            check_cancel(gate)?;
            let (x, y, z) = demorton3(key);
            if [x, y, z]
                .into_iter()
                .zip(source_tiles)
                .any(|(coordinate, extent)| {
                    coordinate as usize >= extent || coordinate >= parent_morton_limit
                })
            {
                return Err(invalid(
                    "dyadic source solid is outside the physical or refinable Morton domain",
                ));
            }
        }
        let sigma = scale(self.sigma, 2.0)?;
        if !crate::admit_explicit_capillary_coefficient(sigma)
            || sigma * Self::CURVATURE_LIMIT / CS2_3 >= 1.0
        {
            return Err(invalid(
                "refined timestep exceeds the native capillary-resolution bound",
            ));
        }
        let geometry_epoch = self
            .geometry_epoch
            .checked_add(1)
            .ok_or_else(|| invalid("dyadic geometry epoch overflows"))?;
        let body_acceleration = self.body_acceleration.map(acceleration).transpose()?;
        if body_acceleration.is_some() && self.local_force.is_none() {
            return Err(invalid(
                "dyadic acceleration lost its accepted force history",
            ));
        }
        let mut targets = Vec::new();
        targets
            .try_reserve_exact(target_tiles)
            .map_err(|_| invalid("dyadic sparse topology allocation failed"))?;
        for &key in self.grid.active_keys() {
            check_cancel(gate)?;
            let (x, y, z) = demorton3(key);
            for child in 0..8u32 {
                targets.push((
                    2 * x + (child & 1),
                    2 * y + ((child >> 1) & 1),
                    2 * z + (child >> 2),
                ));
            }
        }
        let mut grid = SparseGrid3::new(
            target_cells[0],
            target_cells[1],
            target_cells[2],
            fine_tau,
            vector(self.grid.collision().1, 0.5)?,
        )?;
        grid.activate_tiles(&targets)?;
        let mut solid = std::collections::BTreeSet::new();
        for &key in &self.solid {
            check_cancel(gate)?;
            let (x, y, z) = demorton3(key);
            for child in 0..8u32 {
                solid.insert(morton3(
                    2 * x + (child & 1),
                    2 * y + ((child >> 1) & 1),
                    2 * z + (child >> 2),
                ));
            }
        }
        let mut tiles = BTreeMap::new();
        for &key in grid.active_keys() {
            check_cancel(gate)?;
            tiles.insert(key, FsTile::all_gas());
        }
        let mut candidate = Self {
            grid,
            tiles,
            solid,
            sigma,
            local_sigma: None,
            local_bulk_relaxation: None,
            redistribution: self.redistribution,
            contact: self.contact,
            planar_contact: self.planar_contact,
            planar_contact_field: None,
            planar_piston: self.planar_piston.map(|p| PlanarPiston3 {
                wall: p.wall,
                depth_cells: p.depth_cells * 2,
            }),
            geometry_epoch,
            carry: 0.0,
            conversions: self.conversions,
            conversion_momentum: None,
            worst_step_violation: self.worst_step_violation,
            fill_smooth: BTreeMap::new(),
            local_force: self.local_force.as_ref().map(|_| BTreeMap::new()),
            body_acceleration,
            material_gravity: self.material_gravity.map(|v| vector(v, 0.5)).transpose()?,
        };
        let factor = fine_tau / (2.0 * tau);
        let old_stress_factor = (tau - 0.5) / tau;
        let new_stress_factor = (fine_tau - 0.5) / fine_tau;
        let mut audit = Audit::new();
        for (&key, tile) in &self.tiles {
            check_cancel(gate)?;
            let slot = self.grid.slot_of(key).expect("original active tile");
            let (tx, ty, tz) = demorton3(key);
            let target_slots: [usize; 8] = core::array::from_fn(|child| {
                candidate
                    .grid
                    .slot_of(morton3(
                        2 * tx + (child as u32 & 1),
                        2 * ty + ((child as u32 >> 1) & 1),
                        2 * tz + (child as u32 >> 2),
                    ))
                    .expect("prepared dyadic child tile")
            });
            for lane in 0..TILE_CELLS {
                if tile.cells[lane] == Cell3::Gas {
                    continue;
                }
                let populations = self.grid.populations(slot, lane);
                let force = match &self.local_force {
                    Some(field) => *field.get(&(key, lane)).ok_or_else(|| {
                        invalid("dyadic source is missing an accepted wet-cell force")
                    })?,
                    None => self.grid.collision().1,
                };
                let (rho, velocity) = state(&populations, force, controls.maximum_speed)?;
                let mass = match tile.cells[lane] {
                    Cell3::Fluid => rho,
                    Cell3::Interface => tile.mass[lane],
                    Cell3::Gas => unreachable!(),
                };
                if !mass.is_finite()
                    || mass < 0.0
                    || mass > rho
                    || (mass > 0.0 && mass * 0.125 == 0.0)
                    || !crate::admit_explicit_capillary_density(sigma, rho)
                {
                    return Err(invalid(
                        "dyadic source has unresolved fill or unrepresentable child mass/capillarity",
                    ));
                }
                let global = Self::global_of(key, lane);
                if self.planar_piston_contains([global.0, global.1, global.2]) {
                    return Err(invalid("dyadic liquid overlaps its native piston"));
                }
                let position = [
                    global.0 as f64 + 0.5,
                    global.1 as f64 + 0.5,
                    global.2 as f64 + 0.5,
                ];
                let physical = forcing::reframe(populations, velocity, force, [0.0; 3])?;
                let equilibrium = equilibrium3(rho, velocity);
                let mut scaled = core::array::from_fn(|q| {
                    factor.mul_add(physical[q] - equilibrium[q], equilibrium[q])
                });
                forcing::preserve_density(&mut scaled, rho)?;
                let stress = force_corrected_stress3(&populations, rho, velocity, force)
                    .map(|v| v * old_stress_factor);
                let (lx, ly, lz) = lane_coords(lane);
                let half_tile = TILE / 2;
                let child_slot =
                    target_slots[(lx / half_tile) + 2 * (ly / half_tile) + 4 * (lz / half_tile)];
                let child_key = candidate.grid.active_keys()[child_slot];
                let mut children = [(0.0, [0.0; 3], [0.0; 3]); 8];
                for child in 0..8 {
                    let cx = 2 * (lx % half_tile) + (child & 1);
                    let cy = 2 * (ly % half_tile) + ((child >> 1) & 1);
                    let cz = 2 * (lz % half_tile) + (child >> 2);
                    let child_lane = lane_of(cx, cy, cz);
                    let next_force = match candidate.body_acceleration {
                        Some(field) => field.at_cell(child_key, child_lane)?.map(|a| a * rho),
                        None => vector(force, 0.5)?,
                    };
                    let fine_populations =
                        forcing::reframe(scaled, velocity, [0.0; 3], next_force)?;
                    candidate
                        .grid
                        .set_populations(child_slot, child_lane, fine_populations);
                    let target = candidate
                        .tiles
                        .get_mut(&child_key)
                        .expect("prepared child tile");
                    target.cells[child_lane] = tile.cells[lane];
                    target.mass[child_lane] = mass;
                    if let Some(field) = &mut candidate.local_force {
                        field.insert((child_key, child_lane), next_force);
                    }
                    // Read actual stored state. Canonical density must retain
                    // the parent's bits; otherwise fill/pressure would change.
                    let actual = candidate.grid.populations(child_slot, child_lane);
                    let (density, u) = state(&actual, next_force, controls.maximum_speed)?;
                    if density != rho {
                        return Err(invalid("dyadic transfer changed canonical liquid density"));
                    }
                    let child_mass = match candidate.tiles[&child_key].cells[child_lane] {
                        Cell3::Fluid => density * 0.125,
                        Cell3::Interface => candidate.tiles[&child_key].mass[child_lane] * 0.125,
                        Cell3::Gas => unreachable!(),
                    };
                    let offset = [
                        (child & 1) as f64 * 0.5 - 0.25,
                        ((child >> 1) & 1) as f64 * 0.5 - 0.25,
                        (child >> 2) as f64 * 0.5 - 0.25,
                    ];
                    children[child] = (child_mass, offset, u);
                    let next_stress = force_corrected_stress3(&actual, density, u, next_force);
                    let defect = (0..6)
                        .map(|a| (next_stress[a] * new_stress_factor - stress[a]).abs())
                        .sum::<f64>();
                    if !defect.is_finite() {
                        return Err(invalid("dyadic physical stress is unrepresentable"));
                    }
                    audit.stress = audit.stress.max(defect);
                }
                audit.cell(mass, position, velocity, &children)?;
            }
        }
        candidate.assert_closure();
        let report = audit.finish(self, &candidate, controls)?;
        check_cancel(gate)?;
        Ok((candidate, report))
    }
}

struct Inventory {
    mass: Compensated,
    first: [Compensated; 3],
    momentum: [Compensated; 3],
    angular: [Compensated; 3],
    kinetic: Compensated,
}
impl Inventory {
    fn new() -> Self {
        Self {
            mass: Compensated::zero(),
            first: [Compensated::zero(); 3],
            momentum: [Compensated::zero(); 3],
            angular: [Compensated::zero(); 3],
            kinetic: Compensated::zero(),
        }
    }
    fn add(&mut self, mass: f64, position: [f64; 3], velocity: [f64; 3]) {
        let momentum = velocity.map(|u| mass * u);
        let angular = cross(position, momentum);
        self.mass = self.mass.accumulate(mass);
        for a in 0..3 {
            self.first[a] = self.first[a].accumulate(mass * position[a]);
            self.momentum[a] = self.momentum[a].accumulate(momentum[a]);
            self.angular[a] = self.angular[a].accumulate(angular[a]);
        }
        self.kinetic = self
            .kinetic
            .accumulate(0.5 * mass * velocity.iter().map(|v| v * v).sum::<f64>());
    }
    fn value(&self) -> DyadicRefinementInventory3 {
        DyadicRefinementInventory3 {
            mass: self.mass.value(),
            first_mass_moment: self.first.map(Compensated::value),
            momentum: self.momentum.map(Compensated::value),
            angular_momentum: self.angular.map(Compensated::value),
            kinetic_energy: self.kinetic.value(),
        }
    }
}

struct Audit {
    before: Inventory,
    after: Inventory,
    delta: [Compensated; 11],
    absolute: [Compensated; 5],
    stress: f64,
}
impl Audit {
    fn new() -> Self {
        Self {
            before: Inventory::new(),
            after: Inventory::new(),
            delta: [Compensated::zero(); 11],
            absolute: [Compensated::zero(); 5],
            stress: 0.0,
        }
    }
    fn cell(
        &mut self,
        mass: f64,
        position: [f64; 3],
        velocity: [f64; 3],
        children: &[(f64, [f64; 3], [f64; 3]); 8],
    ) -> Result<(), FreeSurfaceError3> {
        self.before.add(mass, position, velocity);
        let mut dm = Compensated::zero().accumulate(-mass);
        let mut relative_first = [Compensated::zero(); 3];
        let mut dp = [Compensated::zero(); 3];
        let mut spin = [Compensated::zero(); 3];
        let mut kinetic = Compensated::zero();
        for &(m, offset, u) in children {
            self.after
                .add(m, core::array::from_fn(|a| position[a] + offset[a]), u);
            dm = dm.accumulate(m);
            let du: [f64; 3] = core::array::from_fn(|a| u[a] - velocity[a]);
            let moment = cross(offset, u.map(|v| m * v));
            for a in 0..3 {
                relative_first[a] = relative_first[a].accumulate(m * offset[a]);
                dp[a] = dp[a].accumulate(m * du[a]);
                spin[a] = spin[a].accumulate(moment[a]);
                kinetic = kinetic.accumulate(m * du[a] * (velocity[a] + 0.5 * du[a]));
            }
        }
        let dm = dm.value();
        let first = core::array::from_fn(|a| dm.mul_add(position[a], relative_first[a].value()));
        let momentum = core::array::from_fn(|a| dm.mul_add(velocity[a], dp[a].value()));
        let orbital = cross(position, momentum);
        let angular = core::array::from_fn(|a| orbital[a] + spin[a].value());
        let kinetic = kinetic.value() + 0.5 * dm * velocity.iter().map(|v| v * v).sum::<f64>();
        let values = [
            dm,
            first[0],
            first[1],
            first[2],
            momentum[0],
            momentum[1],
            momentum[2],
            angular[0],
            angular[1],
            angular[2],
            kinetic,
        ];
        if values.into_iter().any(|v| !v.is_finite()) {
            return Err(invalid("dyadic moment readback is unrepresentable"));
        }
        for (sum, value) in self.delta.iter_mut().zip(values) {
            *sum = sum.accumulate(value);
        }
        for (sum, value) in self.absolute.iter_mut().zip([
            dm.abs(),
            l1(first),
            l1(momentum),
            l1(angular),
            kinetic.abs(),
        ]) {
            *sum = sum.accumulate(value);
        }
        Ok(())
    }
    fn finish(
        self,
        source: &FreeSurface3,
        target: &FreeSurface3,
        limits: DyadicRefinementControls3,
    ) -> Result<DyadicRefinementReport3, FreeSurfaceError3> {
        let before = self.before.value();
        let after = self.after.value();
        let d = self.delta.map(Compensated::value);
        let absolute = self.absolute.map(Compensated::value);
        let defects = [
            absolute[0] + d[0].abs(),
            absolute[1] + l1([d[1], d[2], d[3]]),
            absolute[2] + l1([d[4], d[5], d[6]]),
            absolute[3] + l1([d[7], d[8], d[9]]),
            absolute[4] + d[10].abs(),
            self.stress,
        ];
        for inventory in [before, after] {
            if [inventory.mass, inventory.kinetic_energy]
                .into_iter()
                .chain(inventory.first_mass_moment)
                .chain(inventory.momentum)
                .chain(inventory.angular_momentum)
                .any(|v| !v.is_finite())
            {
                return Err(invalid("dyadic inventory is unrepresentable"));
            }
        }
        if defects
            .into_iter()
            .zip([
                limits.maximum_mass_defect,
                limits.maximum_first_moment_defect,
                limits.maximum_momentum_defect,
                limits.maximum_angular_momentum_defect,
                limits.maximum_kinetic_energy_defect,
                limits.maximum_viscous_stress_defect,
            ])
            .any(|(value, limit)| !value.is_finite() || value > limit)
        {
            return Err(invalid(
                "dyadic physical inventory/stress exceeds its error budget",
            ));
        }
        Ok(DyadicRefinementReport3 {
            source_cells: source.domain_cells(),
            target_cells: target.domain_cells(),
            source_active_tiles: source.tiles.len(),
            target_active_tiles: target.tiles.len(),
            before,
            after,
            absolute_mass_defect: defects[0],
            absolute_first_moment_defect: defects[1],
            absolute_momentum_defect: defects[2],
            absolute_angular_momentum_defect: defects[3],
            kinetic_energy_change: d[10],
            absolute_kinetic_energy_defect: defects[4],
            maximum_viscous_stress_defect: defects[5],
            source_geometry_epoch: source.geometry_epoch,
            target_geometry_epoch: target.geometry_epoch,
        })
    }
}

#[cfg(test)]
mod tests;
