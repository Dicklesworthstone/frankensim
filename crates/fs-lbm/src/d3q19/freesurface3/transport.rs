//! Exact native mass routes for a conservative extensive companion field.
//!
//! These records are emitted by collision/streaming and conversion themselves.
//! They do not reconstruct a second flux from cell velocities. Conversion
//! transfers in the default Local mode go directly to adjacent surviving
//! interfaces. The carry remains
//! a real node only for excess or debt without a local geometric recipient;
//! its missing position is never used to mix disconnected liquid components.
//! Explicit LegacyGlobal replay retains the historical shared-pool routes.

use super::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum MassLocation3 {
    Cell(u64, usize),
    Carry,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct MassRoute3 {
    pub from: MassLocation3,
    pub to: MassLocation3,
    /// Nonnegative transferred lattice mass, not a rate.
    pub mass: f64,
}

pub(crate) struct MassTransport3 {
    pub before: BTreeMap<MassLocation3, f64>,
    pub after: BTreeMap<MassLocation3, f64>,
    pub routes: Vec<MassRoute3>,
    /// Measured residuals already admitted against native donor densities.
    /// They belong to the converting source, not to unrelated route nodes.
    conversion_roundoff: BTreeMap<MassLocation3, f64>,
    max_routes: usize,
}

impl MassTransport3 {
    pub(crate) fn new(before: BTreeMap<MassLocation3, f64>, max_routes: usize) -> Self {
        Self {
            before,
            after: BTreeMap::new(),
            routes: Vec::new(),
            conversion_roundoff: BTreeMap::new(),
            max_routes,
        }
    }

    pub(crate) fn record_conversion_roundoff(
        &mut self,
        source: MassLocation3,
        residual: f64,
    ) -> Result<(), FreeSurfaceError3> {
        let accumulated = self.conversion_roundoff_at(source) + residual.abs();
        if !residual.is_finite() || !accumulated.is_finite() {
            return Err(FreeSurfaceError3::Transport(
                "nonfinite native conversion roundoff",
            ));
        }
        if residual != 0.0 {
            self.conversion_roundoff.insert(source, accumulated);
        }
        Ok(())
    }

    pub(crate) fn conversion_roundoff_at(&self, source: MassLocation3) -> f64 {
        self.conversion_roundoff
            .get(&source)
            .copied()
            .unwrap_or(0.0)
    }

    pub(super) fn exchange(
        &mut self,
        from: MassLocation3,
        to: MassLocation3,
        mass: f64,
    ) -> Result<(), FreeSurfaceError3> {
        if !mass.is_finite() {
            return Err(FreeSurfaceError3::Transport(
                "nonfinite native mass exchange",
            ));
        }
        if mass == 0.0 {
            return Ok(());
        }
        if self.routes.len() >= self.max_routes {
            return Err(FreeSurfaceError3::Transport("mass-route budget exhausted"));
        }
        self.routes.push(if mass > 0.0 {
            MassRoute3 { from, to, mass }
        } else {
            MassRoute3 {
                from: to,
                to: from,
                mass: -mass,
            }
        });
        Ok(())
    }
}

#[cfg(test)]
mod viscous_tests {
    use super::*;
    use crate::d3q19::{W3, force_corrected_shear_rate3};

    #[test]
    fn g0_located_liquid_mass_retains_interface_inventory_beyond_clamped_fill() {
        let mut liquid = FreeSurface3::new(
            8, 8, 8, 0.8, [0.0; 3], 0.0, ContactModel3::Neutral, |x, _, _| x < 4,
        ).unwrap();
        let location = FreeSurface3::transport_location([3, 3, 3]).unwrap();
        let MassLocation3::Cell(key, lane) = location else { unreachable!() };
        assert_eq!(liquid.tiles[&key].cells[lane], Cell3::Interface);
        let density = liquid.sample_at(3, 3, 3).unwrap().density;
        let retained = density.next_up();
        liquid.tiles.get_mut(&key).unwrap().mass[lane] = retained;
        let sample = liquid.sample_at(3, 3, 3).unwrap();
        assert_eq!(sample.fill_fraction, 1.0);
        assert!(sample.density * sample.fill_fraction < retained);
        assert_eq!(liquid.transport_cell_mass(location).unwrap(), retained);
        assert!(liquid.transport_cell_mass(MassLocation3::Carry).is_err());
        assert!(liquid.transport_cell_mass(
            FreeSurface3::transport_location([6, 3, 3]).unwrap(),
        ).is_err());
    }

    #[test]
    fn g1_pure_shear_power_is_invariant_under_guo_force_reframing() {
        let rho = 1.0;
        let tau = 0.8;
        let rate = 0.003;
        let velocity = [0.02, -0.01, 0.005];
        let mut physical = equilibrium3(rho, velocity);
        // Pure xy shear has Pi_xy = -rho*tau*gamma/3. The D3Q19 fourth
        // moment sum(w*ex^2*ey^2)=1/9 supplies this independent field.
        for q in 0..Q3 {
            physical[q] -= 3.0 * rho * tau * rate * W3[q] * f64::from(E3[q].0) * f64::from(E3[q].1);
        }
        let expected = force_corrected_shear_rate3(&physical, rho, velocity, [0.0; 3], tau);
        assert!((expected - rate).abs() < 1e-14);
        let first = [0.0003, -0.0004, 0.0002];
        let second = [-0.0008, 0.0001, 0.0005];
        let f = super::super::forcing::reframe(physical, velocity, [0.0; 3], first).unwrap();
        let g = super::super::forcing::reframe(f, velocity, first, second).unwrap();
        for (populations, force) in [(f, first), (g, second)] {
            let measured = force_corrected_shear_rate3(&populations, rho, velocity, force, tau);
            assert!((measured - expected).abs() < 1e-14);
        }
    }
}

impl FreeSurface3 {
    /// Original located liquid inventory, without reconstructing it from the
    /// rounded/clamped fill fraction. Full cells own population mass; interface
    /// cells own their separately tracked mass. The unlocated carry refuses.
    pub(crate) fn transport_cell_mass(
        &self,
        location: MassLocation3,
    ) -> Result<f64, FreeSurfaceError3> {
        let MassLocation3::Cell(key, lane) = location else {
            return Err(FreeSurfaceError3::Transport("carry has no located liquid mass"));
        };
        let tile = self.tiles.get(&key).ok_or(FreeSurfaceError3::Transport(
            "liquid mass requires an active cell",
        ))?;
        let mass = match tile.cells.get(lane) {
            Some(Cell3::Fluid) => {
                let slot = self.grid.slot_of(key).ok_or(FreeSurfaceError3::Transport(
                    "liquid mass requires active populations",
                ))?;
                self.grid.populations(slot, lane).iter().sum::<f64>()
            }
            Some(Cell3::Interface) => tile.mass[lane],
            _ => return Err(FreeSurfaceError3::Transport("liquid mass requires a wet cell")),
        };
        if !mass.is_finite() || mass <= 0.0 {
            return Err(FreeSurfaceError3::Transport("invalid located liquid mass"));
        }
        Ok(mass)
    }

    /// Read the tracked interface inventory and its population density without
    /// reconstructing mass from the rounded fill fraction or scanning the grid.
    pub(crate) fn transport_interface_mass(
        &self,
        location: MassLocation3,
    ) -> Result<(f64, f64), FreeSurfaceError3> {
        let MassLocation3::Cell(key, lane) = location else {
            return Err(FreeSurfaceError3::Transport("phase transfer needs a located interface"));
        };
        let slot = self.grid.slot_of(key).ok_or(FreeSurfaceError3::Transport(
            "phase transfer requires an active interface",
        ))?;
        let tile = self.tiles.get(&key).ok_or(FreeSurfaceError3::Transport(
            "phase transfer requires an active interface",
        ))?;
        let mass = tile.mass[lane];
        let density = self.grid.populations(slot, lane).iter().sum::<f64>();
        if tile.cells[lane] != Cell3::Interface
            || !density.is_finite()
            || density <= 0.0
            || !mass.is_finite()
            || mass <= 0.0
            || mass > density
        {
            return Err(FreeSurfaceError3::Transport(
                "phase transfer needs positive bounded interface fill",
            ));
        }
        Ok((mass, density))
    }

    /// Replace only the represented liquid inventory of an existing interface.
    /// The caller owns a private candidate and the physical mass/energy audit.
    /// Populations and their physical velocity are unchanged: material crosses
    /// at the local liquid velocity. Exhaustion/overfill need a phase-front
    /// remap and are refused rather than hidden in the conversion carry.
    pub(crate) fn transport_set_interface_mass(
        &mut self,
        location: MassLocation3,
        mass: f64,
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        check_cancel(gate)?;
        let (_, density) = self.transport_interface_mass(location)?;
        let MassLocation3::Cell(key, lane) = location else {
            return Err(FreeSurfaceError3::Transport("phase transfer needs a located interface"));
        };
        if !mass.is_finite()
            || mass <= 0.0
            || mass > density
        {
            return Err(FreeSurfaceError3::Transport(
                "phase transfer must retain positive bounded interface fill; refine or remap",
            ));
        }
        self.tiles.get_mut(&key).expect("admitted interface").mass[lane] = mass;
        self.fill_smooth.clear();
        Ok(())
    }

    /// Retarget one occupied interface in its actual force frame. Its native
    /// density, tracked mass and physical nonequilibrium stress remain intact;
    /// the enclosing inflow owner accounts for momentum and kinetic energy.
    pub(crate) fn transport_set_interface_velocity(
        &mut self,
        location: MassLocation3,
        velocity: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        check_cancel(gate)?;
        let (_, density) = self.transport_interface_mass(location)?;
        let speed = velocity[0].hypot(velocity[1]).hypot(velocity[2]);
        if !speed.is_finite() || speed / CS2_3.sqrt() >= crate::MACH_LIMIT {
            return Err(FreeSurfaceError3::Transport("inflow exceeds the lattice Mach domain"));
        }
        let MassLocation3::Cell(key, lane) = location else {
            return Err(FreeSurfaceError3::Transport("inflow needs a located interface"));
        };
        let slot = self.grid.slot_of(key).expect("admitted interface");
        let (_, raw) = self.grid.cell_macros(slot, lane);
        let force = self.force_at_lane(key, lane);
        let previous = core::array::from_fn(|d| raw[d] + 0.5 * force[d] / density);
        let populations = forcing::retarget_velocity(
            self.grid.populations(slot, lane), previous, velocity, force,
        )?;
        check_cancel(gate)?;
        self.grid.set_populations(slot, lane, populations);
        Ok(())
    }

    /// Permanent mechanical wall, not an inactive tile that can later be gas.
    /// Thermal contacts must not turn an absent liquid neighbor into a solid.
    pub(crate) fn transport_wall_at(&self, coordinates: [i64; 3]) -> bool {
        let (nx, ny, nz) = self.grid.tile_dims();
        let extents = [nx * TILE, ny * TILE, nz * TILE];
        if coordinates
            .iter()
            .zip(extents)
            .any(|(&p, n)| p < 0 || p >= n as i64)
        {
            return true;
        }
        let [x, y, z] = coordinates.map(|p| (p as usize / TILE) as u32);
        self.solid.contains(&morton3(x, y, z))
    }

    /// Accepted molecular strain rate in lattice inverse-time units. Carry
    /// has no resolved velocity/stress and must never be thermalized as a cell.
    pub(crate) fn transport_shear_rate(
        &self,
        location: MassLocation3,
        tau: f64,
    ) -> Result<f64, FreeSurfaceError3> {
        if !tau.is_finite() || tau <= 0.5 {
            return Err(FreeSurfaceError3::Transport("invalid stress relaxation"));
        }
        let rate = self.transport_stress_numerators(location)?[0] / tau;
        if !rate.is_finite() {
            return Err(FreeSurfaceError3::Transport("unrepresentable liquid strain rate"));
        }
        Ok(rate)
    }

    /// Deviatoric shear and sqrt(2/3)|div(u)| in lattice units. Tied models
    /// use the supplied local shear tau for both modes; independent bulk
    /// collision retains its own trace relaxation. The material shear query
    /// never includes dilation, while molecular heating retains both modes.
    pub(crate) fn transport_strain_rates(
        &self,
        location: MassLocation3,
        tau: f64,
    ) -> Result<[f64; 2], FreeSurfaceError3> {
        if !tau.is_finite() || tau <= 0.5 {
            return Err(FreeSurfaceError3::Transport("invalid stress relaxation"));
        }
        let coordinates = Self::transport_coordinates(location)
            .ok_or(FreeSurfaceError3::Transport("carry has no resolved stress"))?;
        let model = self.collision_model_at(coordinates)?;
        let components = self.transport_stress_numerators(location)?;
        let dilation = match model {
            crate::d3q19::CollisionModel3::CentralMomentBulk { bulk_rate, .. } => components[1] * bulk_rate,
            _ => components[1] / tau,
        };
        let rates = [components[0] / tau, dilation];
        if rates.iter().any(|rate| !rate.is_finite() || *rate < 0.0) {
            return Err(FreeSurfaceError3::Transport("unrepresentable liquid strain rate"));
        }
        Ok(rates)
    }

    /// Guo-corrected deviatoric and spherical stress numerators before either
    /// molecular relaxation is applied. Newly wet source queries may use these
    /// moments before their local bulk coefficient has been resolved; no
    /// constitutive rate or collision coefficient is fabricated by this read.
    pub(crate) fn transport_stress_numerators(
        &self,
        location: MassLocation3,
    ) -> Result<[f64; 2], FreeSurfaceError3> {
        let MassLocation3::Cell(key, lane) = location else {
            return Err(FreeSurfaceError3::Transport("carry has no resolved stress"));
        };
        if self.tiles.get(&key).is_none_or(|tile| tile.cells[lane] == Cell3::Gas) {
            return Err(FreeSurfaceError3::Transport("stress requires an active liquid cell"));
        }
        let slot = self.grid.slot_of(key).ok_or(FreeSurfaceError3::Transport(
            "stress requires an active liquid cell",
        ))?;
        let (rho, raw) = self.grid.cell_macros(slot, lane);
        if !rho.is_finite() || rho <= 0.0 {
            return Err(FreeSurfaceError3::Transport("invalid stress density"));
        }
        let force = self.force_at_lane(key, lane);
        let u = core::array::from_fn(|d| raw[d] + 0.5 * force[d] / rho);
        let components = crate::d3q19::rheology_stress_components3(
            &self.grid.populations(slot, lane),
            rho,
            u,
            force,
            crate::d3q19::CollisionModel3::Bgk { tau: 1.0 },
        );
        if components.iter().any(|rate| !rate.is_finite() || *rate < 0.0) {
            return Err(FreeSurfaceError3::Transport(
                "unrepresentable liquid strain rate",
            ));
        }
        Ok(components)
    }

    pub(crate) fn mass_snapshot(&self) -> BTreeMap<MassLocation3, f64> {
        let mut masses = BTreeMap::new();
        for (&key, tile) in &self.tiles {
            let slot = self.grid.slot_of(key).expect("active");
            for lane in 0..TILE_CELLS {
                let mass = match tile.cells[lane] {
                    Cell3::Gas => continue,
                    Cell3::Fluid => self.grid.populations(slot, lane).iter().sum(),
                    Cell3::Interface => tile.mass[lane],
                };
                masses.insert(MassLocation3::Cell(key, lane), mass);
            }
        }
        masses.insert(MassLocation3::Carry, self.carry);
        masses
    }

    pub(crate) fn transport_coordinates(location: MassLocation3) -> Option<[usize; 3]> {
        let MassLocation3::Cell(key, lane) = location else {
            return None;
        };
        let (x, y, z) = Self::global_of(key, lane);
        Some([x as usize, y as usize, z as usize])
    }

    pub(crate) fn transport_location(coordinates: [usize; 3]) -> Option<MassLocation3> {
        let tile = coordinates.map(|v| v / TILE);
        if tile
            .iter()
            .any(|v| *v >= 1usize << super::super::sparse::MORTON_COORD_BITS)
        {
            return None;
        }
        let key = morton3(tile[0] as u32, tile[1] as u32, tile[2] as u32);
        let lane = lane_of(
            coordinates[0] % TILE,
            coordinates[1] % TILE,
            coordinates[2] % TILE,
        );
        Some(MassLocation3::Cell(key, lane))
    }

    /// The caller owns candidate publication: conversion or companion refusal
    /// may occur after native populations have advanced on this private clone.
    #[cfg(test)]
    pub(crate) fn step_with_mass_transport(
        &mut self,
        gate: &CancelGate,
        tau: &BTreeMap<(u64, usize), f64>,
        forces: Option<BTreeMap<(u64, usize), [f64; 3]>>,
        max_routes: usize,
    ) -> Result<MassTransport3, FreeSurfaceError3> {
        self.step_with_mass_transport_and_pressure(gate, tau, forces, max_routes, None)
    }

    /// Prescribed pressure and thermal transport share the exact same native
    /// mass routes. The enclosing thermal owner publishes the complete candidate.
    pub(crate) fn step_with_mass_transport_and_pressure(
        &mut self,
        gate: &CancelGate,
        tau: &BTreeMap<(u64, usize), f64>,
        forces: Option<BTreeMap<(u64, usize), [f64; 3]>>,
        max_routes: usize,
        gas_pressure: Option<(f64, super::pressure::GasPressureField3<'_>)>,
    ) -> Result<MassTransport3, FreeSurfaceError3> {
        check_cancel(gate)?;
        let boundary = gas_pressure
            .map(|(maximum, field)| super::pressure::GasPressure3::new(maximum, field))
            .transpose()?;
        if let Some(forces) = forces {
            self.rebind_force_field(forces, gate)?;
        }
        let mut trace = MassTransport3::new(self.mass_snapshot(), max_routes);
        self.step_inner(gate, Some(&mut trace), Some(tau), boundary.as_ref())?;
        trace.after = self.mass_snapshot();
        Ok(trace)
    }

    /// The enclosing constitutive owner admits and publishes a private candidate.
    /// Reuse native collision, mass exchange and conversion without recording
    /// unused enthalpy routes for an isothermal trajectory.
    pub(crate) fn step_with_material_relaxation_and_pressure(
        &mut self,
        gate: &CancelGate,
        tau: &BTreeMap<(u64, usize), f64>,
        gas_pressure: Option<(f64, super::pressure::GasPressureField3<'_>)>,
    ) -> Result<(), FreeSurfaceError3> {
        check_cancel(gate)?;
        let boundary = gas_pressure
            .map(|(maximum, field)| super::pressure::GasPressure3::new(maximum, field))
            .transpose()?;
        self.step_inner(gate, None, Some(tau), boundary.as_ref())
    }
}

#[cfg(test)]
mod inflow_tests {
    use super::*;

    #[test]
    fn g1_native_interface_inflow_preserves_pressure_and_carried_momentum() {
        let gate = CancelGate::new_clock_free();
        fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
            let cx = fs_exec::Cx::new(&gate, arena,
                fs_exec::StreamKey { seed: 71, kernel_id: 29, tile: 0, iteration: 0 },
                fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic);
            for rho in [0.4, 0.8, 1.3, 2.5] {
                let velocity = [0.02, -0.01, 0.03];
                let incoming = [-0.04, 0.05, -0.02];
                let coordinates = [3, 3, 3];
                let mut surface = FreeSurface3::from_cells([8; 3], 0.8,
                    [0.0003, -0.0002, 0.0001], 0.0, ContactModel3::Neutral, &[],
                    &[FreeSurfaceInitialCell3 { coordinates, fill_fraction: 0.4,
                        density: rho, velocity }], &cx).unwrap();
                let location = FreeSurface3::transport_location(coordinates).unwrap();
                let before = surface.sample_at(3, 3, 3).unwrap();
                let (mass, density) = surface.transport_interface_mass(location).unwrap();
                let next_mass = mass + 0.1 * density;
                let dm = next_mass - mass;
                let next_velocity = core::array::from_fn(|d|
                    (mass * before.velocity[d] + dm * incoming[d]) / next_mass);
                surface.transport_set_interface_mass(location, next_mass, &gate).unwrap();
                surface.transport_set_interface_velocity(location, next_velocity, &gate).unwrap();
                let after = surface.sample_at(3, 3, 3).unwrap();
                assert_eq!(before.density.to_bits(), after.density.to_bits());
                assert!((after.fill_fraction - 0.5).abs() < 2e-14);
                assert_eq!(surface.ledger_mass(), next_mass);
                for d in 0..3 {
                    assert!((next_mass * after.velocity[d] - mass * before.velocity[d]
                        - dm * incoming[d]).abs() < 2e-14 * rho);
                }
                let accepted = after;
                assert!(surface.transport_set_interface_velocity(location, [1.0; 3], &gate).is_err());
                assert_eq!(surface.sample_at(3, 3, 3).unwrap(), accepted);
            }
        });
    }
}


/// Exact, compact native state retained by the planar thermal source split.
/// All transverse cells in one layer must share populations and actual mass.
#[derive(Clone)]
pub(crate) struct VolumeSlabProfile3 {
    axis: usize,
    upper: bool,
    size: [usize; 3],
    count: usize,
    layers: BTreeMap<usize, ([f64; Q3], f64)>,
}

impl FreeSurface3 {
    /// Admit rest or the exact mapped native state from a previous source step.
    /// An arbitrary velocity/stress field cannot enter this affine split.
    pub(crate) fn admit_volume_slab(
        &self,
        axis: usize,
        upper: bool,
        previous: Option<&VolumeSlabProfile3>,
        gate: &CancelGate,
    ) -> Result<ConversionRedistribution3, FreeSurfaceError3> {
        check_cancel(gate)?;
        if !self.solid.is_empty()
            || self.planar_piston.is_some()
            || self.geometry_epoch != 0
            || self.carry != 0.0
            || self.sigma != 0.0
            || self.local_sigma.is_some()
            || self.local_bulk_relaxation.is_some()
            || self.planar_contact.is_some()
            || self.planar_contact_field.is_some()
            || !matches!(self.contact, ContactModel3::Neutral)
            || !matches!(
                self.collision_model(),
                crate::d3q19::CollisionModel3::Bgk { .. }
            )
            || self.body_acceleration != Some(BodyAccelerationField3::Uniform([0.0; 3]))
            || self.material_gravity != Some([0.0; 3])
            || self.reference_force() != [0.0; 3]
        {
            return Err(FreeSurfaceError3::Transport(
                "thermal slab bridge requires stationary neutral walls, zero acceleration/capillarity and ordinary BGK",
            ));
        }
        let actual = self.volume_slab_profile(axis, upper, gate)?;
        if let Some(previous) = previous {
            if previous.axis != axis
                || previous.upper != upper
                || previous.size != actual.size
                || previous.count != actual.count
                || previous.layers != actual.layers
            {
                return Err(FreeSurfaceError3::Transport(
                    "thermal slab continuation lost its exact mapped native state",
                ));
            }
        } else {
            let mut equilibrium = equilibrium3(1.0, [0.0; 3]);
            forcing::preserve_density(&mut equilibrium, 1.0)?;
            if actual
                .layers
                .values()
                .any(|(populations, mass)| *populations != equilibrium || *mass > 1.0)
            {
                return Err(FreeSurfaceError3::Transport(
                    "initial thermal slab requires canonical pressure-density one at exact rest",
                ));
            }
        }
        Ok(self.redistribution)
    }

    /// Retain full native moments and actual inventory, never a reconstructed
    /// equilibrium approximation. Storage is one population vector per layer.
    pub(crate) fn volume_slab_profile(
        &self,
        axis: usize,
        upper: bool,
        gate: &CancelGate,
    ) -> Result<VolumeSlabProfile3, FreeSurfaceError3> {
        check_cancel(gate)?;
        let size = self.domain_cells();
        if axis >= 3 {
            return Err(FreeSurfaceError3::Transport("invalid thermal slab axis"));
        }
        let mut profile = VolumeSlabProfile3 {
            axis,
            upper,
            size,
            count: 0,
            layers: BTreeMap::new(),
        };
        for (&key, tile) in &self.tiles {
            check_cancel(gate)?;
            let slot = self.grid.slot_of(key).expect("active");
            for lane in 0..TILE_CELLS {
                if tile.cells[lane] == Cell3::Gas
                    || (tile.cells[lane] == Cell3::Interface && tile.mass[lane] == 0.0)
                {
                    continue;
                }
                let location = MassLocation3::Cell(key, lane);
                let mass = self.transport_cell_mass(location)?;
                if !mass.is_finite() || mass <= 0.0 || self.force_at_lane(key, lane) != [0.0; 3] {
                    return Err(FreeSurfaceError3::Transport(
                        "thermal slab has invalid inventory or an unmatched force",
                    ));
                }
                let p = Self::transport_coordinates(location).ok_or(
                    FreeSurfaceError3::Transport("thermal slab has unlocated inventory"),
                )?;
                let index = if upper {
                    size[axis] - 1 - p[axis]
                } else {
                    p[axis]
                };
                let state = (self.grid.populations(slot, lane), mass);
                if profile
                    .layers
                    .get(&index)
                    .is_some_and(|previous| *previous != state)
                {
                    return Err(FreeSurfaceError3::Transport(
                        "thermal slab native moments differ across a planar layer",
                    ));
                }
                profile.layers.insert(index, state);
                profile.count += 1;
            }
        }
        check_cancel(gate)?;
        Ok(profile)
    }
}
