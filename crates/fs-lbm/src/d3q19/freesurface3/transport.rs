//! Exact native mass routes for a conservative extensive companion field.
//!
//! These records are emitted by collision/streaming and conversion themselves.
//! They do not reconstruct a second flux from cell velocities. The carry is a
//! real node because native conversion redistributes mass through that pool.

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
    max_routes: usize,
}

impl MassTransport3 {
    pub(crate) fn new(before: BTreeMap<MassLocation3, f64>, max_routes: usize) -> Self {
        Self {
            before,
            after: BTreeMap::new(),
            routes: Vec::new(),
            max_routes,
        }
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
        let MassLocation3::Cell(key, lane) = location else {
            return Err(FreeSurfaceError3::Transport("carry has no resolved stress"));
        };
        if !tau.is_finite() || tau <= 0.5 {
            return Err(FreeSurfaceError3::Transport("invalid stress relaxation"));
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
        let rate = crate::d3q19::force_corrected_shear_rate3(
            &self.grid.populations(slot, lane),
            rho,
            u,
            force,
            tau,
        );
        if !rate.is_finite() || rate < 0.0 {
            return Err(FreeSurfaceError3::Transport(
                "unrepresentable liquid strain rate",
            ));
        }
        Ok(rate)
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
    pub(crate) fn step_with_mass_transport(
        &mut self,
        gate: &CancelGate,
        tau: &BTreeMap<(u64, usize), f64>,
        forces: Option<BTreeMap<(u64, usize), [f64; 3]>>,
        max_routes: usize,
    ) -> Result<MassTransport3, FreeSurfaceError3> {
        if let Some(forces) = forces {
            self.rebind_force_field(forces, gate)?;
        }
        let mut trace = MassTransport3::new(self.mass_snapshot(), max_routes);
        self.step_inner(gate, Some(&mut trace), Some(tau))?;
        trace.after = self.mass_snapshot();
        Ok(trace)
    }
}
