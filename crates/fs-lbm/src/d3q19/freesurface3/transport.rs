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

impl FreeSurface3 {
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
