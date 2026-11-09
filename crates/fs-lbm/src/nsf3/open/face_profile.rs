//! Dense prescribed profiles on original complete reservoir faces.

use super::super::body_force::checkpoint;
use super::{NsfError, NsfGrid3, NsfOpenState, NsfReservoirTileState, OpenGhost};
use crate::lattice3::Lattice3;
use std::sync::Arc;

#[derive(Clone, Copy)]
struct FaceTile {
    index: usize,
    state: NsfOpenState,
    ghost: OpenGhost,
}

/// Immutable dense face payload. Trial owners share this admitted snapshot;
/// a later update stages a new one and never mutates an existing snapshot.
pub(in crate::nsf3) struct FaceProfile {
    tiles: Vec<FaceTile>,
}

fn axes(normal: usize) -> [usize; 2] {
    match normal {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    }
}

fn face_cell(shape: [usize; 3], normal: usize, high: bool, position: usize) -> [usize; 3] {
    let [first, second] = axes(normal);
    let mut cell = [0; 3];
    cell[normal] = if high { shape[normal] - 1 } else { 0 };
    cell[first] = position % shape[first];
    cell[second] = position / shape[first];
    cell
}

impl FaceProfile {
    pub(super) fn ghost(&self, shape: [usize; 3], normal: usize, cell: [usize; 3]) -> OpenGhost {
        let [first, second] = axes(normal);
        self.tiles[cell[first] + shape[first] * cell[second]].ghost
    }
}

impl<L: Lattice3<27>> NsfGrid3<L> {
    // A diagonal incoming population originates on a transversely displaced
    // exterior tile. Uniform controls keep their original face-wide cache.
    pub(in crate::nsf3) fn profiled_reservoir_link_ghost(
        &self,
        cell: [usize; 3],
        axis: usize,
        side: usize,
        velocity: [i32; 3],
    ) -> Option<OpenGhost> {
        let profile = self.reservoir_face_profiles[axis][side].as_ref()?;
        let mut donor = cell;
        for a in 0..3 {
            if a != axis {
                donor[a] = (cell[a] as i64 - i64::from(velocity[a]))
                    .rem_euclid(self.spec.dims[a] as i64) as usize;
            }
        }
        Some(profile.ghost(self.spec.dims, axis, donor))
    }

    /// Accepted states of a complete reservoir face in original z/y/x order
    /// (lowest remaining Cartesian axis varies fastest), without allocation.
    /// `None` means the axis is invalid or this complete face is not open.
    /// Uniform states are converted from the original admitted lattice chart.
    pub fn open_face_tile_states(
        &self,
        axis: usize,
        high: bool,
    ) -> Option<impl ExactSizeIterator<Item = NsfReservoirTileState> + '_> {
        let ghost = *self
            .reservoir_faces
            .get(axis)?
            .get(usize::from(high))?
            .as_ref()?;
        let profile = self.reservoir_face_profiles[axis][usize::from(high)].as_ref();
        Some(
            (0..self.n / self.spec.dims[axis]).map(move |position| NsfReservoirTileState {
                cell: face_cell(self.spec.dims, axis, high, position),
                state: profile.map_or_else(
                    || self.physical_open_state(ghost),
                    |profile| profile.tiles[position].state,
                ),
            }),
        )
    }

    fn physical_open_state(&self, ghost: OpenGhost) -> NsfOpenState {
        NsfOpenState {
            density_kg_m3: ghost.rho * self.law.rho_ref_kg_m3,
            velocity_m_s: ghost.u.map(|v| v * self.law.dx_m / self.law.dt_s),
            temperature_k: self.law.kelvin_of_theta(ghost.theta),
        }
    }

    /// Atomically replace every ghost on an already installed complete face.
    ///
    /// Input order is arbitrary; each adjacent face cell must occur exactly
    /// once. `maximum_tiles` admits the dense profile before allocation. Every
    /// state passes the original material, temperature and Mach gates. Profiles
    /// use constant-time dense lookup in the original population, pressure,
    /// full-spacing energy and SGS-gradient consumers. A uniform profile takes
    /// the original uniform path exactly, including its historical state hash.
    /// Diagonal incoming links wrap on the original transverse periodic axes;
    /// their ghost donor is distinct from the cell-normal finite-volume trace.
    ///
    /// Current physical density, velocity and total energy are preserved while
    /// the original Guo force chart and energy-flux history are rebound. All
    /// candidate transport/material checks and final cancellation precede
    /// publication. Refusal retains the prior profile and every accepted field.
    /// Cancellation is polled during tile admission/population rebase and at
    /// native force, closure, stress and flux stage boundaries; each stage is
    /// bounded by the already admitted grid. No source callback is evaluated.
    ///
    /// This is a held resolved primitive profile, not a spectrum generator,
    /// temporal interpolation rule or exact imposed streamed mass flux. The
    /// native kinetic and total-energy boundary fluxes remain authoritative.
    pub fn set_open_face_profile(
        &mut self,
        axis: usize,
        high: bool,
        states: &[NsfReservoirTileState],
        maximum_tiles: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), NsfError> {
        checkpoint(&mut cancelled)?;
        let side = usize::from(high);
        let ghost = self
            .reservoir_faces
            .get(axis)
            .and_then(|faces| faces[side])
            .ok_or(NsfError::Input(
                "a face profile requires an existing complete reservoir",
            ))?;
        if self.dynamic_sgs.is_some() && !self.dynamic_reservoir_face_update_allowed(axis, true) {
            return Err(NsfError::Input(
                "dynamic SGS requires its admitted reservoir topology",
            ));
        }
        let required = self.n / self.spec.dims[axis];
        if required > maximum_tiles {
            return Err(NsfError::Cap {
                requested: required,
                maximum: maximum_tiles,
            });
        }
        if states.len() != required {
            return Err(NsfError::Input(
                "a complete reservoir profile must cover every face tile once",
            ));
        }
        let empty = FaceTile {
            index: usize::MAX,
            state: self.physical_open_state(ghost),
            ghost,
        };
        let mut tiles = Vec::new();
        tiles
            .try_reserve_exact(required)
            .map_err(|_| NsfError::Input("complete reservoir profile allocation failed"))?;
        if tiles.capacity() > maximum_tiles {
            return Err(NsfError::Cap {
                requested: tiles.capacity(),
                maximum: maximum_tiles,
            });
        }
        tiles.resize(required, empty);
        let [first, second] = axes(axis);
        let plane = if high { self.spec.dims[axis] - 1 } else { 0 };
        for (count, tile) in states.iter().enumerate() {
            if count % 64 == 0 {
                checkpoint(&mut cancelled)?;
            }
            if tile.cell.iter().zip(self.spec.dims).any(|(&x, n)| x >= n)
                || tile.cell[axis] != plane
            {
                return Err(NsfError::Input(
                    "reservoir profile cell lies outside its complete face",
                ));
            }
            let position = tile.cell[first] + self.spec.dims[first] * tile.cell[second];
            if tiles[position].index != usize::MAX {
                return Err(NsfError::Input("duplicate complete reservoir profile tile"));
            }
            tiles[position] = FaceTile {
                index: self.index(tile.cell),
                state: tile.state,
                ghost: self.admit_open_state(tile.state)?,
            };
        }
        checkpoint(&mut cancelled)?;
        let first = tiles[0].state; // Original face dimensions are positive.
        if tiles.iter().all(|tile| tile.state == first) {
            return self.rebind_reservoir_face(
                axis,
                side,
                Some(tiles[0].ghost),
                None,
                &mut cancelled,
            );
        }
        if self.reservoir_face_profiles[axis][side]
            .as_ref()
            .is_some_and(|previous| {
                previous
                    .tiles
                    .iter()
                    .zip(&tiles)
                    .all(|(old, new)| old.state == new.state)
            })
        {
            return Ok(());
        }
        self.rebind_reservoir_face(
            axis,
            side,
            Some(ghost),
            Some(Arc::new(FaceProfile { tiles })),
            &mut cancelled,
        )
    }

    pub(in crate::nsf3) fn feed_face_profile_hash(&self, feed: &mut impl FnMut(f64)) {
        for (axis, faces) in self.reservoir_face_profiles.iter().enumerate() {
            for (side, profile) in faces.iter().enumerate() {
                if let Some(profile) = profile {
                    feed(12.0);
                    feed(axis as f64);
                    feed(side as f64);
                    feed(profile.tiles.len() as f64);
                    for tile in &profile.tiles {
                        feed(tile.index as f64);
                        feed(tile.state.density_kg_m3);
                        for v in tile.state.velocity_m_s {
                            feed(v);
                        }
                        feed(tile.state.temperature_k);
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests;
