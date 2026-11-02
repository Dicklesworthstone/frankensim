//! Spatially resolved prescribed states on one existing reservoir owner.
//! The original kinetic, thermal-force and conservative energy paths all use
//! the same admitted ghost; a profile never creates seams between its tiles.

use super::{BoundReservoirAperture, NsfError, NsfGrid3, NsfOpenState, OpenGhost};
use crate::lattice3::Lattice3;

/// Prescribed physical state at the exterior ghost centre of one aperture tile.
///
/// The cell identifies the adjacent fluid cell, not a new boundary owner. The
/// ghost centre lies one grid spacing beyond that cell along the face normal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfReservoirTileState {
    /// Original adjacent fluid-cell coordinates.
    pub cell: [usize; 3],
    /// Density, Cartesian velocity and temperature in physical units.
    pub state: NsfOpenState,
}

#[derive(Clone)]
struct ReservoirTile {
    index: usize,
    state: NsfOpenState,
    ghost: OpenGhost,
}

#[derive(Clone)]
pub(super) struct ReservoirProfile {
    // Canonical original-cell order, exactly matching the owner's indices.
    tiles: Vec<ReservoirTile>,
}

impl BoundReservoirAperture {
    pub(super) fn ghost_at(&self, index: usize) -> OpenGhost {
        self.profile.as_ref().map_or(self.ghost, |profile| {
            let tile = profile
                .tiles
                .binary_search_by_key(&index, |tile| tile.index)
                .expect("admitted reservoir tile belongs to its owner");
            profile.tiles[tile].ghost
        })
    }

    pub(in crate::nsf3) fn feed_profile_hash(&self, feed: &mut impl FnMut(f64)) {
        let Some(profile) = &self.profile else {
            return;
        };
        // Uniform controls contribute no new bits. Bind the spatial profile
        // and its geometry even for the historical single-owner hash path.
        feed(6.0);
        feed(self.source.axis as f64);
        feed(f64::from(u8::from(self.source.high)));
        feed(profile.tiles.len() as f64);
        for tile in &profile.tiles {
            feed(tile.index as f64);
            feed(tile.state.density_kg_m3);
            for velocity in tile.state.velocity_m_s {
                feed(velocity);
            }
            feed(tile.state.temperature_k);
        }
    }
}

impl<L: Lattice3<27>> NsfGrid3<L> {
    /// Accepted states of one aperture in canonical z/y/x order, without an
    /// allocation. Uniform and spatial controls have the same inspection path.
    /// `None` means that the aperture index is not attached.
    pub fn reservoir_aperture_tile_states(
        &self,
        aperture: usize,
    ) -> Option<impl ExactSizeIterator<Item = NsfReservoirTileState> + '_> {
        let bound = self.reservoir_apertures.get(aperture)?;
        Some(
            bound
                .indices
                .iter()
                .enumerate()
                .map(move |(position, &index)| NsfReservoirTileState {
                    cell: self.at(index),
                    state: bound
                        .profile
                        .as_ref()
                        .map_or(bound.source.state, |profile| profile.tiles[position].state),
                }),
        )
    }

    /// Atomically replace every prescribed tile state of one existing aperture.
    ///
    /// Every owned tile must appear exactly once. Input order is irrelevant;
    /// duplicates, missing tiles, other owners and out-of-grid coordinates are
    /// refused. The original gas source admits each state and includes its
    /// resolved kinetic energy. No density, speed or temperature is clipped.
    ///
    /// The same original boundary supplies incoming populations, thermal force,
    /// SGS boundary gradients, conservative energy and all exchange accounts.
    /// A warm update rebases the Guo force and energy-flux history together,
    /// preserving accepted physical fluid fields. Refusal retains the complete
    /// previous profile and all accepted solver state. The additional cache is
    /// one bounded vector with exactly one record per existing aperture tile.
    ///
    /// A spatially uniform profile delegates to the existing scalar control,
    /// including its exact arithmetic and replay hash. Scalar controls also
    /// clear any profile on the named owner. Nonuniform controls retain the
    /// aperture's uniform reference state for existing source queries; inspect
    /// `reservoir_aperture_tile_states` for the actual prescribed states.
    pub fn set_reservoir_aperture_profile(
        &mut self,
        aperture: usize,
        states: &[NsfReservoirTileState],
    ) -> Result<(), NsfError> {
        let bound = self
            .reservoir_apertures
            .get(aperture)
            .ok_or(NsfError::Input(
                "reservoir profile requires a valid aperture index",
            ))?;
        if states.len() != bound.indices.len() {
            return Err(NsfError::Input(
                "reservoir profile must cover every owned tile exactly once",
            ));
        }
        let mut tiles = Vec::new();
        tiles
            .try_reserve_exact(states.len())
            .map_err(|_| NsfError::Input("reservoir profile staging allocation failed"))?;
        for tile in states {
            if tile
                .cell
                .iter()
                .zip(self.spec.dims)
                .any(|(&x, size)| x >= size)
            {
                return Err(NsfError::Input(
                    "reservoir profile tile lies outside the grid",
                ));
            }
            let index = self.index(tile.cell);
            if bound.indices.binary_search(&index).is_err() {
                return Err(NsfError::Input(
                    "reservoir profile tile is not owned by this aperture",
                ));
            }
            tiles.push(ReservoirTile {
                index,
                state: tile.state,
                ghost: self.admit_open_state(tile.state)?,
            });
        }
        tiles.sort_unstable_by_key(|tile| tile.index);
        if !tiles
            .iter()
            .zip(&bound.indices)
            .all(|(tile, &index)| tile.index == index)
        {
            return Err(NsfError::Input(
                "reservoir profile has duplicate or missing tiles",
            ));
        }
        let first = tiles[0].state; // Existing apertures are nonempty at admission.
        if tiles.iter().all(|tile| tile.state == first) {
            return self.set_reservoir_aperture_states(&[(aperture, first)]);
        }
        if bound.profile.as_ref().is_some_and(|previous| {
            previous
                .tiles
                .iter()
                .zip(&tiles)
                .all(|(old, new)| old.state == new.state)
        }) {
            return Ok(());
        }
        let drives = self.boundary_drives().ok_or(NsfError::Input(
            "reservoir profile requires its original characteristic axis",
        ))?;
        let previous = self.reservoir_apertures[aperture]
            .profile
            .replace(ReservoirProfile { tiles });
        if let Err(error) = self.set_boundary_drives(drives) {
            // The original owner publishes no fields on rebind failure.
            self.reservoir_apertures[aperture].profile = previous;
            return Err(error);
        }
        Ok(())
    }
}
