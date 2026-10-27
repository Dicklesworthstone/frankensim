//! Conservative initial mass and velocity ingress into the original VOF owner.

use fs_exec::Cx;
use std::collections::{BTreeMap, BTreeSet};

use super::*;
use crate::MACH_LIMIT;
use crate::d3q19::sparse::MORTON_COORD_BITS;

/// One initially occupied Cartesian cell. Omitted cells contain no liquid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceInitialCell3 {
    /// Global cell coordinates, with no duplicate entries.
    pub coordinates: [usize; 3],
    /// Liquid volume fraction in (0, 1], retained at boundary cells as well.
    pub fill_fraction: f64,
    /// Positive lattice density. This is a numerical pressure state, not an EOS.
    pub density: f64,
    /// Physical lattice velocity, including the Guo half-force convention.
    pub velocity: [f64; 3],
}

fn refused(coordinates: Option<[usize; 3]>, reason: &'static str) -> FreeSurfaceError3 {
    FreeSurfaceError3::InvalidInitialState {
        coordinates,
        reason,
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), FreeSurfaceError3> {
    cx.checkpoint().map_err(|_| FreeSurfaceError3::Cancelled)
}

impl FreeSurface3 {
    /// Construct the original sparse free surface from transferred cell state.
    ///
    /// Unlike the binary occupancy constructor, this path never invents a
    /// half-filled shell or resets transferred velocity to rest. A full cell
    /// bordering gas is an interface cell with fill ONE; interior partial cells
    /// also remain interface cells. Initial mass is rho*fill (to floating-point
    /// precision), with zero redistribution carry. Populations are equilibrium
    /// at `velocity - force/(2*rho)` so native sampling returns the supplied
    /// Guo-corrected velocity. This projects stress to equilibrium; it does not
    /// transfer unresolved stress, caloric state or subcell velocity moments.
    ///
    /// Input order does not affect allocation, initialization or later sweeps.
    /// Only supplied occupied tiles and their one-tile gas margin are activated.
    /// The caller owns a domain/allocation budget, as for the original native
    /// constructor; the material adapter enforces its explicit cell budget.
    ///
    /// # Errors
    /// Refuses invalid dimensions/solid tiles, duplicate or out-of-domain
    /// cells, wall intersections, nonpositive/nonfinite density, invalid fill,
    /// nonfinite/excessive initial Mach, inadmissible relaxation/capillarity,
    /// nonpositive equilibrium populations, or observed cancellation. The
    /// candidate is returned only after all state and closure checks succeed.
    #[allow(clippy::too_many_arguments)]
    pub fn from_cells(
        cells: [usize; 3],
        tau: f64,
        gravity: [f64; 3],
        sigma: f64,
        contact: ContactModel3,
        solid_tiles: &[(u32, u32, u32)],
        initial: &[FreeSurfaceInitialCell3],
        cx: &Cx<'_>,
    ) -> Result<Self, FreeSurfaceError3> {
        checkpoint(cx)?;
        let [nx, ny, nz] = cells;
        let mut grid = SparseGrid3::new(nx, ny, nz, tau, gravity)?;
        if cells
            .iter()
            .any(|&n| n / TILE > (1usize << MORTON_COORD_BITS))
            || cells
                .iter()
                .try_fold(1usize, |n, &m| n.checked_mul(m))
                .is_none()
        {
            return Err(refused(None, "domain exceeds addressable Cartesian cells"));
        }
        if !tau.is_finite() || tau <= 0.5 || gravity.iter().any(|g| !g.is_finite()) {
            return Err(refused(None, "relaxation or force is inadmissible"));
        }
        if !sigma.is_finite() || sigma < 0.0 || sigma * Self::CURVATURE_LIMIT / CS2_3 >= 1.0 {
            return Err(refused(
                None,
                "capillary boundary density can be nonpositive",
            ));
        }
        let tile_dims = cells.map(|n| n / TILE);
        let mut solid = BTreeSet::new();
        for &(tx, ty, tz) in solid_tiles {
            checkpoint(cx)?;
            if [tx, ty, tz]
                .iter()
                .enumerate()
                .any(|(axis, &t)| t as usize >= tile_dims[axis])
            {
                return Err(FreeSurfaceError3::Grid(SparseError3::TileOutOfDomain {
                    tx,
                    ty,
                    tz,
                }));
            }
            solid.insert(morton3(tx, ty, tz));
        }

        let mut states = BTreeMap::new();
        let mut active = BTreeSet::new();
        for state in initial {
            checkpoint(cx)?;
            let c = state.coordinates;
            if (0..3).any(|axis| c[axis] >= cells[axis]) {
                return Err(refused(Some(c), "cell lies outside the domain"));
            }
            let tile = c.map(|v| v / TILE);
            let key = morton3(tile[0] as u32, tile[1] as u32, tile[2] as u32);
            let lane = c[0] % TILE + TILE * (c[1] % TILE + TILE * (c[2] % TILE));
            if solid.contains(&key) {
                return Err(refused(Some(c), "liquid intersects a permanent solid tile"));
            }
            if !state.fill_fraction.is_finite()
                || state.fill_fraction <= 0.0
                || state.fill_fraction > 1.0
                || !state.density.is_finite()
                || state.density <= 0.0
                || !(state.fill_fraction * state.density).is_finite()
                || state.fill_fraction * state.density <= 0.0
            {
                return Err(refused(
                    Some(c),
                    "fill, density or cell mass is inadmissible",
                ));
            }
            let speed = state.velocity[0]
                .hypot(state.velocity[1])
                .hypot(state.velocity[2]);
            if !speed.is_finite() || speed / CS2_3.sqrt() >= MACH_LIMIT {
                return Err(refused(
                    Some(c),
                    "initial velocity exceeds the lattice Mach domain",
                ));
            }
            let raw_velocity = core::array::from_fn(|axis| {
                state.velocity[axis] - 0.5 * gravity[axis] / state.density
            });
            let populations = equilibrium3(state.density, raw_velocity);
            if populations.iter().any(|f| !f.is_finite() || *f <= 0.0)
                || !populations.iter().sum::<f64>().is_finite()
            {
                return Err(refused(
                    Some(c),
                    "initial equilibrium populations are not positive",
                ));
            }
            if states.insert((key, lane), (*state, populations)).is_some() {
                return Err(refused(Some(c), "duplicate occupied cell"));
            }
            for dz in -1i64..=1 {
                for dy in -1i64..=1 {
                    for dx in -1i64..=1 {
                        let a = [
                            tile[0] as i64 + dx,
                            tile[1] as i64 + dy,
                            tile[2] as i64 + dz,
                        ];
                        if (0..3).all(|axis| a[axis] >= 0 && (a[axis] as usize) < tile_dims[axis]) {
                            let coords = (a[0] as u32, a[1] as u32, a[2] as u32);
                            if !solid.contains(&morton3(coords.0, coords.1, coords.2)) {
                                active.insert(coords);
                            }
                        }
                    }
                }
            }
        }
        if states.is_empty() {
            return Err(FreeSurfaceError3::NoFluid);
        }
        checkpoint(cx)?;
        grid.activate_tiles(&active.into_iter().collect::<Vec<_>>())?;
        let mut tiles = BTreeMap::new();
        for &key in grid.active_keys() {
            checkpoint(cx)?;
            tiles.insert(key, FsTile::all_gas());
        }
        for ((key, lane), (state, populations)) in states {
            checkpoint(cx)?;
            let slot = grid.slot_of(key).expect("supplied occupied tile is active");
            grid.set_populations(slot, lane, populations);
            let tile = tiles.get_mut(&key).expect("active initial tile");
            tile.cells[lane] = if state.fill_fraction == 1.0 {
                Cell3::Fluid
            } else {
                Cell3::Interface
            };
            tile.mass[lane] = state.fill_fraction * state.density;
        }
        let mut result = Self {
            grid,
            tiles,
            solid,
            sigma,
            local_sigma: None,
            contact,
            planar_contact: None,
            planar_contact_field: None,
            carry: 0.0,
            conversions: ConversionStats3::default(),
            worst_step_violation: 0.0,
            fill_smooth: BTreeMap::new(),
            local_force: None,
        };
        let mut interfaces = Vec::new();
        for (&key, tile) in &result.tiles {
            checkpoint(cx)?;
            for lane in 0..TILE_CELLS {
                if tile.cells[lane] == Cell3::Fluid
                    && (1..Q3).any(|q| {
                        matches!(
                            result.neighbor_state(key, lane, q),
                            Neighbor::Active(_, _, Cell3::Gas)
                        )
                    })
                {
                    interfaces.push((key, lane));
                }
            }
        }
        for (key, lane) in interfaces {
            result
                .tiles
                .get_mut(&key)
                .expect("active boundary tile")
                .cells[lane] = Cell3::Interface;
        }
        checkpoint(cx)?;
        result.assert_closure();
        if !result.ledger_mass().is_finite() {
            return Err(refused(None, "total initial mass is unrepresentable"));
        }
        checkpoint(cx)?;
        Ok(result)
    }
}
