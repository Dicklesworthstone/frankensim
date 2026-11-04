//! Conservative initial mass and velocity ingress into the original VOF owner.

use fs_exec::Cx;
use std::collections::{BTreeMap, BTreeSet};

use super::*;
use crate::MACH_LIMIT;
use crate::d3q19::sparse::MORTON_COORD_BITS;

#[cfg(test)]
mod stress_tests;

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
    /// Seed one cell of an unpublished transfer candidate with its admitted
    /// symmetric strain, packed xx, yy, zz, xy, xz, yz, in lattice time units.
    /// The supplied molecular tau belongs to that SAME strain/source query.
    /// This changes no represented mass, velocity, heat or elapsed time. It
    /// supplies Newtonian/generalized-Newtonian stress, never elastic memory.
    /// Callers own cancellation between these fixed-size cell operations and
    /// re-audit actual native momenta before returning their complete candidate.
    pub(crate) fn initialize_strain_at(
        &mut self,
        coordinates: [usize; 3],
        strain: [f64; 6],
        tau: f64,
    ) -> Result<(), FreeSurfaceError3> {
        let fail = |reason| refused(Some(coordinates), reason);
        if !tau.is_finite() || tau <= 0.5 || strain.iter().any(|v| !v.is_finite()) {
            return Err(fail("initial viscous strain or relaxation is inadmissible"));
        }
        let [x, y, z] = coordinates;
        let sample = self.sample_at(x, y, z)
            .filter(|state| matches!(state.cell, Cell3::Fluid | Cell3::Interface)
                && state.fill_fraction > 0.0)
            .ok_or_else(|| fail("initial viscous strain requires an occupied liquid cell"))?;
        // Retain the exact existing equilibrium/Guo representation at zero
        // strain, including rigid-motion input and its density reduction.
        if strain == [0.0; 6] {
            return Ok(());
        }
        let key = morton3((x / TILE) as u32, (y / TILE) as u32, (z / TILE) as u32);
        let lane = x % TILE + TILE * (y % TILE + TILE * (z % TILE));
        let slot = self.grid.slot_of(key)
            .ok_or_else(|| fail("initial viscous strain lost its native tile"))?;
        let rho = sample.density;
        let velocity = sample.velocity;
        let force = self.force_at_lane(key, lane);
        // Pi_neq = -2 rho cs^2 tau D. The physical viscous Cauchy stress is
        // -(1 - 1/(2 tau))*Pi_neq; using tau-1/2 here would seed the wrong rate.
        let expected = strain.map(|rate| -2.0 * CS2_3 * rho * tau * rate);
        if expected.iter().any(|v| !v.is_finite()) {
            return Err(fail("initial viscous stress overflows lattice units"));
        }
        let [xx, yy, zz, xy, xz, yz] = expected;
        let tensor = [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]];
        let mut populations = crate::d3q19::boundary::regularized_populations(rho, velocity, tensor);
        forcing::preserve_density(&mut populations, rho)?;
        let populations = forcing::reframe(populations, velocity, [0.0; 3], force)?;
        if populations.iter().any(|v| !v.is_finite() || *v <= 0.0)
            || populations.iter().sum::<f64>() != rho
        {
            return Err(fail("initial viscous stress needs positive density-preserving populations"));
        }
        let actual_velocity: [f64; 3] = core::array::from_fn(|axis| {
            let momentum = populations.iter().enumerate().map(|(q, &f)| {
                let direction = [E3[q].0, E3[q].1, E3[q].2];
                f * f64::from(direction[axis])
            }).sum::<f64>();
            (momentum + 0.5 * force[axis]) / rho
        });
        let actual = crate::d3q19::force_corrected_stress3(
            &populations, rho, actual_velocity, force,
        );
        let allowance = 256.0 * f64::EPSILON
            * (rho + expected.iter().map(|v| v.abs()).sum::<f64>());
        if !allowance.is_finite()
            || actual_velocity.iter().zip(velocity).any(|(actual, expected)| {
                !actual.is_finite() || ((actual - expected) * rho).abs() > allowance
            })
            || actual.iter().zip(expected).any(|(actual, expected)| {
                !actual.is_finite() || (actual - expected).abs() > allowance
            })
        {
            return Err(fail("initial viscous populations do not reproduce physical moments"));
        }
        self.grid.set_populations(slot, lane, populations);
        Ok(())
    }

    /// Construct the original sparse free surface from transferred cell state.
    ///
    /// Unlike the binary occupancy constructor, this path never invents a
    /// half-filled shell or resets transferred velocity to rest. A full cell
    /// bordering gas is an interface cell with fill ONE; interior partial cells
    /// also remain interface cells. Initial mass is rho*fill (to floating-point
    /// precision), with zero redistribution carry. Populations are physical
    /// equilibrium minus the Guo half source `S(velocity, force)/2`, so native
    /// sampling returns the supplied velocity without an artificial force-squared
    /// pressure/stress. This projects physical stress to equilibrium; it does not
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
        if !crate::admit_explicit_capillary_coefficient(sigma) {
            return Err(refused(
                None,
                "surface tension exceeds conservative explicit capillary-wave bound",
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
            if !crate::admit_explicit_capillary_density(sigma, state.density) {
                return Err(refused(
                    Some(c),
                    "surface tension exceeds conservative explicit capillary-wave bound at initial density",
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
            let populations = forcing::reframe(
                equilibrium3(state.density, state.velocity),
                state.velocity,
                [0.0; 3],
                gravity,
            )?;
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
            redistribution: ConversionRedistribution3::Local,
            contact,
            planar_contact: None,
            planar_contact_field: None,
            planar_piston: None,
            geometry_epoch: 0,
            carry: 0.0,
            conversions: ConversionStats3::default(),
            conversion_momentum: None,
            worst_step_violation: 0.0,
            fill_smooth: BTreeMap::new(),
            local_force: None,
            body_acceleration: None,
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
