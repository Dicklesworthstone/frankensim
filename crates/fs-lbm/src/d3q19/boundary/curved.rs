//! Sparse off-lattice no-slip boundaries on the original D3Q19 pull operator.
//!
//! Linear Bouzidi-Firdaouss-Lallemand interpolation (Phys. Fluids 13,
//! 3452, 2001, doi:10.1063/1.1399290) uses the actual fluid-to-wall fraction.
//! Its returned population is also the population used for momentum exchange.
//! Geometry and occupancy are fixed during a step; this is not a fresh-cell
//! remapper or a displacement/volume-conservation model.

use super::partial_saturation::{cross, low_mach};
use super::{
    BoundaryGrid3, BoundaryLink3, BoundaryStepError3, CancelGate, E3, Face3, FaceBoundary3,
    NO_FLUID_SOURCE, OPP3, PeriodicMomentumTransport3, PullSource, Q3, StreamView, TILE_CELLS,
    WallMomentumExchange3, cell_force, effective_wall_velocity, wall_population,
};

type Result<T> = core::result::Result<T, BoundaryStepError3>;

/// An actual intersection on a native wall link. Pull direction points from
/// wall to fluid, so the intersection is cell_center - fraction * c_q.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurvedWallLink3 {
    /// Original destination cell and incoming direction of the wall link.
    pub link: BoundaryLink3,
    /// Fraction of the lattice link lying in fluid: strictly positive, <= 1.
    pub fraction: f64,
}

impl CurvedWallLink3 {
    /// Intersection in the grid's inertial coordinates (dx = 1).
    #[must_use]
    pub fn intersection(self) -> [f64; 3] {
        let c = E3[self.link.direction];
        let e = [c.0, c.1, c.2];
        core::array::from_fn(|i| self.link.cell[i] as f64 + 0.5 - self.fraction * f64::from(e[i]))
    }
}

#[derive(Debug, Clone, Copy)]
struct CompiledLink3 {
    geometry: CurvedWallLink3,
    /// Opposite pull's native fluid source for the fraction < 1/2 branch.
    far: Option<u32>,
}

/// Reusable ordered boundary geometry, with cached short-link fluid donors.
///
/// Only the selected links are replaced; all other voxel/exterior boundaries
/// retain their original native rule. One binding can include several surfaces
/// with independently prescribed intersection velocities. Its exchange is their
/// sum. Use the same geometry every coupled step; ordinary BoundaryGrid3::step does not
/// retain an implicit boundary override.
///
/// Construction costs O(selected links), apart from the optional level-set
/// generator's one-time occupancy scan. Stepping uses the original post/next
/// buffers and a sparse cursor, with no dense wall-velocity or wall-distance
/// field. Geometry is checked against the current pull table before collision.
#[derive(Debug, Clone)]
pub struct CurvedWall3 {
    dims: [usize; 3],
    links: Vec<CompiledLink3>,
}

impl CurvedWall3 {
    /// Bind a unique canonical (tile, lane, direction) selection of native wall
    /// links. An empty selection is valid and recovers the ordinary closed step.
    /// A short link requires a real opposite-side fluid donor; no wall, solid,
    /// open reservoir or unresolved one-cell gap is substituted for that donor.
    pub fn new(grid: &BoundaryGrid3, links: &[CurvedWallLink3], gate: &CancelGate) -> Result<Self> {
        checkpoint(gate)?;
        let mut prepared = Vec::new();
        prepared
            .try_reserve_exact(links.len())
            .map_err(|_| invalid("curved wall geometry allocation refused"))?;
        let mut previous = None;
        for &geometry in links {
            checkpoint(gate)?;
            let key = key(geometry.link);
            if previous.is_some_and(|old| old >= key) {
                return Err(invalid(
                    "curved wall links must be unique and in canonical order",
                ));
            }
            let far = admit_link(grid, geometry)?;
            prepared.push(CompiledLink3 { geometry, far });
            previous = Some(key);
        }
        checkpoint(gate)?;
        Ok(Self {
            dims: grid.dimensions(),
            links: prepared,
        })
    }

    /// Generate every voxel-solid intersection from the same signed level set
    /// that defines occupancy: positive in fluid, nonpositive in solid.
    ///
    /// The callback must be deterministic, finite, continuous on each selected
    /// link and have one interface crossing on it. Those geometric regularity
    /// assumptions belong to the caller; bisection is not an SDF certificate.
    /// All cell-center signs are checked against the grid before binding.
    /// Exterior domain walls keep their native rule. For periodic solid images
    /// the level set must also be periodic at queried unwrapped link endpoints.
    ///
    /// Each root is bracketed within fractional_tolerance in lattice-link
    /// parameter, using at most max_iterations (1..=64). Insufficient budget,
    /// inconsistent occupancy, no short-link donor or cancellation refuses.
    pub fn from_level_set(
        grid: &BoundaryGrid3,
        mut level_set: impl FnMut([f64; 3]) -> f64,
        max_iterations: usize,
        fractional_tolerance: f64,
        gate: &CancelGate,
    ) -> Result<Self> {
        checkpoint(gate)?;
        if !(1..=64).contains(&max_iterations)
            || !fractional_tolerance.is_finite()
            || fractional_tolerance <= 0.0
            || fractional_tolerance >= 0.5
        {
            return Err(invalid(
                "bounded root iterations and a positive sub-half-link tolerance required",
            ));
        }
        for tile in 0..grid.solid.len() {
            checkpoint(gate)?;
            for lane in 0..TILE_CELLS {
                let cell = grid.coords(tile, lane);
                let value = level_set(cell.map(|x| x as f64 + 0.5));
                let solid = grid.solid[tile] & (1u64 << lane) != 0;
                if !value.is_finite() || (value <= 0.0) != solid {
                    return Err(invalid(
                        "level-set cell signs differ from native solid occupancy",
                    ));
                }
            }
        }
        let mut links = Vec::new();
        for tile in 0..grid.solid.len() {
            checkpoint(gate)?;
            for lane in 0..TILE_CELLS {
                if grid.solid[tile] & (1u64 << lane) != 0 {
                    continue;
                }
                let cell = grid.coords(tile, lane);
                let mask = grid.link_masks[tile].0[lane];
                for q in 1..Q3 {
                    if mask & (1u32 << q) == 0
                        || !matches!(
                            grid.classify_pull(&grid.solid, cell[0], cell[1], cell[2], q),
                            PullSource::Wall { voxel: true }
                        )
                    {
                        continue;
                    }
                    checkpoint(gate)?;
                    let link = BoundaryLink3 {
                        tile,
                        lane,
                        cell,
                        direction: q,
                    };
                    let start = cell.map(|x| x as f64 + 0.5);
                    let c = E3[q];
                    let e = [f64::from(c.0), f64::from(c.1), f64::from(c.2)];
                    let point = |s: f64| core::array::from_fn(|i| start[i] - s * e[i]);
                    let start_value = level_set(start);
                    let end_value = level_set(point(1.0));
                    if !start_value.is_finite()
                        || start_value <= 0.0
                        || !end_value.is_finite()
                        || end_value > 0.0
                    {
                        return Err(invalid(
                            "curved wall root needs a finite fluid/solid bracket",
                        ));
                    }
                    let fraction = if end_value == 0.0 {
                        1.0
                    } else {
                        let (mut lo, mut hi) = (0.0, 1.0);
                        for _ in 0..max_iterations {
                            checkpoint(gate)?;
                            if hi - lo <= fractional_tolerance {
                                break;
                            }
                            let mid = 0.5 * lo + 0.5 * hi;
                            let value = level_set(point(mid));
                            if !value.is_finite() {
                                return Err(invalid("curved wall level set became non-finite"));
                            }
                            if value == 0.0 {
                                lo = mid;
                                hi = mid;
                                break;
                            } else if value > 0.0 {
                                lo = mid;
                            } else {
                                hi = mid;
                            }
                        }
                        if hi - lo > fractional_tolerance {
                            return Err(invalid("curved wall root iteration budget exhausted"));
                        }
                        0.5 * lo + 0.5 * hi
                    };
                    links
                        .try_reserve(1)
                        .map_err(|_| invalid("curved wall geometry allocation refused"))?;
                    links.push(CurvedWallLink3 { link, fraction });
                }
            }
        }
        Self::new(grid, &links, gate)
    }

    /// Bound links, in native canonical order; no allocation.
    pub fn links(&self) -> impl ExactSizeIterator<Item = CurvedWallLink3> + '_ {
        self.links.iter().map(|link| link.geometry)
    }

    fn validate(&self, grid: &BoundaryGrid3, gate: &CancelGate) -> Result<()> {
        checkpoint(gate)?;
        if self.dims != grid.dimensions() {
            return Err(invalid("curved wall binding has different grid dimensions"));
        }
        for link in &self.links {
            checkpoint(gate)?;
            if admit_link(grid, link.geometry)? != link.far {
                return Err(invalid("curved wall donor topology changed after binding"));
            }
        }
        Ok(())
    }
}

/// Actual momentum exchange of the accepted native collide/stream step.
#[derive(Debug, Clone, Copy, PartialEq)]
#[must_use]
pub struct CurvedWallStep3 {
    /// Selected curved links, at their actual off-lattice intersections.
    pub curved: WallMomentumExchange3,
    /// Every remaining native wall link, counted once at its halfway point.
    pub other_walls: WallMomentumExchange3,
    /// Held Guo body-force input, distinct from all wall reactions.
    pub body_force_impulse: [f64; 3],
    /// Moment of the held fluid force about the same requested origin.
    pub body_force_angular_impulse: [f64; 3],
    /// Actual periodic coordinate transport, distinct from surface torque.
    pub periodic: PeriodicMomentumTransport3,
}

impl BoundaryGrid3 {
    /// Advance using native collision and pull with linear interpolated walls.
    ///
    /// velocity is evaluated once per selected link at its actual intersection.
    /// Supply a physical inertial-frame surface velocity, e.g. V + omega cross
    /// (intersection - center). Each speed must satisfy the existing low-Mach
    /// envelope. Geometry, occupancy and the supplied kinematics are frozen for
    /// this step. Evolving solid occupancy requires the existing moving-body
    /// owner or an explicit conservative topology transition; this method does
    /// not insert/remove fluid cells or advance a rigid body's pose.
    ///
    /// All four native collision models and held Guo forces retain their own
    /// original law. Wall/periodic domains are admitted; regularized open faces
    /// refuse because their separate reservoir exchange is not owned here.
    ///
    /// BFL interpolation can change population mass. The returned measured mass
    /// change is not repaired or reclassified as physical swept volume. The
    /// boundary-relative force, actual population impulse and velocity-times-mass
    /// correction are distinct, closing J_wall + dP_fluid = velocity*dM.
    /// Torque uses the same actual intersection. No full-solver Galilean,
    /// geometric-conservation, drag-calibration or fluid-energy claim is made.
    ///
    /// Candidate populations and every load are checked before the final swap.
    /// Any geometry, collision, arithmetic, velocity or cancellation refusal
    /// leaves populations and topology admission unchanged; scratch may change.
    pub fn step_with_curved_walls(
        &mut self,
        wall: &CurvedWall3,
        mut velocity: impl FnMut(CurvedWallLink3) -> [f64; 3],
        moment_origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<CurvedWallStep3> {
        wall.validate(self, gate)?;
        if !moment_origin.into_iter().all(f64::is_finite) {
            return Err(invalid("curved wall moment origin must be finite"));
        }
        if Face3::ALL.into_iter().any(|face| {
            matches!(
                self.boundaries.face(face),
                FaceBoundary3::Velocity { .. } | FaceBoundary3::Pressure { .. }
            )
        }) {
            return Err(invalid("curved wall stepping requires wall/periodic faces"));
        }
        self.try_collide_with_gate(Some(gate))?;
        let periodic = self.periodic_momentum_transport(gate).map_err(|_| {
            if gate.is_requested() {
                BoundaryStepError3::Cancelled
            } else {
                invalid("periodic momentum transport overflowed")
            }
        })?;
        let mut next = self.next.take().unwrap_or_else(|| self.f.clone());
        let outcome = (|| {
            let view = StreamView {
                post: &self.post,
                pull: &self.pull,
                link_masks: &self.link_masks,
                open_link_masks: &self.open_link_masks,
                stationary_link_masks: &self.stationary_link_masks,
                boundaries: self.boundaries,
                dims: self.dimensions(),
            };
            let mut curved = WallMomentumExchange3::default();
            let mut other_walls = WallMomentumExchange3::default();
            let mut body_force_impulse = [0.0; 3];
            let mut body_force_angular_impulse = [0.0; 3];
            let mut cursor = 0;
            for tile in 0..self.f[0].len() {
                checkpoint(gate)?;
                for lane in 0..TILE_CELLS {
                    if self.solid[tile] & (1u64 << lane) != 0 {
                        continue;
                    }
                    let mut values = view.pull_cell(tile, lane);
                    let mask = self.link_masks[tile].0[lane];
                    if mask != 0 {
                        let cell = self.coords(tile, lane);
                        let rho = (0..Q3).map(|q| self.post[q][tile].0[lane]).sum::<f64>();
                        let exterior_velocity =
                            effective_wall_velocity(self.boundaries, self.dimensions(), cell);
                        for (q, value) in values.iter_mut().enumerate() {
                            if mask & (1u32 << q) == 0 {
                                continue;
                            }
                            let link = BoundaryLink3 {
                                tile,
                                lane,
                                cell,
                                direction: q,
                            };
                            let outgoing = self.post[OPP3[q]][tile].0[lane];
                            if let Some(compiled) = wall
                                .links
                                .get(cursor)
                                .filter(|l| key(l.geometry.link) == key(link))
                            {
                                checkpoint(gate)?;
                                let geometry = compiled.geometry;
                                let u = velocity(geometry);
                                if !low_mach(u) {
                                    return Err(invalid(
                                        "curved wall velocity is non-finite or outside low Mach",
                                    ));
                                }
                                let half = wall_population(outgoing, rho, q, u);
                                *value = if geometry.fraction == 0.5 {
                                    half
                                } else if geometry.fraction < 0.5 {
                                    let source =
                                        compiled.far.expect("admitted short-link donor") as usize;
                                    let far = self.post[OPP3[q]][source / TILE_CELLS].0
                                        [source % TILE_CELLS];
                                    // Keep the moving correction separate: subtracting it from
                                    // half would lose a small imposed velocity beside outgoing.
                                    let correction = wall_population(0.0, rho, q, u);
                                    2.0 * geometry.fraction * outgoing
                                        + (1.0 - 2.0 * geometry.fraction) * far
                                        + correction
                                } else {
                                    let width = 2.0 * geometry.fraction;
                                    half / width
                                        + ((width - 1.0) / width) * self.post[q][tile].0[lane]
                                };
                                curved.add_population_link(
                                    link,
                                    outgoing,
                                    *value,
                                    u,
                                    geometry.fraction,
                                    moment_origin,
                                );
                                cursor += 1;
                            } else {
                                let u = if self.stationary_link_masks[tile].0[lane] & (1u32 << q)
                                    != 0
                                {
                                    [0.0; 3]
                                } else {
                                    exterior_velocity
                                };
                                other_walls.add_population_link(
                                    link,
                                    outgoing,
                                    *value,
                                    u,
                                    0.5,
                                    moment_origin,
                                );
                            }
                        }
                    }
                    let rho = values.iter().sum::<f64>();
                    if !values.iter().all(|x| x.is_finite()) || !(rho.is_finite() && rho > 0.0) {
                        return Err(BoundaryStepError3::Unphysical {
                            cell: self.coords(tile, lane),
                        });
                    }
                    if self.force != [0.0; 3] || self.force_field.is_some() {
                        let force = cell_force(self.force, self.force_field.as_deref(), tile, lane);
                        let cell = self.coords(tile, lane);
                        let arm = core::array::from_fn(|i| cell[i] as f64 + 0.5 - moment_origin[i]);
                        let torque = cross(arm, force);
                        for i in 0..3 {
                            body_force_impulse[i] += force[i];
                            body_force_angular_impulse[i] += torque[i];
                        }
                    }
                    for (field, value) in next.iter_mut().zip(values) {
                        field[tile].0[lane] = value;
                    }
                }
            }
            debug_assert_eq!(cursor, wall.links.len());
            let curved = curved.checked()?;
            let other_walls = other_walls.checked()?;
            if !body_force_impulse
                .into_iter()
                .chain(body_force_angular_impulse)
                .all(f64::is_finite)
            {
                return Err(invalid("curved wall body-force account overflowed"));
            }
            checkpoint(gate)?;
            Ok(CurvedWallStep3 {
                curved,
                other_walls,
                body_force_impulse,
                body_force_angular_impulse,
                periodic,
            })
        })();
        if outcome.is_ok() {
            core::mem::swap(&mut self.f, &mut next);
            self.topology_locked = true;
        }
        self.next = Some(next);
        outcome
    }
}

fn key(link: BoundaryLink3) -> (usize, usize, usize) {
    (link.tile, link.lane, link.direction)
}

fn admit_link(grid: &BoundaryGrid3, geometry: CurvedWallLink3) -> Result<Option<u32>> {
    let link = geometry.link;
    if !geometry.fraction.is_finite()
        || geometry.fraction <= 0.0
        || geometry.fraction > 1.0
        || link.tile >= grid.solid.len()
        || link.lane >= TILE_CELLS
        || link.direction == 0
        || link.direction >= Q3
        || grid.coords(link.tile, link.lane) != link.cell
        || grid.solid[link.tile] & (1u64 << link.lane) != 0
        || grid.link_masks[link.tile].0[link.lane] & (1u32 << link.direction) == 0
    {
        return Err(invalid(
            "curved wall requires a native wall link and fraction in (0,1]",
        ));
    }
    if geometry.fraction < 0.5 {
        let donor = grid.pull[(link.tile * TILE_CELLS + link.lane) * Q3 + OPP3[link.direction]];
        if donor == NO_FLUID_SOURCE {
            return Err(invalid(
                "short curved wall link requires an opposite-side fluid donor",
            ));
        }
        Ok(Some(donor))
    } else {
        Ok(None)
    }
}

fn checkpoint(gate: &CancelGate) -> Result<()> {
    if gate.is_requested() {
        Err(BoundaryStepError3::Cancelled)
    } else {
        Ok(())
    }
}

fn invalid(reason: &'static str) -> BoundaryStepError3 {
    BoundaryStepError3::WallMomentum { reason }
}
