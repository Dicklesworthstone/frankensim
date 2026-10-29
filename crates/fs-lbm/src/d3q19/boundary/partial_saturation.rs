//! Spatial rigid-body response for the native partially saturated BGK law.
//!
//! With `G(r) v = V + omega cross r`, the solid impulse at one cell is
//! `B rho (u-Gv)`. Thus `K = sum(B rho G^T G)` is positive semidefinite and
//! one midpoint impulse solve `(M+K/2) delta_v = sum(B rho G^T(u-Gv_n))+J`
//! couples all six velocities. The inertia and sampled geometry are frozen
//! in this substep: no Euler/gyroscopic or pose evolution is hidden here.
//!
//! The Noble-Torczynski non-equilibrium reflection law is the same one used
//! by D2Q9 (`doi:10.1142/S0129183198001084`). Fluid populations remain inside
//! the body, so its inertia is ADDITIONAL to the retained fluid. This is not
//! a displaced-fluid/physical-density model.

use super::{
    BoundaryGrid3, BoundaryLink3, CancelGate, CollisionError3, CollisionModel3, E3, Face3,
    FaceBoundary3, OPP3, Q3, StreamView, TILE_CELLS, WallMomentumExchange3, collide_cell3,
    equilibrium3,
};

/// Frozen spatial rigid-body state in lattice units (`dx = dt = 1`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PartialSaturationBody3 {
    /// Center of mass; lattice cell centers are `[x+0.5,y+0.5,z+0.5]`.
    pub center: [f64; 3],
    /// Positive mass of the additional solid velocity degree of freedom.
    pub mass: f64,
    /// Symmetric positive-definite inertia about the center, in WORLD axes.
    /// The caller must rotate this tensor when advancing body orientation.
    pub inertia_world: [[f64; 3]; 3],
    /// World translational velocity.
    pub velocity: [f64; 3],
    /// World angular velocity.
    pub angular_velocity: [f64; 3],
}

/// The measured correction to ordinary BGK, before streaming.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct PartialSaturationExchange3 {
    /// Population mass introduced by the immersed collision (roundoff only).
    pub fluid_population_mass_change: f64,
    /// Actual population momentum correction, not an analytic force estimate.
    pub fluid_population_impulse: [f64; 3],
    /// Actual population angular correction about the requested origin.
    pub fluid_population_angular_impulse: [f64; 3],
    /// Equal-and-opposite hydrodynamic impulse on the body.
    pub body_impulse: [f64; 3],
    /// Body angular impulse about the requested origin, including its arm.
    pub body_angular_impulse: [f64; 3],
    /// Same measured impulse reduced about the center without subtracting
    /// nearly equal torques about a distant origin.
    pub body_angular_impulse_about_center: [f64; 3],
    /// Hydrodynamic impulse work at the solved rigid midpoint velocity.
    pub body_work: f64,
    /// Number of cells with positive solid fraction.
    pub coupled_cells: usize,
    /// Sum of the Noble-Torczynski coupling weights `B`.
    pub coupling_weight_sum: f64,
}

/// Accepted body update plus separately accounted stationary-wall exchange.
#[derive(Debug, Clone, Copy, PartialEq)]
#[must_use]
pub struct PartialSaturationBodyExchange3 {
    /// Immersed-body collision exchange.
    pub hydrodynamic: PartialSaturationExchange3,
    /// Exchange with all stationary voxel and exterior wall links, measured
    /// from the actual blended populations about the same requested origin.
    /// Periodic transport angular flux is not included in this wall receipt.
    pub walls: WallMomentumExchange3,
    /// Translational velocity used by the fluid collision.
    pub velocity_midpoint: [f64; 3],
    /// Angular velocity used by the fluid collision.
    pub angular_velocity_midpoint: [f64; 3],
    /// Body kinetic-energy change using the accepted, rounded velocities.
    pub kinetic_energy_change: f64,
    /// Work of the supplied linear impulse and center couple at midpoint.
    pub external_work: f64,
    /// Actual body momentum change minus external and measured hydrodynamic
    /// impulses: translation first, then angular momentum about its center.
    pub momentum_residual: [f64; 6],
}

/// A rejected spatial coupling publishes neither fluid nor body velocities.
#[derive(Debug, Clone, PartialEq)]
pub enum PartialSaturationError3 {
    /// Cancellation checked between lattice tiles and before publication.
    Cancelled,
    /// Unsupported or malformed model/body/measurement input.
    Invalid {
        /// Stable diagnostic.
        reason: &'static str,
    },
    /// Invalid sampled fraction, point speed, density, or arithmetic.
    Cell {
        /// Global cell coordinate.
        cell: [usize; 3],
        /// Stable diagnostic.
        reason: &'static str,
    },
    /// Refusal from the unchanged native BGK cell owner.
    Collision {
        /// Global cell coordinate.
        cell: [usize; 3],
        /// Native diagnostic.
        source: CollisionError3,
    },
    /// Refusal from the native body-frame rotation owner.
    Motion {
        /// Native group or time-integration diagnostic.
        source: fs_time::lie::LieStepError,
    },
}

impl core::fmt::Display for PartialSaturationError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Cancelled => write!(f, "D3Q19 immersed-body step cancelled"),
            Self::Invalid { reason } => write!(f, "D3Q19 immersed-body refusal: {reason}"),
            Self::Cell { cell, reason } => write!(f, "D3Q19 immersed-body cell {cell:?}: {reason}"),
            Self::Collision { cell, source } => {
                write!(f, "D3Q19 collision at {cell:?}: {source:?}")
            }
            Self::Motion { source } => write!(f, "D3Q19 body motion: {source}"),
        }
    }
}

impl std::error::Error for PartialSaturationError3 {}

type Result<T> = core::result::Result<T, PartialSaturationError3>;

pub(super) fn invalid(reason: &'static str) -> PartialSaturationError3 {
    PartialSaturationError3::Invalid { reason }
}

pub(super) fn checkpoint(gate: &CancelGate) -> Result<()> {
    if gate.is_requested() {
        Err(PartialSaturationError3::Cancelled)
    } else {
        Ok(())
    }
}

pub(super) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1].mul_add(b[2], -a[2] * b[1]),
        a[2].mul_add(b[0], -a[0] * b[2]),
        a[0].mul_add(b[1], -a[1] * b[0]),
    ]
}

pub(super) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}

pub(super) fn point_velocity(velocity: [f64; 3], omega: [f64; 3], r: [f64; 3]) -> [f64; 3] {
    let rotation = cross(omega, r);
    core::array::from_fn(|axis| velocity[axis] + rotation[axis])
}

pub(super) fn low_mach(velocity: [f64; 3]) -> bool {
    velocity.into_iter().all(f64::is_finite) && dot(velocity, velocity) < 0.03
}

fn moments(populations: [f64; Q3], cell: [usize; 3]) -> Result<(f64, [f64; 3])> {
    let mut rho = 0.0;
    let mut momentum = [0.0; 3];
    for (q, population) in populations.into_iter().enumerate() {
        rho += population;
        momentum[0] += f64::from(E3[q].0) * population;
        momentum[1] += f64::from(E3[q].1) * population;
        momentum[2] += f64::from(E3[q].2) * population;
    }
    let velocity = momentum.map(|p| p / rho);
    if !populations.into_iter().all(f64::is_finite)
        || !(rho.is_finite() && rho > 0.0)
        || !velocity.into_iter().all(f64::is_finite)
    {
        return Err(PartialSaturationError3::Cell {
            cell,
            reason: "positive finite density and finite populations/velocity required",
        });
    }
    Ok((rho, velocity))
}

impl BoundaryGrid3 {
    /// Couple one six-DOF rigid velocity to the native D3Q19 fluid, then stream.
    ///
    /// `solid_fractions[(z*ny+y)*nx+x]` lies in `[0,1]`; positive fractions
    /// may only select fluid cells. `external_impulse` contains three linear
    /// components applied at the body's center, then three couple components.
    /// Geometry, center and world inertia stay fixed; this API advances only
    /// velocities. Pose drift and gyroscopic evolution belong to the caller's
    /// rigid-body integrator. Inertia is additional to retained interior fluid.
    ///
    /// Uses one 6-by-6 SPD solve, no nonlinear fluid/body iteration, and the
    /// native tile-major collision and pull operators. The existing `post`
    /// and `next` buffers are reused; only first-use streaming scratch and two
    /// fixed-size Cholesky factors allocate. Zero fractions recover the
    /// ordinary force-free BGK step bit for bit.
    ///
    /// # Errors
    /// Requires BGK, zero fluid forcing, periodic or stationary-wall faces,
    /// symmetric positive-definite inertia, finite inputs, and body center/
    /// sampled point speeds squared below `0.03` before/midpoint/after.
    /// Invalid inputs, numerical overflow and cancellation leave all fluid
    /// populations, topology admission and body state unchanged. Scratch may
    /// change. Fractions are immutable for the step. No moving geometry,
    /// physical-density correction, overlapping bodies, or free surface claim.
    #[allow(clippy::too_many_lines)] // One atomic coupled transition and publication.
    pub fn step_with_partial_saturation_body(
        &mut self,
        solid_fractions: &[f64],
        body: &mut PartialSaturationBody3,
        external_impulse: [f64; 6],
        moment_origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<PartialSaturationBodyExchange3> {
        self.step_with_partial_saturation_body_checked(
            solid_fractions,
            body,
            external_impulse,
            moment_origin,
            gate,
            |_, _, _| Ok(()),
        )
        .map(|(exchange, ())| exchange)
    }

    // The motion owner completes its final drift, admission and diagnostics
    // here, while every population still belongs to the previous state.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(super) fn step_with_partial_saturation_body_checked<T>(
        &mut self,
        solid_fractions: &[f64],
        body: &mut PartialSaturationBody3,
        external_impulse: [f64; 6],
        moment_origin: [f64; 3],
        gate: &CancelGate,
        finish: impl FnOnce(
            &Self,
            &PartialSaturationBody3,
            &PartialSaturationBodyExchange3,
        ) -> Result<T>,
    ) -> Result<(PartialSaturationBodyExchange3, T)> {
        checkpoint(gate)?;
        let CollisionModel3::Bgk { tau } = self.collision_model else {
            return Err(invalid("finite-inertia immersed coupling requires BGK"));
        };
        if solid_fractions.len() != self.nx * self.ny * self.nz {
            return Err(invalid("solid fractions must cover the whole grid"));
        }
        if self.force != [0.0; 3]
            || Face3::ALL.into_iter().any(|face| {
                !matches!(self.boundaries.face(face), FaceBoundary3::Periodic)
                    && !matches!(self.boundaries.face(face), FaceBoundary3::Wall { velocity } if velocity == [0.0; 3])
            })
        {
            return Err(invalid("requires force-free fluid and periodic/stationary-wall faces"));
        }
        if !(body.mass.is_finite() && body.mass > 0.0)
            || !body
                .center
                .into_iter()
                .chain(body.velocity)
                .chain(body.angular_velocity)
                .chain(body.inertia_world.into_iter().flatten())
                .chain(external_impulse)
                .chain(moment_origin)
                .all(f64::is_finite)
            || !low_mach(body.velocity)
        {
            return Err(invalid(
                "body state, impulse, origin or mass is inadmissible",
            ));
        }
        let inertia: [f64; 9] = core::array::from_fn(|i| body.inertia_world[i / 3][i % 3]);
        if (0..3).any(|i| (0..i).any(|j| inertia[i * 3 + j] != inertia[j * 3 + i])) {
            return Err(invalid("world inertia must be symmetric"));
        }
        fs_la::factor::cholesky(&inertia, 3)
            .map_err(|_| invalid("world inertia must be numerically positive definite"))?;

        let mut increment = external_impulse;
        let mut weight_sum = 0.0;
        let mut first = [0.0; 3];
        let mut rotary = [[0.0; 3]; 3];
        for tile in 0..self.f[0].len() {
            checkpoint(gate)?;
            for lane in 0..TILE_CELLS {
                let cell = self.coords(tile, lane);
                if self
                    .force_field
                    .as_ref()
                    .is_some_and(|f| f[tile * TILE_CELLS + lane] != [0.0; 3])
                {
                    return Err(PartialSaturationError3::Cell {
                        cell,
                        reason: "immersed coupling requires zero fluid forcing",
                    });
                }
                let epsilon = solid_fractions[(cell[2] * self.ny + cell[1]) * self.nx + cell[0]];
                if !(epsilon.is_finite() && (0.0..=1.0).contains(&epsilon))
                    || (epsilon > 0.0 && self.solid[tile] & (1u64 << lane) != 0)
                {
                    return Err(PartialSaturationError3::Cell {
                        cell,
                        reason: "fraction must be in [0,1] and positive only on fluid",
                    });
                }
                if epsilon == 0.0 {
                    continue;
                }
                let r = core::array::from_fn(|axis| cell[axis] as f64 + 0.5 - body.center[axis]);
                let solid_velocity = point_velocity(body.velocity, body.angular_velocity, r);
                if !low_mach(solid_velocity) {
                    return Err(PartialSaturationError3::Cell {
                        cell,
                        reason: "initial body point exceeds low-Mach envelope",
                    });
                }
                let (rho, velocity) = moments(self.f.each_ref().map(|f| f[tile].0[lane]), cell)?;
                let weight = epsilon * (tau - 0.5) / (1.0 - epsilon + tau - 0.5) * rho;
                let slip = core::array::from_fn(|axis| velocity[axis] - solid_velocity[axis]);
                let torque = cross(r, slip);
                weight_sum += weight;
                for axis in 0..3 {
                    first[axis] = weight.mul_add(r[axis], first[axis]);
                    increment[axis] = weight.mul_add(slip[axis], increment[axis]);
                    increment[axis + 3] = weight.mul_add(torque[axis], increment[axis + 3]);
                    let b = (axis + 1) % 3;
                    let c = (axis + 2) % 3;
                    rotary[axis][axis] =
                        weight.mul_add(r[b].mul_add(r[b], r[c] * r[c]), rotary[axis][axis]);
                    for other in 0..axis {
                        rotary[axis][other] =
                            (-weight * r[axis]).mul_add(r[other], rotary[axis][other]);
                    }
                }
            }
        }
        let mut matrix = [0.0; 36];
        for axis in 0..3 {
            matrix[axis * 6 + axis] = body.mass + 0.5 * weight_sum;
            for other in 0..=axis {
                matrix[(axis + 3) * 6 + other + 3] =
                    body.inertia_world[axis][other] + 0.5 * rotary[axis][other];
            }
        }
        // Lower-left block is sum(B rho [r]x)/2; upper-right is its transpose.
        matrix[3 * 6 + 1] = -0.5 * first[2];
        matrix[3 * 6 + 2] = 0.5 * first[1];
        matrix[4 * 6] = 0.5 * first[2];
        matrix[4 * 6 + 2] = -0.5 * first[0];
        matrix[5 * 6] = -0.5 * first[1];
        matrix[5 * 6 + 1] = 0.5 * first[0];
        for i in 0..6 {
            for j in 0..i {
                matrix[j * 6 + i] = matrix[i * 6 + j];
            }
        }
        if !matrix.into_iter().chain(increment).all(f64::is_finite) {
            return Err(invalid("spatial inertia/response overflowed"));
        }
        fs_la::factor::cholesky(&matrix, 6)
            .map_err(|_| invalid("coupled inertia is not numerically positive definite"))?
            .solve(&mut increment);
        let midpoint = core::array::from_fn(|i| body.velocity[i] + 0.5 * increment[i]);
        let after = core::array::from_fn(|i| body.velocity[i] + increment[i]);
        let omega_mid = core::array::from_fn(|i| body.angular_velocity[i] + 0.5 * increment[i + 3]);
        let omega_after = core::array::from_fn(|i| body.angular_velocity[i] + increment[i + 3]);
        if !low_mach(midpoint)
            || !low_mach(after)
            || !omega_mid.into_iter().chain(omega_after).all(f64::is_finite)
        {
            return Err(invalid(
                "solved body velocity exceeds finite low-Mach envelope",
            ));
        }

        let mut exchange = PartialSaturationExchange3::default();
        let mut walls = WallMomentumExchange3::default();
        for tile in 0..self.f[0].len() {
            checkpoint(gate)?;
            for lane in 0..TILE_CELLS {
                if self.solid[tile] & (1u64 << lane) != 0 {
                    continue;
                }
                let cell = self.coords(tile, lane);
                let f = self.f.each_ref().map(|f| f[tile].0[lane]);
                let ordinary = collide_cell3(f, self.collision_model, [0.0; 3])
                    .map_err(|source| PartialSaturationError3::Collision { cell, source })?;
                let mut post = ordinary;
                let epsilon = solid_fractions[(cell[2] * self.ny + cell[1]) * self.nx + cell[0]];
                if epsilon > 0.0 {
                    let r =
                        core::array::from_fn(|axis| cell[axis] as f64 + 0.5 - body.center[axis]);
                    let solid_velocity = point_velocity(midpoint, omega_mid, r);
                    if !low_mach(solid_velocity) || !low_mach(point_velocity(after, omega_after, r))
                    {
                        return Err(PartialSaturationError3::Cell {
                            cell,
                            reason: "solved body point exceeds low-Mach envelope",
                        });
                    }
                    let (rho, velocity) = moments(f, cell)?;
                    let fluid_eq = equilibrium3(rho, velocity);
                    let solid_eq = equilibrium3(rho, solid_velocity);
                    let b = epsilon * (tau - 0.5) / (1.0 - epsilon + tau - 0.5);
                    let mut fluid_impulse = [0.0; 3];
                    for q in 0..Q3 {
                        let fluid_collision = (fluid_eq[q] - f[q]) / tau;
                        let solid_collision = f[OPP3[q]] - f[q] + solid_eq[q] - fluid_eq[OPP3[q]];
                        post[q] = f[q] + (1.0 - b) * fluid_collision + b * solid_collision;
                        let correction = post[q] - ordinary[q];
                        exchange.fluid_population_mass_change += correction;
                        fluid_impulse[0] += f64::from(E3[q].0) * correction;
                        fluid_impulse[1] += f64::from(E3[q].1) * correction;
                        fluid_impulse[2] += f64::from(E3[q].2) * correction;
                    }
                    let offset =
                        core::array::from_fn(|axis| cell[axis] as f64 + 0.5 - moment_origin[axis]);
                    let angular = cross(offset, fluid_impulse);
                    let spin = cross(r, fluid_impulse);
                    for axis in 0..3 {
                        exchange.fluid_population_impulse[axis] += fluid_impulse[axis];
                        exchange.body_impulse[axis] -= fluid_impulse[axis];
                        exchange.fluid_population_angular_impulse[axis] += angular[axis];
                        exchange.body_angular_impulse[axis] -= angular[axis];
                        exchange.body_angular_impulse_about_center[axis] -= spin[axis];
                    }
                    exchange.body_work -= dot(fluid_impulse, solid_velocity);
                    exchange.coupled_cells += 1;
                    exchange.coupling_weight_sum += b;
                }
                if !post.into_iter().all(f64::is_finite) {
                    return Err(PartialSaturationError3::Cell {
                        cell,
                        reason: "blended collision overflowed",
                    });
                }
                for (field, value) in self.post.iter_mut().zip(post) {
                    field[tile].0[lane] = value;
                }
                let mask = self.link_masks[tile].0[lane];
                if mask != 0 {
                    let rho = post.iter().sum();
                    for q in 1..Q3 {
                        if mask & (1u32 << q) != 0 {
                            walls.add_link(
                                BoundaryLink3 {
                                    tile,
                                    lane,
                                    cell,
                                    direction: q,
                                },
                                post[OPP3[q]],
                                rho,
                                [0.0; 3],
                                moment_origin,
                            );
                        }
                    }
                }
            }
        }
        let walls = walls
            .checked()
            .map_err(|_| invalid("wall exchange overflowed"))?;
        let dv = core::array::from_fn(|i| after[i] - body.velocity[i]);
        let dw = core::array::from_fn(|i| omega_after[i] - body.angular_velocity[i]);
        let dp = dv.map(|v| body.mass * v);
        let dl = body.inertia_world.map(|row| dot(row, dw));
        let mean_v = core::array::from_fn(|i| 0.5 * (after[i] + body.velocity[i]));
        let mean_w = core::array::from_fn(|i| 0.5 * (omega_after[i] + body.angular_velocity[i]));
        let kinetic_energy_change = dot(dp, mean_v) + dot(dl, mean_w);
        let external_work = dot(external_impulse[..3].try_into().unwrap(), midpoint)
            + dot(external_impulse[3..].try_into().unwrap(), omega_mid);
        let momentum_residual = core::array::from_fn(|i| {
            if i < 3 {
                dp[i] - external_impulse[i] - exchange.body_impulse[i]
            } else {
                dl[i - 3] - external_impulse[i] - exchange.body_angular_impulse_about_center[i - 3]
            }
        });
        if ![
            exchange.fluid_population_mass_change,
            exchange.body_work,
            exchange.coupling_weight_sum,
            kinetic_energy_change,
            external_work,
        ]
        .into_iter()
        .chain(exchange.fluid_population_impulse)
        .chain(exchange.fluid_population_angular_impulse)
        .chain(exchange.body_angular_impulse_about_center)
        .chain(momentum_residual)
        .all(f64::is_finite)
        {
            return Err(invalid("body energy or momentum exchange overflowed"));
        }

        let mut next = match self.next.take() {
            Some(next) => next,
            None => {
                let mut next = core::array::from_fn(|_| Vec::new());
                for field in &mut next {
                    field
                        .try_reserve_exact(self.f[0].len())
                        .map_err(|_| invalid("streaming scratch allocation refused"))?;
                }
                // First-use scratch preserves solid slots too, with the same
                // tile-sized cancellation bound as collision and streaming.
                for tile in 0..self.f[0].len() {
                    checkpoint(gate)?;
                    for (q, field) in next.iter_mut().enumerate() {
                        field.push(self.f[q][tile].clone());
                    }
                }
                next
            }
        };
        let view = StreamView {
            post: &self.post,
            pull: &self.pull,
            link_masks: &self.link_masks,
            open_link_masks: &self.open_link_masks,
            stationary_link_masks: &self.stationary_link_masks,
            boundaries: self.boundaries,
            dims: [self.nx, self.ny, self.nz],
        };
        for tile in 0..self.f[0].len() {
            if let Err(error) = checkpoint(gate) {
                self.next = Some(next);
                return Err(error);
            }
            for lane in 0..TILE_CELLS {
                if self.solid[tile] & (1u64 << lane) == 0 {
                    let post = view.pull_cell(tile, lane);
                    for (field, value) in next.iter_mut().zip(post) {
                        field[tile].0[lane] = value;
                    }
                }
            }
        }
        if let Err(error) = checkpoint(gate) {
            self.next = Some(next);
            return Err(error);
        }
        let proposed_body = PartialSaturationBody3 {
            velocity: after,
            angular_velocity: omega_after,
            ..*body
        };
        let receipt = PartialSaturationBodyExchange3 {
            hydrodynamic: exchange,
            walls,
            velocity_midpoint: midpoint,
            angular_velocity_midpoint: omega_mid,
            kinetic_energy_change,
            external_work,
            momentum_residual,
        };
        let final_value = match finish(self, &proposed_body, &receipt) {
            Ok(value) => value,
            Err(error) => {
                self.next = Some(next);
                return Err(error);
            }
        };
        if let Err(error) = checkpoint(gate) {
            self.next = Some(next);
            return Err(error);
        }
        core::mem::swap(&mut self.f, &mut next);
        self.next = Some(next);
        self.topology_locked = true;
        *body = proposed_body;
        Ok((receipt, final_value))
    }
}
