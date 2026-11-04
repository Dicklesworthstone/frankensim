//! Prescribed moving surfaces on the original, fixed halfway-wall topology.
//!
//! The native pull stream uses `f_in = f_out + 2 w_q rho (e_q.u_wall)/cs²`,
//! with post-collision density and `e_q` directed into the fluid. This is the
//! usual low-Mach, second-order moving bounce-back law, shared in convention
//! with the specialized D3Q19 wall owner. It does not change collision physics.
//! Surface velocity is held over one lattice step and is supplied in probe
//! order, so an evolving drive needs no topology recompilation or step allocation.
//!
//! Geometry remains fixed. Exterior motion must be tangential, and the full
//! selected link set at each fluid cell must have zero moving-wall mass source
//! to roundoff. Incomplete moving faces and normal/swept motion refuse; no
//! velocity is projected, clipped or silently replaced at corners. Interior
//! solid links obey the same discrete impermeability condition. This is not a
//! fresh-cell/remap, moving body, curved-wall or thermal-energy solver.
//!
//! The load uses relative incident/reflected velocities, following Wen et al.,
//! <https://arxiv.org/abs/1303.0625>. Separate population and moving-mass terms
//! retain the identity `wall + fluid_population = wall_velocity_mass` for both
//! linear and angular impulse. Wall work is sum(J_wall.u_wall); a prescribed
//! actuator supplies its negative to the fluid. This mechanical boundary work
//! is not the kinetic-energy change of a dissipative isothermal fluid.

use super::wall_loads::WallMomentumProbe3;
use super::{Grid3, GridError3, Lattice3, StepData, macro_cell};
use fs_exec::{CancelGate, Cx, KernelRunner};

/// The same physical receipt used by the specialized D3Q19 wall owner, now
/// supplied by shared D3Q19/D3Q27 transport. All values use lattice units.
/// Impulse scales by rho_scale*dx^4/dt; angular impulse by rho_scale*dx^5/dt;
/// work by rho_scale*dx^5/dt². Divide the impulses by dt for mean force/torque.
pub type MovingWallMomentum3 = crate::d3q19::WallMomentumExchange3;

pub(super) struct WallMotion3<'a, L: Lattice3<Q>, const Q: usize> {
    probe: &'a WallMomentumProbe3<L, Q>,
    velocities: &'a [[f64; 3]],
}

fn dot_direction(e: [i32; 3], u: [f64; 3]) -> f64 {
    f64::from(e[0]).mul_add(u[0], f64::from(e[1]).mul_add(u[1], f64::from(e[2]) * u[2]))
}

fn cross(r: [f64; 3], p: [f64; 3]) -> [f64; 3] {
    [
        r[1] * p[2] - r[2] * p[1],
        r[2] * p[0] - r[0] * p[2],
        r[0] * p[1] - r[1] * p[0],
    ]
}

impl<L: Lattice3<Q>, const Q: usize> WallMotion3<'_, L, Q> {
    fn admit(
        &self,
        grid: &Grid3<L, Q>,
        mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<bool, GridError3<L::Error>> {
        self.probe.validate(grid, &mut checkpoint)?;
        if self.velocities.len() != self.probe.links.len() {
            return Err(GridError3::Input("moving-wall velocity field length"));
        }
        let mut flux = 0.0;
        let mut magnitude = 0.0;
        let mut moving = false;
        for (i, (link, descriptor)) in self
            .probe
            .links
            .iter()
            .zip(self.probe.wall_links())
            .enumerate()
        {
            if i % 64 == 0 {
                checkpoint()?;
            }
            let u = self.velocities[i];
            moving |= u != [0.0; 3];
            if u.iter().any(|v| !v.is_finite()) {
                return Err(GridError3::Input("nonfinite moving-wall velocity"));
            }
            let mach = (u.iter().map(|v| v * v).sum::<f64>() / L::CS2).sqrt();
            if !mach.is_finite() || mach > grid.budget.maximum_mach {
                return Err(GridError3::Mach {
                    cell: link.cell,
                    observed: mach,
                    maximum: grid.budget.maximum_mach,
                });
            }
            for (face, &crossed) in descriptor.exterior_faces.iter().enumerate() {
                if crossed && u[face / 2] != 0.0 {
                    return Err(GridError3::Input(
                        "fixed exterior wall velocity must be tangential",
                    ));
                }
            }
            // All corrections at this cell share its one post-collision rho.
            // Admission is independent of density and permits no mass source.
            let source = 2.0 * L::W[link.incoming] * dot_direction(L::E[link.incoming], u) / L::CS2;
            flux += source;
            magnitude += source.abs();
            if i + 1 == self.probe.links.len() || self.probe.links[i + 1].cell != link.cell {
                if !flux.is_finite()
                    || !magnitude.is_finite()
                    || flux.abs() > 64.0 * f64::EPSILON * magnitude
                {
                    return Err(GridError3::Input(
                        "fixed moving-wall links must conserve local mass",
                    ));
                }
                flux = 0.0;
                magnitude = 0.0;
            }
        }
        checkpoint()?;
        Ok(moving)
    }

    fn measure(
        &self,
        grid: &Grid3<L, Q>,
        mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<MovingWallMomentum3, GridError3<L::Error>> {
        let mut sum = [0.0_f64; 20];
        let mut correction = [0.0_f64; 20];
        for (i, link) in self.probe.links.iter().enumerate() {
            if i % 64 == 0 {
                checkpoint()?;
            }
            let q = link.incoming;
            let outgoing = grid.post[link.cell][L::OPP[q]];
            let incoming = grid.next[link.cell][q];
            let velocity = self.velocities[i];
            let e = L::E[q].map(f64::from);
            let dm = incoming - outgoing;
            let fluid = e.map(|v| v * (outgoing + incoming));
            let moving_mass = velocity.map(|v| v * dm);
            let wall = if velocity == [0.0; 3] {
                // Preserve the original stationary report arithmetic exactly.
                e.map(|v| -(outgoing + incoming) * v)
            } else {
                core::array::from_fn(|a| {
                    (-e[a] - velocity[a]) * outgoing - (e[a] - velocity[a]) * incoming
                })
            };
            let work = wall[0] * velocity[0] + wall[1] * velocity[1] + wall[2] * velocity[2];
            let values = [
                wall,
                fluid,
                moving_mass,
                cross(link.arm, wall),
                cross(link.arm, fluid),
                cross(link.arm, moving_mass),
            ];
            for (a, value) in values.into_iter().flatten().chain([dm, work]).enumerate() {
                let next = sum[a] + value;
                correction[a] += if sum[a].abs() >= value.abs() {
                    (sum[a] - next) + value
                } else {
                    (value - next) + sum[a]
                };
                sum[a] = next;
                if !next.is_finite() || !correction[a].is_finite() {
                    return Err(GridError3::NonFiniteCell { cell: link.cell });
                }
            }
        }
        for a in 0..20 {
            sum[a] += correction[a];
            if !sum[a].is_finite() {
                return Err(GridError3::Input("nonfinite moving-wall receipt"));
            }
        }
        checkpoint()?;
        Ok(MovingWallMomentum3 {
            wall_impulse: [sum[0], sum[1], sum[2]],
            fluid_population_impulse: [sum[3], sum[4], sum[5]],
            wall_velocity_mass_impulse: [sum[6], sum[7], sum[8]],
            wall_angular_impulse: [sum[9], sum[10], sum[11]],
            fluid_population_angular_impulse: [sum[12], sum[13], sum[14]],
            wall_velocity_mass_angular_impulse: [sum[15], sum[16], sum[17]],
            fluid_mass_change: sum[18],
            wall_work: sum[19],
            measured_links: self.probe.links.len(),
        })
    }
}

impl<L: Lattice3<Q>, const Q: usize> StepData<'_, L, Q> {
    pub(super) fn process_moving_walls(
        &self,
        src: &[[f64; Q]],
        dst: &mut [[f64; Q]],
        base: usize,
        motion: &WallMotion3<'_, L, Q>,
    ) -> Result<(), GridError3<L::Error>> {
        let links = &motion.probe.links;
        let mut selected = links.partition_point(|link| link.cell < base);
        for (offset, output) in dst.iter_mut().enumerate() {
            let cell = base + offset;
            if self.solid[cell] {
                *output = src[cell];
                continue;
            }
            self.pull_cell(src, output, cell);
            if selected < links.len() && links[selected].cell == cell {
                let rho: f64 = src[cell].iter().sum();
                while selected < links.len() && links[selected].cell == cell {
                    let q = links[selected].incoming;
                    let u = motion.velocities[selected];
                    if u != [0.0; 3] {
                        // This is the actual streamed population, also read by
                        // the post-pass receipt; there is no second wall model.
                        output[q] = src[cell][L::OPP[q]]
                            + 2.0 * L::W[q] * rho * dot_direction(L::E[q], u) / L::CS2;
                    }
                    selected += 1;
                }
            }
            let force = self.forces.map_or(self.force, |field| field[cell]);
            macro_cell::<L, Q>(output, force, cell, self.budget.maximum_mach)?;
        }
        Ok(())
    }
}

impl<L: Lattice3<Q>, const Q: usize> Grid3<L, Q> {
    /// Advance with prescribed velocities on the probe's fixed wall links.
    /// `velocities` has exactly one vector per `probe.wall_links()` item. Other
    /// reflected links remain stationary. Include zero velocities to measure
    /// stationary surfaces in the same receipt. Surface speed, tangency and
    /// discrete local impermeability are admitted before collision.
    ///
    /// The optional material/SGS collision field has the ordinary step's
    /// meaning. Velocities and coefficients are borrowed for this step only;
    /// there is no latched drive, extra population buffer or accepted partial
    /// update. The final corrected stream and complete load receipt must both
    /// pass before populations and step count publish together.
    pub fn step_serial_with_moving_wall_momentum(
        &mut self,
        cx: &Cx<'_>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        collisions: Option<&[L::Collision]>,
    ) -> Result<MovingWallMomentum3, GridError3<L::Error>> {
        let poll = || cx.checkpoint().map_err(|_| GridError3::Cancelled);
        let motion = WallMotion3 { probe, velocities };
        let moving = motion.admit(self, poll)?;
        if let Some(field) = collisions {
            self.admit_collisions(field, poll)?;
        }
        self.step_serial_with_wall_motion(cx, collisions, moving.then_some(&motion), |grid| {
            motion.measure(grid, poll)
        })
    }

    /// Executor counterpart, with the original collision and pull-stream
    /// tiles. Sparse corrections are owned by their destination tile; the
    /// compensated receipt uses canonical probe order for worker invariance.
    pub fn step_pooled_with_moving_wall_momentum<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        collisions: Option<&[L::Collision]>,
    ) -> Result<MovingWallMomentum3, GridError3<L::Error>> {
        let poll = || {
            if gate.is_requested() {
                Err(GridError3::Cancelled)
            } else {
                Ok(())
            }
        };
        let motion = WallMotion3 { probe, velocities };
        let moving = motion.admit(self, poll)?;
        if let Some(field) = collisions {
            self.admit_collisions(field, poll)?;
        }
        self.step_pooled_with_wall_motion(
            pool,
            gate,
            collisions,
            moving.then_some(&motion),
            |grid| motion.measure(grid, poll),
        )
    }
}
