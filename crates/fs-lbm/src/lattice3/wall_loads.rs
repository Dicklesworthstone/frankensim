//! Native momentum-exchange loads on stationary halfway walls.
//!
//! A probe compiles its selected links once. Each measured step reads only
//! those links from the collision and stream candidates, before publication;
//! it neither reconstructs a stress field nor runs a second collision. The
//! returned impulse is ON THE WALL, including pressure and viscous traction.
//! For lattice spacing/time step one it is also the mean force over the step.
//! Physical force scales by rho_scale*dx^4/dt^2 and torque by rho_scale*dx^5/dt^2,
//! where rho_scale converts lattice density to kg/m^3.

use super::{Grid3, GridError3, Lattice3, coordinates};
use core::marker::PhantomData;
use fs_exec::{CancelGate, Cx, KernelRunner};

/// A single native reflected link, offered to the probe's selection predicate.
/// Each link is offered exactly once, including edges/corners of exterior walls.
#[derive(Clone, Copy, Debug)]
pub struct StationaryWallLink3 {
    pub fluid_cell: [usize; 3],
    /// The population direction pointing FROM the wall INTO the fluid cell.
    pub incoming_direction: [i32; 3],
    /// Wrapped solid source cell, or None for a nonperiodic domain boundary.
    pub solid_cell: Option<[usize; 3]>,
    /// Crossed nonperiodic faces in x-, x+, y-, y+, z-, z+ order.
    /// A diagonal corner link can cross several faces but is never duplicated.
    pub exterior_faces: [bool; 6],
    /// Halfway link intersection; cell centres are at integer coordinates + 1/2.
    pub position: [f64; 3],
}

struct MeasuredLink {
    cell: usize,
    incoming: usize,
    solid: Option<usize>,
    arm: [f64; 3],
}

/// Immutable sparse selection on a grid's fixed stationary wall topology.
/// Selections may be reused for later steps or a grid with compatible selected
/// links. Dimensions, periodicity and each selected fluid/solid link are checked.
pub struct WallMomentumProbe3<L: Lattice3<Q>, const Q: usize> {
    dims: [usize; 3],
    periodic: [bool; 3],
    links: Vec<MeasuredLink>,
    marker: PhantomData<L>,
}

/// Actual reflected-population impulse and its moment about the probe origin.
/// Stationary walls exchange momentum while doing zero mechanical work.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct StationaryWallMomentum3 {
    pub wall_impulse: [f64; 3],
    pub wall_angular_impulse: [f64; 3],
    pub measured_links: usize,
}

impl<L: Lattice3<Q>, const Q: usize> WallMomentumProbe3<L, Q> {
    #[must_use]
    pub fn measured_links(&self) -> usize {
        self.links.len()
    }

    fn validate(
        &self,
        grid: &Grid3<L, Q>,
        mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<(), GridError3<L::Error>> {
        checkpoint()?;
        if self.dims != grid.dims || self.periodic != grid.periodic {
            return Err(GridError3::Input("wall probe grid topology"));
        }
        for (index, link) in self.links.iter().enumerate() {
            if index % 64 == 0 {
                checkpoint()?;
            }
            if grid.solid[link.cell] || link.solid.is_some_and(|cell| !grid.solid[cell]) {
                return Err(GridError3::Input("wall probe selected link topology"));
            }
        }
        checkpoint()
    }

    fn measure(
        &self,
        grid: &Grid3<L, Q>,
        mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<StationaryWallMomentum3, GridError3<L::Error>> {
        let mut sum = [0.0; 6];
        let mut correction = [0.0; 6];
        for (index, link) in self.links.iter().enumerate() {
            if index % 64 == 0 {
                checkpoint()?;
            }
            let q = link.incoming;
            let outgoing = grid.post[link.cell][L::OPP[q]];
            let incoming = grid.next[link.cell][q];
            let impulse = L::E[q].map(|e| -(outgoing + incoming) * f64::from(e));
            let r = link.arm;
            let torque = [
                r[1] * impulse[2] - r[2] * impulse[1],
                r[2] * impulse[0] - r[0] * impulse[2],
                r[0] * impulse[1] - r[1] * impulse[0],
            ];
            // Compensated canonical link-order sums resolve small net drag
            // against much larger cancelling equilibrium pressure loads.
            for (axis, value) in impulse.into_iter().chain(torque).enumerate() {
                let next = sum[axis] + value;
                correction[axis] += if sum[axis].abs() >= value.abs() {
                    (sum[axis] - next) + value
                } else {
                    (value - next) + sum[axis]
                };
                sum[axis] = next;
                if !next.is_finite() || !correction[axis].is_finite() {
                    return Err(GridError3::NonFiniteCell { cell: link.cell });
                }
            }
        }
        for axis in 0..6 {
            sum[axis] += correction[axis];
            if !sum[axis].is_finite() {
                return Err(GridError3::Input("nonfinite total wall momentum"));
            }
        }
        checkpoint()?;
        Ok(StationaryWallMomentum3 {
            wall_impulse: [sum[0], sum[1], sum[2]],
            wall_angular_impulse: [sum[3], sum[4], sum[5]],
            measured_links: self.links.len(),
        })
    }
}

impl<L: Lattice3<Q>, const Q: usize> Grid3<L, Q> {
    /// Compile a sparse stationary wall selection in x-fastest cell, then
    /// lattice direction order. Select all solid links for an immersed body's
    /// total drag, or specific exterior faces for channel pressure/shear loads.
    /// The origin and wall positions use the domain's cell-coordinate image.
    /// A body crossing a periodic seam needs a continuous-image torque assembly;
    /// this probe does not implicitly wrap its moment arms around the origin.
    pub fn wall_momentum_probe(
        &self,
        cx: &Cx<'_>,
        origin: [f64; 3],
        mut select: impl FnMut(StationaryWallLink3) -> bool,
    ) -> Result<WallMomentumProbe3<L, Q>, GridError3<L::Error>> {
        if !origin.iter().all(|v| v.is_finite()) {
            return Err(GridError3::Input("wall moment origin"));
        }
        let mut links = Vec::new();
        for cell in 0..self.f.len() {
            if cell % 64 == 0 {
                cx.checkpoint().map_err(|_| GridError3::Cancelled)?;
            }
            if self.solid[cell] {
                continue;
            }
            let p = coordinates(cell, self.dims);
            for q in 0..Q {
                let mut from = [0; 3];
                let mut exterior_faces = [false; 6];
                for a in 0..3 {
                    let v = p[a] as isize - L::E[q][a] as isize;
                    if (v < 0 || v >= self.dims[a] as isize) && !self.periodic[a] {
                        exterior_faces[2 * a + usize::from(v >= 0)] = true;
                    }
                    from[a] = v.rem_euclid(self.dims[a] as isize) as usize;
                }
                let exterior = exterior_faces.iter().any(|&face| face);
                let source = (from[2] * self.dims[1] + from[1]) * self.dims[0] + from[0];
                if !exterior && !self.solid[source] {
                    continue;
                }
                let position =
                    core::array::from_fn(|a| p[a] as f64 + 0.5 - 0.5 * f64::from(L::E[q][a]));
                if select(StationaryWallLink3 {
                    fluid_cell: p,
                    incoming_direction: L::E[q],
                    solid_cell: (!exterior).then_some(from),
                    exterior_faces,
                    position,
                }) {
                    links.try_reserve(1).map_err(|_| GridError3::Allocation)?;
                    links.push(MeasuredLink {
                        cell,
                        incoming: q,
                        solid: (!exterior).then_some(source),
                        arm: core::array::from_fn(|a| position[a] - origin[a]),
                    });
                }
            }
        }
        cx.checkpoint().map_err(|_| GridError3::Cancelled)?;
        Ok(WallMomentumProbe3 {
            dims: self.dims,
            periodic: self.periodic,
            links,
            marker: PhantomData,
        })
    }

    /// Advance the original serial solver and return actual wall loads only
    /// when the step is accepted. An optional collision field binds the same
    /// material/SGS viscosity as `step_serial_with_collisions`; None uses the
    /// constructor law. Measurement failure/cancellation publishes neither.
    pub fn step_serial_with_wall_momentum(
        &mut self,
        cx: &Cx<'_>,
        probe: &WallMomentumProbe3<L, Q>,
        collisions: Option<&[L::Collision]>,
    ) -> Result<StationaryWallMomentum3, GridError3<L::Error>> {
        probe.validate(self, || cx.checkpoint().map_err(|_| GridError3::Cancelled))?;
        if let Some(field) = collisions {
            self.admit_collisions(field, || cx.checkpoint().map_err(|_| GridError3::Cancelled))?;
        }
        self.step_serial_observed(cx, collisions, |grid| {
            probe.measure(grid, || cx.checkpoint().map_err(|_| GridError3::Cancelled))
        })
    }

    /// Pooled counterpart. Collision and streaming use the original tiles;
    /// sparse compensated load reduction has the serial canonical order,
    /// independently of worker count. Ordinary steps allocate no probe storage.
    pub fn step_pooled_with_wall_momentum<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        probe: &WallMomentumProbe3<L, Q>,
        collisions: Option<&[L::Collision]>,
    ) -> Result<StationaryWallMomentum3, GridError3<L::Error>> {
        let poll = || {
            if gate.is_requested() {
                Err(GridError3::Cancelled)
            } else {
                Ok(())
            }
        };
        probe.validate(self, poll)?;
        if let Some(field) = collisions {
            self.admit_collisions(field, poll)?;
        }
        self.step_pooled_observed(pool, gate, collisions, |grid| probe.measure(grid, poll))
    }
}
