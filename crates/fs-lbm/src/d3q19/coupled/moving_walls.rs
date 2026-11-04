//! Prescribed tangential plate motion on the original D3Q19 pull stream.
//!
//! The ordinary moving halfway reflection is
//! f_in = f_out + 2 w_q rho (e_q.u_wall)/cs², with post-collision rho.
//! Relative-velocity momentum exchange follows the same Wen et al. convention
//! as `WallMomentumExchange3` (<https://arxiv.org/abs/1303.0625>).
//! Geometry and represented volume remain fixed; neither normal motion nor
//! swept-cell remapping is represented. Mechanical wall work is not a heat source.

use super::{CS2_F, E3, OPP3, PlateStepError3, PlatesGrid3, Q3, admit_force_frame3};
use crate::d3q19::{W3, WallMomentumExchange3};

#[derive(Clone, Copy, Debug)]
pub(super) struct PlateMotion3 {
    velocity: [[f64; 3]; 2],
    origin: [f64; 3],
}

fn cross(r: [f64; 3], p: [f64; 3]) -> [f64; 3] {
    [
        r[1] * p[2] - r[2] * p[1],
        r[2] * p[0] - r[0] * p[2],
        r[0] * p[1] - r[1] * p[0],
    ]
}

#[derive(Clone, Copy, Default)]
struct Exchange {
    sum: [f64; 20],
    correction: [f64; 20],
    links: usize,
}

impl Exchange {
    fn add(
        &mut self,
        outgoing: f64,
        incoming: f64,
        e: [f64; 3],
        velocity: [f64; 3],
        arm: [f64; 3],
    ) -> Result<(), PlateStepError3> {
        let dm = incoming - outgoing;
        let fluid = e.map(|v| v * (outgoing + incoming));
        let moving_mass = velocity.map(|v| v * dm);
        let wall = if velocity == [0.0; 3] {
            e.map(|v| -(outgoing + incoming) * v)
        } else {
            core::array::from_fn(|a| {
                (-e[a] - velocity[a]) * outgoing - (e[a] - velocity[a]) * incoming
            })
        };
        let work = wall[0] * velocity[0] + wall[1] * velocity[1] + wall[2] * velocity[2];
        let vectors = [
            wall,
            fluid,
            moving_mass,
            cross(arm, wall),
            cross(arm, fluid),
            cross(arm, moving_mass),
        ];
        for (a, value) in vectors.into_iter().flatten().chain([dm, work]).enumerate() {
            let next = self.sum[a] + value;
            self.correction[a] += if self.sum[a].abs() >= value.abs() {
                (self.sum[a] - next) + value
            } else {
                (value - next) + self.sum[a]
            };
            if !next.is_finite() || !self.correction[a].is_finite() {
                return Err(PlateStepError3::WallMotion(
                    "nonfinite plate momentum/work account",
                ));
            }
            self.sum[a] = next;
        }
        self.links += 1;
        Ok(())
    }

    fn finish(self) -> Result<WallMomentumExchange3, PlateStepError3> {
        let value: [f64; 20] = core::array::from_fn(|a| self.sum[a] + self.correction[a]);
        if value.iter().any(|v| !v.is_finite()) {
            return Err(PlateStepError3::WallMotion(
                "nonfinite summed plate exchange",
            ));
        }
        Ok(WallMomentumExchange3 {
            wall_impulse: [value[0], value[1], value[2]],
            fluid_population_impulse: [value[3], value[4], value[5]],
            wall_velocity_mass_impulse: [value[6], value[7], value[8]],
            wall_angular_impulse: [value[9], value[10], value[11]],
            fluid_population_angular_impulse: [value[12], value[13], value[14]],
            wall_velocity_mass_angular_impulse: [value[15], value[16], value[17]],
            fluid_mass_change: value[18],
            wall_work: value[19],
            measured_links: self.links,
        })
    }
}

impl PlatesGrid3 {
    /// Retain bottom/top uniform tangential velocities in lattice units.
    /// The fixed plates are z=1/2 and z=nz-3/2 in cell-index coordinates.
    /// `moment_origin` uses the same coordinates and defines both wall torques.
    /// With periodic tangents, the receipt uses this unwrapped grid image;
    /// periodic seam transport is a separate angular-momentum account.
    ///
    /// Velocity must be finite, normal velocity exactly zero and speed within
    /// the existing lattice Mach envelope. All admission precedes mutation.
    /// Changing a drive changes no population, viscosity, force or thermal
    /// history. The new drive applies from the next actual native step.
    /// A successful change invalidates the previous plate receipt.
    pub fn set_plate_motion(
        &mut self,
        velocities: [[f64; 3]; 2],
        moment_origin: [f64; 3],
    ) -> Result<(), PlateStepError3> {
        self.checked_layout()?;
        if moment_origin.iter().any(|v| !v.is_finite()) {
            return Err(PlateStepError3::WallMotion(
                "plate moment origin must be finite",
            ));
        }
        for velocity in velocities {
            let mach = (velocity.iter().map(|v| v * v).sum::<f64>() / CS2_F).sqrt();
            if velocity.iter().any(|v| !v.is_finite()) || velocity[2] != 0.0 {
                return Err(PlateStepError3::WallMotion(
                    "plate velocity must be finite and tangential",
                ));
            }
            if !mach.is_finite() || mach > crate::MACH_LIMIT {
                return Err(PlateStepError3::WallMotion(
                    "plate speed exceeds the lattice Mach envelope",
                ));
            }
        }
        self.plate_motion = Some(PlateMotion3 {
            velocity: velocities,
            origin: moment_origin,
        });
        self.plate_exchange = None;
        Ok(())
    }

    /// Restore stationary reflection without resetting flow or thermal history.
    /// Use an explicitly zero drive to continue measuring stationary wall loads.
    pub fn reset_plate_motion(&mut self) {
        self.plate_motion = None;
        self.plate_exchange = None;
    }

    /// Bottom/top prescribed lattice velocities; absent drive means rest.
    #[must_use]
    pub fn plate_velocities(&self) -> [[f64; 3]; 2] {
        self.plate_motion
            .map_or([[0.0; 3]; 2], |motion| motion.velocity)
    }

    /// Actual native-stream wall impulse, angular impulse and work, in lattice
    /// units, bottom then top. `wall + fluid_population = wall_velocity_mass`.
    /// Negative wall work is energy supplied by the prescribed wall actuator.
    /// This is not an internal-energy increment or a transient dissipation audit.
    /// Refused native steps retain the last accepted receipt. No receipt exists
    /// before a driven step or after changing/resetting the drive.
    #[must_use]
    pub const fn last_plate_exchange(&self) -> Option<[WallMomentumExchange3; 2]> {
        self.plate_exchange
    }

    // Preserve the physical drive and accepted receipt in the original source
    // owner's full candidate snapshot, without exposing internal drive state.
    pub(crate) fn copy_plate_state_from(&mut self, original: &Self) {
        self.plate_motion = original.plate_motion;
        self.plate_exchange = original.plate_exchange;
    }

    #[inline]
    pub(super) fn reflect_plate(&self, outgoing: f64, rho: f64, q: usize, side: usize) -> f64 {
        let Some(motion) = self.plate_motion else {
            return outgoing;
        };
        let u = motion.velocity[side];
        if u == [0.0; 3] {
            return outgoing;
        }
        let e = E3[q];
        let eu = f64::from(e.0).mul_add(u[0], f64::from(e.1).mul_add(u[1], f64::from(e.2) * u[2]));
        outgoing + 2.0 * W3[q] * rho * eu / CS2_F
    }

    // Boundary-only preflight. The same reflection function is used by the
    // publishing stream. Stack storage suffices; interior work is unchanged.
    pub(super) fn prepare_plate_exchange(
        &self,
    ) -> Result<Option<[WallMomentumExchange3; 2]>, PlateStepError3> {
        let Some(motion) = self.plate_motion else {
            return Ok(None);
        };
        let mut exchange = [Exchange::default(); 2];
        for z in [1, self.nz - 2] {
            for y in 0..self.ny {
                for x in 0..self.nx {
                    let cell = self.idx(x, y, z);
                    let rho = self.post[cell].iter().sum::<f64>();
                    if !rho.is_finite() || rho <= 0.0 {
                        return Err(PlateStepError3::WallState { cell });
                    }
                    let mut candidate = [0.0; Q3];
                    for (q, incoming) in candidate.iter_mut().enumerate() {
                        let (ex, ey, ez) = E3[q];
                        let sx = (x as i64 - i64::from(ex)).rem_euclid(self.nx as i64) as usize;
                        let sy = (y as i64 - i64::from(ey)).rem_euclid(self.ny as i64) as usize;
                        let sz = (z as i64 - i64::from(ez)) as usize;
                        if self.is_wall_layer(sz) {
                            let side = usize::from(sz != 0);
                            let outgoing = self.post[cell][OPP3[q]];
                            *incoming = self.reflect_plate(outgoing, rho, q, side);
                            let e = [f64::from(ex), f64::from(ey), f64::from(ez)];
                            let p = [x as f64, y as f64, z as f64];
                            let arm =
                                core::array::from_fn(|a| p[a] - 0.5 * e[a] - motion.origin[a]);
                            exchange[side].add(
                                outgoing,
                                *incoming,
                                e,
                                motion.velocity[side],
                                arm,
                            )?;
                        } else {
                            *incoming = self.post[self.idx(sx, sy, sz)][q];
                        }
                    }
                    let velocity = admit_force_frame3(&candidate, self.fext[cell], cell)
                        .map_err(|_| PlateStepError3::WallState { cell })?;
                    let mach = (velocity.iter().map(|v| v * v).sum::<f64>() / CS2_F).sqrt();
                    if !mach.is_finite() || mach > crate::MACH_LIMIT {
                        return Err(PlateStepError3::WallState { cell });
                    }
                }
            }
        }
        Ok(Some([exchange[0].finish()?, exchange[1].finish()?]))
    }
}
