//! Tangential exterior-wall motion on the original open-grid pull stencil.
//!
//! Wall and reservoir links retain disjoint ownership at their intersections.
//! The complete moving-wall correction must conserve mass at each fluid cell;
//! incompatible wall junctions refuse instead of averaging or clipping motion.

use super::*;
use crate::lattice3::moving_walls::{MovingWallMomentum3, WallExchange, reflect};

mod profiles;
pub use profiles::WallVelocityProfileReport3;

pub(super) struct WallVelocityProfile3 {
    values: Vec<[f64; 3]>,
    pending: Vec<[f64; 3]>,
    moving: bool,
}
pub(super) type WallProfileView3<'a> = [Option<&'a [[f64; 3]]>; 6];

pub(super) fn wall_profile_values(
    profiles: &[Option<WallVelocityProfile3>; 6],
) -> WallProfileView3<'_> {
    core::array::from_fn(|face| profiles[face].as_ref().map(|p| p.values.as_slice()))
}

fn plane_index(dims: [usize; 3], face: Face3, cell: usize) -> usize {
    let axis = face as usize / 2;
    let stride = dims[..axis].iter().product::<usize>();
    cell / (stride * dims[axis]) * stride + cell % stride
}

impl<L: EquilibriumLattice3<Q>, const Q: usize> OpenBoundaryGrid3<L, Q> {
    /// Actual exterior-wall impulse, angular impulse and work from the last
    /// accepted interval with at least one moving wall. Interior solid links
    /// retain their stationary law and are outside this exterior-wall account.
    ///
    /// The relative-velocity convention closes `wall + fluid = moving_mass`.
    /// Wall work is delivered to the wall; its negative is prescribed actuator
    /// work on the fluid, not a heat source or an isothermal kinetic-energy
    /// balance. Units and moment origin match the native moving-wall owner.
    /// No receipt exists before a moving step or after a wholly stationary step.
    /// Refused/cancelled steps and boundary rebindings retain the last receipt.
    #[must_use]
    pub fn last_wall_exchange(&self) -> Option<MovingWallMomentum3> {
        self.wall_exchange
    }

    // All coupled owners call this only at their final publication point.
    pub(super) fn publish_wall_exchange(&mut self) {
        self.wall_exchange = self.pending_wall_exchange;
    }

    pub(super) fn wall_face_moves(&self, face: Face3) -> bool {
        self.wall_profiles[face as usize].as_ref().map_or_else(
            || matches!(self.plan.boundaries.face(face), FaceBoundary3::Wall { velocity } if velocity != [0.0; 3]),
            |profile| profile.moving,
        )
    }
}

impl FacePlan {
    fn touches_wall(&self, p: [usize; 3]) -> bool {
        (0..3).any(|a| {
            (p[a] == 0
                && matches!(
                    self.boundaries.face(Face3::ALL[2 * a]),
                    FaceBoundary3::Wall { .. }
                ))
                || (p[a] + 1 == self.dims[a]
                    && matches!(
                        self.boundaries.face(Face3::ALL[2 * a + 1]),
                        FaceBoundary3::Wall { .. }
                    ))
        })
    }

    // Surface cells in canonical x-fast order, with no duplicate edges or
    // cell-sized map. Interior x rows contribute no iteration or allocation.
    fn wall_cells(self) -> impl Iterator<Item = usize> {
        let wall =
            Face3::ALL.map(|face| matches!(self.boundaries.face(face), FaceBoundary3::Wall { .. }));
        let [nx, ny, nz] = self.dims;
        (0..nz).flat_map(move |z| {
            (0..ny).flat_map(move |y| {
                let full = (y == 0 && wall[2])
                    || (y + 1 == ny && wall[3])
                    || (z == 0 && wall[4])
                    || (z + 1 == nz && wall[5]);
                let count = if full {
                    nx
                } else if nx == 1 {
                    usize::from(wall[0] || wall[1])
                } else {
                    usize::from(wall[0]) + usize::from(wall[1])
                };
                (0..count).map(move |i| {
                    let x = if full {
                        i
                    } else if wall[0] && i == 0 {
                        0
                    } else {
                        nx - 1
                    };
                    (z * ny + y) * nx + x
                })
            })
        })
    }

    pub(super) fn exterior_wall_velocity<E, const Q: usize>(
        &self,
        profiles: &WallProfileView3<'_>,
        p: [usize; 3],
        direction: [i32; 3],
        q: usize,
    ) -> Result<Option<[f64; 3]>, GridError3<E>> {
        let mut wall = None;
        for a in 0..3 {
            let from = p[a] as isize - direction[a] as isize;
            if from < 0 || from >= self.dims[a] as isize {
                let face = 2 * a + usize::from(from >= 0);
                if let FaceBoundary3::Wall { velocity: uniform } =
                    self.boundaries.face(Face3::ALL[face])
                {
                    let velocity = profiles[face].map_or(uniform, |values| {
                        let cell = (p[2] * self.dims[1] + p[1]) * self.dims[0] + p[0];
                        values[plane_index(self.dims, Face3::ALL[face], cell) * Q + q]
                    });
                    if wall.is_some_and(|previous| previous != velocity) {
                        return Err(GridError3::Input(
                            "intersecting fixed walls require compatible velocities",
                        ));
                    }
                    wall = Some(velocity);
                }
            }
        }
        Ok(wall)
    }

    pub(super) fn admit_wall_motion<L: Lattice3<Q>, const Q: usize>(
        &self,
        grid: &Grid3<L, Q>,
        profiles: WallProfileView3<'_>,
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<(), GridError3<L::Error>> {
        if !self.moving {
            return Ok(());
        }
        for (index, cell) in self.wall_cells().enumerate() {
            if index % 64 == 0 {
                poll()?;
            }
            if grid.solid[cell] {
                continue;
            }
            let p = coordinates(cell, self.dims);
            let mut mass = 0.0;
            let mut magnitude = 0.0;
            for q in 1..Q {
                if let Some(velocity) =
                    self.exterior_wall_velocity::<L::Error, Q>(&profiles, p, L::E[q], q)?
                {
                    // Every correction at this cell shares its post-collision
                    // density, so the impermeability check is state independent.
                    let correction = reflect::<L, Q>(0.0, 1.0, q, velocity);
                    mass += correction;
                    magnitude += correction.abs();
                }
            }
            if !mass.is_finite()
                || !magnitude.is_finite()
                || mass.abs() > 64.0 * f64::EPSILON * magnitude
            {
                return Err(GridError3::Input(
                    "fixed moving-wall links must conserve local mass",
                ));
            }
        }
        poll()
    }

    pub(super) fn gradient_wall_velocity<L: Lattice3<Q>, const Q: usize>(
        &self,
        profiles: &[Option<WallVelocityProfile3>; 6],
        cell: usize,
        axis: usize,
        upper: bool,
    ) -> [f64; 3] {
        let p = coordinates(cell, self.dims);
        if (upper && p[axis] + 1 == self.dims[axis]) || (!upper && p[axis] == 0) {
            if let FaceBoundary3::Wall { velocity } = self
                .boundaries
                .face(Face3::ALL[2 * axis + usize::from(upper)])
            {
                if let Some(profile) = &profiles[2 * axis + usize::from(upper)] {
                    let mut normal = [0; 3];
                    normal[axis] = if upper { -1 } else { 1 };
                    let q = L::E
                        .iter()
                        .position(|&e| e == normal)
                        .expect("native stencil has straight normal links");
                    return profile.values[plane_index(
                        self.dims,
                        Face3::ALL[2 * axis + usize::from(upper)],
                        cell,
                    ) * Q
                        + q];
                }
                return velocity;
            }
        }
        // Interior voxels and the provisional open-normal stencil are static.
        // The open rheology owner replaces that latter derivative separately.
        [0.0; 3]
    }
}

impl<L: EquilibriumLattice3<Q>, const Q: usize> OpenStep<'_, L, Q> {
    pub(super) fn reflect_walls(
        &self,
        src: &[[f64; Q]],
        output: &mut [f64; Q],
        cell: usize,
    ) -> Result<(), GridError3<L::Error>> {
        if !self.plan.moving {
            return Ok(());
        }
        let p = coordinates(cell, self.data.dims);
        if !self.plan.touches_wall(p) {
            return Ok(());
        }
        let rho = src[cell].iter().sum();
        for (q, value) in output.iter_mut().enumerate().skip(1) {
            if let Some(velocity) = self.plan.exterior_wall_velocity::<L::Error, Q>(
                &self.wall_profiles,
                p,
                L::E[q],
                q,
            )? && velocity != [0.0; 3]
            {
                *value = reflect::<L, Q>(src[cell][L::OPP[q]], rho, q, velocity);
            }
        }
        Ok(())
    }

    pub(super) fn measure_walls(
        &self,
        post: &[[f64; Q]],
        next: &[[f64; Q]],
        origin: [f64; 3],
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<Option<MovingWallMomentum3>, GridError3<L::Error>> {
        if !self.plan.moving {
            return Ok(None);
        }
        let mut exchange = WallExchange::default();
        for (index, cell) in self.plan.wall_cells().enumerate() {
            if index % 64 == 0 {
                poll()?;
            }
            if self.data.solid[cell] {
                continue;
            }
            let p = coordinates(cell, self.data.dims);
            for q in 1..Q {
                if let Some(velocity) = self.plan.exterior_wall_velocity::<L::Error, Q>(
                    &self.wall_profiles,
                    p,
                    L::E[q],
                    q,
                )? {
                    let e = L::E[q].map(f64::from);
                    let arm = core::array::from_fn(|a| p[a] as f64 + 0.5 - 0.5 * e[a] - origin[a]);
                    exchange.add(cell, post[cell][L::OPP[q]], next[cell][q], e, velocity, arm)?;
                }
            }
        }
        exchange.finish().map(Some)
    }
}
