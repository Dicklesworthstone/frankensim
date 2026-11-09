//! Spatial velocity at actual halfway-wall links, held over each native tick.

use super::*;
use crate::lattice3::open_boundaries::forcing::ForcedOpenBoundaryGrid3;

/// Prescribed surface diagnostics; actual impulse and work come from the step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallVelocityProfileReport3 {
    pub face: Face3,
    pub sampled_links: usize,
    pub maximum_speed_lu: f64,
}

fn check<E>(gate: &CancelGate) -> Result<(), GridError3<E>> {
    if gate.is_requested() {
        Err(GridError3::Cancelled)
    } else {
        Ok(())
    }
}

impl<L: EquilibriumLattice3<Q>, const Q: usize> OpenBoundaryGrid3<L, Q> {
    /// Sample a tangential wall drive at each actual reflected-link midpoint.
    ///
    /// Coordinates are `p + 1/2 - e_q/2` in lattice cell units, including
    /// diagonal links and wall/open seams. The callback supplies grid-frame
    /// velocity, not slip relative to the fluid. It runs only here, in x-fast
    /// fluid-cell/native-direction order; the admitted trace is held over the
    /// next intervals until explicitly updated. Caller owns physical time,
    /// SI conversion and resolution of an oscillating or spatially varying drive.
    ///
    /// Native tangency, Mach, intersecting-wall compatibility and local mass
    /// conservation are checked on the complete reflected set. A tangential
    /// field with unresolved surface divergence can therefore refuse; no wall
    /// mass is compensated through a reservoir. No normal geometry motion or
    /// interior-voxel motion is introduced. Impulse, torque, work and gradient
    /// consumers use these same link samples. Stationary finite-wall heat
    /// storage still refuses a nonzero trace because it has no solid advection.
    ///
    /// Two reused buffers contain Q vectors per face-plane cell, including
    /// inactive zero padding. `maximum_storage_bytes` bounds their combined
    /// payload over ALL attached wall profiles, excluding allocator metadata.
    /// There is no source allocation or callback during a fluid step. Refusal
    /// or cancellation preserves the old complete trace and all physical state;
    /// a later refused step retains an already accepted control update.
    pub fn set_wall_velocity_profile(
        &mut self,
        face: Face3,
        maximum_storage_bytes: usize,
        gate: &CancelGate,
        mut sample: impl FnMut([f64; 3]) -> Result<[f64; 3], GridError3<L::Error>>,
    ) -> Result<WallVelocityProfileReport3, GridError3<L::Error>> {
        check(gate)?;
        if !matches!(self.plan.boundaries.face(face), FaceBoundary3::Wall { .. }) {
            return Err(GridError3::Input("wall profile requires a wall face"));
        }
        let axis = face as usize / 2;
        let count = (self.grid.f.len() / self.grid.dims[axis])
            .checked_mul(Q)
            .ok_or(GridError3::Allocation)?;
        let vectors = self
            .wall_profiles
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != face as usize)
            .try_fold(count, |n, (_, profile)| {
                n.checked_add(profile.as_ref().map_or(0, |p| p.values.len()))
            })
            .ok_or(GridError3::Allocation)?;
        let bytes = vectors
            .checked_mul(2 * core::mem::size_of::<[f64; 3]>())
            .ok_or(GridError3::Allocation)?;
        if bytes > maximum_storage_bytes {
            return Err(GridError3::Input(
                "wall profiles exceed their storage budget",
            ));
        }
        let index = face as usize;
        let had_profile = self.wall_profiles[index].is_some();
        let mut profile = match self.wall_profiles[index].take() {
            Some(profile) => profile,
            None => {
                let mut values = Vec::new();
                let mut pending = Vec::new();
                values
                    .try_reserve_exact(count)
                    .map_err(|_| GridError3::Allocation)?;
                pending
                    .try_reserve_exact(count)
                    .map_err(|_| GridError3::Allocation)?;
                values.resize(count, [0.0; 3]);
                pending.resize(count, [0.0; 3]);
                WallVelocityProfile3 {
                    values,
                    pending,
                    moving: false,
                }
            }
        };
        let outcome = (|| {
            let mut report = WallVelocityProfileReport3 {
                face,
                sampled_links: 0,
                maximum_speed_lu: 0.0,
            };
            let normal = if index % 2 == 0 { 1 } else { -1 };
            for cell in super::super::velocity::plane_cells(self.grid.dims, face) {
                check(gate)?;
                if self.grid.solid[cell] {
                    continue;
                }
                let p = coordinates(cell, self.grid.dims);
                for (q, e) in L::E.iter().enumerate().skip(1) {
                    if e[axis] != normal {
                        continue;
                    }
                    let point = core::array::from_fn(|a| p[a] as f64 + 0.5 - 0.5 * f64::from(e[a]));
                    let velocity = sample(point)?;
                    if velocity[axis] != 0.0 {
                        return Err(GridError3::Input(
                            "fixed exterior wall velocity must be tangential",
                        ));
                    }
                    initial::admit_velocity::<L, Q>(
                        self.grid.collision,
                        velocity,
                        cell,
                        self.grid.budget.maximum_mach,
                    )?;
                    profile.pending[plane_index(self.grid.dims, face, cell) * Q + q] = velocity;
                    report.sampled_links += 1;
                    report.maximum_speed_lu = report
                        .maximum_speed_lu
                        .max(velocity[0].hypot(velocity[1]).hypot(velocity[2]));
                }
            }
            if report.sampled_links == 0 {
                return Err(GridError3::Input("wall profile has no fluid links"));
            }
            let mut values = wall_profile_values(&self.wall_profiles);
            values[index] = Some(&profile.pending);
            let plan = FacePlan::admit(&self.grid, self.plan.boundaries, values, &mut || {
                check(gate)
            })?;
            check(gate)?;
            Ok((report, plan))
        })();
        match outcome {
            Ok((report, plan)) => {
                core::mem::swap(&mut profile.values, &mut profile.pending);
                profile.moving = report.maximum_speed_lu != 0.0;
                self.wall_profiles[index] = Some(profile);
                self.plan = plan;
                Ok(report)
            }
            Err(error) => {
                if had_profile {
                    self.wall_profiles[index] = Some(profile);
                }
                Err(error)
            }
        }
    }

    /// Restore a held uniform wall only if the complete remaining wall trace
    /// still has compatible corners and zero local mass source. No impulse,
    /// population reset or physical time advance occurs during this update.
    pub fn clear_wall_velocity_profile(
        &mut self,
        face: Face3,
        gate: &CancelGate,
    ) -> Result<(), GridError3<L::Error>> {
        check(gate)?;
        let mut values = wall_profile_values(&self.wall_profiles);
        values[face as usize] = None;
        let plan = FacePlan::admit(&self.grid, self.plan.boundaries, values, &mut || {
            check(gate)
        })?;
        check(gate)?;
        self.wall_profiles[face as usize] = None;
        self.plan = plan;
        Ok(())
    }

    /// Plane cells in x-fast order, each followed by Q native direction slots.
    /// Only incoming exterior links at fluid cells are active; other slots are
    /// zero padding. Readback never samples or advances the prescribed drive.
    pub fn wall_velocity_profile(&self, face: Face3) -> Option<&[[f64; 3]]> {
        self.wall_profiles[face as usize]
            .as_ref()
            .map(|p| p.values.as_slice())
    }
}

impl<L: EquilibriumLattice3<Q>, const Q: usize> ForcedOpenBoundaryGrid3<L, Q> {
    /// Bind the same held halfway-link trace without changing force history.
    pub fn set_wall_velocity_profile(
        &mut self,
        face: Face3,
        maximum_storage_bytes: usize,
        gate: &CancelGate,
        sample: impl FnMut([f64; 3]) -> Result<[f64; 3], GridError3<L::Error>>,
    ) -> Result<WallVelocityProfileReport3, GridError3<L::Error>> {
        self.inner
            .set_wall_velocity_profile(face, maximum_storage_bytes, gate, sample)
    }

    pub fn clear_wall_velocity_profile(
        &mut self,
        face: Face3,
        gate: &CancelGate,
    ) -> Result<(), GridError3<L::Error>> {
        self.inner.clear_wall_velocity_profile(face, gate)
    }

    pub fn wall_velocity_profile(&self, face: Face3) -> Option<&[[f64; 3]]> {
        self.inner.wall_velocity_profile(face)
    }
}
