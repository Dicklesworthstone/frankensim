//! Resolved velocity traces on the original native reservoir operator.

use super::*;
use crate::d3q19::{BoundaryProfileError3, InletFlux3, VelocityProfileReport3};

pub(super) struct VelocityProfile3 {
    pub(super) values: Vec<[f64; 3]>,
    pending: Vec<[f64; 3]>,
}

fn check(gate: &CancelGate) -> Result<(), BoundaryProfileError3> {
    if gate.is_requested() {
        Err(BoundaryProfileError3::Cancelled)
    } else {
        Ok(())
    }
}

// Only the requested plane, in x-fast order with its fixed axis omitted.
pub(super) fn plane_cells(dims: [usize; 3], face: Face3) -> impl Iterator<Item = usize> {
    let axis = face as usize / 2;
    let stride = dims[..axis].iter().product::<usize>();
    let block = stride * dims[axis];
    let offset = if face as usize % 2 == 0 {
        0
    } else {
        (dims[axis] - 1) * stride
    };
    (0..dims.iter().product::<usize>() / block).flat_map(move |outer| {
        let start = outer * block + offset;
        start..start + stride
    })
}

impl<L: EquilibriumLattice3<Q>, const Q: usize> OpenBoundaryGrid3<L, Q> {
    /// Atomically replace one native velocity face's held spatial trace.
    ///
    /// `sample` receives fluid-node lattice indices in x-fast plane order.
    /// Sample once at the caller's chosen physical time, then hold these values
    /// during collision/streaming; no callback or profile allocation occurs in
    /// a step. Repeated updates reuse two face-sized vector buffers. Sampling
    /// does not change populations, forces, wall reactions or accepted time.
    /// A later refused fluid step retains this already accepted control.
    ///
    /// The existing force-corrected donor stress, native equilibrium, reservoir
    /// receipt and wall/open link ownership are unchanged. At mixed seams the
    /// actual complete moments need not equal the prescribed trace. Mean-normal
    /// control refuses such seams, as on the specialized D3Q19 owner. Its mean
    /// is an authored volume-velocity constraint, not measured mass transport.
    /// Physical-unit conversion and temporal/spectral resolution remain with
    /// the caller. No turbulent-spectrum or nonreflecting-outlet claim is added.
    pub fn set_velocity_profile(
        &mut self,
        face: Face3,
        flux: InletFlux3,
        gate: &CancelGate,
        sample: impl FnMut([usize; 3]) -> Result<[f64; 3], BoundaryProfileError3>,
    ) -> Result<VelocityProfileReport3, BoundaryProfileError3> {
        self.set_velocity_profile_inner(face, flux, None, gate, sample)
    }

    /// Also bound the FINAL vector speed after optional mean-normal control.
    /// Samples and corrected targets must satisfy the native grid/collision
    /// Mach envelope. The additional inclusive limit is in lattice speed units;
    /// zero admits only a quiescent corrected trace. Failure preserves the old
    /// complete profile, including cancellation after the last source sample.
    pub fn set_velocity_profile_with_speed_limit(
        &mut self,
        face: Face3,
        flux: InletFlux3,
        maximum_speed_lu: f64,
        gate: &CancelGate,
        sample: impl FnMut([usize; 3]) -> Result<[f64; 3], BoundaryProfileError3>,
    ) -> Result<VelocityProfileReport3, BoundaryProfileError3> {
        check(gate)?;
        if !maximum_speed_lu.is_finite() || maximum_speed_lu < 0.0 {
            return Err(BoundaryProfileError3::InvalidSpeedLimit);
        }
        self.set_velocity_profile_inner(face, flux, Some(maximum_speed_lu), gate, sample)
    }

    fn set_velocity_profile_inner(
        &mut self,
        face: Face3,
        flux: InletFlux3,
        maximum_speed_lu: Option<f64>,
        gate: &CancelGate,
        mut sample: impl FnMut([usize; 3]) -> Result<[f64; 3], BoundaryProfileError3>,
    ) -> Result<VelocityProfileReport3, BoundaryProfileError3> {
        check(gate)?;
        let FaceBoundary3::Velocity { velocity: uniform } = self.plan.boundaries.face(face) else {
            return Err(BoundaryProfileError3::NotVelocityFace);
        };
        if let InletFlux3::MeanNormal(mean) = flux
            && (!mean.is_finite() || mean * mean >= 0.03)
        {
            return Err(BoundaryProfileError3::InvalidMean);
        }
        let dims = self.grid.dims;
        let axis = face as usize / 2;
        let sign = if face as usize % 2 == 0 { -1.0 } else { 1.0 };
        let plane_len = self.grid.f.len() / dims[axis];
        let index = face as usize;
        let had_profile = self.velocity_profiles[index].is_some();
        let mut profile = match self.velocity_profiles[index].take() {
            Some(profile) => profile,
            None => {
                let mut values = Vec::new();
                let mut pending = Vec::new();
                values
                    .try_reserve_exact(plane_len)
                    .map_err(|_| BoundaryProfileError3::Allocation)?;
                pending
                    .try_reserve_exact(plane_len)
                    .map_err(|_| BoundaryProfileError3::Allocation)?;
                values.resize(plane_len, uniform);
                pending.resize(plane_len, uniform);
                VelocityProfile3 { values, pending }
            }
        };
        let outcome = (|| {
            let data = StepData::<L, Q> {
                dims,
                periodic: self.grid.periodic,
                solid: &self.grid.solid,
                collision: self.grid.collision,
                collisions: None,
                force: self.grid.force,
                forces: None,
                budget: self.grid.budget,
                marker: PhantomData,
            };
            let admit = |cell: usize, value: [f64; 3]| {
                let speed_squared = value.iter().map(|v| v * v).sum::<f64>();
                if value.iter().any(|v| !v.is_finite())
                    || speed_squared >= 0.03
                    || initial::admit_velocity::<L, Q>(
                        self.grid.collision,
                        value,
                        cell,
                        self.grid.budget.maximum_mach,
                    )
                    .is_err()
                {
                    Err(BoundaryProfileError3::InvalidVelocity {
                        cell: coordinates(cell, dims),
                    })
                } else {
                    Ok(speed_squared)
                }
            };
            let mut count = 0;
            let mut total = [Dd::ZERO; 3];
            let mut offset = Dd::ZERO;
            for (ordinal, cell) in plane_cells(dims, face).enumerate() {
                if ordinal.is_multiple_of(64) {
                    check(gate)?;
                }
                if self.grid.solid[cell] {
                    continue;
                }
                if matches!(flux, InletFlux3::MeanNormal(_)) && self.plan.masks(&data, cell).0 != 0
                {
                    return Err(BoundaryProfileError3::MixedWallSeam);
                }
                let value = sample(coordinates(cell, dims))?;
                admit(cell, value)?;
                profile.pending[self.plan.plane_index(cell)] = value;
                count += 1;
                for a in 0..3 {
                    total[a] = total[a] + Dd::from_f64(value[a]);
                }
                if let InletFlux3::MeanNormal(target) = flux {
                    // Uniform samples equal to the target retain exact zero
                    // correction; compensate cancellation of varying samples.
                    offset = offset + Dd::from_f64(sign * value[axis]) - Dd::from_f64(target);
                }
            }
            if count == 0 {
                return Err(BoundaryProfileError3::EmptyFace);
            }
            let correction = -offset.to_f64() / count as f64;
            let mut corrected = [Dd::ZERO; 3];
            let mut maximum_speed_squared: f64 = 0.0;
            for (ordinal, cell) in plane_cells(dims, face).enumerate() {
                if ordinal.is_multiple_of(64) {
                    check(gate)?;
                }
                if self.grid.solid[cell] {
                    continue;
                }
                let value = &mut profile.pending[self.plan.plane_index(cell)];
                if correction != 0.0 {
                    value[axis] += sign * correction;
                }
                let speed_squared = admit(cell, *value)?;
                if maximum_speed_lu
                    .is_some_and(|limit| value[0].hypot(value[1]).hypot(value[2]) > limit)
                {
                    return Err(BoundaryProfileError3::InvalidVelocity {
                        cell: coordinates(cell, dims),
                    });
                }
                maximum_speed_squared = maximum_speed_squared.max(speed_squared);
                for a in 0..3 {
                    corrected[a] = corrected[a] + Dd::from_f64(value[a]);
                }
            }
            check(gate)?;
            Ok(VelocityProfileReport3 {
                face,
                sampled_cells: count,
                sampled_mean_velocity_lu: total.map(|v| v.to_f64() / count as f64),
                prescribed_mean_velocity_lu: corrected.map(|v| v.to_f64() / count as f64),
                normal_velocity_correction_lu: correction,
                maximum_speed_squared_lu: maximum_speed_squared,
            })
        })();
        if outcome.is_ok() {
            core::mem::swap(&mut profile.values, &mut profile.pending);
            self.velocity_profiles[index] = Some(profile);
        } else if had_profile {
            self.velocity_profiles[index] = Some(profile);
        }
        outcome
    }

    /// Restore this face's held uniform target; populations and time are unchanged.
    pub fn clear_velocity_profile(&mut self, face: Face3) {
        self.velocity_profiles[face as usize] = None;
    }

    /// Held vector plane in x-fast order with the normal coordinate omitted.
    /// Solid slots are inactive. Readback does not sample or advance any source.
    #[must_use]
    pub fn velocity_profile(&self, face: Face3) -> Option<&[[f64; 3]]> {
        self.velocity_profiles[face as usize]
            .as_ref()
            .map(|profile| profile.values.as_slice())
    }
}
