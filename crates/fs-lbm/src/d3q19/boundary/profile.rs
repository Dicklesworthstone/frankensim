//! Spatial inlet traces through the existing regularized population boundary.
//! The caller owns physical time and the coherent field; sampling occurs once
//! per active face node, never inside a collision worker or observer query.

use super::{BoundaryGrid3, Face3, FaceBoundary3, planar_face_cells};
use fs_exec::CancelGate;

/// Control of the prescribed normal velocity, in lattice units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum InletFlux3 {
    /// Keep the supplied profile, including its time-dependent bulk flow.
    Unconstrained,
    /// Remove the plane-mean normal fluctuation and impose this signed
    /// outward-normal mean velocity. Inflow is negative on either face.
    /// This is a volumetric constraint on the prescribed trace, not a
    /// density-weighted mass-flux constraint or an exact population flux.
    /// Mixed wall/open nodes are refused for this mode because their partial
    /// reconstruction does not impose all three velocity moments.
    MeanNormal(f64),
}

/// Admission or cancellation leaves the previous complete profile in place.
#[derive(Clone, Debug, PartialEq)]
pub enum BoundaryProfileError3 {
    /// The requested face is not a declared velocity boundary.
    NotVelocityFace,
    /// No fluid open nodes are available on this face.
    EmptyFace,
    /// Mean-normal control is undefined at partially reconstructed wall seams.
    MixedWallSeam,
    /// The requested mean is nonfinite or outside the existing low-Mach bound.
    InvalidMean,
    /// A sampled or corrected velocity violates the existing low-Mach bound.
    InvalidVelocity { cell: [usize; 3] },
    /// An open node has a solid first-interior neighbor.
    InvalidInterior { cell: [usize; 3] },
    /// Callback refusal, retaining the original field's reason.
    Source(String),
    /// A reusable plane buffer could not be allocated.
    Allocation,
    /// The caller cancelled before publication.
    Cancelled,
}

impl core::fmt::Display for BoundaryProfileError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-LBM-BOUNDARY-PROFILE3: {self:?}")
    }
}
impl std::error::Error for BoundaryProfileError3 {}

/// Diagnostics of the prescribed velocity trace, not measured lattice fluxes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VelocityProfileReport3 {
    pub face: Face3,
    /// Number of fluid nodes with an open link on the selected face.
    pub sampled_cells: usize,
    /// Plane mean before optional normal correction, in Cartesian axes.
    pub sampled_mean_velocity_lu: [f64; 3],
    /// Plane mean after correction, in Cartesian axes.
    pub prescribed_mean_velocity_lu: [f64; 3],
    /// Correction in the outward-normal direction, applied to every node.
    pub normal_velocity_correction_lu: f64,
    /// Largest prescribed speed squared (the admitted bound is strictly 0.03).
    pub maximum_speed_squared_lu: f64,
}

pub(super) struct ProfileState3 {
    values: Vec<[f64; 3]>,
    scratch: Vec<[f64; 3]>,
}

fn check(gate: &CancelGate) -> Result<(), BoundaryProfileError3> {
    if gate.is_requested() {
        Err(BoundaryProfileError3::Cancelled)
    } else {
        Ok(())
    }
}

fn speed_squared(value: [f64; 3]) -> f64 {
    value.iter().map(|v| v * v).sum()
}

impl BoundaryGrid3 {
    /// Atomically replace one velocity face's spatial trace. Call after
    /// voxelization and before the step whose boundary endpoint it supplies.
    /// The method does not advance populations, the field's clock, or RNGs.
    /// A successful first profile locks solid topology, as stepping does.
    ///
    /// `sample` is called in deterministic z/y/x order only for fluid nodes
    /// with open links. Coordinates are lattice node indices; physical-unit
    /// conversion, temporal evolution, wall compatibility, and spectral
    /// resolution belong to the caller. Repeated calls reuse both plane
    /// buffers; no volume-sized buffer or per-worker random stream is added.
    ///
    /// The original regularized boundary and wall/open link ownership remain
    /// in force. At mixed seams only open-owned populations are reconstructed;
    /// their actual moments need not equal the complete prescribed trace.
    /// This is not a new nonreflecting outlet or a qualified turbulent-inflow
    /// spectrum. The same profile is consumed by ordinary and SGS stepping.
    ///
    /// # Errors
    /// Typed geometry, source, low-Mach, allocation, or cancellation refusal.
    pub fn set_velocity_profile(
        &mut self,
        face: Face3,
        flux: InletFlux3,
        gate: &CancelGate,
        mut sample: impl FnMut([usize; 3]) -> Result<[f64; 3], BoundaryProfileError3>,
    ) -> Result<VelocityProfileReport3, BoundaryProfileError3> {
        check(gate)?;
        let FaceBoundary3::Velocity { velocity: uniform } = self.boundaries.face(face) else {
            return Err(BoundaryProfileError3::NotVelocityFace);
        };
        if let InletFlux3::MeanNormal(mean) = flux {
            if !mean.is_finite() || mean * mean >= 0.03 {
                return Err(BoundaryProfileError3::InvalidMean);
            }
        }
        let dims = self.dimensions();
        let axis = face.axis();
        let sign = f64::from(face.normal()[axis]);
        let plane_len = dims.iter().product::<usize>() / dims[axis];
        let had_profile = self.velocity_profiles[face.index()].is_some();
        let mut state = match self.velocity_profiles[face.index()].take() {
            Some(state) => state,
            None => {
                let mut values = Vec::new();
                let mut scratch = Vec::new();
                values
                    .try_reserve_exact(plane_len)
                    .map_err(|_| BoundaryProfileError3::Allocation)?;
                scratch
                    .try_reserve_exact(plane_len)
                    .map_err(|_| BoundaryProfileError3::Allocation)?;
                values.resize(plane_len, uniform);
                scratch.resize(plane_len, uniform);
                ProfileState3 { values, scratch }
            }
        };
        let outcome = (|| {
            let mut count = 0;
            let mut total = [0.0; 3];
            let mut normal_offset = 0.0;
            for (index, cell) in planar_face_cells(dims, face).enumerate() {
                if index.is_multiple_of(64) {
                    check(gate)?;
                }
                if self.is_solid(cell[0], cell[1], cell[2])
                    || self.open_link_mask(cell[0], cell[1], cell[2]) == 0
                {
                    continue;
                }
                if matches!(flux, InletFlux3::MeanNormal(_))
                    && self.link_mask(cell[0], cell[1], cell[2]) != 0
                {
                    return Err(BoundaryProfileError3::MixedWallSeam);
                }
                let interior = Self::interior_neighbor(cell[0], cell[1], cell[2], face);
                if self.is_solid(interior[0], interior[1], interior[2]) {
                    return Err(BoundaryProfileError3::InvalidInterior { cell });
                }
                let value = sample(cell)?;
                if !value.iter().all(|v| v.is_finite()) || speed_squared(value) >= 0.03 {
                    return Err(BoundaryProfileError3::InvalidVelocity { cell });
                }
                state.scratch[index] = value;
                count += 1;
                for a in 0..3 {
                    total[a] += value[a];
                }
                if let InletFlux3::MeanNormal(target) = flux {
                    // Sum fluctuations relative to target: a uniform control
                    // identical to target gets an exact zero correction.
                    normal_offset += sign * value[axis] - target;
                }
            }
            if count == 0 {
                return Err(BoundaryProfileError3::EmptyFace);
            }
            let correction = -normal_offset / count as f64;
            let mut corrected = [0.0; 3];
            let mut max_speed2 = 0.0_f64;
            for (index, cell) in planar_face_cells(dims, face).enumerate() {
                if index.is_multiple_of(64) {
                    check(gate)?;
                }
                if self.is_solid(cell[0], cell[1], cell[2])
                    || self.open_link_mask(cell[0], cell[1], cell[2]) == 0
                {
                    continue;
                }
                let value = &mut state.scratch[index];
                if correction != 0.0 {
                    value[axis] += sign * correction;
                }
                let speed2 = speed_squared(*value);
                if !value.iter().all(|v| v.is_finite()) || speed2 >= 0.03 {
                    return Err(BoundaryProfileError3::InvalidVelocity { cell });
                }
                max_speed2 = max_speed2.max(speed2);
                for a in 0..3 {
                    corrected[a] += value[a];
                }
            }
            check(gate)?;
            Ok(VelocityProfileReport3 {
                face,
                sampled_cells: count,
                sampled_mean_velocity_lu: total.map(|v| v / count as f64),
                prescribed_mean_velocity_lu: corrected.map(|v| v / count as f64),
                normal_velocity_correction_lu: correction,
                maximum_speed_squared_lu: max_speed2,
            })
        })();
        if outcome.is_ok() {
            core::mem::swap(&mut state.values, &mut state.scratch);
            self.velocity_profiles[face.index()] = Some(state);
            self.topology_locked = true;
        } else if had_profile {
            self.velocity_profiles[face.index()] = Some(state);
        }
        outcome
    }

    /// Restore the uniform velocity in the original face specification.
    /// Existing populations change only when the next lattice step applies it.
    pub fn clear_velocity_profile(&mut self, face: Face3) {
        self.velocity_profiles[face.index()] = None;
    }

    /// Full plane buffer in increasing z/y/x order, omitting the fixed normal
    /// coordinate. Inactive solid nodes carry unused values. Queries are pure.
    #[must_use]
    pub fn velocity_profile(&self, face: Face3) -> Option<&[[f64; 3]]> {
        self.velocity_profiles[face.index()]
            .as_ref()
            .map(|s| s.values.as_slice())
    }

    pub(super) fn profile_velocity(&self, face: Face3, cell: [usize; 3]) -> Option<[f64; 3]> {
        let index = match face.axis() {
            0 => cell[2] * self.ny + cell[1],
            1 => cell[2] * self.nx + cell[0],
            _ => cell[1] * self.nx + cell[0],
        };
        self.velocity_profiles[face.index()]
            .as_ref()
            .map(|s| s.values[index])
    }
}
