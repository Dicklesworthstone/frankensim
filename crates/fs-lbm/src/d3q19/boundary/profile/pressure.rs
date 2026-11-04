//! Spatial reservoir pressure targets through the original density boundary.

use super::*;

/// Diagnostics of the held pressure trace, not measured transport or work.
/// Native pressure is represented by positive density in lattice units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PressureProfileReport3 {
    pub face: Face3,
    /// Number of active fluid nodes sampled on the selected face.
    pub sampled_cells: usize,
    pub minimum_density_lu: f64,
    pub maximum_density_lu: f64,
}

impl BoundaryGrid3 {
    /// Atomically replace a pressure face's held spatial density trace.
    ///
    /// Sampling uses deterministic z/y/x order at active fluid node indices.
    /// Each sample must be positive and finite, exactly as for a uniform
    /// pressure target. The first successful profile locks solid topology.
    /// Repeated calls reuse two scalar plane buffers; refusal or cancellation
    /// retains the entire previous trace and accepted populations.
    ///
    /// The original pressure reconstruction retains the first-interior
    /// velocity/stress and Guo frame. At wall/open seams only open-owned
    /// populations change, so actual density need not equal the full trace.
    /// Actual mass and momentum belong to the next step's shared open receipt.
    /// This neither advances time nor supplies a pressure-work balance.
    pub fn set_pressure_profile(
        &mut self,
        face: Face3,
        gate: &CancelGate,
        mut sample: impl FnMut([usize; 3]) -> Result<f64, BoundaryProfileError3>,
    ) -> Result<PressureProfileReport3, BoundaryProfileError3> {
        check(gate)?;
        let FaceBoundary3::Pressure { density: uniform } = self.boundaries.face(face) else {
            return Err(BoundaryProfileError3::NotPressureFace);
        };
        let dims = self.dimensions();
        let plane_len = dims.iter().product::<usize>() / dims[face.axis()];
        let had_profile = self.pressure_profiles[face.index()].is_some();
        let mut state = match self.pressure_profiles[face.index()].take() {
            Some(state) => state,
            None => {
                let mut values = Vec::new();
                let mut scratch = Vec::new();
                values.try_reserve_exact(plane_len).map_err(|_| BoundaryProfileError3::Allocation)?;
                scratch.try_reserve_exact(plane_len).map_err(|_| BoundaryProfileError3::Allocation)?;
                values.resize(plane_len, uniform);
                scratch.resize(plane_len, uniform);
                ProfileState3 { values, scratch }
            }
        };
        let outcome = (|| {
            let mut count = 0;
            let mut minimum = f64::INFINITY;
            let mut maximum = 0.0_f64;
            for (index, cell) in planar_face_cells(dims, face).enumerate() {
                if index.is_multiple_of(64) {
                    check(gate)?;
                }
                if self.is_solid(cell[0], cell[1], cell[2])
                    || self.open_link_mask(cell[0], cell[1], cell[2]) == 0
                {
                    continue;
                }
                let interior = Self::interior_neighbor(cell[0], cell[1], cell[2], face);
                if self.is_solid(interior[0], interior[1], interior[2]) {
                    return Err(BoundaryProfileError3::InvalidInterior { cell });
                }
                let value = sample(cell)?;
                if !value.is_finite() || value <= 0.0 {
                    return Err(BoundaryProfileError3::InvalidDensity { cell });
                }
                state.scratch[index] = value;
                count += 1;
                minimum = minimum.min(value);
                maximum = maximum.max(value);
            }
            if count == 0 {
                return Err(BoundaryProfileError3::EmptyFace);
            }
            check(gate)?;
            Ok(PressureProfileReport3 {
                face,
                sampled_cells: count,
                minimum_density_lu: minimum,
                maximum_density_lu: maximum,
            })
        })();
        if outcome.is_ok() {
            core::mem::swap(&mut state.values, &mut state.scratch);
            self.pressure_profiles[face.index()] = Some(state);
            self.topology_locked = true;
        } else if had_profile {
            self.pressure_profiles[face.index()] = Some(state);
        }
        outcome
    }

    /// Restore the held uniform density target without advancing populations.
    pub fn clear_pressure_profile(&mut self, face: Face3) {
        self.pressure_profiles[face.index()] = None;
    }

    /// Full plane in z/y/x order, omitting the fixed normal coordinate.
    /// Inactive solid slots carry unused values. This observation is pure.
    #[must_use]
    pub fn pressure_profile(&self, face: Face3) -> Option<&[f64]> {
        self.pressure_profiles[face.index()].as_ref().map(|s| s.values.as_slice())
    }

    pub(in super::super) fn profile_density(&self, face: Face3, cell: [usize; 3]) -> Option<f64> {
        let index = match face.axis() {
            0 => cell[2] * self.ny + cell[1],
            1 => cell[2] * self.nx + cell[0],
            _ => cell[1] * self.nx + cell[0],
        };
        self.pressure_profiles[face.index()].as_ref().map(|s| s.values[index])
    }
}
