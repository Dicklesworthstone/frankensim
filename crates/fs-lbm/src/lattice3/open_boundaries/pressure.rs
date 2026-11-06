//! Held spatial pressure on the original native D3Q19/D3Q27 reservoirs.

use super::*;
use crate::d3q19::{BoundaryProfileError3, PressureProfileReport3};

pub(super) struct PressureProfile3 {
    pub(super) values: Vec<f64>,
    pending: Vec<f64>,
}

fn check(gate: &CancelGate) -> Result<(), BoundaryProfileError3> {
    if gate.is_requested() {
        Err(BoundaryProfileError3::Cancelled)
    } else {
        Ok(())
    }
}

impl<L: EquilibriumLattice3<Q>, const Q: usize> OpenBoundaryGrid3<L, Q> {
    /// Replace one reservoir's held spatial pressure-density field atomically.
    ///
    /// `sample` receives native cell indices in x-fast order, once per fluid
    /// cell on the selected pressure plane. Pressure is `L::CS2 * density` in
    /// lattice units. SI pressure perturbations use the existing material
    /// scaling; these densities are not absolute thermodynamic pressures.
    ///
    /// This admits transverse hydrostatic heads and localized pressure loads
    /// with the original D3Q19/D3Q27 collision, force frame, and reservoir
    /// mass/momentum/angular-momentum accounting. No callback runs in a step.
    /// Two reusable scalar plane buffers keep storage proportional to face
    /// area. Fixed topology means solid slots are never sampled or consumed.
    ///
    /// Both serial and pooled steps use this same held field. The original
    /// first-interior velocity/stress and wall/open ownership are retained;
    /// only open-owned populations change at a wall seam, where the actual
    /// density need not equal the complete prescribed target. The step receipt
    /// accounts for actual transport, not an inferred pressure-work balance.
    ///
    /// A source refusal, invalid density or cancellation leaves the previous
    /// complete trace, populations and time unchanged. This update itself
    /// advances neither fluid nor source time; callers own temporal sampling.
    pub fn set_pressure_profile(
        &mut self,
        face: Face3,
        gate: &CancelGate,
        mut sample: impl FnMut([usize; 3]) -> Result<f64, BoundaryProfileError3>,
    ) -> Result<PressureProfileReport3, BoundaryProfileError3> {
        check(gate)?;
        let FaceBoundary3::Pressure { density } = self.plan.boundaries.face(face) else {
            return Err(BoundaryProfileError3::NotPressureFace);
        };
        let plane_len = self.grid.f.len() / self.grid.dims[face as usize / 2];
        let index = face as usize;
        let had_profile = self.pressure_profiles[index].is_some();
        let mut profile = match self.pressure_profiles[index].take() {
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
                values.resize(plane_len, density);
                pending.resize(plane_len, density);
                PressureProfile3 { values, pending }
            }
        };
        let outcome = (|| {
            let mut count = 0;
            let mut minimum = f64::INFINITY;
            let mut maximum = 0.0_f64;
            for (ordinal, cell) in self.plan.cells().enumerate() {
                if ordinal.is_multiple_of(64) {
                    check(gate)?;
                }
                if self.grid.solid[cell]
                    || self
                        .plan
                        .locate(cell)
                        .is_none_or(|(candidate, _)| candidate != face)
                {
                    continue;
                }
                // Every native fluid face cell owns its straight normal
                // incoming link, including corners. Donors were admitted at
                // construction, and this owner never changes solid topology.
                let point = coordinates(cell, self.grid.dims);
                let value = sample(point)?;
                if !value.is_finite() || value <= 0.0 {
                    return Err(BoundaryProfileError3::InvalidDensity { cell: point });
                }
                profile.pending[self.plan.plane_index(cell)] = value;
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
            core::mem::swap(&mut profile.values, &mut profile.pending);
            self.pressure_profiles[index] = Some(profile);
        } else if had_profile {
            self.pressure_profiles[index] = Some(profile);
        }
        outcome
    }

    /// Return to the face's held uniform pressure without advancing fluid time.
    pub fn clear_pressure_profile(&mut self, face: Face3) {
        self.pressure_profiles[face as usize] = None;
    }

    /// Held plane in x-fast order with the normal coordinate omitted.
    /// Inactive solid slots carry unused values; this observation is pure.
    #[must_use]
    pub fn pressure_profile(&self, face: Face3) -> Option<&[f64]> {
        self.pressure_profiles[face as usize]
            .as_ref()
            .map(|profile| profile.values.as_slice())
    }
}
