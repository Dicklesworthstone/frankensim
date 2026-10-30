//! A separately driven aperture within one original characteristic face.
//! The complementary face cells retain their explicit whole-face drive. All
//! reconstruction, incoming populations and FV fluxes still have one owner.

use super::*;

/// Half-open cell bounds on the two increasing Cartesian axes other than the
/// characteristic axis. Bounds select original boundary quadrature cells;
/// there is no inferred area, interpolated edge or geometry repair.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NsfCharacteristicPatch {
    pub start: [usize; 2],
    pub end: [usize; 2],
}

impl NsfCharacteristicPatch {
    pub(super) fn contains(self, cell: [usize; 3], normal_axis: usize) -> bool {
        let mut tangent = 0;
        for (axis, coordinate) in cell.into_iter().enumerate() {
            if axis == normal_axis {
                continue;
            }
            if coordinate < self.start[tangent] || coordinate >= self.end[tangent] {
                return false;
            }
            tangent += 1;
        }
        true
    }

    fn validate(self, dims: [usize; 3], normal_axis: usize) -> Result<(), NsfError> {
        let mut tangent = 0;
        for (axis, size) in dims.into_iter().enumerate() {
            if axis == normal_axis {
                continue;
            }
            if self.start[tangent] >= self.end[tangent] || self.end[tangent] > size {
                return Err(NsfError::Input(
                    "characteristic patch needs nonempty in-grid cell bounds",
                ));
            }
            tangent += 1;
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug)]
pub(in crate::nsf3) struct CharacteristicPatchBinding {
    pub high: bool,
    pub shape: NsfCharacteristicPatch,
    pub drive: NsfBoundaryDrive,
}

impl OpenFaces {
    pub(in crate::nsf3) fn contains_patch(&self, cell: [usize; 3], side: usize) -> bool {
        self.patch.is_some_and(|patch| {
            usize::from(patch.high) == side && patch.shape.contains(cell, self.axis)
        })
    }
}

impl<L: Lattice3<27>> NsfGrid3<L> {
    /// Bind one rectangular aperture before the first step. Other cells keep
    /// the original face control, including an explicitly selected ambient
    /// pressure. The patch shares that face's material/reference state and
    /// original two-cell fluid collar. This does not introduce a side wall,
    /// a second characteristic axis, nonlinear NSCBC or a turbulence closure.
    /// Invalid geometry/control leaves every owner unpublished.
    pub fn with_characteristic_patch(
        mut self,
        high: bool,
        patch: NsfCharacteristicPatch,
        drive: NsfBoundaryDrive,
    ) -> Result<Self, NsfError> {
        let open = self.open.ok_or(NsfError::Input(
            "a patch requires an original characteristic axis",
        ))?;
        if self.steps != 0 || self.characteristic_patch.is_some() || !drive.is_finite() {
            return Err(NsfError::Input(
                "one cold characteristic patch requires a finite drive",
            ));
        }
        patch.validate(self.spec.dims, open.axis)?;
        self.rebind_characteristics(
            open,
            self.characteristic_pressure_offset_pa,
            Some(CharacteristicPatchBinding {
                high,
                shape: patch,
                drive,
            }),
        )?;
        Ok(self)
    }

    /// Selected face and original quadrature-cell bounds.
    pub fn characteristic_patch(&self) -> Option<(bool, NsfCharacteristicPatch)> {
        self.characteristic_patch
            .map(|patch| (patch.high, patch.shape))
    }

    /// Accepted aperture control. Whole-face controls describe the complement.
    pub fn characteristic_patch_drive(&self) -> Option<NsfBoundaryDrive> {
        self.characteristic_patch.map(|patch| patch.drive)
    }

    /// Atomically change both whole-face drives and the existing patch drive.
    /// Physical density, velocity and total energy remain unchanged; the
    /// original force chart and boundary flux history are rebased together.
    pub fn set_boundary_drives_with_patch(
        &mut self,
        drives: [NsfBoundaryDrive; 2],
        drive: NsfBoundaryDrive,
    ) -> Result<(), NsfError> {
        let mut patch = self
            .characteristic_patch
            .ok_or(NsfError::Input("no characteristic patch is attached"))?;
        if !drive.is_finite() {
            return Err(NsfError::Input("characteristic patch drive must be finite"));
        }
        patch.drive = drive;
        self.set_boundary_drives_with_binding(drives, Some(patch))
    }

    /// Original accepted traces restricted to the aperture, in original order.
    pub fn open_patch_traces(&self) -> Option<Vec<NsfBoundaryTrace>> {
        let patch = self.characteristic_patch?;
        let axis = self.open.as_ref()?.axis;
        Some(
            self.open_boundary_traces(patch.high)?
                .into_iter()
                .filter(|trace| patch.shape.contains(trace.cell, axis))
                .collect(),
        )
    }

    /// Aperture exchange from the last accepted step. Whole-face accounts
    /// already include this exchange; adding it again would double count it.
    pub const fn last_patch_exchange(&self) -> Option<&NsfBoundaryExchange> {
        self.last_patch_exchange.as_ref()
    }

    pub(in crate::nsf3) fn finish_patch_exchange(
        &self,
        links: Option<([LinkExchange; 2], Option<LinkExchange>)>,
        before: Option<&OpenFaces>,
        after: Option<&OpenFaces>,
        flux_after: BoundaryFlux,
    ) -> Option<NsfBoundaryExchange> {
        let patch = self.characteristic_patch?;
        let link = links?.1?;
        let (before, after) = (before?, after?);
        let axis = before.axis;
        let side = usize::from(patch.high);
        let sign = if patch.high { 1.0 } else { -1.0 };
        let edge = if patch.high {
            self.spec.dims[axis] - 1
        } else {
            0
        };
        let reference = self.open.as_ref()?.faces[side].reference_pressure_pa;
        let mass_unit = self.law.rho_ref_kg_m3 * self.law.dx_m.powi(3);
        let velocity_unit = self.law.dx_m / self.law.dt_s;
        let area = self.law.dx_m * self.law.dx_m;
        let mut out = NsfBoundaryExchange {
            mass_out_kg: link.mass * mass_unit,
            momentum_out_kg_m_s: link.momentum.map(|p| p * mass_unit * velocity_unit),
            total_energy_out_j: 0.5
                * (self.open_out.patch_total_out + flux_after.patch_total_out)
                * self.law.energy_unit_j(),
            ..NsfBoundaryExchange::default()
        };
        let mut count = 0;
        for i in 0..self.n {
            let x = self.at(i);
            if x[axis] != edge || !patch.shape.contains(x, axis) {
                continue;
            }
            count += 1;
            let index = plane_index(x, self.spec.dims, axis);
            let old = before.values[side][index];
            let new = after.values[side][index];
            out.momentum_out_kg_m_s[axis] += 0.5
                * sign
                * (old.rho * (old.theta - THETA0) + new.rho * (new.theta - THETA0))
                * mass_unit
                * velocity_unit;
            out.mean_pressure_pa += 0.5 * (old.pressure_pa + new.pressure_pa);
            let old_acoustic = old.pressure_pa - reference;
            let new_acoustic = new.pressure_pa - reference;
            out.mean_acoustic_pressure_pa += 0.5 * (old_acoustic + new_acoustic);
            out.acoustic_pressure_work_out_j += 0.5
                * sign
                * (old_acoustic * old.u[axis] + new_acoustic * new.u[axis])
                * velocity_unit
                * area
                * self.law.dt_s;
            out.normal_volume_flow_m3_s +=
                0.5 * sign * (old.u[axis] + new.u[axis]) * velocity_unit * area;
            out.pressure_work_out_j += 0.5
                * sign
                * (old.pressure_pa * old.u[axis] + new.pressure_pa * new.u[axis])
                * velocity_unit
                * area
                * self.law.dt_s;
        }
        out.mean_pressure_pa /= count as f64;
        out.mean_acoustic_pressure_pa /= count as f64;
        Some(out)
    }
}
