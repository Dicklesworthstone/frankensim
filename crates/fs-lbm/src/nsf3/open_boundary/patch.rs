//! Independently driven apertures within the original characteristic faces.
//! Complementary face cells retain their explicit whole-face drives. All
//! reconstruction, incoming populations and FV fluxes still have one owner.

use super::*;

#[cfg(test)]
mod tests;

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
        self.patch[side].is_some_and(|patch| patch.shape.contains(cell, self.axis))
    }
}

impl<L: Lattice3<27>> NsfGrid3<L> {
    /// Bind one rectangular aperture on the selected face before the first
    /// step. Each face may own one independent aperture. Other cells keep
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
        let side = usize::from(high);
        if self.steps != 0 || self.characteristic_patch[side].is_some() || !drive.is_finite() {
            return Err(NsfError::Input(
                "each face admits one cold characteristic patch with a finite drive",
            ));
        }
        patch.validate(self.spec.dims, open.axis)?;
        let mut patches = self.characteristic_patch;
        patches[side] = Some(CharacteristicPatchBinding { high, shape: patch, drive });
        self.rebind_characteristics(
            open,
            self.characteristic_pressure_offset_pa,
            self.characteristic_inflow_temperature_k,
            patches,
        )?;
        Ok(self)
    }

    /// Selected face and original quadrature-cell bounds when exactly one
    /// aperture is attached. Returns `None` for zero or two apertures; use
    /// [`Self::characteristic_patch_on`] to address both faces explicitly.
    pub fn characteristic_patch(&self) -> Option<(bool, NsfCharacteristicPatch)> {
        let side = self.sole_patch_side()?;
        self.characteristic_patch[side].map(|patch| (patch.high, patch.shape))
    }

    /// Original quadrature-cell bounds of the selected face's aperture.
    pub fn characteristic_patch_on(&self, high: bool) -> Option<NsfCharacteristicPatch> {
        self.characteristic_patch[usize::from(high)].map(|patch| patch.shape)
    }

    fn sole_patch_side(&self) -> Option<usize> {
        match self.characteristic_patch {
            [Some(_), None] => Some(0),
            [None, Some(_)] => Some(1),
            _ => None,
        }
    }

    pub(in crate::nsf3) fn single_patch_controls(
        &self, drive: NsfBoundaryDrive,
    ) -> Result<[Option<NsfBoundaryDrive>; 2], NsfError> {
        let side = self.sole_patch_side().ok_or(NsfError::Input(
            "a singular patch control requires exactly one attached aperture",
        ))?;
        let mut controls = [None; 2];
        controls[side] = Some(drive);
        Ok(controls)
    }

    /// Accepted control when exactly one aperture is attached; `None` for zero
    /// or two. Whole-face controls describe the complementary cells.
    pub fn characteristic_patch_drive(&self) -> Option<NsfBoundaryDrive> {
        self.characteristic_patch_drive_on(self.sole_patch_side()? == 1)
    }

    /// Accepted aperture control on the selected original face.
    pub fn characteristic_patch_drive_on(&self, high: bool) -> Option<NsfBoundaryDrive> {
        self.characteristic_patch[usize::from(high)].map(|patch| patch.drive)
    }

    /// Atomically change both whole-face drives and the existing patch drive.
    /// Physical density, velocity and total energy remain unchanged; the
    /// original force chart and boundary flux history are rebased together.
    pub fn set_boundary_drives_with_patch(
        &mut self,
        drives: [NsfBoundaryDrive; 2],
        drive: NsfBoundaryDrive,
    ) -> Result<(), NsfError> {
        let controls = self.single_patch_controls(drive)?;
        self.set_boundary_drives_with_patches(drives, controls)
    }

    /// Atomically change both complementary face drives and explicitly selected
    /// aperture drives. Entries are low then high; `None` retains the accepted
    /// aperture control and `Some` requires an aperture on that face. Invalid
    /// data on either face publishes no boundary, population or energy history.
    pub fn set_boundary_drives_with_patches(
        &mut self,
        drives: [NsfBoundaryDrive; 2],
        patch_drives: [Option<NsfBoundaryDrive>; 2],
    ) -> Result<(), NsfError> {
        let mut patches = self.characteristic_patch;
        for (patch, drive) in patches.iter_mut().zip(patch_drives) {
            if let Some(drive) = drive {
                if !drive.is_finite() {
                    return Err(NsfError::Input("characteristic patch drive must be finite"));
                }
                let patch = patch.as_mut().ok_or(NsfError::Input(
                    "a selected patch drive requires an aperture on the same face",
                ))?;
                patch.drive = drive;
            }
        }
        self.set_boundary_drives_with_binding(drives, patches)
    }

    /// Original traces when exactly one aperture is attached; `None` for zero
    /// or two. Use [`Self::open_patch_traces_on`] for explicit face selection.
    pub fn open_patch_traces(&self) -> Option<Vec<NsfBoundaryTrace>> {
        self.open_patch_traces_on(self.sole_patch_side()? == 1)
    }

    /// Original accepted traces restricted to one face's aperture, in original order.
    pub fn open_patch_traces_on(&self, high: bool) -> Option<Vec<NsfBoundaryTrace>> {
        let patch = self.characteristic_patch[usize::from(high)]?;
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
    /// Returns `None` if both faces carry apertures; use the explicit face accessor.
    pub const fn last_patch_exchange(&self) -> Option<&NsfBoundaryExchange> {
        match &self.last_patch_exchange {
            [Some(exchange), None] | [None, Some(exchange)] => Some(exchange),
            _ => None,
        }
    }

    /// Last accepted aperture exchange on the selected face, already included
    /// in that face's complete boundary account.
    pub const fn last_patch_exchange_on(&self, high: bool) -> Option<&NsfBoundaryExchange> {
        self.last_patch_exchange[high as usize].as_ref()
    }

    pub(in crate::nsf3) fn finish_patch_exchange(
        &self,
        links: Option<([LinkExchange; 2], [Option<LinkExchange>; 2])>,
        before: Option<&OpenFaces>,
        after: Option<&OpenFaces>,
        flux_after: BoundaryFlux,
    ) -> [Option<NsfBoundaryExchange>; 2] {
        std::array::from_fn(|side| self.finish_patch_exchange_on(
            side, links, before, after, flux_after,
        ))
    }

    fn finish_patch_exchange_on(
        &self,
        side: usize,
        links: Option<([LinkExchange; 2], [Option<LinkExchange>; 2])>,
        before: Option<&OpenFaces>,
        after: Option<&OpenFaces>,
        flux_after: BoundaryFlux,
    ) -> Option<NsfBoundaryExchange> {
        let patch = self.characteristic_patch[side]?;
        let link = links?.1[side]?;
        let (before, after) = (before?, after?);
        let axis = before.axis;
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
                * (self.open_out.patch_total_out[side] + flux_after.patch_total_out[side])
                * self.law.energy_unit_j(),
            advected_stagnation_enthalpy_out_j: 0.5
                * (self.open_out.patch_advected_out[side] + flux_after.patch_advected_out[side])
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
        out.nonadvective_energy_out_j =
            out.total_energy_out_j - out.advected_stagnation_enthalpy_out_j;
        Some(out)
    }
}
