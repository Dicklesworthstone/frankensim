//! Full Lighthill tensor from the original configured dissipative operator.

use super::*;
use fs_flux::compressible::spatial::{mixture::molecular as mixture, molecular as frozen};
use mixture::{DissipationControl, DissipationError, TransportLaw};

impl PreparedCompressibleVolume3 {
    /// Bind a frozen-gas molecular flow and publish its current accepted sample.
    ///
    /// The tensor includes `-tau_ij` from the original molecular/bulk transport
    /// and any explicitly configured SGS closure. It is evaluated at the
    /// accepted flow time, not averaged over the previous time step. The native
    /// operator reuses its allocated workspace; no acoustic transport law or
    /// independent velocity-gradient approximation is introduced.
    pub fn from_molecular_frozen<L: TransportLaw>(
        flow: &mut frozen::MolecularEuler3d<L>,
        cells: &[usize],
        config: CompressibleVolumeConfig3,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CompressibleVolumeError3> {
        Self::admit(
            flow.flow().grid(),
            cells,
            config,
            maximum_workspace_bytes,
            &mut cancelled,
        )?;
        let mut result = Self::prepare(
            Source::Frozen(flow.flow().gas().clone()),
            *flow.flow().grid(),
            flow.flow().boundaries(),
            cells,
            config,
            maximum_workspace_bytes,
            &mut cancelled,
        )?;
        result.with_viscous_stress = true;
        result.sample_molecular_frozen(flow, &mut cancelled)?;
        Ok(result)
    }

    /// Bind a full-composition molecular flow and its actual native stress.
    ///
    /// The transport source retains the original species/caloric owner and
    /// composition-dependent coefficients. The supplied controls are passed to
    /// that native operator unchanged; absent species are not filled in.
    pub fn from_molecular_mixture<L: TransportLaw>(
        flow: &mut mixture::MolecularEuler3d<L>,
        cells: &[usize],
        config: CompressibleVolumeConfig3,
        maximum_workspace_bytes: usize,
        controls: DissipationControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CompressibleVolumeError3> {
        Self::admit(
            flow.flow().grid(),
            cells,
            config,
            maximum_workspace_bytes,
            &mut cancelled,
        )?;
        let mut result = Self::prepare(
            Source::Mixture(flow.flow().gas().clone()),
            *flow.flow().grid(),
            flow.flow().boundaries(),
            cells,
            config,
            maximum_workspace_bytes,
            &mut cancelled,
        )?;
        result.with_viscous_stress = true;
        result.sample_molecular_mixture(flow, controls, &mut cancelled)?;
        Ok(result)
    }

    /// Publish the full tensor at a later accepted frozen-gas source time.
    pub fn sample_molecular_frozen<L: TransportLaw>(
        &mut self,
        flow: &mut frozen::MolecularEuler3d<L>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), CompressibleVolumeError3> {
        poll(&mut cancelled)?;
        if !self.with_viscous_stress || self.frozen_source() != Some(flow.flow().gas()) {
            return Err(CompressibleVolumeError3::SourceMismatch);
        }
        self.admit_sample(
            flow.flow().grid(),
            flow.flow().boundaries(),
            flow.flow().time().value(),
        )?;
        let field = flow.stress_field(&mut cancelled).map_err(molecular_error)?;
        self.sample_stress(&field, &mut cancelled)
    }

    /// Publish the full tensor at a later accepted mixture source time.
    pub fn sample_molecular_mixture<L: TransportLaw>(
        &mut self,
        flow: &mut mixture::MolecularEuler3d<L>,
        controls: DissipationControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), CompressibleVolumeError3> {
        poll(&mut cancelled)?;
        if !self.with_viscous_stress || self.mixture_source() != Some(flow.flow().gas()) {
            return Err(CompressibleVolumeError3::SourceMismatch);
        }
        self.admit_sample(
            flow.flow().grid(),
            flow.flow().boundaries(),
            flow.flow().time().value(),
        )?;
        let field = flow
            .stress_field(controls, &mut cancelled)
            .map_err(dissipation_error)?;
        self.sample_stress(&field, &mut cancelled)
    }

    /// Shared with the joint exterior recorder after it checks the original
    /// source identity. Both surface and volume can read one native endpoint
    /// evaluation, retaining their original source order and accepted clock.
    pub(super) fn sample_stress(
        &mut self,
        field: &frozen::StressField3d<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), CompressibleVolumeError3> {
        poll(cancelled)?;
        if !self.with_viscous_stress {
            return Err(CompressibleVolumeError3::SourceMismatch);
        }
        self.admit_sample(field.grid(), field.boundaries(), field.time_s())?;
        for (i, &cell) in self.cells.iter().enumerate() {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            let p = field
                .primitive(cell)
                .ok_or(CompressibleVolumeError3::GeometryMismatch)?;
            let stress = field
                .cell_stress_pa(cell)
                .ok_or(CompressibleVolumeError3::GeometryMismatch)?;
            self.trial[i] = kernel::molecular_moment(
                self.config.acoustic.density_kg_m3(),
                self.config.acoustic.sound_speed_m_s(),
                self.config.reference_pressure_pa,
                p.density,
                p.pressure,
                self.grid.cell_volume(),
                kernel::rotate(self.config.grid_to_acoustic, p.velocity),
                stress,
                self.config.grid_to_acoustic,
            )
            .ok_or(CompressibleVolumeError3::Invalid(
                "unrepresentable full molecular stress moment",
            ))?;
        }
        poll(cancelled)?;
        std::mem::swap(&mut self.moments, &mut self.trial);
        self.time_s = Some(field.time_s());
        Ok(())
    }
}

pub(super) fn molecular_error(
    error: fs_flux::compressible::molecular::MolecularError,
) -> CompressibleVolumeError3 {
    if error == fs_flux::compressible::molecular::MolecularError::Cancelled {
        CompressibleVolumeError3::Cancelled
    } else {
        CompressibleVolumeError3::Molecular(error)
    }
}

pub(super) fn dissipation_error(error: DissipationError) -> CompressibleVolumeError3 {
    if error == DissipationError::Cancelled {
        CompressibleVolumeError3::Cancelled
    } else {
        CompressibleVolumeError3::Dissipation(error)
    }
}
