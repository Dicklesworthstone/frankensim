//! Coherent original molecular surface/volume sampling from one stress field.
use super::*;
use super::super::super::molecular::{dissipation_error, molecular_error};
use fs_flux::compressible::spatial::{mixture::molecular as mixture, molecular as frozen};
use mixture::{DissipationControl, TransportLaw};

impl ExteriorSurfaceWelchStream {
    /// Record the current native molecular endpoint once, retaining both
    /// -tau.n on the surface and -tau in the disjoint exterior volume.
    /// One source-owned stress evaluation supplies both accepted samplers.
    #[allow(clippy::too_many_arguments)]
    pub fn from_molecular_frozen<L: TransportLaw>(
        flow: &mut frozen::MolecularEuler3d<L>,
        control_box: PermeableBox3,
        exterior_cells: &[usize],
        source_config: CompressibleVolumeConfig3,
        clock: SourceSampleClock,
        expected_samples: usize,
        config: WelchConfig,
        budget: ExteriorSurfaceRecordingBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ExteriorSurfaceRecordingError> {
        poll(&mut cancelled)?;
        clock.admit(0, flow.flow().time().value())?;
        admit_exterior(
            flow.flow().grid(),
            control_box,
            exterior_cells,
            source_config,
            budget.recording.maximum_compact_sources,
            budget.recording.maximum_quadrupole_sources,
            budget.maximum_volume_workspace_bytes,
            &mut cancelled,
        )?;
        // Reuse the original complete ambient/source/geometry admission.
        // These private provisional samples never enter the recording.
        let mut surface = PreparedPermeableSurface3::from_frozen(
            flow.flow(), control_box, source_config, budget.surface, &mut cancelled,
        )?;
        let mut volume = PreparedCompressibleVolume3::from_frozen(
            flow.flow(), exterior_cells, source_config, budget.maximum_volume_workspace_bytes,
            &mut cancelled,
        )?;
        surface.with_viscous_stress = true;
        surface.time_s = None;
        volume.with_viscous_stress = true;
        volume.time_s = None;
        let field = flow.stress_field(&mut cancelled).map_err(molecular_error)?;
        surface.sample_stress(&field, &mut cancelled)?;
        volume.sample_stress(&field, &mut cancelled)?;
        Self::from_sources(
            surface, volume, clock, expected_samples, config, budget.recording, cancelled,
        )
    }

    /// The same one-evaluation path with the actual evolving species source
    /// and explicit ambient composition, without a frozen-composition surrogate.
    #[allow(clippy::too_many_arguments)]
    pub fn from_molecular_mixture<L: TransportLaw>(
        flow: &mut mixture::MolecularEuler3d<L>,
        reference_partial_densities_kg_m3: &[f64],
        control_box: PermeableBox3,
        exterior_cells: &[usize],
        source_config: CompressibleVolumeConfig3,
        clock: SourceSampleClock,
        expected_samples: usize,
        config: WelchConfig,
        budget: ExteriorSurfaceRecordingBudget,
        controls: DissipationControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ExteriorSurfaceRecordingError> {
        poll(&mut cancelled)?;
        clock.admit(0, flow.flow().time().value())?;
        admit_exterior(
            flow.flow().grid(),
            control_box,
            exterior_cells,
            source_config,
            budget.recording.maximum_compact_sources,
            budget.recording.maximum_quadrupole_sources,
            budget.maximum_volume_workspace_bytes,
            &mut cancelled,
        )?;
        let mut surface = PreparedPermeableSurface3::from_mixture(
            flow.flow(), reference_partial_densities_kg_m3, control_box,
            source_config, budget.surface, &mut cancelled,
        )?;
        let mut volume = PreparedCompressibleVolume3::from_mixture(
            flow.flow(), exterior_cells, source_config, budget.maximum_volume_workspace_bytes,
            &mut cancelled,
        )?;
        surface.with_viscous_stress = true;
        surface.time_s = None;
        volume.with_viscous_stress = true;
        volume.time_s = None;
        let field = flow.stress_field(controls, &mut cancelled).map_err(dissipation_error)?;
        surface.sample_stress(&field, &mut cancelled)?;
        volume.sample_stress(&field, &mut cancelled)?;
        Self::from_sources(
            surface, volume, clock, expected_samples, config, budget.recording, cancelled,
        )
    }

    /// Acquire surface and exterior stress at the same accepted native time,
    /// using one full-grid constitutive/gradient evaluation. A failed volume
    /// extraction or joint transform restores BOTH preceding source samples.
    pub fn sample_molecular_frozen<L: TransportLaw>(
        &mut self,
        flow: &mut frozen::MolecularEuler3d<L>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), ExteriorSurfaceRecordingError> {
        self.admit_time(flow.flow().time().value(), &mut cancelled)?;
        if !self.surface.with_viscous_stress
            || !self.volume.with_viscous_stress
            || self.surface.frozen_source() != Some(flow.flow().gas())
            || self.volume.frozen_source() != Some(flow.flow().gas())
        {
            return Err(CompressibleVolumeError3::SourceMismatch.into());
        }
        self.surface.admit_sample(
            flow.flow().grid(), flow.flow().boundaries(), flow.flow().time().value(),
        )?;
        self.volume.admit_sample(
            flow.flow().grid(), flow.flow().boundaries(), flow.flow().time().value(),
        )?;
        let field = flow.stress_field(&mut cancelled).map_err(molecular_error)?;
        self.sample_stress(&field, &mut cancelled)
    }

    pub fn sample_molecular_mixture<L: TransportLaw>(
        &mut self,
        flow: &mut mixture::MolecularEuler3d<L>,
        controls: DissipationControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), ExteriorSurfaceRecordingError> {
        self.admit_time(flow.flow().time().value(), &mut cancelled)?;
        if !self.surface.with_viscous_stress
            || !self.volume.with_viscous_stress
            || self.surface.mixture_source() != Some(flow.flow().gas())
            || self.volume.mixture_source() != Some(flow.flow().gas())
        {
            return Err(CompressibleVolumeError3::SourceMismatch.into());
        }
        self.surface.admit_sample(
            flow.flow().grid(), flow.flow().boundaries(), flow.flow().time().value(),
        )?;
        self.volume.admit_sample(
            flow.flow().grid(), flow.flow().boundaries(), flow.flow().time().value(),
        )?;
        let field = flow.stress_field(controls, &mut cancelled).map_err(dissipation_error)?;
        self.sample_stress(&field, &mut cancelled)
    }

    fn sample_stress(
        &mut self,
        field: &frozen::StressField3d<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ExteriorSurfaceRecordingError> {
        let previous = (self.surface.time_s, self.volume.time_s);
        self.surface.sample_stress(field, cancelled)?;
        let volume_result = self.volume.sample_stress(field, cancelled);
        self.finish_sample(previous, volume_result, cancelled)
    }
}
