//! Original full-stress permeable-surface source extraction.
use super::*;
use super::super::molecular::{dissipation_error, molecular_error};
use fs_flux::compressible::spatial::{mixture::molecular as mixture, molecular as frozen};
use mixture::{DissipationControl, TransportLaw};

impl PreparedPermeableSurface3 {
    /// Prepare the same source/ambient/closed-box admission as Euler sampling,
    /// retaining the original viscous/SGS face traction for every observation.
    pub fn from_molecular_frozen<L: TransportLaw>(
        flow: &mut frozen::MolecularEuler3d<L>,
        control_box: PermeableBox3,
        config: CompressibleVolumeConfig3,
        controls: PermeableSurfaceControls3,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CompressibleVolumeError3> {
        let mut result = Self::from_frozen(flow.flow(), control_box, config, controls, &mut cancelled)?;
        result.with_viscous_stress = true;
        result.time_s = None;
        result.sample_molecular_frozen(flow, &mut cancelled)?;
        Ok(result)
    }

    /// Retain the complete original mixture, explicit ambient composition and
    /// actual configured molecular/SGS traction, including bulk viscosity.
    pub fn from_molecular_mixture<L: TransportLaw>(
        flow: &mut mixture::MolecularEuler3d<L>,
        reference_partial_densities_kg_m3: &[f64],
        control_box: PermeableBox3,
        config: CompressibleVolumeConfig3,
        controls: PermeableSurfaceControls3,
        transport_controls: DissipationControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CompressibleVolumeError3> {
        let mut result = Self::from_mixture(
            flow.flow(),
            reference_partial_densities_kg_m3,
            control_box,
            config,
            controls,
            &mut cancelled,
        )?;
        result.with_viscous_stress = true;
        result.time_s = None;
        result.sample_molecular_mixture(flow, transport_controls, &mut cancelled)?;
        Ok(result)
    }

    /// Publish a later accepted endpoint without replacing the original law,
    /// face stress or physical clock. Native reusable constitutive caches may
    /// change; accepted flow, closure history and interval loads do not.
    pub fn sample_molecular_frozen<L: TransportLaw>(
        &mut self,
        flow: &mut frozen::MolecularEuler3d<L>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), CompressibleVolumeError3> {
        poll(&mut cancelled)?;
        if !self.with_viscous_stress || self.frozen_source() != Some(flow.flow().gas()) {
            return Err(CompressibleVolumeError3::SourceMismatch);
        }
        self.admit_sample(flow.flow().grid(), flow.flow().boundaries(), flow.flow().time().value())?;
        let field = flow.stress_field(&mut cancelled).map_err(molecular_error)?;
        self.sample_stress(&field, &mut cancelled)
    }

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
        self.admit_sample(flow.flow().grid(), flow.flow().boundaries(), flow.flow().time().value())?;
        let field = flow.stress_field(controls, &mut cancelled).map_err(dissipation_error)?;
        self.sample_stress(&field, &mut cancelled)
    }

    /// The joint exterior owner checks gas identity before borrowing one field
    /// for both source mechanisms. This method stages the whole face sample.
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
        let [nx, ny, _] = self.grid.shape();
        for (i, face) in self.faces.iter().enumerate() {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            let state = |cell| -> Result<(Primitive, [f64; 3]), CompressibleVolumeError3> {
                Ok((
                    field.primitive(cell).ok_or(CompressibleVolumeError3::GeometryMismatch)?,
                    field.momentum(cell).ok_or(CompressibleVolumeError3::GeometryMismatch)?,
                ))
            };
            let states = [state(face.adjacent_cells[0])?, state(face.adjacent_cells[1])?];
            // The positive-axis adjacent cell has the same native coordinates
            // as this intervening face: no midpoint or ghost index is inferred.
            let right = face.adjacent_cells[1];
            let coordinates = [right % nx, (right / nx) % ny, right / nx / ny];
            let traction = field.face_traction_pa(face.side.axis(), coordinates)
                .ok_or(CompressibleVolumeError3::GeometryMismatch)?;
            self.trial[i] = surface_flux(face, self.config, states, Some(traction))?;
        }
        poll(cancelled)?;
        std::mem::swap(&mut self.sources, &mut self.trial);
        self.time_s = Some(field.time_s());
        Ok(())
    }
}
