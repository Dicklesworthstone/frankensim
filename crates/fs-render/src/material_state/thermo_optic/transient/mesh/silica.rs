//! Sourced silica optics on the original accepted transient conduction mesh.
//!
//! The NASA CHARMS law reads each original P1 temperature reconstruction.
//! Shared nodal temperatures preserve the phase index across internal faces;
//! unsmoothed gradient jumps alter curvature without inventing Fresnel events.
//! The native mesh walker retains face-crossing variations, complete phase/group
//! lengths and Beer attenuation, including its two-element transition collars.
//!
//! Spatial preparation allocates once in the original mesh owner. Wavelength
//! preparation and subsequent propagation add no heap workspace. Absorption,
//! the accepted snapshot identity and source temperature/wavelength bounds are
//! retained. Exterior faces and duplicated contact nodes stop at their native
//! boundaries; this adapter does not infer stress optics or thermal expansion.

use super::bundle::{P1MeshBundleOptions, P1MeshBundlePropagation};
use super::*;
use crate::material_state::silica_refraction::{MaterialSilicaRefractiveField, SilicaFieldError};

impl<'a> P1ConductionMesh<'a> {
    /// Follow the original silica thermal binding across this exact mesh.
    /// The binding's element selects the launch cell; its mesh and nodal slice
    /// must be the original accepted snapshot, including pointer identity.
    #[allow(clippy::too_many_arguments)]
    pub fn propagate_silica(
        &self,
        cx: &Cx<'_>,
        binding: &MaterialSilicaRefractiveField<'_, P1ConductionElement<'a>>,
        ray: Ray,
        transverse_axis: Vec3,
        wavelength_nm: f64,
        max_distance_m: f64,
        options: P1MeshOptions,
    ) -> Result<P1MeshPropagation, SilicaMeshError> {
        options.event.validate().map_err(P1MeshError::from)?;
        if options.max_internal_faces > MAX_INTERNAL_FACES {
            return Err(P1MeshError::Propagation(GradientIndexError::InvalidOptions).into());
        }
        let cursor = self.silica_cursor(binding)?;
        let original = MaterialSilicaRefractiveField::try_new(
            &cursor,
            binding.absorption_source(),
            binding.snapshot_identity(),
        )?;
        let prepared = original.prepare(wavelength_nm)?;
        self.propagate_prepared(
            cx,
            &cursor,
            &prepared,
            ray,
            transverse_axis,
            wavelength_nm,
            max_distance_m,
            options,
        )
        .map_err(SilicaMeshError::Mesh)
    }

    /// Carry the same launch-angle family across unsmoothed thermal faces.
    /// The source's nonlinear temperature curvature remains in the index
    /// Hessian even though each P1 element has zero temperature Hessian.
    #[allow(clippy::too_many_arguments)]
    pub fn propagate_silica_bundle(
        &self,
        cx: &Cx<'_>,
        binding: &MaterialSilicaRefractiveField<'_, P1ConductionElement<'a>>,
        ray: Ray,
        transverse_axis: Vec3,
        wavelength_nm: f64,
        max_distance_m: f64,
        options: P1MeshBundleOptions,
    ) -> Result<P1MeshBundlePropagation, SilicaMeshError> {
        options.mesh.event.validate().map_err(P1MeshError::from)?;
        if !max_distance_m.is_finite()
            || max_distance_m < 0.0
            || options.mesh.max_internal_faces > MAX_INTERNAL_FACES
        {
            return Err(P1MeshError::Propagation(GradientIndexError::InvalidOptions).into());
        }
        let cursor = self.silica_cursor(binding)?;
        let original = MaterialSilicaRefractiveField::try_new(
            &cursor,
            binding.absorption_source(),
            binding.snapshot_identity(),
        )?;
        let prepared = original.prepare(wavelength_nm)?;
        self.propagate_bundle_prepared(
            cx,
            &cursor,
            &prepared,
            ray,
            transverse_axis,
            wavelength_nm,
            max_distance_m,
            options,
            |_, _| Ok(false),
        )
        .map_err(SilicaMeshError::Mesh)
    }

    fn silica_cursor(
        &self,
        binding: &MaterialSilicaRefractiveField<'_, P1ConductionElement<'a>>,
    ) -> Result<MeshTemperature<'a>, SilicaMeshError> {
        let start = *binding.temperature_source();
        let source = start.snapshot();
        if !core::ptr::eq(source.mesh(), self.snapshot.mesh())
            || !core::ptr::eq(
                source.nodal_temperature_k(),
                self.snapshot.nodal_temperature_k(),
            )
        {
            return Err(P1MeshError::ForeignSnapshot.into());
        }
        Ok(MeshTemperature {
            current: Cell::new(start),
            previous: Cell::new(None),
        })
    }
}

/// Original geometric/transport failures and original silica-source refusals.
#[derive(Clone, Debug, PartialEq)]
pub enum SilicaMeshError {
    Mesh(P1MeshError),
    Optical(SilicaFieldError),
}
impl From<P1MeshError> for SilicaMeshError {
    fn from(error: P1MeshError) -> Self {
        Self::Mesh(error)
    }
}
impl From<SilicaFieldError> for SilicaMeshError {
    fn from(error: SilicaFieldError) -> Self {
        Self::Optical(error)
    }
}
impl core::fmt::Display for SilicaMeshError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "transient P1 silica refraction: {self:?}")
    }
}
impl core::error::Error for SilicaMeshError {}

#[cfg(test)]
mod tests;
