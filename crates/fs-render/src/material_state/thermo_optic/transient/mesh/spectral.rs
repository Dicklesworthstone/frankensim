//! Original tabulated thermal optics on an accepted P1 conduction mesh.
//!
//! The original spectral/temperature source remains borrowed while the native
//! cursor selects each one-sided P1 reconstruction. A wavelength is prepared
//! once per propagation. The existing mesh walker owns geometric admission,
//! transition collars, path budgets and the bundle's gradient-jump correction.
//! This adapter introduces no solver, temperature smoothing or heap workspace.

use super::bundle::{P1MeshBundleOptions, P1MeshBundlePropagation};
use super::*;
use crate::material_state::spectral_refraction::{
    MaterialSpectralRefractiveField, SpectralRefractiveFieldError,
};

impl<'a> P1ConductionMesh<'a> {
    /// Follow the original phase, group and absorption law across this exact
    /// accepted mesh/nodal snapshot. The binding's original element selects
    /// the launch cell; matching nodal values in a different slice do not
    /// authorize using another physical snapshot.
    #[allow(clippy::too_many_arguments)]
    pub fn propagate_spectral(
        &self,
        cx: &Cx<'_>,
        binding: &MaterialSpectralRefractiveField<'_, P1ConductionElement<'a>>,
        ray: Ray,
        transverse_axis: Vec3,
        wavelength_nm: f64,
        max_distance_m: f64,
        options: P1MeshOptions,
    ) -> Result<P1MeshPropagation, SpectralMeshError> {
        options.event.validate().map_err(P1MeshError::from)?;
        if options.max_internal_faces > MAX_INTERNAL_FACES {
            return Err(P1MeshError::Propagation(GradientIndexError::InvalidOptions).into());
        }
        let cursor = self.spectral_cursor(binding)?;
        let original = binding.with_temperature_source(&cursor)?;
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
        .map_err(SpectralMeshError::Mesh)
    }

    /// Carry native angular variations with the source's complete analytic
    /// thermal chain rule. P1 temperature Hessians vanish inside each element;
    /// nonlinear source curvature and face-crossing saltation do not.
    #[allow(clippy::too_many_arguments)]
    pub fn propagate_spectral_bundle(
        &self,
        cx: &Cx<'_>,
        binding: &MaterialSpectralRefractiveField<'_, P1ConductionElement<'a>>,
        ray: Ray,
        transverse_axis: Vec3,
        wavelength_nm: f64,
        max_distance_m: f64,
        options: P1MeshBundleOptions,
    ) -> Result<P1MeshBundlePropagation, SpectralMeshError> {
        options.mesh.event.validate().map_err(P1MeshError::from)?;
        if !max_distance_m.is_finite()
            || max_distance_m < 0.0
            || options.mesh.max_internal_faces > MAX_INTERNAL_FACES
        {
            return Err(P1MeshError::Propagation(GradientIndexError::InvalidOptions).into());
        }
        let cursor = self.spectral_cursor(binding)?;
        let original = binding.with_temperature_source(&cursor)?;
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
        .map_err(SpectralMeshError::Mesh)
    }

    fn spectral_cursor(
        &self,
        binding: &MaterialSpectralRefractiveField<'_, P1ConductionElement<'a>>,
    ) -> Result<MeshTemperature<'a>, SpectralMeshError> {
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

/// Original mesh/transport failures and original spectral-source refusals.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpectralMeshError {
    Mesh(P1MeshError),
    Optical(SpectralRefractiveFieldError),
}
impl From<P1MeshError> for SpectralMeshError {
    fn from(error: P1MeshError) -> Self {
        Self::Mesh(error)
    }
}
impl From<SpectralRefractiveFieldError> for SpectralMeshError {
    fn from(error: SpectralRefractiveFieldError) -> Self {
        Self::Optical(error)
    }
}
impl core::fmt::Display for SpectralMeshError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "transient P1 tabulated spectral refraction: {self:?}")
    }
}
impl core::error::Error for SpectralMeshError {}
