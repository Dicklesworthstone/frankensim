//! IF97/IAPWS water refraction through the original conduction mesh.
//!
//! One fixed absolute pressure and explicit IF97 branch apply to the same
//! accepted P1 temperature snapshot on every element. Shared nodal temperature
//! gives continuous density and phase index across an internal face, while
//! their gradients can jump. Those faces change ray curvature without an
//! artificial Fresnel event or direction jump.
//!
//! The existing mesh event and bundle owners retain the complete transition
//! collar, crossing-time variations, optical integrals and whole-path budgets.
//! The original EOS and optical chain rule supply the nonlinear index Hessian;
//! a zero in-element temperature Hessian does not imply a zero index Hessian.
//! Preparation is per wavelength; propagation adds no workspace allocation.
//!
//! This is a frozen single-phase water thermal field at fixed pressure. Native
//! domain admission refuses saturation/regional crossings, and duplicated
//! contact nodes or exterior faces terminate at their original boundary.
//! No fluid-flow closure, mixture interpolation or free-surface optics is added.

use super::*;
use super::bundle::{P1MeshBundleOptions, P1MeshBundlePropagation};
use crate::material_state::water_refraction::{
    MaterialWaterRefractiveField, WaterRefractiveFieldError,
};
use crate::material_state::water_refraction::isobaric::{
    IsobaricWaterError, IsobaricWaterField,
};

impl<'a> P1ConductionMesh<'a> {
    /// Propagate the original water binding through this exact thermal snapshot.
    ///
    /// The supplied element selects the launch cell. Both its mesh and nodal
    /// buffer must be the original snapshot borrowed by this prepared mesh.
    /// Pressure, selected IF97 branch, source and optical snapshot identities,
    /// and the independent absorption law are retained without substitution.
    #[allow(clippy::too_many_arguments)]
    pub fn propagate_water(
        &self,
        cx: &Cx<'_>,
        binding: &MaterialWaterRefractiveField<'_, IsobaricWaterField<'_, P1ConductionElement<'a>>>,
        ray: Ray,
        transverse_axis: Vec3,
        wavelength_nm: f64,
        max_distance_m: f64,
        options: P1MeshOptions,
    ) -> Result<P1MeshPropagation, WaterMeshError> {
        options.event.validate().map_err(P1MeshError::from)?;
        if options.max_internal_faces > MAX_INTERNAL_FACES {
            return Err(P1MeshError::Propagation(GradientIndexError::InvalidOptions).into());
        }
        let cursor = self.water_cursor(binding)?;
        let source = binding.state_source();
        let state = IsobaricWaterField::try_new(
            &cursor, source.pressure(), source.region(), source.snapshot_identity(),
        )?;
        // An outer optical binding can carry a distinct accepted identity.
        // Preserve it rather than replacing it with state.content_identity().
        let original = MaterialWaterRefractiveField::try_new(
            &state, binding.absorption_source(), binding.snapshot_identity(),
        )?;
        let prepared = original.prepare(wavelength_nm)?;
        self.propagate_prepared(
            cx, &cursor, &prepared, ray, transverse_axis, wavelength_nm, max_distance_m, options,
        ).map_err(WaterMeshError::Mesh)
    }

    /// Carry the original launch-angle family through the same water field.
    ///
    /// Full IF97 density curvature and the IAPWS temperature/density chain rule
    /// feed the native variational equations. Original gradient jumps use the
    /// existing continuous-index face-coordinate continuation, including their
    /// crossing-time correction; receiving-plane projection stays with the
    /// returned native bundle.
    #[allow(clippy::too_many_arguments)]
    pub fn propagate_water_bundle(
        &self,
        cx: &Cx<'_>,
        binding: &MaterialWaterRefractiveField<'_, IsobaricWaterField<'_, P1ConductionElement<'a>>>,
        ray: Ray,
        transverse_axis: Vec3,
        wavelength_nm: f64,
        max_distance_m: f64,
        options: P1MeshBundleOptions,
    ) -> Result<P1MeshBundlePropagation, WaterMeshError> {
        options.mesh.event.validate().map_err(P1MeshError::from)?;
        if !max_distance_m.is_finite()
            || max_distance_m < 0.0
            || options.mesh.max_internal_faces > MAX_INTERNAL_FACES
        {
            return Err(P1MeshError::Propagation(GradientIndexError::InvalidOptions).into());
        }
        let cursor = self.water_cursor(binding)?;
        let source = binding.state_source();
        let state = IsobaricWaterField::try_new(
            &cursor, source.pressure(), source.region(), source.snapshot_identity(),
        )?;
        let original = MaterialWaterRefractiveField::try_new(
            &state, binding.absorption_source(), binding.snapshot_identity(),
        )?;
        let prepared = original.prepare(wavelength_nm)?;
        self.propagate_bundle_prepared(
            cx, &cursor, &prepared, ray, transverse_axis, wavelength_nm, max_distance_m, options,
        ).map_err(WaterMeshError::Mesh)
    }

    fn water_cursor(
        &self,
        binding: &MaterialWaterRefractiveField<'_, IsobaricWaterField<'_, P1ConductionElement<'a>>>,
    ) -> Result<MeshTemperature<'a>, WaterMeshError> {
        let start = *binding.state_source().temperature_source();
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

/// Native mesh/transport failures and original water-source preparation failures.
#[derive(Clone, Debug, PartialEq)]
pub enum WaterMeshError {
    Mesh(P1MeshError),
    State(IsobaricWaterError),
    Optical(WaterRefractiveFieldError),
}
impl From<P1MeshError> for WaterMeshError {
    fn from(error: P1MeshError) -> Self {
        Self::Mesh(error)
    }
}
impl From<IsobaricWaterError> for WaterMeshError {
    fn from(error: IsobaricWaterError) -> Self {
        Self::State(error)
    }
}
impl From<WaterRefractiveFieldError> for WaterMeshError {
    fn from(error: WaterRefractiveFieldError) -> Self {
        Self::Optical(error)
    }
}
impl core::fmt::Display for WaterMeshError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "transient P1 water refraction: {self:?}")
    }
}
impl core::error::Error for WaterMeshError {}

#[cfg(test)]
mod tests;
