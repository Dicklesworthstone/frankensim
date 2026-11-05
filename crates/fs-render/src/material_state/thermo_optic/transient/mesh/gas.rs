//! Molecular gas refraction through the original continuous P1 mesh.
//!
//! The native mesh owns element selection, face events and launch-angle
//! continuation. The original gas owner supplies the local EOS, Lorentz--Lorenz
//! index and exact density Hessian. Shared nodal temperature and, when supplied,
//! pressure give continuous index while allowing unsmoothed gradient jumps.
//! No Fresnel interface or instantaneous direction jump is invented there.
//!
//! Original field buffers and identities remain attached. The joint p/T path
//! reuses the same barycentric sample and two-element pressure cache as water
//! transport. Duplicated contact topology and exterior faces retain native
//! boundary termination; no extrapolation or replacement pressure solve occurs.

use super::bundle::{P1MeshBundleOptions, P1MeshBundlePropagation};
use super::water::pressure_temperature::P1PressureTemperatureElement;
use super::water::pressure_temperature::transport::MeshPressureTemperature;
use super::*;
use crate::material_state::gas_scattering::GasRefractiveFieldError;
use crate::material_state::gas_scattering::isobaric::{IsobaricGasError, IsobaricGasField};
use crate::material_state::gas_scattering::pressure_temperature::{
    PressureTemperatureGasError, PressureTemperatureGasField,
};
use fs_qty::{Pressure, Temperature};

impl<'a> P1ConductionMesh<'a> {
    /// Propagate an original isobaric gas through its exact thermal mesh.
    /// Reference temperature prepares the original molecular coefficients;
    /// every stage still samples the actual native temperature and EOS.
    #[allow(clippy::too_many_arguments)]
    pub fn propagate_gas(
        &self,
        cx: &Cx<'_>,
        state: &IsobaricGasField<'_, P1ConductionElement<'a>>,
        ray: Ray,
        transverse_axis: Vec3,
        wavelength_nm: f64,
        max_distance_m: f64,
        options: P1MeshOptions,
    ) -> Result<P1MeshPropagation, GasMeshError> {
        admit_gas_mesh_options(options, max_distance_m)?;
        let cursor = self.gas_cursor(*state.temperature_source())?;
        let reference = state.reference_binding().state().gas();
        let original = IsobaricGasField::try_new(
            &cursor,
            state.model(),
            state.pressure(),
            Temperature::new(reference.temperature),
            state.snapshot_identity(),
        )?;
        let optics = original.bind_refraction()?;
        self.propagate_prepared(
            cx,
            &cursor,
            &optics,
            ray,
            transverse_axis,
            wavelength_nm,
            max_distance_m,
            options,
        )
        .map_err(GasMeshError::Mesh)
    }

    /// Carry the native launch-angle family through thermal gas focusing.
    /// Density curvature and the nonlinear optical Hessian stay active even
    /// within a P1 element whose temperature Hessian is zero.
    #[allow(clippy::too_many_arguments)]
    pub fn propagate_gas_bundle(
        &self,
        cx: &Cx<'_>,
        state: &IsobaricGasField<'_, P1ConductionElement<'a>>,
        ray: Ray,
        transverse_axis: Vec3,
        wavelength_nm: f64,
        max_distance_m: f64,
        options: P1MeshBundleOptions,
    ) -> Result<P1MeshBundlePropagation, GasMeshError> {
        admit_gas_mesh_options(options.mesh, max_distance_m)?;
        let cursor = self.gas_cursor(*state.temperature_source())?;
        let reference = state.reference_binding().state().gas();
        let original = IsobaricGasField::try_new(
            &cursor,
            state.model(),
            state.pressure(),
            Temperature::new(reference.temperature),
            state.snapshot_identity(),
        )?;
        let optics = original.bind_refraction()?;
        self.propagate_bundle_prepared(
            cx,
            &cursor,
            &optics,
            ray,
            transverse_axis,
            wavelength_nm,
            max_distance_m,
            options,
            |_, _| Ok(false),
        )
        .map_err(GasMeshError::Mesh)
    }

    /// Propagate the original joint absolute-pressure/temperature gas state.
    /// The state must identify the native joint snapshot's content identity.
    /// No acoustic perturbation is implicitly treated as absolute pressure.
    #[allow(clippy::too_many_arguments)]
    pub fn propagate_pressure_temperature_gas(
        &self,
        cx: &Cx<'_>,
        state: &PressureTemperatureGasField<'_, P1PressureTemperatureElement<'a>>,
        ray: Ray,
        transverse_axis: Vec3,
        wavelength_nm: f64,
        max_distance_m: f64,
        options: P1MeshOptions,
    ) -> Result<P1MeshPropagation, GasMeshError> {
        admit_gas_mesh_options(options, max_distance_m)?;
        let (cursor, element) = self.joint_gas_cursor(state)?;
        let joint = MeshPressureTemperature::new(&cursor, element);
        let reference = state.reference_binding().state().gas();
        let original = PressureTemperatureGasField::try_new(
            &joint,
            state.model(),
            Temperature::new(reference.temperature),
            Pressure::new(reference.pressure),
            state.snapshot_identity(),
        )?;
        let optics = original.bind_refraction()?;
        self.propagate_prepared(
            cx,
            &cursor,
            &optics,
            ray,
            transverse_axis,
            wavelength_nm,
            max_distance_m,
            options,
        )
        .map_err(GasMeshError::Mesh)
    }

    /// Propagate rays and their full angular derivatives across native p/T
    /// gradient jumps. Pressure-only jumps retain the same crossing-time
    /// correction and one-sided endpoint requirement as temperature jumps.
    #[allow(clippy::too_many_arguments)]
    pub fn propagate_pressure_temperature_gas_bundle(
        &self,
        cx: &Cx<'_>,
        state: &PressureTemperatureGasField<'_, P1PressureTemperatureElement<'a>>,
        ray: Ray,
        transverse_axis: Vec3,
        wavelength_nm: f64,
        max_distance_m: f64,
        options: P1MeshBundleOptions,
    ) -> Result<P1MeshBundlePropagation, GasMeshError> {
        admit_gas_mesh_options(options.mesh, max_distance_m)?;
        let (cursor, element) = self.joint_gas_cursor(state)?;
        let joint = MeshPressureTemperature::new(&cursor, element);
        let reference = state.reference_binding().state().gas();
        let original = PressureTemperatureGasField::try_new(
            &joint,
            state.model(),
            Temperature::new(reference.temperature),
            Pressure::new(reference.pressure),
            state.snapshot_identity(),
        )?;
        let optics = original.bind_refraction()?;
        self.propagate_bundle_prepared(
            cx,
            &cursor,
            &optics,
            ray,
            transverse_axis,
            wavelength_nm,
            max_distance_m,
            options,
            |previous, next| joint.pressure_gradient_jump(previous, next),
        )
        .map_err(GasMeshError::Mesh)
    }

    fn gas_cursor(
        &self,
        element: P1ConductionElement<'a>,
    ) -> Result<MeshTemperature<'a>, GasMeshError> {
        let source = element.snapshot();
        if !core::ptr::eq(source.mesh(), self.snapshot.mesh())
            || !core::ptr::eq(
                source.nodal_temperature_k(),
                self.snapshot.nodal_temperature_k(),
            )
        {
            return Err(GasMeshError::Mesh(P1MeshError::ForeignSnapshot));
        }
        Ok(MeshTemperature {
            current: Cell::new(element),
            previous: Cell::new(None),
        })
    }

    fn joint_gas_cursor(
        &self,
        state: &PressureTemperatureGasField<'_, P1PressureTemperatureElement<'a>>,
    ) -> Result<(MeshTemperature<'a>, P1PressureTemperatureElement<'a>), GasMeshError> {
        let element = *state.state_source();
        if state.snapshot_identity() != element.snapshot().content_identity() {
            return Err(GasMeshError::SourceIdentityMismatch);
        }
        Ok((self.gas_cursor(element.temperature_element())?, element))
    }
}

fn admit_gas_mesh_options(options: P1MeshOptions, distance_m: f64) -> Result<(), GasMeshError> {
    options.event.validate().map_err(P1MeshError::from)?;
    if !distance_m.is_finite()
        || distance_m < 0.0
        || options.max_internal_faces > MAX_INTERNAL_FACES
    {
        return Err(P1MeshError::Propagation(GradientIndexError::InvalidOptions).into());
    }
    Ok(())
}

#[derive(Clone, Debug, PartialEq)]
pub enum GasMeshError {
    Mesh(P1MeshError),
    Isobaric(IsobaricGasError),
    PressureTemperature(PressureTemperatureGasError),
    Optical(GasRefractiveFieldError),
    SourceIdentityMismatch,
}
impl From<P1MeshError> for GasMeshError {
    fn from(error: P1MeshError) -> Self {
        Self::Mesh(error)
    }
}
impl From<IsobaricGasError> for GasMeshError {
    fn from(error: IsobaricGasError) -> Self {
        Self::Isobaric(error)
    }
}
impl From<PressureTemperatureGasError> for GasMeshError {
    fn from(error: PressureTemperatureGasError) -> Self {
        Self::PressureTemperature(error)
    }
}
impl From<GasRefractiveFieldError> for GasMeshError {
    fn from(error: GasRefractiveFieldError) -> Self {
        Self::Optical(error)
    }
}
impl core::fmt::Display for GasMeshError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "original P1 gas refraction transport: {self:?}")
    }
}
impl core::error::Error for GasMeshError {}

#[cfg(test)]
mod tests;
