//! Native joint pressure/temperature water transport through the original mesh.
//!
//! Both fields share one native element choice and one barycentric query at
//! every stage. Two small pressure polynomials cache the current/previous
//! elements; no field is resampled, averaged or converted into temperature.
//! The original joint water owner supplies all EOS and optical derivatives.
//!
//! Shared nodal p and T make the phase index continuous, so the existing
//! face-coordinate bundle continuation applies without a Fresnel event. An
//! unresolved-face endpoint is admitted only when both original source
//! gradients agree. Pressure jumps represented by duplicated contact topology
//! still terminate at the original boundary rather than being smoothed.

use super::*;
use crate::charts::Ray;
use crate::material_state::thermo_optic::transient::mesh::{
    MAX_INTERNAL_FACES, MeshTemperature, P1ConductionMesh, P1MeshError, P1MeshOptions,
    P1MeshPropagation,
};
use crate::material_state::thermo_optic::transient::mesh::bundle::{
    P1MeshBundleOptions, P1MeshBundlePropagation,
};
use crate::material_state::water_refraction::{
    MaterialWaterRefractiveField, WaterRefractiveFieldError,
};
use core::cell::Cell;

impl<'a> P1ConductionMesh<'a> {
    /// Propagate a joint original P1 water state using native mesh boundaries.
    ///
    /// Construct the inner binding with `P1PressureTemperatureElement::bind_water`
    /// or the same derived joint identity. The outer optical identity may be
    /// distinct and is retained separately. The original temperature mesh and
    /// nodal buffer must be exactly those of this prepared traversal.
    #[allow(clippy::too_many_arguments)]
    pub fn propagate_pressure_temperature_water(
        &self,
        cx: &Cx<'_>,
        binding: &MaterialWaterRefractiveField<'_, PressureTemperatureWaterField<'_, P1PressureTemperatureElement<'a>>>,
        ray: Ray,
        transverse_axis: Vec3,
        wavelength_nm: f64,
        max_distance_m: f64,
        options: P1MeshOptions,
    ) -> Result<P1MeshPropagation, PressureTemperatureMeshError> {
        options.event.validate().map_err(P1MeshError::from)?;
        if options.max_internal_faces > MAX_INTERNAL_FACES {
            return Err(P1MeshError::Propagation(GradientIndexError::InvalidOptions).into());
        }
        let (cursor, element) = self.joint_water_cursor(binding)?;
        let joint = MeshPressureTemperature::new(&cursor, element);
        let original_state = binding.state_source();
        let state = PressureTemperatureWaterField::try_new(
            &joint, original_state.region(), original_state.snapshot_identity(),
        )?;
        let original = MaterialWaterRefractiveField::try_new(
            &state, binding.absorption_source(), binding.snapshot_identity(),
        )?;
        let prepared = original.prepare(wavelength_nm)?;
        self.propagate_prepared(
            cx, &cursor, &prepared, ray, transverse_axis, wavelength_nm, max_distance_m, options,
        ).map_err(PressureTemperatureMeshError::Mesh)
    }

    /// Carry the original launch-angle family through coupled p/T refraction.
    ///
    /// The full nonlinear mixed Hessian and native gradient-jump crossing-time
    /// corrections remain active. A pressure-only gradient jump has the same
    /// one-sided endpoint requirement as an original temperature-gradient jump.
    #[allow(clippy::too_many_arguments)]
    pub fn propagate_pressure_temperature_water_bundle(
        &self,
        cx: &Cx<'_>,
        binding: &MaterialWaterRefractiveField<'_, PressureTemperatureWaterField<'_, P1PressureTemperatureElement<'a>>>,
        ray: Ray,
        transverse_axis: Vec3,
        wavelength_nm: f64,
        max_distance_m: f64,
        options: P1MeshBundleOptions,
    ) -> Result<P1MeshBundlePropagation, PressureTemperatureMeshError> {
        options.mesh.event.validate().map_err(P1MeshError::from)?;
        if !max_distance_m.is_finite()
            || max_distance_m < 0.0
            || options.mesh.max_internal_faces > MAX_INTERNAL_FACES
        {
            return Err(P1MeshError::Propagation(GradientIndexError::InvalidOptions).into());
        }
        let (cursor, element) = self.joint_water_cursor(binding)?;
        let joint = MeshPressureTemperature::new(&cursor, element);
        let original_state = binding.state_source();
        let state = PressureTemperatureWaterField::try_new(
            &joint, original_state.region(), original_state.snapshot_identity(),
        )?;
        let original = MaterialWaterRefractiveField::try_new(
            &state, binding.absorption_source(), binding.snapshot_identity(),
        )?;
        let prepared = original.prepare(wavelength_nm)?;
        self.propagate_bundle_prepared(
            cx, &cursor, &prepared, ray, transverse_axis, wavelength_nm, max_distance_m, options,
            |previous, next| joint.pressure_gradient_jump(previous, next),
        ).map_err(PressureTemperatureMeshError::Mesh)
    }

    fn joint_water_cursor(
        &self,
        binding: &MaterialWaterRefractiveField<'_, PressureTemperatureWaterField<'_, P1PressureTemperatureElement<'a>>>,
    ) -> Result<(MeshTemperature<'a>, P1PressureTemperatureElement<'a>), PressureTemperatureMeshError> {
        let state = binding.state_source();
        let element = *state.state_source();
        if state.snapshot_identity() != element.snapshot().content_identity() {
            return Err(PressureTemperatureMeshError::SourceIdentityMismatch);
        }
        let temperature = element.temperature_element();
        let source = temperature.snapshot();
        if !core::ptr::eq(source.mesh(), self.snapshot.mesh())
            || !core::ptr::eq(source.nodal_temperature_k(), self.snapshot.nodal_temperature_k())
        {
            return Err(P1MeshError::ForeignSnapshot.into());
        }
        Ok((MeshTemperature {
            current: Cell::new(temperature),
            previous: Cell::new(None),
        }, element))
    }
}

// This cache is private per propagation and never chooses a geometric side.
// Original native barycentric admission supplies the one selected element.
// Only its pressure polynomial is looked up or prepared, with no allocation.
struct MeshPressureTemperature<'cursor, 'mesh> {
    cursor: &'cursor MeshTemperature<'mesh>,
    snapshot: P1PressureTemperatureSnapshot<'mesh>,
    current: Cell<(usize, PressurePolynomial)>,
    previous: Cell<Option<(usize, PressurePolynomial)>>,
}
impl<'cursor, 'mesh> MeshPressureTemperature<'cursor, 'mesh> {
    fn new(
        cursor: &'cursor MeshTemperature<'mesh>,
        element: P1PressureTemperatureElement<'mesh>,
    ) -> Self {
        Self {
            cursor,
            snapshot: element.snapshot(),
            current: Cell::new((element.temperature_element().element_index(), element.pressure)),
            previous: Cell::new(None),
        }
    }
    fn polynomial(&self, element: usize) -> Result<PressurePolynomial, GradientIndexError> {
        let current = self.current.get();
        if current.0 == element {
            return Ok(current.1);
        }
        if let Some(previous) = self.previous.get() {
            if previous.0 == element {
                self.current.set(previous);
                self.previous.set(Some(current));
                return Ok(previous.1);
            }
        }
        let pressure = self.snapshot.pressure_snapshot().polynomial(element)
            .map_err(PressureTemperatureSnapshotError::into_gradient_error)?;
        self.current.set((element, pressure));
        self.previous.set(Some(current));
        Ok(pressure)
    }
    fn sample_jet(&self, position_m: Point3) -> Result<PressureTemperatureJet, GradientIndexError> {
        let (temperature, lambda) = self.cursor.sample_barycentric(position_m)?;
        let pressure = self.polynomial(temperature.element_index())?;
        joint_jet(&temperature, pressure, lambda)
            .map_err(PressureTemperatureSnapshotError::into_gradient_error)
    }
    fn pressure_gradient_jump(
        &self,
        previous: P1ConductionElement<'mesh>,
        next: P1ConductionElement<'mesh>,
    ) -> Result<bool, GradientIndexError> {
        Ok(self.polynomial(previous.element_index())?.gradient_pa_per_m
            != self.polynomial(next.element_index())?.gradient_pa_per_m)
    }
}
impl PressureTemperatureGradientField for MeshPressureTemperature<'_, '_> {
    fn pressure_temperature_gradient(&self, p: Point3) -> Result<PressureTemperatureGradientSample, GradientIndexError> {
        self.sample_jet(p).map(|jet| jet.sample())
    }
}
impl PressureTemperatureSecondOrderField for MeshPressureTemperature<'_, '_> {
    fn pressure_temperature_jet(&self, p: Point3) -> Result<PressureTemperatureJet, GradientIndexError> {
        self.sample_jet(p)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum PressureTemperatureMeshError {
    Mesh(P1MeshError),
    State(PressureTemperatureWaterError),
    Optical(WaterRefractiveFieldError),
    SourceIdentityMismatch,
}
impl From<P1MeshError> for PressureTemperatureMeshError {
    fn from(error: P1MeshError) -> Self { Self::Mesh(error) }
}
impl From<PressureTemperatureWaterError> for PressureTemperatureMeshError {
    fn from(error: PressureTemperatureWaterError) -> Self { Self::State(error) }
}
impl From<WaterRefractiveFieldError> for PressureTemperatureMeshError {
    fn from(error: WaterRefractiveFieldError) -> Self { Self::Optical(error) }
}
impl core::fmt::Display for PressureTemperatureMeshError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "original P1 pressure/temperature water transport: {self:?}")
    }
}
impl core::error::Error for PressureTemperatureMeshError {}

#[cfg(test)]
mod tests;
