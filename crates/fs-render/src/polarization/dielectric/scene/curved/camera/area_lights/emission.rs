//! Original source laws shared by area proposals and complementary hits.
//!
//! A thermal source retains the complete native snapshot. It is never fitted
//! into a calibrated spectrum; Planck and emissivity are evaluated only at
//! the actual path wavelength. Camera-owned thermal snapshots are values,
//! while calibrated sources retain their original borrowed identity.

use super::super::{CoatingThermalEmitter, DiffuseThermalEmitter, phase_index};
use super::{
    AreaEmitterError, CameraRadianceError, CurvedCameraScene, EmissionSides, SurfaceEmitter,
};
use crate::polarization::dielectric::scene::curved::CurvedPolarizedScene;
use crate::instances::GeometryInstance;
use crate::polarization::Stokes;
use crate::gradient_index::event::BoundarySide;
use fs_geom::{Point3, Vec3};

#[derive(Clone, Debug)]
pub(super) enum AreaEmission<'e> {
    Calibrated(&'e SurfaceEmitter),
    Thermal(DiffuseThermalEmitter),
    Coating(CoatingThermalEmitter),
}

impl<'e> AreaEmission<'e> {
    pub(super) const fn object_id(&self) -> u64 {
        match self {
            Self::Calibrated(emitter) => emitter.object_id(),
            Self::Thermal(emitter) => emitter.object_id(),
            Self::Coating(emitter) => emitter.object_id(),
        }
    }

    pub(super) fn supports_side(&self, side: BoundarySide) -> bool {
        if side == BoundarySide::OnSurface { return false; }
        match self {
            Self::Calibrated(emitter) => emitter.sides() == EmissionSides::Both
                || side == BoundarySide::Positive,
            Self::Thermal(_) => true,
            Self::Coating(emitter) => emitter.side() == side,
        }
    }

    /// A direction-independent source can be evaluated before image discovery.
    /// A coating retains its actual finite band here, but its radiance and
    /// polarization require the converged outgoing direction of each image.
    pub(super) fn uniform_radiance(
        &self, wavelength_nm: f64,
    ) -> Result<Option<f64>, CameraRadianceError> {
        match self {
            Self::Calibrated(emitter) => emitter.spectrum().evaluate(wavelength_nm).map(Some),
            Self::Thermal(emitter) => emitter
                .surface()
                .radiance_w_m2_sr_nm(wavelength_nm)
                .map(Some)
                .map_err(CameraRadianceError::Thermal),
            Self::Coating(emitter) => {
                let (low, high) = emitter.source().wavelength_bounds_nm();
                if !wavelength_nm.is_finite() || !(low..=high).contains(&wavelength_nm) {
                    return Err(CameraRadianceError::CoatingThermal(
                        crate::thermal_surface::ThermalSurfaceError::Wavelength.into()));
                }
                Ok(None)
            }
        }
    }

    /// The declared body axis owns every area-bundle frame, including a
    /// spectral knot or interpolated wavelength whose Stokes state is
    /// unpolarized. Preparation and shooting therefore use the same chart
    /// without a wavelength-dependent switch to an unrelated caller axis.
    /// Sources with no declaration retain their original arbitrary basis.
    pub(super) fn polarization_axis(
        &self,
        instance: &GeometryInstance,
        unpolarized_axis: Vec3,
    ) -> Vec3 {
        if let Self::Calibrated(emitter) = self
            && let Some(polarization) = emitter.polarization()
        {
            instance.transform().transform_vector(polarization.reference_axis())
        } else {
            unpolarized_axis
        }
    }

    /// The same original spectral source law as camera hits, expressed in the
    /// source-owned area-bundle frame admitted during preparation.
    pub(super) fn polarization(
        &self,
        instance: &GeometryInstance,
        wavelength_nm: f64,
        unpolarized_axis: Vec3,
    ) -> Result<(Vec3, Stokes), CameraRadianceError> {
        if let Self::Calibrated(emitter) = self
            && let Some(polarization) = emitter.polarization()
        {
            let state = polarization.normalized_stokes(wavelength_nm)?;
            return Ok((self.polarization_axis(instance, unpolarized_axis), state));
        }
        Ok((unpolarized_axis, Stokes::unpolarized(1.0)))
    }

    pub(super) fn admit_local_medium(
        &self,
        scene: &CurvedPolarizedScene<'_>,
        position: Point3,
        medium: usize,
        wavelength_nm: f64,
    ) -> Result<(), CameraRadianceError> {
        if matches!(self, Self::Thermal(_) | Self::Coating(_)) {
            let local = scene.medium_at(medium, position, wavelength_nm)?;
            self.admit_phase_index(phase_index(local, wavelength_nm)?)?;
        }
        Ok(())
    }

    pub(super) fn admit_phase_index(&self, index: f64) -> Result<(), CameraRadianceError> {
        if matches!(self, Self::Thermal(_) | Self::Coating(_)) && index != 1.0 {
            return Err(CameraRadianceError::UnsupportedSurface {
                object_id: self.object_id(),
            });
        }
        Ok(())
    }

    pub(super) fn from_camera(
        camera: &CurvedCameraScene<'e, '_>,
        object_id: u64,
        side: Option<BoundarySide>,
    ) -> Result<Self, AreaEmitterError> {
        let invalid = || AreaEmitterError::InvalidEmitter { object_id };
        let slot = camera
            .scene
            .bindings
            .binary_search_by_key(&object_id, |binding| binding.object_id)
            .map_err(|_| invalid())?;
        let coating = camera.coating_thermal_emitters.get(slot);
        let has_coating = coating.is_some_and(|faces| faces.iter().any(Option::is_some));
        let emitter = match (
            camera.emitters.get(slot).and_then(|entry| *entry),
            camera
                .diffuse_thermal_emitters
                .get(slot)
                .and_then(Option::as_ref),
            has_coating,
        ) {
            (Some(emitter), None, false) => Self::Calibrated(emitter),
            (None, Some(emitter), false) => Self::Thermal(*emitter),
            (None, None, true) => {
                let faces = coating.ok_or_else(invalid)?;
                let source = match side {
                    Some(BoundarySide::Positive) => faces[0].as_ref(),
                    Some(BoundarySide::Negative) => faces[1].as_ref(),
                    Some(BoundarySide::OnSurface) => None,
                    None => match (&faces[0], &faces[1]) {
                        (Some(source), None) | (None, Some(source)) => Some(source),
                        _ => None,
                    },
                }.ok_or_else(invalid)?;
                Self::Coating(source.clone())
            }
            _ => return Err(invalid()),
        };
        if side.is_some_and(|side| !emitter.supports_side(side)) { return Err(invalid()); }
        Ok(emitter)
    }

    pub(super) fn matches_camera(&self, camera: &CurvedCameraScene<'_, '_>, slot: usize) -> bool {
        let coating = camera.coating_thermal_emitters.get(slot);
        let has_coating = coating.is_some_and(|faces| faces.iter().any(Option::is_some));
        match (
            self,
            camera.emitters.get(slot).and_then(|entry| *entry),
            camera
                .diffuse_thermal_emitters
                .get(slot)
                .and_then(Option::as_ref),
            has_coating,
        ) {
            (Self::Calibrated(original), Some(emitter), None, false) => core::ptr::eq(*original, emitter),
            (Self::Thermal(original), None, Some(emitter), false) => original == emitter,
            (Self::Coating(original), None, None, true) => coating
                .and_then(|faces| faces[usize::from(original.side() == BoundarySide::Negative)].as_ref())
                .is_some_and(|emitter| emitter.object_id() == original.object_id()
                    && emitter.side() == original.side()
                    && emitter.source().identity() == original.source().identity()),
            _ => false,
        }
    }
}
