//! Original source laws shared by area proposals and complementary hits.
//!
//! A thermal source retains the complete native snapshot. It is never fitted
//! into a calibrated spectrum; Planck and emissivity are evaluated only at
//! the actual path wavelength. Camera-owned thermal snapshots are values,
//! while calibrated sources retain their original borrowed identity.

use super::super::{DiffuseThermalEmitter, phase_index};
use super::{
    AreaEmitterError, CameraRadianceError, CurvedCameraScene, EmissionSides, SurfaceEmitter,
};
use crate::polarization::dielectric::scene::curved::CurvedPolarizedScene;
use fs_geom::Point3;

#[derive(Clone, Copy, Debug)]
pub(super) enum AreaEmission<'e> {
    Calibrated(&'e SurfaceEmitter),
    Thermal(DiffuseThermalEmitter),
}

impl<'e> AreaEmission<'e> {
    pub(super) const fn object_id(&self) -> u64 {
        match self {
            Self::Calibrated(emitter) => emitter.object_id(),
            Self::Thermal(emitter) => emitter.object_id(),
        }
    }

    pub(super) const fn sides(&self) -> EmissionSides {
        match self {
            Self::Calibrated(emitter) => emitter.sides(),
            Self::Thermal(_) => EmissionSides::Both,
        }
    }

    pub(super) fn evaluate(&self, wavelength_nm: f64) -> Result<f64, CameraRadianceError> {
        match self {
            Self::Calibrated(emitter) => emitter.spectrum().evaluate(wavelength_nm),
            Self::Thermal(emitter) => emitter
                .surface()
                .radiance_w_m2_sr_nm(wavelength_nm)
                .map_err(CameraRadianceError::Thermal),
        }
    }

    pub(super) fn admit_local_medium(
        &self,
        scene: &CurvedPolarizedScene<'_>,
        position: Point3,
        medium: usize,
        wavelength_nm: f64,
    ) -> Result<(), CameraRadianceError> {
        if let Self::Thermal(_) = self {
            let local = scene.medium_at(medium, position, wavelength_nm)?;
            self.admit_phase_index(phase_index(local, wavelength_nm)?)?;
        }
        Ok(())
    }

    pub(super) fn admit_phase_index(&self, index: f64) -> Result<(), CameraRadianceError> {
        if matches!(self, Self::Thermal(_)) && index != 1.0 {
            return Err(CameraRadianceError::UnsupportedSurface {
                object_id: self.object_id(),
            });
        }
        Ok(())
    }

    pub(super) fn from_camera(
        camera: &CurvedCameraScene<'e, '_>,
        object_id: u64,
    ) -> Result<Self, AreaEmitterError> {
        let invalid = || AreaEmitterError::InvalidEmitter { object_id };
        let slot = camera
            .scene
            .bindings
            .binary_search_by_key(&object_id, |binding| binding.object_id)
            .map_err(|_| invalid())?;
        match (
            camera.emitters.get(slot).and_then(|entry| *entry),
            camera
                .diffuse_thermal_emitters
                .get(slot)
                .and_then(Option::as_ref),
        ) {
            (Some(emitter), None) => Ok(Self::Calibrated(emitter)),
            (None, Some(emitter)) => Ok(Self::Thermal(*emitter)),
            _ => Err(invalid()),
        }
    }

    pub(super) fn matches_camera(&self, camera: &CurvedCameraScene<'_, '_>, slot: usize) -> bool {
        match (
            self,
            camera.emitters.get(slot).and_then(|entry| *entry),
            camera
                .diffuse_thermal_emitters
                .get(slot)
                .and_then(Option::as_ref),
        ) {
            (Self::Calibrated(original), Some(emitter), None) => core::ptr::eq(*original, emitter),
            (Self::Thermal(original), None, Some(emitter)) => original == emitter,
            _ => false,
        }
    }
}
