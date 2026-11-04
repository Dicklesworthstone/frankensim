//! The original finite-band diffuse Kirchhoff source on native camera paths.
//!
//! Binding verifies the complete native reflectance curve against the same
//! emissivity that emits. Planck radiance is evaluated at the actual path
//! wavelength, without resampling it into a tabulated light spectrum. The
//! accepted thermal state is observed only; this binding never advances heat.

use super::super::{CameraRadianceError, CurvedCameraScene, Interaction};
use crate::polarization::dielectric::scene::curved::CurvedPolarizedScene;
use crate::spectral::SpectralReflectance;
use crate::thermal_surface::SpectralDiffuseThermalSurface;
use fs_exec::Cx;

/// An isothermal opaque diffuse surface with one original spectral law.
/// Both faces emit into their own hemisphere, as in the original diffuse
/// thermal owner. Every emitting contact must be in an index-one medium;
/// a vacuum Planck law is not silently promoted into a refractive host.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiffuseThermalEmitter {
    object_id: u64,
    surface: SpectralDiffuseThermalSurface,
}

impl DiffuseThermalEmitter {
    pub const fn new(object_id: u64, surface: SpectralDiffuseThermalSurface) -> Self {
        Self { object_id, surface }
    }

    pub const fn object_id(&self) -> u64 {
        self.object_id
    }

    pub const fn surface(&self) -> &SpectralDiffuseThermalSurface {
        &self.surface
    }

    /// Prepare the exact opaque reciprocal curve for the native scene's
    /// `DepolarizingDiffuse` material. This copies at most 16 inline knots;
    /// no light/ray stores an independently fitted emissivity or RGB value.
    pub fn reflectance_spectrum(&self) -> Result<SpectralReflectance, CameraRadianceError> {
        let emissivity = self.surface.emissivity();
        let mut knots = [(0.0, 0.0); 16];
        for (slot, &(wavelength, epsilon)) in knots.iter_mut().zip(emissivity.knots()) {
            *slot = (wavelength, 1.0 - epsilon);
        }
        SpectralReflectance::table_in_band(&knots[..emissivity.knots().len()]).map_err(|_| {
            CameraRadianceError::InvalidEmitter {
                object_id: self.object_id,
            }
        })
    }

    // Camera hits and native area proposals must admit the same complete
    // Kirchhoff material; neither consumer may substitute a calibrated lamp.
    pub(in super::super) fn admit_scene(
        &self,
        scene: &CurvedPolarizedScene<'_>,
    ) -> Result<usize, CameraRadianceError> {
        let invalid = || CameraRadianceError::InvalidEmitter {
            object_id: self.object_id,
        };
        let index = scene
            .bindings
            .binary_search_by_key(&self.object_id, |binding| binding.object_id)
            .map_err(|_| invalid())?;
        match scene.bindings[index].interaction {
            Interaction::DepolarizingDiffuse { reflectance, .. } => {
                if scene.reflectances[reflectance] != self.reflectance_spectrum()? {
                    return Err(invalid());
                }
            }
            Interaction::Opaque
                if self
                    .surface
                    .emissivity()
                    .knots()
                    .iter()
                    .all(|&(_, epsilon)| epsilon == 1.0) => {}
            _ => return Err(invalid()),
        }
        Ok(index)
    }
}

impl CurvedCameraScene<'_, '_> {
    /// Attach original thermal surfaces to existing opaque diffuse geometry.
    /// The native reflectance must exactly match `reflectance_spectrum()` in
    /// both knots and values; an `Opaque` termination is admitted only when
    /// every source knot is a blackbody (emissivity one). A calibrated emitter
    /// on the same object or any duplicate thermal binding refuses.
    ///
    /// Preparation allocates one indexed table bounded by native surface
    /// count. Paths and detector lanes query it without allocation. Empty
    /// input removes these sources; the native material remains its owner.
    /// Surface-only and molecular/particle-scattering paths share this source
    /// accumulator, including after earlier emission or optical boundaries.
    pub fn with_diffuse_thermal_emitters(
        mut self,
        cx: &Cx<'_>,
        emitters: &[DiffuseThermalEmitter],
    ) -> Result<Self, CameraRadianceError> {
        cx.checkpoint()?;
        let mut slots = Vec::new();
        if !emitters.is_empty() {
            if emitters.len() > self.scene.bindings.len() {
                return Err(CameraRadianceError::InvalidEmitter {
                    object_id: emitters[self.scene.bindings.len()].object_id,
                });
            }
            slots
                .try_reserve_exact(self.scene.bindings.len())
                .map_err(|_| CameraRadianceError::AllocationRefused)?;
            slots.resize(self.scene.bindings.len(), None);
            for emitter in emitters {
                cx.checkpoint()?;
                let invalid = CameraRadianceError::InvalidEmitter {
                    object_id: emitter.object_id,
                };
                let index = emitter.admit_scene(self.scene)?;
                if slots[index].is_some() || self.emitters[index].is_some() {
                    return Err(invalid);
                }
                slots[index] = Some(*emitter);
            }
        }
        cx.checkpoint()?;
        self.diffuse_thermal_emitters = slots;
        Ok(self)
    }
}
