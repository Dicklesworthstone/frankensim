//! Directional polarized Kirchhoff emission on the original native camera.
//!
//! A source is one oriented face of the same smooth conductor used by the
//! traversal. Exact optical-law matching happens once, before publication.
//! Each encounter evaluates the original layer temperatures, wavelength and
//! outgoing cosine; no scalar emissivity or averaged temperature is fitted.

use super::super::{CameraRadianceError, CurvedCameraScene, CurvedSceneError, Interaction};
use crate::gradient_index::event::BoundarySide;
use crate::polarization::dielectric::scene::curved::CurvedPolarizedScene;
use crate::polarization::dielectric::transport::PolarizationFrame;
use crate::thermal_surface::CoatingThermalEmission;
use fs_exec::Cx;
use fs_geom::Vec3;

/// One physical side of a native smooth conductor, observed at its retained
/// thermal state. Positive side follows the original geometric normal;
/// negative side opposes it. Two bindings may give the two sides independent
/// temperatures, while both must match their original native optical law.
///
/// The existing Kirchhoff owner admits a bare substrate or one exact lossless
/// or absorbing coating here. Rough surfaces, compact periodic coatings and
/// unsupported multilayer/transmitting scene representations refuse. Each
/// emitting hit must be in the original vacuum/exterior-air approximation
/// (local phase index one). No inferred cooling or thermal advance occurs.
///
/// A native accepted-state observation can be attached by cloning
/// `NativeCoatingEmission::emitter()`. Immutable film storage remains shared.
#[derive(Clone, Debug)]
pub struct CoatingThermalEmitter {
    object_id: u64,
    side: BoundarySide,
    source: CoatingThermalEmission,
}

impl CoatingThermalEmitter {
    /// Retain a previously admitted thermal source and select one physical
    /// emitting side. `OnSurface` is not an emission direction.
    pub fn try_new(
        object_id: u64,
        side: BoundarySide,
        source: CoatingThermalEmission,
    ) -> Result<Self, CameraRadianceError> {
        if side == BoundarySide::OnSurface {
            return Err(CameraRadianceError::InvalidEmitter { object_id });
        }
        Ok(Self {
            object_id,
            side,
            source,
        })
    }

    pub const fn object_id(&self) -> u64 {
        self.object_id
    }

    pub const fn side(&self) -> BoundarySide {
        self.side
    }

    /// Original layer optics, individual temperatures and source identities.
    pub const fn source(&self) -> &CoatingThermalEmission {
        &self.source
    }

    fn side_slot(&self) -> usize {
        usize::from(self.side == BoundarySide::Negative)
    }

    fn admit_scene(&self, scene: &CurvedPolarizedScene<'_>) -> Result<usize, CameraRadianceError> {
        let invalid = || CameraRadianceError::InvalidEmitter {
            object_id: self.object_id,
        };
        let index = scene
            .bindings
            .binary_search_by_key(&self.object_id, |binding| binding.object_id)
            .map_err(|_| invalid())?;
        match scene.bindings[index].interaction {
            Interaction::SmoothConductor { optics, .. }
                if scene
                    .conductors
                    .get(optics)
                    .is_some_and(|optics| self.source.matches_smooth_conductor(optics)) =>
            {
                Ok(index)
            }
            _ => Err(invalid()),
        }
    }

    /// Native unit ray/normal and the arriving camera's valid transverse axis.
    /// At exact normal incidence, the existing transport frame resolves the
    /// degenerate s/p plane; away from it, positive Q is physical s polarization.
    pub(in super::super) fn sample(
        &self,
        wavelength_nm: f64,
        outgoing: Vec3,
        geometric_normal: Vec3,
        fallback_axis: Vec3,
    ) -> Result<([f64; 2], PolarizationFrame), CameraRadianceError> {
        let cosine = geometric_normal.dot(outgoing)
            * if self.side == BoundarySide::Positive {
                1.0
            } else {
                -1.0
            };
        // Both native inputs are unit vectors. Clamp only their positive
        // normal-incidence roundoff; zero, wrong-side and nonfinite modes
        // retain the original source owner's explicit refusal.
        let cosine = if cosine > 1.0 && cosine <= 1.0 + 16.0 * f64::EPSILON {
            1.0
        } else {
            cosine
        };
        let radiance = self
            .source
            .evaluate(wavelength_nm, cosine)
            .map_err(CameraRadianceError::CoatingThermal)?
            .polarized_radiance_w_m2_sr_nm();
        let s = Vec3::new(
            outgoing.y * geometric_normal.z - outgoing.z * geometric_normal.y,
            outgoing.z * geometric_normal.x - outgoing.x * geometric_normal.z,
            outgoing.x * geometric_normal.y - outgoing.y * geometric_normal.x,
        );
        let axis = if (s.x == 0.0 && s.y == 0.0 && s.z == 0.0) || radiance[0] == radiance[1] {
            fallback_axis
        } else {
            s
        };
        let frame = PolarizationFrame::try_new(outgoing, axis).map_err(CurvedSceneError::from)?;
        Ok((radiance, frame))
    }
}

impl CurvedCameraScene<'_, '_> {
    /// Bind the original polarized coating source to matching native smooth
    /// conductor faces. Temperature snapshots can be refreshed by rebinding;
    /// the source shares its original immutable films and never re-solves heat.
    ///
    /// Opposite sides may have distinct thermal states. Duplicate same-side
    /// bindings, a different bulk/film law, or a calibrated/diffuse thermal
    /// emitter on the same object refuse transactionally. Empty input removes
    /// these sources. Allocation occurs only during binding, bounded by twice
    /// the original surface count; camera queries allocate nothing.
    ///
    /// Every camera path using the native source accumulator observes this
    /// law, including analyzed, molecular, raster and transient traversals.
    /// These directional sources are currently sampled by camera hits; they
    /// are not coerced into the uniform-radiance finite-area proposal API.
    pub fn with_coating_thermal_emitters(
        mut self,
        cx: &Cx<'_>,
        emitters: &[CoatingThermalEmitter],
    ) -> Result<Self, CameraRadianceError> {
        cx.checkpoint()?;
        let mut slots = Vec::new();
        if !emitters.is_empty() {
            let maximum = self.scene.bindings.len().saturating_mul(2);
            if emitters.len() > maximum {
                return Err(CameraRadianceError::InvalidEmitter {
                    object_id: emitters[maximum].object_id,
                });
            }
            slots
                .try_reserve_exact(self.scene.bindings.len())
                .map_err(|_| CameraRadianceError::AllocationRefused)?;
            slots.resize_with(self.scene.bindings.len(), || [None, None]);
            for emitter in emitters {
                cx.checkpoint()?;
                let index = emitter.admit_scene(self.scene)?;
                let side = emitter.side_slot();
                if slots[index][side].is_some()
                    || self.emitters[index].is_some()
                    || self
                        .diffuse_thermal_emitters
                        .get(index)
                        .is_some_and(Option::is_some)
                {
                    return Err(CameraRadianceError::InvalidEmitter {
                        object_id: emitter.object_id,
                    });
                }
                slots[index][side] = Some(emitter.clone());
            }
        }
        cx.checkpoint()?;
        self.coating_thermal_emitters = slots;
        Ok(self)
    }
}

#[cfg(test)]
mod tests;
