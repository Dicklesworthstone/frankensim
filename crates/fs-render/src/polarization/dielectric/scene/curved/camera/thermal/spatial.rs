//! Explicit source bindings for native adaptive nonisothermal LTE transport.

use super::*;
#[cfg(feature = "material-state")]
pub use crate::gradient_index::thermal::AxialGasTemperatureField;
pub use crate::gradient_index::thermal::TemperatureField;

/// An immutable temperature field for one original optical medium. Kirchhoff
/// emission uses that medium's true absorption at the native optical stages.
/// The field must be supported throughout all reached stages and numerical
/// interface collars, with spatial structure resolved by the optical step cap.
/// No index, opacity, temperature interpolation or material state is inferred.
///
/// The source tolerance is in reduced spectral radiance W/(m² sr nm), per
/// accepted step before the camera prefix and n_camera². It joins the native
/// optical relative tolerance and accept/reject decision; it is not a bound
/// on total image error. Unknown infinite source tails explicitly refuse.
#[derive(Clone, Copy)]
pub struct SpatialVolumeEmitter<'a> {
    medium_index: usize,
    temperature: &'a (dyn TemperatureField + Sync),
    source_identity: ContentHash,
    absolute_tolerance: f64,
}

impl<'a> SpatialVolumeEmitter<'a> {
    /// Bind the original frozen source with a nonzero identity and a positive
    /// finite local reduced-radiance tolerance. Medium ownership is admitted
    /// by `CurvedCameraScene::with_spatial_volume_emitters`.
    pub fn try_new(
        medium_index: usize,
        temperature: &'a (dyn TemperatureField + Sync),
        source_identity: ContentHash,
        absolute_tolerance: f64,
    ) -> Result<Self, CameraRadianceError> {
        if source_identity == ContentHash([0; 32]) {
            return Err(CameraRadianceError::Thermal(
                ThermalSurfaceError::SourceIdentity,
            ));
        }
        if !absolute_tolerance.is_finite() || absolute_tolerance <= 0.0 {
            return Err(CameraRadianceError::InvalidVolumeEmitter { medium_index });
        }
        Ok(Self {
            medium_index,
            temperature,
            source_identity,
            absolute_tolerance,
        })
    }

    /// Index in the original scene's medium table.
    pub const fn medium_index(&self) -> usize {
        self.medium_index
    }
    /// Borrowed physical temperature reconstruction.
    pub const fn temperature_field(&self) -> &'a (dyn TemperatureField + Sync) {
        self.temperature
    }
    /// Identity of the caller's frozen physical temperature state.
    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }
    /// Local absolute tolerance before camera throughput and index scaling.
    pub const fn absolute_radiance_tolerance(&self) -> f64 {
        self.absolute_tolerance
    }
}

impl<'s, 'm> CurvedCameraScene<'s, 'm> {
    /// Add nonisothermal LTE sources to actual camera, molecular and raster
    /// consumers. Existing isothermal bindings remain exact and may coexist
    /// in different media. Duplicate or overlapping source laws refuse.
    /// Empty input disables only spatial sources. Preparation allocates once;
    /// the native path carries stack-only adaptive source state thereafter.
    pub fn with_spatial_volume_emitters(
        mut self,
        cx: &Cx<'_>,
        emitters: &[SpatialVolumeEmitter<'s>],
    ) -> Result<Self, CameraRadianceError> {
        cx.checkpoint()?;
        let mut slots = Vec::new();
        if !emitters.is_empty() {
            if emitters.len() > self.scene.media.len() {
                return Err(CameraRadianceError::InvalidVolumeEmitter {
                    medium_index: emitters[self.scene.media.len()].medium_index,
                });
            }
            slots
                .try_reserve_exact(self.scene.media.len())
                .map_err(|_| CameraRadianceError::AllocationRefused)?;
            slots.resize(self.scene.media.len(), None);
            for emitter in emitters {
                cx.checkpoint()?;
                let invalid = CameraRadianceError::InvalidVolumeEmitter {
                    medium_index: emitter.medium_index,
                };
                let slot = slots.get_mut(emitter.medium_index).ok_or(invalid.clone())?;
                if slot.is_some()
                    || self
                        .volume_emitters
                        .get(emitter.medium_index)
                        .is_some_and(Option::is_some)
                {
                    return Err(invalid);
                }
                *slot = Some(*emitter);
            }
        }
        cx.checkpoint()?;
        self.spatial_volume_emitters = slots;
        Ok(self)
    }
}
