//! Isothermal LTE emission on the original curved absorption integral.
//!
//! With vacuum wavelength fixed, the refractive transfer equation evolves
//! reduced radiance L/n². Kirchhoff's law supplies local source n² alpha B(T),
//! so a leg contributes n_camera² I_start B(T) (1-exp(-Delta tau_abs)). The
//! source's local phase index cancels; a group index never enters radiance.
//! Temperature is constant within each explicitly bound medium, while its
//! original phase index and true absorption may vary spatially. Molecular
//! free-flight probability already supplies scattering survival; it is not
//! added to the optical depth or multiplied into this contribution again.
//!
//! References: Zhao, Tan and Liu (2011), Eq. 21,
//! <https://arxiv.org/abs/1110.5134>; Liu, Shi and Liu (2018), Eq. 1,
//! <https://doi.org/10.1177/1687814018809613>.

use super::{
    BeamMedium, CameraRadianceError, CurvedCameraScene, CurvedPathResult, PathTermination,
    finite_nonnegative,
};
use crate::thermal_surface::{ThermalSurfaceError, planck_radiance_w_m2_sr_nm};
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_math::det;

mod spatial;
pub use spatial::{SpatialVolumeEmitter, TemperatureField};
mod surface;
pub use surface::DiffuseThermalEmitter;
mod coating;
pub use coating::CoatingThermalEmitter;
/// Kirchhoff emission consistent with complete rough conductor reflection.
pub mod smith;
#[cfg(feature = "material-state")]
pub use spatial::AxialGasTemperatureField;

/// An explicit isothermal, isotropic local-thermodynamic-equilibrium source.
/// The selected native medium owns refraction and true absorption; this
/// binding uses that same absorption in Kirchhoff's law. It adds no opacity,
/// scattering extinction, surface emissivity or independent thermal marcher.
///
/// The source identity names the frozen physical temperature state supplied
/// by the caller. This declaration does not infer temperature from optical
/// data or update the thermal state's energy. Nonuniform temperature, gain
/// and polarized thermal emission are outside this source law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsothermalVolumeEmitter {
    medium_index: usize,
    temperature_k: f64,
    source_identity: ContentHash,
}

impl IsothermalVolumeEmitter {
    /// Declare LTE for one original scene medium at a positive temperature.
    /// The camera validates medium ownership and duplicate bindings.
    pub fn try_new(
        medium_index: usize,
        temperature_k: f64,
        source_identity: ContentHash,
    ) -> Result<Self, ThermalSurfaceError> {
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(ThermalSurfaceError::Temperature);
        }
        if source_identity == ContentHash([0; 32]) {
            return Err(ThermalSurfaceError::SourceIdentity);
        }
        Ok(Self {
            medium_index,
            temperature_k,
            source_identity,
        })
    }

    /// Index in the original curved scene's medium table.
    pub const fn medium_index(&self) -> usize {
        self.medium_index
    }

    /// Fixed thermodynamic temperature of this medium, in kelvin.
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }

    /// Identity of the caller's frozen physical temperature state.
    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }
}

impl CurvedCameraScene<'_, '_> {
    /// Bind optional LTE emission to original media. Unlisted media retain
    /// their existing absorption without an inferred thermal source. Empty
    /// input disables volume emission and retains the ordinary camera bits.
    /// Binding allocates once; rays, molecular paths and pixels use stack state.
    pub fn with_isothermal_volume_emitters(
        mut self,
        cx: &Cx<'_>,
        emitters: &[IsothermalVolumeEmitter],
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
                    || self.spatial_volume_emitters.get(emitter.medium_index).is_some_and(Option::is_some)
                {
                    return Err(invalid);
                }
                *slot = Some(*emitter);
            }
        }
        cx.checkpoint()?;
        self.volume_emitters = slots;
        Ok(self)
    }
}

/// One current leg. Cache Planck evaluation across molecular collisions in
/// the same medium; all spectral and source ownership stays with the camera.
pub(super) struct ThermalEmission {
    optical_depth: f64,
    departure_power: f64,
    source: Option<(usize, f64)>,
}

impl ThermalEmission {
    pub(super) const fn new() -> Self {
        Self {
            optical_depth: 0.0,
            departure_power: 1.0,
            source: None,
        }
    }

    pub(super) fn depart(&mut self, path: &CurvedPathResult) {
        self.optical_depth = path.optical_depth();
        self.departure_power = path.beam().power().i;
    }

    pub(super) fn arrive(
        &mut self,
        camera: &CurvedCameraScene<'_, '_>,
        initial_index: f64,
        path: &CurvedPathResult,
        terminal: bool,
    ) -> Result<f64, CameraRadianceError> {
        let depth = finite_nonnegative(path.optical_depth() - self.optical_depth)?;
        let wavelength = path.beam().wavelength_nm();
        let tail = terminal
            && path.termination() == PathTermination::Escaped
            && match path.beam().medium() {
                BeamMedium::Homogeneous(glass) => {
                    glass.absorption().extinction_per_m(wavelength)? > 0.0
                }
                _ => false,
            };
        if depth == 0.0 && !tail {
            self.depart(path);
            return Ok(0.0);
        }
        let medium = path.medium_index();
        let source = if let Some((cached, value)) = self.source
            && cached == medium
        {
            value
        } else {
            let value = match camera.volume_emitters.get(medium).and_then(Option::as_ref) {
                Some(emitter) => planck_radiance_w_m2_sr_nm(emitter.temperature_k, wavelength)
                    .map_err(CameraRadianceError::Thermal)?,
                None => 0.0,
            };
            self.source = Some((medium, value));
            value
        };
        let mut emitted = weighted_source(
            source,
            initial_index,
            self.departure_power,
            -det::expm1(-depth),
        )?;
        self.depart(path);
        if tail {
            emitted = finite_nonnegative(
                emitted + weighted_source(source, initial_index, self.departure_power, 1.0)?,
            )?;
        }
        Ok(emitted)
    }
}

pub(super) fn weighted_source(
    source: f64,
    initial_index: f64,
    power: f64,
    absorbed: f64,
) -> Result<f64, CameraRadianceError> {
    if source == 0.0 || power == 0.0 || absorbed == 0.0 {
        return Ok(0.0);
    }
    let direct = ((initial_index * initial_index) * source) * (power * absorbed);
    if direct.is_finite() && direct > 0.0 {
        return Ok(direct);
    }
    // A tiny absorption fraction can rescue a large source, or a bright
    // source can rescue a tiny prefix; reject only the final numeric range.
    let value = det::exp(
        2.0 * det::ln(initial_index) + det::ln(source) + det::ln(power) + det::ln(absorbed),
    );
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(CameraRadianceError::NumericalRange)
    }
}
