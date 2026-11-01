//! Molecular multiple scattering along actual curved physical paths.
//!
//! Each declared medium has a spatially uniform mass-scattering spectrum and
//! density. Its original phase-index field may vary, and that original field
//! continues to own absorption and group delay. An exponential scattering-only
//! flight is measured in arc length, not endpoint separation or a tangent ray.
//! Native geometry decides whether a surface interrupts that flight. At a
//! collision the original Rayleigh Mueller/PDF sampler acts on the locally
//! transported polarization frame and current physical refractive index.
//!
//! Sampling scattering separately from continuous Beer absorption is the
//! survival-biased power estimator: scattering survival is already represented
//! by flight probability, so no extra scattering extinction or albedo weight
//! multiplies the packet. Source power receives no radiance/index-squared
//! factor. A volume collision may redirect a ray toward a previously missed
//! detector; an active scattering medium cannot use a straight escape proof.
//!
//! This owner does not turn a heterogeneous density grid into constant cells or
//! trace its chords. Spatially varying scattering requires a curved optical-
//! depth event owner. Only molecular Rayleigh/depolarized-Rayleigh phase laws
//! have the requisite native Stokes response; scalar particle laws refuse.
//! The medium snapshot must remain immutable for the complete path.

use super::{
    CurvedPathResult, CurvedPolarizedScene, CurvedSceneError, CurvedTraceSettings,
    CurvedVolumeFlight, MAX_POLARIZED_PATH_EVENTS, PolarizedBeam, SceneError,
};
use crate::volumes::transport::VolumeTransportError;
use crate::volumes::transport::scattering::polarized::PolarizedScatteringError;
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_geom::Point3;
use fs_math::det;
use fs_rand::philox::philox4x32_10;

/// Immutable density snapshot queried only at delta-tracking candidates.
/// Implementations must return kg/m³ for one frozen physical state.
pub trait CurvedDensityField: Sync {
    fn density_kg_m3(&self, position: Point3) -> Result<f64, VolumeTransportError>;
}

#[derive(Clone, Copy)]
enum CurvedDensityProfile<'a> {
    Uniform(f64),
    Field {
        field: &'a (dyn CurvedDensityField + Sync),
        majorant_kg_m3: f64,
        identity: ContentHash,
    },
}

/// Immutable scattering law for one explicit native medium. Absorption remains
/// owned by the original refractive field. Heterogeneous density uses unbiased
/// delta tracking against an explicit global density majorant.
#[derive(Clone, Copy)]
pub struct CurvedRayleighMedium<'a> {
    scattering: &'a MassScatteringSpectrum,
    density: CurvedDensityProfile<'a>,
    phase: &'a PhaseFunction,
}

impl<'a> CurvedRayleighMedium<'a> {
    /// Admit kg/m³ and m²/kg without inferring either from refractive index.
    /// Every declared phase is checked, including a zero-density medium.
    pub fn try_new(
        scattering: &'a MassScatteringSpectrum,
        density_kg_m3: f64,
        phase: &'a PhaseFunction,
    ) -> Result<Self, CurvedSceneError> {
        if !density_kg_m3.is_finite() || density_kg_m3 < 0.0 {
            return Err(volume_error(VolumeTransportError::Density { index: 0 }));
        }
        phase.rayleigh_mueller(0.0).map_err(scattering_error)?;
        Ok(Self {
            scattering,
            density: CurvedDensityProfile::Uniform(density_kg_m3),
            phase,
        })
    }

    /// Bind a frozen heterogeneous density field. The positive finite majorant
    /// must bound the field everywhere this medium can be reached.
    pub fn try_new_heterogeneous(
        scattering: &'a MassScatteringSpectrum,
        density: &'a (dyn CurvedDensityField + Sync),
        density_majorant_kg_m3: f64,
        density_identity: ContentHash,
        phase: &'a PhaseFunction,
    ) -> Result<Self, CurvedSceneError> {
        if !density_majorant_kg_m3.is_finite()
            || density_majorant_kg_m3 <= 0.0
            || density_identity == ContentHash([0; 32])
        {
            return Err(volume_error(VolumeTransportError::Density { index: 0 }));
        }
        phase.rayleigh_mueller(0.0).map_err(scattering_error)?;
        Ok(Self {
            scattering,
            density: CurvedDensityProfile::Field {
                field: density,
                majorant_kg_m3: density_majorant_kg_m3,
                identity: density_identity,
            },
            phase,
        })
    }

    pub const fn scattering_spectrum(&self) -> &'a MassScatteringSpectrum {
        self.scattering
    }
    /// Uniform density, or the global majorant for a heterogeneous field.
    pub const fn density_kg_m3(&self) -> f64 {
        match self.density {
            CurvedDensityProfile::Uniform(value) => value,
            CurvedDensityProfile::Field { majorant_kg_m3, .. } => majorant_kg_m3,
        }
    }
    pub const fn is_heterogeneous(&self) -> bool {
        matches!(self.density, CurvedDensityProfile::Field { .. })
    }
    pub const fn density_identity(&self) -> Option<ContentHash> {
        match self.density {
            CurvedDensityProfile::Uniform(_) => None,
            CurvedDensityProfile::Field { identity, .. } => Some(identity),
        }
    }
    pub const fn phase(&self) -> &'a PhaseFunction {
        self.phase
    }

    /// Majorant scattering rate per metre of ray arc length.
    pub fn scattering_per_m(&self, wavelength_nm: f64) -> Result<f64, CurvedSceneError> {
        let mass = self.scattering.at(wavelength_nm).map_err(volume_error)?;
        let density = self.density_kg_m3();
        let rate = mass * density;
        if !rate.is_finite() || (mass > 0.0 && density > 0.0 && rate == 0.0) {
            return Err(volume_error(VolumeTransportError::Unrepresentable));
        }
        Ok(rate)
    }

    fn accepts_candidate(
        &self,
        position: Point3,
        acceptance_sample: f64,
    ) -> Result<bool, CurvedSceneError> {
        match self.density {
            CurvedDensityProfile::Uniform(_) => Ok(true),
            CurvedDensityProfile::Field {
                field,
                majorant_kg_m3,
                ..
            } => {
                if !acceptance_sample.is_finite() || !(0.0..1.0).contains(&acceptance_sample) {
                    return Err(volume_error(VolumeTransportError::Sample));
                }
                let density = field.density_kg_m3(position).map_err(volume_error)?;
                if !density.is_finite() || density < 0.0 || density > majorant_kg_m3 {
                    return Err(volume_error(VolumeTransportError::Density { index: 0 }));
                }
                if density == 0.0 {
                    return Ok(false);
                }
                let probability = density / majorant_kg_m3;
                if !probability.is_finite() || probability <= 0.0 {
                    return Err(volume_error(VolumeTransportError::Unrepresentable));
                }
                Ok(acceptance_sample < probability)
            }
        }
    }
}

/// Independent physical-distance, surface-event and volume-collision limits.
#[derive(Clone, Copy, Debug)]
pub struct CurvedRayleighTraceSettings {
    pub surfaces: CurvedTraceSettings,
    pub max_volume_events: u32,
}

/// The original detector-consumable path plus actual volume-work counts.
/// Surviving power at a volume cap has `PathTermination::VolumeLimit`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurvedRayleighPath {
    path: CurvedPathResult,
    volume_events: u32,
    free_flight_draws: u32,
}
impl CurvedRayleighPath {
    pub const fn path(&self) -> &CurvedPathResult {
        &self.path
    }
    pub const fn volume_events(&self) -> u32 {
        self.volume_events
    }
    pub const fn free_flight_draws(&self) -> u32 {
        self.free_flight_draws
    }
}

/// Borrowed original geometry and one optional scattering law per medium.
/// Binding and tracing allocate no heap workspace or copied optical tables.
pub struct CurvedRayleighScene<'a, 'scene, 'optics> {
    scene: &'a CurvedPolarizedScene<'scene>,
    media: &'a [Option<CurvedRayleighMedium<'optics>>],
}
impl<'a, 'scene, 'optics> CurvedRayleighScene<'a, 'scene, 'optics> {
    pub fn try_new(
        cx: &Cx<'_>,
        scene: &'a CurvedPolarizedScene<'scene>,
        media: &'a [Option<CurvedRayleighMedium<'optics>>],
    ) -> Result<Self, CurvedSceneError> {
        cx.checkpoint()?;
        if media.len() != scene.media.len() {
            return Err(SceneError::InvalidBindings.into());
        }
        Ok(Self { scene, media })
    }

    pub const fn surface_scene(&self) -> &'a CurvedPolarizedScene<'scene> {
        self.scene
    }

    /// Surface draws keep the native domain and surface ordinal. Flights use a
    /// disjoint Philox domain and their own ordinal, including interrupted
    /// flights. Sharing path IDs across wavelengths intentionally correlates
    /// spectral estimates. Zero scattering consumes no random flight draw.
    #[allow(clippy::too_many_arguments)]
    pub fn trace(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<CurvedRayleighPath, CurvedSceneError> {
        let draws = |ordinal, domain| {
            let block = philox4x32_10(
                [path_id as u32, (path_id >> 32) as u32, ordinal, domain],
                [seed as u32, (seed >> 32) as u32],
            );
            core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
        };
        self.trace_with_tracking_samples_impl(
            cx,
            beam,
            initial_medium,
            settings,
            |event| draws(event, 0x706f_6c32),
            |flight| draws(flight, 0x706f_6c34),
            true,
        )
    }

    /// Flight draws are [CDF variate, phase cosine variate, phase azimuth
    /// variate], all in [0,1). The first gives `-ln(1-u)/sigma_s` metres of
    /// actual arc length. A surface winning the race starts a fresh flight in
    /// its selected medium, using exponential memorylessness. Interfaces and
    /// colliding packets retain their original physical power exactly once.
    /// Deterministic custom sampler for the original uniform-density path.\n    /// Heterogeneous density needs an independent fourth acceptance variate\n    /// and therefore refuses here rather than correlating acceptance and phase.\n    #[allow(clippy::too_many_arguments)]\n    pub fn trace_with_samples(\n        &self,\n        cx: &Cx<'_>,\n        beam: PolarizedBeam,\n        initial_medium: usize,\n        settings: CurvedRayleighTraceSettings,\n        surface_samples: impl FnMut(u32) -> [f64; 3],\n        mut flight_samples: impl FnMut(u32) -> [f64; 3],\n    ) -> Result<CurvedRayleighPath, CurvedSceneError> {\n        self.trace_with_tracking_samples_impl(\n            cx, beam, initial_medium, settings, surface_samples,\n            |flight| { let d = flight_samples(flight); [d[0], d[1], d[2], 0.0] },\n            false,\n        )\n    }\n\n    #[allow(clippy::too_many_arguments)]
    fn trace_with_tracking_samples_impl(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        mut flight_samples: impl FnMut(u32) -> [f64; 4],
        heterogeneous_samples: bool,
    ) -> Result<CurvedRayleighPath, CurvedSceneError> {
        cx.checkpoint()?;
        if settings.max_volume_events == 0 || settings.max_volume_events > MAX_POLARIZED_PATH_EVENTS
        {
            return Err(SceneError::InvalidSettings.into());
        }
        let mut free_flight_draws = 0_u32;
        let (path, volume_events) = self.scene.trace_with_flights(
            cx,
            beam,
            initial_medium,
            settings.surfaces,
            surface_samples,
            settings.max_volume_events,
            |medium, wavelength| {
                let Some(law) = self.media[medium] else {
                    return Ok(None);
                };
                if law.is_heterogeneous() && !heterogeneous_samples {
                    return Err(volume_error(VolumeTransportError::Sample));
                }
                let rate = law.scattering_per_m(wavelength)?;
                if rate == 0.0 {
                    return Ok(None);
                }
                let draws = flight_samples(free_flight_draws);
                free_flight_draws = free_flight_draws
                    .checked_add(1)
                    .ok_or(CurvedSceneError::NumericalRange)?;
                if !draws
                    .into_iter()
                    .all(|u| u.is_finite() && (0.0..1.0).contains(&u))
                {
                    return Err(volume_error(VolumeTransportError::Sample));
                }
                // Compensate the rounded subtraction in 1-u. If it rounds to
                // one, u itself is the correctly resolved first-order depth.
                let complement = 1.0 - draws[0];
                let depth = if complement == 1.0 {
                    draws[0]
                } else {
                    -det::ln(complement) * (draws[0] / (1.0 - complement))
                };
                let distance_m = depth / rate;
                if distance_m.is_nan() || distance_m < 0.0 || (depth > 0.0 && distance_m == 0.0) {
                    return Err(volume_error(VolumeTransportError::Unrepresentable));
                }
                Ok(Some(CurvedVolumeFlight {
                    distance_m,
                    phase_samples: [draws[1], draws[2]],
                    acceptance_sample: draws[3],
                }))
            },
            |medium, beam, samples, acceptance_sample| {
                let law = self.media[medium].ok_or(CurvedSceneError::NumericalRange)?;
                if !law.accepts_candidate(beam.ray().origin, acceptance_sample)? {
                    return Ok(None);
                }
                Ok(Some(law.phase.sample_polarized(beam, samples).map_err(scattering_error)?.beam))
            },
        )?;
        cx.checkpoint()?;
        Ok(CurvedRayleighPath {
            path,
            volume_events,
            free_flight_draws,
        })
    }
}

fn scattering_error(error: PolarizedScatteringError) -> CurvedSceneError {
    SceneError::from(error).into()
}
fn volume_error(error: VolumeTransportError) -> CurvedSceneError {
    scattering_error(PolarizedScatteringError::Volume(error))
}

#[cfg(test)]
mod tests;
