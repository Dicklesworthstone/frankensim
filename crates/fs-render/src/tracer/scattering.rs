//! Participating media share the surface path integrator and light proposals.
//!
//! A sampled free flight supplies spectral likelihood weights; phase sampling
//! supplies the next direction. Direct connections integrate exact grid
//! transmittance. Both MIS directions include the probability of escaping the
//! medium without another collision, in the same solid-angle measure.

use super::*;
use crate::volumes::transport::scattering::{FreeFlightSample, PhaseFunction, ScatteringGrid};
use crate::volumes::transport::{AbsorbingGrid, VolumeTransportError};

mod execution;
pub use execution::{
    render_absorbing_range_with_execution, render_absorbing_with_execution,
    render_scattering_range_with_execution, render_scattering_with_execution,
};

/// Bit-affecting physical free-flight, phase and survival-aware MIS semantics.
// v2 admits ambient grids around dielectric bodies and independently samples
// rough Fresnel events; volume-only and vacuum streams remain unchanged.
pub const SCATTERING_TRACER_SEMANTICS_VERSION: u32 = 2;

/// Progressive film bound to one immutable participating-medium snapshot.
///
/// All settings and the physical field identity are checked before every
/// append. The caller must retain the same geometry, camera, materials and
/// lights between appends, as with `AbsorbingFilm`; this is not a serialized
/// scene checkpoint.
#[derive(Debug, Clone, PartialEq)]
pub struct ScatteringFilm {
    film: Film,
    medium_identity: ContentHash,
    settings: Settings,
}

impl ScatteringFilm {
    /// Allocate a film bound to the field and complete render policy.
    pub fn try_new(medium: &ScatteringGrid<'_>, settings: &Settings) -> Result<Self, TracerError> {
        if settings.spp == 0 || settings.max_depth == 0 {
            return Err(TracerError::InvalidInput);
        }
        Ok(Self {
            film: Film::try_new(settings.width, settings.height)?,
            medium_identity: medium.content_identity(),
            settings: *settings,
        })
    }

    /// Read accumulated XYZ sums and the completed sample prefix.
    #[must_use]
    pub const fn film(&self) -> &Film {
        &self.film
    }

    /// Finish accumulation and take the ordinary output film.
    #[must_use]
    pub fn into_film(self) -> Film {
        self.film
    }

    fn checked_append(
        &mut self,
        medium: &ScatteringGrid<'_>,
        settings: &Settings,
        from: u32,
        to: u32,
    ) -> Result<&mut Film, TracerError> {
        if to < from {
            return Err(TracerError::InvalidRange { from, to });
        }
        if *settings != self.settings
            || medium.content_identity() != self.medium_identity
            || self.film.spp_done != from
            || to > settings.spp
        {
            return Err(TracerError::InvalidInput);
        }
        Ok(&mut self.film)
    }
}

/// Render coupled surface and multiple-scattering transport through a physical grid.
///
/// Density is kg/m³, mass absorption/scattering are m²/kg and world distances
/// are metres. The medium has index one and one reciprocal phase law. Camera,
/// reflected, scattered and shadow rays all read the same borrowed field.
/// The grid is the ambient index-one medium: closed dielectric interiors
/// displace it and retain their own homogeneous spectral absorption. Refraction
/// and reflection continue through that same accepted ambient field after
/// every exit. Specialized refractive NEE connectors are disabled for nonzero
/// physical grids; the complete interface random walk supplies their paths.
/// A zero-scattering field uses deterministic absorption; a vacuum preserves
/// the ordinary tracer.
/// `max_depth` bounds the total number of surface and volume vertices together.
pub fn render_scattering(
    scene: &Scene,
    medium: &ScatteringGrid<'_>,
    cx: &Cx<'_>,
    settings: &Settings,
) -> Result<Film, TracerError> {
    cx.checkpoint()?;
    let mut film = ScatteringFilm::try_new(medium, settings)?;
    render_scattering_range(scene, medium, cx, settings, &mut film, 0, settings.spp)?;
    Ok(film.into_film())
}

/// Append a complete sample prefix transactionally, preserving exact replay.
///
/// Any changed field/settings, invalid range, exhausted traversal budget,
/// cancelled operation or nonfinite contribution leaves the film unchanged.
#[allow(clippy::too_many_arguments)]
pub fn render_scattering_range(
    scene: &Scene,
    medium: &ScatteringGrid<'_>,
    cx: &Cx<'_>,
    settings: &Settings,
    film: &mut ScatteringFilm,
    from: u32,
    to: u32,
) -> Result<(), TracerError> {
    cx.checkpoint()?;
    let film = film.checked_append(medium, settings, from, to)?;
    render_medium_range(
        scene,
        cx,
        settings,
        film,
        scattering_medium(medium),
        medium.has_thermal_emission(),
        from,
        to,
    )
}

#[derive(Clone, Copy)]
pub(super) enum PathMedium<'view, 'field> {
    Absorbing(&'view AbsorbingGrid<'field>),
    Scattering(&'view ScatteringGrid<'field>),
}

fn scattering_medium<'view, 'field>(
    medium: &'view ScatteringGrid<'field>,
) -> Option<PathMedium<'view, 'field>> {
    if medium.is_vacuum() {
        None
    } else if medium.has_scattering() {
        Some(PathMedium::Scattering(medium))
    } else {
        Some(PathMedium::Absorbing(medium.absorption_grid()))
    }
}

impl PathMedium<'_, '_> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn sample_segment(
        self,
        cx: &Cx<'_>,
        origin: Point3,
        direction: Vec3,
        distance: f64,
        wavelengths: &[f64; PACKET],
        rng: &mut PathRng,
        sobol: Option<&Sobol>,
        depth: u32,
    ) -> Result<FreeFlightSample<PACKET>, VolumeTransportError> {
        match self {
            Self::Absorbing(medium) => {
                let transport = medium.transport(cx, origin, direction, distance, wavelengths)?;
                Ok(FreeFlightSample {
                    event: None,
                    weight: transport.transmittance,
                    emission_w_m2_sr_nm: transport.emission_w_m2_sr_nm,
                    escape_probability: 1.0,
                })
            }
            Self::Scattering(medium) => {
                // The third independent bounce pair was previously reserved
                // for dielectric events. Dielectrics are not admitted here,
                // so light and phase dimensions remain disjoint from flight.
                let (lane, distance_sample) = rng.path_pair(sobol, depth, 2);
                medium.sample_segment(
                    cx,
                    origin,
                    direction,
                    distance,
                    wavelengths,
                    lane,
                    distance_sample,
                )
            }
        }
    }

    pub(super) fn transmittance(
        self,
        cx: &Cx<'_>,
        origin: Point3,
        direction: Vec3,
        distance: f64,
        wavelengths: &[f64; PACKET],
    ) -> Result<[f64; PACKET], VolumeTransportError> {
        match self {
            Self::Absorbing(medium) => {
                medium.transmittance(cx, origin, direction, distance, wavelengths)
            }
            Self::Scattering(medium) => {
                medium.transmittance(cx, origin, direction, distance, wavelengths)
            }
        }
    }

    pub(super) fn phase(self) -> Option<PhaseFunction> {
        match self {
            Self::Absorbing(_) => None,
            Self::Scattering(medium) => Some(medium.phase()),
        }
    }

    /// Density of the no-collision outcome under the packet-mixture proposal.
    /// Deterministic absorbing segments have no random escape event.
    pub(super) fn escape_probability(self, transmittance: &[f64; PACKET]) -> f64 {
        match self {
            Self::Absorbing(_) => 1.0,
            Self::Scattering(_) => transmittance.iter().sum::<f64>() / PACKET as f64,
        }
    }
}

/// Shared serial accumulator for physical media; staging publishes only once.
#[allow(clippy::too_many_arguments)]
pub(super) fn render_medium_range(
    scene: &Scene,
    cx: &Cx<'_>,
    settings: &Settings,
    film: &mut Film,
    medium: Option<PathMedium<'_, '_>>,
    has_emission: bool,
    from: u32,
    to: u32,
) -> Result<(), TracerError> {
    validate_medium_scene(scene, cx, medium)?;
    let lighting = validate_scene_with_medium_emission(scene, cx, None, has_emission)?;
    if from == to {
        return Ok(());
    }
    let mut staged = Vec::new();
    staged
        .try_reserve_exact(film.xyz.len())
        .map_err(|_| TracerError::InvalidInput)?;
    for chunk in film.xyz.chunks(4096) {
        cx.checkpoint()?;
        staged.extend_from_slice(chunk);
    }
    let key = [
        (settings.seed & 0xffff_ffff) as u32,
        (settings.seed >> 32) as u32,
    ];
    let sobol = pixel_sobol(settings.sampler, settings.seed);
    let path_sobol = path_sobol(settings.sampler);
    let kn = 1.0 / y_integral();
    let render_cx = cx.with_stream_seed(settings.seed);
    for (pixel, slot) in staged.iter_mut().enumerate() {
        let pixel = u32::try_from(pixel).map_err(|_| TracerError::InvalidInput)?;
        for sample in from..to {
            cx.checkpoint()?;
            let (jx, jy, ul) = pixel_dims(settings, sobol.as_ref(), key, pixel, sample)?;
            let xyz = trace_path_with_medium(
                scene,
                &lighting,
                &render_cx,
                settings,
                kn,
                pixel,
                sample,
                jx,
                jy,
                ul,
                None,
                CameraPath::Legacy,
                false,
                path_sobol.as_ref(),
                None,
                medium,
            )?
            .xyz;
            for channel in 0..3 {
                slot[channel] += xyz[channel];
                if !slot[channel].is_finite() {
                    return Err(TracerError::InvalidInput);
                }
            }
        }
    }
    cx.checkpoint()?;
    film.xyz = staged;
    film.spp_done = to;
    film.time_mode = FilmTimeMode::Static;
    Ok(())
}

fn validate_medium_scene(
    _scene: &Scene,
    cx: &Cx<'_>,
    _medium: Option<PathMedium<'_, '_>>,
) -> Result<(), TracerError> {
    // Closed-interface nesting and unclosed media are validated by the path
    // walk. The ambient physical grid is displaced by those interiors, so
    // admitting a scene no longer requires a blanket dielectric prohibition.
    cx.checkpoint()?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub(super) fn direct_contribution(
    scene: &Scene,
    cx: &Cx<'_>,
    ray_time: Option<&PathTime>,
    settings: &Settings,
    depth: u32,
    medium: PathMedium<'_, '_>,
    phase: PhaseFunction,
    point: Point3,
    incoming_travel: Vec3,
    light: LightSample,
    wavelengths: &[f64; PACKET],
) -> Result<[f64; PACKET], TracerError> {
    let mut result = [0.0; PACKET];
    let (direction, distance, pdf, endpoint) = match light {
        LightSample::ThermalSurface(light) => {
            let displacement = light.point.delta_from(point);
            let distance = displacement.norm();
            if !(distance > 0.0 && distance.is_finite()) {
                return Ok(result);
            }
            (
                displacement.scale(1.0 / distance),
                distance,
                light.pdf_solid_angle,
                Some((light.point, light.primitive_index)),
            )
        }
        _ => {
            let Some(direct) = prepare_direct_light(point, light) else {
                return Ok(result);
            };
            let (distance, endpoint) = match direct.target {
                DirectLightTarget::Environment => (f64::INFINITY, None),
                DirectLightTarget::Rectangle {
                    primitive_index,
                    distance_m,
                    point,
                    ..
                } => (distance_m, Some((point, primitive_index))),
            };
            (direct.direction, distance, direct.pdf_solid_angle, endpoint)
        }
    };
    if !(pdf > 0.0 && pdf.is_finite()) {
        return Ok(result);
    }
    let shadow = Ray {
        origin: point,
        dir: direction,
    };
    let visible = if let Some((target, primitive)) = endpoint {
        thermal_lighting::visible_endpoint(scene, cx, ray_time, point, point, target, primitive)?
    } else {
        intersect(scene, cx, &shadow, ray_time)?.is_none()
    };
    if !visible {
        return Ok(result);
    }
    let transmittance = medium.transmittance(cx, point, direction, distance, wavelengths)?;
    let phase_value = phase.eval_cos(incoming_travel.dot(direction))?;
    let competing_pdf = phase_value * medium.escape_probability(&transmittance);
    let weight = if depth + 1 == settings.max_depth {
        1.0
    } else {
        match settings.strategy {
            DirectStrategy::Mis => balance_heuristic(1, pdf, 1, competing_pdf),
            DirectStrategy::PowerMis => power_heuristic(1, pdf, 1, competing_pdf),
            DirectStrategy::NeeOnly => 1.0,
            DirectStrategy::BsdfOnly => return Err(TracerError::InvalidInput),
        }
    };
    for (lane, &wavelength) in wavelengths.iter().enumerate() {
        let emitted = match light {
            LightSample::ThermalSurface(light) => light
                .surface
                .radiance_w_m2_sr_nm(wavelength)
                .map_err(|_| TracerError::InvalidInput)?,
            LightSample::Rectangle(light) => light.emission.0.eval(wavelength) * light.emission.1,
            LightSample::Environment(light) => light.emission.0.eval(wavelength) * light.emission.1,
        };
        result[lane] = phase_value * emitted * transmittance[lane] / pdf * weight;
        if !result[lane].is_finite() || result[lane] < 0.0 {
            return Err(TracerError::InvalidInput);
        }
    }
    Ok(result)
}
