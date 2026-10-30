//! Calibrated point fixtures feed the existing spectral surface/volume paths.
//!
//! Every admitted fixture is evaluated at each scattering vertex. These delta
//! position sources have no competing continuous BSDF/phase hit technique, so
//! their direct contribution has unit MIS weight. The established area-light
//! mixture and its random-number consumption are unchanged.

use super::*;
use crate::camera::{Aperture, PhysicalCamera};
use crate::photometry::area::PhotometricEmitters;
use crate::photometry::medium::HomogeneousMedium;
use crate::photometry::source::{LUMINOUS_EFFICACY_LM_W, PhotometricSource};
use crate::photometry::transport::PhotometricMedium;
use crate::photometry::{
    IlluminanceEstimate, IlluminanceProbe, PhotometricScene, PhotometryError, TransportOptions,
};
use crate::volumes::transport::AbsorbingGrid;
use crate::volumes::transport::scattering::{PhaseFunction, ScatteringGrid};

pub(super) mod area;

#[cfg(test)]
mod thermal_range_tests;

/// Borrowed physical medium used by the calibrated path estimator.
#[derive(Clone, Copy)]
pub enum FixtureMedium<'view, 'field> {
    /// No participating medium.
    Vacuum,
    /// Borrowed SI absorption and optional local-equilibrium emission.
    Absorbing(&'view AbsorbingGrid<'field>),
    /// Borrowed SI absorption, scattering and optional thermal emission.
    Scattering(&'view ScatteringGrid<'field>),
}

/// Immutable fixture state and finite transport support for one render job.
pub(super) struct FixtureLighting<'a> {
    sources: &'a [PhotometricSource],
    pub(super) areas: area::AreaLighting<'a>,
    options: TransportOptions,
    medium_bounds_m: Option<([f64; 3], [f64; 3])>,
    scattering: bool,
}

struct FixtureSegment {
    direction: Vec3,
    distance: f64,
    medium_interval: Option<(f64, f64)>,
}

/// One admitted scene, fixture set, medium and sample policy. All retained
/// assets are borrowed, including the original mesh BVHs and physical cells.
pub(crate) struct FixtureTransport<'scene, 'field> {
    scene: &'scene Scene,
    direct_scene: PhotometricScene<'scene>,
    fixtures: FixtureLighting<'scene>,
    lighting: SceneLighting<'scene>,
    medium: Option<PathMedium<'scene, 'field>>,
    settings: Settings,
    key: [u32; 2],
    sobol: Option<Sobol>,
    path_sobol: Option<Sobol>,
}

impl<'scene, 'field> FixtureTransport<'scene, 'field> {
    pub(crate) fn try_new_with_emitters(
        scene: &'scene Scene,
        emitters: PhotometricEmitters<'scene>,
        medium: FixtureMedium<'scene, 'field>,
        options: TransportOptions,
        settings: Settings,
        cx: &Cx<'_>,
    ) -> Result<Self, PhotometryError> {
        Self::try_admit(scene, emitters, medium, options, settings, true, cx)
    }

    /// Admit probe addresses without reading or validating the unused camera.
    /// The receiver index occupies the existing deterministic pixel address.
    pub(crate) fn try_for_probes_with_emitters(
        scene: &'scene Scene,
        emitters: PhotometricEmitters<'scene>,
        medium: FixtureMedium<'scene, 'field>,
        options: TransportOptions,
        settings: Settings,
        probe_count: u32,
        cx: &Cx<'_>,
    ) -> Result<Self, PhotometryError> {
        cx.checkpoint()?;
        if probe_count == 0 {
            return Err(PhotometryError::InvalidOptions);
        }
        let settings = Settings {
            width: probe_count,
            height: 1,
            ..settings
        };
        Self::try_admit(scene, emitters, medium, options, settings, false, cx)
    }

    fn try_admit(
        scene: &'scene Scene,
        emitters: PhotometricEmitters<'scene>,
        medium: FixtureMedium<'scene, 'field>,
        options: TransportOptions,
        settings: Settings,
        admit_camera: bool,
        cx: &Cx<'_>,
    ) -> Result<Self, PhotometryError> {
        cx.checkpoint()?;
        let sources = emitters.points;
        if settings.spp == 0
            || settings.max_depth == 0
            || (settings.strategy == DirectStrategy::BsdfOnly
                && sources.iter().any(|source| source.flux_lm() > 0.0))
            || !matches!(settings.sampler, Sampler::Iid | Sampler::OwenSobol)
        {
            return Err(PhotometryError::InvalidOptions);
        }
        let count = checked_pixel_len(settings.width, settings.height)?;
        if count > options.max_outputs {
            return Err(PhotometryError::OutputBudgetExceeded);
        }
        // This admission boundary must never assign SI units to legacy RGB/
        // relative emitters. Physical cell emission already has SI radiance.
        if !scene.lights.is_empty() || scene.environment.is_some() {
            return Err(PhotometryError::InvalidOptions);
        }
        let direct_medium = match medium {
            FixtureMedium::Vacuum => PhotometricMedium::Homogeneous(HomogeneousMedium::vacuum()),
            FixtureMedium::Absorbing(grid) => PhotometricMedium::Absorbing(grid),
            FixtureMedium::Scattering(grid) => PhotometricMedium::Scattering(grid),
        };
        let direct_scene = PhotometricScene::try_with_path_medium(
            &scene.primitives,
            sources,
            direct_medium,
            options,
            cx,
        )?;
        if admit_camera {
            PhysicalCamera::try_legacy_compatible(
                scene.camera.eye,
                scene.camera.forward,
                scene.camera.up,
                scene.camera.half_tan,
                1.0,
                Aperture::try_circular(0.0)?,
            )?;
        }
        let areas = area::AreaLighting::try_new(scene, emitters.surfaces, cx)?;
        let (medium, medium_bounds_m, scattering) = match medium {
            FixtureMedium::Vacuum => (None, None, false),
            FixtureMedium::Absorbing(grid) if grid.is_vacuum() => (None, None, false),
            FixtureMedium::Absorbing(grid) => (
                Some(PathMedium::Absorbing(grid)),
                Some(grid.bounds_m()),
                false,
            ),
            FixtureMedium::Scattering(grid) if grid.is_vacuum() => (None, None, false),
            FixtureMedium::Scattering(grid) => (
                if grid.has_scattering() {
                    Some(PathMedium::Scattering(grid))
                } else {
                    Some(PathMedium::Absorbing(grid.absorption_grid()))
                },
                Some(grid.absorption_grid().bounds_m()),
                grid.has_scattering(),
            ),
        };
        // Point fixtures are sampled separately, so an empty continuous light
        // mixture is intentional and does not mean the calibrated scene is dark.
        let lighting = validate_scene_with_medium_emission(scene, cx, None, true)?;
        let key = [
            (settings.seed & 0xffff_ffff) as u32,
            (settings.seed >> 32) as u32,
        ];
        let sobol = pixel_sobol(settings.sampler, settings.seed);
        let path_sobol = path_sobol(settings.sampler);
        cx.checkpoint()?;
        Ok(Self {
            scene,
            direct_scene,
            fixtures: FixtureLighting {
                sources,
                areas,
                options,
                medium_bounds_m,
                scattering,
            },
            lighting,
            medium,
            settings,
            key,
            sobol,
            path_sobol,
        })
    }

    /// Return one unaccumulated absolute photometric XYZ sample. Y is cd/m²;
    /// there is no display white normalization, exposure or tone mapping.
    pub(crate) fn sample(
        &self,
        cx: &Cx<'_>,
        pixel: u32,
        sample: u32,
    ) -> Result<[f64; 3], TracerError> {
        cx.checkpoint()?;
        if pixel >= self.settings.width * self.settings.height || sample >= self.settings.spp {
            return Err(TracerError::InvalidInput);
        }
        let (jx, jy, ul) =
            pixel_dims(&self.settings, self.sobol.as_ref(), self.key, pixel, sample)?;
        let render_cx = cx.with_stream_seed(self.settings.seed);
        let result = trace_path_with_emitters(
            self.scene,
            &self.lighting,
            &render_cx,
            &self.settings,
            LUMINOUS_EFFICACY_LM_W,
            pixel,
            sample,
            jx,
            jy,
            ul,
            None,
            CameraPath::Legacy,
            false,
            self.path_sobol.as_ref(),
            None,
            self.medium,
            Some(&self.fixtures),
        )?;
        if result.xyz.iter().any(|value| !value.is_finite()) {
            return Err(TracerError::InvalidInput);
        }
        cx.checkpoint()?;
        Ok(result.xyz)
    }

    /// The delta-position fixture contribution is integrated separately from
    /// hemisphere samples, which have zero probability of finding a point.
    pub(crate) fn direct_illuminance(
        &self,
        cx: &Cx<'_>,
        probe: &IlluminanceProbe,
    ) -> Result<IlluminanceEstimate, PhotometryError> {
        self.direct_scene.illuminance(cx, probe)
    }

    /// One cosine-weighted sample of continuous incoming irradiance. Y is lux;
    /// the direct delta-fixture term is supplied by `direct_illuminance`.
    /// Reflected/scattered fixture light and physical medium emission use the
    /// original spectral path estimator, including its finite depth ceiling.
    pub(crate) fn sample_probe(
        &self,
        cx: &Cx<'_>,
        probe: &IlluminanceProbe,
        probe_index: u32,
        sample: u32,
    ) -> Result<[f64; 3], TracerError> {
        cx.checkpoint()?;
        if probe_index >= self.settings.width * self.settings.height || sample >= self.settings.spp
        {
            return Err(TracerError::InvalidInput);
        }
        let (u1, u2, ul) = pixel_dims(
            &self.settings,
            self.sobol.as_ref(),
            self.key,
            probe_index,
            sample,
        )?;
        let position = probe.position_m();
        let normal = probe.normal();
        let radius = u1.sqrt();
        let phi = 2.0 * PI * u2;
        let direction = to_world_all_sphere(
            normal,
            [
                radius * det::cos(phi),
                radius * det::sin(phi),
                (1.0 - u1).max(0.0).sqrt(),
            ],
        );
        // Advancing along the sampled direction keeps the geometric ray and
        // physical receiver-to-hit segment collinear. A normal offset would
        // change that segment's direction without changing the cosine PDF.
        let origin = position.offset(direction.scale(self.fixtures.options.ray_epsilon_m));
        let offset = origin.delta_from(position);
        if [origin.x, origin.y, origin.z]
            .iter()
            .any(|value| !value.is_finite())
            || !offset.dot(normal).is_finite()
            || offset.dot(normal) <= 0.0
        {
            // An unrepresentable offset must not silently reuse the unshifted
            // receiver or invent a larger physical exclusion distance.
            return Err(TracerError::InvalidInput);
        }
        let previous_bsdf = if self.fixtures.areas.has_emission() {
            let pdf = normal.dot(direction).max(0.0) / PI;
            Some(PreviousBsdf {
                pdf: DirectionalPdfPair::continuous(pdf, pdf).ok_or(TracerError::InvalidInput)?,
                // The primary area-light proposal does not compete with the
                // independently integrated physical volume-emission term.
                volume_emission_competitor: false,
                opaque_source_geometric_normal: Some(normal),
                smooth_slab: None,
            })
        } else {
            None
        };
        let state = SpectralPathState {
            ray: Ray {
                origin,
                dir: direction,
            },
            throughput: [1.0; PACKET],
            previous_bsdf,
            prev_origin: position,
            // Only geometry uses the offset. Medium transport starts at the
            // exact receiver, retaining even cells thinner than ray epsilon.
            segment_origin: position,
            medium_stack: MediumStack::new(),
            rng: PathRng {
                pixel: probe_index,
                sample,
                dim: 1,
                key: self.key,
                qmc_depth: 0,
                qmc_values: [0.0; PATH_QMC_DIMENSIONS],
                qmc_valid: false,
            },
            // The receiver is a measurement, not a scattering vertex. Depth
            // one therefore includes one actual surface or volume collision.
            next_depth: 0,
            active_lane: None,
            remaining_smooth_dielectric_splits: SMOOTH_DIELECTRIC_SPLIT_BUDGET,
        };
        let render_cx = cx.with_stream_seed(self.settings.seed);
        let area_xyz = if self.fixtures.areas.has_emission() {
            let wavelengths = wavelength_packet(ul);
            let direct = self.fixtures.areas.direct_contribution(
                self.scene,
                &self.fixtures,
                &render_cx,
                &self.settings,
                0,
                self.medium,
                position,
                area::Receiver::Probe { normal },
                &wavelengths,
                &state.rng,
            )?;
            Some(packet_radiance_to_xyz(
                direct,
                &wavelengths,
                LUMINOUS_EFFICACY_LM_W,
            ))
        } else {
            None
        };
        let result = trace_path_with_emitters(
            self.scene,
            &self.lighting,
            &render_cx,
            &self.settings,
            LUMINOUS_EFFICACY_LM_W,
            probe_index,
            sample,
            u1,
            u2,
            ul,
            None,
            CameraPath::Legacy,
            false,
            self.path_sobol.as_ref(),
            Some(state),
            self.medium,
            Some(&self.fixtures),
        )?;
        // cos(theta)/pdf = pi for the sampled receiving hemisphere. This is
        // the irradiance integral, with no receiver BRDF or additional cosine.
        let mut xyz = result.xyz.map(|value| PI * value);
        if let Some(area_xyz) = area_xyz {
            for channel in 0..3 {
                xyz[channel] += area_xyz[channel];
            }
        }
        if xyz.iter().any(|value| !value.is_finite()) {
            return Err(TracerError::InvalidInput);
        }
        cx.checkpoint()?;
        Ok(xyz)
    }
}

/// Shared calibrated packet for incoming paths and probe-primary area NEE.
/// Its 20-bin lane spacing retains the source's represented photopic integral.
pub(super) fn wavelength_packet(ul: f64) -> [f64; PACKET] {
    let hero = crate::photometry::source::wavelength_nm(
        ((ul * crate::photometry::source::SPECTRAL_BINS as f64) as usize)
            .min(crate::photometry::source::SPECTRAL_BINS - 1),
    );
    hero_wavelength_packet::<PACKET>(hero, LAMBDA_MIN, LAMBDA_MAX)
}

impl FixtureLighting<'_> {
    /// Admit a transport segment, including exclusion of singular point-source
    /// lines when its sampled volume collisions could connect to those points.
    pub(super) fn admit_segment(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        direction: Vec3,
        surface: Option<Point3>,
    ) -> Result<f64, TracerError> {
        let segment = self.physical_segment(origin, direction, surface)?;
        if self.scattering
            && let Some((start, end)) = segment.medium_interval
        {
            for (index, source) in self.sources.iter().enumerate() {
                cx.checkpoint()?;
                if source.flux_lm() == 0.0 {
                    continue;
                }
                let parameter = source
                    .position_m()
                    .delta_from(origin)
                    .dot(segment.direction);
                if !parameter.is_finite() {
                    return Err(connection_error(
                        index,
                        "unrepresentable ray/source distance",
                    ));
                }
                let closest = origin.offset(segment.direction.scale(parameter.clamp(start, end)));
                let delta = source.position_m().delta_from(closest);
                let distance = if delta == Vec3::new(0.0, 0.0, 0.0) {
                    0.0
                } else {
                    norm(delta)?
                };
                if distance <= self.options.ray_epsilon_m {
                    return Err(connection_error(
                        index,
                        "ray enters point-source exclusion radius",
                    ));
                }
            }
        }
        Ok(segment.distance)
    }

    /// A thermal source line integral has no sampled fixture connection along
    /// the segment, so only physical geometry and range admission apply.
    pub(super) fn admit_source_segment(
        &self,
        origin: Point3,
        direction: Vec3,
        surface: Option<Point3>,
    ) -> Result<(Vec3, f64), TracerError> {
        let segment = self.physical_segment(origin, direction, surface)?;
        Ok((segment.direction, segment.distance))
    }

    /// Preserve the actual surface endpoint, or the finite grid exit for a
    /// geometric miss. Empty vacuum tails carry no emission or boundary light.
    fn physical_segment(
        &self,
        origin: Point3,
        direction: Vec3,
        surface: Option<Point3>,
    ) -> Result<FixtureSegment, TracerError> {
        if [
            origin.x,
            origin.y,
            origin.z,
            direction.x,
            direction.y,
            direction.z,
        ]
        .iter()
        .any(|value| !value.is_finite())
            || (direction.dot(direction) - 1.0).abs() > 1.0e-10
        {
            return Err(TracerError::InvalidInput);
        }
        let (direction, surface_distance) = if let Some(endpoint) = surface {
            let delta = endpoint.delta_from(origin);
            let distance = norm(delta)?;
            (delta.scale(1.0 / distance), Some(distance))
        } else {
            (direction, None)
        };
        let interval = self
            .medium_bounds_m
            .map(|bounds| clip_bounds(bounds, origin, direction, surface_distance))
            .transpose()?
            .flatten();
        let distance = surface_distance
            .or_else(|| interval.map(|(_, end)| end))
            .unwrap_or(0.0);
        if !distance.is_finite() || distance > self.options.max_distance_m {
            return Err(TracerError::FixturePathRange);
        }
        Ok(FixtureSegment {
            direction,
            distance,
            medium_interval: interval,
        })
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn surface_contribution(
    scene: &Scene,
    fixtures: &FixtureLighting<'_>,
    cx: &Cx<'_>,
    medium: Option<PathMedium<'_, '_>>,
    material: &Material,
    point: Point3,
    normal: Vec3,
    outgoing: Vec3,
    brushed: Option<&crate::surface_detail::ShadingFrame>,
    wavelengths: &[f64; PACKET],
) -> Result<[f64; PACKET], TracerError> {
    contribution(
        scene,
        fixtures,
        cx,
        medium,
        point,
        wavelengths,
        |incoming, lambda| {
            let cosine = normal.dot(incoming).max(0.0);
            if cosine == 0.0 {
                Ok(0.0)
            } else {
                Ok(cosine
                    * opaque_bsdf_eval(
                        material, normal, outgoing, incoming, lambda, None, brushed,
                    )?)
            }
        },
    )
}

#[allow(clippy::too_many_arguments)]
pub(super) fn volume_contribution(
    scene: &Scene,
    fixtures: &FixtureLighting<'_>,
    cx: &Cx<'_>,
    medium: PathMedium<'_, '_>,
    phase: PhaseFunction,
    point: Point3,
    incoming_travel: Vec3,
    wavelengths: &[f64; PACKET],
) -> Result<[f64; PACKET], TracerError> {
    contribution(
        scene,
        fixtures,
        cx,
        Some(medium),
        point,
        wavelengths,
        |towards_source, _| Ok(phase.eval_cos(incoming_travel.dot(towards_source))?),
    )
}

fn contribution(
    scene: &Scene,
    fixtures: &FixtureLighting<'_>,
    cx: &Cx<'_>,
    medium: Option<PathMedium<'_, '_>>,
    point: Point3,
    wavelengths: &[f64; PACKET],
    mut response: impl FnMut(Vec3, f64) -> Result<f64, TracerError>,
) -> Result<[f64; PACKET], TracerError> {
    let mut result = [0.0; PACKET];
    let mut correction = [0.0; PACKET];
    for (index, source) in fixtures.sources.iter().enumerate() {
        cx.checkpoint()?;
        if source.flux_lm() == 0.0 {
            continue;
        }
        let displacement = source.position_m().delta_from(point);
        let distance =
            norm(displacement).map_err(|_| connection_error(index, "invalid source distance"))?;
        if distance <= fixtures.options.ray_epsilon_m {
            return Err(connection_error(
                index,
                "receiver enters point-source exclusion radius",
            ));
        }
        if distance > fixtures.options.max_distance_m {
            return Err(connection_error(index, "source exceeds distance budget"));
        }
        let towards_source = displacement.scale(1.0 / distance);
        let from_source = towards_source.scale(-1.0);
        let intensity = source
            .intensity_cd(from_source)
            .map_err(TracerError::FixtureSource)?;
        if intensity == 0.0 {
            continue;
        }
        let ray = Ray {
            origin: source.position_m(),
            dir: from_source,
        };
        if photometric_intersect(
            &scene.primitives,
            cx,
            &ray,
            distance - fixtures.options.ray_epsilon_m,
        )?
        .is_some()
        {
            continue;
        }
        let transmittance = if let Some(medium) = medium {
            medium.transmittance(cx, ray.origin, ray.dir, distance, wavelengths)?
        } else {
            [1.0; PACKET]
        };
        let incident_lux = intensity / distance / distance;
        for (lane, &lambda) in wavelengths.iter().enumerate() {
            let value = incident_lux
                * source.spectrum().value_w_per_lumen_nm(lambda)
                * response(towards_source, lambda)?
                * transmittance[lane];
            if !value.is_finite() || value < 0.0 {
                return Err(connection_error(index, "nonfinite spectral contribution"));
            }
            let sum = result[lane] + value;
            correction[lane] += if result[lane].abs() >= value.abs() {
                (result[lane] - sum) + value
            } else {
                (value - sum) + result[lane]
            };
            result[lane] = sum;
        }
    }
    for lane in 0..PACKET {
        result[lane] += correction[lane];
        if !result[lane].is_finite() || result[lane] < 0.0 {
            return Err(TracerError::InvalidInput);
        }
    }
    Ok(result)
}

fn connection_error(source: usize, reason: &'static str) -> TracerError {
    TracerError::FixtureConnection { source, reason }
}

fn norm(vector: Vec3) -> Result<f64, TracerError> {
    let scale = vector.x.abs().max(vector.y.abs()).max(vector.z.abs());
    if !scale.is_finite() || scale == 0.0 {
        return Err(TracerError::InvalidInput);
    }
    let value = Vec3::new(vector.x / scale, vector.y / scale, vector.z / scale).norm() * scale;
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(TracerError::InvalidInput)
    }
}

fn clip_bounds(
    (lo, hi): ([f64; 3], [f64; 3]),
    origin: Point3,
    direction: Vec3,
    endpoint: Option<f64>,
) -> Result<Option<(f64, f64)>, TracerError> {
    let o = [origin.x, origin.y, origin.z];
    let d = [direction.x, direction.y, direction.z];
    let (mut enter, mut exit) = (0.0_f64, endpoint.unwrap_or(f64::INFINITY));
    for axis in 0..3 {
        if d[axis] == 0.0 {
            if o[axis] < lo[axis] || o[axis] >= hi[axis] {
                return Ok(None);
            }
        } else {
            let a = (lo[axis] - o[axis]) / d[axis];
            let b = (hi[axis] - o[axis]) / d[axis];
            enter = enter.max(a.min(b));
            exit = exit.min(a.max(b));
        }
    }
    if exit <= enter {
        return Ok(None);
    }
    if !enter.is_finite() || !exit.is_finite() {
        return Err(TracerError::FixturePathRange);
    }
    Ok(Some((enter, exit)))
}
