//! Calibrated point fixtures feed the existing spectral surface/volume paths.
//!
//! Every admitted fixture is evaluated at each scattering vertex. These delta
//! position sources have no competing continuous BSDF/phase hit technique, so
//! their direct contribution has unit MIS weight. The established area-light
//! mixture and its random-number consumption are unchanged.

use super::*;
use crate::camera::{Aperture, PhysicalCamera};
use crate::photometry::medium::HomogeneousMedium;
use crate::photometry::source::{LUMINOUS_EFFICACY_LM_W, PhotometricSource};
use crate::photometry::{PhotometricScene, PhotometryError, TransportOptions};
use crate::volumes::transport::AbsorbingGrid;
use crate::volumes::transport::scattering::{PhaseFunction, ScatteringGrid};

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
    options: TransportOptions,
    medium_bounds_m: Option<([f64; 3], [f64; 3])>,
    scattering: bool,
}

/// One admitted scene, fixture set, medium and sample policy. All retained
/// assets are borrowed, including the original mesh BVHs and physical cells.
pub(crate) struct FixtureTransport<'scene, 'field> {
    scene: &'scene Scene,
    fixtures: FixtureLighting<'scene>,
    lighting: AdmittedLighting<'scene>,
    medium: Option<PathMedium<'scene, 'field>>,
    settings: Settings,
    key: [u32; 2],
    sobol: Option<Sobol>,
    path_sobol: Option<Sobol>,
}

impl<'scene, 'field> FixtureTransport<'scene, 'field> {
    pub(crate) fn try_new(
        scene: &'scene Scene,
        sources: &'scene [PhotometricSource],
        medium: FixtureMedium<'scene, 'field>,
        options: TransportOptions,
        settings: Settings,
        cx: &Cx<'_>,
    ) -> Result<Self, PhotometryError> {
        cx.checkpoint()?;
        if settings.spp == 0
            || settings.max_depth == 0
            || settings.strategy == DirectStrategy::BsdfOnly
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
        PhotometricScene::try_new(
            &scene.primitives,
            sources,
            HomogeneousMedium::vacuum(),
            options,
            cx,
        )?;
        PhysicalCamera::try_legacy_compatible(
            scene.camera.eye,
            scene.camera.forward,
            scene.camera.up,
            scene.camera.half_tan,
            1.0,
            Aperture::try_circular(0.0)?,
        )?;
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
            fixtures: FixtureLighting {
                sources,
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
}

impl FixtureLighting<'_> {
    /// Preserve the actual surface endpoint, or the finite grid exit for a
    /// geometric miss. Empty vacuum tails carry no emission or boundary light.
    pub(super) fn admit_segment(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        direction: Vec3,
        surface: Option<Point3>,
    ) -> Result<f64, TracerError> {
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
        if self.scattering
            && let Some((start, end)) = interval
        {
            for (index, source) in self.sources.iter().enumerate() {
                cx.checkpoint()?;
                if source.flux_lm() == 0.0 {
                    continue;
                }
                let parameter = source.position_m().delta_from(origin).dot(direction);
                if !parameter.is_finite() {
                    return Err(connection_error(
                        index,
                        "unrepresentable ray/source distance",
                    ));
                }
                let closest = origin.offset(direction.scale(parameter.clamp(start, end)));
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
        Ok(distance)
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
                    * opaque_bsdf_eval(material, normal, outgoing, incoming, lambda, None, None)?)
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
