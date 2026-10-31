//! Physical thermal emission in the existing surface-path estimator.
//!
//! The proposal uses actual mesh faces, including rigidly placed body
//! boundaries. Moving meshes retain one local hierarchy and evaluate rigid
//! placement at the same shutter time as every visibility ray. Thermal charts
//! retain camera-path support.
//! No spectrum is converted through RGB, and observing light never updates
//! the thermal owner's stored energy.

use super::*;
use crate::lighting::{ThermalLightSample, ThermalTriangleLight};

mod moving;
pub(super) use moving::MovingThermalLights;

pub(super) fn triangles<'scene>(
    scene: &'scene Scene,
    cx: &Cx<'_>,
) -> Result<Vec<ThermalTriangleLight<'scene>>, TracerError> {
    let mut lights = Vec::new();
    for (primitive_index, primitive) in scene.primitives.iter().enumerate() {
        cx.checkpoint()?;
        let Some(surface) = primitive.material.thermal_surface() else {
            continue;
        };
        if !surface.is_emitting() {
            continue;
        }
        let (mesh, placement) = match &primitive.shape {
            Shape::Mesh(mesh) => (mesh, None),
            Shape::Instance(instance) => match instance.geometry() {
                SharedGeometry::Mesh(mesh) => (mesh.as_ref(), Some(instance.transform())),
                SharedGeometry::Chart(_) => continue,
            },
            Shape::Chart(_) | Shape::AnimatedInstance(_) => continue,
        };
        lights
            .try_reserve(mesh.triangles.len())
            .map_err(|_| TracerError::Lighting(LightingError::EnvironmentTooLarge))?;
        for triangle in &mesh.triangles {
            cx.checkpoint()?;
            let mut vertices = [Point3::new(0.0, 0.0, 0.0); 3];
            for (corner, &index) in vertices.iter_mut().zip(triangle) {
                let vertex = mesh
                    .vertices
                    .get(index as usize)
                    .ok_or(TracerError::InvalidInput)?;
                let local = Point3::new(vertex[0], vertex[1], vertex[2]);
                *corner = placement.map_or(local, |transform| transform.transform_point(local));
            }
            lights.push(ThermalTriangleLight {
                vertices,
                primitive_index,
                surface,
            });
        }
    }
    Ok(lights)
}

pub(super) fn hit_weight(
    lighting: &AdmittedLighting<'_>,
    strategy: DirectStrategy,
    previous: Option<PreviousBsdf>,
    primitive_index: usize,
    origin: Point3,
    point: Point3,
) -> f64 {
    let Some(previous) = previous else { return 1.0 };
    if previous.pdf.is_delta() || previous.smooth_slab.is_some() {
        return 1.0;
    }
    let pdf = lighting.thermal_mixture_pdf(primitive_index, origin, point);
    if pdf > 0.0 {
        emissive_hit_weight(strategy, Some(previous), Some(pdf))
    } else {
        // An unsupported shape or direction has no competing proposal,
        // including in the NEE-only diagnostic mode.
        1.0
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn direct_contribution(
    scene: &Scene,
    cx: &Cx<'_>,
    ray_time: Option<&PathTime>,
    settings: &Settings,
    depth: u32,
    light: ThermalLightSample<'_>,
    material: &Material,
    point: Point3,
    frame: SurfaceFrame,
    outgoing: Vec3,
    media: &MediumStack<'_>,
    boundary: Option<BoundaryMedia<'_>>,
    brushed: Option<&crate::surface_detail::ShadingFrame>,
    wavelengths: &[f64; PACKET],
    #[cfg(feature = "volumes")] medium: Option<PathMedium<'_, '_>>,
) -> Result<[f64; PACKET], TracerError> {
    let mut result = [0.0; PACKET];
    let displacement = light.point.delta_from(point);
    let distance = displacement.norm();
    if !(distance > 0.0 && distance.is_finite() && light.pdf_solid_angle > 0.0) {
        return Ok(result);
    }
    let incoming = displacement.scale(1.0 / distance);
    let cosine = frame.oriented.dot(incoming);
    if cosine == 0.0 || (boundary.is_none() && cosine < 0.0) {
        return Ok(result);
    }
    let shadow_origin = dielectric_spawn_origin(point, frame.geometric, incoming);
    if !visible_endpoint(
        scene,
        cx,
        ray_time,
        point,
        shadow_origin,
        light.point,
        light.primitive_index,
    )? {
        return Ok(result);
    }
    let shadow_medium = boundary.map_or_else(
        || media.last().map(|entry| *entry.glass),
        |boundary| {
            if cosine > 0.0 {
                boundary.incident
            } else {
                boundary.transmitted
            }
        },
    );
    if shadow_medium.is_some() {
        return Err(TracerError::ThermalSurfaceMediumUnsupported);
    }
    #[cfg(feature = "volumes")]
    let volume_transmittance = if let Some(medium) = medium {
        medium.transmittance(cx, point, incoming, distance, wavelengths)?
    } else {
        [1.0; PACKET]
    };
    for (lane, &wavelength) in wavelengths.iter().enumerate() {
        let (value, competing_pdf) = match *material {
            Material::Dielectric { surface, .. } => {
                let Some(alpha) = surface.roughness_alpha() else {
                    continue;
                };
                let boundary = boundary.ok_or(TracerError::InvalidInput)?;
                let evaluation = surface.interface(wavelength, boundary.coating_reversed())?.evaluate_rough(
                    frame.oriented,
                    outgoing,
                    incoming,
                    medium_ior(boundary.incident, wavelength)?,
                    medium_ior(boundary.transmitted, wavelength)?,
                    alpha,
                )?;
                (evaluation.value, evaluation.pdf)
            }
            _ => (
                opaque_bsdf_eval(
                    material,
                    frame.oriented,
                    outgoing,
                    incoming,
                    wavelength,
                    None,
                    brushed,
                )?,
                bsdf_pdf(material, frame.oriented, outgoing, incoming, brushed),
            ),
        };
        if value <= 0.0 {
            continue;
        }
        #[cfg(feature = "volumes")]
        let competing_pdf = competing_pdf
            * medium.map_or(1.0, |medium| {
                medium.escape_probability(&volume_transmittance)
            });
        let weight = if depth + 1 == settings.max_depth {
            // No BSDF ray will be traced beyond the final admitted vertex.
            1.0
        } else {
            match settings.strategy {
                DirectStrategy::Mis => {
                    balance_heuristic(1, light.pdf_solid_angle, 1, competing_pdf)
                }
                DirectStrategy::PowerMis => {
                    power_heuristic(1, light.pdf_solid_angle, 1, competing_pdf)
                }
                DirectStrategy::NeeOnly => 1.0,
                DirectStrategy::BsdfOnly => return Err(TracerError::InvalidInput),
            }
        };
        let emitted = light
            .surface
            .radiance_w_m2_sr_nm(wavelength)
            .map_err(|_| TracerError::InvalidInput)?;
        result[lane] = value * cosine.abs() * emitted / light.pdf_solid_angle * weight;
        #[cfg(feature = "volumes")]
        {
            result[lane] *= volume_transmittance[lane];
        }
        if !result[lane].is_finite() || result[lane] < 0.0 {
            return Err(TracerError::InvalidInput);
        }
    }
    Ok(result)
}

/// Trace an exact sampled endpoint, rejecting hidden faces of the same body.
#[allow(clippy::too_many_arguments)]
pub(super) fn visible_endpoint(
    scene: &Scene,
    cx: &Cx<'_>,
    ray_time: Option<&PathTime>,
    source: Point3,
    shadow_origin: Point3,
    target: Point3,
    primitive_index: usize,
) -> Result<bool, TracerError> {
    let displacement = target.delta_from(shadow_origin);
    let distance = displacement.norm();
    if !(distance > 0.0 && distance.is_finite()) {
        return Ok(false);
    }
    let ray = Ray {
        origin: shadow_origin,
        dir: displacement.scale(1.0 / distance),
    };
    let Some(hit) = intersect(scene, cx, &ray, ray_time)? else {
        return Ok(false);
    };
    let coordinate_scale = [source.x, source.y, source.z, target.x, target.y, target.z]
        .into_iter()
        .map(f64::abs)
        .fold(target.delta_from(source).norm(), f64::max);
    // Aimed-at-target origin offsets cause no endpoint bias. A metric
    // RAY_EPS allowance here would leak a nearby hidden face of the same body.
    let tolerance = 128.0 * f64::EPSILON * coordinate_scale;
    Ok(hit.primitive_index == primitive_index
        && hit.hit.point.delta_from(target).norm() <= tolerance)
}
