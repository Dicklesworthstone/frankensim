//! Time-correct direct proposals for accepted rigid thermal trajectories.
//!
//! The source set is disjoint from the static light mixture. A body is selected
//! by cached area-times-visible-luminance, then its immutable local hierarchy
//! proposes an actual face endpoint. Proper-rigid transforms preserve area and
//! solid angle, so the same local-space PDF applies in world space. No triangle
//! positions, source spectra, or acceleration structures are rebuilt per ray.

use std::collections::BTreeMap;

use fs_math::det;

use super::*;

/// A small uniform component preserves every positive body's proposal support
/// when physical temperatures or sizes span more than binary64's dynamic range.
const UNIFORM_BODY_FRACTION: f64 = 1.0 / 16.0;

struct MovingSource<'scene> {
    primitive_index: usize,
    instance: &'scene AnimatedGeometryInstance,
    local: AdmittedLighting<'scene>,
    log_power: f64,
}

/// Prepared once per render admission, borrowed unchanged by all path workers.
/// Body selection and source lookup are logarithmic in emitting-body count;
/// each body's existing light hierarchy handles its own triangle count.
pub(in crate::tracer) struct MovingThermalLights<'scene> {
    sources: Vec<MovingSource<'scene>>,
    /// Stored interval differences define both selection and reverse density.
    cdf: Vec<f64>,
    source_for_primitive: BTreeMap<usize, usize>,
}

impl<'scene> MovingThermalLights<'scene> {
    pub(in crate::tracer) fn try_new(
        scene: &'scene Scene,
        cx: &Cx<'_>,
    ) -> Result<Self, TracerError> {
        // Stable object identities, rather than caller insertion order, fix
        // the body-selection stream. Scene admission checks their uniqueness.
        let mut ordered = BTreeMap::new();
        for (primitive_index, primitive) in scene.primitives.iter().enumerate() {
            cx.checkpoint()?;
            let Shape::AnimatedInstance(instance) = &primitive.shape else {
                continue;
            };
            let Some(surface) = primitive.material.thermal_surface() else {
                continue;
            };
            if !surface.is_emitting() {
                continue;
            }
            let SharedGeometry::Mesh(mesh) = instance.geometry() else {
                continue;
            };
            let mut triangles = Vec::new();
            triangles
                .try_reserve_exact(mesh.triangles.len())
                .map_err(|_| LightingError::EnvironmentTooLarge)?;
            for triangle in &mesh.triangles {
                cx.checkpoint()?;
                let mut vertices = [Point3::new(0.0, 0.0, 0.0); 3];
                for (vertex, &index) in vertices.iter_mut().zip(triangle) {
                    let point = mesh
                        .vertices
                        .get(index as usize)
                        .ok_or(TracerError::InvalidInput)?;
                    *vertex = Point3::new(point[0], point[1], point[2]);
                }
                triangles.push(ThermalTriangleLight {
                    vertices,
                    primitive_index,
                    surface,
                });
            }
            let local =
                AdmittedLighting::with_thermal_surface_emission(cx, &[], None, triangles, true)?;
            let Some(log_power) = local.thermal_log_power(cx)? else {
                // The gray source's entire visible spectrum can underflow.
                // No total-spectrum power is inferred from a spectral table.
                continue;
            };
            if ordered
                .insert(
                    instance.object_id(),
                    MovingSource {
                        primitive_index,
                        instance,
                        local,
                        log_power,
                    },
                )
                .is_some()
            {
                return Err(TracerError::InvalidInstance);
            }
        }
        let mut result = Self {
            sources: Vec::new(),
            cdf: Vec::new(),
            source_for_primitive: BTreeMap::new(),
        };
        if ordered.is_empty() {
            return Ok(result);
        }
        result
            .sources
            .try_reserve_exact(ordered.len())
            .map_err(|_| LightingError::EnvironmentTooLarge)?;
        result
            .cdf
            .try_reserve_exact(
                ordered
                    .len()
                    .checked_add(1)
                    .ok_or(LightingError::EnvironmentTooLarge)?,
            )
            .map_err(|_| LightingError::EnvironmentTooLarge)?;
        let mut maximum = f64::NEG_INFINITY;
        for source in ordered.into_values() {
            cx.checkpoint()?;
            maximum = maximum.max(source.log_power);
            result
                .source_for_primitive
                .insert(source.primitive_index, result.sources.len());
            result.sources.push(source);
        }
        let mut total = 0.0;
        for source in &result.sources {
            cx.checkpoint()?;
            total += det::exp(source.log_power - maximum);
        }
        if !(total > 0.0 && total.is_finite()) {
            return Err(TracerError::InvalidInput);
        }
        let uniform = UNIFORM_BODY_FRACTION / result.sources.len() as f64;
        let mut cumulative = 0.0;
        result.cdf.push(0.0);
        for source in &result.sources {
            cx.checkpoint()?;
            cumulative += uniform
                + (1.0 - UNIFORM_BODY_FRACTION) * det::exp(source.log_power - maximum) / total;
            result.cdf.push(cumulative);
        }
        let last = result.cdf.len() - 1;
        result.cdf[last] = 1.0;
        for pair in result.cdf.windows(2) {
            cx.checkpoint()?;
            if pair[1] <= pair[0] {
                return Err(TracerError::InvalidInput);
            }
        }
        cx.checkpoint()?;
        Ok(result)
    }

    pub(in crate::tracer) fn is_empty(&self) -> bool {
        self.sources.is_empty()
    }

    /// The first variate selects a body; the other two retain the existing
    /// triangle-selection residual and area parameterization in its local frame.
    pub(in crate::tracer) fn sample(
        &self,
        cx: &Cx<'_>,
        ray_time: Option<&PathTime>,
        origin: Point3,
        variates: [f64; 3],
    ) -> Result<Option<ThermalLightSample<'scene>>, TracerError> {
        cx.checkpoint()?;
        if self.is_empty() {
            return Ok(None);
        }
        if variates
            .iter()
            .any(|&value| !value.is_finite() || !(0.0..1.0).contains(&value))
        {
            return Err(TracerError::InvalidInput);
        }
        let index = self.cdf.partition_point(|&edge| edge <= variates[0]) - 1;
        let source = &self.sources[index];
        let transform = source.transform(ray_time)?;
        let inverse = transform.inverse()?;
        let local_origin = inverse.transform_point(origin);
        if !finite_point(local_origin) {
            return Err(TracerError::InvalidInstance);
        }
        let Some(sample) = source.local.sample(local_origin, variates[1], variates[2]) else {
            return Ok(None);
        };
        let LightSample::ThermalSurface(mut sample) = sample else {
            return Err(TracerError::InvalidInput);
        };
        sample.point = transform.transform_point(sample.point);
        sample.normal = transform.transform_vector(sample.normal);
        // Classify the rounded world endpoint through the same inverse and
        // uncertainty envelope as a BSDF hit. Forward local-only support can
        // disagree after a large translation, especially at shared mesh edges.
        sample.pdf_solid_angle = self.probability(index)
            * source.world_endpoint_pdf(transform, inverse, local_origin, sample.point)?;
        if !finite_point(sample.point)
            || !sample.pdf_solid_angle.is_finite()
            || sample.pdf_solid_angle <= 0.0
        {
            return Err(TracerError::InvalidInstance);
        }
        cx.checkpoint()?;
        Ok(Some(sample))
    }

    /// Reverse density in exactly the source proposal and path-time measure.
    /// A static or unsupported emitter has no competitor from this technique.
    pub(in crate::tracer) fn hit_pdf(
        &self,
        cx: &Cx<'_>,
        ray_time: Option<&PathTime>,
        primitive_index: usize,
        origin: Point3,
        point: Point3,
    ) -> Result<f64, TracerError> {
        cx.checkpoint()?;
        let Some(&index) = self.source_for_primitive.get(&primitive_index) else {
            return Ok(0.0);
        };
        let source = &self.sources[index];
        let transform = source.transform(ray_time)?;
        let inverse = transform.inverse()?;
        let local_origin = inverse.transform_point(origin);
        if !finite_point(local_origin) {
            return Err(TracerError::InvalidInstance);
        }
        let pdf = self.probability(index)
            * source.world_endpoint_pdf(transform, inverse, local_origin, point)?;
        if !pdf.is_finite() || pdf < 0.0 {
            return Err(TracerError::InvalidInput);
        }
        cx.checkpoint()?;
        Ok(pdf)
    }

    fn probability(&self, index: usize) -> f64 {
        self.cdf[index + 1] - self.cdf[index]
    }
}

impl MovingSource<'_> {
    fn world_endpoint_pdf(
        &self,
        transform: RigidTransform,
        inverse: RigidTransform,
        local_origin: Point3,
        point: Point3,
    ) -> Result<f64, TracerError> {
        let local_point = inverse.transform_point(point);
        if !finite_point(point) || !finite_point(local_point) {
            return Err(TracerError::InvalidInstance);
        }
        // Both quaternion rotations and the inverse translation are rounded.
        // Their fixed operation count gives this conservative absolute error
        // envelope in local metres. The local triangle's own arithmetic error
        // remains separate. Use world and translation magnitudes so cancellation
        // cannot hide the uncertainty after a world-to-body transformation.
        let scale = [point.x, point.y, point.z]
            .into_iter()
            .chain(transform.translation_m())
            .chain(inverse.translation_m())
            .map(f64::abs)
            .fold(0.0_f64, f64::max);
        let support_tolerance_m = (256.0 * f64::EPSILON) * scale;
        Ok(self.local.thermal_mixture_pdf_with_tolerance(
            self.primitive_index,
            local_origin,
            local_point,
            support_tolerance_m,
        ))
    }

    fn transform(&self, ray_time: Option<&PathTime>) -> Result<RigidTransform, TracerError> {
        let time = ray_time.ok_or(TracerError::MissingRayTime)?;
        if let Some(cached) = time
            .cached_animated
            .iter()
            .flatten()
            .find(|cached| cached.primitive_index == self.primitive_index)
        {
            return Ok(cached.instance.transform());
        }
        Ok(self
            .instance
            .trajectory()
            .evaluate(time.interval.time_at(time.normalized))?
            .transform)
    }
}

fn finite_point(point: Point3) -> bool {
    [point.x, point.y, point.z].into_iter().all(f64::is_finite)
}

#[cfg(test)]
mod tests;
