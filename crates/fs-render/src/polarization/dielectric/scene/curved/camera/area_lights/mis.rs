//! Reciprocal area/direction densities for native finite-source camera MIS.
//!
//! For a regular source-to-receiver image, the differentiated bundle gives
//! J = dA_receiver / dOmega_source. Optical reciprocity gives
//! dA_source / dOmega_receiver = (n_receiver/n_source)^2
//! * cos_receiver/cos_source * J. This uses the actual endpoint phase indices
//! and normals, including continuous refraction. See Pediredla et al. (2020),
//! "Path Tracing Estimators for Refractive Radiative Transfer", and PBRT 4e
//! section 13.4 for complementary emitter-hit/light-sampling weights.
//!
//! One area draw enumerates each retained image once. Its per-image density
//! is p_area times that reciprocal Jacobian; there is no image-count divisor.
//! The competing diffuse continuation uses cos_receiver/pi times the exact
//! reverse native scalar probabilities of its selected smooth interfaces.
//! Absorption and Mueller power belong to radiance, not these proposals.
//! These weights do not prove that finite authored seeds discover every image.
//!
//! Participating consumers use the same deterministic angular/branch scores
//! at both ends of a connection, without scattering survival in the heuristic.
//! The full continuation density would also contain no-collision survival;
//! these scores do not claim to be that marginalized path density. A shared
//! partition of unity still preserves the two estimators' expectation (PBRT 4e
//! section 2.2.3). Source survival remains an independent radiance estimator.
//! Inserting a noisy ratio-tracking estimate into this nonlinear heuristic
//! would correlate its weight with source attenuation and bias the result.

use super::{AreaEmitterError, AreaEmitterSample, CurvedAreaEmitter, SourceSupport, cross, delta, norm};
use crate::gradient_index::bundle::ProjectionOptions;
use crate::instances::SharedGeometry;
use crate::polarization::dielectric::scene::curved::CurvedContact;
use crate::polarization::dielectric::scene::curved::bundle::connection::{
    BundleConnection, BundleConnectionTarget,
};
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_math::det;

/// Complementary power-heuristic weights for one regular optical image.
/// Angular/branch scores use the same receiver solid-angle measure. They omit
/// participating free-flight survival. An image absent from a finite light
/// catalog must instead retain emitter-hit weight one.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AreaImageMisWeights {
    light: f64,
    continuation: f64,
    light_log_pdf: f64,
    continuation_log_pdf: f64,
}

impl AreaImageMisWeights {
    pub const fn light_weight(&self) -> f64 {
        self.light
    }

    pub const fn continuation_weight(&self) -> f64 {
        self.continuation
    }

    pub const fn light_log_pdf_solid_angle(&self) -> f64 {
        self.light_log_pdf
    }

    /// Log of angular PDF times reverse smooth-interface branch probability.
    /// In participating media this is a deterministic weighting score, not
    /// the full suffix probability including no-collision survival.
    pub const fn continuation_log_pdf_solid_angle(&self) -> f64 {
        self.continuation_log_pdf
    }
}

impl CurvedAreaEmitter<'_, '_, '_> {
    /// Evaluate the two proposals on an existing native diffuse-receiver
    /// connection. This does not discover images or establish membership of
    /// an independently sampled camera path in a light's finite seed catalog.
    pub fn mis_weights(
        &self,
        sample: AreaEmitterSample<'_>,
        connection: &BundleConnection,
        target: BundleConnectionTarget,
        projection: ProjectionOptions,
    ) -> Result<AreaImageMisWeights, AreaEmitterError> {
        let path = connection.path();
        let receiver = path.last_contact().ok_or(AreaEmitterError::SourceContactMismatch)?;
        if !target.is_diffuse_receiver()
            || !path.is_diffuse_receiver()
            || receiver.object_id != target.object_id()
        {
            return Err(AreaEmitterError::SourceContactMismatch);
        }
        let cosine_receiver = -target.normal().dot(path.bundle().propagation().propagation().ray.dir);
        if !cosine_receiver.is_finite() || cosine_receiver <= 0.0 {
            return Err(AreaEmitterError::SourceContactMismatch);
        }
        let cosine_receiver = cosine_receiver.min(1.0);
        let (light_log_pdf, reverse_log_probability) = self.image_log_densities(
            sample, connection, target.normal(), cosine_receiver, projection,
        )?;
        let continuation_log_pdf = det::ln(cosine_receiver)
            - det::ln(core::f64::consts::PI) + reverse_log_probability;
        power_weights(light_log_pdf, continuation_log_pdf)
    }

    /// Deterministic angular partition at an original volume receiver.
    /// `continuation_pdf_solid_angle` must be the native angular sampler's PDF
    /// evaluated for the camera incoming direction and the reverse of this
    /// connection's final direction. For a molecular vertex use the original
    /// local `PhaseFunction::eval_cos`; polarized power is a separate factor.
    ///
    /// The physical receiving section is perpendicular to the actual arriving
    /// ray. The target's arbitrary Newton section therefore cannot alter the
    /// returned density. Neither collision rate nor scattering survival enters
    /// these scores; participating consumers must use the same deterministic
    /// weights for area sampling and a matching native emitter hit.
    pub fn volume_mis_weights(
        &self,
        sample: AreaEmitterSample<'_>,
        connection: &BundleConnection,
        target: BundleConnectionTarget,
        projection: ProjectionOptions,
        continuation_pdf_solid_angle: f64,
    ) -> Result<AreaImageMisWeights, AreaEmitterError> {
        let path = connection.path();
        let endpoint = path.bundle().propagation().propagation();
        if !path.is_volume_receiver()
            || target.volume_medium_index() != Some(path.medium_index())
            || norm(endpoint.ray.origin.delta_from(target.position_m()))
                > target.position_tolerance_m()
        {
            return Err(AreaEmitterError::SourceContactMismatch);
        }
        if !continuation_pdf_solid_angle.is_finite() || continuation_pdf_solid_angle < 0.0 {
            return Err(AreaEmitterError::NumericalRange);
        }
        let (light_log_pdf, reverse_log_probability) = self.image_log_densities(
            sample, connection, endpoint.ray.dir.scale(-1.0), 1.0, projection,
        )?;
        let continuation_log_pdf = if continuation_pdf_solid_angle == 0.0 {
            f64::NEG_INFINITY
        } else {
            det::ln(continuation_pdf_solid_angle) + reverse_log_probability
        };
        power_weights(light_log_pdf, continuation_log_pdf)
    }

    /// Reciprocal directional weights at an original continuous surface BSDF.
    /// The supplied PDF is the native camera sampler evaluated toward this
    /// source image, including its local reflection/transmission selection.
    /// It is independent of the physical polarized BSDF contraction.
    ///
    /// The actual native contact normal determines the receiving footprint
    /// and cosine, on either geometric side. The target normal is only the
    /// numerical shooting section and need not be the surface tangent. The
    /// endpoint phase index belongs to the light's arriving medium, which may
    /// differ from the camera's medium across a transmitting rough interface.
    pub fn surface_mis_weights(
        &self,
        sample: AreaEmitterSample<'_>,
        connection: &BundleConnection,
        target: BundleConnectionTarget,
        projection: ProjectionOptions,
        continuation_pdf_solid_angle: f64,
    ) -> Result<AreaImageMisWeights, AreaEmitterError> {
        let path = connection.path();
        let endpoint = path.bundle().propagation().propagation();
        let contact = path.last_contact().ok_or(AreaEmitterError::SourceContactMismatch)?;
        if !target.is_surface_receiver() || !path.is_surface_receiver()
            || contact.object_id != target.object_id()
            || norm(endpoint.ray.origin.delta_from(target.position_m()))
                > target.position_tolerance_m()
        {
            return Err(AreaEmitterError::SourceContactMismatch);
        }
        if !continuation_pdf_solid_angle.is_finite() || continuation_pdf_solid_angle < 0.0 {
            return Err(AreaEmitterError::NumericalRange);
        }
        let cosine = contact.witness.normal.dot(endpoint.ray.dir);
        if !cosine.is_finite() || cosine == 0.0 {
            return Err(AreaEmitterError::SourceContactMismatch);
        }
        let normal = if cosine > 0.0 { contact.witness.normal.scale(-1.0) }
            else { contact.witness.normal };
        let (light_log_pdf, reverse_log_probability) = self.image_log_densities(
            sample, connection, normal, cosine.abs().min(1.0), projection,
        )?;
        let continuation_log_pdf = if continuation_pdf_solid_angle == 0.0 {
            f64::NEG_INFINITY
        } else {
            det::ln(continuation_pdf_solid_angle) + reverse_log_probability
        };
        power_weights(light_log_pdf, continuation_log_pdf)
    }

    // All endpoint kinds share source ownership, reciprocal etendue and the
    // exact reverse scalar interface proposal. Only their physical receiving
    // section and local continuous angular law differ.
    fn image_log_densities(
        &self,
        sample: AreaEmitterSample<'_>,
        connection: &BundleConnection,
        receiver_normal: Vec3,
        cosine_receiver: f64,
        projection: ProjectionOptions,
    ) -> Result<(f64, f64), AreaEmitterError> {
        if !core::ptr::eq(sample.emitter, &self.emitter)
            || !core::ptr::eq(sample.instance, self.instance)
            || !core::ptr::eq(sample.support, &self.support)
        {
            return Err(AreaEmitterError::SampleOwnerMismatch);
        }
        let path = connection.path();
        let contact = path.source_contact().ok_or(AreaEmitterError::SourceContactMismatch)?;
        let bundle = path.bundle();
        let geometric = bundle.propagation();
        if contact.object_id != self.emitter.object_id()
            || contact.geometry_identity != self.instance.geometry_identity()
            || contact.frame_identity != self.instance.frame_identity()
            || geometric.initial_ray().origin != sample.position_m
            || contact.witness.normal.dot(sample.normal) < 1.0 - 1e-8
        {
            return Err(AreaEmitterError::SourceContactMismatch);
        }
        let cosine_source = contact.witness.normal.dot(geometric.initial_ray().dir).abs().min(1.0);
        if !cosine_source.is_finite() || cosine_source <= 0.0 {
            return Err(AreaEmitterError::SourceContactMismatch);
        }
        let area = geometric.project_at_endpoint(receiver_normal, projection)
            .map_err(|error| AreaEmitterError::Refractive(error.into()))?
            .area_per_solid_angle_m2();
        let source_index = geometric.propagation().initial_sample.phase_index();
        self.emitter.admit_phase_index(source_index).map_err(AreaEmitterError::Camera)?;
        let receiver_index = geometric.propagation().final_sample.phase_index();
        let reverse_log_probability = bundle.reverse_branch_log_probability();
        if !(source_index > 0.0 && source_index.is_finite()
            && receiver_index > 0.0 && receiver_index.is_finite()
            && area > 0.0 && area.is_finite())
            || reverse_log_probability.is_nan()
            || reverse_log_probability > 0.0
        {
            return Err(AreaEmitterError::NumericalRange);
        }
        let light_log_pdf = det::ln(sample.pdf_area_per_m2)
            + 2.0 * (det::ln(receiver_index) - det::ln(source_index))
            + det::ln(cosine_receiver) - det::ln(cosine_source) + det::ln(area);
        Ok((light_log_pdf, reverse_log_probability))
    }

    /// Recover the original source's area density at a native emitter hit.
    /// Mesh membership uses its original BVH; analytic charts supply their own
    /// boundary density and normal. No chart is replaced by a proxy surface.
    /// The launch position stays on the integrated camera ray and is unsnapped.
    /// An original face outside this source's declared support returns `None`;
    /// edges, unresolved membership and work limits still refuse.
    pub(super) fn sample_at_contact_if_supported(
        &self,
        cx: &Cx<'_>,
        contact: &CurvedContact,
        position_m: Point3,
    ) -> Result<Option<AreaEmitterSample<'_>>, AreaEmitterError> {
        cx.checkpoint()?;
        if contact.object_id != self.emitter.object_id()
            || contact.geometry_identity != self.instance.geometry_identity()
            || contact.frame_identity != self.instance.frame_identity()
            || ![position_m.x, position_m.y, position_m.z].into_iter().all(f64::is_finite)
        {
            return Err(AreaEmitterError::SourceContactMismatch);
        }
        let local = self.instance.transform().inverse()
            .map_err(|_| AreaEmitterError::NumericalRange)?
            .transform_point(contact.witness.surface_point);
        // The native witness is in world coordinates. Its inverse placement
        // adds roundoff even when the original contact distance was zero.
        // Include that allowance in the BVH query as well as the plane test.
        let point = contact.witness.surface_point;
        let coordinate_scale = [point.x, point.y, point.z, local.x, local.y, local.z]
            .into_iter().fold(0.0_f64, |scale, x| scale.max(x.abs()));
        let tolerance = (contact.resolution_m + 64.0 * f64::EPSILON * coordinate_scale).next_up();
        let (triangle, local_normal, pdf_area_per_m2) = match (&self.support, self.instance.geometry()) {
            (SourceSupport::Mesh(triangles), SharedGeometry::Mesh(mesh)) => {
                let triangle = mesh.area_source_face_at(cx, local, tolerance,
                    self.limits.max_lookup_work).map_err(AreaEmitterError::Mesh)?
                    .ok_or(AreaEmitterError::SourceContactMismatch)?;
                let Ok(slot) = triangles.binary_search_by_key(&triangle, |entry| entry.original_triangle)
                    else { return Ok(None); };
                let [a, b, c] = mesh.triangles[triangle].map(|index| mesh.vertices[index as usize]);
                let entry = &triangles[slot];
                let lower = if slot == 0 { 0.0 } else { triangles[slot - 1].cumulative_m2 };
                let pdf = ((entry.cumulative_m2 - lower) / self.total_area_m2) / entry.area_m2;
                (Some(triangle), cross(delta(b, a), delta(c, a)), pdf)
            }
            (SourceSupport::Chart, SharedGeometry::Chart(chart)) => {
                if self.limits.max_lookup_work == 0 {
                    return Err(AreaEmitterError::SourceWorkLimit { limit: 0 });
                }
                let density = chart.surface_area_pdf_at(cx, local, tolerance)
                    .map_err(AreaEmitterError::Chart)?;
                let length = norm(density.normal);
                if !length.is_finite() || (length - 1.0).abs() > 1e-8 {
                    return Err(AreaEmitterError::SourceContactMismatch);
                }
                (None, density.normal, density.pdf_area_per_m2)
            }
            _ => return Err(AreaEmitterError::SourceContactMismatch),
        };
        let raw = self.instance.transform().transform_vector(local_normal);
        let length = norm(raw);
        let normal = Vec3::new(raw.x / length, raw.y / length, raw.z / length);
        if !normal.dot(contact.witness.normal).is_finite()
            || normal.dot(contact.witness.normal) < 1.0 - 1e-8
        {
            return Err(AreaEmitterError::SourceContactMismatch);
        }
        if !pdf_area_per_m2.is_finite() || pdf_area_per_m2 <= 0.0 {
            return Err(AreaEmitterError::NumericalRange);
        }
        Ok(Some(AreaEmitterSample {
            emitter: &self.emitter,
            instance: self.instance,
            support: &self.support,
            triangle,
            position_m,
            normal,
            pdf_area_per_m2,
        }))
    }

}

fn power_weights(
    light_log_pdf: f64,
    continuation_log_pdf: f64,
) -> Result<AreaImageMisWeights, AreaEmitterError> {
    if !light_log_pdf.is_finite() || continuation_log_pdf.is_nan()
        || continuation_log_pdf == f64::INFINITY
    {
        return Err(AreaEmitterError::NumericalRange);
    }
    let (light, continuation) = if continuation_log_pdf == f64::NEG_INFINITY {
        (1.0, 0.0)
    } else if light_log_pdf >= continuation_log_pdf {
        let ratio_squared = det::exp(2.0 * (continuation_log_pdf - light_log_pdf));
        (1.0 / (1.0 + ratio_squared), ratio_squared / (1.0 + ratio_squared))
    } else {
        let ratio_squared = det::exp(2.0 * (light_log_pdf - continuation_log_pdf));
        (ratio_squared / (1.0 + ratio_squared), 1.0 / (1.0 + ratio_squared))
    };
    // A positive but unrepresentable weight is not exact darkness: later
    // radiance factors could amplify it. Refuse instead of silently dropping it.
    if !light.is_finite() || light <= 0.0 || !continuation.is_finite()
        || (continuation_log_pdf.is_finite() && continuation <= 0.0)
    {
        return Err(AreaEmitterError::NumericalRange);
    }
    Ok(AreaImageMisWeights { light, continuation, light_log_pdf, continuation_log_pdf })
}


