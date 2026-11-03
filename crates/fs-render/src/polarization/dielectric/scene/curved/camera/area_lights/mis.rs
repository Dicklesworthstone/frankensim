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

use super::{AreaEmitterError, AreaEmitterSample, CurvedAreaEmitter, cross, delta, norm};
use crate::gradient_index::bundle::ProjectionOptions;
use crate::instances::SharedGeometry;
use crate::polarization::dielectric::scene::curved::CurvedContact;
use crate::polarization::dielectric::scene::curved::bundle::connection::{
    BundleConnection, BundleConnectionTarget,
};
use fs_geom::{Point3, Vec3};
use fs_math::det;

/// Complementary power-heuristic weights for one regular optical image.
/// Densities use the same receiver solid-angle measure. An image absent from
/// a finite light catalog must instead retain emitter-hit weight one.
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
        if !core::ptr::eq(sample.emitter, self.emitter)
            || !core::ptr::eq(sample.instance, self.instance)
        {
            return Err(AreaEmitterError::SampleOwnerMismatch);
        }
        let path = connection.path();
        let contact = path.source_contact().ok_or(AreaEmitterError::SourceContactMismatch)?;
        let receiver = path.last_contact().ok_or(AreaEmitterError::SourceContactMismatch)?;
        let bundle = path.bundle();
        let geometric = bundle.propagation();
        if !target.is_diffuse_receiver()
            || !path.is_diffuse_receiver()
            || receiver.object_id != target.object_id()
            || contact.object_id != self.emitter.object_id()
            || contact.geometry_identity != self.instance.geometry_identity()
            || contact.frame_identity != self.instance.frame_identity()
            || geometric.initial_ray().origin != sample.position_m
            || contact.witness.normal.dot(sample.normal) < 1.0 - 1e-8
        {
            return Err(AreaEmitterError::SourceContactMismatch);
        }
        let cosine_source = contact.witness.normal.dot(geometric.initial_ray().dir).abs().min(1.0);
        let cosine_receiver = -target.normal().dot(geometric.propagation().ray.dir);
        if !cosine_source.is_finite() || cosine_source <= 0.0
            || !cosine_receiver.is_finite() || cosine_receiver <= 0.0
        {
            return Err(AreaEmitterError::SourceContactMismatch);
        }
        let cosine_receiver = cosine_receiver.min(1.0);
        let area = geometric.project_at_endpoint(target.normal(), projection)
            .map_err(|error| AreaEmitterError::Refractive(error.into()))?
            .area_per_solid_angle_m2();
        let source_index = geometric.propagation().initial_sample.phase_index();
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
        let continuation_log_pdf = det::ln(cosine_receiver)
            - det::ln(core::f64::consts::PI) + reverse_log_probability;
        power_weights(light_log_pdf, continuation_log_pdf)
    }

    /// Recover the original represented area density at an actual native
    /// emitter hit. The witness locates its original face; the returned launch
    /// position stays on the integrated camera ray and is never snapped.
    /// Triangle edges and ambiguous membership refuse, matching sample_area.
    pub(super) fn sample_at_contact(
        &self,
        contact: &CurvedContact,
        position_m: Point3,
    ) -> Result<AreaEmitterSample<'_>, AreaEmitterError> {
        if contact.object_id != self.emitter.object_id()
            || contact.geometry_identity != self.instance.geometry_identity()
            || contact.frame_identity != self.instance.frame_identity()
            || ![position_m.x, position_m.y, position_m.z].into_iter().all(f64::is_finite)
        {
            return Err(AreaEmitterError::SourceContactMismatch);
        }
        let SharedGeometry::Mesh(mesh) = self.instance.geometry() else {
            return Err(AreaEmitterError::UnsupportedGeometry);
        };
        let local = self.instance.transform().inverse()
            .map_err(|_| AreaEmitterError::NumericalRange)?
            .transform_point(contact.witness.surface_point);
        let p = [local.x, local.y, local.z];
        let mut found = None;
        for (triangle, indices) in mesh.triangles.iter().enumerate() {
            let [a, b, c] = indices.map(|index| mesh.vertices[index as usize]);
            if !interior_barycentric(p, a, b, c)? {
                continue;
            }
            if found.is_some() {
                return Err(AreaEmitterError::SourceContactMismatch);
            }
            let raw = self.instance.transform().transform_vector(cross(delta(b, a), delta(c, a)));
            let length = norm(raw);
            let normal = Vec3::new(raw.x / length, raw.y / length, raw.z / length);
            if !normal.dot(contact.witness.normal).is_finite()
                || normal.dot(contact.witness.normal) < 1.0 - 1e-8
            {
                return Err(AreaEmitterError::SourceContactMismatch);
            }
            let entry = &self.triangles[triangle];
            let lower = if triangle == 0 { 0.0 } else { self.triangles[triangle - 1].cumulative_m2 };
            let pdf_area_per_m2 = ((entry.cumulative_m2 - lower) / self.total_area_m2) / entry.area_m2;
            if !pdf_area_per_m2.is_finite() || pdf_area_per_m2 <= 0.0 {
                return Err(AreaEmitterError::NumericalRange);
            }
            found = Some(AreaEmitterSample {
                emitter: self.emitter,
                instance: self.instance,
                triangle,
                position_m,
                normal,
                pdf_area_per_m2,
            });
        }
        found.ok_or(AreaEmitterError::SourceContactMismatch)
    }
}

// Work in the dominant coordinate plane and normalize local differences
// before products. The original native witness supplies the surface point;
// this lookup identifies its original triangle, not a prepared patch ID.
fn interior_barycentric(
    point: [f64; 3],
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
) -> Result<bool, AreaEmitterError> {
    let e = delta(b, a);
    let f = delta(c, a);
    let v = delta(point, a);
    let n = cross(e, f);
    let (u, w) = if n.x.abs() >= n.y.abs() && n.x.abs() >= n.z.abs() {
        (1, 2)
    } else if n.y.abs() >= n.z.abs() {
        (0, 2)
    } else {
        (0, 1)
    };
    let e = [e.x, e.y, e.z];
    let f = [f.x, f.y, f.z];
    let v = [v.x, v.y, v.z];
    let scale = e[u].abs().max(e[w].abs()).max(f[u].abs()).max(f[w].abs());
    if !scale.is_finite() || scale <= 0.0 {
        return Err(AreaEmitterError::NumericalRange);
    }
    let ex = e[u] / scale;
    let ey = e[w] / scale;
    let fx = f[u] / scale;
    let fy = f[w] / scale;
    let vx = v[u] / scale;
    let vy = v[w] / scale;
    let determinant = ex * fy - ey * fx;
    let first = (vx * fy - vy * fx) / determinant;
    let second = (ex * vy - ey * vx) / determinant;
    let third = 1.0 - first - second;
    if ![first, second, third].into_iter().all(f64::is_finite) {
        return Err(AreaEmitterError::NumericalRange);
    }
    Ok(first > 0.0 && second > 0.0 && third > 0.0)
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

