//! Area proposals for the actual triangles of physical thermal surfaces.
//!
//! The Planck law stays attached to each sample. Visible CIE luminance only
//! controls source selection; it never replaces the transported spectrum.

use fs_blake3::{ContentHash, DomainHasher};
use fs_geom::{Point3, Vec3};

use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN, QUAD_BINS, cie_y};
use crate::thermal_surface::GrayDiffuseThermalSurface;

use super::{LightingError, TAU, canonical_f64_bits, cross, triangle_solid_angle};

/// Source-bound triangle from accepted static or rigid-instance geometry.
#[derive(Clone, Copy, Debug)]
pub(crate) struct ThermalTriangleLight {
    pub(crate) vertices: [Point3; 3],
    pub(crate) primitive_index: usize,
    pub(crate) surface: GrayDiffuseThermalSurface,
}

/// One physical thermal-source NEE sample, with the complete mixture PDF.
#[derive(Clone, Copy, Debug)]
pub struct ThermalLightSample {
    /// Primitive carrying the same gray/diffuse surface in the traced scene.
    pub primitive_index: usize,
    /// Sampled world-space triangle point, in metres.
    pub point: Point3,
    /// Triangle unit normal. Physical thermal emission is two-sided.
    pub normal: Vec3,
    /// The original source-bound Planck law, without an RGB approximation.
    pub surface: GrayDiffuseThermalSurface,
    /// Full light-mixture density at the shading point, in inverse steradians.
    pub pdf_solid_angle: f64,
}

#[derive(Clone, Copy, Debug)]
pub(super) struct AdmittedThermalTriangle {
    pub(super) source: ThermalTriangleLight,
    pub(super) identity: ContentHash,
    pub(super) normal: Vec3,
    pub(super) luminance: f64,
    pub(super) area: f64,
    edge_u: Vec3,
    edge_v: Vec3,
    coordinate_scale: f64,
}

impl AdmittedThermalTriangle {
    pub(super) fn try_new(mut source: ThermalTriangleLight, luminance: f64) -> Option<Self> {
        if source
            .vertices
            .iter()
            .any(|point| !point.x.is_finite() || !point.y.is_finite() || !point.z.is_finite())
            || !luminance.is_finite()
            || luminance < 0.0
            || !(luminance * TAU).is_finite()
        {
            return None;
        }
        // A triangle's vertex numbering and winding do not change this
        // two-sided source. Canonicalize the sampling map as well as its hash.
        for point in &mut source.vertices {
            point.x = f64::from_bits(canonical_f64_bits(point.x));
            point.y = f64::from_bits(canonical_f64_bits(point.y));
            point.z = f64::from_bits(canonical_f64_bits(point.z));
        }
        source.vertices.sort_by(|left, right| {
            left.x
                .total_cmp(&right.x)
                .then_with(|| left.y.total_cmp(&right.y))
                .then_with(|| left.z.total_cmp(&right.z))
        });
        let edge_u = source.vertices[1].delta_from(source.vertices[0]);
        let edge_v = source.vertices[2].delta_from(source.vertices[0]);
        let unnormalized = cross(edge_u, edge_v);
        let (normal, twice_area) = normalized(unnormalized)?;
        let area = 0.5 * twice_area;
        if !(area > 0.0 && area.is_finite()) {
            return None;
        }
        let coordinate_scale = source.vertices.iter().fold(0.0_f64, |scale, point| {
            scale
                .max(point.x.abs())
                .max(point.y.abs())
                .max(point.z.abs())
        });
        let mut hasher = DomainHasher::new("org.frankensim.render.thermal-triangle-light.v1");
        for point in source.vertices {
            for value in [point.x, point.y, point.z] {
                hasher.update(&canonical_f64_bits(value).to_le_bytes());
            }
        }
        hasher.update(source.surface.identity().as_bytes());
        Some(Self {
            source,
            identity: hasher.finalize(),
            normal,
            luminance,
            area,
            edge_u,
            edge_v,
            coordinate_scale,
        })
    }

    pub(super) fn sample_point(self, u1: f64, u2: f64) -> Point3 {
        let radius = u1.sqrt();
        self.source.vertices[0]
            .offset(self.edge_u.scale(radius * (1.0 - u2)))
            .offset(self.edge_v.scale(radius * u2))
    }

    pub(super) fn weight(self, origin: Point3) -> f64 {
        if self.luminance == 0.0 {
            return 0.0;
        }
        let a = self.source.vertices[0].delta_from(origin);
        // A same-plane triangle has zero projected solid angle, including a
        // receiver that lies inside it. Do not let atan2(0, negative) invent a
        // 2*pi source at a self-hit. No arbitrary grazing cutoff removes
        // otherwise representable non-coplanar sampling support.
        let plane_distance = a.dot(self.normal);
        if !plane_distance.is_finite() || plane_distance == 0.0 {
            return 0.0;
        }
        let Some((a, _)) = normalized(a) else {
            return 0.0;
        };
        let Some((b, _)) = normalized(self.source.vertices[1].delta_from(origin)) else {
            return 0.0;
        };
        let Some((c, _)) = normalized(self.source.vertices[2].delta_from(origin)) else {
            return 0.0;
        };
        let angle = triangle_solid_angle(a, b, c);
        if !angle.is_finite() {
            return 0.0;
        }
        // A very small angular triangle can lose its triple product through
        // cancellation. The floor affects variance only and preserves the
        // proposal's support; it is not substituted physical radiance.
        (self.luminance * angle.max(f64::MIN_POSITIVE)).max(f64::MIN_POSITIVE)
    }

    pub(super) fn directional_pdf(self, origin: Point3, point: Point3) -> f64 {
        if !self.contains(point) {
            return 0.0;
        }
        let Some((direction, distance)) = normalized(point.delta_from(origin)) else {
            return 0.0;
        };
        let cosine = self.normal.dot(direction).abs();
        if cosine == 0.0 {
            return 0.0;
        }
        let pdf = (distance / self.area) * (distance / cosine);
        if pdf.is_finite() && pdf > 0.0 {
            pdf
        } else {
            0.0
        }
    }

    fn contains(self, point: Point3) -> bool {
        let relative = point.delta_from(self.source.vertices[0]);
        let plane_error = relative.dot(self.normal).abs();
        let scale = self
            .coordinate_scale
            .max(point.x.abs())
            .max(point.y.abs())
            .max(point.z.abs());
        if !plane_error.is_finite() || plane_error > 128.0 * f64::EPSILON * scale {
            return false;
        }
        // Project away the largest normal coordinate. Scale all coordinates
        // before the two-dimensional oriented-area ratios so metre units and
        // ordinary coordinate magnitude cannot overflow the barycentrics.
        let project = |vector: Vec3| {
            if self.normal.x.abs() >= self.normal.y.abs()
                && self.normal.x.abs() >= self.normal.z.abs()
            {
                [vector.y, vector.z]
            } else if self.normal.y.abs() >= self.normal.z.abs() {
                [vector.x, vector.z]
            } else {
                [vector.x, vector.y]
            }
        };
        let a = project(self.edge_u);
        let b = project(self.edge_v);
        let p = project(relative);
        let edge_scale = a[0].abs().max(a[1].abs()).max(b[0].abs()).max(b[1].abs());
        let a = a.map(|value| value / edge_scale);
        let b = b.map(|value| value / edge_scale);
        let p = p.map(|value| value / edge_scale);
        let determinant = a[0] * b[1] - a[1] * b[0];
        let u = (p[0] * b[1] - p[1] * b[0]) / determinant;
        let v = (a[0] * p[1] - a[1] * p[0]) / determinant;
        let tolerance = 128.0 * f64::EPSILON;
        u.is_finite()
            && v.is_finite()
            && u >= -tolerance
            && v >= -tolerance
            && u + v <= 1.0 + tolerance
    }
}

pub(super) fn visible_luminance(
    surface: GrayDiffuseThermalSurface,
    checkpoint: &mut impl FnMut() -> Result<(), LightingError>,
) -> Result<f64, LightingError> {
    let dl = (LAMBDA_MAX - LAMBDA_MIN) / QUAD_BINS as f64;
    let mut luminance = 0.0;
    let mut normalization = 0.0;
    for bin in 0..QUAD_BINS {
        if bin % 8 == 0 {
            checkpoint()?;
        }
        let lambda = LAMBDA_MIN + (bin as f64 + 0.5) * dl;
        let observer_weight = cie_y(lambda) * dl;
        let radiance = surface
            .radiance_w_m2_sr_nm(lambda)
            .map_err(|_| LightingError::InvalidThermalTriangle { light_index: 0 })?;
        luminance += radiance * observer_weight;
        normalization += observer_weight;
    }
    luminance /= normalization;
    if !luminance.is_finite() || luminance < 0.0 {
        return Err(LightingError::InvalidThermalTriangle { light_index: 0 });
    }
    // At the extreme Wien tail, an endpoint may still have representable
    // radiance after every fixed quadrature product has underflowed. Give it
    // nonzero proposal support without altering emitted spectral values.
    if luminance == 0.0
        && surface
            .radiance_w_m2_sr_nm(LAMBDA_MAX)
            .map_err(|_| LightingError::InvalidThermalTriangle { light_index: 0 })?
            > 0.0
    {
        luminance = f64::MIN_POSITIVE;
    }
    Ok(luminance)
}

fn normalized(vector: Vec3) -> Option<(Vec3, f64)> {
    let scale = vector.x.abs().max(vector.y.abs()).max(vector.z.abs());
    if !(scale > 0.0 && scale.is_finite()) {
        return None;
    }
    let scaled = vector.scale(1.0 / scale);
    let unit_length = scaled.dot(scaled).sqrt();
    let length = scale * unit_length;
    (length > 0.0 && length.is_finite()).then_some((scaled.scale(1.0 / unit_length), length))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lighting::{AdmittedLighting, LightSample};

    fn triangle(x: f64, primitive_index: usize) -> ThermalTriangleLight {
        ThermalTriangleLight {
            vertices: [
                Point3::new(x, 1.0, -0.5),
                Point3::new(x + 1.0, 1.0, -0.5),
                Point3::new(x, 1.0, 0.5),
            ],
            primitive_index,
            surface: GrayDiffuseThermalSurface::try_new(2400.0, 0.7, ContentHash([71; 32]))
                .unwrap(),
        }
    }

    fn admit(triangles: Vec<ThermalTriangleLight>) -> AdmittedLighting<'static> {
        AdmittedLighting::try_new_controlled(&[], None, triangles, false, || Ok(())).unwrap()
    }

    #[test]
    fn g0_thermal_forward_reverse_and_overlapping_support_agree() {
        let light = triangle(-0.5, 4);
        let single = admit(vec![light]);
        let repeated = admit(vec![light, light]);
        for origin in [Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 2.0, 0.0)] {
            for index in 0..128 {
                let u1 = (f64::from(index) + 0.5) / 128.0;
                let u2 = f64::from((index * 53) % 128) / 128.0;
                for lights in [&single, &repeated] {
                    let LightSample::ThermalSurface(sample) =
                        lights.sample(origin, u1, u2).unwrap()
                    else {
                        panic!("thermal source expected");
                    };
                    let reverse = lights.thermal_mixture_pdf(4, origin, sample.point);
                    assert_eq!(sample.pdf_solid_angle.to_bits(), reverse.to_bits());
                    let unnormalized = cross(
                        light.vertices[1].delta_from(light.vertices[0]),
                        light.vertices[2].delta_from(light.vertices[0]),
                    );
                    let area = 0.5 * unnormalized.norm();
                    let delta = sample.point.delta_from(origin);
                    let expected = delta.dot(delta)
                        / (area * sample.normal.dot(delta.scale(1.0 / delta.norm())).abs());
                    assert!((reverse / expected - 1.0).abs() < 2.0e-14);
                    assert_eq!(sample.surface, light.surface);
                }
            }
        }
    }

    #[test]
    fn g5_thermal_source_order_and_vertex_numbering_do_not_change_samples() {
        let a = triangle(-2.0, 1);
        let mut b = triangle(1.0, 2);
        b.surface = GrayDiffuseThermalSurface::try_new(2700.0, 0.9, ContentHash([72; 32])).unwrap();
        let forward = admit(vec![a, b]);
        let mut a_reversed = a;
        a_reversed.vertices.swap(0, 2);
        a_reversed.primitive_index = 22;
        let mut b_reversed = b;
        b_reversed.vertices.swap(0, 1);
        b_reversed.primitive_index = 11;
        let reversed = admit(vec![b_reversed, a_reversed]);
        let origin = Point3::new(0.0, 0.0, 0.0);
        for index in 0..128 {
            let u1 = (f64::from(index) + 0.5) / 128.0;
            let u2 = f64::from((index * 53) % 128) / 128.0;
            let LightSample::ThermalSurface(a) = forward.sample(origin, u1, u2).unwrap() else {
                panic!("thermal source expected");
            };
            let LightSample::ThermalSurface(b) = reversed.sample(origin, u1, u2).unwrap() else {
                panic!("thermal source expected");
            };
            assert_eq!(
                [a.point.x, a.point.y, a.point.z],
                [b.point.x, b.point.y, b.point.z]
            );
            assert_eq!(a.pdf_solid_angle.to_bits(), b.pdf_solid_angle.to_bits());
            assert_eq!(a.surface, b.surface);
        }
    }

    #[test]
    fn g0_thermal_proposal_refuses_degenerate_and_excludes_coplanar_or_off_surface_hits() {
        let source = triangle(-0.5, 4);
        let lights = admit(vec![source]);
        let in_plane = Point3::new(0.0, 1.0, 0.0);
        assert!(lights.sample(in_plane, 0.5, 0.5).is_none());
        assert_eq!(
            lights.thermal_mixture_pdf(4, in_plane, source.vertices[0]),
            0.0
        );
        let origin = Point3::new(0.0, 0.0, 0.0);
        assert_eq!(
            lights.thermal_mixture_pdf(4, origin, Point3::new(5.0, 1.0, 0.0)),
            0.0
        );
        assert_eq!(
            lights.thermal_mixture_pdf(4, origin, Point3::new(0.0, 1.1, 0.0)),
            0.0
        );
        let mut degenerate = source;
        degenerate.vertices[2] = degenerate.vertices[1];
        assert!(matches!(
            AdmittedLighting::try_new_controlled(&[], None, vec![degenerate], false, || Ok(())),
            Err(LightingError::InvalidThermalTriangle { light_index: 0 })
        ));
    }

    #[test]
    fn g4_thermal_admission_cancels_inside_spectral_quadrature() {
        let mut polls = 0;
        let result =
            AdmittedLighting::try_new_controlled(&[], None, vec![triangle(-0.5, 4)], false, || {
                polls += 1;
                if polls == 5 {
                    Err(LightingError::Cancelled)
                } else {
                    Ok(())
                }
            });
        assert!(matches!(result, Err(LightingError::Cancelled)));
        assert_eq!(polls, 5);
    }
}
