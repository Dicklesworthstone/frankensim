//! Material-independent volume moments; no named shape implies a material.
use fs_ga::{Mat3, Se3, Vec3};
use fs_mbd::articulated::SpatialInertia;
use crate::{MaterialInertiaError, positive};

/// Integral geometry of a uniform region, expressed in an explicit body frame.
/// The tensor is integral (|r-c|² 1 - (r-c)(r-c)^T) dV [m5], NOT kg m2.
/// Multiplication by source density yields the existing centroidal mass inertia.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VolumeMoments {
    volume_m3: f64,
    centroid_m: Vec3,
    inertia_m5: Mat3,
}
impl VolumeMoments {
    /// Admit geometry supplied by an analytic shape or a mesh integration owner.
    /// Numeric inertia consistency is checked by the existing dynamics owner;
    /// this does not prove mesh coverage, absence of overlap or geometric accuracy.
    pub fn new(volume_m3: f64, centroid_m: Vec3, inertia_m5: Mat3)
        -> Result<Self, MaterialInertiaError> {
        let checked = SpatialInertia::new(volume_m3, centroid_m, inertia_m5)?;
        Ok(Self { volume_m3, centroid_m, inertia_m5: checked.inertia_com() })
    }
    /// Volume [m3].
    pub const fn volume_m3(self) -> f64 { self.volume_m3 }
    /// Geometric centroid in the body frame [m].
    pub const fn centroid_m(self) -> Vec3 { self.centroid_m }
    /// Centroidal volume-inertia tensor in the body frame [m5].
    pub const fn inertia_m5(self) -> Mat3 { self.inertia_m5 }

    /// Centered solid sphere with explicit radius [m]. The volume and all three
    /// centroidal moments come from this SAME geometry, before material density.
    pub fn sphere(radius_m: f64) -> Result<Self, MaterialInertiaError> {
        if !positive(radius_m) {
            return Err(MaterialInertiaError::Invalid("sphere radius must be finite and positive"));
        }
        let volume = (4.0 * core::f64::consts::PI / 3.0) * radius_m * radius_m * radius_m;
        let moment = (2.0 / 5.0) * volume * radius_m * radius_m;
        Self::new(volume, Vec3::new(0.0, 0.0, 0.0), Mat3 { m: [
            moment, 0.0, 0.0, 0.0, moment, 0.0, 0.0, 0.0, moment,
        ] })
    }

    /// Centered solid ellipsoid with positive semiaxes [m]. Principal volume
    /// moments are V*(b²+c²)/5 and cyclic permutations. Use `transformed` to
    /// express the SAME geometry and its complete tensor in another body frame.
    pub fn ellipsoid(semiaxes_m: Vec3) -> Result<Self, MaterialInertiaError> {
        let [a, b, c] = [semiaxes_m.x, semiaxes_m.y, semiaxes_m.z];
        if ![a, b, c].into_iter().all(positive) {
            return Err(MaterialInertiaError::Invalid(
                "ellipsoid semiaxes must be finite and positive",
            ));
        }
        let volume = (4.0 * core::f64::consts::PI / 3.0) * a * b * c;
        Self::new(
            volume,
            Vec3::new(0.0, 0.0, 0.0),
            Mat3 {
                m: [
                    (volume / 5.0) * (b * b + c * c),
                    0.0,
                    0.0,
                    0.0,
                    (volume / 5.0) * (a * a + c * c),
                    0.0,
                    0.0,
                    0.0,
                    (volume / 5.0) * (a * a + b * b),
                ],
            },
        )
    }

    /// Centered axis-aligned solid cuboid with explicit full side lengths [m].
    pub fn cuboid(sides_m: Vec3) -> Result<Self, MaterialInertiaError> {
        let [x, y, z] = [sides_m.x, sides_m.y, sides_m.z];
        if ![x, y, z].into_iter().all(positive) {
            return Err(MaterialInertiaError::Invalid("cuboid side lengths must be finite and positive"));
        }
        let v = x*y*z;
        Self::new(v, Vec3::new(0.0, 0.0, 0.0), Mat3 { m: [
            (v/12.0)*(y*y+z*z), 0.0, 0.0,
            0.0, (v/12.0)*(x*x+z*z), 0.0,
            0.0, 0.0, (v/12.0)*(x*x+y*y),
        ] })
    }

    /// Exact polynomial volume moments of a nondegenerate straight tetrahedron.
    /// Vertex order does not choose mass sign. Vertices are in metres in the
    /// destination body frame. A tetrahedral partition may assign different
    /// material cards to different cells without changing a dynamics equation.
    /// Integration is analytic; floating-point geometry is not certified exact.
    pub fn tetrahedron(vertices_m: [Vec3; 4]) -> Result<Self, MaterialInertiaError> {
        if vertices_m.iter().any(|v| ![v.x, v.y, v.z].into_iter().all(f64::is_finite)) {
            return Err(MaterialInertiaError::Invalid("tetrahedron vertices must be finite"));
        }
        let origin = vertices_m[0];
        let offsets = vertices_m.map(|v| v - origin);
        let volume = offsets[1].dot(offsets[2].cross(offsets[3])).abs() / 6.0;
        let offset_center = (offsets[1] + offsets[2] + offsets[3]).scale(0.25);
        let center = origin + offset_center;
        let mut inertia = Mat3 { m: [0.0; 9] };
        // E[(r-c)(r-c)^T] = sum_i (v_i-c)(v_i-c)^T / 20.
        for offset in offsets {
            let d = offset - offset_center;
            let [x, y, z] = [d.x, d.y, d.z];
            let terms = [y*y+z*z, -x*y, -x*z, -x*y, x*x+z*z, -y*z, -x*z, -y*z, x*x+y*y];
            for (entry, term) in inertia.m.iter_mut().zip(terms) { *entry += (volume/20.0)*term; }
        }
        Self::new(volume, center, inertia)
    }

    /// Place the geometry using the same validated SE(3) owner as dynamics.
    /// Rotation acts on the complete tensor; translation acts on the centroid.
    /// No parallel-axis shift is added to a tensor that is still centroidal.
    pub fn transformed(self, body_from_region: Se3) -> Result<Self, MaterialInertiaError> {
        let r = body_from_region.rotation().matrix();
        Self::new(self.volume_m3,
            r.apply(self.centroid_m) + body_from_region.translation(),
            r.compose(self.inertia_m5).compose(r.transpose()))
    }

    pub(crate) fn with_density(self, density: f64) -> Result<SpatialInertia, MaterialInertiaError> {
        if !positive(density) { return Err(MaterialInertiaError::Invalid("density must be finite and positive")); }
        let mut tensor = self.inertia_m5;
        for entry in &mut tensor.m { *entry *= density; }
        Ok(SpatialInertia::new(self.volume_m3*density, self.centroid_m, tensor)?)
    }
}
