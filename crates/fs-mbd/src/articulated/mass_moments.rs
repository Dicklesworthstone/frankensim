//! Analytic weighted geometry; source/material resolution belongs upstream.

use super::{ArticulatedError, SpatialInertia, validate_slice, validate_vec3};
use fs_ga::{Mat3, Vec3};

impl SpatialInertia {
    /// Exact polynomial mass moments of a straight tetrahedron with positive
    /// linearly interpolated vertex weights and an explicitly supplied mass.
    ///
    /// Weights set the relative density distribution, not the total mass; their
    /// common scale has no physical effect. This permits a reference inventory
    /// to retain its exact mass during an affine deformation. Coordinates are
    /// metres in one body frame. The returned centroidal tensor is kg m².
    ///
    /// Integration uses the analytic barycentric moments, with fixed work and
    /// no quadrature allocation. Centred positive dyads avoid subtraction of
    /// origin moments. This is ordinary floating-point integration, not a
    /// geometric enclosure, material measurement or overlap certificate.
    ///
    /// # Errors
    /// Refuses nonpositive/nonfinite mass or weights, nonfinite/collapsed
    /// geometry, and any derived inertia refused by the original owner.
    pub fn from_weighted_tetrahedron(
        mass: f64,
        vertices_m: [Vec3; 4],
        weights: [f64; 4],
    ) -> Result<Self, ArticulatedError> {
        if !mass.is_finite() || mass <= 0.0 {
            return Err(ArticulatedError::InvalidMass);
        }
        validate_slice(&weights, "tetrahedron.weights")?;
        if weights.iter().any(|&w| w <= 0.0) {
            return Err(ArticulatedError::InvalidMass);
        }
        for vertex in vertices_m {
            validate_vec3(vertex, "tetrahedron.vertex")?;
        }
        let origin = vertices_m[0];
        let offsets = vertices_m.map(|v| v - origin);
        let determinant = offsets[1].dot(offsets[2].cross(offsets[3]));
        validate_slice(&[determinant], "tetrahedron.determinant")?;
        if determinant == 0.0 {
            return Err(ArticulatedError::NonPhysicalInertia { measure: 0.0 });
        }
        let scale = weights.into_iter().fold(0.0_f64, f64::max);
        let normalized = weights.map(|w| w / scale);
        let total: f64 = normalized.iter().sum();
        let fractions = normalized.map(|w| w / total);
        let geometric_center = (offsets[1] + offsets[2] + offsets[3]).scale(0.25);
        let centered = offsets.map(|v| v - geometric_center);
        let shift = centered
            .into_iter()
            .zip(fractions)
            .fold(Vec3::new(0.0, 0.0, 0.0), |sum, (v, w)| sum + v.scale(w));
        let center = origin + geometric_center + shift.scale(0.2);

        // For q_i = vertex_i - geometric centroid and a = sum(w_i*q_i),
        // C/M = sum(q_i*q_i')/30 + sum(w_i*(q_i-a)*(q_i-a)')/15
        //       + 2*a*a'/75. Every dyad has nonnegative weight. The density
        // centroid is the geometric centroid + a/5. These are the exact
        // degree-three barycentric integrals, including all cross moments.
        let mut tensor = Mat3 { m: [0.0; 9] };
        for (offset, fraction) in centered.into_iter().zip(fractions) {
            add_dyad_inertia(&mut tensor, mass / 30.0, offset);
            add_dyad_inertia(&mut tensor, (mass / 15.0) * fraction, offset - shift);
        }
        add_dyad_inertia(&mut tensor, mass * (2.0 / 75.0), shift);
        Self::new(mass, center, tensor)
    }
}

fn add_dyad_inertia(tensor: &mut Mat3, weight: f64, offset: Vec3) {
    let [x, y, z] = [offset.x, offset.y, offset.z];
    let terms = [
        y * y + z * z,
        -x * y,
        -x * z,
        -x * y,
        x * x + z * z,
        -y * z,
        -x * z,
        -y * z,
        x * x + y * y,
    ];
    for (entry, term) in tensor.m.iter_mut().zip(terms) {
        *entry += weight * term;
    }
}
