//! Whole-source-patch extinction bounds for the existing area integrator.
//!
//! For a nonnegative geometric kernel K and homogeneous extinction sigma,
//! T(r_max) integral(K) <= integral(K T(r)) <= T(r_min) integral(K).
//! No centroid ray or quadrature discrepancy is used as an error bound. Reuse
//! HomogeneousMedium's retained coefficients and fs-ivl's proved exponential.

use super::*;
use fs_ivl::decay::exp_neg;

/// Prepared numerical view, not a new physical-medium owner. Keep both original
/// coefficients for provenance even when their sums happen to be identical.
pub(super) struct PreparedAttenuation {
    absorption_per_m: f64,
    scattering_per_m: f64,
    extinction: Interval,
}

impl PreparedAttenuation {
    pub(super) fn vacuum() -> Self {
        Self { absorption_per_m: 0.0, scattering_per_m: 0.0, extinction: zero() }
    }

    pub(super) fn new(medium: &HomogeneousMedium) -> Self {
        let absorption_per_m = medium.absorption_per_m();
        let scattering_per_m = medium.scattering_per_m();
        // Same real function as the original direct point-fixture enclosure:
        // sum the exact retained coefficients outward, not their rounded sum.
        let extinction = nonnegative(point(absorption_per_m) + point(scattering_per_m));
        Self { absorption_per_m, scattering_per_m, extinction }
    }

    pub(super) fn absorption_per_m(&self) -> f64 { self.absorption_per_m }
    pub(super) fn scattering_per_m(&self) -> f64 { self.scattering_per_m }
    pub(super) fn is_vacuum(&self) -> bool {
        self.absorption_per_m == 0.0 && self.scattering_per_m == 0.0
    }
    pub(super) fn hash_into(&self, hash: &mut DomainHasher) {
        // The phase law changes only omitted in-scattering, not this direct
        // numerical model. Do not claim this identity describes all transport.
        for value in [self.absorption_per_m, self.scattering_per_m] {
            hash_float(hash, if value == 0.0 { 0.0 } else { value });
        }
    }

    pub(super) fn patch_transmittance(
        &self, vectors: [Vector; 3], perpendicular: Vector, box_distance: Interval,
    ) -> Option<Interval> {
        if self.is_vacuum() { return Some(point(1.0)); }
        let distance = patch_distance(vectors, perpendicular, box_distance)?;
        // Finite distance/coefficient products may overflow their upper endpoint;
        // exp_neg accepts positive infinity and retains positive underflow bounds.
        exp_neg(nonnegative(self.extinction * distance))
    }
}

/// The norm is convex: its maximum over a triangle is at a corner. A lower
/// bound is the larger of the enclosing box's distance and the supporting
/// plane's distance. The latter avoids a false zero for broad tilted patches.
/// All candidates are outward enclosures, including rigid-placement roundoff.
fn patch_distance(vectors: [Vector; 3], perpendicular: Vector, box_distance: Interval) -> Option<Interval> {
    if !finite(box_distance) { return None; }
    let normal_length = norm(perpendicular)?;
    if normal_length.lo() <= 0.0 { return None; }
    let plane = dot(perpendicular, vectors[0]).abs() / normal_length;
    if !finite(plane) { return None; }
    let lower = box_distance.lo().max(plane.lo()).max(0.0);
    let mut upper = 0.0_f64;
    for corner in vectors {
        upper = upper.max(norm(corner)?.hi());
    }
    if !upper.is_finite() || lower > upper { return None; }
    Some(Interval::new(lower, upper))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::photometry::medium::PhaseFunction;

    #[test]
    fn supporting_plane_bound_covers_a_patch_whose_box_contains_the_receiver() {
        // Triangle x+y+z=1: minimum distance is 1/sqrt(3), maximum is 1.
        let vectors = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]].map(|v| v.map(point));
        let distance = patch_distance(vectors, [point(1.0);3], Interval::new(0.0,2.0)).unwrap();
        assert!(distance.lo() > 0.57 && distance.lo() <= 1.0/3.0_f64.sqrt());
        assert!(distance.hi() >= 1.0 && distance.hi() < 1.0+1e-12);
        let medium = HomogeneousMedium::try_new(0.25,0.25,PhaseFunction::Isotropic).unwrap();
        let transfer = PreparedAttenuation::new(&medium)
            .patch_transmittance(vectors,[point(1.0);3],Interval::new(0.0,2.0)).unwrap();
        assert!(transfer.lo() <= (-0.5_f64).exp());
        assert!(transfer.hi() >= (-0.5/3.0_f64.sqrt()).exp());
        assert!(transfer.hi() < 0.76);
    }

    #[test]
    fn corner_maximum_avoids_a_nonexistent_far_corner_of_the_box() {
        let vectors = [[1.0,0.0,2.0],[0.0,1.0,2.0],[0.0,0.0,2.0]].map(|v| v.map(point));
        let distance = patch_distance(vectors,[point(0.0),point(0.0),point(1.0)],Interval::new(2.0,3.0)).unwrap();
        assert!(distance.lo() <= 2.0 && distance.lo() > 1.99);
        assert!(distance.hi() >= 5.0_f64.sqrt() && distance.hi() < 2.237);
    }
}
