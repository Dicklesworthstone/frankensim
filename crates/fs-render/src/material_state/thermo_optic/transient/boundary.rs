//! Four original mesh triangles through the native interval proximity owner.

use super::P1ConductionElement;
use crate::charts::proximity::{ProximitySide, segment_clearance, triangle_proximity};
use crate::gradient_index::GradientIndexError;
use crate::gradient_index::event::{
    BoundaryOracle, BoundaryQuery, BoundarySide, BoundaryWitness, DepartureNeighborhood,
    DepartureRequest,
};
use fs_geom::Point3;
use fs_ivl::Interval;

/// Native first-face geometry for one original tetrahedron. It reuses the
/// existing interval lower clearance and independent convex upper witness.
/// At most four triangles are queried normally; a near contact additionally
/// checks three edges and the three other faces. No heap allocation occurs.
///
/// This oracle covers these four faces only. Other physical scene surfaces
/// need the existing scene oracle. Face events identify the local opposite
/// vertex (0..4); the caller supplies the enclosing boundary identity.
/// This owner supports arrival only. Its proven planar radius guards the
/// unresolved incident tail at shallow incidence; it does not authorize a
/// reflected/transmitted continuation or accept caller-authored witnesses.
pub struct ConductionElementBoundary {
    faces_m: [[[f64; 3]; 3]; 4],
    boundary_id: usize,
    contact_tolerance_m: f64,
}

impl P1ConductionElement<'_> {
    /// Geometry for `gradient_index::event::propagate_to_boundary` or the
    /// native bundle event owner. Use the same metric contact tolerance in
    /// that request. An interior launch stops before leaving this thermal
    /// element; the oracle never skips an interpolation face.
    pub fn boundary_oracle(
        &self,
        boundary_id: usize,
        contact_tolerance_m: f64,
    ) -> Result<ConductionElementBoundary, GradientIndexError> {
        if !contact_tolerance_m.is_finite() || contact_tolerance_m <= 0.0 {
            return Err(GradientIndexError::InvalidOptions);
        }
        Ok(ConductionElementBoundary {
            faces_m: self.faces_m(),
            boundary_id,
            contact_tolerance_m,
        })
    }
}

impl BoundaryOracle for ConductionElementBoundary {
    fn query(&mut self, position: Point3) -> Result<BoundaryQuery, GradientIndexError> {
        let invalid = |_| GradientIndexError::InvalidBoundaryQuery;
        let first = triangle_proximity(position, self.faces_m[0], 0).map_err(invalid)?;
        let mut candidates = [first; 4];
        let mut selected = 0;
        let mut clearance = first.0;
        for (face, candidate) in candidates.iter_mut().enumerate().skip(1) {
            *candidate =
                triangle_proximity(position, self.faces_m[face], face as u32).map_err(invalid)?;
            clearance = clearance.min(candidate.0);
        }
        for face in 1..4 {
            if candidates[face].1.distance_upper_m < candidates[selected].1.distance_upper_m {
                selected = face;
            }
        }
        let witness = candidates[selected].1;
        if witness.distance_upper_m > self.contact_tolerance_m {
            return Ok(BoundaryQuery {
                clearance_m: clearance,
                witness: None,
            });
        }
        // A rounded representative differs from a true point on the original
        // face by its existing interval error. The following bound includes
        // that collar without promoting a tiny plane residual into a hit.
        let collar = (Interval::point(witness.distance_upper_m)
            + Interval::point(2.0) * Interval::point(witness.surface_point_error_m))
        .hi();
        if collar > self.contact_tolerance_m {
            return Err(GradientIndexError::BoundaryResolutionLimit);
        }
        let mut radius = f64::INFINITY;
        let vertices = self.faces_m[selected];
        for edge in 0..3 {
            radius = radius.min(
                segment_clearance(
                    witness.surface_point,
                    vertices[edge],
                    vertices[(edge + 1) % 3],
                )
                .map_err(invalid)?,
            );
        }
        for face in 0..4 {
            if face != selected {
                radius = radius.min(
                    triangle_proximity(witness.surface_point, self.faces_m[face], face as u32)
                        .map_err(invalid)?
                        .0,
                );
            }
        }
        radius = ((Interval::point(radius) - Interval::point(witness.surface_point_error_m))
            * Interval::point(0.5))
        .lo();
        // Corners, competing faces and unresolved translated representatives
        // do not acquire a unique normal by choosing the first triangle.
        if !radius.is_finite() || radius <= 4.0 * self.contact_tolerance_m {
            return Err(GradientIndexError::BoundaryResolutionLimit);
        }
        let normal = witness
            .normal
            .ok_or(GradientIndexError::InvalidBoundaryQuery)?;
        Ok(BoundaryQuery {
            clearance_m: clearance,
            witness: Some(BoundaryWitness {
                boundary_id: self.boundary_id,
                feature_id: selected,
                surface_point: witness.surface_point,
                normal,
                distance_upper_m: collar,
                side: match witness.side {
                    ProximitySide::Positive => BoundarySide::Positive,
                    ProximitySide::Negative => BoundarySide::Negative,
                    ProximitySide::Unresolved => BoundarySide::OnSurface,
                },
                departure_radius_m: Some(radius),
            }),
        })
    }

    fn prepare_departure(
        &mut self,
        _request: DepartureRequest,
    ) -> Result<DepartureNeighborhood, GradientIndexError> {
        // Do not inherit the default planar departure: an external witness
        // has not been rebound to this original mesh and local face. Crossing
        // also requires the next original thermal field; no exterior collar
        // or automatic temperature/material selection is supplied here.
        Err(GradientIndexError::UnsupportedBoundaryDeparture)
    }
}
