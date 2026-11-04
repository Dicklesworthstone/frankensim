//! A private departure is rebound to two original adjacent tetrahedra.

use super::*;
use crate::gradient_index::event::{
    BoundaryOracle, BoundaryQuery, BoundarySide, BoundaryWitness, DepartureNeighborhood,
    DepartureRequest,
};
use fs_ivl::Interval;

pub(super) struct Transition<'a> {
    previous: P1ConductionElement<'a>,
    next: P1ConductionElement<'a>,
    next_face: usize,
    origin: Point3,
    witness: BoundaryWitness,
}

impl<'a> Transition<'a> {
    pub(super) fn new(
        previous: P1ConductionElement<'a>,
        next: P1ConductionElement<'a>,
        next_face: usize,
        origin: Point3,
        witness: BoundaryWitness,
    ) -> Self {
        Self {
            previous,
            next,
            next_face,
            origin,
            witness,
        }
    }

    fn prepare(
        &self,
        request: DepartureRequest,
    ) -> Result<DepartureNeighborhood, GradientIndexError> {
        if request.origin != self.origin
            || request.contact.witness != self.witness
            || request.contact.side != BoundarySide::Positive
        {
            return Err(GradientIndexError::UnsupportedBoundaryDeparture);
        }
        let witness = self.witness;
        let invalid = |_| GradientIndexError::InvalidBoundaryQuery;
        let (_, original) = triangle_proximity(
            self.origin,
            self.previous.faces_m()[witness.feature_id],
            witness.feature_id as u32,
        )
        .map_err(invalid)?;
        if original.surface_point != witness.surface_point
            || original.normal != Some(witness.normal)
            || original.side != ProximitySide::Negative
            || witness.side != BoundarySide::Negative
            || original.distance_upper_m > witness.distance_upper_m
        {
            return Err(GradientIndexError::UnsupportedBoundaryDeparture);
        }
        // The arrival owner already excludes the original face's edges and
        // the previous tetrahedron's other faces. Also exclude every other
        // face of the destination tetrahedron, including representation error.
        let mut clearance = f64::INFINITY;
        for (face, vertices) in self.next.faces_m().into_iter().enumerate() {
            if face != self.next_face {
                clearance = clearance.min(
                    triangle_proximity(witness.surface_point, vertices, face as u32)
                        .map_err(invalid)?
                        .0,
                );
            }
        }
        let radius_m = ((Interval::point(clearance)
            - Interval::point(original.surface_point_error_m))
            * Interval::point(0.5))
        .lo()
        .min(
            witness
                .departure_radius_m
                .ok_or(GradientIndexError::UnsupportedBoundaryDeparture)?,
        );
        if !radius_m.is_finite() || radius_m <= 4.0 * request.surface_tolerance_m {
            return Err(GradientIndexError::BoundaryResolutionLimit);
        }
        // This proof runs inside the native charged prepare_departure call,
        // so it shares the whole-path geometry budget and cancellation point.
        Ok(DepartureNeighborhood::Planar { radius_m })
    }

    pub(super) fn departure(&self) -> GradientDeparture {
        GradientDeparture {
            witness: self.witness,
            side: BoundarySide::Positive,
        }
    }
}

pub(super) struct MeshBoundary<'a> {
    current: ConductionElementBoundary,
    transition: Option<Transition<'a>>,
}

impl<'a> MeshBoundary<'a> {
    pub(super) fn new(current: ConductionElementBoundary) -> Self {
        Self {
            current,
            transition: None,
        }
    }
    pub(super) fn with_transition(
        current: ConductionElementBoundary,
        transition: Transition<'a>,
    ) -> Self {
        Self {
            current,
            transition: Some(transition),
        }
    }
}

impl BoundaryOracle for MeshBoundary<'_> {
    fn query(&mut self, position: Point3) -> Result<BoundaryQuery, GradientIndexError> {
        self.current.query(position)
    }

    fn prepare_departure(
        &mut self,
        request: DepartureRequest,
    ) -> Result<DepartureNeighborhood, GradientIndexError> {
        let transition = self
            .transition
            .as_ref()
            .ok_or(GradientIndexError::UnsupportedBoundaryDeparture)?;
        // The native owner independently checks initial proximity, forward
        // direction, all stages and progress inside this proven finite ball.
        transition.prepare(request)
    }
}
