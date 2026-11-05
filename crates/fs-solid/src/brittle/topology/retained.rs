//! Admission of retained material-point geometry through additional cuts.

use super::{FractureTopologyError, TetFractureRefinement, TetFractureTopology};
use crate::cohesive::p1::{P1Interfaces, frame};
use fs_exec::Cx;

impl TetFractureRefinement<'_, '_> {
    /// Require identical retained physical points, frames, ordered surfaces and
    /// original law/source owners on the complete matching P1 seam maps.
    /// Construct both maps with the ordered constructors and `SplitFace` order.
    ///
    /// This is linear in the admitted seam count, uses no mesh-sized temporary
    /// allocation and checks reference geometry with the native P1 frame owner.
    /// Free/prescribed DOFs may change: boundary/load reassignment remains the
    /// caller's responsibility. This does not rebind an opaque solver state or
    /// initialize new cohesive histories, nor does it admit a material change.
    pub fn require_retained_interfaces(
        &self,
        previous: &P1Interfaces<'_>,
        next: &P1Interfaces<'_>,
        cx: &Cx<'_>,
    ) -> Result<(), FractureTopologyError> {
        self.require_retained_interfaces_controlled(previous, next, &mut || {
            cx.checkpoint().is_err()
        })
    }

    /// Predicate form with identical geometry, ownership and cancellation gates.
    pub fn require_retained_interfaces_controlled<C: FnMut() -> bool>(
        &self,
        previous: &P1Interfaces<'_>,
        next: &P1Interfaces<'_>,
        cancelled: &mut C,
    ) -> Result<(), FractureTopologyError> {
        require_reference_map(self.previous(), previous, cancelled)?;
        require_reference_map(self.topology(), next, cancelled)?;
        for (index, old) in previous.points().iter().enumerate() {
            poll(cancelled)?;
            let new = &next.points()[index];
            if old.interface_id != new.interface_id
                || old.point_in_face != new.point_in_face
                || old.location_m != new.location_m
                || old.area_weight_m2 != new.area_weight_m2
                || old.axes != new.axes
                || old.shape_values != new.shape_values
                || old.side_b_shape_values != new.side_b_shape_values
                || (0..2).any(|side| {
                    (0..3).any(|i| {
                        self.previous_nodes()[new.paired_nodes[side][i]]
                            != old.paired_nodes[side][i]
                    })
                })
            {
                return Err(FractureTopologyError::Invalid(
                    "retained interface physical point or local frame changed",
                ));
            }
            let old_patch = &previous.patches()[index];
            let new_patch = &next.patches()[index];
            let same_source = match (old_patch.source(), new_patch.source()) {
                (None, None) => true,
                (Some(a), Some(b)) => core::ptr::eq(a, b),
                _ => false,
            };
            if !same_source || !core::ptr::eq(old_patch.law(), new_patch.law()) {
                return Err(FractureTopologyError::Invalid(
                    "retained interface material/source owner changed",
                ));
            }
        }
        poll(cancelled)
    }
}

fn require_reference_map<C: FnMut() -> bool>(
    topology: &TetFractureTopology<'_>,
    map: &P1Interfaces<'_>,
    cancelled: &mut C,
) -> Result<(), FractureTopologyError> {
    poll(cancelled)?;
    let count = topology
        .seams()
        .len()
        .checked_mul(3)
        .ok_or(FractureTopologyError::Budget("retained P1 point count"))?;
    if map.full_dimension() != 3 * topology.nodes_m().len()
        || map.points().len() != count
        || map.patches().len() != count
    {
        return Err(FractureTopologyError::Invalid(
            "P1 map must cover the complete corresponding split topology",
        ));
    }
    for (index, seam) in topology.seams().iter().enumerate() {
        poll(cancelled)?;
        let opposite = |side: usize, face: [usize; 3]| {
            topology.tetrahedra()[seam.elements[side]]
                .iter()
                .copied()
                .find(|node| !face.contains(node))
                .ok_or(FractureTopologyError::Invalid(
                    "missing opposite seam vertex",
                ))
        };
        // Element/face quality was already admitted by both geometry owners.
        // Reuse their exact P1 frame arithmetic with only nondegeneracy here.
        let (axes, area) = frame(
            topology.nodes_m(),
            seam.side_a,
            opposite(0, seam.side_a)?,
            opposite(1, seam.side_b)?,
            0.0,
        )
        .map_err(|_| FractureTopologyError::Invalid("P1 seam reference frame"))?;
        for q in 0..3 {
            let point = &map.points()[3 * index + q];
            let patch = &map.patches()[3 * index + q];
            let shape = core::array::from_fn(|i| if i == q { 2.0 / 3.0 } else { 1.0 / 6.0 });
            let mut location = [0.0; 3];
            for c in 0..3 {
                for i in 0..3 {
                    location[c] += shape[i] * topology.nodes_m()[seam.side_a[i]][c];
                }
            }
            if point.interface_id != seam.id
                || point.point_in_face != q
                || point.paired_nodes != [seam.side_a, seam.side_b]
                || point.shape_values != shape
                || point.side_b_shape_values != shape
                || point.location_m != location
                || point.axes != axes
                || point.area_weight_m2 != area / 3.0
                || patch.id() != 3 * seam.id + q as u64
                || patch.area_m2() != point.area_weight_m2
            {
                return Err(FractureTopologyError::Invalid(
                    "P1 map does not retain the seam's ordered reference geometry",
                ));
            }
        }
    }
    poll(cancelled)
}

fn poll<C: FnMut() -> bool>(cancelled: &mut C) -> Result<(), FractureTopologyError> {
    if cancelled() {
        Err(FractureTopologyError::Cancelled)
    } else {
        Ok(())
    }
}
