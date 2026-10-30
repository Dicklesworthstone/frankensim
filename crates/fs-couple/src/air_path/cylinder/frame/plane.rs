//! Prescribed exterior regions from the original reference geometry.
use super::*;

impl<F> FrameFluidRegion<F> {
    /// Partition one original member at `normal . reference_position = offset_m`.
    /// The normal is a finite unit vector in the SAME planar coordinates as the
    /// source frame; offset is metres. Negative/positive sides receive exactly
    /// the caller's respective material objects. Sources are neither evaluated
    /// nor cloned. Returns at most two regions in original i-to-j order.
    ///
    /// A boundary endpoint contributes zero measure, so an otherwise one-sided
    /// member produces one whole span. A member lying wholly on the plane refuses
    /// as ambiguous. Invalid geometry, non-unit normals and an unrepresentable
    /// interior crossing refuse; no epsilon thickens or deletes a physical span.
    /// Final model construction separately admits quadrature and resource limits.
    ///
    /// This partitions a REFERENCE CENTERLINE, not a circular surface intersection
    /// or a solved liquid interface. Linear reference geometry is required.
    /// Surface waves, buoyancy, capillary forces and interface hydrodynamics are
    /// not added. This is cold preparation, not moving-boundary history transfer.
    #[allow(clippy::too_many_arguments)]
    pub fn split_at_reference_plane(
        frame: &Frame2d,
        element: usize,
        normal: [f64; 2],
        offset_m: f64,
        negative: F,
        positive: F,
        gate: &CancelGate,
    ) -> Result<Vec<Self>> {
        poll(gate)?;
        if normal.iter().any(|x| !x.is_finite()) || !offset_m.is_finite()
            || (normal[0].hypot(normal[1])-1.0).abs() > 1e-12 {
            return Err(invalid("fluid plane requires a finite unit normal and offset in metres"));
        }
        let points = frame.reference_member_segment(element)?;
        // Admit the SAME linear-reference approximation as the actual loading.
        frame.reference_transverse_trace(element, 0.5)?;
        let distance = |point: [f64; 2]| finite(normal[0].mul_add(point[0], normal[1]*point[1])-offset_m);
        let (a, b) = (distance(points[0])?, distance(points[1])?);
        let regions = if a == 0.0 && b == 0.0 {
            return Err(invalid("member lies on the fluid boundary; explicit source coverage required"));
        } else if a <= 0.0 && b <= 0.0 {
            vec![Self { element, span: [0.0, 1.0], fluid: negative }]
        } else if a >= 0.0 && b >= 0.0 {
            vec![Self { element, span: [0.0, 1.0], fluid: positive }]
        } else {
            // Opposite signs: scale before summing magnitudes. No overflowing
            // a-b, clipping, artificial gap or overlapping side interval.
            let scale = a.abs().max(b.abs());
            let left = a.abs()/scale;
            let cut = left/(left+b.abs()/scale);
            if !cut.is_finite() || cut <= 0.0 || cut >= 1.0 {
                return Err(invalid("fluid-plane crossing is not a representable interior coordinate"));
            }
            let (first, second) = if a < 0.0 { (negative, positive) } else { (positive, negative) };
            vec![Self { element, span: [0.0, cut], fluid: first },
                Self { element, span: [cut, 1.0], fluid: second }]
        };
        poll(gate)?;
        Ok(regions)
    }
}

#[cfg(test)]
mod tests;
