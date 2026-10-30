//! Original reference geometry and power-conjugate transverse trace.
use super::{Frame2d,SolidError,invalid};
use crate::frame2d::Geometry;

impl Frame2d {
    /// Original reference endpoints [m], in the element's own i-to-j order.
    /// This is the SAME geometry used by reference inertia and transverse
    /// traces, not a second authored segment or the current deformed chord.
    /// The read-only query supports corotational frames too; consumers must
    /// independently admit their reference/current-geometry approximation.
    pub fn reference_member_segment(&self, element: usize) -> Result<[[f64; 2]; 2], SolidError> {
        let member = self.elements.get(element)
            .ok_or_else(|| invalid("reference member does not exist"))?;
        let points = [self.nodes[member.i], self.nodes[member.j]];
        if points.iter().flatten().any(|v| !v.is_finite())
            || !member.l0.is_finite() || member.l0 <= 0.0 {
            return Err(invalid("reference member geometry is not finite and positive"));
        }
        Ok(points)
    }

    /// Work-conjugate reference transverse trace at member coordinate xi in [0,1].
    ///
    /// Returns original node-major DOFs, their six displacement weights and
    /// reference member length [m]. The dot product with displacement gives
    /// transverse displacement [m]; with velocity it gives m/s. Multiplying
    /// these SAME weights by a transverse point force returns nodal forces and
    /// moments, so virtual work is unchanged. Rotation weights have metres.
    ///
    /// This is the cubic Hermite interpolation of `ConsistentReference` inertia,
    /// along the member's reference normal (-sin,cos). The selected dry inertia
    /// formulation is not changed. Quadrature weights/areas belong to the caller;
    /// this method supplies neither fluid forces nor a material law. Only linear
    /// reference geometry is admitted; finite-rotation follower loads need their
    /// own current-geometry trace. No allocation or accepted-state change occurs.
    pub fn reference_transverse_trace(
        &self,
        element: usize,
        xi: f64,
    ) -> Result<([usize; 6], [f64; 6], f64), SolidError> {
        if self.geometry != Geometry::Linear || !xi.is_finite()
            || !(0.0..=1.0).contains(&xi)
        {
            return Err(invalid("transverse trace requires linear reference geometry and xi in [0,1]"));
        }
        let e = self.elements.get(element)
            .ok_or_else(|| invalid("transverse trace member does not exist"))?;
        if !e.l0.is_finite() || e.l0 <= 0.0 || !e.cos0.is_finite() || !e.sin0.is_finite() {
            return Err(invalid("transverse trace reference geometry is not finite"));
        }
        let one_minus = 1.0 - xi;
        let h0 = one_minus * one_minus * (1.0 + 2.0 * xi);
        let h1 = e.l0 * xi * one_minus * one_minus;
        let h2 = xi * xi * (3.0 - 2.0 * xi);
        let h3 = -e.l0 * xi * xi * one_minus;
        let weights = [-e.sin0*h0, e.cos0*h0, h1, -e.sin0*h2, e.cos0*h2, h3];
        if weights.iter().any(|x| !x.is_finite()) {
            return Err(invalid("transverse trace weights are unrepresentable"));
        }
        Ok(([3*e.i, 3*e.i+1, 3*e.i+2, 3*e.j, 3*e.j+1, 3*e.j+2], weights, e.l0))
    }
}
