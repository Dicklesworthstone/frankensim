//! Actual physical shell rotation, conjugate to torque in the common body basis.
use super::*;
use super::super::torsion::ModalRotationAttachment;
use fs_plate::shell::point::{triangle_point_rows,triangle_rotation_rows};

/// A local material rotation in the original shell/assembly geometry.
#[derive(Clone, Debug, PartialEq)]
pub struct ShellRotationGeometry {
    /// Original assembly component, or zero for one material shell.
    pub shell_component: usize,
    /// Original triangle index within that component.
    pub triangle: usize,
    /// Nonnegative unit-sum weights in the original triangle's node order.
    pub barycentric: [f64; 3],
    /// Signed unit global XYZ axis; the connected string twist uses this sign.
    pub axis: [f64; 3],
}

impl ShellRotationGeometry {
    /// Admit the explicit rotation coordinate without normalizing its inputs.
    /// Component/triangle existence belongs to the original source query.
    pub fn validate(&self) -> Result<(), ReciprocalStringError> {
        triangle_point_rows(self.barycentric,self.axis,[[0.0;3];3])
            .map(|_|()).map_err(|e|invalid(e.reason()))
    }
}

impl MaterialStringBankSource {
    /// Project local physical orientation at one original shell triangle.
    /// The returned row observes angular velocity and returns its conjugate
    /// torque through `StringBodyRotationConnection` in the existing solve.
    ///
    /// Every assembly component is evaluated on the SAME assembled body basis
    /// at original body index zero. Material tensors, inertia, supports and
    /// mortar ties already determine director motion and membrane spin. This query
    /// introduces no new modes, baseline, angular gain or joint constitutive law.
    /// Tangential orientation comes from the original interpolated directors;
    /// normal spin comes from the original P1 membrane curl. Every independent
    /// director-axis drilling gauge is annihilated. Numerical drilling inertia
    /// and penalty therefore cannot act as this port's physical angular motion.
    /// This P1 interface query does not reconstruct the condensed interior bubble.
    ///
    /// # Errors
    /// Absent source/component/triangle, invalid axis or interpolation, detached
    /// retained rotation, unrepresentable modal projection, or cancellation.
    pub fn shell_rotation(&self, geometry: &ShellRotationGeometry, gate: &CancelGate)
        -> Result<ModalRotationAttachment, ReciprocalStringError>
    {
        poll(gate)?;
        geometry.validate()?;
        let (specimen,directors,modes) = self.physical_shell_component(geometry.shell_component)?;
        let triangle = specimen.mesh().tris.get(geometry.triangle)
            .ok_or_else(||invalid("shell rotation names an absent original triangle"))?;
        let facet = specimen.mesh().facet(geometry.triangle)
            .map_err(|_|invalid("shell rotation lost its original admitted facet geometry"))?;
        let rows = triangle_rotation_rows(geometry.barycentric,geometry.axis,facet.frame,facet.gradient,
            triangle.map(|node|directors.vectors()[node])).map_err(|e|invalid(e.reason()))?;
        let angular_shapes = modes.iter().map(|mode| {
            poll(gate)?;
            let mut angle = 0.0;
            for a in 0..3 {
                let start = triangle[a].checked_mul(6)
                    .ok_or_else(||invalid("shell rotation node extent overflows"))?;
                let nodal = mode.get(start..start+6)
                    .ok_or_else(||invalid("shell rotation basis lost its original physical node"))?;
                angle += nodal.iter().zip(rows[a]).map(|(q,torque)|q*torque).sum::<f64>();
            }
            if !angle.is_finite() {
                return Err(invalid("physical shell rotation projection overflows"));
            }
            Ok(angle)
        }).collect::<Result<Vec<_>,ReciprocalStringError>>()?;
        if angular_shapes.is_empty() || angular_shapes.iter().all(|v|*v == 0.0) {
            return Err(invalid("physical shell rotation is detached from the retained body basis"));
        }
        poll(gate)?;
        Ok(ModalRotationAttachment { component: 0, angular_shapes })
    }
}
