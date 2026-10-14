//! Attach the exact in-crate source-law correction to a checked endpoint.
use super::*;
use crate::interface::temperature::tangent::ContactTemperatureLinearization;

impl super::super::RobinResponse {
    pub(crate) fn with_contact_temperature_feedback(
        mut self, cx: &Cx<'_>, contact: &ContactTemperatureLinearization,
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        if self.temperature.iter().map(|v| v.to_bits())
            .ne(contact.temperature().iter().map(|v| v.to_bits()))
            || self.linear.restart == 0
        {
            return Err(invalid("contact feedback must match the checked endpoint and admit a positive FGMRES restart"));
        }
        let (extra, _) = reduce_matrix_and_lift(contact.feedback(), &self.dofs);
        let n = self.matrix.nrows();
        if extra.nrows() != n { return Err(invalid("contact feedback has a different free-node layout")); }
        let mut coo = fs_sparse::Coo::new(n, n);
        for row in 0..n {
            poll(cx, row)?;
            for matrix in [&self.matrix, &extra] {
                let (columns, values) = matrix.row(row);
                for (&column, &value) in columns.iter().zip(values) {
                    coo.push(row, column, checked(value)?);
                }
            }
        }
        let matrix = coo.assemble();
        for row in 0..n {
            poll(cx, row)?;
            for &value in matrix.row(row).1 { checked(value)?; }
        }
        // R'(T) makes contact nonsymmetric even with a constant bulk tensor.
        // Preserve the old material smoothness gate; do not fall back to PCG.
        let smooth = self.has_smooth_material_tangent();
        self.nonlinear = Some(TangentSystem { transpose: ops::transpose(&matrix), smooth });
        self.matrix = matrix;
        poll(cx, n)?;
        Ok(self)
    }
}
