//! Exact discrete virtual work from original tetrahedral dielectric fields.

use super::{Result, SpatialDielectricLaw, SpatialElectricField, add, check, finite, invalid};
use crate::electrostatic::virtual_work::kernel;

impl SpatialElectricField<'_> {
    /// Electric mechanical virtual work [J] along `x(t) = x + t * displacement_m`,
    /// with dimensionless `t` and one displacement per original mesh vertex.
    ///
    /// The analytic contraction is
    /// `sum_cells(U_cell div(u) - volume E_i D_j d_j(u_i))`, using the accepted
    /// field, original tetrahedral shape gradients and physical cell volumes.
    /// It needs no field re-solve, stress recovery or interface averaging. Local
    /// displacement differences make a common translation cancel exactly.
    ///
    /// For fixed electrode voltages this equals `dU/dt`; the voltage sources
    /// supply the complementary work. For fixed electrode charges it equals
    /// `-dU/dt` after electrostatic re-equilibration. A unit-voltage two-terminal
    /// field therefore supplies `dC/dt = 2 * virtual_work`.
    ///
    /// Scalar dielectric coefficients and their complete material states stay
    /// fixed in each transported cell. Explicitly anisotropic laws and nonzero
    /// free volume/sheet charge, imposed normal displacement or impressed
    /// polarization refuse until their deformation/source work is supplied.
    /// Electrostriction and mechanical pressure are not inferred from epsilon.
    ///
    /// The caller owns the physical displacement map, including the exterior
    /// field domain. Mesh-coordinate work is not automatically a solid nodal
    /// traction or finite mechanical trajectory work. Accuracy inherits the
    /// accepted field residual and spatial approximation. Do not add another
    /// overlapping electric load. Refusal/cancellation leaves the field intact.
    pub fn virtual_work(
        &self,
        displacement_m: &[[f64; 3]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64> {
        check(cancelled)?;
        let g = &self.problem.geometry;
        if displacement_m.len() != g.positions_m.len() {
            return Err(invalid(
                "virtual work requires complete original mesh displacements",
            ));
        }
        for displacement in displacement_m {
            check(cancelled)?;
            for &value in displacement {
                finite(value)?;
            }
        }
        if self.has_imposed_displacement {
            return Err(invalid(
                "electric virtual work requires zero imposed displacement",
            ));
        }
        for sheet in &self.sheets {
            check(cancelled)?;
            if sheet.free_charge_density_c_m2 != 0.0 {
                return Err(invalid(
                    "electric virtual work requires zero free sheet charge",
                ));
            }
        }
        let (mut work, mut correction) = (0.0, 0.0);
        for (cell, tet) in g.mesh.tets.iter().enumerate() {
            check(cancelled)?;
            if !matches!(self.problem.laws[cell], SpatialDielectricLaw::Isotropic)
                || self.free_charge_density_c_m3[cell] != 0.0
                || self
                    .impressed_polarization_c_m2
                    .as_ref()
                    .is_some_and(|p| p[cell].iter().any(|&value| value != 0.0))
            {
                return Err(invalid(
                    "electric virtual work requires explicitly isotropic source-free cells",
                ));
            }
            for local in 1..4 {
                let force = kernel::shape_force(
                    self.cell_energy_j[cell],
                    g.geometry.vol_signed[cell].abs(),
                    self.electric_field_v_m[cell],
                    self.displacement_c_m2[cell],
                    g.geometry.grads[cell][local],
                )
                .ok_or_else(|| invalid("unrepresentable electric shape force"))?;
                for axis in 0..3 {
                    let difference = finite(
                        displacement_m[tet[local] as usize][axis]
                            - displacement_m[tet[0] as usize][axis],
                    )?;
                    add(
                        &mut work,
                        &mut correction,
                        finite(force[axis] * difference)?,
                    )?;
                }
            }
        }
        check(cancelled)?;
        finite(work)
    }
}
