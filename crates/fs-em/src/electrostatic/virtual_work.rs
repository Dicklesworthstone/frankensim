//! Geometry derivatives of the original source-free dielectric field energy.

pub(crate) mod kernel;

use super::{
    ElectricField, ElectrostaticError, PlanarDielectricLaw, accumulate, finite, geometry, poll,
};
use crate::planar::gradients;
use fs_rep_mesh::Metric2;

impl ElectricField<'_> {
    /// Electric mechanical virtual work [J] along `x(t) = x + t * displacement_m`,
    /// where `t` is dimensionless and the original extrusion thickness is fixed.
    /// The complete displacement follows the original mesh vertex order.
    ///
    /// This analytically differentiates the SAME discrete field potential,
    /// using its original cell fields, shape gradients and physical volumes:
    /// `sum_cells(U_cell div(u) - volume E_i D_j d_j(u_i))`.
    /// No perturbed field solve, stress reconstruction or nodal averaging is
    /// used. Local displacement differences preserve exact invariance under a
    /// common translation, including a large displacement offset.
    ///
    /// With held electrode voltages this equals the stored-energy derivative
    /// `dU/dt`; the electrical sources supply both the storage change and this
    /// mechanical work. With held electrode charges the same mechanical work equals
    /// `-dU/dt` after the field is re-equilibrated. In particular a unit-voltage
    /// two-terminal field gives `dC/dt = 2 * virtual_work`.
    ///
    /// Every cell must have an explicitly isotropic law with its resolved
    /// coefficient and material state held fixed. Nonzero free volume/sheet
    /// charge, imposed normal displacement and impressed polarization refuse:
    /// their transport/source work is not inferred. Neither electrostriction
    /// nor a constitutive mechanical pressure follows from permittivity alone.
    ///
    /// This is an infinitesimal derivative of the represented mesh, not finite
    /// trajectory work or a solid nodal-traction map. The caller must define
    /// the physical geometry path, including the surrounding field domain.
    /// Accuracy inherits the accepted field residual and spatial discretization.
    /// Do not add an overlapping Maxwell or other electric load to this work.
    /// Failure/cancellation publishes no work and leaves the field unchanged.
    pub fn virtual_work(
        &self,
        displacement_m: &[[f64; 2]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, ElectrostaticError> {
        poll(cancelled)?;
        let mesh = self.problem.mesh();
        if !matches!(mesh.metric(), Metric2::Planar { .. })
            || displacement_m.len() != mesh.vertices().len()
        {
            return Err(ElectrostaticError::Invalid(
                "virtual work requires complete displacements on the original planar mesh",
            ));
        }
        for displacement in displacement_m {
            poll(cancelled)?;
            for &value in displacement {
                finite(value)?;
            }
        }
        for sheet in &self.free_charge_sheets {
            poll(cancelled)?;
            if sheet.free_charge_density_c_m2 != 0.0 {
                return Err(ElectrostaticError::Invalid(
                    "electric virtual work requires zero free sheet charge",
                ));
            }
        }
        for trace in &self.displacement_traces {
            poll(cancelled)?;
            if trace.outward_displacement_c_m2 != 0.0 {
                return Err(ElectrostaticError::Invalid(
                    "electric virtual work requires zero imposed displacement",
                ));
            }
        }
        let (mut work, mut correction) = (0.0, 0.0);
        for (cell, face) in mesh.faces().iter().enumerate() {
            poll(cancelled)?;
            if !matches!(
                self.problem.permittivity_laws[cell],
                PlanarDielectricLaw::Isotropic
            ) || self.free_charge_density_c_m3[cell] != 0.0
                || self
                    .impressed_polarization_c_m2
                    .as_ref()
                    .is_some_and(|p| p[cell].iter().any(|&value| value != 0.0))
            {
                return Err(ElectrostaticError::Invalid(
                    "electric virtual work requires explicitly isotropic source-free cells",
                ));
            }
            let gradients =
                gradients(face.map(|v| mesh.vertices()[v as usize])).map_err(geometry)?;
            let volume = mesh.face_measure(cell).ok_or(ElectrostaticError::Invalid(
                "missing original dielectric cell volume",
            ))?;
            for local in 1..3 {
                let force = kernel::shape_force(
                    self.cell_energy_j[cell],
                    volume,
                    self.electric_field_v_m[cell],
                    self.displacement_c_m2[cell],
                    gradients[local],
                )
                .ok_or(ElectrostaticError::Invalid(
                    "unrepresentable electric shape force",
                ))?;
                for axis in 0..2 {
                    let difference = finite(
                        displacement_m[face[local] as usize][axis]
                            - displacement_m[face[0] as usize][axis],
                    )?;
                    accumulate(
                        &mut work,
                        &mut correction,
                        finite(force[axis] * difference)?,
                    )?;
                }
            }
        }
        poll(cancelled)?;
        finite(work + correction)
    }
}
