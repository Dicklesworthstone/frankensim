//! Geometry derivatives of the original static planar magnetic potential.
//!
//! Under an in-plane P1 deformation F, the transported A_z gives
//! B -> F B / det(F). Thus -dW at fixed potential coordinates is
//! `W'_cell div(u) - volume H dot (grad(u) B)`, with W'=B.H-W.
//! Stationarity removes the field-state derivative when integrated impressed
//! currents and original Dirichlet A_z values are held fixed. This is a
//! derivative of the native discrete potential, not a reconstructed stress.
//!
//! Energy/control conventions: MIT 6.685, Class Notes 1, sections 4.1–4.2,
//! https://web.mit.edu/6.685/www/chapter1.pdf .

use super::shape::shape_force;
use crate::planar::{
    MagneticError, PlanarField, gradients,
    nonlinear::{IsotropicMagneticMaterial, NonlinearMagneticError, PlanarNonlinearField},
    poll,
};
use fs_math::dd::Dd;
use fs_rep_mesh::{Metric2, TriComplex2};

fn finite(value: f64) -> Result<f64, MagneticError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(MagneticError::Invalid(
            "unrepresentable magnetic virtual work",
        ))
    }
}

fn integrate(
    mesh: &TriComplex2,
    displacement_m: &[[f64; 2]],
    cell_coenergy_j: &[f64],
    magnetic_field_a_m: &[[f64; 2]],
    magnetic_flux_density_t: &[[f64; 2]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, MagneticError> {
    poll(cancelled)?;
    if !matches!(mesh.metric(), Metric2::Planar { .. })
        || displacement_m.len() != mesh.vertices().len()
        || [
            cell_coenergy_j.len(),
            magnetic_field_a_m.len(),
            magnetic_flux_density_t.len(),
        ]
        .into_iter()
        .any(|length| length != mesh.faces().len())
    {
        return Err(MagneticError::Invalid(
            "magnetic virtual work requires the complete original planar field and displacement",
        ));
    }
    for displacement in displacement_m {
        poll(cancelled)?;
        for &value in displacement {
            finite(value)?;
        }
    }
    let mut work = Dd::from_f64(0.0);
    for (cell, face) in mesh.faces().iter().enumerate() {
        poll(cancelled)?;
        let gradients = gradients(face.map(|v| mesh.vertices()[v as usize]))?;
        let volume = mesh.face_measure(cell).ok_or(MagneticError::Invalid(
            "missing original magnetic cell volume",
        ))?;
        // Shape gradients sum to zero. Contract only local displacement
        // differences so rigid translation vanishes before summation.
        for local in 1..3 {
            let force = shape_force(
                cell_coenergy_j[cell],
                volume,
                magnetic_field_a_m[cell],
                magnetic_flux_density_t[cell],
                gradients[local],
            )
            .ok_or(MagneticError::Invalid(
                "unrepresentable magnetic shape force",
            ))?;
            for axis in 0..2 {
                let difference = finite(
                    displacement_m[face[local] as usize][axis]
                        - displacement_m[face[0] as usize][axis],
                )?;
                work = work + Dd::from_f64(force[axis]) * Dd::from_f64(difference);
                finite(work.to_f64())?;
            }
        }
    }
    poll(cancelled)?;
    finite(work.to_f64())
}

impl PlanarField<'_> {
    /// Magnetic mechanical virtual work [J] along `x(s)=x+s*displacement_m`.
    /// The dimensionless path parameter, complete P1 displacement in original
    /// vertex order, and fixed extrusion thickness define the geometry change.
    /// A rigid in-plane rotation tangent on a body, extended through the
    /// surrounding mesh, gives its axial torque times the virtual angle.
    ///
    /// Hold each integrated cell current `J_z * area` and the original
    /// Dirichlet A_z values fixed while transporting the cells. Material
    /// coefficients/state are fixed. The result is `-d(Pi)/ds`, where
    /// `Pi = sum(W_cell) - sum_cells(J_z*volume*mean(A_z))` is the actual
    /// stationary discrete magnetic potential. At zero Dirichlet A_z this
    /// is `+d(coenergy)/ds` at fixed currents; without impressed currents it
    /// is `-d(energy)/ds` at fixed boundary flux. No field re-solve is needed.
    ///
    /// Only original static induced-linear fields are admitted. Recoil
    /// magnets and fields recovered from an induction history refuse; their
    /// material/source motion and work are not inferred. The caller owns the
    /// physical displacement extension and any finite geometry evolution.
    /// This supplies no magnetostriction, finite-step energy integrator or
    /// field-error enclosure. Do not add overlapping Maxwell/Lorentz work.
    /// The query allocates no heap storage, polls original vertices/cells,
    /// and publishes no work on failure or cancellation.
    pub fn virtual_work(
        &self,
        displacement_m: &[[f64; 2]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, MagneticError> {
        poll(cancelled)?;
        if self.impressed_current_density_a_m2().is_none() || self.problem().has_permanent_magnets()
        {
            return Err(MagneticError::Invalid(
                "magnetic virtual work requires original static induced-linear source controls",
            ));
        }
        integrate(
            self.problem().mesh(),
            displacement_m,
            self.cell_energy_j(),
            self.magnetic_field_a_m(),
            self.magnetic_flux_density_t(),
            cancelled,
        )
    }
}

impl PlanarNonlinearField<'_> {
    /// Native magnetic virtual work for static scalar reversible materials.
    /// Geometry, held integrated currents and Dirichlet A_z controls have the
    /// same meaning as [`PlanarField::virtual_work`]. This uses the original
    /// nonlinear cell coenergy `B.H-W`, which generally differs from storage,
    /// and the original constitutive H; no secant-linear replacement is made.
    /// Recoil magnets, hysteretic and induction-recovered histories refuse.
    pub fn virtual_work(
        &self,
        displacement_m: &[[f64; 2]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, NonlinearMagneticError> {
        poll(cancelled)?;
        if self.impressed_current_density_a_m2().is_none() {
            return Err(MagneticError::Invalid(
                "nonlinear magnetic virtual work requires original static source controls",
            )
            .into());
        }
        for material in self.problem().materials() {
            poll(cancelled)?;
            if matches!(material, IsotropicMagneticMaterial::Recoil(_)) {
                return Err(MagneticError::Invalid(
                    "recoil-material motion is outside magnetic virtual work",
                )
                .into());
            }
        }
        Ok(integrate(
            self.problem().mesh(),
            displacement_m,
            self.cell_coenergy_j(),
            self.magnetic_field_a_m(),
            self.magnetic_flux_density_t(),
            cancelled,
        )?)
    }
}
