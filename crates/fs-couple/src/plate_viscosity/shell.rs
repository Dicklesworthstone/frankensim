//! Material membrane/bending loss of the original CST/DKT shell.
//!
//! The spatial owner contracts physical strain rates before projecting its
//! passive resultant tensor. This adapter retains every cross-modal term and
//! realizes the result through the existing modal dashpot owner. It neither
//! damps the shell's artificial drilling stiffness nor changes its modes.

use super::{ModalPlateViscosity, PlateViscosityLimits, invalid, owner, poll};
use crate::acoustic_realize::AcousticRealizeError;
use fs_exec::CancelGate;
use fs_plate::shell::viscosity::{
    ShellViscosity, ShellViscosityBudget, project_viscosity, projection_terms,
};
use fs_plate::shell::{ShellMesh, ShellModel};

impl ModalPlateViscosity {
    /// Project physical CST membrane and DKT bending viscosity on the caller's
    /// unchanged reduced shell basis, including membrane/bending coupling.
    ///
    /// `modes` must be finite, mass-orthonormal vectors in `model.dof_map`
    /// order. Their Gram matrix is checked against the ORIGINAL sparse mass,
    /// including consistent and translation/director coupling terms, to an
    /// absolute entry tolerance of 1e-8. No normalization or new eigensolve is
    /// performed. The returned matrix has units 1/s and its passive channels
    /// can enter [`Self::connections`] or the existing certified block path.
    ///
    /// Supply one viscosity section per original facet, in its local tangent
    /// axes. The mesh, pencil and basis must come from the same CST/DKT model;
    /// the spatial owner's map checks do not establish geometry identity for a
    /// public `ShellModel`. MITC, tied and obliquely constrained shell pencils
    /// have their own projection owners. Constitutive applicability, source
    /// bands and modal completeness remain caller obligations. Support,
    /// radiation, and other material losses must not be counted again here.
    ///
    /// `max_projection_terms` includes spatial projection, the sparse mass
    /// products, all Gram entries and passive factorization. Preparation polls
    /// cancellation; accepted propagation uses the existing prepared solver.
    /// There is no material state or heat update in this cold projection.
    ///
    /// # Errors
    /// Invalid dimensions, nonfinite or non-orthonormal modes, exhausted work,
    /// cancellation, or the original constitutive/projection/factor refusals.
    pub fn from_cst_dkt_shell(
        mesh: &ShellMesh,
        viscosities: &[ShellViscosity],
        model: &ShellModel,
        modes: &[Vec<f64>],
        limits: PlateViscosityLimits,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        let n = modes.len();
        let full_dofs = mesh
            .nodes
            .len()
            .checked_mul(6)
            .ok_or_else(|| invalid("shell viscosity coordinate count overflows"))?;
        if !(1..=64).contains(&limits.max_modes)
            || n == 0
            || n > limits.max_modes
            || model.free == 0
            || model.free > full_dofs
            || model.m.nrows() != model.free
            || model.m.ncols() != model.free
            || modes.iter().any(|mode| mode.len() != model.free)
            || !limits.reference_resistance_n_s_m.is_finite()
            || limits.reference_resistance_n_s_m <= 0.0
        {
            return Err(invalid(
                "shell viscosity needs its original bounded mass-normalized basis",
            ));
        }
        let spatial_work = projection_terms(full_dofs, mesh.tris.len(), n)
            .ok_or_else(|| invalid("shell viscosity projection work overflows"))?;
        let basis_work = model
            .m
            .nnz()
            .checked_mul(n)
            .and_then(|work| {
                model
                    .free
                    .checked_mul(n * n)
                    .and_then(|gram| work.checked_add(gram))
            })
            .and_then(|work| work.checked_add(4 * n * n * n))
            .ok_or_else(|| invalid("shell viscosity basis work overflows"))?;
        if spatial_work
            .checked_add(basis_work)
            .is_none_or(|work| work > limits.max_projection_terms)
        {
            return Err(invalid(
                "shell viscosity exceeds its declared projection work",
            ));
        }
        for mode in modes {
            for chunk in mode.chunks(256) {
                poll(gate)?;
                if chunk.iter().any(|value| !value.is_finite()) {
                    return Err(invalid("shell viscosity requires finite modal vectors"));
                }
            }
        }
        let mut mass_mode = vec![0.0; model.free];
        for (i, mode) in modes.iter().enumerate() {
            poll(gate)?;
            model.m.spmv(mode, &mut mass_mode);
            for (j, other) in modes.iter().enumerate() {
                let mut sum = 0.0;
                let mut correction = 0.0;
                for (chunk, weighted) in other.chunks(256).zip(mass_mode.chunks(256)) {
                    poll(gate)?;
                    for (&value, &mass_value) in chunk.iter().zip(weighted) {
                        let term = value * mass_value;
                        let next = sum + term;
                        correction += if sum.abs() >= term.abs() {
                            (sum - next) + term
                        } else {
                            (term - next) + sum
                        };
                        sum = next;
                    }
                }
                let gram = sum + correction;
                let expected = if i == j { 1.0 } else { 0.0 };
                if !gram.is_finite() || (gram - expected).abs() > 1e-8 {
                    return Err(invalid(
                        "shell viscosity basis is not orthonormal in the original mass",
                    ));
                }
            }
        }
        let matrix = project_viscosity(
            mesh,
            viscosities,
            model,
            modes,
            ShellViscosityBudget {
                max_elements: mesh.tris.len(),
                max_full_dofs: full_dofs,
                max_dofs: model.free,
                max_modes: limits.max_modes,
                max_work: limits.max_projection_terms - basis_work,
            },
            || gate.is_requested(),
        )
        .map_err(|error| match error {
            fs_plate::PlateError::Cancelled => AcousticRealizeError::Cancelled,
            other => owner(other),
        })?;
        Self::factor(matrix, n, limits.reference_resistance_n_s_m, gate)
    }
}

#[cfg(test)]
mod tests;
