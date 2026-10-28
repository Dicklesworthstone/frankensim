//! Reuse immutable reference geometry throughout a finite-deformation trajectory.
//!
//! Shape gradients and reference volumes depend on the reference mesh, not on
//! its current displacement. Preparation admits them once; every evaluation
//! still checks the current deformation and the original material domain.
//! Borrowed mesh and law data cannot change while this preparation is alive.
//! Constitutive evaluation and assembly remain owned by the parent hyper3
//! kernel, with the same deterministic element order and cancellation points.

use fs_exec::Cx;
use fs_material::hyper::Hyperelastic;
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_sparse::Csr;

use super::{Element, Evaluated, HyperTetError, HyperTetProblem, checkpoint};

/// An admitted immutable mesh and material binding for repeated evaluations.
///
/// Retains O(number of tetrahedra) reference gradients and volumes. No tangent,
/// current force, current deformation or constitutive history is cached.
#[derive(Debug)]
pub struct PreparedHyperTetProblem<'a, L: HyperelasticMaterial = Hyperelastic> {
    problem: HyperTetProblem<'a, L>,
    elements: Vec<Element>,
}

impl<'a, L: HyperelasticMaterial> HyperTetProblem<'a, L> {
    /// Admit the reference problem once for repeated dynamic evaluations.
    ///
    /// # Errors
    /// The same mesh, material-field, load, support and budget refusals as the
    /// native problem, or cancellation before/after reference preparation.
    pub fn prepare(&self, cx: &Cx<'_>) -> Result<PreparedHyperTetProblem<'a, L>, HyperTetError> {
        checkpoint(cx)?;
        let elements = self.validate()?;
        checkpoint(cx)?;
        Ok(PreparedHyperTetProblem {
            problem: *self,
            elements,
        })
    }
}

impl<L: HyperelasticMaterial> PreparedHyperTetProblem<'_, L> {
    fn evaluate(&self, displacement_m: &[f64], cx: &Cx<'_>) -> Result<Evaluated, HyperTetError> {
        checkpoint(cx)?;
        self.problem.check_len(displacement_m)?;
        let state = self
            .problem
            .evaluate(&self.elements, displacement_m, cx)?
            .map_err(|(element, what)| HyperTetError::MaterialRefused { element, what })?;
        if !state.energy.is_finite() || state.internal.iter().any(|value| !value.is_finite()) {
            return Err(HyperTetError::InvalidInput {
                what: "stored energy or internal force is not finite".into(),
            });
        }
        checkpoint(cx)?;
        Ok(state)
    }

    /// Full native internal force, N, at the current displacement.
    ///
    /// # Errors
    /// Invalid displacement, material-domain refusal, overflow or cancellation.
    pub fn internal_force(
        &self,
        displacement_m: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, HyperTetError> {
        Ok(self.evaluate(displacement_m, cx)?.internal)
    }

    /// Native stored strain energy, J, with full current-state admission.
    ///
    /// # Errors
    /// Invalid displacement, material-domain refusal, overflow or cancellation.
    pub fn stored_energy_j(
        &self,
        displacement_m: &[f64],
        cx: &Cx<'_>,
    ) -> Result<f64, HyperTetError> {
        Ok(self.evaluate(displacement_m, cx)?.energy)
    }

    /// Native consistent tangent at the current displacement.
    ///
    /// # Errors
    /// Invalid displacement, material-domain refusal, overflow or cancellation.
    pub fn tangent_matrix(
        &self,
        displacement_m: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Csr, HyperTetError> {
        checkpoint(cx)?;
        self.problem.check_len(displacement_m)?;
        let tangent = self
            .problem
            .tangent(&self.elements, displacement_m, cx)?
            .map_err(|(element, what)| HyperTetError::MaterialRefused { element, what })?;
        for row in 0..tangent.nrows() {
            if row % 4096 == 0 {
                checkpoint(cx)?;
            }
            if tangent.row(row).1.iter().any(|value| !value.is_finite()) {
                return Err(HyperTetError::InvalidInput {
                    what: "assembled consistent tangent is not finite".into(),
                });
            }
        }
        checkpoint(cx)?;
        Ok(tangent)
    }

    /// Assemble new held dead loads using the admitted reference volumes.
    ///
    /// Loads are validated on EVERY call and never retained in this object.
    /// The body-force density is per reference volume, N/m^3; nodal loads are N.
    ///
    /// # Errors
    /// Out-of-range nodal index, non-finite load, overflow or cancellation.
    pub fn external_force_n(
        &self,
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, HyperTetError> {
        checkpoint(cx)?;
        let dofs = 3 * self.problem.nodes_m.len();
        if nodal_forces_n
            .iter()
            .any(|&(dof, value)| dof >= dofs || !value.is_finite())
            || body_force_n_m3.iter().any(|value| !value.is_finite())
        {
            return Err(HyperTetError::InvalidInput {
                what: "finite dead loads on admitted nodal DOFs required".into(),
            });
        }
        let loading = HyperTetProblem {
            nodal_forces_n,
            body_force_n_m3,
            ..self.problem
        };
        let force = loading.external(&self.elements, 1.0);
        if force.iter().any(|value| !value.is_finite()) {
            return Err(HyperTetError::InvalidInput {
                what: "assembled dead force is not finite".into(),
            });
        }
        checkpoint(cx)?;
        Ok(force)
    }
}
