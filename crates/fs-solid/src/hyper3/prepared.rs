//! Reuse immutable reference geometry throughout a finite-deformation trajectory.
//!
//! Shape gradients and reference volumes depend on the reference mesh, not on
//! its current displacement. Preparation admits them once; every evaluation
//! still checks the current deformation and the original material domain.
//! Borrowed mesh and law data cannot change while this preparation is alive.
//! Constitutive evaluation and assembly remain owned by the parent hyper3
//! kernel, with the same deterministic element order and cancellation points.

mod energy_momentum;

use fs_exec::Cx;
use fs_material::hyper::Hyperelastic;
use fs_material::hyper::cards::FiniteStrainPotential;
use fs_sparse::Csr;

use super::{Element, Evaluated, HyperTetError, HyperTetProblem, checkpoint};

/// An admitted immutable mesh and material binding for repeated evaluations.
///
/// Retains O(number of tetrahedra) reference gradients and volumes. No tangent,
/// current force, current deformation or constitutive history is cached.
#[derive(Debug)]
pub struct PreparedHyperTetProblem<'a, L: FiniteStrainPotential = Hyperelastic> {
    problem: HyperTetProblem<'a, L>,
    elements: Vec<Element>,
}

impl<'a, L: FiniteStrainPotential> HyperTetProblem<'a, L> {
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

impl<L: FiniteStrainPotential> PreparedHyperTetProblem<'_, L> {
    /// Current TOTAL P1 deformation gradient of one original element, row-major.
    ///
    /// The callback supplies displacement in metres for each of this element's
    /// four ORIGINAL global node IDs, once each in original local-node order.
    /// This sparse observation uses the same cached reference gradients and
    /// translation-invariant arithmetic as native force/energy/tangent assembly;
    /// it requires no full-mesh displacement buffer or reference preparation.
    ///
    /// This admits finite displacement and positive finite det F, not the
    /// constitutive domain or an elastic/thermal deformation substitution.
    /// Original material evaluation remains the force/energy owner's work.
    /// Callback-owned side effects are outside this read-only geometry query.
    ///
    /// # Errors
    /// Missing element, nonfinite displacement, invalid current total F, or
    /// cancellation before/after the bounded four-node observation.
    pub fn deformation_gradient_with(
        &self,
        element: usize,
        mut displacement_at_node: impl FnMut(usize) -> [f64; 3],
        cx: &Cx<'_>,
    ) -> Result<[f64; 9], HyperTetError> {
        checkpoint(cx)?;
        let geometry = self
            .elements
            .get(element)
            .ok_or_else(|| HyperTetError::InvalidInput {
                what: "deformation element is outside the original prepared mesh".into(),
            })?;
        let mut finite = true;
        let f = HyperTetProblem::<L>::deformation_gradient_with(geometry, |node| {
            let displacement = displacement_at_node(node);
            finite &= displacement.iter().all(|value| value.is_finite());
            displacement
        });
        if !finite {
            return Err(HyperTetError::InvalidInput {
                what: "finite nodal displacement required for deformation observation".into(),
            });
        }
        HyperTetProblem::<L>::deformation_jacobian(&f)
            .map_err(|what| HyperTetError::MaterialRefused { element, what })?;
        checkpoint(cx)?;
        Ok(f)
    }

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

    /// Native material potential, J, with full current-state admission.
    ///
    /// Equals stored strain energy for hyperelastic laws. Frozen inelastic
    /// adapters own the distinction between incremental and stored energy.
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

    /// Conservative gravity on current quadrature positions, retaining the
    /// original reference geometry and element mass inventory.
    ///
    /// The density field is the same one supplied to the native reference mass
    /// matrix. No deformed-volume mass update or nodal density average occurs.
    /// The callback supplies the derivatives of one stationary specific potential.
    ///
    /// # Errors
    /// Invalid density/displacement/field, overflow, or cancellation.
    pub fn gravity_load(
        &self,
        reference_density_kg_m3: &[f64],
        displacement_m: &[f64],
        field: &super::gravity::GravityField3<'_>,
        cx: &Cx<'_>,
    ) -> Result<super::gravity::TetGravityLoad3, HyperTetError> {
        self.gravity_load_with_density(
            super::gravity::ReferenceGravityDensity::PerElement(reference_density_kg_m3),
            displacement_m,
            field,
            cx,
        )
    }

    /// Current conservative gravity for positive P1 density on the original
    /// mesh nodes. The cubic rule integrates the graded consistent mass exactly
    /// for affine acceleration and retains reference inventory under deformation.
    /// The callback is sampled at current vertices and current face centroids;
    /// general fields still require spatial refinement.
    ///
    /// # Errors
    /// Invalid density/displacement/field, overflow, or cancellation.
    pub fn gravity_load_with_nodal_density(
        &self,
        reference_nodal_density_kg_m3: &[f64],
        displacement_m: &[f64],
        field: &super::gravity::GravityField3<'_>,
        cx: &Cx<'_>,
    ) -> Result<super::gravity::TetGravityLoad3, HyperTetError> {
        self.gravity_load_with_density(
            super::gravity::ReferenceGravityDensity::Nodal(reference_nodal_density_kg_m3),
            displacement_m,
            field,
            cx,
        )
    }

    fn gravity_load_with_density(
        &self,
        density: super::gravity::ReferenceGravityDensity<'_>,
        displacement_m: &[f64],
        field: &super::gravity::GravityField3<'_>,
        cx: &Cx<'_>,
    ) -> Result<super::gravity::TetGravityLoad3, HyperTetError> {
        let result =
            self.gravity_evaluate_with_density(density, displacement_m, field, true, cx)?;
        Ok(super::gravity::TetGravityLoad3 {
            force_n: result.force_n,
            force_tangent_n_m: result.tangent.expect("requested gravity tangent"),
            potential_energy_j: result.potential_energy_j,
        })
    }

    /// The identical native gravity force and potential without Newton-tangent
    /// allocation or assembly. Intended for force/impulse splitting consumers.
    /// Returns `(full_nodal_force_n, potential_energy_j)`.
    pub fn gravity_force_and_potential(
        &self,
        reference_density_kg_m3: &[f64],
        displacement_m: &[f64],
        field: &super::gravity::GravityField3<'_>,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, f64), HyperTetError> {
        let result = self.gravity_evaluate_with_density(
            super::gravity::ReferenceGravityDensity::PerElement(reference_density_kg_m3),
            displacement_m,
            field,
            false,
            cx,
        )?;
        Ok((result.force_n, result.potential_energy_j))
    }

    /// Force/potential-only counterpart retaining the original positive nodal
    /// density quadrature. No cell-average or lumped-mass substitution is made.
    pub fn gravity_force_and_potential_with_nodal_density(
        &self,
        reference_nodal_density_kg_m3: &[f64],
        displacement_m: &[f64],
        field: &super::gravity::GravityField3<'_>,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, f64), HyperTetError> {
        let result = self.gravity_evaluate_with_density(
            super::gravity::ReferenceGravityDensity::Nodal(reference_nodal_density_kg_m3),
            displacement_m,
            field,
            false,
            cx,
        )?;
        Ok((result.force_n, result.potential_energy_j))
    }

    fn gravity_evaluate_with_density(
        &self,
        density: super::gravity::ReferenceGravityDensity<'_>,
        displacement_m: &[f64],
        field: &super::gravity::GravityField3<'_>,
        tangent: bool,
        cx: &Cx<'_>,
    ) -> Result<super::gravity::GravityEvaluation, HyperTetError> {
        let weights = super::gravity::bind_reference_weights(
            &self.elements,
            self.problem.nodes_m.len(),
            density,
            cx,
        )?;
        super::gravity::evaluate(
            self.problem.nodes_m,
            &self.elements,
            &weights,
            displacement_m,
            field,
            tangent,
            cx,
        )
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
