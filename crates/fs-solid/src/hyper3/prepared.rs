//! Reuse immutable reference geometry throughout a finite-deformation trajectory.
//!
//! Shape gradients and reference volumes depend on the reference mesh, not on
//! its current displacement. Preparation admits them once; every evaluation
//! still checks the current deformation and the original material domain.
//! Borrowed mesh and law data cannot change while this preparation is alive.
//! Constitutive evaluation and assembly remain owned by the parent hyper3
//! kernel, with the same deterministic element order and cancellation points.

mod energy_momentum;
/// Source-bound tensile assessment of original finite-strain Cauchy stresses.
pub mod tensile;

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

/// Original material storage at two configurations and its direct change [J].
///
/// `change_j` is accumulated from each law's own potential increment, not from
/// the difference of the rounded totals. Frozen inelastic adapters retain their
/// original incremental-potential semantics; this does not relabel dissipation
/// as stored elastic energy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HyperTetEnergyInterval {
    pub before_j: f64,
    pub after_j: f64,
    pub change_j: f64,
}

/// Exact material tangent at one fixed displacement and constitutive history.
///
/// Borrows the original prepared reference geometry and retains one 9-by-9
/// constitutive Hessian per P1 cell, without a global sparse matrix or element
/// 12-by-12 matrices. Repeated actions do not reevaluate a material law. Rebuild
/// this operator when displacement or the held constitutive history changes.
/// This is the complete derivative of native internal force; the original
/// Piola Hessian already includes the finite-strain geometric contribution.
/// No support elimination, follower-load tangent, gravity or inertia is included.
///
/// Actions remain explicitly fallible: solver callbacks must propagate refusal
/// and cancellation rather than interpret an incomplete result as a linear solve.
#[derive(Debug)]
pub struct PreparedHyperTetTangent<'a> {
    elements: &'a [Element],
    tangents: Vec<[[f64; 9]; 9]>,
    dofs: usize,
}

impl PreparedHyperTetTangent<'_> {
    /// Full node-major xyz dimension; output and scratch each need this length.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.dofs
    }

    /// Apply the original tangent to a displacement increment, returning force.
    ///
    /// Scratch may change on any call; output changes only on success. Both
    /// buffers are supplied by the caller and reused without allocation.
    ///
    /// # Errors
    /// Incorrect buffer lengths, nonfinite direction, unrepresentable action,
    /// or cancellation. No partial output is published.
    pub fn apply(
        &self,
        increment_m: &[f64],
        force_increment_n: &mut [f64],
        scratch_n: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), HyperTetError> {
        self.apply_impl::<false>(increment_m, force_increment_n, scratch_n, cx)
    }

    /// Apply the actual transpose, retaining every original tangent component.
    ///
    /// The constitutive Hessian of a potential is symmetric analytically, but
    /// this operation transposes its supplied numerical entries explicitly;
    /// it never silently assumes or repairs symmetry for an adjoint client.
    /// Buffer and refusal semantics are identical to [`Self::apply`].
    ///
    /// # Errors
    /// Incorrect buffer lengths, nonfinite direction, unrepresentable action,
    /// or cancellation.
    pub fn apply_transpose(
        &self,
        increment_m: &[f64],
        force_increment_n: &mut [f64],
        scratch_n: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), HyperTetError> {
        self.apply_impl::<true>(increment_m, force_increment_n, scratch_n, cx)
    }

    fn apply_impl<const TRANSPOSE: bool>(
        &self,
        increment_m: &[f64],
        force_increment_n: &mut [f64],
        scratch_n: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), HyperTetError> {
        begin_action(self.dofs, increment_m, force_increment_n, scratch_n, cx)?;
        for (element, tangent) in self.elements.iter().zip(&self.tangents) {
            checkpoint(cx)?;
            element_action::<TRANSPOSE>(element, tangent, increment_m, scratch_n);
        }
        finish_action(force_increment_n, scratch_n, cx)
    }
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

    /// Constant stress, total strain, volume and potential on every original P1 cell.
    ///
    /// Uses the original cached geometry and material field in connectivity
    /// order. Material interfaces remain discontinuous; no nodal averaging or
    /// constitutive tangent is needed. Potential semantics remain with the law.
    ///
    /// # Errors
    /// Malformed displacement, material-domain refusal, unrepresentable fields
    /// or cancellation. No partial field or constitutive update is published.
    pub fn element_responses(
        &self,
        displacement_m: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<super::response::HyperTetElementResponse>, HyperTetError> {
        super::response::recover_all(&self.problem, &self.elements, displacement_m, cx)
    }

    /// Recover one original element using only its four nodal displacements.
    ///
    /// Calls the supplied observation once for each original node, in original
    /// local-node order. This avoids a full-mesh displacement/result buffer.
    /// Material/source admission uses the total F and original law together;
    /// callback-owned side effects are outside this read-only observation.
    ///
    /// # Errors
    /// Missing element, nonfinite displacement, invalid total deformation,
    /// source/material refusal, unrepresentable response or cancellation.
    pub fn element_response_with(
        &self,
        element: usize,
        displacement_at_node: impl FnMut(usize) -> [f64; 3],
        cx: &Cx<'_>,
    ) -> Result<super::response::HyperTetElementResponse, HyperTetError> {
        let deformation = self.deformation_gradient_with(element, displacement_at_node, cx)?;
        let response =
            super::response::recover(&self.problem, element, &self.elements[element], deformation)?;
        checkpoint(cx)?;
        Ok(response)
    }

    pub(super) fn evaluate(
        &self,
        displacement_m: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Evaluated, HyperTetError> {
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

    /// Evaluate both original material potentials and their direct change [J].
    ///
    /// Streams the cached reference cells once in connectivity order, retaining
    /// each heterogeneous law and reference volume. Compensated products and
    /// sums preserve small local work beside large stored-energy baselines.
    /// Both total deformations and material endpoints are admitted before the
    /// source-owned increment is used. No force vector or tangent is allocated.
    /// Custom laws without an increment override retain their own endpoint
    /// difference; no stronger cancellation guarantee is inferred for them.
    ///
    /// # Errors
    /// Malformed displacement, invalid total deformation, original material
    /// refusal, nonfinite stress/potential/increment, overflow or cancellation.
    pub fn stored_energy_interval_j(
        &self,
        before_displacement_m: &[f64],
        after_displacement_m: &[f64],
        cx: &Cx<'_>,
    ) -> Result<HyperTetEnergyInterval, HyperTetError> {
        use fs_math::dd::Dd;

        checkpoint(cx)?;
        self.problem.check_len(before_displacement_m)?;
        self.problem.check_len(after_displacement_m)?;
        let mut total = [Dd::from_f64(0.0); 3];
        for (index, element) in self.elements.iter().enumerate() {
            if index % 4096 == 0 {
                checkpoint(cx)?;
            }
            let before = HyperTetProblem::<L>::deformation_gradient(element, before_displacement_m);
            let after = HyperTetProblem::<L>::deformation_gradient(element, after_displacement_m);
            let refused = |what| HyperTetError::MaterialRefused {
                element: index,
                what,
            };
            HyperTetProblem::<L>::deformation_jacobian(&before).map_err(refused)?;
            HyperTetProblem::<L>::deformation_jacobian(&after).map_err(refused)?;
            let material = self.problem.material.at(index);
            let (stress_before, energy_before) = material
                .stress_and_potential(&before)
                .map_err(|error| refused(error.to_string()))?;
            let (stress_after, energy_after) = material
                .stress_and_potential(&after)
                .map_err(|error| refused(error.to_string()))?;
            if stress_before
                .iter()
                .chain(&stress_after)
                .any(|v| !v.is_finite())
                || !energy_before.is_finite()
                || !energy_after.is_finite()
            {
                return Err(refused(
                    "nonfinite material endpoint stress or potential".into(),
                ));
            }
            let change = material.energy_density_increment(&before, &after);
            if !change.is_finite() {
                return Err(refused("nonfinite material potential increment".into()));
            }
            let volume = Dd::from_f64(element.volume);
            for (sum, value) in total.iter_mut().zip([energy_before, energy_after, change]) {
                *sum = *sum + volume * Dd::from_f64(value);
            }
        }
        let [before_j, after_j, change_j] = total.map(Dd::to_f64);
        if [before_j, after_j, change_j].iter().any(|v| !v.is_finite()) {
            return Err(HyperTetError::InvalidInput {
                what: "material potential interval is not finite".into(),
            });
        }
        checkpoint(cx)?;
        Ok(HyperTetEnergyInterval {
            before_j,
            after_j,
            change_j,
        })
    }

    /// Direct change in the original material potential, with both endpoints admitted.
    ///
    /// # Errors
    /// The same geometry, material, arithmetic and cancellation refusals as
    /// [`Self::stored_energy_interval_j`].
    pub fn stored_energy_change_j(
        &self,
        before_displacement_m: &[f64],
        after_displacement_m: &[f64],
        cx: &Cx<'_>,
    ) -> Result<f64, HyperTetError> {
        Ok(self
            .stored_energy_interval_j(before_displacement_m, after_displacement_m, cx)?
            .change_j)
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

    /// Apply the exact current material tangent without assembling a matrix.
    ///
    /// Streams the original constitutive Hessian once per cell, computes
    /// `delta F = sum_a delta u_a tensor grad N_a`, then scatters
    /// `V (A : delta F) grad N_a`. Relative nodal increments remove rigid
    /// translation before multiplication; no stiffness, bulk law or material
    /// field is averaged. All original DOFs remain, including support rows.
    ///
    /// The caller supplies two full node-major xyz buffers. Scratch may change
    /// after refusal; output changes only on success. No workspace or matrix is
    /// allocated by this action. The original law owns its evaluation costs.
    /// Use [`Self::linearized_tangent`] to retain local Hessians when many
    /// actions use the same displacement and held constitutive history.
    ///
    /// # Errors
    /// Invalid displacement/direction or buffer sizes, material-domain exit,
    /// nonfinite tangent/action, or cancellation. No partial output is published.
    pub fn tangent_apply(
        &self,
        displacement_m: &[f64],
        increment_m: &[f64],
        force_increment_n: &mut [f64],
        scratch_n: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), HyperTetError> {
        checkpoint(cx)?;
        self.problem.check_len(displacement_m)?;
        begin_action(
            displacement_m.len(),
            increment_m,
            force_increment_n,
            scratch_n,
            cx,
        )?;
        for (index, element) in self.elements.iter().enumerate() {
            checkpoint(cx)?;
            let tangent = self
                .problem
                .element_tangent(index, element, displacement_m)
                .map_err(|what| HyperTetError::MaterialRefused {
                    element: index,
                    what,
                })?;
            element_action::<false>(element, &tangent, increment_m, scratch_n);
        }
        finish_action(force_increment_n, scratch_n, cx)
    }

    /// Freeze original point Hessians for repeated matrix-free forward/adjoint
    /// actions at this displacement and held material history.
    ///
    /// `maximum_tangent_bytes` bounds the retained 81-f64-per-element payload;
    /// existing borrowed geometry and caller-owned action buffers are separate.
    /// This is a payload allowance, not an OS memory/RSS bound. No global or
    /// element stiffness matrix is assembled. The original law supplies every
    /// Hessian through the same seam as [`Self::tangent_matrix`].
    ///
    /// # Errors
    /// Invalid displacement, insufficient payload allowance, material-domain
    /// refusal, nonfinite tangent or cancellation. No partial operator is returned.
    pub fn linearized_tangent(
        &self,
        displacement_m: &[f64],
        maximum_tangent_bytes: usize,
        cx: &Cx<'_>,
    ) -> Result<PreparedHyperTetTangent<'_>, HyperTetError> {
        checkpoint(cx)?;
        self.problem.check_len(displacement_m)?;
        let bytes = self
            .elements
            .len()
            .checked_mul(core::mem::size_of::<[[f64; 9]; 9]>())
            .ok_or_else(|| action_error("local tangent payload size overflow"))?;
        if bytes > maximum_tangent_bytes {
            return Err(action_error(
                "local tangent payload exceeds its byte allowance",
            ));
        }
        let mut tangents = Vec::with_capacity(self.elements.len());
        for (index, element) in self.elements.iter().enumerate() {
            checkpoint(cx)?;
            tangents.push(
                self.problem
                    .element_tangent(index, element, displacement_m)
                    .map_err(|what| HyperTetError::MaterialRefused {
                        element: index,
                        what,
                    })?,
            );
        }
        checkpoint(cx)?;
        Ok(PreparedHyperTetTangent {
            elements: &self.elements,
            tangents,
            dofs: displacement_m.len(),
        })
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

    /// Change in the original reference-mass gravitational potential [J].
    ///
    /// The callback supplies `phi(after)-phi(before)` at the same material
    /// quadrature points as the force owner. Analytic field differences cancel
    /// arbitrary potential datums before rounding; subtracting two independently
    /// rounded total energies cannot recover small physical work. Both cellwise
    /// and positive nodal P1 inventories retain their original quadrature.
    /// Explicit endpoint source snapshots may include reciprocal source motion;
    /// their original field parameters and domains belong to the callback owner.
    /// This evaluates storage, not a midpoint-force approximation to work, and
    /// allocates neither a nodal force vector nor a Newton tangent.
    ///
    /// # Errors
    /// Invalid density/displacement, field-domain refusal, overflow or cancellation.
    pub fn gravity_potential_energy_change(
        &self,
        density: super::gravity::ReferenceGravityDensity<'_>,
        before_displacement_m: &[f64],
        after_displacement_m: &[f64],
        field_increment: &super::gravity::GravityPotentialIncrement3<'_>,
        cx: &Cx<'_>,
    ) -> Result<f64, HyperTetError> {
        let weights = self.bind_gravity_weights(density, cx)?;
        super::gravity::potential_increment(
            self.problem.nodes_m,
            &self.elements,
            &weights,
            before_displacement_m,
            after_displacement_m,
            field_increment,
            cx,
        )
    }

    fn gravity_evaluate_with_density(
        &self,
        density: super::gravity::ReferenceGravityDensity<'_>,
        displacement_m: &[f64],
        field: &super::gravity::GravityField3<'_>,
        tangent: bool,
        cx: &Cx<'_>,
    ) -> Result<super::gravity::GravityEvaluation, HyperTetError> {
        let weights = self.bind_gravity_weights(density, cx)?;
        self.gravity_evaluate_with_weights(&weights, displacement_m, field, tangent, cx)
    }

    /// Internal load owners retain these weights with this exact prepared mesh.
    /// Repeated field evaluations must not rebind or average its mass inventory.
    pub(super) fn bind_gravity_weights(
        &self,
        density: super::gravity::ReferenceGravityDensity<'_>,
        cx: &Cx<'_>,
    ) -> Result<super::gravity::ReferenceGravityWeights, HyperTetError> {
        super::gravity::bind_reference_weights(
            &self.elements,
            self.problem.nodes_m.len(),
            density,
            cx,
        )
    }

    pub(super) fn gravity_evaluate_with_weights(
        &self,
        weights: &super::gravity::ReferenceGravityWeights,
        displacement_m: &[f64],
        field: &super::gravity::GravityField3<'_>,
        tangent: bool,
        cx: &Cx<'_>,
    ) -> Result<super::gravity::GravityEvaluation, HyperTetError> {
        super::gravity::evaluate(
            self.problem.nodes_m,
            &self.elements,
            weights,
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

fn action_error(what: &str) -> HyperTetError {
    HyperTetError::InvalidInput {
        what: format!("matrix-free finite-strain tangent: {what}"),
    }
}

fn begin_action(
    dofs: usize,
    increment: &[f64],
    output: &[f64],
    scratch: &mut [f64],
    cx: &Cx<'_>,
) -> Result<(), HyperTetError> {
    checkpoint(cx)?;
    if [increment.len(), output.len(), scratch.len()]
        .iter()
        .any(|&length| length != dofs)
    {
        return Err(action_error(
            "direction, output and scratch must have the full DOF count",
        ));
    }
    for chunk in increment.chunks(4096) {
        checkpoint(cx)?;
        if chunk.iter().any(|value| !value.is_finite()) {
            return Err(action_error("direction must be finite"));
        }
    }
    for chunk in scratch.chunks_mut(4096) {
        checkpoint(cx)?;
        chunk.fill(0.0);
    }
    Ok(())
}

fn finish_action(output: &mut [f64], scratch: &[f64], cx: &Cx<'_>) -> Result<(), HyperTetError> {
    for chunk in scratch.chunks(4096) {
        checkpoint(cx)?;
        if chunk.iter().any(|value| !value.is_finite()) {
            return Err(action_error("force increment is not representable"));
        }
    }
    checkpoint(cx)?;
    output.copy_from_slice(scratch);
    Ok(())
}

fn element_action<const TRANSPOSE: bool>(
    element: &Element,
    tangent: &[[f64; 9]; 9],
    increment: &[f64],
    output: &mut [f64],
) {
    // Start at zero, not deformation_gradient(increment)-I: the latter loses
    // increments smaller than the identity's ULP. Relative motion also makes
    // a constant increment an EXACT null vector before constitutive work.
    let mut delta_f = [0.0; 9];
    for (a, &node) in element.nodes.iter().enumerate().skip(1) {
        for i in 0..3 {
            let delta_u = increment[3 * node + i] - increment[3 * element.nodes[0] + i];
            for j in 0..3 {
                delta_f[3 * i + j] = delta_u.mul_add(element.gradients[a][j], delta_f[3 * i + j]);
            }
        }
    }
    let delta_p: [f64; 9] = core::array::from_fn(|ij| {
        (0..9).fold(0.0, |sum, kl| {
            let a = if TRANSPOSE {
                tangent[kl][ij]
            } else {
                tangent[ij][kl]
            };
            a.mul_add(delta_f[kl], sum)
        })
    });
    for (a, &node) in element.nodes.iter().enumerate() {
        for i in 0..3 {
            let force = (0..3).fold(0.0, |sum, j| {
                delta_p[3 * i + j].mul_add(element.gradients[a][j], sum)
            });
            output[3 * node + i] = element.volume.mul_add(force, output[3 * node + i]);
        }
    }
}
