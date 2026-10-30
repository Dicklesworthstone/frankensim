//! Initial-stress linearization of the existing three-dimensional solid operator.
//!
//! On the SAME reference mesh and supports as the accepted small-displacement
//! state, form K_t = K_e + K_g with
//! K_g[(a,i),(b,j)] = delta_ij integral grad(N_a)^T sigma_0 grad(N_b) dV.
//! Tensile stress stiffens and compressive stress softens the relevant motions.
//! K_g and K_t need not be positive definite: negative eigenvalues must not be
//! clamped, absolutized or converted into real vibration frequencies.
//!
//! This is the conservative, small-displacement initial-stress approximation,
//! NOT the full finite-strain tangent. It deliberately keeps reference geometry
//! and mass. It includes neither follower-load stiffness nor geometry-dependent
//! boundary forces, material memory, live mode/state transfer or post-buckling.
//! See Bleyer/Pierre, Computational Mechanics Numerical Tours, “Linear Buckling
//! Analysis of a 3D solid”, equation (32), for the weak-form approximation.
//!
//! Elasticity, mass, source/solution binding, stress recovery and quality checks
//! remain in linear3; inversion uses fs-la. No material or eigensolver is copied.

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_la::batched::{BatchMat, batch_inv};
use fs_sparse::{Coo, Csr};
use crate::linear3::{TetElasticAssembly, TetElasticError, TetLinearElasticProblem,
    TetThermalDisplacementSolution, TetThermalStrainField, TetThermalStressReport,
    TetMaterialField};
use crate::linear3::{TetEigenstrainDisplacementSolution, eigenstrain::TetEigenstrainField};

/// The original consistent-mass pencil with its source-bound stress correction.
/// Matrices and stress are read-only together; this is not a stability certificate.
#[derive(Clone, Debug)]
pub struct PrestressedTetAssembly {
    pencil: TetElasticAssembly,
    elastic_stiffness: Csr,
    geometric_stiffness: Csr,
    stress: TetThermalStressReport,
    identity: ContentHash,
}
impl PrestressedTetAssembly {
    /// Existing (K_t, M) carrier, directly consumable by fs-modal.
    pub const fn pencil(&self) -> &TetElasticAssembly { &self.pencil }
    /// Selected incremental material stiffness K_e, without a stress multiplier.
    pub const fn elastic_stiffness(&self) -> &Csr { &self.elastic_stiffness }
    /// Signed initial-stress contribution K_g, in the same free-DOF order.
    pub const fn geometric_stiffness(&self) -> &Csr { &self.geometric_stiffness }
    /// Actual recovered stress, checked against the original accepted state.
    pub const fn stress(&self) -> &TetThermalStressReport { &self.stress }
    /// Binds the approximation to its original geometry/material/support/solution.
    pub const fn identity(&self) -> ContentHash { self.identity }
}

impl TetLinearElasticProblem<'_> {
    /// Initial-stress pencil from an accepted nonthermal eigenstrain solve.
    /// Uses the same reference-geometry approximation and signed stress kernel
    /// as thermal prestress; the physical source remains explicitly nonthermal.
    pub fn assemble_eigenstrain_tangent(
        &self, solution: &TetEigenstrainDisplacementSolution,
        strains: TetEigenstrainField<'_>, mesh_frame: ContentHash,
        maximum_absolute_mandel_strain: f64, cx: &Cx<'_>,
    ) -> Result<PrestressedTetAssembly, TetElasticError> {
        let stress = self.recover_eigenstrain_stress(solution, strains, mesh_frame,
            maximum_absolute_mandel_strain, cx)?;
        self.assemble_tangent_from_stress(stress, cx)
    }
    /// Assemble a thermal-prestressed vibration/stability pencil.
    ///
    /// The original stress recovery checks the EXACT mesh, material, thermal
    /// field, constraints and accepted solution. Neither a caller-authored stress
    /// array nor a solution from an old temperature can silently enter this path.
    /// Incremental vibration retains these same zero-displacement supports.
    /// All-fixed systems correctly refuse with NoFreeDofs.
    ///
    /// The existing assembly budget bounds allocations. Extra work is linear in
    /// element/DOF count, with at most 64 Jacobian inversions between checkpoints;
    /// sparse sorting inherits the existing assembler's bound. No partial pencil
    /// is returned on material, geometry, arithmetic or cancellation failure.
    ///
    /// This does NOT establish dynamic stability. Feed signed eigenvalues to
    /// the existing spectral owner; a nonpositive value is not an audible mode.
    pub fn assemble_thermal_tangent(
        &self, solution: &TetThermalDisplacementSolution,
        thermal: TetThermalStrainField<'_>, mesh_frame: ContentHash,
        maximum_absolute_mandel_strain: f64, cx: &Cx<'_>,
    ) -> Result<PrestressedTetAssembly, TetElasticError> {
        let stress = self.recover_thermal_stress(solution, thermal, mesh_frame,
            maximum_absolute_mandel_strain, cx)?;
        self.assemble_tangent_from_stress(stress, cx)
    }

    /// The same initial-stress pencil for equilibrium with explicit dead nodal
    /// forces, including distributed loads already integrated by loads3.
    ///
    /// Recover stress against the exact original force vector before assembling
    /// the shared correction. A changed force, material, thermal field, mesh or
    /// support cannot reuse a stale displacement. Dead forces add no follower
    /// stiffness; the original thermal-only entrypoint retains its input contract.
    #[allow(clippy::too_many_arguments)]
    pub fn assemble_thermal_tangent_with_nodal_forces(
        &self, solution: &TetThermalDisplacementSolution,
        thermal: TetThermalStrainField<'_>, nodal_forces_n: &[(usize, f64)],
        mesh_frame: ContentHash, maximum_absolute_mandel_strain: f64, cx: &Cx<'_>,
    ) -> Result<PrestressedTetAssembly, TetElasticError> {
        let stress = self.recover_thermal_stress_with_nodal_forces(solution, thermal,
            nodal_forces_n, mesh_frame, maximum_absolute_mandel_strain, cx)?;
        self.assemble_tangent_from_stress(stress, cx)
    }

    /// Preserve the original recovered stress while using an explicitly
    /// supplied INCREMENTAL material tangent. This is needed when eliminating
    /// a coupled state changes d(sigma)/d(strain), as for local adiabatic
    /// thermoelasticity, without changing the accepted finite stress.
    ///
    /// Recover stress against the exact original solution, thermal field and
    /// dead loads FIRST. The replacement tangent is never used to reconstruct
    /// that stress. Reference geometry, supports and each inertial density
    /// remain unchanged. The caller owns the additional constitutive law and
    /// its source identity; this API does not infer adiabaticity by itself.
    ///
    /// The signed geometric stiffness still comes from the single native
    /// initial-stress kernel. Negative eigenvalues remain signed, and this
    /// frozen conservative pencil is not a transient or stability certificate.
    #[allow(clippy::too_many_arguments)]
    pub fn assemble_thermal_tangent_with_constitutive_materials(
        &self, solution: &TetThermalDisplacementSolution,
        thermal: TetThermalStrainField<'_>, nodal_forces_n: &[(usize, f64)],
        constitutive_materials: TetMaterialField<'_>, mesh_frame: ContentHash,
        maximum_absolute_mandel_strain: f64, cx: &Cx<'_>,
    ) -> Result<PrestressedTetAssembly, TetElasticError> {
        let stress = self.recover_thermal_stress_with_nodal_forces(
            solution, thermal, nodal_forces_n, mesh_frame,
            maximum_absolute_mandel_strain, cx)?;
        self.assemble_constitutive_tangent_from_stress(stress, constitutive_materials, cx)
    }

    // Only native stress-recovery owners call this seam. In particular the
    // transient owner recovers from its accepted motion; arbitrary public
    // stress arrays cannot bypass original geometry/material/history binding.
    pub(crate) fn assemble_constitutive_tangent_from_stress(
        &self, stress: TetThermalStressReport,
        constitutive_materials: TetMaterialField<'_>, cx: &Cx<'_>,
    ) -> Result<PrestressedTetAssembly, TetElasticError> {
        checkpoint(cx)?;
        let count = self.tetrahedra.len();
        for field in [self.materials, constitutive_materials] {
            if let TetMaterialField::PerElement(values) = field {
                if values.len() != count {
                    return Err(TetElasticError::MaterialCountMismatch {
                        expected: count, got: values.len(),
                    });
                }
            }
        }
        let mut hash = DomainHasher::new(
            "org.frankensim.fs-solid.reference-constitutive-initial-stress-pencil.v1");
        hash.update(&(count as u64).to_le_bytes());
        for element in 0..count {
            if element % 64 == 0 { checkpoint(cx)?; }
            let original = match self.materials {
                TetMaterialField::Uniform(value) => value,
                TetMaterialField::PerElement(values) => &values[element],
            };
            let tangent = match constitutive_materials {
                TetMaterialField::Uniform(value) => value,
                TetMaterialField::PerElement(values) => &values[element],
            };
            if original.density_kg_m3.to_bits() != tangent.density_kg_m3.to_bits() {
                return Err(invalid(
                    "an incremental constitutive tangent cannot replace reference inertial density"));
            }
            hash.update(&tangent.material_state_identity);
            hash.update(&tangent.density_kg_m3.to_bits().to_le_bytes());
            for value in tangent.stiffness_mandel_pa().iter().flatten() {
                hash.update(&value.to_bits().to_le_bytes());
            }
        }
        let problem = TetLinearElasticProblem {
            nodes_m: self.nodes_m, tetrahedra: self.tetrahedra,
            materials: constitutive_materials, fixed_dofs: self.fixed_dofs,
            budget: self.budget,
        };
        let mut result = problem.assemble_tangent_from_stress(stress, cx)?;
        // The original identity binds the accepted stress; the new one ALSO
        // binds the selected incremental constitutive law and exact coefficients.
        hash.update(result.identity().as_bytes());
        // Dynamic stress recovery binds geometry and motion independently of
        // supports. Bind the actual reduced pencil's canonical DOF order too.
        hash.update(&(result.pencil.free_dofs.len() as u64).to_le_bytes());
        for (i, dof) in result.pencil.free_dofs.iter().enumerate() {
            if i % 4096 == 0 { checkpoint(cx)?; }
            hash.update(&(*dof as u64).to_le_bytes());
        }
        result.identity = hash.finalize();
        checkpoint(cx)?;
        Ok(result)
    }

    fn assemble_tangent_from_stress(
        &self, stress: TetThermalStressReport, cx: &Cx<'_>,
    ) -> Result<PrestressedTetAssembly, TetElasticError> {
        let mut pencil = self.assemble(cx)?;
        checkpoint(cx)?;
        let mut reduced = vec![None; self.nodes_m.len().checked_mul(3)
            .ok_or(TetElasticError::SizeOverflow)?];
        for (i, &dof) in pencil.free_dofs.iter().enumerate() {
            if i % 4096 == 0 { checkpoint(cx)?; }
            reduced[dof] = Some(i);
        }
        let n = pencil.free_dofs.len();
        let mut geometric = Coo::new(n, n);
        // The original solid owner has already checked every element and its
        // quality. Use the existing batched inverse, not another inverse kernel.
        for (batch, cells) in self.tetrahedra.chunks(64).enumerate() {
            checkpoint(cx)?;
            let jac = BatchMat::from_fn(3, cells.len(), |e, i, j| {
                self.nodes_m[cells[e][j + 1]][i] - self.nodes_m[cells[e][0]][i]
            });
            let mut inverse = BatchMat::zeros(3, cells.len());
            if !batch_inv(&jac, &mut inverse).is_empty() {
                return Err(invalid("admitted initial-stress Jacobian could not be inverted"));
            }
            for (local, tet) in cells.iter().enumerate() {
                let e = batch * 64 + local;
                let mut grad = [[0.0; 3]; 4];
                for a in 1..4 {
                    for i in 0..3 { grad[a][i] = inverse.get(local, a - 1, i); }
                }
                for i in 0..3 { grad[0][i] = -(grad[1][i] + grad[2][i] + grad[3][i]); }
                if grad.iter().flatten().any(|v| !v.is_finite()) {
                    return Err(invalid("initial-stress basis gradient is nonfinite"));
                }
                let s = stress.elements()[e].stress_mandel_pa();
                let h = core::f64::consts::FRAC_1_SQRT_2;
                let sigma = [[s[0], h*s[3], h*s[5]],
                    [h*s[3], s[1], h*s[4]], [h*s[5], h*s[4], s[2]]];
                for a in 0..4 {
                    for b in a..4 {
                        let mut value = 0.0;
                        for i in 0..3 {
                            for j in 0..3 { value += grad[a][i] * sigma[i][j] * grad[b][j]; }
                        }
                        value *= pencil.element_volumes_m3[e];
                        if !value.is_finite() { return Err(invalid("initial-stress stiffness is nonfinite")); }
                        if value == 0.0 { continue; }
                        for axis in 0..3 {
                            if let (Some(row), Some(col)) = (reduced[3*tet[a]+axis], reduced[3*tet[b]+axis]) {
                                geometric.push(row, col, value);
                                // Copy one contraction to its transpose exactly.
                                if a != b { geometric.push(col, row, value); }
                            }
                        }
                    }
                }
            }
        }
        checkpoint(cx)?;
        let geometric_stiffness = geometric.assemble();
        let elastic_stiffness = pencil.stiffness.clone();
        let mut tangent = Coo::new(n, n);
        for row in 0..n {
            if row % 256 == 0 { checkpoint(cx)?; }
            for matrix in [&elastic_stiffness, &geometric_stiffness] {
                let (columns, values) = matrix.row(row);
                for (&column, &value) in columns.iter().zip(values) {
                    if !value.is_finite() { return Err(invalid("initial-stress pencil contains nonfinite stiffness")); }
                    tangent.push(row, column, value);
                }
            }
        }
        pencil.stiffness = tangent.assemble();
        for row in 0..n {
            if row % 256 == 0 { checkpoint(cx)?; }
            if pencil.stiffness.row(row).1.iter().chain(pencil.mass.row(row).1).any(|v| !v.is_finite()) {
                return Err(invalid("initial-stress assembled pencil is unrepresentable"));
            }
        }
        if !pencil.total_mass_kg.is_finite() {
            return Err(invalid("initial-stress mass is unrepresentable"));
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-solid.reference-initial-stress-pencil.v1");
        hash.update(stress.identity().as_bytes());
        let identity = hash.finalize();
        checkpoint(cx)?;
        Ok(PrestressedTetAssembly { pencil, elastic_stiffness, geometric_stiffness, stress, identity })
    }
}
fn invalid(what: &'static str) -> TetElasticError { TetElasticError::InvalidStressObservation { what } }
fn checkpoint(cx: &Cx<'_>) -> Result<(), TetElasticError> {
    cx.checkpoint().map_err(|_| TetElasticError::Cancelled)
}

