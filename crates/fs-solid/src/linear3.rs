//! Three-dimensional small-strain elasticity on body-fitted linear
//! tetrahedra.
//!
//! This module owns the generic volume operator needed by solid vibration,
//! transient deformation, and coupled acoustics. It assembles the standard
//! Galerkin stiffness and consistent mass matrices
//!
//! `K = sum_e integral B^T C B dV`, `M = sum_e integral rho N^T N dV`
//!
//! for isotropic material state points. No object class or material name is
//! encoded: a squat metal disc, a wooden body, or a glass support are merely
//! geometry plus resolved material fields. The output is the ordinary
//! `(K,M)` pencil consumed by `fs-modal`.
//!
//! The v1 applicability boundary is intentionally narrow and explicit:
//! infinitesimal strain, isotropic elasticity, undeformed geometry, and
//! conforming non-degenerate P1 tetrahedra. Plasticity, finite strain,
//! thermoelastic prestress, phase change, and evolving topology must select a
//! different constitutive/kinematic rung rather than silently passing through
//! this operator.

use std::collections::BTreeSet;

use fs_exec::Cx;
use fs_material::state_point::IsotropicSolidStatePoint;
use fs_sparse::{Coo, Csr};

/// Bounded assembly request. Limits are checked before matrix allocation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TetAssemblyBudget {
    /// Maximum admitted mesh vertices.
    pub maximum_nodes: usize,
    /// Maximum admitted tetrahedra.
    pub maximum_tetrahedra: usize,
    /// Maximum admitted free displacement degrees of freedom.
    pub maximum_free_dofs: usize,
    /// Minimum `abs(det(J)) / h_max^3` admitted for any tetrahedron.
    pub minimum_scaled_jacobian: f64,
}

impl TetAssemblyBudget {
    /// A conservative interactive/default envelope.
    #[must_use]
    pub const fn standard() -> Self {
        Self {
            maximum_nodes: 1_000_000,
            maximum_tetrahedra: 6_000_000,
            maximum_free_dofs: 3_000_000,
            minimum_scaled_jacobian: 1.0e-12,
        }
    }
}

impl Default for TetAssemblyBudget {
    fn default() -> Self {
        Self::standard()
    }
}

/// One uniform isotropic elastic state carried into the element operator.
#[derive(Debug, Clone, PartialEq)]
pub struct TetElasticMaterial {
    /// Density [kg/m^3].
    pub density_kg_m3: f64,
    /// Young's modulus [Pa].
    pub young_modulus_pa: f64,
    /// Poisson ratio [-].
    pub poisson_ratio: f64,
    /// Raw bytes of the evidence-bearing `fs-material` state identity.
    pub material_state_identity: [u8; 32],
}

impl TetElasticMaterial {
    /// Construct from scalar properties and an upstream material-state
    /// identity. The identity is mandatory because otherwise changing a
    /// temperature-dependent property could leave no trace in downstream
    /// modal or acoustic artifacts.
    pub fn try_new(
        density_kg_m3: f64,
        young_modulus_pa: f64,
        poisson_ratio: f64,
        material_state_identity: [u8; 32],
    ) -> Result<Self, TetElasticError> {
        let material = Self {
            density_kg_m3,
            young_modulus_pa,
            poisson_ratio,
            material_state_identity,
        };
        material.validate()?;
        Ok(material)
    }

    /// Bind directly to the canonical temperature/environment-resolved
    /// isotropic material state.
    #[must_use]
    pub fn from_resolved_state(state: &IsotropicSolidStatePoint) -> Self {
        Self {
            density_kg_m3: state.density_kg_m3(),
            young_modulus_pa: state.young_modulus_pa(),
            poisson_ratio: state.poisson_ratio(),
            material_state_identity: *state.resolved().identity().as_bytes(),
        }
    }

    fn validate(&self) -> Result<(), TetElasticError> {
        if !(self.density_kg_m3.is_finite() && self.density_kg_m3 > 0.0) {
            return Err(TetElasticError::InvalidMaterial {
                what: "density_kg_m3 must be finite and positive",
            });
        }
        if !(self.young_modulus_pa.is_finite() && self.young_modulus_pa > 0.0) {
            return Err(TetElasticError::InvalidMaterial {
                what: "young_modulus_pa must be finite and positive",
            });
        }
        if !(self.poisson_ratio.is_finite()
            && self.poisson_ratio > -1.0
            && self.poisson_ratio < 0.5)
        {
            return Err(TetElasticError::InvalidMaterial {
                what: "poisson_ratio must lie in (-1, 0.5)",
            });
        }
        if self.material_state_identity == [0; 32] {
            return Err(TetElasticError::InvalidMaterial {
                what: "material_state_identity must not be the zero identity",
            });
        }
        Ok(())
    }
}

/// Material field over tetrahedra. A one-card solid uses `Uniform`; a
/// multi-material or spatially varying state uses one entry per element.
#[derive(Debug, Clone, Copy)]
pub enum TetMaterialField<'a> {
    /// One state for every tetrahedron.
    Uniform(&'a TetElasticMaterial),
    /// One state per tetrahedron, in connectivity order.
    PerElement(&'a [TetElasticMaterial]),
}

impl TetMaterialField<'_> {
    fn validate(&self, element_count: usize) -> Result<(), TetElasticError> {
        match self {
            Self::Uniform(material) => material.validate(),
            Self::PerElement(materials) => {
                if materials.len() != element_count {
                    return Err(TetElasticError::MaterialCountMismatch {
                        expected: element_count,
                        got: materials.len(),
                    });
                }
                for material in *materials {
                    material.validate()?;
                }
                Ok(())
            }
        }
    }

    fn at(&self, element: usize) -> &TetElasticMaterial {
        match self {
            Self::Uniform(material) => material,
            Self::PerElement(materials) => &materials[element],
        }
    }
}

/// A body-fitted linear-tetrahedron elasticity/mass assembly request.
pub struct TetLinearElasticProblem<'a> {
    /// Vertex coordinates [m].
    pub nodes_m: &'a [[f64; 3]],
    /// Four vertex indices per conforming tetrahedron.
    pub tetrahedra: &'a [[usize; 4]],
    /// Uniform or per-element resolved material field.
    pub materials: TetMaterialField<'a>,
    /// Strongly fixed global displacement DOFs (`3*node + component`).
    pub fixed_dofs: &'a [usize],
    /// Explicit work/quality envelope.
    pub budget: TetAssemblyBudget,
}

/// Canonical reduced matrices and the map back to full displacement DOFs.
#[derive(Debug, Clone)]
pub struct TetElasticAssembly {
    /// Reduced symmetric stiffness matrix [N/m].
    pub stiffness: Csr,
    /// Reduced symmetric consistent mass matrix [kg].
    pub mass: Csr,
    /// `free_dofs[i]` is the full displacement DOF represented by reduced row i.
    pub free_dofs: Vec<usize>,
    /// Total physical mass `sum rho_e V_e` [kg].
    pub total_mass_kg: f64,
    /// Smallest admitted element quality `abs(det(J))/h_max^3`.
    pub minimum_scaled_jacobian: f64,
    /// Positive element volumes [m^3], in connectivity order.
    pub element_volumes_m3: Vec<f64>,
}

/// Typed, stable refusal surface for 3-D elastic assembly.
#[derive(Debug, Clone, PartialEq)]
pub enum TetElasticError {
    /// Empty node or element set.
    EmptyMesh,
    /// A coordinate is NaN or infinite.
    NonFiniteNode {
        /// Zero-based node index.
        node: usize,
    },
    /// A connectivity entry is outside the node array.
    NodeIndexOutOfRange {
        /// Zero-based tetrahedron index.
        element: usize,
        /// Invalid node index.
        node: usize,
    },
    /// A tetrahedron repeats one of its four vertices.
    RepeatedTetNode {
        /// Zero-based tetrahedron index.
        element: usize,
    },
    /// A mesh node belongs to no tetrahedron.
    UnreferencedNode {
        /// Zero-based node index.
        node: usize,
    },
    /// A tetrahedron is singular or below the declared quality floor.
    DegenerateTet {
        /// Zero-based tetrahedron index.
        element: usize,
        /// Measured dimensionless quality.
        scaled_jacobian: f64,
        /// Declared quality floor.
        minimum: f64,
    },
    /// A material scalar or identity is inadmissible.
    InvalidMaterial {
        /// Failed material invariant.
        what: &'static str,
    },
    /// A per-element material array has the wrong length.
    MaterialCountMismatch {
        /// Required number of cards.
        expected: usize,
        /// Supplied number of cards.
        got: usize,
    },
    /// A fixed displacement DOF is outside the full vector.
    FixedDofOutOfRange {
        /// Invalid full displacement DOF.
        dof: usize,
        /// Total full displacement DOFs.
        full_dofs: usize,
    },
    /// A fixed displacement DOF was declared more than once.
    DuplicateFixedDof {
        /// Repeated full displacement DOF.
        dof: usize,
    },
    /// Every displacement DOF was constrained.
    NoFreeDofs,
    /// A checked integer size operation overflowed.
    SizeOverflow,
    /// A declared work or quality envelope was exceeded.
    BudgetExceeded {
        /// Bounded resource.
        what: &'static str,
        /// Requested count.
        requested: usize,
        /// Admitted count.
        maximum: usize,
    },
    /// Cancellation was observed at an element boundary.
    Cancelled,
}

impl core::fmt::Display for TetElasticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyMesh => write!(f, "FS-SOLID-TET-EMPTY-MESH"),
            Self::NonFiniteNode { node } => {
                write!(f, "FS-SOLID-TET-NONFINITE-NODE: node {node}")
            }
            Self::NodeIndexOutOfRange { element, node } => write!(
                f,
                "FS-SOLID-TET-NODE-RANGE: element {element} names node {node}"
            ),
            Self::RepeatedTetNode { element } => {
                write!(f, "FS-SOLID-TET-REPEATED-NODE: element {element}")
            }
            Self::UnreferencedNode { node } => {
                write!(f, "FS-SOLID-TET-UNREFERENCED-NODE: node {node}")
            }
            Self::DegenerateTet {
                element,
                scaled_jacobian,
                minimum,
            } => write!(
                f,
                "FS-SOLID-TET-DEGENERATE: element {element} scaled Jacobian \
                 {scaled_jacobian:.3e} below {minimum:.3e}"
            ),
            Self::InvalidMaterial { what } => {
                write!(f, "FS-SOLID-TET-MATERIAL: {what}")
            }
            Self::MaterialCountMismatch { expected, got } => write!(
                f,
                "FS-SOLID-TET-MATERIAL-COUNT: expected {expected}, got {got}"
            ),
            Self::FixedDofOutOfRange { dof, full_dofs } => write!(
                f,
                "FS-SOLID-TET-FIXED-DOF-RANGE: {dof} outside 0..{full_dofs}"
            ),
            Self::DuplicateFixedDof { dof } => {
                write!(f, "FS-SOLID-TET-DUPLICATE-FIXED-DOF: {dof}")
            }
            Self::NoFreeDofs => write!(f, "FS-SOLID-TET-NO-FREE-DOFS"),
            Self::SizeOverflow => write!(f, "FS-SOLID-TET-SIZE-OVERFLOW"),
            Self::BudgetExceeded {
                what,
                requested,
                maximum,
            } => write!(
                f,
                "FS-SOLID-TET-BUDGET: {what} requested {requested}, maximum {maximum}"
            ),
            Self::Cancelled => write!(f, "FS-SOLID-TET-CANCELLED"),
        }
    }
}

impl std::error::Error for TetElasticError {}

impl TetLinearElasticProblem<'_> {
    /// Assemble reduced stiffness and consistent mass matrices.
    ///
    /// Cancellation is polled before allocation and at every element
    /// boundary. No partial matrix is published on refusal.
    pub fn assemble(&self, cx: &Cx<'_>) -> Result<TetElasticAssembly, TetElasticError> {
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        self.validate_budget()?;
        self.materials.validate(self.tetrahedra.len())?;

        let full_dofs = self
            .nodes_m
            .len()
            .checked_mul(3)
            .ok_or(TetElasticError::SizeOverflow)?;
        let mut fixed = BTreeSet::new();
        for &dof in self.fixed_dofs {
            if dof >= full_dofs {
                return Err(TetElasticError::FixedDofOutOfRange { dof, full_dofs });
            }
            if !fixed.insert(dof) {
                return Err(TetElasticError::DuplicateFixedDof { dof });
            }
        }
        let free_dofs: Vec<usize> = (0..full_dofs).filter(|dof| !fixed.contains(dof)).collect();
        if free_dofs.is_empty() {
            return Err(TetElasticError::NoFreeDofs);
        }
        if free_dofs.len() > self.budget.maximum_free_dofs {
            return Err(TetElasticError::BudgetExceeded {
                what: "free_dofs",
                requested: free_dofs.len(),
                maximum: self.budget.maximum_free_dofs,
            });
        }
        let mut reduced_of = vec![None; full_dofs];
        for (reduced, &full) in free_dofs.iter().enumerate() {
            reduced_of[full] = Some(reduced);
        }

        let mut incident = vec![false; self.nodes_m.len()];
        let n = free_dofs.len();
        let mut stiffness = Coo::new(n, n);
        let mut mass = Coo::new(n, n);
        let mut total_mass_kg = 0.0;
        let mut minimum_scaled_jacobian = f64::INFINITY;
        let mut element_volumes_m3 = Vec::with_capacity(self.tetrahedra.len());

        for (element, tet) in self.tetrahedra.iter().enumerate() {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            validate_connectivity(element, tet, self.nodes_m.len())?;
            for &node in tet {
                incident[node] = true;
            }
            let points = tet.map(|node| self.nodes_m[node]);
            let geometry =
                TetGeometry::try_new(element, points, self.budget.minimum_scaled_jacobian)?;
            minimum_scaled_jacobian = minimum_scaled_jacobian.min(geometry.scaled_jacobian);
            element_volumes_m3.push(geometry.volume);

            let material = self.materials.at(element);
            let (ke, me) = element_matrices(&geometry, material);
            total_mass_kg += material.density_kg_m3 * geometry.volume;
            for local_row in 0..12 {
                let full_row = 3 * tet[local_row / 3] + local_row % 3;
                let Some(row) = reduced_of[full_row] else {
                    continue;
                };
                for local_col in 0..12 {
                    let full_col = 3 * tet[local_col / 3] + local_col % 3;
                    let Some(col) = reduced_of[full_col] else {
                        continue;
                    };
                    let stiffness_value = ke[local_row][local_col];
                    if stiffness_value != 0.0 {
                        stiffness.push(row, col, stiffness_value);
                    }
                    let mass_value = me[local_row][local_col];
                    if mass_value != 0.0 {
                        mass.push(row, col, mass_value);
                    }
                }
            }
        }
        if let Some(node) = incident.iter().position(|seen| !seen) {
            return Err(TetElasticError::UnreferencedNode { node });
        }
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        Ok(TetElasticAssembly {
            stiffness: stiffness.assemble(),
            mass: mass.assemble(),
            free_dofs,
            total_mass_kg,
            minimum_scaled_jacobian,
            element_volumes_m3,
        })
    }

    fn validate_budget(&self) -> Result<(), TetElasticError> {
        if self.nodes_m.is_empty() || self.tetrahedra.is_empty() {
            return Err(TetElasticError::EmptyMesh);
        }
        if self.nodes_m.len() > self.budget.maximum_nodes {
            return Err(TetElasticError::BudgetExceeded {
                what: "nodes",
                requested: self.nodes_m.len(),
                maximum: self.budget.maximum_nodes,
            });
        }
        if self.tetrahedra.len() > self.budget.maximum_tetrahedra {
            return Err(TetElasticError::BudgetExceeded {
                what: "tetrahedra",
                requested: self.tetrahedra.len(),
                maximum: self.budget.maximum_tetrahedra,
            });
        }
        if !(self.budget.minimum_scaled_jacobian.is_finite()
            && self.budget.minimum_scaled_jacobian > 0.0)
        {
            return Err(TetElasticError::DegenerateTet {
                element: 0,
                scaled_jacobian: f64::NAN,
                minimum: self.budget.minimum_scaled_jacobian,
            });
        }
        for (node, point) in self.nodes_m.iter().enumerate() {
            if point.iter().any(|value| !value.is_finite()) {
                return Err(TetElasticError::NonFiniteNode { node });
            }
        }
        Ok(())
    }
}

fn validate_connectivity(
    element: usize,
    tet: &[usize; 4],
    node_count: usize,
) -> Result<(), TetElasticError> {
    for &node in tet {
        if node >= node_count {
            return Err(TetElasticError::NodeIndexOutOfRange { element, node });
        }
    }
    let unique: BTreeSet<usize> = tet.iter().copied().collect();
    if unique.len() != 4 {
        return Err(TetElasticError::RepeatedTetNode { element });
    }
    Ok(())
}

struct TetGeometry {
    gradients: [[f64; 3]; 4],
    volume: f64,
    scaled_jacobian: f64,
}

impl TetGeometry {
    fn try_new(element: usize, p: [[f64; 3]; 4], minimum: f64) -> Result<Self, TetElasticError> {
        let j = [sub(p[1], p[0]), sub(p[2], p[0]), sub(p[3], p[0])];
        // `j` stores columns; spell the determinant and inverse in coordinate
        // rows to keep the gradient convention auditable.
        let a = [
            [j[0][0], j[1][0], j[2][0]],
            [j[0][1], j[1][1], j[2][1]],
            [j[0][2], j[1][2], j[2][2]],
        ];
        let det = determinant3(a);
        let h_max = maximum_edge_length(p);
        let scaled_jacobian = det.abs() / (h_max * h_max * h_max);
        if !scaled_jacobian.is_finite() || scaled_jacobian < minimum {
            return Err(TetElasticError::DegenerateTet {
                element,
                scaled_jacobian,
                minimum,
            });
        }
        let inv = inverse3(a, det);
        // grad_x N_i = J^{-T} grad_ref N_i. Rows of J^{-1} are the
        // physical gradients of reference coordinates xi_i.
        let gradients = [
            [
                -inv[0][0] - inv[1][0] - inv[2][0],
                -inv[0][1] - inv[1][1] - inv[2][1],
                -inv[0][2] - inv[1][2] - inv[2][2],
            ],
            inv[0],
            inv[1],
            inv[2],
        ];
        Ok(Self {
            gradients,
            volume: det.abs() / 6.0,
            scaled_jacobian,
        })
    }
}

fn element_matrices(
    geometry: &TetGeometry,
    material: &TetElasticMaterial,
) -> ([[f64; 12]; 12], [[f64; 12]; 12]) {
    let e = material.young_modulus_pa;
    let nu = material.poisson_ratio;
    let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e / (2.0 * (1.0 + nu));
    // Engineering-shear Voigt order: xx, yy, zz, xy, yz, zx.
    let mut d = [[0.0; 6]; 6];
    for (row, drow) in d.iter_mut().enumerate().take(3) {
        for slot in drow.iter_mut().take(3) {
            *slot = lambda;
        }
        drow[row] += 2.0 * mu;
    }
    d[3][3] = mu;
    d[4][4] = mu;
    d[5][5] = mu;

    let mut b = [[0.0; 12]; 6];
    for (node, g) in geometry.gradients.iter().enumerate() {
        let c = 3 * node;
        b[0][c] = g[0];
        b[1][c + 1] = g[1];
        b[2][c + 2] = g[2];
        b[3][c] = g[1];
        b[3][c + 1] = g[0];
        b[4][c + 1] = g[2];
        b[4][c + 2] = g[1];
        b[5][c] = g[2];
        b[5][c + 2] = g[0];
    }
    let mut db = [[0.0; 12]; 6];
    for i in 0..6 {
        for j in 0..12 {
            for k in 0..6 {
                db[i][j] = d[i][k].mul_add(b[k][j], db[i][j]);
            }
        }
    }
    let mut ke = [[0.0; 12]; 12];
    for i in 0..12 {
        for j in 0..12 {
            for k in 0..6 {
                ke[i][j] = b[k][i].mul_add(db[k][j], ke[i][j]);
            }
            ke[i][j] *= geometry.volume;
        }
    }

    let mut me = [[0.0; 12]; 12];
    let mass_scale = material.density_kg_m3 * geometry.volume / 20.0;
    for a in 0..4 {
        for b_node in 0..4 {
            let value = mass_scale * if a == b_node { 2.0 } else { 1.0 };
            for component in 0..3 {
                me[3 * a + component][3 * b_node + component] = value;
            }
        }
    }
    (ke, me)
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn determinant3(a: [[f64; 3]; 3]) -> f64 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}

fn inverse3(a: [[f64; 3]; 3], det: f64) -> [[f64; 3]; 3] {
    let inv_det = det.recip();
    [
        [
            (a[1][1] * a[2][2] - a[1][2] * a[2][1]) * inv_det,
            (a[0][2] * a[2][1] - a[0][1] * a[2][2]) * inv_det,
            (a[0][1] * a[1][2] - a[0][2] * a[1][1]) * inv_det,
        ],
        [
            (a[1][2] * a[2][0] - a[1][0] * a[2][2]) * inv_det,
            (a[0][0] * a[2][2] - a[0][2] * a[2][0]) * inv_det,
            (a[0][2] * a[1][0] - a[0][0] * a[1][2]) * inv_det,
        ],
        [
            (a[1][0] * a[2][1] - a[1][1] * a[2][0]) * inv_det,
            (a[0][1] * a[2][0] - a[0][0] * a[2][1]) * inv_det,
            (a[0][0] * a[1][1] - a[0][1] * a[1][0]) * inv_det,
        ],
    ]
}

fn maximum_edge_length(p: [[f64; 3]; 4]) -> f64 {
    let mut max_squared: f64 = 0.0;
    for i in 0..4 {
        for j in i + 1..4 {
            let d = sub(p[i], p[j]);
            max_squared = max_squared.max(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        }
    }
    max_squared.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    use fs_mesh::{RoundedCylinderMeshSpec, rounded_cylinder_tet_mesh};
    use fs_modal::eigh_gen_dense;

    fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new();
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 1,
                    kernel_id: 3,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    fn material(e: f64, rho: f64) -> TetElasticMaterial {
        TetElasticMaterial::try_new(rho, e, 0.25, [7; 32]).unwrap()
    }

    fn reference_problem<'a>(
        nodes: &'a [[f64; 3]],
        tetrahedra: &'a [[usize; 4]],
        mat: &'a TetElasticMaterial,
        fixed: &'a [usize],
    ) -> TetLinearElasticProblem<'a> {
        TetLinearElasticProblem {
            nodes_m: nodes,
            tetrahedra,
            materials: TetMaterialField::Uniform(mat),
            fixed_dofs: fixed,
            budget: TetAssemblyBudget::standard(),
        }
    }

    #[test]
    fn g0_single_tet_has_rigid_translation_nullspace_and_exact_mass() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tets = [[0, 1, 2, 3]];
        let mat = material(120.0, 6.0);
        let assembly =
            with_cx(|cx| reference_problem(&nodes, &tets, &mat, &[]).assemble(cx)).unwrap();
        assert_eq!(assembly.total_mass_kg, 1.0);
        let k = assembly.stiffness.to_dense();
        let m = assembly.mass.to_dense();
        for component in 0..3 {
            let u: Vec<f64> = (0..12).map(|dof| f64::from(dof % 3 == component)).collect();
            for row in 0..12 {
                let ku: f64 = (0..12).map(|col| k[row * 12 + col] * u[col]).sum();
                assert!(ku.abs() < 1.0e-12, "row {row}: {ku}");
            }
            let directional_mass: f64 = (0..12)
                .map(|row| {
                    (0..12)
                        .map(|col| u[row] * m[row * 12 + col] * u[col])
                        .sum::<f64>()
                })
                .sum();
            assert!((directional_mass - 1.0).abs() < 1.0e-14);
        }
    }

    #[test]
    fn g0_affine_strain_energy_and_rigid_rotation_are_exact() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tets = [[0, 1, 2, 3]];
        let mat = material(120.0, 6.0);
        let assembly =
            with_cx(|cx| reference_problem(&nodes, &tets, &mat, &[]).assemble(cx)).unwrap();
        let k = assembly.stiffness.to_dense();

        // Infinitesimal rotation about z: u=(-y,x,0) has zero symmetric strain.
        let rotation: Vec<f64> = nodes.iter().flat_map(|p| [-p[1], p[0], 0.0]).collect();
        let rotation_energy = quadratic(&k, &rotation);
        assert!(rotation_energy.abs() < 1.0e-12, "{rotation_energy}");

        // u=(alpha*x,0,0): eps_xx=alpha. Analytic energy is
        // 1/2 * (lambda+2mu) * alpha^2 * V.
        let alpha = 0.125;
        let extension: Vec<f64> = nodes
            .iter()
            .flat_map(|p| [alpha * p[0], 0.0, 0.0])
            .collect();
        let nu = mat.poisson_ratio;
        let lambda = mat.young_modulus_pa * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let mu = mat.young_modulus_pa / (2.0 * (1.0 + nu));
        let expected = 0.5 * (lambda + 2.0 * mu) * alpha * alpha / 6.0;
        let actual = 0.5 * quadratic(&k, &extension);
        assert!(
            (actual - expected).abs() < 1.0e-13,
            "{actual} != {expected}"
        );
    }

    #[test]
    fn g1_frequency_scaling_follows_geometry_modulus_and_density() {
        let base = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tets = [[0, 1, 2, 3]];
        // Minimal exact rigid-body suppression: node 0 xyz, node 1 yz,
        // node 2 z. Six elastic coordinates remain.
        let fixed = [0, 1, 2, 4, 5, 8];
        let lowest = |nodes: &[[f64; 3]], e: f64, rho: f64| {
            let mat = material(e, rho);
            let a =
                with_cx(|cx| reference_problem(nodes, &tets, &mat, &fixed).assemble(cx)).unwrap();
            let modes = eigh_gen_dense(
                &a.stiffness.to_dense(),
                &a.mass.to_dense(),
                a.free_dofs.len(),
            )
            .unwrap();
            modes[0].lambda
        };
        let lambda = lowest(&base, 120.0, 6.0);
        let scaled = base.map(|p| p.map(|x| 2.0 * x));
        let lambda_size = lowest(&scaled, 120.0, 6.0);
        let lambda_e = lowest(&base, 480.0, 6.0);
        let lambda_rho = lowest(&base, 120.0, 24.0);
        assert!((lambda_size / lambda - 0.25).abs() < 1.0e-10);
        assert!((lambda_e / lambda - 4.0).abs() < 1.0e-10);
        assert!((lambda_rho / lambda - 0.25).abs() < 1.0e-10);
    }

    #[test]
    fn malformed_topology_materials_and_quality_refuse() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0e-16, 1.0e-16, 0.0],
        ];
        let tets = [[0, 1, 2, 3]];
        let mat = material(1.0, 1.0);
        let error =
            with_cx(|cx| reference_problem(&nodes, &tets, &mat, &[]).assemble(cx)).unwrap_err();
        assert!(matches!(error, TetElasticError::DegenerateTet { .. }));

        let bad = TetElasticMaterial::try_new(1.0, 1.0, 0.25, [0; 32]).unwrap_err();
        assert_eq!(
            bad,
            TetElasticError::InvalidMaterial {
                what: "material_state_identity must not be the zero identity"
            }
        );
    }

    #[test]
    fn g1_rounded_cylinder_geometry_composes_with_elastic_mass() {
        let mesh = with_cx(|cx| {
            rounded_cylinder_tet_mesh(
                RoundedCylinderMeshSpec {
                    outer_radius_m: 0.038,
                    thickness_m: 0.006,
                    fillet_radius_m: 0.001,
                    core_radial_segments: 3,
                    fillet_radial_segments: 2,
                    azimuthal_segments: 12,
                    axial_segments: 2,
                    maximum_vertices: 10_000,
                    maximum_tetrahedra: 50_000,
                },
                cx,
            )
        })
        .unwrap();
        let density = 7_800.0;
        let mat = material(200.0e9, density);
        let assembly = with_cx(|cx| {
            TetLinearElasticProblem {
                nodes_m: &mesh.nodes_m,
                tetrahedra: &mesh.tetrahedra,
                materials: TetMaterialField::Uniform(&mat),
                fixed_dofs: &[],
                budget: TetAssemblyBudget::standard(),
            }
            .assemble(cx)
        })
        .unwrap();
        let geometric_volume: f64 = mesh
            .tetrahedra
            .iter()
            .map(|tet| {
                let [a, b, c, d] = tet.map(|vertex| mesh.nodes_m[vertex]);
                determinant_columns([sub(b, a), sub(c, a), sub(d, a)]).abs() / 6.0
            })
            .sum();
        assert!((assembly.total_mass_kg - density * geometric_volume).abs() < 1.0e-12);
        assert_eq!(assembly.free_dofs.len(), 3 * mesh.nodes_m.len());
        assert!(!mesh.boundary.triangles.is_empty());
    }

    fn quadratic(matrix: &[f64], vector: &[f64]) -> f64 {
        let n = vector.len();
        (0..n)
            .map(|row| {
                (0..n)
                    .map(|col| vector[row] * matrix[row * n + col] * vector[col])
                    .sum::<f64>()
            })
            .sum()
    }

    fn determinant_columns(columns: [[f64; 3]; 3]) -> f64 {
        columns[0][0] * (columns[1][1] * columns[2][2] - columns[1][2] * columns[2][1])
            - columns[1][0] * (columns[0][1] * columns[2][2] - columns[0][2] * columns[2][1])
            + columns[2][0] * (columns[0][1] * columns[1][2] - columns[0][2] * columns[1][1])
    }
}
