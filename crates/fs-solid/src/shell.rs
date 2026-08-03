//! Deterministic, spatial-only flat-plate surrogate assembly.
//!
//! This fixture-scale operator assembles lumped inertia, edge-spring in-plane
//! resistance, and normal-slope/rotation bending penalties on one coplanar,
//! consistently oriented triangular patch. It is estimate-only: it is neither
//! a CST plane-stress membrane nor a Kirchhoff--Love/isotropic continuum
//! discretisation. It is the spatial base used by the Euler-disc ladder; it
//! does not integrate in time or claim curved, multi-patch IGA shell behavior.

use core::marker::PhantomData;

/// A point of a plate mid-surface, in the declared Cartesian SI frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShellNode {
    /// Cartesian position in metres.
    pub position_m: [f64; 3],
}

/// Isotropic material data for the thin-plate spatial law.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShellMaterial {
    /// Young's modulus in Pa.
    pub youngs_modulus_pa: f64,
    /// Poisson ratio, strictly between -1 and 0.5.
    pub poisson_ratio: f64,
    /// Volumetric density in kg / m³.
    pub density_kg_m3: f64,
}

/// Required identity and interpretation fields for an assembly request.
///
/// The present operator always uses SI metres, kilograms, and seconds in a
/// right-handed Cartesian frame.  The IDs let a caller bind the returned
/// spatial operator to its model/source/state ledger records without this
/// synchronous leaf claiming ledger authority itself.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ShellIdentity {
    /// Caller-owned stable model identity.
    pub model_id: String,
    /// Caller-owned stable geometry/source identity.
    pub source_id: String,
    /// Caller-owned state identity; spatial assembly does not advance it.
    pub state_id: String,
}

/// Explicit bounds for this dense, fixture-scale spatial assembler.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AssemblyBudget {
    /// Largest admissible node count.
    pub max_nodes: usize,
    /// Largest admissible triangle count.
    pub max_triangles: usize,
    /// Largest admissible dense matrix entry count per operator.
    pub max_matrix_entries: usize,
    /// Maximum triangle visits; bounds assembly work before allocation.
    pub max_work_units: usize,
    /// Largest operator dimension for exact Jacobi conditioning diagnostics.
    pub max_conditioning_dofs: usize,
}

impl Default for AssemblyBudget {
    fn default() -> Self {
        Self {
            max_nodes: 64,
            max_triangles: 128,
            max_matrix_entries: 64 * 6 * 64 * 6,
            max_work_units: 128,
            max_conditioning_dofs: 48,
        }
    }
}

/// Optional Rayleigh damping, applicable only after mass and stiffness exist.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DampingModel {
    /// No damping operator is assembled.
    None,
    /// `C = αM + βK`, with both coefficients finite and non-negative.
    Rayleigh {
        /// Mass-proportional coefficient in 1/s.
        mass_proportional_per_s: f64,
        /// Stiffness-proportional coefficient in s.
        stiffness_proportional_s: f64,
    },
}

/// Three pointwise pinned supports on the plate's declared normal side.
///
/// Each support removes its node's three translations while retaining all
/// rotations. The normal must point to the same side as the flat plate normal;
/// it intentionally does not invent rotational clamps or tangential friction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShellSupport {
    /// Exactly three distinct node indices.
    pub node_indices: [usize; 3],
    /// Unit normal in the declared Cartesian frame.
    pub normal: [f64; 3],
}

/// A flat, consistently oriented triangular estimate-only plate surrogate.
#[derive(Debug, Clone, PartialEq)]
pub struct ShellPlate {
    /// Mid-surface nodes in deterministic DOF order.
    pub nodes: Vec<ShellNode>,
    /// Oriented triangle connectivity.  Input order is retained as the
    /// deterministic element tie-breaker.
    pub triangles: Vec<[usize; 3]>,
    /// Uniform thickness in metres.
    pub thickness_m: f64,
    /// Isotropic material card.
    pub material: ShellMaterial,
    /// Identity and SI/frame declaration.
    pub identity: ShellIdentity,
    /// Optional three-point pinned support.
    pub support: Option<ShellSupport>,
    /// Damping request, if applicable.
    pub damping: DampingModel,
    /// Work, memory, and conditioning bounds.
    pub budget: AssemblyBudget,
}

/// Marker for a typed mass matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mass {}
/// Marker for a typed stiffness matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stiffness {}
/// Marker for a typed damping matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Damping {}

/// Dense symmetric operator whose physical kind is carried in its type.
#[derive(Debug, Clone, PartialEq)]
pub struct ShellMatrix<K> {
    dimension: usize,
    values: Vec<f64>,
    kind: PhantomData<K>,
}

impl<K> ShellMatrix<K> {
    /// Matrix order.
    #[must_use]
    pub const fn dimension(&self) -> usize {
        self.dimension
    }

    /// Row-major coefficients, including both symmetric triangles.
    #[must_use]
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    /// Apply the operator to a compatible generalized-displacement vector.
    #[must_use]
    pub fn apply(&self, vector: &[f64]) -> Vec<f64> {
        assert_eq!(
            vector.len(),
            self.dimension,
            "operator/vector size mismatch"
        );
        let mut output = vec![0.0; self.dimension];
        for (row, output_value) in output.iter_mut().enumerate() {
            let begin = row * self.dimension;
            *output_value = self.values[begin..begin + self.dimension]
                .iter()
                .zip(vector)
                .map(|(a, x)| a * x)
                .sum();
        }
        output
    }

    /// Quadratic energy/work form `1/2 xᵀAx`.
    #[must_use]
    pub fn quadratic_energy(&self, vector: &[f64]) -> f64 {
        0.5 * vector
            .iter()
            .zip(self.apply(vector))
            .map(|(x, ax)| x * ax)
            .sum::<f64>()
    }
}

/// Raw algebraic diagnostics for a small mixed-unit operator, or a bounded refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum OperatorDiagnostics {
    /// Jacobi algebraic diagnostics were completed within the declared bound.
    Computed {
        /// Smallest raw mass-array eigenvalue; a positive value is a structural
        /// check for the assembled lumped array, not a physical eigenvalue.
        raw_mass_min_eigenvalue: f64,
        /// Raw stiffness-array nullity under the deterministic tolerance.
        raw_stiffness_nullity: usize,
        /// `λ_max / λ_min_positive` of the raw mixed-unit stiffness array.
        /// It is an algebraic spread only, not a dimensionless physical
        /// condition number or continuum-quality certificate.
        raw_stiffness_eigenvalue_spread: f64,
        /// Largest absolute antisymmetry residue before symmetrization.
        symmetry_residual: f64,
    },
    /// Assembly remains valid, but raw algebraic diagnostics exceeded its budget.
    NotComputed {
        /// Why the declared conditioning budget intentionally skipped Jacobi.
        reason: String,
    },
}

/// Completed spatial assembly. `free_dofs` maps reduced rows to the full,
/// deterministic node-major order: x, y, z, rx, ry, rz for every node.
#[derive(Debug, Clone, PartialEq)]
pub struct ShellAssembly {
    /// Full, unconstrained typed mass operator.
    pub full_mass: ShellMatrix<Mass>,
    /// Full, unconstrained typed stiffness operator.
    pub full_stiffness: ShellMatrix<Stiffness>,
    /// Full Rayleigh damping operator when requested.
    pub full_damping: Option<ShellMatrix<Damping>>,
    /// Reduced free-DOF mass operator.
    pub mass: ShellMatrix<Mass>,
    /// Reduced free-DOF stiffness operator.
    pub stiffness: ShellMatrix<Stiffness>,
    /// Reduced damping operator when requested.
    pub damping: Option<ShellMatrix<Damping>>,
    /// Full DOF indices retained in the reduced operators.
    pub free_dofs: Vec<usize>,
    /// Raw algebraic report for the full mixed-unit spatial operator.
    pub diagnostics: OperatorDiagnostics,
}

/// Refusals are structured so callers never receive partial operators.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ShellError {
    /// A field is non-finite, empty, or violates a material/geometry domain.
    InvalidInput {
        /// Stable explanation of the invalid input.
        what: String,
    },
    /// The requested geometry is outside the flat single-patch plate law.
    UnsupportedGeometry {
        /// Stable explanation of the unsupported geometry class.
        what: String,
    },
    /// The requested boundary is malformed or outside normal point supports.
    UnsupportedBoundary {
        /// Stable explanation of the unsupported boundary condition.
        what: String,
    },
    /// Explicit work, memory, or diagnostic input limit was exceeded.
    BudgetExceeded {
        /// Stable explanation of the exceeded explicit budget.
        what: String,
    },
}

impl core::fmt::Display for ShellError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput { what } => write!(f, "invalid shell input: {what}"),
            Self::UnsupportedGeometry { what } => write!(f, "unsupported shell geometry: {what}"),
            Self::UnsupportedBoundary { what } => write!(f, "unsupported shell boundary: {what}"),
            Self::BudgetExceeded { what } => write!(f, "shell assembly budget exceeded: {what}"),
        }
    }
}

impl std::error::Error for ShellError {}

impl ShellPlate {
    /// Assemble complete and optionally three-point-supported spatial operators.
    ///
    /// The routine validates all input and budgets before allocation, then
    /// assembles private matrices and returns only a complete result. It is a
    /// synchronous spatial leaf: no time advancement or partial authority is
    /// exposed.
    pub fn assemble(&self) -> Result<ShellAssembly, ShellError> {
        self.validate()?;
        let dofs = self.nodes.len() * 6;
        let entries = dofs * dofs;
        let mut mass = vec![0.0; entries];
        let mut stiffness = vec![0.0; entries];
        let reference_normal = triangle_geometry(self, self.triangles[0])?.2;

        for &triangle in &self.triangles {
            let (area, edges, normal) = triangle_geometry(self, triangle)?;
            if dot(normal, reference_normal) < 1.0 - 1e-10 {
                return Err(ShellError::UnsupportedGeometry {
                    what: "curved, folded, or inconsistently oriented multi-patch triangles require the fs-iga shell frontend".into(),
                });
            }
            self.assemble_triangle(triangle, area, edges, normal, &mut mass, &mut stiffness);
        }

        let full_mass = matrix::<Mass>(dofs, mass);
        let full_stiffness = matrix::<Stiffness>(dofs, stiffness);
        validate_derived_matrix("mass", full_mass.values())?;
        validate_derived_matrix("stiffness", full_stiffness.values())?;
        let full_damping = match self.damping {
            DampingModel::None => None,
            DampingModel::Rayleigh {
                mass_proportional_per_s,
                stiffness_proportional_s,
            } => Some(matrix::<Damping>(
                dofs,
                full_mass
                    .values
                    .iter()
                    .zip(&full_stiffness.values)
                    .map(|(m, k)| mass_proportional_per_s * m + stiffness_proportional_s * k)
                    .collect(),
            )),
        };
        if let Some(damping) = &full_damping {
            validate_derived_matrix("damping", damping.values())?;
        }
        let diagnostics = diagnostics(
            &full_mass,
            &full_stiffness,
            self.budget.max_conditioning_dofs,
        );
        let free_dofs = self.free_dofs(dofs)?;
        let mass = reduce::<Mass>(&full_mass, &free_dofs);
        let stiffness = reduce::<Stiffness>(&full_stiffness, &free_dofs);
        let damping = full_damping
            .as_ref()
            .map(|operator| reduce::<Damping>(operator, &free_dofs));
        Ok(ShellAssembly {
            full_mass,
            full_stiffness,
            full_damping,
            mass,
            stiffness,
            damping,
            free_dofs,
            diagnostics,
        })
    }

    fn validate(&self) -> Result<(), ShellError> {
        if self.nodes.len() < 3 || self.triangles.is_empty() {
            return Err(ShellError::InvalidInput {
                what: "at least three nodes and one triangle are required".into(),
            });
        }
        if self.nodes.len() > self.budget.max_nodes
            || self.triangles.len() > self.budget.max_triangles
        {
            return Err(ShellError::BudgetExceeded {
                what: "node or triangle limit".into(),
            });
        }
        let dofs = self
            .nodes
            .len()
            .checked_mul(6)
            .ok_or_else(|| ShellError::BudgetExceeded {
                what: "DOF count overflow".into(),
            })?;
        let entries = dofs
            .checked_mul(dofs)
            .ok_or_else(|| ShellError::BudgetExceeded {
                what: "matrix size overflow".into(),
            })?;
        if entries > self.budget.max_matrix_entries
            || self.triangles.len() > self.budget.max_work_units
        {
            return Err(ShellError::BudgetExceeded {
                what: "matrix memory or element-work limit".into(),
            });
        }
        if !self.thickness_m.is_finite() || self.thickness_m <= 0.0 {
            return Err(ShellError::InvalidInput {
                what: "thickness_m must be finite and positive".into(),
            });
        }
        let material = self.material;
        if !material.youngs_modulus_pa.is_finite()
            || material.youngs_modulus_pa <= 0.0
            || !material.density_kg_m3.is_finite()
            || material.density_kg_m3 <= 0.0
            || !material.poisson_ratio.is_finite()
            || !(-1.0 < material.poisson_ratio && material.poisson_ratio < 0.5)
        {
            return Err(ShellError::InvalidInput {
                what: "material requires E > 0, density > 0, and -1 < nu < 0.5".into(),
            });
        }
        if self.identity.model_id.is_empty()
            || self.identity.source_id.is_empty()
            || self.identity.state_id.is_empty()
        {
            return Err(ShellError::InvalidInput {
                what: "model_id, source_id, and state_id are required".into(),
            });
        }
        if self
            .nodes
            .iter()
            .any(|node| node.position_m.iter().any(|value| !value.is_finite()))
        {
            return Err(ShellError::InvalidInput {
                what: "node position must be finite".into(),
            });
        }
        for triangle in &self.triangles {
            if triangle.iter().any(|&node| node >= self.nodes.len())
                || triangle[0] == triangle[1]
                || triangle[1] == triangle[2]
                || triangle[0] == triangle[2]
            {
                return Err(ShellError::InvalidInput {
                    what: "triangle indices must be distinct in range".into(),
                });
            }
        }
        for (index, triangle) in self.triangles.iter().enumerate() {
            let mut canonical = *triangle;
            canonical.sort_unstable();
            if self.triangles[..index].iter().any(|other| {
                let mut other_canonical = *other;
                other_canonical.sort_unstable();
                other_canonical == canonical
            }) {
                return Err(ShellError::InvalidInput {
                    what: "duplicate triangles are not admitted by the flat plate surrogate".into(),
                });
            }
        }
        let mut connected = vec![false; self.triangles.len()];
        connected[0] = true;
        let mut changed = true;
        while changed {
            changed = false;
            for (index, triangle) in self.triangles.iter().enumerate() {
                if connected[index] {
                    continue;
                }
                if self
                    .triangles
                    .iter()
                    .enumerate()
                    .any(|(other_index, other)| {
                        connected[other_index] && triangles_share_edge(*triangle, *other)
                    })
                {
                    connected[index] = true;
                    changed = true;
                }
            }
        }
        if connected.iter().any(|connected| !connected) {
            return Err(ShellError::UnsupportedGeometry {
                what: "disconnected triangles require separate plate assemblies".into(),
            });
        }
        let mut used_nodes = vec![false; self.nodes.len()];
        for triangle in &self.triangles {
            for &node in triangle {
                used_nodes[node] = true;
            }
        }
        if used_nodes.iter().any(|used| !used) {
            return Err(ShellError::InvalidInput {
                what: "every plate node must be used by a triangle".into(),
            });
        }
        let (_, _, reference_normal) = triangle_geometry(self, self.triangles[0])?;
        let reference_point = self.nodes[self.triangles[0][0]].position_m;
        let geometric_scale = self
            .nodes
            .iter()
            .map(|node| norm(sub(node.position_m, reference_point)))
            .fold(1.0_f64, f64::max);
        let coplanarity_tolerance = 1.0e-10 * geometric_scale;
        if self.nodes.iter().any(|node| {
            dot(sub(node.position_m, reference_point), reference_normal).abs()
                > coplanarity_tolerance
        }) {
            return Err(ShellError::UnsupportedGeometry {
                what: "all nodes must be coplanar for the flat plate surrogate".into(),
            });
        }
        if let Some(support) = self.support {
            validate_support_geometry(self, support, reference_normal, geometric_scale)?;
        }
        if let DampingModel::Rayleigh {
            mass_proportional_per_s,
            stiffness_proportional_s,
        } = self.damping
        {
            if !mass_proportional_per_s.is_finite()
                || !stiffness_proportional_s.is_finite()
                || mass_proportional_per_s < 0.0
                || stiffness_proportional_s < 0.0
            {
                return Err(ShellError::InvalidInput {
                    what: "Rayleigh damping coefficients must be finite and non-negative".into(),
                });
            }
        }
        Ok(())
    }

    fn assemble_triangle(
        &self,
        triangle: [usize; 3],
        area: f64,
        edges: [([f64; 3], f64); 3],
        normal: [f64; 3],
        mass: &mut [f64],
        stiffness: &mut [f64],
    ) {
        let density = self.material.density_kg_m3;
        let thickness = self.thickness_m;
        let translational_mass = density * thickness * area / 3.0;
        let rotary_mass = density * thickness.powi(3) * area / 36.0;
        for node in triangle {
            for component in 0..3 {
                add(
                    mass,
                    self.nodes.len() * 6,
                    dof(node, component),
                    dof(node, component),
                    translational_mass,
                );
            }
            for component in 3..6 {
                add(
                    mass,
                    self.nodes.len() * 6,
                    dof(node, component),
                    dof(node, component),
                    rotary_mass,
                );
            }
        }
        // This is an edge-spring surrogate, not a CST plane-stress membrane.
        let membrane_scale = self.material.youngs_modulus_pa * thickness * area / 3.0;
        let bending = self.material.youngs_modulus_pa * thickness.powi(3)
            / (12.0 * (1.0 - self.material.poisson_ratio.powi(2)));
        for ((from, to), (edge, length)) in [
            (triangle[0], triangle[1]),
            (triangle[1], triangle[2]),
            (triangle[2], triangle[0]),
        ]
        .into_iter()
        .zip(edges)
        {
            let tangent = scale(edge, 1.0 / length);
            let cross_edge = cross(normal, tangent);
            let membrane = membrane_scale / (length * length);
            add_rank_one(
                stiffness,
                self.nodes.len() * 6,
                &[
                    (dof(from, 0), -tangent[0]),
                    (dof(from, 1), -tangent[1]),
                    (dof(from, 2), -tangent[2]),
                    (dof(to, 0), tangent[0]),
                    (dof(to, 1), tangent[1]),
                    (dof(to, 2), tangent[2]),
                ],
                membrane,
            );
            let slope = bending / 3.0;
            add_rank_one(
                stiffness,
                self.nodes.len() * 6,
                &[
                    (dof(from, 0), -normal[0] / length),
                    (dof(from, 1), -normal[1] / length),
                    (dof(from, 2), -normal[2] / length),
                    (dof(to, 0), normal[0] / length),
                    (dof(to, 1), normal[1] / length),
                    (dof(to, 2), normal[2] / length),
                    (dof(from, 3), cross_edge[0] * 0.5),
                    (dof(from, 4), cross_edge[1] * 0.5),
                    (dof(from, 5), cross_edge[2] * 0.5),
                    (dof(to, 3), cross_edge[0] * 0.5),
                    (dof(to, 4), cross_edge[1] * 0.5),
                    (dof(to, 5), cross_edge[2] * 0.5),
                ],
                slope,
            );
            for (direction, weight) in [(cross_edge, 1.0), (tangent, 1.0), (normal, 0.01)] {
                add_rank_one(
                    stiffness,
                    self.nodes.len() * 6,
                    &[
                        (dof(from, 3), -direction[0]),
                        (dof(from, 4), -direction[1]),
                        (dof(from, 5), -direction[2]),
                        (dof(to, 3), direction[0]),
                        (dof(to, 4), direction[1]),
                        (dof(to, 5), direction[2]),
                    ],
                    bending * weight / 3.0,
                );
            }
        }
        // A triangle-level normal-compatibility term closes the edge-slope
        // rank deficiency without penalizing a rigid motion:
        // grad(w) = n cross omega for w = n dot (omega cross x).
        let positions = triangle.map(|node| self.nodes[node].position_m);
        let gradients = [
            scale(
                cross(normal, sub(positions[2], positions[1])),
                1.0 / (2.0 * area),
            ),
            scale(
                cross(normal, sub(positions[0], positions[2])),
                1.0 / (2.0 * area),
            ),
            scale(
                cross(normal, sub(positions[1], positions[0])),
                1.0 / (2.0 * area),
            ),
        ];
        let tangent = scale(edges[0].0, 1.0 / edges[0].1);
        for direction in [tangent, cross(normal, tangent)] {
            let mut compatibility = Vec::with_capacity(12);
            for (node, gradient) in triangle.into_iter().zip(gradients) {
                let normal_gradient = dot(gradient, direction);
                compatibility.push((dof(node, 0), normal[0] * normal_gradient));
                compatibility.push((dof(node, 1), normal[1] * normal_gradient));
                compatibility.push((dof(node, 2), normal[2] * normal_gradient));
                let rotation = scale(cross(direction, normal), -1.0 / 3.0);
                compatibility.push((dof(node, 3), rotation[0]));
                compatibility.push((dof(node, 4), rotation[1]));
                compatibility.push((dof(node, 5), rotation[2]));
            }
            add_rank_one(
                stiffness,
                self.nodes.len() * 6,
                &compatibility,
                bending * area / (3.0 * thickness.powi(2)),
            );
        }
        // Couple the drilling rotation to the in-plane displacement curl.
        // For a rigid spin omega parallel to n, curl(u) = 2 omega dot n,
        // so this removes only the independent drilling mechanism.
        let bitangent = cross(normal, tangent);
        let mut drilling = Vec::with_capacity(12);
        for (node, gradient) in triangle.into_iter().zip(gradients) {
            let along_tangent = dot(gradient, tangent);
            let along_bitangent = dot(gradient, bitangent);
            drilling.push((dof(node, 0), 0.5 * along_bitangent * tangent[0]));
            drilling.push((dof(node, 1), 0.5 * along_bitangent * tangent[1]));
            drilling.push((dof(node, 2), 0.5 * along_bitangent * tangent[2]));
            drilling.push((dof(node, 0), -0.5 * along_tangent * bitangent[0]));
            drilling.push((dof(node, 1), -0.5 * along_tangent * bitangent[1]));
            drilling.push((dof(node, 2), -0.5 * along_tangent * bitangent[2]));
            drilling.push((dof(node, 3), normal[0] / 3.0));
            drilling.push((dof(node, 4), normal[1] / 3.0));
            drilling.push((dof(node, 5), normal[2] / 3.0));
        }
        add_rank_one(stiffness, self.nodes.len() * 6, &drilling, membrane_scale);
    }

    fn free_dofs(&self, dofs: usize) -> Result<Vec<usize>, ShellError> {
        let mut constrained = vec![false; dofs];
        if let Some(support) = self.support {
            for node in support.node_indices {
                constrained[dof(node, 0)] = true;
                constrained[dof(node, 1)] = true;
                constrained[dof(node, 2)] = true;
            }
        }
        Ok(constrained
            .into_iter()
            .enumerate()
            .filter_map(|(index, fixed)| (!fixed).then_some(index))
            .collect())
    }
}

fn matrix<K>(dimension: usize, values: Vec<f64>) -> ShellMatrix<K> {
    ShellMatrix {
        dimension,
        values,
        kind: PhantomData,
    }
}
fn dof(node: usize, component: usize) -> usize {
    node * 6 + component
}
fn add(values: &mut [f64], dimension: usize, row: usize, column: usize, value: f64) {
    values[row * dimension + column] += value;
}
fn add_rank_one(values: &mut [f64], dimension: usize, entries: &[(usize, f64)], scale: f64) {
    for &(row, row_value) in entries {
        for &(column, column_value) in entries {
            add(
                values,
                dimension,
                row,
                column,
                scale * row_value * column_value,
            );
        }
    }
}
fn reduce<K>(operator: &ShellMatrix<K>, free_dofs: &[usize]) -> ShellMatrix<K> {
    matrix(
        free_dofs.len(),
        free_dofs
            .iter()
            .flat_map(|&row| {
                free_dofs
                    .iter()
                    .map(move |&column| operator.values[row * operator.dimension + column])
            })
            .collect(),
    )
}

fn validate_derived_matrix(name: &str, values: &[f64]) -> Result<(), ShellError> {
    if values.iter().all(|value| value.is_finite()) {
        Ok(())
    } else {
        Err(ShellError::InvalidInput {
            what: format!("derived {name} matrix is non-finite"),
        })
    }
}

fn triangles_share_edge(left: [usize; 3], right: [usize; 3]) -> bool {
    left.iter().filter(|node| right.contains(*node)).count() >= 2
}

fn validate_support_geometry(
    plate: &ShellPlate,
    support: ShellSupport,
    reference_normal: [f64; 3],
    geometric_scale: f64,
) -> Result<(), ShellError> {
    if support
        .node_indices
        .iter()
        .any(|&index| index >= plate.nodes.len())
        || support.node_indices[0] == support.node_indices[1]
        || support.node_indices[1] == support.node_indices[2]
        || support.node_indices[0] == support.node_indices[2]
    {
        return Err(ShellError::UnsupportedBoundary {
            what: "three-point support requires three distinct in-range nodes".into(),
        });
    }
    let length = norm(support.normal);
    if !length.is_finite() || (length - 1.0).abs() > 1.0e-10 {
        return Err(ShellError::UnsupportedBoundary {
            what: "support normal must be unit length; oblique constraints need a multiplier formulation".into(),
        });
    }
    if dot(support.normal, reference_normal) < 1.0 - 1.0e-10 {
        return Err(ShellError::UnsupportedBoundary {
            what: "support normal must point to the declared flat plate normal side".into(),
        });
    }
    let [first, second, third] = support.node_indices;
    let a = plate.nodes[first].position_m;
    let b = plate.nodes[second].position_m;
    let c = plate.nodes[third].position_m;
    let distance_tolerance = 1.0e-12 * geometric_scale;
    if norm(sub(b, a)) <= distance_tolerance
        || norm(sub(c, b)) <= distance_tolerance
        || norm(sub(a, c)) <= distance_tolerance
        || norm(cross(sub(b, a), sub(c, a))) <= distance_tolerance * geometric_scale
    {
        return Err(ShellError::UnsupportedBoundary {
            what: "three-point support positions must be geometrically distinct and non-collinear"
                .into(),
        });
    }
    Ok(())
}

fn triangle_geometry(
    plate: &ShellPlate,
    triangle: [usize; 3],
) -> Result<(f64, [([f64; 3], f64); 3], [f64; 3]), ShellError> {
    let a = plate.nodes[triangle[0]].position_m;
    let b = plate.nodes[triangle[1]].position_m;
    let c = plate.nodes[triangle[2]].position_m;
    let ab = sub(b, a);
    let bc = sub(c, b);
    let ca = sub(a, c);
    let cross_product = cross(ab, sub(c, a));
    let twice_area = norm(cross_product);
    if !twice_area.is_finite() || twice_area <= 1e-12 {
        return Err(ShellError::InvalidInput {
            what: "triangle area must be finite and non-zero".into(),
        });
    }
    let edges = [(ab, norm(ab)), (bc, norm(bc)), (ca, norm(ca))];
    if edges
        .iter()
        .any(|(_, length)| !length.is_finite() || *length <= 1e-12)
    {
        return Err(ShellError::InvalidInput {
            what: "triangle edge length must be finite and positive".into(),
        });
    }
    Ok((
        twice_area * 0.5,
        edges,
        scale(cross_product, 1.0 / twice_area),
    ))
}

fn diagnostics(
    mass: &ShellMatrix<Mass>,
    stiffness: &ShellMatrix<Stiffness>,
    max_dofs: usize,
) -> OperatorDiagnostics {
    if mass.dimension > max_dofs {
        return OperatorDiagnostics::NotComputed {
            reason: format!(
                "{} DOFs exceeds max_conditioning_dofs {max_dofs}",
                mass.dimension
            ),
        };
    }
    let symmetry_residual = stiffness
        .values
        .iter()
        .enumerate()
        .map(|(index, value)| {
            (value
                - stiffness.values[(index % stiffness.dimension) * stiffness.dimension
                    + index / stiffness.dimension])
                .abs()
        })
        .fold(0.0_f64, f64::max);
    let mass_eigenvalues = jacobi_eigenvalues(&mass.values, mass.dimension);
    let stiffness_eigenvalues = jacobi_eigenvalues(&stiffness.values, stiffness.dimension);
    let scale = stiffness_eigenvalues
        .iter()
        .map(|value| value.abs())
        .fold(1.0_f64, f64::max);
    let tolerance = 1e-9 * scale;
    let mut positive: Vec<f64> = stiffness_eigenvalues
        .into_iter()
        .filter(|value| *value > tolerance)
        .collect();
    positive.sort_by(|left, right| left.total_cmp(right));
    OperatorDiagnostics::Computed {
        raw_mass_min_eigenvalue: mass_eigenvalues.into_iter().fold(f64::INFINITY, f64::min),
        raw_stiffness_nullity: stiffness.values.len() / stiffness.dimension - positive.len(),
        raw_stiffness_eigenvalue_spread: positive.first().map_or(f64::INFINITY, |minimum| {
            positive.last().copied().unwrap_or(*minimum) / minimum
        }),
        symmetry_residual,
    }
}

fn jacobi_eigenvalues(values: &[f64], dimension: usize) -> Vec<f64> {
    let mut matrix = values.to_vec();
    for _ in 0..(dimension * dimension * 32).max(1) {
        let mut pivot = (0, 0, 0.0_f64);
        for row in 0..dimension {
            for column in (row + 1)..dimension {
                let value = matrix[row * dimension + column].abs();
                if value > pivot.2 {
                    pivot = (row, column, value);
                }
            }
        }
        if pivot.2 <= 1e-12 {
            break;
        }
        let (p, q, _) = pivot;
        let app = matrix[p * dimension + p];
        let aqq = matrix[q * dimension + q];
        let apq = matrix[p * dimension + q];
        let angle = 0.5 * (2.0 * apq).atan2(aqq - app);
        let (sine, cosine) = angle.sin_cos();
        for index in 0..dimension {
            if index != p && index != q {
                let aip = matrix[index * dimension + p];
                let aiq = matrix[index * dimension + q];
                matrix[index * dimension + p] = cosine * aip - sine * aiq;
                matrix[p * dimension + index] = matrix[index * dimension + p];
                matrix[index * dimension + q] = sine * aip + cosine * aiq;
                matrix[q * dimension + index] = matrix[index * dimension + q];
            }
        }
        matrix[p * dimension + p] =
            cosine * cosine * app - 2.0 * sine * cosine * apq + sine * sine * aqq;
        matrix[q * dimension + q] =
            sine * sine * app + 2.0 * sine * cosine * apq + cosine * cosine * aqq;
        matrix[p * dimension + q] = 0.0;
        matrix[q * dimension + p] = 0.0;
    }
    (0..dimension)
        .map(|index| matrix[index * dimension + index])
        .collect()
}

fn sub(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
}
fn scale(value: [f64; 3], factor: f64) -> [f64; 3] {
    [value[0] * factor, value[1] * factor, value[2] * factor]
}
fn dot(left: [f64; 3], right: [f64; 3]) -> f64 {
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
}
fn norm(value: [f64; 3]) -> f64 {
    dot(value, value).sqrt()
}
fn cross(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]
}
