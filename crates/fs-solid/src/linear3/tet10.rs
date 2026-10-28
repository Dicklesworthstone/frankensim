//! Quadratic Lagrange displacement and geometry on tetrahedra (IRP M1.3-W1a).
//!
//! The canonical edge order is shared by promotion, assembly and importers.
//! Geometry is isoparametric; curved elements require a separate quadrature
//! ladder. Positive sampled Jacobians are admission checks, not a certificate
//! of injectivity between samples. Consistent mass is retained: row-sum lumping
//! would give negative corner masses. Thermal eigenstrain uses the shared
//! isoparametric load kernel; P1 arithmetic remains unchanged.

use std::collections::{BTreeMap, BTreeSet};

use fs_exec::Cx;
use fs_feec::highorder::simplex::duffy_quadrature;

use super::{
    QuadraticThermalStrainField, TetAssemblyBudget, TetElasticAssembly, TetElasticError,
    TetMaterialField, TetThermalLoad, determinant3, inverse3, maximum_edge_length,
    signed_tet_determinant, validate_connectivity,
};

/// `v0..v3`, then `(0,1),(1,2),(0,2),(0,3),(1,3),(2,3)`.
pub const TET10_EDGES: [[usize; 2]; 6] = [[0, 1], [1, 2], [0, 2], [0, 3], [1, 3], [2, 3]];

/// Bounded Duffy rule: increasing this is a declared quadrature axis.
pub const MAX_QUADRATURE_POINTS_PER_DIR: usize = 12;

/// Required quadrature choices, with no implicit underintegration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TetQuadrature {
    /// At least 3: affine stiffness has polynomial degree 2.
    pub stiffness_points_per_dir: usize,
    /// At least 4: affine consistent mass has polynomial degree 4.
    pub mass_points_per_dir: usize,
}

impl TetQuadrature {
    /// Duffy adds degree two in its first collapsed coordinate.
    pub fn validate(self) -> Result<(), TetElasticError> {
        if !(3..=MAX_QUADRATURE_POINTS_PER_DIR).contains(&self.stiffness_points_per_dir) {
            return Err(TetElasticError::QuadratureInadmissible {
                what: "stiffness Duffy points per direction must be in 3..=12",
            });
        }
        if !(4..=MAX_QUADRATURE_POINTS_PER_DIR).contains(&self.mass_points_per_dir) {
            return Err(TetElasticError::QuadratureInadmissible {
                what: "mass Duffy points per direction must be in 4..=12",
            });
        }
        Ok(())
    }
}

/// Fully specified quadratic elasticity problem, in SI units.
pub struct QuadraticTetProblem<'a> {
    /// Physical nodes [m], including shared midside nodes.
    pub nodes_m: &'a [[f64; 3]],
    /// Positive reference orientation, in [`TET10_EDGES`] order.
    pub elements: &'a [[usize; 10]],
    /// One resolved tensor per element, or one uniform tensor.
    pub materials: TetMaterialField<'a>,
    /// Homogeneous fixed node-major xyz displacement DOFs.
    pub fixed_dofs: &'a [usize],
    /// Caller-owned mesh, free-DOF and Jacobian quality limits.
    pub budget: TetAssemblyBudget,
    /// Explicit stiffness and consistent-mass rules.
    pub quadrature: TetQuadrature,
}

/// Reference shape functions and derivatives with respect to `(λ1,λ2,λ3)`.
/// The caller supplies barycentric coordinates summing to one.
#[must_use]
pub fn tet10_shapes(lambda: [f64; 4]) -> ([f64; 10], [[f64; 3]; 10]) {
    let dl = [[-1.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let mut n = [0.0; 10];
    let mut d = [[0.0; 3]; 10];
    for a in 0..4 {
        n[a] = lambda[a] * (2.0 * lambda[a] - 1.0);
        d[a] = dl[a].map(|v| (4.0 * lambda[a] - 1.0) * v);
    }
    for (edge, &[a, b]) in TET10_EDGES.iter().enumerate() {
        n[4 + edge] = 4.0 * lambda[a] * lambda[b];
        d[4 + edge] = core::array::from_fn(|k| 4.0 * (lambda[a] * dl[b][k] + lambda[b] * dl[a][k]));
    }
    (n, d)
}

/// Physical isoparametric evaluation, also usable for consistent loads and norms.
#[derive(Debug, Clone)]
pub struct Tet10Point {
    /// Interpolated position [m].
    pub position_m: [f64; 3],
    /// Nodal shape functions.
    pub shapes: [f64; 10],
    /// Physical shape gradients [1/m].
    pub gradients: [[f64; 3]; 10],
    /// Positive signed `det J` [m³]; multiply reference weights by this.
    pub determinant_m3: f64,
}

/// Evaluate a map, refusing folded/nonfinite geometry and the declared quality
/// floor. The indices are retained in any diagnostic.
pub fn tet10_point(
    nodes_m: [[f64; 3]; 10],
    lambda: [f64; 4],
    minimum_scaled_jacobian: f64,
    element: usize,
    point: usize,
) -> Result<Tet10Point, TetElasticError> {
    let (shapes, reference) = tet10_shapes(lambda);
    let mut position_m = nodes_m[0];
    let mut jac = [[0.0; 3]; 3];
    // Subtract the first vertex so translations do not contaminate J by
    // multiplying a rounded partition-of-unity residual by a large origin.
    for a in 0..10 {
        for i in 0..3 {
            let offset = nodes_m[a][i] - nodes_m[0][i];
            position_m[i] += shapes[a] * offset;
            for j in 0..3 {
                jac[i][j] += offset * reference[a][j];
            }
        }
    }
    let det = determinant3(jac);
    if !(det.is_finite() && det > 0.0) {
        return Err(TetElasticError::NonPositiveJacobian {
            element,
            point,
            det,
        });
    }
    let h = maximum_edge_length([nodes_m[0], nodes_m[1], nodes_m[2], nodes_m[3]]);
    let scaled = det / (h * h * h);
    if !(minimum_scaled_jacobian.is_finite() && minimum_scaled_jacobian > 0.0)
        || !scaled.is_finite()
        || scaled < minimum_scaled_jacobian
    {
        return Err(TetElasticError::DegenerateTet {
            element,
            scaled_jacobian: scaled,
            minimum: minimum_scaled_jacobian,
        });
    }
    let inv = inverse3(jac, det);
    let gradients = reference
        .map(|d| core::array::from_fn(|i| d[0] * inv[0][i] + d[1] * inv[1][i] + d[2] * inv[2][i]));
    if position_m
        .iter()
        .chain(gradients.iter().flatten())
        .any(|v| !v.is_finite())
    {
        return Err(TetElasticError::NonPositiveJacobian {
            element,
            point,
            det,
        });
    }
    Ok(Tet10Point {
        position_m,
        shapes,
        gradients,
        determinant_m3: det,
    })
}

impl QuadraticTetProblem<'_> {
    fn validate(&self, require_free_dofs: bool, cx: &Cx<'_>) -> Result<(), TetElasticError> {
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        self.quadrature.validate()?;
        if self.nodes_m.is_empty() || self.elements.is_empty() {
            return Err(TetElasticError::EmptyMesh);
        }
        for (what, requested, maximum) in [
            ("nodes", self.nodes_m.len(), self.budget.maximum_nodes),
            (
                "tetrahedra",
                self.elements.len(),
                self.budget.maximum_tetrahedra,
            ),
        ] {
            if requested > maximum {
                return Err(TetElasticError::BudgetExceeded {
                    what,
                    requested,
                    maximum,
                });
            }
        }
        let full = self
            .nodes_m
            .len()
            .checked_mul(3)
            .ok_or(TetElasticError::SizeOverflow)?;
        // Validate and count fixed DOFs before allocating either the maps or
        // matrix storage. O(n²) in constraints but no allocation on cap refusal.
        for (i, &dof) in self.fixed_dofs.iter().enumerate() {
            if i % 256 == 0 {
                cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            }
            if dof >= full {
                return Err(TetElasticError::FixedDofOutOfRange {
                    dof,
                    full_dofs: full,
                });
            }
            if self.fixed_dofs[..i].contains(&dof) {
                return Err(TetElasticError::DuplicateFixedDof { dof });
            }
        }
        let requested = full - self.fixed_dofs.len();
        if requested > self.budget.maximum_free_dofs {
            return Err(TetElasticError::BudgetExceeded {
                what: "free_dofs",
                requested,
                maximum: self.budget.maximum_free_dofs,
            });
        }
        if require_free_dofs && requested == 0 {
            return Err(TetElasticError::NoFreeDofs);
        }
        for (node, p) in self.nodes_m.iter().enumerate() {
            if node % 256 == 0 {
                cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            }
            if p.iter().any(|v| !v.is_finite()) {
                return Err(TetElasticError::NonFiniteNode { node });
            }
        }
        self.materials.validate(self.elements.len())?;
        let mut incident = vec![false; self.nodes_m.len()];
        let mut edges = BTreeMap::new();
        let mut midsides = BTreeMap::new();
        let mut vertices = BTreeSet::new();
        for tet in self.elements {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            vertices.extend(tet[..4].iter().copied());
        }
        for (element, tet) in self.elements.iter().enumerate() {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            for (a, &node) in tet.iter().enumerate() {
                if node >= self.nodes_m.len() {
                    return Err(TetElasticError::NodeIndexOutOfRange { element, node });
                }
                if tet[..a].contains(&node) {
                    return Err(TetElasticError::RepeatedTetNode { element });
                }
                incident[node] = true;
            }
            for (edge, &[a, b]) in TET10_EDGES.iter().enumerate() {
                let key = sorted_edge(tet[a], tet[b]);
                let mid = tet[4 + edge];
                if vertices.contains(&mid)
                    || edges.insert(key, mid).is_some_and(|old| old != mid)
                    || midsides.insert(mid, key).is_some_and(|old| old != key)
                {
                    return Err(TetElasticError::MidsideMismatch { element, edge });
                }
            }
        }
        if let Some(node) = incident.iter().position(|seen| !seen) {
            return Err(TetElasticError::UnreferencedNode { node });
        }
        Ok(())
    }

    /// Interpolate the original reference and current nodal temperatures at
    /// this problem's native stiffness points, then integrate each pinned
    /// source expansion law there. The returned borrow feeds the existing
    /// thermal load, response and tangent methods directly. Elastic tensors and
    /// reference mass remain those explicitly supplied in this problem.
    pub fn resolve_nodal_thermal_strains(
        &self,
        temperatures: super::source::NodalThermalTemperatures<'_>,
        sources: super::source::QuadraticThermalSources<'_>,
        mesh_frame: fs_blake3::ContentHash,
        strain_limit: f64,
        cx: &Cx<'_>,
    ) -> Result<super::source::ResolvedQuadraticThermalStrains, super::source::SolidMaterialError>
    {
        self.validate(false, cx)?;
        super::source::quadratic_thermal::resolve(
            self.nodes_m,
            self.elements,
            self.materials,
            temperatures,
            sources,
            &duffy_quadrature(self.quadrature.stiffness_points_per_dir),
            mesh_frame,
            strain_limit,
            |nodes, lambda, element, index| {
                tet10_point(
                    nodes,
                    lambda,
                    self.budget.minimum_scaled_jacobian,
                    element,
                    index,
                )
                .map(|point| point.shapes)
            },
            cx,
        )
    }

    /// Assemble into the same `(K,M)` type as P1. Each element and each
    /// quadrature point is a cancellation boundary. No partial result escapes.
    pub fn assemble(&self, cx: &Cx<'_>) -> Result<TetElasticAssembly, TetElasticError> {
        self.validate(true, cx)?;
        let kq = duffy_quadrature(self.quadrature.stiffness_points_per_dir);
        let mq = duffy_quadrature(self.quadrature.mass_points_per_dir);
        super::isoparametric::assemble::<10, 30, _>(
            self.nodes_m,
            self.elements,
            self.materials,
            self.fixed_dofs,
            self.budget,
            &kq,
            &mq,
            |nodes, lambda, element, point| {
                let p = tet10_point(
                    nodes,
                    lambda,
                    self.budget.minimum_scaled_jacobian,
                    element,
                    point,
                )?;
                let h = maximum_edge_length([nodes[0], nodes[1], nodes[2], nodes[3]]);
                Ok(super::isoparametric::Point {
                    shapes: p.shapes,
                    gradients: p.gradients,
                    determinant: p.determinant_m3,
                    scaled_jacobian: p.determinant_m3 / (h * h * h),
                })
            },
            cx,
        )
    }

    /// Integrate `B^T C epsilon_free` on the same quadratic geometry and
    /// stiffness quadrature as [`Self::assemble`]. Eigenstrain is an already
    /// integrated, source-bound thermal state, constant per element or supplied
    /// at each stiffness quadrature point. No temperature interpolation or
    /// expansion coefficient is inferred.
    ///
    /// The result is directly usable with the existing stiffness and solver
    /// owners: `K u = f_external + reduced_force_n`. Full force rows retain
    /// constrained contributions, including when every degree of freedom is fixed.
    /// Neither geometry nor material history changes,
    /// and a refusal/cancellation publishes no partial load. Quadratic response
    /// is available through solve_thermal_response; geometric prestress remains
    /// a separate capability.
    pub fn assemble_thermal_load<'t>(
        &self,
        thermal_strains: impl Into<QuadraticThermalStrainField<'t>>,
        cx: &Cx<'_>,
    ) -> Result<TetThermalLoad, TetElasticError> {
        self.validate(false, cx)?;
        super::isoparametric::thermal_load::<10, 30, _>(
            self.nodes_m,
            self.elements,
            self.materials,
            thermal_strains.into(),
            self.fixed_dofs,
            self.budget,
            &duffy_quadrature(self.quadrature.stiffness_points_per_dir),
            &duffy_quadrature(self.quadrature.mass_points_per_dir),
            |nodes, lambda, element, point| {
                let p = tet10_point(
                    nodes,
                    lambda,
                    self.budget.minimum_scaled_jacobian,
                    element,
                    point,
                )?;
                let h = maximum_edge_length([nodes[0], nodes[1], nodes[2], nodes[3]]);
                Ok(super::isoparametric::Point {
                    shapes: p.shapes,
                    gradients: p.gradients,
                    determinant: p.determinant_m3,
                    scaled_jacobian: p.determinant_m3 / (h * h * h),
                })
            },
            cx,
        )
    }

    /// Assemble thermal eigenstrain and dead nodal forces into one load.
    ///
    /// Forces are strictly increasing, unique full node-major xyz DOFs paired
    /// with newtons. Loads at fixed DOFs are retained for support reactions.
    pub fn assemble_thermal_load_with_nodal_forces<'t>(
        &self,
        thermal_strains: impl Into<QuadraticThermalStrainField<'t>>,
        nodal_forces_n: &[(usize, f64)],
        cx: &Cx<'_>,
    ) -> Result<TetThermalLoad, TetElasticError> {
        super::add_nodal_forces_to_thermal_load(
            self.assemble_thermal_load(thermal_strains, cx)?,
            self.nodes_m.len(),
            nodal_forces_n,
            cx,
        )
    }

    /// Solve the source-bound thermal equilibrium and recover physical stress,
    /// elastic energy and support forces on the same reference mesh.
    ///
    /// Uses the existing SPD solver with a true-residual gate. Every connected
    /// component must have its six rigid motions restrained. An entirely fixed
    /// body returns zero displacement with its thermal stress and reactions.
    /// Both total and elastic strains must satisfy the caller's small-strain
    /// ceiling at every stiffness quadrature point. Reported element stresses
    /// are volume averages; no pointwise or peak-stress accuracy is implied.
    /// The geometry, density and material state remain fixed through the solve.
    pub fn solve_thermal_response<'t>(
        &self,
        thermal_strains: impl Into<QuadraticThermalStrainField<'t>>,
        config: super::TetStaticSolveConfig,
        mesh_frame: fs_blake3::ContentHash,
        maximum_absolute_mandel_strain: f64,
        cx: &Cx<'_>,
    ) -> Result<super::QuadraticThermalResponse, TetElasticError> {
        self.solve_thermal_response_impl(
            thermal_strains.into(),
            None,
            config,
            mesh_frame,
            maximum_absolute_mandel_strain,
            cx,
        )
    }

    /// Solve thermal eigenstrain and applied dead nodal forces together.
    ///
    /// Forces use strictly increasing, unique full node-major xyz DOFs and
    /// finite newtons. Fixed-DOF loads change support reactions without creating
    /// artificial displacement. The response preserves actual internal forces
    /// and separately reports internal minus applied force: support reactions
    /// on fixed DOFs and unmodified equilibrium residuals on free DOFs.
    pub fn solve_thermal_response_with_nodal_forces<'t>(
        &self,
        thermal_strains: impl Into<QuadraticThermalStrainField<'t>>,
        nodal_forces_n: &[(usize, f64)],
        config: super::TetStaticSolveConfig,
        mesh_frame: fs_blake3::ContentHash,
        maximum_absolute_mandel_strain: f64,
        cx: &Cx<'_>,
    ) -> Result<super::QuadraticThermalResponse, TetElasticError> {
        self.solve_thermal_response_impl(
            thermal_strains.into(),
            Some(nodal_forces_n),
            config,
            mesh_frame,
            maximum_absolute_mandel_strain,
            cx,
        )
    }

    /// Solve thermal equilibrium and assemble its initial-stress vibration pencil.
    ///
    /// Actual stress is integrated at each original stiffness quadrature point,
    /// including its within-element variation. Reference mass, geometry and
    /// homogeneous incremental supports are unchanged. The tangent may be
    /// indefinite; signed eigenvalues belong to the existing spectral owner.
    /// This is the conservative initial-stress approximation, not a finite-strain
    /// tangent, follower-load model, stability certificate or live mode transfer.
    /// Fully fixed systems refuse because no vibration DOFs remain.
    pub fn solve_thermal_tangent<'t>(
        &self,
        thermal_strains: impl Into<QuadraticThermalStrainField<'t>>,
        config: super::TetStaticSolveConfig,
        mesh_frame: fs_blake3::ContentHash,
        maximum_absolute_mandel_strain: f64,
        cx: &Cx<'_>,
    ) -> Result<super::QuadraticPrestressedAssembly, TetElasticError> {
        self.solve_thermal_tangent_impl(
            thermal_strains.into(),
            None,
            config,
            mesh_frame,
            maximum_absolute_mandel_strain,
            cx,
        )
    }

    /// The same signed tangent about equilibrium with explicit dead nodal forces.
    ///
    /// The complete equilibrium is solved from these exact inputs in this call;
    /// an old response or independently authored stress cannot enter the pencil.
    /// Fixed-DOF forces retain their reaction semantics and add no follower-load
    /// stiffness. All material, strain, geometry and solve checks remain active.
    pub fn solve_thermal_tangent_with_nodal_forces<'t>(
        &self,
        thermal_strains: impl Into<QuadraticThermalStrainField<'t>>,
        nodal_forces_n: &[(usize, f64)],
        config: super::TetStaticSolveConfig,
        mesh_frame: fs_blake3::ContentHash,
        maximum_absolute_mandel_strain: f64,
        cx: &Cx<'_>,
    ) -> Result<super::QuadraticPrestressedAssembly, TetElasticError> {
        self.solve_thermal_tangent_impl(
            thermal_strains.into(),
            Some(nodal_forces_n),
            config,
            mesh_frame,
            maximum_absolute_mandel_strain,
            cx,
        )
    }

    fn solve_thermal_tangent_impl(
        &self,
        thermal_strains: QuadraticThermalStrainField<'_>,
        nodal_forces_n: Option<&[(usize, f64)]>,
        config: super::TetStaticSolveConfig,
        mesh_frame: fs_blake3::ContentHash,
        maximum_absolute_mandel_strain: f64,
        cx: &Cx<'_>,
    ) -> Result<super::QuadraticPrestressedAssembly, TetElasticError> {
        let response = self.solve_thermal_response_impl(
            thermal_strains,
            nodal_forces_n,
            config,
            mesh_frame,
            maximum_absolute_mandel_strain,
            cx,
        )?;
        let pencil = self.assemble(cx)?;
        super::isoparametric::thermal_response::prestress::assemble::<10, 30, _>(
            self.nodes_m,
            self.elements,
            self.materials,
            thermal_strains,
            pencil,
            response,
            &duffy_quadrature(self.quadrature.stiffness_points_per_dir),
            |nodes, lambda, element, point| {
                let p = tet10_point(
                    nodes,
                    lambda,
                    self.budget.minimum_scaled_jacobian,
                    element,
                    point,
                )?;
                let h = maximum_edge_length([nodes[0], nodes[1], nodes[2], nodes[3]]);
                Ok(super::isoparametric::Point {
                    shapes: p.shapes,
                    gradients: p.gradients,
                    determinant: p.determinant_m3,
                    scaled_jacobian: p.determinant_m3 / (h * h * h),
                })
            },
            maximum_absolute_mandel_strain,
            cx,
        )
    }

    fn solve_thermal_response_impl(
        &self,
        thermal_strains: QuadraticThermalStrainField<'_>,
        nodal_forces_n: Option<&[(usize, f64)]>,
        config: super::TetStaticSolveConfig,
        mesh_frame: fs_blake3::ContentHash,
        maximum_absolute_mandel_strain: f64,
        cx: &Cx<'_>,
    ) -> Result<super::QuadraticThermalResponse, TetElasticError> {
        self.validate(false, cx)?;
        super::isoparametric::thermal_response::solve::<10, 30, _>(
            self.nodes_m,
            self.elements,
            self.materials,
            thermal_strains,
            nodal_forces_n,
            self.fixed_dofs,
            self.budget,
            &duffy_quadrature(self.quadrature.stiffness_points_per_dir),
            &duffy_quadrature(self.quadrature.mass_points_per_dir),
            |nodes, lambda, element, point| {
                let p = tet10_point(
                    nodes,
                    lambda,
                    self.budget.minimum_scaled_jacobian,
                    element,
                    point,
                )?;
                let h = maximum_edge_length([nodes[0], nodes[1], nodes[2], nodes[3]]);
                Ok(super::isoparametric::Point {
                    shapes: p.shapes,
                    gradients: p.gradients,
                    determinant: p.determinant_m3,
                    scaled_jacobian: p.determinant_m3 / (h * h * h),
                })
            },
            config,
            mesh_frame,
            maximum_absolute_mandel_strain,
            cx,
        )
    }
}

/// A caller's source-surface projection; no surface or material is guessed.
pub trait BoundaryProjection {
    /// Project a boundary-edge midpoint [m], or refuse.
    fn project(&self, p: [f64; 3]) -> Result<[f64; 3], TetElasticError>;
}

/// A generated part's exact-surface midside table (IRP M1.3-W4). A point that
/// is not a boundary-edge midpoint of the tabulated mesh refuses.
impl BoundaryProjection for fs_mesh::solid_parts::MidsideProjection {
    fn project(&self, p: [f64; 3]) -> Result<[f64; 3], TetElasticError> {
        self.get(p).ok_or(TetElasticError::InvalidGeometryUpdate {
            what: "boundary midpoint is not tabulated by the exact-surface projection",
        })
    }
}

fn sorted_edge(a: usize, b: usize) -> [usize; 2] {
    [a.min(b), a.max(b)]
}

/// Promote conforming P1 topology with one globally shared node per edge.
/// New IDs follow sorted global endpoint pairs, independent of element order.
/// Only edges on boundary faces are projected. Negative input orientations are
/// canonicalized to positive; the original node IDs are preserved.
pub fn promote_p1_to_p2(
    nodes_m: &[[f64; 3]],
    tets: &[[usize; 4]],
    boundary: Option<&dyn BoundaryProjection>,
    cx: &Cx<'_>,
) -> Result<(Vec<[f64; 3]>, Vec<[usize; 10]>), TetElasticError> {
    cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
    if nodes_m.is_empty() || tets.is_empty() {
        return Err(TetElasticError::EmptyMesh);
    }
    for (node, p) in nodes_m.iter().enumerate() {
        if p.iter().any(|v| !v.is_finite()) {
            return Err(TetElasticError::NonFiniteNode { node });
        }
    }
    let mut edges = BTreeSet::new();
    let mut faces = BTreeMap::<[usize; 3], usize>::new();
    for (element, tet) in tets.iter().enumerate() {
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        validate_connectivity(element, tet, nodes_m.len())?;
        let det = signed_tet_determinant(tet.map(|n| nodes_m[n]));
        if !det.is_finite() || det == 0.0 {
            return Err(TetElasticError::NonPositiveJacobian {
                element,
                point: 0,
                det,
            });
        }
        for [a, b] in TET10_EDGES {
            edges.insert(sorted_edge(tet[a], tet[b]));
        }
        for omit in 0..4 {
            let mut face = [0; 3];
            let mut j = 0;
            for (i, &v) in tet.iter().enumerate() {
                if i != omit {
                    face[j] = v;
                    j += 1;
                }
            }
            face.sort_unstable();
            *faces.entry(face).or_default() += 1;
        }
    }
    let mut boundary_edges = BTreeSet::new();
    for (face, count) in faces {
        if count == 1 {
            for [a, b] in [[0, 1], [1, 2], [0, 2]] {
                boundary_edges.insert(sorted_edge(face[a], face[b]));
            }
        } else if count > 2 {
            return Err(TetElasticError::InvalidGeometryUpdate {
                what: "nonmanifold P1 face during quadratic promotion",
            });
        }
    }
    let extent = nodes_m
        .len()
        .checked_add(edges.len())
        .ok_or(TetElasticError::SizeOverflow)?;
    let mut nodes = Vec::with_capacity(extent);
    nodes.extend_from_slice(nodes_m);
    let mut edge_nodes = BTreeMap::new();
    for key in edges {
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        let mut p = core::array::from_fn(|k| f64::midpoint(nodes_m[key[0]][k], nodes_m[key[1]][k]));
        if boundary_edges.contains(&key) {
            if let Some(projector) = boundary {
                p = projector.project(p)?;
            }
        }
        if p.iter().any(|v| !v.is_finite()) {
            return Err(TetElasticError::NonFiniteNode { node: nodes.len() });
        }
        edge_nodes.insert(key, nodes.len());
        nodes.push(p);
    }
    let mut elements = Vec::with_capacity(tets.len());
    for tet in tets {
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        let mut tet = *tet;
        if signed_tet_determinant(tet.map(|n| nodes_m[n])) < 0.0 {
            tet.swap(1, 2);
        }
        let mut quadratic = [0; 10];
        quadratic[..4].copy_from_slice(&tet);
        for (edge, [a, b]) in TET10_EDGES.into_iter().enumerate() {
            quadratic[4 + edge] = edge_nodes[&sorted_edge(tet[a], tet[b])];
        }
        elements.push(quadratic);
    }
    Ok((nodes, elements))
}
