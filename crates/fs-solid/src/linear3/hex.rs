//! Fully integrated isoparametric HEX20/HEX27. The element-specific pieces are
//! shapes, geometry and topology admission; TET10 owns the shared numerical
//! assembly path in isoparametric.rs. Positive sampled Jacobians do not prove
//! global injectivity. Materials are constant per element; mass is consistent.
use super::{
    QuadraticThermalStrainField, TetAssemblyBudget, TetElasticAssembly, TetElasticError,
    TetMaterialField, TetThermalLoad, determinant3, inverse3,
};
use fs_exec::Cx;
use fs_feec::highorder::quad1d::gauss_legendre;
pub use fs_mesh::structured_hex::{
    GMSH20_TO_CANONICAL, GMSH27_TO_CANONICAL, HEX_EDGES, HEX20_REFERENCE_NODES,
    HEX27_REFERENCE_NODES, HexFamily,
};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Copy)]
pub enum HexConnectivity<'a> {
    Hex20(&'a [[usize; 20]]),
    Hex27(&'a [[usize; 27]]),
}
#[derive(Debug, Clone, Copy)]
pub struct HexQuadrature {
    pub points_per_dir: usize,
}
impl HexQuadrature {
    pub fn validate(self) -> Result<(), TetElasticError> {
        if self.points_per_dir < 3 {
            return Err(TetElasticError::UnderIntegrated {
                points_per_dir: self.points_per_dir,
            });
        }
        if self.points_per_dir > 12 {
            return Err(TetElasticError::QuadratureInadmissible {
                what: "hex Gauss points per direction must be in 3..=12",
            });
        }
        Ok(())
    }
    pub fn points(self) -> Result<Vec<([f64; 3], f64)>, TetElasticError> {
        self.validate()?;
        let (x, w) = gauss_legendre(self.points_per_dir);
        let mut rule = Vec::with_capacity(x.len() * x.len() * x.len());
        for i in 0..x.len() {
            for j in 0..x.len() {
                for k in 0..x.len() {
                    rule.push(([x[i], x[j], x[k]], w[i] * w[j] * w[k]));
                }
            }
        }
        Ok(rule)
    }
}
pub struct HexElasticProblem<'a> {
    pub nodes_m: &'a [[f64; 3]],
    pub family: HexFamily,
    pub elements: HexConnectivity<'a>,
    pub materials: TetMaterialField<'a>,
    pub fixed_dofs: &'a [usize],
    pub budget: TetAssemblyBudget,
    pub quadrature: HexQuadrature,
}

/// Quadratic serendipity shapes and reference gradients on [-1,1]^3.
pub fn hex20_shapes(p: [f64; 3]) -> ([f64; 20], [[f64; 3]; 20]) {
    let mut n = [0.0; 20];
    let mut d = [[0.0; 3]; 20];
    for (a, q) in HEX20_REFERENCE_NODES.iter().enumerate() {
        if a < 8 {
            let s = q.map(f64::from);
            let f = core::array::from_fn::<_, 3, _>(|i| 1.0 + s[i] * p[i]);
            let g = (0..3).map(|i| s[i] * p[i]).sum::<f64>() - 2.0;
            n[a] = f[0] * f[1] * f[2] * g / 8.0;
            for i in 0..3 {
                let j = (i + 1) % 3;
                let k = (i + 2) % 3;
                d[a][i] = s[i] * f[j] * f[k] * (g + f[i]) / 8.0;
            }
        } else {
            let axis = q
                .iter()
                .position(|v| *v == 0)
                .expect("canonical edge has one zero coordinate");
            let j = (axis + 1) % 3;
            let k = (axis + 2) % 3;
            let sj = f64::from(q[j]);
            let sk = f64::from(q[k]);
            let fj = 1.0 + sj * p[j];
            let fk = 1.0 + sk * p[k];
            let bubble = 1.0 - p[axis] * p[axis];
            n[a] = bubble * fj * fk / 4.0;
            d[a][axis] = -p[axis] * fj * fk / 2.0;
            d[a][j] = bubble * sj * fk / 4.0;
            d[a][k] = bubble * fj * sk / 4.0;
        }
    }
    (n, d)
}
/// Tensor Q2 Lagrange shapes, in the same corner/edge order as HEX20.
pub fn hex27_shapes(p: [f64; 3]) -> ([f64; 27], [[f64; 3]; 27]) {
    let mut n = [0.0; 27];
    let mut d = [[0.0; 3]; 27];
    for (a, q) in HEX27_REFERENCE_NODES.iter().enumerate() {
        let mut l = [0.0; 3];
        let mut dl = [0.0; 3];
        for i in 0..3 {
            (l[i], dl[i]) = match q[i] {
                -1 => (p[i] * (p[i] - 1.0) / 2.0, p[i] - 0.5),
                0 => (1.0 - p[i] * p[i], -2.0 * p[i]),
                1 => (p[i] * (p[i] + 1.0) / 2.0, p[i] + 0.5),
                _ => unreachable!("canonical Q2 coordinate"),
            };
        }
        n[a] = l[0] * l[1] * l[2];
        for i in 0..3 {
            d[a][i] = dl[i] * l[(i + 1) % 3] * l[(i + 2) % 3];
        }
    }
    (n, d)
}
#[derive(Debug, Clone)]
pub struct HexPoint<const N: usize> {
    pub position_m: [f64; 3],
    pub shapes: [f64; N],
    pub gradients: [[f64; 3]; N],
    pub determinant_m3: f64,
    pub scaled_jacobian: f64,
}
fn point<const N: usize>(
    nodes: [[f64; 3]; N],
    shapes: [f64; N],
    reference: [[f64; 3]; N],
    minimum: f64,
    element: usize,
    point: usize,
) -> Result<HexPoint<N>, TetElasticError> {
    let mut position = nodes[0];
    let mut jac = [[0.0; 3]; 3];
    for a in 0..N {
        for i in 0..3 {
            let offset = nodes[a][i] - nodes[0][i];
            position[i] += shapes[a] * offset;
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
    let mut h2 = 0.0f64;
    for a in 0..8 {
        for b in a + 1..8 {
            h2 = h2.max(
                (0..3)
                    .map(|k| (nodes[a][k] - nodes[b][k]) * (nodes[a][k] - nodes[b][k]))
                    .sum(),
            );
        }
    }
    let h = fs_math::det::sqrt(h2);
    let scaled = det / (h * h * h);
    if !(minimum.is_finite() && minimum > 0.0 && scaled.is_finite() && scaled >= minimum) {
        return Err(TetElasticError::DegenerateTet {
            element,
            scaled_jacobian: scaled,
            minimum,
        });
    }
    let inv = inverse3(jac, det);
    let gradients = reference
        .map(|g| core::array::from_fn(|i| g[0] * inv[0][i] + g[1] * inv[1][i] + g[2] * inv[2][i]));
    if position
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
    Ok(HexPoint {
        position_m: position,
        shapes,
        gradients,
        determinant_m3: det,
        scaled_jacobian: scaled,
    })
}
pub fn hex20_point(
    nodes: [[f64; 3]; 20],
    p: [f64; 3],
    minimum: f64,
    element: usize,
    index: usize,
) -> Result<HexPoint<20>, TetElasticError> {
    let (n, d) = hex20_shapes(p);
    point(nodes, n, d, minimum, element, index)
}
pub fn hex27_point(
    nodes: [[f64; 3]; 27],
    p: [f64; 3],
    minimum: f64,
    element: usize,
    index: usize,
) -> Result<HexPoint<27>, TetElasticError> {
    let (n, d) = hex27_shapes(p);
    point(nodes, n, d, minimum, element, index)
}

impl HexElasticProblem<'_> {
    fn validate<const N: usize>(
        &self,
        elements: &[[usize; N]],
        require_free_dofs: bool,
        cx: &Cx<'_>,
    ) -> Result<(), TetElasticError> {
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        self.quadrature.validate()?;
        let expected = match self.family {
            HexFamily::Serendipity20 => 20,
            HexFamily::Lagrange27 => 27,
        };
        if expected != N {
            return Err(TetElasticError::HexConnectivityMismatch { expected, got: N });
        }
        if self.nodes_m.is_empty() || elements.is_empty() {
            return Err(TetElasticError::EmptyMesh);
        }
        for (what, requested, maximum) in [
            ("nodes", self.nodes_m.len(), self.budget.maximum_nodes),
            ("tetrahedra", elements.len(), self.budget.maximum_tetrahedra),
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
        self.materials.validate(elements.len())?;
        let mut corners = BTreeSet::new();
        for e in elements {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            corners.extend(e[..8].iter().copied());
        }
        let mut entities = BTreeMap::<Vec<usize>, usize>::new();
        let mut roles = BTreeMap::<usize, Vec<usize>>::new();
        let mut incident = vec![false; self.nodes_m.len()];
        let mut cells = BTreeSet::new();
        const FACES: [[usize; 4]; 6] = [
            [0, 1, 2, 3],
            [4, 5, 6, 7],
            [0, 1, 5, 4],
            [1, 2, 6, 5],
            [2, 3, 7, 6],
            [3, 0, 4, 7],
        ];
        for (element, e) in elements.iter().enumerate() {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            for (a, &node) in e.iter().enumerate() {
                if node >= self.nodes_m.len() {
                    return Err(TetElasticError::NodeIndexOutOfRange { element, node });
                }
                if e[..a].contains(&node) {
                    return Err(TetElasticError::RepeatedTetNode { element });
                }
                incident[node] = true;
            }
            let mut cell = e[..8].to_vec();
            cell.sort_unstable();
            if !cells.insert(cell.clone()) {
                return Err(TetElasticError::InvalidGeometryUpdate {
                    what: "duplicate hexahedron",
                });
            }
            for local in 8..N {
                let mut key = if local < 20 {
                    HEX_EDGES[local - 8].map(|i| e[i]).to_vec()
                } else if local < 26 {
                    FACES[local - 20].map(|i| e[i]).to_vec()
                } else {
                    cell.clone()
                };
                key.sort_unstable();
                let node = e[local];
                if corners.contains(&node)
                    || entities
                        .insert(key.clone(), node)
                        .is_some_and(|old| old != node)
                    || roles
                        .insert(node, key.clone())
                        .is_some_and(|old| old != key)
                {
                    return Err(TetElasticError::MidsideMismatch {
                        element,
                        edge: local - 8,
                    });
                }
            }
        }
        if let Some(node) = incident.iter().position(|v| !*v) {
            return Err(TetElasticError::UnreferencedNode { node });
        }
        Ok(())
    }
    pub fn assemble(&self, cx: &Cx<'_>) -> Result<TetElasticAssembly, TetElasticError> {
        match self.elements {
            HexConnectivity::Hex20(elements) => {
                self.validate(elements, true, cx)?;
                let q = self.quadrature.points()?;
                super::isoparametric::assemble::<20, 60, _>(
                    self.nodes_m,
                    elements,
                    self.materials,
                    self.fixed_dofs,
                    self.budget,
                    &q,
                    &q,
                    |nodes, r, e, i| {
                        let p = hex20_point(nodes, r, self.budget.minimum_scaled_jacobian, e, i)?;
                        Ok(super::isoparametric::Point {
                            shapes: p.shapes,
                            gradients: p.gradients,
                            determinant: p.determinant_m3,
                            scaled_jacobian: p.scaled_jacobian,
                        })
                    },
                    cx,
                )
            }
            HexConnectivity::Hex27(elements) => {
                self.validate(elements, true, cx)?;
                let q = self.quadrature.points()?;
                super::isoparametric::assemble::<27, 81, _>(
                    self.nodes_m,
                    elements,
                    self.materials,
                    self.fixed_dofs,
                    self.budget,
                    &q,
                    &q,
                    |nodes, r, e, i| {
                        let p = hex27_point(nodes, r, self.budget.minimum_scaled_jacobian, e, i)?;
                        Ok(super::isoparametric::Point {
                            shapes: p.shapes,
                            gradients: p.gradients,
                            determinant: p.determinant_m3,
                            scaled_jacobian: p.scaled_jacobian,
                        })
                    },
                    cx,
                )
            }
        }
    }

    /// Thermal eigenstrain load on the existing HEX20/HEX27 geometry. Uses
    /// the same material tensor, fully integrated rule, constraints and quality
    /// admission as [`Self::assemble`]. Each element's integrated symmetric
    /// eigenstrain must refer to that element's exact elastic material state.
    /// It may be constant per element or resolved at every stiffness point;
    /// gradient heating is retained through stress, energy and prestress recovery.
    ///
    /// The shared load carrier supplies `K u = f_external + f_thermal` to the
    /// existing solver and retains full rows at supports, including fully fixed
    /// meshes with an empty reduced force vector. This does not update
    /// geometry or imply an error bound on curved-element quadrature. Physical
    /// response is available through solve_thermal_response. Cancellation never
    /// publishes partial forces.
    pub fn assemble_thermal_load<'t>(
        &self,
        thermal_strains: impl Into<QuadraticThermalStrainField<'t>>,
        cx: &Cx<'_>,
    ) -> Result<TetThermalLoad, TetElasticError> {
        let thermal_strains = thermal_strains.into();
        match self.elements {
            HexConnectivity::Hex20(elements) => {
                self.validate(elements, false, cx)?;
                let q = self.quadrature.points()?;
                super::isoparametric::thermal_load::<20, 60, _>(
                    self.nodes_m,
                    elements,
                    self.materials,
                    thermal_strains,
                    self.fixed_dofs,
                    self.budget,
                    &q,
                    &q,
                    |nodes, r, e, i| {
                        let p = hex20_point(nodes, r, self.budget.minimum_scaled_jacobian, e, i)?;
                        Ok(super::isoparametric::Point {
                            shapes: p.shapes,
                            gradients: p.gradients,
                            determinant: p.determinant_m3,
                            scaled_jacobian: p.scaled_jacobian,
                        })
                    },
                    cx,
                )
            }
            HexConnectivity::Hex27(elements) => {
                self.validate(elements, false, cx)?;
                let q = self.quadrature.points()?;
                super::isoparametric::thermal_load::<27, 81, _>(
                    self.nodes_m,
                    elements,
                    self.materials,
                    thermal_strains,
                    self.fixed_dofs,
                    self.budget,
                    &q,
                    &q,
                    |nodes, r, e, i| {
                        let p = hex27_point(nodes, r, self.budget.minimum_scaled_jacobian, e, i)?;
                        Ok(super::isoparametric::Point {
                            shapes: p.shapes,
                            gradients: p.gradients,
                            determinant: p.determinant_m3,
                            scaled_jacobian: p.scaled_jacobian,
                        })
                    },
                    cx,
                )
            }
        }
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

    /// Solve the source-bound thermal equilibrium with the existing SPD solver,
    /// then integrate actual stress, stored energy and full support forces.
    ///
    /// Every connected component must have its six rigid motions restrained.
    /// Fully fixed bodies retain zero displacement and nonzero thermal reactions.
    /// Total and elastic strains are checked against the caller's ceiling at
    /// each stiffness quadrature point; element stresses are volume averages.
    /// This frozen-reference solve retains geometry and density and does not
    /// add geometric stiffness or certify curved-map integration accuracy.
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
        let q = self.quadrature.points()?;
        match self.elements {
            HexConnectivity::Hex20(elements) => {
                super::isoparametric::thermal_response::prestress::assemble::<20, 60, _>(
                    self.nodes_m,
                    elements,
                    self.materials,
                    thermal_strains,
                    pencil,
                    response,
                    &q,
                    |nodes, r, e, i| {
                        let p = hex20_point(nodes, r, self.budget.minimum_scaled_jacobian, e, i)?;
                        Ok(super::isoparametric::Point {
                            shapes: p.shapes,
                            gradients: p.gradients,
                            determinant: p.determinant_m3,
                            scaled_jacobian: p.scaled_jacobian,
                        })
                    },
                    maximum_absolute_mandel_strain,
                    cx,
                )
            }
            HexConnectivity::Hex27(elements) => {
                super::isoparametric::thermal_response::prestress::assemble::<27, 81, _>(
                    self.nodes_m,
                    elements,
                    self.materials,
                    thermal_strains,
                    pencil,
                    response,
                    &q,
                    |nodes, r, e, i| {
                        let p = hex27_point(nodes, r, self.budget.minimum_scaled_jacobian, e, i)?;
                        Ok(super::isoparametric::Point {
                            shapes: p.shapes,
                            gradients: p.gradients,
                            determinant: p.determinant_m3,
                            scaled_jacobian: p.scaled_jacobian,
                        })
                    },
                    maximum_absolute_mandel_strain,
                    cx,
                )
            }
        }
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
        match self.elements {
            HexConnectivity::Hex20(elements) => {
                self.validate(elements, false, cx)?;
                let q = self.quadrature.points()?;
                super::isoparametric::thermal_response::solve::<20, 60, _>(
                    self.nodes_m,
                    elements,
                    self.materials,
                    thermal_strains,
                    nodal_forces_n,
                    self.fixed_dofs,
                    self.budget,
                    &q,
                    &q,
                    |nodes, r, e, i| {
                        let p = hex20_point(nodes, r, self.budget.minimum_scaled_jacobian, e, i)?;
                        Ok(super::isoparametric::Point {
                            shapes: p.shapes,
                            gradients: p.gradients,
                            determinant: p.determinant_m3,
                            scaled_jacobian: p.scaled_jacobian,
                        })
                    },
                    config,
                    mesh_frame,
                    maximum_absolute_mandel_strain,
                    cx,
                )
            }
            HexConnectivity::Hex27(elements) => {
                self.validate(elements, false, cx)?;
                let q = self.quadrature.points()?;
                super::isoparametric::thermal_response::solve::<27, 81, _>(
                    self.nodes_m,
                    elements,
                    self.materials,
                    thermal_strains,
                    nodal_forces_n,
                    self.fixed_dofs,
                    self.budget,
                    &q,
                    &q,
                    |nodes, r, e, i| {
                        let p = hex27_point(nodes, r, self.budget.minimum_scaled_jacobian, e, i)?;
                        Ok(super::isoparametric::Point {
                            shapes: p.shapes,
                            gradients: p.gradients,
                            determinant: p.determinant_m3,
                            scaled_jacobian: p.scaled_jacobian,
                        })
                    },
                    config,
                    mesh_frame,
                    maximum_absolute_mandel_strain,
                    cx,
                )
            }
        }
    }
}
