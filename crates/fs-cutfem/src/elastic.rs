//! Vector Q1 small-strain elasticity on certified SDF cuts.
//!
//! This is the vector sibling of [`crate::fem::Space`]. It reuses the
//! certified [`crate::quad::cut_cell_rules`] surface verbatim, assembles
//! the symmetric Nitsche form for displacement data, and stabilizes cut
//! faces with a first-normal-derivative ghost penalty. The Nitsche
//! constant is deliberately **independent of cut fraction**:
//! `beta * mu / h`. Degenerating cuts are controlled by the ghost
//! penalty, not hidden behind an exploding boundary penalty.
//!
//! The current vector lift is restricted to a uniform active quadtree
//! level. That refusal is explicit because silently ignoring scalar
//! hanging-node expansions would make a nonconforming elasticity space.

use crate::CutFemError;
use crate::fem::{JacobiPrecond, q1};
use crate::grid::{CellKey, NodeKey, Quadtree};
use crate::quad::{CutRules, cut_cell_rules, tensor_gauss};
use crate::sdf::CutSdf;
use fs_material::IsotropicElastic;
use fs_solver::krylov::CgState;
use fs_solver::op::LinearOp;
use fs_sparse::{Coo, Csr};
use std::collections::{BTreeMap, BTreeSet};

/// Stable prefix for content-addressed linear displacement apply VJP keys.
#[cfg(feature = "adjoint-vjp")]
pub const ELASTICITY_APPLY_VJP_OP: &str = "fs-cutfem.elasticity-apply.v1";

/// Largest plane-strain constitutive stiffness ratio certified by the vector
/// frontend, `(lambda + 2*mu) / mu`.
///
/// The first-generation Nitsche and ghost terms scale with `mu`. Capping the
/// ratio at four keeps the admitted material family inside the compressible
/// regime exercised by the coercivity battery instead of silently extending
/// that evidence to the nearly incompressible limit.
pub const MAX_PLANE_STRAIN_STIFFNESS_RATIO: f64 = 4.0;

type FaceKey = (CellKey, CellKey);

/// Vector Q1 CutFEM problem on `Omega = {phi < 0}`.
///
/// The constitutive parameters come from [`IsotropicElastic`], so the
/// material's admissibility checks and model-card identity are shared
/// with the rest of FLUX rather than duplicated here.
pub struct CutElasticity<'a> {
    /// Uniform background quadtree.
    pub grid: &'a Quadtree,
    /// Certified negative-inside level set.
    pub sdf: &'a dyn CutSdf,
    /// Isotropic small-strain material (plane-strain restriction). The v1
    /// certified regime requires `(lambda + 2*mu) / mu <= 4`.
    pub material: &'a IsotropicElastic,
    /// Dimensionless symmetric-Nitsche constant. The applied penalty is
    /// `nitsche_beta * mu / h`, never divided by cut fraction.
    pub nitsche_beta: f64,
    /// First-derivative ghost-penalty constant. Zero disables it.
    pub ghost_gamma: f64,
    /// Certified cut-quadrature subdivision depth.
    pub quad_depth: u32,
    /// Optional zero-displacement clamp on active design-box boundary nodes.
    pub clamp: Option<&'a dyn Fn(f64, f64) -> bool>,
    /// Optional dead traction on active design-box boundary edges.
    pub boundary_traction: Option<&'a dyn Fn(f64, f64) -> [f64; 2]>,
    /// Use a natural traction-free embedded interface instead of Nitsche
    /// displacement data. A clamp is then normally required to remove
    /// rigid-body modes.
    pub traction_free_interface: bool,
    /// CG relative-residual target.
    pub solver_tol: f64,
    /// CG iteration cap.
    pub solver_max_iters: usize,
}

impl core::fmt::Debug for CutElasticity<'_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("CutElasticity")
            .field("material", self.material)
            .field("nitsche_beta", &self.nitsche_beta)
            .field("ghost_gamma", &self.ghost_gamma)
            .field("quad_depth", &self.quad_depth)
            .field("has_clamp", &self.clamp.is_some())
            .field("has_boundary_traction", &self.boundary_traction.is_some())
            .field("traction_free_interface", &self.traction_free_interface)
            .field("solver_tol", &self.solver_tol)
            .field("solver_max_iters", &self.solver_max_iters)
            .finish_non_exhaustive()
    }
}

/// Assembled vector elasticity operator and its deterministic topology.
///
/// The operator exposes apply/transpose-apply separately from the solve.
/// That keeps adjoints on the operator boundary rather than differentiating
/// through CG iterations.
#[derive(Debug, Clone)]
pub struct CutElasticityOperator {
    matrix: Csr,
    rhs: Vec<f64>,
    node_ids: BTreeMap<NodeKey, usize>,
    clamped: Vec<bool>,
    active: Vec<CellKey>,
    rules: BTreeMap<CellKey, CutRules>,
    dropped_cut_cells: usize,
}

impl CutElasticityOperator {
    /// Canonical symmetric CSR matrix.
    #[must_use]
    pub fn matrix(&self) -> &Csr {
        &self.matrix
    }

    /// Assembled load vector.
    #[must_use]
    pub fn rhs(&self) -> &[f64] {
        &self.rhs
    }

    /// Vector displacement DOF count.
    #[must_use]
    pub fn dof_count(&self) -> usize {
        self.matrix.nrows()
    }

    /// Deterministic node-to-block map. Node `id` owns displacement DOFs
    /// `2*id` and `2*id + 1`.
    #[must_use]
    pub fn node_ids(&self) -> &BTreeMap<NodeKey, usize> {
        &self.node_ids
    }

    /// Per-DOF zero-clamp mask. Clamped rows remain in the operator as unit
    /// identity rows rather than being eliminated from the coefficient vector.
    #[must_use]
    pub fn clamped_dofs(&self) -> &[bool] {
        &self.clamped
    }

    /// Conservatively classified cut cells whose quadrature retained less
    /// than `1e-12` of the full-cell area and were therefore omitted.
    #[must_use]
    pub fn dropped_cut_cells(&self) -> usize {
        self.dropped_cut_cells
    }

    /// Apply `y = K x` in canonical CSR column order.
    #[must_use]
    pub fn apply_vec(&self, x: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0; self.matrix.nrows()];
        self.matrix.spmv(x, &mut y);
        y
    }

    /// Apply `y = K^T x`. The assembled symmetric form makes this
    /// bit-identical to [`Self::apply_vec`].
    #[must_use]
    pub fn apply_transpose_vec(&self, x: &[f64]) -> Vec<f64> {
        self.apply_vec(x)
    }

    /// Expand one coefficient vector into deterministic nodal displacements.
    #[must_use]
    pub fn nodal_values(&self, x: &[f64]) -> BTreeMap<NodeKey, [f64; 2]> {
        assert_eq!(x.len(), self.dof_count(), "elasticity coefficient length");
        self.node_ids
            .iter()
            .map(|(&node, &id)| (node, [x[2 * id], x[2 * id + 1]]))
            .collect()
    }
}

impl LinearOp for CutElasticityOperator {
    fn n(&self) -> usize {
        self.dof_count()
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.matrix.spmv(x, y);
    }

    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        // `symmetrize_local` makes every element pair bit-identical before
        // canonical COO accumulation, so the assembled CSR is exactly
        // symmetric rather than merely symmetric up to roundoff.
        self.matrix.spmv(x, y);
    }
}

/// Solved vector field plus convergence and integration metadata.
#[derive(Debug, Clone)]
pub struct CutElasticitySolution {
    coefficients: Vec<f64>,
    nodal: BTreeMap<NodeKey, [f64; 2]>,
    active: Vec<CellKey>,
    rules: BTreeMap<CellKey, CutRules>,
    dropped_cut_cells: usize,
    /// CG iterations.
    pub iters: usize,
    /// Final relative residual.
    pub rel_residual: f64,
}

impl CutElasticitySolution {
    /// Displacement coefficients, two components per active node. Zero-clamped
    /// DOFs remain present as unit-identity rows.
    #[must_use]
    pub fn coefficients(&self) -> &[f64] {
        &self.coefficients
    }

    /// Nodal displacements in deterministic node-key order.
    #[must_use]
    pub fn nodal(&self) -> &BTreeMap<NodeKey, [f64; 2]> {
        &self.nodal
    }

    /// Active cells used by the solve.
    #[must_use]
    pub fn active_cells(&self) -> &[CellKey] {
        &self.active
    }

    /// Conservatively classified cut cells omitted below the documented area
    /// threshold during assembly.
    #[must_use]
    pub fn dropped_cut_cells(&self) -> usize {
        self.dropped_cut_cells
    }
}

impl CutElasticity<'_> {
    fn validate_assembly(&self) -> Result<(f64, f64), CutFemError> {
        if !self.traction_free_interface
            && !(self.nitsche_beta > 0.0 && self.nitsche_beta.is_finite())
        {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!(
                    "nitsche_beta {} must be finite and positive",
                    self.nitsche_beta
                ),
            });
        }
        if !(self.ghost_gamma >= 0.0 && self.ghost_gamma.is_finite()) {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!(
                    "ghost_gamma {} must be finite and non-negative",
                    self.ghost_gamma
                ),
            });
        }
        if self.quad_depth > 12 {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!("quad_depth {} exceeds the bounded cap 12", self.quad_depth),
            });
        }
        if !(self.material.youngs > 0.0 && self.material.youngs.is_finite()) {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!(
                    "material Young's modulus {} must be finite and positive",
                    self.material.youngs
                ),
            });
        }
        if !(self.material.poisson.is_finite()
            && self.material.poisson > -1.0
            && self.material.poisson < 0.5)
        {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!(
                    "material Poisson ratio {} must lie in (-1, 0.5)",
                    self.material.poisson
                ),
            });
        }
        if !(self.material.strain_limit > 0.0 && self.material.strain_limit.is_finite()) {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!(
                    "material strain_limit {} must be finite and positive",
                    self.material.strain_limit
                ),
            });
        }
        let (lambda, mu) = self.material.lame();
        let bulk_2d = lambda + mu;
        if !(lambda.is_finite()
            && mu > 0.0
            && mu.is_finite()
            && bulk_2d > 0.0
            && bulk_2d.is_finite())
        {
            return Err(CutFemError::InvalidElasticityInput {
                what: "material Lamé parameters do not define a finite coercive plane-strain law"
                    .to_string(),
            });
        }
        let stiffness_ratio = (lambda + 2.0 * mu) / mu;
        if !(stiffness_ratio.is_finite() && stiffness_ratio <= MAX_PLANE_STRAIN_STIFFNESS_RATIO) {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!(
                    "plane-strain stiffness ratio (lambda + 2*mu)/mu = {stiffness_ratio} exceeds the certified compressible-regime limit {MAX_PLANE_STRAIN_STIFFNESS_RATIO}; near-incompressible stabilization is not claimed"
                ),
            });
        }
        Ok((lambda, mu))
    }

    fn validate_solver(&self) -> Result<(), CutFemError> {
        if !(self.solver_tol > 0.0 && self.solver_tol.is_finite()) {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!("solver_tol {} must be finite and positive", self.solver_tol),
            });
        }
        if self.solver_max_iters == 0 {
            return Err(CutFemError::InvalidElasticityInput {
                what: "solver_max_iters must be positive".to_string(),
            });
        }
        Ok(())
    }

    /// Assemble `K u = b` for `-div sigma(u) = f`.
    ///
    /// On the embedded interface, `g` is imposed by symmetric Nitsche unless
    /// [`Self::traction_free_interface`] selects the natural boundary.
    ///
    /// # Errors
    /// Returns a structured refusal for invalid parameters/callback values,
    /// empty domains, or non-uniform active levels.
    #[allow(clippy::too_many_lines)]
    pub fn assemble(
        &self,
        f: &dyn Fn(f64, f64) -> [f64; 2],
        g: &dyn Fn(f64, f64) -> [f64; 2],
    ) -> Result<CutElasticityOperator, CutFemError> {
        let (lambda, mu) = self.validate_assembly()?;
        let mut active = Vec::new();
        let mut cut = BTreeSet::new();
        let mut rules = BTreeMap::new();
        let mut dropped_cut_cells = 0usize;
        for cell in self.grid.leaves() {
            let (lo, hi) = self.grid.rect(cell);
            let enclosure = self.sdf.enclose(lo, hi);
            if enclosure.hi() < 0.0 {
                active.push(cell);
            } else if enclosure.lo() <= 0.0 {
                let rule = cut_cell_rules(self.sdf, lo, hi, self.quad_depth);
                let area = (hi[0] - lo[0]) * (hi[1] - lo[1]);
                let inside_area: f64 = rule.bulk.iter().map(|&(_, weight)| weight).sum();
                if inside_area >= 1e-12 * area {
                    active.push(cell);
                    cut.insert(cell);
                    rules.insert(cell, rule);
                } else {
                    dropped_cut_cells += 1;
                }
            }
        }
        let Some(expected_level) = active.first().map(|cell| cell.0) else {
            return Err(CutFemError::EmptyDomain);
        };
        if let Some(&cell) = active.iter().find(|cell| cell.0 != expected_level) {
            return Err(CutFemError::ElasticityGridNotUniform {
                cell,
                expected_level,
            });
        }
        let active_set: BTreeSet<CellKey> = active.iter().copied().collect();
        let mut node_ids = BTreeMap::new();
        for &cell in &active {
            for node in self.grid.corner_nodes(cell) {
                let next = node_ids.len();
                node_ids.entry(node).or_insert(next);
            }
        }
        let ndof = 2 * node_ids.len();
        let extent = self.grid.node_extent();
        let mut clamped = vec![false; ndof];
        if let Some(predicate) = self.clamp {
            for (&node, &id) in &node_ids {
                let on_box = node.0 == 0 || node.0 == extent || node.1 == 0 || node.1 == extent;
                if on_box {
                    let point = self.grid.node_pos(node);
                    if predicate(point[0], point[1]) {
                        clamped[2 * id] = true;
                        clamped[2 * id + 1] = true;
                    }
                }
            }
        }

        let mut coo = Coo::new(ndof, ndof);
        let mut rhs = vec![0.0; ndof];
        for &cell in &active {
            let (lo, hi) = self.grid.rect(cell);
            let corners = self.grid.corner_nodes(cell);
            let ids = corners.map(|node| node_ids[&node]);
            let h = self.grid.cell_h(cell);
            let mut local_k = [[0.0; 8]; 8];
            let mut local_f = [0.0; 8];
            let full_rule;
            let bulk: &[([f64; 2], f64)] = if cut.contains(&cell) {
                &rules[&cell].bulk
            } else {
                full_rule = {
                    let mut points = Vec::with_capacity(9);
                    tensor_gauss(lo, hi, &mut points);
                    points
                };
                &full_rule
            };
            for &(point, weight) in bulk {
                let (shape, gradients) = q1(lo, hi, point);
                let body = f(point[0], point[1]);
                if body.iter().any(|value| !value.is_finite()) {
                    return Err(CutFemError::InvalidElasticityInput {
                        what: format!("body force is non-finite at {point:?}"),
                    });
                }
                for a in 0..4 {
                    for ca in 0..2 {
                        let ba = strain_row(gradients[a], ca);
                        let dba = constitutive_mul(lambda, mu, ba);
                        for b in 0..4 {
                            for cb in 0..2 {
                                let bb = strain_row(gradients[b], cb);
                                local_k[2 * a + ca][2 * b + cb] += weight * dot3(dba, bb);
                            }
                        }
                        local_f[2 * a + ca] += weight * shape[a] * body[ca];
                    }
                }
            }
            if cut.contains(&cell) && !self.traction_free_interface {
                let penalty = self.nitsche_beta * mu / h;
                if !penalty.is_finite() {
                    return Err(CutFemError::InvalidElasticityInput {
                        what: format!("derived Nitsche penalty is non-finite on cell {cell:?}"),
                    });
                }
                for &(point, weight, normal) in &rules[&cell].iface {
                    let (shape, gradients) = q1(lo, hi, point);
                    let data = g(point[0], point[1]);
                    if data.iter().any(|value| !value.is_finite()) {
                        return Err(CutFemError::InvalidElasticityInput {
                            what: format!("embedded Dirichlet data is non-finite at {point:?}"),
                        });
                    }
                    for a in 0..4 {
                        for ca in 0..2 {
                            let traction_a = shape_traction(lambda, mu, gradients[a], ca, normal);
                            for b in 0..4 {
                                for cb in 0..2 {
                                    let traction_b =
                                        shape_traction(lambda, mu, gradients[b], cb, normal);
                                    let diagonal = f64::from(ca == cb);
                                    let value = penalty * shape[a] * shape[b] * diagonal
                                        - traction_a[cb] * shape[b]
                                        - shape[a] * traction_b[ca];
                                    local_k[2 * a + ca][2 * b + cb] += weight * value;
                                }
                            }
                            local_f[2 * a + ca] += weight
                                * (penalty * shape[a] * data[ca]
                                    - traction_a[0] * data[0]
                                    - traction_a[1] * data[1]);
                        }
                    }
                }
            }
            // The weak form is symmetric analytically. Evaluate each local
            // entry independently for clarity, then canonicalize each pair to
            // one bit pattern so CG and the registered K^T VJP do not rest on
            // an "equal within roundoff" assumption.
            symmetrize_local(&mut local_k);
            scatter_local(&mut coo, &mut rhs, &clamped, &ids, &local_k, &local_f);
        }

        self.assemble_outer_traction(&active, &node_ids, &clamped, &mut rhs)?;
        for (dof, is_clamped) in clamped.iter().enumerate() {
            if *is_clamped {
                coo.push(dof, dof, 1.0);
            }
        }
        if self.ghost_gamma > 0.0 {
            let mut seen = BTreeSet::<FaceKey>::new();
            for &cell in &cut {
                for direction in 0..4u8 {
                    let Some(neighbor) = self.grid.covering_neighbor(cell, direction) else {
                        continue;
                    };
                    if !active_set.contains(&neighbor) {
                        continue;
                    }
                    if neighbor.0 != cell.0 {
                        return Err(CutFemError::CutBandNotUniform { cell, neighbor });
                    }
                    let face = if cell < neighbor {
                        (cell, neighbor)
                    } else {
                        (neighbor, cell)
                    };
                    if seen.insert(face) {
                        self.assemble_ghost_face(face, mu, &node_ids, &clamped, &mut coo)?;
                    }
                }
            }
        }
        Ok(CutElasticityOperator {
            matrix: coo.assemble(),
            rhs,
            node_ids,
            clamped,
            active,
            rules,
            dropped_cut_cells,
        })
    }

    fn assemble_outer_traction(
        &self,
        active: &[CellKey],
        node_ids: &BTreeMap<NodeKey, usize>,
        clamped: &[bool],
        rhs: &mut [f64],
    ) -> Result<(), CutFemError> {
        let Some(traction) = self.boundary_traction else {
            return Ok(());
        };
        for &cell in active {
            let (level, i, j) = cell;
            let nmax = 1u32 << level;
            let corners = self.grid.corner_nodes(cell);
            let edges = [
                (j == 0, [0usize, 1usize]),
                (i + 1 == nmax, [1, 2]),
                (j + 1 == nmax, [2, 3]),
                (i == 0, [3, 0]),
            ];
            for (on_boundary, corner_indices) in edges {
                if !on_boundary {
                    continue;
                }
                let pa = self.grid.node_pos(corners[corner_indices[0]]);
                let pb = self.grid.node_pos(corners[corner_indices[1]]);
                let dx = pb[0] - pa[0];
                let dy = pb[1] - pa[1];
                let edge_lo = [pa[0].min(pb[0]), pa[1].min(pb[1])];
                let edge_hi = [pa[0].max(pb[0]), pa[1].max(pb[1])];
                let enclosure = self.sdf.enclose(edge_lo, edge_hi);
                if !(enclosure.lo().is_finite() && enclosure.hi().is_finite()) {
                    return Err(CutFemError::InvalidElasticityInput {
                        what: format!(
                            "SDF enclosure is non-finite on loaded design-box edge {pa:?}--{pb:?}"
                        ),
                    });
                }
                if enclosure.lo() <= 0.0 && enclosure.hi() >= 0.0 {
                    return Err(CutFemError::InvalidElasticityInput {
                        what: format!(
                            "boundary traction edge {pa:?}--{pb:?} is cut by the SDF; \
                             cut-edge traction quadrature is not yet certified"
                        ),
                    });
                }
                if enclosure.lo() > 0.0 {
                    continue;
                }
                let length = dx.hypot(dy);
                let gauss = 0.5 / 3.0f64.sqrt();
                for t in [0.5 - gauss, 0.5 + gauss] {
                    let point = [pa[0] + t * dx, pa[1] + t * dy];
                    let value = traction(point[0], point[1]);
                    if value.iter().any(|component| !component.is_finite()) {
                        return Err(CutFemError::InvalidElasticityInput {
                            what: format!("boundary traction is non-finite at {point:?}"),
                        });
                    }
                    let weight = 0.5 * length;
                    for (corner_index, shape) in
                        [(corner_indices[0], 1.0 - t), (corner_indices[1], t)]
                    {
                        let id = node_ids[&corners[corner_index]];
                        for (component, traction_component) in value.iter().enumerate() {
                            let dof = 2 * id + component;
                            if !clamped[dof] {
                                rhs[dof] += weight * shape * traction_component;
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }

    fn assemble_ghost_face(
        &self,
        face: FaceKey,
        mu: f64,
        node_ids: &BTreeMap<NodeKey, usize>,
        clamped: &[bool],
        coo: &mut Coo,
    ) -> Result<(), CutFemError> {
        let (cell_a, cell_b) = face;
        let (lo_a, hi_a) = self.grid.rect(cell_a);
        let (lo_b, hi_b) = self.grid.rect(cell_b);
        let h = self.grid.cell_h(cell_a);
        let axis = usize::from(cell_a.1 == cell_b.1);
        let (t0, t1) = if axis == 0 {
            (lo_a[1], hi_a[1])
        } else {
            (lo_a[0], hi_a[0])
        };
        let normal = if axis == 0 { [1.0, 0.0] } else { [0.0, 1.0] };
        let face_coordinate = if axis == 0 { hi_a[0] } else { hi_a[1] };
        let corners_a = self.grid.corner_nodes(cell_a);
        let corners_b = self.grid.corner_nodes(cell_b);
        let gauss = 0.5 / 3.0f64.sqrt();
        let weight = 0.5 * (t1 - t0);
        let mut jump = BTreeMap::<NodeKey, [f64; 2]>::new();
        for (quadrature_index, t) in [0.5 - gauss, 0.5 + gauss].into_iter().enumerate() {
            let varying = t0 + t * (t1 - t0);
            let point = if axis == 0 {
                [face_coordinate, varying]
            } else {
                [varying, face_coordinate]
            };
            let (_, gradients_a) = q1(lo_a, hi_a, point);
            let (_, gradients_b) = q1(lo_b, hi_b, point);
            for a in 0..4 {
                jump.entry(corners_a[a]).or_default()[quadrature_index] +=
                    dot2(gradients_a[a], normal);
                jump.entry(corners_b[a]).or_default()[quadrature_index] -=
                    dot2(gradients_b[a], normal);
            }
        }
        let scale = self.ghost_gamma * mu * h * weight;
        if !scale.is_finite() {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!("derived ghost penalty is non-finite on face {face:?}"),
            });
        }
        let entries: Vec<(NodeKey, [f64; 2])> = jump.into_iter().collect();
        for (node_a, jump_a) in &entries {
            for (node_b, jump_b) in &entries {
                let value = scale * (jump_a[0] * jump_b[0] + jump_a[1] * jump_b[1]);
                if value == 0.0 {
                    continue;
                }
                for component in 0..2 {
                    let row = 2 * node_ids[node_a] + component;
                    let col = 2 * node_ids[node_b] + component;
                    if !clamped[row] && !clamped[col] {
                        coo.push(row, col, value);
                    }
                }
            }
        }
        Ok(())
    }

    /// Assemble and solve the vector problem with deterministic CG.
    ///
    /// # Errors
    /// Returns [`CutFemError::SolveNotConverged`] when the configured
    /// residual gate is missed, plus all assembly refusals.
    pub fn solve(
        &self,
        f: &dyn Fn(f64, f64) -> [f64; 2],
        g: &dyn Fn(f64, f64) -> [f64; 2],
    ) -> Result<CutElasticitySolution, CutFemError> {
        self.validate_solver()?;
        let operator = self.assemble(f, g)?;
        let preconditioner = JacobiPrecond::new(operator.matrix());
        let mut state = CgState::new(&operator, &preconditioner, operator.rhs());
        let report = state.run(
            &operator,
            &preconditioner,
            self.solver_tol,
            self.solver_max_iters,
        );
        if !report.converged {
            return Err(CutFemError::SolveNotConverged {
                iters: report.iters,
                rel_residual: report.rel_residual,
            });
        }
        let nodal = operator.nodal_values(&state.x);
        Ok(CutElasticitySolution {
            coefficients: state.x,
            nodal,
            active: operator.active,
            rules: operator.rules,
            dropped_cut_cells: operator.dropped_cut_cells,
            iters: report.iters,
            rel_residual: report.rel_residual,
        })
    }

    /// L2 and H1-seminorm displacement errors, integrated with one deeper
    /// cut rule than the solve.
    #[must_use]
    pub fn l2_h1_error(
        &self,
        solution: &CutElasticitySolution,
        exact: &dyn Fn(f64, f64) -> [f64; 2],
        exact_gradient: &dyn Fn(f64, f64) -> [[f64; 2]; 2],
    ) -> (f64, f64) {
        let mut l2 = 0.0f64;
        let mut h1 = 0.0f64;
        for &cell in &solution.active {
            let (lo, hi) = self.grid.rect(cell);
            let corners = self.grid.corner_nodes(cell);
            let values = corners.map(|node| solution.nodal[&node]);
            let quadrature;
            let rule: &[([f64; 2], f64)] = if solution.rules.contains_key(&cell) {
                quadrature = cut_cell_rules(self.sdf, lo, hi, self.quad_depth + 1).bulk;
                &quadrature
            } else {
                quadrature = {
                    let mut points = Vec::with_capacity(9);
                    tensor_gauss(lo, hi, &mut points);
                    points
                };
                &quadrature
            };
            for &(point, weight) in rule {
                let (shape, gradients) = q1(lo, hi, point);
                let mut computed = [0.0; 2];
                let mut computed_gradient = [[0.0; 2]; 2];
                for a in 0..4 {
                    for component in 0..2 {
                        computed[component] += shape[a] * values[a][component];
                        computed_gradient[component][0] += gradients[a][0] * values[a][component];
                        computed_gradient[component][1] += gradients[a][1] * values[a][component];
                    }
                }
                let expected = exact(point[0], point[1]);
                let expected_gradient = exact_gradient(point[0], point[1]);
                for component in 0..2 {
                    let value_error = expected[component] - computed[component];
                    l2 += weight * value_error * value_error;
                    for axis in 0..2 {
                        let gradient_error =
                            expected_gradient[component][axis] - computed_gradient[component][axis];
                        h1 += weight * gradient_error * gradient_error;
                    }
                }
            }
        }
        // Every contribution above is non-negative.  Do not clamp here:
        // `f64::max` treats a single NaN as the other operand, which could
        // turn a poisoned error integral into a false zero-error certificate.
        // A negative or non-finite accumulator must remain non-finite so the
        // acceptance gates fail closed.
        (l2.sqrt(), h1.sqrt())
    }
}

fn symmetrize_local(matrix: &mut [[f64; 8]; 8]) {
    for row in 0..8 {
        for column in (row + 1)..8 {
            let value = f64::midpoint(matrix[row][column], matrix[column][row]);
            matrix[row][column] = value;
            matrix[column][row] = value;
        }
    }
}

fn scatter_local(
    coo: &mut Coo,
    rhs: &mut [f64],
    clamped: &[bool],
    node_ids: &[usize; 4],
    local_k: &[[f64; 8]; 8],
    local_f: &[f64; 8],
) {
    for a in 0..8 {
        let row = 2 * node_ids[a / 2] + a % 2;
        if clamped[row] {
            continue;
        }
        rhs[row] += local_f[a];
        for b in 0..8 {
            let col = 2 * node_ids[b / 2] + b % 2;
            let value = local_k[a][b];
            if value != 0.0 && !clamped[col] {
                coo.push(row, col, value);
            }
        }
    }
}

fn strain_row(gradient: [f64; 2], component: usize) -> [f64; 3] {
    if component == 0 {
        [gradient[0], 0.0, gradient[1]]
    } else {
        [0.0, gradient[1], gradient[0]]
    }
}

fn constitutive_mul(lambda: f64, mu: f64, strain: [f64; 3]) -> [f64; 3] {
    [
        (lambda + 2.0 * mu) * strain[0] + lambda * strain[1],
        lambda * strain[0] + (lambda + 2.0 * mu) * strain[1],
        mu * strain[2],
    ]
}

fn shape_traction(
    lambda: f64,
    mu: f64,
    gradient: [f64; 2],
    component: usize,
    normal: [f64; 2],
) -> [f64; 2] {
    let gradient_normal = dot2(gradient, normal);
    let mut traction = [0.0; 2];
    for axis in 0..2 {
        traction[axis] =
            lambda * gradient[component] * normal[axis] + mu * normal[component] * gradient[axis];
    }
    traction[component] += mu * gradient_normal;
    traction
}

fn dot2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[cfg(feature = "adjoint-vjp")]
#[derive(Clone)]
struct ElasticityApplyVjp {
    matrix: Csr,
}

#[cfg(feature = "adjoint-vjp")]
impl fs_adjoint::transpose::Vjp for ElasticityApplyVjp {
    fn vjp(&self, primal_inputs: &[&[f64]], out_cotangent: &[f64]) -> Vec<Vec<f64>> {
        assert_eq!(primal_inputs.len(), 1, "elasticity apply VJP arity");
        assert_eq!(
            primal_inputs[0].len(),
            self.matrix.ncols(),
            "elasticity primal length"
        );
        assert_eq!(
            out_cotangent.len(),
            self.matrix.nrows(),
            "elasticity cotangent length"
        );
        let mut input_cotangent = vec![0.0; self.matrix.ncols()];
        // The symmetric Nitsche + ghost form has K^T = K. Keeping the
        // explicit VJP object makes that fact testable at the registry seam.
        self.matrix.spmv(out_cotangent, &mut input_cotangent);
        vec![input_cotangent]
    }
}

/// Content-addressed registry key for this exact discrete operator.
///
/// Distinct matrices must not share a registry entry: the registry is keyed by
/// strings and deliberately replaces duplicate keys. Hashing the canonical CSR
/// shape, sparsity, and value bits makes multi-operator tapes deterministic and
/// prevents a later registration from silently changing an earlier node's VJP.
#[cfg(feature = "adjoint-vjp")]
#[must_use]
pub fn elasticity_apply_vjp_key(operator: &CutElasticityOperator) -> String {
    const DOMAIN: &str = "frankensim/fs-cutfem/elasticity-apply-key/v1";
    let matrix = operator.matrix();
    let mut payload = Vec::with_capacity(matrix.nnz().saturating_mul(16));
    push_usize(&mut payload, matrix.nrows());
    push_usize(&mut payload, matrix.ncols());
    for row in 0..matrix.nrows() {
        let (columns, values) = matrix.row(row);
        push_usize(&mut payload, columns.len());
        for (&column, &value) in columns.iter().zip(values) {
            push_usize(&mut payload, column);
            payload.extend_from_slice(&value.to_bits().to_le_bytes());
        }
    }
    format!(
        "{ELASTICITY_APPLY_VJP_OP}:{}",
        fs_blake3::hash_domain(DOMAIN, &payload)
    )
}

#[cfg(feature = "adjoint-vjp")]
fn push_usize(payload: &mut Vec<u8>, value: usize) {
    let encoded = u64::try_from(value).expect("CSR index must fit the portable u64 key encoding");
    payload.extend_from_slice(&encoded.to_le_bytes());
}

/// Register this exact vector-elasticity apply with fs-adjoint's ledger-DAG
/// VJP registry and return the content-addressed op key to record on the tape.
#[cfg(feature = "adjoint-vjp")]
#[must_use = "record the returned content-addressed key on the adjoint tape"]
pub fn register_elasticity_apply_vjp(
    registry: &mut fs_adjoint::transpose::VjpRegistry,
    operator: &CutElasticityOperator,
) -> String {
    use std::sync::Arc;
    let key = elasticity_apply_vjp_key(operator);
    registry.register(
        &key,
        Arc::new(ElasticityApplyVjp {
            matrix: operator.matrix.clone(),
        }),
    );
    key
}
