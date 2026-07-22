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
//! On graded active trees, scalar 2:1 hanging-node constraints are lifted
//! componentwise and eliminated by the explicit transform `T`: element and
//! face terms scatter as `T^T K T`, loads as `T^T f`, and solved terminal
//! coefficients reconstruct every active mesh node. Uniform trees retain the
//! literal legacy scatter path so their operator and topology bits do not move.

use crate::CutFemError;
use crate::fem::{JacobiPrecond, q1, recomputed_euclidean_residual_claim};
use crate::grid::{CellKey, FaceDirection, NodeKey, Quadtree, SharedFacePatch};
use crate::quad::{CutRules, cut_cell_rules, tensor_gauss};
use crate::sdf::CutSdf;
use fs_material::IsotropicElastic;
use fs_solver::krylov::{CgState, ResidualClaim};
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
type NodeExpansion = Vec<(usize, f64)>;
type TerminalIds = BTreeMap<NodeKey, usize>;
type NodeExpansions = BTreeMap<NodeKey, NodeExpansion>;

/// One named edge of the unit-square CutFEM design box.
///
/// [`EdgeBand`] coordinates increase with the corresponding Cartesian
/// coordinate: `x` on [`Self::Bottom`] and [`Self::Top`], and `y` on
/// [`Self::Left`] and [`Self::Right`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DesignBoxEdge {
    /// The edge `y = 0`, parameterized by increasing `x`.
    Bottom,
    /// The edge `x = 1`, parameterized by increasing `y`.
    Right,
    /// The edge `y = 1`, parameterized by increasing `x`.
    Top,
    /// The edge `x = 0`, parameterized by increasing `y`.
    Left,
}

/// Checked closed support band on one named design-box edge.
///
/// Both endpoints are normalized design-box coordinates in `[0, 1]`. The
/// fields are private so an invalid, non-finite, or reversed band cannot enter
/// certified boundary assembly.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EdgeBand {
    edge: DesignBoxEdge,
    start: f64,
    end: f64,
}

impl EdgeBand {
    /// Construct a closed support band `start..=end` on `edge`.
    ///
    /// A zero-length band is valid and represents measure-zero support. Its
    /// endpoint is still classified rather than assumed to carry zero data.
    ///
    /// # Errors
    /// Refuses non-finite endpoints, endpoints outside `[0, 1]`, or a reversed
    /// interval.
    pub fn new(edge: DesignBoxEdge, start: f64, end: f64) -> Result<Self, CutFemError> {
        if !(start.is_finite() && end.is_finite()) {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!(
                    "boundary traction edge-band endpoints {start}..={end} must be finite"
                ),
            });
        }
        if !(0.0..=1.0).contains(&start) || !(0.0..=1.0).contains(&end) {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!(
                    "boundary traction edge-band endpoints {start}..={end} must lie in [0, 1]"
                ),
            });
        }
        if start > end {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!("boundary traction edge-band start {start} exceeds end {end}"),
            });
        }
        Ok(Self { edge, start, end })
    }

    /// Named design-box edge carrying this support.
    #[must_use]
    pub fn edge(self) -> DesignBoxEdge {
        self.edge
    }

    /// Inclusive normalized start coordinate.
    #[must_use]
    pub fn start(self) -> f64 {
        self.start
    }

    /// Inclusive normalized end coordinate.
    #[must_use]
    pub fn end(self) -> f64 {
        self.end
    }
}

/// Boundary-traction data supplied explicitly to typed assembly.
#[derive(Clone, Copy)]
pub enum BoundaryTraction<'a> {
    /// A legacy callback with uncertified, potentially nonzero support on every
    /// design-box edge. Every active box-edge segment is therefore classified.
    Uncertified(&'a dyn Fn(f64, f64) -> [f64; 2]),
    /// A callback whose applied value is defined to be zero outside `support`.
    /// Only the exact intersection with that named band is classified and
    /// integrated.
    EdgeBand {
        /// Checked closed support on one named design-box edge.
        support: EdgeBand,
        /// Traction value evaluated only on the supported subsegment.
        value: &'a dyn Fn(f64, f64) -> [f64; 2],
    },
}

#[derive(Clone, Copy)]
struct SupportedTractionSegment<'a> {
    value: &'a dyn Fn(f64, f64) -> [f64; 2],
    local_start: f64,
    local_end: f64,
    start: [f64; 2],
    end: [f64; 2],
}

fn supported_traction_segment(
    boundary_traction: BoundaryTraction<'_>,
    edge: DesignBoxEdge,
    point_a: [f64; 2],
    point_b: [f64; 2],
) -> Option<SupportedTractionSegment<'_>> {
    match boundary_traction {
        BoundaryTraction::Uncertified(value) => Some(SupportedTractionSegment {
            value,
            local_start: 0.0,
            local_end: 1.0,
            start: point_a,
            end: point_b,
        }),
        BoundaryTraction::EdgeBand { support, value } => {
            if support.edge != edge {
                return None;
            }
            let (coordinate_a, coordinate_b) = match edge {
                DesignBoxEdge::Bottom | DesignBoxEdge::Top => (point_a[0], point_b[0]),
                DesignBoxEdge::Left | DesignBoxEdge::Right => (point_a[1], point_b[1]),
            };
            let edge_start = coordinate_a.min(coordinate_b);
            let edge_end = coordinate_a.max(coordinate_b);
            // Closed intervals deliberately use strict disjointness: endpoint
            // contact remains potentially nonzero support and is classified.
            if support.end < edge_start || support.start > edge_end {
                return None;
            }
            let supported_start = support.start.max(edge_start);
            let supported_end = support.end.min(edge_end);
            let first = (supported_start - coordinate_a) / (coordinate_b - coordinate_a);
            let second = (supported_end - coordinate_a) / (coordinate_b - coordinate_a);
            let (start, end) = match edge {
                DesignBoxEdge::Bottom | DesignBoxEdge::Top => {
                    ([supported_start, point_a[1]], [supported_end, point_a[1]])
                }
                DesignBoxEdge::Left | DesignBoxEdge::Right => {
                    ([point_a[0], supported_start], [point_a[0], supported_end])
                }
            };
            Some(SupportedTractionSegment {
                value,
                local_start: first.min(second),
                local_end: first.max(second),
                start,
                end,
            })
        }
    }
}

/// Vector Q1 CutFEM problem on `Omega = {phi < 0}`.
///
/// The constitutive parameters come from [`IsotropicElastic`], so the
/// material's admissibility checks and model-card identity are shared
/// with the rest of FLUX rather than duplicated here.
#[derive(Clone, Copy)]
pub struct CutElasticity<'a> {
    /// Uniform or 2:1-balanced graded background quadtree.
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
    /// Optional legacy dead traction with uncertified, potentially nonzero
    /// support on every active design-box boundary edge. Prefer
    /// [`Self::assemble_with_boundary_traction`] for checked named support.
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
    node_ids: TerminalIds,
    /// Present only when hanging nodes exist. Keys cover every active mesh
    /// node; values resolve it into unconstrained terminal blocks.
    node_expansions: Option<NodeExpansions>,
    clamped: Vec<bool>,
    active: Vec<CellKey>,
    rules: BTreeMap<CellKey, CutRules>,
    ghost_faces: Vec<FaceKey>,
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

    /// Deterministic terminal-node-to-block map. Node `id` owns displacement
    /// DOFs `2*id` and `2*id + 1`. On a uniform tree every active node is a
    /// terminal, preserving the original map exactly.
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

    /// Active cells in canonical quadtree-key order.
    #[must_use]
    pub fn active_cells(&self) -> &[CellKey] {
        &self.active
    }

    /// Certified cut-cell quadrature retained by the assembly.
    #[must_use]
    pub fn cut_rules(&self) -> &BTreeMap<CellKey, CutRules> {
        &self.rules
    }

    /// Canonically ordered leaf pairs whose exact shared patches carry ghost
    /// stabilization. One coarse face can therefore contribute two pairs,
    /// one for each fine neighbor. The slice is empty when
    /// `ghost_gamma == 0`.
    #[must_use]
    pub fn ghost_faces(&self) -> &[(CellKey, CellKey)] {
        &self.ghost_faces
    }

    /// Algebraic compliance `b^T x` for this assembled load vector.
    ///
    /// This is the discrete assembled-load functional. In particular, when
    /// nonzero embedded displacement data contribute to `b`, it is not by
    /// itself a claim about physical external work.
    ///
    /// # Errors
    /// Refuses a coefficient-length mismatch or any non-finite input/result.
    pub fn algebraic_compliance(&self, x: &[f64]) -> Result<f64, CutFemError> {
        if x.len() != self.dof_count() {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!(
                    "elasticity coefficient length {} does not match DOF count {}",
                    x.len(),
                    self.dof_count()
                ),
            });
        }
        if x.iter().any(|value| !value.is_finite()) {
            return Err(CutFemError::InvalidElasticityInput {
                what: "elasticity coefficients must be finite".to_string(),
            });
        }
        let compliance: f64 = self
            .rhs
            .iter()
            .zip(x)
            .map(|(load, coefficient)| load * coefficient)
            .sum();
        if !compliance.is_finite() {
            return Err(CutFemError::InvalidElasticityInput {
                what: "algebraic compliance b^T x is non-finite".to_string(),
            });
        }
        Ok(compliance)
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
        let Some(expansions) = &self.node_expansions else {
            return self
                .node_ids
                .iter()
                .map(|(&node, &id)| (node, [x[2 * id], x[2 * id + 1]]))
                .collect();
        };
        expansions
            .iter()
            .map(|(&node, expansion)| {
                let mut value = [0.0; 2];
                for &(id, weight) in expansion {
                    value[0] += weight * x[2 * id];
                    value[1] += weight * x[2 * id + 1];
                }
                (node, value)
            })
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
    ghost_faces: Vec<FaceKey>,
    compliance: f64,
    dropped_cut_cells: usize,
    /// CG iterations.
    pub iters: usize,
    /// Final recomputed Euclidean relative residual.
    pub rel_residual: f64,
    residual_claim: ResidualClaim,
}

impl CutElasticitySolution {
    /// Typed provenance for [`Self::rel_residual`].
    #[must_use]
    pub fn residual_claim(&self) -> ResidualClaim {
        self.residual_claim
    }

    /// The recomputed Euclidean relative residual.
    ///
    /// This accessor refuses any future non-Euclidean residual claim rather
    /// than letting a recurrence estimate enter a certificate as a bound.
    #[must_use]
    pub fn euclidean_rel_residual(&self) -> Option<f64> {
        self.residual_claim.euclidean()
    }

    /// Displacement coefficients, two components per algebraic terminal node.
    /// Zero-clamped DOFs remain present as unit-identity rows; hanging values
    /// are available through [`Self::nodal`].
    #[must_use]
    pub fn coefficients(&self) -> &[f64] {
        &self.coefficients
    }

    /// Number of vector displacement DOFs, including zero-clamped identity
    /// rows retained in the coefficient vector.
    #[must_use]
    pub fn dof_count(&self) -> usize {
        self.coefficients.len()
    }

    /// Exact deterministic dot product `b^T x` for the assembled discrete
    /// load and converged coefficient vector.
    #[must_use]
    pub fn compliance(&self) -> f64 {
        self.compliance
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

    /// Certified cut-cell quadrature retained by the solve's assembly.
    #[must_use]
    pub fn cut_rules(&self) -> &BTreeMap<CellKey, CutRules> {
        &self.rules
    }

    /// Canonically ordered leaf pairs whose exact shared patches carry ghost
    /// stabilization. The slice is empty when `ghost_gamma == 0`.
    #[must_use]
    pub fn ghost_faces(&self) -> &[(CellKey, CellKey)] {
        &self.ghost_faces
    }

    /// Evaluate the vector Q1 displacement and its componentwise gradient on
    /// one active cell without inventing values for absent or poisoned nodes.
    ///
    /// The gradient is indexed `[component][axis]`.
    ///
    /// # Errors
    /// Refuses inactive cells, points outside the cell, missing corners, or
    /// any non-finite point, nodal value, shape quantity, or result.
    pub fn value_gradient(
        &self,
        grid: &Quadtree,
        cell: CellKey,
        point: [f64; 2],
    ) -> Result<([f64; 2], [[f64; 2]; 2]), CutFemError> {
        if self.active.binary_search(&cell).is_err() {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!("cannot evaluate inactive elasticity cell {cell:?}"),
            });
        }
        if point.iter().any(|coordinate| !coordinate.is_finite()) {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!("elasticity evaluation point {point:?} must be finite"),
            });
        }
        let (lo, hi) = grid.rect(cell);
        if (0..2).any(|axis| point[axis] < lo[axis] || point[axis] > hi[axis]) {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!(
                    "elasticity evaluation point {point:?} lies outside cell {cell:?} rectangle {lo:?}--{hi:?}"
                ),
            });
        }
        let corners = grid.corner_nodes(cell);
        let mut values = [[0.0; 2]; 4];
        for (corner, value) in corners.iter().zip(&mut values) {
            let Some(nodal) = self.nodal.get(corner) else {
                return Err(CutFemError::InvalidElasticityInput {
                    what: format!("elasticity cell {cell:?} is missing nodal corner {corner:?}"),
                });
            };
            if nodal.iter().any(|component| !component.is_finite()) {
                return Err(CutFemError::InvalidElasticityInput {
                    what: format!("elasticity nodal value at corner {corner:?} is non-finite"),
                });
            }
            *value = *nodal;
        }
        let (shapes, gradients) = q1(lo, hi, point);
        if shapes.iter().any(|shape| !shape.is_finite())
            || gradients
                .iter()
                .flatten()
                .any(|gradient| !gradient.is_finite())
        {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!("Q1 evaluation is non-finite on cell {cell:?} at {point:?}"),
            });
        }
        let mut displacement = [0.0; 2];
        let mut gradient = [[0.0; 2]; 2];
        for corner in 0..4 {
            for component in 0..2 {
                displacement[component] += shapes[corner] * values[corner][component];
                for axis in 0..2 {
                    gradient[component][axis] +=
                        gradients[corner][axis] * values[corner][component];
                }
            }
        }
        if displacement
            .iter()
            .chain(gradient.iter().flatten())
            .any(|value| !value.is_finite())
        {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!("elasticity field evaluation is non-finite on cell {cell:?}"),
            });
        }
        Ok((displacement, gradient))
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
        if !(self.traction_free_interface
            || (self.nitsche_beta > 0.0 && self.nitsche_beta.is_finite()))
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

    fn require_explicit_boundary_traction_slot(&self) -> Result<(), CutFemError> {
        if self.boundary_traction.is_some() {
            return Err(CutFemError::InvalidElasticityInput {
                what: "assemble_with_boundary_traction/solve_with_boundary_traction require CutElasticity::boundary_traction to be None"
                    .to_string(),
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
    /// empty domains, or malformed hanging-node transforms.
    #[allow(clippy::too_many_lines)]
    pub fn assemble(
        &self,
        f: &dyn Fn(f64, f64) -> [f64; 2],
        g: &dyn Fn(f64, f64) -> [f64; 2],
    ) -> Result<CutElasticityOperator, CutFemError> {
        self.assemble_impl(
            f,
            g,
            self.boundary_traction.map(BoundaryTraction::Uncertified),
        )
    }

    /// Assemble with an explicit boundary-traction support descriptor.
    ///
    /// [`BoundaryTraction::Uncertified`] preserves the full-edge fail-closed
    /// callback semantics of [`Self::boundary_traction`]. A checked
    /// [`BoundaryTraction::EdgeBand`] instead defines the applied traction as
    /// zero outside its named support, so only the exact supported subsegment
    /// is classified and integrated.
    ///
    /// # Errors
    /// In addition to the ordinary assembly refusals, this method requires
    /// [`Self::boundary_traction`] to be `None`; two competing sources are
    /// never merged implicitly.
    pub fn assemble_with_boundary_traction(
        &self,
        f: &dyn Fn(f64, f64) -> [f64; 2],
        g: &dyn Fn(f64, f64) -> [f64; 2],
        boundary_traction: BoundaryTraction<'_>,
    ) -> Result<CutElasticityOperator, CutFemError> {
        self.require_explicit_boundary_traction_slot()?;
        self.assemble_impl(f, g, Some(boundary_traction))
    }

    #[allow(clippy::too_many_lines)]
    fn assemble_impl(
        &self,
        f: &dyn Fn(f64, f64) -> [f64; 2],
        g: &dyn Fn(f64, f64) -> [f64; 2],
        boundary_traction: Option<BoundaryTraction<'_>>,
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
        if active.is_empty() {
            return Err(CutFemError::EmptyDomain);
        }
        let active_set: BTreeSet<CellKey> = active.iter().copied().collect();
        let mut nodes = BTreeSet::new();
        let mut legacy_node_ids = BTreeMap::new();
        for &cell in &active {
            for node in self.grid.corner_nodes(cell) {
                nodes.insert(node);
                let next = legacy_node_ids.len();
                legacy_node_ids.entry(node).or_insert(next);
            }
        }
        let constraints: BTreeMap<NodeKey, Vec<(NodeKey, f64)>> = self
            .grid
            .hanging_constraints(&active_set, &nodes)
            .into_iter()
            .map(|(node, terms)| (node, terms.to_vec()))
            .collect();
        let (node_ids, node_expansions) = if constraints.is_empty() {
            // Do not route a uniform operator through even an identity
            // transform: this is the literal pre-graded numbering/scatter
            // path, preserving CSR/RHS/map/clamp/topology bits.
            (legacy_node_ids, None)
        } else {
            let (terminal_ids, expansions) = build_terminal_expansions(&nodes, &constraints)?;
            (terminal_ids, Some(expansions))
        };
        let ndof = 2 * node_ids.len();
        let clamped = self.build_clamp_mask(&node_ids, node_expansions.as_ref())?;

        let mut coo = Coo::new(ndof, ndof);
        let mut rhs = vec![0.0; ndof];
        for &cell in &active {
            let (lo, hi) = self.grid.rect(cell);
            let corners = self.grid.corner_nodes(cell);
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
            if let Some(expansions) = &node_expansions {
                scatter_local_reduced(
                    &mut coo, &mut rhs, &clamped, &corners, expansions, &local_k, &local_f,
                )?;
            } else {
                let ids = corners.map(|node| node_ids[&node]);
                scatter_local(&mut coo, &mut rhs, &clamped, &ids, &local_k, &local_f);
            }
        }

        self.assemble_outer_traction(boundary_traction, &active, &node_ids, &clamped, &mut rhs)?;
        for (dof, is_clamped) in clamped.iter().enumerate() {
            if *is_clamped {
                coo.push(dof, dof, 1.0);
            }
        }
        let mut ghost_faces = BTreeSet::<FaceKey>::new();
        if self.ghost_gamma > 0.0 {
            for &cell in &cut {
                for direction in FaceDirection::ALL {
                    let patches = self.grid.face_neighbors(cell, direction).map_err(|error| {
                        CutFemError::InvalidElasticityInput {
                            what: format!(
                                "cannot derive 2:1 ghost-face topology from cut cell {cell:?}: {error}"
                            ),
                        }
                    })?;
                    for patch in patches {
                        let neighbor = patch.other(cell).ok_or_else(|| {
                            CutFemError::InvalidElasticityInput {
                                what: format!(
                                    "shared ghost patch {patch:?} does not contain source cell {cell:?}"
                                ),
                            }
                        })?;
                        if !active_set.contains(&neighbor) {
                            continue;
                        }
                        let face = patch.canonical_cells();
                        if ghost_faces.insert(face) {
                            self.assemble_ghost_face(
                                patch,
                                mu,
                                &node_ids,
                                node_expansions.as_ref(),
                                &clamped,
                                &mut coo,
                            )?;
                        }
                    }
                }
            }
        }
        Ok(CutElasticityOperator {
            matrix: coo.assemble(),
            rhs,
            node_ids,
            node_expansions,
            clamped,
            active,
            rules,
            ghost_faces: ghost_faces.into_iter().collect(),
            dropped_cut_cells,
        })
    }

    fn build_clamp_mask(
        &self,
        node_ids: &TerminalIds,
        expansions: Option<&NodeExpansions>,
    ) -> Result<Vec<bool>, CutFemError> {
        let mut selected = BTreeSet::new();
        let Some(predicate) = self.clamp else {
            return Ok(vec![false; 2 * node_ids.len()]);
        };
        let extent = self.grid.node_extent();
        let mut select = |node: NodeKey| {
            let on_box = node.0 == 0 || node.0 == extent || node.1 == 0 || node.1 == extent;
            if on_box {
                let point = self.grid.node_pos(node);
                if predicate(point[0], point[1]) {
                    selected.insert(node);
                }
            }
        };
        if let Some(expansions) = expansions {
            for node in expansions.keys().copied() {
                select(node);
            }
        } else {
            for node in node_ids.keys().copied() {
                select(node);
            }
        }
        clamp_terminal_blocks(node_ids, expansions, &selected)
    }

    fn assemble_outer_traction(
        &self,
        boundary_traction: Option<BoundaryTraction<'_>>,
        active: &[CellKey],
        node_ids: &TerminalIds,
        clamped: &[bool],
        rhs: &mut [f64],
    ) -> Result<(), CutFemError> {
        let Some(boundary_traction) = boundary_traction else {
            return Ok(());
        };
        for &cell in active {
            let (level, i, j) = cell;
            let nmax = 1u32 << level;
            let corners = self.grid.corner_nodes(cell);
            let edges = [
                (j == 0, DesignBoxEdge::Bottom, [0usize, 1usize]),
                (i + 1 == nmax, DesignBoxEdge::Right, [1, 2]),
                (j + 1 == nmax, DesignBoxEdge::Top, [2, 3]),
                (i == 0, DesignBoxEdge::Left, [3, 0]),
            ];
            for (on_boundary, edge, corner_indices) in edges {
                if !on_boundary {
                    continue;
                }
                let pa = self.grid.node_pos(corners[corner_indices[0]]);
                let pb = self.grid.node_pos(corners[corner_indices[1]]);
                let dx = pb[0] - pa[0];
                let dy = pb[1] - pa[1];
                let Some(segment) = supported_traction_segment(boundary_traction, edge, pa, pb)
                else {
                    continue;
                };
                let segment_lo = [
                    segment.start[0].min(segment.end[0]),
                    segment.start[1].min(segment.end[1]),
                ];
                let segment_hi = [
                    segment.start[0].max(segment.end[0]),
                    segment.start[1].max(segment.end[1]),
                ];
                let enclosure = self.sdf.enclose(segment_lo, segment_hi);
                if !(enclosure.lo().is_finite() && enclosure.hi().is_finite()) {
                    return Err(CutFemError::InvalidElasticityInput {
                        what: format!(
                            "SDF enclosure is non-finite on supported {edge:?} design-box segment {:?}--{:?}",
                            segment.start, segment.end
                        ),
                    });
                }
                if enclosure.lo() <= 0.0 && enclosure.hi() >= 0.0 {
                    return Err(CutFemError::InvalidElasticityInput {
                        what: format!(
                            "supported {edge:?} boundary traction segment {:?}--{:?} is cut by the SDF; \
                             supported cut-edge traction quadrature is not yet certified",
                            segment.start, segment.end
                        ),
                    });
                }
                if enclosure.lo() > 0.0 {
                    continue;
                }
                let length = dx.hypot(dy);
                let midpoint = f64::midpoint(segment.local_start, segment.local_end);
                let half_span = 0.5 * (segment.local_end - segment.local_start);
                let gauss = half_span / 3.0f64.sqrt();
                for t in [midpoint - gauss, midpoint + gauss] {
                    let point = [pa[0] + t * dx, pa[1] + t * dy];
                    let value = (segment.value)(point[0], point[1]);
                    if value.iter().any(|component| !component.is_finite()) {
                        return Err(CutFemError::InvalidElasticityInput {
                            what: format!("boundary traction is non-finite at {point:?}"),
                        });
                    }
                    let weight = half_span * length;
                    for (corner_index, shape) in
                        [(corner_indices[0], 1.0 - t), (corner_indices[1], t)]
                    {
                        let node = corners[corner_index];
                        let id = node_ids.get(&node).copied().ok_or_else(|| {
                            CutFemError::InvalidElasticityInput {
                                what: format!(
                                    "design-box boundary traction node {node:?} is not an unconstrained terminal"
                                ),
                            }
                        })?;
                        for (component, traction_component) in value.iter().enumerate() {
                            let load = weight * shape * traction_component;
                            let dof = 2 * id + component;
                            if !clamped[dof] {
                                rhs[dof] += load;
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
        patch: SharedFacePatch,
        mu: f64,
        node_ids: &TerminalIds,
        expansions: Option<&NodeExpansions>,
        clamped: &[bool],
        coo: &mut Coo,
    ) -> Result<(), CutFemError> {
        let (cell_a, cell_b) = patch.oriented_cells();
        let (lo_a, hi_a) = self.grid.rect(cell_a);
        let (lo_b, hi_b) = self.grid.rect(cell_b);
        let h = patch.h_f();
        let axis = patch.axis().index();
        let (t0, t1) = patch.tangent_interval();
        let normal = patch.axis().normal();
        let face_coordinate = patch.coordinate();
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
                what: format!("derived ghost penalty is non-finite on patch {patch:?}"),
            });
        }
        if let Some(expansions) = expansions {
            let mut terminal_jump = BTreeMap::<usize, [f64; 2]>::new();
            for (node, node_jump) in jump {
                for &(id, terminal_weight) in &expansions[&node] {
                    let entry = terminal_jump.entry(id).or_default();
                    entry[0] += terminal_weight * node_jump[0];
                    entry[1] += terminal_weight * node_jump[1];
                }
            }
            let entries: Vec<(usize, [f64; 2])> = terminal_jump.into_iter().collect();
            for (id_a, jump_a) in &entries {
                for (id_b, jump_b) in &entries {
                    let value = scale * (jump_a[0] * jump_b[0] + jump_a[1] * jump_b[1]);
                    if value == 0.0 {
                        continue;
                    }
                    for component in 0..2 {
                        let row = 2 * *id_a + component;
                        let col = 2 * *id_b + component;
                        if !clamped[row] && !clamped[col] {
                            coo.push(row, col, value);
                        }
                    }
                }
            }
            return Ok(());
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
        self.solve_impl(
            f,
            g,
            self.boundary_traction.map(BoundaryTraction::Uncertified),
        )
    }

    /// Assemble and solve with an explicit boundary-traction support
    /// descriptor.
    ///
    /// # Errors
    /// In addition to the ordinary solve and assembly refusals, this method
    /// requires [`Self::boundary_traction`] to be `None`.
    pub fn solve_with_boundary_traction(
        &self,
        f: &dyn Fn(f64, f64) -> [f64; 2],
        g: &dyn Fn(f64, f64) -> [f64; 2],
        boundary_traction: BoundaryTraction<'_>,
    ) -> Result<CutElasticitySolution, CutFemError> {
        self.require_explicit_boundary_traction_slot()?;
        self.solve_impl(f, g, Some(boundary_traction))
    }

    fn solve_impl(
        &self,
        f: &dyn Fn(f64, f64) -> [f64; 2],
        g: &dyn Fn(f64, f64) -> [f64; 2],
        boundary_traction: Option<BoundaryTraction<'_>>,
    ) -> Result<CutElasticitySolution, CutFemError> {
        self.validate_solver()?;
        let operator = self.assemble_impl(f, g, boundary_traction)?;
        let preconditioner = JacobiPrecond::new(operator.matrix());
        let mut state = CgState::new(&operator, &preconditioner, operator.rhs());
        let _ = state.run(
            &operator,
            &preconditioner,
            self.solver_tol,
            self.solver_max_iters,
        );
        let residual_claim =
            recomputed_euclidean_residual_claim(&operator, &state.x, operator.rhs());
        let rel_residual = residual_claim
            .euclidean()
            .expect("the elasticity solve stores an explicitly recomputed Euclidean residual");
        if !rel_residual.is_finite() || rel_residual >= self.solver_tol {
            return Err(CutFemError::SolveNotConverged {
                iters: state.iters,
                rel_residual,
            });
        }
        let nodal = operator.nodal_values(&state.x);
        let compliance = operator.algebraic_compliance(&state.x)?;
        Ok(CutElasticitySolution {
            coefficients: state.x,
            nodal,
            active: operator.active,
            rules: operator.rules,
            ghost_faces: operator.ghost_faces,
            compliance,
            dropped_cut_cells: operator.dropped_cut_cells,
            iters: state.iters,
            rel_residual,
            residual_claim,
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

fn build_terminal_expansions(
    nodes: &BTreeSet<NodeKey>,
    constraints: &BTreeMap<NodeKey, Vec<(NodeKey, f64)>>,
) -> Result<(TerminalIds, NodeExpansions), CutFemError> {
    for (&node, terms) in constraints {
        if !nodes.contains(&node) {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!("hanging constraint target {node:?} is not an active mesh node"),
            });
        }
        validate_affine_terms(node, terms, "raw")?;
        if let Some((child, _)) = terms.iter().find(|(child, _)| !nodes.contains(child)) {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!(
                    "hanging constraint target {node:?} references absent mesh node {child:?}"
                ),
            });
        }
    }
    let node_ids: TerminalIds = nodes
        .iter()
        .filter(|node| !constraints.contains_key(node))
        .copied()
        .enumerate()
        .map(|(id, node)| (node, id))
        .collect();
    if node_ids.is_empty() {
        return Err(CutFemError::InvalidElasticityInput {
            what: "hanging constraint graph has no unconstrained terminal node".to_string(),
        });
    }
    let mut memo = BTreeMap::new();
    for &node in nodes {
        let mut stack = BTreeSet::new();
        expand_terminal_node(node, constraints, &node_ids, &mut memo, &mut stack)?;
    }
    Ok((node_ids, memo))
}

fn expand_terminal_node(
    node: NodeKey,
    constraints: &BTreeMap<NodeKey, Vec<(NodeKey, f64)>>,
    node_ids: &TerminalIds,
    memo: &mut NodeExpansions,
    stack: &mut BTreeSet<NodeKey>,
) -> Result<(), CutFemError> {
    if memo.contains_key(&node) {
        return Ok(());
    }
    if !stack.insert(node) {
        return Err(CutFemError::ConstraintCycle { node });
    }
    let expansion = if let Some(terms) = constraints.get(&node) {
        let mut composed = BTreeMap::<usize, f64>::new();
        for &(child, raw_weight) in terms {
            expand_terminal_node(child, constraints, node_ids, memo, stack)?;
            let child_expansion =
                memo.get(&child)
                    .ok_or_else(|| CutFemError::InvalidElasticityInput {
                        what: format!("hanging child {child:?} did not resolve to terminals"),
                    })?;
            for &(id, child_weight) in child_expansion {
                let contribution = raw_weight * child_weight;
                if !contribution.is_finite() {
                    return Err(CutFemError::InvalidElasticityInput {
                        what: format!(
                            "hanging transform contribution through {node:?} and {child:?} is non-finite"
                        ),
                    });
                }
                let weight = composed.entry(id).or_insert(0.0);
                *weight += contribution;
                if !weight.is_finite() {
                    return Err(CutFemError::InvalidElasticityInput {
                        what: format!(
                            "composed hanging transform for node {node:?} has a non-finite terminal weight"
                        ),
                    });
                }
            }
        }
        let expansion: NodeExpansion = composed.into_iter().collect();
        validate_affine_terms(node, &expansion, "composed")?;
        expansion
    } else {
        let id =
            node_ids
                .get(&node)
                .copied()
                .ok_or_else(|| CutFemError::InvalidElasticityInput {
                    what: format!(
                        "active node {node:?} has neither a constraint nor a terminal block"
                    ),
                })?;
        vec![(id, 1.0)]
    };
    stack.remove(&node);
    memo.insert(node, expansion);
    Ok(())
}

fn validate_affine_terms<K: core::fmt::Debug>(
    node: NodeKey,
    terms: &[(K, f64)],
    stage: &str,
) -> Result<(), CutFemError> {
    if terms.is_empty() {
        return Err(CutFemError::InvalidElasticityInput {
            what: format!("{stage} hanging transform for node {node:?} is empty"),
        });
    }
    if let Some((key, weight)) = terms.iter().find(|(_, weight)| !weight.is_finite()) {
        return Err(CutFemError::InvalidElasticityInput {
            what: format!(
                "{stage} hanging transform for node {node:?} has non-finite weight {weight} at {key:?}"
            ),
        });
    }
    let sum: f64 = terms.iter().map(|(_, weight)| weight).sum();
    if sum.to_bits() != 1.0f64.to_bits() {
        return Err(CutFemError::InvalidElasticityInput {
            what: format!(
                "{stage} hanging transform for node {node:?} has affine weight sum {sum}, expected exactly 1"
            ),
        });
    }
    Ok(())
}

fn clamp_terminal_blocks(
    node_ids: &TerminalIds,
    expansions: Option<&NodeExpansions>,
    selected: &BTreeSet<NodeKey>,
) -> Result<Vec<bool>, CutFemError> {
    let mut clamped = vec![false; 2 * node_ids.len()];
    for node in selected {
        if let Some(&id) = node_ids.get(node) {
            clamped[2 * id] = true;
            clamped[2 * id + 1] = true;
        }
    }
    if let Some(expansions) = expansions {
        for node in selected {
            if node_ids.contains_key(node) {
                continue;
            }
            let expansion =
                expansions
                    .get(node)
                    .ok_or_else(|| CutFemError::InvalidElasticityInput {
                        what: format!("selected clamp node {node:?} has no terminal expansion"),
                    })?;
            let unclamped: Vec<usize> = expansion
                .iter()
                .filter(|&&(id, _)| !clamped[2 * id])
                .map(|&(id, _)| id)
                .collect();
            if !unclamped.is_empty() {
                return Err(CutFemError::InvalidElasticityInput {
                    what: format!(
                        "clamp selects hanging node {node:?}, but terminal blocks {unclamped:?} are not all clamped"
                    ),
                });
            }
        }
    }
    Ok(clamped)
}

fn symmetrize_local(matrix: &mut [[f64; 8]; 8]) {
    let mut remaining_rows = matrix.as_mut_slice();
    let mut row = 0;
    while let Some((upper_row, lower_rows)) = remaining_rows.split_first_mut() {
        for (offset, lower_row) in lower_rows.iter_mut().enumerate() {
            let column = row + offset + 1;
            let value = f64::midpoint(upper_row[column], lower_row[row]);
            upper_row[column] = value;
            lower_row[row] = value;
        }
        remaining_rows = lower_rows;
        row += 1;
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

fn scatter_local_reduced(
    coo: &mut Coo,
    rhs: &mut [f64],
    clamped: &[bool],
    corners: &[NodeKey; 4],
    expansions: &NodeExpansions,
    local_k: &[[f64; 8]; 8],
    local_f: &[f64; 8],
) -> Result<(), CutFemError> {
    let mut reduced_rhs = BTreeMap::<usize, f64>::new();
    let mut reduced_matrix = BTreeMap::<(usize, usize), f64>::new();
    for a in 0..8 {
        let row_expansion = &expansions[&corners[a / 2]];
        for &(row_id, row_weight) in row_expansion {
            let row = 2 * row_id + a % 2;
            if clamped[row] {
                continue;
            }
            *reduced_rhs.entry(row).or_insert(0.0) += row_weight * local_f[a];
            for b in 0..8 {
                let value = local_k[a][b];
                if value == 0.0 {
                    continue;
                }
                for &(column_id, column_weight) in &expansions[&corners[b / 2]] {
                    let column = 2 * column_id + b % 2;
                    if !clamped[column] {
                        *reduced_matrix.entry((row, column)).or_insert(0.0) +=
                            value * (row_weight * column_weight);
                    }
                }
            }
        }
    }
    // Constraint collisions can make the transpose pair accumulate the same
    // analytical terms in a different order. Canonicalize each reduced cell
    // pair before COO insertion, just as `symmetrize_local` does before the
    // legacy scatter. Cell traversal is shared by both orientations, so the
    // final duplicate-accumulation sequences remain bit-identical.
    let off_diagonal: BTreeSet<(usize, usize)> = reduced_matrix
        .keys()
        .copied()
        .filter(|(row, column)| row != column)
        .map(|(row, column)| (row.min(column), row.max(column)))
        .collect();
    for (row, column) in off_diagonal {
        let Some(&upper) = reduced_matrix.get(&(row, column)) else {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!(
                    "reduced elasticity cell contribution is missing transpose entry ({row}, {column})"
                ),
            });
        };
        let Some(&lower) = reduced_matrix.get(&(column, row)) else {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!(
                    "reduced elasticity cell contribution is missing transpose entry ({column}, {row})"
                ),
            });
        };
        let value = f64::midpoint(upper, lower);
        if !value.is_finite() {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!(
                    "reduced elasticity cell contribution is non-finite at ({row}, {column})"
                ),
            });
        }
        reduced_matrix.insert((row, column), value);
        reduced_matrix.insert((column, row), value);
    }
    for (row, value) in reduced_rhs {
        if !value.is_finite() {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!("reduced elasticity load is non-finite at terminal DOF {row}"),
            });
        }
        rhs[row] += value;
    }
    for ((row, column), value) in reduced_matrix {
        if !value.is_finite() {
            return Err(CutFemError::InvalidElasticityInput {
                what: format!(
                    "reduced elasticity cell contribution is non-finite at ({row}, {column})"
                ),
            });
        }
        if value != 0.0 {
            coo.push(row, column, value);
        }
    }
    Ok(())
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hanging_transform_refuses_cycles() {
        let nodes = BTreeSet::from([(0, 0), (1, 0), (2, 0)]);
        let constraints =
            BTreeMap::from([((0, 0), vec![((1, 0), 1.0)]), ((1, 0), vec![((0, 0), 1.0)])]);
        assert!(matches!(
            build_terminal_expansions(&nodes, &constraints),
            Err(CutFemError::ConstraintCycle { .. })
        ));
    }

    #[test]
    fn hanging_transform_refuses_nonfinite_and_nonunit_weights() {
        let nodes = BTreeSet::from([(0, 0), (1, 0), (2, 0)]);
        for weights in [[f64::NAN, 0.5], [0.25, 0.5]] {
            let constraints =
                BTreeMap::from([((1, 0), vec![((0, 0), weights[0]), ((2, 0), weights[1])])]);
            assert!(matches!(
                build_terminal_expansions(&nodes, &constraints),
                Err(CutFemError::InvalidElasticityInput { .. })
            ));
        }
    }

    #[test]
    fn hanging_transform_refuses_empty_and_absent_children() {
        let nodes = BTreeSet::from([(0, 0), (1, 0), (2, 0)]);
        for constraints in [
            BTreeMap::from([((1, 0), Vec::new())]),
            BTreeMap::from([((1, 0), vec![((0, 0), 0.5), ((3, 0), 0.5)])]),
        ] {
            assert!(matches!(
                build_terminal_expansions(&nodes, &constraints),
                Err(CutFemError::InvalidElasticityInput { .. })
            ));
        }
    }

    #[test]
    fn hanging_clamp_refuses_a_partially_clamped_terminal_trace() {
        // A valid leaf partition cannot place a hanging midpoint on the
        // exterior design-box edge, which is the public clamp predicate's
        // current scope. Exercise the generic transform/clamp invariant here
        // with the same synthetic corruption seam used for cycle coverage.
        let nodes = BTreeSet::from([(0, 0), (1, 0), (2, 0)]);
        let constraints = BTreeMap::from([((1, 0), vec![((0, 0), 0.5), ((2, 0), 0.5)])]);
        let (node_ids, expansions) =
            build_terminal_expansions(&nodes, &constraints).expect("valid midpoint transform");
        let selected = BTreeSet::from([(0, 0), (1, 0)]);
        let error = clamp_terminal_blocks(&node_ids, Some(&expansions), &selected)
            .expect_err("one unclamped representative must refuse");
        assert!(
            matches!(&error, CutFemError::InvalidElasticityInput { what } if what.contains("not all clamped")),
            "unexpected refusal: {error}"
        );

        let all_selected = BTreeSet::from([(0, 0), (1, 0), (2, 0)]);
        let all_clamped = clamp_terminal_blocks(&node_ids, Some(&expansions), &all_selected)
            .expect("a wholly clamped midpoint trace is compatible");
        assert!(all_clamped.iter().all(|value| *value));

        let terminal_only = BTreeSet::from([(0, 0)]);
        let one_clamped = clamp_terminal_blocks(&node_ids, Some(&expansions), &terminal_only)
            .expect("clamping one terminal alone is compatible");
        let selected_id = node_ids[&(0, 0)];
        assert!(one_clamped[2 * selected_id] && one_clamped[2 * selected_id + 1]);
        assert_eq!(one_clamped.iter().filter(|&&value| value).count(), 2);
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn equal_level_ghost_face_applies_nonidentity_terminal_transform() {
        let mut grid = Quadtree::with_room(1, 2);
        grid.split((1, 1, 0));
        let active: BTreeSet<_> = grid.leaves().collect();
        let nodes: BTreeSet<_> = active
            .iter()
            .flat_map(|&cell| grid.corner_nodes(cell))
            .collect();
        let constraints: BTreeMap<_, Vec<_>> = grid
            .hanging_constraints(&active, &nodes)
            .into_iter()
            .map(|(node, terms)| (node, terms.into()))
            .collect();
        let (terminal_ids, expansions) =
            build_terminal_expansions(&nodes, &constraints).expect("valid 2:1 transform");

        let midpoint = (2, 1);
        assert_eq!(
            expansions[&midpoint],
            vec![(terminal_ids[&(2, 0)], 0.5), (terminal_ids[&(2, 2)], 0.5),],
            "fixture midpoint must use the real quadtree hanging constraint"
        );

        let face_cells = ((2, 2, 0), (2, 3, 0));
        assert!(grid.is_leaf(face_cells.0) && grid.is_leaf(face_cells.1));
        assert_eq!(
            face_cells.0.0, face_cells.1.0,
            "ghost face must be equal-level"
        );
        let face = grid
            .shared_face_patch(face_cells.0, face_cells.1)
            .expect("exact equal-level face patch");
        assert!(
            grid.corner_nodes(face_cells.0).contains(&midpoint),
            "constrained midpoint must participate in the face kernel"
        );

        let physical_ids: TerminalIds = nodes
            .iter()
            .copied()
            .enumerate()
            .map(|(id, node)| (node, id))
            .collect();
        let material = IsotropicElastic::new(1.0, 0.3, 10.0).expect("valid fixture material");
        let sdf = crate::sdf::HalfPlane {
            normal: [1.0, 0.0],
            offset: 0.75,
        };
        let problem = CutElasticity {
            grid: &grid,
            sdf: &sdf,
            material: &material,
            nitsche_beta: 100.0,
            ghost_gamma: 0.5,
            quad_depth: 4,
            clamp: None,
            boundary_traction: None,
            traction_free_interface: false,
            solver_tol: 1e-10,
            solver_max_iters: 100,
        };
        let (_, mu) = material.lame();

        let physical_dofs = 2 * physical_ids.len();
        let physical_clamped = vec![false; physical_dofs];
        let mut physical_coo = Coo::new(physical_dofs, physical_dofs);
        problem
            .assemble_ghost_face(
                face,
                mu,
                &physical_ids,
                None,
                &physical_clamped,
                &mut physical_coo,
            )
            .expect("physical face assembly");
        let physical = physical_coo.assemble();
        let (_, midpoint_row) = physical.row(2 * physical_ids[&midpoint]);
        assert!(
            midpoint_row.iter().any(|value| *value != 0.0),
            "constrained midpoint must carry nonzero physical ghost stiffness"
        );

        let terminal_dofs = 2 * terminal_ids.len();
        let terminal_clamped = vec![false; terminal_dofs];
        let mut reduced_coo = Coo::new(terminal_dofs, terminal_dofs);
        problem
            .assemble_ghost_face(
                face,
                mu,
                &terminal_ids,
                Some(&expansions),
                &terminal_clamped,
                &mut reduced_coo,
            )
            .expect("reduced face assembly");
        let reduced = reduced_coo.assemble();

        let mut saw_nonidentity_effect = false;
        for basis_index in 0..terminal_dofs {
            let mut basis = vec![0.0; terminal_dofs];
            basis[basis_index] = 1.0;

            let mut transformed = vec![0.0; physical_dofs];
            for (&node, &physical_id) in &physical_ids {
                for &(terminal_id, weight) in &expansions[&node] {
                    for component in 0..2 {
                        transformed[2 * physical_id + component] +=
                            weight * basis[2 * terminal_id + component];
                    }
                }
            }
            let mut physical_applied = vec![0.0; physical_dofs];
            physical.spmv(&transformed, &mut physical_applied);
            let mut expected = vec![0.0; terminal_dofs];
            for (&node, &physical_id) in &physical_ids {
                for &(terminal_id, weight) in &expansions[&node] {
                    for component in 0..2 {
                        expected[2 * terminal_id + component] +=
                            weight * physical_applied[2 * physical_id + component];
                    }
                }
            }

            let mut actual = vec![0.0; terminal_dofs];
            reduced.spmv(&basis, &mut actual);
            for (row, (&actual_value, &expected_value)) in actual.iter().zip(&expected).enumerate()
            {
                let tolerance =
                    256.0 * f64::EPSILON * actual_value.abs().max(expected_value.abs()).max(1.0);
                assert!(
                    (actual_value - expected_value).abs() <= tolerance,
                    "T^T K T mismatch at ({row}, {basis_index}): \
                     actual={actual_value:e}, expected={expected_value:e}, tolerance={tolerance:e}"
                );
            }

            let mut untransformed = vec![0.0; physical_dofs];
            for (&node, &terminal_id) in &terminal_ids {
                let physical_id = physical_ids[&node];
                for component in 0..2 {
                    untransformed[2 * physical_id + component] = basis[2 * terminal_id + component];
                }
            }
            let mut untransformed_applied = vec![0.0; physical_dofs];
            physical.spmv(&untransformed, &mut untransformed_applied);
            for (&node, &terminal_id) in &terminal_ids {
                let physical_id = physical_ids[&node];
                for component in 0..2 {
                    let naive = untransformed_applied[2 * physical_id + component];
                    saw_nonidentity_effect |=
                        (expected[2 * terminal_id + component] - naive).abs() > 1e-12;
                }
            }
        }
        assert!(
            saw_nonidentity_effect,
            "hanging-node transform must differ from deleting the constrained block"
        );
    }
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
