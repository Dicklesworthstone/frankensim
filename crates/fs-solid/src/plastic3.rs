//! Incremental small-strain inelasticity on body-fitted linear tetrahedra
//! (plan §8.2 "J2 plasticity with return mapping"), generic over any
//! `fs-material` [`SmallStrainLaw`]: the law's algorithmic stress update and
//! consistent tangent drive a load-path-following Newton solve, and the
//! internal variables are committed only at converged increments. Uniform and
//! per-element material fields use the same assembly and Newton solve. Each
//! tetrahedron selects its own law for initialization, stress, tangent, and
//! history update; material regions are never averaged into one effective law.
//!
//! # Discretization
//!
//! P1 tetrahedra have one constant strain, so one law evaluation per element
//! is exact quadrature. With tensor-component Voigt strain
//! `eps = B u` (`[xx, yy, zz, xy, yz, zx]`, `eps_xy = (u_x,y + u_y,x)/2`)
//! and `W = diag(1, 1, 1, 2, 2, 2)` (the tensor double contraction),
//!
//! ```text
//! internal force  f = sum_e V B^T W sigma(eps, state_committed)
//! tangent         K = sum_e V B^T W C_alg B
//! ```
//!
//! # Load path
//!
//! `load_path` lists load factors applied in sequence; each scales the
//! reference prescribed displacements and dead nodal forces. Factors may
//! decrease (unloading) or change sign (load reversal). Every increment is
//! one Newton solve from the previous converged state with the prescribed
//! increment taken through the consistent linearization `K_ff d = -r_f -
//! K_fp dp`; a step that increases the residual is halved. The gate is
//! `||r_f|| <= relative_tolerance * force_scale`. The scale bounds actual
//! element stress forces on free rows, with the external free-force norm and
//! the largest accepted scale in this increment retained through unloading.
//! Reactions on unrelated fixed bodies and rigid translations do not weaken
//! the gate. After convergence each element history is committed.
//!
//! # No-claim boundaries
//!
//! Small strain and rotation (no finite-strain plasticity); no implicit clock.
//! Rate-dependent laws must supply a physical interval for each load increment.
//! The [`hex`] child supplies a locking-resistant eight-node B-bar alternative
//! using this same load-path/Newton owner.
//! Constant-strain tetrahedra (plastic incompressibility locks
//! P1 tetrahedra in confined flow: limit loads of confined problems are
//! over-predicted; no mixed/B-bar formulation), dead loads only, no
//! limit-load continuation (an increment past a limit load stalls and
//! refuses; drive by displacement), no mesh-convergence claim. Linear solves
//! as in `hyper3` (ILU(0)-PCG, |diag|-Jacobi PMINRES fallback, recomputed
//! residual acceptance). Deterministic sequential assembly.

use fs_exec::Cx;
use fs_material::{SmallStrainLaw, Tangent6, Voigt};
use fs_solver::krylov::{PminresState, StallDiagnosis};
use fs_solver::nonlinear::{FgmresState, SolverRunProgress};
use fs_solver::op::CsrOp;
use fs_sparse::precond::{Ilu0, Precond, ilu0, pcg};
use fs_sparse::{Coo, Csr};

use crate::linear3::TetAssemblyBudget;

/// Physical-time creep increments with independent work and heat accounting.
pub mod creep;
/// Inertial increments using the same constitutive assembly and Newton owner.
pub mod dynamics;

mod increment;

/// Fully integrated eight-node hexahedra with optional volumetric B-bar projection.
pub mod hex;

/// Fixed-plane strain-life observation of the native cyclic J2 load path.
pub mod fatigue;

/// Newton controls for each increment.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IncrementSettings {
    /// Newton iterations allowed per increment.
    pub max_newton_iterations: usize,
    /// Relative force-residual gate (see module docs).
    pub relative_tolerance: f64,
    /// Residual-increase halvings allowed per iteration.
    pub max_backtracks: usize,
    /// Relative residual required of each linear solve.
    pub linear_tolerance: f64,
    /// Iteration cap per Krylov path. A symmetric correction may try PCG and
    /// then MINRES; any true-residual MINRES correction shares that path's
    /// original allowance. The general FGMRES path uses one such allowance.
    pub max_linear_iterations: usize,
}

impl Default for IncrementSettings {
    fn default() -> Self {
        Self {
            max_newton_iterations: 30,
            relative_tolerance: 1e-10,
            max_backtracks: 20,
            linear_tolerance: 1e-12,
            max_linear_iterations: 20_000,
        }
    }
}

/// Constitutive data over a continuum mesh.
///
/// Use `Uniform` for one material and `PerElement` for different material
/// cards or spatially varying states of the same constitutive family. Entries
/// follow the element connectivity order. All coefficients remain fixed
/// over this problem's load path; each element evolves its own internal state.
/// Retain source receipts alongside these laws when lowering material cards.
#[derive(Debug, Clone, Copy)]
pub enum SmallStrainLawField<'a, L: SmallStrainLaw> {
    /// One constitutive law for every element.
    Uniform(&'a L),
    /// One constitutive law per element, without homogenization.
    PerElement(&'a [L]),
}

impl<L: SmallStrainLaw> SmallStrainLawField<'_, L> {
    fn validate_count(&self, expected: usize) -> Result<(), SmallStrainTetError> {
        if let Self::PerElement(laws) = self
            && laws.len() != expected
        {
            return Err(SmallStrainTetError::LawCountMismatch {
                expected,
                got: laws.len(),
            });
        }
        Ok(())
    }

    fn at(&self, element: usize) -> &L {
        match self {
            Self::Uniform(law) => law,
            Self::PerElement(laws) => &laws[element],
        }
    }
}

/// A load-path problem on a tetrahedral mesh.
#[derive(Debug, Clone, Copy)]
pub struct SmallStrainTetProblem<'a, L: SmallStrainLaw> {
    /// Vertex coordinates, m.
    pub nodes_m: &'a [[f64; 3]],
    /// Four vertex indices per conforming tetrahedron.
    pub tetrahedra: &'a [[usize; 4]],
    /// Uniform or per-element constitutive laws, in connectivity order.
    pub laws: SmallStrainLawField<'a, L>,
    /// Reference prescribed displacements: (`3 node + component`, m).
    pub prescribed_m: &'a [(usize, f64)],
    /// Reference dead nodal forces: (`3 node + component`, N).
    pub nodal_forces_n: &'a [(usize, f64)],
    /// Load factors applied in sequence.
    pub load_path: &'a [f64],
    /// Size and element-quality envelope.
    pub budget: TetAssemblyBudget,
    /// Newton controls.
    pub settings: IncrementSettings,
}

/// Evidence of one converged increment.
#[derive(Debug, Clone, PartialEq)]
pub struct IncrementRecord {
    /// Load factor reached.
    pub load_factor: f64,
    /// Relative residual before each Newton update.
    pub residual_history: Vec<f64>,
    /// Reactions at the prescribed DOFs (internal minus external), N, in
    /// `prescribed_m` order.
    pub reactions_n: Vec<f64>,
    /// Elements whose committed state changed in this increment.
    pub evolving_elements: usize,
    /// Linear solves that fell back to PMINRES.
    pub indefinite_fallbacks: usize,
}

/// State at the end of the load path.
#[derive(Debug, Clone, PartialEq)]
pub struct SmallStrainTetSolution<S> {
    /// Nodal displacements, m.
    pub displacement_m: Vec<[f64; 3]>,
    /// Element strains (tensor-component Voigt).
    pub element_strain: Vec<Voigt>,
    /// Element stresses (Voigt), Pa.
    pub element_stress: Vec<Voigt>,
    /// Committed element internal variables.
    pub element_states: Vec<S>,
    /// One record per load-path entry.
    pub increments: Vec<IncrementRecord>,
}

/// Structured refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum SmallStrainTetError {
    /// Inconsistent, out-of-budget, or non-finite input.
    InvalidInput {
        /// What was refused.
        what: String,
    },
    /// A per-element field does not cover exactly the mesh connectivity.
    LawCountMismatch {
        /// Number of elements.
        expected: usize,
        /// Number of supplied constitutive laws.
        got: usize,
    },
    /// A reference tetrahedron is degenerate or inverted.
    DegenerateElement {
        /// Element index.
        element: usize,
        /// `det(J) / h_max^3`.
        scaled_jacobian: f64,
    },
    /// Newton did not reach the gate for this increment. Repair: refine the
    /// load path; past a limit load, drive by displacement.
    NewtonStalled {
        /// Index into `load_path`.
        increment: usize,
        /// Relative residual history.
        history: Vec<f64>,
    },
    /// A linear solve failed on both paths.
    LinearSolveFailed {
        /// Index into `load_path`.
        increment: usize,
        /// Recomputed relative residual.
        relative_residual: f64,
    },
    /// The context was cancelled; nothing was committed past the last
    /// converged increment.
    Cancelled,
}

impl core::fmt::Display for SmallStrainTetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput { what } => write!(f, "FS-SOLID-PLASTIC3-INPUT: {what}"),
            Self::LawCountMismatch { expected, got } => write!(
                f,
                "FS-SOLID-PLASTIC3-LAW-COUNT: expected {expected} laws, got {got}"
            ),
            Self::DegenerateElement {
                element,
                scaled_jacobian,
            } => write!(
                f,
                "FS-SOLID-PLASTIC3-DEGENERATE: element {element} scaled Jacobian {scaled_jacobian:e}"
            ),
            Self::NewtonStalled { increment, history } => write!(
                f,
                "FS-SOLID-PLASTIC3-NEWTON: increment {increment} stalled after {} iterations; refine the load path",
                history.len()
            ),
            Self::LinearSolveFailed {
                increment,
                relative_residual,
            } => write!(
                f,
                "FS-SOLID-PLASTIC3-LINEAR: increment {increment}: residual {relative_residual:e}"
            ),
            Self::Cancelled => write!(f, "FS-SOLID-PLASTIC3-CANCELLED"),
        }
    }
}

impl std::error::Error for SmallStrainTetError {}

#[derive(Debug, Clone, Copy)]
struct Element {
    nodes: [usize; 4],
    gradients: [[f64; 3]; 4],
    volume: f64,
}

const WEIGHT: [f64; 6] = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0];

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).fold(0.0, |s, (x, y)| x.mul_add(*y, s))
}

fn norm(a: &[f64]) -> f64 {
    fs_math::det::sqrt(dot(a, a))
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), SmallStrainTetError> {
    cx.checkpoint().map_err(|_| SmallStrainTetError::Cancelled)
}

/// `B` column of DOF `(a, i)`: d eps / d u_(a,i) in tensor-component Voigt.
fn b_column(g: [f64; 3], i: usize) -> Voigt {
    let mut col = [0.0; 6];
    col[i] = g[i];
    // Shear slots: 3 = xy, 4 = yz, 5 = zx.
    match i {
        0 => {
            col[3] = 0.5 * g[1];
            col[5] = 0.5 * g[2];
        }
        1 => {
            col[3] = 0.5 * g[0];
            col[4] = 0.5 * g[2];
        }
        _ => {
            col[4] = 0.5 * g[1];
            col[5] = 0.5 * g[0];
        }
    }
    col
}

impl<L: SmallStrainLaw> SmallStrainTetProblem<'_, L> {
    fn validate(&self) -> Result<Vec<Element>, SmallStrainTetError> {
        let invalid = |what: String| Err(SmallStrainTetError::InvalidInput { what });
        let n = self.nodes_m.len();
        if n == 0 || self.tetrahedra.is_empty() || self.load_path.is_empty() {
            return invalid("empty mesh or load path".into());
        }
        self.laws.validate_count(self.tetrahedra.len())?;
        if n > self.budget.maximum_nodes || self.tetrahedra.len() > self.budget.maximum_tetrahedra {
            return invalid(format!("mesh exceeds the budget {:?}", self.budget));
        }
        if self
            .nodes_m
            .iter()
            .flatten()
            .chain(self.load_path)
            .any(|v| !v.is_finite())
        {
            return invalid("non-finite node coordinate or load factor".into());
        }
        let s = self.settings;
        if s.max_newton_iterations == 0
            || !(s.relative_tolerance.is_finite() && s.relative_tolerance > 0.0)
            || !(s.linear_tolerance.is_finite() && s.linear_tolerance > 0.0)
        {
            return invalid(format!("inadmissible settings {s:?}"));
        }
        let dofs = 3 * n;
        let mut seen = vec![false; dofs];
        for &(dof, value) in self.prescribed_m {
            if dof >= dofs || !value.is_finite() || seen[dof] {
                return invalid(format!(
                    "prescribed DOF {dof} = {value} invalid or repeated"
                ));
            }
            seen[dof] = true;
        }
        for &(dof, value) in self.nodal_forces_n {
            if dof >= dofs || !value.is_finite() {
                return invalid(format!("nodal force on DOF {dof} = {value} invalid"));
            }
        }
        if dofs - self.prescribed_m.len() > self.budget.maximum_free_dofs {
            return invalid("free DOFs exceed the budget".into());
        }
        let mut elements = Vec::with_capacity(self.tetrahedra.len());
        for (element, tet) in self.tetrahedra.iter().enumerate() {
            if tet.iter().any(|&v| v >= n) {
                return invalid(format!("element {element} names a missing vertex"));
            }
            let x = tet.map(|v| self.nodes_m[v]);
            let e = |k: usize| [0, 1, 2].map(|d| x[k][d] - x[0][d]);
            let (e1, e2, e3) = (e(1), e(2), e(3));
            let cross = |a: [f64; 3], b: [f64; 3]| {
                [
                    a[1] * b[2] - a[2] * b[1],
                    a[2] * b[0] - a[0] * b[2],
                    a[0] * b[1] - a[1] * b[0],
                ]
            };
            let c23 = cross(e2, e3);
            let det = dot(&c23, &e1);
            let mut h_max = 0.0f64;
            for a in 0..4 {
                for b in a + 1..4 {
                    let d = [0, 1, 2].map(|k| x[b][k] - x[a][k]);
                    h_max = h_max.max(norm(&d));
                }
            }
            let scaled = det / (h_max * h_max * h_max);
            if !(scaled.is_finite() && scaled >= self.budget.minimum_scaled_jacobian) {
                return Err(SmallStrainTetError::DegenerateElement {
                    element,
                    scaled_jacobian: scaled,
                });
            }
            let g1 = c23.map(|v| v / det);
            let g2 = cross(e3, e1).map(|v| v / det);
            let g3 = cross(e1, e2).map(|v| v / det);
            let g0 = [0, 1, 2].map(|k| -(g1[k] + g2[k] + g3[k]));
            elements.push(Element {
                nodes: *tet,
                gradients: [g0, g1, g2, g3],
                volume: det / 6.0,
            });
        }
        Ok(elements)
    }

    fn strain(element: &Element, u: &[f64]) -> Voigt {
        let mut eps = [0.0; 6];
        // grad(N_0) = -sum_{a>0} grad(N_a): subtract translation before
        // contraction rather than cancelling large absolute displacements.
        for (a, &node) in element.nodes.iter().enumerate().skip(1) {
            for i in 0..3 {
                let col = b_column(element.gradients[a], i);
                let ui = u[3 * node + i] - u[3 * element.nodes[0] + i];
                for (e, c) in eps.iter_mut().zip(col) {
                    *e = c.mul_add(ui, *e);
                }
            }
        }
        eps
    }

    fn internal(
        &self,
        elements: &[Element],
        states: &[L::State],
        u: &[f64],
        prescribed: Option<&[bool]>,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, f64), SmallStrainTetError> {
        let mut f = vec![0.0f64; u.len()];
        let mut force_bounds = prescribed.map(|_| vec![0.0; u.len()]);
        for (index, element) in elements.iter().enumerate() {
            if index % 4096 == 0 {
                checkpoint(cx)?;
            }
            let strain = Self::strain(element, u);
            let sigma = self.laws.at(index).stress(&strain, &states[index]);
            if strain.iter().chain(&sigma).any(|v| !v.is_finite()) {
                return Err(SmallStrainTetError::InvalidInput {
                    what: format!("element {index} has non-finite strain or constitutive stress"),
                });
            }
            let stress_norm = if prescribed.is_some() {
                sigma.iter().enumerate().fold(0.0_f64, |n, (v, s)| {
                    let n = n.hypot(*s);
                    if v >= 3 { n.hypot(*s) } else { n }
                })
            } else {
                0.0
            };
            for (a, &node) in element.nodes.iter().enumerate() {
                for i in 0..3 {
                    let col = b_column(element.gradients[a], i);
                    let mut s = 0.0;
                    for v in 0..6 {
                        s = (col[v] * WEIGHT[v]).mul_add(sigma[v], s);
                    }
                    f[3 * node + i] = element.volume.mul_add(s, f[3 * node + i]);
                    if let (Some(prescribed), Some(bounds)) = (prescribed, &mut force_bounds)
                        && !prescribed[3 * node + i]
                    {
                        // Cauchy-Schwarz for B_column : sigma in the tensor
                        // metric. It remains nonzero at a traction-free row
                        // when individual stress components cancel there.
                        let b_norm =
                            fs_math::det::sqrt((0..6).map(|v| WEIGHT[v] * col[v] * col[v]).sum());
                        bounds[3 * node + i] += element.volume * b_norm * stress_norm;
                    }
                }
            }
        }
        let scale = force_bounds.as_ref().map_or(0.0, |bounds| {
            bounds.iter().fold(0.0_f64, |n, v| n.hypot(*v))
        });
        if !scale.is_finite() || f.iter().any(|v| !v.is_finite()) {
            return Err(SmallStrainTetError::InvalidInput {
                what: "non-finite assembled constitutive force or force scale".into(),
            });
        }
        Ok((f, scale))
    }

    fn tangent(
        &self,
        elements: &[Element],
        states: &[L::State],
        u: &[f64],
        inertia: Option<&dynamics::Inertia<'_>>,
        cx: &Cx<'_>,
    ) -> Result<Csr, SmallStrainTetError> {
        let mut coo = Coo::new(u.len(), u.len());
        for (index, element) in elements.iter().enumerate() {
            if index % 4096 == 0 {
                checkpoint(cx)?;
            }
            let c: Tangent6 = self
                .laws
                .at(index)
                .tangent(&Self::strain(element, u), &states[index]);
            if c.iter().flatten().any(|v| !v.is_finite()) {
                return Err(SmallStrainTetError::InvalidInput {
                    what: format!("element {index} has a non-finite constitutive tangent"),
                });
            }
            let cols: [Voigt; 12] =
                std::array::from_fn(|k| b_column(element.gradients[k / 3], k % 3));
            // C B_q is shared by all twelve test columns. Keep both small
            // arrays on the stack instead of allocating once per element.
            let cb: [Voigt; 12] = cols.map(|bq| {
                std::array::from_fn(|v| (0..6).fold(0.0, |sum, w| c[v][w].mul_add(bq[w], sum)))
            });
            for (p, bp) in cols.iter().enumerate() {
                let row = 3 * element.nodes[p / 3] + p % 3;
                for (q, cbq) in cb.iter().enumerate() {
                    let mut s = 0.0;
                    for v in 0..6 {
                        s = (bp[v] * WEIGHT[v]).mul_add(cbq[v], s);
                    }
                    coo.push(row, 3 * element.nodes[q / 3] + q % 3, element.volume * s);
                }
            }
        }
        if let Some(inertia) = inertia {
            inertia.add_tangent(&mut coo, cx)?;
        }
        let matrix = coo.assemble();
        if (0..matrix.nrows()).any(|row| matrix.row(row).1.iter().any(|v| !v.is_finite())) {
            return Err(SmallStrainTetError::InvalidInput {
                what: "non-finite assembled constitutive tangent".into(),
            });
        }
        Ok(matrix)
    }

    /// Internal force at `displacement_m` (all DOFs) with the given committed
    /// element states, exposed for verification.
    ///
    /// # Errors
    /// Input/budget refusals or cancellation.
    pub fn internal_force(
        &self,
        displacement_m: &[f64],
        states: &[L::State],
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, SmallStrainTetError> {
        let elements = self.validate()?;
        self.check(displacement_m, states)?;
        self.internal(&elements, states, displacement_m, None, cx)
            .map(|(force, _)| force)
    }

    /// Algorithmic tangent at `displacement_m` with the given committed
    /// element states, exposed for verification.
    ///
    /// # Errors
    /// Input/budget refusals or cancellation.
    pub fn tangent_matrix(
        &self,
        displacement_m: &[f64],
        states: &[L::State],
        cx: &Cx<'_>,
    ) -> Result<Csr, SmallStrainTetError> {
        let elements = self.validate()?;
        self.check(displacement_m, states)?;
        self.tangent(&elements, states, displacement_m, None, cx)
    }

    /// Reference displacement gradients from the same P1 shape gradients as
    /// the mechanical weak form. Coupled owners can bound small rotations as
    /// well as symmetric strain, without reconstructing a second FE geometry.
    pub fn displacement_gradients(
        &self,
        displacement_m: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<[[f64; 3]; 3]>, SmallStrainTetError> {
        checkpoint(cx)?;
        let elements = self.validate()?;
        if displacement_m.len() != 3 * self.nodes_m.len()
            || displacement_m.iter().any(|v| !v.is_finite())
        {
            return Err(SmallStrainTetError::InvalidInput {
                what: "invalid displacement field for element kinematics".into(),
            });
        }
        let mut gradients = Vec::with_capacity(elements.len());
        for element in &elements {
            checkpoint(cx)?;
            let mut gradient = [[0.0; 3]; 3];
            for (a, &node) in element.nodes.iter().enumerate().skip(1) {
                for i in 0..3 {
                    for j in 0..3 {
                        gradient[i][j] = (displacement_m[3 * node + i]
                            - displacement_m[3 * element.nodes[0] + i])
                            .mul_add(element.gradients[a][j], gradient[i][j]);
                    }
                }
            }
            gradients.push(gradient);
        }
        Ok(gradients)
    }

    fn check(&self, u: &[f64], states: &[L::State]) -> Result<(), SmallStrainTetError> {
        if u.len() != 3 * self.nodes_m.len()
            || states.len() != self.tetrahedra.len()
            || u.iter().any(|v| !v.is_finite())
        {
            return Err(SmallStrainTetError::InvalidInput {
                what: "displacement or state arrays have the wrong length or are non-finite".into(),
            });
        }
        Ok(())
    }

    /// Follow the load path.
    ///
    /// # Errors
    /// See [`SmallStrainTetError`].
    pub fn solve(
        &self,
        cx: &Cx<'_>,
    ) -> Result<SmallStrainTetSolution<L::State>, SmallStrainTetError> {
        self.solve_candidate(None, None, cx)
    }

    /// Continue the same native Newton solve from accepted displacement and
    /// element history. Load factors still specify absolute targets, not
    /// increments added to the accepted displacement. All proposals borrow the
    /// same old history until their load increment converges; refusal publishes
    /// none of the candidate state.
    ///
    /// This is the ingress for thermal/material recoupling between increments.
    /// The constitutive owner must admit the retained history, coefficient
    /// changes and energy transfer before accepting the resulting solution.
    pub fn solve_from(
        &self,
        displacement_m: &[f64],
        element_states: &[L::State],
        cx: &Cx<'_>,
    ) -> Result<SmallStrainTetSolution<L::State>, SmallStrainTetError> {
        self.solve_candidate(Some((displacement_m, element_states)), None, cx)
    }

    fn solve_candidate(
        &self,
        initial: Option<(&[f64], &[L::State])>,
        inertia: Option<&dynamics::Inertia<'_>>,
        cx: &Cx<'_>,
    ) -> Result<SmallStrainTetSolution<L::State>, SmallStrainTetError> {
        checkpoint(cx)?;
        let elements = self.validate()?;
        if let Some((displacement, states)) = initial {
            self.check(displacement, states)?;
        }
        let operator = TetIncrement {
            problem: self,
            elements: &elements,
            inertia,
        };
        let result = increment::solve(
            &operator,
            increment::IncrementControls {
                dofs: 3 * self.nodes_m.len(),
                prescribed_m: self.prescribed_m,
                nodal_forces_n: self.nodal_forces_n,
                load_path: self.load_path,
                settings: self.settings,
            },
            initial,
            inertia,
            cx,
        )?;
        Ok(SmallStrainTetSolution {
            displacement_m: result.displacement_m.as_chunks::<3>().0.to_vec(),
            element_strain: elements
                .iter()
                .map(|e| Self::strain(e, &result.displacement_m))
                .collect(),
            element_stress: result.point_stress,
            element_states: result.point_states,
            increments: result.increments,
        })
    }
}

struct TetIncrement<'a, 'b, L: SmallStrainLaw> {
    problem: &'a SmallStrainTetProblem<'b, L>,
    elements: &'a [Element],
    inertia: Option<&'a dynamics::Inertia<'a>>,
}

impl<L: SmallStrainLaw> increment::IncrementOperator for TetIncrement<'_, '_, L> {
    type State = L::State;

    fn initial_states(&self) -> Vec<Self::State> {
        (0..self.elements.len())
            .map(|index| self.problem.laws.at(index).initial_state())
            .collect()
    }

    fn internal(
        &self,
        states: &[Self::State],
        u: &[f64],
        prescribed: Option<&[bool]>,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, f64), SmallStrainTetError> {
        self.problem
            .internal(self.elements, states, u, prescribed, cx)
    }

    fn tangent(
        &self,
        states: &[Self::State],
        u: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Csr, SmallStrainTetError> {
        self.problem
            .tangent(self.elements, states, u, self.inertia, cx)
    }

    fn tangent_symmetric(&self) -> bool {
        match self.problem.laws {
            SmallStrainLawField::Uniform(law) => law.admissibility().tangent_symmetric,
            SmallStrainLawField::PerElement(laws) => {
                laws.iter().all(|law| law.admissibility().tangent_symmetric)
            }
        }
    }

    fn commit(
        &self,
        states: &mut [Self::State],
        u: &[f64],
        mut final_stress: Option<&mut Vec<Voigt>>,
        cx: &Cx<'_>,
    ) -> Result<usize, SmallStrainTetError> {
        let mut evolving = 0;
        for (index, element) in self.elements.iter().enumerate() {
            if index % 4096 == 0 {
                checkpoint(cx)?;
            }
            let strain = SmallStrainTetProblem::<L>::strain(element, u);
            let law = self.problem.laws.at(index);
            if let Some(stresses) = &mut final_stress {
                stresses.push(law.stress(&strain, &states[index]));
            }
            let next = law.update_state(&strain, &states[index]);
            evolving += usize::from(next != states[index]);
            states[index] = next;
        }
        Ok(evolving)
    }
}

/// Shared sparse linear correction for small-strain continuum increments.
pub(crate) fn solve_increment_linear(
    kff: &Csr,
    rhs: &[f64],
    increment: usize,
    settings: IncrementSettings,
    symmetric: bool,
    cx: &Cx<'_>,
) -> Result<(Vec<f64>, bool), SmallStrainTetError> {
    checkpoint(cx)?;
    let tol = settings.linear_tolerance;
    if norm(rhs) == 0.0 {
        return Ok((vec![0.0; rhs.len()], false));
    }
    let true_relative = |x: &[f64]| {
        let mut ax = vec![0.0f64; rhs.len()];
        kff.spmv(x, &mut ax);
        let r: Vec<f64> = rhs.iter().zip(&ax).map(|(b, a)| b - a).collect();
        norm(&r) / norm(rhs)
    };
    // Dynamic-recovery plasticity has a genuinely nonsymmetric consistent
    // tangent under nonproportional loading. Preserve it and select the existing
    // general Krylov owner; neither PCG nor MINRES is admitted for this operator.
    if !symmetric {
        if settings.max_linear_iterations == 0 {
            return Err(SmallStrainTetError::LinearSolveFailed {
                increment,
                relative_residual: 1.0,
            });
        }
        let preconditioner = match ilu0(kff) {
            Ok(ilu) => GeneralPreconditioner::Ilu(ilu),
            Err(_) => GeneralPreconditioner::Jacobi(AbsJacobi::new(kff)),
        };
        let op = CsrOp::general(kff.clone());
        let restart = rhs.len().min(32).min(settings.max_linear_iterations);
        let mut state = FgmresState::new(rhs, restart);
        while state.iters < settings.max_linear_iterations {
            state.restart = restart.min(settings.max_linear_iterations - state.iters);
            let before = state.iters;
            let report = state.run_cancellable(&op, &preconditioner, rhs, tol, 1, cx);
            match report.progress {
                SolverRunProgress::Complete => {}
                SolverRunProgress::Paused => return Err(SmallStrainTetError::Cancelled),
                _ => {
                    return Err(SmallStrainTetError::LinearSolveFailed {
                        increment,
                        relative_residual: f64::INFINITY,
                    });
                }
            }
            if report.report.converged
                || report.report.diagnosis == Some(StallDiagnosis::Breakdown)
                || state.iters == before
            {
                break;
            }
        }
        let relative_residual = true_relative(&state.x);
        checkpoint(cx)?;
        return if relative_residual.is_finite()
            && relative_residual <= 10.0 * tol
            && state.x.iter().all(|value| value.is_finite())
        {
            Ok((state.x, false))
        } else {
            Err(SmallStrainTetError::LinearSolveFailed {
                increment,
                relative_residual,
            })
        };
    }
    if let Ok(m) = ilu0(kff) {
        let mut x = vec![0.0f64; rhs.len()];
        let report = pcg(kff, rhs, &mut x, &m, tol, settings.max_linear_iterations);
        if report.converged && x.iter().all(|v| v.is_finite()) && true_relative(&x) <= 10.0 * tol {
            return Ok((x, false));
        }
    }
    let jacobi = AbsJacobi::new(kff);
    let op = CsrOp::symmetric(kff.clone());
    let mut state = PminresState::new(&op, &jacobi, rhs);
    let report = state.run(&op, &jacobi, tol, settings.max_linear_iterations);
    checkpoint(cx)?;
    let mut solution = state.x;
    let mut rel = true_relative(&solution);
    let gate = 1e3 * tol.max(1e-14);
    let remaining = settings.max_linear_iterations.saturating_sub(report.iters);
    // MINRES stops on an M-norm recurrence estimate, which need not meet this
    // owner's independently recomputed Euclidean gate. On a nearly converged
    // ill-conditioned system, solve ONE correction equation from the true
    // residual. Restarting the already-converged recurrence would do no work.
    // Charge every correction iteration to the SAME MINRES-path allowance.
    if rel.is_finite()
        && rel > gate
        && remaining > 0
        && solution.iter().all(|value| value.is_finite())
    {
        let mut residual = vec![0.0; rhs.len()];
        kff.spmv(&solution, &mut residual);
        for (value, &source) in residual.iter_mut().zip(rhs) {
            *value = source - *value;
        }
        checkpoint(cx)?;
        let mut weighted = vec![0.0; rhs.len()];
        jacobi.apply(&residual, &mut weighted);
        let quadratic_form = fs_solver::dot(&residual, &weighted);
        // A finite Euclidean residual may still underflow in the M-norm.
        // Refuse normally instead of triggering the constructor's assertion.
        if quadratic_form.is_finite() && quadratic_form > 0.0 {
            let mut correction = PminresState::new(&op, &jacobi, &residual);
            let _ = correction.run(&op, &jacobi, tol, remaining);
            checkpoint(cx)?;
            let candidate: Vec<f64> = solution
                .iter()
                .zip(&correction.x)
                .map(|(value, change)| value + change)
                .collect();
            // The correction RHS has a different norm. Acceptance must still
            // use the ORIGINAL rhs, never the correction's residual report.
            let improved = true_relative(&candidate);
            if improved.is_finite()
                && improved < rel
                && candidate.iter().all(|value| value.is_finite())
            {
                solution = candidate;
                rel = improved;
            }
        }
    }
    checkpoint(cx)?;
    if rel.is_finite() && rel <= gate && solution.iter().all(|value| value.is_finite()) {
        Ok((solution, true))
    } else {
        Err(SmallStrainTetError::LinearSolveFailed {
            increment,
            relative_residual: rel,
        })
    }
}

struct AbsJacobi {
    inv: Vec<f64>,
}

impl AbsJacobi {
    fn new(matrix: &Csr) -> Self {
        Self {
            inv: (0..matrix.nrows())
                .map(|i| {
                    let diagonal = matrix.get(i, i).abs();
                    if diagonal > 0.0 { 1.0 / diagonal } else { 1.0 }
                })
                .collect(),
        }
    }
}

enum GeneralPreconditioner {
    Ilu(Ilu0),
    Jacobi(AbsJacobi),
}

impl Precond for GeneralPreconditioner {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        match self {
            Self::Ilu(ilu) => ilu.apply(r, z),
            Self::Jacobi(jacobi) => jacobi.apply(r, z),
        }
    }
}

impl Precond for AbsJacobi {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        for (zi, (ri, di)) in z.iter_mut().zip(r.iter().zip(&self.inv)) {
            *zi = ri * di;
        }
    }
}
