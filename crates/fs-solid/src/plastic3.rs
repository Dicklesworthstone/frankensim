//! Incremental small-strain inelasticity on body-fitted linear tetrahedra
//! (plan §8.2 "J2 plasticity with return mapping"), generic over any
//! `fs-material` [`SmallStrainLaw`]: the law's algorithmic stress update and
//! consistent tangent drive a load-path-following Newton solve, and the
//! internal variables are committed only at converged increments.
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
//! `||r_f|| <= relative_tolerance * max(||f_int||, ||f_ext||,
//! max|K_dd| max|u|)`. After convergence each element state is updated by
//! the law (`update_state`) and committed.
//!
//! # No-claim boundaries
//!
//! Small strain and rotation (no finite-strain plasticity), rate
//! independence, constant-strain tetrahedra (plastic incompressibility locks
//! P1 tetrahedra in confined flow: limit loads of confined problems are
//! over-predicted; no mixed/B-bar formulation), dead loads only, no
//! limit-load continuation (an increment past a limit load stalls and
//! refuses; drive by displacement), no mesh-convergence claim. Linear solves
//! as in `hyper3` (ILU(0)-PCG, |diag|-Jacobi PMINRES fallback, recomputed
//! residual acceptance). Deterministic sequential assembly.

use fs_exec::Cx;
use fs_material::{SmallStrainLaw, Tangent6, Voigt};
use fs_solver::krylov::PminresState;
use fs_solver::op::CsrOp;
use fs_sparse::precond::{Precond, ilu0, pcg};
use fs_sparse::{Coo, Csr};

use crate::linear3::TetAssemblyBudget;

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
    /// Iteration cap of each linear solve.
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

/// A load-path problem on a tetrahedral mesh.
#[derive(Debug, Clone, Copy)]
pub struct SmallStrainTetProblem<'a, L: SmallStrainLaw> {
    /// Vertex coordinates, m.
    pub nodes_m: &'a [[f64; 3]],
    /// Four vertex indices per conforming tetrahedron.
    pub tetrahedra: &'a [[usize; 4]],
    /// The constitutive law (uniform over the mesh).
    pub law: &'a L,
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
        for (a, &node) in element.nodes.iter().enumerate() {
            for i in 0..3 {
                let col = b_column(element.gradients[a], i);
                let ui = u[3 * node + i];
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
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, SmallStrainTetError> {
        let mut f = vec![0.0f64; u.len()];
        for (index, element) in elements.iter().enumerate() {
            if index % 4096 == 0 {
                checkpoint(cx)?;
            }
            let sigma = self.law.stress(&Self::strain(element, u), &states[index]);
            for (a, &node) in element.nodes.iter().enumerate() {
                for i in 0..3 {
                    let col = b_column(element.gradients[a], i);
                    let mut s = 0.0;
                    for v in 0..6 {
                        s = (col[v] * WEIGHT[v]).mul_add(sigma[v], s);
                    }
                    f[3 * node + i] = element.volume.mul_add(s, f[3 * node + i]);
                }
            }
        }
        Ok(f)
    }

    fn tangent(
        &self,
        elements: &[Element],
        states: &[L::State],
        u: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Csr, SmallStrainTetError> {
        let mut coo = Coo::new(u.len(), u.len());
        for (index, element) in elements.iter().enumerate() {
            if index % 4096 == 0 {
                checkpoint(cx)?;
            }
            let c: Tangent6 = self.law.tangent(&Self::strain(element, u), &states[index]);
            let cols: Vec<Voigt> = (0..12)
                .map(|k| b_column(element.gradients[k / 3], k % 3))
                .collect();
            for (p, bp) in cols.iter().enumerate() {
                let row = 3 * element.nodes[p / 3] + p % 3;
                for (q, bq) in cols.iter().enumerate() {
                    let mut s = 0.0;
                    for v in 0..6 {
                        let mut cb = 0.0;
                        for w in 0..6 {
                            cb = c[v][w].mul_add(bq[w], cb);
                        }
                        s = (bp[v] * WEIGHT[v]).mul_add(cb, s);
                    }
                    coo.push(row, 3 * element.nodes[q / 3] + q % 3, element.volume * s);
                }
            }
        }
        Ok(coo.assemble())
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
        self.internal(&elements, states, displacement_m, cx)
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
        self.tangent(&elements, states, displacement_m, cx)
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
    #[allow(clippy::too_many_lines)] // one increment / Newton / backtracking loop
    pub fn solve(
        &self,
        cx: &Cx<'_>,
    ) -> Result<SmallStrainTetSolution<L::State>, SmallStrainTetError> {
        let elements = self.validate()?;
        let dofs = 3 * self.nodes_m.len();
        let mut is_prescribed = vec![false; dofs];
        for &(dof, _) in self.prescribed_m {
            is_prescribed[dof] = true;
        }
        let free: Vec<usize> = (0..dofs).filter(|&d| !is_prescribed[d]).collect();
        let mut row_of = vec![usize::MAX; dofs];
        for (row, &dof) in free.iter().enumerate() {
            row_of[dof] = row;
        }
        let mut states: Vec<L::State> = vec![self.law.initial_state(); elements.len()];
        let mut u = vec![0.0f64; dofs];
        let mut increments = Vec::with_capacity(self.load_path.len());
        for (increment, &load) in self.load_path.iter().enumerate() {
            let mut f_ext = vec![0.0f64; dofs];
            for &(dof, value) in self.nodal_forces_n {
                f_ext[dof] += load * value;
            }
            let target: Vec<(usize, f64)> = self
                .prescribed_m
                .iter()
                .map(|&(dof, value)| (dof, load * value))
                .collect();
            let mut history = Vec::new();
            let mut fallbacks = 0usize;
            let mut converged = false;
            let mut f_int = self.internal(&elements, &states, &u, cx)?;
            for _ in 0..self.settings.max_newton_iterations {
                checkpoint(cx)?;
                let full = self.tangent(&elements, &states, &u, cx)?;
                let residual: Vec<f64> = free.iter().map(|&d| f_int[d] - f_ext[d]).collect();
                let pending: Vec<(usize, f64)> = target
                    .iter()
                    .map(|&(dof, value)| (dof, value - u[dof]))
                    .collect();
                let pending_any = pending.iter().any(|&(_, dp)| dp != 0.0);
                let stiffness = (0..dofs).map(|d| full.get(d, d).abs()).fold(0.0, f64::max);
                let reach = u.iter().map(|v| v.abs()).fold(0.0, f64::max);
                let scale = norm(&f_int)
                    .max(norm(&f_ext))
                    .max(stiffness * reach)
                    .max(f64::MIN_POSITIVE);
                let relative = norm(&residual) / scale;
                history.push(relative);
                if !pending_any && relative <= self.settings.relative_tolerance {
                    converged = true;
                    break;
                }
                let mut dp_full = vec![0.0f64; dofs];
                for &(dof, dp) in &pending {
                    dp_full[dof] = dp;
                }
                let mut coo = Coo::new(free.len(), free.len());
                let mut rhs: Vec<f64> = residual.iter().map(|r| -r).collect();
                for (row, &dof) in free.iter().enumerate() {
                    let (cols, vals) = full.row(dof);
                    for (&c, &v) in cols.iter().zip(vals) {
                        if is_prescribed[c] {
                            rhs[row] = (-v).mul_add(dp_full[c], rhs[row]);
                        } else {
                            coo.push(row, row_of[c], v);
                        }
                    }
                }
                let (direction, fallback) = self.linear_solve(&coo.assemble(), &rhs, increment)?;
                fallbacks += usize::from(fallback);
                let before = norm(&residual);
                let mut alpha = 1.0f64;
                let mut accepted = None;
                for _ in 0..=self.settings.max_backtracks {
                    let mut trial = u.clone();
                    for (row, &dof) in free.iter().enumerate() {
                        trial[dof] = alpha.mul_add(direction[row], trial[dof]);
                    }
                    for &(dof, dp) in &pending {
                        trial[dof] = alpha.mul_add(dp, trial[dof]);
                    }
                    let f_trial = self.internal(&elements, &states, &trial, cx)?;
                    let after: Vec<f64> = free.iter().map(|&d| f_trial[d] - f_ext[d]).collect();
                    if pending_any || norm(&after) < before || alpha < 1e-3 {
                        accepted = Some((trial, f_trial));
                        break;
                    }
                    alpha *= 0.5;
                }
                let Some((trial, f_trial)) = accepted else {
                    return Err(SmallStrainTetError::NewtonStalled { increment, history });
                };
                u = trial;
                f_int = f_trial;
            }
            if !converged {
                return Err(SmallStrainTetError::NewtonStalled { increment, history });
            }
            let mut evolving = 0usize;
            for (index, element) in elements.iter().enumerate() {
                let next = self
                    .law
                    .update_state(&Self::strain(element, &u), &states[index]);
                evolving += usize::from(next != states[index]);
                states[index] = next;
            }
            increments.push(IncrementRecord {
                load_factor: load,
                residual_history: history,
                reactions_n: self
                    .prescribed_m
                    .iter()
                    .map(|&(dof, _)| f_int[dof] - f_ext[dof])
                    .collect(),
                evolving_elements: evolving,
                indefinite_fallbacks: fallbacks,
            });
        }
        let element_strain: Vec<Voigt> = elements.iter().map(|e| Self::strain(e, &u)).collect();
        let element_stress = element_strain
            .iter()
            .zip(&states)
            .map(|(eps, state)| self.law.stress(eps, state))
            .collect();
        Ok(SmallStrainTetSolution {
            displacement_m: u.as_chunks::<3>().0.to_vec(),
            element_strain,
            element_stress,
            element_states: states,
            increments,
        })
    }

    fn linear_solve(
        &self,
        kff: &Csr,
        rhs: &[f64],
        increment: usize,
    ) -> Result<(Vec<f64>, bool), SmallStrainTetError> {
        let tol = self.settings.linear_tolerance;
        if norm(rhs) == 0.0 {
            return Ok((vec![0.0; rhs.len()], false));
        }
        let true_relative = |x: &[f64]| {
            let mut ax = vec![0.0f64; rhs.len()];
            kff.spmv(x, &mut ax);
            let r: Vec<f64> = rhs.iter().zip(&ax).map(|(b, a)| b - a).collect();
            norm(&r) / norm(rhs)
        };
        if let Ok(m) = ilu0(kff) {
            let mut x = vec![0.0f64; rhs.len()];
            let report = pcg(
                kff,
                rhs,
                &mut x,
                &m,
                tol,
                self.settings.max_linear_iterations,
            );
            if report.converged
                && x.iter().all(|v| v.is_finite())
                && true_relative(&x) <= 10.0 * tol
            {
                return Ok((x, false));
            }
        }
        let jacobi = AbsJacobi {
            inv: (0..kff.nrows())
                .map(|i| {
                    let d = kff.get(i, i).abs();
                    if d > 0.0 { 1.0 / d } else { 1.0 }
                })
                .collect(),
        };
        let op = CsrOp::symmetric(kff.clone());
        let mut state = PminresState::new(&op, &jacobi, rhs);
        let _ = state.run(&op, &jacobi, tol, self.settings.max_linear_iterations);
        let rel = true_relative(&state.x);
        if rel.is_finite() && rel <= 1e3 * tol.max(1e-14) {
            Ok((state.x.clone(), true))
        } else {
            Err(SmallStrainTetError::LinearSolveFailed {
                increment,
                relative_residual: rel,
            })
        }
    }
}

struct AbsJacobi {
    inv: Vec<f64>,
}

impl Precond for AbsJacobi {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        for (zi, (ri, di)) in z.iter_mut().zip(r.iter().zip(&self.inv)) {
            *zi = ri * di;
        }
    }
}
