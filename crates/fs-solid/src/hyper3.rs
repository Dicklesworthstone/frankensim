//! Three-dimensional finite-strain hyperelasticity on body-fitted linear
//! tetrahedra (plan §8.2): Total-Lagrangian kinematics, the fs-material
//! energy cards (exact AD first Piola–Kirchhoff stress and consistent 9x9
//! tangent), load stepping, and a globalized Newton method.
//!
//! # Discretization
//!
//! On each P1 tetrahedron the reference shape-function gradients `g_a`
//! (a = 0..4) are constant, so the deformation gradient
//! `F = I + sum_a u_a (x) g_a` is constant per element and one stress /
//! tangent evaluation is exact quadrature. With `V` the reference volume,
//!
//! ```text
//! internal force   r_(a,i)       = V sum_J P_iJ g_a,J
//! tangent          K_(a,i),(b,k) = V sum_(J,L) A_iJkL g_a,J g_b,L
//! potential        Pi(u) = sum_e V W(F_e) - f_ext . u
//! ```
//!
//! External loads are DEAD (referential): nodal forces and a body force per
//! reference volume, lumped equally to the four vertices. Prescribed
//! displacements and loads ramp linearly over `load_steps`.
//!
//! # Newton globalization
//!
//! Each step solves `K_ff d_f = -r_f - K_fp dp` where `dp` is the prescribed
//! increment still pending, i.e. the boundary motion enters through the
//! consistent linearization instead of a distorting jump in one element
//! layer. A trial state that inverts any element (`det F <= 0`, refused by
//! the material card) halves the step. While a prescribed increment is
//! pending the first admissible step is taken; once the boundary data are
//! fully applied an Armijo condition on the potential is enforced (or, on a
//! non-descent direction near a limit point, a decrease of the residual
//! norm). Convergence: `||r_f|| <= relative_tolerance * max(||f_int||,
//! ||f_ext||, max|K_dd| max|u|)` with the full internal-force vector
//! (reactions included) and no pending boundary increment; the stiffness
//! term keeps the gate meaningful at stress-free states such as a rigid
//! rotation, where every force vanishes.
//!
//! The tangent is SPD for stable states and is solved with ILU(0)-PCG; an
//! indefinite tangent (detected by a stalled CG) falls back to PMINRES with
//! an absolute-diagonal Jacobi preconditioner. The accepted linear residual
//! is the recomputed true residual.
//!
//! # No-claim boundaries
//!
//! Static, isothermal, rate-independent hyperelasticity on conforming P1
//! tetrahedra (constant strain: volumetric locking for nearly incompressible
//! cards is NOT mitigated; no mixed/F-bar formulation). Dead loads only (no
//! follower pressure). No contact, no limit-point continuation (a solve that
//! passes a limit point refuses or must be driven by displacement), no
//! mesh-convergence claim from a single run. Determinism: sequential
//! assembly and solves in canonical order, bit-reproducible on one
//! ISA/toolchain profile.

use fs_exec::Cx;
use fs_material::hyper::Hyperelastic;
use fs_solver::krylov::PminresState;
use fs_solver::op::CsrOp;
use fs_sparse::precond::{Precond, ilu0, pcg};
use fs_sparse::{Coo, Csr};

use crate::linear3::TetAssemblyBudget;

/// Newton / load-stepping controls.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HyperTetSettings {
    /// Equal load increments from zero to the full load.
    pub load_steps: usize,
    /// Newton iterations allowed per load step.
    pub max_newton_iterations: usize,
    /// Relative force-residual gate (see module docs).
    pub relative_tolerance: f64,
    /// Step halvings allowed per Newton iteration.
    pub max_backtracks: usize,
    /// Relative residual required of each linear solve.
    pub linear_tolerance: f64,
    /// Iteration cap of each linear solve.
    pub max_linear_iterations: usize,
}

impl Default for HyperTetSettings {
    fn default() -> Self {
        Self {
            load_steps: 1,
            max_newton_iterations: 40,
            relative_tolerance: 1e-10,
            max_backtracks: 30,
            linear_tolerance: 1e-12,
            max_linear_iterations: 20_000,
        }
    }
}

/// A static finite-strain problem on a tetrahedral mesh.
#[derive(Debug, Clone, Copy)]
pub struct HyperTetProblem<'a> {
    /// Reference vertex coordinates, m.
    pub nodes_m: &'a [[f64; 3]],
    /// Four vertex indices per conforming tetrahedron.
    pub tetrahedra: &'a [[usize; 4]],
    /// The hyperelastic card (uniform over the mesh).
    pub material: &'a Hyperelastic,
    /// Prescribed displacements at full load: (`3 node + component`, m).
    pub prescribed_m: &'a [(usize, f64)],
    /// Dead nodal forces at full load: (`3 node + component`, N).
    pub nodal_forces_n: &'a [(usize, f64)],
    /// Dead body force per reference volume at full load, N/m^3.
    pub body_force_n_m3: [f64; 3],
    /// Size and element-quality envelope.
    pub budget: TetAssemblyBudget,
    /// Newton controls.
    pub settings: HyperTetSettings,
}

/// Evidence of one load step.
#[derive(Debug, Clone, PartialEq)]
pub struct HyperLoadStep {
    /// Load factor in `(0, 1]`.
    pub load_factor: f64,
    /// Relative residual after each Newton iteration (first entry: before
    /// the first update).
    pub residual_history: Vec<f64>,
    /// Total step halvings.
    pub backtracks: usize,
    /// Total linear-solver iterations.
    pub linear_iterations: usize,
    /// Linear solves that fell back to PMINRES (indefinite tangent).
    pub indefinite_fallbacks: usize,
}

/// Converged finite-strain state.
#[derive(Debug, Clone, PartialEq)]
pub struct HyperTetSolution {
    /// Nodal displacements, m.
    pub displacement_m: Vec<[f64; 3]>,
    /// Reaction forces at the prescribed DOFs (internal minus external), N,
    /// in the order of `prescribed_m`.
    pub reactions_n: Vec<(usize, f64)>,
    /// Stored strain energy `sum V W(F)`, J.
    pub strain_energy_j: f64,
    /// Smallest element `det F` in the converged state.
    pub min_det_f: f64,
    /// Per-load-step evidence.
    pub steps: Vec<HyperLoadStep>,
}

/// Structured refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum HyperTetError {
    /// Inconsistent, out-of-budget, or non-finite input.
    InvalidInput {
        /// What was refused.
        what: String,
    },
    /// A reference tetrahedron is degenerate or inverted.
    DegenerateElement {
        /// Element index.
        element: usize,
        /// `det(J) / h_max^3` of the reference element.
        scaled_jacobian: f64,
    },
    /// The material refused the state of an element (in the converged or
    /// starting state; trial states are absorbed by step halving).
    MaterialRefused {
        /// Element index.
        element: usize,
        /// The card's message.
        what: String,
    },
    /// Newton did not reach the gate. Repair: more load steps.
    NewtonStalled {
        /// Load factor of the failing step.
        load_factor: f64,
        /// Relative residual history of that step.
        history: Vec<f64>,
    },
    /// A linear solve missed its gate on both the PCG and PMINRES paths.
    LinearSolveFailed {
        /// Load factor of the failing step.
        load_factor: f64,
        /// Recomputed relative residual reached.
        relative_residual: f64,
    },
    /// The context was cancelled.
    Cancelled,
}

impl core::fmt::Display for HyperTetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput { what } => write!(f, "FS-SOLID-HYPER3-INPUT: {what}"),
            Self::DegenerateElement {
                element,
                scaled_jacobian,
            } => write!(
                f,
                "FS-SOLID-HYPER3-DEGENERATE: element {element} scaled Jacobian {scaled_jacobian:e}"
            ),
            Self::MaterialRefused { element, what } => {
                write!(f, "FS-SOLID-HYPER3-MATERIAL: element {element}: {what}")
            }
            Self::NewtonStalled {
                load_factor,
                history,
            } => write!(
                f,
                "FS-SOLID-HYPER3-NEWTON: load factor {load_factor} stalled after {} iterations (last relative residual {:e}); add load steps",
                history.len(),
                history.last().copied().unwrap_or(f64::NAN)
            ),
            Self::LinearSolveFailed {
                load_factor,
                relative_residual,
            } => write!(
                f,
                "FS-SOLID-HYPER3-LINEAR: load factor {load_factor}: linear residual {relative_residual:e}"
            ),
            Self::Cancelled => write!(f, "FS-SOLID-HYPER3-CANCELLED"),
        }
    }
}

impl std::error::Error for HyperTetError {}

/// Reference geometry of one element.
#[derive(Debug, Clone, Copy)]
struct Element {
    nodes: [usize; 4],
    gradients: [[f64; 3]; 4],
    volume: f64,
}

/// Per-element state for one displacement field.
struct Evaluated {
    internal: Vec<f64>,
    energy: f64,
    min_det_f: f64,
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).fold(0.0, |s, (x, y)| x.mul_add(*y, s))
}

fn norm(a: &[f64]) -> f64 {
    fs_math::det::sqrt(dot(a, a))
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), HyperTetError> {
    cx.checkpoint().map_err(|_| HyperTetError::Cancelled)
}

/// `(a, b, c) -> (a x b) . c`.
fn triple(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    (a[1] * b[2] - a[2] * b[1]) * c[0]
        + (a[2] * b[0] - a[0] * b[2]) * c[1]
        + (a[0] * b[1] - a[1] * b[0]) * c[2]
}

impl HyperTetProblem<'_> {
    fn validate(&self) -> Result<Vec<Element>, HyperTetError> {
        let invalid = |what: String| Err(HyperTetError::InvalidInput { what });
        let n = self.nodes_m.len();
        if n == 0 || self.tetrahedra.is_empty() {
            return invalid("empty mesh".into());
        }
        if n > self.budget.maximum_nodes || self.tetrahedra.len() > self.budget.maximum_tetrahedra {
            return invalid(format!(
                "{n} nodes / {} tetrahedra exceed the budget {:?}",
                self.tetrahedra.len(),
                self.budget
            ));
        }
        if self.nodes_m.iter().flatten().any(|v| !v.is_finite()) {
            return invalid("non-finite node coordinate".into());
        }
        let s = self.settings;
        if s.load_steps == 0
            || s.max_newton_iterations == 0
            || !(s.relative_tolerance.is_finite() && s.relative_tolerance > 0.0)
            || !(s.linear_tolerance.is_finite() && s.linear_tolerance > 0.0)
        {
            return invalid(format!("inadmissible settings {s:?}"));
        }
        let dofs = 3 * n;
        let mut prescribed = vec![false; dofs];
        for &(dof, value) in self.prescribed_m {
            if dof >= dofs || !value.is_finite() || prescribed[dof] {
                return invalid(format!(
                    "prescribed DOF {dof} = {value} is out of range, non-finite, or repeated"
                ));
            }
            prescribed[dof] = true;
        }
        for &(dof, value) in self.nodal_forces_n {
            if dof >= dofs || !value.is_finite() {
                return invalid(format!(
                    "nodal force on DOF {dof} = {value} is out of range or non-finite"
                ));
            }
        }
        if self.body_force_n_m3.iter().any(|v| !v.is_finite()) {
            return invalid("non-finite body force".into());
        }
        let free = dofs - self.prescribed_m.len();
        if free > self.budget.maximum_free_dofs {
            return invalid(format!("{free} free DOFs exceed the budget"));
        }
        let mut elements = Vec::with_capacity(self.tetrahedra.len());
        for (element, tet) in self.tetrahedra.iter().enumerate() {
            if tet.iter().any(|&v| v >= n) {
                return invalid(format!("element {element} names a missing vertex"));
            }
            let x = tet.map(|v| self.nodes_m[v]);
            let e = |k: usize| [x[k][0] - x[0][0], x[k][1] - x[0][1], x[k][2] - x[0][2]];
            let (e1, e2, e3) = (e(1), e(2), e(3));
            let det = triple(e1, e2, e3);
            let mut h_max = 0.0f64;
            for a in 0..4 {
                for b in a + 1..4 {
                    let d = [x[b][0] - x[a][0], x[b][1] - x[a][1], x[b][2] - x[a][2]];
                    h_max = h_max.max(fs_math::det::sqrt(dot(&d, &d)));
                }
            }
            let scaled = det / (h_max * h_max * h_max);
            if !(scaled.is_finite() && scaled >= self.budget.minimum_scaled_jacobian) {
                return Err(HyperTetError::DegenerateElement {
                    element,
                    scaled_jacobian: scaled,
                });
            }
            // Rows of D^-1 (D = [e1 e2 e3] as columns) are the gradients of
            // the barycentric coordinates of vertices 1..3.
            let cross = |a: [f64; 3], b: [f64; 3]| {
                [
                    a[1] * b[2] - a[2] * b[1],
                    a[2] * b[0] - a[0] * b[2],
                    a[0] * b[1] - a[1] * b[0],
                ]
            };
            let g1 = cross(e2, e3).map(|v| v / det);
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

    fn deformation_gradient(element: &Element, u: &[f64]) -> [f64; 9] {
        let mut f = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        for (a, &node) in element.nodes.iter().enumerate() {
            for i in 0..3 {
                let ui = u[3 * node + i];
                for j in 0..3 {
                    f[3 * i + j] = ui.mul_add(element.gradients[a][j], f[3 * i + j]);
                }
            }
        }
        f
    }

    /// Internal force (all DOFs) and stored energy; `Err(element)` when the
    /// card refuses a state.
    fn evaluate(
        &self,
        elements: &[Element],
        u: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Result<Evaluated, (usize, String)>, HyperTetError> {
        let mut internal = vec![0.0f64; u.len()];
        let mut energy = 0.0f64;
        let mut min_det_f = f64::INFINITY;
        for (index, element) in elements.iter().enumerate() {
            if index % 4096 == 0 {
                checkpoint(cx)?;
            }
            let f = Self::deformation_gradient(element, u);
            let j = f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
                + f[2] * (f[3] * f[7] - f[4] * f[6]);
            min_det_f = min_det_f.min(j);
            let p = match self.material.piola(&f) {
                Ok(p) => p,
                Err(error) => return Ok(Err((index, error.to_string()))),
            };
            let w = self.material.energy(&f);
            if !w.is_finite() {
                return Ok(Err((
                    index,
                    format!("non-finite stored energy at det F = {j}"),
                )));
            }
            energy = element.volume.mul_add(w, energy);
            for (a, &node) in element.nodes.iter().enumerate() {
                for i in 0..3 {
                    let mut s = 0.0;
                    for jj in 0..3 {
                        s = p[3 * i + jj].mul_add(element.gradients[a][jj], s);
                    }
                    internal[3 * node + i] = element.volume.mul_add(s, internal[3 * node + i]);
                }
            }
        }
        Ok(Ok(Evaluated {
            internal,
            energy,
            min_det_f,
        }))
    }

    fn tangent(
        &self,
        elements: &[Element],
        u: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Result<Csr, (usize, String)>, HyperTetError> {
        let dofs = u.len();
        let mut coo = Coo::new(dofs, dofs);
        for (index, element) in elements.iter().enumerate() {
            if index % 4096 == 0 {
                checkpoint(cx)?;
            }
            let f = Self::deformation_gradient(element, u);
            let a4 = match self.material.tangent(&f) {
                Ok(a4) => a4,
                Err(error) => return Ok(Err((index, error.to_string()))),
            };
            let g = &element.gradients;
            for a in 0..4 {
                for i in 0..3 {
                    let row = 3 * element.nodes[a] + i;
                    for b in 0..4 {
                        for k in 0..3 {
                            let mut s = 0.0;
                            for jj in 0..3 {
                                for l in 0..3 {
                                    s = (a4[3 * i + jj][3 * k + l] * g[a][jj]).mul_add(g[b][l], s);
                                }
                            }
                            coo.push(row, 3 * element.nodes[b] + k, element.volume * s);
                        }
                    }
                }
            }
        }
        Ok(Ok(coo.assemble()))
    }

    fn external(&self, elements: &[Element], load: f64) -> Vec<f64> {
        let mut f = vec![0.0f64; 3 * self.nodes_m.len()];
        for &(dof, value) in self.nodal_forces_n {
            f[dof] += load * value;
        }
        if self.body_force_n_m3.iter().any(|&b| b != 0.0) {
            for element in elements {
                for &node in &element.nodes {
                    for i in 0..3 {
                        f[3 * node + i] += load * self.body_force_n_m3[i] * element.volume / 4.0;
                    }
                }
            }
        }
        f
    }

    /// Internal force vector `sum_e V P : grad N` over all DOFs (node-major
    /// xyz) at a full displacement field — the exact gradient of the stored
    /// energy, exposed for verification.
    ///
    /// # Errors
    /// Input/budget refusals, [`HyperTetError::MaterialRefused`].
    pub fn internal_force(
        &self,
        displacement_m: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, HyperTetError> {
        let elements = self.validate()?;
        self.check_len(displacement_m)?;
        match self.evaluate(&elements, displacement_m, cx)? {
            Ok(evaluated) => Ok(evaluated.internal),
            Err((element, what)) => Err(HyperTetError::MaterialRefused { element, what }),
        }
    }

    /// Assembled consistent tangent over all DOFs at a full displacement
    /// field, exposed for verification.
    ///
    /// # Errors
    /// Input/budget refusals, [`HyperTetError::MaterialRefused`].
    pub fn tangent_matrix(
        &self,
        displacement_m: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Csr, HyperTetError> {
        let elements = self.validate()?;
        self.check_len(displacement_m)?;
        match self.tangent(&elements, displacement_m, cx)? {
            Ok(k) => Ok(k),
            Err((element, what)) => Err(HyperTetError::MaterialRefused { element, what }),
        }
    }

    fn check_len(&self, u: &[f64]) -> Result<(), HyperTetError> {
        if u.len() == 3 * self.nodes_m.len() && u.iter().all(|v| v.is_finite()) {
            Ok(())
        } else {
            Err(HyperTetError::InvalidInput {
                what: "displacement must hold 3 finite values per node".into(),
            })
        }
    }

    /// Solve by load stepping and globalized Newton.
    ///
    /// # Errors
    /// See [`HyperTetError`].
    #[allow(clippy::too_many_lines)] // one load-step / Newton / line-search loop
    pub fn solve(&self, cx: &Cx<'_>) -> Result<HyperTetSolution, HyperTetError> {
        let elements = self.validate()?;
        let n = self.nodes_m.len();
        let dofs = 3 * n;
        let mut is_prescribed = vec![false; dofs];
        for &(dof, _) in self.prescribed_m {
            is_prescribed[dof] = true;
        }
        let free: Vec<usize> = (0..dofs).filter(|&d| !is_prescribed[d]).collect();
        let mut row_of = vec![usize::MAX; dofs];
        for (row, &dof) in free.iter().enumerate() {
            row_of[dof] = row;
        }
        let mut u = vec![0.0f64; dofs];
        let mut steps = Vec::with_capacity(self.settings.load_steps);
        let refused =
            |(element, what): (usize, String)| HyperTetError::MaterialRefused { element, what };
        let mut state = self.evaluate(&elements, &u, cx)?.map_err(refused)?;
        for step in 1..=self.settings.load_steps {
            let load = step as f64 / self.settings.load_steps as f64;
            let f_ext = self.external(&elements, load);
            let target: Vec<(usize, f64)> = self
                .prescribed_m
                .iter()
                .map(|&(dof, value)| (dof, load * value))
                .collect();
            let mut record = HyperLoadStep {
                load_factor: load,
                residual_history: Vec::new(),
                backtracks: 0,
                linear_iterations: 0,
                indefinite_fallbacks: 0,
            };
            let mut converged = false;
            for _ in 0..self.settings.max_newton_iterations {
                checkpoint(cx)?;
                let residual: Vec<f64> =
                    free.iter().map(|&d| state.internal[d] - f_ext[d]).collect();
                let pending: Vec<(usize, f64)> = target
                    .iter()
                    .map(|&(dof, value)| (dof, value - u[dof]))
                    .collect();
                let pending_any = pending.iter().any(|&(_, dp)| dp != 0.0);
                // Tangent first: its diagonal times the displacement
                // magnitude is the force scale of the current deformation,
                // which stays meaningful at stress-free states (a rigid
                // motion) where the internal force itself vanishes.
                let full = self.tangent(&elements, &u, cx)?.map_err(refused)?;
                let stiffness = (0..dofs).map(|d| full.get(d, d).abs()).fold(0.0, f64::max);
                let reach = u.iter().map(|v| v.abs()).fold(0.0, f64::max);
                let scale = norm(&state.internal)
                    .max(norm(&f_ext))
                    .max(stiffness * reach)
                    .max(f64::MIN_POSITIVE);
                let relative = norm(&residual) / scale;
                record.residual_history.push(relative);
                if !pending_any && relative <= self.settings.relative_tolerance {
                    converged = true;
                    break;
                }
                // K_ff d = -r_f - K_fp dp.
                let mut coo = Coo::new(free.len(), free.len());
                let mut rhs: Vec<f64> = residual.iter().map(|r| -r).collect();
                let mut dp_full = vec![0.0f64; dofs];
                for &(dof, dp) in &pending {
                    dp_full[dof] = dp;
                }
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
                let kff = coo.assemble();
                let (direction, iterations, fallback) = self.linear_solve(&kff, &rhs, load)?;
                record.linear_iterations += iterations;
                record.indefinite_fallbacks += usize::from(fallback);
                // Line search on the joint (free, prescribed) update.
                let slope = dot(&residual, &direction);
                let residual_norm = norm(&residual);
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
                    if let Ok(evaluated) = self.evaluate(&elements, &trial, cx)? {
                        let acceptable = if pending_any {
                            true
                        } else if slope < 0.0 {
                            let potential = |e: &Evaluated, x: &[f64]| e.energy - dot(&f_ext, x);
                            let (before, after) =
                                (potential(&state, &u), potential(&evaluated, &trial));
                            after
                                <= (1e-4 * alpha).mul_add(slope, before)
                                    + 1e-12 * before.abs().max(after.abs())
                        } else {
                            let trial_residual: Vec<f64> = free
                                .iter()
                                .map(|&d| evaluated.internal[d] - f_ext[d])
                                .collect();
                            norm(&trial_residual) < residual_norm
                        };
                        if acceptable {
                            accepted = Some((trial, evaluated));
                            break;
                        }
                    }
                    alpha *= 0.5;
                    record.backtracks += 1;
                }
                let Some((trial, evaluated)) = accepted else {
                    return Err(HyperTetError::NewtonStalled {
                        load_factor: load,
                        history: record.residual_history,
                    });
                };
                u = trial;
                state = evaluated;
            }
            if !converged {
                return Err(HyperTetError::NewtonStalled {
                    load_factor: load,
                    history: record.residual_history,
                });
            }
            steps.push(record);
        }
        let f_ext = self.external(&elements, 1.0);
        let reactions_n = self
            .prescribed_m
            .iter()
            .map(|&(dof, _)| (dof, state.internal[dof] - f_ext[dof]))
            .collect();
        Ok(HyperTetSolution {
            displacement_m: u.as_chunks::<3>().0.to_vec(),
            reactions_n,
            strain_energy_j: state.energy,
            min_det_f: state.min_det_f,
            steps,
        })
    }

    /// PCG with ILU(0); on a stalled or non-finite CG (indefinite tangent),
    /// PMINRES with |diag| Jacobi. Returns (solution, iterations, fallback).
    fn linear_solve(
        &self,
        kff: &Csr,
        rhs: &[f64],
        load: f64,
    ) -> Result<(Vec<f64>, usize, bool), HyperTetError> {
        let tol = self.settings.linear_tolerance;
        let true_relative = |x: &[f64]| {
            let mut ax = vec![0.0f64; rhs.len()];
            kff.spmv(x, &mut ax);
            let r: Vec<f64> = rhs.iter().zip(&ax).map(|(b, a)| b - a).collect();
            norm(&r) / norm(rhs).max(f64::MIN_POSITIVE)
        };
        if norm(rhs) == 0.0 {
            return Ok((vec![0.0; rhs.len()], 0, false));
        }
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
            if report.converged && x.iter().all(|v| v.is_finite()) {
                let rel = true_relative(&x);
                if rel <= 10.0 * tol {
                    return Ok((x, report.iters, false));
                }
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
            Ok((state.x.clone(), state.iters, true))
        } else {
            Err(HyperTetError::LinearSolveFailed {
                load_factor: load,
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
