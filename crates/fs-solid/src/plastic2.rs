//! Incremental plane-stress and plane-strain inelasticity on the existing
//! [`Mesh2`] P1/Q1 elements. The three-dimensional [`SmallStrainLaw`] remains
//! the constitutive owner: plane stress uses its checked transverse reduction
//! and Schur-complement tangent; plane strain retains the constrained stress
//! and every component of the original material history.
//!
//! Each `mesh2::quad_points` integration point has its own committed history.
//! Newton and line-search proposals always evaluate from the OLD history;
//! only a converged physical load entry advances it. A timed law therefore
//! advances its declared interval once per entry, including repeated load
//! factors. Endpoint stresses are retained from that same update.
//!
//! Work uses tensor shear: `B = [epsilon_xx, epsilon_yy, epsilon_zz, epsilon_xy]`
//! and `W = diag(1, 1, 1, 2)`, so `f = integral B^T W sigma` and
//! `K = integral B^T W C_alg B`. Positive reference thickness converts
//! stresses into forces; geometry and quadrature are prepared once per solve.
//! Full compatible strain is bounded, including the solved thickness strain.
//! The independent small-rotation bound checks the skew displacement gradient
//! at the same integration points; a large rigid spin cannot pass as zero strain.
//!
//! For Q1 plane strain, opt-in [`Formulation::BBar`] replaces only the
//! three-dimensional volumetric strain by its quadrature-weighted element
//! mean: `eps_bar = dev(eps) + mean(trace(eps)) I / 3`. The same projected B
//! drives strain, force and tangent, including the out-of-plane normal stress.
//! This retains local deviatoric strain and full quadrature histories. The
//! returned material strain is then an ASSUMED strain: its zz component can
//! be nonzero even though the underlying displacement has zero zz strain.
//! Both raw compatible and projected strains are checked against the domain.
//! See the original INL MOOSE small-strain implementation:
//! <https://mooseframework.inl.gov/docs/doxygen/modules/Compute2DSmallStrain_8C_source.html>.
//!
//! This is small-strain, small-rotation FE with dead nodal loads. Standard
//! displacement elements can lock in confined plastic flow. B-bar is admitted
//! only for Q1 plane strain; there is no B-bar plane-stress, shear-locking-free,
//! finite-strain forming, shell bending or plastic limit-load convergence
//! claim. Use prescribed displacement past a load
//! maximum. Material coefficients and time intervals stay fixed within a
//! path; `solve_from` supports a subsequent law/interval with retained history.

use crate::linear::{Formulation, PlaneKind};
use crate::mesh2::{Mesh2, quad_points, shapes_at};
use crate::plastic3::{
    IncrementRecord, IncrementSettings, SmallStrainLawField, SmallStrainTetError,
    solve_increment_linear,
};
use fs_exec::Cx;
use fs_material::plastic::reduction::{
    PlaneStressReduction, PlaneStressWorkspace, ReducedState, ReductionConfig, ReductionError,
};
use fs_material::{SmallStrainLaw, Voigt};
use fs_sparse::{Coo, Csr};

/// Caller-owned allocation and geometry bounds, checked before preparation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlaneAssemblyBudget {
    /// Maximum mesh vertices.
    pub maximum_nodes: usize,
    /// Maximum P1/Q1 elements.
    pub maximum_elements: usize,
    /// Maximum retained integration-point histories.
    pub maximum_quadrature_points: usize,
    /// Maximum unconstrained scalar displacement unknowns.
    pub maximum_free_dofs: usize,
    /// Maximum element-matrix entries before sparse duplicate assembly.
    pub maximum_tangent_entries: usize,
    /// Maximum physical increments and retained increment records.
    pub maximum_load_steps: usize,
    /// Minimum positive `det(J) / longest_element_edge^2`, also at Q1 corners.
    pub minimum_scaled_jacobian: f64,
}

impl Default for PlaneAssemblyBudget {
    fn default() -> Self {
        Self {
            maximum_nodes: 100_000,
            maximum_elements: 100_000,
            maximum_quadrature_points: 400_000,
            maximum_free_dofs: 200_000,
            maximum_tangent_entries: 6_400_000,
            maximum_load_steps: 10_000,
            minimum_scaled_jacobian: 1e-12,
        }
    }
}

/// One homogeneous or heterogeneous, history-dependent plane continuum.
#[derive(Debug, Clone, Copy)]
pub struct SmallStrainPlaneProblem<'a, L: SmallStrainLaw> {
    /// Existing body-fitted mesh. Coordinates are in metres.
    pub mesh: &'a Mesh2,
    /// Native laws in element connectivity order; no region averaging.
    pub laws: SmallStrainLawField<'a, L>,
    /// Zero transverse strain or zero transverse stress.
    pub plane: PlaneKind,
    /// Standard displacement, or full 3-D volumetric B-bar projection for
    /// Q1 plane strain. B-bar plane stress and P1 triangles are refused.
    pub formulation: Formulation,
    /// Positive reference thickness in metres, including plane-strain slices.
    pub thickness_m: f64,
    /// Full tensor strain admission and plane-stress local-solve controls.
    /// Plane strain uses the same strain bound without a transverse solve.
    pub reduction: ReductionConfig,
    /// Maximum absolute infinitesimal in-plane rotation [rad], in `(0, 0.05]`.
    /// This geometric allowance is independent of the material strain domain.
    pub rotation_limit_rad: f64,
    /// Reference prescribed displacements `(2 * node + component, m)`.
    pub prescribed_m: &'a [(usize, f64)],
    /// Reference dead nodal forces `(2 * node + component, N)`.
    pub nodal_forces_n: &'a [(usize, f64)],
    /// Absolute factors scaling both reference loads, in chronological order.
    pub load_path: &'a [f64],
    /// Explicit geometry, memory and load-path envelope.
    pub budget: PlaneAssemblyBudget,
    /// Global Newton, line-search and shared sparse linear-solve controls.
    pub settings: IncrementSettings,
}

/// Accepted endpoint and independently retained three-dimensional histories.
#[derive(Debug, Clone, PartialEq)]
pub struct SmallStrainPlaneSolution<S> {
    /// Nodal displacement in metres.
    pub displacement_m: Vec<[f64; 2]>,
    /// Committed full material strain/history, element-major and then
    /// `quad_points` order. B-bar stores its projected assumed strain; the
    /// underlying compatible plane-strain displacement still has epsilon_zz=0.
    pub quadrature_states: Vec<ReducedState<S>>,
    /// Full endpoint stresses [Pa] from the same frozen-history update.
    pub quadrature_stress_pa: Vec<Voigt>,
    /// Point range of element `e` is `element_offsets[e]..element_offsets[e+1]`.
    pub element_offsets: Vec<usize>,
    /// Physical step records; `evolving_elements` counts elements with any
    /// changed material history, not merely a changed total strain.
    pub increments: Vec<IncrementRecord>,
}

/// Failure leaves every supplied displacement and material history untouched.
#[derive(Debug, Clone, PartialEq)]
pub enum SmallStrainPlaneError {
    /// Invalid, nonfinite, unsupported or out-of-budget input.
    InvalidInput { what: String },
    /// Refusal of the original material's checked local reduction.
    Material {
        /// Element connectivity index.
        element: usize,
        /// Local quadrature index.
        point: usize,
        /// Original material reduction failure.
        source: ReductionError,
    },
    /// No accepted equilibrium at this load-path entry.
    NewtonStalled {
        /// Index into `load_path`.
        increment: usize,
        /// Accepted relative force residuals.
        history: Vec<f64>,
    },
    /// The shared sparse solve missed its recomputed residual gate.
    LinearSolveFailed {
        /// Index into `load_path`.
        increment: usize,
        /// Original shared-solver failure.
        source: SmallStrainTetError,
    },
    /// Cancellation; no candidate histories are published.
    Cancelled,
}

impl core::fmt::Display for SmallStrainPlaneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-SOLID-PLASTIC2: {self:?}")
    }
}
impl std::error::Error for SmallStrainPlaneError {}

fn invalid(what: impl Into<String>) -> SmallStrainPlaneError {
    SmallStrainPlaneError::InvalidInput { what: what.into() }
}
fn checkpoint(cx: &Cx<'_>) -> Result<(), SmallStrainPlaneError> {
    cx.checkpoint()
        .map_err(|_| SmallStrainPlaneError::Cancelled)
}
fn norm(values: impl IntoIterator<Item = f64>) -> f64 {
    values.into_iter().fold(0.0_f64, f64::hypot)
}
fn tensor_norm(strain: &Voigt) -> f64 {
    norm(strain.iter().copied().chain(strain[3..].iter().copied()))
}
const IN_PLANE: [usize; 3] = [0, 1, 3];
const WEIGHT: [f64; 4] = [1.0, 1.0, 1.0, 2.0];

#[derive(Debug)]
struct Point {
    /// Original compatible [xx, yy, xy] columns, also used for rotation admission.
    columns: [[f64; 3]; 8],
    /// One third of (mean divergence - local divergence), per nodal DOF.
    volumetric_correction: [f64; 8],
    volume: f64,
}
impl Point {
    fn assumed_columns(&self) -> [[f64; 4]; 8] {
        std::array::from_fn(|i| {
            let b = self.columns[i];
            let correction = self.volumetric_correction[i];
            [b[0] + correction, b[1] + correction, correction, b[2]]
        })
    }
}
#[derive(Debug)]
struct Element {
    nodes: [usize; 4],
    count: usize,
    points: Vec<Point>,
}
struct Trial<S> {
    states: Vec<ReducedState<S>>,
    stress: Vec<Voigt>,
    tangent: Vec<[[f64; 4]; 4]>,
    force: Vec<f64>,
    force_scale: f64,
}

impl<L: SmallStrainLaw> SmallStrainPlaneProblem<'_, L> {
    fn law(&self, element: usize) -> &L {
        match &self.laws {
            SmallStrainLawField::Uniform(law) => law,
            SmallStrainLawField::PerElement(laws) => &laws[element],
        }
    }

    fn prepare(&self, cx: &Cx<'_>) -> Result<Vec<Element>, SmallStrainPlaneError> {
        checkpoint(cx)?;
        let n = self.mesh.nodes.len();
        let ne = self.mesh.elems.len();
        if n == 0
            || ne == 0
            || self.load_path.is_empty()
            || n > self.budget.maximum_nodes
            || ne > self.budget.maximum_elements
            || self.load_path.len() > self.budget.maximum_load_steps
            || !(self.thickness_m.is_finite() && self.thickness_m > 0.0)
            || !(self.rotation_limit_rad.is_finite()
                && self.rotation_limit_rad > 0.0
                && self.rotation_limit_rad <= 0.05)
            || !(self.budget.minimum_scaled_jacobian.is_finite()
                && self.budget.minimum_scaled_jacobian > 0.0)
        {
            return Err(invalid("empty, oversized or inadmissible plane problem"));
        }
        if let SmallStrainLawField::PerElement(laws) = &self.laws
            && laws.len() != ne
        {
            return Err(invalid("material count must equal the element count"));
        }
        if self.formulation == Formulation::BBar
            && (self.plane != PlaneKind::Strain || self.mesh.elems.iter().any(|e| e.len() != 4))
        {
            return Err(invalid("B-bar requires Q1 quadrilaterals in plane strain"));
        }
        // Validate the shared controls once, even when no transverse solve is needed.
        PlaneStressReduction::new(self.law(0), self.reduction)
            .map_err(|e| invalid(e.to_string()))?;
        let s = self.settings;
        if s.max_newton_iterations == 0
            || s.max_linear_iterations == 0
            || !(s.relative_tolerance.is_finite()
                && s.relative_tolerance > 0.0
                && s.relative_tolerance < 1.0)
            || !(s.linear_tolerance.is_finite()
                && s.linear_tolerance > 0.0
                && s.linear_tolerance < 1.0)
            || self
                .mesh
                .nodes
                .iter()
                .flatten()
                .chain(self.load_path)
                .any(|v| !v.is_finite())
        {
            return Err(invalid("invalid solve controls, coordinates or load path"));
        }
        let dofs = n
            .checked_mul(2)
            .ok_or_else(|| invalid("DOF count overflow"))?;
        let mut fixed = vec![false; dofs];
        for &(dof, value) in self.prescribed_m {
            if dof >= dofs || !value.is_finite() || fixed[dof] {
                return Err(invalid("invalid or repeated prescribed displacement"));
            }
            fixed[dof] = true;
        }
        if dofs - self.prescribed_m.len() > self.budget.maximum_free_dofs
            || self
                .nodal_forces_n
                .iter()
                .any(|&(d, v)| d >= dofs || !v.is_finite())
        {
            return Err(invalid(
                "invalid forces or oversized free displacement space",
            ));
        }
        let mut count = 0usize;
        let mut entries = 0usize;
        for conn in &self.mesh.elems {
            if !(conn.len() == 3 || conn.len() == 4) || conn.iter().any(|&v| v >= n) {
                return Err(invalid(
                    "only valid P1 triangles and Q1 quadrilaterals are supported",
                ));
            }
            count = count
                .checked_add(conn.len())
                .ok_or_else(|| invalid("point count overflow"))?;
            entries = entries
                .checked_add((2 * conn.len()).pow(2))
                .ok_or_else(|| invalid("matrix count overflow"))?;
        }
        if count > self.budget.maximum_quadrature_points
            || entries > self.budget.maximum_tangent_entries
        {
            return Err(invalid(
                "quadrature histories or tangent entries exceed the budget",
            ));
        }
        let mut elements = Vec::with_capacity(ne);
        for (e, conn) in self.mesh.elems.iter().enumerate() {
            checkpoint(cx)?;
            let mut nodes = [0; 4];
            nodes[..conn.len()].copy_from_slice(conn);
            // The original geometry owner receives translated coordinates to
            // avoid cancellation when a small element is far from the origin.
            let x: Vec<[f64; 2]> = conn
                .iter()
                .map(|&node| [0, 1].map(|c| self.mesh.nodes[node][c] - self.mesh.nodes[conn[0]][c]))
                .collect();
            let h = (0..conn.len())
                .flat_map(|a| (0..a).map(move |b| (a, b)))
                .map(|(a, b)| (x[a][0] - x[b][0]).hypot(x[a][1] - x[b][1]))
                .fold(0.0_f64, f64::max);
            if !h.is_finite() || h == 0.0 {
                return Err(invalid(format!("element {e} has invalid extent")));
            }
            let local = &[0, 1, 2, 3][..conn.len()];
            let corners = if conn.len() == 4 {
                &[(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)][..]
            } else {
                &[(0.0, 0.0)][..]
            };
            for &(xi, eta) in corners {
                let (_, _, det) = shapes_at(&x, local, xi, eta);
                let scaled = det / h / h;
                if !scaled.is_finite() || scaled < self.budget.minimum_scaled_jacobian {
                    return Err(invalid(format!(
                        "element {e} is inverted, folded or degenerate"
                    )));
                }
            }
            let mut points = Vec::with_capacity(conn.len());
            for (xi, eta, weight) in quad_points(conn.len()) {
                let (_, gradients, det) = shapes_at(&x, local, xi, eta);
                let volume = weight * det * self.thickness_m;
                if !(volume.is_finite() && volume > 0.0)
                    || gradients.iter().flatten().any(|v| !v.is_finite())
                {
                    return Err(invalid(format!(
                        "element {e} has invalid quadrature geometry"
                    )));
                }
                let mut columns = [[0.0; 3]; 8];
                for (a, g) in gradients.iter().enumerate() {
                    columns[2 * a] = [g[0], 0.0, 0.5 * g[1]];
                    columns[2 * a + 1] = [0.0, g[1], 0.5 * g[0]];
                }
                points.push(Point {
                    columns,
                    volumetric_correction: [0.0; 8],
                    volume,
                });
            }
            if self.formulation == Formulation::BBar {
                let volume: f64 = points.iter().map(|p| p.volume).sum();
                if !volume.is_finite() || volume <= 0.0 {
                    return Err(invalid(format!("element {e} has invalid B-bar volume")));
                }
                let mut mean = [0.0; 8];
                for point in &points {
                    let fraction = point.volume / volume;
                    for (d, value) in mean.iter_mut().enumerate() {
                        *value =
                            fraction.mul_add(point.columns[d][0] + point.columns[d][1], *value);
                    }
                }
                for point in &mut points {
                    for (d, value) in point.volumetric_correction.iter_mut().enumerate() {
                        *value = (mean[d] - (point.columns[d][0] + point.columns[d][1])) / 3.0;
                    }
                    if point.volumetric_correction.iter().any(|v| !v.is_finite()) {
                        return Err(invalid(format!(
                            "element {e} has nonfinite B-bar gradients"
                        )));
                    }
                }
            }
            elements.push(Element {
                nodes,
                count: conn.len(),
                points,
            });
        }
        Ok(elements)
    }

    fn check(
        &self,
        elements: &[Element],
        u: &[f64],
        states: &[ReducedState<L::State>],
    ) -> Result<(), SmallStrainPlaneError> {
        let count: usize = elements.iter().map(|e| e.points.len()).sum();
        if u.len() != 2 * self.mesh.nodes.len()
            || states.len() != count
            || u.iter().any(|v| !v.is_finite())
            || states.iter().any(|s| {
                !tensor_norm(&s.strain).is_finite()
                    || tensor_norm(&s.strain) > self.reduction.strain_limit
            })
        {
            return Err(invalid("invalid displacement or quadrature history field"));
        }
        Ok(())
    }

    fn strain_and_rotation(element: &Element, point: &Point, u: &[f64]) -> ([f64; 3], f64) {
        let mut gradient = [[0.0; 2]; 2];
        for a in 1..element.count {
            let g = [point.columns[2 * a][0], point.columns[2 * a + 1][1]];
            for (c, row) in gradient.iter_mut().enumerate() {
                let du = u[2 * element.nodes[a] + c] - u[2 * element.nodes[0] + c];
                for (j, value) in row.iter_mut().enumerate() {
                    *value = g[j].mul_add(du, *value);
                }
            }
        }
        (
            [
                gradient[0][0],
                gradient[1][1],
                0.5 * (gradient[0][1] + gradient[1][0]),
            ],
            0.5 * (gradient[1][0] - gradient[0][1]),
        )
    }

    fn evaluate(
        &self,
        elements: &[Element],
        u: &[f64],
        old: &[ReducedState<L::State>],
        fixed: &[bool],
        workspace: &mut PlaneStressWorkspace,
        cx: &Cx<'_>,
    ) -> Result<Trial<L::State>, SmallStrainPlaneError> {
        let mut trial = Trial {
            states: Vec::with_capacity(old.len()),
            stress: Vec::with_capacity(old.len()),
            tangent: Vec::with_capacity(old.len()),
            force: vec![0.0; u.len()],
            force_scale: 0.0,
        };
        let mut bounds = vec![0.0; u.len()];
        for (e, element) in elements.iter().enumerate() {
            checkpoint(cx)?;
            let law = self.law(e);
            let reduction = PlaneStressReduction::new(law, self.reduction)
                .map_err(|e| invalid(e.to_string()))?;
            for (q, point) in element.points.iter().enumerate() {
                let (strain, rotation) = Self::strain_and_rotation(element, point, u);
                if !rotation.is_finite() || rotation.abs() > self.rotation_limit_rad {
                    return Err(invalid(format!(
                        "element {e} point {q} exceeds the small-rotation domain"
                    )));
                }
                let incoming = &old[trial.states.len()];
                let (next, stress, tangent) = match self.plane {
                    PlaneKind::Stress => {
                        let response = reduction
                            .evaluate_with_workspace(strain, incoming, workspace)
                            .map_err(|source| SmallStrainPlaneError::Material {
                                element: e,
                                point: q,
                                source,
                            })?;
                        let mut tangent = [[0.0; 4]; 4];
                        for (a, &i) in IN_PLANE.iter().enumerate() {
                            for (b, &j) in IN_PLANE.iter().enumerate() {
                                tangent[i][j] = response.tangent_pa[a][b];
                            }
                        }
                        (response.state, response.full_stress_pa, tangent)
                    }
                    PlaneKind::Strain => {
                        let mut full = [strain[0], strain[1], 0.0, strain[2], 0.0, 0.0];
                        if !tensor_norm(&full).is_finite()
                            || tensor_norm(&full) > self.reduction.strain_limit
                        {
                            return Err(invalid(format!(
                                "element {e} point {q} exceeds the small-strain domain"
                            )));
                        }
                        if self.formulation == Formulation::BBar {
                            let mut correction = 0.0;
                            for a in 1..element.count {
                                for c in 0..2 {
                                    let du =
                                        u[2 * element.nodes[a] + c] - u[2 * element.nodes[0] + c];
                                    correction = point.volumetric_correction[2 * a + c]
                                        .mul_add(du, correction);
                                }
                            }
                            for value in &mut full[..3] {
                                *value += correction;
                            }
                            if !tensor_norm(&full).is_finite()
                                || tensor_norm(&full) > self.reduction.strain_limit
                            {
                                return Err(invalid(format!(
                                    "element {e} point {q} projected strain exceeds the domain"
                                )));
                            }
                        }
                        let stress = law.stress(&full, &incoming.history);
                        let tangent = law.tangent(&full, &incoming.history);
                        if tangent.iter().flatten().any(|v| !v.is_finite()) {
                            return Err(invalid(format!(
                                "element {e} point {q} has a nonfinite native tangent"
                            )));
                        }
                        (
                            ReducedState {
                                strain: full,
                                history: law.update_state(&full, &incoming.history),
                            },
                            stress,
                            std::array::from_fn(|i| std::array::from_fn(|j| tangent[i][j])),
                        )
                    }
                };
                if stress.iter().any(|v| !v.is_finite())
                    || tangent.iter().flatten().any(|v| !v.is_finite())
                {
                    return Err(invalid(format!(
                        "element {e} point {q} has a nonfinite constitutive response"
                    )));
                }
                let stress_norm = if self.formulation == Formulation::BBar {
                    norm([stress[0], stress[1], stress[2], stress[3], stress[3]])
                } else {
                    norm([stress[0], stress[1], stress[3], stress[3]])
                };
                let columns = point.assumed_columns();
                for a in 0..2 * element.count {
                    let dof = 2 * element.nodes[a / 2] + a % 2;
                    let b = columns[a];
                    let work: f64 = (0..4).map(|i| WEIGHT[i] * b[i] * stress[i]).sum();
                    trial.force[dof] = point.volume.mul_add(work, trial.force[dof]);
                    if !fixed[dof] {
                        bounds[dof] +=
                            point.volume * norm([b[0], b[1], b[2], b[3], b[3]]) * stress_norm;
                    }
                }
                trial.states.push(next);
                trial.stress.push(stress);
                trial.tangent.push(tangent);
            }
        }
        trial.force_scale = norm(bounds);
        if !trial.force_scale.is_finite() || trial.force.iter().any(|v| !v.is_finite()) {
            return Err(invalid("nonfinite assembled force or force bound"));
        }
        Ok(trial)
    }

    fn assemble(
        &self,
        elements: &[Element],
        trial: &Trial<L::State>,
        cx: &Cx<'_>,
    ) -> Result<Csr, SmallStrainPlaneError> {
        let n = 2 * self.mesh.nodes.len();
        let mut coo = Coo::new(n, n);
        let mut index = 0;
        for element in elements {
            checkpoint(cx)?;
            let mut stiffness = [[0.0; 8]; 8];
            for point in &element.points {
                let c = trial.tangent[index];
                index += 1;
                // Reuse C B for every test column and accumulate quadrature
                // locally: at most one 8-by-8 matrix enters the sparse COO.
                let columns = point.assumed_columns();
                let cb = columns.map(|b| {
                    std::array::from_fn::<_, 4, _>(|i| (0..4).map(|j| c[i][j] * b[j]).sum::<f64>())
                });
                for (p, row) in stiffness.iter_mut().enumerate().take(2 * element.count) {
                    for (q, entry) in row.iter_mut().enumerate().take(2 * element.count) {
                        let work: f64 = (0..4).map(|i| WEIGHT[i] * columns[p][i] * cb[q][i]).sum();
                        *entry = point.volume.mul_add(work, *entry);
                    }
                }
            }
            for (p, row) in stiffness.iter().enumerate().take(2 * element.count) {
                for (q, &value) in row.iter().enumerate().take(2 * element.count) {
                    coo.push(
                        2 * element.nodes[p / 2] + p % 2,
                        2 * element.nodes[q / 2] + q % 2,
                        value,
                    );
                }
            }
        }
        let matrix = coo.assemble();
        if (0..n).any(|r| matrix.row(r).1.iter().any(|v| !v.is_finite())) {
            return Err(invalid("nonfinite assembled tangent"));
        }
        Ok(matrix)
    }

    /// Follow the path from zero displacement and each law's virgin history.
    pub fn solve(
        &self,
        cx: &Cx<'_>,
    ) -> Result<SmallStrainPlaneSolution<L::State>, SmallStrainPlaneError> {
        self.solve_candidate(None, cx)
    }

    /// Continue with absolute load factors from supplied committed history.
    /// The inputs are borrowed and remain unchanged on success or refusal.
    /// Integration-point ordering must match this mesh, plane reduction and
    /// formulation; histories are not transferred between discretizations.
    pub fn solve_from(
        &self,
        displacement_m: &[f64],
        states: &[ReducedState<L::State>],
        cx: &Cx<'_>,
    ) -> Result<SmallStrainPlaneSolution<L::State>, SmallStrainPlaneError> {
        self.solve_candidate(Some((displacement_m, states)), cx)
    }

    /// Candidate internal forces from the frozen incoming histories; no commit.
    pub fn internal_force(
        &self,
        displacement_m: &[f64],
        states: &[ReducedState<L::State>],
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, SmallStrainPlaneError> {
        let elements = self.prepare(cx)?;
        self.check(&elements, displacement_m, states)?;
        let fixed = vec![false; displacement_m.len()];
        let mut workspace = PlaneStressWorkspace::new().map_err(|e| invalid(e.to_string()))?;
        Ok(self
            .evaluate(
                &elements,
                displacement_m,
                states,
                &fixed,
                &mut workspace,
                cx,
            )?
            .force)
    }

    /// Consistent global derivative of `internal_force` from the same history.
    pub fn tangent_matrix(
        &self,
        displacement_m: &[f64],
        states: &[ReducedState<L::State>],
        cx: &Cx<'_>,
    ) -> Result<Csr, SmallStrainPlaneError> {
        let elements = self.prepare(cx)?;
        self.check(&elements, displacement_m, states)?;
        let fixed = vec![false; displacement_m.len()];
        let mut workspace = PlaneStressWorkspace::new().map_err(|e| invalid(e.to_string()))?;
        let trial = self.evaluate(
            &elements,
            displacement_m,
            states,
            &fixed,
            &mut workspace,
            cx,
        )?;
        self.assemble(&elements, &trial, cx)
    }

    fn solve_candidate(
        &self,
        initial: Option<(&[f64], &[ReducedState<L::State>])>,
        cx: &Cx<'_>,
    ) -> Result<SmallStrainPlaneSolution<L::State>, SmallStrainPlaneError> {
        let elements = self.prepare(cx)?;
        let symmetric = (0..elements.len()).all(|e| self.law(e).admissibility().tangent_symmetric);
        let n = 2 * self.mesh.nodes.len();
        let mut fixed = vec![false; n];
        for &(d, _) in self.prescribed_m {
            fixed[d] = true;
        }
        let free: Vec<usize> = (0..n).filter(|&d| !fixed[d]).collect();
        let mut row_of = vec![usize::MAX; n];
        for (row, &d) in free.iter().enumerate() {
            row_of[d] = row;
        }
        let mut offsets = vec![0];
        for element in &elements {
            offsets.push(offsets.last().unwrap() + element.points.len());
        }
        let (mut u, mut states) = if let Some((u, states)) = initial {
            self.check(&elements, u, states)?;
            (u.to_vec(), states.to_vec())
        } else {
            (
                vec![0.0; n],
                elements
                    .iter()
                    .enumerate()
                    .flat_map(|(e, element)| {
                        (0..element.points.len()).map(move |_| ReducedState {
                            strain: [0.0; 6],
                            history: self.law(e).initial_state(),
                        })
                    })
                    .collect(),
            )
        };
        let mut workspace = PlaneStressWorkspace::new().map_err(|e| invalid(e.to_string()))?;
        let mut increments = Vec::with_capacity(self.load_path.len());
        let mut stresses = Vec::new();
        for (increment, &load) in self.load_path.iter().enumerate() {
            checkpoint(cx)?;
            let mut external = vec![0.0; n];
            for &(d, f) in self.nodal_forces_n {
                external[d] += load * f;
            }
            let targets: Vec<(usize, f64)> = self
                .prescribed_m
                .iter()
                .map(|&(d, v)| (d, load * v))
                .collect();
            if external
                .iter()
                .chain(targets.iter().map(|(_, v)| v))
                .any(|v| !v.is_finite())
            {
                return Err(invalid(
                    "load factor produces nonfinite forces or displacements",
                ));
            }
            let mut trial = self.evaluate(&elements, &u, &states, &fixed, &mut workspace, cx)?;
            let mut scale = norm(free.iter().map(|&d| external[d]))
                .max(trial.force_scale)
                .max(f64::MIN_POSITIVE);
            let mut history = Vec::new();
            let mut fallbacks = 0;
            let mut converged = false;
            for iteration in 0..=self.settings.max_newton_iterations {
                checkpoint(cx)?;
                let residual: Vec<f64> =
                    free.iter().map(|&d| trial.force[d] - external[d]).collect();
                let before = norm(residual.iter().copied());
                let pending = targets.iter().any(|&(d, v)| u[d] != v);
                history.push(before / scale);
                if !pending && before <= self.settings.relative_tolerance * scale {
                    converged = true;
                    break;
                }
                if iteration == self.settings.max_newton_iterations {
                    break;
                }
                let mut dp = vec![0.0; n];
                for &(d, v) in &targets {
                    dp[d] = v - u[d];
                }
                let mut rhs: Vec<f64> = residual.iter().map(|r| -r).collect();
                let mut coo = Coo::new(free.len(), free.len());
                // A fully prescribed body needs only a checked material update.
                if !free.is_empty() {
                    let full = self.assemble(&elements, &trial, cx)?;
                    for (row, &d) in free.iter().enumerate() {
                        let (cols, values) = full.row(d);
                        for (&c, &v) in cols.iter().zip(values) {
                            if fixed[c] {
                                rhs[row] = (-v).mul_add(dp[c], rhs[row]);
                            } else {
                                coo.push(row, row_of[c], v);
                            }
                        }
                    }
                }
                let (direction, fallback) = solve_increment_linear(
                    &coo.assemble(),
                    &rhs,
                    increment,
                    self.settings,
                    symmetric,
                    cx,
                )
                .map_err(|source| match source {
                    SmallStrainTetError::Cancelled => SmallStrainPlaneError::Cancelled,
                    source => SmallStrainPlaneError::LinearSolveFailed { increment, source },
                })?;
                fallbacks += usize::from(fallback);
                let mut alpha = 1.0_f64;
                let mut accepted = None;
                for _ in 0..=self.settings.max_backtracks {
                    checkpoint(cx)?;
                    let mut candidate = u.clone();
                    for (row, &d) in free.iter().enumerate() {
                        candidate[d] = alpha.mul_add(direction[row], u[d]);
                    }
                    for &(d, v) in &targets {
                        candidate[d] = if alpha == 1.0 {
                            v
                        } else {
                            alpha.mul_add(dp[d], u[d])
                        };
                    }
                    match self.evaluate(&elements, &candidate, &states, &fixed, &mut workspace, cx)
                    {
                        Ok(next) => {
                            let after = norm(free.iter().map(|&d| next.force[d] - external[d]));
                            if pending
                                || after < before
                                || after
                                    <= self.settings.relative_tolerance
                                        * scale.max(next.force_scale)
                            {
                                accepted = Some((candidate, next));
                                break;
                            }
                        }
                        Err(SmallStrainPlaneError::Cancelled) => {
                            return Err(SmallStrainPlaneError::Cancelled);
                        }
                        Err(_) => {} // A refused material-domain proposal may be shortened.
                    }
                    alpha *= 0.5;
                }
                let Some((candidate, next)) = accepted else {
                    return Err(SmallStrainPlaneError::NewtonStalled { increment, history });
                };
                u = candidate;
                scale = scale.max(next.force_scale); // Rejected proposals cannot weaken acceptance.
                trial = next;
            }
            if !converged {
                return Err(SmallStrainPlaneError::NewtonStalled { increment, history });
            }
            let evolving_elements = offsets
                .windows(2)
                .filter(|range| {
                    (range[0]..range[1]).any(|q| trial.states[q].history != states[q].history)
                })
                .count();
            increments.push(IncrementRecord {
                load_factor: load,
                residual_history: history,
                reactions_n: self
                    .prescribed_m
                    .iter()
                    .map(|&(d, _)| trial.force[d] - external[d])
                    .collect(),
                evolving_elements,
                indefinite_fallbacks: fallbacks,
            });
            // Keep stress and history from ONE accepted frozen-history update.
            stresses = trial.stress;
            states = trial.states;
        }
        checkpoint(cx)?;
        Ok(SmallStrainPlaneSolution {
            displacement_m: u.as_chunks::<2>().0.to_vec(),
            quadrature_states: states,
            quadrature_stress_pa: stresses,
            element_offsets: offsets,
            increments,
        })
    }
}
