//! Finite-strain Taylor–Hood solids: quadratic displacement, continuous P1 pressure.
//!
//! The material owner supplies the original mixed potential and all derivatives.
//! Reference Tet10 geometry is prepared once. Element matrices are accumulated
//! locally before sparse assembly; the native sparse LDLᵀ analysis is reused
//! across Newton iterates. Native Newton–Krylov supplies residual globalization.
//! A positive work-conjugate scaling preserves the saddle-point inertia and
//! prevents forces [N] from being compared directly with volumes [m³].
//!
//! This is a bounded spatial discretization, not a general inf-sup or mesh
//! convergence certificate. The retained motion owner adds the native consistent
//! reference inertia to displacement rows; pressure has no fictitious mass.
//! Exact incompressibility is imposed weakly in the
//! continuous P1 pressure space. Current Jacobians are checked at quadrature
//! points; that does not certify injectivity between samples. The sparse direct
//! factorization is one bounded, non-interruptible kernel between checkpoints.

use std::cell::RefCell;
use std::sync::Arc;

use fs_exec::Cx;
use fs_material::MaterialError;
use fs_material::hyper::Hyperelastic;
use fs_material::hyper::mixed::{MixedFiniteStrainPotential, VolumetricTreatment};
use fs_solver::nonlinear::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, NonlinearProblem,
    SolverRunProgress,
};
use fs_sparse::{
    Coo, Csr, DirectOrdering, Inertia, LdltError, LdltFactor, LdltOptions, SymbolicLdlt,
};

use crate::linear3::{TetAssemblyBudget, TetElasticError};

mod element;
use element::PreparedMixedMesh;

/// Retained finite-strain viscoelastic histories on the mixed solid.
pub mod visco;

/// Objective discrete-gradient elasticity with conservative reference inertia.
pub mod conservative;

/// Constitutive laws retain their spatial identity; no coefficients are averaged.
#[derive(Debug)]
pub enum MixedLawField<'a, L: MixedFiniteStrainPotential = Hyperelastic> {
    Uniform(&'a L),
    PerElement(&'a [L]),
    /// Element-major, then the original deterministic Duffy quadrature order.
    PerQuadrature(&'a [L]),
}

impl<L: MixedFiniteStrainPotential> Copy for MixedLawField<'_, L> {}
impl<L: MixedFiniteStrainPotential> Clone for MixedLawField<'_, L> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<L: MixedFiniteStrainPotential> MixedLawField<'_, L> {
    fn validate(&self, elements: usize, points: usize, cx: &Cx<'_>) -> Result<(), MixedUpError> {
        let (laws, expected): (&[L], usize) = match self {
            Self::Uniform(law) => (core::slice::from_ref(*law), 1),
            Self::PerElement(laws) => (laws, elements),
            Self::PerQuadrature(laws) => (
                laws,
                elements
                    .checked_mul(points)
                    .ok_or_else(|| invalid("material count overflow"))?,
            ),
        };
        if laws.len() != expected {
            return Err(invalid(
                "mixed law field does not cover every element/quadrature point",
            ));
        }
        for (index, law) in laws.iter().enumerate() {
            if index % 128 == 0 {
                checkpoint(cx)?;
            }
            law.validate_mixed()
                .map_err(|error| MixedUpError::Material {
                    element: index,
                    point: 0,
                    error,
                })?;
        }
        Ok(())
    }

    fn at(&self, element: usize, point: usize, points: usize) -> &L {
        match self {
            Self::Uniform(law) => law,
            Self::PerElement(laws) => &laws[element],
            Self::PerQuadrature(laws) => &laws[element * points + point],
        }
    }
}

/// Explicit numerical scales and bounded native solver controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixedUpSettings {
    pub load_steps: usize,
    pub maximum_newton_iterations: usize,
    pub newton: NewtonKrylovConfig,
    /// Characteristic displacement/geometry scale, in m.
    pub length_scale_m: f64,
    /// Characteristic shear/load stress, in Pa; do not scale by a huge bulk penalty.
    pub stress_scale_pa: f64,
    /// Conservative sparse assembly payload cap, separate from reference points.
    pub maximum_assembly_bytes: usize,
    /// Native numeric factor and scratch payload cap, not process RSS.
    pub maximum_factor_bytes: usize,
}

/// Static finite-strain problem with an independent continuous pressure field.
#[derive(Debug)]
pub struct MixedUpProblem<'a, L: MixedFiniteStrainPotential = Hyperelastic> {
    pub nodes_m: &'a [[f64; 3]],
    pub elements: &'a [[usize; 10]],
    pub material: MixedLawField<'a, L>,
    pub volumetric: VolumetricTreatment,
    /// Node-major xyz DOF and absolute displacement at full load [m].
    pub prescribed_m: &'a [(usize, f64)],
    pub nodal_forces_n: &'a [(usize, f64)],
    /// Dead force per original reference volume [N/m³].
    pub body_force_n_m3: [f64; 3],
    /// Optional constant-pressure gauge: original corner-node ID and value [Pa].
    /// Only exact incompressibility admits this, and the selected component's
    /// constant-pressure mode must exert zero force on every free displacement.
    /// Traction boundaries normally determine pressure and require `None`.
    pub pressure_gauge: Option<(usize, f64)>,
    pub budget: TetAssemblyBudget,
    pub quadrature_points_per_dir: usize,
    pub maximum_quadrature_points: usize,
    pub settings: MixedUpSettings,
}

impl<L: MixedFiniteStrainPotential> Copy for MixedUpProblem<'_, L> {}
impl<L: MixedFiniteStrainPotential> Clone for MixedUpProblem<'_, L> {
    fn clone(&self) -> Self {
        *self
    }
}

/// Raw physical residual and consistent tangent, including constrained rows.
#[derive(Clone, Debug)]
pub struct MixedUpEvaluation {
    pub internal_forces_n: Vec<f64>,
    pub pressure_residual_m3: Vec<f64>,
    /// All displacement DOFs followed by pressures in `pressure_nodes()` order.
    pub tangent: Option<Csr>,
    /// Mixed saddle potential [J], not necessarily physical stored energy.
    pub potential_j: f64,
    pub min_det_f: f64,
}

#[derive(Clone, Debug)]
pub struct MixedUpSolution {
    pub displacement_m: Vec<[f64; 3]>,
    pub pressure_nodes: Vec<usize>,
    /// Tensile bulk multiplier; see the material's original split convention.
    pub pressure_pa: Vec<f64>,
    pub reactions_n: Vec<(usize, f64)>,
    pub potential_j: f64,
    pub min_det_f: f64,
    /// Final constrained, congruently scaled KKT inertia.
    pub kkt_inertia: Inertia,
    pub steps: Vec<NewtonReport>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum MixedUpError {
    InvalidInput {
        what: String,
    },
    Geometry(TetElasticError),
    Material {
        element: usize,
        point: usize,
        error: MaterialError,
    },
    Direct(LdltError),
    KktInertia {
        expected: Inertia,
        got: Inertia,
    },
    NewtonConfiguration(NewtonError),
    NewtonStalled {
        load_factor: f64,
        report: NewtonReport,
    },
    IncompatiblePressureGauge {
        node: usize,
    },
    Cancelled,
}

impl core::fmt::Display for MixedUpError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput { what } => write!(f, "mixed solid input: {what}"),
            Self::Geometry(error) => write!(f, "mixed reference geometry: {error}"),
            Self::Material {
                element,
                point,
                error,
            } => write!(
                f,
                "mixed material element {element}, point {point}: {error}"
            ),
            Self::Direct(error) => write!(f, "mixed KKT factorization: {error}"),
            Self::KktInertia { expected, got } => {
                write!(f, "mixed KKT inertia {got:?}, expected {expected:?}")
            }
            Self::NewtonConfiguration(error) => write!(f, "mixed Newton configuration: {error}"),
            Self::NewtonStalled {
                load_factor,
                report,
            } => write!(
                f,
                "mixed Newton stalled at load {load_factor}: {:?}, residual {}",
                report.diagnosis, report.residual_norm
            ),
            Self::IncompatiblePressureGauge { node } => write!(
                f,
                "pressure pin at node {node} is not a free constant-pressure gauge, or prescribed volume is incompatible"
            ),
            Self::Cancelled => write!(f, "mixed solid cancelled"),
        }
    }
}
impl core::error::Error for MixedUpError {}

fn invalid(what: impl Into<String>) -> MixedUpError {
    MixedUpError::InvalidInput { what: what.into() }
}
fn checkpoint(cx: &Cx<'_>) -> Result<(), MixedUpError> {
    cx.checkpoint().map_err(|_| MixedUpError::Cancelled)
}

fn inertia_error(error: crate::inertia::InertiaError) -> MixedUpError {
    match error {
        crate::inertia::InertiaError::Cancelled => MixedUpError::Cancelled,
        crate::inertia::InertiaError::InvalidInput { what } => {
            invalid(format!("mixed reference inertia: {what}"))
        }
    }
}

/// Immutable shared geometry; constitutive histories remain borrowed and fixed.
pub struct PreparedMixedUpProblem<'a, L: MixedFiniteStrainPotential = Hyperelastic> {
    problem: MixedUpProblem<'a, L>,
    mesh: Arc<PreparedMixedMesh>,
    external_n: Vec<f64>,
    free: Vec<usize>,
    free_index: Vec<Option<usize>>,
    free_displacements: usize,
    gauge_pressure: Option<usize>,
    gauge_component: Vec<bool>,
    work_scale_j: f64,
}

impl<'a, L: MixedFiniteStrainPotential> MixedUpProblem<'a, L> {
    fn validate_controls(&self, cx: &Cx<'_>) -> Result<(), MixedUpError> {
        checkpoint(cx)?;
        let s = self.settings;
        if s.load_steps == 0
            || s.maximum_newton_iterations == 0
            || !(s.length_scale_m.is_finite() && s.length_scale_m > 0.0)
            || !(s.stress_scale_pa.is_finite() && s.stress_scale_pa > 0.0)
            || s.maximum_factor_bytes == 0
            || self.body_force_n_m3.iter().any(|x| !x.is_finite())
        {
            return Err(invalid(
                "invalid mixed solver scales, iteration budgets or body force",
            ));
        }
        s.newton
            .validate()
            .map_err(MixedUpError::NewtonConfiguration)?;
        self.validate_assembly_budget(false)
    }

    fn validate_assembly_budget(&self, with_inertia: bool) -> Result<(), MixedUpError> {
        // COO, stable-sort scratch, CSR and vectors. Count element blocks once,
        // not once per quadrature point. Inertia adds three scalar 10x10 blocks.
        // Native factor storage has its own cap.
        let entries = 34 * 34 + if with_inertia { 3 * 10 * 10 } else { 0 };
        let estimated = self
            .elements
            .len()
            .checked_mul(entries * 128)
            .and_then(|v| {
                self.nodes_m
                    .len()
                    .checked_mul(4 * 128)
                    .and_then(|n| v.checked_add(n))
            })
            .ok_or_else(|| invalid("mixed assembly size overflow"))?;
        if estimated > self.settings.maximum_assembly_bytes {
            return Err(invalid(format!(
                "mixed assembly needs at most {estimated} payload bytes; cap {}",
                self.settings.maximum_assembly_bytes
            )));
        }
        Ok(())
    }

    pub fn prepare(&self, cx: &Cx<'_>) -> Result<PreparedMixedUpProblem<'a, L>, MixedUpError> {
        self.validate_controls(cx)?;
        let mesh = Arc::new(
            PreparedMixedMesh::new(
                self.nodes_m,
                self.elements,
                self.budget,
                self.quadrature_points_per_dir,
                self.maximum_quadrature_points,
                cx,
            )
            .map_err(MixedUpError::Geometry)?,
        );
        self.prepare_with_mesh(mesh, cx)
    }

    fn prepare_with_mesh(
        &self,
        mesh: Arc<PreparedMixedMesh>,
        cx: &Cx<'_>,
    ) -> Result<PreparedMixedUpProblem<'a, L>, MixedUpError> {
        self.validate_controls(cx)?;
        if mesh.nodes_count != self.nodes_m.len()
            || mesh.elements.len() != self.elements.len()
            || mesh.nodes_count > self.budget.maximum_nodes
            || mesh.elements.len() > self.budget.maximum_tetrahedra
            || mesh
                .elements
                .len()
                .checked_mul(mesh.points_per_element)
                .is_none_or(|n| n > self.maximum_quadrature_points)
            || mesh.points_per_element != self.quadrature_points_per_dir.checked_pow(3).unwrap_or(0)
        {
            return Err(invalid(
                "cached mixed geometry differs from the admitted mesh/budget",
            ));
        }
        self.material
            .validate(mesh.elements.len(), mesh.points_per_element, cx)?;
        let nu = mesh.nodes_count * 3;
        let np = mesh.pressure_nodes.len();
        let mut constrained = vec![false; nu + np];
        for &(d, value) in self.prescribed_m {
            if d >= nu || !value.is_finite() || constrained[d] {
                return Err(invalid("invalid or duplicate prescribed displacement"));
            }
            constrained[d] = true;
        }
        let gauge_pressure = match self.pressure_gauge {
            Some((node, value)) => {
                if self.volumetric != VolumetricTreatment::ExactIncompressible || !value.is_finite()
                {
                    return Err(invalid(
                        "pressure gauges require exact incompressibility and a finite value",
                    ));
                }
                let p = mesh
                    .pressure_nodes
                    .binary_search(&node)
                    .map_err(|_| invalid("pressure gauge must name a corner node"))?;
                constrained[nu + p] = true;
                Some(p)
            }
            None => None,
        };
        let free_displacements = constrained[..nu].iter().filter(|&&v| !v).count();
        if free_displacements > self.budget.maximum_free_dofs {
            return Err(invalid("free displacement DOF budget exceeded"));
        }
        let free: Vec<_> = (0..nu + np).filter(|&d| !constrained[d]).collect();
        if free.is_empty() {
            return Err(invalid("mixed problem has no free unknowns"));
        }
        let mut free_index = vec![None; nu + np];
        for (i, &d) in free.iter().enumerate() {
            free_index[d] = Some(i);
        }
        let mut external_n = vec![0.0; nu];
        for &(d, value) in self.nodal_forces_n {
            if d >= nu || !value.is_finite() {
                return Err(invalid("invalid nodal force"));
            }
            external_n[d] += value;
        }
        for element in &mesh.elements {
            checkpoint(cx)?;
            for q in &element.points {
                for a in 0..10 {
                    for i in 0..3 {
                        external_n[3 * element.connectivity[a] + i] +=
                            q.weight_m3 * q.shapes[a] * self.body_force_n_m3[i];
                    }
                }
            }
        }
        let work_scale_j = self.settings.stress_scale_pa * mesh.total_volume_m3;
        if !work_scale_j.is_finite()
            || work_scale_j <= 0.0
            || external_n.iter().any(|x| !x.is_finite())
        {
            return Err(invalid("mixed force/work scale is unrepresentable"));
        }
        let mut gauge_component = vec![false; np];
        if let Some(pin) = gauge_pressure {
            let mut roots: Vec<_> = (0..np).collect();
            fn root(roots: &[usize], mut p: usize) -> usize {
                while roots[p] != p {
                    p = roots[p];
                }
                p
            }
            for element in &mesh.elements {
                for &p in &element.pressure_indices[1..] {
                    let a = root(&roots, element.pressure_indices[0]);
                    let b = root(&roots, p);
                    roots[a.max(b)] = a.min(b);
                }
            }
            let selected = root(&roots, pin);
            for (p, selected_component) in gauge_component.iter_mut().enumerate() {
                *selected_component = root(&roots, p) == selected;
            }
        }
        checkpoint(cx)?;
        Ok(PreparedMixedUpProblem {
            problem: *self,
            mesh,
            external_n,
            free,
            free_index,
            free_displacements,
            gauge_pressure,
            gauge_component,
            work_scale_j,
        })
    }

    pub fn solve(&self, cx: &Cx<'_>) -> Result<MixedUpSolution, MixedUpError> {
        self.prepare(cx)?.solve(cx)
    }

    pub fn solve_from(
        &self,
        displacement_m: &[[f64; 3]],
        pressure_pa: &[f64],
        cx: &Cx<'_>,
    ) -> Result<MixedUpSolution, MixedUpError> {
        self.prepare(cx)?
            .solve_from(displacement_m, pressure_pa, cx)
    }
}

impl<L: MixedFiniteStrainPotential> PreparedMixedUpProblem<'_, L> {
    pub fn pressure_nodes(&self) -> &[usize] {
        &self.mesh.pressure_nodes
    }
    pub fn reference_volume_m3(&self) -> f64 {
        self.mesh.total_volume_m3
    }

    fn variable_scale(&self, dof: usize) -> f64 {
        if dof < 3 * self.mesh.nodes_count {
            self.problem.settings.length_scale_m
        } else {
            self.problem.settings.stress_scale_pa
        }
    }

    fn validate_state(&self, u: &[[f64; 3]], p: &[f64]) -> Result<(), MixedUpError> {
        if u.len() != self.mesh.nodes_count
            || p.len() != self.mesh.pressure_nodes.len()
            || u.iter().flatten().chain(p.iter()).any(|v| !v.is_finite())
        {
            return Err(invalid(
                "mixed state dimensions or finite entries are invalid",
            ));
        }
        Ok(())
    }

    pub fn evaluate(
        &self,
        displacement_m: &[[f64; 3]],
        pressure_pa: &[f64],
        with_tangent: bool,
        cx: &Cx<'_>,
    ) -> Result<MixedUpEvaluation, MixedUpError> {
        self.evaluate_with_inertia(displacement_m, pressure_pa, with_tangent, None, cx)
    }

    // The motion owner solves the same mixed system, adding M*a and M/dt^2
    // only to displacement rows. The reported potential remains material-only.
    fn evaluate_with_inertia(
        &self,
        displacement_m: &[[f64; 3]],
        pressure_pa: &[f64],
        with_tangent: bool,
        inertia: Option<&crate::inertia::BackwardEulerInertia<'_>>,
        cx: &Cx<'_>,
    ) -> Result<MixedUpEvaluation, MixedUpError> {
        checkpoint(cx)?;
        self.validate_state(displacement_m, pressure_pa)?;
        let nu = self.mesh.nodes_count * 3;
        let np = self.mesh.pressure_nodes.len();
        let mut residual = vec![0.0; nu + np];
        let mut matrix = with_tangent.then(|| Coo::new(nu + np, nu + np));
        let mut potential_j = 0.0;
        let mut min_det_f = f64::INFINITY;
        for (e, element) in self.mesh.elements.iter().enumerate() {
            checkpoint(cx)?;
            let dofs: [usize; 34] = core::array::from_fn(|d| {
                if d < 30 {
                    3 * element.connectivity[d / 3] + d % 3
                } else {
                    nu + element.pressure_indices[d - 30]
                }
            });
            let mut local_r = [0.0; 34];
            let mut local_k = [[0.0; 34]; 34];
            for (q, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let mut f = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
                for a in 1..10 {
                    for i in 0..3 {
                        for j in 0..3 {
                            f[3 * i + j] += (displacement_m[element.connectivity[a]][i]
                                - displacement_m[element.connectivity[0]][i])
                                * point.gradients[a][j];
                        }
                    }
                }
                let p = (0..4)
                    .map(|a| pressure_pa[element.pressure_indices[a]] * point.pressure_shapes[a])
                    .sum();
                let law = self.problem.material.at(e, q, self.mesh.points_per_element);
                let response = if with_tangent {
                    Some(
                        law.mixed_response(&f, p, self.problem.volumetric)
                            .map_err(|error| MixedUpError::Material {
                                element: e,
                                point: q,
                                error,
                            })?,
                    )
                } else {
                    None
                };
                let state = match &response {
                    Some(r) => r.state.clone(),
                    None => law
                        .mixed_state(&f, p, self.problem.volumetric)
                        .map_err(|error| MixedUpError::Material {
                            element: e,
                            point: q,
                            error,
                        })?,
                };
                potential_j += point.weight_m3 * state.potential_pa;
                min_det_f = min_det_f.min(state.volume_ratio);
                for a in 0..10 {
                    for i in 0..3 {
                        local_r[3 * a + i] += point.weight_m3
                            * (0..3)
                                .map(|j| state.first_piola_pa[3 * i + j] * point.gradients[a][j])
                                .sum::<f64>();
                    }
                }
                for a in 0..4 {
                    local_r[30 + a] +=
                        point.weight_m3 * point.pressure_shapes[a] * state.pressure_residual;
                }
                if let Some(r) = response {
                    // Contract A with one shape gradient once, then reuse it
                    // for all test functions (3510 rather than 8100 products).
                    let a_gradient: [[f64; 30]; 9] = core::array::from_fn(|ij| {
                        core::array::from_fn(|bk| {
                            (0..3)
                                .map(|l| {
                                    r.tangent_pa[ij][3 * (bk % 3) + l] * point.gradients[bk / 3][l]
                                })
                                .sum()
                        })
                    });
                    for a in 0..10 {
                        for i in 0..3 {
                            let row = 3 * a + i;
                            for b in 0..10 {
                                for k in 0..3 {
                                    let value: f64 = (0..3)
                                        .map(|j| {
                                            point.gradients[a][j] * a_gradient[3 * i + j][3 * b + k]
                                        })
                                        .sum();
                                    local_k[row][3 * b + k] += point.weight_m3 * value;
                                }
                            }
                            let coupling = point.weight_m3
                                * (0..3)
                                    .map(|j| r.pressure_coupling[3 * i + j] * point.gradients[a][j])
                                    .sum::<f64>();
                            for b in 0..4 {
                                let value = coupling * point.pressure_shapes[b];
                                local_k[row][30 + b] += value;
                                local_k[30 + b][row] += value;
                            }
                        }
                    }
                    for a in 0..4 {
                        for b in 0..4 {
                            local_k[30 + a][30 + b] += point.weight_m3
                                * point.pressure_shapes[a]
                                * point.pressure_shapes[b]
                                * r.pressure_tangent_per_pa;
                        }
                    }
                }
            }
            for a in 0..34 {
                residual[dofs[a]] += local_r[a];
                if let Some(coo) = matrix.as_mut() {
                    for b in 0..34 {
                        coo.push(dofs[a], dofs[b], local_k[a][b]);
                    }
                }
            }
        }
        if let Some(inertia) = inertia {
            inertia
                .add_force(displacement_m.as_flattened(), &mut residual[..nu], cx)
                .map_err(inertia_error)?;
            if let Some(coo) = matrix.as_mut() {
                inertia.add_tangent(coo, cx).map_err(inertia_error)?;
            }
        }
        if !potential_j.is_finite() || residual.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "assembled mixed residual or potential is unrepresentable",
            ));
        }
        checkpoint(cx)?;
        let tangent = matrix.map(|m| m.assemble());
        if let Some(k) = &tangent {
            for row in 0..k.nrows() {
                if row % 128 == 0 {
                    checkpoint(cx)?;
                }
                if k.row(row).1.iter().any(|v| !v.is_finite()) {
                    return Err(invalid("assembled mixed tangent is unrepresentable"));
                }
            }
        }
        Ok(MixedUpEvaluation {
            internal_forces_n: residual[..nu].to_vec(),
            pressure_residual_m3: residual[nu..].to_vec(),
            tangent,
            potential_j,
            min_det_f,
        })
    }

    pub fn solve(&self, cx: &Cx<'_>) -> Result<MixedUpSolution, MixedUpError> {
        self.solve_from(
            &vec![[0.0; 3]; self.mesh.nodes_count],
            &vec![0.0; self.mesh.pressure_nodes.len()],
            cx,
        )
    }

    /// Continue from supplied displacements and pressures; each endpoint applies
    /// the problem's absolute support targets. No constitutive history changes.
    pub fn solve_from(
        &self,
        displacement_m: &[[f64; 3]],
        pressure_pa: &[f64],
        cx: &Cx<'_>,
    ) -> Result<MixedUpSolution, MixedUpError> {
        self.solve_loaded_from(displacement_m, pressure_pa, None, cx)
    }

    fn solve_loaded_from(
        &self,
        displacement_m: &[[f64; 3]],
        pressure_pa: &[f64],
        inertia: Option<&crate::inertia::BackwardEulerInertia<'_>>,
        cx: &Cx<'_>,
    ) -> Result<MixedUpSolution, MixedUpError> {
        checkpoint(cx)?;
        self.validate_state(displacement_m, pressure_pa)?;
        if inertia.is_some() {
            if self.problem.settings.load_steps != 1 {
                return Err(invalid(
                    "one physical endpoint per mixed momentum solve required",
                ));
            }
            self.problem.validate_assembly_budget(true)?;
        }
        let nu = 3 * self.mesh.nodes_count;
        let mut u = displacement_m.to_vec();
        let mut p = pressure_pa.to_vec();
        let symbolic = RefCell::new(None);
        let mut steps = Vec::new();
        let mut last_inertia = None;
        for step in 1..=self.problem.settings.load_steps {
            checkpoint(cx)?;
            let load = step as f64 / self.problem.settings.load_steps as f64;
            for &(d, value) in self.problem.prescribed_m {
                u[d / 3][d % 3] = load * value;
            }
            if let (Some(pin), Some((_, value))) =
                (self.gauge_pressure, self.problem.pressure_gauge)
            {
                p[pin] = load * value;
            }
            let initial: Vec<_> = self
                .free
                .iter()
                .map(|&d| {
                    (if d < nu { u[d / 3][d % 3] } else { p[d - nu] }) / self.variable_scale(d)
                })
                .collect();
            let endpoint = Endpoint {
                prepared: self,
                base_u: &u,
                base_p: &p,
                load,
                inertia,
                cx,
                symbolic: &symbolic,
                linearization: RefCell::new(None),
                fatal_error: RefCell::new(None),
                trial_error: RefCell::new(None),
            };
            let mut newton =
                match NewtonKrylovState::new(&endpoint, initial, self.problem.settings.newton) {
                    Ok(state) => state,
                    Err(error) => {
                        return Err(endpoint
                            .trial_error
                            .borrow_mut()
                            .take()
                            .unwrap_or(MixedUpError::NewtonConfiguration(error)));
                    }
                };
            let target = self
                .problem
                .settings
                .newton
                .absolute_tolerance
                .max(self.problem.settings.newton.relative_tolerance * newton.residual_norm());
            let result = newton.run_controlled(
                &endpoint,
                self.problem.settings.maximum_newton_iterations,
                &mut || cx.checkpoint().is_err(),
            );
            checkpoint(cx)?;
            if let Some(error) = endpoint.fatal_error.borrow_mut().take() {
                return Err(error);
            }
            if result.progress != SolverRunProgress::Complete {
                return Err(invalid(format!(
                    "native mixed Newton invocation did not complete: {:?}",
                    result.progress
                )));
            }
            if !result.report.converged {
                return Err(MixedUpError::NewtonStalled {
                    load_factor: load,
                    report: result.report,
                });
            }
            // Zero-residual initial states also need the stability/gauge gate.
            endpoint.ensure_linearization(&newton.x)?;
            let (next_u, next_p) = endpoint.physical_state(&newton.x)?;
            let state = self.evaluate_with_inertia(&next_u, &next_p, false, inertia, cx)?;
            let force_residual: Vec<_> = self
                .free
                .iter()
                .copied()
                .take(self.free_displacements)
                .map(|d| {
                    (state.internal_forces_n[d] - load * self.external_n[d])
                        * self.problem.settings.length_scale_m
                        / self.work_scale_j
                })
                .collect();
            let volume_residual: Vec<_> = state
                .pressure_residual_m3
                .iter()
                .map(|r| r * self.problem.settings.stress_scale_pa / self.work_scale_j)
                .collect();
            if fs_solver::norm2(&force_residual) > target
                || fs_solver::norm2(&volume_residual) > target
            {
                if let Some((node, _)) = self.problem.pressure_gauge {
                    return Err(MixedUpError::IncompatiblePressureGauge { node });
                }
                return Err(MixedUpError::NewtonStalled {
                    load_factor: load,
                    report: result.report,
                });
            }
            last_inertia = Some(
                endpoint
                    .linearization
                    .borrow()
                    .as_ref()
                    .expect("prepared final tangent")
                    .factor
                    .inertia(),
            );
            steps.push(result.report);
            u = next_u;
            p = next_p;
        }
        let state = self.evaluate_with_inertia(&u, &p, false, inertia, cx)?;
        let reactions_n = self
            .problem
            .prescribed_m
            .iter()
            .map(|&(d, _)| (d, state.internal_forces_n[d] - self.external_n[d]))
            .collect();
        checkpoint(cx)?;
        Ok(MixedUpSolution {
            displacement_m: u,
            pressure_nodes: self.mesh.pressure_nodes.clone(),
            pressure_pa: p,
            reactions_n,
            potential_j: state.potential_j,
            min_det_f: state.min_det_f,
            kkt_inertia: last_inertia.expect("positive load-step count"),
            steps,
        })
    }
}

struct Linearization {
    x: Vec<f64>,
    matrix: Csr,
    factor: LdltFactor,
}

struct Endpoint<'p, 'a, 'c, 'scope, L: MixedFiniteStrainPotential> {
    prepared: &'p PreparedMixedUpProblem<'a, L>,
    base_u: &'p [[f64; 3]],
    base_p: &'p [f64],
    load: f64,
    inertia: Option<&'p crate::inertia::BackwardEulerInertia<'p>>,
    cx: &'c Cx<'scope>,
    symbolic: &'p RefCell<Option<SymbolicLdlt>>,
    linearization: RefCell<Option<Linearization>>,
    fatal_error: RefCell<Option<MixedUpError>>,
    trial_error: RefCell<Option<MixedUpError>>,
}

impl<L: MixedFiniteStrainPotential> Endpoint<'_, '_, '_, '_, L> {
    fn physical_state(&self, x: &[f64]) -> Result<(Vec<[f64; 3]>, Vec<f64>), MixedUpError> {
        if x.len() != self.prepared.free.len() || x.iter().any(|v| !v.is_finite()) {
            return Err(invalid("native mixed Newton state is invalid"));
        }
        let mut u = self.base_u.to_vec();
        let mut p = self.base_p.to_vec();
        let nu = 3 * u.len();
        for (&d, &value) in self.prepared.free.iter().zip(x) {
            let physical = value * self.prepared.variable_scale(d);
            if d < nu {
                u[d / 3][d % 3] = physical;
            } else {
                p[d - nu] = physical;
            }
        }
        self.prepared.validate_state(&u, &p)?;
        Ok((u, p))
    }

    fn ensure_linearization(&self, x: &[f64]) -> Result<(), MixedUpError> {
        checkpoint(self.cx)?;
        if self
            .linearization
            .borrow()
            .as_ref()
            .is_some_and(|cached| cached.x == x)
        {
            return Ok(());
        }
        // The previous Newton point will not be queried again. Release its
        // matrix/factor before allocating the next numeric factor and scratch.
        *self.linearization.borrow_mut() = None;
        let (u, p) = self.physical_state(x)?;
        let evaluation =
            self.prepared
                .evaluate_with_inertia(&u, &p, true, self.inertia, self.cx)?;
        let full = evaluation.tangent.expect("requested original tangent");
        let nu = 3 * u.len();
        if let Some((node, _)) = self.prepared.problem.pressure_gauge {
            // A pressure pin is only a gauge when its constant component is a
            // true null mode on the free displacement space, at this actual F.
            for &d in &self.prepared.free[..self.prepared.free_displacements] {
                let (cols, values) = full.row(d);
                let mut signed = 0.0;
                let mut absolute = 0.0;
                for (&col, &value) in cols.iter().zip(values) {
                    if col >= nu && self.prepared.gauge_component[col - nu] {
                        signed += value;
                        absolute += value.abs();
                    }
                }
                if signed.abs() > 4096.0 * f64::EPSILON * absolute.max(f64::MIN_POSITIVE) {
                    return Err(MixedUpError::IncompatiblePressureGauge { node });
                }
            }
        }
        let n = self.prepared.free.len();
        let mut coo = Coo::new(n, n);
        for (i, &d) in self.prepared.free.iter().enumerate() {
            if i % 128 == 0 {
                checkpoint(self.cx)?;
            }
            let row_scale = self.prepared.variable_scale(d) / self.prepared.work_scale_j;
            let (cols, values) = full.row(d);
            for (&col, &value) in cols.iter().zip(values) {
                if let Some(j) = self.prepared.free_index[col] {
                    coo.push(i, j, row_scale * value * self.prepared.variable_scale(col));
                }
            }
        }
        let matrix = coo.assemble();
        if self.symbolic.borrow().is_none() {
            let analyzed = SymbolicLdlt::analyze(&matrix, DirectOrdering::Amd)
                .map_err(MixedUpError::Direct)?;
            *self.symbolic.borrow_mut() = Some(analyzed);
        }
        checkpoint(self.cx)?;
        let options = LdltOptions {
            max_factor_bytes: Some(self.prepared.problem.settings.maximum_factor_bytes),
            ..LdltOptions::default()
        };
        let factor = self
            .symbolic
            .borrow()
            .as_ref()
            .expect("analyzed fixed pattern")
            .factor(&matrix, &options)
            .map_err(MixedUpError::Direct)?;
        checkpoint(self.cx)?;
        let expected = Inertia {
            positive: self.prepared.free_displacements,
            negative: n - self.prepared.free_displacements,
        };
        if factor.inertia() != expected {
            return Err(MixedUpError::KktInertia {
                expected,
                got: factor.inertia(),
            });
        }
        *self.linearization.borrow_mut() = Some(Linearization {
            x: x.to_vec(),
            matrix,
            factor,
        });
        Ok(())
    }

    fn refuse_callback(&self, output: &mut [f64], error: MixedUpError) {
        output.fill(f64::NAN);
        *self.fatal_error.borrow_mut() = Some(error);
    }
}

impl<L: MixedFiniteStrainPotential> NonlinearProblem for Endpoint<'_, '_, '_, '_, L> {
    fn dimension(&self) -> usize {
        self.prepared.free.len()
    }

    fn residual(&self, x: &[f64], residual: &mut [f64]) {
        let result = self.physical_state(x).and_then(|(u, p)| {
            self.prepared
                .evaluate_with_inertia(&u, &p, false, self.inertia, self.cx)
        });
        match result {
            Ok(state) => {
                let nu = state.internal_forces_n.len();
                for (out, &d) in residual.iter_mut().zip(&self.prepared.free) {
                    let physical = if d < nu {
                        state.internal_forces_n[d] - self.load * self.prepared.external_n[d]
                    } else {
                        state.pressure_residual_m3[d - nu]
                    };
                    *out = physical * self.prepared.variable_scale(d) / self.prepared.work_scale_j;
                }
                *self.trial_error.borrow_mut() = None;
            }
            Err(error) => {
                residual.fill(f64::NAN);
                *self.trial_error.borrow_mut() = Some(error);
            }
        }
    }

    fn jacobian_apply(&self, x: &[f64], direction: &[f64], output: &mut [f64]) {
        match self.ensure_linearization(x) {
            Ok(()) => self
                .linearization
                .borrow()
                .as_ref()
                .expect("cached tangent")
                .matrix
                .spmv(direction, output),
            Err(error) => self.refuse_callback(output, error),
        }
    }

    fn preconditioner_apply(
        &self,
        x: &[f64],
        _outer: usize,
        _inner: usize,
        residual: &[f64],
        output: &mut [f64],
    ) {
        match self.ensure_linearization(x) {
            Ok(()) => output.copy_from_slice(
                &self
                    .linearization
                    .borrow()
                    .as_ref()
                    .expect("cached factor")
                    .factor
                    .solve(residual),
            ),
            Err(error) => self.refuse_callback(output, error),
        }
    }
}
