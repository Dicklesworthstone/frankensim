//! Retained small-strain membrane and bending plasticity on CST/DKT shells.
//!
//! Geometry, physical axial rotations and the three midside curvature points
//! are the original `fs-plate` operators. Each surface point owns a complete
//! through-thickness `fs-material` section. Its native plane-stress return,
//! transverse history and consistent tangent drive both force and stiffness;
//! yielding is not represented by a scalar reduction of the elastic matrix.
//!
//! The common solid increment solver owns Newton, line search and publication.
//! Every proposal starts from frozen accepted history, including timed laws.
//! Coordinates are internally `(u, length * theta)` and conjugate loads are
//! `(force, moment / length)`, so residuals compare forces in newtons. This is
//! a work-preserving change of coordinates, not a physical stiffness change.
//!
//! This lane admits small strains and small rotations, dead nodal loading and
//! facet-frame material axes. It does not model finite-strain forming, resolved
//! transverse shear or evolving contact. The generic solve alone has no material
//! heat account; [`energy`] audits the native energetic laws. The drilling penalty
//! remains a separate numerical elastic energy and never changes plastic state.

use crate::plastic3::increment::{self, IncrementControls, IncrementOperator};
use crate::plastic3::{IncrementSettings, SmallStrainTetError};
use fs_exec::Cx;
use fs_material::plastic::reduction::PlaneStressWorkspace;
use fs_material::plastic::section::{LayeredPlaneStressSection, LayeredPlaneStressState};
use fs_material::{SmallStrainLaw, Voigt};
use fs_plate::shell::{FacetGeometry, ShellMesh};
use fs_sparse::{Coo, Csr};
use std::cell::RefCell;

/// Native physical/numerical energy and drilling work for a single accepted step.
pub mod energy;

/// Retained reference inertia and native momentum evolution of plastic shells.
pub mod dynamics;

/// Immutable numerical drilling stiffness across changing material coefficients.
pub mod drilling;

/// Original immutable section assignments; no material averaging across facets.
#[derive(Debug, Clone, Copy)]
pub enum PlasticShellSections<'a, L: SmallStrainLaw> {
    /// One thickness stack for the whole midsurface.
    Uniform(&'a LayeredPlaneStressSection<L>),
    /// Exactly one thickness stack per triangle, in connectivity order.
    PerElement(&'a [LayeredPlaneStressSection<L>]),
    /// Borrow independently owned stacks without cloning their material laws.
    PerElementBorrowed(&'a [&'a LayeredPlaneStressSection<L>]),
}
impl<L: SmallStrainLaw> PlasticShellSections<'_, L> {
    fn at(&self, element: usize) -> &LayeredPlaneStressSection<L> {
        match self {
            Self::Uniform(section) => section,
            Self::PerElement(sections) => &sections[element],
            Self::PerElementBorrowed(sections) => sections[element],
        }
    }
}

/// Explicit allocation and physical-path limits, checked before preparation.
#[derive(Debug, Clone, Copy)]
pub struct PlasticShellBudget {
    /// Midsurface vertices.
    pub maximum_nodes: usize,
    /// CST/DKT triangles.
    pub maximum_elements: usize,
    /// Total material histories, including all three surface points per facet.
    pub maximum_material_points: usize,
    /// Unconstrained displacement and rotation coordinates.
    pub maximum_free_dofs: usize,
    /// Tangent entries before sparse duplicate assembly. Dynamic steps also
    /// count the retained mass entries added by the original inertia owner.
    pub maximum_tangent_entries: usize,
    /// Chronological load entries and their retained records.
    pub maximum_load_steps: usize,
}
impl Default for PlasticShellBudget {
    fn default() -> Self {
        Self {
            maximum_nodes: 100_000,
            maximum_elements: 100_000,
            maximum_material_points: 2_400_000,
            maximum_free_dofs: 600_000,
            maximum_tangent_entries: 32_400_000,
            maximum_load_steps: 10_000,
        }
    }
}

/// One quasistatic load path in the original six physical coordinates per node.
#[derive(Debug)]
pub struct SmallStrainShellProblem<'a, L: SmallStrainLaw> {
    /// Reference midsurface and oriented triangle connectivity, in metres.
    pub mesh: &'a ShellMesh,
    /// Native plane-stress material stacks in each facet's local axes.
    pub sections: PlasticShellSections<'a, L>,
    /// Positive numerical length converting rotations to solver coordinates [m].
    pub coordinate_length_m: f64,
    /// Maximum physical rotation and transverse displacement slope [rad],
    /// strictly positive and at most 0.05. Material strain limits are separate.
    pub rotation_limit_rad: f64,
    /// Nonnegative numerical drilling-spin penalty relative to virgin A66.
    /// The original elastic shell uses `fs_plate::shell::DRILLING_ALPHA`.
    /// Retain [`Self::prepare_reference_drilling`] when material coefficients evolve.
    pub drilling_alpha: f64,
    /// Prescribed `(6 * node + component, value)`: translations [m], rotations [rad].
    pub prescribed: &'a [(usize, f64)],
    /// Dead nodal `(6 * node + component, value)`: forces [N], moments [N m].
    pub nodal_loads: &'a [(usize, f64)],
    /// Absolute factors scaling both reference fields, in chronological order.
    pub load_path: &'a [f64],
    /// Geometry, memory and path envelope.
    pub budget: PlasticShellBudget,
    /// Original solid Newton/backtracking and sparse solve controls.
    pub settings: IncrementSettings,
}

/// Physical reactions and convergence evidence at one accepted load entry.
#[derive(Debug, Clone, PartialEq)]
pub struct PlasticShellIncrement {
    /// Absolute load factor reached.
    pub load_factor: f64,
    /// Dimensionless force residuals in the declared scaled coordinates.
    pub residual_history: Vec<f64>,
    /// Prescribed-coordinate reactions: forces [N] and moments [N m], in
    /// `SmallStrainShellProblem::prescribed` order.
    pub reactions: Vec<f64>,
    /// Facets whose underlying material history changed.
    pub evolving_elements: usize,
    /// Corrections using the original solid solver's indefinite fallback.
    pub indefinite_fallbacks: usize,
}

/// Accepted endpoint. Returning an error never mutates a supplied checkpoint.
#[derive(Debug, Clone, PartialEq)]
pub struct SmallStrainShellSolution<S> {
    /// Original nodal `(ux,uy,uz,theta_x,theta_y,theta_z)`, in metres/radians.
    pub coordinates: Vec<[f64; 6]>,
    /// Facet-major, then original DKT midside order; each section retains all
    /// layer histories, including transverse plastic strain and contraction.
    pub section_states: Vec<LayeredPlaneStressState<S>>,
    /// Same endpoint return as the committed histories: `(Nx,Ny,Nxy,Mx,My,Mxy)`
    /// in local facet axes, with forces per length [N/m] and moments [N].
    pub section_resultants: Vec<[f64; 6]>,
    /// Physical path records, retaining unloading and reversal entries.
    pub increments: Vec<PlasticShellIncrement>,
}

fn invalid(what: impl Into<String>) -> SmallStrainTetError {
    SmallStrainTetError::InvalidInput { what: what.into() }
}
fn checkpoint(cx: &Cx<'_>) -> Result<(), SmallStrainTetError> {
    cx.checkpoint().map_err(|_| SmallStrainTetError::Cancelled)
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).fold(0.0, |s, (x, y)| x.mul_add(*y, s))
}

struct Element {
    nodes: [usize; 3],
    geometry: FacetGeometry,
    /// Generalized engineering strain from scaled global nodal coordinates.
    strain: [[[f64; 18]; 6]; 3],
    drilling: [[f64; 18]; 3],
    drilling_weight: f64,
}

/// Transform one local element row into scaled global nodal coordinates.
fn transform(row: &[f64; 18], frame: &[[f64; 3]; 3], length: f64) -> [f64; 18] {
    core::array::from_fn(|d| {
        let node = d / 6;
        let component = d % 6;
        let offset = if component < 3 { 0 } else { 3 };
        let value: f64 = (0..3)
            .map(|axis| row[6 * node + offset + axis] * frame[axis][component % 3])
            .sum();
        value / if component < 3 { 1.0 } else { length }
    })
}

impl<L: SmallStrainLaw> SmallStrainShellProblem<'_, L> {
    fn prepare(&self, cx: &Cx<'_>) -> Result<Vec<Element>, SmallStrainTetError> {
        self.prepare_with_reference_drilling(None, cx)
    }

    fn prepare_with_reference_drilling(
        &self,
        reference: Option<&drilling::PlasticShellReferenceDrilling<'_>>,
        cx: &Cx<'_>,
    ) -> Result<Vec<Element>, SmallStrainTetError> {
        checkpoint(cx)?;
        let n = self.mesh.nodes.len();
        let ne = self.mesh.tris.len();
        let b = self.budget;
        let s = self.settings;
        if n < 3
            || ne == 0
            || n > b.maximum_nodes
            || ne > b.maximum_elements
            || self.load_path.is_empty()
            || self.load_path.len() > b.maximum_load_steps
            || !self.coordinate_length_m.is_finite()
            || self.coordinate_length_m <= 0.0
            || !self.rotation_limit_rad.is_finite()
            || self.rotation_limit_rad <= 0.0
            || self.rotation_limit_rad > 0.05
            || !self.drilling_alpha.is_finite()
            || self.drilling_alpha < 0.0
            || s.max_newton_iterations == 0
            || s.max_linear_iterations == 0
            || !(s.relative_tolerance > 0.0 && s.relative_tolerance < 1.0)
            || !(s.linear_tolerance > 0.0 && s.linear_tolerance < 1.0)
            || self
                .load_path
                .iter()
                .chain(self.mesh.nodes.iter().flatten())
                .any(|x| !x.is_finite())
        {
            return Err(invalid(
                "invalid shell geometry, controls, load path or budget",
            ));
        }
        let section_count = match &self.sections {
            PlasticShellSections::Uniform(_) => ne,
            PlasticShellSections::PerElement(sections) => sections.len(),
            PlasticShellSections::PerElementBorrowed(sections) => sections.len(),
        };
        if section_count != ne {
            return Err(invalid("shell section count must equal triangle count"));
        }
        let dofs = n
            .checked_mul(6)
            .ok_or_else(|| invalid("shell DOF count overflow"))?;
        let entries = ne
            .checked_mul(324)
            .ok_or_else(|| invalid("shell tangent count overflow"))?;
        if entries > b.maximum_tangent_entries {
            return Err(invalid("shell tangent exceeds allocation budget"));
        }
        let mut fixed = vec![false; dofs];
        for &(d, v) in self.prescribed {
            if d >= dofs || !v.is_finite() || fixed[d] {
                return Err(invalid("invalid or repeated prescribed shell coordinate"));
            }
            fixed[d] = true;
        }
        if dofs - self.prescribed.len() > b.maximum_free_dofs
            || self
                .nodal_loads
                .iter()
                .any(|&(d, v)| d >= dofs || !v.is_finite())
        {
            return Err(invalid("invalid shell loads or free-coordinate budget"));
        }
        let mut material_points = 0usize;
        for e in 0..ne {
            material_points = self
                .sections
                .at(e)
                .point_count()
                .checked_mul(3)
                .and_then(|q| material_points.checked_add(q))
                .ok_or_else(|| invalid("shell history count overflow"))?;
        }
        if material_points > b.maximum_material_points {
            return Err(invalid(
                "shell thickness histories exceed allocation budget",
            ));
        }
        if let Some(reference) = reference {
            reference.admit(self, cx)?;
        }
        let mut workspace = if reference.is_none() && self.drilling_alpha != 0.0 {
            Some(PlaneStressWorkspace::new().map_err(|e| invalid(e.to_string()))?)
        } else {
            None
        };
        let mut elements = Vec::with_capacity(ne);
        for e in 0..ne {
            checkpoint(cx)?;
            let g = self.mesh.facet(e).map_err(|e| invalid(e.to_string()))?;
            let (curvature, _) = fs_plate::dkt_curvature_quadrature(&g.x, &g.y, e)
                .map_err(|e| invalid(e.to_string()))?;
            let drilling_weight = if let Some(reference) = reference {
                reference.weight(e)
            } else if let Some(workspace) = &mut workspace {
                let section = self.sections.at(e);
                let virgin = section
                    .evaluate_with_workspace([0.0; 6], &section.initial_state(), workspace)
                    .map_err(|error| {
                        invalid(format!("shell element {e} virgin section: {error}"))
                    })?;
                self.drilling_alpha * virgin.tangent[2][2] * g.area_m2 / 3.0
            } else {
                0.0
            };
            if !drilling_weight.is_finite() || drilling_weight < 0.0 {
                return Err(invalid("invalid virgin shell drilling stiffness"));
            }
            let mut strain = [[[0.0; 18]; 6]; 3];
            let slots = [2, 4, 3, 8, 10, 9, 14, 16, 15];
            let signs = [1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0];
            for (q, rows) in strain.iter_mut().enumerate() {
                for a in 0..3 {
                    let [dx, dy] = g.gradient[a];
                    rows[0][6 * a] = dx;
                    rows[1][6 * a + 1] = dy;
                    rows[2][6 * a] = dy;
                    rows[2][6 * a + 1] = dx;
                }
                for p in 0..3 {
                    for j in 0..9 {
                        rows[3 + p][slots[j]] = curvature[q][p][j] * signs[j];
                    }
                }
                for row in rows {
                    *row = transform(row, &g.frame, self.coordinate_length_m);
                }
            }
            let mut drilling = [[0.0; 18]; 3];
            for (a, row) in drilling.iter_mut().enumerate() {
                row[6 * a + 5] = 1.0;
                for j in 0..3 {
                    row[6 * j] = 0.5 * g.gradient[j][1];
                    row[6 * j + 1] = -0.5 * g.gradient[j][0];
                }
                *row = transform(row, &g.frame, self.coordinate_length_m);
            }
            if strain
                .iter()
                .flatten()
                .flatten()
                .chain(drilling.iter().flatten())
                .any(|x| !x.is_finite())
            {
                return Err(invalid("shell strain operator is not representable"));
            }
            elements.push(Element {
                nodes: self.mesh.tris[e],
                geometry: g,
                strain,
                drilling,
                drilling_weight,
            });
        }
        Ok(elements)
    }

    /// Solve from virgin material history. Each load entry advances physical
    /// time once if the selected material law carries a time interval.
    pub fn solve(
        &self,
        cx: &Cx<'_>,
    ) -> Result<SmallStrainShellSolution<L::State>, SmallStrainTetError> {
        self.solve_initial(None, cx)
    }

    /// Continue from the same reference mesh, section assignments and history
    /// ordering. Callers own any material/coefficient change and its energy
    /// implications. Inputs remain unchanged even if a later path entry fails.
    /// [`Self::solve_from_with_reference_drilling`] keeps the original numerical
    /// gauge unchanged when refreshing temperature-dependent material laws.
    pub fn solve_from(
        &self,
        coordinates: &[[f64; 6]],
        states: &[LayeredPlaneStressState<L::State>],
        cx: &Cx<'_>,
    ) -> Result<SmallStrainShellSolution<L::State>, SmallStrainTetError> {
        self.solve_initial(Some((coordinates, states)), cx)
    }

    fn solve_initial(
        &self,
        initial: Option<(&[[f64; 6]], &[LayeredPlaneStressState<L::State>])>,
        cx: &Cx<'_>,
    ) -> Result<SmallStrainShellSolution<L::State>, SmallStrainTetError> {
        let elements = self.prepare(cx)?;
        self.solve_prepared(&elements, initial, None, cx)
    }

    fn solve_prepared(
        &self,
        elements: &[Element],
        initial: Option<(&[[f64; 6]], &[LayeredPlaneStressState<L::State>])>,
        inertia: Option<&crate::inertia::BackwardEulerInertia<'_>>,
        cx: &Cx<'_>,
    ) -> Result<SmallStrainShellSolution<L::State>, SmallStrainTetError> {
        let length = self.coordinate_length_m;
        let prescribed: Vec<_> = self
            .prescribed
            .iter()
            .map(|&(d, v)| (d, v * if d % 6 < 3 { 1.0 } else { length }))
            .collect();
        let loads: Vec<_> = self
            .nodal_loads
            .iter()
            .map(|&(d, v)| (d, v / if d % 6 < 3 { 1.0 } else { length }))
            .collect();
        let coordinates = if let Some((coordinates, states)) = initial {
            if coordinates.len() != self.mesh.nodes.len() || states.len() != 3 * elements.len() {
                return Err(invalid(
                    "shell checkpoint has incompatible node or section count",
                ));
            }
            coordinates
                .iter()
                .flat_map(|q| {
                    q.iter()
                        .enumerate()
                        .map(move |(c, v)| v * if c < 3 { 1.0 } else { length })
                })
                .collect::<Vec<_>>()
        } else {
            Vec::new()
        };
        if prescribed.iter().chain(&loads).any(|(_, x)| !x.is_finite())
            || coordinates.iter().any(|x| !x.is_finite())
        {
            return Err(invalid("scaled shell coordinates or forces overflow"));
        }
        let operator = Operator {
            problem: self,
            elements,
            inertia,
            workspace: RefCell::new(
                PlaneStressWorkspace::new().map_err(|e| invalid(e.to_string()))?,
            ),
        };
        if let Some((_, states)) = initial {
            for (e, element) in elements.iter().enumerate() {
                checkpoint(cx)?;
                let local = operator.local_coordinates(element, &coordinates)?;
                for (q, rows) in element.strain.iter().enumerate() {
                    let state = &states[3 * e + q];
                    self.sections.at(e).validate_state(state).map_err(|error| {
                        invalid(format!("shell checkpoint element {e}, point {q}: {error}"))
                    })?;
                    for (p, row) in rows.iter().enumerate() {
                        let expected = dot(row, &local);
                        let previous = state.generalized_strain[p];
                        let scale = if p < 3 { 1.0 } else { length };
                        let allowance = 256.0
                            * f64::EPSILON
                            * (scale * expected.abs())
                                .max(scale * previous.abs())
                                .max(1.0);
                        if scale * (expected - previous).abs() > allowance {
                            return Err(invalid(
                                "shell checkpoint strain does not belong to its geometry",
                            ));
                        }
                    }
                }
            }
        }
        let result = increment::solve(
            &operator,
            IncrementControls {
                dofs: 6 * self.mesh.nodes.len(),
                prescribed_m: &prescribed,
                nodal_forces_n: &loads,
                load_path: self.load_path,
                settings: self.settings,
            },
            initial.map(|(_, states)| (coordinates.as_slice(), states)),
            inertia,
            cx,
        )?;
        let coordinates = result
            .displacement_m
            .chunks_exact(6)
            .map(|q| core::array::from_fn(|c| q[c] / if c < 3 { 1.0 } else { length }))
            .collect();
        let increments = result
            .increments
            .into_iter()
            .map(|r| PlasticShellIncrement {
                load_factor: r.load_factor,
                residual_history: r.residual_history,
                reactions: r
                    .reactions_n
                    .into_iter()
                    .zip(self.prescribed)
                    .map(|(v, &(d, _))| v * if d % 6 < 3 { 1.0 } else { length })
                    .collect(),
                evolving_elements: r.evolving_elements,
                indefinite_fallbacks: r.indefinite_fallbacks,
            })
            .collect();
        Ok(SmallStrainShellSolution {
            coordinates,
            section_states: result.point_states,
            section_resultants: result.point_stress,
            increments,
        })
    }
}

struct Operator<'a, 'm, L: SmallStrainLaw> {
    problem: &'a SmallStrainShellProblem<'m, L>,
    elements: &'a [Element],
    inertia: Option<&'a crate::inertia::BackwardEulerInertia<'a>>,
    workspace: RefCell<PlaneStressWorkspace>,
}

impl<L: SmallStrainLaw> Operator<'_, '_, L> {
    fn local_coordinates(
        &self,
        element: &Element,
        u: &[f64],
    ) -> Result<[f64; 18], SmallStrainTetError> {
        // Remove a common translation before differentiating. Rotations are
        // physical axial vectors; a constant rotation is not removed alone.
        let values = core::array::from_fn(|d| {
            let c = d % 6;
            u[6 * element.nodes[d / 6] + c]
                - if c < 3 {
                    u[6 * element.nodes[0] + c]
                } else {
                    0.0
                }
        });
        let mut gradient = [[0.0; 2]; 3];
        for a in 0..3 {
            let node = element.nodes[a];
            let rotation = (3..6).fold(0.0_f64, |n, c| {
                n.hypot(u[6 * node + c] / self.problem.coordinate_length_m)
            });
            if !rotation.is_finite() || rotation > self.problem.rotation_limit_rad {
                return Err(invalid("shell rotation exceeds its small-rotation domain"));
            }
            for c in 0..3 {
                let local = dot(&element.geometry.frame[c], &values[6 * a..6 * a + 3]);
                for j in 0..2 {
                    gradient[c][j] += local * element.geometry.gradient[a][j];
                }
            }
        }
        let spin = 0.5 * (gradient[1][0] - gradient[0][1]);
        let slope = gradient[2][0].hypot(gradient[2][1]);
        if !spin.is_finite()
            || !slope.is_finite()
            || spin.abs() > self.problem.rotation_limit_rad
            || slope > self.problem.rotation_limit_rad
        {
            return Err(invalid(
                "shell displacement gradient exceeds its small-rotation domain",
            ));
        }
        Ok(values)
    }

    fn evaluate(
        &self,
        states: &[LayeredPlaneStressState<L::State>],
        u: &[f64],
        prescribed: Option<&[bool]>,
        tangent: bool,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, f64, Option<Csr>), SmallStrainTetError> {
        let mut force = vec![0.0; u.len()];
        let mut bounds = vec![0.0; u.len()];
        let mut matrix = tangent.then(|| Coo::new(u.len(), u.len()));
        let mut workspace = self.workspace.borrow_mut();
        for (e, element) in self.elements.iter().enumerate() {
            checkpoint(cx)?;
            let local = self.local_coordinates(element, u)?;
            let mut stiffness = [0.0; 324];
            let mut element_force = [0.0; 18];
            let mut element_bound = [0.0; 18];
            for (q, b) in element.strain.iter().enumerate() {
                checkpoint(cx)?;
                let strain = core::array::from_fn(|i| dot(&b[i], &local));
                let response = self
                    .problem
                    .sections
                    .at(e)
                    .evaluate_with_workspace(strain, &states[3 * e + q], &mut workspace)
                    .map_err(|error| {
                        invalid(format!("shell element {e}, surface point {q}: {error}"))
                    })?;
                let weight = element.geometry.area_m2 / 3.0;
                for i in 0..18 {
                    for p in 0..6 {
                        let term = weight * b[p][i] * response.resultants[p];
                        element_force[i] += term;
                        element_bound[i] += term.abs();
                    }
                }
                if tangent {
                    let cb: [[f64; 18]; 6] = core::array::from_fn(|p| {
                        core::array::from_fn(|j| {
                            (0..6).map(|r| response.tangent[p][r] * b[r][j]).sum()
                        })
                    });
                    for i in 0..18 {
                        for j in 0..18 {
                            stiffness[18 * i + j] +=
                                weight * (0..6).map(|p| b[p][i] * cb[p][j]).sum::<f64>();
                        }
                    }
                }
            }
            for spin in &element.drilling {
                let effort = element.drilling_weight * dot(spin, &local);
                for i in 0..18 {
                    let term = spin[i] * effort;
                    element_force[i] += term;
                    element_bound[i] += term.abs();
                    if tangent {
                        for j in 0..18 {
                            stiffness[18 * i + j] += element.drilling_weight * spin[i] * spin[j];
                        }
                    }
                }
            }
            for i in 0..18 {
                let row = 6 * element.nodes[i / 6] + i % 6;
                force[row] += element_force[i];
                bounds[row] += element_bound[i];
                if let Some(matrix) = matrix.as_mut() {
                    for j in 0..18 {
                        let col = 6 * element.nodes[j / 6] + j % 6;
                        matrix.push(row, col, stiffness[18 * i + j]);
                    }
                }
            }
            if stiffness
                .iter()
                .chain(&element_force)
                .chain(&element_bound)
                .any(|x| !x.is_finite())
            {
                return Err(invalid("shell force or tangent is not representable"));
            }
        }
        let force_scale = bounds
            .iter()
            .enumerate()
            .filter(|(d, _)| prescribed.is_none_or(|p| !p[*d]))
            .fold(0.0_f64, |n, (_, &v)| n.hypot(v));
        if !force_scale.is_finite() || force.iter().any(|x| !x.is_finite()) {
            return Err(invalid("assembled shell force is not representable"));
        }
        if let (Some(inertia), Some(matrix)) = (self.inertia, matrix.as_mut()) {
            inertia.add_tangent(matrix, cx)?;
        }
        Ok((force, force_scale, matrix.map(|coo| coo.assemble())))
    }
}

impl<L: SmallStrainLaw> IncrementOperator for Operator<'_, '_, L> {
    type State = LayeredPlaneStressState<L::State>;
    fn initial_states(&self) -> Vec<Self::State> {
        (0..self.elements.len())
            .flat_map(|e| (0..3).map(move |_| self.problem.sections.at(e).initial_state()))
            .collect()
    }
    fn components_per_node(&self) -> usize {
        6
    }
    fn deformation_coordinate(&self, u: &[f64], row: usize, column: usize) -> f64 {
        let c = column % 6;
        u[column] - if c < 3 { u[6 * (row / 6) + c] } else { 0.0 }
    }
    fn internal(
        &self,
        states: &[Self::State],
        u: &[f64],
        prescribed: Option<&[bool]>,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, f64), SmallStrainTetError> {
        let (force, scale, _) = self.evaluate(states, u, prescribed, false, cx)?;
        Ok((force, scale))
    }
    fn tangent(
        &self,
        states: &[Self::State],
        u: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Csr, SmallStrainTetError> {
        let (_, _, matrix) = self.evaluate(states, u, None, true, cx)?;
        Ok(matrix.expect("requested shell tangent"))
    }
    fn tangent_symmetric(&self) -> bool {
        (0..self.elements.len()).all(|e| self.problem.sections.at(e).tangent_symmetric())
    }
    fn commit(
        &self,
        states: &mut [Self::State],
        u: &[f64],
        final_stress: Option<&mut Vec<Voigt>>,
        cx: &Cx<'_>,
    ) -> Result<usize, SmallStrainTetError> {
        let mut workspace = self.workspace.borrow_mut();
        let mut resultants = final_stress;
        if let Some(out) = resultants.as_mut() {
            out.clear();
        }
        let mut changed = 0;
        for (e, element) in self.elements.iter().enumerate() {
            checkpoint(cx)?;
            let local = self.local_coordinates(element, u)?;
            let mut evolving = false;
            for (q, b) in element.strain.iter().enumerate() {
                checkpoint(cx)?;
                let index = 3 * e + q;
                let strain = core::array::from_fn(|i| dot(&b[i], &local));
                let response = self
                    .problem
                    .sections
                    .at(e)
                    .evaluate_with_workspace(strain, &states[index], &mut workspace)
                    .map_err(|error| {
                        invalid(format!("shell element {e}, surface point {q}: {error}"))
                    })?;
                evolving |= states[index]
                    .points
                    .iter()
                    .zip(&response.state.points)
                    .any(|(old, new)| old.history != new.history);
                if let Some(out) = resultants.as_mut() {
                    out.push(response.resultants);
                }
                states[index] = response.state;
            }
            changed += usize::from(evolving);
        }
        Ok(changed)
    }
    fn recoverable_trial_error(&self, error: &SmallStrainTetError) -> bool {
        matches!(error, SmallStrainTetError::InvalidInput { .. })
    }
}
