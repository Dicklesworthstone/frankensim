//! The transient finite-element body on the composition contract.
//!
//! [`ConductionBody`] puts one [`BackwardEuler`] owner, with its mesh, its
//! conductivity and its boundary partition, behind the Component contract.
//! Each named boundary patch of the body is one heat port, of one of two
//! kinds ([`PatchPort`]).
//!
//! A *mean* port is scalar: it receives the heat rate `Q` [W] into the body
//! and applies it on its faces as the uniform outward flux `-Q / A`, through
//! the owner's own Neumann rows. It returns the patch's area mean of the END
//! temperature of the step, the temperature at which a uniform flux does its
//! work (`∫ T q dA = Q T̄`) and backward Euler's own reading of it.
//!
//! A *trace* port is a field heat port over the patch's vertices: it
//! receives one heat load `Q_i` [W] per vertex, already integrated against
//! that vertex's P1 function, and returns the END temperature at each. The
//! loads reach the owner as the nodal outward flux `q` with `M q = -Q`, `M`
//! the patch's consistent face mass, so the owner's own Neumann rows put
//! `Q_i` on vertex `i` again. Two bodies meet through
//! [`ConductionBody::matching_contact`] where their traces match and through
//! [`ConductionBody::mapped_contact`] where they do not: either states the
//! owner's contact between them as the quadrature points of the traces'
//! common refinement, each with its weights on both traces and its
//! conductance.
//!
//! The step is the owner's unchanged `advance`; the tangent is the owner's
//! `linearize_step` with no Robin port, driven by the same nodal loads.
//!
//! No claim: the entropy generated inside the body (conduction between nodes
//! at different temperatures) is reported as zero, a lower bound. Robin,
//! Dirichlet and contact traces, volumetric sources, temperature-dependent
//! conductivity and sheet conduction are refused, not wrapped. Traces that
//! are not coplanar, or do not cover each other, are refused. The state has
//! no snapshot codec.

use std::collections::{BTreeMap, BTreeSet};
use std::num::NonZeroUsize;

use fs_exec::Cx;
use fs_ports::StableId;
use fs_ports::component::{Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt};

use crate::interface::InterfaceResistance;
use crate::interface::nonmatching::{
    NonmatchingOptions, TraceFace, degree_two_weights, trace_overlaps,
};
use crate::transient::backward_euler::{BackwardEuler, StepConfig, StepLinearization};
use crate::{ConductionError, ConductionProblem, ScalarField, ThermalBc, ThermalBoundary};

/// How one named boundary patch of a [`ConductionBody`] meets the kernel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PatchPort<'a> {
    /// A scalar heat port on the named region: a heat rate in, spread as a
    /// uniform flux, and the patch's area-mean end temperature out.
    Mean(&'a str),
    /// A field heat port on the named region: one integrated heat load per
    /// vertex of the patch in, the end temperature at each vertex out, the
    /// vertices in ascending mesh order.
    Trace(&'a str),
}

impl<'a> PatchPort<'a> {
    /// The boundary region this port sits on.
    #[must_use]
    pub const fn region(self) -> &'a str {
        match self {
            Self::Mean(name) | Self::Trace(name) => name,
        }
    }
}

/// The points of a trace port and the solve that turns their loads into the
/// owner's nodal flux.
#[derive(Debug)]
struct Trace {
    /// Mesh vertex of each point, ascending.
    vertices: Vec<usize>,
    /// The patch's faces over the points (indices into `vertices`), in
    /// boundary-slot order.
    faces: Vec<[usize; 3]>,
    /// The unit normal of each face, out of the body.
    normals: Vec<[f64; 3]>,
    /// Lower Cholesky factor of the patch's consistent P1 face mass [m²],
    /// row-major over the points.
    mass_factor: Vec<f64>,
}

impl Trace {
    /// The outward flux `q` [W/m²] at the points with `M q = -loads`.
    fn outward_flux(&self, loads: &[f64]) -> Vec<f64> {
        let n = self.vertices.len();
        let l = &self.mass_factor;
        let mut q: Vec<f64> = loads.iter().map(|load| -load).collect();
        for i in 0..n {
            let mut sum = q[i];
            for j in 0..i {
                sum -= l[i * n + j] * q[j];
            }
            q[i] = sum / l[i * n + i];
        }
        for i in (0..n).rev() {
            let mut sum = q[i];
            for j in (i + 1)..n {
                sum -= l[j * n + i] * q[j];
            }
            q[i] = sum / l[i * n + i];
        }
        q
    }
}

/// One port's boundary patch: its region, its area and its faces.
#[derive(Debug)]
struct Patch {
    region: usize,
    area_m2: f64,
    /// Vertices and area of every face of the patch, in boundary-slot order.
    faces: Vec<([usize; 3], f64)>,
    /// First entry of the port in the body's inputs and in its outputs.
    offset: usize,
    /// The points of a trace port; `None` for a mean port.
    trace: Option<Trace>,
}

impl Patch {
    /// Entries of the port in the body's inputs and in its outputs.
    fn width(&self) -> usize {
        self.trace.as_ref().map_or(1, |trace| trace.vertices.len())
    }
}

/// One quadrature point of a contact between two traces: where it lies on
/// either trace, as up to three points of the trace and their P1 weights
/// there, and its conductance, its quadrature weight over `R''`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContactPoint {
    /// Points of side A's trace and their weights at this point.
    pub a: [(usize, f64); 3],
    /// Points of side B's trace and their weights at this point.
    pub b: [(usize, f64); 3],
    /// The conductance of the contact at this point [W/K].
    pub conductance_w_k: f64,
}

/// The owner's contact between two trace ports, stated as the interface
/// points a `ConductanceMatrix` of fs-system couples two field heat ports
/// through: the degree-two quadrature of `(T_a - T_b) (v_a - v_b) / R''` over
/// the traces' common refinement.
///
/// Built by [`ConductionBody::matching_contact`] and
/// [`ConductionBody::mapped_contact`].
#[derive(Debug, Clone, PartialEq)]
pub struct TraceContact {
    /// Points of side A's trace.
    pub width_a: usize,
    /// Points of side B's trace.
    pub width_b: usize,
    /// Three quadrature points for every triangle of the common refinement.
    pub points: Vec<ContactPoint>,
}

/// A backward-Euler finite-element body whose ports are named boundary patches.
///
/// Port `i` has the id of region `ports[i]` and [`PortRole::FlowIn`]. A mean
/// port is a [`PortSpec::heat`]: its input is the heat rate into the body
/// [W], its output the patch's area-mean end temperature [K]. A trace port
/// is a [`PortSpec::heat_field`] over the patch's vertices: its inputs are
/// the heat loads into the body at the vertices [W], its outputs their end
/// temperatures [K]. The state is the nodal temperature field [K] in vertex
/// order. The stored energy is `1ᵀ C T`, with the owner's prepared
/// capacitance `C`, so its change between two states is the owner's own
/// `stored_energy_change_j` up to rounding.
#[derive(Debug)]
pub struct ConductionBody<'m> {
    solver: BackwardEuler<'m>,
    problem: ConductionProblem<'m>,
    config: StepConfig,
    ports: Vec<PortSpec>,
    patches: Vec<Patch>,
    /// Entries of all ports together, in the inputs and in the outputs.
    width: usize,
    /// `C 1`: the stored energy of a field `T` is `Σ_v weight_v T_v`.
    weights: Vec<f64>,
}

impl<'m> ConductionBody<'m> {
    /// Put `solver` and `problem` behind heat ports on the regions `ports`.
    ///
    /// Every region of `problem.boundary` must be a uniform zero-flux Neumann
    /// row: a port's region because the port's heat replaces that zero, any
    /// other because heat through it would cross no port. The source must be
    /// `ScalarField::Uniform(0.0)` and the conductivity independent of
    /// temperature. `config` is the owner's step configuration, used for
    /// every step and tangent.
    ///
    /// # Errors
    /// [`ConductionError::Config`] naming the refused port, region, source or
    /// material; a solver prepared on another mesh.
    pub fn new(
        solver: BackwardEuler<'m>,
        problem: ConductionProblem<'m>,
        ports: &[PatchPort<'_>],
        config: StepConfig,
    ) -> Result<Self, ConductionError> {
        if !std::ptr::eq(solver.mesh(), problem.mesh) {
            return Err(invalid(
                "the solver and the problem must use the same mesh".into(),
            ));
        }
        let boundary = problem.boundary;
        for (name, condition) in boundary.region_names().iter().zip(boundary.conditions()) {
            if !insulated(condition) {
                return Err(invalid(format!(
                    "region {name} is a {} row that is not a uniform zero flux; heat through it would cross no port",
                    condition.tag()
                )));
            }
        }
        if boundary.has_temperature_dependent_surface_conduction() {
            return Err(invalid(
                "temperature-dependent sheet conduction is not wrapped".into(),
            ));
        }
        if !matches!(problem.source, ScalarField::Uniform(s) if *s == 0.0) {
            return Err(invalid(
                "a volumetric source would add heat through no port".into(),
            ));
        }
        let nonlinear = match problem.element_materials {
            Some(materials) => (0..problem.mesh.element_count()).try_fold(false, |any, e| {
                Ok::<_, ConductionError>(any || materials.model_for(e)?.is_temperature_dependent())
            })?,
            None => problem.material.is_temperature_dependent(),
        };
        if nonlinear {
            return Err(invalid(
                "temperature-dependent conductivity is not wrapped".into(),
            ));
        }
        let mut specs = Vec::with_capacity(ports.len());
        let mut patches = Vec::with_capacity(ports.len());
        let mut width = 0;
        for (index, &port) in ports.iter().enumerate() {
            let name = port.region();
            if ports[..index]
                .iter()
                .any(|earlier| earlier.region() == name)
            {
                return Err(invalid(format!("port {name} is named twice")));
            }
            let (spec, patch) = patch_of(&problem, port, width)?;
            width += patch.width();
            specs.push(spec);
            patches.push(patch);
        }
        if specs.is_empty() {
            return Err(invalid("a body needs at least one port".into()));
        }
        let ones = vec![1.0; problem.mesh.vertex_count()];
        let mut weights = vec![0.0; ones.len()];
        solver.capacitance().spmv(&ones, &mut weights);
        Ok(Self {
            solver,
            problem,
            config,
            ports: specs,
            patches,
            width,
            weights,
        })
    }

    /// Area [m²] of port `port`'s patch: the flux of a heat rate `Q` on it is
    /// `Q / patch_area_m2(port)`.
    #[must_use]
    pub fn patch_area_m2(&self, port: usize) -> f64 {
        self.patches[port].area_m2
    }

    /// The mesh vertices of trace port `port`, in the order of its inputs and
    /// outputs; `None` for a mean port.
    #[must_use]
    pub fn trace_vertices(&self, port: usize) -> Option<&[usize]> {
        self.patches[port]
            .trace
            .as_ref()
            .map(|trace| trace.vertices.as_slice())
    }

    /// The owner's matching-P1 contact between trace port `port` of this
    /// body (side A) and trace port `other_port` of `other` (side B), through
    /// the contact resistance `resistance`.
    ///
    /// `offset_m` places `other` in this body's frame: a point of B lies at
    /// its mesh position plus `offset_m`. The traces match when every point
    /// of A has the coordinates of exactly one point of B, bit for bit, and
    /// every face of A's patch is a face of B's with the opposite normal, as
    /// the owner's contact requires of its face pairs. Each face is its own
    /// overlap triangle: three quadrature points by the owner's rule for
    /// exact pairs, each a third of A's face area over the resistance. Bonded
    /// to the two ports through a `ConductanceMatrix`, they make the pair
    /// the owner's contact problem on one mesh.
    ///
    /// # Errors
    /// [`ConductionError::Config`] when a port is a mean port, the traces do
    /// not match point for point and face for face, or two matched faces do
    /// not face each other.
    pub fn matching_contact(
        &self,
        port: usize,
        other: &ConductionBody<'_>,
        other_port: usize,
        offset_m: [f64; 3],
        resistance: &InterfaceResistance,
    ) -> Result<TraceContact, ConductionError> {
        let (Some(a), Some(b)) = (
            self.patches[port].trace.as_ref(),
            other.patches[other_port].trace.as_ref(),
        ) else {
            return Err(invalid(
                "a matching contact joins two trace ports; a mean port has one temperature".into(),
            ));
        };
        let (n_a, n_b) = (a.vertices.len(), b.vertices.len());
        let key = |point: [f64; 3]| point.map(f64::to_bits);
        let mut places = BTreeMap::new();
        for (index, &vertex) in b.vertices.iter().enumerate() {
            let [x, y, z] = other.problem.mesh.positions()[vertex];
            let placed = [x + offset_m[0], y + offset_m[1], z + offset_m[2]];
            if places.insert(key(placed), index).is_some() {
                return Err(invalid(format!(
                    "two points of the second trace lie at {placed:?}"
                )));
            }
        }
        let mut partner = Vec::with_capacity(n_a);
        for &vertex in &a.vertices {
            let point = self.problem.mesh.positions()[vertex];
            let index = places.get(&key(point)).copied().ok_or_else(|| {
                invalid(format!(
                    "the traces do not match: no point of the second trace lies at {point:?}; the mapped transfer between traces that do not match is not built"
                ))
            })?;
            partner.push(index);
        }
        if n_a != n_b || partner.iter().collect::<BTreeSet<_>>().len() != n_a {
            return Err(invalid(format!(
                "the traces do not match: {n_a} points against {n_b}"
            )));
        }
        let sorted = |mut face: [usize; 3]| {
            face.sort_unstable();
            face
        };
        let faces_b: BTreeMap<[usize; 3], [f64; 3]> = b
            .faces
            .iter()
            .zip(&b.normals)
            .map(|(&face, &normal)| (sorted(face), normal))
            .collect();
        let mut shared = 0;
        for (face, normal_a) in a.faces.iter().zip(&a.normals) {
            let Some(normal_b) = faces_b.get(&sorted(face.map(|point| partner[point]))) else {
                continue;
            };
            let dot = normal_a[0].mul_add(
                normal_b[0],
                normal_a[1].mul_add(normal_b[1], normal_a[2] * normal_b[2]),
            );
            // The owner's own band for the opposing normals of a face pair.
            if dot > -1.0 + 64.0 * f64::EPSILON {
                return Err(invalid(format!(
                    "the traces do not face each other: two matched faces have normals {normal_a:?} and {normal_b:?}"
                )));
            }
            shared += 1;
        }
        if shared != a.faces.len() || a.faces.len() != b.faces.len() {
            return Err(invalid(format!(
                "the traces share their points but only {shared} of {} and {} faces",
                a.faces.len(),
                b.faces.len()
            )));
        }
        let conductance_density = 1.0 / resistance.value_m2_k_per_w();
        let mut points = Vec::with_capacity(3 * a.faces.len());
        for (face, &(_, area)) in a.faces.iter().zip(&self.patches[port].faces) {
            for q in 0..3 {
                let weights = degree_two_weights(q);
                points.push(ContactPoint {
                    a: [0, 1, 2].map(|i| (face[i], weights[i])),
                    b: [0, 1, 2].map(|i| (partner[face[i]], weights[i])),
                    conductance_w_k: conductance_density * (area / 3.0),
                });
            }
        }
        Ok(TraceContact {
            width_a: n_a,
            width_b: n_b,
            points,
        })
    }

    /// The owner's contact between trace port `port` of this body (side A)
    /// and trace port `other_port` of `other` (side B) where the traces need
    /// not match: the quadrature points of their common refinement, through
    /// the contact resistance `resistance`.
    ///
    /// `offset_m` places `other` in this body's frame, as for
    /// [`Self::matching_contact`]. The refinement is the owner's own, the
    /// routine that binds a `NonmatchingSurface` on one mesh, under
    /// `options` and with its refusals: traces that leave one plane, do not
    /// face each other, overlap themselves or do not cover each other, and
    /// a budget that is exceeded.
    ///
    /// # Errors
    /// [`ConductionError::Config`] when a port is a mean port;
    /// [`ConductionError::Interface`] for the owner's refusals;
    /// [`ConductionError::Cancelled`].
    #[allow(clippy::too_many_arguments)]
    pub fn mapped_contact(
        &self,
        cx: &Cx<'_>,
        port: usize,
        other: &ConductionBody<'_>,
        other_port: usize,
        offset_m: [f64; 3],
        resistance: &InterfaceResistance,
        options: NonmatchingOptions,
    ) -> Result<TraceContact, ConductionError> {
        let (Some(a), Some(b)) = (
            self.patches[port].trace.as_ref(),
            other.patches[other_port].trace.as_ref(),
        ) else {
            return Err(invalid(
                "a mapped contact joins two trace ports; a mean port has one temperature".into(),
            ));
        };
        let faces_a = self.trace_faces(port, a, [0.0; 3]);
        let faces_b = other.trace_faces(other_port, b, offset_m);
        let conductance_density = 1.0 / resistance.value_m2_k_per_w();
        let points = trace_overlaps(cx, "trace contact", &faces_a, &faces_b, &options)?
            .into_iter()
            .map(|(a, b, weight_m2)| ContactPoint {
                a,
                b,
                conductance_w_k: conductance_density * weight_m2,
            })
            .collect();
        Ok(TraceContact {
            width_a: a.vertices.len(),
            width_b: b.vertices.len(),
            points,
        })
    }

    /// The faces of trace port `port` over the trace's own points, placed
    /// `offset_m` from where the mesh has them.
    fn trace_faces(&self, port: usize, trace: &Trace, offset_m: [f64; 3]) -> Vec<TraceFace> {
        let positions = self.problem.mesh.positions();
        trace
            .faces
            .iter()
            .zip(&trace.normals)
            .zip(&self.patches[port].faces)
            .map(|((&vertices, &normal), &(mesh_vertices, area))| TraceFace {
                vertices,
                points: mesh_vertices.map(|vertex| {
                    let [x, y, z] = positions[vertex];
                    [x + offset_m[0], y + offset_m[1], z + offset_m[2]]
                }),
                normal,
                area,
            })
            .collect()
    }

    /// The owner's boundary with each port's heat applied on its patch: a
    /// mean port's heat rate as a uniform flux, a trace port's loads as the
    /// nodal flux whose consistent loads they are.
    fn boundary_for(&self, inputs: &[f64]) -> Result<ThermalBoundary, ConductionError> {
        let mut uniform = Vec::with_capacity(self.patches.len());
        let mut nodal = Vec::new();
        for patch in &self.patches {
            match &patch.trace {
                None => uniform.push((patch.region, -inputs[patch.offset] / patch.area_m2)),
                Some(trace) => {
                    let loads = &inputs[patch.offset..patch.offset + trace.vertices.len()];
                    let mut field = vec![0.0; self.weights.len()];
                    for (&vertex, flux) in trace.vertices.iter().zip(trace.outward_flux(loads)) {
                        field[vertex] = flux;
                    }
                    nodal.push((patch.region, field));
                }
            }
        }
        self.problem
            .boundary
            .with_uniform_neumann_additions(&uniform)?
            .with_nodal_neumann_replacements(nodal)
    }

    /// The owner's tangent at `old` under the held port inputs.
    fn linearize<'a>(
        &'a self,
        old: &[f64],
        dt_s: f64,
        inputs: &[f64],
        cx: &Cx<'_>,
    ) -> Result<StepLinearization<'a>, ConductionError> {
        self.check_width("inputs", inputs.len())?;
        let boundary = self.boundary_for(inputs)?;
        let problem = ConductionProblem {
            boundary: &boundary,
            ..self.problem
        };
        self.solver
            .linearize_step(cx, problem, None, old, dt_s, self.config, None, &[])
    }

    /// What every port returns of a nodal field, into `outputs`: a mean
    /// port its patch's area mean, a trace port the values at its points.
    fn port_readings(&self, field: &[f64], outputs: &mut [f64]) {
        for patch in &self.patches {
            match &patch.trace {
                None => {
                    let mut sum = 0.0;
                    for &(vertices, area) in &patch.faces {
                        for vertex in vertices {
                            sum += (area / patch.area_m2 / 3.0) * field[vertex];
                        }
                    }
                    outputs[patch.offset] = sum;
                }
                Some(trace) => {
                    for (output, &vertex) in outputs[patch.offset..].iter_mut().zip(&trace.vertices)
                    {
                        *output = field[vertex];
                    }
                }
            }
        }
    }

    /// The nodal loads [W] of the port inputs `inputs`: a mean port's heat
    /// rate spread uniformly over its patch, `A_f / 3` of each face's flux on
    /// each of its vertices, and a trace port's loads on their own vertices.
    fn nodal_loads(&self, inputs: &[f64], loads: &mut [f64]) {
        for patch in &self.patches {
            match &patch.trace {
                None => {
                    let q = inputs[patch.offset];
                    for &(vertices, area) in &patch.faces {
                        for vertex in vertices {
                            loads[vertex] += (area / 3.0) * (q / patch.area_m2);
                        }
                    }
                }
                Some(trace) => {
                    for (&vertex, &load) in trace.vertices.iter().zip(&inputs[patch.offset..]) {
                        loads[vertex] += load;
                    }
                }
            }
        }
    }

    fn check_width(&self, vector: &str, length: usize) -> Result<(), ConductionError> {
        if length == self.width {
            Ok(())
        } else {
            Err(invalid(format!(
                "{vector} have {length} entries for the {} of {} ports",
                self.width,
                self.ports.len()
            )))
        }
    }
}

/// The port `port` declares on its region of `problem`'s boundary, and its
/// patch, whose first entry in the body's inputs and outputs is `offset`.
fn patch_of(
    problem: &ConductionProblem<'_>,
    port: PatchPort<'_>,
    offset: usize,
) -> Result<(PortSpec, Patch), ConductionError> {
    let name = port.region();
    let boundary = problem.boundary;
    let region = boundary
        .region_names()
        .iter()
        .position(|region| region == name)
        .ok_or_else(|| invalid(format!("port {name} names no boundary region")))?;
    let faces: Vec<([usize; 3], f64)> = problem
        .mesh
        .boundary()
        .iter()
        .enumerate()
        .filter(|&(slot, _)| boundary.region_for(slot) == Some(region))
        .map(|(_, face)| (face.vertices.map(|v| v as usize), face.area))
        .collect();
    let area_m2: f64 = faces.iter().map(|&(_, area)| area).sum();
    if !(area_m2.is_finite() && area_m2 > 0.0) {
        return Err(invalid(format!("port {name} has no boundary area")));
    }
    let id = StableId::new(name).map_err(|error| invalid(format!("port {name}: {error:?}")))?;
    let trace = match port {
        PatchPort::Mean(_) => None,
        PatchPort::Trace(_) => {
            let normals = problem
                .mesh
                .boundary()
                .iter()
                .enumerate()
                .filter(|&(slot, _)| boundary.region_for(slot) == Some(region))
                .map(|(_, face)| face.outward_normal)
                .collect();
            Some(trace_of(name, &faces, normals)?)
        }
    };
    let spec = match &trace {
        None => PortSpec::heat(id, PortRole::FlowIn),
        Some(trace) => {
            let points = NonZeroUsize::new(trace.vertices.len())
                .ok_or_else(|| invalid(format!("port {name} has no vertex")))?;
            PortSpec::heat_field(id, PortRole::FlowIn, points)
        }
    }
    .map_err(|error| invalid(format!("port {name}: {error:?}")))?;
    let patch = Patch {
        region,
        area_m2,
        faces,
        offset,
        trace,
    };
    Ok((spec, patch))
}

/// The points of a patch and the Cholesky factor of its consistent P1 face
/// mass, `A_f (1 + δ_ab) / 12` over each face, the weights of the owner's
/// Neumann rows.
fn trace_of(
    name: &str,
    faces: &[([usize; 3], f64)],
    normals: Vec<[f64; 3]>,
) -> Result<Trace, ConductionError> {
    let vertices: Vec<usize> = faces
        .iter()
        .flat_map(|&(vertices, _)| vertices)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let n = vertices.len();
    let point = |vertex: usize| {
        vertices
            .binary_search(&vertex)
            .expect("a face vertex is a point of its patch")
    };
    let local: Vec<[usize; 3]> = faces
        .iter()
        .map(|&(vertices, _)| vertices.map(point))
        .collect();
    let mut factor = vec![0.0; n * n];
    for (face, &(_, area)) in local.iter().zip(faces) {
        for (a, &row) in face.iter().enumerate() {
            for (b, &column) in face.iter().enumerate() {
                factor[row * n + column] += if a == b { area / 6.0 } else { area / 12.0 };
            }
        }
    }
    for i in 0..n {
        for j in 0..=i {
            let mut sum = factor[i * n + j];
            for k in 0..j {
                sum -= factor[i * n + k] * factor[j * n + k];
            }
            if i == j {
                if !(sum.is_finite() && sum > 0.0) {
                    return Err(invalid(format!(
                        "port {name}: the patch's face mass is not positive definite at point {i}"
                    )));
                }
                factor[i * n + i] = sum.sqrt();
            } else {
                factor[i * n + j] = sum / factor[j * n + j];
            }
        }
        for j in (i + 1)..n {
            factor[i * n + j] = 0.0;
        }
    }
    Ok(Trace {
        vertices,
        faces: local,
        normals,
        mass_factor: factor,
    })
}

fn insulated(condition: &ThermalBc) -> bool {
    matches!(condition, ThermalBc::Neumann { outward_flux: ScalarField::Uniform(flux) } if *flux == 0.0)
}

fn invalid(what: String) -> ConductionError {
    ConductionError::Config {
        parameter: "conduction body",
        what,
    }
}

impl Component for ConductionBody<'_> {
    type State = Vec<f64>;
    type Error = ConductionError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    /// Affine in the port inputs with slopes that depend on the step alone,
    /// but only to the tolerance of the owner's Krylov solve.
    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, state: &Vec<f64>) -> Result<f64, ConductionError> {
        if state.len() != self.weights.len() {
            return Err(ConductionError::FieldLength {
                field: "body temperature",
                expected: self.weights.len(),
                found: state.len(),
            });
        }
        let energy: f64 = self.weights.iter().zip(state).map(|(w, t)| w * t).sum();
        if energy.is_finite() {
            Ok(energy)
        } else {
            Err(ConductionError::NonFinite {
                field: "body stored energy",
                bits: energy.to_bits(),
            })
        }
    }

    fn propose(
        &self,
        old: &Vec<f64>,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Vec<f64>,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ConductionError> {
        self.check_width("inputs", inputs.len())?;
        self.check_width("outputs", outputs.len())?;
        let boundary = self.boundary_for(inputs)?;
        let problem = ConductionProblem {
            boundary: &boundary,
            ..self.problem
        };
        let step = self
            .solver
            .advance(cx, problem, None, old, dt_s, self.config)?;
        self.port_readings(&step.temperature, outputs);
        // The heat in, summed as the kernel sums the port energies: a trace
        // port's loads first, then the ports.
        let mut heat_j = 0.0;
        for patch in &self.patches {
            heat_j += match &patch.trace {
                None => dt_s * inputs[patch.offset],
                Some(trace) => {
                    let loads = &inputs[patch.offset..patch.offset + trace.vertices.len()];
                    dt_s * loads.iter().sum::<f64>()
                }
            };
        }
        let stored_old_j = self.stored_energy_j(old)?;
        let stored_new_j = self.stored_energy_j(&step.temperature)?;
        // Each stored energy is a sum of `n` positive terms, and the owner's
        // change one of `n` terms no larger: `(n + 16) eps` of both ends
        // covers the rounding of all three. What the solve for a trace's
        // flux rounds away of its loads is in the numerical loss with the
        // owner's own residual.
        let units = (self.weights.len() + 16) as f64 * f64::EPSILON;
        let receipt = Receipt {
            stored_change_j: step.stored_energy_change_j,
            dissipated_j: 0.0,
            numerical_loss_j: heat_j - step.stored_energy_change_j,
            entropy_generated_j_per_k: 0.0,
            energy_tolerance_j: units * (stored_old_j.abs() + stored_new_j.abs() + heat_j.abs()),
            residual: step.relative_residual,
            tolerance: self.config.linear.tolerance,
        };
        new.clear();
        new.extend_from_slice(&step.temperature);
        Ok(Proposal::Step(receipt))
    }
}

impl PortTangent for ConductionBody<'_> {
    fn jvp(
        &self,
        old: &Vec<f64>,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ConductionError> {
        self.check_width("tangent inputs", d_inputs.len())?;
        self.check_width("tangent outputs", d_outputs.len())?;
        let step = self.linearize(old, dt_s, inputs, cx)?;
        let mut direction = step.zero_direction();
        self.nodal_loads(d_inputs, &mut direction.nodal_load_w);
        let change = step.apply(cx, &direction)?;
        self.port_readings(&change.temperature_k, d_outputs);
        Ok(())
    }

    fn vjp(
        &self,
        old: &Vec<f64>,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ConductionError> {
        self.check_width("adjoint outputs", bar_outputs.len())?;
        self.check_width("adjoint inputs", bar_inputs.len())?;
        let step = self.linearize(old, dt_s, inputs, cx)?;
        // The readings are the transpose of the nodal loads: a mean port's
        // up to its patch's area, weights `A_f / (3 A_p)` on its vertices,
        // a trace port's exactly, one on each of its points.
        let mut weights = vec![0.0; self.weights.len()];
        for patch in &self.patches {
            match &patch.trace {
                None => {
                    let bar = bar_outputs[patch.offset];
                    for &(vertices, area) in &patch.faces {
                        for vertex in vertices {
                            weights[vertex] += (area / patch.area_m2 / 3.0) * bar;
                        }
                    }
                }
                Some(trace) => {
                    for (&vertex, &bar) in trace.vertices.iter().zip(&bar_outputs[patch.offset..]) {
                        weights[vertex] += bar;
                    }
                }
            }
        }
        let gradient = step.pullback(cx, &weights, &[], &[])?;
        for patch in &self.patches {
            match &patch.trace {
                None => {
                    let mut sum = 0.0;
                    for &(vertices, area) in &patch.faces {
                        for vertex in vertices {
                            sum += (area / 3.0) * (gradient.nodal_load[vertex] / patch.area_m2);
                        }
                    }
                    bar_inputs[patch.offset] = sum;
                }
                Some(trace) => {
                    for (bar, &vertex) in bar_inputs[patch.offset..].iter_mut().zip(&trace.vertices)
                    {
                        *bar = gradient.nodal_load[vertex];
                    }
                }
            }
        }
        Ok(())
    }
}
