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
//! conductance. Where the card's resistance follows the temperatures of the
//! two traces, [`TraceResistanceLaw`] is the factor on those conductances
//! ([`resistance_law`]).
//!
//! The step is the owner's unchanged `advance`; the tangent is the owner's
//! `linearize_step` with no Robin port, driven by the same nodal loads. A
//! conductivity that depends on temperature is put behind the same ports by
//! [`ConductionBody::temperature_dependent`] with the owner's Newton policy:
//! the step is then the owner's `advance_nonlinear` and the tangent
//! `linearize_step` of that endpoint. Heat capacity stays independent of
//! temperature, as in the owner. Either way the body also offers the step as
//! the owner's residual in its end field (`StageResidual`), for a monolithic
//! solve.
//!
//! [`EnthalpyBody`] puts the nonlinear enthalpy owner, whose capacity,
//! reference mass and conductivity can all come from one material card,
//! behind the same ports ([`enthalpy`]).
//!
//! The entropy generated inside the body is the conduction term alone,
//! `dt sum_{i<j} -K_ij (T_i - T_j)^2 / (T_i T_j)` over the owner's operator at
//! the endpoint: a lower bound away from steady state, since backward Euler's
//! temporal term is not added, and a negative sum refuses the step. Robin,
//! Dirichlet and contact traces, volumetric sources and temperature-dependent
//! sheet conduction are refused, not wrapped. Traces that
//! are not coplanar, or do not cover each other, are refused.
//!
//! Both bodies have a snapshot codec: the state is written as the vertex
//! count and every nodal value's bits, and the identity a snapshot is bound
//! to is one hash of the mesh, every element's conductivity, the solver's
//! storage, the ports and the step policies. A conductivity sampled from a
//! material card is bound by its receipts, a declared one by its knots. A
//! body with a part the identity cannot bind has no codec: a boundary with
//! sheet conduction, or a solver whose storage is more than its matrix.

use std::collections::{BTreeMap, BTreeSet};
use std::num::NonZeroUsize;

use fs_blake3::DomainHasher;
use fs_exec::Cx;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_ports::StableId;
use fs_ports::component::{
    Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt, StageResidual,
};

use crate::assemble::assemble_operator_scaled_with_interfaces;
use crate::interface::InterfaceResistance;
use crate::interface::nonmatching::{
    NonmatchingOptions, TraceFace, degree_two_weights, trace_overlaps,
};
use crate::transient::backward_euler::{
    BackwardEuler, NonlinearStepConfig, StepConfig, StepLinearization, StepSolution,
};
use crate::transient::enthalpy::entropy::body_generation_j_k;
use crate::{ConductionError, ConductionProblem, ScalarField, ThermalBc, ThermalBoundary};

pub mod enthalpy;
pub use enthalpy::EnthalpyBody;
pub mod resistance_law;
pub use resistance_law::TraceResistanceLaw;

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
    /// The pairs `[a, b]` of trace points that face each other: every point
    /// of a matching contact with its partner, and none of a mapped one.
    ///
    /// A contact whose resistance follows the trace temperatures
    /// ([`TraceResistanceLaw`]) is guarded at these pairs, as the owner's
    /// pointwise contact asks its source at the paired corners of every
    /// face: a corner the card has no answer for is refused even where every
    /// quadrature point lies inside the card.
    pub facing: Vec<[usize; 2]>,
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
    /// The owner's Newton policy, present exactly when the conductivity
    /// depends on temperature: every step and tangent is then the owner's
    /// `advance_nonlinear`.
    nonlinear: Option<NonlinearStepConfig>,
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
    /// temperature (a conductivity that depends on it is put behind ports by
    /// [`Self::temperature_dependent`]). `config` is the owner's step
    /// configuration, used for every step and tangent.
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
        Self::admit(solver, problem, ports, config, None)
    }

    /// Put a body whose conductivity depends on temperature behind heat
    /// ports, as [`Self::new`] puts one whose conductivity does not. Every
    /// step is the owner's `advance_nonlinear` under `nonlinear`, its Newton
    /// policy, and every tangent the owner's `linearize_step` of that same
    /// endpoint. The owner's refusals stand: an endpoint that leaves the
    /// conductivity's temperature span, or a Newton solve that does not meet
    /// the policy, refuses the step; a tangent at a vertex on a knot of the
    /// table refuses the tangent.
    ///
    /// # Errors
    /// Those of [`Self::new`]; a conductivity that does not depend on
    /// temperature (it steps on the owner's linear path, through
    /// [`Self::new`]); a policy the owner refuses.
    pub fn temperature_dependent(
        solver: BackwardEuler<'m>,
        problem: ConductionProblem<'m>,
        ports: &[PatchPort<'_>],
        config: StepConfig,
        nonlinear: NonlinearStepConfig,
    ) -> Result<Self, ConductionError> {
        nonlinear.validate()?;
        Self::admit(solver, problem, ports, config, Some(nonlinear))
    }

    fn admit(
        solver: BackwardEuler<'m>,
        problem: ConductionProblem<'m>,
        ports: &[PatchPort<'_>],
        config: StepConfig,
        nonlinear: Option<NonlinearStepConfig>,
    ) -> Result<Self, ConductionError> {
        if !std::ptr::eq(solver.mesh(), problem.mesh) {
            return Err(invalid(
                "the solver and the problem must use the same mesh".into(),
            ));
        }
        admit_boundary(&problem)?;
        let dependent = match problem.element_materials {
            Some(materials) => (0..problem.mesh.element_count()).try_fold(false, |any, e| {
                Ok::<_, ConductionError>(any || materials.model_for(e)?.is_temperature_dependent())
            })?,
            None => problem.material.is_temperature_dependent(),
        };
        match (dependent, nonlinear.is_some()) {
            (true, false) => {
                return Err(invalid(
                    "a temperature-dependent conductivity steps on the owner's Newton path: \
                     ConductionBody::temperature_dependent with its policy"
                        .into(),
                ));
            }
            (false, true) => {
                return Err(invalid(
                    "a conductivity independent of temperature steps on the owner's linear path: \
                     ConductionBody::new"
                        .into(),
                ));
            }
            _ => {}
        }
        let (specs, patches, width) = admit_ports(&problem, ports)?;
        let ones = vec![1.0; problem.mesh.vertex_count()];
        let mut weights = vec![0.0; ones.len()];
        solver.capacitance().spmv(&ones, &mut weights);
        Ok(Self {
            solver,
            problem,
            config,
            nonlinear,
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
                    "the traces do not match: no point of the second trace lies at {point:?}; `mapped_contact` joins traces that do not match"
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
            facing: partner
                .iter()
                .enumerate()
                .map(|(point, &other)| [point, other])
                .collect(),
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
            // The owner guards the paired corners of exact pairs only.
            facing: Vec::new(),
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
        boundary_for(&self.patches, &self.problem, inputs)
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
        self.solver.linearize_step(
            cx,
            problem,
            None,
            old,
            dt_s,
            self.config,
            self.nonlinear,
            &[],
        )
    }

    /// What every port returns of a nodal field, into `outputs`: a mean
    /// port its patch's area mean, a trace port the values at its points.
    fn port_readings(&self, field: &[f64], outputs: &mut [f64]) {
        port_readings(&self.patches, field, outputs);
    }

    /// The nodal loads [W] of the port inputs `inputs`: a mean port's heat
    /// rate spread uniformly over its patch, `A_f / 3` of each face's flux on
    /// each of its vertices, and a trace port's loads on their own vertices.
    fn nodal_loads(&self, inputs: &[f64], loads: &mut [f64]) {
        nodal_loads(&self.patches, inputs, loads);
    }

    fn check_width(&self, vector: &str, length: usize) -> Result<(), ConductionError> {
        check_width(self.width, self.ports.len(), vector, length)
    }

    /// Whether the identity binds all of this body: its storage is the
    /// capacitance matrix alone and its boundary carries no sheet.
    fn has_codec(&self) -> bool {
        self.solver.is_plain() && unsheeted(&self.problem)
    }
}

/// The checks every body behind heat ports makes of its boundary and source:
/// every region a uniform zero-flux Neumann row (a port's heat replaces that
/// zero; heat through any other region would cross no port), no
/// temperature-dependent sheet conduction, no volumetric source.
fn admit_boundary(problem: &ConductionProblem<'_>) -> Result<(), ConductionError> {
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
    Ok(())
}

/// The ports `ports` declare on `problem`'s boundary, their patches and
/// their entries in all; each region once, and at least one.
fn admit_ports(
    problem: &ConductionProblem<'_>,
    ports: &[PatchPort<'_>],
) -> Result<(Vec<PortSpec>, Vec<Patch>, usize), ConductionError> {
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
        let (spec, patch) = patch_of(problem, port, width)?;
        width += patch.width();
        specs.push(spec);
        patches.push(patch);
    }
    if specs.is_empty() {
        return Err(invalid("a body needs at least one port".into()));
    }
    Ok((specs, patches, width))
}

/// `problem`'s boundary with each port's heat in `inputs` applied on its
/// patch: a mean port's heat rate as a uniform flux, a trace port's loads as
/// the nodal flux whose consistent loads they are.
fn boundary_for(
    patches: &[Patch],
    problem: &ConductionProblem<'_>,
    inputs: &[f64],
) -> Result<ThermalBoundary, ConductionError> {
    let mut uniform = Vec::with_capacity(patches.len());
    let mut nodal = Vec::new();
    for patch in patches {
        match &patch.trace {
            None => uniform.push((patch.region, -inputs[patch.offset] / patch.area_m2)),
            Some(trace) => {
                let loads = &inputs[patch.offset..patch.offset + trace.vertices.len()];
                let mut field = vec![0.0; problem.mesh.vertex_count()];
                for (&vertex, flux) in trace.vertices.iter().zip(trace.outward_flux(loads)) {
                    field[vertex] = flux;
                }
                nodal.push((patch.region, field));
            }
        }
    }
    problem
        .boundary
        .with_uniform_neumann_additions(&uniform)?
        .with_nodal_neumann_replacements(nodal)
}

/// What every port returns of a nodal field, into `outputs`: a mean port its
/// patch's area mean, a trace port the values at its points.
fn port_readings(patches: &[Patch], field: &[f64], outputs: &mut [f64]) {
    for patch in patches {
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
                for (output, &vertex) in outputs[patch.offset..].iter_mut().zip(&trace.vertices) {
                    *output = field[vertex];
                }
            }
        }
    }
}

/// The nodal loads [W] of the port inputs `inputs`, added into `loads`.
fn nodal_loads(patches: &[Patch], inputs: &[f64], loads: &mut [f64]) {
    for patch in patches {
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

/// The transpose of [`port_readings`]: the nodal weights, added into
/// `weights`, whose sum against a field is `bar_outputs` against its
/// readings. A mean port's are `A_f / (3 A_p)` on its vertices, a trace
/// port's one on each of its points.
fn reading_weights(patches: &[Patch], bar_outputs: &[f64], weights: &mut [f64]) {
    for patch in patches {
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
}

/// The transpose of [`nodal_loads`]: each port input's derivative, into
/// `bar_inputs`, from the derivatives `nodal_bar` with respect to the nodal
/// loads.
fn input_bars(patches: &[Patch], nodal_bar: &[f64], bar_inputs: &mut [f64]) {
    for patch in patches {
        match &patch.trace {
            None => {
                let mut sum = 0.0;
                for &(vertices, area) in &patch.faces {
                    for vertex in vertices {
                        sum += (area / 3.0) * (nodal_bar[vertex] / patch.area_m2);
                    }
                }
                bar_inputs[patch.offset] = sum;
            }
            Some(trace) => {
                for (bar, &vertex) in bar_inputs[patch.offset..].iter_mut().zip(&trace.vertices) {
                    *bar = nodal_bar[vertex];
                }
            }
        }
    }
}

/// The heat the port inputs put in over `dt_s` [J], summed as the kernel
/// sums the port energies: a trace port's loads first, then the ports.
fn heat_in_j(patches: &[Patch], dt_s: f64, inputs: &[f64]) -> f64 {
    let mut heat_j = 0.0;
    for patch in patches {
        heat_j += match &patch.trace {
            None => dt_s * inputs[patch.offset],
            Some(trace) => {
                let loads = &inputs[patch.offset..patch.offset + trace.vertices.len()];
                dt_s * loads.iter().sum::<f64>()
            }
        };
    }
    heat_j
}

fn check_width(
    width: usize,
    ports: usize,
    vector: &str,
    length: usize,
) -> Result<(), ConductionError> {
    if length == width {
        Ok(())
    } else {
        Err(invalid(format!(
            "{vector} have {length} entries for the {width} of {ports} ports"
        )))
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

/// A nodal field as the words of a durable snapshot: the vertex count, then
/// every value's bits in vertex order.
fn encode_field(field: &[f64], words: &mut Vec<u64>) {
    words.push(field.len() as u64);
    words.extend(field.iter().map(|value| value.to_bits()));
}

/// The nodal field [`encode_field`] wrote for a mesh of `vertices` vertices.
/// Refused: another count, a value that is not finite, a value `admitted`
/// refuses.
fn decode_field(words: &[u64], vertices: usize, admitted: fn(f64) -> bool) -> Option<Vec<f64>> {
    let (&count, values) = words.split_first()?;
    if count != vertices as u64 || values.len() != vertices {
        return None;
    }
    let field: Vec<f64> = values.iter().map(|&bits| f64::from_bits(bits)).collect();
    field
        .iter()
        .all(|&value| value.is_finite() && admitted(value))
        .then_some(field)
}

/// Bind to `hasher` what a body's problem and ports fix of its physics, for
/// the identity of a durable snapshot: the mesh (every vertex position and
/// every element's vertices), every element's conductivity, and each port's
/// region, faces and kind. The regions' conditions and the source are not
/// bound: a body admits only insulated regions and a zero source. Sheet
/// conduction on the boundary is not bound either, and a body with any has
/// no codec ([`unsheeted`]).
fn bind_problem(hasher: &mut DomainHasher, problem: &ConductionProblem<'_>, patches: &[Patch]) {
    let count = |hasher: &mut DomainHasher, n: usize| hasher.update(&(n as u64).to_le_bytes());
    let mesh = problem.mesh;
    count(hasher, mesh.vertex_count());
    for coordinate in mesh.positions().iter().flatten() {
        hasher.update(&coordinate.to_bits().to_le_bytes());
    }
    let tets = &mesh.complex().tets;
    count(hasher, tets.len());
    for vertex in tets.iter().flatten() {
        hasher.update(&vertex.to_le_bytes());
    }
    problem.material.bind_identity(hasher);
    match problem.element_materials {
        None => hasher.update(&[0]),
        Some(materials) => {
            hasher.update(&[1]);
            materials.bind_identity(hasher);
        }
    }
    count(hasher, patches.len());
    for patch in patches {
        count(hasher, patch.region);
        count(hasher, patch.faces.len());
        for &vertex in patch.faces.iter().flat_map(|(vertices, _)| vertices) {
            count(hasher, vertex);
        }
        hasher.update(&[u8::from(patch.trace.is_some())]);
    }
}

/// Whether `problem`'s boundary carries no sheet conduction. A body admits a
/// sheet that does not depend on temperature, and no identity binds one.
fn unsheeted(problem: &ConductionProblem<'_>) -> bool {
    problem.boundary.surface_heat_conduction().is_empty()
}

/// Bind the solver's capacitance matrix to `hasher`: every row's columns
/// and entries. Its row sums alone would not tell two storages apart.
fn bind_capacitance(hasher: &mut DomainHasher, solver: &BackwardEuler<'_>) {
    let capacity = solver.capacitance();
    hasher.update(&(capacity.nrows() as u64).to_le_bytes());
    for row in 0..capacity.nrows() {
        let (columns, entries) = capacity.row(row);
        hasher.update(&(columns.len() as u64).to_le_bytes());
        for (&column, entry) in columns.iter().zip(entries) {
            hasher.update(&(column as u64).to_le_bytes());
            hasher.update(&entry.to_bits().to_le_bytes());
        }
    }
}

/// Bind numbers a body holds to `hasher`: their count, then their bits.
fn bind_values(hasher: &mut DomainHasher, values: &[f64]) {
    hasher.update(&(values.len() as u64).to_le_bytes());
    for value in values {
        hasher.update(&value.to_bits().to_le_bytes());
    }
}

/// Bind a step policy to `hasher` as its `Debug` text. A derived `Debug`
/// prints every field, each number as the shortest text that reads back to
/// the same bits, so a policy that gains a field is bound with no edit
/// here.
fn bind_policy(hasher: &mut DomainHasher, policy: &dyn core::fmt::Debug) {
    let text = format!("{policy:?}");
    hasher.update(&(text.len() as u64).to_le_bytes());
    hasher.update(text.as_bytes());
}

/// The identity a body appends: its hash as four words.
fn identity_of(hasher: &DomainHasher, words: &mut Vec<u64>) {
    let hash = hasher.finalize();
    for word in hash.as_bytes().as_chunks::<8>().0 {
        words.push(u64::from_le_bytes(*word));
    }
}

impl Component for ConductionBody<'_> {
    type State = Vec<f64>;
    type Error = ConductionError;

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-conduction::component::ConductionBody",
        state_family: "nodal-temperature-field",
        schema_grammar: "the temperature [K] at every vertex of the body's mesh, in vertex order",
        codec_grammar: "the vertex count, then each temperature as the bits of an f64",
        codec_version: 1,
    });

    /// Refuses, writing nothing, for a body with a part its identity cannot
    /// bind: a boundary with sheet conduction, or a solver whose storage is
    /// more than its capacitance matrix.
    fn encode_state(&self, state: &Vec<f64>, words: &mut Vec<u64>) -> bool {
        if !self.has_codec() {
            return false;
        }
        encode_field(state, words);
        true
    }

    /// Refuses a count that is not the mesh's, a temperature that is not
    /// positive and finite, and everything for a body `encode_state`
    /// refuses.
    fn decode_state(&self, words: &[u64]) -> Option<Vec<f64>> {
        if !self.has_codec() {
            return None;
        }
        decode_field(words, self.weights.len(), |kelvin| kelvin > 0.0)
    }

    /// One hash, as four words, of what fixes the body's steps: the mesh,
    /// every element's conductivity and the ports ([`bind_problem`]), the
    /// solver's capacitance matrix, the step configuration and the Newton
    /// policy of a body whose conductivity depends on temperature. The body
    /// steps itself, so a resume under another tolerance could not give the
    /// uninterrupted run's state.
    fn identity_words(&self, words: &mut Vec<u64>) {
        let mut hasher = DomainHasher::new("fs-conduction/conduction-body/v1");
        bind_problem(&mut hasher, &self.problem, &self.patches);
        bind_capacitance(&mut hasher, &self.solver);
        bind_policy(&mut hasher, &self.config);
        bind_policy(&mut hasher, &self.nonlinear);
        identity_of(&hasher, words);
    }

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
        let step = match self.nonlinear {
            None => self
                .solver
                .advance(cx, problem, None, old, dt_s, self.config)?,
            Some(policy) => {
                self.solver
                    .advance_nonlinear(cx, problem, None, old, dt_s, self.config, policy)?
                    .step
            }
        };
        self.port_readings(&step.temperature, outputs);
        let receipt = self.step_receipt(cx, problem, old, dt_s, inputs, &step)?;
        new.clear();
        new.extend_from_slice(&step.temperature);
        Ok(Proposal::Step(receipt))
    }
}

impl ConductionBody<'_> {
    /// The receipt of the owner's step `step` from `old` under the held port
    /// inputs, whichever solver found its endpoint.
    fn step_receipt(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        old: &Vec<f64>,
        dt_s: f64,
        inputs: &[f64],
        step: &StepSolution,
    ) -> Result<Receipt, ConductionError> {
        // The operator the owner's energy equation uses at the endpoint.
        let end = &step.temperature;
        let system = assemble_operator_scaled_with_interfaces(
            cx,
            problem.mesh,
            problem.boundary,
            problem.material,
            problem.source,
            end,
            None,
            None,
            problem.element_materials,
        )?;
        let generated_j_per_k = body_generation_j_k(cx, &system.operator, end, dt_s)?;
        let heat_j = heat_in_j(&self.patches, dt_s, inputs);
        let stored_old_j = self.stored_energy_j(old)?;
        let stored_new_j = self.stored_energy_j(&step.temperature)?;
        // Each stored energy is a sum of `n` positive terms, and the owner's
        // change one of `n` terms no larger: `(n + 16) eps` of both ends
        // covers the rounding of all three. What the solve for a trace's
        // flux rounds away of its loads is in the numerical loss with the
        // owner's own residual.
        let units = (self.weights.len() + 16) as f64 * f64::EPSILON;
        Ok(Receipt {
            stored_change_j: step.stored_energy_change_j,
            dissipated_j: 0.0,
            numerical_loss_j: heat_j - step.stored_energy_change_j,
            entropy_generated_j_per_k: generated_j_per_k,
            energy_tolerance_j: units * (stored_old_j.abs() + stored_new_j.abs() + heat_j.abs()),
            residual: step.relative_residual,
            tolerance: self.config.linear.tolerance,
        })
    }
}

/// The step as the owner's residual in the nodal temperatures `z`: the
/// unknowns are the end field, the residual `C (z - T_old) + dt (K(z) z - b)`
/// [J] that `propose` solves, its derivative `C + dt J(z)` the owner's Newton
/// matrix. `accept` books the owner's accounts on the supplied field; the
/// owner solves nothing there, so its receipt's residual is zero, as every
/// closed-form receipt's.
impl StageResidual for ConductionBody<'_> {
    fn unknowns(&self) -> usize {
        self.weights.len()
    }

    fn initial_guess(
        &self,
        old: &Vec<f64>,
        _dt_s: f64,
        _inputs: &[f64],
        z: &mut [f64],
    ) -> Result<(), ConductionError> {
        if old.len() != z.len() || z.len() != self.weights.len() {
            return Err(invalid("a stage guess needs one entry per vertex".into()));
        }
        z.copy_from_slice(old);
        Ok(())
    }

    fn residual(
        &self,
        old: &Vec<f64>,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        r: &mut [f64],
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ConductionError> {
        self.check_width("inputs", inputs.len())?;
        self.check_width("outputs", outputs.len())?;
        let boundary = self.boundary_for(inputs)?;
        let problem = ConductionProblem {
            boundary: &boundary,
            ..self.problem
        };
        let residual = self
            .solver
            .stage_residual(cx, problem, old, z, dt_s, self.config)?;
        if r.len() != residual.len() {
            return Err(invalid("a residual needs one entry per vertex".into()));
        }
        r.copy_from_slice(&residual);
        self.port_readings(z, outputs);
        Ok(())
    }

    fn jvp_z(
        &self,
        _old: &Vec<f64>,
        dt_s: f64,
        _inputs: &[f64],
        z: &[f64],
        dz: &[f64],
        dr: &mut [f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ConductionError> {
        self.check_width("tangent outputs", d_outputs.len())?;
        // The port loads enter `b` alone: the derivative holds no input.
        let solver = &self.solver;
        solver.stage_jacobian_apply(cx, self.problem, z, dt_s, false, dz, dr)?;
        self.port_readings(dz, d_outputs);
        Ok(())
    }

    fn vjp_z(
        &self,
        _old: &Vec<f64>,
        dt_s: f64,
        _inputs: &[f64],
        z: &[f64],
        bar_r: &[f64],
        bar_outputs: &[f64],
        bar_z: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ConductionError> {
        self.check_width("adjoint outputs", bar_outputs.len())?;
        let solver = &self.solver;
        solver.stage_jacobian_apply(cx, self.problem, z, dt_s, true, bar_r, bar_z)?;
        reading_weights(&self.patches, bar_outputs, bar_z);
        Ok(())
    }

    fn accept(
        &self,
        old: &Vec<f64>,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        new: &mut Vec<f64>,
        cx: &Cx<'_>,
    ) -> Result<Receipt, ConductionError> {
        self.check_width("inputs", inputs.len())?;
        let boundary = self.boundary_for(inputs)?;
        let problem = ConductionProblem {
            boundary: &boundary,
            ..self.problem
        };
        let step = self
            .solver
            .stage_endpoint(cx, problem, old, z, dt_s, self.config)?;
        let receipt = self.step_receipt(cx, problem, old, dt_s, inputs, &step)?;
        new.clear();
        new.extend_from_slice(&step.temperature);
        Ok(receipt)
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
        reading_weights(&self.patches, bar_outputs, &mut weights);
        let gradient = step.pullback(cx, &weights, &[], &[])?;
        input_bars(&self.patches, &gradient.nodal_load, bar_inputs);
        Ok(())
    }
}
