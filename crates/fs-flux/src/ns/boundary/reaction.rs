//! Discrete boundary reactions from the original mixed momentum equations.

use super::{
    AdmittedBoundary, OpenBoundaryPolicy, Samples, StokesBoundaryError, StokesPressureReference,
};
use super::transport::{BoundaryTransport, Sums, TransportAccumulator, TransportQuadrature};
use crate::bdm::{edge_gauss_pub, eval_basis};
use crate::ns::{
    FluxSolution, FluxSystem, convection_boundary_trace, convection_edge_quadrature, viscous_normal,
};
use fs_sparse::{Coo, Csr};

/// What the flow boundary actually prescribes; velocity does not imply a wall.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoundaryReactionOwner {
    /// Strong normal velocity and weak SIP velocity, possibly an inlet/outlet.
    VelocityConstraint,
    /// Prescribed total stress, with an optional convective open-boundary law.
    OpenTraction {
        /// None for Stokes; otherwise the explicitly admitted backflow policy.
        backflow: Option<OpenBoundaryPolicy>,
    },
    /// Literal supplied stress with a separately prescribed incoming velocity.
    /// Advective momentum uses that reservoir trace only on inflow; it does
    /// not constrain volume flow or add a mechanical backflow resistance.
    OpenReservoir,
}

/// Integrated two-dimensional effort, divided by the constant density.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BoundaryWrench {
    /// Force per unit out-of-plane depth / density [m³/s²].
    pub force: [f64; 2],
    /// Moment per unit depth / density [m⁴/s²], about the declared origin.
    pub moment: f64,
}

impl BoundaryWrench {
    fn add(&mut self, other: Self) {
        for axis in 0..2 {
            self.force[axis] += other.force[axis];
        }
        self.moment += other.moment;
    }
    fn sample(&mut self, point: [f64; 2], force: [f64; 2]) {
        self.add(Self {
            force,
            moment: cross(point, force),
        });
    }
    fn about(mut self, origin: [f64; 2]) -> Self {
        self.moment -= cross(origin, self.force);
        self
    }
    fn finite(self) -> bool {
        self.force.iter().all(|x| x.is_finite()) && self.moment.is_finite()
    }
    /// Work rate per unit depth / density [m⁴/s³] on an explicitly supplied
    /// rigid virtual motion. Translation is the velocity of the moment origin
    /// [m/s]; angular velocity [1/s] is positive counterclockwise. This does not
    /// advance body/mesh motion or include advective kinetic-energy transport.
    ///
    /// # Errors
    /// Refuses nonfinite motion or unrepresentable work.
    pub fn rigid_power(
        self,
        translation: [f64; 2],
        angular_velocity: f64,
    ) -> Result<f64, StokesBoundaryError> {
        let power = self.force[0] * translation[0]
            + self.force[1] * translation[1]
            + self.moment * angular_velocity;
        if !translation.iter().all(|v| v.is_finite())
            || !angular_velocity.is_finite()
            || !power.is_finite()
        {
            return Err(invalid("nonfinite rigid boundary work"));
        }
        Ok(power)
    }
}

/// Selected-edge exchange at an explicit pressure reference and moment origin.
/// Every effort is divided by the same constant density. Stress and advective
/// transport have distinct meanings even though their integrated units agree.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BoundaryLoad {
    /// Origin [m] for every reported moment.
    pub origin_m: [f64; 2],
    /// Total stress exerted BY the fluid ON the exterior, including recovered
    /// velocity-constraint reaction or the declared open-boundary resistance.
    pub stress: BoundaryWrench,
    /// Supplied open traction, with the opposite sign (fluid on exterior).
    /// Velocity constraints contribute zero to this separately retained term.
    pub imposed: BoundaryWrench,
    /// Added open-boundary backflow stress, also fluid on exterior.
    pub backflow: BoundaryWrench,
    /// Outward conservative momentum transport and its moment, using the
    /// original upwind trace. This is NOT a mechanical wall force. It is zero
    /// when convection is absent from the selected flow model.
    pub outward_momentum: BoundaryWrench,
    /// Actual reference for each selected component: None retains its traction
    /// datum; Some(offset) adds that explicit [m²/s²] offset to its mean-zero
    /// pressure. Unselected components are omitted.
    pub pressure_offsets: Vec<(usize, Option<f64>)>,
}

impl BoundaryLoad {
    fn add(&mut self, other: Self) {
        self.stress.add(other.stress);
        self.imposed.add(other.imposed);
        self.backflow.add(other.backflow);
        self.outward_momentum.add(other.outward_momentum);
    }
    fn finite(&self) -> bool {
        self.origin_m.iter().all(|x| x.is_finite())
            && self.stress.finite()
            && self.imposed.finite()
            && self.backflow.finite()
            && self.outward_momentum.finite()
    }
}

/// One discrete reaction functional, NOT a pointwise Cauchy-traction estimate.
/// Normal constraint equations and the FULL numerical viscous flux restore
/// the original rigid virtual work, including symmetric Nitsche consistency.
#[derive(Clone, Debug)]
pub struct DiscreteBoundaryReaction {
    /// Original mesh boundary edge.
    pub edge: usize,
    /// Independently connected fluid component.
    pub component: usize,
    /// Actual boundary owner; a prescribed velocity may carry normal throughflow.
    pub owner: BoundaryReactionOwner,
    /// The SAME pressure reference as the returned cell pressure.
    pub pressure_reference: StokesPressureReference,
    load: BoundaryLoad,
    normal_measure: [f64; 2],
    normal_moment_measure: f64,
    transport: TransportQuadrature,
}

impl DiscreteBoundaryReaction {
    /// Advected volume, momentum and kinetic power from this edge's original
    /// quadrature, with a uniform carried-velocity offset. Pressure work and
    /// mechanical stress are separate; no pressure datum is needed.
    ///
    /// # Errors
    /// Refuses nonfinite origins/offsets or unrepresentable transport.
    pub fn transport(
        &self,
        origin_m: [f64; 2],
        velocity_offset_m_s: [f64; 2],
    ) -> Result<BoundaryTransport, StokesBoundaryError> {
        self.transport.evaluate(origin_m, velocity_offset_m_s)
    }

    /// Evaluate this edge's functional about origin_m. A MeanZero component
    /// REQUIRES Some(offset), the authored pressure offset [m²/s²] to add to
    /// its reported mean-zero field; Some(0) is an explicit gauge choice.
    /// TractionDatum requires None and keeps its actual imposed stress.
    ///
    /// # Errors
    /// Refuses missing/incompatible pressure references or nonfinite outputs.
    pub fn load(
        &self,
        origin_m: [f64; 2],
        pressure_offset: Option<f64>,
    ) -> Result<BoundaryLoad, StokesBoundaryError> {
        let mut load = self.load_unrecorded(origin_m, pressure_offset)?;
        load.pressure_offsets
            .push((self.component, pressure_offset));
        Ok(load)
    }

    fn load_unrecorded(
        &self,
        origin_m: [f64; 2],
        pressure_offset: Option<f64>,
    ) -> Result<BoundaryLoad, StokesBoundaryError> {
        let offset = match (self.pressure_reference, pressure_offset) {
            (StokesPressureReference::MeanZero, Some(value)) if value.is_finite() => value,
            (StokesPressureReference::TractionDatum, None) => 0.0,
            _ => {
                return Err(invalid(
                    "explicit mean-zero offset or retained traction datum required",
                ));
            }
        };
        let mut load = self.load.clone();
        load.origin_m = origin_m;
        load.stress.force[0] += offset * self.normal_measure[0];
        load.stress.force[1] += offset * self.normal_measure[1];
        load.stress.moment += offset * self.normal_moment_measure;
        load.stress = load.stress.about(origin_m);
        load.imposed = load.imposed.about(origin_m);
        load.backflow = load.backflow.about(origin_m);
        load.outward_momentum = load.outward_momentum.about(origin_m);
        if !load.finite() {
            return Err(invalid("unrepresentable boundary resultant"));
        }
        Ok(load)
    }
}

/// Accepted endpoint reactions, retaining boundary ownership and pressure axes.
#[derive(Clone, Debug)]
pub struct BoundaryReactions {
    edges: Vec<DiscreteBoundaryReaction>,
    components: usize,
}

impl BoundaryReactions {
    /// Original edge order. Evaluation requires each component's pressure datum.
    #[must_use]
    pub fn edges(&self) -> &[DiscreteBoundaryReaction] {
        &self.edges
    }

    /// Signed outward volume flux [m²/s] and its first spatial moment [m³/s]
    /// about origin_m, from the SAME retained advector quadrature as the
    /// reported outward momentum. Zero when convection is absent.
    ///
    /// This is deliberately not reconstructed from a later endpoint velocity:
    /// a retained Oseen/Picard advector need not equal that velocity exactly.
    ///
    /// # Errors
    /// Refuses empty, duplicate or nonboundary selections and nonfinite output.
    pub fn volume_flux_moments(
        &self,
        edges: &[usize],
        origin_m: [f64; 2],
    ) -> Result<(f64, [f64; 2]), StokesBoundaryError> {
        if !origin_m.iter().all(|x| x.is_finite()) {
            return Err(invalid("finite flux origin required"));
        }
        let mut sums = Sums::<3>::new();
        for index in self.selected_transport_edges(edges)? {
            let (flux, moment) = self.edges[index].transport.volume_flux_moments(origin_m);
            sums.add([flux, moment[0], moment[1]]);
        }
        let values = sums.values();
        if !values.iter().all(|v| v.is_finite()) {
            return Err(invalid("unrepresentable boundary volume flux moments"));
        }
        Ok((values[0], [values[1], values[2]]))
    }

    /// Selected-edge volume, linear/angular momentum and kinetic transport,
    /// retaining the original upwind traces and advector. The velocity offset
    /// changes carried velocities while preserving flow relative to the control
    /// surface, as in translation from an ALE frame to the world frame.
    ///
    /// Boundary upwind dissipation is reported separately. Pressure/stress
    /// work, internal energy, heat and species are not included. This is not
    /// a full discrete kinetic-energy balance or passivity guarantee.
    ///
    /// # Errors
    /// Refuses empty/duplicate/nonboundary selections, nonfinite frame values,
    /// and unrepresentable output; never resamples boundary or material data.
    pub fn transport(
        &self,
        edges: &[usize],
        origin_m: [f64; 2],
        velocity_offset_m_s: [f64; 2],
    ) -> Result<BoundaryTransport, StokesBoundaryError> {
        if !origin_m.iter().chain(&velocity_offset_m_s).all(|v| v.is_finite()) {
            return Err(invalid("finite transport origin and velocity offset required"));
        }
        let mut sums = Sums::<6>::new();
        for index in self.selected_transport_edges(edges)? {
            let item = self.edges[index].transport(origin_m, velocity_offset_m_s)?;
            sums.add([
                item.outward_volume_flux_m2_s,
                item.outward_momentum.force[0],
                item.outward_momentum.force[1],
                item.outward_momentum.moment,
                item.outward_kinetic_energy_flux,
                item.upwind_kinetic_dissipation,
            ]);
        }
        let values = sums.values();
        let transport = BoundaryTransport {
            origin_m,
            velocity_offset_m_s,
            outward_volume_flux_m2_s: values[0],
            outward_momentum: BoundaryWrench {
                force: [values[1], values[2]],
                moment: values[3],
            },
            outward_kinetic_energy_flux: values[4],
            upwind_kinetic_dissipation: values[5],
        };
        if !transport.finite() {
            return Err(invalid("unrepresentable boundary transport sum"));
        }
        Ok(transport)
    }

    fn selected_transport_edges(&self, edges: &[usize]) -> Result<Vec<usize>, StokesBoundaryError> {
        if edges.is_empty() {
            return Err(invalid("nonempty transport edge selection required"));
        }
        let mut selected = edges.to_vec();
        selected.sort_unstable();
        if selected.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(invalid("duplicate flux edge"));
        }
        for edge in &mut selected {
            *edge = self.edges.binary_search_by_key(edge, |item| item.edge)
                .map_err(|_| invalid("flux selection must contain original boundary edges"))?;
        }
        Ok(selected)
    }

    /// Sum selected edges in canonical order without reassembly or source
    /// queries. pressure_offsets is indexed by the original connected fluid
    /// component. Selected MeanZero components need an explicit Some(offset);
    /// selected TractionDatum components need None. Output stress acts on the
    /// exterior; momentum transported through an open boundary stays separate.
    ///
    /// # Errors
    /// Refuses empty/duplicate/nonboundary selections, mismatched datum axes,
    /// missing references and nonfinite/unrepresentable resultants.
    pub fn resultant(
        &self,
        edges: &[usize],
        origin_m: [f64; 2],
        pressure_offsets: &[Option<f64>],
    ) -> Result<BoundaryLoad, StokesBoundaryError> {
        if edges.is_empty() || pressure_offsets.len() != self.components {
            return Err(invalid(
                "nonempty boundary selection and component pressure offsets required",
            ));
        }
        let mut selected = edges.to_vec();
        selected.sort_unstable();
        if selected.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(invalid("duplicate reaction edge"));
        }
        let mut total = BoundaryLoad {
            origin_m,
            ..BoundaryLoad::default()
        };
        let mut referenced = vec![false; self.components];
        for edge in selected {
            let index = self
                .edges
                .binary_search_by_key(&edge, |item| item.edge)
                .map_err(|_| invalid("reaction selection must contain original boundary edges"))?;
            let item = &self.edges[index];
            total.add(item.load_unrecorded(origin_m, pressure_offsets[item.component])?);
            referenced[item.component] = true;
        }
        total.pressure_offsets = referenced
            .iter()
            .enumerate()
            .filter_map(|(component, &selected)| {
                selected.then_some((component, pressure_offsets[component]))
            })
            .collect();
        if !total.finite() {
            return Err(invalid("unrepresentable boundary sum"));
        }
        Ok(total)
    }
}

struct Face {
    edge: usize,
    samples: Samples,
    policy: Option<OpenBoundaryPolicy>,
    reservoir: Option<[[f64; 2]; 2]>,
    viscosity: [f64; 2],
    penalty: [f64; 2],
    advection: Option<([f64; 2], [f64; 6])>,
}

/// Only rows removed by strong normal velocity are retained, in a rectangular
/// boundary-sized sparse matrix. No second full operator or source call.
pub(in crate::ns) struct ReactionAssembly {
    global_rows: Vec<usize>,
    matrix: Coo,
    rhs: Vec<f64>,
    faces: Vec<Face>,
}

impl ReactionAssembly {
    pub(super) fn new(
        system: &FluxSystem<'_>,
        boundary: &AdmittedBoundary,
        advection: Option<&[f64]>,
    ) -> Self {
        let mut global_rows = Vec::new();
        let faces = boundary
            .samples
            .iter()
            .enumerate()
            .filter_map(|(edge, samples)| {
                let samples = (*samples)?;
                if matches!(samples, Samples::Velocity(_)) {
                    global_rows.extend([2 * edge, 2 * edge + 1]);
                }
                let cell = system.mesh.edges[edge].tris.0;
                Some(Face {
                    edge,
                    samples,
                    policy: boundary.open_policy(edge),
                    reservoir: boundary.reservoir_velocity(edge),
                    viscosity: [0.0; 2],
                    penalty: [0.0; 2],
                    advection: advection.map(|v| {
                        (
                            [v[2 * edge], v[2 * edge + 1]],
                            core::array::from_fn(|i| v[system.dof(cell, i)]),
                        )
                    }),
                })
            })
            .collect();
        Self {
            matrix: Coo::new(global_rows.len(), system.n),
            rhs: vec![0.0; global_rows.len()],
            global_rows,
            faces,
        }
    }
    pub(in crate::ns) fn push(&mut self, row: usize, column: usize, value: f64) {
        let local = self
            .global_rows
            .binary_search(&row)
            .expect("captured normal velocity row");
        self.matrix.push(local, column, value);
    }
    pub(in crate::ns) fn viscosity(&mut self, edge: usize, q: usize, viscosity: f64, penalty: f64) {
        let face = self
            .faces
            .binary_search_by_key(&edge, |face| face.edge)
            .expect("boundary face");
        self.faces[face].viscosity[q] = viscosity;
        self.faces[face].penalty[q] = penalty;
    }
    pub(in crate::ns) fn capture_rhs(&mut self, rhs: &[f64]) {
        for (local, &global) in self.global_rows.iter().enumerate() {
            self.rhs[local] = rhs[global];
        }
    }
    pub(super) fn finish(self) -> ReactionOperator {
        ReactionOperator {
            matrix: self.matrix.assemble(),
            rhs: self.rhs,
            faces: self.faces,
        }
    }
}

pub(super) struct ReactionOperator {
    matrix: Csr,
    rhs: Vec<f64>,
    faces: Vec<Face>,
}

impl ReactionOperator {
    #[allow(clippy::too_many_lines)] // one boundary's constraint, stress and advective exchange
    pub(super) fn evaluate(
        &self,
        system: &FluxSystem<'_>,
        field: &FluxSolution,
        pressure: &[f64],
        has_traction: &[bool],
    ) -> Result<BoundaryReactions, StokesBoundaryError> {
        let mut row = 0;
        let mut edges = Vec::with_capacity(self.faces.len());
        for face in &self.faces {
            let edge = &system.mesh.edges[face.edge];
            let cell = edge.tris.0;
            let component = system.pressure_component[cell];
            let a = system.mesh.verts[edge.verts.0];
            let b = system.mesh.verts[edge.verts.1];
            let mut load = BoundaryLoad::default();
            let owner = match face.samples {
                Samples::Velocity(prescribed) => {
                    let residual: [f64; 2] = core::array::from_fn(|k| {
                        let (columns, values) = self.matrix.row(row + k);
                        columns.iter().zip(values).fold(
                            -self.rhs[row + k],
                            |r, (&column, &value)| {
                                // Match the RETURNED component datum, including mean subtraction.
                                r + value
                                    * if column < system.n_u {
                                        field.x[column]
                                    } else {
                                        pressure[column - system.n_u]
                                    }
                            },
                        )
                    });
                    row += 2;
                    // m0(v)=mean(v.n), m1(v)=mean(s*v.n), s in [-1/2,1/2].
                    // For rigid rotation the second moment is (b-a) cross n / 12.
                    load.stress.force =
                        [-residual[0] * edge.normal[0], -residual[0] * edge.normal[1]];
                    load.stress.moment = -residual[0] * cross(edge.mid, edge.normal)
                        - residual[1] * cross([b[0] - a[0], b[1] - a[1]], edge.normal) / 12.0;
                    let mut gradient = [[0.0; 2]; 2];
                    for i in 0..6 {
                        for (r, values) in gradient.iter_mut().enumerate() {
                            for (c, value) in values.iter_mut().enumerate() {
                                *value +=
                                    field.x[system.dof(cell, i)] * system.bases[cell].grad[i][r][c];
                            }
                        }
                    }
                    let normal = viscous_normal(gradient, edge.normal, true);
                    for (q, (point, weight)) in edge_gauss_pub(a, b).into_iter().enumerate() {
                        let u = velocity(system, cell, point, &field.x);
                        // Restore FULL vector flux: the stored normal row also
                        // included its Nitsche consistency term. Restoring only
                        // the tangent would omit viscous normal stress.
                        let force = core::array::from_fn(|axis| {
                            -weight
                                * (face.viscosity[q] * normal[axis]
                                    - face.penalty[q] * (u[axis] - prescribed[q][axis]))
                        });
                        load.stress.sample(point, force);
                    }
                    BoundaryReactionOwner::VelocityConstraint
                }
                Samples::Traction(prescribed) => {
                    for (q, (point, weight)) in edge_gauss_pub(a, b).into_iter().enumerate() {
                        load.imposed
                            .sample(point, prescribed[q].map(|t| -weight * t));
                    }
                    load.stress = load.imposed;
                    if face.reservoir.is_some() {
                        BoundaryReactionOwner::OpenReservoir
                    } else {
                        BoundaryReactionOwner::OpenTraction {
                            backflow: face.policy,
                        }
                    }
                }
            };
            let mut transport = TransportAccumulator::new();
            if let Some((moments, advect)) = face.advection {
                for (point, weight, fraction) in convection_edge_quadrature(a, b, moments) {
                    let u = velocity(system, cell, point, &field.x);
                    let mut w = [0.0; 2];
                    for (i, coefficient) in advect.iter().enumerate() {
                        let basis = eval_basis(&system.bases[cell], i, point);
                        for axis in 0..2 {
                            w[axis] += coefficient * basis[axis];
                        }
                    }
                    let wn = w[0] * edge.normal[0] + w[1] * edge.normal[1];
                    let upwind = if wn < 0.0 {
                        match (face.reservoir, face.samples) {
                            (Some(prescribed), _) | (None, Samples::Velocity(prescribed)) => {
                                convection_boundary_trace(prescribed, fraction)
                            }
                            _ => u,
                        }
                    } else {
                        u
                    };
                    let volume_flux = weight * wn;
                    transport.sample(point, volume_flux, upwind, u)?;
                    load.outward_momentum
                        .sample(point, upwind.map(|v| weight * wn * v));
                    if face.policy == Some(OpenBoundaryPolicy::DissipativeBackflow) {
                        load.backflow
                            .sample(point, u.map(|v| -weight * wn.min(0.0) * v));
                    }
                }
                load.stress.add(load.backflow);
            }
            let transport = transport.finish()?;
            if !load.finite() {
                return Err(invalid("nonfinite discrete boundary reaction"));
            }
            edges.push(DiscreteBoundaryReaction {
                edge: face.edge,
                component,
                owner,
                pressure_reference: if has_traction[component] {
                    StokesPressureReference::TractionDatum
                } else {
                    StokesPressureReference::MeanZero
                },
                load,
                normal_measure: edge.normal.map(|n| edge.len * n),
                normal_moment_measure: edge.len * cross(edge.mid, edge.normal),
                transport,
            });
        }
        Ok(BoundaryReactions {
            edges,
            components: has_traction.len(),
        })
    }
}

fn velocity(system: &FluxSystem<'_>, cell: usize, point: [f64; 2], x: &[f64]) -> [f64; 2] {
    let mut u = [0.0; 2];
    for i in 0..6 {
        let basis = eval_basis(&system.bases[cell], i, point);
        for axis in 0..2 {
            u[axis] += x[system.dof(cell, i)] * basis[axis];
        }
    }
    u
}
fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
fn invalid(reason: &'static str) -> StokesBoundaryError {
    StokesBoundaryError::InvalidInput(reason)
}
