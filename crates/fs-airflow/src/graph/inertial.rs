//! Inertial passive graphs on the original port-Hamiltonian descriptor owner.
//!
//! Every branch retains impulse L Q and kinetic energy L Q^2/2. Free-node
//! pressures are zero-storage multipliers: the skew incidence stamps enforce
//! continuity without a pressure penalty, fictitious capacitance or steady root.
//! Coefficients may be re-resolved at a thermal trial; all accepted flows and
//! work come from the same native midpoint step. This is an incompressible,
//! fixed-reference-density approximation, not acoustic or compressible storage.

use std::sync::Arc;

use fs_exec::Cx;
use fs_phs::{DescriptorPortHamiltonian, PhsError, QuadraticStorage, StepRecord};
use fs_qty::Density;

use super::{FixedPressure, GraphError, LossGraph, Sum, checkpoint};

/// Vector reaction loads on the original fixed pipe geometry.
pub mod loads;
/// Explicit total-pressure junctions, with static traction and kinetic streams.
pub mod total_pressure;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PressureChart {
    Static,
    Total,
}

/// The existing descriptor owner uses bounded dense matrices and a native
/// analytic quadratic tangent. Large spatial CFD belongs to its field owner.
pub const MAX_INERTIAL_GRAPH_STATES: usize = 256;

/// Authored branch inertia geometry; loss geometry remains with the original
/// source-bearing graph. Curved paths use their actual fluid path length.
#[derive(Clone, Copy, Debug)]
pub struct InertialBranch {
    pub length_m: f64,
    pub flow_area_m2: f64,
}

/// Physical approximation caps, separate from nonlinear acceptance tolerances.
#[derive(Clone, Copy, Debug)]
pub struct InertialGraphControls {
    /// Per-branch dt R_secant/L in (0,2]. This bounds the constrained damping
    /// spectrum; individual branches can still reverse through graph coupling.
    pub maximum_drag_step: f64,
    /// Loss plus absolute LOCAL anomalous-gravity work at the largest
    /// old/new/average speed, with the same held constitutive state.
    pub maximum_neglected_power_w: f64,
    pub maximum_neglected_energy_j: f64,
}

#[derive(Debug)]
pub enum InertialGraphError {
    Invalid(&'static str),
    Graph(GraphError),
    Momentum(PhsError),
    PipeLoad(loads::PipeLoadError),
    Continuity {
        node: usize,
        residual_m3_s: f64,
        allowance_m3_s: f64,
    },
    DragStep {
        branch: usize,
        value: f64,
        maximum: f64,
    },
    NeglectedWork {
        power_w: f64,
        energy_j: f64,
    },
    MomentumResidual {
        impulse_pa_s: f64,
        energy_j: f64,
    },
}
impl From<GraphError> for InertialGraphError {
    fn from(error: GraphError) -> Self {
        Self::Graph(error)
    }
}
impl From<PhsError> for InertialGraphError {
    fn from(error: PhsError) -> Self {
        Self::Momentum(error)
    }
}
impl From<loads::PipeLoadError> for InertialGraphError {
    fn from(error: loads::PipeLoadError) -> Self {
        Self::PipeLoad(error)
    }
}
impl core::fmt::Display for InertialGraphError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "inertial hydraulic graph: {self:?}")
    }
}
impl std::error::Error for InertialGraphError {}

#[derive(Clone, Debug)]
struct Compiled {
    pressure_chart: PressureChart,
    topology: LossGraph,
    reference_density: Density,
    branches: Vec<InertialBranch>,
    inertances: Vec<f64>,
    free: Vec<Option<usize>>,
    states: usize,
    j: Vec<f64>,
    g: Vec<f64>,
    storage: QuadraticStorage,
}

/// One retained branch-momentum state. No thermal history is owned here.
/// Clone is cheap for the immutable compiled topology; proposals never mutate
/// this state, so a thermal owner can publish both accepted states atomically.
#[derive(Clone, Debug)]
pub struct InertialGraph {
    compiled: Arc<Compiled>,
    impulses_pa_s: Vec<f64>,
    accepted_time_s: f64,
}

/// Original-law diagnostics and actual native step-average transport.
/// Pressures are interval multipliers, not a new steady solution. Original
/// A Q+B Q|Q| loss, inertia and explicit driving heads complete each branch
/// equation. Endpoint flows must not replace average-flow advection.
#[derive(Clone, Debug)]
pub struct InertialGraphReport {
    pub duration_s: f64,
    pub average_volume_flows_m3_s: Vec<f64>,
    pub average_pressures_pa: Vec<f64>,
    pub node_outflows_m3_s: Vec<f64>,
    pub maximum_node_imbalance_m3_s: f64,
    pub boundary_imbalance_m3_s: f64,
    pub pressure_boundaries: Vec<FixedPressure>,
    pub dissipative_pressure_drops_pa: Vec<f64>,
    pub gravity_pressure_heads_pa: Vec<f64>,
    pub driving_pressure_pa: Vec<f64>,
    pub initial_flows_m3_s: Vec<f64>,
    pub final_flows_m3_s: Vec<f64>,
    /// Candidate minus actual native average, for an enclosing nonlinear solve.
    pub average_flow_residual_m3_s: Vec<f64>,
    pub original_impulse_residuals_pa_s: Vec<f64>,
    pub maximum_impulse_residual_pa_s: f64,
    /// Original kinetic-energy change + original nonlinear loss - supplied work.
    pub mechanical_balance_residual_j: f64,
    pub hydraulic_loss_j: f64,
    pub driving_work_j: f64,
    /// Pressure-port work in the original hydrostatically reduced pressure.
    /// Together with anomalous gravity this equals full static-pressure plus
    /// full gravity work, up to the separately retained continuity residual.
    pub reduced_boundary_work_j: f64,
    pub anomalous_gravity_work_j: f64,
    pub speed_envelope_m3_s: Vec<f64>,
    pub neglected_hydraulic_energy_bound_j: f64,
    pub neglected_buoyancy_energy_bound_j: f64,
    pub neglected_energy_bound_j: f64,
    /// Original native proposal, including algebraic multiplier endpoints.
    /// Rechecking a constitutive response leaves this trajectory and its native
    /// solve receipt unchanged; the physical residuals above use the refreshed
    /// law. Only the branch-impulse prefix is retained as physical history.
    pub momentum: StepRecord,
}

/// An uncommitted physical proposal. Access is read-only until the independent
/// original-law impulse/energy gates accept; the input owner stays untouched.
#[derive(Clone)]
pub struct InertialGraphTrial {
    initial: InertialGraph,
    next: InertialGraph,
    report: InertialGraphReport,
    // Retain the actual refreshed material samples, not the initial topology's
    // gravity binding. Vector pipe loads must use the same accepted body force.
    gravity: Option<Arc<super::GraphGravity>>,
}
impl InertialGraphTrial {
    pub fn report(&self) -> &InertialGraphReport {
        &self.report
    }

    /// Audit this SAME proposed trajectory against a freshly resolved law.
    ///
    /// Neither impulse nor native average flow is stepped again. The exact
    /// original impulses remain available for the physical momentum residual;
    /// refreshed losses, drives, boundary pressures and gravity enter the same
    /// original-law audit used by `propose`. Source support and approximation
    /// caps are re-admitted before the independent acceptance gates are used.
    /// The original candidate-to-average residual is retained unchanged.
    pub fn recheck(
        self,
        cx: &Cx<'_>,
        graph: &LossGraph,
        boundaries: &[FixedPressure],
        driving_pressure_pa: &[f64],
        controls: InertialGraphControls,
    ) -> Result<Self, InertialGraphError> {
        let Self {
            initial,
            next,
            report,
            ..
        } = self;
        let (_, reduced_pressures) = initial.admit_step(
            cx,
            graph,
            boundaries,
            driving_pressure_pa,
            &report.average_volume_flows_m3_s,
            next.accepted_time_s,
            controls,
        )?;
        initial.finish(
            cx,
            graph,
            boundaries,
            driving_pressure_pa,
            next.accepted_time_s,
            controls,
            report.momentum,
            report.average_flow_residual_m3_s,
            reduced_pressures,
        )
    }

    pub fn accept(
        self,
        cx: &Cx<'_>,
        impulse_tolerance_pa_s: f64,
        mechanical_tolerance_j: f64,
    ) -> Result<(InertialGraph, InertialGraphReport), InertialGraphError> {
        checkpoint(cx)?;
        if !nonnegative(impulse_tolerance_pa_s) || !nonnegative(mechanical_tolerance_j) {
            return Err(InertialGraphError::Invalid(
                "finite nonnegative physical tolerances required",
            ));
        }
        if self.report.maximum_impulse_residual_pa_s > impulse_tolerance_pa_s
            || self.report.mechanical_balance_residual_j.abs() > mechanical_tolerance_j
        {
            return Err(InertialGraphError::MomentumResidual {
                impulse_pa_s: self.report.maximum_impulse_residual_pa_s,
                energy_j: self.report.mechanical_balance_residual_j,
            });
        }
        Ok((self.next, self.report))
    }
}

impl InertialGraph {
    /// Bind finite positive inertia, original topology and fixed pressure-node
    /// identities. Each component needs a datum; multiple pressure boundaries
    /// describe physical external ports, whose mass/heat streams remain explicit
    /// in the thermal owner. Initial free-node flow must close within roundoff.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        cx: &Cx<'_>,
        graph: &LossGraph,
        reference_density: Density,
        branches: &[InertialBranch],
        boundaries: &[FixedPressure],
        initial_flows_m3_s: &[f64],
        accepted_time_s: f64,
        maximum_states: usize,
    ) -> Result<Self, InertialGraphError> {
        Self::try_new_in_chart(
            cx,
            graph,
            reference_density,
            branches,
            boundaries,
            initial_flows_m3_s,
            accepted_time_s,
            maximum_states,
            PressureChart::Static,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn try_new_in_chart(
        cx: &Cx<'_>,
        graph: &LossGraph,
        reference_density: Density,
        branches: &[InertialBranch],
        boundaries: &[FixedPressure],
        initial_flows_m3_s: &[f64],
        accepted_time_s: f64,
        maximum_states: usize,
        pressure_chart: PressureChart,
    ) -> Result<Self, InertialGraphError> {
        checkpoint(cx)?;
        let edges = graph.branches().len();
        if !(1..=MAX_INERTIAL_GRAPH_STATES).contains(&maximum_states)
            || graph.node_count() > maximum_states
            || edges == 0
            || edges > maximum_states
            || branches.len() != edges
            || initial_flows_m3_s.len() != edges
            || !positive(reference_density.value())
            || !nonnegative(accepted_time_s)
            || initial_flows_m3_s.iter().any(|v| !v.is_finite())
        {
            return Err(InertialGraphError::Invalid(
                "bounded graph, complete inertia and finite initial state required",
            ));
        }
        let fixed = graph.fixed_pressures(boundaries, cx)?;
        let mut free = vec![None; graph.node_count()];
        let mut states = edges;
        for (node, value) in fixed.iter().enumerate() {
            if value.is_none() {
                free[node] = Some(states);
                states += 1;
            }
        }
        if states > maximum_states {
            return Err(InertialGraphError::Invalid(
                "native descriptor state budget exceeded",
            ));
        }
        admit_reference(graph, reference_density.value())?;
        let mut inertances = Vec::with_capacity(edges);
        let mut impulses = Vec::with_capacity(edges);
        let mut q = vec![0.0; edges * edges];
        let mut j = vec![0.0; states * states];
        let mut g = vec![0.0; states * edges];
        for (edge, (geometry, branch)) in branches.iter().zip(graph.branches()).enumerate() {
            checkpoint(cx)?;
            let inertance = reference_density.value() * (geometry.length_m / geometry.flow_area_m2);
            let impulse = inertance * initial_flows_m3_s[edge];
            if !positive(geometry.length_m)
                || !positive(geometry.flow_area_m2)
                || !positive(inertance)
                || !positive(1.0 / inertance)
                || !impulse.is_finite()
                || !(0.5 * impulse * initial_flows_m3_s[edge]).is_finite()
                || (initial_flows_m3_s[edge] != 0.0 && impulse == 0.0)
            {
                return Err(InertialGraphError::Invalid(
                    "finite positive branch inertance required",
                ));
            }
            inertances.push(inertance);
            impulses.push(impulse);
            q[edge * edges + edge] = 1.0 / inertance;
            g[edge * edges + edge] = 1.0;
            for (node, sign) in [(branch.from, 1.0), (branch.to, -1.0)] {
                if let Some(multiplier) = free[node] {
                    j[edge * states + multiplier] = sign;
                    j[multiplier * states + edge] = -sign;
                }
            }
        }
        let compiled = Compiled {
            pressure_chart,
            topology: graph.clone(),
            reference_density,
            branches: branches.to_vec(),
            inertances,
            free,
            states,
            j,
            g,
            storage: QuadraticStorage::new(q, edges)?,
        };
        continuity(cx, &compiled, initial_flows_m3_s)?;
        let owner = Self {
            compiled: Arc::new(compiled),
            impulses_pa_s: impulses,
            accepted_time_s,
        };
        if !owner.kinetic_energy_j().is_finite() {
            return Err(InertialGraphError::Invalid(
                "initial kinetic energy is unrepresentable",
            ));
        }
        checkpoint(cx)?;
        Ok(owner)
    }

    pub fn node_count(&self) -> usize {
        self.compiled.free.len()
    }
    pub fn branch_count(&self) -> usize {
        self.compiled.inertances.len()
    }
    pub const fn accepted_time_s(&self) -> f64 {
        self.accepted_time_s
    }
    pub fn reference_density(&self) -> Density {
        self.compiled.reference_density
    }
    pub fn inertances_pa_s2_m3(&self) -> &[f64] {
        &self.compiled.inertances
    }
    pub fn volume_flows_m3_s(&self) -> Vec<f64> {
        self.impulses_pa_s
            .iter()
            .zip(&self.compiled.inertances)
            .map(|(p, l)| p / l)
            .collect()
    }
    pub fn kinetic_energy_j(&self) -> f64 {
        self.impulses_pa_s
            .iter()
            .zip(&self.compiled.inertances)
            .map(|(p, l)| 0.5 * p * (p / l))
            .sum()
    }

    /// Evaluate one native descriptor step at held trial coefficients. The
    /// caller closes candidate-average flow against actual native average flow,
    /// preferably in the SAME residual as native thermal energy increments.
    ///
    /// Every trial is pure and repeatable. Source/correlation adapters may use
    /// candidate flow to resolve drag, but must retain their original source
    /// identities and applicability checks. This adapter independently checks
    /// the supplied graph law at actual flow, its Darcy binding at the speed
    /// envelope, and original impulse/work residuals before acceptance.
    ///
    /// Native dense stepping is bounded by maximum_states and its existing
    /// nonlinear budget. Checkpoints bracket it and every graph traversal.
    /// No promise of resolving the within-step reversal instant is made.
    #[allow(clippy::too_many_arguments)]
    pub fn propose(
        &self,
        cx: &Cx<'_>,
        graph: &LossGraph,
        boundaries: &[FixedPressure],
        driving_pressure_pa: &[f64],
        candidate_average_flows_m3_s: &[f64],
        target_s: f64,
        controls: InertialGraphControls,
    ) -> Result<InertialGraphTrial, InertialGraphError> {
        let c = &self.compiled;
        let edges = c.inertances.len();
        let (dt, reduced_pressures) = self.admit_step(
            cx,
            graph,
            boundaries,
            driving_pressure_pa,
            candidate_average_flows_m3_s,
            target_s,
            controls,
        )?;
        let mut r = vec![0.0; c.states * c.states];
        let mut forcing = Vec::with_capacity(edges);
        for (edge, branch) in graph.branches().iter().enumerate() {
            checkpoint(cx)?;
            let resistance = secant(graph, edge, candidate_average_flows_m3_s[edge])?;
            admit_drag(c, edge, resistance, dt, controls.maximum_drag_step)?;
            r[edge * c.states + edge] = resistance;
            forcing.push(finite(
                driving_pressure_pa[edge]
                    + graph.driving_pressure(
                        edge,
                        reduced_pressures[branch.from],
                        reduced_pressures[branch.to],
                    )?,
                "held branch forcing",
            )?);
        }
        let system = DescriptorPortHamiltonian::new(
            c.states,
            edges,
            edges,
            c.j.clone(),
            r,
            c.g.clone(),
            Box::new(c.storage.clone()),
        )?;
        // Algebraic pressure has no physical history. Zero initial multiplier
        // is a coordinate choice; native x1/2 recovers its interval effort.
        let mut old = self.impulses_pa_s.clone();
        old.resize(c.states, 0.0);
        checkpoint(cx)?;
        let momentum = fs_phs::step_descriptor(&system, &old, &forcing, dt)?;
        checkpoint(cx)?;
        let flow_residuals = candidate_average_flows_m3_s
            .iter()
            .zip(&momentum.y)
            .map(|(candidate, actual)| finite(candidate - actual, "average-flow closure"))
            .collect::<Result<Vec<_>, _>>()?;
        self.finish(
            cx,
            graph,
            boundaries,
            driving_pressure_pa,
            target_s,
            controls,
            momentum,
            flow_residuals,
            reduced_pressures,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn admit_step(
        &self,
        cx: &Cx<'_>,
        graph: &LossGraph,
        boundaries: &[FixedPressure],
        driving_pressure_pa: &[f64],
        average_flows_m3_s: &[f64],
        target_s: f64,
        controls: InertialGraphControls,
    ) -> Result<(f64, Vec<f64>), InertialGraphError> {
        checkpoint(cx)?;
        let c = &self.compiled;
        let edges = c.inertances.len();
        let dt = target_s - self.accepted_time_s;
        if !positive(dt)
            || !target_s.is_finite()
            || driving_pressure_pa.len() != edges
            || average_flows_m3_s.len() != edges
            || driving_pressure_pa
                .iter()
                .chain(average_flows_m3_s)
                .any(|v| !v.is_finite())
            || !positive(controls.maximum_drag_step)
            || controls.maximum_drag_step > 2.0
            || !nonnegative(controls.maximum_neglected_power_w)
            || !nonnegative(controls.maximum_neglected_energy_j)
        {
            return Err(InertialGraphError::Invalid(
                "finite complete trial, positive time and explicit physical caps required",
            ));
        }
        admit_topology(c, graph)?;
        let fixed = graph.fixed_pressures(boundaries, cx)?;
        if fixed
            .iter()
            .zip(&c.free)
            .any(|(p, index)| p.is_none() != index.is_some())
        {
            return Err(InertialGraphError::Invalid(
                "pressure-boundary node identities changed",
            ));
        }
        let reduced_boundaries = graph.reduce_boundaries(boundaries)?;
        let reduced_boundaries = reduced_boundaries.as_deref().unwrap_or(boundaries);
        let mut reduced_pressures = vec![0.0; graph.node_count()];
        for boundary in reduced_boundaries {
            reduced_pressures[boundary.node] = boundary.pressure.value();
        }
        Ok((dt, reduced_pressures))
    }

    // The only original-law impulse, work and source-envelope audit. A fresh
    // constitutive response reuses the exact native trajectory through here.
    #[allow(clippy::too_many_arguments)]
    fn finish(
        &self,
        cx: &Cx<'_>,
        graph: &LossGraph,
        boundaries: &[FixedPressure],
        driving_pressure_pa: &[f64],
        target_s: f64,
        controls: InertialGraphControls,
        momentum: StepRecord,
        flow_residuals: Vec<f64>,
        mut reduced_pressures: Vec<f64>,
    ) -> Result<InertialGraphTrial, InertialGraphError> {
        checkpoint(cx)?;
        let c = &self.compiled;
        let edges = c.inertances.len();
        let dt = target_s - self.accepted_time_s;
        let average = &momentum.y;
        if average.len() != edges || momentum.x.len() != c.states || flow_residuals.len() != edges {
            return Err(InertialGraphError::Invalid(
                "native descriptor returned mismatched dimensions",
            ));
        }
        let initial = self.volume_flows_m3_s();
        let final_flows: Vec<_> = momentum.x[..edges]
            .iter()
            .zip(&c.inertances)
            .map(|(p, l)| p / l)
            .collect();
        let (balances, maximum_imbalance, boundary_imbalance) = continuity(cx, c, average)?;
        continuity(cx, c, &final_flows)?;
        for (node, index) in c.free.iter().enumerate() {
            if let Some(index) = index {
                reduced_pressures[node] = 0.5 * momentum.x[*index];
            }
        }
        let mut physical = reduced_pressures.clone();
        if let Some(gravity) = graph.gravity() {
            for (node, p) in physical.iter_mut().enumerate() {
                *p = finite(
                    *p + gravity.reference_pressure(node),
                    "static pressure reconstruction",
                )?;
            }
        }
        for boundary in boundaries {
            physical[boundary.node] = boundary.pressure.value();
        }
        let mut pressure_losses = Vec::with_capacity(edges);
        let mut residuals = Vec::with_capacity(edges);
        let mut envelope = Vec::with_capacity(edges);
        let mut loss_power = Sum::default();
        let mut driving_power = Sum::default();
        let mut boundary_power = Sum::default();
        let mut buoyancy_power = Sum::default();
        let mut envelope_loss_power = Sum::default();
        let mut envelope_buoyancy_power = Sum::default();
        for (edge, branch) in graph.branches().iter().enumerate() {
            checkpoint(cx)?;
            let flow = finite(average[edge], "native average flow")?;
            let resistance = secant(graph, edge, flow)?;
            admit_drag(c, edge, resistance, dt, controls.maximum_drag_step)?;
            let loss = finite(resistance * flow, "original branch loss")?;
            let head = graph.outward_head(edge, branch.from);
            let impulse_residual = finite(
                momentum.x[edge]
                    - self.impulses_pa_s[edge]
                    - dt * (graph.driving_pressure(
                        edge,
                        reduced_pressures[branch.from],
                        reduced_pressures[branch.to],
                    )? + driving_pressure_pa[edge]
                        - loss),
                "original branch impulse residual",
            )?;
            residuals.push(impulse_residual);
            let speed = initial[edge]
                .abs()
                .max(final_flows[edge].abs())
                .max(flow.abs());
            envelope.push(speed);
            if let Some(linear) = &graph.linear_losses()[edge] {
                let pressures = [physical[branch.from], physical[branch.to]];
                let bounds = c.pressure_chart.source_pressure_bounds(
                    pressures,
                    [initial[edge], final_flows[edge]],
                    c.reference_density.value(),
                    c.branches[edge].flow_area_m2,
                )?;
                let samples = if c.pressure_chart == PressureChart::Static {
                    1
                } else {
                    2
                };
                for [from, to] in bounds.into_iter().take(samples) {
                    linear.admit_flow(from, to, speed).map_err(|reason| {
                        GraphError::InvalidBranch {
                            branch: edge,
                            reason,
                        }
                    })?;
                }
            }
            loss_power.add(finite(loss * flow, "original dissipation")?)?;
            driving_power.add(finite(driving_pressure_pa[edge] * flow, "driving power")?)?;
            buoyancy_power.add(finite(head * flow, "anomalous gravity power")?)?;
            envelope_loss_power.add(finite(
                secant(graph, edge, speed)? * speed * speed,
                "loss envelope",
            )?)?;
            let absolute_head = graph
                .gravity()
                .map_or(0.0, |gravity| gravity.absolute_departure_heads_pa()[edge]);
            envelope_buoyancy_power.add(finite(absolute_head * speed, "gravity envelope")?)?;
            pressure_losses.push(loss);
        }
        for boundary in boundaries {
            boundary_power.add(finite(
                reduced_pressures[boundary.node] * balances[boundary.node],
                "boundary power",
            )?)?;
        }
        let hydraulic_loss_j = finite(dt * loss_power.value()?, "hydraulic loss energy")?;
        let driving_work_j = finite(dt * driving_power.value()?, "driving work")?;
        let reduced_boundary_work_j = finite(dt * boundary_power.value()?, "boundary work")?;
        let anomalous_gravity_work_j = finite(dt * buoyancy_power.value()?, "gravity work")?;
        let loss_bound_power = envelope_loss_power.value()?;
        let buoyancy_bound_power = envelope_buoyancy_power.value()?;
        let power_bound = finite(
            loss_bound_power + buoyancy_bound_power,
            "omitted power bound",
        )?;
        let neglected_hydraulic_energy_bound_j = finite(dt * loss_bound_power, "loss work bound")?;
        let neglected_buoyancy_energy_bound_j =
            finite(dt * buoyancy_bound_power, "gravity work bound")?;
        let neglected_energy_bound_j = finite(dt * power_bound, "omitted energy bound")?;
        if power_bound > controls.maximum_neglected_power_w
            || neglected_energy_bound_j > controls.maximum_neglected_energy_j
        {
            return Err(InertialGraphError::NeglectedWork {
                power_w: power_bound,
                energy_j: neglected_energy_bound_j,
            });
        }
        let mechanical_balance_residual_j = finite(
            momentum.delta_h + hydraulic_loss_j
                - driving_work_j
                - reduced_boundary_work_j
                - anomalous_gravity_work_j,
            "mechanical energy balance",
        )?;
        let maximum_impulse_residual_pa_s = residuals.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
        let next = Self {
            compiled: self.compiled.clone(),
            impulses_pa_s: momentum.x[..edges].to_vec(),
            accepted_time_s: target_s,
        };
        checkpoint(cx)?;
        Ok(InertialGraphTrial {
            initial: self.clone(),
            next,
            gravity: graph.gravity.clone(),
            report: InertialGraphReport {
                duration_s: dt,
                average_volume_flows_m3_s: average.clone(),
                average_pressures_pa: physical,
                node_outflows_m3_s: balances,
                maximum_node_imbalance_m3_s: maximum_imbalance,
                boundary_imbalance_m3_s: boundary_imbalance,
                pressure_boundaries: boundaries.to_vec(),
                dissipative_pressure_drops_pa: pressure_losses,
                gravity_pressure_heads_pa: graph
                    .gravity()
                    .map_or_else(|| vec![0.0; edges], |g| g.pressure_heads_pa().to_vec()),
                driving_pressure_pa: driving_pressure_pa.to_vec(),
                initial_flows_m3_s: initial,
                final_flows_m3_s: final_flows,
                average_flow_residual_m3_s: flow_residuals,
                original_impulse_residuals_pa_s: residuals,
                maximum_impulse_residual_pa_s,
                mechanical_balance_residual_j,
                hydraulic_loss_j,
                driving_work_j,
                reduced_boundary_work_j,
                anomalous_gravity_work_j,
                speed_envelope_m3_s: envelope,
                neglected_hydraulic_energy_bound_j,
                neglected_buoyancy_energy_bound_j,
                neglected_energy_bound_j,
                momentum,
            },
        })
    }
}

fn admit_topology(c: &Compiled, graph: &LossGraph) -> Result<(), InertialGraphError> {
    if graph.node_count() != c.topology.node_count()
        || graph.branches().len() != c.inertances.len()
        || graph
            .branches()
            .iter()
            .zip(c.topology.branches())
            .any(|(a, b)| a.from != b.from || a.to != b.to || a.loss.name != b.loss.name)
    {
        return Err(InertialGraphError::Invalid(
            "original branch topology or identities changed",
        ));
    }
    admit_reference(graph, c.reference_density.value())?;
    match (graph.gravity(), c.topology.gravity()) {
        (None, None) => {}
        (Some(a), Some(b))
            if a.same_geometry_as(b)
                && a.maximum_relative_density_departure()
                    <= b.maximum_relative_density_departure() => {}
        _ => {
            return Err(InertialGraphError::Invalid(
                "original gravity geometry or density allowance changed",
            ));
        }
    }
    Ok(())
}
fn admit_reference(graph: &LossGraph, rho: f64) -> Result<(), InertialGraphError> {
    if graph
        .gravity()
        .is_some_and(|g| g.reference_transport().density_kg_m3() != rho)
        || graph
            .linear_losses()
            .iter()
            .flatten()
            .any(|l| l.material().is_some_and(|m| m.gas.density != rho))
    {
        return Err(InertialGraphError::Invalid(
            "hydraulic sources and inertia require the same reference density",
        ));
    }
    Ok(())
}
fn secant(graph: &LossGraph, edge: usize, flow: f64) -> Result<f64, InertialGraphError> {
    let value = graph.linear_losses()[edge]
        .as_ref()
        .map_or(0.0, |l| l.resistance.value())
        + graph.branches()[edge].loss.resistance.value() * flow.abs();
    if nonnegative(value) {
        Ok(value)
    } else {
        Err(InertialGraphError::Invalid(
            "original secant resistance is unrepresentable",
        ))
    }
}
fn admit_drag(
    c: &Compiled,
    edge: usize,
    resistance: f64,
    dt: f64,
    maximum: f64,
) -> Result<(), InertialGraphError> {
    let value = dt * (resistance / c.inertances[edge]);
    if !value.is_finite() || value > maximum {
        return Err(InertialGraphError::DragStep {
            branch: edge,
            value,
            maximum,
        });
    }
    Ok(())
}
fn continuity(
    cx: &Cx<'_>,
    c: &Compiled,
    flows: &[f64],
) -> Result<(Vec<f64>, f64, f64), InertialGraphError> {
    let mut balances = vec![Sum::default(); c.free.len()];
    let mut magnitudes = vec![Sum::default(); c.free.len()];
    for (edge, branch) in c.topology.branches().iter().enumerate() {
        checkpoint(cx)?;
        let flow = finite(flows[edge], "branch continuity flow")?;
        for (node, sign) in [(branch.from, 1.0), (branch.to, -1.0)] {
            balances[node].add(sign * flow)?;
            magnitudes[node].add(flow.abs())?;
        }
    }
    let mut values = Vec::with_capacity(balances.len());
    let mut maximum = 0.0_f64;
    let mut boundary = Sum::default();
    for (node, (balance, magnitude)) in balances.into_iter().zip(magnitudes).enumerate() {
        checkpoint(cx)?;
        let value = balance.value()?;
        if c.free[node].is_some() {
            let allowance = 64.0 * f64::EPSILON * magnitude.value()?;
            if value.abs() > allowance {
                return Err(InertialGraphError::Continuity {
                    node,
                    residual_m3_s: value,
                    allowance_m3_s: allowance,
                });
            }
            maximum = maximum.max(value.abs());
        } else {
            boundary.add(value)?;
        }
        values.push(value);
    }
    Ok((values, maximum, boundary.value()?))
}
fn finite(value: f64, what: &'static str) -> Result<f64, InertialGraphError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(InertialGraphError::Invalid(what))
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}

#[cfg(test)]
mod tests;
