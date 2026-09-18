//! Steady heat transport on a solved passive airflow graph.
//!
//! Incoming streams mix by heat-capacity flow at each zero-volume junction.
//! Every outgoing branch receives that mixture, then applies the existing
//! exponential `AirPath` law. A passive nonzero flow strictly descends pressure,
//! so its directed graph is acyclic even when the undirected loss graph loops.
//! Without return links a deterministic topological pass needs no thermal
//! matrix solve. Optional imposed return links close a bounded supply system.
//!
//! Density, specific heat, hydraulics and heat-transfer coefficients are frozen.
//! This is a nominal dry-air network model, not CFD, solved return hydraulics,
//! buoyancy feedback, an uncertainty enclosure, or experimental validation.
//! The node sensible-enthalpy balance follows EnergyPlus Engineering Reference,
//! "AirflowNetwork Model / Node Temperature Calculations"; no EnergyPlus code
//! or runtime dependency is used.

use std::collections::BTreeSet;
use std::fmt;
use std::ops::Range;

use fs_exec::Cx;
use fs_qty::{Density, Temperature, VolumetricFlowRate};

use crate::AirflowError;
use crate::conjugate::{AirMarch, AirPath, AirSegment};
use crate::graph::GraphSolution;

pub mod sensitivity;
pub mod recirculation;
use recirculation::{Feedback, RecirculationReport};

/// One explicit thermal model for each hydraulic branch, in graph order.
#[derive(Debug, Clone, PartialEq)]
pub enum BranchThermalModel {
    /// No heat exchange, including an explicitly unheated leakage bypass.
    Adiabatic,
    /// Segments ordered from the hydraulic edge's declared `from` to `to`.
    /// Construction reverses this order when the solved flow is negative.
    Exchange(Vec<AirSegment>),
}

/// Frozen coherent-SI dry-air properties; no temperature dependence is inferred.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransportAir {
    /// Constant density in kg/m^3.
    pub density: Density,
    /// Constant specific heat in J/(kg K).
    pub specific_heat_j_kg_k: f64,
}

/// Temperature of external air supplied at a pressure-reservoir node.
/// This prescribes the entering stream, not the mixture with incoming branches.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransportInlet {
    /// Index in the hydraulic solution's pressure array.
    pub node: usize,
    /// Absolute temperature in kelvin.
    pub temperature: Temperature,
}

/// Explicit hydraulic-admission and whole-network sensible-heat tolerances.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransportConfig {
    /// Positive absolute free-node conservation floor in m^3/s.
    pub absolute_flow_tolerance: VolumetricFlowRate,
    /// Relative free-node tolerance, in [0, 1), applied to the larger throughflow.
    pub relative_flow_tolerance: f64,
    /// Positive absolute floor on external enthalpy gain minus wall heat, W.
    pub absolute_heat_tolerance_w: f64,
    /// Relative heat tolerance in [0, 1), applied to the larger gross heat rate.
    pub relative_heat_tolerance: f64,
}

/// A transported branch in actual flow direction. Zero-flow adiabatic branches
/// have no temperature or march: stagnant air is not silently assigned a value.
#[derive(Debug, Clone, PartialEq)]
pub struct BranchTransport {
    /// Temperature entering the branch in actual flow direction, K.
    pub inlet_temperature_k: Option<f64>,
    /// Temperature leaving the branch in actual flow direction, K.
    pub outlet_temperature_k: Option<f64>,
    /// Heated segments in actual flow order; absent for adiabatic branches.
    pub march: Option<AirMarch>,
}

/// A full nominal transport pass, not an uncertainty or passivity certificate.
#[derive(Debug, Clone, PartialEq)]
pub struct TransportMarch {
    /// Perfectly mixed node temperatures; stagnant nodes are `None`.
    pub node_temperatures_k: Vec<Option<f64>>,
    /// Branches in the original hydraulic graph's order.
    pub branches: Vec<BranchTransport>,
    /// Effective Robin references in [`TransportNetwork::regions`] order.
    pub reference_temperatures_k: Vec<f64>,
    /// Heat transferred from all walls into the air, W.
    pub wall_heat_rate_w: f64,
    /// External outlet sensible enthalpy minus external inlet enthalpy, W.
    pub external_heat_gain_w: f64,
    /// `external_heat_gain_w - wall_heat_rate_w`, W. Gated without correction.
    pub heat_imbalance_w: f64,
    /// Contribution predicted from admitted free-node hydraulic residuals, W.
    /// Disclosed for diagnosis, NEVER subtracted to pass the heat-balance gate.
    pub hydraulic_energy_defect_w: f64,
    /// Common reference for sensible enthalpy, K; not a boundary condition.
    pub enthalpy_reference_k: f64,
    /// Fresh/return mixtures and the independent outer balance, when requested.
    pub recirculation: Option<RecirculationReport>,
}

/// Input, transport, thermal-balance and cancellation failures.
#[derive(Debug)]
pub enum TransportError {
    /// Invalid shape, quantity, or tolerance.
    InvalidInput(&'static str),
    /// Invalid or ambiguous node boundary.
    Node { node: usize, reason: &'static str },
    /// Invalid branch or an unsupported zero-flow heat exchanger.
    Branch { branch: usize, reason: &'static str },
    /// A free hydraulic node is not adequately balanced for transport.
    FlowImbalance { node: usize, residual_m3_s: f64, tolerance_m3_s: f64 },
    /// A global heat budget was not met; no march is returned.
    HeatImbalance { residual_w: f64, tolerance_w: f64, hydraulic_defect_w: f64 },
    /// A bounded work checkpoint refused more work.
    Interrupted,
    /// Existing channel or conjugate solver refusal, preserved without promotion.
    Airflow(AirflowError),
}

impl fmt::Display for TransportError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "airflow heat transport: {self:?}")
    }
}
impl std::error::Error for TransportError {}
impl From<AirflowError> for TransportError {
    fn from(error: AirflowError) -> Self { Self::Airflow(error) }
}

/// Thermal transport bound to one immutable nominal hydraulic solution.
///
/// Volume and capacity balances are recomputed from branch flows rather than
/// trusted from the solution's public summary fields. Every external source
/// needs a temperature; every branch needs an explicit model. Regions have one
/// owner across the whole network. After admission, row order is branch-major,
/// then stream-wise in SOLVED flow direction, not declared edge direction.
#[derive(Debug)]
pub struct TransportNetwork<'a> {
    flow: &'a GraphSolution,
    air: TransportAir,
    models: Vec<BranchThermalModel>,
    offsets: Vec<usize>,
    outgoing: Vec<Vec<usize>>,
    downstream: Vec<usize>,
    order: Vec<usize>,
    mass_flows: Vec<f64>,
    capacities: Vec<f64>,
    external_capacity: Vec<f64>,
    capacity_residual: Vec<f64>,
    inlet_temperatures: Vec<Option<f64>>,
    reference_k: f64,
    config: TransportConfig,
    feedback: Option<Feedback>,
}

impl<'a> TransportNetwork<'a> {
    /// Admit frozen flows, explicit branch laws, source temperatures and budgets.
    /// Heat exchange on a zero-flow branch refuses rather than inventing a
    /// stagnant-air convection model. Sink-only and interior inlet declarations
    /// also refuse. Independent anchored components and adiabatic dead legs work.
    pub fn new(
        cx: &Cx<'_>, flow: &'a GraphSolution, air: TransportAir,
        mut models: Vec<BranchThermalModel>, inlets: &[TransportInlet],
        config: TransportConfig,
    ) -> Result<Self, TransportError> {
        poll(cx)?;
        positive(air.density.value(), "air density")?;
        positive(air.specific_heat_j_kg_k, "air specific heat")?;
        positive(config.absolute_flow_tolerance.value(), "absolute flow tolerance")?;
        positive(config.absolute_heat_tolerance_w, "absolute heat tolerance")?;
        for relative in [config.relative_flow_tolerance, config.relative_heat_tolerance] {
            if !(0.0..1.0).contains(&relative) {
                return Err(TransportError::InvalidInput("relative tolerances must be in [0, 1)"));
            }
        }
        let n = flow.pressures.len();
        if n == 0 || models.len() != flow.branches.len() {
            return Err(TransportError::InvalidInput("one thermal model per branch and at least one node are required"));
        }
        let mut fixed = vec![false; n];
        for boundary in &flow.boundaries {
            poll(cx)?;
            if boundary.node >= n || !boundary.pressure.value().is_finite()
                || fixed[boundary.node] || boundary.pressure != flow.pressures[boundary.node]
            {
                return Err(TransportError::Node { node: boundary.node, reason: "invalid or duplicate hydraulic boundary" });
            }
            fixed[boundary.node] = true;
        }
        for pressure in &flow.pressures {
            finite(pressure.value(), "nodal pressure")?;
        }
        let mut outgoing = vec![Vec::new(); n];
        let mut downstream = Vec::with_capacity(models.len());
        let mut indegree = vec![0_usize; n];
        let mut volume_in = vec![0.0; n];
        let mut volume_out = vec![0.0; n];
        let mut capacity_in = vec![0.0; n];
        let mut capacity_out = vec![0.0; n];
        let mut mass_flows = Vec::with_capacity(models.len());
        let mut capacities = Vec::with_capacity(models.len());
        let mut offsets = vec![0_usize];
        let mut names = BTreeSet::new();
        for (index, (branch, model)) in flow.branches.iter().zip(&mut models).enumerate() {
            poll(cx)?;
            let bad = |reason| TransportError::Branch { branch: index, reason };
            if branch.from >= n || branch.to >= n || branch.from == branch.to {
                return Err(bad("invalid endpoints"));
            }
            let q = finite(branch.flow.value(), "branch volume flow")?;
            let (up, down) = if q < 0.0 { (branch.to, branch.from) } else { (branch.from, branch.to) };
            if q != 0.0 && flow.pressures[up].value() <= flow.pressures[down].value() {
                return Err(bad("nonzero passive flow must descend pressure"));
            }
            let mass = finite(q.abs() * air.density.value(), "branch mass flow")?;
            let capacity = finite(mass * air.specific_heat_j_kg_k, "branch capacity rate")?;
            if q != 0.0 && (mass <= 0.0 || capacity <= 0.0) {
                return Err(bad("nonzero mass or capacity flow underflowed"));
            }
            mass_flows.push(mass);
            capacities.push(capacity);
            downstream.push(down);
            if q != 0.0 {
                outgoing[up].push(index);
                indegree[down] += 1;
                volume_out[up] = finite(volume_out[up] + q.abs(), "outgoing volume flow")?;
                volume_in[down] = finite(volume_in[down] + q.abs(), "incoming volume flow")?;
                capacity_out[up] = finite(capacity_out[up] + capacity, "outgoing capacity flow")?;
                capacity_in[down] = finite(capacity_in[down] + capacity, "incoming capacity flow")?;
            }
            let count = match model {
                BranchThermalModel::Adiabatic => 0,
                BranchThermalModel::Exchange(segments) => {
                    if segments.is_empty() || q == 0.0 {
                        return Err(bad("a heat exchanger requires segments and nonzero flow"));
                    }
                    if q < 0.0 { segments.reverse(); }
                    for segment in segments.iter() {
                        poll(cx)?;
                        if !names.insert(segment.region().to_string()) {
                            return Err(bad("a solid region has more than one thermal owner"));
                        }
                    }
                    segments.len()
                }
            };
            let end = offsets.last().copied().unwrap_or(0).checked_add(count)
                .ok_or(TransportError::InvalidInput("thermal row count overflow"))?;
            offsets.push(end);
        }
        let mut external_capacity = vec![0.0; n];
        let mut capacity_residual = vec![0.0; n];
        for node in 0..n {
            poll(cx)?;
            let residual = finite(volume_out[node] - volume_in[node], "node volume residual")?;
            let capacity_difference = finite(capacity_out[node] - capacity_in[node], "node capacity residual")?;
            if fixed[node] {
                external_capacity[node] = capacity_difference;
            } else {
                let tolerance = finite(config.absolute_flow_tolerance.value()
                    + config.relative_flow_tolerance * volume_in[node].max(volume_out[node]), "flow threshold")?;
                if residual.abs() > tolerance {
                    return Err(TransportError::FlowImbalance { node, residual_m3_s: residual, tolerance_m3_s: tolerance });
                }
                capacity_residual[node] = capacity_difference;
            }
        }
        let mut inlet_temperatures = vec![None; n];
        for inlet in inlets {
            poll(cx)?;
            let node = inlet.node;
            if node >= n || external_capacity[node] <= 0.0 || inlet_temperatures[node].is_some() {
                return Err(TransportError::Node { node, reason: "temperature declarations require a unique external supply" });
            }
            inlet_temperatures[node] = Some(positive(inlet.temperature.value(), "inlet temperature")?);
        }
        for node in 0..n {
            if external_capacity[node] > 0.0 && inlet_temperatures[node].is_none() {
                return Err(TransportError::Node { node, reason: "external supply has no inlet temperature" });
            }
        }
        let mut ready: BTreeSet<usize> = indegree.iter().enumerate()
            .filter_map(|(node, &degree)| (degree == 0).then_some(node)).collect();
        let mut order = Vec::with_capacity(n);
        while let Some(node) = ready.pop_first() {
            poll(cx)?;
            order.push(node);
            for &branch in &outgoing[node] {
                poll(cx)?;
                let next = downstream[branch];
                indegree[next] -= 1;
                if indegree[next] == 0 { ready.insert(next); }
            }
        }
        if order.len() != n {
            return Err(TransportError::InvalidInput("directed recirculation is not a passive transport DAG"));
        }
        let reference_k = inlet_temperatures.iter().flatten().next().copied().unwrap_or(0.0);
        let network = Self { flow, air, models, offsets, outgoing, downstream, order,
            mass_flows, capacities, external_capacity, capacity_residual, inlet_temperatures,
            reference_k, config, feedback: None };
        // Admit the very same channel arithmetic before any caller solid work.
        for branch in 0..network.models.len() {
            poll(cx)?;
            if matches!(network.models[branch], BranchThermalModel::Exchange(_)) {
                network.path(branch, reference_k)?;
            }
        }
        Ok(network)
    }

    /// The immutable hydraulic solution supplying this network's branch flows.
    #[must_use]
    pub const fn hydraulics(&self) -> &GraphSolution { self.flow }

    /// Solid regions in branch-major, actual stream-wise order. This is the
    /// ordering for wall inputs, reference outputs, and the coupled callback.
    #[must_use]
    pub fn regions(&self) -> Vec<&str> {
        self.models.iter().flat_map(|model| match model {
            BranchThermalModel::Adiabatic => &[],
            BranchThermalModel::Exchange(segments) => segments.as_slice(),
        }).map(AirSegment::region).collect()
    }

    /// Seed all solid references with zero-heat transport of the declared inlets.
    /// Downstream branches already see any upstream mixing, not an invented
    /// uniform ambient. This is an initial guess, not a solved temperature field.
    pub fn initial_references(&self, cx: &Cx<'_>) -> Result<Vec<f64>, TransportError> {
        Ok(self.march_inner(cx, None)?.reference_temperatures_k)
    }

    /// Transport prescribed positive wall temperatures, in [`Self::regions`]
    /// order. Returns nothing on interruption or a violated global heat budget.
    pub fn march(&self, cx: &Cx<'_>, walls_k: &[f64]) -> Result<TransportMarch, TransportError> {
        if walls_k.len() != self.offsets.last().copied().unwrap_or(0) {
            return Err(TransportError::InvalidInput("one wall temperature per thermal region is required"));
        }
        for &wall in walls_k { poll(cx)?; positive(wall, "wall temperature")?; }
        self.march_inner(cx, Some(walls_k))
    }

    pub(super) fn row_range(&self, branch: usize) -> Range<usize> {
        self.offsets[branch]..self.offsets[branch + 1]
    }

    pub(super) fn path(&self, branch: usize, inlet: f64) -> Result<AirPath, TransportError> {
        let BranchThermalModel::Exchange(segments) = &self.models[branch] else {
            return Err(TransportError::InvalidInput("an adiabatic branch has no air path"));
        };
        Ok(AirPath::new(inlet, self.mass_flows[branch], self.air.specific_heat_j_kg_k, segments.clone())?)
    }

    fn march_inner(&self, cx: &Cx<'_>, walls: Option<&[f64]>) -> Result<TransportMarch, TransportError> {
        match &self.feedback {
            Some(feedback) => feedback.march(self, cx, walls),
            None => self.march_once(cx, walls, &self.inlet_temperatures),
        }
    }

    fn march_once(&self, cx: &Cx<'_>, walls: Option<&[f64]>,
        inlets: &[Option<f64>]) -> Result<TransportMarch, TransportError> {
        poll(cx)?;
        let reference_k = inlets.iter().flatten().next().copied().unwrap_or(0.0);
        let n = self.order.len();
        let mut mixers = vec![Mixer::default(); n];
        for (node, &temperature) in inlets.iter().enumerate() {
            poll(cx)?;
            if let Some(temperature) = temperature {
                mixers[node].add(self.external_capacity[node], temperature)?;
            }
        }
        let mut result = TransportMarch {
            node_temperatures_k: vec![None; n],
            branches: vec![BranchTransport { inlet_temperature_k: None, outlet_temperature_k: None, march: None }; self.models.len()],
            reference_temperatures_k: vec![0.0; self.offsets.last().copied().unwrap_or(0)],
            wall_heat_rate_w: 0.0, external_heat_gain_w: 0.0, heat_imbalance_w: 0.0,
            hydraulic_energy_defect_w: 0.0, enthalpy_reference_k: reference_k, recirculation: None,
        };
        let mut gross_wall = 0.0;
        let mut gross_external = 0.0;
        for &node in &self.order {
            poll(cx)?;
            let Some(temperature) = mixers[node].temperature else {
                if !self.outgoing[node].is_empty() {
                    return Err(TransportError::Node { node, reason: "active node has no thermally supplied incoming flow" });
                }
                continue;
            };
            result.node_temperatures_k[node] = Some(temperature);
            let external = self.external_capacity[node];
            let entering_temperature = inlets[node].unwrap_or(temperature);
            let external_heat = finite(-external * (if external > 0.0 { entering_temperature } else { temperature } - reference_k), "external sensible enthalpy")?;
            result.external_heat_gain_w = finite(result.external_heat_gain_w + external_heat, "external heat gain")?;
            gross_external = finite(gross_external + external_heat.abs(), "gross external enthalpy")?;
            result.hydraulic_energy_defect_w = finite(result.hydraulic_energy_defect_w
                + self.capacity_residual[node] * (temperature - reference_k), "hydraulic energy defect")?;
            for &branch in &self.outgoing[node] {
                poll(cx)?;
                let mut outlet = temperature;
                let mut air_march = None;
                if let BranchThermalModel::Exchange(_) = &self.models[branch] {
                    let range = self.row_range(branch);
                    if let Some(walls) = walls {
                        let marched = self.path(branch, temperature)?.march(&walls[range.clone()])?;
                        result.reference_temperatures_k[range].copy_from_slice(&marched.reference_temperatures_k());
                        outlet = marched.outlet_temperature_k;
                        result.wall_heat_rate_w = finite(result.wall_heat_rate_w + marched.total_heat_rate_w, "total wall heat")?;
                        for segment in &marched.segments {
                            poll(cx)?;
                            gross_wall = finite(gross_wall + segment.heat_rate_w.abs(), "gross wall heat")?;
                        }
                        air_march = Some(marched);
                    } else {
                        result.reference_temperatures_k[range].fill(temperature);
                    }
                }
                positive(outlet, "branch outlet temperature")?;
                mixers[self.downstream[branch]].add(self.capacities[branch], outlet)?;
                result.branches[branch] = BranchTransport {
                    inlet_temperature_k: Some(temperature), outlet_temperature_k: Some(outlet), march: air_march,
                };
            }
        }
        result.heat_imbalance_w = finite(result.external_heat_gain_w - result.wall_heat_rate_w, "global heat imbalance")?;
        let threshold = finite(self.config.absolute_heat_tolerance_w
            + self.config.relative_heat_tolerance * gross_wall.max(gross_external), "global heat threshold")?;
        if result.heat_imbalance_w.abs() > threshold {
            return Err(TransportError::HeatImbalance { residual_w: result.heat_imbalance_w,
                tolerance_w: threshold, hydraulic_defect_w: result.hydraulic_energy_defect_w });
        }
        poll(cx)?;
        Ok(result)
    }
}

#[derive(Debug, Clone, Copy, Default)]
struct Mixer { capacity: f64, temperature: Option<f64> }

impl Mixer {
    fn add(&mut self, capacity: f64, temperature: f64) -> Result<(), TransportError> {
        positive(capacity, "mixing capacity")?;
        positive(temperature, "mixed stream temperature")?;
        let total = positive(self.capacity + capacity, "total mixing capacity")?;
        // Incremental convex averaging avoids overflowing capacity * absolute K.
        self.temperature = Some(positive(match self.temperature {
            Some(previous) => previous + (temperature - previous) * (capacity / total),
            None => temperature,
        }, "mixed node temperature")?);
        self.capacity = total;
        Ok(())
    }
}

fn poll(cx: &Cx<'_>) -> Result<(), TransportError> {
    cx.checkpoint().map_err(|_| TransportError::Interrupted)
}
fn finite(value: f64, name: &'static str) -> Result<f64, TransportError> {
    if value.is_finite() { Ok(value) } else { Err(TransportError::InvalidInput(name)) }
}
fn positive(value: f64, name: &'static str) -> Result<f64, TransportError> {
    if value.is_finite() && value > 0.0 { Ok(value) } else { Err(TransportError::InvalidInput(name)) }
}

#[cfg(test)]
mod tests;
