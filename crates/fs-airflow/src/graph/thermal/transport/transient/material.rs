//! Original liquid/gas sources on branched finite-inventory convection networks.
//!
//! Native cells own M*h and T(h); original sources own rho, mu, k and Cp. This
//! binding changes transport coefficients and anomalous gravity. Reference mass,
//! continuity and inertia retain rho0. Existing PHS and thermal owners still
//! solve momentum, junction mixing, pipe residence time and reciprocal wall heat.
//!
//! Every hydraulic branch needs a nonempty, straight, monotone native pipe path.
//! Its length-weighted density then represents the integrated head exactly for
//! held section states. Curved paths refuse; unresolved buoyancy is never
//! disguised as an external pressure drive.
//!
//! Low-Eckert admission independently sums ABSOLUTE local anomalous-gravity work
//! and hydraulic dissipation at the old/new/average speed envelope. A branch-head
//! bound alone cannot certify hot/cold sections whose heads cancel. These bounds
//! concern the admitted discrete, endpoint-held model, not unresolved spatial or
//! continuous-time temperature extrema.

pub mod convection;
pub mod gas;
pub mod liquid;
pub use gas::GasGraphSources;
pub use liquid::LiquidGraphSources;

use std::collections::BTreeMap;
use std::sync::Arc;

use fs_conduction::transient::deposition::{
    CellId, DepositionError, DepositionSimulation, ThermalLaw,
};
use fs_exec::Cx;
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::state_point::MaterialStatePointError;
use fs_material::thermochemical_gas::energy::GasEnergyError;
use fs_qty::{MassFlowRate, Temperature};

use super::{TransientTransportError, distributed::DistributedGraphInventory};
use crate::graph::inertial::InertialGraphReport;
use crate::graph::thermal::medium::series::MAX_MATERIAL_FLOW_SECTIONS;
use crate::graph::thermal::medium::series::enthalpy::LiquidEnthalpyError;
use crate::graph::thermal::medium::series::hydraulic::{
    ResistanceError, ResistanceEvaluation,
    pressure::{BoussinesqDensity, SeriesPressureError},
};
use crate::graph::thermal::natural_circulation::{
    distributed_response::{
        DistributedMaterialSection, DistributedResponseError, evaluate_distributed_material_graph,
    },
    material_response::MaterialBranchLaw,
};
use crate::graph::{GraphError, LossGraph};

/// Same native geometry and declared section order as the distributed inventory.
/// Each law belongs to exactly one pipe cell and must have an empty exchange:
/// reciprocal finite-cell links, not steady exponential segments, own wall heat.
pub struct SourceGraphBinding<'a> {
    pub graph: &'a LossGraph,
    pub inventory: &'a DistributedGraphInventory,
    pub laws: &'a [Vec<MaterialBranchLaw>],
    pub acceleration_m_s2: [f64; 3],
}

/// Physical source/model limits, independent of nonlinear residual tolerances.
#[derive(Clone, Copy, Debug)]
pub struct SourceGraphLimits {
    /// Includes finite wall cells in the same native thermal problem.
    pub maximum_cells: usize,
    /// Applies to EVERY original node, pipe and evaluated film source state.
    pub maximum_relative_density_departure: f64,
    /// Hydraulic dissipation plus absolute LOCAL anomalous-gravity work.
    pub maximum_neglected_section_power_w: f64,
    pub maximum_neglected_section_energy_j: f64,
}

/// Original refusal or a disclosed source/geometry/low-Eckert limit.
#[derive(Debug)]
pub enum SourceGraphError {
    Invalid(&'static str),
    Cell {
        cell: CellId,
        reason: &'static str,
    },
    Section {
        branch: usize,
        section: usize,
        reason: &'static str,
    },
    LiquidCaloric(LiquidEnthalpyError),
    Liquid(MaterialStatePointError),
    GasCaloric(GasEnergyError),
    Material(MaterialError),
    Thermal(DepositionError),
    Inventory(TransientTransportError),
    Graph(GraphError),
    Density(SeriesPressureError),
    Compiler(DistributedResponseError),
    Resistance {
        branch: usize,
        section: usize,
        source: ResistanceError,
    },
    DensityDeparture {
        cell: CellId,
        relative: f64,
        maximum: f64,
    },
    GasMach {
        branch: usize,
        section: usize,
        value: f64,
        maximum: f64,
    },
    NeglectedWork {
        power_w: f64,
        energy_j: f64,
    },
    Cancelled,
}
impl core::fmt::Display for SourceGraphError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite graph material source: {self:?}")
    }
}
impl std::error::Error for SourceGraphError {}

/// Move sources into the inertial driver's payload and pass that same payload
/// back to admit_actual before evaluating heat.
#[derive(Debug)]
pub struct ResolvedSourceGraph<M> {
    pub graph: LossGraph,
    pub sources: GraphSourceState<M>,
}

/// Source state for one exact thermal candidate and one prepared source binding.
/// Private identity/temperature fields prevent a foreign or stale payload from
/// masquerading as the actual accepted candidate.
#[derive(Debug)]
pub struct GraphSourceState<M> {
    binding: Arc<()>,
    temperatures_k: Vec<f64>,
    node_materials: Vec<M>,
    sections: Vec<Vec<DistributedMaterialSection<M>>>,
    branch_densities_kg_m3: Vec<f64>,
    section_anomalous_heads_pa: Vec<Vec<f64>>,
    gravity_heads_pa: Vec<f64>,
}
impl<M> GraphSourceState<M> {
    /// Native junction order; original owners retain their source receipts.
    pub fn node_materials(&self) -> &[M] {
        &self.node_materials
    }
    /// Native branch declaration order, then native pipe-cell spatial order.
    pub fn sections(&self) -> &[Vec<DistributedMaterialSection<M>>] {
        &self.sections
    }
    /// Length-weighted numbers used ONLY to bind straight-branch gravity.
    pub fn branch_densities_kg_m3(&self) -> &[f64] {
        &self.branch_densities_kg_m3
    }
    /// Original local (rho_i-rho0)*g.dx_i, before branch cancellation.
    pub fn section_anomalous_heads_pa(&self) -> &[Vec<f64>] {
        &self.section_anomalous_heads_pa
    }
    pub fn gravity_heads_pa(&self) -> &[f64] {
        &self.gravity_heads_pa
    }
}

/// Independent source admission at the ACTUAL PHS average and speed envelope.
/// No mechanical work here is silently deposited into native thermal storage.
#[derive(Debug)]
pub struct GraphSourceAdmission {
    pub actual_losses: Vec<Vec<ResistanceEvaluation>>,
    pub envelope_losses: Vec<Vec<ResistanceEvaluation>>,
    pub speed_envelope_m3_s: Vec<f64>,
    pub neglected_hydraulic_power_w: f64,
    pub neglected_buoyancy_power_w: f64,
    pub neglected_power_w: f64,
    pub neglected_hydraulic_energy_j: f64,
    pub neglected_buoyancy_energy_j: f64,
    pub neglected_energy_j: f64,
}

struct BoundGraphCore<'a> {
    caloric: &'a dyn ThermalLaw,
    inventory: &'a DistributedGraphInventory,
    template: LossGraph,
    laws: Vec<Vec<MaterialBranchLaw>>,
    reference: FluidTransportState,
    positions: Vec<[f64; 3]>,
    acceleration: [f64; 3],
    limits: SourceGraphLimits,
    cells: Vec<CellId>,
    locations: BTreeMap<CellId, Option<(usize, usize)>>,
    section_lengths: Vec<Vec<f64>>,
    section_areas: Vec<Vec<f64>>,
    specific_gravity_work: Vec<Vec<f64>>,
    binding: Arc<()>,
}

impl<'a> BoundGraphCore<'a> {
    fn try_new(
        cx: &Cx<'_>,
        caloric: &'a dyn ThermalLaw,
        reference: FluidTransportState,
        thermal: &DepositionSimulation<'_>,
        binding: SourceGraphBinding<'a>,
        limits: SourceGraphLimits,
    ) -> Result<Self, SourceGraphError> {
        poll(cx)?;
        let graph = binding.graph;
        let inventory = binding.inventory;
        if limits.maximum_cells == 0
            || thermal.cells().len() > limits.maximum_cells
            || !limits.maximum_relative_density_departure.is_finite()
            || !(0.0..1.0).contains(&limits.maximum_relative_density_departure)
            || !nonnegative(limits.maximum_neglected_section_power_w)
            || !nonnegative(limits.maximum_neglected_section_energy_j)
            || binding.acceleration_m_s2.iter().any(|x| !x.is_finite())
            || graph.branches().is_empty()
            || binding.laws.len() != graph.branches().len()
            || inventory.branch_cells().len() != graph.branches().len()
            || inventory.node_cells().len() != graph.node_count()
            || reference.density_kg_m3() != inventory.reference_density().value()
        {
            return Err(SourceGraphError::Invalid(
                "bounded native inventory, complete laws and identical source reference density required",
            ));
        }
        inventory
            .admit_storage(cx, graph, thermal)
            .map_err(SourceGraphError::Inventory)?;
        let positions: Vec<_> = inventory
            .node_positions_m()
            .iter()
            .copied()
            .collect::<Option<_>>()
            .ok_or(SourceGraphError::Invalid(
                "every source-bound junction needs its original pipe endpoint position",
            ))?;
        if positions.len() != graph.node_count() {
            return Err(SourceGraphError::Invalid(
                "incomplete native junction positions",
            ));
        }
        if let Some(gravity) = graph.gravity() {
            if gravity.positions_m() != positions.as_slice()
                || gravity.acceleration_m_s2() != binding.acceleration_m_s2
                || gravity.reference_transport().density_kg_m3() != reference.density_kg_m3()
                || gravity.maximum_relative_density_departure()
                    < limits.maximum_relative_density_departure
            {
                return Err(SourceGraphError::Invalid(
                    "retain the original gravity frame, reference density and admission cap",
                ));
            }
        }
        let mut cells = inventory.node_cells().to_vec();
        let mut locations: BTreeMap<_, _> = cells.iter().map(|&cell| (cell, None)).collect();
        let mut lengths = Vec::with_capacity(binding.laws.len());
        let mut areas = Vec::with_capacity(binding.laws.len());
        let mut gravity_work = Vec::with_capacity(binding.laws.len());
        let mut section_count = 0_usize;
        for (branch, (path, laws)) in inventory
            .branch_cells()
            .iter()
            .zip(binding.laws)
            .enumerate()
        {
            poll(cx)?;
            section_count = section_count
                .checked_add(path.len())
                .ok_or(SourceGraphError::Invalid("pipe section count overflow"))?;
            if path.is_empty()
                || path.len() != laws.len()
                || section_count > MAX_MATERIAL_FLOW_SECTIONS
                || cells
                    .len()
                    .checked_add(path.len())
                    .is_none_or(|n| n > limits.maximum_cells)
            {
                return Err(SourceGraphError::Invalid(
                    "one bounded, nonempty law sequence per native pipe path required",
                ));
            }
            let edge = &graph.branches()[branch];
            let delta =
                std::array::from_fn::<_, 3, _>(|a| positions[edge.to][a] - positions[edge.from][a]);
            let total_length = norm(delta);
            if !positive(total_length) {
                return Err(SourceGraphError::Section {
                    branch,
                    section: 0,
                    reason: "straight pipe endpoints must differ",
                });
            }
            let mut branch_lengths = Vec::with_capacity(path.len());
            let mut branch_areas = Vec::with_capacity(path.len());
            let mut branch_work = Vec::with_capacity(path.len());
            for (section, (&id, law)) in path.iter().zip(laws).enumerate() {
                poll(cx)?;
                let cell = thermal.cells().get(&id).ok_or(SourceGraphError::Cell {
                    cell: id,
                    reason: "source pipe cell is not active",
                })?;
                let dx = std::array::from_fn::<_, 3, _>(|a| {
                    cell.geometry.end_m[a] - cell.geometry.start_m[a]
                });
                let length = norm(dx);
                let area = cell.geometry.volume_m3 / length;
                // Parallel AND same orientation. Bent gravity is not silently
                // averaged or attributed to an external drive.
                if !positive(length)
                    || !positive(area)
                    || !law.exchange.is_empty()
                    || !near_positive(area, law.resistance.flow_area().value())
                    || (0..3).any(|a| {
                        (dx[a] - delta[a] * (length / total_length)).abs()
                            > 256.0 * f64::EPSILON * length
                    })
                {
                    return Err(SourceGraphError::Section {
                        branch,
                        section,
                        reason: "native area, adiabatic drag law and straight monotone path required",
                    });
                }
                let work: f64 = dx
                    .iter()
                    .zip(binding.acceleration_m_s2)
                    .map(|(x, g)| x * g)
                    .sum();
                if !work.is_finite() {
                    return Err(SourceGraphError::Section {
                        branch,
                        section,
                        reason: "local gravity work is unrepresentable",
                    });
                }
                cells.push(id);
                locations.insert(id, Some((branch, section)));
                branch_lengths.push(length);
                branch_areas.push(area);
                branch_work.push(work);
            }
            lengths.push(branch_lengths);
            areas.push(branch_areas);
            gravity_work.push(branch_work);
        }
        let core = Self {
            caloric,
            inventory,
            template: graph.clone(),
            laws: binding.laws.to_vec(),
            reference,
            positions,
            acceleration: binding.acceleration_m_s2,
            limits,
            cells,
            locations,
            section_lengths: lengths,
            section_areas: areas,
            specific_gravity_work: gravity_work,
            binding: Arc::new(()),
        };
        core.temperatures(cx, thermal)?;
        Ok(core)
    }

    fn temperatures(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<Vec<f64>, SourceGraphError> {
        poll(cx)?;
        if thermal.cells().len() > self.limits.maximum_cells {
            return Err(SourceGraphError::Invalid(
                "native thermal cell budget exceeded",
            ));
        }
        self.inventory
            .admit_storage(cx, &self.template, thermal)
            .map_err(SourceGraphError::Inventory)?;
        let mut temperatures = Vec::with_capacity(self.cells.len());
        for &id in &self.cells {
            poll(cx)?;
            let cell = thermal.cells().get(&id).ok_or(SourceGraphError::Cell {
                cell: id,
                reason: "bound source cell is not active",
            })?;
            if !thermal.uses_material_law(cell.material_id, self.caloric)
                || !near_positive(
                    cell.mass_kg,
                    self.reference.density_kg_m3() * cell.geometry.volume_m3,
                )
            {
                return Err(SourceGraphError::Cell {
                    cell: id,
                    reason: "exact original caloric law and fixed reference mass rho0*V required",
                });
            }
            let t = thermal
                .temperature_k(id)
                .map_err(SourceGraphError::Thermal)?;
            if !positive(t) {
                return Err(SourceGraphError::Cell {
                    cell: id,
                    reason: "invalid native caloric temperature",
                });
            }
            temperatures.push(t);
        }
        Ok(temperatures)
    }

    fn admit_density(
        &self,
        cell: CellId,
        state: FluidTransportState,
    ) -> Result<(), SourceGraphError> {
        let departure = ((state.density_kg_m3() - self.reference.density_kg_m3())
            / self.reference.density_kg_m3())
        .abs();
        if !departure.is_finite() || departure > self.limits.maximum_relative_density_departure {
            return Err(SourceGraphError::DensityDeparture {
                cell,
                relative: departure,
                maximum: self.limits.maximum_relative_density_departure,
            });
        }
        Ok(())
    }

    fn resolve<M: FluidTransport>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        flows: &[f64],
        mut resolve: impl FnMut(Temperature) -> Result<M, SourceGraphError>,
    ) -> Result<ResolvedSourceGraph<M>, SourceGraphError> {
        let temperatures = self.temperatures(cx, thermal)?;
        if flows.len() != self.laws.len() || flows.iter().any(|q| !q.is_finite()) {
            return Err(SourceGraphError::Invalid(
                "one finite signed trial flow per native branch required",
            ));
        }
        let mut materials = Vec::with_capacity(self.cells.len());
        for (&id, &t) in self.cells.iter().zip(&temperatures) {
            poll(cx)?;
            let material = resolve(Temperature::new(t))?;
            let transport = material
                .transport_state()
                .map_err(SourceGraphError::Material)?;
            self.admit_density(id, transport)?;
            materials.push(material);
        }
        let mut all = materials.into_iter();
        let node_materials = all
            .by_ref()
            .take(self.inventory.node_cells().len())
            .collect();
        let mut bulk = Vec::with_capacity(self.laws.len());
        let mut empty_heat = Vec::with_capacity(self.laws.len());
        for path in self.inventory.branch_cells() {
            bulk.push(all.by_ref().take(path.len()).collect::<Vec<_>>());
            empty_heat.push((0..path.len()).map(|_| Vec::new()).collect::<Vec<_>>());
        }
        let density = BoussinesqDensity::try_new(
            self.reference.density_kg_m3(),
            self.limits.maximum_relative_density_departure,
        )
        .map_err(SourceGraphError::Density)?;
        let compiled = evaluate_distributed_material_graph(
            cx,
            &self.template,
            density,
            &self.laws,
            flows,
            bulk,
            empty_heat,
        )
        .map_err(SourceGraphError::Compiler)?;
        let mut branch_densities = Vec::with_capacity(self.laws.len());
        let mut heads = Vec::with_capacity(self.laws.len());
        let mut branch_snapshots = Vec::with_capacity(self.laws.len());
        for (branch, sections) in compiled.sections.iter().enumerate() {
            poll(cx)?;
            let total: f64 = self.section_lengths[branch].iter().sum();
            if !positive(total) {
                return Err(SourceGraphError::Invalid(
                    "native straight branch length is unrepresentable",
                ));
            }
            let mut departure = 0.0;
            let mut local_heads = Vec::with_capacity(sections.len());
            for (section, source) in sections.iter().enumerate() {
                poll(cx)?;
                let rho = source.momentum.evaluation.transport.density_kg_m3();
                let delta = rho - self.reference.density_kg_m3();
                departure += delta * (self.section_lengths[branch][section] / total);
                let head = delta * self.specific_gravity_work[branch][section];
                if !head.is_finite() {
                    return Err(SourceGraphError::Invalid(
                        "local anomalous head is unrepresentable",
                    ));
                }
                local_heads.push(head);
            }
            let rho = self.reference.density_kg_m3() + departure;
            // Numerical gravity-only tuple. Original material owners remain
            // authoritative; this tuple never supplies calorics or drag.
            branch_snapshots.push(
                FluidTransportState::try_new(
                    rho,
                    self.reference.specific_heat_j_kg_k(),
                    self.reference.thermal_conductivity_w_m_k(),
                    self.reference.dynamic_viscosity_pa_s(),
                )
                .map_err(SourceGraphError::Material)?,
            );
            branch_densities.push(rho);
            heads.push(local_heads);
        }
        let graph = compiled
            .graph
            .with_boussinesq_gravity(
                &self.reference,
                &branch_snapshots,
                &self.positions,
                self.acceleration,
                self.limits.maximum_relative_density_departure,
                cx,
            )
            .map_err(SourceGraphError::Graph)?;
        let gravity_heads_pa = graph
            .gravity()
            .ok_or(SourceGraphError::Invalid("missing compiled source gravity"))?
            .pressure_heads_pa()
            .to_vec();
        Ok(ResolvedSourceGraph {
            graph,
            sources: GraphSourceState {
                binding: self.binding.clone(),
                temperatures_k: temperatures,
                node_materials,
                sections: compiled.sections,
                branch_densities_kg_m3: branch_densities,
                section_anomalous_heads_pa: heads,
                gravity_heads_pa,
            },
        })
    }

    fn admit_actual<M>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        sources: &GraphSourceState<M>,
        report: &InertialGraphReport,
    ) -> Result<GraphSourceAdmission, SourceGraphError> {
        let temperatures = self.temperatures(cx, thermal)?;
        if !Arc::ptr_eq(&self.binding, &sources.binding) || temperatures != sources.temperatures_k {
            return Err(SourceGraphError::Invalid(
                "foreign prepared source binding or stale thermal source payload",
            ));
        }
        let count = self.laws.len();
        if !positive(report.duration_s)
            || [
                report.average_volume_flows_m3_s.len(),
                report.initial_flows_m3_s.len(),
                report.final_flows_m3_s.len(),
                report.speed_envelope_m3_s.len(),
                report.dissipative_pressure_drops_pa.len(),
                report.gravity_pressure_heads_pa.len(),
            ]
            .iter()
            .any(|&n| n != count)
        {
            return Err(SourceGraphError::Invalid(
                "complete native inertial report and positive duration required",
            ));
        }
        let mut actual_losses = Vec::with_capacity(count);
        let mut envelope_losses = Vec::with_capacity(count);
        let mut speed_envelope = Vec::with_capacity(count);
        let mut hydraulic_power = 0.0;
        let mut buoyancy_power = 0.0;
        for branch in 0..count {
            poll(cx)?;
            let actual = report.average_volume_flows_m3_s[branch];
            let initial = report.initial_flows_m3_s[branch];
            let final_flow = report.final_flows_m3_s[branch];
            if [actual, initial, final_flow].iter().any(|v| !v.is_finite()) {
                return Err(SourceGraphError::Invalid(
                    "nonfinite actual native momentum flow",
                ));
            }
            let envelope = actual.abs().max(initial.abs()).max(final_flow.abs());
            if !roundoff_equal(envelope, report.speed_envelope_m3_s[branch], 1)
                || !roundoff_equal(
                    sources.gravity_heads_pa[branch],
                    report.gravity_pressure_heads_pa[branch],
                    1,
                )
            {
                return Err(SourceGraphError::Invalid(
                    "native report differs from source gravity or speed envelope",
                ));
            }
            let mut actual_branch = Vec::with_capacity(self.laws[branch].len());
            let mut envelope_branch = Vec::with_capacity(self.laws[branch].len());
            let mut drop = 0.0;
            for (section, source) in sources.sections[branch].iter().enumerate() {
                poll(cx)?;
                let law = self.laws[branch][section].resistance;
                let transport = source.momentum.evaluation.transport;
                let rho0 = self.reference.density_kg_m3();
                let evaluate = |q: f64| {
                    law.evaluate_snapshot(transport, MassFlowRate::new(rho0 * q), rho0)
                        .map_err(|source| SourceGraphError::Resistance {
                            branch,
                            section,
                            source,
                        })
                };
                let actual_loss = evaluate(actual)?;
                // Support may have a nonzero lower Reynolds limit. Maximum
                // speed alone cannot qualify an accelerating or reversing step.
                evaluate(initial)?;
                evaluate(final_flow)?;
                if (initial < 0.0 && final_flow > 0.0) || (initial > 0.0 && final_flow < 0.0) {
                    evaluate(0.0)?;
                }
                let envelope_loss = evaluate(envelope)?;
                drop += actual_loss.pressure_drop_pa;
                hydraulic_power += envelope_loss.dissipated_power_w;
                buoyancy_power +=
                    sources.section_anomalous_heads_pa[branch][section].abs() * envelope;
                actual_branch.push(actual_loss);
                envelope_branch.push(envelope_loss);
            }
            if !roundoff_equal(
                drop,
                report.dissipative_pressure_drops_pa[branch],
                self.laws[branch].len(),
            ) {
                return Err(SourceGraphError::Invalid(
                    "actual native drag differs from retained original section laws",
                ));
            }
            actual_losses.push(actual_branch);
            envelope_losses.push(envelope_branch);
            speed_envelope.push(envelope);
        }
        let power = hydraulic_power + buoyancy_power;
        let hydraulic_energy = report.duration_s * hydraulic_power;
        let buoyancy_energy = report.duration_s * buoyancy_power;
        let energy = report.duration_s * power;
        if [
            hydraulic_power,
            buoyancy_power,
            power,
            hydraulic_energy,
            buoyancy_energy,
            energy,
        ]
        .into_iter()
        .any(|v| !nonnegative(v))
        {
            return Err(SourceGraphError::Invalid(
                "non-cancelling local omitted-work bound is unrepresentable",
            ));
        }
        if power > self.limits.maximum_neglected_section_power_w
            || energy > self.limits.maximum_neglected_section_energy_j
        {
            return Err(SourceGraphError::NeglectedWork {
                power_w: power,
                energy_j: energy,
            });
        }
        Ok(GraphSourceAdmission {
            actual_losses,
            envelope_losses,
            speed_envelope_m3_s: speed_envelope,
            neglected_hydraulic_power_w: hydraulic_power,
            neglected_buoyancy_power_w: buoyancy_power,
            neglected_power_w: power,
            neglected_hydraulic_energy_j: hydraulic_energy,
            neglected_buoyancy_energy_j: buoyancy_energy,
            neglected_energy_j: energy,
        })
    }
}

fn poll(cx: &Cx<'_>) -> Result<(), SourceGraphError> {
    cx.checkpoint().map_err(|_| SourceGraphError::Cancelled)
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn norm(v: [f64; 3]) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}
fn near_positive(a: f64, b: f64) -> bool {
    positive(a) && positive(b) && (a - b).abs() <= 128.0 * f64::EPSILON * a.abs().max(b.abs())
}
fn roundoff_equal(a: f64, b: f64, terms: usize) -> bool {
    a.is_finite()
        && b.is_finite()
        && (a - b).abs()
            <= 256.0
                * f64::EPSILON
                * (1.0 + terms as f64)
                * a.abs().max(b.abs()).max(f64::MIN_POSITIVE)
}
