//! Regional finite-inventory heat transport through original Carreau pipes.
//!
//! One native thermal cell represents each physical developed pipe section;
//! separate native cells own the finite junction volumes. Every implicit trial
//! re-prepares each original pipe at its own cell temperature and resolves the
//! actual nonlinear pressure network. The native upwind enthalpy owner carries
//! that signed flow through junction -> pipe -> junction, and reverses donors
//! without discarding fluid already stored in the pipe. Only each original
//! steady profile's viscous dissipation is deposited as local fluid heat.
//! Explicit straight series splices transfer directly between real pipe cells;
//! they retain hydraulic pressure nodes but add no mixing or thermal storage.

mod series;
pub mod transient;

use fs_exec::Cx;
use fs_qty::Density;
use fs_thermal::enthalpy::advection::MassBoundary;
use fs_thermal::enthalpy::backward_euler::{
    BoundaryStepCoefficients, StepCoefficients, TransportStepCoefficients,
};
use fs_thermal::enthalpy::{CellId, DepositionSimulation};

use super::distributed::MAX_DISTRIBUTED_FLUID_CELLS;
use super::nonlinear::{
    ImplicitTransportControls, ImplicitTransportError, ImplicitTransportReport,
    advance_endpoint_transport, combine_endpoint_heat,
};
use super::{GraphEnthalpyTransport, thermal_checkpoint};
use crate::graph::FixedPressure;
use crate::graph::rheology::{
    CarreauGraphControls, CarreauGraphError, CarreauGraphSolution, CarreauPipeGraph, FlowInjection,
};

/// Original hydraulic or caller-owned exchange refusal at a thermal trial.
#[derive(Debug)]
pub enum CarreauTransportError<E> {
    Hydraulic(CarreauGraphError),
    Exchange(E),
}
impl<E: core::fmt::Debug> core::fmt::Display for CarreauTransportError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "regional Carreau transport: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for CarreauTransportError<E> {}

/// Explicit open streams and source-owned heat at the SAME thermal/flow state.
pub struct CarreauTransportExchange<T> {
    /// Streams attach only to original junction cells. Their native rates,
    /// caloric identities and donor enthalpies must balance actual solved Q.
    /// A prescribed pressure does not implicitly supply any mass or enthalpy.
    pub mass_boundaries: Vec<MassBoundary>,
    /// Additional sources, reciprocal wall links and ambient coefficients.
    /// Do not include advection, pressure work or pipe viscous dissipation:
    /// this adapter already installs the original transport and viscous heat.
    pub heat: BoundaryStepCoefficients<T>,
}

/// Original hydraulic state and distinct heat contributions at acceptance.
#[derive(Debug)]
pub struct CarreauTransportEndpoint<T> {
    pub hydraulics: CarreauGraphSolution,
    pub exchange: T,
    /// Original accepted-profile losses, in original pipe order, W.
    pub branch_viscous_heat_w: Vec<f64>,
    /// Native datum-reduced NET external stream enthalpy, W. Internal
    /// advection is conservative and is excluded from this number.
    pub external_stream_enthalpy_w: f64,
}

/// Immutable original topology and finite thermal storage binding.
///
/// The caloric law belongs to the original native thermal timeline. Every
/// fluid cell shares one registered enthalpy material and declared caloric
/// identity; the caller associates that original caloric law with the original
/// Carreau fluid. Density is an explicit fixed reference, never inferred from
/// rheology or a heat capacity. Native walls/solids may coexist outside this map.
#[derive(Clone)]
pub struct RegionalCarreauTransport {
    graph: CarreauPipeGraph,
    series: series::SeriesInventory,
    pipe_cells: Vec<CellId>,
    reference_density: Density,
    inventory: GraphEnthalpyTransport,
    maximum_connections: usize,
}

impl RegionalCarreauTransport {
    /// Bind one complete physical pipe cell per original branch and one
    /// distinct, authored finite junction cell per node. Pipe volume must be
    /// pi R²L and its straight centreline length must be L; incident pipe
    /// endpoints must agree at each hydraulic junction. Junction storage is
    /// additional real inventory, never a copy of the adjacent pipe volume.
    ///
    /// All fluid masses must equal rho0 V. The original native storage owner
    /// pins mass, full geometry, material and caloric identity. Initial local
    /// temperatures must admit the original unchanged Carreau rate ceilings.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        cx: &Cx<'_>,
        graph: CarreauPipeGraph,
        reference_density: Density,
        node_cells: &[CellId],
        pipe_cells: &[CellId],
        thermal: &DepositionSimulation<'_>,
        maximum_connections: usize,
    ) -> Result<Self, ImplicitTransportError<CarreauGraphError>> {
        Self::new_with_series_splices(
            cx,
            graph,
            reference_density,
            &node_cells.iter().copied().map(Some).collect::<Vec<_>>(),
            pipe_cells,
            thermal,
            maximum_connections,
        )
    }

    /// Bind real pipe and junction inventory, with explicit zero-volume series
    /// splices. `Some(cell)` retains an authored finite junction; `None` declares
    /// a straight continuation of exactly two same-radius, same-original-law
    /// pipe sections. Their native cell centrelines must meet and continue in
    /// opposite outward directions. No extra volume or mixing cell is created.
    ///
    /// Splices remain original hydraulic pressure unknowns and may not carry
    /// prescribed pressure or injection ports. All external mass streams still
    /// attach to finite junction cells. The original per-section developed-flow
    /// approximation and caloric/source domains remain in force; this is not
    /// an entrance-flow or axial-momentum-convection model.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_series_splices(
        cx: &Cx<'_>,
        graph: CarreauPipeGraph,
        reference_density: Density,
        node_cells: &[Option<CellId>],
        pipe_cells: &[CellId],
        thermal: &DepositionSimulation<'_>,
        maximum_connections: usize,
    ) -> Result<Self, ImplicitTransportError<CarreauGraphError>> {
        thermal_checkpoint(cx)?;
        let count = node_cells
            .iter()
            .flatten()
            .count()
            .checked_add(pipe_cells.len());
        if node_cells.len() != graph.node_count()
            || pipe_cells.len() != graph.branches().len()
            || count.is_none_or(|n| n > MAX_DISTRIBUTED_FLUID_CELLS)
            || !(1..=1_000_000).contains(&maximum_connections)
        {
            return Err(ImplicitTransportError::Invalid(
                "bounded complete Carreau pipe and junction inventories required",
            ));
        }
        let mut positions = vec![None; node_cells.len()];
        let mut temperatures = Vec::with_capacity(pipe_cells.len());
        for (branch, &id) in graph.branches().iter().zip(pipe_cells) {
            thermal_checkpoint(cx)?;
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(ImplicitTransportError::Invalid(
                    "Carreau pipe requires an active native thermal cell",
                ))?;
            let geometry = cell.geometry;
            let length = fs_math::det::hypot(
                fs_math::det::hypot(
                    geometry.end_m[0] - geometry.start_m[0],
                    geometry.end_m[1] - geometry.start_m[1],
                ),
                geometry.end_m[2] - geometry.start_m[2],
            );
            let pipe = branch.pipe;
            let volume =
                core::f64::consts::PI * pipe.radius_m() * pipe.radius_m() * pipe.length_m();
            if !same_positive(length, pipe.length_m()) || !same_positive(volume, geometry.volume_m3)
            {
                return Err(ImplicitTransportError::Invalid(
                    "native Carreau pipe length and volume must match original hydraulic geometry",
                ));
            }
            for (node, point) in [(branch.from, geometry.start_m), (branch.to, geometry.end_m)] {
                if let Some(previous) = positions[node] {
                    if previous != point {
                        return Err(ImplicitTransportError::Invalid(
                            "incident Carreau pipe endpoints must share their original junction",
                        ));
                    }
                } else {
                    positions[node] = Some(point);
                }
            }
            temperatures.push(thermal.temperature_k(id)?);
        }
        let series = series::SeriesInventory::new(cx, &graph, node_cells, pipe_cells, thermal)?;
        if series.connections.len() > maximum_connections {
            return Err(ImplicitTransportError::Invalid(
                "bounded Carreau inventory connections required",
            ));
        }
        let mut all_cells = series.junctions.clone();
        all_cells.extend_from_slice(pipe_cells);
        // This also admits distinct cells, enthalpy storage, reference density
        // and the actual immutable caloric identities, even at initial rest.
        let inventory = GraphEnthalpyTransport::bind(
            cx,
            all_cells.len(),
            None,
            series
                .connections
                .iter()
                .map(|link| (link.from, link.to, 0.0)),
            reference_density,
            &all_cells,
            thermal,
            vec![],
            maximum_connections,
        )?;
        let material = thermal.cells()[&all_cells[0]].material_id;
        if thermal.material_caloric_identity(material)?.is_none()
            || all_cells
                .iter()
                .any(|id| thermal.cells()[id].material_id != material)
        {
            return Err(ImplicitTransportError::Invalid(
                "one original identity-bearing enthalpy material required for all Carreau fluid cells",
            ));
        }
        let graph = graph
            .at_temperatures(&temperatures, cx)
            .map_err(ImplicitTransportError::Resolve)?;
        thermal_checkpoint(cx)?;
        Ok(Self {
            graph,
            series,
            pipe_cells: pipe_cells.to_vec(),
            reference_density,
            inventory,
            maximum_connections,
        })
    }

    /// Immutable source/topology template at the binding's initial temperatures.
    /// The accepted current temperatures and profiles belong to the step report.
    pub fn graph(&self) -> &CarreauPipeGraph {
        &self.graph
    }
    /// Finite junction cells in original hydraulic node order, omitting declared
    /// splices. For `new` this remains the complete original node-to-cell map.
    /// Use `node_cell` or `node_cell_bindings` when hydraulic indices are needed.
    pub fn node_cells(&self) -> &[CellId] {
        &self.series.junctions
    }
    /// Explicit inventory binding for every original hydraulic node.
    pub fn node_cell_bindings(&self) -> &[Option<CellId>] {
        &self.series.nodes
    }
    /// A node's actual inventory, or `None` for a splice or an absent node.
    pub fn node_cell(&self, node: usize) -> Option<CellId> {
        self.series.nodes.get(node).copied().flatten()
    }
    pub fn pipe_cells(&self) -> &[CellId] {
        &self.pipe_cells
    }
    pub fn reference_density(&self) -> Density {
        self.reference_density
    }

    /// Advance native regional enthalpy and quasi-steady original rheology
    /// together. Every residual, derivative probe and final acceptance resolves
    /// local source viscosity, pressure/flow, upwind transport and wall exchange
    /// at that SAME candidate. The exchange callback must be pure and admit
    /// its source laws on every queried state, including derivative probes.
    /// Callback refusals and work exhaustion are terminal.
    ///
    /// The original nonlinear thermal driver owns all stored energies, native
    /// boundary events, local/global energy gates and publication. Inputs stay
    /// unchanged on failure/cancellation. The returned payload retains actual
    /// original hydraulic residual work separately from physical viscous heat.
    /// Native mass closure is stricter than a loose hydraulic stopping criterion;
    /// choose hydraulic tolerances that satisfy the original transport gate.
    ///
    /// This is fixed-density, incompressible, axially developed laminar flow.
    /// Each section is thermally well mixed; author multiple real sections to
    /// resolve residence and spatial viscosity. Endpoint flow represents this
    /// entire backward-Euler interval. No momentum history, expansion, radial
    /// temperature field, entrance loss, buoyancy or reversal-time location is
    /// inferred. Refine space/time and retain the original laminar/source bounds.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<'m, T, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        boundaries: &[FixedPressure],
        injections: &[FlowInjection],
        hydraulic_controls: CarreauGraphControls,
        controls: ImplicitTransportControls,
        energy_guess_j: &[f64],
        exchange: impl Fn(
            &Cx<'_>,
            &DepositionSimulation<'m>,
            &CarreauGraphSolution,
        ) -> Result<CarreauTransportExchange<T>, E>,
    ) -> Result<
        (
            DepositionSimulation<'m>,
            ImplicitTransportReport<CarreauTransportEndpoint<T>>,
        ),
        ImplicitTransportError<CarreauTransportError<E>>,
    > {
        thermal_checkpoint(cx)?;
        self.inventory.validate_storage(cx, thermal)?;
        self.series.validate_ports(cx, boundaries, injections)?;
        advance_endpoint_transport(
            cx,
            thermal,
            target_s,
            controls,
            energy_guess_j,
            |cx, candidate| {
                self.inventory.validate_storage(cx, candidate)?;
                let mut temperatures = Vec::with_capacity(self.pipe_cells.len());
                for &cell in &self.pipe_cells {
                    thermal_checkpoint(cx)?;
                    temperatures.push(candidate.temperature_k(cell)?);
                }
                let hydraulic_error = |error| {
                    ImplicitTransportError::Resolve(CarreauTransportError::Hydraulic(error))
                };
                let local = self
                    .graph
                    .at_temperatures(&temperatures, cx)
                    .map_err(hydraulic_error)?;
                let hydraulics = local
                    .solve_with_injections(boundaries, injections, hydraulic_controls, cx)
                    .map_err(hydraulic_error)?;
                let response = exchange(cx, candidate, &hydraulics).map_err(|error| {
                    ImplicitTransportError::Resolve(CarreauTransportError::Exchange(error))
                })?;
                thermal_checkpoint(cx)?;
                let rates: Vec<_> = hydraulics
                    .branches
                    .iter()
                    .map(|branch| branch.profile.flow().volumetric_flow_m3_s)
                    .collect();
                let viscous_heat_w: Vec<_> = hydraulics
                    .branches
                    .iter()
                    .map(|branch| branch.profile.flow().dissipated_power_w)
                    .collect();
                self.endpoint_coefficients(
                    cx,
                    candidate,
                    &rates,
                    viscous_heat_w,
                    response,
                    |payload, external_stream_enthalpy_w, branch_viscous_heat_w| {
                        CarreauTransportEndpoint {
                            hydraulics,
                            exchange: payload,
                            branch_viscous_heat_w,
                            external_stream_enthalpy_w,
                        }
                    },
                )
            },
        )
    }
    /// Shared original inventory, signed mass binding and local heat deposit.
    /// The momentum owner supplies the actual rates and irreversible losses;
    /// this adapter does not assume a steady pressure-flow identity.
    #[allow(clippy::too_many_arguments)]
    fn endpoint_coefficients<T, P, E>(
        &self,
        cx: &Cx<'_>,
        candidate: &DepositionSimulation<'_>,
        signed_volume_flows_m3_s: &[f64],
        branch_viscous_heat_w: Vec<f64>,
        response: CarreauTransportExchange<T>,
        finish: impl FnOnce(T, f64, Vec<f64>) -> P,
    ) -> Result<TransportStepCoefficients<P>, ImplicitTransportError<E>> {
        self.inventory.validate_storage(cx, candidate)?;
        if signed_volume_flows_m3_s.len() != self.pipe_cells.len()
            || branch_viscous_heat_w.len() != self.pipe_cells.len()
        {
            return Err(ImplicitTransportError::Invalid(
                "one actual flow and viscous heat per original Carreau pipe required",
            ));
        }
        for boundary in &response.mass_boundaries {
            thermal_checkpoint(cx)?;
            let cell = match boundary {
                MassBoundary::Inlet { cell, .. } | MassBoundary::Outlet { cell, .. } => cell,
            };
            if !self.node_cells().contains(cell) {
                return Err(ImplicitTransportError::Invalid(
                    "external Carreau streams must attach to original junction cells",
                ));
            }
        }
        self.series
            .validate_flows(cx, self.reference_density, signed_volume_flows_m3_s)?;
        let binding = GraphEnthalpyTransport::bind(
            cx,
            self.inventory.node_cells().len(),
            None,
            self.series.connections.iter().map(|link| {
                (
                    link.from,
                    link.to,
                    link.sign * signed_volume_flows_m3_s[link.branch],
                )
            }),
            self.reference_density,
            self.inventory.node_cells(),
            candidate,
            response.mass_boundaries,
            self.maximum_connections,
        )?;
        let advected = binding.endpoint_transport(cx, candidate)?;
        let mut stream = crate::graph::Sum::default();
        let arithmetic =
            |_| ImplicitTransportError::Invalid("unrepresentable native external-stream enthalpy");
        for &watts in advected.heat.bulk.sources_w.values() {
            stream.add(watts).map_err(arithmetic)?;
        }
        let external_stream_enthalpy_w = stream.value().map_err(arithmetic)?;
        let BoundaryStepCoefficients { bulk, boundaries } = response.heat;
        let StepCoefficients {
            mut sources_w,
            conductances_w_k,
            payload,
        } = bulk;
        for (&cell, &heat) in self.pipe_cells.iter().zip(&branch_viscous_heat_w) {
            thermal_checkpoint(cx)?;
            let total = sources_w.entry(cell).or_insert(0.0);
            *total += heat;
            if !heat.is_finite() || heat < 0.0 || !total.is_finite() {
                return Err(ImplicitTransportError::Invalid(
                    "unrepresentable original Carreau viscous heat",
                ));
            }
        }
        combine_endpoint_heat(
            cx,
            advected,
            BoundaryStepCoefficients {
                bulk: StepCoefficients {
                    sources_w,
                    conductances_w_k,
                    payload: finish(payload, external_stream_enthalpy_w, branch_viscous_heat_w),
                },
                boundaries,
            },
        )
    }
}

fn same_positive(a: f64, b: f64) -> bool {
    a.is_finite()
        && b.is_finite()
        && a > 0.0
        && b > 0.0
        && (a - b).abs() <= 64.0 * f64::EPSILON * a.max(b)
}

#[cfg(test)]
mod tests;
