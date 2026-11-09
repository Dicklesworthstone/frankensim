//! Source-backed storage and graph resistance inside nonlinear spatial closure.
//!
//! The native DAE, field thermodynamics and finite-heat residual remain their
//! original owners. Each thermal trial refreshes R/C/L and contact states once;
//! every nested bulk constitutive solve sees those same coefficients and starts
//! from the same accepted charge and flux. Reversible field transfer is kept
//! separate from bulk and interface Joule heating.

use super::*;
use crate::circuit_enthalpy::field_storage::{self, FieldThermalCandidate};
use crate::circuit_enthalpy::{
    conduction::heat::ConductionHeatError,
    contacts::heat::{HeatAssembly, HeatNetworkError, closed_balance},
};
use crate::material_circuit::storage::lower;
use fs_em::conduction::motion::ConductionMotion;
use fs_exec::Cx;

/// Nonlinear spatial response, source-backed field storage and one finite heat
/// history. The physical energy account uses U = psi + T s at both endpoints.
pub type ConductionResponseFieldThermalCandidate<'t, S = MaterialNetworkSample> =
    FieldThermalCandidate<'t, ConductionResponseCircuit, ConductionResponseCircuitStep, S>;

impl ConductionResponseCircuit {
    /// Bind source-derived R/C/L to the same native circuit as the nonlinear
    /// spatial body. Every nonspatial resistor and every storage branch requires
    /// its original complete material/geometry owner. Graph numeric values are
    /// replaced by those owners before admitting the initial operating point.
    ///
    /// Conductivity is queried at the actual initial electric field/current,
    /// including native surface contacts. Initial conserved charge/flux are not
    /// changed. Coefficients stay fixed in an ordinary fixed-context step;
    /// [`ConductionThermalNetwork::step_response_field_circuit_candidate`]
    /// refreshes their thermal state and accounts for reversible field entropy.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_material_storage(
        problem: PlanarConductionResponse,
        graph: CircuitGraph,
        terminal: ConductionCircuitTerminal,
        mut material_resistors: Vec<MaterialResistor>,
        mut storage: Vec<MaterialStorage>,
        time_s: f64,
        fluxes_wb: Vec<f64>,
        charges_c: Vec<f64>,
        ports: &[f64],
        assembly: WeightedAssemblyLimits,
        field: ConductionControls,
        response: ConductionResponseControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        let graph = lower(&graph, &mut storage)?;
        let graph = stamp_material_resistors(&graph, &terminal, &mut material_resistors)?;
        let mut circuit = Self::new(
            problem,
            graph,
            terminal,
            time_s,
            fluxes_wb,
            charges_c,
            ports,
            assembly,
            field,
            response,
            &mut cancelled,
        )?;
        circuit.circuit.material_resistors = material_resistors;
        circuit.circuit.storage = storage;
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        Ok(circuit)
    }

    /// Original whole-cell graph-resistor sources at this accepted endpoint.
    pub fn material_resistors(&self) -> &[MaterialResistor] {
        self.circuit.material_resistors()
    }

    /// Original dielectric/core material sources and fixed field geometry.
    pub fn storage(&self) -> &[MaterialStorage] {
        self.circuit.storage()
    }
}

impl<'t, S> ConductionResponseCircuitEnthalpyCandidate<'t, FieldThermalSource<S>> {
    // Plain and complete heat networks close physical field U through the
    // original energy owner. Internal enclosure roundoff is never external heat.
    pub(in crate::circuit_enthalpy::conduction) fn close_field_report(
        self,
        internal_heat_j: f64,
        tolerance_j: f64,
    ) -> Result<ConductionResponseFieldThermalCandidate<'t, S>, Error> {
        let (thermal, circuit, coupled) = self.into_parts();
        let supplied = finite(
            coupled.electrical.circuit.supplied_work_j - coupled.electrical.mechanical_work_j,
        )?;
        let report = field_storage::close_report_with_internal_heat(
            coupled,
            supplied,
            internal_heat_j,
            tolerance_j,
        )
        .map_err(super::super::thermal::field_thermal_error)?;
        Ok(FieldThermalCandidate::new(thermal, circuit, report))
    }

    pub(super) fn without_thermal_fields(
        self,
    ) -> ConductionResponseCircuitEnthalpyCandidate<'t, S> {
        let report = self.report;
        ConductionResponseCircuitEnthalpyCandidate {
            thermal: self.thermal,
            circuit: self.circuit,
            report: CircuitEnthalpyReport {
                electrical: report.electrical,
                thermal: report.thermal,
                source: report.source.network,
                nonlinear: report.nonlinear,
                maximum_cell_residual_j: report.maximum_cell_residual_j,
                circuit_trials: report.circuit_trials,
                deposited_joule_j: report.deposited_joule_j,
                non_electrical_heat_j: report.non_electrical_heat_j,
                total_energy_residual_j: report.total_energy_residual_j,
            },
        }
    }
}

impl ConductionThermalNetwork<'_> {
    /// Evolve nonlinear spatial sigma(T, |E|, |J|), electrical contacts, source
    /// graph R(T), source C(T)/L(T), and finite caloric storage in one candidate.
    ///
    /// All source coefficients are refreshed at the common thermal endpoint.
    /// Nested constitutive solves use actual midpoint fields and restart from
    /// accepted charge/flux. Actual bulk/surface Joule losses enter their named
    /// cells once; signed -delta(psi_parameter)-delta(T s) goes to each original
    /// field material's caloric cell. The final enclosing account uses physical
    /// field U, rather than treating coefficient work as external supply.
    ///
    /// Reuses the existing whole-volume field bindings, affine source laws,
    /// constant-reference-density solid caloric charts and positive total heat
    /// capacity/coordinate envelope. No field-dependent caloric law, nonlinear
    /// magnetic constitutive law, dielectric relaxation or changing geometry is
    /// inferred. Thermal backward Euler plus endpoint-coefficient circuit
    /// midpoint is generally first order; constitutive convergence is bounded.
    /// Missing contact or graph-resistor heat bindings refuse atomically.
    #[allow(clippy::too_many_arguments)]
    pub fn step_response_field_circuit_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &ConductionResponseCircuit,
        ports: &[f64],
        graph_resistors: &[CellResistor],
        contact_bindings: &[contacts::ConductionContactThermal<'_>],
        fields: &[FieldThermalBinding],
        assembly: WeightedAssemblyLimits,
        field: ConductionControls,
        response: ConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ConductionResponseFieldThermalCandidate<'t>, Error> {
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        let destinations = graph_resistors
            .iter()
            .map(|r| r.destination)
            .collect::<Vec<_>>();
        let next = self
            .step_response_circuit_bound(
                thermal,
                target_s,
                circuit,
                ports,
                &destinations,
                graph_resistors,
                contact_bindings,
                Some(fields),
                assembly,
                field,
                response,
                controls,
                initial_energy_increments_j,
                &mut cancelled,
            )?
            .close_field_report(0.0, controls.total_energy_tolerance_j)?;
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        Ok(next)
    }
}

impl ConductionResponseCircuit {
    /// Initialize source-backed R/C/L and the original nonlinear conductor at
    /// the actual moving operating point. Conserved charge and flux stay fixed;
    /// the native motional EMF is part of the circuit, including at zero voltage.
    ///
    /// Source-storage and fixed-context semantics match new_with_material_storage.
    /// The work cap includes every native field solve and original descriptor
    /// evaluation during source closure.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_motion_and_material_storage(
        problem: PlanarConductionResponse,
        graph: CircuitGraph,
        terminal: ConductionCircuitTerminal,
        mut material_resistors: Vec<MaterialResistor>,
        mut storage: Vec<MaterialStorage>,
        time_s: f64,
        fluxes_wb: Vec<f64>,
        charges_c: Vec<f64>,
        ports: &[f64],
        motion: ConductionMotion<'_>,
        assembly: WeightedAssemblyLimits,
        field: ConductionControls,
        response: ConductionResponseControls,
        maximum_field_solves: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        let graph = lower(&graph, &mut storage)?;
        let graph = stamp_material_resistors(&graph, &terminal, &mut material_resistors)?;
        let mut circuit = Self::new_with_motion(
            problem,
            graph,
            terminal,
            time_s,
            fluxes_wb,
            charges_c,
            ports,
            motion,
            assembly,
            field,
            response,
            maximum_field_solves,
            &mut cancelled,
        )?;
        circuit.circuit.material_resistors = material_resistors;
        circuit.circuit.storage = storage;
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        Ok(circuit)
    }
}

impl ConductionThermalNetwork<'_> {
    /// Close moving nonlinear current, source-backed R(T)/C(T)/L(T) and the
    /// original finite caloric cells on one interval. Only conductor-frame
    /// Joule and stationary electrical-contact loss enter as dissipation.
    /// Signed Lorentz work leaves the electrical/thermal account separately;
    /// field entropy and parameter transfer retain their original owner.
    ///
    /// The original geometry, fixed-domain motion, source support, constant
    /// reference density, affine storage laws and explicit contact heat routes
    /// remain binding. All inner work shares the thermal circuit-trial budget.
    #[allow(clippy::too_many_arguments)]
    pub fn step_moving_response_field_circuit_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &ConductionResponseCircuit,
        ports: &[f64],
        motion: ConductionMotion<'_>,
        graph_resistors: &[CellResistor],
        contact_bindings: &[contacts::ConductionContactThermal<'_>],
        fields: &[FieldThermalBinding],
        assembly: WeightedAssemblyLimits,
        field: ConductionControls,
        response: ConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ConductionResponseFieldThermalCandidate<'t>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: graph_resistors,
            links: self.links,
        };
        let bases = circuit
            .material_resistors()
            .iter()
            .map(|r| (r.branch, r.conductor.source_basis()))
            .collect::<Vec<_>>();
        let prepared = network
            .prepare_with_bases(thermal, &bases, || (*cancel.borrow_mut())())
            .map_err(material)?;
        let destinations = graph_resistors
            .iter()
            .map(|r| r.destination)
            .collect::<Vec<_>>();
        let next = self
            .step_response_circuit_resolved_bound(
                thermal,
                target_s,
                circuit,
                ports,
                &destinations,
                graph_resistors,
                contact_bindings,
                Some(fields),
                Some(motion),
                assembly,
                field,
                response,
                controls,
                initial_energy_increments_j,
                |state| {
                    prepared
                        .resolve(state, || (*cancel.borrow_mut())())
                        .map_err(material)
                },
                |source| source,
                &cancel,
            )?
            .close_field_report(0.0, controls.total_energy_tolerance_j)?;
        check(&cancel)?;
        Ok(next)
    }
}

fn complete_heat_error<E: From<HeatNetworkError>>(
    error: HeatNetworkError,
) -> CircuitEnthalpyError<ConductionHeatError<E>> {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(ConductionHeatError::Heat(other.into())),
    }
}

// Complete heat providers share source/boundary admission and the existing
// circuit endpoint. Frame selection never introduces another thermal loop.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub(in crate::circuit_enthalpy::conduction) fn step_response_circuit_with_heat<
    't,
    A: HeatAssembly + ?Sized,
>(
    heat: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    circuit: &ConductionResponseCircuit,
    ports: &[f64],
    motion: Option<ConductionMotion<'_>>,
    graph_resistors: &[ResistorCell],
    contact_bindings: &[contacts::ConductionContactThermal<'_>],
    fields: Option<&[FieldThermalBinding]>,
    regions: &[ConductionRegion],
    volume_tolerance_m3: f64,
    assembly: WeightedAssemblyLimits,
    field_controls: ConductionControls,
    response: ConductionResponseControls,
    controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
) -> Result<
    ConductionResponseCircuitEnthalpyCandidate<'t, FieldThermalSource<A::Sample>>,
    CircuitEnthalpyError<ConductionHeatError<A::Error>>,
> {
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    check(&cancel).map_err(map_conduction_error)?;
    let prepared = heat.prepared_heat_network();
    let bulk = prepared
        .field_material_network()
        .map_err(complete_heat_error)?;
    if regions.len() != circuit.sources().len() || regions.len() > assembly.max_cells {
        return Err(map_conduction_error(invalid(
            "one bounded heat route per original circuit triangle required",
        )));
    }
    let mut routes = BTreeMap::new();
    let mut field_axes = Vec::new();
    for region in regions {
        check(&cancel).map_err(map_conduction_error)?;
        let source = circuit.sources().get(region.triangle).ok_or_else(|| {
            map_conduction_error(invalid("thermal route names an absent circuit conductor"))
        })?;
        if routes.insert(region.triangle, region.cell).is_some() {
            return Err(map_conduction_error(invalid(
                "original circuit heat routes must be unique",
            )));
        }
        for (axis, _) in source.axes() {
            field_axes.push((region.cell, axis.as_str()));
        }
    }
    if let Some(sources) = circuit.sources.contacts() {
        let trace = circuit
            .circuit
            .problem
            .mesh()
            .boundary_trace()
            .map_err(|_| {
                map_conduction_error(invalid("circuit contact lost its original exterior trace"))
            })?;
        let edge_key = |mut vertices: [usize; 2]| {
            vertices.sort_unstable();
            vertices
        };
        let mut adjacent = BTreeMap::new();
        for edge in &trace.edges {
            check(&cancel).map_err(map_conduction_error)?;
            adjacent.insert(
                edge_key(edge.oriented_vertices.map(|vertex| vertex as usize)),
                edge.source_face,
            );
        }
        for (law, contact) in sources.laws().iter().zip(circuit.contacts()) {
            check(&cancel).map_err(map_conduction_error)?;
            for vertices in [contact.first_vertices, contact.second_vertices] {
                let cell = adjacent
                    .get(&edge_key(vertices))
                    .and_then(|face| routes.get(face))
                    .ok_or_else(|| {
                        map_conduction_error(invalid(
                            "nonlinear circuit contact lacks its original heat route",
                        ))
                    })?;
                for (axis, _) in law.axes() {
                    field_axes.push((*cell, axis.as_str()));
                }
            }
        }
    }
    prepared
        .admit_field_independent_boundaries(cx, &field_axes)
        .map_err(complete_heat_error)?;
    let layout = ConductionThermalNetwork {
        cells: bulk.cells,
        links: bulk.links,
        regions,
        volume_tolerance_m3,
    };
    let mut destinations = prepared.resistor_destinations();
    for route in graph_resistors {
        destinations.push(ResistorCell {
            branch: route.branch,
            cell: route.cell,
        });
    }
    let next = layout.step_response_circuit_resolved_bound(
        thermal,
        target_s,
        circuit,
        ports,
        &destinations,
        bulk.resistors,
        contact_bindings,
        fields,
        motion,
        assembly,
        field_controls,
        response,
        controls,
        initial_energy_increments_j,
        |trial| {
            let coefficients = heat
                .resolve_heat(cx, trial)
                .map_err(|error| error.map_source(ConductionHeatError::Heat))?;
            if !coefficients.contacts.is_empty() {
                return Err(map_conduction_error(invalid(
                    "spatial circuit cannot omit series electrical-contact work",
                )));
            }
            Ok(coefficients.bulk)
        },
        |sample: &A::Sample| &sample.as_ref().material,
        &cancel,
    )?;
    check(&cancel).map_err(map_conduction_error)?;
    Ok(next)
}

/// Advance moving nonlinear conduction and native fixed-C/L circuit storage
/// with any complete heat assembly. The original descriptor receives motional
/// EMF at every actual source iteration, including when the body generates
/// electrical power. Source graph resistance, thermal contacts, layers,
/// radiation, reservoirs and finite links share the same caloric endpoint.
///
/// The full typed heat sample and actual moving field report are retained.
/// Only local conductor-frame Joule and stationary contact loss enter heat;
/// signed Lorentz output is subtracted from independently authored port work.
/// Fixed-domain prescribed motion adds neither induction nor body dynamics.
/// The field-storage entry additionally refreshes source C/L and entropy heat.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn step_moving_response_circuit_candidate<'t, A: HeatAssembly + ?Sized>(
    heat: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    circuit: &ConductionResponseCircuit,
    ports: &[f64],
    motion: ConductionMotion<'_>,
    graph_resistors: &[ResistorCell],
    contact_bindings: &[contacts::ConductionContactThermal<'_>],
    regions: &[ConductionRegion],
    volume_tolerance_m3: f64,
    assembly: WeightedAssemblyLimits,
    field_controls: ConductionControls,
    response: ConductionResponseControls,
    controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
) -> Result<
    ConductionResponseCircuitEnthalpyCandidate<'t, A::Sample>,
    CircuitEnthalpyError<ConductionHeatError<A::Error>>,
> {
    let next = step_response_circuit_with_heat(
        heat,
        cx,
        thermal,
        target_s,
        circuit,
        ports,
        Some(motion),
        graph_resistors,
        contact_bindings,
        None,
        regions,
        volume_tolerance_m3,
        assembly,
        field_controls,
        response,
        controls,
        initial_energy_increments_j,
    )?
    .without_thermal_fields();
    closed_balance(
        target_s - thermal.time_s(),
        next.report.total_energy_residual_j,
        next.report.source.as_ref(),
        controls.total_energy_tolerance_j,
    )
    .map_err(complete_heat_error)?;
    cx.checkpoint()
        .map_err(|_| CircuitEnthalpyError::Cancelled)?;
    Ok(next)
}

/// The same moving circuit/complete heat endpoint with original source-bound
/// R(T)/C(T)/L(T) and physical field energy U = psi + T s. Each storage branch
/// keeps its original whole-volume caloric binding, affine material envelope
/// and independently accounted reversible transfer. External supply is native
/// circuit port work minus signed Lorentz output; mechanical transfer is not
/// deposited as Joule heat or included in source parameter work.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn step_moving_response_field_circuit_candidate<'t, A: HeatAssembly + ?Sized>(
    heat: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    circuit: &ConductionResponseCircuit,
    ports: &[f64],
    motion: ConductionMotion<'_>,
    contact_bindings: &[contacts::ConductionContactThermal<'_>],
    fields: &[FieldThermalBinding],
    regions: &[ConductionRegion],
    volume_tolerance_m3: f64,
    assembly: WeightedAssemblyLimits,
    field_controls: ConductionControls,
    response: ConductionResponseControls,
    controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
) -> Result<
    ConductionResponseFieldThermalCandidate<'t, A::Sample>,
    CircuitEnthalpyError<ConductionHeatError<A::Error>>,
> {
    let next = step_response_circuit_with_heat(
        heat,
        cx,
        thermal,
        target_s,
        circuit,
        ports,
        Some(motion),
        &[],
        contact_bindings,
        Some(fields),
        regions,
        volume_tolerance_m3,
        assembly,
        field_controls,
        response,
        controls,
        initial_energy_increments_j,
    )?;
    let internal_heat_j = finite(
        (target_s - thermal.time_s()) * next.report.source.network.as_ref().radiation_net_power_w(),
    )?;
    let next = next
        .close_field_report(internal_heat_j, controls.total_energy_tolerance_j)
        .map_err(map_conduction_error)?;
    cx.checkpoint()
        .map_err(|_| CircuitEnthalpyError::Cancelled)?;
    Ok(next)
}
