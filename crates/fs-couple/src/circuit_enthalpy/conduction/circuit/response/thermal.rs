//! Source-backed storage and graph resistance inside nonlinear spatial closure.
//!
//! The native DAE, field thermodynamics and finite-heat residual remain their
//! original owners. Each thermal trial refreshes R/C/L and contact states once;
//! every nested bulk constitutive solve sees those same coefficients and starts
//! from the same accepted charge and flux. Reversible field transfer is kept
//! separate from bulk and interface Joule heating.

use super::*;
use crate::circuit_enthalpy::field_storage::{self, FieldThermalCandidate};
use crate::material_circuit::storage::lower;

/// Nonlinear spatial response, source-backed field storage and one finite heat
/// history. The physical energy account uses U = psi + T s at both endpoints.
pub type ConductionResponseFieldThermalCandidate<'t> =
    FieldThermalCandidate<'t, ConductionResponseCircuit, ConductionResponseCircuitStep>;

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
        let (thermal, circuit, coupled) = self
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
            .into_parts();
        let supplied = coupled.electrical.circuit.native.supplied;
        let report =
            field_storage::close_report(coupled, supplied, controls.total_energy_tolerance_j)
                .map_err(super::super::thermal::field_thermal_error)?;
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        Ok(FieldThermalCandidate::new(thermal, circuit, report))
    }
}
