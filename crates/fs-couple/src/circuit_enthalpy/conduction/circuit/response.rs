//! Source-dependent spatial conductors on the original circuit/thermal timeline.
//!
//! Every constitutive trial reuses the original circuit's conserved charge and
//! flux, then evaluates material data at its actual midpoint E/J. Iterating a
//! constitutive law never advances electrical time or deposits trial heat.

/// Source-backed R/C/L and reversible field heat inside the nonlinear solve.
pub mod thermal;

use super::*;
use fs_em::conduction::response::{
    ConductionResponseControls, ConductionResponseError, ConductionResponseReport,
    ConductionResponseSample, MAX_RESPONSE_ITERATIONS, PlanarConductionResponse,
};
use fs_material::conductor::response::ScalarConductorResponse;
use std::sync::Arc;

/// One continuing circuit with an original nonlinear spatial material law.
/// The retained native operator is private: ordinary linear circuit steps must
/// not bypass the response closure. Other graph elements may be explicit ideal
/// data or original source-backed R/C/L. Their thermal refresh uses the shared
/// field/circuit/enthalpy entry point, including reversible field entropy.
pub struct ConductionResponseCircuit {
    circuit: ConductionCircuit,
    sources: Arc<[ScalarConductorResponse]>,
}

/// Actual circuit/field work plus separately reported constitutive defects.
/// Native bulk E.J and exact surface-contact losses are the spatial heat
/// account; neither constitutive defects nor circuit reduction add more heat.
#[derive(Clone, Debug)]
pub struct ConductionResponseCircuitStep {
    /// Native charge/flux storage, port work and local field/graph heat.
    pub circuit: ConductionCircuitStep,
    /// Original source evaluations and defects at the actual midpoint E/J.
    pub evaluated: Vec<ConductionResponseSample>,
    /// Complete native circuit/field trials in this constitutive solve.
    pub field_solves: usize,
    /// Original material evaluations across all constitutive trials.
    pub source_evaluations: usize,
}

/// Admitted joint electrical/material proposal; the accepted input is untouched.
pub struct ConductionResponseCircuitCandidate {
    circuit: ConductionResponseCircuit,
    step: ConductionResponseCircuitStep,
}
impl ConductionResponseCircuitCandidate {
    /// Proposed continuing electrical history and source state.
    pub fn circuit(&self) -> &ConductionResponseCircuit {
        &self.circuit
    }
    /// Original electrical balance and separately retained source defects.
    pub fn step(&self) -> &ConductionResponseCircuitStep {
        &self.step
    }
    /// Publish the history and its balance together.
    pub fn into_parts(self) -> (ConductionResponseCircuit, ConductionResponseCircuitStep) {
        (self.circuit, self.step)
    }
}
impl ThermalDriver for ConductionResponseCircuitCandidate {
    type Model = ConductionResponseCircuit;
    type Step = ConductionResponseCircuitStep;
    fn account(&self) -> DriverAccount<'_> {
        let step = &self.step.circuit;
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: step.heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: step.stored_after_j - step.stored_before_j,
            supplied_j: step.native.supplied,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        self.into_parts()
    }
}

impl ConductionResponseCircuit {
    /// Bind the original source laws to a two-electrode circuit branch and
    /// close its initial electrical operating point. Initial charge/flux remain
    /// exactly authored. Unit-field extraction only lowers each FROZEN trial;
    /// material response is always queried at the actual initial voltage/E/J.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        problem: PlanarConductionResponse,
        graph: CircuitGraph,
        terminal: ConductionCircuitTerminal,
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
        admit_graph(&graph)?;
        if !time_s.is_finite() || time_s < 0.0 {
            return Err(invalid("finite nonnegative initial circuit time required"));
        }
        let (candidate, graph) = problem.solve_coupled_with_materials_candidate(
            problem.conductors().to_vec(),
            assembly,
            field,
            response,
            &mut cancelled,
            |native, cancelled| {
                let resistance = native
                    .resistance(&terminal.positive, &terminal.negative, field, &mut || {
                        cancelled()
                    })
                    .map_err(field_error)?;
                let graph = stamp(&graph, &terminal, resistance.resistance_ohm())?;
                let dae = assemble_circuit(&graph).map_err(electrical)?;
                let x = consistent(&dae, &fluxes_wb, &charges_c, ports)?;
                let voltage = |node| {
                    if node == 0 {
                        0.0
                    } else {
                        x[dae.node_potential_index[node - 1]]
                    }
                };
                let (a, b, _) = graph.branches[terminal.branch];
                let field = resistance
                    .field_at_voltage(
                        finite(voltage(a) - voltage(b))?,
                        field.current_tolerance_a,
                        field.power_tolerance_w,
                        &mut || cancelled(),
                    )
                    .map_err(field_error)?;
                Ok::<_, Error>((field, graph))
            },
        )?;
        let (problem, _) = candidate.into_parts();
        let (problem, sources) = problem.into_frozen_parts();
        Ok(Self {
            circuit: ConductionCircuit {
                problem,
                graph,
                terminal,
                fluxes_wb,
                charges_c,
                time_s,
                material_resistors: Vec::new(),
                storage: Vec::new(),
            },
            sources,
        })
    }

    /// Graph snapshot at the accepted operating point, for inspection.
    /// Future steps must reclose the spatial source response.
    pub fn graph(&self) -> &CircuitGraph {
        self.circuit.graph()
    }
    /// Original graph branch and equipotential exterior electrodes.
    pub fn terminal(&self) -> &ConductionCircuitTerminal {
        self.circuit.terminal()
    }
    /// Complete source coefficients used by the accepted spatial operator.
    pub fn conductors(&self) -> &[ResolvedConductorState] {
        self.circuit.problem.conductors()
    }
    /// Original complete bulk laws; changing circuit amplitude never replaces
    /// their selected source claims or their explicitly owned E/J coordinates.
    pub fn sources(&self) -> &[ScalarConductorResponse] {
        &self.sources
    }
    /// Original source-bound contacts in canonical native trace order.
    pub fn contacts(&self) -> &[fs_em::conduction::ConductionContact] {
        self.circuit.problem.contacts()
    }
    /// Conserved magnetic coordinates in original inductor branch order.
    pub fn fluxes_wb(&self) -> &[f64] {
        self.circuit.fluxes_wb()
    }
    /// Conserved electric coordinates in original capacitor branch order.
    pub fn charges_c(&self) -> &[f64] {
        self.circuit.charges_c()
    }
    /// Accepted physical time; source iterations never alter it.
    pub fn time_s(&self) -> f64 {
        self.circuit.time_s()
    }

    fn response(
        &self,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<PlanarConductionResponse, Error> {
        PlanarConductionResponse::from_frozen_parts(
            self.circuit.problem.clone(),
            self.sources.clone(),
            &mut || (*cancel.borrow_mut())(),
        )
        .map_err(Error::from)
    }

    /// Propose one native midpoint circuit interval at fixed thermal context.
    /// All source-dependent spatial resistance is solved inside that SAME
    /// interval from the original charge/flux. Local Joule heat is reported,
    /// not deposited; use `step_response_circuit_candidate` for finite heating.
    /// Contact-bearing circuits require `step_response_contact_circuit_candidate`
    /// so their original source refresh and heat destinations remain explicit.
    /// This Picard closure is bounded and may refuse instead of converging.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate(
        &self,
        target_s: f64,
        ports: &[f64],
        assembly: WeightedAssemblyLimits,
        field: ConductionControls,
        response: ConductionResponseControls,
        circuit: MaterialCircuitSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ConductionResponseCircuitCandidate, Error> {
        let cancel = RefCell::new(cancelled);
        let problem = self.response(&cancel)?;
        self.trial(
            &self.circuit,
            &problem,
            problem.conductors().to_vec(),
            problem.contacts().to_vec(),
            self.circuit.material_resistors.clone(),
            target_s,
            ports,
            assembly,
            field,
            response,
            circuit,
            None,
            &cancel,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn trial(
        &self,
        driver: &ConductionCircuit,
        problem: &PlanarConductionResponse,
        sources: Vec<ResolvedConductorState>,
        contact_sources: Vec<fs_em::conduction::ConductionContact>,
        material_resistors: Vec<MaterialResistor>,
        target_s: f64,
        ports: &[f64],
        assembly: WeightedAssemblyLimits,
        field: ConductionControls,
        response: ConductionResponseControls,
        controls: MaterialCircuitSolve,
        contact_heat: Option<&contacts::PreparedContacts<'_, '_>>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<ConductionResponseCircuitCandidate, Error> {
        check(cancel)?;
        if !target_s.is_finite()
            || target_s <= self.time_s()
            || !controls.energy_tolerance_j.is_finite()
            || controls.energy_tolerance_j <= 0.0
            || !controls.power_tolerance_w.is_finite()
            || controls.power_tolerance_w <= 0.0
        {
            return Err(invalid(
                "positive circuit interval and finite positive balance allowances required",
            ));
        }
        let (candidate, trial) = problem.solve_coupled_with_materials_and_contacts_candidate(
            sources,
            contact_sources,
            assembly,
            field,
            response,
            &mut || (*cancel.borrow_mut())(),
            |native, _| {
                driver.trial_field(
                    native,
                    material_resistors.clone(),
                    target_s,
                    ports,
                    field,
                    controls,
                    contact_heat,
                    cancel,
                )
            },
        )?;
        let (problem, report) = candidate.into_parts();
        let (problem, sources) = problem.into_frozen_parts();
        let native = trial.finish(driver, problem);
        let ConductionResponseReport {
            evaluated,
            field_solves,
            source_evaluations,
            ..
        } = report;
        Ok(ConductionResponseCircuitCandidate {
            circuit: Self {
                circuit: native.circuit,
                sources,
            },
            step: ConductionResponseCircuitStep {
                circuit: native.step,
                evaluated,
                field_solves,
                source_evaluations,
            },
        })
    }
}

/// A single proposal for the nonlinear spatial law, electrical storage, local
/// thermal material state and finite enthalpy. Refusal publishes none of them.
pub struct ConductionResponseCircuitEnthalpyCandidate<'t, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    circuit: ConductionResponseCircuit,
    report: CircuitEnthalpyReport<S, ConductionResponseCircuitStep>,
}
impl<'t, S> ConductionResponseCircuitEnthalpyCandidate<'t, S> {
    /// Proposed original finite-enthalpy state.
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    /// Proposed native circuit history and complete material sources.
    pub fn circuit(&self) -> &ConductionResponseCircuit {
        &self.circuit
    }
    /// Independent native electrical, thermal and combined balance accounts.
    pub fn report(&self) -> &CircuitEnthalpyReport<S, ConductionResponseCircuitStep> {
        &self.report
    }
    /// Publish all coupled state with its original balance report.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        ConductionResponseCircuit,
        CircuitEnthalpyReport<S, ConductionResponseCircuitStep>,
    ) {
        (self.thermal, self.circuit, self.report)
    }
}

impl ConductionThermalNetwork<'_> {
    /// Close the native circuit midpoint, original sigma/rho(T, phase, |E|, |J|)
    /// and finite backward-Euler heat at one common endpoint. Actual E.J heat
    /// enters each original full-volume thermal cell once. Every other ideal
    /// graph resistor needs an explicit receiving cell. All discarded response
    /// trials count against the common circuit work budget.
    ///
    /// The thermal/caloric source restrictions are the existing nonlinear
    /// conductor restrictions: no silently frozen field-dependent Cp or k, and
    /// no reversible field caloric work. Graph R/L/C, geometry and topology are
    /// fixed; spatial sources and their original claims remain pinned. Native
    /// circuit work, constitutive defects and combined energy remain separate.
    #[allow(clippy::too_many_arguments)]
    pub fn step_response_circuit_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &ConductionResponseCircuit,
        ports: &[f64],
        graph_resistors: &[ResistorCell],
        assembly: WeightedAssemblyLimits,
        field: ConductionControls,
        response: ConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ConductionResponseCircuitEnthalpyCandidate<'t>, Error> {
        self.step_response_contact_circuit_candidate(
            thermal,
            target_s,
            circuit,
            ports,
            graph_resistors,
            &[],
            assembly,
            field,
            response,
            controls,
            initial_energy_increments_j,
            cancelled,
        )
    }

    /// Solve native circuit storage, nonlinear bulk E/J response, original
    /// contact resistance and finite two-sided heating at one endpoint. All
    /// constitutive trials restart from the same conserved charge/flux and
    /// include the native contact operator. Only accepted physical bulk/surface
    /// losses enter the original enthalpy residual, under explicit heat splits.
    #[allow(clippy::too_many_arguments)]
    pub fn step_response_contact_circuit_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &ConductionResponseCircuit,
        ports: &[f64],
        graph_resistors: &[ResistorCell],
        contact_bindings: &[contacts::ConductionContactThermal<'_>],
        assembly: WeightedAssemblyLimits,
        field: ConductionControls,
        response: ConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ConductionResponseCircuitEnthalpyCandidate<'t>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        self.step_response_circuit_resolved_candidate(
            thermal,
            target_s,
            circuit,
            ports,
            graph_resistors,
            contact_bindings,
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
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_response_circuit_bound<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &ConductionResponseCircuit,
        ports: &[f64],
        graph_resistors: &[ResistorCell],
        material_resistors: &[CellResistor],
        contact_bindings: &[contacts::ConductionContactThermal<'_>],
        field_bindings: Option<&[FieldThermalBinding]>,
        assembly: WeightedAssemblyLimits,
        field: ConductionControls,
        response: ConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ConductionResponseCircuitEnthalpyCandidate<'t, FieldThermalSource>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: material_resistors,
            links: self.links,
        };
        let bases = circuit
            .circuit
            .material_resistors
            .iter()
            .map(|r| (r.branch, r.conductor.source_basis()))
            .collect::<Vec<_>>();
        let prepared = network
            .prepare_with_bases(thermal, &bases, || (*cancel.borrow_mut())())
            .map_err(material)?;
        self.step_response_circuit_resolved_bound(
            thermal,
            target_s,
            circuit,
            ports,
            graph_resistors,
            material_resistors,
            contact_bindings,
            field_bindings,
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
        )
    }

    // Original plain and complete heat networks enter this SAME continuing
    // circuit/material/enthalpy residual. Source observations remain typed and
    // retained; no contact or radiation owner is flattened into fitted heat.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::circuit_enthalpy::conduction) fn step_response_circuit_resolved_candidate<
        't,
        S,
    >(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &ConductionResponseCircuit,
        ports: &[f64],
        graph_resistors: &[ResistorCell],
        contact_bindings: &[contacts::ConductionContactThermal<'_>],
        assembly: WeightedAssemblyLimits,
        field: ConductionControls,
        response: ConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        resolve_sources: impl Fn(
            &DepositionSimulation<'_>,
        ) -> Result<CircuitThermalCoefficients<S>, Error>,
        sample: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<ConductionResponseCircuitEnthalpyCandidate<'t, S>, Error> {
        self.step_response_circuit_resolved_bound(
            thermal,
            target_s,
            circuit,
            ports,
            graph_resistors,
            &[],
            contact_bindings,
            None,
            assembly,
            field,
            response,
            controls,
            initial_energy_increments_j,
            resolve_sources,
            sample,
            cancel,
        )
        .map(ConductionResponseCircuitEnthalpyCandidate::without_thermal_fields)
    }

    // One constitutive/thermal residual retains the original typed heat
    // observations, including callers that resolve contacts and radiation.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::circuit_enthalpy::conduction) fn step_response_circuit_resolved_bound<
        't,
        S,
        E: From<ConductionThermalError>,
    >(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &ConductionResponseCircuit,
        ports: &[f64],
        graph_resistors: &[ResistorCell],
        material_resistors: &[CellResistor],
        contact_bindings: &[contacts::ConductionContactThermal<'_>],
        field_bindings: Option<&[FieldThermalBinding]>,
        assembly: WeightedAssemblyLimits,
        field: ConductionControls,
        response: ConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        resolve_sources: impl Fn(
            &DepositionSimulation<'_>,
        )
            -> Result<CircuitThermalCoefficients<S>, CircuitEnthalpyError<E>>,
        sample: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<
        ConductionResponseCircuitEnthalpyCandidate<'t, FieldThermalSource<S>>,
        CircuitEnthalpyError<E>,
    > {
        check(cancel).map_err(map_conduction_error)?;
        if !(1..=MAX_RESPONSE_ITERATIONS).contains(&response.maximum_iterations)
            || circuit.time_s().to_bits() != thermal.time_s().to_bits()
            || material_resistors.len() != circuit.circuit.material_resistors.len()
        {
            return Err(map_conduction_error(invalid(
                "matching clocks, bounded response and complete material graph-resistor bindings required",
            )));
        }
        let problem = circuit.response(cancel).map_err(map_conduction_error)?;
        let dt = self
            .admit_mesh_step(
                thermal,
                target_s,
                problem.mesh(),
                assembly,
                controls,
                initial_energy_increments_j,
            )
            .map_err(map_conduction_error)?;
        let routes = self
            .admit_mesh_routes(thermal, problem.mesh(), cancel)
            .map_err(map_conduction_error)?;
        let graph_routes = graph_heat_routes(thermal, &circuit.circuit, graph_resistors, cancel)
            .map_err(map_conduction_error)?;
        self.admit_response_caloric(&problem, cancel)
            .map_err(map_conduction_error)?;
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: material_resistors,
            links: self.links,
        };
        for binding in material_resistors {
            check(cancel).map_err(map_conduction_error)?;
            if routes.contains(&binding.destination.cell) {
                return Err(map_conduction_error(invalid(
                    "a whole-cell material resistor cannot reuse the spatial conductor's volume",
                )));
            }
        }
        let resolve = |state: &DepositionSimulation<'_>| {
            let coefficients = resolve_sources(state)?;
            self.admit_response_heat(sample(&coefficients.source), cancel)
                .map_err(map_conduction_error)?;
            Ok::<_, CircuitEnthalpyError<E>>(coefficients)
        };
        let accepted = resolve(thermal)?;
        for original in &circuit.circuit.material_resistors {
            check(cancel).map_err(map_conduction_error)?;
            let next = accepted
                .resistors
                .iter()
                .find(|r| r.branch == original.branch)
                .ok_or_else(|| invalid("accepted material-resistor branch is absent"))
                .map_err(map_conduction_error)?;
            if !same_conductor(&original.conductor, &next.conductor) {
                return Err(map_conduction_error(invalid(
                    "accepted graph-resistor geometry and sources must match the actual thermal state",
                )));
            }
        }
        let contacts = contacts::PreparedContacts::new_for_contacts(
            self,
            problem.mesh(),
            problem.contacts(),
            problem.contact_liquid_mass_fraction_axes(),
            contact_bindings,
            &routes,
            sample(&accepted.source),
            cancel,
        )
        .map_err(map_conduction_error)?;
        let initial = self
            .response_materials_at(&problem, sample(&accepted.source), &routes, cancel)
            .map_err(map_conduction_error)?;
        if initial.iter().zip(problem.conductors()).any(|(a, b)| {
            !a.has_same_electrical_source(b) || !same_material(a.material(), b.material())
        }) {
            return Err(map_conduction_error(invalid(
                "accepted response-circuit sources must match the actual thermal state",
            )));
        }
        let fields = field_bindings
            .map(|bindings| {
                PreparedFieldThermal::new(
                    &network,
                    thermal,
                    circuit.circuit.field_view(),
                    bindings,
                    cancel,
                )
                .map_err(super::thermal::field_thermal_error)
            })
            .transpose()
            .map_err(map_conduction_error)?;
        let (thermal, circuit, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |state, remaining| {
                let mut coefficients = resolve(state)?;
                let sources = self
                    .response_materials_at(&problem, sample(&coefficients.source), &routes, cancel)
                    .map_err(map_conduction_error)?;
                let contact_sources = contacts
                    .refresh_contacts(problem.contacts(), sample(&coefficients.source), cancel)
                    .map_err(map_conduction_error)?;
                let bounded = ConductionResponseControls {
                    maximum_iterations: response.maximum_iterations.min(remaining),
                    ..response
                };
                let refreshed = fields
                    .as_ref()
                    .map(|fields| {
                        let (storage, samples) = fields
                            .refresh(state, cancel)
                            .map_err(super::thermal::field_thermal_error)?;
                        let rebound =
                            circuit
                                .circuit
                                .rebind_storage_candidate(storage, ports, || {
                                    (*cancel.borrow_mut())()
                                })?;
                        Ok::<_, Error>((rebound, samples))
                    })
                    .transpose()
                    .map_err(map_conduction_error)?;
                let driver = refreshed
                    .as_ref()
                    .map_or(&circuit.circuit, |(rebound, _)| rebound.circuit());
                let candidate = circuit
                    .trial(
                        driver,
                        &problem,
                        sources,
                        contact_sources,
                        coefficients.resistors,
                        target_s,
                        ports,
                        assembly,
                        field,
                        bounded,
                        controls.circuit,
                        Some(&contacts),
                        cancel,
                    )
                    .map_err(|error| match error {
                        CircuitEnthalpyError::Source(ConductionThermalError::Response(
                            ConductionResponseError::NotConverged { .. },
                        )) if bounded.maximum_iterations < response.maximum_iterations => {
                            CircuitEnthalpyError::CircuitBudget
                        }
                        other => other,
                    })
                    .map_err(map_conduction_error)?;
                let circuit_trials = candidate.step.field_solves;
                let source = if let (Some(fields), Some((rebound, samples))) = (&fields, refreshed)
                {
                    fields
                        .account(
                            samples,
                            candidate.circuit.circuit.field_view(),
                            rebound.parameter_work_j(),
                            dt,
                            coefficients.source,
                            &mut coefficients.sources_w,
                            controls.total_energy_tolerance_j,
                            cancel,
                        )
                        .map_err(super::thermal::field_thermal_error)
                        .map_err(map_conduction_error)?
                } else {
                    FieldThermalSource {
                        network: coefficients.source,
                        fields: Vec::new(),
                        field_internal_energy_j: [0.0; 2],
                        reversible_heat_j: 0.0,
                    }
                };
                Ok(DrivenCoefficients {
                    circuit: candidate,
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source,
                    circuit_trials,
                })
            },
            |candidate, _dt, checkpoint| {
                (|| {
                    let mut powers = self.deposit_power(
                        &candidate.step.circuit.spatial.cell_power_w,
                        &routes,
                        checkpoint,
                    )?;
                    for heat in &candidate.step.circuit.graph_resistors {
                        checkpoint()?;
                        let cell = graph_routes[&heat.branch];
                        let previous = *powers.get(&cell).unwrap_or(&0.0);
                        let total = finite(previous + heat.power_w)?;
                        if heat.power_w != 0.0 && total == previous {
                            return Err(invalid(
                                "response-circuit heat disappeared in receiving-cell sum",
                            ));
                        }
                        powers.insert(cell, total);
                    }
                    contacts::deposit(
                        &candidate.step.circuit.spatial.contacts,
                        &mut powers,
                        checkpoint,
                    )?;
                    Ok::<_, Error>(powers)
                })()
                .map_err(map_conduction_error)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(cancel).map_err(map_conduction_error)?;
        Ok(ConductionResponseCircuitEnthalpyCandidate {
            thermal,
            circuit,
            report,
        })
    }
}
