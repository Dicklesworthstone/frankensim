//! Nonlinear tetrahedral resistance on the original circuit and heat timeline.
//!
//! Every source trial extracts the full native terminal matrix and advances
//! the SAME accepted charge/flux once. Only the admitted proposal contributes
//! local bulk, contact and material-branch heat to finite enthalpy.

use super::*;
use fs_em::{
    conduction::response::{ConductionResponseSources, MAX_RESPONSE_ITERATIONS},
    spatial::conduction::response::{
        SpatialConductionResponse, SpatialConductionResponseControls,
        SpatialConductionResponseError, SpatialConductionResponseSample,
    },
};
use fs_material::conductor::response::ScalarConductorResponse;

impl From<SpatialConductionResponseError> for Error {
    fn from(error: SpatialConductionResponseError) -> Self {
        match error {
            SpatialConductionResponseError::Cancelled => Self::Cancelled,
            other => Self::Source(ConductionThermalError::SpatialResponse(other)),
        }
    }
}

/// Continuing original bulk/contact laws and conserved electrical coordinates.
/// The frozen native circuit stays private so future steps cannot omit source
/// closure. C/L, mutual inductance, ideal transformers and geometry are fixed;
/// their temperature-dependent parameter work needs its own reversible owner.
pub struct SpatialConductionResponseCircuit {
    circuit: SpatialConductionCircuit,
    sources: ConductionResponseSources,
}

#[derive(Clone, Debug)]
pub struct SpatialConductionResponseCircuitStep {
    /// Native midpoint storage, source work and all local physical heat.
    pub circuit: SpatialConductionCircuitStep,
    pub evaluated: Vec<SpatialConductionResponseSample>,
    pub evaluated_contacts: Vec<fs_em::conduction::response::contact::ContactResponseSample>,
    pub iterations: usize,
    /// N-1 native basis solves for EVERY accepted or discarded source trial.
    pub field_solves: usize,
    pub source_evaluations: usize,
}

pub struct SpatialConductionResponseCircuitCandidate {
    circuit: SpatialConductionResponseCircuit,
    step: SpatialConductionResponseCircuitStep,
}
impl SpatialConductionResponseCircuitCandidate {
    pub fn circuit(&self) -> &SpatialConductionResponseCircuit {
        &self.circuit
    }
    pub fn step(&self) -> &SpatialConductionResponseCircuitStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        SpatialConductionResponseCircuit,
        SpatialConductionResponseCircuitStep,
    ) {
        (self.circuit, self.step)
    }
}
impl ThermalDriver for SpatialConductionResponseCircuitCandidate {
    type Model = SpatialConductionResponseCircuit;
    type Step = SpatialConductionResponseCircuitStep;
    fn account(&self) -> DriverAccount<'_> {
        let step = &self.step.circuit;
        DriverAccount {
            resistors: &step.resistors,
            contacts: &[],
            heat_j: step.heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: step.stored_change_j,
            supplied_j: step.supplied_work_j,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        self.into_parts()
    }
}

impl SpatialConductionResponseCircuit {
    /// Close the initial actual operating point without changing authored C/L
    /// coordinates. All electrodes retain their full native matrix and order.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        problem: SpatialConductionResponse,
        graph: CircuitGraph,
        terminals: Vec<Vec<usize>>,
        nodes: Vec<usize>,
        time_s: f64,
        fluxes_wb: Vec<f64>,
        charges_c: Vec<f64>,
        ports: &[f64],
        matrix: SpatialConductanceControls,
        assembly: WeightedAssemblyLimits,
        field: ConductionControls,
        response: SpatialConductionResponseControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        Self::new_with_couplings(
            problem,
            graph,
            terminals,
            nodes,
            Vec::new(),
            Vec::new(),
            time_s,
            fluxes_wb,
            charges_c,
            ports,
            matrix,
            assembly,
            field,
            response,
            cancelled,
        )
    }

    /// Original material graph resistors and reciprocal mutual links remain
    /// explicit. Their ordinary source admission and complete branch coverage
    /// are the same as the existing native 3D circuit owner.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_couplings(
        problem: SpatialConductionResponse,
        graph: CircuitGraph,
        terminals: Vec<Vec<usize>>,
        nodes: Vec<usize>,
        material_resistors: Vec<MaterialResistor>,
        mutual: Vec<MutualInductance>,
        time_s: f64,
        fluxes_wb: Vec<f64>,
        charges_c: Vec<f64>,
        ports: &[f64],
        matrix: SpatialConductanceControls,
        assembly: WeightedAssemblyLimits,
        field: ConductionControls,
        response: SpatialConductionResponseControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        if !time_s.is_finite() || time_s < 0.0 {
            return Err(invalid("finite nonnegative initial circuit time required"));
        }
        if nodes.len() != terminals.len()
            || nodes.len() < 2
            || nodes.len() > graph.node_count
            || nodes.iter().any(|&node| node >= graph.node_count)
            || nodes.iter().copied().collect::<BTreeSet<_>>().len() != nodes.len()
        {
            return Err(invalid(
                "complete unique original terminal-to-circuit-node mapping required",
            ));
        }
        let (graph, material_resistors) =
            lower_resistors(&graph, material_resistors).map_err(material_circuit)?;
        let faces = terminals.iter().map(Vec::as_slice).collect::<Vec<_>>();
        let (candidate, terminals) = problem.solve_coupled_terminal_candidate(
            problem.conductors().to_vec(),
            problem.contacts().to_vec(),
            &faces,
            matrix,
            assembly,
            field,
            response,
            &mut cancelled,
            |_native, matrix, cancelled| {
                let cancel = RefCell::new(cancelled);
                let (_, dae) = prepare_matrix(matrix, &graph, &mutual, &nodes, &cancel)?;
                let initial = consistent(&dae, &fluxes_wb, &charges_c, ports)?;
                let voltages = nodes
                    .iter()
                    .map(|&node| {
                        if node == 0 {
                            0.0
                        } else {
                            initial[dae.node_potential_index[node - 1]]
                        }
                    })
                    .collect::<Vec<_>>();
                Ok::<_, Error>((voltages, matrix.terminal_faces().to_vec()))
            },
        )?;
        let (problem, _) = candidate.into_parts();
        let (problem, sources) = problem.into_frozen_parts();
        Ok(Self {
            circuit: SpatialConductionCircuit {
                problem,
                graph,
                mutual,
                terminals,
                nodes,
                fluxes_wb,
                charges_c,
                time_s,
                material_resistors,
            },
            sources,
        })
    }

    pub fn graph(&self) -> &CircuitGraph {
        self.circuit.graph()
    }
    pub fn terminals(&self) -> &[Vec<usize>] {
        &self.circuit.terminals
    }
    pub fn terminal_nodes(&self) -> &[usize] {
        &self.circuit.nodes
    }
    pub fn conductors(&self) -> &[ResolvedConductorState] {
        self.circuit.problem.conductors()
    }
    pub fn sources(&self) -> &[ScalarConductorResponse] {
        self.sources.bulk()
    }
    pub fn contacts(&self) -> &[SpatialConductionContact] {
        self.circuit.problem.contacts()
    }
    pub fn fluxes_wb(&self) -> &[f64] {
        &self.circuit.fluxes_wb
    }
    pub fn charges_c(&self) -> &[f64] {
        &self.circuit.charges_c
    }
    pub fn time_s(&self) -> f64 {
        self.circuit.time_s
    }
    pub fn relaxation(&self) -> f64 {
        self.sources.relaxation()
    }
    pub fn adaptive_relaxation(
        &self,
    ) -> Option<fs_em::conduction::response::AdaptiveConductionRelaxation> {
        self.sources.adaptive_relaxation()
    }
    pub fn material_resistors(&self) -> &[MaterialResistor] {
        &self.circuit.material_resistors
    }

    fn response(
        &self,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<SpatialConductionResponse, Error> {
        SpatialConductionResponse::from_frozen_parts(
            self.circuit.problem.clone(),
            self.sources.clone(),
            &mut || (*cancel.borrow_mut())(),
        )
        .map_err(Into::into)
    }

    /// Propose one interval at fixed external thermal context. Every source
    /// iteration restarts from accepted C/L; rejected trials deposit no heat.
    /// Contact-bearing circuits require the finite-enthalpy entry point so both
    /// source temperatures and heat receivers remain explicitly bound.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate(
        &self,
        target_s: f64,
        ports: &[f64],
        matrix: SpatialConductanceControls,
        assembly: WeightedAssemblyLimits,
        field: ConductionControls,
        response: SpatialConductionResponseControls,
        controls: MaterialCircuitSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialConductionResponseCircuitCandidate, Error> {
        let cancel = RefCell::new(cancelled);
        let problem = self.response(&cancel)?;
        self.trial(
            &problem,
            problem.conductors().to_vec(),
            problem.contacts().to_vec(),
            self.circuit.material_resistors.clone(),
            target_s,
            ports,
            matrix,
            assembly,
            field,
            response,
            controls,
            None,
            usize::MAX,
            &cancel,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn trial(
        &self,
        problem: &SpatialConductionResponse,
        seeds: Vec<ResolvedConductorState>,
        contacts: Vec<SpatialConductionContact>,
        material_resistors: Vec<MaterialResistor>,
        target_s: f64,
        ports: &[f64],
        matrix: SpatialConductanceControls,
        assembly: WeightedAssemblyLimits,
        field: ConductionControls,
        response: SpatialConductionResponseControls,
        controls: MaterialCircuitSolve,
        contact_heat: Option<&contacts::PreparedContacts<'_, '_>>,
        remaining: usize,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<SpatialConductionResponseCircuitCandidate, Error> {
        check(cancel)?;
        if finite(target_s - self.time_s())? <= 0.0
            || !controls.energy_tolerance_j.is_finite()
            || controls.energy_tolerance_j <= 0.0
            || !controls.power_tolerance_w.is_finite()
            || controls.power_tolerance_w <= 0.0
        {
            return Err(invalid(
                "positive circuit interval and finite positive work tolerances required",
            ));
        }
        let cost = self.circuit.terminals.len();
        if cost < 2 {
            return Err(invalid("at least two circuit electrodes required"));
        }
        if remaining < cost {
            return Err(CircuitEnthalpyError::CircuitBudget);
        }
        // Each source iteration costs N-1 field solves PLUS one descriptor
        // step. Reserve that entire cost before any part of the next trial.
        let bounded = SpatialConductionResponseControls {
            maximum_field_solves: response
                .maximum_field_solves
                .min((remaining / cost) * (cost - 1)),
            ..response
        };
        let faces = self
            .circuit
            .terminals
            .iter()
            .map(Vec::as_slice)
            .collect::<Vec<_>>();
        let (candidate, stage) = problem
            .solve_coupled_terminal_candidate(
                seeds,
                contacts,
                &faces,
                matrix,
                assembly,
                field,
                bounded,
                &mut || (*cancel.borrow_mut())(),
                |native, matrix, cancelled| {
                    let cancel = RefCell::new(cancelled);
                    let stage = self.circuit.trial_with_matrix(
                        native,
                        matrix,
                        material_resistors.clone(),
                        target_s,
                        ports,
                        field,
                        controls,
                        contact_heat,
                        &cancel,
                    )?;
                    Ok::<_, Error>((stage.step.terminal_voltage_v.clone(), stage))
                },
            )
            .map_err(|error| match error {
                CircuitEnthalpyError::Source(ConductionThermalError::SpatialResponse(
                    SpatialConductionResponseError::FieldSolveBudget { .. },
                )) if bounded.maximum_field_solves < response.maximum_field_solves => {
                    CircuitEnthalpyError::CircuitBudget
                }
                other => other,
            })?;
        let (problem, report) = candidate.into_parts();
        let (problem, sources) = problem.into_frozen_parts();
        let next = stage.finish(&self.circuit, problem);
        Ok(SpatialConductionResponseCircuitCandidate {
            circuit: Self {
                circuit: next.circuit,
                sources,
            },
            step: SpatialConductionResponseCircuitStep {
                circuit: next.step,
                evaluated: report.evaluated,
                evaluated_contacts: report.evaluated_contacts,
                iterations: report.iterations,
                field_solves: report.field_solves,
                source_evaluations: report.source_evaluations,
            },
        })
    }
}

pub struct SpatialConductionResponseCircuitEnthalpyCandidate<'t> {
    thermal: DepositionSimulation<'t>,
    circuit: SpatialConductionResponseCircuit,
    report: CircuitEnthalpyReport<MaterialNetworkSample, SpatialConductionResponseCircuitStep>,
}
impl<'t> SpatialConductionResponseCircuitEnthalpyCandidate<'t> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn circuit(&self) -> &SpatialConductionResponseCircuit {
        &self.circuit
    }
    pub fn report(
        &self,
    ) -> &CircuitEnthalpyReport<MaterialNetworkSample, SpatialConductionResponseCircuitStep> {
        &self.report
    }
    /// Publish source coefficients, electrical history and finite heat jointly.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        SpatialConductionResponseCircuit,
        CircuitEnthalpyReport<MaterialNetworkSample, SpatialConductionResponseCircuitStep>,
    ) {
        (self.thermal, self.circuit, self.report)
    }
}

impl SpatialConductionThermalNetwork<'_> {
    /// Close actual midpoint E/J, original source-dependent resistance, native
    /// C/L storage and finite backward-Euler heat on one accepted endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn step_response_circuit_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &SpatialConductionResponseCircuit,
        ports: &[f64],
        contact_bindings: &[ConductionContactThermal<'_>],
        matrix: SpatialConductanceControls,
        assembly: WeightedAssemblyLimits,
        field: ConductionControls,
        response: SpatialConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialConductionResponseCircuitEnthalpyCandidate<'t>, Error> {
        self.step_response_material_circuit_candidate(
            thermal,
            target_s,
            circuit,
            ports,
            &[],
            contact_bindings,
            matrix,
            assembly,
            field,
            response,
            controls,
            initial_energy_increments_j,
            cancelled,
        )
    }

    /// Include original whole-cell graph resistors and explicitly partitioned
    /// surface-contact losses. Each native basis solve AND descriptor step uses
    /// the common work allowance; discarded source/thermal trials publish no
    /// history and deposit no trial heat. Caloric and heat-transport sources
    /// must remain independent of E/J, as in the existing 3D response owner.
    #[allow(clippy::too_many_arguments)]
    pub fn step_response_material_circuit_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &SpatialConductionResponseCircuit,
        ports: &[f64],
        material_resistors: &[CellResistor],
        contact_bindings: &[ConductionContactThermal<'_>],
        matrix: SpatialConductanceControls,
        assembly: WeightedAssemblyLimits,
        field: ConductionControls,
        response: SpatialConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialConductionResponseCircuitEnthalpyCandidate<'t>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        if !(1..=MAX_RESPONSE_ITERATIONS).contains(&response.maximum_iterations)
            || response.maximum_field_solves == 0
        {
            return Err(invalid("invalid bounded 3D circuit source response limits"));
        }
        self.admit_circuit_step(
            thermal,
            target_s,
            &circuit.circuit,
            material_resistors.len(),
            assembly,
            controls,
            initial_energy_increments_j,
            &cancel,
        )?;
        let problem = circuit.response(&cancel)?;
        let regions = self
            .regions
            .iter()
            .map(|region| ConductionRegion {
                triangle: region.tetrahedron,
                cell: region.cell,
            })
            .collect::<Vec<_>>();
        let route_owner = ConductionThermalNetwork {
            cells: self.cells,
            links: self.links,
            regions: &regions,
            volume_tolerance_m3: self.volume_tolerance_m3,
        };
        let volumes = (0..problem.mesh().tets.len())
            .map(|tet| {
                check(&cancel)?;
                problem
                    .cell_volume_m3(tet)
                    .ok_or_else(|| invalid("original tetrahedral volume is absent"))
            })
            .collect::<Result<Vec<_>, Error>>()?;
        let routes = route_owner.admit_volume_routes(thermal, &volumes, &cancel)?;
        route_owner.admit_response_caloric_sources(problem.sources(), &cancel)?;
        let mut resistor_routes = BTreeMap::new();
        let mut resistor_receivers = BTreeSet::new();
        for binding in material_resistors {
            check(&cancel)?;
            let destination = binding.destination;
            if !thermal.cells().contains_key(&destination.cell)
                || routes.contains(&destination.cell)
                || !resistor_receivers.insert(destination.cell)
                || resistor_routes
                    .insert(destination.branch, destination.cell)
                    .is_some()
            {
                return Err(invalid(
                    "each whole-cell resistor needs a unique active receiver outside the spatial field volume",
                ));
            }
        }
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: material_resistors,
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
        let resolve = |state: &DepositionSimulation<'t>| {
            let coefficients = prepared
                .resolve(state, || (*cancel.borrow_mut())())
                .map_err(material)?;
            route_owner.admit_response_heat(&coefficients.source, &cancel)?;
            let branches = coefficients
                .resistors
                .iter()
                .map(|r| r.branch)
                .collect::<BTreeSet<_>>();
            if branches.len() != coefficients.resistors.len()
                || branches.len() != resistor_routes.len()
                || branches
                    .iter()
                    .any(|branch| !resistor_routes.contains_key(branch))
            {
                return Err(invalid(
                    "original graph resistor sources and receivers must cover the same branches",
                ));
            }
            Ok::<_, Error>(coefficients)
        };
        let refresh = |sample: &MaterialNetworkSample| {
            route_owner.response_materials_from(
                problem.conductors(),
                problem.sources(),
                sample,
                &routes,
                &cancel,
            )
        };
        let accepted = resolve(thermal)?;
        for original in circuit.material_resistors() {
            let next = accepted
                .resistors
                .iter()
                .find(|r| r.branch == original.branch)
                .ok_or_else(|| invalid("accepted material-resistor branch is absent"))?;
            let (a, b) = (&original.conductor, &next.conductor);
            if a.source_basis() != b.source_basis()
                || a.source_key() != b.source_key()
                || a.length_m().to_bits() != b.length_m().to_bits()
                || a.area_m2().to_bits() != b.area_m2().to_bits()
                || !same_material(a.material(), b.material())
            {
                return Err(invalid(
                    "accepted graph-resistor geometry and sources must match its actual caloric state",
                ));
            }
        }
        let initial = refresh(&accepted.source)?;
        if initial.iter().zip(problem.conductors()).any(|(a, b)| {
            !a.has_same_electrical_source(b) || !same_material(a.material(), b.material())
        }) {
            return Err(invalid(
                "accepted tetrahedral response sources must match the actual thermal state",
            ));
        }
        let contact_sources = problem
            .contacts()
            .iter()
            .map(|c| &c.source)
            .collect::<Vec<_>>();
        let contact_axes = contacts::response_axes(problem.contact_response_sources());
        let contact_owner = contacts::PreparedContacts::new_for_response_sources(
            &route_owner,
            &contact_sources,
            problem.contact_cells(),
            problem.contact_liquid_mass_fraction_axes(),
            contact_bindings,
            &routes,
            &accepted.source,
            &[],
            &contact_axes,
            &cancel,
        )?;
        let (thermal, circuit, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |trial, remaining| {
                check(&cancel)?;
                let coefficients = resolve(trial)?;
                let seeds = refresh(&coefficients.source)?;
                let refreshed_contacts = contact_owner.refresh_sources(
                    &contact_sources,
                    &coefficients.source,
                    &cancel,
                )?;
                let contacts = problem
                    .contacts()
                    .iter()
                    .zip(refreshed_contacts)
                    .map(|(old, source)| SpatialConductionContact {
                        first_vertices: old.first_vertices,
                        second_vertices: old.second_vertices,
                        source,
                    })
                    .collect();
                let candidate = circuit.trial(
                    &problem,
                    seeds,
                    contacts,
                    coefficients.resistors,
                    target_s,
                    ports,
                    matrix,
                    assembly,
                    field,
                    response,
                    controls.circuit,
                    Some(&contact_owner),
                    remaining,
                    &cancel,
                )?;
                let circuit_trials = candidate
                    .step
                    .field_solves
                    .checked_add(candidate.step.iterations)
                    .ok_or_else(|| invalid("native field/circuit work count overflows"))?;
                Ok(DrivenCoefficients {
                    circuit: candidate,
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: coefficients.source,
                    circuit_trials,
                })
            },
            |candidate, _dt, checkpoint| {
                let step = &candidate.step.circuit;
                let mut powers = route_owner.deposit_power(
                    &step.spatial.field.cell_power_w,
                    &routes,
                    checkpoint,
                )?;
                for contact in &step.spatial.contacts {
                    checkpoint()?;
                    contacts::deposit_partition(
                        contact.cells,
                        contact.heat_fraction_to_first,
                        contact.response.power_w,
                        &mut powers,
                    )?;
                }
                for resistor in &step.resistors {
                    checkpoint()?;
                    let cell = *resistor_routes
                        .get(&resistor.material.branch)
                        .ok_or_else(|| invalid("graph resistor heat lost its original receiver"))?;
                    let previous = *powers.get(&cell).unwrap_or(&0.0);
                    let next = finite(previous + resistor.power_w)?;
                    if resistor.power_w != 0.0 && next == previous {
                        return Err(invalid(
                            "material resistor heat disappeared in the receiving-cell sum",
                        ));
                    }
                    powers.insert(cell, next);
                }
                Ok(powers)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(&cancel)?;
        Ok(SpatialConductionResponseCircuitEnthalpyCandidate {
            thermal,
            circuit,
            report,
        })
    }
}
