//! Native multi-electrode conduction inside the original finite C/L descriptor.
//! Every thermal trial starts from the accepted charges/fluxes, refreshes the
//! original field/contact sources, and recovers actual midpoint E/J and heat.

use super::super::circuit::consistent;
use super::*;
use crate::circuit_enthalpy::material::CellResistor;
use crate::material_circuit::{
    MaterialCircuitError, MaterialCircuitSolve, MaterialResistor, ResistorHeat,
    account_resistor_heat, lower_resistors,
};
use fs_em::spatial::conduction::SpatialConductanceMatrix;
use fs_phs::{
    StepRecord,
    circuit::{
        Branch, CircuitDae, CircuitError, CircuitGraph, MutualInductance,
        assemble_circuit_with_couplings_polled, conductance::ConductanceBlock,
    },
};
use std::collections::BTreeSet;

/// Original electrical storage and spatial source history. The graph contains
/// fixed authored C/L, reciprocal mutual inductance, source ports and optional
/// ideal transformers. Each
/// ordinary graph resistor requires its original complete material source and
/// a separate whole-cell thermal receiver; the spatial field retains its own
/// bulk/contact heat. Source-dependent storage requires a reversible-work owner.
pub struct SpatialConductionCircuit<C = ResolvedConductorState> {
    problem: SpatialConduction<C>,
    graph: CircuitGraph,
    mutual: Vec<MutualInductance>,
    terminals: Vec<Vec<usize>>,
    nodes: Vec<usize>,
    fluxes_wb: Vec<f64>,
    charges_c: Vec<f64>,
    time_s: f64,
    material_resistors: Vec<MaterialResistor>,
}

#[derive(Clone, Debug)]
pub struct SpatialConductionCircuitStep {
    pub native: StepRecord,
    pub spatial: SpatialConductionThermalStep,
    /// Actual native midpoint voltages in the original terminal order.
    pub terminal_voltage_v: Vec<f64>,
    /// Native conductance-block currents INTO the spatial body.
    pub terminal_current_a: Vec<f64>,
    /// Original lumped material losses, separate from spatial bulk/contact heat.
    pub resistors: Vec<ResistorHeat>,
    /// Complete field/contact plus lumped material heat, once per physical loss.
    pub heat_j: f64,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Independently evaluated fixed-C/L energy increment, including every
    /// mutual-inductance term. Midpoint constitutive efforts times coordinate
    /// increments avoid subtraction of large stored-energy endpoints.
    pub stored_change_j: f64,
    /// Independent dt*sum(authored source input * native conjugate output).
    pub supplied_work_j: f64,
    pub energy_residual_j: f64,
}

pub struct SpatialConductionCircuitEnthalpyCandidate<
    't,
    C = ResolvedConductorState,
    S = MaterialNetworkSample,
> {
    thermal: DepositionSimulation<'t>,
    circuit: SpatialConductionCircuit<C>,
    report: CircuitEnthalpyReport<S, SpatialConductionCircuitStep>,
}
impl<'t, C: SpatialConductorMaterial, S> SpatialConductionCircuitEnthalpyCandidate<'t, C, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn circuit(&self) -> &SpatialConductionCircuit<C> {
        &self.circuit
    }
    pub fn report(&self) -> &CircuitEnthalpyReport<S, SpatialConductionCircuitStep> {
        &self.report
    }
    /// Publish electrical storage, original source field and finite heat jointly.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        SpatialConductionCircuit<C>,
        CircuitEnthalpyReport<S, SpatialConductionCircuitStep>,
    ) {
        (self.thermal, self.circuit, self.report)
    }
}

struct CircuitCandidate<C> {
    circuit: SpatialConductionCircuit<C>,
    step: SpatialConductionCircuitStep,
}
impl<C: SpatialConductorMaterial> ThermalDriver for CircuitCandidate<C> {
    type Model = SpatialConductionCircuit<C>;
    type Step = SpatialConductionCircuitStep;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &self.step.resistors,
            contacts: &[],
            heat_j: self.step.heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: self.step.stored_change_j,
            supplied_j: self.step.supplied_work_j,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.circuit, self.step)
    }
}

fn electrical(error: CircuitError) -> Error {
    CircuitEnthalpyError::Circuit(MaterialCircuitError::Circuit(error))
}

fn material_circuit(error: MaterialCircuitError) -> Error {
    match error {
        MaterialCircuitError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Circuit(other),
    }
}

fn native_polled<T>(
    cancel: &RefCell<impl FnMut() -> bool>,
    operation: impl FnOnce(&mut dyn FnMut() -> Result<(), CircuitError>) -> Result<T, CircuitError>,
) -> Result<T, Error> {
    check(cancel)?;
    let mut aborted = false;
    let mut checkpoint = || {
        if (*cancel.borrow_mut())() {
            aborted = true;
            Err(CircuitError::Invalid {
                what: "spatial circuit preparation cancelled",
            })
        } else {
            Ok(())
        }
    };
    let result = operation(&mut checkpoint);
    if aborted {
        return Err(CircuitEnthalpyError::Cancelled);
    }
    check(cancel)?;
    result.map_err(electrical)
}

#[allow(clippy::too_many_arguments)]
fn prepare<'p, C: SpatialConductorMaterial>(
    problem: &'p SpatialConduction<C>,
    graph: &CircuitGraph,
    mutual: &[MutualInductance],
    terminals: &[Vec<usize>],
    nodes: &[usize],
    response: SpatialConductanceControls,
    field_controls: ConductionControls,
    cancel: &RefCell<impl FnMut() -> bool>,
) -> Result<
    (
        SpatialConductanceMatrix<'p, C>,
        ConductanceBlock,
        CircuitDae,
    ),
    Error,
> {
    check(cancel)?;
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
    let faces: Vec<_> = terminals.iter().map(Vec::as_slice).collect();
    let matrix = problem
        .conductance_matrix(&faces, response, field_controls, &mut || {
            (*cancel.borrow_mut())()
        })
        .map_err(field_error)?;
    let n = nodes.len();
    let mut grounded = Vec::with_capacity((n - 1) * (n - 1));
    for row in 1..n {
        for column in 1..n {
            grounded.push(matrix.conductance_s()[row * n + column]);
        }
    }
    let block = native_polled(cancel, |checkpoint| {
        ConductanceBlock::try_new_polled(nodes.to_vec(), grounded, checkpoint)
    })?;
    let dae = native_polled(cancel, |checkpoint| {
        assemble_circuit_with_couplings_polled(
            graph,
            mutual,
            std::slice::from_ref(&block),
            checkpoint,
        )
    })?;
    Ok((matrix, block, dae))
}

impl<C: SpatialConductorMaterial> SpatialConductionCircuit<C> {
    /// Bind a full native terminal matrix to original circuit nodes. Nodes and
    /// face sets have the same order; terminal zero is only the FIELD reference
    /// and need not map to circuit ground. No pair-resistor reduction is used.
    /// Charges/fluxes retain graph branch order and original physical units.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        problem: SpatialConduction<C>,
        graph: CircuitGraph,
        terminals: Vec<Vec<usize>>,
        nodes: Vec<usize>,
        time_s: f64,
        fluxes_wb: Vec<f64>,
        charges_c: Vec<f64>,
        ports: &[f64],
        response: SpatialConductanceControls,
        field_controls: ConductionControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        Self::new_with_material_resistors(
            problem,
            graph,
            terminals,
            nodes,
            Vec::new(),
            time_s,
            fluxes_wb,
            charges_c,
            ports,
            response,
            field_controls,
            cancelled,
        )
    }

    /// Bind a complete original source to every ordinary resistor branch.
    /// Authored numeric ohms are replaced by these material/geometry values.
    /// Sources may be obtained from `MaterialNetwork::resolve_with_bases` at
    /// the accepted thermal state. Advance with `step_material_circuit_candidate`
    /// so every source follows its own receiving cell's temperature and phase.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_material_resistors(
        problem: SpatialConduction<C>,
        graph: CircuitGraph,
        terminals: Vec<Vec<usize>>,
        nodes: Vec<usize>,
        material_resistors: Vec<MaterialResistor>,
        time_s: f64,
        fluxes_wb: Vec<f64>,
        charges_c: Vec<f64>,
        ports: &[f64],
        response: SpatialConductanceControls,
        field_controls: ConductionControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        Self::new_with_couplings(
            problem,
            graph,
            terminals,
            nodes,
            material_resistors,
            Vec::new(),
            time_s,
            fluxes_wb,
            charges_c,
            ports,
            response,
            field_controls,
            cancelled,
        )
    }

    /// Couple fixed reciprocal windings to the complete spatial conductor and
    /// original material resistors. Magnetic links use graph branch indices,
    /// including their winding polarity, and persist across every thermal
    /// trial and accepted endpoint. Initial fluxes are TOTAL linkages L i;
    /// current is L^-1 flux, including all mutual terms.
    ///
    /// The full inductance matrix must have positive leakage energy. Its
    /// coefficients are fixed here: changing permeability or geometry needs
    /// a separate reversible parameter-work owner and is not inferred from
    /// the spatial conductor's temperature.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_couplings(
        problem: SpatialConduction<C>,
        graph: CircuitGraph,
        terminals: Vec<Vec<usize>>,
        nodes: Vec<usize>,
        material_resistors: Vec<MaterialResistor>,
        mutual: Vec<MutualInductance>,
        time_s: f64,
        fluxes_wb: Vec<f64>,
        charges_c: Vec<f64>,
        ports: &[f64],
        response: SpatialConductanceControls,
        field_controls: ConductionControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        if !time_s.is_finite() || time_s < 0.0 {
            return Err(invalid("finite nonnegative initial circuit time required"));
        }
        let (graph, material_resistors) =
            lower_resistors(&graph, material_resistors).map_err(material_circuit)?;
        let (matrix, _, dae) = prepare(
            &problem,
            &graph,
            &mutual,
            &terminals,
            &nodes,
            response,
            field_controls,
            &cancel,
        )?;
        consistent(&dae, &fluxes_wb, &charges_c, ports)?;
        let terminals = matrix.terminal_faces().to_vec();
        drop(matrix);
        check(&cancel)?;
        Ok(Self {
            problem,
            graph,
            mutual,
            terminals,
            nodes,
            fluxes_wb,
            charges_c,
            time_s,
            material_resistors,
        })
    }
    pub fn problem(&self) -> &SpatialConduction<C> {
        &self.problem
    }
    pub fn graph(&self) -> &CircuitGraph {
        &self.graph
    }
    /// Fixed reciprocal links in original graph-branch order.
    pub fn mutual_inductance(&self) -> &[MutualInductance] {
        &self.mutual
    }
    pub fn terminal_faces(&self) -> &[Vec<usize>] {
        &self.terminals
    }
    pub fn terminal_nodes(&self) -> &[usize] {
        &self.nodes
    }
    pub fn fluxes_wb(&self) -> &[f64] {
        &self.fluxes_wb
    }
    pub fn charges_c(&self) -> &[f64] {
        &self.charges_c
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    /// Complete branch-ordered conductor sources at the accepted thermal state.
    pub fn material_resistors(&self) -> &[MaterialResistor] {
        &self.material_resistors
    }

    #[allow(clippy::too_many_arguments)]
    fn trial(
        &self,
        problem: SpatialConduction<C>,
        material_resistors: Vec<MaterialResistor>,
        target_s: f64,
        ports: &[f64],
        response: SpatialConductanceControls,
        field_controls: ConductionControls,
        controls: MaterialCircuitSolve,
        contacts: &contacts::PreparedContacts<'_, '_>,
        remaining: usize,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<(CircuitCandidate<C>, usize), Error> {
        check(cancel)?;
        // N-1 FE basis solves plus one native finite circuit step. Admission is
        // before any basis extraction; the original driver accounts this work.
        let trials = self.terminals.len();
        if trials > remaining {
            return Err(CircuitEnthalpyError::CircuitBudget);
        }
        let dt = finite(target_s - self.time_s)?;
        if dt <= 0.0 {
            return Err(invalid("positive spatial circuit interval required"));
        }
        let (graph, material_resistors) =
            lower_resistors(&self.graph, material_resistors).map_err(material_circuit)?;
        let (matrix, block, dae) = prepare(
            &problem,
            &graph,
            &self.mutual,
            &self.terminals,
            &self.nodes,
            response,
            field_controls,
            cancel,
        )?;
        let x0 = consistent(&dae, &self.fluxes_wb, &self.charges_c, ports)?;
        check(cancel)?;
        let (native, _) = dae.step_audited(&x0, ports, dt).map_err(electrical)?;
        check(cancel)?;
        let stored_before_j = dae.system.hamiltonian(&x0);
        let stored_after_j = dae.system.hamiltonian(&native.x);
        if native.x.len() != x0.len()
            || native.y.len() != ports.len()
            || native.x.iter().chain(&native.y).any(|v| !v.is_finite())
            || [
                stored_before_j,
                stored_after_j,
                native.delta_h,
                native.dissipated,
                native.supplied,
                native.solver_residual,
            ]
            .iter()
            .any(|v| !v.is_finite())
            || stored_before_j < 0.0
            || stored_after_j < 0.0
            || native.dissipated < 0.0
        {
            return Err(invalid(
                "native spatial circuit storage/work is not representable",
            ));
        }
        let terminal_voltage_v = self
            .nodes
            .iter()
            .map(|&node| {
                if node == 0 {
                    0.0
                } else {
                    let v = dae.node_potential_index[node - 1];
                    f64::midpoint(x0[v], native.x[v])
                }
            })
            .collect::<Vec<_>>();
        let terminal_current_a = block
            .terminal_currents_a(&terminal_voltage_v)
            .map_err(electrical)?;
        let factored_power = block.power_w(&terminal_voltage_v).map_err(electrical)?;
        // Independent native node work. Individual electrode contributions
        // may have either sign; retain their cancellation in the field gauge.
        let (mut terminal_power, mut correction) = (0.0, 0.0);
        for (&voltage, &current) in terminal_voltage_v.iter().zip(&terminal_current_a) {
            check(cancel)?;
            let work = finite(current * finite(voltage - terminal_voltage_v[0])?)?;
            let adjusted = finite(work - correction)?;
            let next = finite(terminal_power + adjusted)?;
            correction = finite((next - terminal_power) - adjusted)?;
            terminal_power = next;
        }
        let field = matrix
            .field_at_voltages(
                &terminal_voltage_v,
                field_controls.current_tolerance_a,
                field_controls.power_tolerance_w,
                &mut || (*cancel.borrow_mut())(),
            )
            .map_err(field_error)?;
        for (faces, &expected) in matrix.terminal_faces().iter().zip(&terminal_current_a) {
            let vertices: BTreeSet<_> = faces
                .iter()
                .flat_map(|&face| problem.mesh().faces[face])
                .collect();
            let mut reaction = 0.0;
            for vertex in vertices {
                check(cancel)?;
                reaction = finite(
                    reaction
                        + field.nodal_boundary_current_a[vertex as usize].ok_or_else(|| {
                            invalid("native circuit electrode lost its field reaction")
                        })?,
                )?;
            }
            if finite(reaction - expected)?.abs() > field_controls.current_tolerance_a {
                return Err(invalid(
                    "native circuit and original field terminal currents disagree",
                ));
            }
        }
        drop(matrix);
        let spatial_heat_j = finite(dt * field.power_w)?;
        let terminal_work = finite(dt * terminal_power)?;
        if (field.power_w > 0.0 && spatial_heat_j == 0.0)
            || finite(terminal_power - factored_power)?.abs() > controls.power_tolerance_w
            || finite(field.power_w - terminal_power)?.abs() > controls.power_tolerance_w
            || finite(spatial_heat_j - terminal_work)?.abs() > controls.energy_tolerance_j
        {
            return Err(invalid(
                "native field/contact heat and terminal work disagree",
            ));
        }
        let (resistors, _, resistor_heat_j) = account_resistor_heat(
            &graph,
            &material_resistors,
            &[],
            &dae.node_potential_index,
            &x0,
            &native,
            dt,
            controls,
            &mut || (*cancel.borrow_mut())(),
        )
        .map_err(material_circuit)?;
        let heat_j = finite(spatial_heat_j + resistor_heat_j)?;
        if finite(heat_j - native.dissipated)?.abs() > controls.energy_tolerance_j {
            return Err(invalid(
                "spatial/contact and material resistor heat do not close native dissipation",
            ));
        }
        let fluxes_wb: Vec<_> = dae.flux_index.iter().map(|&i| native.x[i]).collect();
        let charges_c: Vec<_> = dae.charge_index.iter().map(|&i| native.x[i]).collect();
        let mean_currents = if self.mutual.is_empty() {
            None
        } else {
            // Independent constitutive work: dH = (L^-1 mean(lambda)) dot
            // d(lambda). Dividing by individual self L discards magnetic
            // exchange and can misclassify it as heat or source power.
            let midpoint = x0
                .iter()
                .zip(&native.x)
                .map(|(&a, &b)| f64::midpoint(a, b))
                .collect::<Vec<_>>();
            Some(dae.inductor_currents(&midpoint).map_err(electrical)?)
        };
        let (mut l, mut c, mut stored_change_j) = (0, 0, 0.0);
        for &(_, _, branch) in &graph.branches {
            check(cancel)?;
            let work = match branch {
                Branch::Inductor { henries } => {
                    let mean_current = mean_currents.as_ref().map_or_else(
                        || f64::midpoint(self.fluxes_wb[l], fluxes_wb[l]) / henries,
                        |currents| currents[l],
                    );
                    let work = finite(mean_current * finite(fluxes_wb[l] - self.fluxes_wb[l])?)?;
                    l += 1;
                    work
                }
                Branch::Capacitor { farads } => {
                    let work = finite(
                        (f64::midpoint(self.charges_c[c], charges_c[c]) / farads)
                            * finite(charges_c[c] - self.charges_c[c])?,
                    )?;
                    c += 1;
                    work
                }
                _ => 0.0,
            };
            stored_change_j = finite(stored_change_j + work)?;
        }
        let mut source_power = 0.0;
        for (&input, &output) in ports.iter().zip(&native.y) {
            source_power = finite(input.mul_add(output, source_power))?;
        }
        let supplied_work_j = finite(dt * source_power)?;
        let energy_residual_j = finite(stored_change_j + heat_j - supplied_work_j)?;
        if finite(stored_change_j - native.delta_h)?.abs() > controls.energy_tolerance_j
            || finite(supplied_work_j - native.supplied)?.abs() > controls.energy_tolerance_j
            || energy_residual_j.abs() > controls.energy_tolerance_j
        {
            return Err(invalid(
                "independent fixed-C/L storage, local heat and authored port work do not close",
            ));
        }
        let contact_heat = contacts.heat(&field.contacts)?;
        let spatial = SpatialConductionThermalStep {
            start_s: self.time_s,
            end_s: target_s,
            field,
            contacts: contact_heat,
            heat_j: spatial_heat_j,
            supplied_work_j: terminal_work,
        };
        check(cancel)?;
        Ok((
            CircuitCandidate {
                circuit: Self {
                    problem,
                    graph,
                    mutual: self.mutual.clone(),
                    terminals: self.terminals.clone(),
                    nodes: self.nodes.clone(),
                    fluxes_wb,
                    charges_c,
                    time_s: target_s,
                    material_resistors,
                },
                step: SpatialConductionCircuitStep {
                    native,
                    spatial,
                    terminal_voltage_v,
                    terminal_current_a,
                    resistors,
                    heat_j,
                    stored_before_j,
                    stored_after_j,
                    stored_change_j,
                    supplied_work_j,
                    energy_residual_j,
                },
            },
            trials,
        ))
    }
}

impl SpatialConductionThermalNetwork<'_> {
    /// Advance the actual multi-electrode field, native capacitor/inductor
    /// history and finite enthalpy on one accepted endpoint. Every trial uses
    /// the SAME accepted C/L coordinates and source/contact/heat ownership.
    /// Native midpoint electrical work enters the existing backward-Euler heat
    /// residual once; source and temporal accuracy still require refinement.
    #[allow(clippy::too_many_arguments)]
    pub fn step_circuit_candidate<'t, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &SpatialConductionCircuit<C>,
        ports: &[f64],
        contact_bindings: &[ConductionContactThermal<'_>],
        response: SpatialConductanceControls,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialConductionCircuitEnthalpyCandidate<'t, C>, Error> {
        self.step_material_circuit_candidate(
            thermal,
            target_s,
            circuit,
            ports,
            &[],
            contact_bindings,
            response,
            assembly,
            field_controls,
            controls,
            initial_energy_increments_j,
            cancelled,
        )
    }

    /// Evolve ordinary material resistors and the multi-electrode field in one
    /// original descriptor/enthalpy step. Each resistor's full source bundle,
    /// electrical basis, geometry and receiving caloric cell remain explicit.
    /// Whole-cell resistor volumes cannot also belong to the spatial field.
    /// Bulk/contact/branch heat stays separate and is deposited exactly once.
    #[allow(clippy::too_many_arguments)]
    pub fn step_material_circuit_candidate<'t, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &SpatialConductionCircuit<C>,
        ports: &[f64],
        material_resistors: &[CellResistor],
        contact_bindings: &[ConductionContactThermal<'_>],
        response: SpatialConductanceControls,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialConductionCircuitEnthalpyCandidate<'t, C>, Error> {
        let cancel = RefCell::new(cancelled);
        self.admit_circuit_step(
            thermal,
            target_s,
            circuit,
            material_resistors.len(),
            assembly,
            controls,
            initial_energy_increments_j,
            &cancel,
        )?;
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: material_resistors,
            links: self.links,
        };
        let bases = circuit
            .material_resistors
            .iter()
            .map(|resistor| (resistor.branch, resistor.conductor.source_basis()))
            .collect::<Vec<_>>();
        let prepared = network
            .prepare_with_bases(thermal, &bases, || (*cancel.borrow_mut())())
            .map_err(material)?;
        let destinations = material_resistors
            .iter()
            .map(|resistor| resistor.destination)
            .collect::<Vec<_>>();
        self.step_circuit_resolved_candidate(
            thermal,
            target_s,
            circuit,
            ports,
            &destinations,
            contact_bindings,
            response,
            assembly,
            field_controls,
            controls,
            initial_energy_increments_j,
            |trial| {
                prepared
                    .resolve(trial, || (*cancel.borrow_mut())())
                    .map_err(material)
            },
            |sample| sample,
            &cancel,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn admit_circuit_step<C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'_>,
        target_s: f64,
        circuit: &SpatialConductionCircuit<C>,
        resistor_count: usize,
        assembly: WeightedAssemblyLimits,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<(), Error> {
        check(cancel)?;
        if thermal.time_s() != circuit.time_s
            || resistor_count != circuit.material_resistors.len()
            || !controls.circuit.energy_tolerance_j.is_finite()
            || controls.circuit.energy_tolerance_j <= 0.0
            || !controls.circuit.power_tolerance_w.is_finite()
            || controls.circuit.power_tolerance_w <= 0.0
        {
            return Err(invalid(
                "common accepted time, complete material resistor bindings and finite positive work tolerances required",
            ));
        }
        self.admit_step(
            thermal,
            target_s,
            &circuit.problem,
            assembly,
            controls,
            initial_energy_increments_j,
        )?;
        Ok(())
    }

    // Bare Fourier and complete heat assemblies share the same native circuit
    // proposal, accepted source checks and exactly-once local heat routing.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::circuit_enthalpy::conduction) fn step_circuit_resolved_candidate<
        't,
        C: SpatialConductorMaterial,
        S,
        E: From<ConductionThermalError>,
    >(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &SpatialConductionCircuit<C>,
        ports: &[f64],
        destinations: &[crate::circuit_enthalpy::ResistorCell],
        contact_bindings: &[ConductionContactThermal<'_>],
        response: SpatialConductanceControls,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        resolve: impl Fn(
            &DepositionSimulation<'t>,
        ) -> Result<CircuitThermalCoefficients<S>, CircuitEnthalpyError<E>>,
        sample: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<SpatialConductionCircuitEnthalpyCandidate<'t, C, S>, CircuitEnthalpyError<E>> {
        self.admit_circuit_step(
            thermal,
            target_s,
            circuit,
            destinations.len(),
            assembly,
            controls,
            initial_energy_increments_j,
            cancel,
        )
        .map_err(map_conduction_error)?;
        if !destinations.is_empty() {
            let accepted = resolve(thermal)?;
            for original in &circuit.material_resistors {
                check(cancel).map_err(map_conduction_error)?;
                let next = accepted
                    .resistors
                    .iter()
                    .find(|resistor| resistor.branch == original.branch)
                    .ok_or_else(|| invalid("accepted material-resistor branch is absent"))
                    .map_err(map_conduction_error)?;
                let (a, b) = (&original.conductor, &next.conductor);
                if a.source_basis() != b.source_basis()
                    || a.source_key() != b.source_key()
                    || a.length_m().to_bits() != b.length_m().to_bits()
                    || a.area_m2().to_bits() != b.area_m2().to_bits()
                    || !same_material(a.material(), b.material())
                {
                    return Err(map_conduction_error(invalid(
                        "accepted graph-resistor geometry and sources must match its actual caloric state",
                    )));
                }
            }
        }
        let (thermal, circuit, report) = self.step_driver_candidate_with_geometry(
            thermal,
            target_s,
            &circuit.problem,
            contact_bindings,
            assembly,
            controls,
            initial_energy_increments_j,
            resolve,
            sample,
            None,
            |_| Ok(None),
            &[],
            |_| Ok(Vec::new()),
            destinations,
            |problem, resistors, contacts, remaining| {
                circuit
                    .trial(
                        problem,
                        resistors,
                        target_s,
                        ports,
                        response,
                        field_controls,
                        controls.circuit,
                        contacts,
                        remaining,
                        cancel,
                    )
                    .map_err(map_conduction_error)
            },
            |candidate| &candidate.step.spatial,
            cancel,
        )?;
        Ok(SpatialConductionCircuitEnthalpyCandidate {
            thermal,
            circuit,
            report,
        })
    }
}
