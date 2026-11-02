//! Native multi-electrode conduction inside the original finite C/L descriptor.
//! Every thermal trial starts from the accepted charges/fluxes, refreshes the
//! original field/contact sources, and recovers actual midpoint E/J and heat.

use super::super::circuit::consistent;
use super::*;
use crate::material_circuit::{MaterialCircuitError, MaterialCircuitSolve, admit_graph};
use fs_em::spatial::conduction::SpatialConductanceMatrix;
use fs_phs::{
    StepRecord,
    circuit::{
        Branch, CircuitDae, CircuitError, CircuitGraph, assemble_circuit_with_conductance_polled,
        conductance::ConductanceBlock,
    },
};
use std::collections::BTreeSet;

/// Original electrical storage and spatial source history. The graph contains
/// fixed authored C/L, source ports and optional ideal transformers. Dissipation
/// belongs solely to this field; separate graph resistors and source-dependent
/// storage require their own heat/work bindings and are explicitly excluded.
pub struct SpatialConductionCircuit<C = ResolvedConductorState> {
    problem: SpatialConduction<C>,
    graph: CircuitGraph,
    terminals: Vec<Vec<usize>>,
    nodes: Vec<usize>,
    fluxes_wb: Vec<f64>,
    charges_c: Vec<f64>,
    time_s: f64,
}

#[derive(Clone, Debug)]
pub struct SpatialConductionCircuitStep {
    pub native: StepRecord,
    pub spatial: SpatialConductionThermalStep,
    /// Actual native midpoint voltages in the original terminal order.
    pub terminal_voltage_v: Vec<f64>,
    /// Native conductance-block currents INTO the spatial body.
    pub terminal_current_a: Vec<f64>,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Independently evaluated fixed-C/L energy increment, with cancellation
    /// avoided by midpoint-coordinate times coordinate-increment products.
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
            resistors: &[],
            contacts: &[],
            heat_j: self.step.spatial.heat_j,
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
        assemble_circuit_with_conductance_polled(graph, std::slice::from_ref(&block), checkpoint)
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
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        if !time_s.is_finite()
            || time_s < 0.0
            || admit_graph(&graph).map_err(CircuitEnthalpyError::Circuit)? != 0
        {
            return Err(invalid(
                "finite initial time and a bounded graph without separate resistor losses required",
            ));
        }
        let (matrix, _, dae) = prepare(
            &problem,
            &graph,
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
            terminals,
            nodes,
            fluxes_wb,
            charges_c,
            time_s,
        })
    }
    pub fn problem(&self) -> &SpatialConduction<C> {
        &self.problem
    }
    pub fn graph(&self) -> &CircuitGraph {
        &self.graph
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

    #[allow(clippy::too_many_arguments)]
    fn trial(
        &self,
        problem: SpatialConduction<C>,
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
        let (matrix, block, dae) = prepare(
            &problem,
            &self.graph,
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
        let heat_j = finite(dt * field.power_w)?;
        let terminal_work = finite(dt * terminal_power)?;
        if (field.power_w > 0.0 && heat_j == 0.0)
            || finite(terminal_power - factored_power)?.abs() > controls.power_tolerance_w
            || finite(field.power_w - terminal_power)?.abs() > controls.power_tolerance_w
            || finite(heat_j - terminal_work)?.abs() > controls.energy_tolerance_j
            || finite(heat_j - native.dissipated)?.abs() > controls.energy_tolerance_j
        {
            return Err(invalid(
                "native field/contact heat and circuit dissipation disagree",
            ));
        }
        let fluxes_wb: Vec<_> = dae.flux_index.iter().map(|&i| native.x[i]).collect();
        let charges_c: Vec<_> = dae.charge_index.iter().map(|&i| native.x[i]).collect();
        let (mut l, mut c, mut stored_change_j) = (0, 0, 0.0);
        for &(_, _, branch) in &self.graph.branches {
            check(cancel)?;
            let work = match branch {
                Branch::Inductor { henries } => {
                    let work = finite(
                        (f64::midpoint(self.fluxes_wb[l], fluxes_wb[l]) / henries)
                            * finite(fluxes_wb[l] - self.fluxes_wb[l])?,
                    )?;
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
            heat_j,
            supplied_work_j: terminal_work,
        };
        check(cancel)?;
        Ok((
            CircuitCandidate {
                circuit: Self {
                    problem,
                    graph: self.graph.clone(),
                    terminals: self.terminals.clone(),
                    nodes: self.nodes.clone(),
                    fluxes_wb,
                    charges_c,
                    time_s: target_s,
                },
                step: SpatialConductionCircuitStep {
                    native,
                    spatial,
                    terminal_voltage_v,
                    terminal_current_a,
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
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        if thermal.time_s() != circuit.time_s
            || !controls.circuit.energy_tolerance_j.is_finite()
            || controls.circuit.energy_tolerance_j <= 0.0
            || !controls.circuit.power_tolerance_w.is_finite()
            || controls.circuit.power_tolerance_w <= 0.0
        {
            return Err(invalid(
                "common accepted time and finite positive circuit work tolerances required",
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
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let (thermal, circuit, report) = self.step_driver_candidate_with_geometry(
            thermal,
            target_s,
            &circuit.problem,
            contact_bindings,
            assembly,
            controls,
            initial_energy_increments_j,
            |trial| {
                prepared
                    .resolve(trial, || (*cancel.borrow_mut())())
                    .map_err(material)
            },
            |sample| sample,
            None,
            |_| Ok(None),
            &[],
            |_| Ok(Vec::new()),
            |problem, contacts, remaining| {
                circuit.trial(
                    problem,
                    target_s,
                    ports,
                    response,
                    field_controls,
                    controls.circuit,
                    contacts,
                    remaining,
                    &cancel,
                )
            },
            |candidate| &candidate.step.spatial,
            &cancel,
        )?;
        Ok(SpatialConductionCircuitEnthalpyCandidate {
            thermal,
            circuit,
            report,
        })
    }
}
