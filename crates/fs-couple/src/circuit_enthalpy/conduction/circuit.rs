//! A spatial Ohmic body inside the original transient Kirchhoff/thermal solve.
//!
//! The field's native two-electrode resistance stamps one graph branch at every
//! proposed thermal endpoint. The original descriptor advances capacitor charge
//! and inductor flux; its actual midpoint branch voltage drives the original
//! spatial E/J solution. Local E.J heat and every other graph-resistor loss enter
//! the existing enthalpy residual once. The common candidate publishes all three
//! states together. No uniform-wire substitute or additional integrator is used.
//! Scalar and oriented tensor materials use that same field extraction; both
//! principal source properties and the fixed frame survive every thermal trial.

use super::*;
use crate::circuit_enthalpy::ResistorCell;
use crate::material_circuit::{MaterialCircuitError, MaterialCircuitSolve, admit_graph, extended};
use fs_phs::{
    StepRecord,
    circuit::{Branch, CircuitDae, CircuitError, CircuitGraph, assemble_circuit},
};

/// One source-resolved two-electrode body replacing an original resistor slot.
/// Positive/negative electrode orientation follows that branch's first/second
/// graph node. Each electrode must be a union of complete exterior mesh edges.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ConductionCircuitTerminal {
    pub branch: usize,
    pub positive: Vec<usize>,
    pub negative: Vec<usize>,
}

/// One accepted electrical history and its complete spatial material operator.
/// Only native fluxes and charges carry electrical storage. The field remains
/// quasistatic at the circuit's common midpoint; C/L and other graph resistors
/// are explicitly authored fixed parameters on this bounded rung.
pub struct ConductionCircuit<C = ResolvedConductorState> {
    problem: PlanarConduction<C>,
    graph: CircuitGraph,
    terminal: ConductionCircuitTerminal,
    fluxes_wb: Vec<f64>,
    charges_c: Vec<f64>,
    time_s: f64,
}

fn electrical(error: CircuitError) -> Error {
    CircuitEnthalpyError::Circuit(MaterialCircuitError::Circuit(error))
}

fn consistent(
    dae: &CircuitDae,
    fluxes: &[f64],
    charges: &[f64],
    ports: &[f64],
) -> Result<Vec<f64>, Error> {
    let seed = extended(dae, fluxes, charges, ports)?;
    let initial = dae
        .consistent_initial_state(&seed, ports)
        .map_err(electrical)?;
    if initial.len() != seed.len()
        || initial.iter().any(|x| !x.is_finite())
        || !dae.system.hamiltonian(&initial).is_finite()
        || dae.system.hamiltonian(&initial) < 0.0
    {
        return Err(invalid(
            "spatial circuit initial storage or constraints are unrepresentable",
        ));
    }
    for i in 0..dae.system.differential_dim() {
        if (initial[i] - seed[i]).abs() > 64.0 * f64::EPSILON * seed[i].abs() {
            return Err(invalid(
                "circuit constraint initialization changed a conserved coordinate",
            ));
        }
    }
    Ok(initial)
}

fn stamp(
    graph: &CircuitGraph,
    terminal: &ConductionCircuitTerminal,
    ohms: f64,
) -> Result<CircuitGraph, Error> {
    let mut next = graph.clone();
    let branch = next
        .branches
        .get_mut(terminal.branch)
        .ok_or_else(|| invalid("spatial terminal branch is absent"))?;
    if !matches!(branch.2, Branch::Resistor { .. })
        || !ohms.is_finite()
        || ohms <= 0.0
        || !ohms.recip().is_finite()
        || ohms.recip() <= 0.0
    {
        return Err(invalid(
            "spatial terminal needs a resistor slot and representable conductance",
        ));
    }
    branch.2 = Branch::Resistor { ohms };
    Ok(next)
}

impl<C: PlanarConductorMaterial> ConductionCircuit<C> {
    /// Bind the native spatial resistance and initialize conserved branch-order
    /// fluxes/charges using the original descriptor. The named resistor slot's
    /// numeric value is replaced by the field extraction; it is not a fallback.
    /// Other graph resistors remain fixed authored components and need explicit
    /// receiving thermal cells when advancing the coupled state.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        problem: PlanarConduction<C>,
        graph: CircuitGraph,
        terminal: ConductionCircuitTerminal,
        time_s: f64,
        fluxes_wb: Vec<f64>,
        charges_c: Vec<f64>,
        ports: &[f64],
        field_controls: ConductionControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        admit_graph(&graph)?;
        if !time_s.is_finite()
            || time_s < 0.0
            || !matches!(
                graph.branches.get(terminal.branch),
                Some((_, _, Branch::Resistor { .. }))
            )
        {
            return Err(invalid(
                "finite initial time and an original resistor terminal are required",
            ));
        }
        let resistance = problem
            .resistance(
                &terminal.positive,
                &terminal.negative,
                field_controls,
                &mut cancelled,
            )
            .map_err(field_error)?;
        let graph = stamp(&graph, &terminal, resistance.resistance_ohm())?;
        let dae = assemble_circuit(&graph).map_err(electrical)?;
        consistent(&dae, &fluxes_wb, &charges_c, ports)?;
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        Ok(Self {
            problem,
            graph,
            terminal,
            fluxes_wb,
            charges_c,
            time_s,
        })
    }

    pub fn problem(&self) -> &PlanarConduction<C> {
        &self.problem
    }
    pub fn graph(&self) -> &CircuitGraph {
        &self.graph
    }
    pub fn terminal(&self) -> &ConductionCircuitTerminal {
        &self.terminal
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
}

/// An explicitly authored fixed resistor outside the represented spatial body.
#[derive(Clone, Debug, PartialEq)]
pub struct ConductionCircuitResistorHeat {
    pub branch: usize,
    pub voltage_v: f64,
    pub current_a: f64,
    pub power_w: f64,
    pub heat_j: f64,
}

/// Original native electrical record, spatial fields and complete heat split.
/// Spatial terminal work is INTERNAL to this circuit/field accounting boundary;
/// only native.supplied represents external electrical source work.
#[derive(Clone, Debug)]
pub struct ConductionCircuitStep {
    pub native: StepRecord,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub spatial: ConductionThermalStep,
    pub spatial_voltage_v: f64,
    pub spatial_current_a: f64,
    pub graph_resistors: Vec<ConductionCircuitResistorHeat>,
    pub heat_j: f64,
    pub energy_residual_j: f64,
}

struct CircuitFieldCandidate<C> {
    circuit: ConductionCircuit<C>,
    step: ConductionCircuitStep,
}
impl<C: PlanarConductorMaterial> ThermalDriver for CircuitFieldCandidate<C> {
    type Model = ConductionCircuit<C>;
    type Step = ConductionCircuitStep;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.step.heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: self.step.stored_after_j - self.step.stored_before_j,
            supplied_j: self.step.native.supplied,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.circuit, self.step)
    }
}

/// Joint proposal for finite heat, conserved circuit storage and spatial sources.
pub struct ConductionCircuitEnthalpyCandidate<'t, C = ResolvedConductorState> {
    thermal: DepositionSimulation<'t>,
    circuit: ConductionCircuit<C>,
    report: CircuitEnthalpyReport<MaterialNetworkSample, ConductionCircuitStep>,
}
impl<'t, C: PlanarConductorMaterial> ConductionCircuitEnthalpyCandidate<'t, C> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn circuit(&self) -> &ConductionCircuit<C> {
        &self.circuit
    }
    pub fn report(&self) -> &CircuitEnthalpyReport<MaterialNetworkSample, ConductionCircuitStep> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        ConductionCircuit<C>,
        CircuitEnthalpyReport<MaterialNetworkSample, ConductionCircuitStep>,
    ) {
        (self.thermal, self.circuit, self.report)
    }
}

impl ConductionThermalNetwork<'_> {
    /// Advance circuit storage, spatial conduction and finite heat together.
    ///
    /// Every trial starts from the original accepted charges/fluxes. Its native
    /// electrical midpoint uses source-pinned ENDPOINT conductivity, while the
    /// existing thermal owner uses backward Euler. This is generally first-order
    /// coupled accuracy and needs time/mesh refinement. Source ports are held
    /// over the interval; topology, geometry and all C/L values remain fixed.
    ///
    /// All triangles have their original full-volume thermal receivers; every
    /// other graph resistor has exactly one additional receiver. Sources and
    /// caloric cards must match the accepted state before any trial. Field and
    /// descriptor current/work, local heat, and combined storage are admitted
    /// independently. No accepted state changes on any refusal or cancellation.
    /// Confined equal-density phase charts use the same explicit liquid-fraction
    /// coordinate and latent-enthalpy admission as the directly driven field;
    /// conductivity changes on an isothermal plateau still restamp the circuit.
    /// This does not add displacement current, induction inside the conductor,
    /// contact resistance, thermoelectric transport or a moving/phase-changing
    /// geometry model to the quasistatic spatial field.
    #[allow(clippy::too_many_arguments)]
    pub fn step_circuit_candidate<'t, C: PlanarConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &ConductionCircuit<C>,
        ports: &[f64],
        graph_resistors: &[ResistorCell],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ConductionCircuitEnthalpyCandidate<'t, C>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        super::super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        let dt = finite(target_s - thermal.time_s())?;
        if dt <= 0.0
            || circuit.time_s.to_bits() != thermal.time_s().to_bits()
            || self.cells.len() != thermal.cells().len()
            || !controls.circuit.energy_tolerance_j.is_finite()
            || controls.circuit.energy_tolerance_j <= 0.0
            || !controls.circuit.power_tolerance_w.is_finite()
            || controls.circuit.power_tolerance_w <= 0.0
        {
            return Err(invalid(
                "matching clocks, positive interval/tolerances and complete thermal bindings required",
            ));
        }
        let problem = &circuit.problem;
        if problem.mesh().faces().len() > assembly.max_cells
            || problem.mesh().vertices().len() > assembly.max_dofs
            || problem.mesh().edges().len() > assembly.max_dofs
            || problem
                .mesh()
                .faces()
                .len()
                .checked_mul(9)
                .is_none_or(|n| n > assembly.max_triplets)
        {
            return Err(invalid(
                "spatial circuit geometry exceeds declared assembly capacity",
            ));
        }
        let mut graph_routes = BTreeMap::new();
        if graph_resistors.len() > circuit.graph.branches.len() {
            return Err(invalid("too many graph-resistor heat destinations"));
        }
        for route in graph_resistors {
            check(&cancel)?;
            if route.branch == circuit.terminal.branch
                || !matches!(
                    circuit.graph.branches.get(route.branch),
                    Some((_, _, Branch::Resistor { .. }))
                )
                || !thermal.cells().contains_key(&route.cell)
                || graph_routes.insert(route.branch, route.cell).is_some()
            {
                return Err(invalid(
                    "each nonspatial graph resistor needs one unique active heat receiver",
                ));
            }
        }
        for (branch, (_, _, law)) in circuit.graph.branches.iter().enumerate() {
            if branch != circuit.terminal.branch
                && matches!(law, Branch::Resistor { .. })
                && !graph_routes.contains_key(&branch)
            {
                return Err(invalid("a graph-resistor heat destination is missing"));
            }
        }
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let routes = self.admit_routes(thermal, problem, &cancel)?;
        let accepted = prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        let initial = self.materials_at(problem, &accepted.source, &routes, &cancel)?;
        if initial.iter().zip(problem.conductors()).any(|(a, b)| {
            !a.has_same_electrical_source(b)
                || a.frame() != b.frame()
                || !same_material(a.material(), b.material())
        }) {
            return Err(invalid(
                "accepted spatial circuit sources must match the actual thermal state",
            ));
        }
        let (thermal, circuit, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |state, _remaining| {
                check(&cancel)?;
                let coefficients = prepared
                    .resolve(state, || (*cancel.borrow_mut())())
                    .map_err(material)?;
                let sources = self.materials_at(problem, &coefficients.source, &routes, &cancel)?;
                let next =
                    PlanarConduction::new(problem.mesh().clone(), sources, assembly, &mut || {
                        (*cancel.borrow_mut())()
                    })
                    .map_err(field_error)?;
                let candidate = circuit.trial(
                    next,
                    target_s,
                    ports,
                    field_controls,
                    controls.circuit,
                    &cancel,
                )?;
                Ok(DrivenCoefficients {
                    circuit: candidate,
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: coefficients.source,
                    circuit_trials: 1,
                })
            },
            |candidate, _dt, checkpoint| {
                let mut powers = BTreeMap::new();
                let mut add = |cell, power| -> Result<(), Error> {
                    let previous = *powers.get(&cell).unwrap_or(&0.0);
                    let total = finite(previous + power)?;
                    if power != 0.0 && total == previous {
                        return Err(invalid(
                            "spatial circuit heat disappeared in receiving-cell sum",
                        ));
                    }
                    powers.insert(cell, total);
                    Ok(())
                };
                for (triangle, &power) in candidate.step.spatial.cell_power_w.iter().enumerate() {
                    checkpoint()?;
                    add(routes[triangle], power)?;
                }
                for heat in &candidate.step.graph_resistors {
                    checkpoint()?;
                    add(graph_routes[&heat.branch], heat.power_w)?;
                }
                Ok(powers)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(&cancel)?;
        Ok(ConductionCircuitEnthalpyCandidate {
            thermal,
            circuit,
            report,
        })
    }
}

impl<C: PlanarConductorMaterial> ConductionCircuit<C> {
    fn trial(
        &self,
        problem: PlanarConduction<C>,
        target_s: f64,
        ports: &[f64],
        field_controls: ConductionControls,
        controls: MaterialCircuitSolve,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<CircuitFieldCandidate<C>, Error> {
        check(cancel)?;
        let dt = finite(target_s - self.time_s)?;
        let resistance = problem
            .resistance(
                &self.terminal.positive,
                &self.terminal.negative,
                field_controls,
                &mut || (*cancel.borrow_mut())(),
            )
            .map_err(field_error)?;
        let ohms = resistance.resistance_ohm();
        let graph = stamp(&self.graph, &self.terminal, ohms)?;
        let dae = assemble_circuit(&graph).map_err(electrical)?;
        let x0 = consistent(&dae, &self.fluxes_wb, &self.charges_c, ports)?;
        check(cancel)?;
        let (native, _) = dae.step_audited(&x0, ports, dt).map_err(electrical)?;
        check(cancel)?;
        let stored_before_j = dae.system.hamiltonian(&x0);
        let stored_after_j = dae.system.hamiltonian(&native.x);
        if native.x.iter().chain(&native.y).any(|x| !x.is_finite())
            || [
                stored_before_j,
                stored_after_j,
                native.delta_h,
                native.dissipated,
                native.supplied,
                native.solver_residual,
            ]
            .iter()
            .any(|x| !x.is_finite())
            || stored_before_j < 0.0
            || stored_after_j < 0.0
            || native.dissipated < 0.0
        {
            return Err(invalid(
                "native spatial circuit storage/loss is unrepresentable",
            ));
        }
        let voltage = |node| {
            if node == 0 {
                0.0
            } else {
                let i = dae.node_potential_index[node - 1];
                f64::midpoint(x0[i], native.x[i])
            }
        };
        let (a, b, _) = graph.branches[self.terminal.branch];
        let spatial_voltage_v = finite(voltage(a) - voltage(b))?;
        let spatial_current_a = finite(spatial_voltage_v / ohms)?;
        if spatial_voltage_v != 0.0 && spatial_current_a == 0.0 {
            return Err(invalid("spatial circuit current underflows"));
        }
        let field = resistance
            .field_at_voltage(
                spatial_voltage_v,
                field_controls.current_tolerance_a,
                field_controls.power_tolerance_w,
                &mut || (*cancel.borrow_mut())(),
            )
            .map_err(field_error)?;
        let mut actual_current = 0.0;
        for &vertex in &self.terminal.positive {
            check(cancel)?;
            actual_current = finite(
                actual_current
                    + field.nodal_boundary_current_a()[vertex].ok_or_else(|| {
                        invalid("spatial circuit electrode lost its native reaction")
                    })?,
            )?;
        }
        let power = field.power_w();
        let terminal_power = finite(spatial_voltage_v * spatial_current_a)?;
        if finite(actual_current - spatial_current_a)?.abs() > field_controls.current_tolerance_a
            || finite(power - terminal_power)?.abs() > controls.power_tolerance_w
        {
            return Err(invalid(
                "native spatial field and circuit terminal current/work disagree",
            ));
        }
        let spatial_heat = finite(dt * power)?;
        if power > 0.0 && spatial_heat == 0.0 {
            return Err(invalid("spatial circuit heat underflows"));
        }
        let spatial = ConductionThermalStep {
            start_s: self.time_s,
            end_s: target_s,
            voltage_v: field.voltage_v().to_vec(),
            electric_field_v_m: field.electric_field_v_m().to_vec(),
            current_density_a_m2: field.current_density_a_m2().to_vec(),
            cell_power_w: field.cell_power_w().to_vec(),
            nodal_boundary_current_a: field.nodal_boundary_current_a().to_vec(),
            nodal_applied_current_a: field.nodal_applied_current_a().to_vec(),
            boundary_current_a: field.boundary_current_a(),
            maximum_free_current_residual_a: field.maximum_free_current_residual_a(),
            power_balance_residual_w: field.power_balance_residual_w(),
            heat_j: spatial_heat,
            supplied_work_j: finite(dt * finite(power - field.power_balance_residual_w())?)?,
            linear: field.report().cloned(),
        };
        let mut graph_resistors = Vec::new();
        let mut heat_j = spatial_heat;
        for (branch, &(a, b, law)) in graph.branches.iter().enumerate() {
            check(cancel)?;
            if branch == self.terminal.branch {
                continue;
            }
            if let Branch::Resistor { ohms } = law {
                let voltage_v = finite(voltage(a) - voltage(b))?;
                let current_a = finite(voltage_v / ohms)?;
                let power_w = finite(voltage_v * current_a)?;
                let q = finite(dt * power_w)?;
                if power_w < 0.0
                    || (voltage_v != 0.0 && power_w == 0.0)
                    || (power_w != 0.0 && q == 0.0)
                    || finite(current_a * ohms * current_a - power_w)?.abs()
                        > controls.power_tolerance_w
                {
                    return Err(invalid("fixed graph-resistor work/heat is unrepresentable"));
                }
                heat_j = finite(heat_j + q)?;
                graph_resistors.push(ConductionCircuitResistorHeat {
                    branch,
                    voltage_v,
                    current_a,
                    power_w,
                    heat_j: q,
                });
            }
        }
        let energy_residual_j =
            finite(stored_after_j - stored_before_j + heat_j - native.supplied)?;
        if finite(heat_j - native.dissipated)?.abs() > controls.energy_tolerance_j
            || finite((stored_after_j - stored_before_j) - native.delta_h)?.abs()
                > controls.energy_tolerance_j
            || energy_residual_j.abs() > controls.energy_tolerance_j
        {
            return Err(invalid(
                "spatial/graph heat, native storage and external circuit work do not close",
            ));
        }
        let fluxes_wb = dae.flux_index.iter().map(|&i| native.x[i]).collect();
        let charges_c = dae.charge_index.iter().map(|&i| native.x[i]).collect();
        let circuit = Self {
            problem,
            graph,
            terminal: self.terminal.clone(),
            fluxes_wb,
            charges_c,
            time_s: target_s,
        };
        let step = ConductionCircuitStep {
            native,
            stored_before_j,
            stored_after_j,
            spatial,
            spatial_voltage_v,
            spatial_current_a,
            graph_resistors,
            heat_j,
            energy_residual_j,
        };
        check(cancel)?;
        Ok(CircuitFieldCandidate { circuit, step })
    }
}
