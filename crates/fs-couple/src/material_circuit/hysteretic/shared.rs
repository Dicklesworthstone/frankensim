//! A source-owned common magnetic history connected to multiple native ports.
//!
//! All selected winding branches become private current ports. The original
//! shared material owner supplies endpoint currents and exact means along the
//! total-linkage chord, including physical leakage, core memory and pinning
//! heat. One vector closure enforces each winding's Faraday equation through
//! the existing native circuit stepper. No independent per-winding hysteresis
//! loop, fitted loss resistance or history-dependent Storage fiction is added.
//!
//! Removing the complete selected magnetic block matters: leaving leakage in
//! the native midpoint DAE while driving only a core's segment mean current
//! violates endpoint current constraints at pinning or saturation crossings.
//! Here each trial freezes the SAME accepted shared state; endpoint currents
//! satisfy the actual common-core MMF after every accepted threshold crossing.
//!
//! Geometry and source coefficients remain fixed. Each interval assumes a
//! monotone common-core flux, with three native instantaneous probes rejecting
//! detectable interior turns. Individual winding linkages may run in different
//! directions. Refinement and caller-supplied waveform/reversal cuts remain
//! necessary; three samples do not certify an unresolved oscillation.

pub mod thermal;

use super::{HystereticCircuitError, HystereticCircuitSolve, finite, invalid, numerical};
use crate::material_circuit::{
    MAX_MATERIAL_CIRCUIT_BRANCHES, MAX_MATERIAL_CIRCUIT_NODES, MAX_PORTS, MAX_TRANSFORMERS,
    MaterialCircuit, MaterialCircuitStep, MaterialResistor, account_step, bind, extended,
};
use fs_matdb::QueryPoint;
use fs_material::field_storage::hysteresis::shared::{
    SharedHystereticCore, SharedHystereticState, SharedHystereticStep,
};
use fs_phs::{
    StepRecord,
    circuit::{Branch, CircuitDae, CircuitGraph},
};
use fs_solver::{
    NewtonKrylovState, NewtonReport, SolverRunProgress, finite_difference::FallibleFiniteDifference,
};
use std::cell::{Cell, RefCell};

type Result<T> = core::result::Result<T, HystereticCircuitError>;

/// Actual shared core memory and all other original conserved circuit states.
pub struct SharedHystereticCircuitInitial {
    /// Total linkages in the declared winding order and ONE magnetic history.
    pub magnetic: SharedHystereticState,
    /// Unselected original inductors, in original branch order [Wb].
    pub remaining_fluxes_wb: Vec<f64>,
    /// All original capacitors, in original branch order [C].
    pub charges_c: Vec<f64>,
    /// Original external source ports, including authored index gaps.
    pub external_inputs: Vec<f64>,
    pub time_s: f64,
}

/// Accepted multiple-winding history with the original circuit identities.
#[derive(Clone, Debug, PartialEq)]
pub struct SharedHystereticCircuit {
    circuit: MaterialCircuit,
    core: SharedHystereticCore,
    state: SharedHystereticState,
    branch_indices: Vec<usize>,
    external_ports: usize,
    reference_flux_linkage_wb: f64,
    time_s: f64,
}

/// One shared magnetic receipt and native electrical work, with no duplicated
/// core stores or losses. Winding arrays use `branch_indices()` order.
pub struct SharedHystereticCircuitStep {
    /// Includes PRIVATE winding ports; the combined external work is supplied_work_j.
    pub circuit: MaterialCircuitStep,
    pub magnetic: SharedHystereticStep,
    /// Native work-conjugate winding terminal voltages [V].
    pub terminal_voltages_v: Vec<f64>,
    pub terminal_voltage_residuals_v: Vec<f64>,
    /// Sum(turns_ratio * terminal_voltage / leakage) at instantaneous probes.
    /// Its sign is the common reference flux rate's sign, by strict monotonicity.
    pub drive_rate_start_mid_end_a_s: [f64; 3],
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub stored_change_j: f64,
    /// Original graph resistor heat plus the ONE shared core's pinning heat.
    pub heat_j: f64,
    pub supplied_work_j: f64,
    /// An accounting defect is never relabeled as physical heat.
    pub numerical_dissipation_j: f64,
    pub energy_residual_j: f64,
    pub circuit_trials: usize,
    pub nonlinear: NewtonReport,
}

pub struct SharedHystereticCircuitCandidate {
    circuit: SharedHystereticCircuit,
    step: SharedHystereticCircuitStep,
}
impl SharedHystereticCircuitCandidate {
    pub const fn circuit(&self) -> &SharedHystereticCircuit {
        &self.circuit
    }
    pub const fn step(&self) -> &SharedHystereticCircuitStep {
        &self.step
    }
    pub fn into_parts(self) -> (SharedHystereticCircuit, SharedHystereticCircuitStep) {
        (self.circuit, self.step)
    }
}

impl SharedHystereticCircuit {
    /// Bind each common-core winding to one original inductor branch, in the
    /// SAME order as `core.windings()`. Selected numeric L values are replaced;
    /// all explicit physical leakage belongs to the source-owned shared core.
    /// Original resistor/source branches, node identities and other L/C remain
    /// native. No virgin or remanent state is guessed. Unsupported ideal-current
    /// cutsets refuse through the original consistent-initial-state owner.
    pub fn new(
        mut graph: CircuitGraph,
        branch_indices: Vec<usize>,
        core: SharedHystereticCore,
        resistors: Vec<MaterialResistor>,
        initial: SharedHystereticCircuitInitial,
    ) -> Result<Self> {
        // Capacity/indices precede both graph cloning and source/numerical work.
        if graph.node_count > MAX_MATERIAL_CIRCUIT_NODES
            || graph.branches.len() > MAX_MATERIAL_CIRCUIT_BRANCHES
            || graph.transformers.len() > MAX_TRANSFORMERS
            || branch_indices.len() != core.windings().len()
            || branch_indices.is_empty()
            || branch_indices.len() > MAX_PORTS
            || !initial.time_s.is_finite()
            || initial.time_s < 0.0
        {
            return Err(invalid(
                "bounded complete shared windings and a nonnegative clock required",
            ));
        }
        let mut external_ports = 0;
        for (_, _, branch) in &graph.branches {
            if let Branch::VoltageSource { port } | Branch::CurrentSource { port } = branch {
                if *port >= MAX_PORTS {
                    return Err(invalid("original shared-core source port exceeds capacity"));
                }
                external_ports = external_ports.max(*port + 1);
            }
        }
        if external_ports + branch_indices.len() > MAX_PORTS
            || initial.external_inputs.len() != external_ports
        {
            return Err(invalid(
                "complete external inputs and private winding-port capacity required",
            ));
        }
        for (winding, &branch) in branch_indices.iter().enumerate() {
            if branch_indices[..winding].contains(&branch)
                || !graph
                    .branches
                    .get(branch)
                    .is_some_and(|entry| matches!(entry.2, Branch::Inductor { .. }))
            {
                return Err(invalid(
                    "shared windings require distinct original inductor branches",
                ));
            }
        }
        let held = core.hold_candidate(&initial.magnetic)?;
        for (winding, &branch) in branch_indices.iter().enumerate() {
            graph.branches[branch].2 = Branch::CurrentSource {
                port: external_ports + winding,
            };
        }
        let mut inputs = initial.external_inputs;
        inputs.extend(held.state.currents_a().iter().map(|&current| -current));
        let circuit = MaterialCircuit::new(
            graph,
            resistors,
            initial.remaining_fluxes_wb,
            initial.charges_c,
            &inputs,
        )?;
        Ok(Self {
            circuit,
            core,
            state: held.state,
            branch_indices,
            external_ports,
            reference_flux_linkage_wb: held.core.flux_linkage_wb,
            time_s: initial.time_s,
        })
    }

    /// Original topology; the selected branches are private winding-current ports.
    pub const fn native_circuit(&self) -> &MaterialCircuit {
        &self.circuit
    }
    pub const fn core(&self) -> &SharedHystereticCore {
        &self.core
    }
    pub const fn magnetic_state(&self) -> &SharedHystereticState {
        &self.state
    }
    pub fn branch_indices(&self) -> &[usize] {
        &self.branch_indices
    }
    /// Linkage in the source core's declared reference winding [Wb].
    pub const fn reference_flux_linkage_wb(&self) -> f64 {
        self.reference_flux_linkage_wb
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }

    /// Recheck the original constant source fit at another supported context,
    /// preserving every winding linkage, shared polarization and native state.
    /// A thermal owner must own the actual temperature and publication.
    pub fn with_core_source_query(&self, query: &QueryPoint) -> Result<Self> {
        let core = self.core.at_source_query(query)?;
        core.hold_candidate(&self.state)?;
        Ok(Self {
            core,
            ..self.clone()
        })
    }

    /// Advance one monotone COMMON-flux interval with held external ports.
    /// Exact source pinning/saturation crossings are admitted inside this
    /// interval; trials never advance the accepted magnetic memory.
    pub fn step_candidate(
        &self,
        target_s: f64,
        resistors: Vec<MaterialResistor>,
        external_inputs: &[f64],
        controls: HystereticCircuitSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SharedHystereticCircuitCandidate> {
        let cancel = RefCell::new(cancelled);
        let check = || {
            if (*cancel.borrow_mut())() {
                Err(HystereticCircuitError::Cancelled)
            } else {
                Ok(())
            }
        };
        check()?;
        controls.newton.validate()?;
        let dt = target_s - self.time_s;
        if !target_s.is_finite()
            || !dt.is_finite()
            || dt <= 0.0
            || controls.maximum_circuit_trials < 4
            || controls.maximum_newton_iterations == 0
            || [
                controls.difference_step_wb,
                controls.terminal_voltage_tolerance_v,
                controls.energy_tolerance_j,
                controls.circuit.energy_tolerance_j,
                controls.circuit.power_tolerance_w,
            ]
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
            || external_inputs.len() != self.external_ports
            || external_inputs.iter().any(|value| !value.is_finite())
        {
            return Err(invalid(
                "positive interval, budgets and complete finite shared-core ports required",
            ));
        }
        let (graph, resistors, dae) = bind(&self.circuit.graph, resistors, &[])?;
        let mut inputs = external_inputs.to_vec();
        inputs.extend(self.state.currents_a().iter().map(|&current| -current));
        let x0 = extended(
            &dae,
            &self.circuit.fluxes_wb,
            &self.circuit.charges_c,
            &inputs,
        )?;
        let circuit_before = finite(dae.system.hamiltonian(&x0))?;
        let trials = Cell::new(0usize);
        let reserve = || {
            check()?;
            if trials.get() >= controls.maximum_circuit_trials {
                return Err(HystereticCircuitError::TrialBudget);
            }
            trials.set(trials.get() + 1);
            Ok(())
        };
        let count = self.branch_indices.len();
        let evaluate = |delta: &[f64]| -> Result<(SharedHystereticStep, StepRecord, Vec<f64>)> {
            reserve()?;
            if delta.len() != count {
                return Err(invalid("shared winding trial has the wrong dimension"));
            }
            let linkages = self
                .state
                .flux_linkages_wb()
                .iter()
                .zip(delta)
                .map(|(&old, &change)| {
                    let linkage = finite(old + change)?;
                    if change != 0.0 && linkage == old {
                        return Err(invalid(
                            "shared linkage increment is below representable resolution",
                        ));
                    }
                    Ok(linkage)
                })
                .collect::<Result<Vec<_>>>()?;
            let magnetic = self.core.candidate_linkages(&self.state, &linkages)?;
            let mut input = external_inputs.to_vec();
            input.extend(magnetic.mean_currents_a.iter().map(|&current| -current));
            let (native, _) = dae.step_audited(&x0, &input, dt)?;
            let residuals = magnetic
                .state
                .flux_linkages_wb()
                .iter()
                .zip(self.state.flux_linkages_wb())
                .enumerate()
                .map(|(i, (&after, &before))| {
                    finite((after - before) - dt * native.y[self.external_ports + i])
                })
                .collect::<Result<Vec<_>>>()?;
            check()?;
            Ok((magnetic, native, residuals))
        };
        let problem =
            FallibleFiniteDifference::new(count, controls.difference_step_wb, |delta: &[f64]| {
                evaluate(delta).map(|(_, _, residuals)| residuals)
            })
            .map_err(numerical)?;
        let construction = NewtonKrylovState::new(&problem, vec![0.0; count], controls.newton);
        if let Some(error) = problem.take_failure() {
            return Err(numerical(error));
        }
        let mut solver = construction?;
        let result =
            solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
                (*cancel.borrow_mut())() || problem.has_failed()
            });
        if let Some(error) = problem.take_failure() {
            return Err(numerical(error));
        }
        check()?;
        if result.progress != SolverRunProgress::Complete {
            return Err(HystereticCircuitError::Interrupted(result.progress));
        }
        if !result.report.converged {
            return Err(HystereticCircuitError::NotConverged(result.report));
        }
        let (magnetic, native, residuals) = evaluate(&solver.x)?;
        let terminal_voltages_v = native.y[self.external_ports..].to_vec();
        let terminal_voltage_residuals_v = residuals
            .into_iter()
            .map(|residual| finite(residual / dt))
            .collect::<Result<Vec<_>>>()?;
        if terminal_voltage_residuals_v
            .iter()
            .any(|residual| residual.abs() > controls.terminal_voltage_tolerance_v)
        {
            return Err(HystereticCircuitError::Balance(
                "native winding voltages and shared total-linkage increments disagree",
            ));
        }
        let next_fluxes: Vec<_> = dae.flux_index.iter().map(|&i| native.x[i]).collect();
        let next_charges: Vec<_> = dae.charge_index.iter().map(|&i| native.x[i]).collect();
        let middle_fluxes: Vec<_> = self
            .circuit
            .fluxes_wb
            .iter()
            .zip(&next_fluxes)
            .map(|(&a, &b)| f64::midpoint(a, b))
            .collect();
        let middle_charges: Vec<_> = self
            .circuit
            .charges_c
            .iter()
            .zip(&next_charges)
            .map(|(&a, &b)| f64::midpoint(a, b))
            .collect();
        let middle_linkages: Vec<_> = self
            .state
            .flux_linkages_wb()
            .iter()
            .zip(magnetic.state.flux_linkages_wb())
            .map(|(&a, &b)| f64::midpoint(a, b))
            .collect();
        let middle = self
            .core
            .candidate_linkages(&self.state, &middle_linkages)?;
        let mut rates = [0.0; 3];
        for (index, (fluxes, charges, currents)) in [
            (
                &self.circuit.fluxes_wb,
                &self.circuit.charges_c,
                self.state.currents_a(),
            ),
            (&middle_fluxes, &middle_charges, middle.state.currents_a()),
            (&next_fluxes, &next_charges, magnetic.state.currents_a()),
        ]
        .into_iter()
        .enumerate()
        {
            reserve()?;
            let voltages = instantaneous_voltages(
                &dae,
                &graph,
                &self.branch_indices,
                fluxes,
                charges,
                external_inputs,
                currents,
            )?;
            for (&voltage, winding) in voltages.iter().zip(self.core.windings()) {
                rates[index] =
                    finite(rates[index] + winding.turns_ratio * (voltage / winding.leakage_h))?;
            }
        }
        let change = finite(magnetic.core.flux_linkage_wb - self.reference_flux_linkage_wb)?;
        let rate_tolerance = finite(
            self.core
                .windings()
                .iter()
                .map(|w| w.turns_ratio.abs() / w.leakage_h * controls.terminal_voltage_tolerance_v)
                .sum::<f64>(),
        )?;
        if rates.iter().any(|&rate| {
            if change == 0.0 {
                rate.abs() > rate_tolerance
            } else {
                rate * change.signum() < -rate_tolerance
            }
        }) {
            return Err(HystereticCircuitError::SharedMagneticReversal {
                reference_linkage_change_wb: change,
                drive_rate_start_mid_end_a_s: rates,
            });
        }
        let circuit_after = finite(dae.system.hamiltonian(&native.x))?;
        let (resistor_heat, contacts, graph_heat, graph_residual) = account_step(
            &graph,
            &resistors,
            &[],
            &dae.node_potential_index,
            &x0,
            &native,
            circuit_before,
            circuit_after,
            dt,
            controls.circuit,
            &mut || (*cancel.borrow_mut())(),
        )?;
        let supplied = finite(
            dt * external_inputs
                .iter()
                .zip(&native.y)
                .map(|(&input, &output)| input * output)
                .sum::<f64>(),
        )?;
        let magnetic_work = finite(magnetic.stored_change_j + magnetic.core.core_heat_j)?;
        let stored_change = finite(native.delta_h + magnetic.stored_change_j)?;
        let total_heat = finite(graph_heat + magnetic.core.core_heat_j)?;
        let balance = finite(stored_change + total_heat - supplied)?;
        if balance.abs() > controls.energy_tolerance_j
            || (supplied - native.supplied - magnetic_work).abs() > controls.energy_tolerance_j
        {
            return Err(HystereticCircuitError::Balance(
                "shared-core storage, physical losses and external native work do not close",
            ));
        }
        let next = Self {
            circuit: MaterialCircuit {
                graph,
                mutual: self.circuit.mutual.clone(),
                magnetic_couplings: self.circuit.magnetic_couplings.clone(),
                capacitor_curves: self.circuit.capacitor_curves.clone(),
                resistors,
                contacts: vec![],
                fluxes_wb: next_fluxes,
                charges_c: next_charges,
                storage: self.circuit.storage.clone(),
            },
            core: self.core.clone(),
            state: magnetic.state.clone(),
            branch_indices: self.branch_indices.clone(),
            external_ports: self.external_ports,
            reference_flux_linkage_wb: magnetic.core.flux_linkage_wb,
            time_s: target_s,
        };
        let step = SharedHystereticCircuitStep {
            terminal_voltages_v,
            terminal_voltage_residuals_v,
            drive_rate_start_mid_end_a_s: rates,
            stored_before_j: finite(circuit_before + self.state.stored_energy_j())?,
            stored_after_j: finite(circuit_after + magnetic.stored_energy_j)?,
            stored_change_j: stored_change,
            heat_j: total_heat,
            supplied_work_j: supplied,
            numerical_dissipation_j: 0.0,
            energy_residual_j: balance,
            circuit_trials: trials.get(),
            nonlinear: result.report,
            circuit: MaterialCircuitStep {
                native,
                stored_before_j: circuit_before,
                stored_after_j: circuit_after,
                resistors: resistor_heat,
                contacts,
                heat_j: graph_heat,
                energy_residual_j: graph_residual,
            },
            magnetic,
        };
        check()?;
        Ok(SharedHystereticCircuitCandidate {
            circuit: next,
            step,
        })
    }
}

#[allow(clippy::too_many_arguments)]
fn instantaneous_voltages(
    dae: &CircuitDae,
    graph: &CircuitGraph,
    branches: &[usize],
    fluxes: &[f64],
    charges: &[f64],
    external: &[f64],
    currents: &[f64],
) -> Result<Vec<f64>> {
    let mut inputs = external.to_vec();
    inputs.extend(currents.iter().map(|&current| -current));
    let seed = extended(dae, fluxes, charges, &inputs)?;
    let solved = dae.consistent_initial_state(&seed, &inputs)?;
    let voltage = |node: usize| {
        if node == 0 {
            0.0
        } else {
            let index = dae.node_potential_index[node - 1];
            f64::midpoint(seed[index], solved[index])
        }
    };
    branches
        .iter()
        .map(|&branch| {
            let (a, b, _) = graph.branches[branch];
            finite(voltage(a) - voltage(b))
        })
        .collect()
}
