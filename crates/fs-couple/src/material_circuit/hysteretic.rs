//! Voltage-driven magnetic memory through a native Kirchhoff current port.
//!
//! One original inductor branch becomes a private current-source port. Every
//! Newton trial evaluates UniformHystereticCore from FROZEN old memory, supplies
//! its work-conjugate segment mean current to the native DAE, and closes
//! delta(lambda) = dt * V. Private-port work cancels the core's exact stored
//! energy increment plus pinning heat. B alone is never falsely advertised as
//! path-independent Storage; no alternate time/material integrator is introduced.
//! Ordinary R/L/C/source branches retain their native owner. Represent coil
//! resistance as a sourced graph resistor, not a hidden magnetic-fit parameter.
//!
//! Each accepted increment assumes monotone magnetic linkage. Split actual
//! field reversals and source discontinuities into separate intervals. Start,
//! midpoint and endpoint native voltage probes refuse detectable turns; they
//! cannot certify the absence of unresolved interior oscillations. The method
//! reduces to native midpoint stepping for linear storage and is second order
//! on a smooth fixed active set. Pinning/saturation crossings need timestep
//! refinement, not an exact-trajectory claim. Thermal coupling is first order.

/// Distinct core and winding heat on the original finite enthalpy timeline.
pub mod thermal;

use super::{
    MaterialCircuit, MaterialCircuitError, MaterialCircuitSolve, MaterialCircuitStep,
    MaterialResistor, account_step, bind, extended,
};
use fs_material::field_storage::{
    FieldStorageError,
    hysteresis::{HystereticCoreStep, MagneticHysteresisState, UniformHystereticCore},
};
use fs_phs::{
    StepRecord,
    circuit::{Branch, CircuitDae, CircuitError, CircuitGraph},
};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
    finite_difference::{FallibleFiniteDifference, FiniteDifferenceError},
};
use std::cell::{Cell, RefCell};

/// Actual prepared magnetic memory and the circuit's OTHER conserved coordinates.
pub struct HystereticCircuitInitial {
    pub core: MagneticHysteresisState,
    /// Unbound inductor branches, in original branch order [Wb].
    pub remaining_fluxes_wb: Vec<f64>,
    /// All original capacitor branches, in original branch order [C].
    pub charges_c: Vec<f64>,
    /// Original external ports only, including authored index gaps.
    pub external_inputs: Vec<f64>,
    pub time_s: f64,
}

/// Bounded existing Newton work and independent physical acceptance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HystereticCircuitSolve {
    /// Native evaluations, including stencils, final admission and three
    /// instantaneous voltage probes. One assembly serves all these trials.
    pub maximum_circuit_trials: usize,
    pub maximum_newton_iterations: usize,
    /// Jacobian perturbation in linkage increments [Wb].
    pub difference_step_wb: f64,
    pub terminal_voltage_tolerance_v: f64,
    pub energy_tolerance_j: f64,
    pub circuit: MaterialCircuitSolve,
    /// Shared Newton policy on a scalar linkage-increment residual [Wb].
    pub newton: NewtonKrylovConfig,
}

#[derive(Debug)]
pub enum HystereticCircuitError {
    Invalid(&'static str),
    Material(FieldStorageError),
    Circuit(MaterialCircuitError),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
    /// Split the interval at the physical magnetic turning point and retry.
    MagneticReversal {
        linkage_change_wb: f64,
        voltage_start_mid_end_v: [f64; 3],
    },
    Balance(&'static str),
    Cancelled,
}
impl core::fmt::Display for HystereticCircuitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "hysteretic circuit: {self:?}")
    }
}
impl std::error::Error for HystereticCircuitError {}
impl From<FieldStorageError> for HystereticCircuitError {
    fn from(e: FieldStorageError) -> Self {
        Self::Material(e)
    }
}
impl From<MaterialCircuitError> for HystereticCircuitError {
    fn from(e: MaterialCircuitError) -> Self {
        match e {
            MaterialCircuitError::Cancelled => Self::Cancelled,
            other => Self::Circuit(other),
        }
    }
}
impl From<CircuitError> for HystereticCircuitError {
    fn from(e: CircuitError) -> Self {
        Self::Circuit(MaterialCircuitError::Circuit(e))
    }
}
impl From<NewtonError> for HystereticCircuitError {
    fn from(e: NewtonError) -> Self {
        Self::Newton(e)
    }
}
type Result<T> = core::result::Result<T, HystereticCircuitError>;
fn invalid(s: &'static str) -> HystereticCircuitError {
    HystereticCircuitError::Invalid(s)
}
fn finite(v: f64) -> Result<f64> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(invalid("unrepresentable hysteretic circuit arithmetic"))
    }
}
fn numerical(e: FiniteDifferenceError<HystereticCircuitError>) -> HystereticCircuitError {
    match e {
        FiniteDifferenceError::Source(e) => e,
        FiniteDifferenceError::Invalid(s) => invalid(s),
    }
}

/// Immutable accepted magnetic memory and native finite circuit history.
#[derive(Clone, Debug, PartialEq)]
pub struct HystereticCircuit {
    circuit: MaterialCircuit,
    core: UniformHystereticCore,
    state: MagneticHysteresisState,
    core_branch: usize,
    external_ports: usize,
    linkage_wb: f64,
    current_a: f64,
    core_energy_j: f64,
    time_s: f64,
}

/// Native circuit work and exact irreversible core increment, kept separately.
pub struct HystereticCircuitStep {
    /// Includes the PRIVATE interface port; use supplied_work_j for the combined
    /// system's external source work.
    pub circuit: MaterialCircuitStep,
    pub core: HystereticCoreStep,
    pub terminal_voltage_v: f64,
    pub terminal_voltage_residual_v: f64,
    pub voltage_start_mid_end_v: [f64; 3],
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Stable physical change of core plus native circuit energy [J].
    pub stored_change_j: f64,
    /// Native graph I^2 R plus core pinning heat, each once [J].
    pub heat_j: f64,
    /// External source ports only; private magnetic work cancels internally.
    pub supplied_work_j: f64,
    /// No BE-like numerical loss is added. The independent residual is never
    /// converted into physical heating.
    pub numerical_dissipation_j: f64,
    pub energy_residual_j: f64,
    pub circuit_trials: usize,
    pub nonlinear: NewtonReport,
}

pub struct HystereticCircuitCandidate {
    circuit: HystereticCircuit,
    step: HystereticCircuitStep,
}
impl HystereticCircuitCandidate {
    pub const fn circuit(&self) -> &HystereticCircuit {
        &self.circuit
    }
    pub const fn step(&self) -> &HystereticCircuitStep {
        &self.step
    }
    pub fn into_parts(self) -> (HystereticCircuit, HystereticCircuitStep) {
        (self.circuit, self.step)
    }
}

impl HystereticCircuit {
    /// Replace one original inductor by its declared uniform hysteretic core.
    /// That numeric L is a placeholder; other L/C remain native stores. All
    /// resistors require material sources. No memory or linkage is guessed.
    /// Native initialization must preserve all other fluxes/charges. Ideal
    /// current-only cutsets may refuse this bounded private-port coupling.
    pub fn new(
        mut graph: CircuitGraph,
        core_branch: usize,
        core: UniformHystereticCore,
        resistors: Vec<MaterialResistor>,
        initial: HystereticCircuitInitial,
    ) -> Result<Self> {
        super::admit_graph(&graph)?;
        if !initial.time_s.is_finite()
            || initial.time_s < 0.0
            || !graph
                .branches
                .get(core_branch)
                .is_some_and(|b| matches!(b.2, Branch::Inductor { .. }))
        {
            return Err(invalid(
                "one original inductor branch and a finite nonnegative initial clock required",
            ));
        }
        let external_ports = graph
            .branches
            .iter()
            .filter_map(|b| match b.2 {
                Branch::VoltageSource { port } | Branch::CurrentSource { port } => Some(port + 1),
                _ => None,
            })
            .max()
            .unwrap_or(0);
        if external_ports >= super::MAX_PORTS || initial.external_inputs.len() != external_ports {
            return Err(invalid(
                "complete external inputs and capacity for one private core port required",
            ));
        }
        let held = core.hold_candidate(&initial.core)?;
        graph.branches[core_branch].2 = Branch::CurrentSource {
            port: external_ports,
        };
        let mut inputs = initial.external_inputs;
        inputs.push(-held.current_a);
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
            state: held.material.state,
            core_branch,
            external_ports,
            linkage_wb: held.flux_linkage_wb,
            current_a: held.current_a,
            core_energy_j: held.stored_energy_j,
            time_s: initial.time_s,
        })
    }

    /// Native graph view; the core branch is a PRIVATE current-source port.
    pub const fn native_circuit(&self) -> &MaterialCircuit {
        &self.circuit
    }
    pub const fn core(&self) -> &UniformHystereticCore {
        &self.core
    }
    pub const fn magnetic_state(&self) -> &MagneticHysteresisState {
        &self.state
    }
    pub const fn core_branch(&self) -> usize {
        self.core_branch
    }
    pub const fn flux_linkage_wb(&self) -> f64 {
        self.linkage_wb
    }
    pub const fn current_a(&self) -> f64 {
        self.current_a
    }
    pub const fn core_stored_energy_j(&self) -> f64 {
        self.core_energy_j
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }

    /// Solve one monotone-linkage interval under held external source ports.
    /// Every trial starts from identical accepted memory, fluxes and charges.
    /// Source domains, waveform cuts and unresolved turns remain the caller's
    /// responsibility; detectable turns return MagneticReversal.
    pub fn step_candidate(
        &self,
        target_s: f64,
        resistors: Vec<MaterialResistor>,
        external_inputs: &[f64],
        controls: HystereticCircuitSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<HystereticCircuitCandidate> {
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
            .any(|x| !x.is_finite() || *x <= 0.0)
            || external_inputs.len() != self.external_ports
            || external_inputs.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "positive interval, budgets, complete finite ports and physical tolerances required",
            ));
        }
        // One admitted assembly of fixed trial resistance/topology serves the
        // SAME original native stepper and final branch-work accountant.
        let (graph, resistors, dae) = bind(&self.circuit.graph, resistors, &[])?;
        let mut inputs = external_inputs.to_vec();
        inputs.push(-self.current_a);
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
        let evaluate = |delta: f64| -> Result<(HystereticCoreStep, StepRecord, f64)> {
            reserve()?;
            let linkage = finite(self.linkage_wb + delta)?;
            if delta != 0.0 && linkage == self.linkage_wb {
                return Err(invalid(
                    "magnetic linkage increment is below representable resolution",
                ));
            }
            let core = self.core.candidate_linkage(&self.state, linkage)?;
            let mut input = external_inputs.to_vec();
            input.push(-core.mean_current_a);
            let (native, _) = dae.step_audited(&x0, &input, dt)?;
            let residual = finite(
                (core.flux_linkage_wb - self.linkage_wb) - dt * native.y[self.external_ports],
            )?;
            check()?;
            Ok((core, native, residual))
        };
        let problem =
            FallibleFiniteDifference::new(1, controls.difference_step_wb, |x: &[f64]| {
                evaluate(x[0]).map(|(_, _, residual)| vec![residual])
            })
            .map_err(numerical)?;
        let construction = NewtonKrylovState::new(&problem, vec![0.0], controls.newton);
        if let Some(e) = problem.take_failure() {
            return Err(numerical(e));
        }
        let mut solver = construction?;
        let result =
            solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
                (*cancel.borrow_mut())() || problem.has_failed()
            });
        if let Some(e) = problem.take_failure() {
            return Err(numerical(e));
        }
        check()?;
        if result.progress != SolverRunProgress::Complete {
            return Err(HystereticCircuitError::Interrupted(result.progress));
        }
        if !result.report.converged {
            return Err(HystereticCircuitError::NotConverged(result.report));
        }
        let (core, native, residual) = evaluate(solver.x[0])?;
        let terminal_voltage = native.y[self.external_ports];
        let voltage_residual = finite(residual / dt)?;
        if voltage_residual.abs() > controls.terminal_voltage_tolerance_v {
            return Err(HystereticCircuitError::Balance(
                "native core voltage and magnetic linkage increment disagree",
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
        let middle_core = self.core.candidate_linkage(
            &self.state,
            f64::midpoint(self.linkage_wb, core.flux_linkage_wb),
        )?;
        let mut voltages = [0.0; 3];
        for (index, (fluxes, charges, current)) in [
            (
                &self.circuit.fluxes_wb,
                &self.circuit.charges_c,
                self.current_a,
            ),
            (&middle_fluxes, &middle_charges, middle_core.current_a),
            (&next_fluxes, &next_charges, core.current_a),
        ]
        .into_iter()
        .enumerate()
        {
            reserve()?;
            voltages[index] = instantaneous_voltage(
                &dae,
                &graph,
                self.core_branch,
                fluxes,
                charges,
                external_inputs,
                current,
            )?;
        }
        let change = core.flux_linkage_wb - self.linkage_wb;
        if voltages.iter().any(|&v| {
            if change == 0.0 {
                v.abs() > controls.terminal_voltage_tolerance_v
            } else {
                v * change.signum() < -controls.terminal_voltage_tolerance_v
            }
        }) {
            return Err(HystereticCircuitError::MagneticReversal {
                linkage_change_wb: change,
                voltage_start_mid_end_v: voltages,
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
                .map(|(u, y)| u * y)
                .sum::<f64>(),
        )?;
        let stored_change =
            finite(native.delta_h + self.core.volume_m3() * core.material.energy_change_j_m3)?;
        let total_heat = finite(graph_heat + core.core_heat_j)?;
        let balance = finite(stored_change + total_heat - supplied)?;
        if balance.abs() > controls.energy_tolerance_j
            || (supplied - native.supplied - core.magnetic_work_j).abs()
                > controls.energy_tolerance_j
        {
            return Err(HystereticCircuitError::Balance(
                "external native work, magnetic storage and distinct losses do not close",
            ));
        }
        let step = HystereticCircuitStep {
            terminal_voltage_v: terminal_voltage,
            terminal_voltage_residual_v: voltage_residual,
            voltage_start_mid_end_v: voltages,
            stored_before_j: finite(circuit_before + self.core_energy_j)?,
            stored_after_j: finite(circuit_after + core.stored_energy_j)?,
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
            core: core.clone(),
        };
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
            state: core.material.state,
            core_branch: self.core_branch,
            external_ports: self.external_ports,
            linkage_wb: core.flux_linkage_wb,
            current_a: core.current_a,
            core_energy_j: core.stored_energy_j,
            time_s: target_s,
        };
        check()?;
        Ok(HystereticCircuitCandidate {
            circuit: next,
            step,
        })
    }
}

fn instantaneous_voltage(
    dae: &CircuitDae,
    graph: &CircuitGraph,
    branch: usize,
    fluxes: &[f64],
    charges: &[f64],
    external: &[f64],
    current: f64,
) -> Result<f64> {
    let mut input = external.to_vec();
    input.push(-current);
    let seed = extended(dae, fluxes, charges, &input)?;
    let solved = dae.consistent_initial_state(&seed, &input)?;
    // Algebraic effort is the start/end MIDPOINT, including the dt=0
    // initializer. Its returned algebraic coordinate is not V itself.
    let voltage = |node: usize| {
        if node == 0 {
            0.0
        } else {
            let i = dae.node_potential_index[node - 1];
            f64::midpoint(seed[i], solved[i])
        }
    };
    let (a, b, _) = graph.branches[branch];
    finite(voltage(a) - voltage(b))
}

#[cfg(test)]
mod tests;
