//! Retained polarization switching on the native Kirchhoff circuit timeline.
//!
//! One original capacitor becomes a private voltage port. Every nonlinear
//! trial reuses the SAME accepted material state and supplies the exact segment
//! mean voltage; charge closes `delta q = dt I` through the native DAE. Private
//! port work cancels the dielectric's stored-energy change and switching heat.
//! Ordinary R/L/C and external sources keep their original solver and accounts.
//!
//! Each interval assumes monotone charge. Supply source discontinuities and
//! turning points explicitly. Three instantaneous current probes reject turns
//! they detect but cannot certify the absence of unresolved interior turns.
//! Linear storage recovers midpoint stepping; smooth fixed active sets have
//! second-order time accuracy. Switching/saturation crossings need refinement.
//! Ideal voltage constraints can refuse native initialization or stepping.
//! Parameters and geometry stay fixed. The returned switching heat is owed
//! once to a separate thermal owner; this module does not itself advance heat,
//! add leakage, infer electromechanical strain or change source coefficients.

use super::{
    MaterialCircuit, MaterialCircuitSolve, MaterialCircuitStep, MaterialResistor,
    account_step_with_storage_change, bind, extended,
    hysteretic::{HystereticCircuitError, finite, invalid, numerical},
};
use fs_material::field_storage::dielectric::hysteresis::{
    DielectricHysteresisState, HystereticCapacitorStep, UniformHystereticCapacitor,
};
use fs_phs::{
    StepRecord,
    circuit::{Branch, CircuitDae, CircuitGraph},
};
use fs_solver::{
    NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
    finite_difference::FallibleFiniteDifference,
};
use std::cell::{Cell, RefCell};

/// Finite enthalpy, source-resolved resistance and retained switching memory.
pub mod thermal;

type Result<T> = core::result::Result<T, HystereticCircuitError>;

/// Explicit prepared polarization plus all other native conserved coordinates.
pub struct DielectricCircuitInitial {
    pub dielectric: DielectricHysteresisState,
    /// Every original inductor, in original branch order [Wb].
    pub fluxes_wb: Vec<f64>,
    /// Unbound original capacitors, in original branch order [C].
    pub remaining_charges_c: Vec<f64>,
    /// Original external source ports, preserving authored index gaps.
    pub external_inputs: Vec<f64>,
    pub time_s: f64,
}

/// Bounded native Newton work with independent physical acceptance tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DielectricCircuitSolve {
    pub maximum_circuit_trials: usize,
    pub maximum_newton_iterations: usize,
    /// Charge-increment Jacobian stencil [C].
    pub difference_step_c: f64,
    pub terminal_current_tolerance_a: f64,
    pub energy_tolerance_j: f64,
    pub circuit: MaterialCircuitSolve,
    /// Existing Newton policy, on a scalar charge residual [C].
    pub newton: NewtonKrylovConfig,
}

/// Accepted polarization and native circuit history. Trials never mutate it.
#[derive(Clone, Debug, PartialEq)]
pub struct DielectricHystereticCircuit {
    circuit: MaterialCircuit,
    capacitor: UniformHystereticCapacitor,
    state: DielectricHysteresisState,
    capacitor_branch: usize,
    external_ports: usize,
    charge_c: f64,
    voltage_v: f64,
    dielectric_energy_j: f64,
    time_s: f64,
}

/// One closed physical work account; private-port work is not external supply.
pub struct DielectricCircuitStep {
    /// Native account includes the PRIVATE port. Use supplied_work_j below for
    /// combined external source work and heat_j for combined physical heating.
    pub circuit: MaterialCircuitStep,
    pub dielectric: HystereticCapacitorStep,
    pub terminal_current_a: f64,
    pub terminal_current_residual_a: f64,
    pub current_start_mid_end_a: [f64; 3],
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub stored_change_j: f64,
    /// Original resistor loss plus polarization switching, each exactly once.
    pub heat_j: f64,
    pub supplied_work_j: f64,
    pub numerical_dissipation_j: f64,
    pub energy_residual_j: f64,
    pub circuit_trials: usize,
    pub nonlinear: NewtonReport,
}

pub struct DielectricCircuitCandidate {
    circuit: DielectricHystereticCircuit,
    step: DielectricCircuitStep,
}
impl DielectricCircuitCandidate {
    pub const fn circuit(&self) -> &DielectricHystereticCircuit {
        &self.circuit
    }
    pub const fn step(&self) -> &DielectricCircuitStep {
        &self.step
    }
    pub fn into_parts(self) -> (DielectricHystereticCircuit, DielectricCircuitStep) {
        (self.circuit, self.step)
    }
}

impl DielectricHystereticCircuit {
    /// Replace one original numeric capacitor placeholder by the explicit fit.
    /// All resistors still require their original material sources; initial
    /// polarization is caller-prepared and is never inferred from charge alone.
    pub fn new(
        mut graph: CircuitGraph,
        capacitor_branch: usize,
        capacitor: UniformHystereticCapacitor,
        resistors: Vec<MaterialResistor>,
        initial: DielectricCircuitInitial,
    ) -> Result<Self> {
        super::admit_graph(&graph)?;
        if !initial.time_s.is_finite()
            || initial.time_s < 0.0
            || !graph
                .branches
                .get(capacitor_branch)
                .is_some_and(|b| matches!(b.2, Branch::Capacitor { .. }))
        {
            return Err(invalid(
                "one original capacitor and a finite nonnegative initial clock required",
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
                "complete external inputs and room for one private dielectric port required",
            ));
        }
        let held = capacitor.hold_candidate(&initial.dielectric)?;
        graph.branches[capacitor_branch].2 = Branch::VoltageSource {
            port: external_ports,
        };
        let mut inputs = initial.external_inputs;
        inputs.push(held.voltage_v);
        let circuit = MaterialCircuit::new(
            graph,
            resistors,
            initial.fluxes_wb,
            initial.remaining_charges_c,
            &inputs,
        )?;
        Ok(Self {
            circuit,
            capacitor,
            state: held.material.state,
            capacitor_branch,
            external_ports,
            charge_c: held.charge_c,
            voltage_v: held.voltage_v,
            dielectric_energy_j: held.stored_energy_j,
            time_s: initial.time_s,
        })
    }

    /// Native graph view contains the private voltage port in place of the C.
    pub const fn native_circuit(&self) -> &MaterialCircuit {
        &self.circuit
    }
    pub const fn capacitor(&self) -> &UniformHystereticCapacitor {
        &self.capacitor
    }
    pub const fn dielectric_state(&self) -> &DielectricHysteresisState {
        &self.state
    }
    pub const fn capacitor_branch(&self) -> usize {
        self.capacitor_branch
    }
    pub const fn charge_c(&self) -> f64 {
        self.charge_c
    }
    pub const fn voltage_v(&self) -> f64 {
        self.voltage_v
    }
    pub const fn dielectric_stored_energy_j(&self) -> f64 {
        self.dielectric_energy_j
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }

    /// Advance one held-source, monotone-charge interval atomically. Newton
    /// probes, refused steps and cancellation never become material history.
    pub fn step_candidate(
        &self,
        target_s: f64,
        resistors: Vec<MaterialResistor>,
        external_inputs: &[f64],
        controls: DielectricCircuitSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DielectricCircuitCandidate> {
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
                controls.difference_step_c,
                controls.terminal_current_tolerance_a,
                controls.energy_tolerance_j,
                controls.circuit.energy_tolerance_j,
                controls.circuit.power_tolerance_w,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
            || external_inputs.len() != self.external_ports
            || external_inputs.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "positive interval, budgets, physical tolerances and complete finite ports required",
            ));
        }
        let (graph, resistors, dae) = bind(&self.circuit.graph, resistors, &[])?;
        let mut inputs = external_inputs.to_vec();
        inputs.push(self.voltage_v);
        let x0 = extended(
            &dae,
            &self.circuit.fluxes_wb,
            &self.circuit.charges_c,
            &inputs,
        )?;
        let before = finite(dae.system.hamiltonian(&x0))?;
        let trials = Cell::new(0usize);
        let reserve = || {
            check()?;
            if trials.get() >= controls.maximum_circuit_trials {
                return Err(HystereticCircuitError::TrialBudget);
            }
            trials.set(trials.get() + 1);
            Ok(())
        };
        let evaluate = |delta: f64| -> Result<(HystereticCapacitorStep, StepRecord, f64)> {
            reserve()?;
            let charge = finite(self.charge_c + delta)?;
            if delta != 0.0 && charge == self.charge_c {
                return Err(invalid(
                    "dielectric charge increment is below representable resolution",
                ));
            }
            let dielectric = self.capacitor.candidate_charge(&self.state, charge)?;
            let mut input = external_inputs.to_vec();
            input.push(dielectric.mean_voltage_v);
            let (native, _) = dae.step_audited(&x0, &input, dt)?;
            // Native voltage-port output is current supplied INTO the network.
            // Capacitor current is the opposite, into its positive terminal.
            let residual =
                finite(dielectric.charge_c - self.charge_c + dt * native.y[self.external_ports])?;
            check()?;
            Ok((dielectric, native, residual))
        };
        let problem =
            FallibleFiniteDifference::new(1, controls.difference_step_c, |x: &[f64]| {
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
        let (dielectric, native, residual) = evaluate(solver.x[0])?;
        let current = -native.y[self.external_ports];
        let current_residual = finite(residual / dt)?;
        if current_residual.abs() > controls.terminal_current_tolerance_a {
            return Err(HystereticCircuitError::Balance(
                "native capacitor current and retained charge increment disagree",
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
        let middle = self.capacitor.candidate_charge(
            &self.state,
            f64::midpoint(self.charge_c, dielectric.charge_c),
        )?;
        let mut currents = [0.0; 3];
        for (index, (fluxes, charges, voltage)) in [
            (
                &self.circuit.fluxes_wb,
                &self.circuit.charges_c,
                self.voltage_v,
            ),
            (&middle_fluxes, &middle_charges, middle.voltage_v),
            (&next_fluxes, &next_charges, dielectric.voltage_v),
        ]
        .into_iter()
        .enumerate()
        {
            reserve()?;
            currents[index] =
                instantaneous_current(&dae, fluxes, charges, external_inputs, voltage)?;
        }
        let change = dielectric.charge_c - self.charge_c;
        if currents.iter().any(|&i| {
            if change == 0.0 {
                i.abs() > controls.terminal_current_tolerance_a
            } else {
                i * change.signum() < -controls.terminal_current_tolerance_a
            }
        }) {
            return Err(HystereticCircuitError::ElectricReversal {
                charge_change_c: change,
                current_start_mid_end_a: currents,
            });
        }
        let after = finite(dae.system.hamiltonian(&native.x))?;
        let (resistor_heat, contacts, graph_heat, graph_residual) =
            account_step_with_storage_change(
                &graph,
                &resistors,
                &[],
                &dae.node_potential_index,
                &x0,
                &native,
                before,
                after,
                native.delta_h,
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
        let stored_change = finite(
            native.delta_h + self.capacitor.volume_m3() * dielectric.material.energy_change_j_m3,
        )?;
        let heat = finite(graph_heat + dielectric.switching_heat_j)?;
        let balance = finite(stored_change + heat - supplied)?;
        if balance.abs() > controls.energy_tolerance_j
            || (supplied - native.supplied - dielectric.electric_work_j).abs()
                > controls.energy_tolerance_j
        {
            return Err(HystereticCircuitError::Balance(
                "native external work, dielectric storage and distinct physical losses do not close",
            ));
        }
        let step = DielectricCircuitStep {
            terminal_current_a: current,
            terminal_current_residual_a: current_residual,
            current_start_mid_end_a: currents,
            stored_before_j: finite(before + self.dielectric_energy_j)?,
            stored_after_j: finite(after + dielectric.stored_energy_j)?,
            stored_change_j: stored_change,
            heat_j: heat,
            supplied_work_j: supplied,
            numerical_dissipation_j: 0.0,
            energy_residual_j: balance,
            circuit_trials: trials.get(),
            nonlinear: result.report,
            circuit: MaterialCircuitStep {
                native,
                stored_before_j: before,
                stored_after_j: after,
                resistors: resistor_heat,
                contacts,
                heat_j: graph_heat,
                energy_residual_j: graph_residual,
            },
            dielectric: dielectric.clone(),
        };
        let mut circuit = self.circuit.clone();
        circuit.graph = graph;
        circuit.resistors = resistors;
        circuit.fluxes_wb = next_fluxes;
        circuit.charges_c = next_charges;
        let next = Self {
            circuit,
            capacitor: self.capacitor.clone(),
            state: dielectric.material.state,
            capacitor_branch: self.capacitor_branch,
            external_ports: self.external_ports,
            charge_c: dielectric.charge_c,
            voltage_v: dielectric.voltage_v,
            dielectric_energy_j: dielectric.stored_energy_j,
            time_s: target_s,
        };
        check()?;
        Ok(DielectricCircuitCandidate {
            circuit: next,
            step,
        })
    }
}

fn instantaneous_current(
    dae: &CircuitDae,
    fluxes: &[f64],
    charges: &[f64],
    external: &[f64],
    voltage: f64,
) -> Result<f64> {
    let mut input = external.to_vec();
    input.push(voltage);
    let seed = extended(dae, fluxes, charges, &input)?;
    let solved = dae.consistent_initial_state(&seed, &input)?;
    // dt=0 returns an algebraic endpoint whose MIDPOINT is the physical effort.
    let midpoint: Vec<_> = seed
        .iter()
        .zip(&solved)
        .map(|(&a, &b)| f64::midpoint(a, b))
        .collect();
    finite(-dae.system.output(&midpoint)[external.len()])
}
