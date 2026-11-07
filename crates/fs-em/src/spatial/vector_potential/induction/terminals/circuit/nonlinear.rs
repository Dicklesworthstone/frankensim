//! Bounded nonlinear electrode closure on the original native circuit graph.
//!
//! Unknowns are physical electrode voltage differences. Every evaluation
//! advances the same accepted field and native flux/charge state, then compares
//! native stage voltage with the field's endpoint drive. The reference current
//! V/Rref-I cancels its numerical resistor exactly at convergence. No affine
//! field response or constitutive state is reused across nonlinear trials.

use super::native::{Allowances, FieldReceipt};
use super::*;
use fs_solver::{
    NewtonKrylovConfig, NewtonKrylovState,
    finite_difference::{FallibleFiniteDifference, FiniteDifferenceError},
};
use std::cell::{Cell, RefCell};

/// Shared physical acceptance and bounded nonlinear interface work controls.
/// Native midpoint stages share electrode V/I with backward-Euler endpoint
/// fields. This coupled discretization has first-order temporal accuracy.
#[derive(Clone, Copy, Debug)]
pub struct SpatialTerminalCircuitNonlinearControls<F> {
    pub field: F,
    /// Voltage perturbation used by the existing finite-difference Jacobian [V].
    pub difference_step_v: f64,
    /// Includes every initial/Jacobian/line-search field solve and final replay.
    pub maximum_field_trials: usize,
    pub maximum_newton_iterations: usize,
    pub voltage_tolerance_v: f64,
    pub current_tolerance_a: f64,
    pub energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
}

#[derive(Clone)]
pub(in crate::spatial::vector_potential::induction) struct Circuit {
    layout: Rc<Layout>,
    state: Vec<f64>,
    graph_materials: Vec<GraphResistorMaterial>,
    /// Accepted field drive, used only as the next nonlinear initial guess.
    relative_voltage_v: Vec<f64>,
}

pub(in crate::spatial::vector_potential::induction) struct Candidate<S, F> {
    pub field: S,
    pub circuit: Circuit,
    pub step: SpatialTerminalCircuitStep<F>,
}

struct Trial<S, F> {
    field: S,
    step: F,
    native: StepRecord,
    input: Vec<f64>,
    residual: Vec<f64>,
    relative_residual: f64,
}

fn numerical(
    error: FiniteDifferenceError<SpatialTerminalCircuitError>,
) -> SpatialTerminalCircuitError {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(message) => invalid(message),
    }
}

impl Circuit {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        reference: &[usize],
        vertex_count: usize,
        graph: CircuitGraph,
        terminal_nodes: Vec<usize>,
        mutual: Vec<MutualInductance>,
        initial: SpatialTerminalCircuitInitial,
        config: SpatialTerminalCircuitConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let (layout, state) = native::prepare(
            reference,
            vertex_count,
            graph,
            terminal_nodes,
            mutual,
            initial,
            config,
            cancelled,
        )?;
        let relative_voltage_v = vec![0.0; layout.coordinates.len()];
        Ok(Self {
            layout,
            state,
            relative_voltage_v,
            graph_materials: Vec::new(),
        })
    }

    pub fn graph(&self) -> &CircuitGraph {
        &self.layout.graph
    }
    pub fn dae(&self) -> &CircuitDae {
        &self.layout.dae
    }
    pub fn mutual_inductance(&self) -> &[MutualInductance] {
        &self.layout.mutual
    }
    pub fn terminal_nodes(&self) -> &[usize] {
        &self.layout.nodes
    }
    pub fn native_state(&self) -> &[f64] {
        &self.state
    }
    pub fn graph_resistor_materials(&self) -> &[GraphResistorMaterial] {
        &self.graph_materials
    }

    pub fn with_graph_resistors(
        &self,
        mut sources: Vec<GraphResistorMaterial>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let count = self
            .layout
            .graph
            .branches
            .iter()
            .filter(|(_, _, branch)| matches!(branch, Branch::Resistor { .. }))
            .count();
        if sources.len() != count || sources.len() > MAX_TERMINAL_CIRCUIT_BRANCHES {
            return Err(invalid(
                "every physical graph resistor requires one bounded source",
            ));
        }
        sources.sort_by_key(|source| source.branch);
        if sources == self.graph_materials {
            poll(cancelled)?;
            return Ok(self.clone());
        }
        let layout = native::rebind_resistors(&self.layout, &self.state, &sources, cancelled)?;
        let next = Self {
            layout,
            state: self.state.clone(),
            graph_materials: sources,
            relative_voltage_v: self.relative_voltage_v.clone(),
        };
        poll(cancelled)?;
        Ok(next)
    }

    /// `advance` owns original source preparation and always uses the same
    /// accepted electrode field. No candidate becomes a subsequent trial's seed.
    pub fn step<S, F: FieldReceipt, FC>(
        &self,
        external: &[f64],
        duration_s: f64,
        controls: SpatialTerminalCircuitNonlinearControls<FC>,
        advance: impl Fn(&[f64], &mut dyn FnMut() -> bool) -> Result<(S, F)>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Candidate<S, F>> {
        poll(cancelled)?;
        admit_inputs(&self.layout, external)?;
        controls.newton.validate()?;
        if !positive(duration_s)
            || !positive(controls.difference_step_v)
            || controls.maximum_field_trials == 0
            || controls.maximum_newton_iterations == 0
            || !positive(controls.voltage_tolerance_v)
            || !positive(controls.current_tolerance_a)
            || !positive(controls.energy_tolerance_j)
        {
            return Err(invalid(
                "positive interval, nonlinear budgets and terminal allowances required",
            ));
        }
        let k = self.layout.coordinates.len();
        let cancel = RefCell::new(cancelled);
        let trials = Cell::new(0_usize);
        let evaluate = |relative: &[f64]| -> Result<Trial<S, F>> {
            let mut check = || (*cancel.borrow_mut())();
            poll(&mut check)?;
            if relative.len() != k || relative.iter().any(|v| !v.is_finite()) {
                return Err(invalid(
                    "complete finite physical electrode voltages required",
                ));
            }
            if trials.get() >= controls.maximum_field_trials {
                return Err(SpatialTerminalCircuitError::FieldTrialBudget);
            }
            trials.set(trials.get() + 1);
            let mut voltage = vec![0.0; self.layout.nodes.len()];
            for (&terminal, &value) in self.layout.coordinates.iter().zip(relative) {
                voltage[terminal] = value;
            }
            let (field, step) = advance(&voltage, &mut check)?;
            let account = step.account();
            if account.terminals.len() != voltage.len() {
                return Err(invalid("field trial lost its original electrode reactions"));
            }
            let mut input = external.to_vec();
            input.resize(self.layout.dae.ports, 0.0);
            for (p, &terminal) in self.layout.coordinates.iter().enumerate() {
                input[self.layout.external_ports + p] = finite(
                    relative[p] / self.layout.reference_ohm
                        - account.terminals[terminal].incoming_current_a,
                )?;
            }
            let native =
                native::advance(&self.layout, &self.state, &input, duration_s, &mut check)?;
            let mut residual = Vec::with_capacity(k);
            let mut norm = 0.0_f64;
            let mut scale = 0.0_f64;
            for (p, &required) in relative.iter().enumerate() {
                let actual = finite(native.y[self.layout.external_ports + p])?;
                let defect = finite(actual - required)?;
                residual.push(defect);
                norm = finite(norm.hypot(defect))?;
                scale = finite(scale.hypot(actual.abs().max(required.abs())))?;
            }
            let relative_residual = finite(norm / scale.max(controls.voltage_tolerance_v))?;
            Ok(Trial {
                field,
                step,
                native,
                input,
                residual,
                relative_residual,
            })
        };
        let (relative, report) = if k == 0 {
            (Vec::new(), None)
        } else {
            let interface =
                FallibleFiniteDifference::new(k, controls.difference_step_v, |x: &[f64]| {
                    evaluate(x).map(|trial| trial.residual)
                })
                .map_err(numerical)?;
            // A k-dimensional interface cannot gain additional independent
            // Krylov directions. Keep its workspace within preparation's O(k²)
            // bound even when the shared solver policy requests a larger restart.
            let mut newton = controls.newton;
            newton.linear_restart = newton.linear_restart.min(k);
            let construction =
                NewtonKrylovState::new(&interface, self.relative_voltage_v.clone(), newton);
            if let Some(error) = interface.take_failure() {
                return Err(numerical(error));
            }
            let mut solver = construction?;
            let result =
                solver.run_controlled(&interface, controls.maximum_newton_iterations, &mut || {
                    (*cancel.borrow_mut())() || interface.has_failed()
                });
            if let Some(error) = interface.take_failure() {
                return Err(numerical(error));
            }
            poll(&mut || (*cancel.borrow_mut())())?;
            if result.progress == SolverRunProgress::Paused {
                return Err(SpatialTerminalCircuitError::Cancelled);
            }
            if result.progress != SolverRunProgress::Complete {
                return Err(SpatialTerminalCircuitError::Interrupted(result.progress));
            }
            if !result.report.converged {
                return Err(SpatialTerminalCircuitError::NotConverged(result.report));
            }
            (solver.x, Some(result.report))
        };
        // Final independent replay, including all original native and physical
        // field gates, cannot be replaced by a cached residual or trial receipt.
        let trial = evaluate(&relative)?;
        for &defect in &trial.residual {
            close(
                defect,
                controls.voltage_tolerance_v,
                "actual nonlinear electrode voltage does not close",
            )?;
        }
        let mut check = || (*cancel.borrow_mut())();
        let mut step = native::account(
            &self.layout,
            &self.state,
            trial.native,
            trial.step,
            &self.graph_materials,
            &trial.input,
            external,
            duration_s,
            Allowances {
                voltage_v: controls.voltage_tolerance_v,
                current_a: controls.current_tolerance_a,
                energy_j: controls.energy_tolerance_j,
            },
            &mut check,
        )?;
        step.field_trials = trials.get();
        step.relative_terminal_residual = trial.relative_residual;
        step.report = report;
        let circuit = Self {
            layout: Rc::clone(&self.layout),
            state: step.native.x.clone(),
            graph_materials: self.graph_materials.clone(),
            relative_voltage_v: relative,
        };
        poll(&mut check)?;
        Ok(Candidate {
            field: trial.field,
            circuit,
            step,
        })
    }
}
