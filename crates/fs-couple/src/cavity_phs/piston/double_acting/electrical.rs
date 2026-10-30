//! Reciprocal electrical work and motion on the original gas/piston and circuit.
//!
//! A reserved voltage-source port represents a frozen linear moving-conductor
//! transducer: e = k*v_mid, F = -k*y, where y is the native current INTO the
//! circuit. The same k supplies force and back emf. KCL/KVL, electrical storage,
//! material resistance, both gases and momentum stay with their original owners.
//! This is a lumped frozen-field approximation, not a magnetic field solver.

/// Finite material heat and reciprocal electrical work on the same shared motion.
pub mod network;

use super::{DoubleActingGasPiston, DoublePistonInput, DoublePistonSolve, DoublePistonStep};
use crate::material_circuit::{MaterialCircuit, MaterialCircuitCandidate, MaterialCircuitError,
    MaterialCircuitSolve, MaterialCircuitStep, MaterialResistor};
use fs_phs::circuit::Branch;
use fs_solver::{NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use std::cell::{Cell, RefCell};

/// Geometry/field coupling to one existing circuit voltage-source branch.
/// Positive branch current is from the first to the second terminal; positive
/// force expands the primary gas. Reversing winding orientation reverses k.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TranslationalCircuitPort {
    /// Original graph branch index, NOT a port index or resistor-only index.
    pub branch: usize,
    /// Signed force/current = back-emf/velocity [N/A = V s/m]. Zero explicitly
    /// decouples mechanics while retaining the electrical branch as a short.
    /// Frozen over this step; no independent force and back-emf constants.
    pub factor_n_a: f64,
}

/// Immutable electrical history and a complete midpoint material evaluation.
#[derive(Clone, Copy)]
pub struct ElectricalMotionInput<'a> {
    /// Original accepted charges, fluxes and graph.
    pub circuit: &'a MaterialCircuit,
    /// Complete source-bound resistors; no material names imply resistance.
    pub resistors: &'a [MaterialResistor],
    /// All held circuit inputs. The reserved mechanical port MUST be zero here:
    /// an applied voltage belongs to a separate source, never a silent overwrite.
    pub inputs: &'a [f64],
    /// The one power-conjugate electrical/mechanical connection.
    pub port: TranslationalCircuitPort,
}

/// Nonlinear velocity residual, absolute physical gates and bounded nested work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElectricalMotionSolve {
    /// Original electrical storage and material branch heat tolerances.
    pub circuit: MaterialCircuitSolve,
    /// Original Newton--Krylov policy on one velocity unknown [m/s].
    pub newton: NewtonKrylovConfig,
    /// Maximum nonlinear iterations; inner linear work is bounded by `newton`.
    pub maximum_iterations: usize,
    /// All circuit/mechanical pairs, including stencils and fresh final admission.
    pub maximum_trials: usize,
    /// Finite-difference velocity scale [m/s] on the SAME accepted histories.
    pub difference_velocity_m_s: f64,
    /// Absolute final |trial velocity - actual midpoint velocity| [m/s].
    pub velocity_tolerance_m_s: f64,
    /// Absolute electrical/mechanical transfer and combined energy allowance [J].
    pub energy_tolerance_j: f64,
}

/// Refusals preserve the original circuit and both gas/mechanical histories.
#[derive(Debug)]
pub enum ElectricalMotionError<E> {
    /// Invalid branch ownership, controls or unrepresentable physical transfer.
    Invalid(&'static str),
    /// Original circuit/material refusal.
    Circuit(MaterialCircuitError),
    /// Original mechanical/thermal owner refusal.
    Mechanical(E),
    /// Original numerical construction refusal.
    Newton(NewtonError),
    /// Nonlinear residual did not converge.
    NotConverged(NewtonReport),
    /// Original numerical solve was interrupted.
    Interrupted(SolverRunProgress),
    /// Joint evaluation budget exhausted, including final admission.
    Budget,
    /// Cancellation before returning a complete candidate.
    Cancelled,
}
impl<E> From<MaterialCircuitError> for ElectricalMotionError<E> {
    fn from(e: MaterialCircuitError) -> Self {
        match e { MaterialCircuitError::Cancelled => Self::Cancelled, other => Self::Circuit(other) }
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for ElectricalMotionError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result { write!(f, "circuit/motion: {self:?}") }
}
impl<E: std::error::Error + 'static> std::error::Error for ElectricalMotionError<E> {}

/// Actual reciprocal work, using the native electrical effort/flow midpoint.
#[derive(Debug)]
pub struct ElectricalMotionStep {
    /// Voltage across the original first/second transducer terminals [V].
    pub back_emf_v: f64,
    /// Current from first to second terminal, negative native source flow [A].
    pub current_a: f64,
    /// Applied electromagnetic force, positive on primary expansion [N].
    pub force_n: f64,
    /// Force times the ACTUAL common displacement [J].
    pub mechanical_work_j: f64,
    /// Native transducer voltage times native port flow times duration [J].
    pub electrical_work_j: f64,
    /// Independently summed work of every OTHER circuit source port [J].
    pub external_electrical_work_j: f64,
    /// Trial back-emf velocity minus actual mechanical midpoint velocity [m/s].
    pub velocity_residual_m_s: f64,
    /// Electrical transducer work + mechanical transducer work [J].
    pub transfer_residual_j: f64,
    /// All circuit/mechanical proposals, including discarded stencils.
    pub trials: usize,
    /// Original nonlinear result, not a temporal-accuracy certificate.
    pub nonlinear: NewtonReport,
}

// Private composition seam shared by the ordinary and thermal-boundary adapters.
// The original mechanical owner, not this solver, computes its actual trajectory.
pub(super) struct MotionTrial<M> {
    pub state: M,
    pub displacement_m: f64,
    pub midpoint_velocity_m_s: f64,
}
pub(super) struct CoupledMotion<M> {
    pub circuit: MaterialCircuitCandidate,
    pub mechanical: M,
    pub step: ElectricalMotionStep,
}
fn positive(v: f64) -> bool { v.is_finite() && v > 0.0 }
fn numerical<E>(e: FiniteDifferenceError<ElectricalMotionError<E>>) -> ElectricalMotionError<E> {
    match e { FiniteDifferenceError::Source(e) => e,
        FiniteDifferenceError::Invalid(s) => ElectricalMotionError::Invalid(s) }
}

pub(super) fn solve<M, E, C: FnMut() -> bool>(
    input: ElectricalMotionInput<'_>, dt: f64, controls: ElectricalMotionSolve,
    initial_velocity_m_s: f64, cancelled: &RefCell<C>,
    mechanical: impl Fn(f64) -> Result<MotionTrial<M>, E>,
) -> Result<CoupledMotion<M>, ElectricalMotionError<E>> {
    let check = || if (*cancelled.borrow_mut())() { Err(ElectricalMotionError::Cancelled) } else { Ok(()) };
    check()?;
    controls.newton.validate().map_err(ElectricalMotionError::Newton)?;
    if !positive(dt) || !initial_velocity_m_s.is_finite() || !input.port.factor_n_a.is_finite()
        || !positive(controls.velocity_tolerance_m_s) || !positive(controls.energy_tolerance_j)
        || controls.maximum_iterations == 0 || controls.maximum_trials == 0 {
        return Err(ElectricalMotionError::Invalid("finite coupling/velocity, positive clock and explicit work controls required"));
    }
    let port = match input.circuit.graph().branches.get(input.port.branch) {
        Some((_, _, Branch::VoltageSource { port })) => *port,
        _ => return Err(ElectricalMotionError::Invalid("motion must own a voltage-source branch")),
    };
    if input.inputs.get(port) != Some(&0.0) || input.inputs.iter().any(|v| !v.is_finite())
        || input.circuit.graph().branches.iter().enumerate().any(|(i, (_, _, b))|
            i != input.port.branch && matches!(b, Branch::VoltageSource { port: p } | Branch::CurrentSource { port: p } if *p == port)) {
        return Err(ElectricalMotionError::Invalid("one exclusive zero-placeholder mechanical source port required"));
    }
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64]| {
        check()?;
        if trials.get() >= controls.maximum_trials { return Err(ElectricalMotionError::Budget); }
        trials.set(trials.get()+1);
        let velocity = x[0];
        let emf = input.port.factor_n_a*velocity;
        if !velocity.is_finite() || !emf.is_finite()
            || (input.port.factor_n_a != 0.0 && velocity != 0.0 && emf == 0.0) {
            return Err(ElectricalMotionError::Invalid("back emf is unrepresentable"));
        }
        let mut inputs = input.inputs.to_vec(); inputs[port] = emf;
        let circuit = input.circuit.step_candidate(input.resistors.to_vec(), &inputs, dt,
            controls.circuit, || (*cancelled.borrow_mut())())?;
        let native = &circuit.step().native;
        let current = -native.y[port];
        let force = input.port.factor_n_a*current;
        if !force.is_finite() || (input.port.factor_n_a != 0.0 && current != 0.0 && force == 0.0) {
            return Err(ElectricalMotionError::Invalid("electromagnetic force is unrepresentable"));
        }
        let motion = mechanical(force).map_err(ElectricalMotionError::Mechanical)?;
        check()?;
        let residual = velocity-motion.midpoint_velocity_m_s;
        let mechanical_work = force*motion.displacement_m;
        let electrical_work = dt*emf*native.y[port];
        let external: f64 = inputs.iter().zip(&native.y).enumerate()
            .filter(|(i, _)| *i != port).map(|(_, (u,y))| dt*u*y).sum();
        let transfer = mechanical_work+electrical_work;
        if ![residual, motion.displacement_m, motion.midpoint_velocity_m_s,
            mechanical_work, electrical_work, external, transfer].iter().all(|v| v.is_finite())
            || ((external+electrical_work)-native.supplied).abs() > controls.energy_tolerance_j {
            return Err(ElectricalMotionError::Invalid("native source work or motion is unrepresentable/inconsistent"));
        }
        Ok((circuit, motion, [emf, current, force, mechanical_work, electrical_work, external, residual, transfer]))
    };
    let problem = FallibleFiniteDifference::new(1, controls.difference_velocity_m_s,
        |x: &[f64]| evaluate(x).map(|(_,_,r)| vec![r[6]])).map_err(numerical)?;
    let construction = NewtonKrylovState::new(&problem, vec![initial_velocity_m_s], controls.newton);
    if let Some(e) = problem.take_failure() { return Err(numerical(e)); }
    let mut solver = construction.map_err(ElectricalMotionError::Newton)?;
    let result = solver.run_controlled(&problem, controls.maximum_iterations,
        &mut || (*cancelled.borrow_mut())() || problem.has_failed());
    if let Some(e) = problem.take_failure() { return Err(numerical(e)); }
    check()?;
    if result.progress != SolverRunProgress::Complete { return Err(ElectricalMotionError::Interrupted(result.progress)); }
    if !result.report.converged { return Err(ElectricalMotionError::NotConverged(result.report)); }
    let (circuit, motion, r) = evaluate(&solver.x)?;
    if r[6].abs() > controls.velocity_tolerance_m_s || r[7].abs() > controls.energy_tolerance_j
        || (motion.displacement_m-dt*motion.midpoint_velocity_m_s).abs() > dt*controls.velocity_tolerance_m_s {
        return Err(ElectricalMotionError::Invalid("actual velocity, displacement and reciprocal work did not close"));
    }
    check()?;
    Ok(CoupledMotion { circuit, mechanical: motion.state, step: ElectricalMotionStep {
        back_emf_v:r[0], current_a:r[1], force_n:r[2], mechanical_work_j:r[3], electrical_work_j:r[4],
        external_electrical_work_j:r[5], velocity_residual_m_s:r[6], transfer_residual_j:r[7],
        trials:trials.get(), nonlinear:result.report,
    } })
}

/// Complete circuit, two-gas and common mechanical endpoint with all work records.
pub struct ElectricalPistonCandidate {
    circuit: MaterialCircuitCandidate,
    piston: DoubleActingGasPiston,
    step: DoublePistonStep,
    motion: ElectricalMotionStep,
    total_energy_residual_j: f64,
}
impl ElectricalPistonCandidate {
    /// Original accepted-candidate electrical state, branch heat and source work.
    pub const fn circuit(&self) -> &MaterialCircuitCandidate { &self.circuit }
    /// Both source-backed gases and the single moving body's state.
    pub const fn piston(&self) -> &DoubleActingGasPiston { &self.piston }
    /// Mechanical report; its applied force includes the transducer force.
    pub const fn piston_step(&self) -> &DoublePistonStep { &self.step }
    /// Actual reciprocal transfer and nonlinear work.
    pub const fn motion(&self) -> &ElectricalMotionStep { &self.motion }
    /// Circuit + gas + kinetic changes + resistor heat - external supplies [J].
    pub const fn total_energy_residual_j(&self) -> f64 { self.total_energy_residual_j }
    /// Publish both histories together. Resistor heat is an output, not yet a
    /// finite thermal body update; this method must not be called twice.
    pub fn into_parts(self) -> (MaterialCircuit, DoubleActingGasPiston, MaterialCircuitStep,
        DoublePistonStep, ElectricalMotionStep, f64) {
        let (circuit, report) = self.circuit.into_parts();
        (circuit, self.piston, report, self.step, self.motion, self.total_energy_residual_j)
    }
}
impl DoubleActingGasPiston {
    /// Solve back emf, native circuit current, electromagnetic force and both
    /// gas trajectories together. `input.external_force_n` excludes the new
    /// transducer; both prescribed gas heats remain independent external inputs.
    ///
    /// The finite-difference Jacobian probes the SAME accepted electrical and
    /// gas histories; no trial consumes charge or advances a piston. No separate
    /// electrical or mechanical integrator is introduced. Missing material data,
    /// failed residual/energy admission or cancellation returns neither endpoint.
    ///
    /// Fixed C/L/topology, frozen linear force factor and uniform closed gases.
    /// No iron saturation, field history, variable inductance, eddy/skin losses,
    /// switching, thermal update or emitted sound is inferred by this adapter.
    pub fn electrical_candidate(
        &self, electrical: ElectricalMotionInput<'_>, input: DoublePistonInput,
        piston_controls: DoublePistonSolve, controls: ElectricalMotionSolve,
        initial_midpoint_velocity_m_s: f64, cancelled: impl FnMut() -> bool,
    ) -> Result<ElectricalPistonCandidate, ElectricalMotionError<super::PistonError>> {
        self.electrical_candidate_with_midpoint_resistance(electrical, input, piston_controls,
            controls, initial_midpoint_velocity_m_s, 0.0, cancelled)
    }

    // Shared by the original electrical path and reciprocal acoustic loads.
    // Only the mechanical callback sees the affine load; all circuit stencils
    // retain the original accepted charge/flux and the same held material laws.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::cavity_phs::piston) fn electrical_candidate_with_midpoint_resistance(
        &self, electrical: ElectricalMotionInput<'_>, input: DoublePistonInput,
        piston_controls: DoublePistonSolve, controls: ElectricalMotionSolve,
        initial_midpoint_velocity_m_s: f64, resistance_n_s_m: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ElectricalPistonCandidate, ElectricalMotionError<super::PistonError>> {
        if !resistance_n_s_m.is_finite() || resistance_n_s_m < 0.0 {
            return Err(ElectricalMotionError::Invalid("midpoint load resistance must be finite and nonnegative"));
        }
        let cancelled = RefCell::new(cancelled);
        let coupled = solve(electrical, input.duration_s, controls, initial_midpoint_velocity_m_s, &cancelled, |force| {
            let total_force = input.external_force_n+force;
            if !total_force.is_finite() || (force != 0.0 && total_force == input.external_force_n) {
                return Err(super::PistonError::Invalid("electromagnetic force vanished or overflowed in mechanical load"));
            }
            let mut next = self.clone();
            let step = next.step_with_midpoint_resistance(
                DoublePistonInput { external_force_n:total_force, ..input }, piston_controls,
                resistance_n_s_m, || (*cancelled.borrow_mut())())?;
            Ok(MotionTrial { displacement_m:step.displacement_m,
                midpoint_velocity_m_s:f64::midpoint(step.velocity_m_s[0], step.velocity_m_s[1]), state:(next,step) })
        });
        let coupled = match coupled {
            Err(ElectricalMotionError::Mechanical(super::PistonError::Cancelled)) => return Err(ElectricalMotionError::Cancelled),
            other => other?,
        };
        let (piston, step) = coupled.mechanical;
        let e = coupled.circuit.step();
        // The actual native force work includes both transducers. Remove only
        // the electrical work when enclosing circuit + gas + wall; the caller
        // owning a velocity-dependent mechanical load must close the remainder.
        // Keep the original arithmetic on the exact zero-resistance path.
        let external_mechanical_work_j = if resistance_n_s_m == 0.0 {
            input.external_force_n*step.displacement_m
        } else {
            step.external_force_work_j-coupled.step.mechanical_work_j
        };
        let total = (e.stored_after_j-e.stored_before_j)+e.heat_j
            +step.gas[0].stored_energy_change_j+step.gas[1].stored_energy_change_j+step.kinetic_energy_change_j
            -coupled.step.external_electrical_work_j-input.heat_j[0]-input.heat_j[1]
            -external_mechanical_work_j;
        if !external_mechanical_work_j.is_finite() || !total.is_finite() || total.abs() > controls.energy_tolerance_j {
            return Err(ElectricalMotionError::Invalid("circuit, both gases and one kinetic store do not close external work"));
        }
        if (*cancelled.borrow_mut())() { return Err(ElectricalMotionError::Cancelled); }
        Ok(ElectricalPistonCandidate { circuit:coupled.circuit, piston, step, motion:coupled.step, total_energy_residual_j:total })
    }
}

#[cfg(test)]
pub(super) mod tests;
