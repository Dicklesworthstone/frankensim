//! One atomic piezoelectric, material-circuit and finite thermal interval.
//!
//! Original rod storage and Kirchhoff motion supply the actual resistor heat.
//! The original enthalpy residual resolves resistance at its proposed thermal
//! endpoint, so load heating feeds back into current and mechanical response.
//! The base entry holds ceramic coefficients fixed. The `thermal` entry adds
//! the same ceramic's source-derived temperature response and field entropy.
//! Neither entry infers depoling, dielectric loss or thermal strain.

use super::*;
use crate::material_circuit::piezoelectric::{
    PiezoelectricCircuit, PiezoelectricCircuitCandidate, PiezoelectricCircuitStep,
};

/// Ceramic coefficient/entropy feedback on this original circuit and heat owner.
pub mod thermal;

/// Accepted physical history, held drives and original resistor heat routes.
#[derive(Clone, Copy)]
pub struct PiezoelectricEnthalpyInput<'a> {
    pub model: &'a PiezoelectricCircuit,
    /// Original circuit port order, excluding the appended mechanical port.
    pub electrical_ports: &'a [f64],
    pub force_n: f64,
    /// Exactly one finite thermal receiver per original material resistor.
    pub destinations: &'a [ResistorCell],
}

pub type PiezoelectricEnthalpyReport<S> = CircuitEnthalpyReport<S, PiezoelectricCircuitStep>;

/// A single publishable mechanical/electrical/thermal candidate.
pub struct PiezoelectricEnthalpyCandidate<'m, S> {
    thermal: DepositionSimulation<'m>,
    model: PiezoelectricCircuit,
    report: PiezoelectricEnthalpyReport<S>,
}
impl<'m, S> PiezoelectricEnthalpyCandidate<'m, S> {
    pub const fn thermal(&self) -> &DepositionSimulation<'m> { &self.thermal }
    pub const fn model(&self) -> &PiezoelectricCircuit { &self.model }
    pub const fn report(&self) -> &PiezoelectricEnthalpyReport<S> { &self.report }
    pub fn into_parts(self) -> (DepositionSimulation<'m>, PiezoelectricCircuit,
        PiezoelectricEnthalpyReport<S>) {
        (self.thermal, self.model, self.report)
    }
}

impl ThermalDriver for PiezoelectricCircuitCandidate {
    type Model = PiezoelectricCircuit;
    type Step = PiezoelectricCircuitStep;
    fn account(&self) -> DriverAccount<'_> {
        let step = self.step();
        DriverAccount {
            resistors: &step.resistors, contacts: &[], heat_j: step.heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: step.stored_after_j-step.stored_before_j,
            supplied_j: step.electrical_supply_j+step.mechanical_supply_j+step.parameter_work_j,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) { PiezoelectricCircuitCandidate::into_parts(self) }
}

/// Advance the original piezoelectric rod, circuit and finite resistor bodies
/// through the SAME native coupled residual. Use a prepared material network's
/// resolve method to keep the source claim selection pinned through trials.
///
/// Every material/enthalpy trial restarts from the accepted charge, flux and
/// mechanical history. Joule heat reaches its original cell exactly once.
/// Combined field/strain/kinetic/thermal storage closes against actual electrical
/// work, force-displacement work and authored heat before returning. Refusal
/// from any original source, solver, budget or cancellation returns no endpoint.
///
/// The spatially uniform small-signal rod and original dense circuit limits
/// apply. Only resistor material and thermal-link states evolve here; the
/// piezoelectric source context and its reference mass stay fixed. Thermal
/// rates use the original backward-Euler endpoint convention, so whole-step
/// temporal accuracy is generally first order; energy closure is not an error
/// bound. Electrical contacts and thermal ceramic feedback need their own laws.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'m, S, E>(
    thermal: &DepositionSimulation<'m>, target_s: f64,
    input: PiezoelectricEnthalpyInput<'_>, controls: CircuitEnthalpySolve,
    guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>) -> Result<CircuitThermalCoefficients<S>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<PiezoelectricEnthalpyCandidate<'m, S>, CircuitEnthalpyError<E>> {
    let cancel = RefCell::new(cancelled);
    let original = CircuitEnthalpyInput {
        circuit: input.model.circuit(), ports: input.electrical_ports,
        destinations: input.destinations,
    };
    let (thermal, model, report) = coupled_driver_candidate(
        thermal, target_s, original, &[], controls, guess_j,
        |state, _remaining| {
            let coefficients = resolve(state).map_err(CircuitEnthalpyError::Source)?;
            let circuit = input.model.step_candidate(
                coefficients.resistors, input.electrical_ports, input.force_n,
                target_s-thermal.time_s(), controls.circuit, || (*cancel.borrow_mut())(),
            )?;
            Ok(DrivenCoefficients {
                circuit, sources_w: coefficients.sources_w,
                conductances_w_k: coefficients.conductances_w_k,
                source: coefficients.source, circuit_trials: 1,
            })
        },
        || (*cancel.borrow_mut())(),
    )?;
    Ok(PiezoelectricEnthalpyCandidate { thermal, model, report })
}
