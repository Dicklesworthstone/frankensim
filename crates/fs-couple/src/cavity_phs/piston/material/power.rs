//! Continue sourced wall inertia through native state-dependent thermal power.
//!
//! The original powered thermal solve owns both feedback equations and every
//! physical budget. This adapter carries its exact endpoints back into the
//! material-bound mechanism, retaining density receipts and accepted port history.

use super::{ClosedGasParcel, MaterialGasPiston, MaterialThermalPistonCandidate, PistonGas};
use crate::cavity_phs::piston::thermal::power::{PowerFeedbackStep, PowerSample};
use crate::cavity_phs::piston::thermal::{
    ConductanceSample, HeatExchangeSolve, ThermalPistonError,
};
use crate::cavity_phs::piston::{PistonInput, PistonSolve};
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStepInput};

/// Joint powered endpoint retaining the original translating material specimen.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialPoweredThermalCandidate<I, P, G: PistonGas = ClosedGasParcel> {
    coupled: MaterialThermalPistonCandidate<I, G>,
    power: PowerSample<P>,
    step: PowerFeedbackStep,
}

impl<I, P, G: PistonGas> MaterialPoweredThermalCandidate<I, P, G> {
    /// Original gas/motion and finite-body endpoint with sourced wall inertia.
    pub const fn coupled(&self) -> &MaterialThermalPistonCandidate<I, G> {
        &self.coupled
    }
    /// Original concrete power source evaluated at the accepted energy midpoints.
    pub const fn power(&self) -> &PowerSample<P> {
        &self.power
    }
    /// Native power consistency, combined energy balance and work counts.
    pub const fn step(&self) -> &PowerFeedbackStep {
        &self.step
    }
    /// Publish the mechanism and thermal endpoint together. Source energy and
    /// interface heat are already included and must not be applied again.
    pub fn into_parts(
        self,
    ) -> (
        MaterialThermalPistonCandidate<I, G>,
        PowerSample<P>,
        PowerFeedbackStep,
    ) {
        (self.coupled, self.power, self.step)
    }
}

impl<G: PistonGas> MaterialGasPiston<G> {
    /// Apply the original state-dependent body-power solve without rebuilding
    /// the material piston, replacing its inertia, or resetting its port account.
    /// The thermal store may be distinct from the translating specimen; its
    /// heating or phase change does not silently change mechanical mass.
    ///
    /// Both callbacks receive the actual gas and body energy midpoints. They
    /// must be repeatable and may not debit live external stores. Original
    /// source values, failure types, work limits and balance gates are retained.
    /// No endpoint is published on refusal or final cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn powered_thermal_candidate<I, P, E>(
        &self,
        input: PistonInput,
        thermal_body: UniformEnthalpyStepInput<'_>,
        controls: PistonSolve,
        exchange: HeatExchangeSolve,
        source: HeatExchangeSolve,
        conductance: impl FnMut(&G::State, EquilibriumPhaseState) -> Result<ConductanceSample<I>, E>,
        power: impl FnMut(&G::State, EquilibriumPhaseState) -> Result<PowerSample<P>, E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialPoweredThermalCandidate<I, P, G>, ThermalPistonError<E>> {
        let candidate = self.piston.powered_thermal_candidate(
            input,
            thermal_body,
            controls,
            exchange,
            source,
            conductance,
            power,
            &mut cancelled,
        )?;
        let (coupled, power, step) = candidate.into_parts();
        let (piston, thermal, thermal_step) = coupled.into_parts();
        let candidate = MaterialPoweredThermalCandidate {
            coupled: MaterialThermalPistonCandidate {
                piston: Self {
                    body: self.body.clone(),
                    piston,
                },
                thermal,
                step: thermal_step,
            },
            power,
            step,
        };
        if cancelled() {
            return Err(ThermalPistonError::Cancelled);
        }
        Ok(candidate)
    }
}

#[cfg(test)]
mod tests;
