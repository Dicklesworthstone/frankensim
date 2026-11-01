//! Metered gas exchange and finite-body heat retaining original solid inertia.
//!
//! Complete native candidates are rebound to their original source-resolved
//! body. Gas inventories, source models, momentum and accepted port histories
//! are retained without reconstructing snapshots or repeating physical solves.

use super::MaterialDoubleActingGasPiston;
use crate::cavity_phs::piston::double_acting::valve::{
    DoublePistonValveStep, MeteredPistonGas, RectangularMeteringValve,
    thermal::DoubleThermalValveStep,
};
use crate::cavity_phs::piston::double_acting::{DoublePistonInput, DoublePistonSolve};
use crate::cavity_phs::piston::thermal::double_acting::DoubleHeatExchangeSolve;
use crate::cavity_phs::piston::thermal::{ConductanceSample, ThermalPistonError};
use crate::cavity_phs::piston::{ClosedGasParcel, PistonError};
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStep, UniformEnthalpyStepInput};

/// Both metered gas inventories and the original material-bound moving body.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialPistonValveCandidate<G: MeteredPistonGas = ClosedGasParcel> {
    piston: MaterialDoubleActingGasPiston<G>,
    step: DoublePistonValveStep<G>,
}

impl<G: MeteredPistonGas> MaterialPistonValveCandidate<G> {
    /// Complete mechanism retaining its density, geometry and inertia records.
    pub const fn piston(&self) -> &MaterialDoubleActingGasPiston<G> {
        &self.piston
    }

    /// Original nozzle, transfer, shared-motion and conservation accounts.
    pub const fn step(&self) -> &DoublePistonValveStep<G> {
        &self.step
    }

    /// Publish the complete material-bound mechanism and its native report.
    pub fn into_parts(self) -> (MaterialDoubleActingGasPiston<G>, DoublePistonValveStep<G>) {
        (self.piston, self.step)
    }
}

/// Material-bound metered mechanism and its jointly accepted finite heat store.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialThermalValveCandidate<S, G: MeteredPistonGas = ClosedGasParcel> {
    piston: MaterialDoubleActingGasPiston<G>,
    thermal: UniformEnthalpyStep<[S; 2]>,
    step: DoubleThermalValveStep<G>,
}

impl<S, G: MeteredPistonGas> MaterialThermalValveCandidate<S, G> {
    /// Both gases and the original source-bound inertial specimen together.
    pub const fn piston(&self) -> &MaterialDoubleActingGasPiston<G> {
        &self.piston
    }

    /// Native enthalpy endpoint retaining both interface-source evaluations.
    pub const fn thermal(&self) -> &UniformEnthalpyStep<[S; 2]> {
        &self.thermal
    }

    /// Original transfer, motion, heat and enclosing conservation accounts.
    pub const fn step(&self) -> &DoubleThermalValveStep<G> {
        &self.step
    }

    /// Publish the mechanism and thermal endpoint together. Their transported
    /// enthalpy and exchanged face heat have already been applied once.
    pub fn into_parts(
        self,
    ) -> (
        MaterialDoubleActingGasPiston<G>,
        UniformEnthalpyStep<[S; 2]>,
        DoubleThermalValveStep<G>,
    ) {
        (self.piston, self.thermal, self.step)
    }
}

impl<G: MeteredPistonGas> MaterialDoubleActingGasPiston<G> {
    /// Meter the original gas inventories through the moving window and retain
    /// the complete accepted native mechanism on its source-resolved solid.
    ///
    /// Ideal and real-gas owners keep their original nozzle, composition and
    /// energy laws. The native first-order transfer-then-motion split and all
    /// source, geometry, clock and conservation admissions remain unchanged.
    /// Neither a refusal nor final cancellation advances the accepted specimen.
    #[allow(clippy::too_many_arguments)]
    pub fn valve_candidate(
        &self,
        geometry: RectangularMeteringValve,
        input: DoublePistonInput,
        exchange: G::NozzleSolve,
        controls: DoublePistonSolve,
        total_energy_tolerance_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialPistonValveCandidate<G>, PistonError> {
        let native = self.piston.valve_candidate(
            geometry,
            input,
            exchange,
            controls,
            total_energy_tolerance_j,
            &mut cancelled,
        )?;
        let (piston, step) = native.into_parts();
        let candidate = MaterialPistonValveCandidate {
            piston: self.bind_piston_candidate(piston)?,
            step,
        };
        if cancelled() {
            return Err(PistonError::Cancelled);
        }
        Ok(candidate)
    }

    /// Retain the original material inertia while native metering, both gas
    /// pressure works and reciprocal finite-body heat advance as one candidate.
    ///
    /// The callback receives both actual gas midpoints and the original thermal
    /// midpoint; it must be repeatable and must not debit live external state.
    /// The heat store may differ from the moving specimen. Its heating never
    /// silently changes that specimen's density, geometry or canonical mass.
    /// A late refusal or cancellation publishes neither endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn thermal_valve_candidate<S, E>(
        &self,
        geometry: RectangularMeteringValve,
        input: DoublePistonInput,
        body: UniformEnthalpyStepInput<'_>,
        controls: DoublePistonSolve,
        flow: G::NozzleSolve,
        exchange: DoubleHeatExchangeSolve,
        total_energy_tolerance_j: f64,
        resolve: impl FnMut(
            &[G::State; 2],
            EquilibriumPhaseState,
        ) -> Result<[ConductanceSample<S>; 2], E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialThermalValveCandidate<S, G>, ThermalPistonError<E>> {
        let native = self.piston.thermal_valve_candidate(
            geometry,
            input,
            body,
            controls,
            flow,
            exchange,
            total_energy_tolerance_j,
            resolve,
            &mut cancelled,
        )?;
        let (piston, thermal, step) = native.into_parts();
        let candidate = MaterialThermalValveCandidate {
            piston: self.bind_piston_candidate(piston)?,
            thermal,
            step,
        };
        if cancelled() {
            return Err(ThermalPistonError::Cancelled);
        }
        Ok(candidate)
    }
}
