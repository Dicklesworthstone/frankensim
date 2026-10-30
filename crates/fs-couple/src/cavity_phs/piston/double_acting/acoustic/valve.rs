//! Reciprocal acoustic loads on the original metered two-gas motion.
//!
//! The valve transports species and total enthalpy between its thermodynamic
//! chambers. Separate fixed-geometry acoustic domains load additional faces of
//! the same body. Their affine mean-pressure reaction enters the original
//! momentum solve, including every finite-body heat trial. Each acoustic load
//! advances once, after all gas, wall and optional thermal states converge.

use super::{
    AcousticMotionError, AcousticMotionFace, AcousticMotionSolve, AcousticMotionStep, finish,
    prepare,
};
use crate::cavity_phs::piston::PistonError;
use crate::cavity_phs::piston::double_acting::valve::thermal::DoubleThermalValveCandidate;
use crate::cavity_phs::piston::double_acting::valve::{
    DoublePistonValveCandidate, RectangularMeteringValve,
};
use crate::cavity_phs::piston::double_acting::{
    DoubleActingGasPiston, DoublePistonInput, DoublePistonSolve,
};
use crate::cavity_phs::piston::thermal::double_acting::DoubleHeatExchangeSolve;
use crate::cavity_phs::piston::thermal::{ConductanceSample, ThermalPistonError};
use crate::dg_port::DgPortLoad;
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStepInput};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::GasNozzleStepSolve;

/// The original finite-gas valve candidate and every continuing acoustic load.
#[derive(Debug)]
pub struct AcousticValveCandidate {
    native: DoublePistonValveCandidate,
    loads: Vec<DgPortLoad>,
    motion: AcousticMotionStep,
    total_energy_residual_j: f64,
}

impl AcousticValveCandidate {
    /// Native gas inventories, one wall and their once-only transport report.
    pub const fn native(&self) -> &DoublePistonValveCandidate {
        &self.native
    }
    /// Accepted-candidate loads in the original face order.
    pub fn loads(&self) -> &[DgPortLoad] {
        &self.loads
    }
    /// Actual displacement flows, pressure work, losses and exterior exchange.
    pub const fn motion(&self) -> &AcousticMotionStep {
        &self.motion
    }
    /// Enclosing gas, kinetic and acoustic balance against external supplies.
    pub const fn total_energy_residual_j(&self) -> f64 {
        self.total_energy_residual_j
    }
    /// Publish these participants together; no stream enthalpy is added again.
    pub fn into_parts(
        self,
    ) -> (
        DoublePistonValveCandidate,
        Vec<DgPortLoad>,
        AcousticMotionStep,
        f64,
    ) {
        (
            self.native,
            self.loads,
            self.motion,
            self.total_energy_residual_j,
        )
    }
}

/// The same gas/valve/wall/acoustic candidate plus one original finite body.
#[derive(Debug)]
pub struct AcousticThermalValveCandidate<S> {
    native: DoubleThermalValveCandidate<S>,
    loads: Vec<DgPortLoad>,
    motion: AcousticMotionStep,
    total_energy_residual_j: f64,
}

impl<S> AcousticThermalValveCandidate<S> {
    /// Original gas, wall, finite enthalpy endpoint and both source evaluations.
    pub const fn native(&self) -> &DoubleThermalValveCandidate<S> {
        &self.native
    }
    pub fn loads(&self) -> &[DgPortLoad] {
        &self.loads
    }
    pub const fn motion(&self) -> &AcousticMotionStep {
        &self.motion
    }
    /// Includes the ONE body's actual enthalpy change; latent heat is not extra.
    pub const fn total_energy_residual_j(&self) -> f64 {
        self.total_energy_residual_j
    }
    /// Consume the complete proposal without separately replaying any transfer.
    pub fn into_parts(
        self,
    ) -> (
        DoubleThermalValveCandidate<S>,
        Vec<DgPortLoad>,
        AcousticMotionStep,
        f64,
    ) {
        (
            self.native,
            self.loads,
            self.motion,
            self.total_energy_residual_j,
        )
    }
}

fn loaded_input<E>(
    input: DoublePistonInput,
    history_force_n: f64,
) -> Result<DoublePistonInput, AcousticMotionError<E>> {
    // Preserve the exact original force, including signed zero, when decoupled.
    let force = if history_force_n == 0.0 {
        input.external_force_n
    } else {
        input.external_force_n + history_force_n
    };
    if !force.is_finite() {
        return Err(AcousticMotionError::Invalid(
            "acoustic history and prescribed force overflowed",
        ));
    }
    Ok(DoublePistonInput {
        external_force_n: force,
        ..input
    })
}

fn enclosing_energy<E>(
    native_residual_j: f64,
    motion: &AcousticMotionStep,
    tolerance_j: f64,
) -> Result<f64, AcousticMotionError<E>> {
    // The native residual subtracts ACTUAL total wall work. Enclosing the DG
    // owners cancels its acoustic part through their measured reciprocal work.
    let total =
        (native_residual_j + motion.acoustic_energy_residual_j) + motion.transfer_residual_j;
    if !total.is_finite() || total.abs() > tolerance_j {
        return Err(AcousticMotionError::Invalid(
            "metered gas, wall, thermal and acoustic stores do not close external work",
        ));
    }
    Ok(total)
}

impl DoubleActingGasPiston {
    /// Meter the two original gas inventories while independent acoustic loads
    /// react on additional faces of their common moving body.
    ///
    /// Transport is sampled once from the accepted valve opening and gas pair.
    /// The original pressure/momentum solve then sees every acoustic history
    /// force and the summed nonnegative mean velocity resistance. It is not a
    /// delayed acoustic-force sample or a second mechanical trajectory.
    ///
    /// The acoustic controls' energy tolerance also gates the original valve's
    /// complete pair balance. Its native donor budget, final pressure-order
    /// gate and first-order frozen-opening approximation remain in force.
    /// Acoustic inlet flows are signed face area times wall velocity; they are
    /// NOT the valve's mass rate divided by either chamber density. No chamber
    /// pressure storage or stream enthalpy is counted as DG energy.
    #[allow(clippy::too_many_arguments)]
    pub fn acoustic_valve_candidate(
        &self,
        faces: &[AcousticMotionFace<'_>],
        geometry: RectangularMeteringValve,
        input: DoublePistonInput,
        flow: GasNozzleStepSolve,
        controls: DoublePistonSolve,
        acoustic: AcousticMotionSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AcousticValveCandidate, AcousticMotionError<PistonError>> {
        let prepared = prepare(
            faces,
            input.duration_s,
            self.ports.next_tick() - 1,
            acoustic,
            &mut cancelled,
        )?;
        let applied = loaded_input(input, prepared.history_force_n)?;
        let native = match self.valve_candidate_with_midpoint_resistance(
            geometry,
            applied,
            flow,
            controls,
            acoustic.energy_tolerance_j,
            prepared.resistance_n_s_m,
            &mut cancelled,
        ) {
            Ok(candidate) => candidate,
            Err(PistonError::Cancelled) => return Err(AcousticMotionError::Cancelled),
            Err(error) => return Err(AcousticMotionError::Mechanical(error)),
        };
        let (loads, motion) = finish(
            prepared,
            &native.step().piston,
            input.external_force_n,
            &mut cancelled,
        )?;
        let total = enclosing_energy(
            native.step().total_energy_residual_j,
            &motion,
            acoustic.energy_tolerance_j,
        )?;
        if cancelled() {
            return Err(AcousticMotionError::Cancelled);
        }
        Ok(AcousticValveCandidate {
            native,
            loads,
            motion,
            total_energy_residual_j: total,
        })
    }

    /// Close native metering, BOTH gas/body heat laws, common wall momentum and
    /// every independent acoustic load in one immutable candidate.
    ///
    /// Every heat trial starts from the same once-transferred gas pair and
    /// original finite body. Its momentum solve includes the same prepared
    /// acoustic reaction; the native body/gas midpoint law therefore sees the
    /// acoustically loaded motion. No trial advances a DG/rational/channel
    /// history or repeats mass transport. Complete acoustic candidates are
    /// prepared only after the bounded original thermal solve has converged.
    ///
    /// The callback must be repeatable and cannot debit live state. Native
    /// source refusals, exhausted heat/impulse budgets, pressure reversal,
    /// acoustic refusal or final cancellation publish no participant.
    ///
    /// The body's original fixed-volume enthalpy model and valve's frozen-
    /// opening split remain. This does not add jet recoil, thermal shape change,
    /// moving-mesh ducts or heterogeneous through-flow acoustic junctions.
    #[allow(clippy::too_many_arguments)]
    pub fn acoustic_thermal_valve_candidate<S, E>(
        &self,
        faces: &[AcousticMotionFace<'_>],
        geometry: RectangularMeteringValve,
        input: DoublePistonInput,
        body: UniformEnthalpyStepInput<'_>,
        controls: DoublePistonSolve,
        flow: GasNozzleStepSolve,
        exchange: DoubleHeatExchangeSolve,
        acoustic: AcousticMotionSolve,
        resolve: impl FnMut(
            &[ResolvedThermochemicalGasState; 2],
            EquilibriumPhaseState,
        ) -> Result<[ConductanceSample<S>; 2], E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AcousticThermalValveCandidate<S>, AcousticMotionError<ThermalPistonError<E>>> {
        let prepared = prepare(
            faces,
            input.duration_s,
            self.ports.next_tick() - 1,
            acoustic,
            &mut cancelled,
        )?;
        let applied = loaded_input(input, prepared.history_force_n)?;
        let native = match self.thermal_valve_candidate_with_midpoint_resistance(
            geometry,
            applied,
            body,
            controls,
            flow,
            exchange,
            acoustic.energy_tolerance_j,
            prepared.resistance_n_s_m,
            resolve,
            &mut cancelled,
        ) {
            Ok(candidate) => candidate,
            Err(ThermalPistonError::Cancelled) => return Err(AcousticMotionError::Cancelled),
            Err(error) => return Err(AcousticMotionError::Mechanical(error)),
        };
        let (loads, motion) = finish(
            prepared,
            &native.step().valve.piston,
            input.external_force_n,
            &mut cancelled,
        )?;
        let total = enclosing_energy(
            native.step().total_energy_residual_j,
            &motion,
            acoustic.energy_tolerance_j,
        )?;
        if cancelled() {
            return Err(AcousticMotionError::Cancelled);
        }
        Ok(AcousticThermalValveCandidate {
            native,
            loads,
            motion,
            total_energy_residual_j: total,
        })
    }
}

#[cfg(test)]
mod tests;
