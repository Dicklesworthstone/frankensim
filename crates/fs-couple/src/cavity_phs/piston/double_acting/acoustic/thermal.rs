//! One finite enthalpy body and both closed gases under reciprocal acoustic load.
//!
//! This uses the original two-face heat solver directly. A caller does not need
//! to invent a closed valve or nozzle just to turn physical heating into loaded
//! wall motion and acoustic output.

use super::{
    AcousticMotionError, AcousticMotionFace, AcousticMotionSolve, AcousticMotionStep, finish,
    prepare,
};
use crate::cavity_phs::piston::double_acting::{
    DoubleActingGasPiston, DoublePistonInput, DoublePistonSolve,
};
use crate::cavity_phs::piston::thermal::double_acting::{
    DoubleHeatExchangeSolve, DoubleThermalPistonCandidate,
};
use crate::cavity_phs::piston::thermal::{ConductanceSample, ThermalPistonError};
use crate::dg_port::DgPortLoad;
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStepInput};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;

/// Original finite-body/two-gas candidate and all continuing acoustic endpoints.
#[derive(Debug)]
pub struct AcousticThermalPistonCandidate<S> {
    native: DoubleThermalPistonCandidate<S>,
    loads: Vec<DgPortLoad>,
    motion: AcousticMotionStep,
    total_energy_residual_j: f64,
}
impl<S> AcousticThermalPistonCandidate<S> {
    pub const fn native(&self) -> &DoubleThermalPistonCandidate<S> {
        &self.native
    }
    pub fn loads(&self) -> &[DgPortLoad] {
        &self.loads
    }
    pub const fn motion(&self) -> &AcousticMotionStep {
        &self.motion
    }
    /// Gas + single wall + single thermal body + original acoustic balance.
    pub const fn total_energy_residual_j(&self) -> f64 {
        self.total_energy_residual_j
    }
    /// Publish every returned state together; face heat is already exchanged.
    pub fn into_parts(
        self,
    ) -> (
        DoubleThermalPistonCandidate<S>,
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

impl DoubleActingGasPiston {
    /// Resolve both gas/body heat laws on the actual acoustically loaded motion.
    ///
    /// All original heat trials retain the same prepared affine acoustic force
    /// in the native momentum equation. The original bounded heat solver owns
    /// phase/enthalpy inversion, conductance evaluation and source admission;
    /// each complete acoustic owner advances once after convergence. No valve,
    /// mass exchange, repeated acoustic marching or separate thermal integrator
    /// is introduced. The callback must not debit live state.
    ///
    /// Source/refinement controls retain their original meaning. This remains
    /// fixed-body enthalpy coupled to small-displacement mean-port acoustics;
    /// sensible/latent heating does not silently change the radiator geometry.
    #[allow(clippy::too_many_arguments)]
    pub fn acoustic_thermal_candidate<S, E>(
        &self,
        faces: &[AcousticMotionFace<'_>],
        input: DoublePistonInput,
        body: UniformEnthalpyStepInput<'_>,
        controls: DoublePistonSolve,
        exchange: DoubleHeatExchangeSolve,
        acoustic: AcousticMotionSolve,
        resolve: impl FnMut(
            &[ResolvedThermochemicalGasState; 2],
            EquilibriumPhaseState,
        ) -> Result<[ConductanceSample<S>; 2], E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AcousticThermalPistonCandidate<S>, AcousticMotionError<ThermalPistonError<E>>> {
        let prepared = prepare(
            faces,
            input.duration_s,
            self.ports.next_tick() - 1,
            acoustic,
            &mut cancelled,
        )?;
        let force = if prepared.history_force_n == 0.0 {
            input.external_force_n
        } else {
            input.external_force_n + prepared.history_force_n
        };
        if !force.is_finite() {
            return Err(AcousticMotionError::Invalid(
                "acoustic history and prescribed force overflowed",
            ));
        }
        let native = match self.thermal_candidate_with_midpoint_resistance(
            DoublePistonInput {
                external_force_n: force,
                ..input
            },
            body,
            controls,
            exchange,
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
            &native.step().piston,
            input.external_force_n,
            &mut cancelled,
        )?;
        let total = (native.step().total_energy_residual_j + motion.acoustic_energy_residual_j)
            + motion.transfer_residual_j;
        if !total.is_finite() || total.abs() > acoustic.energy_tolerance_j {
            return Err(AcousticMotionError::Invalid(
                "both gases, one thermal body, common wall and acoustic stores do not close external work",
            ));
        }
        if cancelled() {
            return Err(AcousticMotionError::Cancelled);
        }
        Ok(AcousticThermalPistonCandidate {
            native,
            loads,
            motion,
            total_energy_residual_j: total,
        })
    }

    /// Resolve both existing source-bound material conduction layers at the
    /// same gas/body midpoints used by the loaded thermal solve.
    #[allow(clippy::too_many_arguments)]
    pub fn acoustic_thermal_layers_candidate(
        &self,
        faces: &[AcousticMotionFace<'_>],
        input: DoublePistonInput,
        body: UniformEnthalpyStepInput<'_>,
        controls: DoublePistonSolve,
        exchange: DoubleHeatExchangeSolve,
        acoustic: AcousticMotionSolve,
        layers: [crate::cavity_phs::piston::thermal::layers::MaterialConductionLayers<'_>; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        AcousticThermalPistonCandidate<fs_material::state_point::ResolvedMaterialStatePoint>,
        AcousticMotionError<
            ThermalPistonError<crate::cavity_phs::piston::thermal::layers::LayerConductionError>,
        >,
    > {
        self.acoustic_thermal_candidate(
            faces,
            input,
            body,
            controls,
            exchange,
            acoustic,
            |gases, thermal| {
                Ok([
                    layers[0].resolve(&gases[0], thermal)?,
                    layers[1].resolve(&gases[1], thermal)?,
                ])
            },
            cancelled,
        )
    }
}
