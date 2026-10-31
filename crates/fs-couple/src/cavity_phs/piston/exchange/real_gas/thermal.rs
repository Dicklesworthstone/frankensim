//! Real-gas mass/enthalpy transfer followed by reciprocal body heat and motion.
//!
//! Transfer is staged ONCE. Every thermal trial starts from that same proposed
//! chamber, retaining real density, calorics and sources without re-debiting the
//! finite reservoir. Stream enthalpy and wall heat each have one physical owner.

use super::{
    GasNozzle, GasPiston, IsochoricRealGasParcel, NozzleDirection, PistonInput, PistonSolve,
    RealGasMassTransferStep, RealGasNozzleFlow, RealGasNozzleStepSolve,
    ResolvedThermochemicalRealGasState, poll,
};
use crate::cavity_phs::piston::thermal::{
    ConductanceSample, HeatExchangeSolve, ThermalPistonError, ThermalPistonStep,
};
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStep, UniformEnthalpyStepInput};

/// Original real-gas flow and thermal/motion accounts enclosing all three stores.
#[derive(Clone, Debug, PartialEq)]
pub struct RealGasNozzleThermalPistonStep {
    /// Initial-pressure-selected upstream side; None for a shut/equilibrated path.
    pub direction: Option<NozzleDirection>,
    /// Original real-gas throat, choking and flow-rate state.
    pub flow: RealGasNozzleFlow,
    /// Original conservative mass/enthalpy transfer in donor order.
    pub transfer: RealGasMassTransferStep,
    /// Reciprocal finite-body heat and momentum from the transferred chamber.
    pub thermal: ThermalPistonStep<ResolvedThermochemicalRealGasState>,
    /// Chamber internal-energy change including stream and wall exchange [J].
    pub chamber_energy_change_j: f64,
    /// Finite-reservoir internal-energy change [J].
    pub reservoir_energy_change_j: f64,
    /// Actual net mass defect across the two accepted gas endpoints [kg].
    pub total_mass_residual_kg: f64,
    /// Original transfer defect plus the complete gas/body/mechanical defect [J].
    /// Stream enthalpy and reciprocal heat cancel inside this enclosing account.
    pub total_energy_residual_j: f64,
}

/// All three endpoints together, retaining a material-bound piston through `P`.
#[derive(Clone, Debug, PartialEq)]
pub struct RealGasNozzleThermalPistonCandidate<S, P = GasPiston<IsochoricRealGasParcel>> {
    piston: P,
    reservoir: IsochoricRealGasParcel,
    thermal: UniformEnthalpyStep<S>,
    step: RealGasNozzleThermalPistonStep,
}
impl<S, P> RealGasNozzleThermalPistonCandidate<S, P> {
    /// Complete tentative gas, momentum, source and accepted port history.
    pub const fn piston(&self) -> &P {
        &self.piston
    }
    /// Tentative finite upstream/downstream reservoir.
    pub const fn reservoir(&self) -> &IsochoricRealGasParcel {
        &self.reservoir
    }
    /// Finite-body endpoint with its actual interface sample.
    pub const fn thermal(&self) -> &UniformEnthalpyStep<S> {
        &self.thermal
    }
    /// Original stage reports and enclosing mass/energy accounting.
    pub const fn step(&self) -> &RealGasNozzleThermalPistonStep {
        &self.step
    }
    /// Publish the gas/mechanism, reservoir and body together in the caller.
    pub fn into_parts(
        self,
    ) -> (
        P,
        IsochoricRealGasParcel,
        UniformEnthalpyStep<S>,
        RealGasNozzleThermalPistonStep,
    ) {
        (self.piston, self.reservoir, self.thermal, self.step)
    }
}
impl<S> RealGasNozzleThermalPistonCandidate<S> {
    pub(in crate::cavity_phs::piston) fn map_piston<P>(
        self,
        bind: impl FnOnce(GasPiston<IsochoricRealGasParcel>) -> P,
    ) -> RealGasNozzleThermalPistonCandidate<S, P> {
        RealGasNozzleThermalPistonCandidate {
            piston: bind(self.piston),
            reservoir: self.reservoir,
            thermal: self.thermal,
            step: self.step,
        }
    }
}

impl GasPiston<IsochoricRealGasParcel> {
    /// Transfer real gas through a nozzle and jointly solve finite-body heat,
    /// source gas energy and piston momentum without mutating any live input.
    ///
    /// This retains the original first-order transfer-then-heat/motion split.
    /// Flow uses accepted initial pressures; subsequent thermal trials share
    /// the one transferred chamber and resolve its actual energy midpoint.
    /// Both stages span `input.duration_s`; temporal refinement is required.
    /// The original PR flow owner admits identical frozen composition and
    /// material sources, a homogeneous stable branch, and its donor/sonic bounds.
    ///
    /// `input.heat_j` and `body.internal_heat_j` are external sources only.
    /// Neither includes transported enthalpy or reciprocal wall heat. The gas
    /// stage receives Q once and the body loses it once. Original-law increment
    /// ledgers close the entire system without subtracting absolute formation
    /// energies. `exchange.energy_tolerance_j` also bounds this enclosing account.
    ///
    /// The same body enthalpy law, negligible body pV work, acoustic clock and
    /// repeatable callback requirements remain. No jet thrust, phase-equilibrium
    /// closure or monolithic flow/heat iteration is inferred. Any refusal or
    /// cancellation returns none of the three tentative endpoints.
    #[allow(clippy::too_many_arguments)]
    pub fn real_gas_nozzle_thermal_candidate<S, E>(
        &self,
        reservoir: &IsochoricRealGasParcel,
        nozzle: GasNozzle,
        input: PistonInput,
        body: UniformEnthalpyStepInput<'_>,
        controls: PistonSolve,
        flow: RealGasNozzleStepSolve,
        exchange: HeatExchangeSolve,
        resolve: impl FnMut(
            &ResolvedThermochemicalRealGasState,
            EquilibriumPhaseState,
        ) -> Result<ConductanceSample<S>, E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RealGasNozzleThermalPistonCandidate<S>, ThermalPistonError<E>> {
        poll(&mut cancelled)?;
        let mut next =
            self.real_nozzle_transfer(reservoir, nozzle, input.duration_s, flow, &mut cancelled)?;
        let coupled = next.piston.thermal_candidate(
            input,
            body,
            controls,
            exchange,
            resolve,
            &mut cancelled,
        )?;
        let (piston, thermal, thermal_step) = coupled.into_parts();
        next.piston = piston;
        let balance = self.real_nozzle_balance(
            reservoir,
            &next,
            &thermal_step.piston,
            flow,
            exchange.energy_tolerance_j,
        )?;
        let total_energy_residual_j =
            next.transfer.energy_residual_j + thermal_step.total_energy_residual_j;
        if !total_energy_residual_j.is_finite()
            || total_energy_residual_j.abs() > exchange.energy_tolerance_j
        {
            return Err(ThermalPistonError::Invalid(
                "real-gas nozzle, finite body and piston do not close total energy",
            ));
        }
        // A local gas heat-port observation includes Q, although Q cancels
        // from the enclosing body/gas account. Stream-port absence stays explicit.
        let gas_input = PistonInput {
            heat_j: thermal_step.piston.gas.supplied_heat_j,
            ..input
        };
        self.record_real_nozzle_account(
            &mut next,
            gas_input,
            &thermal_step.piston,
            &balance,
            exchange.energy_tolerance_j,
        );
        let step = RealGasNozzleThermalPistonStep {
            direction: next.direction,
            flow: next.flow,
            transfer: next.transfer,
            thermal: thermal_step,
            chamber_energy_change_j: balance.chamber_energy_change_j,
            reservoir_energy_change_j: balance.reservoir_energy_change_j,
            total_mass_residual_kg: balance.total_mass_residual_kg,
            total_energy_residual_j,
        };
        poll(&mut cancelled)?;
        Ok(RealGasNozzleThermalPistonCandidate {
            piston: next.piston,
            reservoir: next.reservoir,
            thermal,
            step,
        })
    }
}
