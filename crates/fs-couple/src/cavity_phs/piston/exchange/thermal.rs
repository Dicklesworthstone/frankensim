//! Conservative species flow followed by the existing coupled wall-heat/motion solve.
//!
//! Nozzle transfer is staged once. Every thermal iteration starts from that same
//! transferred chamber, so iteration never withdraws the reservoir a second time.

use super::{ClosedGasParcel, GasNozzle, GasNozzleFlow, GasNozzleStepSolve,
    GasMassTransferStep, GasPiston, NozzleDirection, PistonInput, PistonSolve, poll};
use crate::cavity_phs::piston::thermal::{ConductanceSample, HeatExchangeSolve,
    ThermalPistonError, ThermalPistonStep};
use crate::cavity_phs::piston::thermal::layers::{LayerConductionError, MaterialConductionLayers};
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStep, UniformEnthalpyStepInput};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;

/// Original flow and thermal/mechanical results with a full-system balance.
#[derive(Clone, Debug, PartialEq)]
pub struct NozzleThermalPistonStep {
    /// Initial-pressure-selected upstream side; None for zero mass transfer.
    pub direction: Option<NozzleDirection>,
    /// Source-resolved nozzle state, including actual choking and total enthalpy.
    pub flow: GasNozzleFlow,
    /// Original fixed-volume species/enthalpy exchange, in donor order.
    pub transfer: GasMassTransferStep,
    /// Coupled heat, phase and momentum solve from the post-transfer chamber.
    /// Its local energy residual covers that stage, not the reservoir.
    pub thermal: ThermalPistonStep,
    /// Chamber internal-energy change including nozzle enthalpy and wall heat [J].
    pub chamber_energy_change_j: f64,
    /// Finite-reservoir internal-energy change [J].
    pub reservoir_energy_change_j: f64,
    /// Independently measured total gas mass discrepancy [kg].
    pub total_mass_residual_kg: f64,
    /// dU_chamber+dU_reservoir+dH_body+dK+outside_work-force_work
    /// minus the two prescribed heat sources [J]. Internal transfers cancel.
    pub total_energy_residual_j: f64,
}

/// All three endpoints together; no live gas, body or port history is mutated.
/// `P` retains a material-bound piston when called through `MaterialGasPiston`.
#[derive(Clone, Debug, PartialEq)]
pub struct NozzleThermalPistonCandidate<S, P = GasPiston> {
    piston: P,
    reservoir: ClosedGasParcel,
    thermal: UniformEnthalpyStep<S>,
    step: NozzleThermalPistonStep,
}
impl<S, P> NozzleThermalPistonCandidate<S, P> {
    /// Joint gas, volume and momentum endpoint with the original piston binding.
    pub const fn piston(&self) -> &P { &self.piston }
    /// Tentative finite-reservoir endpoint, including accepted species sources.
    pub const fn reservoir(&self) -> &ClosedGasParcel { &self.reservoir }
    /// Finite-body endpoint with its original enthalpy/phase and interface source.
    pub const fn thermal(&self) -> &UniformEnthalpyStep<S> { &self.thermal }
    /// Source-stage accounts and independent full-system energy/mass balance.
    pub const fn step(&self) -> &NozzleThermalPistonStep { &self.step }
    /// Publish piston, reservoir and thermal state together in the caller's owner.
    pub fn into_parts(self) -> (P, ClosedGasParcel, UniformEnthalpyStep<S>, NozzleThermalPistonStep) {
        (self.piston, self.reservoir, self.thermal, self.step)
    }
}
impl<S> NozzleThermalPistonCandidate<S> {
    pub(in crate::cavity_phs::piston) fn map_piston<P>(self, bind: impl FnOnce(GasPiston) -> P)
        -> NozzleThermalPistonCandidate<S, P>
    {
        NozzleThermalPistonCandidate { piston: bind(self.piston), reservoir: self.reservoir,
            thermal: self.thermal, step: self.step }
    }
}

impl GasPiston {
    /// Exchange species/enthalpy through a nozzle and solve reciprocal wall heat,
    /// gas energy and momentum, returning one complete physical candidate.
    ///
    /// The first-order split transfers once at the accepted initial pressures,
    /// then calls `thermal_candidate` with the resulting gas composition. The
    /// conductance callback sees that gas's actual energy midpoint and the body's
    /// actual enthalpy midpoint, including phase. Both stages span the SAME dt.
    /// Material properties are reevaluated at each thermal/motion trial's energy
    /// state. The nozzle flow itself is split, not monolithically
    /// iterated with motion; temporal refinement remains necessary.
    ///
    /// `input.heat_j` and `body.internal_heat_j` are independent external sources.
    /// Never include nozzle enthalpy or reciprocal wall heat in either input.
    /// `exchange.energy_tolerance_j` bounds the full reservoir/chamber/body/
    /// mechanical balance as well as the original thermal stage. The nozzle's
    /// donor, per-species and source checks retain their original tolerances.
    ///
    /// The body uses the existing fixed-volume, negligible-pV-work enthalpy law;
    /// it may be a chamber wall separate from the translating rigid specimen.
    /// No jet thrust, deformation, reaction chemistry or phase-front model is
    /// inferred. Callbacks must be repeatable and may not debit live state.
    /// Failures and cancellation return none of the three tentative endpoints.
    /// Nonzero flow retains real wall/heat observations with the existing explicit
    /// missing-stream audit boundary; this report closes the full system.
    #[allow(clippy::too_many_arguments)]
    pub fn nozzle_thermal_candidate<S, E>(
        &self, reservoir: &ClosedGasParcel, nozzle: GasNozzle,
        input: PistonInput, body: UniformEnthalpyStepInput<'_>,
        controls: PistonSolve, flow: GasNozzleStepSolve, exchange: HeatExchangeSolve,
        resolve: impl FnMut(&ResolvedThermochemicalGasState, EquilibriumPhaseState)
            -> Result<ConductanceSample<S>, E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NozzleThermalPistonCandidate<S>, ThermalPistonError<E>> {
        poll(&mut cancelled)?;
        let mut next = self.nozzle_transfer(reservoir, nozzle, input.duration_s, flow, &mut cancelled)?;
        let coupled = next.piston.thermal_candidate(input, body, controls, exchange,
            resolve, &mut cancelled)?;
        let (piston, thermal, thermal_step) = coupled.into_parts();
        next.piston = piston;
        // Q enters the gas stage once and leaves the body once. It remains a
        // real chamber heat-port transfer but cancels in the enclosing balance.
        let gas_input = PistonInput { heat_j: thermal_step.piston.gas.supplied_heat_j, ..input };
        let balance = self.nozzle_balance(reservoir, &next, gas_input, flow, exchange.energy_tolerance_j)?;
        let body_energy = body.mass_kg * (thermal.state.specific_enthalpy_j_kg()
            - body.initial.specific_enthalpy_j_kg());
        let total = (((balance.chamber_energy_j + balance.reservoir_energy_j)
            + body_energy) + balance.kinetic_energy_j)
            + balance.exterior_work_j - balance.force_work_j
            - input.heat_j - body.internal_heat_j;
        if !body_energy.is_finite() || !total.is_finite() || total.abs() > exchange.energy_tolerance_j {
            return Err(ThermalPistonError::Invalid("nozzle, finite-body and piston endpoints do not close total energy"));
        }
        self.record_nozzle_account(&mut next, gas_input, thermal_step.piston.midpoint.state().temperature,
            &balance, exchange.energy_tolerance_j)?;
        let step = NozzleThermalPistonStep { direction: next.direction, flow: next.flow,
            transfer: next.transfer, thermal: thermal_step,
            chamber_energy_change_j: balance.chamber_energy_j,
            reservoir_energy_change_j: balance.reservoir_energy_j,
            total_mass_residual_kg: balance.mass_residual_kg, total_energy_residual_j: total };
        poll(&mut cancelled)?;
        Ok(NozzleThermalPistonCandidate { piston: next.piston, reservoir: next.reservoir, thermal, step })
    }

    /// Use the existing source-bound series gas/body conductivity and enthalpy
    /// charts. Composition and material-card swaps reach both heat and motion.
    #[allow(clippy::too_many_arguments)]
    pub fn nozzle_thermal_layers_candidate(
        &self, reservoir: &ClosedGasParcel, nozzle: GasNozzle,
        input: PistonInput, body: UniformEnthalpyStepInput<'_>,
        controls: PistonSolve, flow: GasNozzleStepSolve, exchange: HeatExchangeSolve,
        layers: MaterialConductionLayers<'_>, cancelled: impl FnMut() -> bool,
    ) -> Result<NozzleThermalPistonCandidate<ResolvedMaterialStatePoint>, ThermalPistonError<LayerConductionError>> {
        self.nozzle_thermal_candidate(reservoir, nozzle, input, body, controls, flow, exchange,
            |gas, thermal| layers.resolve(gas, thermal), cancelled)
    }
}

#[cfg(test)]
mod tests;
