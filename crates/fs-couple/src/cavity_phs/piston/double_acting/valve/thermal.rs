//! The original metered gas pair and one shared finite thermal body.
//!
//! Sample transport once at the accepted opening, then solve both face heat
//! flows with the existing enthalpy/momentum owner. Every heat trial sees the
//! same post-transfer gases; it cannot move more mass or debit the body twice.

use super::{DoublePistonValveStep, RectangularMeteringValve};
use crate::cavity_phs::piston::double_acting::{
    DoubleActingGasPiston, DoublePistonInput, DoublePistonSolve,
};
use crate::cavity_phs::piston::thermal::double_acting::{
    DoubleHeatExchangeSolve, DoubleThermalPistonStep,
};
use crate::cavity_phs::piston::thermal::{ConductanceSample, ThermalPistonError};
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStep, UniformEnthalpyStepInput};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::GasNozzleStepSolve;

/// One accepted transport, thermal and shared-wall proposal. Native gas and
/// mechanical accounts appear once, inside `valve.piston`.
#[derive(Clone, Debug, PartialEq)]
pub struct DoubleThermalValveStep {
    /// Original valve/motion report; its gas heats include the two accepted
    /// body exchanges. Those exchanges are internal to this enclosing report.
    pub valve: DoublePistonValveStep,
    /// Body phase at midpoint enthalpy, evaluated by its original source chart.
    pub thermal_midpoint: EquilibriumPhaseState,
    /// Actual source conductance at each accepted gas/body midpoint [W/K].
    pub conductance_w_k: [f64; 2],
    /// Heat OUT of the one finite body and INTO [primary, opposing] gas [J].
    pub exchanged_heat_j: [f64; 2],
    /// Each accepted Q minus dt*G*(T_body_mid-T_gas_mid) [J].
    pub source_residual_j: [f64; 2],
    /// Actual body mass times its specific enthalpy change [J].
    pub thermal_energy_change_j: f64,
    /// Body storage + both face heats - prescribed body heat [J].
    pub thermal_energy_residual_j: f64,
    /// Both gases + one wall kinetic + one body enthalpy change, minus actual
    /// external force work and the three prescribed heats [J]. Transported
    /// enthalpy and the two face heats cancel internally. Native conserved
    /// species datums cancel before arithmetic, as in the original valve owner.
    pub total_energy_residual_j: f64,
    /// Original bounded heat/momentum evaluations, excluding the single native
    /// transport stage, whose bounded work is reported by `valve.flow`.
    pub joint_trials: usize,
    pub material_trials: usize,
}

/// Return the gas pair, common wall and finite thermal endpoint together.
/// Neither accepted piston nor external body is mutated by construction.
#[derive(Clone, Debug, PartialEq)]
pub struct DoubleThermalValveCandidate<S> {
    piston: DoubleActingGasPiston,
    thermal: UniformEnthalpyStep<[S; 2]>,
    step: DoubleThermalValveStep,
}

impl<S> DoubleThermalValveCandidate<S> {
    pub const fn piston(&self) -> &DoubleActingGasPiston {
        &self.piston
    }

    /// Native enthalpy endpoint and both original interface-source evaluations.
    /// Its external heat is already exchanged with the returned gases.
    pub const fn thermal(&self) -> &UniformEnthalpyStep<[S; 2]> {
        &self.thermal
    }

    pub const fn step(&self) -> &DoubleThermalValveStep {
        &self.step
    }

    pub fn into_parts(
        self,
    ) -> (
        DoubleActingGasPiston,
        UniformEnthalpyStep<[S; 2]>,
        DoubleThermalValveStep,
    ) {
        (self.piston, self.thermal, self.step)
    }
}

impl DoubleActingGasPiston {
    /// Exchange native species/enthalpy through the moving metering edge, then
    /// jointly solve the existing two-face finite-body heat and wall momentum.
    ///
    /// Different gas compositions and caloric laws keep their original owners.
    /// Flow is sampled ONCE from the accepted position and reservoirs; nested
    /// heat trials borrow the same post-transfer pair and original body. The
    /// returned valve report is finalized with actual accepted gas heats, then
    /// enclosed by the single body's independently measured enthalpy change.
    /// Pressure ordering is rechecked after both heat exchange and motion.
    ///
    /// `resolve` has the original repeatable thermal-candidate contract: it
    /// returns conductances and source evaluations at the supplied midpoints,
    /// and must not mutate live external state. Native heat brackets and joint
    /// trial limits bound its calls. Source failure, exhausted work or final
    /// cancellation returns no partially advanced gas, wall or body endpoint.
    /// Commit all three returned participants together without depositing the
    /// native transported enthalpy or accepted face heat a second time.
    ///
    /// The original fixed-body enthalpy assumptions and first-order frozen-
    /// opening split remain. No jet recoil, contact, thermal deformation of the
    /// metering geometry or new heterogeneous acoustic junction is introduced.
    #[allow(clippy::too_many_arguments)]
    pub fn thermal_valve_candidate<S, E>(
        &self,
        geometry: RectangularMeteringValve,
        input: DoublePistonInput,
        body: UniformEnthalpyStepInput<'_>,
        controls: DoublePistonSolve,
        flow: GasNozzleStepSolve,
        exchange: DoubleHeatExchangeSolve,
        total_energy_tolerance_j: f64,
        resolve: impl FnMut(
            &[ResolvedThermochemicalGasState; 2],
            EquilibriumPhaseState,
        ) -> Result<[ConductanceSample<S>; 2], E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DoubleThermalValveCandidate<S>, ThermalPistonError<E>> {
        self.thermal_valve_candidate_with_midpoint_resistance(
            geometry, input, body, controls, flow, exchange,
            total_energy_tolerance_j, 0.0, resolve, cancelled,
        )
    }

    // The original post-transfer heat solver includes the affine acoustic
    // reaction without replaying transport or marching an acoustic history.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::cavity_phs::piston) fn thermal_valve_candidate_with_midpoint_resistance<S, E>(
        &self,
        geometry: RectangularMeteringValve,
        input: DoublePistonInput,
        body: UniformEnthalpyStepInput<'_>,
        controls: DoublePistonSolve,
        flow: GasNozzleStepSolve,
        exchange: DoubleHeatExchangeSolve,
        total_energy_tolerance_j: f64,
        resistance_n_s_m: f64,
        resolve: impl FnMut(
            &[ResolvedThermochemicalGasState; 2],
            EquilibriumPhaseState,
        ) -> Result<[ConductanceSample<S>; 2], E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DoubleThermalValveCandidate<S>, ThermalPistonError<E>> {
        if cancelled() {
            return Err(ThermalPistonError::Cancelled);
        }
        if !total_energy_tolerance_j.is_finite()
            || total_energy_tolerance_j <= 0.0
            || !resistance_n_s_m.is_finite()
            || resistance_n_s_m < 0.0
            || body.duration_s != input.duration_s
        {
            return Err(ThermalPistonError::Invalid(
                "metered thermal coupling requires one physical clock and a positive finite total energy tolerance",
            ));
        }
        let (transported, transfer) =
            self.prepare_valve(geometry, input.duration_s, flow, &mut cancelled)?;
        let candidate = transported.thermal_candidate_with_midpoint_resistance(
            input,
            body,
            controls,
            exchange,
            resistance_n_s_m,
            resolve,
            &mut cancelled,
        )?;
        let (next, thermal, native) = candidate.into_parts();
        let DoubleThermalPistonStep {
            piston,
            thermal_midpoint,
            conductance_w_k,
            exchanged_heat_j,
            source_residual_j,
            joint_trials,
            material_trials,
            ..
        } = native;
        let actual_input = DoublePistonInput {
            heat_j: std::array::from_fn(|i| input.heat_j[i] + exchanged_heat_j[i]),
            external_force_n: super::super::midpoint_force_n(
                input, resistance_n_s_m, piston.displacement_m,
            )?,
            ..input
        };
        let accepted = self.finish_valve(
            geometry,
            actual_input,
            flow,
            total_energy_tolerance_j,
            transfer,
            next,
            piston,
            &mut cancelled,
        )?;
        let (piston, valve) = accepted.into_parts();
        let thermal_change = body.mass_kg
            * (thermal.state.specific_enthalpy_j_kg() - body.initial.specific_enthalpy_j_kg());
        let net_heat = exchanged_heat_j[0] + exchanged_heat_j[1];
        let thermal_residual = (thermal_change + net_heat) - body.internal_heat_j;
        // `valve` re-evaluated actual endpoint caloric and kinetic changes.
        // Its gas-only boundary subtracts both body heats; enclosing the body
        // cancels precisely those heats while retaining native datum handling.
        let total = valve.total_energy_residual_j + thermal_residual;
        if ![
            thermal_change,
            net_heat,
            thermal_residual,
            total,
            thermal.external_heat_j,
        ]
        .iter()
        .all(|value| value.is_finite())
            || (thermal.external_heat_j + net_heat).abs() > total_energy_tolerance_j
            || thermal_residual.abs() > total_energy_tolerance_j
            || total.abs() > total_energy_tolerance_j
        {
            return Err(ThermalPistonError::Invalid(
                "metered gas, shared wall and finite thermal body do not close their common energy account",
            ));
        }
        if cancelled() {
            return Err(ThermalPistonError::Cancelled);
        }
        Ok(DoubleThermalValveCandidate {
            piston,
            thermal,
            step: DoubleThermalValveStep {
                valve,
                thermal_midpoint,
                conductance_w_k,
                exchanged_heat_j,
                source_residual_j,
                thermal_energy_change_j: thermal_change,
                thermal_energy_residual_j: thermal_residual,
                total_energy_residual_j: total,
                joint_trials,
                material_trials,
            },
        })
    }
}

#[cfg(test)]
pub(in crate::cavity_phs::piston::double_acting) mod tests;
