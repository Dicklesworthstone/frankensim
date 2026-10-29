//! One accepted chemical/caloric/conduction endpoint around continuing sound.

use super::super::exchange::{
    FiniteGasStringExchangeError, FiniteGasStringExchangeFrame, resolve_conduction,
};
use super::{
    FiniteGasStringError, FiniteGasStringHeat, FiniteGasStringSolve, GasReaction,
    HeatedGasReactionStep, IsochoricGasParcel, ReactiveFiniteGasStringSolve,
    ThermalMaterialStringRuntime,
};
use crate::cavity_phs::piston::thermal::{
    ConductanceSample, HeatExchangeSolve,
    layers::{LayerConductionError, MaterialConductionLayers},
};
use fs_material::phase::EquilibriumPhaseState;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;

/// The original reciprocal exchange and chemistry reports at one endpoint.
#[derive(Debug)]
pub struct ReactiveFiniteGasStringExchangeFrame<R> {
    /// The internal solid-to-gas transfer cancels from this combined balance.
    pub coupled: FiniteGasStringExchangeFrame<R>,
    /// Chemistry uses the heat from this SAME accepted exchange proposal.
    pub chemistry: HeatedGasReactionStep,
}

impl ThermalMaterialStringRuntime {
    /// Solve reciprocal thermal exchange and gas chemistry at one endpoint.
    ///
    /// Each proposed transfer Q removes Q from the original solid and supplies
    /// Q plus aerodynamic/external heat to the original heated implicit reaction
    /// owner. Every trial restarts from the accepted species inventory. The
    /// existing bounded exchange residual uses those reacting endpoint states:
    /// Q = dt G(T_s,end,T_g,end,n_end)(T_s,end-T_g,end).
    ///
    /// The SAME thermal proposal owner also closes any admitted reversible
    /// thermoelastic heat. Its source/temperature/error/work bounds remain in
    /// force. The native chemical work ceiling applies per endpoint evaluation;
    /// exchange and reversible-heat work retain their own bounded trial counts.
    /// Formation energy is inside gas storage and never an extra source. Final
    /// composition drives both interface conductance and drag/compact pressure.
    ///
    /// This remains first-order frozen-initial mechanics with a well-mixed
    /// fixed-volume gas, not a flame, chamber or resolved boundary-layer model.
    /// `resolve` is repeatable, proposal-only and sees actual endpoint states.
    /// Any source, chemical, thermal, observer, work or cancellation refusal
    /// publishes neither participant nor its next sound sample.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_reacting_isochoric_gas_exchange<R, E>(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        gas: &mut IsochoricGasParcel,
        heat: FiniteGasStringHeat,
        reaction: &GasReaction,
        controls: ReactiveFiniteGasStringSolve,
        exchange: HeatExchangeSolve,
        resolve: impl FnMut(
            &ResolvedThermochemicalGasState,
            EquilibriumPhaseState,
        ) -> Result<ConductanceSample<R>, E>,
    ) -> Result<ReactiveFiniteGasStringExchangeFrame<R>, FiniteGasStringExchangeError<E>> {
        let finite = FiniteGasStringSolve {
            gas: controls.reaction.gas,
            maximum_gas_temperature_increment_k: controls.maximum_gas_temperature_increment_k,
        };
        let (coupled, chemistry) = self.step_with_isochoric_gas_exchange_endpoint(
            cx,
            expected_epoch,
            forces,
            gas,
            heat,
            finite,
            exchange,
            resolve,
            |initial, heat_j, duration_s, cancelled| {
                let (gas, report) = reaction
                    .heated_implicit_candidate(
                        &ClosedGasParcel::from_isochoric(initial.clone()),
                        duration_s,
                        heat_j,
                        controls.reaction,
                        cancelled,
                    )
                    .map_err(FiniteGasStringError::from)?
                    .into_parts();
                Ok((gas.into_isochoric(), report.heat.clone(), report))
            },
        )?;
        Ok(ReactiveFiniteGasStringExchangeFrame { coupled, chemistry })
    }

    /// Source-bound series conduction on the reacting endpoint above.
    /// Uses the original solid-card/condition, surface-area, source-selection
    /// and Biot admissions. Gas conductivity comes from the final reacting
    /// mixture at EVERY trial; it is not sampled before chemical evolution.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_reacting_isochoric_gas_conduction(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        gas: &mut IsochoricGasParcel,
        heat: FiniteGasStringHeat,
        reaction: &GasReaction,
        controls: ReactiveFiniteGasStringSolve,
        exchange: HeatExchangeSolve,
        layers: &MaterialConductionLayers<'_>,
    ) -> Result<
        ReactiveFiniteGasStringExchangeFrame<ResolvedMaterialStatePoint>,
        FiniteGasStringExchangeError<LayerConductionError>,
    > {
        let volume = self.admit_finite_gas_conduction(cx, layers)?;
        self.step_with_reacting_isochoric_gas_exchange(
            cx,
            expected_epoch,
            forces,
            gas,
            heat,
            reaction,
            controls,
            exchange,
            |gas_end, solid_end| resolve_conduction(layers, volume, gas_end, solid_end),
        )
    }
}
