//! One continuing material string, reacting atmosphere and finite gas reservoir.
//!
//! Chemistry and reciprocal heat exchange retain their original endpoint
//! owners. The shared nozzle transaction transports the resulting species and
//! upstream enthalpy once, over the same sample interval, then rebuilds drag
//! and diagnostic pressure from that final mixture.

use super::exchange::ReactiveFiniteGasStringExchangeFrame;
use super::{
    FiniteGasStringError, FiniteGasStringHeat, FiniteGasStringSolve, GasReaction,
    ReactiveFiniteGasStringFrame, ReactiveFiniteGasStringSolve, ThermalMaterialStringRuntime,
};
use crate::acoustic_realize::finite_gas::exchange::FiniteGasStringExchangeError;
use crate::acoustic_realize::finite_gas::flow::FiniteGasStringFlowFrame;
use crate::cavity_phs::piston::thermal::{
    HeatExchangeSolve,
    layers::{LayerConductionError, MaterialConductionLayers},
};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::thermochemical_gas::energy::IsochoricGasParcel;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::{
    GasNozzle, GasNozzleStepSolve,
};

impl ThermalMaterialStringRuntime {
    /// Evolve chemistry, sound and conservative finite-volume gas exchange.
    ///
    /// The original heated implicit reaction owner first advances the local
    /// composition and caloric state with actual aerodynamic and prescribed
    /// heat. The original nozzle then transfers species and upstream enthalpy
    /// between two finite volumes over that SAME sample interval. No species,
    /// rate, conductivity, caloric or nozzle law is reconstructed here.
    ///
    /// The reaction catalog must include all species and original source
    /// selections present in the final local mixture, including incoming inert
    /// gas. That condition is checked BEFORE publication, so the next sample
    /// can continue chemistry from the transferred inventory. A reservoir need
    /// not react; its complete species and energy accounts are still retained.
    ///
    /// The nested chemistry/heat report describes the pre-flow endpoint. The
    /// nested string observer describes the final conserved mixture. Formation
    /// energy is already inside gas storage, while transported stream enthalpy
    /// is internal to the two-volume account. Neither is another heat source.
    ///
    /// This is first-order reaction/heat-then-flow splitting with frozen-initial
    /// mechanics, two well-mixed volumes and diagnostic compact pressure. It
    /// supplies no resolved flame, jet momentum or acoustic back-reaction.
    /// Failed chemistry, flow, source admission, energy accounting or
    /// cancellation preserves all three accepted owners and the audio cursor.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_reacting_isochoric_gas_nozzle(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        gas: &mut IsochoricGasParcel,
        reservoir: &mut ClosedGasParcel,
        heat: FiniteGasStringHeat,
        reaction: &GasReaction,
        controls: ReactiveFiniteGasStringSolve,
        nozzle: GasNozzle,
        flow_controls: GasNozzleStepSolve,
    ) -> Result<FiniteGasStringFlowFrame<ReactiveFiniteGasStringFrame>, FiniteGasStringError> {
        self.step_with_finite_gas_nozzle(
            cx,
            expected_epoch,
            forces,
            gas,
            reservoir,
            heat,
            FiniteGasStringSolve {
                gas: controls.reaction.gas,
                maximum_gas_temperature_increment_k: controls.maximum_gas_temperature_increment_k,
            },
            nozzle,
            flow_controls,
            0.0,
            0.0,
            |candidate, gas, heat| {
                candidate.step_with_reacting_isochoric_gas(
                    cx,
                    expected_epoch,
                    forces,
                    gas,
                    heat,
                    reaction,
                    controls,
                )
            },
            |thermal| &mut thermal.coupled,
            |gas, _| admit_reacting_inventory(reaction, gas),
        )
    }

    /// The same transaction with source-bound solid/gas thermal conduction.
    ///
    /// Each original conduction proposal solves temperature and chemistry
    /// together; the final reacting mixture determines gas conductivity.
    /// Transfer occurs only after that solve succeeds. The source-bound solid
    /// layer, geometry, Biot and thermal work limits remain authoritative.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_reacting_isochoric_gas_nozzle_conduction(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        gas: &mut IsochoricGasParcel,
        reservoir: &mut ClosedGasParcel,
        heat: FiniteGasStringHeat,
        reaction: &GasReaction,
        controls: ReactiveFiniteGasStringSolve,
        exchange: HeatExchangeSolve,
        layers: &MaterialConductionLayers<'_>,
        nozzle: GasNozzle,
        flow_controls: GasNozzleStepSolve,
    ) -> Result<
        FiniteGasStringFlowFrame<ReactiveFiniteGasStringExchangeFrame<ResolvedMaterialStatePoint>>,
        FiniteGasStringExchangeError<LayerConductionError>,
    > {
        self.step_with_finite_gas_nozzle(
            cx,
            expected_epoch,
            forces,
            gas,
            reservoir,
            heat,
            FiniteGasStringSolve {
                gas: controls.reaction.gas,
                maximum_gas_temperature_increment_k: controls.maximum_gas_temperature_increment_k,
            },
            nozzle,
            flow_controls,
            0.0,
            0.0,
            |candidate, gas, heat| {
                candidate.step_with_reacting_isochoric_gas_conduction(
                    cx,
                    expected_epoch,
                    forces,
                    gas,
                    heat,
                    reaction,
                    controls,
                    exchange,
                    layers,
                )
            },
            |thermal| &mut thermal.coupled.coupled,
            |gas, _| admit_reacting_inventory(reaction, gas).map_err(Into::into),
        )
    }
}

// Reuse the reaction owner's exact species, caloric/transport source and water
// identity admission. Transport may introduce species absent before this step;
// publishing them outside this catalog would strand the next chemical step.
pub(in crate::acoustic_realize::finite_gas) fn admit_reacting_inventory(
    reaction: &GasReaction,
    gas: &IsochoricGasParcel,
) -> Result<(), FiniteGasStringError> {
    reaction.inventory(&ClosedGasParcel::from_isochoric(gas.clone()))?;
    Ok(())
}
