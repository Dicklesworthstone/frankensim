//! Endpoint electrical work, reciprocal conduction and original gas chemistry.

use super::{
    Candidate, ConductionError, FiniteGasElectricalFeedbackFrame, FiniteGasElectricalFeedbackSolve,
    FiniteGasStringHeat, FiniteGasStringSolve, HeatExchangeSolve, MaterialConductionLayers,
    ThermalMaterialStringRuntime,
};
use crate::acoustic_realize::finite_gas::flow::FiniteGasStringFlowFrame;
use crate::acoustic_realize::finite_gas::reactive::{
    ReactiveFiniteGasStringSolve, exchange::ReactiveFiniteGasStringExchangeFrame,
    flow::admit_reacting_inventory,
};
use crate::cavity_phs::piston::reactive::GasReaction;
use fs_material::conductor::OhmicDrive;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::thermochemical_gas::energy::IsochoricGasParcel;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::{
    GasNozzle, GasNozzleStepSolve,
};

type ReactiveConductionFrame = ReactiveFiniteGasStringExchangeFrame<ResolvedMaterialStatePoint>;

impl ThermalMaterialStringRuntime {
    /// Close electrical work, solid/gas conduction and heated chemistry on the
    /// same accepted endpoint. The original reaction owner supplies composition
    /// and temperature for every conduction proposal, and the original
    /// conductor supplies endpoint current/resistance for every Joule trial.
    ///
    /// All trials restart from the accepted string and species inventory. The
    /// bounded source equation remains J = dt I_end^2 R_end; chemical formation
    /// energy is already in gas storage and is never another heat input. The
    /// original chemical, caloric, conduction, source and work admissions remain
    /// authoritative. Failure publishes neither owner nor its sound sample.
    ///
    /// This remains a fixed-volume, well-mixed gas with frozen-initial mechanics,
    /// a fully solid conductor and diagnostic compact pressure. No flame,
    /// circuit storage, surface reaction or acoustic back-reaction is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_implicit_drive_and_reacting_isochoric_gas_conduction(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        drive: OhmicDrive,
        gas: &mut IsochoricGasParcel,
        heat: FiniteGasStringHeat,
        reaction: &GasReaction,
        controls: ReactiveFiniteGasStringSolve,
        exchange: HeatExchangeSolve,
        layers: &MaterialConductionLayers<'_>,
        electrical: FiniteGasElectricalFeedbackSolve,
    ) -> Result<FiniteGasElectricalFeedbackFrame<ReactiveConductionFrame>, ConductionError> {
        let candidate = self.propose_implicit_reacting_conduction(
            cx,
            expected_epoch,
            forces,
            drive,
            gas,
            heat,
            reaction,
            controls,
            exchange,
            layers,
            electrical,
        )?;
        self.publish_implicit_thermal(cx, gas, candidate)
    }

    /// Advance that same converged electrical/chemical/conduction endpoint and
    /// one conservative nozzle transfer over the original sample interval.
    /// Flow runs once AFTER feedback convergence, so its species and upstream
    /// enthalpy are never withdrawn once per nonlinear trial.
    ///
    /// The common flow owner audits actual source work against both extensive
    /// gas inventories. Before any owner publishes, the final surrounding
    /// mixture must still use the original reaction catalog, caloric/transport
    /// sources and explicit water identity. Incoming unlisted or incompatible
    /// species therefore refuse atomically instead of stranding the next step.
    ///
    /// Chemistry/conduction then flow is first-order splitting, not a
    /// simultaneous chemical transport solve. The chemistry receipt describes
    /// the pre-flow endpoint; accepted gas and sound use the final mixture.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_implicit_drive_and_reacting_isochoric_gas_nozzle_conduction(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        drive: OhmicDrive,
        gas: &mut IsochoricGasParcel,
        reservoir: &mut ClosedGasParcel,
        heat: FiniteGasStringHeat,
        reaction: &GasReaction,
        controls: ReactiveFiniteGasStringSolve,
        exchange: HeatExchangeSolve,
        layers: &MaterialConductionLayers<'_>,
        nozzle: GasNozzle,
        flow_controls: GasNozzleStepSolve,
        electrical: FiniteGasElectricalFeedbackSolve,
    ) -> Result<
        FiniteGasElectricalFeedbackFrame<FiniteGasStringFlowFrame<ReactiveConductionFrame>>,
        ConductionError,
    > {
        let candidate = self.propose_implicit_reacting_conduction(
            cx,
            expected_epoch,
            forces,
            drive,
            gas,
            heat,
            reaction,
            controls,
            exchange,
            layers,
            electrical,
        )?;
        self.publish_implicit_nozzle(
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
            electrical,
            candidate,
            |thermal| &mut thermal.coupled.coupled,
            |gas, _| admit_reacting_inventory(reaction, gas).map_err(ConductionError::from),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn propose_implicit_reacting_conduction(
        &self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        drive: OhmicDrive,
        gas: &IsochoricGasParcel,
        heat: FiniteGasStringHeat,
        reaction: &GasReaction,
        controls: ReactiveFiniteGasStringSolve,
        exchange: HeatExchangeSolve,
        layers: &MaterialConductionLayers<'_>,
        electrical: FiniteGasElectricalFeedbackSolve,
    ) -> Result<Candidate<ReactiveConductionFrame>, ConductionError> {
        self.propose_implicit_thermal(
            cx,
            expected_epoch,
            drive,
            gas,
            heat,
            electrical,
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
            |thermal| &thermal.coupled.coupled,
        )
    }
}
