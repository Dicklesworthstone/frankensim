//! Electrical work, material vibration and finite-gas heating at one epoch.
//!
//! The existing conductor supplies initial resistance/current and Joule work.
//! The existing finite-gas step owns vibration, dissipation, enthalpy, gas EOS,
//! heat exchange and sound. Both proposals remain private until the endpoint
//! electrical law is also admitted. No circuit or transport law is duplicated.

/// Bounded endpoint Joule-work feedback on the same thermal/flow transactions.
pub mod feedback;

use super::exchange::{FiniteGasStringExchangeError, FiniteGasStringExchangeFrame};
use super::flow::FiniteGasStringFlowFrame;
use super::{
    FiniteGasStringError, FiniteGasStringFrame, FiniteGasStringHeat, FiniteGasStringSolve,
};
use crate::acoustic_realize::{ThermalMaterialStringRuntime, ThermalStringError};
use crate::cavity_phs::piston::thermal::{
    HeatExchangeSolve,
    layers::{LayerConductionError, MaterialConductionLayers},
};
use fs_material::conductor::{OhmicDrive, ResolvedConductor};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::thermochemical_gas::energy::IsochoricGasParcel;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::{
    GasNozzle, GasNozzleStepSolve,
};

/// Original coupled report and the electrical source work it already includes.
#[derive(Debug)]
pub struct FiniteGasElectrothermalStringFrame<R = FiniteGasStringFrame> {
    /// Original finite-gas or reciprocal-conduction report. Its solid source
    /// includes `external_string_heat_j + joule_heat_j`; its combined energy
    /// balance therefore already counts electrical source work exactly once.
    pub coupled: R,
    /// Ideal source prescription, retained separately from resolved current.
    pub drive: OhmicDrive,
    /// Initial current [A]. The explicit route freezes this value for source
    /// work; the feedback route uses `final_current_a` for its boundary work.
    pub current_a: f64,
    /// Current demanded by the same source at accepted endpoint resistance [A].
    pub final_current_a: f64,
    /// Electrical work deposited once as solid heat [J]. Endpoint feedback
    /// separately reports actual source work and its bounded residual.
    pub joule_heat_j: f64,
    /// Caller-prescribed solid heat only [J], excluding Joule/dissipation/exchange.
    pub external_string_heat_j: f64,
    /// Original source-resolved conductor at the initial accepted temperature.
    pub conductor_before: ResolvedConductor,
    /// Original source-resolved conductor at the final accepted temperature.
    pub conductor_after: ResolvedConductor,
}

impl<R> FiniteGasElectrothermalStringFrame<R> {
    fn try_map_coupled<T, E>(
        self,
        map: impl FnOnce(R) -> Result<T, E>,
    ) -> Result<FiniteGasElectrothermalStringFrame<T>, E> {
        Ok(FiniteGasElectrothermalStringFrame {
            coupled: map(self.coupled)?,
            drive: self.drive,
            current_a: self.current_a,
            final_current_a: self.final_current_a,
            joule_heat_j: self.joule_heat_j,
            external_string_heat_j: self.external_string_heat_j,
            conductor_before: self.conductor_before,
            conductor_after: self.conductor_after,
        })
    }
}

impl ThermalMaterialStringRuntime {
    /// Evolve an electrically heated, vibrating string and finite surrounding gas.
    ///
    /// Joule work is `I_initial² R_initial dt` from the original conductor owner;
    /// a voltage source resolves `I_initial = V/R_initial`. It is added once to
    /// prescribed solid heat. Material viscous loss heats the solid; aerodynamic
    /// loss heats the gas. The original finite-gas report retains those internal
    /// transfers and its independent combined energy check.
    ///
    /// Endpoint resistance/current, thermal/material state, gas and observer all
    /// publish at one epoch. A late source refusal or cancellation publishes
    /// neither owner. Inputs are not consumed on failure.
    ///
    /// This insulated solid/gas interface has no reciprocal thermal exchange;
    /// use `step_with_drive_and_isochoric_gas_conduction` for a sourced layer.
    /// Electrical and mechanical coefficients remain frozen at the initial
    /// state during one sample: this is first-order partitioned coupling. A
    /// thermoelastic string then updates stiffness/prestress through its original
    /// material laws and reversible heat account, with invariant reference mass.
    /// No circuit dynamics, contact resistance, skin effect, Lorentz force,
    /// pressure work, radiation loading or acoustic chamber is supplied.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_drive_and_isochoric_gas(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        drive: OhmicDrive,
        gas: &mut IsochoricGasParcel,
        heat: FiniteGasStringHeat,
        controls: FiniteGasStringSolve,
    ) -> Result<FiniteGasElectrothermalStringFrame, FiniteGasStringError> {
        self.step_with_drive_and_finite_gas(
            cx,
            expected_epoch,
            drive,
            gas,
            heat,
            |candidate, gas, heat| {
                candidate.step_with_isochoric_gas(cx, expected_epoch, forces, gas, heat, controls)
            },
        )
    }

    /// Electrically heated counterpart of the existing source-bound reciprocal
    /// conduction step. The SAME solid card supplies mechanics, resistivity,
    /// caloric storage and conductivity; the finite source gas supplies its own
    /// energy, EOS and conductivity. Actual endpoint temperatures determine the
    /// original layer's signed transfer, and its geometry/Biot/source gates are
    /// unchanged. Numerical loss and diagnostic radiation are never Joule heat.
    ///
    /// The nested exchange report's `external_string_heat_j` includes Joule work
    /// plus the caller's prescribed solid heat. This wrapper separates them for
    /// reporting; do not add `joule_heat_j` again to its combined energy balance.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_drive_and_isochoric_gas_conduction(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        drive: OhmicDrive,
        gas: &mut IsochoricGasParcel,
        heat: FiniteGasStringHeat,
        controls: FiniteGasStringSolve,
        exchange: HeatExchangeSolve,
        layers: &MaterialConductionLayers<'_>,
    ) -> Result<
        FiniteGasElectrothermalStringFrame<
            FiniteGasStringExchangeFrame<ResolvedMaterialStatePoint>,
        >,
        FiniteGasStringExchangeError<LayerConductionError>,
    > {
        self.step_with_drive_and_finite_gas(
            cx,
            expected_epoch,
            drive,
            gas,
            heat,
            |candidate, gas, heat| {
                candidate.step_with_isochoric_gas_conduction(
                    cx,
                    expected_epoch,
                    forces,
                    gas,
                    heat,
                    controls,
                    exchange,
                    layers,
                )
            },
        )
    }

    /// Electrical work, reciprocal material conduction and conservative
    /// atmospheric exchange on one continuing string. The source-resolved
    /// conductor supplies Joule work once; the original conduction and nozzle
    /// owners advance over the same sample duration. The nested flow report
    /// closes the combined solid/mechanical/two-gas energy balance.
    ///
    /// Endpoint electrical admission happens AFTER both thermal and transport
    /// proposals. A late conductor refusal or cancellation leaves the string,
    /// local gas, reservoir and observer history at their previous accepted
    /// states. No source or transferred inventory is consumed on retry.
    /// Electrical/mechanical coefficients are frozen initially and heat precedes
    /// flow; this retains the original first-order partitioned model.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_drive_and_isochoric_gas_nozzle_conduction(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        drive: OhmicDrive,
        gas: &mut IsochoricGasParcel,
        reservoir: &mut ClosedGasParcel,
        heat: FiniteGasStringHeat,
        controls: FiniteGasStringSolve,
        exchange: HeatExchangeSolve,
        layers: &MaterialConductionLayers<'_>,
        nozzle: GasNozzle,
        flow_controls: GasNozzleStepSolve,
    ) -> Result<
        FiniteGasElectrothermalStringFrame<
            FiniteGasStringFlowFrame<FiniteGasStringExchangeFrame<ResolvedMaterialStatePoint>>,
        >,
        FiniteGasStringExchangeError<LayerConductionError>,
    > {
        let mut next_reservoir = reservoir.clone();
        let (candidate, next_gas, frame) = self.propose_drive_and_finite_gas(
            cx,
            expected_epoch,
            drive,
            gas,
            heat,
            |candidate, gas, heat| {
                candidate.step_with_isochoric_gas_nozzle_conduction(
                    cx,
                    expected_epoch,
                    forces,
                    gas,
                    &mut next_reservoir,
                    heat,
                    controls,
                    exchange,
                    layers,
                    nozzle,
                    flow_controls,
                )
            },
        )?;
        super::check(cx)?;
        // Retain all old owners until the complete exclusive-state publication.
        let _old_string = core::mem::replace(self, candidate);
        let _old_gas = core::mem::replace(gas, next_gas);
        let _old_reservoir = core::mem::replace(reservoir, next_reservoir);
        Ok(frame)
    }

    fn step_with_drive_and_finite_gas<R, E>(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        drive: OhmicDrive,
        gas: &mut IsochoricGasParcel,
        heat: FiniteGasStringHeat,
        advance: impl FnOnce(&mut Self, &mut IsochoricGasParcel, FiniteGasStringHeat) -> Result<R, E>,
    ) -> Result<FiniteGasElectrothermalStringFrame<R>, E>
    where
        E: From<ThermalStringError> + From<FiniteGasStringError>,
    {
        let (candidate, next_gas, frame) =
            self.propose_drive_and_finite_gas(cx, expected_epoch, drive, gas, heat, advance)?;
        super::check(cx)?;
        let _old_string = core::mem::replace(self, candidate);
        let _old_gas = core::mem::replace(gas, next_gas);
        Ok(frame)
    }

    // The same electrical proposal can accompany two- or three-owner thermal
    // transactions. Nothing accepted is published before the final conductor
    // law and current have been resolved successfully.
    fn propose_drive_and_finite_gas<R, E>(
        &self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        drive: OhmicDrive,
        gas: &IsochoricGasParcel,
        heat: FiniteGasStringHeat,
        advance: impl FnOnce(&mut Self, &mut IsochoricGasParcel, FiniteGasStringHeat) -> Result<R, E>,
    ) -> Result<
        (
            Self,
            IsochoricGasParcel,
            FiniteGasElectrothermalStringFrame<R>,
        ),
        E,
    >
    where
        E: From<ThermalStringError> + From<FiniteGasStringError>,
    {
        let (conductor_before, current_a, joule_heat_j) =
            self.prepare_drive(cx, expected_epoch, drive)?;
        self.propose_joule_work_and_finite_gas(
            drive,
            gas,
            heat,
            conductor_before,
            current_a,
            joule_heat_j,
            advance,
        )
    }

    // Initial-law and endpoint-law sources use the same deposition, thermal
    // proposal and final conductor admission. Only work selection differs.
    #[allow(clippy::too_many_arguments)]
    fn propose_joule_work_and_finite_gas<R, E>(
        &self,
        drive: OhmicDrive,
        gas: &IsochoricGasParcel,
        heat: FiniteGasStringHeat,
        conductor_before: ResolvedConductor,
        current_a: f64,
        joule_heat_j: f64,
        advance: impl FnOnce(&mut Self, &mut IsochoricGasParcel, FiniteGasStringHeat) -> Result<R, E>,
    ) -> Result<
        (
            Self,
            IsochoricGasParcel,
            FiniteGasElectrothermalStringFrame<R>,
        ),
        E,
    >
    where
        E: From<ThermalStringError> + From<FiniteGasStringError>,
    {
        let string_heat_j = heat.string_heat_j + joule_heat_j;
        if !heat.string_heat_j.is_finite()
            || !string_heat_j.is_finite()
            || (joule_heat_j != 0.0 && string_heat_j == heat.string_heat_j)
            || (heat.string_heat_j != 0.0 && string_heat_j == joule_heat_j)
        {
            return Err(FiniteGasStringError::Admission(
                "Joule and prescribed solid heat overflow or disappear when combined",
            )
            .into());
        }
        let mut candidate = self.clone();
        let mut next_gas = gas.clone();
        let coupled = advance(
            &mut candidate,
            &mut next_gas,
            FiniteGasStringHeat {
                string_heat_j,
                gas_heat_j: heat.gas_heat_j,
            },
        )?;
        let conductor_after = candidate.conductor()?;
        let final_current_a = conductor_after
            .current_a_for_drive(drive)
            .map_err(ThermalStringError::Electrical)?;
        Ok((
            candidate,
            next_gas,
            FiniteGasElectrothermalStringFrame {
                coupled,
                drive,
                current_a,
                final_current_a,
                joule_heat_j,
                external_string_heat_j: heat.string_heat_j,
                conductor_before,
                conductor_after,
            },
        ))
    }
}
