//! A changing source-backed atmosphere on one continuing material string.
//!
//! The existing acoustic/solid/gas step first advances once. The original
//! quasi-steady nozzle then exchanges species and upstream enthalpy with a
//! second finite, well-mixed volume. Final drag and pressure observation use
//! that same conserved mixture. No new gas, flow or caloric law is supplied.

use super::exchange::{FiniteGasStringExchangeError, FiniteGasStringExchangeFrame};
use super::{
    FiniteGasStringError, FiniteGasStringFrame, FiniteGasStringHeat, FiniteGasStringSolve,
    ThermalMaterialStringRuntime, check,
};
use crate::cavity_phs::piston::thermal::{
    HeatExchangeSolve,
    layers::{LayerConductionError, MaterialConductionLayers},
};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::thermochemical_gas::energy::IsochoricGasParcel;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::GasMassTransferStep;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::{
    GasNozzle, GasNozzleFlow, GasNozzleStepSolve,
};

/// Direction selected from the actual pressures after the heating substep.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FiniteGasStringFlowDirection {
    /// The other finite volume donates to the string's surroundings.
    IntoSurroundings,
    /// The string's surroundings donate to the other finite volume.
    OutOfSurroundings,
}

/// One sample of solid motion, heat and conservative atmospheric exchange.
#[derive(Debug)]
pub struct FiniteGasStringFlowFrame<R = FiniteGasStringFrame> {
    /// Original thermal proposal, including any reciprocal exchange/source
    /// receipt. Gas heating and conductance describe the PRE-transfer endpoint;
    /// the contained string observer uses the FINAL conserved mixture. Thermal
    /// and transported enthalpy accounts are not additional external sources.
    pub thermal: R,
    /// Upstream/downstream ownership for `flow` and `transfer`. Equal pressures
    /// select `IntoSurroundings` with zero transferred mass.
    pub direction: FiniteGasStringFlowDirection,
    /// Original source-resolved flow, evaluated before inventory withdrawal.
    pub flow: GasNozzleFlow,
    /// Original species and enthalpy transfer. The jet's total enthalpy already
    /// includes flow work and thermalized kinetic energy; neither is extra heat.
    pub transfer: GasMassTransferStep,
    /// Actual change of local gas `m*u`, including heat AND exchange [J].
    pub surrounding_energy_change_j: f64,
    /// Actual change of the other gas volume's `m*u` [J].
    pub reservoir_energy_change_j: f64,
    /// External work/heats minus BOTH gas, solid and mechanical storage changes,
    /// authored loss and explicit damping-partition roundoff [J].
    pub energy_balance_residual_j: f64,
    /// Original thermal/transfer allowances plus arithmetic roundoff [J].
    pub energy_tolerance_j: f64,
}

impl ThermalMaterialStringRuntime {
    /// Advance one string sample while its finite atmosphere exchanges gas.
    ///
    /// The initial surrounding snapshot must match the runtime exactly. The
    /// existing finite-gas step advances motion and solid/gas heat once at the
    /// initial coefficients. The higher-pressure gas then supplies the native
    /// nozzle's upstream state for exactly that sample duration. Its original
    /// species/source join and extensive energy inversions reconstruct both
    /// endpoints; donor/receiver solve controls follow the selected direction.
    ///
    /// Drag and the compact observer are rebuilt from the final local mixture,
    /// retaining modal displacement/momentum and the PREVIOUS accepted sample's
    /// acceleration. The next sample therefore evolves with the new density,
    /// viscosity and sound speed. One epoch and one audio sample are consumed.
    /// Source disagreement, excessive withdrawal, pressure-order reversal,
    /// final temperature bounds, observer failure or cancellation publish NONE
    /// of the three owners. No source input is consumed by a rejected proposal.
    ///
    /// This is a first-order partitioned model with two well-mixed fixed gas
    /// volumes and a quasi-steady converging nozzle. It needs physical timestep
    /// refinement; the withdrawal/temperature bounds are admission limits, not
    /// an accuracy estimate. It supplies no resolved jet momentum, wind forcing,
    /// solid-gas conductive exchange, acoustic chamber, radiation loading or
    /// propagation through a spatially/time-varying medium. The compact pressure
    /// signal remains diagnostic. Transported enthalpy is internal to the pair,
    /// never an additional source of heat or acoustic energy.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_isochoric_gas_nozzle(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        gas: &mut IsochoricGasParcel,
        reservoir: &mut ClosedGasParcel,
        heat: FiniteGasStringHeat,
        thermal_controls: FiniteGasStringSolve,
        nozzle: GasNozzle,
        flow_controls: GasNozzleStepSolve,
    ) -> Result<FiniteGasStringFlowFrame, FiniteGasStringError> {
        self.step_with_finite_gas_nozzle(
            cx,
            expected_epoch,
            forces,
            gas,
            reservoir,
            heat,
            thermal_controls,
            nozzle,
            flow_controls,
            0.0,
            0.0,
            |candidate, gas, heat| {
                candidate.step_with_isochoric_gas(
                    cx,
                    expected_epoch,
                    forces,
                    gas,
                    heat,
                    thermal_controls,
                )
            },
            |thermal| thermal,
        )
    }

    /// Evolve reciprocal material/gas conduction and pressure-driven exchange
    /// in one accepted sample. The original source-bound thermal solve runs
    /// first; the original nozzle then transports species and upstream enthalpy
    /// over the SAME interval. Solid heat, gas heat and mechanical loss retain
    /// their original owners, while the final medium drives drag and pressure.
    ///
    /// Layer source/geometry/Biot admissions are unchanged. The typed thermal
    /// report retains its pre-flow endpoint and signed internal transfer. This
    /// is first-order heat-then-flow splitting; timestep refinement is required
    /// for simultaneous heat/transport dynamics. It supplies no resolved jet
    /// momentum, chamber resonance or time-varying propagation path.
    /// Failure after either proposal preserves all three accepted owners.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_isochoric_gas_nozzle_conduction(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        gas: &mut IsochoricGasParcel,
        reservoir: &mut ClosedGasParcel,
        heat: FiniteGasStringHeat,
        thermal_controls: FiniteGasStringSolve,
        exchange: HeatExchangeSolve,
        layers: &MaterialConductionLayers<'_>,
        nozzle: GasNozzle,
        flow_controls: GasNozzleStepSolve,
    ) -> Result<
        FiniteGasStringFlowFrame<FiniteGasStringExchangeFrame<ResolvedMaterialStatePoint>>,
        FiniteGasStringExchangeError<LayerConductionError>,
    > {
        self.step_with_finite_gas_nozzle(
            cx,
            expected_epoch,
            forces,
            gas,
            reservoir,
            heat,
            thermal_controls,
            nozzle,
            flow_controls,
            0.0,
            0.0,
            |candidate, gas, heat| {
                candidate.step_with_isochoric_gas_conduction(
                    cx,
                    expected_epoch,
                    forces,
                    gas,
                    heat,
                    thermal_controls,
                    exchange,
                    layers,
                )
            },
            |thermal| &mut thermal.coupled,
        )
    }

    // One flow/accounting/publication owner for every thermal proposal. The
    // accessor only exposes the original shared finite-gas account; extension
    // reports and their source receipts stay attached to the same proposal.
    // An implicit source may distinguish applied heat from actual boundary
    // work by its admitted residual. That correction is audited here, before
    // publication, instead of renormalizing either gas or solid energy.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn step_with_finite_gas_nozzle<R, E>(
        &mut self,
        cx: &fs_exec::Cx<'_>,
        expected_epoch: u64,
        forces: &[f64],
        gas: &mut IsochoricGasParcel,
        reservoir: &mut ClosedGasParcel,
        heat: FiniteGasStringHeat,
        thermal_controls: FiniteGasStringSolve,
        nozzle: GasNozzle,
        flow_controls: GasNozzleStepSolve,
        source_work_correction_j: f64,
        source_work_tolerance_j: f64,
        advance: impl FnOnce(&mut Self, &mut IsochoricGasParcel, FiniteGasStringHeat) -> Result<R, E>,
        account: impl FnOnce(&mut R) -> &mut FiniteGasStringFrame,
    ) -> Result<FiniteGasStringFlowFrame<R>, E>
    where
        E: From<FiniteGasStringError>,
    {
        if !source_work_correction_j.is_finite()
            || !source_work_tolerance_j.is_finite()
            || source_work_tolerance_j < 0.0
            || source_work_correction_j.abs() > source_work_tolerance_j
        {
            return Err(FiniteGasStringError::Admission(
                "source work differs from deposited heat outside its finite residual allowance",
            )
            .into());
        }
        self.admit_finite_gas_step(cx, expected_epoch, gas, heat, thermal_controls)?;
        let mut candidate = self.clone();
        let mut heated_gas = gas.clone();
        let mut thermal = advance(&mut candidate, &mut heated_gas, heat)?;
        let heated = account(&mut thermal);
        let local = ClosedGasParcel::from_isochoric(heated_gas);
        let into = reservoir.state().state().pressure >= local.state().state().pressure;
        let (donor, receiver) = if into {
            (&*reservoir, &local)
        } else {
            (&local, &*reservoir)
        };
        let (proposal, flow) = donor
            .through_nozzle_to(
                receiver,
                nozzle,
                self.mechanical.member.model.sample_period_s(),
                flow_controls,
                || cx.checkpoint().is_err(),
            )
            .map_err(FiniteGasStringError::from)?;
        let (next_donor, next_receiver, transfer) = proposal.into_parts();
        let (next_gas, next_reservoir, direction) = if into {
            (
                next_receiver.into_isochoric(),
                next_donor,
                FiniteGasStringFlowDirection::IntoSurroundings,
            )
        } else {
            (
                next_donor.into_isochoric(),
                next_receiver,
                FiniteGasStringFlowDirection::OutOfSurroundings,
            )
        };
        if (next_gas.state().state().temperature - gas.state().state().temperature).abs()
            > thermal_controls.maximum_gas_temperature_increment_k
        {
            return Err(FiniteGasStringError::Admission(
                "gas heating and species exchange exceeded the per-sample temperature-increment limit",
            ).into());
        }
        self.refresh_finite_gas_medium(
            cx,
            forces,
            &mut candidate,
            &mut heated.string,
            next_gas.state().state(),
        )?;

        // Measure extensive endpoint storage, including the changed masses.
        // Differencing specific energy first avoids discarding small transfers
        // against a large, valid thermochemical formation-energy reference.
        let gas_u = gas.state().specific_internal_energy_j_kg();
        let next_gas_u = next_gas.state().specific_internal_energy_j_kg();
        let reservoir_u = reservoir.state().specific_internal_energy_j_kg();
        let next_reservoir_u = next_reservoir.state().specific_internal_energy_j_kg();
        let surrounding_energy_change_j = next_gas.mass_kg() * (next_gas_u - gas_u)
            + (next_gas.mass_kg() - gas.mass_kg()) * gas_u;
        let reservoir_energy_change_j = next_reservoir.mass_kg() * (next_reservoir_u - reservoir_u)
            + (next_reservoir.mass_kg() - reservoir.mass_kg()) * reservoir_u;
        let mechanical_before = self
            .mechanical_internal_energy_j()
            .map_err(FiniteGasStringError::from)?;
        let mechanical_after = candidate
            .mechanical_internal_energy_j()
            .map_err(FiniteGasStringError::from)?;
        let mechanical_storage = mechanical_after - mechanical_before;
        let solid_storage = self.mechanical.specimen().mass_kg()
            * (candidate.thermal.specific_enthalpy_j_kg() - self.thermal.specific_enthalpy_j_kg());
        let input =
            heated.string.vibration.acoustic.input_work_j + heat.string_heat_j + heat.gas_heat_j;
        let input = if source_work_correction_j == 0.0 {
            input
        } else {
            input + source_work_correction_j
        };
        let dissipation = &heated.string.vibration.dissipation;
        let residual = input
            - (mechanical_storage
                + solid_storage
                + surrounding_energy_change_j
                + reservoir_energy_change_j
                + dissipation.authored_loss_j
                + dissipation.roundoff_residual_j);
        let scale = input.abs()
            + mechanical_before.abs()
            + mechanical_after.abs()
            + solid_storage.abs()
            + surrounding_energy_change_j.abs()
            + reservoir_energy_change_j.abs()
            + gas.mass_kg() * gas_u.abs()
            + next_gas.mass_kg() * next_gas_u.abs()
            + reservoir.mass_kg() * reservoir_u.abs()
            + next_reservoir.mass_kg() * next_reservoir_u.abs();
        let tolerance = heated.energy_tolerance_j
            + flow_controls.transfer.energy_tolerance_j
            + source_work_tolerance_j
            + 32.0 * f64::EPSILON * scale;
        if !residual.is_finite() || !tolerance.is_finite() || residual.abs() > tolerance {
            return Err(FiniteGasStringError::Admission(
                "independently measured string/solid/two-volume gas energy balance did not close",
            )
            .into());
        }
        check(cx)?;
        // All three exclusive owners publish before any old owner is dropped.
        let _old_string = core::mem::replace(self, candidate);
        let _old_gas = core::mem::replace(gas, next_gas);
        let _old_reservoir = core::mem::replace(reservoir, next_reservoir);
        Ok(FiniteGasStringFlowFrame {
            thermal,
            direction,
            flow,
            transfer,
            surrounding_energy_change_j,
            reservoir_energy_change_j,
            energy_balance_residual_j: residual,
            energy_tolerance_j: tolerance,
        })
    }
}
