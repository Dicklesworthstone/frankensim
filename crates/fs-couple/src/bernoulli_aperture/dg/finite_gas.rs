//! Finite thermochemical feedback through the original moving aperture and load.
//!
//! The original nonlinear source and complete acoustic load advance once in
//! the accepted initial medium. The existing finite-background work solve then
//! refreshes that same acoustic endpoint. Every trial starts from the original
//! gas and the same marched wave, branch, channel and exterior histories.

use super::{ApertureDg, DgApertureFrame, accept_acoustic_candidate};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::dynamic::ApertureDrive;
use crate::bernoulli_aperture::modal::pressure::MaterialPressureFrame;
use crate::bernoulli_aperture::tube::TubeDrive;
use crate::bernoulli_aperture::two_sided::LoadedAperture;
use crate::bernoulli_aperture::wall::material::heated_finite_gas;
use crate::dg_port::finite_gas::{
    FiniteGasDgBackground, FiniteGasDgError, FiniteGasDgPortCandidate, FiniteGasDgPortStep,
    FiniteGasDgSolve, admit_initial, prepare_after_step,
};
use crate::dg_port::{DgPortLoad, sum_terminal_work};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::composition::GasCompositionStep;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};

/// One accepted interval of the original mechanics, acoustic load and finite gas.
#[derive(Debug)]
pub struct FiniteGasApertureDgStep<S> {
    /// Held-force/flow observations retain their actual interval integrals;
    /// endpoint pressure, storage and parameter work use the final gas.
    pub coupled: DgApertureFrame,
    /// Original caloric/transport report and complete bore/channel intervention.
    /// `gas.acoustic` and `coupled.duct` describe the SAME accepted acoustic
    /// interval; their storage, losses and work must not be added together.
    pub gas: FiniteGasDgPortStep<S>,
    /// Gas plus mechanical/acoustic storage, original losses and signed
    /// exterior exchange, less heat, transported enthalpy and actual drives [J].
    pub combined_energy_residual_j: f64,
}

fn acoustic<E>(error: AcousticRealizeError) -> FiniteGasDgError<E> {
    match error {
        AcousticRealizeError::Cancelled => FiniteGasDgError::Cancelled,
        error => FiniteGasDgError::Acoustic(error),
    }
}

fn admit<A, E, C: FnMut() -> bool>(
    model: &ApertureDg<A>,
    gas: &ClosedGasParcel,
    controls: FiniteGasDgSolve,
    cancelled: &mut C,
) -> Result<(), FiniteGasDgError<E>> {
    if model.upstream.is_some() {
        return Err(FiniteGasDgError::Invalid(
            "finite aperture feedback needs an explicit inventory for each tract; one parcel cannot replace two independent gas owners",
        ));
    }
    admit_initial(&model.load, gas, controls, cancelled)
}

#[allow(clippy::too_many_arguments)]
fn finish<S, E, C: FnMut() -> bool>(
    original: &DgPortLoad,
    marched: &DgPortLoad,
    gas: &ClosedGasParcel,
    mut frame: DgApertureFrame,
    initial_energy_j: f64,
    controls: FiniteGasDgSolve,
    evolve: impl FnMut(&ClosedGasParcel, f64, f64, &mut C) -> Result<FiniteGasDgBackground<S>, E>,
    cancelled: &mut C,
) -> Result<(DgPortLoad, ClosedGasParcel, FiniteGasApertureDgStep<S>), FiniteGasDgError<E>> {
    let FiniteGasDgPortCandidate { load, gas, step } = prepare_after_step(
        original,
        marched,
        frame.duct.clone(),
        gas,
        controls,
        evolve,
        cancelled,
    )?;
    frame.duct = step.acoustic.clone();
    frame.parameter_work_j = step.refresh.parameter_work_j;
    frame.stored_energy_j = frame.aperture.stored_energy_j + load.stored_energy_j();
    frame.storage_change_j = frame.stored_energy_j - initial_energy_j;
    let combined_energy_residual_j = ((step.gas_energy_change_j
        + frame.storage_change_j
        + frame.inferred_total_dissipation_j
        + sum_terminal_work(&frame.terminals))
        - frame.upstream_work_j
        - frame.body_work_j
        - frame.aperture.mechanical_work_j
        - step.supplied_heat_j)
        - step.transported_enthalpy_j;
    if !frame.stored_energy_j.is_finite()
        || !frame.storage_change_j.is_finite()
        || !frame.balance_residual_j().is_finite()
        || frame.balance_residual_j().abs() > controls.energy_tolerance_j
        || !combined_energy_residual_j.is_finite()
        || combined_energy_residual_j.abs() > controls.energy_tolerance_j
    {
        return Err(FiniteGasDgError::Invalid(
            "finite gas, moving aperture and complete acoustic load do not close their shared energy balance",
        ));
    }
    if cancelled() {
        return Err(FiniteGasDgError::Cancelled);
    }
    Ok((
        load,
        gas,
        FiniteGasApertureDgStep {
            coupled: frame,
            gas: step,
            combined_energy_residual_j,
        },
    ))
}

impl ApertureDg {
    /// Advance original scalar mechanics, all attached acoustic owners and the
    /// finite thermochemical background as one transaction.
    ///
    /// The callback receives (original gas, dt, work OUT of gas, cancellation).
    /// It must return a candidate from that original snapshot on EVERY trial,
    /// including its actual heat and separately transported enthalpy. Return
    /// uncommitted companion states in `source`; publish them only on success.
    /// The native finite-gas owner accounts for bore AND terminal-channel work.
    /// The moving source then consumes that same final density, viscosity and
    /// exact load impedance without resetting its mechanical/material history.
    ///
    /// This is the existing first-order frozen-medium partition, with explicit
    /// temperature/work limits; refine its physical time step independently.
    /// Acoustic volume flow does not silently become thermodynamic mass flow.
    /// A mass/species source belongs to the supplied native background process.
    /// No mean-pressure prestress, moving thermodynamic volume or deposition of
    /// acoustic/mechanical dissipation into gas heat is inferred. Two independent
    /// finite tracts require their own joint inventory solve and are refused.
    /// Every refusal leaves gas, valve, load, terminal memory and clocks intact.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_finite_gas_controlled<S, E, C: FnMut() -> bool>(
        &mut self,
        gas: &mut ClosedGasParcel,
        drive: TubeDrive,
        force_n: f64,
        controls: FiniteGasDgSolve,
        evolve: impl FnMut(&ClosedGasParcel, f64, f64, &mut C) -> Result<FiniteGasDgBackground<S>, E>,
        cancelled: &mut C,
    ) -> Result<FiniteGasApertureDgStep<S>, FiniteGasDgError<E>> {
        admit(self, gas, controls, cancelled)?;
        let initial_energy_j = self.stored_energy_j();
        let relation = self.load.port_response().map_err(acoustic)?;
        let mechanical = self
            .aperture
            .preview_step_with_force(
                ApertureDrive {
                    upstream_pressure_pa: drive.upstream_pressure_pa,
                    incoming_pressure_pa: relation.history_pressure_pa / 2.0,
                    body_flow_m3_s: drive.body_flow_m3_s,
                },
                force_n,
            )
            .map_err(acoustic)?;
        let mut marched = self.load.clone();
        let frame = accept_acoustic_candidate(
            &mut marched,
            relation,
            mechanical.frame,
            [MaterialPressureFrame::default(); 2],
            drive,
        )
        .map_err(acoustic)?;
        let (load, next_gas, report) = finish(
            &self.load,
            &marched,
            gas,
            frame,
            initial_energy_j,
            controls,
            evolve,
            cancelled,
        )?;
        // This source update validates everything before its first mutation.
        // No fallible operation or cancellation poll follows publication.
        self.aperture
            .rebind_exact_gas_load(load.duct().inlet_gas(), load.impedance_pa_s_m3())
            .map_err(acoustic)?;
        self.aperture.accept_frame(mechanical);
        self.load = load;
        *gas = next_gas;
        Ok(report)
    }

    /// Original fixed-inventory caloric heating, paying the complete native
    /// acoustic parameter work before any participant publishes its endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_finite_gas_heat_controlled<C: FnMut() -> bool>(
        &mut self,
        gas: &mut ClosedGasParcel,
        drive: TubeDrive,
        force_n: f64,
        heat_j: f64,
        gas_solve: GasEnergySolve,
        controls: FiniteGasDgSolve,
        cancelled: &mut C,
    ) -> Result<FiniteGasApertureDgStep<GasCompositionStep>, FiniteGasDgError<GasEnergyError>> {
        self.step_with_finite_gas_controlled(
            gas,
            drive,
            force_n,
            controls,
            |old, _, work, cancel| heated_finite_gas(old, heat_j, work, gas_solve, cancel),
            cancelled,
        )
    }
}

impl ApertureDg<LoadedAperture> {
    /// The same finite-background transaction for single, paired or material-modal
    /// mechanics. Independent blade forces and modal pressure projections remain
    /// in their original owner and enter the work account exactly once.
    /// See `ApertureDg::step_with_finite_gas_controlled` for the physical partition
    /// and callback contract. An independent upstream finite tract is refused.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_finite_gas_controlled<S, E, C: FnMut() -> bool>(
        &mut self,
        gas: &mut ClosedGasParcel,
        drive: TubeDrive,
        forces_n: [f64; 2],
        pressures_pa: [f64; 2],
        controls: FiniteGasDgSolve,
        evolve: impl FnMut(&ClosedGasParcel, f64, f64, &mut C) -> Result<FiniteGasDgBackground<S>, E>,
        cancelled: &mut C,
    ) -> Result<FiniteGasApertureDgStep<S>, FiniteGasDgError<E>> {
        admit(self, gas, controls, cancelled)?;
        let initial_energy_j = self.stored_energy_j();
        let relation = self.load.port_response().map_err(acoustic)?;
        let mechanical = self
            .aperture
            .preview(
                ApertureDrive {
                    upstream_pressure_pa: drive.upstream_pressure_pa,
                    incoming_pressure_pa: relation.history_pressure_pa / 2.0,
                    body_flow_m3_s: drive.body_flow_m3_s,
                },
                forces_n,
                pressures_pa,
            )
            .map_err(acoustic)?;
        let mut marched = self.load.clone();
        let frame = accept_acoustic_candidate(
            &mut marched,
            relation,
            mechanical.frame(),
            mechanical.spatial_pressure(),
            drive,
        )
        .map_err(acoustic)?;
        let (load, next_gas, report) = finish(
            &self.load,
            &marched,
            gas,
            frame,
            initial_energy_j,
            controls,
            evolve,
            cancelled,
        )?;
        self.aperture
            .rebind_exact_gas_load(load.duct().inlet_gas(), load.impedance_pa_s_m3())
            .map_err(acoustic)?;
        self.aperture.accept(mechanical);
        self.load = load;
        *gas = next_gas;
        Ok(report)
    }

    /// Sourced fixed-volume heating with original independent mechanical drives.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_finite_gas_heat_controlled<C: FnMut() -> bool>(
        &mut self,
        gas: &mut ClosedGasParcel,
        drive: TubeDrive,
        forces_n: [f64; 2],
        pressures_pa: [f64; 2],
        heat_j: f64,
        gas_solve: GasEnergySolve,
        controls: FiniteGasDgSolve,
        cancelled: &mut C,
    ) -> Result<FiniteGasApertureDgStep<GasCompositionStep>, FiniteGasDgError<GasEnergyError>> {
        self.step_with_finite_gas_controlled(
            gas,
            drive,
            forces_n,
            pressures_pa,
            controls,
            |old, _, work, cancel| heated_finite_gas(old, heat_j, work, gas_solve, cancel),
            cancelled,
        )
    }
}
