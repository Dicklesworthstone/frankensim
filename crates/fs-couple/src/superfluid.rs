//! Closed-channel mutual-friction heat on the original finite thermal timeline.
//!
//! The channel prepares its wave endpoint and per-cell dissipated joules without
//! publication. The original thermal owner advances its authored boundaries and
//! contacts to that same endpoint, then accepts those joules through its own
//! caloric laws. Both histories publish only after every source, energy and
//! cancellation check succeeds. Retrying a refused step deposits no heat.
//!
//! This is a first-order endpoint heat split. It does not change the channel's
//! frozen thermodynamic reference or fold finite thermal-cell temperatures into
//! its linear temperature perturbation. The caller declares physical heat
//! destinations; geometry matching and interface transport are not inferred.
//! Closed ends exclude unhandled bulk mass/entropy/counterflow imports. Native
//! thermal storage conventions remain authoritative: the combined account is
//! represented wave plus thermal storage, not an unconditional total-internal-
//! energy theorem for a thermal law that stores enthalpy.

use std::collections::{BTreeMap, BTreeSet};

use fs_duct::superfluid::{
    BodyForce, BoundaryVelocities, ChannelError, ChannelStep, ForcedChannelStep, HeatStepError,
    MutualFriction, TwoFluidChannel,
};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionSimulation, EnergyAudit, StepConfig,
};

/// Reciprocal pressure impulse and motion on the original native rigid-body owner.
pub mod piston;

/// Two freely moving native bodies with reciprocal channel pressure and forcing.
pub mod piston_pair;

/// Original thermal stepping controls plus explicit heat and combined-account gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrictionHeatControls {
    /// Original conductive/source time integration and conservation admission.
    pub thermal: StepConfig,
    /// Absolute permitted error in each receiving cell's heat increment [J].
    pub cell_heat_tolerance_j: f64,
    /// Sum of absolute local heat-increment errors, and transfer mismatch [J].
    pub heat_energy_tolerance_j: f64,
    /// Whole interval represented wave-plus-thermal storage defect [J].
    pub combined_energy_tolerance_j: f64,
}

/// Both unchanged owner reports and the one internal energy transfer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrictionHeatReport {
    /// Original channel step, including the physical friction debit.
    pub channel: ChannelStep,
    /// Original thermal march, excluding the endpoint friction deposit.
    pub thermal_march: EnergyAudit,
    /// Original endpoint caloric update; its external heat is INTERNAL to this
    /// composed system and is not added to the thermal march's source account.
    pub heat_intake: EnergyAudit,
    /// Deposited joules minus the original channel friction debit [J].
    pub transfer_residual_j: f64,
    /// Sum of the original wave, march, heat-intake and transfer residuals [J].
    pub combined_storage_residual_j: f64,
}

/// Failed channel, mapping or caloric proposals preserve both accepted owners.
#[derive(Debug)]
pub enum FrictionHeatError {
    /// Incompatible timeline, incomplete/duplicate destination or invalid control.
    Invalid(&'static str),
    /// Original two-fluid wave/source refusal.
    Channel(ChannelError),
    /// Original caloric, thermal work or composed-account refusal.
    Thermal(DepositionError),
    /// Caller cancellation before joint publication.
    Cancelled,
}
impl core::fmt::Display for FrictionHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "superfluid friction heat refused: {self:?}")
    }
}
impl std::error::Error for FrictionHeatError {}
impl From<DepositionError> for FrictionHeatError {
    fn from(value: DepositionError) -> Self {
        match value {
            DepositionError::Cancelled => Self::Cancelled,
            other => Self::Thermal(other),
        }
    }
}
impl From<HeatStepError<DepositionError>> for FrictionHeatError {
    fn from(value: HeatStepError<DepositionError>) -> Self {
        match value {
            HeatStepError::Channel(ChannelError::Cancelled) => Self::Cancelled,
            HeatStepError::Channel(error) => Self::Channel(error),
            HeatStepError::Heat(error) => error.into(),
        }
    }
}

/// Advance a closed two-fluid channel and its original finite heat destinations.
///
/// `destinations[i]` receives channel cell i's actual dissipated joules exactly
/// once. Every destination must already be active and appear once; unrelated
/// thermal cells may participate in their original conduction network. Thermal
/// events and authored sources are advanced by the original `advance_to`, then
/// the original `apply_heat_j` applies only the new interval's heat.
///
/// The clocks must match exactly on entry. No fallible work follows channel
/// publication: assignment publishes the prepared thermal candidate. The gate
/// is shared immutably by both original owners (e.g. `|| gate.is_requested()`).
/// No second damping pass, cumulative-heat replay, fake source power, copied
/// heat-capacity law, or feedback into frozen wave coefficients is introduced.
#[allow(clippy::too_many_arguments)]
pub fn advance_closed_friction_heat(
    channel: &mut TwoFluidChannel,
    thermal: &mut DepositionSimulation<'_>,
    destinations: &[CellId],
    dt_s: f64,
    friction: MutualFriction,
    controls: FrictionHeatControls,
    cancelled: impl Fn() -> bool,
) -> Result<FrictionHeatReport, FrictionHeatError> {
    validate_heat_route(channel, thermal, destinations, controls, &cancelled)?;
    let (_, (next_thermal, report)) = channel.step_with_heat(
        dt_s,
        BoundaryVelocities::CLOSED,
        friction,
        || cancelled(),
        |wave, per_cell_heat| {
            prepare_heat(
                thermal,
                destinations,
                wave,
                per_cell_heat,
                controls,
                &cancelled,
            )
        },
    )?;
    // Native channel publication has succeeded; this assignment cannot refuse.
    *thermal = next_thermal;
    Ok(report)
}

/// Drive a closed two-fluid channel and deposit its actual friction heat once.
///
/// The original physical source owner advances the component momenta and
/// records external body-force work. The original thermal owner advances its
/// authored timeline and accepts only the resulting per-cell friction joules.
/// Force work is external to this composition; friction is an internal transfer.
/// The channel's original wave balance already debits force work, so it must
/// not be added again to the heat intake or combined residual.
///
/// Mapping, clocks, caloric laws and frozen-reference semantics are identical to
/// [`advance_closed_friction_heat`]. A thermal refusal or cancellation preserves
/// the accepted force response, wave fields, momentum, heat and thermal cells.
#[allow(clippy::too_many_arguments)]
pub fn advance_forced_closed_friction_heat(
    channel: &mut TwoFluidChannel,
    thermal: &mut DepositionSimulation<'_>,
    destinations: &[CellId],
    dt_s: f64,
    force: BodyForce<'_>,
    friction: MutualFriction,
    controls: FrictionHeatControls,
    cancelled: impl Fn() -> bool,
) -> Result<(ForcedChannelStep, FrictionHeatReport), FrictionHeatError> {
    validate_heat_route(channel, thermal, destinations, controls, &cancelled)?;
    let (forced, (next_thermal, report)) = channel.step_forced_closed_with_heat(
        dt_s,
        force,
        friction,
        || cancelled(),
        |forced, per_cell_heat| {
            prepare_heat(
                thermal,
                destinations,
                &forced.channel,
                per_cell_heat,
                controls,
                &cancelled,
            )
        },
    )?;
    *thermal = next_thermal;
    Ok((forced, report))
}

fn validate_heat_route(
    channel: &TwoFluidChannel,
    thermal: &DepositionSimulation<'_>,
    destinations: &[CellId],
    controls: FrictionHeatControls,
    cancelled: &impl Fn() -> bool,
) -> Result<(), FrictionHeatError> {
    check(cancelled)?;
    if channel.time_s() != thermal.time_s() {
        return Err(FrictionHeatError::Invalid(
            "channel and thermal clocks must match",
        ));
    }
    if destinations.len() != channel.cells() {
        return Err(FrictionHeatError::Invalid(
            "one active thermal destination per channel cell is required",
        ));
    }
    if [
        controls.cell_heat_tolerance_j,
        controls.heat_energy_tolerance_j,
        controls.combined_energy_tolerance_j,
    ]
    .into_iter()
    .any(|v| !v.is_finite() || v < 0.0)
    {
        return Err(FrictionHeatError::Invalid(
            "heat and combined-energy tolerances must be finite and nonnegative",
        ));
    }
    let mut distinct = BTreeSet::new();
    for (i, &id) in destinations.iter().enumerate() {
        if i % 128 == 0 {
            check(cancelled)?;
        }
        if !thermal.cells().contains_key(&id) || !distinct.insert(id) {
            return Err(FrictionHeatError::Invalid(
                "thermal destinations must be active and unique",
            ));
        }
    }
    Ok(())
}

fn prepare_heat<'a>(
    thermal: &DepositionSimulation<'a>,
    destinations: &[CellId],
    wave: &ChannelStep,
    per_cell_heat: &[f64],
    controls: FrictionHeatControls,
    cancelled: &impl Fn() -> bool,
) -> Result<(DepositionSimulation<'a>, FrictionHeatReport), DepositionError> {
    check(cancelled)?;
    if per_cell_heat.len() != destinations.len() {
        return Err(DepositionError::InvalidInput(
            "channel heat field and destination map differ",
        ));
    }
    let mut heat = BTreeMap::new();
    for (i, (&id, &joules)) in destinations.iter().zip(per_cell_heat).enumerate() {
        if i % 128 == 0 {
            check(cancelled)?;
        }
        if !joules.is_finite() || joules < 0.0 {
            return Err(DepositionError::InvalidInput(
                "channel returned invalid physical friction heat",
            ));
        }
        heat.insert(id, joules);
    }
    let (marched, thermal_march) =
        thermal.advance_to(wave.time_s, controls.thermal, || check(cancelled))?;
    let (next, heat_intake) = marched.apply_heat_j(
        &heat,
        controls.cell_heat_tolerance_j,
        controls.heat_energy_tolerance_j,
        || check(cancelled),
    )?;
    let transfer_residual_j = heat_intake.external_heat_j - wave.mutual_friction_heat_j;
    // Preserve the owners' cell-relative residuals. The wave balance already
    // accounts for any external force work; friction crosses this interface once.
    // Subtracting large absolute thermal energies can erase a small heat deposit.
    let combined_storage_residual_j = wave.wave_balance_defect_j
        + thermal_march.residual_j
        + heat_intake.residual_j
        + transfer_residual_j;
    if !transfer_residual_j.is_finite() || !combined_storage_residual_j.is_finite() {
        return Err(DepositionError::Arithmetic);
    }
    if transfer_residual_j.abs() > controls.heat_energy_tolerance_j {
        return Err(DepositionError::EnergyBalance(transfer_residual_j));
    }
    if combined_storage_residual_j.abs() > controls.combined_energy_tolerance_j {
        return Err(DepositionError::EnergyBalance(combined_storage_residual_j));
    }
    Ok((
        next,
        FrictionHeatReport {
            channel: *wave,
            thermal_march,
            heat_intake,
            transfer_residual_j,
            combined_storage_residual_j,
        },
    ))
}

fn check(cancelled: &impl Fn() -> bool) -> Result<(), DepositionError> {
    if cancelled() {
        Err(DepositionError::Cancelled)
    } else {
        Ok(())
    }
}
