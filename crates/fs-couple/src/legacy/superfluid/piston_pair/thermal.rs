//! Atomic fluid, two-body and original finite caloric-cell composition.

use fs_duct::superfluid::{BodyForce, MutualFriction, TwoFluidChannel};
use fs_mbd::{MassProperties, RigidBodyState};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};

use super::super::{
    FrictionHeatControls, FrictionHeatError, FrictionHeatReport, prepare_heat, validate_heat_route,
};
use super::{
    BalancedPistonBoundary, PistonControls, PistonError, PistonPairStep, advance_pair_with, sum,
};

/// Every original owner report and the whole represented-storage account.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PistonPairHeatReport {
    /// Original reciprocal fluid and native mechanical updates.
    pub pair: PistonPairStep,
    /// Original thermal march, per-cell intake and friction-transfer account.
    pub heat: FrictionHeatReport,
    /// Fluid, both bodies, endpoint fluid and thermal storage residual [J].
    pub combined_storage_residual_j: f64,
}

/// A refused proposal preserves fluid, both bodies and all thermal cells.
#[derive(Debug)]
pub enum PistonPairHeatError {
    /// Original boundary, wave, momentum or native mechanical refusal.
    Piston(PistonError),
    /// Original mapping, timeline, caloric or combined-storage refusal.
    Heat(FrictionHeatError),
    /// Cancellation before all original owners publish.
    Cancelled,
}

impl core::fmt::Display for PistonPairHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "superfluid piston-pair heat step refused: {self:?}")
    }
}
impl std::error::Error for PistonPairHeatError {}
impl From<PistonError> for PistonPairHeatError {
    fn from(error: PistonError) -> Self {
        match error {
            PistonError::Cancelled => Self::Cancelled,
            other => Self::Piston(other),
        }
    }
}
impl From<FrictionHeatError> for PistonPairHeatError {
    fn from(error: FrictionHeatError) -> Self {
        match error {
            FrictionHeatError::Cancelled => Self::Cancelled,
            other => Self::Heat(other),
        }
    }
}
impl From<DepositionError> for PistonPairHeatError {
    fn from(error: DepositionError) -> Self {
        FrictionHeatError::from(error).into()
    }
}

/// Advance the driven two-fluid channel, two native pistons and original heat cells.
///
/// Each `destinations[i]` is an active, unique original caloric cell receiving
/// channel cell i's current-interval friction heat. The thermal owner advances
/// its own authored sources, events and conductive contacts to the channel's
/// proposed endpoint before accepting those joules. Fluid and thermal clocks
/// must agree on entry. The native mechanical candidate and original caloric
/// candidate remain private until the channel's final acceptance succeeds.
///
/// Friction is an internal transfer in this composition. Fluid/body forces and
/// the original thermal march's sources remain external. The combined gate
/// adds the pair's mechanical-fluid residual to the original thermal march,
/// intake and transfer residuals, without counting the wave residual twice.
/// `heat_controls.combined_energy_tolerance_j` gates this full account.
///
/// This retains the original first-order endpoint heat split and the original
/// thermal storage basis. It does not feed finite caloric temperatures back into
/// frozen wave coefficients or claim a nonlinear total-internal-energy theorem
/// for a thermal law that stores enthalpy.
#[allow(clippy::too_many_arguments)]
pub fn advance_balanced_piston_pair_with_heat(
    channel: &mut TwoFluidChannel,
    bodies: &mut [RigidBodyState; 2],
    properties: [MassProperties; 2],
    boundaries: [BalancedPistonBoundary; 2],
    thermal: &mut DepositionSimulation<'_>,
    destinations: &[CellId],
    dt_s: f64,
    force: BodyForce<'_>,
    external_force_outward_n: [f64; 2],
    friction: MutualFriction,
    controls: PistonControls,
    heat_controls: FrictionHeatControls,
    cancelled: impl Fn() -> bool,
) -> Result<PistonPairHeatReport, PistonPairHeatError> {
    validate_heat_route(channel, thermal, destinations, heat_controls, &cancelled)?;
    let (_, (candidate, report)) = advance_pair_with(
        channel,
        bodies,
        properties,
        boundaries,
        dt_s,
        force,
        external_force_outward_n,
        friction,
        controls,
        &cancelled,
        |pair, per_cell_heat| {
            let (candidate, heat) = prepare_heat(
                thermal,
                destinations,
                &pair.channel,
                per_cell_heat,
                heat_controls,
                &cancelled,
            )?;
            let residual = sum([
                pair.combined_energy_residual_j,
                heat.thermal_march.residual_j,
                heat.heat_intake.residual_j,
                heat.transfer_residual_j,
            ]);
            if !residual.is_finite() {
                return Err(PistonPairHeatError::from(DepositionError::Arithmetic));
            }
            if residual.abs() > heat_controls.combined_energy_tolerance_j {
                return Err(PistonPairHeatError::from(DepositionError::EnergyBalance(
                    residual,
                )));
            }
            Ok((
                candidate,
                PistonPairHeatReport {
                    pair: *pair,
                    heat,
                    combined_storage_residual_j: residual,
                },
            ))
        },
    )?;
    // Fluid and both bodies were admitted together; this assignment cannot fail.
    *thermal = candidate;
    Ok(report)
}
