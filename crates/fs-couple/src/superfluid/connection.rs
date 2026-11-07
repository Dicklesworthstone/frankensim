//! Dissipative heat from two connected channels on one native thermal timeline.
//!
//! The connection prepares both wave endpoints and its physical junction
//! inertia. Both heat routes then use the original thermal owner and caloric
//! laws before any of the four accepted owners publishes. External thermal
//! sources advance once; each channel's new dissipative debit is deposited once.

use std::collections::BTreeSet;

use super::{
    FrictionHeatControls, FrictionHeatError, FrictionHeatReport, check, prepare_heat,
    validate_heat_route,
};
use fs_duct::superfluid::TwoFluidChannel;
use fs_duct::superfluid::connection::{ConnectionDrive, ConnectionStep, TwoFluidConnection};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};

/// Original connection and thermal reports for one jointly accepted interval.
#[derive(Debug)]
pub struct ConnectedFrictionHeatReport {
    /// Both channel reports, junction storage, and the original connection audit.
    pub connection: ConnectionStep,
    /// Native thermal march and heat intake for the first and second routes.
    /// Only the first report advances the thermal clock; the second march is
    /// at that same endpoint and therefore applies no additional source work.
    pub heat: [FrictionHeatReport; 2],
    /// Connection, native thermal, and internal heat-transfer residuals [J].
    pub combined_storage_residual_j: f64,
}

/// Advance two serially connected channels and deposit their actual losses.
///
/// Each route names one active, distinct thermal destination per channel cell.
/// Destinations must also be distinct across the two routes. The outer channel
/// ends are closed: bulk mass or heat exchange with an external reservoir
/// requires a separate physical transport composition. The internal connection
/// may be open or a superleak and retains its own inertial energy.
///
/// The first native thermal candidate advances authored sources and contacts to
/// the common endpoint. The second deposits its channel's joules on that same
/// candidate at the same time. Channel and junction work remain in their original
/// mechanical account and are never deposited again as heat. Frozen wave
/// coefficients and each native law's declared thermal storage basis are retained.
///
/// All four owners preserve their accepted state on caloric refusal, an energy
/// rejection, or cancellation. No fallible operation follows joint mechanical
/// publication; assignment publishes the already accepted thermal candidate.
#[allow(clippy::too_many_arguments)]
pub fn advance_connected_friction_heat(
    first: &mut TwoFluidChannel,
    second: &mut TwoFluidChannel,
    connection: &mut TwoFluidConnection,
    thermal: &mut DepositionSimulation<'_>,
    destinations: [&[CellId]; 2],
    dt_s: f64,
    drive: ConnectionDrive<'_>,
    controls: FrictionHeatControls,
    cancelled: impl Fn() -> bool,
) -> Result<ConnectedFrictionHeatReport, FrictionHeatError> {
    validate_heat_route(first, thermal, destinations[0], controls, &cancelled)?;
    validate_heat_route(second, thermal, destinations[1], controls, &cancelled)?;
    if [
        drive.outer.left.normal_m_s,
        drive.outer.left.superfluid_m_s,
        drive.outer.right.normal_m_s,
        drive.outer.right.superfluid_m_s,
    ]
    .into_iter()
    .any(|velocity| velocity != 0.0)
    {
        return Err(FrictionHeatError::Invalid(
            "connected friction heat requires closed outer channel ends",
        ));
    }
    let mut distinct = BTreeSet::new();
    for (i, &id) in destinations.into_iter().flatten().enumerate() {
        if i % 128 == 0 {
            check(&cancelled)?;
        }
        if !distinct.insert(id) {
            return Err(FrictionHeatError::Invalid(
                "connected channel heat destinations must be distinct across both routes",
            ));
        }
    }
    let (report, (next, heat, combined_storage_residual_j)) = connection.step_with_heat(
        first,
        second,
        dt_s,
        drive,
        || cancelled(),
        |report, first_heat, second_heat| {
            let (candidate, first_report) = prepare_heat(
                thermal,
                destinations[0],
                &report.channels[0],
                first_heat,
                controls,
                &cancelled,
            )?;
            let (next, second_report) = prepare_heat(
                &candidate,
                destinations[1],
                &report.channels[1],
                second_heat,
                controls,
                &cancelled,
            )?;
            // The connection residual already includes both original channel
            // wave residuals. Add only native thermal and transfer residuals.
            let mut combined = report.combined_energy_defect_j;
            for heat in [&first_report, &second_report] {
                combined += heat.thermal_march.residual_j
                    + heat.heat_intake.residual_j
                    + heat.transfer_residual_j;
            }
            if !combined.is_finite() {
                return Err(DepositionError::Arithmetic);
            }
            if combined.abs() > controls.combined_energy_tolerance_j {
                return Err(DepositionError::EnergyBalance(combined));
            }
            check(&cancelled)?;
            Ok((next, [first_report, second_report], combined))
        },
    )?;
    *thermal = next;
    Ok(ConnectedFrictionHeatReport {
        connection: report,
        heat,
        combined_storage_residual_j,
    })
}
