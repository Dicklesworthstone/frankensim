//! Distributed wall/vortex heating during finite heater-to-sink transport.
//!
//! This composes the original paired thermal contact with an explicit spatial
//! caloric destination for every channel cell. Both endpoint exchanges and the
//! channel's actual dissipated joules publish on one native thermal timeline.
//! The destinations must differ from the endpoint cells: depositing additional
//! heat into an endpoint after solving its contact would invalidate that
//! contact's midpoint temperature. No geometry or heat route is inferred.
//!
//! Native sources and conduction march once, then the contact and distributed
//! dissipation are applied. This remains a first-order thermal-network split.
//! Deposited heat does not relinearize the channel's frozen helium properties.
//! Reference counterflow heat, wave availability and irreversible heat retain
//! their separate original accounts; no nonlinear total-energy theorem is added.

use std::collections::{BTreeMap, BTreeSet};

use fs_duct::superfluid::{MutualFriction, TwoFluidChannel};
use fs_thermal::enthalpy::{CellId, DepositionSimulation, EnergyAudit};

use super::super::counterflow::{
    CounterflowHeatControls, CounterflowHeatError, CounterflowThermalContact, check, thermal_check,
};
use super::{CounterflowPairHeatStep, advance_pair_prepared};

/// Endpoint transport and separately routed irreversible heat, on one interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CounterflowDissipationStep {
    /// Original paired-contact reports. Its caloric intake contains endpoint
    /// heat only, so all reference heat and entropy identities retain meaning.
    pub contacts: CounterflowPairHeatStep,
    /// Original caloric intake of distributed normal-wall and mutual losses.
    pub dissipation_intake: EnergyAudit,
    pub dissipated_heat_j: f64,
    /// Actual distributed caloric intake minus the channel's dissipation [J].
    pub dissipation_transfer_residual_j: f64,
}

/// Advance both finite contacts and deposit each channel cell's loss once.
///
/// `destinations[i]` receives cell i's accepted normal-viscous plus mutual-
/// friction heat. Destinations must be active, unique, and distinct from both
/// endpoint contact cells. Their native caloric laws and any phase plateaux
/// remain authoritative. The original thermal network advances exactly once;
/// there is no second wave solve or damping pass to estimate the heat.
///
/// Every proposal remains private until the physical channel and both caloric
/// intakes pass. A material refusal or cancellation preserves all wave/profile,
/// thermal and heat-history state. Retry never replays previously accepted heat.
#[allow(clippy::too_many_arguments)]
pub fn advance_counterflow_heat_pair_with_dissipation(
    channel: &mut TwoFluidChannel,
    thermal: &mut DepositionSimulation<'_>,
    contacts: [CounterflowThermalContact; 2],
    destinations: &[CellId],
    dt_s: f64,
    friction: MutualFriction,
    controls: CounterflowHeatControls,
    cancelled: impl Fn() -> bool,
) -> Result<CounterflowDissipationStep, CounterflowHeatError> {
    check(&cancelled)?;
    if destinations.len() != channel.cells() {
        return Err(CounterflowHeatError::Invalid(
            "distributed dissipation requires one destination per channel cell",
        ));
    }
    let mut seen = BTreeSet::new();
    for (i, &id) in destinations.iter().enumerate() {
        if i % 128 == 0 {
            check(&cancelled)?;
        }
        if contacts.iter().any(|contact| contact.cell == id)
            || !seen.insert(id)
            || !thermal.cells().contains_key(&id)
        {
            return Err(CounterflowHeatError::Invalid(
                "dissipation destinations must be active, unique and disjoint from both contacts",
            ));
        }
    }
    advance_pair_prepared(
        channel,
        thermal,
        contacts,
        dt_s,
        friction,
        controls,
        &cancelled,
        |candidate, contacts, per_cell_heat| {
            if per_cell_heat.len() != destinations.len() {
                return Err(CounterflowHeatError::Invalid(
                    "dissipation field and destination mapping differ",
                ));
            }
            let mut heat = BTreeMap::new();
            for (i, (&id, &joules)) in destinations.iter().zip(per_cell_heat).enumerate() {
                if i % 128 == 0 {
                    check(&cancelled)?;
                }
                if !joules.is_finite() || joules < 0.0 {
                    return Err(CounterflowHeatError::Invalid(
                        "distributed dissipation must be finite nonnegative heat",
                    ));
                }
                heat.insert(id, joules);
            }
            let (next, intake) = candidate.apply_heat_j(
                &heat,
                controls.cell_energy_tolerance_j,
                controls.transfer_energy_tolerance_j,
                || thermal_check(&cancelled),
            )?;
            let dissipated = contacts.channel.dissipated_heat_j();
            let residual = intake.external_heat_j - dissipated;
            if !residual.is_finite() || residual.abs() > controls.transfer_energy_tolerance_j {
                return Err(CounterflowHeatError::Invalid(
                    "distributed caloric intake does not match channel dissipation",
                ));
            }
            check(&cancelled)?;
            Ok((
                next,
                CounterflowDissipationStep {
                    contacts,
                    dissipation_intake: intake,
                    dissipated_heat_j: dissipated,
                    dissipation_transfer_residual_j: residual,
                },
            ))
        },
    )
}
