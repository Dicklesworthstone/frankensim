//! Finite-cell heat from both adhesive rolling hysteresis and junction release.
//!
//! The original source regenerates each increment. Only its irreversible
//! rolling loss uses the independently declared rolling partition; normal snap
//! and released spring energy retain the existing release partition. All three
//! channels enter one caloric impulse, so a finite cell admits their sum before
//! any history is published. Signed normal/stiffness work and recoverable
//! rolling storage remain mechanical quantities.
//!
//! As with the snap owner, this is a frozen-source, same-time impulse split.
//! The partition is supplied explicitly; it is not a solved contact conduction
//! law. Temperature feedback and transport between physical events are separate.

use super::{
    JkrRollingSnapHeat, JkrRollingSnapHeatCandidate, JkrRollingSnapHeatError,
    JkrRollingSnapHeatReport, JkrRollingSnapHeatStamp, JkrRollingSnapSegment, JkrSnapHeatControls,
    JkrSnapHeatError,
};
use fs_contact::normal_patch::adhesive::rolling::events::{
    JkrRollingEventReport, JkrRollingEventState, SourceJkrRollingEvents,
};
use fs_math::dd::Dd;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation, EnergyAudit};
use fs_tribo::{DissipationStep, HeatPartition, WorkLedger};
use std::collections::BTreeMap;

/// Atomic native-body motion and all three source-owned heat channels.
pub mod native;

/// Original snap receiver binding plus an explicit destination for plastic
/// rolling loss. The new ledger begins at this binding; previously accumulated
/// source loss is retained and is never heated a second time.
#[derive(Clone)]
pub struct JkrRollingDissipationHeat<'s> {
    snap: JkrRollingSnapHeat<'s>,
    rolling_partition: HeatPartition,
    rolling_outside_destination: Option<String>,
    rolling_loss_ledger: WorkLedger,
}

/// Distinct source channels and their single actual caloric deposition.
#[derive(Clone, Debug, PartialEq)]
pub struct JkrRollingDissipationHeatReport {
    pub initial_state: JkrRollingEventState,
    pub final_state: JkrRollingEventState,
    pub segments: Vec<JkrRollingEventReport>,
    pub normal_release: WorkLedger,
    pub rolling_release: WorkLedger,
    pub rolling_loss: WorkLedger,
    /// Sum of the three channels, excluding all recoverable mechanical work.
    pub work: WorkLedger,
    /// Destination for the outside shares of normal and spring release.
    pub release_outside_destination: Option<String>,
    /// Independent destination for the outside share of plastic rolling loss.
    pub rolling_outside_destination: Option<String>,
    pub thermal: EnergyAudit,
    /// Actual finite-cell energy change plus outside work minus all three
    /// source losses. This is a transfer check, not a total-system certificate.
    pub transfer_residual_j: f64,
}

pub struct JkrRollingDissipationHeatCandidate<'s, 'm> {
    history: JkrRollingDissipationHeat<'s>,
    thermal: DepositionSimulation<'m>,
    report: JkrRollingDissipationHeatReport,
}
impl<'s, 'm> JkrRollingDissipationHeatCandidate<'s, 'm> {
    pub const fn history(&self) -> &JkrRollingDissipationHeat<'s> {
        &self.history
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub const fn report(&self) -> &JkrRollingDissipationHeatReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        JkrRollingDissipationHeat<'s>,
        DepositionSimulation<'m>,
        JkrRollingDissipationHeatReport,
    ) {
        (self.history, self.thermal, self.report)
    }
}

impl<'s> JkrRollingDissipationHeat<'s> {
    /// Preserve the original source, exact caloric objects, current spring and
    /// snap ledgers. No default conversion fraction or destination is assumed.
    pub fn new(
        snap: JkrRollingSnapHeat<'s>,
        rolling_partition: HeatPartition,
        rolling_outside_destination: Option<String>,
    ) -> Result<Self, JkrRollingSnapHeatError> {
        if rolling_outside_destination
            .as_ref()
            .is_some_and(|name| name.trim().is_empty())
            || (rolling_partition.other() > 0.0 && rolling_outside_destination.is_none())
        {
            return Err(JkrSnapHeatError::Invalid(
                "rolling loss needs a named outside destination",
            )
            .into());
        }
        DissipationStep::from_work(0.0, rolling_partition)?;
        Ok(Self {
            snap,
            rolling_partition,
            rolling_outside_destination,
            rolling_loss_ledger: WorkLedger::default(),
        })
    }

    pub fn source(&self) -> &SourceJkrRollingEvents {
        self.snap.source()
    }
    pub const fn state(&self) -> &JkrRollingEventState {
        self.snap.state()
    }
    pub fn stamp(&self) -> JkrRollingSnapHeatStamp {
        self.snap.stamp()
    }
    pub const fn normal_release_ledger(&self) -> &WorkLedger {
        self.snap.normal_release_ledger()
    }
    pub const fn rolling_release_ledger(&self) -> &WorkLedger {
        self.snap.rolling_release_ledger()
    }
    pub const fn rolling_loss_ledger(&self) -> &WorkLedger {
        &self.rolling_loss_ledger
    }

    /// Advance original source segments, then apply all declared heat in one
    /// finite-cell proposal. Both release and plastic loss are derived from
    /// sealed source history; authored heat totals or edited reports are inputs
    /// to neither this method nor its native-body composition.
    pub fn propose<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        expected: JkrRollingSnapHeatStamp,
        segments: &[JkrRollingSnapSegment],
        controls: JkrSnapHeatControls,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<JkrRollingDissipationHeatCandidate<'s, 'm>, JkrRollingSnapHeatError> {
        let staged = stage(
            &self.snap,
            thermal,
            expected,
            segments,
            controls,
            Some(RollingConversion {
                partition: self.rolling_partition,
                outside_destination: &self.rolling_outside_destination,
                ledger: self.rolling_loss_ledger,
            }),
            checkpoint,
        )?;
        Ok(JkrRollingDissipationHeatCandidate {
            history: Self {
                snap: staged.history,
                rolling_partition: self.rolling_partition,
                rolling_outside_destination: self.rolling_outside_destination.clone(),
                rolling_loss_ledger: staged.rolling_loss_ledger,
            },
            thermal: staged.thermal,
            report: staged.report,
        })
    }

    /// Publish the complete source, all channel ledgers and finite calorics.
    pub fn apply<'m>(
        &mut self,
        thermal: &mut DepositionSimulation<'m>,
        expected: JkrRollingSnapHeatStamp,
        segments: &[JkrRollingSnapSegment],
        controls: JkrSnapHeatControls,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<JkrRollingDissipationHeatReport, JkrRollingSnapHeatError> {
        let (history, next, report) = self
            .propose(thermal, expected, segments, controls, checkpoint)?
            .into_parts();
        *self = history;
        *thermal = next;
        Ok(report)
    }
}

struct RollingConversion<'a> {
    partition: HeatPartition,
    outside_destination: &'a Option<String>,
    ledger: WorkLedger,
}

struct StagedHeat<'s, 'm> {
    history: JkrRollingSnapHeat<'s>,
    thermal: DepositionSimulation<'m>,
    rolling_loss_ledger: WorkLedger,
    report: JkrRollingDissipationHeatReport,
}

/// The legacy snap-only entry point deliberately supplies no rolling conversion.
/// Its work ledger and caloric audit therefore retain their original semantics.
pub(super) fn propose_snap<'s, 'm>(
    owner: &JkrRollingSnapHeat<'s>,
    thermal: &DepositionSimulation<'m>,
    expected: JkrRollingSnapHeatStamp,
    segments: &[JkrRollingSnapSegment],
    controls: JkrSnapHeatControls,
    checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<JkrRollingSnapHeatCandidate<'s, 'm>, JkrRollingSnapHeatError> {
    let staged = stage(
        owner, thermal, expected, segments, controls, None, checkpoint,
    )?;
    let r = staged.report;
    Ok(JkrRollingSnapHeatCandidate {
        history: staged.history,
        thermal: staged.thermal,
        report: JkrRollingSnapHeatReport {
            initial_state: r.initial_state,
            final_state: r.final_state,
            segments: r.segments,
            normal_release: r.normal_release,
            rolling_release: r.rolling_release,
            work: r.work,
            outside_destination: r.release_outside_destination,
            thermal: r.thermal,
            transfer_residual_j: r.transfer_residual_j,
        },
    })
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn stage<'s, 'm>(
    owner: &JkrRollingSnapHeat<'s>,
    thermal: &DepositionSimulation<'m>,
    expected: JkrRollingSnapHeatStamp,
    segments: &[JkrRollingSnapSegment],
    controls: JkrSnapHeatControls,
    conversion: Option<RollingConversion<'_>>,
    mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<StagedHeat<'s, 'm>, JkrRollingSnapHeatError> {
    checkpoint()?;
    if expected != owner.stamp {
        return Err(JkrSnapHeatError::HistoryMismatch.into());
    }
    if segments.is_empty()
        || segments.len() > controls.maximum_segments
        || thermal.time_s() < owner.accepted_time_s
        || [
            controls.cell_energy_tolerance_j,
            controls.thermal_energy_tolerance_j,
            controls.transfer_energy_tolerance_j,
        ]
        .iter()
        .any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err(JkrSnapHeatError::Invalid("invalid segment, clock or heat budget").into());
    }
    let revision = owner
        .stamp
        .revision
        .checked_add(1)
        .ok_or(JkrSnapHeatError::Invalid("history version exhausted"))?;
    owner.admit_receivers(thermal)?;
    let mut state = owner.stamp.state.clone();
    let mut reports = Vec::with_capacity(segments.len());
    let mut normal_release = WorkLedger::default();
    let mut rolling_release = WorkLedger::default();
    let mut rolling_loss = WorkLedger::default();
    let mut work = WorkLedger::default();
    // Retain weak rolling loss beside a much larger normal snap. The rounded
    // work summary is not the source of truth for transfer admission.
    let mut source_energy = Dd::ZERO;
    let mut normal_release_ledger = owner.normal_release_ledger;
    let mut rolling_release_ledger = owner.rolling_release_ledger;
    let mut rolling_loss_ledger = conversion
        .as_ref()
        .map_or(WorkLedger::default(), |c| c.ledger);
    for segment in segments {
        checkpoint()?;
        let step = match *segment {
            JkrRollingSnapSegment::NormalIndentation(delta) => {
                owner.source.advance_normal(&state, delta)?
            }
            JkrRollingSnapSegment::RollingAngleIncrement(angle) => {
                owner.source.advance_rolling(&state, angle)?
            }
        };
        let (next, report) = step.into_parts();
        let normal = DissipationStep::from_work(report.normal_released_energy_j, owner.partition)?;
        let release =
            DissipationStep::from_work(report.rolling_released_energy_j, owner.partition)?;
        source_energy = source_energy
            + Dd::from_f64(report.normal_released_energy_j)
            + Dd::from_f64(report.rolling_released_energy_j);
        normal_release.record(normal)?;
        rolling_release.record(release)?;
        normal_release_ledger.record(normal)?;
        rolling_release_ledger.record(release)?;
        work.record(normal)?;
        work.record(release)?;
        if let Some(c) = &conversion {
            let loss = DissipationStep::from_work(report.irreversible_rolling_loss_j, c.partition)?;
            source_energy = source_energy + Dd::from_f64(report.irreversible_rolling_loss_j);
            rolling_loss.record(loss)?;
            rolling_loss_ledger.record(loss)?;
            work.record(loss)?;
        }
        reports.push(report);
        state = next;
    }
    let heat = BTreeMap::from([
        (owner.receivers[0].cell, work.surface_a_heat_j()),
        (owner.receivers[1].cell, work.surface_b_heat_j()),
    ]);
    let (thermal, audit) = thermal.apply_heat_j(
        &heat,
        controls.cell_energy_tolerance_j,
        controls.thermal_energy_tolerance_j,
        &mut checkpoint,
    )?;
    let residual = (Dd::from_f64(work.surface_a_heat_j())
        + Dd::from_f64(work.surface_b_heat_j())
        + Dd::from_f64(work.other_work_j())
        - source_energy
        + Dd::from_f64(audit.residual_j))
    .to_f64();
    if !residual.is_finite() || residual.abs() > controls.transfer_energy_tolerance_j {
        return Err(JkrSnapHeatError::TransferBalance(residual).into());
    }
    let report = JkrRollingDissipationHeatReport {
        initial_state: owner.stamp.state.clone(),
        final_state: state.clone(),
        segments: reports,
        normal_release,
        rolling_release,
        rolling_loss,
        work,
        release_outside_destination: owner.outside_destination.clone(),
        rolling_outside_destination: conversion
            .as_ref()
            .and_then(|c| c.outside_destination.clone()),
        thermal: audit,
        transfer_residual_j: residual,
    };
    let history = JkrRollingSnapHeat {
        source: owner.source.clone(),
        caloric: owner.caloric,
        receivers: owner.receivers,
        partition: owner.partition,
        outside_destination: owner.outside_destination.clone(),
        stamp: JkrRollingSnapHeatStamp {
            source: owner.stamp.source,
            revision,
            state,
        },
        accepted_time_s: thermal.time_s(),
        normal_release_ledger,
        rolling_release_ledger,
    };
    checkpoint()?;
    Ok(StagedHeat {
        history,
        thermal,
        rolling_loss_ledger,
        report,
    })
}
