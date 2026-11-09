//! Finite caloric destinations for the original calibrated torsion lifecycle.
//!
//! Normal snap, remaining twist-spring release and plastic torsion loss keep
//! distinct source and cumulative accounts. The caller declares separate
//! release and plastic partitions; only their caloric shares enter the two
//! original finite receivers. Retained elastic storage and signed normal/area
//! work stay mechanical. Prior source losses are never deposited a second time.
//!
//! Ordered segments share the current thermal time. This frozen-source impulse
//! split adds no conducting link, temperature feedback, mixed-mode instability
//! or body impulse. Intervening transport belongs to the original thermal owner.
//! [`timed`] carries plastic heat through that original physical-time transport.
//! [`temperature`] continues loaded twist at actual receiver-cell temperatures,
//! preserving its history and assigning signed parameter work explicitly.

pub mod timed;

/// Actual thermal-cell temperatures continue loaded torsional material response.
pub mod temperature;

use super::{JkrSnapHeat, JkrSnapHeatControls, JkrSnapHeatError, Receiver};
use fs_blake3::ContentHash;
use fs_contact::normal_patch::adhesive::torsion::{
    SourceJkrTorsionError,
    events::{JkrTorsionEventReport, JkrTorsionEventState, SourceJkrTorsionEvents},
};
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_math::dd::Dd;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation, EnergyAudit};
use fs_tribo::{DissipationStep, HeatPartition, TriboError, WorkLedger};
use std::{collections::BTreeMap, sync::Arc};

/// Ordered monotone source segments; simultaneous paths or hidden reversals
/// must be resolved by the caller rather than inferred from endpoints.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum JkrTorsionHeatSegment {
    NormalIndentation(f64),
    TorsionAngleIncrement(f64),
}

#[derive(Clone, Debug, PartialEq)]
pub enum JkrTorsionHeatError {
    Source(SourceJkrTorsionError),
    Heat(JkrSnapHeatError),
}
impl From<SourceJkrTorsionError> for JkrTorsionHeatError {
    fn from(value: SourceJkrTorsionError) -> Self {
        Self::Source(value)
    }
}
impl From<JkrSnapHeatError> for JkrTorsionHeatError {
    fn from(value: JkrSnapHeatError) -> Self {
        Self::Heat(value)
    }
}
impl From<DepositionError> for JkrTorsionHeatError {
    fn from(value: DepositionError) -> Self {
        Self::Heat(value.into())
    }
}
impl From<TriboError> for JkrTorsionHeatError {
    fn from(value: TriboError) -> Self {
        Self::Heat(value.into())
    }
}
impl core::fmt::Display for JkrTorsionHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "JKR torsion heat: {self:?}")
    }
}
impl std::error::Error for JkrTorsionHeatError {}

/// Sealed complete mechanical history and heat-binding revision prevent replay.
#[derive(Clone, Debug, PartialEq)]
pub struct JkrTorsionHeatStamp {
    source: ContentHash,
    revision: u64,
    state: JkrTorsionEventState,
}

#[derive(Clone)]
pub struct JkrTorsionHeat<'s> {
    source: Arc<SourceJkrTorsionEvents>,
    caloric: [&'s EquilibriumEnthalpyPhaseCurve; 2],
    receivers: [Receiver; 2],
    release_partition: HeatPartition,
    release_outside_destination: Option<String>,
    torsion_partition: HeatPartition,
    torsion_outside_destination: Option<String>,
    stamp: JkrTorsionHeatStamp,
    accepted_time_s: f64,
    normal_release_ledger: WorkLedger,
    torsion_release_ledger: WorkLedger,
    torsion_loss_ledger: WorkLedger,
}

/// Original source observations and actual finite-cell transfer. Writable
/// reports are never accepted as an input to the next mechanical/heat step.
#[derive(Clone, Debug, PartialEq)]
pub struct JkrTorsionHeatReport {
    pub initial_state: JkrTorsionEventState,
    pub final_state: JkrTorsionEventState,
    pub segments: Vec<JkrTorsionEventReport>,
    pub normal_release: WorkLedger,
    pub torsion_release: WorkLedger,
    pub torsion_loss: WorkLedger,
    /// All three outgoing channels; excludes stored energy and signed work.
    pub work: WorkLedger,
    pub release_outside_destination: Option<String>,
    pub torsion_outside_destination: Option<String>,
    pub thermal: EnergyAudit,
    /// Actual caloric change + outside work - original outgoing source energy.
    /// A transfer audit, not a complete mechanical/thermal system certificate.
    pub transfer_residual_j: f64,
    /// Partition, accumulation and caloric defects cannot cancel one another.
    pub absolute_transfer_residual_j: f64,
}

pub struct JkrTorsionHeatCandidate<'s, 'm> {
    history: JkrTorsionHeat<'s>,
    thermal: DepositionSimulation<'m>,
    report: JkrTorsionHeatReport,
}
impl<'s, 'm> JkrTorsionHeatCandidate<'s, 'm> {
    pub const fn history(&self) -> &JkrTorsionHeat<'s> {
        &self.history
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub const fn report(&self) -> &JkrTorsionHeatReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        JkrTorsionHeat<'s>,
        DepositionSimulation<'m>,
        JkrTorsionHeatReport,
    ) {
        (self.history, self.thermal, self.report)
    }
}

impl<'s> JkrTorsionHeat<'s> {
    /// Preserve an already loaded or separated original cycle and bind exact
    /// registered caloric laws. Heat ledgers start here; earlier mechanical
    /// losses remain in source history and cannot be heated retroactively.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        source: SourceJkrTorsionEvents,
        initial: JkrTorsionEventState,
        thermal: &DepositionSimulation<'_>,
        cells: [CellId; 2],
        caloric: [&'s EquilibriumEnthalpyPhaseCurve; 2],
        release_partition: HeatPartition,
        release_outside_destination: Option<String>,
        torsion_partition: HeatPartition,
        torsion_outside_destination: Option<String>,
    ) -> Result<Self, JkrTorsionHeatError> {
        // Zero twist admits the complete sealed history without snapping an
        // attached junction sitting exactly at its marginal normal fold.
        if source.advance_torsion(&initial, 0.0)?.state() != &initial {
            return Err(JkrSnapHeatError::HistoryMismatch.into());
        }
        if torsion_outside_destination
            .as_ref()
            .is_some_and(|name| name.trim().is_empty())
            || (torsion_partition.other() > 0.0 && torsion_outside_destination.is_none())
        {
            return Err(JkrSnapHeatError::Invalid(
                "torsion loss needs a named outside destination",
            )
            .into());
        }
        DissipationStep::from_work(0.0, torsion_partition)?;
        let receivers = JkrSnapHeat::new(
            source.normal_source(),
            initial.normal_response().state(),
            thermal,
            cells,
            caloric,
            release_partition,
            release_outside_destination.clone(),
        )?
        .receivers;
        let stamp = JkrTorsionHeatStamp {
            source: source.normal_source().identity(),
            revision: 0,
            state: initial,
        };
        Ok(Self {
            source: Arc::new(source),
            caloric,
            receivers,
            release_partition,
            release_outside_destination,
            torsion_partition,
            torsion_outside_destination,
            stamp,
            accepted_time_s: thermal.time_s(),
            normal_release_ledger: WorkLedger::default(),
            torsion_release_ledger: WorkLedger::default(),
            torsion_loss_ledger: WorkLedger::default(),
        })
    }

    pub fn source(&self) -> &SourceJkrTorsionEvents {
        &self.source
    }
    pub const fn state(&self) -> &JkrTorsionEventState {
        &self.stamp.state
    }
    pub fn stamp(&self) -> JkrTorsionHeatStamp {
        self.stamp.clone()
    }
    pub const fn normal_release_ledger(&self) -> &WorkLedger {
        &self.normal_release_ledger
    }
    pub const fn torsion_release_ledger(&self) -> &WorkLedger {
        &self.torsion_release_ledger
    }
    pub const fn torsion_loss_ledger(&self) -> &WorkLedger {
        &self.torsion_loss_ledger
    }

    /// Regenerate the ordered original normal/twist segments and deposit all
    /// declared heat together at the current thermal time. Source, receiver,
    /// caloric, budget and cancellation refusals leave every owner unchanged.
    #[allow(clippy::too_many_lines)]
    pub fn propose<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        expected: JkrTorsionHeatStamp,
        segments: &[JkrTorsionHeatSegment],
        controls: JkrSnapHeatControls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<JkrTorsionHeatCandidate<'s, 'm>, JkrTorsionHeatError> {
        checkpoint()?;
        if expected != self.stamp {
            return Err(JkrSnapHeatError::HistoryMismatch.into());
        }
        if segments.is_empty()
            || segments.len() > controls.maximum_segments
            || thermal.time_s() < self.accepted_time_s
            || [
                controls.cell_energy_tolerance_j,
                controls.thermal_energy_tolerance_j,
                controls.transfer_energy_tolerance_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(JkrSnapHeatError::Invalid("invalid segment, clock or heat budget").into());
        }
        let revision = self
            .stamp
            .revision
            .checked_add(1)
            .ok_or(JkrSnapHeatError::Invalid("history version exhausted"))?;
        self.admit_receivers(thermal)?;
        let mut state = self.state().clone();
        let mut reports = Vec::with_capacity(segments.len());
        let mut channels = [WorkLedger::default(); 3];
        let mut cumulative = [
            self.normal_release_ledger,
            self.torsion_release_ledger,
            self.torsion_loss_ledger,
        ];
        let mut work = WorkLedger::default();
        let mut source_energy = Dd::ZERO;
        let mut partition_defect = Dd::ZERO;
        let mut absolute_partition = Dd::ZERO;
        for segment in segments {
            checkpoint()?;
            let step = match *segment {
                JkrTorsionHeatSegment::NormalIndentation(delta) => {
                    self.source.advance_normal(&state, delta)?
                }
                JkrTorsionHeatSegment::TorsionAngleIncrement(angle) => {
                    self.source.advance_torsion(&state, angle)?
                }
            };
            let (next, report) = step.into_parts();
            let losses = [
                (report.normal_released_energy_j, self.release_partition),
                (report.torsion_released_energy_j, self.release_partition),
                (report.irreversible_torsion_loss_j, self.torsion_partition),
            ];
            for (channel, (energy, partition)) in losses.into_iter().enumerate() {
                let loss = DissipationStep::from_work(energy, partition)?;
                let defect = Dd::from_f64(loss.surface_a_heat_j())
                    + Dd::from_f64(loss.surface_b_heat_j())
                    + Dd::from_f64(loss.other_work_j())
                    - Dd::from_f64(energy);
                partition_defect = partition_defect + defect;
                absolute_partition = absolute_partition + defect.abs();
                source_energy = source_energy + Dd::from_f64(energy);
                channels[channel].record(loss)?;
                cumulative[channel].record(loss)?;
                work.record(loss)?;
            }
            state = next;
            reports.push(report);
        }
        let heat = BTreeMap::from([
            (self.receivers[0].cell, work.surface_a_heat_j()),
            (self.receivers[1].cell, work.surface_b_heat_j()),
        ]);
        let (thermal, audit) = thermal.apply_heat_j(
            &heat,
            controls.cell_energy_tolerance_j,
            controls.thermal_energy_tolerance_j,
            &mut checkpoint,
        )?;
        // The original relative caloric audit retains weak increments beside
        // large background enthalpy. Check partition and accumulation defects
        // independently so opposite errors cannot hide in the signed total.
        let partition = Dd::from_f64(work.surface_a_heat_j())
            + Dd::from_f64(work.surface_b_heat_j())
            + Dd::from_f64(work.other_work_j())
            - source_energy;
        let residual = (partition + Dd::from_f64(audit.residual_j)).to_f64();
        let absolute = (absolute_partition
            + (partition - partition_defect).abs()
            + Dd::from_f64(audit.residual_j.abs()))
        .to_f64();
        if !residual.is_finite()
            || !absolute.is_finite()
            || residual.abs() > controls.transfer_energy_tolerance_j
            || absolute > controls.transfer_energy_tolerance_j
        {
            return Err(JkrSnapHeatError::TransferBalance(residual.abs().max(absolute)).into());
        }
        let report = JkrTorsionHeatReport {
            initial_state: self.state().clone(),
            final_state: state.clone(),
            segments: reports,
            normal_release: channels[0],
            torsion_release: channels[1],
            torsion_loss: channels[2],
            work,
            release_outside_destination: self.release_outside_destination.clone(),
            torsion_outside_destination: self.torsion_outside_destination.clone(),
            thermal: audit,
            transfer_residual_j: residual,
            absolute_transfer_residual_j: absolute,
        };
        let mut history = self.clone();
        history.stamp.revision = revision;
        history.stamp.state = state;
        history.accepted_time_s = thermal.time_s();
        history.normal_release_ledger = cumulative[0];
        history.torsion_release_ledger = cumulative[1];
        history.torsion_loss_ledger = cumulative[2];
        checkpoint()?;
        Ok(JkrTorsionHeatCandidate {
            history,
            thermal,
            report,
        })
    }

    /// Publish the whole original cycle and caloric impulse after the final
    /// caller checkpoint. A late refusal cannot publish partial source history.
    pub fn apply<'m>(
        &mut self,
        thermal: &mut DepositionSimulation<'m>,
        expected: JkrTorsionHeatStamp,
        segments: &[JkrTorsionHeatSegment],
        controls: JkrSnapHeatControls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<JkrTorsionHeatReport, JkrTorsionHeatError> {
        let (history, next, report) = self
            .propose(thermal, expected, segments, controls, &mut checkpoint)?
            .into_parts();
        checkpoint()?;
        *self = history;
        *thermal = next;
        Ok(report)
    }

    fn admit_receivers(&self, thermal: &DepositionSimulation<'_>) -> Result<(), JkrSnapHeatError> {
        for (side, receiver) in self.receivers.iter().enumerate() {
            let cell = thermal
                .cells()
                .get(&receiver.cell)
                .ok_or(JkrSnapHeatError::Invalid(
                    "original heat receiver is no longer active",
                ))?;
            if cell.segment_id != receiver.segment
                || cell.body_id != receiver.body
                || cell.material_id != receiver.material
                || cell.mass_kg != receiver.mass_kg
                || cell.geometry != receiver.geometry
                || !thermal.uses_material_law(cell.material_id, self.caloric[side])
            {
                return Err(JkrSnapHeatError::Invalid(
                    "heat receiver changed its original material, body, mass or geometry",
                ));
            }
        }
        Ok(())
    }
}
