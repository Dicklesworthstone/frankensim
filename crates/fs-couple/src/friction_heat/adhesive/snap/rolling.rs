//! Retained adhesive rolling cycles and original finite-cell snap heating.
//!
//! The original source regenerates every normal or angular segment. Normal
//! snap release and the rolling spring released at detachment remain distinct
//! accounts and enter only explicitly declared caloric/outside destinations.
//! Signed normal/rolling parameter work and irreversible rolling loss never
//! become snap heat. Previously accumulated source losses survive recontact.
//!
//! All releases in a proposal enter at the CURRENT thermal time. This is a
//! frozen-source impulse split, with no transport between its ordered segments,
//! no temperature feedback and no inferred conducting contact link.
//! Call separately at physical event times when intervening transport matters.
//! The [`native`] composition stages actual body recoil and this same heat
//! transfer together, checking that both regenerate one original lifecycle.

use super::{JkrSnapHeat, JkrSnapHeatControls, JkrSnapHeatError, Receiver};
use fs_blake3::ContentHash;
use fs_contact::normal_patch::adhesive::rolling::{
    SourceJkrRollingError,
    events::{JkrRollingEventReport, JkrRollingEventState, SourceJkrRollingEvents},
};
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_math::dd::Dd;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation, EnergyAudit};
use fs_tribo::{DissipationStep, HeatPartition, TriboError, WorkLedger};
use std::{collections::BTreeMap, sync::Arc};

/// Atomic native-body recoil, original lifecycle and finite-cell snap transfer.
pub mod native;

/// One physical monotone segment. Order is part of the prescribed path;
/// simultaneous normal/rolling motion and hidden reversals are not inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum JkrRollingSnapSegment {
    /// Prescribed signed normal approach at held rolling angle [m].
    NormalIndentation(f64),
    /// Prescribed rolling increment at held normal approach [rad].
    /// Detached travel cannot preload a later contact.
    RollingAngleIncrement(f64),
}

#[derive(Clone, Debug, PartialEq)]
pub enum JkrRollingSnapHeatError {
    Source(SourceJkrRollingError),
    Heat(JkrSnapHeatError),
}
impl From<SourceJkrRollingError> for JkrRollingSnapHeatError {
    fn from(value: SourceJkrRollingError) -> Self {
        Self::Source(value)
    }
}
impl From<JkrSnapHeatError> for JkrRollingSnapHeatError {
    fn from(value: JkrSnapHeatError) -> Self {
        Self::Heat(value)
    }
}
impl From<DepositionError> for JkrRollingSnapHeatError {
    fn from(value: DepositionError) -> Self {
        Self::Heat(value.into())
    }
}
impl From<TriboError> for JkrRollingSnapHeatError {
    fn from(value: TriboError) -> Self {
        Self::Heat(value.into())
    }
}
impl core::fmt::Display for JkrRollingSnapHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "JKR rolling snap heat: {self:?}")
    }
}
impl std::error::Error for JkrRollingSnapHeatError {}

/// Complete sealed source history and heat-binding revision. A stale stamp
/// cannot redeposit release, including after a full return to the same approach.
#[derive(Clone, Debug, PartialEq)]
pub struct JkrRollingSnapHeatStamp {
    source: ContentHash,
    revision: u64,
    state: JkrRollingEventState,
}

/// Original retained cycle and pinned finite caloric receivers. Release ledgers
/// begin at this binding; pre-existing source releases are never deposited again.
#[derive(Clone)]
pub struct JkrRollingSnapHeat<'s> {
    source: Arc<SourceJkrRollingEvents>,
    caloric: [&'s EquilibriumEnthalpyPhaseCurve; 2],
    receivers: [Receiver; 2],
    partition: HeatPartition,
    outside_destination: Option<String>,
    stamp: JkrRollingSnapHeatStamp,
    accepted_time_s: f64,
    normal_release_ledger: WorkLedger,
    rolling_release_ledger: WorkLedger,
}

/// Source-generated mechanics and actual finite-cell transfer. The source
/// reports retain the canonical normal JkrJump and separate rolling release;
/// neither an edited report nor an authored energy is accepted as input.
#[derive(Clone, Debug, PartialEq)]
pub struct JkrRollingSnapHeatReport {
    pub initial_state: JkrRollingEventState,
    pub final_state: JkrRollingEventState,
    pub segments: Vec<JkrRollingEventReport>,
    pub normal_release: WorkLedger,
    pub rolling_release: WorkLedger,
    /// Only outgoing release, excluding plastic rolling loss and signed work.
    pub work: WorkLedger,
    pub outside_destination: Option<String>,
    pub thermal: EnergyAudit,
    /// Actual finite-cell storage change + outside work - outgoing release.
    /// This audits transfer, not a complete dynamical system.
    pub transfer_residual_j: f64,
}

/// Cycle state, its release ledgers and caloric cells publish together.
pub struct JkrRollingSnapHeatCandidate<'s, 'm> {
    history: JkrRollingSnapHeat<'s>,
    thermal: DepositionSimulation<'m>,
    report: JkrRollingSnapHeatReport,
}
impl<'s, 'm> JkrRollingSnapHeatCandidate<'s, 'm> {
    pub const fn history(&self) -> &JkrRollingSnapHeat<'s> {
        &self.history
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub const fn report(&self) -> &JkrRollingSnapHeatReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        JkrRollingSnapHeat<'s>,
        DepositionSimulation<'m>,
        JkrRollingSnapHeatReport,
    ) {
        (self.history, self.thermal, self.report)
    }
}

impl<'s> JkrRollingSnapHeat<'s> {
    /// Bind an already loaded or detached original cycle without resetting its
    /// spring or prior loss. Caloric objects must be the exact registered laws
    /// of the normal source's ordered material cards. The declared partition
    /// applies to both release channels; outside work needs a physical name.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        source: SourceJkrRollingEvents,
        initial: JkrRollingEventState,
        thermal: &DepositionSimulation<'_>,
        cells: [CellId; 2],
        caloric: [&'s EquilibriumEnthalpyPhaseCurve; 2],
        partition: HeatPartition,
        outside_destination: Option<String>,
    ) -> Result<Self, JkrRollingSnapHeatError> {
        // A zero angular increment validates the sealed history without moving
        // a marginal attached normal state through its detachment event.
        if source.advance_rolling(&initial, 0.0)?.state() != &initial {
            return Err(JkrSnapHeatError::HistoryMismatch.into());
        }
        let receivers = {
            let admitted = JkrSnapHeat::new(
                source.normal_source(),
                initial.normal_response().state(),
                thermal,
                cells,
                caloric,
                partition,
                outside_destination.clone(),
            )?;
            admitted.receivers
        };
        let stamp = JkrRollingSnapHeatStamp {
            source: source.normal_source().identity(),
            revision: 0,
            state: initial,
        };
        Ok(Self {
            source: Arc::new(source),
            caloric,
            receivers,
            partition,
            outside_destination,
            stamp,
            accepted_time_s: thermal.time_s(),
            normal_release_ledger: WorkLedger::default(),
            rolling_release_ledger: WorkLedger::default(),
        })
    }

    pub fn source(&self) -> &SourceJkrRollingEvents {
        &self.source
    }
    pub const fn state(&self) -> &JkrRollingEventState {
        &self.stamp.state
    }
    pub fn stamp(&self) -> JkrRollingSnapHeatStamp {
        self.stamp.clone()
    }
    pub const fn normal_release_ledger(&self) -> &WorkLedger {
        &self.normal_release_ledger
    }
    pub const fn rolling_release_ledger(&self) -> &WorkLedger {
        &self.rolling_release_ledger
    }

    /// Regenerate ordered physical segments, then stage one same-time native
    /// heat impulse. `maximum_segments` counts both normal and angular segments.
    /// Rolling plastic loss stays with its explicit source account, unconverted;
    /// only outgoing snap release uses this binding's declared partition.
    /// Every source, receiver, caloric, budget or cancellation refusal is atomic.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn propose<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        expected: JkrRollingSnapHeatStamp,
        segments: &[JkrRollingSnapSegment],
        controls: JkrSnapHeatControls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<JkrRollingSnapHeatCandidate<'s, 'm>, JkrRollingSnapHeatError> {
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
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(JkrSnapHeatError::Invalid("invalid segment, clock or heat budget").into());
        }
        let revision = self
            .stamp
            .revision
            .checked_add(1)
            .ok_or(JkrSnapHeatError::Invalid("history version exhausted"))?;
        self.admit_receivers(thermal)?;
        let mut state = self.stamp.state.clone();
        let mut reports = Vec::with_capacity(segments.len());
        let mut normal_release = WorkLedger::default();
        let mut rolling_release = WorkLedger::default();
        let mut work = WorkLedger::default();
        let mut normal_release_ledger = self.normal_release_ledger;
        let mut rolling_release_ledger = self.rolling_release_ledger;
        for segment in segments {
            checkpoint()?;
            let step = match *segment {
                JkrRollingSnapSegment::NormalIndentation(delta) => {
                    self.source.advance_normal(&state, delta)?
                }
                JkrRollingSnapSegment::RollingAngleIncrement(angle) => {
                    self.source.advance_rolling(&state, angle)?
                }
            };
            let (next, report) = step.into_parts();
            let normal =
                DissipationStep::from_work(report.normal_released_energy_j, self.partition)?;
            let rolling =
                DissipationStep::from_work(report.rolling_released_energy_j, self.partition)?;
            normal_release.record(normal)?;
            rolling_release.record(rolling)?;
            normal_release_ledger.record(normal)?;
            rolling_release_ledger.record(rolling)?;
            work.record(normal)?;
            work.record(rolling)?;
            reports.push(report);
            state = next;
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
        let residual = (Dd::from_f64(work.surface_a_heat_j())
            + Dd::from_f64(work.surface_b_heat_j())
            + Dd::from_f64(work.other_work_j())
            - Dd::from_f64(work.dissipated_work_j())
            + Dd::from_f64(audit.residual_j))
        .to_f64();
        if !residual.is_finite() || residual.abs() > controls.transfer_energy_tolerance_j {
            return Err(JkrSnapHeatError::TransferBalance(residual).into());
        }
        let report = JkrRollingSnapHeatReport {
            initial_state: self.stamp.state.clone(),
            final_state: state.clone(),
            segments: reports,
            normal_release,
            rolling_release,
            work,
            outside_destination: self.outside_destination.clone(),
            thermal: audit,
            transfer_residual_j: residual,
        };
        let history = Self {
            source: self.source.clone(),
            caloric: self.caloric,
            receivers: self.receivers,
            partition: self.partition,
            outside_destination: self.outside_destination.clone(),
            stamp: JkrRollingSnapHeatStamp {
                source: self.stamp.source,
                revision,
                state,
            },
            accepted_time_s: thermal.time_s(),
            normal_release_ledger,
            rolling_release_ledger,
        };
        checkpoint()?;
        Ok(JkrRollingSnapHeatCandidate {
            history,
            thermal,
            report,
        })
    }

    /// Publish the source cycle, separate release ledgers and original finite
    /// caloric cells once, after the complete proposal has been admitted.
    #[allow(clippy::too_many_arguments)]
    pub fn apply<'m>(
        &mut self,
        thermal: &mut DepositionSimulation<'m>,
        expected: JkrRollingSnapHeatStamp,
        segments: &[JkrRollingSnapSegment],
        controls: JkrSnapHeatControls,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<JkrRollingSnapHeatReport, JkrRollingSnapHeatError> {
        let (history, next, report) = self
            .propose(thermal, expected, segments, controls, checkpoint)?
            .into_parts();
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
