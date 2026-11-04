//! Original JKR snap work deposited once on the native thermal timeline.
//!
//! The original normal owner regenerates every accepted monotone segment from
//! retained, sealed branch history. Only its released jump energy is partitioned;
//! reversible normal work never becomes heat. The declared surface shares enter
//! the two original material caloric cells through `apply_heat_j`. Any remainder
//! has an explicit named destination outside those cells.
//!
//! The default proposals are a same-time frozen-normal-source split. The opt-in
//! [`timed`] composition advances the original thermal timeline between sealed
//! normal events. Neither path invents radiation, dynamic fracture or a
//! temperature-dependent normal potential.

use fs_blake3::ContentHash;
use fs_contact::normal_patch::adhesive::{
    JkrBranch, JkrError, JkrJump, JkrState, SourceJkrSphereContact, motion::JkrNormalMotionState,
};
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_thermal::enthalpy::{
    CellGeometry, CellId, DepositionError, DepositionSimulation, EnergyAudit, MaterialId,
};
use fs_tribo::{DissipationStep, HeatPartition, TriboError, WorkLedger};
use std::collections::BTreeMap;

mod normal_motion;
/// Original thermal transport between the sealed normal trajectory's events.
pub mod timed;
/// Explicit viscous-support losses on their original smooth mechanical intervals.
pub mod viscous;
/// Preserve original snap and caloric history across prescribed source temperature changes.
pub mod temperature;

/// Finite work and heat accuracy budgets, all chosen by the caller.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrSnapHeatControls {
    /// Maximum number of original monotone path segments in one transaction.
    pub maximum_segments: usize,
    /// Maximum represented heat defect in each original caloric cell [J].
    pub cell_energy_tolerance_j: f64,
    /// Maximum sum of absolute cell heat defects [J].
    pub thermal_energy_tolerance_j: f64,
    /// Maximum defect in heat plus named outside work minus snap release [J].
    pub transfer_energy_tolerance_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum JkrSnapHeatError {
    Invalid(&'static str),
    HistoryMismatch,
    Normal(JkrError),
    Partition(TriboError),
    Thermal(DepositionError),
    TransferBalance(f64),
}
impl From<JkrError> for JkrSnapHeatError {
    fn from(error: JkrError) -> Self {
        Self::Normal(error)
    }
}
impl From<TriboError> for JkrSnapHeatError {
    fn from(error: TriboError) -> Self {
        Self::Partition(error)
    }
}
impl From<DepositionError> for JkrSnapHeatError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl core::fmt::Display for JkrSnapHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "JKR snap heat: {self:?}")
    }
}
impl std::error::Error for JkrSnapHeatError {}

/// Sealed source/history version. An accepted version cannot be applied twice.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrSnapHeatStamp {
    source: ContentHash,
    revision: u64,
    state: JkrState,
}
impl JkrSnapHeatStamp {
    /// Original source-bound branch state at this accepted version.
    #[must_use]
    pub const fn state(self) -> JkrState {
        self.state
    }
}

#[derive(Clone, Copy)]
struct Receiver {
    cell: CellId,
    segment: u64,
    body: u64,
    material: MaterialId,
    mass_kg: f64,
    geometry: CellGeometry,
}

/// Accepted original branch and its explicit surface/outside work destinations.
///
/// The caloric objects must belong to the same ordered material cards as the
/// original normal source and be the exact objects registered in the timeline.
/// Fixed cell identity, mass and geometry remain pinned across proposals. This
/// is a heat-only mapping: no fictitious conducting area or contact link is
/// created when a detached surface receives its snap-off heat.
pub struct JkrSnapHeat<'s> {
    source: &'s SourceJkrSphereContact,
    caloric: [&'s EquilibriumEnthalpyPhaseCurve; 2],
    receivers: [Receiver; 2],
    partition: HeatPartition,
    outside_destination: Option<String>,
    stamp: JkrSnapHeatStamp,
    normal_motion_state: Option<JkrNormalMotionState>,
    accepted_time_s: f64,
    ledger: WorkLedger,
}

/// Observations from the original release and actual native thermal update.
#[derive(Clone, Debug, PartialEq)]
pub struct JkrSnapHeatReport {
    pub source_identity: ContentHash,
    pub initial_state: JkrState,
    pub final_state: JkrState,
    /// Original jumps, in physical traversal order, never caller-authored work.
    pub jumps: Vec<JkrJump>,
    /// This transaction's release and declared surface/outside shares [J].
    pub work: WorkLedger,
    /// Named destination of `work.other_work_j()`, outside represented cells.
    pub outside_destination: Option<String>,
    /// Native same-time heat audit. Its external heat is INTERNAL to the
    /// combined normal-contact/thermal/outside energy account.
    pub thermal: EnergyAudit,
    /// Actual caloric storage increase plus outside work minus original release.
    /// This is a transfer check, not a total-system thermodynamic certificate.
    pub transfer_residual_j: f64,
}

/// Both candidate owners publish together after every source and heat gate.
pub struct JkrSnapHeatCandidate<'s, 'm> {
    history: JkrSnapHeat<'s>,
    thermal: DepositionSimulation<'m>,
    report: JkrSnapHeatReport,
}
impl<'s, 'm> JkrSnapHeatCandidate<'s, 'm> {
    #[must_use]
    pub fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    #[must_use]
    pub const fn report(&self) -> &JkrSnapHeatReport {
        &self.report
    }
    /// Consume the complete candidate for a larger native mechanical transaction.
    #[must_use]
    pub fn into_parts(self) -> (JkrSnapHeat<'s>, DepositionSimulation<'m>, JkrSnapHeatReport) {
        (self.history, self.thermal, self.report)
    }
}

impl<'s> JkrSnapHeat<'s> {
    /// Pin the original branch, ordered material objects and finite receivers.
    ///
    /// A nonzero outside share requires a nonempty physical destination name.
    /// No default heat fraction or all-heat assumption is supplied. The initial
    /// state may be attached or detached; it must come from this exact source.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        source: &'s SourceJkrSphereContact,
        initial: JkrState,
        thermal: &DepositionSimulation<'_>,
        cells: [CellId; 2],
        caloric: [&'s EquilibriumEnthalpyPhaseCurve; 2],
        partition: HeatPartition,
        outside_destination: Option<String>,
    ) -> Result<Self, JkrSnapHeatError> {
        if cells[0] == cells[1]
            || outside_destination
                .as_ref()
                .is_some_and(|name| name.trim().is_empty())
            || (partition.other() > 0.0 && outside_destination.is_none())
        {
            return Err(JkrSnapHeatError::Invalid(
                "distinct cells and named outside destination required",
            ));
        }
        // A zero-work validation does not invent a time or physical dissipation.
        DissipationStep::from_work(0.0, partition)?;
        let regenerated = match initial.branch() {
            JkrBranch::Detached => source.detached_state(initial.indentation_m())?,
            JkrBranch::Attached => source.evaluate_attached(initial.indentation_m())?.state(),
        };
        if regenerated != initial {
            return Err(JkrSnapHeatError::HistoryMismatch);
        }
        let elastic = source.normal_state().elastic();
        let card = [
            elastic.surface_a_state().resolved().card_identity(),
            elastic.surface_b_state().resolved().card_identity(),
        ];
        let receiver = |side: usize| -> Result<Receiver, JkrSnapHeatError> {
            let cell = thermal
                .cells()
                .get(&cells[side])
                .ok_or(JkrSnapHeatError::Invalid("inactive original heat receiver"))?;
            if caloric[side].material_card_identity() != card[side]
                || !thermal.uses_material_law(cell.material_id, caloric[side])
            {
                return Err(JkrSnapHeatError::Invalid(
                    "normal and caloric owners need the same ordered original material cards",
                ));
            }
            Ok(Receiver {
                cell: cell.id,
                segment: cell.segment_id,
                body: cell.body_id,
                material: cell.material_id,
                mass_kg: cell.mass_kg,
                geometry: cell.geometry,
            })
        };
        let receivers = [receiver(0)?, receiver(1)?];
        Ok(Self {
            source,
            caloric,
            receivers,
            partition,
            outside_destination,
            stamp: JkrSnapHeatStamp {
                source: source.identity(),
                revision: 0,
                state: initial,
            },
            normal_motion_state: None,
            accepted_time_s: thermal.time_s(),
            ledger: WorkLedger::default(),
        })
    }

    #[must_use]
    pub const fn stamp(&self) -> JkrSnapHeatStamp {
        self.stamp
    }
    #[must_use]
    pub const fn state(&self) -> JkrState {
        self.stamp.state
    }
    #[must_use]
    pub const fn ledger(&self) -> &WorkLedger {
        &self.ledger
    }

    /// Regenerate an ordered accepted path and stage one native heat impulse.
    ///
    /// Each entry ends one monotone segment. The last regenerated sealed state
    /// must equal `expected_final`, allowing the native mechanical owner to bind
    /// its accepted endpoint to exactly the same branch traversal. Reversals need
    /// separate endpoints; an endpoint pair alone cannot describe a hidden loop.
    ///
    /// Every jump is taken from the original source, including an exact fold.
    /// Smooth elastic segments deposit exactly zero. The thermal state may have
    /// advanced since the previous accepted transfer, but cannot move backward.
    /// Once `propose_normal_motion` binds a sealed motion continuation, use that
    /// method for later intervals; prescribed endpoints cannot erase its clock.
    /// All releases in this batch enter at its CURRENT thermal time; no thermal
    /// evolution between batch events is implied. Call once per event when that
    /// intervening transport matters. Cost is O(segments + native cell clone),
    /// with one caloric application and bounded cancellation checks per segment.
    #[allow(clippy::too_many_arguments)]
    pub fn propose<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        expected: JkrSnapHeatStamp,
        endpoints_m: &[f64],
        expected_final: JkrState,
        controls: JkrSnapHeatControls,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<JkrSnapHeatCandidate<'s, 'm>, JkrSnapHeatError> {
        if self.normal_motion_state.is_some() {
            return Err(JkrSnapHeatError::Invalid(
                "normal-motion heat history requires its sealed continuation",
            ));
        }
        self.propose_segments(
            thermal,
            expected,
            endpoints_m,
            expected_final,
            controls,
            checkpoint,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn propose_segments<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        expected: JkrSnapHeatStamp,
        endpoints_m: &[f64],
        expected_final: JkrState,
        controls: JkrSnapHeatControls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<JkrSnapHeatCandidate<'s, 'm>, JkrSnapHeatError> {
        checkpoint()?;
        if expected != self.stamp {
            return Err(JkrSnapHeatError::HistoryMismatch);
        }
        if endpoints_m.is_empty()
            || endpoints_m.len() > controls.maximum_segments
            || [
                controls.cell_energy_tolerance_j,
                controls.thermal_energy_tolerance_j,
                controls.transfer_energy_tolerance_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
            || thermal.time_s() < self.accepted_time_s
        {
            return Err(JkrSnapHeatError::Invalid(
                "invalid segment, clock or heat budget",
            ));
        }
        let revision = self
            .stamp
            .revision
            .checked_add(1)
            .ok_or(JkrSnapHeatError::Invalid("history version exhausted"))?;
        self.admit_receivers(thermal)?;
        let mut state = self.stamp.state;
        let mut work = WorkLedger::default();
        let mut ledger = self.ledger;
        let mut jumps = Vec::new();
        for &end in endpoints_m {
            checkpoint()?;
            let transition = self.source.advance_monotone(state, end)?;
            state = transition.response.state();
            if let Some(jump) = transition.jump {
                let step =
                    DissipationStep::from_work(transition.irreversible_work_j, self.partition)?;
                work.record(step)?;
                ledger.record(step)?;
                jumps.push(jump);
            }
        }
        if state != expected_final {
            return Err(JkrSnapHeatError::HistoryMismatch);
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
        // Keep the small cell representation error separate from the larger
        // release, rather than subtracting absolute thermal reference energies.
        let partition_residual = (work.surface_a_heat_j() + work.surface_b_heat_j())
            + (work.other_work_j() - work.dissipated_work_j());
        let residual = partition_residual + audit.residual_j;
        if !residual.is_finite() || residual.abs() > controls.transfer_energy_tolerance_j {
            return Err(JkrSnapHeatError::TransferBalance(residual));
        }
        let history = Self {
            source: self.source,
            caloric: self.caloric,
            receivers: self.receivers,
            partition: self.partition,
            outside_destination: self.outside_destination.clone(),
            stamp: JkrSnapHeatStamp {
                source: self.source.identity(),
                revision,
                state,
            },
            normal_motion_state: self.normal_motion_state,
            accepted_time_s: thermal.time_s(),
            ledger,
        };
        let report = JkrSnapHeatReport {
            source_identity: self.source.identity(),
            initial_state: self.stamp.state,
            final_state: state,
            jumps,
            work,
            outside_destination: self.outside_destination.clone(),
            thermal: audit,
            transfer_residual_j: residual,
        };
        checkpoint()?;
        Ok(JkrSnapHeatCandidate {
            history,
            thermal,
            report,
        })
    }

    /// Publish accepted branch history and native caloric cells atomically.
    /// Any source failure, stale version, material-domain failure, energy refusal
    /// or cancellation leaves both accepted owners unchanged and permits retry.
    #[allow(clippy::too_many_arguments)]
    pub fn apply<'m>(
        &mut self,
        thermal: &mut DepositionSimulation<'m>,
        expected: JkrSnapHeatStamp,
        endpoints_m: &[f64],
        expected_final: JkrState,
        controls: JkrSnapHeatControls,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<JkrSnapHeatReport, JkrSnapHeatError> {
        let (history, next, report) = self
            .propose(
                thermal,
                expected,
                endpoints_m,
                expected_final,
                controls,
                checkpoint,
            )?
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
