//! Event-timed JKR snap heat on the original native thermal timeline.
//!
//! Ordinary conduction, reservoirs, sources and scheduled events remain owned
//! by `DepositionSimulation::advance_to`. A source-regenerated snap is applied
//! after all native scheduled events at its timestamp, before any later heat
//! transport. No interpolation of temperatures or replacement heat solver is
//! introduced. Normal coefficients and the declared heat partition stay frozen.

use super::{
    JkrSnapHeat, JkrSnapHeatControls, JkrSnapHeatError, JkrSnapHeatReport, JkrSnapHeatStamp,
};
use fs_contact::normal_patch::adhesive::{JkrJump, motion::JkrNormalMotionStep};
use fs_math::dd::Dd;
use fs_thermal::enthalpy::{
    DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit, StepConfig,
};
use fs_tribo::{DissipationStep, WorkLedger};
use std::collections::BTreeMap;

/// Original heat and transport admissions plus whole-interval acceptance gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrTimedSnapHeatControls {
    /// Original per-impulse cell/heat gates and aggregate snap-transfer gate.
    /// `maximum_segments` bounds original events plus the final branch endpoint.
    pub heat: JkrSnapHeatControls,
    /// Native explicit thermal stepping and its per-advance energy gate.
    pub transport: StepConfig,
    /// Shared native substep allowance for the complete proposal. Each march
    /// receives only the remainder and refuses before further heat evaluation.
    /// The timeline also retains its own original per-advance work ceiling.
    pub maximum_transport_substeps: usize,
    /// Sum of absolute transport, impulse and partition residuals [J], also
    /// bounding the combined transfer/transport residual for this interval.
    pub total_energy_tolerance_j: f64,
}

/// Original transport followed by one original snap, at its physical timestamp.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrTimedSnapHeatEvent {
    pub motion_offset_s: f64,
    pub thermal_time_s: f64,
    pub jump: JkrJump,
    /// Transport since the preceding snap (or the initial accepted state).
    /// Scheduled events at this timestamp have already been processed.
    pub transport: EnergyAudit,
    /// Native same-time caloric application of this snap's surface shares.
    pub impulse: EnergyAudit,
}

#[derive(Clone, Debug, PartialEq)]
pub struct JkrTimedSnapHeatReport {
    pub initial_time_s: f64,
    pub final_time_s: f64,
    /// Explicit clock mapping: thermal time = this origin + motion time.
    pub thermal_time_at_motion_zero_s: f64,
    pub events: Vec<JkrTimedSnapHeatEvent>,
    /// Native march from the last snap to the interval endpoint. With no snap,
    /// this is the complete ordinary thermal advance.
    pub final_transport: EnergyAudit,
    /// Sum of native ordinary transport audits, excluding snap deposits.
    pub transport: EnergyAudit,
    /// Original snap work, destinations and the sum of caloric impulse audits.
    /// `snap.thermal.external_heat_j` is internal contact-to-thermal transfer.
    pub snap: JkrSnapHeatReport,
    /// Complete native thermal account, including transport, births and snaps.
    pub thermal: EnergyAudit,
    /// Thermal storage change + outside snap work - original snap release
    /// - transport heat - incoming material energy. It uses native relative
    /// storage audits, avoiding subtraction of absolute enthalpy datums.
    pub combined_energy_residual_j: f64,
    pub absolute_stage_energy_residual_j: f64,
}

/// Original thermal state and full motion/heat history accept together.
pub struct JkrTimedSnapHeatCandidate<'s, 'm> {
    history: JkrSnapHeat<'s>,
    thermal: DepositionSimulation<'m>,
    report: JkrTimedSnapHeatReport,
}
impl<'s, 'm> JkrTimedSnapHeatCandidate<'s, 'm> {
    #[must_use]
    pub fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    #[must_use]
    pub const fn report(&self) -> &JkrTimedSnapHeatReport {
        &self.report
    }
    /// Consume only after every other coupled native owner admits its candidate.
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        JkrSnapHeat<'s>,
        DepositionSimulation<'m>,
        JkrTimedSnapHeatReport,
    ) {
        (self.history, self.thermal, self.report)
    }
}

impl<'s> JkrSnapHeat<'s> {
    /// Advance native heat transport between the sealed motion's actual snaps.
    ///
    /// The caller explicitly maps the mechanical clock: `thermal_time =
    /// thermal_time_at_motion_zero_s + motion_time`. Pass zero for identical
    /// clocks. The supplied thermal state must be exactly at the mapped initial
    /// time. Once motion is bound, its full source/mass/policy/time/step token
    /// and the preceding accepted thermal time must both continue exactly.
    /// Nonfinite, reversed or unrepresentable event clocks refuse; this method
    /// never shifts an event to the current or final thermal timestamp.
    ///
    /// Every native scheduled birth, boundary or contact event at a snap time
    /// executes first, using the existing schedule order. Then original released
    /// work enters the same pinned caloric cells. Remaining transport uses those
    /// heated states. Scheduled topology changes remain caller-authored; a snap
    /// does not manufacture conducting area or remove a thermal link by itself.
    ///
    /// Only sealed source jumps create heat; edited public mechanical reports
    /// and smooth reversible work are ignored. Normal coefficients and the heat
    /// partition are held over this split. This is event-timed one-way heat
    /// coupling, not temperature-dependent normal feedback or a wave solver.
    /// Native explicit transport retains its spatial and timestep limitations.
    ///
    /// All work is private until the complete candidate is returned. A rejected
    /// source, caloric/transport domain, budget or cancellation changes neither
    /// accepted timeline nor branch history. This candidate can be staged in
    /// the native mechanical `step_with_admission` callback.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn propose_timed_normal_motion<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        expected: JkrSnapHeatStamp,
        motion: &JkrNormalMotionStep,
        thermal_time_at_motion_zero_s: f64,
        controls: JkrTimedSnapHeatControls,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<JkrTimedSnapHeatCandidate<'s, 'm>, JkrSnapHeatError> {
        self.propose_timed_normal_motion_observed(
            thermal,
            expected,
            motion,
            thermal_time_at_motion_zero_s,
            controls,
            &(),
            checkpoint,
        )
        .map(|(candidate, ())| candidate)
    }

    // A narrowly shared transport path for original, source-owned smooth losses.
    // Public snap-only proposals keep the original no-op observer and arithmetic.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(super) fn propose_timed_normal_motion_observed<'m, O>(
        &self,
        thermal: &DepositionSimulation<'m>,
        expected: JkrSnapHeatStamp,
        motion: &JkrNormalMotionStep,
        thermal_time_at_motion_zero_s: f64,
        controls: JkrTimedSnapHeatControls,
        observer: &O,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(JkrTimedSnapHeatCandidate<'s, 'm>, O), JkrSnapHeatError>
    where
        O: DepositionObserver,
        JkrSnapHeatError: From<O::Error>,
    {
        checkpoint()?;
        if expected != self.stamp
            || motion.initial_state().contact_state() != self.state()
            || self
                .normal_motion_state
                .is_some_and(|previous| previous != motion.initial_state())
        {
            return Err(JkrSnapHeatError::HistoryMismatch);
        }
        let initial_motion_time = motion.initial_state().time_s();
        let final_motion_time = motion.state().time_s();
        let start = finite(thermal_time_at_motion_zero_s + initial_motion_time)?;
        let end = finite(thermal_time_at_motion_zero_s + final_motion_time)?;
        if !thermal_time_at_motion_zero_s.is_finite()
            || start != thermal.time_s()
            || start < self.accepted_time_s
            || start < 0.0
            || end <= start
            || (self.normal_motion_state.is_some() && start != self.accepted_time_s)
            || motion.event_count() >= controls.heat.maximum_segments
            || controls.maximum_transport_substeps == 0
            || !controls.transport.max_dt_s.is_finite()
            || controls.transport.max_dt_s <= 0.0
            || [
                controls.heat.cell_energy_tolerance_j,
                controls.heat.thermal_energy_tolerance_j,
                controls.heat.transfer_energy_tolerance_j,
                controls.transport.energy_tolerance_j,
                controls.total_energy_tolerance_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(invalid(
                "timed snap heat needs matching explicit clocks and finite positive work budgets",
            ));
        }
        let revision = self
            .stamp
            .revision
            .checked_add(1)
            .ok_or_else(|| invalid("history version exhausted"))?;
        self.admit_receivers(thermal)?;

        // Regenerate the sealed event path before cloning a thermal candidate.
        // Within-branch reversals affect mechanical work, never the snap debit.
        let mut state = self.state();
        let mut work = WorkLedger::default();
        let mut ledger = self.ledger;
        let mut planned = Vec::with_capacity(motion.event_count());
        let mut previous_offset = 0.0;
        let mut previous_time = start;
        let mut previous_motion_time = initial_motion_time;
        for event in motion.events() {
            checkpoint()?;
            let event_motion_time = finite(initial_motion_time + event.offset_s)?;
            let time = finite(thermal_time_at_motion_zero_s + event_motion_time)?;
            if !event.offset_s.is_finite()
                || event.offset_s < previous_offset
                || time < previous_time
                || time > end
                || event_motion_time > final_motion_time
                || (event.offset_s > previous_offset && time <= previous_time)
            {
                return Err(invalid(
                    "snap times must remain ordered and representable on the explicit thermal clock",
                ));
            }
            let transition = self
                .source
                .advance_monotone(state, event.before().indentation_m())?;
            if transition.jump != Some(event.jump) || transition.response.state() != event.after() {
                return Err(JkrSnapHeatError::HistoryMismatch);
            }
            let partition =
                DissipationStep::from_work(transition.irreversible_work_j, self.partition)?;
            work.record(partition)?;
            ledger.record(partition)?;
            planned.push((event.offset_s, time, event.jump, partition));
            state = event.after();
            previous_offset = event.offset_s;
            previous_time = time;
            previous_motion_time = event_motion_time;
        }
        if final_motion_time > previous_motion_time && end <= previous_time {
            return Err(invalid(
                "final transport interval is not representable on the explicit thermal clock",
            ));
        }
        let final_branch = self
            .source
            .advance_monotone(state, motion.state().contact_state().indentation_m())?;
        if final_branch.jump.is_some()
            || final_branch.response.state() != motion.state().contact_state()
        {
            return Err(JkrSnapHeatError::HistoryMismatch);
        }

        let mut current = None;
        let mut observer = observer.clone();
        let mut observer_begun = false;
        let mut events = Vec::with_capacity(planned.len());
        let mut marches = AuditSum::default();
        let mut impulses = AuditSum::default();
        let mut absolute_residual = Dd::ZERO;
        for (offset, time, jump, partition) in planned {
            checkpoint()?;
            let (advanced, next_observer, march) = advance_transport(
                current.as_ref().unwrap_or(thermal),
                &observer,
                &mut observer_begun,
                time,
                controls.transport,
                transport_budget(&marches, controls.maximum_transport_substeps)?,
                &mut checkpoint,
            )?;
            observer = next_observer;
            self.admit_receivers(&advanced)?;
            marches.add(march)?;
            transport_budget(&marches, controls.maximum_transport_substeps)?;
            let heat = BTreeMap::from([
                (self.receivers[0].cell, partition.surface_a_heat_j()),
                (self.receivers[1].cell, partition.surface_b_heat_j()),
            ]);
            let (heated, impulse) = advanced.apply_heat_j(
                &heat,
                controls.heat.cell_energy_tolerance_j,
                controls.heat.thermal_energy_tolerance_j,
                &mut checkpoint,
            )?;
            let partition_residual = partition_defect(partition);
            let transfer_residual =
                finite((partition_residual + Dd::from_f64(impulse.residual_j)).to_f64())?;
            if transfer_residual.abs() > controls.heat.transfer_energy_tolerance_j {
                return Err(JkrSnapHeatError::TransferBalance(transfer_residual));
            }
            absolute_residual = absolute_residual
                + Dd::from_f64(march.residual_j.abs())
                + Dd::from_f64(impulse.residual_j.abs())
                + partition_residual.abs();
            impulses.add(impulse)?;
            events.push(JkrTimedSnapHeatEvent {
                motion_offset_s: offset,
                thermal_time_s: time,
                jump,
                transport: march,
                impulse,
            });
            current = Some(heated);
        }
        let (thermal, observer, final_transport) = advance_transport(
            current.as_ref().unwrap_or(thermal),
            &observer,
            &mut observer_begun,
            end,
            controls.transport,
            transport_budget(&marches, controls.maximum_transport_substeps)?,
            &mut checkpoint,
        )?;
        self.admit_receivers(&thermal)?;
        marches.add(final_transport)?;
        transport_budget(&marches, controls.maximum_transport_substeps)?;
        absolute_residual = absolute_residual + Dd::from_f64(final_transport.residual_j.abs());
        let transport = marches.audit()?;
        let impulse = impulses.audit()?;
        let partition_residual = Dd::from_f64(work.surface_a_heat_j())
            + Dd::from_f64(work.surface_b_heat_j())
            + Dd::from_f64(work.other_work_j())
            - Dd::from_f64(work.dissipated_work_j());
        let snap_residual =
            finite((partition_residual + Dd::from_f64(impulse.residual_j)).to_f64())?;
        if snap_residual.abs() > controls.heat.transfer_energy_tolerance_j {
            return Err(JkrSnapHeatError::TransferBalance(snap_residual));
        }
        let combined =
            finite((Dd::from_f64(transport.residual_j) + Dd::from_f64(snap_residual)).to_f64())?;
        let absolute = finite(absolute_residual.to_f64())?;
        if combined.abs() > controls.total_energy_tolerance_j
            || absolute > controls.total_energy_tolerance_j
        {
            return Err(DepositionError::EnergyBalance(absolute.max(combined.abs())).into());
        }
        let mut total = marches;
        total.add(impulse)?;
        let snap = JkrSnapHeatReport {
            source_identity: self.source.identity(),
            initial_state: self.state(),
            final_state: motion.state().contact_state(),
            jumps: events.iter().map(|e| e.jump).collect(),
            work,
            outside_destination: self.outside_destination.clone(),
            thermal: impulse,
            transfer_residual_j: snap_residual,
        };
        let history = JkrSnapHeat {
            source: self.source,
            caloric: self.caloric,
            receivers: self.receivers,
            partition: self.partition,
            outside_destination: self.outside_destination.clone(),
            stamp: JkrSnapHeatStamp {
                source: self.source.identity(),
                revision,
                state: snap.final_state,
            },
            normal_motion_state: Some(motion.state()),
            accepted_time_s: thermal.time_s(),
            ledger,
        };
        let report = JkrTimedSnapHeatReport {
            initial_time_s: start,
            final_time_s: end,
            thermal_time_at_motion_zero_s,
            events,
            final_transport,
            transport,
            snap,
            thermal: total.audit()?,
            combined_energy_residual_j: combined,
            absolute_stage_energy_residual_j: absolute,
        };
        checkpoint()?;
        Ok((
            JkrTimedSnapHeatCandidate {
                history,
                thermal,
                report,
            },
            observer,
        ))
    }
}

#[allow(clippy::too_many_arguments)]
fn advance_transport<'m, O>(
    thermal: &DepositionSimulation<'m>,
    observer: &O,
    begun: &mut bool,
    target_s: f64,
    config: StepConfig,
    remaining_substeps: usize,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<(DepositionSimulation<'m>, O, EnergyAudit), JkrSnapHeatError>
where
    O: DepositionObserver,
    JkrSnapHeatError: From<O::Error>,
{
    let result = if *begun {
        thermal.advance_observed_continuing_limited(
            observer,
            target_s,
            config,
            remaining_substeps,
            checkpoint,
        )
    } else {
        thermal.advance_observed_limited(observer, target_s, config, remaining_substeps, checkpoint)
    }?;
    *begun = true;
    Ok(result)
}

fn invalid(message: &'static str) -> JkrSnapHeatError {
    JkrSnapHeatError::Invalid(message)
}
fn finite(value: f64) -> Result<f64, JkrSnapHeatError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic.into())
    }
}
fn partition_defect(step: DissipationStep) -> Dd {
    Dd::from_f64(step.surface_a_heat_j())
        + Dd::from_f64(step.surface_b_heat_j())
        + Dd::from_f64(step.other_work_j())
        - Dd::from_f64(step.total_work_j())
}

#[derive(Default)]
struct AuditSum {
    // Native incoming mass, incoming energy, external heat, storage and residual.
    values: [Dd; 5],
    substeps: usize,
    events: usize,
}
impl AuditSum {
    fn add(&mut self, audit: EnergyAudit) -> Result<(), JkrSnapHeatError> {
        for (sum, value) in self.values.iter_mut().zip([
            audit.incoming_mass_kg,
            audit.incoming_energy_j,
            audit.external_heat_j,
            audit.stored_energy_change_j,
            audit.residual_j,
        ]) {
            *sum = *sum + Dd::from_f64(finite(value)?);
            finite(sum.to_f64())?;
        }
        self.substeps = self
            .substeps
            .checked_add(audit.substeps)
            .ok_or_else(|| invalid("thermal substep count overflow"))?;
        self.events = self
            .events
            .checked_add(audit.events)
            .ok_or_else(|| invalid("thermal event count overflow"))?;
        Ok(())
    }
    fn audit(&self) -> Result<EnergyAudit, JkrSnapHeatError> {
        Ok(EnergyAudit {
            incoming_mass_kg: finite(self.values[0].to_f64())?,
            incoming_energy_j: finite(self.values[1].to_f64())?,
            external_heat_j: finite(self.values[2].to_f64())?,
            stored_energy_change_j: finite(self.values[3].to_f64())?,
            residual_j: finite(self.values[4].to_f64())?,
            substeps: self.substeps,
            events: self.events,
        })
    }
}
fn transport_budget(sum: &AuditSum, maximum: usize) -> Result<usize, JkrSnapHeatError> {
    maximum
        .checked_sub(sum.substeps)
        .ok_or_else(|| DepositionError::Budget("timed snap heat total thermal substeps").into())
}
