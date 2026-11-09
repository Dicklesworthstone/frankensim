//! Rolling plastic heat transported over an actual finite thermal interval.
//!
//! The prescribed normal/area update and its original releases occur at the
//! interval start, after native scheduled events at that time. The original
//! rolling increment then supplies an interval-average plastic heat source.
//! Native conduction, reservoirs, ordinary sources and scheduled changes act
//! during that interval. There is no endpoint heat correction or guessed loss.
//!
//! This is a frozen-source time split with average-power reconstruction, not
//! a resolved within-step rolling trajectory or temperature feedback. Refine
//! the mechanical and thermal intervals for transient accuracy. Contact heat
//! partition and thermal topology remain explicit caller-selected physics.

/// Actual driven finite-body motion and thermal transport on one accepted clock.
pub mod native;

use super::super::super::power::{SurfacePower, sum_audits};
use super::{
    JkrRollingDissipationHeat, JkrRollingDissipationHeatReport, JkrRollingSnapHeatError,
    JkrRollingSnapHeatStamp, JkrRollingSnapSegment, JkrSnapHeatError,
};
use fs_contact::normal_patch::adhesive::rolling::events::JkrRollingEventReport;
use fs_math::dd::Dd;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation, EnergyAudit};
use fs_tribo::{DissipationStep, WorkLedger};

/// Reuses the original impulse, explicit-transport and whole-interval gates.
/// `heat.maximum_segments` must admit the normal and rolling source segments.
pub use super::super::super::timed::JkrTimedSnapHeatControls as JkrTimedRollingHeatControls;

/// Original normal impulse, distributed rolling heat, and native thermal work.
#[derive(Clone, Debug, PartialEq)]
pub struct JkrTimedRollingHeatReport {
    pub initial_time_s: f64,
    pub final_time_s: f64,
    /// Scheduled events at the initial time, before the normal release.
    pub initial_events: EnergyAudit,
    /// One original normal/area segment and its same-time caloric impulse.
    pub normal: JkrRollingDissipationHeatReport,
    /// Original held-area angular segment; recoverable work remains mechanical.
    pub rolling: JkrRollingEventReport,
    pub rolling_loss: WorkLedger,
    /// Actual time-integrated native source heat into the original A/B cells.
    pub injected_surface_heat_j: [f64; 2],
    /// Positive-duration native transport, including the rolling source.
    pub transport: EnergyAudit,
    /// Initial events, normal impulse and positive-duration transport together.
    pub thermal: EnergyAudit,
    /// Ordinary boundary/reservoir/source heat, excluding both contact channels.
    pub external_transport_heat_j: f64,
    pub rolling_transfer_residual_j: f64,
    /// Thermal storage + outside work - original contact losses - ordinary
    /// external heat - incoming material energy. Mechanical storage is separate.
    pub combined_energy_residual_j: f64,
    pub absolute_stage_energy_residual_j: f64,
}

/// Original source history and thermal clock advance as one accepted proposal.
pub struct JkrTimedRollingHeatCandidate<'s, 'm> {
    history: JkrRollingDissipationHeat<'s>,
    thermal: DepositionSimulation<'m>,
    report: JkrTimedRollingHeatReport,
}
impl<'s, 'm> JkrTimedRollingHeatCandidate<'s, 'm> {
    pub const fn history(&self) -> &JkrRollingDissipationHeat<'s> {
        &self.history
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub const fn report(&self) -> &JkrTimedRollingHeatReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        JkrRollingDissipationHeat<'s>,
        DepositionSimulation<'m>,
        JkrTimedRollingHeatReport,
    ) {
        (self.history, self.thermal, self.report)
    }
}

impl<'s> JkrRollingDissipationHeat<'s> {
    /// Apply a normal/area change at the start, then transport the source-owned
    /// rolling loss over `dt_s`. The complete monotone angular increment is
    /// prescribed; only its original irreversible loss supplies heat.
    ///
    /// The accepted heat and thermal clocks must match exactly. Native events
    /// at the start execute before release; later events retain their original
    /// order and source power is integrated over actual accepted durations.
    /// Rounding of the physical clock or average watts must satisfy the explicit
    /// transfer gate; no endpoint impulse repairs a discrepancy. All failures
    /// discard the proposed history, transport and caloric changes together.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn propose_interval<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        expected: JkrRollingSnapHeatStamp,
        indentation_m: f64,
        angle_increment_rad: f64,
        dt_s: f64,
        controls: JkrTimedRollingHeatControls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<JkrTimedRollingHeatCandidate<'s, 'm>, JkrRollingSnapHeatError> {
        checkpoint()?;
        let start = thermal.time_s();
        let end = finite(start + dt_s)?;
        if expected != self.snap.stamp || start != self.snap.accepted_time_s {
            return Err(JkrSnapHeatError::HistoryMismatch.into());
        }
        if !dt_s.is_finite()
            || dt_s <= 0.0
            || end <= start
            || controls.heat.maximum_segments < 2
            || controls.maximum_transport_substeps == 0
            || !controls.total_energy_tolerance_j.is_finite()
            || controls.total_energy_tolerance_j < 0.0
        {
            return Err(invalid(
                "rolling heat needs synchronized positive physical intervals",
            ));
        }
        self.snap.admit_receivers(thermal)?;
        let (at_start, (), initial_events) = thermal.advance_observed_limited(
            &(),
            start,
            controls.transport,
            controls.maximum_transport_substeps,
            &mut checkpoint,
        )?;
        self.snap.admit_receivers(&at_start)?;
        let normal = self.propose(
            &at_start,
            expected,
            &[JkrRollingSnapSegment::NormalIndentation(indentation_m)],
            controls.heat,
            &mut checkpoint,
        )?;
        let (mut history, at_start, normal) = normal.into_parts();
        let (source_after, rolling) = history
            .source()
            .advance_rolling(history.state(), angle_increment_rad)?
            .into_parts();
        let loss = DissipationStep::from_work(
            rolling.irreversible_rolling_loss_j,
            history.rolling_partition,
        )?;
        let mut rolling_loss = WorkLedger::default();
        rolling_loss.record(loss)?;
        history.rolling_loss_ledger.record(loss)?;
        let admit = |thermal: &DepositionSimulation<'_>| history.snap.admit_receivers(thermal);
        let observer = SurfacePower::new(
            history.snap.receivers.map(|receiver| receiver.cell),
            start,
            end,
            dt_s,
            loss,
            &admit,
        )?;
        let remaining = controls
            .maximum_transport_substeps
            .checked_sub(initial_events.substeps)
            .ok_or(DepositionError::Budget("rolling heat transport substeps"))?;
        let (thermal, observer, transport) = at_start.advance_observed_limited(
            &observer,
            end,
            controls.transport,
            remaining,
            &mut checkpoint,
        )?;
        let (injected, transfer, transfer_absolute) =
            observer.finish(loss, controls.heat.transfer_energy_tolerance_j)?;
        history.snap.admit_receivers(&thermal)?;
        let total = sum_audits([initial_events, normal.thermal, transport])?;
        let external = Dd::from_f64(initial_events.external_heat_j)
            + Dd::from_f64(transport.external_heat_j)
            - injected[0]
            - injected[1];
        // The original audits retain actual relative-storage defects. Summing
        // those and the independently checked transfer defects is the complete
        // thermal/contact account without re-subtracting rounded aggregate
        // storage, incoming material energy or large background heat totals.
        let combined = (Dd::from_f64(initial_events.residual_j)
            + Dd::from_f64(normal.transfer_residual_j)
            + Dd::from_f64(transport.residual_j)
            + Dd::from_f64(transfer))
        .to_f64();
        // Keep normal partition and actual caloric defects independent in the
        // absolute gate, even if their signs cancel in the signed account.
        let normal_source = normal.segments[0];
        let normal_partition = Dd::from_f64(normal.work.surface_a_heat_j())
            + Dd::from_f64(normal.work.surface_b_heat_j())
            + Dd::from_f64(normal.work.other_work_j())
            - Dd::from_f64(normal_source.normal_released_energy_j)
            - Dd::from_f64(normal_source.rolling_released_energy_j)
            - Dd::from_f64(normal_source.irreversible_rolling_loss_j);
        let absolute = (Dd::from_f64(initial_events.residual_j.abs())
            + Dd::from_f64(normal.thermal.residual_j.abs())
            + normal_partition.abs()
            + Dd::from_f64(transport.residual_j.abs())
            + transfer_absolute)
            .to_f64();
        if !combined.is_finite()
            || !absolute.is_finite()
            || combined.abs() > controls.total_energy_tolerance_j
            || absolute > controls.total_energy_tolerance_j
        {
            return Err(DepositionError::EnergyBalance(combined.abs().max(absolute)).into());
        }
        // The normal proposal's revision denotes this one enclosing accepted
        // interval; its intermediate source/clock has never been published.
        history.snap.stamp.state = source_after;
        history.snap.accepted_time_s = end;
        checkpoint()?;
        Ok(JkrTimedRollingHeatCandidate {
            history,
            thermal,
            report: JkrTimedRollingHeatReport {
                initial_time_s: start,
                final_time_s: end,
                initial_events,
                normal,
                rolling,
                rolling_loss,
                injected_surface_heat_j: [
                    finite(injected[0].to_f64())?,
                    finite(injected[1].to_f64())?,
                ],
                transport,
                thermal: total,
                external_transport_heat_j: finite(external.to_f64())?,
                rolling_transfer_residual_j: transfer,
                combined_energy_residual_j: combined,
                absolute_stage_energy_residual_j: absolute,
            },
        })
    }

    /// Publish source history, cumulative losses and the actual thermal timeline.
    #[allow(clippy::too_many_arguments)]
    pub fn apply_interval<'m>(
        &mut self,
        thermal: &mut DepositionSimulation<'m>,
        expected: JkrRollingSnapHeatStamp,
        indentation_m: f64,
        angle_increment_rad: f64,
        dt_s: f64,
        controls: JkrTimedRollingHeatControls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<JkrTimedRollingHeatReport, JkrRollingSnapHeatError> {
        let (history, next, report) = self
            .propose_interval(
                thermal,
                expected,
                indentation_m,
                angle_increment_rad,
                dt_s,
                controls,
                &mut checkpoint,
            )?
            .into_parts();
        checkpoint()?;
        *self = history;
        *thermal = next;
        Ok(report)
    }
}

fn invalid(message: &'static str) -> JkrRollingSnapHeatError {
    JkrSnapHeatError::Invalid(message).into()
}
fn finite(value: f64) -> Result<f64, JkrRollingSnapHeatError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic.into())
    }
}
