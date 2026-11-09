//! Calibrated torsional loss and junction release on physical thermal time.
//!
//! Original thermal events at the start precede the prescribed normal/area
//! stage and its normal and retained-spring release. The original monotone
//! twist increment supplies only its plastic loss as interval-average power.
//! Native conduction, reservoirs, ordinary sources and scheduled changes then
//! evolve the same finite receivers. Source and thermal clocks accept together.
//!
//! This retains the original homogeneous bonded-twist model and its mode-III
//! applicability screen. Normal loading occurs at interval start and material
//! properties remain frozen. Average-power reconstruction requires temporal
//! refinement; it is not a resolved within-step twist or temperature feedback.

/// Original driven body motion and caloric transport on one accepted interval.
pub mod native;

use super::super::power::{SurfacePower, sum_audits};
use super::{
    JkrSnapHeatError, JkrTorsionHeat, JkrTorsionHeatError, JkrTorsionHeatReport,
    JkrTorsionHeatSegment, JkrTorsionHeatStamp,
};
use fs_contact::normal_patch::adhesive::torsion::events::JkrTorsionEventReport;
use fs_math::dd::Dd;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation, EnergyAudit};
use fs_tribo::{DissipationStep, WorkLedger};

/// Original impulse gates, finite transport budget and whole-interval balance.
/// `heat.maximum_segments` must admit the normal and angular source stages.
pub use super::super::timed::JkrTimedSnapHeatControls as JkrTimedTorsionHeatControls;

#[derive(Clone, Debug, PartialEq)]
pub struct JkrTimedTorsionHeatReport {
    pub initial_time_s: f64,
    pub final_time_s: f64,
    /// Native scheduled events at the initial time, before source release.
    pub initial_events: EnergyAudit,
    /// Original normal/area stage and both release channels at the start.
    pub normal: JkrTorsionHeatReport,
    /// Original angular stage; recoverable twist work remains mechanical.
    pub torsion: JkrTorsionEventReport,
    pub torsion_loss: WorkLedger,
    /// Actual time-integrated plastic heat into the ordered original cells.
    pub injected_surface_heat_j: [f64; 2],
    /// Positive-duration transport, including the distributed plastic source.
    pub transport: EnergyAudit,
    /// Initial native events, release impulse and positive-duration transport.
    pub thermal: EnergyAudit,
    /// Ordinary thermal supply only, excluding release and plastic heat.
    pub external_transport_heat_j: f64,
    pub torsion_transfer_residual_j: f64,
    /// Delta H + outside loss - original contact losses - ordinary heat and
    /// incoming material energy. Mechanical storage and work remain separate.
    pub combined_energy_residual_j: f64,
    /// Absolute transfer and original thermal defects cannot hide by cancelling.
    pub absolute_stage_energy_residual_j: f64,
}

pub struct JkrTimedTorsionHeatCandidate<'s, 'm> {
    history: JkrTorsionHeat<'s>,
    thermal: DepositionSimulation<'m>,
    report: JkrTimedTorsionHeatReport,
}
impl<'s, 'm> JkrTimedTorsionHeatCandidate<'s, 'm> {
    pub const fn history(&self) -> &JkrTorsionHeat<'s> {
        &self.history
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub const fn report(&self) -> &JkrTimedTorsionHeatReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        JkrTorsionHeat<'s>,
        DepositionSimulation<'m>,
        JkrTimedTorsionHeatReport,
    ) {
        (self.history, self.thermal, self.report)
    }
}

impl<'s> JkrTorsionHeat<'s> {
    /// Stage a normal event at interval start, then transport the original
    /// angular step's plastic heat over a positive physical duration.
    ///
    /// Both accepted clocks must match exactly. The shared thermal source
    /// integrates power over actual native step durations and checks each
    /// receiving-cell total. Clock or average-power rounding must satisfy the
    /// transfer gate; no endpoint heat correction repairs a discrepancy.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn propose_interval<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        expected: JkrTorsionHeatStamp,
        indentation_m: f64,
        angle_increment_rad: f64,
        dt_s: f64,
        controls: JkrTimedTorsionHeatControls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<JkrTimedTorsionHeatCandidate<'s, 'm>, JkrTorsionHeatError> {
        checkpoint()?;
        let start = thermal.time_s();
        let end = finite(start + dt_s)?;
        if expected != self.stamp || start != self.accepted_time_s {
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
            return Err(JkrSnapHeatError::Invalid(
                "torsion heat needs synchronized positive physical intervals",
            )
            .into());
        }
        self.admit_receivers(thermal)?;
        let (at_start, (), initial_events) = thermal.advance_observed_limited(
            &(),
            start,
            controls.transport,
            controls.maximum_transport_substeps,
            &mut checkpoint,
        )?;
        self.admit_receivers(&at_start)?;
        let normal = self.propose(
            &at_start,
            expected,
            &[JkrTorsionHeatSegment::NormalIndentation(indentation_m)],
            controls.heat,
            &mut checkpoint,
        )?;
        let (mut history, at_start, normal) = normal.into_parts();
        let (source_after, torsion) = history
            .source()
            .advance_torsion(history.state(), angle_increment_rad)?
            .into_parts();
        let loss = DissipationStep::from_work(
            torsion.irreversible_torsion_loss_j,
            history.torsion_partition,
        )?;
        let mut torsion_loss = WorkLedger::default();
        torsion_loss.record(loss)?;
        history.torsion_loss_ledger.record(loss)?;
        let admit = |thermal: &DepositionSimulation<'_>| history.admit_receivers(thermal);
        let observer = SurfacePower::new(
            history.receivers.map(|r| r.cell),
            start,
            end,
            dt_s,
            loss,
            &admit,
        )?;
        let remaining = controls
            .maximum_transport_substeps
            .checked_sub(initial_events.substeps)
            .ok_or(DepositionError::Budget("torsion heat transport substeps"))?;
        let (thermal, observer, transport) = at_start.advance_observed_limited(
            &observer,
            end,
            controls.transport,
            remaining,
            &mut checkpoint,
        )?;
        let (injected, transfer, transfer_absolute) =
            observer.finish(loss, controls.heat.transfer_energy_tolerance_j)?;
        history.admit_receivers(&thermal)?;
        let total = sum_audits([initial_events, normal.thermal, transport])?;
        let external = Dd::from_f64(initial_events.external_heat_j)
            + Dd::from_f64(transport.external_heat_j)
            - injected[0]
            - injected[1];
        // Preserve the original relative native storage defects. Subtracting
        // rounded absolute formation enthalpies would lose small contact heat.
        let combined = (Dd::from_f64(initial_events.residual_j)
            + Dd::from_f64(normal.transfer_residual_j)
            + Dd::from_f64(transport.residual_j)
            + Dd::from_f64(transfer))
        .to_f64();
        let absolute = (Dd::from_f64(initial_events.residual_j.abs())
            + Dd::from_f64(normal.absolute_transfer_residual_j)
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
        // One enclosing interval owns the normal proposal's new revision.
        // Its intermediate branch, heat and clock have never been published.
        history.stamp.state = source_after;
        history.accepted_time_s = end;
        checkpoint()?;
        Ok(JkrTimedTorsionHeatCandidate {
            history,
            thermal,
            report: JkrTimedTorsionHeatReport {
                initial_time_s: start,
                final_time_s: end,
                initial_events,
                normal,
                torsion,
                torsion_loss,
                injected_surface_heat_j: [
                    finite(injected[0].to_f64())?,
                    finite(injected[1].to_f64())?,
                ],
                transport,
                thermal: total,
                external_transport_heat_j: finite(external.to_f64())?,
                torsion_transfer_residual_j: transfer,
                combined_energy_residual_j: combined,
                absolute_stage_energy_residual_j: absolute,
            },
        })
    }

    /// Publish source history, all loss channels and the advanced thermal clock.
    #[allow(clippy::too_many_arguments)]
    pub fn apply_interval<'m>(
        &mut self,
        thermal: &mut DepositionSimulation<'m>,
        expected: JkrTorsionHeatStamp,
        indentation_m: f64,
        angle_increment_rad: f64,
        dt_s: f64,
        controls: JkrTimedTorsionHeatControls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<JkrTimedTorsionHeatReport, JkrTorsionHeatError> {
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

fn finite(value: f64) -> Result<f64, JkrTorsionHeatError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic.into())
    }
}
