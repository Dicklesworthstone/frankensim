//! Smooth viscous-support loss and adhesion snaps on one native heat timeline.
//!
//! The original mechanical solver retains a bounded, immutable trace of its
//! positive dissipative segments. Their constant average powers are separate
//! from original snap release, reversible contact work and elastic support
//! storage. Only an explicitly supplied support partition routes these powers
//! to the already-bound caloric cells. There is no inferred heat fraction,
//! additional thermal store, changed boundary or manufactured contact link.

use super::timed::{JkrTimedSnapHeatControls, JkrTimedSnapHeatReport};
use super::{JkrSnapHeat, JkrSnapHeatError, JkrSnapHeatStamp};
use fs_contact::normal_patch::adhesive::motion::{
    JkrNormalDissipationSegment, JkrNormalMotionTrace,
};
use fs_math::dd::Dd;
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, Event,
};
use fs_tribo::{DissipationStep, HeatPartition, WorkLedger};
use std::collections::BTreeMap;

/// Caller-selected destinations for the mechanical support's viscous loss.
///
/// Surface shares use the two original cells already bound by `JkrSnapHeat`.
/// Selecting them explicitly asserts that those cells receive this support's
/// loss; the adhesive snap partition provides no such material closure. A
/// nonzero outside share requires a named destination. An all-outside route
/// injects no heat and leaves native transport stepping unchanged.
#[derive(Clone, Debug, PartialEq)]
pub struct JkrViscousHeatRoute {
    partition: HeatPartition,
    outside_destination: Option<String>,
}
impl JkrViscousHeatRoute {
    pub fn new(
        partition: HeatPartition,
        outside_destination: Option<String>,
    ) -> Result<Self, JkrSnapHeatError> {
        if outside_destination
            .as_ref()
            .is_some_and(|name| name.trim().is_empty())
            || (partition.other() > 0.0 && outside_destination.is_none())
        {
            return Err(invalid(
                "viscous support needs a named nonempty outside-work destination",
            ));
        }
        Ok(Self {
            partition,
            outside_destination,
        })
    }
    #[must_use]
    pub const fn partition(&self) -> HeatPartition {
        self.partition
    }
    #[must_use]
    pub fn outside_destination(&self) -> Option<&str> {
        self.outside_destination.as_deref()
    }
}

/// Original passive work and its independently declared support destinations.
#[derive(Clone, Debug, PartialEq)]
pub struct JkrViscousHeatReport {
    /// Original sealed support dissipation, separate from snap release.
    pub work: WorkLedger,
    /// Actual time-integrated native source energy at each original cell [J].
    pub injected_surface_heat_j: [f64; 2],
    pub outside_destination: Option<String>,
    /// Actual source heat + declared outside work - original support loss [J].
    pub transfer_residual_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct JkrTimedViscousHeatReport {
    /// Original native audits. `timed.transport.external_heat_j` includes the
    /// support source, which is internal to the combined mechanical/heat system.
    pub timed: JkrTimedSnapHeatReport,
    pub support: JkrViscousHeatReport,
    /// Authored sources/reservoirs/boundaries, excluding support and snap heat.
    pub external_transport_heat_j: f64,
    /// Original relative transport/snap residual plus support transfer defect.
    pub combined_energy_residual_j: f64,
    pub absolute_stage_energy_residual_j: f64,
}

pub struct JkrTimedViscousHeatCandidate<'s, 'm> {
    history: JkrSnapHeat<'s>,
    thermal: DepositionSimulation<'m>,
    report: JkrTimedViscousHeatReport,
}
impl<'s, 'm> JkrTimedViscousHeatCandidate<'s, 'm> {
    #[must_use]
    pub fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    #[must_use]
    pub const fn report(&self) -> &JkrTimedViscousHeatReport {
        &self.report
    }
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        JkrSnapHeat<'s>,
        DepositionSimulation<'m>,
        JkrTimedViscousHeatReport,
    ) {
        (self.history, self.thermal, self.report)
    }
}

impl<'s> JkrSnapHeat<'s> {
    /// Transport original viscous-support heat throughout its smooth intervals,
    /// and original adhesion snap heat at the source's branch-event times.
    ///
    /// Only the original immutable trace supplies support work and timing.
    /// Editable public frame observations cannot create heat. The same full
    /// motion continuation, explicit clock mapping, receiver identities, shared
    /// thermal work budget and transactional publication as `timed` apply.
    ///
    /// The reconstruction is piecewise constant at each original mechanical
    /// segment's average dissipative power. The native thermal solver clips at
    /// every power discontinuity and conducts/cools the generated heat during
    /// that segment. It adds these sources to the original prepared heat rates;
    /// authored boundary sources, reservoirs and scheduled changes retain their
    /// own meanings. Mechanical and thermal timestep refinement remain necessary.
    ///
    /// Actual source energy uses the native accepted durations. Its difference
    /// from original mechanical work must satisfy the explicit transfer gate;
    /// clock rounding is never repaired by an invented endpoint heat impulse.
    /// Source coefficients and both separately declared partitions stay frozen.
    #[allow(clippy::too_many_arguments)]
    pub fn propose_timed_damped_normal_motion<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        expected: JkrSnapHeatStamp,
        trace: &JkrNormalMotionTrace,
        thermal_time_at_motion_zero_s: f64,
        route: &JkrViscousHeatRoute,
        controls: JkrTimedSnapHeatControls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<JkrTimedViscousHeatCandidate<'s, 'm>, JkrSnapHeatError> {
        checkpoint()?;
        let source = SupportSources::new(
            trace,
            thermal_time_at_motion_zero_s,
            route,
            [self.receivers[0].cell, self.receivers[1].cell],
            controls.heat.transfer_energy_tolerance_j,
            &mut checkpoint,
        )?;
        let (candidate, source) = self.propose_timed_normal_motion_observed(
            thermal,
            expected,
            trace.step(),
            thermal_time_at_motion_zero_s,
            controls,
            &source,
            &mut checkpoint,
        )?;
        let (support, support_absolute) = source.finish(route)?;
        let (history, thermal, timed) = candidate.into_parts();
        let injected = Dd::from_f64(support.injected_surface_heat_j[0])
            + Dd::from_f64(support.injected_surface_heat_j[1]);
        let external_transport_heat_j =
            finite((Dd::from_f64(timed.transport.external_heat_j) - injected).to_f64())?;
        let combined_energy_residual_j = finite(
            (Dd::from_f64(timed.combined_energy_residual_j)
                + Dd::from_f64(support.transfer_residual_j))
            .to_f64(),
        )?;
        let absolute_stage_energy_residual_j = finite(
            (Dd::from_f64(timed.absolute_stage_energy_residual_j) + support_absolute).to_f64(),
        )?;
        if combined_energy_residual_j.abs() > controls.total_energy_tolerance_j
            || absolute_stage_energy_residual_j > controls.total_energy_tolerance_j
        {
            return Err(DepositionError::EnergyBalance(
                combined_energy_residual_j
                    .abs()
                    .max(absolute_stage_energy_residual_j),
            )
            .into());
        }
        checkpoint()?;
        Ok(JkrTimedViscousHeatCandidate {
            history,
            thermal,
            report: JkrTimedViscousHeatReport {
                timed,
                support,
                external_transport_heat_j,
                combined_energy_residual_j,
                absolute_stage_energy_residual_j,
            },
        })
    }
}

#[derive(Clone)]
struct SupportSources<'a> {
    segments: &'a [JkrNormalDissipationSegment],
    initial_motion_time_s: f64,
    origin_s: f64,
    partition: HeatPartition,
    cells: [CellId; 2],
    cursor: usize,
    work: WorkLedger,
    heat_enabled: bool,
    transfer_tolerance_j: f64,
    injected: [Dd; 2],
    segment_injected: [Dd; 2],
    absolute_transfer: Dd,
}
impl<'a> SupportSources<'a> {
    fn new(
        trace: &'a JkrNormalMotionTrace,
        origin_s: f64,
        route: &JkrViscousHeatRoute,
        cells: [CellId; 2],
        transfer_tolerance_j: f64,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, JkrSnapHeatError> {
        if !origin_s.is_finite() || !transfer_tolerance_j.is_finite() || transfer_tolerance_j < 0.0
        {
            return Err(invalid(
                "finite viscous-source clock and transfer tolerance required",
            ));
        }
        let mut work = WorkLedger::default();
        // This Step is held privately by the original trace, so its work is not
        // a mutable observation copied out of a caller's mechanical frame.
        work.record(DissipationStep::from_work(
            trace.step().damping_dissipation_j,
            route.partition,
        )?)?;
        let heat_enabled = work.surface_a_heat_j() > 0.0 || work.surface_b_heat_j() > 0.0;
        let initial_motion_time_s = trace.step().initial_state().time_s();
        let source = Self {
            segments: trace.segments(),
            initial_motion_time_s,
            origin_s,
            partition: route.partition,
            cells,
            cursor: 0,
            work,
            heat_enabled,
            transfer_tolerance_j,
            injected: [Dd::ZERO; 2],
            segment_injected: [Dd::ZERO; 2],
            absolute_transfer: Dd::ZERO,
        };
        let final_motion_time = trace.step().state().time_s();
        let final_time = finite(origin_s + final_motion_time)?;
        let mut previous_offset = 0.0;
        let mut previous_time = finite(origin_s + initial_motion_time_s)?;
        let mut previous_motion_time = initial_motion_time_s;
        let mut original_work = Dd::ZERO;
        for segment in source.segments {
            checkpoint()?;
            let start = segment.start_offset_s();
            let end = segment.end_offset_s();
            if !start.is_finite()
                || !end.is_finite()
                || start < previous_offset
                || end <= start
                || !segment.duration_s().is_finite()
                || segment.duration_s() <= 0.0
                || !segment.power_w().is_finite()
                || segment.power_w() <= 0.0
                || !segment.dissipated_work_j().is_finite()
                || segment.dissipated_work_j() <= 0.0
            {
                return Err(invalid(
                    "original viscous segments must be finite, positive and ordered",
                ));
            }
            original_work = original_work + Dd::from_f64(segment.dissipated_work_j());
            if heat_enabled {
                let (start_time, end_time) = source.times(segment)?;
                let end_motion = finite(initial_motion_time_s + end)?;
                if start_time < previous_time
                    || end_time <= start_time
                    || end_time > final_time
                    || end_motion > final_motion_time
                    || (start > previous_offset && start_time <= previous_time)
                {
                    return Err(invalid(
                        "viscous-source intervals must remain ordered and representable on the thermal clock",
                    ));
                }
                previous_time = end_time;
                previous_motion_time = end_motion;
            }
            previous_offset = end;
        }
        if heat_enabled && final_motion_time > previous_motion_time && final_time <= previous_time {
            return Err(invalid(
                "final viscous-source interval is not representable on the thermal clock",
            ));
        }
        let defect = finite((original_work - Dd::from_f64(work.dissipated_work_j())).to_f64())?;
        if defect.abs() > transfer_tolerance_j {
            return Err(JkrSnapHeatError::TransferBalance(defect));
        }
        Ok(source)
    }

    fn times(&self, segment: &JkrNormalDissipationSegment) -> Result<(f64, f64), JkrSnapHeatError> {
        Ok((
            finite(self.origin_s + finite(self.initial_motion_time_s + segment.start_offset_s())?)?,
            finite(self.origin_s + finite(self.initial_motion_time_s + segment.end_offset_s())?)?,
        ))
    }

    fn finish_segment(&mut self) -> Result<(), JkrSnapHeatError> {
        let segment = &self.segments[self.cursor];
        let expected = DissipationStep::from_work(segment.dissipated_work_j(), self.partition)?;
        let differences = [
            self.segment_injected[0] - Dd::from_f64(expected.surface_a_heat_j()),
            self.segment_injected[1] - Dd::from_f64(expected.surface_b_heat_j()),
        ];
        let partition_defect = Dd::from_f64(expected.surface_a_heat_j())
            + Dd::from_f64(expected.surface_b_heat_j())
            + Dd::from_f64(expected.other_work_j())
            - Dd::from_f64(expected.total_work_j());
        let absolute = differences[0].abs() + differences[1].abs() + partition_defect.abs();
        let defect = finite((differences[0] + differences[1] + partition_defect).to_f64())?;
        if defect.abs() > self.transfer_tolerance_j
            || finite(absolute.to_f64())? > self.transfer_tolerance_j
        {
            return Err(JkrSnapHeatError::TransferBalance(
                defect.abs().max(absolute.to_f64()),
            ));
        }
        self.absolute_transfer = self.absolute_transfer + absolute;
        self.segment_injected = [Dd::ZERO; 2];
        self.cursor += 1;
        Ok(())
    }

    fn finish(
        self,
        route: &JkrViscousHeatRoute,
    ) -> Result<(JkrViscousHeatReport, Dd), JkrSnapHeatError> {
        if self.heat_enabled && self.cursor != self.segments.len() {
            return Err(JkrSnapHeatError::HistoryMismatch);
        }
        let injected_surface_heat_j = [
            finite(self.injected[0].to_f64())?,
            finite(self.injected[1].to_f64())?,
        ];
        let residual = self.injected[0] + self.injected[1] + Dd::from_f64(self.work.other_work_j())
            - Dd::from_f64(self.work.dissipated_work_j());
        let transfer_residual_j = finite(residual.to_f64())?;
        if transfer_residual_j.abs() > self.transfer_tolerance_j {
            return Err(JkrSnapHeatError::TransferBalance(transfer_residual_j));
        }
        Ok((
            JkrViscousHeatReport {
                work: self.work,
                injected_surface_heat_j,
                outside_destination: route.outside_destination.clone(),
                transfer_residual_j,
            },
            self.absolute_transfer,
        ))
    }
}

impl DepositionObserver for SupportSources<'_> {
    type Error = JkrSnapHeatError;
    fn requires_step_interval(&self) -> bool {
        self.heat_enabled
    }
    fn next_transition_time_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, Self::Error> {
        checkpoint()?;
        if !self.heat_enabled {
            return Ok(None);
        }
        let Some(segment) = self.segments.get(self.cursor) else {
            return Ok(None);
        };
        let (start, end) = self.times(segment)?;
        if thermal.time_s() >= end {
            return Err(JkrSnapHeatError::HistoryMismatch);
        }
        Ok(Some(if thermal.time_s() < start { start } else { end }))
    }

    fn step_sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        start_s: f64,
        dt_s: f64,
        end_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        checkpoint()?;
        if !self.heat_enabled {
            return Ok(BTreeMap::new());
        }
        let Some(segment) = self.segments.get(self.cursor) else {
            return Ok(BTreeMap::new());
        };
        let (start, end) = self.times(segment)?;
        if start_s < start {
            return if end_s <= start {
                Ok(BTreeMap::new())
            } else {
                Err(JkrSnapHeatError::HistoryMismatch)
            };
        }
        if start_s >= end || end_s > end {
            return Err(JkrSnapHeatError::HistoryMismatch);
        }
        let powers = [
            finite(segment.power_w() * self.partition.surface_a())?,
            finite(segment.power_w() * self.partition.surface_b())?,
        ];
        let mut sources = BTreeMap::new();
        for (side, &power) in powers.iter().enumerate() {
            if power == 0.0 {
                continue;
            }
            let heat = finite(power * dt_s)?;
            self.injected[side] = self.injected[side] + Dd::from_f64(heat);
            self.segment_injected[side] = self.segment_injected[side] + Dd::from_f64(heat);
            sources.insert(self.cells[side], power);
        }
        Ok(sources)
    }

    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        checkpoint()?;
        if self.heat_enabled {
            if let Some(segment) = self.segments.get(self.cursor) {
                let (_, end) = self.times(segment)?;
                if thermal.time_s() == end {
                    self.finish_segment()?;
                }
            }
        }
        Ok(())
    }
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
