//! Prescribed combined loading and original finite thermal transport.
//!
//! The sealed source ramp owns its affine physical prescription and first
//! boundary time. Native thermal events at that time run before an explicitly
//! admitted fold release. Separation removes only the owned contact patch,
//! restoring the original ambient exposure before subsequent transport.
//!
//! Attached coverage and conductance use a bounded left-endpoint footprint
//! split. Original source samples update the actual patch at interval ends;
//! they do not advance accepted mechanical history. This is a first-order
//! footprint approximation requiring timestep refinement, not an exact swept
//! area integral. Constitutive properties and conductance per area stay frozen.
//! No inertial motion, recoil, rebonding or temperature feedback is inferred.

/// Original thermal contact ownership, evolving coverage and separation.
pub mod footprint;

use self::footprint::JkrCombinedHeatFootprint;
use super::{
    JkrCombinedHeat, JkrCombinedHeatControls, JkrCombinedHeatError, JkrCombinedHeatReport,
    JkrCombinedHeatSegment, JkrCombinedHeatStamp, JkrSnapHeatControls, JkrSnapHeatError, dd,
    finite,
};
use crate::friction_heat::adhesive::snap::timed::AuditSum;
use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::mixed_mode::combined::{
    JkrCombinedBoundary,
    loading::{JkrCombinedLoadingError, timed::JkrCombinedLoadingRamp},
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation, EnergyAudit, StepConfig};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum JkrCombinedBoundaryPolicy {
    /// Accept only the source's first-boundary prefix, leaving it attached.
    Stop,
    /// Release a new separating radial-fold arrival, then traverse the remaining
    /// prescription detached. Other boundaries retain the stopped prefix.
    ReleaseIfSeparatingFold,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrTimedCombinedHeatControls {
    /// Aggregate original mechanical-stage count and per-impulse heat budgets.
    pub heat: JkrSnapHeatControls,
    /// Original native explicit transport controls; its CFL remains authoritative.
    pub transport: StepConfig,
    /// Shared allowance across every native transport interval in this proposal.
    pub maximum_transport_substeps: usize,
    /// Bound original source area samples and footprint updates independently.
    pub maximum_footprint_intervals: usize,
    /// Positive maximum duration over which an attached footprint is held fixed.
    pub maximum_footprint_dt_s: f64,
    /// Whole-window mechanical, partition and native thermal defect allowance.
    pub total_energy_tolerance_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct JkrTimedCombinedHeatReport {
    pub initial_time_s: f64,
    pub requested_final_time_s: f64,
    /// May precede the request when a boundary remains attached.
    pub final_time_s: f64,
    pub untraversed_time_s: f64,
    pub boundary_time_s: Option<f64>,
    pub release_time_s: Option<f64>,
    pub clock_residual_s: f64,
    /// Reversals split at the represented zero-load time, including the actual
    /// maximum source footprint rather than only the two smaller endpoints.
    pub zero_crossing_time_s: Option<f64>,
    pub zero_crossing_clock_residual_s: f64,
    pub footprint_intervals: usize,
    pub maximum_footprint_interval_s: f64,
    pub maximum_admitted_area_m2: f64,
    /// Complete original load, optional explicit release and detached-tail reports.
    pub mechanics: Vec<JkrCombinedHeatReport>,
    /// Ordinary native transport/events, excluding same-time contact heat.
    pub transport: EnergyAudit,
    /// All native transport, births and original contact heat exactly once.
    pub thermal: EnergyAudit,
    /// Delta U_contact + Delta H + outside work - prescribed mechanical work
    /// - ordinary thermal supply - incoming material energy. Relative native
    /// audits retain small contact heat beside large formation enthalpies.
    pub combined_energy_residual_j: f64,
    pub absolute_stage_energy_residual_j: f64,
}

pub struct JkrTimedCombinedHeatCandidate<'s, 'm> {
    history: JkrCombinedHeat<'s>,
    thermal: DepositionSimulation<'m>,
    footprint: JkrCombinedHeatFootprint,
    report: JkrTimedCombinedHeatReport,
}
impl<'s, 'm> JkrTimedCombinedHeatCandidate<'s, 'm> {
    pub const fn history(&self) -> &JkrCombinedHeat<'s> {
        &self.history
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub const fn footprint(&self) -> &JkrCombinedHeatFootprint {
        &self.footprint
    }
    pub const fn report(&self) -> &JkrTimedCombinedHeatReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        JkrCombinedHeat<'s>,
        DepositionSimulation<'m>,
        JkrCombinedHeatFootprint,
        JkrTimedCombinedHeatReport,
    ) {
        (self.history, self.thermal, self.footprint, self.report)
    }
}

impl<'s> JkrCombinedHeat<'s> {
    /// Stage one source-authenticated physical ramp and finite-cell transport.
    ///
    /// Incoming source history, thermal time and last accepted heat time must
    /// match exactly. Original native events execute first at each timestamp;
    /// an event that changes this owner's patch or receivers refuses the whole
    /// candidate. Other contacts, boundaries, reservoirs and births retain their
    /// original thermal ownership and order.
    ///
    /// The whole source radius envelope is reserved against both actual cell
    /// faces and competing coverage, including between-sample reversal maxima.
    /// Pure geometry samples never publish additional mechanical increments.
    /// A release, its heat, patch removal and any detached tail accept together.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn propose_timed_ramp<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        expected: JkrCombinedHeatStamp,
        ramp: &JkrCombinedLoadingRamp,
        footprint: &JkrCombinedHeatFootprint,
        policy: JkrCombinedBoundaryPolicy,
        controls: JkrTimedCombinedHeatControls,
        gate: &CancelGate,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<JkrTimedCombinedHeatCandidate<'s, 'm>, JkrCombinedHeatError> {
        let mut check = || {
            if gate.is_requested() {
                return Err(DepositionError::Cancelled);
            }
            checkpoint()?;
            if gate.is_requested() {
                return Err(DepositionError::Cancelled);
            }
            Ok(())
        };
        check()?;
        if expected != self.stamp
            || ramp.previous_state() != self.state()
            || ramp.initial_time_s() != thermal.time_s()
            || thermal.time_s() != self.accepted_time_s
        {
            return Err(JkrSnapHeatError::HistoryMismatch.into());
        }
        validate_controls(controls)?;
        footprint.admit(self, thermal)?;
        let release = policy == JkrCombinedBoundaryPolicy::ReleaseIfSeparatingFold
            && ramp.boundary_time_s().is_some()
            && ramp.state().contact().is_some_and(|point| {
                point.boundary == Some(JkrCombinedBoundary::RadialFold)
                    && point.coordinates[0] < 0.0
            });
        let tail = release && ramp.accepted_final_time_s() < ramp.requested_final_time_s();
        let stages = 1 + usize::from(release) + usize::from(tail);
        if stages > controls.heat.maximum_segments {
            return Err(DepositionError::Budget("combined timed mechanical segments").into());
        }
        let maximum_area_m2 = if self.state().contact().is_some() {
            let radius = ramp.step().report().maximum_contact_radius_m;
            finite(PI * radius * radius)?
        } else {
            0.0
        };
        let zero = zero_crossing(ramp)?;
        let heat_controls = JkrCombinedHeatControls {
            heat: controls.heat,
            total_energy_tolerance_j: controls.total_energy_tolerance_j,
        };
        let mut accounts = Accounts::default();
        let mut history = self.clone();
        let mut footprint = footprint.clone();
        // Flush due native events even when this source stops at its initial
        // boundary. Patch changes always occur after that native event cursor.
        let mut current = march(
            &history,
            thermal,
            &footprint,
            ramp.initial_time_s(),
            maximum_area_m2,
            controls,
            &mut accounts,
            gate,
            &mut check,
        )?;
        footprint.reserve_maximum_area(&current, maximum_area_m2, &mut check)?;
        let mut footprint_intervals = 0usize;
        let mut maximum_footprint_interval_s: f64 = 0.0;
        if history.state().contact().is_some() {
            while current.time_s() < ramp.accepted_final_time_s() {
                check()?;
                if footprint_intervals >= controls.maximum_footprint_intervals {
                    return Err(DepositionError::Budget("combined footprint intervals").into());
                }
                let start = current.time_s();
                let end = next_footprint_time(
                    start,
                    ramp.accepted_final_time_s(),
                    controls.maximum_footprint_dt_s,
                    zero.map(|crossing| crossing.0),
                )?;
                let point = if end == ramp.accepted_final_time_s() {
                    ramp.state().contact()
                } else {
                    let load = if zero.is_some_and(|crossing| end == crossing.0) {
                        0.0
                    } else {
                        ramp.load_factor_at_time(end)?
                    };
                    history.loading.initial_state(load, gate)?.contact()
                }
                .ok_or_else(|| invalid("attached footprint sample needs original contact"))?;
                current = march(
                    &history,
                    &current,
                    &footprint,
                    end,
                    maximum_area_m2,
                    controls,
                    &mut accounts,
                    gate,
                    &mut check,
                )?;
                if footprint.area_m2() != Some(point.contact_area_m2) {
                    (current, footprint) =
                        footprint.set_area(&current, point.contact_area_m2, &mut check)?;
                }
                footprint_intervals += 1;
                maximum_footprint_interval_s = maximum_footprint_interval_s.max(end - start);
            }
        } else {
            current = march(
                &history,
                &current,
                &footprint,
                ramp.accepted_final_time_s(),
                0.0,
                controls,
                &mut accounts,
                gate,
                &mut check,
            )?;
        }
        let mut mechanics = Vec::with_capacity(stages);
        let (next_history, next_thermal, report) = history
            .propose(
                &current,
                history.stamp(),
                &[JkrCombinedHeatSegment::LoadFactor(
                    ramp.step().report().requested_load_factor,
                )],
                heat_controls,
                gate,
                &mut check,
            )?
            .into_parts();
        if next_history.state() != ramp.state()
            || report.segments.as_slice() != [*ramp.step().report()]
        {
            return Err(JkrSnapHeatError::HistoryMismatch.into());
        }
        accounts.add_heat(&report)?;
        mechanics.push(report);
        history = next_history;
        current = next_thermal;
        footprint.admit(&history, &current)?;
        let release_time_s = if release {
            let (next_history, next_thermal, report) = history
                .propose(
                    &current,
                    history.stamp(),
                    &[JkrCombinedHeatSegment::DetachAtFold],
                    heat_controls,
                    gate,
                    &mut check,
                )?
                .into_parts();
            accounts.add_heat(&report)?;
            mechanics.push(report);
            history = next_history;
            (current, footprint) = footprint.detach(&next_thermal, &mut check)?;
            Some(current.time_s())
        } else {
            None
        };
        if tail {
            current = march(
                &history,
                &current,
                &footprint,
                ramp.requested_final_time_s(),
                0.0,
                controls,
                &mut accounts,
                gate,
                &mut check,
            )?;
            let (next_history, next_thermal, report) = history
                .propose(
                    &current,
                    history.stamp(),
                    &[JkrCombinedHeatSegment::LoadFactor(
                        ramp.step().report().requested_load_factor,
                    )],
                    heat_controls,
                    gate,
                    &mut check,
                )?
                .into_parts();
            accounts.add_heat(&report)?;
            mechanics.push(report);
            history = next_history;
            current = next_thermal;
        }
        footprint.admit(&history, &current)?;
        let combined_energy_residual_j = finite(accounts.combined.to_f64())?;
        let absolute_stage_energy_residual_j = finite(accounts.absolute.to_f64())?;
        if combined_energy_residual_j.abs() > controls.total_energy_tolerance_j
            || absolute_stage_energy_residual_j > controls.total_energy_tolerance_j
        {
            return Err(JkrCombinedHeatError::EnergyBalance {
                residual_j: combined_energy_residual_j,
                absolute_stage_residual_j: absolute_stage_energy_residual_j,
            });
        }
        let report = JkrTimedCombinedHeatReport {
            initial_time_s: ramp.initial_time_s(),
            requested_final_time_s: ramp.requested_final_time_s(),
            final_time_s: current.time_s(),
            untraversed_time_s: finite(
                (dd(ramp.requested_final_time_s()) - dd(current.time_s())).to_f64(),
            )?,
            boundary_time_s: ramp.boundary_time_s(),
            release_time_s,
            clock_residual_s: ramp.clock_residual_s(),
            zero_crossing_time_s: zero.map(|crossing| crossing.0),
            zero_crossing_clock_residual_s: zero.map_or(0.0, |crossing| crossing.1),
            footprint_intervals,
            maximum_footprint_interval_s,
            maximum_admitted_area_m2: maximum_area_m2,
            mechanics,
            transport: accounts.transport.audit()?,
            thermal: accounts.thermal.audit()?,
            combined_energy_residual_j,
            absolute_stage_energy_residual_j,
        };
        check()?;
        Ok(JkrTimedCombinedHeatCandidate {
            history,
            thermal: current,
            footprint,
            report,
        })
    }

    /// Publish loading history, thermal time/enthalpy and the owned footprint
    /// only after the entire prescribed interval and final checkpoint accept.
    #[allow(clippy::too_many_arguments)]
    pub fn apply_timed_ramp<'m>(
        &mut self,
        thermal: &mut DepositionSimulation<'m>,
        expected: JkrCombinedHeatStamp,
        ramp: &JkrCombinedLoadingRamp,
        footprint: &mut JkrCombinedHeatFootprint,
        policy: JkrCombinedBoundaryPolicy,
        controls: JkrTimedCombinedHeatControls,
        gate: &CancelGate,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<JkrTimedCombinedHeatReport, JkrCombinedHeatError> {
        let (history, next, next_footprint, report) = self
            .propose_timed_ramp(
                thermal,
                expected,
                ramp,
                footprint,
                policy,
                controls,
                gate,
                &mut checkpoint,
            )?
            .into_parts();
        checkpoint()?;
        if gate.is_requested() {
            return Err(JkrCombinedLoadingError::Cancelled.into());
        }
        *self = history;
        *thermal = next;
        *footprint = next_footprint;
        Ok(report)
    }
}

#[derive(Default)]
struct Accounts {
    transport: AuditSum,
    thermal: AuditSum,
    combined: Dd,
    absolute: Dd,
}
impl Accounts {
    fn add_transport(&mut self, audit: EnergyAudit) -> Result<(), JkrCombinedHeatError> {
        self.transport.add(audit)?;
        self.thermal.add(audit)?;
        self.combined = self.combined + dd(audit.residual_j);
        self.absolute = self.absolute + dd(audit.residual_j.abs());
        Ok(())
    }

    fn add_heat(&mut self, report: &JkrCombinedHeatReport) -> Result<(), JkrCombinedHeatError> {
        self.thermal.add(report.thermal)?;
        self.combined = self.combined + dd(report.combined_energy_residual_j);
        self.absolute = self.absolute + dd(report.absolute_stage_energy_residual_j);
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn march<'m>(
    history: &JkrCombinedHeat<'_>,
    thermal: &DepositionSimulation<'m>,
    footprint: &JkrCombinedHeatFootprint,
    target_s: f64,
    maximum_area_m2: f64,
    controls: JkrTimedCombinedHeatControls,
    accounts: &mut Accounts,
    gate: &CancelGate,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<DepositionSimulation<'m>, JkrCombinedHeatError> {
    checkpoint()?;
    let remaining = controls
        .maximum_transport_substeps
        .checked_sub(accounts.transport.audit()?.substeps)
        .ok_or(DepositionError::Budget(
            "combined timed total thermal substeps",
        ))?;
    let guard = footprint.guard(history, maximum_area_m2, gate);
    let (thermal, _, audit) = thermal.advance_observed_limited(
        &guard,
        target_s,
        controls.transport,
        remaining,
        &mut *checkpoint,
    )?;
    footprint.admit_transport(history, &thermal)?;
    accounts.add_transport(audit)?;
    checkpoint()?;
    Ok(thermal)
}

fn validate_controls(controls: JkrTimedCombinedHeatControls) -> Result<(), JkrCombinedHeatError> {
    if !controls.maximum_footprint_dt_s.is_finite()
        || controls.maximum_footprint_dt_s <= 0.0
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
        .any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err(invalid(
            "positive footprint/thermal intervals and finite energy budgets required",
        ));
    }
    Ok(())
}

fn next_footprint_time(
    start: f64,
    end: f64,
    maximum_dt_s: f64,
    zero: Option<f64>,
) -> Result<f64, JkrCombinedHeatError> {
    let mut next = if end - start <= maximum_dt_s {
        end
    } else {
        finite(start + maximum_dt_s)?.min(end)
    };
    // An internal footprint cut may round down; a physical source/event time
    // must never be moved. At a large clock origin, rounded addition can make
    // the represented elapsed interval longer than its declared maximum.
    if next - start > maximum_dt_s {
        next = next.next_down();
    }
    if let Some(zero) = zero.filter(|time| *time > start && *time < next) {
        next = zero;
    }
    if next <= start || next - start > maximum_dt_s {
        return Err(invalid(
            "unrepresentable positive combined footprint interval",
        ));
    }
    Ok(next)
}

fn zero_crossing(
    ramp: &JkrCombinedLoadingRamp,
) -> Result<Option<(f64, f64)>, JkrCombinedHeatError> {
    let initial = ramp.previous_state().load_factor();
    let accepted = ramp.state().load_factor();
    if ramp.previous_state().contact().is_none()
        || !((initial > 0.0 && accepted < 0.0) || (initial < 0.0 && accepted > 0.0))
    {
        return Ok(None);
    }
    let fraction = -dd(initial) / (dd(ramp.step().report().requested_load_factor) - dd(initial));
    let exact = dd(ramp.initial_time_s())
        + fraction * (dd(ramp.requested_final_time_s()) - dd(ramp.initial_time_s()));
    let time = finite(exact.to_f64())?;
    if time <= ramp.initial_time_s() || time >= ramp.accepted_final_time_s() {
        return Err(invalid(
            "interior zero-load crossing collapses onto a ramp endpoint",
        ));
    }
    let residual = dd(time) - exact;
    if dd(ramp.maximum_clock_error_s()).lt(residual.abs()) {
        return Err(JkrCombinedLoadingError::ClockRoundoff {
            residual_s: residual.to_f64(),
            tolerance_s: ramp.maximum_clock_error_s(),
        }
        .into());
    }
    Ok(Some((time, finite(residual.to_f64())?)))
}

fn invalid(message: &'static str) -> JkrCombinedHeatError {
    JkrSnapHeatError::Invalid(message).into()
}

#[cfg(test)]
mod tests {
    use super::next_footprint_time;

    #[test]
    fn internal_footprint_cuts_respect_duration_at_large_clock_origins() {
        let start = 9_007_199_254_740_992.0;
        let end = start + 8.0;
        // The nominal 3-second sum rounds upward to a 4-second interval.
        let next = next_footprint_time(start, end, 3.0, None).unwrap();
        assert_eq!(next, start + 2.0);
        assert!(next - start <= 3.0);
        assert!(next_footprint_time(start, end, 1.0, None).is_err());
        // The actual source endpoint is retained when it is reached later.
        assert_eq!(
            next_footprint_time(start + 6.0, end, 3.0, None).unwrap(),
            end
        );
    }

    #[test]
    fn physical_zero_load_and_source_endpoint_times_are_not_shifted() {
        let zero = 0.35;
        assert_eq!(
            next_footprint_time(0.2, 0.6, 0.2, Some(zero)).unwrap(),
            zero
        );
        assert_eq!(next_footprint_time(0.4, 0.6, 0.2, Some(zero)).unwrap(), 0.6);
    }
}
