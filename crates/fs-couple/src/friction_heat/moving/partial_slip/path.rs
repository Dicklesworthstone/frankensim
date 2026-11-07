//! Accelerated, reversing motion through the original material brush and heat.
//!
//! A prescribed quadratic path is split at every relative-slip and material
//! residence reversal, then at the requested time resolution. Each interval
//! retains the actual quadratic swept areas; it is never replaced by the chord
//! joining its endpoints. The original transport/return-map and endpoint-source
//! split still requires time refinement. Each native regional heat increment
//! is held as constant mean power within that subinterval; spatially correct
//! swept area alone does not resolve the within-step plastic-power history.
//! This owner integrates no body forces.

use super::*;
use crate::friction_heat::moving::trajectory::{ContactCurves, QuadraticContactPath};

/// Time resolution and whole-path work/energy limits, in addition to the
/// existing per-interval constitutive, thermal and nonlinear controls.
#[derive(Clone, Copy, Debug)]
pub struct MovingPartialSlipPathControls {
    pub step: MovingPartialSlipControls,
    pub maximum_step_s: f64,
    pub maximum_steps: usize,
    pub maximum_trials: usize,
    pub total_energy_tolerance_j: f64,
}

/// One monotone relative-motion interval of the actual supplied trajectory.
/// Sources may read these endpoint velocities instead of inferring a constant
/// speed from total displacement. Thermal/source evaluation remains at the
/// proposed endpoint, not at an invented mean material temperature.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingPartialSlipPathPhase {
    pub index: usize,
    pub start_s: f64,
    pub end_s: f64,
    pub path: QuadraticContactPath,
    pub relative_velocity_m_s: [f64; 2],
}

pub struct MovingPartialSlipPathStep<P> {
    pub phase: MovingPartialSlipPathPhase,
    pub report: MovingPartialSlipReport<P>,
}

/// Accepted subintervals and their one complete mechanical/thermal account.
pub struct MovingPartialSlipPathReport<P> {
    pub steps: Vec<MovingPartialSlipPathStep<P>>,
    pub trials: usize,
    /// Integral of absolute relative speed, including reversals [m].
    pub relative_travel_m: f64,
    pub mechanical_work_j: f64,
    pub physical_heat_j: f64,
    pub numerical_dissipation_j: f64,
    pub elastic_outflow_j: f64,
    pub thermal_energy_change_j: f64,
    pub external_heat_j: f64,
    pub storage_before_j: f64,
    pub storage_after_j: f64,
    pub total_energy_residual_j: f64,
}

/// No prefix is published on a later source, budget, thermal or cancel failure.
pub struct MovingPartialSlipPathCandidate<'m, P> {
    contact: MovingPartialSlip,
    thermal: DepositionSimulation<'m>,
    report: MovingPartialSlipPathReport<P>,
}
impl<'m, P> MovingPartialSlipPathCandidate<'m, P> {
    pub const fn contact(&self) -> &MovingPartialSlip {
        &self.contact
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub const fn report(&self) -> &MovingPartialSlipPathReport<P> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        MovingPartialSlip,
        DepositionSimulation<'m>,
        MovingPartialSlipPathReport<P>,
    ) {
        (self.contact, self.thermal, self.report)
    }
}

fn position(
    curves: ContactCurves,
    path: QuadraticContactPath,
    t: f64,
) -> Result<MovingContactPosition, DepositionError> {
    // Retain the caller's exact accepted and terminal placements, avoiding a
    // reconstruction roundoff at the two ownership boundaries.
    if t == 0.0 {
        return Ok(path.start);
    }
    if t == 1.0 {
        return Ok(path.end);
    }
    Ok(MovingContactPosition {
        body_offsets_m: [curves.bodies[0].at(t)?, curves.bodies[1].at(t)?],
        footprint_offset_m: curves.footprint.at(t)?,
    })
}

impl MovingPartialSlip {
    /// Continue one complete prescribed quadratic path, including returns to
    /// its initial placement. Slip and footprint/material reversals are exact
    /// split events; the remaining monotone spans obey `maximum_step_s`.
    ///
    /// The current original pressure field and footprint remain fixed. Existing
    /// material strains, kinetic flags and finite-cell heat survive each step.
    /// Fresh material still enters unloaded; expelled memory is not resurrected
    /// when the footprint returns. Each nonlinear source callback is a pure
    /// trial; its payload is retained only with the accepted complete path.
    ///
    /// Actual quadratic swept residence determines local deposition and the
    /// time-mean shear moment. No numerical return loss or elastic outflow is
    /// deposited as heat. Unique subinterval work IDs derive from the supplied
    /// interval owner. Shared step/trial budgets and a final cumulative energy
    /// gate apply before any candidate can be returned.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_quadratic_candidate<'m, P, E>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        path: QuadraticContactPath,
        ownership: &GeneralizedWorkOwnership,
        controls: MovingPartialSlipPathControls,
        interface_temperature_guess_k: f64,
        resolve: impl Fn(
            &MovingPartialSlipPathPhase,
            &DepositionSimulation<'m>,
            &[MovingSlipRegion],
            &[SweptContactPatch],
            &[f64],
        ) -> Result<MovingPartialSlipSource<P>, E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MovingPartialSlipPathCandidate<'m, P>, MovingPartialSlipError<E>> {
        let duration = target_s - self.time_s;
        if cancelled() {
            return Err(DepositionError::Cancelled.into());
        }
        if path.start != self.position
            || thermal.time_s() != self.time_s
            || !duration.is_finite()
            || duration <= 0.0
            || !controls.maximum_step_s.is_finite()
            || controls.maximum_step_s <= 0.0
            || controls.maximum_steps == 0
            || controls.maximum_trials == 0
            || !controls.total_energy_tolerance_j.is_finite()
            || controls.total_energy_tolerance_j <= 0.0
            || ownership.patch_id() != self.reference.normal_patch().patch_id()
        {
            return Err(MovingPartialSlipError::Invalid(
                "synchronized quadratic path, original work owner and finite shared budgets required",
            ));
        }
        let curves = ContactCurves::quadratic(path).map_err(DepositionError::from)?;
        for side in 0..2 {
            let faces = &self.route.surfaces[side];
            if !curves
                .contained(
                    side,
                    self.route.footprint_m,
                    [
                        faces[0].face.span_m[0],
                        faces.last().unwrap().face.span_m[1],
                    ],
                )
                .map_err(DepositionError::from)?
            {
                return Err(MovingPartialSlipError::Invalid(
                    "complete quadratic footprint leaves the original material surface",
                ));
            }
        }
        let slip = curves.bodies[0]
            .subtract(curves.bodies[1])
            .map_err(DepositionError::from)?;
        let mut events = vec![0.0, 1.0];
        for curve in [
            slip,
            curves
                .footprint
                .subtract(curves.bodies[0])
                .map_err(DepositionError::from)?,
            curves
                .footprint
                .subtract(curves.bodies[1])
                .map_err(DepositionError::from)?,
        ] {
            if let Some(t) = curve.turning_point() {
                events.push(t);
            }
        }
        events.sort_by(f64::total_cmp);
        events.dedup();
        // Preflight the complete bounded schedule before cloning histories or
        // invoking any material source. No loop is driven by a rounded dt sum.
        let mut times = vec![0.0];
        for event in events.windows(2) {
            let count = ((event[1] - event[0]) * duration / controls.maximum_step_s).ceil();
            let remaining = controls.maximum_steps - (times.len() - 1);
            if !count.is_finite() || count > remaining as f64 {
                return Err(MovingPartialSlipError::StepBudget);
            }
            let count = (count as usize).max(1);
            if count > remaining {
                return Err(MovingPartialSlipError::StepBudget);
            }
            for i in 1..=count {
                if cancelled() {
                    return Err(DepositionError::Cancelled.into());
                }
                times.push(if i == count {
                    event[1]
                } else {
                    (event[1] - event[0]).mul_add(i as f64 / count as f64, event[0])
                });
            }
        }
        let storage_before_j = self.stored_energy_j()?;
        let mut contact = self.clone();
        let mut thermal = thermal.clone();
        let mut reports = Vec::with_capacity(times.len() - 1);
        let mut trials = 0usize;
        let mut guess_k = interface_temperature_guess_k;
        // Work, physical heat, numerical loss, outflow, thermal storage,
        // external heat, and unsigned slip. Compensated totals keep small
        // intervals visible beside a large accumulated account.
        let mut totals = [Dd::ZERO; 7];
        for (index, pair) in times.windows(2).enumerate() {
            if cancelled() {
                return Err(DepositionError::Cancelled.into());
            }
            let end_s = if pair[1] == 1.0 {
                target_s
            } else {
                finite(pair[1].mul_add(duration, self.time_s))?
            };
            if end_s <= contact.time_s {
                return Err(MovingPartialSlipError::Invalid(
                    "quadratic subinterval is not representable on the original clock",
                ));
            }
            let phase = MovingPartialSlipPathPhase {
                index,
                start_s: contact.time_s,
                end_s,
                path: QuadraticContactPath {
                    start: contact.position,
                    midpoint: position(curves, path, f64::midpoint(pair[0], pair[1]))?,
                    end: position(curves, path, pair[1])?,
                },
                relative_velocity_m_s: [
                    finite(slip.derivative(pair[0]).map_err(DepositionError::from)? / duration)?,
                    finite(slip.derivative(pair[1]).map_err(DepositionError::from)? / duration)?,
                ],
            };
            let remaining = controls.maximum_trials - trials;
            if remaining == 0 {
                return Err(MovingPartialSlipError::TrialBudget);
            }
            let mut step_controls = controls.step;
            step_controls.maximum_trials = step_controls.maximum_trials.min(remaining);
            let interval_owner = GeneralizedWorkOwnership::new(
                ownership.patch_id(),
                format!("{}/segment/{index}", ownership.interval_id()),
                ownership.longitudinal_coordinate_id(),
                ownership.lateral_coordinate_id(),
                ownership.torsional_coordinate_id(),
            )?;
            let proposal = contact.step_candidate_for_pressure(
                &thermal,
                end_s,
                phase.path.end,
                contact.pressure.as_ref(),
                Some(phase.path.midpoint),
                &interval_owner,
                step_controls,
                &vec![0.0; thermal.cells().len()],
                guess_k,
                |t, r, p, ti| resolve(&phase, t, r, p, ti),
                &mut cancelled,
            )?;
            let (next_contact, next_thermal, report) = proposal.into_parts();
            trials += report.trials;
            for (sum, value) in totals.iter_mut().zip([
                report.mechanical_work_j,
                report.heat.dissipated_j,
                report.numerical_dissipation_j,
                report.elastic_outflow_j,
                report.heat.thermal.stored_energy_change_j,
                report.heat.external_heat_j,
                report.relative_slip_m.abs(),
            ]) {
                *sum = *sum + Dd::from_f64(value);
            }
            // A numerical initial iterate only: every new interface retains
            // its own independently solved temperature coordinate.
            guess_k = report.heat.patches[0].interface_temperature_k;
            reports.push(MovingPartialSlipPathStep { phase, report });
            contact = next_contact;
            thermal = next_thermal;
        }
        let storage_after_j = contact.stored_energy_j()?;
        let total_energy_residual_j = finite(
            (Dd::from_f64(storage_after_j) - Dd::from_f64(storage_before_j)
                + totals[4]
                + totals[2]
                + totals[3]
                - totals[0]
                - totals[5])
                .to_f64(),
        )?;
        if total_energy_residual_j.abs() > controls.total_energy_tolerance_j {
            return Err(MovingPartialSlipError::Invalid(
                "complete moving path exceeds its cumulative energy budget",
            ));
        }
        if cancelled() {
            return Err(DepositionError::Cancelled.into());
        }
        let mut values = [0.0; 7];
        for (value, total) in values.iter_mut().zip(totals) {
            *value = finite(total.to_f64())?;
        }
        Ok(MovingPartialSlipPathCandidate {
            contact,
            thermal,
            report: MovingPartialSlipPathReport {
                steps: reports,
                trials,
                relative_travel_m: values[6],
                mechanical_work_j: values[0],
                physical_heat_j: values[1],
                numerical_dissipation_j: values[2],
                elastic_outflow_j: values[3],
                thermal_energy_change_j: values[4],
                external_heat_j: values[5],
                storage_before_j,
                storage_after_j,
                total_energy_residual_j,
            },
        })
    }
}
