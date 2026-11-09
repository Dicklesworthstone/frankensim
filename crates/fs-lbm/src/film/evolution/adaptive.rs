//! Adaptive scheduled evolution around the original conservative film owner.
//!
//! Compare one backward-Euler interval with its two actual representable
//! halves; publish only the native fine trajectory. Accuracy includes each
//! face's signed transported volume and each separate physical heat/work
//! channel, so constant height or counterflow cannot hide unresolved forcing.
//! The existing thermal controller supplies first-order error scaling; no
//! extrapolated film, replacement mobility or thermal-storage law is created.

use super::*;
use fs_conduction::enthalpy::conduction::implicit::adaptive::ConductionTimeOrder;

/// A query belongs to one authored continuous segment. Segment i ends at
/// events[i]; a query exactly at that event still uses segment i (left limit).
/// The next interval uses i+1, without perturbing an event's floating time.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThinFilmDriveQuery {
    pub time_s: f64,
    pub segment: usize,
}

/// Endpoint data held only over one original backward-Euler interval.
///
/// Newtonian source owners require None for temperatures. Carreau owners
/// require every original source temperature on EVERY query, including both
/// open ends. Density, surface tension and normal gravity remain unchanged.
/// These are prescribed temperatures, not evolved enthalpy or thermal feedback.
#[derive(Clone, Debug, PartialEq)]
pub struct ThinFilmDrive {
    pub acceleration_m_s2: f64,
    pub surface_shear_pa: Vec<f64>,
    pub open_ends: Option<ThinFilmOpenEnds>,
    pub face_temperatures_k: Option<Vec<f64>>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThinFilmAdaptiveControls {
    /// Bounds on proposed full spans; accepted native halves are shorter.
    /// Endpoint addition is rounded; a one-ULP tail closes at the exact event.
    /// Every solve, account and error estimate uses the actual representable spans.
    pub initial_dt_s: f64,
    pub minimum_dt_s: f64,
    pub maximum_dt_s: f64,
    /// Includes accuracy rejections and failed original nonlinear solves.
    pub maximum_trials: usize,
    pub maximum_accepted_steps: usize,
    pub height_absolute_tolerance_m: f64,
    /// Per original face, including each open half-cell.
    pub transport_absolute_tolerance_m2: f64,
    /// Per face AND per separate viscous/body/surface/reservoir channel.
    pub energy_absolute_tolerance_j_per_m: f64,
    /// Scales height CHANGES and individual interval transports/work, never
    /// cancelled net work, total liquid inventory or an arbitrary datum.
    pub relative_tolerance: f64,
}
impl Default for ThinFilmAdaptiveControls {
    fn default() -> Self {
        Self {
            initial_dt_s: 0.1,
            minimum_dt_s: 1e-8,
            maximum_dt_s: 1.0,
            maximum_trials: 65_536,
            maximum_accepted_steps: 131_072,
            height_absolute_tolerance_m: 1e-8,
            transport_absolute_tolerance_m2: 1e-10,
            energy_absolute_tolerance_j_per_m: 1e-8,
            relative_tolerance: 1e-3,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ThinFilmTimeError {
    pub maximum_height_difference_m: f64,
    pub maximum_transport_difference_m2: f64,
    pub maximum_energy_difference_j_per_m: f64,
    /// Estimated fine local error with the actual unequal-half correction.
    /// This is not a global, spatial, constitutive or model-error bound.
    pub normalized: f64,
}

/// Accepted time integral at an original profile source. Open ordering is
/// left end, interior links, right end; other topologies use interior links.
/// Positive volume is +x. The right reservoir's work includes its minus sign.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ThinFilmFaceIntegral {
    pub transported_volume_m2: f64,
    pub viscous_heat_j_per_m: f64,
    pub body_work_j_per_m: f64,
    pub surface_work_j_per_m: f64,
    pub boundary_work_j_per_m: f64,
}

#[derive(Debug)]
pub struct ThinFilmAdaptiveReport {
    pub time_s: f64,
    pub duration_s: f64,
    pub comparisons: usize,
    pub accepted_comparisons: usize,
    pub accuracy_rejections: usize,
    pub nonlinear_rejections: usize,
    pub events_reached: usize,
    pub accepted_steps: usize,
    pub maximum_accepted_normalized_error: f64,
    pub last_error: Option<ThinFilmTimeError>,
    /// Explicit continuation suggestion; no hidden PI state is retained.
    pub suggested_dt_s: f64,
    /// Original work includes coarse, rejected and failed native trials.
    pub endpoint_evaluations: usize,
    pub profile_evaluations: u64,
    pub face_integrals: Vec<ThinFilmFaceIntegral>,
    /// None only for an exactly zero-duration call.
    pub last_step: Option<ThinFilmStep>,
    pub volume_change_m2: f64,
    pub volume_balance_residual_m2: f64,
    pub stored_energy_change_j_per_m: f64,
    /// Sum of ONLY accepted native numerical defects, never physical heat.
    pub backward_euler_energy_defect_j_per_m: f64,
    pub viscous_heat_j_per_m: f64,
    pub body_work_j_per_m: f64,
    pub surface_work_j_per_m: f64,
    pub boundary_work_j_per_m: f64,
    pub energy_residual_j_per_m: f64,
}

#[derive(Debug)]
pub enum ThinFilmAdaptiveError<E> {
    Film(ThinFilmError),
    Source(E),
    Invalid(&'static str),
    TrialBudget,
    StepBudget,
    StepUnderflow {
        time_s: f64,
        span_s: f64,
    },
    TimeAccuracy {
        time_s: f64,
        span_s: f64,
        error: ThinFilmTimeError,
    },
}
impl<E> From<ThinFilmError> for ThinFilmAdaptiveError<E> {
    fn from(error: ThinFilmError) -> Self {
        Self::Film(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for ThinFilmAdaptiveError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive thin-film evolution refused: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for ThinFilmAdaptiveError<E> {}

impl ThinFilmLine {
    /// Original accepted Carreau temperature; Newtonian sources return None.
    pub fn face_temperature_k(&self, source: usize) -> Option<f64> {
        match &self.source {
            Source::Carreau(sources) => sources.get(source).map(|s| s.temperature_k()),
            Source::Newtonian(_) => None,
        }
    }

    /// Atomically prepare a prescribed temperature on every original Carreau
    /// source. Retains the exact law, support, density, tension and histories.
    /// No thermal storage, expansion or Marangoni coefficient is inferred.
    pub fn set_face_temperatures(
        &mut self,
        temperatures_k: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<(), ThinFilmError> {
        let check = RefCell::new(cancelled);
        let stopped = Cell::new(false);
        let checkpoint = || {
            if stopped.get() || (check.borrow_mut())() {
                stopped.set(true);
                Err(ThinFilmError::Cancelled)
            } else {
                Ok(())
            }
        };
        let source = self.rebound_source(temperatures_k, &checkpoint)?;
        checkpoint()?;
        self.source = source;
        Ok(())
    }

    fn rebound_source(
        &self,
        temperatures: &[f64],
        check: &dyn Fn() -> Result<(), ThinFilmError>,
    ) -> Result<Source, ThinFilmError> {
        check()?;
        let Source::Carreau(sources) = &self.source else {
            return Err(ThinFilmError::Invalid(
                "temperature rebinding requires original Carreau sources",
            ));
        };
        if temperatures.len() != sources.len() {
            return Err(ThinFilmError::Invalid(
                "one explicit temperature per original profile source required",
            ));
        }
        let mut prepared = Vec::with_capacity(sources.len());
        for (original, &temperature) in sources.iter().zip(temperatures) {
            check()?;
            prepared.push(
                if temperature.to_bits() == original.temperature_k().to_bits() {
                    *original
                } else {
                    LaminarCarreauSlit::new(
                        *original.law(),
                        temperature,
                        original.width_m(),
                        original.gap_m(),
                        original.length_m(),
                        original.maximum_laminar_shear_rate_per_s(),
                    )?
                },
            );
        }
        Ok(Source::Carreau(prepared))
    }

    fn scheduled_source(
        &self,
        drive: &ThinFilmDrive,
        check: &dyn Fn() -> Result<(), ThinFilmError>,
    ) -> Result<Source, ThinFilmError> {
        match (&self.source, &drive.face_temperatures_k) {
            (Source::Newtonian(material), None) => Ok(Source::Newtonian(*material)),
            (Source::Carreau(_), Some(temperatures)) => self.rebound_source(temperatures, check),
            _ => Err(ThinFilmError::Invalid(
                "scheduled Carreau sources require explicit temperatures; frozen Newtonian sources require none",
            )),
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn scheduled_interval(
        &mut self,
        end: f64,
        drive: &ThinFilmDrive,
        source: &Source,
        limits: ThinFilmLimits,
        controls: ThinFilmControls,
        work: &Work,
        check: &dyn Fn() -> Result<(), ThinFilmError>,
    ) -> Result<ThinFilmStep, ThinFilmError> {
        self.source = source.clone();
        let guess = self.heights.clone();
        self.advance_metered(
            end,
            drive.acceleration_m_s2,
            &drive.surface_shear_pa,
            drive.open_ends,
            limits,
            controls,
            &guess,
            work,
            || check().is_err(),
        )
    }

    /// Adapt to a target with explicit continuous drive segments.
    ///
    /// Events are finite, nonnegative, strictly increasing absolute times.
    /// Past/future entries may remain in the list for deterministic resumption.
    /// No comparison straddles an event. The resolver supplies the left limit
    /// at a segment's exact closing event; it must be repeatable and must not
    /// advance external physical history. A resolved endpoint source is shared
    /// by the coarse and fine endpoint, not queried separately for each solve.
    ///
    /// Only excessive time error and original Newton nonconvergence retry.
    /// Source, regime, conservation, cancellation and work refusals propagate.
    /// Native endpoint/profile ceilings apply to the ENTIRE adaptive call.
    /// Every accepted native half retains its original mass/energy gates; the
    /// final whole-interval account independently checks opening/final storage.
    /// On any refusal the complete accepted film, source temperatures and clock
    /// remain unchanged, even after earlier private fine pairs succeeded.
    ///
    /// A final/event-clipped comparison below the floor may be tried once,
    /// never accepted merely because it is small. Normal gravity, density and
    /// surface tension stay fixed. Fast momentum transients, contact lines,
    /// evolving thermal storage and constitutive uncertainty remain outside
    /// this quasi-steady model and its local time estimator.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_adaptive<E>(
        &mut self,
        target_s: f64,
        event_times_s: &[f64],
        limits: ThinFilmLimits,
        controls: ThinFilmControls,
        adaptive: ThinFilmAdaptiveControls,
        mut resolve: impl FnMut(ThinFilmDriveQuery) -> Result<ThinFilmDrive, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ThinFilmAdaptiveReport, ThinFilmAdaptiveError<E>> {
        let checkpoint = RefCell::new(cancelled);
        let stopped = Cell::new(false);
        let check = || {
            if stopped.get() || (checkpoint.borrow_mut())() {
                stopped.set(true);
                Err(ThinFilmError::Cancelled)
            } else {
                Ok(())
            }
        };
        check()?;
        admit(limits, controls)?;
        if !target_s.is_finite()
            || target_s < self.time
            || ![
                adaptive.initial_dt_s,
                adaptive.minimum_dt_s,
                adaptive.maximum_dt_s,
                adaptive.height_absolute_tolerance_m,
                adaptive.transport_absolute_tolerance_m2,
                adaptive.energy_absolute_tolerance_j_per_m,
            ]
            .into_iter()
            .all(positive)
            || adaptive.initial_dt_s < adaptive.minimum_dt_s
            || adaptive.maximum_dt_s < adaptive.minimum_dt_s
            || !(1..=1_000_000).contains(&adaptive.maximum_trials)
            || !(2..=2_000_000).contains(&adaptive.maximum_accepted_steps)
            || !adaptive.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&adaptive.relative_tolerance)
        {
            return Err(ThinFilmAdaptiveError::Invalid(
                "finite forward interval and bounded adaptive film controls required",
            ));
        }
        for (i, &event) in event_times_s.iter().enumerate() {
            if i % 128 == 0 {
                check()?;
            }
            if !nonnegative(event) || i > 0 && event <= event_times_s[i - 1] {
                return Err(ThinFilmAdaptiveError::Invalid(
                    "drive events must be finite, nonnegative and strictly increasing",
                ));
            }
        }
        self.admit_heights(&self.heights, limits, &check)?;
        let order = ConductionTimeOrder::First;
        let work = Work::default();
        let mut next = self.clone();
        let mut segment = event_times_s.partition_point(|&time| time <= self.time);
        let mut span = adaptive.initial_dt_s.min(adaptive.maximum_dt_s);
        let mut sums = vec![[Dd::from_f64(0.0); 5]; self.profile_source_count()];
        let mut defect = Dd::from_f64(0.0);
        let mut report = ThinFilmAdaptiveReport {
            time_s: target_s,
            duration_s: target_s - self.time,
            comparisons: 0,
            accepted_comparisons: 0,
            accuracy_rejections: 0,
            nonlinear_rejections: 0,
            events_reached: 0,
            accepted_steps: 0,
            maximum_accepted_normalized_error: 0.0,
            last_error: None,
            suggested_dt_s: span,
            endpoint_evaluations: 0,
            profile_evaluations: 0,
            face_integrals: Vec::new(),
            last_step: None,
            volume_change_m2: 0.0,
            volume_balance_residual_m2: 0.0,
            stored_energy_change_j_per_m: 0.0,
            backward_euler_energy_defect_j_per_m: 0.0,
            viscous_heat_j_per_m: 0.0,
            body_work_j_per_m: 0.0,
            surface_work_j_per_m: 0.0,
            boundary_work_j_per_m: 0.0,
            energy_residual_j_per_m: 0.0,
        };
        while next.time < target_s {
            check()?;
            if report.comparisons >= adaptive.maximum_trials {
                return Err(ThinFilmAdaptiveError::TrialBudget);
            }
            if adaptive.maximum_accepted_steps - report.accepted_steps < 2 {
                return Err(ThinFilmAdaptiveError::StepBudget);
            }
            let start = next.time;
            let event = event_times_s.get(segment).copied();
            let stop = event.map_or(target_s, |time| time.min(target_s));
            let proposal = span;
            let proposed_end = if span >= stop - start {
                stop
            } else {
                (start + span).min(stop)
            };
            // An exact decimal event can lie one ULP above the rounded sum.
            // Consume that representable endpoint in this actual interval,
            // rather than creating a spurious zero-interior tail comparison.
            let end = if proposed_end < stop && proposed_end.next_up() == stop {
                stop
            } else {
                proposed_end
            };
            let middle = f64::midpoint(start, end);
            let dt = end - start;
            if !(start < middle && middle < end && end <= stop) {
                return Err(ThinFilmAdaptiveError::StepUnderflow {
                    time_s: start,
                    span_s: dt,
                });
            }
            report.comparisons += 1;
            let first_drive = resolve(ThinFilmDriveQuery {
                time_s: middle,
                segment,
            })
            .map_err(ThinFilmAdaptiveError::Source)?;
            check()?;
            let end_drive = resolve(ThinFilmDriveQuery {
                time_s: end,
                segment,
            })
            .map_err(ThinFilmAdaptiveError::Source)?;
            check()?;
            let first_source = next.scheduled_source(&first_drive, &check)?;
            let end_source = next.scheduled_source(&end_drive, &check)?;
            let candidate = (|| -> Result<_, ThinFilmError> {
                let mut coarse = next.clone();
                let coarse_step = coarse.scheduled_interval(
                    end,
                    &end_drive,
                    &end_source,
                    limits,
                    controls,
                    &work,
                    &check,
                )?;
                let mut fine = next.clone();
                let first_step = fine.scheduled_interval(
                    middle,
                    &first_drive,
                    &first_source,
                    limits,
                    controls,
                    &work,
                    &check,
                )?;
                let last_step = fine.scheduled_interval(
                    end,
                    &end_drive,
                    &end_source,
                    limits,
                    controls,
                    &work,
                    &check,
                )?;
                Ok((coarse, fine, coarse_step, first_step, last_step))
            })();
            let (coarse, fine, coarse_step, first_step, last_step) = match candidate {
                Ok(value) => value,
                Err(error @ ThinFilmError::NotConverged(_)) => {
                    report.nonlinear_rejections += 1;
                    if dt <= adaptive.minimum_dt_s {
                        return Err(error.into());
                    }
                    span = (0.5 * dt).max(adaptive.minimum_dt_s);
                    if span >= dt {
                        return Err(error.into());
                    }
                    continue;
                }
                Err(error) => return Err(error.into()),
            };
            let coarse_integrals = interval_integrals(&coarse_step, &end_drive, self)?;
            let first_integrals = interval_integrals(&first_step, &first_drive, self)?;
            let last_integrals = interval_integrals(&last_step, &end_drive, self)?;
            let correction = order
                .fine_error_correction(middle - start, end - middle)
                .ok_or(ThinFilmAdaptiveError::Invalid(
                    "unrepresentable fine time-error correction",
                ))?;
            let error = comparison_error(
                &next.heights,
                &coarse.heights,
                &fine.heights,
                &coarse_integrals,
                &first_integrals,
                &last_integrals,
                correction,
                adaptive,
                &check,
            )?;
            report.last_error = Some(error);
            if error.normalized > 1.0 {
                report.accuracy_rejections += 1;
                if dt <= adaptive.minimum_dt_s {
                    return Err(ThinFilmAdaptiveError::TimeAccuracy {
                        time_s: start,
                        span_s: dt,
                        error,
                    });
                }
                let factor = order
                    .rejected_reduction_factor(error.normalized)
                    .ok_or(ThinFilmAdaptiveError::Invalid("invalid film error ratio"))?;
                span = (dt * factor).max(adaptive.minimum_dt_s);
                if span >= dt {
                    return Err(ThinFilmAdaptiveError::TimeAccuracy {
                        time_s: start,
                        span_s: dt,
                        error,
                    });
                }
                continue;
            }
            check()?;
            for (i, sum) in sums.iter_mut().enumerate() {
                if i % 128 == 0 {
                    check()?;
                }
                for k in 0..5 {
                    sum[k] = sum[k]
                        + Dd::from_f64(first_integrals[i][k])
                        + Dd::from_f64(last_integrals[i][k]);
                }
            }
            defect = defect
                + Dd::from_f64(first_step.backward_euler_energy_defect_j_per_m)
                + Dd::from_f64(last_step.backward_euler_energy_defect_j_per_m);
            report.accepted_comparisons += 1;
            report.accepted_steps += 2;
            report.maximum_accepted_normalized_error = report
                .maximum_accepted_normalized_error
                .max(error.normalized);
            report.last_step = Some(last_step);
            next = fine;
            let factor = order.accepted_growth_factor(error.normalized).ok_or(
                ThinFilmAdaptiveError::Invalid("invalid accepted film error ratio"),
            )?;
            span = (dt * factor)
                .max(adaptive.minimum_dt_s)
                .min(adaptive.maximum_dt_s);
            if event == Some(end) {
                // Restore the unclipped proposal at decimal event boundaries;
                // a tiny final segment must not become the next nominal step.
                span = span.max(proposal).min(adaptive.maximum_dt_s);
                segment += 1;
                report.events_reached += 1;
            }
        }
        let mut totals = [Dd::from_f64(0.0); 5];
        for (i, sum) in sums.iter().enumerate() {
            if i % 128 == 0 {
                check()?;
            }
            let values: [f64; 5] = core::array::from_fn(|k| sum[k].to_f64());
            for k in 0..5 {
                finite(values[k])?;
                totals[k] = totals[k] + sum[k];
            }
            report.face_integrals.push(ThinFilmFaceIntegral {
                transported_volume_m2: values[0],
                viscous_heat_j_per_m: values[1],
                body_work_j_per_m: values[2],
                surface_work_j_per_m: values[3],
                boundary_work_j_per_m: values[4],
            });
        }
        let delta: Vec<_> = next
            .heights
            .iter()
            .zip(&self.heights)
            .map(|(a, b)| a - b)
            .collect();
        report.volume_change_m2 = total(delta.iter().map(|h| h * self.geometry.spacing_m))?;
        report.volume_balance_residual_m2 = if self.geometry.boundary == ThinFilmBoundary::Open {
            total([
                report.volume_change_m2,
                -report.face_integrals[0].transported_volume_m2,
                report.face_integrals.last().unwrap().transported_volume_m2,
            ])?
        } else {
            report.volume_change_m2
        };
        report.stored_energy_change_j_per_m = self
            .energy_terms(&self.heights, &next.heights, &delta, &check)?
            .0;
        report.backward_euler_energy_defect_j_per_m = finite(defect.to_f64())?;
        report.viscous_heat_j_per_m = finite(totals[1].to_f64())?;
        report.body_work_j_per_m = finite(totals[2].to_f64())?;
        report.surface_work_j_per_m = finite(totals[3].to_f64())?;
        report.boundary_work_j_per_m = finite(totals[4].to_f64())?;
        report.energy_residual_j_per_m = total([
            report.stored_energy_change_j_per_m,
            report.backward_euler_energy_defect_j_per_m,
            report.viscous_heat_j_per_m,
            -report.body_work_j_per_m,
            -report.surface_work_j_per_m,
            -report.boundary_work_j_per_m,
        ])?;
        if report.volume_balance_residual_m2.abs() > controls.volume_tolerance_m2 {
            return Err(ThinFilmError::VolumeResidual {
                residual_m2: report.volume_balance_residual_m2,
                tolerance_m2: controls.volume_tolerance_m2,
            }
            .into());
        }
        if report.energy_residual_j_per_m.abs() > controls.energy_tolerance_j_per_m {
            return Err(ThinFilmError::EnergyResidual {
                residual_j_per_m: report.energy_residual_j_per_m,
                tolerance_j_per_m: controls.energy_tolerance_j_per_m,
            }
            .into());
        }
        report.suggested_dt_s = span;
        report.endpoint_evaluations = work.endpoints.get();
        report.profile_evaluations = work.profiles.get();
        check()?;
        *self = next;
        Ok(report)
    }
}

// Original accepted profiles define EVERY integral. Channels are compared
// separately before any sums, including signed counterflow and reservoir work.
fn interval_integrals(
    step: &ThinFilmStep,
    drive: &ThinFilmDrive,
    line: &ThinFilmLine,
) -> Result<Vec<[f64; 5]>, ThinFilmError> {
    let mut result = Vec::with_capacity(line.profile_source_count());
    let dt = step.duration_s;
    let dx = line.geometry.spacing_m;
    let body = finite(line.density * drive.acceleration_m_s2)?;
    let mut append = |profile: &ThinFilmProfile, length: f64, reservoir_pressure: f64| {
        let values = [
            dt * profile.flow_per_width_m2_s,
            dt * length * profile.viscous_dissipation_w_m2,
            dt * length * body * profile.flow_per_width_m2_s,
            dt * length * profile.surface_shear_pa * profile.surface_velocity_m_s,
            dt * reservoir_pressure * profile.flow_per_width_m2_s,
        ];
        for value in values {
            finite(value)?;
        }
        result.push(values);
        Ok::<_, ThinFilmError>(())
    };
    if let Some(open) = &step.open_ends {
        append(&open.left, 0.5 * dx, open.ends.left.pressure_pa)?;
    }
    for face in &step.faces {
        append(&face.profile, dx, 0.0)?;
    }
    if let Some(open) = &step.open_ends {
        append(&open.right, 0.5 * dx, -open.ends.right.pressure_pa)?;
    }
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
fn comparison_error(
    opening: &[f64],
    coarse: &[f64],
    fine: &[f64],
    coarse_flux: &[[f64; 5]],
    first_flux: &[[f64; 5]],
    last_flux: &[[f64; 5]],
    correction: f64,
    controls: ThinFilmAdaptiveControls,
    check: &dyn Fn() -> Result<(), ThinFilmError>,
) -> Result<ThinFilmTimeError, ThinFilmError> {
    let mut error = ThinFilmTimeError::default();
    let normalized = |a: f64, b: f64, absolute: f64| -> Result<f64, ThinFilmError> {
        let scale = finite(absolute + controls.relative_tolerance * a.abs().max(b.abs()))?;
        Ok(correction * ((a - b).abs() / scale))
    };
    for (i, ((&old, &a), &b)) in opening.iter().zip(coarse).zip(fine).enumerate() {
        if i % 128 == 0 {
            check()?;
        }
        error.maximum_height_difference_m = error.maximum_height_difference_m.max((a - b).abs());
        error.normalized = error.normalized.max(normalized(
            a - old,
            b - old,
            controls.height_absolute_tolerance_m,
        )?);
    }
    for i in 0..coarse_flux.len() {
        if i % 128 == 0 {
            check()?;
        }
        for k in 0..5 {
            let a = coarse_flux[i][k];
            let b = total([first_flux[i][k], last_flux[i][k]])?;
            let difference = (a - b).abs();
            let absolute = if k == 0 {
                error.maximum_transport_difference_m2 =
                    error.maximum_transport_difference_m2.max(difference);
                controls.transport_absolute_tolerance_m2
            } else {
                error.maximum_energy_difference_j_per_m =
                    error.maximum_energy_difference_j_per_m.max(difference);
                controls.energy_absolute_tolerance_j_per_m
            };
            error.normalized = error.normalized.max(normalized(a, b, absolute)?);
        }
    }
    Ok(error)
}
