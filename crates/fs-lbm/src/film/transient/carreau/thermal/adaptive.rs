//! Adaptive scheduled evolution of the original coupled caloric/film owners.
//!
//! The accepted trajectory is always two genuine native intervals. The
//! endpoint-temperature/midpoint-momentum composition remains first order;
//! no extrapolated velocity, enthalpy, viscosity or heat is published.

use super::*;
use fs_conduction::enthalpy::StepConfig;
use fs_conduction::enthalpy::conduction::implicit::adaptive::ConductionTimeOrder;

/// One actual native interval within an authored continuous source segment.
/// A closing event belongs to this segment's left limit. Drive queries use
/// midpoint_s; heat queries see the actual endpoint candidate and endpoint_s.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalFilmQuery {
    pub start_s: f64,
    pub midpoint_s: f64,
    pub endpoint_s: f64,
    pub segment: usize,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalFilmAdaptiveControls {
    /// Full comparison spans; accepted native halves are shorter.
    pub initial_dt_s: f64,
    pub minimum_dt_s: f64,
    pub maximum_dt_s: f64,
    pub maximum_trials: usize,
    /// Further capped by the original DepositionBudget::max_substeps for the
    /// complete window; individual native half-steps do not renew that limit.
    pub maximum_accepted_steps: usize,
    pub velocity_absolute_tolerance_m_s: f64,
    pub temperature_absolute_tolerance_k: f64,
    /// Per original cell and separate physical heat/work channel.
    pub energy_absolute_tolerance_j: f64,
    pub transport_absolute_tolerance_m2: f64,
    /// Scales changes and interval exchanges, never absolute enthalpy or
    /// temperature datums, total stored energy, or cancelled net actuator work.
    pub relative_tolerance: f64,
    /// Independent full-window storage/work and heat-partition gates. This
    /// does not replace any individual native algebraic acceptance tolerance.
    pub window_energy_tolerance_j: f64,
}
impl Default for ThermalFilmAdaptiveControls {
    fn default() -> Self {
        Self {
            initial_dt_s: 0.01,
            minimum_dt_s: 1e-8,
            maximum_dt_s: 1.0,
            maximum_trials: 65_536,
            maximum_accepted_steps: 131_072,
            velocity_absolute_tolerance_m_s: 1e-6,
            temperature_absolute_tolerance_k: 1e-5,
            energy_absolute_tolerance_j: 1e-7,
            transport_absolute_tolerance_m2: 1e-8,
            relative_tolerance: 1e-3,
            window_energy_tolerance_j: 1e-6,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ThermalFilmTimeError {
    pub maximum_velocity_difference_m_s: f64,
    pub maximum_temperature_difference_k: f64,
    pub maximum_energy_difference_j: f64,
    pub transport_difference_m2: f64,
    /// Fine local estimate with actual unequal-half correction. Not a global,
    /// spatial, material, phase or turbulent-transition accuracy certificate.
    pub normalized: f64,
}

#[derive(Debug)]
pub struct ThermalFilmAdaptiveReport<T> {
    pub time_s: f64,
    pub duration_s: f64,
    pub comparisons: usize,
    pub accuracy_rejections: usize,
    pub nonlinear_rejections: usize,
    pub accepted_steps: usize,
    pub native_events_processed: usize,
    pub maximum_accepted_normalized_error: f64,
    pub last_error: Option<ThermalFilmTimeError>,
    pub suggested_dt_s: f64,
    /// Includes coarse, rejected and failed original thermal trials.
    pub endpoint_evaluations: usize,
    pub momentum_material_evaluations: u64,
    pub last_exchange: Option<T>,
    pub transported_volume_per_width_m2: f64,
    /// In original film_cells() order; only accepted fine-step deposits.
    pub deposited_cell_heat_j: Vec<f64>,
    /// Original net cell heat before viscosity, including internal conduction.
    /// Distinct thermal mechanisms within the SAME cell can still cancel;
    /// no separate per-boundary or per-mechanism accuracy claim is made.
    pub other_cell_heat_j: BTreeMap<CellId, f64>,
    pub deposited_viscous_heat_j: f64,
    pub viscous_heat_j: f64,
    pub body_work_j: f64,
    pub wall_work_j: f64,
    pub surface_work_j: f64,
    pub other_heat_j: f64,
    pub stored_energy_change_j: f64,
    pub kinetic_energy_change_j: f64,
    pub thermal_energy_residual_j: f64,
    pub heat_partition_residual_j: f64,
    pub combined_energy_residual_j: f64,
}

#[derive(Debug)]
pub enum ThermalFilmAdaptiveError<E> {
    Film(ThermalFilmError<E>),
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
        error: ThermalFilmTimeError,
    },
}
impl<E> From<ThermalFilmError<E>> for ThermalFilmAdaptiveError<E> {
    fn from(error: ThermalFilmError<E>) -> Self {
        Self::Film(error)
    }
}
impl<E> From<DepositionError> for ThermalFilmAdaptiveError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Film(error.into())
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for ThermalFilmAdaptiveError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive thermal film: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for ThermalFilmAdaptiveError<E> {}

impl<'m> ThermalCarreauFilm<'m> {
    fn fork_adaptive<E>(
        &self,
        check: &dyn Fn() -> Result<(), DepositionError>,
    ) -> Result<Self, ThermalFilmError<E>> {
        check()?;
        let make_flow = || -> Result<TransientCarreauFilm, ThermalFilmError<E>> {
            Ok(TransientCarreauFilm::from_profile(
                *self.flow.material().law(),
                self.flow.material().temperature_k(),
                TransientCarreauFilmSpec {
                    density_kg_m3: self.flow.density_kg_m3(),
                    thickness_m: self.flow.thickness_m(),
                    maximum_laminar_shear_rate_per_s: self
                        .flow
                        .material()
                        .maximum_shear_rate_per_s(),
                },
                self.flow.velocities_m_s().to_vec(),
                self.flow.accepted_time_s(),
            )?)
        };
        let mut result = Self {
            flow: make_flow()?,
            trial: make_flow()?,
            thermal: self.thermal.clone(),
            cells: self.cells.clone(),
            area: self.area,
            face_temperatures: filled(self.cells.len() + 1)?,
            cell_heat: filled(self.cells.len())?,
        };
        result.restore_adaptive(self, check)?;
        Ok(result)
    }

    fn restore_adaptive<E>(
        &mut self,
        original: &Self,
        check: &dyn Fn() -> Result<(), DepositionError>,
    ) -> Result<(), ThermalFilmError<E>> {
        check()?;
        // Copy the original already-admitted immutable source preparations,
        // never scratch from a preceding failed thermal endpoint query.
        self.flow.material = original.flow.material;
        self.flow
            .face_materials
            .clone_from(&original.flow.face_materials);
        self.flow.flow.material.kinematic_viscosity_m2_s =
            original.flow.flow.material.kinematic_viscosity_m2_s;
        self.flow.flow.shear_conductance = original.flow.flow.shear_conductance;
        self.flow.flow.diffusion_time = original.flow.flow.diffusion_time;
        self.flow.flow.factored_duration = None;
        self.flow
            .flow
            .velocity
            .copy_from_slice(&original.flow.flow.velocity);
        self.flow.flow.time = original.flow.flow.time;
        self.thermal = original.thermal.clone();
        check()?;
        Ok(())
    }

    fn apply_same_time_events<E>(
        &mut self,
        tolerance: f64,
        check: &dyn Fn() -> Result<(), DepositionError>,
    ) -> Result<usize, ThermalFilmError<E>> {
        if self.thermal.next_event_time_s() != Some(self.thermal.time_s()) {
            return Ok(0);
        }
        let (next, audit) = self.thermal.advance_to(
            self.thermal.time_s(),
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: tolerance,
            },
            || check(),
        )?;
        let same = self.thermal.cells().len() == next.cells().len()
            && self
                .thermal
                .cells()
                .iter()
                .zip(next.cells())
                .all(|((id, a), (next_id, b))| {
                    id == next_id
                        && a.segment_id == b.segment_id
                        && a.body_id == b.body_id
                        && a.material_id == b.material_id
                        && a.mass_kg.to_bits() == b.mass_kg.to_bits()
                        && a.specific_h_j_kg.to_bits() == b.specific_h_j_kg.to_bits()
                        && a.geometry == b.geometry
                });
        if !same
            || audit.incoming_mass_kg != 0.0
            || audit.incoming_energy_j != 0.0
            || audit.external_heat_j != 0.0
            || audit.stored_energy_change_j != 0.0
        {
            return Err(ThermalFilmError::Invalid(
                "adaptive film events must retain EVERY native cell, material, mass, geometry and stored energy",
            ));
        }
        check()?;
        self.thermal = next;
        Ok(audit.events)
    }

    /// Advance a complete scheduled window atomically with the original
    /// caloric and momentum solvers. Drive samples use each actual interval's
    /// midpoint; heat sources use its endpoint candidate and explicit segment.
    /// Both callbacks must be repeatable and free of external side effects.
    ///
    /// The source event list is finite, nonnegative and strictly increasing;
    /// past events establish the current segment. Native boundary/contact
    /// events split intervals and use advance_to(current) for their handoff.
    /// All native cells, materials, reference masses and geometry must remain
    /// unchanged. A target-time event stays pending at its original left limit.
    /// Birth/remap, changing film thickness, phase occupancy and advection are
    /// not supplied by this fixed-inventory, wall-normal model.
    ///
    /// Error compares velocity, caloric energy CHANGES and temperature CHANGES,
    /// signed transport, every cell's physical heat, and separate body/wall/
    /// surface work. Signed exchanges AND their two-half magnitudes are checked
    /// to expose opposite work. Distinct heat mechanisms inside one cell can
    /// still cancel. Only accepted fine exchanges enter returned accounts.
    /// Original nonlinear/algebraic gates remain independent of time accuracy.
    ///
    /// Thermal controls' endpoint and material ceilings apply to this ENTIRE
    /// call, including discarded comparisons. Only outer nonlinear convergence
    /// and temporal error are retryable; sources, support, conservation,
    /// cancellation, work ceilings and representable step floors are terminal.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_adaptive<T, E>(
        &mut self,
        target_s: f64,
        source_events_s: &[f64],
        limits: FilmStepLimits,
        momentum: TransientCarreauControls,
        controls: ThermalFilmControls,
        adaptive: ThermalFilmAdaptiveControls,
        drive: impl Fn(ThermalFilmQuery) -> Result<FilmDrive, E>,
        resolve: impl Fn(
            ThermalFilmQuery,
            &DepositionSimulation<'m>,
            &[CellId],
        ) -> Result<ThermalFilmExchange<T>, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ThermalFilmAdaptiveReport<T>, ThermalFilmAdaptiveError<E>> {
        let checkpoint = RefCell::new(cancelled);
        let stopped = Cell::new(false);
        let check = || {
            if stopped.get() || (checkpoint.borrow_mut())() {
                stopped.set(true);
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        };
        check()?;
        let start_s = self.thermal.time_s();
        if !target_s.is_finite()
            || target_s < start_s
            || ![
                adaptive.initial_dt_s,
                adaptive.minimum_dt_s,
                adaptive.maximum_dt_s,
                adaptive.velocity_absolute_tolerance_m_s,
                adaptive.temperature_absolute_tolerance_k,
                adaptive.energy_absolute_tolerance_j,
                adaptive.transport_absolute_tolerance_m2,
            ]
            .iter()
            .all(|v| v.is_finite() && *v > 0.0)
            || adaptive.minimum_dt_s > adaptive.initial_dt_s
            || adaptive.initial_dt_s > adaptive.maximum_dt_s
            || !adaptive.relative_tolerance.is_finite()
            || adaptive.relative_tolerance < 0.0
            || !adaptive.window_energy_tolerance_j.is_finite()
            || adaptive.window_energy_tolerance_j < 0.0
            || !(1..=1_000_000).contains(&adaptive.maximum_trials)
            || !(2..=2_000_000).contains(&adaptive.maximum_accepted_steps)
            || self.thermal.cells().len() > controls.maximum_thermal_cells
            || !(1..=1_000_000).contains(&controls.maximum_thermal_cells)
            || source_events_s.len() > adaptive.maximum_trials
            || source_events_s.iter().any(|t| !t.is_finite() || *t < 0.0)
            || source_events_s.windows(2).any(|w| w[0] >= w[1])
        {
            return Err(ThermalFilmAdaptiveError::Invalid(
                "finite ordered time/source schedule, positive accuracy scales and bounded work required",
            ));
        }
        let mut next = self.fork_adaptive(&check)?;
        let mut coarse = self.fork_adaptive(&check)?;
        let mut fine = self.fork_adaptive(&check)?;
        let work = ThermalFilmWork::default();
        let mut report = ThermalFilmAdaptiveReport {
            time_s: start_s,
            duration_s: target_s - start_s,
            comparisons: 0,
            accuracy_rejections: 0,
            nonlinear_rejections: 0,
            accepted_steps: 0,
            native_events_processed: 0,
            maximum_accepted_normalized_error: 0.0,
            last_error: None,
            suggested_dt_s: adaptive.initial_dt_s,
            endpoint_evaluations: 0,
            momentum_material_evaluations: 0,
            last_exchange: None,
            transported_volume_per_width_m2: 0.0,
            deposited_cell_heat_j: vec![0.0; self.cells.len()],
            other_cell_heat_j: self.thermal.cells().keys().map(|id| (*id, 0.0)).collect(),
            deposited_viscous_heat_j: 0.0,
            viscous_heat_j: 0.0,
            body_work_j: 0.0,
            wall_work_j: 0.0,
            surface_work_j: 0.0,
            other_heat_j: 0.0,
            stored_energy_change_j: 0.0,
            kinetic_energy_change_j: 0.0,
            thermal_energy_residual_j: 0.0,
            heat_partition_residual_j: 0.0,
            combined_energy_residual_j: 0.0,
        };
        let mut integrals = [Dd::from_f64(0.0); 7];
        let mut deposited = vec![Dd::from_f64(0.0); self.cells.len()];
        let mut other: BTreeMap<_, _> = self
            .thermal
            .cells()
            .keys()
            .map(|id| (*id, Dd::from_f64(0.0)))
            .collect();
        let guess = vec![0.0; self.thermal.cells().len()];
        let step_budget = adaptive
            .maximum_accepted_steps
            .min(self.thermal.budget().max_substeps);
        let mut span = adaptive.initial_dt_s;
        while next.thermal.time_s() < target_s {
            check()?;
            report.native_events_processed +=
                next.apply_same_time_events(controls.total_energy_tolerance_j, &check)?;
            let start = next.thermal.time_s();
            let segment = source_events_s.partition_point(|t| *t <= start);
            let source_stop = source_events_s.get(segment).copied().unwrap_or(target_s);
            let native_stop = next.thermal.next_event_time_s().unwrap_or(target_s);
            let stop = target_s.min(source_stop).min(native_stop);
            if stop <= start {
                return Err(ThermalFilmAdaptiveError::Invalid(
                    "native event timeline did not advance",
                ));
            }
            if report.comparisons == adaptive.maximum_trials {
                return Err(ThermalFilmAdaptiveError::TrialBudget);
            }
            if report.accepted_steps + 2 > step_budget {
                return Err(ThermalFilmAdaptiveError::StepBudget);
            }
            let unclipped = span;
            let proposed = if span >= stop - start {
                stop
            } else {
                (start + span).min(stop)
            };
            let end = if proposed < stop && proposed.next_up() == stop {
                stop
            } else {
                proposed
            };
            let middle = start + 0.5 * (end - start);
            if !(start < middle && middle < end) {
                return Err(ThermalFilmAdaptiveError::StepUnderflow {
                    time_s: start,
                    span_s: end - start,
                });
            }
            let correction = ConductionTimeOrder::First
                .fine_error_correction(middle - start, end - middle)
                .ok_or(ThermalFilmAdaptiveError::StepUnderflow {
                    time_s: start,
                    span_s: end - start,
                })?;
            report.comparisons += 1;
            coarse.restore_adaptive(&next, &check)?;
            fine.restore_adaptive(&next, &check)?;
            let run = |candidate: &mut Self, endpoint: f64| {
                let begin = candidate.thermal.time_s();
                let query = ThermalFilmQuery {
                    start_s: begin,
                    midpoint_s: begin + 0.5 * (endpoint - begin),
                    endpoint_s: endpoint,
                    segment,
                };
                check().map_err(ThermalFilmError::from)?;
                let forcing = drive(query).map_err(ThermalFilmError::Source)?;
                candidate.advance_metered(
                    endpoint,
                    forcing,
                    limits,
                    momentum,
                    controls,
                    &guess,
                    |thermal, cells| resolve(query, thermal, cells),
                    || check().is_err(),
                    &work,
                )
            };
            let attempted = (|| {
                let c = run(&mut coarse, end)?;
                let a = run(&mut fine, middle)?;
                let b = run(&mut fine, end)?;
                Ok::<_, ThermalFilmError<E>>((c, a, b))
            })();
            let (c, a, b) = match attempted {
                Ok(value) => value,
                Err(ThermalFilmError::NotConverged(reason)) => {
                    if end - start <= adaptive.minimum_dt_s {
                        return Err(ThermalFilmError::NotConverged(reason).into());
                    }
                    report.nonlinear_rejections += 1;
                    span = (0.5 * (end - start)).max(adaptive.minimum_dt_s);
                    continue;
                }
                Err(error) => return Err(error.into()),
            };
            let error = compare(&next, &coarse, &fine, &c, &a, &b, correction, adaptive)?;
            report.last_error = Some(error);
            if error.normalized > 1.0 {
                if end - start <= adaptive.minimum_dt_s {
                    return Err(ThermalFilmAdaptiveError::TimeAccuracy {
                        time_s: start,
                        span_s: end - start,
                        error,
                    });
                }
                report.accuracy_rejections += 1;
                span = ((end - start)
                    * ConductionTimeOrder::First
                        .rejected_reduction_factor(error.normalized)
                        .unwrap())
                .max(adaptive.minimum_dt_s);
                continue;
            }
            report.maximum_accepted_normalized_error = report
                .maximum_accepted_normalized_error
                .max(error.normalized);
            for accepted in [&a, &b] {
                let channels = channels(accepted, self.area);
                for (sum, value) in integrals.iter_mut().zip(channels) {
                    *sum = *sum + Dd::from_f64(value);
                }
                for (sum, value) in deposited.iter_mut().zip(&accepted.deposited_cell_heat_j) {
                    *sum = *sum + Dd::from_f64(*value);
                }
                for (id, value) in &accepted.other_cell_heat_j {
                    let sum = other.get_mut(id).expect("unchanged native inventory");
                    *sum = *sum + Dd::from_f64(*value);
                }
            }
            report.last_exchange = Some(b.exchange);
            report.accepted_steps += 2;
            core::mem::swap(&mut next, &mut fine);
            span = ((end - start)
                * ConductionTimeOrder::First
                    .accepted_growth_factor(error.normalized)
                    .unwrap())
            .clamp(adaptive.minimum_dt_s, adaptive.maximum_dt_s);
            // Restore the uncut proposal after an exact event/window tail.
            // Its small clipped width must not create repeated decimal tails.
            if end == stop {
                span = span.max(unclipped).min(adaptive.maximum_dt_s);
            }
        }
        check()?;
        let values = integrals
            .map(|v| finite(v.to_f64()))
            .into_iter()
            .collect::<Result<Vec<_>, _>>()?;
        report.body_work_j = values[0];
        report.wall_work_j = values[1];
        report.surface_work_j = values[2];
        report.other_heat_j = values[3];
        report.viscous_heat_j = values[4];
        report.transported_volume_per_width_m2 = values[5];
        report.deposited_viscous_heat_j = values[6];
        report.deposited_cell_heat_j = deposited
            .into_iter()
            .map(|v| finite(v.to_f64()))
            .collect::<Result<_, _>>()?;
        report.other_cell_heat_j = other
            .into_iter()
            .map(|(id, v)| Ok((id, finite(v.to_f64())?)))
            .collect::<Result<_, DepositionError>>()?;
        let mut storage = Dd::from_f64(0.0);
        for (id, old) in self.thermal.cells() {
            let new = &next.thermal.cells()[id];
            storage = storage
                + Dd::from_f64(old.mass_kg)
                    * (Dd::from_f64(new.specific_h_j_kg) - Dd::from_f64(old.specific_h_j_kg));
        }
        let mut kinetic = Dd::from_f64(0.0);
        let mass = Dd::from_f64(self.area)
            * Dd::from_f64(self.flow.density_kg_m3())
            * Dd::from_f64(self.flow.spacing_m());
        for (&old, &new) in self
            .flow
            .velocities_m_s()
            .iter()
            .zip(next.flow.velocities_m_s())
        {
            kinetic = kinetic
                + mass
                    * Dd::from_f64(0.5 * old + 0.5 * new)
                    * (Dd::from_f64(new) - Dd::from_f64(old));
        }
        report.stored_energy_change_j = finite(storage.to_f64())?;
        report.kinetic_energy_change_j = finite(kinetic.to_f64())?;
        report.thermal_energy_residual_j =
            finite((storage - integrals[3] - integrals[6]).to_f64())?;
        report.heat_partition_residual_j = finite((integrals[6] - integrals[4]).to_f64())?;
        report.combined_energy_residual_j = finite(
            (storage + kinetic - integrals[3] - integrals[0] - integrals[1] - integrals[2])
                .to_f64(),
        )?;
        for residual in [
            report.thermal_energy_residual_j,
            report.heat_partition_residual_j,
            report.combined_energy_residual_j,
        ] {
            if residual.abs() > adaptive.window_energy_tolerance_j {
                return Err(ThermalFilmError::CombinedEnergy {
                    residual_j: residual,
                    tolerance_j: adaptive.window_energy_tolerance_j,
                }
                .into());
            }
        }
        report.time_s = next.thermal.time_s();
        report.suggested_dt_s = span;
        report.endpoint_evaluations = work.endpoints.get();
        report.momentum_material_evaluations = work.materials.get();
        check()?;
        *self = next;
        Ok(report)
    }
}

fn channels<T>(report: &ThermalFilmReport<T>, area: f64) -> [f64; 7] {
    let flow = report.flow.flow;
    [
        area * flow.body_work_j_m2,
        area * flow.wall_work_j_m2,
        area * flow.surface_work_j_m2,
        report.other_heat_j,
        area * flow.viscous_heat_j_m2,
        flow.duration_s * flow.midpoint_flow_per_width_m2_s,
        report.deposited_viscous_heat_j,
    ]
}

#[allow(clippy::too_many_arguments)]
fn compare<T>(
    old: &ThermalCarreauFilm<'_>,
    coarse: &ThermalCarreauFilm<'_>,
    fine: &ThermalCarreauFilm<'_>,
    c: &ThermalFilmReport<T>,
    a: &ThermalFilmReport<T>,
    b: &ThermalFilmReport<T>,
    correction: f64,
    controls: ThermalFilmAdaptiveControls,
) -> Result<ThermalFilmTimeError, DepositionError> {
    let mut result = ThermalFilmTimeError::default();
    let mut normalized: f64 = 0.0;
    let mut measure = |coarse: f64, fine: f64, absolute: f64| -> Result<f64, DepositionError> {
        let difference = finite(correction * (coarse - fine).abs())?;
        let scale = finite(absolute + controls.relative_tolerance * coarse.abs().max(fine.abs()))?;
        normalized = normalized.max(difference / scale);
        Ok(difference)
    };
    for ((&initial, &c), &f) in old
        .flow
        .velocities_m_s()
        .iter()
        .zip(coarse.flow.velocities_m_s())
        .zip(fine.flow.velocities_m_s())
    {
        result.maximum_velocity_difference_m_s =
            result.maximum_velocity_difference_m_s.max(measure(
                c - initial,
                f - initial,
                controls.velocity_absolute_tolerance_m_s,
            )?);
    }
    for (id, initial) in old.thermal.cells() {
        let c = &coarse.thermal.cells()[id];
        let f = &fine.thermal.cells()[id];
        let c_energy = initial.mass_kg * (c.specific_h_j_kg - initial.specific_h_j_kg);
        let f_energy = initial.mass_kg * (f.specific_h_j_kg - initial.specific_h_j_kg);
        result.maximum_energy_difference_j = result.maximum_energy_difference_j.max(measure(
            c_energy,
            f_energy,
            controls.energy_absolute_tolerance_j,
        )?);
        let initial_t = old.thermal.temperature_k(*id)?;
        result.maximum_temperature_difference_k =
            result.maximum_temperature_difference_k.max(measure(
                coarse.thermal.temperature_k(*id)? - initial_t,
                fine.thermal.temperature_k(*id)? - initial_t,
                controls.temperature_absolute_tolerance_k,
            )?);
    }
    for (index, ((c, a), b)) in channels(c, old.area)
        .into_iter()
        .zip(channels(a, old.area))
        .zip(channels(b, old.area))
        .enumerate()
    {
        let absolute = if index == 5 {
            controls.transport_absolute_tolerance_m2
        } else {
            controls.energy_absolute_tolerance_j
        };
        let difference =
            measure(c, a + b, absolute)?.max(measure(c.abs(), a.abs() + b.abs(), absolute)?);
        if index == 5 {
            result.transport_difference_m2 = difference;
        } else {
            result.maximum_energy_difference_j = result.maximum_energy_difference_j.max(difference);
        }
    }
    for ((&c, &a), &b) in c
        .deposited_cell_heat_j
        .iter()
        .zip(&a.deposited_cell_heat_j)
        .zip(&b.deposited_cell_heat_j)
    {
        result.maximum_energy_difference_j = result.maximum_energy_difference_j.max(measure(
            c,
            a + b,
            controls.energy_absolute_tolerance_j,
        )?);
    }
    for (id, &c) in &c.other_cell_heat_j {
        let a = a.other_cell_heat_j[id];
        let b = b.other_cell_heat_j[id];
        result.maximum_energy_difference_j = result
            .maximum_energy_difference_j
            .max(measure(c, a + b, controls.energy_absolute_tolerance_j)?)
            .max(measure(
                c.abs(),
                a.abs() + b.abs(),
                controls.energy_absolute_tolerance_j,
            )?);
    }
    result.normalized = normalized;
    Ok(result)
}
