//! Adaptive backward Euler without replacing the physical enthalpy equation.
//!
//! Each trial compares one full step with two half steps from the SAME accepted
//! history. The fine endpoint is accepted without Richardson extrapolation:
//! extrapolation can leave a latent plateau or violate the step's heat balance.
//! The three error tests cover specific enthalpy, temperature and liquid fraction.
//! In particular, equal temperatures do not hide different latent storage.
//!
//! Relative scales use changes over the trial, not absolute enthalpy or kelvin.
//! Thus an enthalpy reference offset cannot loosen the requested accuracy. These
//! are local time-error estimates, not bounds on accumulated or spatial error.
//! Only the two accepted half steps contribute heat to the running balance.
//!
//! `advance_to` lands exactly on a requested event time. Source, boundary and
//! contact laws are fixed during that call; change them only between event calls.
//! A trial budget, minimum-step stop, cancellation or solver refusal preserves
//! the last accepted field. Calling again continues the same controller state.
//! This full-window controller accepts directly assembled physical laws;
//! authored boundary/source schedules keep their existing history stepper.

use fs_exec::Cx;

use super::{
    EnthalpyBackwardEuler, EnthalpyError, EnthalpyStepConfig, EnthalpyStepSolution, finite, poll,
    validate_work,
};
use crate::radiation::{PointwiseAmbientRadiation, PointwiseEnclosureRadiation};
use crate::{ConductionProblem, ThermalInterfaces, assemble::ASSEMBLY_TILE};

/// Explicit local time-error tolerances and physical endpoint work limits.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AdaptiveEnthalpyConfig {
    /// Absolute specific-enthalpy error allowance, J/kg.
    pub enthalpy_absolute_tolerance_j_kg: f64,
    /// Absolute temperature error allowance, K.
    pub temperature_absolute_tolerance_k: f64,
    /// Absolute liquid mass-fraction error allowance, dimensionless.
    pub liquid_fraction_absolute_tolerance: f64,
    /// Relative allowance applied to each field's change over the trial.
    pub relative_tolerance: f64,
    /// Smallest requested full trial, s. An event may clip it further, but
    /// a rejected event-clipped trial cannot silently relax this floor.
    pub minimum_step_s: f64,
    /// Largest requested full trial, s.
    pub maximum_step_s: f64,
    /// Unchanged nonlinear, Krylov and per-endpoint energy gates.
    pub endpoint: EnthalpyStepConfig,
}

/// The measured full-step/two-half-step discrepancy at one common endpoint.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnthalpyTimeError {
    pub maximum_enthalpy_difference_j_kg: f64,
    pub maximum_temperature_difference_k: f64,
    pub maximum_liquid_fraction_difference: f64,
    /// Maximum componentwise difference divided by its declared allowance.
    /// Acceptance requires this value to be at most one.
    pub normalized: f64,
}

/// A partial trajectory is distinguishable from reaching the requested time.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdaptiveEnthalpyStatus {
    ReachedEnd,
    TrialLimit,
    MinimumStep,
}

/// Per-call work and the whole accepted trajectory's physical heat balance.
#[derive(Debug, Clone, PartialEq)]
pub struct AdaptiveEnthalpyReport {
    pub status: AdaptiveEnthalpyStatus,
    pub time_s: f64,
    /// Accepted full trials; each contains two physical backward-Euler steps.
    pub accepted: usize,
    pub rejected: usize,
    /// Actual attempted endpoint solves, including rejected comparison work.
    pub endpoint_solves: usize,
    pub last_error: Option<EnthalpyTimeError>,
    /// A recoverable nonlinear/energy refusal on the latest rejected trial.
    pub last_endpoint_refusal: Option<EnthalpyError>,
    /// Direct reference-mass enthalpy change from the ORIGINAL initial field.
    pub stored_energy_change_j: f64,
    /// Time-integrated external heat from accepted fine endpoints only.
    pub integrated_net_input_j: f64,
    /// Measured whole-window storage minus integrated external heat, J.
    /// This is not bounded by one endpoint's energy tolerance.
    pub energy_residual_j: f64,
}

/// In-memory resumable transport state bound to one immutable spatial owner.
/// Cloning forks the accepted history and controller; it does not change mesh,
/// reference masses, phase charts or conductivity assignments on the owner.
#[derive(Debug, Clone)]
pub struct AdaptiveEnthalpyState<'s, 'm, 'c> {
    storage: &'s EnthalpyBackwardEuler<'m, 'c>,
    config: AdaptiveEnthalpyConfig,
    initial_h: Vec<f64>,
    initial_temperature: Vec<f64>,
    initial_liquid: Vec<f64>,
    endpoint: Option<EnthalpyStepSolution>,
    time_s: f64,
    next_step_s: f64,
    accepted: usize,
    rejected: usize,
    endpoint_solves: usize,
    input_sum_j: f64,
    input_correction_j: f64,
    stored_energy_change_j: f64,
    last_error: Option<EnthalpyTimeError>,
    last_endpoint_refusal: Option<EnthalpyError>,
}

impl<'m, 'c> EnthalpyBackwardEuler<'m, 'c> {
    /// Prepare adaptive history, resolving every initial node through its actual
    /// phase chart. No transport solve or boundary-temperature projection occurs.
    pub fn adaptive_state<'s>(
        &'s self,
        cx: &Cx<'_>,
        initial_specific_h: &[f64],
        initial_time_s: f64,
        initial_step_s: f64,
        config: AdaptiveEnthalpyConfig,
    ) -> Result<AdaptiveEnthalpyState<'s, 'm, 'c>, EnthalpyError> {
        poll(cx, 0)?;
        config.validate(initial_step_s)?;
        validate_work(self.masses.len(), initial_step_s, config.endpoint)?;
        if !initial_time_s.is_finite() || initial_specific_h.len() != self.masses.len() {
            return Err(EnthalpyError::InvalidInput(
                "adaptive history needs a finite time and one enthalpy per reference vertex",
            ));
        }
        let mut temperature = Vec::with_capacity(initial_specific_h.len());
        let mut liquid = Vec::with_capacity(initial_specific_h.len());
        for (vertex, &h) in initial_specific_h.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            let state = self
                .curve_for_vertex(vertex)
                .state_at_specific_enthalpy(h)
                .map_err(|source| EnthalpyError::Phase { vertex, source })?;
            temperature.push(state.temperature_k());
            liquid.push(state.liquid_mass_fraction());
        }
        poll(cx, initial_specific_h.len())?;
        Ok(AdaptiveEnthalpyState {
            storage: self,
            config,
            initial_h: initial_specific_h.to_vec(),
            initial_temperature: temperature,
            initial_liquid: liquid,
            endpoint: None,
            time_s: initial_time_s,
            next_step_s: initial_step_s,
            accepted: 0,
            rejected: 0,
            endpoint_solves: 0,
            input_sum_j: 0.0,
            input_correction_j: 0.0,
            stored_energy_change_j: 0.0,
            last_error: None,
            last_endpoint_refusal: None,
        })
    }
}

impl AdaptiveEnthalpyConfig {
    fn validate(self, initial_step_s: f64) -> Result<(), EnthalpyError> {
        if ![
            self.enthalpy_absolute_tolerance_j_kg,
            self.temperature_absolute_tolerance_k,
            self.liquid_fraction_absolute_tolerance,
            self.minimum_step_s,
            self.maximum_step_s,
        ]
        .iter()
        .all(|v| v.is_finite() && *v > 0.0)
            || self.liquid_fraction_absolute_tolerance > 1.0
            || !self.relative_tolerance.is_finite()
            || self.relative_tolerance < 0.0
            || self.relative_tolerance > 1.0
            || self.minimum_step_s > self.maximum_step_s
            || !initial_step_s.is_finite()
            || !(self.minimum_step_s..=self.maximum_step_s).contains(&initial_step_s)
        {
            return Err(EnthalpyError::InvalidInput(
                "adaptive tolerances and step bounds must be finite, positive and ordered",
            ));
        }
        Ok(())
    }
}

impl AdaptiveEnthalpyState<'_, '_, '_> {
    #[must_use]
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    #[must_use]
    pub const fn next_step_s(&self) -> f64 {
        self.next_step_s
    }
    #[must_use]
    pub const fn accepted_steps(&self) -> usize {
        self.accepted
    }
    #[must_use]
    pub const fn rejected_steps(&self) -> usize {
        self.rejected
    }
    #[must_use]
    pub const fn endpoint_solves(&self) -> usize {
        self.endpoint_solves
    }
    #[must_use]
    pub fn specific_enthalpy_j_kg(&self) -> &[f64] {
        self.endpoint
            .as_ref()
            .map_or(&self.initial_h, |s| &s.specific_enthalpy_j_kg)
    }
    #[must_use]
    pub fn temperature(&self) -> &[f64] {
        self.endpoint
            .as_ref()
            .map_or(&self.initial_temperature, |s| &s.temperature)
    }
    #[must_use]
    pub fn liquid_mass_fraction(&self) -> &[f64] {
        self.endpoint
            .as_ref()
            .map_or(&self.initial_liquid, |s| &s.liquid_mass_fraction)
    }
    /// Last accepted HALF-step response. Its storage and energy fields cover
    /// only that half step; use the returned report for the whole trajectory.
    #[must_use]
    pub fn last_endpoint(&self) -> Option<&EnthalpyStepSolution> {
        self.endpoint.as_ref()
    }

    /// Advance under fixed source/boundary/contact laws to a declared event.
    /// `max_trials` bounds accepted plus rejected comparisons for THIS call.
    /// The final trial is clipped to the event, never advanced past it.
    pub fn advance_to(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        end_time_s: f64,
        max_trials: usize,
    ) -> Result<AdaptiveEnthalpyReport, EnthalpyError> {
        let storage = self.storage;
        let config = self.config.endpoint;
        self.advance_with(cx, end_time_s, max_trials, |h, dt| {
            storage.advance(cx, problem, interfaces, h, dt, config)
        })
    }

    /// The same controller with local P1 fourth-power ambient radiation.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to_with_pointwise_ambient_radiation(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        radiation: &PointwiseAmbientRadiation<'_>,
        end_time_s: f64,
        max_trials: usize,
    ) -> Result<AdaptiveEnthalpyReport, EnthalpyError> {
        let storage = self.storage;
        let config = self.config.endpoint;
        self.advance_with(cx, end_time_s, max_trials, |h, dt| {
            storage.advance_with_pointwise_ambient_radiation(
                cx, problem, interfaces, h, dt, config, radiation,
            )
        })
    }

    /// The same controller with reflected, spatially resolved enclosure radiation.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to_with_pointwise_enclosure_radiation(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        radiation: &PointwiseEnclosureRadiation<'_>,
        end_time_s: f64,
        max_trials: usize,
    ) -> Result<AdaptiveEnthalpyReport, EnthalpyError> {
        let storage = self.storage;
        let config = self.config.endpoint;
        self.advance_with(cx, end_time_s, max_trials, |h, dt| {
            storage.advance_with_pointwise_enclosure_radiation(
                cx, problem, interfaces, h, dt, config, radiation,
            )
        })
    }

    fn advance_with(
        &mut self,
        cx: &Cx<'_>,
        end: f64,
        max_trials: usize,
        mut advance: impl FnMut(&[f64], f64) -> Result<EnthalpyStepSolution, EnthalpyError>,
    ) -> Result<AdaptiveEnthalpyReport, EnthalpyError> {
        poll(cx, 0)?;
        if !end.is_finite() || end < self.time_s || !(end - self.time_s).is_finite() {
            return Err(EnthalpyError::InvalidInput(
                "adaptive end time must be finite and not precede history",
            ));
        }
        let before = (self.accepted, self.rejected, self.endpoint_solves);
        self.accepted
            .checked_add(max_trials)
            .and_then(|_| self.rejected.checked_add(max_trials))
            .and_then(|_| max_trials.checked_mul(3))
            .and_then(|n| self.endpoint_solves.checked_add(n))
            .ok_or(EnthalpyError::InvalidInput(
                "adaptive work counter overflow",
            ))?;
        let mut status = AdaptiveEnthalpyStatus::TrialLimit;
        for trial in 0..max_trials {
            poll(cx, trial)?;
            if self.time_s == end {
                break;
            }
            let remaining = end - self.time_s;
            let trial_end = if self.next_step_s >= remaining {
                end
            } else {
                finite(self.time_s + self.next_step_s)?
            };
            let dt = trial_end - self.time_s;
            let middle = self.time_s + 0.5 * dt;
            let first_dt = middle - self.time_s;
            let second_dt = trial_end - middle;
            if dt <= 0.0 || first_dt <= 0.0 || second_dt <= 0.0 {
                status = AdaptiveEnthalpyStatus::MinimumStep;
                break;
            }
            let candidate = (|| {
                self.endpoint_solves += 1;
                let coarse = advance(self.specific_enthalpy_j_kg(), dt)?;
                poll(cx, trial)?;
                self.endpoint_solves += 1;
                let first = advance(self.specific_enthalpy_j_kg(), first_dt)?;
                poll(cx, trial)?;
                self.endpoint_solves += 1;
                let fine = advance(&first.specific_enthalpy_j_kg, second_dt)?;
                let input_first = finite(first_dt * external_heat_w(&first)?)?;
                let input_second = finite(second_dt * external_heat_w(&fine)?)?;
                let error = estimate(cx, self, &coarse, &fine)?;
                Ok::<_, EnthalpyError>((fine, input_first, input_second, error))
            })();
            let (fine, input_first, input_second, error) = match candidate {
                Ok(candidate) => candidate,
                Err(
                    error @ (EnthalpyError::NotConverged(_) | EnthalpyError::EnergyBalance { .. }),
                ) => {
                    poll(cx, trial)?;
                    self.rejected += 1;
                    self.last_error = None;
                    self.last_endpoint_refusal = Some(error);
                    if !self.reduce_step(dt, 0.5) {
                        status = AdaptiveEnthalpyStatus::MinimumStep;
                        break;
                    }
                    continue;
                }
                Err(error) => return Err(error),
            };
            poll(cx, trial)?;
            if error.normalized > 1.0 {
                self.rejected += 1;
                self.last_error = Some(error);
                self.last_endpoint_refusal = None;
                let factor = (0.9 / error.normalized.sqrt()).clamp(0.1, 0.5);
                if !self.reduce_step(dt, factor) {
                    status = AdaptiveEnthalpyStatus::MinimumStep;
                    break;
                }
                continue;
            }
            // Prepare all accounting before publishing any physical state.
            // Coarse/rejected trial heat NEVER enters this sum.
            let (sum, correction) =
                compensated_add(self.input_sum_j, self.input_correction_j, input_first)?;
            let (sum, correction) = compensated_add(sum, correction, input_second)?;
            let mut stored = 0.0;
            let mut stored_correction = 0.0;
            for (vertex, ((&h, &initial), &mass)) in fine
                .specific_enthalpy_j_kg
                .iter()
                .zip(&self.initial_h)
                .zip(&self.storage.masses)
                .enumerate()
            {
                if vertex % ASSEMBLY_TILE == 0 {
                    poll(cx, vertex)?;
                }
                (stored, stored_correction) = compensated_add(
                    stored,
                    stored_correction,
                    finite(mass * finite(h - initial)?)?,
                )?;
            }
            finite(stored - sum)?;
            let factor = if error.normalized == 0.0 {
                2.0
            } else {
                (0.9 / error.normalized.sqrt()).clamp(0.5, 2.0)
            };
            let next = (dt * factor)
                .min(self.config.maximum_step_s)
                .max(self.config.minimum_step_s);
            poll(cx, trial)?;
            self.input_sum_j = sum;
            self.input_correction_j = correction;
            self.stored_energy_change_j = stored;
            self.endpoint = Some(fine);
            self.time_s = trial_end;
            self.next_step_s = next;
            self.accepted += 1;
            self.last_error = Some(error);
            self.last_endpoint_refusal = None;
        }
        if self.time_s == end {
            status = AdaptiveEnthalpyStatus::ReachedEnd;
        }
        Ok(AdaptiveEnthalpyReport {
            status,
            time_s: self.time_s,
            accepted: self.accepted - before.0,
            rejected: self.rejected - before.1,
            endpoint_solves: self.endpoint_solves - before.2,
            last_error: self.last_error,
            last_endpoint_refusal: self.last_endpoint_refusal.clone(),
            stored_energy_change_j: self.stored_energy_change_j,
            integrated_net_input_j: self.input_sum_j,
            energy_residual_j: self.stored_energy_change_j - self.input_sum_j,
        })
    }

    fn reduce_step(&mut self, attempted: f64, factor: f64) -> bool {
        if attempted <= self.config.minimum_step_s {
            return false;
        }
        let next = (attempted * factor).max(self.config.minimum_step_s);
        if next >= attempted || self.time_s + next == self.time_s {
            return false;
        }
        self.next_step_s = next;
        true
    }
}

fn estimate(
    cx: &Cx<'_>,
    state: &AdaptiveEnthalpyState<'_, '_, '_>,
    coarse: &EnthalpyStepSolution,
    fine: &EnthalpyStepSolution,
) -> Result<EnthalpyTimeError, EnthalpyError> {
    let mut error = EnthalpyTimeError {
        maximum_enthalpy_difference_j_kg: 0.0,
        maximum_temperature_difference_k: 0.0,
        maximum_liquid_fraction_difference: 0.0,
        normalized: 0.0,
    };
    let fields = [
        (
            &coarse.specific_enthalpy_j_kg,
            &fine.specific_enthalpy_j_kg,
            state.specific_enthalpy_j_kg(),
            state.config.enthalpy_absolute_tolerance_j_kg,
        ),
        (
            &coarse.temperature,
            &fine.temperature,
            state.temperature(),
            state.config.temperature_absolute_tolerance_k,
        ),
        (
            &coarse.liquid_mass_fraction,
            &fine.liquid_mass_fraction,
            state.liquid_mass_fraction(),
            state.config.liquid_fraction_absolute_tolerance,
        ),
    ];
    for (field, (full, half, old, atol)) in fields.into_iter().enumerate() {
        let mut maximum = 0.0_f64;
        for (vertex, ((&full, &half), &old)) in full.iter().zip(half).zip(old).enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            let difference = finite(half - full)?.abs();
            let change = finite(full - old)?.abs().max(finite(half - old)?.abs());
            let allowance = finite(atol + finite(state.config.relative_tolerance * change)?)?;
            let ratio = difference / allowance;
            // An infinite ratio means rejection, not a fabricated finite error.
            if ratio.is_nan() {
                return Err(EnthalpyError::NonFiniteArithmetic);
            }
            maximum = maximum.max(difference);
            error.normalized = error.normalized.max(ratio);
        }
        match field {
            0 => error.maximum_enthalpy_difference_j_kg = maximum,
            1 => error.maximum_temperature_difference_k = maximum,
            _ => error.maximum_liquid_fraction_difference = maximum,
        }
    }
    Ok(error)
}

fn external_heat_w(step: &EnthalpyStepSolution) -> Result<f64, EnthalpyError> {
    let mut net =
        finite(step.source_w + step.dirichlet_in_w - step.neumann_out_w - step.robin_out_w)?;
    if let Some(radiation) = &step.pointwise_radiation {
        net = finite(net - radiation.outward_heat_w)?;
    }
    if let Some(radiation) = &step.pointwise_enclosure_radiation {
        net = finite(net - radiation.outward_heat_w)?;
    }
    Ok(net)
}

fn compensated_add(sum: f64, correction: f64, value: f64) -> Result<(f64, f64), EnthalpyError> {
    let increment = finite(value - correction)?;
    let next = finite(sum + increment)?;
    Ok((next, finite((next - sum) - increment)?))
}
