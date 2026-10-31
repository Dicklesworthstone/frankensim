//! Bounded continuation of the original finite wall, film and gas transaction.
//!
//! Every accepted interval advances all three stores through the same duration.
//! The caller's body and gas heat inputs are uniform in time over this window;
//! cumulative source endpoints apportion them without renewing their totals on
//! retries. Phase step limits discard the entire proposed interval. The nominal
//! time grid resumes after a shortened interval reaches its original endpoint.
//! This is admission-driven subdivision, not a temporal-error estimator: refine
//! `maximum_step_s` to assess the complete coupled time discretization.

use super::{WallFilmThermalEnvironment, WallFilmThermalError};
use fs_exec::Cx;
use fs_material::{
    evaporation::gas::{GasEvaporationError, wall_film::GasWallFilmInventory},
    phase::{UniformEnthalpyStep, UniformEnthalpyStepInput},
    thermochemical_gas::energy::{GasEnergyError, parcel::ClosedGasParcel},
};

/// One physical-window work and conservation allowance, never renewed per step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallFilmWindowSolve {
    pub maximum_step_s: f64,
    /// Smallest retried interval. A final nominal-grid remainder may be shorter.
    pub minimum_step_s: f64,
    /// Includes rejected complete body/film/gas proposals. Each proposal also
    /// retains the original wall, sensible-heat and gas inversion work bounds.
    pub maximum_attempts: usize,
    pub energy_tolerance_j: f64,
    pub mass_tolerance_kg: f64,
}

#[derive(Debug)]
pub enum WallFilmWindowError {
    Invalid(&'static str),
    Step(WallFilmThermalError),
    Gas(GasEnergyError),
    Cancelled,
    WorkLimit {
        attempted_intervals: usize,
        accepted_time_s: f64,
    },
    StepSizeLimit {
        accepted_time_s: f64,
        requested_duration_s: f64,
        maximum_duration_s: f64,
        reason: &'static str,
    },
    EnergyBalance {
        residual_j: f64,
        tolerance_j: f64,
    },
    MassBalance {
        residual_kg: f64,
        tolerance_kg: f64,
    },
}
impl core::fmt::Display for WallFilmWindowError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite wall-film window: {self:?}")
    }
}
impl std::error::Error for WallFilmWindowError {}
impl From<GasEnergyError> for WallFilmWindowError {
    fn from(error: GasEnergyError) -> Self {
        Self::Gas(error)
    }
}

/// Compact accepted physical history; no cloned gas catalog per sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallFilmWindowSample {
    /// Relative to the beginning of this requested window [s].
    pub time_s: f64,
    pub wall_temperature_k: f64,
    pub gas_temperature_k: f64,
    pub film_mass_kg: f64,
    pub transferred_mass_kg: f64,
    pub film_formed: bool,
    /// Within-window event time for the accepted interval's frozen rate.
    pub dry_out_time_s: Option<f64>,
}

/// Original candidate states plus accounts spanning the entire requested time.
#[derive(Clone, Debug, PartialEq)]
pub struct WallFilmWindowReport {
    pub next_gas: ClosedGasParcel,
    pub next_film: GasWallFilmInventory,
    /// Includes the initial state and every accepted endpoint, in time order.
    pub samples: Vec<WallFilmWindowSample>,
    pub accepted_intervals: usize,
    pub attempted_intervals: usize,
    pub rejected_intervals: usize,
    /// Work inside accepted proposals only; rejected proposals still spend the
    /// `maximum_attempts` allowance and their own nested solve bounds.
    pub accepted_wall_trials: usize,
    pub accepted_sensible_trials: usize,
    pub evaporated_mass_kg: f64,
    pub condensed_mass_kg: f64,
    pub vapor_energy_into_gas_j: f64,
    pub latent_energy_j: f64,
    pub film_sensible_energy_j: f64,
    pub solid_to_gas_heat_j: f64,
    pub wall_heat_into_phase_system_j: f64,
    pub internal_body_heat_j: f64,
    pub external_gas_heat_j: f64,
    /// Accepted interval source sums minus the requested full-window totals.
    pub body_heat_apportionment_residual_j: f64,
    pub gas_heat_apportionment_residual_j: f64,
    pub radiation_heat_j: f64,
    pub wall_enthalpy_change_j: f64,
    /// Original owners' source-relative storage sums. Represented mass-reference
    /// rounding is bounded separately by actual endpoint inventory checks.
    pub film_enthalpy_change_j: f64,
    pub gas_internal_energy_change_j: f64,
    pub combined_energy_residual_j: f64,
    pub sum_absolute_energy_residual_j: f64,
    pub mass_residual_kg: f64,
    pub volatile_mass_residual_kg: f64,
    pub maximum_carrier_mass_residual_kg: f64,
}

#[derive(Clone, Copy, Default)]
struct Sum {
    main: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), WallFilmWindowError> {
        let (main, residual) = fs_math::eft::two_sum(self.main, value);
        self.main = main;
        self.correction += residual;
        finite(self.value())?;
        Ok(())
    }
    fn value(self) -> f64 {
        self.main + self.correction
    }
    fn difference(self, reference: f64) -> Result<f64, WallFilmWindowError> {
        let (difference, roundoff) = fs_math::eft::two_sum(self.main, -reference);
        finite(difference + (self.correction + roundoff))
    }
}
#[derive(Default)]
struct Totals {
    evaporated: Sum,
    condensed: Sum,
    vapor_energy: Sum,
    latent: Sum,
    film_sensible: Sum,
    sensible_heat: Sum,
    phase_wall_heat: Sum,
    external_wall_heat: Sum,
    internal_body_heat: Sum,
    external_gas_heat: Sum,
    radiation: Sum,
    wall_storage: Sum,
    film_storage: Sum,
    gas_storage: Sum,
    energy_residual: Sum,
    absolute_energy_residual: Sum,
}
fn finite(value: f64) -> Result<f64, WallFilmWindowError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(WallFilmWindowError::Invalid(
            "window arithmetic is unrepresentable",
        ))
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cx: &Cx<'_>) -> Result<(), WallFilmWindowError> {
    cx.checkpoint().map_err(|_| WallFilmWindowError::Cancelled)
}
fn increment(total: f64, start: f64, end: f64, duration: f64) -> Result<f64, WallFilmWindowError> {
    let at = |time| {
        if time == duration {
            total
        } else {
            total * (time / duration)
        }
    };
    let value = finite(at(end) - at(start))?;
    if total != 0.0 && value == 0.0 {
        return Err(WallFilmWindowError::Invalid(
            "nonzero interval heat disappears in the window source",
        ));
    }
    Ok(value)
}
fn count(total: usize, increment: usize) -> Result<usize, WallFilmWindowError> {
    total
        .checked_add(increment)
        .ok_or(WallFilmWindowError::Invalid("window work count overflows"))
}

impl WallFilmThermalEnvironment<'_> {
    /// Advance one uniform-source window, retaining real onset/dry-out history.
    ///
    /// `input.internal_heat_j` and `self.gas_heat_j` are TOTAL energies over
    /// `input.duration_s`, not energies to repeat at every accepted interval.
    /// The original one-step owner supplies all wall, phase and gas physics.
    /// Only its explicit phase `StepSizeLimit` is retryable here; source domains,
    /// failed nonlinear brackets, energy defects and cancellation still refuse.
    /// Failed work never publishes a partial endpoint or a partial physical sum.
    ///
    /// A separate call can change sources/ambient conditions at a known event.
    /// Film geometry stays with the wall across Dry/Wet transitions. No nominal
    /// output timestep is substituted for an accepted interval's physical clock.
    pub fn advance_window(
        &self,
        cx: &Cx<'_>,
        input: UniformEnthalpyStepInput<'_>,
        controls: WallFilmWindowSolve,
    ) -> Result<UniformEnthalpyStep<WallFilmWindowReport>, WallFilmWindowError> {
        poll(cx)?;
        if !positive(input.duration_s)
            || !input.internal_heat_j.is_finite()
            || !self.gas_heat_j.is_finite()
            || !positive(controls.maximum_step_s)
            || !positive(controls.minimum_step_s)
            || controls.minimum_step_s > controls.maximum_step_s
            || controls.maximum_attempts == 0
            || !positive(controls.energy_tolerance_j)
            || !positive(controls.mass_tolerance_kg)
        {
            return Err(WallFilmWindowError::Invalid(
                "positive finite window, steps, work and balance budgets required",
            ));
        }
        let mut state = input.initial;
        let mut gas = self.gas.clone();
        let mut film = self.film;
        let mut samples = vec![WallFilmWindowSample {
            time_s: 0.0,
            wall_temperature_k: state.temperature_k(),
            gas_temperature_k: gas.state().state().temperature,
            film_mass_kg: film.mass_kg(),
            transferred_mass_kg: 0.0,
            film_formed: false,
            dry_out_time_s: None,
        }];
        let mut totals = Totals::default();
        let mut time = 0.0;
        let mut grid_index = 1_usize;
        let mut grid_end = controls.maximum_step_s.min(input.duration_s);
        let mut trial_end = grid_end;
        let mut attempted = 0;
        let mut rejected = 0;
        let mut accepted_wall_trials = 0;
        let mut accepted_sensible_trials = 0;
        while time < input.duration_s {
            poll(cx)?;
            if attempted >= controls.maximum_attempts {
                return Err(WallFilmWindowError::WorkLimit {
                    attempted_intervals: attempted,
                    accepted_time_s: time,
                });
            }
            let dt = finite(trial_end - time)?;
            if dt <= 0.0 {
                return Err(WallFilmWindowError::Invalid(
                    "window clock cannot represent a forward interval",
                ));
            }
            let internal_heat =
                increment(input.internal_heat_j, time, trial_end, input.duration_s)?;
            let gas_heat = increment(self.gas_heat_j, time, trial_end, input.duration_s)?;
            let environment = WallFilmThermalEnvironment {
                gas: &gas,
                catalog: self.catalog,
                film,
                boundary: self.boundary,
                body: self.body,
                gas_heat_j: gas_heat,
                solve: self.solve,
            };
            attempted += 1;
            let candidate = match environment.advance(
                cx,
                UniformEnthalpyStepInput {
                    initial: state,
                    internal_heat_j: internal_heat,
                    duration_s: dt,
                    ..input
                },
            ) {
                Ok(candidate) => candidate,
                Err(WallFilmThermalError::Phase(GasEvaporationError::StepSizeLimit {
                    maximum_duration_s,
                    reason,
                })) => {
                    rejected += 1;
                    let next_dt = maximum_duration_s.min(0.5 * dt);
                    let next_end = time + next_dt;
                    if !positive(maximum_duration_s)
                        || !positive(next_dt)
                        || next_dt < controls.minimum_step_s
                        || next_end <= time
                        || next_end >= trial_end
                    {
                        return Err(WallFilmWindowError::StepSizeLimit {
                            accepted_time_s: time,
                            requested_duration_s: dt,
                            maximum_duration_s,
                            reason,
                        });
                    }
                    trial_end = next_end;
                    continue;
                }
                Err(WallFilmThermalError::Cancelled) => return Err(WallFilmWindowError::Cancelled),
                Err(error) => return Err(WallFilmWindowError::Step(error)),
            };
            let report = candidate.report;
            let phase = &report.phase;
            totals.evaporated.add(phase.transferred_mass_kg.max(0.0))?;
            totals
                .condensed
                .add((-phase.transferred_mass_kg).max(0.0))?;
            totals.vapor_energy.add(phase.energy_into_gas_j)?;
            totals.latent.add(phase.latent_energy_j)?;
            totals.film_sensible.add(phase.film_sensible_energy_j)?;
            totals.sensible_heat.add(report.solid_to_gas_heat_j)?;
            totals
                .phase_wall_heat
                .add(report.wall_heat_into_phase_system_j)?;
            totals.external_wall_heat.add(candidate.external_heat_j)?;
            totals.radiation.add(report.radiation_heat_j)?;
            totals.internal_body_heat.add(internal_heat)?;
            totals.external_gas_heat.add(gas_heat)?;
            totals.wall_storage.add(
                input.mass_kg
                    * (candidate.state.specific_enthalpy_j_kg() - state.specific_enthalpy_j_kg()),
            )?;
            totals.film_storage.add(phase.film_enthalpy_change_j)?;
            totals
                .gas_storage
                .add(report.gas_heat.stored_energy_change_j)?;
            totals.gas_storage.add(phase.gas_internal_energy_change_j)?;
            totals
                .energy_residual
                .add(report.combined_energy_residual_j)?;
            totals
                .absolute_energy_residual
                .add(report.sum_absolute_energy_residual_j)?;
            if totals.absolute_energy_residual.value() > controls.energy_tolerance_j {
                return Err(WallFilmWindowError::EnergyBalance {
                    residual_j: totals.absolute_energy_residual.value(),
                    tolerance_j: controls.energy_tolerance_j,
                });
            }
            accepted_wall_trials = count(accepted_wall_trials, report.trials)?;
            accepted_sensible_trials =
                count(accepted_sensible_trials, report.total_sensible_trials)?;
            let dry_out_time_s = phase
                .dry_out_time_s
                .map(|event| if event == dt { trial_end } else { time + event });
            if let Some(event) = dry_out_time_s {
                if !event.is_finite() || event <= time || event > trial_end {
                    return Err(WallFilmWindowError::Invalid(
                        "dry-out event is not representable on the window clock",
                    ));
                }
            }
            samples.push(WallFilmWindowSample {
                time_s: trial_end,
                wall_temperature_k: candidate.state.temperature_k(),
                gas_temperature_k: report.next_gas.state().state().temperature,
                film_mass_kg: report.next_film.mass_kg(),
                transferred_mass_kg: phase.transferred_mass_kg,
                film_formed: phase.film_formed,
                dry_out_time_s,
            });
            state = candidate.state;
            gas = report.next_gas;
            film = report.next_film;
            time = trial_end;
            if time == grid_end && time < input.duration_s {
                grid_index = count(grid_index, 1)?;
                grid_end = (grid_index as f64 * controls.maximum_step_s).min(input.duration_s);
                if !grid_end.is_finite() || grid_end <= time {
                    return Err(WallFilmWindowError::Invalid(
                        "nominal window grid cannot advance",
                    ));
                }
            }
            trial_end = grid_end;
        }
        poll(cx)?;
        let wall_storage = finite(
            input.mass_kg
                * (state.specific_enthalpy_j_kg() - input.initial.specific_enthalpy_j_kg()),
        )?;
        // Audit the actual complete wall endpoint independently of accumulated
        // interval caloric increments, without subtracting formation energies.
        let telescoping_residual = -totals.wall_storage.difference(wall_storage)?;
        let body_heat_apportionment_residual_j = totals
            .internal_body_heat
            .difference(input.internal_heat_j)?;
        let gas_heat_apportionment_residual_j =
            totals.external_gas_heat.difference(self.gas_heat_j)?;
        let energy_residual = finite(
            totals.energy_residual.value()
                + telescoping_residual
                + body_heat_apportionment_residual_j
                + gas_heat_apportionment_residual_j,
        )?;
        let absolute_energy_residual = finite(
            totals.absolute_energy_residual.value()
                + telescoping_residual.abs()
                + body_heat_apportionment_residual_j.abs()
                + gas_heat_apportionment_residual_j.abs(),
        )?;
        let wall_residual =
            finite(wall_storage - input.internal_heat_j - totals.external_wall_heat.value())?;
        if absolute_energy_residual > controls.energy_tolerance_j
            || wall_residual.abs() > controls.energy_tolerance_j
            || energy_residual.abs() > controls.energy_tolerance_j
        {
            return Err(WallFilmWindowError::EnergyBalance {
                residual_j: absolute_energy_residual
                    .max(wall_residual.abs())
                    .max(energy_residual.abs()),
                tolerance_j: controls.energy_tolerance_j,
            });
        }
        let film_mass_change = film.mass_kg() - self.film.mass_kg();
        let mass_residual = finite(gas.mass_kg() - self.gas.mass_kg() + film_mass_change)?;
        let volatile = &self.boundary.binding.species;
        let volatile_residual = finite(
            gas.species_mass_kg(volatile)? - self.gas.species_mass_kg(volatile)? + film_mass_change,
        )?;
        let mut maximum_carrier_residual: f64 = 0.0;
        for source in self.gas.model().thermodynamics().components() {
            poll(cx)?;
            if source.species() != volatile {
                let residual = finite(
                    gas.species_mass_kg(source.species())?
                        - self.gas.species_mass_kg(source.species())?,
                )?;
                maximum_carrier_residual = maximum_carrier_residual.max(residual.abs());
            }
        }
        let maximum_mass_residual = mass_residual
            .abs()
            .max(volatile_residual.abs())
            .max(maximum_carrier_residual);
        if maximum_mass_residual > controls.mass_tolerance_kg {
            return Err(WallFilmWindowError::MassBalance {
                residual_kg: maximum_mass_residual,
                tolerance_kg: controls.mass_tolerance_kg,
            });
        }
        poll(cx)?;
        Ok(UniformEnthalpyStep {
            state,
            external_heat_j: totals.external_wall_heat.value(),
            energy_residual_tolerance_j: controls.energy_tolerance_j,
            report: WallFilmWindowReport {
                next_gas: gas,
                next_film: film,
                samples,
                accepted_intervals: attempted - rejected,
                attempted_intervals: attempted,
                rejected_intervals: rejected,
                accepted_wall_trials,
                accepted_sensible_trials,
                evaporated_mass_kg: totals.evaporated.value(),
                condensed_mass_kg: totals.condensed.value(),
                vapor_energy_into_gas_j: totals.vapor_energy.value(),
                latent_energy_j: totals.latent.value(),
                film_sensible_energy_j: totals.film_sensible.value(),
                solid_to_gas_heat_j: totals.sensible_heat.value(),
                wall_heat_into_phase_system_j: totals.phase_wall_heat.value(),
                internal_body_heat_j: input.internal_heat_j,
                external_gas_heat_j: self.gas_heat_j,
                body_heat_apportionment_residual_j,
                gas_heat_apportionment_residual_j,
                radiation_heat_j: totals.radiation.value(),
                wall_enthalpy_change_j: wall_storage,
                film_enthalpy_change_j: totals.film_storage.value(),
                gas_internal_energy_change_j: totals.gas_storage.value(),
                combined_energy_residual_j: energy_residual,
                sum_absolute_energy_residual_j: absolute_energy_residual,
                mass_residual_kg: mass_residual,
                volatile_mass_residual_kg: volatile_residual,
                maximum_carrier_mass_residual_kg: maximum_carrier_residual,
            },
        })
    }
}
