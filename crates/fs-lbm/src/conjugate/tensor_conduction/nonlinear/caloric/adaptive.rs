//! Embedded temporal control around the original caloric tensor SDIRK2 owner.
//!
//! A trial still performs only the two native spatial stages. Its accepted
//! state and physical heat are the original second-order result; the embedded
//! first-order estimate exists only as a datum-free energy error coordinate.

use fs_conduction::enthalpy::conduction::implicit::adaptive::ConductionTimeOrder;
use fs_exec::CancelGate;
use std::borrow::Cow;

use super::super::super::sdirk::{TENSOR_SDIRK2_GAMMA, stage_clock};
use super::super::{TensorConductivityLaw, TensorNonlinearConfig};
use super::radiation::{TensorCaloricRadiationReport, advance_radiating_sdirk2_prepared};
use super::{
    TensorCaloricEndpoint, TensorCaloricSdirk2Step, TensorCaloricStorage,
    advance_caloric_tensor_conduction_sdirk2,
};
use crate::conjugate::radiation::{RadiationConfig, Radiators};
use crate::conjugate::transient::TransientRecord;
use crate::conjugate::{ChtError, SolidMaterial, ThermalSetup, VoxelDomain, finite, poll};

mod error;
use error::Sum;
pub(in crate::conjugate) use error::embedded_ratio;

/// Temporal accuracy and bounded work for one positive physical time window.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TensorCaloricAdaptiveConfig {
    /// One attempt's original Newton, Krylov and backtracking ceilings, shared
    /// across its two native SDIRK2 stages. Rejected work also consumes an
    /// attempt, so whole-window work is bounded by the products of these caps
    /// and maximum_attempts; failed nonlinear work does not escape that bound.
    pub step: TensorNonlinearConfig,
    pub initial_dt_s: f64,
    pub minimum_dt_s: f64,
    pub maximum_dt_s: f64,
    /// Includes accuracy rejections and failed native solves. Retained records
    /// are bounded by this same count. Zero permits no attempts.
    pub maximum_attempts: usize,
    /// Mandatory error allowance on original specific thermal energy, J/kg.
    /// This remains active on an exactly isothermal latent plateau.
    pub energy_absolute_tolerance_j_kg: f64,
    /// Scales the two energy INCREMENTS from the accepted start, never their
    /// absolute h/u datum. Must be finite and in [0,1).
    pub relative_tolerance: f64,
    /// Whole-window joule allowance for BOTH the final independently computed
    /// closure and the sum of absolute accepted native step closures.
    pub energy_tolerance_j: f64,
}

impl Default for TensorCaloricAdaptiveConfig {
    fn default() -> Self {
        Self {
            step: TensorNonlinearConfig::default(),
            initial_dt_s: 0.1,
            minimum_dt_s: 1e-9,
            maximum_dt_s: 1.0,
            maximum_attempts: 65_536,
            energy_absolute_tolerance_j_kg: 1e-3,
            relative_tolerance: 1e-3,
            energy_tolerance_j: 1e-6,
        }
    }
}

/// One actual native stage-time forcing query.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TensorCaloricForcingSample {
    pub time_s: f64,
    /// Attempted physical interval, including rejected attempts.
    pub interval_s: [f64; 2],
    /// True only at the ending sample of an explicitly declared transition.
    /// Return that interval's LEFT limiting data there; the next interval
    /// samples the new branch. A jump must not contaminate the preceding step.
    pub left_limit: bool,
}

/// Existing native source/boundary/contact/sink and thermostat data for a stage.
/// Constant samples may borrow their original arrays without per-stage clones.
#[derive(Clone)]
pub struct TensorCaloricForcing<'a> {
    pub setup: Cow<'a, ThermalSetup>,
    /// Explicit original energy for exactly the stage's prescribed cells.
    pub prescribed_specific_energy_j_kg: Cow<'a, [(usize, f64)]>,
}

impl<'a> TensorCaloricForcing<'a> {
    pub fn borrowed(setup: &'a ThermalSetup, prescribed: &'a [(usize, f64)]) -> Self {
        Self {
            setup: Cow::Borrowed(setup),
            prescribed_specific_energy_j_kg: Cow::Borrowed(prescribed),
        }
    }

    pub fn owned(setup: ThermalSetup, prescribed: Vec<(usize, f64)>) -> Self {
        Self {
            setup: Cow::Owned(setup),
            prescribed_specific_energy_j_kg: Cow::Owned(prescribed),
        }
    }
}

/// Borrow the original thermal part of a coupled physical-stage sample.
pub(in crate::conjugate) trait CaloricAdaptiveForcing {
    fn setup(&self) -> &ThermalSetup;
    fn prescribed(&self) -> &[(usize, f64)];
}

impl CaloricAdaptiveForcing for TensorCaloricForcing<'_> {
    fn setup(&self) -> &ThermalSetup {
        self.setup.as_ref()
    }

    fn prescribed(&self) -> &[(usize, f64)] {
        self.prescribed_specific_energy_j_kg.as_ref()
    }
}

/// One completed original spatial trial. Coupled proposals retain their actual
/// complete work separately from the accepted physical heat record.
pub(in crate::conjugate) struct CaloricAdaptiveTrial<Receipt> {
    pub thermal: TensorCaloricSdirk2Step,
    pub radiation: Option<[TensorCaloricRadiationReport; 2]>,
    pub receipt: Receipt,
    pub completed_iterations: usize,
    pub completed_radiation_solves: usize,
    /// Additional original embedded coordinate, such as local Joule energy.
    /// Zero retains the pure thermal controller's original arithmetic.
    pub additional_normalized_error: f64,
}

pub(in crate::conjugate) struct CaloricAdaptiveRun<Receipt> {
    pub thermal: TensorCaloricRadiatingAdaptiveSolution,
    pub receipts: Vec<Receipt>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TensorCaloricAdaptiveReport {
    pub attempts: usize,
    pub accepted_steps: usize,
    pub accuracy_rejections: usize,
    pub solver_rejections: usize,
    /// Krylov work in completed native trials, including accuracy rejections.
    /// Incomplete failed trials are bounded by step's caps, not misreported as
    /// zero work. They are counted separately in solver_rejections.
    pub completed_trial_iterations: usize,
    /// The embedded FIRST-order discrepancy, not a certified error bound on
    /// the accepted second-order solution or on individual heat channels.
    pub maximum_accepted_normalized_error: f64,
    pub absolute_accepted_closure_j: f64,
    /// Explicit stateless controller output for a subsequent physical window.
    pub suggested_dt_s: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct TensorCaloricAdaptiveSolution {
    pub endpoint: TensorCaloricEndpoint,
    /// Only accepted native intervals appear here; rejected heat is absent.
    pub records: Vec<TransientRecord>,
    /// Whole-window heat. Storage is recomputed from original initial/final
    /// energies and fixed masses. Heat channels are compensated sums of the
    /// accepted native records; iterations sums their native Krylov work.
    /// Closure uses separately compensated signed net heat, preserving small
    /// balances invisible in the independently rounded display channels.
    pub record: TransientRecord,
    pub report: TensorCaloricAdaptiveReport,
}

/// Accepted adaptive radiation on the same original-caloric thermal timeline.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorCaloricRadiatingAdaptiveSolution {
    /// Original endpoint, accepted records, whole-window heat and time control.
    pub thermal: TensorCaloricAdaptiveSolution,
    /// Native physical radiation reports for each accepted interval, in the
    /// same order as thermal.records. Rejected stage heat is never retained.
    /// Geometry counts identify the one shared prepared ray geometry; summing
    /// them would not measure ray-tracing work.
    pub radiation: Vec<[TensorCaloricRadiationReport; 2]>,
    /// Compensated sum of actual SDIRK stage radiation [J], positive outward.
    /// This is an informational component of thermal.record.boundary_outflow_j;
    /// it must not be deposited or debited from that account again.
    pub radiated_energy_j: f64,
    /// Native radiative heat solves in completed trials, including accuracy
    /// rejections. Incomplete failed work is bounded by maximum_attempts times
    /// the per-attempt limits and counted in thermal.report.solver_rejections.
    pub completed_trial_radiation_solves: usize,
}

#[derive(Debug)]
pub enum TensorCaloricAdaptiveError {
    Solve(ChtError),
    Invalid(&'static str),
    AttemptBudget {
        time_s: f64,
        attempts: usize,
    },
    TimeResolution {
        time_s: f64,
        boundary_s: f64,
    },
    TimeAccuracy {
        time_s: f64,
        dt_s: f64,
        normalized_error: f64,
    },
    EnergyBalance {
        closure_j: f64,
        absolute_accepted_closure_j: f64,
        tolerance_j: f64,
    },
}

impl From<ChtError> for TensorCaloricAdaptiveError {
    fn from(value: ChtError) -> Self {
        Self::Solve(value)
    }
}
impl core::fmt::Display for TensorCaloricAdaptiveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive caloric tensor conduction: {self:?}")
    }
}
impl std::error::Error for TensorCaloricAdaptiveError {}

/// Advance original caloric tensor conduction with an embedded SDIRK2(1) pair.
///
/// The native first stage has a=gamma*dt and e1=e0+a*f1. Its first-order
/// companion would have increment dt*f1=(e1-e0)/(a/dt). Compare that increment
/// with e2-e0, using only free cells and original energy CHANGE scales. The
/// native second-order e2 and its ledger remain the accepted physical state.
/// No extra spatial comparison solves, extrapolated physical states, apparent
/// capacities, or T queries at the numerical embedded energy are introduced.
///
/// The O(dt^2) embedded discrepancy uses the existing thermal controller's
/// square-root growth/reduction factors. It is an asymptotic LOCAL estimate,
/// not a certified error bound on the second-order endpoint, a global bound,
/// or a quadrature bound on cancelling external heat channels. Native spatial
/// and nonlinear residual tolerances still need to resolve the requested time
/// accuracy. Sharp material corners can reduce order and require refinement.
///
/// transition_times_s is a finite, strictly increasing nonnegative schedule.
/// Past and later transitions may remain in that schedule for repeated calls.
/// No interval crosses an interior transition. At an exact transition ending,
/// forcing MUST return the declared left limit; its next interval uses the new
/// branch. A terminal/event remainder may be shorter than minimum_dt_s, but
/// its temporal error gate is never waived and it cannot be reduced further.
///
/// Forcing is sampled at the two actual native stage times on every attempt.
/// It must be deterministic and read-only with respect to physical history.
/// FnMut permits caches/diagnostics, not publication of rejected trial state.
/// Long callbacks must respect the same cancellation scope. Prescribed masks
/// must agree across all attempts in a smooth segment; declare a transition
/// to change the mask.
/// Explicit prescribed energies may vary continuously with stage time.
///
/// Accuracy failures and SolverNotConverged may reduce dt. Original source,
/// invalid-data, allocation, preconditioner and cancellation failures remain
/// fatal. Accepted prefixes stay private: any late refusal returns no partial
/// window and leaves all input thermal/material data unchanged.
///
/// # Errors
/// Invalid controls/schedules, unresolved stage/event times, exhausted work,
/// temporal error at the step floor, native source/solve/cancellation refusals,
/// or a failed whole-window joule gate return no partial physical result.
#[allow(clippy::too_many_arguments)]
pub fn advance_caloric_tensor_conduction_adaptive<'f>(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    conductivity_laws: &[&dyn TensorConductivityLaw],
    storage: TensorCaloricStorage<'_>,
    start_time_s: f64,
    target_time_s: f64,
    transition_times_s: &[f64],
    config: &TensorCaloricAdaptiveConfig,
    forcing: impl FnMut(TensorCaloricForcingSample) -> Result<TensorCaloricForcing<'f>, ChtError>,
    gate: &CancelGate,
) -> Result<TensorCaloricAdaptiveSolution, TensorCaloricAdaptiveError> {
    advance_adaptive(
        domain,
        solids,
        conductivity_laws,
        storage,
        start_time_s,
        target_time_s,
        transition_times_s,
        config,
        None,
        forcing,
        gate,
    )
    .map(|value| value.thermal)
}

/// Advance original enthalpy and nonlinear radiation through an adaptive window.
///
/// Uses the SAME embedded SDIRK2(1) controller, original-energy discrepancy,
/// forcing samples, transition left limits, prescribed-mask admission and
/// compensated physical heat account as the nonradiating entry. Radiation
/// converges at both physical stages on every attempt; an accuracy rejection
/// publishes neither its state nor any of its heat.
///
/// The original ray geometry is prepared ONCE for this window and reused by
/// every accepted and rejected attempt. Geometry, emissivities and radiative
/// surroundings remain fixed. Time-varying authored sources, boundary/contact
/// and thermostat conditions use the original forcing callback and transition
/// schedule. Its read-only/deterministic/cancellation requirements also apply.
///
/// Native Newton/Krylov/backtrack and radiation.max_iterations limits are
/// shared across both stages and all radiation retries of each attempt.
/// maximum_attempts bounds their whole-window products, including incomplete
/// failed attempts; caller-owned source objects retain their own shared query
/// budgets across all attempts. Reducing dt never resets a source budget.
///
/// Only physical accepted stage powers enter radiated_energy_j. The final
/// thermal.record.radiated_w remains endpoint power; neither is another heat
/// source. The embedded estimate controls original stored energy, including
/// latent plateaus, not a certified radiation-channel quadrature error.
///
/// Cell-centered gray radiation and the original smooth-segment time-accuracy
/// limits remain. Any late error or cancellation discards the complete private
/// prefix and returns no gas, solid, phase or heat-history partial result.
///
/// # Errors
/// The original adaptive and radiating caloric admissions apply. Invalid source
/// support, controls, geometry, masks and cancellation are fatal; bounded native
/// nonconvergence or embedded accuracy rejection may reduce the next trial dt.
#[allow(clippy::too_many_arguments)]
pub fn advance_caloric_tensor_conduction_radiating_adaptive<'f>(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    conductivity_laws: &[&dyn TensorConductivityLaw],
    storage: TensorCaloricStorage<'_>,
    start_time_s: f64,
    target_time_s: f64,
    transition_times_s: &[f64],
    config: &TensorCaloricAdaptiveConfig,
    radiation: &RadiationConfig,
    forcing: impl FnMut(TensorCaloricForcingSample) -> Result<TensorCaloricForcing<'f>, ChtError>,
    gate: &CancelGate,
) -> Result<TensorCaloricRadiatingAdaptiveSolution, TensorCaloricAdaptiveError> {
    advance_adaptive(
        domain,
        solids,
        conductivity_laws,
        storage,
        start_time_s,
        target_time_s,
        transition_times_s,
        config,
        Some(radiation),
        forcing,
        gate,
    )
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn advance_adaptive<'f>(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    conductivity_laws: &[&dyn TensorConductivityLaw],
    storage: TensorCaloricStorage<'_>,
    start_time_s: f64,
    target_time_s: f64,
    transition_times_s: &[f64],
    config: &TensorCaloricAdaptiveConfig,
    radiation: Option<&RadiationConfig>,
    forcing: impl FnMut(TensorCaloricForcingSample) -> Result<TensorCaloricForcing<'f>, ChtError>,
    gate: &CancelGate,
) -> Result<TensorCaloricRadiatingAdaptiveSolution, TensorCaloricAdaptiveError> {
    advance_caloric_adaptive_with(
        domain,
        solids,
        conductivity_laws,
        storage,
        start_time_s,
        target_time_s,
        transition_times_s,
        config,
        radiation,
        forcing,
        |data, _, trial_storage, time, dt, radiation| {
            let setups = [data[0].setup(), data[1].setup()];
            let prescribed = [data[0].prescribed(), data[1].prescribed()];
            let (thermal, radiation) = match radiation {
                Some((radiation, geometry)) => advance_radiating_sdirk2_prepared(
                    domain,
                    solids,
                    setups,
                    conductivity_laws,
                    trial_storage,
                    prescribed,
                    time,
                    dt,
                    &config.step,
                    radiation,
                    Some(geometry),
                    gate,
                )
                .map(|(step, reports)| (step, Some(reports))),
                None => advance_caloric_tensor_conduction_sdirk2(
                    domain,
                    solids,
                    setups,
                    conductivity_laws,
                    trial_storage,
                    prescribed,
                    time,
                    dt,
                    &config.step,
                    gate,
                )
                .map(|step| (step, None)),
            }?;
            let mut completed_radiation_solves = 0usize;
            if let Some(reports) = &radiation {
                for stage in reports {
                    completed_radiation_solves = completed_radiation_solves
                        .checked_add(stage.iterations)
                        .ok_or(TensorCaloricAdaptiveError::Invalid(
                            "adaptive completed-trial radiation work counter overflow",
                        ))?;
                }
            }
            Ok(CaloricAdaptiveTrial {
                completed_iterations: thermal.record.iterations,
                completed_radiation_solves,
                thermal,
                radiation,
                receipt: (),
                additional_normalized_error: 0.0,
            })
        },
        gate,
    )
    .map(|run| run.thermal)
}

/// Original event, embedded-error and physical-ledger owner with a coupled
/// spatial trial seam. Only accepted small receipts are retained; any refusal
/// discards the entire private prefix, including all receipt state.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(in crate::conjugate) fn advance_caloric_adaptive_with<Forcing, Receipt>(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    conductivity_laws: &[&dyn TensorConductivityLaw],
    storage: TensorCaloricStorage<'_>,
    start_time_s: f64,
    target_time_s: f64,
    transition_times_s: &[f64],
    config: &TensorCaloricAdaptiveConfig,
    radiation: Option<&RadiationConfig>,
    mut forcing: impl FnMut(TensorCaloricForcingSample) -> Result<Forcing, ChtError>,
    mut solve: impl FnMut(
        [&Forcing; 2],
        [TensorCaloricForcingSample; 2],
        TensorCaloricStorage<'_>,
        f64,
        f64,
        Option<(&RadiationConfig, &Radiators)>,
    ) -> Result<CaloricAdaptiveTrial<Receipt>, TensorCaloricAdaptiveError>,
    gate: &CancelGate,
) -> Result<CaloricAdaptiveRun<Receipt>, TensorCaloricAdaptiveError>
where
    Forcing: CaloricAdaptiveForcing,
{
    poll(gate)?;
    validate(
        start_time_s,
        target_time_s,
        transition_times_s,
        config,
        gate,
    )?;
    let cells = domain.cell_count();
    if storage.laws.len() != cells
        || storage.reference_mass_kg.len() != cells
        || storage.previous_specific_energy_j_kg.len() != cells
        || conductivity_laws.len() != cells
    {
        return Err(TensorCaloricAdaptiveError::Invalid(
            "adaptive tensor storage needs one original law, mass, energy and tensor law per cell",
        ));
    }
    let radiators = match radiation {
        Some(config) => {
            domain.check_materials(solids.len())?;
            Some(Radiators::build(domain, solids, config, gate)?)
        }
        None => None,
    };
    let mut radiation_records = Vec::new();
    let mut radiated_energy = Sum::default();
    let mut completed_trial_radiation_solves = 0usize;
    let mut fixed = Vec::new();
    fixed.try_reserve_exact(cells).map_err(|_| {
        TensorCaloricAdaptiveError::Invalid("adaptive prescribed-cell field allocation refused")
    })?;
    fixed.resize(cells, false);
    let mut mask_count = None;
    let mut endpoint: Option<TensorCaloricEndpoint> = None;
    let mut records: Vec<TransientRecord> = Vec::new();
    let mut receipts = Vec::new();
    // Source heat, boundary heat, absolute native closure.
    let mut heat = [Sum::default(); 3];
    // Keep signed heat before independently rounding the displayed channels.
    let mut net_heat = Sum::default();
    let mut report = TensorCaloricAdaptiveReport {
        attempts: 0,
        accepted_steps: 0,
        accuracy_rejections: 0,
        solver_rejections: 0,
        completed_trial_iterations: 0,
        maximum_accepted_normalized_error: 0.0,
        absolute_accepted_closure_j: 0.0,
        suggested_dt_s: config.initial_dt_s,
    };
    let mut time = start_time_s;
    let mut event = transition_times_s.partition_point(|&next| next <= time);
    let mut accepted_iterations = 0usize;
    while time < target_time_s {
        poll(gate)?;
        if report.attempts == config.maximum_attempts {
            return Err(TensorCaloricAdaptiveError::AttemptBudget {
                time_s: time,
                attempts: report.attempts,
            });
        }
        let transition = transition_times_s.get(event).copied();
        let boundary = transition.unwrap_or(target_time_s).min(target_time_s);
        let remaining = boundary - time;
        let dt = report.suggested_dt_s.min(remaining);
        let (a, stage_times) =
            stage_clock(time, dt).ok_or(TensorCaloricAdaptiveError::TimeResolution {
                time_s: time,
                boundary_s: boundary,
            })?;
        let end = stage_times[1];
        // A numerical clock that cannot land exactly on an authored event or
        // target refuses; never label a nearby state with the requested time.
        if end > boundary || (dt == remaining && end != boundary) {
            return Err(TensorCaloricAdaptiveError::TimeResolution {
                time_s: time,
                boundary_s: boundary,
            });
        }
        let minimum = config.minimum_dt_s.min(remaining);
        report.attempts += 1;
        let samples: [TensorCaloricForcingSample; 2] =
            std::array::from_fn(|stage| TensorCaloricForcingSample {
                time_s: stage_times[stage],
                interval_s: [time, end],
                left_limit: stage == 1 && transition == Some(end),
            });
        let mut sample = |stage: usize| {
            poll(gate)?;
            let data = forcing(samples[stage])?;
            poll(gate)?;
            Ok::<_, ChtError>(data)
        };
        let first_data = sample(0)?;
        let second_data = sample(1)?;
        admit_masks(
            first_data.setup(),
            second_data.setup(),
            &mut fixed,
            mask_count,
            gate,
        )?;
        mask_count = Some(first_data.setup().fixed_temperature.len());
        let previous = endpoint
            .as_ref()
            .map_or(storage.previous_specific_energy_j_kg, |value| {
                value.specific_energy_j_kg.as_slice()
            });
        let trial_storage = TensorCaloricStorage {
            previous_specific_energy_j_kg: previous,
            ..storage
        };
        let trial = solve(
            [&first_data, &second_data],
            samples,
            trial_storage,
            time,
            dt,
            radiation.zip(radiators.as_ref()),
        );
        let trial = match trial {
            Ok(trial) => trial,
            Err(TensorCaloricAdaptiveError::Solve(error @ ChtError::SolverNotConverged { .. })) => {
                report.solver_rejections += 1;
                report.suggested_dt_s = reduced_dt(dt, minimum, f64::INFINITY)
                    .ok_or(TensorCaloricAdaptiveError::Solve(error))?;
                continue;
            }
            Err(error) => return Err(error),
        };
        let CaloricAdaptiveTrial {
            thermal: trial,
            radiation: radiation_trial,
            receipt,
            completed_iterations,
            completed_radiation_solves,
            additional_normalized_error,
        } = trial;
        completed_trial_radiation_solves = completed_trial_radiation_solves
            .checked_add(completed_radiation_solves)
            .ok_or(TensorCaloricAdaptiveError::Invalid(
                "adaptive completed-trial radiation work counter overflow",
            ))?;
        report.completed_trial_iterations = report
            .completed_trial_iterations
            .checked_add(completed_iterations)
            .ok_or(TensorCaloricAdaptiveError::Invalid(
                "adaptive completed-trial work counter overflow",
            ))?;
        if additional_normalized_error.is_nan() || additional_normalized_error < 0.0 {
            return Err(TensorCaloricAdaptiveError::Invalid(
                "invalid coupled embedded error coordinate",
            ));
        }
        let mut normalized = additional_normalized_error;
        for cell in 0..cells {
            if cell.is_multiple_of(256) {
                poll(gate)?;
            }
            if !fixed[cell] {
                let first = trial.stages[0].specific_energy_j_kg[cell] - previous[cell];
                let second = trial.stages[1].specific_energy_j_kg[cell] - previous[cell];
                normalized = normalized.max(
                    embedded_ratio(
                        first,
                        second,
                        a / dt,
                        config.energy_absolute_tolerance_j_kg,
                        config.relative_tolerance,
                    )
                    .ok_or(TensorCaloricAdaptiveError::Invalid(
                        "unrepresentable original-energy error comparison",
                    ))?,
                );
            }
        }
        if normalized > 1.0 {
            report.accuracy_rejections += 1;
            report.suggested_dt_s = reduced_dt(dt, minimum, normalized).ok_or(
                TensorCaloricAdaptiveError::TimeAccuracy {
                    time_s: time,
                    dt_s: dt,
                    normalized_error: normalized,
                },
            )?;
            continue;
        }
        poll(gate)?;
        records.try_reserve(1).map_err(|_| {
            TensorCaloricAdaptiveError::Invalid("adaptive accepted-record allocation refused")
        })?;
        receipts.try_reserve(1).map_err(|_| {
            TensorCaloricAdaptiveError::Invalid("adaptive accepted-receipt allocation refused")
        })?;
        if let Some(reports) = radiation_trial {
            radiation_records.try_reserve(1).map_err(|_| {
                TensorCaloricAdaptiveError::Invalid("adaptive radiation-record allocation refused")
            })?;
            // Match the original tableau: the first physical stage record
            // spans a=gamma*dt and its contribution is weighted by (1-gamma)/gamma.
            let carry_weight = (1.0 - TENSOR_SDIRK2_GAMMA) / TENSOR_SDIRK2_GAMMA;
            for heat_j in [
                carry_weight * (a * reports[0].radiated_w),
                a * reports[1].radiated_w,
            ] {
                radiated_energy
                    .add(heat_j)
                    .ok_or(TensorCaloricAdaptiveError::Invalid(
                        "unrepresentable accepted radiation accumulation",
                    ))?;
            }
            radiation_records.push(reports);
        }
        for (sum, value) in heat.iter_mut().zip([
            trial.record.source_j,
            trial.record.boundary_outflow_j,
            trial.record.closure_j.abs(),
        ]) {
            sum.add(value).ok_or(TensorCaloricAdaptiveError::Invalid(
                "unrepresentable accepted heat accumulation",
            ))?;
        }
        for value in [trial.record.source_j, -trial.record.boundary_outflow_j] {
            net_heat
                .add(value)
                .ok_or(TensorCaloricAdaptiveError::Invalid(
                    "unrepresentable accepted net heat accumulation",
                ))?;
        }
        accepted_iterations = accepted_iterations
            .checked_add(trial.record.iterations)
            .ok_or(TensorCaloricAdaptiveError::Invalid(
                "accepted physical-step work counter overflow",
            ))?;
        report.accepted_steps += 1;
        report.maximum_accepted_normalized_error =
            report.maximum_accepted_normalized_error.max(normalized);
        let growth = ConductionTimeOrder::First
            .accepted_growth_factor(normalized)
            .ok_or(TensorCaloricAdaptiveError::Invalid(
                "invalid embedded controller error",
            ))?;
        report.suggested_dt_s = (dt * growth).clamp(config.minimum_dt_s, config.maximum_dt_s);
        records.push(trial.record);
        receipts.push(receipt);
        let [_, accepted] = trial.stages;
        endpoint = Some(accepted);
        time = end;
        if transition == Some(time) {
            event += 1;
            mask_count = None;
        }
    }
    poll(gate)?;
    let endpoint = endpoint.ok_or(TensorCaloricAdaptiveError::Invalid(
        "positive adaptive interval produced no endpoint",
    ))?;
    let mut stored = Sum::default();
    for cell in 0..cells {
        if cell.is_multiple_of(256) {
            poll(gate)?;
        }
        stored
            .add(
                storage.reference_mass_kg[cell]
                    * (endpoint.specific_energy_j_kg[cell]
                        - storage.previous_specific_energy_j_kg[cell]),
            )
            .ok_or(TensorCaloricAdaptiveError::Invalid(
                "unrepresentable whole-window original storage",
            ))?;
    }
    let total = |sum: Sum| {
        sum.total().ok_or(TensorCaloricAdaptiveError::Invalid(
            "unrepresentable whole-window heat account",
        ))
    };
    let mut record = *records.last().ok_or(TensorCaloricAdaptiveError::Invalid(
        "adaptive endpoint has no accepted physical record",
    ))?;
    record.stored_energy_change_j = total(stored)?;
    record.source_j = total(heat[0])?;
    record.boundary_outflow_j = total(heat[1])?;
    record.closure_j = stored
        .difference(net_heat)
        .ok_or(TensorCaloricAdaptiveError::Invalid(
            "unrepresentable whole-window heat closure",
        ))?;
    finite(
        "thermal.tensor_caloric.adaptive.closure_j",
        record.closure_j,
    )?;
    record.iterations = accepted_iterations;
    report.absolute_accepted_closure_j = total(heat[2])?;
    if record.closure_j.abs() > config.energy_tolerance_j
        || report.absolute_accepted_closure_j > config.energy_tolerance_j
    {
        return Err(TensorCaloricAdaptiveError::EnergyBalance {
            closure_j: record.closure_j,
            absolute_accepted_closure_j: report.absolute_accepted_closure_j,
            tolerance_j: config.energy_tolerance_j,
        });
    }
    poll(gate)?;
    Ok(CaloricAdaptiveRun {
        thermal: TensorCaloricRadiatingAdaptiveSolution {
            thermal: TensorCaloricAdaptiveSolution {
                endpoint,
                records,
                record,
                report,
            },
            radiation: radiation_records,
            radiated_energy_j: total(radiated_energy)?,
            completed_trial_radiation_solves,
        },
        receipts,
    })
}

fn reduced_dt(dt: f64, minimum: f64, error: f64) -> Option<f64> {
    let factor = ConductionTimeOrder::First.rejected_reduction_factor(error)?;
    let reduced = (dt * factor).max(minimum);
    (reduced.is_finite() && reduced > 0.0 && reduced < dt).then_some(reduced)
}

fn admit_masks(
    first: &ThermalSetup,
    second: &ThermalSetup,
    fixed: &mut [bool],
    previous_count: Option<usize>,
    gate: &CancelGate,
) -> Result<(), TensorCaloricAdaptiveError> {
    if first.fixed_temperature.len() != second.fixed_temperature.len() {
        return Err(TensorCaloricAdaptiveError::Invalid(
            "prescribed-cell masks change inside an adaptive interval; declare a transition",
        ));
    }
    if let Some(count) = previous_count {
        if first.fixed_temperature.len() != count {
            return Err(TensorCaloricAdaptiveError::Invalid(
                "prescribed-cell masks change between attempts without a declared transition",
            ));
        }
        for (index, &(cell, _)) in first.fixed_temperature.iter().enumerate() {
            if index.is_multiple_of(4096) {
                poll(gate)?;
            }
            if !fixed.get(cell).copied().unwrap_or(false) {
                return Err(TensorCaloricAdaptiveError::Invalid(
                    "prescribed-cell masks change between attempts without a declared transition",
                ));
            }
        }
    }
    for tile in fixed.chunks_mut(4096) {
        poll(gate)?;
        tile.fill(false);
    }
    for (index, &(cell, _)) in first.fixed_temperature.iter().enumerate() {
        if index.is_multiple_of(4096) {
            poll(gate)?;
        }
        let entry = fixed
            .get_mut(cell)
            .ok_or(TensorCaloricAdaptiveError::Invalid(
                "prescribed cell is outside the tensor domain",
            ))?;
        if *entry {
            return Err(TensorCaloricAdaptiveError::Invalid(
                "prescribed cells must be unique",
            ));
        }
        *entry = true;
    }
    for (index, &(cell, _)) in second.fixed_temperature.iter().enumerate() {
        if index.is_multiple_of(4096) {
            poll(gate)?;
        }
        let entry = fixed
            .get_mut(cell)
            .ok_or(TensorCaloricAdaptiveError::Invalid(
                "prescribed cell is outside the tensor domain",
            ))?;
        if !*entry {
            return Err(TensorCaloricAdaptiveError::Invalid(
                "prescribed-cell masks change or duplicate inside an adaptive interval",
            ));
        }
        *entry = false;
    }
    for (index, &(cell, _)) in first.fixed_temperature.iter().enumerate() {
        if index.is_multiple_of(4096) {
            poll(gate)?;
        }
        fixed[cell] = true;
    }
    Ok(())
}

fn validate(
    start: f64,
    target: f64,
    transitions: &[f64],
    config: &TensorCaloricAdaptiveConfig,
    gate: &CancelGate,
) -> Result<(), TensorCaloricAdaptiveError> {
    let positive = |x: f64| x.is_finite() && x > 0.0;
    if !start.is_finite()
        || start < 0.0
        || !target.is_finite()
        || target <= start
        || !(target - start).is_finite()
        || !positive(config.minimum_dt_s)
        || !positive(config.initial_dt_s)
        || !positive(config.maximum_dt_s)
        || config.minimum_dt_s > config.initial_dt_s
        || config.initial_dt_s > config.maximum_dt_s
        || !positive(config.energy_absolute_tolerance_j_kg)
        || !config.relative_tolerance.is_finite()
        || !(0.0..1.0).contains(&config.relative_tolerance)
        || !config.energy_tolerance_j.is_finite()
        || config.energy_tolerance_j < 0.0
    {
        return Err(TensorCaloricAdaptiveError::Invalid(
            "invalid adaptive target, step bounds or original-energy tolerances",
        ));
    }
    for cap in [
        2, // Forcing callback evaluations.
        config.step.max_newton_iterations,
        config.step.max_backtracks,
        config.step.linear.max_iterations,
    ] {
        config
            .maximum_attempts
            .checked_mul(cap)
            .ok_or(TensorCaloricAdaptiveError::Invalid(
                "adaptive work bound is unrepresentable",
            ))?;
    }
    config
        .maximum_attempts
        .checked_mul(std::mem::size_of::<TransientRecord>())
        .filter(|&bytes| isize::try_from(bytes).is_ok())
        .ok_or(TensorCaloricAdaptiveError::Invalid(
            "adaptive accepted-record bound exceeds addressable storage",
        ))?;
    let mut previous = None;
    for (index, &next) in transitions.iter().enumerate() {
        if index.is_multiple_of(4096) {
            poll(gate)?;
        }
        if !next.is_finite() || next < 0.0 || previous.is_some_and(|old| next <= old) {
            return Err(TensorCaloricAdaptiveError::Invalid(
                "transition times must be finite, nonnegative and strictly increasing",
            ));
        }
        previous = Some(next);
    }
    poll(gate)?;
    Ok(())
}
