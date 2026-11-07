//! Reacting volume histories with an original bath or finite caloric body.
//!
//! Knots share the adiabatic history admission. Each segment delegates to the
//! prepared thermal compression owner and its fixed/adaptive native stages.
//! Original gas composition and finite-body enthalpy/phase continue through every
//! turn and dwell; neither physical endpoint is published before the whole
//! history succeeds. Only reusable numerical workspace changes on refusal.

use super::super::thermal::adaptive::{
    AdaptiveMaterialThermalCompressionError, AdaptiveMaterialThermalCompressionSolve,
    AdaptiveMaterialThermalCompressionStep,
};
use super::super::thermal::{
    MaterialThermalCompressionBodyStep, MaterialThermalCompressionReactor,
    MaterialThermalCompressionSolve, MaterialThermalCompressionStep,
};
use super::{
    ClosedGasParcel, CompressionKnot, CompressionTotals, MAX_COMPRESSION_EVALUATIONS,
    MAX_COMPRESSION_SUBSTEPS, MaterialCompressionError, MaterialCompressionHistoryError,
    accumulate, add_count, poll, segment_error, validate_knots,
};
use crate::chemistry::heat::{HeatBath, material::FiniteHeatBody};
use fs_math::dd::Dd;
use fs_qty::{Time, Volume};

type Error = MaterialCompressionHistoryError;

/// Original thermal controls with work allowances shared by every segment.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MaterialThermalCompressionHistorySolve {
    /// Moving segments use the original fixed subdivision; dwells cost one
    /// native interval. Both fixed bounds are in 1..=4096.
    Fixed {
        controls: MaterialThermalCompressionSolve,
        maximum_substeps: usize,
    },
    /// Native full/two-half comparisons. Maximum intervals/evaluations apply
    /// globally, including coarse and rejected probes; dwells cost one each.
    Adaptive(AdaptiveMaterialThermalCompressionSolve),
}

/// Complete accepted history with original source heat, work and entropy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialThermalCompressionHistoryStep<P = MaterialThermalCompressionStep> {
    /// Compensated original balances. The finite-body payload also retains the
    /// accepted body, original enthalpy/internal-energy/expansion-work changes
    /// and separately accumulated coupled residuals.
    pub physical: P,
    pub completed_segments: usize,
    /// All native intervals, including coarse/rejected adaptive probes.
    pub evaluations: usize,
    /// Native accepted comparisons plus dwells; zero in fixed mode.
    pub accepted_adaptive_intervals: usize,
    pub rejected_adaptive_intervals: usize,
    /// Native temperature/species/work/heat/body-enthalpy maxima; None in fixed
    /// mode. These local comparisons are not a global accuracy certificate.
    pub maximum_accepted_error_ratios: Option<[f64; 5]>,
    pub suggested_step_s: Option<f64>,
    /// Chemical work on all probes; physical counts include retained paths only.
    pub total_chemical_accepted_steps: usize,
    pub total_chemical_rejected_steps: usize,
    /// Compensated segment durations minus the exact final elapsed knot [s].
    /// Exposes absolute-clock subtraction rounding, not integration error.
    pub elapsed_time_error_s: f64,
}

/// A gas/body history publishes one accepted original body with its gas balances.
pub type MaterialThermalCompressionHistoryBodyStep<'m> =
    MaterialThermalCompressionHistoryStep<MaterialThermalCompressionBodyStep<'m>>;

impl<P> MaterialThermalCompressionHistoryStep<P> {
    fn with_physical<Q>(self, physical: Q) -> MaterialThermalCompressionHistoryStep<Q> {
        MaterialThermalCompressionHistoryStep {
            physical,
            completed_segments: self.completed_segments,
            evaluations: self.evaluations,
            accepted_adaptive_intervals: self.accepted_adaptive_intervals,
            rejected_adaptive_intervals: self.rejected_adaptive_intervals,
            maximum_accepted_error_ratios: self.maximum_accepted_error_ratios,
            suggested_step_s: self.suggested_step_s,
            total_chemical_accepted_steps: self.total_chemical_accepted_steps,
            total_chemical_rejected_steps: self.total_chemical_rejected_steps,
            elapsed_time_error_s: self.elapsed_time_error_s,
        }
    }
}

impl MaterialThermalCompressionReactor {
    /// Follow every prescribed turn and dwell with a held finite-conductance bath.
    ///
    /// The original coupled chemistry/heat stage owns reaction and signed
    /// G(Tbath-Tgas). Pressure work remains the native volume-map account.
    /// Equal first/last volumes never bypass the intermediate trajectory.
    /// Knots, global work caps, original domains and the final cancellation
    /// checkpoint must all pass before publishing the complete gas parcel.
    pub fn advance_history_with_heat_bath(
        &mut self,
        parcel: &mut ClosedGasParcel,
        knots: &[CompressionKnot],
        bath: HeatBath,
        solve: MaterialThermalCompressionHistorySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialThermalCompressionHistoryStep, Error> {
        let (candidate, report) = propose_history(
            self, parcel, knots, solve, &mut cancelled,
            |reactor, candidate, volume, duration, solve, cancelled| match solve {
                MaterialThermalCompressionHistorySolve::Fixed { controls, .. } => {
                    let physical = reactor.advance_with_heat_bath(
                        candidate, volume, duration, bath, controls, cancelled,
                    ).map_err(NativeError::Fixed)?;
                    Ok((physical, SegmentStatistics::fixed(&physical)))
                }
                MaterialThermalCompressionHistorySolve::Adaptive(controls) => {
                    let report = reactor.advance_with_heat_bath_adaptive(
                        candidate, volume, duration, bath, controls, cancelled,
                    ).map_err(NativeError::Adaptive)?;
                    let statistics = SegmentStatistics::adaptive(&report);
                    Ok((report.physical, statistics))
                }
            },
        )?;
        poll(&mut cancelled)?;
        *parcel = candidate;
        Ok(report)
    }

    /// Retain one original finite body's enthalpy, phase and source throughout
    /// compression, expansion and physical-time dwells.
    ///
    /// The body is copied by value as in the original endpoint API. Gas and
    /// returned body publish only after the entire history succeeds. A refusal
    /// after accepted internal segments leaves both caller endpoints intact.
    /// Source mass, pressure, curve and conductance remain original; gas volume
    /// follows the knots and never the body's density-derived volume.
    ///
    /// Heat cancels once in combined balances. Body outward pressure work stays
    /// separate from work on the gas, with no latent/reaction heat counted twice.
    pub fn advance_history_with_finite_body<'m>(
        &mut self,
        parcel: &mut ClosedGasParcel,
        knots: &[CompressionKnot],
        mut body: FiniteHeatBody<'m>,
        solve: MaterialThermalCompressionHistorySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialThermalCompressionHistoryBodyStep<'m>, Error> {
        let mut body_totals = [(0.0, 0.0); 5];
        let (candidate, gas) = propose_history(
            self, parcel, knots, solve, &mut cancelled,
            |reactor, candidate, volume, duration, solve, cancelled| {
                let (physical, statistics) = match solve {
                    MaterialThermalCompressionHistorySolve::Fixed { controls, .. } => {
                        let physical = reactor.advance_with_finite_body(
                            candidate, volume, duration, body, controls, cancelled,
                        ).map_err(NativeError::Fixed)?;
                        let statistics = SegmentStatistics::fixed(&physical.gas);
                        (physical, statistics)
                    }
                    MaterialThermalCompressionHistorySolve::Adaptive(controls) => {
                        let report = reactor.advance_with_finite_body_adaptive(
                            candidate, volume, duration, body, controls, cancelled,
                        ).map_err(NativeError::Adaptive)?;
                        let statistics = SegmentStatistics::adaptive(&report);
                        (report.physical, statistics)
                    }
                };
                for (total, value) in body_totals.iter_mut().zip([
                    physical.body_enthalpy_change_j,
                    physical.body_internal_energy_change_j,
                    physical.body_boundary_work_j,
                    physical.combined_heat_residual_j,
                    physical.combined_energy_residual_j,
                ]) {
                    accumulate(total, value);
                }
                body = physical.body;
                Ok((physical.gas, statistics))
            },
        )?;
        if body_totals.iter().any(|total| !total.0.is_finite()) {
            return Err(Error::Invalid(
                "whole-history finite-body work, storage or residual is unrepresentable",
            ));
        }
        let physical = MaterialThermalCompressionBodyStep {
            gas: gas.physical,
            body,
            body_enthalpy_change_j: body_totals[0].0,
            body_internal_energy_change_j: body_totals[1].0,
            body_boundary_work_j: body_totals[2].0,
            body_entropy_change_j_k: gas.physical.thermal_partner_entropy_change_j_k,
            combined_heat_residual_j: body_totals[3].0,
            combined_energy_residual_j: body_totals[4].0,
        };
        let report = gas.with_physical(physical);
        poll(&mut cancelled)?;
        *parcel = candidate;
        Ok(report)
    }
}

struct SegmentStatistics {
    evaluations: usize,
    accepted: usize,
    rejected: usize,
    errors: Option<[f64; 5]>,
    suggested_step: Option<f64>,
    chemical: [usize; 2],
}

impl SegmentStatistics {
    fn fixed(physical: &MaterialThermalCompressionStep) -> Self {
        Self {
            evaluations: physical.substeps,
            accepted: 0,
            rejected: 0,
            errors: None,
            suggested_step: None,
            chemical: [physical.chemical_accepted_steps, physical.chemical_rejected_steps],
        }
    }

    fn adaptive<P>(step: &AdaptiveMaterialThermalCompressionStep<P>) -> Self {
        Self {
            evaluations: step.evaluations,
            accepted: step.accepted_intervals,
            rejected: step.rejected_intervals,
            errors: Some(step.maximum_accepted_error_ratios),
            suggested_step: Some(step.suggested_step_s),
            chemical: [step.total_chemical_accepted_steps, step.total_chemical_rejected_steps],
        }
    }
}

enum NativeError {
    Fixed(MaterialCompressionError),
    Adaptive(AdaptiveMaterialThermalCompressionError),
}

fn adaptive_error(
    segment: usize,
    start_time_s: f64,
    end_time_s: f64,
    prior: [usize; 3],
    source: AdaptiveMaterialThermalCompressionError,
) -> Error {
    use AdaptiveMaterialThermalCompressionError as Native;
    let source = match source {
        Native::Cancelled => return Error::Cancelled,
        Native::Accuracy { elapsed_time_s, step_s, error_ratios } => Native::Accuracy {
            elapsed_time_s: start_time_s + elapsed_time_s, step_s, error_ratios,
        },
        Native::WorkLimit { elapsed_time_s, accepted_intervals, rejected_intervals, evaluations } => Native::WorkLimit {
            elapsed_time_s: start_time_s + elapsed_time_s,
            accepted_intervals: prior[0] + accepted_intervals,
            rejected_intervals: prior[1] + rejected_intervals,
            evaluations: prior[2] + evaluations,
        },
        Native::TimeResolution { elapsed_time_s, target_time_s } => Native::TimeResolution {
            elapsed_time_s: start_time_s + elapsed_time_s,
            target_time_s: start_time_s + target_time_s,
        },
        other => other,
    };
    Error::ThermalAdaptive { segment, start_time_s, end_time_s, source }
}

fn validate_controls(solve: MaterialThermalCompressionHistorySolve) -> Result<(), Error> {
    match solve {
        MaterialThermalCompressionHistorySolve::Fixed { controls, maximum_substeps } => {
            if !(1..=MAX_COMPRESSION_SUBSTEPS).contains(&controls.substeps)
                || !(1..=MAX_COMPRESSION_SUBSTEPS).contains(&maximum_substeps)
            {
                return Err(Error::Invalid(
                    "fixed thermal per-segment and whole-history substeps must be in 1..=4096",
                ));
            }
        }
        MaterialThermalCompressionHistorySolve::Adaptive(controls) => {
            if !(1..=MAX_COMPRESSION_SUBSTEPS).contains(&controls.maximum_intervals)
                || !(3..=MAX_COMPRESSION_EVALUATIONS).contains(&controls.maximum_evaluations)
            {
                return Err(Error::Invalid(
                    "thermal adaptive whole-history intervals/evaluations exceed original bounds",
                ));
            }
        }
    }
    Ok(())
}

/// Only original endpoint solvers evaluate physics; this loop shares schedule,
/// work accounting and atomic publication for the two thermal partners.
#[allow(clippy::too_many_lines)]
fn propose_history<F: FnMut() -> bool>(
    reactor: &mut MaterialThermalCompressionReactor,
    parcel: &ClosedGasParcel,
    knots: &[CompressionKnot],
    solve: MaterialThermalCompressionHistorySolve,
    cancelled: &mut F,
    mut advance: impl FnMut(
        &mut MaterialThermalCompressionReactor,
        &mut ClosedGasParcel,
        Volume,
        Time,
        MaterialThermalCompressionHistorySolve,
        &mut F,
    ) -> Result<(MaterialThermalCompressionStep, SegmentStatistics), NativeError>,
) -> Result<(ClosedGasParcel, MaterialThermalCompressionHistoryStep), Error> {
    validate_knots(parcel, knots, cancelled)?;
    validate_controls(solve)?;
    let initial_volume = parcel.volume_m3();
    let final_knot = knots[knots.len() - 1];
    let mut candidate = parcel.clone();
    let mut totals = CompressionTotals::default();
    let mut thermal = [(0.0, 0.0); 3];
    let mut elapsed = Dd::ZERO;
    let (mut evaluations, mut accepted, mut rejected) = (0, 0, 0);
    let mut chemical = [0_usize; 2];
    let mut maxima = match solve {
        MaterialThermalCompressionHistorySolve::Fixed { .. } => None,
        MaterialThermalCompressionHistorySolve::Adaptive(_) => Some([0.0_f64; 5]),
    };
    let mut suggested_step = None;
    for (segment, pair) in knots.windows(2).enumerate() {
        poll(cancelled)?;
        let start = pair[0].elapsed_time.value();
        let end = pair[1].elapsed_time.value();
        let moving = pair[1].volume.value() != candidate.volume_m3();
        let work_limit = Error::WorkLimit {
            completed_segments: segment, elapsed_time_s: start, evaluations,
        };
        let local = match solve {
            MaterialThermalCompressionHistorySolve::Fixed { controls, maximum_substeps } => {
                let needed = if moving { controls.substeps } else { 1 };
                if needed > maximum_substeps - evaluations {
                    return Err(work_limit);
                }
                solve
            }
            MaterialThermalCompressionHistorySolve::Adaptive(mut controls) => {
                let needed = if moving { 3 } else { 1 };
                if accepted == controls.maximum_intervals
                    || needed > controls.maximum_evaluations - evaluations
                {
                    return Err(work_limit);
                }
                controls.maximum_intervals -= accepted;
                controls.maximum_evaluations -= evaluations;
                // A known dwell calls the native rigid path exactly once. It
                // may use the final 1/2 evaluations despite native admission3.
                if !moving {
                    controls.maximum_evaluations = controls.maximum_evaluations.max(3);
                }
                if let Some(step) = suggested_step {
                    controls.initial_step_s = f64::clamp(
                        step, controls.minimum_step_s, controls.maximum_step_s,
                    );
                }
                MaterialThermalCompressionHistorySolve::Adaptive(controls)
            }
        };
        let (physical, statistics) = advance(
            reactor, &mut candidate, pair[1].volume, Time::new(end - start), local, cancelled,
        ).map_err(|source| match source {
            NativeError::Fixed(source) => segment_error(segment, start, end, source),
            NativeError::Adaptive(source) => adaptive_error(
                segment, start, end, [accepted, rejected, evaluations], source,
            ),
        })?;
        evaluations += statistics.evaluations;
        accepted += statistics.accepted;
        rejected += statistics.rejected;
        for (total, next) in chemical.iter_mut().zip(statistics.chemical) {
            add_count(total, next)?;
        }
        if let (Some(maxima), Some(errors)) = (&mut maxima, statistics.errors) {
            for (maximum, error) in maxima.iter_mut().zip(errors) {
                *maximum = (*maximum).max(error);
            }
        }
        suggested_step = statistics.suggested_step;
        accumulate(&mut totals.work, physical.boundary_work_on_gas_j);
        accumulate(&mut totals.energy, physical.internal_energy_change_j);
        accumulate(&mut totals.residual, physical.energy_residual_j);
        accumulate(&mut totals.entropy, physical.gas_entropy_change_j_k);
        totals.add_counts(
            physical.substeps, physical.chemical_accepted_steps, physical.chemical_rejected_steps,
        ).map_err(|source| segment_error(segment, start, end, source))?;
        for (total, value) in thermal.iter_mut().zip([
            physical.heat_input_j,
            physical.thermal_partner_entropy_change_j_k,
            physical.entropy_production_j_k,
        ]) {
            accumulate(total, value);
        }
        elapsed = elapsed + Dd::from_f64(physical.elapsed_time.value());
    }
    let gas = totals.finish(initial_volume, final_knot.volume.value(), final_knot.elapsed_time)
        .map_err(|_| Error::Invalid("whole-history thermal gas work/energy/entropy is unrepresentable"))?;
    let elapsed_time_error_s = (elapsed - Dd::from_f64(final_knot.elapsed_time.value())).to_f64();
    if !elapsed_time_error_s.is_finite() || thermal.iter().any(|total| !total.0.is_finite()) {
        return Err(Error::Invalid("whole-history thermal heat/entropy/time is unrepresentable"));
    }
    Ok((
        candidate,
        MaterialThermalCompressionHistoryStep {
            physical: MaterialThermalCompressionStep {
                elapsed_time: final_knot.elapsed_time,
                initial_volume_m3: initial_volume,
                final_volume_m3: final_knot.volume.value(),
                boundary_work_on_gas_j: gas.boundary_work_on_gas_j,
                heat_input_j: thermal[0].0,
                internal_energy_change_j: gas.internal_energy_change_j,
                energy_residual_j: gas.energy_residual_j,
                gas_entropy_change_j_k: gas.entropy_change_j_k,
                thermal_partner_entropy_change_j_k: thermal[1].0,
                entropy_production_j_k: thermal[2].0,
                substeps: gas.substeps,
                chemical_accepted_steps: gas.chemical_accepted_steps,
                chemical_rejected_steps: gas.chemical_rejected_steps,
            },
            completed_segments: knots.len() - 1,
            evaluations,
            accepted_adaptive_intervals: accepted,
            rejected_adaptive_intervals: rejected,
            maximum_accepted_error_ratios: maxima,
            suggested_step_s: suggested_step,
            total_chemical_accepted_steps: chemical[0],
            total_chemical_rejected_steps: chemical[1],
            elapsed_time_error_s,
        },
    ))
}
