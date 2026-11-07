//! Prescribed reacting compression histories with turns, dwells and closed cycles.
//!
//! Every adjacent pair of elapsed-time/volume knots is one original stroke;
//! volume is log-linear within that segment. Matching the first and last volume
//! never collapses a cycle into a rigid reactor. The same candidate inventory
//! passes through every segment, with one atomic publication at the end.
//!
//! This module orchestrates the existing fixed/adaptive owners. It does not
//! duplicate kinetics, caloric inversion, volume work or error comparisons.
//! Adaptive local comparisons are not a global accuracy certificate.

/// Bath and original finite-body heat exchange throughout prescribed histories.
pub mod thermal;

use super::adaptive::{
    AdaptiveMaterialCompressionError, AdaptiveMaterialCompressionSolve,
    MAX_COMPRESSION_EVALUATIONS, advance_closed_parcel_compression_adaptive,
};
use super::{
    ClosedGasParcel, CompressionTotals, MAX_COMPRESSION_SUBSTEPS, MaterialCompressionError,
    MaterialCompressionSolve, MaterialCompressionStep, ThermochemicalGasModel, accumulate,
    advance_closed_parcel_compression, positive,
};
use fs_math::dd::Dd;
use fs_qty::{Time, Volume};
use fs_thermochem::kinetics::IdealGasKineticsV1;

/// At most 4096 finite positive-duration segments in one atomic update.
pub const MAX_COMPRESSION_HISTORY_KNOTS: usize = MAX_COMPRESSION_SUBSTEPS + 1;

/// Prescribed physical volume at an elapsed time measured from this update.
///
/// The first knot must be at zero and at the parcel's exact current volume.
/// Times then increase strictly; equal adjacent volumes explicitly declare a
/// dwell. Intermediate volumes and turning points are never removed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompressionKnot {
    pub elapsed_time: Time,
    pub volume: Volume,
}

/// Original physical controls, with work budgets applying to the WHOLE history.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MaterialCompressionHistorySolve {
    /// `controls.substeps` per moving segment; one native interval per dwell.
    /// `maximum_substeps` bounds their total, in 1..=4096.
    Fixed {
        controls: MaterialCompressionSolve,
        maximum_substeps: usize,
    },
    /// Original comparison tolerances and step limits. The original
    /// `maximum_intervals` and `maximum_evaluations` are global history caps,
    /// not fresh allowances at every knot. Dwells use one interval/evaluation.
    Adaptive(AdaptiveMaterialCompressionSolve),
}

/// Original physical balances and actual work of the entire accepted history.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialCompressionHistoryStep {
    /// Compensated sum of accepted source work, storage, residual and entropy.
    /// `substeps` counts actual retained physical intervals, including both fine
    /// halves of an accepted adaptive comparison and one interval per dwell.
    pub physical: MaterialCompressionStep,
    pub completed_segments: usize,
    /// All original interval evaluations, including coarse/rejected probes.
    pub evaluations: usize,
    /// Accepted adaptive comparisons or dwell intervals; zero in fixed mode.
    pub accepted_adaptive_intervals: usize,
    pub rejected_adaptive_intervals: usize,
    /// Original temperature/species/work comparison maxima; None in fixed mode.
    pub maximum_accepted_error_ratios: Option<[f64; 3]>,
    /// Final native step suggestion; retained across intervening dwells.
    /// None in fixed mode.
    pub suggested_step_s: Option<f64>,
    /// Chemical work across every successful probe, including rejected/coarse
    /// trajectories. Physical counts above retain accepted trajectories only.
    pub total_chemical_accepted_steps: usize,
    pub total_chemical_rejected_steps: usize,
    /// Compensated sum of segment-reported durations minus the final elapsed
    /// knot [s]. This exposes subtraction rounding between absolute knots; it
    /// is not a bound on chemical or splitting time-integration error.
    pub elapsed_time_error_s: f64,
}

/// Every refusal leaves the complete caller parcel and its source receipts intact.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialCompressionHistoryError {
    Invalid(&'static str),
    Knot {
        index: usize,
        requirement: &'static str,
    },
    /// No further original interval can be started under the global allowance.
    WorkLimit {
        completed_segments: usize,
        elapsed_time_s: f64,
        evaluations: usize,
    },
    /// Fixed/native physical refusal; inner owners do not report a partial clock.
    Segment {
        segment: usize,
        start_time_s: f64,
        end_time_s: f64,
        source: MaterialCompressionError,
    },
    /// Native accuracy/time/work errors use the GLOBAL history clock and
    /// cumulative interval/evaluation counts. Other source errors retain the
    /// enclosing segment bounds without inventing partial progress.
    Adaptive {
        segment: usize,
        start_time_s: f64,
        end_time_s: f64,
        source: AdaptiveMaterialCompressionError,
    },
    /// Thermal native comparisons use the same global clock/count convention,
    /// including original heat and finite-body enthalpy error components.
    ThermalAdaptive {
        segment: usize,
        start_time_s: f64,
        end_time_s: f64,
        source: super::thermal::adaptive::AdaptiveMaterialThermalCompressionError,
    },
    Cancelled,
}
impl core::fmt::Display for MaterialCompressionHistoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reacting material compression history: {self:?}")
    }
}
impl std::error::Error for MaterialCompressionHistoryError {}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MaterialCompressionHistoryError> {
    if cancelled() {
        Err(MaterialCompressionHistoryError::Cancelled)
    } else {
        Ok(())
    }
}

fn validate(
    parcel: &ClosedGasParcel,
    knots: &[CompressionKnot],
    solve: MaterialCompressionHistorySolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), MaterialCompressionHistoryError> {
    validate_knots(parcel, knots, cancelled)?;
    match solve {
        MaterialCompressionHistorySolve::Fixed {
            controls,
            maximum_substeps,
        } => {
            if !(1..=MAX_COMPRESSION_SUBSTEPS).contains(&maximum_substeps)
                || !(1..=MAX_COMPRESSION_SUBSTEPS).contains(&controls.substeps)
            {
                return Err(MaterialCompressionHistoryError::Invalid(
                    "fixed per-segment and whole-history substeps must be in 1..=4096",
                ));
            }
        }
        MaterialCompressionHistorySolve::Adaptive(controls) => {
            if !(1..=MAX_COMPRESSION_SUBSTEPS).contains(&controls.maximum_intervals)
                || !(3..=MAX_COMPRESSION_EVALUATIONS).contains(&controls.maximum_evaluations)
            {
                return Err(MaterialCompressionHistoryError::Invalid(
                    "adaptive whole-history intervals/evaluations exceed original bounds",
                ));
            }
        }
    }
    Ok(())
}

/// Common schedule admission for adiabatic and thermal original-source owners.
fn validate_knots(
    parcel: &ClosedGasParcel,
    knots: &[CompressionKnot],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), MaterialCompressionHistoryError> {
    poll(cancelled)?;
    if !(2..=MAX_COMPRESSION_HISTORY_KNOTS).contains(&knots.len()) {
        return Err(MaterialCompressionHistoryError::Invalid(
            "history requires 2..=4097 knots",
        ));
    }
    for (index, knot) in knots.iter().enumerate() {
        poll(cancelled)?;
        let time = knot.elapsed_time.value();
        if !time.is_finite() || time < 0.0 || !positive(knot.volume.value()) {
            return Err(MaterialCompressionHistoryError::Knot {
                index,
                requirement: "finite nonnegative elapsed seconds and positive finite volume required",
            });
        }
        if index == 0 {
            if time != 0.0 || knot.volume.value().to_bits() != parcel.volume_m3().to_bits() {
                return Err(MaterialCompressionHistoryError::Knot {
                    index,
                    requirement: "first knot must bind time zero and the exact current parcel volume",
                });
            }
        } else if time <= knots[index - 1].elapsed_time.value()
            || !positive(knot.volume.value() / knots[index - 1].volume.value())
        {
            return Err(MaterialCompressionHistoryError::Knot {
                index,
                requirement: "elapsed times must increase strictly and volume ratios remain representable",
            });
        }
    }
    Ok(())
}

fn segment_error(
    segment: usize,
    start_time_s: f64,
    end_time_s: f64,
    source: MaterialCompressionError,
) -> MaterialCompressionHistoryError {
    match source {
        MaterialCompressionError::Cancelled => MaterialCompressionHistoryError::Cancelled,
        source => MaterialCompressionHistoryError::Segment {
            segment,
            start_time_s,
            end_time_s,
            source,
        },
    }
}

fn adaptive_error(
    segment: usize,
    start_time_s: f64,
    end_time_s: f64,
    prior: [usize; 3],
    source: AdaptiveMaterialCompressionError,
) -> MaterialCompressionHistoryError {
    let source = match source {
        AdaptiveMaterialCompressionError::Cancelled => {
            return MaterialCompressionHistoryError::Cancelled;
        }
        AdaptiveMaterialCompressionError::Accuracy {
            elapsed_time_s,
            step_s,
            error_ratios,
        } => AdaptiveMaterialCompressionError::Accuracy {
            elapsed_time_s: start_time_s + elapsed_time_s,
            step_s,
            error_ratios,
        },
        AdaptiveMaterialCompressionError::WorkLimit {
            elapsed_time_s,
            accepted_intervals,
            rejected_intervals,
            evaluations,
        } => AdaptiveMaterialCompressionError::WorkLimit {
            elapsed_time_s: start_time_s + elapsed_time_s,
            accepted_intervals: prior[0] + accepted_intervals,
            rejected_intervals: prior[1] + rejected_intervals,
            evaluations: prior[2] + evaluations,
        },
        AdaptiveMaterialCompressionError::TimeResolution {
            elapsed_time_s,
            target_time_s,
        } => AdaptiveMaterialCompressionError::TimeResolution {
            elapsed_time_s: start_time_s + elapsed_time_s,
            target_time_s: start_time_s + target_time_s,
        },
        other => other,
    };
    MaterialCompressionHistoryError::Adaptive {
        segment,
        start_time_s,
        end_time_s,
        source,
    }
}

fn add_count(total: &mut usize, next: usize) -> Result<(), MaterialCompressionHistoryError> {
    *total = total
        .checked_add(next)
        .ok_or(MaterialCompressionHistoryError::Invalid(
            "whole-history chemical work count overflow",
        ))?;
    Ok(())
}

/// Advance a sealed reacting parcel through every prescribed segment atomically.
///
/// Volume is log-linear between each consecutive pair of knots; a constant
/// adjacent volume is a physical-time dwell using the original rigid reactor.
/// Turns and dwells remain explicit even if the final volume equals the first.
/// The original mass, species bank, kinetics and material sources remain the
/// owners throughout; product inventories are never reset at a knot.
///
/// All knots are admitted before evolution. Fixed work and adaptive comparison
/// budgets apply globally, including rejected/coarse probes. The previous
/// adaptive step suggestion seeds the next segment, subject to the same native
/// step limits and error comparisons. Knots are exact integration boundaries,
/// so a comparison never crosses a prescribed rate discontinuity.
///
/// Work, energy, native reconciliation residual and stage entropy are reduced
/// in chronological compensated sums. The entropy and local-error limitations
/// of the original owners remain unchanged. No piston inertia, heat exchange,
/// reaction-rate calibration, extrapolation or global error certificate is added.
#[allow(clippy::too_many_lines)]
pub fn advance_closed_parcel_compression_history(
    parcel: &mut ClosedGasParcel,
    bank: &ThermochemicalGasModel,
    mechanism: &IdealGasKineticsV1,
    knots: &[CompressionKnot],
    solve: MaterialCompressionHistorySolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MaterialCompressionHistoryStep, MaterialCompressionHistoryError> {
    validate(parcel, knots, solve, &mut cancelled)?;
    let initial_volume = parcel.volume_m3();
    let final_knot = knots[knots.len() - 1];
    let mut candidate = parcel.clone();
    let mut totals = CompressionTotals::default();
    let mut elapsed = Dd::ZERO;
    let (mut evaluations, mut accepted, mut rejected) = (0, 0, 0);
    let mut chemical = [0_usize; 2];
    let mut maximum_errors = match solve {
        MaterialCompressionHistorySolve::Fixed { .. } => None,
        MaterialCompressionHistorySolve::Adaptive(_) => Some([0.0_f64; 3]),
    };
    let mut suggested_step = None;
    for (segment, pair) in knots.windows(2).enumerate() {
        poll(&mut cancelled)?;
        let start = pair[0].elapsed_time.value();
        let end = pair[1].elapsed_time.value();
        let duration = Time::new(end - start);
        let moving = pair[1].volume.value() != candidate.volume_m3();
        let work_limit = || MaterialCompressionHistoryError::WorkLimit {
            completed_segments: segment,
            elapsed_time_s: start,
            evaluations,
        };
        let physical = match solve {
            MaterialCompressionHistorySolve::Fixed {
                controls,
                maximum_substeps,
            } => {
                let needed = if moving { controls.substeps } else { 1 };
                if needed > maximum_substeps - evaluations {
                    return Err(work_limit());
                }
                let physical = advance_closed_parcel_compression(
                    &mut candidate,
                    bank,
                    mechanism,
                    pair[1].volume,
                    duration,
                    controls,
                    &mut cancelled,
                )
                .map_err(|source| segment_error(segment, start, end, source))?;
                evaluations += physical.substeps;
                add_count(&mut chemical[0], physical.chemical_accepted_steps)?;
                add_count(&mut chemical[1], physical.chemical_rejected_steps)?;
                physical
            }
            MaterialCompressionHistorySolve::Adaptive(mut controls) => {
                let needed = if moving { 3 } else { 1 };
                if accepted == controls.maximum_intervals
                    || needed > controls.maximum_evaluations - evaluations
                {
                    return Err(work_limit());
                }
                controls.maximum_intervals -= accepted;
                controls.maximum_evaluations -= evaluations;
                // Native adaptive admission requires >=3 even on its exactly
                // one-evaluation rigid path. A known dwell may consume the last
                // one or two global evaluations; moving strokes never get this.
                if !moving {
                    controls.maximum_evaluations = controls.maximum_evaluations.max(3);
                }
                if let Some(step) = suggested_step {
                    controls.initial_step_s =
                        f64::clamp(step, controls.minimum_step_s, controls.maximum_step_s);
                }
                let step = advance_closed_parcel_compression_adaptive(
                    &mut candidate,
                    bank,
                    mechanism,
                    pair[1].volume,
                    duration,
                    controls,
                    &mut cancelled,
                )
                .map_err(|source| {
                    adaptive_error(
                        segment,
                        start,
                        end,
                        [accepted, rejected, evaluations],
                        source,
                    )
                })?;
                accepted += step.accepted_intervals;
                rejected += step.rejected_intervals;
                evaluations += step.evaluations;
                add_count(&mut chemical[0], step.total_chemical_accepted_steps)?;
                add_count(&mut chemical[1], step.total_chemical_rejected_steps)?;
                if let Some(maxima) = &mut maximum_errors {
                    for (maximum, error) in
                        maxima.iter_mut().zip(step.maximum_accepted_error_ratios)
                    {
                        *maximum = (*maximum).max(error);
                    }
                }
                suggested_step = Some(step.suggested_step_s);
                step.physical
            }
        };
        // These are the original stage balances, not differences of absolute
        // stored energies or accumulated work. Small residuals stay separate.
        accumulate(&mut totals.work, physical.boundary_work_on_gas_j);
        accumulate(&mut totals.energy, physical.internal_energy_change_j);
        accumulate(&mut totals.residual, physical.energy_residual_j);
        accumulate(&mut totals.entropy, physical.entropy_change_j_k);
        totals
            .add_counts(
                physical.substeps,
                physical.chemical_accepted_steps,
                physical.chemical_rejected_steps,
            )
            .map_err(|source| segment_error(segment, start, end, source))?;
        elapsed = elapsed + Dd::from_f64(physical.elapsed_time.value());
    }
    let physical = totals
        .finish(
            initial_volume,
            final_knot.volume.value(),
            final_knot.elapsed_time,
        )
        .map_err(|_| {
            MaterialCompressionHistoryError::Invalid(
                "whole-history original work/energy/entropy totals are unrepresentable",
            )
        })?;
    let elapsed_time_error_s = (elapsed - Dd::from_f64(final_knot.elapsed_time.value())).to_f64();
    if !elapsed_time_error_s.is_finite() {
        return Err(MaterialCompressionHistoryError::Invalid(
            "whole-history elapsed-time sum is unrepresentable",
        ));
    }
    let report = MaterialCompressionHistoryStep {
        physical,
        completed_segments: knots.len() - 1,
        evaluations,
        accepted_adaptive_intervals: accepted,
        rejected_adaptive_intervals: rejected,
        maximum_accepted_error_ratios: maximum_errors,
        suggested_step_s: suggested_step,
        total_chemical_accepted_steps: chemical[0],
        total_chemical_rejected_steps: chemical[1],
        elapsed_time_error_s,
    };
    poll(&mut cancelled)?;
    *parcel = candidate;
    Ok(report)
}
