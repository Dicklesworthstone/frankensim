//! Adaptive time intervals for the original reacting compression owner.
//!
//! Every trial compares one full Strang interval with two half intervals along
//! the SAME prescribed logarithmic volume trajectory. Temperature, original
//! species amounts and interval pressure work must all satisfy the declared
//! local comparison tolerances. Only the two-half native states are retained;
//! no extrapolation, species clipping or replacement thermodynamics is used.
//! This is local temporal error control, not a global accuracy certificate.

use super::{
    ClosedGasParcel, CompressionInterval, CompressionTotals, GasEntropySolve,
    MAX_COMPRESSION_SUBSTEPS, MaterialCompressionError, MaterialCompressionSolve,
    MaterialCompressionStep, MaterialKineticsSolve, ThermochemicalGasModel,
    advance_closed_parcel_compression, advance_interval, positive,
};
use crate::chemistry::IsochoricKineticsMethod;
use fs_qty::{Time, Volume};
use fs_thermochem::{SpeciesId, kinetics::IdealGasKineticsV1};

/// Three original interval evaluations per full/two-half comparison.
pub const MAX_COMPRESSION_EVALUATIONS: usize = 3 * MAX_COMPRESSION_SUBSTEPS;

/// Original physical controls and explicit local time-comparison tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveMaterialCompressionSolve {
    pub kinetics: MaterialKineticsSolve,
    pub volume: GasEntropySolve,
    /// Positive absolute full/fine temperature difference allowance [K].
    pub absolute_temperature_tolerance_k: f64,
    /// Positive absolute allowance for EACH original species inventory [mol].
    /// Initially absent species are compared against zero, never trace seeded.
    pub absolute_amount_tolerance_mol: f64,
    /// Positive absolute allowance for interval pressure work [J].
    pub absolute_work_tolerance_j: f64,
    /// Nonnegative and below one. For temperature/species, multiplies the
    /// largest coarse/fine CHANGE from the accepted starting state. For work,
    /// multiplies the largest coarse/fine interval work magnitude. Large carrier
    /// inventories or accumulated work cannot dilute the local comparison.
    pub relative_tolerance: f64,
    pub initial_step_s: f64,
    /// Smallest retry interval. The exact final remainder may be smaller, but
    /// cannot be further reduced after a failed comparison.
    pub minimum_step_s: f64,
    pub maximum_step_s: f64,
    /// Accepted full/two-half comparisons, 1..=4096; each retains two intervals.
    pub maximum_intervals: usize,
    /// Counts ALL coarse/fine original interval calls, including rejected
    /// comparisons: 3..=12288. Native chemistry has its own per-call work bound.
    pub maximum_evaluations: usize,
}

/// Physical totals come only from the accepted fine trajectory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveMaterialCompressionStep {
    /// Same work/energy/residual/stage-entropy semantics as fixed-step evolution.
    /// External heat is zero: this remains an adiabatic stroke.
    pub physical: MaterialCompressionStep,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub evaluations: usize,
    /// Temperature, maximum original-species inventory, interval work. These
    /// are raw full/fine differences divided by declared comparison allowances,
    /// without an assumed-order Richardson divisor. All accepted ratios <= 1.
    /// A fixed-volume call bypasses splitting and reports zero ratios.
    pub maximum_accepted_error_ratios: [f64; 3],
    pub minimum_accepted_step_s: f64,
    pub maximum_accepted_step_s: f64,
    pub suggested_step_s: f64,
    /// Native chemical work across ALL successful probes, including coarse
    /// and rejected trajectories. Physical counts include accepted fine only.
    pub total_chemical_accepted_steps: usize,
    pub total_chemical_rejected_steps: usize,
}

/// Any refusal preserves the complete original parcel, even after an accepted prefix.
#[derive(Clone, Debug, PartialEq)]
pub enum AdaptiveMaterialCompressionError {
    Invalid(&'static str),
    /// Native source, kinetics or volume inversion refusal. These failures
    /// are terminal; the controller only retries an actual error comparison.
    Physical(MaterialCompressionError),
    Accuracy {
        elapsed_time_s: f64,
        step_s: f64,
        error_ratios: [f64; 3],
    },
    WorkLimit {
        elapsed_time_s: f64,
        accepted_intervals: usize,
        rejected_intervals: usize,
        evaluations: usize,
    },
    TimeResolution {
        elapsed_time_s: f64,
        target_time_s: f64,
    },
    Cancelled,
}
impl core::fmt::Display for AdaptiveMaterialCompressionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive reacting compression: {self:?}")
    }
}
impl std::error::Error for AdaptiveMaterialCompressionError {}
impl From<MaterialCompressionError> for AdaptiveMaterialCompressionError {
    fn from(error: MaterialCompressionError) -> Self {
        match error {
            MaterialCompressionError::Cancelled => Self::Cancelled,
            other => Self::Physical(other),
        }
    }
}

impl AdaptiveMaterialCompressionSolve {
    fn validate(self) -> Result<(), AdaptiveMaterialCompressionError> {
        if ![
            self.absolute_temperature_tolerance_k,
            self.absolute_amount_tolerance_mol,
            self.absolute_work_tolerance_j,
            self.initial_step_s,
            self.minimum_step_s,
            self.maximum_step_s,
        ]
        .into_iter()
        .all(positive)
            || !self.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_tolerance)
            || self.minimum_step_s > self.initial_step_s
            || self.initial_step_s > self.maximum_step_s
            || !(1..=MAX_COMPRESSION_SUBSTEPS).contains(&self.maximum_intervals)
            || !(3..=MAX_COMPRESSION_EVALUATIONS).contains(&self.maximum_evaluations)
        {
            return Err(AdaptiveMaterialCompressionError::Invalid(
                "positive absolute tolerances, 0<=relative<1, ordered positive step bounds and bounded work required",
            ));
        }
        Ok(())
    }
}

/// Complete a prescribed reacting stroke with bounded full/two-half comparisons.
///
/// The original native chemistry integrator retains its own error, conservation,
/// source-domain and work controls. The outer comparison includes its numerical
/// differences; tolerances below that native accuracy may refuse at the declared
/// minimum step. There is no independent total entropy certificate (see the
/// physical report), global time-error bound, ignition or empirical rate claim.
/// Fixed volume, including zero duration, delegates to the original one-reactor
/// path after control admission. Moving strokes require positive duration.
/// The exact requested final volume and elapsed time are retained.
pub fn advance_closed_parcel_compression_adaptive(
    parcel: &mut ClosedGasParcel,
    bank: &ThermochemicalGasModel,
    mechanism: &IdealGasKineticsV1,
    final_volume: Volume,
    duration: Time,
    controls: AdaptiveMaterialCompressionSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<AdaptiveMaterialCompressionStep, AdaptiveMaterialCompressionError> {
    poll(&mut cancelled)?;
    controls.validate()?;
    let fixed = MaterialCompressionSolve {
        kinetics: controls.kinetics,
        volume: controls.volume,
        substeps: 1,
    };
    let initial_volume = parcel.volume_m3();
    let final_volume_m3 = final_volume.value();
    let duration_s = duration.value();
    if !positive(final_volume_m3) || !duration_s.is_finite() || duration_s < 0.0 {
        return Err(AdaptiveMaterialCompressionError::Invalid(
            "positive finite final volume and nonnegative finite duration required",
        ));
    }
    if final_volume_m3 == initial_volume {
        let physical = advance_closed_parcel_compression(
            parcel,
            bank,
            mechanism,
            final_volume,
            duration,
            fixed,
            &mut cancelled,
        )?;
        return Ok(AdaptiveMaterialCompressionStep {
            physical,
            accepted_intervals: usize::from(duration_s > 0.0),
            rejected_intervals: 0,
            evaluations: 1,
            maximum_accepted_error_ratios: [0.0; 3],
            minimum_accepted_step_s: duration_s,
            maximum_accepted_step_s: duration_s,
            suggested_step_s: controls.initial_step_s,
            total_chemical_accepted_steps: physical.chemical_accepted_steps,
            total_chemical_rejected_steps: physical.chemical_rejected_steps,
        });
    }
    let ratio = final_volume_m3 / initial_volume;
    if !positive(duration_s) || !positive(ratio) {
        return Err(AdaptiveMaterialCompressionError::Invalid(
            "moving strokes require positive elapsed time and a representable volume ratio",
        ));
    }
    let log_ratio = fs_math::det::ln(ratio);
    let volume_at = |time_s: f64| {
        if time_s == duration_s {
            final_volume_m3
        } else {
            initial_volume * fs_math::det::exp(log_ratio * (time_s / duration_s))
        }
    };
    // Canonical original ordering is bound once. Error comparisons merge the
    // active sorted species with this ordering; no per-retry keyed maps exist.
    let species = bank.thermodynamics().components();
    let mut candidate = parcel.clone();
    candidate
        .set_volume_isentropic(initial_volume, controls.volume, &mut cancelled)
        .map_err(MaterialCompressionError::from)?;
    let mut totals = CompressionTotals::default();
    let (mut elapsed, mut step) = (0.0, controls.initial_step_s);
    let (mut accepted, mut rejected, mut evaluations) = (0, 0, 0);
    let mut chemical_work = [0_usize; 2];
    let mut maxima = [0.0_f64; 3];
    let (mut minimum, mut maximum) = (f64::INFINITY, 0.0_f64);

    while elapsed < duration_s {
        poll(&mut cancelled)?;
        if accepted == controls.maximum_intervals || evaluations + 3 > controls.maximum_evaluations
        {
            return Err(AdaptiveMaterialCompressionError::WorkLimit {
                elapsed_time_s: elapsed,
                accepted_intervals: accepted,
                rejected_intervals: rejected,
                evaluations,
            });
        }
        let end = if step >= duration_s - elapsed {
            duration_s
        } else {
            elapsed + step
        };
        let span = end - elapsed;
        let middle = elapsed + 0.5 * span;
        let first_midpoint = elapsed + 0.5 * (middle - elapsed);
        let second_midpoint = middle + 0.5 * (end - middle);
        if !(elapsed < first_midpoint
            && first_midpoint < middle
            && middle < second_midpoint
            && second_midpoint < end)
        {
            return Err(AdaptiveMaterialCompressionError::TimeResolution {
                elapsed_time_s: elapsed,
                target_time_s: end,
            });
        }
        let endpoint = volume_at(end);
        let midpoint = volume_at(middle);
        let mut coarse = candidate.clone();
        let mut coarse_totals = CompressionTotals::default();
        evaluations += 1;
        advance_interval(
            &mut coarse,
            bank,
            mechanism,
            CompressionInterval {
                duration_s: span,
                midpoint_volume_m3: midpoint,
                final_volume_m3: endpoint,
            },
            fixed,
            &mut coarse_totals,
            &mut cancelled,
        )?;
        let mut fine = candidate.clone();
        let mut fine_totals = CompressionTotals::default();
        evaluations += 1;
        advance_interval(
            &mut fine,
            bank,
            mechanism,
            CompressionInterval {
                duration_s: middle - elapsed,
                midpoint_volume_m3: volume_at(first_midpoint),
                final_volume_m3: midpoint,
            },
            fixed,
            &mut fine_totals,
            &mut cancelled,
        )?;
        evaluations += 1;
        advance_interval(
            &mut fine,
            bank,
            mechanism,
            CompressionInterval {
                duration_s: end - middle,
                midpoint_volume_m3: volume_at(second_midpoint),
                final_volume_m3: endpoint,
            },
            fixed,
            &mut fine_totals,
            &mut cancelled,
        )?;
        // Validate native totals before deriving any error or publishing state.
        let coarse_report =
            coarse_totals.finish(candidate.volume_m3(), endpoint, Time::new(span))?;
        let fine_report = fine_totals.finish(candidate.volume_m3(), endpoint, Time::new(span))?;
        for (i, count) in [
            coarse_report
                .chemical_accepted_steps
                .checked_add(fine_report.chemical_accepted_steps),
            coarse_report
                .chemical_rejected_steps
                .checked_add(fine_report.chemical_rejected_steps),
        ]
        .into_iter()
        .enumerate()
        {
            chemical_work[i] = count.and_then(|n| chemical_work[i].checked_add(n)).ok_or(
                AdaptiveMaterialCompressionError::Invalid("chemical work count overflow"),
            )?;
        }
        let errors = compare(
            species.iter().map(|source| source.species()),
            [&candidate, &coarse, &fine],
            [
                coarse_report.boundary_work_on_gas_j,
                fine_report.boundary_work_on_gas_j,
            ],
            controls,
            &mut cancelled,
        )?;
        let error = errors.into_iter().fold(0.0_f64, f64::max);
        if error <= 1.0 {
            totals.merge(fine_totals)?;
            candidate = fine;
            elapsed = end;
            accepted += 1;
            minimum = minimum.min(span);
            maximum = maximum.max(span);
            for i in 0..3 {
                maxima[i] = maxima[i].max(errors[i]);
            }
            let factor = if error == 0.0 {
                2.0
            } else {
                (0.9 / error_root(error, controls.kinetics.chemistry.method)).clamp(0.5, 2.0)
            };
            step = if span >= controls.maximum_step_s / factor {
                controls.maximum_step_s
            } else {
                (span * factor).max(controls.minimum_step_s)
            };
        } else {
            rejected += 1;
            let factor =
                (0.8 / error_root(error, controls.kinetics.chemistry.method)).clamp(0.1, 0.5);
            let next = (span * factor).max(controls.minimum_step_s);
            if span <= controls.minimum_step_s || !positive(next) || next >= span {
                return Err(AdaptiveMaterialCompressionError::Accuracy {
                    elapsed_time_s: elapsed,
                    step_s: span,
                    error_ratios: errors,
                });
            }
            step = next;
        }
    }
    let physical = totals.finish(initial_volume, final_volume_m3, duration)?;
    let report = AdaptiveMaterialCompressionStep {
        physical,
        accepted_intervals: accepted,
        rejected_intervals: rejected,
        evaluations,
        maximum_accepted_error_ratios: maxima,
        minimum_accepted_step_s: minimum,
        maximum_accepted_step_s: maximum,
        suggested_step_s: step,
        total_chemical_accepted_steps: chemical_work[0],
        total_chemical_rejected_steps: chemical_work[1],
    };
    poll(&mut cancelled)?;
    *parcel = candidate;
    Ok(report)
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AdaptiveMaterialCompressionError> {
    if cancelled() {
        Err(AdaptiveMaterialCompressionError::Cancelled)
    } else {
        Ok(())
    }
}

/// Only a step-size heuristic: acceptance always tests the raw comparison.
pub(super) fn error_root(error: f64, method: IsochoricKineticsMethod) -> f64 {
    match method {
        IsochoricKineticsMethod::BackwardEuler => fs_math::det::sqrt(error),
        IsochoricKineticsMethod::Sdirk2 | IsochoricKineticsMethod::ExactFirstOrder { .. } => {
            fs_math::det::pow(error, 1.0 / 3.0)
        }
    }
}

pub(super) fn scaled_difference(
    initial: f64,
    coarse: f64,
    fine: f64,
    absolute: f64,
    relative: f64,
) -> Result<f64, AdaptiveMaterialCompressionError> {
    let change = (coarse - initial).abs().max((fine - initial).abs());
    let scale = absolute + relative * change;
    let difference = (fine - coarse).abs();
    if !positive(scale) || !difference.is_finite() {
        return Err(AdaptiveMaterialCompressionError::Invalid(
            "local comparison scale or difference is unrepresentable",
        ));
    }
    // An overflowing ratio is an accuracy rejection, never a passing estimate.
    Ok(difference / scale)
}

pub(super) fn compare<'s>(
    sources: impl IntoIterator<Item = &'s SpeciesId>,
    states: [&ClosedGasParcel; 3],
    work: [f64; 2],
    controls: AdaptiveMaterialCompressionSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; 3], AdaptiveMaterialCompressionError> {
    let temperatures = states.map(|p| p.state().state().temperature);
    let mut errors = [
        scaled_difference(
            temperatures[0],
            temperatures[1],
            temperatures[2],
            controls.absolute_temperature_tolerance_k,
            controls.relative_tolerance,
        )?,
        0.0_f64,
        scaled_difference(
            0.0,
            work[0],
            work[1],
            controls.absolute_work_tolerance_j,
            controls.relative_tolerance,
        )?,
    ];
    let total_moles =
        states.map(|p| p.mass_kg() / p.state().thermodynamics().properties().molar_mass().value());
    if !total_moles.into_iter().all(positive) {
        return Err(AdaptiveMaterialCompressionError::Invalid(
            "unrepresentable parcel inventory",
        ));
    }
    let mut indices = [0; 3];
    for source in sources {
        poll(cancelled)?;
        let mut amounts = [0.0; 3];
        for i in 0..3 {
            let mixture = states[i].model().thermodynamics();
            if mixture
                .components()
                .get(indices[i])
                .is_some_and(|model| model.species() == source)
            {
                amounts[i] = total_moles[i] * mixture.mole_composition().fractions()[indices[i]];
                if !positive(amounts[i]) {
                    return Err(AdaptiveMaterialCompressionError::Invalid(
                        "positive species inventory is unrepresentable",
                    ));
                }
                indices[i] += 1;
            }
        }
        errors[1] = errors[1].max(scaled_difference(
            amounts[0],
            amounts[1],
            amounts[2],
            controls.absolute_amount_tolerance_mol,
            controls.relative_tolerance,
        )?);
    }
    if (0..3).any(|i| indices[i] != states[i].model().thermodynamics().components().len()) {
        return Err(AdaptiveMaterialCompressionError::Invalid(
            "unmapped original species inventory",
        ));
    }
    Ok(errors)
}
