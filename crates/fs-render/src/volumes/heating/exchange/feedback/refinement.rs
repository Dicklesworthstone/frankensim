//! Complete native thermal trajectories with independently sampled optical
//! response fields. Within a replica, photon streams are shared across evolving
//! temperatures and paired time/history resolutions. Between replicas they are
//! independent. Final-state variation therefore includes radiative feedback,
//! caloric nonlinearities and the covariance of every earlier heat transfer.
//!
//! The gates below are empirical time, sampling and finite-history sensitivity
//! assessments. They are not confidence sequences, certified bias bounds or
//! guarantees under adaptive selection. Rare optical paths may remain unseen.

use super::super::super::{
    CheckedSum,
    ensemble_statistics::{self, StatisticsError},
};
use super::{RadiativeFeedbackError, TemperatureRadiativeExchange};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionSimulation, EnergyAudit, StepConfig,
};
use std::collections::BTreeMap;

/// Physical allowances with relative scales measured from the input state.
/// An arbitrary formation enthalpy never increases the relative allowance.
#[derive(Clone, Copy, Debug)]
pub struct RadiativeEnsembleTolerance {
    pub enthalpy_absolute_j_kg: f64,
    pub temperature_absolute_k: f64,
    pub radiative_energy_absolute_j: f64,
    pub relative: f64,
}

/// Independent gates for temporal discrepancy, final sampling variability and
/// H/2H photon sensitivity, with shared work and retained-storage ceilings.
#[derive(Clone, Copy, Debug)]
pub struct RadiativeFeedbackRefinement {
    pub initial_max_dt_s: f64,
    pub minimum_max_dt_s: f64,
    /// At least two: two successive real native time refinements must pass.
    pub maximum_halvings: usize,
    /// Fixed independent replica count at every resolution, at least two.
    pub replicas: u32,
    /// Replica i uses this key plus i with wrapping u64 arithmetic.
    pub ensemble_seed: u64,
    /// At least one. Volume emission is sampled even without scattering, so
    /// every accepted study includes a paired H/2H comparison.
    pub maximum_photon_doublings: u32,
    /// Multiplies empirical SE; it does not specify a coverage probability.
    pub standard_error_multiplier: f64,
    pub temporal: RadiativeEnsembleTolerance,
    pub sampling: RadiativeEnsembleTolerance,
    pub photon_sensitivity: RadiativeEnsembleTolerance,
    /// The original native joule-conservation gate, separate from accuracy.
    pub conservation_tolerance_j: f64,
    pub maximum_total_substeps: usize,
    pub maximum_total_source_evaluations: usize,
    /// Includes all operator rebuilds in every accepted/discarded trajectory.
    /// Zero remains legal for an identically nonemitting optical field.
    pub maximum_total_traced_segments: u64,
    /// At least 2*replicas, including the native in-flight candidate. No full
    /// time history is retained; original per-state native limits also apply.
    pub maximum_retained_trajectories: usize,
    pub maximum_observable_cells: usize,
    /// Conservative peak storage admission: (2*replicas+4)*max_pairs. Extra
    /// capacity covers the incoming shared snapshot, configured/native clones,
    /// replacement operators and construction tally/output overlap. A slot
    /// admits the pair's optical/mapped representation, not bytes; original
    /// per-operator and mapped-pair limits still apply.
    pub maximum_retained_operator_pair_slots: usize,
    /// Same peak factor times (geometry cells+max_material_cells+max_emitter_cells).
    /// A geometry-cell slot admits its κa/κs packet and associated material/
    /// temperature scratch, including vacuum cells; report slots retain their
    /// metadata/estimates. This is a cell-count cap, not bytes, and does not
    /// shrink under many-to-one thermal maps.
    pub maximum_retained_operator_cell_slots: usize,
}

/// Mean and empirical standard error of that mean across whole trajectories.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RadiativeScalarEstimate {
    pub mean: f64,
    pub standard_error: f64,
}

/// Native enthalpy and actual caloric temperature are averaged separately.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadiativeCellEstimate {
    pub specific_h_j_kg: RadiativeScalarEstimate,
    pub temperature_k: RadiativeScalarEstimate,
}

/// Interval increments over this call, preserving correlations across steps.
/// Internal exchange counts an unordered transfer once, never external heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadiativeEnsembleEnergy<const N: usize> {
    pub net_bath_to_material_j: [RadiativeScalarEstimate; N],
    pub internal_transfer_magnitude_j: [RadiativeScalarEstimate; N],
    pub thermal_external_j: RadiativeScalarEstimate,
    pub unresolved_exchange_scale_j: RadiativeScalarEstimate,
}

/// Paired widths are |mean(fine-coarse)| + k*SE(fine-coarse); final sampling
/// widths are k*SE(fine). Their interpretation stays empirical after selection.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RadiativeEnsembleAssessment {
    pub maximum_enthalpy_width_j_kg: f64,
    pub maximum_temperature_width_k: f64,
    pub maximum_radiative_energy_width_j: f64,
    pub maximum_normalized_width: f64,
    /// Paired h/T gates also refuse requested accuracy below native coordinate
    /// resolution. This floor never widens the caller's error allowance.
    pub maximum_coordinate_resolution_ratio: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RadiativeEnsembleWork {
    pub substeps: usize,
    pub source_evaluations: usize,
    pub traced_segments: u64,
    pub operator_builds: usize,
}

/// Last evaluated resolution and independent assessment gates. A missing gate
/// was not evaluated because an earlier gate refused. No observed order or
/// probabilistic confidence is inferred from these finite sample diagnostics.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RadiativeRefinementDiagnostics {
    pub temporal: Option<RadiativeEnsembleAssessment>,
    pub temporal_refinement_observed: bool,
    pub sampling: Option<RadiativeEnsembleAssessment>,
    pub photon_sensitivity: Option<RadiativeEnsembleAssessment>,
    pub maximum_dt_s: f64,
    pub halvings: usize,
    pub photon_doublings: u32,
    pub histories_per_emitter_lane: u32,
    pub work: RadiativeEnsembleWork,
    pub admitted_peak_operator_pair_slots: usize,
    pub admitted_peak_operator_cell_slots: usize,
}

/// One actual accepted trajectory. No ensemble-mean physical state is made.
pub struct RadiativeFeedbackTrajectory<'t, 'g, const N: usize> {
    pub thermal: DepositionSimulation<'t>,
    pub optical: TemperatureRadiativeExchange<'g, N>,
    pub audit: EnergyAudit,
}

pub struct RadiativeFeedbackEnsemble<'t, 'g, const N: usize> {
    pub trajectories: Vec<RadiativeFeedbackTrajectory<'t, 'g, N>>,
    pub cells: BTreeMap<CellId, RadiativeCellEstimate>,
    pub radiative_energy: RadiativeEnsembleEnergy<N>,
    pub diagnostics: RadiativeRefinementDiagnostics,
}

#[derive(Clone, Debug, PartialEq)]
pub enum RadiativeFeedbackRefinementError {
    Invalid(&'static str),
    Budget(&'static str),
    Physics(RadiativeFeedbackError),
    /// Available refinements did not meet every independent assessment gate.
    Unresolved(RadiativeRefinementDiagnostics),
}
impl From<RadiativeFeedbackError> for RadiativeFeedbackRefinementError {
    fn from(error: RadiativeFeedbackError) -> Self {
        Self::Physics(error)
    }
}
impl From<DepositionError> for RadiativeFeedbackRefinementError {
    fn from(error: DepositionError) -> Self {
        Self::Physics(error.into())
    }
}
impl From<StatisticsError> for RadiativeFeedbackRefinementError {
    fn from(error: StatisticsError) -> Self {
        match error {
            StatisticsError::Invalid(message) => Self::Invalid(message),
            StatisticsError::Heating(error) => {
                Self::Physics(RadiativeFeedbackError::Exchange(error.into()))
            }
        }
    }
}
impl core::fmt::Display for RadiativeFeedbackRefinementError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "radiative feedback refinement refused: {self:?}")
    }
}
impl core::error::Error for RadiativeFeedbackRefinementError {}

type Result<T> = core::result::Result<T, RadiativeFeedbackRefinementError>;
type Checkpoint<'a> = &'a mut dyn FnMut() -> core::result::Result<(), DepositionError>;

impl<'g, const N: usize> TemperatureRadiativeExchange<'g, N> {
    /// Assess complete native thermal trajectories with independently sampled
    /// temperature-dependent reciprocal operators. All trials start from the
    /// same accepted input, conditional on that fixed state; existing input
    /// uncertainty is not resampled. A replica's photon field stays common
    /// across temperatures and paired resolutions, including H/2H prefixes.
    ///
    /// Two adjacent real time refinements must pass, then final sampling and
    /// photon-sensitivity gates pass separately. More photons repeat the time
    /// assessment at the new history count. Every trial shares the same total
    /// work allowance. Native source/material/domain/conservation/cancellation
    /// failures refuse the call transactionally; only unresolved assessment
    /// gates lead to further bounded refinement. Empirical acceptance does not
    /// establish confidence coverage or a bound on nonlinear finite-H bias.
    pub fn advance_ensemble_refined<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        controls: RadiativeFeedbackRefinement,
        mut checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<RadiativeFeedbackEnsemble<'t, 'g, N>> {
        checkpoint()?;
        let (pair_slots, cell_slots) = validate(thermal, self, target_s, controls)?;
        let mut work = RadiativeEnsembleWork::default();
        let mut diagnostics = RadiativeRefinementDiagnostics::default();
        let mut start_halving = 0;
        for doubling in 1..=controls.maximum_photon_doublings {
            checkpoint()?;
            let histories = self
                .histories_per_emitter_lane()
                .checked_mul(
                    1_u32
                        .checked_shl(doubling)
                        .ok_or(budget("photon-count overflow"))?,
                )
                .ok_or(budget("photon-count overflow"))?;
            let mut previous: Option<Vec<RadiativeFeedbackTrajectory<'t, 'g, N>>> = None;
            let mut previous_time: Option<RadiativeEnsembleAssessment> = None;
            let mut dt = controls.initial_max_dt_s * 0.5_f64.powi(start_halving as i32);
            for halving in start_halving..=controls.maximum_halvings {
                if dt < controls.minimum_max_dt_s {
                    break;
                }
                let mut fine = run_ensemble(
                    self,
                    thermal,
                    target_s,
                    dt,
                    histories,
                    controls,
                    &mut work,
                    &mut checkpoint,
                )?;
                diagnostics = RadiativeRefinementDiagnostics {
                    maximum_dt_s: dt,
                    halvings: halving,
                    photon_doublings: doubling,
                    histories_per_emitter_lane: histories,
                    work,
                    admitted_peak_operator_pair_slots: pair_slots,
                    admitted_peak_operator_cell_slots: cell_slots,
                    ..RadiativeRefinementDiagnostics::default()
                };
                if let Some(coarse) = previous.take() {
                    let (_, _, time) = measure_ensemble(
                        thermal,
                        self,
                        &fine,
                        Some(&coarse),
                        controls.temporal,
                        controls,
                        false,
                        &mut checkpoint,
                    )?;
                    let refined = fine
                        .iter()
                        .zip(&coarse)
                        .all(|(a, b)| a.audit.substeps > b.audit.substeps);
                    diagnostics.temporal = Some(time);
                    diagnostics.temporal_refinement_observed = refined && previous_time.is_some();
                    // Release the prior trajectories before allocating the
                    // low-history ensemble; never retain three ensembles.
                    drop(coarse);
                    if let Some(prior) = previous_time.filter(|_| refined) {
                        if prior.maximum_normalized_width <= 1.0
                            && time.maximum_normalized_width <= 1.0
                        {
                            diagnostics.temporal = Some(combine(prior, time));
                            let (cells, energy, sampling) = measure_ensemble(
                                thermal,
                                self,
                                &fine,
                                None,
                                controls.sampling,
                                controls,
                                true,
                                &mut checkpoint,
                            )?;
                            diagnostics.sampling = Some(sampling);
                            if sampling.maximum_normalized_width <= 1.0 {
                                let low = run_ensemble(
                                    self,
                                    thermal,
                                    target_s,
                                    dt,
                                    histories / 2,
                                    controls,
                                    &mut work,
                                    &mut checkpoint,
                                )?;
                                let (_, _, photon) = measure_ensemble(
                                    thermal,
                                    self,
                                    &fine,
                                    Some(&low),
                                    controls.photon_sensitivity,
                                    controls,
                                    false,
                                    &mut checkpoint,
                                )?;
                                diagnostics.photon_sensitivity = Some(photon);
                                diagnostics.work = work;
                                if photon.maximum_normalized_width <= 1.0 {
                                    for candidate in &mut fine {
                                        candidate.optical.budget = self.budget;
                                        // Flat constitutive laws keep this fixed
                                        // owner on later advances. Its effective
                                        // ceiling must leave the trial allowance
                                        // together with the outer observer's.
                                        if let Some(current) = &mut candidate.optical.current {
                                            current.evaluation_limit =
                                                self.budget.max_source_evaluations;
                                        }
                                    }
                                    checkpoint()?;
                                    return Ok(RadiativeFeedbackEnsemble {
                                        trajectories: fine,
                                        cells,
                                        radiative_energy: energy,
                                        diagnostics,
                                    });
                                }
                            }
                            // A new history resolution repeats two time pairs,
                            // starting near the just-assessed physical grid.
                            start_halving = halving.saturating_sub(2);
                            break;
                        }
                    }
                    previous_time = refined.then_some(time);
                }
                previous = Some(fine);
                let next_dt = dt * 0.5;
                if next_dt <= 0.0 || next_dt == dt {
                    break;
                }
                dt = next_dt;
            }
        }
        Err(RadiativeFeedbackRefinementError::Unresolved(diagnostics))
    }
}

fn validate<const N: usize>(
    thermal: &DepositionSimulation<'_>,
    optical: &TemperatureRadiativeExchange<'_, N>,
    target: f64,
    c: RadiativeFeedbackRefinement,
) -> Result<(usize, usize)> {
    optical.check_clock(thermal)?;
    if !target.is_finite()
        || target <= thermal.time_s()
        || !positive(c.initial_max_dt_s)
        || !positive(c.minimum_max_dt_s)
        || c.minimum_max_dt_s > c.initial_max_dt_s * 0.25
        || !(2..=64).contains(&c.maximum_halvings)
        || c.replicas < 2
        || !(1..=31).contains(&c.maximum_photon_doublings)
        || !positive(c.standard_error_multiplier)
        || !nonnegative(c.conservation_tolerance_j)
        || ![c.temporal, c.sampling, c.photon_sensitivity]
            .into_iter()
            .all(|t| {
                positive(t.enthalpy_absolute_j_kg)
                    && positive(t.temperature_absolute_k)
                    && positive(t.radiative_energy_absolute_j)
                    && t.relative.is_finite()
                    && (0.0..1.0).contains(&t.relative)
            })
    {
        return Err(invalid(
            "invalid radiative ensemble, timeline or accuracy controls",
        ));
    }
    let retained = (c.replicas as usize)
        .checked_mul(2)
        .ok_or(budget("trajectory-count overflow"))?;
    let peak_operators = retained
        .checked_add(4)
        .ok_or(budget("operator-count overflow"))?;
    let config = optical.sampling_config();
    let pair_slots = peak_operators
        .checked_mul(config.max_pairs)
        .ok_or(budget("operator pair-slot overflow"))?;
    let cells_per_operator = config
        .max_material_cells
        .checked_add(config.max_emitter_cells)
        .and_then(|cells| cells.checked_add(optical.field().geometry().cell_count()))
        .ok_or(budget("operator cell-slot overflow"))?;
    let cell_slots = peak_operators
        .checked_mul(cells_per_operator)
        .ok_or(budget("operator cell-slot overflow"))?;
    if retained > c.maximum_retained_trajectories
        || pair_slots > c.maximum_retained_operator_pair_slots
        || cell_slots > c.maximum_retained_operator_cell_slots
        || c.maximum_observable_cells == 0
        || thermal.cells().len() > c.maximum_observable_cells
        || c.maximum_total_substeps == 0
        || c.maximum_total_source_evaluations == 0
    {
        return Err(budget("radiative ensemble retained-storage/work allowance"));
    }
    Ok((pair_slots, cell_slots))
}

#[allow(clippy::too_many_arguments)]
fn run_ensemble<'t, 'g, const N: usize>(
    optical: &TemperatureRadiativeExchange<'g, N>,
    thermal: &DepositionSimulation<'t>,
    target: f64,
    dt: f64,
    histories: u32,
    c: RadiativeFeedbackRefinement,
    work: &mut RadiativeEnsembleWork,
    checkpoint: Checkpoint<'_>,
) -> Result<Vec<RadiativeFeedbackTrajectory<'t, 'g, N>>> {
    let mut trajectories = Vec::new();
    trajectories
        .try_reserve_exact(c.replicas as usize)
        .map_err(|_| budget("ensemble allocation"))?;
    for replica in 0..c.replicas {
        checkpoint()?;
        let steps = c
            .maximum_total_substeps
            .checked_sub(work.substeps)
            .ok_or(budget("ensemble substeps"))?;
        let sources = c
            .maximum_total_source_evaluations
            .checked_sub(work.source_evaluations)
            .ok_or(budget("ensemble source evaluations"))?;
        let segments = c
            .maximum_total_traced_segments
            .checked_sub(work.traced_segments)
            .ok_or(budget("ensemble transport segments"))?;
        if steps == 0 || sources == 0 {
            return Err(budget("complete radiative ensemble work exhausted"));
        }
        let trial = optical.configured_trial(
            c.ensemble_seed.wrapping_add(u64::from(replica)),
            histories,
            sources,
            segments,
        )?;
        let (state, source, audit) = thermal.advance_observed_limited(
            &trial,
            target,
            StepConfig {
                max_dt_s: dt,
                energy_tolerance_j: c.conservation_tolerance_j,
            },
            steps,
            &mut *checkpoint,
        )?;
        drop(trial);
        work.substeps = work
            .substeps
            .checked_add(audit.substeps)
            .ok_or(budget("substep-count overflow"))?;
        work.source_evaluations = work
            .source_evaluations
            .checked_add(source.source_evaluations())
            .ok_or(budget("source-count overflow"))?;
        work.traced_segments = work
            .traced_segments
            .checked_add(source.traced_segments())
            .ok_or(budget("segment-count overflow"))?;
        work.operator_builds = work
            .operator_builds
            .checked_add(source.operator_builds())
            .ok_or(budget("operator-count overflow"))?;
        if work.substeps > c.maximum_total_substeps
            || work.source_evaluations > c.maximum_total_source_evaluations
            || work.traced_segments > c.maximum_total_traced_segments
        {
            return Err(budget("complete radiative ensemble work exceeded"));
        }
        if state.cells().len() > c.maximum_observable_cells {
            return Err(budget("ensemble observable cells"));
        }
        if state.time_s() != target || source.time_s() != target {
            return Err(invalid("native radiative ensemble endpoint mismatch"));
        }
        trajectories.push(RadiativeFeedbackTrajectory {
            thermal: state,
            optical: source,
            audit,
        });
    }
    Ok(trajectories)
}

#[allow(clippy::too_many_arguments)]
fn measure_ensemble<const N: usize>(
    initial: &DepositionSimulation<'_>,
    initial_optical: &TemperatureRadiativeExchange<'_, N>,
    fine: &[RadiativeFeedbackTrajectory<'_, '_, N>],
    coarse: Option<&[RadiativeFeedbackTrajectory<'_, '_, N>]>,
    tolerance: RadiativeEnsembleTolerance,
    controls: RadiativeFeedbackRefinement,
    retain_cells: bool,
    checkpoint: Checkpoint<'_>,
) -> Result<(
    BTreeMap<CellId, RadiativeCellEstimate>,
    RadiativeEnsembleEnergy<N>,
    RadiativeEnsembleAssessment,
)> {
    checkpoint()?;
    if fine.len() != controls.replicas as usize || coarse.is_some_and(|c| c.len() != fine.len()) {
        return Err(invalid("radiative replica-count mismatch"));
    }
    let reference = &fine[0].thermal;
    for candidate in fine.iter().chain(coarse.into_iter().flatten()) {
        checkpoint()?;
        if candidate.thermal.time_s() != reference.time_s()
            || candidate
                .thermal
                .cells()
                .keys()
                .ne(reference.cells().keys())
        {
            return Err(invalid("radiative ensemble timeline/topology mismatch"));
        }
    }
    let paired = coarse.is_some();
    let mut assessment = RadiativeEnsembleAssessment::default();
    let mut cells = BTreeMap::new();
    for &id in reference.cells().keys() {
        checkpoint()?;
        let initial_h = initial.cells().get(&id).map(|c| c.specific_h_j_kg);
        let initial_t = if initial_h.is_some() {
            Some(initial.temperature_k(id)?)
        } else {
            None
        };
        let (h, h_width, h_scale, h_magnitude) = measure_scalar(
            fine.len(),
            initial_h,
            paired,
            controls.standard_error_multiplier,
            |r| {
                Ok((
                    fine[r].thermal.cells()[&id].specific_h_j_kg,
                    coarse.map(|c| c[r].thermal.cells()[&id].specific_h_j_kg),
                    None,
                ))
            },
            checkpoint,
        )?;
        let (t, t_width, t_scale, t_magnitude) = measure_scalar(
            fine.len(),
            initial_t,
            paired,
            controls.standard_error_multiplier,
            |r| {
                Ok((
                    fine[r].thermal.temperature_k(id)?,
                    coarse.map(|c| c[r].thermal.temperature_k(id)).transpose()?,
                    None,
                ))
            },
            checkpoint,
        )?;
        assessment.maximum_enthalpy_width_j_kg =
            assessment.maximum_enthalpy_width_j_kg.max(h_width);
        assessment.maximum_temperature_width_k =
            assessment.maximum_temperature_width_k.max(t_width);
        for (width, scale, magnitude, absolute) in [
            (
                h_width,
                h_scale,
                h_magnitude,
                tolerance.enthalpy_absolute_j_kg,
            ),
            (
                t_width,
                t_scale,
                t_magnitude,
                tolerance.temperature_absolute_k,
            ),
        ] {
            assessment.maximum_normalized_width = assessment
                .maximum_normalized_width
                .max(normalize(width, scale, absolute, tolerance.relative)?);
            if paired {
                let floor = (f64::EPSILON * magnitude).max(f64::from_bits(1));
                assessment.maximum_coordinate_resolution_ratio = assessment
                    .maximum_coordinate_resolution_ratio
                    .max(normalize(floor, scale, absolute, tolerance.relative)?);
            }
        }
        assessment.maximum_normalized_width = assessment
            .maximum_normalized_width
            .max(assessment.maximum_coordinate_resolution_ratio);
        if retain_cells {
            cells.insert(
                id,
                RadiativeCellEstimate {
                    specific_h_j_kg: h,
                    temperature_k: t,
                },
            );
        }
    }
    let mut energy = RadiativeEnsembleEnergy {
        net_bath_to_material_j: [RadiativeScalarEstimate::default(); N],
        internal_transfer_magnitude_j: [RadiativeScalarEstimate::default(); N],
        thermal_external_j: RadiativeScalarEstimate::default(),
        unresolved_exchange_scale_j: RadiativeScalarEstimate::default(),
    };
    for channel in 0..4 {
        let lanes = if channel < 2 { N } else { 1 };
        for lane in 0..lanes {
            let base = radiative_sum(initial_optical, channel, lane);
            let (estimate, width, scale, _) = measure_scalar(
                fine.len(),
                Some(0.0),
                paired,
                controls.standard_error_multiplier,
                |r| {
                    Ok((
                        ensemble_statistics::sum_difference(
                            radiative_sum(&fine[r].optical, channel, lane),
                            base,
                        )?,
                        coarse
                            .map(|c| {
                                ensemble_statistics::sum_difference(
                                    radiative_sum(&c[r].optical, channel, lane),
                                    base,
                                )
                            })
                            .transpose()?,
                        coarse
                            .map(|c| {
                                ensemble_statistics::sum_difference(
                                    radiative_sum(&fine[r].optical, channel, lane),
                                    radiative_sum(&c[r].optical, channel, lane),
                                )
                            })
                            .transpose()?,
                    ))
                },
                checkpoint,
            )?;
            match channel {
                0 => energy.net_bath_to_material_j[lane] = estimate,
                1 => energy.internal_transfer_magnitude_j[lane] = estimate,
                2 => energy.thermal_external_j = estimate,
                _ => energy.unresolved_exchange_scale_j = estimate,
            }
            assessment.maximum_radiative_energy_width_j =
                assessment.maximum_radiative_energy_width_j.max(width);
            assessment.maximum_normalized_width =
                assessment.maximum_normalized_width.max(normalize(
                    width,
                    scale,
                    tolerance.radiative_energy_absolute_j,
                    tolerance.relative,
                )?);
        }
    }
    Ok((cells, energy, assessment))
}

fn radiative_sum<const N: usize>(
    optical: &TemperatureRadiativeExchange<'_, N>,
    channel: usize,
    lane: usize,
) -> CheckedSum {
    let Some(current) = &optical.current else {
        return CheckedSum::default();
    };
    match channel {
        0 => current.sums.bath[lane],
        1 => current.sums.internal[lane],
        2 => current.sums.external,
        _ => current.sums.unresolved,
    }
}

fn measure_scalar(
    count: usize,
    base: Option<f64>,
    paired: bool,
    multiplier: f64,
    values: impl FnMut(usize) -> Result<(f64, Option<f64>, Option<f64>)>,
    checkpoint: Checkpoint<'_>,
) -> Result<(RadiativeScalarEstimate, f64, f64, f64)> {
    let (estimate, width, scale, magnitude) =
        ensemble_statistics::measure_scalar(count, base, paired, multiplier, values, checkpoint)?;
    Ok((
        RadiativeScalarEstimate {
            mean: estimate.mean,
            standard_error: estimate.standard_error,
        },
        width,
        scale,
        magnitude,
    ))
}
fn combine(
    a: RadiativeEnsembleAssessment,
    b: RadiativeEnsembleAssessment,
) -> RadiativeEnsembleAssessment {
    RadiativeEnsembleAssessment {
        maximum_enthalpy_width_j_kg: a
            .maximum_enthalpy_width_j_kg
            .max(b.maximum_enthalpy_width_j_kg),
        maximum_temperature_width_k: a
            .maximum_temperature_width_k
            .max(b.maximum_temperature_width_k),
        maximum_radiative_energy_width_j: a
            .maximum_radiative_energy_width_j
            .max(b.maximum_radiative_energy_width_j),
        maximum_normalized_width: a.maximum_normalized_width.max(b.maximum_normalized_width),
        maximum_coordinate_resolution_ratio: a
            .maximum_coordinate_resolution_ratio
            .max(b.maximum_coordinate_resolution_ratio),
    }
}
fn normalize(width: f64, scale: f64, absolute: f64, relative: f64) -> Result<f64> {
    let allowance = absolute + relative * scale;
    if !allowance.is_finite() {
        return Err(invalid("unrepresentable ensemble allowance"));
    }
    // An overflowing width/allowance ratio fails the gate instead of admitting.
    Ok(width / allowance)
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn invalid(message: &'static str) -> RadiativeFeedbackRefinementError {
    RadiativeFeedbackRefinementError::Invalid(message)
}
fn budget(message: &'static str) -> RadiativeFeedbackRefinementError {
    RadiativeFeedbackRefinementError::Budget(message)
}
