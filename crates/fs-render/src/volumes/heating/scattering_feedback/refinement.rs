//! Bounded complete-trajectory ensembles for stochastic thermal feedback.
//!
//! Sampling uncertainty belongs to whole trajectories: earlier random heat
//! changes later coefficients, caloric slopes and transport. Per-step packet
//! errors therefore cannot be propagated as independent thermal errors. Replica
//! streams are independent; matching replicas share streams across timestep and
//! photon-count comparisons. These are empirical sensitivity/standard-error
//! diagnostics, not confidence sequences, certified bias bounds or guarantees
//! under optional stopping. Finite samples can still miss rare transport paths.

use super::super::CheckedSum;
use super::super::feedback::OpticalFeedbackError;
use super::TemperatureScatteringHeating;
use fs_math::det;
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionSimulation, EnergyAudit, StepConfig,
};
use std::collections::BTreeMap;

/// Absolute physical tolerances plus a scale based on change from the original
/// state. Arbitrary enthalpy reference energies never enlarge the allowance.
#[derive(Clone, Copy, Debug)]
pub struct ScatteringEnsembleTolerance {
    pub enthalpy_absolute_j_kg: f64,
    pub temperature_absolute_k: f64,
    pub optical_energy_absolute_j: f64,
    pub relative: f64,
}

/// Independent budgets for time error, final sampling noise and photon-count
/// sensitivity. Each is tested for every thermal cell and optical energy lane.
#[derive(Clone, Copy, Debug)]
pub struct ScatteringFeedbackRefinement {
    pub initial_max_dt_s: f64,
    pub minimum_max_dt_s: f64,
    /// At least two; two successive coarse/fine comparisons must pass.
    pub maximum_halvings: usize,
    /// Fixed replica count at each resolution, at least two.
    pub replicas: u32,
    /// Replica i uses this seed plus i, with wrapping u64 arithmetic. Distinct
    /// replicas have distinct stream identities; refinements retain identity.
    pub ensemble_seed: u64,
    /// A scattering calculation starts with H versus 2H. At least one doubling
    /// is required unless every compiled scattering coefficient is zero.
    pub maximum_photon_doublings: u32,
    /// Multiplies the empirical SE, without assigning a coverage probability.
    pub standard_error_multiplier: f64,
    pub temporal: ScatteringEnsembleTolerance,
    pub sampling: ScatteringEnsembleTolerance,
    pub photon_sensitivity: ScatteringEnsembleTolerance,
    pub conservation_tolerance_j: f64,
    /// Includes all discarded replicas and all time/photon resolutions.
    pub maximum_total_substeps: usize,
    pub maximum_total_source_evaluations: usize,
    pub maximum_total_traced_segments: u64,
    /// At least twice `replicas`. Counts additional native thermal snapshots,
    /// including an in-flight native candidate. Each retains the original
    /// thermal/grid allocation ceilings; no complete timestep history is kept.
    pub maximum_retained_trajectories: usize,
    /// Bounds per-cell summary storage, including any newly activated cells.
    pub maximum_observable_cells: usize,
}

/// Mean and empirical standard error of that mean across complete replicas.
/// This is not a confidence interval or a bound on model/discretization error.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ScatteringScalarEstimate {
    pub mean: f64,
    pub standard_error: f64,
}

/// Both estimates use actual native coordinates. Temperature is averaged after
/// evaluating each replica's caloric law, never by inverting mean enthalpy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScatteringCellEstimate {
    pub specific_h_j_kg: ScatteringScalarEstimate,
    pub temperature_k: ScatteringScalarEstimate,
}

/// Optical energy INCREMENTS over this invocation, not historic cumulative
/// energy. Per-lane moments retain the covariance of all steps of a trajectory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScatteringEnsembleEnergy<const N: usize> {
    pub incident_j: [ScatteringScalarEstimate; N],
    pub ballistic_escape_j: [ScatteringScalarEstimate; N],
    pub diffuse_escape_j: [ScatteringScalarEstimate; N],
    pub absorbed_j: [ScatteringScalarEstimate; N],
    pub unresolved_j: [ScatteringScalarEstimate; N],
    pub residual_j: [ScatteringScalarEstimate; N],
    pub thermal_input_j: ScatteringScalarEstimate,
}

/// For paired comparisons the tested width is |mean(fine-coarse)| + k SE of
/// paired differences. For final sampling it is k SE of the finest observable.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ScatteringEnsembleAssessment {
    pub maximum_enthalpy_width_j_kg: f64,
    pub maximum_temperature_width_k: f64,
    pub maximum_optical_energy_width_j: f64,
    /// Largest empirical width/allowance or paired h/T representation-floor
    /// ratio. A sub-resolution user allowance is refused, never enlarged.
    pub maximum_normalized_width: f64,
    pub maximum_coordinate_resolution_ratio: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ScatteringEnsembleWork {
    pub substeps: usize,
    pub source_evaluations: usize,
    pub traced_segments: u64,
}

/// Empirical diagnostics at the last evaluated resolution. Absence means a
/// preceding gate failed before that assessment; photon sensitivity is also
/// absent for the exact no-scattering path. No convergence order is inferred.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ScatteringRefinementDiagnostics {
    pub temporal: Option<ScatteringEnsembleAssessment>,
    /// Both final comparisons actually increased native substep counts in
    /// every replica. Merely lowering an inactive CFL/event ceiling is not
    /// observed time refinement, even if resulting trajectories are identical.
    pub temporal_refinement_observed: bool,
    pub sampling: Option<ScatteringEnsembleAssessment>,
    pub photon_sensitivity: Option<ScatteringEnsembleAssessment>,
    pub maximum_dt_s: f64,
    pub halvings: usize,
    pub photon_doublings: u32,
    pub histories_per_lane: u32,
    pub work: ScatteringEnsembleWork,
}

/// One genuine accepted native trajectory, with its original caloric state,
/// radiation accounts and conservation audit. No mean physical state is made.
pub struct ScatteringFeedbackTrajectory<'t, 'g, const N: usize> {
    pub thermal: DepositionSimulation<'t>,
    pub optical: TemperatureScatteringHeating<'g, N>,
    pub audit: EnergyAudit,
}

pub struct ScatteringFeedbackEnsemble<'t, 'g, const N: usize> {
    pub trajectories: Vec<ScatteringFeedbackTrajectory<'t, 'g, N>>,
    pub cells: BTreeMap<CellId, ScatteringCellEstimate>,
    pub optical_energy: ScatteringEnsembleEnergy<N>,
    pub diagnostics: ScatteringRefinementDiagnostics,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ScatteringRefinementError {
    Invalid(&'static str),
    Budget(&'static str),
    Physics(OpticalFeedbackError),
    /// The bounded calculation could not meet every independent accuracy gate.
    Unresolved(ScatteringRefinementDiagnostics),
}
impl From<OpticalFeedbackError> for ScatteringRefinementError {
    fn from(error: OpticalFeedbackError) -> Self {
        Self::Physics(error)
    }
}
impl From<DepositionError> for ScatteringRefinementError {
    fn from(error: DepositionError) -> Self {
        Self::Physics(error.into())
    }
}
impl std::fmt::Display for ScatteringRefinementError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "scattering feedback refinement refused: {self:?}")
    }
}
impl std::error::Error for ScatteringRefinementError {}

type Result<T> = std::result::Result<T, ScatteringRefinementError>;
type Checkpoint<'a> = &'a mut dyn FnMut() -> std::result::Result<(), DepositionError>;

impl<'g, const N: usize> TemperatureScatteringHeating<'g, N> {
    /// Advance bounded independent whole-trajectory replicas through the native
    /// thermal/transport owners. Every trial starts from the supplied accepted
    /// snapshot, conditional on that fixed starting state; uncertainty already
    /// present in the input state is not resampled. Corresponding replicas share
    /// only their deterministic stream
    /// identities across resolutions; thermal steps use fresh physical-time
    /// streams. Refusal/cancellation never mutates either input.
    ///
    /// Two adjacent time comparisons must pass at the accepted photon count.
    /// Final ensemble sampling noise is tested independently of paired time
    /// differences. For scattering, a paired H/2H comparison then tests finite
    /// photon sensitivity at that same finest time resolution. Failure can
    /// trigger another bounded photon doubling, which rechecks time accuracy.
    /// Selection uses empirical fixed-replica diagnostics and does NOT confer
    /// anytime-valid confidence or certify nonlinear finite-photon bias.
    pub fn advance_ensemble_refined<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        controls: ScatteringFeedbackRefinement,
        mut checkpoint: impl FnMut() -> std::result::Result<(), DepositionError>,
    ) -> Result<ScatteringFeedbackEnsemble<'t, 'g, N>> {
        checkpoint()?;
        validate(thermal, self, target_s, controls)?;
        let scattering = self.has_scattering();
        let first_doubling = u32::from(scattering);
        let last_doubling = if scattering {
            controls.maximum_photon_doublings
        } else {
            0
        };
        let mut work = ScatteringEnsembleWork::default();
        let mut diagnostics = ScatteringRefinementDiagnostics::default();
        let mut start_halving = 0;
        for doubling in first_doubling..=last_doubling {
            checkpoint()?;
            let histories = self
                .histories_per_lane()
                .checked_mul(
                    1_u32
                        .checked_shl(doubling)
                        .ok_or(ScatteringRefinementError::Budget("photon-count overflow"))?,
                )
                .ok_or(ScatteringRefinementError::Budget("photon-count overflow"))?;
            let mut previous: Option<Vec<ScatteringFeedbackTrajectory<'t, 'g, N>>> = None;
            let mut previous_time: Option<ScatteringEnsembleAssessment> = None;
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
                diagnostics = ScatteringRefinementDiagnostics {
                    maximum_dt_s: dt,
                    halvings: halving,
                    photon_doublings: doubling,
                    histories_per_lane: histories,
                    work,
                    ..ScatteringRefinementDiagnostics::default()
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
                    // Drop the coarse physical ensemble before allocating the
                    // low-photon comparison; never retain three ensembles.
                    drop(coarse);
                    if let Some(prior) = previous_time.filter(|_| refined) {
                        if prior.maximum_normalized_width <= 1.0
                            && time.maximum_normalized_width <= 1.0
                        {
                            diagnostics.temporal = Some(combine(prior, time));
                            let (cells, optical_energy, sampling) = measure_ensemble(
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
                                if scattering {
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
                                    if photon.maximum_normalized_width > 1.0 {
                                        start_halving = halving.saturating_sub(2);
                                        break;
                                    }
                                }
                                // Shared work has already been charged. Keep the
                                // caller's original per-advance ceilings for a
                                // future continuation of each genuine endpoint.
                                for candidate in &mut fine {
                                    candidate.optical.budget = self.budget;
                                }
                                checkpoint()?;
                                return Ok(ScatteringFeedbackEnsemble {
                                    trajectories: fine,
                                    cells,
                                    optical_energy,
                                    diagnostics,
                                });
                            }
                            // More photons address propagated sampling noise;
                            // time comparisons are repeated at the new count.
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
        Err(ScatteringRefinementError::Unresolved(diagnostics))
    }
}

fn validate<const N: usize>(
    thermal: &DepositionSimulation<'_>,
    optical: &TemperatureScatteringHeating<'_, N>,
    target: f64,
    c: ScatteringFeedbackRefinement,
) -> Result<()> {
    if !target.is_finite()
        || target <= thermal.time_s()
        || thermal.time_s() != optical.time_s()
        || !positive(c.initial_max_dt_s)
        || !positive(c.minimum_max_dt_s)
        || c.minimum_max_dt_s > c.initial_max_dt_s * 0.25
        || !(2..=64).contains(&c.maximum_halvings)
        || c.replicas < 2
        || c.maximum_photon_doublings > 31
        || (optical.has_scattering() && c.maximum_photon_doublings == 0)
        || !positive(c.standard_error_multiplier)
        || !nonnegative(c.conservation_tolerance_j)
        || ![c.temporal, c.sampling, c.photon_sensitivity]
            .into_iter()
            .all(|t| {
                positive(t.enthalpy_absolute_j_kg)
                    && positive(t.temperature_absolute_k)
                    && positive(t.optical_energy_absolute_j)
                    && t.relative.is_finite()
                    && (0.0..1.0).contains(&t.relative)
            })
    {
        return Err(ScatteringRefinementError::Invalid(
            "invalid ensemble, timeline or accuracy controls",
        ));
    }
    let retained =
        (c.replicas as usize)
            .checked_mul(2)
            .ok_or(ScatteringRefinementError::Budget(
                "trajectory-count overflow",
            ))?;
    if retained > c.maximum_retained_trajectories
        || c.maximum_observable_cells == 0
        || thermal.cells().len() > c.maximum_observable_cells
        || c.maximum_total_substeps == 0
        || c.maximum_total_source_evaluations == 0
        || c.maximum_total_traced_segments == 0
    {
        return Err(ScatteringRefinementError::Budget(
            "ensemble storage/work allowance",
        ));
    }
    Ok(())
}

fn run_ensemble<'t, 'g, const N: usize>(
    optical: &TemperatureScatteringHeating<'g, N>,
    thermal: &DepositionSimulation<'t>,
    target: f64,
    dt: f64,
    histories: u32,
    c: ScatteringFeedbackRefinement,
    work: &mut ScatteringEnsembleWork,
    checkpoint: Checkpoint<'_>,
) -> Result<Vec<ScatteringFeedbackTrajectory<'t, 'g, N>>> {
    let mut trajectories = Vec::new();
    trajectories
        .try_reserve_exact(c.replicas as usize)
        .map_err(|_| ScatteringRefinementError::Budget("ensemble allocation"))?;
    for replica in 0..c.replicas {
        checkpoint()?;
        let remaining_steps = c.maximum_total_substeps - work.substeps;
        let remaining_sources = c.maximum_total_source_evaluations - work.source_evaluations;
        let remaining_segments = c.maximum_total_traced_segments - work.traced_segments;
        if remaining_steps == 0 || remaining_sources == 0 || remaining_segments == 0 {
            return Err(ScatteringRefinementError::Budget(
                "complete-ensemble work exhausted",
            ));
        }
        // Configure one clone: chaining three clone-based public setters would
        // copy the mutable cell-coefficient arrays three times per replica.
        let mut source = optical.with_histories_per_lane(histories)?;
        source.replica_seed = c.ensemble_seed.wrapping_add(u64::from(replica));
        source.budget.max_source_evaluations =
            source.budget.max_source_evaluations.min(remaining_sources);
        source.budget.max_traced_segments =
            source.budget.max_traced_segments.min(remaining_segments);
        let (state, source, audit) = thermal.advance_observed_limited(
            &source,
            target,
            StepConfig {
                max_dt_s: dt,
                energy_tolerance_j: c.conservation_tolerance_j,
            },
            remaining_steps,
            &mut *checkpoint,
        )?;
        work.substeps = work
            .substeps
            .checked_add(audit.substeps)
            .ok_or(ScatteringRefinementError::Budget("substep-count overflow"))?;
        work.source_evaluations = work
            .source_evaluations
            .checked_add(source.source_evaluations())
            .ok_or(ScatteringRefinementError::Budget("source-count overflow"))?;
        work.traced_segments = work
            .traced_segments
            .checked_add(source.traced_segments())
            .ok_or(ScatteringRefinementError::Budget(
                "transport-count overflow",
            ))?;
        if work.substeps > c.maximum_total_substeps
            || work.source_evaluations > c.maximum_total_source_evaluations
            || work.traced_segments > c.maximum_total_traced_segments
        {
            return Err(ScatteringRefinementError::Budget(
                "complete-ensemble work exceeded",
            ));
        }
        if state.cells().len() > c.maximum_observable_cells {
            return Err(ScatteringRefinementError::Budget(
                "ensemble observable cells",
            ));
        }
        if state.time_s() != target || source.time_s() != target {
            return Err(ScatteringRefinementError::Invalid(
                "native ensemble endpoint mismatch",
            ));
        }
        trajectories.push(ScatteringFeedbackTrajectory {
            thermal: state,
            optical: source,
            audit,
        });
    }
    Ok(trajectories)
}

fn measure_ensemble<const N: usize>(
    initial: &DepositionSimulation<'_>,
    initial_optical: &TemperatureScatteringHeating<'_, N>,
    fine: &[ScatteringFeedbackTrajectory<'_, '_, N>],
    coarse: Option<&[ScatteringFeedbackTrajectory<'_, '_, N>]>,
    tolerance: ScatteringEnsembleTolerance,
    controls: ScatteringFeedbackRefinement,
    retain_cells: bool,
    checkpoint: Checkpoint<'_>,
) -> Result<(
    BTreeMap<CellId, ScatteringCellEstimate>,
    ScatteringEnsembleEnergy<N>,
    ScatteringEnsembleAssessment,
)> {
    checkpoint()?;
    if fine.len() != controls.replicas as usize || coarse.is_some_and(|c| c.len() != fine.len()) {
        return Err(ScatteringRefinementError::Invalid("replica-count mismatch"));
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
            return Err(ScatteringRefinementError::Invalid(
                "ensemble timeline/topology mismatch",
            ));
        }
    }
    let paired = coarse.is_some();
    let mut assessment = ScatteringEnsembleAssessment::default();
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
        assessment.maximum_normalized_width = assessment
            .maximum_normalized_width
            .max(normalize(
                h_width,
                h_scale,
                tolerance.enthalpy_absolute_j_kg,
                tolerance.relative,
            )?)
            .max(normalize(
                t_width,
                t_scale,
                tolerance.temperature_absolute_k,
                tolerance.relative,
            )?);
        if paired {
            // One binary64 epsilon of the represented native coordinate is a
            // conservative resolution floor, not a rounding-error certificate.
            // It prevents equal rounded h/T from claiming sub-ULP accuracy.
            let h_floor = (f64::EPSILON * h_magnitude).max(f64::from_bits(1));
            let t_floor = (f64::EPSILON * t_magnitude).max(f64::from_bits(1));
            assessment.maximum_coordinate_resolution_ratio = assessment
                .maximum_coordinate_resolution_ratio
                .max(normalize(
                    h_floor,
                    h_scale,
                    tolerance.enthalpy_absolute_j_kg,
                    tolerance.relative,
                )?)
                .max(normalize(
                    t_floor,
                    t_scale,
                    tolerance.temperature_absolute_k,
                    tolerance.relative,
                )?);
            assessment.maximum_normalized_width = assessment
                .maximum_normalized_width
                .max(assessment.maximum_coordinate_resolution_ratio);
        }
        if retain_cells {
            cells.insert(
                id,
                ScatteringCellEstimate {
                    specific_h_j_kg: h,
                    temperature_k: t,
                },
            );
        }
    }
    let mut energy = ScatteringEnsembleEnergy {
        incident_j: [ScatteringScalarEstimate::default(); N],
        ballistic_escape_j: [ScatteringScalarEstimate::default(); N],
        diffuse_escape_j: [ScatteringScalarEstimate::default(); N],
        absorbed_j: [ScatteringScalarEstimate::default(); N],
        unresolved_j: [ScatteringScalarEstimate::default(); N],
        residual_j: [ScatteringScalarEstimate::default(); N],
        thermal_input_j: ScatteringScalarEstimate::default(),
    };
    // All optical exit accounts matter, including untransported collision tails.
    // Compensated interval differences retain tiny new heat after large history.
    for lane in 0..N {
        for channel in 0..6 {
            let base = optical_sum(initial_optical, channel, lane);
            let (estimate, width, scale, _) = measure_scalar(
                fine.len(),
                Some(0.0),
                paired,
                controls.standard_error_multiplier,
                |r| {
                    Ok((
                        sum_difference(optical_sum(&fine[r].optical, channel, lane), base)?,
                        coarse
                            .map(|c| {
                                sum_difference(optical_sum(&c[r].optical, channel, lane), base)
                            })
                            .transpose()?,
                        coarse
                            .map(|c| {
                                sum_difference(
                                    optical_sum(&fine[r].optical, channel, lane),
                                    optical_sum(&c[r].optical, channel, lane),
                                )
                            })
                            .transpose()?,
                    ))
                },
                checkpoint,
            )?;
            match channel {
                0 => energy.incident_j[lane] = estimate,
                1 => energy.ballistic_escape_j[lane] = estimate,
                2 => energy.diffuse_escape_j[lane] = estimate,
                3 => energy.absorbed_j[lane] = estimate,
                4 => energy.unresolved_j[lane] = estimate,
                _ => energy.residual_j[lane] = estimate,
            }
            record_optical(&mut assessment, width, scale, tolerance)?;
        }
    }
    let (heat, width, scale, _) = measure_scalar(
        fine.len(),
        Some(0.0),
        paired,
        controls.standard_error_multiplier,
        |r| {
            Ok((
                sum_difference(fine[r].optical.sums.thermal, initial_optical.sums.thermal)?,
                coarse
                    .map(|c| {
                        sum_difference(c[r].optical.sums.thermal, initial_optical.sums.thermal)
                    })
                    .transpose()?,
                coarse
                    .map(|c| {
                        sum_difference(fine[r].optical.sums.thermal, c[r].optical.sums.thermal)
                    })
                    .transpose()?,
            ))
        },
        checkpoint,
    )?;
    energy.thermal_input_j = heat;
    record_optical(&mut assessment, width, scale, tolerance)?;
    Ok((cells, energy, assessment))
}

fn optical_sum<const N: usize>(
    o: &TemperatureScatteringHeating<'_, N>,
    channel: usize,
    lane: usize,
) -> CheckedSum {
    match channel {
        0 => o.sums.incident[lane],
        1 => o.sums.ballistic_escape[lane],
        2 => o.sums.diffuse_escape[lane],
        3 => o.sums.absorbed[lane],
        4 => o.sums.unresolved[lane],
        _ => o.sums.residual[lane],
    }
}

fn measure_scalar(
    count: usize,
    base: Option<f64>,
    paired: bool,
    multiplier: f64,
    mut values: impl FnMut(usize) -> Result<(f64, Option<f64>, Option<f64>)>,
    checkpoint: Checkpoint<'_>,
) -> Result<(ScatteringScalarEstimate, f64, f64, f64)> {
    let mut fine = Moments::default();
    let mut coarse = Moments::default();
    let mut difference = Moments::default();
    let origin = base.unwrap_or(0.0);
    let mut magnitude = origin.abs();
    for replica in 0..count {
        checkpoint()?;
        let (a, b, paired_difference) = values(replica)?;
        magnitude = magnitude.max(a.abs());
        fine.observe(finite(a - origin)?)?;
        if let Some(b) = b {
            magnitude = magnitude.max(b.abs());
            coarse.observe(finite(b - origin)?)?;
            difference.observe(match paired_difference {
                Some(value) => value,
                None => finite(a - b)?,
            })?;
        }
    }
    let mut estimate = fine.estimate()?;
    let scale = if base.is_some() {
        estimate.mean.abs().max(if paired {
            coarse.estimate()?.mean.abs()
        } else {
            0.0
        })
    } else {
        0.0
    };
    let width = if paired {
        let delta = difference.estimate()?;
        finite(delta.mean.abs() + multiplier * delta.standard_error)?
    } else {
        finite(multiplier * estimate.standard_error)?
    };
    estimate.mean = finite(origin + estimate.mean)?;
    Ok((estimate, width, scale, magnitude))
}

/// Signed, centered and scaled Welford moments. Centering avoids losing thermal
/// variability against an arbitrary caloric datum; scaling avoids squaring SI
/// magnitudes that would overflow or underflow the second central moment.
#[derive(Clone, Copy, Default)]
struct Moments {
    count: usize,
    anchor: f64,
    scale: f64,
    mean: f64,
    m2: f64,
}
impl Moments {
    fn observe(&mut self, value: f64) -> Result<()> {
        finite(value)?;
        if self.count == 0 {
            self.anchor = value;
        }
        let shifted = finite(value - self.anchor)?;
        if shifted.abs() > self.scale {
            let factor = self.scale / shifted.abs();
            self.mean *= factor;
            self.m2 *= factor * factor;
            self.scale = shifted.abs();
        }
        let normalized = if self.scale == 0.0 {
            0.0
        } else {
            shifted / self.scale
        };
        self.count += 1;
        let delta = normalized - self.mean;
        self.mean += delta / self.count as f64;
        self.m2 += delta * (normalized - self.mean);
        Ok(())
    }
    fn estimate(self) -> Result<ScatteringScalarEstimate> {
        if self.count < 2 {
            return Err(ScatteringRefinementError::Invalid(
                "independent replicas required",
            ));
        }
        Ok(ScatteringScalarEstimate {
            mean: finite(self.anchor + self.mean * self.scale)?,
            standard_error: finite(
                self.scale
                    * det::sqrt(self.m2.max(0.0) / self.count as f64 / (self.count - 1) as f64),
            )?,
        })
    }
}
fn sum_difference(a: CheckedSum, b: CheckedSum) -> Result<f64> {
    let mut difference = CheckedSum::default();
    difference
        .add(a.sum - b.sum)
        .map_err(OpticalFeedbackError::from)?;
    difference
        .add(a.correction)
        .map_err(OpticalFeedbackError::from)?;
    difference
        .add(-b.correction)
        .map_err(OpticalFeedbackError::from)?;
    Ok(difference.value().map_err(OpticalFeedbackError::from)?)
}
fn record_optical(
    a: &mut ScatteringEnsembleAssessment,
    width: f64,
    scale: f64,
    t: ScatteringEnsembleTolerance,
) -> Result<()> {
    a.maximum_optical_energy_width_j = a.maximum_optical_energy_width_j.max(width);
    a.maximum_normalized_width = a.maximum_normalized_width.max(normalize(
        width,
        scale,
        t.optical_energy_absolute_j,
        t.relative,
    )?);
    Ok(())
}
fn combine(
    a: ScatteringEnsembleAssessment,
    b: ScatteringEnsembleAssessment,
) -> ScatteringEnsembleAssessment {
    ScatteringEnsembleAssessment {
        maximum_enthalpy_width_j_kg: a
            .maximum_enthalpy_width_j_kg
            .max(b.maximum_enthalpy_width_j_kg),
        maximum_temperature_width_k: a
            .maximum_temperature_width_k
            .max(b.maximum_temperature_width_k),
        maximum_optical_energy_width_j: a
            .maximum_optical_energy_width_j
            .max(b.maximum_optical_energy_width_j),
        maximum_normalized_width: a.maximum_normalized_width.max(b.maximum_normalized_width),
        maximum_coordinate_resolution_ratio: a
            .maximum_coordinate_resolution_ratio
            .max(b.maximum_coordinate_resolution_ratio),
    }
}
fn normalize(width: f64, scale: f64, absolute: f64, relative: f64) -> Result<f64> {
    // An overflowing ratio means failure of the gate, not successful admission.
    let allowance = finite(absolute + relative * scale)?;
    Ok(width / allowance)
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ScatteringRefinementError::Invalid(
            "unrepresentable ensemble statistic",
        ))
    }
}
