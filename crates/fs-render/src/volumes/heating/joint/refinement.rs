//! Complete native trajectories driven jointly by incident illumination and
//! temperature-dependent thermal radiation, with independent spectral packets.
//! Photon response fields are common within a replica and independent between
//! replicas. Each photon-sensitivity comparison lowers only one owner's history
//! count, so opposing illumination and cooling errors cannot cancel that gate.
//!
//! These are empirical time, sampling and finite-history diagnostics, conditional
//! on the accepted input state. They do not certify confidence coverage, rare
//! paths, optional stopping, finite-history bias, spatial source quadrature or
//! constitutive/spectral accuracy.

use super::super::{
    CheckedSum,
    ensemble_statistics::{self, StatisticsError},
    exchange::feedback::refinement::{RadiativeEnsembleEnergy, RadiativeScalarEstimate},
    scattering_feedback::{ScatteringEnsembleEnergy, ScatteringScalarEstimate},
};
use super::{JointOpticalError, JointOpticalHeating};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionSimulation, EnergyAudit, StepConfig,
};
use fs_thermal::observers::CoupledObservers;
use std::collections::BTreeMap;

/// Absolute physical allowances and a relative scale measured from the input.
/// Formation enthalpy never enlarges an allowance. Separate component and net
/// heat tolerances preserve meaningful checks close to a radiative balance.
#[derive(Clone, Copy, Debug)]
pub struct JointOpticalTolerance {
    pub enthalpy_absolute_j_kg: f64,
    pub temperature_absolute_k: f64,
    pub illumination_energy_absolute_j: f64,
    pub radiation_energy_absolute_j: f64,
    pub net_optical_heat_absolute_j: f64,
    pub relative: f64,
}

/// Work and storage ceilings apply to the complete study, including discarded
/// trials. Each physical owner has its own work allowance; one remaining photon
/// budget is never granted twice to the coupled source evaluation.
#[derive(Clone, Copy, Debug)]
pub struct JointOpticalRefinement {
    pub initial_max_dt_s: f64,
    pub minimum_max_dt_s: f64,
    /// At least two: two consecutive actual time refinements must pass.
    pub maximum_halvings: usize,
    /// Fixed independent replicas at every resolution, at least two.
    pub replicas: u32,
    pub ensemble_seed: u64,
    /// At least one. Emitted radiation remains sampled without scattering.
    pub maximum_photon_doublings: u32,
    /// Empirical SE multiplier; it does not specify a coverage probability.
    pub standard_error_multiplier: f64,
    pub temporal: JointOpticalTolerance,
    pub sampling: JointOpticalTolerance,
    pub illumination_photon_sensitivity: JointOpticalTolerance,
    pub radiation_photon_sensitivity: JointOpticalTolerance,
    pub conservation_tolerance_j: f64,
    pub maximum_total_substeps: usize,
    pub maximum_total_illumination_source_evaluations: usize,
    pub maximum_total_radiation_source_evaluations: usize,
    pub maximum_total_illumination_traced_segments: u64,
    pub maximum_total_radiation_traced_segments: u64,
    /// At least 2*replicas, including an in-flight native thermal candidate.
    pub maximum_retained_trajectories: usize,
    pub maximum_observable_cells: usize,
    /// Conservative (2*replicas+4)*radiation.max_pairs. Slots cover optical and
    /// mapped pair representations, including replacement/construction overlap.
    pub maximum_retained_operator_pair_slots: usize,
    /// The same peak factor times BOTH geometry cell counts plus radiation's
    /// material/emitter limits. Vacuum cells retain coefficient packets too.
    pub maximum_retained_operator_cell_slots: usize,
    /// The same factor times (illumination.max_retained_deposition_cells +
    /// illumination.source_count). Source headers count even when no cell heats.
    /// All storage caps count typed slots, not allocator bytes.
    pub maximum_retained_illumination_report_slots: usize,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct JointScalarEstimate {
    pub mean: f64,
    pub standard_error: f64,
}

/// Temperature is measured on each native trajectory before averaging T(h).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JointCellEstimate {
    pub specific_h_j_kg: JointScalarEstimate,
    pub temperature_k: JointScalarEstimate,
}

/// Interval increments, retaining the covariance of both radiation owners on
/// each complete trajectory. Internal exchange is not external heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JointOpticalEnergy<const NI: usize, const NR: usize> {
    pub illumination: ScatteringEnsembleEnergy<NI>,
    pub radiation: RadiativeEnsembleEnergy<NR>,
    /// Absorbed incident heat plus signed bath-to-material heat, combined in
    /// compensated arithmetic before rounding or estimating cross-replica SE.
    pub net_optical_heat_j: JointScalarEstimate,
}

/// Paired widths are |mean(fine-coarse)| + k*SE(fine-coarse); the final sampling
/// gate uses k*SE(fine). All interpretations remain empirical after selection.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct JointOpticalAssessment {
    pub maximum_enthalpy_width_j_kg: f64,
    pub maximum_temperature_width_k: f64,
    pub maximum_illumination_energy_width_j: f64,
    pub maximum_radiation_energy_width_j: f64,
    pub maximum_net_optical_heat_width_j: f64,
    pub maximum_normalized_width: f64,
    /// Paired h/T comparisons refuse accuracy below native coordinate resolution.
    pub maximum_coordinate_resolution_ratio: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct JointOpticalWork {
    pub substeps: usize,
    pub illumination_source_evaluations: usize,
    pub radiation_source_evaluations: usize,
    pub illumination_traced_segments: u64,
    pub radiation_traced_segments: u64,
    pub illumination_transport_builds: usize,
    pub radiation_operator_builds: usize,
}

/// A missing gate was skipped because an earlier gate failed, or because
/// illumination uses the exact no-scattering path. Radiation H/2H is mandatory.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct JointOpticalDiagnostics {
    pub temporal: Option<JointOpticalAssessment>,
    pub temporal_refinement_observed: bool,
    pub sampling: Option<JointOpticalAssessment>,
    pub illumination_photon_sensitivity: Option<JointOpticalAssessment>,
    pub radiation_photon_sensitivity: Option<JointOpticalAssessment>,
    pub maximum_dt_s: f64,
    pub halvings: usize,
    pub photon_doublings: u32,
    pub illumination_histories_per_lane: u32,
    pub radiation_histories_per_emitter_lane: u32,
    pub work: JointOpticalWork,
    pub admitted_peak_operator_pair_slots: usize,
    pub admitted_peak_operator_cell_slots: usize,
    pub admitted_peak_illumination_report_slots: usize,
}

/// An actual accepted native trajectory, including both optical histories.
pub struct JointOpticalTrajectory<'t, 'g, const NI: usize, const NR: usize> {
    pub thermal: DepositionSimulation<'t>,
    pub optical: JointOpticalHeating<'g, NI, NR>,
    pub audit: EnergyAudit,
}

pub struct JointOpticalEnsemble<'t, 'g, const NI: usize, const NR: usize> {
    pub trajectories: Vec<JointOpticalTrajectory<'t, 'g, NI, NR>>,
    pub cells: BTreeMap<CellId, JointCellEstimate>,
    pub energy: JointOpticalEnergy<NI, NR>,
    pub diagnostics: JointOpticalDiagnostics,
}

#[derive(Clone, Debug, PartialEq)]
pub enum JointOpticalRefinementError {
    Invalid(&'static str),
    Budget(&'static str),
    Physics(JointOpticalError),
    Unresolved(JointOpticalDiagnostics),
}
impl From<JointOpticalError> for JointOpticalRefinementError {
    fn from(error: JointOpticalError) -> Self {
        Self::Physics(error)
    }
}
impl From<DepositionError> for JointOpticalRefinementError {
    fn from(error: DepositionError) -> Self {
        Self::Physics(error.into())
    }
}
impl From<StatisticsError> for JointOpticalRefinementError {
    fn from(error: StatisticsError) -> Self {
        match error {
            StatisticsError::Invalid(message) => Self::Invalid(message),
            StatisticsError::Heating(error) => {
                Self::Physics(JointOpticalError::First(error.into()))
            }
        }
    }
}
impl core::fmt::Display for JointOpticalRefinementError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "joint optical refinement refused: {self:?}")
    }
}
impl core::error::Error for JointOpticalRefinementError {}

type Result<T> = core::result::Result<T, JointOpticalRefinementError>;
type Checkpoint<'a> = &'a mut dyn FnMut() -> core::result::Result<(), DepositionError>;

impl<'g, const NI: usize, const NR: usize> JointOpticalHeating<'g, NI, NR> {
    /// Assess the SAME coupled native evolution with independent optical
    /// replicas. Incident and thermal wavelength packets remain independent;
    /// supported visible data never creates missing infrared material data.
    ///
    /// All trajectories start from the accepted input. Photon streams remain
    /// common within a replica across time/temperature/history resolutions.
    /// Two real time refinements, final sampling, and separate one-owner-at-a-
    /// time photon sensitivities must pass. Physics/domain/work/cancellation
    /// failures abort transactionally; only unresolved gates trigger refinement.
    pub fn advance_ensemble_refined<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        controls: JointOpticalRefinement,
        mut checkpoint: impl FnMut() -> core::result::Result<(), DepositionError>,
    ) -> Result<JointOpticalEnsemble<'t, 'g, NI, NR>> {
        checkpoint()?;
        let (pair_slots, cell_slots, report_slots) = validate(thermal, self, target_s, controls)?;
        let mut work = JointOpticalWork::default();
        let mut diagnostics = JointOpticalDiagnostics::default();
        let mut start_halving = 0;
        let stochastic_illumination = self.illumination().has_scattering();
        for doubling in 1..=controls.maximum_photon_doublings {
            checkpoint()?;
            let factor = 1_u32
                .checked_shl(doubling)
                .ok_or(budget("photon-count overflow"))?;
            let illumination_histories = self
                .illumination()
                .histories_per_lane()
                .checked_mul(if stochastic_illumination { factor } else { 1 })
                .ok_or(budget("illumination photon-count overflow"))?;
            let radiation_histories = self
                .radiation()
                .histories_per_emitter_lane()
                .checked_mul(factor)
                .ok_or(budget("radiation photon-count overflow"))?;
            let histories = [illumination_histories, radiation_histories];
            let mut previous: Option<Vec<JointOpticalTrajectory<'t, 'g, NI, NR>>> = None;
            let mut previous_time: Option<JointOpticalAssessment> = None;
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
                diagnostics = JointOpticalDiagnostics {
                    maximum_dt_s: dt,
                    halvings: halving,
                    photon_doublings: doubling,
                    illumination_histories_per_lane: illumination_histories,
                    radiation_histories_per_emitter_lane: radiation_histories,
                    work,
                    admitted_peak_operator_pair_slots: pair_slots,
                    admitted_peak_operator_cell_slots: cell_slots,
                    admitted_peak_illumination_report_slots: report_slots,
                    ..JointOpticalDiagnostics::default()
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
                                let mut photon_passed = true;
                                for component in 0..2 {
                                    if component == 0 && !stochastic_illumination {
                                        continue;
                                    }
                                    let mut low_histories = histories;
                                    low_histories[component] /= 2;
                                    let low = run_ensemble(
                                        self,
                                        thermal,
                                        target_s,
                                        dt,
                                        low_histories,
                                        controls,
                                        &mut work,
                                        &mut checkpoint,
                                    )?;
                                    let tolerance = if component == 0 {
                                        controls.illumination_photon_sensitivity
                                    } else {
                                        controls.radiation_photon_sensitivity
                                    };
                                    let (_, _, photon) = measure_ensemble(
                                        thermal,
                                        self,
                                        &fine,
                                        Some(&low),
                                        tolerance,
                                        controls,
                                        false,
                                        &mut checkpoint,
                                    )?;
                                    // Release low illumination before allocating low radiation.
                                    drop(low);
                                    if component == 0 {
                                        diagnostics.illumination_photon_sensitivity = Some(photon);
                                    } else {
                                        diagnostics.radiation_photon_sensitivity = Some(photon);
                                    }
                                    diagnostics.work = work;
                                    if photon.maximum_normalized_width > 1.0 {
                                        photon_passed = false;
                                        break;
                                    }
                                }
                                if photon_passed {
                                    for candidate in &mut fine {
                                        candidate
                                            .optical
                                            .observers
                                            .first
                                            .restore_budget(self.illumination().budget());
                                        candidate
                                            .optical
                                            .observers
                                            .second
                                            .restore_budget(self.radiation().budget());
                                    }
                                    checkpoint()?;
                                    return Ok(JointOpticalEnsemble {
                                        trajectories: fine,
                                        cells,
                                        energy,
                                        diagnostics,
                                    });
                                }
                            }
                            start_halving = halving.saturating_sub(2);
                            break;
                        }
                    }
                    previous_time = refined.then_some(time);
                }
                previous = Some(fine);
                let next = dt * 0.5;
                if next <= 0.0 || next == dt {
                    break;
                }
                dt = next;
            }
        }
        Err(JointOpticalRefinementError::Unresolved(diagnostics))
    }
}

fn validate<const NI: usize, const NR: usize>(
    thermal: &DepositionSimulation<'_>,
    optical: &JointOpticalHeating<'_, NI, NR>,
    target: f64,
    c: JointOpticalRefinement,
) -> Result<(usize, usize, usize)> {
    if optical.illumination().time_s() != thermal.time_s()
        || optical.radiation().time_s() != thermal.time_s()
    {
        return Err(invalid("joint optical and native thermal clocks differ"));
    }
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
        || ![
            c.temporal,
            c.sampling,
            c.illumination_photon_sensitivity,
            c.radiation_photon_sensitivity,
        ]
        .into_iter()
        .all(|t| {
            positive(t.enthalpy_absolute_j_kg)
                && positive(t.temperature_absolute_k)
                && positive(t.illumination_energy_absolute_j)
                && positive(t.radiation_energy_absolute_j)
                && positive(t.net_optical_heat_absolute_j)
                && t.relative.is_finite()
                && (0.0..1.0).contains(&t.relative)
        })
    {
        return Err(invalid(
            "invalid joint optical ensemble, timeline or accuracy controls",
        ));
    }
    let retained = (c.replicas as usize)
        .checked_mul(2)
        .ok_or(budget("trajectory-count overflow"))?;
    let peak = retained
        .checked_add(4)
        .ok_or(budget("operator-count overflow"))?;
    let config = optical.radiation().sampling_config();
    let pairs = peak
        .checked_mul(config.max_pairs)
        .ok_or(budget("operator pair-slot overflow"))?;
    let cells = config
        .max_material_cells
        .checked_add(config.max_emitter_cells)
        .and_then(|v| v.checked_add(optical.radiation().field().geometry().cell_count()))
        .and_then(|v| {
            v.checked_add(
                optical
                    .illumination()
                    .optical_field()
                    .geometry()
                    .cell_count(),
            )
        })
        .and_then(|v| peak.checked_mul(v))
        .ok_or(budget("operator cell-slot overflow"))?;
    let reports = optical
        .illumination()
        .budget()
        .max_retained_deposition_cells
        .checked_add(optical.illumination().source_count())
        .and_then(|v| peak.checked_mul(v))
        .ok_or(budget("illumination report-slot overflow"))?;
    if retained > c.maximum_retained_trajectories
        || pairs > c.maximum_retained_operator_pair_slots
        || cells > c.maximum_retained_operator_cell_slots
        || reports > c.maximum_retained_illumination_report_slots
        || c.maximum_observable_cells == 0
        || thermal.cells().len() > c.maximum_observable_cells
        || c.maximum_total_substeps == 0
        || c.maximum_total_illumination_source_evaluations == 0
        || c.maximum_total_radiation_source_evaluations == 0
    {
        return Err(budget(
            "joint optical ensemble retained-storage/work allowance",
        ));
    }
    Ok((pairs, cells, reports))
}

#[allow(clippy::too_many_arguments)]
fn run_ensemble<'t, 'g, const NI: usize, const NR: usize>(
    optical: &JointOpticalHeating<'g, NI, NR>,
    thermal: &DepositionSimulation<'t>,
    target: f64,
    dt: f64,
    histories: [u32; 2],
    c: JointOpticalRefinement,
    work: &mut JointOpticalWork,
    checkpoint: Checkpoint<'_>,
) -> Result<Vec<JointOpticalTrajectory<'t, 'g, NI, NR>>> {
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
        let illumination_sources = c
            .maximum_total_illumination_source_evaluations
            .checked_sub(work.illumination_source_evaluations)
            .ok_or(budget("illumination source evaluations"))?;
        let radiation_sources = c
            .maximum_total_radiation_source_evaluations
            .checked_sub(work.radiation_source_evaluations)
            .ok_or(budget("radiation source evaluations"))?;
        let illumination_segments = c
            .maximum_total_illumination_traced_segments
            .checked_sub(work.illumination_traced_segments)
            .ok_or(budget("illumination transport segments"))?;
        let radiation_segments = c
            .maximum_total_radiation_traced_segments
            .checked_sub(work.radiation_traced_segments)
            .ok_or(budget("radiation transport segments"))?;
        if steps == 0 || illumination_sources == 0 || radiation_sources == 0 {
            return Err(budget("complete joint optical ensemble work exhausted"));
        }
        let seed = c.ensemble_seed.wrapping_add(u64::from(replica));
        let illumination = optical
            .illumination()
            .configured_trial(
                seed,
                histories[0],
                illumination_sources,
                illumination_segments,
            )
            .map_err(JointOpticalError::First)?;
        let radiation = optical
            .radiation()
            .configured_trial(seed, histories[1], radiation_sources, radiation_segments)
            .map_err(JointOpticalError::Second)?;
        let trial = JointOpticalHeating {
            observers: CoupledObservers::new(illumination, radiation),
        };
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
        work.substeps = add_count(work.substeps, audit.substeps)?;
        work.illumination_source_evaluations = add_count(
            work.illumination_source_evaluations,
            source.illumination().source_evaluations(),
        )?;
        work.radiation_source_evaluations = add_count(
            work.radiation_source_evaluations,
            source.radiation().source_evaluations(),
        )?;
        work.radiation_operator_builds = add_count(
            work.radiation_operator_builds,
            source.radiation().operator_builds(),
        )?;
        work.illumination_transport_builds = add_count(
            work.illumination_transport_builds,
            source.illumination().transport_builds(),
        )?;
        work.illumination_traced_segments = work
            .illumination_traced_segments
            .checked_add(source.illumination().traced_segments())
            .ok_or(budget("segment-count overflow"))?;
        work.radiation_traced_segments = work
            .radiation_traced_segments
            .checked_add(source.radiation().traced_segments())
            .ok_or(budget("segment-count overflow"))?;
        if work.substeps > c.maximum_total_substeps
            || work.illumination_source_evaluations
                > c.maximum_total_illumination_source_evaluations
            || work.radiation_source_evaluations > c.maximum_total_radiation_source_evaluations
            || work.illumination_traced_segments > c.maximum_total_illumination_traced_segments
            || work.radiation_traced_segments > c.maximum_total_radiation_traced_segments
        {
            return Err(budget("complete joint optical ensemble work exceeded"));
        }
        if state.cells().len() > c.maximum_observable_cells {
            return Err(budget("ensemble observable cells"));
        }
        if state.time_s() != target
            || source.illumination().time_s() != target
            || source.radiation().time_s() != target
        {
            return Err(invalid("native joint optical ensemble endpoint mismatch"));
        }
        trajectories.push(JointOpticalTrajectory {
            thermal: state,
            optical: source,
            audit,
        });
    }
    Ok(trajectories)
}

#[allow(clippy::too_many_arguments)]
fn measure_ensemble<const NI: usize, const NR: usize>(
    initial: &DepositionSimulation<'_>,
    initial_optical: &JointOpticalHeating<'_, NI, NR>,
    fine: &[JointOpticalTrajectory<'_, '_, NI, NR>],
    coarse: Option<&[JointOpticalTrajectory<'_, '_, NI, NR>]>,
    tolerance: JointOpticalTolerance,
    controls: JointOpticalRefinement,
    retain_cells: bool,
    checkpoint: Checkpoint<'_>,
) -> Result<(
    BTreeMap<CellId, JointCellEstimate>,
    JointOpticalEnergy<NI, NR>,
    JointOpticalAssessment,
)> {
    checkpoint()?;
    if fine.len() != controls.replicas as usize || coarse.is_some_and(|c| c.len() != fine.len()) {
        return Err(invalid("joint optical replica-count mismatch"));
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
            return Err(invalid("joint optical ensemble timeline/topology mismatch"));
        }
    }
    let paired = coarse.is_some();
    let mut assessment = JointOpticalAssessment::default();
    let mut cells = BTreeMap::new();
    for &id in reference.cells().keys() {
        checkpoint()?;
        let initial_h = initial.cells().get(&id).map(|c| c.specific_h_j_kg);
        let initial_t = if initial_h.is_some() {
            Some(initial.temperature_k(id)?)
        } else {
            None
        };
        let (h, hw, hs, hm) = measure_scalar(
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
        let (t, tw, ts, tm) = measure_scalar(
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
        assessment.maximum_enthalpy_width_j_kg = assessment.maximum_enthalpy_width_j_kg.max(hw);
        assessment.maximum_temperature_width_k = assessment.maximum_temperature_width_k.max(tw);
        for (width, scale, magnitude, absolute) in [
            (hw, hs, hm, tolerance.enthalpy_absolute_j_kg),
            (tw, ts, tm, tolerance.temperature_absolute_k),
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
                JointCellEstimate {
                    specific_h_j_kg: h,
                    temperature_k: t,
                },
            );
        }
    }
    let mut energy = JointOpticalEnergy {
        illumination: ScatteringEnsembleEnergy {
            incident_j: [ScatteringScalarEstimate::default(); NI],
            ballistic_escape_j: [ScatteringScalarEstimate::default(); NI],
            diffuse_escape_j: [ScatteringScalarEstimate::default(); NI],
            absorbed_j: [ScatteringScalarEstimate::default(); NI],
            unresolved_j: [ScatteringScalarEstimate::default(); NI],
            residual_j: [ScatteringScalarEstimate::default(); NI],
            thermal_input_j: ScatteringScalarEstimate::default(),
        },
        radiation: RadiativeEnsembleEnergy {
            net_bath_to_material_j: [RadiativeScalarEstimate::default(); NR],
            internal_transfer_magnitude_j: [RadiativeScalarEstimate::default(); NR],
            thermal_external_j: RadiativeScalarEstimate::default(),
            unresolved_exchange_scale_j: RadiativeScalarEstimate::default(),
        },
        net_optical_heat_j: JointScalarEstimate::default(),
    };
    // Component accounts are gated separately even when their net heat cancels.
    for channel in 0..12 {
        let lanes = if channel < 6 {
            NI
        } else if (7..9).contains(&channel) {
            NR
        } else {
            1
        };
        for lane in 0..lanes {
            let (estimate, width, scale, _) = measure_scalar(
                fine.len(),
                Some(0.0),
                paired,
                controls.standard_error_multiplier,
                |r| {
                    Ok((
                        energy_difference(&fine[r].optical, initial_optical, channel, lane)?,
                        coarse
                            .map(|c| {
                                energy_difference(&c[r].optical, initial_optical, channel, lane)
                            })
                            .transpose()?,
                        coarse
                            .map(|c| {
                                energy_difference(&fine[r].optical, &c[r].optical, channel, lane)
                            })
                            .transpose()?,
                    ))
                },
                checkpoint,
            )?;
            let absolute = if channel < 7 {
                let value = ScatteringScalarEstimate {
                    mean: estimate.mean,
                    standard_error: estimate.standard_error,
                };
                match channel {
                    0 => energy.illumination.incident_j[lane] = value,
                    1 => energy.illumination.ballistic_escape_j[lane] = value,
                    2 => energy.illumination.diffuse_escape_j[lane] = value,
                    3 => energy.illumination.absorbed_j[lane] = value,
                    4 => energy.illumination.unresolved_j[lane] = value,
                    5 => energy.illumination.residual_j[lane] = value,
                    _ => energy.illumination.thermal_input_j = value,
                }
                assessment.maximum_illumination_energy_width_j =
                    assessment.maximum_illumination_energy_width_j.max(width);
                tolerance.illumination_energy_absolute_j
            } else if channel < 11 {
                let value = RadiativeScalarEstimate {
                    mean: estimate.mean,
                    standard_error: estimate.standard_error,
                };
                match channel {
                    7 => energy.radiation.net_bath_to_material_j[lane] = value,
                    8 => energy.radiation.internal_transfer_magnitude_j[lane] = value,
                    9 => energy.radiation.thermal_external_j = value,
                    _ => energy.radiation.unresolved_exchange_scale_j = value,
                }
                assessment.maximum_radiation_energy_width_j =
                    assessment.maximum_radiation_energy_width_j.max(width);
                tolerance.radiation_energy_absolute_j
            } else {
                energy.net_optical_heat_j = estimate;
                assessment.maximum_net_optical_heat_width_j = width;
                tolerance.net_optical_heat_absolute_j
            };
            assessment.maximum_normalized_width = assessment
                .maximum_normalized_width
                .max(normalize(width, scale, absolute, tolerance.relative)?);
        }
    }
    Ok((cells, energy, assessment))
}

fn energy_difference<const NI: usize, const NR: usize>(
    a: &JointOpticalHeating<'_, NI, NR>,
    b: &JointOpticalHeating<'_, NI, NR>,
    channel: usize,
    lane: usize,
) -> Result<f64> {
    if channel < 7 {
        Ok(ensemble_statistics::sum_difference(
            a.illumination().energy_sum(channel, lane),
            b.illumination().energy_sum(channel, lane),
        )?)
    } else if channel < 11 {
        Ok(ensemble_statistics::sum_difference(
            a.radiation().energy_sum(channel - 7, lane),
            b.radiation().energy_sum(channel - 7, lane),
        )?)
    } else {
        combined_difference([
            (
                a.illumination().energy_sum(6, 0),
                b.illumination().energy_sum(6, 0),
            ),
            (
                a.radiation().energy_sum(2, 0),
                b.radiation().energy_sum(2, 0),
            ),
        ])
    }
}

// Retain opposing component corrections until the complete net interval is
// reduced; rounding two large component intervals first would erase net heat.
fn combined_difference(accounts: [(CheckedSum, CheckedSum); 2]) -> Result<f64> {
    let mut difference = CheckedSum::default();
    for (a, b) in accounts {
        difference
            .add(a.sum - b.sum)
            .map_err(StatisticsError::from)?;
        difference
            .add(a.correction)
            .map_err(StatisticsError::from)?;
        difference
            .add(-b.correction)
            .map_err(StatisticsError::from)?;
    }
    Ok(difference.value().map_err(StatisticsError::from)?)
}

fn measure_scalar(
    count: usize,
    base: Option<f64>,
    paired: bool,
    multiplier: f64,
    values: impl FnMut(usize) -> Result<(f64, Option<f64>, Option<f64>)>,
    checkpoint: Checkpoint<'_>,
) -> Result<(JointScalarEstimate, f64, f64, f64)> {
    let (value, width, scale, magnitude) =
        ensemble_statistics::measure_scalar(count, base, paired, multiplier, values, checkpoint)?;
    Ok((
        JointScalarEstimate {
            mean: value.mean,
            standard_error: value.standard_error,
        },
        width,
        scale,
        magnitude,
    ))
}
fn combine(a: JointOpticalAssessment, b: JointOpticalAssessment) -> JointOpticalAssessment {
    JointOpticalAssessment {
        maximum_enthalpy_width_j_kg: a
            .maximum_enthalpy_width_j_kg
            .max(b.maximum_enthalpy_width_j_kg),
        maximum_temperature_width_k: a
            .maximum_temperature_width_k
            .max(b.maximum_temperature_width_k),
        maximum_illumination_energy_width_j: a
            .maximum_illumination_energy_width_j
            .max(b.maximum_illumination_energy_width_j),
        maximum_radiation_energy_width_j: a
            .maximum_radiation_energy_width_j
            .max(b.maximum_radiation_energy_width_j),
        maximum_net_optical_heat_width_j: a
            .maximum_net_optical_heat_width_j
            .max(b.maximum_net_optical_heat_width_j),
        maximum_normalized_width: a.maximum_normalized_width.max(b.maximum_normalized_width),
        maximum_coordinate_resolution_ratio: a
            .maximum_coordinate_resolution_ratio
            .max(b.maximum_coordinate_resolution_ratio),
    }
}
fn normalize(width: f64, scale: f64, absolute: f64, relative: f64) -> Result<f64> {
    let allowance = absolute + relative * scale;
    if !allowance.is_finite() {
        return Err(invalid("unrepresentable joint ensemble allowance"));
    }
    Ok(width / allowance)
}
fn add_count(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b).ok_or(budget("work-count overflow"))
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn invalid(message: &'static str) -> JointOpticalRefinementError {
    JointOpticalRefinementError::Invalid(message)
}
fn budget(message: &'static str) -> JointOpticalRefinementError {
    JointOpticalRefinementError::Budget(message)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn opposing_heat_retains_small_interval_and_paired_corrections() {
        let mut old_light = CheckedSum::default();
        let mut old_radiation = CheckedSum::default();
        old_light.add(1e25).unwrap();
        old_radiation.add(-1e25).unwrap();
        let mut light = old_light;
        let mut radiation = old_radiation;
        light.add(1e16).unwrap();
        light.add(0.25).unwrap();
        radiation.add(-1e16).unwrap();
        radiation.add(-0.125).unwrap();
        assert_eq!(
            combined_difference([(light, old_light), (radiation, old_radiation)]).unwrap(),
            0.125
        );
        let mut coarse_light = light;
        let mut coarse_radiation = radiation;
        coarse_light.add(-0.0625).unwrap();
        coarse_radiation.add(0.03125).unwrap();
        assert_eq!(
            combined_difference([(light, coarse_light), (radiation, coarse_radiation)]).unwrap(),
            0.03125
        );
    }
}
