//! Local temporal accuracy control for the original conservative surface maps.
//!
//! One palindromic sweep is compared with two half sweeps. The second-order
//! fine-path error estimate is one third of their concentration discrepancy
//! for equal halves; actual unequal representable halves use the corresponding
//! cubic correction. Only the unextrapolated fine state can be accepted, so
//! positivity, partition equilibrium and the original free-energy law remain
//! those of the existing pair maps. This is a local splitting-error estimate,
//! not a certified global or spatial error bound.
//!
//! Coarse and rejected work consume the same explicit budgets. Three amount
//! buffers are reused across trials; the original bounded diffusion kernel
//! owns each sweep and its physical admission. Whole-window amount and energy
//! allowances are checked again against the original inventory before atomic
//! publication. Source, constitutive and physical-balance refusals propagate;
//! only excessive temporal error is retried. Geometry and laws remain held.

use super::{
    Material, Sum, SurfaceDiffusionControls, SurfaceDiffusionEdge, SurfaceDiffusionError,
    SurfaceDiffusionReport, SurfacePartitionDiffusionEdge, SurfactantFace, SurfactantLaw, det,
    diffuse, finite, poll, poll_at,
};

/// Explicit time-error, physical, storage and whole-call work limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveSurfaceDiffusionControls {
    pub initial_step_s: f64,
    /// A final remainder below this floor may be tried once; it still must
    /// satisfy the same error and physical gates.
    pub minimum_step_s: f64,
    /// Upper proposed comparison span; each accepted span has two fine sweeps.
    pub maximum_step_s: f64,
    /// Positive per-face absolute tolerance in physical concentration [mol/m²].
    pub absolute_concentration_tolerance_mol_m2: f64,
    /// Relative tolerance in [0,1), scaled by max(start, fine) concentration.
    pub relative_tolerance: f64,
    /// All full-versus-two-half comparisons, including accuracy rejections.
    pub maximum_trials: usize,
    /// All pair applications, including coarse and rejected trajectories.
    /// A nonempty E-edge comparison costs exactly 3*(2*E-1) applications.
    pub maximum_pair_exchanges: usize,
    pub maximum_faces: usize,
    pub maximum_edges: usize,
    /// The original amount allowance, applied to each sweep AND the whole call.
    pub amount_tolerance_mol: f64,
    /// The original free-energy allowance, applied per sweep AND to the call.
    pub energy_increase_tolerance_j: f64,
}

impl AdaptiveSurfaceDiffusionControls {
    fn step(self, duration_s: f64) -> SurfaceDiffusionControls {
        SurfaceDiffusionControls {
            maximum_step_s: duration_s,
            maximum_substeps: 1,
            maximum_faces: self.maximum_faces,
            maximum_edges: self.maximum_edges,
            amount_tolerance_mol: self.amount_tolerance_mol,
            energy_increase_tolerance_j: self.energy_increase_tolerance_j,
        }
    }
}

/// Original whole-window physical account and separately measured solver work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveSurfaceDiffusionReport {
    /// Only the accepted fine trajectory contributes physical storage/loss and
    /// accepted substeps/pair counts. Coarse/rejected losses are never deposited.
    pub diffusion: SurfaceDiffusionReport,
    pub trials: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub evaluated_pair_exchanges: usize,
    pub maximum_accepted_normalized_error: f64,
    /// Bounded proposal reusable as the initial step for another held interval.
    pub suggested_next_step_s: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum AdaptiveSurfaceDiffusionError {
    Diffusion(SurfaceDiffusionError),
    TrialBudget { completed_trials: usize },
    PairExchangeBudget { completed: usize, required: usize },
    TemporalAccuracy { step_s: f64, normalized_error: f64 },
    TimeResolution { elapsed_s: f64 },
}
impl From<SurfaceDiffusionError> for AdaptiveSurfaceDiffusionError {
    fn from(error: SurfaceDiffusionError) -> Self {
        Self::Diffusion(error)
    }
}
impl core::fmt::Display for AdaptiveSurfaceDiffusionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive surface diffusion: {self:?}")
    }
}
impl std::error::Error for AdaptiveSurfaceDiffusionError {}

/// Advance a homogeneous coating with mixed concentration-error control.
#[allow(clippy::too_many_arguments)]
pub fn diffuse_surface_surfactant_adaptive(
    law: SurfactantLaw,
    areas_m2: &[f64],
    amounts_mol: &mut [f64],
    edges: &[SurfaceDiffusionEdge],
    duration_s: f64,
    controls: AdaptiveSurfaceDiffusionControls,
    cancelled: impl FnMut() -> bool,
) -> Result<AdaptiveSurfaceDiffusionReport, AdaptiveSurfaceDiffusionError> {
    advance(
        Material::Uniform(law, edges),
        areas_m2,
        amounts_mol,
        duration_s,
        controls,
        cancelled,
    )
}

/// Advance heterogeneous original laws with adaptive conservative partitioning.
///
/// Each trial uses the same law-family, held-temperature and chemical-datum
/// requirements as the original partitioned operator. The infinity norm of
/// the estimated fine-path concentration error must be at most one. No amount
/// is extrapolated, clipped or published before the entire duration succeeds.
#[allow(clippy::too_many_arguments)]
pub fn diffuse_partitioned_surface_surfactant_adaptive(
    laws: &[SurfactantLaw],
    areas_m2: &[f64],
    amounts_mol: &mut [f64],
    edges: &[SurfacePartitionDiffusionEdge],
    duration_s: f64,
    controls: AdaptiveSurfaceDiffusionControls,
    cancelled: impl FnMut() -> bool,
) -> Result<AdaptiveSurfaceDiffusionReport, AdaptiveSurfaceDiffusionError> {
    if laws.len() != areas_m2.len() {
        return Err(SurfaceDiffusionError::InvalidInput(
            "one original constitutive law required per face",
        )
        .into());
    }
    advance(
        Material::Partitioned(laws, edges),
        areas_m2,
        amounts_mol,
        duration_s,
        controls,
        cancelled,
    )
}

fn advance(
    material: Material<'_>,
    areas: &[f64],
    amounts: &mut [f64],
    duration_s: f64,
    controls: AdaptiveSurfaceDiffusionControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<AdaptiveSurfaceDiffusionReport, AdaptiveSurfaceDiffusionError> {
    poll(&mut cancelled)?;
    if areas.is_empty() || areas.len() != amounts.len() {
        return Err(SurfaceDiffusionError::InvalidInput(
            "nonempty matching areas and amounts required",
        )
        .into());
    }
    if !duration_s.is_finite()
        || duration_s < 0.0
        || [
            controls.minimum_step_s,
            controls.initial_step_s,
            controls.maximum_step_s,
            controls.absolute_concentration_tolerance_mol_m2,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        || controls.minimum_step_s > controls.initial_step_s
        || controls.initial_step_s > controls.maximum_step_s
        || !controls.relative_tolerance.is_finite()
        || !(0.0..1.0).contains(&controls.relative_tolerance)
        || controls.maximum_trials == 0
    {
        return Err(SurfaceDiffusionError::InvalidInput(
            "finite duration and positive ordered adaptive time/accuracy/work controls required",
        )
        .into());
    }
    if areas.len() > controls.maximum_faces || material.edges() > controls.maximum_edges {
        return Err(SurfaceDiffusionError::WorkBudget.into());
    }
    let mut current = copied(amounts, &mut cancelled)?;
    // Zero-duration admission retains every original source/geometry/edge gate.
    let initial = diffuse(
        material,
        areas,
        &mut current,
        0.0,
        controls.step(controls.maximum_step_s),
        &mut cancelled,
    )?;
    let mut report = AdaptiveSurfaceDiffusionReport {
        diffusion: SurfaceDiffusionReport {
            duration_s,
            ..initial
        },
        trials: 0,
        accepted_intervals: 0,
        rejected_intervals: 0,
        evaluated_pair_exchanges: 0,
        maximum_accepted_normalized_error: 0.0,
        suggested_next_step_s: controls.initial_step_s,
    };
    if duration_s == 0.0 || material.edges() == 0 {
        poll(&mut cancelled)?;
        return Ok(report);
    }
    let sweep_work = material
        .edges()
        .checked_mul(2)
        .and_then(|n| n.checked_sub(1))
        .ok_or(SurfaceDiffusionError::WorkBudget)?;
    let trial_work = sweep_work
        .checked_mul(3)
        .ok_or(SurfaceDiffusionError::WorkBudget)?;
    let mut coarse = copied(&current, &mut cancelled)?;
    let mut fine = copied(&current, &mut cancelled)?;
    let mut elapsed = 0.0;
    while elapsed < duration_s {
        poll(&mut cancelled)?;
        if report.trials >= controls.maximum_trials {
            return Err(AdaptiveSurfaceDiffusionError::TrialBudget {
                completed_trials: report.trials,
            });
        }
        if trial_work
            > controls
                .maximum_pair_exchanges
                .saturating_sub(report.evaluated_pair_exchanges)
        {
            return Err(AdaptiveSurfaceDiffusionError::PairExchangeBudget {
                completed: report.evaluated_pair_exchanges,
                required: trial_work,
            });
        }
        let remaining = duration_s - elapsed;
        let end = if report.suggested_next_step_s >= remaining {
            duration_s
        } else {
            elapsed + report.suggested_next_step_s
        };
        let span = end - elapsed;
        let first_half = 0.5 * span;
        let second_half = span - first_half;
        if end <= elapsed || end > duration_s || first_half <= 0.0 || second_half <= 0.0 {
            return Err(AdaptiveSurfaceDiffusionError::TimeResolution { elapsed_s: elapsed });
        }
        copy_into(&current, &mut coarse, &mut cancelled)?;
        copy_into(&current, &mut fine, &mut cancelled)?;
        report.trials += 1;
        report.evaluated_pair_exchanges += trial_work;
        diffuse(
            material,
            areas,
            &mut coarse,
            span,
            controls.step(span),
            &mut cancelled,
        )?;
        for half in [first_half, second_half] {
            diffuse(
                material,
                areas,
                &mut fine,
                half,
                controls.step(half),
                &mut cancelled,
            )?;
        }
        // (a³+b³)/(3ab(a+b)), scaled without cubing dimensional seconds.
        let ratio = first_half.min(second_half) / first_half.max(second_half);
        let correction = (1.0 - ratio + ratio * ratio) / (3.0 * ratio);
        let error = error_norm(
            &current,
            &coarse,
            &fine,
            areas,
            correction,
            controls,
            &mut cancelled,
        )?;
        if error > 1.0 {
            report.rejected_intervals += 1;
            let smaller = (span * factor(error).min(0.8)).max(controls.minimum_step_s);
            if span <= controls.minimum_step_s
                || smaller >= span
                || smaller >= report.suggested_next_step_s
            {
                return Err(AdaptiveSurfaceDiffusionError::TemporalAccuracy {
                    step_s: span,
                    normalized_error: error,
                });
            }
            report.suggested_next_step_s = smaller;
            continue;
        }
        std::mem::swap(&mut current, &mut fine);
        report.accepted_intervals += 1;
        report.maximum_accepted_normalized_error =
            report.maximum_accepted_normalized_error.max(error);
        elapsed = end;
        // Clamp before multiplication so a finite maximum cannot overflow.
        let growth = factor(error);
        report.suggested_next_step_s = if span >= controls.maximum_step_s / growth {
            controls.maximum_step_s
        } else {
            (span * growth).max(controls.minimum_step_s)
        };
    }
    report.diffusion = account(
        material,
        areas,
        amounts,
        &current,
        duration_s,
        controls,
        &mut cancelled,
    )?;
    report.diffusion.substeps = report
        .accepted_intervals
        .checked_mul(2)
        .ok_or(SurfaceDiffusionError::WorkBudget)?;
    report.diffusion.pair_exchanges = report
        .diffusion
        .substeps
        .checked_mul(sweep_work)
        .ok_or(SurfaceDiffusionError::WorkBudget)?;
    poll(&mut cancelled)?;
    amounts.copy_from_slice(&current);
    Ok(report)
}

fn copied(
    source: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, SurfaceDiffusionError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(source.len())
        .map_err(|_| SurfaceDiffusionError::Allocation)?;
    for (i, &value) in source.iter().enumerate() {
        poll_at(cancelled, i)?;
        values.push(value);
    }
    Ok(values)
}

fn copy_into(
    source: &[f64],
    target: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), SurfaceDiffusionError> {
    for (i, (&value, out)) in source.iter().zip(target).enumerate() {
        poll_at(cancelled, i)?;
        *out = value;
    }
    Ok(())
}

fn error_norm(
    initial: &[f64],
    coarse: &[f64],
    fine: &[f64],
    areas: &[f64],
    correction: f64,
    controls: AdaptiveSurfaceDiffusionControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, SurfaceDiffusionError> {
    let mut error = 0.0_f64;
    for i in 0..initial.len() {
        poll_at(cancelled, i)?;
        // Original constitutive admission already guarantees finite physical
        // concentrations. Scale the tolerance sum, not the inventories, to
        // retain small errors after a large initial concentration has decayed.
        let difference = (fine[i] / areas[i] - coarse[i] / areas[i]).abs();
        let relative =
            controls.relative_tolerance * (initial[i] / areas[i]).max(fine[i] / areas[i]);
        let absolute = controls.absolute_concentration_tolerance_mol_m2;
        let scale = absolute.max(relative);
        error =
            error.max((difference / scale) * correction / (absolute / scale + relative / scale));
    }
    Ok(error)
}

fn factor(error: f64) -> f64 {
    if error == 0.0 {
        2.0
    } else if error.is_infinite() {
        0.2
    } else {
        (0.9 * det::pow(error, -1.0 / 3.0)).clamp(0.2, 2.0)
    }
}

#[allow(clippy::too_many_arguments)]
fn account(
    material: Material<'_>,
    areas: &[f64],
    before: &[f64],
    after: &[f64],
    duration_s: f64,
    controls: AdaptiveSurfaceDiffusionControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SurfaceDiffusionReport, SurfaceDiffusionError> {
    let (mut amount_before, mut amount_after, mut defect) =
        (Sum::default(), Sum::default(), Sum::default());
    let (mut energy_before, mut energy_after, mut change) =
        (Sum::default(), Sum::default(), Sum::default());
    for i in 0..before.len() {
        poll_at(cancelled, i)?;
        let initial = SurfactantFace {
            law: material.law(i),
            amount_mol: before[i],
        };
        let final_face = SurfactantFace {
            law: material.law(i),
            amount_mol: after[i],
        };
        let failure = |error| SurfaceDiffusionError::Constitutive { face: i, error };
        energy_before.add(initial.at_area(areas[i]).map_err(failure)?.free_energy_j);
        energy_after.add(final_face.at_area(areas[i]).map_err(failure)?.free_energy_j);
        change.add(
            initial
                .energy_change_at_fixed_area(areas[i], after[i])
                .map_err(failure)?,
        );
        amount_before.add(before[i]);
        amount_after.add(after[i]);
        defect.add(after[i] - before[i]);
    }
    let defect = finite(defect.value())?;
    let change = finite(change.value())?;
    if defect.abs() > controls.amount_tolerance_mol {
        return Err(SurfaceDiffusionError::AmountBalance { defect_mol: defect });
    }
    if change > controls.energy_increase_tolerance_j {
        return Err(SurfaceDiffusionError::FreeEnergyIncrease { change_j: change });
    }
    Ok(SurfaceDiffusionReport {
        duration_s,
        substeps: 0,
        pair_exchanges: 0,
        amount_before_mol: finite(amount_before.value())?,
        amount_after_mol: finite(amount_after.value())?,
        amount_defect_mol: defect,
        free_energy_before_j: finite(energy_before.value())?,
        free_energy_after_j: finite(energy_after.value())?,
        free_energy_change_j: change,
        free_energy_dissipation_j: -change,
    })
}

#[cfg(test)]
mod tests;
