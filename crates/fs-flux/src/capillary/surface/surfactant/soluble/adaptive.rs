//! Local time-error control for the original finite-bath soluble surface map.
//!
//! Each comparison uses one original backward-Euler exchange step versus two
//! half steps, with the original lateral sweeps in each map. Its fine-path
//! estimate is FIRST order: the equal-half discrepancy has factor one, not the
//! one-third factor of the second-order surface-only diffusion composition.
//! No extrapolation, frozen bath or extra chemical reference is introduced.
//!
//! Areas, temperature, solvent volume, topology and coefficients are held.
//! Accuracy acceptance checks every actual face amount, the selected finite
//! bulk amount, and that species' activity from the original nonideal model.
//! Relative scales measure changes since the last accepted state, so a large
//! background inventory does not hide a small transfer. These are local error
//! estimates, not a global or spatial accuracy guarantee.

use super::super::time_control::{self, Limits, Order, Schedule};
use super::{
    AqueousSolution, MAX_OSMOTIC_SOLUTES, OsmosisError, SolubleSurfaceControls,
    SolubleSurfaceError, SolubleSurfaceKinetics, SolubleSurfaceReport, Sum,
    SurfacePartitionDiffusionEdge, SurfactantLaw, admit_energy, advance_soluble_surface, copied,
    det, energy_change, finite, poll, poll_at,
};

/// Time-error limits in addition to the original physical/work controls.
/// The original `maximum_step_s` caps each full comparison span. Every full
/// or half proposal executes exactly one original soluble substep.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveSolubleSurfaceControls {
    pub initial_step_s: f64,
    /// A shorter final remainder may be tried, but must satisfy the same gates.
    pub minimum_step_s: f64,
    pub absolute_surface_amount_tolerance_mol: f64,
    pub absolute_bulk_amount_tolerance_mol: f64,
    /// Dimensionless selected-solute activity on the original source datum.
    /// Rescale this allowance if the source's activity reference is changed.
    pub absolute_bulk_activity_tolerance: f64,
    /// In [0,1), relative to max(|fine-start|, |coarse-start|).
    pub relative_tolerance: f64,
    /// Full-versus-two-half comparisons, including rejected comparisons.
    pub maximum_trials: usize,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveSolubleSurfaceReport {
    /// Only accepted fine steps supply transfers, dissipation, substeps and
    /// root/pair counts. The final amount/energy account is reconstructed
    /// separately from the original inventories and the final candidate.
    pub transport: SolubleSurfaceReport,
    pub trials: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    /// All completed original substeps, including coarse and rejected work.
    pub evaluated_substeps: usize,
    pub evaluated_root_evaluations: usize,
    pub evaluated_surface_pair_exchanges: usize,
    pub maximum_accepted_normalized_error: f64,
    pub suggested_next_step_s: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum AdaptiveSolubleSurfaceError {
    /// A root-budget error reports the whole-call completed evaluation count.
    Transport(SolubleSurfaceError),
    TrialBudget {
        completed_trials: usize,
    },
    SubstepBudget {
        completed: usize,
        required: usize,
    },
    PairExchangeBudget {
        completed: usize,
        required: usize,
    },
    TemporalAccuracy {
        step_s: f64,
        normalized_error: f64,
    },
    TimeResolution {
        elapsed_s: f64,
    },
}
impl From<SolubleSurfaceError> for AdaptiveSolubleSurfaceError {
    fn from(error: SolubleSurfaceError) -> Self {
        Self::Transport(error)
    }
}
impl core::fmt::Display for AdaptiveSolubleSurfaceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive soluble surface transport: {self:?}")
    }
}
impl std::error::Error for AdaptiveSolubleSurfaceError {}
type Result<T> = core::result::Result<T, AdaptiveSolubleSurfaceError>;

#[derive(Default)]
struct Work {
    substeps: usize,
    roots: usize,
    pairs: usize,
}

/// Complete a held finite-bath duration and publish both original inventories
/// only after every interval and the independent whole-duration account pass.
///
/// Original substep, root-evaluation and pair-exchange caps are shared by the
/// full, both half, and every rejected trajectory. Three substeps and all their
/// pair applications are reserved before a comparison. Root work is charged
/// as the original monotone solver consumes it. Physical/source refusals and
/// cancellation propagate; only excessive temporal error shortens a proposal.
/// The original packing, convex-source, chemical-datum and fixed-bath domain
/// of [`advance_soluble_surface`] is retained, including nonideal co-solutes.
#[allow(clippy::too_many_arguments)]
pub fn advance_soluble_surface_adaptive(
    bulk: &mut AqueousSolution,
    laws: &[SurfactantLaw],
    areas_m2: &[f64],
    amounts_mol: &mut [f64],
    kinetics: SolubleSurfaceKinetics<'_>,
    edges: &[SurfacePartitionDiffusionEdge],
    duration_s: f64,
    physical: SolubleSurfaceControls,
    accuracy: AdaptiveSolubleSurfaceControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<AdaptiveSolubleSurfaceReport> {
    poll(&mut cancelled)?;
    let mut schedule = Schedule::new(
        duration_s,
        Limits {
            initial_step_s: accuracy.initial_step_s,
            minimum_step_s: accuracy.minimum_step_s,
            maximum_step_s: physical.maximum_step_s,
            maximum_trials: accuracy.maximum_trials,
        },
        Order::First,
    )
    .map_err(schedule_error)?;
    if [
        accuracy.absolute_surface_amount_tolerance_mol,
        accuracy.absolute_bulk_amount_tolerance_mol,
        accuracy.absolute_bulk_activity_tolerance,
    ]
    .iter()
    .any(|value| !value.is_finite() || *value <= 0.0)
        || !accuracy.relative_tolerance.is_finite()
        || !(0.0..1.0).contains(&accuracy.relative_tolerance)
    {
        return Err(SolubleSurfaceError::Invalid(
            "positive finite amount/activity tolerances and relative tolerance in [0,1) required",
        )
        .into());
    }
    if areas_m2.is_empty()
        || laws.len() != areas_m2.len()
        || amounts_mol.len() != areas_m2.len()
        || kinetics.desorption_rates_s_inverse.len() != areas_m2.len()
    {
        return Err(
            SolubleSurfaceError::Invalid("matching nonempty surface layouts required").into(),
        );
    }
    if areas_m2.len() > physical.maximum_faces || edges.len() > physical.maximum_edges {
        return Err(SolubleSurfaceError::WorkBudget.into());
    }
    let mut current = copied(amounts_mol, &mut cancelled)?;
    let mut current_bulk = bulk.clone();
    // Keep every original layout, thermodynamic, packing and graph admission,
    // even when the requested duration is zero.
    let mut transport = advance_soluble_surface(
        &mut current_bulk,
        laws,
        areas_m2,
        &mut current,
        kinetics,
        edges,
        0.0,
        physical,
        &mut cancelled,
    )?;
    transport.duration_s = duration_s;
    let mut work = Work::default();
    let mut adsorbed = Sum::default();
    let mut desorbed = Sum::default();
    let mut exchange_loss = Sum::default();
    let mut diffusion_loss = Sum::default();
    if !schedule.complete() {
        let pairs_per_step = if edges.is_empty() {
            0
        } else {
            edges
                .len()
                .checked_mul(2)
                .and_then(|value| value.checked_sub(1))
                .and_then(|value| value.checked_mul(2))
                .ok_or(SolubleSurfaceError::WorkBudget)?
        };
        let pairs_per_trial = pairs_per_step
            .checked_mul(3)
            .ok_or(SolubleSurfaceError::WorkBudget)?;
        let mut coarse = copied(&current, &mut cancelled)?;
        let mut fine = copied(&current, &mut cancelled)?;
        let mut coarse_bulk = bulk.clone();
        let mut fine_bulk = bulk.clone();
        let mut current_activity = log_activity(&current_bulk, kinetics.bulk_species)?;
        while !schedule.complete() {
            poll(&mut cancelled)?;
            schedule.check_trial_budget().map_err(schedule_error)?;
            if 3 > physical.maximum_substeps.saturating_sub(work.substeps) {
                return Err(AdaptiveSolubleSurfaceError::SubstepBudget {
                    completed: work.substeps,
                    required: 3,
                });
            }
            if pairs_per_trial
                > physical
                    .maximum_surface_pair_exchanges
                    .saturating_sub(work.pairs)
            {
                return Err(AdaptiveSolubleSurfaceError::PairExchangeBudget {
                    completed: work.pairs,
                    required: pairs_per_trial,
                });
            }
            let trial = schedule.begin().map_err(schedule_error)?;
            for (face, &amount) in current.iter().enumerate() {
                poll_at(&mut cancelled, face)?;
                coarse[face] = amount;
                fine[face] = amount;
            }
            for candidate in [&mut coarse_bulk, &mut fine_bulk] {
                candidate
                    .replace_inventories(
                        current_bulk.solvent_volume_m3(),
                        current_bulk.solute_moles(),
                    )
                    .map_err(SolubleSurfaceError::from)?;
            }
            step(
                &mut coarse_bulk,
                laws,
                areas_m2,
                &mut coarse,
                kinetics,
                edges,
                trial.span_s,
                physical,
                &mut work,
                &mut cancelled,
            )?;
            let first = step(
                &mut fine_bulk,
                laws,
                areas_m2,
                &mut fine,
                kinetics,
                edges,
                trial.halves_s[0],
                physical,
                &mut work,
                &mut cancelled,
            )?;
            let second = step(
                &mut fine_bulk,
                laws,
                areas_m2,
                &mut fine,
                kinetics,
                edges,
                trial.halves_s[1],
                physical,
                &mut work,
                &mut cancelled,
            )?;
            let coarse_activity = log_activity(&coarse_bulk, kinetics.bulk_species)?;
            let fine_activity = log_activity(&fine_bulk, kinetics.bulk_species)?;
            let mut error = 0.0_f64;
            for face in 0..current.len() {
                poll_at(&mut cancelled, face)?;
                error = error.max(amount_error(
                    current[face],
                    coarse[face],
                    fine[face],
                    accuracy.absolute_surface_amount_tolerance_mol,
                    accuracy.relative_tolerance,
                    trial.fine_error_factor,
                ));
            }
            let species = kinetics.bulk_species;
            error = error.max(amount_error(
                current_bulk.solute_moles()[species],
                coarse_bulk.solute_moles()[species],
                fine_bulk.solute_moles()[species],
                accuracy.absolute_bulk_amount_tolerance_mol,
                accuracy.relative_tolerance,
                trial.fine_error_factor,
            ));
            error = error.max(activity_error(
                current_activity,
                coarse_activity,
                fine_activity,
                accuracy.absolute_bulk_activity_tolerance,
                accuracy.relative_tolerance,
                trial.fine_error_factor,
            ));
            if !schedule.consider(trial, error).map_err(schedule_error)? {
                continue;
            }
            for accepted in [first, second] {
                transport.substeps += accepted.substeps;
                transport.root_evaluations += accepted.root_evaluations;
                transport.surface_pair_exchanges += accepted.surface_pair_exchanges;
                transport.maximum_root_residual_mol = transport
                    .maximum_root_residual_mol
                    .max(accepted.maximum_root_residual_mol);
                transport.maximum_site_kinetic_residual_mol = transport
                    .maximum_site_kinetic_residual_mol
                    .max(accepted.maximum_site_kinetic_residual_mol);
                adsorbed.add(accepted.gross_adsorption_mol);
                desorbed.add(accepted.gross_desorption_mol);
                exchange_loss.add(accepted.exchange_free_energy_dissipation_j);
                diffusion_loss.add(accepted.surface_diffusion_free_energy_dissipation_j);
            }
            core::mem::swap(&mut current, &mut fine);
            core::mem::swap(&mut current_bulk, &mut fine_bulk);
            current_activity = fine_activity;
        }
    }
    // Whole-duration gates use the original endpoints, not an accumulation of
    // each stage's permitted amount/energy drift.
    let species = kinetics.bulk_species;
    let mut defect = Sum::default();
    let mut surface_after = Sum::default();
    defect.add(current_bulk.solute_moles()[species] - bulk.solute_moles()[species]);
    for (face, (&before, &after)) in amounts_mol.iter().zip(&current).enumerate() {
        poll_at(&mut cancelled, face)?;
        defect.add(after - before);
        surface_after.add(after);
    }
    transport.amount_defect_mol = finite(defect.value())?;
    if transport.amount_defect_mol.abs() > physical.amount_tolerance_mol {
        return Err(SolubleSurfaceError::AmountBalance {
            defect_mol: transport.amount_defect_mol,
        }
        .into());
    }
    transport.free_energy_change_j = energy_change(
        bulk,
        species,
        bulk.solute_moles(),
        current_bulk.solute_moles(),
        laws,
        areas_m2,
        areas_m2,
        amounts_mol,
        &current,
        kinetics.standard_potential_offset_j_mol,
        &mut cancelled,
    )?;
    admit_energy(transport.free_energy_change_j, physical)?;
    transport.exchange_free_energy_dissipation_j = finite(exchange_loss.value())?;
    transport.surface_diffusion_free_energy_dissipation_j = finite(diffusion_loss.value())?;
    let mut account = Sum::default();
    account.add(transport.free_energy_change_j);
    account.add(transport.exchange_free_energy_dissipation_j);
    account.add(transport.surface_diffusion_free_energy_dissipation_j);
    transport.energy_account_residual_j = finite(account.value())?;
    if transport.energy_account_residual_j.abs() > physical.energy_tolerance_j {
        return Err(SolubleSurfaceError::EnergyAccount {
            residual_j: transport.energy_account_residual_j,
        }
        .into());
    }
    transport.bulk_amount_after_mol = current_bulk.solute_moles()[species];
    transport.surface_amount_after_mol = finite(surface_after.value())?;
    transport.gross_adsorption_mol = finite(adsorbed.value())?;
    transport.gross_desorption_mol = finite(desorbed.value())?;
    let statistics = schedule.statistics;
    let report = AdaptiveSolubleSurfaceReport {
        transport,
        trials: statistics.trials,
        accepted_intervals: statistics.accepted_intervals,
        rejected_intervals: statistics.rejected_intervals,
        evaluated_substeps: work.substeps,
        evaluated_root_evaluations: work.roots,
        evaluated_surface_pair_exchanges: work.pairs,
        maximum_accepted_normalized_error: statistics.maximum_accepted_normalized_error,
        suggested_next_step_s: statistics.suggested_next_step_s,
    };
    poll(&mut cancelled)?;
    *bulk = current_bulk;
    amounts_mol.copy_from_slice(&current);
    Ok(report)
}

#[allow(clippy::too_many_arguments)]
fn step(
    bulk: &mut AqueousSolution,
    laws: &[SurfactantLaw],
    areas: &[f64],
    amounts: &mut [f64],
    kinetics: SolubleSurfaceKinetics<'_>,
    edges: &[SurfacePartitionDiffusionEdge],
    duration_s: f64,
    physical: SolubleSurfaceControls,
    work: &mut Work,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SolubleSurfaceReport> {
    let roots_remaining = physical.maximum_root_evaluations.saturating_sub(work.roots);
    if roots_remaining == 0 {
        return Err(SolubleSurfaceError::RootBudget {
            completed_evaluations: work.roots,
        }
        .into());
    }
    let report = advance_soluble_surface(
        bulk,
        laws,
        areas,
        amounts,
        kinetics,
        edges,
        duration_s,
        SolubleSurfaceControls {
            maximum_step_s: duration_s,
            maximum_substeps: 1,
            maximum_root_evaluations: roots_remaining,
            maximum_surface_pair_exchanges: physical
                .maximum_surface_pair_exchanges
                .saturating_sub(work.pairs),
            ..physical
        },
        cancelled,
    )
    .map_err(|error| match error {
        SolubleSurfaceError::RootBudget {
            completed_evaluations,
        } => AdaptiveSolubleSurfaceError::Transport(SolubleSurfaceError::RootBudget {
            completed_evaluations: work.roots + completed_evaluations,
        }),
        other => other.into(),
    })?;
    work.substeps += report.substeps;
    work.roots += report.root_evaluations;
    work.pairs += report.surface_pair_exchanges;
    Ok(report)
}

fn amount_error(
    initial: f64,
    coarse: f64,
    fine: f64,
    absolute: f64,
    relative: f64,
    factor: f64,
) -> f64 {
    let relative = relative * (fine - initial).abs().max((coarse - initial).abs());
    let scale = absolute.max(relative);
    ((fine - coarse).abs() / scale) * factor / (absolute / scale + relative / scale)
}

// The exact same concentration/activity-coefficient law as the original root.
// A vacant bath is activity zero, represented by -infinity only in this
// comparison helper. Co-solute interactions retain their original inventories.
fn log_activity(bulk: &AqueousSolution, species: usize) -> super::Result<f64> {
    if bulk.solute_moles()[species] == 0.0 {
        return Ok(f64::NEG_INFINITY);
    }
    let model = bulk.nonideal_model().ok_or(SolubleSurfaceError::Invalid(
        "explicit bulk activity model required",
    ))?;
    let count = bulk.solute_moles().len();
    let mut concentration = [0.0; MAX_OSMOTIC_SOLUTES];
    let mut excess = [0.0; MAX_OSMOTIC_SOLUTES];
    bulk.concentrations_into(&mut concentration[..count])?;
    model
        .log_activity_coefficients_into(&concentration[..count], &mut excess[..count])
        .map_err(OsmosisError::from)?;
    finite(
        det::ln(concentration[species]) - det::ln(model.reference_concentration_mol_m3())
            + excess[species],
    )
}

// Normalize activities in log space, avoiding overflow/underflow from an
// authored reference while preserving small differences with expm1. No rate,
// reference, clipping or material limit is changed by the error estimator.
fn activity_error(
    initial: f64,
    coarse: f64,
    fine: f64,
    absolute: f64,
    relative: f64,
    factor: f64,
) -> f64 {
    let absolute_log = det::ln(absolute);
    let mut scale = coarse.max(fine).max(absolute_log);
    if relative > 0.0 {
        scale = scale.max(initial);
    }
    let difference = scaled_activity_difference(fine, coarse, scale);
    if difference == 0.0 {
        return 0.0;
    }
    let change = if relative == 0.0 {
        0.0
    } else {
        scaled_activity_difference(fine, initial, scale)
            .max(scaled_activity_difference(coarse, initial, scale))
    };
    difference * factor / (det::exp(absolute_log - scale) + relative * change)
}

fn scaled_activity_difference(a: f64, b: f64, scale: f64) -> f64 {
    if a == b {
        return 0.0;
    }
    let upper = a.max(b);
    let lower = a.min(b);
    let fraction = if lower == f64::NEG_INFINITY {
        1.0
    } else {
        -det::expm1(lower - upper)
    };
    det::exp(upper - scale) * fraction
}

fn schedule_error(error: time_control::Error) -> AdaptiveSolubleSurfaceError {
    match error {
        time_control::Error::Invalid => SolubleSurfaceError::Invalid(
            "finite duration and positive ordered adaptive time/work controls required",
        )
        .into(),
        time_control::Error::TrialBudget { completed_trials } => {
            AdaptiveSolubleSurfaceError::TrialBudget { completed_trials }
        }
        time_control::Error::TemporalAccuracy {
            step_s,
            normalized_error,
        } => AdaptiveSolubleSurfaceError::TemporalAccuracy {
            step_s,
            normalized_error,
        },
        time_control::Error::TimeResolution { elapsed_s } => {
            AdaptiveSolubleSurfaceError::TimeResolution { elapsed_s }
        }
    }
}

#[cfg(test)]
mod tests;
