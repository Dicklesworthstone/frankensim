//! Conjugate cooling of multiple independently supplied airflow branches.
//!
//! Every iteration runs ONE solid solve with all Robin reference temperatures.
//! Each branch then marches its own inlet, capacity rate, and ordered segments.
//! Branches exchange heat through the solid, not through an invented air mixer.
//! The hydraulic solution is frozen: heating does not redistribute flow.
//!
//! The existing single-path driver remains the admission authority for each
//! response, temperature criterion, and watt-balance gate. A one-step probe
//! reuses those checks without rerunning the solid. In particular a high-power
//! branch cannot dilute a smaller branch's balance threshold. This composes
//! the existing nominal model; it does not certify uncertainty or passivity.

use std::collections::BTreeSet;

use fs_couple::AitkenRelaxation;
use fs_exec::Cx;

use crate::AirflowError;
use crate::conjugate::{
    AirPath, ConjugateConfig, ConjugateIteration, ConjugateSolution, Relaxation,
    SolidRegionState, solve_conjugate_from,
};

/// A common-solid fixed point, retaining each branch's independent air state.
#[derive(Debug, Clone, PartialEq)]
pub struct ConjugateBranchesSolution {
    /// Branch-major, then stream-wise references, in the caller's path order.
    pub reference_temperatures_k: Vec<f64>,
    /// Branch solutions from the SAME final solid response. Each history spans
    /// all shared iterations, including iterations when that branch alone met
    /// its stopping criterion but another branch did not.
    pub branches: Vec<ConjugateSolution>,
    /// Number of shared solid solves, not the sum of per-branch probes.
    pub iterations: usize,
}

/// Solve independently supplied air paths coupled through one solid response.
///
/// The callback receives all references in branch-major, stream-wise order
/// and must return all regions in that same order. Each path needs a distinct
/// set of solid Robin regions; one surface cannot have two reference rows.
/// Inputs and iteration limits are admitted before the first solid solve.
///
/// # Errors
/// Returns the existing typed single-path refusals, plus `EmptyAirPath` for
/// no branches and `DuplicateAirSegment` for a region shared between paths.
/// Cancellation retains the COMPLETE reference vector, never just one branch.
pub fn solve_conjugate_branches<F>(
    cx: &Cx<'_>,
    paths: &[AirPath],
    config: &ConjugateConfig,
    solid: F,
) -> Result<ConjugateBranchesSolution, AirflowError>
where
    F: FnMut(&Cx<'_>, &[f64]) -> Result<Vec<SolidRegionState>, AirflowError>,
{
    let initial: Vec<f64> = paths
        .iter()
        .flat_map(|path| vec![path.inlet_temperature_k(); path.segments().len()])
        .collect();
    solve_conjugate_branches_from(cx, paths, config, &initial, solid)
}

/// Resume from a complete branch-major reference vector.
///
/// Fixed relaxation reproduces the uninterrupted numerical tail for a
/// deterministic solid callback. As in the single-path driver, Aitken history
/// restarts: its resumed tail is a continuation, not a bitwise replay.
/// Scalar Aitken relaxation is independent per branch; it is not IQN-ILS.
///
/// # Errors
/// The same refusals as [`solve_conjugate_branches`], including arity and
/// finiteness checks on the resumed vector before invoking the callback.
pub fn solve_conjugate_branches_from<F>(
    cx: &Cx<'_>,
    paths: &[AirPath],
    config: &ConjugateConfig,
    initial_references_k: &[f64],
    mut solid: F,
) -> Result<ConjugateBranchesSolution, AirflowError>
where
    F: FnMut(&Cx<'_>, &[f64]) -> Result<Vec<SolidRegionState>, AirflowError>,
{
    let offsets = admit(paths, config, initial_references_k)?;
    if paths.len() == 1 {
        let solution = solve_conjugate_from(cx, &paths[0], config, initial_references_k, solid)?;
        return Ok(ConjugateBranchesSolution {
            reference_temperatures_k: solution.reference_temperatures_k.clone(),
            iterations: solution.iterations,
            branches: vec![solution],
        });
    }
    let probe_config = ConjugateConfig {
        max_iterations: 1,
        relaxation: Relaxation::Fixed { omega: 1.0 },
        ..*config
    };
    let mut relaxers: Vec<Option<AitkenRelaxation>> = paths
        .iter()
        .map(|_| match config.relaxation {
            Relaxation::Fixed { .. } => None,
            Relaxation::Aitken { omega_init, omega_max } => {
                Some(AitkenRelaxation::new(omega_init, omega_max))
            }
        })
        .collect();
    let mut histories: Vec<Vec<ConjugateIteration>> = vec![Vec::new(); paths.len()];
    let mut worst = vec![0.0_f64; paths.len()];
    let mut reference = initial_references_k.to_vec();
    let mut last_change = 0.0_f64;
    for iteration in 0..config.max_iterations {
        checkpoint(cx, iteration, &reference)?;
        let response = solid(cx, &reference)?;
        checkpoint(cx, iteration, &reference)?;
        if response.len() != reference.len() {
            return Err(AirflowError::SolidResponseArity {
                expected: reference.len(),
                found: response.len(),
            });
        }
        let mut admitted = Vec::with_capacity(paths.len());
        let mut updated = Vec::with_capacity(reference.len());
        let mut omegas = Vec::with_capacity(paths.len());
        last_change = 0.0;
        for (index, path) in paths.iter().enumerate() {
            checkpoint(cx, iteration, &reference)?;
            let range = offsets[index]..offsets[index + 1];
            let current = &reference[range.clone()];
            let states = &response[range];
            // This probe executes no solid work. It reuses every production
            // wiring and per-path balance check instead of forking their rules.
            let accepted = match solve_conjugate_from(
                cx, path, &probe_config, current, |_, _| Ok(states.to_vec()),
            ) {
                Ok(solution) => Some(solution),
                Err(AirflowError::ConjugateNotConverged { .. }) => None,
                Err(AirflowError::Cancelled { .. }) => {
                    return Err(AirflowError::Cancelled {
                        iteration,
                        references_k: reference,
                    });
                }
                Err(error) => return Err(error),
            };
            let march = match accepted.as_ref() {
                Some(solution) => solution.march.clone(),
                None => path.march(
                    &states.iter().map(|state| state.mean_wall_temperature_k).collect::<Vec<_>>(),
                )?,
            };
            let next = march.reference_temperatures_k();
            let total_area = checked_sum("branch wetted area", path.segments().iter().map(|s| s.area_m2()))?;
            let mut max_change = 0.0_f64;
            let mut scalar = 0.0;
            for ((&old, &new), segment) in current.iter().zip(&next).zip(path.segments()) {
                let delta = checked("branch reference residual", new - old)?;
                max_change = max_change.max(delta.abs());
                scalar = checked("branch weighted reference residual", scalar + delta * (segment.area_m2() / total_area))?;
            }
            last_change = last_change.max(max_change);
            let omega = checked("branch relaxation omega", match config.relaxation {
                Relaxation::Fixed { omega } => omega,
                Relaxation::Aitken { omega_init, .. } => relaxers[index]
                    .as_mut().map_or(omega_init, |relaxer| relaxer.next_omega(scalar)),
            })?;
            let solid_total = checked_sum("branch solid heat rate", states.iter().map(|s| s.heat_rate_w))?;
            let imbalance = checked("branch interface heat rate", solid_total - march.total_heat_rate_w)?;
            worst[index] = worst[index].max(imbalance.abs());
            histories[index].push(ConjugateIteration {
                iteration,
                max_reference_change_k: max_change,
                scalar_residual_k: scalar,
                relaxation_omega: omega,
                air_outlet_temperature_k: march.outlet_temperature_k,
                solid_heat_rate_w: solid_total,
                air_heat_rate_w: march.total_heat_rate_w,
                interface_imbalance_w: imbalance,
                reference_temperatures_k: current.to_vec(),
            });
            updated.extend(next);
            omegas.push(omega);
            admitted.push(accepted);
        }
        checkpoint(cx, iteration, &reference)?;
        if admitted.iter().all(Option::is_some) {
            let mut branches = Vec::with_capacity(paths.len());
            for (index, accepted) in admitted.into_iter().enumerate() {
                let mut solution = accepted.expect("all branch probes admitted");
                solution.iterations = iteration + 1;
                solution.history = std::mem::take(&mut histories[index]);
                solution.worst_recorded_imbalance_w = worst[index];
                branches.push(solution);
            }
            return Ok(ConjugateBranchesSolution {
                reference_temperatures_k: updated,
                branches,
                iterations: iteration + 1,
            });
        }
        // Even an individually converged branch moves again: another branch
        // can change its wall temperature through the shared solid next time.
        for (index, &omega) in omegas.iter().enumerate() {
            for slot in offsets[index]..offsets[index + 1] {
                reference[slot] = checked("relaxed branch reference temperature",
                    reference[slot] + omega * (updated[slot] - reference[slot]))?;
            }
        }
    }
    Err(AirflowError::ConjugateNotConverged {
        iterations: config.max_iterations,
        max_change_bits: last_change.to_bits(),
        tolerance_bits: config.temperature_tolerance_k.to_bits(),
    })
}

fn admit(paths: &[AirPath], config: &ConjugateConfig, initial: &[f64]) -> Result<Vec<usize>, AirflowError> {
    if paths.is_empty() {
        return Err(AirflowError::EmptyAirPath);
    }
    for (field, value) in [
        ("conjugate temperature tolerance", config.temperature_tolerance_k),
        ("conjugate balance tolerance", config.balance_tolerance_w),
    ] {
        if !(value.is_finite() && value > 0.0) {
            return Err(AirflowError::InvalidConjugateInput { field, value_bits: value.to_bits() });
        }
    }
    if config.max_iterations == 0 {
        return Err(AirflowError::InvalidConjugateInput { field: "conjugate max iterations", value_bits: 0 });
    }
    if !(0.0..1.0).contains(&config.balance_relative_tolerance) {
        return Err(AirflowError::InvalidConjugateInput {
            field: "conjugate relative balance tolerance",
            value_bits: config.balance_relative_tolerance.to_bits(),
        });
    }
    let (omega, cap) = match config.relaxation {
        Relaxation::Fixed { omega } => (omega, 2.0),
        Relaxation::Aitken { omega_init, omega_max } => (omega_init, omega_max),
    };
    if !(omega.is_finite() && cap.is_finite() && omega > 0.0 && omega <= cap) {
        return Err(AirflowError::InvalidConjugateInput {
            field: "conjugate relaxation omega", value_bits: omega.to_bits(),
        });
    }
    let mut names = BTreeSet::new();
    let mut offsets = vec![0_usize];
    for path in paths {
        for segment in path.segments() {
            if !names.insert(segment.region()) {
                return Err(AirflowError::DuplicateAirSegment { region: segment.region().to_string() });
            }
        }
        let end = offsets.last().copied().unwrap_or(0).checked_add(path.segments().len())
            .ok_or(AirflowError::InvalidConjugateInput { field: "branch segment count", value_bits: u64::MAX })?;
        offsets.push(end);
    }
    let expected = offsets.last().copied().unwrap_or(0);
    if expected != initial.len() {
        return Err(AirflowError::SolidResponseArity { expected, found: initial.len() });
    }
    for &value in initial {
        if !value.is_finite() {
            return Err(AirflowError::InvalidConjugateInput {
                field: "resumed reference temperature", value_bits: value.to_bits(),
            });
        }
    }
    Ok(offsets)
}

fn checkpoint(cx: &Cx<'_>, iteration: usize, references: &[f64]) -> Result<(), AirflowError> {
    cx.checkpoint().map_err(|_| AirflowError::Cancelled {
        iteration, references_k: references.to_vec(),
    })
}

fn checked(stage: &'static str, value: f64) -> Result<f64, AirflowError> {
    if value.is_finite() { Ok(value) }
    else { Err(AirflowError::NonFiniteCoupling { stage, value_bits: value.to_bits() }) }
}

fn checked_sum(stage: &'static str, values: impl IntoIterator<Item = f64>) -> Result<f64, AirflowError> {
    values.into_iter().try_fold(0.0, |sum, value| checked(stage, sum + value))
}

#[cfg(test)]
mod tests;
