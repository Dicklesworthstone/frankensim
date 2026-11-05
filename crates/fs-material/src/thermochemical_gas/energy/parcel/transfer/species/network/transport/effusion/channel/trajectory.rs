//! Adaptive finite-reservoir molecular channels with explicitly held walls.
//!
//! Every trial compares one original full step with two resampled half steps;
//! only the half steps enter the accepted trajectory. Error control covers gas
//! temperature, pressure, species inventories, directional transport, reflected
//! withdrawals and each wall's heat. Original material/source refusals remain
//! hard errors; only typed time-step refusals can retry.
//!
//! This is a local time-discretization estimate, not a global accuracy certificate.
//! Geometry, original sources and thermostat temperatures remain fixed over the
//! window. A separate original-to-final audit admits gas plus wall energy and
//! entropy with ONE physical budget, independent of the accepted step count.

mod accounts;
pub use super::super::trajectory::{EffusionAdvanceControls, EffusionAdvanceError};
pub use accounts::{DiffuseChannelIntegratedWall, DiffuseChannelWindowBalance};

use super::super::trajectory::{MAX_EFFUSION_RETAINED_ROWS, retryable_native, state_error};
use super::{
    ClosedGasParcel, DiffuseChannelNetworkCandidate, DiffuseChannelNetworkLink,
    DiffuseChannelNetworkStep, GasEnergySolve, GasTransportError, GasTransportPaths,
    GasTransportSolve, MAX_SNAPSHOT_PARCELS, step_candidate,
};
use EffusionAdvanceError as Error;
use std::collections::BTreeMap;

#[derive(Clone, Debug, PartialEq)]
pub struct DiffuseChannelAdvanceStep {
    pub duration_s: f64,
    pub attempted_steps: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub maximum_accepted_error: f64,
    /// Chronological original half steps; rejected/coarse probes are excluded.
    pub accepted_steps: Vec<DiffuseChannelNetworkStep>,
    /// Canonical per-channel thermostat accounts over the complete window.
    pub walls: Vec<DiffuseChannelIntegratedWall>,
    pub balance: DiffuseChannelWindowBalance,
}

/// Complete immutable proposal; a refused window publishes no partial gas state.
#[derive(Clone, Debug, PartialEq)]
pub struct DiffuseChannelAdvanceCandidate {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    step: DiffuseChannelAdvanceStep,
}
impl DiffuseChannelAdvanceCandidate {
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub const fn step(&self) -> &DiffuseChannelAdvanceStep {
        &self.step
    }
    pub fn into_parts(self) -> (BTreeMap<u64, ClosedGasParcel>, DiffuseChannelAdvanceStep) {
        (self.parcels, self.step)
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn physical(error: GasTransportError, elapsed_s: f64, duration_s: f64) -> Error {
    if error == GasTransportError::Cancelled {
        Error::Cancelled
    } else {
        Error::Step {
            elapsed_s,
            duration_s,
            error,
        }
    }
}
fn retryable(error: &Error) -> bool {
    matches!(error, Error::Step { error, .. } if retryable_native(error))
}

/// Advance channels and all original internal transport paths together.
///
/// Native gas reconstruction, gross withdrawals, wall accommodation, collision
/// qualification and endpoint source admission remain owned by step_candidate.
/// The shared effusion controls preserve dimensional tolerances and finite work.
/// Stream energy tolerance also applies separately to every wall heat account;
/// cancellation between hot and cold walls never hides a local error.
#[allow(clippy::too_many_arguments)]
pub fn advance_candidate<C: FnMut() -> bool>(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    paths: GasTransportPaths<'_>,
    links: &[DiffuseChannelNetworkLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    step_controls: GasTransportSolve,
    controls: EffusionAdvanceControls,
    cancelled: C,
) -> Result<DiffuseChannelAdvanceCandidate, Error> {
    advance_candidate_with_row_limit(
        parcels,
        paths,
        links,
        duration_s,
        endpoints,
        step_controls,
        controls,
        MAX_EFFUSION_RETAINED_ROWS,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
pub(super) fn advance_candidate_with_row_limit<C: FnMut() -> bool>(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    paths: GasTransportPaths<'_>,
    links: &[DiffuseChannelNetworkLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    step_controls: GasTransportSolve,
    controls: EffusionAdvanceControls,
    maximum_retained_rows: usize,
    mut cancelled: C,
) -> Result<DiffuseChannelAdvanceCandidate, Error> {
    poll(&mut cancelled)?;
    controls.validate()?;
    if !positive(duration_s)
        || parcels.is_empty()
        || parcels.len() > MAX_SNAPSHOT_PARCELS
        || !(1..=MAX_EFFUSION_RETAINED_ROWS).contains(&maximum_retained_rows)
    {
        return Err(Error::Invalid(
            "channel advance needs a positive duration and bounded gases/history",
        ));
    }
    let mut current = parcels.clone();
    let mut elapsed = 0.0;
    let mut ceiling = controls.initial_step_s;
    let (mut attempted, mut accepted, mut rejected, mut retained_rows) = (0, 0, 0, 0usize);
    let mut maximum_error = 0.0_f64;
    let mut history = Vec::new();
    let mut accounts = accounts::Accepted::default();

    let evaluate = |state: &BTreeMap<u64, ClosedGasParcel>,
                    span: f64,
                    elapsed: f64,
                    attempted: &mut usize,
                    accepted: usize,
                    cancel: &mut C|
     -> Result<DiffuseChannelNetworkCandidate, Error> {
        poll(cancel)?;
        if *attempted == controls.maximum_attempts {
            return Err(Error::WorkBudget {
                elapsed_s: elapsed,
                attempted_steps: *attempted,
                accepted_intervals: accepted,
            });
        }
        *attempted += 1;
        step_candidate(state, paths, links, span, endpoints, step_controls, cancel)
            .map_err(|error| physical(error, elapsed, span))
    };

    while elapsed < duration_s {
        poll(&mut cancelled)?;
        if accepted == controls.maximum_accepted_intervals {
            return Err(Error::WorkBudget {
                elapsed_s: elapsed,
                attempted_steps: attempted,
                accepted_intervals: accepted,
            });
        }
        let remaining = duration_s - elapsed;
        let next_time = if ceiling >= remaining {
            duration_s
        } else {
            elapsed + ceiling
        };
        let span = next_time - elapsed;
        let half = 0.5 * span;
        if !positive(span) || !positive(half) || half >= span || half + half != span {
            return Err(Error::Invalid(
                "channel clock or half step is unrepresentable",
            ));
        }
        let trial = (|| {
            let coarse = evaluate(
                &current,
                span,
                elapsed,
                &mut attempted,
                accepted,
                &mut cancelled,
            )?;
            let first = evaluate(
                &current,
                half,
                elapsed,
                &mut attempted,
                accepted,
                &mut cancelled,
            )?;
            let fine = evaluate(
                first.parcels(),
                half,
                elapsed + half,
                &mut attempted,
                accepted,
                &mut cancelled,
            )?;
            let error = state_error(
                &current,
                coarse.parcels(),
                fine.parcels(),
                controls,
                &mut cancelled,
            )?
            .max(accounts::stream_error(
                coarse.step(),
                first.step(),
                fine.step(),
                controls,
                &mut cancelled,
            )?);
            Ok::<_, Error>((first, fine, error))
        })();
        match trial {
            Ok((first, fine, error)) if error <= 1.0 => {
                let rows = accounts::row_count(first.step())
                    .checked_add(accounts::row_count(fine.step()))
                    .and_then(|rows| retained_rows.checked_add(rows));
                if rows.is_none_or(|rows| rows > maximum_retained_rows) {
                    return Err(Error::WorkBudget {
                        elapsed_s: elapsed,
                        attempted_steps: attempted,
                        accepted_intervals: accepted,
                    });
                }
                retained_rows = rows.unwrap();
                accounts.record(first.step(), &mut cancelled)?;
                accounts.record(fine.step(), &mut cancelled)?;
                let (_, first_step) = first.into_parts();
                let (next, fine_step) = fine.into_parts();
                history.push(first_step);
                history.push(fine_step);
                current = next;
                elapsed = next_time;
                accepted += 1;
                maximum_error = maximum_error.max(error);
                let factor = if error == 0.0 {
                    2.0
                } else {
                    (0.9 / fs_math::det::sqrt(error)).clamp(0.5, 2.0)
                };
                ceiling = (span.min(controls.maximum_step_s / factor) * factor)
                    .clamp(controls.minimum_step_s, controls.maximum_step_s);
            }
            Ok((_, _, error)) => {
                rejected += 1;
                let smaller = (0.5 * span).max(controls.minimum_step_s);
                if smaller >= span {
                    return Err(Error::MinimumStep {
                        elapsed_s: elapsed,
                        duration_s: span,
                        normalized_error: error,
                    });
                }
                ceiling = smaller;
            }
            Err(error) if retryable(&error) => {
                rejected += 1;
                let smaller = (0.5 * span).max(controls.minimum_step_s);
                if smaller >= span {
                    return Err(error);
                }
                ceiling = smaller;
            }
            Err(error) => return Err(error),
        }
    }
    let (walls, balance) =
        accounts.finish(parcels, &current, duration_s, step_controls, &mut cancelled)?;
    poll(&mut cancelled)?;
    Ok(DiffuseChannelAdvanceCandidate {
        parcels: current,
        step: DiffuseChannelAdvanceStep {
            duration_s,
            attempted_steps: attempted,
            accepted_intervals: accepted,
            rejected_intervals: rejected,
            maximum_accepted_error: maximum_error,
            accepted_steps: history,
            walls,
            balance,
        },
    })
}

#[cfg(test)]
mod tests;
