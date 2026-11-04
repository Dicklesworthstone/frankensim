//! Adaptive effusion on the original simultaneous species/energy transaction.
//!
//! Compare one full step with two half steps and accept only the halves. Each
//! native call resamples ALL internal paths and external molecular beams from
//! its actual original state. There is no second gas law, clipped withdrawal,
//! interpolated composition, or imposed constant heat capacity.
//!
//! The estimate controls local temperature, pressure, every species inventory,
//! directional species throughput and source-relative boundary energy. It is
//! not a rigorous global error bound. Geometry, held reservoirs and relative
//! source velocities remain fixed over this window; moving boundaries, heat,
//! reactions and time-dependent controls require their enclosing physical owner.

mod accounts;

use super::super::{
    DiffusionNetworkError, MembraneNetworkError, PorousNetworkError, SnapshotNozzleError,
};
use super::{
    reservoir::{
        self, EffusionDriftingReservoirLink, EffusionReservoirLink, ReservoirEffusionCandidate,
        ReservoirEffusionStep,
    },
    vacuum::VacuumEffusionLink,
};
use crate::thermochemical_gas::energy::{
    GasEnergyError, GasEnergySolve,
    parcel::{
        ClosedGasParcel,
        transfer::{
            snapshot::MAX_SNAPSHOT_PARCELS,
            species::{
                SpeciesId,
                network::{
                    SpeciesNetworkError,
                    transport::{self, GasTransportError, GasTransportPaths, GasTransportSolve},
                },
            },
        },
    },
};
use std::collections::BTreeMap;

pub use accounts::{EffusionIntegratedBoundary, EffusionWindowBalance};

/// Hard work/history bounds, in addition to the caller's smaller budgets.
pub const MAX_EFFUSION_ATTEMPTS: usize = 65_536;
pub const MAX_EFFUSION_INTERVALS: usize = 4_096;
pub const MAX_EFFUSION_RETAINED_ROWS: usize = 131_072;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EffusionAdvanceControls {
    pub initial_step_s: f64,
    /// Minimum retry interval; native accepted steps are half an interval.
    /// A shorter final remainder is permitted only if it needs no retry.
    pub minimum_step_s: f64,
    pub maximum_step_s: f64,
    /// ALL native candidate evaluations, including failed and coarse probes.
    /// An accepted interval requires three evaluations.
    pub maximum_attempts: usize,
    pub maximum_accepted_intervals: usize,
    pub relative_tolerance: f64,
    pub temperature_tolerance_k: f64,
    pub pressure_tolerance_pa: f64,
    /// Applied separately to each inventory and each directional stream.
    pub species_mass_tolerance_kg: f64,
    /// Boundary energy uses one common partial-enthalpy reference per trial.
    pub stream_energy_tolerance_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum EffusionAdvanceError {
    Invalid(&'static str),
    Step {
        elapsed_s: f64,
        duration_s: f64,
        error: GasTransportError,
    },
    MinimumStep {
        elapsed_s: f64,
        duration_s: f64,
        normalized_error: f64,
    },
    WorkBudget {
        elapsed_s: f64,
        attempted_steps: usize,
        accepted_intervals: usize,
    },
    Cancelled,
}
impl core::fmt::Display for EffusionAdvanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive material effusion: {self:?}")
    }
}
impl std::error::Error for EffusionAdvanceError {}
use EffusionAdvanceError as Error;

#[derive(Clone, Debug, PartialEq)]
pub struct EffusionAdvanceStep {
    pub duration_s: f64,
    pub attempted_steps: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub maximum_accepted_error: f64,
    /// Chronological ORIGINAL half-step reports. Rejected/coarse probes never
    /// enter this history or the integrated physical accounts.
    pub accepted_steps: Vec<ReservoirEffusionStep>,
    /// Canonical boundary/species sums; gross energy is a reporting coordinate.
    pub boundaries: Vec<EffusionIntegratedBoundary>,
    /// One fixed physical budget over the complete requested duration.
    pub balance: EffusionWindowBalance,
}

/// Either every endpoint reaches the requested duration, or none is published.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusionAdvanceCandidate {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    step: EffusionAdvanceStep,
}
impl EffusionAdvanceCandidate {
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub const fn step(&self) -> &EffusionAdvanceStep {
        &self.step
    }
    pub fn into_parts(self) -> (BTreeMap<u64, ClosedGasParcel>, EffusionAdvanceStep) {
        (self.parcels, self.step)
    }
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
    // Source, geometry, model-domain and unrepresentable-arithmetic errors
    // remain hard refusals. No classification by diagnostic message text.
    matches!(
        error,
        Error::Step {
            error: GasTransportError::Exchange(
                SpeciesNetworkError::Withdrawal { .. }
                    | SpeciesNetworkError::Parcel {
                        error: GasEnergyError::OutsideBracket { .. }
                            | GasEnergyError::NotConverged { .. },
                        ..
                    }
            ) | GasTransportError::Nozzle(SnapshotNozzleError::EquilibriumCrossed { .. })
                | GasTransportError::Porous(PorousNetworkError::EquilibriumCrossed { .. })
                | GasTransportError::Membrane(MembraneNetworkError::EquilibriumCrossed { .. })
                | GasTransportError::Diffusion(DiffusionNetworkError::EntropyDecrease { .. }),
            ..
        }
    )
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

fn ratio(a: f64, b: f64, scale: f64, absolute: f64, relative: f64) -> Result<f64, Error> {
    let allowance = relative.mul_add(scale, absolute);
    let error = (a - b).abs() / allowance;
    if !positive(allowance) || !error.is_finite() {
        Err(Error::Invalid(
            "unrepresentable dimensional effusion error estimate",
        ))
    } else {
        Ok(error)
    }
}

fn inventory(gas: &ClosedGasParcel) -> Result<BTreeMap<SpeciesId, f64>, Error> {
    gas.species_inventory_kg()
        .map_err(|_| Error::Invalid("unrepresentable original species inventory"))
}

fn state_error(
    before: &BTreeMap<u64, ClosedGasParcel>,
    coarse: &BTreeMap<u64, ClosedGasParcel>,
    fine: &BTreeMap<u64, ClosedGasParcel>,
    controls: EffusionAdvanceControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    let mut error = 0.0_f64;
    for (&id, old) in before {
        poll(cancelled)?;
        let (a, b) = (&coarse[&id], &fine[&id]);
        let (s0, s1, s2) = (old.state().state(), a.state().state(), b.state().state());
        for (x, y, z, absolute) in [
            (
                s0.temperature,
                s1.temperature,
                s2.temperature,
                controls.temperature_tolerance_k,
            ),
            (
                s0.pressure,
                s1.pressure,
                s2.pressure,
                controls.pressure_tolerance_pa,
            ),
        ] {
            error = error.max(ratio(
                y,
                z,
                x.abs().max(y.abs()).max(z.abs()),
                absolute,
                controls.relative_tolerance,
            )?);
        }
        // Query each original/native inventory once. Newly supplied species
        // enter the union immediately with exact zero in the original state.
        let (old, a, b) = (inventory(old)?, inventory(a)?, inventory(b)?);
        for species in old.keys().chain(a.keys()).chain(b.keys()) {
            poll(cancelled)?;
            let x = old.get(species).copied().unwrap_or(0.0);
            let y = a.get(species).copied().unwrap_or(0.0);
            let z = b.get(species).copied().unwrap_or(0.0);
            error = error.max(ratio(
                y,
                z,
                x.max(y).max(z),
                controls.species_mass_tolerance_kg,
                controls.relative_tolerance,
            )?);
        }
    }
    Ok(error)
}

/// Advance arbitrary original finite transport paths and molecular boundaries.
///
/// Vacuum, reciprocal prescribed sources and qualified drifting sources share
/// every original gross-withdrawal limit. A large initial step can retry on an
/// inventory or caloric-solve limit; unsupported source states never silently
/// become admissible by retrying. All tolerances retain their physical units.
///
/// Actual original-to-final species inventories and source-relative internal
/// energy are independently checked against accepted streams. Native caloric
/// residuals are also accumulated. Both energy checks must fit ONE original
/// budget, never a step-count multiple. Formation-energy datums never scale an
/// error allowance.
#[allow(clippy::too_many_arguments)]
pub fn advance_candidate<C: FnMut() -> bool>(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    paths: GasTransportPaths<'_>,
    vacuum: &[VacuumEffusionLink<'_>],
    reservoirs: &[EffusionReservoirLink<'_>],
    drifting: &[EffusionDriftingReservoirLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    step_controls: GasTransportSolve,
    controls: EffusionAdvanceControls,
    mut cancelled: C,
) -> Result<EffusionAdvanceCandidate, Error> {
    poll(&mut cancelled)?;
    if !positive(duration_s)
        || parcels.is_empty()
        || parcels.len() > MAX_SNAPSHOT_PARCELS
        || ![
            controls.initial_step_s,
            controls.minimum_step_s,
            controls.maximum_step_s,
            controls.temperature_tolerance_k,
            controls.pressure_tolerance_pa,
            controls.species_mass_tolerance_kg,
            controls.stream_energy_tolerance_j,
        ]
        .into_iter()
        .all(positive)
        || controls.minimum_step_s > controls.initial_step_s
        || controls.initial_step_s > controls.maximum_step_s
        || !controls.relative_tolerance.is_finite()
        || !(0.0..=1.0).contains(&controls.relative_tolerance)
        || !(3..=MAX_EFFUSION_ATTEMPTS).contains(&controls.maximum_attempts)
        || !(1..=MAX_EFFUSION_INTERVALS).contains(&controls.maximum_accepted_intervals)
    {
        return Err(Error::Invalid(
            "positive duration, dimensional tolerances and bounded adaptive effort required",
        ));
    }
    let mut current = parcels.clone();
    let mut elapsed = 0.0;
    let mut ceiling = controls.initial_step_s;
    let mut attempted = 0;
    let mut accepted = 0;
    let mut rejected = 0;
    let mut maximum_error = 0.0_f64;
    let mut history = Vec::new();
    let mut accounts = accounts::Accepted::default();
    let mut retained_rows = 0usize;

    let evaluate = |state: &BTreeMap<u64, ClosedGasParcel>,
                    span: f64,
                    elapsed: f64,
                    attempted: &mut usize,
                    accepted: usize,
                    cancel: &mut C|
     -> Result<ReservoirEffusionCandidate, Error> {
        poll(cancel)?;
        if *attempted == controls.maximum_attempts {
            return Err(Error::WorkBudget {
                elapsed_s: elapsed,
                attempted_steps: *attempted,
                accepted_intervals: accepted,
            });
        }
        *attempted += 1;
        reservoir::step_candidate_with_drift(
            state,
            paths,
            vacuum,
            reservoirs,
            drifting,
            span,
            endpoints,
            step_controls,
            cancel,
        )
        .map_err(|e| physical(e, elapsed, span))
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
                "effusion clock or half step is unrepresentable",
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
                &current,
                &coarse,
                &first,
                &fine,
                controls,
                &mut cancelled,
            )?);
            Ok::<_, Error>((first, fine, error))
        })();
        match trial {
            Ok((first, fine, error)) if error <= 1.0 => {
                let rows = accounts::row_count(first.step())
                    .checked_add(accounts::row_count(fine.step()))
                    .and_then(|n| retained_rows.checked_add(n));
                if rows.is_none_or(|n| n > MAX_EFFUSION_RETAINED_ROWS) {
                    return Err(Error::WorkBudget {
                        elapsed_s: elapsed,
                        attempted_steps: attempted,
                        accepted_intervals: accepted,
                    });
                }
                retained_rows = rows.unwrap();
                accounts.record(&current, &first, &mut cancelled)?;
                accounts.record(first.parcels(), &fine, &mut cancelled)?;
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
                // Clamp before multiplication, so a large finite caller ceiling
                // cannot overflow during harmless step-size growth.
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
    let (boundaries, balance) =
        accounts.finish(parcels, &current, step_controls, &mut cancelled)?;
    if boundaries.is_empty() {
        transport::audit_closed_window(parcels, &current, step_controls, &mut cancelled)
            .map_err(|e| physical(e, duration_s, duration_s))?;
    }
    poll(&mut cancelled)?;
    Ok(EffusionAdvanceCandidate {
        parcels: current,
        step: EffusionAdvanceStep {
            duration_s,
            attempted_steps: attempted,
            accepted_intervals: accepted,
            rejected_intervals: rejected,
            maximum_accepted_error: maximum_error,
            accepted_steps: history,
            boundaries,
            balance,
        },
    })
}

#[cfg(test)]
mod tests;
