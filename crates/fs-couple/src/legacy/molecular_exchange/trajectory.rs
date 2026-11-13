//! Adaptive finite molecular exchange through the original transactional owner.
//!
//! One complete step is compared with two original half steps. Only the latter
//! enter accepted history. Both finite gases, species and relative drifts are
//! resampled at every trial; no continuous RHS or second caloric law is used.
//! This controls a first-order local estimate, not a global accuracy bound.
//! Refinement resolves this aperture operator. Other boundary traction and gas
//! transport must be composed before claiming isolated-vessel equilibration,
//! entropy production or passivity; no compensating heat is invented here.

use super::{
    EffusiveDriftInterface, GasEnergyError, MAX_EFFUSIVE_SPECIES, MolecularExchangeControls,
    MolecularExchangeError, MolecularExchangeStep, MolecularExchangeSurface, MovingGasParcel,
    step_candidate,
};
use fs_math::dd::Dd;
use fs_thermochem::Nasa9StandardStateModelV1;

/// Bounded trial intervals, work and dimensional local error allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularExchangeAdaptiveControls {
    pub initial_step_s: f64,
    /// Minimum retry interval. Accepted native steps are half a trial interval;
    /// a shorter final remainder is permitted only if it needs no retry.
    pub minimum_step_s: f64,
    pub maximum_step_s: f64,
    pub maximum_attempts: usize,
    /// Each accepted interval stores exactly TWO original half-step reports.
    pub maximum_accepted_intervals: usize,
    pub relative_tolerance: f64,
    pub temperature_tolerance_k: f64,
    pub pressure_tolerance_pa: f64,
    pub species_mass_tolerance_kg: f64,
    /// Relative velocity scaling uses each trial's INCREMENT, not common drift.
    pub velocity_tolerance_m_s: f64,
}

#[derive(Debug)]
pub enum MolecularExchangeTrajectoryError {
    Invalid(&'static str),
    Inventory {
        elapsed_s: f64,
        error: GasEnergyError,
    },
    Step {
        elapsed_s: f64,
        duration_s: f64,
        error: MolecularExchangeError,
    },
    MinimumStep {
        elapsed_s: f64,
        duration_s: f64,
        normalized_error: f64,
    },
    WorkBudget {
        elapsed_s: f64,
        attempted_intervals: usize,
        accepted_intervals: usize,
    },
    Unrepresentable {
        elapsed_s: f64,
        detail: &'static str,
    },
    Cancelled {
        elapsed_s: f64,
    },
}
impl core::fmt::Display for MolecularExchangeTrajectoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive finite molecular exchange refused: {self:?}")
    }
}
impl std::error::Error for MolecularExchangeTrajectoryError {}
type Error = MolecularExchangeTrajectoryError;

/// Compensated sums from accepted original reports only. Gross directions are
/// [left-to-right, right-to-left]; relative energies use each step's ORIGINAL
/// local enthalpy reference. Gross formation energies are reporting coordinates,
/// never the operands of the small conservation residuals below.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularExchangeTrajectoryTotals {
    pub gross_mass_kg: [f64; 2],
    pub gross_momentum_relative_kg_m_s: [[f64; 3]; 2],
    pub gross_energy_j: [f64; 2],
    pub net_mass_left_to_right_kg: f64,
    pub net_momentum_left_to_right_n_s: [f64; 3],
    pub relative_energy_into_gas_j: [f64; 2],
    pub kinetic_energy_change_j: [f64; 2],
    pub mass_residual_kg: f64,
    pub momentum_residual_n_s: [f64; 3],
    pub stream_energy_residual_j: f64,
    pub energy_residual_j: f64,
    pub laboratory_energy_residual_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MolecularExchangeTrajectoryStep {
    pub duration_s: f64,
    pub attempted_intervals: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub maximum_accepted_error: f64,
    /// Chronological native half steps. No coarse or rejected report survives.
    pub accepted_steps: Vec<MolecularExchangeStep>,
    pub totals: MolecularExchangeTrajectoryTotals,
}

/// The complete requested window publishes both finite endpoints together.
#[derive(Clone, Debug, PartialEq)]
pub struct MolecularExchangeTrajectoryCandidate {
    gases: [MovingGasParcel; 2],
    step: MolecularExchangeTrajectoryStep,
}
impl MolecularExchangeTrajectoryCandidate {
    pub const fn gases(&self) -> &[MovingGasParcel; 2] {
        &self.gases
    }
    pub const fn step(&self) -> &MolecularExchangeTrajectoryStep {
        &self.step
    }
    pub fn into_parts(self) -> ([MovingGasParcel; 2], MolecularExchangeTrajectoryStep) {
        (self.gases, self.step)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool, elapsed_s: f64) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled { elapsed_s })
    } else {
        Ok(())
    }
}
fn unresolved(elapsed_s: f64, detail: &'static str) -> Error {
    Error::Unrepresentable { elapsed_s, detail }
}
fn physical(error: MolecularExchangeError, elapsed_s: f64, duration_s: f64) -> Error {
    match error {
        MolecularExchangeError::Cancelled => Error::Cancelled { elapsed_s },
        error => Error::Step {
            elapsed_s,
            duration_s,
            error,
        },
    }
}
fn ratio(
    full: f64,
    half: f64,
    scale: f64,
    absolute: f64,
    relative: f64,
    at: f64,
) -> Result<f64, Error> {
    let limit = relative.mul_add(scale, absolute);
    let error = (half - full).abs() / limit;
    if !limit.is_finite() || limit <= 0.0 || !error.is_finite() {
        return Err(unresolved(at, "local error or dimensional allowance"));
    }
    Ok(error)
}

// One ORIGINAL inventory query per gas/state. Missing catalog species retain
// zero mass and enter error control immediately; no per-species map rebuilding.
fn masses(
    gas: &MovingGasParcel,
    axis: &[Nasa9StandardStateModelV1],
    at: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; MAX_EFFUSIVE_SPECIES], Error> {
    let inventory = gas
        .parcel()
        .species_inventory_kg()
        .map_err(|error| Error::Inventory {
            elapsed_s: at,
            error,
        })?;
    if axis.len() > MAX_EFFUSIVE_SPECIES || inventory.len() > MAX_EFFUSIVE_SPECIES {
        return Err(unresolved(at, "original molecular species work bound"));
    }
    let mut result = [0.0; MAX_EFFUSIVE_SPECIES];
    for (species, mass) in inventory {
        poll(cancelled, at)?;
        let index = axis
            .binary_search_by(|s| s.species().cmp(&species))
            .map_err(|_| {
                unresolved(at, "endpoint species absent from original aperture catalog")
            })?;
        result[index] = mass;
    }
    Ok(result)
}

fn normalized_error(
    before: &[MovingGasParcel; 2],
    coarse: &[MovingGasParcel; 2],
    fine: &[MovingGasParcel; 2],
    axis: &[Nasa9StandardStateModelV1],
    controls: MolecularExchangeAdaptiveControls,
    at: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    let mut error = 0.0_f64;
    for side in 0..2 {
        poll(cancelled, at)?;
        let states = [&before[side], &coarse[side], &fine[side]];
        let [old, full, half] = states.map(|s| s.parcel().state().state());
        for (a, b, c, tolerance) in [
            (
                old.temperature,
                full.temperature,
                half.temperature,
                controls.temperature_tolerance_k,
            ),
            (
                old.pressure,
                full.pressure,
                half.pressure,
                controls.pressure_tolerance_pa,
            ),
        ] {
            error = error.max(ratio(
                b,
                c,
                a.abs().max(b.abs()).max(c.abs()),
                tolerance,
                controls.relative_tolerance,
                at,
            )?);
        }
        let [old, full, half] = states.map(MovingGasParcel::velocity_m_s);
        for d in 0..3 {
            let a = full[d] - old[d];
            let b = half[d] - old[d];
            error = error.max(ratio(
                a,
                b,
                a.abs().max(b.abs()),
                controls.velocity_tolerance_m_s,
                controls.relative_tolerance,
                at,
            )?);
        }
        // A closed geometric path must not inspect the inactive interface catalog.
        if !axis.is_empty() {
            let [old, full, half] = [
                masses(states[0], axis, at, cancelled)?,
                masses(states[1], axis, at, cancelled)?,
                masses(states[2], axis, at, cancelled)?,
            ];
            for i in 0..axis.len() {
                poll(cancelled, at)?;
                error = error.max(ratio(
                    full[i],
                    half[i],
                    old[i].max(full[i]).max(half[i]),
                    controls.species_mass_tolerance_kg,
                    controls.relative_tolerance,
                    at,
                )?);
            }
        }
    }
    Ok(error)
}

// The coordinates are summed before rounding into the public receipt. No
// absolute formation-energy subtraction is used for the residual coordinates.
fn accumulate(sum: &mut [Dd; 25], step: &MolecularExchangeStep, at: f64) -> Result<(), Error> {
    for row in &step.exchange.species {
        for side in 0..2 {
            sum[side] = sum[side] + Dd::from_f64(row.mass_kg[side]);
            for d in 0..3 {
                sum[2 + 3 * side + d] =
                    sum[2 + 3 * side + d] + Dd::from_f64(row.momentum_relative_kg_m_s[side][d]);
            }
            sum[8 + side] = sum[8 + side] + Dd::from_f64(row.energy_j[side]);
            sum[14 + side] = sum[14 + side] + Dd::from_f64(row.relative_energy_into_gas_j[side]);
        }
    }
    sum[10] = sum[10] + Dd::from_f64(step.net_mass_left_to_right_kg);
    for d in 0..3 {
        sum[11 + d] = sum[11 + d] + Dd::from_f64(step.net_momentum_left_to_right_n_s[d]);
        sum[19 + d] = sum[19 + d] + Dd::from_f64(step.momentum_residual_n_s[d]);
    }
    for side in 0..2 {
        sum[16 + side] = sum[16 + side] + Dd::from_f64(step.gas[side].kinetic_energy_change_j);
    }
    for (index, value) in [
        (18, step.mass_residual_kg),
        (22, step.stream_energy_residual_j),
        (23, step.energy_residual_j),
        (24, step.laboratory_energy_residual_j),
    ] {
        sum[index] = sum[index] + Dd::from_f64(value);
    }
    if sum
        .iter()
        .any(|s| !s.hi.is_finite() || !s.lo.is_finite() || !s.to_f64().is_finite())
    {
        return Err(unresolved(at, "accepted molecular account total"));
    }
    Ok(())
}
fn totals(sum: [Dd; 25]) -> MolecularExchangeTrajectoryTotals {
    let s = sum.map(Dd::to_f64);
    MolecularExchangeTrajectoryTotals {
        gross_mass_kg: [s[0], s[1]],
        gross_momentum_relative_kg_m_s: [[s[2], s[3], s[4]], [s[5], s[6], s[7]]],
        gross_energy_j: [s[8], s[9]],
        net_mass_left_to_right_kg: s[10],
        net_momentum_left_to_right_n_s: [s[11], s[12], s[13]],
        relative_energy_into_gas_j: [s[14], s[15]],
        kinetic_energy_change_j: [s[16], s[17]],
        mass_residual_kg: s[18],
        momentum_residual_n_s: [s[19], s[20], s[21]],
        stream_energy_residual_j: s[22],
        energy_residual_j: s[23],
        laboratory_energy_residual_j: s[24],
    }
}

/// Evolve one complete positive interval with a held aperture and original
/// sources. Every full/half trial refreshes BOTH finite gas distributions.
/// Only typed Duration/Withdrawal refusals and local error permit halving;
/// source, conservation and caloric failures propagate with their physical time.
/// A small accepted error (<=1/4) permits doubling. All accepted native half
/// steps retain the original physical budgets and endpoint source admission.
/// Cancellation or exhausted bounds publishes no partial gas state or history.
#[allow(clippy::too_many_arguments)]
pub fn advance_candidate(
    gases: [&MovingGasParcel; 2],
    interface: &EffusiveDriftInterface,
    surface: MolecularExchangeSurface,
    duration_s: f64,
    controls: MolecularExchangeControls,
    adaptive: MolecularExchangeAdaptiveControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MolecularExchangeTrajectoryCandidate, Error> {
    poll(&mut cancelled, 0.0)?;
    if !duration_s.is_finite()
        || duration_s <= 0.0
        || ![
            adaptive.initial_step_s,
            adaptive.minimum_step_s,
            adaptive.maximum_step_s,
            adaptive.temperature_tolerance_k,
            adaptive.pressure_tolerance_pa,
            adaptive.species_mass_tolerance_kg,
            adaptive.velocity_tolerance_m_s,
        ]
        .into_iter()
        .all(|v| v.is_finite() && v > 0.0)
        || adaptive.minimum_step_s > adaptive.initial_step_s
        || adaptive.initial_step_s > adaptive.maximum_step_s
        || !adaptive.relative_tolerance.is_finite()
        || adaptive.relative_tolerance < 0.0
        || adaptive.maximum_attempts == 0
        || adaptive.maximum_accepted_intervals == 0
        || adaptive.maximum_accepted_intervals.checked_mul(2).is_none()
    {
        return Err(Error::Invalid(
            "positive bounded time, dimensional tolerances and work limits required",
        ));
    }
    let axis = if surface.area_m2 > 0.0 {
        interface
            .momentum_interface()
            .interface()
            .gas_catalog()
            .thermodynamics()
            .components()
    } else {
        &[]
    };
    let mut current = gases.map(Clone::clone);
    let (mut elapsed, mut next_step, mut maximum_error) = (0.0, adaptive.initial_step_s, 0.0_f64);
    let (mut attempts, mut rejected, mut accepted) = (0, 0, 0);
    let mut history = Vec::new();
    let mut sum = [Dd::ZERO; 25];
    while elapsed < duration_s {
        poll(&mut cancelled, elapsed)?;
        if attempts == adaptive.maximum_attempts || accepted == adaptive.maximum_accepted_intervals
        {
            return Err(Error::WorkBudget {
                elapsed_s: elapsed,
                attempted_intervals: attempts,
                accepted_intervals: accepted,
            });
        }
        let remaining = duration_s - elapsed;
        let mut h = next_step.min(remaining);
        let tail = remaining - h;
        if tail > 0.0
            && tail < adaptive.minimum_step_s
            && 0.5 * remaining >= adaptive.minimum_step_s
        {
            h = 0.5 * remaining;
        }
        let half = 0.5 * h;
        let end = elapsed + h;
        let middle_time = elapsed + half;
        if half <= 0.0
            || half + half != h
            || end <= elapsed
            || end > duration_s
            || middle_time <= elapsed
            || middle_time >= end
            || (h == remaining && end != duration_s)
        {
            return Err(unresolved(
                elapsed,
                "trial halves or final physical time cannot advance exactly",
            ));
        }
        attempts += 1;
        let trial = (|| {
            let mut step = |state: &[MovingGasParcel; 2], dt, at| {
                step_candidate(
                    [&state[0], &state[1]],
                    interface,
                    surface,
                    dt,
                    controls,
                    &mut cancelled,
                )
                .map_err(|error| physical(error, at, dt))
            };
            let coarse = step(&current, h, elapsed)?;
            let middle = step(&current, half, elapsed)?;
            let fine = step(middle.gases(), half, middle_time)?;
            let error = normalized_error(
                &current,
                coarse.gases(),
                fine.gases(),
                axis,
                adaptive,
                elapsed,
                &mut cancelled,
            )?;
            let (_, first) = middle.into_parts();
            let (fine, second) = fine.into_parts();
            Ok::<_, Error>((fine, [first, second], error))
        })();
        let (fine, steps, error) = match trial {
            Ok(result) => result,
            Err(Error::Step {
                error:
                    MolecularExchangeError::Duration { .. } | MolecularExchangeError::Withdrawal { .. },
                ..
            }) if half >= adaptive.minimum_step_s => {
                rejected += 1;
                next_step = half;
                continue;
            }
            Err(error) => return Err(error),
        };
        if error > 1.0 {
            rejected += 1;
            if half < adaptive.minimum_step_s {
                return Err(Error::MinimumStep {
                    elapsed_s: elapsed,
                    duration_s: h,
                    normalized_error: error,
                });
            }
            next_step = half;
            continue;
        }
        for step in steps {
            poll(&mut cancelled, elapsed)?;
            accumulate(&mut sum, &step, end)?;
            history.push(step);
        }
        current = fine;
        elapsed = end;
        accepted += 1;
        maximum_error = maximum_error.max(error);
        next_step = if error <= 0.25 {
            (2.0 * h).min(adaptive.maximum_step_s)
        } else {
            h
        };
    }
    poll(&mut cancelled, elapsed)?;
    Ok(MolecularExchangeTrajectoryCandidate {
        gases: current,
        step: MolecularExchangeTrajectoryStep {
            duration_s: elapsed,
            attempted_intervals: attempts,
            accepted_intervals: accepted,
            rejected_intervals: rejected,
            maximum_accepted_error: maximum_error,
            accepted_steps: history,
            totals: totals(sum),
        },
    })
}
