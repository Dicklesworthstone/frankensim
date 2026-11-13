//! Adaptive simultaneous molecular exchange on finite moving-gas networks.
//!
//! Every trial calls the original network owner. Only two accepted half steps
//! survive; changing composition, temperature and drift feed the next sample.
//! Each link's two gross streams and authoritative signed net are controlled
//! separately, in addition to the original gas-state estimates. Fixed window
//! budgets audit actual endpoints and accumulated absolute numerical defects.
//! This refines a first-order interface operator, not the omitted wall forces,
//! container motion, molecular flight time or spatial gas transport.

pub(super) mod accounts;

use super::super::{GasEnergyError, MAX_EFFUSIVE_SPECIES};
use super::{
    MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS, MolecularNetworkControls, MolecularNetworkError,
    MolecularNetworkLink, MolecularNetworkStep, MovingGasBoundarySolve, MovingGasParcel,
    step_candidate,
};
use accounts::{Sources, Totals, window_balance};
use fs_thermochem::{SpeciesId, ThermochemErrorV1};
use std::collections::BTreeMap;

pub use super::super::trajectory::MolecularExchangeAdaptiveControls;

/// Existing gas-state controls plus stream estimates and fixed physical budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularNetworkAdaptiveControls {
    pub states: MolecularExchangeAdaptiveControls,
    /// Per component, in `energy_frame_velocity_m_s`'s translating frame.
    pub stream_momentum_tolerance_n_s: f64,
    /// Original species formation baselines are removed before arithmetic.
    pub stream_energy_tolerance_j: f64,
    /// Applies independently to total and each species' actual mass closure.
    pub window_mass_tolerance_kg: f64,
    /// Per component in the original inertial momentum frame.
    pub window_momentum_tolerance_n_s: f64,
    /// Applies to both common-frame and laboratory-frame energy closure.
    pub window_energy_tolerance_j: f64,
}

/// A beam or signed net in the explicitly declared common energy frame.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MolecularNetworkStreamMoments {
    pub mass_kg: f64,
    pub momentum_n_s: [f64; 3],
    /// Complete molecular energy with one fixed original u_i(T_reference)
    /// removed per species; no absolute formation-energy subtraction is used.
    pub reference_energy_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MolecularNetworkLinkTotals {
    pub id: u64,
    pub gases: [u64; 2],
    /// [left-to-right, right-to-left], preserving each beam's vector direction.
    pub gross: [MolecularNetworkStreamMoments; 2],
    /// Original small-drift net mass and energy coordinates remain authoritative.
    pub net_left_to_right: MolecularNetworkStreamMoments,
}

/// Signed actual closure, or nonnegative accumulated numerical defects.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MolecularNetworkWindowBalance {
    pub mass_kg: f64,
    pub maximum_species_mass_kg: f64,
    pub momentum_n_s: [f64; 3],
    pub energy_j: f64,
    pub laboratory_energy_j: f64,
}

#[derive(Debug)]
pub enum MolecularNetworkTrajectoryError {
    Invalid(&'static str),
    Inventory {
        gas: u64,
        elapsed_s: f64,
        error: GasEnergyError,
    },
    Caloric {
        species: SpeciesId,
        elapsed_s: f64,
        error: ThermochemErrorV1,
    },
    Step {
        elapsed_s: f64,
        duration_s: f64,
        error: MolecularNetworkError,
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
    WindowBalance {
        elapsed_s: f64,
        residual: MolecularNetworkWindowBalance,
        debit: MolecularNetworkWindowBalance,
    },
    Unrepresentable {
        elapsed_s: f64,
        detail: &'static str,
    },
    Cancelled {
        elapsed_s: f64,
    },
}
impl core::fmt::Display for MolecularNetworkTrajectoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive molecular network refused: {self:?}")
    }
}
impl std::error::Error for MolecularNetworkTrajectoryError {}
type Error = MolecularNetworkTrajectoryError;

#[derive(Clone, Debug, PartialEq)]
pub struct MolecularNetworkTrajectoryStep {
    pub duration_s: f64,
    pub energy_frame_velocity_m_s: [f64; 3],
    pub attempted_intervals: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub maximum_accepted_error: f64,
    /// Exactly two chronological native reports per accepted interval.
    pub accepted_steps: Vec<MolecularNetworkStep>,
    /// Canonical link order; no coarse or rejected stream contributes.
    pub links: Vec<MolecularNetworkLinkTotals>,
    pub residual: MolecularNetworkWindowBalance,
    /// Absolute native defects and actual step closure, not a step-count-scaled
    /// allowance. Opposite signs cannot refund this fixed physical budget.
    pub debit: MolecularNetworkWindowBalance,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MolecularNetworkTrajectoryCandidate {
    gases: BTreeMap<u64, MovingGasParcel>,
    step: MolecularNetworkTrajectoryStep,
}
impl MolecularNetworkTrajectoryCandidate {
    pub fn gases(&self) -> &BTreeMap<u64, MovingGasParcel> {
        &self.gases
    }
    pub const fn step(&self) -> &MolecularNetworkTrajectoryStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        BTreeMap<u64, MovingGasParcel>,
        MolecularNetworkTrajectoryStep,
    ) {
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
fn physical(error: MolecularNetworkError, elapsed_s: f64, duration_s: f64) -> Error {
    match error {
        MolecularNetworkError::Cancelled => Error::Cancelled { elapsed_s },
        error => Error::Step {
            elapsed_s,
            duration_s,
            error,
        },
    }
}
pub(super) fn ratio(
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

pub(super) fn state_error(
    before: &BTreeMap<u64, MovingGasParcel>,
    coarse: &BTreeMap<u64, MovingGasParcel>,
    fine: &BTreeMap<u64, MovingGasParcel>,
    sources: &Sources<'_>,
    controls: MolecularExchangeAdaptiveControls,
    at: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    let mut error = 0.0_f64;
    for (&id, initial) in before {
        poll(cancelled, at)?;
        let states = [initial, &coarse[&id], &fine[&id]];
        let [old, full, half] = states.map(|g| g.parcel().state().state());
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
        if let Some(axis) = sources.axes.get(&id) {
            let [old, full, half] = [
                sources.masses(id, states[0], at, cancelled)?,
                sources.masses(id, states[1], at, cancelled)?,
                sources.masses(id, states[2], at, cancelled)?,
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

pub(super) fn admit(
    gases: &BTreeMap<u64, MovingGasParcel>,
    links: &[MolecularNetworkLink<'_>],
    endpoints: &BTreeMap<u64, MovingGasBoundarySolve>,
    duration: f64,
    physical: MolecularNetworkControls,
    adaptive: MolecularNetworkAdaptiveControls,
) -> Result<(), Error> {
    let a = adaptive.states;
    if gases.is_empty()
        || gases.len() > MAX_SNAPSHOT_PARCELS
        || links.len() > MAX_SNAPSHOT_TRANSFERS
        || !gases.keys().eq(endpoints.keys())
        || ![
            duration,
            a.initial_step_s,
            a.minimum_step_s,
            a.maximum_step_s,
            a.temperature_tolerance_k,
            a.pressure_tolerance_pa,
            a.species_mass_tolerance_kg,
            a.velocity_tolerance_m_s,
            adaptive.stream_momentum_tolerance_n_s,
            adaptive.stream_energy_tolerance_j,
            adaptive.window_mass_tolerance_kg,
            adaptive.window_momentum_tolerance_n_s,
            adaptive.window_energy_tolerance_j,
        ]
        .into_iter()
        .all(|v| v.is_finite() && v > 0.0)
        || a.minimum_step_s > a.initial_step_s
        || a.initial_step_s > a.maximum_step_s
        || !a.relative_tolerance.is_finite()
        || a.relative_tolerance < 0.0
        || !physical.maximum_duration_s.is_finite()
        || physical.maximum_duration_s < a.minimum_step_s
        || a.maximum_attempts == 0
        || a.maximum_accepted_intervals == 0
        || a.maximum_accepted_intervals
            .checked_mul(2)
            .and_then(|n| n.checked_mul(links.len().max(gases.len())))
            .is_none()
    {
        return Err(Error::Invalid(
            "bounded original gas/link IDs, positive times/tolerances and finite work/history limits required",
        ));
    }
    Ok(())
}

pub(super) fn admit_window(
    residual: MolecularNetworkWindowBalance,
    debit: MolecularNetworkWindowBalance,
    controls: MolecularNetworkAdaptiveControls,
    elapsed_s: f64,
) -> Result<(), Error> {
    let exceeds = |row: MolecularNetworkWindowBalance| {
        row.mass_kg.abs() > controls.window_mass_tolerance_kg
            || row.maximum_species_mass_kg > controls.window_mass_tolerance_kg
            || row
                .momentum_n_s
                .iter()
                .any(|v| v.abs() > controls.window_momentum_tolerance_n_s)
            || row.energy_j.abs() > controls.window_energy_tolerance_j
            || row.laboratory_energy_j.abs() > controls.window_energy_tolerance_j
    };
    if exceeds(residual) || exceeds(debit) {
        Err(Error::WindowBalance {
            elapsed_s,
            residual,
            debit,
        })
    } else {
        Ok(())
    }
}

/// Evolve the complete window with held topology, areas, normals and aperture
/// velocities. Every native trial refreshes all original finite gas states.
/// Incoming species enter error control from exact zero. Only the typed gross
/// withdrawal error permits a physical retry; source/caloric/closure failures
/// remain hard refusals. The original maximum duration caps each trial.
///
/// Time estimates control local first-order error, not a certified global bound.
/// All native inventory, source and drift support checks remain authoritative.
/// History is bounded by twice the accepted-interval cap times the original
/// bounded network size. No caller-owned state changes, even on a late refusal.
#[allow(clippy::too_many_arguments)]
pub fn advance_candidate(
    gases: &BTreeMap<u64, MovingGasParcel>,
    links: &[MolecularNetworkLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, MovingGasBoundarySolve>,
    controls: MolecularNetworkControls,
    adaptive: MolecularNetworkAdaptiveControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MolecularNetworkTrajectoryCandidate, Error> {
    poll(&mut cancelled, 0.0)?;
    admit(gases, links, endpoints, duration_s, controls, adaptive)?;
    let sources = Sources::new(gases, links, &mut cancelled)?;
    let a = adaptive.states;
    let mut current = gases.clone();
    let (mut elapsed, mut next_step, mut maximum_error) = (0.0, a.initial_step_s, 0.0_f64);
    let (mut attempts, mut rejected, mut accepted) = (0, 0, 0);
    let mut history = Vec::new();
    let mut totals = Totals::default();
    let mut residual = MolecularNetworkWindowBalance::default();
    while elapsed < duration_s {
        poll(&mut cancelled, elapsed)?;
        if attempts == a.maximum_attempts || accepted == a.maximum_accepted_intervals {
            return Err(Error::WorkBudget {
                elapsed_s: elapsed,
                attempted_intervals: attempts,
                accepted_intervals: accepted,
            });
        }
        let remaining = duration_s - elapsed;
        let mut h = next_step.min(controls.maximum_duration_s).min(remaining);
        let tail = remaining - h;
        if tail > 0.0 && tail < a.minimum_step_s && 0.5 * remaining >= a.minimum_step_s {
            h = 0.5 * remaining;
        }
        let (half, end, middle_time) = (0.5 * h, elapsed + h, elapsed + 0.5 * h);
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
            let mut step = |state: &BTreeMap<u64, MovingGasParcel>, dt, at| {
                step_candidate(state, links, dt, endpoints, controls, &mut cancelled)
                    .map_err(|error| physical(error, at, dt))
            };
            let coarse = step(&current, h, elapsed)?;
            let middle = step(&current, half, elapsed)?;
            let fine = step(middle.gases(), half, middle_time)?;
            let error = state_error(
                &current,
                coarse.gases(),
                fine.gases(),
                &sources,
                a,
                elapsed,
                &mut cancelled,
            )?
            .max(accounts::stream_error(
                &current,
                middle.gases(),
                coarse.step(),
                middle.step(),
                fine.step(),
                &sources,
                adaptive,
                elapsed,
                &mut cancelled,
            )?);
            let (middle, first) = middle.into_parts();
            let (fine, second) = fine.into_parts();
            Ok::<_, Error>((middle, fine, [first, second], error))
        })();
        let (middle, fine, steps, error) = match trial {
            Ok(result) => result,
            Err(Error::Step {
                error: MolecularNetworkError::Withdrawal { .. },
                ..
            }) if half >= a.minimum_step_s => {
                rejected += 1;
                next_step = half;
                continue;
            }
            Err(error) => return Err(error),
        };
        if error > 1.0 {
            rejected += 1;
            if half < a.minimum_step_s {
                return Err(Error::MinimumStep {
                    elapsed_s: elapsed,
                    duration_s: h,
                    normalized_error: error,
                });
            }
            next_step = half;
            continue;
        }
        totals.accept(
            &current,
            &middle,
            &steps[0],
            &sources,
            middle_time,
            &mut cancelled,
        )?;
        totals.accept(&middle, &fine, &steps[1], &sources, end, &mut cancelled)?;
        residual = window_balance(
            gases,
            &fine,
            &sources,
            controls.energy_frame_velocity_m_s,
            end,
            &mut cancelled,
        )?;
        admit_window(residual, totals.debit(end)?, adaptive, end)?;
        poll(&mut cancelled, end)?;
        history.extend(steps);
        current = fine;
        elapsed = end;
        accepted += 1;
        maximum_error = maximum_error.max(error);
        next_step = if error <= 0.25 {
            (2.0 * h).min(a.maximum_step_s)
        } else {
            h
        };
    }
    let debit = totals.debit(elapsed)?;
    let links = totals.links(elapsed, &mut cancelled)?;
    poll(&mut cancelled, elapsed)?;
    Ok(MolecularNetworkTrajectoryCandidate {
        gases: current,
        step: MolecularNetworkTrajectoryStep {
            duration_s: elapsed,
            energy_frame_velocity_m_s: controls.energy_frame_velocity_m_s,
            attempted_intervals: attempts,
            accepted_intervals: accepted,
            rejected_intervals: rejected,
            maximum_accepted_error: maximum_error,
            accepted_steps: history,
            links,
            residual,
            debit,
        },
    })
}
