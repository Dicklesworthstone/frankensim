//! Adaptive finite material walls coupled to finite molecular-channel gases.
//!
//! Each trial compares one original full step with two resampled half steps and
//! accepts only the half steps. Gas and wall inventories are advanced together:
//! no wall temperature is interpolated, no heat capacity is frozen, and every
//! accepted half step reconstructs the original material law at its actual state.
//!
//! A complete-window audit uses the original gas and wall inventories and one
//! fixed mass/energy/entropy budget. Rejected/coarse probes never enter that
//! audit. This remains a local first-order time-discretization estimate, not a
//! continuum or constitutive accuracy certificate.

use super::{
    FiniteChannelError, FiniteChannelNetworkLink, FiniteChannelNetworkStep, FixedVolumeChannelWall,
    step_candidate,
};
use crate::thermochemical_gas::energy::{
    GasEnergySolve,
    parcel::{
        ClosedGasParcel,
        transfer::species::{
            entropy,
            network::transport::{GasTransportError, GasTransportPaths, GasTransportSolve},
        },
    },
};
use fs_math::dd::Dd;
use fs_thermal::enthalpy::ThermalLaw;
use fs_thermochem::SpeciesId;
use std::collections::{BTreeMap, BTreeSet};

use super::super::super::trajectory::retryable_native;
pub use super::super::trajectory::EffusionAdvanceControls;

#[derive(Clone, Debug, PartialEq)]
pub enum FiniteChannelAdvanceError {
    Invalid(&'static str),
    Step {
        elapsed_s: f64,
        duration_s: f64,
        error: FiniteChannelError,
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
impl core::fmt::Display for FiniteChannelAdvanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive finite molecular channel refusal: {self:?}")
    }
}
impl std::error::Error for FiniteChannelAdvanceError {}
type Error = FiniteChannelAdvanceError;

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteChannelWindowBalance {
    pub total_mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    pub gas_energy_change_j: f64,
    pub wall_energy_change_j: f64,
    pub energy_residual_j: f64,
    pub gas_entropy_change_j_k: f64,
    pub wall_entropy_change_j_k: f64,
    pub total_entropy_change_j_k: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteChannelAdvanceStep {
    pub duration_s: f64,
    pub attempted_steps: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub maximum_accepted_error: f64,
    /// Chronological accepted half-step reports; coarse/rejected probes omitted.
    pub accepted_steps: Vec<FiniteChannelNetworkStep>,
    pub balance: FiniteChannelWindowBalance,
}

#[derive(Clone, Debug)]
pub struct FiniteChannelAdvanceCandidate<'l, 'p> {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    walls: BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    step: FiniteChannelAdvanceStep,
}
impl<'l, 'p> FiniteChannelAdvanceCandidate<'l, 'p> {
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub fn walls(&self) -> &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>> {
        &self.walls
    }
    pub const fn step(&self) -> &FiniteChannelAdvanceStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        BTreeMap<u64, ClosedGasParcel>,
        BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
        FiniteChannelAdvanceStep,
    ) {
        (self.parcels, self.walls, self.step)
    }
}

fn positive(v: f64) -> bool {
    v.is_finite() && v > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn physical(error: FiniteChannelError, elapsed_s: f64, duration_s: f64) -> Error {
    if error == FiniteChannelError::Cancelled {
        Error::Cancelled
    } else {
        Error::Step {
            elapsed_s,
            duration_s,
            error,
        }
    }
}
pub(in super::super) fn retryable(error: &Error) -> bool {
    match error {
        Error::Step {
            error: FiniteChannelError::Transport(e),
            ..
        } => retryable_native(e),
        Error::Step {
            error: FiniteChannelError::WallOvershoot { .. } | FiniteChannelError::WallDomain { .. },
            ..
        } => true,
        _ => false,
    }
}
pub(in super::super) fn ratio(
    a: f64,
    b: f64,
    scale: f64,
    absolute: f64,
    relative: f64,
) -> Result<f64, Error> {
    let allowance = relative.mul_add(scale, absolute);
    let value = (a - b).abs() / allowance;
    if !positive(allowance) || !value.is_finite() {
        Err(Error::Invalid(
            "finite-wall adaptive error is unrepresentable",
        ))
    } else {
        Ok(value)
    }
}

pub(in super::super) fn gas_error(
    coarse: &BTreeMap<u64, ClosedGasParcel>,
    fine: &BTreeMap<u64, ClosedGasParcel>,
    controls: EffusionAdvanceControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    if !coarse.keys().eq(fine.keys()) {
        return Err(Error::Invalid(
            "finite-wall adaptive gas identities changed",
        ));
    }
    let mut error = 0.0_f64;
    for (&id, a) in coarse {
        poll(cancelled)?;
        let b = &fine[&id];
        let sa = a.state().state();
        let sb = b.state().state();
        error = error.max(ratio(
            sa.temperature,
            sb.temperature,
            sa.temperature.abs().max(sb.temperature.abs()),
            controls.temperature_tolerance_k,
            controls.relative_tolerance,
        )?);
        error = error.max(ratio(
            sa.pressure,
            sb.pressure,
            sa.pressure.abs().max(sb.pressure.abs()),
            controls.pressure_tolerance_pa,
            controls.relative_tolerance,
        )?);
        let ia: BTreeMap<SpeciesId, f64> = a
            .species_inventory_kg()
            .map_err(|_| Error::Invalid("finite-wall gas inventory is unrepresentable"))?
            .into_iter()
            .collect();
        let ib: BTreeMap<SpeciesId, f64> = b
            .species_inventory_kg()
            .map_err(|_| Error::Invalid("finite-wall gas inventory is unrepresentable"))?
            .into_iter()
            .collect();
        let keys: BTreeSet<_> = ia.keys().chain(ib.keys()).cloned().collect();
        for species in keys {
            poll(cancelled)?;
            let x = ia.get(&species).copied().unwrap_or(0.0);
            let y = ib.get(&species).copied().unwrap_or(0.0);
            error = error.max(ratio(
                x,
                y,
                x.abs().max(y.abs()),
                controls.species_mass_tolerance_kg,
                controls.relative_tolerance,
            )?);
        }
    }
    poll(cancelled)?;
    Ok(error)
}

pub(in super::super) fn wall_error<'l, 'p>(
    coarse: &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    fine: &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    controls: EffusionAdvanceControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    if !coarse.keys().eq(fine.keys()) {
        return Err(Error::Invalid(
            "finite-wall adaptive wall identities changed",
        ));
    }
    let mut error = 0.0_f64;
    for (&id, a) in coarse {
        poll(cancelled)?;
        let b = &fine[&id];
        if a.volume_m3().to_bits() != b.volume_m3().to_bits()
            || a.mass_kg().to_bits() != b.mass_kg().to_bits()
            || a.law().identity() != b.law().identity()
        {
            return Err(Error::Invalid(
                "finite-wall adaptive material identity changed",
            ));
        }
        error = error.max(ratio(
            a.temperature_k(),
            b.temperature_k(),
            a.temperature_k().abs().max(b.temperature_k().abs()),
            controls.temperature_tolerance_k,
            controls.relative_tolerance,
        )?);
        let absolute_specific = controls.stream_energy_tolerance_j / a.mass_kg();
        // Use the original caloric integral and a source-domain anchor, so
        // changing the arbitrary wall energy datum cannot loosen accuracy.
        let law = a.law();
        let lower = law.path().temperature_span_k()[0];
        let relative_storage = |from, to| {
            law.path()
                .storage_change_j_m3(from, to)
                .map(|value| value / law.reference_density_kg_m3())
                .map_err(|_| Error::Invalid("finite-wall error left its original caloric support"))
        };
        let scale = relative_storage(lower, a.temperature_k())?
            .abs()
            .max(relative_storage(lower, b.temperature_k())?.abs());
        error = error.max(ratio(
            0.0,
            relative_storage(a.temperature_k(), b.temperature_k())?,
            scale,
            absolute_specific,
            controls.relative_tolerance,
        )?);
        let _ = id;
    }
    poll(cancelled)?;
    Ok(error)
}

fn audit_window<'l, 'p>(
    before_gas: &BTreeMap<u64, ClosedGasParcel>,
    after_gas: &BTreeMap<u64, ClosedGasParcel>,
    before_walls: &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    after_walls: &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    controls: GasTransportSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<FiniteChannelWindowBalance, Error> {
    poll(cancelled)?;
    if !before_walls.keys().eq(after_walls.keys()) {
        return Err(Error::Invalid("finite-wall window lost wall identities"));
    }
    let gas = entropy::closed_inventory_change(
        before_gas
            .iter()
            .map(|(&id, old)| (id, old, &after_gas[&id])),
        controls.exchange.mass_tolerance_kg,
        &mut *cancelled,
    )
    .map_err(|error| match error {
        entropy::Error::Cancelled => Error::Cancelled,
        entropy::Error::Invalid(_) | entropy::Error::Parcel { .. } => {
            Error::Invalid("finite-wall window gas audit failed")
        }
    })?;
    let mut wall_energy = Dd::ZERO;
    let mut wall_entropy = Dd::ZERO;
    for (&id, old) in before_walls {
        poll(cancelled)?;
        let next = &after_walls[&id];
        if old.volume_m3().to_bits() != next.volume_m3().to_bits()
            || old.mass_kg().to_bits() != next.mass_kg().to_bits()
            || old.law().identity() != next.law().identity()
        {
            return Err(Error::Invalid(
                "finite-wall window changed a material inventory",
            ));
        }
        wall_energy = wall_energy
            + (Dd::from_f64(next.specific_energy_j_kg())
                - Dd::from_f64(old.specific_energy_j_kg()))
                * Dd::from_f64(old.mass_kg());
        let ds = old
            .law()
            .specific_entropy_change_j_kg_k(old.specific_energy_j_kg(), next.specific_energy_j_kg())
            .map_err(|_| Error::Invalid("finite-wall entropy audit left its material domain"))?
            .ok_or(Error::Invalid("finite-wall entropy law became unavailable"))?;
        wall_entropy = wall_entropy + Dd::from_f64(ds) * Dd::from_f64(old.mass_kg());
    }
    let energy_residual = (Dd::from_f64(gas.energy_change_j) + wall_energy).to_f64();
    let total_entropy = (Dd::from_f64(gas.entropy_change_j_k) + wall_entropy).to_f64();
    if !energy_residual.is_finite()
        || !total_entropy.is_finite()
        || gas.total_mass_residual_kg.abs() > controls.exchange.mass_tolerance_kg
        || gas.maximum_species_mass_residual_kg > controls.exchange.mass_tolerance_kg
        || energy_residual.abs() > controls.exchange.energy_tolerance_j
        || total_entropy < -controls.entropy_tolerance_j_k
    {
        return Err(Error::Invalid(
            "finite-wall window exceeded its fixed physical budget",
        ));
    }
    poll(cancelled)?;
    Ok(FiniteChannelWindowBalance {
        total_mass_residual_kg: gas.total_mass_residual_kg,
        maximum_species_mass_residual_kg: gas.maximum_species_mass_residual_kg,
        gas_energy_change_j: gas.energy_change_j,
        wall_energy_change_j: wall_energy.to_f64(),
        energy_residual_j: energy_residual,
        gas_entropy_change_j_k: gas.entropy_change_j_k,
        wall_entropy_change_j_k: wall_entropy.to_f64(),
        total_entropy_change_j_k: total_entropy,
    })
}

/// Adapt gas plus finite material walls over one requested physical interval.
///
/// A successful result publishes only accepted resampled half steps. Every trial
/// invokes the original finite-wall transaction, so shared-wall pooling, wall
/// source inversion, gas reconstruction, endpoint caloric admission and all
/// gross withdrawal limits remain authoritative at each actual state.
#[allow(clippy::too_many_arguments)]
pub fn advance_candidate<'l, 'p, C: FnMut() -> bool>(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    walls: &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    paths: GasTransportPaths<'_>,
    links: &[FiniteChannelNetworkLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    step_controls: GasTransportSolve,
    controls: EffusionAdvanceControls,
    mut cancelled: C,
) -> Result<FiniteChannelAdvanceCandidate<'l, 'p>, Error> {
    let history = advance_history(
        parcels,
        walls,
        duration_s,
        controls,
        super::super::super::trajectory::MAX_EFFUSION_RETAINED_ROWS,
        |gas, wall, dt, cancelled| {
            step_candidate(
                gas,
                wall,
                paths,
                links,
                dt,
                endpoints,
                step_controls,
                cancelled,
            )
            .map(|candidate| candidate.into_parts())
        },
        |_, _, _, _| Ok(0.0),
        |step| {
            step.walls.len()
                + step.exchange.fluxes.len()
                + step.exchange.heat.len()
                + step.exchange.external_heat.len()
                + step.exchange.boundaries.len()
                + step.exchange.parcels.len()
                + step
                    .channels
                    .iter()
                    .map(|row| 1 + row.transmitted.species.len() + row.wall.reflected_species.len())
                    .sum::<usize>()
        },
        &mut cancelled,
    )?;
    let balance = audit_window(
        parcels,
        &history.parcels,
        walls,
        &history.walls,
        step_controls,
        &mut cancelled,
    )?;
    poll(&mut cancelled)?;
    Ok(FiniteChannelAdvanceCandidate {
        parcels: history.parcels,
        walls: history.walls,
        step: FiniteChannelAdvanceStep {
            duration_s,
            attempted_steps: history.attempted_steps,
            accepted_intervals: history.accepted_intervals,
            rejected_intervals: history.rejected_intervals,
            maximum_accepted_error: history.maximum_accepted_error,
            accepted_steps: history.steps,
            balance,
        },
    })
}

pub(in super::super) type ThermalStepParts<'l, 'p, S> = (
    BTreeMap<u64, ClosedGasParcel>,
    BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    S,
);

/// Shared time integration for closed channels and finite-wall vacuum outlets.
/// The physical step and complete-window audit retain their original owners.
pub(in super::super) struct ThermalHistory<'l, 'p, S> {
    pub parcels: BTreeMap<u64, ClosedGasParcel>,
    pub walls: BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    pub steps: Vec<S>,
    pub attempted_steps: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub maximum_accepted_error: f64,
}

#[allow(clippy::too_many_arguments)]
pub(in super::super) fn advance_history<'l, 'p, S, C: FnMut() -> bool>(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    walls: &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    duration_s: f64,
    controls: EffusionAdvanceControls,
    maximum_retained_rows: usize,
    mut evaluate: impl FnMut(
        &BTreeMap<u64, ClosedGasParcel>,
        &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
        f64,
        &mut C,
    ) -> Result<ThermalStepParts<'l, 'p, S>, FiniteChannelError>,
    mut stream_error: impl FnMut(&S, &S, &S, &mut C) -> Result<f64, Error>,
    retained_rows: impl Fn(&S) -> usize,
    cancelled: &mut C,
) -> Result<ThermalHistory<'l, 'p, S>, Error> {
    poll(cancelled)?;
    controls
        .validate()
        .map_err(|_| Error::Invalid("invalid finite-wall adaptive controls"))?;
    if !positive(duration_s) || parcels.is_empty() || walls.is_empty() || maximum_retained_rows == 0
    {
        return Err(Error::Invalid(
            "finite-wall advance needs positive duration and inventories",
        ));
    }
    let mut gas = parcels.clone();
    let mut wall = walls.clone();
    let mut elapsed = 0.0;
    let mut ceiling = controls.initial_step_s;
    let (mut attempted, mut accepted, mut rejected) = (0usize, 0usize, 0usize);
    let mut maximum_error = 0.0_f64;
    let mut history = Vec::new();
    let mut history_rows = 0usize;

    while elapsed < duration_s {
        poll(cancelled)?;
        if accepted == controls.maximum_accepted_intervals || attempted >= controls.maximum_attempts
        {
            return Err(Error::WorkBudget {
                elapsed_s: elapsed,
                attempted_steps: attempted,
                accepted_intervals: accepted,
            });
        }
        let remaining = duration_s - elapsed;
        let span = ceiling.min(remaining);
        let half = 0.5 * span;
        if !positive(span)
            || !positive(half)
            || half + half != span
            || elapsed + half <= elapsed
            || elapsed + span <= elapsed
        {
            return Err(Error::Invalid(
                "finite-wall adaptive clock is unrepresentable",
            ));
        }
        // The cancel callback is a parameter, not a capture: the error estimates
        // below borrow the same callback inside the trial.
        let mut eval = |g: &BTreeMap<u64, ClosedGasParcel>,
                        w: &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
                        dt: f64,
                        at: f64,
                        cancelled: &mut C|
         -> Result<ThermalStepParts<'l, 'p, S>, Error> {
            if attempted == controls.maximum_attempts {
                return Err(Error::WorkBudget {
                    elapsed_s: at,
                    attempted_steps: attempted,
                    accepted_intervals: accepted,
                });
            }
            attempted += 1;
            evaluate(g, w, dt, cancelled).map_err(|e| physical(e, at, dt))
        };
        let trial = (|| {
            let coarse = eval(&gas, &wall, span, elapsed, cancelled)?;
            let first = eval(&gas, &wall, half, elapsed, cancelled)?;
            let fine = eval(&first.0, &first.1, half, elapsed + half, cancelled)?;
            let e = gas_error(&coarse.0, &fine.0, controls, cancelled)?
                .max(wall_error(&coarse.1, &fine.1, controls, cancelled)?)
                .max(stream_error(&coarse.2, &first.2, &fine.2, cancelled)?);
            if !e.is_finite() || e < 0.0 {
                return Err(Error::Invalid("finite-wall adaptive estimate is invalid"));
            }
            Ok::<_, Error>((first, fine, e))
        })();
        match trial {
            Ok((first, fine, error)) if error <= 1.0 => {
                history_rows = history_rows
                    .checked_add(retained_rows(&first.2))
                    .and_then(|n| n.checked_add(retained_rows(&fine.2)))
                    .filter(|&n| n <= maximum_retained_rows)
                    .ok_or(Error::WorkBudget {
                        elapsed_s: elapsed,
                        attempted_steps: attempted,
                        accepted_intervals: accepted,
                    })?;
                let (_, _, first_step) = first;
                let (next_gas, next_wall, fine_step) = fine;
                history.push(first_step);
                history.push(fine_step);
                gas = next_gas;
                wall = next_wall;
                elapsed = if span == remaining {
                    duration_s
                } else {
                    elapsed + span
                };
                accepted += 1;
                maximum_error = maximum_error.max(error);
                let factor = if error == 0.0 {
                    2.0
                } else {
                    (0.9 / fs_math::det::sqrt(error)).clamp(0.5, 2.0)
                };
                ceiling = (span * factor).clamp(controls.minimum_step_s, controls.maximum_step_s);
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

    poll(cancelled)?;
    Ok(ThermalHistory {
        parcels: gas,
        walls: wall,
        attempted_steps: attempted,
        accepted_intervals: accepted,
        rejected_intervals: rejected,
        maximum_accepted_error: maximum_error,
        steps: history,
    })
}

#[cfg(test)]
mod tests;
