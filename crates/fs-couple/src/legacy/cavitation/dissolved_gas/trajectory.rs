//! Temperature-dependent absorption histories through the original finite pair.
//!
//! One native full step estimates the local error of two native half steps.
//! Only the two-half state is accepted; no extrapolated species, temperature or
//! material source is constructed. For first-order splitting the fine error is
//! estimated by `(fine - coarse) * (h1²+h2²)/(2 h1 h2)`. Actual representable
//! durations determine that factor, which is one for equal halves.
//!
//! This is an asymptotic local estimate, not a bound on global trajectory error
//! or uncertainty in the material laws. Held external power, fixed films and
//! the original common-temperature dilute model apply throughout one call.
//! Cut a call at changes in forcing or geometry. Energy and species conservation
//! have separate, fixed whole-window budgets, never multiplied by step count.

use super::{
    ClosedGasParcel, DissolutionFilm, DissolvedGasInventory, HenryDissolution,
    IncompressibleSolvent, ThermalDissolutionControls, ThermalDissolutionError,
    ThermalDissolutionStep, advance_thermal_mixture,
};
use fs_material::cavitation::dissolved_gas::mixture::MAX_DISSOLVING_SPECIES;
use fs_material::thermochemical_gas::energy::GasEnergyError;
use fs_math::dd::Dd;
use fs_thermochem::SpeciesId;
use std::collections::BTreeMap;

/// Temporal accuracy and whole-window work/conservation limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalDissolutionTimeControls {
    pub initial_step_s: f64,
    /// Minimum comparison span. Its accepted native steps are approximately
    /// half as long; a smaller final remainder must pass without a retry.
    pub minimum_step_s: f64,
    /// Must not exceed the original mass-transfer owner's maximum interval.
    pub maximum_step_s: f64,
    /// Counts every native call, including rejected/failed trials; 3..=100_000.
    pub maximum_evaluations: usize,
    /// Multiplies each coordinate's CHANGE from its last accepted state.
    /// Large carrier inventories and absolute temperatures cannot mask error.
    pub relative_tolerance: f64,
    pub temperature_tolerance_k: f64,
    pub pressure_tolerance_pa: f64,
    /// Applied separately to each gas and dissolved species, without averaging.
    pub species_tolerance_mol: f64,
    /// Fixed complete-window energy budget, independent of local root budgets.
    pub energy_tolerance_j: f64,
    /// Fixed original-to-final per-species AND total mass budget [kg].
    pub mass_tolerance_kg: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ThermalDissolutionTrajectoryError {
    Invalid(&'static str),
    Inventory(GasEnergyError),
    Step {
        time_s: f64,
        step_s: f64,
        error: ThermalDissolutionError,
    },
    Accuracy {
        time_s: f64,
        step_s: f64,
        error_ratio: f64,
    },
    WorkLimit {
        time_s: f64,
        evaluations: usize,
        accepted_intervals: usize,
    },
    TimeResolution {
        time_s: f64,
        target_s: f64,
    },
    Balance {
        quantity: &'static str,
        residual: f64,
        tolerance: f64,
    },
    Cancelled,
}
impl core::fmt::Display for ThermalDissolutionTrajectoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermal absorption trajectory: {self:?}")
    }
}
impl std::error::Error for ThermalDissolutionTrajectoryError {}
type Error = ThermalDissolutionTrajectoryError;

/// Accepted gross transfers retain reversals that a net endpoint cannot show.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalDissolutionSpeciesHistory {
    pub species: SpeciesId,
    pub absorbed_moles: f64,
    pub exsolved_moles: f64,
    pub dissolved_moles_change: f64,
    pub mass_residual_kg: f64,
}

/// Bounded history: canonical species totals and the last original receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalDissolutionTrajectory {
    pub interval_s: [f64; 2],
    pub evaluations: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub minimum_accepted_step_s: f64,
    pub maximum_accepted_step_s: f64,
    pub maximum_accepted_error_ratio: f64,
    pub suggested_step_s: f64,
    pub species: Vec<ThermalDissolutionSpeciesHistory>,
    /// Native transfer storage increments, including their actual residuals.
    pub isothermal_transfer_energy_change_j: f64,
    pub gas_sensible_energy_change_j: f64,
    pub dissolved_sensible_energy_change_j: f64,
    pub solvent_energy_change_j: f64,
    /// Prescribed held power times the complete requested duration.
    pub external_heat_j: f64,
    /// Accepted native heat inputs only; rejected probes never enter accounts.
    pub integrated_external_heat_j: f64,
    pub energy_residual_j: f64,
    pub mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    pub last_substep: ThermalDissolutionStep,
}

#[derive(Clone)]
struct State {
    laws: Vec<HenryDissolution>,
    gas: ClosedGasParcel,
    liquids: Vec<DissolvedGasInventory>,
}
#[derive(Clone)]
struct Snapshot {
    temperature: f64,
    pressure: f64,
    gas_mass: BTreeMap<SpeciesId, f64>,
    dissolved_moles: Vec<f64>,
}
impl Snapshot {
    fn new(state: &State, order: &[usize]) -> Result<Self, Error> {
        Ok(Self {
            temperature: state.gas.state().state().temperature,
            pressure: state.gas.state().state().pressure,
            gas_mass: state.gas.species_inventory_kg().map_err(Error::Inventory)?,
            dissolved_moles: order
                .iter()
                .map(|&i| state.liquids[i].dissolved_moles())
                .collect(),
        })
    }
    fn mass(&self, species: &SpeciesId) -> f64 {
        self.gas_mass.get(species).copied().unwrap_or(0.0)
    }
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn finite(x: f64) -> Result<f64, Error> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(Error::Invalid("nonfinite trajectory arithmetic"))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn physical(error: ThermalDissolutionError, time_s: f64, step_s: f64) -> Error {
    if error == ThermalDissolutionError::Cancelled {
        Error::Cancelled
    } else {
        Error::Step {
            time_s,
            step_s,
            error,
        }
    }
}
fn ratio(old: f64, coarse: f64, fine: f64, absolute: f64, relative: f64) -> Result<f64, Error> {
    let scale = (coarse - old).abs().max((fine - old).abs());
    let allowance = finite(relative.mul_add(scale, absolute))?;
    finite((fine - coarse).abs() / allowance)
}
fn error_ratio(
    old: &Snapshot,
    coarse: &Snapshot,
    fine: &Snapshot,
    axis: &[(SpeciesId, f64)],
    controls: ThermalDissolutionTimeControls,
    factor: f64,
) -> Result<f64, Error> {
    let mut error = ratio(
        old.temperature,
        coarse.temperature,
        fine.temperature,
        controls.temperature_tolerance_k,
        controls.relative_tolerance,
    )?;
    error = error.max(ratio(
        old.pressure,
        coarse.pressure,
        fine.pressure,
        controls.pressure_tolerance_pa,
        controls.relative_tolerance,
    )?);
    for (species, molar_mass) in axis {
        error = error.max(ratio(
            old.mass(species) / molar_mass,
            coarse.mass(species) / molar_mass,
            fine.mass(species) / molar_mass,
            controls.species_tolerance_mol,
            controls.relative_tolerance,
        )?);
    }
    for i in 0..old.dissolved_moles.len() {
        error = error.max(ratio(
            old.dissolved_moles[i],
            coarse.dissolved_moles[i],
            fine.dissolved_moles[i],
            controls.species_tolerance_mol,
            controls.relative_tolerance,
        )?);
    }
    finite(factor * error)
}

/// Resolve a held-power thermal absorption window to explicit local tolerances.
///
/// Sources, finite inventories and heat always evolve through the original
/// conservative endpoint. Only temporal-error or temperature-increment failures
/// permit shorter retries; source, regime, caloric and conservation failures
/// remain physical refusals. Every unsuccessful call preserves ALL input state,
/// even after an internally accepted prefix. Memory is independent of duration.
#[allow(clippy::too_many_arguments)]
pub fn advance_thermal_to(
    laws: &mut [HenryDissolution],
    gas: &mut ClosedGasParcel,
    liquids: &mut [DissolvedGasInventory],
    films: &[DissolutionFilm],
    solvent: &IncompressibleSolvent,
    end_time_s: f64,
    external_power_w: f64,
    step: ThermalDissolutionControls,
    time: ThermalDissolutionTimeControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ThermalDissolutionTrajectory, Error> {
    poll(&mut cancelled)?;
    if laws.is_empty()
        || laws.len() > MAX_DISSOLVING_SPECIES
        || laws.len() != liquids.len()
        || laws.len() != films.len()
    {
        return Err(Error::Invalid(
            "matching nonempty absorption layouts required",
        ));
    }
    let start = liquids[0].time_s();
    let duration = end_time_s - start;
    if !positive(duration)
        || !end_time_s.is_finite()
        || !external_power_w.is_finite()
        || ![
            time.initial_step_s,
            time.minimum_step_s,
            time.maximum_step_s,
            time.temperature_tolerance_k,
            time.pressure_tolerance_pa,
            time.species_tolerance_mol,
            time.energy_tolerance_j,
            time.mass_tolerance_kg,
        ]
        .into_iter()
        .all(positive)
        || time.minimum_step_s > time.initial_step_s
        || time.initial_step_s > time.maximum_step_s
        || time.maximum_step_s > step.exchange.maximum_step_s
        || !(3..=100_000).contains(&time.maximum_evaluations)
        || !time.relative_tolerance.is_finite()
        || !(0.0..1.0).contains(&time.relative_tolerance)
    {
        return Err(Error::Invalid(
            "forward finite window, bounded time/work and positive dimensional tolerances required",
        ));
    }
    let external_heat = finite(external_power_w * duration)?;
    let mut order: Vec<_> = (0..laws.len()).collect();
    order.sort_by(|&a, &b| laws[a].species().cmp(laws[b].species()));
    let mut axis: Vec<_> = laws[0]
        .catalog()
        .thermodynamics()
        .components()
        .iter()
        .map(|s| (s.species().clone(), s.molar_mass().value()))
        .collect();
    axis.sort_by(|a, b| a.0.cmp(&b.0));
    let mut state = State {
        laws: laws.to_vec(),
        gas: gas.clone(),
        liquids: liquids.to_vec(),
    };
    let initial = Snapshot::new(&state, &order)?;
    let mut previous = initial.clone();
    let (mut current, mut proposed) = (start, time.initial_step_s);
    let (mut evaluations, mut accepted, mut rejected) = (0, 0, 0);
    let (mut minimum_step, mut maximum_step, mut maximum_error) = (f64::INFINITY, 0.0_f64, 0.0_f64);
    let mut sum = [Dd::ZERO; 5];
    let mut gross = vec![[Dd::ZERO; 2]; laws.len()];
    let mut last = None;
    while current < end_time_s {
        poll(&mut cancelled)?;
        let remaining = end_time_s - current;
        let mut h = proposed.min(remaining);
        let tail = remaining - h;
        if tail > 0.0 && tail < time.minimum_step_s && 0.5 * remaining >= time.minimum_step_s {
            h = 0.5 * remaining;
        }
        let end = if h == remaining {
            end_time_s
        } else {
            current + h
        };
        let span = end - current;
        let middle = current + 0.5 * span;
        let first_h = middle - current;
        let second_h = end - middle;
        if !end.is_finite()
            || end > end_time_s
            || span > time.maximum_step_s
            || middle <= current
            || middle >= end
        {
            return Err(Error::TimeResolution {
                time_s: current,
                target_s: end,
            });
        }
        let half_ratio = first_h / second_h;
        let factor = finite(0.5 * (half_ratio + half_ratio.recip()))?;
        let trial = (|| {
            let mut advance = |candidate: &mut State, target| {
                poll(&mut cancelled)?;
                if evaluations == time.maximum_evaluations {
                    return Err(Error::WorkLimit {
                        time_s: current,
                        evaluations,
                        accepted_intervals: accepted,
                    });
                }
                evaluations += 1;
                let at = candidate.liquids[0].time_s();
                let dt = target - at;
                let heat = finite(external_power_w * dt)?;
                advance_thermal_mixture(
                    &mut candidate.laws,
                    &mut candidate.gas,
                    &mut candidate.liquids,
                    films,
                    solvent,
                    target,
                    heat,
                    step,
                    &mut cancelled,
                )
                .map_err(|error| physical(error, at, dt))
            };
            let coarse = {
                let mut candidate = state.clone();
                advance(&mut candidate, end)?;
                Snapshot::new(&candidate, &order)?
            };
            let mut fine = state.clone();
            let first = advance(&mut fine, middle)?;
            let second = advance(&mut fine, end)?;
            let snapshot = Snapshot::new(&fine, &order)?;
            let error = error_ratio(&previous, &coarse, &snapshot, &axis, time, factor)?;
            Ok::<_, Error>((fine, snapshot, [first, second], error))
        })();
        let (fine, snapshot, reports, error) = match trial {
            Ok(trial) => trial,
            Err(Error::Step {
                error: ThermalDissolutionError::TemperatureIncrement { .. },
                ..
            }) if 0.5 * span >= time.minimum_step_s => {
                rejected += 1;
                proposed = 0.5 * span;
                continue;
            }
            Err(error) => return Err(error),
        };
        if error > 1.0 {
            if 0.5 * span < time.minimum_step_s {
                return Err(Error::Accuracy {
                    time_s: current,
                    step_s: span,
                    error_ratio: error,
                });
            }
            rejected += 1;
            proposed = 0.5 * span;
            continue;
        }
        for report in reports {
            poll(&mut cancelled)?;
            let values = [
                report.isothermal_exchange.thermostat_heat_into_system_j,
                report.gas_sensible_energy_change_j,
                report.dissolved_sensible_energy_change_j,
                report.solvent_energy_change_j,
                report.external_heat_j,
            ];
            for (total, value) in sum.iter_mut().zip(values) {
                *total = *total + Dd::from_f64(value);
            }
            sum[0] = sum[0] + Dd::from_f64(report.isothermal_exchange.energy_residual_j);
            for (total, exchange) in gross.iter_mut().zip(&report.isothermal_exchange.exchanges) {
                let dn = exchange.dissolved_moles_change;
                let direction = usize::from(dn < 0.0);
                total[direction] = total[direction] + Dd::from_f64(dn.abs());
            }
            last = Some(report);
        }
        state = fine;
        previous = snapshot;
        current = end;
        accepted += 1;
        minimum_step = minimum_step.min(first_h).min(second_h);
        maximum_step = maximum_step.max(first_h).max(second_h);
        maximum_error = maximum_error.max(error);
        proposed = if error <= 0.25 {
            (2.0 * span).min(time.maximum_step_s)
        } else {
            span
        };
    }

    let selected: BTreeMap<_, _> = order
        .iter()
        .enumerate()
        .map(|(i, &input)| (laws[input].species().clone(), i))
        .collect();
    let mut mass_sum = Dd::ZERO;
    let mut maximum_mass_residual = 0.0_f64;
    let mut species_residuals = vec![0.0; laws.len()];
    for (species, molar_mass) in &axis {
        poll(&mut cancelled)?;
        let mut residual =
            Dd::from_f64(previous.mass(species)) - Dd::from_f64(initial.mass(species));
        if let Some(&i) = selected.get(species) {
            residual = residual
                + (Dd::from_f64(previous.dissolved_moles[i])
                    - Dd::from_f64(initial.dissolved_moles[i]))
                    * Dd::from_f64(*molar_mass);
            species_residuals[i] = finite(residual.to_f64())?;
        }
        maximum_mass_residual = maximum_mass_residual.max(finite(residual.to_f64())?.abs());
        mass_sum = mass_sum + residual;
    }
    let mass_residual = finite(mass_sum.to_f64())?;
    let energy_residual =
        finite((sum[0] + sum[1] + sum[2] + sum[3] - Dd::from_f64(external_heat)).to_f64())?;
    for (quantity, residual, tolerance) in [
        (
            "whole-window species mass",
            maximum_mass_residual,
            time.mass_tolerance_kg,
        ),
        (
            "whole-window total mass",
            mass_residual,
            time.mass_tolerance_kg,
        ),
        (
            "whole-window internal energy",
            energy_residual,
            time.energy_tolerance_j,
        ),
    ] {
        if residual.abs() > tolerance {
            return Err(Error::Balance {
                quantity,
                residual,
                tolerance,
            });
        }
    }
    let species = order
        .iter()
        .enumerate()
        .map(|(i, &input)| {
            Ok(ThermalDissolutionSpeciesHistory {
                species: laws[input].species().clone(),
                absorbed_moles: finite(gross[i][0].to_f64())?,
                exsolved_moles: finite(gross[i][1].to_f64())?,
                dissolved_moles_change: previous.dissolved_moles[i] - initial.dissolved_moles[i],
                mass_residual_kg: species_residuals[i],
            })
        })
        .collect::<Result<_, Error>>()?;
    let report = ThermalDissolutionTrajectory {
        interval_s: [start, end_time_s],
        evaluations,
        accepted_intervals: accepted,
        rejected_intervals: rejected,
        minimum_accepted_step_s: minimum_step,
        maximum_accepted_step_s: maximum_step,
        maximum_accepted_error_ratio: maximum_error,
        suggested_step_s: proposed,
        species,
        isothermal_transfer_energy_change_j: finite(sum[0].to_f64())?,
        gas_sensible_energy_change_j: finite(sum[1].to_f64())?,
        dissolved_sensible_energy_change_j: finite(sum[2].to_f64())?,
        solvent_energy_change_j: finite(sum[3].to_f64())?,
        external_heat_j: external_heat,
        integrated_external_heat_j: finite(sum[4].to_f64())?,
        energy_residual_j: energy_residual,
        mass_residual_kg: mass_residual,
        maximum_species_mass_residual_kg: maximum_mass_residual,
        last_substep: last.expect("positive duration accepted at least one native pair"),
    };
    poll(&mut cancelled)?;
    laws.clone_from_slice(&state.laws);
    liquids.copy_from_slice(&state.liquids);
    *gas = state.gas;
    Ok(report)
}
