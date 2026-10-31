//! Complete molecular-diffusion windows with local temporal error control.
//!
//! One full explicit step estimates the error of two half steps; only the fine
//! path can be accepted. Each stage uses the original simultaneous face/source,
//! species inventory, caloric reconstruction and entropy owners. This remains
//! first order in time. The local estimate is not a global error certificate or
//! a license to leave the original isothermal/isobaric constitutive domain.

use super::super::super::super::snapshot::{MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS};
use super::super::super::{SpeciesId, Sum, inventory};
use super::{
    ClosedGasParcel, DiffusionNetworkCandidate, DiffusionNetworkError, DiffusionNetworkLink,
    DiffusionNetworkSolve, DiffusionNetworkStep, GasEnergySolve, SpeciesNetworkError,
    candidate_before_entropy_gate, check, entropy_change,
};
use std::collections::BTreeMap;

/// Local accuracy and finite work; these never replace physical exchange gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiffusionNetworkAdvanceControls {
    /// Largest coarse comparison interval [s]. Accepted states use two half steps.
    pub maximum_step_s: f64,
    /// Smallest retry interval [s]. A shorter final remainder may be attempted.
    pub minimum_step_s: f64,
    /// Complete snapshot evaluations, including rejected coarse and fine probes.
    /// At least three are required; the hard bound is 100,000.
    pub maximum_evaluations: usize,
    /// Per-parcel, per-species absolute inventory error allowance [kg].
    pub absolute_species_tolerance_kg: f64,
    /// Relative allowance times the largest original/coarse/fine species mass.
    /// Zero selects the absolute allowance alone.
    pub relative_species_tolerance: f64,
    /// Absolute difference between original-source endpoint temperatures [K].
    pub temperature_tolerance_k: f64,
}

/// Accepted gross transport only, accumulated in canonical physical-edge order.
/// Reverse transfer has its own row; arrivals may move onward in later substeps.
#[derive(Clone, Debug, PartialEq)]
pub struct DiffusionNetworkIntegratedTransfer {
    pub donor: u64,
    pub receiver: u64,
    pub species: SpeciesId,
    pub mass_kg: f64,
    /// Actual upstream partial enthalpy, including its declared reference [J].
    /// This gross account is not used to subtract large formation energies.
    pub carried_enthalpy_j: f64,
}

/// Bounded whole-window report. Rejected probes spend work but no inventory.
#[derive(Clone, Debug, PartialEq)]
pub struct DiffusionNetworkAdvanceReport {
    pub duration_s: f64,
    pub accepted_intervals: usize,
    /// Two conservative half steps per accepted comparison interval.
    pub accepted_substeps: usize,
    pub rejected_intervals: usize,
    pub evaluations: usize,
    pub minimum_accepted_substep_s: f64,
    pub maximum_accepted_substep_s: f64,
    pub maximum_species_error_ratio: f64,
    pub maximum_temperature_error_ratio: f64,
    pub transfers: Vec<DiffusionNetworkIntegratedTransfer>,
    /// Actual final-minus-original total mass [kg].
    pub total_mass_residual_kg: f64,
    /// Largest actual final-minus-original species mass over the closed graph.
    pub maximum_species_mass_residual_kg: f64,
    /// Compensated sum of accepted original-caloric residuals after reference
    /// cancellation, as in SpeciesNetworkStep. Inventory roundoff is separate.
    pub total_energy_residual_j: f64,
    /// Original-to-final source-qualified closed-network entropy change [J/K].
    pub entropy_change_j_k: f64,
    /// Source receipts for the last accepted half step; history is not retained.
    pub last_substep: DiffusionNetworkStep,
}

/// Complete proposed network and accounting, published together after all gates.
#[derive(Clone, Debug, PartialEq)]
pub struct DiffusionNetworkAdvanceCandidate {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    report: DiffusionNetworkAdvanceReport,
}
impl DiffusionNetworkAdvanceCandidate {
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub const fn report(&self) -> &DiffusionNetworkAdvanceReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        BTreeMap<u64, ClosedGasParcel>,
        DiffusionNetworkAdvanceReport,
    ) {
        (self.parcels, self.report)
    }
}

/// No error returns a partially advanced network.
#[derive(Clone, Debug, PartialEq)]
pub enum DiffusionNetworkAdvanceError {
    Invalid(&'static str),
    Step {
        elapsed_s: f64,
        step_s: f64,
        error: DiffusionNetworkError,
    },
    Accuracy {
        elapsed_s: f64,
        step_s: f64,
        species_error_ratio: f64,
        temperature_error_ratio: f64,
    },
    Entropy {
        elapsed_s: f64,
        step_s: f64,
        change_j_k: f64,
    },
    WorkLimit {
        elapsed_s: f64,
        accepted_intervals: usize,
        evaluations: usize,
    },
    Cancelled,
}
impl core::fmt::Display for DiffusionNetworkAdvanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive molecular network: {self:?}")
    }
}
impl std::error::Error for DiffusionNetworkAdvanceError {}
use DiffusionNetworkAdvanceError as Error;

enum TrialFailure {
    Physical(DiffusionNetworkError),
    Entropy(f64),
    Work,
}

#[allow(clippy::too_many_arguments)]
fn trial(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    links: &[DiffusionNetworkLink<'_>],
    dt: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    physical: DiffusionNetworkSolve,
    evaluations: &mut usize,
    maximum_evaluations: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<DiffusionNetworkCandidate, TrialFailure> {
    check(cancelled).map_err(TrialFailure::Physical)?;
    if *evaluations == maximum_evaluations {
        return Err(TrialFailure::Work);
    }
    *evaluations += 1;
    let candidate =
        candidate_before_entropy_gate(parcels, links, dt, endpoints, physical, cancelled)
            .map_err(TrialFailure::Physical)?;
    if candidate.step().entropy_change_j_k < -physical.entropy_tolerance_j_k {
        return Err(TrialFailure::Entropy(candidate.step().entropy_change_j_k));
    }
    Ok(candidate)
}

fn physical_error(error: DiffusionNetworkError, elapsed_s: f64, step_s: f64) -> Error {
    if error == DiffusionNetworkError::Cancelled {
        Error::Cancelled
    } else {
        Error::Step {
            elapsed_s,
            step_s,
            error,
        }
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn arithmetic(_: super::GasEnergyError) -> Error {
    Error::Invalid("unrepresentable whole-window diffusion account")
}

fn error_ratios(
    initial: &BTreeMap<u64, ClosedGasParcel>,
    coarse: &BTreeMap<u64, ClosedGasParcel>,
    fine: &BTreeMap<u64, ClosedGasParcel>,
    controls: DiffusionNetworkAdvanceControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, f64), Error> {
    let (mut species, mut temperature) = (0.0_f64, 0.0_f64);
    for (&id, old) in initial {
        poll(cancelled)?;
        let before = inventory(old.model(), old.mass_kg()).map_err(arithmetic)?;
        let a = inventory(coarse[&id].model(), coarse[&id].mass_kg()).map_err(arithmetic)?;
        let b = inventory(fine[&id].model(), fine[&id].mass_kg()).map_err(arithmetic)?;
        // The union retains species that were absent before the first stage or
        // have not yet reached the coarse endpoint. No trace is seeded/clipped.
        for name in before.keys().chain(a.keys()).chain(b.keys()) {
            let (m0, ma, mb) = (
                before.get(name).copied().unwrap_or(0.0),
                a.get(name).copied().unwrap_or(0.0),
                b.get(name).copied().unwrap_or(0.0),
            );
            let allowance = controls.absolute_species_tolerance_kg
                + controls.relative_species_tolerance * m0.max(ma).max(mb);
            let ratio = (mb - ma).abs() / allowance;
            if !allowance.is_finite() || !ratio.is_finite() {
                return Err(Error::Invalid(
                    "unrepresentable species temporal error scale",
                ));
            }
            species = species.max(ratio);
        }
        let ratio =
            (coarse[&id].state().state().temperature - fine[&id].state().state().temperature).abs()
                / controls.temperature_tolerance_k;
        if !ratio.is_finite() {
            return Err(Error::Invalid("unrepresentable temperature temporal error"));
        }
        temperature = temperature.max(ratio);
    }
    Ok((species, temperature))
}

#[derive(Default)]
struct Ledger {
    energy: Sum,
    transfers: BTreeMap<(u64, u64, SpeciesId), (Sum, Sum)>,
}
impl Ledger {
    fn record(
        &mut self,
        step: &DiffusionNetworkStep,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        self.energy
            .add(step.exchange.total_energy_residual_j)
            .map_err(arithmetic)?;
        for row in &step.exchange.fluxes {
            poll(cancelled)?;
            let key = (
                row.transfer.donor,
                row.transfer.receiver,
                row.transfer.species.clone(),
            );
            if !self.transfers.contains_key(&key) && self.transfers.len() == MAX_SNAPSHOT_TRANSFERS
            {
                return Err(Error::Invalid(
                    "whole-window diffusion transfer report exceeds its work bound",
                ));
            }
            let totals = self.transfers.entry(key).or_default();
            totals.0.add(row.transfer.mass_kg).map_err(arithmetic)?;
            totals.1.add(row.carried_enthalpy_j).map_err(arithmetic)?;
        }
        Ok(())
    }
}

fn mass_balance(
    initial: &BTreeMap<u64, ClosedGasParcel>,
    final_state: &BTreeMap<u64, ClosedGasParcel>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, f64), Error> {
    let mut total = Sum::default();
    let mut species: BTreeMap<SpeciesId, Sum> = BTreeMap::new();
    for (&id, old) in initial {
        poll(cancelled)?;
        let next = &final_state[&id];
        total
            .add(next.mass_kg() - old.mass_kg())
            .map_err(arithmetic)?;
        let a = inventory(old.model(), old.mass_kg()).map_err(arithmetic)?;
        let b = inventory(next.model(), next.mass_kg()).map_err(arithmetic)?;
        for (name, mass) in &a {
            species
                .entry(name.clone())
                .or_default()
                .add(b.get(name).copied().unwrap_or(0.0) - mass)
                .map_err(arithmetic)?;
        }
        for (name, mass) in &b {
            if !a.contains_key(name) {
                species
                    .entry(name.clone())
                    .or_default()
                    .add(*mass)
                    .map_err(arithmetic)?;
            }
        }
    }
    Ok((
        total.value,
        species.values().map(|s| s.value.abs()).fold(0.0, f64::max),
    ))
}

/// Advance an arbitrary branching/cyclic molecular network through exactly the
/// requested positive duration. Inputs remain unchanged even after late failure.
///
/// Every attempted full/half step samples the actual original state and material
/// laws. Only typed withdrawal limits, local error excess and endpoint entropy
/// decrease can shorten an interval; source/domain failures propagate directly.
/// Gross withdrawal gates apply separately to each actual accepted half step;
/// arrivals may then diffuse onward during later physical time. A whole-window
/// entropy and actual-inventory audit plus the compensated energy ledger must fit
/// the SAME absolute physical budgets, independent of the number of steps.
///
/// The full-minus-two-half difference estimates the fine first-order local error
/// (denominator 2^1-1=1). No Richardson state extrapolation, clipping, artificial
/// transport coefficients, pressure work or Soret/Dufour terms are introduced.
#[allow(clippy::too_many_arguments)]
pub fn advance_candidate(
    original: &BTreeMap<u64, ClosedGasParcel>,
    links: &[DiffusionNetworkLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    physical: DiffusionNetworkSolve,
    controls: DiffusionNetworkAdvanceControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<DiffusionNetworkAdvanceCandidate, Error> {
    poll(&mut cancelled)?;
    if ![
        duration_s,
        controls.minimum_step_s,
        controls.maximum_step_s,
        controls.absolute_species_tolerance_kg,
        controls.temperature_tolerance_k,
    ]
    .iter()
    .all(|x| super::positive(*x))
        || controls.minimum_step_s > controls.maximum_step_s
        || !(3..=100_000).contains(&controls.maximum_evaluations)
        || !controls.relative_species_tolerance.is_finite()
        || !(0.0..=1.0).contains(&controls.relative_species_tolerance)
        || original.is_empty()
        || original.len() > MAX_SNAPSHOT_PARCELS
        || links.len() > MAX_SNAPSHOT_TRANSFERS
        || !original.keys().eq(endpoints.keys())
    {
        return Err(Error::Invalid(
            "bounded nonempty graph, positive adaptive controls and exact endpoint map required",
        ));
    }
    for (&id, solve) in endpoints {
        poll(&mut cancelled)?;
        solve.validate().map_err(|error| {
            physical_error(DiffusionNetworkError::Parcel { id, error }, 0.0, 0.0)
        })?;
        if !(solve.minimum_temperature_k..=solve.maximum_temperature_k)
            .contains(&original[&id].state().state().temperature)
        {
            return Err(Error::Invalid(
                "original gas temperature outside endpoint caloric bracket",
            ));
        }
    }
    let mut state = original.clone();
    let mut elapsed = 0.0;
    let mut ceiling = controls.maximum_step_s;
    let (mut evaluations, mut accepted, mut rejected) = (0usize, 0usize, 0usize);
    let (mut minimum, mut maximum) = (f64::INFINITY, 0.0_f64);
    let (mut max_species_error, mut max_temperature_error) = (0.0_f64, 0.0_f64);
    let mut ledger = Ledger::default();
    let mut last = None;
    while elapsed < duration_s {
        poll(&mut cancelled)?;
        let remaining = duration_s - elapsed;
        let dt = remaining.min(ceiling);
        let half = 0.5 * dt;
        if !super::positive(half) || half + half != dt || elapsed + dt <= elapsed {
            return Err(Error::Invalid(
                "diffusion interval cannot advance representable physical time",
            ));
        }
        let attempted = (|| {
            let coarse = trial(
                &state,
                links,
                dt,
                endpoints,
                physical,
                &mut evaluations,
                controls.maximum_evaluations,
                &mut cancelled,
            )?;
            let first = trial(
                &state,
                links,
                half,
                endpoints,
                physical,
                &mut evaluations,
                controls.maximum_evaluations,
                &mut cancelled,
            )?;
            let second = trial(
                first.parcels(),
                links,
                half,
                endpoints,
                physical,
                &mut evaluations,
                controls.maximum_evaluations,
                &mut cancelled,
            )?;
            Ok((coarse, first, second))
        })();
        let refusal = match attempted {
            Err(TrialFailure::Work) => {
                return Err(Error::WorkLimit {
                    elapsed_s: elapsed,
                    accepted_intervals: accepted,
                    evaluations,
                });
            }
            Err(TrialFailure::Physical(error)) => {
                let retry = matches!(
                    error,
                    DiffusionNetworkError::Exchange(SpeciesNetworkError::Withdrawal { .. })
                );
                let error = physical_error(error, elapsed, dt);
                if !retry {
                    return Err(error);
                }
                error
            }
            Err(TrialFailure::Entropy(change_j_k)) => Error::Entropy {
                elapsed_s: elapsed,
                step_s: dt,
                change_j_k,
            },
            Ok((coarse, first, second)) => {
                let (species_error, temperature_error) = error_ratios(
                    &state,
                    coarse.parcels(),
                    second.parcels(),
                    controls,
                    &mut cancelled,
                )?;
                if species_error <= 1.0 && temperature_error <= 1.0 {
                    ledger.record(first.step(), &mut cancelled)?;
                    ledger.record(second.step(), &mut cancelled)?;
                    let (next, step) = second.into_parts();
                    state = next;
                    last = Some(step);
                    accepted += 1;
                    minimum = minimum.min(half);
                    maximum = maximum.max(half);
                    max_species_error = max_species_error.max(species_error);
                    max_temperature_error = max_temperature_error.max(temperature_error);
                    elapsed = if dt == remaining {
                        duration_s
                    } else {
                        elapsed + dt
                    };
                    // First-order local error scales quadratically with dt.
                    ceiling = if species_error.max(temperature_error) <= 0.25 {
                        if dt < 0.5 * controls.maximum_step_s {
                            2.0 * dt
                        } else {
                            controls.maximum_step_s
                        }
                    } else {
                        dt
                    };
                    continue;
                }
                Error::Accuracy {
                    elapsed_s: elapsed,
                    step_s: dt,
                    species_error_ratio: species_error,
                    temperature_error_ratio: temperature_error,
                }
            }
        };
        rejected += 1;
        if half < controls.minimum_step_s || elapsed + half <= elapsed {
            return Err(refusal);
        }
        ceiling = half;
    }
    let (mass, species_mass) = mass_balance(original, &state, &mut cancelled)?;
    let entropy = entropy_change(
        original,
        &state,
        physical.exchange.mass_tolerance_kg,
        &mut cancelled,
    )
    .map_err(|e| physical_error(e, elapsed, 0.0))?;
    if mass.abs() > physical.exchange.mass_tolerance_kg
        || species_mass > physical.exchange.mass_tolerance_kg
        || ledger.energy.value.abs() > physical.exchange.energy_tolerance_j
    {
        return Err(Error::Invalid(
            "complete diffusion window exceeds its species/mass/energy budget",
        ));
    }
    if entropy < -physical.entropy_tolerance_j_k {
        return Err(Error::Entropy {
            elapsed_s: elapsed,
            step_s: 0.0,
            change_j_k: entropy,
        });
    }
    poll(&mut cancelled)?;
    let transfers = ledger
        .transfers
        .into_iter()
        .map(
            |((donor, receiver, species), (mass, enthalpy))| DiffusionNetworkIntegratedTransfer {
                donor,
                receiver,
                species,
                mass_kg: mass.value,
                carried_enthalpy_j: enthalpy.value,
            },
        )
        .collect();
    Ok(DiffusionNetworkAdvanceCandidate {
        parcels: state,
        report: DiffusionNetworkAdvanceReport {
            duration_s,
            accepted_intervals: accepted,
            accepted_substeps: 2 * accepted,
            rejected_intervals: rejected,
            evaluations,
            minimum_accepted_substep_s: minimum,
            maximum_accepted_substep_s: maximum,
            maximum_species_error_ratio: max_species_error,
            maximum_temperature_error_ratio: max_temperature_error,
            transfers,
            total_mass_residual_kg: mass,
            maximum_species_mass_residual_kg: species_mass,
            total_energy_residual_j: ledger.energy.value,
            entropy_change_j_k: entropy,
            last_substep: last.ok_or(Error::Invalid("diffusion window accepted no interval"))?,
        },
    })
}
