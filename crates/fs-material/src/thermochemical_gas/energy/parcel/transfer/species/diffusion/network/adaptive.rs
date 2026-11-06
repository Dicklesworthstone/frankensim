//! Complete molecular-diffusion windows with local temporal error control.
//!
//! One full step estimates the error of two half steps; only the fine path can
//! be accepted. The same controller, physical inventory comparisons, accepted
//! transfer ledger and final audit serve the original explicit snapshot and the
//! exponential binary sweep. Conservative first-order step-doubling scaling is
//! used for both. The local estimate is not a global error certificate or a
//! license to leave the selected molecular closure's constitutive domain.

use super::super::super::super::snapshot::{MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS};
use super::super::super::network::transport::{
    GasTransportError, GasTransportSolve, audit_closed_window,
};
use super::super::super::{SPECIES_WITHDRAWAL_REFUSAL, SpeciesId, Sum, inventory};
use super::super::DIFFUSION_ENTROPY_REFUSAL;
use super::binary::{
    self, BinaryDiffusionNetworkAdvanceCandidate, BinaryDiffusionNetworkAdvanceControls,
    BinaryDiffusionNetworkCandidate, BinaryDiffusionNetworkControls, BinaryDiffusionNetworkError,
    BinaryDiffusionNetworkFaceStep,
};
use super::{
    ClosedGasParcel, DiffusionNetworkCandidate, DiffusionNetworkError, DiffusionNetworkLink,
    DiffusionNetworkSolve, DiffusionNetworkStep, GasDiffusionDriving, GasEnergyError,
    GasEnergySolve, SpeciesNetworkError, candidate_before_entropy_gate, check,
};
use std::collections::BTreeMap;

/// Local accuracy and finite work; these never replace physical exchange gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiffusionNetworkAdvanceControls {
    /// Largest coarse comparison interval [s]. Accepted states use two half steps.
    pub maximum_step_s: f64,
    /// Smallest retry interval [s]. A shorter final remainder may be attempted.
    pub minimum_step_s: f64,
    /// Complete candidate evaluations, including rejected coarse and fine probes.
    /// The binary owner additionally bounds actual pair calls separately.
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
pub struct DiffusionNetworkAdvanceReport<Step = DiffusionNetworkStep> {
    pub duration_s: f64,
    pub accepted_intervals: usize,
    /// Two conservative half steps per accepted comparison interval.
    pub accepted_substeps: usize,
    pub rejected_intervals: usize,
    pub evaluations: usize,
    /// Actual exponential pair calls, including rejected/failed probes.
    /// Zero for the explicit snapshot entry points.
    pub pair_evaluations: usize,
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
    /// Actual original-to-final source-relative internal-energy change [J].
    /// Independently bounded by the same fixed energy tolerance as the ledger.
    pub endpoint_energy_change_j: f64,
    /// Original-to-final source-qualified closed-network entropy change [J/K].
    pub entropy_change_j_k: f64,
    /// Source receipts for the last accepted half step; history is not retained.
    pub last_substep: Step,
}

/// Complete proposed network and accounting, published together after all gates.
#[derive(Clone, Debug, PartialEq)]
pub struct DiffusionNetworkAdvanceCandidate<Step = DiffusionNetworkStep> {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    report: DiffusionNetworkAdvanceReport<Step>,
}
impl<Step> DiffusionNetworkAdvanceCandidate<Step> {
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub const fn report(&self) -> &DiffusionNetworkAdvanceReport<Step> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        BTreeMap<u64, ClosedGasParcel>,
        DiffusionNetworkAdvanceReport<Step>,
    ) {
        (self.parcels, self.report)
    }
}

/// No error returns a partially advanced network.
#[derive(Clone, Debug, PartialEq)]
pub enum DiffusionNetworkAdvanceError {
    Invalid(&'static str),
    /// Direct original-to-final native conservation/thermodynamic refusal.
    Window(GasTransportError),
    Step {
        elapsed_s: f64,
        step_s: f64,
        error: DiffusionNetworkError,
    },
    BinaryStep {
        elapsed_s: f64,
        step_s: f64,
        error: BinaryDiffusionNetworkError,
    },
    PairWorkLimit {
        elapsed_s: f64,
        accepted_intervals: usize,
        pair_evaluations: usize,
        maximum_pair_evaluations: usize,
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
    Binary(BinaryDiffusionNetworkError),
    Entropy(f64),
    Work,
    PairWork,
}

#[derive(Default)]
struct Work {
    evaluations: usize,
    pair_evaluations: usize,
}

// Borrowed transfer views let both original owners share accounting without
// allocating or cloning an intermediate stream for rejected probes.
struct Transfer<'a> {
    donor: u64,
    receiver: u64,
    species: &'a SpeciesId,
    mass_kg: f64,
    enthalpy_j: f64,
}
trait Candidate {
    type Receipt;
    fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel>;
    fn entropy_change(&self) -> f64;
    fn energy_residual(&self) -> f64;
    fn transfers(&self) -> impl Iterator<Item = Transfer<'_>>;
    fn into_parts(self) -> (BTreeMap<u64, ClosedGasParcel>, Self::Receipt);
}
impl Candidate for DiffusionNetworkCandidate {
    type Receipt = DiffusionNetworkStep;
    fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        DiffusionNetworkCandidate::parcels(self)
    }
    fn entropy_change(&self) -> f64 {
        self.step().entropy_change_j_k
    }
    fn energy_residual(&self) -> f64 {
        self.step().exchange.total_energy_residual_j
    }
    fn transfers(&self) -> impl Iterator<Item = Transfer<'_>> {
        self.step().exchange.fluxes.iter().map(|row| Transfer {
            donor: row.transfer.donor,
            receiver: row.transfer.receiver,
            species: &row.transfer.species,
            mass_kg: row.transfer.mass_kg,
            enthalpy_j: row.carried_enthalpy_j,
        })
    }
    fn into_parts(self) -> (BTreeMap<u64, ClosedGasParcel>, Self::Receipt) {
        DiffusionNetworkCandidate::into_parts(self)
    }
}
impl Candidate for BinaryDiffusionNetworkCandidate {
    type Receipt = Vec<BinaryDiffusionNetworkFaceStep>;
    fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        BinaryDiffusionNetworkCandidate::parcels(self)
    }
    fn entropy_change(&self) -> f64 {
        self.report().entropy_change_j_k
    }
    fn energy_residual(&self) -> f64 {
        self.report().total_energy_residual_j
    }
    fn transfers(&self) -> impl Iterator<Item = Transfer<'_>> {
        self.report().transfers.iter().map(|row| Transfer {
            donor: row.donor,
            receiver: row.receiver,
            species: &row.species,
            mass_kg: row.mass_kg,
            enthalpy_j: row.carried_enthalpy_j,
        })
    }
    fn into_parts(self) -> (BTreeMap<u64, ClosedGasParcel>, Self::Receipt) {
        let (parcels, report) = BinaryDiffusionNetworkCandidate::into_parts(self);
        (parcels, report.last_sweep)
    }
}

#[allow(clippy::too_many_arguments)]
fn trial<C: Candidate>(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    dt: f64,
    physical: DiffusionNetworkSolve,
    work: &mut Work,
    maximum_evaluations: usize,
    evaluate: &mut impl FnMut(
        &BTreeMap<u64, ClosedGasParcel>,
        f64,
        &mut Work,
        &mut dyn FnMut() -> bool,
    ) -> Result<C, TrialFailure>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<C, TrialFailure> {
    check(cancelled).map_err(TrialFailure::Physical)?;
    if work.evaluations == maximum_evaluations {
        return Err(TrialFailure::Work);
    }
    work.evaluations += 1;
    let candidate = evaluate(parcels, dt, work, cancelled)?;
    if candidate.entropy_change() < -physical.entropy_tolerance_j_k {
        return Err(TrialFailure::Entropy(candidate.entropy_change()));
    }
    Ok(candidate)
}

fn retryable_binary(error: &BinaryDiffusionNetworkError) -> bool {
    match error {
        BinaryDiffusionNetworkError::Face { error, .. } => {
            *error == SPECIES_WITHDRAWAL_REFUSAL
                || *error == DIFFUSION_ENTROPY_REFUSAL
                || matches!(
                    error,
                    GasEnergyError::OutsideBracket { .. } | GasEnergyError::NotConverged { .. }
                )
        }
        BinaryDiffusionNetworkError::Window(GasTransportError::Diffusion(
            DiffusionNetworkError::EntropyDecrease { .. },
        )) => true,
        _ => false,
    }
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
    fn record<C: Candidate>(
        &mut self,
        step: &C,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        self.energy
            .add(step.energy_residual())
            .map_err(arithmetic)?;
        for row in step.transfers() {
            poll(cancelled)?;
            let key = (row.donor, row.receiver, row.species.clone());
            if !self.transfers.contains_key(&key) && self.transfers.len() == MAX_SNAPSHOT_TRANSFERS
            {
                return Err(Error::Invalid(
                    "whole-window diffusion transfer report exceeds its work bound",
                ));
            }
            let totals = self.transfers.entry(key).or_default();
            totals.0.add(row.mass_kg).map_err(arithmetic)?;
            totals.1.add(row.enthalpy_j).map_err(arithmetic)?;
        }
        Ok(())
    }
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
    cancelled: impl FnMut() -> bool,
) -> Result<DiffusionNetworkAdvanceCandidate, Error> {
    advance_candidate_with_driving(
        original,
        links,
        duration_s,
        endpoints,
        physical,
        controls,
        GasDiffusionDriving::IsobaricMolar,
        cancelled,
    )
}

/// Advance the selected molecular force/frame through the complete window.
/// Each coarse and fine probe resamples pressure driving and qualifies actual
/// source endpoints. Mass-frame pressure diffusion preserves every original
/// species, energy and closed-network entropy gate, including the independent
/// original-to-final audit. No bulk-flow or thermal-diffusion law is introduced.
#[allow(clippy::too_many_arguments)]
pub fn advance_candidate_with_driving(
    original: &BTreeMap<u64, ClosedGasParcel>,
    links: &[DiffusionNetworkLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    physical: DiffusionNetworkSolve,
    controls: DiffusionNetworkAdvanceControls,
    driving: GasDiffusionDriving,
    cancelled: impl FnMut() -> bool,
) -> Result<DiffusionNetworkAdvanceCandidate, Error> {
    advance_with_kernel(
        original,
        links,
        duration_s,
        endpoints,
        physical,
        controls,
        None,
        |state, dt, _, cancelled| {
            candidate_before_entropy_gate(state, links, dt, endpoints, physical, driving, cancelled)
                .map_err(TrialFailure::Physical)
        },
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
pub(super) fn advance_binary_candidate(
    original: &BTreeMap<u64, ClosedGasParcel>,
    links: &[DiffusionNetworkLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    physical: DiffusionNetworkSolve,
    controls: BinaryDiffusionNetworkAdvanceControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<BinaryDiffusionNetworkAdvanceCandidate, Error> {
    poll(&mut cancelled)?;
    if !(1..=1_000_000).contains(&controls.maximum_pair_evaluations) {
        return Err(Error::Invalid(
            "bounded positive binary pair-work budget required",
        ));
    }
    advance_with_kernel(
        original,
        links,
        duration_s,
        endpoints,
        physical,
        controls.temporal,
        Some(controls.maximum_pair_evaluations),
        |state, dt, work, cancelled| {
            let remaining = controls.maximum_pair_evaluations - work.pair_evaluations;
            let required = if links.is_empty() {
                0
            } else {
                2 * links.len() - 1
            };
            if required > remaining {
                return Err(TrialFailure::PairWork);
            }
            binary::advance_counted_candidate(
                state,
                links,
                dt,
                endpoints,
                physical,
                BinaryDiffusionNetworkControls {
                    substeps: 1,
                    maximum_pair_evaluations: remaining.clamp(1, 100_000),
                },
                &mut work.pair_evaluations,
                cancelled,
            )
            .map_err(TrialFailure::Binary)
        },
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn advance_with_kernel<C: Candidate>(
    original: &BTreeMap<u64, ClosedGasParcel>,
    links: &[DiffusionNetworkLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    physical: DiffusionNetworkSolve,
    controls: DiffusionNetworkAdvanceControls,
    maximum_pair_evaluations: Option<usize>,
    mut evaluate: impl FnMut(
        &BTreeMap<u64, ClosedGasParcel>,
        f64,
        &mut Work,
        &mut dyn FnMut() -> bool,
    ) -> Result<C, TrialFailure>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<DiffusionNetworkAdvanceCandidate<C::Receipt>, Error> {
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
    let mut work = Work::default();
    let (mut accepted, mut rejected) = (0usize, 0usize);
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
                dt,
                physical,
                &mut work,
                controls.maximum_evaluations,
                &mut evaluate,
                &mut cancelled,
            )?;
            let first = trial(
                &state,
                half,
                physical,
                &mut work,
                controls.maximum_evaluations,
                &mut evaluate,
                &mut cancelled,
            )?;
            let second = trial(
                first.parcels(),
                half,
                physical,
                &mut work,
                controls.maximum_evaluations,
                &mut evaluate,
                &mut cancelled,
            )?;
            Ok((coarse, first, second))
        })();
        let refusal = match attempted {
            Err(TrialFailure::Work) => {
                return Err(Error::WorkLimit {
                    elapsed_s: elapsed,
                    accepted_intervals: accepted,
                    evaluations: work.evaluations,
                });
            }
            Err(TrialFailure::PairWork) => {
                return Err(Error::PairWorkLimit {
                    elapsed_s: elapsed,
                    accepted_intervals: accepted,
                    pair_evaluations: work.pair_evaluations,
                    maximum_pair_evaluations: maximum_pair_evaluations
                        .expect("binary trials own a pair-work budget"),
                });
            }
            Err(TrialFailure::Binary(error)) => {
                if error == BinaryDiffusionNetworkError::Cancelled {
                    return Err(Error::Cancelled);
                }
                let retry = retryable_binary(&error);
                let error = Error::BinaryStep {
                    elapsed_s: elapsed,
                    step_s: dt,
                    error,
                };
                if !retry {
                    return Err(error);
                }
                error
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
                    ledger.record(&first, &mut cancelled)?;
                    ledger.record(&second, &mut cancelled)?;
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
    let audit = audit_closed_window(
        original,
        &state,
        GasTransportSolve {
            exchange: physical.exchange,
            entropy_tolerance_j_k: physical.entropy_tolerance_j_k,
        },
        &mut cancelled,
    )
    .map_err(|error| match error {
        GasTransportError::Cancelled => Error::Cancelled,
        GasTransportError::Diffusion(DiffusionNetworkError::EntropyDecrease {
            change_j_k, ..
        }) => Error::Entropy {
            elapsed_s: elapsed,
            step_s: 0.0,
            change_j_k,
        },
        other => Error::Window(other),
    })?;
    if ledger.energy.value.abs() > physical.exchange.energy_tolerance_j {
        return Err(Error::Invalid(
            "complete diffusion window exceeds its species/mass/energy budget",
        ));
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
            evaluations: work.evaluations,
            pair_evaluations: work.pair_evaluations,
            minimum_accepted_substep_s: minimum,
            maximum_accepted_substep_s: maximum,
            maximum_species_error_ratio: max_species_error,
            maximum_temperature_error_ratio: max_temperature_error,
            transfers,
            total_mass_residual_kg: audit.total_mass_residual_kg,
            maximum_species_mass_residual_kg: audit.maximum_species_mass_residual_kg,
            total_energy_residual_j: ledger.energy.value,
            endpoint_energy_change_j: audit.energy_change_j,
            entropy_change_j_k: audit.entropy_change_j_k,
            last_substep: last.ok_or(Error::Invalid("diffusion window accepted no interval"))?,
        },
    })
}
