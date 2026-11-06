//! Conservative exponential pair splitting on finite binary gas networks.
//!
//! Every edge calls the original binary molar-frame relaxation and reconstructs
//! both material states before another edge uses them. A forward half sweep,
//! one full central edge, and the reverse half sweep give symmetric composition.
//! The work has no explicit diffusion CFL restriction. Splitting error still
//! requires refinement; state-dependent frozen pair laws remain first order.

use super::super::super::network::transport::{
    GasTransportError, GasTransportSolve, audit_closed_window,
};
use super::super::super::{GasSpeciesExchangeStep, SpeciesId, Sum};
use super::super::binary::BinaryDiffusionStep;
use super::adaptive::DiffusionNetworkIntegratedTransfer;
use super::{
    ClosedGasParcel, DiffusionNetworkLink, DiffusionNetworkSolve, GasDiffusionFlux,
    GasDiffusionSolve, GasEnergyError, GasEnergySolve, GasSpeciesExchangeSolve,
    MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS, positive,
};
use std::collections::{BTreeMap, BTreeSet};

/// Explicit time resolution and a bound on actual source-aware pair updates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BinaryDiffusionNetworkControls {
    /// Number of symmetric sweeps across the requested physical duration.
    /// Must be in 1..=100,000. No local or global time-error estimate is implied.
    pub substeps: usize,
    /// Maximum pair calls, including closed faces and both half sweeps.
    /// Must be in 1..=100,000. For E>0 links, each sweep costs 2*E-1 calls;
    /// an empty graph costs zero. Closed faces count as pair calls.
    pub maximum_pair_evaluations: usize,
}

/// One actual pair update from the last symmetric sweep, in execution order.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryDiffusionNetworkFaceStep {
    pub id: u64,
    /// Canonical increasing endpoint IDs, independent of authored orientation.
    pub parcels: [u64; 2],
    /// Source-qualified instantaneous flux BEFORE this particular pair update.
    /// Its product with duration is not the integrated exponential transfer.
    pub initial_flux: GasDiffusionFlux,
    pub relaxation: BinaryDiffusionStep,
    pub exchange: GasSpeciesExchangeStep,
    pub entropy_change_j_k: f64,
}

/// Complete-window balances; detailed pair receipts retain only the last sweep.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryDiffusionNetworkReport {
    pub duration_s: f64,
    pub substeps: usize,
    pub pair_evaluations: usize,
    /// Actual gross species/partial-enthalpy transfers, not initial flux times dt.
    pub transfers: Vec<DiffusionNetworkIntegratedTransfer>,
    pub total_mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    /// Compensated sum of original-caloric pair residuals after reference cancellation.
    pub total_energy_residual_j: f64,
    /// Independent original-to-final source-relative internal-energy change.
    pub endpoint_energy_change_j: f64,
    pub entropy_change_j_k: f64,
    pub last_sweep: Vec<BinaryDiffusionNetworkFaceStep>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BinaryDiffusionNetworkCandidate {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    report: BinaryDiffusionNetworkReport,
}
impl BinaryDiffusionNetworkCandidate {
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub const fn report(&self) -> &BinaryDiffusionNetworkReport {
        &self.report
    }
    pub fn into_parts(self) -> (BTreeMap<u64, ClosedGasParcel>, BinaryDiffusionNetworkReport) {
        (self.parcels, self.report)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum BinaryDiffusionNetworkError {
    Invalid(&'static str),
    WorkLimit {
        required_pair_evaluations: usize,
        maximum_pair_evaluations: usize,
    },
    Parcel {
        id: u64,
        error: GasEnergyError,
    },
    Face {
        id: u64,
        /// Zero-based symmetric sweep containing the failed pair update.
        substep: usize,
        error: GasEnergyError,
    },
    Window(GasTransportError),
    Cancelled,
}
impl core::fmt::Display for BinaryDiffusionNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "exponential binary diffusion network: {self:?}")
    }
}
impl std::error::Error for BinaryDiffusionNetworkError {}
use BinaryDiffusionNetworkError as Error;

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn arithmetic(_: GasEnergyError) -> Error {
    Error::Invalid("unrepresentable integrated binary diffusion account")
}

#[derive(Default)]
struct Ledger {
    energy: Sum,
    transfers: BTreeMap<(u64, u64, SpeciesId), (Sum, Sum)>,
}
impl Ledger {
    fn record(
        &mut self,
        parcels: [u64; 2],
        step: &GasSpeciesExchangeStep,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        self.energy
            .add(step.energy_residual_j)
            .map_err(arithmetic)?;
        for row in &step.species {
            poll(cancelled)?;
            let mass = row.transfer.mass_left_to_right_kg;
            if mass == 0.0 {
                continue;
            }
            let (donor, receiver, sign) = if mass > 0.0 {
                (parcels[0], parcels[1], 1.0)
            } else {
                (parcels[1], parcels[0], -1.0)
            };
            let key = (donor, receiver, row.transfer.species.clone());
            if !self.transfers.contains_key(&key) && self.transfers.len() == MAX_SNAPSHOT_TRANSFERS
            {
                return Err(Error::Invalid(
                    "binary diffusion transfer report exceeds its bound",
                ));
            }
            let totals = self.transfers.entry(key).or_default();
            totals.0.add(mass.abs()).map_err(arithmetic)?;
            totals
                .1
                .add(sign * row.enthalpy_left_to_right_j)
                .map_err(arithmetic)?;
        }
        Ok(())
    }
}

/// Advance a closed binary gas network through the complete positive duration.
///
/// Each open edge must join exactly two original species and supply their one
/// binary material. Branches, cycles, disconnected components and unequal finite
/// volumes are supported. There must be at most one edge per unordered endpoint
/// pair; distinct parallel source laws are not combined. IDs are unique receipt
/// labels only. Canonical endpoint ordering fixes the arithmetic under input
/// permutations, ID renaming and edge-orientation reversal.
///
/// Each sweep composes half updates in increasing physical-edge order, a full
/// update of the last edge, then the reverse half updates. Every invocation
/// resamples actual face T/p/D, carries actual upstream species enthalpy, and
/// retains the original source, withdrawal, endpoint and pair-entropy gates.
/// Arriving species can leave in a later split operation. No concentration is
/// clipped or seeded and no extra thermal reservoir or pressure work is added.
///
/// A single isothermal/isobaric constant-D edge has the original analytic
/// exponential endpoint. Composition of exact constant-state pair flows is
/// second order under time refinement. General source-dependent/nonisothermal
/// pair updates remain frozen-coefficient FIRST ORDER; symmetry alone does not
/// raise their order. Large stiff sweeps can have large network splitting error.
/// Choose substeps by refinement; this API supplies no temporal-error estimate.
///
/// Half sweeps count toward the explicit work bound. Only the last sweep retains
/// detailed source receipts. Actual original-to-final species, mass, energy and
/// entropy changes, and the compensated pair-energy ledger, must fit the SAME
/// fixed absolute budgets over the whole window. Late failure/cancellation
/// returns no partial network. Bulk flow, pressure diffusion, reactions and
/// moving volumes require their own physical owners.
#[allow(clippy::too_many_arguments)]
pub fn advance_candidate(
    original: &BTreeMap<u64, ClosedGasParcel>,
    links: &[DiffusionNetworkLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    physical: DiffusionNetworkSolve,
    controls: BinaryDiffusionNetworkControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<BinaryDiffusionNetworkCandidate, Error> {
    poll(&mut cancelled)?;
    if !positive(duration_s)
        || !(1..=100_000).contains(&controls.substeps)
        || !(1..=100_000).contains(&controls.maximum_pair_evaluations)
        || original.is_empty()
        || original.len() > MAX_SNAPSHOT_PARCELS
        || links.len() > MAX_SNAPSHOT_TRANSFERS
        || !original.keys().eq(endpoints.keys())
        || !positive(physical.entropy_tolerance_j_k)
        || !positive(physical.exchange.mass_tolerance_kg)
        || !positive(physical.exchange.energy_tolerance_j)
        || !positive(physical.exchange.maximum_species_fraction)
        || physical.exchange.maximum_species_fraction >= 1.0
    {
        return Err(Error::Invalid(
            "bounded graph, exact endpoint map and positive physical/time controls required",
        ));
    }
    for (&id, solve) in endpoints {
        poll(&mut cancelled)?;
        solve
            .validate()
            .map_err(|error| Error::Parcel { id, error })?;
        if !(solve.minimum_temperature_k..=solve.maximum_temperature_k)
            .contains(&original[&id].state().state().temperature)
        {
            return Err(Error::Invalid(
                "initial gas temperature outside endpoint bracket",
            ));
        }
    }
    let mut ids = BTreeSet::new();
    let mut ordered = BTreeMap::new();
    for link in links {
        poll(&mut cancelled)?;
        let [a, b] = link.parcels;
        if a == b
            || !original.contains_key(&a)
            || !original.contains_key(&b)
            || !ids.insert(link.id)
            || ordered.insert((a.min(b), a.max(b)), link).is_some()
        {
            return Err(Error::Invalid(
                "unique face IDs and unique unordered pairs of existing endpoints required",
            ));
        }
    }
    let ordered: Vec<_> = ordered.into_iter().collect();
    let per_sweep = if ordered.is_empty() {
        0
    } else {
        2 * ordered.len() - 1
    };
    let required = per_sweep
        .checked_mul(controls.substeps)
        .ok_or(Error::Invalid("binary diffusion work count overflow"))?;
    if required > controls.maximum_pair_evaluations {
        return Err(Error::WorkLimit {
            required_pair_evaluations: required,
            maximum_pair_evaluations: controls.maximum_pair_evaluations,
        });
    }
    let nominal = duration_s / controls.substeps as f64;
    if !positive(nominal) {
        return Err(Error::Invalid(
            "binary diffusion substep is unrepresentable",
        ));
    }
    let mut state = original.clone();
    let mut ledger = Ledger::default();
    let mut last_sweep = Vec::with_capacity(per_sweep);
    let mut elapsed = 0.0;
    let mut evaluations = 0;
    if !ordered.is_empty() {
        for substep in 0..controls.substeps {
            poll(&mut cancelled)?;
            let dt = if substep + 1 == controls.substeps {
                duration_s - elapsed
            } else {
                nominal
            };
            let half = 0.5 * dt;
            if !positive(dt)
                || elapsed + dt <= elapsed
                || (ordered.len() > 1 && (!positive(half) || half + half != dt))
            {
                return Err(Error::Invalid(
                    "binary diffusion cannot advance physical time",
                ));
            }
            let last = ordered.len() - 1;
            for slot in 0..per_sweep {
                poll(&mut cancelled)?;
                let edge = if slot <= last {
                    slot
                } else {
                    per_sweep - 1 - slot
                };
                let ((left, right), link) = ordered[edge];
                let pair_duration = if edge == last { dt } else { half };
                let solve = GasDiffusionSolve {
                    maximum_temperature_difference_k: link.solve.maximum_temperature_difference_k,
                    maximum_pressure_difference_pa: link.solve.maximum_pressure_difference_pa,
                    gradient_tolerance_per_m: link.solve.gradient_tolerance_per_m,
                    molar_flux_tolerance_mol_m2_s: link.solve.molar_flux_tolerance_mol_m2_s,
                    entropy_tolerance_j_k: physical.entropy_tolerance_j_k,
                    exchange: GasSpeciesExchangeSolve {
                        gas: [endpoints[&left], endpoints[&right]],
                        maximum_species_fraction: physical.exchange.maximum_species_fraction,
                        mass_tolerance_kg: physical.exchange.mass_tolerance_kg,
                        energy_tolerance_j: physical.exchange.energy_tolerance_j,
                    },
                };
                evaluations += 1;
                let candidate = state[&left]
                    .diffuse_binary_exponential_with(
                        &state[&right],
                        link.geometry,
                        pair_duration,
                        link.pairs,
                        solve,
                        &mut cancelled,
                    )
                    .map_err(|error| {
                        if error == GasEnergyError::Cancelled {
                            Error::Cancelled
                        } else {
                            Error::Face {
                                id: link.id,
                                substep,
                                error,
                            }
                        }
                    })?;
                let (diffusion, relaxation) = candidate.into_parts();
                let (exchange, initial_flux, entropy_change_j_k) = diffusion.into_parts();
                let (next_left, next_right, exchange) = exchange.into_parts();
                ledger.record([left, right], &exchange, &mut cancelled)?;
                state.insert(left, next_left);
                state.insert(right, next_right);
                if substep + 1 == controls.substeps {
                    last_sweep.push(BinaryDiffusionNetworkFaceStep {
                        id: link.id,
                        parcels: [left, right],
                        initial_flux,
                        relaxation,
                        exchange,
                        entropy_change_j_k,
                    });
                }
            }
            elapsed = if substep + 1 == controls.substeps {
                duration_s
            } else {
                elapsed + dt
            };
        }
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
    .map_err(|error| {
        if error == GasTransportError::Cancelled {
            Error::Cancelled
        } else {
            Error::Window(error)
        }
    })?;
    if ledger.energy.value.abs() > physical.exchange.energy_tolerance_j {
        return Err(Error::Invalid(
            "complete binary diffusion window exceeds its energy budget",
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
    poll(&mut cancelled)?;
    Ok(BinaryDiffusionNetworkCandidate {
        parcels: state,
        report: BinaryDiffusionNetworkReport {
            duration_s,
            substeps: controls.substeps,
            pair_evaluations: evaluations,
            transfers,
            total_mass_residual_kg: audit.total_mass_residual_kg,
            maximum_species_mass_residual_kg: audit.maximum_species_mass_residual_kg,
            total_energy_residual_j: ledger.energy.value,
            endpoint_energy_change_j: audit.energy_change_j,
            entropy_change_j_k: audit.entropy_change_j_k,
            last_sweep,
        },
    })
}
