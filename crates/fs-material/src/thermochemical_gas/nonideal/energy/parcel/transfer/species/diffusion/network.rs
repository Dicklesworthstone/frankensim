//! One source snapshot and one conservative real-gas species/enthalpy update.
//!
//! All faces sample original inventories; incoming gas cannot fund withdrawal
//! or travel through a second face in this explicit interval. The original
//! partial PR enthalpies already contain flow work. Local friction entropy is
//! not added as heat and does not replace the actual endpoint entropy audit.

/// Complete real-gas windows on the shared conservative temporal controller.
pub mod adaptive;

use super::super::network::{
    self as exchange, RealGasSpeciesNetworkStep, SpeciesNetworkError, SpeciesNetworkSolve,
    SpeciesNetworkTransfer,
};
use super::super::{GasEnergyError, GasEnergySolve, IsochoricRealGasParcel, positive};
use super::{
    RealGasDiffusionFlux, RealGasDiffusionModel, RealGasDiffusionSolve, RealGasDiffusionWorkspace,
};
use crate::thermochemical_gas::energy::parcel::transfer::{
    snapshot::{MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS},
    species::{
        diffusion::GasDiffusionGeometry, entropy, network::transport::GasTransportWindowAudit,
    },
};
use fs_math::dd::Dd;
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug)]
pub struct RealGasDiffusionNetworkLink<'a> {
    pub id: u64,
    pub parcels: [u64; 2],
    pub model: &'a RealGasDiffusionModel<'a>,
    /// A zero area closes this path without querying inactive material data.
    pub geometry: GasDiffusionGeometry,
    pub solve: RealGasDiffusionSolve,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RealGasDiffusionNetworkSolve {
    pub exchange: SpeciesNetworkSolve,
    /// Absolute allowed rounding decrease in the entire finite network [J/K].
    pub entropy_tolerance_j_k: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RealGasDiffusionNetworkLinkStep {
    pub id: u64,
    pub parcels: [u64; 2],
    /// None for a geometrically closed face.
    pub flux: Option<RealGasDiffusionFlux>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RealGasDiffusionNetworkStep {
    pub duration_s: f64,
    pub links: Vec<RealGasDiffusionNetworkLinkStep>,
    pub exchange: RealGasSpeciesNetworkStep,
    /// Independent actual original-to-final U change with species references
    /// removed before arithmetic, including the original PR departures.
    pub energy_change_j: f64,
    /// Original-to-final source-relative ideal mixing plus PR density entropy.
    pub entropy_change_j_k: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RealGasDiffusionNetworkCandidate {
    parcels: BTreeMap<u64, IsochoricRealGasParcel>,
    step: RealGasDiffusionNetworkStep,
}
impl RealGasDiffusionNetworkCandidate {
    pub fn parcels(&self) -> &BTreeMap<u64, IsochoricRealGasParcel> {
        &self.parcels
    }
    pub const fn step(&self) -> &RealGasDiffusionNetworkStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        BTreeMap<u64, IsochoricRealGasParcel>,
        RealGasDiffusionNetworkStep,
    ) {
        (self.parcels, self.step)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum RealGasDiffusionNetworkError {
    Invalid(&'static str),
    Face {
        id: u64,
        error: GasEnergyError,
    },
    Parcel {
        id: u64,
        error: GasEnergyError,
    },
    Exchange(SpeciesNetworkError),
    Energy {
        change_j: f64,
        tolerance_j: f64,
    },
    Entropy {
        change_j_k: f64,
        tolerance_j_k: f64,
    },
    /// Counted adaptive probes exhausted their original face-call budget.
    FaceWorkLimit {
        face_evaluations: usize,
        maximum_face_evaluations: usize,
    },
    Cancelled,
}
impl core::fmt::Display for RealGasDiffusionNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "real-gas diffusion network: {self:?}")
    }
}
impl std::error::Error for RealGasDiffusionNetworkError {}
use RealGasDiffusionNetworkError as Error;
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn face_error(id: u64, error: GasEnergyError) -> Error {
    if error == GasEnergyError::Cancelled {
        Error::Cancelled
    } else {
        Error::Face { id, error }
    }
}

/// Advance all molecular faces together, then admit actual source endpoints,
/// each species, total mass, original-law energy and actual whole-network S.
/// Source refusal, overdraw, nonconvergence or cancellation returns no endpoint.
/// This supplies one first-order explicit interval, not an implicit stability,
/// momentum, phase-equilibrium or global temporal-error guarantee.
pub fn step_candidate(
    parcels: &BTreeMap<u64, IsochoricRealGasParcel>,
    links: &[RealGasDiffusionNetworkLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: RealGasDiffusionNetworkSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<RealGasDiffusionNetworkCandidate, Error> {
    step_counted_candidate(
        parcels,
        links,
        duration_s,
        endpoints,
        controls,
        &mut BTreeMap::new(),
        &mut 0,
        usize::MAX,
        cancelled,
    )
}

// One original face/update path serves finite and adaptive transport. Work is
// spent immediately before every actual constitutive call, including failures.
// Scratch survives coarse/fine/rejected probes, never their physical states.
#[allow(clippy::too_many_arguments)]
fn step_counted_candidate(
    parcels: &BTreeMap<u64, IsochoricRealGasParcel>,
    links: &[RealGasDiffusionNetworkLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: RealGasDiffusionNetworkSolve,
    work: &mut BTreeMap<usize, RealGasDiffusionWorkspace>,
    face_evaluations: &mut usize,
    maximum_face_evaluations: usize,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RealGasDiffusionNetworkCandidate, Error> {
    poll(&mut cancelled)?;
    if !positive(duration_s)
        || !positive(controls.entropy_tolerance_j_k)
        || !positive(controls.exchange.mass_tolerance_kg)
        || !positive(controls.exchange.energy_tolerance_j)
        || !positive(controls.exchange.maximum_species_fraction)
        || controls.exchange.maximum_species_fraction >= 1.0
        || parcels.is_empty()
        || parcels.len() > MAX_SNAPSHOT_PARCELS
        || links.len() > MAX_SNAPSHOT_TRANSFERS
        || !parcels.keys().eq(endpoints.keys())
    {
        return Err(Error::Invalid(
            "bounded real-gas diffusion graph and positive physical controls required",
        ));
    }
    for (&id, solve) in endpoints {
        poll(&mut cancelled)?;
        solve
            .validate()
            .map_err(|error| Error::Parcel { id, error })?;
    }
    let mut ordered = BTreeMap::new();
    for link in links {
        poll(&mut cancelled)?;
        if link.parcels[0] == link.parcels[1]
            || link.parcels.iter().any(|id| !parcels.contains_key(id))
            || !link.geometry.area_m2.is_finite()
            || link.geometry.area_m2 < 0.0
            || !positive(link.geometry.distance_m)
            || ordered.insert(link.id, link).is_some()
        {
            return Err(Error::Invalid(
                "unique real-gas diffusion IDs, distinct owners and physical geometry required",
            ));
        }
        link.solve
            .validate()
            .map_err(|error| face_error(link.id, error))?;
    }
    let mut transfers = Vec::new();
    let mut steps = Vec::with_capacity(links.len());
    for (&id, link) in &ordered {
        poll(&mut cancelled)?;
        let [left, right] = link.parcels;
        let flux = if link.geometry.area_m2 == 0.0 {
            None
        } else {
            if *face_evaluations == maximum_face_evaluations {
                return Err(Error::FaceWorkLimit {
                    face_evaluations: *face_evaluations,
                    maximum_face_evaluations,
                });
            }
            let n = link.model.species().len();
            if let std::collections::btree_map::Entry::Vacant(entry) = work.entry(n) {
                entry.insert(
                    link.model
                        .workspace()
                        .map_err(|error| face_error(id, error))?,
                );
            }
            *face_evaluations += 1;
            let flux = link
                .model
                .face_flux(
                    &parcels[&left],
                    &parcels[&right],
                    link.geometry,
                    link.solve,
                    work.get_mut(&n).unwrap(),
                    &mut cancelled,
                )
                .map_err(|error| face_error(id, error))?;
            for species in &flux.species {
                poll(&mut cancelled)?;
                let rate = species.mass_rate_left_to_right_kg_s;
                let mass = rate * duration_s;
                if !mass.is_finite() || (rate != 0.0 && mass == 0.0) {
                    return Err(Error::Invalid(
                        "integrated real-gas diffusion mass is unrepresentable",
                    ));
                }
                if mass == 0.0 {
                    continue;
                }
                if transfers.len() == MAX_SNAPSHOT_TRANSFERS {
                    return Err(Error::Invalid(
                        "real-gas diffusion exceeds the shared species-transfer bound",
                    ));
                }
                let (donor, receiver) = if mass > 0.0 {
                    (left, right)
                } else {
                    (right, left)
                };
                transfers.push(SpeciesNetworkTransfer {
                    donor,
                    receiver,
                    species: species.species.clone(),
                    mass_kg: mass.abs(),
                });
            }
            Some(flux)
        };
        steps.push(RealGasDiffusionNetworkLinkStep {
            id,
            parcels: link.parcels,
            flux,
        });
    }
    let candidate = exchange::exchange_candidate(
        parcels,
        &transfers,
        endpoints,
        controls.exchange,
        &mut cancelled,
    )
    .map_err(|error| {
        if error == SpeciesNetworkError::Cancelled {
            Error::Cancelled
        } else {
            Error::Exchange(error)
        }
    })?;
    for link in ordered.values() {
        poll(&mut cancelled)?;
        if link.geometry.area_m2 > 0.0 {
            link.model
                .admit_pair(
                    &candidate.parcels()[&link.parcels[0]],
                    &candidate.parcels()[&link.parcels[1]],
                    link.solve,
                    &mut cancelled,
                )
                .map_err(|error| face_error(link.id, error))?;
        }
    }
    let audit = audit(parcels, candidate.parcels(), controls, &mut cancelled)?;
    poll(&mut cancelled)?;
    let (parcels, exchange) = candidate.into_parts();
    Ok(RealGasDiffusionNetworkCandidate {
        parcels,
        step: RealGasDiffusionNetworkStep {
            duration_s,
            links: steps,
            exchange,
            energy_change_j: audit.energy_change_j,
            entropy_change_j_k: audit.entropy_change_j_k,
        },
    })
}

fn view(parcel: &IsochoricRealGasParcel) -> entropy::ClosedGasInventory<'_> {
    entropy::ClosedGasInventory {
        model: parcel.model().source(),
        mass_kg: parcel.mass_kg(),
        molar_mass_kg_mol: parcel
            .state()
            .source()
            .thermodynamics()
            .properties()
            .molar_mass()
            .value(),
        volume_m3: parcel.volume_m3(),
        temperature_k: parcel.state().temperature_k(),
    }
}

fn audit(
    original: &BTreeMap<u64, IsochoricRealGasParcel>,
    next: &BTreeMap<u64, IsochoricRealGasParcel>,
    controls: RealGasDiffusionNetworkSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<GasTransportWindowAudit, Error> {
    if !original.keys().eq(next.keys()) {
        return Err(Error::Invalid(
            "real-gas diffusion changed physical owner identities",
        ));
    }
    let ideal = entropy::closed_inventory_views::<true>(
        original
            .iter()
            .map(|(&id, old)| (id, view(old), view(&next[&id]))),
        controls.exchange.mass_tolerance_kg,
        &mut *cancelled,
    )
    .map_err(|error| match error {
        entropy::Error::Invalid(message) => Error::Invalid(message),
        entropy::Error::Parcel { id, error } => Error::Parcel { id, error },
        entropy::Error::Cancelled => Error::Cancelled,
    })?;
    if ideal.total_mass_residual_kg.abs() > controls.exchange.mass_tolerance_kg
        || ideal.maximum_species_mass_residual_kg > controls.exchange.mass_tolerance_kg
    {
        return Err(Error::Invalid(
            "real-gas diffusion final inventory did not conserve mass/species",
        ));
    }
    let mut energy = Dd::from_f64(ideal.energy_change_j);
    let mut entropy = Dd::from_f64(ideal.entropy_change_j_k);
    for (&id, old) in original {
        poll(cancelled)?;
        let new = &next[&id];
        if old == new {
            continue;
        }
        let n0 = old.mass_kg() / view(old).molar_mass_kg_mol;
        let n1 = new.mass_kg() / view(new).molar_mass_kg_mol;
        energy = energy
            + Dd::from_f64(n1) * Dd::from_f64(new.state().root().residual_internal_energy_j_mol())
            - Dd::from_f64(n0) * Dd::from_f64(old.state().root().residual_internal_energy_j_mol());
        // The shared ideal audit uses nRT/V. Its density-chart PR departure
        // excludes R ln Z, which would otherwise count the pressure correction twice.
        entropy = entropy
            + Dd::from_f64(n1)
                * Dd::from_f64(new.state().root().residual_entropy_at_density_j_mol_k())
            - Dd::from_f64(n0)
                * Dd::from_f64(old.state().root().residual_entropy_at_density_j_mol_k());
    }
    let energy_change = energy.to_f64();
    let entropy_change = entropy.to_f64();
    if ![
        energy.hi,
        energy.lo,
        entropy.hi,
        entropy.lo,
        energy_change,
        entropy_change,
    ]
    .iter()
    .all(|value| value.is_finite())
    {
        return Err(Error::Invalid(
            "real-gas endpoint energy/entropy audit is unrepresentable",
        ));
    }
    if energy_change.abs() > controls.exchange.energy_tolerance_j {
        return Err(Error::Energy {
            change_j: energy_change,
            tolerance_j: controls.exchange.energy_tolerance_j,
        });
    }
    if entropy_change < -controls.entropy_tolerance_j_k {
        return Err(Error::Entropy {
            change_j_k: entropy_change,
            tolerance_j_k: controls.entropy_tolerance_j_k,
        });
    }
    Ok(GasTransportWindowAudit {
        total_mass_residual_kg: ideal.total_mass_residual_kg,
        maximum_species_mass_residual_kg: ideal.maximum_species_mass_residual_kg,
        energy_change_j: energy_change,
        entropy_change_j_k: entropy_change,
    })
}
