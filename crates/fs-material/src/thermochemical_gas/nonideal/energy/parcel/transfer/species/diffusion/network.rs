//! One source snapshot and one conservative real-gas species/enthalpy update.
//!
//! All faces sample original inventories; incoming gas cannot fund withdrawal
//! or travel through a second face in this explicit interval. The original
//! partial PR enthalpies already contain flow work. Local friction entropy is
//! not added as heat and does not replace the actual endpoint entropy audit.

/// Complete real-gas windows on the shared conservative temporal controller.
pub mod adaptive;

use super::super::network::{
    self as exchange, RealGasSpeciesNetworkStep, SpeciesNetworkError, SpeciesNetworkHeatTransfer,
    SpeciesNetworkSolve, SpeciesNetworkTransfer,
};
use super::super::{GasEnergyError, GasEnergySolve, IsochoricRealGasParcel, positive};
use super::{
    RealGasDiffusionFlux, RealGasDiffusionModel, RealGasDiffusionSolve, RealGasDiffusionWorkspace,
    SpeciesDiffusionFlux, SpeciesId,
};
use crate::thermochemical_gas::energy::parcel::transfer::{
    snapshot::{MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS},
    species::{
        diffusion::GasDiffusionGeometry, entropy, network::transport::GasTransportWindowAudit,
    },
};
use fs_math::dd::Dd;
use std::collections::BTreeMap;

#[derive(Debug)]
pub struct RealGasDiffusionNetworkLink<'a, Model = RealGasDiffusionModel<'a>> {
    pub id: u64,
    pub parcels: [u64; 2],
    pub model: &'a Model,
    /// A zero area closes this path without querying inactive material data.
    pub geometry: GasDiffusionGeometry,
    pub solve: RealGasDiffusionSolve,
}
impl<Model> Copy for RealGasDiffusionNetworkLink<'_, Model> {}
impl<Model> Clone for RealGasDiffusionNetworkLink<'_, Model> {
    fn clone(&self) -> Self {
        *self
    }
}

// Both constitutive laws lower to the SAME simultaneous species/heat owner.
// The interface carries no state or source substitutions into its transactions.
pub(super) trait MolecularFace {
    type Flux;
    fn species(&self) -> &[SpeciesId];
    fn workspace(&self) -> Result<RealGasDiffusionWorkspace, GasEnergyError>;
    fn sample(
        &self,
        left: &IsochoricRealGasParcel,
        right: &IsochoricRealGasParcel,
        geometry: GasDiffusionGeometry,
        controls: RealGasDiffusionSolve,
        work: &mut RealGasDiffusionWorkspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self::Flux, GasEnergyError>;
    fn admit(
        &self,
        left: &IsochoricRealGasParcel,
        right: &IsochoricRealGasParcel,
        controls: RealGasDiffusionSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), GasEnergyError>;
    fn species_flux(flux: &Self::Flux) -> &[SpeciesDiffusionFlux];
    fn heat_rate(_flux: &Self::Flux) -> f64 {
        0.0
    }
    /// The mass-frame residual the constitutive owner published with this
    /// flux [kg/(m2 s)], its compensated sum of the published species mass
    /// fluxes. No default: a face that publishes none does not lower.
    fn frame_residual(flux: &Self::Flux) -> f64;
}
impl MolecularFace for RealGasDiffusionModel<'_> {
    type Flux = RealGasDiffusionFlux;
    fn species(&self) -> &[SpeciesId] {
        RealGasDiffusionModel::species(self)
    }
    fn workspace(&self) -> Result<RealGasDiffusionWorkspace, GasEnergyError> {
        RealGasDiffusionModel::workspace(self)
    }
    fn sample(
        &self,
        left: &IsochoricRealGasParcel,
        right: &IsochoricRealGasParcel,
        geometry: GasDiffusionGeometry,
        controls: RealGasDiffusionSolve,
        work: &mut RealGasDiffusionWorkspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self::Flux, GasEnergyError> {
        self.face_flux(left, right, geometry, controls, work, cancelled)
    }
    fn admit(
        &self,
        left: &IsochoricRealGasParcel,
        right: &IsochoricRealGasParcel,
        controls: RealGasDiffusionSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), GasEnergyError> {
        self.admit_pair(left, right, controls, cancelled)
    }
    fn species_flux(flux: &Self::Flux) -> &[SpeciesDiffusionFlux] {
        &flux.species
    }
    fn frame_residual(flux: &Self::Flux) -> f64 {
        flux.mass_frame_residual_kg_m2_s
    }
}

/// Unit roundoff of binary64.
const U: f64 = f64::EPSILON / 2.0;
/// `gamma_3 = 3u / (1 - 3u)` (Higham, Accuracy and Stability of Numerical
/// Algorithms, 2nd ed., Lemma 3.1).
const GAMMA_3: f64 = 3.0 * U / (1.0 - 3.0 * U);

// Both molecular laws are mass-frame laws: the owner's published species mass
// fluxes J_i sum to zero within its published frame residual r, which is its
// compensated sum of those J_i. A face therefore moves no net mass by its law.
// This lowering rounds m_i = fl(fl(J_i A) dt), so (Higham 2nd ed., Lemma 3.1,
// and section 4.3, (4.8), for the compensated sum)
//   |sum_i m_i| <= |r| A dt + (4u + O(n u^2)) sum_i |m_i|,
// 2u from the two roundings of each m_i and 2u from the owner's compensated
// sum. The audit sums the m_i in double-double and rounds once (u); 1 + gamma_3
// covers that rounding and the two of |r| A dt, and the fifth u covers the
// O(n u^2) terms and this bound's own evaluation. A face outside the bound is
// refused. Inside it, each parcel of the face is allowed the bound plus
// 6u sum_i |m_i| for the network's three compensated sums of these rows (the
// per-key merge, the per-species budget and the parcel's net), and the network
// keeps a parcel's mass exactly when its computed net lies within the sum of
// its faces' allowances. On the general Maxwell-Stefan path r is the frame row
// residual of its LU solve, which the owner has already admitted against its
// own tolerance; it is never recomputed here.
pub(super) fn face_allowance(
    masses: &[f64],
    frame_residual: f64,
    area_m2: f64,
    duration_s: f64,
) -> Result<f64, GasEnergyError> {
    let (mut net, mut gross) = (Dd::ZERO, Dd::ZERO);
    for &mass in masses {
        net = net + Dd::from_f64(mass);
        gross = gross + Dd::from_f64(mass.abs());
    }
    let (net, gross) = (net.to_f64(), gross.to_f64());
    let rate = (1.0 + GAMMA_3) * (frame_residual.abs() * area_m2 * duration_s);
    let bound = rate + 5.0 * U * gross;
    if ![net, gross, rate, bound].iter().all(|v| v.is_finite()) || net.abs() > bound {
        return Err(GasEnergyError::Invalid(
            "real-gas mass-frame face nets more mass than its rounding bound",
        ));
    }
    Ok(bound + 6.0 * U * gross)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RealGasDiffusionNetworkSolve {
    pub exchange: SpeciesNetworkSolve,
    /// Absolute allowed rounding decrease in the entire finite network [J/K].
    pub entropy_tolerance_j_k: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RealGasDiffusionNetworkLinkStep<Flux = RealGasDiffusionFlux> {
    pub id: u64,
    pub parcels: [u64; 2],
    /// None for a geometrically closed face.
    pub flux: Option<Flux>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RealGasDiffusionNetworkStep<Flux = RealGasDiffusionFlux> {
    pub duration_s: f64,
    pub links: Vec<RealGasDiffusionNetworkLinkStep<Flux>>,
    pub exchange: RealGasSpeciesNetworkStep,
    /// Independent actual original-to-final U change with species references
    /// removed before arithmetic, including the original PR departures.
    pub energy_change_j: f64,
    /// Original-to-final source-relative ideal mixing plus PR density entropy.
    pub entropy_change_j_k: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RealGasDiffusionNetworkCandidate<Flux = RealGasDiffusionFlux> {
    parcels: BTreeMap<u64, IsochoricRealGasParcel>,
    step: RealGasDiffusionNetworkStep<Flux>,
}
impl<Flux> RealGasDiffusionNetworkCandidate<Flux> {
    pub fn parcels(&self) -> &BTreeMap<u64, IsochoricRealGasParcel> {
        &self.parcels
    }
    pub const fn step(&self) -> &RealGasDiffusionNetworkStep<Flux> {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        BTreeMap<u64, IsochoricRealGasParcel>,
        RealGasDiffusionNetworkStep<Flux>,
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
pub(super) fn step_counted_candidate<Model: MolecularFace>(
    parcels: &BTreeMap<u64, IsochoricRealGasParcel>,
    links: &[RealGasDiffusionNetworkLink<'_, Model>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: RealGasDiffusionNetworkSolve,
    work: &mut BTreeMap<usize, RealGasDiffusionWorkspace>,
    face_evaluations: &mut usize,
    maximum_face_evaluations: usize,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RealGasDiffusionNetworkCandidate<Model::Flux>, Error> {
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
    let mut heat = Vec::new();
    let mut steps = Vec::with_capacity(links.len());
    let mut face_masses = Vec::new();
    let mut allowances: BTreeMap<u64, Dd> = BTreeMap::new();
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
                .sample(
                    &parcels[&left],
                    &parcels[&right],
                    link.geometry,
                    link.solve,
                    work.get_mut(&n).unwrap(),
                    &mut cancelled,
                )
                .map_err(|error| face_error(id, error))?;
            face_masses.clear();
            for species in Model::species_flux(&flux) {
                poll(&mut cancelled)?;
                let rate = species.mass_rate_left_to_right_kg_s;
                let mass = rate * duration_s;
                if !mass.is_finite() || (rate != 0.0 && mass == 0.0) {
                    return Err(Error::Invalid(
                        "integrated real-gas diffusion mass is unrepresentable",
                    ));
                }
                face_masses.push(mass);
                if mass == 0.0 {
                    continue;
                }
                if transfers.len() + heat.len() == MAX_SNAPSHOT_TRANSFERS {
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
            let allowance = face_allowance(
                &face_masses,
                Model::frame_residual(&flux),
                link.geometry.area_m2,
                duration_s,
            )
            .map_err(|error| face_error(id, error))?;
            for parcel in link.parcels {
                let total = allowances.entry(parcel).or_insert(Dd::ZERO);
                *total = *total + Dd::from_f64(allowance);
            }
            let rate = Model::heat_rate(&flux);
            let energy = rate * duration_s;
            if !energy.is_finite() || (rate != 0.0 && energy == 0.0) {
                return Err(Error::Invalid(
                    "integrated real-gas molecular heat is unrepresentable",
                ));
            }
            if energy != 0.0 {
                if transfers.len() + heat.len() == MAX_SNAPSHOT_TRANSFERS {
                    return Err(Error::Invalid(
                        "real-gas species/heat exceed their shared work bound",
                    ));
                }
                heat.push(SpeciesNetworkHeatTransfer {
                    parcels: link.parcels,
                    heat_left_to_right_j: energy,
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
    // Summed in double-double and rounded once: within the fifth u of each
    // face's allowance for every face of a bounded snapshot.
    let allowances = allowances
        .into_iter()
        .map(|(parcel, total)| (parcel, total.to_f64()))
        .collect();
    let candidate = exchange::exchange_candidate_with_heat(
        parcels,
        &transfers,
        &heat,
        &allowances,
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
                .admit(
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
