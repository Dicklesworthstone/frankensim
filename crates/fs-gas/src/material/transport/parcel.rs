//! Source-backed molecular faces on finite, rigid material gas inventories.
//!
//! The original continuum face owns mass-frame Maxwell--Stefan, Soret, Fourier
//! and Dufour fluxes at zero barycentric velocity. The material snapshot owner
//! owns gross withdrawals, source catalog joins, caloric reconstruction and the
//! final closed energy/entropy gates. This is a first-order candidate operator,
//! not a second timeline or a thermal-transpiration model for dusty-gas pores.

use crate::multicomponent::dissipation::{
    DissipationError, MaxwellStefanSolveV1, PlanarTransport, TransportFaceReport, TransportLaw,
};
use crate::multicomponent::{Error as FlowError, NormalPrimitive};
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::snapshot::{
    MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::{
    SpeciesNetworkHeatTransfer, SpeciesNetworkTransfer,
    transport::{self, GasTransportError, GasTransportPaths, GasTransportSolve, GasTransportStep},
};
use fs_qty::Temperature;
use fs_thermochem::SpeciesId;
use std::collections::BTreeMap;

/// Independently exposed continuum diffusion face; orientation is left to right.
/// This area must not also be supplied as a dusty-gas or duplicate diffusion path.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParcelMolecularLink {
    pub id: u64,
    pub parcels: [u64; 2],
    pub area_m2: f64,
    pub distance_m: f64,
}

/// Original numerical owners and fixed absolute physical acceptance budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParcelMolecularSolve {
    pub face: MaxwellStefanSolveV1,
    pub physical: GasTransportSolve,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ParcelMolecularError {
    Invalid(&'static str),
    SourceMismatch(SpeciesId),
    Composition(fs_qty::semantic::SemanticError),
    Face(DissipationError),
    Transport(GasTransportError),
    Cancelled,
}
impl core::fmt::Display for ParcelMolecularError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite molecular gas transport: {self:?}")
    }
}
impl std::error::Error for ParcelMolecularError {}
impl From<DissipationError> for ParcelMolecularError {
    fn from(error: DissipationError) -> Self {
        match error {
            DissipationError::Cancelled => Self::Cancelled,
            other => Self::Face(other),
        }
    }
}
impl From<GasTransportError> for ParcelMolecularError {
    fn from(error: GasTransportError) -> Self {
        match error {
            GasTransportError::Cancelled => Self::Cancelled,
            other => Self::Transport(other),
        }
    }
}
type Error = ParcelMolecularError;

/// Signed original-snapshot transfers; source identities remain in the prepared law.
#[derive(Clone, Debug, PartialEq)]
pub struct ParcelMolecularLinkStep {
    pub id: u64,
    pub parcels: [u64; 2],
    /// Integrated species mass in `ParcelMolecularStep::species` order [kg].
    pub mass_left_to_right_kg: Vec<f64>,
    /// Additional energy beyond original donor species enthalpy [J].
    pub heat_correction_left_to_right_j: f64,
    pub face: TransportFaceReport,
}
#[derive(Clone, Debug, PartialEq)]
pub struct ParcelMolecularStep {
    pub species: Vec<SpeciesId>,
    pub links: Vec<ParcelMolecularLinkStep>,
    /// Actual combined native endpoint balances, including all other paths.
    pub transport: GasTransportStep,
}
#[derive(Clone, Debug, PartialEq)]
pub struct ParcelMolecularCandidate {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    step: ParcelMolecularStep,
}
impl ParcelMolecularCandidate {
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub const fn step(&self) -> &ParcelMolecularStep {
        &self.step
    }
    pub fn into_parts(self) -> (BTreeMap<u64, ClosedGasParcel>, ParcelMolecularStep) {
        (self.parcels, self.step)
    }
}

/// Prepared original molecular owner usable with `MaterialSoretDufourTransport`.
///
/// Physical state is immutable until the caller accepts a complete candidate.
/// The law retains its normal last-successful-query receipts even during a
/// refused proposal; these may describe a state qualification rather than the
/// last transported face. Accepted gases are never replaced by that cache.
pub struct ParcelMolecularTransport<L: TransportLaw> {
    face: PlanarTransport<L>,
}
impl<L: TransportLaw> ParcelMolecularTransport<L> {
    pub const fn new(face: PlanarTransport<L>) -> Self {
        Self { face }
    }
    /// Prepared face, exact source bank and last successful constitutive query.
    pub const fn face(&self) -> &PlanarTransport<L> {
        &self.face
    }

    /// Sample every continuum and native path from one unchanged inventory.
    ///
    /// Continuum faces use zero bulk velocity and no viscous work. Native paths
    /// may describe separately exposed pressure/porous/membrane connections.
    /// Species always carry their original donor enthalpy once; the analytic
    /// face-to-donor correction makes the net energy equal the original shared
    /// molecular face flux without subtracting large formation references.
    /// The native closed entropy gate can refuse an overlarge explicit step;
    /// callers must refine time, never clip species or heat to force acceptance.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate(
        &mut self,
        parcels: &BTreeMap<u64, ClosedGasParcel>,
        links: &[ParcelMolecularLink],
        paths: GasTransportPaths<'_>,
        duration_s: f64,
        endpoints: &BTreeMap<u64, GasEnergySolve>,
        controls: ParcelMolecularSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ParcelMolecularCandidate, Error> {
        poll(&mut cancelled)?;
        let n = self.face.gas().species_count();
        if !positive(duration_s)
            || parcels.is_empty()
            || parcels.len() > MAX_SNAPSHOT_PARCELS
            || !parcels.keys().eq(endpoints.keys())
            || links
                .len()
                .checked_mul(n + 1)
                .is_none_or(|m| m > MAX_SNAPSHOT_TRANSFERS)
        {
            return Err(Error::Invalid(
                "bounded graph, exact endpoint map and positive duration required",
            ));
        }
        let mut ordered = BTreeMap::new();
        let mut traces = BTreeMap::new();
        for link in links {
            poll(&mut cancelled)?;
            if link.parcels[0] == link.parcels[1]
                || !positive(link.area_m2)
                || !positive(link.distance_m)
                || ordered.insert(link.id, link).is_some()
            {
                return Err(Error::Invalid(
                    "unique continuum face IDs, distinct gas owners and positive geometry required",
                ));
            }
            for &id in &link.parcels {
                let parcel = parcels
                    .get(&id)
                    .ok_or(Error::Invalid("missing continuum gas owner"))?;
                if let std::collections::btree_map::Entry::Vacant(entry) = traces.entry(id) {
                    entry.insert(self.partial_densities(parcel)?);
                }
            }
        }
        for (&id, partial) in &traces {
            poll(&mut cancelled)?;
            self.face.qualify_primitive(
                primitive(&parcels[&id], partial),
                controls.face,
                &mut cancelled,
            )?;
        }
        let mut transfers = Vec::new();
        let mut heat = Vec::with_capacity(links.len());
        let mut reports = Vec::with_capacity(links.len());
        let mut flux = vec![0.0; n + 3];
        for (id, link) in ordered {
            poll(&mut cancelled)?;
            let [left, right] = link.parcels;
            let trace = |id| primitive(&parcels[&id], &traces[&id]);
            let face = self.face.face_from_primitives(
                trace(left),
                trace(right),
                link.distance_m,
                [0.0; 2],
                controls.face,
                &mut flux,
                &mut cancelled,
            )?;
            let area_time = link.area_m2 * duration_s;
            if !positive(area_time) {
                return Err(Error::Invalid("face area-time is unrepresentable"));
            }
            let mut correction = vec![face.reduced_heat_flux_w_m2 * area_time];
            let mut masses = Vec::with_capacity(n);
            for (source, &rate) in self.face.gas().source().models().iter().zip(&flux[..n]) {
                poll(&mut cancelled)?;
                let mass = rate * area_time;
                if !mass.is_finite() || (rate != 0.0 && mass == 0.0) {
                    return Err(Error::Invalid(
                        "integrated continuum species flux is unrepresentable",
                    ));
                }
                masses.push(mass);
                if mass == 0.0 {
                    continue;
                }
                let (donor, receiver) = if mass > 0.0 {
                    (left, right)
                } else {
                    (right, left)
                };
                let change = source
                    .molar_enthalpy_change(
                        Temperature::new(parcels[&donor].state().state().temperature),
                        Temperature::new(face.enthalpy_temperature_k),
                    )
                    .map_err(|e| {
                        Error::Face(DissipationError::Flow(FlowError::Thermochemistry(e)))
                    })?;
                correction.push(mass * (change.value() / source.molar_mass().value()));
                transfers.push(SpeciesNetworkTransfer {
                    donor,
                    receiver,
                    species: source.species().clone(),
                    mass_kg: mass.abs(),
                });
            }
            let correction = crate::material::compensated_sum(correction.into_iter());
            if !correction.is_finite() {
                return Err(Error::Invalid(
                    "integrated continuum heat correction is unrepresentable",
                ));
            }
            heat.push(SpeciesNetworkHeatTransfer {
                parcels: [left, right],
                heat_left_to_right_j: correction,
            });
            reports.push(ParcelMolecularLinkStep {
                id,
                parcels: [left, right],
                mass_left_to_right_kg: masses,
                heat_correction_left_to_right_j: correction,
                face,
            });
        }
        let candidate = transport::step_candidate_with_paths_and_fluxes(
            parcels,
            paths,
            duration_s,
            &transfers,
            &heat,
            endpoints,
            controls.physical,
            &mut cancelled,
        )?;
        for &id in traces.keys() {
            poll(&mut cancelled)?;
            let parcel = &candidate.parcels()[&id];
            let partial = self.partial_densities(parcel)?;
            self.face.qualify_primitive(
                primitive(parcel, &partial),
                controls.face,
                &mut cancelled,
            )?;
        }
        poll(&mut cancelled)?;
        let (parcels, transport) = candidate.into_parts();
        Ok(ParcelMolecularCandidate {
            parcels,
            step: ParcelMolecularStep {
                species: self
                    .face
                    .gas()
                    .source()
                    .models()
                    .iter()
                    .map(|s| s.species().clone())
                    .collect(),
                links: reports,
                transport,
            },
        })
    }

    fn partial_densities(&self, parcel: &ClosedGasParcel) -> Result<Vec<f64>, Error> {
        let active = parcel.model().thermodynamics();
        let sources = self.face.gas().source().models();
        let masses: Vec<_> = active.components().iter().map(|s| s.molar_mass()).collect();
        let fractions = active
            .mole_composition()
            .to_mass_fractions(&masses)
            .map_err(Error::Composition)?;
        let density = parcel.mass_kg() / parcel.volume_m3();
        let mut partial = vec![0.0; sources.len()];
        for (source, &fraction) in active.components().iter().zip(fractions.fractions()) {
            let Some(i) = sources.iter().position(|s| s.species() == source.species()) else {
                return Err(Error::SourceMismatch(source.species().clone()));
            };
            if source != &sources[i] {
                return Err(Error::SourceMismatch(source.species().clone()));
            }
            partial[i] = density * fraction;
            if !positive(partial[i]) {
                return Err(Error::Invalid(
                    "positive original species density is unrepresentable",
                ));
            }
        }
        Ok(partial)
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn primitive<'a>(parcel: &ClosedGasParcel, partial: &'a [f64]) -> NormalPrimitive<'a> {
    NormalPrimitive {
        partial_densities: partial,
        normal_velocity: 0.0,
        transverse_velocity: 0.0,
        pressure: parcel.state().state().pressure,
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
