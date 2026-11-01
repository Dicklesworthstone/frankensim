//! Held material baths on the original finite-gas species/enthalpy snapshot.
//!
//! The shared reservoir face owns mass-frame Maxwell--Stefan, pressure driving,
//! Soret and reciprocal Dufour/Fourier heat. Its bath enthalpy chart is lowered
//! analytically into the original material boundary owner. Every exposed bath
//! reads the same finite inventory; incoming mass cannot fund another outflow.
//! This is a first-order, zero-barycentric-velocity candidate, with no invented
//! finite bath capacity, bulk-flow law or independent caloric reconstruction.

use super::*;
use fs_matdb::QueryPoint;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::{
    self, SpeciesNetworkExternalHeat, SpeciesNetworkStep,
    boundary::{BoundarySpeciesTransfer, SpeciesNetworkBoundary},
    transport::{GasReservoirExchange, GasReservoirWindowAudit, ReservoirSpeciesTransfer},
};
use fs_material::thermochemical_gas::{ResolvedThermochemicalGasState, ThermochemicalGasModel};

/// An immutable, externally held material state with its actual source receipts.
/// Its thermodynamic composition belongs to `model`; its T, p and fixed external
/// context belong to `point`. No finite reservoir mass or volume is fabricated.
#[derive(Clone, Debug, PartialEq)]
pub struct MolecularReservoir {
    model: ThermochemicalGasModel,
    point: QueryPoint,
    state: ResolvedThermochemicalGasState,
}
impl MolecularReservoir {
    pub fn new(
        model: ThermochemicalGasModel,
        point: QueryPoint,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ParcelMolecularError> {
        poll(&mut cancelled)?;
        let state = model.evaluate(&point).map_err(DissipationError::Material)?;
        poll(&mut cancelled)?;
        Ok(Self {
            model,
            point,
            state,
        })
    }
    pub const fn model(&self) -> &ThermochemicalGasModel {
        &self.model
    }
    pub const fn point(&self) -> &QueryPoint {
        &self.point
    }
    pub const fn state(&self) -> &ResolvedThermochemicalGasState {
        &self.state
    }
}

/// Independently exposed bath boundary, oriented from reservoir into finite gas.
/// Geometry is the actual exposed area and cell-center-to-boundary distance.
#[derive(Clone, Copy, Debug)]
pub struct ReservoirMolecularLink<'a> {
    pub id: u64,
    pub parcel: u64,
    pub reservoir: &'a MolecularReservoir,
    pub area_m2: f64,
    pub distance_m: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ReservoirMolecularLinkStep {
    pub id: u64,
    pub parcel: u64,
    /// Signed mass into finite gas in `ReservoirMolecularStep::species` order.
    pub mass_into_gas_kg: Vec<f64>,
    /// Additional heat beyond the bath's original partial-enthalpy streams [J].
    pub reduced_heat_into_gas_j: f64,
    /// Actual held bath material/caloric receipt used for this face and audit.
    pub reservoir: ResolvedThermochemicalGasState,
    pub face: TransportFaceReport,
}
#[derive(Clone, Debug, PartialEq)]
pub struct ReservoirMolecularStep {
    pub duration_s: f64,
    pub species: Vec<SpeciesId>,
    pub links: Vec<ReservoirMolecularLinkStep>,
    /// Original gross withdrawal, source join and one-shot caloric reconstruction.
    pub exchange: SpeciesNetworkStep,
    /// Actual gas endpoint changes minus all held-bath mass, enthalpy and heat.
    pub audit: GasReservoirWindowAudit,
}
#[derive(Clone, Debug, PartialEq)]
pub struct ReservoirMolecularCandidate {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    step: ReservoirMolecularStep,
}
impl ReservoirMolecularCandidate {
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub const fn step(&self) -> &ReservoirMolecularStep {
        &self.step
    }
    pub fn into_parts(self) -> (BTreeMap<u64, ClosedGasParcel>, ReservoirMolecularStep) {
        (self.parcels, self.step)
    }
}

impl<L: TransportLaw> ParcelMolecularTransport<L> {
    /// Sample all held-bath faces from one unchanged finite-gas inventory.
    ///
    /// Both directions of a species stream carry the boundary bath enthalpy,
    /// as prescribed by the original reservoir face. The analytic sensible
    /// correction and separate reduced heat reach the SAME native inversion.
    /// The final audit includes bath entropy flow; gas entropy need not rise.
    /// A nonzero stream at zero bath activity refuses, without adding a trace.
    /// Initial gas/bath, face and final gas states must all satisfy the original
    /// constitutive law. Cancellation or any refusal publishes no partial state.
    #[allow(clippy::too_many_arguments)]
    pub fn reservoir_candidate(
        &mut self,
        parcels: &BTreeMap<u64, ClosedGasParcel>,
        links: &[ReservoirMolecularLink<'_>],
        duration_s: f64,
        endpoints: &BTreeMap<u64, GasEnergySolve>,
        controls: ParcelMolecularSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReservoirMolecularCandidate, ParcelMolecularError> {
        poll(&mut cancelled)?;
        let n = self.face.gas().species_count();
        if !positive(duration_s)
            || parcels.is_empty()
            || parcels.len() > MAX_SNAPSHOT_PARCELS
            || !parcels.keys().eq(endpoints.keys())
            || links
                .len()
                .checked_mul(n + 2)
                .is_none_or(|m| m > MAX_SNAPSHOT_TRANSFERS)
        {
            return Err(Error::Invalid(
                "bounded reservoir graph, exact endpoint map and positive duration required",
            ));
        }
        let mut ordered = BTreeMap::new();
        let mut traces = BTreeMap::new();
        // Local immutable identities avoid re-resolving a shared bath per face.
        let mut baths: BTreeMap<*const MolecularReservoir, Vec<f64>> = BTreeMap::new();
        for link in links {
            poll(&mut cancelled)?;
            if !positive(link.area_m2)
                || !positive(link.distance_m)
                || ordered.insert(link.id, link).is_some()
            {
                return Err(Error::Invalid(
                    "reservoir faces need unique IDs and positive geometry",
                ));
            }
            let gas = parcels
                .get(&link.parcel)
                .ok_or(Error::Invalid("missing reservoir-connected gas owner"))?;
            if let std::collections::btree_map::Entry::Vacant(entry) = traces.entry(link.parcel) {
                entry.insert(self.partial_densities(gas)?);
            }
            let key = core::ptr::from_ref(link.reservoir);
            if let std::collections::btree_map::Entry::Vacant(entry) = baths.entry(key) {
                let partial = self.reservoir_partial_densities(link.reservoir)?;
                self.face.qualify_primitive(
                    reservoir_primitive(link.reservoir, &partial),
                    controls.face,
                    &mut cancelled,
                )?;
                entry.insert(partial);
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
        let mut species_rows = Vec::with_capacity(links.len());
        let mut audit_rows = Vec::with_capacity(links.len());
        let mut heat = Vec::with_capacity(links.len());
        let mut reports = Vec::with_capacity(links.len());
        let mut flux = vec![0.0; n + 3];
        for (&id, link) in &ordered {
            poll(&mut cancelled)?;
            let bath = link.reservoir;
            let partial = &baths[&core::ptr::from_ref(bath)];
            let face = self.face.reservoir_face_from_primitives(
                reservoir_primitive(bath, partial),
                primitive(&parcels[&link.parcel], &traces[&link.parcel]),
                link.distance_m,
                [0.0; 2],
                true,
                controls.face,
                &mut flux,
                &mut cancelled,
            )?;
            let area_time = link.area_m2 * duration_s;
            if !positive(area_time) {
                return Err(Error::Invalid(
                    "reservoir face area-time is unrepresentable",
                ));
            }
            let q = face.reduced_heat_flux_w_m2 * area_time;
            if !q.is_finite() || (q == 0.0 && face.reduced_heat_flux_w_m2 != 0.0) {
                return Err(Error::Invalid(
                    "integrated reservoir heat is unrepresentable",
                ));
            }
            let mut rows = Vec::new();
            let mut audit = Vec::new();
            let mut masses = Vec::with_capacity(n);
            for (i, source) in self.face.gas().source().models().iter().enumerate() {
                poll(&mut cancelled)?;
                let dm = flux[i] * area_time;
                if !dm.is_finite() || (flux[i] != 0.0 && dm == 0.0) {
                    return Err(Error::Invalid(
                        "integrated reservoir species flux is unrepresentable",
                    ));
                }
                masses.push(dm);
                if dm == 0.0 {
                    continue;
                }
                if partial[i] == 0.0 {
                    return Err(Error::Invalid(
                        "nonzero reservoir stream requires positive bath species activity",
                    ));
                }
                let dh = source
                    .molar_enthalpy_change(
                        Temperature::new(parcels[&link.parcel].state().state().temperature),
                        Temperature::new(face.enthalpy_temperature_k),
                    )
                    .map_err(|e| {
                        Error::Face(DissipationError::Flow(FlowError::Thermochemistry(e)))
                    })?
                    .value()
                    / source.molar_mass().value();
                let relative_energy = dm * dh;
                if !relative_energy.is_finite() {
                    return Err(Error::Invalid(
                        "reservoir sensible enthalpy is unrepresentable",
                    ));
                }
                rows.push(BoundarySpeciesTransfer {
                    species: source.species().clone(),
                    incoming_mass_kg: dm.max(0.0),
                    outgoing_mass_kg: (-dm).max(0.0),
                    net_mass_into_gas_kg: dm,
                    relative_energy_into_gas_j: relative_energy,
                });
                audit.push(ReservoirSpeciesTransfer {
                    species: source.species().clone(),
                    mass_into_gas_kg: dm,
                });
            }
            species_rows.push(rows);
            audit_rows.push(audit);
            heat.push(SpeciesNetworkExternalHeat {
                id,
                parcel: link.parcel,
                heat_into_gas_j: q,
            });
            reports.push(ReservoirMolecularLinkStep {
                id,
                parcel: link.parcel,
                mass_into_gas_kg: masses,
                reduced_heat_into_gas_j: q,
                reservoir: bath.state.clone(),
                face,
            });
        }
        let boundaries: Vec<_> = ordered
            .values()
            .zip(&species_rows)
            .filter(|(_, rows)| !rows.is_empty())
            .map(|(link, rows)| SpeciesNetworkBoundary {
                id: link.id,
                parcel: link.parcel,
                catalog: &link.reservoir.model,
                species: rows,
            })
            .collect();
        let candidate = network::exchange_candidate_with_external_heat(
            parcels,
            &[],
            &[],
            &boundaries,
            &heat,
            endpoints,
            controls.physical.exchange,
            &mut cancelled,
        )
        .map_err(GasTransportError::from)?;
        for &id in traces.keys() {
            poll(&mut cancelled)?;
            let gas = &candidate.parcels()[&id];
            let partial = self.partial_densities(gas)?;
            self.face
                .qualify_primitive(primitive(gas, &partial), controls.face, &mut cancelled)?;
        }
        let exchanges: Vec<_> = ordered
            .values()
            .zip(&audit_rows)
            .zip(&heat)
            .map(|((link, rows), q)| GasReservoirExchange {
                parcel: link.parcel,
                catalog: &link.reservoir.model,
                state: &link.reservoir.state,
                species: rows,
                reduced_heat_into_gas_j: q.heat_into_gas_j,
            })
            .collect();
        let audit = transport::audit_reservoir_window(
            parcels,
            candidate.parcels(),
            &exchanges,
            controls.physical,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        let (parcels, exchange) = candidate.into_parts();
        Ok(ReservoirMolecularCandidate {
            parcels,
            step: ReservoirMolecularStep {
                duration_s,
                species: self
                    .face
                    .gas()
                    .source()
                    .models()
                    .iter()
                    .map(|s| s.species().clone())
                    .collect(),
                links: reports,
                exchange,
                audit,
            },
        })
    }

    fn reservoir_partial_densities(&self, bath: &MolecularReservoir) -> Result<Vec<f64>, Error> {
        let active = bath.model.thermodynamics();
        let sources = self.face.gas().source().models();
        let masses: Vec<_> = active.components().iter().map(|s| s.molar_mass()).collect();
        let fractions = active
            .mole_composition()
            .to_mass_fractions(&masses)
            .map_err(Error::Composition)?;
        let mut partial = vec![0.0; sources.len()];
        for (source, &fraction) in active.components().iter().zip(fractions.fractions()) {
            let Some(i) = sources.iter().position(|s| s.species() == source.species()) else {
                return Err(Error::SourceMismatch(source.species().clone()));
            };
            if source != &sources[i] {
                return Err(Error::SourceMismatch(source.species().clone()));
            }
            partial[i] = bath.state.state().density * fraction;
            if !positive(partial[i]) {
                return Err(Error::Invalid(
                    "positive reservoir species density is unrepresentable",
                ));
            }
        }
        Ok(partial)
    }
}

fn reservoir_primitive<'a>(bath: &MolecularReservoir, partial: &'a [f64]) -> NormalPrimitive<'a> {
    NormalPrimitive {
        partial_densities: partial,
        normal_velocity: 0.0,
        transverse_velocity: 0.0,
        pressure: bath.state.state().pressure,
    }
}
