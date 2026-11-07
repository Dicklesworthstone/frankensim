//! Original conductivity and convection on the retained gas/wall heat timeline.
//!
//! Patch area, effective distances and correlation geometry are explicitly held
//! thermal-network geometry. They are not inferred from the deformed mesh.
//! Every flux uses the actual heated gas endpoint and current native cell phase.
//! The existing contact owner supplies source admission, series conduction,
//! surface balance, integrated k(T), correlation regimes and their evidence.

use super::{
    CavityHeatError, CavityHeatObserverError, ThermalCavityGas, ThermalPorousGasSystem, number,
    poll,
};
use crate::solid_enthalpy::SolidEnthalpyError;
use crate::solid_enthalpy::finite_plastic::FiniteThermoplasticQuery;
use crate::solid_enthalpy::gas::convection::endpoint_law;
use crate::solid_enthalpy::gas::{
    PreparedContact, SolidGasContact, SolidGasConvectionSample, SolidGasError, StructuralSource,
};
use fs_exec::Cx;
use fs_material::gas_species::input::temperature::IsobaricGasLaw;
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint};
use fs_material::thermal::conduction::source::ResolvedConductionLayer;
use fs_material::thermal::{PureSolidThermalLaw, source::ThermalMaterialPoint};
use fs_material::thermochemical_gas::{
    ResolvedThermochemicalGasState, energy::parcel::ClosedGasParcel,
};
use fs_solid::hyper3::pressure::cavity::CavityGasEpoch;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

/// One original source-bound patch and the finite gas that exchanges its heat.
///
/// The complete `SolidGasContact` retains its original layer/convection scope.
/// All patches sharing a cell must fit its exposed reference-network area.
/// Fixed W/K links may coexist as separately authored parallel heat paths;
/// their physical coverage is not inferred from their conductance.
#[derive(Clone, Debug)]
pub struct CavityWallSourceContact {
    pub cavity: usize,
    pub contact: SolidGasContact,
}

/// Actual integrated transfer and last applied native source query of a patch.
#[derive(Clone, Debug)]
pub struct CavityWallSourceReceipt {
    pub id: u64,
    pub cavity: usize,
    pub cell: CellId,
    pub heat_into_solid_j: f64,
    pub last_cell_temperature_k: f64,
    pub last_power_into_solid_w: f64,
    /// Original series diagnostic; actual film power remains authoritative.
    pub last_conductance_w_k: f64,
    /// Nonnegative explicit flux secant, not an inferred nonlinear tangent.
    pub last_explicit_decay_conductance_w_k: f64,
    pub source: ResolvedMaterialStatePoint,
    pub convection: Option<SolidGasConvectionSample>,
}

struct SourceCell<'m> {
    cell: CellId,
    query: FiniteThermoplasticQuery<'m>,
    reference_h_j_kg: f64,
}
impl SourceCell<'_> {
    fn point(&self) -> ThermalMaterialPoint<'_> {
        ThermalMaterialPoint {
            material: self.query.card,
            context: &self.query.context,
            temperature_axis: &self.query.temperature_axis,
            // Finite caloric admission already requires the original all-solid
            // chart. No extra phase coordinate or mixture rule is introduced.
            liquid_fraction_axis: None,
        }
    }
    fn binding(&self) -> StructuralSource<'_> {
        StructuralSource {
            cell: self.cell,
            curve: self.query.caloric.curve(),
            point: self.point(),
            reference_h_j_kg: self.reference_h_j_kg,
        }
    }
}

struct Film<'m> {
    route: CavityWallSourceContact,
    element: usize,
    source: usize,
    integrated_solid: Option<Arc<ResolvedConductionLayer<'m>>>,
}
impl<'m> Film<'m> {
    fn prepared(&self) -> PreparedContact<'_, 'm> {
        PreparedContact::new(
            &self.route.contact,
            self.element,
            self.integrated_solid.clone(),
        )
    }
}

#[derive(Clone)]
struct GasEndpoint<'g> {
    gas: &'g ResolvedThermochemicalGasState,
    law: Option<Arc<IsobaricGasLaw>>,
}

// The callback only borrows an original physical endpoint. Keeping it with the
// immutable preparation leaves the generic ideal/real fixed-link path intact;
// this source-contact constructor admits only the existing ideal film owner.
pub(super) struct PreparedContacts<'m, G> {
    cells: Vec<SourceCell<'m>>,
    films: Vec<Film<'m>>,
    areas: BTreeMap<CellId, f64>,
    needs_convection: Vec<bool>,
    endpoint: for<'g> fn(&'g G, bool) -> Result<GasEndpoint<'g>, SolidGasError>,
}

fn ideal_endpoint(
    gas: &ClosedGasParcel,
    convection: bool,
) -> Result<GasEndpoint<'_>, SolidGasError> {
    Ok(GasEndpoint {
        gas: gas.state(),
        law: if convection {
            Some(Arc::new(endpoint_law(
                gas.model(),
                gas.point(),
                gas.state(),
            )?))
        } else {
            None
        },
    })
}

impl<'m> ThermalPorousGasSystem<'m, ClosedGasParcel> {
    /// Add source-owned gas/solid conduction or convection to this retained owner.
    ///
    /// Conductivity is selected at the original stress-free material reference
    /// and remains pinned through all accepted increments. Integrated pure-phase
    /// profiles are compiled once. Current gas pressure/composition and current
    /// cell temperature feed the original source law on every heat proposal.
    ///
    /// Configure once; unique patch IDs share the 16,384-link budget with fixed
    /// W/K links. Geometry stays explicitly held, and ambient cell convection
    /// must be absent on source-contact cells. Existing correlation/phase and
    /// homogeneous wall assumptions remain unchanged. Dense real-gas source
    /// convection needs its own original property-law support and is not inferred.
    pub fn with_source_contacts(
        mut self,
        thermal: &DepositionSimulation<'_>,
        mut contacts: Vec<CavityWallSourceContact>,
        cx: &Cx<'_>,
    ) -> Result<Self, CavityHeatError> {
        poll(cx)?;
        if self.source_contacts.is_some()
            || contacts.len() > 16_384 - self.links.len()
            || thermal.time_s().to_bits() != self.time_s().to_bits()
        {
            return Err(CavityHeatError::Invalid(
                "source contacts require one bounded configuration at the accepted clock",
            ));
        }
        self.solid
            .binding()
            .enthalpies::<std::convert::Infallible>(thermal, true, cx)
            .map_err(CavityHeatError::Binding)?;
        let elements: BTreeMap<_, _> = self
            .solid
            .cells()
            .iter()
            .enumerate()
            .map(|(element, cell)| (cell.id, element))
            .collect();
        let mut used: BTreeSet<_> = self.links.iter().map(|link| link.id).collect();
        let mut cell_indices = BTreeMap::new();
        let mut prepared = PreparedContacts {
            cells: Vec::new(),
            films: Vec::with_capacity(contacts.len()),
            areas: BTreeMap::new(),
            needs_convection: vec![false; self.cavities.len()],
            endpoint: ideal_endpoint,
        };
        contacts.sort_by_key(|entry| entry.contact.id);
        for route in contacts {
            poll(cx)?;
            let contact = &route.contact;
            let element = *elements.get(&contact.cell).ok_or(CavityHeatError::Invalid(
                "source contact must name an original structural cell",
            ))?;
            if route.cavity >= self.cavities.len()
                || !used.insert(contact.id)
                || !contact.geometry.area_m2.is_finite()
                || contact.geometry.area_m2 < 0.0
            {
                return Err(CavityHeatError::Invalid(
                    "source contacts need unique IDs, an actual gas and finite nonnegative area",
                ));
            }
            let area = prepared.areas.entry(contact.cell).or_insert(0.0);
            *area = number(*area + contact.geometry.area_m2)?;
            let cell = &thermal.cells()[&contact.cell];
            if *area
                > thermal
                    .exposed_area_m2(contact.cell)
                    .map_err(SolidGasError::from)
                    .map_err(CavityHeatError::Source)?
                || cell.boundary.htc_w_m2_k != 0.0
            {
                return Err(CavityHeatError::Invalid(
                    "source contacts exceed exposed area or duplicate ambient convection",
                ));
            }
            let source = match cell_indices.get(&contact.cell) {
                Some(&index) => index,
                None => {
                    let query = self.solid.mechanics().families()[element].query().clone();
                    let temperature = query.context.axes()[&query.temperature_axis];
                    let caloric = PureSolidThermalLaw::new(query.caloric.curve())
                        .map_err(SolidEnthalpyError::from)
                        .map_err(SolidGasError::from)
                        .map_err(CavityHeatError::Source)?;
                    let reference_h_j_kg = caloric
                        .specific_enthalpy_at_temperature(temperature)
                        .map_err(SolidEnthalpyError::from)
                        .map_err(SolidGasError::from)
                        .map_err(CavityHeatError::Source)?;
                    let index = prepared.cells.len();
                    prepared.cells.push(SourceCell {
                        cell: contact.cell,
                        query,
                        reference_h_j_kg,
                    });
                    cell_indices.insert(contact.cell, index);
                    index
                }
            };
            prepared.needs_convection[route.cavity] |= contact.convection.is_some();
            let mut film = Film {
                route,
                element,
                source,
                integrated_solid: None,
            };
            let bound = &prepared.cells[source];
            let reference = bound
                .query
                .caloric
                .curve()
                .state_at_specific_enthalpy(bound.reference_h_j_kg)
                .map_err(SolidEnthalpyError::from)
                .map_err(SolidGasError::from)
                .map_err(CavityHeatError::Source)?;
            let selected = film
                .prepared()
                .reference_material(
                    &bound.binding(),
                    self.cavities[film.route.cavity].gas.state(),
                    reference,
                )
                .map_err(CavityHeatError::Source)?;
            let pin = selected
                .property_by_key(&film.route.contact.conductivity)
                .ok_or(CavityHeatError::Invalid(
                    "source contact lacks selected conductivity",
                ))?
                .answer()
                .receipt
                .selected;
            film.route.contact.selection = MaterialPropertySelection::PinnedByProperty(vec![(
                film.route.contact.conductivity.name().to_owned(),
                pin,
            )]);
            film.integrated_solid = film
                .prepared()
                .prepare_solid_profile_at(bound.query.caloric.curve(), bound.point(), &mut || {
                    cx.checkpoint().map_err(|_| DepositionError::Cancelled)
                })
                .map_err(CavityHeatError::Source)?;
            prepared.films.push(film);
        }
        // Admit each CURRENT gas/cell endpoint, not a convection regime at the
        // unrelated stress-free solid temperature used only for source pinning.
        let endpoints = self
            .cavities
            .iter()
            .zip(&prepared.needs_convection)
            .map(|(entry, &needed)| {
                cx.checkpoint().map_err(|_| SolidGasError::Cancelled)?;
                ideal_endpoint(&entry.gas, needed)
            })
            .collect::<Result<Vec<_>, _>>()
            .map_err(CavityHeatError::Source)?;
        SourceFlux {
            prepared: &prepared,
            endpoints,
        }
        .samples(thermal, &mut || {
            cx.checkpoint().map_err(|_| DepositionError::Cancelled)
        })
        .map_err(CavityHeatError::Source)?;
        poll(cx)?;
        self.source_contacts = Some(Arc::new(prepared));
        Ok(self)
    }
}

// Borrow each actual physical endpoint and share its property law through all
// cloned thermal probes. Source preparation remains immutable.
#[derive(Clone)]
pub(super) struct SourceFlux<'s, 'm, G> {
    prepared: &'s PreparedContacts<'m, G>,
    endpoints: Vec<GasEndpoint<'s>>,
}

impl<'m, G: ThermalCavityGas> PreparedContacts<'m, G> {
    pub(super) fn flux<'s>(
        &'s self,
        gases: &'s [CavityGasEpoch<G>],
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<SourceFlux<'s, 'm, G>, SolidGasError> {
        let endpoints = gases
            .iter()
            .zip(&self.needs_convection)
            .map(|(gas, &needed)| {
                checkpoint()?;
                (self.endpoint)(gas.gas(), needed)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(SourceFlux {
            prepared: self,
            endpoints,
        })
    }
}

#[derive(Clone)]
pub(super) struct SourceSample {
    pub receipt: CavityWallSourceReceipt,
}

impl<G> SourceFlux<'_, '_, G> {
    pub(super) fn len(&self) -> usize {
        self.prepared.films.len()
    }

    pub(super) fn samples(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<SourceSample>, SolidGasError> {
        for (&cell, &area) in &self.prepared.areas {
            checkpoint()?;
            let bound = thermal
                .cells()
                .get(&cell)
                .ok_or(SolidGasError::Invalid("inactive source-contact cell"))?;
            if area > thermal.exposed_area_m2(cell)? || bound.boundary.htc_w_m2_k != 0.0 {
                return Err(SolidGasError::Invalid(
                    "source-contact coverage changed or duplicates ambient convection",
                ));
            }
        }
        let mut samples = Vec::with_capacity(self.len());
        for film in &self.prepared.films {
            checkpoint()?;
            let source = &self.prepared.cells[film.source];
            let phase = source
                .query
                .caloric
                .curve()
                .state_at_specific_enthalpy(thermal.cells()[&source.cell].specific_h_j_kg)
                .map_err(SolidEnthalpyError::from)?;
            let endpoint = &self.endpoints[film.route.cavity];
            let sample = film.prepared().resolve_with_law(
                &source.binding(),
                endpoint.gas,
                phase,
                None,
                endpoint.law.as_deref(),
                checkpoint,
            )?;
            let power = -sample.power_to_gas_w;
            let contrast = endpoint.gas.state().temperature - phase.temperature_k();
            let decay = if sample.convection.is_some() {
                if contrast == 0.0 {
                    if power != 0.0 {
                        return Err(SolidGasError::Invalid(
                            "equilibrium source film carries heat",
                        ));
                    }
                    0.0
                } else {
                    let value = power / contrast;
                    if !value.is_finite() || value < 0.0 || (power != 0.0 && value == 0.0) {
                        return Err(SolidGasError::Invalid(
                            "source-film explicit decay is unrepresentable",
                        ));
                    }
                    value
                }
            } else {
                sample.conductance_w_k
            };
            if !power.is_finite() || !decay.is_finite() || decay < 0.0 {
                return Err(SolidGasError::Invalid(
                    "source-contact flux is unrepresentable",
                ));
            }
            if sample.convection.is_none() && decay > 0.0 && contrast != 0.0 && power == 0.0 {
                return Err(SolidGasError::Invalid(
                    "source-contact heat power underflows",
                ));
            }
            samples.push(SourceSample {
                receipt: CavityWallSourceReceipt {
                    id: film.route.contact.id,
                    cavity: film.route.cavity,
                    cell: source.cell,
                    heat_into_solid_j: 0.0,
                    last_cell_temperature_k: phase.temperature_k(),
                    last_power_into_solid_w: power,
                    last_conductance_w_k: sample.conductance_w_k,
                    last_explicit_decay_conductance_w_k: decay,
                    source: sample.source,
                    convection: sample.convection,
                },
            });
        }
        Ok(samples)
    }
}

pub(super) fn observer_error<E>(error: SolidGasError) -> CavityHeatObserverError<E> {
    if matches!(error, SolidGasError::Cancelled) {
        CavityHeatObserverError::Thermal(DepositionError::Cancelled)
    } else {
        CavityHeatObserverError::Source(error)
    }
}
