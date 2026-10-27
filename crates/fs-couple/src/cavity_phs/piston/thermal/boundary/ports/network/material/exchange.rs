//! Electrical losses, contact conduction and enclosure radiation on one state.
//! This composes existing material/heat owners; it introduces no transfer law,
//! nonlinear solver, thermal history, or material-name-based default.
/// Pluggable gas heat laws on this same complete physical assembly.
pub mod gas_law;
/// Same-endpoint reacting chambers on this original electrical/thermal assembly.
pub mod reactive;

use super::MaterialContact;
use super::super::{Chamber, DrivenNetworkError, NetworkCandidate, NetworkError, NetworkSolve};
use crate::cavity_phs::piston::thermal::boundary::network::material::resolve_contact;
use crate::circuit_enthalpy::{CircuitEnthalpyError, CircuitEnthalpyInput, CircuitThermalCoefficients};
use crate::circuit_enthalpy::contacts::{ContactThermalCoefficients, ContactThermalSample};
use crate::circuit_enthalpy::contacts::material::{
    ElectricalContactNetwork, ElectricalContactNetworkError, PreparedElectricalContactNetwork,
};
use crate::circuit_enthalpy::material::{MaterialNetworkError, MaterialNetworkSample};
use crate::material_circuit::MaterialCircuit;
use fs_conduction::ConductionError;
use fs_conduction::transient::deposition::contact::{
    apply_material_contacts, ContactLinkSample, MaterialContactLink, MAX_MATERIAL_NETWORK_CONTACTS,
};
use fs_conduction::transient::deposition::radiation::{
    MaterialRadiation, RadiatingCell, RadiationCellError, RadiationCellSources,
};
use fs_exec::Cx;
use fs_matdb::ClaimSelection;
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};

/// Data-only physical assembly. Empty contact lists or absent radiation mean
/// explicitly omitted mechanisms, not inferred perfect contact or black surfaces.
pub struct MaterialExchangeNetwork<'a> {
    /// Existing bulk transport, electrical contacts and their heat partitions.
    pub electrical: &'a ElectricalContactNetwork<'a>,
    /// Contact resistance added in SERIES to the original bulk thermal links.
    pub thermal_contacts: &'a [MaterialContactLink<'a>],
    /// Optional closed opaque gray-diffuse enclosure on those same finite cells.
    /// Gas is assumed nonparticipating in this radiation model.
    pub radiation: Option<&'a MaterialRadiation<'a>>,
}

/// Original source failures, without substituting missing material properties.
#[derive(Debug)]
pub enum MaterialExchangeError {
    /// Inconsistent physical source bindings or unrepresentable transfer.
    Invalid(&'static str),
    /// Original electrical/bulk material owner refusal.
    Electrical(ElectricalContactNetworkError),
    /// Original contact-resistance/series-law refusal.
    Contact(ConductionError),
    /// Original emissivity/radiosity/thermal refusal.
    Radiation(RadiationCellError),
    /// Cancellation of the caller's execution scope.
    Cancelled,
}
impl core::fmt::Display for MaterialExchangeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material exchange network: {self:?}")
    }
}
impl std::error::Error for MaterialExchangeError {}
impl From<ElectricalContactNetworkError> for MaterialExchangeError {
    fn from(e: ElectricalContactNetworkError) -> Self {
        match e {
            ElectricalContactNetworkError::Thermal(DepositionError::Cancelled)
            | ElectricalContactNetworkError::Bulk(MaterialNetworkError::Thermal(DepositionError::Cancelled)) => Self::Cancelled,
            other => Self::Electrical(other),
        }
    }
}
impl From<ConductionError> for MaterialExchangeError {
    fn from(e: ConductionError) -> Self {
        match e { ConductionError::Cancelled { .. } => Self::Cancelled, other => Self::Contact(other) }
    }
}
impl From<RadiationCellError> for MaterialExchangeError {
    fn from(e: RadiationCellError) -> Self {
        match e { RadiationCellError::Cancelled => Self::Cancelled, other => Self::Radiation(other) }
    }
}
fn checkpoint(cx: &Cx<'_>) -> Result<(), MaterialExchangeError> {
    cx.checkpoint().map_err(|_| MaterialExchangeError::Cancelled)
}

/// Actual trial material payloads. Bulk link samples remain BULK-ONLY; the
/// contact samples retain the series-corrected conductances used by the residual.
#[derive(Debug)]
pub struct MaterialExchangeSample {
    /// Original complete electrical/thermal snapshots and bulk link sources.
    pub bulk: MaterialNetworkSample,
    /// Original contact-property receipts, phases, area and series reductions.
    pub contacts: Vec<ContactLinkSample>,
    /// Actual surface emissivities, signed cell heat and original radiosity solve.
    pub radiation: Option<RadiationCellSources>,
}

/// Selected source identities, never cached physical coefficients or history.
/// Keep across accepted steps when continuity of the selected laws is required.
pub struct PreparedMaterialExchangeNetwork<'n, 'a> {
    binding: &'n MaterialExchangeNetwork<'a>,
    electrical: PreparedElectricalContactNetwork<'n, 'a>,
    thermal_contacts: Vec<MaterialContactLink<'a>>,
    radiating_cells: Vec<RadiatingCell<'a>>,
}

impl<'a> MaterialExchangeNetwork<'a> {
    /// Prepare ALL active mechanisms before any physical step. Each original
    /// resolver performs its own source, geometry and phase admission. Contacts
    /// shared by electrical and thermal transport must describe the same ordered
    /// interface card and conditions; their distinct property claims stay distinct.
    pub fn prepare<'n>(&'n self, cx: &Cx<'_>, thermal: &DepositionSimulation<'_>)
        -> Result<PreparedMaterialExchangeNetwork<'n, 'a>, MaterialExchangeError>
    {
        checkpoint(cx)?;
        let electrical = self.electrical.prepare(thermal, || cx.checkpoint().is_err())?;
        let mut initial = electrical.resolve(thermal, || cx.checkpoint().is_err())?;
        self.admit_common_states(cx, &initial.bulk.source)?;
        let contacts = apply_material_contacts(thermal, self.thermal_contacts,
            &mut initial.bulk.conductances_w_k, || cx.checkpoint().is_err())?;
        // The original owner has admitted bounded, unique active links. Preserve
        // caller order; pin by physical link ID, not incidental vector position.
        let mut thermal_contacts = Vec::with_capacity(self.thermal_contacts.len());
        for &binding in self.thermal_contacts {
            checkpoint(cx)?;
            let sample = contacts.iter().find(|s| s.link == binding.link)
                .ok_or(MaterialExchangeError::Invalid("admitted thermal contact sample missing"))?;
            thermal_contacts.push(MaterialContactLink {
                selection: ClaimSelection::Pinned(sample.interface.receipt().selected), ..binding
            });
        }
        let mut radiating_cells = Vec::new();
        if let Some(radiation) = self.radiation {
            let initial = radiation.resolve(cx, thermal)?;
            for binding in radiation.cells {
                checkpoint(cx)?;
                let sample = &initial.surfaces()[binding.surface];
                let pins = sample.material.properties().iter().map(|p|
                    (p.requirement().name().to_owned(), p.answer().receipt.selected)).collect();
                radiating_cells.push(RadiatingCell {
                    surface: binding.surface, cell: binding.cell, material_id: binding.material_id,
                    curve: binding.curve, material: binding.material,
                    selection: MaterialPropertySelection::PinnedByProperty(pins),
                });
            }
        }
        checkpoint(cx)?;
        Ok(PreparedMaterialExchangeNetwork { binding: self, electrical, thermal_contacts, radiating_cells })
    }

    fn admit_common_states(&self, cx: &Cx<'_>, samples: &MaterialNetworkSample) -> Result<(), MaterialExchangeError> {
        let invalid = MaterialExchangeError::Invalid;
        checkpoint(cx)?;
        if self.thermal_contacts.len() > MAX_MATERIAL_NETWORK_CONTACTS
            || self.radiation.is_some_and(|r| r.cells.len() != r.enclosure.surfaces().len()) {
            return Err(invalid("complete bounded thermal contact and radiative surface bindings required"));
        }
        for contact in self.thermal_contacts {
            checkpoint(cx)?;
            if let Some(electrical) = self.electrical.contacts.iter().find(|e| e.link == contact.link) {
                if contact.point.interface.content_hash() != electrical.point.interface.content_hash() {
                    return Err(invalid("electrical and thermal contact need the same ordered interface card"));
                }
                let a = samples.cells.iter().find(|c| c.cell == contact.sides[0].cell)
                    .ok_or(invalid("thermal contact first cell is absent"))?.phase;
                let b = samples.cells.iter().find(|c| c.cell == contact.sides[1].cell)
                    .ok_or(invalid("thermal contact second cell is absent"))?.phase;
                let left = contact.point.at([a,b]).map_err(|e|
                    MaterialExchangeError::Electrical(ElectricalContactNetworkError::State(e)))?;
                let right = electrical.point.at([a,b]).map_err(|e|
                    MaterialExchangeError::Electrical(ElectricalContactNetworkError::State(e)))?;
                if left != right { return Err(invalid("electrical and thermal interface conditions differ")); }
            }
        }
        if let Some(radiation) = self.radiation {
            for surface in radiation.cells {
                checkpoint(cx)?;
                let bulk = self.electrical.bulk.cells.iter().find(|c| c.cell == surface.cell)
                    .ok_or(invalid("radiating surface lacks its bulk material binding"))?;
                let phase = samples.cells.iter().find(|c| c.cell == surface.cell)
                    .ok_or(invalid("radiating surface lacks its actual bulk phase"))?.phase;
                let state_error = |e| MaterialExchangeError::Electrical(ElectricalContactNetworkError::State(e));
                if surface.material.material.content_hash() != bulk.point.material.content_hash()
                    || surface.material.at(phase).map_err(state_error)? != bulk.point.at(phase).map_err(state_error)? {
                    return Err(invalid("radiation must retain the exact bulk card and physical conditions"));
                }
            }
        }
        Ok(())
    }
}

impl PreparedMaterialExchangeNetwork<'_, '_> {
    /// Resolve fresh values on ONE trial. Start with fresh bulk coefficients;
    /// add thermal contact resistance once and radiation heat once. Electrical
    /// loss deposition remains exclusively with the existing circuit owner.
    pub fn resolve(&self, cx: &Cx<'_>, thermal: &DepositionSimulation<'_>)
        -> Result<ContactThermalCoefficients<MaterialExchangeSample>, MaterialExchangeError>
    {
        checkpoint(cx)?;
        let mut resolved = self.electrical.resolve(thermal, || cx.checkpoint().is_err())?;
        self.binding.admit_common_states(cx, &resolved.bulk.source)?;
        let contacts = apply_material_contacts(thermal, &self.thermal_contacts,
            &mut resolved.bulk.conductances_w_k, || cx.checkpoint().is_err())?;
        let radiation = match self.binding.radiation {
            None => None,
            Some(binding) => {
                let value = MaterialRadiation { enclosure: binding.enclosure, cells: &self.radiating_cells,
                    power_tolerance_w: binding.power_tolerance_w,
                    radiosity_tolerance_w_m2: binding.radiosity_tolerance_w_m2 }.resolve(cx, thermal)?;
                value.add_to(&mut resolved.bulk.sources_w)?;
                Some(value)
            }
        };
        checkpoint(cx)?;
        Ok(ContactThermalCoefficients {
            contacts: resolved.contacts,
            bulk: CircuitThermalCoefficients {
                resistors: resolved.bulk.resistors, sources_w: resolved.bulk.sources_w,
                conductances_w_k: resolved.bulk.conductances_w_k,
                source: MaterialExchangeSample { bulk: resolved.bulk.source, contacts, radiation },
            },
        })
    }

    /// One simultaneous electrical/contact/radiative/thermal/gas-motion interval.
    /// Empty chamber/patch slices give the finite electrical/thermal case. The
    /// original network residual, per-cell balances and signed gas transfers are
    /// reused. Closed-enclosure roundoff is NOT counted as external heat supply.
    ///
    /// Independent walls, fixed geometry/topology/C/L, zero interface storage,
    /// uniform thermal cells and nonparticipating closed gases remain assumptions.
    /// Radiation and gas heat may coexist on the same exposed surface. Each
    /// mechanism's trace-area constraints and source regimes still apply.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn step_candidate<'m>(
        &self, cx: &Cx<'_>, thermal: &DepositionSimulation<'m>, target_s: f64,
        circuit: &MaterialCircuit, ports: &[f64], chambers: &[Chamber<'_>],
        gas_patches: &[MaterialContact], controls: NetworkSolve, guess_j: &[f64],
    ) -> Result<NetworkCandidate<'m, ContactThermalSample<MaterialExchangeSample>, ResolvedMaterialStatePoint>,
        CircuitEnthalpyError<NetworkError<MaterialExchangeError, MaterialNetworkError>>>
    {
        let (surfaces, mapped) = conductive_surfaces(cx, gas_patches, controls)?;
        self.step_with_gas_law_candidate(cx, thermal, target_s, circuit, ports,
            chambers, &surfaces, controls, guess_j, |contact, state, source, gas| {
                self.resolve_conductive_gas_contact(mapped[&contact.id], state, source, gas)
            })
    }
}
// Material geometry admission is shared by both native chamber time policies.
#[allow(clippy::type_complexity)]
fn conductive_surfaces<B>(cx: &Cx<'_>, patches: &[MaterialContact], controls: NetworkSolve)
    -> Result<(Vec<gas_law::GasSurface>, std::collections::BTreeMap<u64, MaterialContact>),
        CircuitEnthalpyError<DrivenNetworkError<MaterialExchangeError, B>>>
{
    if patches.len() > controls.maximum_contacts {
        return Err(CircuitEnthalpyError::Invalid("material gas contact capacity exceeded"));
    }
    let mut surfaces = Vec::with_capacity(patches.len());
    let mut mapped = std::collections::BTreeMap::new();
    for &patch in patches {
        checkpoint(cx).map_err(coupled_error)?;
        let geometry = patch.geometry;
        if ![geometry.gas_distance_m, geometry.cell_distance_m].iter()
            .all(|length| length.is_finite() && *length > 0.0)
            || mapped.insert(patch.contact.id, patch).is_some() {
            return Err(CircuitEnthalpyError::Invalid("unique gas contacts and positive finite layer distances required"));
        }
        surfaces.push(gas_law::GasSurface { contact: patch.contact, area_m2: geometry.area_m2 });
    }
    Ok((surfaces, mapped))
}
fn coupled_error<B>(e: MaterialExchangeError)
    -> CircuitEnthalpyError<DrivenNetworkError<MaterialExchangeError, B>>
{
    match e {
        MaterialExchangeError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(DrivenNetworkError::Material(other)),
    }
}
