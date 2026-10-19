//! Data-only several-chamber coupling with shared source and surface admission.
/// Compose electrical heating, thermal contacts and radiation on the shared gas network.
pub mod exchange;

use super::*;
use crate::circuit_enthalpy::material::{MaterialNetwork,MaterialNetworkError,MaterialNetworkSample};
use crate::circuit_enthalpy::contacts::material::{ElectricalContactNetwork,ElectricalContactNetworkError};
use crate::cavity_phs::piston::thermal::boundary::network::material::{GasConductionGeometry,resolve_contact};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::conductor::ElectricalPropertyBasis;

/// One physical gas/cell patch using the original material-layer geometry.
#[derive(Clone,Copy,Debug,PartialEq)]
pub struct MaterialContact {
    /// Stable contact, chamber and existing thermal-cell identities.
    pub contact: Contact,
    /// Effective layer lengths and physical contact area, not electrical area.
    pub geometry: GasConductionGeometry,
}

/// Simultaneously evolve the finite circuit, material cells and ALL chambers.
/// No material callback, gas-name preset or second heat store is needed.
///
/// Every contact reuses the existing single-boundary material/conduction owner
/// and exact network conductivity snapshot. All patch areas on a cell are summed
/// across ALL chambers before comparing against its exposed surface. Gas patches
/// cannot overlap authored links/reservoirs or be exposed to ambient convection
/// again. Geometric non-overlap within that area budget remains a caller claim;
/// area arithmetic alone is not a surface-intersection test.
///
/// Each chamber keeps its original gas composition and independent translating
/// wall. This does not connect gas masses or impose shared mechanical motion.
/// The legacy electrical basis remains explicitly resistivity-only.
#[allow(clippy::too_many_arguments,clippy::type_complexity)]
pub fn step_candidate<'m>(
    binding:&MaterialNetwork<'_>,thermal:&DepositionSimulation<'m>,target_s:f64,
    circuit:&MaterialCircuit,ports:&[f64],chambers:&[Chamber<'_>],contacts:&[MaterialContact],
    controls:NetworkSolve,guess_j:&[f64],mut cancelled:impl FnMut()->bool,
) -> Result<NetworkCandidate<'m,MaterialNetworkSample,ResolvedMaterialStatePoint>,
    CircuitEnthalpyError<NetworkError<MaterialNetworkError,MaterialNetworkError>>> {
    if cancelled() {return Err(CircuitEnthalpyError::Cancelled);}
    if binding.resistors.len()>crate::material_circuit::MAX_MATERIAL_CIRCUIT_BRANCHES {
        return Err(CircuitEnthalpyError::Invalid("material circuit capacity exceeded"));
    }
    let bases=binding.resistors.iter().map(|r|(r.destination.branch,ElectricalPropertyBasis::Resistivity)).collect::<Vec<_>>();
    step_candidate_with_bases(binding,thermal,target_s,circuit,ports,&bases,chambers,contacts,controls,guess_j,cancelled)
}

/// Explicit conductivity/resistivity choices for every original branch.
/// Prepare transport claims once, then resolve live values at every thermal trial.
/// All chambers share that same selection and exact complete material snapshots.
#[allow(clippy::too_many_arguments,clippy::type_complexity)]
pub fn step_candidate_with_bases<'m>(
    binding:&MaterialNetwork<'_>,thermal:&DepositionSimulation<'m>,target_s:f64,
    circuit:&MaterialCircuit,ports:&[f64],bases:&[(usize,ElectricalPropertyBasis)],
    chambers:&[Chamber<'_>],contacts:&[MaterialContact],controls:NetworkSolve,guess_j:&[f64],cancelled:impl FnMut()->bool,
) -> Result<NetworkCandidate<'m,MaterialNetworkSample,ResolvedMaterialStatePoint>,
    CircuitEnthalpyError<NetworkError<MaterialNetworkError,MaterialNetworkError>>> {
    let cancel=RefCell::new(cancelled);
    let mapped=admit_surfaces(thermal,chambers,contacts,controls,||(*cancel.borrow_mut())())?;
    let prepared=match binding.prepare_with_bases(thermal,bases,||(*cancel.borrow_mut())()) {
        Err(MaterialNetworkError::Thermal(DepositionError::Cancelled))=>return Err(CircuitEnthalpyError::Cancelled),
        other=>other.map_err(|e|CircuitEnthalpyError::Source(NetworkError::Material(e)))?,
    };
    let destinations:Vec<_>=binding.resistors.iter().map(|r|r.destination).collect();
    let patches:Vec<_>=mapped.values().map(|p|p.contact).collect();
    let result=super::step_candidate(thermal,target_s,CircuitEnthalpyInput {circuit,ports,destinations:&destinations},
        chambers,&patches,controls,guess_j,|state|prepared.resolve(state,||(*cancel.borrow_mut())()),
        |contact,state,source,gas|resolve_contact(binding,contact.cell,state,source,gas,mapped[&contact.id].geometry),
        ||(*cancel.borrow_mut())());
    match result {
        Err(CircuitEnthalpyError::Source(NetworkError::Material(MaterialNetworkError::Thermal(DepositionError::Cancelled)))) =>
            Err(CircuitEnthalpyError::Cancelled),
        other=>other,
    }
}

/// Data-only electrical contacts, finite material bodies and all gas chambers.
/// The electrical binding owns heat partitions and thermal-link contact area;
/// the gas patches own their separate exposed areas and layer distances.
/// Original bulk and interface claims are pinned, not their evolving values.
/// Source snapshots go directly from this network into the shared gas contact.
///
/// All electrical losses and gas transfers are already included in the returned
/// endpoints. The original fixed-geometry, uniform-state, independent-wall limits
/// apply. No conductivity, material, heat split or mechanical linkage is guessed.
#[allow(clippy::too_many_arguments,clippy::type_complexity)]
pub fn step_with_electrical_contacts_candidate<'m>(
    binding:&ElectricalContactNetwork<'_>,thermal:&DepositionSimulation<'m>,target_s:f64,
    circuit:&MaterialCircuit,ports:&[f64],chambers:&[Chamber<'_>],contacts:&[MaterialContact],
    controls:NetworkSolve,guess_j:&[f64],cancelled:impl FnMut()->bool,
) -> Result<NetworkCandidate<'m,ContactThermalSample<MaterialNetworkSample>,ResolvedMaterialStatePoint>,
    CircuitEnthalpyError<NetworkError<ElectricalContactNetworkError,MaterialNetworkError>>> {
    let cancel=RefCell::new(cancelled);
    let mapped=admit_surfaces(thermal,chambers,contacts,controls,||(*cancel.borrow_mut())())?;
    let prepared=binding.prepare(thermal,||(*cancel.borrow_mut())()).map_err(contact_error)?;
    let routes=binding.heat_destinations(thermal).map_err(contact_error)?;
    let destinations:Vec<_>=binding.bulk.resistors.iter().map(|r|r.destination).collect();
    let patches:Vec<_>=mapped.values().map(|p|p.contact).collect();
    let result=super::step_with_contacts_candidate(thermal,target_s,
        CircuitEnthalpyInput {circuit,ports,destinations:&destinations},&routes,
        chambers,&patches,controls,guess_j,|state|prepared.resolve(state,||(*cancel.borrow_mut())()),
        |contact,state,source,gas|resolve_contact(binding.bulk,contact.cell,state,source,gas,mapped[&contact.id].geometry),
        ||(*cancel.borrow_mut())());
    match result {
        Err(CircuitEnthalpyError::Source(NetworkError::Material(error)))=>Err(contact_error(error)),
        other=>other,
    }
}

fn contact_error(error:ElectricalContactNetworkError)
    ->CircuitEnthalpyError<NetworkError<ElectricalContactNetworkError,MaterialNetworkError>> {
    match error {
        ElectricalContactNetworkError::Thermal(DepositionError::Cancelled)
        |ElectricalContactNetworkError::Bulk(MaterialNetworkError::Thermal(DepositionError::Cancelled))=>CircuitEnthalpyError::Cancelled,
        other=>CircuitEnthalpyError::Source(NetworkError::Material(other)),
    }
}

// One aggregate surface check for ordinary and contact-heated circuit networks.
fn admit_surfaces<'a,E>(thermal:&DepositionSimulation<'_>,chambers:&[Chamber<'_>],
    contacts:&'a [MaterialContact],controls:NetworkSolve,mut cancelled:impl FnMut()->bool)
    ->Result<BTreeMap<u64,&'a MaterialContact>,CircuitEnthalpyError<E>> {
    if cancelled() {return Err(CircuitEnthalpyError::Cancelled);}
    if contacts.len()>controls.maximum_contacts || chambers.len()>controls.maximum_chambers {
        return Err(CircuitEnthalpyError::Invalid("material chamber/contact capacity exceeded"));
    }
    let mut mapped=BTreeMap::new();let mut areas=BTreeMap::<CellId,f64>::new();
    for contact in contacts {
        if cancelled() {return Err(CircuitEnthalpyError::Cancelled);}
        let geometry=contact.geometry;
        if !geometry.area_m2.is_finite() || geometry.area_m2<0.0
            || !positive(geometry.gas_distance_m) || !positive(geometry.cell_distance_m)
            || mapped.insert(contact.contact.id,contact).is_some() {
            return Err(CircuitEnthalpyError::Invalid("unique material contact IDs and valid layer geometry required"));
        }
        let area=areas.entry(contact.contact.cell).or_insert(0.0);*area+=geometry.area_m2;
        let cell=thermal.cells().get(&contact.contact.cell)
            .ok_or(CircuitEnthalpyError::Invalid("inactive gas contact cell"))?;
        if !area.is_finite() || *area>thermal.exposed_area_m2(cell.id)? || cell.boundary.htc_w_m2_k!=0.0 {
            return Err(CircuitEnthalpyError::Invalid("combined gas patches exceed exposed surface or duplicate ambient convection"));
        }
    }
    Ok(mapped)
}
