//! Data-only several-chamber coupling with shared source and surface admission.
use super::*;
use crate::circuit_enthalpy::material::{MaterialNetwork,MaterialNetworkError,MaterialNetworkSample};
use crate::cavity_phs::piston::thermal::boundary::network::material::{GasConductionGeometry,resolve_contact};
use fs_material::state_point::ResolvedMaterialStatePoint;

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
#[allow(clippy::too_many_arguments,clippy::type_complexity)]
pub fn step_candidate<'m>(
    binding:&MaterialNetwork<'_>,thermal:&DepositionSimulation<'m>,target_s:f64,
    circuit:&MaterialCircuit,ports:&[f64],chambers:&[Chamber<'_>],contacts:&[MaterialContact],
    controls:NetworkSolve,guess_j:&[f64],cancelled:impl FnMut()->bool,
) -> Result<NetworkCandidate<'m,MaterialNetworkSample,ResolvedMaterialStatePoint>,
    CircuitEnthalpyError<NetworkError<MaterialNetworkError,MaterialNetworkError>>> {
    let cancel=RefCell::new(cancelled);
    if (*cancel.borrow_mut())() { return Err(CircuitEnthalpyError::Cancelled); }
    if contacts.len()>controls.maximum_contacts || chambers.len()>controls.maximum_chambers
        || binding.resistors.len()>crate::material_circuit::MAX_MATERIAL_CIRCUIT_BRANCHES {
        return Err(CircuitEnthalpyError::Invalid("material chamber/contact/circuit capacity exceeded"));
    }
    let mut mapped=BTreeMap::new(); let mut areas=BTreeMap::<CellId,f64>::new();
    for contact in contacts {
        if (*cancel.borrow_mut())() { return Err(CircuitEnthalpyError::Cancelled); }
        let geometry=contact.geometry;
        if !geometry.area_m2.is_finite() || geometry.area_m2<0.0
            || !positive(geometry.gas_distance_m) || !positive(geometry.cell_distance_m)
            || mapped.insert(contact.contact.id,contact).is_some() {
            return Err(CircuitEnthalpyError::Invalid("unique material contact IDs and valid layer geometry required"));
        }
        let area=areas.entry(contact.contact.cell).or_insert(0.0); *area+=geometry.area_m2;
        let cell=thermal.cells().get(&contact.contact.cell)
            .ok_or(CircuitEnthalpyError::Invalid("inactive gas contact cell"))?;
        if !area.is_finite() || *area>thermal.exposed_area_m2(cell.id)? || cell.boundary.htc_w_m2_k!=0.0 {
            return Err(CircuitEnthalpyError::Invalid("combined gas patches exceed exposed surface or duplicate ambient convection"));
        }
    }
    let destinations:Vec<_>=binding.resistors.iter().map(|r|r.destination).collect();
    let patches:Vec<_>=mapped.values().map(|p|p.contact).collect();
    let result=super::step_candidate(thermal,target_s,CircuitEnthalpyInput {circuit,ports,destinations:&destinations},
        chambers,&patches,controls,guess_j,|state|binding.resolve(state,||(*cancel.borrow_mut())()),
        |contact,state,source,gas|resolve_contact(binding,contact.cell,state,source,gas,mapped[&contact.id].geometry),
        ||(*cancel.borrow_mut())());
    match result {
        Err(CircuitEnthalpyError::Source(NetworkError::Material(MaterialNetworkError::Thermal(DepositionError::Cancelled)))) =>
            Err(CircuitEnthalpyError::Cancelled),
        other=>other,
    }
}
