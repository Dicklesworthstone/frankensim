//! Material convection on the original finite circuit/cell/multi-chamber solve.
//!
//! Every patch receives its chamber's actual gas-energy midpoint and its cell's
//! actual caloric endpoint. One original nonlinear assembly debits cell heat,
//! credits gas heat and evolves charge/flux and wall momentum. No ambient gas
//! cache, extra Joule deposition, second thermal state or integration law.

use super::{FiniteGasConvectionSample, ForcedGasSurface, MaterialConvectionError, invalid};
use fs_conduction::transient::deposition::{CellId, DepositionError, DepositionSimulation};
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::{
    self as coupled, Chamber, Contact, NetworkCandidate, NetworkError, NetworkSolve,
};
use fs_couple::circuit_enthalpy::{CircuitEnthalpyError, CircuitEnthalpyInput};
use fs_couple::circuit_enthalpy::material::{MaterialNetwork, MaterialNetworkError, MaterialNetworkSample};
use fs_couple::material_circuit::MaterialCircuit;
use fs_material::conductor::ElectricalPropertyBasis;
use std::cell::RefCell;
use std::collections::BTreeMap;

/// Distinct exposed patch on a real thermal cell and finite gas chamber.
/// Several patches may share either owner; geometry, not IDs, defines their area.
#[derive(Clone, Copy, Debug)]
pub struct ConvectiveContact {
    pub contact: Contact,
    pub surface: ForcedGasSurface,
}

/// Original material, thermal, electrical or gas/correlation refusal.
pub type ConvectionNetworkError = CircuitEnthalpyError<
    NetworkError<MaterialNetworkError, MaterialConvectionError>,
>;

/// Propose ALL finite material cells, circuit storage and chamber momenta together.
///
/// Uses the original material network's complete caloric cards and pinned
/// conductivity/resistivity selections. `electrical_bases` explicitly covers
/// every resistor; it is empty only when the binding has no resistors. No wall
/// conductivity is demanded for a lumped convective patch unless the original
/// thermal/electrical network itself needs it for another claimed response.
///
/// The sum of patches from ALL chambers must fit each cell's remaining exposed
/// area after original contacts/reservoirs. Nonzero authored ambient convection
/// on a bound cell refuses to prevent double counting. Area sums cannot prove
/// geometric non-overlap: the caller must supply disjoint physical patches.
///
/// Each chamber is advanced ONCE per original outer trial with its TOTAL signed
/// boundary heat, while each cell receives its OWN signed contact debit. Thus
/// simultaneous heating and cooling contacts do not disappear when net chamber
/// heat happens to cancel. Electrical dissipation is deposited by the original
/// circuit owner, never again by this adapter.
///
/// This is the existing first-order thermal-endpoint/gas-midpoint coupling with
/// uniform cells/gases, fixed thermal geometry, independent rigid walls and
/// prescribed relative speeds. It does not add gas mass transport, chemistry,
/// spatial boundary layers, fan power or viscous heating. Use physical event
/// boundaries and temporal refinement. Return contains all endpoints or none.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn step_candidate<'m>(
    binding: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'m>, target_s: f64,
    circuit: &MaterialCircuit, ports: &[f64], electrical_bases: &[(usize, ElectricalPropertyBasis)],
    chambers: &[Chamber<'_>], contacts: &[ConvectiveContact],
    controls: NetworkSolve, guess_j: &[f64], cancelled: impl FnMut() -> bool,
) -> Result<NetworkCandidate<'m, MaterialNetworkSample, FiniteGasConvectionSample>, ConvectionNetworkError> {
    let cancel = RefCell::new(cancelled);
    let check = || -> Result<(), ConvectionNetworkError> {
        if (*cancel.borrow_mut())() { Err(CircuitEnthalpyError::Cancelled) } else { Ok(()) }
    };
    check()?;
    if contacts.len() > controls.maximum_contacts || chambers.len() > controls.maximum_chambers {
        return Err(CircuitEnthalpyError::Invalid("convective chamber/contact capacity exceeded"));
    }
    // Canonical IDs order reports. Validate their full map before any solve.
    let mut mapped = BTreeMap::new();
    for patch in contacts {
        check()?;
        if patch.surface.validate().is_err()
            || !chambers.iter().any(|chamber| chamber.id == patch.contact.chamber)
            || mapped.insert(patch.contact.id, patch).is_some()
        {
            return Err(CircuitEnthalpyError::Invalid("unique existing chamber/patch IDs and valid convection geometry required"));
        }
    }
    let mut areas = BTreeMap::<CellId, f64>::new();
    for patch in mapped.values() {
        check()?;
        let cell = thermal.cells().get(&patch.contact.cell)
            .ok_or(CircuitEnthalpyError::Invalid("inactive convective cell"))?;
        let area = areas.entry(cell.id).or_insert(0.0);
        *area += patch.surface.area_m2;
        if !area.is_finite() || *area > thermal.exposed_area_m2(cell.id)?
            || cell.boundary.htc_w_m2_k != 0.0
        {
            return Err(CircuitEnthalpyError::Invalid("combined convective patches exceed exposed surface or duplicate ambient convection"));
        }
    }
    let prepared = match binding.prepare_with_bases(thermal, electrical_bases, || (*cancel.borrow_mut())()) {
        Err(MaterialNetworkError::Thermal(DepositionError::Cancelled)) => return Err(CircuitEnthalpyError::Cancelled),
        other => other.map_err(|error| CircuitEnthalpyError::Source(NetworkError::Material(error)))?,
    };
    let destinations = binding.resistors.iter().map(|r| r.destination).collect::<Vec<_>>();
    let patches = mapped.values().map(|p| p.contact).collect::<Vec<_>>();
    let result = coupled::step_candidate(
        thermal, target_s, CircuitEnthalpyInput { circuit, ports, destinations: &destinations },
        chambers, &patches, controls, guess_j,
        |state| prepared.resolve(state, || (*cancel.borrow_mut())()),
        |contact, state, source, gas| {
            let patch = mapped.get(&contact.id).ok_or_else(|| invalid("convective patch disappeared"))?;
            let cell = state.cells().get(&contact.cell).ok_or_else(|| invalid("convective trial omitted its cell"))?;
            let exposed = state.exposed_area_m2(contact.cell)
                .map_err(|error| invalid(&format!("convective trial area: {error}")))?;
            if cell.boundary.htc_w_m2_k != 0.0 || areas[&contact.cell] > exposed {
                return Err(invalid("trial convection overlaps another surface boundary"));
            }
            let sample = source.cells.iter().find(|s| s.cell == contact.cell)
                .ok_or_else(|| invalid("convective trial omitted its original material phase"))?;
            // Preserve the exact complete material solve's endpoint. Do not
            // re-query its chart, replace its temperature, or demand a spare k.
            patch.surface.resolve(gas, sample.phase)
        },
        || (*cancel.borrow_mut())(),
    );
    match result {
        Err(CircuitEnthalpyError::Source(NetworkError::Material(MaterialNetworkError::Thermal(DepositionError::Cancelled)))) =>
            Err(CircuitEnthalpyError::Cancelled),
        other => other,
    }
}

/// Mixed gas heat laws on the complete reacting/contact/radiative assembly.
pub mod exchange;
