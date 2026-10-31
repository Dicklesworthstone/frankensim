//! Contact Joule heat enters BOTH finite material bodies on the shared solve.
//! No guessed 50/50 partition, separate thermal step, or lost contact energy.
pub mod material;
/// Compose material radiation and thermal-contact resistance on the same cells.
pub mod heat;
use super::*;
use crate::material_circuit::contacts::{SeriesElectricalContact,ContactHeat,MAX_SERIES_CONTACTS};
use fs_material::electrical_interface::ELECTRICAL_CONTACT_RESISTANCE_PROPERTY;

/// Fixed physical heat partition for one zero-storage electrical interface.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactHeatCells {
    /// Original electrical contact ID.
    pub id: u64,
    /// Active thermal cells in interface source A/B order.
    pub cells: [CellId;2],
    /// Fraction of contact heat entering A, explicitly supplied in [0,1].
    /// The remainder enters B. This is a model input, not inferred effusivity.
    pub heat_fraction_to_first: f64,
}
/// One coherent bulk/interface material evaluation of the same thermal trial.
pub struct ContactThermalCoefficients<S> {
    /// Original bulk electrical/thermal coefficients and additional heat sources.
    pub bulk: CircuitThermalCoefficients<S>,
    /// Fresh contacts, with the accepted IDs, branch assignments and source pins.
    pub contacts: Vec<SeriesElectricalContact>,
}
/// Retain physical heat destinations with the original material payload.
pub struct ContactThermalSample<S> {
    /// Original bulk and optional extra thermal source records.
    pub material: S,
    /// Exact partition inputs used by every residual and final admission.
    pub destinations: Vec<ContactHeatCells>,
}

pub(crate) fn admit_routes<E>(thermal:&DepositionSimulation<'_>,circuit:&MaterialCircuit,
    routes:&[ContactHeatCells],mut checkpoint:impl FnMut()->Result<(),DepositionError>)
    ->Result<BTreeMap<u64,ContactHeatCells>,CircuitEnthalpyError<E>> {
    if routes.len()>MAX_SERIES_CONTACTS || routes.len()!=circuit.contacts().len() {
        return Err(CircuitEnthalpyError::Invalid("every electrical contact needs one complete heat partition"));
    }
    let mut result=BTreeMap::new();
    for &route in routes {
        checkpoint()?;
        if !route.heat_fraction_to_first.is_finite() || !(0.0..=1.0).contains(&route.heat_fraction_to_first)
            || route.cells.iter().any(|c|!thermal.cells().contains_key(c))
            || !circuit.contacts().iter().any(|c|c.id()==route.id)
            || result.insert(route.id,route).is_some() {
            return Err(CircuitEnthalpyError::Invalid("invalid, duplicate or inactive contact heat destination"));
        }
    }
    Ok(result)
}

pub(crate) fn deposit<E>(losses:&[ContactHeat],routes:&BTreeMap<u64,ContactHeatCells>,dt:f64,
    powers:&mut BTreeMap<CellId,f64>,mut checkpoint:impl FnMut()->Result<(),DepositionError>)
    ->Result<(),CircuitEnthalpyError<E>> {
    if losses.len()!=routes.len() {return Err(CircuitEnthalpyError::Invalid("contact heat coverage changed during trial"));}
    for loss in losses {
        checkpoint()?;
        let route=routes.get(&loss.contact.id()).ok_or(CircuitEnthalpyError::Invalid("unrouted electrical contact heat"))?;
        let fraction=route.heat_fraction_to_first;
        let first=loss.heat_j*fraction;
        let second=loss.heat_j-first;
        if !first.is_finite() || !second.is_finite() || first<0.0 || second<0.0
            || (loss.heat_j!=0.0 && ((fraction>0.0 && first==0.0) || (fraction<1.0 && second==0.0))) {
            return Err(CircuitEnthalpyError::Invalid("electrical contact heat partition is unrepresentable"));
        }
        for (cell,heat) in route.cells.into_iter().zip([first,second]) {
            let power=finite(heat/dt)?;
            let previous=powers.get(&cell).copied().unwrap_or(0.0);
            let next=finite(previous+power)?;
            if heat!=0.0 && (power==0.0 || previous==next) {
                return Err(CircuitEnthalpyError::Invalid("contact heat disappears in receiving cell power"));
            }
            powers.insert(cell,next);
        }
    }
    Ok(())
}

/// Simultaneously solve bulk/contact currents, both receiving enthalpies and
/// all existing thermal links/boundaries. Unknowns and acceptance remain with
/// the original circuit_enthalpy driver. Each trial restarts from the same
/// charges/fluxes, and every joule of bulk AND interface loss is deposited once.
///
/// Source cards and selected contact claims must match the accepted circuit;
/// reconfiguration is explicit, not fallback during a nonlinear iteration.
/// Contact-free callers retain the existing entry point and payload shape.
/// Geometry/topology and partition fractions are fixed over this interval.
/// Zero interface storage and uniform cell temperatures are model assumptions;
/// no opening, arc, contact-mechanics or spreading-temperature law is inferred.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'m,S,E>(thermal:&DepositionSimulation<'m>,target_s:f64,
    input:CircuitEnthalpyInput<'_>,routes:&[ContactHeatCells],controls:CircuitEnthalpySolve,guess_j:&[f64],
    resolve:impl Fn(&DepositionSimulation<'m>)->Result<ContactThermalCoefficients<S>,E>,
    cancelled:impl FnMut()->bool)
    ->Result<CircuitEnthalpyCandidate<'m,ContactThermalSample<S>>,CircuitEnthalpyError<E>> {
    resolved_candidate(thermal,target_s,input,routes,controls,guess_j,
        |state|resolve(state).map_err(CircuitEnthalpyError::Source),cancelled)
}

// Complete assemblies already distinguish source, thermal and cancellation
// failures. Keep those native classes while sharing the exact contact solver.
#[allow(clippy::too_many_arguments)]
pub(crate) fn resolved_candidate<'m,S,E>(thermal:&DepositionSimulation<'m>,target_s:f64,
    input:CircuitEnthalpyInput<'_>,routes:&[ContactHeatCells],controls:CircuitEnthalpySolve,guess_j:&[f64],
    resolve:impl Fn(&DepositionSimulation<'m>)->Result<ContactThermalCoefficients<S>,CircuitEnthalpyError<E>>,
    cancelled:impl FnMut()->bool)
    ->Result<CircuitEnthalpyCandidate<'m,ContactThermalSample<S>>,CircuitEnthalpyError<E>> {
    let cancel=RefCell::new(cancelled);
    let check=||if (*cancel.borrow_mut())() {Err(DepositionError::Cancelled)} else {Ok(())};
    check()?;
    let routes=admit_routes(thermal,input.circuit,routes,check)?.into_values().collect::<Vec<_>>();
    driven_with_contacts_candidate(thermal,target_s,input,&routes,controls,guess_j,|state,_| {
        let coefficients=resolve(state)?;
        let bulk=coefficients.bulk;
        let circuit=advance_circuit::<E>(input.circuit,bulk.resistors,coefficients.contacts,input.ports,
            target_s-thermal.time_s(),controls.circuit,||(*cancel.borrow_mut())())?;
        Ok(DrivenCoefficients {circuit,sources_w:bulk.sources_w,conductances_w_k:bulk.conductances_w_k,
            source:ContactThermalSample {material:bulk.source,destinations:routes.clone()},circuit_trials:1})
    },||(*cancel.borrow_mut())())
}

// Shared by ordinary contact heating and gas/mechanical network composition.
// This only proposes electrical history: the thermal owner alone deposits loss.
// Retain one source-continuity check instead of copying it into every coupling.
#[allow(clippy::too_many_arguments)]
pub(crate) fn advance_circuit<E>(
    circuit:&MaterialCircuit, resistors:Vec<MaterialResistor>,
    contacts:Vec<SeriesElectricalContact>, ports:&[f64], dt:f64,
    controls:MaterialCircuitSolve, mut cancelled:impl FnMut()->bool,
) -> Result<MaterialCircuitCandidate,CircuitEnthalpyError<E>> {
    admit_sources(circuit,&contacts,&mut cancelled)?;
    // Keep the contact-free circuit's exact existing arithmetic and entry point.
    if contacts.is_empty() {
        Ok(circuit.step_candidate(resistors,ports,dt,controls,cancelled)?)
    } else {
        Ok(circuit.step_with_contacts_candidate(resistors,contacts,ports,dt,controls,cancelled)?)
    }
}

// The original contact-source continuity law also applies when the circuit
// owns a moving field gap. Neither driver may change a source during feedback.
pub(crate) fn admit_sources<E>(circuit:&MaterialCircuit,contacts:&[SeriesElectricalContact],
    mut cancelled:impl FnMut()->bool) -> Result<(),CircuitEnthalpyError<E>> {
    if cancelled() {return Err(CircuitEnthalpyError::Cancelled);}
    if contacts.len()!=circuit.contacts().len() {
        return Err(CircuitEnthalpyError::Invalid("contact source coverage changed during feedback"));
    }
    for contact in contacts {
        if cancelled() {return Err(CircuitEnthalpyError::Cancelled);}
        let old=circuit.contacts().iter().find(|c|c.id()==contact.id())
            .ok_or(CircuitEnthalpyError::Invalid("unbound electrical contact source"))?;
        let before=old.source().resolved();let after=contact.source().resolved();
        let claim=|s:&fs_material::state_point::ResolvedInterfaceStatePoint|s.property(ELECTRICAL_CONTACT_RESISTANCE_PROPERTY)
            .map(|p|p.answer().receipt.selected);
        if before.card_identity()!=after.card_identity() || claim(before)!=claim(after) {
            return Err(CircuitEnthalpyError::Invalid("electrical contact source changed during thermal feedback"));
        }
    }
    Ok(())
}
