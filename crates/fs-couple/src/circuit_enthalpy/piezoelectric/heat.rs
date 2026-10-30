//! Complete source-defined heat transfer on the original piezoelectric driver.
//!
//! The assembly resolves transport; it never advances a second thermal or
//! mechanical history. The rod retains its independently held ceramic material.
//! Temperature-dependent ceramic storage uses the separate thermal entry.

use super::{PiezoelectricEnthalpyCandidate, PiezoelectricEnthalpyInput, resolved_candidate};
use crate::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpySolve, admit_driver_unknowns,
    contacts::heat::{HeatAssembly, HeatNetworkError, closed_balance},
    field_storage::admit_circuit_sources,
};
use crate::material_circuit::piezoelectric::PiezoelectricCircuit;
use fs_exec::Cx;
use fs_thermal::enthalpy::DepositionSimulation;

fn heat_error<E: From<HeatNetworkError>>(error: HeatNetworkError) -> CircuitEnthalpyError<E> {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(other.into()),
    }
}

/// Advance the original reciprocal rod, current and finite material heat with
/// one complete source assembly at each proposed thermal endpoint.
///
/// Material convection, thermal contacts, integrated conduction and radiation
/// use the same proposed temperatures as the original load resistance and
/// Joule heat. A finite coolant stays in the native thermal inventory; a held
/// convection reservoir contributes its actual external heat. The full accepted
/// assembly sample is retained in the existing report.
///
/// Original conductor sources, branches and receiving cells must agree before
/// any trial. Electrical contacts need a piezoelectric contact-work owner and
/// are refused. Every trial starts from the same accepted charge, flux, strain
/// and momentum, with each physical resistor loss deposited exactly once.
///
/// Closed radiation remains internal, including its retained numerical net.
/// The same energy tolerance gates that correction; true external heat is
/// `non_electrical_heat_j - dt * source.as_ref().radiation_net_power_w()`.
/// No second solve or correction changes a cell's energy.
///
/// Ceramic coefficients, reference inertia and source context stay fixed under
/// the base small-signal model. This entry does not infer ceramic temperature
/// feedback, thermal strain, dielectric loss or depoling. Heat geometry and
/// contact/view-factor topology remain fixed; the original first-order thermal
/// composition requires time refinement. Any source, budget, native solver or
/// cancellation refusal returns without changing accepted histories.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'t, A: HeatAssembly + ?Sized>(
    assembly: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    model: &PiezoelectricCircuit,
    electrical_ports: &[f64],
    force_n: f64,
    controls: CircuitEnthalpySolve,
    guess_j: &[f64],
) -> Result<PiezoelectricEnthalpyCandidate<'t, A::Sample>, CircuitEnthalpyError<A::Error>> {
    cx.checkpoint()
        .map_err(|_| CircuitEnthalpyError::Cancelled)?;
    admit_driver_unknowns(thermal, controls, guess_j)?;
    let prepared = assembly.prepared_heat_network();
    // This accessor refuses electrical contacts before their work can be lost.
    prepared
        .field_material_network()
        .map_err(heat_error::<A::Error>)?;
    prepared
        .admit_field_independent_boundaries(cx, &[])
        .map_err(heat_error::<A::Error>)?;
    let destinations = prepared.resistor_destinations();
    let mut initial = assembly.resolve_heat(cx, thermal)?;
    if !initial.contacts.is_empty() {
        return Err(CircuitEnthalpyError::Invalid(
            "piezoelectric circuit cannot discard electrical contact work or heat",
        ));
    }
    admit_circuit_sources(
        model.circuit(),
        &mut initial.bulk.resistors,
        &mut initial.contacts,
        || cx.checkpoint().is_err(),
    )?;
    drop(initial);
    let candidate = resolved_candidate(
        thermal,
        target_s,
        PiezoelectricEnthalpyInput {
            model,
            electrical_ports,
            force_n,
            destinations: &destinations,
        },
        controls,
        guess_j,
        |state| {
            let coefficients = assembly.resolve_heat(cx, state)?;
            if !coefficients.contacts.is_empty() {
                return Err(CircuitEnthalpyError::Invalid(
                    "piezoelectric circuit cannot discard electrical contact work or heat",
                ));
            }
            Ok(coefficients.bulk)
        },
        || cx.checkpoint().is_err(),
    )?;
    closed_balance(
        target_s - thermal.time_s(),
        candidate.report().total_energy_residual_j,
        candidate.report().source.as_ref(),
        controls.total_energy_tolerance_j,
    )
    .map_err(heat_error::<A::Error>)?;
    cx.checkpoint()
        .map_err(|_| CircuitEnthalpyError::Cancelled)?;
    Ok(candidate)
}
