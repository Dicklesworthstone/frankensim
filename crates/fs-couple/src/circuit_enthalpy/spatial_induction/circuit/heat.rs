//! Complete original heat transport around the retained electrode/RLC history.

use super::super::terminals::circuit::{close_heat_balance, heat_error};
use super::*;
use crate::circuit_enthalpy::contacts::heat::HeatAssembly;
use fs_exec::Cx;

/// Advance the retained three-dimensional electrode circuit and every finite
/// caloric receiver through the complete original material heat assembly.
///
/// The same proposed endpoint temperatures determine conductor conductivity,
/// magnetic entropy and internal energy, sourced thermal contacts, integrated
/// layers, radiation and finite coolant exchange. Explicit original graph
/// resistor bindings also refresh their own R(T); otherwise authored graph
/// resistance stays fixed. Every physical field or graph loss reaches its
/// original receiver exactly once. The native interface reference resistance
/// remains numerical and cannot supply heat.
///
/// The assembly retains its complete accepted `A::Sample`, including provider
/// payloads and finite coolant sources. Provider failures keep their original
/// typed error. Closed-enclosure radiation and finite thermal exchange remain
/// internal to the same enthalpy solve; no separate cooling step is advanced.
///
/// Graph material resistors occupy distinct, whole-volume finite bodies outside
/// the magnetic mesh. Thermal boundaries are field independent. Electrical
/// interface contacts require their own native electrical work owner and refuse.
/// Geometry, reactive coefficients and source selections stay fixed; the
/// original first-order thermal coupling and linear-magnetic source admission
/// apply. Prescribed motion exchanges signed mechanical work without advancing
/// mechanics. Failure publishes none of the field, circuit or thermal histories.
#[allow(clippy::too_many_arguments)]
pub fn heat_assembly_candidate<'t, 'm, A: HeatAssembly + ?Sized, C: SpatialConductorMaterial>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    circuit: &SpatialTerminalCircuit<'m, C>,
    external_inputs: &[f64],
    endpoint_impressed_current_a_m2: &[[f64; 3]],
    velocity_m_s: &[[f64; 3]],
    regions: &[SpatialMagneticRegion],
    graph_resistors: &[ResistorCell],
    volume_tolerance_m3: f64,
    magnetic_laws: &[SpatialMagneticThermalLaw],
    controls: CircuitEnthalpySolve,
    circuit_controls: SpatialTerminalCircuitControls,
    assembly: WeightedAssemblyLimits,
    initial_energy_increments_j: &[f64],
) -> Result<SpatialCircuitThermalCandidate<'t, 'm, C, A::Sample>, TerminalError<A::Error>> {
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    check(&cancel).map_err(native_error)?;
    crate::circuit_enthalpy::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
    let prepared = provider.prepared_heat_network();
    let bulk = prepared.field_material_network().map_err(heat_error)?;
    prepared
        .admit_field_independent_boundaries(cx, &[])
        .map_err(heat_error)?;
    let layout = SpatialTerminalCircuitThermalNetwork {
        magnetic: SpatialMagneticThermalNetwork {
            cells: bulk.cells,
            links: bulk.links,
            magnetic_regions: regions,
            volume_tolerance_m3,
        },
        graph_resistors,
    };
    let mut next = layout.step_linear_resolved_candidate(
        thermal,
        target_s,
        circuit,
        external_inputs,
        endpoint_impressed_current_a_m2,
        velocity_m_s,
        (!bulk.resistors.is_empty()).then_some(bulk.resistors),
        magnetic_laws,
        controls,
        circuit_controls,
        assembly,
        initial_energy_increments_j,
        |state| {
            let coefficients = provider
                .resolve_heat(cx, state)
                .map_err(|error| error.map_source(SpatialTerminalCircuitThermalError::Heat))?;
            if !coefficients.contacts.is_empty() {
                return Err(native_error(invalid(
                    "electrode circuit heat cannot discard electrical contact work",
                )));
            }
            Ok(coefficients.bulk)
        },
        |sample: &A::Sample| &sample.as_ref().material,
        &cancel,
    )?;
    close_heat_balance(
        &mut next.report,
        target_s - thermal.time_s(),
        controls.total_energy_tolerance_j,
    )?;
    check(&cancel).map_err(native_error)?;
    Ok(next)
}
