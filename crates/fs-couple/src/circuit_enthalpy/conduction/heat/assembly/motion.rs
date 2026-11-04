//! Prescribed motion and resistive joints through the existing complete heat owner.
use super::*;
use crate::circuit_enthalpy::conduction::spatial::motion::SpatialMovingConductionEnthalpyCandidate;
use fs_em::spatial::conduction::{SpatialConductionBoundary, motion::SpatialConductionMotion};

/// Close moving bulk current, stationary resistive interfaces and a complete
/// sourced heat assembly at one finite-enthalpy endpoint.
///
/// Every trial refreshes bulk conductivity and both sides of each electrical
/// contact from the same receiving temperature/phase states used by the
/// original thermal interfaces, conductive layers, reservoirs and radiation.
/// Native conductor-frame loss and explicitly partitioned contact loss enter
/// those finite cells once. Laboratory E.J and Lorentz output are separate:
/// terminal work minus mechanical output supplies the deposited Joule heat.
///
/// Contact-adjacent tetrahedra must have zero prescribed velocity. Moving or
/// sliding interfaces require a surface-work law and remain a native refusal.
/// The held geometry, velocities and magnetic field must preserve the original
/// material occupancy throughout this quasistatic step; no advective inventory,
/// induction, momentum or magnet-field storage is inferred.
///
/// Any complete [`HeatAssembly`] can supply the other thermal transfers. Its
/// source sample and typed errors are retained. Reciprocal links and enclosure
/// radiation remain internal; external non-electrical heat is
/// `report.non_electrical_heat_j - dt * report.source.as_ref().radiation_net_power_w()`.
/// Unsupported source domains, missing contact bindings, graph electrical work,
/// event cuts, budgets and cancellation refuse without publishing either state.
/// Independent mesh and timestep refinement remain necessary.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn step_spatial_moving_conduction_candidate<
    't,
    A: HeatAssembly + ?Sized,
    C: SpatialConductorMaterial,
>(
    heat: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    problem: &SpatialConduction<C>,
    boundary: &SpatialConductionBoundary,
    motion: SpatialConductionMotion<'_>,
    contacts: &[contacts::ConductionContactThermal<'_>],
    regions: &[SpatialConductionRegion],
    volume_tolerance_m3: f64,
    assembly: WeightedAssemblyLimits,
    field_controls: ConductionControls,
    controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
) -> Result<
    SpatialMovingConductionEnthalpyCandidate<'t, C, A::Sample>,
    CircuitEnthalpyError<ConductionHeatError<A::Error>>,
> {
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    check(&cancel).map_err(map_conduction_error)?;
    let prepared = heat.prepared_heat_network();
    let bulk = prepared
        .thermal_material_network()
        .map_err(assembly_heat_error)?;
    prepared
        .admit_field_independent_boundaries(cx, &[])
        .map_err(assembly_heat_error)?;
    let layout = SpatialConductionThermalNetwork {
        cells: bulk.cells,
        links: bulk.links,
        regions,
        volume_tolerance_m3,
    };
    let next = layout.step_moving_resolved_candidate(
        thermal,
        target_s,
        problem,
        boundary,
        motion,
        contacts,
        assembly,
        field_controls,
        controls,
        initial_energy_increments_j,
        |trial| {
            let coefficients = heat
                .resolve_heat(cx, trial)
                .map_err(|error| error.map_source(ConductionHeatError::Heat))?;
            if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
                return Err(map_conduction_error(invalid(
                    "moving 3D current cannot silently omit graph electrical work",
                )));
            }
            Ok(coefficients.bulk)
        },
        |sample: &A::Sample| &sample.as_ref().material,
        &cancel,
    )?;
    let report = next.report();
    closed_balance(
        target_s - thermal.time_s(),
        report.total_energy_residual_j,
        report.source.as_ref(),
        controls.total_energy_tolerance_j,
    )
    .map_err(assembly_heat_error)?;
    check(&cancel).map_err(map_conduction_error)?;
    Ok(next)
}
