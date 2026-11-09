//! Complete sourced heat transport around the original nonlinear 3D field.

use super::*;
use crate::circuit_enthalpy::{
    conduction::heat::ConductionHeatError,
    contacts::heat::{HeatAssembly, HeatNetworkError, closed_balance},
};
use fs_exec::Cx;

fn heat_error<E: From<HeatNetworkError>>(
    error: HeatNetworkError,
) -> CircuitEnthalpyError<ConductionHeatError<E>> {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(ConductionHeatError::Heat(other.into())),
    }
}

/// Close source-dependent tetrahedral current, native bulk/contact Joule heat,
/// finite caloric storage and any complete [`HeatAssembly`] at one endpoint.
///
/// Thermal contacts, conduction layers, radiation and finite coolant links
/// retain their original source/routing owners. They see the same trial
/// temperatures as the electrical laws and feed back into the actual field.
/// Each loss enters storage once; enclosure radiation and reciprocal finite
/// links remain internal exchange. The complete typed heat sample is retained.
///
/// Original static geometry, source-domain, confined-phase, work-budget and
/// cancellation rules apply. Thermal laws must be independent of electrical
/// response and applicability; this owner supplies no reversible field work.
/// Other graph electrical sources refuse because their heat is not owned by
/// this standalone field. No rejected trial publishes either history.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn step_spatial_conduction_response_candidate<'t, 'd, A: HeatAssembly + ?Sized>(
    heat: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    problem: &SpatialConductionResponse,
    drive: impl Into<SpatialConductionResponseDrive<'d>>,
    contacts: &[ConductionContactThermal<'_>],
    regions: &[SpatialConductionRegion],
    volume_tolerance_m3: f64,
    assembly: WeightedAssemblyLimits,
    field_controls: ConductionControls,
    response: SpatialConductionResponseControls,
    controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
) -> Result<
    SpatialConductionResponseEnthalpyCandidate<'t, A::Sample>,
    CircuitEnthalpyError<ConductionHeatError<A::Error>>,
> {
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    step_with_heat(
        heat,
        cx,
        thermal,
        target_s,
        problem,
        regions,
        volume_tolerance_m3,
        assembly,
        controls,
        &cancel,
        |layout, resolve| {
            layout.step_response_resolved_candidate(
                thermal,
                target_s,
                problem,
                drive,
                contacts,
                assembly,
                field_controls,
                response,
                controls,
                initial_energy_increments_j,
                resolve,
                |sample: &A::Sample| &sample.as_ref().material,
                &cancel,
            )
        },
    )
}

/// Couple nonlinear prescribed-motion current to the complete original heat
/// assembly at one endpoint. Source laws see the actual conductor-frame E' and
/// J; thermal contacts, layers, radiation and finite coolant links see the same
/// trial temperature/phase. Native Joule heat enters the original cells once.
///
/// The typed heat sample and native moving/source-closure report are retained.
/// Signed terminal input minus signed Lorentz output is the thermal supply,
/// including generation and braking. Motion/B remain prescribed on the fixed
/// domain; the native stationary electrical-contact velocity guard still binds.
/// All source-domain, field-independent caloric/thermal, budget, cancellation
/// and atomic-continuation rules of the stationary shared owner remain active.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn step_spatial_moving_conduction_response_candidate<'t, 'd, A: HeatAssembly + ?Sized>(
    heat: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    problem: &SpatialConductionResponse,
    drive: impl Into<
        fs_em::spatial::conduction::response::motion::SpatialMovingConductionResponseDrive<'d>,
    >,
    motion: fs_em::spatial::conduction::motion::SpatialConductionMotion<'_>,
    contacts: &[ConductionContactThermal<'_>],
    regions: &[SpatialConductionRegion],
    volume_tolerance_m3: f64,
    assembly: WeightedAssemblyLimits,
    field_controls: ConductionControls,
    response: SpatialConductionResponseControls,
    controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
) -> Result<
    motion::SpatialMovingConductionResponseEnthalpyCandidate<'t, A::Sample>,
    CircuitEnthalpyError<ConductionHeatError<A::Error>>,
> {
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    step_with_heat(
        heat,
        cx,
        thermal,
        target_s,
        problem,
        regions,
        volume_tolerance_m3,
        assembly,
        controls,
        &cancel,
        |layout, resolve| {
            layout.step_moving_response_resolved_candidate(
                thermal,
                target_s,
                problem,
                drive,
                motion,
                contacts,
                assembly,
                field_controls,
                response,
                controls,
                initial_energy_increments_j,
                resolve,
                |sample: &A::Sample| &sample.as_ref().material,
                &cancel,
            )
        },
    )
}

// One complete-assembly source/axis admission and closed-balance owner for
// both native field frames. The callback only selects the electrical adapter.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub(super) fn step_with_heat<'t, A: HeatAssembly + ?Sized, F>(
    heat: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    problem: &SpatialConductionResponse,
    regions: &[SpatialConductionRegion],
    volume_tolerance_m3: f64,
    assembly: WeightedAssemblyLimits,
    controls: CircuitEnthalpySolve,
    cancel: &RefCell<impl FnMut() -> bool>,
    step: impl FnOnce(
        &SpatialConductionThermalNetwork<'_>,
        &dyn Fn(
            &DepositionSimulation<'t>,
        ) -> Result<
            CircuitThermalCoefficients<A::Sample>,
            CircuitEnthalpyError<ConductionHeatError<A::Error>>,
        >,
    ) -> Result<
        SpatialConductionResponseEnthalpyCandidate<'t, A::Sample, F>,
        CircuitEnthalpyError<ConductionHeatError<A::Error>>,
    >,
) -> Result<
    SpatialConductionResponseEnthalpyCandidate<'t, A::Sample, F>,
    CircuitEnthalpyError<ConductionHeatError<A::Error>>,
> {
    check(cancel).map_err(map_conduction_error)?;
    let prepared = heat.prepared_heat_network();
    let bulk = prepared.thermal_material_network().map_err(heat_error)?;
    if regions.len() != problem.sources().len() || regions.len() > assembly.max_cells {
        return Err(map_conduction_error(invalid(
            "one bounded heat route per original nonlinear tetrahedron required",
        )));
    }
    let mut field_axes = Vec::new();
    for region in regions {
        check(cancel).map_err(map_conduction_error)?;
        let source = problem.sources().get(region.tetrahedron).ok_or_else(|| {
            map_conduction_error(invalid("heat route names an absent tetrahedron"))
        })?;
        for (axis, _) in source.axes() {
            field_axes.push((region.cell, axis.as_str()));
        }
    }
    if let Some(sources) = problem.contact_response_sources() {
        for (law, adjacent) in sources.laws().iter().zip(problem.contact_cells()) {
            for &tetrahedron in adjacent {
                let region = regions
                    .iter()
                    .find(|r| r.tetrahedron == tetrahedron)
                    .ok_or_else(|| {
                        map_conduction_error(invalid("nonlinear contact lacks a thermal route"))
                    })?;
                for (axis, _) in law.axes() {
                    field_axes.push((region.cell, axis.as_str()));
                }
            }
        }
    }
    prepared
        .admit_field_independent_boundaries(cx, &field_axes)
        .map_err(heat_error)?;
    let layout = SpatialConductionThermalNetwork {
        cells: bulk.cells,
        links: bulk.links,
        regions,
        volume_tolerance_m3,
    };
    let resolve = |trial: &DepositionSimulation<'t>| {
        let coefficients = heat
            .resolve_heat(cx, trial)
            .map_err(|error| error.map_source(ConductionHeatError::Heat))?;
        if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
            return Err(map_conduction_error(invalid(
                "nonlinear 3D current cannot omit other graph electrical work",
            )));
        }
        Ok(coefficients.bulk)
    };
    let next = step(&layout, &resolve)?;
    let report = next.report();
    closed_balance(
        target_s - thermal.time_s(),
        report.total_energy_residual_j,
        report.source.as_ref(),
        controls.total_energy_tolerance_j,
    )
    .map_err(heat_error)?;
    check(cancel).map_err(map_conduction_error)?;
    Ok(next)
}
