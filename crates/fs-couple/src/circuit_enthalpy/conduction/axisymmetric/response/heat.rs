//! Complete sourced heat transport around the original meridional current field.

use super::*;
use crate::circuit_enthalpy::{
    conduction::{ConductionRegion, heat::ConductionHeatError},
    contacts::heat::{HeatAssembly, HeatNetworkError, closed_balance},
};
use fs_exec::Cx;
use std::collections::BTreeMap;

/// Held boundary or grouped-electrode drive on one axisymmetric heat interval.
/// The original native solver determines its actual nonlinear operating point.
#[derive(Clone, Copy)]
pub enum AxisymmetricConductionResponseDrive<'a> {
    Field(ConductionDrive<'a>),
    Terminals(ConductionTerminalDrive<'a>),
}
impl<'a> From<ConductionDrive<'a>> for AxisymmetricConductionResponseDrive<'a> {
    fn from(drive: ConductionDrive<'a>) -> Self {
        Self::Field(drive)
    }
}
impl<'a> From<ConductionTerminalDrive<'a>> for AxisymmetricConductionResponseDrive<'a> {
    fn from(drive: ConductionTerminalDrive<'a>) -> Self {
        Self::Terminals(drive)
    }
}

fn heat_error<E: From<HeatNetworkError>>(
    error: HeatNetworkError,
) -> CircuitEnthalpyError<ConductionHeatError<E>> {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(ConductionHeatError::Heat(other.into())),
    }
}

/// Close nonlinear swept bulk/contact current and any complete heat assembly
/// at the SAME finite-enthalpy endpoint.
///
/// Original thermal contacts, conduction layers, enclosure radiation and finite
/// coolant links see the same proposed temperatures as the original electrical
/// laws. Native swept-volume and swept-surface losses enter their original
/// receiving cells once. Internal Fourier and enclosure transfers remain
/// internal; the report retains the complete typed heat source sample.
///
/// Existing voltage, current-density, total-current, regulated-power and
/// grouped voltage/current/mixed terminal drives all use the original native
/// nonlinear field solver. Its actual field solves consume the shared thermal
/// work budget, including rejected constitutive trials and final admission.
///
/// Source identities, phase charts, fixed geometry, complete swept-volume
/// routes, sampling policies and electrical-contact partitions remain owned by
/// the existing field and heat providers. Thermal source response/applicability
/// must be independent of electrical operating coordinates. Graph electrical
/// bindings require their own work owner and are refused. Neither endpoint is
/// published on source failure, cancellation or exhausted numerical budgets.
///
/// This first-order quasistatic coupling still requires independent time, mesh
/// and contact-integration refinement. It adds no reversible electrical storage
/// or surface-motion law.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn step_axisymmetric_conduction_response_candidate<'t, 'd, A: HeatAssembly + ?Sized>(
    heat: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    problem: &AxisymmetricConductionResponse,
    drive: impl Into<AxisymmetricConductionResponseDrive<'d>>,
    contacts: &[contacts::ConductionContactThermal<'_>],
    regions: &[ConductionRegion],
    volume_tolerance_m3: f64,
    assembly: WeightedAssemblyLimits,
    field_controls: ConductionControls,
    response: ConductionResponseControls,
    controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
) -> Result<
    AxisymmetricConductionResponseEnthalpyCandidate<'t, A::Sample>,
    CircuitEnthalpyError<ConductionHeatError<A::Error>>,
> {
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    check(&cancel).map_err(map_conduction_error)?;
    let prepared = heat.prepared_heat_network();
    let bulk = prepared.thermal_material_network().map_err(heat_error)?;
    if regions.len() != problem.sources().len() || regions.len() > assembly.max_cells {
        return Err(map_conduction_error(invalid(
            "one bounded heat route per original meridional triangle required",
        )));
    }
    let mut routes = BTreeMap::new();
    let mut field_axes = Vec::new();
    for region in regions {
        check(&cancel).map_err(map_conduction_error)?;
        let source = problem.sources().get(region.triangle).ok_or_else(|| {
            map_conduction_error(invalid("heat route names an absent meridional conductor"))
        })?;
        if routes.insert(region.triangle, region.cell).is_some() {
            return Err(map_conduction_error(invalid(
                "original meridional heat routes must be unique",
            )));
        }
        for (axis, _) in source.axes() {
            field_axes.push((region.cell, axis.as_str()));
        }
    }
    if let Some(sources) = problem.contact_response_sources() {
        let trace = problem.mesh().boundary_trace().map_err(|_| {
            map_conduction_error(invalid(
                "meridional contact lost its original exterior trace",
            ))
        })?;
        let edge_key = |mut vertices: [usize; 2]| {
            vertices.sort_unstable();
            vertices
        };
        let mut adjacent = BTreeMap::new();
        for edge in &trace.edges {
            check(&cancel).map_err(map_conduction_error)?;
            adjacent.insert(
                edge_key(edge.oriented_vertices.map(|vertex| vertex as usize)),
                edge.source_face,
            );
        }
        for (law, contact) in sources.laws().iter().zip(problem.contacts()) {
            check(&cancel).map_err(map_conduction_error)?;
            for vertices in [contact.first_vertices, contact.second_vertices] {
                let cell = adjacent
                    .get(&edge_key(vertices))
                    .and_then(|face| routes.get(face))
                    .ok_or_else(|| {
                        map_conduction_error(invalid(
                            "nonlinear meridional contact lacks its original heat route",
                        ))
                    })?;
                for (axis, _) in law.axes() {
                    field_axes.push((*cell, axis.as_str()));
                }
            }
        }
    }
    prepared
        .admit_field_independent_boundaries(cx, &field_axes)
        .map_err(heat_error)?;
    let layout = ConductionThermalNetwork {
        cells: bulk.cells,
        links: bulk.links,
        regions,
        volume_tolerance_m3,
    };
    let drive = match drive.into() {
        AxisymmetricConductionResponseDrive::Field(drive) => Drive::Field(drive),
        AxisymmetricConductionResponseDrive::Terminals(drive) => Drive::Terminals(drive),
    };
    let next = layout.step_axisymmetric_response_resolved_candidate(
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
        |trial| {
            let coefficients = heat
                .resolve_heat(cx, trial)
                .map_err(|error| error.map_source(ConductionHeatError::Heat))?;
            if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
                return Err(map_conduction_error(invalid(
                    "nonlinear meridional current cannot omit graph electrical work",
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
    .map_err(heat_error)?;
    check(&cancel).map_err(map_conduction_error)?;
    Ok(next)
}
