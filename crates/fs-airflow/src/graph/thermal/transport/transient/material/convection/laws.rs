//! Shared native-flow admission for finite walls and prescribed wall reservoirs.
//! Geometry, source state and signed interval rates keep the same physical owner.

use super::*;
use fs_material::fluid_transport::ThermalFluidTransportState;

#[derive(Debug)]
pub(super) enum FlowConvectionError {
    Source(SourceGraphError),
    Geometry(MaterialConvectionError),
    Envelope(MaterialConvectionError),
}

pub(super) struct FlowConvectionPoint<'a, M> {
    pub fluid_cell: CellId,
    pub characteristic_length: Length,
    pub law: GraphConvectionLaw,
    pub wall_temperature_k: f64,
    pub fluid_temperature_k: f64,
    pub material: &'a M,
    pub thermal: &'a ThermalFluidTransportState,
}

pub(super) fn gravity(core: &BoundGraphCore<'_>) -> Result<f64, SourceGraphError> {
    let gravity = norm(core.acceleration);
    if !nonnegative(gravity) {
        return Err(SourceGraphError::Invalid(
            "convection gravity magnitude is unrepresentable",
        ));
    }
    Ok(gravity)
}

// Use this native cell's tangent, including each local section of a curved
// path. A whole-branch endpoint displacement cannot stand in for wall flow.
fn mixed_geometry(
    core: &BoundGraphCore<'_>,
    thermal: &DepositionSimulation<'_>,
    fluid_cell: CellId,
    branch: usize,
    plate_up: [f64; 3],
) -> Result<(AidingVerticalPlateFrame, [f64; 3], f64), FlowConvectionError> {
    let invalid = |reason| FlowConvectionError::Source(SourceGraphError::Invalid(reason));
    let Some(Some((actual_branch, section))) = core.locations.get(&fluid_cell) else {
        return Err(invalid(
            "mixed graph convection requires a native pipe cell",
        ));
    };
    if *actual_branch != branch {
        return Err(invalid(
            "mixed convection branch does not own the fluid cell",
        ));
    }
    let cell = thermal
        .cells()
        .get(&fluid_cell)
        .ok_or_else(|| invalid("missing original mixed-convection cell geometry"))?;
    let displacement =
        core::array::from_fn(|axis| cell.geometry.end_m[axis] - cell.geometry.start_m[axis]);
    let length = norm(displacement);
    if !positive(length) {
        return Err(invalid("mixed-convection cell tangent is unrepresentable"));
    }
    let frame = AidingVerticalPlateFrame::try_new(plate_up, core.acceleration)
        .map_err(|error| FlowConvectionError::Geometry(error.into()))?;
    Ok((
        frame,
        displacement.map(|value| value / length),
        core.section_areas[branch][*section],
    ))
}

pub(super) fn resolve(
    core: &BoundGraphCore<'_>,
    thermal: &DepositionSimulation<'_>,
    flow: GraphFlow<'_>,
    fluid_cell: CellId,
    authored: GraphConvectionLaw,
    gravity: f64,
) -> Result<SourceCellConvectionLaw, FlowConvectionError> {
    let source_error = FlowConvectionError::Source;
    let location = core.locations.get(&fluid_cell).ok_or_else(|| {
        source_error(SourceGraphError::Cell {
            cell: fluid_cell,
            reason: "convective fluid cell must belong to the bound original source",
        })
    })?;
    let law = match authored {
        GraphConvectionLaw::Forced {
            branch,
            correlation,
            settings,
        } => {
            let &(actual_branch, section) = location.as_ref().ok_or_else(|| {
                source_error(SourceGraphError::Invalid(
                    "forced graph convection requires a pipe cell with native flow area",
                ))
            })?;
            if actual_branch != branch {
                return Err(source_error(SourceGraphError::Invalid(
                    "forced convection branch does not own the fluid cell",
                )));
            }
            let speed = flow.average(branch).abs() / core.section_areas[branch][section];
            if !nonnegative(speed) {
                return Err(source_error(SourceGraphError::Invalid(
                    "actual graph convection speed is unrepresentable",
                )));
            }
            SourceCellConvectionLaw::Forced {
                correlation,
                speed_m_s: speed,
                settings,
            }
        }
        GraphConvectionLaw::Natural { geometry } => {
            let quiescent = |branch: usize| {
                flow.initial(branch) == 0.0
                    && flow.final_flow(branch) == 0.0
                    && flow.average(branch) == 0.0
            };
            let admitted = if let Some((branch, _)) = location {
                quiescent(*branch)
            } else {
                let node = core
                    .inventory
                    .node_cells()
                    .iter()
                    .position(|&id| id == fluid_cell)
                    .ok_or_else(|| {
                        source_error(SourceGraphError::Invalid(
                            "missing original junction identity",
                        ))
                    })?;
                core.template
                    .branches()
                    .iter()
                    .enumerate()
                    .filter(|(_, edge)| edge.from == node || edge.to == node)
                    .all(|(branch, _)| quiescent(branch))
            };
            if !admitted {
                return Err(source_error(SourceGraphError::Invalid(
                    "pure natural wall convection requires a quiescent native speed envelope",
                )));
            }
            SourceCellConvectionLaw::Natural {
                geometry,
                gravity_m_s2: gravity,
            }
        }
        GraphConvectionLaw::AidingVerticalPlate {
            branch,
            plate_up,
            maximum_relative_density_departure,
        } => {
            let (frame, tangent, area) =
                mixed_geometry(core, thermal, fluid_cell, branch, plate_up)?;
            let speed = flow.average(branch) / area;
            let settings = frame
                .settings(
                    tangent.map(|value| speed * value),
                    maximum_relative_density_departure,
                )
                .map_err(|error| FlowConvectionError::Geometry(error.into()))?;
            SourceCellConvectionLaw::AidingVerticalPlate { settings }
        }
    };
    Ok(law)
}

pub(super) fn admit_envelope<M: ThermalFluidTransport>(
    core: &BoundGraphCore<'_>,
    thermal: &DepositionSimulation<'_>,
    flow: GraphFlow<'_>,
    point: FlowConvectionPoint<'_, M>,
) -> Result<(), FlowConvectionError> {
    let source_error = FlowConvectionError::Source;
    match point.law {
        GraphConvectionLaw::Forced {
            branch,
            correlation,
            mut settings,
        } => {
            let Some(Some((_, section))) = core.locations.get(&point.fluid_cell) else {
                return Err(source_error(SourceGraphError::Invalid(
                    "missing admitted pipe heat geometry",
                )));
            };
            let initial = flow.initial(branch);
            let final_flow = flow.final_flow(branch);
            let area = core.section_areas[branch][*section];
            // Same actual thermal-direction convention as the native link.
            if point.wall_temperature_k > point.fluid_temperature_k {
                settings.direction = ThermalDirection::HeatingFluid;
            } else if point.wall_temperature_k < point.fluid_temperature_k {
                settings.direction = ThermalDirection::CoolingFluid;
            }
            let admit = |flow: f64| {
                evaluate_forced_from_material(
                    correlation,
                    point.material,
                    flow.abs() / area,
                    point.characteristic_length,
                    settings,
                )
                .map_err(FlowConvectionError::Envelope)
            };
            // Native evaluation already admitted the actual average. Check
            // both endpoints, including lower Reynolds limits, and the minimum
            // speed within a signed reversal of the held momentum trajectory.
            admit(initial)?;
            admit(final_flow)?;
            if (initial < 0.0 && final_flow > 0.0) || (initial > 0.0 && final_flow < 0.0) {
                admit(0.0)?;
            }
        }
        GraphConvectionLaw::AidingVerticalPlate {
            branch,
            plate_up,
            maximum_relative_density_departure,
        } => {
            let (frame, tangent, area) =
                mixed_geometry(core, thermal, point.fluid_cell, branch, plate_up)?;
            let expansion_per_k =
                point.thermal.volumetric_expansion_per_k().ok_or_else(|| {
                    source_error(SourceGraphError::Invalid(
                        "mixed convection requires the retained original source expansion",
                    ))
                })?;
            let admit = |volume_flow: f64| {
                let speed = volume_flow / area;
                let admission = frame
                    .settings(
                        tangent.map(|value| speed * value),
                        maximum_relative_density_departure,
                    )
                    .map_err(MaterialConvectionError::from)
                    .and_then(|settings| {
                        // This is the native finite-cell owner's admission,
                        // including exact resting thermal equilibrium. Reuse the
                        // same film state, not a forced-only surrogate or a second
                        // property query with different pressure/composition.
                        evaluate_convection_state(
                            &point.thermal.transport(),
                            point.wall_temperature_k,
                            Temperature::new(point.fluid_temperature_k),
                            point.characteristic_length,
                            CellConvectionLaw::AidingVerticalPlate {
                                settings,
                                expansion_per_k,
                            },
                        )
                    });
                admission.map_err(FlowConvectionError::Envelope)
            };
            let initial = flow.initial(branch);
            let final_flow = flow.final_flow(branch);
            admit(initial)?;
            admit(final_flow)?;
            if (initial < 0.0 && final_flow > 0.0) || (initial > 0.0 && final_flow < 0.0) {
                admit(0.0)?;
            }
        }
        GraphConvectionLaw::Natural { .. } => {}
    }
    Ok(())
}
