//! Existing reciprocal finite-cell correlations with solved graph velocity.
//!
//! Geometry and bulk/film convention stay explicit. The source graph fixes
//! pressure/composition and binds the exact native fluid caloric object. Forced
//! velocity is |Q_reference|/A_native; original bulk/film rho, mu, k and Cp supply
//! the correlation. The Boussinesq reference density is not substituted into the
//! original material tuple. Natural and aiding mixed convection retain sourced
//! expansion and the selected geometry's original regime. The mixed plate
//! binds explicit upward orientation to actual gravity and signed native cell
//! velocity. No generic blend, speed floor or empirical entrance length is used.

use super::*;
use crate::graph::thermal::medium::convection::{
    AidingVerticalPlateFrame, ForcedConvectionSettings, MaterialConvectionError,
    NaturalConvectionGeometry,
    deposition::{
        CellConvectionLaw, evaluate_convection_state,
        source::{
            ConvectionPropertyTemperature, SourceCellConvectionError, SourceCellConvectionLaw,
            links::{
                ResolvedLinkConvection, SourceConvectiveLink, evaluate_source_link_convection,
            },
        },
    },
    evaluate_forced_from_material,
};
use fs_convection::{CorrelationId, ThermalDirection};
use fs_material::fluid_transport::ThermalFluidTransport;
use fs_qty::Length;

/// Authored wall geometry. A forced link's fluid cell must be a native pipe
/// section of this branch; a natural link may use a bound junction or pipe cell.
/// Pure natural cards require a zero old/new/average bulk-flow envelope on that
/// pipe (or every incident branch of a junction). Nonzero imposed flow needs
/// an admitted mixed relation; no arbitrary small-speed threshold is used.
#[derive(Clone, Copy, Debug)]
pub enum GraphConvectionLaw {
    Forced {
        branch: usize,
        correlation: CorrelationId,
        settings: ForcedConvectionSettings,
    },
    Natural {
        geometry: NaturalConvectionGeometry,
    },
    /// Original Manning-Qureshi unconfined isothermal laminar vertical plate.
    /// The cell must belong to `branch`; actual signed Q/A follows that native
    /// cell's start/end tangent. `plate_up` is authored in the geometry/gravity
    /// world frame. Neither inclined flow nor confinement is admitted by this
    /// card. Characteristic length remains the authored full plate height.
    /// Mean-film properties and sourced expansion are required. The caller
    /// must admit the original constant-property and free-stream assumptions.
    AidingVerticalPlate {
        branch: usize,
        plate_up: [f64; 3],
        maximum_relative_density_departure: f64,
    },
}

/// A native zero-conductance contact; native link area/coverage remain the
/// authority. Characteristic/developing lengths stay authored correlation input.
#[derive(Clone, Copy, Debug)]
pub struct GraphConvectiveLink {
    pub link: u64,
    pub fluid_cell: CellId,
    pub characteristic_length: Length,
    pub property_temperature: ConvectionPropertyTemperature,
    pub law: GraphConvectionLaw,
}

#[derive(Debug)]
pub enum GraphConvectionError {
    Source(SourceGraphError),
    Convection(SourceCellConvectionError<SourceGraphError>),
    Geometry {
        link: u64,
        source: MaterialConvectionError,
    },
    Envelope {
        link: u64,
        source: MaterialConvectionError,
    },
}
impl core::fmt::Display for GraphConvectionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite graph source convection: {self:?}")
    }
}
impl std::error::Error for GraphConvectionError {}

/// Hydraulic source admission and source-retaining reciprocal conductances.
/// Apply convection to the native step's link-conductance map. Do not also add
/// its equal/opposite powers as cell sources: the native link accounts for both.
#[derive(Debug)]
pub struct GraphSourceHeat<M> {
    pub admission: GraphSourceAdmission,
    pub convection: ResolvedLinkConvection<M>,
}

// Use this native cell's tangent, including each local section of a curved
// path. A whole-branch endpoint displacement cannot stand in for wall flow.
fn mixed_geometry(
    core: &BoundGraphCore<'_>,
    thermal: &DepositionSimulation<'_>,
    binding: &GraphConvectiveLink,
    branch: usize,
    plate_up: [f64; 3],
) -> Result<(AidingVerticalPlateFrame, [f64; 3], f64), GraphConvectionError> {
    let invalid = |reason| GraphConvectionError::Source(SourceGraphError::Invalid(reason));
    let Some(Some((actual_branch, section))) = core.locations.get(&binding.fluid_cell) else {
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
        .get(&binding.fluid_cell)
        .ok_or_else(|| invalid("missing original mixed-convection cell geometry"))?;
    let displacement =
        core::array::from_fn(|axis| cell.geometry.end_m[axis] - cell.geometry.start_m[axis]);
    let length = norm(displacement);
    if !positive(length) {
        return Err(invalid("mixed-convection cell tangent is unrepresentable"));
    }
    let frame =
        AidingVerticalPlateFrame::try_new(plate_up, core.acceleration).map_err(|source| {
            GraphConvectionError::Geometry {
                link: binding.link,
                source,
            }
        })?;
    Ok((
        frame,
        displacement.map(|value| value / length),
        core.section_areas[branch][*section],
    ))
}

#[allow(clippy::too_many_arguments)]
pub(super) fn evaluate_links<M: ThermalFluidTransport>(
    core: &BoundGraphCore<'_>,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'_>,
    flow: GraphFlow<'_>,
    bindings: &[GraphConvectiveLink],
    maximum_links: usize,
    mut resolve: impl FnMut(Temperature) -> Result<M, SourceGraphError>,
) -> Result<ResolvedLinkConvection<M>, GraphConvectionError> {
    let source_error = GraphConvectionError::Source;
    poll(cx).map_err(source_error)?;
    if !(1..=1_000_000).contains(&maximum_links) || bindings.len() > maximum_links {
        return Err(source_error(SourceGraphError::Invalid(
            "bounded finite convection link count required",
        )));
    }
    let gravity = norm(core.acceleration);
    if !nonnegative(gravity) {
        return Err(source_error(SourceGraphError::Invalid(
            "convection gravity magnitude is unrepresentable",
        )));
    }
    let mut native = Vec::with_capacity(bindings.len());
    for binding in bindings {
        poll(cx).map_err(source_error)?;
        let location = core.locations.get(&binding.fluid_cell).ok_or_else(|| {
            source_error(SourceGraphError::Cell {
                cell: binding.fluid_cell,
                reason: "convective fluid cell must belong to the bound original source",
            })
        })?;
        let law = match binding.law {
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
                        .position(|&id| id == binding.fluid_cell)
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
                    mixed_geometry(core, thermal, binding, branch, plate_up)?;
                let speed = flow.average(branch) / area;
                let settings = frame
                    .settings(
                        tangent.map(|value| speed * value),
                        maximum_relative_density_departure,
                    )
                    .map_err(|source| GraphConvectionError::Geometry {
                        link: binding.link,
                        source,
                    })?;
                SourceCellConvectionLaw::AidingVerticalPlate { settings }
            }
        };
        native.push(SourceConvectiveLink {
            link: binding.link,
            fluid_cell: binding.fluid_cell,
            characteristic_length: binding.characteristic_length,
            property_temperature: binding.property_temperature,
            law,
        });
    }
    let result = evaluate_source_link_convection(
        thermal,
        &native,
        maximum_links,
        |binding, _, temperature| {
            poll(cx)?;
            let material = resolve(temperature)?;
            core.admit_density(
                binding.fluid_cell,
                material
                    .transport_state()
                    .map_err(SourceGraphError::Material)?,
            )?;
            Ok(material)
        },
        || cx.checkpoint().is_err(),
    )
    .map_err(GraphConvectionError::Convection)?;
    // Heat conductance uses the actual signed-average velocity. Independently
    // qualify each moving-flow card at the same source's interval envelope; a
    // reversal must not conceal an out-of-regime or opposing-flow endpoint.
    // The native owner returns samples in link-ID order, independently of the
    // caller's binding order. Keep each card paired with its own source sample.
    let mut ordered: Vec<_> = bindings.iter().collect();
    ordered.sort_by_key(|binding| binding.link);
    if ordered.len() != result.samples().len() {
        return Err(source_error(SourceGraphError::Invalid(
            "native convection samples do not retain every bound link",
        )));
    }
    for (binding, sample) in ordered.into_iter().zip(result.samples()) {
        if binding.link != sample.link {
            return Err(source_error(SourceGraphError::Invalid(
                "native convection sample and binding identities differ",
            )));
        }
        poll(cx).map_err(source_error)?;
        match binding.law {
            GraphConvectionLaw::Forced {
                branch,
                correlation,
                mut settings,
            } => {
                let Some(Some((_, section))) = core.locations.get(&binding.fluid_cell) else {
                    return Err(source_error(SourceGraphError::Invalid(
                        "missing admitted pipe heat geometry",
                    )));
                };
                let initial = flow.initial(branch);
                let final_flow = flow.final_flow(branch);
                let area = core.section_areas[branch][*section];
                // Same actual thermal-direction convention as the native link.
                if sample.wall_temperature_k > sample.fluid_temperature_k {
                    settings.direction = ThermalDirection::HeatingFluid;
                } else if sample.wall_temperature_k < sample.fluid_temperature_k {
                    settings.direction = ThermalDirection::CoolingFluid;
                }
                let admit = |flow: f64| {
                    evaluate_forced_from_material(
                        correlation,
                        &sample.material,
                        flow.abs() / area,
                        binding.characteristic_length,
                        settings,
                    )
                    .map_err(|source| GraphConvectionError::Envelope {
                        link: binding.link,
                        source,
                    })
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
                    mixed_geometry(core, thermal, binding, branch, plate_up)?;
                let expansion_per_k =
                    sample.thermal.volumetric_expansion_per_k().ok_or_else(|| {
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
                        .and_then(|settings| {
                            // This is the native finite-cell owner's admission,
                            // including exact resting thermal equilibrium. Reuse the
                            // same film state, not a forced-only surrogate or a second
                            // property query with different pressure/composition.
                            evaluate_convection_state(
                                &sample.thermal.transport(),
                                sample.wall_temperature_k,
                                Temperature::new(sample.fluid_temperature_k),
                                binding.characteristic_length,
                                CellConvectionLaw::AidingVerticalPlate {
                                    settings,
                                    expansion_per_k,
                                },
                            )
                        });
                    admission.map_err(|source| GraphConvectionError::Envelope {
                        link: binding.link,
                        source,
                    })
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
    }
    Ok(result)
}
