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

mod laws;
/// Original flow and material laws on prescribed wall reservoirs.
pub mod reservoir;
/// Section-integrated entrance heat on original finite walls and reservoirs.
pub mod developing;
/// Source-integrated solid paths and actual native wall/fluid surface films.
pub mod conjugate;

use super::*;
use crate::graph::thermal::medium::convection::{
    AidingVerticalPlateFrame, ForcedConvectionSettings, MaterialConvectionError,
    NaturalConvectionGeometry,
    deposition::{
        CellConvectionLaw, evaluate_convection_state,
        source::{
            ConvectionPropertyTemperature, SourceCellConvectionError, SourceCellConvectionLaw,
            links::{
                ConvectiveLinkSample, ResolvedLinkConvection, SourceConvectiveLink,
                evaluate_source_link_at_surface, evaluate_source_link_convection,
            },
        },
    },
    evaluate_flat_plate_section_from_material, evaluate_forced_from_material,
    evaluate_sphere_from_material,
};
use fs_convection::{CorrelationId, ThermalDirection};
use fs_material::fluid_transport::ThermalFluidTransport;
use fs_qty::Length;

/// Direction away from an external plate's physical leading edge, relative to
/// the native branch's `from` -> `to` path. Reversing flow requires rebinding
/// both this direction and the section distances from the opposite edge.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PlateFlowDirection {
    AlongBranch,
    AgainstBranch,
}

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
    /// Original Whitaker isolated stationary isothermal sphere in external flow.
    /// The binding's characteristic length is the whole sphere diameter.
    /// Select `Reservoir` properties: Re, Pr and k use the actual bulk state;
    /// the same source is also queried at the wall solely for its viscosity.
    /// Original Reynolds, Prandtl and viscosity-ratio bounds apply at every
    /// native flow endpoint and any interior zero crossing. No stagnant-flow,
    /// natural-convection or confinement correction is inferred.
    ForcedSphere {
        branch: usize,
    },
    /// An external isothermal plate strip using the original laminar or
    /// leading-edge-corrected turbulent mean. `start` and the binding's
    /// `characteristic_length` are upstream/downstream distances from ONE
    /// physical leading edge, not cell widths or hydraulic diameters.
    /// Mean-film properties are required. The initial, final and average native
    /// flow must follow `direction`; an interval reversal refuses atomically.
    /// Each positive endpoint retains the same card's regime admission.
    /// The caller admits attached, unconfined, zero-pressure-gradient flow;
    /// locally updated properties do not resolve a variable-property boundary layer.
    ForcedPlateSection {
        branch: usize,
        correlation: CorrelationId,
        start: Length,
        direction: PlateFlowDirection,
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
    let gravity = laws::gravity(core).map_err(source_error)?;
    let mut native = Vec::with_capacity(bindings.len());
    for binding in bindings {
        poll(cx).map_err(source_error)?;
        native.push(native_link(core, thermal, flow, binding, gravity)?);
    }
    let result = evaluate_source_link_convection(
        thermal,
        &native,
        maximum_links,
        |binding, _, temperature| link_material(core, cx, binding, temperature, &mut resolve),
        || cx.checkpoint().is_err(),
    )
    .map_err(GraphConvectionError::Convection)?;
    // Native samples are in link-ID order, independent of caller order.
    let mut ordered: Vec<_> = bindings.iter().collect();
    ordered.sort_by_key(|binding| binding.link);
    if ordered.len() != result.samples().len() {
        return Err(source_error(SourceGraphError::Invalid(
            "native convection samples do not retain every bound link",
        )));
    }
    for (binding, sample) in ordered.into_iter().zip(result.samples()) {
        admit_link_sample(core, cx, thermal, flow, binding, sample)?;
    }
    Ok(result)
}

fn native_link(
    core: &BoundGraphCore<'_>,
    thermal: &DepositionSimulation<'_>,
    flow: GraphFlow<'_>,
    binding: &GraphConvectiveLink,
    gravity: f64,
) -> Result<SourceConvectiveLink, GraphConvectionError> {
    let law = laws::resolve(
        core,
        thermal,
        flow,
        binding.fluid_cell,
        binding.law,
        gravity,
    )
    .map_err(|error| link_law_error(binding.link, error))?;
    Ok(SourceConvectiveLink {
        link: binding.link,
        fluid_cell: binding.fluid_cell,
        characteristic_length: binding.characteristic_length,
        property_temperature: binding.property_temperature,
        law,
    })
}

fn link_material<M: ThermalFluidTransport>(
    core: &BoundGraphCore<'_>,
    cx: &Cx<'_>,
    binding: &SourceConvectiveLink,
    temperature: Temperature,
    resolve: &mut impl FnMut(Temperature) -> Result<M, SourceGraphError>,
) -> Result<M, SourceGraphError> {
    poll(cx)?;
    let material = resolve(temperature)?;
    if !matches!(binding.law, SourceCellConvectionLaw::ForcedSphere { .. }) {
        core.admit_density(
            binding.fluid_cell,
            material
                .transport_state()
                .map_err(SourceGraphError::Material)?,
        )?;
    }
    // A sphere's additional query supplies only surface viscosity. Its density
    // is not another stored fluid state; admit the primary bulk snapshot below.
    Ok(material)
}

fn admit_link_sample<M: ThermalFluidTransport>(
    core: &BoundGraphCore<'_>,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'_>,
    flow: GraphFlow<'_>,
    binding: &GraphConvectiveLink,
    sample: &ConvectiveLinkSample<M>,
) -> Result<(), GraphConvectionError> {
    if binding.link != sample.link {
        return Err(GraphConvectionError::Source(SourceGraphError::Invalid(
            "native convection sample and binding identities differ",
        )));
    }
    poll(cx).map_err(GraphConvectionError::Source)?;
    if matches!(binding.law, GraphConvectionLaw::ForcedSphere { .. }) {
        core.admit_density(binding.fluid_cell, sample.thermal.transport())
            .map_err(GraphConvectionError::Source)?;
    }
    laws::admit_envelope(
        core,
        thermal,
        flow,
        laws::FlowConvectionPoint {
            fluid_cell: binding.fluid_cell,
            characteristic_length: binding.characteristic_length,
            law: binding.law,
            wall_temperature_k: sample.wall_temperature_k,
            fluid_temperature_k: sample.fluid_temperature_k,
            material: &sample.material,
            thermal: &sample.thermal,
            surface: sample.surface.as_ref().map(|surface| &surface.thermal),
        },
    )
    .map_err(|error| link_law_error(binding.link, error))
}

// A real material surface uses the SAME fluid source and flow-envelope owner.
// The result is one diagnostic sample; only the surface solver publishes the
// combined reciprocal coefficient to the original thermal link.
#[allow(clippy::too_many_arguments)]
fn evaluate_link_surface<M: ThermalFluidTransport>(
    core: &BoundGraphCore<'_>,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'_>,
    flow: GraphFlow<'_>,
    binding: &GraphConvectiveLink,
    native: &SourceConvectiveLink,
    surface_temperature_k: f64,
    resolve: &mut impl FnMut(Temperature) -> Result<M, SourceGraphError>,
) -> Result<ConvectiveLinkSample<M>, GraphConvectionError> {
    let sample = evaluate_source_link_at_surface(
        thermal,
        native,
        surface_temperature_k,
        |binding, _, temperature| link_material(core, cx, binding, temperature, resolve),
        || cx.checkpoint().is_err(),
    )
    .map_err(GraphConvectionError::Convection)?;
    admit_link_sample(core, cx, thermal, flow, binding, &sample)?;
    Ok(sample)
}

fn link_law_error(link: u64, error: laws::FlowConvectionError) -> GraphConvectionError {
    match error {
        laws::FlowConvectionError::Source(error) => GraphConvectionError::Source(error),
        laws::FlowConvectionError::Geometry(source) => GraphConvectionError::Geometry { link, source },
        laws::FlowConvectionError::Envelope(source) => GraphConvectionError::Envelope { link, source },
    }
}
