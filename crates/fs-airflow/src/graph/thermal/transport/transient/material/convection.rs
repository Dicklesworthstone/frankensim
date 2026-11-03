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
        let law = laws::resolve(
            core, thermal, flow, binding.fluid_cell, binding.law, gravity,
        ).map_err(|error| link_law_error(binding.link, error))?;
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
    // Native samples are in link-ID order, independent of caller order.
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
            },
        ).map_err(|error| link_law_error(binding.link, error))?;
    }
    Ok(result)
}

fn link_law_error(link: u64, error: laws::FlowConvectionError) -> GraphConvectionError {
    match error {
        laws::FlowConvectionError::Source(error) => GraphConvectionError::Source(error),
        laws::FlowConvectionError::Geometry(source) => GraphConvectionError::Geometry { link, source },
        laws::FlowConvectionError::Envelope(source) => GraphConvectionError::Envelope { link, source },
    }
}
