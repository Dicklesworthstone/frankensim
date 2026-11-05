//! Original source films evaluated at the actual conducting-wall surface.

use super::super::super::{BoundGraphCore, GraphFlow, SourceGraphError, poll};
use super::super::{evaluate_link_surface, laws, native_link};
use super::*;
use fs_exec::Cx;
use fs_material::fluid_transport::ThermalFluidTransport;
use fs_material::thermal::conduction::{ConvectiveLayerError, ConvectiveLayerResponse};
use fs_qty::Temperature;

/// One caloric trace, one source-integrated material path and the actual film.
/// `layer.film.wall_temperature_k` is the solved surface temperature, distinct
/// from the original stored wall-cell temperature below. The native link uses
/// only the combined conductance; film hA remains source/correlation evidence.
#[derive(Debug)]
pub struct GraphConjugateSurfaceSample<M> {
    pub link: u64,
    pub wall_cell: CellId,
    pub fluid_cell: CellId,
    pub wall_temperature_k: f64,
    pub fluid_temperature_k: f64,
    pub geometry: ReservoirLayerGeometry,
    pub source_line: Arc<PropertyLineSlice>,
    pub layer: ConvectiveLayerResponse<ConvectiveLinkSample<M>>,
    pub applied_heat_into_fluid_w: f64,
    pub representation_error_w: f64,
    pub representation_allowance_w: f64,
}

/// Complete native hydraulic admission and one reciprocal heat coefficient
/// per original link. No bare-film conductance map is exposed for application.
#[derive(Debug)]
pub struct ResolvedGraphSurfaceConjugateHeat<M> {
    pub admission: GraphSourceAdmission,
    samples: Vec<GraphConjugateSurfaceSample<M>>,
    conductances_w_k: BTreeMap<u64, f64>,
}
impl<M> ResolvedGraphSurfaceConjugateHeat<M> {
    pub fn samples(&self) -> &[GraphConjugateSurfaceSample<M>] {
        &self.samples
    }
    pub fn conductances_w_k(&self) -> &BTreeMap<u64, f64> {
        &self.conductances_w_k
    }

    /// Apply the combined material/film link exactly once. Conflicting owners
    /// and cancellation refuse before changing the destination map.
    pub fn apply_to(
        &self,
        destination: &mut BTreeMap<u64, f64>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), GraphConjugateError> {
        apply_conductances(&self.conductances_w_k, destination, checkpoint)
    }
}

impl PreparedGraphConjugateLinks<'_, '_> {
    #[allow(clippy::too_many_arguments)]
    pub(in crate::graph::thermal::transport::transient::material) fn resolve_surface<
        M: ThermalFluidTransport,
    >(
        &self,
        core: &BoundGraphCore<'_>,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        flow: GraphFlow<'_>,
        admission: GraphSourceAdmission,
        config: LayerConductionConfig,
        mut resolve: impl FnMut(Temperature) -> Result<M, SourceGraphError>,
    ) -> Result<ResolvedGraphSurfaceConjugateHeat<M>, GraphConjugateError> {
        let source_error =
            |error| GraphConjugateError::Convection(GraphConvectionError::Source(error));
        poll(cx).map_err(source_error)?;
        let gravity = laws::gravity(core).map_err(source_error)?;
        // Validate all prepared storage/contact identities before fluid queries.
        for link in &self.links {
            poll(cx).map_err(source_error)?;
            link.temperatures(thermal)?;
        }
        let mut samples = Vec::with_capacity(self.links.len());
        let mut conductances_w_k = BTreeMap::new();
        for prepared in &self.links {
            poll(cx).map_err(source_error)?;
            let binding = &prepared.binding;
            let [wall, bulk] = prepared.temperatures(thermal)?;
            let native = native_link(core, thermal, flow, &binding.convection, gravity)
                .map_err(GraphConjugateError::Convection)?;
            let layer = binding
                .geometry
                .evaluate_convective_material_to_bath(
                    binding.wall_layer.conductivity_curve(),
                    [wall, bulk],
                    config,
                    |surface| {
                        let film = evaluate_link_surface(
                            core,
                            cx,
                            thermal,
                            flow,
                            &binding.convection,
                            &native,
                            surface,
                            &mut resolve,
                        )?;
                        Ok((film.heat_into_fluid_w, film.conductance_w_k, film))
                    },
                    || cx.checkpoint().map_err(|_| DepositionError::Cancelled),
                )
                .map_err(|error| surface_error(binding.convection.link, error))?;
            let difference = wall - bulk;
            let conductance = layer.conductance_w_k;
            let applied = conductance * difference;
            let error = applied - layer.power_to_bath_w;
            let representation =
                8.0 * f64::EPSILON * applied.abs().max(layer.power_to_bath_w.abs());
            if !conductance.is_finite()
                || conductance < 0.0
                || !applied.is_finite()
                || !error.is_finite()
                || error.abs() > representation
                || (difference != 0.0 && applied == 0.0)
                || (difference == 0.0 && layer.power_to_bath_w != 0.0)
            {
                return Err(DepositionError::Arithmetic.into());
            }
            conductances_w_k.insert(binding.convection.link, conductance);
            samples.push(GraphConjugateSurfaceSample {
                link: binding.convection.link,
                wall_cell: prepared.wall,
                fluid_cell: binding.convection.fluid_cell,
                wall_temperature_k: wall,
                fluid_temperature_k: bulk,
                geometry: binding.geometry,
                source_line: Arc::clone(binding.wall_layer.source_line()),
                layer,
                applied_heat_into_fluid_w: applied,
                representation_error_w: error,
                representation_allowance_w: representation,
            });
        }
        poll(cx).map_err(source_error)?;
        Ok(ResolvedGraphSurfaceConjugateHeat {
            admission,
            samples,
            conductances_w_k,
        })
    }
}

fn surface_error(
    link: u64,
    error: ConvectiveLayerError<GraphConvectionError>,
) -> GraphConjugateError {
    match error {
        ConvectiveLayerError::Thermal(error) => GraphConjugateError::Thermal(error),
        ConvectiveLayerError::Film(error) => GraphConjugateError::Convection(error),
        ConvectiveLayerError::SurfaceBalance {
            residual_w,
            iterations,
        } => GraphConjugateError::SurfaceBalance {
            link,
            residual_w,
            iterations,
        },
    }
}
