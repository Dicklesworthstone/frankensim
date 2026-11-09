//! Conducting walls in series with the original developing-duct sections.
//!
//! This is a locally held, quasi-steady normal path. The original section card
//! retains its isothermal-wall and upstream entrance assumptions; no axial
//! wall spreading or nonisothermal boundary-layer history is inferred.

use super::super::conjugate::{GraphConjugateError, PreparedConjugateWall, apply_conductances};
use super::*;
use fs_matdb::PropertyLineSlice;
use fs_material::thermal::conduction::{
    ConvectiveLayerError, ConvectiveLayerResponse, LayerConductionConfig, ReservoirLayerGeometry,
    source::ResolvedConductionLayer,
};
use std::sync::Arc;

/// One whole developing section and its original finite wall's normal path.
/// `section.wall` must name a native `Link`; a prescribed reservoir has no
/// finite solid caloric trace. The same source, pure phase and exact contacted
/// area are retained by the common wall preparation.
pub struct GraphDevelopingConjugateLinkBinding<'a, 'm> {
    pub section: GraphDevelopingSection,
    pub wall_layer: &'a ResolvedConductionLayer<'m>,
    pub geometry: ReservoirLayerGeometry,
}

struct PreparedSection<'a, 'm> {
    binding: GraphDevelopingConjugateLinkBinding<'a, 'm>,
    link: u64,
    wall: PreparedConjugateWall,
}

/// Immutable source/contact identities; neither film nor solid k is frozen.
/// Native liquid/gas owners admit the original duct geometry and full signed
/// flow envelope at every proposed endpoint. Geometry is resolved once before
/// surface trials, and each trial uses the SAME section evaluator as the
/// ordinary developing path. An explicitly selected liquid wall correction
/// therefore queries the same fluid at every actual surface trial; its bulk
/// properties and native storage remain at the original fluid-cell state.
pub struct PreparedGraphDevelopingConjugateLinks<'a, 'm> {
    links: Vec<PreparedSection<'a, 'm>>,
    sections: Vec<GraphDevelopingSection>,
    maximum_sections: usize,
}

/// Actual surface/source evidence and one reciprocal finite-wall heat owner.
#[derive(Debug)]
pub struct GraphDevelopingConjugateSample<M> {
    pub link: u64,
    pub wall_cell: CellId,
    pub fluid_cell: CellId,
    pub wall_temperature_k: f64,
    pub fluid_temperature_k: f64,
    pub geometry: ReservoirLayerGeometry,
    pub source_line: Arc<PropertyLineSlice>,
    /// The retained section film uses the solved surface, not the stored wall
    /// temperature above. Its bare hA is diagnostic, never an extra heat path.
    pub layer: ConvectiveLayerResponse<DevelopingConvectionSample<M>>,
    pub applied_heat_into_fluid_w: f64,
    pub representation_error_w: f64,
    pub representation_allowance_w: f64,
}

#[derive(Debug)]
pub struct ResolvedGraphDevelopingConjugateHeat<M> {
    pub admission: GraphSourceAdmission,
    samples: Vec<GraphDevelopingConjugateSample<M>>,
    conductances_w_k: BTreeMap<u64, f64>,
}
impl<M> ResolvedGraphDevelopingConjugateHeat<M> {
    pub fn samples(&self) -> &[GraphDevelopingConjugateSample<M>] {
        &self.samples
    }
    pub fn conductances_w_k(&self) -> &BTreeMap<u64, f64> {
        &self.conductances_w_k
    }

    /// Publish only the combined reciprocal coefficients. Conflicting owners
    /// or cancellation refuse before changing any destination entry.
    pub fn apply_to(
        &self,
        destination: &mut BTreeMap<u64, f64>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), GraphConjugateError> {
        apply_conductances(&self.conductances_w_k, destination, checkpoint)
    }
}

impl<'a, 'm> PreparedGraphDevelopingConjugateLinks<'a, 'm> {
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        mut bindings: Vec<GraphDevelopingConjugateLinkBinding<'a, 'm>>,
        maximum_sections: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, GraphConjugateError> {
        checkpoint()?;
        if !(1..=1_000_000).contains(&maximum_sections)
            || bindings.len() > maximum_sections
            || bindings.len() > thermal.links().len()
        {
            return Err(GraphConjugateError::Invalid(
                "bounded native developing conjugate link count required",
            ));
        }
        bindings.sort_unstable_by_key(|b| b.section.fluid_cell);
        let mut links: Vec<PreparedSection<'a, 'm>> = Vec::with_capacity(bindings.len());
        let mut sections = Vec::with_capacity(bindings.len());
        let mut ids = BTreeSet::new();
        for binding in bindings {
            checkpoint()?;
            let DevelopingWall::Link(link) = binding.section.wall else {
                return Err(GraphConjugateError::Invalid(
                    "developing conjugacy requires an original finite wall link, not a reservoir",
                ));
            };
            if !ids.insert(link)
                || links
                    .last()
                    .is_some_and(|p| p.binding.section.fluid_cell == binding.section.fluid_cell)
            {
                return Err(GraphConjugateError::Invalid(
                    "a developing section and its whole wall need one conjugate owner",
                ));
            }
            let wall = PreparedConjugateWall::new(
                thermal,
                link,
                binding.section.fluid_cell,
                binding.geometry,
                binding.wall_layer,
            )?;
            sections.push(binding.section);
            links.push(PreparedSection {
                binding,
                link,
                wall,
            });
        }
        checkpoint()?;
        Ok(Self {
            links,
            sections,
            maximum_sections,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::graph::thermal::transport::transient::material) fn resolve<
        M: ThermalFluidTransport,
    >(
        &self,
        core: &BoundGraphCore<'_>,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        flow: GraphFlow<'_>,
        admission: GraphSourceAdmission,
        config: LayerConductionConfig,
        fluid: DevelopingFluid,
        mut resolve: impl FnMut(Temperature) -> Result<M, SourceGraphError>,
    ) -> Result<ResolvedGraphDevelopingConjugateHeat<M>, GraphConjugateError> {
        let source =
            |error| GraphConjugateError::Developing(GraphDevelopingConvectionError::Source(error));
        poll(cx).map_err(source)?;
        // Every original wall and full duct patch is admitted before the first
        // fluid source call. Stored and geometry samples share fluid-cell order.
        for prepared in &self.links {
            poll(cx).map_err(source)?;
            prepared
                .wall
                .temperatures(thermal, prepared.binding.wall_layer)?;
        }
        let geometry = admit_sections(core, cx, thermal, &self.sections, self.maximum_sections, fluid)
            .map_err(GraphConjugateError::Developing)?;
        let mut samples = Vec::with_capacity(self.links.len());
        let mut conductances_w_k = BTreeMap::new();
        for (prepared, g) in self.links.iter().zip(geometry) {
            poll(cx).map_err(source)?;
            let binding = &prepared.binding;
            let [wall, bulk] = prepared.wall.temperatures(thermal, binding.wall_layer)?;
            let layer = binding
                .geometry
                .evaluate_convective_material_to_bath(
                    binding.wall_layer.conductivity_curve(),
                    [wall, bulk],
                    config,
                    |surface| {
                        let film = evaluate_surface(core, cx, flow, &g, surface, &mut resolve)?;
                        Ok((film.heat_into_fluid_w, film.conductance_w_k, film))
                    },
                    || cx.checkpoint().map_err(|_| DepositionError::Cancelled),
                )
                .map_err(|error| match error {
                    ConvectiveLayerError::Thermal(error) => GraphConjugateError::Thermal(error),
                    ConvectiveLayerError::Film(error) => GraphConjugateError::Developing(error),
                    ConvectiveLayerError::SurfaceBalance {
                        residual_w,
                        iterations,
                    } => GraphConjugateError::SurfaceBalance {
                        link: prepared.link,
                        residual_w,
                        iterations,
                    },
                })?;
            let difference = wall - bulk;
            let conductance = layer.conductance_w_k;
            let applied = conductance * difference;
            let error = applied - layer.power_to_bath_w;
            let representation =
                8.0 * f64::EPSILON * applied.abs().max(layer.power_to_bath_w.abs());
            if !positive(conductance)
                || !applied.is_finite()
                || !error.is_finite()
                || error.abs() > representation
                || (difference != 0.0 && applied == 0.0)
                || (difference == 0.0 && layer.power_to_bath_w != 0.0)
            {
                return Err(DepositionError::Arithmetic.into());
            }
            conductances_w_k.insert(prepared.link, conductance);
            samples.push(GraphDevelopingConjugateSample {
                link: prepared.link,
                wall_cell: prepared.wall.cell,
                fluid_cell: binding.section.fluid_cell,
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
        poll(cx).map_err(source)?;
        Ok(ResolvedGraphDevelopingConjugateHeat {
            admission,
            samples,
            conductances_w_k,
        })
    }
}
