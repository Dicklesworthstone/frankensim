//! One solid trace, its actual surface, and the original finite fluid cell.
//!
//! Bulk-property forced ducts use the original fixed-film layer closure.
//! Surface-dependent laws use the same source and geometry through the shared
//! nonlinear material/film balance, without another cell or heat owner.

mod surface;
pub use surface::{GraphConjugateSurfaceSample, ResolvedGraphSurfaceConjugateHeat};

use super::{GraphConvectionError, GraphConvectionLaw, GraphConvectiveLink, GraphSourceHeat};
use crate::graph::thermal::medium::convection::deposition::source::{
    ConvectionPropertyTemperature,
    links::{ConvectiveLinkSample, ResolvedLinkConvection},
};
use crate::graph::thermal::transport::transient::material::GraphSourceAdmission;
use fs_convection::CorrelationId;
use fs_matdb::PropertyLineSlice;
use fs_material::thermal::conduction::{
    LayerConductionConfig, ReservoirLayerGeometry, ReservoirLayerResponse,
    source::{LayerSourceError, ResolvedConductionLayer},
};
use fs_thermal::enthalpy::{
    CellGeometry, CellId, ContactPatch, DepositionError, DepositionSimulation, MaterialId,
    ThermalLaw,
};
use std::collections::BTreeMap;
use std::sync::Arc;

/// An original zero-conductance contact and its homogeneous solid path.
/// The layer's caloric object must be the exact wall-cell law. Geometry spans
/// that cell's temperature trace to the fluid-facing surface; its actual face
/// area must equal the original contact patch. No extra solid storage or
/// lateral spreading is inferred. The caller admits this quasi-steady path.
pub struct GraphConjugateLinkBinding<'a, 'm> {
    pub convection: GraphConvectiveLink,
    pub wall_layer: &'a ResolvedConductionLayer<'m>,
    pub geometry: ReservoirLayerGeometry,
}

#[derive(Debug)]
pub enum GraphConjugateError {
    Invalid(&'static str),
    Convection(GraphConvectionError),
    Developing(super::developing::GraphDevelopingConvectionError),
    Layer(LayerSourceError),
    Thermal(DepositionError),
    SurfaceBalance {
        link: u64,
        residual_w: f64,
        iterations: usize,
    },
}
impl core::fmt::Display for GraphConjugateError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native conjugate convection: {self:?}")
    }
}
impl std::error::Error for GraphConjugateError {}
impl From<LayerSourceError> for GraphConjugateError {
    fn from(error: LayerSourceError) -> Self {
        Self::Layer(error)
    }
}
impl From<DepositionError> for GraphConjugateError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}

struct PreparedLink<'a, 'm> {
    binding: GraphConjugateLinkBinding<'a, 'm>,
    wall: PreparedConjugateWall,
}

/// One original contact and caloric trace, shared by all native film owners.
/// It contains no heat coefficient, source law or additional thermal storage.
pub(super) struct PreparedConjugateWall {
    link: u64,
    fluid: CellId,
    patch: ContactPatch,
    pub(super) cell: CellId,
    wall_geometry: CellGeometry,
    wall_material: MaterialId,
}

/// Prepared geometry and source identities, never a frozen coefficient.
/// `new`/`resolve` select the bulk-property forced-duct fast path. The explicit
/// `new_with_surface_feedback` preparation also admits original film-property,
/// sphere, plate, natural and mixed laws through the native source graph's
/// `evaluate_conjugate_links` and quasi-steady counterpart.
/// The original convection owner still admits every actual flow endpoint and
/// the interval average, all source limits, and the selected duct assumptions.
pub struct PreparedGraphConjugateLinks<'a, 'm> {
    links: Vec<PreparedLink<'a, 'm>>,
    convection: Vec<GraphConvectiveLink>,
}

/// Source and physical surface response for one accepted endpoint coefficient.
#[derive(Debug)]
pub struct GraphConjugateLinkSample {
    pub link: u64,
    pub wall_cell: CellId,
    pub fluid_cell: CellId,
    pub wall_temperature_k: f64,
    pub fluid_temperature_k: f64,
    /// Exact original layer response, including solved surface T and heat residual.
    /// Its diagnostic secant is not the coefficient applied to the native link.
    pub layer: ReservoirLayerResponse,
    pub geometry: ReservoirLayerGeometry,
    pub source_line: Arc<PropertyLineSlice>,
    pub conductance_w_k: f64,
    /// Actual G*(T_wall_cell - T_fluid), used once by the reciprocal native link.
    pub applied_heat_into_fluid_w: f64,
    /// Applied power minus the layer owner's accepted signed power.
    pub representation_error_w: f64,
    pub representation_allowance_w: f64,
}

/// One applicable heat owner. The bare hA map is retained privately only for
/// source evidence; it cannot be applied alongside the combined conductance.
#[derive(Debug)]
pub struct ResolvedGraphConjugateHeat<M> {
    pub admission: GraphSourceAdmission,
    fluid: ResolvedLinkConvection<M>,
    samples: Vec<GraphConjugateLinkSample>,
    conductances_w_k: BTreeMap<u64, f64>,
}
impl<M> ResolvedGraphConjugateHeat<M> {
    pub fn samples(&self) -> &[GraphConjugateLinkSample] {
        &self.samples
    }
    /// Original bulk-fluid sources and correlation evidence in the same link
    /// order. Their bare film powers are diagnostics, not another heat transfer.
    pub fn fluid_samples(&self) -> &[ConvectiveLinkSample<M>] {
        self.fluid.samples()
    }
    pub fn conductances_w_k(&self) -> &BTreeMap<u64, f64> {
        &self.conductances_w_k
    }

    /// Add exactly one reciprocal coefficient per original link. Any competing
    /// override refuses before modifying the destination, even if numerically equal.
    pub fn apply_to(
        &self,
        destination: &mut BTreeMap<u64, f64>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), GraphConjugateError> {
        apply_conductances(&self.conductances_w_k, destination, checkpoint)
    }
}

pub(super) fn apply_conductances(
    conductances: &BTreeMap<u64, f64>,
    destination: &mut BTreeMap<u64, f64>,
    mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<(), GraphConjugateError> {
    checkpoint()?;
    for link in conductances.keys() {
        checkpoint()?;
        if destination.contains_key(link) {
            return Err(GraphConjugateError::Invalid(
                "conjugate link already has a heat owner",
            ));
        }
    }
    checkpoint()?;
    destination.extend(conductances.iter().map(|(&id, &g)| (id, g)));
    Ok(())
}

impl<'a, 'm> PreparedGraphConjugateLinks<'a, 'm> {
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        bindings: Vec<GraphConjugateLinkBinding<'a, 'm>>,
        maximum_links: usize,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, GraphConjugateError> {
        Self::prepare(thermal, bindings, maximum_links, true, checkpoint)
    }

    /// Prepare exact solid/contact identities for a nonlinear source film.
    /// The original liquid/gas graph owner admits correlation, signed flow,
    /// property convention and source support at every actual trial surface.
    /// No new enclosed-flow or nonisothermal boundary-layer model is inferred.
    pub fn new_with_surface_feedback(
        thermal: &DepositionSimulation<'_>,
        bindings: Vec<GraphConjugateLinkBinding<'a, 'm>>,
        maximum_links: usize,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, GraphConjugateError> {
        Self::prepare(thermal, bindings, maximum_links, false, checkpoint)
    }

    fn prepare(
        thermal: &DepositionSimulation<'_>,
        mut bindings: Vec<GraphConjugateLinkBinding<'a, 'm>>,
        maximum_links: usize,
        bulk_only: bool,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, GraphConjugateError> {
        checkpoint()?;
        if !(1..=1_000_000).contains(&maximum_links)
            || bindings.len() > maximum_links
            || bindings.len() > thermal.links().len()
        {
            return Err(GraphConjugateError::Invalid(
                "bounded native conjugate link count required",
            ));
        }
        bindings.sort_unstable_by_key(|b| b.convection.link);
        let mut links: Vec<PreparedLink<'a, 'm>> = Vec::with_capacity(bindings.len());
        let mut convection = Vec::with_capacity(bindings.len());
        for binding in bindings {
            checkpoint()?;
            if bulk_only {
                admit_law(binding.convection)?;
            }
            if links
                .last()
                .is_some_and(|p| p.binding.convection.link == binding.convection.link)
            {
                return Err(GraphConjugateError::Invalid(
                    "duplicate native conjugate link",
                ));
            }
            let wall = PreparedConjugateWall::new(
                thermal, binding.convection.link, binding.convection.fluid_cell,
                binding.geometry, binding.wall_layer,
            )?;
            convection.push(binding.convection);
            links.push(PreparedLink { binding, wall });
        }
        checkpoint()?;
        Ok(Self { links, convection })
    }

    /// Resolve at the actual proposed native endpoint. The callback must invoke
    /// the bound liquid/gas graph owner's existing `evaluate_links` or
    /// `evaluate_quasi_steady_links` with the supplied endpoint and bindings
    /// and the SAME thermal/source/flow candidate. Its original endpoint-flow
    /// admission and source evidence are preserved; bare film heat stays private.
    /// No layer conductivity is sampled at a mean temperature or copied.
    pub fn resolve<M>(
        &self,
        thermal: &DepositionSimulation<'_>,
        config: LayerConductionConfig,
        evaluate_convection: impl FnOnce(
            &DepositionSimulation<'_>,
            &[GraphConvectiveLink],
        ) -> Result<GraphSourceHeat<M>, GraphConvectionError>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ResolvedGraphConjugateHeat<M>, GraphConjugateError> {
        checkpoint()?;
        // Refuse changed wall/contact ownership before any fluid-source work.
        for link in &self.links {
            checkpoint()?;
            admit_law(link.binding.convection)?;
            link.temperatures(thermal)?;
        }
        let GraphSourceHeat {
            admission,
            convection: fluid,
        } = evaluate_convection(thermal, &self.convection)
            .map_err(GraphConjugateError::Convection)?;
        if fluid.samples().len() != self.links.len() {
            return Err(GraphConjugateError::Invalid(
                "conjugate callback changed its link set",
            ));
        }
        let mut samples = Vec::with_capacity(self.links.len());
        let mut conductances_w_k = BTreeMap::new();
        for (prepared, sample) in self.links.iter().zip(fluid.samples()) {
            checkpoint()?;
            let binding = &prepared.binding;
            let [wall, bulk] = prepared.temperatures(thermal)?;
            let correlation = admit_law(binding.convection)?;
            if sample.link != binding.convection.link
                || sample.wall_cell != prepared.wall.cell
                || sample.fluid_cell != binding.convection.fluid_cell
                || sample.wall_temperature_k != wall
                || sample.fluid_temperature_k != bulk
                || sample.property_temperature_k != bulk
                || sample.area_m2 != prepared.wall.patch.link.area_m2
                || sample.surface.is_some()
                || sample
                    .evaluation
                    .as_ref()
                    .is_none_or(|e| e.nusselt.card().id != correlation || e.section.is_some())
            {
                return Err(GraphConjugateError::Invalid(
                    "conjugate callback returned foreign or stale film evidence",
                ));
            }
            let film_resistance = 1.0 / sample.conductance_w_k;
            if !positive(sample.conductance_w_k) || !positive(film_resistance) {
                return Err(DepositionError::Arithmetic.into());
            }
            let layer = binding.geometry.evaluate_material_to_bath(
                binding.wall_layer.conductivity_curve(),
                [wall, bulk],
                film_resistance,
                config,
                &mut checkpoint,
            )?;
            let difference = wall - bulk;
            let conductance = if difference == 0.0 {
                layer.conductance_w_k
            } else {
                layer.power_to_reservoir_w / difference
            };
            let applied = conductance * difference;
            let error = applied - layer.power_to_reservoir_w;
            let representation =
                8.0 * f64::EPSILON * applied.abs().max(layer.power_to_reservoir_w.abs());
            let tolerance = config.absolute_heat_tolerance_w
                + config.relative_heat_tolerance
                    * applied.abs().max(layer.power_to_reservoir_w.abs());
            if !positive(conductance)
                || conductance > sample.conductance_w_k
                || !applied.is_finite()
                || !error.is_finite()
                || !tolerance.is_finite()
                || error.abs() > tolerance + representation
                || (difference != 0.0 && applied == 0.0)
                || (difference == 0.0 && layer.power_to_reservoir_w != 0.0)
            {
                return Err(DepositionError::Arithmetic.into());
            }
            conductances_w_k.insert(sample.link, conductance);
            samples.push(GraphConjugateLinkSample {
                link: sample.link,
                wall_cell: prepared.wall.cell,
                fluid_cell: sample.fluid_cell,
                wall_temperature_k: wall,
                fluid_temperature_k: bulk,
                layer,
                geometry: binding.geometry,
                source_line: Arc::clone(binding.wall_layer.source_line()),
                conductance_w_k: conductance,
                applied_heat_into_fluid_w: applied,
                representation_error_w: error,
                representation_allowance_w: representation,
            });
        }
        checkpoint()?;
        Ok(ResolvedGraphConjugateHeat {
            admission,
            fluid,
            samples,
            conductances_w_k,
        })
    }
}

impl PreparedLink<'_, '_> {
    fn temperatures(
        &self,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<[f64; 2], GraphConjugateError> {
        self.wall.temperatures(thermal, self.binding.wall_layer)
    }
}

impl PreparedConjugateWall {
    pub(super) fn new(
        thermal: &DepositionSimulation<'_>,
        link: u64,
        fluid: CellId,
        geometry: ReservoirLayerGeometry,
        layer: &ResolvedConductionLayer<'_>,
    ) -> Result<Self, GraphConjugateError> {
        let patch = thermal.contact_patch(link).ok_or(
            GraphConjugateError::Invalid("inactive native conjugate link"),
        )?;
        if patch.link.conductance_w_k != 0.0
            || patch.link.area_m2 != geometry.surface_area_m2()
        {
            return Err(GraphConjugateError::Invalid(
                "conjugate geometry needs the exact original area and zero authored conductance",
            ));
        }
        let cell = if fluid == patch.link.a { patch.link.b }
            else if fluid == patch.link.b { patch.link.a }
            else { return Err(GraphConjugateError::Invalid(
                "conjugate fluid is not a link endpoint",
            )); };
        let original = &thermal.cells()[&cell];
        let prepared = Self {
            link, fluid, patch, cell,
            wall_geometry: original.geometry,
            wall_material: original.material_id,
        };
        prepared.temperatures(thermal, layer)?;
        Ok(prepared)
    }

    pub(super) fn temperatures(
        &self,
        thermal: &DepositionSimulation<'_>,
        layer: &ResolvedConductionLayer<'_>,
    ) -> Result<[f64; 2], GraphConjugateError> {
        if thermal.contact_patch(self.link) != Some(self.patch) {
            return Err(GraphConjugateError::Invalid(
                "conjugate contact changed; reprepare its geometry",
            ));
        }
        let cell = thermal
            .cells()
            .get(&self.cell)
            .ok_or(GraphConjugateError::Invalid("inactive conjugate wall"))?;
        if cell.geometry != self.wall_geometry
            || cell.material_id != self.wall_material
            || !thermal.uses_material_law(
                cell.material_id,
                layer.caloric_curve() as &dyn ThermalLaw,
            )
        {
            return Err(GraphConjugateError::Invalid(
                "conjugate wall needs its original geometry and caloric object",
            ));
        }
        Ok([
            layer.temperature_at_enthalpy(cell.specific_h_j_kg)?,
            thermal.temperature_k(self.fluid)?,
        ])
    }
}

fn admit_law(binding: GraphConvectiveLink) -> Result<CorrelationId, GraphConjugateError> {
    if matches!(
        binding.property_temperature,
        ConvectionPropertyTemperature::Reservoir
    ) && let GraphConvectionLaw::Forced { correlation, .. } = binding.law
        && matches!(
            correlation,
            CorrelationId::CircularDuctLaminarCwt
                | CorrelationId::CircularDuctLaminarChf
                | CorrelationId::RectangularDuctLaminarCwt
                | CorrelationId::RectangularDuctLaminarChf
                | CorrelationId::DittusBoelter
                | CorrelationId::Gnielinski
        )
    {
        return Ok(correlation);
    }
    Err(GraphConjugateError::Invalid(
        "conjugate solid layers require a supported bulk-property forced duct law",
    ))
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
