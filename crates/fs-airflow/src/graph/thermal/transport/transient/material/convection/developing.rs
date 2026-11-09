//! Section-integrated developing duct heat on original graph storage.
//!
//! Explicit entrance distances follow signed flow. Native fluid cells retain
//! mass, calorics, momentum geometry and time; native wall links or prescribed
//! reservoirs own one heat transfer. Local property updates are a held-property
//! approximation to the original constant-property mean duct relation.

use super::*;
use crate::graph::thermal::medium::convection::deposition::{
    ConvectionSourceError,
    source::{SurfaceConvectionMaterialSample, material_at},
};
use fs_conduction::transient::deposition::backward_euler::BoundaryConductances;
use fs_convection::developing::{
    CircularDuctSection, RectangularDuctSection, TurbulentCircularDuctSection,
    TurbulentTubePropertyCorrection,
    evaluate_circular_duct_section, evaluate_rectangular_duct_section,
    evaluate_turbulent_circular_duct_section,
};
use fs_convection::{CorrelationError, HeatTransferCoefficient, NusseltEvaluation, ThermalConductivity};
use fs_evidence::Evidence;
use fs_material::fluid_transport::ThermalFluidTransportState;
use std::collections::{BTreeMap, BTreeSet};

/// Source-integrated solid walls on the original finite developing sections.
pub mod conjugate;

/// Existing owner of a whole section's wetted-wall patch.
#[derive(Clone, Copy, Debug)]
pub enum DevelopingWall {
    Link(u64),
    Reservoir,
}

/// Physical cross-section of the original uniform heated passage.
#[derive(Clone, Copy, Debug)]
pub enum DevelopingDuct {
    /// Hydrodynamically developed circular inlet, using the Hausen relation.
    Circular {
        /// Physical inner diameter.
        diameter: Length,
    },
    /// Smooth turbulent circular passage using the corrected Gnielinski
    /// finite-heated-length relation. Bulk properties and K=1 are explicit:
    /// select Reservoir, with no film substitution or wall-viscosity correction.
    /// The original card must admit each positive entrance distance and the
    /// full signed native flow envelope. Roughness, transition, secondary flow
    /// and strong property variation require a different physical model.
    TurbulentCircular {
        diameter: Length,
    },
    /// The same developing turbulent circular section, with the explicit
    /// Gnielinski liquid correction (Pr_bulk/Pr_wall)^0.11. Only the native
    /// liquid graph owner admits this selection. Reservoir bulk properties
    /// still supply Re, Pr and conductivity; the SAME source supplies wall
    /// Pr at the actual fluid-facing surface, including a conjugate solution.
    /// A held section correction does not resolve an axial property profile.
    TurbulentCircularLiquidWall {
        diameter: Length,
    },
    /// Simultaneously developing rectangular inlet, using the original
    /// Shah-London CWT table slice. Its aspect ratio must be 0.5 and the
    /// actual source Prandtl number must be 0.72; neither is adjusted to fit.
    /// All four walls share the declared isothermal wall condition.
    Rectangular {
        /// One full internal side length.
        width: Length,
        /// The perpendicular full internal side length.
        height: Length,
    },
}

impl DevelopingDuct {
    fn metrics(self) -> Result<DuctMetrics, GraphDevelopingConvectionError> {
        let (diameter, area, perimeter, aspect) = match self {
            Self::Circular { diameter } | Self::TurbulentCircular { diameter }
                | Self::TurbulentCircularLiquidWall { diameter } => {
                let d = diameter.value();
                (
                    d,
                    (core::f64::consts::PI / 4.0) * d * d,
                    core::f64::consts::PI * d,
                    None,
                )
            }
            Self::Rectangular { width, height } => {
                let (w, h) = (width.value(), height.value());
                if !positive(w) || !positive(h) {
                    return Err(invalid(
                        "positive finite rectangular passage dimensions required",
                    ));
                }
                let area = w * h;
                let perimeter = 2.0 * w + 2.0 * h;
                (
                    4.0 * (area / perimeter),
                    area,
                    perimeter,
                    Some(w.min(h) / w.max(h)),
                )
            }
        };
        if !positive(diameter) || !positive(area) || !positive(perimeter) {
            return Err(invalid("developing duct cross-section is unrepresentable"));
        }
        Ok(DuctMetrics {
            diameter: Length::new(diameter),
            area,
            perimeter,
            aspect,
        })
    }
}

// The original source owner supplies this capability, never a material name
// or a gas-mixture heuristic. It cannot be selected by a duct binding.
#[derive(Clone, Copy)]
pub(in super::super) enum DevelopingFluid {
    Liquid,
    Gas,
}

#[derive(Clone, Copy)]
struct DuctMetrics {
    diameter: Length,
    area: f64,
    perimeter: f64,
    aspect: Option<f64>,
}

/// Original cumulative evaluations and derived transfer for the selected duct.
/// Rectangular sections retain the table/engineering-bridge distinction and
/// their unbounded local model discrepancy; circular sections retain Hausen's
/// original mean-card support without inventing a local validation band.
#[derive(Debug)]
pub enum DevelopingSectionEvaluation {
    /// Circular Hausen mean and its conservative finite section.
    Circular(CircularDuctSection),
    /// Original Gnielinski K=1 means and their actual section, including any
    /// explicitly selected liquid wall-property correction on that section.
    TurbulentCircular(TurbulentCircularDuctSection),
    /// Original rectangular table means and their conservative finite section.
    Rectangular(RectangularDuctSection),
}

impl DevelopingSectionEvaluation {
    /// Section-average Nusselt number based on hydraulic diameter.
    pub fn nusselt(&self) -> f64 {
        match self {
            Self::Circular(s) => s.nusselt(),
            Self::TurbulentCircular(s) => s.nusselt(),
            Self::Rectangular(s) => s.nusselt(),
        }
    }
    /// Original cumulative mean at the section start, absent at the entrance.
    pub fn upstream(&self) -> Option<&NusseltEvaluation> {
        match self {
            Self::Circular(s) => s.upstream(),
            Self::TurbulentCircular(s) => s.upstream(),
            Self::Rectangular(s) => s.upstream(),
        }
    }
    /// Original cumulative mean at the section end, retaining its source domain.
    pub fn downstream(&self) -> &NusseltEvaluation {
        match self {
            Self::Circular(s) => s.downstream(),
            Self::TurbulentCircular(s) => s.downstream(),
            Self::Rectangular(s) => s.downstream(),
        }
    }
    fn coefficient(&self, conductivity: f64, diameter: Length)
        -> Result<Evidence<HeatTransferCoefficient>, CorrelationError>
    {
        let conductivity = ThermalConductivity::new(conductivity);
        match self {
            Self::Circular(s) => s.heat_transfer_coefficient(conductivity, diameter),
            Self::TurbulentCircular(s) if s.property_correction().is_some() => {
                s.heat_transfer_coefficient_with_sensitivities(conductivity, diameter)
            }
            Self::TurbulentCircular(s) => s.heat_transfer_coefficient(conductivity, diameter),
            Self::Rectangular(s) => s.heat_transfer_coefficient(conductivity, diameter),
        }
    }
}

/// One native pipe cell within an authored uniform heated duct.
///
/// Entrance distances include upstream heated length, not just this cell.
/// `forward_start` measures from the entrance for positive branch flow;
/// `reverse_start` measures from the opposite entrance for negative flow.
/// A mixing junction, change of cross-section, bend-induced secondary flow or changed
/// wall condition needs its own physical admission; this card cannot infer one.
#[derive(Clone, Copy, Debug)]
pub struct GraphDevelopingSection {
    pub fluid_cell: CellId,
    pub wall: DevelopingWall,
    pub duct: DevelopingDuct,
    pub forward_start: Length,
    pub reverse_start: Length,
    pub property_temperature: ConvectionPropertyTemperature,
}

#[derive(Debug)]
pub enum GraphDevelopingConvectionError {
    Source(SourceGraphError),
    Convection(SourceCellConvectionError<SourceGraphError>),
}
impl core::fmt::Display for GraphDevelopingConvectionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "developing graph convection: {self:?}")
    }
}
impl std::error::Error for GraphDevelopingConvectionError {}

#[derive(Debug)]
pub struct GraphDevelopingHeat<M> {
    pub admission: GraphSourceAdmission,
    pub convection: ResolvedDevelopingConvection<M>,
}

#[derive(Debug)]
pub struct DevelopingConvectionSample<M> {
    pub fluid_cell: CellId,
    pub wall: DevelopingWall,
    /// Actual fluid-facing temperature; a conjugate wall retains its separate
    /// stored cell trace in the surrounding surface-balance report.
    pub wall_temperature_k: f64,
    pub fluid_temperature_k: f64,
    pub property_temperature_k: f64,
    pub length_m: f64,
    pub area_m2: f64,
    pub conductance_w_k: f64,
    /// Diagnostic only: native links/boundaries apply this heat exactly once.
    pub heat_into_fluid_w: f64,
    pub material: M,
    pub thermal: ThermalFluidTransportState,
    /// The original liquid queried at the ACTUAL fluid-facing wall when the
    /// explicit correction is selected. This diagnostic source has no mass,
    /// caloric state or density admission as a second stored fluid cell.
    pub surface: Option<SurfaceConvectionMaterialSample<M>>,
    /// Actual section coefficient and model evidence. A corrected turbulent
    /// section includes its original analytic input partials, including
    /// Pr_wall; these are not total derivatives through the material source.
    pub coefficient: Evidence<HeatTransferCoefficient>,
    pub section: DevelopingSectionEvaluation,
}

#[derive(Debug)]
pub struct ResolvedDevelopingConvection<M> {
    links: BTreeMap<u64, f64>,
    reservoirs: BTreeMap<CellId, f64>,
    samples: Vec<DevelopingConvectionSample<M>>,
}
impl<M> ResolvedDevelopingConvection<M> {
    pub fn link_conductances_w_k(&self) -> &BTreeMap<u64, f64> { &self.links }
    pub fn reservoir_conductances_w_k(&self) -> &BTreeMap<CellId, f64> { &self.reservoirs }
    pub fn samples(&self) -> &[DevelopingConvectionSample<M>] { &self.samples }

    /// Apply to a single native thermal trial, checking both destinations before
    /// mutation. Never also add the diagnostic heat powers as cell sources.
    pub fn apply_to(
        &self,
        links: &mut BTreeMap<u64, f64>,
        boundaries: &mut BoundaryConductances,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), ConvectionSourceError> {
        for id in self.links.keys() {
            if cancelled() { return Err(ConvectionSourceError::Cancelled); }
            if links.contains_key(id) {
                return Err(ConvectionSourceError::Invalid("developing wall link already has an owner"));
            }
        }
        for id in self.reservoirs.keys() {
            if cancelled() { return Err(ConvectionSourceError::Cancelled); }
            if boundaries.reservoir_w_k.contains_key(id) {
                return Err(ConvectionSourceError::Invalid("developing wall reservoir already has an owner"));
            }
        }
        if cancelled() { return Err(ConvectionSourceError::Cancelled); }
        links.extend(self.links.iter().map(|(&id, &g)| (id, g)));
        boundaries.reservoir_w_k.extend(self.reservoirs.iter().map(|(&id, &g)| (id, g)));
        Ok(())
    }
}

fn native(error: ConvectionSourceError) -> GraphDevelopingConvectionError {
    GraphDevelopingConvectionError::Convection(SourceCellConvectionError::Convection(error))
}
fn invalid(reason: &'static str) -> GraphDevelopingConvectionError {
    native(ConvectionSourceError::Invalid(reason))
}

struct Geometry<'a> {
    binding: &'a GraphDevelopingSection,
    duct: DuctMetrics,
    branch: usize,
    flow_area: f64,
    length: f64,
    wall_temperature: f64,
    fluid_temperature: f64,
    wall_area: f64,
}

#[allow(clippy::too_many_arguments)]
pub(in super::super) fn evaluate<M: ThermalFluidTransport>(
    core: &BoundGraphCore<'_>,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'_>,
    flow: GraphFlow<'_>,
    bindings: &[GraphDevelopingSection],
    maximum_sections: usize,
    fluid: DevelopingFluid,
    mut resolve: impl FnMut(Temperature) -> Result<M, SourceGraphError>,
) -> Result<ResolvedDevelopingConvection<M>, GraphDevelopingConvectionError> {
    let geometry = admit_sections(core, cx, thermal, bindings, maximum_sections, fluid)?;
    let mut links = BTreeMap::new();
    let mut reservoirs = BTreeMap::new();
    let mut samples = Vec::with_capacity(geometry.len());
    for g in geometry {
        let sample = evaluate_surface(core, cx, flow, &g, g.wall_temperature, &mut resolve)?;
        match g.binding.wall {
            DevelopingWall::Link(id) => { links.insert(id, sample.conductance_w_k); }
            DevelopingWall::Reservoir => {
                reservoirs.insert(g.binding.fluid_cell, sample.conductance_w_k);
            }
        }
        samples.push(sample);
    }
    poll(cx).map_err(GraphDevelopingConvectionError::Source)?;
    Ok(ResolvedDevelopingConvection { links, reservoirs, samples })
}

// This admission is shared by ordinary and source-integrated finite walls.
// It runs once for the whole endpoint before any material query or root trial.
fn admit_sections<'a>(
    core: &BoundGraphCore<'_>,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'_>,
    bindings: &'a [GraphDevelopingSection],
    maximum_sections: usize,
    fluid: DevelopingFluid,
) -> Result<Vec<Geometry<'a>>, GraphDevelopingConvectionError> {
    let source = GraphDevelopingConvectionError::Source;
    poll(cx).map_err(source)?;
    if !(1..=1_000_000).contains(&maximum_sections)
        || bindings.len() > maximum_sections || bindings.len() > thermal.cells().len()
    {
        return Err(invalid("bounded native developing section count required"));
    }
    let mut cells = BTreeSet::new();
    let mut link_ids = BTreeSet::new();
    let mut geometry = Vec::with_capacity(bindings.len());
    // Admit every whole-wall patch before calling any external material source.
    for binding in bindings {
        poll(cx).map_err(source)?;
        if matches!(binding.duct, DevelopingDuct::TurbulentCircularLiquidWall { .. })
            && !matches!(fluid, DevelopingFluid::Liquid)
        {
            return Err(invalid("liquid wall-property correction requires the original liquid graph owner"));
        }
        if matches!(binding.duct, DevelopingDuct::TurbulentCircular { .. }
            | DevelopingDuct::TurbulentCircularLiquidWall { .. })
            && binding.property_temperature != ConvectionPropertyTemperature::Reservoir
        {
            return Err(invalid(
                "developing turbulent tubes require actual bulk properties with their explicit wall-property selection",
            ));
        }
        if !cells.insert(binding.fluid_cell) {
            return Err(invalid("a developing section must have one whole-wall convection owner"));
        }
        let Some(Some((branch, section))) = core.locations.get(&binding.fluid_cell) else {
            return Err(invalid("developing convection requires an original bound pipe cell"));
        };
        let duct = binding.duct.metrics()?;
        let forward = binding.forward_start.value();
        let reverse = binding.reverse_start.value();
        let length = core.section_lengths[*branch][*section];
        let flow_area = core.section_areas[*branch][*section];
        let wall_area = duct.perimeter * length;
        if !nonnegative(forward) || !nonnegative(reverse)
            || !near_positive(flow_area, duct.area) || !positive(wall_area)
            || !positive(forward + length) || !positive(reverse + length)
            || forward + length <= forward || reverse + length <= reverse
        {
            return Err(invalid("finite entrance distances and original uniform duct section geometry required"));
        }
        let cell = thermal.cells().get(&binding.fluid_cell)
            .ok_or_else(|| invalid("missing original developing fluid cell"))?;
        let (area, wall_temperature) = match binding.wall {
            DevelopingWall::Link(id) => {
                if !link_ids.insert(id) {
                    return Err(invalid("duplicate developing wall link"));
                }
                let link = thermal.links().get(&id)
                    .ok_or_else(|| invalid("missing original developing wall link"))?;
                if link.conductance_w_k != 0.0 {
                    return Err(invalid("developing wall link requires zero authored conductance"));
                }
                let wall = if link.a == binding.fluid_cell { link.b }
                    else if link.b == binding.fluid_cell { link.a }
                    else { return Err(invalid("developing fluid must be an original wall-link endpoint")); };
                (link.area_m2, thermal.temperature_k(wall)
                    .map_err(|error| native(error.into()))?)
            }
            DevelopingWall::Reservoir => {
                let reservoir = cell.boundary.reservoir
                    .ok_or_else(|| invalid("developing wall requires an existing native reservoir"))?;
                if reservoir.conductance_w_k != 0.0 {
                    return Err(invalid("developing reservoir requires zero authored conductance"));
                }
                (reservoir.area_m2, reservoir.temperature_k)
            }
        };
        if !near_positive(area, wall_area) {
            return Err(invalid("developing wall must cover the original section's full wetted perimeter and length"));
        }
        geometry.push(Geometry {
            binding, duct, branch: *branch, flow_area, length, wall_temperature,
            fluid_temperature: thermal.temperature_k(binding.fluid_cell)
                .map_err(|error| native(error.into()))?,
            wall_area: area,
        });
    }
    geometry.sort_by_key(|g| g.binding.fluid_cell);
    poll(cx).map_err(source)?;
    Ok(geometry)
}

// Evaluate the SAME section/source and signed native flow envelope at an
// actual passive surface, returning evidence only, without an applicable map.
fn evaluate_surface<M: ThermalFluidTransport>(
    core: &BoundGraphCore<'_>,
    cx: &Cx<'_>,
    flow: GraphFlow<'_>,
    g: &Geometry<'_>,
    surface_temperature: f64,
    resolve: &mut impl FnMut(Temperature) -> Result<M, SourceGraphError>,
) -> Result<DevelopingConvectionSample<M>, GraphDevelopingConvectionError> {
    let source = GraphDevelopingConvectionError::Source;
    poll(cx).map_err(source)?;
    if !positive(surface_temperature)
        || surface_temperature < g.wall_temperature.min(g.fluid_temperature)
        || surface_temperature > g.wall_temperature.max(g.fluid_temperature)
    {
        return Err(invalid("developing surface must lie between its actual wall and fluid temperatures"));
    }
    let binding = g.binding;
    let property_temperature_k = match binding.property_temperature {
        ConvectionPropertyTemperature::Reservoir => g.fluid_temperature,
        ConvectionPropertyTemperature::FilmMean => {
            0.5 * surface_temperature + 0.5 * g.fluid_temperature
        }
    };
    let material = resolve(Temperature::new(property_temperature_k)).map_err(|error| {
        GraphDevelopingConvectionError::Convection(SourceCellConvectionError::Source {
            cell: binding.fluid_cell, property_temperature_k, source: error,
        })
    })?;
    let snapshot = material.thermal_transport_state()
        .map_err(|error| native(ConvectionSourceError::Material(error.into())))?;
    let transport = material.transport_state()
        .map_err(|error| native(ConvectionSourceError::Material(error.into())))?;
    if snapshot.temperature_k() != property_temperature_k || snapshot.transport() != transport {
        return Err(invalid("developing source must match requested temperature and transport state"));
    }
    core.admit_density(binding.fluid_cell, transport).map_err(source)?;
    let surface = if matches!(binding.duct, DevelopingDuct::TurbulentCircularLiquidWall { .. }) {
        poll(cx).map_err(source)?;
        Some(material_at(binding.fluid_cell, surface_temperature, || {
            resolve(Temperature::new(surface_temperature))
        }).map_err(GraphDevelopingConvectionError::Convection)?)
    } else { None };
    let correction = surface.as_ref().map(|wall| TurbulentTubePropertyCorrection::Liquid {
        wall_prandtl: wall.thermal.transport().prandtl(),
    });
    let evaluate_flow = |rate: f64| {
        let start = if rate < 0.0 { binding.reverse_start.value() }
            else { binding.forward_start.value() };
        let reynolds = transport.reynolds(rate.abs() / g.flow_area, g.duct.diameter.value())
            .map_err(|error| native(ConvectionSourceError::Material(error.into())))?;
        let start_ratio = start / g.duct.diameter.value();
        let end_ratio = (start + g.length) / g.duct.diameter.value();
        if start > 0.0 && start_ratio == 0.0 {
            return Err(invalid("positive developing distance cannot underflow to a new entrance"));
        }
        match g.duct.aspect {
            None if matches!(binding.duct, DevelopingDuct::TurbulentCircular { .. }
                | DevelopingDuct::TurbulentCircularLiquidWall { .. }) => {
                evaluate_turbulent_circular_duct_section(
                    reynolds, transport.prandtl(), start_ratio, end_ratio,
                ).and_then(|section| match correction {
                    Some(correction) => section.with_property_correction(correction),
                    None => Ok(section),
                }).map(DevelopingSectionEvaluation::TurbulentCircular)
            }
            None => evaluate_circular_duct_section(
                reynolds, transport.prandtl(), start_ratio, end_ratio,
            ).map(DevelopingSectionEvaluation::Circular),
            Some(aspect) => evaluate_rectangular_duct_section(
                reynolds, transport.prandtl(), aspect, start_ratio, end_ratio,
            ).map(DevelopingSectionEvaluation::Rectangular),
        }.map_err(|error| native(ConvectionSourceError::Material(error.into())))
    };
    let average = flow.average(g.branch);
    let section = evaluate_flow(average)?;
    // Every actually traversed endpoint retains its own entrance direction.
    // A signed reversal crosses Re=0, outside the original forced card.
    let initial = flow.initial(g.branch);
    let final_flow = flow.final_flow(g.branch);
    if initial != average { evaluate_flow(initial)?; }
    if final_flow != average && final_flow != initial { evaluate_flow(final_flow)?; }
    if (initial < 0.0 && final_flow > 0.0) || (initial > 0.0 && final_flow < 0.0) {
        evaluate_flow(0.0)?;
    }
    let coefficient = section.coefficient(
        transport.thermal_conductivity_w_m_k(), g.duct.diameter,
    ).map_err(|error| native(ConvectionSourceError::Material(error.into())))?;
    let conductance_w_k = g.wall_area * coefficient.value.value();
    let difference = surface_temperature - g.fluid_temperature;
    let heat_into_fluid_w = conductance_w_k * difference;
    if !positive(conductance_w_k) || !difference.is_finite() || !heat_into_fluid_w.is_finite()
        || (difference != 0.0 && heat_into_fluid_w == 0.0)
    {
        return Err(invalid("developing wall transfer is unrepresentable"));
    }
    poll(cx).map_err(source)?;
    Ok(DevelopingConvectionSample {
        fluid_cell: binding.fluid_cell, wall: binding.wall,
        wall_temperature_k: surface_temperature, fluid_temperature_k: g.fluid_temperature,
        property_temperature_k, length_m: g.length, area_m2: g.wall_area,
        conductance_w_k, heat_into_fluid_w, material, thermal: snapshot, surface,
        coefficient, section,
    })
}
