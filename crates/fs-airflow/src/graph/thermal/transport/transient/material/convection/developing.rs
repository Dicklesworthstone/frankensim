//! Section-integrated developing duct heat on original graph storage.
//!
//! Explicit entrance distances follow signed flow. Native fluid cells retain
//! mass, calorics, momentum geometry and time; native wall links or prescribed
//! reservoirs own one heat transfer. Local property updates are a held-property
//! approximation to the original constant-property Hausen mean relation.

use super::*;
use crate::graph::thermal::medium::convection::deposition::ConvectionSourceError;
use fs_conduction::transient::deposition::backward_euler::BoundaryConductances;
use fs_convection::developing::{CircularDuctSection, evaluate_circular_duct_section};
use fs_convection::ThermalConductivity;
use fs_material::fluid_transport::ThermalFluidTransportState;
use std::collections::{BTreeMap, BTreeSet};

/// Existing owner of a whole section's circular wetted-wall patch.
#[derive(Clone, Copy, Debug)]
pub enum DevelopingWall {
    Link(u64),
    Reservoir,
}

/// One native pipe cell within an authored uniform circular heated duct.
///
/// Entrance distances include upstream heated length, not just this cell.
/// `forward_start` measures from the entrance for positive branch flow;
/// `reverse_start` measures from the opposite entrance for negative flow.
/// A mixing junction, change of diameter, bend-induced secondary flow or changed
/// wall condition needs its own physical admission; this card cannot infer one.
#[derive(Clone, Copy, Debug)]
pub struct GraphDevelopingSection {
    pub fluid_cell: CellId,
    pub wall: DevelopingWall,
    pub diameter: Length,
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
    pub section: CircularDuctSection,
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
    mut resolve: impl FnMut(Temperature) -> Result<M, SourceGraphError>,
) -> Result<ResolvedDevelopingConvection<M>, GraphDevelopingConvectionError> {
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
        if !cells.insert(binding.fluid_cell) {
            return Err(invalid("a developing section must have one whole-wall convection owner"));
        }
        let Some(Some((branch, section))) = core.locations.get(&binding.fluid_cell) else {
            return Err(invalid("developing convection requires an original bound pipe cell"));
        };
        let diameter = binding.diameter.value();
        let forward = binding.forward_start.value();
        let reverse = binding.reverse_start.value();
        let length = core.section_lengths[*branch][*section];
        let flow_area = core.section_areas[*branch][*section];
        let circular_area = (core::f64::consts::PI / 4.0) * diameter * diameter;
        let wall_area = core::f64::consts::PI * diameter * length;
        if !positive(diameter) || !nonnegative(forward) || !nonnegative(reverse)
            || !near_positive(flow_area, circular_area) || !positive(wall_area)
            || !positive(forward + length) || !positive(reverse + length)
            || forward + length <= forward || reverse + length <= reverse
        {
            return Err(invalid("finite entrance distances and original uniform circular section geometry required"));
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
            return Err(invalid("developing wall must cover the original section's full circular perimeter and length"));
        }
        geometry.push(Geometry {
            binding, branch: *branch, flow_area, length, wall_temperature,
            fluid_temperature: thermal.temperature_k(binding.fluid_cell)
                .map_err(|error| native(error.into()))?,
            wall_area: area,
        });
    }
    geometry.sort_by_key(|g| g.binding.fluid_cell);
    let mut links = BTreeMap::new();
    let mut reservoirs = BTreeMap::new();
    let mut samples = Vec::with_capacity(geometry.len());
    for g in geometry {
        poll(cx).map_err(source)?;
        let binding = g.binding;
        let property_temperature_k = match binding.property_temperature {
            ConvectionPropertyTemperature::Reservoir => g.fluid_temperature,
            ConvectionPropertyTemperature::FilmMean => {
                0.5 * g.wall_temperature + 0.5 * g.fluid_temperature
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
        let evaluate_flow = |rate: f64| {
            let start = if rate < 0.0 { binding.reverse_start.value() }
                else { binding.forward_start.value() };
            let reynolds = transport.reynolds(rate.abs() / g.flow_area, binding.diameter.value())
                .map_err(|error| native(ConvectionSourceError::Material(error.into())))?;
            evaluate_circular_duct_section(
                reynolds, transport.prandtl(),
                start / binding.diameter.value(),
                (start + g.length) / binding.diameter.value(),
            ).map_err(|error| native(ConvectionSourceError::Material(error.into())))
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
        let coefficient = section.heat_transfer_coefficient(
            ThermalConductivity::new(transport.thermal_conductivity_w_m_k()), binding.diameter,
        ).map_err(|error| native(ConvectionSourceError::Material(error.into())))?;
        let conductance_w_k = g.wall_area * coefficient.value();
        let difference = g.wall_temperature - g.fluid_temperature;
        let heat_into_fluid_w = conductance_w_k * difference;
        if !positive(conductance_w_k) || !difference.is_finite() || !heat_into_fluid_w.is_finite()
            || (difference != 0.0 && heat_into_fluid_w == 0.0)
        {
            return Err(invalid("developing wall transfer is unrepresentable"));
        }
        match binding.wall {
            DevelopingWall::Link(id) => { links.insert(id, conductance_w_k); }
            DevelopingWall::Reservoir => { reservoirs.insert(binding.fluid_cell, conductance_w_k); }
        }
        samples.push(DevelopingConvectionSample {
            fluid_cell: binding.fluid_cell, wall: binding.wall,
            wall_temperature_k: g.wall_temperature, fluid_temperature_k: g.fluid_temperature,
            property_temperature_k, length_m: g.length, area_m2: g.wall_area,
            conductance_w_k, heat_into_fluid_w, material, thermal: snapshot, section,
        });
    }
    poll(cx).map_err(source)?;
    Ok(ResolvedDevelopingConvection { links, reservoirs, samples })
}
