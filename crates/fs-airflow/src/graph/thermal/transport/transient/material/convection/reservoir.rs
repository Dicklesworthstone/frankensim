//! Material convection between an original graph fluid cell and prescribed wall.
//!
//! The native reservoir retains its wall temperature and covered reference area.
//! This supplies only its endpoint conductance. Native thermal storage, heat and
//! entropy accounting, scheduled boundary changes and coupled acceptance remain
//! authoritative. Never add the diagnostic wall-to-fluid power as another source.

use super::*;
use super::laws::{FlowConvectionError, FlowConvectionPoint};
use crate::graph::thermal::medium::convection::deposition::{
    CellConvectionEvaluation, ConvectionSourceError,
};
use fs_conduction::transient::deposition::backward_euler::BoundaryConductances;
use fs_material::fluid_transport::ThermalFluidTransportState;
use std::collections::{BTreeMap, BTreeSet};

/// Convection on one already-present prescribed wall reservoir. Its authored
/// conductance must be zero; its area and temperature are never synthesized.
#[derive(Clone, Copy, Debug)]
pub struct GraphConvectiveReservoir {
    pub fluid_cell: CellId,
    pub characteristic_length: Length,
    /// Reservoir means the actual finite fluid's BULK temperature here;
    /// FilmMean uses the mean of that temperature and the prescribed wall.
    pub property_temperature: ConvectionPropertyTemperature,
    pub law: GraphConvectionLaw,
}

#[derive(Debug)]
pub enum GraphReservoirConvectionError {
    Source(SourceGraphError),
    Convection(SourceCellConvectionError<SourceGraphError>),
    Geometry { cell: CellId, source: MaterialConvectionError },
    Envelope { cell: CellId, source: MaterialConvectionError },
}
impl core::fmt::Display for GraphReservoirConvectionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite graph reservoir convection: {self:?}")
    }
}
impl std::error::Error for GraphReservoirConvectionError {}

/// Same source/flow admission as finite-wall convection, and the native
/// prescribed-boundary conductances from that proposed thermal endpoint.
#[derive(Debug)]
pub struct GraphReservoirHeat<M> {
    pub admission: GraphSourceAdmission,
    pub convection: ResolvedReservoirConvection<M>,
}

#[derive(Debug)]
pub struct ReservoirConvectionSample<M> {
    pub fluid_cell: CellId,
    pub wall_temperature_k: f64,
    pub fluid_temperature_k: f64,
    pub property_temperature_k: f64,
    pub area_m2: f64,
    pub conductance_w_k: f64,
    /// Diagnostic power, positive from the prescribed wall into the fluid.
    /// Native reservoir accounting supplies the opposite external debit.
    pub heat_into_fluid_w: f64,
    pub material: M,
    pub thermal: ThermalFluidTransportState,
    /// None only at exact natural-convection thermal equilibrium.
    pub evaluation: Option<CellConvectionEvaluation>,
}

/// One bounded set of original prescribed-wall conductances and material owners.
/// These are total W/K coefficients, not additional heat-source watts.
#[derive(Debug)]
pub struct ResolvedReservoirConvection<M> {
    reservoir_w_k: BTreeMap<CellId, f64>,
    samples: Vec<ReservoirConvectionSample<M>>,
}
impl<M> ResolvedReservoirConvection<M> {
    pub fn reservoir_conductances_w_k(&self) -> &BTreeMap<CellId, f64> {
        &self.reservoir_w_k
    }
    pub fn samples(&self) -> &[ReservoirConvectionSample<M>] {
        &self.samples
    }

    /// Apply once to the native boundary replacements. All conflicts and
    /// cancellation are checked before mutating the destination. Independent
    /// ambient patches, authored sources and finite links retain their owners.
    pub fn apply_to(
        &self,
        destination: &mut BoundaryConductances,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), ConvectionSourceError> {
        for cell in self.reservoir_w_k.keys() {
            if cancelled() {
                return Err(ConvectionSourceError::Cancelled);
            }
            if destination.reservoir_w_k.contains_key(cell) {
                return Err(ConvectionSourceError::Invalid(
                    "convective reservoir already has a conductance owner",
                ));
            }
        }
        if cancelled() {
            return Err(ConvectionSourceError::Cancelled);
        }
        destination.reservoir_w_k.extend(
            self.reservoir_w_k.iter().map(|(&cell, &conductance)| (cell, conductance)),
        );
        Ok(())
    }
}

fn native_error(error: ConvectionSourceError) -> GraphReservoirConvectionError {
    GraphReservoirConvectionError::Convection(SourceCellConvectionError::Convection(error))
}

fn law_error(cell: CellId, error: FlowConvectionError) -> GraphReservoirConvectionError {
    match error {
        FlowConvectionError::Source(error) => GraphReservoirConvectionError::Source(error),
        FlowConvectionError::Geometry(source) => GraphReservoirConvectionError::Geometry { cell, source },
        FlowConvectionError::Envelope(source) => GraphReservoirConvectionError::Envelope { cell, source },
    }
}

#[allow(clippy::too_many_arguments)]
pub(in super::super) fn evaluate<M: ThermalFluidTransport>(
    core: &BoundGraphCore<'_>,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'_>,
    flow: GraphFlow<'_>,
    bindings: &[GraphConvectiveReservoir],
    maximum_reservoirs: usize,
    mut resolve: impl FnMut(Temperature) -> Result<M, SourceGraphError>,
) -> Result<ResolvedReservoirConvection<M>, GraphReservoirConvectionError> {
    let source_error = GraphReservoirConvectionError::Source;
    let invalid = |reason| native_error(ConvectionSourceError::Invalid(reason));
    poll(cx).map_err(source_error)?;
    if !(1..=1_000_000).contains(&maximum_reservoirs)
        || bindings.len() > maximum_reservoirs
        || bindings.len() > thermal.cells().len()
    {
        return Err(invalid("bounded active prescribed reservoir count required"));
    }
    let gravity = laws::gravity(core).map_err(source_error)?;
    let mut ids = BTreeSet::new();
    let mut ordered = Vec::with_capacity(bindings.len());
    // Qualify all original geometry and ownership before resolving film sources.
    for binding in bindings {
        poll(cx).map_err(source_error)?;
        if !ids.insert(binding.fluid_cell) {
            return Err(invalid("duplicate convective reservoir binding"));
        }
        let reservoir = thermal.cells()
            .get(&binding.fluid_cell)
            .and_then(|cell| cell.boundary.reservoir)
            .ok_or_else(|| invalid("convective wall requires an existing native reservoir"))?;
        if reservoir.conductance_w_k != 0.0 {
            return Err(invalid("convective reservoir requires zero authored conductance"));
        }
        if !positive(binding.characteristic_length.value()) {
            return Err(invalid("convective reservoir requires positive finite characteristic length"));
        }
        let law = laws::resolve(
            core, thermal, flow, binding.fluid_cell, binding.law, gravity,
        ).map_err(|error| law_error(binding.fluid_cell, error))?;
        law.admit_property_temperature(binding.property_temperature).map_err(native_error)?;
        ordered.push((binding, law));
    }
    ordered.sort_by_key(|(binding, _)| binding.fluid_cell);
    let mut reservoir_w_k = BTreeMap::new();
    let mut samples = Vec::with_capacity(ordered.len());
    for (binding, law) in ordered {
        poll(cx).map_err(source_error)?;
        let fluid_cell = binding.fluid_cell;
        let reservoir = thermal.cells()
            .get(&fluid_cell)
            .and_then(|cell| cell.boundary.reservoir)
            .ok_or_else(|| invalid("missing admitted prescribed wall reservoir"))?;
        let wall_temperature_k = reservoir.temperature_k;
        let fluid_temperature_k = thermal.temperature_k(fluid_cell)
            .map_err(|error| native_error(ConvectionSourceError::from(error)))?;
        let property_temperature_k = match binding.property_temperature {
            ConvectionPropertyTemperature::Reservoir => fluid_temperature_k,
            ConvectionPropertyTemperature::FilmMean => {
                0.5 * wall_temperature_k + 0.5 * fluid_temperature_k
            }
        };
        let material = resolve(Temperature::new(property_temperature_k))
            .map_err(|source| GraphReservoirConvectionError::Convection(
                SourceCellConvectionError::Source {
                    cell: fluid_cell,
                    property_temperature_k,
                    source,
                },
            ))?;
        let snapshot = material.thermal_transport_state()
            .map_err(|error| native_error(ConvectionSourceError::Material(error.into())))?;
        let transport = material.transport_state()
            .map_err(|error| native_error(ConvectionSourceError::Material(error.into())))?;
        if snapshot.temperature_k() != property_temperature_k || snapshot.transport() != transport {
            return Err(invalid(
                "fluid source does not match the requested temperature and transport state",
            ));
        }
        core.admit_density(fluid_cell, transport).map_err(source_error)?;
        let evaluation = evaluate_convection_state(
            &transport,
            wall_temperature_k,
            Temperature::new(fluid_temperature_k),
            binding.characteristic_length,
            law.resolve(&snapshot).map_err(native_error)?,
        ).map_err(|error| native_error(ConvectionSourceError::Material(error)))?;
        // Identical actual-source and signed old/new/crossing admission to
        // finite walls; maximum-speed envelopes are never substituted for Q.
        laws::admit_envelope(
            core,
            thermal,
            flow,
            FlowConvectionPoint {
                fluid_cell,
                characteristic_length: binding.characteristic_length,
                law: binding.law,
                wall_temperature_k,
                fluid_temperature_k,
                material: &material,
                thermal: &snapshot,
            },
        ).map_err(|error| law_error(fluid_cell, error))?;
        let conductance_w_k = match &evaluation {
            Some(sample) => reservoir.area_m2 * sample.coefficient.value.value(),
            None => 0.0,
        };
        let difference = wall_temperature_k - fluid_temperature_k;
        let heat_into_fluid_w = conductance_w_k * difference;
        if !nonnegative(conductance_w_k) || !difference.is_finite() || !heat_into_fluid_w.is_finite()
            || (evaluation.is_some() && reservoir.area_m2 > 0.0 && conductance_w_k == 0.0)
            || (conductance_w_k != 0.0 && difference != 0.0 && heat_into_fluid_w == 0.0)
        {
            return Err(invalid("convective reservoir transfer is unrepresentable"));
        }
        reservoir_w_k.insert(fluid_cell, conductance_w_k);
        samples.push(ReservoirConvectionSample {
            fluid_cell,
            wall_temperature_k,
            fluid_temperature_k,
            property_temperature_k,
            area_m2: reservoir.area_m2,
            conductance_w_k,
            heat_into_fluid_w,
            material,
            thermal: snapshot,
            evaluation,
        });
    }
    poll(cx).map_err(source_error)?;
    Ok(ResolvedReservoirConvection { reservoir_w_k, samples })
}
