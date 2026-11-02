//! One material-card binding for stationary conducting walls. Parsing and file
//! access stay with callers; caloric integration and k(T) remain material owners.

use super::{MaterialIngressError, MaterialRecord};
use fs_matdb::{
    ClaimSelection, NormalizedMaterialCardPack, PropertyKey, QueryPoint, SelectionPolicy,
};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::thermal::conduction::source::ResolvedConductionLayer;
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_thermal::enthalpy::DepositionError;
use std::sync::Arc;

/// Explicit source coordinates and property names; no material constants or
/// temperature extrapolation are supplied by a performance reader.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialWallSource {
    pub material: MaterialRecord,
    pub temperature_axis: String,
    pub phase: SolidLiquidPhase,
    /// Density, specific heat capacity and thermal conductivity, in that order.
    pub property_names: [String; 3],
    pub temperature_span_k: [f64; 2],
    pub caloric_interpolation_tolerance_j_kg: f64,
}

/// The immutable original pack and its admitted pure-phase caloric chart.
/// Keeping this owner alive lets several cells borrow the same conductivity
/// law without self-references or a second thermal solver.
#[derive(Clone, Debug)]
pub struct ResolvedMaterialWallSource {
    source: MaterialWallSource,
    pack: Arc<NormalizedMaterialCardPack>,
    context: QueryPoint,
    caloric: ResolvedHeatCapacityCurve,
    conductivity: PropertyKey,
}

impl ResolvedMaterialWallSource {
    pub fn new(
        source: MaterialWallSource,
        pack: Arc<NormalizedMaterialCardPack>,
    ) -> Result<Self, MaterialIngressError> {
        let fail = |error: &dyn core::fmt::Display| MaterialIngressError(error.to_string());
        let valid_token = |name: &str| !name.is_empty() && name.len() <= 128;
        if !valid_token(&source.temperature_axis)
            || source
                .material
                .axes
                .iter()
                .any(|(axis, _)| axis == &source.temperature_axis)
            || !source.property_names.iter().all(|name| valid_token(name))
            || !source
                .temperature_span_k
                .iter()
                .all(|t| t.is_finite() && *t > 0.0)
            || source.temperature_span_k[0] >= source.temperature_span_k[1]
            || !source.caloric_interpolation_tolerance_j_kg.is_finite()
            || source.caloric_interpolation_tolerance_j_kg <= 0.0
        {
            return Err(MaterialIngressError("wall source needs bounded names, a distinct temperature axis, increasing positive kelvin and a positive caloric tolerance".into()));
        }
        let card = source.material.card(&pack)?;
        let phase_name = match source.phase {
            SolidLiquidPhase::Solid => "solid",
            SolidLiquidPhase::Liquid => "liquid",
            SolidLiquidPhase::SolidLiquid => {
                return Err(MaterialIngressError(
                    "stationary conducting walls require an explicit pure phase".into(),
                ));
            }
        };
        if card.id().phase != phase_name {
            return Err(MaterialIngressError(
                "wall phase must match the compiled material-card phase".into(),
            ));
        }
        let mut lower_record = source.material.clone();
        lower_record.axes.push((
            source.temperature_axis.clone(),
            source.temperature_span_k[0],
        ));
        let context = lower_record.query_point(card)?;
        let upper = context
            .clone()
            .with(&source.temperature_axis, source.temperature_span_k[1])
            .map_err(|error| fail(&error))?;
        let key = |name: &str| {
            card.claims()
                .registered_quantity(name)
                .map(|quantity| PropertyKey::with_quantity(name, quantity))
                .ok_or_else(|| {
                    MaterialIngressError(format!("wall material lacks property {name:?}"))
                })
        };
        let density = key(&source.property_names[0])?;
        let heat_capacity = key(&source.property_names[1])?;
        let conductivity = key(&source.property_names[2])?;
        let caloric = ResolvedHeatCapacityCurve::try_new(
            card,
            HeatCapacitySourceQuery {
                heat_capacity: &heat_capacity,
                density: &density,
                temperature_axis: &source.temperature_axis,
                lower: &context,
                upper: &upper,
                selection: SelectionPolicy::SingleClaimOnly,
            },
            source.phase,
            0.0,
            source.caloric_interpolation_tolerance_j_kg,
        )
        .map_err(|error| fail(&error))?;
        Ok(Self {
            source,
            pack,
            context,
            caloric,
            conductivity,
        })
    }

    pub fn source(&self) -> &MaterialWallSource {
        &self.source
    }
    pub fn pack(&self) -> &NormalizedMaterialCardPack {
        &self.pack
    }

    /// Resolve the original conductivity line on exactly the retained caloric
    /// range. No k(T), Cp(T), source quantity or fixed coordinate is replaced.
    pub fn layer(
        &self,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ResolvedConductionLayer<'_>, MaterialIngressError> {
        ResolvedConductionLayer::new(
            self.caloric.curve(),
            ThermalMaterialPoint {
                material: self.pack.card(),
                context: &self.context,
                temperature_axis: &self.source.temperature_axis,
                liquid_fraction_axis: None,
            },
            &self.conductivity,
            ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            self.source.phase,
            self.source.temperature_span_k,
            checkpoint,
        )
        .map_err(|error| MaterialIngressError(error.to_string()))
    }
}
