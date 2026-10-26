//! Source-bound constant heat storage for the existing temperature steppers.
//!
//! This reduction is admitted only when BOTH density and specific heat are
//! validity-wide scalar laws. Temperature-varying data belong to the existing
//! reference-mass enthalpy owner; they are never silently frozen here.

use fs_matdb::{MaterialCard, PropertyKey, QueryPoint, SelectionPolicy};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, ResolvedMaterialStatePoint,
    SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY, ScalarAdmissibility,
    ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
    resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, QuantitySpec, Temperature};

use super::VolumetricHeatCapacity;
use crate::{ConductionError, ConductivityModel, ConductivityTable, material::TemperatureSpan};

/// Density and specific heat from one immutable material card, with an explicit
/// temperature interval and unchanged external pressure/composition conditions.
/// The original source decisions at the reference and both endpoints survive
/// preparation; no material name, conductivity, or representative preset supplies
/// an absent caloric coefficient.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialHeatCapacity {
    reference: ResolvedMaterialStatePoint,
    endpoints: [ResolvedMaterialStatePoint; 2],
    capacity: VolumetricHeatCapacity,
    span: TemperatureSpan,
}

impl MaterialHeatCapacity {
    /// Resolve canonical density [kg/m³] and specific heat [J/(kg K)] at `point`
    /// and pin those exact claims over `[low_k, high_k]`. The source coordinate
    /// `T` must denote absolute temperature; all other conditions stay fixed.
    ///
    /// Only scalar, static source laws are admitted. A curve is not made
    /// constant by sampling it once, even if its endpoint values happen to agree.
    /// Nonconstant caloric data require `EnthalpyBackwardEuler` and its original
    /// `ResolvedHeatCapacityCurve`. No phase boundary is inferred from a card's
    /// name: the caller must select a single-phase source validity interval.
    pub fn from_card(
        card: &MaterialCard,
        point: &QueryPoint,
        low_k: f64,
        high_k: f64,
        selection: MaterialPropertySelection,
    ) -> Result<Self, ConductionError> {
        if !low_k.is_finite() || !high_k.is_finite() || low_k <= 0.0 || low_k >= high_k {
            return Err(invalid(
                "a finite positive temperature interval is required",
            ));
        }
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        if point.axis_quantities().get("T").is_some_and(|quantity| {
            *quantity != absolute && *quantity != QuantitySpec::dimensional(Temperature::DIMS)
        }) || !point
            .axes()
            .get("T")
            .is_some_and(|value| (low_k..=high_k).contains(value))
            || point
                .axis_quantities()
                .iter()
                .any(|(axis, quantity)| axis != "T" && quantity.dims() == Temperature::DIMS)
        {
            return Err(invalid(
                "the reference requires one absolute-temperature coordinate T inside the interval",
            ));
        }
        let requirements = [
            (DENSITY_PROPERTY, Density::DIMS),
            (SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS),
        ]
        .into_iter()
        .map(|(name, dims)| {
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
        })
        .collect::<Result<Vec<_>, _>>()
        .map_err(|error| invalid(&error.to_string()))?;
        let reference = resolve_material_state_point(card, point, &requirements, selection)
            .map_err(|error| invalid(&error.to_string()))?;
        for property in reference.properties() {
            if !property.source_axes().is_empty()
                || !property.is_static_at_fixed_external_conditions()
            {
                return Err(invalid(
                    "constant storage requires static scalar density and specific heat; use source-bound enthalpy for varying laws",
                ));
            }
        }
        let resolve_endpoint = |temperature| {
            let query = point
                .clone()
                .with("T", temperature)
                .map_err(|error| invalid(&error.to_string()))?;
            reference
                .at_query_point(card, &query)
                .map_err(|error| invalid(&error.to_string()))
        };
        let endpoints = [resolve_endpoint(low_k)?, resolve_endpoint(high_k)?];
        let density = reference
            .property(DENSITY_PROPERTY)
            .ok_or_else(|| invalid("density source is missing"))?
            .value_si();
        let specific_heat = reference
            .property(SPECIFIC_HEAT_CAPACITY_PROPERTY)
            .ok_or_else(|| invalid("specific-heat source is missing"))?
            .value_si();
        let capacity = VolumetricHeatCapacity::declared(density * specific_heat)?;
        Ok(Self {
            reference,
            endpoints,
            capacity,
            span: TemperatureSpan::Sampled {
                low: low_k,
                high: high_k,
            },
        })
    }

    /// Exact reference material state, including density and specific-heat receipts.
    #[must_use]
    pub const fn reference(&self) -> &ResolvedMaterialStatePoint {
        &self.reference
    }

    /// Pinned source decisions at the two temperature endpoints, lower first.
    #[must_use]
    pub const fn endpoints(&self) -> &[ResolvedMaterialStatePoint; 2] {
        &self.endpoints
    }

    /// Constant rho*c_p [J/(m³ K)] used by the original P1 capacity assembler.
    #[must_use]
    pub const fn value_j_per_m3_k(&self) -> f64 {
        self.capacity.value_j_per_m3_k()
    }

    /// Closed single-phase temperature interval supplied at preparation.
    #[must_use]
    pub const fn temperature_span(&self) -> TemperatureSpan {
        self.span
    }

    /// Build the existing isotropic k(T) transport model from THIS card and the
    /// SAME external conditions/temperature interval as heat storage. Source
    /// breakpoints remain owned by `ConductivityTable`; k need not be constant.
    /// This explicit pairing prevents a material swap retaining another card's
    /// conductivity. Manually assembled conductivity models remain independent.
    pub fn isotropic_conductivity(
        &self,
        card: &MaterialCard,
        selection: SelectionPolicy,
    ) -> Result<ConductivityModel, ConductionError> {
        if card.content_hash() != self.reference.card_identity() {
            return Err(invalid(
                "paired conductivity requires the original heat-capacity material card",
            ));
        }
        let quantity = if card
            .claims()
            .registered_quantity(THERMAL_CONDUCTIVITY_PROPERTY)
            .is_some_and(|quantity| quantity.semantic_type().is_some())
        {
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            ))
        } else {
            QuantitySpec::dimensional(THERMAL_CONDUCTIVITY_DIMS)
        };
        let key = PropertyKey::with_quantity(THERMAL_CONDUCTIVITY_PROPERTY, quantity);
        let points = [self.endpoints[0].query(), self.endpoints[1].query()];
        Ok(ConductivityModel::isotropic(
            ConductivityTable::from_claims_at_query_points(
                card.claims(),
                &key,
                "T",
                &points,
                selection,
            )?,
        ))
    }
}

fn invalid(what: &str) -> ConductionError {
    ConductionError::Config {
        parameter: "material heat capacity",
        what: what.to_owned(),
    }
}
