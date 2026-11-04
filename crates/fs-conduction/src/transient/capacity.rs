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
use super::backward_euler::BackwardEuler;
use crate::{ConductionError, ConductivityModel, ConductivityTable, material::TemperatureSpan};
use crate::{ConductionMesh, ElementMaterials, MaterialId};
use fs_exec::Cx;
use std::collections::BTreeMap;

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

/// Explicit constant heat-capacity distribution, in J/(m^3 K).
#[derive(Debug, Clone, Copy)]
pub enum HeatCapacity<'a> {
    /// One declared value for the whole mesh, retaining the original path.
    Uniform(VolumetricHeatCapacity),
    /// Exactly one declared value per tetrahedron, in the original mesh order.
    PerElement(&'a [VolumetricHeatCapacity]),
    /// Exactly one value per material actually used by `element_materials`.
    /// Input ordering is irrelevant. No missing value inherits a fallback, and
    /// duplicate or unused IDs refuse instead of silently changing a model.
    ByMaterial(&'a [(MaterialId, VolumetricHeatCapacity)]),
    /// Uniform density and specific heat selected from one material card.
    /// Source validity is checked against initial and accepted nodal fields.
    MaterialUniform(&'a MaterialHeatCapacity),
    /// Source-derived rho*c_p using the SAME material IDs as conductivity.
    /// Each source remains distinct; no representative material is averaged in.
    MaterialById(&'a [(MaterialId, MaterialHeatCapacity)]),
}

impl HeatCapacity<'_> {
    pub(crate) fn prepare<'m>(
        self,
        cx: &Cx<'_>,
        mesh: &'m ConductionMesh,
        materials: Option<&ElementMaterials>,
    ) -> Result<BackwardEuler<'m>, ConductionError> {
        poll(cx, 0)?;
        match self {
            Self::Uniform(value) => BackwardEuler::uniform(cx, mesh, value),
            Self::PerElement(values) => BackwardEuler::per_element(cx, mesh, values),
            Self::MaterialUniform(source) => BackwardEuler::uniform_material(cx, mesh, source),
            Self::MaterialById(rows) => {
                let materials = materials.ok_or_else(|| {
                    invalid_assignment(
                        "material-bound capacity requires an explicit element material assignment"
                            .into(),
                    )
                })?;
                materials.validate_for(mesh)?;
                let mut indexed = BTreeMap::new();
                let mut sources = Vec::with_capacity(rows.len());
                for (slot, (id, source)) in rows.iter().enumerate() {
                    if slot % 512 == 0 {
                        poll(cx, slot)?;
                    }
                    if indexed.insert(*id, slot).is_some() {
                        return Err(invalid_assignment(format!(
                            "duplicate heat-capacity source for material {id:?}"
                        )));
                    }
                    sources.push(source.clone());
                }
                let mut assigned = Vec::with_capacity(mesh.element_count());
                for (element, id) in materials.of_element().iter().enumerate() {
                    if element % 512 == 0 {
                        poll(cx, element)?;
                    }
                    assigned.push(*indexed.get(id).ok_or_else(|| invalid_assignment(format!(
                        "element {element} uses material {id:?}, which has no heat-capacity source")))?);
                }
                BackwardEuler::material_assigned(cx, mesh, &sources, &assigned)
            }
            Self::ByMaterial(rows) => {
                let materials = materials.ok_or_else(|| {
                    invalid_assignment(
                        "material-bound capacity requires an explicit element material assignment"
                            .into(),
                    )
                })?;
                materials.validate_for(mesh)?;
                let mut indexed = BTreeMap::new();
                for (slot, &(id, value)) in rows.iter().enumerate() {
                    if slot % 512 == 0 {
                        poll(cx, slot)?;
                    }
                    if indexed.insert(id, (value, false)).is_some() {
                        return Err(invalid_assignment(format!(
                            "duplicate heat-capacity declaration for material {id:?}"
                        )));
                    }
                }
                let mut values = Vec::with_capacity(mesh.element_count());
                for (element, id) in materials.of_element().iter().enumerate() {
                    if element % 512 == 0 {
                        poll(cx, element)?;
                    }
                    let (value, used) = indexed.get_mut(id).ok_or_else(|| invalid_assignment(format!(
                        "element {element} uses material {id:?}, which has no declared heat capacity")))?;
                    *used = true;
                    values.push(*value);
                }
                for (slot, (id, (_, used))) in indexed.iter().enumerate() {
                    if slot % 512 == 0 {
                        poll(cx, slot)?;
                    }
                    if !*used {
                        return Err(invalid_assignment(format!(
                            "heat-capacity declaration names unused material {id:?}"
                        )));
                    }
                }
                poll(cx, values.len())?;
                BackwardEuler::per_element(cx, mesh, &values)
            }
        }
    }
}

fn invalid_assignment(what: String) -> ConductionError {
    ConductionError::Config {
        parameter: "heat-capacity assignment",
        what,
    }
}

fn poll(cx: &Cx<'_>, at: usize) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "heat-capacity-prepare",
        at,
    })
}
