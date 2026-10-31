//! Continue original surface sources onto an accepted, locally heated solid.

use fs_conduction::ConductionError;
use fs_conduction::radiation::{SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity};
use fs_matdb::MaterialCard;
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, VisibleOpticalStatePoint,
};
use fs_material::thermal::source::ThermalMaterialPoint;

use super::{
    Cx, MaterialGeometryError, MaterialRegionAppearance, MaterialRegionBinding,
    MaterialRegionBindingError, Primitive, TetAssemblyBudget, ThermalBody, bind_ordered_regions,
    deformed_body_regions, poll,
};
use crate::material_state::{
    MaterialOpticalBinding, MaterialOpticalRefreshError, MaterialThermalBinding,
    MaterialThermalError,
};

/// Original source and independent surface state for an exterior material region.
/// The source is retained after every accepted observation; no endpoint performs
/// fresh claim discovery or switches the declared appearance model.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialRegionSource {
    /// Bulk optical claims paired with the exact binding they produced.
    Optical {
        /// Original accepted optical observation and independent surface finish.
        binding: MaterialOpticalBinding,
        /// Original optical family and all selected material claim receipts.
        source: VisibleOpticalStatePoint,
    },
    /// Explicit opaque gray/diffuse observation of total emissivity.
    Thermal {
        /// Retains the exact emissivity claim and independent surface identity.
        binding: MaterialThermalBinding,
        /// Existing absolute-temperature coordinate of the original source.
        temperature_axis: String,
    },
}

impl MaterialRegionSource {
    fn appearance(&self) -> MaterialRegionAppearance<'_> {
        match self {
            Self::Optical { binding, .. } => MaterialRegionAppearance::Optical(binding),
            Self::Thermal { binding, .. } => MaterialRegionAppearance::Thermal(binding),
        }
    }

    fn at_material_state(
        &self,
        region: &str,
        card: &MaterialCard,
        accepted: &ResolvedMaterialStatePoint,
    ) -> Result<Self, MaterialRegionContinuationError> {
        match self {
            Self::Optical { binding, source } => {
                let (binding, source) =
                    binding
                        .at_material_state(source, card, accepted)
                        .map_err(|error| MaterialRegionContinuationError::Optical {
                            region: region.to_owned(),
                            error,
                        })?;
                Ok(Self::Optical { binding, source })
            }
            Self::Thermal {
                binding,
                temperature_axis,
            } => {
                let original = binding.emissivity();
                let receipt = original.receipt();
                let original_temperature = receipt
                    .query_point
                    .iter()
                    .find(|(axis, _)| axis == temperature_axis)
                    .map(|(_, value)| *value);
                if original.card_identity() != card.content_hash()
                    || accepted.card_identity() != original.card_identity()
                    || accepted.axis_quantities() != &receipt.axis_quantities
                    || !accepted
                        .query_point()
                        .iter()
                        .map(|(axis, _)| axis)
                        .eq(receipt.query_point.iter().map(|(axis, _)| axis))
                    || original_temperature != Some(original.temperature_k())
                {
                    return Err(MaterialRegionBindingError::StateMismatch {
                        region: region.to_owned(),
                    }
                    .into());
                }
                // The accepted solid supplies the complete physical condition;
                // the original radiation owner retains units/static admission.
                let point = accepted.query();
                let (emissivity, _) = SurfaceEmissivity::from_condition(
                    region,
                    ThermalMaterialPoint {
                        material: card,
                        context: &point,
                        temperature_axis,
                        liquid_fraction_axis: None,
                    },
                    MaterialPropertySelection::PinnedByProperty(vec![(
                        SURFACE_EMISSIVITY_PROPERTY.to_owned(),
                        receipt.selected,
                    )]),
                )
                .map_err(|error| {
                    MaterialRegionContinuationError::ThermalSource {
                        region: region.to_owned(),
                        error,
                    }
                })?;
                let binding = MaterialThermalBinding::try_gray_diffuse(
                    accepted,
                    &emissivity,
                    binding.surface_state_identity(),
                )
                .map_err(|error| {
                    MaterialRegionContinuationError::ThermalBinding {
                        region: region.to_owned(),
                        error,
                    }
                })?;
                Ok(Self::Thermal {
                    binding,
                    temperature_axis: temperature_axis.clone(),
                })
            }
        }
    }
}

/// One exterior region's retained source for the next accepted solid endpoint.
#[derive(Clone, Debug)]
pub struct MaterialRegionContinuation<'a> {
    /// Exact reference-region key, independent of input/element order.
    pub key: &'a str,
    /// The immutable original material card, including competing claims.
    pub card: &'a MaterialCard,
    /// Previously accepted observation and its exact original sources.
    pub source: MaterialRegionSource,
}

/// Refusal before any successor sources or replacement primitives are exposed.
#[derive(Debug)]
pub enum MaterialRegionContinuationError {
    /// Original geometry, regional assignment or accepted-state admission.
    Binding(MaterialRegionBindingError),
    /// The region's original optical law cannot support the accepted condition.
    Optical {
        /// Exterior region whose original source refused.
        region: String,
        /// Original optical owner refusal.
        error: MaterialOpticalRefreshError,
    },
    /// The original total-emissivity claim cannot support the accepted condition.
    ThermalSource {
        /// Exterior region whose original source refused.
        region: String,
        /// Original radiation owner refusal.
        error: ConductionError,
    },
    /// The refreshed thermal observation does not match the accepted state.
    ThermalBinding {
        /// Exterior region whose observation refused.
        region: String,
        /// Original renderer admission refusal.
        error: MaterialThermalError,
    },
}

impl core::fmt::Display for MaterialRegionContinuationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Binding(error) => error.fmt(f),
            Self::Optical { region, error } => {
                write!(f, "region {region} optical continuation: {error}")
            }
            Self::ThermalSource { region, error } => {
                write!(f, "region {region} emissivity continuation: {error}")
            }
            Self::ThermalBinding { region, error } => {
                write!(f, "region {region} thermal continuation: {error}")
            }
        }
    }
}

impl core::error::Error for MaterialRegionContinuationError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        Some(match self {
            Self::Binding(error) => error,
            Self::Optical { error, .. } => error,
            Self::ThermalSource { error, .. } => error,
            Self::ThermalBinding { error, .. } => error,
        })
    }
}

impl From<MaterialRegionBindingError> for MaterialRegionContinuationError {
    fn from(error: MaterialRegionBindingError) -> Self {
        Self::Binding(error)
    }
}

impl From<MaterialGeometryError> for MaterialRegionContinuationError {
    fn from(error: MaterialGeometryError) -> Self {
        Self::Binding(error.into())
    }
}

/// Continue a complete native renderer observation onto an accepted thermal body.
///
/// Every exterior region evaluates its original optical/emissivity claims once
/// at that region's complete accepted constitutive condition. The existing
/// exterior owner extracts geometry once and validates all returned primitives.
/// Sources and primitives return together in canonical exterior-region order;
/// callers publish them together only on success and retain sources for the next
/// endpoint. No input scene, physical state or original binding is mutated.
///
/// This observes a fixed-topology opaque solid. It preserves each independently
/// declared surface finish/identity and does not infer surface evolution, a new
/// heat debit, or support along the intervening physical path. Optical families
/// remain intact; the existing regional owner refuses transparent interfaces.
/// Wholly interior regions do not acquire optical or emissivity requirements.
///
/// # Errors
/// Refuses stale/foreign source pairs, schema changes, unsupported original
/// claims, incomplete regional assignments and original geometry/cancellation
/// failures without returning any partial successor or scene replacement.
pub fn continue_thermal_body_primitives<'a>(
    body: &ThermalBody<'_>,
    bindings: &[MaterialRegionContinuation<'a>],
    budget: TetAssemblyBudget,
    maximum_face_entries: usize,
    cx: &Cx<'_>,
) -> Result<(Vec<Primitive>, Vec<MaterialRegionContinuation<'a>>), MaterialRegionContinuationError>
{
    let regions = deformed_body_regions(body.body(), budget, maximum_face_entries, cx)?;
    if bindings.len() != regions.len() {
        return Err(MaterialRegionBindingError::Bindings(
            "exactly one source per exterior region is required",
        )
        .into());
    }
    let mut ordered: Vec<_> = bindings.iter().collect();
    ordered.sort_by_key(|binding| binding.key);
    // Validate the whole assignment before resolving any source. Source input
    // order cannot swap ownership, and duplicate keys cannot hide missing cells.
    for (region, binding) in regions.iter().zip(&ordered) {
        poll(cx)?;
        if binding.key != body.body().reference().regions()[region.region_index].key() {
            return Err(MaterialRegionBindingError::Bindings(
                "source keys must match the exact exterior region inventory",
            )
            .into());
        }
    }
    let mut successor = Vec::with_capacity(regions.len());
    for (region, binding) in regions.iter().zip(ordered) {
        poll(cx)?;
        let accepted = body.sources()[region.region_index].stiffness();
        successor.push(MaterialRegionContinuation {
            key: binding.key,
            card: binding.card,
            source: binding
                .source
                .at_material_state(binding.key, binding.card, accepted)?,
        });
    }
    let appearances: Vec<_> = successor
        .iter()
        .map(|binding| MaterialRegionBinding {
            key: binding.key,
            appearance: binding.source.appearance(),
        })
        .collect();
    let primitives = bind_ordered_regions(body, regions, &appearances, cx)?;
    Ok((primitives, successor))
}
