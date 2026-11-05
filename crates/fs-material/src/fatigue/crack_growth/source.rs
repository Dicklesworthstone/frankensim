//! Original material-card Paris fits with explicit full-cycle and range units.
//!
//! Stress intensity has a half-integer length exponent, which `Dims` cannot
//! represent. Source properties and the fit coordinate therefore explicitly
//! store its SQUARE in Pa² m; only the native adapter takes the square root.
//! Original values, uncertainty, observations and query receipts stay intact.
//! No static strength, toughness or reversed-cycle datum supplies a Paris fit.

use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_matdb::{
    ClaimSelection, MatDbError, MaterialAnswer, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, QueryPoint,
};
use fs_math::det;
use fs_qty::semantic::{QuantityKind, ValueForm};
use fs_qty::{Dims, Length, QuantitySpec, Temperature};

use super::{CrackGrowthError, ParisDomain, ParisLaw, ParisParameters};
use crate::fatigue::FatigueScope;

/// Pa² m = kg² m⁻¹ s⁻⁴. This is squared intensity, never pressure or energy.
pub const STRESS_INTENSITY_SQUARED_DIMS: Dims = Dims([-1, 2, -4, 0, 0, 0]);

/// Four original properties from one fit, in selection/receipt order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParisProperty {
    /// Square of the reference FULL stress-intensity RANGE, in Pa² m.
    ReferenceDeltaKSquared,
    /// Metres per FULL cycle at the reference range; not per reversal.
    ReferenceGrowthRate,
    /// Exponent of Delta-K, NOT the exponent of its squared source coordinate.
    Exponent,
    /// Square of the exclusive maximum-K validity cutoff, NOT inferred toughness.
    MaximumKSquared,
}

impl ParisProperty {
    pub const ALL: [Self; 4] = [
        Self::ReferenceDeltaKSquared,
        Self::ReferenceGrowthRate,
        Self::Exponent,
        Self::MaximumKSquared,
    ];
    pub const fn index(self) -> usize {
        self as usize
    }
    pub fn property_key(self) -> PropertyKey {
        let (name, dims) = match self {
            Self::ReferenceDeltaKSquared => (
                "paris-reference-delta-k-squared",
                STRESS_INTENSITY_SQUARED_DIMS,
            ),
            Self::ReferenceGrowthRate => ("paris-reference-growth-per-full-cycle", Length::DIMS),
            Self::Exponent => ("paris-exponent", Dims::NONE),
            Self::MaximumKSquared => ("paris-maximum-k-squared", STRESS_INTENSITY_SQUARED_DIMS),
        };
        PropertyKey::new(name, dims)
    }
}

/// The source must explicitly declare common typed T, R and squared-Delta-K
/// support. Labels describe applicability; they do not validate transfer to a
/// different specimen, process, environment, load spectrum or geometry.
#[derive(Debug, Clone, PartialEq)]
pub struct ParisSourceRequest {
    pub temperature_axis: String,
    pub stress_ratio_axis: String,
    pub delta_k_squared_axis: String,
    pub surface_and_geometry: String,
    pub environment: String,
    pub load_basis: String,
}

/// Original source claims and answers retained inside the immutable Paris law.
#[derive(Debug, Clone, PartialEq)]
pub struct ParisCalibrationSource {
    card_identity: ContentHash,
    material_state: MaterialStateId,
    request: ParisSourceRequest,
    query: QueryPoint,
    selections: [ClaimSelection; 4],
    properties: Vec<(PropertyClaim, MaterialAnswer)>,
}

impl ParisCalibrationSource {
    pub const fn card_identity(&self) -> ContentHash {
        self.card_identity
    }
    pub fn material_state(&self) -> &MaterialStateId {
        &self.material_state
    }
    pub fn request(&self) -> &ParisSourceRequest {
        &self.request
    }
    pub fn query(&self) -> &QueryPoint {
        &self.query
    }
    pub fn selections(&self) -> &[ClaimSelection; 4] {
        &self.selections
    }
    pub fn property(&self, property: ParisProperty) -> &MaterialAnswer {
        &self.properties[property.index()].1
    }
    pub fn claim(&self, property: ParisProperty) -> &PropertyClaim {
        &self.properties[property.index()].0
    }

    /// Require the original card, named material state and every fixed physical
    /// coordinate with its exact quantity convention. Only calibration R and
    /// squared-Delta-K coordinates are excluded: the runtime law checks those.
    /// Additional native coordinates confer no additional calibration support.
    pub fn matches_material_state(
        &self,
        card_identity: ContentHash,
        material: &MaterialStateId,
        point: &QueryPoint,
    ) -> bool {
        self.card_identity == card_identity
            && &self.material_state == material
            && self.query.axes().iter().all(|(axis, value)| {
                axis == &self.request.stress_ratio_axis
                    || axis == &self.request.delta_k_squared_axis
                    || (point.axes().get(axis) == Some(value)
                        && point.axis_quantities().get(axis)
                            == self.query.axis_quantities().get(axis))
            })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum ParisSourceError {
    Source(MatDbError),
    Model(CrackGrowthError),
    Invalid(&'static str),
}
impl core::fmt::Display for ParisSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound Paris calibration refused: {self:?}")
    }
}
impl std::error::Error for ParisSourceError {}
impl From<MatDbError> for ParisSourceError {
    fn from(error: MatDbError) -> Self {
        Self::Source(error)
    }
}
impl From<CrackGrowthError> for ParisSourceError {
    fn from(error: CrackGrowthError) -> Self {
        Self::Model(error)
    }
}
fn require(condition: bool, reason: &'static str) -> Result<(), ParisSourceError> {
    if condition {
        Ok(())
    } else {
        Err(ParisSourceError::Invalid(reason))
    }
}

/// Resolve one ORIGINAL normalized Paris fit and its source-declared support.
/// Selections follow `ParisProperty::ALL`; pins may disambiguate each claim.
/// All four claims must retain one provenance, observation population and
/// identical typed validity domain. Scalars explicitly valid over that box and
/// temperature-only curves are admitted. R/range-dependent coefficients and
/// surfaces refuse; their bounding boxes cannot replace actual joint support.
///
/// Temperature is fixed at the queried semantic absolute-kelvin coordinate.
/// Source interpolation is used unchanged, with no fit, inferred fatigue data,
/// extrapolation, new uncertainty model or claim of component life. The upper
/// K cutoff is a source-declared model limit, not a fracture-toughness claim.
pub fn resolve_paris(
    card: &MaterialCard,
    point: &QueryPoint,
    request: &ParisSourceRequest,
    selections: [ClaimSelection; 4],
) -> Result<ParisLaw, ParisSourceError> {
    let axes = [
        &request.temperature_axis,
        &request.stress_ratio_axis,
        &request.delta_k_squared_axis,
    ];
    require(
        axes.iter().all(|axis| !axis.trim().is_empty())
            && axes[0] != axes[1]
            && axes[0] != axes[2]
            && axes[1] != axes[2],
        "declare distinct temperature, stress-ratio and squared-range axes",
    )?;
    let temperature = point
        .axis_quantities()
        .get(&request.temperature_axis)
        .ok_or(ParisSourceError::Invalid(
            "typed source temperature is absent",
        ))?;
    require(
        temperature.semantic_type().is_some_and(|quantity| {
            quantity.kind() == QuantityKind::AbsoluteTemperature
                && quantity.form() == ValueForm::Static
        }),
        "Paris source temperature must be semantic static absolute kelvin",
    )?;
    let quantities = [
        *temperature,
        QuantitySpec::dimensional(Dims::NONE),
        QuantitySpec::dimensional(STRESS_INTENSITY_SQUARED_DIMS),
    ];
    require(
        axes.iter().zip(quantities).all(|(axis, quantity)| {
            point.axis_quantities().get(*axis) == Some(&quantity)
                && point.axes().contains_key(*axis)
        }),
        "query T, tensile R and squared-Delta-K need their explicit quantity conventions",
    )?;
    let mut properties: Vec<(PropertyClaim, MaterialAnswer)> = Vec::with_capacity(4);
    for property in ParisProperty::ALL {
        let key = property.property_key();
        let answer = match selections[property.index()] {
            ClaimSelection::Policy(policy) => card.claims().query_typed(&key, point, policy)?,
            ClaimSelection::Pinned(pin) => card.claims().query_pinned_typed(&key, point, pin)?,
        };
        card.claims().verify_receipt(&answer.receipt)?;
        let claim =
            card.claims()
                .claim(answer.receipt.selected)
                .ok_or(ParisSourceError::Invalid(
                    "original Paris fit claim is absent",
                ))?;
        require(
            axes.iter().zip(quantities).all(|(axis, quantity)| {
                claim.validity.bounds().contains_key(*axis)
                    && claim.validity.axis_quantities().get(*axis) == Some(&quantity)
            }),
            "every Paris fit property needs explicit matching typed T, R and squared-range support",
        )?;
        require(
            match &claim.value {
                PropertyValue::Scalar { .. } => true,
                PropertyValue::Curve {
                    abscissa,
                    abscissa_dims,
                    ..
                } => abscissa == &request.temperature_axis && *abscissa_dims == Temperature::DIMS,
                _ => false,
            },
            "Paris parameters must be fixed scalars or functions of the declared temperature",
        )?;
        if let Some((first, _)) = properties.first() {
            require(
                first.provenance == claim.provenance
                    && first.observations == claim.observations
                    && first.validity == claim.validity,
                "Paris parameters must retain one source fit, population and complete typed domain",
            )?;
        }
        properties.push((claim.clone(), answer));
    }
    let first = &properties[0].0;
    let range_squared = first.validity.bounds()[&request.delta_k_squared_axis];
    let ratio = first.validity.bounds()[&request.stress_ratio_axis];
    let value = |property: ParisProperty| properties[property.index()].1.evidence.value.value;
    require(
        range_squared.0 > 0.0
            && range_squared.0 < range_squared.1
            && value(ParisProperty::ReferenceDeltaKSquared) > 0.0
            && value(ParisProperty::MaximumKSquared) > 0.0,
        "squared reference, fit interval and upper-K limit must be strictly positive",
    )?;
    let mut law = ParisLaw::new(
        FatigueScope {
            source_reference: first.provenance.source.clone(),
            material: card.id().chemistry.clone(),
            process: card.id().process.clone(),
            surface_and_geometry: request.surface_and_geometry.clone(),
            environment: request.environment.clone(),
            load_basis: request.load_basis.clone(),
            temperature_k: point.axes()[&request.temperature_axis],
        },
        ParisParameters {
            reference_delta_k_pa_sqrt_m: det::sqrt(value(ParisProperty::ReferenceDeltaKSquared)),
            reference_rate_m_per_cycle: value(ParisProperty::ReferenceGrowthRate),
            exponent: value(ParisProperty::Exponent),
        },
        ParisDomain {
            minimum_delta_k_pa_sqrt_m: det::sqrt(range_squared.0),
            maximum_delta_k_pa_sqrt_m: det::sqrt(range_squared.1),
            minimum_load_ratio: ratio.0,
            maximum_load_ratio: ratio.1,
            maximum_k_pa_sqrt_m: det::sqrt(value(ParisProperty::MaximumKSquared)),
        },
    )?;
    law.source = Some(Arc::new(ParisCalibrationSource {
        card_identity: card.content_hash(),
        material_state: card.id().clone(),
        request: request.clone(),
        query: point.clone(),
        selections,
        properties,
    }));
    Ok(law)
}
