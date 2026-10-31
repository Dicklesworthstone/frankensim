//! Original material-card coefficients for one bounded CMB calibration.
//!
//! The four fitted fatigue parameters must name one common fit, observation
//! population, fully reversed stress baseline and finite life interval. The
//! elastic modulus is queried from the same immutable material card and fixed
//! physical context. No monotonic strength, hardness, ductility or hardening
//! parameter is substituted for a missing cyclic fatigue measurement.
//!
//! Coefficients fitted against full cycles N and reversals 2N are distinct:
//! a full-cycle coefficient A is converted to A*2^(-exponent) for the native
//! reversal-based CMB law. Original values and all receipts remain unchanged.

use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_matdb::{
    ClaimSelection, MatDbError, MaterialAnswer, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, QueryPoint,
};
use fs_math::det;
use fs_qty::semantic::{QuantityKind, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

use super::{CmbParameters, StrainLifeCurve};
use crate::fatigue::{FatigueError, FatigueScope};

/// The order of the five explicit per-property source selections.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CmbProperty {
    YoungsModulus,
    FatigueStrengthCoefficient,
    FatigueStrengthExponent,
    FatigueDuctilityCoefficient,
    FatigueDuctilityExponent,
}
impl CmbProperty {
    pub const ALL: [Self; 5] = [
        Self::YoungsModulus,
        Self::FatigueStrengthCoefficient,
        Self::FatigueStrengthExponent,
        Self::FatigueDuctilityCoefficient,
        Self::FatigueDuctilityExponent,
    ];
    pub const fn index(self) -> usize {
        self as usize
    }
    pub fn property_key(self) -> PropertyKey {
        let (name, dims) = match self {
            Self::YoungsModulus => ("young_modulus", Pressure::DIMS),
            Self::FatigueStrengthCoefficient => ("fatigue_strength_coefficient", Pressure::DIMS),
            Self::FatigueStrengthExponent => ("fatigue_strength_exponent", Dims::NONE),
            Self::FatigueDuctilityCoefficient => ("fatigue_ductility_coefficient", Dims::NONE),
            Self::FatigueDuctilityExponent => ("fatigue_ductility_exponent", Dims::NONE),
        };
        PropertyKey::new(name, dims)
    }
}

/// Independent variable in BOTH original fitted power terms, not just the
/// reporting convention for life. Total mechanical strain is always amplitude.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CmbPowerBasis {
    /// epsilon_a = sigma_c/E * N^b + epsilon_c * N^c.
    FullCycles,
    /// epsilon_a = sigma_f'/E * (2N)^b + epsilon_f' * (2N)^c.
    Reversals,
}
impl CmbPowerBasis {
    pub const fn life_axis(self) -> &'static str {
        match self {
            Self::FullCycles => "cycles-to-failure",
            Self::Reversals => "reversals-to-failure",
        }
    }
}

/// Applicability labels remain caller declarations. The actual temperature,
/// fully reversed baseline and life interval must also be in the source data.
#[derive(Debug, Clone, PartialEq)]
pub struct CmbSourceRequest {
    pub power_basis: CmbPowerBasis,
    pub temperature_axis: String,
    pub stress_ratio_axis: String,
    pub surface_and_geometry: String,
    pub environment: String,
    pub load_basis: String,
}

/// Complete original coefficients, fit context and query evidence. Stored in
/// the immutable curve, so counters and native observers retain the source.
#[derive(Debug, Clone, PartialEq)]
pub struct CmbCalibrationSource {
    card_identity: ContentHash,
    material_state: MaterialStateId,
    request: CmbSourceRequest,
    query: QueryPoint,
    selections: [ClaimSelection; 5],
    properties: Vec<(PropertyClaim, MaterialAnswer)>,
}
impl CmbCalibrationSource {
    pub const fn card_identity(&self) -> ContentHash {
        self.card_identity
    }
    pub fn material_state(&self) -> &MaterialStateId {
        &self.material_state
    }
    pub fn request(&self) -> &CmbSourceRequest {
        &self.request
    }
    pub fn query(&self) -> &QueryPoint {
        &self.query
    }
    pub fn selections(&self) -> &[ClaimSelection; 5] {
        &self.selections
    }
    pub fn property(&self, property: CmbProperty) -> &MaterialAnswer {
        &self.properties[property.index()].1
    }
    pub fn claim(&self, property: CmbProperty) -> &PropertyClaim {
        &self.properties[property.index()].0
    }

    /// Check against the native owner's ORIGINAL resolved state and modulus
    /// answer. Calibration R and life coordinates are not operating conditions;
    /// every other source-query coordinate and its quantity convention must
    /// agree. This admits no temperature/environment transfer or source switch.
    pub fn matches_material_state(
        &self,
        card_identity: ContentHash,
        material: &MaterialStateId,
        point: &QueryPoint,
        youngs: &MaterialAnswer,
    ) -> bool {
        let original = self.property(CmbProperty::YoungsModulus);
        self.card_identity == card_identity
            && &self.material_state == material
            && original.receipt.selected == youngs.receipt.selected
            && original.evidence.value == youngs.evidence.value
            && self.query.axes().iter().all(|(axis, value)| {
                axis == &self.request.stress_ratio_axis
                    || axis == self.request.power_basis.life_axis()
                    || point.axes().get(axis) == Some(value)
                        && point.axis_quantities().get(axis)
                            == self.query.axis_quantities().get(axis)
            })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum CmbSourceError {
    Source(MatDbError),
    Fatigue(FatigueError),
    Invalid(&'static str),
}
impl core::fmt::Display for CmbSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound CMB calibration refused: {self:?}")
    }
}
impl std::error::Error for CmbSourceError {}
impl From<MatDbError> for CmbSourceError {
    fn from(error: MatDbError) -> Self {
        Self::Source(error)
    }
}
impl From<FatigueError> for CmbSourceError {
    fn from(error: FatigueError) -> Self {
        Self::Fatigue(error)
    }
}
fn require(accepted: bool, reason: &'static str) -> Result<(), CmbSourceError> {
    if accepted {
        Ok(())
    } else {
        Err(CmbSourceError::Invalid(reason))
    }
}

/// Resolve a fully reversed total-mechanical-strain AMPLITUDE fit. Selections
/// follow `CmbProperty::ALL`; each may independently pin an original claim.
/// Four fatigue coefficients must have the same provenance, observations and
/// life bounds. Coefficients may be scalar or depend on the declared fixed T;
/// life-dependent fitted coefficients and unsupported source forms refuse.
/// The model's finite life interval is obtained from those original bounds.
pub fn resolve_cmb(
    card: &MaterialCard,
    point: &QueryPoint,
    request: &CmbSourceRequest,
    selections: [ClaimSelection; 5],
) -> Result<StrainLifeCurve, CmbSourceError> {
    let life_axis = request.power_basis.life_axis();
    require(
        !request.temperature_axis.trim().is_empty()
            && !request.stress_ratio_axis.trim().is_empty()
            && request.temperature_axis != request.stress_ratio_axis
            && request.temperature_axis != life_axis
            && request.stress_ratio_axis != life_axis,
        "declare distinct temperature, baseline stress-ratio and life axes",
    )?;
    let dimensionless = QuantitySpec::dimensional(Dims::NONE);
    let temperature = point
        .axis_quantities()
        .get(&request.temperature_axis)
        .ok_or(CmbSourceError::Invalid(
            "CMB source needs typed absolute temperature",
        ))?;
    require(
        *temperature == QuantitySpec::dimensional(Temperature::DIMS)
            || temperature.semantic_type().is_some_and(|q| {
                q.kind() == QuantityKind::AbsoluteTemperature && q.form() == ValueForm::Static
            }),
        "CMB source temperature must be absolute kelvin",
    )?;
    require(
        point.axis_quantities().get(life_axis) == Some(&dimensionless)
            && point.axis_quantities().get(&request.stress_ratio_axis) == Some(&dimensionless)
            && point.axes().get(&request.stress_ratio_axis) == Some(&-1.0),
        "declare typed source life and a fully reversed R = -1 calibration",
    )?;
    let temperature_k = *point
        .axes()
        .get(&request.temperature_axis)
        .ok_or(CmbSourceError::Invalid("CMB source temperature is absent"))?;
    let mut properties = Vec::with_capacity(5);
    for property in CmbProperty::ALL {
        let key = property.property_key();
        let answer = match selections[property.index()] {
            ClaimSelection::Policy(policy) => card.claims().query_typed(&key, point, policy)?,
            ClaimSelection::Pinned(pin) => card.claims().query_pinned_typed(&key, point, pin)?,
        };
        card.claims().verify_receipt(&answer.receipt)?;
        let claim = card
            .claims()
            .claim(answer.receipt.selected)
            .ok_or(CmbSourceError::Invalid(
                "original CMB source claim is absent",
            ))?;
        require(
            claim
                .validity
                .axis_quantities()
                .get(&request.temperature_axis)
                == Some(temperature)
                && claim
                    .validity
                    .bounds()
                    .contains_key(&request.temperature_axis),
            "every CMB parameter and modulus needs the same typed source temperature",
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
            "CMB coefficients must be fixed scalars or functions of the declared temperature",
        )?;
        if property != CmbProperty::YoungsModulus {
            require(
                claim
                    .validity
                    .axis_quantities()
                    .get(&request.stress_ratio_axis)
                    == Some(&dimensionless)
                    && claim.validity.bounds().get(&request.stress_ratio_axis)
                        == Some(&(-1.0, -1.0))
                    && claim.validity.axis_quantities().get(life_axis) == Some(&dimensionless)
                    && claim.validity.bounds().contains_key(life_axis),
                "all fitted coefficients need source-declared R = -1 and typed finite life support",
            )?;
        } else {
            require(
                !claim.validity.bounds().contains_key(life_axis)
                    && !claim
                        .validity
                        .bounds()
                        .contains_key(&request.stress_ratio_axis),
                "the elastic modulus must be a physical-state property, not a fatigue-life fit",
            )?;
        }
        properties.push((claim.clone(), answer));
    }
    let fit = &properties[CmbProperty::FatigueStrengthCoefficient.index()].0;
    let source_domain = fit.validity.bounds()[life_axis];
    for (claim, _) in properties.iter().skip(2) {
        require(
            claim.provenance == fit.provenance
                && claim.observations == fit.observations
                && claim.validity.bounds()[life_axis] == source_domain,
            "fatigue coefficients must retain one source fit, population and life interval",
        )?;
    }
    let value = |p: CmbProperty| properties[p.index()].1.evidence.value.value;
    let b = value(CmbProperty::FatigueStrengthExponent);
    let c = value(CmbProperty::FatigueDuctilityExponent);
    let convert = |coefficient: f64, exponent: f64| {
        if request.power_basis == CmbPowerBasis::Reversals || coefficient == 0.0 {
            coefficient
        } else {
            coefficient * det::exp(-exponent * core::f64::consts::LN_2)
        }
    };
    let parameters = CmbParameters {
        youngs_modulus_pa: value(CmbProperty::YoungsModulus),
        fatigue_strength_coefficient_pa: convert(value(CmbProperty::FatigueStrengthCoefficient), b),
        fatigue_strength_exponent: b,
        fatigue_ductility_coefficient: convert(value(CmbProperty::FatigueDuctilityCoefficient), c),
        fatigue_ductility_exponent: c,
    };
    let life_scale = if request.power_basis == CmbPowerBasis::Reversals {
        0.5
    } else {
        1.0
    };
    let id = card.id();
    let scope = FatigueScope {
        source_reference: fit.provenance.source.clone(),
        material: id.chemistry.clone(),
        process: id.process.clone(),
        surface_and_geometry: request.surface_and_geometry.clone(),
        environment: request.environment.clone(),
        load_basis: request.load_basis.clone(),
        temperature_k,
    };
    let mut curve = StrainLifeCurve::new(
        scope,
        parameters,
        source_domain.0 * life_scale,
        source_domain.1 * life_scale,
    )?;
    curve.source = Some(Arc::new(CmbCalibrationSource {
        card_identity: card.content_hash(),
        material_state: id.clone(),
        request: request.clone(),
        query: point.clone(),
        selections,
        properties,
    }));
    Ok(curve)
}
