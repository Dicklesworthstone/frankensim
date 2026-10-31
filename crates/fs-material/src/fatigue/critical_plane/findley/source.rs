//! Original material-card Findley calibration, with no inferred fatigue data.
//!
//! The two distinct properties are the already calibrated normal-stress
//! sensitivity k and Findley limit factor f. Neither ordinary strength nor a
//! single axial S-N datum supplies this pair. The original fit, observations,
//! uncertainty, source interpolation and selected claims remain inspectable.
//! This owner does not fit parameters, manufacture an endurance limit or turn
//! a limit-factor comparison into finite life or accumulated Miner damage.

use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_matdb::{
    ClaimSelection, MatDbError, MaterialAnswer, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, QueryPoint,
};
use fs_qty::semantic::{QuantityKind, ValueForm};
use fs_qty::{Dims, Pressure, Temperature};

use super::{FindleyCalibration, FindleyError};
use crate::fatigue::FatigueScope;

/// Source meanings are separate even though k shares the dimensions of many
/// unrelated ratios and f shares the dimensions of ordinary static strength.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FindleyProperty {
    NormalStressSensitivity,
    LimitFactor,
}

impl FindleyProperty {
    pub const ALL: [Self; 2] = [Self::NormalStressSensitivity, Self::LimitFactor];

    pub const fn index(self) -> usize {
        match self {
            Self::NormalStressSensitivity => 0,
            Self::LimitFactor => 1,
        }
    }

    pub fn property_key(self) -> PropertyKey {
        match self {
            Self::NormalStressSensitivity => {
                PropertyKey::new("findley-normal-stress-sensitivity", Dims::NONE)
            }
            Self::LimitFactor => PropertyKey::new("findley-limit-factor", Pressure::DIMS),
        }
    }
}

/// Fixed external context for the original fitted k/f pair. Temperature is
/// explicit semantic AbsoluteTemperature/Static in SI kelvin. These labels
/// declare applicability; they do not independently validate its transfer to
/// another specimen, process, load spectrum or component.
#[derive(Debug, Clone, PartialEq)]
pub struct FindleySourceRequest {
    pub temperature_axis: String,
    pub surface_and_geometry: String,
    pub environment: String,
    pub load_basis: String,
}

/// Kept inside the immutable calibration and shared by native cycle gauges.
#[derive(Debug, Clone, PartialEq)]
pub struct FindleyCalibrationSource {
    card_identity: ContentHash,
    material_state: MaterialStateId,
    request: FindleySourceRequest,
    query: QueryPoint,
    selections: [ClaimSelection; 2],
    properties: Vec<(PropertyClaim, MaterialAnswer)>,
}

impl FindleyCalibrationSource {
    pub const fn card_identity(&self) -> ContentHash {
        self.card_identity
    }
    pub fn material_state(&self) -> &MaterialStateId {
        &self.material_state
    }
    pub fn request(&self) -> &FindleySourceRequest {
        &self.request
    }
    pub fn query(&self) -> &QueryPoint {
        &self.query
    }
    pub fn selections(&self) -> &[ClaimSelection; 2] {
        &self.selections
    }
    pub fn property(&self, property: FindleyProperty) -> &MaterialAnswer {
        &self.properties[property.index()].1
    }
    pub fn claim(&self, property: FindleyProperty) -> &PropertyClaim {
        &self.properties[property.index()].0
    }

    /// Match the ORIGINAL native material bundle's card, named state and every
    /// source-query coordinate and quantity convention. Additional native axes
    /// do not extend a calibration; the original claims and their domains are
    /// retained. No source switch or temperature/environment transfer occurs.
    pub fn matches_material_state(
        &self,
        card_identity: ContentHash,
        material: &MaterialStateId,
        point: &QueryPoint,
    ) -> bool {
        self.card_identity == card_identity
            && &self.material_state == material
            && self.query.axes().iter().all(|(axis, value)| {
                point.axes().get(axis) == Some(value)
                    && point.axis_quantities().get(axis) == self.query.axis_quantities().get(axis)
            })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum FindleySourceError {
    Source(MatDbError),
    Model(FindleyError),
    Invalid(&'static str),
}

impl core::fmt::Display for FindleySourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound Findley calibration refused: {self:?}")
    }
}
impl std::error::Error for FindleySourceError {}
impl From<MatDbError> for FindleySourceError {
    fn from(error: MatDbError) -> Self {
        Self::Source(error)
    }
}
impl From<FindleyError> for FindleySourceError {
    fn from(error: FindleyError) -> Self {
        Self::Model(error)
    }
}

fn require(condition: bool, reason: &'static str) -> Result<(), FindleySourceError> {
    if condition {
        Ok(())
    } else {
        Err(FindleySourceError::Invalid(reason))
    }
}

/// Resolve the ORIGINAL fitted sensitivity and limit factor. Selections follow
/// `FindleyProperty::ALL`; explicit individual pins resolve ambiguous sources.
/// Both properties must belong to the same provenance/observation fit and each
/// must admit the actual typed temperature. Fixed scalars and source-supported
/// temperature curves are admitted; other parameter dependencies refuse.
///
/// The source's own query/interpolation policy is used without inventing a fit
/// or discarding source knots. Existing source validity bounds intersect at the
/// query; their original uncertainty and domains remain attached unchanged.
pub fn resolve_findley(
    card: &MaterialCard,
    point: &QueryPoint,
    request: &FindleySourceRequest,
    selections: [ClaimSelection; 2],
) -> Result<FindleyCalibration, FindleySourceError> {
    require(
        !request.temperature_axis.trim().is_empty(),
        "declare the source temperature axis",
    )?;
    let temperature_quantity = point
        .axis_quantities()
        .get(&request.temperature_axis)
        .ok_or(FindleySourceError::Invalid(
            "typed source temperature is missing",
        ))?;
    require(
        temperature_quantity.semantic_type().is_some_and(|q| {
            q.kind() == QuantityKind::AbsoluteTemperature && q.form() == ValueForm::Static
        }),
        "source temperature must be semantic static absolute temperature in kelvin",
    )?;
    let temperature_k =
        *point
            .axes()
            .get(&request.temperature_axis)
            .ok_or(FindleySourceError::Invalid(
                "source temperature value is missing",
            ))?;
    require(
        temperature_k.is_finite() && temperature_k > 0.0,
        "source temperature must be finite positive kelvin",
    )?;

    let mut properties = Vec::with_capacity(2);
    for property in FindleyProperty::ALL {
        let key = property.property_key();
        let answer = match selections[property.index()] {
            ClaimSelection::Policy(policy) => card.claims().query_typed(&key, point, policy)?,
            ClaimSelection::Pinned(pin) => card.claims().query_pinned_typed(&key, point, pin)?,
        };
        card.claims().verify_receipt(&answer.receipt)?;
        let claim =
            card.claims()
                .claim(answer.receipt.selected)
                .ok_or(FindleySourceError::Invalid(
                    "selected original Findley claim is absent",
                ))?;
        require(
            claim
                .validity
                .axis_quantities()
                .get(&request.temperature_axis)
                == Some(temperature_quantity)
                && claim
                    .validity
                    .bounds()
                    .contains_key(&request.temperature_axis),
            "both Findley parameters need an explicit matching typed source temperature domain",
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
            "Findley parameters must be fixed scalars or functions of the declared temperature",
        )?;
        properties.push((claim.clone(), answer));
    }

    let first = &properties[0].0;
    let second = &properties[1].0;
    require(
        first.provenance == second.provenance && first.observations == second.observations,
        "Findley sensitivity and limit must retain one source fit and observation population",
    )?;
    let scope = FatigueScope {
        source_reference: first.provenance.source.clone(),
        material: card.id().chemistry.clone(),
        process: card.id().process.clone(),
        surface_and_geometry: request.surface_and_geometry.clone(),
        environment: request.environment.clone(),
        load_basis: request.load_basis.clone(),
        temperature_k,
    };
    let mut calibration = FindleyCalibration::new(
        scope,
        properties[0].1.evidence.value.value,
        properties[1].1.evidence.value.value,
    )?;
    calibration.source = Some(Arc::new(FindleyCalibrationSource {
        card_identity: card.content_hash(),
        material_state: card.id().clone(),
        request: request.clone(),
        query: point.clone(),
        selections,
        properties,
    }));
    Ok(calibration)
}
