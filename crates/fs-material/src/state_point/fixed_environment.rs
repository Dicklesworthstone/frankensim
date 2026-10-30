//! Shared admission of held environmental coordinates on original material
//! and ordered-interface claims. This builds a physical query; it does not
//! evolve transport or replace source-domain/property resolution.

use core::fmt;
use fs_matdb::{ClaimSet, MatDbError, MaterialCard, QueryPoint};
use fs_qty::{QuantitySpec, semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm}};

/// Maximum supplied coordinates in one bounded base or held-environment stage.
pub const MAX_FIXED_ENVIRONMENT_COORDINATES: usize = 8;

/// A held coordinate cannot be represented with the original source semantics.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialEnvironmentError {
    /// The caller exceeded the explicit coordinate budget.
    TooManyCoordinates {
        /// Offered number of additional coordinates.
        observed: usize,
        /// Admitted maximum.
        maximum: usize,
    },
    /// Missing, conflicting or unsupported physical coordinate meaning.
    Coordinate {
        /// Original source axis.
        axis: String,
        /// Stable explanation of the refusal.
        reason: &'static str,
    },
    /// The original quantity/query owner refused the coordinate value.
    Query {
        /// Original source axis.
        axis: String,
        /// Original refusal, including semantic value-domain checks.
        source: MatDbError,
    },
}
impl fmt::Display for MaterialEnvironmentError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TooManyCoordinates { observed, maximum } =>
                write!(f, "material query has {observed} supplied coordinates; maximum {maximum}"),
            Self::Coordinate { axis, reason } => write!(f, "material coordinate {axis}: {reason}"),
            Self::Query { axis, source } => write!(f, "material coordinate {axis}: {source}"),
        }
    }
}
impl std::error::Error for MaterialEnvironmentError {}

fn coordinate(axis: &str, reason: &'static str) -> MaterialEnvironmentError {
    MaterialEnvironmentError::Coordinate { axis: axis.into(), reason }
}

/// Add explicitly held environmental coordinates without changing an existing
/// physical query. The caller retains ownership of base coordinates such as
/// temperature, normal pressure, slip speed or angular frequency.
///
/// Every claim declaring an extra axis must supply the same explicit static
/// pressure, composition, moisture or concentration meaning. A dry-mass
/// moisture fraction is never inferred from atmospheric relative humidity.
/// Source bounds and property selection remain the subsequent resolver's job.
///
/// # Errors
/// Refuses too many coordinates, base-coordinate overwrites, duplicates,
/// nonfinite values, missing/untyped/contradictory source axes, dynamic forms,
/// or quantities outside the supported external-state roles.
pub fn append_held_environment(
    claims: &ClaimSet,
    mut point: QueryPoint,
    environment: &[(String, f64)],
) -> Result<QueryPoint, MaterialEnvironmentError> {
    if environment.len() > MAX_FIXED_ENVIRONMENT_COORDINATES {
        return Err(MaterialEnvironmentError::TooManyCoordinates {
            observed: environment.len(), maximum: MAX_FIXED_ENVIRONMENT_COORDINATES,
        });
    }
    for (axis, value) in environment {
        if axis.is_empty() || axis.len() > 256 || axis.chars().any(char::is_control)
            || axis == "T" || axis == "omega" || point.axes().contains_key(axis)
            || !value.is_finite()
        {
            return Err(coordinate(axis,
                "held coordinates must be distinct finite values outside reserved or existing physical axes"));
        }
        let mut quantity = None;
        for (_, claim) in claims.claims_ordered() {
            if claim.validity.bound(axis).is_none() { continue; }
            let declared = claim.validity.axis_quantities().get(axis)
                .copied().ok_or_else(|| coordinate(axis,
                    "held coordinates require explicit source quantity semantics"))?;
            if quantity.is_some_and(|prior| prior != declared) {
                return Err(coordinate(axis, "source claims disagree on coordinate semantics"));
            }
            quantity = Some(declared);
        }
        let quantity = quantity.ok_or_else(|| coordinate(axis,
            "held coordinate is absent from the original source"))?;
        if !quantity.semantic_type().is_some_and(|semantic|
            semantic.form() == ValueForm::Static && matches!(semantic.kind(),
                QuantityKind::Pressure | QuantityKind::Composition(_) | QuantityKind::Moisture(_)
                    | QuantityKind::MassConcentration | QuantityKind::AmountConcentration))
        {
            return Err(coordinate(axis,
                "held coordinates require static pressure, composition, moisture or concentration"));
        }
        point = point.with_quantity(axis.clone(), quantity, *value)
            .map_err(|source| MaterialEnvironmentError::Query { axis: axis.clone(), source })?;
    }
    Ok(point)
}

/// Build a bounded base query in explicitly owned physical conventions.
/// Every source declaration must equal the expected convention or its legacy
/// dimension-only form. Bare axes remain bare; existing source descriptors are
/// retained. Incompatible aliases (stress/pressure, delta/absolute temperature,
/// cyclic/angular frequency, or waveform forms) never undergo implicit conversion.
///
/// # Errors
/// Refuses repeated or malformed axes, excess coordinates, conflicting source
/// declarations, and any value refused by the original quantity/query owner.
pub fn query_with_source_conventions(
    claims: &ClaimSet,
    coordinates: &[(&str, QuantitySpec, f64)],
) -> Result<QueryPoint, MaterialEnvironmentError> {
    if coordinates.len() > MAX_FIXED_ENVIRONMENT_COORDINATES {
        return Err(MaterialEnvironmentError::TooManyCoordinates {
            observed: coordinates.len(), maximum: MAX_FIXED_ENVIRONMENT_COORDINATES,
        });
    }
    let mut point = QueryPoint::new();
    for &(axis, expected, value) in coordinates {
        if axis.is_empty() || axis.len() > 256 || axis.chars().any(char::is_control)
            || point.axes().contains_key(axis)
        {
            return Err(coordinate(axis, "owned query axes must be distinct bounded source names"));
        }
        let mut quantity = None;
        for (_, claim) in claims.claims_ordered() {
            if let Some(&declared) = claim.validity.axis_quantities().get(axis) {
                if declared != expected && declared != QuantitySpec::dimensional(expected.dims()) {
                    return Err(coordinate(axis, "source convention disagrees with the physical coordinate owner"));
                }
                if quantity.is_none() { quantity = Some(declared); }
            }
        }
        point = match quantity {
            Some(quantity) => point.with_quantity(axis, quantity, value),
            None => point.with(axis, value),
        }.map_err(|source| MaterialEnvironmentError::Query { axis: axis.into(), source })?;
    }
    Ok(point)
}

/// Build the shared kelvin/radians-per-second material query and append the
/// explicitly held original environment. This is the convention used by the
/// native string, plate, reed and duct-wall structural consumers.
///
/// Existing bare or dimension-only T/omega declarations retain their query
/// convention. Explicit semantic declarations must mean absolute temperature
/// and angular frequency. No material property or missing state is inferred.
///
/// # Errors
/// Refuses incompatible temperature/frequency semantics or any refusal from
/// [`append_held_environment`] and the original query-value owner.
pub fn query_at_fixed_environment(
    card: &MaterialCard,
    temperature: f64,
    omega: f64,
    environment: &[(String, f64)],
) -> Result<QueryPoint, MaterialEnvironmentError> {
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature, ValueForm::Static));
    let angular = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Frequency(FrequencyConvention::Angular), ValueForm::Static));
    let point = query_with_source_conventions(card.claims(), &[
        ("T", absolute, temperature), ("omega", angular, omega),
    ])?;
    append_held_environment(card.claims(), point, environment)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
        Provenance, UncertaintyModel};
    use fs_qty::{Dims, semantic::MoistureBasis};

    fn source(name: &str, axis: &str, quantity: Option<QuantitySpec>) -> PropertyClaim {
        let (low, high) = match axis { "T" => (250.0, 350.0), "p" => (0.0, 20000.0), _ => (0.0, 1.0) };
        let validity = match quantity {
            Some(quantity) => ValidityDomain::unconstrained().with_quantity(axis, quantity, low, high),
            None => ValidityDomain::unconstrained().with(axis, low, high),
        };
        PropertyClaim {
            key: PropertyKey::new(name, Dims::NONE),
            value: PropertyValue::Scalar { value: 0.3, dims: Dims::NONE },
            validity, interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated, observations: vec![],
            provenance: Provenance { source: "manufactured environmental-coordinate fixture".into(),
                license: "CC0-1.0".into(), artifact: None },
        }
    }
    fn moisture(basis: MoistureBasis) -> QuantitySpec {
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Moisture(basis), ValueForm::Static))
    }
    fn card(claims: ClaimSet) -> MaterialCard {
        MaterialCard::assemble(MaterialStateId { chemistry: "fixture".into(), phase: "solid".into(),
            process: "held environment".into(), revision: 0 }, claims, vec![]).unwrap()
    }

    #[test]
    fn held_material_and_interface_context_retains_original_coordinates() {
        let q = moisture(MoistureBasis::DryMass);
        let mut claims = ClaimSet::new();
        claims.insert_claim(source("friction", "moisture_content", Some(q))).unwrap();
        let material = card(claims.clone());
        let extras = [("moisture_content".into(), 0.12)];
        let query = query_at_fixed_environment(&material, 293.15, 1000.0, &extras).unwrap();
        assert_eq!(query.axes()["T"], 293.15);
        assert_eq!(query.axes()["omega"], 1000.0);
        assert_eq!(query.axes()["moisture_content"], 0.12);
        assert_eq!(query.axis_quantities()["moisture_content"], q);
        assert!(!query.axis_quantities().contains_key("T"), "legacy bare temperature stays bare");
        let base = QueryPoint::new().with("T", 293.15).unwrap()
            .with("p", 12000.0).unwrap().with("v", 0.2).unwrap();
        let interface = append_held_environment(&claims, base.clone(), &extras).unwrap();
        for (axis, value) in base.axes() { assert_eq!(interface.axes()[axis], *value); }
        assert_eq!(interface.axis_quantities()["moisture_content"], q);
        assert_eq!(append_held_environment(&claims, base.clone(), &[]).unwrap(), base);
        for axis in ["T", "omega", "p", "v"] {
            assert!(append_held_environment(&claims, base.clone(), &[(axis.into(), 0.1)]).is_err());
        }
        assert!(append_held_environment(&claims, base.clone(),
            &[extras[0].clone(), extras[0].clone()]).is_err());
        assert!(append_held_environment(&claims, base.clone(), &[("absent".into(), 0.1)]).is_err());
        assert!(append_held_environment(&claims, base, &[("moisture_content".into(), f64::NAN)]).is_err());
    }

    #[test]
    fn held_context_requires_every_source_to_agree_on_exact_physical_meaning() {
        let q = moisture(MoistureBasis::DryMass);
        for other in [None, Some(moisture(MoistureBasis::WetMass)),
            Some(moisture(MoistureBasis::RelativeHumidity))] {
            let mut claims = ClaimSet::new();
            claims.insert_claim(source("first", "moisture", Some(q))).unwrap();
            claims.insert_claim(source("second", "moisture", other)).unwrap();
            assert!(append_held_environment(&claims, QueryPoint::new(),
                &[("moisture".into(), 0.12)]).is_err());
        }
        let mut claims = ClaimSet::new();
        claims.insert_claim(source("untyped", "moisture", None)).unwrap();
        assert!(append_held_environment(&claims, QueryPoint::new(),
            &[("moisture".into(), 0.12)]).is_err());
        let cyclic = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Frequency(FrequencyConvention::Cyclic), ValueForm::Static));
        let mut claims = ClaimSet::new();
        claims.insert_claim(source("cyclic", "omega", Some(cyclic))).unwrap();
        assert!(query_at_fixed_environment(&card(claims), 293.15, 1.0, &[]).is_err());
    }

    #[test]
    fn owned_source_conventions_retain_typed_contact_axes_and_refuse_aliases() {
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature, ValueForm::Static));
        let pressure = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Pressure, ValueForm::Static));
        let speed = QuantitySpec::dimensional(Dims([1, 0, -1, 0, 0, 0]));
        let coordinates = [("T", absolute, 293.15), ("p", pressure, 12000.0), ("v", speed, 0.2)];
        let mut claims = ClaimSet::new();
        for (i, &(axis, quantity, _)) in coordinates.iter().enumerate() {
            claims.insert_claim(source(&format!("base-{i}"), axis, Some(quantity))).unwrap();
        }
        let point = query_with_source_conventions(&claims, &coordinates).unwrap();
        for (axis, quantity, value) in coordinates {
            assert_eq!(point.axes()[axis], value);
            assert_eq!(point.axis_quantities()[axis], quantity);
        }
        let stress = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Stress, ValueForm::Static));
        claims.insert_claim(source("different-pressure-meaning", "p", Some(stress))).unwrap();
        assert!(query_with_source_conventions(&claims, &coordinates).is_err());
        assert!(query_with_source_conventions(&ClaimSet::new(),
            &[coordinates[0], coordinates[0]]).is_err());
        assert!(query_with_source_conventions(&ClaimSet::new(),
            &[("T", absolute, f64::NAN)]).is_err());
        assert!(matches!(append_held_environment(&ClaimSet::new(), QueryPoint::new(),
            &vec![("moisture".into(), 0.1); 9]), Err(MaterialEnvironmentError::TooManyCoordinates { .. })));
    }
}
