//! Shared material-coordinate grammar for file-driven physics commands.
//! Units and semantic kinds survive admission into the canonical QueryPoint;
//! neither a command name nor a material name selects their interpretation.

pub(crate) use crate::json_read::JsonValue;

use fs_matdb::QueryPoint;
use fs_qty::parse::parse_qty;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec};

/// Explicit legacy, dimensional or semantic coordinate convention. Equal
/// dimensions alone do not equate pressure with stress or energy with torque.
pub(crate) fn quantity(name: &str, dims: Dims) -> Result<Option<QuantitySpec>, String> {
    if name == "legacy" { return Ok(None); }
    if name == "dimensional" { return Ok(Some(QuantitySpec::dimensional(dims))); }
    let kind = match name {
        "absolute-temperature" => QuantityKind::AbsoluteTemperature,
        "temperature-difference" => QuantityKind::TemperatureDifference,
        "pressure" => QuantityKind::Pressure,
        "stress" => QuantityKind::Stress,
        "energy" => QuantityKind::Energy,
        "torque" => QuantityKind::Torque,
        "entropy" => QuantityKind::Entropy,
        "heat-capacity" => QuantityKind::HeatCapacity,
        name => QuantityKind::from_material_convention_name(name)
            .ok_or_else(|| format!("unsupported quantity kind {name:?}"))?,
    };
    if kind.expected_dims() != dims { return Err("quantity kind and unit dimensions disagree".into()); }
    Ok(Some(QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))))
}

/// The same bounded fixed-coordinate input for bulk, interface and fluid
/// source queries. Commands retain their complete surrounding JSON document.
pub(crate) fn state_point(definition: &JsonValue) -> Result<QueryPoint, String> {
    let axes = definition.get("state").ok_or("state is required")?
        .as_array().ok_or("state must be an array")?;
    if axes.len() > 32 { return Err("state exceeds 32 entries".into()); }
    let mut point = QueryPoint::new();
    for axis in axes {
        let object = axis.as_object().ok_or("expected an object")?;
        if let Some((key, _)) = object.iter().find(|(key, _)|
            !["name", "value", "kind"].contains(&key.as_str()))
        {
            return Err(format!("unsupported field {key:?}"));
        }
        let text = |key: &str| -> Result<&str, String> {
            axis.get(key).ok_or_else(|| format!("{key} is required"))?
                .as_str().filter(|s| !s.trim().is_empty())
                .ok_or_else(|| format!("{key} must be a nonempty string"))
        };
        let name = text("name")?;
        if point.axes().contains_key(name) { return Err(format!("duplicate state axis {name}")); }
        let v = parse_qty(text("value")?).map_err(|e| e.to_string())?;
        point = match quantity(text("kind")?, v.dims)? {
            Some(q) => point.with_quantity(name, q, v.value),
            None => point.with(name, v.value),
        }.map_err(|e| e.to_string())?;
    }
    Ok(point)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn physical_coordinates_keep_units_and_distinct_semantic_kinds() {
        let input = JsonValue::parse(r#"{"state":[
            {"name":"temperature","value":"300 K","kind":"absolute-temperature"},
            {"name":"pressure","value":"100 kPa","kind":"pressure"},
            {"name":"stress","value":"100 kPa","kind":"stress"}
        ]}"#).unwrap();
        let point = state_point(&input).unwrap();
        assert_eq!(point.axes()["pressure"], 100_000.0);
        assert_eq!(point.axes()["pressure"], point.axes()["stress"]);
        assert_ne!(point.axis_quantities()["pressure"], point.axis_quantities()["stress"]);
        assert_eq!(point.axis_quantities()["temperature"], quantity("absolute-temperature", fs_qty::Temperature::DIMS).unwrap().unwrap());
    }

    #[test]
    fn malformed_and_ambiguous_material_coordinates_refuse() {
        for input in [
            r#"{"state":[{"name":"T","value":"300 K","kind":"pressure"}]}"#,
            r#"{"state":[{"name":"T","value":"300 K","kind":"legacy"},{"name":"T","value":"310 K","kind":"legacy"}]}"#,
            r#"{"state":[{"name":"T","value":"300 K","kind":"legacy","ignored":true}]}"#,
            r#"{"state":[{"name":"T","value":"300 K"}]}"#,
        ] {
            assert!(state_point(&JsonValue::parse(input).unwrap()).is_err(), "{input}");
        }
    }

    #[test]
    fn legacy_and_dimensional_coordinates_remain_explicitly_distinct() {
        let input = JsonValue::parse(r#"{"state":[
            {"name":"T","value":"300 K","kind":"legacy"},
            {"name":"temperature","value":"300 K","kind":"dimensional"}
        ]}"#).unwrap();
        let point = state_point(&input).unwrap();
        assert!(!point.axis_quantities().contains_key("T"));
        assert_eq!(point.axis_quantities()["temperature"], QuantitySpec::dimensional(fs_qty::Temperature::DIMS));
    }
}
