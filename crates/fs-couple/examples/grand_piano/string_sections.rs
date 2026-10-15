//! Data-driven bonded concentric speaking-string materials.
//!
//! One or more layer records per course, ordered from centre outward:
//! `layer,key,outer_radius_m,density_kg_m3,young_pa,source,license`.
//! The source and license are explicit input claims, NOT measurements verified
//! by this parser. No preset or representative metal is silently selected.
//! Helically wound strings must use their own geometry/slip model instead.
use std::collections::{BTreeMap, BTreeSet};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY,
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point};
use fs_qty::{Density, Pressure};

use super::geometry::Course;

pub const HEADER: &str = "frankensim-string-sections-v1";
pub const ROW: &str = "layer,key,outer_radius_m,density_kg_m3,young_pa,source,license";

fn layer(key: u8, ordinal: usize, radius: f64, density: f64, young: f64,
    source: &str, license: &str) -> Result<ConcentricMaterialLayer, String> {
    if source.is_empty() || license.is_empty() || source.contains(',')
        || license.contains(',') || !density.is_finite() || density <= 0.0
        || !young.is_finite() || young <= 0.0 {
        return Err("layer needs positive finite properties, source and license".into());
    }
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (DENSITY_PROPERTY, Density::DIMS, density),
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, young),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: source.into(), license: license.into(), artifact: None,
            },
        }).map_err(|e| format!("invalid property claim: {e}"))?;
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: format!("declared-string-{key}-layer-{ordinal}"),
        phase: "solid".into(), process: "authored-source-claim".into(), revision: 0,
    }, claims, vec![]).map_err(|e| format!("invalid material card: {e}"))?;
    let material = resolve_material_state_point(&card, &QueryPoint::new(), &[
        ScalarPropertyRequirement::try_new(DENSITY_PROPERTY, Density::DIMS,
            ScalarAdmissibility::StrictlyPositive).map_err(|e| e.to_string())?,
        ScalarPropertyRequirement::try_new(YOUNG_MODULUS_PROPERTY, Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive).map_err(|e| e.to_string())?,
    ], MaterialPropertySelection::SingleClaimOnly).map_err(|e| e.to_string())?;
    Ok(ConcentricMaterialLayer { outer_radius_m: radius, material })
}

/// Replace physical string coefficients for every admitted key, all-or-nothing.
/// Does not change tension or the equations in the string/bridge/contact solver.
/// The source model is perfectly bonded concentric layers; winding is NOT bonded.
///
/// # Errors
/// Missing keys, repeated or out-of-order layers, invalid material claims,
/// hollow/non-solid geometry, and unrepresentable derived dynamics are refused.
pub fn apply(text: &str, courses: &[Course]) -> Result<Vec<Course>, String> {
    if courses.is_empty() { return Err("string material override needs a nonempty scale".into()); }
    let keys: BTreeSet<u8> = courses.iter().map(|c| c.midi).collect();
    if keys.len() != courses.len() { return Err("string scale has duplicate keys".into()); }
    let mut rows: BTreeMap<u8, Vec<ConcentricMaterialLayer>> = BTreeMap::new();
    let mut header_seen = false;
    let mut column_seen = false;
    for (index, raw) in text.lines().enumerate() {
        let row = raw.split('#').next().unwrap_or("").trim();
        if row.is_empty() { continue; }
        let error = |why: &str| format!("string sections line {}: {why}", index + 1);
        if !header_seen {
            if row != HEADER { return Err(error("expected string-sections header")); }
            header_seen = true;
            continue;
        }
        if !column_seen {
            if row != ROW { return Err(error("expected explicit SI column header")); }
            column_seen = true;
            continue;
        }
        let fields: Vec<_> = row.split(',').map(str::trim).collect();
        if fields.len() != 7 || fields[0] != "layer" { return Err(error("expected seven layer fields")); }
        let key: u8 = fields[1].parse().map_err(|_| error("invalid key"))?;
        if !keys.contains(&key) { return Err(error("key absent from the selected scale")); }
        let radius: f64 = fields[2].parse().map_err(|_| error("invalid outer radius"))?;
        let rho: f64 = fields[3].parse().map_err(|_| error("invalid density"))?;
        let young: f64 = fields[4].parse().map_err(|_| error("invalid Young modulus"))?;
        let layers = rows.entry(key).or_default();
        if layers.len() >= fs_material::section::MAX_CONCENTRIC_LAYERS {
            return Err(error("too many layers"));
        }
        if !radius.is_finite() || radius <= 0.0 ||
            layers.last().is_some_and(|last| radius <= last.outer_radius_m) {
            return Err(error("solid concentric layers need strictly increasing positive radii"));
        }
        layers.push(layer(key, layers.len(), radius, rho, young, fields[5], fields[6])
            .map_err(|e| error(&e))?);
    }
    if !header_seen || !column_seen { return Err("missing string section header or columns".into()); }
    courses.iter().map(|course| {
        let layers = rows.remove(&course.midi)
            .ok_or_else(|| format!("missing bonded section for key {}; no fallback", course.midi))?;
        let section = ConcentricElasticSection::try_new(0.0, layers)
            .map_err(|e| format!("key {}: {e}", course.midi))?;
        course.with_bonded_material_section(&section)
            .map_err(|e| format!("key {}: {e}", course.midi))
    }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn material_change_reaches_string_spectrum_without_tension_retuning() {
        let source = [super::super::geometry::demonstration_scale().unwrap()[48]];
        let row = |rho: f64, young: f64| format!(
            "{HEADER}\n{ROW}\nlayer,69,0.0005,{rho},{young},synthetic,CC0-1.0\n");
        let light = apply(&row(7800.0, 200e9), &source).unwrap();
        let heavy = apply(&row(10500.0, 80e9), &source).unwrap();
        assert_eq!(light[0].tension_n, heavy[0].tension_n);
        assert!(heavy[0].linear_density_kg_m > light[0].linear_density_kg_m);
        assert!(heavy[0].flexural_rigidity_nm2 < light[0].flexural_rigidity_nm2);
        assert!(heavy[0].partial_hz(1, heavy[0].tension_n) <
            light[0].partial_hz(1, light[0].tension_n));
        assert_eq!(light[0].midi, source[0].midi);
        assert!(apply(&row(7800.0, 200e9), &[]).is_err());
        assert!(apply(&row(-1.0, 200e9), &source).is_err());
        assert!(apply(&row(7800.0, f64::NAN), &source).is_err());
        assert!(apply(&row(7800.0, 200e9).replace("synthetic", ""), &source).is_err());
        assert!(apply(&format!("{}layer,69,0.0004,7800,200000000000,synthetic,CC0\n",
            row(7800.0, 200e9)), &source).is_err());
    }
}
