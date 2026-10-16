//! Synthetic source data, not measured properties of a named material.
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_qty::{Density, Dims, Temperature};

pub fn density_card(label: &str, low: f64, high: f64) -> MaterialCard {
    card(label, "density", Density::DIMS, &[(low, high)])
}
pub fn card(label: &str, key: &str, dims: Dims, values: &[(f64, f64)]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (i, &(low, high)) in values.iter().enumerate() {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(key, dims),
            value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(200.0, low), (400.0, high)], dims },
            validity: ValidityDomain::unconstrained().with("T", 200.0, 400.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside, observations: vec![],
            provenance: Provenance { source: format!("synthetic density source {label}-{i}"),
                license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId { chemistry: label.into(), phase: "solid".into(),
        process: "synthetic fixture".into(), revision: 0 }, claims, vec![]).unwrap()
}
pub fn point(t: f64) -> QueryPoint { QueryPoint::new().with("T", t).unwrap() }
pub fn close(a: f64, b: f64) {
    assert!((a-b).abs() <= 2e-12*a.abs().max(b.abs()).max(1.0), "{a} != {b}");
}
