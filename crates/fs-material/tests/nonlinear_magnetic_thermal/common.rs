#![allow(dead_code)]
//! Manufactured constitutive data only; no measured magnetic material claims.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::ElectricalPropertyBasis,
    field_storage::{
        MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY,
        nonlinear::{
            REVERSIBLE_FIELD_STRENGTH_PROPERTY, ReversibleMagneticMaterial,
            thermal::ThermalReversibleMagneticMaterial,
        },
        permanent::{FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS},
    },
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, SPECIFIC_HEAT_CAPACITY_DIMS,
        ScalarAdmissibility, ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS,
        resolve_material_state_point,
    },
};
use fs_qty::{
    Density, Dims, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

pub const SCALE: &str = "reversible_magnetic_reference_scale";
pub fn absolute() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
pub fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", absolute(), t)
        .unwrap()
        .with_quantity("B", QuantitySpec::dimensional(FLUX_DENSITY_DIMS), 0.0)
        .unwrap()
}
pub fn scale_key() -> PropertyKey {
    PropertyKey::new(SCALE, Dims([0; 6]))
}
pub fn card_with_scale(knots: Vec<(f64, f64)>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let domain = ValidityDomain::unconstrained().with_quantity("T", absolute(), 200.0, 500.0);
    let mut insert = |name: &str, dims: Dims, value: PropertyValue, validity| {
        let interpolation = if matches!(value, PropertyValue::Scalar { .. }) {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value,
                validity,
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured separable magnetic thermal fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    insert(
        REVERSIBLE_FIELD_STRENGTH_PROPERTY,
        FIELD_STRENGTH_DIMS,
        PropertyValue::Curve {
            abscissa: "B".into(),
            abscissa_dims: FLUX_DENSITY_DIMS,
            knots: vec![(0.0, 0.0), (1.0, 1.0), (4.0, 13.0)],
            dims: FIELD_STRENGTH_DIMS,
        },
        domain
            .clone()
            .with_quantity("B", QuantitySpec::dimensional(FLUX_DENSITY_DIMS), 0.0, 4.0),
    );
    insert(
        SCALE,
        Dims([0; 6]),
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots,
            dims: Dims([0; 6]),
        },
        domain.clone(),
    );
    let basis = ElectricalPropertyBasis::Resistivity;
    insert(
        basis.property_name(),
        basis.dimensions(),
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![
                (200.0, 0.5 / 18.0),
                (300.0, 1.0 / 18.0),
                (500.0, 2.0 / 18.0),
            ],
            dims: basis.dimensions(),
        },
        domain.clone(),
    );
    for (name, dims, value) in [
        ("density", Density::DIMS, 2.0),
        ("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS, 3.0),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 1.0),
        (
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            1.0,
        ),
    ] {
        insert(
            name,
            dims,
            PropertyValue::Scalar { value, dims },
            domain.clone(),
        );
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic mixed magnetic conductor".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
pub fn card(slope: f64) -> MaterialCard {
    card_with_scale(vec![
        (200.0, 1.0 - 100.0 * slope),
        (300.0, 1.0),
        (500.0, 1.0 + 200.0 * slope),
    ])
}
pub fn bundle(card: &MaterialCard, t: f64) -> ResolvedMaterialStatePoint {
    let positive = ScalarAdmissibility::StrictlyPositive;
    let mut requirements = vec![
        ElectricalPropertyBasis::Resistivity.requirement().unwrap(),
        ScalarPropertyRequirement::try_new(
            REVERSIBLE_FIELD_STRENGTH_PROPERTY,
            FIELD_STRENGTH_DIMS,
            ScalarAdmissibility::NonNegative,
        )
        .unwrap(),
    ];
    for (name, dims) in [
        (SCALE, Dims([0; 6])),
        ("density", Density::DIMS),
        ("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
        (MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS),
    ] {
        requirements.push(ScalarPropertyRequirement::try_new(name, dims, positive).unwrap());
    }
    resolve_material_state_point(
        card,
        &point(t),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
pub fn reference(card: &MaterialCard) -> ReversibleMagneticMaterial {
    ReversibleMagneticMaterial::resolve(
        card,
        &point(300.0),
        "B",
        4.0,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
pub fn law(card: &MaterialCard) -> ThermalReversibleMagneticMaterial {
    ThermalReversibleMagneticMaterial::new(
        &reference(card),
        &bundle(card, 300.0),
        &scale_key(),
        "T",
        [290.0, 360.0],
    )
    .unwrap()
}
pub fn h(b: f64) -> f64 {
    if b > 1.0 { 1.0 + 4.0 * (b - 1.0) } else { b }
}
pub fn psi(b: f64) -> f64 {
    if b > 1.0 {
        0.5 + (b - 1.0) + 2.0 * (b - 1.0).powi(2)
    } else {
        0.5 * b * b
    }
}
pub fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-9 * (1.0 + a.abs() + b.abs()),
        "{a:e} != {b:e}"
    );
}
