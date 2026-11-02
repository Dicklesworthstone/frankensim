//! Authored equilibrium/storage oracles, not measured sorption datasets.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::*;
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{MoistureBasis, QuantityKind, SemanticType, ValueForm},
};

const NAME: &str = "equilibrium_dry_basis_moisture";
fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn moisture_quantity(basis: MoistureBasis) -> QuantitySpec {
    quantity(QuantityKind::Moisture(basis))
}
fn point(activity: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            300.0,
        )
        .unwrap()
        .with_quantity(
            "relative_humidity",
            moisture_quantity(MoistureBasis::RelativeHumidity),
            activity,
        )
        .unwrap()
}
fn card(value: PropertyValue, output: MoistureBasis) -> MaterialCard {
    let interpolation = match value {
        PropertyValue::Scalar { .. } => InterpolationPolicy::ConstantWithinValidity,
        _ => InterpolationPolicy::LinearInside,
    };
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::with_quantity(NAME, moisture_quantity(output)),
            value,
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "temperature",
                    quantity(QuantityKind::AbsoluteTemperature),
                    290.0,
                    310.0,
                )
                .with_quantity(
                    "relative_humidity",
                    moisture_quantity(MoistureBasis::RelativeHumidity),
                    0.0,
                    1.0,
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation,
            observations: vec![],
            provenance: Provenance {
                source: "authored monotone sorption verification law".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured hygroscopic solid".into(),
            phase: "solid".into(),
            process: "verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn curve(knots: Vec<(f64, f64)>) -> MaterialCard {
    card(
        PropertyValue::Curve {
            abscissa: "relative_humidity".into(),
            abscissa_dims: Dims::NONE,
            knots,
            dims: Dims::NONE,
        },
        MoistureBasis::DryMass,
    )
}
fn binding(
    card: &MaterialCard,
    point: &QueryPoint,
    output: MoistureBasis,
) -> ResolvedMaterialStatePoint {
    resolve_material_state_point(
        card,
        point,
        &[ScalarPropertyRequirement::try_with_quantity(
            NAME,
            moisture_quantity(output),
            ScalarAdmissibility::NonNegative,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn chart(card: &MaterialCard, bounds: [f64; 2]) -> MonotoneMoistureSorption {
    MonotoneMoistureSorption::from_binding(
        card,
        &binding(card, &point(bounds[0]), MoistureBasis::DryMass),
        NAME,
        bounds,
    )
    .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

#[test]
fn nonlinear_isotherm_inverts_and_integrates_every_storage_segment() {
    let c = curve(vec![(0.0, 0.0), (0.25, 0.05), (0.5, 0.15), (1.0, 0.65)]);
    let s = chart(&c, [0.0, 1.0]);
    for (activity, moisture) in [
        (0.0, 0.0),
        (0.125, 0.025),
        (0.25, 0.05),
        (0.375, 0.1),
        (0.5, 0.15),
        (0.75, 0.4),
        (1.0, 0.65),
    ] {
        near(s.moisture(activity).unwrap(), moisture, 1e-15);
        near(s.activity(moisture).unwrap(), activity, 1e-15);
    }
    near(s.capacity(0.0).unwrap(), 0.2, 1e-15);
    near(s.capacity(0.25).unwrap(), 0.4, 1e-15); // Defined right derivative.
    near(s.capacity(1.0).unwrap(), 1.0, 1e-15);
    near(s.mean_capacity(0.125, 0.75).unwrap(), 0.6, 1e-15);
    near(s.moisture_change(0.125, 0.75).unwrap(), 0.375, 1e-15);
    near(s.moisture_change(0.75, 0.125).unwrap(), -0.375, 1e-15);
    assert_eq!(s.moisture_change(0.25, 0.25).unwrap(), 0.0);
    assert_eq!(s.mean_capacity(0.125, 0.75), s.mean_capacity(0.75, 0.125));
    assert_eq!(s.card_identity(), c.content_hash());
    assert_eq!(s.temperature_kelvin(), 300.0);
    assert_eq!(s.source_point(), &point(0.0));
    assert_eq!(s.knot_count(), 4);
    assert_eq!(s, chart(&c, [0.0, 1.0]));
    assert_ne!(s.identity(), chart(&c, [0.25, 0.75]).identity());
}

#[test]
fn distinct_materials_partition_storage_at_a_common_activity() {
    let a = chart(&curve(vec![(0.0, 0.0), (1.0, 0.2)]), [0.0, 1.0]);
    let b = chart(&curve(vec![(0.0, 0.0), (0.5, 0.3), (1.0, 0.8)]), [0.0, 1.0]);
    let common = 0.5;
    let ma = a.moisture(common).unwrap();
    let mb = b.moisture(common).unwrap();
    near(mb / ma, 3.0, 1e-15);
    assert_eq!(a.activity(ma).unwrap(), b.activity(mb).unwrap());
    // A fixed total water store uses material-local dry masses and capacities.
    near(
        2.0 * a.moisture_change(0.25, 0.5).unwrap() + 3.0 * b.moisture_change(0.25, 0.5).unwrap(),
        0.55,
        2e-15,
    );
}

#[test]
fn narrow_changes_are_not_lost_to_a_large_moisture_datum() {
    let c = curve(vec![(0.0, 1e12), (1.0, 1e12 + 2.0)]);
    let s = chart(&c, [0.0, 1.0]);
    let a = 0.9_f64;
    let b = f64::from_bits(a.to_bits() + 1);
    assert_eq!(s.moisture(a).unwrap(), s.moisture(b).unwrap());
    assert_eq!(s.mean_capacity(a, b).unwrap(), 2.0);
    assert_eq!(s.moisture_change(a, b).unwrap(), 2.0 * (b - a));
    let bounded = chart(&c, [0.25, 0.75]);
    for a in [0.24, 0.76, f64::NAN, f64::INFINITY] {
        assert!(bounded.moisture(a).is_err());
        assert!(bounded.capacity(a).is_err());
        assert!(bounded.mean_capacity(0.5, a).is_err());
    }
    for m in [1e12, 1e12 + 2.0, f64::NAN] {
        assert!(bounded.activity(m).is_err());
    }
}

#[test]
fn noninvertible_sources_wrong_basis_and_duplicate_coordinates_refuse() {
    for c in [
        curve(vec![(0.0, 0.0), (0.5, 0.5), (1.0, 0.25)]),
        curve(vec![(0.0, 0.0), (0.5, 0.0), (1.0, 0.25)]),
        card(
            PropertyValue::Scalar {
                value: 0.2,
                dims: Dims::NONE,
            },
            MoistureBasis::DryMass,
        ),
    ] {
        let bind = binding(&c, &point(0.0), MoistureBasis::DryMass);
        assert!(MonotoneMoistureSorption::from_binding(&c, &bind, NAME, [0.0, 1.0]).is_err());
    }
    let c = curve(vec![(0.0, 0.0), (1.0, 0.2)]);
    let bind = binding(&c, &point(0.5), MoistureBasis::DryMass);
    let other = curve(vec![(0.0, 0.0), (1.0, 0.3)]);
    assert!(MonotoneMoistureSorption::from_binding(&other, &bind, NAME, [0.0, 1.0]).is_err());
    for bounds in [[0.5, 0.5], [1.0, 0.0], [-0.1, 0.5], [0.0, 1.1]] {
        assert!(MonotoneMoistureSorption::from_binding(&c, &bind, NAME, bounds).is_err());
    }
    let wet = card(
        PropertyValue::Curve {
            abscissa: "relative_humidity".into(),
            abscissa_dims: Dims::NONE,
            knots: vec![(0.0, 0.0), (1.0, 0.2)],
            dims: Dims::NONE,
        },
        MoistureBasis::WetMass,
    );
    let wet_bind = binding(&wet, &point(0.5), MoistureBasis::WetMass);
    assert!(MonotoneMoistureSorption::from_binding(&wet, &wet_bind, NAME, [0.0, 1.0]).is_err());
    let duplicate = point(0.5)
        .with_quantity(
            "second_activity",
            moisture_quantity(MoistureBasis::RelativeHumidity),
            0.5,
        )
        .unwrap();
    let duplicate_bind = binding(&c, &duplicate, MoistureBasis::DryMass);
    assert!(MonotoneMoistureSorption::from_binding(&c, &duplicate_bind, NAME, [0.0, 1.0]).is_err());
}

fn surface(gap: bool) -> MaterialCard {
    let (samples, triangles) = if gap {
        (
            vec![
                ([0.0, 290.0], 0.0),
                ([0.25, 290.0], 0.1),
                ([0.0, 310.0], 0.0),
                ([0.25, 310.0], 0.1),
                ([0.75, 290.0], 0.3),
                ([1.0, 290.0], 0.4),
                ([0.75, 310.0], 0.3),
                ([1.0, 310.0], 0.4),
            ],
            vec![[0, 1, 2], [1, 3, 2], [4, 5, 6], [5, 7, 6]],
        )
    } else {
        (
            vec![
                ([0.0, 290.0], 0.0),
                ([1.0, 290.0], 0.4),
                ([0.0, 310.0], 0.2),
                ([1.0, 310.0], 0.8),
            ],
            vec![[0, 1, 2], [1, 3, 2]],
        )
    };
    card(
        PropertyValue::Surface2D {
            axes: ["relative_humidity".into(), "temperature".into()],
            axis_dims: [Dims::NONE, Temperature::DIMS],
            samples,
            triangles,
            dims: Dims::NONE,
        },
        MoistureBasis::DryMass,
    )
}

#[test]
fn fixed_temperature_surface_keeps_crossed_cells_and_refuses_source_holes() {
    let complete = surface(false);
    let s = chart(&complete, [0.0, 1.0]);
    assert_eq!(s.knot_count(), 3);
    for a in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let exact = binding(&complete, &point(a), MoistureBasis::DryMass)
            .property(NAME)
            .unwrap()
            .value_si();
        near(s.moisture(a).unwrap(), exact, 1e-15);
        near(s.activity(exact).unwrap(), a, 1e-15);
    }
    near(s.moisture_change(0.0, 1.0).unwrap(), 0.5, 1e-15);
    let gap = surface(true);
    let lower = binding(&gap, &point(0.0), MoistureBasis::DryMass);
    binding(&gap, &point(1.0), MoistureBasis::DryMass);
    assert!(matches!(
        MonotoneMoistureSorption::from_binding(&gap, &lower, NAME, [0.0, 1.0]),
        Err(MaterialStatePointError::Support { .. })
    ));
}
