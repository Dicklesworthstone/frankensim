//! Independent nonlinear dielectric potentials, source support and stable work.
//! The curves are authored numerical fixtures, not measured material data.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::dielectric::reversible::{
        ELECTRIC_DISPLACEMENT_DIMS, ELECTRIC_FIELD_DIMS, REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY,
        ReversibleDielectricMaterial,
    },
    state_point::MaterialPropertySelection,
};
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("E", QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS), 0.0)
        .unwrap()
        .with_quantity("T", temperature(), t)
        .unwrap()
}
fn card(value: PropertyValue, interpolation: InterpolationPolicy) -> MaterialCard {
    let maximum_field = match &value {
        PropertyValue::Curve { knots, .. } => knots.last().unwrap().0,
        PropertyValue::Surface2D { samples, .. } => samples
            .iter()
            .map(|(point, _)| point[0])
            .fold(0.0_f64, f64::max),
        _ => 4.0,
    };
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY,
                ELECTRIC_DISPLACEMENT_DIMS,
            ),
            value,
            validity: ValidityDomain::unconstrained()
                .with_quantity("T", temperature(), 250.0, 450.0)
                .with_quantity(
                    "E",
                    QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS),
                    0.0,
                    maximum_field,
                ),
            interpolation,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "authored reversible D(E) fixture".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric".into(),
            phase: "solid".into(),
            process: "reversible numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn curve(knots: &[(f64, f64)]) -> MaterialCard {
    card(
        PropertyValue::Curve {
            abscissa: "E".into(),
            abscissa_dims: ELECTRIC_FIELD_DIMS,
            knots: knots.to_vec(),
            dims: ELECTRIC_DISPLACEMENT_DIMS,
        },
        InterpolationPolicy::LinearInside,
    )
}
fn law(source: &MaterialCard, maximum: f64) -> ReversibleDielectricMaterial {
    ReversibleDielectricMaterial::resolve(
        source,
        &point(300.0),
        "E",
        maximum,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}
fn nonlinear() -> MaterialCard {
    curve(&[(0.0, 0.0), (1.0, 2.0), (3.0, 3.0), (4.0, 7.0)])
}
// Independent closed forms for the fixture, in displacement coordinates.
fn energy(d: f64) -> f64 {
    let d = d.abs();
    if d <= 2.0 {
        d * d / 4.0
    } else if d <= 3.0 {
        let x = d - 2.0;
        1.0 + x + x * x
    } else {
        let x = d - 3.0;
        3.0 + 3.0 * x + x * x / 8.0
    }
}
fn field(d: f64) -> f64 {
    let magnitude = d.abs();
    d.signum()
        * if magnitude <= 2.0 {
            magnitude / 2.0
        } else if magnitude <= 3.0 {
            1.0 + 2.0 * (magnitude - 2.0)
        } else {
            3.0 + (magnitude - 3.0) / 4.0
        }
}
fn coenergy(e: f64) -> f64 {
    if e <= 1.0 {
        e * e
    } else if e <= 3.0 {
        let x = e - 1.0;
        1.0 + 2.0 * x + 0.25 * x * x
    } else {
        let x = e - 3.0;
        6.0 + 3.0 * x + 2.0 * x * x
    }
}

#[test]
fn source_path_supplies_independent_energy_coenergy_inverse_and_tangents() {
    let source = nonlinear();
    let law = law(&source, 4.0);
    let cloned = law.clone();
    assert!(std::ptr::eq(law.source_path(), cloned.source_path()));
    assert_eq!(law.card(), &source);
    assert_eq!(law.material().card_identity(), source.content_hash());
    assert_eq!(law.maximum_displacement_c_m2(), 7.0);
    assert_eq!(law.maximum_field_strength_v_m(), 4.0);
    for i in 0..=32 {
        let e = f64::from(i) / 8.0;
        let response = law.evaluate(e).unwrap();
        near(
            response.energy_density_j_m3,
            energy(response.displacement_c_m2),
            2e-14,
        );
        near(response.coenergy_density_j_m3, coenergy(e), 2e-14);
        near(
            response.energy_density_j_m3 + response.coenergy_density_j_m3,
            e * response.displacement_c_m2,
            2e-14,
        );
        near(
            law.field_at_displacement(response.displacement_c_m2)
                .unwrap(),
            e,
            2e-14,
        );
        let inverse = law
            .evaluate_displacement(response.displacement_c_m2)
            .unwrap();
        near(
            inverse.energy_density_j_m3,
            response.energy_density_j_m3,
            2e-14,
        );
        let canonical = law
            .material()
            .at_query_point(&source, &point(300.0).with("E", e).unwrap())
            .unwrap();
        near(
            response.displacement_c_m2,
            canonical
                .property(REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY)
                .unwrap()
                .value_si(),
            2e-14,
        );
    }
    for (e, answer) in &law.source_path().knots {
        assert_eq!(
            law.evaluate(*e).unwrap().displacement_c_m2,
            answer.evidence.value.value
        );
        assert_eq!(
            law.field_at_displacement(answer.evidence.value.value)
                .unwrap(),
            *e
        );
        assert_eq!(
            answer.receipt.selected,
            law.material()
                .property(REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY)
                .unwrap()
                .answer()
                .receipt
                .selected
        );
    }
    assert_eq!(law.evaluate(0.0).unwrap().secant_permittivity_f_m, 2.0);
    assert_eq!(
        law.evaluate(1.0).unwrap().differential_permittivity_f_m,
        0.5
    );
    assert_eq!(
        law.evaluate(4.0).unwrap().differential_permittivity_f_m,
        4.0
    );
    let step = 1e-5;
    let e = 2.0;
    let d = 2.5;
    near(
        (law.evaluate(e + step).unwrap().coenergy_density_j_m3
            - law.evaluate(e - step).unwrap().coenergy_density_j_m3)
            / (2.0 * step),
        d,
        1e-9,
    );
    near(
        (law.evaluate_displacement(d + step)
            .unwrap()
            .energy_density_j_m3
            - law
                .evaluate_displacement(d - step)
                .unwrap()
                .energy_density_j_m3)
            / (2.0 * step),
        e,
        1e-9,
    );
    near(
        (law.evaluate(e + step).unwrap().displacement_c_m2
            - law.evaluate(e - step).unwrap().displacement_c_m2)
            / (2.0 * step),
        0.5,
        1e-9,
    );
    let low = law.evaluate(1e-8).unwrap();
    assert!((low.energy_density_j_m3 / 1e-16 - 1.0).abs() < 1e-14);
    assert_eq!(low.energy_density_j_m3, low.coenergy_density_j_m3);
}

#[test]
fn signed_work_and_native_remainders_cross_original_knots_and_both_polarities() {
    let law = law(&nonlinear(), 4.0);
    for (before, after) in [
        (-7.0, 7.0),
        (7.0, -2.0),
        (-2.0, -3.0),
        (2.5, 4.0),
        (4.0, 2.5),
        (0.0, -7.0),
        (7.0, 0.0),
    ] {
        near(
            law.energy_change_density_j_m3(before, after).unwrap(),
            energy(after) - energy(before),
            2e-14,
        );
        for reference in [-3.0, 0.0, 2.0] {
            let expected = energy(after) - energy(before) - field(reference) * (after - before);
            near(
                law.energy_change_remainder_density_j_m3(before, after, reference)
                    .unwrap(),
                expected,
                2e-14,
            );
        }
    }
    assert_eq!(law.energy_change_density_j_m3(-7.0, 7.0).unwrap(), 0.0);
    for (before, after) in [(0.0, 4.0), (4.0, 0.0), (0.5, 2.0), (2.0, 3.5), (3.5, 1.0)] {
        near(
            law.coenergy_change_density_j_m3(before, after).unwrap(),
            coenergy(after) - coenergy(before),
            2e-14,
        );
    }
}

#[test]
fn direct_integrals_and_work_retain_small_values_under_large_bias() {
    let almost_flat = law(
        &curve(&[(0.0, 0.0), (1e-30, 1.0), (1.0, 1.0 + f64::EPSILON)]),
        1.0,
    );
    let response = almost_flat.evaluate(0.5).unwrap();
    assert_eq!(
        0.5 * response.displacement_c_m2 - response.coenergy_density_j_m3,
        0.0
    );
    let expected = 0.125 * f64::EPSILON + 0.5e-30;
    assert!((response.energy_density_j_m3 / expected - 1.0).abs() < 1e-14);
    let biased = law(
        &curve(&[(0.0, 0.0), (1e16, 1e16), (1e16 + 128.0, 1e16 + 256.0)]),
        1e16 + 128.0,
    );
    let before = 1e16 + 32.0;
    let after = 1e16 + 36.0;
    near(
        biased.energy_change_density_j_m3(before, after).unwrap(),
        (1e16 + 17.0) * 4.0,
        4.0 * f64::EPSILON,
    );
    assert_eq!(
        biased
            .energy_change_remainder_density_j_m3(before, after, before)
            .unwrap(),
        4.0
    );
    assert_eq!(
        biased
            .energy_change_remainder_density_j_m3(before, after, 1e16 + 34.0)
            .unwrap(),
        0.0
    );
    assert_eq!(
        biased
            .energy_change_remainder_density_j_m3(-before, -after, -before)
            .unwrap(),
        4.0
    );
    near(
        biased
            .coenergy_change_density_j_m3(1e16 + 16.0, 1e16 + 18.0)
            .unwrap(),
        (1e16 + 34.0) * 2.0,
        2.0 * f64::EPSILON,
    );
}

#[test]
fn vector_coenergy_retains_transverse_work_and_midpoint_corrections_under_bias() {
    let law = law(&curve(&[(0.0, 0.0), (2e16, 2e16)]), 2e16);
    for sign in [-1.0, 1.0] {
        let before = [1e16, 0.0];
        let after = [1e16, sign];
        assert_eq!(
            law.coenergy_change_vector_density_j_m3(before, after)
                .unwrap(),
            0.5
        );
        assert_eq!(
            law.coenergy_change_vector_density_j_m3(after, before)
                .unwrap(),
            -0.5
        );
    }
    // Each large component contributes approximately +/-2e16. Their corrected
    // midpoints retain the exact surviving sum: (2²+(-2)²)/2 = 4.
    let before = [1e16, 1e16];
    let after = [1e16 + 2.0, 1e16 - 2.0];
    assert_eq!(
        law.coenergy_change_vector_density_j_m3(before, after)
            .unwrap(),
        4.0
    );
    assert_eq!(
        law.coenergy_change_vector_density_j_m3(after, before)
            .unwrap(),
        -4.0
    );
    assert_eq!(
        law.coenergy_change_vector_density_j_m3([1e16, 1.0], [1e16, -1.0])
            .unwrap(),
        0.0
    );
}

#[test]
fn vector_coenergy_corrects_unequal_rounded_radii_on_the_original_nonlinear_curve() {
    let source = curve(&[(0.0, 0.0), (5e7, 1e8), (2e8, 1.75e8)]);
    let material = law(&source, 2e8);
    let before = [1e8_f64, 1.0];
    let after = [1e8_f64, 2.0];
    let a = before[0].hypot(before[1]);
    let b = after[0].hypot(after[1]);
    assert_ne!(a, b);
    // In this segment D(E)=E/2+7.5e7. The exact squared-radius increment
    // is 3/2; the radial increment is 3/(sqrt(1e16+4)+sqrt(1e16+1)).
    // To double precision, delta W*=0.5*(3/2)+7.5e7*(3/(2e8))=1.875.
    near(
        material
            .coenergy_change_vector_density_j_m3(before, after)
            .unwrap(),
        1.875,
        2e-14,
    );
    near(
        material
            .coenergy_change_vector_density_j_m3(after, before)
            .unwrap(),
        -1.875,
        2e-14,
    );
    assert!((material.coenergy_change_density_j_m3(a, b).unwrap() - 1.875).abs() > 1e-3);
    // The correction also retains the original exact crossed-segment integral.
    let material = law(&nonlinear(), 4.0);
    near(
        material
            .coenergy_change_vector_density_j_m3([0.6, 0.8], [1.8, 2.4])
            .unwrap(),
        5.0,
        2e-14,
    );
    assert!(
        material
            .coenergy_change_vector_density_j_m3([0.0, 0.0], [4.0, 4.0])
            .is_err()
    );
    assert!(
        material
            .coenergy_change_vector_density_j_m3([f64::NAN, 0.0], [0.0, 0.0])
            .is_err()
    );
}

#[test]
fn vector_coenergy_weights_before_unrepresentable_squares_and_differences() {
    for (field, displacement) in [(1e200, 1.0), (1e-200, 1.0)] {
        let material = law(
            &curve(&[(0.0, 0.0), (2.0 * field, 2.0 * displacement)]),
            2.0 * field,
        );
        let actual = material
            .coenergy_change_vector_density_j_m3([field, 0.0], [field, field])
            .unwrap();
        let expected = 0.5 * field * displacement;
        assert!(
            (actual / expected - 1.0).abs() < 2e-14,
            "{actual} != {expected}"
        );
    }
    let material = law(&curve(&[(0.0, 0.0), (1e308, 1.0)]), 1e308);
    let actual = material
        .coenergy_change_vector_density_j_m3([-1e308], [9e307])
        .unwrap();
    let expected = material.coenergy_change_density_j_m3(1e308, 9e307).unwrap();
    assert!((actual / expected - 1.0).abs() < 2e-14);
    assert_eq!(
        material
            .coenergy_change_vector_density_j_m3([-1e308], [1e308])
            .unwrap(),
        0.0
    );
}

#[test]
fn fixed_temperature_surface_retains_original_breakpoints_and_receipts() {
    let source = card(
        PropertyValue::Surface2D {
            axes: ["E".into(), "T".into()],
            axis_dims: [ELECTRIC_FIELD_DIMS, Temperature::DIMS],
            samples: vec![
                ([0.0, 250.0], 0.0),
                ([1.0, 250.0], 2.0),
                ([4.0, 250.0], 5.0),
                ([0.0, 450.0], 0.0),
                ([1.0, 450.0], 3.0),
                ([4.0, 450.0], 8.0),
            ],
            triangles: vec![[0, 1, 4], [0, 4, 3], [1, 2, 5], [1, 5, 4]],
            dims: ELECTRIC_DISPLACEMENT_DIMS,
        },
        InterpolationPolicy::LinearInside,
    );
    let law = law(&source, 4.0);
    assert!(law.source_path().knots.len() > 3);
    assert_eq!(
        law.material()
            .property(REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY)
            .unwrap()
            .source_axes()
            .len(),
        2
    );
    for i in 0..=32 {
        let e = f64::from(i) / 8.0;
        let canonical = law
            .material()
            .at_query_point(&source, &point(300.0).with("E", e).unwrap())
            .unwrap();
        let d = canonical
            .property(REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY)
            .unwrap()
            .value_si();
        let response = law.evaluate(e).unwrap();
        near(response.displacement_c_m2, d, 2e-14);
        near(law.field_at_displacement(d).unwrap(), e, 2e-14);
        near(
            response.energy_density_j_m3 + response.coenergy_density_j_m3,
            e * d,
            2e-14,
        );
    }
}

#[test]
fn source_domain_monotonicity_and_unowned_response_conditions_refuse() {
    let resolve = |source: &MaterialCard, point: &QueryPoint, maximum| {
        ReversibleDielectricMaterial::resolve(
            source,
            point,
            "E",
            maximum,
            MaterialPropertySelection::SingleClaimOnly,
        )
    };
    for knots in [
        vec![(0.0, 0.1), (4.0, 4.0)],
        vec![(0.0, 0.0), (2.0, 2.0), (4.0, 2.0)],
        vec![(0.0, 0.0), (2.0, 2.0), (4.0, 1.0)],
    ] {
        assert!(resolve(&curve(&knots), &point(300.0), 4.0).is_err());
    }
    let source = nonlinear();
    for maximum in [0.0, -1.0, 4.01, f64::INFINITY] {
        assert!(resolve(&source, &point(300.0), maximum).is_err());
    }
    assert!(resolve(&source, &point(451.0), 4.0).is_err());
    assert!(resolve(&source, &QueryPoint::new().with("E", 0.0).unwrap(), 4.0).is_err());
    let wrong_units = QueryPoint::new()
        .with_quantity(
            "E",
            QuantitySpec::dimensional(ELECTRIC_DISPLACEMENT_DIMS),
            0.0,
        )
        .unwrap();
    assert!(resolve(&source, &wrong_units, 4.0).is_err());
    let frequency = point(300.0)
        .with_quantity(
            "frequency",
            QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
            1.0,
        )
        .unwrap();
    assert!(resolve(&source, &frequency, 4.0).is_err());
    let constant = card(
        PropertyValue::Scalar {
            value: 0.0,
            dims: ELECTRIC_DISPLACEMENT_DIMS,
        },
        InterpolationPolicy::ConstantWithinValidity,
    );
    assert!(resolve(&constant, &point(300.0), 4.0).is_err());
    let discrete = card(
        PropertyValue::Curve {
            abscissa: "E".into(),
            abscissa_dims: ELECTRIC_FIELD_DIMS,
            knots: vec![(0.0, 0.0), (4.0, 7.0)],
            dims: ELECTRIC_DISPLACEMENT_DIMS,
        },
        InterpolationPolicy::TabulatedOnly,
    );
    assert!(resolve(&discrete, &point(300.0), 4.0).is_err());
    let law = law(&source, 4.0);
    for e in [-0.1, 4.01, f64::NAN] {
        assert!(law.evaluate(e).is_err());
        assert!(law.coenergy_change_density_j_m3(0.0, e).is_err());
    }
    for d in [-0.1, 7.01, f64::INFINITY] {
        assert!(law.evaluate_displacement(d).is_err());
    }
    assert!(law.energy_change_density_j_m3(0.0, -7.01).is_err());
    assert!(
        law.energy_change_remainder_density_j_m3(0.0, 1.0, 7.01)
            .is_err()
    );
}
