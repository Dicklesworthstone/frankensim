//! Authored numerical sources; these fixtures make no measured-material claims.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::dielectric::reversible::{
        ELECTRIC_DISPLACEMENT_DIMS, ELECTRIC_FIELD_DIMS, REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY,
        ReversibleDielectricMaterial,
        thermal::{
            ResolvedThermalReversibleDielectricMaterial, ThermalReversibleDielectricMaterial,
        },
    },
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_qty::{
    Density, Dims, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const SCALE: &str = "reversible_dielectric_reference_scale";
fn absolute() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("E", QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS), 0.0)
        .unwrap()
        .with_quantity("T", absolute(), t)
        .unwrap()
}
fn scale_key() -> PropertyKey {
    PropertyKey::new(SCALE, Dims([0; 6]))
}
fn card(field: PropertyValue, maximum: f64, scale: Vec<(f64, f64)>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let domain = ValidityDomain::unconstrained().with_quantity("T", absolute(), 200.0, 500.0);
    for (name, dims, value, validity) in [
        (
            REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY,
            ELECTRIC_DISPLACEMENT_DIMS,
            field,
            domain.clone().with_quantity(
                "E",
                QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS),
                0.0,
                maximum,
            ),
        ),
        (
            SCALE,
            Dims([0; 6]),
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: scale,
                dims: Dims([0; 6]),
            },
            domain.clone(),
        ),
        (
            "density",
            Density::DIMS,
            PropertyValue::Scalar {
                value: 2.0,
                dims: Density::DIMS,
            },
            domain,
        ),
    ] {
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
                    source: "authored separable dielectric thermal fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic nonlinear dielectric".into(),
            phase: "solid".into(),
            process: "thermal numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn curve(knots: &[(f64, f64)]) -> PropertyValue {
    PropertyValue::Curve {
        abscissa: "E".into(),
        abscissa_dims: ELECTRIC_FIELD_DIMS,
        knots: knots.to_vec(),
        dims: ELECTRIC_DISPLACEMENT_DIMS,
    }
}
fn nonlinear(slope: f64) -> MaterialCard {
    card(
        curve(&[(0.0, 0.0), (1.0, 2.0), (3.0, 3.0), (4.0, 7.0)]),
        4.0,
        vec![
            (200.0, 1.0 - 100.0 * slope),
            (300.0, 1.0),
            (500.0, 1.0 + 200.0 * slope),
        ],
    )
}
fn reference(card: &MaterialCard, t: f64, maximum: f64) -> ReversibleDielectricMaterial {
    ReversibleDielectricMaterial::resolve(
        card,
        &point(t),
        "E",
        maximum,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn bundle(card: &MaterialCard, t: f64) -> ResolvedMaterialStatePoint {
    let requirements = [
        ScalarPropertyRequirement::try_new(
            REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY,
            ELECTRIC_DISPLACEMENT_DIMS,
            ScalarAdmissibility::NonNegative,
        )
        .unwrap(),
        ScalarPropertyRequirement::try_new(
            SCALE,
            Dims([0; 6]),
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
        ScalarPropertyRequirement::try_new(
            "density",
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    ];
    resolve_material_state_point(
        card,
        &point(t),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn law(card: &MaterialCard) -> ThermalReversibleDielectricMaterial {
    ThermalReversibleDielectricMaterial::new(
        &reference(card, 300.0, 4.0),
        &bundle(card, 300.0),
        &scale_key(),
        "T",
        [280.0, 340.0],
    )
    .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-30),
        "{a:e} != {b:e}"
    );
}
// Independent fixed-D potential of the original fixture.
fn w(d: f64) -> f64 {
    if d <= 2.0 {
        d * d / 4.0
    } else if d <= 3.0 {
        d * d - 3.0 * d + 3.0
    } else {
        let x = d - 3.0;
        3.0 + 3.0 * x + x * x / 8.0
    }
}
fn e(d: f64) -> f64 {
    if d <= 2.0 {
        d / 2.0
    } else if d <= 3.0 {
        2.0 * d - 3.0
    } else {
        3.0 + (d - 3.0) / 4.0
    }
}
fn epsilon(d: f64) -> f64 {
    if d < 2.0 {
        2.0
    } else if d < 3.0 {
        0.5
    } else {
        4.0
    }
}
fn g(d: f64) -> f64 {
    d * e(d) - w(d)
}
fn offset(at: &ResolvedThermalReversibleDielectricMaterial, d: f64, volume: f64) -> f64 {
    at.temperature_k() * at.state(d.abs(), volume).unwrap().entropy_j_k
}

#[test]
fn source_scale_supplies_reciprocal_entropy_internal_energy_and_negative_field_capacity() {
    for a in [-0.002, 0.0, 0.002] {
        let source = nonlinear(a);
        let law = law(&source);
        near(law.scale_slope_per_k(), a, 2e-15);
        assert_eq!(law.material().properties().len(), 3);
        assert!(!law.reference().is_temperature_scaled());
        for t in [280.0, 300.0, 337.0] {
            let at = law.at_temperature(t).unwrap();
            let s = 1.0 + a * (t - 300.0);
            assert!(at.dielectric().is_temperature_scaled());
            near(at.displacement_scale(), s, 2e-15);
            assert!(std::ptr::eq(
                at.dielectric().source_path(),
                law.reference().source_path()
            ));
            assert_eq!(at.material().properties().len(), 3);
            for (old, new) in law
                .material()
                .properties()
                .iter()
                .zip(at.material().properties())
            {
                assert_eq!(old.requirement(), new.requirement());
                assert_eq!(old.answer().receipt.selected, new.answer().receipt.selected);
            }
            for z in [0.0, 0.4, 1.7, 2.4, 2.9, 4.5, 7.0] {
                let d = at.displacement_scale() * z;
                let actual = at.state(d, 0.7).unwrap();
                let response = at.dielectric().evaluate_displacement(d).unwrap();
                let energy = 0.7 * s * w(z);
                let entropy = 0.7 * a * g(z);
                let capacity = -0.7 * t * a * a * z * z / (s * epsilon(z));
                near(response.field_strength_v_m, e(z), 2e-14);
                near(actual.free_energy_j, energy, 3e-14);
                near(actual.entropy_j_k, entropy, 3e-14);
                near(actual.internal_energy_j, energy + t * entropy, 3e-14);
                near(actual.heat_capacity_j_k, capacity, 3e-14);
            }
        }
        if a != 0.0 {
            let t = 320.0;
            let d = 2.55;
            let h = 1e-3;
            let hd = 1e-5;
            let low = law.at_temperature(t - h).unwrap();
            let high = law.at_temperature(t + h).unwrap();
            let at = law.at_temperature(t).unwrap();
            near(
                -(high.state(d, 1.0).unwrap().free_energy_j
                    - low.state(d, 1.0).unwrap().free_energy_j)
                    / (2.0 * h),
                at.state(d, 1.0).unwrap().entropy_j_k,
                2e-8,
            );
            near(
                (high.state(d, 1.0).unwrap().internal_energy_j
                    - low.state(d, 1.0).unwrap().internal_energy_j)
                    / (2.0 * h),
                at.state(d, 1.0).unwrap().heat_capacity_j_k,
                2e-8,
            );
            near(
                (at.state(d + hd, 1.0).unwrap().entropy_j_k
                    - at.state(d - hd, 1.0).unwrap().entropy_j_k)
                    / (2.0 * hd),
                -(high.dielectric().field_at_displacement(d).unwrap()
                    - low.dielectric().field_at_displacement(d).unwrap())
                    / (2.0 * h),
                2e-8,
            );
        }
    }
}

#[test]
fn every_existing_constitutive_kernel_uses_actual_scale_and_original_knots() {
    let source = nonlinear(0.002);
    let law = law(&source);
    let at = law.at_temperature(335.0).unwrap();
    let physical = at.dielectric();
    let s = at.displacement_scale();
    let reference = law.reference();
    for field in [0.0, 0.25, 1.0, 2.0, 3.0, 3.8, 4.0] {
        let a = physical.evaluate(field).unwrap();
        let b = reference.evaluate(field).unwrap();
        near(a.displacement_c_m2, s * b.displacement_c_m2, 2e-15);
        near(
            a.differential_permittivity_f_m,
            s * b.differential_permittivity_f_m,
            2e-15,
        );
        near(
            a.secant_permittivity_f_m,
            s * b.secant_permittivity_f_m,
            2e-15,
        );
        near(a.energy_density_j_m3, s * b.energy_density_j_m3, 3e-15);
        near(a.coenergy_density_j_m3, s * b.coenergy_density_j_m3, 3e-15);
        near(
            physical.field_at_displacement(a.displacement_c_m2).unwrap(),
            field,
            3e-15,
        );
    }
    for (before, after) in [(-5.5, 2.7), (2.6, 2.2), (-1.0, 1.0), (2.0, 6.3)] {
        near(
            physical
                .energy_change_density_j_m3(s * before, s * after)
                .unwrap(),
            s * (w(after.abs()) - w(before.abs())),
            4e-14,
        );
        near(
            physical
                .coenergy_change_displacement_density_j_m3(s * before, s * after)
                .unwrap(),
            s * (g(after.abs()) - g(before.abs())),
            4e-14,
        );
        near(
            physical
                .energy_change_remainder_density_j_m3(s * before, s * after, s * 1.3)
                .unwrap(),
            s * (w(after.abs()) - w(before.abs()) - e(1.3) * (after - before)),
            4e-14,
        );
    }
    near(
        physical.coenergy_change_density_j_m3(0.5, 3.5).unwrap(),
        s * reference.coenergy_change_density_j_m3(0.5, 3.5).unwrap(),
        4e-15,
    );
    near(
        physical
            .coenergy_change_vector_density_j_m3([1.0, 0.2], [0.8, 1.1])
            .unwrap(),
        s * reference
            .coenergy_change_vector_density_j_m3([1.0, 0.2], [0.8, 1.1])
            .unwrap(),
        4e-15,
    );
    assert!(
        physical
            .evaluate_displacement(physical.maximum_displacement_c_m2().next_up())
            .is_err()
    );
}

#[test]
fn stable_thermal_offset_handles_crossed_knots_signed_charge_and_source_boundaries() {
    for a in [-0.002, 0.002] {
        let source = nonlinear(a);
        let law = law(&source);
        for (ta, tb, da, db) in [
            (280.0, 340.0, 2.1, -3.2),
            (340.0, 280.0, -6.0, 1.4),
            (300.0, 300.0, 2.4, 2.4_f64.next_up()),
            (340.0, 280.0, 6.2, 0.0),
        ] {
            let before = law.at_temperature(ta).unwrap();
            let after = law.at_temperature(tb).unwrap();
            let actual = before.thermal_offset_change_j(&after, da, db, 0.3).unwrap();
            if ta == tb {
                let expected = 0.3 * ta * a * (db - da) * (db + da); // middle segment g(D)=D²-3 at s=1
                near(actual, expected, 3e-15);
            } else {
                near(
                    actual,
                    offset(&after, db, 0.3) - offset(&before, da, 0.3),
                    2e-13,
                );
            }
            let reverse = after.thermal_offset_change_j(&before, db, da, 0.3).unwrap();
            near(actual, -reverse, 4e-13);
        }
        // Both endpoints are valid although the old charge does not fit the
        // final source range: an artificial fixed-old-D split must not refuse.
        let before = law
            .at_temperature(if a > 0.0 { 340.0 } else { 280.0 })
            .unwrap();
        let after = law
            .at_temperature(if a > 0.0 { 280.0 } else { 340.0 })
            .unwrap();
        let da = before.dielectric().maximum_displacement_c_m2();
        assert!(da > after.dielectric().maximum_displacement_c_m2());
        near(
            before
                .thermal_offset_change_j(&after, da, 0.0, 0.7)
                .unwrap(),
            -offset(&before, da, 0.7),
            3e-14,
        );
    }
    // At an original upper endpoint, a sub-ulp scale increase must still
    // retain the constitutive derivative in the thermal work account.
    let source = nonlinear(0.001);
    let law = law(&source);
    let before = law.at_temperature(300.0).unwrap();
    let after = law.at_temperature(300.0_f64.next_up()).unwrap();
    assert_eq!(before.displacement_scale(), after.displacement_scale());
    let d = before.dielectric().maximum_displacement_c_m2();
    let expected =
        (300.0_f64.next_up() - 300.0) * 0.001 * (g(d) - 300.0 * 0.001 * d * d / epsilon(d));
    near(
        before.thermal_offset_change_j(&after, d, d, 1.0).unwrap(),
        expected,
        2e-14,
    );
}

#[test]
fn entropy_work_below_large_unchanged_bias_survives_charge_temperature_cancellation() {
    // s=T/256 and linear epsilon_ref=2. For D1=D0(1+x), T1=T0(1+2x),
    // Delta(TS)=(D0²/4)*x²/(1+2x), far below one ulp of endpoint TS.
    let d = 2.0_f64.powi(40);
    let x = 2.0_f64.powi(-30);
    let maximum = 2.0 * d;
    let source = card(
        curve(&[(0.0, 0.0), (maximum, 2.0 * maximum)]),
        maximum,
        vec![(200.0, 200.0 / 256.0), (256.0, 1.0), (500.0, 500.0 / 256.0)],
    );
    let law = ThermalReversibleDielectricMaterial::new(
        &reference(&source, 256.0, maximum),
        &bundle(&source, 256.0),
        &scale_key(),
        "T",
        [250.0, 400.0],
    )
    .unwrap();
    let before = law.at_temperature(256.0).unwrap();
    let after = law.at_temperature(256.0 * (1.0 + 2.0 * x)).unwrap();
    let expected = (d * d / 4.0) * x * x / (1.0 + 2.0 * x);
    let actual = before
        .thermal_offset_change_j(&after, d, d * (1.0 + x), 1.0)
        .unwrap();
    near(actual, expected, 2e-14);
    let bias = offset(&before, d, 1.0);
    assert!(actual < bias.next_up() - bias);
}

#[test]
fn complete_strip_capacity_bound_includes_both_tangents_and_matches_affine_linear_limit() {
    let source = nonlinear(0.002);
    let law = law(&source);
    let volume = 0.7;
    let minimum = law.minimum_heat_capacity_j_k(4.8, volume).unwrap();
    near(
        minimum,
        -volume * 340.0 * 0.002_f64.powi(2) / 0.96 * 18.0,
        2e-14,
    );
    for t in [280.0, 300.0, 320.0, 340.0] {
        let at = law.at_temperature(t).unwrap();
        for i in 0..=80 {
            assert!(
                at.state(4.8 * f64::from(i) / 80.0, volume)
                    .unwrap()
                    .heat_capacity_j_k
                    >= minimum
            );
        }
    }
    near(
        law.minimum_heat_capacity_j_k(1.92, volume).unwrap(),
        -volume * 340.0 * 0.002_f64.powi(2) / 0.96 * 8.0,
        2e-14,
    );
    assert!(law.minimum_heat_capacity_j_k(7.0, volume).is_err());
    assert_eq!(law.minimum_heat_capacity_j_k(0.0, volume).unwrap(), 0.0);
    assert!(law.minimum_heat_capacity_j_k(1.0, 0.0).is_err());
    let endpoint_source = nonlinear(-0.001);
    let endpoint_law = ThermalReversibleDielectricMaterial::new(
        &reference(&endpoint_source, 300.0, 4.0),
        &bundle(&endpoint_source, 300.0),
        &scale_key(),
        "T",
        [280.0, 500.0],
    )
    .unwrap();
    let envelope = 1.491218878564746;
    assert!(
        endpoint_law
            .minimum_heat_capacity_j_k(envelope, 1.0)
            .unwrap()
            <= endpoint_law
                .at_temperature(500.0)
                .unwrap()
                .state(envelope, 1.0)
                .unwrap()
                .heat_capacity_j_k
    );
    // The nonlinear owner reduces exactly to q²/(2C(T)) in its linear segment.
    let at = law.at_temperature(330.0).unwrap();
    let s = at.displacement_scale();
    let d = 0.5;
    let state = at.state(d, volume).unwrap();
    let free = volume * d * d / (4.0 * s);
    let ratio = 0.002 / s;
    near(state.free_energy_j, free, 2e-15);
    near(state.entropy_j_k, free * ratio, 2e-15);
    near(
        state.heat_capacity_j_k,
        -2.0 * 330.0 * free * ratio * ratio,
        2e-15,
    );
}

#[test]
fn source_admission_refuses_missing_scale_nonaffinity_temperature_surfaces_and_substitution() {
    let source = nonlinear(0.002);
    let original = reference(&source, 300.0, 4.0);
    assert!(
        ThermalReversibleDielectricMaterial::new(
            &original,
            original.material(),
            &scale_key(),
            "T",
            [280.0, 340.0]
        )
        .is_err()
    );
    assert!(
        ThermalReversibleDielectricMaterial::new(
            &original,
            &bundle(&source, 300.0),
            &scale_key(),
            "T",
            [190.0, 340.0]
        )
        .is_err()
    );
    assert!(
        ThermalReversibleDielectricMaterial::new(
            &original,
            &bundle(&source, 300.0),
            &scale_key(),
            "missing_T",
            [280.0, 340.0]
        )
        .is_err()
    );
    for scale in [
        vec![(200.0, 0.9), (300.0, 1.0), (330.0, 1.2), (500.0, 1.4)],
        vec![(200.0, 2.0), (300.0, 1.0), (330.0, -0.1), (500.0, 0.5)],
    ] {
        let other = card(curve(&[(0.0, 0.0), (4.0, 4.0)]), 4.0, scale);
        assert!(
            ThermalReversibleDielectricMaterial::new(
                &reference(&other, 300.0, 4.0),
                &bundle(&other, 300.0),
                &scale_key(),
                "T",
                [280.0, 340.0]
            )
            .is_err()
        );
    }
    let surface = card(
        PropertyValue::Surface2D {
            axes: ["E".into(), "T".into()],
            axis_dims: [ELECTRIC_FIELD_DIMS, Temperature::DIMS],
            samples: vec![
                ([0.0, 200.0], 0.0),
                ([4.0, 200.0], 4.0),
                ([0.0, 500.0], 0.0),
                ([4.0, 500.0], 4.0),
            ],
            triangles: vec![[0, 1, 2], [1, 3, 2]],
            dims: ELECTRIC_DISPLACEMENT_DIMS,
        },
        4.0,
        vec![(200.0, 0.8), (300.0, 1.0), (500.0, 1.4)],
    );
    assert!(
        ThermalReversibleDielectricMaterial::new(
            &reference(&surface, 300.0, 4.0),
            &bundle(&surface, 300.0),
            &scale_key(),
            "T",
            [280.0, 340.0]
        )
        .is_err()
    );
    let a = law(&source).at_temperature(300.0).unwrap();
    let other = nonlinear(-0.002);
    let b = law(&other).at_temperature(300.0).unwrap();
    assert!(a.admit_source_refresh(&b).is_err());
    assert!(a.thermal_offset_change_j(&b, 1.0, 1.0, 1.0).is_err());
    assert!(
        ThermalReversibleDielectricMaterial::new(
            a.dielectric(),
            a.material(),
            &scale_key(),
            "T",
            [280.0, 340.0]
        )
        .is_err()
    );
    assert!(a.state(-1.0, 1.0).is_err());
    assert!(law(&source).at_temperature(340.0_f64.next_up()).is_err());
}
