//! Synthetic affine-profile integration references; no material calibration claim.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, HardnessLoadStep, HardnessTestContext, InterpolationPolicy, MaterialCard,
    MaterialStateId, ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::hardness::INDENTATION_HARDNESS_QUANTITY;
use fs_material::hardness::depth::{
    DepthHardnessError, EXPOSED_MATERIAL_DEPTH_AXIS as Z, ResolvedDepthHardness,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QtyAny, QuantitySpec, Temperature};

fn profile(knots: Vec<(f64, f64)>) -> ResolvedDepthHardness {
    let maximum = knots.last().unwrap().0;
    let length = Dims([1, 0, 0, 0, 0, 0]);
    let provenance = Provenance {
        source: "synthetic hardness-depth integral fixture".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    };
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "synthetic graded specimen".into(),
            method: "synthetic force/projected-area indentation".into(),
            artifact: fs_blake3::hash_domain("fixture", b"depth integral"),
            caveats: "numerical references only".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    let key = PropertyKey::with_quantity("H", INDENTATION_HARDNESS_QUANTITY)
        .with_hardness_test(
            HardnessTestContext::new(
                "synthetic sphere radius 1 mm",
                vec![
                    HardnessLoadStep::new(
                        QtyAny::new(10.0, Dims([1, 1, -2, 0, 0, 0])),
                        QtyAny::new(1.0, Dims([0, 0, 1, 0, 0, 0])),
                    )
                    .unwrap(),
                ],
                "synthetic depth profile revision 1",
                observation,
            )
            .unwrap(),
        )
        .unwrap();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Curve {
                abscissa: Z.into(),
                abscissa_dims: length,
                knots,
                dims: Pressure::DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                Z,
                QuantitySpec::dimensional(length),
                0.0,
                maximum,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![observation],
            provenance,
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic solid".into(),
            phase: "solid".into(),
            process: "synthetic depth profile".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    ResolvedDepthHardness::new(
        &card,
        &QueryPoint::new(),
        key,
        MaterialPropertySelection::SingleClaimOnly,
        maximum,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(f64::MIN_POSITIVE),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn constant_hardening_and_softening_profiles_match_independent_affine_integrals() {
    for (h0, h1) in [(1.0e9, 1.0e9), (1.0e8, 1.0e9), (1.0e9, 1.0e8)] {
        let source = profile(vec![(0.0, h0), (1.0, h1)]);
        for (start, end) in [(0.0, 0.25), (0.125, 0.5), (0.25, 1.0)] {
            let expected = h0 * (end - start) + 0.5 * (h1 - h0) * (end * end - start * start);
            close(
                source.hardness_integral_between(start, end).unwrap(),
                expected,
                3.0e-15,
            );
            let next = source.advance_by_integral(start, expected).unwrap();
            close(next.depth_m, end, 3.0e-15);
            close(next.increment_m, end - start, 5.0e-15);
            close(
                next.path_mean_hardness_pa,
                expected / (end - start),
                5.0e-15,
            );
        }
    }
}

#[test]
fn crossed_knots_and_path_partition_retain_the_same_original_profile() {
    let source = profile(vec![
        (0.0, 100.0),
        (1.0, 100.0),
        (2.0, 300.0),
        (3.0, 50.0),
        (4.0, 50.0),
    ]);
    let refined = profile(vec![
        (0.0, 100.0),
        (0.5, 100.0),
        (1.0, 100.0),
        (1.5, 200.0),
        (2.0, 300.0),
        (2.5, 175.0),
        (3.0, 50.0),
        (3.5, 50.0),
        (4.0, 50.0),
    ]);
    assert_eq!(source.hardness_integral_between(0.0, 4.0).unwrap(), 525.0);
    let advance = source.advance_by_integral(0.5, 450.0).unwrap();
    assert_eq!(advance.depth_m, 3.5);
    assert_eq!(advance.path_mean_hardness_pa, 150.0);
    for endpoint in [1.0, 2.0, 3.0, 4.0] {
        let integral = source.hardness_integral_between(0.0, endpoint).unwrap();
        assert_eq!(
            source.advance_by_integral(0.0, integral).unwrap().depth_m,
            endpoint
        );
    }
    let first = source.hardness_integral_between(0.37, 1.2).unwrap();
    let second = source.hardness_integral_between(1.2, 3.4).unwrap();
    let split = source
        .advance_by_integral(
            source.advance_by_integral(0.37, first).unwrap().depth_m,
            second,
        )
        .unwrap();
    let whole = source.advance_by_integral(0.37, first + second).unwrap();
    let refined_whole = refined.advance_by_integral(0.37, first + second).unwrap();
    close(split.depth_m, 3.4, 3.0e-15);
    close(whole.depth_m, split.depth_m, 3.0e-15);
    close(refined_whole.depth_m, whole.depth_m, 3.0e-15);
}

#[test]
fn small_late_increments_do_not_subtract_large_reference_integrals() {
    let source = profile(vec![(0.0, 1.0e12), (1.0, 2.0e12)]);
    let start = 0.875;
    let increment = 2.0_f64.powi(-40);
    let end = start + increment;
    let expected = (1.0e12 * (1.0 + start) + 0.5e12 * increment) * increment;
    close(
        source.hardness_integral_between(start, end).unwrap(),
        expected,
        2.0e-15,
    );
    let next = source.advance_by_integral(start, expected).unwrap();
    assert_eq!(next.depth_m.to_bits(), end.to_bits());
    assert_eq!(next.increment_m.to_bits(), increment.to_bits());
}

#[test]
fn scaled_inverse_avoids_squaring_dimensional_hardness() {
    for scale in [1.0e-200, 1.0, 1.0e200] {
        let source = profile(vec![(0.0, scale), (0.1, 2.0 * scale)]);
        let expected = scale * 0.0725;
        close(
            source.hardness_integral_between(0.02, 0.07).unwrap(),
            expected,
            4.0e-15,
        );
        close(
            source.advance_by_integral(0.02, expected).unwrap().depth_m,
            0.07,
            4.0e-15,
        );
    }
}

#[test]
fn source_limits_invalid_measures_and_lost_depth_resolution_refuse() {
    let source = profile(vec![(0.0, 100.0), (1.0, 200.0)]);
    let before = source.clone();
    for measure in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(source.advance_by_integral(0.0, measure).is_err());
    }
    for depth in [-1.0, 1.01, f64::NAN, f64::INFINITY] {
        assert!(source.advance_by_integral(depth, 0.0).is_err());
    }
    assert!(source.advance_by_integral(0.0, 150.000001).is_err());
    assert!(source.advance_by_integral(1.0, 1.0).is_err());
    assert!(source.advance_by_integral(0.5, f64::MIN_POSITIVE).is_err());
    assert!(source.hardness_integral_between(0.5, 0.25).is_err());
    assert!(source.hardness_integral_between(0.0, 1.01).is_err());
    assert_eq!(source.hardness_integral_between(1.0, 1.0).unwrap(), 0.0);
    let unchanged = source.advance_by_integral(1.0, 0.0).unwrap();
    assert_eq!(unchanged.depth_m, 1.0);
    assert_eq!(unchanged.increment_m, 0.0);
    assert_eq!(unchanged.path_mean_hardness_pa, 200.0);
    assert_eq!(source, before);
}

fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn thermal_surface(negative_interior: bool, hole: bool) -> PropertyValue {
    let depths = [0.0, 0.4, 0.6, 1.0];
    let temperatures = [300.0, 350.0, 400.0];
    let mut samples = Vec::new();
    for temperature in temperatures {
        for depth in depths {
            let h = if negative_interior && temperature == 400.0 && depth > 0.0 && depth < 1.0 {
                -50.0
            } else {
                900.0 - 2.0 * (temperature - 300.0) + 400.0 * depth
            };
            samples.push(([depth, temperature], h));
        }
    }
    let mut triangles = Vec::new();
    for row in 0..2_u32 {
        for column in 0..3_u32 {
            if hole && row == 1 && column == 1 {
                continue;
            }
            let lower = 4 * row + column;
            triangles.push([lower, lower + 1, lower + 5]);
            triangles.push([lower, lower + 5, lower + 4]);
        }
    }
    PropertyValue::Surface2D {
        axes: [Z.into(), "T".into()],
        axis_dims: [Dims([1, 0, 0, 0, 0, 0]), Temperature::DIMS],
        samples,
        triangles,
        dims: Pressure::DIMS,
    }
}

fn thermal_profile(
    value: PropertyValue,
    temperature_quantity: Option<QuantitySpec>,
    competing_hot_claim: bool,
) -> (MaterialCard, ResolvedDepthHardness) {
    let pressure =
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
    let environment = ValidityDomain::unconstrained()
        .with_quantity("P", pressure, 101_325.0, 101_325.0)
        .with_quantity(
            Z,
            QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
            0.0,
            1.0,
        );
    let fixed = QueryPoint::new()
        .with_quantity("P", pressure, 101_325.0)
        .unwrap();
    let (validity, fixed) = match temperature_quantity {
        Some(quantity) => (
            environment
                .clone()
                .with_quantity("T", quantity, 300.0, 400.0),
            fixed.with_quantity("T", quantity, 300.0).unwrap(),
        ),
        None => (
            environment.clone().with("T", 300.0, 400.0),
            fixed.with("T", 300.0).unwrap(),
        ),
    };
    let provenance = Provenance {
        source: "synthetic temperature/depth hardness; no measured qualification".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    };
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "synthetic graded specimen at declared temperature".into(),
            method: "synthetic force/projected-area indentation".into(),
            artifact: fs_blake3::hash_domain("fixture", b"thermal depth hardness"),
            caveats: "numerical references only".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    let key = PropertyKey::with_quantity("H", INDENTATION_HARDNESS_QUANTITY)
        .with_hardness_test(
            HardnessTestContext::new(
                "synthetic sphere radius 1 mm",
                vec![
                    HardnessLoadStep::new(
                        QtyAny::new(10.0, Dims([1, 1, -2, 0, 0, 0])),
                        QtyAny::new(1.0, Dims([0, 0, 1, 0, 0, 0])),
                    )
                    .unwrap(),
                ],
                "synthetic thermal depth protocol revision 1",
                observation,
            )
            .unwrap(),
        )
        .unwrap();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value,
            validity,
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![observation],
            provenance: provenance.clone(),
        })
        .unwrap();
    if competing_hot_claim {
        let validity = match temperature_quantity {
            Some(quantity) => environment.with_quantity("T", quantity, 350.0, 450.0),
            None => environment.with("T", 350.0, 450.0),
        };
        claims
            .insert_claim(PropertyClaim {
                key: key.clone(),
                value: PropertyValue::Scalar {
                    value: 10_000.0,
                    dims: Pressure::DIMS,
                },
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![observation],
                provenance,
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermal graded solid".into(),
            phase: "solid".into(),
            process: "synthetic hardness depth/temperature observations".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let profile = ResolvedDepthHardness::new(
        &card,
        &fixed,
        key,
        MaterialPropertySelection::SingleClaimOnly,
        1.0,
    )
    .unwrap();
    (card, profile)
}

#[test]
fn original_thermal_profile_refresh_preserves_test_context_and_exact_depth_integrals() {
    for convention in [
        None,
        Some(QuantitySpec::dimensional(Temperature::DIMS)),
        Some(absolute_temperature()),
    ] {
        let (card, cold) = thermal_profile(thermal_surface(false, false), convention, true);
        let before = cold.clone();
        let noop = cold.at_temperature(&card, 300.0).unwrap();
        assert_eq!(noop, cold);
        let original = cold.reference_source().resolved();
        let fixed = original
            .query_point()
            .iter()
            .filter(|(axis, _)| axis != Z)
            .try_fold(QueryPoint::new(), |point, (axis, value)| {
                match original.axis_quantities().get(axis) {
                    Some(quantity) => point.with_quantity(axis, *quantity, *value),
                    None => point.with(axis, *value),
                }
            })
            .unwrap();
        let pinned_cold = ResolvedDepthHardness::new(
            &card,
            &fixed,
            cold.reference_source().key().clone(),
            MaterialPropertySelection::PinnedByProperty(vec![(
                cold.reference_source().key().name().to_owned(),
                cold.reference_source().property().answer().receipt.selected,
            )]),
            cold.maximum_depth_m(),
        )
        .unwrap();
        for temperature in [325.0, 350.0, 400.0] {
            let hot = cold.at_temperature(&card, temperature).unwrap();
            assert_eq!(hot.maximum_depth_m(), cold.maximum_depth_m());
            assert_ne!(hot.identity(), cold.identity());
            assert_eq!(hot.reference_source().key(), cold.reference_source().key());
            assert_eq!(
                hot.reference_source().property().answer().receipt.selected,
                cold.reference_source().property().answer().receipt.selected
            );
            assert_eq!(
                hot.reference_source().resolved().axis_quantities(),
                cold.reference_source().resolved().axis_quantities()
            );
            assert_eq!(
                hot.reference_source().resolved().state_coordinate("T"),
                Some(temperature)
            );
            assert_eq!(
                hot.reference_source().resolved().state_coordinate("P"),
                Some(101_325.0)
            );
            assert_eq!(
                hot.reference_source().resolved().state_coordinate(Z),
                Some(0.0)
            );
            assert_eq!(
                hot.reference_source()
                    .resolved()
                    .query_point()
                    .iter()
                    .filter(|(axis, _)| axis != "T")
                    .collect::<Vec<_>>(),
                cold.reference_source()
                    .resolved()
                    .query_point()
                    .iter()
                    .filter(|(axis, _)| axis != "T")
                    .collect::<Vec<_>>()
            );
            let base = 900.0 - 2.0 * (temperature - 300.0);
            for depth in [0.0, 0.2, 0.4, 0.55, 0.8, 1.0] {
                close(
                    hot.hardness_pa_at(depth).unwrap(),
                    base + 400.0 * depth,
                    2e-15,
                );
                let typed = hot.resolve_at_depth(&card, depth).unwrap();
                close(typed.hardness_pa(), base + 400.0 * depth, 2e-15);
                assert_eq!(typed.key(), cold.reference_source().key());
            }
            let expected = base * (0.8 - 0.2) + 200.0 * (0.8_f64.powi(2) - 0.2_f64.powi(2));
            close(
                hot.hardness_integral_between(0.2, 0.8).unwrap(),
                expected,
                3e-15,
            );
            close(
                hot.advance_by_integral(0.2, expected).unwrap().depth_m,
                0.8,
                3e-15,
            );
            assert!(
                hot.advance_by_integral(0.2, 100.0).unwrap().depth_m
                    > cold.advance_by_integral(0.2, 100.0).unwrap().depth_m
            );
            // A real thermal revisit truthfully records explicit source pinning;
            // compare with the same original cold profile prepared with that pin.
            assert_eq!(hot.at_temperature(&card, 300.0).unwrap(), pinned_cold);
        }
        assert_eq!(cold, before);
    }
}

#[test]
fn thermal_refresh_checks_positive_continuous_full_depth_support_not_only_endpoints() {
    for (negative_interior, hole) in [(true, false), (false, true)] {
        let (card, cold) = thermal_profile(
            thermal_surface(negative_interior, hole),
            Some(absolute_temperature()),
            true,
        );
        let before = cold.clone();
        for depth in [0.0, 1.0] {
            let query = cold
                .reference_source()
                .resolved()
                .query()
                .with("T", 400.0)
                .unwrap()
                .with(Z, depth)
                .unwrap();
            let endpoint = cold
                .reference_source()
                .at_query_point(&card, &query)
                .unwrap();
            assert!(endpoint.hardness_pa() > 0.0);
        }
        let refused = cold.at_temperature(&card, 400.0);
        if hole {
            assert!(matches!(refused, Err(DepthHardnessError::Support(_))));
        } else {
            assert!(matches!(refused, Err(DepthHardnessError::Invalid(_))));
        }
        // The competing positive hot claim cannot fill a hole or replace a
        // nonpositive interior of the exact original indentation profile.
        assert_eq!(cold, before);
        close(
            cold.hardness_integral_between(0.0, 1.0).unwrap(),
            1100.0,
            2e-15,
        );
    }
}

#[test]
fn temperature_refuses_changed_cards_unowned_temperature_and_outside_original_claim() {
    let (card, cold) = thermal_profile(
        thermal_surface(false, false),
        Some(absolute_temperature()),
        true,
    );
    let before = cold.clone();
    let foreign = MaterialCard::supersede(&card, card.claims().clone(), vec![]).unwrap();
    for temperature in [300.0, 350.0] {
        assert!(cold.at_temperature(&foreign, temperature).is_err());
    }
    for temperature in [0.0, -1.0, f64::NAN, f64::INFINITY, 299.0, 450.0] {
        assert!(cold.at_temperature(&card, temperature).is_err());
    }
    assert_eq!(cold, before);
    let no_temperature = profile(vec![(0.0, 100.0), (1.0, 200.0)]);
    assert!(matches!(
        no_temperature.at_temperature(&card, 300.0),
        Err(DepthHardnessError::Invalid(_))
    ));
    let pressure =
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
    let (misnamed_card, misnamed) = thermal_profile(
        PropertyValue::Curve {
            abscissa: Z.into(),
            abscissa_dims: Dims([1, 0, 0, 0, 0, 0]),
            knots: vec![(0.0, 100.0), (1.0, 200.0)],
            dims: Pressure::DIMS,
        },
        Some(pressure),
        false,
    );
    assert!(matches!(
        misnamed.at_temperature(&misnamed_card, 350.0),
        Err(DepthHardnessError::Invalid(_))
    ));
}
