//! Native pump spectra and thermal refresh retain original source decisions.

use super::*;
use fs_material::phosphorescence::pump::*;

const WAVELENGTH: QuantitySpec = QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0]));
const ABSOLUTE_TEMPERATURE: QuantitySpec = QuantitySpec::semantic(SemanticType::new(
    QuantityKind::AbsoluteTemperature,
    ValueForm::Static,
));

fn point(wavelength_m: f64, temperature_k: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("lambda", WAVELENGTH, wavelength_m)
        .unwrap()
        .with_quantity("T", ABSOLUTE_TEMPERATURE, temperature_k)
        .unwrap()
}

fn pump_claims() -> Vec<PropertyClaim> {
    GROUND_STATE_PUMP_CROSS_SECTION_PROPERTIES
        .into_iter()
        .zip([[1e-20, 2e-20, 2e-20, 3e-20], [2e-21, 4e-21, 5e-21, 7e-21]])
        .map(|(name, values)| PropertyClaim {
            key: PropertyKey::new(name, GROUND_STATE_PUMP_CROSS_SECTION_DIMS),
            value: PropertyValue::Surface2D {
                axes: ["lambda".into(), "T".into()],
                axis_dims: [WAVELENGTH.dims(), Temperature::DIMS],
                samples: [
                    [400e-9, 300.0],
                    [600e-9, 300.0],
                    [400e-9, 600.0],
                    [600e-9, 600.0],
                ]
                .into_iter()
                .zip(values)
                .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: GROUND_STATE_PUMP_CROSS_SECTION_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("lambda", WAVELENGTH, 400e-9, 600e-9)
                .with_quantity("T", ABSOLUTE_TEMPERATURE, 300.0, 600.0),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: Vec::new(),
            provenance: Provenance {
                source: format!("manufactured molecular pump surface: {name}"),
                license: "internal-test-use".into(),
                artifact: None,
            },
        })
        .collect()
}

fn resolve_pump(card: &MaterialCard) -> GroundStatePumpStatePoint {
    GroundStatePumpStatePoint::resolve(
        card,
        &point(400e-9, 300.0),
        MaterialPropertySelection::SingleClaimOnly,
        "lambda",
        Some("T"),
    )
    .unwrap()
}

#[test]
fn g1_ground_state_pump_queries_exact_native_wavelength_temperature_surfaces() {
    let card = assemble(pump_claims());
    let cold = resolve_pump(&card);
    for (wavelength, expected) in [(450e-9, [1.75e-20, 4e-21]), (550e-9, [2.25e-20, 5e-21])] {
        let warm = cold
            .at_wavelength_temperature(&card, wavelength, Some(450.0))
            .unwrap();
        for (actual, expected) in warm.cross_sections_m2().into_iter().zip(expected) {
            close(actual, expected);
        }
        assert_eq!(warm.wavelength_m(), wavelength);
        assert_eq!(warm.temperature_k(), Some(450.0));
        assert_ne!(warm.identity(), cold.identity());
        for (a, b) in cold
            .resolved()
            .properties()
            .iter()
            .zip(warm.resolved().properties())
        {
            assert_eq!(a.requirement(), b.requirement());
            assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
        }
    }
    assert_eq!(cold.cross_sections_m2(), [1e-20, 2e-21]);
    assert_eq!(cold.wavelength_m(), 400e-9);
    assert_eq!(cold.temperature_k(), Some(300.0));
    assert_eq!(cold.wavelength_axis(), "lambda");
    assert_eq!(cold.temperature_axis(), Some("T"));
}

#[test]
fn g0_ground_state_pump_keeps_complete_bundle_and_pins_claims_when_support_overlaps() {
    let mut input = pump_claims();
    let mut competing = input[0].clone();
    competing.value = PropertyValue::Scalar {
        value: 9e-20,
        dims: GROUND_STATE_PUMP_CROSS_SECTION_DIMS,
    };
    competing.validity = competing
        .validity
        .with_quantity("T", ABSOLUTE_TEMPERATURE, 450.0, 600.0);
    competing.interpolation = InterpolationPolicy::ConstantWithinValidity;
    competing.provenance.source.push_str(" alternative source");
    input.push(competing);
    let mut unrelated = input[0].clone();
    unrelated.key = PropertyKey::new("density", Dims([-3, 1, 0, 0, 0, 0]));
    unrelated.value = PropertyValue::Scalar {
        value: 1200.0,
        dims: unrelated.key.dims(),
    };
    unrelated.interpolation = InterpolationPolicy::ConstantWithinValidity;
    input.push(unrelated);
    let card = assemble(input);
    let original = resolve_pump(&card);
    let expanded = original
        .resolved()
        .with_additional_requirements(
            &card,
            &[ScalarPropertyRequirement::try_new(
                "density",
                Dims([-3, 1, 0, 0, 0, 0]),
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap()],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
    let complete =
        GroundStatePumpStatePoint::try_from_material(expanded.clone(), "lambda", Some("T"))
            .unwrap();
    assert_eq!(complete.resolved(), &expanded);
    let warm = complete
        .at_wavelength_temperature(&card, 500e-9, Some(450.0))
        .unwrap();
    close(warm.cross_sections_m2()[0], 2e-20);
    close(warm.cross_sections_m2()[1], 4.5e-21);
    assert_eq!(
        warm.resolved().property("density").unwrap().value_si(),
        1200.0
    );
    assert!(
        GroundStatePumpStatePoint::resolve(
            &card,
            &point(500e-9, 450.0),
            MaterialPropertySelection::SingleClaimOnly,
            "lambda",
            Some("T"),
        )
        .is_err()
    );
}

#[test]
fn g0_ground_state_pump_refuses_source_escape_and_unowned_or_mistyped_axes() {
    let card = assemble(pump_claims());
    let source = resolve_pump(&card);
    for (wavelength, temperature) in [
        (399e-9, Some(300.0)),
        (601e-9, Some(300.0)),
        (500e-9, Some(299.0)),
        (500e-9, Some(601.0)),
        (0.0, Some(300.0)),
        (f64::NAN, Some(300.0)),
        (500e-9, None),
    ] {
        assert!(
            source
                .at_wavelength_temperature(&card, wavelength, temperature)
                .is_err()
        );
    }
    for point in [
        QueryPoint::new()
            .with("lambda", 500e-9)
            .unwrap()
            .with_quantity("T", ABSOLUTE_TEMPERATURE, 300.0)
            .unwrap(),
        QueryPoint::new()
            .with_quantity("lambda", WAVELENGTH, 500e-9)
            .unwrap()
            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
            .unwrap(),
    ] {
        assert!(
            GroundStatePumpStatePoint::resolve(
                &card,
                &point,
                MaterialPropertySelection::SingleClaimOnly,
                "lambda",
                Some("T"),
            )
            .is_err()
        );
    }
    let mut different = pump_claims();
    different[0].provenance.source.push_str(" changed card");
    assert!(
        source
            .at_wavelength_temperature(&assemble(different), 500e-9, Some(300.0))
            .is_err()
    );

    let mut input = pump_claims();
    let elapsed = QuantitySpec::dimensional(Dims([0, 0, 1, 0, 0, 0]));
    input[0].validity = input[0]
        .validity
        .clone()
        .with_quantity("elapsed", elapsed, 0.0, 10.0);
    let history = point(500e-9, 300.0)
        .with_quantity("elapsed", elapsed, 1.0)
        .unwrap();
    assert!(
        GroundStatePumpStatePoint::resolve(
            &assemble(input),
            &history,
            MaterialPropertySelection::SingleClaimOnly,
            "lambda",
            Some("T"),
        )
        .is_err()
    );
    assert_eq!(source, resolve_pump(&card));
}

#[test]
fn g0_ground_state_pump_rechecks_permissive_bundles_and_preserves_explicit_zero() {
    for wrong_units in [false, true] {
        let mut input = pump_claims();
        let dims = if wrong_units {
            Dims([2, -1, 0, 0, 0, 0])
        } else {
            GROUND_STATE_PUMP_CROSS_SECTION_DIMS
        };
        input[1].key = PropertyKey::new(GROUND_STATE_PUMP_CROSS_SECTION_PROPERTIES[1], dims);
        input[1].value = PropertyValue::Scalar {
            value: -1e-21,
            dims,
        };
        input[1].interpolation = InterpolationPolicy::ConstantWithinValidity;
        let requirements: Vec<_> = input
            .iter()
            .map(|claim| {
                ScalarPropertyRequirement::try_with_key(&claim.key, ScalarAdmissibility::Finite)
                    .unwrap()
            })
            .collect();
        let resolved = resolve_material_state_point(
            &assemble(input),
            &point(500e-9, 300.0),
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert!(
            GroundStatePumpStatePoint::try_from_material(resolved, "lambda", Some("T")).is_err()
        );
    }
    let mut input = pump_claims();
    input[1].value = PropertyValue::Scalar {
        value: 0.0,
        dims: GROUND_STATE_PUMP_CROSS_SECTION_DIMS,
    };
    input[1].interpolation = InterpolationPolicy::ConstantWithinValidity;
    let card = assemble(input);
    let unowned_temperature = GroundStatePumpStatePoint::resolve(
        &card,
        &point(500e-9, 300.0),
        MaterialPropertySelection::SingleClaimOnly,
        "lambda",
        None,
    )
    .unwrap();
    assert_eq!(unowned_temperature.cross_sections_m2()[1], 0.0);
    assert_eq!(unowned_temperature.temperature_k(), None);
    assert!(
        unowned_temperature
            .at_wavelength_temperature(&card, 550e-9, Some(300.0))
            .is_err()
    );
    let refreshed = unowned_temperature
        .at_wavelength_temperature(&card, 550e-9, None)
        .unwrap();
    assert_eq!(refreshed.resolved().state_coordinate("T"), Some(300.0));
    close(refreshed.cross_sections_m2()[0], 1.75e-20);
}
