//! Original binary Maxwell–Stefan source kinetics; all data are manufactured.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimId, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    liquid::solution::{
        BinaryLiquidError, MaterialBinaryDiffusivity,
        diffusion::BINARY_LIQUID_DIFFUSIVITY_DIMS as DIMS,
    },
    state_point::MaterialPropertySelection,
};
use fs_qty::{
    Dimensionless, QuantitySpec, Temperature,
    semantic::{CompositionBasis, QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::SpeciesId;
use std::sync::Arc;

fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn species() -> [SpeciesId; 2] {
    [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()]
}
fn property() -> PropertyKey {
    PropertyKey::new("binary_ms_diffusivity", DIMS)
}
fn query(temperature: f64, fraction: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            typed(QuantityKind::AbsoluteTemperature),
            temperature,
        )
        .unwrap()
        .with_quantity(
            "x_a",
            typed(QuantityKind::Composition(CompositionBasis::MoleFraction)),
            fraction,
        )
        .unwrap()
        .with_quantity("pressure", typed(QuantityKind::Pressure), 100_000.0)
        .unwrap()
}
fn claim(value: PropertyValue) -> PropertyClaim {
    PropertyClaim {
        key: property(),
        interpolation: if matches!(&value, PropertyValue::Scalar { .. }) {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        },
        value,
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                "temperature",
                typed(QuantityKind::AbsoluteTemperature),
                300.0,
                400.0,
            )
            .with_quantity(
                "x_a",
                typed(QuantityKind::Composition(CompositionBasis::MoleFraction)),
                0.0,
                1.0,
            )
            .with_quantity(
                "pressure",
                typed(QuantityKind::Pressure),
                90_000.0,
                110_000.0,
            ),
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured binary Maxwell-Stefan kinetic data".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn card(sources: Vec<PropertyClaim>) -> (Arc<MaterialCard>, Vec<ClaimId>) {
    let mut set = ClaimSet::new();
    let ids = sources
        .into_iter()
        .map(|source| set.insert_claim(source).unwrap())
        .collect();
    let material = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured binary A+B liquid".into(),
            phase: "liquid".into(),
            process: "source-kinetics fixture".into(),
            revision: 0,
        },
        set,
        vec![],
    )
    .unwrap();
    (Arc::new(material), ids)
}
fn bind(material: Arc<MaterialCard>, reference: QueryPoint) -> MaterialBinaryDiffusivity {
    MaterialBinaryDiffusivity::new(
        species(),
        material,
        property(),
        MaterialPropertySelection::SingleClaimOnly,
        reference,
        "temperature",
        "x_a",
    )
    .unwrap()
}
fn curve(values: [f64; 3]) -> PropertyValue {
    PropertyValue::Curve {
        abscissa: "x_a".into(),
        abscissa_dims: Dimensionless::DIMS,
        knots: vec![(0.0, values[0]), (0.5, values[1]), (1.0, values[2])],
        dims: DIMS,
    }
}
fn surface() -> PropertyValue {
    PropertyValue::Surface2D {
        axes: ["temperature".into(), "x_a".into()],
        axis_dims: [Temperature::DIMS, Dimensionless::DIMS],
        samples: vec![
            ([300.0, 0.0], 1e-9),
            ([300.0, 0.5], 4e-9),
            ([300.0, 1.0], 2e-9),
            ([400.0, 0.0], 3e-9),
            ([400.0, 0.5], 8e-9),
            ([400.0, 1.0], 5e-9),
        ],
        triangles: vec![[0, 3, 4], [0, 4, 1], [1, 4, 5], [1, 5, 2]],
        dims: DIMS,
    }
}
fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() <= 2e-14 * expected.abs().max(1e-30));
}

#[test]
fn g1_original_surface_tracks_temperature_and_composition_without_averaging_diffusivities() {
    let (material, ids) = card(vec![claim(surface())]);
    let law = bind(material.clone(), query(300.0, 0.5));
    let endpoints = [[0.25, 0.75], [0.75, 0.25]];
    let cold = law
        .resolve_face(Temperature::new(300.0), endpoints, || false)
        .unwrap();
    close(cold.left.properties()[0].value_si(), 2.5e-9);
    close(cold.right.properties()[0].value_si(), 3e-9);
    close(cold.diffusivity_m2_s, 4e-9);
    assert!(cold.diffusivity_m2_s > 0.5 * (2.5e-9 + 3e-9));
    let warm = law
        .resolve_face(Temperature::new(400.0), endpoints, || false)
        .unwrap();
    close(warm.left.properties()[0].value_si(), 5.5e-9);
    close(warm.right.properties()[0].value_si(), 6.5e-9);
    close(warm.diffusivity_m2_s, 8e-9);
    let shifted = law
        .resolve_face(Temperature::new(300.0), [[0.0, 1.0], [0.5, 0.5]], || false)
        .unwrap();
    close(shifted.diffusivity_m2_s, 2.5e-9);
    for result in [&cold, &warm, &shifted] {
        for receipt in [&result.left, &result.face, &result.right] {
            assert_eq!(receipt.card_identity(), material.content_hash());
            assert_eq!(receipt.state_coordinate("pressure"), Some(100_000.0));
            assert_eq!(receipt.properties()[0].answer().receipt.selected, ids[0]);
        }
    }
    assert_eq!(cold.face.state_coordinate("x_a"), Some(0.5));
    assert_eq!(warm.face.state_coordinate("temperature"), Some(400.0));
    assert_eq!(
        law.resolve_face(Temperature::new(300.0), endpoints, || false)
            .unwrap(),
        cold
    );
    assert_eq!(law.selected_claim(), ids[0]);
    assert_eq!(law.species(), &species());
}

#[test]
fn g0_zero_is_an_explicit_closed_face_and_negative_endpoint_kinetics_refuse() {
    let (zero, _) = card(vec![claim(curve([1e-9, 0.0, 1e-9]))]);
    let zero = bind(zero, query(300.0, 0.5));
    let endpoints = [[0.0, 1.0], [1.0, 0.0]];
    let resolved = zero
        .resolve_face(Temperature::new(300.0), endpoints, || false)
        .unwrap();
    assert_eq!(resolved.diffusivity_m2_s, 0.0);
    assert!(resolved.left.properties()[0].value_si() > 0.0);
    assert!(resolved.right.properties()[0].value_si() > 0.0);

    let (negative, _) = card(vec![claim(curve([-1e-9, 1e-9, 1e-9]))]);
    let negative = bind(negative, query(300.0, 0.5));
    assert!(matches!(
        negative.resolve_face(Temperature::new(300.0), endpoints, || false),
        Err(BinaryLiquidError::Source(_))
    ));
}

#[test]
fn g0_endpoint_and_joint_face_support_both_remain_authoritative() {
    let mut bounded = claim(curve([1e-9, 2e-9, 3e-9]));
    bounded.validity = bounded.validity.with_quantity(
        "x_a",
        typed(QuantityKind::Composition(CompositionBasis::MoleFraction)),
        0.1,
        0.9,
    );
    let (bounded, _) = card(vec![bounded]);
    let bounded = bind(bounded, query(300.0, 0.5));
    // The mean is in range, but the left endpoint is not.
    assert!(
        bounded
            .resolve_face(Temperature::new(300.0), [[0.0, 1.0], [0.8, 0.2]], || false)
            .is_err()
    );

    let (isolated, _) = card(vec![claim(PropertyValue::Surface2D {
        axes: ["temperature".into(), "x_a".into()],
        axis_dims: [Temperature::DIMS, Dimensionless::DIMS],
        samples: vec![([300.0, 0.25], 1e-9), ([300.0, 0.75], 2e-9)],
        triangles: vec![],
        dims: DIMS,
    })]);
    let isolated = bind(isolated, query(300.0, 0.25));
    // Both endpoints are original samples; marginal axis bounds cannot invent
    // interpolation support for their absent joint face state.
    assert!(
        isolated
            .resolve_face(
                Temperature::new(300.0),
                [[0.25, 0.75], [0.75, 0.25]],
                || false
            )
            .is_err()
    );
}

#[test]
fn g0_temperature_updates_cannot_switch_to_another_in_domain_claim() {
    let mut first = claim(curve([1e-9, 2e-9, 3e-9]));
    first.validity = first.validity.with_quantity(
        "temperature",
        typed(QuantityKind::AbsoluteTemperature),
        300.0,
        350.0,
    );
    let mut second = claim(curve([4e-9, 5e-9, 6e-9]));
    second.validity = second.validity.with_quantity(
        "temperature",
        typed(QuantityKind::AbsoluteTemperature),
        350.0,
        400.0,
    );
    let (material, ids) = card(vec![first, second]);
    let cold = bind(material.clone(), query(325.0, 0.5));
    let warm = bind(material, query(375.0, 0.5));
    assert_eq!(cold.selected_claim(), ids[0]);
    assert_eq!(warm.selected_claim(), ids[1]);
    let endpoints = [[0.25, 0.75], [0.75, 0.25]];
    assert!(
        cold.resolve_face(Temperature::new(375.0), endpoints, || false)
            .is_err()
    );
    assert!(
        warm.resolve_face(Temperature::new(325.0), endpoints, || false)
            .is_err()
    );
}

#[test]
fn g0_typed_coordinates_and_explicit_binary_identity_are_required() {
    let (material, _) = card(vec![claim(surface())]);
    let wrong_points = [
        QueryPoint::new()
            .with("temperature", 300.0)
            .unwrap()
            .with_quantity(
                "x_a",
                typed(QuantityKind::Composition(CompositionBasis::MoleFraction)),
                0.5,
            )
            .unwrap(),
        query(300.0, 0.5)
            .with_quantity(
                "temperature",
                typed(QuantityKind::TemperatureDifference),
                300.0,
            )
            .unwrap(),
        query(300.0, 0.5)
            .with_quantity(
                "x_a",
                typed(QuantityKind::Composition(CompositionBasis::MassFraction)),
                0.5,
            )
            .unwrap(),
        query(300.0, 0.5)
            .with_quantity(
                "other_temperature",
                typed(QuantityKind::AbsoluteTemperature),
                300.0,
            )
            .unwrap(),
        query(300.0, 0.5)
            .with_quantity("stress", typed(QuantityKind::Stress), 1e5)
            .unwrap(),
    ];
    for point in wrong_points {
        assert!(
            MaterialBinaryDiffusivity::new(
                species(),
                material.clone(),
                property(),
                MaterialPropertySelection::SingleClaimOnly,
                point,
                "temperature",
                "x_a",
            )
            .is_err()
        );
    }
    let duplicate = [species()[0].clone(), species()[0].clone()];
    assert!(
        MaterialBinaryDiffusivity::new(
            duplicate,
            material.clone(),
            property(),
            MaterialPropertySelection::SingleClaimOnly,
            query(300.0, 0.5),
            "temperature",
            "x_a",
        )
        .is_err()
    );
    assert!(
        MaterialBinaryDiffusivity::new(
            species(),
            material,
            PropertyKey::new("binary_ms_diffusivity", Dimensionless::DIMS),
            MaterialPropertySelection::SingleClaimOnly,
            query(300.0, 0.5),
            "temperature",
            "x_a",
        )
        .is_err()
    );
}

#[test]
fn g0_binding_equality_retains_ordered_species_original_card_and_held_pressure() {
    let (material, _) = card(vec![claim(surface())]);
    let a = bind(material.clone(), query(300.0, 0.5));
    let (reconstructed, _) = card(vec![claim(surface())]);
    assert_eq!(a, bind(reconstructed, query(300.0, 0.5)));
    assert_ne!(
        a,
        bind(
            material.clone(),
            query(300.0, 0.5).with("pressure", 101_000.0).unwrap()
        )
    );
    let mut reversed = species();
    reversed.swap(0, 1);
    let reversed = MaterialBinaryDiffusivity::new(
        reversed,
        material,
        property(),
        MaterialPropertySelection::SingleClaimOnly,
        query(300.0, 0.5),
        "temperature",
        "x_a",
    )
    .unwrap();
    assert_ne!(a, reversed);
}

#[test]
fn g4_invalid_simplex_temperature_and_cancellation_publish_no_partial_face() {
    let (material, _) = card(vec![claim(curve([1e-9; 3]))]);
    let law = bind(material, query(300.0, 0.5));
    let before = law.clone();
    let endpoints = [[0.25, 0.75], [0.75, 0.25]];
    let mut polls = 0;
    let expected = law
        .resolve_face(Temperature::new(300.0), endpoints, || {
            polls += 1;
            false
        })
        .unwrap();
    assert_eq!(expected.diffusivity_m2_s, 1e-9);
    for stop in 1..=polls {
        let mut count = 0;
        assert_eq!(
            law.resolve_face(Temperature::new(300.0), endpoints, || {
                count += 1;
                count == stop
            }),
            Err(BinaryLiquidError::Cancelled)
        );
    }
    for endpoint in [
        [0.8, 0.4],
        [-0.1, 1.1],
        [f64::NAN, 0.5],
        [0.0, 0.0],
        [1.0, 1e-20],
    ] {
        assert!(
            law.resolve_face(Temperature::new(300.0), [endpoint, endpoints[1]], || false)
                .is_err()
        );
    }
    for temperature in [0.0, -300.0, f64::INFINITY, f64::NAN, 299.0, 401.0] {
        assert!(
            law.resolve_face(Temperature::new(temperature), endpoints, || false)
                .is_err()
        );
    }
    let trace = law
        .resolve_face(Temperature::new(300.0), [[1e-20, 1.0]; 2], || false)
        .unwrap();
    assert_eq!(trace.face.state_coordinate("x_a"), Some(1e-20));
    assert_eq!(law, before);
    assert_eq!(
        law.resolve_face(Temperature::new(300.0), endpoints, || false)
            .unwrap(),
        expected
    );
}
