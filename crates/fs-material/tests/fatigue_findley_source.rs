//! Original Findley fits, typed applicability and native material binding.
//! Synthetic coefficients exercise the implementation; they are not measured data.
#![deny(unsafe_code)]

use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimId, ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    ObservationDataset, ObservationId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::fatigue::critical_plane::findley::source::{
    FindleyProperty, FindleySourceRequest, resolve_findley,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn request() -> FindleySourceRequest {
    FindleySourceRequest {
        temperature_axis: "specimen-temperature".into(),
        surface_and_geometry: "synthetic smooth torsion coupon".into(),
        environment: "synthetic dry atmosphere".into(),
        load_basis: "synthetic common axial/torsional Findley calibration".into(),
    }
}

fn point(request: &FindleySourceRequest, temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            &request.temperature_axis,
            absolute_temperature(),
            temperature,
        )
        .unwrap()
}

fn rows(request: &FindleySourceRequest) -> Vec<PropertyClaim> {
    FindleyProperty::ALL
        .into_iter()
        .map(|property| {
            let key = property.property_key();
            PropertyClaim {
                value: PropertyValue::Scalar {
                    value: [0.2, 200e6][property.index()],
                    dims: key.dims(),
                },
                key,
                validity: ValidityDomain::unconstrained().with_quantity(
                    &request.temperature_axis,
                    absolute_temperature(),
                    300.0,
                    400.0,
                ),
                uncertainty: UncertaintyModel::RelativeHalfWidth {
                    fraction: 0.05,
                    confidence: 0.9,
                },
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic common Findley fit".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            }
        })
        .collect()
}

fn card_with_observations(
    rows: Vec<PropertyClaim>,
    observations: Vec<ObservationDataset>,
) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for observation in observations {
        claims.register_observation(observation).unwrap();
    }
    for row in rows {
        claims.insert_claim(row).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic alloy".into(),
            phase: "solid".into(),
            process: "synthetic coupon process".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn card(rows: Vec<PropertyClaim>) -> MaterialCard {
    card_with_observations(rows, vec![])
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-12 * expected.abs().max(1e-10),
        "{actual:e} vs {expected:e}"
    );
}

#[test]
fn original_findley_claims_uncertainty_and_receipts_survive() {
    let request = request();
    let original = rows(&request);
    let material = card(original.clone());
    let at = point(&request, 300.0);
    let calibration = resolve_findley(&material, &at, &request, [SINGLE; 2]).unwrap();
    close(calibration.normal_stress_sensitivity(), 0.2);
    close(calibration.limit_factor_pa(), 200e6);
    assert_eq!(
        FindleyProperty::NormalStressSensitivity
            .property_key()
            .dims(),
        Dims::NONE
    );
    assert_eq!(
        FindleyProperty::LimitFactor.property_key().dims(),
        Pressure::DIMS
    );
    assert_ne!(
        FindleyProperty::NormalStressSensitivity.property_key(),
        FindleyProperty::LimitFactor.property_key()
    );
    let source = calibration.source().unwrap();
    assert_eq!(source.card_identity(), material.content_hash());
    assert_eq!(source.material_state(), material.id());
    assert_eq!(source.request(), &request);
    assert_eq!(source.query(), &at);
    assert_eq!(source.selections(), &[SINGLE; 2]);
    for property in FindleyProperty::ALL {
        assert_eq!(source.claim(property), &original[property.index()]);
        let answer = material
            .claims()
            .query_typed(
                &property.property_key(),
                &at,
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
        assert_eq!(source.property(property), &answer);
        material
            .claims()
            .verify_receipt(&source.property(property).receipt)
            .unwrap();
    }
    assert_eq!(
        calibration.scope().source_reference,
        original[0].provenance.source
    );
    assert_eq!(calibration.scope().material, material.id().chemistry);
    assert_eq!(calibration.scope().process, material.id().process);
    assert_eq!(calibration.scope().temperature_k, 300.0);
    assert_eq!(
        calibration.scope().surface_and_geometry,
        request.surface_and_geometry
    );
    assert_eq!(calibration.scope().environment, request.environment);
    assert_eq!(calibration.scope().load_basis, request.load_basis);
}

#[test]
fn source_temperature_curves_keep_their_knots_and_interpolation() {
    let request = request();
    let mut original = rows(&request);
    for property in FindleyProperty::ALL {
        let row = &mut original[property.index()];
        row.value = PropertyValue::Curve {
            abscissa: request.temperature_axis.clone(),
            abscissa_dims: Temperature::DIMS,
            knots: if property == FindleyProperty::NormalStressSensitivity {
                vec![(300.0, 0.2), (400.0, 0.4)]
            } else {
                vec![(300.0, 200e6), (400.0, 100e6)]
            },
            dims: row.key.dims(),
        };
        row.interpolation = InterpolationPolicy::LinearInside;
    }
    let material = card(original.clone());
    let calibration =
        resolve_findley(&material, &point(&request, 350.0), &request, [SINGLE; 2]).unwrap();
    close(calibration.normal_stress_sensitivity(), 0.3);
    close(calibration.limit_factor_pa(), 150e6);
    for property in FindleyProperty::ALL {
        let source = calibration.source().unwrap();
        assert_eq!(source.claim(property), &original[property.index()]);
        material
            .claims()
            .verify_receipt(&source.property(property).receipt)
            .unwrap();
    }
    for outside in [299.0, 401.0] {
        assert!(
            resolve_findley(&material, &point(&request, outside), &request, [SINGLE; 2]).is_err()
        );
    }
}

#[test]
fn missing_parameters_are_not_inferred_and_ambiguity_requires_pins() {
    let request = request();
    let at = point(&request, 300.0);
    for missing in FindleyProperty::ALL {
        let mut incomplete = rows(&request);
        let mut ordinary_strength = incomplete.remove(missing.index());
        ordinary_strength.key = PropertyKey::new("ultimate-tensile-strength", Pressure::DIMS);
        ordinary_strength.value = PropertyValue::Scalar {
            value: 500e6,
            dims: Pressure::DIMS,
        };
        incomplete.push(ordinary_strength);
        assert!(resolve_findley(&card(incomplete), &at, &request, [SINGLE; 2]).is_err());
    }
    for ambiguous in FindleyProperty::ALL {
        let mut original = rows(&request);
        let keep = ClaimId(original[ambiguous.index()].content_hash());
        let mut alternative = original[ambiguous.index()].clone();
        alternative.value = PropertyValue::Scalar {
            value: [0.3, 250e6][ambiguous.index()],
            dims: alternative.key.dims(),
        };
        original.push(alternative);
        let material = card(original);
        assert!(resolve_findley(&material, &at, &request, [SINGLE; 2]).is_err());
        let mut selections = [SINGLE; 2];
        selections[ambiguous.index()] = ClaimSelection::Pinned(keep);
        let calibration = resolve_findley(&material, &at, &request, selections).unwrap();
        close(calibration.normal_stress_sensitivity(), 0.2);
        close(calibration.limit_factor_pa(), 200e6);
        assert_eq!(
            calibration
                .source()
                .unwrap()
                .property(ambiguous)
                .receipt
                .selected,
            keep
        );
    }
}

#[test]
fn parameters_from_distinct_source_fits_or_observation_populations_refuse() {
    let request = request();
    let at = point(&request, 300.0);
    let mut mixed = rows(&request);
    mixed[1].provenance.source = "another synthetic Findley fit".into();
    assert!(resolve_findley(&card(mixed), &at, &request, [SINGLE; 2]).is_err());

    let mut original = rows(&request);
    let observations: Vec<_> = [1_u8, 2]
        .into_iter()
        .map(|number| ObservationDataset {
            specimen: format!("synthetic population {number}"),
            method: "synthetic Findley fixture; not a measured test".into(),
            artifact: ContentHash([number; 32]),
            caveats: "test-only".into(),
            provenance: original[0].provenance.clone(),
        })
        .collect();
    original[0].observations = vec![ObservationId(observations[0].content_hash())];
    original[1].observations = vec![ObservationId(observations[1].content_hash())];
    let mixed = card_with_observations(original.clone(), observations.clone());
    assert!(resolve_findley(&mixed, &at, &request, [SINGLE; 2]).is_err());
    original[1].observations = original[0].observations.clone();
    let common = card_with_observations(original.clone(), observations);
    let calibration = resolve_findley(&common, &at, &request, [SINGLE; 2]).unwrap();
    for property in FindleyProperty::ALL {
        assert_eq!(
            calibration.source().unwrap().claim(property).observations,
            original[property.index()].observations
        );
    }
}

#[test]
fn source_and_query_both_need_explicit_static_absolute_temperature() {
    let request = request();
    let material = card(rows(&request));
    let untyped = QueryPoint::new()
        .with(&request.temperature_axis, 300.0)
        .unwrap();
    assert!(resolve_findley(&material, &untyped, &request, [SINGLE; 2]).is_err());
    assert!(resolve_findley(&material, &QueryPoint::new(), &request, [SINGLE; 2]).is_err());
    for quantity in [
        QuantitySpec::dimensional(Temperature::DIMS),
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::TemperatureDifference,
            ValueForm::Static,
        )),
    ] {
        let at = QueryPoint::new()
            .with_quantity(&request.temperature_axis, quantity, 300.0)
            .unwrap();
        assert!(resolve_findley(&material, &at, &request, [SINGLE; 2]).is_err());
        // Equal dimensions do not make either source convention absolute temperature.
        let mut original = rows(&request);
        original[0].validity = ValidityDomain::unconstrained().with_quantity(
            &request.temperature_axis,
            quantity,
            300.0,
            400.0,
        );
        assert!(
            resolve_findley(
                &card(original),
                &point(&request, 300.0),
                &request,
                [SINGLE; 2]
            )
            .is_err()
        );
    }
    for missing_domain in FindleyProperty::ALL {
        let mut original = rows(&request);
        original[missing_domain.index()].validity = ValidityDomain::unconstrained();
        assert!(
            resolve_findley(
                &card(original),
                &point(&request, 300.0),
                &request,
                [SINGLE; 2]
            )
            .is_err()
        );
    }
}

#[test]
fn unsupported_parameter_dependencies_and_nonphysical_fits_refuse() {
    let request = request();
    let mut original = rows(&request);
    let row = &mut original[FindleyProperty::LimitFactor.index()];
    row.value = PropertyValue::Curve {
        abscissa: "fitted-cycles".into(),
        abscissa_dims: Dims::NONE,
        knots: vec![(1.0, 200e6), (100.0, 100e6)],
        dims: Pressure::DIMS,
    };
    row.interpolation = InterpolationPolicy::LinearInside;
    row.validity = row.validity.clone().with_quantity(
        "fitted-cycles",
        QuantitySpec::dimensional(Dims::NONE),
        1.0,
        100.0,
    );
    let at = point(&request, 300.0)
        .with_quantity("fitted-cycles", QuantitySpec::dimensional(Dims::NONE), 50.0)
        .unwrap();
    let material = card(original);
    // This is a valid matdb query, but it does not declare one Findley limit factor.
    assert!(
        material
            .claims()
            .query_typed(
                &FindleyProperty::LimitFactor.property_key(),
                &at,
                SelectionPolicy::SingleClaimOnly
            )
            .is_ok()
    );
    assert!(resolve_findley(&material, &at, &request, [SINGLE; 2]).is_err());
    for (property, invalid) in [
        (FindleyProperty::NormalStressSensitivity, -0.1),
        (FindleyProperty::LimitFactor, 0.0),
        (FindleyProperty::LimitFactor, -1.0),
    ] {
        let mut original = rows(&request);
        let row = &mut original[property.index()];
        row.value = PropertyValue::Scalar {
            value: invalid,
            dims: row.key.dims(),
        };
        assert!(
            resolve_findley(
                &card(original),
                &point(&request, 300.0),
                &request,
                [SINGLE; 2]
            )
            .is_err()
        );
    }
}

#[test]
fn native_binding_requires_the_original_card_state_and_all_source_coordinates() {
    let request = request();
    let mut original = rows(&request);
    for row in &mut original {
        row.validity = row.validity.clone().with_quantity(
            "confining-pressure",
            QuantitySpec::dimensional(Pressure::DIMS),
            1e5,
            2e5,
        );
    }
    let at = point(&request, 300.0)
        .with_quantity(
            "confining-pressure",
            QuantitySpec::dimensional(Pressure::DIMS),
            1e5,
        )
        .unwrap();
    let material = card(original);
    let calibration = resolve_findley(&material, &at, &request, [SINGLE; 2]).unwrap();
    let source = calibration.source().unwrap();
    assert!(source.matches_material_state(material.content_hash(), material.id(), &at));
    let extra_native_axis = at
        .clone()
        .with_quantity(
            "native-only-coordinate",
            QuantitySpec::dimensional(Dims::NONE),
            0.25,
        )
        .unwrap();
    assert!(source.matches_material_state(
        material.content_hash(),
        material.id(),
        &extra_native_axis
    ));
    assert!(!source.matches_material_state(ContentHash([77; 32]), material.id(), &at));
    let mut other_state = material.id().clone();
    other_state.revision += 1;
    assert!(!source.matches_material_state(material.content_hash(), &other_state, &at));
    for changed in [
        at.clone().with(&request.temperature_axis, 350.0).unwrap(),
        at.clone().with("confining-pressure", 2e5).unwrap(),
        point(&request, 300.0),
        at.clone()
            .with_quantity(
                &request.temperature_axis,
                QuantitySpec::dimensional(Temperature::DIMS),
                300.0,
            )
            .unwrap(),
        at.clone()
            .with_quantity(
                "confining-pressure",
                QuantitySpec::dimensional(Dims::NONE),
                1e5,
            )
            .unwrap(),
    ] {
        assert!(!source.matches_material_state(material.content_hash(), material.id(), &changed));
    }
}
