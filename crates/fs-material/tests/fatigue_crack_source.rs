//! G0/G1/G3: original Paris claims, squared-intensity units and native binding.
//! Synthetic fits verify implementation, not measured material or component life.
#![deny(unsafe_code)]

use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimId, ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    ObservationDataset, ObservationId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::fatigue::crack_growth::source::{
    ParisProperty, ParisSourceRequest, STRESS_INTENSITY_SQUARED_DIMS, resolve_paris,
};
use fs_material::fatigue::crack_growth::{
    ConstantGeometryCrack, CrackGrowthError, CrackGrowthLimits, ParisCrackGrowth, StressCycle,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Length, Pressure, QuantitySpec, Temperature};

const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
const VALUES: [f64; 4] = [1e12, 1e-9, 2.0, 1e14];

fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn request() -> ParisSourceRequest {
    ParisSourceRequest {
        temperature_axis: "T".into(),
        stress_ratio_axis: "calibration-R".into(),
        delta_k_squared_axis: "calibration-delta-K-squared".into(),
        surface_and_geometry: "synthetic through-cracked coupon".into(),
        environment: "synthetic dry atmosphere".into(),
        load_basis: "tensile Mode I; full Delta-K range; growth per full cycle".into(),
    }
}
fn physical_point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", absolute_temperature(), temperature)
        .unwrap()
        .with_quantity(
            "ambient-pressure",
            QuantitySpec::dimensional(Pressure::DIMS),
            1e5,
        )
        .unwrap()
}
fn point(request: &ParisSourceRequest, temperature: f64) -> QueryPoint {
    physical_point(temperature)
        .with_quantity(
            &request.stress_ratio_axis,
            QuantitySpec::dimensional(Dims::NONE),
            0.25,
        )
        .unwrap()
        .with_quantity(
            &request.delta_k_squared_axis,
            QuantitySpec::dimensional(STRESS_INTENSITY_SQUARED_DIMS),
            1e12,
        )
        .unwrap()
}
fn rows(request: &ParisSourceRequest) -> Vec<PropertyClaim> {
    let validity = ValidityDomain::unconstrained()
        .with_quantity("T", absolute_temperature(), 300.0, 400.0)
        .with_quantity(
            &request.stress_ratio_axis,
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
            0.8,
        )
        .with_quantity(
            &request.delta_k_squared_axis,
            QuantitySpec::dimensional(STRESS_INTENSITY_SQUARED_DIMS),
            1e8,
            9e12,
        )
        .with_quantity(
            "ambient-pressure",
            QuantitySpec::dimensional(Pressure::DIMS),
            1e5,
            2e5,
        );
    ParisProperty::ALL
        .into_iter()
        .map(|property| {
            let key = property.property_key();
            PropertyClaim {
                value: PropertyValue::Scalar {
                    value: VALUES[property.index()],
                    dims: key.dims(),
                },
                key,
                validity: validity.clone(),
                uncertainty: UncertaintyModel::RelativeHalfWidth {
                    fraction: 0.05,
                    confidence: 0.9,
                },
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic common normalized Paris fit".into(),
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
        (actual - expected).abs() <= 2e-12 * expected.abs(),
        "{actual:e} vs {expected:e}"
    );
}

#[test]
fn g1_original_fit_drives_growth_and_keeps_all_source_evidence() {
    let request = request();
    let original = rows(&request);
    let material = card(original.clone());
    let at = point(&request, 300.0);
    let law = resolve_paris(&material, &at, &request, [SINGLE; 4]).unwrap();
    assert_eq!(law.parameters().reference_delta_k_pa_sqrt_m, 1e6);
    assert_eq!(law.parameters().reference_rate_m_per_cycle, 1e-9);
    assert_eq!(law.parameters().exponent, 2.0);
    assert_eq!(law.domain().minimum_delta_k_pa_sqrt_m, 1e4);
    assert_eq!(law.domain().maximum_delta_k_pa_sqrt_m, 3e6);
    assert_eq!(law.domain().maximum_k_pa_sqrt_m, 1e7);
    assert_eq!(
        (
            law.domain().minimum_load_ratio,
            law.domain().maximum_load_ratio
        ),
        (0.0, 0.8)
    );
    close(law.growth_rate(2e6, 0.25).unwrap(), 4e-9);
    let source = law.source().unwrap();
    assert_eq!(source.card_identity(), material.content_hash());
    assert_eq!(source.material_state(), material.id());
    assert_eq!(source.request(), &request);
    assert_eq!(source.query(), &at);
    assert_eq!(source.selections(), &[SINGLE; 4]);
    for property in ParisProperty::ALL {
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
    assert_eq!(law.scope().source_reference, original[0].provenance.source);
    assert_eq!(law.scope().material, material.id().chemistry);
    assert_eq!(law.scope().process, material.id().process);
    assert_eq!(law.scope().temperature_k, 300.0);
    let mut growth = ParisCrackGrowth::new(
        Arc::new(law),
        ConstantGeometryCrack {
            geometry_factor: 1.0,
            minimum_length_m: 1e-5,
            maximum_length_m: 0.01,
        },
        1e-3,
        CrackGrowthLimits {
            maximum_blocks: 4,
            maximum_cycles: 1e6,
        },
    )
    .unwrap();
    let step = growth
        .advance(StressCycle {
            minimum_pa: 2e6,
            maximum_pa: 8e6,
            cycles: 1000.0,
        })
        .unwrap();
    let expected = 1e-3 * (1e-9 * 36.0 * std::f64::consts::PI * 1000.0).exp();
    close(step.after.crack_length_m, expected);
    assert_eq!(
        growth.law().source().unwrap().card_identity(),
        material.content_hash()
    );
}

#[test]
fn g1_temperature_interpolation_keeps_the_original_knots_and_fixed_context() {
    let request = request();
    let mut original = rows(&request);
    for property in ParisProperty::ALL {
        let row = &mut original[property.index()];
        let start = VALUES[property.index()];
        let end = if property == ParisProperty::ReferenceGrowthRate {
            3e-9
        } else {
            start
        };
        row.value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, start), (400.0, end)],
            dims: row.key.dims(),
        };
        row.interpolation = InterpolationPolicy::LinearInside;
    }
    let material = card(original.clone());
    let law = resolve_paris(&material, &point(&request, 350.0), &request, [SINGLE; 4]).unwrap();
    close(law.parameters().reference_rate_m_per_cycle, 2e-9);
    assert_eq!(law.scope().temperature_k, 350.0);
    for property in ParisProperty::ALL {
        assert_eq!(
            law.source().unwrap().claim(property),
            &original[property.index()]
        );
    }
    for outside in [299.0, 401.0] {
        assert!(
            resolve_paris(&material, &point(&request, outside), &request, [SINGLE; 4]).is_err()
        );
    }
    assert!(!law.source().unwrap().matches_material_state(
        material.content_hash(),
        material.id(),
        &physical_point(300.0),
    ));
}

#[test]
fn g0_missing_growth_data_and_ambiguous_fits_are_not_inferred() {
    let request = request();
    let at = point(&request, 300.0);
    let mut incomplete = rows(&request);
    let rate = &mut incomplete[ParisProperty::ReferenceGrowthRate.index()];
    rate.key = PropertyKey::new("paris-reference-growth-per-reversal", Length::DIMS);
    assert!(resolve_paris(&card(incomplete), &at, &request, [SINGLE; 4]).is_err());
    let mut wrong_strength = rows(&request);
    let exponent = &mut wrong_strength[ParisProperty::Exponent.index()];
    exponent.key = PropertyKey::new("ultimate-tensile-strength", Pressure::DIMS);
    exponent.value = PropertyValue::Scalar {
        value: 500e6,
        dims: Pressure::DIMS,
    };
    assert!(resolve_paris(&card(wrong_strength), &at, &request, [SINGLE; 4]).is_err());
    let mut original = rows(&request);
    let index = ParisProperty::Exponent.index();
    let keep = ClaimId(original[index].content_hash());
    let mut alternative = original[index].clone();
    alternative.value = PropertyValue::Scalar {
        value: 3.0,
        dims: Dims::NONE,
    };
    original.push(alternative);
    let material = card(original);
    assert!(resolve_paris(&material, &at, &request, [SINGLE; 4]).is_err());
    let mut selections = [SINGLE; 4];
    selections[index] = ClaimSelection::Pinned(keep);
    let law = resolve_paris(&material, &at, &request, selections).unwrap();
    assert_eq!(law.parameters().exponent, 2.0);
    assert_eq!(law.source().unwrap().selections(), &selections);
}

#[test]
fn g0_fit_population_and_complete_support_cannot_be_mixed() {
    let request = request();
    let at = point(&request, 300.0);
    let mut different_fit = rows(&request);
    different_fit[2].provenance.source = "different synthetic fit".into();
    assert!(resolve_paris(&card(different_fit), &at, &request, [SINGLE; 4]).is_err());
    for (axis, low, high) in [
        ("T", 300.0, 450.0),
        (request.stress_ratio_axis.as_str(), 0.0, 0.9),
        (request.delta_k_squared_axis.as_str(), 1e8, 16e12),
        ("ambient-pressure", 1e5, 3e5),
    ] {
        let mut different_domain = rows(&request);
        different_domain[2].validity = different_domain[2].validity.clone().with(axis, low, high);
        assert!(resolve_paris(&card(different_domain), &at, &request, [SINGLE; 4]).is_err());
    }
    let mut original = rows(&request);
    let observations: Vec<_> = [11, 12]
        .into_iter()
        .map(|tag| ObservationDataset {
            specimen: format!("synthetic population {tag}"),
            method: "synthetic full-cycle Paris fit".into(),
            artifact: ContentHash([tag; 32]),
            caveats: "test-only".into(),
            provenance: original[0].provenance.clone(),
        })
        .collect();
    for row in &mut original {
        row.observations = vec![ObservationId(observations[0].content_hash())];
    }
    original[2].observations = vec![ObservationId(observations[1].content_hash())];
    assert!(
        resolve_paris(
            &card_with_observations(original.clone(), observations.clone()),
            &at,
            &request,
            [SINGLE; 4],
        )
        .is_err()
    );
    original[2].observations = original[0].observations.clone();
    let material = card_with_observations(original.clone(), observations);
    let law = resolve_paris(&material, &at, &request, [SINGLE; 4]).unwrap();
    assert_eq!(
        law.source()
            .unwrap()
            .claim(ParisProperty::Exponent)
            .observations,
        original[0].observations
    );
}

#[test]
fn g0_squared_intensity_and_absolute_temperature_are_explicit_in_source_and_query() {
    let request = request();
    let material = card(rows(&request));
    let at = point(&request, 300.0);
    assert_eq!(
        ParisProperty::ReferenceDeltaKSquared.property_key().dims(),
        STRESS_INTENSITY_SQUARED_DIMS
    );
    assert_eq!(
        ParisProperty::ReferenceGrowthRate.property_key().dims(),
        Length::DIMS
    );
    assert_ne!(STRESS_INTENSITY_SQUARED_DIMS, Pressure::DIMS);
    for (axis, quantity, value) in [
        ("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0),
        (
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::TemperatureDifference,
                ValueForm::Static,
            )),
            300.0,
        ),
        (
            request.stress_ratio_axis.as_str(),
            QuantitySpec::dimensional(Pressure::DIMS),
            0.25,
        ),
        (
            request.delta_k_squared_axis.as_str(),
            QuantitySpec::dimensional(Pressure::DIMS),
            1e12,
        ),
    ] {
        let wrong = at.clone().with_quantity(axis, quantity, value).unwrap();
        assert!(resolve_paris(&material, &wrong, &request, [SINGLE; 4]).is_err());
    }
    let mut missing = rows(&request);
    for row in &mut missing {
        row.validity = ValidityDomain::unconstrained()
            .with_quantity("T", absolute_temperature(), 300.0, 400.0)
            .with_quantity(
                &request.stress_ratio_axis,
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
                0.8,
            );
    }
    assert!(resolve_paris(&card(missing), &at, &request, [SINGLE; 4]).is_err());
    let mut wrong_source = rows(&request);
    for row in &mut wrong_source {
        row.validity = row.validity.clone().with_quantity(
            &request.delta_k_squared_axis,
            QuantitySpec::dimensional(Pressure::DIMS),
            1e8,
            9e12,
        );
    }
    let wrong_at = at
        .clone()
        .with_quantity(
            &request.delta_k_squared_axis,
            QuantitySpec::dimensional(Pressure::DIMS),
            1e12,
        )
        .unwrap();
    assert!(resolve_paris(&card(wrong_source), &wrong_at, &request, [SINGLE; 4]).is_err());
    let mut aliased = request.clone();
    aliased.delta_k_squared_axis = aliased.stress_ratio_axis.clone();
    assert!(resolve_paris(&material, &at, &aliased, [SINGLE; 4]).is_err());
}

#[test]
fn g0_crossdependent_fits_and_surface_bounding_boxes_do_not_become_a_paris_domain() {
    let request = request();
    let at = point(&request, 300.0);
    for surface in [false, true] {
        let mut original = rows(&request);
        let row = &mut original[ParisProperty::Exponent.index()];
        row.value = if surface {
            PropertyValue::Surface2D {
                axes: [
                    request.stress_ratio_axis.clone(),
                    request.delta_k_squared_axis.clone(),
                ],
                axis_dims: [Dims::NONE, STRESS_INTENSITY_SQUARED_DIMS],
                samples: vec![([0.0, 1e8], 2.0), ([0.8, 1e8], 3.0), ([0.0, 9e12], 4.0)],
                triangles: vec![[0, 1, 2]],
                dims: Dims::NONE,
            }
        } else {
            PropertyValue::Curve {
                abscissa: request.stress_ratio_axis.clone(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0.0, 2.0), (0.8, 3.0)],
                dims: Dims::NONE,
            }
        };
        row.interpolation = InterpolationPolicy::LinearInside;
        let material = card(original);
        assert!(
            material
                .claims()
                .query_typed(
                    &ParisProperty::Exponent.property_key(),
                    &at,
                    SelectionPolicy::SingleClaimOnly,
                )
                .is_ok()
        );
        assert!(resolve_paris(&material, &at, &request, [SINGLE; 4]).is_err());
    }
}

#[test]
fn g0_source_bounds_and_exclusive_k_cutoff_reach_the_executable_law() {
    let request = request();
    let at = point(&request, 300.0);
    let law = resolve_paris(&card(rows(&request)), &at, &request, [SINGLE; 4]).unwrap();
    assert!(matches!(
        law.growth_rate(1e3, 0.25),
        Err(CrackGrowthError::DeltaKOutsideDomain { .. })
    ));
    assert!(matches!(
        law.growth_rate(4e6, 0.25),
        Err(CrackGrowthError::DeltaKOutsideDomain { .. })
    ));
    assert!(matches!(
        law.growth_rate(1e6, 0.9),
        Err(CrackGrowthError::LoadRatioOutsideDomain { .. })
    ));
    assert!(matches!(
        law.growth_rate(2e6, 0.8),
        Err(CrackGrowthError::MaximumKExceeded { .. })
    ));
    for (property, bad) in [
        (ParisProperty::ReferenceDeltaKSquared, -1.0),
        (ParisProperty::ReferenceGrowthRate, 0.0),
        (ParisProperty::Exponent, -2.0),
        (ParisProperty::MaximumKSquared, 1e10),
    ] {
        let mut original = rows(&request);
        original[property.index()].value = PropertyValue::Scalar {
            value: bad,
            dims: property.property_key().dims(),
        };
        assert!(resolve_paris(&card(original), &at, &request, [SINGLE; 4]).is_err());
    }
}

#[test]
fn g3_native_binding_retains_every_physical_axis_but_excludes_fit_coordinates() {
    let request = request();
    let material = card(rows(&request));
    let at = point(&request, 300.0);
    let law = resolve_paris(&material, &at, &request, [SINGLE; 4]).unwrap();
    let source = law.source().unwrap();
    let native = physical_point(300.0);
    assert!(source.matches_material_state(material.content_hash(), material.id(), &native));
    let extra = native
        .clone()
        .with_quantity("native-only", QuantitySpec::dimensional(Dims::NONE), 1.0)
        .unwrap();
    assert!(source.matches_material_state(material.content_hash(), material.id(), &extra));
    assert!(!source.matches_material_state(ContentHash([99; 32]), material.id(), &native));
    let mut another_state = material.id().clone();
    another_state.revision += 1;
    assert!(!source.matches_material_state(material.content_hash(), &another_state, &native));
    for mismatch in [
        physical_point(350.0),
        native.clone().with("ambient-pressure", 2e5).unwrap(),
        QueryPoint::new()
            .with_quantity("T", absolute_temperature(), 300.0)
            .unwrap(),
        native
            .clone()
            .with_quantity(
                "ambient-pressure",
                QuantitySpec::dimensional(Dims::NONE),
                1e5,
            )
            .unwrap(),
    ] {
        assert!(!source.matches_material_state(material.content_hash(), material.id(), &mismatch));
    }
}
