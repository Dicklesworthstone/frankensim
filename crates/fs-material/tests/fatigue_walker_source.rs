//! Original Walker calibration drives the shared crack-growth state.
//! Synthetic fits verify implementation, not material or component lifetime.
#![deny(unsafe_code)]

use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimId, ClaimSelection, ClaimSet, InterpolationPolicy, MaterialAnswer, MaterialCard,
    MaterialStateId, ObservationDataset, ObservationId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::fatigue::crack_growth::source::{
    ParisCalibrationSource, ParisProperty, ParisSourceRequest, STRESS_INTENSITY_SQUARED_DIMS,
    WalkerProperty, resolve_paris, resolve_walker,
};
use fs_material::fatigue::crack_growth::{
    ConstantGeometryCrack, CrackGrowthLimits, ParisCrackGrowth, StressCycle, WalkerParameters,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
const VALUES: [f64; 6] = [1e12, 1e-9, 2.0, 1e14, 0.25, 0.5];

fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn request() -> ParisSourceRequest {
    ParisSourceRequest {
        temperature_axis: "absolute-temperature".into(),
        stress_ratio_axis: "calibration-R".into(),
        delta_k_squared_axis: "calibration-delta-K-squared".into(),
        surface_and_geometry: "synthetic through-cracked coupon".into(),
        environment: "synthetic dry atmosphere".into(),
        load_basis: "tensile Mode I; full Delta-K range; growth per full cycle".into(),
    }
}

fn physical_point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("absolute-temperature", absolute_temperature(), temperature)
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
            0.1, // Query selection R is deliberately not the fit's reference R.
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
        .with_quantity(
            &request.temperature_axis,
            absolute_temperature(),
            300.0,
            400.0,
        )
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
    WalkerProperty::ALL
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
                    source: "synthetic common normalized Walker fit".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            }
        })
        .collect()
}

fn material(rows: Vec<PropertyClaim>, observations: Vec<ObservationDataset>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for observation in observations {
        claims.register_observation(observation).unwrap();
    }
    for row in rows {
        claims.insert_claim(row).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic Walker alloy".into(),
            phase: "solid".into(),
            process: "synthetic coupon process".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn original(
    source: &ParisCalibrationSource,
    property: WalkerProperty,
) -> (&PropertyClaim, &MaterialAnswer) {
    let walker = source.walker().unwrap();
    match property {
        WalkerProperty::ReferenceLoadRatio => (
            walker.reference_load_ratio_claim(),
            walker.reference_load_ratio_answer(),
        ),
        WalkerProperty::Gamma => (walker.gamma_claim(), walker.gamma_answer()),
        common => {
            let role = ParisProperty::ALL[common.index()];
            (source.claim(role), source.property(role))
        }
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs(),
        "{actual:e} vs {expected:e}"
    );
}

#[test]
fn original_walker_fit_drives_growth_retains_receipts_and_isolates_paris_keys() {
    let request = request();
    let claims = rows(&request);
    let card = material(claims.clone(), vec![]);
    let at = point(&request, 300.0);
    let law = resolve_walker(&card, &at, &request, [SINGLE; 6]).unwrap();
    assert_eq!(
        law.walker(),
        Some(WalkerParameters {
            reference_load_ratio: 0.25,
            gamma: 0.5,
        })
    );
    close(law.growth_rate(1e6, 0.25).unwrap(), 1e-9);
    close(law.growth_rate(2e6, 0.5).unwrap(), 6e-9);
    let source = law.source().unwrap();
    assert_eq!(source.selections(), &[SINGLE; 4]);
    assert_eq!(source.walker().unwrap().selections(), &[SINGLE; 2]);
    assert_eq!(source.query(), &at);
    for property in WalkerProperty::ALL {
        let (claim, answer) = original(source, property);
        assert_eq!(claim, &claims[property.index()]);
        let direct = card
            .claims()
            .query_typed(
                &property.property_key(),
                &at,
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
        assert_eq!(answer, &direct);
        card.claims().verify_receipt(&answer.receipt).unwrap();
    }
    let native = physical_point(300.0);
    assert!(source.matches_material_state(card.content_hash(), card.id(), &native));
    assert!(!source.matches_material_state(ContentHash([99; 32]), card.id(), &native));
    assert!(!source.matches_material_state(
        card.content_hash(),
        card.id(),
        &native.with("ambient-pressure", 2e5).unwrap(),
    ));
    assert!(resolve_paris(&card, &at, &request, [SINGLE; 4]).is_err());
    let mut paris_rows = claims[..4].to_vec();
    for property in ParisProperty::ALL {
        assert_ne!(property.property_key(), paris_rows[property.index()].key);
        paris_rows[property.index()].key = property.property_key();
    }
    let paris_card = material(paris_rows, vec![]);
    assert!(resolve_walker(&paris_card, &at, &request, [SINGLE; 6]).is_err());
    let paris = resolve_paris(&paris_card, &at, &request, [SINGLE; 4]).unwrap();
    assert!(paris.walker().is_none());
    assert!(paris.source().unwrap().walker().is_none());

    let mut growth = ParisCrackGrowth::new(
        Arc::new(law),
        ConstantGeometryCrack {
            geometry_factor: 1.0,
            minimum_length_m: 1e-5,
            maximum_length_m: 0.01,
        },
        1e-3,
        CrackGrowthLimits {
            maximum_blocks: 2,
            maximum_cycles: 1e6,
        },
    )
    .unwrap();
    let step = growth
        .advance(StressCycle {
            minimum_pa: 6e6,
            maximum_pa: 12e6,
            cycles: 1000.0,
        })
        .unwrap();
    // n=2, gamma=1/2, Rref=1/4 and R=1/2: independent multiplier = 3/2.
    let expected = 1e-3 * (1e-9 * 36.0 * std::f64::consts::PI * 1.5 * 1000.0).exp();
    close(step.after.crack_length_m, expected);
    assert_eq!(
        growth.law().source().unwrap().card_identity(),
        card.content_hash()
    );
}

#[test]
fn temperature_interpolation_retains_both_original_walker_coordinates() {
    let request = request();
    let at = point(&request, 350.0);
    let mut claims = rows(&request);
    for (property, start, end) in [
        (WalkerProperty::ReferenceGrowthRate, 1e-9, 3e-9),
        (WalkerProperty::ReferenceLoadRatio, 0.1, 0.3),
        (WalkerProperty::Gamma, 0.2, 0.8),
    ] {
        let row = &mut claims[property.index()];
        row.value = PropertyValue::Curve {
            abscissa: request.temperature_axis.clone(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, start), (400.0, end)],
            dims: row.key.dims(),
        };
        row.interpolation = InterpolationPolicy::LinearInside;
    }
    let card = material(claims.clone(), vec![]);
    let law = resolve_walker(&card, &at, &request, [SINGLE; 6]).unwrap();
    close(law.parameters().reference_rate_m_per_cycle, 2e-9);
    close(law.walker().unwrap().reference_load_ratio, 0.2);
    close(law.walker().unwrap().gamma, 0.5);
    close(law.growth_rate(1e6, 0.5).unwrap(), 3.2e-9);
    assert_eq!(law.source().unwrap().query(), &at);
    for property in [WalkerProperty::ReferenceLoadRatio, WalkerProperty::Gamma] {
        let (claim, answer) = original(law.source().unwrap(), property);
        assert_eq!(claim, &claims[property.index()]);
        card.claims().verify_receipt(&answer.receipt).unwrap();
    }
    assert!(resolve_walker(&card, &point(&request, 450.0), &request, [SINGLE; 6]).is_err());
}

#[test]
fn extra_parameters_require_original_dimensionless_claims_and_explicit_selection() {
    let request = request();
    let at = point(&request, 300.0);
    for property in [WalkerProperty::ReferenceLoadRatio, WalkerProperty::Gamma] {
        let mut missing = rows(&request);
        missing.remove(property.index());
        assert!(resolve_walker(&material(missing, vec![]), &at, &request, [SINGLE; 6]).is_err());
        let mut wrong_units = rows(&request);
        let row = &mut wrong_units[property.index()];
        row.key = PropertyKey::new(row.key.name(), Pressure::DIMS);
        row.value = PropertyValue::Scalar {
            value: VALUES[property.index()],
            dims: Pressure::DIMS,
        };
        assert!(
            resolve_walker(&material(wrong_units, vec![]), &at, &request, [SINGLE; 6]).is_err()
        );
    }
    let mut claims = rows(&request);
    let gamma_index = WalkerProperty::Gamma.index();
    let keep = ClaimId(claims[gamma_index].content_hash());
    let mut alternative = claims[gamma_index].clone();
    alternative.value = PropertyValue::Scalar {
        value: 0.75,
        dims: Dims::NONE,
    };
    claims.push(alternative);
    let card = material(claims, vec![]);
    assert!(resolve_walker(&card, &at, &request, [SINGLE; 6]).is_err());
    let mut selections = [SINGLE; 6];
    selections[gamma_index] = ClaimSelection::Pinned(keep);
    let law = resolve_walker(&card, &at, &request, selections).unwrap();
    close(law.walker().unwrap().gamma, 0.5);
    assert_eq!(
        law.source().unwrap().walker().unwrap().selections(),
        &[SINGLE, selections[gamma_index]]
    );
    let wrong_temperature = at
        .with_quantity(
            &request.temperature_axis,
            QuantitySpec::dimensional(Temperature::DIMS),
            300.0,
        )
        .unwrap();
    assert!(resolve_walker(&card, &wrong_temperature, &request, selections).is_err());
}

#[test]
fn all_six_parameters_require_one_original_fit_population_and_complete_domain() {
    let request = request();
    let at = point(&request, 300.0);
    for property in [WalkerProperty::ReferenceLoadRatio, WalkerProperty::Gamma] {
        let mut other_fit = rows(&request);
        other_fit[property.index()].provenance.source = "different Walker fit".into();
        assert!(resolve_walker(&material(other_fit, vec![]), &at, &request, [SINGLE; 6]).is_err());
        let mut other_domain = rows(&request);
        let row = &mut other_domain[property.index()];
        row.validity = row.validity.clone().with("ambient-pressure", 1e5, 3e5);
        assert!(
            resolve_walker(&material(other_domain, vec![]), &at, &request, [SINGLE; 6]).is_err()
        );
    }
    let mut claims = rows(&request);
    let observations: Vec<_> = [41, 42]
        .into_iter()
        .map(|tag| ObservationDataset {
            specimen: format!("synthetic Walker population {tag}"),
            method: "synthetic full-cycle Walker fit".into(),
            artifact: ContentHash([tag; 32]),
            caveats: "test-only".into(),
            provenance: claims[0].provenance.clone(),
        })
        .collect();
    for row in &mut claims {
        row.observations = vec![ObservationId(observations[0].content_hash())];
    }
    claims[WalkerProperty::Gamma.index()].observations =
        vec![ObservationId(observations[1].content_hash())];
    assert!(
        resolve_walker(
            &material(claims.clone(), observations.clone()),
            &at,
            &request,
            [SINGLE; 6]
        )
        .is_err()
    );
    claims[WalkerProperty::Gamma.index()].observations = claims[0].observations.clone();
    let card = material(claims.clone(), observations);
    let law = resolve_walker(&card, &at, &request, [SINGLE; 6]).unwrap();
    assert_eq!(
        law.source()
            .unwrap()
            .walker()
            .unwrap()
            .gamma_claim()
            .observations,
        claims[0].observations
    );
}

#[test]
fn walker_coefficients_cannot_hide_crossdependent_fits_or_invalid_reference_support() {
    let request = request();
    let at = point(&request, 300.0);
    for surface in [false, true] {
        let mut claims = rows(&request);
        let row = &mut claims[WalkerProperty::Gamma.index()];
        row.value = if surface {
            PropertyValue::Surface2D {
                axes: [
                    request.stress_ratio_axis.clone(),
                    request.delta_k_squared_axis.clone(),
                ],
                axis_dims: [Dims::NONE, STRESS_INTENSITY_SQUARED_DIMS],
                samples: vec![([0.0, 1e8], 0.25), ([0.8, 1e8], 0.5), ([0.0, 9e12], 0.75)],
                triangles: vec![[0, 1, 2]],
                dims: Dims::NONE,
            }
        } else {
            PropertyValue::Curve {
                abscissa: request.stress_ratio_axis.clone(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0.0, 0.25), (0.8, 0.75)],
                dims: Dims::NONE,
            }
        };
        row.interpolation = InterpolationPolicy::LinearInside;
        let card = material(claims, vec![]);
        assert!(
            card.claims()
                .query_typed(
                    &WalkerProperty::Gamma.property_key(),
                    &at,
                    SelectionPolicy::SingleClaimOnly
                )
                .is_ok()
        );
        assert!(resolve_walker(&card, &at, &request, [SINGLE; 6]).is_err());
    }
    for (property, invalid) in [
        (WalkerProperty::Gamma, -0.1),
        (WalkerProperty::Gamma, 1.1),
        (WalkerProperty::ReferenceLoadRatio, -0.1),
        (WalkerProperty::ReferenceLoadRatio, 0.9),
        (WalkerProperty::MaximumKSquared, 1.44e12),
    ] {
        let mut claims = rows(&request);
        claims[property.index()].value = PropertyValue::Scalar {
            value: invalid,
            dims: property.property_key().dims(),
        };
        assert!(resolve_walker(&material(claims, vec![]), &at, &request, [SINGLE; 6]).is_err());
    }
}
