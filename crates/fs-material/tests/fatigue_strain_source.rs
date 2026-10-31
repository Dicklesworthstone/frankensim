//! G0/G1/G3: original CMB calibration, count conventions and native binding.
//! Synthetic coefficients test implementation; they are not measured life data.
#![deny(unsafe_code)]

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimId, ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::fatigue::strain_life::source::{
    CmbPowerBasis, CmbProperty, CmbSourceRequest, resolve_cmb,
};
use fs_qty::{Dims, QuantitySpec, Temperature};

const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

fn request(basis: CmbPowerBasis) -> CmbSourceRequest {
    CmbSourceRequest {
        power_basis: basis,
        temperature_axis: "T".into(),
        stress_ratio_axis: "calibration-stress-ratio".into(),
        surface_and_geometry: "smooth axial coupon".into(),
        environment: "synthetic dry atmosphere".into(),
        load_basis: "total mechanical strain amplitude; fully reversed calibration".into(),
    }
}
fn state(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            temperature,
        )
        .unwrap()
}
fn point(request: &CmbSourceRequest) -> QueryPoint {
    state(300.0)
        .with_quantity(
            &request.stress_ratio_axis,
            QuantitySpec::dimensional(Dims::NONE),
            -1.0,
        )
        .unwrap()
        .with_quantity(
            request.power_basis.life_axis(),
            QuantitySpec::dimensional(Dims::NONE),
            100.0,
        )
        .unwrap()
}
fn rows(request: &CmbSourceRequest) -> Vec<PropertyClaim> {
    let source_factor = if request.power_basis == CmbPowerBasis::FullCycles {
        2.0_f64.powf(-0.5)
    } else {
        1.0
    };
    let lives = if request.power_basis == CmbPowerBasis::FullCycles {
        (0.5, 10_000.0)
    } else {
        (1.0, 20_000.0)
    };
    let values = [
        1000.0,
        10.0 * source_factor,
        -0.5,
        0.04 * source_factor,
        -0.5,
    ];
    CmbProperty::ALL
        .into_iter()
        .map(|property| {
            let key = property.property_key();
            let mut validity = ValidityDomain::unconstrained().with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                300.0,
                400.0,
            );
            if property != CmbProperty::YoungsModulus {
                validity = validity
                    .with_quantity(
                        &request.stress_ratio_axis,
                        QuantitySpec::dimensional(Dims::NONE),
                        -1.0,
                        -1.0,
                    )
                    .with_quantity(
                        request.power_basis.life_axis(),
                        QuantitySpec::dimensional(Dims::NONE),
                        lives.0,
                        lives.1,
                    );
            }
            PropertyClaim {
                value: PropertyValue::Scalar {
                    value: values[property.index()],
                    dims: key.dims(),
                },
                key,
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: if property == CmbProperty::YoungsModulus {
                        "synthetic elastic modulus"
                    } else {
                        "synthetic common CMB fit"
                    }
                    .into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            }
        })
        .collect()
}
fn card(rows: Vec<PropertyClaim>) -> MaterialCard {
    let mut claims = ClaimSet::new();
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
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-12 * expected.abs().max(1e-10),
        "{actual:e} vs {expected:e}"
    );
}

#[test]
fn g1_full_cycle_and_reversal_fits_produce_the_same_physical_curve() {
    let cycle_request = request(CmbPowerBasis::FullCycles);
    let reversal_request = request(CmbPowerBasis::Reversals);
    let cycles = card(rows(&cycle_request));
    let reversals = card(rows(&reversal_request));
    let a = resolve_cmb(&cycles, &point(&cycle_request), &cycle_request, [SINGLE; 5]).unwrap();
    let b = resolve_cmb(
        &reversals,
        &point(&reversal_request),
        &reversal_request,
        [SINGLE; 5],
    )
    .unwrap();
    assert_eq!(a.cycle_domain(), (0.5, 10_000.0));
    assert_eq!(a.cycle_domain(), b.cycle_domain());
    close(a.parameters().fatigue_strength_coefficient_pa, 10.0);
    close(a.parameters().fatigue_ductility_coefficient, 0.04);
    for n in [0.5, 10.0, 50.0, 1000.0, 10_000.0] {
        let expected = 0.05 / (2.0_f64 * n).sqrt();
        close(a.amplitude_at_cycles(n).unwrap(), expected);
        close(b.amplitude_at_cycles(n).unwrap(), expected);
    }
    // Original coefficients retain the source convention; only the native
    // derived curve has reversal-based coefficients.
    close(
        a.source()
            .unwrap()
            .property(CmbProperty::FatigueStrengthCoefficient)
            .evidence
            .value
            .value,
        10.0 / 2.0_f64.sqrt(),
    );
    assert!(a.cycles_to_failure(0.06).is_err());
    assert!(a.cycles_to_failure(1e-5).is_err());
}

#[test]
fn g0_original_fit_claims_uncertainty_and_receipts_survive() {
    let request = request(CmbPowerBasis::Reversals);
    let original = rows(&request);
    let material = card(original.clone());
    let at = point(&request);
    let curve = resolve_cmb(&material, &at, &request, [SINGLE; 5]).unwrap();
    let source = curve.source().unwrap();
    assert_eq!(source.card_identity(), material.content_hash());
    assert_eq!(source.material_state(), material.id());
    assert_eq!(source.query(), &at);
    assert_eq!(source.request(), &request);
    assert_eq!(source.selections(), &[SINGLE; 5]);
    for property in CmbProperty::ALL {
        assert_eq!(source.claim(property), &original[property.index()]);
        let expected = material
            .claims()
            .query_typed(
                &property.property_key(),
                &at,
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
        assert_eq!(source.property(property), &expected);
        material
            .claims()
            .verify_receipt(&source.property(property).receipt)
            .unwrap();
    }
    assert_eq!(curve.scope().material, material.id().chemistry);
    assert_eq!(curve.scope().temperature_k, 300.0);
}

#[test]
fn g0_ambiguous_coefficients_need_explicit_individual_pins() {
    let request = request(CmbPowerBasis::Reversals);
    let mut original = rows(&request);
    let keep = ClaimId(original[1].content_hash());
    let mut alternative = original[1].clone();
    alternative.value = PropertyValue::Scalar {
        value: 15.0,
        dims: alternative.key.dims(),
    };
    original.push(alternative);
    let material = card(original);
    assert!(resolve_cmb(&material, &point(&request), &request, [SINGLE; 5]).is_err());
    let mut selections = [SINGLE; 5];
    selections[CmbProperty::FatigueStrengthCoefficient.index()] = ClaimSelection::Pinned(keep);
    let curve = resolve_cmb(&material, &point(&request), &request, selections).unwrap();
    assert_eq!(curve.parameters().fatigue_strength_coefficient_pa, 10.0);
    assert_eq!(
        curve
            .source()
            .unwrap()
            .property(CmbProperty::FatigueStrengthCoefficient)
            .receipt
            .selected,
        keep
    );
}

#[test]
fn g0_mixed_fits_missing_coefficients_and_mismatched_life_domains_refuse() {
    let request = request(CmbPowerBasis::Reversals);
    let mut different_fit = rows(&request);
    different_fit[4].provenance.source = "unrelated population fit".into();
    assert!(
        resolve_cmb(
            &card(different_fit),
            &point(&request),
            &request,
            [SINGLE; 5]
        )
        .is_err()
    );
    let mut incomplete = rows(&request);
    incomplete.remove(3);
    assert!(resolve_cmb(&card(incomplete), &point(&request), &request, [SINGLE; 5]).is_err());
    let mut different_domain = rows(&request);
    different_domain[3].validity = different_domain[3].validity.clone().with_quantity(
        request.power_basis.life_axis(),
        QuantitySpec::dimensional(Dims::NONE),
        1.0,
        10_000.0,
    );
    assert!(
        resolve_cmb(
            &card(different_domain),
            &point(&request),
            &request,
            [SINGLE; 5]
        )
        .is_err()
    );
}

#[test]
fn g0_baseline_temperature_and_count_quantity_are_actual_source_requirements() {
    let request = request(CmbPowerBasis::Reversals);
    let material = card(rows(&request));
    let wrong_ratio = point(&request)
        .with(&request.stress_ratio_axis, 0.0)
        .unwrap();
    assert!(resolve_cmb(&material, &wrong_ratio, &request, [SINGLE; 5]).is_err());
    let hot = point(&request).with("T", 401.0).unwrap();
    assert!(resolve_cmb(&material, &hot, &request, [SINGLE; 5]).is_err());
    let untyped = QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with(&request.stress_ratio_axis, -1.0)
        .unwrap()
        .with(request.power_basis.life_axis(), 100.0)
        .unwrap();
    assert!(resolve_cmb(&material, &untyped, &request, [SINGLE; 5]).is_err());
    let mut undeclared_baseline = rows(&request);
    undeclared_baseline[1].validity = ValidityDomain::unconstrained()
        .with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            300.0,
            400.0,
        )
        .with_quantity(
            request.power_basis.life_axis(),
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
            20_000.0,
        );
    assert!(
        resolve_cmb(
            &card(undeclared_baseline),
            &point(&request),
            &request,
            [SINGLE; 5]
        )
        .is_err()
    );
}

#[test]
fn g1_temperature_curve_preserves_source_interpolation_and_life_fit_bounds() {
    let request = request(CmbPowerBasis::Reversals);
    let mut original = rows(&request);
    original[1].value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, 10.0), (400.0, 8.0)],
        dims: original[1].key.dims(),
    };
    original[1].interpolation = InterpolationPolicy::LinearInside;
    let material = card(original);
    let at = point(&request).with("T", 350.0).unwrap();
    let curve = resolve_cmb(&material, &at, &request, [SINGLE; 5]).unwrap();
    assert_eq!(curve.parameters().fatigue_strength_coefficient_pa, 9.0);
    assert_eq!(curve.cycle_domain(), (0.5, 10_000.0));
    assert_eq!(curve.scope().temperature_k, 350.0);
    material
        .claims()
        .verify_receipt(
            &curve
                .source()
                .unwrap()
                .property(CmbProperty::FatigueStrengthCoefficient)
                .receipt,
        )
        .unwrap();
}

#[test]
fn g0_life_dependent_parameters_and_nonphysical_domains_refuse() {
    let request = request(CmbPowerBasis::Reversals);
    let mut wrong_dependency = rows(&request);
    wrong_dependency[1].value = PropertyValue::Curve {
        abscissa: request.power_basis.life_axis().into(),
        abscissa_dims: Dims::NONE,
        knots: vec![(1.0, 10.0), (20_000.0, 8.0)],
        dims: wrong_dependency[1].key.dims(),
    };
    wrong_dependency[1].interpolation = InterpolationPolicy::LinearInside;
    assert!(
        resolve_cmb(
            &card(wrong_dependency),
            &point(&request),
            &request,
            [SINGLE; 5]
        )
        .is_err()
    );
    let mut wrong_life = rows(&request);
    for row in wrong_life.iter_mut().skip(1) {
        row.validity = row.validity.clone().with_quantity(
            request.power_basis.life_axis(),
            QuantitySpec::dimensional(Dims::NONE),
            0.1,
            20_000.0,
        );
    }
    assert!(resolve_cmb(&card(wrong_life), &point(&request), &request, [SINGLE; 5]).is_err());
}

#[test]
fn g3_native_binding_checks_original_card_context_and_modulus_source() {
    let request = request(CmbPowerBasis::Reversals);
    let material = card(rows(&request));
    let curve = resolve_cmb(&material, &point(&request), &request, [SINGLE; 5]).unwrap();
    let source = curve.source().unwrap();
    let at = state(300.0);
    let youngs = material
        .claims()
        .query_typed(
            &CmbProperty::YoungsModulus.property_key(),
            &at,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    assert!(source.matches_material_state(material.content_hash(), material.id(), &at, &youngs));
    assert!(!source.matches_material_state(
        material.content_hash(),
        material.id(),
        &state(301.0),
        &youngs
    ));
    let untyped = QueryPoint::new().with("T", 300.0).unwrap();
    assert!(!source.matches_material_state(
        material.content_hash(),
        material.id(),
        &untyped,
        &youngs
    ));
    let mut changed = rows(&request);
    changed[0].provenance.source = "different elastic source with coincident value".into();
    let other = card(changed);
    let other_youngs = other
        .claims()
        .query_typed(
            &CmbProperty::YoungsModulus.property_key(),
            &at,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    assert!(!source.matches_material_state(other.content_hash(), other.id(), &at, &other_youngs));
    assert!(!source.matches_material_state(
        material.content_hash(),
        material.id(),
        &at,
        &other_youngs
    ));
}
