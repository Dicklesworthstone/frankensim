//! Source-preserving affine thermal hardening, with no inferred caloric kinks.
//! Synthetic coefficients test implementation, not an alloy calibration.
#![deny(unsafe_code)]

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimId, ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
    UncertaintyModel,
};
use fs_material::plastic::tabulated::source::{
    ACCUMULATED_PLASTIC_STRAIN_AXIS as ALPHA, CAUCHY_YIELD_CURVE_PROPERTY,
    KIRCHHOFF_YIELD_CURVE_PROPERTY,
};
use fs_material::plastic::tabulated::thermal_source::{
    ThermalHardeningCurveRequest, ThermalHardeningSourceError, resolve_thermal_tabulated_hardening,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

const TEMPERATURE: &str = "specimen-temperature";
const AMBIENT_PRESSURE: &str = "ambient-pressure";
const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
const STRAINS: [f64; 4] = [0.0, 0.01, 0.02, 0.04];
const TEMPERATURES: [f64; 3] = [280.0, 300.0, 320.0];
const REFERENCE_YIELD: [f64; 4] = [3.0, 4.0, 4.0, 8.0];
const SLOPES: [f64; 4] = [-1.0 / 64.0, -1.0 / 32.0, -1.0 / 32.0, -1.0 / 16.0];

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn context(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            TEMPERATURE,
            semantic(QuantityKind::AbsoluteTemperature),
            temperature,
        )
        .unwrap()
        .with_quantity(
            AMBIENT_PRESSURE,
            semantic(QuantityKind::Pressure),
            100_000.0,
        )
        .unwrap()
}
fn request(reference_temperature_k: f64) -> ThermalHardeningCurveRequest {
    ThermalHardeningCurveRequest {
        temperature_axis: TEMPERATURE.into(),
        reference_temperature_k,
        temperature_range_k: [280.0, 320.0],
        maximum_samples: STRAINS.len() * TEMPERATURES.len(),
        maximum_knots: STRAINS.len(),
    }
}
fn curve() -> PropertyClaim {
    let mut samples = Vec::new();
    for temperature in TEMPERATURES {
        for (index, alpha) in STRAINS.into_iter().enumerate() {
            samples.push((
                [alpha, temperature],
                SLOPES[index].mul_add(temperature - 300.0, REFERENCE_YIELD[index]),
            ));
        }
    }
    PropertyClaim {
        key: PropertyKey::new(KIRCHHOFF_YIELD_CURVE_PROPERTY, Pressure::DIMS),
        value: PropertyValue::Surface2D {
            axes: [ALPHA.into(), TEMPERATURE.into()],
            axis_dims: [Dims::NONE, Temperature::DIMS],
            samples,
            triangles: vec![],
            dims: Pressure::DIMS,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity(ALPHA, QuantitySpec::dimensional(Dims::NONE), 0.0, 0.04)
            .with_quantity(
                TEMPERATURE,
                semantic(QuantityKind::AbsoluteTemperature),
                280.0,
                320.0,
            )
            .with_quantity(
                AMBIENT_PRESSURE,
                semantic(QuantityKind::Pressure),
                99_000.0,
                101_000.0,
            ),
        interpolation: InterpolationPolicy::TabulatedOnly,
        uncertainty: UncertaintyModel::HalfWidth {
            half_width: 0.05,
            confidence: 0.95,
        },
        observations: vec![],
        provenance: Provenance {
            source: "synthetic thermal hardening measurements".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn with_original_cells(mut claim: PropertyClaim) -> PropertyClaim {
    let PropertyValue::Surface2D { triangles, .. } = &mut claim.value else {
        panic!("surface fixture");
    };
    for row in 0..TEMPERATURES.len() - 1 {
        for column in 0..STRAINS.len() - 1 {
            let lower = (row * STRAINS.len() + column) as u32;
            let upper = lower + STRAINS.len() as u32;
            triangles.extend([[lower, lower + 1, upper + 1], [lower, upper + 1, upper]]);
        }
    }
    claim.interpolation = InterpolationPolicy::LinearInside;
    claim
}
fn plane() -> PropertyClaim {
    let mut claim = curve();
    let PropertyValue::Surface2D { samples, .. } = &mut claim.value else {
        panic!("surface fixture");
    };
    for (point, value) in samples {
        *value = (-1.0_f64 / 64.0).mul_add(point[1] - 300.0, 50.0_f64.mul_add(point[0], 3.0));
    }
    with_original_cells(claim)
}
fn card(properties: Vec<PropertyClaim>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for property in properties {
        claims.insert_claim(property).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic alloy".into(),
            phase: "solid".into(),
            process: "explicit affine thermal hardening fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn original_grid_retains_its_receipts_and_supplies_signed_hardening_entropy() {
    let original = curve();
    let material = card(vec![original.clone()]);
    let resolved =
        resolve_thermal_tabulated_hardening(&material, &context(300.0), &request(300.0), SINGLE)
            .unwrap();
    let source = resolved.source();
    assert_eq!(source.card_identity(), material.content_hash());
    assert_eq!(source.material_state(), material.id());
    assert_eq!(source.context(), &context(300.0));
    assert_eq!(source.request(), &request(300.0));
    assert_eq!(source.selection(), SINGLE);
    assert_eq!(source.claim(), &original);
    assert!(source.continuous_support().is_none());
    material
        .claims()
        .verify_receipt(&source.selection_answer().receipt)
        .unwrap();
    let PropertyValue::Surface2D { samples, .. } = &original.value else {
        panic!("surface fixture");
    };
    assert_eq!(source.observations().len(), samples.len());
    for ((point, answer), (coordinates, stress)) in source.observations().iter().zip(samples) {
        assert_eq!(point.axes()[ALPHA], coordinates[0]);
        assert_eq!(point.axes()[TEMPERATURE], coordinates[1]);
        assert_eq!(point.axes()[AMBIENT_PRESSURE], 100_000.0);
        assert_eq!(
            point.axis_quantities()[TEMPERATURE],
            semantic(QuantityKind::AbsoluteTemperature)
        );
        assert_eq!(answer.evidence.value.value, *stress);
        assert_eq!(
            answer.receipt.selected,
            source.selection_answer().receipt.selected
        );
        material.claims().verify_receipt(&answer.receipt).unwrap();
    }
    assert_eq!(resolved.hardening().yield_stress_slopes_pa_k(), &SLOPES);
    assert_eq!(
        resolved.hardening().reference().maximum_plastic_strain(),
        0.04
    );
    near(
        resolved
            .hardening()
            .reference()
            .yield_stress_pa(0.03)
            .unwrap(),
        6.0,
        1e-14,
    );
    near(
        resolved
            .hardening()
            .reference()
            .hardening_energy_j_m3(0.03)
            .unwrap(),
        0.035,
        1e-16,
    );
    // Three original segments contribute -1/12800, -1/6400 and -3/12800.
    let derivative = resolved
        .hardening()
        .hardening_energy_temperature_derivative_j_m3_k(0.03)
        .unwrap();
    near(derivative, -0.00046875, 2e-18);
    let warmer =
        resolve_thermal_tabulated_hardening(&material, &context(310.0), &request(310.0), SINGLE)
            .unwrap();
    near(
        warmer
            .hardening()
            .reference()
            .hardening_energy_j_m3(0.03)
            .unwrap()
            - resolved
                .hardening()
                .reference()
                .hardening_energy_j_m3(0.03)
                .unwrap(),
        10.0 * derivative,
        1e-16,
    );
    assert!(
        resolved
            .hardening()
            .hardening_energy_temperature_derivative_j_m3_k(0.040001)
            .is_err()
    );
    let interpolated = context(310.0)
        .with_quantity(ALPHA, QuantitySpec::dimensional(Dims::NONE), 0.03)
        .unwrap();
    assert!(
        material
            .claims()
            .query_typed(
                &original.key,
                &interpolated,
                SelectionPolicy::SingleClaimOnly
            )
            .is_err(),
        "native constitutive interpolation must preserve the original discrete source policy"
    );
}

#[test]
fn swapped_axes_and_narrower_temperature_span_preserve_every_original_vertex() {
    let ordinary = card(vec![curve()]);
    let baseline =
        resolve_thermal_tabulated_hardening(&ordinary, &context(310.0), &request(310.0), SINGLE)
            .unwrap();
    let mut swapped = curve();
    let PropertyValue::Surface2D {
        axes,
        axis_dims,
        samples,
        ..
    } = &mut swapped.value
    else {
        panic!("surface fixture");
    };
    axes.swap(0, 1);
    axis_dims.swap(0, 1);
    for (point, _) in samples.iter_mut() {
        point.swap(0, 1);
    }
    samples.reverse();
    let material = card(vec![swapped.clone()]);
    let mut narrower = request(310.0);
    narrower.temperature_range_k = [300.0, 320.0];
    let resolved =
        resolve_thermal_tabulated_hardening(&material, &context(310.0), &narrower, SINGLE).unwrap();
    assert_eq!(resolved.hardening(), baseline.hardening());
    assert_eq!(resolved.source().claim(), &swapped);
    assert_eq!(resolved.source().request(), &narrower);
    let PropertyValue::Surface2D { samples, .. } = &swapped.value else {
        panic!("surface fixture");
    };
    for ((point, answer), (original, stress)) in
        resolved.source().observations().iter().zip(samples)
    {
        assert_eq!(point.axes()[TEMPERATURE], original[0]);
        assert_eq!(point.axes()[ALPHA], original[1]);
        assert_eq!(answer.evidence.value.value, *stress);
    }
    assert_eq!(resolved.source().observations().len(), 12);
    assert!(
        resolved
            .source()
            .observations()
            .iter()
            .any(|(point, _)| point.axes()[TEMPERATURE] == 280.0)
    );
}

#[test]
fn every_original_temperature_row_and_full_grid_are_required_with_explicit_budgets() {
    let material = card(vec![curve()]);
    for (samples, knots) in [(11, 4), (12, 3), (0, 0)] {
        let mut bounded = request(300.0);
        bounded.maximum_samples = samples;
        bounded.maximum_knots = knots;
        assert!(
            resolve_thermal_tabulated_hardening(&material, &context(300.0), &bounded, SINGLE)
                .is_err()
        );
    }
    for interval in [[290.0, 320.0], [280.0, 330.0], [300.0, 300.0]] {
        let mut unsupported = request(300.0);
        unsupported.temperature_range_k = interval;
        assert!(
            resolve_thermal_tabulated_hardening(&material, &context(300.0), &unsupported, SINGLE)
                .is_err()
        );
    }
    let mut missing = curve();
    let PropertyValue::Surface2D { samples, .. } = &mut missing.value else {
        panic!("surface fixture");
    };
    samples.remove(5);
    assert!(
        resolve_thermal_tabulated_hardening(
            &card(vec![missing]),
            &context(300.0),
            &request(300.0),
            SINGLE
        )
        .is_err()
    );

    // The endpoint grid is unchanged. An interior-temperature kink cannot be
    // discarded: it would change the hardening-storage temperature derivative.
    let mut kinked = curve();
    let PropertyValue::Surface2D { samples, .. } = &mut kinked.value else {
        panic!("surface fixture");
    };
    samples[5].1 += 0.125;
    assert!(matches!(
        resolve_thermal_tabulated_hardening(
            &card(vec![kinked]),
            &context(300.0),
            &request(300.0),
            SINGLE
        ),
        Err(ThermalHardeningSourceError::Invalid(_))
    ));

    let mut single_temperature = curve();
    let PropertyValue::Surface2D { samples, .. } = &mut single_temperature.value else {
        panic!("surface fixture");
    };
    samples.retain(|(point, _)| point[1] == 300.0);
    let mut no_derivative = request(300.0);
    no_derivative.temperature_range_k = [300.0, 320.0];
    assert!(
        resolve_thermal_tabulated_hardening(
            &card(vec![single_temperature]),
            &context(300.0),
            &no_derivative,
            SINGLE
        )
        .is_err()
    );

    let mut constant = curve();
    let PropertyValue::Surface2D { samples, .. } = &mut constant.value else {
        panic!("surface fixture");
    };
    for (index, (_, value)) in samples.iter_mut().enumerate() {
        *value = REFERENCE_YIELD[index % STRAINS.len()];
    }
    let measured_constant = resolve_thermal_tabulated_hardening(
        &card(vec![constant]),
        &context(300.0),
        &request(300.0),
        SINGLE,
    )
    .unwrap();
    assert_eq!(
        measured_constant.hardening().yield_stress_slopes_pa_k(),
        &[0.0; 4]
    );
    assert_eq!(
        measured_constant
            .hardening()
            .hardening_energy_temperature_derivative_j_m3_k(0.03)
            .unwrap(),
        0.0
    );
}

#[test]
fn original_linear_cells_must_match_the_model_and_cover_the_whole_rectangle() {
    let original = plane();
    let material = card(vec![original.clone()]);
    let resolved =
        resolve_thermal_tabulated_hardening(&material, &context(300.0), &request(300.0), SINGLE)
            .unwrap();
    let envelope = resolved
        .source()
        .continuous_support()
        .expect("original cell support");
    material
        .claims()
        .verify_receipt(&envelope.lower.receipt)
        .unwrap();
    material
        .claims()
        .verify_receipt(&envelope.upper.receipt)
        .unwrap();
    near(
        resolved
            .hardening()
            .reference()
            .yield_stress_pa(0.015)
            .unwrap(),
        3.75,
        1e-14,
    );
    near(
        resolved
            .hardening()
            .hardening_energy_temperature_derivative_j_m3_k(0.03)
            .unwrap(),
        0.0,
        1e-18,
    );

    let mut bilinear = curve();
    let PropertyValue::Surface2D { samples, .. } = &mut bilinear.value else {
        panic!("surface fixture");
    };
    for (index, (point, value)) in samples.iter_mut().enumerate() {
        *value = SLOPES[index % STRAINS.len()]
            .mul_add(point[1] - 300.0, 50.0_f64.mul_add(point[0], 3.0));
    }
    assert!(
        resolve_thermal_tabulated_hardening(
            &card(vec![with_original_cells(bilinear)]),
            &context(300.0),
            &request(300.0),
            SINGLE
        )
        .is_err(),
        "triangular source cells must not become a different bilinear law"
    );

    let mut hole = plane();
    let PropertyValue::Surface2D { triangles, .. } = &mut hole.value else {
        panic!("surface fixture");
    };
    triangles.remove(4);
    assert!(matches!(
        resolve_thermal_tabulated_hardening(
            &card(vec![hole]),
            &context(300.0),
            &request(300.0),
            SINGLE
        ),
        Err(ThermalHardeningSourceError::Support(_))
    ));

    let mut conflicting_policy = original;
    conflicting_policy.interpolation = InterpolationPolicy::TabulatedOnly;
    assert!(
        resolve_thermal_tabulated_hardening(
            &card(vec![conflicting_policy]),
            &context(300.0),
            &request(300.0),
            SINGLE
        )
        .is_err(),
        "a discrete observation grid cannot also declare competing triangular interpolation"
    );
}

#[test]
fn explicit_claim_pins_preserve_fixed_conditions_and_exact_kinematic_meaning() {
    let first = curve();
    let mut competing = first.clone();
    competing.provenance.source = "a different synthetic experiment".into();
    let pin = ClaimSelection::Pinned(ClaimId(first.content_hash()));
    let material = card(vec![first.clone(), competing]);
    assert!(
        resolve_thermal_tabulated_hardening(&material, &context(300.0), &request(300.0), SINGLE)
            .is_err()
    );
    let resolved =
        resolve_thermal_tabulated_hardening(&material, &context(300.0), &request(300.0), pin)
            .unwrap();
    assert_eq!(resolved.source().claim(), &first);
    assert_eq!(resolved.source().selection(), pin);
    assert!(
        resolved
            .source()
            .observations()
            .iter()
            .all(|(_, answer)| answer.receipt.selected == ClaimId(first.content_hash()))
    );

    let wrong_fixed_pressure = context(300.0).with(AMBIENT_PRESSURE, 102_000.0).unwrap();
    assert!(
        resolve_thermal_tabulated_hardening(&material, &wrong_fixed_pressure, &request(300.0), pin)
            .is_err()
    );
    let frozen_rate = context(300.0).with("strain-rate", 1.0).unwrap();
    assert!(
        resolve_thermal_tabulated_hardening(&material, &frozen_rate, &request(300.0), pin).is_err()
    );
    let duplicate_temperature = context(300.0).with("T", 300.0).unwrap();
    assert!(
        resolve_thermal_tabulated_hardening(
            &material,
            &duplicate_temperature,
            &request(300.0),
            pin
        )
        .is_err()
    );
    let wrong_temperature_meaning = QueryPoint::new()
        .with_quantity(
            TEMPERATURE,
            QuantitySpec::dimensional(Temperature::DIMS),
            300.0,
        )
        .unwrap();
    assert!(
        resolve_thermal_tabulated_hardening(
            &material,
            &wrong_temperature_meaning,
            &request(300.0),
            pin
        )
        .is_err()
    );

    let mut cauchy = first.clone();
    cauchy.key = PropertyKey::new(CAUCHY_YIELD_CURVE_PROPERTY, Pressure::DIMS);
    assert!(
        resolve_thermal_tabulated_hardening(
            &card(vec![cauchy]),
            &context(300.0),
            &request(300.0),
            SINGLE
        )
        .is_err()
    );
    let mut total_strain = first;
    let PropertyValue::Surface2D { axes, .. } = &mut total_strain.value else {
        panic!("surface fixture");
    };
    axes[0] = "total-strain".into();
    total_strain.validity = total_strain.validity.with_quantity(
        "total-strain",
        QuantitySpec::dimensional(Dims::NONE),
        0.0,
        0.04,
    );
    assert!(
        resolve_thermal_tabulated_hardening(
            &card(vec![total_strain]),
            &context(300.0),
            &request(300.0),
            SINGLE
        )
        .is_err()
    );
}
