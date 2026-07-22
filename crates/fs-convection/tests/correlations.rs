//! G0/G3 correlation-card, refusal, typing, and metamorphic coverage.

use fs_convection::{
    CorrelationError, CorrelationId, CorrelationInputs, DiscrepancyBasis, HeatTransferCoefficient,
    ThermalConductivity, ThermalDirection, correlation_catalog, evaluate,
};
use fs_evidence::{CertifyError, NumericalKind};
use fs_qty::Length;
use fs_vvreg::thermal_level_a::{
    ThermalLevelAAcceptance, ThermalLevelACase, ThermalLevelAFamily, ThermalLevelAKind,
    thermal_level_a_cases,
};

const LEVEL_A_CONVECTION_BINDINGS: [(&str, &str); 2] = [
    (
        "thermal-a-duct-nu-cwt",
        "tests/correlations.rs::level_a_fully_developed_limits_and_rectangular_endpoints_are_frozen",
    ),
    (
        "thermal-a-duct-nu-chf",
        "tests/correlations.rs::level_a_fully_developed_limits_and_rectangular_endpoints_are_frozen",
    ),
];

fn close(observed: f64, expected: f64, relative: f64) {
    let scale = expected.abs().max(1.0);
    assert!(
        (observed - expected).abs() <= relative * scale,
        "observed={observed:.17e} expected={expected:.17e} tolerance={relative:.3e}"
    );
}

fn duct_inputs() -> CorrelationInputs {
    CorrelationInputs::forced(1_000.0, 0.7).with_length_ratio(100.0)
}

fn level_a_convection_case(case_id: &str) -> &'static ThermalLevelACase {
    let case = thermal_level_a_cases()
        .iter()
        .find(|case| case.id == case_id)
        .unwrap_or_else(|| panic!("missing Level-A case {case_id}"));
    assert_eq!(case.family, ThermalLevelAFamily::ConvectionLimit);
    assert_eq!(case.kind, ThermalLevelAKind::AnalyticReference);
    assert_eq!(case.metric, "nusselt-number");
    assert!(
        LEVEL_A_CONVECTION_BINDINGS
            .iter()
            .any(|(id, _)| *id == case_id),
        "{case_id} is not declared as an executing fs-convection binding"
    );
    let reynolds = case
        .context
        .iter()
        .find(|axis| axis.name == "reynolds-number")
        .unwrap_or_else(|| panic!("{case_id} must declare its Reynolds context"));
    assert!(reynolds.lo <= 1_000.0 && 1_000.0 <= reynolds.hi);
    case
}

fn assert_level_a_convection_limit(case_id: &str, observed: f64) {
    let case = level_a_convection_case(case_id);
    let ThermalLevelAAcceptance::Tolerance { atol, rtol } = case.acceptance else {
        panic!("{case_id}: analytic Level-A row must carry a scalar tolerance");
    };
    let absolute_error = (observed - case.reference_value_si).abs();
    let envelope = atol + rtol * case.reference_value_si.abs();
    assert!(
        absolute_error <= envelope,
        "{case_id}: Nu={observed:.17e}, reference={:.17e}, error={absolute_error:.3e}, \
         envelope={envelope:.3e}",
        case.reference_value_si
    );
    assert_eq!(observed.to_bits(), case.reference_value_si.to_bits());
    println!(
        "{{\"suite\":\"fs-convection/level-a\",\"case_id\":\"{case_id}\",\
         \"computed\":{observed:.17e},\"reference\":{:.17e},\
         \"absolute_error\":{absolute_error:.17e},\"envelope\":{envelope:.17e},\
         \"authority\":\"executed-formula-limit-not-registry-receipt\",\
         \"verdict\":\"pass\"}}",
        case.reference_value_si
    );
}

#[test]
fn catalog_has_eleven_sourced_cards_and_no_unlabeled_discrepancy() {
    let catalog = correlation_catalog();
    assert_eq!(catalog.len(), 11);
    assert_eq!(catalog.len(), CorrelationId::ALL.len());

    let mut names = std::collections::BTreeSet::new();
    for card in catalog {
        assert!(
            names.insert(card.id.name()),
            "duplicate card {}",
            card.id.name()
        );
        assert_eq!(card.model.name, card.id.name());
        assert!(!card.model.validity.bounds().is_empty());
        assert!(!card.source.citation.trim().is_empty());
        assert!(!card.source.identifier.trim().is_empty());
        assert!(!card.model.assumptions.is_empty());
        assert!(!card.model.known_failures.is_empty());
        match card.discrepancy_basis {
            DiscrepancyBasis::AnalyticIdealLimit => {
                assert_eq!(card.model.discrepancy_rel.to_bits(), 0);
            }
            DiscrepancyBasis::EngineeringAllowance => {
                assert!(card.model.discrepancy_rel >= 0.10 && card.model.discrepancy_rel <= 0.25);
            }
        }
    }
}

#[test]
fn level_a_fully_developed_limits_and_rectangular_endpoints_are_frozen() {
    let cwt = evaluate(CorrelationId::CircularDuctLaminarCwt, duct_inputs()).expect("CWT");
    let chf = evaluate(CorrelationId::CircularDuctLaminarChf, duct_inputs()).expect("CHF");
    assert_level_a_convection_limit("thermal-a-duct-nu-cwt", cwt.evidence().value);
    assert_level_a_convection_limit("thermal-a-duct-nu-chf", chf.evidence().value);

    let square = CorrelationInputs::forced(1_000.0, 0.7)
        .with_length_ratio(100.0)
        .with_aspect_ratio(1.0);
    close(
        evaluate(CorrelationId::RectangularDuctLaminarCwt, square)
            .expect("square CWT")
            .evidence()
            .value,
        2.978_695,
        2.0e-15,
    );
    close(
        evaluate(CorrelationId::RectangularDuctLaminarChf, square)
            .expect("square CHF")
            .evidence()
            .value,
        3.610_224,
        2.0e-15,
    );
}

#[test]
fn level_a_convection_binding_partition_is_complete() {
    let catalog_ids = thermal_level_a_cases()
        .iter()
        .filter(|case| case.family == ThermalLevelAFamily::ConvectionLimit)
        .map(|case| case.id)
        .collect::<std::collections::BTreeSet<_>>();
    let binding_ids = LEVEL_A_CONVECTION_BINDINGS
        .iter()
        .map(|(id, _)| *id)
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(binding_ids, catalog_ids);
    assert_eq!(binding_ids.len(), 2);
    for (id, test) in LEVEL_A_CONVECTION_BINDINGS {
        assert!(
            test.starts_with("tests/correlations.rs::"),
            "{id}: executing test path must be stable"
        );
    }
}

#[test]
fn every_nonconstant_formula_has_a_frozen_source_formula_spot_value() {
    let cases = [
        (
            CorrelationId::CircularDuctHausen,
            CorrelationInputs::forced(1_000.0, 7.0).with_length_ratio(20.0),
            11.488_360_610_697_356,
        ),
        (
            CorrelationId::RectangularDuctLaminarCwt,
            CorrelationInputs::forced(1_000.0, 0.7)
                .with_length_ratio(100.0)
                .with_aspect_ratio(0.5),
            3.388_736_875_000_000_6,
        ),
        (
            CorrelationId::RectangularDuctLaminarChf,
            CorrelationInputs::forced(1_000.0, 0.7)
                .with_length_ratio(100.0)
                .with_aspect_ratio(0.5),
            4.125_812_203_124_999,
        ),
        (
            CorrelationId::DittusBoelter,
            CorrelationInputs::forced(100_000.0, 0.7).with_length_ratio(100.0),
            199.419_237_807_658_48,
        ),
        (
            CorrelationId::Gnielinski,
            CorrelationInputs::forced(100_000.0, 0.7).with_length_ratio(100.0),
            178.622_951_779_291_2,
        ),
        (
            CorrelationId::FlatPlateLaminarAverage,
            CorrelationInputs::forced(100_000.0, 0.7),
            186.437_852_875_226_2,
        ),
        (
            CorrelationId::FlatPlateTurbulentAverage,
            CorrelationInputs::forced(1_000_000.0, 0.7),
            1_299.484_953_525_734_2,
        ),
        (
            CorrelationId::ChurchillBernsteinCylinder,
            CorrelationInputs::forced(10_000.0, 0.7),
            53.327_788_670_209_97,
        ),
        (
            CorrelationId::ChurchillChuVerticalPlate,
            CorrelationInputs::natural(1.0e6, 0.7),
            16.530_366_876_407_225,
        ),
    ];

    for (id, inputs, expected) in cases {
        let observed = evaluate(id, inputs).unwrap_or_else(|error| panic!("{id:?}: {error}"));
        close(observed.evidence().value, expected, 3.0e-13);
    }
}

#[test]
fn validity_edges_are_inclusive_and_missing_or_outside_axes_refuse() {
    evaluate(
        CorrelationId::DittusBoelter,
        CorrelationInputs::forced(10_000.0, 0.7).with_length_ratio(10.0),
    )
    .expect("inclusive lower edge");
    evaluate(
        CorrelationId::DittusBoelter,
        CorrelationInputs::forced(120_000.0, 120.0).with_length_ratio(1.0e6),
    )
    .expect("inclusive upper edge");

    let missing = evaluate(
        CorrelationId::DittusBoelter,
        CorrelationInputs::forced(100_000.0, 0.7),
    )
    .expect_err("L/Dh is mandatory");
    match missing {
        CorrelationError::OutOfDomain { violations, .. } => {
            assert_eq!(violations.len(), 1);
            assert_eq!(violations[0].axis, "L_over_Dh");
            assert_eq!(violations[0].value, None);
        }
        other => panic!("unexpected refusal: {other}"),
    }

    let outside = evaluate(
        CorrelationId::DittusBoelter,
        CorrelationInputs::forced(9_999.0, 0.69).with_length_ratio(9.0),
    )
    .expect_err("three axes are outside");
    match outside {
        CorrelationError::OutOfDomain { violations, .. } => {
            let axes = violations
                .iter()
                .map(|violation| violation.axis.as_str())
                .collect::<Vec<_>>();
            assert_eq!(axes, ["L_over_Dh", "Pr", "Re"]);
        }
        other => panic!("unexpected refusal: {other}"),
    }

    assert!(matches!(
        evaluate(
            CorrelationId::FlatPlateLaminarAverage,
            CorrelationInputs::forced(f64::NAN, 0.7)
        ),
        Err(CorrelationError::InvalidGroup { axis: "Re", .. })
    ));
}

#[test]
fn cylinder_card_checks_the_product_constraint_through_a_named_peclet_axis() {
    let accepted = evaluate(
        CorrelationId::ChurchillBernsteinCylinder,
        CorrelationInputs::forced(1.0, 0.2),
    )
    .expect("Pe=0.2 is inclusive");
    assert_eq!(accepted.groups().get("Pe"), Some(&0.2));

    let refused = evaluate(
        CorrelationId::ChurchillBernsteinCylinder,
        CorrelationInputs::forced(1.0, 0.199),
    )
    .expect_err("Pe and Pr are outside");
    let CorrelationError::OutOfDomain { violations, .. } = refused else {
        panic!("expected domain refusal");
    };
    assert!(violations.iter().any(|violation| violation.axis == "Pe"));
    assert!(violations.iter().any(|violation| violation.axis == "Pr"));
}

#[test]
fn typed_h_retains_model_evidence_and_empirical_result_cannot_certify() {
    fn require_htc(_: &fs_evidence::Evidence<HeatTransferCoefficient>) {}

    let evaluated = evaluate(
        CorrelationId::Gnielinski,
        CorrelationInputs::forced(100_000.0, 0.7).with_length_ratio(100.0),
    )
    .expect("in domain");
    let coefficient = evaluated
        .heat_transfer_coefficient(ThermalConductivity::new(0.026), Length::new(0.010))
        .expect("typed h");
    require_htc(&coefficient);
    close(
        coefficient.value.value(),
        evaluated.evidence().value * 0.026 / 0.010,
        1.0e-15,
    );
    assert_eq!(coefficient.numerical.kind, NumericalKind::Estimate);
    assert_eq!(coefficient.model.cards, [CorrelationId::Gnielinski.name()]);
    assert!(coefficient.model.in_domain);
    assert_eq!(coefficient.model.discrepancy_rel, 0.10);
    assert!(matches!(
        coefficient.certified(),
        Err(CertifyError::NotRigorous {
            kind: NumericalKind::Estimate
        })
    ));
}

#[test]
fn g3_coherent_unit_rescaling_leaves_nu_and_h_invariant() {
    let evaluated = evaluate(
        CorrelationId::FlatPlateLaminarAverage,
        CorrelationInputs::forced(100_000.0, 0.7),
    )
    .expect("in domain");
    let si = evaluated
        .heat_transfer_coefficient(ThermalConductivity::new(0.026), Length::new(0.1))
        .expect("SI");
    // The same k supplied as 26 mW/(m K), and L as 100 mm, normalized
    // explicitly at the boundary before entering the coherent-SI types.
    let rescaled = evaluated
        .heat_transfer_coefficient(
            ThermalConductivity::new(26.0 * 1.0e-3),
            Length::new(100.0 * 1.0e-3),
        )
        .expect("rescaled");
    close(si.value.value(), rescaled.value.value(), 1.0e-15);
    assert_eq!(si.model, rescaled.model);
}

#[test]
fn dittus_boelter_direction_is_semantic_and_provenance_bearing() {
    let heating = evaluate(
        CorrelationId::DittusBoelter,
        CorrelationInputs::forced(100_000.0, 0.7)
            .with_length_ratio(100.0)
            .with_direction(ThermalDirection::HeatingFluid),
    )
    .expect("heating");
    let cooling = evaluate(
        CorrelationId::DittusBoelter,
        CorrelationInputs::forced(100_000.0, 0.7)
            .with_length_ratio(100.0)
            .with_direction(ThermalDirection::CoolingFluid),
    )
    .expect("cooling");
    close(cooling.evidence().value, 206.660_391_611_847_25, 3.0e-13);
    assert_ne!(
        heating.evidence().value.to_bits(),
        cooling.evidence().value.to_bits()
    );
    assert_ne!(heating.evidence().provenance, cooling.evidence().provenance);
}

#[test]
fn dimensional_inputs_refuse_zero_negative_nan_and_overflow() {
    let evaluated = evaluate(
        CorrelationId::FlatPlateLaminarAverage,
        CorrelationInputs::forced(100_000.0, 0.7),
    )
    .expect("in domain");
    for conductivity in [0.0, -1.0, f64::NAN] {
        assert!(matches!(
            evaluated.heat_transfer_coefficient(
                ThermalConductivity::new(conductivity),
                Length::new(0.1)
            ),
            Err(CorrelationError::InvalidDimensionalInput {
                field: "fluid thermal conductivity",
                ..
            })
        ));
    }
    for length in [0.0, -1.0, f64::INFINITY] {
        assert!(matches!(
            evaluated
                .heat_transfer_coefficient(ThermalConductivity::new(0.026), Length::new(length)),
            Err(CorrelationError::InvalidDimensionalInput {
                field: "characteristic length",
                ..
            })
        ));
    }
    assert!(matches!(
        evaluated.heat_transfer_coefficient(
            ThermalConductivity::new(f64::MAX),
            Length::new(f64::MIN_POSITIVE)
        ),
        Err(CorrelationError::NonFiniteResult {
            stage: "Nu-to-h conversion",
            ..
        })
    ));
}
