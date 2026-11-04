//! Implementation checks against independent high-precision evaluation of
//! Wright (1996), Eqs. (9a-c) and (7); these are not experimental validation.

use fs_convection::{
    CorrelationError, CorrelationId, CorrelationInputs, DiscrepancyBasis, EvaluationSourceRegion,
    ThermalConductivity, evaluate,
};
use fs_qty::Length;

const CARD: CorrelationId = CorrelationId::WrightVerticalAirLayer;

fn inputs(rayleigh: f64, height_over_gap: f64) -> CorrelationInputs {
    CorrelationInputs::vertical_air_layer(rayleigh, 0.71, height_over_gap)
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual / expected - 1.0).abs() < 3e-13,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

fn nu(rayleigh: f64, height_over_gap: f64) -> f64 {
    evaluate(CARD, inputs(rayleigh, height_over_gap))
        .unwrap()
        .evidence()
        .value
}

#[test]
fn independent_reference_values_cover_each_source_branch_and_finite_height() {
    // Independent 80-digit Decimal arithmetic of the published equations.
    for (ra, aspect, expected) in [
        (100.0, 40.0, 1.000_006_956_362_895_7),
        (1000.0, 40.0, 1.001_383_113_212_262_8),
        (5000.0, 40.0, 1.055_901_461_222_033_2),
        (10_000.0, 40.0, 1.275_000_339_485_781_6),
        (20_000.0, 40.0, 1.688_829_888_799_044_3),
        (50_000.0, 40.0, 2.466_574_890_260_911),
        (100_000.0, 40.0, 3.127_678_936_463_967_7),
        (1_000_000.0, 40.0, 6.738_38),
        (5000.0, 5.0, 1.584_219_541_242_805_5),
        (20_000.0, 5.0, 2.309_806_949_236_508_5),
    ] {
        near(nu(ra, aspect), expected);
    }
    // A vertical cavity does not have the horizontal layer's Ra=1708 onset.
    assert!(nu(1000.0, 40.0) > 1.0);
    let sample = evaluate(CARD, inputs(20_000.0, 40.0)).unwrap();
    assert!(sample.card().source.citation.contains("Wright"));
    assert!(sample.card().source.identifier.contains("AT-96-16-2"));
    assert_eq!(sample.source_region(), EvaluationSourceRegion::CitedFormula);
    assert_eq!(
        sample.card().discrepancy_basis,
        DiscrepancyBasis::EngineeringAllowance
    );
    assert_eq!(sample.groups()["aspect_ratio"], 40.0);
}

#[test]
fn zero_driving_is_total_conduction_and_typed_lowering_counts_the_gap_once() {
    for ra in [0.0, 1e-300] {
        let sample = evaluate(CARD, inputs(ra, 40.0)).unwrap();
        assert_eq!(sample.evidence().value, 1.0);
        assert_eq!(
            sample
                .heat_transfer_coefficient(ThermalConductivity::new(0.028), Length::new(0.02))
                .unwrap()
                .value
                .value(),
            0.028 / 0.02
        );
    }
    let sample = evaluate(CARD, inputs(20_000.0, 40.0)).unwrap();
    let h = sample
        .heat_transfer_coefficient(ThermalConductivity::new(0.028), Length::new(0.02))
        .unwrap();
    let wider_gap = sample
        .heat_transfer_coefficient(ThermalConductivity::new(0.028), Length::new(0.04))
        .unwrap();
    near(h.value.value(), 1.688_829_888_799_044_3 * 0.028 / 0.02);
    assert_eq!(wider_gap.value.value(), h.value.value() * 0.5);
    assert_eq!(h.model.cards, [CARD.name()]);
}

#[test]
fn height_correction_selects_the_maximum_without_adding_a_second_heat_path() {
    for ra in [100.0, 1000.0, 5000.0, 20_000.0, 100_000.0, 1_000_000.0] {
        assert!(nu(ra, 5.0) >= nu(ra, 10.0));
        assert!(nu(ra, 10.0) >= nu(ra, 40.0));
        assert_eq!(nu(ra, 40.0).to_bits(), nu(ra, 80.0).to_bits());
        assert_eq!(nu(ra, 40.0).to_bits(), nu(ra, 110.0).to_bits());
    }
    assert!(nu(5000.0, 5.0) > 1.4 * nu(5000.0, 40.0));
    assert_ne!(
        evaluate(CARD, inputs(5000.0, 5.0)).unwrap().evidence().provenance,
        evaluate(CARD, inputs(5000.0, 40.0)).unwrap().evidence().provenance
    );
}

#[test]
fn analytic_partials_cover_tall_and_finite_height_branches_without_low_ra_cancellation() {
    for (ra, expected) in [
        (100.0, 1.598_902_968_514_361_2e-7),
        (5000.0, 2.569_762_780_660_87e-5),
        (20_000.0, 3.490_811_380_147_624_3e-5),
        (100_000.0, 1.042_559_645_487_989_3e-5),
    ] {
        let sample = inputs(ra, 40.0).evaluate_with_sensitivities(CARD).unwrap();
        near(sample.evidence().sensitivity.d_qoi["Ra"], expected);
        assert_eq!(sample.evidence().sensitivity.d_qoi["aspect_ratio"], 0.0);
        assert_eq!(sample.evidence().sensitivity.d_qoi["Pr"], 0.0);
    }
    let tiny = inputs(1e-10, 40.0).evaluate_with_sensitivities(CARD).unwrap();
    assert_eq!(tiny.evidence().value, 1.0);
    assert!(tiny.evidence().sensitivity.d_qoi["Ra"] > 0.0);

    let sample = inputs(5000.0, 6.0).evaluate_with_sensitivities(CARD).unwrap();
    for (axis, step, plus, minus) in [
        ("Ra", 0.05, nu(5000.05, 6.0), nu(4999.95, 6.0)),
        ("aspect_ratio", 0.00006, nu(5000.0, 6.00006), nu(5000.0, 5.99994)),
    ] {
        let finite_difference = (plus - minus) / (2.0 * step);
        assert!(
            (sample.evidence().sensitivity.d_qoi[axis] / finite_difference - 1.0).abs()
                < 2e-8
        );
    }
    let h = sample
        .heat_transfer_coefficient_with_sensitivities(
            ThermalConductivity::new(0.028),
            Length::new(0.02),
        )
        .unwrap();
    near(
        h.sensitivity.d_qoi["aspect_ratio"],
        (0.028 / 0.02) * sample.evidence().sensitivity.d_qoi["aspect_ratio"],
    );
}

#[test]
fn source_breaks_and_domain_edges_do_not_publish_cross_branch_derivatives() {
    for ra in [0.0, 10_000.0, 50_000.0, 1_000_000.0] {
        let sample = inputs(ra, 40.0).evaluate_with_sensitivities(CARD).unwrap();
        assert!(!sample.evidence().sensitivity.d_qoi.contains_key("Ra"));
        assert_eq!(sample.evidence().sensitivity.d_qoi["Pr"], 0.0);
    }
    for aspect in [5.0, 110.0] {
        let sample = inputs(5000.0, aspect).evaluate_with_sensitivities(CARD).unwrap();
        assert!(!sample.evidence().sensitivity.d_qoi.contains_key("aspect_ratio"));
        assert!(sample.evidence().sensitivity.d_qoi.contains_key("Ra"));
    }
}

#[test]
fn missing_geometry_invalid_groups_and_unsupported_domains_refuse() {
    for (ra, pr, aspect) in [
        (1_000_000.0_f64.next_up(), 0.71, 40.0),
        (10_000.0, 0.70_f64.next_down(), 40.0),
        (10_000.0, 0.72_f64.next_up(), 40.0),
        (10_000.0, 0.71, 5.0_f64.next_down()),
        (10_000.0, 0.71, 110.0_f64.next_up()),
    ] {
        assert!(matches!(
            evaluate(CARD, CorrelationInputs::vertical_air_layer(ra, pr, aspect)),
            Err(CorrelationError::OutOfDomain { .. })
        ));
    }
    for ra in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            evaluate(CARD, inputs(ra, 40.0)),
            Err(CorrelationError::InvalidGroup { .. })
        ));
    }
    assert!(matches!(
        evaluate(CARD, CorrelationInputs::natural(5000.0, 0.71)),
        Err(CorrelationError::OutOfDomain { .. })
    ));
    assert!(evaluate(CARD, inputs(0.0, 5.0)).is_ok());
    assert!(evaluate(CARD, inputs(1_000_000.0, 110.0)).is_ok());
}
