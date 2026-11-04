//! G0/G3 horizontal gas-layer source, stability, domain and derivative checks.
//! Formula arithmetic verifies implementation; it is not experimental validation.

use fs_convection::{
    CorrelationError, CorrelationId, CorrelationInputs, DiscrepancyBasis, EvaluationSourceRegion,
    ThermalConductivity, evaluate,
};
use fs_qty::Length;

const CARD: CorrelationId = CorrelationId::HollandsHorizontalAirLayer;

fn inputs(rayleigh: f64) -> CorrelationInputs {
    CorrelationInputs::horizontal_air_layer(rayleigh, 0.71, 0.05)
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual / expected - 1.0).abs() < 3e-13,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn horizontal_source_equation_matches_independent_decimal_arithmetic() {
    // Independent 70-digit Python Decimal evaluation of the horizontal
    // reduction of report COO-2597-2, Eq. (4-1), doi:10.2172/7125623.
    // Decimal powers use exponent Decimal(1)/Decimal(3), not fs-math::det.
    for (rayleigh, expected) in [
        (1708.0, 1.0),
        (2000.0, 1.210_24),
        (3416.0, 1.72),
        (5830.0, 2.018_126_929_674_099_5),
        (10_000.0, 2.391_093_013_195_09),
        (50_000.0, 3.437_727_779_652_893_2),
        (100_000.0, 3.994_360_101_957_176),
    ] {
        near(
            evaluate(CARD, inputs(rayleigh)).unwrap().evidence().value,
            expected,
        );
    }
    let sample = evaluate(CARD, inputs(10_000.0)).unwrap();
    assert_eq!(
        sample.card().source.identifier,
        "doi:10.2172/7125623; doi:10.1115/1.3450517"
    );
    assert!(sample.card().source.citation.contains("Eq. (4-1)"));
    assert_eq!(sample.source_region(), EvaluationSourceRegion::CitedFormula);
    assert_eq!(
        sample.card().discrepancy_basis,
        DiscrepancyBasis::EngineeringAllowance
    );
    assert!(
        sample
            .card()
            .model
            .assumptions
            .iter()
            .any(|line| line.contains("not a source-published validation range"))
    );
}

#[test]
fn stable_zero_gravity_and_subcritical_layers_keep_conduction_exactly_once() {
    for rayleigh in [-f64::MAX, -1e12, -10_000.0, -0.0, 0.0, 1e-300, 1708.0] {
        let sample = evaluate(CARD, inputs(rayleigh)).unwrap();
        assert_eq!(sample.groups()["Ra"].to_bits(), rayleigh.to_bits());
        assert_eq!(sample.evidence().value, 1.0);
        let h = sample
            .heat_transfer_coefficient(ThermalConductivity::new(0.028), Length::new(0.02))
            .unwrap();
        assert_eq!(h.value.value(), 0.028 / 0.02);
        // Stable hot-above configuration: a positive total conductance and
        // one complete temperature contrast send heat down into the cold wall.
        let heat_into_lower = h.value.value() * 0.25 * (320.0 - 300.0);
        near(heat_into_lower, 7.0);
    }
    let unstable = evaluate(CARD, inputs(10_000.0)).unwrap();
    assert!(unstable.evidence().value > 1.0);
    assert_ne!(
        unstable.evidence().provenance,
        evaluate(CARD, inputs(-10_000.0))
            .unwrap()
            .evidence()
            .provenance
    );
}

#[test]
fn admission_rejects_missing_geometry_and_unsupported_air_layer_domains() {
    let Err(CorrelationError::OutOfDomain { violations, .. }) =
        evaluate(CARD, CorrelationInputs::natural(10_000.0, 0.71))
    else {
        panic!("whole layer aspect ratio must be supplied");
    };
    assert!(
        violations
            .iter()
            .any(|v| v.axis == "aspect_ratio" && v.value.is_none())
    );
    for (rayleigh, prandtl, aspect) in [
        (100_000.0_f64.next_up(), 0.71, 0.05),
        (0.0, 0.70_f64.next_down(), 0.05),
        (-10_000.0, 0.72_f64.next_up(), 0.05),
        (10_000.0, 0.71, 0.1),
    ] {
        assert!(matches!(
            evaluate(
                CARD,
                CorrelationInputs::horizontal_air_layer(rayleigh, prandtl, aspect)
            ),
            Err(CorrelationError::OutOfDomain { .. })
        ));
    }
    for (rayleigh, prandtl, aspect) in [
        (f64::NAN, 0.71, 0.05),
        (f64::NEG_INFINITY, 0.71, 0.05),
        (f64::INFINITY, 0.71, 0.05),
        (0.0, f64::NAN, 0.05),
        (0.0, 0.71, 0.0),
        (0.0, 0.71, f64::INFINITY),
    ] {
        assert!(matches!(
            evaluate(
                CARD,
                CorrelationInputs::horizontal_air_layer(rayleigh, prandtl, aspect)
            ),
            Err(CorrelationError::InvalidGroup { .. })
        ));
    }
    for prandtl in [0.70, 0.72] {
        assert!(
            evaluate(
                CARD,
                CorrelationInputs::horizontal_air_layer(100_000.0, prandtl, 0.1_f64.next_down(),)
            )
            .is_ok()
        );
    }
    // The new signed axis must not admit negative Rayleigh on external cards.
    assert!(matches!(
        evaluate(
            CorrelationId::ChurchillChuVerticalPlate,
            CorrelationInputs::natural(-10_000.0, 0.71)
        ),
        Err(CorrelationError::InvalidGroup { .. })
    ));
}

#[test]
fn analytic_partials_keep_stable_zero_rayleigh_and_source_kinks_honest() {
    for rayleigh in [-1e12, -10_000.0, -0.0, 0.0, 1000.0] {
        let sample = inputs(rayleigh).evaluate_with_sensitivities(CARD).unwrap();
        for axis in ["Ra", "Pr", "aspect_ratio"] {
            assert_eq!(sample.evidence().sensitivity.d_qoi[axis], 0.0);
        }
    }
    for (rayleigh, expected) in [
        (2000.0, 0.000_614_88),
        (3416.0, 0.000_210_772_833_723_653_4),
        (10_000.0, 0.000_064_496_700_439_836_33),
        (50_000.0, 0.000_014_629_929_197_685_955),
    ] {
        let sample = inputs(rayleigh).evaluate_with_sensitivities(CARD).unwrap();
        near(sample.evidence().sensitivity.d_qoi["Ra"], expected);
        assert_eq!(sample.evidence().sensitivity.d_qoi["Pr"], 0.0);
        assert_eq!(sample.evidence().sensitivity.d_qoi["aspect_ratio"], 0.0);
        assert_eq!(
            sample.evidence().value.to_bits(),
            evaluate(CARD, inputs(rayleigh))
                .unwrap()
                .evidence()
                .value
                .to_bits()
        );
    }
    for rayleigh in [1708.0, 5830.0, 100_000.0] {
        let sample = inputs(rayleigh).evaluate_with_sensitivities(CARD).unwrap();
        assert!(!sample.evidence().sensitivity.d_qoi.contains_key("Ra"));
        assert_eq!(sample.evidence().sensitivity.d_qoi["Pr"], 0.0);
        assert_eq!(sample.evidence().sensitivity.d_qoi["aspect_ratio"], 0.0);
    }
}

#[test]
fn source_branches_are_continuous_monotone_and_preserve_typed_scaling() {
    let mut previous = 1.0;
    for rayleigh in [
        1708.0,
        1708.0_f64.next_up(),
        2000.0,
        5830.0_f64.next_down(),
        5830.0,
        5830.0_f64.next_up(),
        10_000.0,
        100_000.0,
    ] {
        let value = evaluate(CARD, inputs(rayleigh)).unwrap().evidence().value;
        assert!(value >= previous);
        previous = value;
    }
    for kink in [1708.0_f64, 5830.0] {
        let value = evaluate(CARD, inputs(kink)).unwrap().evidence().value;
        near(
            evaluate(CARD, inputs(kink.next_down()))
                .unwrap()
                .evidence()
                .value,
            value,
        );
        near(
            evaluate(CARD, inputs(kink.next_up()))
                .unwrap()
                .evidence()
                .value,
            value,
        );
    }
    let sample = inputs(10_000.0).evaluate_with_sensitivities(CARD).unwrap();
    let convert = |k, gap| {
        sample
            .heat_transfer_coefficient_with_sensitivities(
                ThermalConductivity::new(k),
                Length::new(gap),
            )
            .unwrap()
    };
    let h = convert(0.028, 0.02);
    near(convert(0.056, 0.02).value.value(), 2.0 * h.value.value());
    near(convert(0.028, 0.04).value.value(), 0.5 * h.value.value());
    near(
        h.sensitivity.d_qoi["Ra"],
        (0.028 / 0.02) * sample.evidence().sensitivity.d_qoi["Ra"],
    );
}
