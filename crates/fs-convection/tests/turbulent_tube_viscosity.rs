//! Source-formula and sensitivity checks, not experimental validation.

use fs_convection::{
    CorrelationError, CorrelationId, CorrelationInputs, DiscrepancyBasis, ThermalConductivity,
    ThermalDirection, evaluate,
};
use fs_qty::{Length, Temperature};

const CARD: CorrelationId = CorrelationId::SiederTateTurbulent;

fn inputs(re: f64, pr: f64, ratio: f64) -> CorrelationInputs {
    CorrelationInputs::forced(re, pr)
        .with_length_ratio(100.0)
        .with_viscosity_ratio(ratio)
}

fn nusselt(re: f64, pr: f64, ratio: f64) -> f64 {
    evaluate(CARD, inputs(re, pr, ratio))
        .unwrap()
        .evidence()
        .value
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual / expected - 1.0).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn source_formula_matches_independent_decimal_anchors() {
    // Independent Python Decimal arithmetic at 80 digits, with exact decimal
    // coefficients and rational exponents 4/5, 1/3 and 14/100. Source:
    // Matuska & Zmrhal, KOLEKTOR 2.2 (2009), Table 2.5, p. 28.
    // This intentionally pins the table's 0.027 prefactor rather than 0.023.
    for (re, pr, ratio, expected) in [
        (10_000.0, 0.7, 0.1, 27.525_155_268_843_22),
        (100_000.0, 1.2, 1.0, 286.917_813_679_305),
        (100_000.0, 7.0, 4.0, 627.121_239_206_940_3),
        (100_000.0, 7.0, 0.25, 425.377_693_471_240_9),
        (50_000.0, 200.0, 2.0, 999.295_577_917_829_7),
        (1_000_000.0, 16_700.0, 10.0, 60_108.504_272_297_26),
    ] {
        near(nusselt(re, pr, ratio), expected, 3e-13);
    }
    near(nusselt(100_000.0, 1.0, 1.0), 270.0, 3e-13);
    let evaluated = evaluate(CARD, inputs(100_000.0, 7.0, 4.0)).unwrap();
    assert!(
        evaluated
            .card()
            .source
            .citation
            .contains("Table 2.5, p. 28")
    );
    assert_eq!(
        evaluated.card().discrepancy_basis,
        DiscrepancyBasis::EngineeringAllowance
    );
    assert_eq!(evaluated.card().model.discrepancy_rel, 0.25);
    assert!(evaluated.card().model.assumptions.iter().any(|a| {
        a.contains("implementation admission limits, not source-published validation bounds")
    }));
}

#[test]
fn viscosity_contrast_not_the_dittus_direction_controls_the_correction() {
    let nominal = inputs(100_000.0, 7.0, 1.0);
    let reference = evaluate(CARD, nominal).unwrap();
    for (ratio, factor) in [
        (4.0, 1.214_194_884_395_047),
        (0.25, 0.823_591_017_267_573_1),
    ] {
        let corrected = evaluate(CARD, nominal.with_viscosity_ratio(ratio)).unwrap();
        near(
            corrected.evidence().value / reference.evidence().value,
            factor,
            3e-13,
        );
        assert_eq!(corrected.groups()["mu_bulk_over_surface"], ratio);
        assert_ne!(
            corrected.evidence().provenance,
            reference.evidence().provenance
        );
        let cooling = evaluate(
            CARD,
            nominal
                .with_viscosity_ratio(ratio)
                .with_direction(ThermalDirection::CoolingFluid),
        )
        .unwrap();
        assert_eq!(
            cooling.evidence().value.to_bits(),
            corrected.evidence().value.to_bits()
        );
    }
}

#[test]
fn missing_viscosity_transition_and_extrapolation_refuse_without_fallback() {
    for point in [
        CorrelationInputs::forced(100_000.0, 7.0).with_length_ratio(100.0),
        CorrelationInputs::forced(100_000.0, 7.0).with_viscosity_ratio(1.0),
    ] {
        let Err(CorrelationError::OutOfDomain {
            correlation,
            violations,
        }) = evaluate(CARD, point)
        else {
            panic!("required geometry and wall viscosity must remain explicit");
        };
        assert_eq!(correlation, CARD);
        assert!(violations.iter().any(|v| v.value.is_none()));
    }
    for point in [
        inputs(2_000.0, 7.0, 1.0),
        inputs(10_000.0_f64.next_down(), 7.0, 1.0),
        inputs(1_000_000.0_f64.next_up(), 7.0, 1.0),
        inputs(100_000.0, 0.7_f64.next_down(), 1.0),
        inputs(100_000.0, 16_700.0_f64.next_up(), 1.0),
        inputs(100_000.0, 7.0, 0.1_f64.next_down()),
        inputs(100_000.0, 7.0, 10.0_f64.next_up()),
        inputs(100_000.0, 7.0, 1.0).with_length_ratio(60.0_f64.next_down()),
        inputs(100_000.0, 7.0, 1.0).with_length_ratio(1e6_f64.next_up()),
    ] {
        assert!(matches!(
            evaluate(CARD, point),
            Err(CorrelationError::OutOfDomain {
                correlation: CARD,
                ..
            })
        ));
        assert!(point.evaluate_with_sensitivities(CARD).is_err());
    }
    for ratio in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            evaluate(CARD, inputs(100_000.0, 7.0, ratio)),
            Err(CorrelationError::InvalidGroup {
                axis: "mu_bulk_over_surface",
                ..
            })
        ));
    }
}

#[test]
fn flow_changes_transfer_but_fully_developed_length_does_not_reset_it() {
    let k = ThermalConductivity::new(0.14);
    let diameter = Length::new(0.01);
    let low = evaluate(CARD, inputs(40_000.0, 200.0, 2.0)).unwrap();
    let high = evaluate(CARD, inputs(80_000.0, 200.0, 2.0)).unwrap();
    let low_h = low.heat_transfer_coefficient(k, diameter).unwrap();
    let high_h = high.heat_transfer_coefficient(k, diameter).unwrap();
    near(
        high_h.value.value() / low_h.value.value(),
        1.741_101_126_592_248_2,
        3e-13,
    );
    assert_ne!(high.evidence().provenance, low.evidence().provenance);
    for length_ratio in [60.0, 400.0, 1e6] {
        let same_flow = evaluate(
            CARD,
            inputs(40_000.0, 200.0, 2.0).with_length_ratio(length_ratio),
        )
        .unwrap();
        let h = same_flow.heat_transfer_coefficient(k, diameter).unwrap();
        assert_eq!(h.value.value().to_bits(), low_h.value.value().to_bits());
        assert_eq!(h.model.cards, [CARD.name()]);
    }
}

#[test]
fn analytic_partials_match_perturbations_and_survive_robin_lowering() {
    for (re, pr, ratio) in [
        (20_000.0, 1.0, 0.2),
        (100_000.0, 7.0, 4.0),
        (500_000.0, 5_000.0, 8.0),
    ] {
        let point = inputs(re, pr, ratio);
        let sensitive = point.evaluate_with_sensitivities(CARD).unwrap();
        let plain = evaluate(CARD, point).unwrap();
        assert_eq!(
            sensitive.evidence().value.to_bits(),
            plain.evidence().value.to_bits()
        );
        assert_eq!(sensitive.evidence().model, plain.evidence().model);
        let step = 1e-5;
        let slopes = [
            (
                "Re",
                (nusselt(re * (1.0 + step), pr, ratio) - nusselt(re * (1.0 - step), pr, ratio))
                    / (2.0 * step * re),
            ),
            (
                "Pr",
                (nusselt(re, pr * (1.0 + step), ratio) - nusselt(re, pr * (1.0 - step), ratio))
                    / (2.0 * step * pr),
            ),
            (
                "mu_bulk_over_surface",
                (nusselt(re, pr, ratio * (1.0 + step)) - nusselt(re, pr, ratio * (1.0 - step)))
                    / (2.0 * step * ratio),
            ),
        ];
        let k = ThermalConductivity::new(0.14);
        let diameter = Length::new(0.01);
        let robin = sensitive
            .robin_boundary_with_sensitivities(k, diameter, Temperature::new(350.0))
            .unwrap();
        for (axis, derivative) in slopes {
            near(
                sensitive.evidence().sensitivity.d_qoi[axis],
                derivative,
                2e-8,
            );
            near(
                robin.coefficient().sensitivity.d_qoi[axis],
                14.0 * derivative,
                2e-8,
            );
        }
        assert_eq!(sensitive.evidence().sensitivity.d_qoi["L_over_Dh"], 0.0);
        near(
            plain.reynolds_elasticity().unwrap().unwrap(),
            slopes[0].1 * re / plain.evidence().value,
            2e-8,
        );
        assert_eq!(robin.coefficient().model.cards, [CARD.name()]);
        assert_eq!(robin.coefficient().model, plain.evidence().model);
    }
}

#[test]
fn inclusive_endpoints_omit_only_the_partial_that_leaves_the_domain() {
    for (axis, points) in [
        ("Re", [inputs(10_000.0, 7.0, 1.0), inputs(1e6, 7.0, 1.0)]),
        (
            "Pr",
            [
                inputs(100_000.0, 0.7, 1.0),
                inputs(100_000.0, 16_700.0, 1.0),
            ],
        ),
        (
            "mu_bulk_over_surface",
            [inputs(100_000.0, 7.0, 0.1), inputs(100_000.0, 7.0, 10.0)],
        ),
        (
            "L_over_Dh",
            [
                inputs(100_000.0, 7.0, 1.0).with_length_ratio(60.0),
                inputs(100_000.0, 7.0, 1.0).with_length_ratio(1e6),
            ],
        ),
    ] {
        for point in points {
            let evaluated = point.evaluate_with_sensitivities(CARD).unwrap();
            for input_axis in ["Re", "Pr", "mu_bulk_over_surface", "L_over_Dh"] {
                assert_eq!(
                    evaluated
                        .evidence()
                        .sensitivity
                        .d_qoi
                        .contains_key(input_axis),
                    input_axis != axis
                );
            }
        }
    }
}
