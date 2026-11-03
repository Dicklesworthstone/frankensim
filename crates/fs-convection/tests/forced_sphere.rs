//! G0/G3 checks of Whitaker's source formula and its declared bulk-property domain.
//! Independent arithmetic checks do not constitute experimental validation.

use fs_convection::{
    CorrelationError, CorrelationId, CorrelationInputs, ThermalConductivity, evaluate,
};
use fs_qty::Length;

const CARD: CorrelationId = CorrelationId::WhitakerSphere;

fn inputs(re: f64, pr: f64, ratio: f64) -> CorrelationInputs {
    CorrelationInputs::forced(re, pr).with_viscosity_ratio(ratio)
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
fn source_equation_matches_independent_high_precision_arithmetic() {
    // Test-owned 70-decimal-digit Python Decimal transcription of Whitaker
    // (1972), Eq. (33), DOI 10.1002/aic.690180219. Decimal source coefficients
    // and rational exponents are evaluated independently of fs-math::det.
    for (re, pr, ratio, expected) in [
        (3.5, 0.71, 1.0, 2.773_130_659_789_787_3),
        (100.0, 1.0, 1.0, 7.292_660_814_019_13),
        (1_000.0, 7.0, 2.0, 50.300_857_440_281_03),
        (76_000.0, 380.0, 3.2, 3_138.961_444_107_458),
    ] {
        near(nusselt(re, pr, ratio), expected, 3e-13);
    }
    let evaluated = evaluate(CARD, inputs(1_000.0, 7.0, 2.0)).unwrap();
    assert_eq!(
        evaluated.card().source.identifier,
        "doi:10.1002/aic.690180219"
    );
    assert!(
        evaluated
            .card()
            .source
            .citation
            .contains("Eq. (33), Table 1")
    );
    assert_eq!(
        evaluated.card().model.validity.bounds()["Re"],
        (3.5, 76_000.0)
    );
    assert_eq!(
        evaluated.card().model.validity.bounds()["Pr"],
        (0.71, 380.0)
    );
    assert_eq!(
        evaluated.card().model.validity.bounds()["mu_bulk_over_surface"],
        (1.0, 3.2)
    );
}

#[test]
fn original_domain_refuses_stagnation_extrapolation_and_missing_viscosity() {
    let missing = evaluate(CARD, CorrelationInputs::forced(100.0, 7.0));
    let Err(CorrelationError::OutOfDomain { violations, .. }) = missing else {
        panic!("surface viscosity must not be silently replaced by bulk viscosity");
    };
    assert!(
        violations
            .iter()
            .any(|v| v.axis == "mu_bulk_over_surface" && v.value.is_none())
    );
    for point in [
        inputs(3.5_f64.next_down(), 7.0, 2.0),
        inputs(76_000.0_f64.next_up(), 7.0, 2.0),
        inputs(100.0, 0.71_f64.next_down(), 2.0),
        inputs(100.0, 380.0_f64.next_up(), 2.0),
        inputs(100.0, 7.0, 1.0_f64.next_down()),
        inputs(100.0, 7.0, 3.2_f64.next_up()),
    ] {
        assert!(matches!(
            evaluate(CARD, point),
            Err(CorrelationError::OutOfDomain { .. })
        ));
    }
    for point in [
        inputs(0.0, 7.0, 2.0),
        inputs(f64::NAN, 7.0, 2.0),
        inputs(100.0, 7.0, 0.0),
        inputs(100.0, 7.0, -1.0),
        inputs(100.0, 7.0, f64::NAN),
        inputs(100.0, 7.0, f64::INFINITY),
    ] {
        assert!(matches!(
            evaluate(CARD, point),
            Err(CorrelationError::InvalidGroup { .. })
        ));
    }
}

#[test]
fn viscosity_correction_changes_only_the_forced_excess_and_retains_provenance() {
    let reference = evaluate(CARD, inputs(100.0, 7.0, 1.0)).unwrap();
    let corrected = evaluate(CARD, inputs(100.0, 7.0, 2.0)).unwrap();
    near(
        corrected.evidence().value - 2.0,
        (reference.evidence().value - 2.0) * 2.0_f64.sqrt().sqrt(),
        3e-14,
    );
    assert_ne!(
        reference.evidence().provenance,
        corrected.evidence().provenance
    );
    assert_eq!(corrected.groups()["mu_bulk_over_surface"], 2.0);
    for diameter in [0.002, 0.02, 0.2] {
        let coefficient = corrected
            .heat_transfer_coefficient(ThermalConductivity::new(0.6), Length::new(diameter))
            .unwrap();
        let area = core::f64::consts::PI * diameter * diameter;
        let radial = 2.0 * core::f64::consts::PI * 0.6 * diameter;
        // The conductive contribution is already in Nu; count it once.
        near(
            coefficient.value.value() * area,
            radial * (corrected.evidence().value / 2.0),
            3e-15,
        );
        assert_eq!(coefficient.model.cards, [CARD.name()]);
    }
}

#[test]
fn analytic_partials_and_dimensional_chain_rule_match_independent_perturbations() {
    for (re, pr, ratio) in [
        (10.0, 0.8, 1.1),
        (1_000.0, 7.0, 2.0),
        (50_000.0, 300.0, 3.0),
    ] {
        let sensitive = inputs(re, pr, ratio)
            .evaluate_with_sensitivities(CARD)
            .unwrap();
        let plain = evaluate(CARD, inputs(re, pr, ratio)).unwrap();
        assert_eq!(
            sensitive.evidence().value.to_bits(),
            plain.evidence().value.to_bits()
        );
        assert_eq!(sensitive.evidence().model, plain.evidence().model);
        let step = 1e-5;
        let d_re = (nusselt(re * (1.0 + step), pr, ratio) - nusselt(re * (1.0 - step), pr, ratio))
            / (2.0 * step * re);
        let d_pr = (nusselt(re, pr * (1.0 + step), ratio) - nusselt(re, pr * (1.0 - step), ratio))
            / (2.0 * step * pr);
        let d_ratio = (nusselt(re, pr, ratio * (1.0 + step))
            - nusselt(re, pr, ratio * (1.0 - step)))
            / (2.0 * step * ratio);
        for (axis, derivative) in [
            ("Re", d_re),
            ("Pr", d_pr),
            ("mu_bulk_over_surface", d_ratio),
        ] {
            near(
                sensitive.evidence().sensitivity.d_qoi[axis],
                derivative,
                2e-8,
            );
        }
        near(
            plain.reynolds_elasticity().unwrap().unwrap(),
            d_re * re / plain.evidence().value,
            2e-8,
        );
        let h = sensitive
            .heat_transfer_coefficient_with_sensitivities(
                ThermalConductivity::new(0.6),
                Length::new(0.01),
            )
            .unwrap();
        near(
            h.sensitivity.d_qoi["mu_bulk_over_surface"],
            60.0 * d_ratio,
            2e-8,
        );
    }
}

#[test]
fn a_fixed_domain_endpoint_does_not_erase_other_input_partials() {
    for ratio in [1.0, 3.2] {
        let point = inputs(100.0, 7.0, ratio)
            .evaluate_with_sensitivities(CARD)
            .unwrap();
        let partials = &point.evidence().sensitivity.d_qoi;
        assert!(partials.contains_key("Re"));
        assert!(partials.contains_key("Pr"));
        assert!(!partials.contains_key("mu_bulk_over_surface"));
    }
    for re in [3.5, 76_000.0] {
        let point = inputs(re, 7.0, 2.0)
            .evaluate_with_sensitivities(CARD)
            .unwrap();
        let partials = &point.evidence().sensitivity.d_qoi;
        assert!(!partials.contains_key("Re"));
        assert!(partials.contains_key("Pr"));
        assert!(partials.contains_key("mu_bulk_over_surface"));
    }
}
