//! G0/G1/G3: corrected source formula, inlet history, section integration and tangents.

use fs_convection::developing::evaluate_turbulent_circular_duct_section;
use fs_convection::{
    CorrelationError, CorrelationId, CorrelationInputs, ThermalConductivity, evaluate,
};
use fs_qty::{Length, Temperature};

const CARD: CorrelationId = CorrelationId::GnielinskiDeveloping;

fn mean(re: f64, pr: f64, length: f64) -> f64 {
    evaluate(
        CARD,
        CorrelationInputs::forced(re, pr).with_length_ratio(length),
    )
    .unwrap()
    .evidence()
    .value
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn corrected_source_formula_matches_independent_eighty_digit_decimal_evaluations() {
    // Python Decimal at precision 80, evaluating Eq. (2) with K=1 and the
    // corrected Eq. (6) directly; these are formula checks, not measurements.
    for (re, pr, length, expected) in [
        (10_000.0, 0.5, 10.0, 29.700_584_552_775_114),
        (20_000.0, 0.7, 20.0, 57.165_062_570_900_855),
        (100_000.0, 7.0, 50.0, 637.891_653_724_976_5),
        (1_000_000.0, 200.0, 1_000_000.0, 16958.709_772_180_937),
    ] {
        close(mean(re, pr, length), expected, 3e-13);
    }
    // At Re=10^5 and Pr=1, z=1.8*5-1.5=7.5, so Nu_infinity=220.
    // This distinguishes the corrected -2 friction exponent and the 2013
    // coefficients from the existing, independently retained legacy card.
    close(mean(100_000.0, 1.0, 1000.0), 222.2, 3e-13);
    let legacy = evaluate(
        CorrelationId::Gnielinski,
        CorrelationInputs::forced(100_000.0, 1.0).with_length_ratio(1000.0),
    )
    .unwrap();
    assert!((legacy.evidence().value - 220.0).abs() > 2.0);
}

#[test]
fn inlet_enhancement_decays_toward_the_correct_developed_limit() {
    let mut previous = f64::INFINITY;
    for length in [10.0_f64, 27.0, 64.0, 125.0, 1000.0, 1.0e6] {
        let nu = mean(100_000.0, 1.0, length);
        assert!(nu < previous && nu > 220.0);
        close(nu / 220.0 - 1.0, 1.0 / length.cbrt().powi(2), 3e-11);
        previous = nu;
    }
    assert!(mean(100_000.0, 1.0, 10.0) > 1.2 * 220.0);
}

#[test]
fn analytic_mean_partials_match_physical_input_perturbations_and_keep_endpoints_honest() {
    for (re, pr, length) in [(20_000.0, 0.7, 20.0), (100_000.0, 7.0, 100.0)] {
        let result = CorrelationInputs::forced(re, pr)
            .with_length_ratio(length)
            .evaluate_with_sensitivities(CARD)
            .unwrap();
        assert_eq!(result.evidence().value, mean(re, pr, length));
        let epsilon = 1e-5;
        for (axis, input, plus, minus) in [
            (
                "Re",
                re,
                mean(re * (1.0 + epsilon), pr, length),
                mean(re * (1.0 - epsilon), pr, length),
            ),
            (
                "Pr",
                pr,
                mean(re, pr * (1.0 + epsilon), length),
                mean(re, pr * (1.0 - epsilon), length),
            ),
            (
                "L_over_Dh",
                length,
                mean(re, pr, length * (1.0 + epsilon)),
                mean(re, pr, length * (1.0 - epsilon)),
            ),
        ] {
            close(
                result.evidence().sensitivity.d_qoi[axis],
                (plus - minus) / (2.0 * epsilon * input),
                2e-8,
            );
        }
        assert!(result.evidence().sensitivity.d_qoi["L_over_Dh"] < 0.0);
        close(
            result.reynolds_elasticity().unwrap().unwrap(),
            result.evidence().sensitivity.d_qoi["Re"] * re / result.evidence().value,
            1e-13,
        );
    }
    let at_inlet_bound = CorrelationInputs::forced(100_000.0, 7.0)
        .with_length_ratio(10.0)
        .evaluate_with_sensitivities(CARD)
        .unwrap();
    assert!(
        !at_inlet_bound
            .evidence()
            .sensitivity
            .d_qoi
            .contains_key("L_over_Dh")
    );
    assert!(at_inlet_bound.evidence().sensitivity.d_qoi["Re"] > 0.0);
    assert!(at_inlet_bound.evidence().sensitivity.d_qoi["Pr"] > 0.0);
}

#[test]
fn section_heat_telescopes_without_restarting_the_entrance_and_lowers_to_robin() {
    let conductivity = ThermalConductivity::new(0.04);
    let diameter = Length::new(0.01);
    let positions = [0.0, 27.0, 64.0, 125.0, 1000.0];
    let mut transfer = 0.0;
    let mut restarted = 0.0;
    for pair in positions.windows(2) {
        let [start, end] = [pair[0], pair[1]];
        let section = evaluate_turbulent_circular_duct_section(100_000.0, 1.0, start, end).unwrap();
        assert_eq!(section.start_over_diameter(), start);
        assert_eq!(section.end_over_diameter(), end);
        assert_eq!(section.downstream().card().id, CARD);
        let coefficient = section
            .heat_transfer_coefficient(conductivity, diameter)
            .unwrap();
        let robin = section
            .robin_boundary(conductivity, diameter, Temperature::new(300.0))
            .unwrap();
        assert_eq!(robin.coefficient(), &coefficient);
        assert_eq!(robin.reference_temperature(), Temperature::new(300.0));
        if start == 0.0 {
            assert!(section.upstream().is_none());
            assert_eq!(section.nusselt(), mean(100_000.0, 1.0, end));
            assert_eq!(
                coefficient,
                section
                    .downstream()
                    .heat_transfer_coefficient(conductivity, diameter)
                    .unwrap()
            );
        } else {
            let upstream = section.upstream().unwrap();
            assert_eq!(upstream.groups()["L_over_Dh"], start);
            close(
                section.nusselt(),
                (end * mean(100_000.0, 1.0, end) - start * mean(100_000.0, 1.0, start))
                    / (end - start),
                3e-13,
            );
            assert_eq!(coefficient.model.discrepancy_rel, f64::INFINITY);
            assert!(coefficient.sensitivity.d_qoi.is_empty());
        }
        transfer += coefficient.value.value() * (end - start) * diameter.value();
        restarted += mean(100_000.0, 1.0, end - start) * (end - start) * conductivity.value();
    }
    close(transfer, 222.2 * 1000.0 * conductivity.value(), 3e-13);
    assert!(
        restarted > transfer * 1.008,
        "a per-cell entrance reset must spuriously add heat"
    );
    let middle = evaluate_turbulent_circular_duct_section(100_000.0, 1.0, 27.0, 64.0).unwrap();
    close(middle.nusselt(), 220.0 * 38.0 / 37.0, 3e-13);
    let (start, end, upstream, downstream, nu) = middle.into_parts();
    assert_eq!((start, end), (27.0, 64.0));
    assert_eq!(upstream.unwrap().groups()["L_over_Dh"], start);
    assert_eq!(downstream.groups()["L_over_Dh"], end);
    close(nu, 220.0 * 38.0 / 37.0, 3e-13);
    assert_ne!(nu, downstream.evidence().value);
}

#[test]
fn adjacent_endpoints_retain_the_finite_local_limit() {
    for start in [27.0_f64, 64.0, 1000.0, 999_999.0] {
        let end = start.next_up();
        let section = evaluate_turbulent_circular_duct_section(100_000.0, 1.0, start, end).unwrap();
        let expected = 220.0 * (1.0 + 1.0 / (3.0 * start.cbrt().powi(2)));
        close(section.nusselt(), expected, 3e-13);
    }
    let a = 27.0_f64;
    let b = a.next_up();
    let naive = (b * mean(100_000.0, 1.0, b) - a * mean(100_000.0, 1.0, a)) / (b - a);
    assert!(
        (naive - 220.0 * (1.0 + 1.0 / 27.0)).abs() > 1.0,
        "fixture must expose cancellation in the raw cumulative difference"
    );
}

#[test]
fn source_admission_and_positive_section_history_are_not_silently_extended() {
    assert!(matches!(
        evaluate(CARD, CorrelationInputs::forced(100_000.0, 7.0)),
        Err(CorrelationError::OutOfDomain { .. })
    ));
    for (re, pr, length) in [
        (9999.0, 7.0, 50.0),
        (1_000_001.0, 7.0, 50.0),
        (100_000.0, 0.49, 50.0),
        (100_000.0, 201.0, 50.0),
        (100_000.0, 7.0, 9.99),
        (100_000.0, 7.0, 1_000_001.0),
    ] {
        assert!(matches!(
            evaluate(
                CARD,
                CorrelationInputs::forced(re, pr).with_length_ratio(length)
            ),
            Err(CorrelationError::OutOfDomain { .. })
        ));
    }
    for (a, b) in [
        (-1.0, 20.0),
        (1.0, 20.0),
        (9.99, 20.0),
        (20.0, 20.0),
        (20.0, 19.0),
        (f64::NAN, 20.0),
        (20.0, f64::INFINITY),
        (0.0, 9.99),
    ] {
        assert!(evaluate_turbulent_circular_duct_section(100_000.0, 7.0, a, b).is_err());
    }
    let section = evaluate_turbulent_circular_duct_section(100_000.0, 7.0, 20.0, 30.0).unwrap();
    assert!(
        section
            .heat_transfer_coefficient(ThermalConductivity::new(0.0), Length::new(0.01))
            .is_err()
    );
    assert!(
        section
            .heat_transfer_coefficient(ThermalConductivity::new(0.04), Length::new(0.0))
            .is_err()
    );
    assert!(
        section
            .robin_boundary(
                ThermalConductivity::new(0.04),
                Length::new(0.01),
                Temperature::new(f64::NAN)
            )
            .is_err()
    );
}

#[test]
fn section_partials_use_the_integrated_coefficient_including_adjacent_endpoints() {
    let k = ThermalConductivity::new(0.04);
    let diameter = Length::new(0.01);
    let epsilon = 1e-5;
    for (re, pr, start, end) in [
        (20_000.0, 0.7, 0.0, 20.0),
        (100_000.0, 7.0, 27.0, 64.0),
        (300_000.0, 80.0, 1000.0, 1000.0_f64.next_up()),
    ] {
        let section = evaluate_turbulent_circular_duct_section(re, pr, start, end).unwrap();
        let primal = section.heat_transfer_coefficient(k, diameter).unwrap();
        let tangent = section.heat_transfer_coefficient_with_sensitivities(k, diameter).unwrap();
        assert_eq!(tangent.value, primal.value);
        assert_eq!(tangent.model, primal.model);
        assert_eq!(tangent.numerical, primal.numerical);
        assert!(tangent.adjoint_ref.is_none());
        let h = |re, pr| evaluate_turbulent_circular_duct_section(re, pr, start, end).unwrap()
            .heat_transfer_coefficient(k, diameter).unwrap().value.value();
        for (axis, input, plus, minus) in [
            ("Re", re, h(re * (1.0 + epsilon), pr), h(re * (1.0 - epsilon), pr)),
            ("Pr", pr, h(re, pr * (1.0 + epsilon)), h(re, pr * (1.0 - epsilon))),
        ] {
            close(tangent.sensitivity.d_qoi[axis], (plus - minus) / (2.0 * epsilon * input), 3e-8);
        }
        close(tangent.sensitivity.d_qoi["k"], section.nusselt() / diameter.value(), 3e-13);
        close(tangent.sensitivity.d_qoi["L"], -primal.value.value() / diameter.value(), 3e-13);
        let downstream = CorrelationInputs::forced(re, pr).with_length_ratio(end)
            .evaluate_with_sensitivities(CARD).unwrap()
            .heat_transfer_coefficient_with_sensitivities(k, diameter).unwrap();
        for axis in ["Re", "Pr"] {
            close(tangent.sensitivity.d_qoi[axis] / primal.value.value(),
                downstream.sensitivity.d_qoi[axis] / downstream.value.value(), 3e-13);
            if start > 0.0 {
                assert!(tangent.sensitivity.d_qoi[axis] < downstream.sensitivity.d_qoi[axis]);
            }
        }
        if start > 0.0 {
            assert_eq!(tangent.model.discrepancy_rel, f64::INFINITY);
        }
    }
}

#[test]
fn section_partials_add_over_a_shared_inlet_and_preserve_boundary_omissions() {
    let k = ThermalConductivity::new(0.04);
    let diameter = Length::new(0.01);
    let positions = [0.0, 27.0, 64.0, 125.0, 1000.0];
    let mut d_re = 0.0;
    let mut d_pr = 0.0;
    for pair in positions.windows(2) {
        let [a, b] = [pair[0], pair[1]];
        let h = evaluate_turbulent_circular_duct_section(100_000.0, 7.0, a, b).unwrap()
            .heat_transfer_coefficient_with_sensitivities(k, diameter).unwrap();
        d_re += h.sensitivity.d_qoi["Re"] * (b - a) * diameter.value();
        d_pr += h.sensitivity.d_qoi["Pr"] * (b - a) * diameter.value();
    }
    let whole = evaluate_turbulent_circular_duct_section(100_000.0, 7.0, 0.0, 1000.0).unwrap()
        .heat_transfer_coefficient_with_sensitivities(k, diameter).unwrap();
    close(d_re, whole.sensitivity.d_qoi["Re"] * 1000.0 * diameter.value(), 3e-13);
    close(d_pr, whole.sensitivity.d_qoi["Pr"] * 1000.0 * diameter.value(), 3e-13);
    // Holding either source length boundary fixed does not move Re or Pr.
    let length_bounds = evaluate_turbulent_circular_duct_section(100_000.0, 7.0, 10.0, 1.0e6)
        .unwrap().heat_transfer_coefficient_with_sensitivities(k, diameter).unwrap();
    assert!(length_bounds.sensitivity.d_qoi["Re"] > 0.0);
    assert!(length_bounds.sensitivity.d_qoi["Pr"] > 0.0);
    for re in [10_000.0, 1.0e6] {
        let h = evaluate_turbulent_circular_duct_section(re, 7.0, 27.0, 64.0).unwrap()
            .heat_transfer_coefficient_with_sensitivities(k, diameter).unwrap();
        assert!(!h.sensitivity.d_qoi.contains_key("Re"));
        assert!(h.sensitivity.d_qoi["Pr"] > 0.0);
    }
    for pr in [0.5, 200.0] {
        let h = evaluate_turbulent_circular_duct_section(100_000.0, pr, 27.0, 64.0).unwrap()
            .heat_transfer_coefficient_with_sensitivities(k, diameter).unwrap();
        assert!(!h.sensitivity.d_qoi.contains_key("Pr"));
        assert!(h.sensitivity.d_qoi["Re"] > 0.0);
    }
}
