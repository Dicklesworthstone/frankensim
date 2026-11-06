//! G0/G1/G3: actual source-card admission, strip heat, and conduction lowering.

use fs_convection::developing::{FlatPlateSectionError, evaluate_flat_plate_section};
use fs_convection::{CorrelationId, CorrelationInputs, ThermalConductivity, evaluate};
use fs_qty::{Length, Temperature};

const LAMINAR: CorrelationId = CorrelationId::FlatPlateLaminarAverage;
const TURBULENT: CorrelationId = CorrelationId::FlatPlateTurbulentAverage;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs(),
        "{actual} != {expected}"
    );
}

#[test]
fn first_strip_preserves_the_original_mean_and_source_bearing_robin_row_exactly() {
    for (card, re) in [(LAMINAR, 20_000.0), (TURBULENT, 2.0e6)] {
        let whole = evaluate(card, CorrelationInputs::forced(re, 0.71)).unwrap();
        let strip = evaluate_flat_plate_section(card, re, 0.71, 0.0).unwrap();
        assert!(strip.upstream().is_none());
        assert_eq!(strip.nusselt(), whole.evidence().value);
        let k = ThermalConductivity::new(0.026);
        let end = Length::new(0.2);
        let ambient = Temperature::new(300.0);
        assert_eq!(
            strip.heat_transfer_coefficient(k, end).unwrap(),
            whole.heat_transfer_coefficient(k, end).unwrap()
        );
        assert_eq!(
            strip.robin_boundary(k, end, ambient).unwrap(),
            whole.robin_boundary(k, end, ambient).unwrap()
        );
    }
}

#[test]
fn physical_strip_heat_adds_and_both_original_source_endpoints_are_retained() {
    let conductivity = ThermalConductivity::new(0.04);
    let full_length = 0.8;
    for (card, re) in [(LAMINAR, 200_000.0), (TURBULENT, 8.0e6)] {
        let positions = [0.0, 0.125, 0.3, 0.5, 0.8, 1.0];
        let mut heat_per_width_and_kelvin = 0.0;
        for pair in positions.windows(2) {
            let [a, b] = [pair[0], pair[1]];
            let strip = evaluate_flat_plate_section(card, re * b, 0.71, a / b).unwrap();
            assert_eq!(strip.start_fraction(), a / b);
            assert_eq!(strip.downstream().card().id, card);
            if a > 0.0 {
                let upstream = strip.upstream().unwrap();
                assert_eq!(upstream.card().id, card);
                close(upstream.groups()["Re"], re * a);
                close(
                    strip.nusselt(),
                    (strip.downstream().evidence().value - upstream.evidence().value)
                        / (1.0 - a / b),
                );
            }
            let h = strip
                .heat_transfer_coefficient(conductivity, Length::new(b * full_length))
                .unwrap();
            heat_per_width_and_kelvin += h.value.value() * (b - a) * full_length;
        }
        let original = evaluate(card, CorrelationInputs::forced(re, 0.71)).unwrap();
        close(
            heat_per_width_and_kelvin,
            original.evidence().value * conductivity.value(),
        );
    }
}

#[test]
fn positive_strip_evidence_does_not_invent_a_local_validation_band() {
    let strip = evaluate_flat_plate_section(LAMINAR, 20_000.0, 0.71, 0.25).unwrap();
    let k = ThermalConductivity::new(0.026);
    let end = Length::new(0.2);
    let h = strip.heat_transfer_coefficient(k, end).unwrap();
    assert_eq!(h.model.cards, strip.downstream().evidence().model.cards);
    assert!(h.model.in_domain);
    assert_eq!(h.model.discrepancy_rel, f64::INFINITY);
    assert!(h.sensitivity.d_qoi.is_empty());
    assert!(h.adjoint_ref.is_none());
    let robin = strip
        .robin_boundary(k, end, Temperature::new(300.0))
        .unwrap();
    assert_eq!(robin.coefficient(), &h);
    assert_eq!(robin.reference_temperature(), Temperature::new(300.0));
    let changed_k = strip
        .heat_transfer_coefficient(ThermalConductivity::new(0.052), end)
        .unwrap();
    let changed_end = strip
        .heat_transfer_coefficient(k, Length::new(0.4))
        .unwrap();
    let changed_strip = evaluate_flat_plate_section(LAMINAR, 20_000.0, 0.71, 0.5)
        .unwrap()
        .heat_transfer_coefficient(k, end)
        .unwrap();
    close(changed_k.value.value(), 2.0 * h.value.value());
    close(changed_end.value.value(), 0.5 * h.value.value());
    for changed in [changed_k, changed_end, changed_strip] {
        assert_ne!(h.provenance, changed.provenance);
    }
}

#[test]
fn original_source_domains_refuse_unsupported_positive_upstream_positions_and_models() {
    assert!(matches!(
        evaluate_flat_plate_section(CorrelationId::Gnielinski, 1.0e6, 0.71, 0.5),
        Err(FlatPlateSectionError::UnsupportedCorrelation(
            CorrelationId::Gnielinski
        ))
    ));
    // An original turbulent whole plate includes its laminar entrance. A
    // positive strip start cannot silently borrow a different laminar card.
    assert!(evaluate_flat_plate_section(TURBULENT, 1.0e7, 0.71, 0.0).is_ok());
    assert!(evaluate_flat_plate_section(TURBULENT, 1.0e7, 0.71, 0.05).is_ok());
    assert!(evaluate_flat_plate_section(TURBULENT, 1.0e7, 0.71, 0.049).is_err());
    assert!(evaluate_flat_plate_section(TURBULENT, 500_000.0, 0.71, 0.5).is_err());
    assert!(evaluate_flat_plate_section(LAMINAR, 1.0, 0.71, 0.0).is_ok());
    assert!(evaluate_flat_plate_section(LAMINAR, 10.0, 0.71, 0.01).is_err());
    for fraction in [-0.1, 1.0, 1.1, f64::NAN, f64::INFINITY] {
        assert!(evaluate_flat_plate_section(LAMINAR, 10_000.0, 0.71, fraction).is_err());
    }
    for (re, pr) in [
        (0.0, 0.71),
        (500_001.0, 0.71),
        (10_000.0, 0.59),
        (10_000.0, 50.1),
        (f64::INFINITY, 0.71),
        (10_000.0, f64::NAN),
    ] {
        assert!(evaluate_flat_plate_section(LAMINAR, re, pr, 0.0).is_err());
    }
    let strip = evaluate_flat_plate_section(LAMINAR, 10_000.0, 0.71, 0.5).unwrap();
    for (k, end) in [
        (0.0, 1.0),
        (-1.0, 1.0),
        (1.0, 0.0),
        (1.0, -1.0),
        (f64::NAN, 1.0),
        (1.0, f64::INFINITY),
    ] {
        assert!(
            strip
                .heat_transfer_coefficient(ThermalConductivity::new(k), Length::new(end))
                .is_err()
        );
    }
    assert!(
        strip
            .robin_boundary(
                ThermalConductivity::new(0.026),
                Length::new(0.2),
                Temperature::new(f64::NAN)
            )
            .is_err()
    );
}

#[test]
fn analytic_section_partials_follow_the_actual_coefficient_including_thin_strips() {
    let k = ThermalConductivity::new(0.026);
    let length = Length::new(0.2);
    let epsilon = 1e-5;
    for (card, re) in [(LAMINAR, 20_000.0), (TURBULENT, 2.0e6)] {
        for fraction in [0.0, 0.5, 1.0_f64.next_down()] {
            let section = evaluate_flat_plate_section(card, re, 0.71, fraction).unwrap();
            let primal = section.heat_transfer_coefficient(k, length).unwrap();
            let tangent = section.heat_transfer_coefficient_with_sensitivities(k, length).unwrap();
            assert_eq!(tangent.value, primal.value);
            assert_eq!(tangent.model, primal.model);
            assert_eq!(tangent.numerical, primal.numerical);
            assert!(tangent.adjoint_ref.is_none());
            let h = |re, pr| evaluate_flat_plate_section(card, re, pr, fraction).unwrap()
                .heat_transfer_coefficient(k, length).unwrap().value.value();
            for (axis, input, plus, minus) in [
                ("Re", re, h(re * (1.0 + epsilon), 0.71), h(re * (1.0 - epsilon), 0.71)),
                ("Pr", 0.71, h(re, 0.71 * (1.0 + epsilon)), h(re, 0.71 * (1.0 - epsilon))),
            ] {
                let fd = (plus - minus) / (2.0 * epsilon * input);
                assert!((tangent.sensitivity.d_qoi[axis] - fd).abs() <= 3e-8 * fd.abs(),
                    "{card:?}, fraction={fraction}, axis={axis}");
            }
            close(tangent.sensitivity.d_qoi["k"], section.nusselt() / length.value());
            close(tangent.sensitivity.d_qoi["L"], -primal.value.value() / length.value());
            close(tangent.sensitivity.d_qoi["Pr"], primal.value.value() / (3.0 * 0.71));
            if fraction > 0.0 {
                let exponent = if card == LAMINAR { 0.5 } else { 0.8 };
                close(tangent.sensitivity.d_qoi["Re"], exponent * primal.value.value() / re);
                if card == TURBULENT {
                    // Reusing a downstream mean tangent would retain the 871
                    // entrance correction that cancels out of this strip.
                    assert!(section.downstream().reynolds_elasticity().unwrap().unwrap() > 0.9);
                }
                assert_eq!(tangent.model.discrepancy_rel, f64::INFINITY);
            } else {
                let mean = CorrelationInputs::forced(re, 0.71)
                    .evaluate_with_sensitivities(card).unwrap()
                    .heat_transfer_coefficient_with_sensitivities(k, length).unwrap();
                for axis in ["Re", "Pr", "k", "L"] {
                    close(tangent.sensitivity.d_qoi[axis], mean.sensitivity.d_qoi[axis]);
                }
            }
        }
    }
}

#[test]
fn section_partials_omit_only_inputs_that_move_a_source_endpoint() {
    let k = ThermalConductivity::new(0.026);
    let length = Length::new(0.2);
    for (card, re, fraction) in [
        (LAMINAR, 1.0, 0.0),
        (LAMINAR, 2.0, 0.5),
        (LAMINAR, 500_000.0, 0.5),
        (TURBULENT, 500_000.0, 0.0),
        (TURBULENT, 2.0e6, 0.25),
        (TURBULENT, 1.0e7, 0.5),
    ] {
        let tangent = evaluate_flat_plate_section(card, re, 0.71, fraction).unwrap()
            .heat_transfer_coefficient_with_sensitivities(k, length).unwrap();
        assert!(!tangent.sensitivity.d_qoi.contains_key("Re"));
        assert!(tangent.sensitivity.d_qoi["Pr"] > 0.0);
        assert!(tangent.sensitivity.d_qoi["k"] > 0.0);
        assert!(tangent.sensitivity.d_qoi["L"] < 0.0);
    }
    // A boundary in the held Pr coordinate must not erase the Re partial.
    for (card, re, high_pr) in [(LAMINAR, 20_000.0, 50.0), (TURBULENT, 2.0e6, 60.0)] {
        for pr in [0.6, high_pr] {
            let tangent = evaluate_flat_plate_section(card, re, pr, 0.5).unwrap()
                .heat_transfer_coefficient_with_sensitivities(k, length).unwrap();
            assert!(!tangent.sensitivity.d_qoi.contains_key("Pr"));
            assert!(tangent.sensitivity.d_qoi["Re"] > 0.0);
        }
    }
    let interior = evaluate_flat_plate_section(TURBULENT, 2.0e6_f64.next_up(), 0.71, 0.25)
        .unwrap().heat_transfer_coefficient_with_sensitivities(k, length).unwrap();
    assert!(interior.sensitivity.d_qoi["Re"] > 0.0);
}

#[test]
fn section_heat_derivatives_add_under_physical_subdivision() {
    let k = ThermalConductivity::new(0.04);
    let full_length = 0.8;
    let positions = [0.0, 0.125, 0.3, 0.5, 0.8, 1.0];
    for (card, re) in [(LAMINAR, 200_000.0), (TURBULENT, 8.0e6)] {
        let mut d_re = 0.0;
        let mut d_pr = 0.0;
        for pair in positions.windows(2) {
            let [a, b] = [pair[0], pair[1]];
            let h = evaluate_flat_plate_section(card, re * b, 0.71, a / b).unwrap()
                .heat_transfer_coefficient_with_sensitivities(k, Length::new(b * full_length))
                .unwrap();
            // dRe_b/dRe_full=b, while every strip shares the same Pr.
            d_re += h.sensitivity.d_qoi["Re"] * b * (b - a) * full_length;
            d_pr += h.sensitivity.d_qoi["Pr"] * (b - a) * full_length;
        }
        let whole = evaluate_flat_plate_section(card, re, 0.71, 0.0).unwrap()
            .heat_transfer_coefficient_with_sensitivities(k, Length::new(full_length)).unwrap();
        close(d_re, whole.sensitivity.d_qoi["Re"] * full_length);
        close(d_pr, whole.sensitivity.d_qoi["Pr"] * full_length);
    }
}
