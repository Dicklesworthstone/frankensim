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
