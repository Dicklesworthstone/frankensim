//! Physical-section admission, conservation and original conduction lowering.

use fs_convection::developing::evaluate_rectangular_duct_section;
use fs_convection::{CorrelationId, CorrelationInputs, ThermalConductivity, evaluate};
use fs_qty::{Length, Temperature};

#[test]
fn first_section_preserves_original_coefficient_and_later_sections_preserve_total_heat() {
    let original = evaluate(
        CorrelationId::RectangularDuctLaminarCwtDevelopingPr072,
        CorrelationInputs::forced(100.0, 0.72)
            .with_aspect_ratio(0.5)
            .with_length_ratio(12.0),
    )
    .unwrap();
    let first = evaluate_rectangular_duct_section(100.0, 0.72, 0.5, 0.0, 12.0).unwrap();
    let k = ThermalConductivity::new(0.03);
    let diameter = Length::new(0.01);
    let ambient = Temperature::new(300.0);
    assert!(first.upstream().is_none());
    assert_eq!(
        first.heat_transfer_coefficient(k, diameter).unwrap(),
        original.heat_transfer_coefficient(k, diameter).unwrap()
    );
    assert_eq!(
        first.robin_boundary(k, diameter, ambient).unwrap(),
        original.robin_boundary(k, diameter, ambient).unwrap()
    );

    let mut integral = 0.0;
    for pair in [0.0, 0.4, 0.63, 1.4, 4.3, 7.8, 12.0].windows(2) {
        let section =
            evaluate_rectangular_duct_section(100.0, 0.72, 0.5, pair[0], pair[1]).unwrap();
        let h = section.heat_transfer_coefficient(k, diameter).unwrap();
        integral += h.qoi * (pair[1] - pair[0]) * diameter.value();
        assert_eq!(
            section.downstream().card().id,
            CorrelationId::RectangularDuctLaminarCwtDevelopingPr072
        );
        if pair[0] > 0.0 {
            assert_eq!(section.upstream().unwrap().groups()["Pr"], 0.72);
            assert_eq!(h.model.discrepancy_rel, f64::INFINITY);
            assert!(h.sensitivity.d_qoi.is_empty());
        }
        assert_eq!(
            section
                .robin_boundary(k, diameter, ambient)
                .unwrap()
                .coefficient(),
            &h
        );
    }
    let expected = 12.0 * k.value() * original.evidence().value;
    assert!((integral - expected).abs() < 3e-12 * expected);
}

#[test]
fn narrow_source_support_and_both_positive_endpoints_remain_authoritative() {
    for (re, pr, aspect, start, end) in [
        (0.0, 0.72, 0.5, 0.0, 1.0),
        (2301.0, 0.72, 0.5, 0.0, 20.0),
        (100.0, 0.71, 0.5, 0.0, 1.0),
        (100.0, 0.72, 0.4, 0.0, 1.0),
        (100.0, 0.72, 0.5, 0.1, 1.0), // upstream Gz=720
        (100.0, 0.72, 0.5, 0.0, 0.1),
        (100.0, 0.72, 0.5, 0.0, 7201.0),
        (100.0, 0.72, 0.5, -1.0, 1.0),
        (100.0, 0.72, 0.5, 1.0, 1.0),
        (100.0, 0.72, 0.5, 1.0, f64::INFINITY),
    ] {
        assert!(evaluate_rectangular_duct_section(re, pr, aspect, start, end).is_err());
    }
    let section = evaluate_rectangular_duct_section(100.0, 0.72, 0.5, 5.0, 10.0).unwrap();
    assert!(
        section
            .heat_transfer_coefficient(ThermalConductivity::new(0.03), Length::new(0.0))
            .is_err()
    );
    assert!(
        section
            .robin_boundary(
                ThermalConductivity::new(0.03),
                Length::new(0.01),
                Temperature::new(f64::NAN)
            )
            .is_err()
    );
}
