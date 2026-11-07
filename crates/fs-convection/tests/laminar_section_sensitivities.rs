//! Actual finite Hausen-section tangents retain the original heated entrance.
use fs_convection::developing::evaluate_circular_duct_section;
use fs_convection::{CorrelationId, CorrelationInputs, ThermalConductivity};
use fs_qty::Length;

fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-10),
        "{actual} != {expected}"
    );
}

#[test]
fn section_transport_partials_match_local_derivatives_and_conservative_subdivision() {
    let reynolds = 1000.0;
    let prandtl = 7.0;
    let conductivity = ThermalConductivity::new(0.6);
    let diameter = Length::new(0.01);
    let sample = |re, pr, start, end| {
        evaluate_circular_duct_section(re, pr, start, end)
            .unwrap()
            .heat_transfer_coefficient_with_sensitivities(conductivity, diameter)
            .unwrap()
    };
    for (start, end) in [(0.0, 1.0), (1.0, 50.0), (50.0, 50.0 + 1e-10)] {
        let actual = sample(reynolds, prandtl, start, end);
        let delta_re = 0.01;
        let delta_pr = 1e-4;
        near(
            actual.sensitivity.d_qoi["Re"],
            (sample(reynolds + delta_re, prandtl, start, end)
                .value
                .value()
                - sample(reynolds - delta_re, prandtl, start, end)
                    .value
                    .value())
                / (2.0 * delta_re),
            3e-8,
        );
        near(
            actual.sensitivity.d_qoi["Pr"],
            (sample(reynolds, prandtl + delta_pr, start, end)
                .value
                .value()
                - sample(reynolds, prandtl - delta_pr, start, end)
                    .value
                    .value())
                / (2.0 * delta_pr),
            3e-8,
        );
        near(
            actual.sensitivity.d_qoi["k"],
            actual.value.value() / 0.6,
            2e-14,
        );
        near(
            actual.sensitivity.d_qoi["L"],
            -actual.value.value() / 0.01,
            2e-14,
        );
        assert_eq!(actual.model.discrepancy_rel.is_infinite(), start > 0.0);
    }
    let whole = sample(reynolds, prandtl, 0.0, 50.0);
    let first = sample(reynolds, prandtl, 0.0, 1.0);
    let second = sample(reynolds, prandtl, 1.0, 50.0);
    for axis in ["Re", "Pr"] {
        near(
            50.0 * whole.sensitivity.d_qoi[axis],
            first.sensitivity.d_qoi[axis] + 49.0 * second.sensitivity.d_qoi[axis],
            3e-13,
        );
    }
    let mean = CorrelationInputs::forced(reynolds, prandtl)
        .with_length_ratio(50.0)
        .evaluate_with_sensitivities(CorrelationId::CircularDuctHausen)
        .unwrap()
        .heat_transfer_coefficient_with_sensitivities(conductivity, diameter)
        .unwrap();
    near(
        whole.sensitivity.d_qoi["Re"],
        mean.sensitivity.d_qoi["Re"],
        2e-14,
    );
    assert!((second.sensitivity.d_qoi["Re"] - mean.sensitivity.d_qoi["Re"]).abs() > 1e-4);
}

#[test]
fn original_moving_domain_endpoints_remove_only_the_unavailable_partial() {
    let sample = |re, pr, start, end| {
        evaluate_circular_duct_section(re, pr, start, end)
            .unwrap()
            .heat_transfer_coefficient_with_sensitivities(
                ThermalConductivity::new(0.6),
                Length::new(0.01),
            )
            .unwrap()
    };
    for re in [1.0, 2300.0] {
        let coefficient = sample(re, 7.0, 0.05, 1000.0);
        assert!(!coefficient.sensitivity.d_qoi.contains_key("Re"));
        assert!(coefficient.sensitivity.d_qoi.contains_key("Pr"));
    }
    for pr in [0.6, 1000.0] {
        let coefficient = sample(1000.0, pr, 0.05, 1000.0);
        assert!(coefficient.sensitivity.d_qoi.contains_key("Re"));
        assert!(!coefficient.sensitivity.d_qoi.contains_key("Pr"));
    }
    // Fixed heated-length endpoints do not invalidate independent Re or Pr.
    let coefficient = sample(1000.0, 7.0, 0.05, 1000.0);
    assert!(coefficient.sensitivity.d_qoi.contains_key("Re"));
    assert!(coefficient.sensitivity.d_qoi.contains_key("Pr"));
}
