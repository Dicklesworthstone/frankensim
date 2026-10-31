//! G0/G3 isolated-sphere formula, conductive limit, domain and derivative checks.
//! Formula arithmetic and physical limits are not experimental validation.

use fs_convection::{
    CorrelationError, CorrelationId, CorrelationInputs, NaturalConvectionGeometry,
    ThermalConductivity, evaluate,
};
use fs_qty::Length;

const CARD: CorrelationId = CorrelationId::ChurchillSphere;

fn nusselt(ra: f64, pr: f64) -> f64 {
    evaluate(CARD, CorrelationInputs::natural(ra, pr))
        .expect("admitted isolated sphere")
        .evidence()
        .value
}

fn relative_close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual / expected - 1.0).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

#[test]
fn source_equation_matches_independent_high_precision_arithmetic() {
    // Test-owned 70-decimal-digit Python Decimal evaluation of COMSOL 6.4
    // Heat Transfer Theory, Eq. (4-194), with decimal source coefficients.
    // These are arithmetic references, not measurements or fitted outputs.
    for (ra, pr, expected) in [
        (1.0e-8, 0.7, 2.004_537_775_894_906_5),
        (1.0e4, 0.7, 6.537_775_894_906_464),
        (1.0e6, 7.0, 19.058_994_303_879_462),
        (1.0e11, 100.0, 326.253_814_441_648_8),
    ] {
        relative_close(nusselt(ra, pr), expected, 2.0e-13);
    }
    let value = evaluate(CARD, CorrelationInputs::natural(1.0e6, 7.0)).unwrap();
    assert_eq!(value.card().model.validity.bounds()["Ra"], (0.0, 1.0e11));
    assert_eq!(value.card().model.validity.bounds()["Pr"], (0.7, f64::MAX));
    assert!(value.card().source.citation.contains("(4-194)"));
}

#[test]
fn stagnant_sphere_recovers_unbounded_radial_conduction_with_diameter_length() {
    let geometry = NaturalConvectionGeometry::Sphere;
    let stagnant = geometry.evaluate(0.0, 7.0).unwrap();
    assert_eq!(stagnant.card().id, CARD);
    assert_eq!(stagnant.evidence().value, 2.0);
    // Solving steady spherical conduction gives G=4*pi*k*r. The Nu=2 term
    // already supplies that entire conductance; no extra conduction link is
    // appropriate. Test several lengths to distinguish diameter from radius.
    let conductivity = 0.6;
    for diameter in [0.002, 0.02, 0.2] {
        let h = stagnant
            .heat_transfer_coefficient(
                ThermalConductivity::new(conductivity),
                Length::new(diameter),
            )
            .unwrap()
            .value
            .value();
        let area = core::f64::consts::PI * diameter * diameter;
        let radial_conductance = 4.0 * core::f64::consts::PI * conductivity * (diameter / 2.0);
        relative_close(h * area, radial_conductance, 2.0e-15);
    }
    // An isolated sphere has no face orientation. Reversing signed buoyancy
    // preserves the average magnitude under this constant-property model.
    assert_eq!(
        geometry.evaluate(1.0e6, 7.0),
        geometry.evaluate(-1.0e6, 7.0)
    );
    for ra in [1.0e-8, 1.0, 1.0e6] {
        let excess = nusselt(ra, 7.0) - 2.0;
        relative_close(nusselt(16.0 * ra, 7.0) - 2.0, 2.0 * excess, 2.0e-12);
    }
}

#[test]
fn published_domain_and_sphere_only_zero_rayleigh_admission_are_enforced() {
    for ra in [0.0, f64::from_bits(1), 1.0e11] {
        for pr in [0.7, 7.0, f64::MAX] {
            assert!(nusselt(ra, pr).is_finite());
        }
    }
    for inputs in [
        CorrelationInputs::natural(1.0e11_f64.next_up(), 0.7),
        CorrelationInputs::natural(1.0e6, 0.7_f64.next_down()),
        CorrelationInputs::forced(1.0e3, 0.7),
    ] {
        assert!(matches!(
            evaluate(CARD, inputs),
            Err(CorrelationError::OutOfDomain { .. })
        ));
    }
    for (ra, pr) in [
        (-1.0, 0.7),
        (f64::NAN, 0.7),
        (1.0, f64::INFINITY),
        (1.0, 0.0),
    ] {
        assert!(matches!(
            evaluate(CARD, CorrelationInputs::natural(ra, pr)),
            Err(CorrelationError::InvalidGroup { .. })
        ));
    }
    for id in [
        CorrelationId::ChurchillChuVerticalPlate,
        CorrelationId::ChurchillChuHorizontalCylinder,
    ] {
        assert!(matches!(
            evaluate(id, CorrelationInputs::natural(0.0, 0.7)),
            Err(CorrelationError::InvalidGroup { axis: "Ra", .. })
        ));
    }
    // The sphere's stagnant limit does not admit zero unrelated groups.
    assert!(matches!(
        evaluate(CARD, CorrelationInputs::mixed(0.0, 1.0e6, 0.7)),
        Err(CorrelationError::InvalidGroup { axis: "Re", .. })
    ));
}

#[test]
fn analytic_partials_match_independent_perturbations_without_changing_primal_evidence() {
    for ra in [1.0, 1.0e6, 1.0e10] {
        for pr in [1.0, 7.0, 100.0] {
            let sensitive = CorrelationInputs::natural(ra, pr)
                .evaluate_with_sensitivities(CARD)
                .unwrap();
            let plain = evaluate(CARD, CorrelationInputs::natural(ra, pr)).unwrap();
            assert_eq!(
                sensitive.evidence().value.to_bits(),
                plain.evidence().value.to_bits()
            );
            assert_eq!(sensitive.evidence().model, plain.evidence().model);
            let step = 1.0e-5;
            let dra = (nusselt(ra * (1.0 + step), pr) - nusselt(ra * (1.0 - step), pr))
                / (2.0 * step * ra);
            let dpr = (nusselt(ra, pr * (1.0 + step)) - nusselt(ra, pr * (1.0 - step)))
                / (2.0 * step * pr);
            relative_close(sensitive.evidence().sensitivity.d_qoi["Ra"], dra, 2.0e-8);
            relative_close(sensitive.evidence().sensitivity.d_qoi["Pr"], dpr, 2.0e-8);
        }
    }
}

#[test]
fn conductive_endpoint_omits_singular_partial_and_tiny_buoyancy_keeps_its_derivative() {
    let stagnant = CorrelationInputs::natural(0.0, 7.0)
        .evaluate_with_sensitivities(CARD)
        .unwrap();
    assert!(!stagnant.evidence().sensitivity.d_qoi.contains_key("Ra"));
    assert_eq!(stagnant.evidence().sensitivity.d_qoi["Pr"], 0.0);

    let tiny = CorrelationInputs::natural(1.0e-80, 7.0)
        .evaluate_with_sensitivities(CARD)
        .unwrap();
    assert_eq!(tiny.evidence().value, 2.0);
    // Independent 70-digit derivatives remain finite when Nu's buoyant excess
    // is smaller than one ulp of the conduction floor.
    relative_close(
        tiny.evidence().sensitivity.d_qoi["Ra"],
        1.348_631_914_802_441_4e59,
        2.0e-13,
    );
    relative_close(
        tiny.evidence().sensitivity.d_qoi["Pr"],
        3.456_195_182_470_212e-23,
        2.0e-13,
    );

    let upper = CorrelationInputs::natural(1.0e11, 7.0)
        .evaluate_with_sensitivities(CARD)
        .unwrap();
    assert!(!upper.evidence().sensitivity.d_qoi.contains_key("Ra"));
    assert!(upper.evidence().sensitivity.d_qoi["Pr"] > 0.0);
    let lower_pr = CorrelationInputs::natural(1.0e6, 0.7)
        .evaluate_with_sensitivities(CARD)
        .unwrap();
    assert!(!lower_pr.evidence().sensitivity.d_qoi.contains_key("Pr"));
    assert!(lower_pr.evidence().sensitivity.d_qoi["Ra"] > 0.0);
}
