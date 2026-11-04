//! G0/G3 checks of the bounded inclined air-gap card and its evidence lowering.
//! Formula checks establish implementation accuracy, not experimental validation.

use fs_convection::{
    CorrelationError, CorrelationId, CorrelationInputs, DiscrepancyBasis, EvaluationSourceRegion,
    ThermalConductivity, evaluate,
};
use fs_qty::Length;

const CARD: CorrelationId = CorrelationId::HollandsInclinedAirLayer;

fn inputs(rayleigh: f64, degrees: f64) -> CorrelationInputs {
    CorrelationInputs::inclined_air_layer(rayleigh, 0.71, 0.05, degrees)
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual / expected - 1.0).abs() < 3e-13,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn explicit_angle_and_source_authority_survive_typed_conductance_lowering() {
    let sample = evaluate(CARD, inputs(10_000.0, 45.0)).unwrap();
    // Independent 80-digit Decimal evaluation of the printed EnergyPlus
    // equation, including the additional sin(1.8 theta)^1.6 correction.
    near(sample.evidence().value, 1.899_983_267_836_469_4);
    assert_eq!(sample.groups()["Ra"], 10_000.0);
    assert_eq!(sample.groups()["inclination_degrees"], 45.0);
    assert!(
        sample
            .card()
            .source
            .identifier
            .contains("doi:10.1115/1.3450517")
    );
    assert!(
        sample
            .card()
            .source
            .identifier
            .contains("solar-collectors.html")
    );
    assert!(sample.card().source.citation.contains("EnergyPlus"));
    assert_eq!(sample.source_region(), EvaluationSourceRegion::CitedFormula);
    assert_eq!(
        sample.card().discrepancy_basis,
        DiscrepancyBasis::EngineeringAllowance
    );
    let h = sample
        .heat_transfer_coefficient(ThermalConductivity::new(0.028), Length::new(0.02))
        .unwrap();
    near(h.value.value(), 1.899_983_267_836_469_4 * 0.028 / 0.02);
    assert_ne!(
        sample.evidence().provenance,
        evaluate(CARD, inputs(10_000.0, 30.0))
            .unwrap()
            .evidence()
            .provenance
    );
    assert!(
        sample
            .card()
            .model
            .assumptions
            .iter()
            .any(|line| line.contains("conservative implementation subset"))
    );
}

#[test]
fn horizontal_limit_is_bitwise_compatible_and_subcritical_tilt_keeps_total_conduction() {
    for rayleigh in [0.0, 1e-300, 1708.0, 2000.0, 5830.0, 10_000.0, 100_000.0] {
        let old = evaluate(
            CorrelationId::HollandsHorizontalAirLayer,
            CorrelationInputs::horizontal_air_layer(rayleigh, 0.71, 0.05),
        )
        .unwrap();
        assert_eq!(
            evaluate(CARD, inputs(rayleigh, 0.0))
                .unwrap()
                .evidence()
                .value
                .to_bits(),
            old.evidence().value.to_bits()
        );
    }
    for (rayleigh, degrees) in [(0.0, 45.0), (1000.0, 30.0), (3000.0, 60.0)] {
        let sample = evaluate(CARD, inputs(rayleigh, degrees)).unwrap();
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
    // Stable horizontal stratification retains its separate, signed authority.
    assert_eq!(
        evaluate(
            CorrelationId::HollandsHorizontalAirLayer,
            CorrelationInputs::horizontal_air_layer(-10_000.0, 0.71, 0.05),
        )
        .unwrap()
        .evidence()
        .value,
        1.0
    );
}

#[test]
fn missing_angle_heated_above_and_unsupported_inclined_domains_refuse() {
    let Err(CorrelationError::OutOfDomain { violations, .. }) = evaluate(
        CARD,
        CorrelationInputs::horizontal_air_layer(10_000.0, 0.71, 0.05),
    ) else {
        panic!("inclination cannot be inferred from missing geometry");
    };
    assert!(
        violations
            .iter()
            .any(|v| v.axis == "inclination_degrees" && v.value.is_none())
    );
    for (rayleigh, prandtl, aspect, degrees) in [
        (100_000.0_f64.next_up(), 0.71, 0.05, 45.0),
        (10_000.0, 0.71, 0.05, 60.0_f64.next_up()),
        (10_000.0, 0.70_f64.next_down(), 0.05, 30.0),
        (10_000.0, 0.72_f64.next_up(), 0.05, 30.0),
        (10_000.0, 0.71, 0.1, 30.0),
    ] {
        assert!(matches!(
            evaluate(
                CARD,
                CorrelationInputs::inclined_air_layer(rayleigh, prandtl, aspect, degrees)
            ),
            Err(CorrelationError::OutOfDomain { .. })
        ));
    }
    for (rayleigh, degrees) in [
        (-1.0, 45.0),
        (10_000.0, -1.0),
        (f64::NAN, 30.0),
        (10_000.0, f64::NAN),
        (10_000.0, f64::INFINITY),
    ] {
        assert!(matches!(
            evaluate(CARD, inputs(rayleigh, degrees)),
            Err(CorrelationError::InvalidGroup { .. })
        ));
    }
    for degrees in [0.0, 60.0] {
        assert!(evaluate(CARD, inputs(100_000.0, degrees)).is_ok());
    }
}

#[test]
fn analytic_partials_preserve_degree_units_domain_edges_and_physical_kinks() {
    let sample = inputs(10_000.0, 45.0)
        .evaluate_with_sensitivities(CARD)
        .unwrap();
    let partials = &sample.evidence().sensitivity.d_qoi;
    near(partials["Ra"], 0.000_087_957_625_570_388_23);
    near(partials["inclination_degrees"], -0.017_410_556_438_046_24);
    assert_eq!(partials["Pr"], 0.0);
    assert_eq!(partials["aspect_ratio"], 0.0);
    let h = sample
        .heat_transfer_coefficient_with_sensitivities(
            ThermalConductivity::new(0.028),
            Length::new(0.02),
        )
        .unwrap();
    near(
        h.sensitivity.d_qoi["inclination_degrees"],
        (0.028 / 0.02) * partials["inclination_degrees"],
    );
    for (rayleigh, degrees, absent) in [
        (10_000.0, 0.0, "inclination_degrees"),
        (10_000.0, 60.0, "inclination_degrees"),
        (0.0, 45.0, "Ra"),
        (100_000.0, 45.0, "Ra"),
        (1708.0, 0.0, "Ra"),
        (5830.0, 0.0, "Ra"),
    ] {
        let sample = inputs(rayleigh, degrees)
            .evaluate_with_sensitivities(CARD)
            .unwrap();
        assert!(!sample.evidence().sensitivity.d_qoi.contains_key(absent));
        assert_eq!(sample.evidence().sensitivity.d_qoi["Pr"], 0.0);
        assert_eq!(sample.evidence().sensitivity.d_qoi["aspect_ratio"], 0.0);
    }
    let rest = inputs(1000.0, 30.0)
        .evaluate_with_sensitivities(CARD)
        .unwrap();
    assert_eq!(rest.evidence().sensitivity.d_qoi["Ra"], 0.0);
    assert_eq!(
        rest.evidence().sensitivity.d_qoi["inclination_degrees"],
        0.0
    );
}
