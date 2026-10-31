//! G0/G3 checks of the aiding mixed-plate formula, limits and local derivatives.
//! Formula arithmetic is verified here; this is not experimental validation.

use fs_convection::{CorrelationError, CorrelationId, CorrelationInputs, evaluate};

const CARD: CorrelationId = CorrelationId::ManningQureshiAidingVerticalPlate;

fn nusselt(re: f64, ra: f64, pr: f64) -> f64 {
    evaluate(CARD, CorrelationInputs::mixed(re, ra, pr))
        .expect("admitted aiding plate")
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
fn source_equations_match_independent_high_precision_arithmetic() {
    // Test-owned 60-decimal-digit evaluation of Manning & Qureshi (1994)
    // Eqs. (4), (8), (10), (11) using Python Decimal ln/exp and the original
    // Fn*(Ra/Pr)^(1/4) expression. These are formula checks, not source-table
    // measurements or observations of this Rust implementation.
    for (pr, expected) in [
        (0.003, 5.989_106_929_891_576),
        (0.72, 21.570_096_902_182_86),
        (1.0, 23.419_283_913_172_51),
        (10.0, 46.330_632_705_354_51),
        (100.0, 99.455_679_605_428_89),
    ] {
        relative_close(nusselt(1_000.0, 1.0e6, pr), expected, 2.0e-13);
    }
    relative_close(nusselt(0.0, 1.0e6, 0.72), 16.308_419_837_353_866, 2.0e-13);
    relative_close(nusselt(1_000.0, 0.0, 0.72), 18.709_186_006_703_454, 2.0e-13);

    let value = evaluate(CARD, CorrelationInputs::mixed(1_000.0, 1.0e6, 0.72)).unwrap();
    assert_eq!(value.card().source.identifier, "doi:10.1115/1.2910940");
    for disclosure in [
        "implementation policy",
        "frozen-property approximation",
        "conduction intercept",
    ] {
        assert!(
            value
                .card()
                .model
                .assumptions
                .iter()
                .any(|text| text.contains(disclosure))
        );
    }
    assert!(
        value
            .card()
            .model
            .known_failures
            .iter()
            .any(|text| text.contains("opposing"))
    );
}

#[test]
fn each_flow_driver_approaches_its_pure_limit_continuously() {
    for pr in [0.003, 0.72, 10.0, 100.0] {
        let natural = nusselt(0.0, 1.0e6, pr);
        let forced = nusselt(1_000.0, 0.0, pr);
        let mixed = nusselt(1_000.0, 1.0e6, pr);
        assert!(mixed > natural.max(forced));
        assert!(mixed < natural + forced);

        let mut previous = mixed;
        for re in [100.0, 1.0, 1.0e-4, 1.0e-12, 0.0] {
            let next = nusselt(re, 1.0e6, pr);
            assert!(next >= natural && next <= previous);
            previous = next;
        }
        relative_close(nusselt(1.0e-12, 1.0e6, pr), natural, 2.0e-13);

        let mut previous = mixed;
        for ra in [1.0e4, 1.0, 1.0e-8, 1.0e-24, 0.0] {
            let next = nusselt(1_000.0, ra, pr);
            assert!(next >= forced && next <= previous);
            previous = next;
        }
        relative_close(nusselt(1_000.0, 1.0e-24, pr), forced, 2.0e-13);
    }
}

#[test]
fn normalized_power_sum_retains_tiny_limits_and_homogeneity() {
    // Multiplying Re by s^2 and Ra by s^4 multiplies both components and Nu
    // by s. Direct component powers underflow for the smallest case at Pr=100.
    for pr in [0.003, 0.72, 100.0] {
        let nominal = nusselt(1_000.0, 1.0e6, pr);
        for scale in [1.0e-2_f64, 1.0e-30, 1.0e-70] {
            let scaled = nusselt(1_000.0 * scale.powi(2), 1.0e6 * scale.powi(4), pr);
            relative_close(scaled / scale, nominal, 3.0e-13);
        }
        let tiny = f64::from_bits(1);
        relative_close(nusselt(tiny, 1.0e6, pr), nusselt(0.0, 1.0e6, pr), 2.0e-13);
        relative_close(
            nusselt(1_000.0, tiny, pr),
            nusselt(1_000.0, 0.0, pr),
            2.0e-13,
        );
        assert!(nusselt(tiny, tiny, pr).is_finite());
    }
}

#[test]
fn strict_laminar_bounds_and_declared_prandtl_policy_refuse_extrapolation() {
    for pr in [0.003, 100.0] {
        assert!(
            evaluate(
                CARD,
                CorrelationInputs::mixed(5.0e5_f64.next_down(), 1.0e9_f64.next_down(), pr)
            )
            .is_ok()
        );
    }
    for inputs in [
        CorrelationInputs::mixed(5.0e5, 1.0e6, 0.72),
        CorrelationInputs::mixed(1_000.0, 1.0e9, 0.72),
        CorrelationInputs::mixed(1_000.0, 1.0e6, 0.003_f64.next_down()),
        CorrelationInputs::mixed(1_000.0, 1.0e6, 100.0_f64.next_up()),
        CorrelationInputs::forced(1_000.0, 0.72),
        CorrelationInputs::natural(1.0e6, 0.72),
    ] {
        assert!(matches!(
            evaluate(CARD, inputs),
            Err(CorrelationError::OutOfDomain { .. })
        ));
    }
    for inputs in [
        CorrelationInputs::mixed(-1.0, 1.0e6, 0.72),
        CorrelationInputs::mixed(1_000.0, -1.0, 0.72),
        CorrelationInputs::mixed(f64::NAN, 1.0e6, 0.72),
        CorrelationInputs::mixed(1_000.0, f64::INFINITY, 0.72),
        CorrelationInputs::mixed(1_000.0, 1.0e6, 0.0),
    ] {
        assert!(matches!(
            evaluate(CARD, inputs),
            Err(CorrelationError::InvalidGroup { .. })
        ));
    }
    assert!(matches!(
        evaluate(CARD, CorrelationInputs::mixed(0.0, 0.0, 0.72)),
        Err(CorrelationError::InvalidGroup {
            axis: "max(Re,Ra)",
            ..
        })
    ));
    // Adding a new card's pure limits must not relax other models' gates.
    assert!(matches!(
        evaluate(
            CorrelationId::FlatPlateLaminarAverage,
            CorrelationInputs::forced(0.0, 0.72)
        ),
        Err(CorrelationError::InvalidGroup { axis: "Re", .. })
    ));
    assert!(matches!(
        evaluate(
            CorrelationId::ChurchillChuVerticalPlate,
            CorrelationInputs::natural(0.0, 0.72)
        ),
        Err(CorrelationError::InvalidGroup { axis: "Ra", .. })
    ));
}

#[test]
fn analytic_partials_include_the_prandtl_dependent_mixing_exponent() {
    let step = 1.0e-5_f64;
    for pr in [0.004, 0.01, 0.72, 1.0, 10.0, 99.0] {
        for (re, ra) in [(1.0, 1.0e6), (1_000.0, 1.0e6), (1.0e5, 1.0)] {
            let sensitive = CorrelationInputs::mixed(re, ra, pr)
                .evaluate_with_sensitivities(CARD)
                .unwrap();
            let plain = evaluate(CARD, CorrelationInputs::mixed(re, ra, pr)).unwrap();
            assert_eq!(
                sensitive.evidence().value.to_bits(),
                plain.evidence().value.to_bits()
            );
            assert_eq!(sensitive.evidence().model, plain.evidence().model);
            for (index, axis) in ["Re", "Ra", "Pr"].into_iter().enumerate() {
                let inputs = [re, ra, pr];
                let evaluate_scaled = |factor: f64| {
                    let mut inputs = inputs;
                    inputs[index] *= factor;
                    nusselt(inputs[0], inputs[1], inputs[2])
                };
                let expected = (evaluate_scaled(step.exp()).ln()
                    - evaluate_scaled((-step).exp()).ln())
                    / (2.0 * step);
                let actual = sensitive.evidence().sensitivity.d_qoi[axis] * inputs[index]
                    / sensitive.evidence().value;
                assert!(
                    (actual - expected).abs() < 1.0e-8,
                    "Pr={pr}, Re={re}, Ra={ra}, {axis}: {actual} != {expected}"
                );
                if axis == "Re" {
                    assert!(
                        (plain.reynolds_elasticity().unwrap().unwrap() - actual).abs() < 2.0e-13
                    );
                }
            }
        }
    }
}

#[test]
fn pure_limits_omit_only_the_endpoint_partial_and_preserve_tiny_derivatives() {
    for (re, ra, missing, retained) in [(0.0, 1.0e6, "Re", "Ra"), (1_000.0, 0.0, "Ra", "Re")] {
        let value = CorrelationInputs::mixed(re, ra, 0.72)
            .evaluate_with_sensitivities(CARD)
            .unwrap();
        let partials = &value.evidence().sensitivity.d_qoi;
        assert!(!partials.contains_key(missing));
        assert!(partials[retained] > 0.0);
        assert!(partials["Pr"] > 0.0);
    }
    // At large Pr, n/4 > 1. A tiny Ra contribution's power weight can
    // underflow even though dNu/dRa stays representable. In the forced limit,
    // that derivative scales as Ra^(n/4-1); test this independently of Nu's
    // rounded-to-forced primal value.
    let derivative = |ra| {
        CorrelationInputs::mixed(1_000.0, ra, 100.0)
            .evaluate_with_sensitivities(CARD)
            .unwrap()
            .evidence()
            .sensitivity
            .d_qoi["Ra"]
    };
    let small = 1.0e-200;
    let tiny = 1.0e-320;
    let exponent = 3.5 * 100.0_f64.powf(0.075);
    let expected = derivative(small) * (tiny / small).powf(exponent / 4.0 - 1.0);
    assert!(derivative(tiny) > 0.0);
    relative_close(derivative(tiny), expected, 2.0e-12);
}
