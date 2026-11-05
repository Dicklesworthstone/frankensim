//! G0/G1/G3: authored interface physics, weak signals and source-domain refusal.
//! Coefficients below are mathematical fixtures, not measured surface data.

use fs_material::superfluid::kapitza::{KapitzaInterface, KapitzaParameters};

fn parameters(exponent: f64) -> KapitzaParameters {
    KapitzaParameters {
        area_m2: 0.125,
        coefficient_w_m2_k_n: 2.0,
        exponent,
        min_temperature_k: 0.25,
        max_temperature_k: 4.0,
    }
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(f64::MIN_POSITIVE),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn quartic_interface_matches_independent_power_and_entropy_balances() {
    let law = KapitzaInterface::try_new(parameters(4.0)).unwrap();
    assert_eq!(law.parameters(), parameters(4.0));
    assert_eq!(law.temperature_bounds_k(), [0.25, 4.0]);
    assert_eq!(law.heat_flow_w(2.0, 1.0).unwrap(), 3.75);
    assert_eq!(law.heat_flow_w(1.0, 2.0).unwrap(), -3.75);
    assert_eq!(law.tangent_from_w_k(2.0).unwrap(), 8.0);
    assert_eq!(law.tangent_from_w_k(1.0).unwrap(), 1.0);
    assert_eq!(law.entropy_production_w_k(2.0, 1.0).unwrap(), 1.875);
    assert_eq!(law.entropy_production_w_k(1.0, 2.0).unwrap(), 1.875);
    assert_eq!(law.heat_flow_w(2.0, 2.0).unwrap(), 0.0);
    assert_eq!(law.entropy_production_w_k(2.0, 2.0).unwrap(), 0.0);
}

#[test]
fn fractional_and_integer_laws_match_independent_wide_separation_oracles() {
    for n in [1.0, 2.0, 2.7, 3.0, 3.5, 4.0, 4.25, 5.0] {
        let law = KapitzaInterface::try_new(parameters(n)).unwrap();
        for (hot, cold) in [(4.0_f64, 0.25_f64), (2.0, 1.0), (1.8, 1.7)] {
            let expected = 0.25 * (hot.powf(n) - cold.powf(n));
            let power = law.heat_flow_w(hot, cold).unwrap();
            close(power, expected, 2e-13);
            assert_eq!(law.heat_flow_w(cold, hot).unwrap(), -power);
            close(
                law.tangent_from_w_k(hot).unwrap(),
                0.25 * n * hot.powf(n - 1.0),
                2e-13,
            );
            assert!(law.entropy_production_w_k(hot, cold).unwrap() > 0.0);
            assert!(law.tangent_from_w_k(cold).unwrap() >= law.minimum_tangent_w_k());
        }
    }
}

#[test]
fn one_ulp_temperature_difference_retains_differential_heat_transfer() {
    let cold = 1.8_f64;
    let hot = f64::from_bits(cold.to_bits() + 1);
    for n in [1.0, 2.0, 2.7, 3.0, 3.5, 4.0, 4.25] {
        let law = KapitzaInterface::try_new(parameters(n)).unwrap();
        let expected = 0.25 * n * cold.powf(n - 1.0) * (hot - cold);
        let power = law.heat_flow_w(hot, cold).unwrap();
        assert!(power > 0.0, "lost weak heat transfer for n={n}");
        close(power, expected, 2e-13);
        assert_eq!(law.heat_flow_w(cold, hot).unwrap(), -power);
    }
}

#[test]
fn reference_offsets_preserve_heat_below_an_absolute_temperature_ulp() {
    let reference = 1.8_f64;
    for n in [1.0, 3.0, 3.5, 4.0] {
        let law = KapitzaInterface::try_new(parameters(n)).unwrap();
        for offset in [1e-18, 1e-100, 1e-300] {
            assert_eq!(reference + offset, reference);
            let expected = 0.25 * n * reference.powf(n - 1.0) * offset;
            let power = law
                .heat_flow_about_reference_w(reference, offset, 0.0)
                .unwrap();
            assert!(power > 0.0);
            close(power, expected, 2e-13);
            assert_eq!(
                law.heat_flow_about_reference_w(reference, 0.0, offset)
                    .unwrap(),
                -power,
            );
        }
        let from = 0.1_f64;
        let to = -0.2_f64;
        close(
            law.heat_flow_about_reference_w(reference, from, to)
                .unwrap(),
            law.heat_flow_w(reference + from, reference + to).unwrap(),
            2e-13,
        );
    }
    let law = KapitzaInterface::try_new(parameters(3.5)).unwrap();
    assert!(
        law.heat_flow_about_reference_w(2.0, f64::from_bits(1), 0.0)
            .unwrap()
            > 0.0
    );
}

#[test]
fn sub_ulp_offsets_cannot_round_an_invalid_temperature_onto_a_source_boundary() {
    let law = KapitzaInterface::try_new(parameters(4.0)).unwrap();
    let tiny = 1e-20;
    assert_eq!(0.25 - tiny, 0.25);
    assert_eq!(4.0 + tiny, 4.0);
    assert!(law.heat_flow_about_reference_w(0.25, -tiny, 0.0).is_err());
    assert!(law.heat_flow_about_reference_w(4.0, 0.0, tiny).is_err());
    assert!(law.heat_flow_about_reference_w(0.25, tiny, 0.0).unwrap() > 0.0);
    assert!(law.heat_flow_about_reference_w(4.0, 0.0, -tiny).unwrap() > 0.0);
    // A common reference outside the source range is allowed when the actual
    // two temperatures, retained as compensated sums, are both inside.
    close(
        law.heat_flow_about_reference_w(5.0, -3.0, -4.0).unwrap(),
        3.75,
        2e-15,
    );
    for reference in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            law.heat_flow_about_reference_w(reference, 0.0, 0.0)
                .is_err()
        );
    }
    assert!(law.heat_flow_about_reference_w(2.0, f64::NAN, 0.0).is_err());
    assert!(
        law.heat_flow_about_reference_w(2.0, 0.0, f64::INFINITY)
            .is_err()
    );
}

#[test]
fn interface_area_and_source_coefficient_scale_actual_power() {
    let p = parameters(3.5);
    let law = KapitzaInterface::try_new(p).unwrap();
    let area = KapitzaInterface::try_new(KapitzaParameters {
        area_m2: 4.0 * p.area_m2,
        ..p
    })
    .unwrap();
    let source = KapitzaInterface::try_new(KapitzaParameters {
        coefficient_w_m2_k_n: 4.0 * p.coefficient_w_m2_k_n,
        ..p
    })
    .unwrap();
    let power = law.heat_flow_w(1.9, 1.7).unwrap();
    assert_eq!(area.heat_flow_w(1.9, 1.7).unwrap(), 4.0 * power);
    assert_eq!(source.heat_flow_w(1.9, 1.7).unwrap(), 4.0 * power);
    close(
        area.minimum_tangent_w_k(),
        4.0 * law.minimum_tangent_w_k(),
        2e-15,
    );
}

#[test]
fn monotone_temperature_response_has_a_positive_lower_tangent() {
    for n in [1.0, 3.0, 3.5, 4.0] {
        let law = KapitzaInterface::try_new(parameters(n)).unwrap();
        let mut previous = law.heat_flow_w(0.25, 2.0).unwrap();
        for k in 1..=60 {
            let t = 0.25 + 3.75 * f64::from(k) / 60.0;
            let next = law.heat_flow_w(t, 2.0).unwrap();
            assert!(next > previous);
            assert!(law.tangent_from_w_k(t).unwrap() >= law.minimum_tangent_w_k());
            previous = next;
        }
    }
}

#[test]
fn source_domains_and_unrepresentable_coefficients_are_refused() {
    for invalid in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        for p in [
            KapitzaParameters {
                area_m2: invalid,
                ..parameters(4.0)
            },
            KapitzaParameters {
                coefficient_w_m2_k_n: invalid,
                ..parameters(4.0)
            },
            KapitzaParameters {
                min_temperature_k: invalid,
                ..parameters(4.0)
            },
            KapitzaParameters {
                max_temperature_k: invalid,
                ..parameters(4.0)
            },
            KapitzaParameters {
                exponent: invalid,
                ..parameters(4.0)
            },
        ] {
            assert!(KapitzaInterface::try_new(p).is_err(), "accepted {p:?}");
        }
    }
    for p in [
        KapitzaParameters {
            exponent: 0.5,
            ..parameters(4.0)
        },
        KapitzaParameters {
            max_temperature_k: 0.25,
            ..parameters(4.0)
        },
        KapitzaParameters {
            area_m2: f64::MAX,
            ..parameters(4.0)
        },
        KapitzaParameters {
            coefficient_w_m2_k_n: f64::MIN_POSITIVE,
            exponent: 1000.0,
            ..parameters(4.0)
        },
    ] {
        assert!(KapitzaInterface::try_new(p).is_err(), "accepted {p:?}");
    }
    let law = KapitzaInterface::try_new(parameters(4.0)).unwrap();
    for invalid in [0.0, 0.249, 4.001, f64::NAN, f64::INFINITY] {
        assert!(law.heat_flow_w(invalid, 2.0).is_err());
        assert!(law.heat_flow_w(2.0, invalid).is_err());
        assert!(law.tangent_from_w_k(invalid).is_err());
        assert!(law.entropy_production_w_k(invalid, 2.0).is_err());
    }
}
