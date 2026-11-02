use super::*;
use fs_qty::{Length, Pressure, VolumetricFlowRate};

fn state() -> SlurryConditions {
    SlurryConditions {
        temperature: Temperature::new(298.15),
        solids_mass_fraction: Dimensionless::new(0.4),
    }
}

fn domain(minimum_rate: f64) -> SlurryDomain {
    SlurryDomain::new(
        [Frequency::new(minimum_rate), Frequency::new(1_000.0)],
        [Temperature::new(293.15), Temperature::new(303.15)],
        [Dimensionless::new(0.35), Dimensionless::new(0.45)],
    )
    .unwrap()
}

fn context() -> SlurryTestContext {
    SlurryTestContext {
        specimen: "authored homogeneous suspension, not a measured ore slurry".into(),
        method: "exact scalar rheology and pipe reference fixture".into(),
        provenance: Provenance {
            source: "authored arithmetic fixture; constitutive form NIST SP 946 §6".into(),
            license: "project test fixture".into(),
            artifact: None,
        },
    }
}

fn card(yield_stress: f64, reference_stress: f64, n: f64) -> HerschelBulkleyCard {
    HerschelBulkleyCard::new(
        HerschelBulkleyParameters {
            yield_stress: Stress::new(yield_stress),
            reference_excess_stress: Stress::new(reference_stress),
            reference_shear_rate: Frequency::new(10.0),
            flow_index: Dimensionless::new(n),
        },
        domain(0.0),
        context(),
    )
    .unwrap()
}

fn pipe() -> LaminarCircularPipe {
    LaminarCircularPipe::new(Length::new(0.05), Length::new(10.0), Frequency::new(500.0)).unwrap()
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-100),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn signed_constitutive_law_yield_interval_and_dissipation() {
    let hb = card(7.0, 12.0, 0.5);
    for sign in [-1.0, 1.0] {
        // gamma=40/s: tau = 7 + 12*sqrt(40/10) = 31 Pa.
        let response = hb
            .stress_at_shear_rate(Frequency::new(sign * 40.0), state())
            .unwrap();
        let ShearResponse::Yielded {
            stress,
            dissipation,
        } = response
        else {
            panic!()
        };
        close(stress.value(), sign * 31.0, 2e-14);
        close(dissipation.value(), 1_240.0, 2e-14);
        close(
            hb.shear_rate_at_stress(stress, state()).unwrap().value(),
            sign * 40.0,
            2e-14,
        );
    }
    assert_eq!(
        hb.stress_at_shear_rate(Frequency::new(0.0), state())
            .unwrap(),
        ShearResponse::Unyielded {
            stress_bounds: [Stress::new(-7.0), Stress::new(7.0)]
        }
    );
    for stress in [-7.0, -2.0, 0.0, 2.0, 7.0] {
        assert_eq!(
            hb.shear_rate_at_stress(Stress::new(stress), state())
                .unwrap()
                .value(),
            0.0
        );
    }
    assert_eq!(hb.test_context(), &context());
    for n in [0.5, 1.0, 1.4] {
        let hb = card(7.0, 12.0, n);
        let mut last = 0.0;
        for gamma in [0.001, 0.1, 1.0, 10.0, 100.0, 1_000.0] {
            let ShearResponse::Yielded {
                stress,
                dissipation,
            } = hb
                .stress_at_shear_rate(Frequency::new(gamma), state())
                .unwrap()
            else {
                panic!()
            };
            assert!(stress.value() > last && dissipation.value() > 0.0);
            close(
                hb.shear_rate_at_stress(stress, state()).unwrap().value(),
                gamma,
                2e-10,
            );
            last = stress.value();
        }
    }
}

#[test]
fn newtonian_pipe_reproduces_poiseuille_and_pressure_inverse() {
    let fluid = card(0.0, 5.0, 1.0); // mu = 5 Pa / (10/s) = 0.5 Pa.s.
    let pipe = pipe();
    let expected = core::f64::consts::PI * 0.05_f64.powi(4) * 1_000.0 / (8.0 * 0.5 * 10.0);
    for direction in [-1.0, 1.0] {
        let flow = pipe
            .flow_at_pressure(&fluid, state(), Pressure::new(direction * 1_000.0))
            .unwrap();
        close(flow.volumetric_flow.value(), direction * expected, 1e-14);
        close(flow.wall_stress.value(), direction * 2.5, 1e-14);
        close(flow.wall_shear_rate.value(), direction * 5.0, 1e-14);
        assert_eq!(flow.plug_radius.value(), 0.0);
        close(flow.dissipated_power.value(), 1_000.0 * expected, 1e-14);
        let inverse = pipe
            .pressure_for_flow(
                &fluid,
                state(),
                VolumetricFlowRate::new(direction * expected),
                VolumetricFlowRate::new(1e-14),
                80,
            )
            .unwrap();
        close(
            inverse.flow.pressure_drop.value(),
            direction * 1_000.0,
            3e-11,
        );
        assert!(inverse.flow_residual.value().abs() <= 1e-14);
    }
}

#[test]
fn bingham_pipe_reproduces_buckingham_reiner_and_plug_threshold() {
    let fluid = card(10.0, 5.0, 1.0); // yield=10 Pa, plastic viscosity=0.5 Pa.s.
    let pipe = pipe();
    for dp in [-4_000.0, -100.0, 0.0, 100.0, 4_000.0] {
        let rest = pipe
            .flow_at_pressure(&fluid, state(), Pressure::new(dp))
            .unwrap();
        assert_eq!(rest.volumetric_flow.value(), 0.0);
        assert_eq!(rest.wall_shear_rate.value(), 0.0);
        assert_eq!(rest.plug_radius.value(), 0.05);
        assert_eq!(rest.dissipated_power.value(), 0.0);
    }
    let dp = 8_000.0;
    let a = 0.5_f64; // yield / wall stress = 10 / 20.
    let expected = core::f64::consts::PI * 0.05_f64.powi(4) * dp / (8.0 * 0.5 * 10.0)
        * (1.0 - 4.0 * a / 3.0 + a.powi(4) / 3.0);
    let flow = pipe
        .flow_at_pressure(&fluid, state(), Pressure::new(dp))
        .unwrap();
    close(flow.volumetric_flow.value(), expected, 2e-14);
    close(flow.plug_radius.value(), 0.025, 1e-14);
    close(flow.wall_shear_rate.value(), 20.0, 1e-14);
    // A one-ULP stress increment above yield must not disappear in 1 - tau_y/tau_w.
    let barely = pipe
        .flow_at_pressure(&fluid, state(), Pressure::new(4_000.0_f64.next_up()))
        .unwrap();
    assert!(barely.volumetric_flow.value() > 0.0);
    let inverse = pipe
        .pressure_for_flow(
            &fluid,
            state(),
            VolumetricFlowRate::new(expected),
            VolumetricFlowRate::new(1e-13),
            80,
        )
        .unwrap();
    close(inverse.flow.pressure_drop.value(), dp, 2e-10);
    let zero = pipe
        .pressure_for_flow(
            &fluid,
            state(),
            VolumetricFlowRate::new(0.0),
            VolumetricFlowRate::new(1e-13),
            80,
        )
        .unwrap();
    assert_eq!(zero.flow.pressure_drop.value(), 0.0);
}

#[test]
fn herschel_bulkley_pipe_matches_independent_radial_quadrature() {
    // Integrate Q = pi integral gamma(r) r² dr directly, without the implementation's
    // integrated polynomial. The plug contributes through the no-slip velocity integral.
    let radius = 0.05_f64;
    let wall_stress = 40.0;
    for n in [0.5_f64, 0.8, 1.4] {
        let fluid = card(7.0, 12.0, n);
        let actual = pipe()
            .flow_at_pressure(&fluid, state(), Pressure::new(16_000.0))
            .unwrap();
        let cells = 100_000;
        let dr = radius / f64::from(cells);
        let mut integral = 0.0;
        for i in 0..cells {
            let r = (f64::from(i) + 0.5) * dr;
            let excess = wall_stress * r / radius - 7.0;
            let gamma = if excess <= 0.0 {
                0.0
            } else {
                10.0 * (excess / 12.0).powf(1.0 / n)
            };
            integral += gamma * r * r * dr;
        }
        close(
            actual.volumetric_flow.value(),
            core::f64::consts::PI * integral,
            2e-8,
        );
        close(actual.plug_radius.value(), radius * 7.0 / 40.0, 1e-14);
    }
}

#[test]
fn source_limits_cover_actual_state_and_every_sheared_radius() {
    let fluid = card(7.0, 12.0, 0.5);
    for state in [
        SlurryConditions {
            temperature: Temperature::new(290.0),
            ..state()
        },
        SlurryConditions {
            solids_mass_fraction: Dimensionless::new(0.5),
            ..state()
        },
    ] {
        assert!(
            fluid
                .stress_at_shear_rate(Frequency::new(10.0), state)
                .is_err()
        );
        assert!(
            pipe()
                .flow_at_pressure(&fluid, state, Pressure::new(16_000.0))
                .is_err()
        );
    }
    assert!(
        fluid
            .stress_at_shear_rate(Frequency::new(1_001.0), state())
            .is_err()
    );
    assert!(
        fluid
            .shear_rate_at_stress(Stress::new(200.0), state())
            .is_err()
    );
    let only_measured_rates =
        HerschelBulkleyCard::new(fluid.parameters(), domain(1.0), context()).unwrap();
    let ShearResponse::Yielded { stress, .. } = only_measured_rates
        .stress_at_shear_rate(Frequency::new(1.0), state())
        .unwrap()
    else {
        panic!()
    };
    assert_eq!(
        only_measured_rates
            .shear_rate_at_stress(stress, state())
            .unwrap()
            .value(),
        1.0
    );
    assert!(
        only_measured_rates
            .shear_rate_at_stress(Stress::new(stress.value().next_down()), state())
            .is_err()
    );
    assert!(
        only_measured_rates
            .stress_at_shear_rate(Frequency::new(10.0), state())
            .is_ok()
    );
    assert!(
        pipe()
            .flow_at_pressure(&only_measured_rates, state(), Pressure::new(16_000.0))
            .is_err()
    );
    let laminar_limit =
        LaminarCircularPipe::new(Length::new(0.05), Length::new(10.0), Frequency::new(20.0))
            .unwrap();
    assert!(
        laminar_limit
            .flow_at_pressure(&fluid, state(), Pressure::new(16_000.0))
            .is_err()
    );
    assert!(
        pipe()
            .pressure_for_flow(
                &fluid,
                state(),
                VolumetricFlowRate::new(1.0),
                VolumetricFlowRate::new(1e-12),
                80
            )
            .is_err()
    );
    assert!(matches!(
        pipe().pressure_for_flow(
            &fluid,
            state(),
            VolumetricFlowRate::new(0.001),
            VolumetricFlowRate::new(1e-30),
            1
        ),
        Err(SlurryError::NoConvergence { .. })
    ));
}

#[test]
fn invalid_declarations_and_nonfinite_queries_refuse() {
    let good = card(7.0, 12.0, 0.5);
    for bad in [
        HerschelBulkleyParameters {
            yield_stress: Stress::new(-1.0),
            ..good.parameters()
        },
        HerschelBulkleyParameters {
            reference_excess_stress: Stress::new(0.0),
            ..good.parameters()
        },
        HerschelBulkleyParameters {
            flow_index: Dimensionless::new(f64::NAN),
            ..good.parameters()
        },
        HerschelBulkleyParameters {
            reference_shear_rate: Frequency::new(0.0),
            ..good.parameters()
        },
    ] {
        assert!(HerschelBulkleyCard::new(bad, domain(0.0), context()).is_err());
    }
    let mut missing = context();
    missing.method.clear();
    assert!(HerschelBulkleyCard::new(good.parameters(), domain(0.0), missing).is_err());
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(
            good.stress_at_shear_rate(Frequency::new(bad), state())
                .is_err()
        );
        assert!(
            good.shear_rate_at_stress(Stress::new(bad), state())
                .is_err()
        );
        assert!(
            pipe()
                .flow_at_pressure(&good, state(), Pressure::new(bad))
                .is_err()
        );
    }
    assert!(
        LaminarCircularPipe::new(Length::new(0.0), Length::new(10.0), Frequency::new(1.0)).is_err()
    );
    // Poiseuille gives a positive, representable flow at this pressure, while
    // dp*Q is about 4.91e-407 W and cannot be represented as positive f64.
    assert!(matches!(
        pipe().flow_at_pressure(&card(0.0, 5.0, 1.0), state(), Pressure::new(1e-200)),
        Err(SlurryError::Invalid(
            "hydraulic dissipation is unrepresentable"
        ))
    ));
}
