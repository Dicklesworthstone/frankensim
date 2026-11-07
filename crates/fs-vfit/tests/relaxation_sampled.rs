//! The sampled transfer must describe the actual power-conjugate runtime port.
use fs_math::c64::C64;
use fs_vfit::impedance::SeriesImpedanceSpec;
use fs_vfit::relaxation::{
    ComplianceRelaxationTerm, RelaxationImpedance, RelaxationImpedanceSpec, RelaxationSampling,
    RelaxationTerm,
};
use fs_vfit::waveguide::network::admittance::{
    AdmittanceTerm, RelaxationAdmittance, RelaxationAdmittanceSpec,
};

fn base() -> SeriesImpedanceSpec {
    SeriesImpedanceSpec {
        resistance_pa_s_m3: 0.2,
        inertance_pa_s2_m3: 0.002,
        compliance_m3_pa: Some(0.3),
    }
}

fn mixed(sampling: RelaxationSampling) -> RelaxationImpedanceSpec {
    RelaxationImpedanceSpec::new(
        base(),
        &[
            RelaxationTerm {
                resistance_pa_s_m3: 0.6,
                rate_per_s: 80.0,
            },
            RelaxationTerm {
                resistance_pa_s_m3: 0.4,
                rate_per_s: 40_000.0,
            },
        ],
    )
    .unwrap()
    .with_compliance_terms(&[
        ComplianceRelaxationTerm {
            resistance_pa_s_m3: 0.7,
            rate_per_s: 300.0,
        },
        ComplianceRelaxationTerm {
            resistance_pa_s_m3: 0.3,
            rate_per_s: 20_000.0,
        },
    ])
    .unwrap()
    .with_sampling(sampling)
}

fn close_complex(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "actual={actual:?}, expected={expected:?}"
    );
}

#[test]
fn sampled_pressure_matches_harmonic_stepping_and_independent_heat_for_both_rules() {
    let dt = 0.001;
    for sampling in [
        RelaxationSampling::BilinearMidpoint,
        RelaxationSampling::ExactHeldFlow,
    ] {
        for period in [128, 16, 4] {
            let spec = mixed(sampling);
            let mut load = RelaxationImpedance::new(spec, 1.3, dt).unwrap();
            let frequency = 1.0 / (f64::from(period) * dt);
            let expected = spec.sampled_impedance(frequency, dt).unwrap();
            assert_eq!(expected, load.sampled_impedance(frequency).unwrap());
            let mut pressure = C64::ZERO;
            let mut flow = C64::ZERO;
            let mut heat = 0.0;
            let measured = 64 * period;
            // The base C=0.3 and static series/port resistance retain a slow
            // charging mode near 0.75 s; settle that physical history as well
            // as the much faster RL/RC arms before measuring a stationary tone.
            let warm = 32_768;
            for n in 0..warm + measured {
                let angle = core::f64::consts::TAU * f64::from(n) / f64::from(period);
                let (sin, cos) = angle.sin_cos();
                let frame = load.step(cos).unwrap();
                if n >= warm {
                    let demodulate = C64::new(cos, -sin);
                    pressure = pressure + demodulate.scale(frame.port.pressure_pa);
                    flow = flow + demodulate.scale(frame.port.flow_m3_s);
                    heat += frame.port.dissipated_energy_j;
                }
            }
            close_complex(pressure / flow, expected, 2e-10);
            let amplitude = flow.scale(2.0 / f64::from(measured)).abs();
            let predicted_heat =
                0.5 * expected.re * amplitude * amplitude * dt * f64::from(measured);
            assert!(
                (heat - predicted_heat).abs() < 3e-10 * heat,
                "sampling={sampling:?}, heat={heat:e}, predicted={predicted_heat:e}"
            );
        }
    }
}

#[test]
fn midpoint_response_is_the_original_bilinear_impedance_and_queries_keep_history() {
    let dt = 0.001;
    let spec = mixed(RelaxationSampling::BilinearMidpoint);
    let mut load = RelaxationImpedance::new(spec, 1.3, dt).unwrap();
    load.step(0.7).unwrap();
    let mut unchanged = load;
    for frequency in [0.03, 5.0, 100.0, 490.0] {
        let omega = 2.0 * (core::f64::consts::PI * frequency * dt).tan() / dt;
        let s = C64::new(0.0, omega);
        let mut expected = C64::from_re(base().resistance_pa_s_m3)
            + s.scale(base().inertance_pa_s2_m3)
            + s.scale(base().compliance_m3_pa.unwrap()).recip();
        for term in spec.terms() {
            expected =
                expected + (s / (s + C64::from_re(term.rate_per_s))).scale(term.resistance_pa_s_m3);
        }
        for term in spec.compliance_terms() {
            expected = expected
                + (s + C64::from_re(term.rate_per_s))
                    .recip()
                    .scale(term.resistance_pa_s_m3 * term.rate_per_s);
        }
        close_complex(load.sampled_impedance(frequency).unwrap(), expected, 5e-13);
    }
    for frequency in [0.0, 500.0, -1.0, 501.0, f64::NAN, f64::INFINITY] {
        assert!(load.sampled_impedance(frequency).is_err());
    }
    for bad_dt in [0.0, -dt, f64::NAN, f64::INFINITY] {
        assert!(spec.sampled_impedance(5.0, bad_dt).is_err());
    }
    assert_eq!(load.step(-0.2).unwrap(), unchanged.step(-0.2).unwrap());
}

#[test]
fn complementary_branches_keep_static_limits_and_weak_positive_nyquist_heat() {
    let none = SeriesImpedanceSpec {
        resistance_pa_s_m3: 0.0,
        inertance_pa_s2_m3: 0.0,
        compliance_m3_pa: None,
    };
    for sampling in [
        RelaxationSampling::BilinearMidpoint,
        RelaxationSampling::ExactHeldFlow,
    ] {
        for rate in [1e-12, 0.01, 1.0, 1000.0] {
            let high = RelaxationImpedanceSpec::new(
                none,
                &[RelaxationTerm {
                    resistance_pa_s_m3: 2.0,
                    rate_per_s: rate,
                }],
            )
            .unwrap()
            .with_sampling(sampling);
            let low = RelaxationImpedanceSpec::new(none, &[])
                .unwrap()
                .with_compliance_terms(&[ComplianceRelaxationTerm {
                    resistance_pa_s_m3: 2.0,
                    rate_per_s: rate,
                }])
                .unwrap()
                .with_sampling(sampling);
            assert_eq!(high.sampled_impedance(0.0, 1.0).unwrap(), C64::ZERO);
            assert_eq!(low.sampled_impedance(0.0, 1.0).unwrap(), C64::from_re(2.0));
            for frequency in [1e-14, 0.001, 0.1, 0.49, 0.5] {
                let a = high.sampled_impedance(frequency, 1.0).unwrap();
                let b = low.sampled_impedance(frequency, 1.0).unwrap();
                assert!(a.re >= 0.0 && b.re >= 0.0);
                close_complex(a + b, C64::from_re(2.0), 5e-15);
            }
        }
    }
    // At Nyquist, the exact mean RC resistance is R*(1-2*tanh(r/2)/r).
    // Its leading R*r^2/12 stays representable although 1-2h would round to 0.
    for (resistance, rate) in [(1.0, 1e-6), (1.0, 1e-12), (1e100, 1e-100)] {
        let low = RelaxationImpedanceSpec::new(none, &[])
            .unwrap()
            .with_compliance_terms(&[ComplianceRelaxationTerm {
                resistance_pa_s_m3: resistance,
                rate_per_s: rate,
            }])
            .unwrap()
            .with_sampling(RelaxationSampling::ExactHeldFlow);
        let value = low.sampled_impedance(0.5, 1.0).unwrap();
        let expected = ((resistance * rate) * rate) / 12.0;
        assert!(value.re > 0.0 && value.im == 0.0);
        assert!((value.re - expected).abs() < 2e-12 * expected);
    }
}

#[test]
fn sampled_admittance_is_the_same_mean_flow_as_the_actual_dual_owner() {
    let dt = 0.001;
    for sampling in [
        RelaxationSampling::BilinearMidpoint,
        RelaxationSampling::ExactHeldFlow,
    ] {
        let spec = RelaxationAdmittanceSpec::new(
            0.2,
            0.003,
            &[
                AdmittanceTerm {
                    conductance_m3_pa_s: 0.6,
                    rate_per_s: 80.0,
                },
                AdmittanceTerm {
                    conductance_m3_pa_s: 0.4,
                    rate_per_s: 40_000.0,
                },
            ],
        )
        .unwrap()
        .with_sampling(sampling);
        let mut load = RelaxationAdmittance::new(spec, 1.3, dt).unwrap();
        let (mut pressure, mut flow) = (C64::ZERO, C64::ZERO);
        for n in 0..8192 {
            let (sin, cos) = (core::f64::consts::TAU * f64::from(n) / 32.0).sin_cos();
            let frame = load.step(cos).unwrap();
            if n >= 4096 {
                let demodulate = C64::new(cos, -sin);
                pressure = pressure + demodulate.scale(frame.pressure_pa);
                flow = flow + demodulate.scale(frame.flow_m3_s);
            }
        }
        let expected = spec.sampled_admittance(1.0 / (32.0 * dt), dt).unwrap();
        assert_eq!(
            expected,
            load.sampled_admittance(1.0 / (32.0 * dt)).unwrap()
        );
        close_complex(flow / pressure, expected, 2e-10);
    }
}
