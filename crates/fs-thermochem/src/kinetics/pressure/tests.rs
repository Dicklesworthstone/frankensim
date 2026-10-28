use super::*;

fn colliders() -> ThirdBodyEfficienciesV1 {
    ThirdBodyEfficienciesV1::new(0.0, vec![(id("X"), 1.0)]).unwrap()
}

fn mechanism(
    released_energy: f64,
    high: ArrheniusRateV1,
    pressure: PressureDependenceV1,
) -> IdealGasKineticsV1 {
    IdealGasKineticsV1::new(
        thermo(
            vec![
                source("A", 0.01, 4.0, 0.0, 0.0),
                source("B", 0.01, 4.0, -released_energy, 0.0),
                source("X", 0.04, 2.5, 0.0, 0.0),
            ],
            vec![1, 1, 1],
        ),
        "pressure-kinetics-fixture",
        1,
        vec![reaction("A-B", &[("A", 1)], &[("B", 1)], high).with_pressure_dependence(pressure)],
        1e-12,
    )
    .unwrap()
}

fn troe() -> TroeFalloffV1 {
    TroeFalloffV1::new(
        0.55,
        Temperature::new(100.0),
        Temperature::new(1200.0),
        Some(Temperature::new(5000.0)),
    )
    .unwrap()
}

fn falloff(low: ArrheniusRateV1, broadening: FalloffBroadeningV1) -> PressureDependenceV1 {
    PressureDependenceV1::Falloff {
        low_pressure_rate: low,
        colliders: colliders(),
        broadening,
    }
}

fn arrhenius(law: ArrheniusRateV1, t: f64) -> f64 {
    let tr = law.reference_temperature().value();
    law.reference_rate().value()
        * (t / tr).powf(law.temperature_exponent())
        * (law.activation_energy().value() / R * (1.0 / tr - 1.0 / t)).exp()
}

fn center(law: TroeFalloffV1, t: f64) -> f64 {
    (1.0 - law.alpha()) * (-t / law.temperature_3().value()).exp()
        + law.alpha() * (-t / law.temperature_1().value()).exp()
        + law
            .temperature_2()
            .map_or(0.0, |t2| (-t2.value() / t).exp())
}

// Independent direct textbook expression, deliberately distinct from the
// production logarithmic/scaled rational evaluation and derivative algebra.
fn troe_factor(law: TroeFalloffV1, t: f64, reduced_pressure: f64) -> f64 {
    let l = center(law, t).log10();
    let c = -0.4 - 0.67 * l;
    let n = 0.75 - 1.27 * l;
    let y = reduced_pressure.log10() + c;
    let f1 = y / (n - 0.14 * y);
    10.0_f64.powf(l / (1.0 + f1 * f1))
}

#[test]
fn g1_third_body_efficiencies_and_all_jacobian_columns_match_mass_action() {
    let law = rate(4.0, 0.3, 1200.0);
    let m = mechanism(
        8000.0,
        law,
        PressureDependenceV1::ThirdBody(
            ThirdBodyEfficienciesV1::new(0.25, vec![(id("X"), 4.0), (id("A"), 0.5)]).unwrap(),
        ),
    );
    let c = [2.0, 0.5, 3.0];
    let collision = 0.5 * c[0] + 0.25 * c[1] + 4.0 * c[2];
    let mut w = KineticWorkspaceV1::new(&m).unwrap();
    for t in [300.0, 500.0, 1400.0] {
        let result = m
            .evaluate(Temperature::new(t), &c, &mut w, || false)
            .unwrap();
        let kf = arrhenius(law, t);
        let kr = kf * (-8000.0 / (R * t)).exp();
        let base = kf * c[0] - kr * c[1];
        close(w.reaction_rates()[0], collision * base, 3e-12 * base.abs());
        close(
            w.rate_jacobian()[0],
            collision * kf + 0.5 * base,
            3e-12 * kf,
        );
        close(
            w.rate_jacobian()[1],
            -collision * kr + 0.25 * base,
            3e-12 * kf,
        );
        close(w.rate_jacobian()[2], 4.0 * base, 3e-12 * kf);
        close(
            w.temperature_derivatives()[0],
            collision
                * (base * (0.3 / t + 1200.0 / (R * t * t)) - kr * c[1] * 8000.0 / (R * t * t)),
            1e-12 * kf,
        );
        close(
            result.entropy_production_w_m3_k,
            R * collision * base * (8000.0 / (R * t) + (c[0] / c[1]).ln()),
            2e-10 * kf,
        );
        assert_eq!(w.species_rates()[2], 0.0);
    }
}

#[test]
fn g1_lindemann_limits_and_logarithmic_extremes_keep_representable_rates() {
    let m = mechanism(
        0.0,
        rate(7.0, 0.0, 0.0),
        falloff(rate(2.0, 0.0, 0.0), FalloffBroadeningV1::Lindemann),
    );
    let mut w = KineticWorkspaceV1::new(&m).unwrap();
    for collider in [1e-200, 1e-12, 0.1, 3.5, 1e6, 1e100] {
        m.evaluate(
            Temperature::new(500.0),
            &[2.0, 1.0, collider],
            &mut w,
            || false,
        )
        .unwrap();
        let pr = 2.0 * collider / 7.0;
        let k = 2.0 * collider / (1.0 + pr);
        let dk = 2.0 / ((1.0 + pr) * (1.0 + pr));
        close(w.reaction_rates()[0], k, k * 3e-13);
        close(w.rate_jacobian()[0], k, k * 3e-13);
        close(w.rate_jacobian()[1], -k, k * 3e-13);
        close(w.rate_jacobian()[2], dk, dk * 5e-13);
    }
    // k0*M overflows in the first case and Pr underflows in the second;
    // both complete physical rate and Jacobian remain representable.
    for (high, low, expected_rate, expected_collider_derivative) in
        [(1e200, 1e300, 1e200, 1e80), (1e200, 1e-300, 1e-290, 1e-300)]
    {
        let m = mechanism(
            0.0,
            rate(high, 0.0, 0.0),
            falloff(rate(low, 0.0, 0.0), FalloffBroadeningV1::Lindemann),
        );
        let mut w = KineticWorkspaceV1::new(&m).unwrap();
        m.evaluate(Temperature::new(500.0), &[2.0, 1.0, 1e10], &mut w, || false)
            .unwrap();
        close(w.reaction_rates()[0], expected_rate, expected_rate * 8e-13);
        close(
            w.rate_jacobian()[2],
            expected_collider_derivative,
            expected_collider_derivative * 8e-13,
        );
    }
}

#[test]
fn g1_troe_rates_and_analytic_derivatives_match_independent_perturbations() {
    let high = rate(12.0, 0.2, 9000.0);
    let low = rate(3.0, -0.4, 3000.0);
    let tr = troe();
    let m = mechanism(
        3000.0,
        high,
        PressureDependenceV1::Falloff {
            low_pressure_rate: low,
            colliders: ThirdBodyEfficienciesV1::new(0.3, vec![(id("X"), 2.0)]).unwrap(),
            broadening: FalloffBroadeningV1::Troe(tr),
        },
    );
    let mut w = KineticWorkspaceV1::new(&m).unwrap();
    for t in [250.0, 700.0, 1800.0] {
        for bath in [1e-5, 1.0, 1e5] {
            let c = [1.8, 0.4, bath];
            m.evaluate(Temperature::new(t), &c, &mut w, || false)
                .unwrap();
            let collider = 0.3 * (c[0] + c[1]) + 2.0 * c[2];
            let kh = arrhenius(high, t);
            let kl = arrhenius(low, t);
            let pr = kl * collider / kh;
            let k = kh * pr / (1.0 + pr) * troe_factor(tr, t, pr);
            let expected = k * (c[0] - c[1] * (-3000.0 / (R * t)).exp());
            close(w.reaction_rates()[0], expected, expected.abs() * 2e-12);
            let jac = w.rate_jacobian().to_vec();
            let dt = w.temperature_derivatives()[0];
            for i in 0..3 {
                let h = c[i] * 2e-5;
                let mut plus = c;
                let mut minus = c;
                plus[i] += h;
                minus[i] -= h;
                m.evaluate(Temperature::new(t), &plus, &mut w, || false)
                    .unwrap();
                let rp = w.reaction_rates()[0];
                m.evaluate(Temperature::new(t), &minus, &mut w, || false)
                    .unwrap();
                let rm = w.reaction_rates()[0];
                close(jac[i], (rp - rm) / (2.0 * h), 8e-7 * jac[i].abs().max(1.0));
            }
            let h = t * 1e-5;
            m.evaluate(Temperature::new(t + h), &c, &mut w, || false)
                .unwrap();
            let rp = w.reaction_rates()[0];
            m.evaluate(Temperature::new(t - h), &c, &mut w, || false)
                .unwrap();
            let rm = w.reaction_rates()[0];
            close(dt, (rp - rm) / (2.0 * h), 1e-7 * dt.abs().max(1e-8));
        }
    }
}

#[test]
fn g0_zero_collider_keeps_its_inward_derivative_without_changing_affinity() {
    let tr = troe();
    let f0 = center(tr, 500.0).powf(0.14 * 0.14 / (1.0 + 0.14 * 0.14));
    let laws = [
        (PressureDependenceV1::ThirdBody(colliders()), 7.0),
        (
            falloff(rate(4.0, 0.0, 0.0), FalloffBroadeningV1::Lindemann),
            4.0,
        ),
        (
            falloff(rate(4.0, 0.0, 0.0), FalloffBroadeningV1::Troe(tr)),
            4.0 * f0,
        ),
    ];
    for (law, slope) in laws {
        let m = mechanism(0.0, rate(7.0, 0.0, 0.0), law);
        let mut w = KineticWorkspaceV1::new(&m).unwrap();
        for b in [0.0, 0.5, 2.0] {
            let state = m
                .evaluate(Temperature::new(500.0), &[2.0, b, 0.0], &mut w, || false)
                .unwrap();
            assert_eq!(w.reaction_rates()[0], 0.0);
            assert_eq!(w.species_rates(), &[0.0, 0.0, 0.0]);
            assert_eq!(w.rate_jacobian()[0], 0.0);
            assert_eq!(w.rate_jacobian()[1], 0.0);
            close(w.rate_jacobian()[2], (2.0 - b) * slope, slope * 2e-14);
            assert_eq!(state.entropy_production_w_m3_k, 0.0);
            assert_eq!(w.temperature_derivatives()[0], 0.0);
            if b == 0.0 {
                assert_eq!(w.dimensionless_affinities()[0], Some(f64::INFINITY));
            } else {
                close(
                    w.dimensionless_affinities()[0].unwrap(),
                    (2.0 / b).ln(),
                    3e-15,
                );
            }
        }
    }
}

#[test]
fn g0_troe_apparent_rational_pole_has_smooth_finite_rate_and_derivatives() {
    let tr = TroeFalloffV1::new(
        0.5,
        Temperature::new(1000.0),
        Temperature::new(1000.0),
        None,
    )
    .unwrap();
    let m = mechanism(
        0.0,
        rate(7.0, 0.0, 0.0),
        falloff(rate(2.0, 0.0, 0.0), FalloffBroadeningV1::Troe(tr)),
    );
    let l = center(tr, 500.0).log10();
    let c = -0.4 - 0.67 * l;
    let n = 0.75 - 1.27 * l;
    let pr = 10.0_f64.powf(n / 0.14 - c);
    let collider = pr * 7.0 / 2.0;
    let mut w = KineticWorkspaceV1::new(&m).unwrap();
    m.evaluate(
        Temperature::new(500.0),
        &[2.0, 1.0, collider],
        &mut w,
        || false,
    )
    .unwrap();
    close(w.reaction_rates()[0], 7.0 * pr / (1.0 + pr), 4e-14);
    close(w.rate_jacobian()[2], 2.0 / ((1.0 + pr) * (1.0 + pr)), 1e-25);
    assert!(w.temperature_derivatives()[0].is_finite());
}

#[test]
fn g3_pressure_multipliers_preserve_near_equilibrium_entropy_and_source_identity() {
    for scale in [1e-150, 1.0, 1e150] {
        let tr = troe();
        for (law, factor) in [
            (PressureDependenceV1::ThirdBody(colliders()), 4.0),
            (
                falloff(rate(scale * 0.5, 0.0, 0.0), FalloffBroadeningV1::Lindemann),
                2.0 / 3.0,
            ),
            (
                falloff(rate(scale * 0.5, 0.0, 0.0), FalloffBroadeningV1::Troe(tr)),
                2.0 / 3.0 * troe_factor(tr, 500.0, 2.0),
            ),
        ] {
            let m = mechanism(0.0, rate(scale, 0.0, 0.0), law);
            let mut w = KineticWorkspaceV1::new(&m).unwrap();
            for sign in [-1.0, 1.0] {
                let a = 1.0 + sign * 1e-12;
                let state = m
                    .evaluate(Temperature::new(500.0), &[a, 1.0, 4.0], &mut w, || false)
                    .unwrap();
                close(w.reaction_rates()[0] / scale, factor * (a - 1.0), 2e-24);
                assert!(state.entropy_production_w_m3_k > 0.0);
                assert_eq!(w.log_equilibrium_constants()[0], 0.0);
                assert_eq!(w.species_rates()[2], 0.0);
                assert!(m.conservation().matches(
                    m.thermodynamics().elemental_matrix(),
                    m.stoichiometry(),
                    m.thermodynamics().charge_vector()
                ));
            }
        }
    }
}

#[test]
fn g0_invalid_pressure_data_refuse_and_exact_shutdown_never_evaluates_falloff() {
    assert!(ThirdBodyEfficienciesV1::new(-1.0, vec![]).is_err());
    assert!(ThirdBodyEfficienciesV1::new(1.0, vec![(id("X"), f64::NAN)]).is_err());
    assert!(ThirdBodyEfficienciesV1::new(1.0, vec![(id("X"), 1.0), (id("X"), 2.0)]).is_err());
    assert!(
        TroeFalloffV1::new(1.1, Temperature::new(100.0), Temperature::new(100.0), None).is_err()
    );
    assert!(TroeFalloffV1::new(0.5, Temperature::new(0.0), Temperature::new(100.0), None).is_err());
    assert!(
        TroeFalloffV1::new(
            0.5,
            Temperature::new(100.0),
            Temperature::new(100.0),
            Some(Temperature::new(f64::INFINITY))
        )
        .is_err()
    );
    let source = mechanism(0.0, rate(1.0, 0.0, 0.0), PressureDependenceV1::Independent);
    let unknown = reaction(
        "unknown-collider",
        &[("A", 1)],
        &[("B", 1)],
        rate(1.0, 0.0, 0.0),
    )
    .with_pressure_dependence(PressureDependenceV1::ThirdBody(
        ThirdBodyEfficienciesV1::new(1.0, vec![(id("absent"), 0.0)]).unwrap(),
    ));
    assert!(
        IdealGasKineticsV1::new(
            source.thermodynamics().clone(),
            "unknown",
            1,
            vec![unknown],
            1e-12
        )
        .is_err()
    );
    let singular = TroeFalloffV1::new(
        0.5,
        Temperature::new(f64::from_bits(1)),
        Temperature::new(f64::from_bits(1)),
        None,
    )
    .unwrap();
    for (high, low) in [(0.0, 1.0), (1.0, 0.0)] {
        let m = mechanism(
            0.0,
            rate(high, 0.0, 0.0),
            falloff(rate(low, 0.0, 0.0), FalloffBroadeningV1::Troe(singular)),
        );
        let mut w = KineticWorkspaceV1::new(&m).unwrap();
        m.evaluate(Temperature::new(500.0), &[2.0, 0.0, 4.0], &mut w, || false)
            .unwrap();
        assert_eq!(w.reaction_rates()[0], 0.0);
        assert_eq!(w.rate_jacobian(), &[0.0, 0.0, 0.0]);
        assert_eq!(w.temperature_derivatives()[0], 0.0);
        assert_eq!(w.dimensionless_affinities()[0], Some(f64::INFINITY));
    }
    let unrepresentable_high = ArrheniusRateV1::new(
        Temperature::new(1.0),
        ReactionProgressRateQuantityV1::new(1.0),
        f64::MAX,
        MolarEnergyQuantityV1::new(0.0),
    )
    .unwrap();
    for law in [
        falloff(rate(0.0, 0.0, 0.0), FalloffBroadeningV1::Troe(singular)),
        PressureDependenceV1::ThirdBody(ThirdBodyEfficienciesV1::new(0.0, vec![]).unwrap()),
    ] {
        let m = mechanism(0.0, unrepresentable_high, law);
        let mut w = KineticWorkspaceV1::new(&m).unwrap();
        m.evaluate(Temperature::new(500.0), &[2.0, 0.0, 4.0], &mut w, || false)
            .unwrap();
        assert_eq!(w.reaction_rates()[0], 0.0);
        assert_eq!(w.rate_jacobian(), &[0.0, 0.0, 0.0]);
    }
}
