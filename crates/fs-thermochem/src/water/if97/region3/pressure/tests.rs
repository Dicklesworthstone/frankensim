use super::*;

fn controls() -> DensePressureSolve {
    DensePressureSolve {
        pressure_tolerance_pa: 1e-3,
        density_tolerance_kg_m3: 1e-8,
        maximum_iterations: 128,
    }
}
fn solve_at(p: f64, t: f64) -> DensePressureSolution {
    If97DenseState::at_pressure_temperature(
        Pressure::new(p),
        Temperature::new(t),
        controls(),
        || false,
    )
    .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        a.is_finite() && (a - b).abs() <= tolerance,
        "{a:.16e} versus {b:.16e}, residual {:.8e}, limit {tolerance:.8e}",
        a - b
    );
}

#[test]
fn table_33_pressure_driven_inversion_preserves_original_source() {
    // G2: original Table 33 printed p,T -> independent published rho,h,s.
    for (p, t, rho, h, s) in [
        (25.5837018e6, 650.0, 500.0, 1863430.19, 4054.27273),
        (22.2930643e6, 650.0, 200.0, 2375124.01, 4854.38792),
        (78.3095639e6, 750.0, 500.0, 2258688.45, 4469.71906),
    ] {
        let solved = solve_at(p, t);
        // Table p is rounded to 0.1 Pa. Propagate its half-last-digit
        // uncertainty through the independent analytic p derivatives before
        // comparing the printed calorics at the unrounded density state.
        let source = solved.state;
        let dv_dt = source.dv_dt_at_pressure().value();
        let dh_dp = source.specific_volume().value() - t * dv_dt;
        near(
            source.density().value(),
            rho,
            0.05 / source.dp_drho_at_temperature().value() + 2e-8,
        );
        near(
            source.specific_enthalpy().value(),
            h,
            0.005 + 0.05 * dh_dp.abs(),
        );
        near(
            source.specific_entropy().value(),
            s,
            0.000005 + 0.05 * dv_dt.abs(),
        );
        assert!(solved.pressure_residual_pa.abs() <= controls().pressure_tolerance_pa);
        assert!(
            solved.density_bracket_kg_m3[1] - solved.density_bracket_kg_m3[0]
                <= controls().density_tolerance_kg_m3
        );
        assert!(solved.iterations <= controls().maximum_iterations);
        assert!(solved.evaluations <= 3 * controls().maximum_iterations + 4);
        assert_eq!(
            solved.state,
            If97DenseState::at_density_temperature(solved.state.density(), Temperature::new(t))
                .unwrap()
        );
        assert_eq!(
            solved.pressure_residual_pa,
            solved.state.pressure().value() - p
        );
        near(
            solved.state.specific_enthalpy().value()
                - solved.state.specific_internal_energy().value(),
            solved.state.pressure().value() / solved.state.density().value(),
            6e-10,
        );
    }
}

#[test]
fn stable_branches_density_roundtrips_and_pressure_derivatives() {
    // G1/G3: independently perturb pressure; the inverse must reproduce the
    // analytic compressibility of the original Helmholtz evaluator.
    for (p, t) in [
        (20e6, 630.0),
        (17e6, 625.0),
        (20e6, 640.0),
        (21e6, 645.0),
        (25e6, 650.0),
        (80e6, 750.0),
        (95e6, 850.0),
    ] {
        let s = solve_at(p, t).state;
        let exact = solve_at(s.pressure().value(), t);
        near(exact.state.density().value(), s.density().value(), 3e-8);
        let step = 100.0;
        let drho = (solve_at(p + step, t).state.density().value()
            - solve_at(p - step, t).state.density().value())
            / (2.0 * step);
        near(drho, 1.0 / s.dp_drho_at_temperature().value(), 2e-9);
        let (lo, hi) = (s.density().value() - 1e-8, s.density().value() + 1e-8);
        let (stiffness, cv) = stability_bounds(lo, hi, t);
        assert!(stiffness > 0.0 && cv > 0.0);
        assert!(stiffness <= s.dp_drho_at_temperature().value());
        assert!(cv <= s.cv().value());
        let raw = Isotherm::new(t).evaluate(s.density().value());
        assert_eq!(raw.pressure, s.pressure().value());
        assert_eq!(raw.stiffness, s.dp_drho_at_temperature().value());
    }
}

#[test]
fn high_temperature_saturation_has_two_original_stable_endpoints() {
    // G0/G3: true separate phase densities, positive latent h, common original
    // R4 T, and explicit finite original R3/R4 Gibbs mismatch.
    for t in [624.0, 630.0, 640.0, 646.0, 647.0] {
        let pair =
            If97DenseSaturationPair::at_temperature(Temperature::new(t), controls(), || false)
                .unwrap();
        let p = pair.pressure().value();
        assert_eq!(pair.temperature().value(), t);
        assert_eq!(pair.saturation_pressure_residual_pa(), 0.0);
        assert!(pair.liquid().density().value() > CRITICAL_DENSITY_KG_M3);
        assert!(pair.vapor().density().value() < CRITICAL_DENSITY_KG_M3);
        assert!(pair.liquid().pressure().value() >= p);
        assert!(pair.vapor().pressure().value() <= p);
        assert!(pair.latent_enthalpy().value() > 0.0);
        assert!(pair.latent_entropy().value() > 0.0);
        assert!(pair.gibbs_energy_gap().value().is_finite());
        near(
            pair.latent_enthalpy().value() - t * pair.latent_entropy().value(),
            pair.gibbs_energy_gap().value(),
            1e-8,
        );
        for endpoint in [pair.liquid_solution(), pair.vapor_solution()] {
            assert!(endpoint.pressure_residual_pa.abs() <= controls().pressure_tolerance_pa);
            assert_eq!(
                endpoint.state,
                If97DenseState::at_density_temperature(
                    endpoint.state.density(),
                    endpoint.state.temperature()
                )
                .unwrap()
            );
        }
        let inverse =
            If97DenseSaturationPair::at_pressure(Pressure::new(p), controls(), || false).unwrap();
        assert_eq!(inverse.pressure().value(), p);
        near(inverse.temperature().value(), t, 1e-8);
        near(
            inverse.liquid().density().value(),
            pair.liquid().density().value(),
            2e-6,
        );
        near(
            inverse.vapor().density().value(),
            pair.vapor().density().value(),
            2e-6,
        );
        assert_eq!(
            inverse.saturation_pressure_residual_pa(),
            saturation_pressure(inverse.temperature()).unwrap().value() - p
        );
        assert_eq!(
            If97DenseState::at_pressure_temperature(
                Pressure::new(p),
                Temperature::new(t),
                controls(),
                || false
            ),
            Err(DensePressureError::AmbiguousSaturation)
        );
    }
    let pair =
        If97DenseSaturationPair::at_temperature(Temperature::new(640.0), controls(), || false)
            .unwrap();
    assert!(
        pair.gibbs_energy_gap().value().abs() > 1e-6,
        "do not silently repair the original regional Gibbs gap"
    );
}

#[test]
fn requested_pressure_and_closed_source_boundaries_remain_explicit() {
    let pair =
        If97DenseSaturationPair::at_pressure(Pressure::new(20e6), controls(), || false).unwrap();
    assert_eq!(pair.pressure().value(), 20e6);
    for endpoint in [pair.liquid_solution(), pair.vapor_solution()] {
        assert_eq!(endpoint.requested_pressure_pa, 20e6);
        assert_eq!(
            endpoint.pressure_residual_pa,
            endpoint.state.pressure().value() - 20e6
        );
    }
    for t in [623.15, 650.0, 750.0, 850.0] {
        let s = solve_at(100e6, t);
        assert!(s.state.pressure().value() <= 100e6);
        assert!(s.pressure_residual_pa.abs() <= controls().pressure_tolerance_pa);
    }
    for t in [650.0, 700.0, 850.0] {
        let p = boundary23_pressure(Temperature::new(t)).unwrap().value();
        let s = solve_at(p, t);
        assert!(s.state.pressure().value() >= p);
        assert!(s.pressure_residual_pa.abs() <= controls().pressure_tolerance_pa);
    }
}

#[test]
fn domain_work_resolution_and_cancellation_refuse_atomically() {
    for (p, t) in [
        (0.0, 650.0),
        (f64::NAN, 650.0),
        (20e6, f64::INFINITY),
        (20e6, 600.0),
        (1e6, 650.0),
        (100e6 + 1.0, 650.0),
        (100e6, 900.0),
    ] {
        assert!(
            If97DenseState::at_pressure_temperature(
                Pressure::new(p),
                Temperature::new(t),
                controls(),
                || false
            )
            .is_err()
        );
    }
    assert_eq!(
        If97DenseState::at_pressure_temperature(
            Pressure::new(CRITICAL_PRESSURE_PA),
            Temperature::new(CRITICAL_TEMPERATURE_K),
            controls(),
            || false
        ),
        Err(DensePressureError::CriticalSingularity)
    );
    for t in [623.15, 647.096, 650.0] {
        assert!(
            If97DenseSaturationPair::at_temperature(Temperature::new(t), controls(), || false)
                .is_err()
        );
    }
    for bad in [
        DensePressureSolve {
            maximum_iterations: 0,
            ..controls()
        },
        DensePressureSolve {
            maximum_iterations: 257,
            ..controls()
        },
        DensePressureSolve {
            pressure_tolerance_pa: -1.0,
            ..controls()
        },
        DensePressureSolve {
            density_tolerance_kg_m3: 0.0,
            ..controls()
        },
    ] {
        assert!(matches!(
            If97DenseState::at_pressure_temperature(
                Pressure::new(25e6),
                Temperature::new(650.0),
                bad,
                || false
            ),
            Err(DensePressureError::Invalid(_))
        ));
    }
    assert!(matches!(
        If97DenseState::at_pressure_temperature(
            Pressure::new(25e6),
            Temperature::new(650.0),
            DensePressureSolve {
                maximum_iterations: 1,
                ..controls()
            },
            || false
        ),
        Err(DensePressureError::NotConverged { .. })
    ));
    let unresolved = DensePressureSolve {
        density_tolerance_kg_m3: f64::MIN_POSITIVE,
        ..controls()
    };
    assert!(
        If97DenseState::at_pressure_temperature(
            Pressure::new(25e6),
            Temperature::new(650.0),
            unresolved,
            || false
        )
        .is_err()
    );
    assert_eq!(
        If97DenseState::at_pressure_temperature(
            Pressure::new(25e6),
            Temperature::new(650.0),
            controls(),
            || true
        ),
        Err(DensePressureError::Cancelled)
    );
    let mut polls = 0;
    assert_eq!(
        If97DenseSaturationPair::at_pressure(Pressure::new(20e6), controls(), || {
            polls += 1;
            polls == 8
        }),
        Err(DensePressureError::Cancelled)
    );
    assert_eq!(polls, 8);
    let complete = solve_at(25e6, 650.0);
    let mut calls = 0;
    let mut total_calls = 0;
    If97DenseState::at_pressure_temperature(
        Pressure::new(25e6),
        Temperature::new(650.0),
        controls(),
        || {
            total_calls += 1;
            false
        },
    )
    .unwrap();
    assert_eq!(
        If97DenseState::at_pressure_temperature(
            Pressure::new(25e6),
            Temperature::new(650.0),
            controls(),
            || {
                calls += 1;
                calls == total_calls
            }
        ),
        Err(DensePressureError::Cancelled)
    );
    assert_eq!(solve_at(25e6, 650.0), complete);
}
