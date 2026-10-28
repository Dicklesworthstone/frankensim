use super::*;

fn state(t: f64) -> SaturatedWaterState {
    SaturatedWaterState::at_temperature(Temperature::new(t)).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.17e}, expected {expected:.17e}, tolerance {tolerance:.3e}"
    );
}

fn solve() -> WaterSaturationSolve {
    WaterSaturationSolve {
        pressure_tolerance_pa: 1e-6,
        maximum_iterations: 96,
    }
}

#[test]
fn released_table_one_values_use_their_published_rounding_precision() {
    // IAPWS SR1-86(1992), Table 1. Columns: p, dp/dT, rho_l, rho_v,
    // alpha, h_l, h_v, phi, s_l, s_v, in coherent SI. Tolerances are half
    // of the printed last digit, except the zero entropy reference, where
    // the last published digit of d_phi sets its finite coefficient precision.
    let rows = [
        (
            273.16,
            [
                611.657, 44.436693, 999.789, 0.00485426, -11.529101, 0.611786, 2500.5e3, -0.04,
                0.0, 9.154e3,
            ],
            [
                0.0005,
                0.0000005,
                0.0005,
                0.000000005,
                0.0000005,
                0.0000005,
                50.0,
                0.005,
                0.00008,
                0.5,
            ],
        ),
        (
            373.1243,
            [
                0.101325e6, 3.616e3, 958.365, 0.597586, 417.65e3, 419.05e3, 2675.7e3, 1.303e3,
                1.307e3, 7.355e3,
            ],
            [0.5, 0.5, 0.0005, 0.0000005, 5.0, 5.0, 50.0, 0.5, 0.5, 0.5],
        ),
        (
            647.096,
            [
                22.064e6, 268e3, 322.0, 322.0, 1548e3, 2086.6e3, 2086.6e3, 3.578e3, 4.410e3,
                4.410e3,
            ],
            [0.0, 500.0, 0.0, 0.0, 500.0, 50.0, 50.0, 0.5, 0.5, 0.5],
        ),
    ];
    for (t, expected, tolerances) in rows {
        let pair = state(t);
        let slope = pair.pressure_derivative().value();
        let rho = pair.liquid_density().value();
        let h = pair.liquid_specific_enthalpy().value();
        let s = pair.liquid_specific_entropy().value();
        let actual = [
            pair.pressure().value(),
            slope,
            rho,
            pair.vapor_density().value(),
            h - t * slope / rho,
            h,
            pair.vapor_specific_enthalpy().value(),
            s - slope / rho,
            s,
            pair.vapor_specific_entropy().value(),
        ];
        for ((actual, expected), tolerance) in actual.into_iter().zip(expected).zip(tolerances) {
            close(actual, expected, tolerance);
        }
    }
    // No post-evaluation datum shift is applied to make rounded coefficients
    // reproduce the convention as bitwise zero.
    let triple = state(TRIPLE_POINT_TEMPERATURE_K);
    close(triple.liquid_specific_internal_energy().value(), 0.0, 5e-7);
    assert_ne!(triple.liquid_specific_internal_energy().value(), 0.0);
    assert_ne!(triple.liquid_specific_entropy().value(), 0.0);
}

#[test]
fn analytic_pressure_slope_matches_independent_central_difference() {
    for t in [280.0, 350.0, 500.0, 620.0] {
        let step = 0.01;
        let p = |at| state(at).pressure().value();
        let numerical = (p(t - 2.0 * step) - 8.0 * p(t - step) + 8.0 * p(t + step)
            - p(t + 2.0 * step))
            / (12.0 * step);
        let analytic = state(t).pressure_derivative().value();
        close(numerical, analytic, 2e-9 * analytic);
    }
}

#[test]
fn latent_heat_obeys_clapeyron_entropy_and_pressure_volume_work() {
    for t in [
        273.16,
        300.0,
        373.1243,
        450.0,
        600.0,
        CRITICAL_TEMPERATURE_K.next_down(),
    ] {
        let pair = state(t);
        let p = pair.pressure().value();
        let vl = 1.0 / pair.liquid_density().value();
        let vv = 1.0 / pair.vapor_density().value();
        let dv = pair.specific_volume_gap().value();
        let latent = pair.latent_enthalpy().value();
        let entropy = pair.latent_entropy().value();
        let work = p * dv;
        close(dv, vv - vl, 16.0 * f64::EPSILON * (vv + vl));
        close(
            latent,
            t * pair.pressure_derivative().value() * (vv - vl),
            32.0 * f64::EPSILON * t * pair.pressure_derivative().value() * (vv + vl),
        );
        close(latent, t * entropy, 2.0 * f64::EPSILON * latent);
        close(
            pair.latent_internal_energy().value() + work,
            latent,
            8.0 * f64::EPSILON * (latent + work),
        );
        for (h, u, v) in [
            (
                pair.liquid_specific_enthalpy().value(),
                pair.liquid_specific_internal_energy().value(),
                vl,
            ),
            (
                pair.vapor_specific_enthalpy().value(),
                pair.vapor_specific_internal_energy().value(),
                vv,
            ),
        ] {
            close(
                u + p * v,
                h,
                16.0 * f64::EPSILON * (h.abs() + u.abs() + p * v),
            );
        }
        let gl =
            pair.liquid_specific_enthalpy().value() - t * pair.liquid_specific_entropy().value();
        let gv = pair.vapor_specific_enthalpy().value() - t * pair.vapor_specific_entropy().value();
        close(
            gl,
            gv,
            16.0 * f64::EPSILON
                * (pair.vapor_specific_enthalpy().value()
                    + t * pair.vapor_specific_entropy().value()),
        );
    }
}

#[test]
fn critical_approach_retains_positive_phase_gap_until_exact_coalescence() {
    let mut last_gap = f64::INFINITY;
    for t in [
        600.0,
        CRITICAL_TEMPERATURE_K - 1.0,
        CRITICAL_TEMPERATURE_K - 0.001,
        CRITICAL_TEMPERATURE_K - 1e-6,
        CRITICAL_TEMPERATURE_K.next_down(),
    ] {
        let pair = state(t);
        assert!(pair.liquid_density().value() > pair.vapor_density().value());
        assert!(pair.latent_enthalpy().value() > 0.0);
        assert!(pair.latent_entropy().value() > 0.0);
        assert!(pair.latent_internal_energy().value() > 0.0);
        assert!(pair.specific_volume_gap().value() > 0.0);
        assert!(pair.latent_enthalpy().value() < last_gap);
        last_gap = pair.latent_enthalpy().value();
    }
    let critical = state(CRITICAL_TEMPERATURE_K);
    assert_eq!(critical.pressure().value(), CRITICAL_PRESSURE_PA);
    assert_eq!(critical.liquid_density().value(), CRITICAL_DENSITY_KG_M3);
    assert_eq!(critical.vapor_density(), critical.liquid_density());
    assert_eq!(
        critical.liquid_specific_enthalpy(),
        critical.vapor_specific_enthalpy()
    );
    assert_eq!(
        critical.liquid_specific_entropy(),
        critical.vapor_specific_entropy()
    );
    assert_eq!(
        critical.liquid_specific_internal_energy(),
        critical.vapor_specific_internal_energy()
    );
    assert_eq!(critical.specific_volume_gap().value(), 0.0);
    assert_eq!(critical.latent_enthalpy().value(), 0.0);
    assert_eq!(critical.latent_entropy().value(), 0.0);
    assert_eq!(critical.latent_internal_energy().value(), 0.0);
}

#[test]
fn original_pressure_inverse_roundtrips_endpoints_and_interior_states() {
    let controls = solve();
    let mut previous_pressure = 0.0;
    for t in [
        TRIPLE_POINT_TEMPERATURE_K,
        TRIPLE_POINT_TEMPERATURE_K.next_up(),
        280.0,
        373.1243,
        450.0,
        600.0,
        CRITICAL_TEMPERATURE_K.next_down(),
        CRITICAL_TEMPERATURE_K,
    ] {
        let original = state(t);
        assert!(original.pressure().value() > previous_pressure);
        previous_pressure = original.pressure().value();
        let solved =
            SaturatedWaterState::at_pressure(original.pressure(), controls, || false).unwrap();
        close(
            solved.state().pressure().value(),
            original.pressure().value(),
            controls.pressure_tolerance_pa,
        );
        close(
            solved.state().temperature().value(),
            t,
            4.0 * controls.pressure_tolerance_pa / original.pressure_derivative().value()
                + 2.0 * (t.next_up() - t),
        );
        assert_eq!(
            solved.pressure_residual().value(),
            solved.state().pressure().value() - original.pressure().value()
        );
        assert_eq!(
            solved.into_state(),
            state(solved.state().temperature().value())
        );
        if t == TRIPLE_POINT_TEMPERATURE_K || t == CRITICAL_TEMPERATURE_K {
            assert_eq!(solved.iterations(), 0);
            assert_eq!(solved.state(), &original);
        }
    }
    let normal_boiling =
        SaturatedWaterState::at_pressure(Pressure::new(101325.0), controls, || false).unwrap();
    close(
        normal_boiling.state().temperature().value(),
        373.1243,
        0.00005,
    );
    // A rounded table temperature is not used as an exact inverse result.
    assert_ne!(normal_boiling.state().temperature().value(), 373.1243);
    assert_eq!(
        normal_boiling,
        SaturatedWaterState::at_pressure(Pressure::new(101325.0), controls, || false).unwrap()
    );
}

#[test]
fn source_domain_and_explicit_work_budgets_are_enforced() {
    for t in [
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        0.0,
        TRIPLE_POINT_TEMPERATURE_K.next_down(),
        CRITICAL_TEMPERATURE_K.next_up(),
    ] {
        assert!(matches!(
            SaturatedWaterState::at_temperature(Temperature::new(t)),
            Err(WaterSaturationError::TemperatureOutOfRange { .. })
        ));
    }
    let minimum = state(TRIPLE_POINT_TEMPERATURE_K).pressure().value();
    for p in [
        f64::NAN,
        f64::INFINITY,
        0.0,
        611.657,
        minimum.next_down(),
        CRITICAL_PRESSURE_PA.next_up(),
    ] {
        assert!(matches!(
            SaturatedWaterState::at_pressure(Pressure::new(p), solve(), || false),
            Err(WaterSaturationError::PressureOutOfRange { .. })
        ));
    }
    for tolerance in [f64::NAN, f64::INFINITY, 0.0, -1.0] {
        let controls = WaterSaturationSolve {
            pressure_tolerance_pa: tolerance,
            ..solve()
        };
        assert!(matches!(
            SaturatedWaterState::at_pressure(Pressure::new(1e6), controls, || false),
            Err(WaterSaturationError::InvalidSolve(_))
        ));
    }
    for maximum_iterations in [0, MAX_PRESSURE_ITERATIONS + 1] {
        let controls = WaterSaturationSolve {
            maximum_iterations,
            ..solve()
        };
        assert!(matches!(
            SaturatedWaterState::at_pressure(Pressure::new(1e6), controls, || false),
            Err(WaterSaturationError::InvalidSolve(_))
        ));
    }
    let short = WaterSaturationSolve {
        maximum_iterations: 1,
        ..solve()
    };
    assert!(matches!(
        SaturatedWaterState::at_pressure(Pressure::new(1e6), short, || false),
        Err(WaterSaturationError::NotConverged { iterations: 1, .. })
    ));
}

#[test]
fn inverse_polls_cancellation_before_work_during_solve_and_before_acceptance() {
    let pressure = Pressure::new(1e6);
    let mut polls = 0;
    let solution = SaturatedWaterState::at_pressure(pressure, solve(), || {
        polls += 1;
        false
    })
    .unwrap();
    assert_eq!(polls, solution.iterations() + 2);
    for stop in [1, 3, polls] {
        let mut count = 0;
        let error = SaturatedWaterState::at_pressure(pressure, solve(), || {
            count += 1;
            count == stop
        });
        assert_eq!(error, Err(WaterSaturationError::Cancelled));
        assert_eq!(count, stop);
    }
    let mut polls = 0;
    assert_eq!(
        SaturatedWaterState::at_pressure(Pressure::new(CRITICAL_PRESSURE_PA), solve(), || {
            polls += 1;
            polls == 2
        }),
        Err(WaterSaturationError::Cancelled)
    );
}
