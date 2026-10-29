use super::*;

fn state(rho: f64, t: f64) -> If97DenseState {
    If97DenseState::at_density_temperature(Density::new(rho), Temperature::new(t)).unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance,
        "{actual:.15e} differs from {expected:.15e} by {:.8e}, limit {tolerance:.8e}",
        actual - expected,
    );
}
fn derivative(f: impl Fn(f64) -> f64, x: f64, step: f64) -> f64 {
    (f(x - 2.0 * step) - 8.0 * f(x - step) + 8.0 * f(x + step) - f(x + 2.0 * step)) / (12.0 * step)
}

#[test]
fn original_table_33_all_three_independent_reference_states() {
    // Original published values and their printed final-digit precision.
    // SI conversion changes no coefficient or caloric reference.
    let rows = [
        (
            500.0,
            650.0,
            [
                25.5837018e6,
                1863.43019e3,
                1812.26279e3,
                4054.27273,
                13893.5717,
                502.005554,
            ],
        ),
        (
            200.0,
            650.0,
            [
                22.2930643e6,
                2375.12401e3,
                2263.65868e3,
                4854.38792,
                44657.9342,
                383.444594,
            ],
        ),
        (
            500.0,
            750.0,
            [
                78.3095639e6,
                2258.68845e3,
                2102.06932e3,
                4469.71906,
                6341.65359,
                760.696041,
            ],
        ),
    ];
    for (rho, t, expected) in rows {
        let source = state(rho, t);
        let actual = [
            source.pressure().value(),
            source.specific_enthalpy().value(),
            source.specific_internal_energy().value(),
            source.specific_entropy().value(),
            source.cp().value(),
            source.sound_speed().value(),
        ];
        for ((actual, expected), tolerance) in actual
            .into_iter()
            .zip(expected)
            .zip([0.05, 0.005, 0.005, 0.000005, 0.00005, 0.0000005])
        {
            near(actual, expected, tolerance);
        }
        assert_eq!(source.density().value(), rho);
        assert_eq!(source.temperature().value(), t);
        assert_eq!(source.phase(), If97DensePhase::Supercritical);
        near(source.specific_volume().value() * rho, 1.0, 2e-16);
        near(
            source.specific_enthalpy().value() - source.specific_internal_energy().value(),
            source.pressure().value() / rho,
            6e-10,
        );
        near(
            source.specific_gibbs_energy().value(),
            source.specific_enthalpy().value() - t * source.specific_entropy().value(),
            8e-10,
        );
        near(
            source.specific_helmholtz_energy().value(),
            source.specific_internal_energy().value() - t * source.specific_entropy().value(),
            8e-10,
        );
    }
}

#[test]
fn analytic_density_temperature_calorics_match_independent_differences() {
    for (rho, t) in [
        (500.0, 650.0),
        (200.0, 650.0),
        (500.0, 750.0),
        (700.0, 630.0),
    ] {
        let source = state(rho, t);
        near(
            source.dp_drho_at_temperature().value(),
            derivative(|r| state(r, t).pressure().value(), rho, 0.01),
            0.02,
        );
        near(
            source.dp_dt_at_density().value(),
            derivative(
                |temperature| state(rho, temperature).pressure().value(),
                t,
                0.01,
            ),
            0.02,
        );
        near(
            source.cv().value(),
            derivative(
                |temperature| state(rho, temperature).specific_internal_energy().value(),
                t,
                0.01,
            ),
            0.0001,
        );
        near(
            source.dh_dt_at_density().value(),
            derivative(
                |temperature| state(rho, temperature).specific_enthalpy().value(),
                t,
                0.01,
            ),
            0.0001,
        );
        near(
            source.du_drho_at_temperature().value(),
            derivative(
                |r| state(r, t).specific_internal_energy().value(),
                rho,
                0.01,
            ),
            0.0001,
        );
        near(
            -source.dp_dt_at_density().value() / (rho * rho),
            derivative(|r| state(r, t).specific_entropy().value(), rho, 0.01),
            0.000001,
        );
        near(
            source.cv().value() / t,
            derivative(
                |temperature| state(rho, temperature).specific_entropy().value(),
                t,
                0.01,
            ),
            0.000001,
        );
        assert!(source.cp().value() > source.cv().value());
        assert!(source.dv_dp_at_temperature().value() < 0.0);
        near(
            source.thermal_expansion().value(),
            source.dv_dt_at_pressure().value() * rho,
            1e-16,
        );
        near(
            source.isothermal_compressibility().value(),
            -source.dv_dp_at_temperature().value() * rho,
            1e-20,
        );
    }
}

#[test]
fn constant_pressure_cp_and_isentropic_sound_use_independent_paths() {
    for (rho, t) in [(500.0, 650.0), (200.0, 650.0), (500.0, 750.0)] {
        let initial = state(rho, t);
        let h_at_pressure = |temperature| {
            let mut density = rho;
            for _ in 0..8 {
                let current = state(density, temperature);
                density -= (current.pressure().value() - initial.pressure().value())
                    / current.dp_drho_at_temperature().value();
            }
            let answer = state(density, temperature);
            near(
                answer.pressure().value(),
                initial.pressure().value(),
                0.00002,
            );
            answer.specific_enthalpy().value()
        };
        near(
            initial.cp().value(),
            derivative(h_at_pressure, t, 0.002),
            0.002,
        );
        let p_at_entropy = |density| {
            let mut temperature = t;
            for _ in 0..8 {
                let current = state(density, temperature);
                temperature -= (current.specific_entropy().value()
                    - initial.specific_entropy().value())
                    / (current.cv().value() / temperature);
            }
            let answer = state(density, temperature);
            near(
                answer.specific_entropy().value(),
                initial.specific_entropy().value(),
                2e-9,
            );
            answer.pressure().value()
        };
        near(
            initial.sound_speed().value().powi(2),
            derivative(p_at_entropy, rho, 0.01),
            0.02,
        );
    }
}

#[test]
fn stable_source_domain_refuses_coexistence_metastability_and_extrapolated_roots() {
    assert_eq!(state(600.0, 640.0).phase(), If97DensePhase::Liquid);
    assert_eq!(state(150.0, 640.0).phase(), If97DensePhase::Vapor);
    assert_eq!(
        state(500.0, CRITICAL_TEMPERATURE_K).phase(),
        If97DensePhase::CriticalIsotherm
    );
    for (rho, t) in [
        (f64::NAN, 650.0),
        (0.0, 650.0),
        (-1.0, 650.0),
        (500.0, f64::INFINITY),
    ] {
        assert!(matches!(
            If97DenseState::at_density_temperature(Density::new(rho), Temperature::new(t)),
            Err(If97DenseError::InvalidInput(_))
        ));
    }
    for (rho, t) in [
        (500.0, MINIMUM_TEMPERATURE_K.next_down()),
        (500.0, MAXIMUM_TEMPERATURE_K.next_up()),
        (DENSITY_CEILING_KG_M3.next_up(), 650.0),
        (10.0, 650.0),
        (800.0, 650.0),
        (2000.0, 650.0),
    ] {
        assert!(matches!(
            If97DenseState::at_density_temperature(Density::new(rho), Temperature::new(t)),
            Err(If97DenseError::OutsideRegion { .. })
        ));
    }
    assert_eq!(
        If97DenseState::at_density_temperature(
            Density::new(CRITICAL_DENSITY_KG_M3),
            Temperature::new(CRITICAL_TEMPERATURE_K),
        ),
        Err(If97DenseError::CriticalSingularity),
    );
    assert!(matches!(
        If97DenseState::at_density_temperature(Density::new(322.0), Temperature::new(640.0)),
        Err(If97DenseError::Unstable { .. })
    ));
    for rho in [470.0, 190.0] {
        assert!(matches!(
            If97DenseState::at_density_temperature(Density::new(rho), Temperature::new(640.0)),
            Err(If97DenseError::UnsupportedPhase { .. })
        ));
    }
}

// A small independent outward-rounded test enclosure, rather than sampling
// temperatures to assert that the upper density guard excludes all real states.
#[derive(Clone, Copy)]
struct Interval(f64, f64);
impl Interval {
    fn point(x: f64) -> Self {
        Self(x, x)
    }
    fn add(self, other: Self) -> Self {
        Self((self.0 + other.0).next_down(), (self.1 + other.1).next_up())
    }
    fn mul(self, other: Self) -> Self {
        let values = [
            self.0 * other.0,
            self.0 * other.1,
            self.1 * other.0,
            self.1 * other.1,
        ];
        Self(
            values.into_iter().fold(f64::INFINITY, f64::min).next_down(),
            values
                .into_iter()
                .fold(f64::NEG_INFINITY, f64::max)
                .next_up(),
        )
    }
    fn reciprocal_positive(self) -> Self {
        assert!(self.0 > 0.0);
        Self((1.0 / self.1).next_down(), (1.0 / self.0).next_up())
    }
}

#[test]
fn upper_density_envelope_is_beyond_100_mpa_over_the_entire_source_interval() {
    // Density is fixed. Collect equal tau powers before interval evaluation;
    // otherwise cancelling density terms needlessly widen the temperature box.
    let mut d = [Interval::point(1.0); 12];
    for i in 1..d.len() {
        d[i] = d[i - 1].mul(Interval::point(2.5));
    }
    let mut coefficients = [Interval::point(0.0); 27];
    coefficients[0] = Interval::point(LOG_COEFFICIENT);
    for &(i, j, n) in &TERMS {
        if i != 0 {
            coefficients[j] =
                coefficients[j].add(Interval::point(n).mul(d[i]).mul(Interval::point(i as f64)));
        }
    }
    let mut pending = vec![(MINIMUM_TEMPERATURE_K, MAXIMUM_TEMPERATURE_K, 0)];
    let mut intervals = 0;
    while let Some((lower, upper, depth)) = pending.pop() {
        intervals += 1;
        assert!(
            intervals <= 32767,
            "density-envelope interval proof exhausted"
        );
        let temperature = Interval(lower, upper);
        let tau = Interval::point(CRITICAL_TEMPERATURE_K).mul(temperature.reciprocal_positive());
        let mut pressure_term = coefficients[26];
        for coefficient in coefficients[..26].iter().rev() {
            pressure_term = pressure_term.mul(tau).add(*coefficient);
        }
        let pressure = Interval::point(DENSITY_CEILING_KG_M3)
            .mul(Interval::point(R_SPECIFIC_J_KG_K))
            .mul(temperature)
            .mul(pressure_term);
        if pressure.0 > MAXIMUM_PRESSURE_PA {
            continue;
        }
        assert!(
            depth < 14,
            "upper-density pressure bound [{},{}] at [{lower},{upper}]",
            pressure.0,
            pressure.1
        );
        let midpoint = lower + 0.5 * (upper - lower);
        pending.push((midpoint, upper, depth + 1));
        pending.push((lower, midpoint, depth + 1));
    }
}
