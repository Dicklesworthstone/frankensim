//! Shifted power-law weakening and the measured rosin fit of Woodhouse &
//! Galluzzo (2025), "Enhanced tribological modelling of violin rosin",
//! Tribology Letters, doi:10.1007/s11249-025-02062-4, eq. (13).
use fs_tribo::FrictionLaw;

/// The printed equation, evaluated with std (not the law's det::pow route).
fn printed_eq13(speed_m_s: f64) -> f64 {
    0.2 * (speed_m_s + 0.011).powf(-0.4)
}

fn power(static_mu: f64, kinetic_mu: f64, characteristic_speed: f64, exponent: f64) -> FrictionLaw {
    FrictionLaw::ShiftedPower { static_mu, kinetic_mu, characteristic_speed, exponent }
}

#[test]
fn eq13_matches_the_printed_formula_at_three_speeds() {
    let law = FrictionLaw::rosin_woodhouse_galluzzo_2025_eq13();
    // 0.01, 0.1 and 0.3 m/s: inside the measured range (about 0.3 m/s).
    // The second column is the printed formula to 5 digits, a units check.
    for (speed, rounded) in [(0.01, 0.93787), (0.1, 0.48184), (0.3, 0.31910)] {
        let mu = law.kinetic_coefficient(speed).unwrap();
        let printed = printed_eq13(speed);
        assert!((mu - printed).abs() <= 1.0e-13 * printed, "{speed} m/s: {mu} vs {printed}");
        assert!((mu - rounded).abs() < 5.0e-6, "{speed} m/s: {mu} vs {rounded}");
    }
    // Continuous at rest: static capacity is the curve's zero-speed limit.
    let rest = printed_eq13(0.0);
    let static_mu = law.static_coefficient().unwrap();
    assert!((static_mu - rest).abs() <= 1.0e-13 * rest, "{static_mu} vs {rest}");
    assert!((static_mu - 1.21471).abs() < 5.0e-6, "{static_mu}");
    assert_eq!(law.zero_speed_kinetic_coefficient().unwrap(), static_mu);
    // Weakens monotonically towards zero, never below it.
    let mut previous = static_mu;
    for speed in [1.0e-6, 1.0e-3, 0.011, 0.05, 0.3, 1.0, 10.0, 1.0e4] {
        let mu = law.kinetic_coefficient(speed).unwrap();
        assert!(mu > 0.0 && mu < previous, "{speed}: {mu} vs {previous}");
        previous = mu;
    }
}

#[test]
fn unit_exponent_is_the_hyperbolic_curve() {
    let hyperbolic = FrictionLaw::Hyperbolic { static_mu: 0.8, kinetic_mu: 0.3, characteristic_speed: 0.04 };
    let unit = power(0.8, 0.3, 0.04, 1.0);
    for speed in [1.0e-9, 1.0e-3, 0.04, 0.2, 3.0, 1.0e6] {
        let a = hyperbolic.kinetic_coefficient(speed).unwrap();
        let b = unit.kinetic_coefficient(speed).unwrap();
        assert!((a - b).abs() <= 4.0 * f64::EPSILON * a, "{speed}: {a} vs {b}");
    }
    assert_eq!(unit.static_coefficient().unwrap(), 0.8);
}

#[test]
fn malformed_shifted_power_laws_are_refused_with_text() {
    for (law, field) in [
        (power(0.8, 0.3, 0.04, 0.0), "exponent"),
        (power(0.8, 0.3, 0.04, -0.4), "exponent"),
        (power(0.8, 0.3, 0.04, f64::NAN), "exponent"),
        (power(0.8, 0.3, 0.04, f64::INFINITY), "exponent"),
        (power(0.8, 0.3, 0.0, 0.4), "characteristic_speed"),
        (power(0.8, 0.3, -0.011, 0.4), "characteristic_speed"),
        (power(0.8, 0.9, 0.04, 0.4), "kinetic_mu > static_mu"),
        (power(-0.8, 0.3, 0.04, 0.4), "static_mu"),
        (power(0.8, f64::NAN, 0.04, 0.4), "kinetic_mu"),
    ] {
        let error = law.validate().unwrap_err();
        assert!(error.to_string().contains(field), "{law:?}: {error}");
        assert!(law.kinetic_coefficient(0.1).is_err(), "{law:?}");
        assert!(law.static_coefficient().is_err(), "{law:?}");
    }
    let law = FrictionLaw::rosin_woodhouse_galluzzo_2025_eq13();
    assert!(law.kinetic_coefficient(0.0).is_err());
    assert!(law.kinetic_coefficient(-0.1).is_err());
}
