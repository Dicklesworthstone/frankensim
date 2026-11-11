use super::*;
use crate::sublimation::molecular::momentum::flux_weighted_normal_speed;

const NORMAL: [f64; 3] = [1.0, 0.0, 0.0];

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs(),
        "{actual:e} != {expected:e} (relative tolerance {tolerance:e})"
    );
}

// Independent Simpson integration of the ORIGINAL incoming Gaussian, without
// erf, moment recurrences, continued fractions or the production energy chart.
fn quadrature(a: f64) -> [f64; 3] {
    let end = if a < 0.0 {
        (40.0 / (1.0 - a)).min(12.0)
    } else {
        a + 12.0
    };
    let count = 32_768;
    let step = end / count as f64;
    let mut moments = [0.0; 3];
    for i in 0..=count {
        let y = step * i as f64;
        let weight = if i == 0 || i == count {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        let density = (-0.5 * (y - a) * (y - a)).exp() / (2.0 * core::f64::consts::PI).sqrt();
        let mut integrand = weight * y * density;
        for value in &mut moments {
            *value += integrand;
            integrand *= y;
        }
    }
    moments.map(|value| value * step / 3.0)
}

#[test]
fn drifting_mass_vector_and_crossing_energy_match_independent_integrals() {
    let (rho, r, temperature): (f64, f64, f64) = (0.73, 287.0, 340.0);
    let sigma = (r * temperature).sqrt();
    for a in [
        -8.0, -2.000_001, -2.0, -1.999_999, -1.0, 0.0, 0.5, 2.0, 8.0, 15.999_999, 16.0, 16.000_001,
    ] {
        let [m1, m2, m3] = quadrature(a);
        let actual = incident_maxwellian(
            rho,
            r,
            temperature,
            [-a * sigma, 0.7 * sigma, -0.3 * sigma],
            NORMAL,
        )
        .unwrap();
        relative(actual.incoming_mass_flux_kg_m2_s, rho * sigma * m1, 2e-11);
        if a == 0.0 {
            assert_eq!(
                actual.incoming_mass_flux_change_from_zero_drift_kg_m2_s,
                0.0
            );
        } else {
            relative(
                actual.incoming_mass_flux_change_from_zero_drift_kg_m2_s,
                rho * sigma * (m1 - 1.0 / (2.0 * core::f64::consts::PI).sqrt()),
                3e-11,
            );
        }
        let conditional = actual.moments.unwrap();
        relative(
            conditional.mean_relative_velocity_m_s[0],
            -sigma * m2 / m1,
            2e-11,
        );
        relative(
            conditional.mean_relative_velocity_m_s[1],
            0.7 * sigma,
            3e-15,
        );
        relative(
            conditional.mean_relative_velocity_m_s[2],
            -0.3 * sigma,
            3e-15,
        );
        // Two tangential thermal modes plus their drift energy, and the
        // independently integrated normal kinetic energy. Subtract 5RT/2
        // because the caller supplies the original species enthalpy.
        let kinetic = r * temperature * (1.0 + 0.5 * (0.7 * 0.7 + 0.3 * 0.3) + 0.5 * m3 / m1);
        relative(
            conditional.crossing_energy_minus_enthalpy_j_kg,
            kinetic - 2.5 * r * temperature,
            3e-10,
        );
    }
}

#[test]
fn tangent_drift_preserves_original_impingement_and_adds_its_crossing_energy() {
    let (rho, r, temperature): (f64, f64, f64) = (0.02, 461.5, 280.0);
    let beam = incident_maxwellian(rho, r, temperature, [0.0, 40.0, -15.0], NORMAL).unwrap();
    let moments = beam.moments.unwrap();
    relative(
        beam.incoming_mass_flux_kg_m2_s,
        rho * (r * temperature / (2.0 * core::f64::consts::PI)).sqrt(),
        4e-15,
    );
    assert_eq!(beam.incoming_mass_flux_change_from_zero_drift_kg_m2_s, 0.0);
    relative(
        moments.mean_relative_velocity_m_s[0],
        -flux_weighted_normal_speed(r, temperature).unwrap(),
        4e-15,
    );
    assert_eq!(moments.mean_relative_velocity_m_s[1..], [40.0, -15.0]);
    relative(
        moments.crossing_energy_minus_enthalpy_j_kg,
        0.5 * (40.0 * 40.0 + 15.0 * 15.0 - r * temperature),
        4e-15,
    );
    // Incident normal pressure remains p/2 despite nonzero tangential momentum.
    relative(
        -beam.incoming_mass_flux_kg_m2_s * moments.mean_relative_velocity_m_s[0],
        0.5 * rho * r * temperature,
        4e-15,
    );
}

#[test]
fn opposite_half_spaces_reconstruct_full_mass_momentum_and_energy_fluxes() {
    let (rho, r, temperature, h): (f64, f64, f64, f64) = (1.7, 2.0, 3.0, 37.0);
    for normal_velocity in [-5.0, -1.0, 0.0, 1.0, 5.0] {
        let velocity = [normal_velocity, 0.4, -0.7];
        let incoming = incident_maxwellian(rho, r, temperature, velocity, NORMAL).unwrap();
        let outgoing =
            incident_maxwellian(rho, r, temperature, velocity, [-1.0, 0.0, 0.0]).unwrap();
        let im = incoming.moments.unwrap();
        let om = outgoing.moments.unwrap();
        let ji = incoming.incoming_mass_flux_kg_m2_s;
        let jo = outgoing.incoming_mass_flux_kg_m2_s;
        let expected_mass = rho * normal_velocity;
        if normal_velocity == 0.0 {
            assert_eq!(jo - ji, 0.0);
        } else {
            relative(jo - ji, expected_mass, 3e-13);
        }
        for i in 0..3 {
            let actual =
                jo * om.mean_relative_velocity_m_s[i] - ji * im.mean_relative_velocity_m_s[i];
            let expected = rho * normal_velocity * velocity[i] + rho * r * temperature * NORMAL[i];
            if expected == 0.0 {
                assert_eq!(actual, 0.0);
            } else {
                relative(actual, expected, 3e-13);
            }
        }
        let actual_energy = jo * (h + om.crossing_energy_minus_enthalpy_j_kg)
            - ji * (h + im.crossing_energy_minus_enthalpy_j_kg);
        let expected_energy =
            expected_mass * (h + 0.5 * velocity.into_iter().map(|x| x * x).sum::<f64>());
        if normal_velocity == 0.0 {
            assert_eq!(actual_energy, 0.0);
        } else {
            relative(actual_energy, expected_energy, 3e-13);
        }
    }
}

#[test]
fn rotation_covariance_preserves_the_incident_beam() {
    let rotate = |v: [f64; 3]| [0.6 * v[0] - 0.8 * v[1], 0.8 * v[0] + 0.6 * v[1], v[2]];
    let velocity = [3.0, -2.0, 4.0];
    let original = incident_maxwellian(1.0, 2.0, 3.0, velocity, NORMAL).unwrap();
    let rotated = incident_maxwellian(1.0, 2.0, 3.0, rotate(velocity), rotate(NORMAL)).unwrap();
    relative(
        rotated.incoming_mass_flux_kg_m2_s,
        original.incoming_mass_flux_kg_m2_s,
        3e-14,
    );
    relative(
        rotated.incoming_mass_flux_change_from_zero_drift_kg_m2_s,
        original.incoming_mass_flux_change_from_zero_drift_kg_m2_s,
        3e-14,
    );
    let original = original.moments.unwrap();
    let rotated = rotated.moments.unwrap();
    for (actual, expected) in rotated
        .mean_relative_velocity_m_s
        .into_iter()
        .zip(rotate(original.mean_relative_velocity_m_s))
    {
        relative(actual, expected, 3e-14);
    }
    relative(
        rotated.crossing_energy_minus_enthalpy_j_kg,
        original.crossing_energy_minus_enthalpy_j_kg,
        3e-14,
    );
}

#[test]
fn tiny_normal_drift_survives_equal_rounded_gross_fluxes() {
    for (rho, r, temperature, speed) in [
        (1e200, 1.0, 1.0, 1e-200),
        // The dimensionless speed itself underflows, while the physical
        // first-order mass increment remains representable.
        (1e150, 1e300, 1.0, 1e-200),
    ] {
        let baseline = incident_maxwellian(rho, r, temperature, [0.0; 3], NORMAL).unwrap();
        for sign in [-1.0, 1.0] {
            let beam =
                incident_maxwellian(rho, r, temperature, [sign * speed, 0.0, 0.0], NORMAL).unwrap();
            assert_eq!(
                beam.incoming_mass_flux_kg_m2_s,
                baseline.incoming_mass_flux_kg_m2_s
            );
            // dJ/dv_n at v_n=0 is exactly -rho/2. Higher terms here are
            // hundreds of orders below the precision of this independent oracle.
            relative(
                beam.incoming_mass_flux_change_from_zero_drift_kg_m2_s,
                -0.5 * sign * (rho * speed),
                4e-15,
            );
        }
    }
}

#[test]
fn scaled_far_tail_retains_representable_flux_and_conditional_moments() {
    let z: f64 = 40.0;
    let rho = 1e200;
    let beam = incident_maxwellian(rho, 1.0, 1.0, [z, 0.0, 0.0], NORMAL).unwrap();
    // Independent scaled integral with t=z*y. Keeping exp(-800) split
    // avoids erasing this representable physical flux before density enters.
    let count = 32_768;
    let step = 60.0 / count as f64;
    let mut integrals = [0.0; 3];
    for i in 0..=count {
        let t = i as f64 * step;
        let weight = if i == 0 || i == count {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        let mut value = weight * t * (-t - 0.5 * (t / z) * (t / z)).exp();
        for integral in &mut integrals {
            *integral += value;
            value *= t;
        }
    }
    let integrals = integrals.map(|x| x * step / 3.0);
    let expected_mass = ((rho * (-400.0_f64).exp()) * (-400.0_f64).exp()) * integrals[0]
        / (z * z * (2.0 * core::f64::consts::PI).sqrt());
    relative(beam.incoming_mass_flux_kg_m2_s, expected_mass, 3e-11);
    let moments = beam.moments.unwrap();
    relative(
        moments.mean_relative_velocity_m_s[0],
        -integrals[1] / (z * integrals[0]),
        3e-11,
    );
    relative(
        moments.crossing_energy_minus_enthalpy_j_kg,
        0.5 * integrals[2] / (z * z * integrals[0]) - 1.5,
        3e-13,
    );
    relative(
        beam.incoming_mass_flux_change_from_zero_drift_kg_m2_s,
        -rho / (2.0 * core::f64::consts::PI).sqrt(),
        3e-15,
    );
}

#[test]
fn cold_stream_and_finite_extreme_scales_do_not_require_squared_thermal_speed() {
    for (r, temperature) in [(1e-200, 1e-200), (f64::from_bits(1), f64::from_bits(1))] {
        let beam = incident_maxwellian(2.0, r, temperature, [-3.0, 2.0, 0.0], NORMAL).unwrap();
        assert_eq!(beam.incoming_mass_flux_kg_m2_s, 6.0);
        assert_eq!(
            beam.moments.unwrap().mean_relative_velocity_m_s,
            [-3.0, 2.0, 0.0]
        );
        assert_eq!(
            beam.moments.unwrap().crossing_energy_minus_enthalpy_j_kg,
            6.5
        );
    }
    let beam = incident_maxwellian(1e-300, 1e300, 1.0, [0.0; 3], NORMAL).unwrap();
    relative(
        beam.incoming_mass_flux_kg_m2_s,
        1e-150 / (2.0 * core::f64::consts::PI).sqrt(),
        4e-15,
    );
    relative(
        beam.moments.unwrap().crossing_energy_minus_enthalpy_j_kg,
        -0.5e300,
        4e-15,
    );
    // rho*sigma alone rounds to zero; the complete incoming product does not.
    let subnormal =
        incident_maxwellian(f64::from_bits(1), 0.25, 1.0, [-1.0, 0.0, 0.0], NORMAL).unwrap();
    assert_eq!(subnormal.incoming_mass_flux_kg_m2_s.to_bits(), 1);
}

#[test]
fn vacuum_stays_empty_but_invalid_or_unrepresentable_active_inputs_refuse() {
    let vacuum = incident_maxwellian(0.0, f64::MAX, f64::MAX, [f64::MAX; 3], NORMAL).unwrap();
    assert_eq!(vacuum.incoming_mass_flux_kg_m2_s, 0.0);
    assert_eq!(
        vacuum.incoming_mass_flux_change_from_zero_drift_kg_m2_s,
        0.0
    );
    assert_eq!(vacuum.moments, None);
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(incident_maxwellian(0.0, bad, 300.0, [0.0; 3], NORMAL).is_err());
        assert!(incident_maxwellian(0.0, 287.0, bad, [0.0; 3], NORMAL).is_err());
    }
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(incident_maxwellian(bad, 287.0, 300.0, [0.0; 3], NORMAL).is_err());
    }
    for normal in [[0.0; 3], [2.0, 0.0, 0.0], [f64::NAN, 0.0, 0.0]] {
        assert!(incident_maxwellian(0.0, 1.0, 1.0, [0.0; 3], normal).is_err());
    }
    assert!(incident_maxwellian(0.0, 1.0, 1.0, [f64::INFINITY, 0.0, 0.0], NORMAL).is_err());
    assert!(incident_maxwellian(f64::MAX, 1e100, 1.0, [0.0; 3], NORMAL).is_err());
    assert!(incident_maxwellian(1.0, 1.0, 1.0, [40.0, 0.0, 0.0], NORMAL).is_err());
    // A nonzero drift increment and a nonzero energy correction each have
    // their own representability requirement; neither may silently disappear.
    assert!(incident_maxwellian(1.0, 1.0, 1.0, [f64::from_bits(1), 0.0, 0.0], NORMAL).is_err());
    assert!(incident_maxwellian(1e200, 1e-200, 1e-200, [0.0; 3], NORMAL).is_err());
}

#[test]
fn three_factor_product_rounds_its_exact_value_once() {
    let quantum = f64::from_bits(1); // 2^-1074
    // Exactly half of 2^-1074: ties to even give 0, so the nonzero product refuses.
    assert!(product(1.0, quantum, 0.5).is_err());
    // 3 * 2^-1076 = 0.75 * 2^-1074 rounds to 2^-1074 and is admitted.
    assert_eq!(product(3.0, quantum, 0.25).unwrap().to_bits(), 1);
    // A normal a*b with a normal result is the earlier (a*b)*c, bit for bit.
    let (a, b, c) = (1.1, 3.3e-5, 0.7);
    assert_eq!(product(a, b, c).unwrap().to_bits(), ((a * b) * c).to_bits());
    // With c = 1 the exact product is a*b, which one IEEE multiply rounds
    // correctly onto the subnormal grid. 2^-1055 is 2^19 quanta; the first two
    // factors put the product exactly half-way, and ties go to the even neighbour.
    let x = f64::from_bits(1 << 19);
    let step = 1.0 / 1_048_576.0; // 2^-20
    for (a, bits) in [
        (1.0 + step, 1 << 19),
        (1.0 + 3.0 * step, (1 << 19) + 2),
        (1.0 + f64::EPSILON, 1 << 19),
    ] {
        let rounded = product(a, x, 1.0).unwrap();
        assert_eq!(rounded.to_bits(), bits);
        assert_eq!(rounded.to_bits(), (a * x).to_bits());
    }
}
