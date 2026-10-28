//! Evaluate small radiative exchanges without subtracting blackbody backgrounds.
//!
//! Write J = B_ref + dJ, G = B_ref * sum(F) + dG and
//! B(T) = B_ref + dB. These are the original gray-diffuse equations in shifted
//! coordinates. An imperfect view-factor row keeps its background leakage;
//! neither geometry nor heat is normalized to manufacture energy closure.

/// B(a) - B(b), retaining relative accuracy as the two temperatures converge.
pub(crate) fn emission_difference(sigma: f64, a: f64, b: f64) -> f64 {
    sigma * (a - b) * (a + b) * a.mul_add(a, b * b)
}

/// Right-hand side for [I - (1-epsilon) F] dJ.
pub(crate) fn shifted_load(
    emissivity: f64,
    emission_delta: f64,
    reference_emission: f64,
    row_sum: f64,
) -> f64 {
    let leakage = reference_emission * (row_sum - 1.0);
    emissivity.mul_add(emission_delta, (1.0 - emissivity) * leakage)
}

/// Net outward heat flux epsilon * (B - G), in W/m².
///
/// Using the emissive balance also avoids subtracting nearly equal J and G on
/// reflective surfaces. The linear-system residual still measures disagreement
/// between this emission law and the independently solved radiosity equation.
pub(crate) fn outward_flux(
    emissivity: f64,
    emission_delta: f64,
    irradiation_delta: f64,
    reference_emission: f64,
    row_sum: f64,
) -> f64 {
    let difference = reference_emission.mul_add(1.0 - row_sum, emission_delta - irradiation_delta);
    emissivity * difference
}

#[cfg(test)]
mod tests {
    use super::*;

    const SIGMA: f64 = 5.670374419e-8;

    #[test]
    fn gray_plate_exchange_survives_near_equilibrium_and_reflective_surfaces() {
        // Independent scalar radiative resistance, with a binomial expansion
        // of (T+dT)^4 - T^4 rather than the production factored identity.
        for (ea, eb) in [(0.8, 0.6), (1e-4, 2e-4), (1e-8, 2e-8)] {
            for delta in [f64::from_bits(300.0_f64.to_bits() + 1) - 300.0, 1e-6, 100.0] {
                let cold = 300.0;
                let hot = cold + delta;
                let dt = hot - cold;
                let exact = SIGMA
                    * dt
                    * (4.0 * cold * cold * cold
                        + dt * (6.0 * cold * cold + dt * (4.0 * cold + dt)))
                    / (1.0 / ea + 1.0 / eb - 1.0);
                for reference_t in [cold, hot] {
                    let reference = SIGMA * reference_t * reference_t * reference_t * reference_t;
                    let emission_hot = emission_difference(SIGMA, hot, reference_t);
                    let emission_cold = emission_difference(SIGMA, cold, reference_t);
                    let load_hot = shifted_load(ea, emission_hot, reference, 1.0);
                    let load_cold = shifted_load(eb, emission_cold, reference, 1.0);
                    // Closed-form solution for the two original radiosity rows.
                    let j_cold = (load_cold + (1.0 - eb) * load_hot) / (ea + eb - ea * eb);
                    let j_hot = load_hot + (1.0 - ea) * j_cold;
                    let hot_flux = outward_flux(ea, emission_hot, j_cold, reference, 1.0);
                    let cold_flux = outward_flux(eb, emission_cold, j_hot, reference, 1.0);
                    assert!(
                        (hot_flux - exact).abs() < 2e-14 * exact,
                        "hot {hot_flux:e}, expected {exact:e}"
                    );
                    assert!(
                        (cold_flux + exact).abs() < 2e-14 * exact,
                        "cold {cold_flux:e}, expected {:e}",
                        -exact
                    );
                }
            }
        }
    }

    #[test]
    fn equal_temperature_has_exactly_zero_load_and_exchange() {
        for temperature in [300.0, 1000.0, 10000.0] {
            let background = SIGMA * temperature * temperature * temperature * temperature;
            assert_eq!(emission_difference(SIGMA, temperature, temperature), 0.0);
            for emissivity in [1.0, 0.6, 1e-8] {
                assert_eq!(shifted_load(emissivity, 0.0, background, 1.0), 0.0);
                assert_eq!(outward_flux(emissivity, 0.0, 0.0, background, 1.0), 0.0);
            }
        }
    }

    #[test]
    fn imperfect_view_factor_rows_retain_physical_equation_defects() {
        let reference: f64 = 1000.0;
        let emissivity: f64 = 0.7;
        for row_sum in [0.999, 1.001] {
            let expected_load =
                emissivity * 10.0 + (1.0 - emissivity) * reference * (row_sum - 1.0);
            assert!(
                (shifted_load(emissivity, 10.0, reference, row_sum) - expected_load).abs() < 1e-14
            );
            let original_flux = emissivity * ((reference + 10.0) - (reference * row_sum + 8.0));
            let flux = outward_flux(emissivity, 10.0, 8.0, reference, row_sum);
            assert!((flux - original_flux).abs() < 1e-12);
            assert_ne!(outward_flux(emissivity, 0.0, 0.0, reference, row_sum), 0.0);
        }
    }
}
