//! Explicit Hoffman–Lauritzen temperature dependence for a Nakamura rate.
//!
//! Within the admitted crystallization window,
//! `K(T) = K0 exp[-U*/(R(T-Tinf)) - Kg/(T(Tm-T)f)]`,
//! where `f = 2T/(Tm+T)`. The first barrier arrests chain transport near Tinf;
//! the second suppresses growth near equilibrium melting. This opt-in closure
//! extends continuously by zero at/below Tinf and at/above Tm, but only inside
//! the caller's declared source domain. That extension is part of this model,
//! not permission to extrapolate a measured curve or a complete glass model.
//!
//! All coefficients require material/process admission. In particular, K0 is
//! the prefactor of macroscopic Nakamura exposure in inverse seconds. A linear
//! crystal-growth prefactor in metres/second, an Avrami coefficient in s^-n or
//! inverse half-time cannot be substituted without an explicit conversion.
//! No conventional U* or Tg-minus-offset value is silently supplied. Separate
//! kinetic regimes, secondary crystallization, pressure and flow dependence
//! need their own admitted laws; this closure does not infer them.
//!
//! Equation: <https://doi.org/10.3390/polym15071603>, equation 17;
//! also <https://doi.org/10.1021/acsomega.5c00067>, equation 2.

use super::{CrystallizationError, CrystallizationRate, nonnegative, positive};

/// Authored parameters of one continuous Hoffman–Lauritzen kinetic regime.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HoffmanLauritzenConfig {
    /// Nakamura exposure-rate prefactor K0 [s^-1], positive and finite.
    pub prefactor_s_inv: f64,
    /// Segmental transport activation energy U* [J/mol], strictly positive.
    pub transport_activation_energy_j_mol: f64,
    /// Hypothetical mobility-arrest temperature Tinf [K], positive and below Tm.
    /// This is supplied directly, not inferred from a glass-transition point.
    pub mobility_limit_temperature_k: f64,
    /// Secondary surface-nucleation barrier Kg [K²], strictly positive.
    pub nucleation_barrier_k2: f64,
    /// Equilibrium melting limit Tm [K]. The rate vanishes at/above this limit.
    pub melting_temperature_k: f64,
    /// Complete supported absolute-temperature domain, including any explicitly
    /// admitted arrested or superheated portion. No query outside it is allowed.
    pub temperature_range_k: [f64; 2],
}

/// Bounded scalar rate evaluation, with no allocations or iterative work.
///
/// Preparation computes invariant logarithms and U*/R once. Queries use the
/// deterministic elementary-function owner and one exponential. Strong
/// attenuation is combined with log(K0) before exponentiation so a large
/// prefactor cannot lose a representable rate to premature underflow.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HoffmanLauritzenRate {
    config: HoffmanLauritzenConfig,
    log_prefactor: f64,
    transport_activation_temperature_k: f64,
}

impl HoffmanLauritzenRate {
    pub fn new(config: HoffmanLauritzenConfig) -> Result<Self, CrystallizationError> {
        let [minimum_k, maximum_k] = config.temperature_range_k;
        // The molar gas constant is N_A*k_B in SI, rather than a fitted kinetic
        // parameter. Keeping the exact SI defining decimals avoids importing a
        // material-dependent or inconsistent energy-unit convention.
        const MOLAR_GAS_CONSTANT_J_MOL_K: f64 = 8.314_462_618_153_24;
        let activation_k = config.transport_activation_energy_j_mol / MOLAR_GAS_CONSTANT_J_MOL_K;
        if !positive(config.prefactor_s_inv)
            || !positive(config.transport_activation_energy_j_mol)
            || !positive(activation_k)
            || !positive(config.nucleation_barrier_k2)
            || !positive(config.mobility_limit_temperature_k)
            || !positive(config.melting_temperature_k)
            || config.mobility_limit_temperature_k >= config.melting_temperature_k
            || !positive(minimum_k)
            || !positive(maximum_k)
            || minimum_k >= maximum_k
        {
            return Err(CrystallizationError::Invalid(
                "invalid Hoffman-Lauritzen parameters or temperature domain",
            ));
        }
        Ok(Self {
            config,
            log_prefactor: fs_math::det::ln(config.prefactor_s_inv),
            transport_activation_temperature_k: activation_k,
        })
    }

    /// Original authored parameters; the combined kinetic-law identity must
    /// bind these and their source applicability, as for every authored rate.
    #[must_use]
    pub const fn config(&self) -> HoffmanLauritzenConfig {
        self.config
    }
}

impl CrystallizationRate for HoffmanLauritzenRate {
    fn rate_s_inv(&self, temperature_k: f64) -> Result<f64, CrystallizationError> {
        let [minimum_k, maximum_k] = self.config.temperature_range_k;
        if !temperature_k.is_finite() || temperature_k < minimum_k || temperature_k > maximum_k {
            return Err(CrystallizationError::TemperatureOutOfRange {
                temperature_k,
                minimum_k,
                maximum_k,
            });
        }
        if temperature_k <= self.config.mobility_limit_temperature_k
            || temperature_k >= self.config.melting_temperature_k
        {
            return Ok(0.0);
        }

        let mobility = self.transport_activation_temperature_k
            / (temperature_k - self.config.mobility_limit_temperature_k);
        let undercooling = self.config.melting_temperature_k - temperature_k;
        // Kg/(T*deltaT*f) = (Kg/T)/T * (1/2 + T/deltaT).
        // Do not form T², T*deltaT or Tm+T: those can overflow even when
        // the resulting barrier is representable. The divided difference is
        // well-defined arbitrarily close to either admitted kinetic endpoint.
        let nucleation = (self.config.nucleation_barrier_k2 / temperature_k) / temperature_k
            * (0.5 + temperature_k / undercooling);
        let attenuation = mobility + nucleation;
        let rate = if attenuation <= 700.0 {
            // Here the exponential is normal: multiplication preserves small
            // attenuation and avoids log/exp round-trip overflow near MAX.
            self.config.prefactor_s_inv * fs_math::det::exp(-attenuation)
        } else {
            // +infinity is a correctly arrested rate, not nonfinite physics.
            fs_math::det::exp(self.log_prefactor - attenuation)
        };
        if !nonnegative(rate) {
            return Err(CrystallizationError::Invalid(
                "Hoffman-Lauritzen rate is not representable",
            ));
        }
        Ok(rate)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crystallization::{CrystallizationConfig, CrystallizationLaw, CrystallizationSolve};

    fn config() -> HoffmanLauritzenConfig {
        HoffmanLauritzenConfig {
            prefactor_s_inv: 100.0,
            transport_activation_energy_j_mol: 4_157.231_309_076_62,
            mobility_limit_temperature_k: 250.0,
            nucleation_barrier_k2: 20_000.0,
            melting_temperature_k: 500.0,
            temperature_range_k: [200.0, 650.0],
        }
    }

    fn near(actual: f64, expected: f64, relative: f64) {
        assert!(
            (actual - expected).abs() <= relative * expected.abs(),
            "{actual:e} != {expected:e} (relative tolerance {relative:e})",
        );
    }

    #[test]
    fn g1_original_equation_and_competing_barriers() {
        let config = config();
        let rate = HoffmanLauritzenRate::new(config).unwrap();
        assert_eq!(rate.config(), config);
        for temperature in [251.0, 270.0, 300.0, 350.0, 400.0, 450.0, 490.0, 499.0] {
            // Independent direct expression and platform libm oracle.
            let correction = 2.0 * temperature / (config.melting_temperature_k + temperature);
            let expected = config.prefactor_s_inv
                * (-config.transport_activation_energy_j_mol
                    / (8.314_462_618_153_24 * (temperature - config.mobility_limit_temperature_k)))
                    .exp()
                * (-config.nucleation_barrier_k2
                    / (temperature * (config.melting_temperature_k - temperature) * correction))
                    .exp();
            near(rate.rate_s_inv(temperature).unwrap(), expected, 1e-12);
        }
        assert!(rate.rate_s_inv(350.0).unwrap() > rate.rate_s_inv(270.0).unwrap());
        assert!(rate.rate_s_inv(400.0).unwrap() > rate.rate_s_inv(490.0).unwrap());
    }

    #[test]
    fn g0_arrest_melting_and_exact_source_domain_are_distinct() {
        let rate = HoffmanLauritzenRate::new(config()).unwrap();
        for temperature in [200.0, 249.0, 250.0, 500.0, 510.0, 650.0] {
            assert_eq!(rate.rate_s_inv(temperature).unwrap(), 0.0);
        }
        assert_eq!(rate.rate_s_inv(250.0_f64.next_up()).unwrap(), 0.0);
        assert_eq!(rate.rate_s_inv(500.0_f64.next_down()).unwrap(), 0.0);
        for temperature in [199.0, 651.0, 0.0, f64::INFINITY, f64::NAN] {
            assert!(matches!(
                rate.rate_s_inv(temperature),
                Err(CrystallizationError::TemperatureOutOfRange { .. })
            ));
        }
    }

    #[test]
    fn g0_invalid_parameters_refuse_instead_of_supplying_material_defaults() {
        let original = config();
        for invalid in [
            HoffmanLauritzenConfig {
                prefactor_s_inv: 0.0,
                ..original
            },
            HoffmanLauritzenConfig {
                prefactor_s_inv: f64::INFINITY,
                ..original
            },
            HoffmanLauritzenConfig {
                transport_activation_energy_j_mol: -1.0,
                ..original
            },
            HoffmanLauritzenConfig {
                transport_activation_energy_j_mol: f64::NAN,
                ..original
            },
            HoffmanLauritzenConfig {
                nucleation_barrier_k2: 0.0,
                ..original
            },
            HoffmanLauritzenConfig {
                mobility_limit_temperature_k: 500.0,
                ..original
            },
            HoffmanLauritzenConfig {
                melting_temperature_k: f64::INFINITY,
                ..original
            },
            HoffmanLauritzenConfig {
                temperature_range_k: [200.0, 200.0],
                ..original
            },
            HoffmanLauritzenConfig {
                temperature_range_k: [0.0, 650.0],
                ..original
            },
        ] {
            assert!(HoffmanLauritzenRate::new(invalid).is_err());
        }
    }

    #[test]
    fn g1_nakamura_half_time_latent_heat_and_arrest_use_the_existing_owner() {
        let rate = HoffmanLauritzenRate::new(config()).unwrap();
        let amorphous = |t: f64| Ok(1000.0 * (t - 300.0) + 50_000.0);
        let crystalline = |t: f64| Ok(1000.0 * (t - 300.0));
        let law = CrystallizationLaw::new(
            &amorphous,
            &crystalline,
            &rate,
            CrystallizationConfig {
                identity: [0x72; 32],
                avrami_exponent: 2.5,
                maximum_crystalline_mass_fraction: 0.8,
                temperature_range_k: config().temperature_range_k,
                melting_temperature_k: config().melting_temperature_k,
            },
        )
        .unwrap();
        let half_time_s = 2.0_f64.ln().powf(1.0 / 2.5) / rate.rate_s_inv(350.0).unwrap();
        let initial = law.state_at_temperature(350.0, 0.0).unwrap();
        let step = law.advance_isothermal(initial, half_time_s).unwrap();
        near(step.state.crystalline_mass_fraction(), 0.4, 2e-14);
        near(step.specific_heat_j_kg, -20_000.0, 2e-14);
        law.validate_state(step.state, CrystallizationSolve::default())
            .unwrap();

        // Retained nuclei/exposure do not create a spurious rate below Tinf;
        // superheated amorphous feed also remains admissible without growth.
        let arrested = law.state_at_temperature(240.0, 0.7).unwrap();
        assert_eq!(
            law.advance_isothermal(arrested, 1e6).unwrap().state,
            arrested
        );
        let hot = law.state_at_temperature(600.0, 0.0).unwrap();
        assert_eq!(law.advance_isothermal(hot, 1e6).unwrap().state, hot);
        assert_eq!(law.advance_isothermal(initial, half_time_s).unwrap(), step);
    }

    #[test]
    fn g0_extreme_prefactors_and_barriers_do_not_create_false_growth_or_arrest() {
        let strongly_attenuated = HoffmanLauritzenRate::new(HoffmanLauritzenConfig {
            prefactor_s_inv: 1e300,
            transport_activation_energy_j_mol: 8.314_462_618_153_24 * 150.0 * 800.0,
            nucleation_barrier_k2: 1.0,
            ..config()
        })
        .unwrap();
        let expected = (1e300_f64.ln() - 800.0 - 1.0 / (400.0 * 100.0 * (800.0 / 900.0))).exp();
        let actual = strongly_attenuated.rate_s_inv(400.0).unwrap();
        assert!(
            actual > 0.0,
            "representable rate must survive exp(-800) underflow"
        );
        near(actual, expected, 3e-13);

        let tiny_barriers = HoffmanLauritzenRate::new(HoffmanLauritzenConfig {
            prefactor_s_inv: f64::MAX,
            transport_activation_energy_j_mol: 1e-200,
            nucleation_barrier_k2: 1e-200,
            ..config()
        })
        .unwrap();
        assert_eq!(tiny_barriers.rate_s_inv(400.0).unwrap(), f64::MAX);

        // A subnormal Kg cannot be halved to zero before the temperature
        // factors amplify it. No practical calibration is claimed here.
        let scaled = HoffmanLauritzenRate::new(HoffmanLauritzenConfig {
            prefactor_s_inv: 1.0,
            transport_activation_energy_j_mol: 1e-320,
            mobility_limit_temperature_k: 1e-171,
            nucleation_barrier_k2: f64::from_bits(1),
            melting_temperature_k: 1e-169,
            temperature_range_k: [1e-171, 2e-169],
        })
        .unwrap();
        assert_eq!(scaled.rate_s_inv(2e-170).unwrap(), 0.0);
    }
}
