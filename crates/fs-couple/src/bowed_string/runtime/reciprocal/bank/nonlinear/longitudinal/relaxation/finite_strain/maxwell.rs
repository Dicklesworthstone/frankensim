//! Source Maxwell strain integrated by the shared exact material-ramp owner.
//! The geometry adapter supplies a work-conjugate midpoint-coordinate pullback.
#[cfg(test)]
use fs_nlmodal::longitudinal::relaxation::maxwell_step as step;

#[cfg(test)]
mod tests {
    use super::step;

    #[test]
    fn g0_step_closes_storage_work_and_positive_heat() {
        for rate in [1e-8, 0.001, 0.2, 1.0, 2.0, 20.0, 1e6] {
            for sign in [-1.0, 1.0] {
                let (a, b, z, c) = (sign * 0.001, sign * 0.003, sign * 0.0007, 2500.0);
                let (next, stress, heat) = step(a, b, z, rate, c).unwrap();
                let initial = 0.5 * c * (a - z).powi(2);
                let final_energy = 0.5 * c * (b - next).powi(2);
                assert!(heat >= 0.0);
                assert!(next >= a.min(b).min(z) && next <= a.max(b).max(z));
                assert!((final_energy - initial + heat - stress * (b - a)).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn g1_fixed_strain_relaxation_has_the_original_time_constant() {
        let (epsilon, z0, c, total, tau) = (0.002, 0.0003, 1700.0, 0.002, 0.001);
        for steps in [1, 4, 16, 64] {
            let mut z = z0;
            let rate = total / (steps as f64 * tau);
            let mut heat = 0.0;
            for _ in 0..steps {
                let next = step(epsilon, epsilon, z, rate, c).unwrap();
                assert!(next.0 >= z && next.0 <= epsilon);
                z = next.0;
                heat += next.2;
            }
            let expected = epsilon + (z0 - epsilon) * (-total / tau).exp();
            assert!((z - expected).abs() < 8e-18);
            assert!(
                (0.5 * c * (epsilon - z).powi(2) + heat - 0.5 * c * (epsilon - z0).powi(2)).abs()
                    < 1e-14
            );
        }
    }

    #[test]
    fn g1_endpoint_strain_tangent_is_the_exact_ramp_modulus() {
        let (a, b, z, rate, c) = (0.001, 0.003, 0.0007, 0.25_f64, 2500.0);
        let h = 1e-7;
        let plus = step(a, b + h, z, rate, c).unwrap();
        let minus = step(a, b - h, z, rate, c).unwrap();
        let relaxation = (1.0 - (-rate).exp()) / rate;
        let gain = (1.0 - relaxation) / rate;
        assert!(((plus.1 - minus.1) / (2.0 * h) - c * gain).abs() < 1e-8);
        assert!(((plus.0 - minus.0) / (2.0 * h) - (1.0 - relaxation)).abs() < 1e-10);
    }

    #[test]
    fn g0_relaxed_state_and_elastic_stiff_limits_are_physical() {
        assert_eq!(step(0.002, 0.002, 0.002, 0.2, 1500.0), Some((0.002, 0.0, 0.0)));
        let elastic = step(0.001, 0.003, 0.0007, 0.0, 2500.0).unwrap();
        assert_eq!(elastic.0, 0.0007);
        assert!((elastic.1 - 2500.0 * (0.002 - 0.0007)).abs() < 1e-14);
        assert_eq!(elastic.2, 0.0);
        let relaxed = step(0.001, 0.003, 0.0007, f64::INFINITY, 2500.0).unwrap();
        assert_eq!(relaxed.0, 0.003);
        assert_eq!(relaxed.1, 0.0);
        assert!((relaxed.2 - 0.5 * 2500.0 * (0.001_f64 - 0.0007).powi(2)).abs() < 1e-18);
        for rate in [-1.0, f64::NEG_INFINITY, f64::NAN] {
            assert!(step(0.0, 0.0, 0.0, rate, 1.0).is_none());
        }
        for stiffness in [0.0, -1.0, f64::INFINITY, f64::NAN] {
            assert!(step(0.0, 0.0, 0.0, 0.2, stiffness).is_none());
        }
        assert!(step(f64::NAN, 0.0, 0.0, 0.2, 1.0).is_none());
    }
}
