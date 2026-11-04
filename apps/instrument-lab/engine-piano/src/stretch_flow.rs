//! Exact position mean of a held-force axial oscillator.
//!
//! Preparation consumes the existing modal owner's displacement responses to
//! unit held and endpoint-ramp forces over the SAME interval. No oscillator
//! transition is reimplemented here. With impulse response B and held response
//! G, G' = B; the ramp response is (integral G)/h. Therefore
//! mean(q) = q0 + G(h)/h * v0 + ramp(h) * (F - omega^2*q0).
//! The extension force is affine in axial displacement at frozen transverse
//! geometry, so this mean gives its exact integrated reciprocal impulse.

#[derive(Clone, Copy, Debug)]
pub(crate) struct AxialMean {
    velocity_gain_s: f64,
    residual_gain_s2: f64,
}

impl AxialMean {
    pub(crate) fn try_new(held_displacement: f64, ramp_displacement: f64,
        duration_s: f64) -> Option<Self> {
        if !duration_s.is_finite() || duration_s <= 0.0
            || [held_displacement, ramp_displacement].iter().any(|x| !x.is_finite() || *x < 0.0) {
            return None;
        }
        let velocity_gain_s = held_displacement / duration_s;
        velocity_gain_s.is_finite().then_some(Self {
            velocity_gain_s, residual_gain_s2: ramp_displacement,
        })
    }

    /// Original mass-normalized displacement, velocity and held acceleration.
    /// The caller checks representability before publishing a physical state.
    #[inline]
    pub(crate) fn displacement(self, q: f64, v: f64, force: f64, omega2: f64) -> f64 {
        q + self.velocity_gain_s * v + self.residual_gain_s2 * (force - omega2 * q)
    }
}

#[cfg(test)]
mod tests {
    use super::AxialMean;
    use std::f64::consts::TAU;

    // Independent analytic scalar motion and moments. Production preparation
    // obtains both moments from fs-modal's admitted held/ramp transition.
    fn free(omega: f64, zeta: f64, q: f64, v: f64, h: f64) -> [f64; 2] {
        let a = zeta * omega;
        let decay = (-a * h).exp();
        if zeta == 1.0 {
            [decay * (q + h * (v + omega * q)),
                decay * (v - omega * h * (v + omega * q))]
        } else {
            let wd = omega * (1.0 - zeta * zeta).sqrt();
            let (s, c) = (wd * h).sin_cos();
            [decay * (q * c + (v + a * q) * s / wd),
                decay * (v * c - (a * v + omega * omega * q) * s / wd)]
        }
    }

    fn moments(omega: f64, zeta: f64, h: f64) -> AxialMean {
        let p = free(omega, zeta, 1.0, 0.0, h)[0];
        let b = free(omega, zeta, 0.0, 1.0, h)[0];
        let held = (1.0 - p) / (omega * omega);
        let ramp = (h - b - 2.0 * zeta * omega * held) / (h * omega * omega);
        AxialMean::try_new(held, ramp, h).unwrap()
    }

    #[test]
    fn axial_mean_matches_independent_damped_trajectory_quadrature() {
        let omega = TAU * 20_000.0;
        let h = 0.5 / 48_000.0;
        let (q, v, equilibrium) = (2e-5, 0.1, -1.2e-5);
        for zeta in [0.0, 1e-10, 0.2, 0.8, 1.0] {
            let mean = moments(omega, zeta, h).displacement(q, v,
                omega * omega * equilibrium, omega * omega);
            let n = 4096;
            let integral = (0..=n).map(|j| {
                let weight = if j == 0 || j == n { 1.0 } else if j % 2 == 0 { 2.0 } else { 4.0 };
                weight * (equilibrium + free(omega, zeta, q - equilibrium, v, h * j as f64 / n as f64)[0])
            }).sum::<f64>() / (3.0 * n as f64);
            assert!((mean - integral).abs() < 2e-12 * (q.abs() + equilibrium.abs()),
                "zeta={zeta}: mean={mean:e}, integral={integral:e}");
        }
        let free_mass = AxialMean::try_new(0.5 * h * h, h * h / 6.0, h).unwrap();
        assert_eq!(free_mass.displacement(0.3, -0.1, 2.0, 0.0),
            0.3 - 0.5 * h * 0.1 + h * h / 3.0);
        assert!(AxialMean::try_new(1.0, 1.0, 0.0).is_none());
        assert!(AxialMean::try_new(f64::NAN, 1.0, h).is_none());
    }

    fn linear_extension(x: &mut [f64; 4], omega: f64, coupling: f64, h: f64, mean: AxialMean) {
        let [q, v, u, w] = *x;
        let force = coupling * q;
        let average = mean.displacement(u, w, force, omega * omega);
        let equilibrium = force / (omega * omega);
        let next = free(omega, 0.0, u - equilibrium, w, h);
        *x = [q, v + h * coupling * average, next[0] + equilibrium, next[1]];
    }

    #[test]
    fn symmetric_reciprocal_subflows_resolve_normal_modes_and_beating_at_second_order() {
        let omega = TAU * 400.0;
        let beta = 0.08_f64;
        let coupling = beta * omega * omega;
        let low = omega * (1.0 - beta).sqrt();
        let high = omega * (1.0 + beta).sqrt();
        let energy0 = 0.5 * omega * omega;
        let (mut errors, mut energy_errors) = (Vec::new(), Vec::new());
        for rate in [32_000, 64_000, 128_000] {
            let h = 1.0 / f64::from(rate);
            let mean = moments(omega, 0.0, 0.5 * h);
            let mut x = [1.0, 0.0, 0.0, 0.0];
            let (mut error, mut energy_error, mut transfer) = (0.0_f64, 0.0_f64, 0.0_f64);
            for j in 1..=rate / 40 {
                linear_extension(&mut x, omega, coupling, 0.5 * h, mean);
                let next = free(omega, 0.0, x[0], x[1], h);
                (x[0], x[1]) = (next[0], next[1]);
                linear_extension(&mut x, omega, coupling, 0.5 * h, mean);
                let t = f64::from(j) * h;
                let (sl, cl) = (low * t).sin_cos();
                let (sh, ch) = (high * t).sin_cos();
                let exact = [0.5 * (cl + ch), -0.5 * (low * sl + high * sh),
                    0.5 * (cl - ch), -0.5 * (low * sl - high * sh)];
                for k in 0..4 {
                    let scale = if k % 2 == 0 { 1.0 } else { omega };
                    error = error.max((x[k] - exact[k]).abs() / scale);
                }
                let [q, v, u, w] = x;
                let energy = 0.5 * (v * v + w * w + omega * omega * (q * q + u * u)) - coupling * q * u;
                energy_error = energy_error.max((energy / energy0 - 1.0).abs());
                transfer = transfer.max(u.abs());
            }
            assert!(transfer > 0.9, "the coupled receiver must acquire the beating energy");
            errors.push(error);
            energy_errors.push(energy_error);
        }
        for values in [&errors, &energy_errors] {
            for pair in values.windows(2) {
                assert!((3.7..4.3).contains(&(pair[0] / pair[1])), "refinement {values:?}");
            }
        }
        assert!(errors[2] < 2e-4, "{errors:?}");
        assert!(energy_errors[2] < 2e-5, "{energy_errors:?}");
    }

    // Positive moderate-slope potential:
    // q^2/2 + 2*u^2 + gamma*u*q^2 + lambda*q^4/4.
    // The axial quadratic can be completed without a negative quartic remainder.
    const GAMMA: f64 = 0.6;
    const LAMBDA: f64 = 1.2;
    fn nonlinear_rhs(x: [f64; 4]) -> [f64; 4] {
        let [q, v, u, w] = x;
        [v, -q - LAMBDA * q.powi(3) - 2.0 * GAMMA * q * u,
            w, -4.0 * u - GAMMA * q * q]
    }
    fn energy(x: [f64; 4]) -> f64 {
        let [q, v, u, w] = x;
        0.5 * (v * v + w * w + q * q) + 2.0 * u * u
            + GAMMA * u * q * q + 0.25 * LAMBDA * q.powi(4)
    }
    fn rk4(x: [f64; 4], h: f64) -> [f64; 4] {
        let a = nonlinear_rhs(x);
        let b = nonlinear_rhs(std::array::from_fn(|i| x[i] + 0.5 * h * a[i]));
        let c = nonlinear_rhs(std::array::from_fn(|i| x[i] + 0.5 * h * b[i]));
        let d = nonlinear_rhs(std::array::from_fn(|i| x[i] + h * c[i]));
        std::array::from_fn(|i| x[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0)
    }
    fn nonlinear_extension(x: &mut [f64; 4], h: f64, mean: AxialMean) {
        let [q, v, u, w] = *x;
        let force = -GAMMA * q * q;
        let average = mean.displacement(u, w, force, 4.0);
        let next = free(2.0, 0.0, u - 0.25 * force, w, h);
        *x = [q, v - h * (LAMBDA * q.powi(3) + 2.0 * GAMMA * q * average),
            next[0] + 0.25 * force, next[1]];
    }

    #[test]
    fn nonlinear_exchange_and_splitting_work_converge_without_added_damping() {
        let initial = [0.6, 0.2, -0.03, 0.1];
        let energy0 = energy(initial);
        let duration = 12.8;
        let mut reference = initial;
        for _ in 0..65_536 { reference = rk4(reference, duration / 65_536.0); }
        assert!((energy(reference) / energy0 - 1.0).abs() < 2e-11);
        let (mut errors, mut energy_errors) = (Vec::new(), Vec::new());
        for steps in [256, 512, 1024] {
            let h = duration / steps as f64;
            let mean = moments(2.0, 0.0, 0.5 * h);
            let mut x = initial;
            let mut work = 0.0_f64;
            for _ in 0..steps {
                nonlinear_extension(&mut x, 0.5 * h, mean);
                let next = free(1.0, 0.0, x[0], x[1], h);
                (x[0], x[1]) = (next[0], next[1]);
                nonlinear_extension(&mut x, 0.5 * h, mean);
                work = work.max((energy(x) / energy0 - 1.0).abs());
            }
            errors.push(x.into_iter().zip(reference).map(|(a, b)| (a - b).powi(2)).sum::<f64>().sqrt());
            energy_errors.push(work);
        }
        for values in [&errors, &energy_errors] {
            for pair in values.windows(2) {
                assert!((3.7..4.3).contains(&(pair[0] / pair[1])), "refinement {values:?}");
            }
        }
        assert!(errors[2] < 5e-5, "{errors:?}");
        assert!(energy_errors[2] < 5e-5, "{energy_errors:?}");
    }
}
