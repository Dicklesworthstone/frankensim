//! Causal uniform-piston radiation in an infinite rigid baffle.
//!
//! `p = rho*c/A * (Q - K*Q)`, where the regular impulse kernel is
//! `K(t) = 4/(pi*T) sqrt(1-(t/T)^2)`, `0 <= t <= T = 2a/c`.
//! This is the inverse transform of the shared Bessel/Struve piston law,
//! not a resistance frozen at one resonance. See Lindemann, JASA 55,
//! 708–717 (1974), https://doi.org/10.1121/1.1914587.
//!
//! Positive linear-interpolation weights give the discrete load
//! `p_n = sum_j g_j (Q_n - Q_{n-j})`, `g_j >= 0`. Each delay branch has
//! storage `dt*g_j/2 * sum_{k=0}^{j-1} Q_{n-k}^2` and dissipation
//! `dt*g_j/2 * (Q_n-Q_{n-j})^2`. Thus the realized load is passive on
//! the actual clock, including its memory, without fitting a pole model.

use std::f64::consts::PI;

#[derive(Clone)]
pub(super) struct BaffledRadiation {
    pub(super) dt: f64,
    pub(super) area: f64,
    pub(super) density: f64,
    gain: f64,
    delayed_gain: Vec<f64>,
    storage_gain: Vec<f64>,
    history: Vec<f64>,
    next: usize,
    pub(super) stored_j: f64,
    pub(super) radiated_j: f64,
}

impl BaffledRadiation {
    pub(super) fn new(area: f64, density: f64, sound_speed: f64, dt: f64) -> Result<Self, String> {
        if ![area, density, sound_speed, dt]
            .iter()
            .all(|v| v.is_finite() && *v > 0.0)
        {
            return Err(
                "baffled vent radiation needs finite positive area, gas and time step".into(),
            );
        }
        let transit = 2.0 * (area / PI).sqrt() / sound_speed;
        let impedance = density * sound_speed / area;
        if !transit.is_finite() || transit <= 0.0 || !impedance.is_finite() || impedance <= 0.0 {
            return Err("baffled vent radiation scale is not representable".into());
        }
        let weights = fs_math::special::semicircle_hat_weights(dt / transit, 4096)?;
        // Express the direct term as the sum of delayed gains. This preserves
        // the kernel's exact DC cancellation instead of subtracting 1-w_0.
        let delayed_gain: Vec<f64> = weights[1..].iter().map(|w| impedance * w).collect();
        let mut storage_gain = vec![0.0; delayed_gain.len()];
        let mut tail = 0.0;
        for j in (0..delayed_gain.len()).rev() {
            tail += delayed_gain[j];
            storage_gain[j] = 0.5 * dt * tail;
        }
        if !tail.is_finite() || tail <= 0.0 || storage_gain.iter().any(|g| !g.is_finite()) {
            return Err("baffled vent radiation coefficients are not representable".into());
        }
        Ok(Self {
            dt,
            area,
            density,
            gain: tail,
            history: vec![0.0; delayed_gain.len()],
            delayed_gain,
            storage_gain,
            next: 0,
            stored_j: 0.0,
            radiated_j: 0.0,
        })
    }

    /// `p_n = gain*Q_n - memory`, with Q centered at the same time as p.
    pub(super) fn pressure_coefficients(&self) -> (f64, f64) {
        let n = self.history.len();
        let memory = self
            .delayed_gain
            .iter()
            .enumerate()
            .map(|(j, g)| g * self.history[(self.next + n - 1 - j) % n])
            .sum();
        (self.gain, memory)
    }

    pub(super) fn advance(&mut self, mean_flow: f64) {
        let n = self.history.len();
        let mut stored = self.storage_gain[0] * mean_flow * mean_flow;
        let mut radiated = 0.0;
        for j in 0..n {
            let past = self.history[(self.next + n - 1 - j) % n];
            radiated += 0.5 * self.dt * self.delayed_gain[j] * (mean_flow - past).powi(2);
            if j + 1 < n {
                stored += self.storage_gain[j + 1] * past * past;
            }
        }
        self.stored_j = stored;
        self.radiated_j += radiated;
        self.history[self.next] = mean_flow;
        self.next = (self.next + 1) % n;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn impulse_kernel_matches_bessel_struve_impedance_and_refines() {
        // G1: independent frequency-domain special functions, including both
        // resistance and added mass; no fitted resonance or calibration.
        let (radius, density, c) = (0.025, 1.225, 343.0);
        let area = PI * radius * radius;
        let transit = 2.0 * radius / c;
        let impedance = density * c / area;
        let mut previous = None;
        for cells in [16, 32, 64] {
            let dt = transit / cells as f64;
            let weights =
                fs_math::special::semicircle_hat_weights(1.0 / cells as f64, 4096).unwrap();
            assert!(weights.iter().all(|w| *w >= 0.0));
            assert!((weights.iter().sum::<f64>() - 1.0).abs() < 2e-13);
            let first: f64 = weights
                .iter()
                .enumerate()
                .map(|(j, w)| j as f64 * dt * w)
                .sum();
            assert!((first / transit - 4.0 / (3.0 * PI)).abs() < 2e-13);
            let load = BaffledRadiation::new(area, density, c, dt).unwrap();
            for ka in [0.02, 0.2, 1.5, 3.0, 6.0] {
                let theta = ka * c / radius * dt;
                let (mut re, mut im) = (0.0, 0.0);
                for (j, gain) in load.delayed_gain.iter().enumerate() {
                    let phase = theta * (j + 1) as f64;
                    re += gain / impedance * (1.0 - phase.cos());
                    im += gain / impedance * phase.sin();
                }
                let expected_re = fs_math::special::one_minus_two_j1_over_x(2.0 * ka);
                let expected_im = fs_math::special::struve_h1(2.0 * ka) / ka;
                let error = (re - expected_re).hypot(im - expected_im);
                // A unit-amplitude exponential's linear interpolation error
                // is bounded by (omega*dt)^2/8; the positive kernel integrates to one.
                assert!(
                    error <= theta * theta / 8.0 + 5e-12,
                    "cells={cells} ka={ka}: impedance error {error}"
                );
                if ka == 1.5 {
                    if let Some(old) = previous {
                        let ratio = old / error;
                        assert!(
                            (3.8..4.2).contains(&ratio),
                            "clock refinement ratio {ratio}"
                        );
                    }
                    previous = Some(error);
                }
            }
        }
        assert!(previous.unwrap() < 3e-4);
        assert!(fs_math::special::semicircle_hat_weights(1e-9, 4096).is_err());
    }

    #[test]
    fn delay_memory_and_emission_close_arbitrary_pressure_flow_work() {
        // G0: prescribed flow with a reversal, a DC interval, and a drained tail.
        let mut load =
            BaffledRadiation::new(PI * 0.025_f64.powi(2), 1.225, 343.0, 1.0 / 96000.0).unwrap();
        let mut work = 0.0;
        let mut scale: f64 = 0.0;
        for i in 0..600 {
            let flow = match i {
                0..=199 => 0.001 * (i as f64 * 0.17).sin(),
                200..=399 => -0.0003,
                _ => 0.0,
            };
            let (gain, memory) = load.pressure_coefficients();
            let pressure = gain * flow - memory;
            let increment = load.dt * pressure * flow;
            work += increment;
            scale += increment.abs();
            load.advance(flow);
            assert!(load.stored_j >= 0.0 && load.radiated_j >= 0.0);
            assert!((load.stored_j + load.radiated_j - work).abs() <= 2e-13 * scale.max(1e-30));
            if i == 399 {
                assert!(pressure.abs() < 1e-12, "DC radiation must vanish");
            }
        }
        assert_eq!(load.stored_j, 0.0);
        assert!(load.radiated_j > 0.0);
    }
}
