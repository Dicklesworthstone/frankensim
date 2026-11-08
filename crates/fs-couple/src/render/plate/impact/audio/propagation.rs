//! Fixed physical travel on the original source clock. The existing Rayleigh
//! path owner supplies the qualified fractional-delay kernel; this leaf only
//! retains and advances its history. No future sample or timing compensation.
use crate::pcm_wav::baffled::path::{HALF, TAPS, taps, window};

pub(super) struct Propagation {
    history: Vec<f64>,
    kernel: [f64; TAPS],
    head: usize,
    whole: usize,
    fraction: f64,
}
impl Propagation {
    /// The history cap includes physical travel and the whole causal FIR.
    pub(super) fn new(delay_samples: f64, maximum_samples: usize) -> Result<Self, &'static str> {
        if !delay_samples.is_finite() || delay_samples < 0.0
            || delay_samples >= usize::MAX as f64
        { return Err("listener propagation exceeds representable delay budget"); }
        let whole = delay_samples.floor() as usize;
        let count = whole.checked_add(TAPS)
            .filter(|&count| count <= maximum_samples)
            .ok_or("listener propagation and reconstruction exceed the delay history budget")?;
        let fraction = delay_samples - whole as f64;
        let kernel = taps(fraction, &window());
        let mut history = Vec::new();
        history.try_reserve_exact(count)
            .map_err(|_| "pressure propagation allocation failed")?;
        history.resize(count, 0.0);
        Ok(Self { history, kernel, head: 0, whole, fraction })
    }
    pub(super) const fn whole_samples(&self) -> usize { self.whole }
    pub(super) const fn fractional_samples(&self) -> f64 { self.fraction }
    pub(super) const fn reconstruction_delay_samples(&self) -> usize { HALF }

    /// Retain one current sample, then evaluate only current/past samples.
    /// Ring traversal uses two slices rather than one modulo per tap.
    pub(super) fn step(&mut self, value: f64) -> f64 {
        self.history[self.head] = value;
        let first = if self.head >= self.whole { self.head - self.whole }
            else { self.head + self.history.len() - self.whole };
        let pressure = if self.fraction == 0.0 {
            let center = if first >= HALF { first - HALF }
                else { first + self.history.len() - HALF };
            self.history[center]
        } else {
            let (before, after) = self.history.split_at(first + 1);
            before.iter().rev().chain(after.iter().rev()).zip(&self.kernel)
                .map(|(sample, weight)| sample * weight).sum()
        };
        self.head += 1;
        if self.head == self.history.len() { self.head = 0; }
        pressure
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::TAU;

    #[test]
    fn fractional_travel_preserves_upper_harmonic_gain_and_phase() {
        for whole in [0, 1, 32, 137] {
            for fraction in [0.0, 0.125, 0.5, 0.875] {
                let physical = whole as f64 + fraction;
                for frequency in [0.05, 0.2, 0.4, 0.45] {
                    let mut path = Propagation::new(physical, 256).unwrap();
                    assert_eq!(path.whole_samples(), whole);
                    assert_eq!(path.fractional_samples(), fraction);
                    assert_eq!(path.reconstruction_delay_samples(), 32);
                    let omega = TAU * frequency;
                    let (mut error, mut signal) = (0.0, 0.0);
                    for n in 0..2048 {
                        let actual = path.step((omega * n as f64).cos());
                        if n >= 256 {
                            let exact = (omega * (n as f64 - physical - 32.0)).cos();
                            error += (actual - exact).powi(2);
                            signal += exact * exact;
                        }
                    }
                    let relative = (error / signal).sqrt();
                    assert!(relative < 4e-4, "delay {physical}, frequency {frequency}: {relative}");
                }
            }
        }
    }

    #[test]
    fn distinct_fractional_paths_keep_coherent_interference() {
        // Rounding both physical paths to 13 samples would erase this entire
        // signed upper-harmonic field, although the true phase differs by 81°.
        let (mut near, mut far) = (Propagation::new(12.125, 128).unwrap(),
            Propagation::new(12.625, 128).unwrap());
        let omega = TAU * 0.45;
        let (mut error, mut signal) = (0.0, 0.0);
        for n in 0..2048 {
            let source = (omega * n as f64).cos();
            let actual = near.step(source) - far.step(source);
            if n >= 256 {
                let exact = (omega * (n as f64 - 44.125)).cos()
                    - (omega * (n as f64 - 44.625)).cos();
                error += (actual - exact).powi(2);
                signal += exact * exact;
            }
        }
        assert!(signal > 1000.0);
        assert!((error / signal).sqrt() < 4e-4);
    }

    #[test]
    fn startup_and_integral_paths_remain_causal_without_compensation() {
        for delay in [0.0, 0.125, 0.875, 1.0, 32.5] {
            let mut path = Propagation::new(delay, 128).unwrap();
            let mut response = [0.0; 128];
            for (n, sample) in response.iter_mut().enumerate() {
                *sample = path.step(if n == 0 { 1.0 } else { 0.0 });
            }
            let first = delay.ceil() as usize;
            assert!(response[..first].iter().all(|&x| x == 0.0));
            assert!((response.iter().sum::<f64>() - 1.0).abs() < 1e-14);
            if delay.fract() == 0.0 {
                assert_eq!(response[delay as usize + 32], 1.0);
                assert_eq!(response.iter().filter(|&&x| x != 0.0).count(), 1);
            }
        }
    }

    #[test]
    fn history_capacity_is_admitted_before_allocating_the_complete_fir() {
        for delay in [f64::NAN, f64::INFINITY, -0.1, usize::MAX as f64] {
            assert!(Propagation::new(delay, 256).is_err());
        }
        assert!(Propagation::new(0.0, TAPS - 1).is_err());
        assert!(Propagation::new(4.75, 4 + TAPS - 1).is_err());
        assert_eq!(Propagation::new(4.75, 4 + TAPS).unwrap().history.len(), 4 + TAPS);
    }
}
