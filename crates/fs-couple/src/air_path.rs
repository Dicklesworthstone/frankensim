//! Classical viscothermal absorption on a compact observer path.
//!
//! `H(ω) = exp(−α(ω) r)` with Stokes–Kirchhoff `α ~ ω²` from
//! [`GasState`]. A 10 m path and a 0.1 m path are the same law.
//! This is not ISO 9613 humidity relaxation.

use fs_fft::{C64 as FftC64, Fft};
use fs_material::gas::GasState;
use fs_math::det;

/// Apply Stokes–Kirchhoff absorption to a finished pressure history.
///
/// The transfer is real and even; the block IFFT is the linear
/// zero-phase map of that transfer on the DFT grid. `range_m ≤ 0`
/// is a no-op.
pub fn absorb_pressure_history(pressure_pa: &mut [f64], dt: f64, range_m: f64, gas: &GasState) {
    if pressure_pa.len() < 4 || !(dt > 0.0 && range_m > 0.0) {
        return;
    }
    let n = pressure_pa.len();
    let n_fft = n.next_power_of_two().clamp(8, 8192);
    let fft = Fft::new(n_fft);
    let mut buf: Vec<FftC64> = (0..n_fft)
        .map(|i| {
            if i < n {
                FftC64::new(pressure_pa[i], 0.0)
            } else {
                FftC64::new(0.0, 0.0)
            }
        })
        .collect();
    let mut scratch = vec![FftC64::new(0.0, 0.0); n_fft];
    fft.forward(&mut buf, &mut scratch);
    for (k, bin) in buf.iter_mut().enumerate().take(n_fft / 2 + 1) {
        if k == 0 {
            continue;
        }
        let omega = core::f64::consts::TAU * k as f64 / (n_fft as f64 * dt);
        let atten = det::exp(-gas.stokes_kirchhoff_absorption(omega) * range_m);
        *bin = FftC64::new(bin.re * atten, bin.im * atten);
    }
    for k in 1..n_fft / 2 {
        let c = buf[k];
        buf[n_fft - k] = FftC64::new(c.re, -c.im);
    }
    fft.inverse(&mut buf, &mut scratch);
    for (i, sample) in pressure_pa.iter_mut().enumerate() {
        *sample = buf[i].re;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_material::gas::{GasSpec, GasState};

    #[test]
    fn farther_path_kills_high_frequency_more() {
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 288.15, 101_325.0).expect("air");
        let n = 256;
        let dt = 1.0 / 8_000.0;
        let mut near: Vec<f64> = (0..n)
            .map(|i| (core::f64::consts::TAU * 2_000.0 * i as f64 * dt).sin())
            .collect();
        let mut far = near.clone();
        absorb_pressure_history(&mut near, dt, 1.0, &gas);
        absorb_pressure_history(&mut far, dt, 200.0, &gas);
        let en: f64 = near.iter().map(|x| x * x).sum();
        let ef: f64 = far.iter().map(|x| x * x).sum();
        assert!(
            ef < en * 0.98,
            "200 m must absorb more 2 kHz energy than 1 m"
        );
    }
}
