//! Gas-related losses: cylinder motion and free-field observer absorption.
//!
//! `H(ω) = exp(−α(ω) r)` with ISO 9613-1 when the `(T, p, humidity)`
//! window admits it. Humidity is an explicit `[0, 1]` argument.
//! ISO already includes a classical term — this path does **not**
//! add Stokes–Kirchhoff on top. Outside the ISO meteorological
//! band the Stokes–Kirchhoff `α ~ ω²` law is the fallback.

use fs_fft::{C64 as FftC64, Fft};
use fs_material::gas::GasState;
use fs_material::iso9613::iso9613_absorption;
use fs_math::det;

/// Mechanical resistance per unit length [N s/m²] of a slender circular
/// cylinder oscillating transversely in stationary gas.
///
/// Uses `R = 2 pi mu (1 + r sqrt(2 rho omega / mu))`, the small-amplitude
/// boundary-layer approximation in Desvages (2018), section 3.2.2, equations
/// 3.8 and 3.12: <https://hdl.handle.net/1842/31273>.
/// Here `mu` is dynamic viscosity [Pa s], not kinematic viscosity [m²/s].
/// For linear mass density `m_l`, modal damping is `zeta = R/(2 m_l omega)`.
///
/// This is an estimate for a long cylinder, including the model's constant
/// low-frequency continuation. It does not solve end effects, fluid added mass,
/// turbulence, confinement, rarefied gas, or large-amplitude motion. Applying it
/// to a nonlinear oscillator retains a linear, modal loss approximation.
///
/// # Errors
/// Refuses nonpositive/nonfinite radius [m], frequency [rad/s], gas viscosity
/// or density, and unrepresentable derived resistance.
pub fn oscillating_cylinder_air_resistance_per_length(
    radius_m: f64,
    omega_rad_s: f64,
    gas: &GasState,
) -> Result<f64, crate::acoustic_realize::AcousticRealizeError> {
    use crate::acoustic_realize::AcousticRealizeError;
    if [radius_m, omega_rad_s, gas.dynamic_viscosity, gas.density]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "cylinder air resistance needs positive finite radius, angular frequency, viscosity and gas density",
        });
    }
    let resistance = core::f64::consts::TAU
        * gas.dynamic_viscosity
        * (1.0 + radius_m * det::sqrt(2.0 * gas.density * omega_rad_s / gas.dynamic_viscosity));
    if !resistance.is_finite() || resistance <= 0.0 {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "cylinder air resistance is unrepresentable",
        });
    }
    Ok(resistance)
}

/// Apply atmospheric absorption to a finished pressure history.
///
/// The transfer is real and even; the block IFFT is the linear
/// zero-phase map of that transfer on the DFT grid. `range_m ≤ 0`
/// is a no-op.
pub fn absorb_pressure_history(
    pressure_pa: &mut [f64],
    dt: f64,
    range_m: f64,
    gas: &GasState,
    relative_humidity: f64,
) {
    if pressure_pa.len() < 4 || !(dt > 0.0 && range_m > 0.0) {
        return;
    }
    let n = pressure_pa.len();
    // The transform covers the complete history. Capping below `n` would
    // truncate the input and make the writeback index beyond the FFT buffer.
    let n_fft = n.next_power_of_two().max(8);
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
        let alpha = iso9613_absorption(gas, relative_humidity, omega)
            .unwrap_or_else(|_| gas.stokes_kirchhoff_absorption(omega));
        let atten = det::exp(-alpha * range_m);
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
            .map(|i| (core::f64::consts::TAU * 2_000.0 * f64::from(i) * dt).sin())
            .collect();
        let mut far = near.clone();
        absorb_pressure_history(&mut near, dt, 1.0, &gas, 0.50);
        absorb_pressure_history(&mut far, dt, 200.0, &gas, 0.50);
        let en: f64 = near.iter().map(|x| x * x).sum();
        let ef: f64 = far.iter().map(|x| x * x).sum();
        assert!(
            ef < en * 0.98,
            "200 m must absorb more 2 kHz energy than 1 m"
        );
    }

    #[test]
    fn g0_history_beyond_legacy_fft_cap_is_transformed_in_full() {
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 288.15, 101_325.0).expect("air");
        let mut pressure = vec![0.0; 8_193];
        pressure[8_192] = 1.0;

        absorb_pressure_history(&mut pressure, 1.0 / 48_000.0, 200.0, &gas, 0.50);

        assert_eq!(pressure.len(), 8_193);
        assert!(pressure.iter().all(|sample| sample.is_finite()));
        assert!(pressure[8_192].abs() < 1.0);
    }
}
