//! Gas-related losses: cylinder motion and free-field observer absorption.
//!
//! A generic gas uses its own Stokes-Kirchhoff viscothermal loss. The
//! air-specific ISO 9613 law requires an explicit choice; temperature and
//! pressure alone do not establish that an arbitrary mixture is air.

/// Refinable causal gas inertia and viscous memory for circular transverse motion.
pub mod cylinder;

use fs_fft::{C64 as FftC64, Fft};
use fs_material::gas::GasState;
use fs_material::iso9613::iso9613_absorption;
use fs_math::det;
use crate::acoustic_realize::AcousticRealizeError;

/// Constitutive choice for observer-path attenuation, never inferred from a
/// material name, gamma, molar mass, or an unrelated humidity coordinate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GasAbsorptionModel {
    /// Classical viscothermal absorption from this gas's rho, c, mu, k and cp.
    /// Molecular relaxation, bulk viscosity and chemical relaxation are not
    /// supplied by this model. No air-specific coefficients are substituted.
    StokesKirchhoff,
    /// Explicit caller assertion of terrestrial air and ISO applicability.
    /// The ISO term already includes classical loss; it must not be added to
    /// Stokes-Kirchhoff. Out-of-window requests refuse rather than switch laws.
    Iso9613Air {
        /// Relative humidity as a fraction in [0, 1], not water mole fraction.
        relative_humidity: f64,
    },
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

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
/// Refuses nonpositive/nonfinite radius [m], frequency [rad/s], incoherent gas
/// states, and unrepresentable derived resistance.
pub fn oscillating_cylinder_air_resistance_per_length(
    radius_m: f64,
    omega_rad_s: f64,
    gas: &GasState,
) -> Result<f64, AcousticRealizeError> {
    gas.validate_coherence().map_err(|_| invalid("cylinder resistance needs a coherent gas state"))?;
    if [radius_m, omega_rad_s, gas.dynamic_viscosity, gas.density]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err(invalid("cylinder air resistance needs positive finite radius, angular frequency, viscosity and gas density"));
    }
    let resistance = core::f64::consts::TAU
        * gas.dynamic_viscosity
        * (1.0 + radius_m * det::sqrt(2.0 * gas.density * omega_rad_s / gas.dynamic_viscosity));
    if !resistance.is_finite() || resistance <= 0.0 {
        return Err(invalid("cylinder air resistance is unrepresentable"));
    }
    Ok(resistance)
}

/// Apply an explicitly chosen absorption law to the entire pressure history.
///
/// The transfer is real and even: this is the existing zero-phase map on the
/// padded DFT grid, not a causal propagation or dispersion solver. A zero range
/// and histories shorter than four samples are no-ops after input admission.
/// Source/material receipts remain with the caller's evaluated gas model.
///
/// # Errors
/// Invalid samples, time/range, incoherent gas, unavailable ISO conditions, or
/// unrepresentable absorption/FFT output refuse WITHOUT changing the history.
pub fn absorb_pressure_history_with_model(
    pressure_pa: &mut [f64],
    dt: f64,
    range_m: f64,
    gas: &GasState,
    model: GasAbsorptionModel,
) -> Result<(), AcousticRealizeError> {
    gas.validate_coherence().map_err(|_| invalid("observer absorption needs a coherent gas state"))?;
    if let GasAbsorptionModel::Iso9613Air { relative_humidity } = model {
        // Admit the declared air conditions even for an empty or zero-range
        // history. A failed air law is not permission to use a different law.
        iso9613_absorption(gas, relative_humidity, 1.0)
            .map_err(|_| invalid("explicit ISO 9613 air conditions are outside the model domain"))?;
    }
    transform_history(pressure_pa, dt, range_m, |omega| match model {
        GasAbsorptionModel::StokesKirchhoff => Ok(gas.stokes_kirchhoff_absorption(omega)),
        GasAbsorptionModel::Iso9613Air { relative_humidity } => {
            iso9613_absorption(gas, relative_humidity, omega)
                .map_err(|_| invalid("explicit ISO 9613 absorption could not be evaluated"))
        }
    })
}

/// Compatibility path for callers that explicitly supply terrestrial-air
/// humidity. Retains the legacy ISO-to-classical fallback, including NaN
/// humidity used by older generic-gas callers. New generic consumers should
/// use [`absorb_pressure_history_with_model`] and [`GasAbsorptionModel`] so
/// physical model selection is not encoded in an invalid numerical value.
///
/// `range_m <= 0` and histories shorter than four samples remain no-ops.
/// Invalid/unrepresentable transformations leave the history unchanged.
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
    let _ = transform_history(pressure_pa, dt, range_m, |omega| {
        Ok(iso9613_absorption(gas, relative_humidity, omega)
            .unwrap_or_else(|_| gas.stokes_kirchhoff_absorption(omega)))
    });
}

// One transform implementation for explicit and legacy model selection.
fn transform_history(
    pressure_pa: &mut [f64],
    dt: f64,
    range_m: f64,
    mut absorption: impl FnMut(f64) -> Result<f64, AcousticRealizeError>,
) -> Result<(), AcousticRealizeError> {
    if !dt.is_finite() || dt <= 0.0 || !range_m.is_finite() || range_m < 0.0
        || pressure_pa.iter().any(|sample| !sample.is_finite())
    {
        return Err(invalid("absorption needs finite pressure, positive dt and nonnegative range"));
    }
    if pressure_pa.len() < 4 || range_m == 0.0 { return Ok(()); }
    let n = pressure_pa.len();
    // Cover the complete history, including samples beyond the old FFT cap.
    let n_fft = n.checked_next_power_of_two()
        .ok_or_else(|| invalid("absorption FFT length is unrepresentable"))?.max(8);
    let fft = Fft::new(n_fft);
    let mut buf: Vec<FftC64> = (0..n_fft)
        .map(|i| FftC64::new(if i < n { pressure_pa[i] } else { 0.0 }, 0.0))
        .collect();
    let mut scratch = vec![FftC64::new(0.0, 0.0); n_fft];
    fft.forward(&mut buf, &mut scratch);
    for (k, bin) in buf.iter_mut().enumerate().take(n_fft / 2 + 1) {
        if k == 0 { continue; }
        let omega = core::f64::consts::TAU * k as f64 / (n_fft as f64 * dt);
        if !omega.is_finite() || omega <= 0.0 {
            return Err(invalid("absorption frequency is unrepresentable"));
        }
        let alpha = absorption(omega)?;
        if !alpha.is_finite() || alpha < 0.0 {
            return Err(invalid("absorption must be finite and nonnegative"));
        }
        let atten = det::exp(-alpha * range_m);
        *bin = FftC64::new(bin.re * atten, bin.im * atten);
    }
    for k in 1..n_fft / 2 {
        let c = buf[k];
        buf[n_fft - k] = FftC64::new(c.re, -c.im);
    }
    fft.inverse(&mut buf, &mut scratch);
    if buf[..n].iter().any(|value| !value.re.is_finite()) {
        return Err(invalid("absorbed pressure is unrepresentable"));
    }
    for (i, sample) in pressure_pa.iter_mut().enumerate() { *sample = buf[i].re; }
    Ok(())
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
        assert!(ef < en * 0.98, "200 m must absorb more 2 kHz energy than 1 m");
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

    #[test]
    fn generic_transfer_uses_the_actual_material_for_air_and_steam() {
        let dt = 1.0 / 8000.0;
        let omega = core::f64::consts::TAU * 1000.0;
        let input: Vec<_> = (0..256).map(|i| (omega * f64::from(i) * dt).cos()).collect();
        let mut outputs = Vec::new();
        for spec in [GasSpec::dry_air_ussa1976(), GasSpec::water_vapor_nist()] {
            let gas = GasState::try_new(&spec, 450.0, 100000.0).unwrap();
            let expected = det::exp(-gas.stokes_kirchhoff_absorption(omega) * 100.0);
            let mut signal = input.clone();
            absorb_pressure_history_with_model(&mut signal, dt, 100.0, &gas,
                GasAbsorptionModel::StokesKirchhoff).unwrap();
            for (&actual, &before) in signal.iter().zip(&input) {
                assert!((actual - expected * before).abs() < 1e-10);
            }
            outputs.push(signal);
        }
        assert_ne!(outputs[0], outputs[1]);
    }

    #[test]
    fn explicit_air_refuses_instead_of_silently_using_another_law() {
        let gas = GasState::try_new(&GasSpec::water_vapor_nist(), 450.0, 100000.0).unwrap();
        let original = vec![1.0, -1.0, 0.5, -0.5, 1.0, 0.0, -1.0, 0.0];
        let mut signal = original.clone();
        assert!(absorb_pressure_history_with_model(&mut signal, 1.0 / 8000.0, 100.0,
            &gas, GasAbsorptionModel::Iso9613Air { relative_humidity: 0.5 }).is_err());
        assert_eq!(signal, original);
        let mut forged = gas;
        forged.density *= 2.0;
        assert!(absorb_pressure_history_with_model(&mut signal, 1.0 / 8000.0, 100.0,
            &forged, GasAbsorptionModel::StokesKirchhoff).is_err());
        assert_eq!(signal, original);
    }

    #[test]
    fn old_nan_selector_and_explicit_generic_share_the_same_operator() {
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100000.0).unwrap();
        for (humidity, model) in [
            (f64::NAN, GasAbsorptionModel::StokesKirchhoff),
            (0.5, GasAbsorptionModel::Iso9613Air { relative_humidity: 0.5 }),
        ] {
            let mut legacy = vec![1.0, -1.0, 0.5, -0.5, 1.0, 0.0, -1.0, 0.0];
            let mut explicit = legacy.clone();
            absorb_pressure_history(&mut legacy, 1.0 / 8000.0, 100.0, &gas, humidity);
            absorb_pressure_history_with_model(&mut explicit, 1.0 / 8000.0, 100.0, &gas, model).unwrap();
            assert_eq!(legacy, explicit);
        }
    }
}
