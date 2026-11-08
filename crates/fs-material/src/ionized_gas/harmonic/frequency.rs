//! Analytic frequency derivatives of the original collisional carrier modes.
//!
//! The prepared mass-normalized drag owns rates nu_j and projected forcing f_j.
//! Its exact discrete constitutive law is sigma=sum f_j^2/(nu_j+i*omega).
//! Under exp(+i*omega*t), the free-carrier permittivity is -i*sigma/omega.
//! Differentiating these original modes retains mutual Coulomb drag; no
//! independent-carrier substitute, fitted plasma pole or new eigensolve enters.
//!
//! Constitutive convention: R. Fitzpatrick, equations 1201--1208 (opposite
//! phasor), https://farside.ph.utexas.edu/teaching/em/lectures/node103.html.

use super::{IonizedGasHarmonic, finite, invalid};
use crate::ionized_gas::{IonizedGasConductivity, IonizedGasError, poll, positive};
use fs_math::{c64::C64, dd::Dd};

/// Local free-carrier derivatives with the complete frozen source retained.
/// A separately declared bound-charge background contributes no derivatives
/// here. Actual field/drift and harmonic momentum admission remain owned by
/// IonizedGasHarmonic::response at the frequency used by the consumer.
#[derive(Clone, Copy, Debug)]
pub struct IonizedGasFrequencyDerivatives<'a> {
    /// Original populations, temperatures, neutral rates and mutual drag.
    pub source: &'a IonizedGasConductivity,
    /// Angular frequency [rad/s], without a cyclic-frequency conversion.
    pub omega_rad_s: f64,
    /// [omega*epsilon', omega^2*epsilon''] [F/m] of free carriers alone.
    pub scaled_permittivity_derivatives_f_m: [C64; 2],
    /// [Im(epsilon+omega*epsilon'), Im(2*omega*epsilon'+omega^2*epsilon'')].
    /// These [F/m] combinations keep small attenuation slopes independent
    /// of cancellation against the much larger Ohmic low-frequency term.
    pub permittivity_loss_moments_f_m: [f64; 2],
}

impl<'a> IonizedGasHarmonic<'a> {
    /// Differentiate the original prepared relaxation spectrum at held gas
    /// state and collision data. O(number of modes) work, no allocation or
    /// additional factorization. Polls cancellation at each original mode.
    ///
    /// This derivative-only capability preserves ordinary harmonic response
    /// admission: consumers retain their already admitted native value and
    /// frequency. It supplies no changing populations, magnetized response,
    /// kinetic/nonlocal transport or broader source validity.
    pub fn permittivity_frequency_derivatives(
        &self,
        omega_rad_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IonizedGasFrequencyDerivatives<'a>, IonizedGasError> {
        poll(&mut cancelled)?;
        if !positive(omega_rad_s) {
            return Err(invalid("carrier dispersion requires positive finite angular frequency"));
        }
        let mut sums = [Dd::from_f64(0.0); 6];
        for (&force, &rate) in self.projected_forcing.iter().zip(&self.rates) {
            poll(&mut cancelled)?;
            let values = mode_derivatives(force, rate, omega_rad_s)?;
            for (sum, value) in sums.iter_mut().zip(values) {
                *sum = *sum + Dd::from_f64(value);
            }
        }
        let values = sums.map(Dd::to_f64);
        if !values.into_iter().all(f64::is_finite) {
            return Err(invalid("carrier frequency derivative sum is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(IonizedGasFrequencyDerivatives {
            source: self.source,
            omega_rad_s,
            scaled_permittivity_derivatives_f_m: [
                C64::new(values[0], values[1]), C64::new(values[2], values[3]),
            ],
            permittivity_loss_moments_f_m: [values[4], values[5]],
        })
    }
}

// For chi=-f^2/[omega*(omega-i*nu)], put u=omega^2/(nu^2+omega^2).
// Re(omega*chi')=-2*u*Re(chi), Im(omega*chi')=-Im(chi)-2*u*Im(chi).
// The second derivative and both loss moments are polynomials in u.
// Evaluate u*Im(chi)=Re(chi)*nu*omega/(nu^2+omega^2) directly: forming u
// first can underflow even when this small absorption combination is finite.
fn mode_derivatives(force: f64, rate: f64, omega: f64) -> Result<[f64; 6], IonizedGasError> {
    if force == 0.0 {
        return Ok([0.0; 6]);
    }
    let scale = rate.max(omega);
    let weighted = force / scale;
    let x = rate.min(omega) / scale;
    let denominator = 1.0 + x * x;
    let real = -(weighted / denominator) * weighted;
    // Keep the native conductivity division for the dissipative component;
    // dividing its real part by omega does not cancel constitutive terms.
    let conductivity = (C64::from_re(force) / C64::new(rate, omega)).scale(force);
    let imaginary = -conductivity.re / omega;
    let cross = x / denominator;
    let u = if omega <= rate { (x * x) / denominator } else { 1.0 / denominator };
    let u_real = if omega <= rate { (real * x) * cross } else { real * u };
    let u_imaginary = real * cross;
    let first = C64::new(-2.0 * u_real, -imaginary - 2.0 * u_imaginary);
    let second = C64::new((8.0 * u - 2.0) * u_real,
        2.0 * imaginary + (2.0 + 8.0 * u) * u_imaginary);
    let moments = [-2.0 * u_imaginary, (8.0 * u - 2.0) * u_imaginary];
    if !finite(conductivity) || !real.is_finite() || !imaginary.is_finite()
        || !finite(first) || !finite(second) || !moments.into_iter().all(f64::is_finite)
        || (weighted == 0.0 && force != 0.0)
    {
        return Err(invalid("carrier frequency derivatives are unrepresentable"));
    }
    Ok([first.re, first.im, second.re, second.im, moments[0], moments[1]])
}

#[cfg(test)]
mod tests;
