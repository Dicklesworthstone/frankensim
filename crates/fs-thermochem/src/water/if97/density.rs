//! Isobaric density jets from the original IF97 Region 1/2 Gibbs equations.
//!
//! These optional derivatives support spatial constitutive laws driven by an
//! accepted temperature field at a declared fixed pressure. No fit, finite
//! difference, phase interpolation or new EOS is introduced. Existing state
//! evaluation and its operation tree remain unchanged.
//!
//! With `q = pi * gamma_pi`, `v = R*T*q/p` and `tau = T_star/T`,
//! `v_TT/v = (tau/T)^2 * (pi*gamma_pi_tau_tau)/q`. Differentiating
//! `rho = 1/v` gives `rho_T = -rho*v_T/v` and
//! `rho_TT = rho*(2*(v_T/v)^2 - v_TT/v)`. The third Gibbs partial is evaluated
//! directly from the original coefficients only when this jet is requested.
//! Source: IAPWS R7-97(2012), equations (7), (15)--(17), Tables 2 and 11:
//! <https://iapws.org/technical-guidance/release/IF97-Rev.download>.
//!
//! Saturation endpoints retain the selected single-phase branch derivative;
//! this is not a derivative of a two-phase mixture or a crossing of Region 4.

use super::{
    GibbsTerm, If97Error, If97Region, If97State, R_SPECIFIC_J_KG_K, R1_PRESSURE_STAR_PA,
    R1_TEMPERATURE_STAR_K, R1_TERMS, R2_PRESSURE_STAR_PA, R2_RESIDUAL, R2_TEMPERATURE_STAR_K, Sum,
    powers,
};
use fs_qty::{Density, Qty};

/// `(d rho / d T)_p`, in kg/(m3 K).
pub type DensityTemperatureDerivative = Qty<-3, 1, 0, -1, 0, 0>;
/// `(d² rho / d T²)_p`, in kg/(m3 K²).
pub type DensityTemperatureSecondDerivative = Qty<-3, 1, 0, -2, 0, 0>;

/// Value and analytic temperature derivatives at the original fixed pressure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsobaricDensityJet {
    pub density: Density,
    pub temperature_derivative: DensityTemperatureDerivative,
    pub temperature_second_derivative: DensityTemperatureSecondDerivative,
}

impl If97State {
    /// Evaluate the density jet of this exact original state without allocation.
    ///
    /// The value and first derivative reuse the stored Gibbs state. Only the
    /// additional mixed third partial is summed, with the native compensated
    /// accumulator. This method does not broaden the state's source domain.
    pub fn isobaric_density_jet(self) -> Result<IsobaricDensityJet, If97Error> {
        let t = self.temperature().value();
        let p = self.pressure().value();
        let v = self.specific_volume().value();
        let rho = self.density().value();
        let (tau, qtt) = match self.region() {
            If97Region::Region1 => {
                let tau = R1_TEMPERATURE_STAR_K / t;
                (tau, region1_mixed_third(p / R1_PRESSURE_STAR_PA, tau))
            }
            If97Region::Region2 => {
                let tau = R2_TEMPERATURE_STAR_K / t;
                (tau, region2_mixed_third(p / R2_PRESSURE_STAR_PA, tau))
            }
        };
        let q = (v * p) / (R_SPECIFIC_J_KG_K * t);
        let slope = self.dv_dt_at_pressure().value() / v;
        let scaled_tau = tau / t;
        let curvature = (qtt / q) * scaled_tau * scaled_tau;
        let first = -rho * slope;
        let second = rho * (2.0 * slope).mul_add(slope, -curvature);
        if ![q, curvature, first, second]
            .into_iter()
            .all(f64::is_finite)
            || q <= 0.0
        {
            return Err(If97Error::Unrepresentable("isobaric density derivatives"));
        }
        Ok(IsobaricDensityJet {
            density: self.density(),
            temperature_derivative: DensityTemperatureDerivative::new(first),
            temperature_second_derivative: DensityTemperatureSecondDerivative::new(second),
        })
    }
}

// Retain the original scaled-pressure convention to avoid 1/pi intermediates
// in the vapor limit. The ideal Region 2 log(pi) has zero mixed third partial.
fn region1_mixed_third(pi: f64, tau: f64) -> f64 {
    let x = powers::<33>(7.1 - pi, 0);
    let y = powers::<61>(tau - 1.222, 43);
    let mut sum = Sum::default();
    for &GibbsTerm { i, j, n } in &R1_TERMS {
        if i == 0 || j == 0 || j == 1 {
            continue;
        }
        let coefficient = -n * f64::from(i) * f64::from(j) * f64::from(j - 1);
        sum.add(coefficient * x[(i - 1) as usize] * y[(j + 41) as usize]);
    }
    pi * sum.value
}

fn region2_mixed_third(pi: f64, tau: f64) -> f64 {
    let x = powers::<25>(pi, 0);
    let y = powers::<59>(tau - 0.5, 0);
    let mut sum = Sum::default();
    for &GibbsTerm { i, j, n } in &R2_RESIDUAL {
        if j < 2 {
            continue;
        }
        let coefficient = n * f64::from(i) * f64::from(j) * f64::from(j - 1);
        sum.add(coefficient * x[i as usize] * y[(j - 2) as usize]);
    }
    sum.value
}

#[cfg(test)]
mod tests;
