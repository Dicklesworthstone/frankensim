//! Density jets from the original IF97 Region 1/2 Gibbs equations.
//!
//! These optional derivatives support spatial constitutive laws driven by
//! accepted pressure and temperature fields. No fit, finite difference, phase
//! interpolation or new EOS is introduced. Existing state evaluation and its
//! operation tree remain unchanged; the isobaric path retains its smaller cost.
//!
//! With `q = pi * gamma_pi`, `v = R*T*q/p` and `tau = T_star/T`,
//! `v_TT/v = (tau/T)^2 * (pi*gamma_pi_tau_tau)/q`. Differentiating
//! `rho = 1/v` gives `rho_T = -rho*v_T/v` and
//! `rho_TT = rho*(2*(v_T/v)^2 - v_TT/v)`. The third Gibbs partial is evaluated
//! directly from the original coefficients only when this jet is requested.
//! The full jet also retains pressure curvature and the mixed derivative.
//! Region 2 pressure curvature is evaluated through the residual derivatives
//! of `q`, avoiding cancellation of the ideal-gas contribution as `p -> 0`.
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
/// `(d rho / d p)_T`, in kg/(m3 Pa).
pub type DensityPressureDerivative = Qty<-2, 0, 2, 0, 0, 0>;
/// `(d² rho / d p²)_T`, in kg/(m3 Pa²).
pub type DensityPressureSecondDerivative = Qty<-1, -1, 4, 0, 0, 0>;
/// `d² rho / (d p d T)`, in kg/(m3 Pa K).
pub type DensityPressureTemperatureDerivative = Qty<-2, 0, 2, -1, 0, 0>;

/// Value and analytic temperature derivatives at the original fixed pressure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsobaricDensityJet {
    pub density: Density,
    pub temperature_derivative: DensityTemperatureDerivative,
    pub temperature_second_derivative: DensityTemperatureSecondDerivative,
}

/// Value and complete analytic Hessian in the original `(p, T)` coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DensityJet {
    pub density: Density,
    pub temperature_derivative: DensityTemperatureDerivative,
    pub pressure_derivative: DensityPressureDerivative,
    pub temperature_second_derivative: DensityTemperatureSecondDerivative,
    pub pressure_second_derivative: DensityPressureSecondDerivative,
    pub pressure_temperature_derivative: DensityPressureTemperatureDerivative,
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
        let (tau, qtt) = match self.region() {
            If97Region::Region1 => {
                let tau = R1_TEMPERATURE_STAR_K / t;
                (
                    tau,
                    region1_thirds::<false>(p / R1_PRESSURE_STAR_PA, tau)[0],
                )
            }
            If97Region::Region2 => {
                let tau = R2_TEMPERATURE_STAR_K / t;
                (
                    tau,
                    region2_thirds::<false>(p / R2_PRESSURE_STAR_PA, tau)[0],
                )
            }
        };
        self.temperature_density_jet(tau, qtt)
    }

    /// Evaluate all pressure/temperature density derivatives without allocation.
    ///
    /// Value and first derivatives reuse this admitted Gibbs state. The extra
    /// third partials share one set of fixed-size power tables and compensated
    /// sums. This is a tangent to the selected single-phase branch, including
    /// at saturation; no mixture or automatic regional continuation is implied.
    pub fn density_jet(self) -> Result<DensityJet, If97Error> {
        let t = self.temperature().value();
        let p = self.pressure().value();
        let v = self.specific_volume().value();
        let rho = self.density().value();
        let q = (v * p) / (R_SPECIFIC_J_KG_K * t);
        let alpha = self.thermal_expansion().value();
        let pressure_slope = self.dv_dp_at_temperature().value() / v;
        let (tau, qtt, ppt, pressure_second) = match self.region() {
            If97Region::Region1 => {
                let tau = R1_TEMPERATURE_STAR_K / t;
                let [qtt, ppt, ppp] = region1_thirds::<true>(p / R1_PRESSURE_STAR_PA, tau);
                let v_pp_over_v = ((ppp / q) / p) / p;
                let second = rho * (2.0 * pressure_slope).mul_add(pressure_slope, -v_pp_over_v);
                (tau, qtt, ppt, second)
            }
            If97Region::Region2 => {
                let tau = R2_TEMPERATURE_STAR_K / t;
                let [qtt, ppt, q_pi, q_pipi] = region2_thirds::<true>(p / R2_PRESSURE_STAR_PA, tau);
                // rho = p/(R*T*q). Differentiate q's residual directly:
                // rho_pp = (rho/p) * ((2*p*a - 2)*a - p*b),
                // a = q_p/q, b = q_pp/q. Neither residual contains log(pi),
                // so dilute vapor retains its small, finite virial curvature.
                let a = (q_pi / R2_PRESSURE_STAR_PA) / q;
                let b = ((q_pipi / R2_PRESSURE_STAR_PA) / R2_PRESSURE_STAR_PA) / q;
                let second = (rho / p) * ((2.0 * p).mul_add(a, -2.0)).mul_add(a, -p * b);
                (tau, qtt, ppt, second)
            }
        };
        let thermal = self.temperature_density_jet(tau, qtt)?;
        // v_pT/v = (pi²*gamma_pipi - tau*pi²*gamma_pipitau)/(p*T*q).
        // Scale by p before combining, retaining the finite ideal-gas limit.
        let beta = pressure_slope * p;
        let mixed = (rho / p) * beta.mul_add(2.0 * alpha - 1.0 / t, (tau / t) * (ppt / q));
        let first = rho * self.isothermal_compressibility().value();
        if ![first, pressure_second, mixed]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(If97Error::Unrepresentable(
                "pressure/temperature density derivatives",
            ));
        }
        Ok(DensityJet {
            density: thermal.density,
            temperature_derivative: thermal.temperature_derivative,
            pressure_derivative: DensityPressureDerivative::new(first),
            temperature_second_derivative: thermal.temperature_second_derivative,
            pressure_second_derivative: DensityPressureSecondDerivative::new(pressure_second),
            pressure_temperature_derivative: DensityPressureTemperatureDerivative::new(mixed),
        })
    }

    fn temperature_density_jet(self, tau: f64, qtt: f64) -> Result<IsobaricDensityJet, If97Error> {
        let t = self.temperature().value();
        let p = self.pressure().value();
        let v = self.specific_volume().value();
        let rho = self.density().value();
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
fn region1_thirds<const PRESSURE: bool>(pi: f64, tau: f64) -> [f64; 3] {
    let x = powers::<33>(7.1 - pi, 0);
    let y = powers::<61>(tau - 1.222, 43);
    let mut sums = [Sum::default(); 3];
    for &GibbsTerm { i, j, n } in &R1_TERMS {
        if i != 0 && j != 0 && j != 1 {
            let coefficient = -n * f64::from(i) * f64::from(j) * f64::from(j - 1);
            sums[0].add(coefficient * x[(i - 1) as usize] * y[(j + 41) as usize]);
        }
        if PRESSURE && i >= 2 && j != 0 {
            let coefficient = n * f64::from(i) * f64::from(i - 1) * f64::from(j);
            sums[1].add(coefficient * x[(i - 2) as usize] * y[(j + 42) as usize]);
        }
        if PRESSURE && i >= 3 {
            let coefficient = -n * f64::from(i) * f64::from(i - 1) * f64::from(i - 2);
            sums[2].add(coefficient * x[(i - 3) as usize] * y[(j + 43) as usize]);
        }
    }
    [
        pi * sums[0].value,
        pi * pi * sums[1].value,
        pi * pi * pi * sums[2].value,
    ]
}

fn region2_thirds<const PRESSURE: bool>(pi: f64, tau: f64) -> [f64; 4] {
    let x = powers::<25>(pi, 0);
    let y = powers::<59>(tau - 0.5, 0);
    let mut sums = [Sum::default(); 4];
    for &GibbsTerm { i, j, n } in &R2_RESIDUAL {
        if j >= 2 {
            let coefficient = n * f64::from(i) * f64::from(j) * f64::from(j - 1);
            sums[0].add(coefficient * x[i as usize] * y[(j - 2) as usize]);
        }
        if PRESSURE {
            if i >= 2 && j != 0 {
                let coefficient = n * f64::from(i) * f64::from(i - 1) * f64::from(j);
                sums[1].add(coefficient * x[i as usize] * y[(j - 1) as usize]);
            }
            let coefficient = n * f64::from(i * i);
            sums[2].add(coefficient * x[(i - 1) as usize] * y[j as usize]);
            if i >= 2 {
                sums[3].add(coefficient * f64::from(i - 1) * x[(i - 2) as usize] * y[j as usize]);
            }
        }
    }
    [sums[0].value, sums[1].value, sums[2].value, sums[3].value]
}

#[cfg(test)]
mod tests;
