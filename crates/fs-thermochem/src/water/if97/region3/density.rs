//! Density gradients and curvature of the original Region 3 Helmholtz source.
//!
//! Implicit differentiation of `p(rho,T)` gives `rho_p = 1/p_rho`,
//! `rho_T = -p_T/p_rho` and the complete Hessian below. The pressure Hessian
//! is summed directly from Eq. (28), Table 30 of IAPWS R7-97(2012):
//! <https://iapws.org/technical-guidance/release/IF97-Rev.download>.
//! No density inversion, finite difference, branch switch or regularization is
//! performed. The admitted state's source domain and uncertainty are unchanged.

use super::super::density::{
    DensityJet, DensityPressureDerivative, DensityPressureSecondDerivative,
    DensityPressureTemperatureDerivative, DensityTemperatureDerivative,
    DensityTemperatureSecondDerivative, IsobaricDensityJet,
};
use super::super::{Sum, powers};
use super::{
    CRITICAL_DENSITY_KG_M3, CRITICAL_TEMPERATURE_K, If97DenseError, If97DenseState,
    LOG_COEFFICIENT, R_SPECIFIC_J_KG_K, TERMS,
};

impl If97DenseState {
    /// Complete density Hessian in `(p,T)` coordinates at this admitted state.
    ///
    /// First responses retain the original state stiffness and thermal pressure
    /// derivative. Three additional compensated sums over fixed-size power
    /// tables supply the analytic pressure Hessian without allocation. Close to
    /// the critical point the physical response may be very large; unrepresentable
    /// responses refuse instead of being clipped to a finite value.
    pub fn density_jet(self) -> Result<DensityJet, If97DenseError> {
        let rho = self.density().value();
        let t = self.temperature().value();
        let stiffness = self.dp_drho_at_temperature().value();
        let rho_p = stiffness.recip();
        let rho_t = -self.dp_dt_at_density().value() / stiffness;
        let [p_rr, p_rt, p_tt] = pressure_hessian(rho, t);
        // Divide successively, rather than forming stiffness^3, which can
        // overflow even when its reciprocal response remains representable.
        let rho_pp = -((p_rr / stiffness) / stiffness) / stiffness;
        let rho_pt = -((p_rt + p_rr * rho_t) / stiffness) / stiffness;
        let rho_tt = -(p_tt + rho_t * (2.0 * p_rt + p_rr * rho_t)) / stiffness;
        if ![rho_p, rho_t, rho_pp, rho_pt, rho_tt]
            .into_iter()
            .all(f64::is_finite)
            || rho_p <= 0.0
        {
            return Err(If97DenseError::Unrepresentable("density response Hessian"));
        }
        Ok(DensityJet {
            density: self.density(),
            temperature_derivative: DensityTemperatureDerivative::new(rho_t),
            pressure_derivative: DensityPressureDerivative::new(rho_p),
            temperature_second_derivative: DensityTemperatureSecondDerivative::new(rho_tt),
            pressure_second_derivative: DensityPressureSecondDerivative::new(rho_pp),
            pressure_temperature_derivative: DensityPressureTemperatureDerivative::new(rho_pt),
        })
    }

    /// Isobaric density value, slope and curvature from the same native source.
    ///
    /// A dense source is expressed in `(rho,T)`, so even the isobaric curvature
    /// requires all three pressure Hessian entries. It performs the same single
    /// coefficient traversal as [`Self::density_jet`].
    pub fn isobaric_density_jet(self) -> Result<IsobaricDensityJet, If97DenseError> {
        let full = self.density_jet()?;
        Ok(IsobaricDensityJet {
            density: full.density,
            temperature_derivative: full.temperature_derivative,
            temperature_second_derivative: full.temperature_second_derivative,
        })
    }
}

// For one Helmholtz monomial q=n*delta^i*tau^j, p=R*T*rho*i*q.
// Differentiating that original monomial directly avoids cancellation between
// separately evaluated first, second and third Helmholtz derivatives. The log
// term contributes only R*n_log to the mixed pressure derivative.
fn pressure_hessian(rho: f64, t: f64) -> [f64; 3] {
    let density = powers::<12>(rho / CRITICAL_DENSITY_KG_M3, 0);
    let temperature = powers::<27>(CRITICAL_TEMPERATURE_K / t, 0);
    let mut sums = [Sum::default(); 3];
    sums[1].add(LOG_COEFFICIENT);
    for &(i, j, n) in &TERMS {
        if i == 0 {
            continue;
        }
        let term = n * density[i] * temperature[j];
        let (i, j) = (i as f64, j as f64);
        sums[0].add(i * i * (i + 1.0) * term);
        sums[1].add(i * (i + 1.0) * (1.0 - j) * term);
        sums[2].add(i * j * (j - 1.0) * term);
    }
    [
        (R_SPECIFIC_J_KG_K * t / rho) * sums[0].value,
        R_SPECIFIC_J_KG_K * sums[1].value,
        (R_SPECIFIC_J_KG_K * rho / t) * sums[2].value,
    ]
}

#[cfg(test)]
mod tests;
