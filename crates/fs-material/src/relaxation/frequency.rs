//! Analytic frequency derivative of the original molecular acoustic branch.
//!
//! The equilibrium gas, composition and evaluated reservoir times are held.
//! Both parts of dk/domega matter: its real part is group slowness, while its
//! imaginary part is the frequency slope of molecular attenuation. Neither
//! may be replaced by an equilibrium sound speed or a separate loss fit.

use super::{RelaxationState, debye, parameters, positive, state_error};
use crate::MaterialError;
use fs_math::c64::C64;

/// Molecular-only propagation for exp(i(k*x - omega*t)). Group slowness is a
/// phase derivative, not an energy velocity or a wavefront-speed certificate.
/// Shear viscosity, conduction and wall layers remain separate owners.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularFrequencyResponse {
    /// Original relaxing wavenumber [1/m], including positive attenuation.
    pub wavenumber_per_m: C64,
    /// Analytic dk/domega [s/m]. The real part is group slowness; the imaginary
    /// part is d(alpha)/domega. No positivity clipping is applied to either.
    pub wavenumber_derivative_s_m: C64,
}

impl RelaxationState {
    /// Original molecular wavenumber and its analytic angular-frequency slope.
    ///
    /// For d_j=1/(1-i*omega*tau_j), omega*cv'=sum C_j*(d_j^2-d_j) and
    /// k'=sqrt(rho/K)*(1+R*omega*cv'/(2*cv*cp)). Evaluation is allocation-free
    /// and bounded by the original reservoir count. The value retains the
    /// exact arithmetic of [`Self::relaxing_wavenumber`].
    ///
    /// The derivative is evaluated as the frozen slowness plus an analytic
    /// excess. Writing I=sum C_j*d_j, J=sum C_j*d_j^2, u=I/cp,
    /// a=(R/cv_inf)*u, t=sqrt(1+a), r=a/(t+1), the excess is
    /// `((R/cv_inf)*(u*u+(1-u)*J/cp)+r*r)/(2*t)`.
    /// This is the same derivative without subtracting two O(1/omega) terms
    /// to recover the O(1/omega^3) attenuation slope in the frozen tail.
    ///
    /// # Errors
    /// Refuses nonpositive/nonfinite omega or an unrepresentable response.
    pub fn relaxing_frequency_response(
        &self,
        omega: f64,
    ) -> Result<MolecularFrequencyResponse, MaterialError> {
        if !positive(omega) {
            return Err(parameters(
                "molecular frequency response requires positive finite omega",
            ));
        }
        let caloric = self.caloric_response(omega)?;
        let slowness = caloric
            .bulk_modulus_pa
            .scale(1.0 / self.gas.density)
            .sqrt()
            .recip();
        let wavenumber = slowness.scale(omega);
        let mut increment = C64::ZERO;
        let mut squared = C64::ZERO;
        for process in &self.processes {
            let (re, im) = debye(omega * process.relaxation_time_s);
            let d = C64::new(re, im);
            increment = increment + d.scale(process.heat_capacity_j_kg_k);
            squared = squared + (d * d).scale(process.heat_capacity_j_kg_k);
        }
        let ratio = self.gas.specific_gas_constant / self.cv_frozen;
        let u = increment / caloric.cp_j_kg_k;
        let a = u.scale(ratio);
        let t = (C64::ONE + a).sqrt();
        let r = a / (t + C64::ONE);
        let summed = (u * u + (C64::ONE - u) * (squared / caloric.cp_j_kg_k)).scale(ratio);
        let frozen_slowness = (C64::ONE
            + C64::from_re(self.cv_frozen)
                .recip()
                .scale(self.gas.specific_gas_constant))
        .scale(self.gas.pressure)
        .scale(1.0 / self.gas.density)
        .sqrt()
        .recip();
        let derivative = frozen_slowness * (C64::ONE + ((summed + r * r) / t).scale(0.5));
        if !positive(wavenumber.re)
            || !wavenumber.im.is_finite()
            || wavenumber.im < 0.0
            || !derivative.re.is_finite()
            || !derivative.im.is_finite()
        {
            return Err(state_error(
                "molecular frequency response is unrepresentable",
            ));
        }
        Ok(MolecularFrequencyResponse {
            wavenumber_per_m: wavenumber,
            wavenumber_derivative_s_m: derivative,
        })
    }
}

#[cfg(test)]
mod tests;
