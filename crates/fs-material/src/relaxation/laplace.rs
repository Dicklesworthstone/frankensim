//! Right-half-plane response of the ORIGINAL molecular caloric reservoirs.
//!
//! `cv(s)=cv_inf+sum Cj/(1+s*tauj)`, `K(s)=p(1+R/cv(s))`.
//! This is an analytic constitutive continuation, not an interpolation of
//! frequency samples or a continuation of the classical omega-squared loss.
//! The frozen-speed delay is separated without subtracting nearly equal
//! wavenumbers. Positive finite reservoirs give the bounded wavefront of a
//! completely monotone relaxation law; see Hanyga (2013), Wave Motion 50,
//! 909–928, <https://doi.org/10.1016/j.wavemoti.2013.03.002>.

use super::{RelaxationState, positive, state_error};
use crate::MaterialError;
use fs_math::c64::C64;

/// Molecular-only Laplace response, with transforms using exp(-s t).
/// On the boundary s=-i omega it matches the original exp(-i omega t)
/// harmonic constitutive and propagation convention. Shear viscosity,
/// conduction and wall layers remain separate physical mechanisms.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularLaplaceResponse {
    pub cv_j_kg_k: C64,
    pub cp_j_kg_k: C64,
    pub bulk_modulus_pa: C64,
    /// Principal causal slowness sqrt(rho/K), [s/m].
    pub slowness_s_m: C64,
    /// kappa(s)-s/c_inf, [1/m], evaluated without catastrophic subtraction.
    pub excess_propagation_rate_per_m: C64,
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

// A scaled Debye reciprocal avoids overflowing s*tau, including on the
// imaginary boundary. For s=0 the equilibrated reservoir survives exactly.
fn debye(s: C64, tau: f64) -> C64 {
    let magnitude = s.re.abs().max(s.im.abs());
    if magnitude == 0.0 {
        return C64::ONE;
    }
    if magnitude <= tau.recip() {
        (C64::ONE + s.scale(tau)).recip()
    } else {
        let reciprocal = (1.0 / magnitude) / tau;
        C64::new(s.re / magnitude + reciprocal, s.im / magnitude)
            .recip()
            .scale(reciprocal)
    }
}

impl RelaxationState {
    /// Evaluate all original reservoirs at finite Re(s)>=0, including s=0.
    ///
    /// No state, rate, species or heat capacity is fitted or discarded. The
    /// ratio `a=R*(cv-cv_inf)/(cv_inf*(cv+R))` gives
    /// `c_inf*slowness=sqrt(1+a)` and
    /// `kappa-s/c_inf=(s/c_inf)*a/(sqrt(1+a)+1)`. Retaining the small excess
    /// explicitly preserves weak absorption after the exact wavefront delay.
    /// The density/EOS and frozen sound speed remain this owner's originals.
    pub fn laplace_response(&self, s: C64) -> Result<MolecularLaplaceResponse, MaterialError> {
        if !finite(s) || s.re < 0.0 {
            return Err(state_error(
                "molecular Laplace response requires finite Re(s)>=0",
            ));
        }
        let gas = self.gas();
        let frozen = self.frozen_cv_j_kg_k();
        let mut increment = C64::ZERO;
        for process in self.processes() {
            increment =
                increment + debye(s, process.relaxation_time_s).scale(process.heat_capacity_j_kg_k);
        }
        let cv = C64::from_re(frozen) + increment;
        let cp = cv + C64::from_re(gas.specific_gas_constant);
        let bulk = (C64::ONE + cv.recip().scale(gas.specific_gas_constant)).scale(gas.pressure);
        let ratio = (increment / cp).scale(gas.specific_gas_constant / frozen);
        let relative_slowness = (C64::ONE + ratio).sqrt();
        let frozen_slowness = self.frozen_sound_speed_m_s().recip();
        let slowness = relative_slowness.scale(frozen_slowness);
        let excess = s * (ratio / (relative_slowness + C64::ONE)).scale(frozen_slowness);
        if !positive(cv.re)
            || !positive(cp.re)
            || !positive(bulk.re)
            || !positive(slowness.re)
            || !positive(frozen_slowness)
            || [cv, cp, bulk, slowness, excess]
                .into_iter()
                .any(|z| !finite(z))
            || excess.re < 0.0
        {
            return Err(state_error("molecular Laplace response is unrepresentable"));
        }
        Ok(MolecularLaplaceResponse {
            cv_j_kg_k: cv,
            cp_j_kg_k: cp,
            bulk_modulus_pa: bulk,
            slowness_s_m: slowness,
            excess_propagation_rate_per_m: excess,
        })
    }
}

#[cfg(test)]
mod tests;
