//! Local Fresnel reflection from inside a passive absorbing medium.
//!
//! The incident field is specified at z=0, with a real conserved kx. An
//! absorbing incident wave has a complex normal wave number; assigning it a
//! real refraction angle would generally invent transverse attenuation.
//! Incident and reflected powers alone are not additive at the interface:
//! their interference contributes a real normal Poynting flux. All four terms
//! are returned in W/m² for the caller's actual peak electric field.
//!
//! The exp(+i omega t), exp(-i k.r) and fixed tangential-E conventions are the
//! parent's. Fields are local boundary values supplied by a finite upstream
//! source; no illumination at infinite depth or launch-distance attenuation
//! is inferred. This stationary, isotropic single-interface owner does not
//! change the lossless-port admission of films or angular spectra.
//!
//! Poynting interference and absorbing incidence: S. J. Byrnes, "Multilayer
//! optical calculations", §§4.3–4.4.1 and Appendix B, arXiv:1603.02720v5,
//! <https://arxiv.org/pdf/1603.02720>. That reference uses the conjugate time
//! convention; here both source permeabilities remain explicitly supplied.

use super::{DielectricResponse, PlanarDielectricInterface, ReflectionError, positive};
use fs_math::c64::C64;

mod kernel;
pub use kernel::{AbsorbingInterfaceReflection, InterfaceBoundaryField};

impl PlanarDielectricInterface<'_> {
    /// Reflect the actual complex peak [Ex,Ey] [V/m] at z=0 with signed real
    /// kx [rad/m], ky=0. Ex is TM and Ey is TE; arbitrary elliptical mixtures
    /// retain their relative phase. The incident source may absorb.
    ///
    /// Returned normal powers are signed along +z, including the independent
    /// incident/reflected interference term. They are not R/T fractions or
    /// interface heat. The incoming wave must carry positive unit-amplitude
    /// normal power; pure evanescent/grazing incidence refuses. Zero authored
    /// illumination returns zero fields and power after the same admission.
    /// Evaluation is fixed-size and allocation-free after source preparation.
    pub fn response_at_transverse_wavenumber(
        &self,
        omega_rad_s: f64,
        transverse_wavenumber_rad_m: f64,
        incident_tangential_e_v_m: [C64; 2],
    ) -> Result<AbsorbingInterfaceReflection, ReflectionError> {
        let omega = positive(omega_rad_s)?;
        // Plasma's actual fields require its original drift-amplitude gate.
        // This local two-wave incident state has no such admission owner yet.
        if matches!(self.incident_response, DielectricResponse::IonizedGas(_))
            || matches!(self.transmitted, DielectricResponse::IonizedGas(_))
        {
            return Err(ReflectionError::Invalid(
                "absorbing-incidence fields require original plasma drift admission",
            ));
        }
        kernel::solve(kernel::Input {
            omega,
            kx: transverse_wavenumber_rad_m,
            epsilon_i: self.incident_response.permittivity(omega)?,
            epsilon_t: self.transmitted.permittivity(omega)?,
            mu_i: self.incident_mu_h_m,
            mu_t: self.transmitted_mu_h_m,
            electric: incident_tangential_e_v_m,
        })
    }
}
