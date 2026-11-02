//! Opaque surface absorption from the original electromagnetic boundary flux.
//!
//! The reflected and absorbed fractions share the native complex-index and
//! coherent-film owners. Absorption uses normal Poynting flux, not `1-R`, so
//! a representable weak loss survives a reflectance that rounds to one.

use super::{
    CoatingLayer, CoherentCoating, CoherentCoatingStack, Complex, ConductorError, ConductorOptics,
};

/// Fractions of incident ray power at a smooth opaque surface, ordered s, p.
/// Each polarization closes reflected + absorbed = 1 up to roundoff. There is
/// no transmitted scene ray: absorption includes every coating and the opaque
/// half space together, with no inferred layer or penetration-depth partition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductorPowerBalance {
    pub reflected: [f64; 2],
    pub absorbed: [f64; 2],
}

impl ConductorOptics {
    /// Resolve reflected and absorbed fractions from the original scaled
    /// boundary fields, including the explicitly attached coherent coating.
    ///
    /// Optical support, incidence, source identities and local planar/coherence
    /// assumptions are unchanged. Substrate eta and k are strictly positive.
    /// Exact grazing has zero absorbed fraction. A numerically unresolved net
    /// boundary flux refuses rather than publishing a negative or invented heat.
    pub fn power_balance(
        &self,
        wavelength_nm: f64,
        incident_eta: f64,
        incidence_cosine: f64,
    ) -> Result<ConductorPowerBalance, ConductorError> {
        match &self.coating {
            Some(CoherentCoating::Dielectric(film)) => {
                CoherentCoatingStack::try_new(&[CoatingLayer::Dielectric(film)])?.power_balance(
                    self,
                    wavelength_nm,
                    incident_eta,
                    incidence_cosine,
                )
            }
            Some(CoherentCoating::Absorbing(film)) => {
                CoherentCoatingStack::try_new(&[CoatingLayer::Absorbing(film)])?.power_balance(
                    self,
                    wavelength_nm,
                    incident_eta,
                    incidence_cosine,
                )
            }
            Some(CoherentCoating::Periodic(coating)) => coating.power_balance(
                wavelength_nm,
                incident_eta,
                self.complex_ior(wavelength_nm)?,
                incidence_cosine,
            ),
            None => CoherentCoatingStack::try_new(&[])?.power_balance(
                self,
                wavelength_nm,
                incident_eta,
                incidence_cosine,
            ),
        }
    }
}

impl ConductorPowerBalance {
    pub(super) fn from_boundary(
        cosine: f64,
        incident_admittance: [f64; 2],
        field: [Complex; 2],
        flux: [Complex; 2],
        amplitudes: [Complex; 2],
    ) -> Result<Self, ConductorError> {
        if cosine == 0.0 {
            return Ok(Self {
                reflected: [1.0; 2],
                absorbed: [0.0; 2],
            });
        }
        let mut result = Self {
            reflected: [0.0; 2],
            absorbed: [0.0; 2],
        };
        for pol in 0..2 {
            let f = field[pol];
            let g = flux[pol];
            let y = incident_admittance[pol];
            // Re(F conj(G)). Two product residuals retain the representable
            // part of a cancellation without constructing complex power.
            let a = f.re * g.re;
            let b = f.im * g.im;
            let normal_flux = (a + b) + (f.re.mul_add(g.re, -a) + f.im.mul_add(g.im, -b));
            let flux_resolution = 128.0 * f64::EPSILON * (a.abs() + b.abs());
            let denominator = (f * y + g).norm_squared();
            if !normal_flux.is_finite()
                || normal_flux <= flux_resolution
                || !denominator.is_finite()
                || denominator <= 0.0
                || !y.is_finite()
                || y <= 0.0
            {
                return Err(ConductorError::NumericalFailure);
            }
            // r=(Y F-G)/(Y F+G) gives A=4Y Re(F conj(G))/|Y F+G|².
            // A common field/flux normalization cancels, so no removed
            // propagation scale or subtraction from rounded reflectance enters.
            let absorbed = 4.0 * y * (normal_flux / denominator);
            let reflected = amplitudes[pol].norm_squared();
            let allowance = 512.0 * f64::EPSILON;
            if !absorbed.is_finite()
                || absorbed <= 0.0
                || absorbed > 1.0 + allowance
                || !reflected.is_finite()
                || !(0.0..=1.0 + allowance).contains(&reflected)
                || ((absorbed + reflected) - 1.0).abs() > allowance
            {
                return Err(ConductorError::NumericalFailure);
            }
            result.absorbed[pol] = absorbed.min(1.0);
            result.reflected[pol] = reflected.min(1.0);
        }
        if cosine == 1.0 {
            result.absorbed[1] = result.absorbed[0];
            result.reflected[1] = result.reflected[0];
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
