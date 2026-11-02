//! Opaque surface absorption from the original electromagnetic boundary flux.
//!
//! The reflected and absorbed fractions share the native complex-index and
//! coherent-film owners. Absorption uses normal Poynting flux, not `1-R`, so
//! a representable weak loss survives a reflectance that rounds to one.

use super::{
    CoatingLayer, CoherentCoating, CoherentCoatingStack, Complex, ConductorError, ConductorOptics,
    MAX_COHERENT_STACK_LAYERS,
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

/// Positive absorbed fractions resolved over the original physical films and
/// opaque substrate. Fractions are ordered s, p and refer to incident ray power.
/// The substrate receives the flux crossing its first face; no penetration
/// depth or thermal control volume is inferred from that semi-infinite law.
#[derive(Clone, Debug, PartialEq)]
pub struct ConductorLayerPowerBalance {
    /// Original reflected and aggregate absorbed fractions.
    pub total: ConductorPowerBalance,
    /// Absorption in the opaque half space after all finite coating layers.
    pub substrate_absorbed: [f64; 2],
    pub(super) layers: [[f64; 2]; MAX_COHERENT_STACK_LAYERS],
    pub(super) layer_count: usize,
}

impl ConductorLayerPowerBalance {
    /// Absorption in authored incident-to-substrate order, including explicit
    /// zero-thickness and lossless films. A periodic coating expands its repeats
    /// in that same order. Entries plus substrate absorption equal the original
    /// aggregate absorption to floating-point rounding.
    #[must_use]
    pub fn layer_absorbed(&self) -> &[[f64; 2]] {
        &self.layers[..self.layer_count]
    }
}

impl ConductorOptics {
    /// Resolve where coherent surface absorption occurs before assigning heat
    /// to thermal cells. The original field solve retains interference, local
    /// incidence, spectral laws and source support. Positive dissipative modes
    /// preserve thin-film loss without subtracting nearly equal power fluxes.
    ///
    /// Work/storage are bounded by 64 physical layers, with no heap allocation.
    /// Existing scalar and Mueller reflection paths perform no extra work.
    pub fn layer_power_balance(
        &self,
        wavelength_nm: f64,
        incident_eta: f64,
        incidence_cosine: f64,
    ) -> Result<ConductorLayerPowerBalance, ConductorError> {
        match &self.coating {
            Some(CoherentCoating::Dielectric(film)) => {
                CoherentCoatingStack::try_new(&[CoatingLayer::Dielectric(film)])?
                    .layer_power_balance(self, wavelength_nm, incident_eta, incidence_cosine)
            }
            Some(CoherentCoating::Absorbing(film)) => {
                CoherentCoatingStack::try_new(&[CoatingLayer::Absorbing(film)])?
                    .layer_power_balance(self, wavelength_nm, incident_eta, incidence_cosine)
            }
            Some(CoherentCoating::Periodic(coating)) => coating.layer_power_balance(
                wavelength_nm,
                incident_eta,
                self.complex_ior(wavelength_nm)?,
                incidence_cosine,
            ),
            None => CoherentCoatingStack::try_new(&[])?.layer_power_balance(
                self,
                wavelength_nm,
                incident_eta,
                incidence_cosine,
            ),
        }
    }

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
