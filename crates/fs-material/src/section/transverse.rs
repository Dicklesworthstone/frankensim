//! Kelvin--Voigt transverse shear loss for the existing circular rod section.
//!
//! The original uniform-shear reduction supplies its declared correction kappa.
//! Applying that SAME reduction to an explicitly selected shear viscosity gives
//! `C_s = kappa sum(eta_s A)` [N s] in each transverse material direction.
//! Its dissipation per reference length is `C_s (gamma_dot_y²+gamma_dot_z²)`.
//! No axial/bending viscosity, fitted modal ratio or torsional polar moment is
//! substituted. This is the section's fixed-state small-strain approximation;
//! it does not resolve cross-sectional warping or interlayer slip.
//! The strain-rate Kelvin rod model follows Linn, Lang and Tuganov (2013),
//! <https://doi.org/10.5194/ms-4-79-2013>; layer integration and the declared
//! shear correction retain this repository's original section approximation.

use super::{ConcentricShearSection, MaterialStatePointError, invalid};
use fs_qty::{DynViscosity, QuantitySpec};

/// Original section with independently sourced transverse shear viscosity.
#[derive(Clone, Debug, PartialEq)]
pub struct ConcentricTransverseViscosity {
    shear: ConcentricShearSection,
    viscosity_property: String,
    layer_viscous_rigidity_n_s: Vec<f64>,
    viscous_rigidity_n_s: f64,
    omega_band_rad_s: (f64, f64),
}

impl ConcentricShearSection {
    /// Bind a nonnegative shear viscosity [Pa s] from every original layer.
    ///
    /// Density, shear modulus and viscosity must be constants at the retained
    /// environment. Viscosity needs an explicit frequency interval; all source
    /// intervals are intersected. Zero is supported, absent data is not zero.
    /// The original shear correction is applied exactly once to each eta*A.
    pub fn with_transverse_viscosity(
        self,
        property_name: &str,
    ) -> Result<ConcentricTransverseViscosity, MaterialStatePointError> {
        let mut band = (0.0_f64, f64::MAX);
        let mut total = 0.0;
        let mut layers = Vec::with_capacity(self.elastic().layers().len());
        for (index, (layer, geometry)) in self
            .elastic()
            .layers()
            .iter()
            .zip(self.elastic().layer_integrals())
            .enumerate()
        {
            let eta = layer.material.property(property_name).ok_or_else(|| {
                invalid("each transverse shear layer needs an explicit shear viscosity")
            })?;
            if eta.requirement().quantity() != QuantitySpec::dimensional(DynViscosity::DIMS)
                || !eta.value_si().is_finite()
                || eta.value_si() < 0.0
            {
                return Err(invalid(
                    "transverse shear viscosity must be nonnegative dimension-only Pa s",
                ));
            }
            let interval = eta.angular_frequency_band_rad_s()?.ok_or_else(|| {
                invalid("transverse shear viscosity needs an explicit frequency validity band")
            })?;
            band = (band.0.max(interval.0), band.1.min(interval.1));
            for coefficient in self.source_properties(index)?.into_iter().chain([eta]) {
                if !coefficient.is_constant_at_fixed_environment() {
                    return Err(invalid(
                        "transverse shear loss requires fixed-state constant density, modulus and viscosity",
                    ));
                }
                if let Some((lo, hi)) = coefficient.angular_frequency_band_rad_s()? {
                    band = (band.0.max(lo), band.1.min(hi));
                }
            }
            let value = (eta.value_si() * geometry.area_m2()) * self.shear_correction();
            if !value.is_finite() || (eta.value_si() > 0.0 && value == 0.0) {
                return Err(invalid(
                    "layer transverse shear viscous rigidity is unrepresentable",
                ));
            }
            total += value;
            if !total.is_finite() {
                return Err(invalid("transverse shear viscous rigidity sum overflow"));
            }
            layers.push(value);
        }
        if !band.0.is_finite() || !band.1.is_finite() || band.0 < 0.0 || band.1 < band.0 {
            return Err(invalid(
                "transverse shear coefficient frequency bands do not intersect",
            ));
        }
        Ok(ConcentricTransverseViscosity {
            shear: self,
            viscosity_property: property_name.into(),
            layer_viscous_rigidity_n_s: layers,
            viscous_rigidity_n_s: total,
            omega_band_rad_s: band,
        })
    }
}

impl ConcentricTransverseViscosity {
    /// Original geometry, correction, elastic sources and retained state.
    #[must_use]
    pub const fn shear(&self) -> &ConcentricShearSection {
        &self.shear
    }
    /// Explicit shear-viscosity selector.
    #[must_use]
    pub fn viscosity_property(&self) -> &str {
        &self.viscosity_property
    }
    /// Each kappa*eta_i*A_i [N s], in original spatial layer order.
    #[must_use]
    pub fn layer_viscous_rigidity_n_s(&self) -> &[f64] {
        &self.layer_viscous_rigidity_n_s
    }
    /// C_s [N s] multiplying either transverse shear-strain rate.
    #[must_use]
    pub const fn viscous_rigidity_n_s(&self) -> f64 {
        self.viscous_rigidity_n_s
    }
    /// Intersection of the original fixed-state source intervals [rad/s].
    #[must_use]
    pub const fn omega_band_rad_s(&self) -> (f64, f64) {
        self.omega_band_rad_s
    }
}
