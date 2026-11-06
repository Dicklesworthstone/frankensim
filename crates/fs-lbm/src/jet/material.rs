//! Original liquid-material binding for cylindrical capillary-jet screening.
//!
//! Density, viscosity and surface tension come from the same admitted source
//! point as the existing free-surface and film owners. The inviscid dispersion
//! omits viscous stresses; its explicit Ohnesorge budget is an applicability
//! screen, not a bound on growth-rate error or a measured breakup certificate.
//! The viscous binding retains the source viscosity in the Newtonian slender
//! model, with an explicit wavelength budget instead of an inviscid allowance.

use super::viscous::ViscousSlenderJet;
use super::{CapillaryJet, JetBandMaximum, JetDispersion, JetError};
use fs_material::liquid::surface::ResolvedLiquidSurfaceState;

/// An inviscid jet together with its original immutable material and viscosity
/// applicability. No source record, receipt or query coordinate is cloned.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialCapillaryJet<'a> {
    material: &'a ResolvedLiquidSurfaceState,
    jet: CapillaryJet,
    ohnesorge: f64,
    maximum_ohnesorge: f64,
}

impl<'a> MaterialCapillaryJet<'a> {
    /// Admit the original liquid for a low-viscosity cylindrical jet.
    ///
    /// The caller supplies the actual radius and a finite nonnegative maximum
    /// Oh = mu/sqrt(rho*sigma*R). There is no default material or viscosity
    /// cutoff. This gate does not establish negligible surrounding-gas inertia,
    /// a uniform unstretched jet, or linear disturbance amplitudes; those remain
    /// assumptions of the native model. A changed temperature/composition needs
    /// a newly resolved source point and a new binding.
    ///
    /// # Errors
    /// Refuses omitted surface tension, an evolving rheology source, invalid
    /// geometry/budget, unrepresentable Oh, or a source exceeding that budget.
    pub fn inviscid(
        material: &'a ResolvedLiquidSurfaceState,
        radius_m: f64,
        maximum_ohnesorge: f64,
    ) -> Result<Self, JetError> {
        if material.rheology_axes().is_some() || material.rheology_pressure_axis().is_some() {
            return Err(JetError::InvalidInput(
                "a constant Newtonian jet cannot freeze an evolving rheology source",
            ));
        }
        if !maximum_ohnesorge.is_finite() || maximum_ohnesorge < 0.0 {
            return Err(JetError::InvalidInput(
                "maximum Ohnesorge number must be finite and nonnegative",
            ));
        }
        let sigma = material.surface_tension_n_m().ok_or(JetError::InvalidInput(
            "a material capillary jet requires explicitly sourced surface tension",
        ))?;
        let jet = CapillaryJet::inviscid(radius_m, material.density_kg_m3(), sigma)?;
        // Logarithmic construction avoids overflowing rho*sigma*R or squaring
        // a small radius. It uses the same deterministic scalar math owner as
        // other material/film evaluations; this is a one-time admission cost.
        let ohnesorge = fs_math::det::exp(
            fs_math::det::ln(material.kinematic_viscosity_m2_s())
                + 0.5
                    * (fs_math::det::ln(material.density_kg_m3())
                        - fs_math::det::ln(sigma)
                        - fs_math::det::ln(radius_m)),
        );
        if !ohnesorge.is_finite() || ohnesorge <= 0.0 {
            return Err(JetError::NonRepresentable);
        }
        if ohnesorge > maximum_ohnesorge {
            return Err(JetError::InvalidInput(
                "source viscosity exceeds the declared inviscid-jet Ohnesorge budget",
            ));
        }
        Ok(Self {
            material,
            jet,
            ohnesorge,
            maximum_ohnesorge,
        })
    }

    /// The exact source state, including its property selections and receipts.
    #[must_use]
    pub const fn material(self) -> &'a ResolvedLiquidSurfaceState {
        self.material
    }

    /// Native allocation-free dispersion/band owner, using these source values.
    #[must_use]
    pub const fn jet(self) -> CapillaryJet {
        self.jet
    }

    /// Actual viscosity-to-capillary-inertia ratio at this radius.
    #[must_use]
    pub const fn ohnesorge(self) -> f64 {
        self.ohnesorge
    }

    /// Explicit scenario allowance retained beside the actual source ratio.
    #[must_use]
    pub const fn maximum_ohnesorge(self) -> f64 {
        self.maximum_ohnesorge
    }
}

impl JetDispersion for MaterialCapillaryJet<'_> {
    fn radius_m(self) -> f64 {
        self.jet.radius_m()
    }

    fn fastest_mode(self, wavenumbers_rad_m: [f64; 2]) -> Result<JetBandMaximum, JetError> {
        self.jet.fastest_mode(wavenumbers_rad_m)
    }
}

/// Newtonian viscous jet retaining the original source state and receipts.
///
/// No viscosity cutoff is inferred. The native long-wave approximation uses
/// the source density, kinematic viscosity and surface tension together.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialViscousSlenderJet<'a> {
    material: &'a ResolvedLiquidSurfaceState,
    jet: ViscousSlenderJet,
}

impl<'a> MaterialViscousSlenderJet<'a> {
    /// Bind an admitted source at the actual radius with an explicit kR budget.
    ///
    /// A changed temperature/composition requires a new source resolution and
    /// binding. The budget must be in (0,1]; it declares the admitted wavelengths
    /// and does not certify the accuracy of the slender approximation.
    ///
    /// # Errors
    /// Refuses omitted surface tension, evolving non-Newtonian rheology,
    /// invalid geometry/budget, or nonrepresentable native coefficients.
    pub fn try_new(
        material: &'a ResolvedLiquidSurfaceState,
        radius_m: f64,
        maximum_kr: f64,
    ) -> Result<Self, JetError> {
        if material.rheology_axes().is_some() || material.rheology_pressure_axis().is_some() {
            return Err(JetError::InvalidInput(
                "a constant Newtonian jet cannot freeze an evolving rheology source",
            ));
        }
        let sigma = material.surface_tension_n_m().ok_or(JetError::InvalidInput(
            "a material capillary jet requires explicitly sourced surface tension",
        ))?;
        let jet = ViscousSlenderJet::new(
            radius_m,
            material.density_kg_m3(),
            sigma,
            material.kinematic_viscosity_m2_s(),
            maximum_kr,
        )?;
        Ok(Self { material, jet })
    }

    /// Exact original material, with its selections, coordinates and receipts.
    #[must_use]
    pub const fn material(self) -> &'a ResolvedLiquidSurfaceState {
        self.material
    }

    /// Native allocation-free viscous dispersion using this source's values.
    #[must_use]
    pub const fn jet(self) -> ViscousSlenderJet {
        self.jet
    }
}

impl JetDispersion for MaterialViscousSlenderJet<'_> {
    fn radius_m(self) -> f64 {
        self.jet.radius_m()
    }

    fn fastest_mode(self, wavenumbers_rad_m: [f64; 2]) -> Result<JetBandMaximum, JetError> {
        self.jet.fastest_mode(wavenumbers_rad_m)
    }
}
