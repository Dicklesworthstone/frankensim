//! Streamwise Boussinesq force for the fully developed channel screen.
//!
//! Subtracting the reference hydrostatic pressure from rho*g with
//! rho/rho_ref = 1 - beta*(T-T_ref) gives a = -g*beta*(T-T_ref).
//! See the OpenFOAM Boussinesq equation of state:
//! https://doc.openfoam.com/2306/tools/processing/models/thermophysical/equation-of-state/rtm/Boussinesq/
//!
//! Gravity is parallel to U; temperature varies only across the channel.
//! The gradient-diffusion buoyancy production of k is proportional to
//! g dot grad(T), which is zero in this restricted geometry. The existing
//! shear production remains active. This is not a stratified or 3-D model.

use super::RansTransportError;

/// Admitted thermal expansion force along the mean-flow coordinate.
///
/// Pressure supplied to the stepping API excludes reference-density
/// hydrostatic pressure. Positive gravity points along positive U, so a
/// hotter parcel accelerates opposite that direction when beta is positive.
/// Density, heat capacity, viscosity and conductivity otherwise stay constant.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RansBoussinesqForce {
    beta_per_k: f64,
    reference_temperature_k: f64,
    gravity_streamwise_m_s2: f64,
    maximum_relative_density_change: f64,
}

impl RansBoussinesqForce {
    /// Declare material expansion, gravity and a small-density-departure limit.
    ///
    /// beta [1/K] is finite and nonnegative, reference temperature [K] is
    /// positive, and signed gravity [m/s²] is finite. The caller's maximum
    /// |beta*(T-T_ref)| must be in (0, 0.1]. This admission limit is a scope
    /// restriction, not a bound on Boussinesq approximation error. Zero beta
    /// or gravity is permitted and supplies exactly zero acceleration.
    pub fn try_new(
        beta_per_k: f64,
        reference_temperature_k: f64,
        gravity_streamwise_m_s2: f64,
        maximum_relative_density_change: f64,
    ) -> Result<Self, RansTransportError> {
        if !beta_per_k.is_finite()
            || beta_per_k < 0.0
            || !reference_temperature_k.is_finite()
            || reference_temperature_k <= 0.0
            || !gravity_streamwise_m_s2.is_finite()
            || !maximum_relative_density_change.is_finite()
            || maximum_relative_density_change <= 0.0
            || maximum_relative_density_change > 0.1
        {
            return Err(RansTransportError::InvalidInput(
                "Boussinesq material, gravity or density-departure limit",
            ));
        }
        Ok(Self {
            beta_per_k,
            reference_temperature_k,
            gravity_streamwise_m_s2,
            maximum_relative_density_change,
        })
    }

    pub(super) fn acceleration(
        self,
        temperature_k: f64,
        node: usize,
    ) -> Result<f64, RansTransportError> {
        let departure = self.beta_per_k * (temperature_k - self.reference_temperature_k);
        if !temperature_k.is_finite()
            || temperature_k <= 0.0
            || !departure.is_finite()
            || departure.abs() > self.maximum_relative_density_change
        {
            return Err(RansTransportError::Arithmetic {
                equation: "Boussinesq density-departure limit",
                node,
            });
        }
        let acceleration = -self.gravity_streamwise_m_s2 * departure;
        if !acceleration.is_finite() {
            return Err(RansTransportError::Arithmetic {
                equation: "Boussinesq acceleration",
                node,
            });
        }
        Ok(if acceleration == 0.0 {
            0.0
        } else {
            acceleration
        })
    }
}
