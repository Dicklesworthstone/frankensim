//! Source response at the actual exterior temperature and original stored energy.

use super::ImplicitGraySurfaceState;
use crate::enthalpy::{DepositionError, ThermalLaw};

/// Passive gray coefficients with independent bulk-energy and surface tangents.
///
/// `coefficients` retains the original bulk h/Cp or u/Cv coordinate. Its energy
/// derivatives are PARTIAL derivatives at fixed exterior and gas temperatures.
/// The two derivatives below hold that original energy and gas temperature
/// fixed. The boundary owner combines both paths after solving the exterior
/// temperature; a source must not include that implicit response itself.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitGrayExteriorState {
    pub coefficients: ImplicitGraySurfaceState,
    /// Partial d(emissivity)/dTs, 1/K.
    pub emissivity_temperature_derivative_per_k: f64,
    /// Partial d(htc)/dTs, W/(m² K²).
    pub convection_temperature_derivative_w_m2_k2: f64,
}

impl ImplicitGrayExteriorState {
    /// Admit passive coefficients and finite, possibly signed partials.
    /// Zero coefficients can retain an inward source derivative at onset.
    pub fn admit(self) -> Result<(), DepositionError> {
        self.coefficients.admit()?;
        if !self.emissivity_temperature_derivative_per_k.is_finite()
            || !self.convection_temperature_derivative_w_m2_k2.is_finite()
        {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(())
    }
}

/// Immutable surface law on the original finite body's caloric object.
///
/// The original specific energy, actual solved exterior temperature and actual
/// gas bath temperature are DISTINCT physical inputs. The energy coordinate is
/// h for the original isobaric law or u for its isochoric law; substituting a
/// surface-derived bulk energy would change the physical model. A coating can
/// separate this exterior face from the supported material face.
///
/// Bounds describe finite nonempty inclusive source intervals. Surface bounds
/// need not contain the bulk or either bath temperature. A law also checks its
/// actual gas-temperature support and any nonrectangular source restrictions
/// on every query. No extrapolation, clipping or material-name inference is
/// authorized by these bounds. Geometry, flow, pressure and other external
/// conditions belong to the immutable source; none evolves implicitly here.
///
/// The same bounded surface solve evaluates this law at each trial Ts. Return
/// coefficients and BOTH partials at that exact input without changing history.
/// HTC includes only the exterior convective film, and emissivity is total
/// hemispherical gray emissivity. Neither includes material/coating resistance,
/// a separate heat source, surface storage or a finite surrounding's energy.
/// Positive coefficients establish passive heat directions; they do not alone
/// establish a monotone flux law or uniqueness of the surface balance.
pub trait ImplicitGrayExteriorLaw {
    fn caloric_law(&self) -> &dyn ThermalLaw;
    fn energy_bounds_j_kg(&self) -> [f64; 2];
    fn surface_temperature_bounds_k(&self) -> [f64; 2];
    fn state_at_surface(
        &self,
        specific_energy_j_kg: f64,
        surface_temperature_k: f64,
        gas_temperature_k: f64,
    ) -> Result<ImplicitGrayExteriorState, DepositionError>;
}
