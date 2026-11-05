//! Source-owned surface coefficients on one original bulk enthalpy chart.

use crate::enthalpy::{DepositionError, ThermalLaw, finite, nonnegative};

// The same rounded CODATA value as the native gray radiation owners.
const SIGMA: f64 = 5.670_374_419e-8;

/// Effective gray surface response at one original bulk enthalpy.
///
/// Derivatives hold geometry, bath temperatures and every external source
/// coordinate fixed. They include the source's original temperature/phase
/// chain rules; latent phase response can remain nonzero when dT/dh is zero.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitGraySurfaceState {
    /// Total hemispherical gray emissivity, dimensionless in [0,1].
    pub emissivity: f64,
    /// Complete d(emissivity)/dh, kg/J.
    pub emissivity_enthalpy_derivative_kg_j: f64,
    /// Exterior convective film coefficient, W/(m² K), nonnegative.
    pub convection_w_m2_k: f64,
    /// Complete d(htc)/dh, kg W/(m² K J).
    pub convection_enthalpy_derivative_kg_w_m2_k_j: f64,
}

impl ImplicitGraySurfaceState {
    /// Admit passive coefficients and finite tangents without clipping them.
    /// A source-domain endpoint can have a nonzero inward derivative at zero
    /// emissivity/htc; its physical onset must not disappear from the Jacobian.
    pub fn admit(self) -> Result<(), DepositionError> {
        if !self.emissivity.is_finite()
            || !(0.0..=1.0).contains(&self.emissivity)
            || !self.emissivity_enthalpy_derivative_kg_j.is_finite()
            || !nonnegative(self.convection_w_m2_k)
            || !self.convection_enthalpy_derivative_kg_w_m2_k_j.is_finite()
        {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(())
    }

    pub(super) fn fixed(emissivity: f64, convection_w_m2_k: f64) -> Self {
        Self {
            emissivity,
            emissivity_enthalpy_derivative_kg_j: 0.0,
            convection_w_m2_k,
            convection_enthalpy_derivative_kg_w_m2_k_j: 0.0,
        }
    }

    pub(super) fn has_enthalpy_response(self) -> bool {
        self.emissivity_enthalpy_derivative_kg_j != 0.0
            || self.convection_enthalpy_derivative_kg_w_m2_k_j != 0.0
    }

    /// Partial derivative of the exterior flux at a FIXED surface temperature.
    /// Deltas come from the scalar solver's relative coordinate, preserving
    /// small fluxes without subtracting rounded Kelvin or fourth powers.
    pub(super) fn flux_enthalpy_derivative(
        self,
        surface_k: f64,
        radiative_bath_k: f64,
        gas_delta_k: f64,
        radiation_delta_k: f64,
    ) -> Result<f64, DepositionError> {
        let convection = finite(self.convection_enthalpy_derivative_kg_w_m2_k_j * gas_delta_k)?;
        let radiation = if self.emissivity_enthalpy_derivative_kg_j == 0.0 {
            0.0
        } else {
            // Do not divide the radiative power by epsilon: epsilon may be
            // exactly zero at the admitted onset of a source/phase branch.
            let coefficient = finite(
                (SIGMA * self.emissivity_enthalpy_derivative_kg_j)
                    * (surface_k + radiative_bath_k)
                    * surface_k.mul_add(surface_k, radiative_bath_k * radiative_bath_k),
            )?;
            finite(coefficient * radiation_delta_k)?
        };
        finite(convection + radiation)
    }
}

/// Immutable effective surface law on the SAME caloric object as the cell.
///
/// This law explicitly maps the original BULK specific enthalpy to emissivity
/// and convection. A material owner supplies its own supported source domain
/// and temperature/phase derivatives. It is a cell-state surface closure:
/// solved exterior temperature is not substituted into the source coordinates,
/// and no distinct surface phase, evolving coating, flow correlation, or bath
/// state is inferred. Sources requiring those additional unknowns need their
/// own closure. Geometry, pressure, roughness and other conditions remain fixed
/// by the immutable source law.
///
/// Coefficients are resolved once per original enthalpy trial, then held fixed
/// during that trial's bounded surface solve. They must include neither bulk
/// half-cell/coating resistance nor an additional heat source. The existing
/// boundary owner retains heat, bath entropy, storage and time integration.
pub trait ImplicitGraySurfaceLaw {
    fn caloric_law(&self) -> &dyn ThermalLaw;
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2];
    fn state_at_specific_enthalpy(
        &self,
        specific_h_j_kg: f64,
    ) -> Result<ImplicitGraySurfaceState, DepositionError>;
}

pub(super) struct BoundGraySurface<'a> {
    law: &'a dyn ImplicitGraySurfaceLaw,
    bounds_j_kg: [f64; 2],
}

impl<'a> BoundGraySurface<'a> {
    pub(super) fn new(
        caloric: &dyn ThermalLaw,
        law: &'a dyn ImplicitGraySurfaceLaw,
    ) -> Result<Self, DepositionError> {
        let [lo, hi] = law.enthalpy_bounds_j_kg();
        if !lo.is_finite()
            || !hi.is_finite()
            || lo >= hi
            || !std::ptr::addr_eq(
                caloric as *const dyn ThermalLaw,
                law.caloric_law() as *const dyn ThermalLaw,
            )
        {
            return Err(DepositionError::InvalidInput(
                "gray surface source requires finite bounds and the original caloric object",
            ));
        }
        Ok(Self {
            law,
            bounds_j_kg: [lo, hi],
        })
    }

    pub(super) fn evaluate(
        &self,
        h: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitGraySurfaceState, DepositionError> {
        let [lo, hi] = self.bounds_j_kg;
        if !h.is_finite() || !(lo..=hi).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        checkpoint()?;
        let state = self.law.state_at_specific_enthalpy(h)?;
        state.admit()?;
        checkpoint()?;
        Ok(state)
    }
}
