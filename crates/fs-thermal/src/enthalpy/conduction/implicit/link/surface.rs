//! Passive two-surface exchange on two original finite caloric coordinates.

use crate::enthalpy::{DepositionError, ThermalLaw};

/// One source flux density and independent surface/energy partial derivatives.
///
/// Positive flux is from original A to B and equals G*(Ts_A-Ts_B). G includes
/// the entire zero-storage interface, with no material half-cell resistance.
/// Surface partials hold BOTH original energies fixed. Energy partials hold
/// BOTH surface temperatures fixed and retain each original h/Cp or u/Cv
/// coordinate. The consuming layer owner supplies the implicit surface response.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitSurfaceExchangeState {
    pub flux_a_to_b_w_m2: f64,
    pub conductance_w_m2_k: f64,
    /// Signed [partial j/partial Ts_A, partial j/partial Ts_B].
    pub surface_derivatives_w_m2_k: [f64; 2],
    /// Signed [partial j/partial e_A, partial j/partial e_B], W kg/(m² J).
    pub energy_derivatives_w_kg_m2_j: [f64; 2],
}

impl ImplicitSurfaceExchangeState {
    /// Construct a passive source from its original conductance and partials.
    ///
    /// The supplied derivatives differentiate G, not j. They may include a
    /// declared conductive and radiative surface response at the same query.
    /// This constructor applies the exact product rule, including G=0 onset.
    /// It infers no geometry, emissivity, view factor or source authority.
    pub fn from_conductance(
        surface_temperatures_k: [f64; 2],
        conductance_w_m2_k: f64,
        conductance_surface_derivatives_w_m2_k2: [f64; 2],
        conductance_energy_derivatives_w_kg_m2_k_j: [f64; 2],
    ) -> Result<Self, DepositionError> {
        if surface_temperatures_k
            .iter()
            .any(|t| !t.is_finite() || *t <= 0.0)
            || !conductance_w_m2_k.is_finite()
            || conductance_w_m2_k < 0.0
            || conductance_surface_derivatives_w_m2_k2
                .iter()
                .chain(&conductance_energy_derivatives_w_kg_m2_k_j)
                .any(|x| !x.is_finite())
        {
            return Err(DepositionError::MaterialDomain);
        }
        let difference = surface_temperatures_k[0] - surface_temperatures_k[1];
        let state = Self {
            flux_a_to_b_w_m2: conductance_w_m2_k * difference,
            conductance_w_m2_k,
            surface_derivatives_w_m2_k: [
                conductance_surface_derivatives_w_m2_k2[0].mul_add(difference, conductance_w_m2_k),
                conductance_surface_derivatives_w_m2_k2[1].mul_add(difference, -conductance_w_m2_k),
            ],
            energy_derivatives_w_kg_m2_j: conductance_energy_derivatives_w_kg_m2_k_j
                .map(|derivative| derivative * difference),
        };
        state.admit()?;
        if difference != 0.0 && conductance_w_m2_k > 0.0 && state.flux_a_to_b_w_m2 == 0.0 {
            return Err(DepositionError::Arithmetic);
        }
        Ok(state)
    }

    /// Admit finite passive coefficients and finite, possibly signed partials.
    /// The consuming layer owner also checks j=G*(Ts_A-Ts_B) at the actual query.
    pub fn admit(self) -> Result<(), DepositionError> {
        if !self.flux_a_to_b_w_m2.is_finite()
            || !self.conductance_w_m2_k.is_finite()
            || self.conductance_w_m2_k < 0.0
            || self
                .surface_derivatives_w_m2_k
                .iter()
                .chain(&self.energy_derivatives_w_kg_m2_j)
                .any(|x| !x.is_finite())
        {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(())
    }
}

/// Immutable interface source on two ordered ORIGINAL caloric objects.
///
/// Each energy is the original body's specific h or u; solved surface
/// temperatures are independent inputs and never substitute for bulk energies.
/// This can express a temperature-dependent contact, gas gap or a declared
/// radiative/conductive exchange between equal-area faces. The source owns
/// physical assumptions, external conditions and coefficient data; a passive
/// scalar flux law does not infer a general enclosure or its view factors.
///
/// Bounds are outer inclusive source intervals. They do not establish joint
/// rectangular support: each query must admit its complete coupled conditions.
/// Surface support need not contain either bulk temperature or overlap the
/// opposite surface's chart. No extrapolation or changing source history occurs.
///
/// The original conductivity-integral owner solves both surface temperatures
/// and chains BOTH sets of partials through the same scalar balance. The finite
/// thermal owner retains mass, stored energy, entropy, accepted time and
/// cancellation. Extra surface storage or internal heat generation requires
/// its explicitly represented physical owner.
pub trait ImplicitSurfaceExchangeLaw {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2];
    fn energy_bounds_j_kg(&self) -> [[f64; 2]; 2];
    fn surface_temperature_bounds_k(&self) -> [[f64; 2]; 2];
    fn state_at_surfaces(
        &self,
        energies_j_kg: [f64; 2],
        surface_temperatures_k: [f64; 2],
    ) -> Result<ImplicitSurfaceExchangeState, DepositionError>;
}
