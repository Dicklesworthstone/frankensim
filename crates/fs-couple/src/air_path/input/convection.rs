//! Source-gas forced convection into an explicitly prescribed reservoir.
//!
//! The original isobaric gas law supplies the film-state rho, mu, Cp and k;
//! fs-convection owns every correlation and its regime/evidence. This adapter
//! supplies signed boundary heat, not a new thermal or fluid time integrator.
//! Reservoir temperature, pressure, composition and velocity are prescribed.
//! No finite gas inventory, buoyancy, rarefaction, radiation or moving exposure
//! is inferred, and no transport coefficient is substituted from another gas.

use fs_convection::{
    CorrelationError, CorrelationId, CorrelationInputs, HeatTransferCoefficient,
    NusseltEvaluation, ThermalConductivity,
};
use fs_evidence::Evidence;
use fs_material::{
    MaterialError,
    fluid_transport::FluidTransport,
    gas_species::input::{GasRecord, temperature::IsobaricGasLaw},
};
use fs_qty::{Length, Temperature};

/// Explicit exterior geometry and flow regime, not an instrument/material name.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ForcedGasGeometry {
    /// Diameter is the characteristic length; uniform crossflow, not vibration.
    CylinderCrossflow,
    /// Streamwise length; the shared laminar average isothermal-plate relation.
    LaminarFlatPlate,
}

/// Fixed exposed geometry. Contact-covered area must be excluded by its owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ForcedGasBoundary {
    pub geometry: ForcedGasGeometry,
    pub characteristic_length_m: f64,
    pub exposed_area_m2: f64,
    /// Prescribed free-stream speed relative to the stationary surface [m/s].
    pub speed_m_s: f64,
}

/// Physical/source/regime refusal. No failed law silently selects another one.
#[derive(Debug)]
pub enum GasConvectionError {
    Invalid(&'static str),
    Material(MaterialError),
    Correlation(CorrelationError),
}
impl core::fmt::Display for GasConvectionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "gas convection: {what}"),
            Self::Material(error) => write!(f, "gas convection material: {error}"),
            Self::Correlation(error) => write!(f, "gas convection correlation: {error}"),
        }
    }
}
impl std::error::Error for GasConvectionError {}
impl From<MaterialError> for GasConvectionError {
    fn from(error: MaterialError) -> Self { Self::Material(error) }
}
impl From<CorrelationError> for GasConvectionError {
    fn from(error: CorrelationError) -> Self { Self::Correlation(error) }
}

/// The same evaluated gas and model evidence that produced the boundary heat.
#[derive(Debug)]
pub struct GasConvectionSample {
    pub film_gas: GasRecord,
    pub nusselt: NusseltEvaluation,
    pub coefficient: Evidence<HeatTransferCoefficient>,
    pub wall_temperature_k: f64,
    pub reservoir_temperature_k: f64,
    /// h A [W/K]; recomputed at every proposed wall temperature.
    pub conductance_w_k: f64,
    /// Positive INTO the wall; the external reservoir receives its negative.
    pub heat_into_wall_w: f64,
    /// Q (1/T_wall - 1/T_reservoir) [W/K], nonnegative for this passive boundary.
    /// This is an instantaneous production, not time-integrated entropy.
    pub entropy_generation_w_k: f64,
}

impl ForcedGasBoundary {
    pub fn validate(self) -> Result<(), GasConvectionError> {
        if [self.characteristic_length_m, self.exposed_area_m2, self.speed_m_s]
            .into_iter().any(|x| !x.is_finite() || x <= 0.0)
        {
            return Err(GasConvectionError::Invalid(
                "forced flow requires positive finite length, exposed area and speed",
            ));
        }
        Ok(())
    }

    /// Evaluate an original source law at (T_wall + T_bulk)/2, holding its
    /// reference pressure/composition. Heat uses the FULL wall/bulk difference.
    ///
    /// A thermochemical law retains the original source model and pinned
    /// claims; a snapshot alone is not sufficient. Humid laws preserve water
    /// fraction, not relative humidity. Their phase/domain refusals propagate.
    /// Correlation model discrepancy is retained, not a propagated source-error
    /// bound. A caller must insert this heat once into its existing thermal
    /// residual, reevaluate each trial, and account the external reservoir.
    pub fn evaluate(
        self, gas: &IsobaricGasLaw, wall_temperature_k: f64,
    ) -> Result<GasConvectionSample, GasConvectionError> {
        self.validate()?;
        if !wall_temperature_k.is_finite() || wall_temperature_k <= 0.0 {
            return Err(GasConvectionError::Invalid("wall temperature must be positive finite kelvin"));
        }
        let reservoir_temperature_k = gas.reference().state().temperature;
        let film_temperature_k = 0.5 * wall_temperature_k + 0.5 * reservoir_temperature_k;
        let film_gas = gas.at_temperature(Temperature::new(film_temperature_k))?;
        let transport = film_gas.transport_state()?;
        let reynolds = transport.reynolds(self.speed_m_s, self.characteristic_length_m)?;
        let id = match self.geometry {
            ForcedGasGeometry::CylinderCrossflow => CorrelationId::ChurchillBernsteinCylinder,
            ForcedGasGeometry::LaminarFlatPlate => CorrelationId::FlatPlateLaminarAverage,
        };
        let nusselt = fs_convection::evaluate(id, CorrelationInputs::forced(reynolds, transport.prandtl()))?;
        let coefficient = nusselt.heat_transfer_coefficient(
            ThermalConductivity::new(transport.thermal_conductivity_w_m_k()),
            Length::new(self.characteristic_length_m),
        )?;
        let conductance_w_k = coefficient.value.value() * self.exposed_area_m2;
        let difference = reservoir_temperature_k - wall_temperature_k;
        let heat_into_wall_w = conductance_w_k * difference;
        // Use the original temperature difference rather than subtracting
        // nearly equal reciprocal temperatures; avoid a potentially huge T*T.
        let entropy_generation_w_k = (heat_into_wall_w / wall_temperature_k)
            * (difference / reservoir_temperature_k);
        if !conductance_w_k.is_finite() || conductance_w_k <= 0.0
            || !heat_into_wall_w.is_finite()
            || !entropy_generation_w_k.is_finite() || entropy_generation_w_k < 0.0
            || (difference != 0.0 && (heat_into_wall_w == 0.0 || entropy_generation_w_k == 0.0))
        {
            return Err(GasConvectionError::Invalid("boundary heat or entropy is unrepresentable"));
        }
        Ok(GasConvectionSample { film_gas, nusselt, coefficient,
            wall_temperature_k, reservoir_temperature_k, conductance_w_k,
            heat_into_wall_w, entropy_generation_w_k })
    }
}

#[cfg(test)]
mod tests;
