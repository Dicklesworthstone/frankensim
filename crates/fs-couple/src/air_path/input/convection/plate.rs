//! Source-gas heat on a physical strip of an external isothermal plate.

use super::*;
use fs_convection::developing::evaluate_flat_plate_section;

/// Explicit source regime; a strip never switches correlations at a mesh edge.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FlatPlateRegime {
    /// Original laminar isothermal-plate mean and its Reynolds/Prandtl domain.
    Laminar,
    /// Original turbulent mean includes its laminar leading-edge correction.
    Turbulent,
}

impl FlatPlateRegime {
    const fn correlation(self) -> CorrelationId {
        match self {
            Self::Laminar => CorrelationId::FlatPlateLaminarAverage,
            Self::Turbulent => CorrelationId::FlatPlateTurbulentAverage,
        }
    }
}

/// One exposed strip `[start_m, end_m]` measured downstream from the physical
/// leading edge. `end_m`, not strip width, is the correlation length.
///
/// Supply disjoint physical areas; contact-covered area is excluded by its
/// original owner. Reversing flow requires distances from the opposite edge.
/// The prescribed positive speed is relative to the stationary wall.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlatPlateGasBoundary {
    /// One declared correlation for both cumulative endpoint evaluations.
    pub regime: FlatPlateRegime,
    /// Distance to the upstream strip edge [m], including an exact zero entrance.
    pub start_m: f64,
    /// Distance to the downstream edge [m]; strictly greater than `start_m`.
    pub end_m: f64,
    /// Actual exposed strip area [m²] after contact masking by its owner.
    pub exposed_area_m2: f64,
    /// Prescribed positive stream speed relative to the stationary surface [m/s].
    pub speed_m_s: f64,
}

impl FlatPlateGasBoundary {
    pub fn validate(self) -> Result<(), GasConvectionError> {
        if !self.start_m.is_finite()
            || self.start_m < 0.0
            || !self.end_m.is_finite()
            || self.end_m <= self.start_m
            || !self.exposed_area_m2.is_finite()
            || self.exposed_area_m2 <= 0.0
            || !self.speed_m_s.is_finite()
            || self.speed_m_s <= 0.0
        {
            return Err(GasConvectionError::Invalid(
                "plate strip needs finite 0 <= start < end, positive area and positive relative speed",
            ));
        }
        if self.start_m > 0.0 && self.start_m / self.end_m == 0.0 {
            return Err(GasConvectionError::Invalid(
                "a positive plate-strip start cannot underflow to a new leading edge",
            ));
        }
        Ok(())
    }

    /// Reevaluate the original gas at mean film temperature and unchanged
    /// reservoir pressure/composition for every thermal trial. The shared
    /// section owner integrates the original mean correlation between both
    /// physical endpoints. Its source domains remain active at both ends.
    ///
    /// At start zero this gives exactly the original whole-plate coefficient.
    /// Subdividing a common isothermal plate of constant transverse width
    /// recovers its total heat while allowing downstream h to decrease. Using
    /// varying local wall temperatures is a held-property approximation, not a resolved boundary
    /// layer. Positive-start strips retain the shared owner's unbounded local
    /// model discrepancy rather than inheriting a whole-plate uncertainty band.
    /// Heat uses the full wall/bulk contrast and enters the existing thermal
    /// residual once; the film is a property query with no separate heat store.
    pub fn evaluate(
        self,
        gas: &IsobaricGasLaw,
        wall_temperature_k: f64,
    ) -> Result<GasConvectionSample, GasConvectionError> {
        self.validate()?;
        if !wall_temperature_k.is_finite() || wall_temperature_k <= 0.0 {
            return Err(GasConvectionError::Invalid(
                "wall temperature must be positive finite kelvin",
            ));
        }
        let reservoir_temperature_k = gas.reference().state().temperature;
        let property_gas = gas.at_temperature(Temperature::new(
            0.5 * wall_temperature_k + 0.5 * reservoir_temperature_k,
        ))?;
        let transport = property_gas.transport_state()?;
        let section = evaluate_flat_plate_section(
            self.regime.correlation(),
            transport.reynolds(self.speed_m_s, self.end_m)?,
            transport.prandtl(),
            self.start_m / self.end_m,
        )?;
        let coefficient = section.heat_transfer_coefficient(
            ThermalConductivity::new(transport.thermal_conductivity_w_m_k()),
            Length::new(self.end_m),
        )?;
        finish_coefficient(
            property_gas,
            None,
            None,
            Some(section),
            Some(coefficient),
            self.exposed_area_m2,
            wall_temperature_k,
            reservoir_temperature_k,
        )
    }
}
