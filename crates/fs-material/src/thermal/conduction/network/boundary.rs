//! Original material-to-bath flux and tangent in native implicit heat networks.

use super::super::ReservoirLayerGeometry;
use super::implicit::temperature_slope;
use super::{LayerConductionConfig, LayerSourceError, ResolvedConductionLayer};
use fs_thermal::enthalpy::conduction::implicit::boundary::{
    ImplicitBoundaryGeometry, ImplicitBoundaryLaw, ImplicitBoundaryState,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

/// Source-qualified planar half-cell followed by its original boundary film.
///
/// Attach this borrowed law with `ImplicitConductionChain::with_boundary_laws`
/// or `ImplicitConductionGraph::with_boundary_laws`. The native owner supplies
/// each current reservoir or ambient temperature, retains the original area,
/// and owns all enthalpy storage, heat, entropy and accepted time evolution.
/// The material source, pure-phase caloric branch, distance and areal film
/// resistance stay fixed. A bath may lie outside the source's temperature span
/// only when a finite film keeps the actual material face inside that span.
///
/// Geometry and total film resistance are prepared once. Trial evaluation uses
/// the existing Kirchhoff integral and bounded scalar face solve, with no
/// source query, allocation or independent boundary heat-capacity account.
pub struct ImplicitLayerBoundary<'a, 'm> {
    layer: &'a ResolvedConductionLayer<'m>,
    geometry: ImplicitBoundaryGeometry,
    flux_geometry: ReservoirLayerGeometry,
    film_resistance_k_w: f64,
    config: LayerConductionConfig,
}

impl<'a, 'm> ImplicitLayerBoundary<'a, 'm> {
    /// Prepare the same centre-to-face distance, face area and areal film
    /// resistance that the original native conduction observer binds.
    /// Zero film selects a prescribed face temperature; no extra film is added
    /// from the cell's authored conductance or convection coefficient.
    pub fn new(
        layer: &'a ResolvedConductionLayer<'m>,
        geometry: ImplicitBoundaryGeometry,
        config: LayerConductionConfig,
    ) -> Result<Self, DepositionError> {
        let flux_geometry =
            ReservoirLayerGeometry::planar(geometry.area_m2, geometry.half_distance_m)?;
        let film = geometry.surface_resistance_m2_k_w;
        if !film.is_finite()
            || film < 0.0
            || config.max_iterations == 0
            || !config.absolute_heat_tolerance_w.is_finite()
            || config.absolute_heat_tolerance_w < 0.0
            || !config.relative_heat_tolerance.is_finite()
            || config.relative_heat_tolerance < 0.0
        {
            return Err(DepositionError::InvalidInput(
                "implicit material boundary needs a finite nonnegative film and bounded solve controls",
            ));
        }
        let film_resistance_k_w = film / geometry.area_m2;
        if !film_resistance_k_w.is_finite() || (film > 0.0 && film_resistance_k_w == 0.0) {
            return Err(DepositionError::Arithmetic);
        }
        Ok(Self {
            layer,
            geometry,
            flux_geometry,
            film_resistance_k_w,
            config,
        })
    }
}

impl ImplicitBoundaryLaw for ImplicitLayerBoundary<'_, '_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.layer.caloric_curve()
    }

    fn state_at_specific_enthalpy(
        &self,
        specific_h_j_kg: f64,
        bath_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitBoundaryState, DepositionError> {
        checkpoint()?;
        if geometry.area_m2 != self.geometry.area_m2
            || geometry.half_distance_m != self.geometry.half_distance_m
            || geometry.surface_resistance_m2_k_w != self.geometry.surface_resistance_m2_k_w
        {
            return Err(DepositionError::InvalidInput(
                "implicit material boundary requires its original area, distance and film",
            ));
        }
        let temperature = self
            .layer
            .temperature_at_enthalpy(specific_h_j_kg)
            .map_err(|error| match error {
                LayerSourceError::Thermal(error) => error,
                _ => DepositionError::MaterialDomain,
            })?;
        let slope = temperature_slope(self.layer, specific_h_j_kg, temperature)?;
        let response = self.flux_geometry.evaluate_material_to_bath(
            self.layer.conductivity_curve(),
            [temperature, bath_temperature_k],
            self.film_resistance_k_w,
            self.config,
            &mut *checkpoint,
        )?;
        // Native heat accounting applies the secant, so admit its actual power
        // against the original scalar solve's tolerance before publication.
        let applied = response.conductance_w_k * (temperature - bath_temperature_k);
        let allowance = self.config.absolute_heat_tolerance_w
            + self.config.relative_heat_tolerance
                * applied.abs().max(response.power_to_reservoir_w.abs());
        let derivative = response.endpoint_derivative_w_k[0] * slope;
        if !applied.is_finite()
            || !allowance.is_finite()
            || (applied - response.power_to_reservoir_w).abs() > allowance
            || !derivative.is_finite()
        {
            return Err(DepositionError::Arithmetic);
        }
        checkpoint()?;
        Ok(ImplicitBoundaryState {
            conductance_w_k: response.conductance_w_k,
            power_enthalpy_derivative_w_kg_j: derivative,
        })
    }
}
