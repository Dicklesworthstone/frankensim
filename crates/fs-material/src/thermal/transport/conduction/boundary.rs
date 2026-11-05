//! Original Cp/Cv storage with source-integrated material-to-bath heat flow.
//!
//! The existing reservoir-layer owner solves the material surface temperature
//! and film jump. The native thermal owner retains the heat, bath entropy and
//! time integration. This adapter supplies the exact derivative with respect
//! to the same original h or u; it adds no surface heat inventory.

use super::CaloricConductionLayer;
use crate::thermal::conduction::{
    LayerConductionConfig, ReservoirLayerGeometry, ReservoirLayerResponse,
};
use fs_thermal::enthalpy::conduction::implicit::boundary::{
    ImplicitBoundaryGeometry, ImplicitBoundaryLaw, ImplicitBoundaryState,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

/// A planar source-resolved half-cell and its one authored boundary film.
///
/// Attach this law to the original reservoir or ambient channel using the
/// native implicit chain/graph's `with_boundary_laws`. Its current bath
/// temperature, including scheduled changes, comes from that channel.
///
/// The material face must stay inside the original conductivity support. A
/// hotter or colder bath outside that support is permitted when a finite
/// film keeps the actual material face supported. The bath temperature is
/// never evaluated as though it were another material sample.
///
/// The source is single-phase and the normal profile is quasi-steady within
/// the half-cell. Geometry, composition and the original Cp/Cv constraint
/// remain fixed. Spatial/time accuracy still requires the original solver's
/// refinement; this constitutive integral does not establish it by itself.
pub struct IntegratedCaloricBoundary<'l, 'a, 'm> {
    layer: &'l CaloricConductionLayer<'a, 'm>,
    geometry: ImplicitBoundaryGeometry,
    flux_geometry: ReservoirLayerGeometry,
    film_resistance_k_w: f64,
    config: LayerConductionConfig,
}

impl<'l, 'a, 'm> IntegratedCaloricBoundary<'l, 'a, 'm> {
    /// Bind the same positive area and centre-to-face distance as the native
    /// observer. The resistance is area-specific, in m² K/W; zero selects a
    /// prescribed material face. No additional film is inferred from htc.
    pub fn new(
        layer: &'l CaloricConductionLayer<'a, 'm>,
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
                "integrated caloric boundary needs a finite nonnegative film and bounded solve controls",
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

    /// Resolve actual material-face temperature and signed outward power.
    /// The returned endpoint derivatives hold the source and geometry fixed.
    /// They are with respect to temperatures; the native trait below chains
    /// the material derivative with the original 1/Cp or 1/Cv.
    pub fn response_at_specific_energy(
        &self,
        specific_energy_j_kg: f64,
        bath_temperature_k: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ReservoirLayerResponse, DepositionError> {
        checkpoint()?;
        let temperature = self
            .layer
            .caloric_law()
            .temperature_k(specific_energy_j_kg)?;
        self.flux_geometry.evaluate_material_to_bath(
            self.layer.conductivity_curve(),
            [temperature, bath_temperature_k],
            self.film_resistance_k_w,
            self.config,
            checkpoint,
        )
    }
}

impl ImplicitBoundaryLaw for IntegratedCaloricBoundary<'_, '_, '_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.layer.caloric_law()
    }

    fn state_at_specific_enthalpy(
        &self,
        specific_energy_j_kg: f64,
        bath_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitBoundaryState, DepositionError> {
        checkpoint()?;
        if geometry != self.geometry {
            return Err(DepositionError::InvalidInput(
                "integrated caloric boundary requires its original area, distance and film",
            ));
        }
        let state = self.layer.state_at_specific_energy(specific_energy_j_kg)?;
        let response = self.flux_geometry.evaluate_material_to_bath(
            self.layer.conductivity_curve(),
            [state.temperature_k, bath_temperature_k],
            self.film_resistance_k_w,
            self.config,
            &mut *checkpoint,
        )?;
        // The original thermal residual applies G*(Tcell-Tbath). Admit that
        // very power against the original layer solve before returning its
        // analytic derivative; never repair the physical energy afterward.
        let applied = response.conductance_w_k * (state.temperature_k - bath_temperature_k);
        let allowance = self.config.absolute_heat_tolerance_w
            + self.config.relative_heat_tolerance
                * applied.abs().max(response.power_to_reservoir_w.abs());
        let derivative = response.endpoint_derivative_w_k[0] * state.temperature_slope_kg_k_j;
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
