//! Original Cp/Cv storage with source-integrated material-to-bath heat flow.
//!
//! The existing reservoir-layer owner solves the material surface temperature
//! and film jump. The native thermal owner retains the heat, bath entropy and
//! time integration. This adapter supplies the exact derivative with respect
//! to the same original h or u; it adds no surface heat inventory.

#[path = "boundary/current.rs"]
mod current;
pub use current::CurrentPlanarCaloricBoundary;

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
/// the half-cell. Area, film, composition and the original Cp/Cv constraint
/// remain fixed. The explicit geometry-aware query admits density-driven normal
/// half-thickness while preserving this same source-integrated surface solve.
/// Spatial/time accuracy still requires the original solver's
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

    fn implicit_state(
        &self,
        specific_energy_j_kg: f64,
        bath_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        relative_distance_slope_kg_j: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitBoundaryState, DepositionError> {
        let state = self.layer.state_at_specific_energy(specific_energy_j_kg)?;
        let flux_geometry = if geometry == self.geometry {
            self.flux_geometry
        } else {
            ReservoirLayerGeometry::planar(geometry.area_m2, geometry.half_distance_m)?
        };
        let response = flux_geometry.evaluate_material_to_bath(
            self.layer.conductivity_curve(),
            [state.temperature_k, bath_temperature_k],
            self.film_resistance_k_w,
            self.config,
            &mut *checkpoint,
        )?;
        // Admit the SAME heat applied by the native residual before publishing
        // its complete derivative. No flux or stored energy is repaired.
        let applied = response.conductance_w_k * (state.temperature_k - bath_temperature_k);
        let allowance = self.config.absolute_heat_tolerance_w
            + self.config.relative_heat_tolerance
                * applied.abs().max(response.power_to_reservoir_w.abs());
        let mut derivative = response.endpoint_derivative_w_k[0] * state.temperature_slope_kg_k_j;
        if relative_distance_slope_kg_j != 0.0 {
            // Q=(A/d)*integral(Ts..T) k dT; Ts=Tbath+Rfilm*Q.
            // The geometry term is -Q*(d'/d)/(1+(A/d)*k(Ts)*Rfilm).
            // Resolve k at the SAME accepted material face, not at the bath.
            let kface = self
                .layer
                .conductivity_curve()
                .conductivity_w_m_k(response.material_face_temperature_k)?;
            let material_resistance = geometry.half_distance_m / kface;
            let total_resistance = material_resistance + geometry.surface_resistance_m2_k_w;
            if !material_resistance.is_finite()
                || material_resistance <= 0.0
                || !total_resistance.is_finite()
                || total_resistance <= 0.0
            {
                return Err(DepositionError::Arithmetic);
            }
            let geometry_term =
                -response.power_to_reservoir_w * (material_resistance / total_resistance);
            derivative = geometry_term.mul_add(relative_distance_slope_kg_j, derivative);
        }
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
        self.implicit_state(
            specific_energy_j_kg,
            bath_temperature_k,
            geometry,
            0.0,
            checkpoint,
        )
    }

    fn state_at_specific_enthalpy_with_geometry_tangent(
        &self,
        specific_energy_j_kg: f64,
        bath_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        half_distance_relative_enthalpy_derivative_kg_j: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<ImplicitBoundaryState>, DepositionError> {
        checkpoint()?;
        if geometry.area_m2 != self.geometry.area_m2
            || geometry.surface_resistance_m2_k_w != self.geometry.surface_resistance_m2_k_w
            || !half_distance_relative_enthalpy_derivative_kg_j.is_finite()
        {
            return Err(DepositionError::InvalidInput(
                "density-driven integrated boundary requires its original area, film and finite length tangent",
            ));
        }
        self.implicit_state(
            specific_energy_j_kg,
            bath_temperature_k,
            geometry,
            half_distance_relative_enthalpy_derivative_kg_j,
            checkpoint,
        )
        .map(Some)
    }
}
