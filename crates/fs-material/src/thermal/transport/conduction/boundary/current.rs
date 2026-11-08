//! One explicitly prepared current face on the original caloric boundary.

use super::IntegratedCaloricBoundary;
use crate::thermal::conduction::ReservoirLayerResponse;
use fs_thermal::enthalpy::conduction::CurrentConductiveBoundary;
use fs_thermal::enthalpy::conduction::implicit::boundary::{
    ImplicitBoundaryGeometry, ImplicitBoundaryLaw, ImplicitBoundaryState,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

/// A current planar area and normal distance with the original material film.
///
/// Preparation rebuilds the original geometric resistance AND converts the
/// original areal film resistance with the CURRENT area, exactly once. Both
/// caches then belong to this explicitly selected mechanical configuration.
/// Every trial still checks its supplied geometry and uses the original
/// conductivity integral, bounded surface solve and exact caloric tangent.
///
/// The original source, bath support, storage basis and fixed external context
/// stay intact. This supplies no density, pressure work, face storage or motion
/// model. Rebind between configurations; geometry stays fixed during the native
/// implicit interval. Density-driven geometry remains a separate query/model.
pub struct CurrentPlanarCaloricBoundary<'r, 'l, 'a, 'm> {
    reference: &'r IntegratedCaloricBoundary<'l, 'a, 'm>,
    prepared: IntegratedCaloricBoundary<'l, 'a, 'm>,
}

impl<'l, 'a, 'm> IntegratedCaloricBoundary<'l, 'a, 'm> {
    /// Prepare current measures while retaining the original areal film.
    ///
    /// This is bounded constant work without source queries or table copies.
    /// The native boundary attachment probes the actual cell/bath state before
    /// returning a solver. The same original caloric object and reference
    /// reservoir/exposed coverage remain owned by the native thermal binding.
    /// Source density/volume constraints require compatible mechanical motion;
    /// changing face measures does not authorize a new material context.
    pub fn with_current_planar_geometry<'r>(
        &'r self,
        current: CurrentConductiveBoundary,
    ) -> Result<CurrentPlanarCaloricBoundary<'r, 'l, 'a, 'm>, DepositionError> {
        let geometry = ImplicitBoundaryGeometry {
            area_m2: current.area_m2,
            half_distance_m: current.half_distance_m,
            surface_resistance_m2_k_w: self.geometry.surface_resistance_m2_k_w,
        };
        // Reuse the original preparation, including representability of R''/A.
        // Retaining the old total film resistance would silently change the
        // physical film whenever the current area differs from the reference.
        let prepared = Self::new(self.layer, geometry, self.config)?;
        Ok(CurrentPlanarCaloricBoundary {
            reference: self,
            prepared,
        })
    }
}

impl CurrentPlanarCaloricBoundary<'_, '_, '_, '_> {
    #[must_use]
    pub fn reference_geometry(&self) -> ImplicitBoundaryGeometry {
        self.reference.geometry
    }

    /// The exact measures to supply to the native current-geometry map.
    #[must_use]
    pub fn current_geometry(&self) -> CurrentConductiveBoundary {
        CurrentConductiveBoundary {
            area_m2: self.prepared.geometry.area_m2,
            half_distance_m: self.prepared.geometry.half_distance_m,
        }
    }

    /// Complete current query geometry, with the unchanged authored areal film.
    #[must_use]
    pub fn native_geometry(&self) -> ImplicitBoundaryGeometry {
        self.prepared.geometry
    }

    /// Actual material-face temperature and signed heat in this configuration.
    pub fn response_at_specific_energy(
        &self,
        specific_energy_j_kg: f64,
        bath_temperature_k: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ReservoirLayerResponse, DepositionError> {
        self.prepared.response_at_specific_energy(
            specific_energy_j_kg,
            bath_temperature_k,
            checkpoint,
        )
    }
}

impl ImplicitBoundaryLaw for CurrentPlanarCaloricBoundary<'_, '_, '_, '_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.prepared.caloric_law()
    }

    fn supports_current_geometry(&self) -> bool {
        true
    }

    fn state_at_specific_enthalpy(
        &self,
        specific_energy_j_kg: f64,
        bath_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitBoundaryState, DepositionError> {
        // The original owner checks this exact prepared geometry, computes the
        // full integrated flux/tangent and admits the same applied secant heat.
        self.prepared.state_at_specific_enthalpy(
            specific_energy_j_kg,
            bath_temperature_k,
            geometry,
            checkpoint,
        )
    }
}
