//! One gray conducting gap on the original finite Cp/Cv thermal network.

use super::{CaloricConductionLayer, admit_original_link};
use crate::thermal::conduction::{
    GrayLayerConductionResponse, GrayLayerGap, LayerConductionConfig, LayerGeometry,
};
use fs_thermal::enthalpy::conduction::{
    ConductiveLink,
    implicit::link::{ImplicitLinkLaw, ImplicitLinkState},
};
use fs_thermal::enthalpy::{ContactPatch, DepositionError, DepositionSimulation, ThermalLaw};

/// Two original caloric stores coupled through their actual material surfaces.
///
/// Each planar half-cell uses its own source-integrated k(T) profile. The
/// existing gray-gap kernel solves both surface temperatures and the ONE heat
/// transfer, with radiation and gap conduction acting in parallel between
/// those surfaces. The native thermal owner retains both finite masses, their
/// original h/Cp or u/Cv storage, heat/entropy accounts and accepted time steps.
/// Neither the material legs nor the gap introduce another thermal inventory.
///
/// `GrayLayerGap::equal_area_full_view` explicitly supplies the equal-area,
/// opaque diffuse-gray, nonparticipating-gap/full-view model. This adapter
/// admits planar geometry only; it does not infer view factors for curved or
/// unequal-area faces. Gap coefficients and physical geometry remain fixed.
/// Material profiles are source-free and quasi-steady within each half-cell;
/// transient accuracy still requires the original spatial/time refinement.
#[derive(Debug)]
pub struct IntegratedCaloricGrayLink<'l, 'a, 'm> {
    layers: [&'l CaloricConductionLayer<'a, 'm>; 2],
    patch: ContactPatch,
    geometry: LayerGeometry,
    gap: GrayLayerGap,
    config: LayerConductionConfig,
}

impl<'l, 'a, 'm> IntegratedCaloricGrayLink<'l, 'a, 'm> {
    /// Bind the original ordered calorics, contact coverage and physical gap.
    ///
    /// Preparation checks the initial flux using this actual gap immediately.
    /// It does not demand an intermediate perfect-contact state: distinct
    /// material source temperature spans may be bridged by a finite gap while
    /// each actual material surface remains on its own supported chart.
    ///
    /// Attach using the native chain/graph's `with_link_laws` and the returned
    /// `native_geometry()`. Its fixed contact resistance is zero because the
    /// gap owns the complete interface. A separate contact law or additional
    /// authored resistance would count the same interface twice and refuses.
    /// Preparation and trial evaluation allocate no tables; every link borrows
    /// the same prepared material layers and their original source receipts.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        link_id: u64,
        layers: [&'l CaloricConductionLayer<'a, 'm>; 2],
        geometry: LayerGeometry,
        gap: GrayLayerGap,
        config: LayerConductionConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if geometry.radial_coordinates_m().is_some() {
            return Err(DepositionError::InvalidInput(
                "integrated caloric gray gap requires equal-area planar material faces",
            ));
        }
        let (patch, energies) =
            admit_original_link(thermal, link_id, layers, geometry, &mut checkpoint)?;
        if patch.covered_area_m2 != [patch.link.area_m2; 2] {
            return Err(DepositionError::InvalidInput(
                "integrated caloric gray gap requires both physical covered areas to equal its area",
            ));
        }
        let link = Self {
            layers,
            patch,
            geometry,
            gap,
            config,
        };
        link.evaluate(energies, &mut checkpoint)?;
        checkpoint()?;
        Ok(link)
    }

    #[must_use]
    pub const fn link_id(&self) -> u64 {
        self.patch.link.id
    }

    /// Original planar lengths and common physical surface area.
    #[must_use]
    pub const fn geometry(&self) -> LayerGeometry {
        self.geometry
    }

    /// Original prepared material layers in the contact's A/B order.
    #[must_use]
    pub const fn layers(&self) -> [&'l CaloricConductionLayer<'a, 'm>; 2] {
        self.layers
    }

    /// Native geometry with zero additional interface resistance.
    #[must_use]
    pub fn native_geometry(&self) -> ConductiveLink {
        ConductiveLink {
            link_id: self.link_id(),
            half_distances_m: self.geometry.equivalent_distances_m(),
            contact_resistance_m2_k_w: 0.0,
        }
    }

    /// Resolved material surfaces and the components of one conservative flux.
    ///
    /// The response's endpoint derivatives are with respect to temperature.
    /// The native trait below chains them through each original caloric
    /// inverse. Gap conduction and radiation are diagnostics of this SAME
    /// transfer; neither is deposited as an additional thermal source.
    pub fn evaluate(
        &self,
        energies_j_kg: [f64; 2],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<GrayLayerConductionResponse, DepositionError> {
        self.response(energies_j_kg, &mut checkpoint)
            .map(|(response, _)| response)
    }

    fn response(
        &self,
        energies_j_kg: [f64; 2],
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(GrayLayerConductionResponse, [f64; 2]), DepositionError> {
        checkpoint()?;
        let a = self.layers[0].state_at_specific_energy(energies_j_kg[0])?;
        let b = self.layers[1].state_at_specific_energy(energies_j_kg[1])?;
        let response = self
            .geometry
            .conduction(
                self.layers.map(CaloricConductionLayer::conductivity_curve),
                0.0,
            )?
            .evaluate_gray_gap(
                [a.temperature_k, b.temperature_k],
                self.gap,
                self.config,
                &mut *checkpoint,
            )?;
        // Apply the same integral flux through the original native secant.
        // The original kernel admits the physical layer and surface balances;
        // this gate admits its lowering without changing the heat afterward.
        let applied = response.conductance_w_k * (a.temperature_k - b.temperature_k);
        let allowance = self.config.absolute_heat_tolerance_w
            + self.config.relative_heat_tolerance
                * applied.abs().max(response.power_a_to_b_w.abs());
        let derivatives = [
            response.endpoint_derivative_w_k[0] * a.temperature_slope_kg_k_j,
            response.endpoint_derivative_w_k[1] * b.temperature_slope_kg_k_j,
        ];
        if !applied.is_finite()
            || !allowance.is_finite()
            || (applied - response.power_a_to_b_w).abs() > allowance
            || derivatives.iter().any(|value| !value.is_finite())
            || derivatives[0] < 0.0
            || derivatives[1] > 0.0
        {
            return Err(DepositionError::Arithmetic);
        }
        checkpoint()?;
        Ok((response, derivatives))
    }
}

impl ImplicitLinkLaw for IntegratedCaloricGrayLink<'_, '_, '_> {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
        self.layers.map(CaloricConductionLayer::caloric_law)
    }

    fn state_at_specific_enthalpies(
        &self,
        energies_j_kg: [f64; 2],
        geometry: ConductiveLink,
        patch: ContactPatch,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitLinkState, DepositionError> {
        checkpoint()?;
        if patch != self.patch || geometry != self.native_geometry() {
            return Err(DepositionError::InvalidInput(
                "integrated caloric gray gap requires its original geometry and contact coverage",
            ));
        }
        let (response, derivatives) = self.response(energies_j_kg, checkpoint)?;
        Ok(ImplicitLinkState {
            conductance_w_k: response.conductance_w_k,
            power_enthalpy_derivatives_w_kg_j: derivatives,
        })
    }
}
