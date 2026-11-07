//! One gray conducting gap on the original finite Cp/Cv thermal network.

use super::{BoundCaloricContact, CaloricConductionLayer, admit_original_link};
use crate::thermal::conduction::{
    GrayLayerConductionResponse, GrayLayerGap, LayerConductionConfig, LayerGeometry,
};
use fs_thermal::enthalpy::conduction::{
    ConductiveLink,
    implicit::{
        contact::ImplicitContactLaw,
        link::{ImplicitLinkLaw, ImplicitLinkState},
    },
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
/// unequal-area faces. Emissivities, area and independently authored gap remain
/// fixed; an optional original source supplies its nonradiative resistance.
/// The geometry-aware native query includes density-driven material thickness.
/// Material profiles are source-free and quasi-steady within each half-cell;
/// transient accuracy still requires the original spatial/time refinement.
#[derive(Debug)]
pub struct IntegratedCaloricGrayLink<'l, 'a, 'm> {
    layers: [&'l CaloricConductionLayer<'a, 'm>; 2],
    patch: ContactPatch,
    geometry: LayerGeometry,
    gap: GrayLayerGap,
    config: LayerConductionConfig,
    gap_resistance: Option<BoundCaloricContact<'l>>,
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
        Self::prepare(
            thermal,
            link_id,
            layers,
            geometry,
            gap,
            config,
            None,
            &mut checkpoint,
        )
    }

    /// Bind one original energy-dependent NONRADIATIVE gap resistance.
    ///
    /// The declared gray gap must have zero fixed HTC. Its emissivities remain
    /// fixed, while this source supplies R''(e_A,e_B) at the actual ordered
    /// original bulk h/u pair, with its complete derivatives. Each resistance
    /// must be strictly positive and finite with a representable reciprocal.
    /// It includes neither radiation nor either material half-cell, coating,
    /// or any additional heat source. External source conditions remain fixed
    /// by that law; solved surface temperatures do not replace its coordinates.
    ///
    /// The same source/identity/bounds owner serves ordinary integrated
    /// contacts. Initial admission and every subsequent trial use the actual
    /// source immediately, once per energy pair. Both dR''/de contributions
    /// enter the complete native tangent after the original gray surface solve.
    /// Valid resistance feedback may reverse a differential heat-flow slope
    /// while positive R'' still gives a passive signed transfer. The source
    /// retains its own evaluation costs and introduces no thermal inventory.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_gap_resistance(
        thermal: &DepositionSimulation<'_>,
        link_id: u64,
        layers: [&'l CaloricConductionLayer<'a, 'm>; 2],
        geometry: LayerGeometry,
        gap: GrayLayerGap,
        source: &'l dyn ImplicitContactLaw,
        config: LayerConductionConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        let bound = BoundCaloricContact::new(source, layers)?;
        Self::prepare(
            thermal,
            link_id,
            layers,
            geometry,
            gap,
            config,
            Some(bound),
            &mut checkpoint,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare(
        thermal: &DepositionSimulation<'_>,
        link_id: u64,
        layers: [&'l CaloricConductionLayer<'a, 'm>; 2],
        geometry: LayerGeometry,
        gap: GrayLayerGap,
        config: LayerConductionConfig,
        gap_resistance: Option<BoundCaloricContact<'l>>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if geometry.radial_coordinates_m().is_some() {
            return Err(DepositionError::InvalidInput(
                "integrated caloric gray gap requires equal-area planar material faces",
            ));
        }
        let (patch, energies) =
            admit_original_link(thermal, link_id, layers, geometry, checkpoint)?;
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
            gap_resistance,
        };
        link.evaluate(energies, &mut *checkpoint)?;
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
    /// The response's endpoint derivatives are with respect to temperature,
    /// holding the evaluated gap resistance fixed. The native trait below
    /// chains them through each original caloric inverse and includes the
    /// original gap-resistance derivatives. Gap conduction and radiation are
    /// diagnostics of this SAME transfer; neither is deposited as an additional
    /// thermal source.
    pub fn evaluate(
        &self,
        energies_j_kg: [f64; 2],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<GrayLayerConductionResponse, DepositionError> {
        self.response(energies_j_kg, self.geometry, [0.0; 2], &mut checkpoint)
            .map(|(response, _)| response)
    }

    fn response(
        &self,
        energies_j_kg: [f64; 2],
        geometry: LayerGeometry,
        relative_distance_slopes_kg_j: [f64; 2],
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(GrayLayerConductionResponse, [f64; 2]), DepositionError> {
        checkpoint()?;
        let a = self.layers[0].state_at_specific_energy(energies_j_kg[0])?;
        let b = self.layers[1].state_at_specific_energy(energies_j_kg[1])?;
        let resistance = if let Some(bound) = &self.gap_resistance {
            Some(bound.evaluate(energies_j_kg, checkpoint)?)
        } else {
            None
        };
        let gap = if let Some(state) = resistance {
            self.gap
                .with_nonradiative_resistance(state.resistance_m2_k_w)?
        } else {
            self.gap
        };
        let response = geometry
            .conduction(
                self.layers.map(CaloricConductionLayer::conductivity_curve),
                0.0,
            )?
            .evaluate_gray_gap(
                [a.temperature_k, b.temperature_k],
                gap,
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
        let mut derivatives = [
            response.endpoint_derivative_w_k[0] * a.temperature_slope_kg_k_j,
            response.endpoint_derivative_w_k[1] * b.temperature_slope_kg_k_j,
        ];
        if let Some(state) = resistance {
            for side in 0..2 {
                derivatives[side] = response.gap_resistance_derivative_w2_m2_k.mul_add(
                    state.resistance_enthalpy_derivatives_m2_kg_k_w_j[side],
                    derivatives[side],
                );
            }
        }
        for side in 0..2 {
            let slope = relative_distance_slopes_kg_j[side];
            if slope != 0.0 {
                // The same accepted surface elimination gives
                // Q_dA=-(Q/A)*Q_TA/k(TA), Q_dB=(Q/A)*Q_TB/k(TB).
                // Endpoint partials hold the original gap source fixed; its
                // explicit energy response above must not enter this factor.
                let k = [a.conductivity_w_m_k, b.conductivity_w_m_k][side];
                let sign = if side == 0 { -1.0 } else { 1.0 };
                let temperature_response = (sign
                    * (response.power_a_to_b_w / geometry.contact_area_m2()))
                    * (geometry.equivalent_distances_m()[side] / k)
                    * slope;
                derivatives[side] = response.endpoint_derivative_w_k[side]
                    .mul_add(temperature_response, derivatives[side]);
            }
        }
        if !applied.is_finite()
            || !allowance.is_finite()
            || (applied - response.power_a_to_b_w).abs() > allowance
            || derivatives.iter().any(|value| !value.is_finite())
            || (resistance.is_none()
                && relative_distance_slopes_kg_j == [0.0; 2]
                && (derivatives[0] < 0.0 || derivatives[1] > 0.0))
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
        let (response, derivatives) =
            self.response(energies_j_kg, self.geometry, [0.0; 2], checkpoint)?;
        Ok(ImplicitLinkState {
            conductance_w_k: response.conductance_w_k,
            power_enthalpy_derivatives_w_kg_j: derivatives,
        })
    }

    fn state_at_specific_enthalpies_with_geometry_tangent(
        &self,
        energies_j_kg: [f64; 2],
        geometry: ConductiveLink,
        patch: ContactPatch,
        half_distance_relative_enthalpy_derivatives_kg_j: [f64; 2],
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<ImplicitLinkState>, DepositionError> {
        checkpoint()?;
        if patch != self.patch
            || geometry.link_id != self.link_id()
            || geometry.contact_resistance_m2_k_w != 0.0
            || half_distance_relative_enthalpy_derivatives_kg_j
                .iter()
                .any(|slope| !slope.is_finite())
        {
            return Err(DepositionError::InvalidInput(
                "density-driven caloric gray gap requires its original patch, gap ownership and finite length tangents",
            ));
        }
        let current =
            LayerGeometry::planar(geometry.half_distances_m, self.geometry.contact_area_m2())?;
        let (response, derivatives) = self.response(
            energies_j_kg,
            current,
            half_distance_relative_enthalpy_derivatives_kg_j,
            checkpoint,
        )?;
        Ok(Some(ImplicitLinkState {
            conductance_w_k: response.conductance_w_k,
            power_enthalpy_derivatives_w_kg_j: derivatives,
        }))
    }
}
