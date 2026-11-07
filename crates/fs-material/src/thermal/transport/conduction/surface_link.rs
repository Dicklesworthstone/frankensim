//! Source-resolved interface temperatures on the original finite Cp/Cv stores.

use super::{CaloricConductionLayer, admit_original_link};
use crate::thermal::conduction::{
    LayerConductionConfig, LayerGeometry, SurfaceExchangeState, SurfaceLayerConductionResponse,
};
use fs_thermal::enthalpy::conduction::{
    ConductiveLink,
    implicit::link::{
        ImplicitLinkLaw, ImplicitLinkState, ImplicitSurfaceExchangeLaw,
        ImplicitSurfaceExchangeState,
    },
};
use fs_thermal::enthalpy::{ContactPatch, DepositionError, DepositionSimulation, ThermalLaw};

/// Two finite original caloric bodies coupled through their solved surfaces.
///
/// The original material layers integrate each supported k(T) profile. A single
/// source evaluates the actual ordered energies and both material-face
/// temperatures inside that SAME bounded scalar solve. Complete source and
/// caloric derivatives reach the original implicit chain/graph. No material
/// chart is resampled, and no extra storage or thermal clock is introduced.
///
/// The interface has equal actual covered areas and is source-free,
/// quasi-steady, passive and zero-storage. Planar, cylindrical and spherical
/// material legs use their original exact resistance coordinates. Radial faces
/// are the two material-side traces at ONE shared interface radius: a cylinder
/// retains one common axial length, and a sphere retains its full solid angle.
/// Both actual covered areas must equal that geometric interface area.
///
/// The source owns the complete interface response; no separate authored
/// contact resistance may coexist. Radiation/view assumptions or gas/contact
/// conditions belong to that source. Radial material geometry does not adapt a
/// planar radiation law to separated concentric radiating faces. Such a gap,
/// partial circumferential contact, axial end effects and general enclosures
/// require their separately resolved geometry and original physical owners.
pub struct IntegratedCaloricSurfaceLink<'l, 'a, 'm> {
    layers: [&'l CaloricConductionLayer<'a, 'm>; 2],
    patch: ContactPatch,
    geometry: LayerGeometry,
    source: &'l dyn ImplicitSurfaceExchangeLaw,
    energy_bounds_j_kg: [[f64; 2]; 2],
    surface_bounds_k: [[f64; 2]; 2],
    config: LayerConductionConfig,
}

impl<'l, 'a, 'm> IntegratedCaloricSurfaceLink<'l, 'a, 'm> {
    /// Bind the original cells, material profiles and one interface source.
    ///
    /// Initial admission immediately uses the actual source, including its own
    /// surface support; no intermediate perfect-contact state is required.
    /// Both source caloric objects must exactly match the original ordered
    /// layers. Explicit energy bounds lie inside each original material chart.
    /// Surface bounds are independent and may exclude the bulk temperatures.
    ///
    /// The original patch and BOTH covered areas must match the geometry's
    /// common interface area, including the cylinder's declared axial length.
    /// `LayerGeometry` admits only positive radial traces straddling that same
    /// radius, in either A/B order. Its equivalent distances are resistance
    /// coordinates; they neither replace physical cell volumes nor infer mass.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        link_id: u64,
        layers: [&'l CaloricConductionLayer<'a, 'm>; 2],
        geometry: LayerGeometry,
        source: &'l dyn ImplicitSurfaceExchangeLaw,
        config: LayerConductionConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        let caloric = source.caloric_laws();
        let energy_bounds = source.energy_bounds_j_kg();
        let surface_bounds = source.surface_temperature_bounds_k();
        for side in 0..2 {
            let [lo, hi] = energy_bounds[side];
            let [original_lo, original_hi] = layers[side].energy_bounds_j_kg();
            let [surface_lo, surface_hi] = surface_bounds[side];
            if !core::ptr::addr_eq(
                caloric[side] as *const dyn ThermalLaw,
                layers[side].caloric_law() as *const dyn ThermalLaw,
            ) || !lo.is_finite()
                || !hi.is_finite()
                || lo >= hi
                || lo < original_lo
                || hi > original_hi
                || !surface_lo.is_finite()
                || surface_lo <= 0.0
                || !surface_hi.is_finite()
                || surface_lo > surface_hi
            {
                return Err(DepositionError::InvalidInput(
                    "surface source needs its original ordered calorics and supported energy/temperature bounds",
                ));
            }
        }
        // LayerGeometry already owns the common interface radius/length and
        // ordered radial resistance integrals. Match its physical area before
        // lowering those integrals to the original scalar surface solve.
        let (patch, energies) =
            admit_original_link(thermal, link_id, layers, geometry, &mut checkpoint)?;
        if patch.covered_area_m2 != [patch.link.area_m2; 2] {
            return Err(DepositionError::InvalidInput(
                "caloric surface exchange requires both physical covered areas to equal its area",
            ));
        }
        let link = Self {
            layers,
            patch,
            geometry,
            source,
            energy_bounds_j_kg: energy_bounds,
            surface_bounds_k: surface_bounds,
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

    #[must_use]
    pub const fn geometry(&self) -> LayerGeometry {
        self.geometry
    }

    #[must_use]
    pub const fn layers(&self) -> [&'l CaloricConductionLayer<'a, 'm>; 2] {
        self.layers
    }

    /// Original native geometry; the source owns the entire interface.
    #[must_use]
    pub fn native_geometry(&self) -> ConductiveLink {
        ConductiveLink {
            link_id: self.link_id(),
            half_distances_m: self.geometry.equivalent_distances_m(),
            contact_resistance_m2_k_w: 0.0,
        }
    }

    /// Resolve both surfaces and one conservative original A-to-B transfer.
    ///
    /// Endpoint derivatives hold source energy coordinates fixed and already
    /// include both source surface partials. The native trait additionally
    /// applies original 1/Cp or 1/Cv and explicit source energy partials.
    /// The returned source is the accepted surface query, without re-evaluation.
    pub fn evaluate(
        &self,
        energies_j_kg: [f64; 2],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<SurfaceLayerConductionResponse<ImplicitSurfaceExchangeState>, DepositionError> {
        self.response(energies_j_kg, &mut checkpoint)
            .map(|(response, _)| response)
    }

    fn response(
        &self,
        energies_j_kg: [f64; 2],
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (
            SurfaceLayerConductionResponse<ImplicitSurfaceExchangeState>,
            [f64; 2],
        ),
        DepositionError,
    > {
        checkpoint()?;
        for side in 0..2 {
            let [lo, hi] = self.energy_bounds_j_kg[side];
            if !energies_j_kg[side].is_finite() || !(lo..=hi).contains(&energies_j_kg[side]) {
                return Err(DepositionError::MaterialDomain);
            }
        }
        let states = [
            self.layers[0].state_at_specific_energy(energies_j_kg[0])?,
            self.layers[1].state_at_specific_energy(energies_j_kg[1])?,
        ];
        let temperatures = states.map(|state| state.temperature_k);
        let response = self
            .geometry
            .conduction(
                self.layers.map(CaloricConductionLayer::conductivity_curve),
                0.0,
            )?
            .evaluate_surface_exchange(
                temperatures,
                self.surface_bounds_k,
                self.config,
                |surfaces| {
                    let source = self.source.state_at_surfaces(energies_j_kg, surfaces)?;
                    source.admit()?;
                    Ok(SurfaceExchangeState {
                        flux_a_to_b_w_m2: source.flux_a_to_b_w_m2,
                        conductance_w_m2_k: source.conductance_w_m2_k,
                        surface_derivatives_w_m2_k: source.surface_derivatives_w_m2_k,
                        source,
                    })
                },
                &mut *checkpoint,
            )?;
        // Admit the original native secant lowering without changing its heat.
        let applied = response.conductance_w_k * (temperatures[0] - temperatures[1]);
        let allowance = self.config.absolute_heat_tolerance_w
            + self.config.relative_heat_tolerance
                * applied.abs().max(response.power_a_to_b_w.abs());
        let derivatives = std::array::from_fn(|side| {
            response.source_flux_derivative_scale_m2.mul_add(
                response.source.energy_derivatives_w_kg_m2_j[side],
                response.endpoint_derivative_w_k[side] * states[side].temperature_slope_kg_k_j,
            )
        });
        if !applied.is_finite()
            || !allowance.is_finite()
            || (applied - response.power_a_to_b_w).abs() > allowance
            || derivatives.iter().any(|value| !value.is_finite())
        {
            return Err(DepositionError::Arithmetic);
        }
        checkpoint()?;
        Ok((response, derivatives))
    }
}

impl ImplicitLinkLaw for IntegratedCaloricSurfaceLink<'_, '_, '_> {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
        self.layers.map(|layer| layer.caloric_law())
    }

    fn state_at_specific_enthalpies(
        &self,
        energies_j_kg: [f64; 2],
        geometry: ConductiveLink,
        patch: ContactPatch,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitLinkState, DepositionError> {
        checkpoint()?;
        if geometry != self.native_geometry() || patch != self.patch {
            return Err(DepositionError::InvalidInput(
                "caloric surface exchange requires its original ordered patch and complete native geometry",
            ));
        }
        let (response, derivatives) = self.response(energies_j_kg, checkpoint)?;
        Ok(ImplicitLinkState {
            conductance_w_k: response.conductance_w_k,
            power_enthalpy_derivatives_w_kg_j: derivatives,
        })
    }
}
