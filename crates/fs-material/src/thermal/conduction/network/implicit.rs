//! Source-integrated link tangents for the original implicit thermal owner.
use super::*;
use fs_thermal::enthalpy::conduction::ConductiveLink;
use fs_thermal::enthalpy::conduction::implicit::link::{ImplicitLinkLaw, ImplicitLinkState};

/// Borrowed original source/geometry binding for an implicit chain or graph.
///
/// The same Kirchhoff/contact solve used by `PreparedLayerConduction::resolve`
/// supplies both heat and its full analytic derivative. No source query,
/// geometry conversion, allocation or material identity construction occurs
/// inside a flux evaluation. The enthalpy solver still owns every stored joule,
/// accepted step, entropy account and nonlinear convergence decision.
pub struct ImplicitLayerLink<'p, 'a, 'm> {
    link: &'p PreparedLink<'a, 'm>,
    config: LayerConductionConfig,
}

impl ImplicitLayerLink<'_, '_, '_> {
    pub fn link_id(&self) -> u64 {
        self.link.binding.link
    }
}

impl<'a, 'm> PreparedLayerConduction<'a, 'm> {
    /// Borrow each original link for `ImplicitConductionChain::with_link_laws`
    /// or `ImplicitConductionGraph::with_link_laws`.
    ///
    /// The native conduction observer must carry the exact same equivalent
    /// distances and fixed contact resistance as `bindings()`. Planar and
    /// radial bindings retain their physical metrics and ordered cell coverage.
    /// Attachment probes source support and scalar-solve controls. Other links
    /// and all reservoir faces keep their native bulk-conductivity semantics.
    pub fn implicit_link_laws(
        &self,
        config: LayerConductionConfig,
    ) -> impl ExactSizeIterator<Item = ImplicitLayerLink<'_, 'a, 'm>> {
        self.links.iter().map(move |link| ImplicitLayerLink { link, config })
    }
}

impl ImplicitLinkLaw for ImplicitLayerLink<'_, '_, '_> {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
        self.link.binding.layers.map(|layer| layer.caloric_curve() as &dyn ThermalLaw)
    }

    fn state_at_specific_enthalpies(
        &self,
        enthalpies_j_kg: [f64; 2],
        geometry: ConductiveLink,
        patch: ContactPatch,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitLinkState, DepositionError> {
        checkpoint()?;
        let binding = &self.link.binding;
        if patch != self.link.patch
            || geometry.link_id != binding.link
            || geometry.half_distances_m != binding.distances_m
            || geometry.contact_resistance_m2_k_w != binding.contact_area_resistance_m2_k_w
        {
            return Err(DepositionError::InvalidInput(
                "implicit layers require the original contact, distances and resistance",
            ));
        }
        let mut temperature = [0.0; 2];
        let mut slope = [0.0; 2];
        for side in 0..2 {
            let layer = binding.layers[side];
            let h = enthalpies_j_kg[side];
            temperature[side] = layer.temperature_at_enthalpy(h).map_err(|error| match error {
                LayerSourceError::Thermal(error) => error,
                _ => DepositionError::MaterialDomain,
            })?;
            slope[side] = temperature_slope(layer, h, temperature[side])?;
        }
        let response = self.link.flux_geometry.conduction(
            binding.layers.map(|layer| layer.conductivity_curve()),
            binding.contact_area_resistance_m2_k_w,
        )?.evaluate(temperature, self.config, &mut *checkpoint)?;
        // The native energy owner consumes the secant. Keep its actual applied
        // power within the SAME constitutive tolerance as the existing adapter.
        let applied = response.conductance_w_k * (temperature[0] - temperature[1]);
        let allowance = self.config.absolute_heat_tolerance_w
            + self.config.relative_heat_tolerance * applied.abs().max(response.power_a_to_b_w.abs());
        let derivatives = core::array::from_fn(|side| {
            response.endpoint_derivative_w_k[side] * slope[side]
        });
        if !applied.is_finite()
            || !allowance.is_finite()
            || (applied - response.power_a_to_b_w).abs() > allowance
            || derivatives.iter().any(|value| !value.is_finite())
        {
            return Err(DepositionError::Arithmetic);
        }
        checkpoint()?;
        Ok(ImplicitLinkState {
            conductance_w_k: response.conductance_w_k,
            power_enthalpy_derivatives_w_kg_j: derivatives,
        })
    }
}

fn temperature_slope(layer: &ResolvedConductionLayer<'_>, h: f64, temperature: f64) -> Result<f64, DepositionError> {
    let caloric = layer.caloric_curve();
    let knots = caloric.knots();
    // At a closed upper source endpoint use the admitted left derivative.
    // In particular a solid-only layer must not use the excluded plateau's
    // zero right derivative. Other points keep the original chart convention.
    if temperature == layer.conductivity_curve().temperature_range_k()[1] {
        let hi = knots.partition_point(|k| k.specific_enthalpy_j_kg < h).min(knots.len() - 1);
        let lo = hi.checked_sub(1).ok_or(DepositionError::MaterialDomain)?;
        let slope = (knots[hi].temperature_k - knots[lo].temperature_k)
            / (knots[hi].specific_enthalpy_j_kg - knots[lo].specific_enthalpy_j_kg);
        if !slope.is_finite() || slope <= 0.0 { return Err(DepositionError::Arithmetic); }
        Ok(slope)
    } else {
        caloric.temperature_derivative_at_specific_enthalpy(h)
            .map_err(|_| DepositionError::MaterialDomain)
    }
}
