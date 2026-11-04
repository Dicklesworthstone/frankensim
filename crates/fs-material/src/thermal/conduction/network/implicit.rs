//! Source-integrated link tangents for the original implicit thermal owner.
use super::*;
use fs_thermal::enthalpy::conduction::ConductiveLink;
use fs_thermal::enthalpy::conduction::implicit::contact::ImplicitContactLaw;
use fs_thermal::enthalpy::conduction::implicit::link::{ImplicitLinkLaw, ImplicitLinkState};

/// Borrowed original source/geometry binding for an implicit chain or graph.
///
/// The same Kirchhoff/contact solve used by `PreparedLayerConduction::resolve`
/// supplies both heat and its full analytic derivative. The bulk kernel performs
/// no source query, geometry conversion, allocation or identity construction
/// during evaluation. An optional contact retains its original source response
/// and evaluation costs. The enthalpy solver still owns every stored joule,
/// accepted step, entropy account and nonlinear convergence decision.
pub struct ImplicitLayerLink<'p, 'a, 'm> {
    link: &'p PreparedLink<'a, 'm>,
    config: LayerConductionConfig,
    contact: Option<BoundLayerContact<'p>>,
}

struct BoundLayerContact<'a> {
    law: &'a dyn ImplicitContactLaw,
    bounds_j_kg: [[f64; 2]; 2],
}

impl<'p, 'a, 'm> ImplicitLayerLink<'p, 'a, 'm> {
    pub fn link_id(&self) -> u64 {
        self.link.binding.link
    }

    /// Compose the original state-dependent contact with integrated bulk k(T).
    ///
    /// The prepared geometry must declare zero fixed contact resistance. This
    /// source owns the entire contact once, and must use the same ordered
    /// caloric objects as the two layers. Declare its enthalpy bounds entirely
    /// inside the supported pure-phase layer domains, so the original contact
    /// owner supplies tangents on the same admitted branch. External pressure
    /// and composition remain fixed by that original immutable law.
    ///
    /// Attach the resulting complete link with `with_link_laws`; do not also
    /// attach the contact separately. Attachment to a thermal solver probes the
    /// actual initial state. No heat or history is changed by this composition.
    pub fn with_contact(
        mut self,
        law: &'p dyn ImplicitContactLaw,
    ) -> Result<Self, LayerSourceError> {
        if self.contact.is_some() || self.link.binding.contact_area_resistance_m2_k_w != 0.0 {
            return Err(LayerSourceError::Invalid(
                "integrated contact needs one source owner and zero authored resistance",
            ));
        }
        let caloric = law.caloric_laws();
        let bounds_j_kg = law.enthalpy_bounds_j_kg();
        for side in 0..2 {
            let [lo, hi] = bounds_j_kg[side];
            let original = self.link.binding.layers[side].caloric_curve() as &dyn ThermalLaw;
            if !lo.is_finite()
                || !hi.is_finite()
                || lo >= hi
                || !std::ptr::addr_eq(
                    original as *const dyn ThermalLaw,
                    caloric[side] as *const dyn ThermalLaw,
                )
            {
                return Err(LayerSourceError::Invalid(
                    "integrated contact requires finite bounds and original ordered caloric objects",
                ));
            }
            // Preparation fixes one common admissible branch. A wider contact
            // chart could select a derivative into an excluded latent interval.
            self.link.binding.layers[side].temperature_at_enthalpy(lo)?;
            self.link.binding.layers[side].temperature_at_enthalpy(hi)?;
        }
        self.contact = Some(BoundLayerContact { law, bounds_j_kg });
        Ok(self)
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
    /// retain their native coefficients; `ImplicitLayerBoundary` independently
    /// attaches source-integrated reservoir and ambient faces.
    pub fn implicit_link_laws(
        &self,
        config: LayerConductionConfig,
    ) -> impl ExactSizeIterator<Item = ImplicitLayerLink<'_, 'a, 'm>> {
        self.links.iter().map(move |link| ImplicitLayerLink {
            link,
            config,
            contact: None,
        })
    }
}

impl ImplicitLinkLaw for ImplicitLayerLink<'_, '_, '_> {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
        self.link
            .binding
            .layers
            .map(|layer| layer.caloric_curve() as &dyn ThermalLaw)
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
            temperature[side] = layer
                .temperature_at_enthalpy(h)
                .map_err(|error| match error {
                    LayerSourceError::Thermal(error) => error,
                    _ => DepositionError::MaterialDomain,
                })?;
            slope[side] = if self
                .contact
                .as_ref()
                .is_some_and(|contact| h == contact.bounds_j_kg[side][1])
            {
                inward_temperature_slope(layer, h)?
            } else {
                temperature_slope(layer, h, temperature[side])?
            };
        }
        let contact = if let Some(bound) = &self.contact {
            for (h, [lo, hi]) in enthalpies_j_kg.into_iter().zip(bound.bounds_j_kg) {
                if !(lo..=hi).contains(&h) {
                    return Err(DepositionError::MaterialDomain);
                }
            }
            checkpoint()?;
            let state = bound.law.state_at_specific_enthalpies(enthalpies_j_kg)?;
            state.admit()?;
            checkpoint()?;
            Some(state)
        } else {
            None
        };
        let resistance = contact.map_or(binding.contact_area_resistance_m2_k_w, |state| {
            state.resistance_m2_k_w
        });
        let response = self
            .link
            .flux_geometry
            .conduction(
                binding.layers.map(|layer| layer.conductivity_curve()),
                resistance,
            )?
            .evaluate(temperature, self.config, &mut *checkpoint)?;
        // The native energy owner consumes the secant. Keep its actual applied
        // power within the SAME constitutive tolerance as the existing adapter.
        let applied = response.conductance_w_k * (temperature[0] - temperature[1]);
        let allowance = self.config.absolute_heat_tolerance_w
            + self.config.relative_heat_tolerance
                * applied.abs().max(response.power_a_to_b_w.abs());
        let derivatives = core::array::from_fn(|side| {
            let bulk = response.endpoint_derivative_w_k[side] * slope[side];
            // Differentiate the SAME integrated flux: dQ/dh = dQ/dT*dT/dh
            // + dQ/dR''*dR''/dh. A secant-G derivative would miss the two
            // interface-temperature responses when k varies through a layer.
            contact.map_or(bulk, |state| {
                response.contact_resistance_derivative_w2_m2_k.mul_add(
                    state.resistance_enthalpy_derivatives_m2_kg_k_w_j[side],
                    bulk,
                )
            })
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

pub(super) fn temperature_slope(
    layer: &ResolvedConductionLayer<'_>,
    h: f64,
    temperature: f64,
) -> Result<f64, DepositionError> {
    if temperature == layer.conductivity_curve().temperature_range_k()[1] {
        inward_temperature_slope(layer, h)
    } else {
        layer
            .caloric_curve()
            .temperature_derivative_at_specific_enthalpy(h)
            .map_err(|_| DepositionError::MaterialDomain)
    }
}

fn inward_temperature_slope(
    layer: &ResolvedConductionLayer<'_>,
    h: f64,
) -> Result<f64, DepositionError> {
    let knots = layer.caloric_curve().knots();
    // Use the exact admitted left segment, including a stricter contact upper
    // bound at an interior caloric knot. No finite-difference tangent or change
    // to the actual constitutive state is introduced.
    let hi = knots
        .partition_point(|k| k.specific_enthalpy_j_kg < h)
        .min(knots.len() - 1);
    let lo = hi.checked_sub(1).ok_or(DepositionError::MaterialDomain)?;
    let slope = (knots[hi].temperature_k - knots[lo].temperature_k)
        / (knots[hi].specific_enthalpy_j_kg - knots[lo].specific_enthalpy_j_kg);
    if !slope.is_finite() || slope <= 0.0 {
        return Err(DepositionError::Arithmetic);
    }
    Ok(slope)
}
