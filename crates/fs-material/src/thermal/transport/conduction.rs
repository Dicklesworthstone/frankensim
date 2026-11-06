//! Source-integrated heat exchange with the original Cp(T) and Cv(T) storage.
//!
//! The existing Kirchhoff layer kernel supplies one signed heat transfer and
//! its analytic temperature derivatives. This adapter chains those derivatives
//! through each ORIGINAL caloric inverse, retaining h for Cp and u for Cv.
//! It neither resamples a caloric chart nor integrates a second thermal state.

/// Source-integrated material surfaces and their original reservoir/ambient film.
#[path = "conduction/boundary.rs"]
pub mod boundary;
pub use boundary::IntegratedCaloricBoundary;

/// One source-integrated material surface exchanging convection and radiation.
#[path = "conduction/gray_boundary.rs"]
pub mod gray_boundary;
pub use gray_boundary::IntegratedCaloricGrayBoundary;

/// Source-integrated material surfaces joined by one gray conducting gap.
#[path = "conduction/gray_link.rs"]
pub mod gray_link;
pub use gray_link::IntegratedCaloricGrayLink;

/// Original finite caloric bodies joined by a source-resolved two-surface law.
#[path = "conduction/surface_link.rs"]
pub mod surface_link;
pub use surface_link::IntegratedCaloricSurfaceLink;

use super::{
    isobaric::conduction::IsobaricConductionLaw, isochoric::conduction::IsochoricConductionLaw,
};
use crate::thermal::conduction::{
    ConductivityCurve, LayerConductionConfig, LayerConductionResponse, LayerGeometry,
};
use fs_matdb::PropertyLineSlice;
use fs_thermal::enthalpy::conduction::ConductiveLink;
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionLaw, ImplicitConductionState,
    contact::{ImplicitContactLaw, ImplicitContactState},
    link::{ImplicitLinkLaw, ImplicitLinkState},
};
use fs_thermal::enthalpy::{ContactPatch, DepositionError, DepositionSimulation, ThermalLaw};

#[derive(Clone, Copy, Debug)]
enum OriginalConduction<'a, 'm> {
    Isobaric(&'a IsobaricConductionLaw<'m>),
    Isochoric(&'a IsochoricConductionLaw<'m>),
}

/// One original caloric/conductivity source prepared for integrated layer fluxes.
///
/// Preparation copies only its original k(T) knots into the existing bounded
/// conductivity-integral owner. It retains the exact Cp or Cv object and source
/// receipts. Many faces can borrow this same prepared layer; face construction
/// never copies its tables. Bulk trial evaluation allocates nothing and queries
/// no material card; a separately bound contact retains its own evaluation costs.
/// The separate integral-knot bound never narrows the original point law.
#[derive(Debug)]
pub struct CaloricConductionLayer<'a, 'm> {
    original: OriginalConduction<'a, 'm>,
    conductivity: ConductivityCurve,
}

impl<'a, 'm> CaloricConductionLayer<'a, 'm> {
    /// Retain source-exact Cp storage at its original fixed pressure/composition.
    pub fn from_isobaric(law: &'a IsobaricConductionLaw<'m>) -> Result<Self, DepositionError> {
        Ok(Self {
            original: OriginalConduction::Isobaric(law),
            conductivity: law.integrated_conductivity_curve()?,
        })
    }

    /// Retain source-exact Cv storage, fixed density and original cell admission.
    pub fn from_isochoric(law: &'a IsochoricConductionLaw<'m>) -> Result<Self, DepositionError> {
        Ok(Self {
            original: OriginalConduction::Isochoric(law),
            conductivity: law.integrated_conductivity_curve()?,
        })
    }

    /// The exact caloric object registered on the original thermal timeline.
    #[must_use]
    pub fn caloric_law(&self) -> &'m dyn ThermalLaw {
        match self.original {
            OriginalConduction::Isobaric(law) => law.caloric(),
            OriginalConduction::Isochoric(law) => law.caloric(),
        }
    }

    /// Original retained conductivity knots, with no caloric-grid resampling.
    #[must_use]
    pub const fn conductivity_curve(&self) -> &ConductivityCurve {
        &self.conductivity
    }

    /// Original source receipts and uncertainty; preparation adds no authority.
    #[must_use]
    pub fn source_line(&self) -> &PropertyLineSlice {
        match self.original {
            OriginalConduction::Isobaric(law) => law.source_line(),
            OriginalConduction::Isochoric(law) => law.source_line(),
        }
    }

    /// Original specific h or u support, selected by `caloric_law().storage_basis()`.
    #[must_use]
    pub fn energy_bounds_j_kg(&self) -> [f64; 2] {
        match self.original {
            OriginalConduction::Isobaric(law) => law.enthalpy_bounds_j_kg(),
            OriginalConduction::Isochoric(law) => law.enthalpy_bounds_j_kg(),
        }
    }

    /// Original inverse and analytic derivatives with respect to that SAME h or u.
    pub fn state_at_specific_energy(
        &self,
        energy: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        match self.original {
            OriginalConduction::Isobaric(law) => law.state_at_specific_enthalpy(energy),
            OriginalConduction::Isochoric(law) => law.state_at_specific_internal_energy(energy),
        }
        .map_err(|_| DepositionError::MaterialDomain)
    }
}

// Shared admission for every integrated flux on an original finite-body link.
// Probe each concrete interface only AFTER this binding, so finite physical
// jumps do not first have to admit an unrelated perfect-contact solution.
fn admit_original_link(
    thermal: &DepositionSimulation<'_>,
    link_id: u64,
    layers: [&CaloricConductionLayer<'_, '_>; 2],
    geometry: LayerGeometry,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(ContactPatch, [f64; 2]), DepositionError> {
    checkpoint()?;
    let patch = thermal
        .contact_patch(link_id)
        .ok_or(DepositionError::InvalidInput(
            "integrated caloric link requires an active original contact",
        ))?;
    if patch.link.area_m2 != geometry.contact_area_m2() {
        return Err(DepositionError::InvalidInput(
            "integrated caloric geometry must retain the original contact area",
        ));
    }
    let mut energies = [0.0; 2];
    for (side, id) in [patch.link.a, patch.link.b].into_iter().enumerate() {
        checkpoint()?;
        let cell = thermal
            .cells()
            .get(&id)
            .ok_or(DepositionError::InvalidInput(
                "integrated caloric endpoint is not active",
            ))?;
        let caloric = layers[side].caloric_law();
        if !thermal.uses_material_law(cell.material_id, caloric) {
            return Err(DepositionError::InvalidInput(
                "integrated caloric link requires the original ordered caloric objects",
            ));
        }
        caloric.validate_cell(cell)?;
        energies[side] = cell.specific_h_j_kg;
    }
    Ok((patch, energies))
}

/// One bound source-integrated face of the original implicit chain or graph.
///
/// The caller declares homogeneous, source-free, quasi-steady layers aligned
/// with each scalar conductivity and fixed planar/radial geometry. The native
/// timeline owns all finite-cell storage, mass, entropy and accepted time steps.
/// Layers add no separate heat inventory, convection or mechanical work.
/// Either ordered side may use Cp or Cv: heat exchange does not equate h with u.
#[derive(Debug)]
pub struct IntegratedCaloricLink<'l, 'a, 'm> {
    layers: [&'l CaloricConductionLayer<'a, 'm>; 2],
    patch: ContactPatch,
    geometry: LayerGeometry,
    resistance_m2_k_w: f64,
    config: LayerConductionConfig,
    contact: Option<BoundCaloricContact<'l>>,
}

struct BoundCaloricContact<'a> {
    law: &'a dyn ImplicitContactLaw,
    bounds_j_kg: [[f64; 2]; 2],
}

impl core::fmt::Debug for BoundCaloricContact<'_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("BoundCaloricContact")
            .field("bounds_j_kg", &self.bounds_j_kg)
            .finish_non_exhaustive()
    }
}

impl<'a> BoundCaloricContact<'a> {
    fn new(
        law: &'a dyn ImplicitContactLaw,
        layers: [&CaloricConductionLayer<'_, '_>; 2],
    ) -> Result<Self, DepositionError> {
        let caloric = law.caloric_laws();
        let bounds_j_kg = law.enthalpy_bounds_j_kg();
        for side in 0..2 {
            let [lo, hi] = bounds_j_kg[side];
            let original = layers[side].caloric_law();
            let [source_lo, source_hi] = layers[side].energy_bounds_j_kg();
            if !lo.is_finite()
                || !hi.is_finite()
                || lo >= hi
                || lo < source_lo
                || hi > source_hi
                || !core::ptr::addr_eq(
                    original as *const dyn ThermalLaw,
                    caloric[side] as *const dyn ThermalLaw,
                )
            {
                return Err(DepositionError::InvalidInput(
                    "integrated contact requires original ordered calorics and supported energy bounds",
                ));
            }
        }
        Ok(Self { law, bounds_j_kg })
    }

    fn evaluate(
        &self,
        energies_j_kg: [f64; 2],
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitContactState, DepositionError> {
        for (energy, [lo, hi]) in energies_j_kg.into_iter().zip(self.bounds_j_kg) {
            if !(lo..=hi).contains(&energy) {
                return Err(DepositionError::MaterialDomain);
            }
        }
        checkpoint()?;
        let state = self.law.state_at_specific_enthalpies(energies_j_kg)?;
        state.admit()?;
        checkpoint()?;
        Ok(state)
    }
}

impl<'l, 'a, 'm> IntegratedCaloricLink<'l, 'a, 'm> {
    /// Bind the actual ordered caloric objects, contact coverage and geometry.
    ///
    /// `geometry` must use the original link's nominal contact area. Its exact
    /// equivalent distances are returned by `native_geometry()` for the existing
    /// conduction observer, and checked again on every native flux query.
    /// Source density/mass constraints, initial integrated flux support, scalar
    /// work limits and cancellation are admitted before returning any binding.
    /// No accepted physical state changes. Reprepare after geometry or contact
    /// changes; source tables are shared by reference, never copied per link.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        link_id: u64,
        layers: [&'l CaloricConductionLayer<'a, 'm>; 2],
        geometry: LayerGeometry,
        contact_resistance_m2_k_w: f64,
        config: LayerConductionConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        Self::prepare(
            thermal,
            link_id,
            layers,
            geometry,
            contact_resistance_m2_k_w,
            config,
            None,
            &mut checkpoint,
        )
    }

    /// Bind one original state-dependent contact to the integrated material legs.
    ///
    /// The contact must retain the same ordered caloric objects and declare
    /// finite h/u bounds inside their original support. Its resistance replaces
    /// the contact once; `native_geometry()` therefore declares zero fixed
    /// resistance. Attach this complete link with `with_link_laws`, without a
    /// separate contact attachment or additional authored resistance.
    ///
    /// Every trial queries the actual resistance and its derivatives, then
    /// applies dQ/de = (dQ/dT)(dT/de) + (dQ/dR'')(dR''/de) on each side.
    /// Temperature-dependent contact can reverse a local differential slope
    /// while its nonnegative resistance still gives passive signed heat flow.
    /// The original contact owns its external conditions and evaluation costs.
    ///
    /// Initial admission uses this actual contact immediately. Disjoint source
    /// temperature spans can be joined by a finite supported contact jump;
    /// construction never first demands an unsupported perfect-contact state.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_contact(
        thermal: &DepositionSimulation<'_>,
        link_id: u64,
        layers: [&'l CaloricConductionLayer<'a, 'm>; 2],
        geometry: LayerGeometry,
        contact: &'l dyn ImplicitContactLaw,
        config: LayerConductionConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        let contact = BoundCaloricContact::new(contact, layers)?;
        Self::prepare(
            thermal,
            link_id,
            layers,
            geometry,
            0.0,
            config,
            Some(contact),
            &mut checkpoint,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare(
        thermal: &DepositionSimulation<'_>,
        link_id: u64,
        layers: [&'l CaloricConductionLayer<'a, 'm>; 2],
        geometry: LayerGeometry,
        contact_resistance_m2_k_w: f64,
        config: LayerConductionConfig,
        contact: Option<BoundCaloricContact<'l>>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        let (patch, energies) =
            admit_original_link(thermal, link_id, layers, geometry, checkpoint)?;
        let link = Self {
            layers,
            patch,
            geometry,
            resistance_m2_k_w: contact_resistance_m2_k_w,
            config,
            contact,
        };
        link.evaluate(energies, &mut *checkpoint)?;
        checkpoint()?;
        Ok(link)
    }

    #[must_use]
    pub const fn link_id(&self) -> u64 {
        self.patch.link.id
    }

    /// Bound physical geometry; radial distances retain their geometric meaning.
    #[must_use]
    pub const fn geometry(&self) -> LayerGeometry {
        self.geometry
    }

    /// Original prepared material layers in the contact's A/B order.
    #[must_use]
    pub const fn layers(&self) -> [&'l CaloricConductionLayer<'a, 'm>; 2] {
        self.layers
    }

    /// Exact geometry supplied to the native `ConductionObserver`. A bound
    /// contact source owns its resistance once, so the fixed entry stays zero.
    #[must_use]
    pub fn native_geometry(&self) -> ConductiveLink {
        ConductiveLink {
            link_id: self.link_id(),
            half_distances_m: self.geometry.equivalent_distances_m(),
            contact_resistance_m2_k_w: self.resistance_m2_k_w,
        }
    }

    /// Full original layer/contact response at supported specific h/u endpoints.
    /// Positive power leaves A and enters B; no thermal history is changed.
    /// Response derivatives hold contact resistance fixed. The `ImplicitLinkLaw`
    /// response additionally includes the original contact's dR''/de chain rule.
    pub fn evaluate(
        &self,
        energies_j_kg: [f64; 2],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<LayerConductionResponse, DepositionError> {
        self.response(energies_j_kg, &mut checkpoint)
            .map(|(response, _)| response)
    }

    fn response(
        &self,
        energies_j_kg: [f64; 2],
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(LayerConductionResponse, [f64; 2]), DepositionError> {
        checkpoint()?;
        let a = self.layers[0].state_at_specific_energy(energies_j_kg[0])?;
        let b = self.layers[1].state_at_specific_energy(energies_j_kg[1])?;
        let contact = if let Some(bound) = &self.contact {
            Some(bound.evaluate(energies_j_kg, checkpoint)?)
        } else {
            None
        };
        let resistance = contact.map_or(self.resistance_m2_k_w, |state| state.resistance_m2_k_w);
        let response = self
            .geometry
            .conduction(
                self.layers.map(CaloricConductionLayer::conductivity_curve),
                resistance,
            )?
            .evaluate(
                [a.temperature_k, b.temperature_k],
                self.config,
                &mut *checkpoint,
            )?;
        // The native residual applies G*(Ta-Tb). Keep that lowering inside the
        // same heat-continuity allowance as the actual integrated layer flux.
        let applied = response.conductance_w_k * (a.temperature_k - b.temperature_k);
        let allowance = self.config.absolute_heat_tolerance_w
            + self.config.relative_heat_tolerance
                * applied.abs().max(response.power_a_to_b_w.abs());
        let mut derivatives = [
            response.endpoint_derivative_w_k[0] * a.temperature_slope_kg_k_j,
            response.endpoint_derivative_w_k[1] * b.temperature_slope_kg_k_j,
        ];
        if let Some(state) = contact {
            for side in 0..2 {
                derivatives[side] = response.contact_resistance_derivative_w2_m2_k.mul_add(
                    state.resistance_enthalpy_derivatives_m2_kg_k_w_j[side],
                    derivatives[side],
                );
            }
        }
        if !applied.is_finite()
            || !allowance.is_finite()
            || (applied - response.power_a_to_b_w).abs() > allowance
            || derivatives.iter().any(|value| !value.is_finite())
            || (contact.is_none()
                && (!super::positive(derivatives[0]) || !super::positive(-derivatives[1])))
        {
            return Err(DepositionError::Arithmetic);
        }
        checkpoint()?;
        Ok((response, derivatives))
    }
}

impl ImplicitLinkLaw for IntegratedCaloricLink<'_, '_, '_> {
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
                "integrated caloric link requires its original geometry and contact coverage",
            ));
        }
        let (response, derivatives) = self.response(energies_j_kg, checkpoint)?;
        Ok(ImplicitLinkState {
            conductance_w_k: response.conductance_w_k,
            power_enthalpy_derivatives_w_kg_j: derivatives,
        })
    }
}

#[cfg(test)]
#[path = "conduction/tests.rs"]
mod tests;
