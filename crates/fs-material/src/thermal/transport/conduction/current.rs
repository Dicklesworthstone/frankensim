//! Explicitly prepared mechanical measures on the original integrated face.
//!
//! Preparation retains the reference contact and original source objects, then
//! resolves the same source-integrated flux at the selected current measures.
//! The native implicit owner holds those measures throughout a physical interval.
//! A new mechanical configuration needs a new binding, not a coefficient mutation.

use super::{IntegratedCaloricLink, admit_original_link};
use crate::thermal::conduction::{LayerConductionResponse, LayerGeometry};
use fs_thermal::enthalpy::conduction::implicit::link::{ImplicitLinkLaw, ImplicitLinkState};
use fs_thermal::enthalpy::conduction::{ConductiveLink, CurrentConductiveLink};
use fs_thermal::enthalpy::{ContactPatch, DepositionError, DepositionSimulation, ThermalLaw};

/// One current planar configuration borrowing its original integrated link.
///
/// The current common area and both physical normal distances are explicit.
/// Each query retains the original ordered caloric objects, reference coverage,
/// authored areal resistance or single source-contact owner, and source support.
/// The original Kirchhoff solve supplies the flux and complete energy tangents;
/// no material table is copied, source is re-resolved, or heat store is added.
///
/// Geometry is fixed during thermal differentiation. This is separate from
/// density-driven thickness: no dA/de or dd/de is inferred from these measures.
/// Source pressure, composition and caloric constraints remain authoritative;
/// fixed-density/volume sources require independently compatible mechanical
/// motion. Face measures alone do not establish a new density or work law.
#[derive(Debug)]
pub struct CurrentPlanarCaloricLink<'r, 'l, 'a, 'm> {
    reference: &'r IntegratedCaloricLink<'l, 'a, 'm>,
    geometry: LayerGeometry,
}

impl<'l, 'a, 'm> IntegratedCaloricLink<'l, 'a, 'm> {
    /// Prepare a selected current planar interval without changing this binding.
    ///
    /// Build the native observer from this original link's `native_geometry()`,
    /// then supply the same current measures through `with_current_geometry`.
    /// Attach the returned law to that observer's implicit chain or graph.
    /// The actual current face is probed at the supplied thermal state before
    /// return; reference identity, coverage and source cell admission are checked.
    ///
    /// Repeating this operation prepares a later mechanical configuration while
    /// retaining the original material and contact owners. Equivalent radial
    /// distances are not physical planar lengths and explicitly refuse.
    pub fn with_current_planar_geometry<'r>(
        &'r self,
        thermal: &DepositionSimulation<'_>,
        current: CurrentConductiveLink,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<CurrentPlanarCaloricLink<'r, 'l, 'a, 'm>, DepositionError> {
        checkpoint()?;
        if self.geometry.radial_coordinates_m().is_some() {
            return Err(DepositionError::InvalidInput(
                "current planar caloric geometry cannot replace a radial source binding",
            ));
        }
        let (patch, energies) = admit_original_link(
            thermal,
            self.link_id(),
            self.layers,
            self.geometry,
            &mut checkpoint,
        )?;
        if patch != self.patch {
            return Err(DepositionError::InvalidInput(
                "current caloric geometry requires the original reference contact and coverage",
            ));
        }
        let geometry = LayerGeometry::planar(current.half_distances_m, current.area_m2)?;
        self.response(energies, geometry, [0.0; 2], &mut checkpoint)?;
        checkpoint()?;
        Ok(CurrentPlanarCaloricLink {
            reference: self,
            geometry,
        })
    }
}

impl CurrentPlanarCaloricLink<'_, '_, '_, '_> {
    /// The immutable reference geometry; it is not used as a current flux metric.
    #[must_use]
    pub fn reference_geometry(&self) -> LayerGeometry {
        self.reference.geometry()
    }

    /// The exact measures to supply to the native current-geometry map.
    #[must_use]
    pub fn current_geometry(&self) -> CurrentConductiveLink {
        CurrentConductiveLink {
            area_m2: self.geometry.contact_area_m2(),
            half_distances_m: self.geometry.equivalent_distances_m(),
        }
    }

    /// Geometry expected on a native query after its current measures refresh.
    /// Source-owned resistance still replaces a zero authored resistance once.
    #[must_use]
    pub fn native_geometry(&self) -> ConductiveLink {
        ConductiveLink {
            half_distances_m: self.geometry.equivalent_distances_m(),
            ..self.reference.native_geometry()
        }
    }

    /// Original full face response in this prepared mechanical configuration.
    /// Temperature derivatives hold the current contact resistance fixed;
    /// the native trait also includes the original contact's energy derivatives.
    pub fn evaluate(
        &self,
        energies_j_kg: [f64; 2],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<LayerConductionResponse, DepositionError> {
        self.reference
            .response(energies_j_kg, self.geometry, [0.0; 2], &mut checkpoint)
            .map(|(response, _)| response)
    }
}

impl ImplicitLinkLaw for CurrentPlanarCaloricLink<'_, '_, '_, '_> {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
        self.reference.caloric_laws()
    }

    fn supports_current_geometry(&self) -> bool {
        true
    }

    fn state_at_specific_enthalpies(
        &self,
        energies_j_kg: [f64; 2],
        geometry: ConductiveLink,
        patch: ContactPatch,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitLinkState, DepositionError> {
        checkpoint()?;
        // The mechanical measure replaces only the common flux area. Per-side
        // reference coverage and every original channel identity remain intact.
        let mut expected = self.reference.patch;
        expected.link.area_m2 = self.geometry.contact_area_m2();
        if patch != expected || geometry != self.native_geometry() {
            return Err(DepositionError::InvalidInput(
                "current caloric link requires its prepared measures and original reference coverage",
            ));
        }
        let (response, derivatives) =
            self.reference
                .response(energies_j_kg, self.geometry, [0.0; 2], checkpoint)?;
        Ok(ImplicitLinkState {
            conductance_w_k: response.conductance_w_k,
            power_enthalpy_derivatives_w_kg_j: derivatives,
        })
    }
}
