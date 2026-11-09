//! Finite material-contact voltage jumps on matching swept physical surfaces.
//!
//! Both ordered material sources and the specific resistance [ohm m2] remain
//! explicit. The same native bulk/contact solve retains exact radial P1
//! surface moments and separates bulk from interface Joule heat. Thermal
//! partition, contact motion and surface electrical storage are not inferred.

use super::{AxisymmetricConduction, AxisymmetricConductionField};
use crate::conduction::{
    ConductionContact, ConductionContactResponse, ConductionError, PlanarConduction,
    admit_contact_source, poll,
};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::conductor::planar::{PlanarConductorMaterial, PlanarOhmicResponse};
use fs_rep_mesh::TriComplex2;

impl<C: PlanarConductorMaterial> AxisymmetricConduction<C> {
    /// Prepare source-resolved bulk and matching electrical interfaces in one
    /// meridional field. R'' multiplies current density, so the physical swept
    /// area and radial surface basis determine terminal resistance and heat.
    /// Material, pressure, temperature or phase changes require preparation
    /// from freshly resolved original sources.
    ///
    /// Positive finite contact resistance is required. Perfect contact uses
    /// shared mesh vertices. The symmetry axis has no physical contact area.
    pub fn new_with_contacts(
        mesh: TriComplex2,
        conductors: Vec<C>,
        contacts: Vec<ConductionContact>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ConductionError> {
        poll(cancelled)?;
        if contacts.len() > limits.max_triplets / 16 {
            return Err(ConductionError::Invalid(
                "electrical contacts exceed the declared triplet capacity",
            ));
        }
        let axes = vec![[None, None]; contacts.len()];
        Self::new_with_phase_contacts(mesh, conductors, contacts, axes, limits, cancelled)
    }

    /// The same physical interfaces with explicitly owned liquid mass-fraction
    /// coordinates on either ordered side. Rows follow the input contacts;
    /// retained rows use the canonical order returned by `contacts()`.
    /// The phase owner resolves the actual fractions before every preparation.
    pub fn new_with_phase_contacts(
        mesh: TriComplex2,
        conductors: Vec<C>,
        contacts: Vec<ConductionContact>,
        liquid_mass_fraction_axes: Vec<[Option<String>; 2]>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ConductionError> {
        poll(cancelled)?;
        if !mesh.metric().is_axisymmetric() {
            return Err(ConductionError::Invalid(
                "axisymmetric conduction requires an explicit radial metric",
            ));
        }
        let core = PlanarConduction::assemble_metric_with_admitted_contacts(
            mesh,
            conductors,
            contacts,
            liquid_mass_fraction_axes,
            limits,
            cancelled,
            |source| source.has_admitted_static_sources(),
            |_, source, axes| admit_contact_source(source, axes),
        )?;
        Ok(Self { core })
    }
}

impl<C: PlanarOhmicResponse> AxisymmetricConduction<C> {
    /// Original ordered material interfaces, in canonical paired-edge order.
    pub fn contacts(&self) -> &[ConductionContact] {
        self.core.contacts()
    }

    /// Explicit phase-owner coordinates, in the same order as `contacts()`.
    pub fn contact_liquid_mass_fraction_axes(&self) -> &[[Option<String>; 2]] {
        self.core.contact_liquid_mass_fraction_axes()
    }
}

impl<C: PlanarOhmicResponse> AxisymmetricConductionField<'_, C> {
    /// Exact swept-surface current and irreversible heat. The order matches
    /// `problem().contacts()`. This heat is already included in `power_w()`
    /// and is separate from the volume heat in `cell_power_w()`.
    /// A thermal owner must declare its partition between the receiving bodies.
    pub fn contact_response(&self) -> &[ConductionContactResponse] {
        self.field.contact_response()
    }
}
