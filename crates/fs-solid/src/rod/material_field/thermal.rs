//! Prescribed spatial thermal states coupled to the existing rod equilibrium.
//!
//! Every segment can have different current temperature, material and layer
//! expansion. The reference length describes the initially unloaded compatible
//! configuration at each retained expansion path's reference state. Current
//! moduli, density and integrated expansion must come from exactly the same
//! layer bundles. No endpoint-alpha approximation or whole-body scaling is used.
//!
//! This is one-way quasistatic coupling to a supplied thermal field, not a heat
//! solver, transient thermodynamic balance or remap of a moving material. Section
//! geometry is frozen. The source's bound controls total and per-layer axial
//! strains only; it is not a 3D yield criterion, bending-fiber strain certificate,
//! or general proof that a large-rotation equilibrium is physically admissible.

use fs_exec::CancelGate;
use fs_material::section::thermal::{ConcentricAxialResponse, ConcentricThermoelasticSection};
use super::{SpatialMaterialRod, SpatialRodError, checkpoint};

/// Temperature-only refresh of all section properties from pinned material cards.
pub mod source;

#[derive(Clone, Debug)]
pub(super) struct ThermalField {
    sections: Vec<ConcentricThermoelasticSection>,
    limit: f64,
}
impl ThermalField {
    pub(super) fn summary(&self, segment: usize, strain: f64) -> Result<(f64, f64), SpatialRodError> {
        Ok(self.sections[segment].axial_force_and_energy_at_strain(strain, self.limit)?)
    }
}

impl SpatialMaterialRod {
    /// Return a new static specimen with a complete source-bound thermal field.
    /// The old pose is an initial guess, not a transported dynamic/thermal history.
    /// Each thermal section must retain EXACTLY the elastic section assigned to
    /// that segment, including source identities, query points and layer geometry.
    /// `linear_strain_limit` must lie in (0,1); every layer is checked even when
    /// opposing thermal stresses make the net section force small.
    ///
    /// Native forces, tangents and support reactions are analytical at the
    /// actual admitted strain; no finite-difference probe needs extra domain
    /// margin. Every physical trial still obeys every original layer bound.
    ///
    /// # Errors
    /// Refuses count/source/state mismatches, invalid limits or strains, and
    /// cancellation. The previous specimen and its accepted pose remain intact.
    pub fn with_thermal_sections(
        &self, sections: Vec<ConcentricThermoelasticSection>,
        linear_strain_limit: f64, gate: &CancelGate,
    ) -> Result<Self, SpatialRodError> {
        checkpoint(gate)?;
        if sections.len() != self.sections.len() {
            return Err(SpatialRodError::Invalid("exactly one thermal section is required per rod segment"));
        }
        for (thermal, mechanical) in sections.iter().zip(&self.sections) {
            checkpoint(gate)?;
            if thermal.section() != mechanical.elastic() {
                return Err(SpatialRodError::Invalid("thermal and mechanical segment must share the exact elastic source bundle and geometry"));
            }
        }
        let field = ThermalField { sections, limit: linear_strain_limit };
        let mut candidate = self.clone();
        candidate.thermal = Some(field);
        candidate.energy_j()?;
        checkpoint(gate)?;
        Ok(candidate)
    }

    /// Retained thermal sections in spatial order, or None for the elastic path.
    #[must_use]
    pub fn thermal_sections(&self) -> Option<&[ConcentricThermoelasticSection]> {
        self.thermal.as_ref().map(|field| field.sections.as_slice())
    }

    /// Centerline axial stress/force in EVERY layer at the accepted nodal state.
    /// A free section can have zero resultant but nonzero opposing layer forces
    /// and stored energy. Bending stresses across a layer are not reported here.
    /// None means that no thermal field was bound, not zero thermal stress.
    pub fn thermal_axial_response(
        &self, segment: usize,
    ) -> Result<Option<ConcentricAxialResponse>, SpatialRodError> {
        let (g, _) = self.segment_strains(segment)?;
        self.thermal.as_ref().map(|field| {
            field.sections[segment].response_at_axial_strain(g[0], field.limit)
                .map_err(SpatialRodError::from)
        }).transpose()
    }
}
