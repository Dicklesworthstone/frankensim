//! Quasistatic axial response observed at every accepted thermal substep.
use super::{HeatedConcentricSection, SectionHeatingError};
use crate::section::{ConcentricElasticSection, thermal::{AxialConstraint,
    ConcentricAxialResponse, ConcentricThermoelasticSection}};
use crate::section_temperature::ConcentricSectionTemperatureLaw;

/// One accepted thermal/mechanical source state and its independent support data.
/// No public field permits substituting stress from another material epoch.
#[derive(Clone, Debug, PartialEq)]
pub struct SectionAxialState {
    constraint: AxialConstraint,
    linear_strain_limit: f64,
    section: ConcentricThermoelasticSection,
    response: ConcentricAxialResponse,
}
impl SectionAxialState {
    pub(super) fn resolve(
        law: &ConcentricSectionTemperatureLaw<'_>, section: ConcentricElasticSection,
        constraint: AxialConstraint, linear_strain_limit: f64,
    ) -> Result<Self, SectionHeatingError> {
        let section = law.bind_expansion_paths(section)?;
        let response = section.response_with_axial_constraint(constraint, linear_strain_limit)?;
        Ok(Self { constraint, linear_strain_limit, section, response })
    }
    /// Original support prescription, never a target-pitch retuning operation.
    pub const fn constraint(&self) -> AxialConstraint { self.constraint }
    /// Explicit layer/total strain bound, not a measured yield envelope.
    pub const fn linear_strain_limit(&self) -> f64 { self.linear_strain_limit }
    /// Same material bundles as the accepted thermal endpoint, with integrated
    /// expansion paths from the original common stress-free temperature.
    pub const fn section(&self) -> &ConcentricThermoelasticSection { &self.section }
    /// Per-layer stress, support resultant and diagnostic elastic energy.
    pub const fn response(&self) -> &ConcentricAxialResponse { &self.response }
}

impl HeatedConcentricSection<'_> {
    /// Attach or replace a quasistatic support prescription at the CURRENT
    /// thermal state, retaining the ORIGINAL stress-free material reference.
    /// The resulting response is refreshed and strain-checked at every actual
    /// thermal substep, not merely at the requested output time. A failed step
    /// publishes neither the thermal nor the mechanical candidate.
    ///
    /// This is deliberately one-way thermal-to-mechanical coupling: heat rates,
    /// reference mass and transport geometry are unchanged. Returned elastic
    /// energy is diagnostic, never added to/debited from the heat audit. Force
    /// and spring control report equilibrium strain but do not move the thermal
    /// mesh. No thermodynamic feedback, buckling or finite-deformation claim.
    /// Missing expansion is an error; ordinary unbound heating remains available.
    pub fn with_axial_constraint(
        mut self, constraint: AxialConstraint, linear_strain_limit: f64,
    ) -> Result<Self, SectionHeatingError> {
        self.observer.axial = Some(SectionAxialState::resolve(
            &self.observer.law, self.observer.section.clone(), constraint, linear_strain_limit,
        )?);
        Ok(self)
    }

    /// Accepted layer stresses and support force, absent only when the caller
    /// did not request thermoelastic response. No stale or partial state escapes.
    pub const fn axial_state(&self) -> Option<&SectionAxialState> { self.observer.axial.as_ref() }

    /// Fixed physical length used by the thermal and electrical geometry [m].
    /// It is not recomputed from a quasistatic force-controlled strain response.
    pub const fn length_m(&self) -> f64 { self.observer.length_m }
}
