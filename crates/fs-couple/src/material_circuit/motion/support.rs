//! A layered axial member loading the original reciprocal field-motion owner.
//!
//! The material section owns force, elastic energy and layer strain admission.
//! Here the authored displacement field u(X) = X * (gap-reference_gap) / L
//! maps that law to the guided coordinate and gives modal mass rho A L / 3.
//! This is a one-coordinate axial reduction with a fixed far end, not a resolved
//! elastic wave, bending/buckling, or temperature-evolution model.

use fs_material::section::thermal::ConcentricThermoelasticSection;
use fs_material::state_point::{DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY};

use super::{EvaluationDecision, MaterialCircuitError, invalid, positive};

/// Immutable source-bound return member, including any frozen thermal preload.
/// Positive gap displacement extends the member. Its mass is distinct from the
/// guided body's regional mass; neither is supplied as an override coefficient.
#[derive(Clone, Debug, PartialEq)]
pub struct AxialSupport {
    section: ConcentricThermoelasticSection,
    length_m: f64,
    reference_gap_m: f64,
    linear_strain_limit: f64,
    modal_mass_kg: f64,
}

impl AxialSupport {
    /// Retain the original section and its complete integrated expansion paths.
    /// The reference gap has zero TOTAL axial strain, not necessarily zero force:
    /// the original per-layer thermal eigenstrains determine any preload.
    ///
    /// E and reference density must be validity-wide constant source claims.
    /// Temperatures and eigenstrains remain frozen during motion; no thermal
    /// modulus work, changing mass, dashpot heat or numerical damping is inferred.
    /// The explicit strain limit is the original section's model-applicability
    /// assumption, not a yield or buckling certificate.
    pub fn try_new(section: ConcentricThermoelasticSection, length_m: f64,
        reference_gap_m: f64, linear_strain_limit: f64) -> Result<Self, MaterialCircuitError>
    {
        if !positive(length_m) || !positive(reference_gap_m) {
            return Err(invalid("axial support requires positive finite length and reference gap"));
        }
        for layer in section.section().layers() {
            for name in [DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY] {
                let property = layer.material.property(name)
                    .ok_or_else(|| invalid("axial support layer is missing its original constitutive property"))?;
                if property.answer().receipt.decision != EvaluationDecision::ConstantWithinValidity {
                    return Err(invalid("conservative axial support requires validity-wide constant modulus and reference density"));
                }
            }
        }
        let modal_mass_kg = section.section().line_density_kg_m() * length_m / 3.0;
        let stiffness = section.section().axial_stiffness_n() / length_m;
        if !positive(modal_mass_kg) || !positive(stiffness) {
            return Err(invalid("source-derived support inertia or stiffness is unrepresentable"));
        }
        // Validate the material's strain policy and all retained eigenstrains.
        // Initial actuator geometry is checked separately by its constructor.
        section.axial_force_and_energy_at_strain(0.0, linear_strain_limit)?;
        Ok(Self { section, length_m, reference_gap_m, linear_strain_limit, modal_mass_kg })
    }

    /// Original layers, material states and integrated thermal expansion paths.
    pub const fn section(&self) -> &ConcentricThermoelasticSection { &self.section }
    /// Fixed reference member length [m].
    pub const fn length_m(&self) -> f64 { self.length_m }
    /// Gap at zero total member strain [m].
    pub const fn reference_gap_m(&self) -> f64 { self.reference_gap_m }
    /// Per-layer total, thermal and elastic strain applicability bound.
    pub const fn linear_strain_limit(&self) -> f64 { self.linear_strain_limit }
    /// Consistent inertia of the linear axial displacement shape [kg].
    pub const fn modal_mass_kg(&self) -> f64 { self.modal_mass_kg }
    /// Original section EA divided by authored reference length [N/m].
    pub fn stiffness_n_m(&self) -> f64 { self.section.section().axial_stiffness_n() / self.length_m }

    /// Signed restoring force on the moving coordinate and elastic energy [N,J].
    /// Every call delegates to the original per-layer constitutive owner, which
    /// rejects excess thermal/elastic strain and unrepresentable layer stresses.
    pub fn force_and_energy(&self, gap_m: f64) -> Result<(f64, f64), MaterialCircuitError> {
        if !positive(gap_m) {
            return Err(invalid("axial support gap must be finite and positive"));
        }
        let displacement = gap_m - self.reference_gap_m;
        let strain = displacement / self.length_m;
        if !strain.is_finite() || (displacement != 0.0 && strain == 0.0) {
            return Err(invalid("axial support strain is unrepresentable"));
        }
        let (force, per_length) = self.section.axial_force_and_energy_at_strain(
            strain, self.linear_strain_limit)?;
        let energy = per_length * self.length_m;
        if !energy.is_finite() || (per_length != 0.0 && energy == 0.0) {
            return Err(invalid("axial support elastic energy is unrepresentable"));
        }
        Ok((-force, energy))
    }
}
