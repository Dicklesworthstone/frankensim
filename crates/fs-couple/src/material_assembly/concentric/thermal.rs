//! Source-bound layered thermal prestress, lowered to the existing realizer.
//!
//! Changing a layer's expansion path or current elastic properties changes its
//! stress and the string's resultant tension without a material-name branch.
//! The common material section owns all constitutive arithmetic; this adapter
//! owns only fixed supports and the existing taut-string admission boundary.

use fs_material::section::thermal::{ConcentricAxialResponse, ConcentricThermoelasticSection};
use fs_scenario::{AcousticAssembly, PrestressedString};

use crate::acoustic_realize::{AcousticRealizeError, RealizedAssembly};
use crate::string_specimen::StringPrestress;
use super::{invalid, positive, with_concentric_section};
use super::super::{AcousticMaterialBindings, CompiledMaterialAssembly, MaterialAssemblyError,
    compile_material_assembly, missing, physical};

/// Bind fixed supports to the shared per-layer thermal response.
///
/// Both endpoints are stationary. `reference_stress_free_length_m` is common
/// to all layers at their retained common reference temperature. Current length
/// and section geometry are explicit; this does not simulate heat diffusion,
/// Poisson contraction, helical winding or the evolution of material history.
/// No implicit retuning or inherited internal damping is permitted.
///
/// # Errors
/// Invalid reference length, moving supports, any inadmissible layer strain,
/// nonpositive net tension, or an ordinary elastic binding refusal. The generic
/// section admits compression; the taut-string realizer deliberately does not.
pub fn with_concentric_thermal_extension(
    input: PrestressedString,
    section: &ConcentricThermoelasticSection,
    reference_stress_free_length_m: f64,
    linear_strain_limit: f64,
) -> Result<(PrestressedString, ConcentricAxialResponse), AcousticRealizeError> {
    if input.moving_end {
        return Err(invalid("concentric thermal extension requires stationary endpoints"));
    }
    let length = positive(input.length_m, "thermal section member length must be positive and finite")?;
    let reference = positive(reference_stress_free_length_m,
        "thermal section reference length must be positive and finite")?;
    let strain = (length - reference) / reference;
    let response = section.response_at_axial_strain(strain, linear_strain_limit)
        .map_err(|_| invalid("concentric thermal response exceeds its declared strain or arithmetic domain"))?;
    let bound = with_concentric_section(input, section.section(),
        StringPrestress::FixedTension(response.force_n()))?;
    Ok((bound, response))
}

/// One existing compiled assembly with all layered thermal source data retained.
#[derive(Clone, Debug)]
pub struct CompiledThermoelasticConcentricAssembly {
    compiled: CompiledMaterialAssembly,
    section: ConcentricThermoelasticSection,
    response: ConcentricAxialResponse,
    reference_stress_free_length_m: f64,
    linear_strain_limit: f64,
}
impl CompiledThermoelasticConcentricAssembly {
    /// Original compiled assembly, including other material bindings.
    #[must_use]
    pub const fn compiled(&self) -> &CompiledMaterialAssembly { &self.compiled }
    /// Exact source-bound thermal and elastic section, not an effective metal.
    #[must_use]
    pub const fn section(&self) -> &ConcentricThermoelasticSection { &self.section }
    /// Per-layer stresses and the resultant passed to the realizer.
    #[must_use]
    pub const fn axial_response(&self) -> &ConcentricAxialResponse { &self.response }
    /// Authored common stress-free reference length [m].
    #[must_use]
    pub const fn reference_stress_free_length_m(&self) -> f64 { self.reference_stress_free_length_m }
    /// Authored small-strain applicability bound; not a material yield limit.
    #[must_use]
    pub const fn linear_strain_limit(&self) -> f64 { self.linear_strain_limit }
    /// Run the ordinary coupled dynamics, medium loading and pressure observer.
    pub fn realize(&self) -> Result<RealizedAssembly, AcousticRealizeError> { self.compiled.realize() }
}

/// Compile a thermally prestressed layered member with the ordinary material
/// assembly compiler. Plate/body and medium behavior stay with their existing
/// owners; this function neither creates a second solver nor retunes the string.
///
/// # Errors
/// Missing string, competing homogeneous string binding, invalid supports or
/// thermal strain, or an existing compiler refusal. Inputs remain unchanged.
pub fn compile_concentric_thermal_material_assembly(
    template: &AcousticAssembly,
    section: ConcentricThermoelasticSection,
    reference_stress_free_length_m: f64,
    linear_strain_limit: f64,
    other_bindings: &AcousticMaterialBindings<'_>,
) -> Result<CompiledThermoelasticConcentricAssembly, MaterialAssemblyError> {
    if other_bindings.string.is_some() {
        return Err(missing("string", "concentric and homogeneous string bindings are mutually exclusive"));
    }
    let input = template.string.as_ref()
        .ok_or_else(|| missing("string", "concentric thermal binding requires a string description"))?;
    let (string, response) = with_concentric_thermal_extension(input.clone(), &section,
        reference_stress_free_length_m, linear_strain_limit)
        .map_err(|error| physical("string", error))?;
    let mut candidate = template.clone();
    candidate.string = Some(string);
    let compiled = compile_material_assembly(&candidate, other_bindings)?;
    Ok(CompiledThermoelasticConcentricAssembly { compiled, section, response,
        reference_stress_free_length_m, linear_strain_limit })
}

/// Acoustic observations of accepted, thermally constrained material histories.
pub mod heating;
