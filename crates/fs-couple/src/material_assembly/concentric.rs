//! Bind source-resolved concentric sections into the existing assembly solver.
//!
//! This is an elastic, perfectly bonded slender member, NOT a helical wrapped
//! string. Geometry, material assignment and support prescription are data.
//! Every downstream mass, axial nonlinear stiffness, bending dispersion and
//! radiating width comes from the same section; no solver is reimplemented.

use fs_material::section::ConcentricElasticSection;
use fs_scenario::{AcousticAssembly, PrestressedString};
use crate::acoustic_realize::{AcousticRealizeError, RealizedAssembly};
use crate::string_specimen::StringPrestress;
use super::{AcousticMaterialBindings, CompiledMaterialAssembly, MaterialAssemblyError,
    compile_material_assembly, missing, physical};

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}
fn positive(value: f64, what: &'static str) -> Result<f64, AcousticRealizeError> {
    if value.is_finite() && value > 0.0 { Ok(value) } else { Err(invalid(what)) }
}

/// Replace all derived elastic coefficients in one transaction.
///
/// The section is stress-free in a common reference configuration. Differential
/// thermal strain requires a layered thermal model, not one invented average.
/// Existing nonzero/internal material damping is refused: elastic material
/// data cannot justify copying another material's dissipation spectrum. Air
/// loading remains the existing realizer's separate medium-dependent model.
///
/// # Errors
/// Refuses invalid supports, slack/compression, incompatible thermal or damping
/// laws, and target-frequency tuning of moving-end modes. Remaining modal and
/// excitation checks stay in the ordinary acoustic realizer.
pub fn with_concentric_section(
    mut input: PrestressedString,
    section: &ConcentricElasticSection,
    prestress: StringPrestress,
) -> Result<PrestressedString, AcousticRealizeError> {
    let length = positive(input.length_m, "section member length must be finite and positive")?;
    if input.damping_ratio != 0.0 || input.rayleigh.is_some()
        || input.kelvin_voigt_bending.is_some() || input.relaxation_bending.is_some()
    {
        return Err(invalid("elastic section binding requires explicitly zero internal loss and no inherited loss law"));
    }
    let axial = section.axial_stiffness_n();
    let tension = match prestress {
        StringPrestress::FixedTension(value) => value,
        StringPrestress::FixedExtension { stress_free_length_m, linear_strain_limit } => {
            let rest = positive(stress_free_length_m, "section stress-free length must be finite and positive")?;
            if !linear_strain_limit.is_finite() || linear_strain_limit <= 0.0 || linear_strain_limit >= 1.0 {
                return Err(invalid("section linear strain limit must lie in (0, 1)"));
            }
            let strain = (length - rest) / rest;
            if !strain.is_finite() || strain <= 0.0 || strain > linear_strain_limit {
                return Err(invalid("section extension is slack, compressed or outside its declared linear strain limit"));
            }
            axial * strain
        }
        StringPrestress::TargetFundamentalHz(frequency) => {
            if input.moving_end {
                return Err(invalid("section target-frequency constraint requires fixed-fixed modes"));
            }
            positive(frequency, "section target frequency must be finite and positive")?;
            // omega_1^2 = (T k^2 + EI k^4) / mu, k = pi/L.
            let k = core::f64::consts::PI / length;
            let speed = 2.0 * length * frequency;
            section.line_density_kg_m() * speed * speed - section.bending_stiffness_n_m2() * k * k
        }
        StringPrestress::FixedThermalExtension { .. } => {
            return Err(invalid("composite thermal prestress requires per-layer free strains, not a shared scalar override"));
        }
    };
    input.tension_n = positive(tension, "section support prescription must produce finite positive tension")?;
    input.lin_density_kg_m = section.line_density_kg_m();
    input.axial_stiffness_n = axial;
    input.bending_stiffness_n_m2 = section.bending_stiffness_n_m2();
    input.width_m = positive(2.0 * section.outer_radius_m(), "section radiating width must be finite and positive")?;
    Ok(input)
}

/// An ordinary material-compiled assembly with its layered section retained.
#[derive(Clone, Debug)]
pub struct CompiledConcentricAssembly {
    compiled: CompiledMaterialAssembly,
    section: ConcentricElasticSection,
    prestress: StringPrestress,
}
impl CompiledConcentricAssembly {
    /// Existing compiled assembly, including any uniform/regional plate binding.
    #[must_use]
    pub const fn compiled(&self) -> &CompiledMaterialAssembly { &self.compiled }
    /// Original source layers and section integrals.
    #[must_use]
    pub const fn section(&self) -> &ConcentricElasticSection { &self.section }
    /// Mechanical prescription used when compiling this immutable specimen.
    #[must_use]
    pub const fn prestress(&self) -> StringPrestress { self.prestress }
    /// Run the same acoustic solver and observer as ordinary material assemblies.
    pub fn realize(&self) -> Result<RealizedAssembly, AcousticRealizeError> { self.compiled.realize() }
}

/// Compile a layered string together with the existing plate/body bindings.
/// A simultaneous homogeneous string binding refuses instead of overwriting
/// one set of material coefficients with another. The input remains unchanged.
///
/// # Errors
/// Propagates section/support, conflicting binding and ordinary material
/// compilation refusals. No candidate is published on failure.
pub fn compile_concentric_material_assembly(
    template: &AcousticAssembly,
    section: ConcentricElasticSection,
    prestress: StringPrestress,
    other_bindings: &AcousticMaterialBindings<'_>,
) -> Result<CompiledConcentricAssembly, MaterialAssemblyError> {
    if other_bindings.string.is_some() {
        return Err(missing("string", "concentric and homogeneous string bindings are mutually exclusive"));
    }
    let input = template.string.as_ref()
        .ok_or_else(|| missing("string", "concentric material binding requires a string description"))?;
    let string = with_concentric_section(input.clone(), &section, prestress)
        .map_err(|error| physical("string", error))?;
    let mut candidate = template.clone();
    candidate.string = Some(string);
    let compiled = compile_material_assembly(&candidate, other_bindings)?;
    Ok(CompiledConcentricAssembly { compiled, section, prestress })
}
