//! Accepted thermal history -> physical fixed-end prestress -> existing sound.
use fs_material::section::thermal::AxialConstraint;
use fs_material::section_temperature::heating::HeatedConcentricSection;
use fs_scenario::AcousticAssembly;
use crate::acoustic_realize::{AcousticRealizeError, RealizedAssembly};
use crate::material_assembly::{AcousticMaterialBindings, MaterialAssemblyError, missing};
use crate::material_assembly::concentric::{CompiledConcentricAssembly, compile_concentric_material_assembly};
use crate::string_specimen::StringPrestress;

/// An ordinary acoustic compilation tied to one accepted heated specimen.
/// Borrowing the source prevents a different thermal epoch being substituted.
#[derive(Clone)]
pub struct CompiledHeatedConcentricAssembly<'s, 'm> {
    source: &'s HeatedConcentricSection<'m>,
    compiled: CompiledConcentricAssembly,
}
impl<'s, 'm> CompiledHeatedConcentricAssembly<'s, 'm> {
    /// Original thermal history, layer stresses and source receipts.
    pub const fn source(&self) -> &'s HeatedConcentricSection<'m> { self.source }
    /// Existing compiled assembly, with mass/stiffness and physical tension
    /// taken from the same accepted state, including other body/plate bindings.
    pub const fn compiled(&self) -> &CompiledConcentricAssembly { &self.compiled }
    /// Run the existing realizer and physical-pressure observer unchanged.
    pub fn realize(&self) -> Result<RealizedAssembly, AcousticRealizeError> { self.compiled.realize() }
}

/// Observe a heated fixed-end member without retensioning it to its old force
/// or tuning it to its old pitch. The accepted axial response supplies the
/// resultant; the same section supplies line mass, bending and axial stiffness.
///
/// The supplied string length must equal the fixed thermal specimen length.
/// Force/spring-controlled supports are deliberately refused on this observer:
/// their equilibrium strain does not move the existing thermal geometry, and
/// they must not be passed off as stationary fixed-end strings. Compression and
/// slack remain valid material outcomes, but refuse at the taut-string boundary.
///
/// This is a frozen-state observation, not audio during a changing-temperature
/// trajectory, coupled acoustic heat/work exchange or a helical winding model.
/// Inherited material damping and duplicate homogeneous bindings remain errors
/// under the original assembly compiler; no fallback material is selected.
pub fn compile_heated_concentric_material_assembly<'s, 'm>(
    template: &AcousticAssembly,
    source: &'s HeatedConcentricSection<'m>,
    other_bindings: &AcousticMaterialBindings<'_>,
) -> Result<CompiledHeatedConcentricAssembly<'s, 'm>, MaterialAssemblyError> {
    let state = source.axial_state()
        .ok_or_else(|| missing("string", "heated acoustic source needs an accepted axial constraint"))?;
    let string = template.string.as_ref()
        .ok_or_else(|| missing("string", "heated acoustic observation needs string geometry"))?;
    if !matches!(state.constraint(), AxialConstraint::FixedStrain(_)) || string.moving_end
        || string.length_m.to_bits() != source.length_m().to_bits()
    {
        return Err(missing("string", "heated sound requires stationary fixed-strain supports and the same physical length"));
    }
    // This is only the realizer's frozen force input. It is NOT a new support
    // prescription: the borrowed thermal source retains the actual fixed strain.
    let compiled = compile_concentric_material_assembly(template, source.section().clone(),
        StringPrestress::FixedTension(state.response().force_n()), other_bindings)?;
    Ok(CompiledHeatedConcentricAssembly { source, compiled })
}
