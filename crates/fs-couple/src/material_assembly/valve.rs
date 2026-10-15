//! Material/geometry-derived pressure valves through the ordinary assembly.
//!
//! A reed's mesh, oriented material field, pressure load and reduction controls
//! determine its mechanics. This adapter does not choose an instrument model or
//! copy a nominal cane stiffness. It admits the actual first-mode image, uses
//! its work-conjugate projected area, then composes the EXISTING material-bound
//! string/plate compiler, duct solver, contact law and acoustic medium.
//!
//! The supported playing image remains one inward-striking degree of freedom.
//! Two independent oboe leaflets, rolling lay contact, material-dependent
//! damping spectra, large strain, and transient moisture evolution require
//! their own admitted models. A successful compilation does not claim them.

use fs_material::gas::GasState;
use fs_scenario::AcousticAssembly;
use fs_scenario::reed_source::BeatingReedSource;

use super::{AcousticMaterialBindings, CompiledMaterialAssembly, MaterialAssemblyError,
    compile_material_assembly, missing};
use crate::acoustic_realize::{AcousticRealizeError, RealizedAssembly};

/// A component-local compiler refusal, retaining the responsible owner's error.
#[derive(Debug)]
pub enum ValveMaterialAssemblyError<E> {
    /// The exact playing image failed its reduction/conversion admission.
    Reduction(E),
    /// Missing assembly members or ordinary solid material binding refusal.
    Assembly(MaterialAssemblyError),
}
impl<E: core::fmt::Display> core::fmt::Display for ValveMaterialAssemblyError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Reduction(e) => write!(f, "FS-COUPLE-VALVE-REDUCTION: {e}"),
            Self::Assembly(e) => write!(f, "{e}"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for ValveMaterialAssemblyError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self { Self::Reduction(e) => Some(e), Self::Assembly(e) => Some(e) }
    }
}
impl<E> From<MaterialAssemblyError> for ValveMaterialAssemblyError<E> {
    fn from(error: MaterialAssemblyError) -> Self { Self::Assembly(error) }
}

/// Immutable numeric assembly retaining its actual valve and adequacy criteria.
/// Source material/mesh identities and pressure projections remain accessible;
/// publishing solver coefficients does not turn authored loss into measurement.
#[derive(Clone, Debug)]
pub struct CompiledValveMaterialAssembly<S: BeatingReedSource> {
    compiled: CompiledMaterialAssembly,
    valve: S,
    admission: S::Admission,
}
impl<S: BeatingReedSource> CompiledValveMaterialAssembly<S> {
    /// The ordinary compiled assembly, including any source-bound string/body.
    #[must_use]
    pub const fn compiled(&self) -> &CompiledMaterialAssembly { &self.compiled }
    /// Original structural reduction, pressure port and exact static responses.
    #[must_use]
    pub const fn valve(&self) -> &S { &self.valve }
    /// Explicit tolerances used for the first-mode image, not a hidden preset.
    #[must_use]
    pub const fn admission(&self) -> &S::Admission { &self.admission }

    /// Render through the ordinary compiler's ambient-gas path.
    ///
    /// # Errors
    /// Existing duct, body, contact, bandwidth, excitation and runtime refusals.
    pub fn realize(&self) -> Result<RealizedAssembly, AcousticRealizeError> {
        self.compiled.realize()
    }

    /// Render the same frozen solid materials in the explicitly supplied gas.
    ///
    /// This uses the common gas-dependent excitation, bore, radiation and
    /// propagation path. It does not silently equilibrate the solid's material
    /// temperature/moisture with the ambient or retune it to a nominal pitch.
    /// A new solid state must be resolved and reduced explicitly.
    ///
    /// # Errors
    /// Common supplied-gas coherence and ordinary realization refusals propagate.
    pub fn realize_with_gas(&self, gas: &GasState) -> Result<RealizedAssembly, AcousticRealizeError> {
        self.compiled.realize_with_gas(gas)
    }
}

/// Compile a pressure-loaded material valve with the existing solid bindings.
///
/// The template supplies a duct and a reed's playing controls (mouth pressure
/// and attack). Rest opening, width, mass, stiffness, damping and closing
/// pressure are ALL replaced by the admitted reduction. The effective area is
/// the actual modal pressure projection, not width times a nominal face length.
/// Other body/solid material bindings still use their original shared compiler.
/// The template and any previously compiled assembly are unchanged on refusal.
///
/// # Errors
/// A missing reed/duct, inadequate first-mode image, nonphysical controls, or
/// an ordinary string/plate material refusal aborts before publication.
pub fn compile_pressure_loaded_reed_assembly<S: BeatingReedSource>(
    template: &AcousticAssembly,
    valve: S,
    admission: S::Admission,
    other_bindings: &AcousticMaterialBindings<'_>,
) -> Result<CompiledValveMaterialAssembly<S>, ValveMaterialAssemblyError<S::Error>> {
    let control = template.reed.ok_or_else(|| missing(
        "reed", "material valve binding requires a reed description with playing controls",
    ))?;
    if template.duct.is_none() {
        return Err(missing("reed", "material valve binding requires an existing duct").into());
    }
    let reed = valve.compile_beating_reed(
        control.blowing_pressure_pa, control.attack_s, &admission,
    ).map_err(ValveMaterialAssemblyError::Reduction)?;
    let mut candidate = template.clone();
    candidate.reed = Some(reed);
    let compiled = compile_material_assembly(&candidate, other_bindings)?;
    Ok(CompiledValveMaterialAssembly { compiled, valve, admission })
}
