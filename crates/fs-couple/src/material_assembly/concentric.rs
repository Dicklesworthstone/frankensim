//! Bind source-resolved concentric sections into the existing assembly solver.
//!
//! This is a perfectly bonded slender member, NOT a helical wrapped string.
//! Geometry, material assignment and support prescription are data. Every
//! downstream mass, axial nonlinear stiffness, bending dispersion and radiating
//! width comes from the same section; no solver is reimplemented. Optional
//! source-bound Kelvin-Voigt bending replaces, rather than inherits, internal
//! loss. The [`thermal`] adapter retains each layer's thermal path and stress.

pub mod thermal;
pub use thermal::{CompiledThermoelasticConcentricAssembly,
    compile_concentric_thermal_material_assembly, with_concentric_thermal_extension};

use fs_material::section::{ConcentricElasticSection, ConcentricKelvinVoigtSection};
use fs_material::section::relaxation::ConcentricPronySection;
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
    input: PrestressedString,
    section: &ConcentricElasticSection,
    prestress: StringPrestress,
) -> Result<PrestressedString, AcousticRealizeError> {
    with_section_coefficients(input, section.line_density_kg_m(), section.axial_stiffness_n(),
        section.bending_stiffness_n_m2(), section.outer_radius_m(), prestress)
}

// One mechanical support law for every source-bound circular section. A new
// constitutive section supplies integrals, never a second tuning prescription.
pub(super) fn with_section_coefficients(
    mut input: PrestressedString, mass_per_length: f64, axial: f64,
    bending: f64, outer_radius: f64, prestress: StringPrestress,
) -> Result<PrestressedString, AcousticRealizeError> {
    let length = positive(input.length_m, "section member length must be finite and positive")?;
    if input.damping_ratio != 0.0 || input.rayleigh.is_some()
        || input.kelvin_voigt_bending.is_some() || input.relaxation_bending.is_some()
        || input.bending_loss_factor.is_some()
    {
        return Err(invalid("elastic section binding requires explicitly zero internal loss and no inherited loss law"));
    }
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
            mass_per_length * speed * speed - bending * k * k
        }
        StringPrestress::FixedThermalExtension { .. } => {
            return Err(invalid("composite thermal prestress requires per-layer free strains; use with_concentric_thermal_extension"));
        }
    };
    input.tension_n = positive(tension, "section support prescription must produce finite positive tension")?;
    input.lin_density_kg_m = mass_per_length;
    input.axial_stiffness_n = axial;
    input.bending_stiffness_n_m2 = bending;
    input.width_m = positive(2.0 * outer_radius, "section radiating width must be finite and positive")?;
    Ok(input)
}

/// Bind the original per-layer bending spectra and the same section/support law.
///
/// Every arm keeps its own integrated stiffness and relaxation time. No
/// Kelvin-Voigt coefficient, equivalent material, or averaged time is inferred.
/// Consumers must retain the material histories or refuse this description.
/// Original layer identities and selectors remain in `section`.
///
/// # Errors
/// The original elastic section/support admission rejects inherited internal
/// loss, unsupported geometry or nonphysical support prescriptions.
pub fn with_concentric_prony_section(
    input: PrestressedString,
    section: &ConcentricPronySection,
    prestress: StringPrestress,
) -> Result<PrestressedString, AcousticRealizeError> {
    let mut string = with_concentric_section(input, section.elastic(), prestress)?;
    string.relaxation_bending = Some(concentric_prony_bending(section));
    Ok(string)
}

/// Shared mechanical lowering of an already admitted source spectrum.
pub(crate) fn concentric_prony_bending(
    section: &ConcentricPronySection,
) -> fs_scenario::acoustic::PronyBending {
    fs_scenario::acoustic::PronyBending {
        branches: section.branches().iter().map(|arm|
            fs_scenario::acoustic::BendingRelaxationBranch {
                relaxing_stiffness_n_m2: arm.relaxing_stiffness_n_m2,
                relaxation_time_s: arm.relaxation_time_s,
            }).collect(),
        omega_band_rad_s: section.omega_band_rad_s(),
        // There is no single equivalent material: the owning section retains
        // every layer's actual source state, spatial assignment and selectors.
        material_state_identity: None,
        source_properties: None,
    }
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

/// Bind all elastic coefficients AND per-layer Kelvin-Voigt bending loss.
/// The ordinary realizer consumes the summed viscous stiffness and checks its
/// intersected source band against every retained reference-mode frequency.
/// Air drag and radiation continue to use the assembly's independent medium.
///
/// # Errors
/// Inherited internal/Rayleigh/Prony loss and invalid support prescriptions
/// refuse through the same elastic binding. No alternative acoustic solver or
/// nominal-material damping preset is selected.
pub fn with_concentric_kelvin_voigt_section(
    input: PrestressedString,
    section: &ConcentricKelvinVoigtSection,
    prestress: StringPrestress,
) -> Result<PrestressedString, AcousticRealizeError> {
    let mut string = with_concentric_section(input, section.elastic(), prestress)?;
    string.kelvin_voigt_bending = Some(fs_scenario::acoustic::KelvinVoigtBending {
        viscous_stiffness_n_m2_s: section.viscous_stiffness_n_m2_s(),
        omega_band_rad_s: section.omega_band_rad_s(),
        // This is NOT one homogeneous state. The compiled wrapper below retains
        // every source bundle rather than mislabeling an aggregate as a source ID.
        material_state_identity: None,
    });
    Ok(string)
}

/// Actual material-compiled acoustic assembly and its source-bound layered loss.
#[derive(Clone, Debug)]
pub struct CompiledConcentricKelvinVoigtAssembly {
    compiled: CompiledMaterialAssembly,
    section: ConcentricKelvinVoigtSection,
    prestress: StringPrestress,
}
impl CompiledConcentricKelvinVoigtAssembly {
    /// Existing assembly compiler result, including other material bindings.
    #[must_use]
    pub const fn compiled(&self) -> &CompiledMaterialAssembly { &self.compiled }
    /// All elastic and viscous source data, geometry and validity bounds.
    #[must_use]
    pub const fn section(&self) -> &ConcentricKelvinVoigtSection { &self.section }
    /// Independent mechanical support prescription.
    #[must_use]
    pub const fn prestress(&self) -> StringPrestress { self.prestress }
    /// Execute the existing pressure solver, not a separately synthesized decay.
    pub fn realize(&self) -> Result<RealizedAssembly, AcousticRealizeError> { self.compiled.realize() }
}

/// Compile a material-swappable, lossy layered member into the ordinary coupled
/// assembly. Use the existing string-specimen bending-viscosity property name
/// when constructing the section. The template is unchanged on success/failure.
///
/// # Errors
/// Conflicting homogeneous bindings, inherited loss, missing string geometry,
/// invalid support prescriptions and ordinary assembly admission failures refuse.
pub fn compile_concentric_kelvin_voigt_material_assembly(
    template: &AcousticAssembly,
    section: ConcentricKelvinVoigtSection,
    prestress: StringPrestress,
    other_bindings: &AcousticMaterialBindings<'_>,
) -> Result<CompiledConcentricKelvinVoigtAssembly, MaterialAssemblyError> {
    if other_bindings.string.is_some() {
        return Err(missing("string", "concentric and homogeneous string bindings are mutually exclusive"));
    }
    let input = template.string.as_ref()
        .ok_or_else(|| missing("string", "concentric material binding requires a string description"))?;
    let string = with_concentric_kelvin_voigt_section(input.clone(), &section, prestress)
        .map_err(|error| physical("string", error))?;
    let mut candidate = template.clone();
    candidate.string = Some(string);
    let compiled = compile_material_assembly(&candidate, other_bindings)?;
    Ok(CompiledConcentricKelvinVoigtAssembly { compiled, section, prestress })
}

#[cfg(test)]
mod viscous_tests {
    use super::*;
    use crate::string_specimen::KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY as ETA;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
        PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
    use fs_material::section::ConcentricMaterialLayer;
    use fs_material::state_point::{MaterialPropertySelection, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point};
    use fs_qty::{Density, DynViscosity, Pressure};
    use fs_scenario::{AmbientGas, Listener, Pluck};

    fn section(sleeve_eta: f64, band: (f64, f64)) -> ConcentricKelvinVoigtSection {
        let layers = [(0.00025, 0.0), (0.0005, sleeve_eta)].into_iter().map(|(radius, eta)| {
            let mut claims = ClaimSet::new();
            let mut requirements = Vec::new();
            for (key, value, dims) in [("density", 8000.0, Density::DIMS),
                ("young_modulus", 2e11, Pressure::DIMS), (ETA, eta, DynViscosity::DIMS)] {
                claims.insert_claim(PropertyClaim {
                    key: PropertyKey::new(key, dims), value: PropertyValue::Scalar { value, dims },
                    validity: ValidityDomain::unconstrained().with("omega", band.0, band.1),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
                    provenance: Provenance { source: "synthetic layered loss".into(), license: "test-only".into(), artifact: None },
                }).unwrap();
                requirements.push(ScalarPropertyRequirement::try_new(key, dims, ScalarAdmissibility::Finite).unwrap());
            }
            let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(),
                phase: "solid".into(), process: "test".into(), revision: 0 }, claims, vec![]).unwrap();
            ConcentricMaterialLayer { outer_radius_m: radius, material: resolve_material_state_point(&card,
                &QueryPoint::new().with("omega", (band.0 + band.1) / 2.0).unwrap(), &requirements,
                MaterialPropertySelection::SingleClaimOnly).unwrap() }
        }).collect();
        ConcentricElasticSection::try_new(0.0, layers).unwrap().with_kelvin_voigt_bending(ETA).unwrap()
    }
    fn assembly() -> AcousticAssembly {
        AcousticAssembly { ambient: AmbientGas::sea_level(), string: Some(PrestressedString {
            length_m: 0.65, tension_n: 10.0, lin_density_kg_m: 1.0, axial_stiffness_n: 0.0,
            width_m: 0.01, n_modes: 4, damping_ratio: 0.0, rayleigh: None,
            bending_stiffness_n_m2: 0.0, kelvin_voigt_bending: None, relaxation_bending: None,
            polarization_detune: 0.0, moving_end: false, tension_loss_factor: None, bending_loss_factor: None }), duct: None,
            pluck: Some(Pluck { station_frac: 0.25, height_m: 1e-5 }), bow: None, blow: None,
            reed: None, soundboard: None, body_modes: vec![], plate: None, cavity: None,
            obstacles: vec![], contact_texture: None, listener: Listener { distance_m: 1.0 },
            sample_rate_hz: 128000, duration_s: 0.008 }
    }

    #[test]
    fn replacing_only_sleeve_viscosity_reaches_actual_pressure_without_retuning() {
        let template = assembly();
        let snapshot = template.clone();
        let bindings = AcousticMaterialBindings::default();
        let a = compile_concentric_kelvin_voigt_material_assembly(&template,
            section(1e7, (1.0, 20000.0)), StringPrestress::FixedTension(10.0), &bindings).unwrap();
        let b = compile_concentric_kelvin_voigt_material_assembly(&template,
            section(1e9, (1.0, 20000.0)), StringPrestress::FixedTension(10.0), &bindings).unwrap();
        let sa = a.compiled().assembly().string.as_ref().unwrap();
        let sb = b.compiled().assembly().string.as_ref().unwrap();
        assert_eq!(sa.lin_density_kg_m, sb.lin_density_kg_m);
        assert_eq!(sa.bending_stiffness_n_m2, sb.bending_stiffness_n_m2);
        assert_eq!(sa.tension_n, sb.tension_n);
        assert!(sb.kelvin_voigt_bending.unwrap().viscous_stiffness_n_m2_s
            > sa.kelvin_voigt_bending.unwrap().viscous_stiffness_n_m2_s);
        let ra = a.realize().unwrap();
        let rb = b.realize().unwrap();
        assert_eq!(ra.gas, rb.gas);
        assert_eq!(ra.pressure_pa.len(), 1024);
        assert!(ra.pressure_pa.iter().chain(&rb.pressure_pa).all(|v| v.is_finite()));
        assert!(ra.pressure_pa.iter().any(|v| *v != 0.0));
        assert_ne!(ra.pressure_pa, rb.pressure_pa);
        assert_eq!(template, snapshot);
    }

    #[test]
    fn all_source_layers_and_explicit_zero_loss_survive_compilation() {
        let loss = section(0.0, (1.0, 20000.0));
        let compiled = compile_concentric_kelvin_voigt_material_assembly(&assembly(), loss.clone(),
            StringPrestress::FixedTension(10.0), &AcousticMaterialBindings::default()).unwrap();
        assert_eq!(compiled.section(), &loss);
        let descriptor = compiled.compiled().assembly().string.as_ref().unwrap().kelvin_voigt_bending.unwrap();
        assert_eq!(descriptor.viscous_stiffness_n_m2_s, 0.0);
        assert_eq!(descriptor.material_state_identity, None);
        assert!(compiled.realize().is_ok());
    }

    #[test]
    fn reference_modes_outside_any_layer_support_refuse_before_sound_is_returned() {
        let result = compile_concentric_kelvin_voigt_material_assembly(&assembly(),
            section(1e7, (1.0, 2.0)), StringPrestress::FixedTension(10.0), &AcousticMaterialBindings::default());
        if let Ok(compiled) = result { assert!(compiled.realize().is_err()); }
    }

    #[test]
    fn inherited_internal_loss_is_not_silently_mixed_with_new_material_loss() {
        for (damping, bending_factor) in [(0.01, None), (0.0, Some(0.01))] {
            let mut template = assembly();
            let string = template.string.as_mut().unwrap();
            string.damping_ratio = damping;
            string.bending_loss_factor = bending_factor;
            let snapshot = template.clone();
            assert!(compile_concentric_kelvin_voigt_material_assembly(&template,
                section(1e7, (1.0, 20000.0)), StringPrestress::FixedTension(10.0),
                &AcousticMaterialBindings::default()).is_err());
            assert_eq!(template, snapshot);
        }
    }
}
