//! Card-driven compilation into the existing acoustic assembly.
//!
//! A source is an exact immutable card, physical query point and explicit claim
//! selection. Only properties needed by the selected mechanical model are
//! resolved. Unbound components remain authored inputs. This compiles frozen
//! specimens; ambient temperature does not silently set material temperature,
//! transport thermal history, or promote raw coefficients to sourced authority.

pub mod concentric;
/// Source-bound close-coiled winding on the same generic string/body solver.
pub mod wound;
/// Direct measured whole-string properties without invented core/winding laws.
pub mod measured_string;
/// Material-reduced reeds composed with the existing solid/gas assembly.
pub mod valve;
/// Source-backed duct walls shared with graph and DG material consumers.
pub mod duct;

use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::gas::GasState;
use fs_material::plane_stress_viscosity::PLANE_STRESS_VISCOSITY_PROPERTIES;
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY,
    integrate_isotropic_thermal_expansion,
    resolve_isotropic_thermoelastic_state_point_with_requirements, resolve_material_state_point,
};
use fs_plate::{PlateMesh, PlateRegion};
use fs_qty::{Density, DynViscosity, Pressure, Time};
use fs_scenario::{AcousticAssembly, PrestressedString};

use crate::acoustic_realize::{
    AcousticRealizeError, RealizeOptions, RealizedAssembly,
    realize_assembly_with_material_plate_chart_and_options, realize_assembly_with_options,
};
use crate::string_specimen::{
    EQUILIBRIUM_YOUNG_MODULUS_PROPERTY, KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,
    ResolvedStringSpecimen, StringGeometryConstraint, StringPrestress,
    with_uniform_circular_material_and_constraints, with_uniform_circular_thermal_extension,
};
use crate::thin_plate::{
    PlateChartRadiation, PlateMaterialModel, PlateRegionMaterial, PlateThicknessConstraint,
    ResolvedPlateChart, ResolvedPlateSpecimen, compile_plate_material_chart,
    with_uniform_isotropic_thermoelastic_material_state, with_uniform_plate_material_state,
};

/// Exact material data and query conditions, independent of the ambient gas.
#[derive(Clone, Debug)]
pub struct MaterialSource<'a> {
    /// Immutable material card; its content identity is retained by resolution.
    pub card: &'a MaterialCard,
    /// Physical state coordinates, with their original quantity schemas.
    pub point: &'a QueryPoint,
    /// Explicit ambiguity policy. Pins must cover exactly the consumer's needs.
    pub selection: MaterialPropertySelection,
}

/// Material and independent mechanical constraints for the described string.
#[derive(Clone, Debug)]
pub struct StringMaterialBinding<'a> {
    /// Source for density, elasticity and any selected sourced bending loss.
    pub source: MaterialSource<'a>,
    /// Current circular radius or total specimen mass.
    pub geometry: StringGeometryConstraint,
    /// Mechanical prescription or expansion from the same card and current state.
    pub prestress: StringPrestressBinding<'a>,
}

/// How the described supports establish string tension at the resolved state.
#[derive(Clone, Debug)]
pub enum StringPrestressBinding<'a> {
    /// Explicit mechanical constraint; any raw thermal strain remains authored.
    Prescribed(StringPrestress),
    /// Fixed supports with free strain integrated from the card's alpha(T).
    /// The source card and current point come from `StringMaterialBinding::source`.
    /// This selects frozen small-strain states, not a heat-transport history.
    ThermalExtension {
        /// State where the reference stress-free length is specified. All other
        /// coordinates and quantity conventions must match the current point.
        reference_point: &'a QueryPoint,
        /// Select only the expansion claim; mechanical/loss pins stay in `source`.
        expansion_selection: MaterialPropertySelection,
        /// Stress-free length at the reference temperature [m].
        reference_stress_free_length_m: f64,
        /// Declared small-strain limit; slack/compression is unsupported.
        linear_strain_limit: f64,
    },
}

/// One plate region before material-property resolution.
#[derive(Clone, Debug)]
pub struct PlateRegionSource<'a> {
    /// Unique region name and its triangles in the supplied mesh.
    pub region: PlateRegion,
    /// Exact card and physical state for this region.
    pub source: MaterialSource<'a>,
    /// In-plane elasticity and principal-axis mapping.
    pub model: PlateMaterialModel,
    /// Regional thickness or mass constraint.
    pub thickness: PlateThicknessConstraint,
    /// Counterclockwise material-axis angle in the chart plane [rad].
    pub material_angle_rad: f64,
}

/// Loss selection for all regions of a source-bound plate.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum RegionalPlateLoss {
    /// Keep the independently authored `PlateChartRadiation::damping_ratio`.
    /// No loss property is inferred from the elastic material or its name.
    #[default]
    Authored,
    /// Resolve Kelvin-Voigt bending from every region's same material state.
    /// Isotropic regions need `kelvin_voigt_bending_viscosity`; orthotropic
    /// regions need all six reduced `plane_stress_viscosity_*` coefficients.
    /// The existing weak-loss modal image retains each mode's curvature-energy
    /// weighting and source band, but not off-diagonal modal damping. Authored
    /// damping must be zero to avoid counting the same loss twice.
    KelvinVoigt,
}

/// Uniform rectangular or regional chart material prescription.
#[derive(Clone, Debug)]
pub enum PlateMaterialBinding<'a> {
    /// Bind the existing uniform plate description; its geometry and controls stay.
    Uniform {
        /// Exact card and physical state.
        source: MaterialSource<'a>,
        /// In-plane elasticity and principal-axis mapping.
        model: PlateMaterialModel,
        /// Uniform thickness or total plate mass constraint.
        thickness: PlateThicknessConstraint,
    },
    /// Compile a linear chart; the assembly's uniform plate must be absent.
    Regional {
        /// Flat plate mid-surface geometry and connectivity.
        mesh: PlateMesh,
        /// Support-node set in the supplied mesh.
        boundary_nodes: Vec<usize>,
        /// Complete disjoint material partition, without a fallback region.
        regions: Vec<PlateRegionSource<'a>>,
        /// Explicit material-loss choice, retained through acoustic realization.
        loss: RegionalPlateLoss,
        /// Existing modal, mechanical and force-footprint controls.
        radiation: PlateChartRadiation,
    },
}

/// Explicit source-backed replacements for selected assembly members.
#[derive(Clone, Debug, Default)]
pub struct AcousticMaterialBindings<'a> {
    /// Requires an existing string description; replaces its derived coefficients.
    pub string: Option<StringMaterialBinding<'a>>,
    /// Uniform or regional plate binding. Other body modes remain authored.
    pub plate: Option<PlateMaterialBinding<'a>>,
    /// Uniform or complete segment wall sources; the receiver remains explicit.
    pub duct_wall: Option<duct::DuctWallMaterialBindings<'a>>,
}

/// Retained plate specimen and any chart-specific realization controls.
#[derive(Clone, Debug)]
pub enum CompiledMaterialPlate {
    /// Uniform source-bound specimen.
    Uniform(ResolvedPlateSpecimen),
    /// Source-bound regional chart and its explicit radiation controls.
    Regional {
        /// Immutable compiled geometry, assignments and material receipts.
        chart: ResolvedPlateChart,
        /// Modal window, support/prestress, damping and force footprint.
        radiation: PlateChartRadiation,
    },
}

/// Immutable numeric assembly with retained source-bound specimens.
/// It does not certify unbound members or the complete acoustic approximation.
#[derive(Clone, Debug)]
pub struct CompiledMaterialAssembly {
    assembly: AcousticAssembly,
    string: Option<ResolvedStringSpecimen>,
    plate: Option<CompiledMaterialPlate>,
    duct_wall: Option<duct::MaterialWallLayout>,
}

impl CompiledMaterialAssembly {
    /// Numeric description; regional plates and material walls are retained in
    /// their respective accessors and must use this compiled realizer.
    #[must_use]
    pub const fn assembly(&self) -> &AcousticAssembly {
        &self.assembly
    }

    /// Bound string and exact material receipts, if explicitly requested.
    #[must_use]
    pub const fn string(&self) -> Option<&ResolvedStringSpecimen> {
        self.string.as_ref()
    }

    /// Bound plate and exact material receipts, if explicitly requested.
    #[must_use]
    pub const fn plate(&self) -> Option<&CompiledMaterialPlate> {
        self.plate.as_ref()
    }

    /// Original wall layout, with each source, spectrum, thickness and applicability.
    #[must_use]
    pub const fn duct_wall(
        &self,
    ) -> Option<&duct::MaterialWallLayout> {
        self.duct_wall.as_ref()
    }

    /// Run the ordinary acoustic realizer with the compiled members.
    ///
    /// # Errors
    /// Existing assembly, modal-window, excitation and runtime refusals propagate.
    /// A source-bound duct wall requires an explicit physical observer through
    /// `realize_with_options`; no receiver geometry is invented here.
    pub fn realize(&self) -> Result<RealizedAssembly, AcousticRealizeError> {
        self.realize_with_options(&RealizeOptions::default())
    }

    /// Render the retained material specimens with explicit medium, observer
    /// and numerical controls. Uniform and regional plates share the same
    /// options path; choosing a material chart cannot reset the requested bore
    /// resolution or silently discard an air-specific absorption declaration.
    /// Source receipts and solid-state constraints remain those of compilation.
    ///
    /// # Errors
    /// Material-chart, gas, absorption and numerical refusals propagate from
    /// the common acoustic realizer without changing the compiled specimen.
    pub fn realize_with_options(
        &self,
        options: &RealizeOptions,
    ) -> Result<RealizedAssembly, AcousticRealizeError> {
        let observer = match (&self.duct_wall, options.baffled_duct) {
            (Some(_), None) => return Err(AcousticRealizeError::InvalidDescription {
                what: "compiled duct wall requires an explicit physical baffled-duct observer",
            }),
            (Some(_), Some(binding)) if binding.material_wall.is_some() => {
                return Err(AcousticRealizeError::InvalidDescription {
                    what: "compiled duct wall cannot also select an observer wall override",
                });
            }
            (Some(wall), Some(binding)) => Some(
                crate::acoustic_realize::duct_observer::BaffledDuctObserver {
                    material_wall: Some(wall),
                    ..binding
                },
            ),
            (None, binding) => binding,
        };
        let options = RealizeOptions { baffled_duct: observer, ..*options };
        match &self.plate {
            Some(CompiledMaterialPlate::Regional { chart, radiation }) => {
                realize_assembly_with_material_plate_chart_and_options(
                    &self.assembly,
                    chart,
                    radiation,
                    &options,
                )
            }
            _ => realize_assembly_with_options(&self.assembly, &options),
        }
    }

    /// Run the compiled solid-material assembly in an explicitly supplied gas.
    ///
    /// This composes the two independent material ingress paths instead of
    /// forcing a caller to choose between source-bound string/plate mechanics
    /// and a source-/mixture-derived acoustic medium. String and plate
    /// coefficients remain exactly those retained by this compilation, while
    /// density, sound speed, viscosity, conductivity, heat capacity and
    /// impedance come from `gas`. The authored scenario ambient is ignored by
    /// the acoustic realizer on this path.
    ///
    /// A regional plate keeps its compiled material chart and force/radiation
    /// controls; only the surrounding acoustic medium is replaced. No material
    /// name inference, coefficient copying, or alternate solver is introduced.
    ///
    /// # Errors
    /// Existing assembly/chart/runtime refusals plus the supplied-gas coherence
    /// checks propagate unchanged.
    pub fn realize_with_gas(
        &self,
        gas: &GasState,
    ) -> Result<RealizedAssembly, AcousticRealizeError> {
        self.realize_with_options(&RealizeOptions {
            gas: Some(*gas),
            ..RealizeOptions::default()
        })
    }
}

/// Component-local refusal, retaining the material or physics owner's error.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialAssemblyError {
    /// The exact card/query/requirements could not be resolved.
    Resolution {
        /// String, uniform plate, or named plate region.
        component: String,
        /// Original material-state resolution error.
        source: MaterialStatePointError,
    },
    /// Geometry, constraints, law selection or assignment failed binding.
    Binding {
        /// String, uniform plate, or named plate region.
        component: String,
        /// Original physical binding error.
        source: AcousticRealizeError,
    },
}

impl core::fmt::Display for MaterialAssemblyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Resolution { component, source } => {
                write!(f, "FS-COUPLE-MATERIAL-RESOLUTION ({component}): {source}")
            }
            Self::Binding { component, source } => {
                write!(f, "FS-COUPLE-MATERIAL-BINDING ({component}): {source}")
            }
        }
    }
}

impl std::error::Error for MaterialAssemblyError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Resolution { source, .. } => Some(source),
            Self::Binding { source, .. } => Some(source),
        }
    }
}

/// Resolve cards and compile selected members through the existing specimen APIs.
/// Source values replace only the bound members' derived mechanical coefficients.
/// Geometry, prescribed constraints, supports, excitations, loss choices and
/// numerical controls remain explicit. No property is inferred from a material
/// name. Failed compilation leaves the input and any previous result unchanged.
///
/// This synchronous compiler carries neither thermal evolution nor cancellation
/// guarantees. It is a Rust description API, not a serialized scenario format.
///
/// # Errors
    /// Missing members, unsupported laws, invalid assignments or card/physical
/// refusals propagate with component context. Realization admission runs later.
pub fn compile_material_assembly(
    template: &AcousticAssembly,
    bindings: &AcousticMaterialBindings<'_>,
) -> Result<CompiledMaterialAssembly, MaterialAssemblyError> {
    let mut assembly = template.clone();
    let string = if let Some(binding) = &bindings.string {
        let input = template
            .string
            .as_ref()
            .ok_or_else(|| missing("string", "material binding requires a string description"))?;
        let requirements =
            string_requirements(input).map_err(|source| resolution("string", source))?;
        let state = resolve(&binding.source, &requirements, "string")?;
        let specimen = match &binding.prestress {
            StringPrestressBinding::Prescribed(prestress) => {
                with_uniform_circular_material_and_constraints(
                    input.clone(),
                    binding.geometry,
                    &state,
                    *prestress,
                )
            }
            StringPrestressBinding::ThermalExtension {
                reference_point,
                expansion_selection,
                reference_stress_free_length_m,
                linear_strain_limit,
            } => {
                let expansion = integrate_isotropic_thermal_expansion(
                    binding.source.card,
                    reference_point,
                    binding.source.point,
                    expansion_selection.clone(),
                )
                .map_err(|source| resolution("string", source))?;
                with_uniform_circular_thermal_extension(
                    input.clone(),
                    binding.geometry,
                    &state,
                    &expansion,
                    *reference_stress_free_length_m,
                    *linear_strain_limit,
                )
            }
        }
        .map_err(|source| physical("string", source))?;
        assembly.string = Some(specimen.string());
        Some(specimen)
    } else {
        None
    };
    let plate = bindings
        .plate
        .as_ref()
        .map(|binding| compile_plate(&mut assembly, binding))
        .transpose()?;
    let duct_wall = bindings.duct_wall.as_ref()
        .map(|binding| duct::compile(template, binding)).transpose()?;
    Ok(CompiledMaterialAssembly {
        assembly,
        string,
        plate,
        duct_wall,
    })
}

fn compile_plate(
    assembly: &mut AcousticAssembly,
    binding: &PlateMaterialBinding<'_>,
) -> Result<CompiledMaterialPlate, MaterialAssemblyError> {
    match binding {
        PlateMaterialBinding::Uniform {
            source,
            model,
            thickness,
        } => {
            let input = assembly.plate.ok_or_else(|| {
                missing(
                    "plate",
                    "uniform material binding requires a plate description",
                )
            })?;
            let specimen = if input.thermoelastic.is_some() {
                if *model != PlateMaterialModel::Isotropic {
                    return Err(missing(
                        "plate",
                        "sourced thermal plate loss requires the original isotropic Young/Poisson material bundle",
                    ));
                }
                let additional = if input.kelvin_voigt_bending.is_some() {
                    vec![
                        ScalarPropertyRequirement::try_new(
                            KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,
                            DynViscosity::DIMS,
                            ScalarAdmissibility::NonNegative,
                        )
                        .map_err(|source| resolution("plate", source))?,
                    ]
                } else {
                    Vec::new()
                };
                let state = resolve_isotropic_thermoelastic_state_point_with_requirements(
                    source.card,
                    source.point,
                    source.selection.clone(),
                    &additional,
                )
                .map_err(|source| resolution("plate", source))?;
                with_uniform_isotropic_thermoelastic_material_state(input, &state, *thickness)
            } else {
                let state = resolve_plate(
                    source,
                    *model,
                    input.kelvin_voigt_bending.is_some(),
                    "plate",
                )?;
                with_uniform_plate_material_state(input, &state, *model, *thickness)
            }
            .map_err(|source| physical("plate", source))?;
            assembly.plate = Some(specimen.plate());
            Ok(CompiledMaterialPlate::Uniform(specimen))
        }
        PlateMaterialBinding::Regional {
            mesh,
            boundary_nodes,
            regions,
            loss,
            radiation,
        } => {
            if assembly.plate.is_some() {
                return Err(missing(
                    "plate",
                    "regional material binding requires the uniform plate to be absent",
                ));
            }
            let bending_viscosity = *loss == RegionalPlateLoss::KelvinVoigt;
            if bending_viscosity && radiation.damping_ratio != 0.0 {
                return Err(missing(
                    "plate",
                    "source-bound regional bending loss requires zero authored damping",
                ));
            }
            let inputs = regions
                .iter()
                .map(|region| {
                    Ok(PlateRegionMaterial {
                        region: region.region.clone(),
                        material: resolve_plate(
                            &region.source,
                            region.model,
                            bending_viscosity,
                            &format!("plate/{}", region.region.name),
                        )?,
                        model: region.model,
                        thickness_constraint: region.thickness,
                        material_angle_rad: region.material_angle_rad,
                    })
                })
                .collect::<Result<Vec<_>, MaterialAssemblyError>>()?;
            let chart = compile_plate_material_chart(mesh.clone(), boundary_nodes.clone(), inputs)
                .map_err(|source| physical("plate", source))?;
            let chart = if bending_viscosity {
                chart
                    .with_kelvin_voigt_bending_loss()
                    .map_err(|source| physical("plate", source))?
            } else {
                chart
            };
            Ok(CompiledMaterialPlate::Regional {
                chart,
                radiation: radiation.clone(),
            })
        }
    }
}

fn string_requirements(
    input: &PrestressedString,
) -> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
    use ScalarAdmissibility::{NonNegative, StrictlyPositive};
    let mut requirements = vec![
        ScalarPropertyRequirement::try_new(DENSITY_PROPERTY, Density::DIMS, StrictlyPositive)?,
        ScalarPropertyRequirement::try_new(
            YOUNG_MODULUS_PROPERTY,
            Pressure::DIMS,
            StrictlyPositive,
        )?,
    ];
    if input.kelvin_voigt_bending.is_some() {
        requirements.push(ScalarPropertyRequirement::try_new(
            KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,
            DynViscosity::DIMS,
            NonNegative,
        )?);
    }
    if let Some(law) = &input.relaxation_bending {
        let pairs = law.source_properties.as_ref().ok_or_else(|| {
            MaterialStatePointError::InvalidRequirement {
                property: "relaxation_bending".into(),
                reason: "rebinding a relaxation spectrum requires explicit source-property pairs",
            }
        })?;
        requirements.push(ScalarPropertyRequirement::try_new(
            EQUILIBRIUM_YOUNG_MODULUS_PROPERTY,
            Pressure::DIMS,
            StrictlyPositive,
        )?);
        for pair in pairs {
            for requirement in [
                ScalarPropertyRequirement::try_new(&pair.modulus, Pressure::DIMS, NonNegative)?,
                ScalarPropertyRequirement::try_new(
                    &pair.relaxation_time,
                    Time::DIMS,
                    StrictlyPositive,
                )?,
            ] {
                if !requirements.contains(&requirement) {
                    requirements.push(requirement);
                }
            }
        }
    }
    Ok(requirements)
}

fn resolve_plate(
    source: &MaterialSource<'_>,
    model: PlateMaterialModel,
    bending_viscosity: bool,
    component: &str,
) -> Result<ResolvedMaterialStatePoint, MaterialAssemblyError> {
    let mut requirements = model.elastic_requirements()
        .map_err(|source| resolution(component, source))?;
    if bending_viscosity {
        if model.is_isotropic() {
            requirements.push(ScalarPropertyRequirement::try_new(
                KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,
                DynViscosity::DIMS,
                ScalarAdmissibility::NonNegative,
            ).map_err(|source| resolution(component, source))?);
        } else {
            // Couplings may be signed. The shared plane-stress owner admits
            // the complete tensor, including explicit zero coefficients.
            for key in PLANE_STRESS_VISCOSITY_PROPERTIES {
                requirements.push(ScalarPropertyRequirement::try_new(key, DynViscosity::DIMS,
                    ScalarAdmissibility::Finite).map_err(|source| resolution(component, source))?);
            }
        }
    }
    resolve(source, &requirements, component)
}

fn resolve(
    source: &MaterialSource<'_>,
    requirements: &[ScalarPropertyRequirement],
    component: &str,
) -> Result<ResolvedMaterialStatePoint, MaterialAssemblyError> {
    resolve_material_state_point(
        source.card,
        source.point,
        requirements,
        source.selection.clone(),
    )
    .map_err(|source| resolution(component, source))
}

fn resolution(component: &str, source: MaterialStatePointError) -> MaterialAssemblyError {
    MaterialAssemblyError::Resolution {
        component: component.into(),
        source,
    }
}

fn physical(component: &str, source: AcousticRealizeError) -> MaterialAssemblyError {
    MaterialAssemblyError::Binding {
        component: component.into(),
        source,
    }
}

fn missing(component: &str, what: &'static str) -> MaterialAssemblyError {
    physical(component, AcousticRealizeError::InvalidDescription { what })
}
