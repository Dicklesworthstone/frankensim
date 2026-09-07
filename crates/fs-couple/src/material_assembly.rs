//! Card-driven compilation into the existing acoustic assembly.
//!
//! A source is an exact immutable card, physical query point and explicit claim
//! selection. Only properties needed by the selected mechanical model are
//! resolved. Unbound components remain authored inputs. This compiles frozen
//! specimens; ambient temperature does not silently set material temperature,
//! transport thermal history, or promote raw coefficients to sourced authority.

use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, MaterialStatePointError,
    ORTHOTROPIC_POISSON_RATIO_PROPERTIES, ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES,
    ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES, POISSON_RATIO_PROPERTY, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY,
    resolve_material_state_point,
};
use fs_plate::{PlateMesh, PlateRegion};
use fs_qty::{Density, Dims, DynViscosity, Pressure, Time};
use fs_scenario::{AcousticAssembly, PrestressedString};

use crate::acoustic_realize::{
    AcousticRealizeError, RealizedAssembly, realize_assembly, realize_assembly_with_plate_chart,
};
use crate::string_specimen::{
    EQUILIBRIUM_YOUNG_MODULUS_PROPERTY, KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,
    ResolvedStringSpecimen, StringGeometryConstraint, StringPrestress,
    with_uniform_circular_material_and_constraints,
};
use crate::thin_plate::{
    PlateChartRadiation, PlateMaterialModel, PlateRegionMaterial, PlateThicknessConstraint,
    ResolvedPlateChart, ResolvedPlateSpecimen, compile_plate_material_chart,
    with_uniform_plate_material_state,
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
    /// Tension, extension or target-frequency prescription.
    pub prestress: StringPrestress,
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
}

impl CompiledMaterialAssembly {
    /// Numeric description; a regional chart is additionally retained in `plate`.
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

    /// Run the ordinary acoustic realizer with the compiled members.
    ///
    /// # Errors
    /// Existing assembly, modal-window, excitation and runtime refusals propagate.
    pub fn realize(&self) -> Result<RealizedAssembly, AcousticRealizeError> {
        match &self.plate {
            Some(CompiledMaterialPlate::Regional { chart, radiation }) => {
                realize_assembly_with_plate_chart(&self.assembly, chart.chart(), radiation)
            }
            _ => realize_assembly(&self.assembly),
        }
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
        let specimen = with_uniform_circular_material_and_constraints(
            input.clone(),
            binding.geometry,
            &state,
            binding.prestress,
        )
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
    Ok(CompiledMaterialAssembly {
        assembly,
        string,
        plate,
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
            let state = resolve_plate(source, *model, "plate")?;
            let specimen = with_uniform_plate_material_state(input, &state, *model, *thickness)
                .map_err(|source| physical("plate", source))?;
            assembly.plate = Some(specimen.plate());
            Ok(CompiledMaterialPlate::Uniform(specimen))
        }
        PlateMaterialBinding::Regional {
            mesh,
            boundary_nodes,
            regions,
            radiation,
        } => {
            if assembly.plate.is_some() {
                return Err(missing(
                    "plate",
                    "regional material binding requires the uniform plate to be absent",
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
    component: &str,
) -> Result<ResolvedMaterialStatePoint, MaterialAssemblyError> {
    use ScalarAdmissibility::{Finite, OpenInterval, StrictlyPositive};
    let mut fields = vec![(DENSITY_PROPERTY, Density::DIMS, StrictlyPositive)];
    match model {
        PlateMaterialModel::Isotropic => fields.extend([
            (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, StrictlyPositive),
            (
                POISSON_RATIO_PROPERTY,
                Dims::NONE,
                OpenInterval {
                    lower: -1.0,
                    upper: 0.5,
                },
            ),
        ]),
        PlateMaterialModel::Orthotropic12 | PlateMaterialModel::Orthotropic21 => fields.extend([
            (
                ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES[0],
                Pressure::DIMS,
                StrictlyPositive,
            ),
            (
                ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES[1],
                Pressure::DIMS,
                StrictlyPositive,
            ),
            (ORTHOTROPIC_POISSON_RATIO_PROPERTIES[0], Dims::NONE, Finite),
            (
                ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES[0],
                Pressure::DIMS,
                StrictlyPositive,
            ),
        ]),
    }
    let requirements = fields
        .into_iter()
        .map(|(key, dims, domain)| ScalarPropertyRequirement::try_new(key, dims, domain))
        .collect::<Result<Vec<_>, _>>()
        .map_err(|source| resolution(component, source))?;
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
