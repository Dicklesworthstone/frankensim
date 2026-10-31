//! Source-backed elastic laminate snapshots at prescribed layer temperatures.
//!
//! Each declared layer (or explicitly subdivided thermal cell) is uniform in
//! temperature. Its original claims supply the current elastic constants and
//! the exact integrated free strain from its own reference temperature. Density
//! is resolved only at that reference: refreshing stiffness cannot create mass.
//! Geometry, material axes and the reference surface stay fixed.
//!
//! These are small-strain, held-temperature mechanical snapshots. They do not
//! account for heat/work exchanged while stiffness changes, evolve internal
//! variables, or certify a solved strain against a source calibration envelope.
//! Original source domains and receipts remain available to the consumer. No
//! missing three-dimensional constants or isotropic expansion are inferred.

use super::{LaminaThermalStrain, LaminateThermalError, PreparedLaminateThermal};
use crate::laminate::{
    LaminaLayer, LaminateError, LaminateSection, LayerLaw, MAX_LAMINA_LAYERS, ReferenceSurface,
};
use crate::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, IntegratedIsotropicThermalExpansion,
    IntegratedOrthotropicThermalExpansion, IsotropicElasticParameters,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection,
    MaterialStatePointError, ORTHOTROPIC_POISSON_RATIO_PROPERTIES,
    ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES, ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES,
    ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES, ResolvedIsotropicElasticModuli,
    ResolvedMaterialStatePoint, ResolvedScalarProperty, ScalarAdmissibility,
    ScalarPropertyRequirement, integrate_isotropic_thermal_expansion,
    integrate_orthotropic_thermal_expansion, is_fixed_external_coordinate,
    resolve_material_state_point,
};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Density, Dims, Pressure, Temperature};
use std::sync::Arc;

/// Explicit constitutive interpretation of the original elastic measurements.
#[derive(Clone, Debug, PartialEq)]
pub enum LaminaElasticSource {
    /// An independently selected E/nu, E/G or K/G isotropic parameter pair.
    Isotropic(IsotropicElasticParameters),
    /// Canonical E1,E2,nu12,G12,G31,G23 in the fixed material 1/2 plane.
    /// Requires both measured transverse shear moduli; no E3 or out-of-plane
    /// Poisson ratios are requested or fabricated.
    OrthotropicPlane12,
}

/// Expansion symmetry is selected independently of elastic symmetry.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LaminaExpansionSource {
    /// Canonical scalar linear CTE, integrated by the isotropic source owner.
    Isotropic,
    /// Three canonical principal CTEs, with axes 1/2 in-plane and 3 normal.
    Orthotropic,
}

/// One physical layer or explicitly declared uniform-temperature sublayer.
#[derive(Clone, Debug)]
pub struct LaminaSourceInput {
    /// Shared original immutable card, never copied into synthetic claims.
    pub card: Arc<MaterialCard>,
    /// Complete original condition point, including positive absolute `T` [K].
    pub reference: QueryPoint,
    /// Undeformed thickness [m].
    pub thickness_m: f64,
    /// Fixed principal axis 1 counterclockwise from section x [rad].
    pub orientation_rad: f64,
    /// Explicit original elastic source roles.
    pub elasticity: LaminaElasticSource,
    /// Explicit original expansion source roles.
    pub expansion: LaminaExpansionSource,
    /// Initial selection for the combined density, elasticity and CTE set.
    /// A pinned plan must cover that complete set exactly. Every later query
    /// pins these original claims, regardless of newly available alternatives.
    pub selection: MaterialPropertySelection,
}

/// Typed refusal of the complete preparation; earlier snapshots stay valid.
#[derive(Clone, Debug, PartialEq)]
pub enum LaminateThermalSourceError {
    /// Original source, unit, selection or path-domain refusal for one layer.
    Material {
        layer: usize,
        source: MaterialStatePointError,
    },
    /// Unsupported physical context or layer geometry.
    Invalid {
        layer: usize,
        quantity: &'static str,
    },
    /// One held temperature is required per declared layer.
    TemperatureCount { expected: usize, actual: usize },
    /// Native laminate geometry or constitutive admission failed.
    Section(LaminateError),
    /// Native thermal loading admission failed.
    Thermal(LaminateThermalError),
    /// Cancellation occurred before a complete result was available.
    Cancelled,
}
impl core::fmt::Display for LaminateThermalSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-MATERIAL-LAMINATE-THERMAL-SOURCE: {self:?}")
    }
}
impl std::error::Error for LaminateThermalSourceError {}

/// Original expansion receipts and the complete signed source integral.
#[derive(Clone, Debug, PartialEq)]
pub enum IntegratedLaminaExpansion {
    /// One isotropic source law.
    Isotropic(IntegratedIsotropicThermalExpansion),
    /// Three independently sourced principal laws in the fixed layer frame.
    Orthotropic(IntegratedOrthotropicThermalExpansion),
}
impl IntegratedLaminaExpansion {
    /// Unmodified source bundle at the original reference temperature.
    pub fn reference(&self) -> &ResolvedMaterialStatePoint {
        match self {
            Self::Isotropic(p) => p.reference().resolved(),
            Self::Orthotropic(p) => p.reference(),
        }
    }
    /// Same original claims at this snapshot's actual layer temperature.
    pub fn current(&self) -> &ResolvedMaterialStatePoint {
        match self {
            Self::Isotropic(p) => p.current().resolved(),
            Self::Orthotropic(p) => p.current(),
        }
    }
    fn strain(&self) -> Result<LaminaThermalStrain, LaminateThermalError> {
        match self {
            Self::Isotropic(p) => LaminaThermalStrain::from_isotropic(p),
            Self::Orthotropic(p) => LaminaThermalStrain::from_orthotropic(p),
        }
    }
}

/// Original layer input and immutable reference receipts, shared by snapshots.
#[derive(Clone, Debug)]
pub struct BoundLaminaSource {
    input: LaminaSourceInput,
    reference_material: ResolvedMaterialStatePoint,
    reference_elasticity: ResolvedMaterialStatePoint,
    reference_expansion: IntegratedLaminaExpansion,
}
impl BoundLaminaSource {
    /// Exact input, including the original card, context, geometry and choices.
    pub const fn input(&self) -> &LaminaSourceInput {
        &self.input
    }
    /// Complete original density, elastic and expansion selection receipts.
    pub const fn reference_material(&self) -> &ResolvedMaterialStatePoint {
        &self.reference_material
    }
    /// Original density receipt; this is never refreshed at the current T.
    pub fn reference_density(&self) -> &ResolvedScalarProperty {
        self.reference_material
            .property(DENSITY_PROPERTY)
            .expect("reference density was admitted")
    }
    /// Reference density [kg/m³] used for the fixed areal mass inventory.
    pub fn reference_density_kg_m3(&self) -> f64 {
        self.reference_density().value_si()
    }
    /// Selected original elastic constants, without density or thermal members.
    pub const fn reference_elasticity(&self) -> &ResolvedMaterialStatePoint {
        &self.reference_elasticity
    }
    /// Original admitted zero-length expansion path and its source receipts.
    pub const fn reference_expansion(&self) -> &IntegratedLaminaExpansion {
        &self.reference_expansion
    }
}

/// Bounded source preparation for a fixed ordered reference laminate.
/// Cloning shares original cards, reference bundles and section geometry.
#[derive(Clone, Debug)]
pub struct LaminateThermalSource {
    layers: Arc<[BoundLaminaSource]>,
    reference_section: Arc<LaminateSection>,
}
impl LaminateThermalSource {
    /// Bind original per-layer sources. At most `MAX_LAMINA_LAYERS` are admitted;
    /// a finer thermal discretization must refuse explicitly, never be averaged.
    /// Only `T` can change; every other source coordinate remains fixed.
    pub fn try_new(
        inputs: Vec<LaminaSourceInput>,
        reference: ReferenceSurface,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LaminateThermalSourceError> {
        if cancelled() {
            return Err(LaminateThermalSourceError::Cancelled);
        }
        if inputs.is_empty() || inputs.len() > MAX_LAMINA_LAYERS {
            return Err(LaminateThermalSourceError::Section(
                LaminateError::LayerCount {
                    count: inputs.len(),
                },
            ));
        }
        let mut layers = Vec::with_capacity(inputs.len());
        let mut native = Vec::with_capacity(inputs.len());
        for (layer, input) in inputs.into_iter().enumerate() {
            if cancelled() {
                return Err(LaminateThermalSourceError::Cancelled);
            }
            if !input.thickness_m.is_finite()
                || input.thickness_m <= 0.0
                || !input.orientation_rad.is_finite()
            {
                return Err(LaminateThermalSourceError::Invalid {
                    layer,
                    quantity: "finite orientation and positive reference thickness",
                });
            }
            admit_context(&input.reference, layer)?;
            let material = |source| LaminateThermalSourceError::Material { layer, source };
            let elastic = elastic_requirements(&input.elasticity).map_err(material)?;
            let expansion = expansion_requirements(input.expansion).map_err(material)?;
            let mut requirements = vec![
                ScalarPropertyRequirement::try_new(
                    DENSITY_PROPERTY,
                    Density::DIMS,
                    ScalarAdmissibility::StrictlyPositive,
                )
                .map_err(material)?,
            ];
            requirements.extend_from_slice(&elastic);
            requirements.extend_from_slice(&expansion);
            let reference_material = resolve_material_state_point(
                &input.card,
                &input.reference,
                &requirements,
                input.selection.clone(),
            )
            .map_err(material)?;
            admit_sources(&reference_material, layer)?;
            // A separate elastic bundle prevents density from being re-queried
            // during temperature updates; its claims are selected only once.
            let reference_elasticity = resolve_material_state_point(
                &input.card,
                &input.reference,
                &elastic,
                pin(&reference_material, &elastic).map_err(material)?,
            )
            .map_err(material)?;
            let reference_expansion = integrate(
                &input,
                &input.reference,
                pin(&reference_material, &expansion).map_err(material)?,
            )
            .map_err(material)?;
            let bound = BoundLaminaSource {
                input,
                reference_material,
                reference_elasticity,
                reference_expansion,
            };
            native.push(native_layer(&bound, &bound.reference_elasticity).map_err(material)?);
            layers.push(bound);
        }
        let section = LaminateSection::try_new(native, reference)
            .map_err(LaminateThermalSourceError::Section)?;
        if cancelled() {
            return Err(LaminateThermalSourceError::Cancelled);
        }
        Ok(Self {
            layers: layers.into(),
            reference_section: Arc::new(section),
        })
    }

    /// Sources in original declaration order; sublayers are never merged.
    pub fn layers(&self) -> &[BoundLaminaSource] {
        &self.layers
    }
    /// Original geometry, elastic stiffness and fixed reference mass moments.
    pub fn reference_section(&self) -> &LaminateSection {
        &self.reference_section
    }

    /// Resolve one current absolute temperature [K] per declared layer, then
    /// prepare exact layer-resolved loading and energy. The returned section
    /// contains the refreshed stiffness and must be used for native assembly.
    /// No prior snapshot is mutated if any source, path or layer refuses.
    pub fn at_temperatures(
        &self,
        temperatures_k: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedLaminateThermal, LaminateThermalSourceError> {
        if cancelled() {
            return Err(LaminateThermalSourceError::Cancelled);
        }
        if temperatures_k.len() != self.layers.len() {
            return Err(LaminateThermalSourceError::TemperatureCount {
                expected: self.layers.len(),
                actual: temperatures_k.len(),
            });
        }
        let mut states = Vec::with_capacity(self.layers.len());
        let mut native = Vec::with_capacity(self.layers.len());
        let mut profiles = Vec::with_capacity(self.layers.len());
        for (layer, (source, &temperature)) in self.layers.iter().zip(temperatures_k).enumerate() {
            if cancelled() {
                return Err(LaminateThermalSourceError::Cancelled);
            }
            if !temperature.is_finite() || temperature <= 0.0 {
                return Err(LaminateThermalSourceError::Invalid {
                    layer,
                    quantity: "positive absolute layer temperature",
                });
            }
            let point = source
                .input
                .reference
                .clone()
                .with("T", temperature)
                .map_err(|_| LaminateThermalSourceError::Invalid {
                    layer,
                    quantity: "original temperature coordinate convention",
                })?;
            let material = |source| LaminateThermalSourceError::Material { layer, source };
            let elasticity = source
                .reference_elasticity
                .at_query_point(&source.input.card, &point)
                .map_err(material)?;
            let expansion = integrate(
                &source.input,
                &point,
                pinned_bundle(source.reference_expansion.reference()),
            )
            .map_err(material)?;
            native.push(native_layer(source, &elasticity).map_err(material)?);
            profiles.push(
                expansion
                    .strain()
                    .map_err(LaminateThermalSourceError::Thermal)?,
            );
            states.push(ResolvedLaminaThermal {
                elasticity,
                expansion,
            });
        }
        let section = Arc::new(
            LaminateSection::try_new(native, self.reference_section.reference())
                .map_err(LaminateThermalSourceError::Section)?,
        );
        let prepared = PreparedLaminateThermal::try_new(section, &profiles, &mut cancelled)
            .map_err(|e| match e {
                LaminateThermalError::Cancelled => LaminateThermalSourceError::Cancelled,
                other => LaminateThermalSourceError::Thermal(other),
            })?;
        Ok(ResolvedLaminateThermal {
            source: self.clone(),
            layers: states,
            prepared: Arc::new(prepared),
        })
    }
}

/// Current layer receipts. Reference mass and input are retained by `source()`
/// on the containing snapshot at this same layer index.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedLaminaThermal {
    elasticity: ResolvedMaterialStatePoint,
    expansion: IntegratedLaminaExpansion,
}
impl ResolvedLaminaThermal {
    /// Original elastic claims evaluated at the actual held layer temperature.
    pub const fn elasticity(&self) -> &ResolvedMaterialStatePoint {
        &self.elasticity
    }
    /// Complete signed expansion integral and endpoint receipts.
    pub const fn expansion(&self) -> &IntegratedLaminaExpansion {
        &self.expansion
    }
}

/// Complete immutable current section, thermal loading and original sources.
#[derive(Clone, Debug)]
pub struct ResolvedLaminateThermal {
    source: LaminateThermalSource,
    layers: Vec<ResolvedLaminaThermal>,
    prepared: Arc<PreparedLaminateThermal>,
}
impl ResolvedLaminateThermal {
    /// Shared original cards, layer context, geometry and density receipts.
    pub const fn source(&self) -> &LaminateThermalSource {
        &self.source
    }
    /// Current receipts in the same order as `source().layers()`.
    pub fn layers(&self) -> &[ResolvedLaminaThermal] {
        &self.layers
    }
    /// Exact current section and integrated loading for native shell assembly.
    pub const fn prepared(&self) -> &Arc<PreparedLaminateThermal> {
        &self.prepared
    }
    /// Transfer the prepared mechanics. Keep this snapshot separately whenever
    /// complete expansion-path receipts are needed by the consuming application.
    pub fn into_prepared(self) -> Arc<PreparedLaminateThermal> {
        self.prepared
    }
}

fn admit_context(point: &QueryPoint, layer: usize) -> Result<(), LaminateThermalSourceError> {
    if !point
        .axes()
        .get("T")
        .is_some_and(|t| t.is_finite() && *t > 0.0)
        || point.axes().iter().any(|(axis, &value)| {
            let quantity = point.axis_quantities().get(axis).copied();
            !is_fixed_external_coordinate(axis, quantity, Some(value))
                || axis != "T" && quantity.is_some_and(|q| q.dims() == Temperature::DIMS)
        })
    {
        return Err(LaminateThermalSourceError::Invalid {
            layer,
            quantity: "one positive absolute T and otherwise fixed external context",
        });
    }
    Ok(())
}

fn admit_sources(
    state: &ResolvedMaterialStatePoint,
    layer: usize,
) -> Result<(), LaminateThermalSourceError> {
    if state
        .properties()
        .iter()
        .any(|p| !p.has_only_fixed_external_conditions(Some("T"), false))
    {
        return Err(LaminateThermalSourceError::Invalid {
            layer,
            quantity: "held-temperature elastic, density and expansion sources",
        });
    }
    Ok(())
}

fn elastic_requirements(
    source: &LaminaElasticSource,
) -> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
    match source {
        LaminaElasticSource::Isotropic(parameters) => {
            Ok(parameters.requirements().into_iter().cloned().collect())
        }
        LaminaElasticSource::OrthotropicPlane12 => [
            (
                ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES[0],
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            ),
            (
                ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES[1],
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            ),
            (
                ORTHOTROPIC_POISSON_RATIO_PROPERTIES[0],
                Dims::NONE,
                ScalarAdmissibility::Finite,
            ),
            (
                ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES[0],
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            ),
            (
                ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES[2],
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            ),
            (
                ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES[1],
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            ),
        ]
        .into_iter()
        .map(|(key, dims, domain)| ScalarPropertyRequirement::try_new(key, dims, domain))
        .collect(),
    }
}

fn expansion_requirements(
    source: LaminaExpansionSource,
) -> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
    let names: &[&str] = match source {
        LaminaExpansionSource::Isotropic => &[LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY],
        LaminaExpansionSource::Orthotropic => &ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES,
    };
    names
        .iter()
        .map(|&name| {
            ScalarPropertyRequirement::try_new(
                name,
                INVERSE_TEMPERATURE_DIMS,
                ScalarAdmissibility::Finite,
            )
        })
        .collect()
}

fn pin(
    state: &ResolvedMaterialStatePoint,
    requirements: &[ScalarPropertyRequirement],
) -> Result<MaterialPropertySelection, MaterialStatePointError> {
    let mut selected = Vec::with_capacity(requirements.len());
    for requirement in requirements {
        let property =
            state
                .property(requirement.name())
                .ok_or(MaterialStatePointError::InvalidDerived {
                    quantity: "distinct scalar laminate source roles",
                })?;
        selected.push((
            requirement.name().to_owned(),
            property.answer().receipt.selected,
        ));
    }
    Ok(MaterialPropertySelection::PinnedByProperty(selected))
}

fn pinned_bundle(state: &ResolvedMaterialStatePoint) -> MaterialPropertySelection {
    MaterialPropertySelection::PinnedByProperty(
        state
            .properties()
            .iter()
            .map(|p| {
                (
                    p.requirement().name().to_owned(),
                    p.answer().receipt.selected,
                )
            })
            .collect(),
    )
}

fn integrate(
    input: &LaminaSourceInput,
    current: &QueryPoint,
    selection: MaterialPropertySelection,
) -> Result<IntegratedLaminaExpansion, MaterialStatePointError> {
    match input.expansion {
        LaminaExpansionSource::Isotropic => {
            integrate_isotropic_thermal_expansion(&input.card, &input.reference, current, selection)
                .map(IntegratedLaminaExpansion::Isotropic)
        }
        LaminaExpansionSource::Orthotropic => integrate_orthotropic_thermal_expansion(
            &input.card,
            &input.reference,
            current,
            selection,
        )
        .map(IntegratedLaminaExpansion::Orthotropic),
    }
}

fn native_layer(
    source: &BoundLaminaSource,
    resolved: &ResolvedMaterialStatePoint,
) -> Result<LaminaLayer, MaterialStatePointError> {
    let (e, nu, g, transverse) = match &source.input.elasticity {
        LaminaElasticSource::Isotropic(parameters) => {
            let m =
                ResolvedIsotropicElasticModuli::try_from_resolved(resolved.clone(), parameters)?;
            (
                [m.young_modulus_pa(); 2],
                m.poisson_ratio(),
                m.shear_modulus_pa(),
                [m.shear_modulus_pa(); 2],
            )
        }
        LaminaElasticSource::OrthotropicPlane12 => {
            let value = |name| {
                resolved
                    .property(name)
                    .expect("plane constants were resolved")
                    .value_si()
            };
            (
                [
                    value(ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES[0]),
                    value(ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES[1]),
                ],
                value(ORTHOTROPIC_POISSON_RATIO_PROPERTIES[0]),
                value(ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES[0]),
                [
                    value(ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES[2]),
                    value(ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES[1]),
                ],
            )
        }
    };
    // Reciprocal plane-stress compliance, in engineering shear. These derived
    // entries are not new source claims or a fabricated 3-D strain envelope.
    let denominator = 1.0 - nu * (nu * e[1] / e[0]);
    if !denominator.is_finite() || denominator <= 0.0 {
        return Err(MaterialStatePointError::InvalidDerived {
            quantity: "stable reciprocal layer compliance",
        });
    }
    let cross = nu * e[1] / denominator;
    Ok(LaminaLayer {
        thickness_m: source.input.thickness_m,
        orientation_rad: source.input.orientation_rad,
        density_kg_m3: source.reference_density_kg_m3(),
        material: Some(resolved.clone()),
        law: LayerLaw::PlaneStress {
            stiffness_pa: [
                e[0] / denominator,
                cross,
                0.0,
                cross,
                e[1] / denominator,
                0.0,
                0.0,
                0.0,
                g,
            ],
            transverse_shear_pa: transverse,
        },
    })
}
