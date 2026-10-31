//! Source-selected RC fiber laws and geometry-consistent section mass.
//!
//! Every constituent keeps one immutable material-state bundle. Selecting source
//! coefficients does not calibrate the Mander or Menegotto-Pinto cyclic model.
//! Temperature and other external conditions are prescribed and fixed; sampled
//! strain, rate, frequency or history dependence cannot be frozen into a law.
//!
//! Reinforcement is represented by two symmetric, full-core-width thin bands.
//! Their specified areas replace concrete at the specified centroids. The steel
//! response is lumped at those centroids, as in the existing fiber formulation.
//! This is a smeared reinforcement model, not resolved round-bar geometry.

use super::{Fiber, FiberLaw, Section};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::MaterialError;
use fs_material::creep::uniaxial::UniaxialNorton;
use fs_material::fiber::{ManderConcrete, MenegottoPintoSteel};
use fs_material::plastic::chaboche::axial::UniaxialChaboche;
use fs_material::plastic::chaboche::source::{
    BackstressProperties, ResolvedChabocheBranches, ResolvedChabocheJ2,
    chaboche_branch_requirements, chaboche_requirements,
};
use fs_material::plastic::reduction::{UniaxialJ2, UniaxialVoce};
use fs_material::plastic::source::{
    HardeningPropertyBasis, J2SourceError, ResolvedJ2Plasticity, j2_requirements,
};
use fs_material::plastic::voce::source::{ResolvedVoceJ2, voce_requirements};
use fs_material::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::{Density, Dims, Pressure};

/// Explicit property roles for a Mander material; no strength-to-modulus or
/// confinement relationship is inferred.
#[derive(Clone, Debug)]
pub struct ConcreteFiberProperties {
    /// Reference mass density in kg/m³.
    pub density: String,
    /// Positive peak compressive stress fcc in Pa.
    pub compressive_strength: String,
    /// Positive compressive strain magnitude at fcc.
    pub peak_strain: String,
    /// Initial concrete modulus Ec in Pa.
    pub initial_modulus: String,
    /// Positive ultimate compressive strain magnitude.
    pub ultimate_strain: String,
}

/// Explicit property roles for every Menegotto-Pinto parameter, including its
/// cyclic curvature law. An elastic/yield pair alone does not identify this law.
#[derive(Clone, Debug)]
pub struct SteelFiberProperties {
    /// Reference mass density in kg/m³.
    pub density: String,
    /// Initial modulus E0 in Pa.
    pub young_modulus: String,
    /// Positive yield stress fy in Pa.
    pub yield_stress: String,
    /// Post-yield to initial tangent ratio b.
    pub hardening_ratio: String,
    /// Initial transition curvature R0.
    pub initial_curvature: String,
    /// Curvature degradation coefficient a1, in [0, R0).
    pub curvature_degradation: String,
    /// Positive curvature evolution scale a2.
    pub curvature_scale: String,
}

/// Admission failures preserve the owning material resolver's diagnostics.
#[derive(Debug)]
pub enum RcSectionError {
    /// Material-state selection or admissibility failed.
    Source(MaterialStatePointError),
    /// The owning constitutive law rejected its parameters.
    Law(MaterialError),
    /// The native J2 material-source owner rejected its constitutive reduction.
    Plasticity(J2SourceError),
    /// Geometry, density or numerical representation was inadmissible.
    Invalid(&'static str),
}
impl core::fmt::Display for RcSectionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Source(error) => error.fmt(f),
            Self::Law(error) => error.fmt(f),
            Self::Plasticity(error) => error.fmt(f),
            Self::Invalid(message) => write!(f, "RC section: {message}"),
        }
    }
}
impl std::error::Error for RcSectionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Source(error) => Some(error),
            Self::Law(error) => Some(error),
            Self::Plasticity(error) => Some(error),
            Self::Invalid(_) => None,
        }
    }
}
impl From<MaterialStatePointError> for RcSectionError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Source(error)
    }
}
impl From<MaterialError> for RcSectionError {
    fn from(error: MaterialError) -> Self {
        Self::Law(error)
    }
}
impl From<J2SourceError> for RcSectionError {
    fn from(error: J2SourceError) -> Self {
        Self::Plasticity(error)
    }
}
type Result<T> = core::result::Result<T, RcSectionError>;

fn require(admitted: bool, message: &'static str) -> Result<()> {
    if admitted {
        Ok(())
    } else {
        Err(RcSectionError::Invalid(message))
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

/// A checked immutable uniaxial law and its matching reference density.
/// Authored parameters carry no source bundle; resolved parameters retain every
/// selected property, its query conditions, and its original usage receipt.
#[derive(Clone, Debug)]
pub struct FiberMaterial {
    law: FiberLaw,
    density_kg_m3: f64,
    source: Option<ResolvedMaterialStatePoint>,
    chaboche_branch_properties: Option<Vec<BackstressProperties>>,
}
impl FiberMaterial {
    /// Explicit authored law and density, without implied source authority.
    pub fn authored(law: FiberLaw, density_kg_m3: f64) -> Result<Self> {
        match law {
            FiberLaw::Norton(law) => Self::norton(law, density_kg_m3),
            FiberLaw::J2(law) => Self::j2(law, density_kg_m3),
            FiberLaw::Voce(law) => Self::voce(law, density_kg_m3),
            FiberLaw::Chaboche(law) => Self::chaboche(law, density_kg_m3),
            FiberLaw::ChabocheBranches(law) => {
                require(
                    positive(density_kg_m3),
                    "Chaboche fiber density must be finite and positive",
                )?;
                Ok(Self {
                    law: FiberLaw::ChabocheBranches(law),
                    density_kg_m3,
                    source: None,
                    chaboche_branch_properties: None,
                })
            }
            FiberLaw::Concrete(law) => Self::concrete(law, density_kg_m3),
            FiberLaw::Steel(law) => Self::steel(law, density_kg_m3),
        }
    }

    /// Native axial creep with explicitly authored reference density. A sourced
    /// law retains its original kinetic/elastic receipts inside `law`; this
    /// constructor does not promote the separately supplied density to data.
    pub fn norton(law: UniaxialNorton, density_kg_m3: f64) -> Result<Self> {
        require(
            positive(density_kg_m3),
            "creep fiber density must be finite and positive",
        )?;
        Ok(Self {
            law: FiberLaw::Norton(law),
            density_kg_m3,
            source: None,
            chaboche_branch_properties: None,
        })
    }

    /// An authored native J2 law for metal fibers, with explicit reference mass.
    /// This uses the exact axial stress reduction, including plastic lateral
    /// contraction; it is not a uniaxial-strain approximation.
    pub fn j2(law: UniaxialJ2, density_kg_m3: f64) -> Result<Self> {
        require(
            positive(density_kg_m3),
            "J2 fiber density must be finite and positive",
        )?;
        Ok(Self {
            law: FiberLaw::J2(law),
            density_kg_m3,
            source: None,
            chaboche_branch_properties: None,
        })
    }

    /// Resolve density and the original J2 coefficients at one fixed source
    /// state. The same hardening-basis conversion and source admission used by
    /// tetrahedral solids now supplies fiber beams and frame reinforcement.
    /// Mixed hardening requires independent isotropic and Prager data.
    /// Temperature remains fixed for this load path; no creep, thermal feedback
    /// or material-specific cyclic calibration is inferred.
    pub fn resolve_j2(
        card: &MaterialCard,
        point: &QueryPoint,
        density_property: &str,
        basis: HardeningPropertyBasis,
        strain_limit: f64,
        selection: MaterialPropertySelection,
    ) -> Result<Self> {
        let mut requirements = j2_requirements(basis)?;
        requirements.push(ScalarPropertyRequirement::try_new(
            density_property,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?);
        let source = resolve_material_state_point(card, point, &requirements, selection)?;
        let density = source
            .property(density_property)
            .expect("required density was resolved");
        if !density.is_static_at_fixed_external_conditions() {
            return Err(MaterialStatePointError::InvalidRequirement {
                property: density_property.into(),
                reason: "J2 fiber density cannot freeze sampled strain, stress, rate or history dependence",
            }.into());
        }
        let density_kg_m3 = density.value_si();
        let resolved = ResolvedJ2Plasticity::try_from_material(source, basis, strain_limit)?;
        let mut material = Self::j2(UniaxialJ2::new(resolved.law().clone())?, density_kg_m3)?;
        material.source = Some(resolved.material().clone());
        Ok(material)
    }

    /// Authored saturating metal law, sharing the native tensor J2 history and
    /// its exact free-sided axial reduction. No empirical calibration is implied.
    pub fn voce(law: UniaxialVoce, density_kg_m3: f64) -> Result<Self> {
        require(
            positive(density_kg_m3),
            "Voce fiber density must be finite and positive",
        )?;
        Ok(Self {
            law: FiberLaw::Voce(law),
            density_kg_m3,
            source: None,
            chaboche_branch_properties: None,
        })
    }

    /// Resolve independent saturation magnitude/rate, linear hardening and mass
    /// density from one fixed-condition material bundle. A post-yield tangent
    /// alone cannot determine the nonlinear hardening curve and is refused.
    pub fn resolve_voce(
        card: &MaterialCard,
        point: &QueryPoint,
        density_property: &str,
        basis: HardeningPropertyBasis,
        strain_limit: f64,
        selection: MaterialPropertySelection,
    ) -> Result<Self> {
        let mut requirements = voce_requirements(basis)?;
        requirements.push(ScalarPropertyRequirement::try_new(
            density_property,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?);
        let source = resolve_material_state_point(card, point, &requirements, selection)?;
        let density = source
            .property(density_property)
            .expect("required density was resolved");
        if !density.is_static_at_fixed_external_conditions() {
            return Err(MaterialStatePointError::InvalidRequirement {
                property: density_property.into(),
                reason: "Voce fiber density cannot freeze sampled strain, stress, rate or history dependence",
            }.into());
        }
        let density_kg_m3 = density.value_si();
        let resolved = ResolvedVoceJ2::try_from_material(source, basis, strain_limit)?;
        let mut material = Self::voce(UniaxialVoce::new(resolved.law().clone())?, density_kg_m3)?;
        material.source = Some(resolved.material().clone());
        Ok(material)
    }

    /// Authored Armstrong-Frederick/Voce metal with independent backstress
    /// recovery. Density remains an explicit physical input.
    pub fn chaboche(law: UniaxialChaboche, density_kg_m3: f64) -> Result<Self> {
        require(
            positive(density_kg_m3),
            "Chaboche fiber density must be finite and positive",
        )?;
        Ok(Self {
            law: FiberLaw::Chaboche(law),
            density_kg_m3,
            source: None,
            chaboche_branch_properties: None,
        })
    }

    /// Resolve one nonlinear backstress branch, Voce/Prager data and density
    /// from one immutable fixed-condition material bundle. C and gamma need
    /// independent source properties; monotonic tangent data cannot supply them.
    pub fn resolve_chaboche(
        card: &MaterialCard,
        point: &QueryPoint,
        density_property: &str,
        basis: HardeningPropertyBasis,
        strain_limit: f64,
        selection: MaterialPropertySelection,
    ) -> Result<Self> {
        let mut requirements = chaboche_requirements(basis)?;
        requirements.push(ScalarPropertyRequirement::try_new(
            density_property,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?);
        let source = resolve_material_state_point(card, point, &requirements, selection)?;
        let density = source
            .property(density_property)
            .expect("required density was resolved");
        if !density.is_static_at_fixed_external_conditions() {
            return Err(MaterialStatePointError::InvalidRequirement {
                property: density_property.into(),
                reason: "Chaboche fiber density cannot freeze sampled strain, stress, rate or history dependence",
            }.into());
        }
        let density_kg_m3 = density.value_si();
        let resolved = ResolvedChabocheJ2::try_from_material(source, basis, strain_limit)?;
        let mut material = Self::chaboche(
            UniaxialChaboche::new(resolved.law().clone())?,
            density_kg_m3,
        )?;
        material.source = Some(resolved.material().clone());
        Ok(material)
    }

    /// Authored one-through-eight-branch metal law, preserving the native
    /// ordered histories and exact axial restriction. No calibration is inferred.
    pub fn chaboche_branches<const N: usize>(
        law: UniaxialChaboche<N>,
        density_kg_m3: f64,
    ) -> Result<Self> {
        Self::authored(FiberLaw::chaboche_branches(law)?, density_kg_m3)
    }

    /// Resolve named recovery branches, Voce/Prager coefficients and reference
    /// density from one fixed-condition bundle. Branch order is caller supplied
    /// and retained; no unused native branches are introduced.
    pub fn resolve_chaboche_branches<const N: usize>(
        card: &MaterialCard,
        point: &QueryPoint,
        density_property: &str,
        basis: HardeningPropertyBasis,
        strain_limit: f64,
        selection: MaterialPropertySelection,
        properties: [BackstressProperties; N],
    ) -> Result<Self> {
        let mut requirements = chaboche_branch_requirements(basis, &properties)?;
        if requirements
            .iter()
            .any(|requirement| requirement.name() == density_property)
        {
            return Err(MaterialStatePointError::InvalidRequirement {
                property: density_property.into(),
                reason: "fiber density must have a distinct name from every constitutive property role",
            }.into());
        }
        requirements.push(ScalarPropertyRequirement::try_new(
            density_property,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?);
        let source = resolve_material_state_point(card, point, &requirements, selection)?;
        let density = source
            .property(density_property)
            .expect("required density was resolved");
        if !density.is_static_at_fixed_external_conditions() {
            return Err(MaterialStatePointError::InvalidRequirement {
                property: density_property.into(),
                reason: "Chaboche fiber density cannot freeze sampled strain, stress, rate or history dependence",
            }.into());
        }
        let density_kg_m3 = density.value_si();
        let resolved =
            ResolvedChabocheBranches::try_from_material(source, basis, strain_limit, properties)?;
        let mut material = Self::chaboche_branches(
            UniaxialChaboche::new(resolved.law().clone())?,
            density_kg_m3,
        )?;
        material.chaboche_branch_properties = Some(resolved.properties().to_vec());
        material.source = Some(resolved.material().clone());
        Ok(material)
    }

    /// Explicit research-model parameters, with no implied catalog authority.
    pub fn concrete(law: ManderConcrete, density_kg_m3: f64) -> Result<Self> {
        require(
            positive(density_kg_m3)
                && [law.fcc, law.eps_cc, law.ec, law.eps_cu]
                    .into_iter()
                    .all(positive),
            "concrete parameters and density must be finite and positive",
        )?;
        ManderConcrete::new(law.fcc, law.eps_cc, law.ec, law.eps_cu)?;
        require(
            (law.ec / (law.ec - law.fcc / law.eps_cc)).is_finite(),
            "Mander shape parameter must be representable",
        )?;
        Ok(Self {
            law: FiberLaw::Concrete(law),
            density_kg_m3,
            source: None,
            chaboche_branch_properties: None,
        })
    }

    /// Explicit research-model parameters, including cyclic curvature.
    pub fn steel(law: MenegottoPintoSteel, density_kg_m3: f64) -> Result<Self> {
        require(
            positive(density_kg_m3)
                && [law.e0, law.fy, law.r0, law.a2].into_iter().all(positive)
                && law.b.is_finite()
                && law.a1.is_finite()
                && law.a1 >= 0.0
                && law.r0 > law.a1
                && positive(law.fy / law.e0),
            "steel parameters need positive finite density/modulus/yield/curvature, 0 <= a1 < R0 and a2 > 0",
        )?;
        MenegottoPintoSteel::new(law.e0, law.fy, law.b)?;
        Ok(Self {
            law: FiberLaw::Steel(law),
            density_kg_m3,
            source: None,
            chaboche_branch_properties: None,
        })
    }

    /// Resolve all Mander coefficients and density from one card at one fixed
    /// state, then admit the same constitutive law used by the section engine.
    pub fn resolve_concrete(
        card: &MaterialCard,
        point: &QueryPoint,
        properties: &ConcreteFiberProperties,
        selection: MaterialPropertySelection,
    ) -> Result<Self> {
        let (source, v) = resolve(
            card,
            point,
            [
                (
                    &properties.density,
                    Density::DIMS,
                    ScalarAdmissibility::StrictlyPositive,
                ),
                (
                    &properties.compressive_strength,
                    Pressure::DIMS,
                    ScalarAdmissibility::StrictlyPositive,
                ),
                (
                    &properties.peak_strain,
                    Dims::NONE,
                    ScalarAdmissibility::StrictlyPositive,
                ),
                (
                    &properties.initial_modulus,
                    Pressure::DIMS,
                    ScalarAdmissibility::StrictlyPositive,
                ),
                (
                    &properties.ultimate_strain,
                    Dims::NONE,
                    ScalarAdmissibility::StrictlyPositive,
                ),
            ],
            selection,
        )?;
        let mut material = Self::concrete(ManderConcrete::new(v[1], v[2], v[3], v[4])?, v[0])?;
        material.source = Some(source);
        Ok(material)
    }

    /// Resolve every steel coefficient; no default hardening or Bauschinger
    /// curvature is substituted for missing source data.
    pub fn resolve_steel(
        card: &MaterialCard,
        point: &QueryPoint,
        properties: &SteelFiberProperties,
        selection: MaterialPropertySelection,
    ) -> Result<Self> {
        let (source, v) = resolve(
            card,
            point,
            [
                (
                    &properties.density,
                    Density::DIMS,
                    ScalarAdmissibility::StrictlyPositive,
                ),
                (
                    &properties.young_modulus,
                    Pressure::DIMS,
                    ScalarAdmissibility::StrictlyPositive,
                ),
                (
                    &properties.yield_stress,
                    Pressure::DIMS,
                    ScalarAdmissibility::StrictlyPositive,
                ),
                (
                    &properties.hardening_ratio,
                    Dims::NONE,
                    ScalarAdmissibility::NonNegative,
                ),
                (
                    &properties.initial_curvature,
                    Dims::NONE,
                    ScalarAdmissibility::StrictlyPositive,
                ),
                (
                    &properties.curvature_degradation,
                    Dims::NONE,
                    ScalarAdmissibility::NonNegative,
                ),
                (
                    &properties.curvature_scale,
                    Dims::NONE,
                    ScalarAdmissibility::StrictlyPositive,
                ),
            ],
            selection,
        )?;
        let mut law = MenegottoPintoSteel::new(v[1], v[2], v[3])?;
        law.r0 = v[4];
        law.a1 = v[5];
        law.a2 = v[6];
        let mut material = Self::steel(law, v[0])?;
        material.source = Some(source);
        Ok(material)
    }

    /// Selected material state and its original usage receipts, when resolved.
    #[must_use]
    pub const fn source(&self) -> Option<&ResolvedMaterialStatePoint> {
        self.source.as_ref()
    }
    /// Ordered source property roles for a named branch binding, including N=1.
    /// Authored laws and legacy single-branch bindings return None. Source names
    /// identify branch roles; the complete original receipt bundle remains in source().
    #[must_use]
    pub fn chaboche_branch_properties(&self) -> Option<&[BackstressProperties]> {
        self.chaboche_branch_properties.as_deref()
    }
    /// Immutable law coefficients; no integration-point history is stored here.
    #[must_use]
    pub const fn law(&self) -> &FiberLaw {
        &self.law
    }
    /// Reference mass density in kg/m³.
    #[must_use]
    pub const fn density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }
}

fn resolve<const N: usize>(
    card: &MaterialCard,
    point: &QueryPoint,
    roles: [(&String, Dims, ScalarAdmissibility); N],
    selection: MaterialPropertySelection,
) -> Result<(ResolvedMaterialStatePoint, [f64; N])> {
    let requirements: Vec<_> = roles
        .iter()
        .map(|(name, dims, admissibility)| {
            ScalarPropertyRequirement::try_new(name.as_str(), *dims, *admissibility)
        })
        .collect::<core::result::Result<_, _>>()?;
    let source = resolve_material_state_point(card, point, &requirements, selection)?;
    let mut values = [0.0; N];
    for (i, (name, _, _)) in roles.iter().enumerate() {
        let property = source
            .property(name)
            .expect("required property was resolved");
        if !property.is_static_at_fixed_external_conditions() {
            return Err(MaterialStatePointError::InvalidRequirement {
                property: name.as_str().to_owned(),
                reason: "RC fiber laws cannot freeze sampled strain, stress, rate, frequency or history dependence",
            }.into());
        }
        values[i] = property.value_si();
    }
    Ok((source, values))
}

/// Original rectangular geometry in SI units. Cover surrounds all four faces.
/// Concrete strips are split at every material boundary before midpoint
/// quadrature; refinement never changes constituent areas or reference mass.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RcSectionGeometry {
    /// Gross depth along the section's fiber coordinate, in metres.
    pub depth_m: f64,
    /// Gross width perpendicular to the fiber coordinate, in metres.
    pub width_m: f64,
    /// Concrete cover thickness on each of the four faces, in metres.
    pub cover_m: f64,
    /// Uniform depth subdivisions, before splitting at material boundaries.
    /// Bounded to 1..=4096.
    pub concrete_layers: usize,
    /// Reinforcement area at each of the two symmetric faces, in m².
    pub steel_area_per_face_m2: f64,
    /// Positive distance of each symmetric steel layer from the section center.
    pub steel_centroid_y_m: f64,
}

/// An immutable virgin section, its constituent bindings, and its coherent mass.
/// Runtime constitutive history lives in a clone of the section, never in cards.
#[derive(Clone, Debug)]
pub struct RcSection {
    geometry: RcSectionGeometry,
    materials: [FiberMaterial; 3],
    section: Section,
    line_mass_kg_m: f64,
    area_scale: f64,
}
impl RcSection {
    /// Build core, four-sided cover and symmetric reinforcement without counting
    /// the displaced concrete in either constitutive response or member mass.
    #[allow(clippy::cast_precision_loss)]
    pub fn new(
        geometry: RcSectionGeometry,
        core: FiberMaterial,
        cover: FiberMaterial,
        steel: FiberMaterial,
    ) -> Result<Self> {
        let g = geometry;
        require(
            [
                g.depth_m,
                g.width_m,
                g.cover_m,
                g.steel_area_per_face_m2,
                g.steel_centroid_y_m,
            ]
            .into_iter()
            .all(positive)
                && (1..=4096).contains(&g.concrete_layers),
            "positive finite dimensions/reinforcement and 1..=4096 concrete layers are required",
        )?;
        require(
            matches!(core.law, FiberLaw::Concrete(_))
                && matches!(cover.law, FiberLaw::Concrete(_))
                && matches!(
                    steel.law,
                    FiberLaw::Steel(_)
                        | FiberLaw::J2(_)
                        | FiberLaw::Voce(_)
                        | FiberLaw::Chaboche(_)
                        | FiberLaw::ChabocheBranches(_)
                ),
            "core and cover require concrete; reinforcement requires a steel or native J2 law",
        )?;
        let core_width = g.width_m - 2.0 * g.cover_m;
        let core_half_depth = 0.5 * g.depth_m - g.cover_m;
        require(
            positive(core_width) && positive(core_half_depth) && positive(g.depth_m * g.width_m),
            "cover must leave a positive representable core and gross area",
        )?;
        let half_band = 0.5 * g.steel_area_per_face_m2 / core_width;
        let inner = g.steel_centroid_y_m - half_band;
        let outer = g.steel_centroid_y_m + half_band;
        require(
            positive(half_band)
                && positive(inner)
                && inner < g.steel_centroid_y_m
                && g.steel_centroid_y_m < outer
                && outer < core_half_depth,
            "reinforcement bands must be representable, disjoint and fit strictly inside the core",
        )?;
        let half_depth = 0.5 * g.depth_m;
        let step = g.depth_m / g.concrete_layers as f64;
        require(positive(step), "concrete strip depth must be representable")?;
        let mut cuts: Vec<f64> = (0..g.concrete_layers)
            .map(|i| -half_depth + i as f64 * step)
            .collect();
        cuts.extend([
            half_depth,
            -core_half_depth,
            -outer,
            -inner,
            inner,
            outer,
            core_half_depth,
        ]);
        cuts.sort_by(f64::total_cmp);
        cuts.dedup();
        let materials = [core, cover, steel];
        let mut fibers = Vec::new();
        let mut line_mass = 0.0;
        let mut add = |y: f64, area: f64, material: &FiberMaterial| -> Result<()> {
            require(
                y.is_finite() && positive(area),
                "fiber position and area must be representable",
            )?;
            let mass = area * material.density_kg_m3;
            require(positive(mass), "fiber line mass must be representable")?;
            line_mass += mass;
            fibers.push(Fiber {
                y,
                area,
                law: material.law.clone(),
            });
            Ok(())
        };
        for interval in cuts.windows(2) {
            let (lo, hi) = (interval[0], interval[1]);
            let y = 0.5 * lo + 0.5 * hi;
            let thickness = hi - lo;
            // Cuts include every material interface. Classify the interval,
            // since a rounded midpoint can coincide with an adjacent interface.
            if hi <= -core_half_depth || lo >= core_half_depth {
                add(y, g.width_m * thickness, &materials[1])?;
            } else {
                add(y, 2.0 * g.cover_m * thickness, &materials[1])?;
                let in_steel_band = (lo >= inner && hi <= outer) || (lo >= -outer && hi <= -inner);
                if !in_steel_band {
                    add(y, core_width * thickness, &materials[0])?;
                }
            }
        }
        for y in [-g.steel_centroid_y_m, g.steel_centroid_y_m] {
            add(y, g.steel_area_per_face_m2, &materials[2])?;
        }
        require(
            positive(line_mass),
            "section line mass must be representable",
        )?;
        Ok(Self {
            geometry,
            materials,
            section: Section::new(fibers),
            line_mass_kg_m: line_mass,
            area_scale: 1.0,
        })
    }

    /// Multiply tributary fiber areas and their mass together. This is an
    /// explicit area-weight scaling, not a geometric depth/width transformation.
    /// Original geometry and all source bindings remain available.
    pub fn scaled(&self, factor: f64) -> Result<Self> {
        require(positive(factor), "area scale must be finite and positive")?;
        let mut scaled = self.clone();
        scaled.line_mass_kg_m *= factor;
        scaled.area_scale *= factor;
        require(
            positive(scaled.line_mass_kg_m) && positive(scaled.area_scale),
            "scaled mass and area must be representable",
        )?;
        for fiber in &mut scaled.section.fibers {
            fiber.area *= factor;
            require(
                positive(fiber.area),
                "scaled fiber area must be representable",
            )?;
        }
        Ok(scaled)
    }

    /// Original geometry, before any tributary area scaling.
    #[must_use]
    pub const fn geometry(&self) -> RcSectionGeometry {
        self.geometry
    }
    /// Virgin section blueprint; clone it for an independent response history.
    #[must_use]
    pub const fn section(&self) -> &Section {
        &self.section
    }
    /// Core, cover and reinforcement, in that order.
    #[must_use]
    pub const fn materials(&self) -> &[FiberMaterial; 3] {
        &self.materials
    }
    /// Constituent-derived reference line mass in kg/m, including area scaling.
    #[must_use]
    pub const fn line_mass_kg_m(&self) -> f64 {
        self.line_mass_kg_m
    }
    /// Accumulated tributary area multiplier, independent of original geometry.
    #[must_use]
    pub const fn area_scale(&self) -> f64 {
        self.area_scale
    }
}
