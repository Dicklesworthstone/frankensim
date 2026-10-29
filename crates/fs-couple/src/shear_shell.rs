//! Regional materials on the existing curved shear-deformable shell owner.
//!
//! Original density, directional elasticity, transverse shear, thickness and
//! full Kelvin--Voigt loss travel together into one immutable element assembly.
//! A source that states `constant_loss_reference_omega` has its loss factors
//! read at every retained mode's own frequency instead ([`constant_loss`]).
//! Substrate materials are homogeneous through their own thickness, at a frozen
//! source state. Explicit source-backed surface films use the original laminate
//! and MITC owners, retaining extension/bending and mass coupling. Geometry,
//! support and numerical resolution remain independent inputs.
//! No body frequency, damping ratio, carved-instrument preset or fitting to sound
//! is introduced. Prestress evolution and finite strain are absent.

use crate::acoustic_realize::AcousticRealizeError;
use crate::thin_plate::{PlateMaterialModel, PlateRegionMaterial, PlateThicknessConstraint,
    plate_elasticity, plate_property};
use fs_exec::CancelGate;
use fs_material::plane_stress_viscosity::{MaterialPlane, PlaneStressViscosity,
    PLANE_STRESS_VISCOSITY_PROPERTIES};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_plate::shell::{ShellMesh, mitc::{MitcBudget, ShearShellSection,
    viscosity::ShearShellViscosity}};
use fs_qty::{Density, DynViscosity, Pressure};

/// Ordered material-interface compliance on physical shell attachments.
pub mod connections;
/// Source-bound component surfaces joined by the original shell mortar owner.
pub mod assembly;
/// Independently owned exterior and cavity pressure boundaries.
pub mod pressure;
mod reduction;
pub use reduction::{
    MAX_SHELL_MODES, MaterialShellImage, ShellProjectionBudget, ShellReductionOptions,
    shell_mode_budget_admits,
};

/// Explicit 13,23 and cross 13/23 shear viscosities [Pa s], in material axes.
/// Missing values do not mean isotropy or zero.
pub const TRANSVERSE_SHEAR_VISCOSITY_PROPERTIES: [&str; 3] = [
    "transverse_shear_viscosity_13", "transverse_shear_viscosity_23",
    "transverse_shear_viscosity_13_23",
];
/// Optional angular frequency [rad/s] at which a source wrote its nine viscosities as `eta x stiffness / omega`: its loss
/// factors are then constant in frequency, and the modal reductions apply them at each mode's own frequency
/// ([`constant_loss`]). A region without it is Kelvin-Voigt, exactly as before.
pub const CONSTANT_LOSS_REFERENCE_OMEGA: &str = "constant_loss_reference_omega";
/// The anchor's physical dimension, 1/s.
pub const CONSTANT_LOSS_REFERENCE_OMEGA_DIMS: fs_qty::Dims = fs_qty::Dims([0, 0, -1, 0, 0, 0]);
/// The constant-loss-factor (generalized Maxwell) reading of anchored shell viscosities.
pub mod constant_loss;
/// Perfectly bonded source-backed films about the substrate reference surface.
pub mod coating;
pub use coating::{CoatingSide, ResolvedShellCoating, ShellCoating};

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}
fn owner(error: impl core::fmt::Display) -> AcousticRealizeError {
    AcousticRealizeError::Nonlinear(error.to_string())
}
fn poll(gate: &CancelGate) -> Result<(), AcousticRealizeError> {
    if gate.is_requested() { Err(AcousticRealizeError::Cancelled) } else { Ok(()) }
}

/// Same regional source, material axes and mass policy as the plate compiler,
/// with an explicit transverse-shear correction and GLOBAL material axis.
/// The region angle is an additional rotation around each facet normal.
#[derive(Clone, Debug)]
pub struct ShellRegionMaterial {
    /// Exact source object and disjoint triangle assignment.
    pub region: PlateRegionMaterial,
    /// Positive transverse-shear correction, applied to stiffness AND viscosity.
    pub shear_correction: f64,
    /// First selected section axis projected into each facet (material axis 2
    /// for Orthotropic21); no grain direction inferred from triangle edges.
    pub axis1_global: [f64;3],
}
/// Immutable source and resulting section for one material region.
#[derive(Clone, Debug)]
pub struct ResolvedShellRegion {
    input: ShellRegionMaterial,
    area_m2: f64,
    thickness_m: f64,
    mass_kg: f64,
    section: ShearShellSection,
    viscosity: ShearShellViscosity,
    plane_viscosity: PlaneStressViscosity,
    shear_pa_s: [f64;3],
    omega_band: (f64, f64),
    loss_reference_omega: Option<f64>,
}
impl ResolvedShellRegion {
    /// Original material state, assignment and comparison policy.
    pub const fn input(&self) -> &ShellRegionMaterial { &self.input }
    /// Actual midsurface area [m²].
    pub const fn area_m2(&self) -> f64 { self.area_m2 }
    /// Substrate thickness [m], independent of arching height and any coating.
    /// Attachment geometry remains relative to the original substrate skins.
    pub const fn thickness_m(&self) -> f64 { self.thickness_m }
    /// Actual mass [kg], including this region only.
    pub const fn mass_kg(&self) -> f64 { self.mass_kg }
    /// Principal-frame shear-shell elastic section, before geometric rotation.
    pub const fn section(&self) -> &ShearShellSection { &self.section }
    /// Substrate-only principal-frame viscosity, before facet rotation.
    /// Coating loss remains separately owned and is never wood-law rescaled.
    pub const fn viscosity(&self) -> &ShearShellViscosity { &self.viscosity }
    /// The original six source-owned reduced plane-stress coefficients.
    pub const fn plane_viscosity(&self) -> &PlaneStressViscosity { &self.plane_viscosity }
    /// Intersected elastic, inertial and viscous support [rad/s].
    pub const fn omega_band_rad_s(&self) -> (f64, f64) { self.omega_band }
    /// The source's [`CONSTANT_LOSS_REFERENCE_OMEGA`] [rad/s], or `None` for Kelvin-Voigt loss.
    pub const fn constant_loss_reference_omega_rad_s(&self) -> Option<f64> { self.loss_reference_omega }
}
/// Bound arbitrary triangular midsurface and every local material section.
/// Mutation requires recompilation; material records cannot drift from operators.
#[derive(Clone, Debug)]
pub struct ResolvedShearShell {
    mesh: ShellMesh,
    regions: Vec<ResolvedShellRegion>,
    element_regions: Vec<usize>,
    sections: Vec<ShearShellSection>,
    viscosity: Vec<ShearShellViscosity>,
    coatings: Vec<ResolvedShellCoating>,
    coating_viscosity: Vec<Vec<ShearShellViscosity>>,
    budget: MitcBudget,
    mass_kg: f64,
    interfaces: Vec<connections::ResolvedShellConnection>,
    interface_budget: fs_plate::shell::mitc::connections::ShellConnectionBudget,
    pressure_surfaces: Option<pressure::ResolvedPressureSurfaces>,
}
impl ResolvedShearShell {
    /// Original finite, oriented physical midsurface.
    pub const fn mesh(&self) -> &ShellMesh { &self.mesh }
    /// Original regional material records and section results.
    pub fn regions(&self) -> &[ResolvedShellRegion] { &self.regions }
    /// Total physical mass integrated over actual facet area [kg].
    pub const fn mass_kg(&self) -> f64 { self.mass_kg }
    /// Original material assignment for an addressed element.
    pub fn material_at_element(&self, element: usize) -> Option<&ResolvedMaterialStatePoint> {
        self.element_regions.get(element).map(|&r| &self.regions[r].input.region.material)
    }
    /// One exact material section per triangle, in original mesh order.
    pub fn sections(&self) -> &[ShearShellSection] { &self.sections }
    /// Original substrate viscosities, in facet coordinates. Coating groups
    /// remain separate because their measured frequency laws differ.
    pub fn substrate_viscosities(&self) -> &[ShearShellViscosity] { &self.viscosity }
    /// One film's complete facet viscosities at its declared storage reference.
    /// Use `project_material_viscosity` to read its loss at actual modal frequencies.
    pub fn coating_viscosities(&self, coating: usize) -> Option<&[ShearShellViscosity]> {
        self.coating_viscosity.get(coating).map(Vec::as_slice)
    }
    /// The substrate viscosities grouped by material loss law, each with its anchor ([`CONSTANT_LOSS_REFERENCE_OMEGA`],
    /// `None` for Kelvin-Voigt). One law borrows the original viscosities unchanged; otherwise each group zeroes the
    /// elements of the other laws, so the groups' projections sum to the single one.
    pub(crate) fn loss_groups(&self) -> Result<Vec<(Option<f64>, std::borrow::Cow<'_, [ShearShellViscosity]>)>,
        AcousticRealizeError>
    {
        let laws: Vec<Option<u64>> = self.regions.iter().map(|r| r.loss_reference_omega.map(f64::to_bits)).collect();
        let mut distinct = laws.clone();
        distinct.sort_unstable();
        distinct.dedup();
        if let [law] = distinct[..] {
            return Ok(vec![(law.map(f64::from_bits), std::borrow::Cow::Borrowed(&self.viscosity[..]))]);
        }
        let zero = ShearShellViscosity::new([0.0; 9], [0.0; 9], [0.0; 3]).map_err(owner)?;
        Ok(distinct.into_iter().map(|law| {
            let only = self.element_regions.iter().zip(&self.viscosity)
                .map(|(&r, v)| if laws[r] == law { *v } else { zero }).collect();
            (law.map(f64::from_bits), std::borrow::Cow::Owned(only))
        }).collect())
    }
    /// Swap one region's complete material data; original geometry, mass policy,
    /// orientation and correction survive. Failed rebind leaves self unchanged.
    pub fn with_region_material(&self, name: &str, material: &ResolvedMaterialStatePoint,
        gate: &CancelGate) -> Result<Self, AcousticRealizeError>
    {
        let index = self.regions.iter().position(|r| r.input.region.region.name == name)
            .ok_or_else(|| invalid("unknown shell material region"))?;
        let mut inputs: Vec<_> = self.regions.iter().map(|r| r.input.clone()).collect();
        inputs[index].region.material = material.clone();
        let rebuilt = compile_material_shear_shell(self.mesh.clone(), inputs, self.budget, gate)?;
        let rebuilt = coating::apply_coatings(rebuilt,
            self.coatings.iter().map(|c| c.input().clone()).collect(), gate)?;
        let mut rebuilt = rebuilt.with_interfaces(self.interfaces.iter().map(|i| i.input().clone()).collect(),
            self.interface_budget, gate)?;
        // Pressure geometry is unchanged by a material rebind. Retain both the
        // authored selection and its original work-conjugate physical rows.
        rebuilt.pressure_surfaces = self.pressure_surfaces.clone();
        Ok(rebuilt)
    }
}
fn restrict_band(source: &ResolvedMaterialStatePoint, key: &str,
    band: &mut (f64,f64), require_frequency: bool) -> Result<(), AcousticRealizeError>
{
    let property = source.property(key).ok_or_else(|| invalid("missing shell material coefficient"))?;
    if !property.is_constant_at_fixed_temperature() {
        return Err(invalid("shell coefficients must be constants at the original frozen temperature"));
    }
    match property.angular_frequency_band_rad_s().map_err(owner)? {
        Some((lo,hi)) => { band.0=band.0.max(lo); band.1=band.1.min(hi); },
        None if require_frequency => return Err(invalid("each shell viscosity requires an explicit frequency band")),
        None => {},
    }
    Ok(())
}

/// Compile complete disjoint regional assignments. Isotropic elasticity derives
/// G from E and nu; orthotropic sources additionally require G13 and G23.
/// In-plane viscosity uses the shared reduced plane-stress owner. The three
/// transverse shear viscosities are independently sourced, including zeros.
/// Neither elastic anisotropy nor loss is inferred from a material's name.
///
/// Every selected coefficient is constant at its source temperature and every
/// viscosity has frequency applicability. Actual retained frequencies are
/// checked by reduction. This is unloaded linear geometry, not a prestress solve.
/// # Errors
/// Incomplete coverage, invalid material/units/geometry, excessive dimensions,
/// nonpositive mass/stiffness/viscosity, nonintersecting bands or cancellation.
pub fn compile_material_shear_shell(mesh: ShellMesh, inputs: Vec<ShellRegionMaterial>,
    budget: MitcBudget, gate: &CancelGate) -> Result<ResolvedShearShell, AcousticRealizeError>
{
    poll(gate)?;
    if mesh.nodes.len() > budget.max_nodes || mesh.tris.len() > budget.max_triangles
        || inputs.is_empty() || inputs.len() > 1024 {
        return Err(invalid("material shell exceeds its explicit mesh/region budget"));
    }
    fs_plate::shell::survey::validate_midsurface(&mesh, fs_plate::shell::profile::ProfileBudget {
        max_nodes:budget.max_nodes, max_triangles:budget.max_triangles, max_feature_evaluations:0,
    }).map_err(owner)?;
    let mut assigned=vec![usize::MAX;mesh.tris.len()];
    let mut names=std::collections::BTreeSet::new();
    for (r,input) in inputs.iter().enumerate() {
        poll(gate)?;
        let region=&input.region.region;
        if region.name.trim().is_empty() || !names.insert(region.name.as_str())
            || region.triangle_indices.is_empty() {
            return Err(invalid("shell regions need unique names and nonempty assignments"));
        }
        for &element in &region.triangle_indices {
            let slot=assigned.get_mut(element).ok_or_else(|| invalid("shell region addresses unknown triangle"))?;
            if *slot!=usize::MAX { return Err(invalid("shell material regions overlap")); }
            *slot=r;
        }
    }
    if assigned.contains(&usize::MAX) { return Err(invalid("shell has an unassigned material triangle")); }
    let mut areas=vec![0.0;inputs.len()];
    let mut element_areas=Vec::with_capacity(mesh.tris.len());
    for (e,&region) in assigned.iter().enumerate() {
        poll(gate)?;
        let area=mesh.facet(e).map_err(owner)?.area_m2;
        areas[region]+=area; element_areas.push(area);
    }
    let mut regions=Vec::with_capacity(inputs.len());
    for (input,area_m2) in inputs.into_iter().zip(areas) {
        poll(gate)?;
        let r=&input.region; let state=&r.material;
        let rho=plate_property(state,"density",Density::DIMS)?;
        let (e1,e2,nu,g12)=plate_elasticity(state,r.model)?;
        let (g13,g23)=if r.model==PlateMaterialModel::Isotropic { (g12,g12) } else {
            let a=plate_property(state,"shear_modulus_13",Pressure::DIMS)?;
            let b=plate_property(state,"shear_modulus_23",Pressure::DIMS)?;
            if r.model==PlateMaterialModel::Orthotropic21 { (b,a) } else { (a,b) }
        };
        let h=match r.thickness_constraint {
            PlateThicknessConstraint::FixedThickness(h)=>h,
            PlateThicknessConstraint::FixedMass(m)=>m/area_m2/rho,
        };
        let section=ShearShellSection::from_orthotropic_solid(e1,e2,nu,g12,g13,g23,
            h,rho,0.0,input.shear_correction).map_err(owner)?;
        let keys=PLANE_STRESS_VISCOSITY_PROPERTIES;
        let (plane,keys)=if r.model==PlateMaterialModel::Orthotropic21 {
            (MaterialPlane::Axes21,[keys[1],keys[0],keys[2],keys[3],keys[5],keys[4]])
        } else { (MaterialPlane::Axes12,keys) };
        let plane_viscosity=PlaneStressViscosity::from_material(state.clone(),plane,keys).map_err(owner)?;
        let mut band=plane_viscosity.omega_band_rad_s();
        let elastic: &[&str]=if r.model==PlateMaterialModel::Isotropic {
            &["density","young_modulus","poisson_ratio"]
        } else {
            &["density","young_modulus_1","young_modulus_2","poisson_ratio_12",
                "shear_modulus_12","shear_modulus_13","shear_modulus_23"]
        };
        for key in elastic { restrict_band(state,key,&mut band,false)?; }
        let mut shear=[0.0;3];
        for (i,key) in TRANSVERSE_SHEAR_VISCOSITY_PROPERTIES.iter().enumerate() {
            shear[i]=plate_property(state,key,DynViscosity::DIMS)?;
            restrict_band(state,key,&mut band,true)?;
        }
        // The anchor is read only when the source states it; its law is claimed only over the spectrum's fit band.
        let loss_reference_omega = if state.property(CONSTANT_LOSS_REFERENCE_OMEGA).is_some() {
            let w = plate_property(state, CONSTANT_LOSS_REFERENCE_OMEGA, CONSTANT_LOSS_REFERENCE_OMEGA_DIMS)?;
            if w <= 0.0 { return Err(invalid("a constant-loss reference angular frequency must be positive")); }
            restrict_band(state, CONSTANT_LOSS_REFERENCE_OMEGA, &mut band, false)?;
            let (lo, hi) = constant_loss::SPECTRUM_BAND_RAD_S;
            band.0 = band.0.max(lo);
            band.1 = band.1.min(hi);
            Some(w)
        } else { None };
        if !band.0.is_finite() || !band.1.is_finite() || band.0<0.0 || band.1<band.0 {
            return Err(invalid("shell material coefficient applicability bands do not intersect"));
        }
        if r.model==PlateMaterialModel::Orthotropic21 { shear.swap(0,1); }
        let viscosity=viscosity_in_frame(&plane_viscosity,shear,h,input.shear_correction,0.0)?;
        let mass_kg=rho*h*area_m2;
        if !area_m2.is_finite() || area_m2<=0.0 || !mass_kg.is_finite() || mass_kg<=0.0 {
            return Err(invalid("shell region area or mass is unrepresentable"));
        }
        regions.push(ResolvedShellRegion { input,area_m2,thickness_m:h,mass_kg,
            section,viscosity,plane_viscosity,shear_pa_s:shear,omega_band:band,loss_reference_omega });
    }
    let mut sections=Vec::with_capacity(assigned.len());
    let mut viscosity=Vec::with_capacity(assigned.len());
    for (e,&r) in assigned.iter().enumerate() {
        poll(gate)?;
        let region=&regions[r];
        let angle=fs_plate::shell::survey::material_axis_angle(
            &mesh.facet(e).map_err(owner)?,region.input.axis1_global).map_err(owner)?
            +region.input.region.material_angle_rad;
        sections.push(region.section.with_material_angle_rad(angle).map_err(owner)?);
        viscosity.push(viscosity_in_frame(&region.plane_viscosity,region.shear_pa_s,
            region.thickness_m,region.input.shear_correction,angle)?);
    }
    let mass_kg:f64=sections.iter().zip(element_areas).map(|(s,a)|s.areal_mass_kg_m2()*a).sum();
    if !mass_kg.is_finite() || mass_kg<=0.0 { return Err(invalid("total shell mass is unrepresentable")); }
    Ok(ResolvedShearShell { mesh,regions,element_regions:assigned,sections,viscosity,
        coatings:Vec::new(),coating_viscosity:Vec::new(),budget,mass_kg,
        interfaces:Vec::new(), interface_budget:fs_plate::shell::mitc::connections::ShellConnectionBudget {
            max_connections:0, max_terms:0, max_work:0,
        }, pressure_surfaces:None })
}

fn viscosity_in_frame(law: &PlaneStressViscosity, shear: [f64;3],
    h: f64, correction: f64, angle: f64) -> Result<ShearShellViscosity,AcousticRealizeError>
{
    // Same R^T V R convention as the original elastic shear section.
    let (s,c)=(fs_math::det::sin(angle),fs_math::det::cos(angle));
    let [a,b,cross]=shear;
    let scale=h*correction;
    let rotated=[scale*(a*c*c-2.0*cross*c*s+b*s*s),
        scale*((a-b)*c*s+cross*(c*c-s*s)),
        scale*(a*s*s+2.0*cross*c*s+b*c*c)];
    let v=law.rotated_matrix_pa_s(angle).map_err(owner)?;
    ShearShellViscosity::new(v.map(|x|x*h),
        law.bending_matrix_n_m_s(h,angle).map_err(owner)?,rotated).map_err(owner)
}
