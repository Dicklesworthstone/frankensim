//! Material temperatures -> constrained FEM displacement/stress -> conserved mass.
//!
//! This composes existing owners. The material database selects E(T), nu(T)
//! and the original integrated alpha path. `fs-solid::linear3` solves compatible
//! displacement and recovers support reactions/stress. The mesh mass adapter
//! transports the same reference matter through those actual nodes, not free
//! isotropic scaling. Temperatures are prescribed per cell, not solved here.

mod vibration;
mod loads;

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::state_point::{
    DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY, POISSON_RATIO_PROPERTY,
    IntegratedIsotropicThermalExpansion, MaterialPropertySelection, MaterialStatePointError,
    ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement,
    integrate_isotropic_thermal_expansion, resolve_material_state_point,
};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_solid::linear3::{TetAssemblyBudget, TetDeformedMesh, TetElasticError, TetElasticMaterial,
    TetLinearElasticProblem, TetMaterialField, TetStaticSolveConfig,
    TetThermalDisplacementSolution, TetThermalStrainField, TetThermalStrainState, TetThermalStressReport};
use fs_solid::loads3::LoadError3;
use super::{DeformedBody, TetrahedralMass};
use crate::MaterialInertiaError;

/// One cell's actual temperature and unchanged source card. All cells are required.
#[derive(Clone, Debug)]
pub struct ThermalCell<'a> {
    /// Exact reference region key, independent of solver element order.
    pub key: &'a str,
    /// Same immutable card that supplied this region's original density.
    pub card: &'a MaterialCard,
    /// Current positive absolute temperature [K].
    pub temperature_k: f64,
    /// Selection of E and nu at the current physical state.
    pub stiffness_selection: MaterialPropertySelection,
    /// Selection of alpha at the original reference state, pinned along its path.
    pub expansion_selection: MaterialPropertySelection,
}

/// Explicit geometry, solver and small-strain applicability limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalSolveConfig {
    /// Existing element/vertex/DOF/quality envelope.
    pub assembly: TetAssemblyBudget,
    /// Existing true-residual PCG controls.
    pub solve: TetStaticSolveConfig,
    /// Caller-declared coordinate-frame identity for physical stress.
    pub mesh_frame: ContentHash,
    /// Bound on free, total and elastic Mandel-strain components.
    pub maximum_absolute_mandel_strain: f64,
    /// Physical maximum nodal displacement [m], checked before geometry publication.
    pub maximum_nodal_displacement_m: f64,
}

/// Original source or physical-solver error; no partly configured body escapes.
#[derive(Clone, Debug, PartialEq)]
pub enum ThermalBodyError {
    /// Invalid material/geometry binding or unrepresentable mass transport.
    Mass(MaterialInertiaError),
    /// Original material resolution/path refusal.
    Material(MaterialStatePointError),
    /// Original solid assembly, equilibrium, stress or geometry refusal.
    Solid(TetElasticError),
    /// Original reference-geometry load integration refusal.
    Loads(LoadError3),
}
impl core::fmt::Display for ThermalBodyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self { Self::Mass(e)=>write!(f,"thermal body mass: {e}"),
            Self::Material(e)=>write!(f,"thermal body material: {e}"),
            Self::Solid(e)=>write!(f,"thermal body solid: {e}"),
            Self::Loads(e)=>write!(f,"thermal body load: {e}") }
    }
}
impl std::error::Error for ThermalBodyError {
    fn source(&self)->Option<&(dyn std::error::Error+'static)> {
        match self {Self::Mass(e)=>Some(e),Self::Material(e)=>Some(e),Self::Solid(e)=>Some(e),Self::Loads(e)=>Some(e)}
    }
}
impl From<MaterialInertiaError> for ThermalBodyError {fn from(e:MaterialInertiaError)->Self{Self::Mass(e)}}
impl From<MaterialStatePointError> for ThermalBodyError {fn from(e:MaterialStatePointError)->Self{Self::Material(e)}}
impl From<TetElasticError> for ThermalBodyError {fn from(e:TetElasticError)->Self{Self::Solid(e)}}
impl From<LoadError3> for ThermalBodyError {fn from(e:LoadError3)->Self{Self::Loads(e)}}

/// Actual source records used by one solver cell; retained in reference region order.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalCellSource {
    stiffness: ResolvedMaterialStatePoint,
    expansion: IntegratedIsotropicThermalExpansion,
    material_state_identity: ContentHash,
}
impl ThermalCellSource {
    /// Current E/nu data and full query/claim records. No current density was requested.
    pub const fn stiffness(&self)->&ResolvedMaterialStatePoint {&self.stiffness}
    /// Original source-integrated stress-free strain, not a displacement field.
    pub const fn expansion(&self)->&IntegratedIsotropicThermalExpansion {&self.expansion}
    /// Identity binding current stiffness to the original reference-density source.
    pub const fn material_state_identity(&self)->ContentHash {self.material_state_identity}
}

/// One accepted physical solution, stress field, mesh and material mass configuration.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalBody<'a> {
    body: DeformedBody<'a>,
    displacement: TetThermalDisplacementSolution,
    stress: TetThermalStressReport,
    support_reactions_n: Vec<[f64; 3]>,
    geometry: TetDeformedMesh,
    sources: Vec<ThermalCellSource>,
    // Exact accepted operator inputs, never reconstructed from material names.
    linearization: vibration::Linearization,
}
impl<'a> ThermalBody<'a> {
    /// Material inventory and full inertia from the ACTUAL constrained displacement.
    pub const fn body(&self)->&DeformedBody<'a> {&self.body}
    /// Original equilibrium result, including independently checked residual.
    pub const fn displacement(&self)->&TetThermalDisplacementSolution {&self.displacement}
    /// Compatible elastic stress, energy and internal-force field on reference cells.
    /// With applied mechanical loads, internal force alone is not support reaction.
    pub const fn stress(&self)->&TetThermalStressReport {&self.stress}
    /// Actual support force: internal minus applied force on fixed DOFs; zero on
    /// free DOFs. Applied loads at supports are retained in this balance.
    pub fn support_reactions_n(&self)->&[[f64; 3]] {&self.support_reactions_n}
    /// Exact accepted dead forces in sorted full node-major xyz DOF order [N].
    pub fn nodal_forces_n(&self)->&[(usize, f64)] {&self.linearization.nodal_forces_n}
    /// Existing admitted geometric update and its exact solution-bound identity.
    pub const fn geometry(&self)->&TetDeformedMesh {&self.geometry}
    /// Current constitutive sources in canonical reference-region order.
    pub fn sources(&self)->&[ThermalCellSource] {&self.sources}
}

impl<'a> TetrahedralMass<'a> {
    /// Solve nonuniform small-strain isotropic thermoelasticity and transport mass.
    ///
    /// Fixed DOFs are the existing solid owner's zero-displacement supports.
    /// Each cell's full reference query is restored; only absolute T changes.
    /// E and nu are queried at current T and alpha is integrated from its own
    /// original T. Reference density supplies the Lagrangian FEM mass measure:
    /// querying rho(current T) on unchanged reference volume would replace matter.
    /// No current density, yield, damping or unrelated material data are required.
    ///
    /// Shared FEM nodes impose compatibility; differently heated bonded cells do
    /// not get independent geometric scales. Stress, reactions, geometry quality,
    /// displacement and strain bounds must all accept before returning the body.
    /// This is static linear elasticity on a fixed reference mesh. Its exact
    /// operator inputs are retained for `ThermalBody::initial_stress_pencil`.
    /// It is not a heat solve, plasticity, finite-strain dynamics or a
    /// momentum/energy-conserving transition of an already moving body.
    pub fn solve_thermal(
        &self, cells:&[ThermalCell<'_>], fixed_dofs:&[usize],
        config:ThermalSolveConfig, cx:&Cx<'_>,
    )->Result<ThermalBody<'a>,ThermalBodyError> {
        self.solve_thermal_with_nodal_forces(cells, fixed_dofs, &[], config, cx)
    }

    /// Solve the same material state with concurrent mechanical dead loads.
    ///
    /// `nodal_forces_n` uses the original solid owner's sorted, unique full xyz
    /// DOFs and newtons. It can consume integrated pressure/traction loads from
    /// `fs_solid::loads3` or [`Self::gravity_nodal_forces`]. Forces are fixed in
    /// the reference frame; follower loads require a different tangent.
    /// The exact force vector remains attached to this snapshot's stress and
    /// initial-stress pencil. Empty forces preserve the thermal-only identity.
    /// No mechanical work is inserted into a heat balance by this static solve.
    pub fn solve_thermal_with_nodal_forces(
        &self, cells:&[ThermalCell<'_>], fixed_dofs:&[usize],
        nodal_forces_n:&[(usize, f64)], config:ThermalSolveConfig, cx:&Cx<'_>,
    )->Result<ThermalBody<'a>,ThermalBodyError> {
        checkpoint(cx)?;
        if cells.len()!=self.reference.regions().len()
            || self.nodes_m.len()>config.assembly.maximum_nodes
            || self.tetrahedra.len()>config.assembly.maximum_tetrahedra
            || nodal_forces_n.len()>self.nodes_m.len().saturating_mul(3) {
            return Err(invalid("complete cell states and adequate original solid budgets required"));
        }
        if !config.maximum_absolute_mandel_strain.is_finite() || config.maximum_absolute_mandel_strain<=0.0
            || !config.maximum_nodal_displacement_m.is_finite() || config.maximum_nodal_displacement_m<=0.0 {
            return Err(invalid("positive finite strain and displacement limits required"));
        }
        let mut ordered:Vec<_>=cells.iter().collect();
        ordered.sort_by_key(|cell|cell.key);
        let requirements=[
            ScalarPropertyRequirement::try_new(YOUNG_MODULUS_PROPERTY,Pressure::DIMS,ScalarAdmissibility::StrictlyPositive)?,
            ScalarPropertyRequirement::try_new(POISSON_RATIO_PROPERTY,Dims([0;6]),
                ScalarAdmissibility::OpenInterval {lower:-1.0,upper:0.5})?,
        ];
        let mut sources=Vec::with_capacity(cells.len());
        let mut operators=Vec::with_capacity(cells.len());
        for (index,(original,cell)) in self.reference.regions().iter().zip(ordered).enumerate() {
            checkpoint(cx)?;
            if cell.key!=original.key() || cell.card.content_hash()!=original.material().card_identity()
                || cell.card.id()!=original.material().material() {
                return Err(invalid("cell keys and material cards must match the complete reference inventory"));
            }
            let reference=reference_point(original.material())?;
            if !cell.temperature_k.is_finite() || cell.temperature_k<=0.0 {
                return Err(invalid("cell temperature must be positive absolute kelvin"));
            }
            let current=reference.clone().with("T",cell.temperature_k)
                .map_err(|_|invalid("current temperature does not match the reference query schema"))?;
            let expansion=integrate_isotropic_thermal_expansion(cell.card,&reference,&current,cell.expansion_selection.clone())?;
            let stiffness=resolve_material_state_point(cell.card,&current,&requirements,cell.stiffness_selection.clone())?;
            let mut hash=DomainHasher::new("org.frankensim.fs-mbd-material.reference-mass-thermal-elastic.v1");
            hash.update(original.material().identity().as_bytes());
            hash.update(stiffness.identity().as_bytes());
            let material_state_identity=hash.finalize();
            let value=|key|stiffness.property(key).expect("required stiffness property resolved").value_si();
            let density=original.material().property(DENSITY_PROPERTY).expect("reference density resolved").value_si();
            let material=TetElasticMaterial::try_new(density,value(YOUNG_MODULUS_PROPERTY),value(POISSON_RATIO_PROPERTY),
                *material_state_identity.as_bytes())?;
            let strain=expansion.free_linear_strain();
            let thermal=TetThermalStrainState::try_new(cell.temperature_k,reference.axes()["T"],
                [strain,strain,strain,0.0,0.0,0.0],*material_state_identity.as_bytes(),expansion.identity(),
                config.maximum_absolute_mandel_strain)?;
            operators.push((self.elements[index],material,thermal));
            sources.push(ThermalCellSource {stiffness,expansion,material_state_identity});
        }
        operators.sort_by_key(|row|row.0);
        let materials:Vec<_>=operators.iter().map(|row|row.1.clone()).collect();
        let strains:Vec<_>=operators.iter().map(|row|row.2).collect();
        let problem=TetLinearElasticProblem {nodes_m:&self.nodes_m,tetrahedra:&self.tetrahedra,
            materials:TetMaterialField::PerElement(&materials),fixed_dofs,budget:config.assembly};
        let field=TetThermalStrainField::PerElement(&strains);
        let (displacement, stress) = if nodal_forces_n.is_empty() {
            let displacement=problem.solve_thermal_displacement(field,config.solve,cx)?;
            let stress=problem.recover_thermal_stress(&displacement,field,config.mesh_frame,
                config.maximum_absolute_mandel_strain,cx)?;
            (displacement, stress)
        } else {
            let displacement=problem.solve_thermal_displacement_with_nodal_forces(
                field,nodal_forces_n,config.solve,cx)?;
            let stress=problem.recover_thermal_stress_with_nodal_forces(&displacement,
                field,nodal_forces_n,config.mesh_frame,config.maximum_absolute_mandel_strain,cx)?;
            (displacement, stress)
        };
        let support_reactions_n=loads::support_reactions(&stress,fixed_dofs,nodal_forces_n,cx)?;
        let geometry=problem.update_geometry_from_displacement(&displacement,config.maximum_nodal_displacement_m,cx)?;
        let body=self.deform(geometry.nodes_m(),||cx.checkpoint().is_err())?;
        checkpoint(cx)?;
        let linearization = vibration::Linearization {
            reference_nodes_m: self.nodes_m.clone(), materials, strains,
            fixed_dofs: fixed_dofs.to_vec(), nodal_forces_n: nodal_forces_n.to_vec(), config,
        };
        checkpoint(cx)?;
        Ok(ThermalBody {body,displacement,stress,support_reactions_n,geometry,sources,linearization})
    }
}

fn reference_point(source:&ResolvedMaterialStatePoint)->Result<QueryPoint,ThermalBodyError> {
    let mut point=QueryPoint::new();
    for (axis,value) in source.query_point() {
        point=match source.axis_quantities().get(axis) {
            Some(&quantity)=>point.with_quantity(axis.clone(),quantity,*value),
            None=>point.with(axis.clone(),*value),
        }.map_err(|_|invalid("reference source query could not be restored"))?;
    }
    let absolute=QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature,ValueForm::Static));
    if !point.axes().contains_key("T") || point.axis_quantities().get("T").is_some_and(|q|
        *q!=absolute && *q!=QuantitySpec::dimensional(Temperature::DIMS)) {
        return Err(invalid("reference T must be an absolute-kelvin coordinate"));
    }
    Ok(point)
}
fn invalid(what:&'static str)->ThermalBodyError {MaterialInertiaError::Invalid(what).into()}
fn checkpoint(cx:&Cx<'_>)->Result<(),ThermalBodyError> {
    cx.checkpoint().map_err(|_|ThermalBodyError::Solid(TetElasticError::Cancelled))
}
