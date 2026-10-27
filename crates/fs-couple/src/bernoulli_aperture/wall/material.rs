//! Geometry + original material receipts -> the existing locally reacting wall.
//!
//! The axisymmetric thin-cylinder membrane limit gives hoop strain x/a and
//! sigma_theta = p*a/h. With zero axial stress, K=E*h/a^2; with zero axial
//! strain, K=E*h/[a^2*(1-nu^2)]. Surface mass is rho*h. An explicitly declared
//! Kelvin--Voigt Young viscosity uses the SAME constant Poisson reduction.
//! No frequency-independent loss factor is relabeled as viscosity.
//!
//! Principal-aligned orthotropic walls use the same local membrane owner as
//! DG ducts, retaining the original directional moduli and material axes.
//! This is not a carved shell, off-axis grain, moving geometry or exterior
//! wall radiation. Axial/bending coupling between rings, end-cap
//! traction, static-pressure prestress and thermal evolution are absent. It is
//! a frozen-state, small-displacement thin-wall image, not a claim
//! that every change of bore metal must produce an audible timbre difference.

/// Hereditary material spectra through the shared passive acoustic load.
pub mod relaxation;
pub use relaxation::WallRelaxation;
pub(crate) mod source;
mod local;
pub use local::{LocalMaterialWall, MaterialWallDg, RegionalMaterialWallDg};
pub use local::finite_gas::{
    FiniteGasDgBackground, FiniteGasDgError, FiniteGasDgLoss, FiniteGasDgSolve,
    FiniteGasDgStep,
};

use fs_vfit::relaxation::RelaxationImpedanceSpec;
use super::{WallPatch, invalid};
use super::super::network::{NetworkNode, TubeNetworkSpec};
use crate::acoustic_realize::AcousticRealizeError;
use fs_exec::CancelGate;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_qty::{Dims, QuantitySpec};

/// Uniaxial solid Young viscosity [Pa s], NOT fluid or bending viscosity.
/// Its isotropic Kelvin--Voigt law uses the elastic Poisson ratio unchanged.
pub const YOUNG_VISCOSITY_PROPERTY: &str = "kelvin_voigt_young_viscosity";

/// Actual axial restraint, not a material-name heuristic or an end-cap model.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AxialConstraint {
    /// Incremental axial stress is zero; axial inertia is neglected.
    StressFree,
    /// Incremental axial strain is zero; the support supplies the axial reaction.
    StrainFixed,
}

/// Absence of damping is explicit, never replaced by a guessed loss factor.
#[derive(Clone, Debug, PartialEq)]
pub enum WallMaterialLoss {
    /// Lossless membrane. Other network loss channels are unchanged.
    Elastic,
    /// Viscous response from the explicit constitutive selection. Isotropic
    /// walls require `kelvin_voigt_young_viscosity`; directional walls retain
    /// their declared proportional time or normal-plane tensor properties.
    KelvinVoigt,
    /// Existing generalized-Maxwell law; retains its original spectrum and conditions.
    /// Isotropic source spectra add no Kelvin damping. Exact normal-plane
    /// elimination may retain its source-derived residual Kelvin resistance.
    Maxwell(WallRelaxation),
}

/// Explicit constitutive image of a local cylindrical membrane.
/// Material names never choose this reduction or its principal directions.
#[derive(Clone, Debug, Default, PartialEq)]
pub enum WallMaterialModel {
    /// Original isotropic Young/Poisson reduction and loss source.
    #[default]
    Isotropic,
    /// Material principal axes follow the local [hoop, axial] directions.
    /// A time property declares C_viscous=tau*C_elastic on this plane.
    /// Off-axis grain and general directional relaxation are not approximated.
    Orthotropic {
        /// Distinct zero-based principal material directions in 0..3.
        axes: [usize; 2],
        /// None selects elastic response; Some requires KelvinVoigt loss.
        proportional_viscosity_time_property: Option<String>,
    },
    /// Independent symmetric normal-plane Kelvin--Voigt viscosity. Eliminating
    /// free axial stress retains its exact Maxwell memory and Kelvin remainder.
    OrthotropicNormal {
        /// Distinct zero-based local [hoop, axial] principal directions.
        axes: [usize; 2],
        /// Original [D_hh,D_ha,D_aa] source properties, each in Pa s.
        properties: [String; 3],
        /// Explicit midpoint-network relaxation resolution, in (0,2].
        max_dt_over_tau: f64,
        /// Explicit midpoint-network angular-step limit, in (0,1].
        max_angular_step: f64,
    },
}

impl WallMaterialModel {
    fn resolve(&self, material: ResolvedMaterialStatePoint, thickness_m: f64,
        axial: AxialConstraint, loss: WallMaterialLoss, band_hz: [f64; 2])
        -> Result<LocalMaterialWall, AcousticRealizeError>
    {
        match self {
            Self::Isotropic => LocalMaterialWall::new(material, thickness_m, axial, loss, band_hz),
            Self::Orthotropic { axes, proportional_viscosity_time_property } => {
                let expected = if proportional_viscosity_time_property.is_some() {
                    WallMaterialLoss::KelvinVoigt
                } else { WallMaterialLoss::Elastic };
                if loss != expected {
                    return Err(invalid("orthotropic wall loss must match its explicit proportional viscosity source"));
                }
                LocalMaterialWall::new_orthotropic(material, thickness_m, axial, *axes,
                    proportional_viscosity_time_property.as_deref(), band_hz)
            }
            Self::OrthotropicNormal { axes, properties, max_dt_over_tau, max_angular_step } => {
                let wall = LocalMaterialWall::new_orthotropic_normal_kelvin_voigt(material,
                    thickness_m, axial, *axes, properties.each_ref().map(String::as_str),
                    band_hz, *max_dt_over_tau, *max_angular_step)?;
                if &loss != wall.loss() {
                    return Err(invalid("normal wall viscosity must retain the exact source-derived relaxation spectrum"));
                }
                Ok(wall)
            }
        }
    }
}

/// One source assignment to a physical, degree-two cylindrical junction.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialWallAssignment {
    /// Original graph node; original section/node indices are never renumbered.
    pub node: usize,
    /// Exact material state with the selected model's dimensioned properties.
    pub material: ResolvedMaterialStatePoint,
    /// Explicit isotropic or principal-aligned membrane reduction.
    pub model: WallMaterialModel,
    /// Uniform radial wall thickness [m] at this already resolved state.
    pub thickness_m: f64,
    /// Explicit structural boundary assumption of the local ring image.
    pub axial: AxialConstraint,
    /// Explicit internal loss model, independent of gas viscothermal loss.
    pub loss: WallMaterialLoss,
}

/// Source assignments and explicit cold-admission limits, shared by source types.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialWallInput {
    /// Original graph addresses and independent material states.
    pub assignments: Vec<MaterialWallAssignment>,
    /// Positive applicability interval [Hz]; not a low-pass or error certificate.
    pub band_hz: [f64; 2],
    /// Maximum admitted assignments, bounded above by 1024.
    pub max_walls: usize,
}

/// Original state/geometry and its work-conjugate pressure/volume-flow law.
#[derive(Clone, Debug, PartialEq)]
pub struct BoundMaterialWall {
    source: MaterialWallAssignment,
    radius_m: f64,
    length_m: f64,
    band_hz: [f64; 2],
    patch: WallPatch,
    load: RelaxationImpedanceSpec,
}
impl BoundMaterialWall {
    /// Original material receipts, structural restraint and loss selection.
    #[must_use]
    pub const fn source(&self) -> &MaterialWallAssignment { &self.source }
    /// Actual bore radius, not an equivalent stiffness radius [m].
    #[must_use]
    pub const fn radius_m(&self) -> f64 { self.radius_m }
    /// Wall's dual-cell length: half of each adjacent physical section [m].
    #[must_use]
    pub const fn length_m(&self) -> f64 { self.length_m }
    /// Declared frequency applicability [Hz], not an antialiasing filter.
    #[must_use]
    pub const fn band_hz(&self) -> [f64; 2] { self.band_hz }
    /// Base inertia/equilibrium/viscous wall law, excluding hereditary arms.
    /// Use `shunt()` or `load()` for the COMPLETE pressure/flow connection.
    #[must_use]
    pub const fn patch(&self) -> WallPatch { self.patch }
}

fn poll(gate: &CancelGate) -> Result<(), AcousticRealizeError> {
    if gate.is_requested() { Err(AcousticRealizeError::Cancelled) } else { Ok(()) }
}

// Only constant coefficients at the original fixed material state are admitted.
// Every selected claim's frequency domain must cover the whole requested band.
fn property(material: &ResolvedMaterialStatePoint, name: &str, dims: Dims,
    band: [f64; 2]) -> Result<f64, AcousticRealizeError>
{
    let p = material.property(name).ok_or_else(|| invalid("wall material is missing a required property"))?;
    if p.requirement().quantity() != QuantitySpec::dimensional(dims)
        || !p.is_constant_at_fixed_temperature() || !p.value_si().is_finite() {
        return Err(invalid("wall material requires dimensioned, fixed-state constant coefficients"));
    }
    if let Some((lo, hi)) = p.angular_frequency_band_rad_s()
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?
        && (core::f64::consts::TAU * band[0] < lo || core::f64::consts::TAU * band[1] > hi) {
        return Err(invalid("wall material does not cover the declared frequency band"));
    }
    Ok(p.value_si())
}

fn reduce(source: MaterialWallAssignment, radius: f64, length: f64, band: [f64; 2])
    -> Result<BoundMaterialWall, AcousticRealizeError>
{
    let law = source.model.resolve(source.material.clone(), source.thickness_m,
        source.axial, source.loss.clone(), band)?;
    let patch = law.base_patch(radius, length)?;
    let load = law.load(patch)?;
    Ok(BoundMaterialWall { source, radius_m: radius, length_m: length, band_hz: band, patch, load })
}

/// Bind material walls without replacing gas, lengths, loads or numerical owners.
///
/// Each selected node must be an EMPTY degree-two junction joining equal-radius
/// cylinders. Its wall strip covers half of each adjacent section. Adjacent
/// selected nodes therefore do not double count wall area. Unselected half-strips
/// remain rigid, including inlet/outlet half-strips. Subdivide the input geometry
/// explicitly to refine the local-wall approximation. Branch junctions, existing
/// loads and radius steps refuse instead of inventing a shell or replacing a load.
///
/// Node and section indices are unchanged. Existing regional gas assignments and
/// viscothermal section selections remain valid. The returned walls retain their
/// sources; callers must keep them with the runtime rather than discard receipts.
/// The given band is a model applicability statement, not a broadband certificate.
///
/// # Errors
/// Invalid/missing material data or geometry, duplicate/over-budget assignments,
/// incompatible topology, unsupported constitutive domain or cancellation.
/// All changes are on owned cold input; no running state is mutated.
pub fn bind_material_walls(
    mut graph: TubeNetworkSpec, assignments: Vec<MaterialWallAssignment>, band_hz: [f64; 2],
    max_walls: usize, gate: &CancelGate,
) -> Result<(TubeNetworkSpec, Vec<BoundMaterialWall>), AcousticRealizeError> {
    poll(gate)?;
    if !band_hz.iter().all(|f| f.is_finite() && *f > 0.0)
        || band_hz[0] >= band_hz[1] || !(core::f64::consts::TAU * band_hz[1]).is_finite()
        || max_walls > 1024 || assignments.len() > max_walls {
        return Err(invalid("material walls require a positive band and at most 1024 bounded assignments"));
    }
    let mut walls = Vec::new();
    walls.try_reserve_exact(assignments.len()).map_err(|_| invalid("material wall allocation failed"))?;
    for source in assignments {
        poll(gate)?;
        if graph.nodes.get(source.node) != Some(&NetworkNode::Junction) {
            return Err(invalid("material wall must bind an unused interior junction, never replace a load"));
        }
        let mut radius = None;
        let mut length = 0.0;
        let mut count = 0;
        for section in &graph.sections {
            poll(gate)?;
            if section.nodes.contains(&source.node) {
                if section.nodes[0] == section.nodes[1] || section.nodes.iter().any(|&n| n >= graph.nodes.len())
                    || !section.length_m.is_finite() || section.length_m <= 0.0
                    || !section.radius_m.is_finite() || section.radius_m <= 0.0 {
                    return Err(invalid("material wall adjacency requires physical positive-length cylinders"));
                }
                if radius.is_some_and(|r: f64| r.to_bits() != section.radius_m.to_bits()) {
                    return Err(invalid("a local cylindrical material wall cannot span a radius discontinuity"));
                }
                radius = Some(section.radius_m);
                length += 0.5 * section.length_m;
                count += 1;
            }
        }
        if count != 2 { return Err(invalid("material wall requires exactly two adjacent tube sections")); }
        let wall = reduce(source, radius.ok_or_else(|| invalid("wall has no radius"))?, length, band_hz)?;
        graph.nodes[wall.source.node] = wall.shunt();
        walls.push(wall);
    }
    poll(gate)?;
    Ok((graph, walls))
}

