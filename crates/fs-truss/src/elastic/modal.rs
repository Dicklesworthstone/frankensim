//! Prestressed vibration of the native compatible planar truss.
//!
//! One ordinary elastic solve recovers the actual axial force from nodal loads,
//! stress-free strain and support motion. The frozen-reference tangent is
//! K_t = K_e + K_g(N), sharing the bar kernel with global buckling. Original
//! density supplies consistent translational mass, independent of prestress:
//! M_e = rho A L / 6 [2I I; I 2I]. Explicit point masses may be added at joints.
//!
//! Signed spectral windows use fs-modal directly. Negative omega^2 denotes a
//! linearized instability; it is never clamped to a vibration frequency.
//! This initial-stress approximation does not follow finite deformation,
//! follower loads, local member bending, damping or changing material state.

use super::{
    ElasticControls, ElasticError, ElasticImposedState, ElasticResponse, MemberTangent, checkpoint,
    invalid,
};
use crate::{GroundStructure, LayoutCase, TrussConstructionError};
use fs_exec::Cx;
use fs_modal::{ModalError, MultiSliceReport, SliceOptions, WindowPlan};
use fs_sparse::{Coo, Csr};

/// Maximum free translational dimension admitted before equilibrium/assembly.
pub const HARD_MAX_MODAL_FREE_DOFS: usize = 100_000;
/// Combined triplet cap bounds native deterministic COO sorting/storage.
pub const HARD_MAX_MODAL_STAGED_ENTRIES: usize = 8_000_000;

/// Reference inertia in complete original member/node ordering.
#[derive(Debug, Clone, Copy)]
pub struct ModalInertia<'a> {
    /// Strictly positive finite reference density for every original member.
    /// Area zero removes that member's mass as well as its stiffness.
    pub member_density_kg_m3: &'a [f64],
    /// Optional isotropic point mass [kg] at every original node.
    /// Entries must be finite and nonnegative. None means no added point mass.
    /// Mass on supported coordinates is retained in the physical inventory.
    pub joint_mass_kg: Option<&'a [f64]>,
}

/// Bounded equilibrium and sparse pencil assembly.
/// Spectral factor/vector/output work is separately admitted by WindowPlan.
#[derive(Debug, Clone, Copy)]
pub struct ModalAssemblyControls {
    /// Existing checked-CG tolerance and total equilibrium iteration allowance.
    pub elastic: ElasticControls,
    /// Maximum free translations, at most HARD_MAX_MODAL_FREE_DOFS.
    pub maximum_free_dofs: usize,
    /// Combined nonzero tangent/mass contributions, before duplicate reduction.
    pub maximum_staged_entries: usize,
}

impl Default for ModalAssemblyControls {
    fn default() -> Self {
        Self {
            elastic: ElasticControls::default(),
            maximum_free_dofs: HARD_MAX_MODAL_FREE_DOFS,
            maximum_staged_entries: 4_000_000,
        }
    }
}

/// Admission, equilibrium and native spectral errors retain their owning cause.
#[derive(Debug)]
pub enum TrussModalError {
    /// Invalid physical data/work, unresolved equilibrium or cancellation.
    Elastic(ElasticError),
    /// Native mass, factorization, window, residual, work or cancellation refusal.
    Modal(ModalError),
}

impl From<ElasticError> for TrussModalError {
    fn from(error: ElasticError) -> Self {
        Self::Elastic(error)
    }
}
impl From<TrussConstructionError> for TrussModalError {
    fn from(error: TrussConstructionError) -> Self {
        Self::Elastic(error.into())
    }
}
impl From<ModalError> for TrussModalError {
    fn from(error: ModalError) -> Self {
        Self::Modal(error)
    }
}
impl core::fmt::Display for TrussModalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "truss prestressed vibration: {self:?}")
    }
}
impl std::error::Error for TrussModalError {}

/// Immutable physical pencil bound to its accepted equilibrium and DOF map.
///
/// No free coordinate is discarded because its axial stiffness vanishes.
/// Positive mass admission is owned by the native modal solve; a massless
/// free subsystem is not repaired with invented inertia or static condensation.
#[derive(Debug, Clone)]
pub struct TrussModalPencil {
    equilibrium: ElasticResponse,
    tangent: Csr,
    mass: Csr,
    dof_map: Vec<Option<usize>>,
    total_reference_mass_kg: f64,
}

impl TrussModalPencil {
    /// Actual load/eigenstrain/settlement equilibrium used to form K_g.
    #[must_use]
    pub fn equilibrium(&self) -> &ElasticResponse {
        &self.equilibrium
    }

    /// Symmetric free-coordinate K_e + K_g(N_actual) [N/m].
    /// This can be indefinite or singular even when the axial solve converged.
    #[must_use]
    pub fn tangent(&self) -> &Csr {
        &self.tangent
    }

    /// Consistent free-coordinate reference translational mass [kg].
    #[must_use]
    pub fn mass(&self) -> &Csr {
        &self.mass
    }

    /// Original node-major [x, y] to reduced coordinate; supported entries None.
    /// Native mode shapes use this reduced ordering. Supported perturbations
    /// vanish even when the base equilibrium has nonzero prescribed motion.
    #[must_use]
    pub fn dof_map(&self) -> &[Option<usize>] {
        &self.dof_map
    }

    /// Sum of original rho A L and explicit joint masses [kg], including
    /// supported material/joints. This is not the trace of the reduced mass.
    #[must_use]
    pub fn total_reference_mass_kg(&self) -> f64 {
        self.total_reference_mass_kg
    }

    /// Reuse this physical pencil for native signed omega^2 windows [s^-2].
    ///
    /// Modes are mass normalized, phi^T M phi = 1, in dof_map ordering.
    /// A positive lambda gives angular frequency sqrt(lambda) [rad/s] only
    /// when its reported interval resolves a positive sign. Negative lambda
    /// gives growth rate sqrt(-lambda) only when its interval resolves negative.
    /// Near-zero unresolved modes remain visible in the pencil; native relative
    /// interval-width admission may refuse rather than declare them stable.
    ///
    /// The explicit native factor and payload caps, endpoint inertia counts,
    /// residual checks, deterministic reduction and cancellation remain intact.
    /// Repeated windows do not repeat equilibrium or sparse physical assembly.
    ///
    /// # Errors
    /// Invalid windows/options, non-SPD mass, singular endpoints, unresolved
    /// eigenpairs, exhausted work/memory and cancellation return no partial
    /// report. A successful requested window is not a full-spectrum stability
    /// assertion; instability may lie outside it.
    pub fn slice_windows(
        &self,
        plan: &WindowPlan,
        options: &SliceOptions,
        cx: &Cx<'_>,
    ) -> Result<MultiSliceReport, TrussModalError> {
        Ok(fs_modal::slice_windows(
            &self.tangent,
            &self.mass,
            plan,
            options,
            cx,
        )?)
    }
}

fn stage(
    tangent: &mut Coo,
    mass: &mut Coo,
    row: usize,
    column: usize,
    stiffness: f64,
    inertia: f64,
    cap: usize,
) -> Result<(), TrussConstructionError> {
    if !stiffness.is_finite() || !inertia.is_finite() {
        return Err(invalid("truss modal pencil", "entries must remain finite"));
    }
    let staged = tangent.staged()
        + mass.staged()
        + usize::from(stiffness != 0.0)
        + usize::from(inertia != 0.0);
    if staged > cap {
        return Err(TrussConstructionError::WorkBudget {
            resource: "truss modal staged entries",
            limit: cap,
            observed: staged,
        });
    }
    if stiffness != 0.0 {
        tangent.push(row, column, stiffness);
    }
    if inertia != 0.0 {
        mass.push(row, column, inertia);
    }
    Ok(())
}

/// Recover actual prestress and assemble consistent reference inertia.
///
/// Inputs retain complete original member/node ordering. Source density is a
/// reference density: thermal/assembly strain and support settlement do not
/// create mass or update rho A L. Both Cartesian translations carry inertia;
/// projecting mass onto the member axis would lose transverse kinetic energy.
///
/// The base displacement is the existing infinitesimal compatible response.
/// Its member force supplies the shared initial-stress tangent, with tension
/// stiffening and compression softening. Original lengths/directions stay
/// frozen. Loads and imposed states remain fixed during modal perturbations;
/// a gravity design must pass its actual converged body-load case explicitly.
///
/// Every original free coordinate enters the pencil, including unloaded
/// mechanisms. Assembly alone proves neither mass positivity nor stability;
/// the native signed window solve admits mass and reports the queried spectrum.
/// Member/node loops poll every 256 entries. Native COO sorting is bounded by
/// maximum_staged_entries and polled immediately before and after.
///
/// # Errors
/// Invalid original-order data, unrepresentable physical mass/stiffness, no
/// free coordinates, exhausted work, unresolved equilibrium and cancellation
/// refuse before returning a bound pencil.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn assemble(
    ground: &GroundStructure,
    case: &LayoutCase,
    areas_m2: &[f64],
    young_pa: &[f64],
    imposed: ElasticImposedState<'_>,
    inertia: ModalInertia<'_>,
    controls: ModalAssemblyControls,
    cx: &Cx<'_>,
) -> Result<TrussModalPencil, TrussModalError> {
    checkpoint(cx)?;
    if case.node_count() != ground.nodes().len() {
        return Err(invalid("truss modal case", "must match original node ordering").into());
    }
    if controls.maximum_free_dofs == 0
        || controls.maximum_free_dofs > HARD_MAX_MODAL_FREE_DOFS
        || controls.maximum_staged_entries == 0
        || controls.maximum_staged_entries > HARD_MAX_MODAL_STAGED_ENTRIES
    {
        return Err(invalid(
            "truss modal controls",
            "require positive free-coordinate/staging allowances within declared bounds",
        )
        .into());
    }
    let full_dofs = 2 * ground.nodes().len();
    if full_dofs > HARD_MAX_MODAL_FREE_DOFS * 8 {
        return Err(TrussConstructionError::WorkBudget {
            resource: "truss modal original DOFs",
            limit: HARD_MAX_MODAL_FREE_DOFS * 8,
            observed: full_dofs,
        }
        .into());
    }
    if inertia.member_density_kg_m3.len() != ground.members().len() {
        return Err(TrussConstructionError::VectorLength {
            field: "truss modal member densities",
            expected: ground.members().len(),
            actual: inertia.member_density_kg_m3.len(),
        }
        .into());
    }
    for (index, &density) in inertia.member_density_kg_m3.iter().enumerate() {
        if index % 256 == 0 {
            checkpoint(cx)?;
        }
        if !density.is_finite() || density <= 0.0 {
            return Err(invalid(
                "truss modal reference density",
                "all original members require finite positive density",
            )
            .into());
        }
    }
    if let Some(masses) = inertia.joint_mass_kg {
        if masses.len() != ground.nodes().len() {
            return Err(TrussConstructionError::VectorLength {
                field: "truss modal joint masses",
                expected: ground.nodes().len(),
                actual: masses.len(),
            }
            .into());
        }
        for (index, &mass) in masses.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            if !mass.is_finite() || mass < 0.0 {
                return Err(
                    invalid("truss modal joint mass", "must be finite and nonnegative").into(),
                );
            }
        }
    }
    let mut dof_map = vec![None; full_dofs];
    let mut free_count = 0;
    for (node, supported) in case.supported().iter().enumerate() {
        if node % 256 == 0 {
            checkpoint(cx)?;
        }
        for (axis, &fixed) in supported.iter().enumerate() {
            if !fixed {
                dof_map[2 * node + axis] = Some(free_count);
                free_count += 1;
            }
        }
    }
    if free_count == 0 {
        return Err(invalid(
            "truss modal free DOFs",
            "require at least one free coordinate",
        )
        .into());
    }
    if free_count > controls.maximum_free_dofs {
        return Err(TrussConstructionError::WorkBudget {
            resource: "truss modal free DOFs",
            limit: controls.maximum_free_dofs,
            observed: free_count,
        }
        .into());
    }
    let equilibrium = super::solve_with_imposed(
        ground,
        case,
        areas_m2,
        young_pa,
        imposed,
        controls.elastic,
        cx,
    )?;
    let mut tangent = Coo::new(free_count, free_count);
    let mut mass = Coo::new(free_count, free_count);
    let mut total_mass = 0.0;
    for index in 0..ground.members().len() {
        if index % 256 == 0 {
            checkpoint(cx)?;
        }
        if areas_m2[index] == 0.0 {
            continue;
        }
        let length = ground.lengths()[index];
        let volume = areas_m2[index] * length;
        let member_mass = inertia.member_density_kg_m3[index] * volume;
        let diagonal_mass = member_mass / 3.0;
        let off_diagonal_mass = member_mass / 6.0;
        total_mass += member_mass;
        if !volume.is_finite()
            || volume <= 0.0
            || !member_mass.is_finite()
            || diagonal_mass <= 0.0
            || off_diagonal_mass <= 0.0
            || !total_mass.is_finite()
        {
            return Err(invalid(
                "truss modal member mass",
                "positive rho A L and its consistent entries must remain representable",
            )
            .into());
        }
        let axial_stiffness = young_pa[index] * (areas_m2[index] / length);
        let geometric = equilibrium.member_force_n[index] / length;
        if !geometric.is_finite() {
            return Err(
                invalid("truss modal geometric stiffness", "N/L must remain finite").into(),
            );
        }
        let member = MemberTangent::new(ground, index);
        for (i, &row_end) in member.ends.iter().enumerate() {
            let Some(row) = dof_map[row_end] else {
                continue;
            };
            for (j, &column_end) in member.ends.iter().enumerate() {
                let Some(column) = dof_map[column_end] else {
                    continue;
                };
                let value = if i % 2 != j % 2 {
                    0.0
                } else if i / 2 == j / 2 {
                    diagonal_mass
                } else {
                    off_diagonal_mass
                };
                stage(
                    &mut tangent,
                    &mut mass,
                    row,
                    column,
                    member.entry(axial_stiffness, geometric, i, j),
                    value,
                    controls.maximum_staged_entries,
                )?;
            }
        }
    }
    if let Some(masses) = inertia.joint_mass_kg {
        for (node, &joint_mass) in masses.iter().enumerate() {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            total_mass += joint_mass;
            if !total_mass.is_finite() {
                return Err(invalid("truss modal mass inventory", "must remain finite").into());
            }
            for axis in 0..2 {
                if let Some(row) = dof_map[2 * node + axis] {
                    stage(
                        &mut tangent,
                        &mut mass,
                        row,
                        row,
                        0.0,
                        joint_mass,
                        controls.maximum_staged_entries,
                    )?;
                }
            }
        }
    }
    checkpoint(cx)?;
    let tangent_csr = tangent.assemble();
    drop(tangent);
    checkpoint(cx)?;
    let mass_csr = mass.assemble();
    drop(mass);
    checkpoint(cx)?;
    for row in 0..free_count {
        if row % 256 == 0 {
            checkpoint(cx)?;
        }
        if tangent_csr
            .row(row)
            .1
            .iter()
            .chain(mass_csr.row(row).1)
            .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "truss modal assembled pencil",
                "reduced sums must remain finite",
            )
            .into());
        }
    }
    checkpoint(cx)?;
    Ok(TrussModalPencil {
        equilibrium,
        tangent: tangent_csr,
        mass: mass_csr,
        dof_map,
        total_reference_mass_kg: total_mass,
    })
}
