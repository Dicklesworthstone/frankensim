//! Original moving-current transport coupled to a sourced rigid-body momentum kick.
//!
//! This is a fixed-configuration implicit-midpoint FORCE substep. The native
//! impulse owner updates linear/angular momentum; the same midpoint velocity
//! changes the electrical current through `u cross B`. Pose, geometry, material
//! conditions, boundary traces and imposed B remain fixed. Free rigid transport,
//! gyroscopic/advective evolution, induction and thermal evolution belong to
//! their existing owners; this substep is not an ALE or full-flight integrator.
//!
//! Preparation retains the original field and exact tetrahedral mass inventory.
//! Six homogeneous native electrical solves produce the complete rigid-velocity
//! response, not fitted motor constants. Every warm kick re-solves the original
//! electrical field at the coupled midpoint before actual native impulse and
//! energy admission. No field, inertia, force law or Joule source is replaced.

mod prepare;

use fs_em::{
    conduction::ConductionControls,
    spatial::conduction::{
        SpatialConduction, SpatialConductionBoundary, SpatialConductionError,
        motion::{
            SpatialConductionMotion, SpatialLorentzLoadError, SpatialLorentzLoadLimits,
            SpatialMovingConductionField, SpatialMovingConductionSnapshot,
        },
    },
};
use fs_ga::{Se3, Twist, Vec3, Wrench};
use fs_la::{LuWorkspace, LuWorkspaceError};
use fs_material::conductor::{ResolvedConductorState, spatial::SpatialConductorMaterial};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, FreeFloatingBaseState, FreeFloatingImpulseResponse,
    JointModel, Link, apply_free_floating_impulses,
};
use fs_mbd_material::deformation::TetrahedralMass;

/// Original electrical work bounds plus independent mechanical acceptance.
#[derive(Clone, Copy, Debug)]
pub struct SpatialConductionKickControls {
    pub field: ConductionControls,
    pub loads: SpatialLorentzLoadLimits,
    /// Preparation requires seven full electrical solves; a warm kick uses one.
    /// This cap applies separately to each call, in addition to native CG limits.
    pub maximum_field_solves: usize,
    /// Maximum mismatch between field and actual mechanical midpoint [rad/s].
    pub angular_velocity_tolerance_rad_s: f64,
    /// Maximum mismatch between field and actual mechanical midpoint [m/s].
    pub linear_velocity_tolerance_m_s: f64,
    /// Absolute Lorentz/kinetic and complete source/heat/storage defect [J].
    pub energy_tolerance_j: f64,
}

/// Refusal publishes no electrical or mechanical candidate.
#[derive(Debug)]
pub enum SpatialConductionKickError {
    Invalid(&'static str),
    Cancelled,
    Electrical(SpatialConductionError),
    Loads(SpatialLorentzLoadError),
    Dynamics(ArticulatedError),
    Linear(LuWorkspaceError),
    Balance {
        angular_velocity_residual_rad_s: f64,
        linear_velocity_residual_m_s: f64,
        mechanical_work_residual_j: f64,
        total_energy_residual_j: f64,
    },
}
impl core::fmt::Display for SpatialConductionKickError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial conduction momentum kick: {self:?}")
    }
}
impl std::error::Error for SpatialConductionKickError {}
impl From<SpatialConductionError> for SpatialConductionKickError {
    fn from(error: SpatialConductionError) -> Self {
        match error {
            SpatialConductionError::Cancelled => Self::Cancelled,
            error => Self::Electrical(error),
        }
    }
}
impl From<SpatialLorentzLoadError> for SpatialConductionKickError {
    fn from(error: SpatialLorentzLoadError) -> Self {
        match error {
            SpatialLorentzLoadError::Cancelled => Self::Cancelled,
            error => Self::Loads(error),
        }
    }
}
impl From<ArticulatedError> for SpatialConductionKickError {
    fn from(error: ArticulatedError) -> Self {
        Self::Dynamics(error)
    }
}
impl From<LuWorkspaceError> for SpatialConductionKickError {
    fn from(error: LuWorkspaceError) -> Self {
        Self::Linear(error)
    }
}
type Result<T> = core::result::Result<T, SpatialConductionKickError>;

/// Complete velocity/current endpoint on the unchanged source configuration.
#[derive(Debug)]
pub struct SpatialConductionKick {
    /// Original native impulse result, including actual kinetic/momentum change.
    pub mechanical: FreeFloatingImpulseResponse,
    /// Same accepted electrical solve; its Joule heat is NOT deposited here.
    pub field: SpatialMovingConductionSnapshot,
    /// Body-frame rigid velocity supplied to the original motional-EMF solve.
    /// Its difference from the actual native endpoint mean is bounded below.
    pub midpoint_twist_body: Twist,
    /// Actual solved-current load about the original body origin, in body axes.
    /// Already included in `mechanical`; do not apply it a second time.
    pub lorentz_wrench_body: Wrench,
    pub terminal_work_j: f64,
    pub joule_heat_j: f64,
    pub mechanical_work_j: f64,
    pub external_work_j: f64,
    pub angular_velocity_residual_rad_s: f64,
    pub linear_velocity_residual_m_s: f64,
    pub mechanical_work_residual_j: f64,
    /// `delta_kinetic + Joule - terminal - external` [J].
    pub total_energy_residual_j: f64,
    pub field_solves: usize,
}

/// Prepared original-mesh momentum coupling for ONE held field configuration.
///
/// The borrowed current/material operator and tetrahedral source mass cannot be
/// replaced while this cache lives. Boundary and B snapshots are owned; the
/// fixed native pose must match every incoming state. Changing geometry, B,
/// conductivity, boundary law, pose or duration requires new preparation.
/// Warm calls alter only private algebra scratch, never accepted history.
pub struct PreparedSpatialConductionKick<'a, C = ResolvedConductorState> {
    problem: &'a SpatialConduction<C>,
    mass: &'a TetrahedralMass<'a>,
    model: ArticulatedModel,
    world_from_body: Se3,
    boundary: SpatialConductionBoundary,
    flux_density_t: Vec<[f64; 3]>,
    selected_cells: Vec<usize>,
    duration_s: f64,
    controls: SpatialConductionKickControls,
    source_wrench: [f64; 6],
    // Rows [torque, force], columns [angular, linear] rigid velocity.
    response: [f64; 36],
    schur: LuWorkspace,
}

impl<C: SpatialConductorMaterial> PreparedSpatialConductionKick<'_, C> {
    pub fn problem(&self) -> &SpatialConduction<C> {
        self.problem
    }
    pub fn mass(&self) -> &TetrahedralMass<'_> {
        self.mass
    }
    pub fn model(&self) -> &ArticulatedModel {
        &self.model
    }
    pub fn duration_s(&self) -> f64 {
        self.duration_s
    }

    /// Apply `M*(v_after-v_before)=dt*(F_Lorentz(v_mid)+F_external)` through
    /// the original native impulse owner, with `v_mid=(v_before+v_after)/2`.
    ///
    /// `external_wrench_body` is a held physical load about the same body
    /// origin, in the same axes, and has its own signed work. There is no
    /// additional Maxwell/induction load. The returned pose is unchanged;
    /// callers must separately own geometry/material transport and reprepare
    /// after a configuration change. Only `joule_heat_j` belongs in a thermal
    /// receiver; the kinetic increment is physical mechanical storage.
    ///
    /// The six-dimensional solve changes only a rigid midpoint velocity. The
    /// final field is recovered with the complete original sparse solve, and
    /// the actual native impulse endpoint must pass independent velocity,
    /// Lorentz-work and combined energy gates before anything is returned.
    pub fn kick_candidate(
        &mut self,
        initial: FreeFloatingBaseState,
        external_wrench_body: Wrench,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialConductionKick> {
        poll(cancelled)?;
        if initial.world_from_base != self.world_from_body {
            return Err(invalid(
                "incoming pose differs from the prepared field configuration",
            ));
        }
        let old = initial.twist_body.to_array();
        let external = external_wrench_body.to_array();
        let mut rhs = [0.0; 6];
        for row in 0..6 {
            finite(old[row])?;
            let mut force = finite(self.source_wrench[row] + finite(external[row])?)?;
            for col in 0..6 {
                force = finite(self.response[row * 6 + col].mul_add(old[col], force))?;
            }
            rhs[row] = finite((0.5 * self.duration_s) * force)?;
        }
        // Solve for the midpoint INCREMENT, preserving small changes beside
        // large unchanged motion instead of subtracting two M*v vectors.
        let mut increment = [0.0; 6];
        self.schur.solve_factored_into(&rhs, &mut increment)?;
        let mut midpoint = [0.0; 6];
        for i in 0..6 {
            midpoint[i] = finite(old[i] + increment[i])?;
        }
        let midpoint_twist_body = twist(midpoint);
        let (field, lorentz_wrench_body) = solve_field(
            self.problem,
            &self.boundary,
            &self.flux_density_t,
            &self.selected_cells,
            midpoint_twist_body,
            self.controls,
            cancelled,
        )?;
        let applied = core::array::from_fn(|i| {
            self.duration_s * (lorentz_wrench_body.to_array()[i] + external[i])
        });
        for value in applied {
            finite(value)?;
        }
        poll(cancelled)?;
        let mechanical =
            apply_free_floating_impulses(&self.model, initial, &[], &[], &[], &[wrench(applied)])?;
        let after = mechanical.base_after.twist_body.to_array();
        let actual_midpoint = core::array::from_fn(|i| 0.5 * old[i] + 0.5 * after[i]);
        let mut angular_residual: f64 = 0.0;
        let mut linear_residual: f64 = 0.0;
        for i in 0..6 {
            let residual = finite(actual_midpoint[i] - midpoint[i])?.abs();
            if i < 3 {
                angular_residual = angular_residual.max(residual);
            } else {
                linear_residual = linear_residual.max(residual);
            }
        }
        let terminal_work_j = finite(self.duration_s * field.terminal_power_w())?;
        let joule_heat_j = finite(self.duration_s * field.joule_power_w())?;
        let mechanical_work_j = finite(self.duration_s * field.mechanical_power_w())?;
        let external_work_j =
            finite(self.duration_s * external_wrench_body.pairing(twist(actual_midpoint)))?;
        let mechanical_work_residual_j = finite(
            finite(mechanical.kinetic_energy_change - external_work_j)? - mechanical_work_j,
        )?;
        let total_energy_residual_j = finite(
            finite(finite(mechanical.kinetic_energy_change + joule_heat_j)? - terminal_work_j)?
                - external_work_j,
        )?;
        if (field.joule_power_w() > 0.0 && joule_heat_j == 0.0)
            || angular_residual > self.controls.angular_velocity_tolerance_rad_s
            || linear_residual > self.controls.linear_velocity_tolerance_m_s
            || mechanical.kinetic_work_residual.abs() > self.controls.energy_tolerance_j
            || mechanical_work_residual_j.abs() > self.controls.energy_tolerance_j
            || total_energy_residual_j.abs() > self.controls.energy_tolerance_j
        {
            return Err(SpatialConductionKickError::Balance {
                angular_velocity_residual_rad_s: angular_residual,
                linear_velocity_residual_m_s: linear_residual,
                mechanical_work_residual_j,
                total_energy_residual_j,
            });
        }
        let field = field.into_snapshot();
        poll(cancelled)?;
        Ok(SpatialConductionKick {
            mechanical,
            field,
            midpoint_twist_body,
            lorentz_wrench_body,
            terminal_work_j,
            joule_heat_j,
            mechanical_work_j,
            external_work_j,
            angular_velocity_residual_rad_s: angular_residual,
            linear_velocity_residual_m_s: linear_residual,
            mechanical_work_residual_j,
            total_energy_residual_j,
            field_solves: 1,
        })
    }
}

fn invalid(reason: &'static str) -> SpatialConductionKickError {
    SpatialConductionKickError::Invalid(reason)
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("nonfinite coupled arithmetic"))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(SpatialConductionKickError::Cancelled)
    } else {
        Ok(())
    }
}
fn vector([x, y, z]: [f64; 3]) -> Vec3 {
    Vec3::new(x, y, z)
}
fn twist(v: [f64; 6]) -> Twist {
    Twist::new(Vec3::new(v[0], v[1], v[2]), Vec3::new(v[3], v[4], v[5]))
}
fn wrench(v: [f64; 6]) -> Wrench {
    Wrench::new(Vec3::new(v[0], v[1], v[2]), Vec3::new(v[3], v[4], v[5]))
}

#[allow(clippy::too_many_arguments)]
fn solve_field<'a, C: SpatialConductorMaterial>(
    problem: &'a SpatialConduction<C>,
    boundary: &SpatialConductionBoundary,
    flux_density_t: &[[f64; 3]],
    selected_cells: &[usize],
    velocity: Twist,
    controls: SpatialConductionKickControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(SpatialMovingConductionField<'a, C>, Wrench)> {
    let mut nodal = Vec::with_capacity(problem.positions_m().len());
    for &position in problem.positions_m() {
        poll(cancelled)?;
        let v = velocity.linear + velocity.angular.cross(vector(position));
        let value = [v.x, v.y, v.z];
        for component in value {
            finite(component)?;
        }
        nodal.push(value);
    }
    let mut cell_velocity = Vec::with_capacity(problem.mesh().tets.len());
    for tet in &problem.mesh().tets {
        poll(cancelled)?;
        let mut mean = [0.0; 3];
        for axis in 0..3 {
            mean[axis] = finite(tet.iter().map(|&v| 0.25 * nodal[v as usize][axis]).sum())?;
        }
        cell_velocity.push(mean);
    }
    let field = problem.solve_with_motion(
        boundary,
        SpatialConductionMotion {
            velocity_m_s: &cell_velocity,
            flux_density_t,
        },
        controls.field,
        cancelled,
    )?;
    let loads =
        field.lorentz_loads_for_velocity(selected_cells, &nodal, 0.0, controls.loads, cancelled)?;
    let wrench = loads.wrench(&Se3::identity(), [0.0; 3])?;
    if finite(wrench.pairing(velocity) - field.mechanical_power_w())?.abs()
        > controls.field.power_tolerance_w
    {
        return Err(invalid(
            "native rigid wrench and original-cell Lorentz work disagree",
        ));
    }
    poll(cancelled)?;
    Ok((field, wrench))
}
