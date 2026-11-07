//! Shared finite-inertia stage for one supported induction rotor.
//!
//! Spatial and planar field owners retain their own material geometry, field
//! evolution and bounded coupling iteration. This stage applies their full
//! native Lorentz loads through the original articulated impulse map. It owns
//! neither a new mechanics integrator nor an approximate scalar torque law.

use fs_em::spatial::vector_potential::induction::mechanical::{
    SpatialLorentzLoadError, SpatialLorentzLoads,
};
use fs_ga::{GaError, Se3, So3, Twist, Vec3, Wrench};
use fs_math::dd::Dd;
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, JointKind, PrescribedBaseImpulseResponse,
    apply_prescribed_base_impulses, forward_kinematics,
};

/// Coupling convergence and physical-accounting allowances, independent of the
/// existing field Krylov/energy/Peclet controls.
#[derive(Clone, Copy, Debug)]
pub struct RotorControls {
    pub maximum_iterations: usize,
    pub relaxation: f64,
    pub angular_velocity_tolerance_rad_s: f64,
    pub energy_tolerance_j: f64,
    pub angular_impulse_tolerance_nms: f64,
    pub linear_impulse_tolerance_ns: f64,
}

/// Diagnostics of the actual accepted joint/field stage, never predictor loads.
#[derive(Clone, Copy, Debug)]
pub struct RotorBalance {
    pub iterations: usize,
    pub field_angular_velocity_rad_s: f64,
    pub mechanical_midpoint_velocity_rad_s: f64,
    pub velocity_residual_rad_s: f64,
    pub electromagnetic_torque_nm: f64,
    /// Actual midpoint speed times the held applied shaft torque and duration.
    pub shaft_work_j: f64,
    /// Native P1 field work minus duration * field speed * integrated torque.
    pub load_work_residual_j: f64,
    /// Actual kinetic increment minus native Lorentz work and shaft work.
    pub mechanical_work_residual_j: f64,
    /// Field energy residual plus the independently reconstructed kinetic work.
    pub total_energy_residual_j: f64,
}

#[derive(Debug)]
pub(crate) enum RotorStageError {
    Invalid(&'static str),
    Cancelled,
    Mechanical(ArticulatedError),
    Load(SpatialLorentzLoadError),
    Geometry(GaError),
    Balance(&'static str),
}
impl From<ArticulatedError> for RotorStageError {
    fn from(error: ArticulatedError) -> Self {
        Self::Mechanical(error)
    }
}
impl From<SpatialLorentzLoadError> for RotorStageError {
    fn from(error: SpatialLorentzLoadError) -> Self {
        Self::Load(error)
    }
}
impl From<GaError> for RotorStageError {
    fn from(error: GaError) -> Self {
        Self::Geometry(error)
    }
}
type Result<T> = core::result::Result<T, RotorStageError>;

#[derive(Clone, Copy, Debug)]
pub(crate) struct RotorEndpoint {
    pub angle_rad: f64,
    pub angular_velocity_rad_s: f64,
}

pub(crate) struct RotorStage {
    pub mechanics: PrescribedBaseImpulseResponse,
    pub balance: RotorBalance,
    /// Present only after velocity, full impulse and all work accounts pass.
    pub endpoint: Option<RotorEndpoint>,
}

fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(RotorStageError::Invalid("nonfinite rotor arithmetic"))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(RotorStageError::Cancelled)
    } else {
        Ok(())
    }
}
fn base() -> BaseState {
    BaseState::stationary(Se3::identity())
}

/// Admit one centered, axisymmetric revolute-z root without prescribing the
/// field's material occupancy. The returned axis origin is in world coordinates.
pub(crate) fn admit_model(
    model: &ArticulatedModel,
    angle_rad: f64,
    angular_velocity_rad_s: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec3> {
    poll(cancelled)?;
    if model.link_count() != 1 || model.dof_count() != 1 {
        return Err(RotorStageError::Invalid(
            "one native revolute rotor required",
        ));
    }
    let link = &model.links()[0];
    let inertia = link.inertia();
    let tensor = inertia.inertia_com().m;
    if link.joint().kind() != JointKind::Revolute
        || link.joint().motion_subspace() != Twist::new(Vec3::new(0.0, 0.0, 1.0), Vec3::default())
        || link.parent_from_child_zero().rotation() != So3::identity()
        || inertia.center_of_mass() != Vec3::default()
        || tensor[0] != tensor[4]
        || [1, 2, 3, 5, 6, 7].into_iter().any(|i| tensor[i] != 0.0)
    {
        return Err(RotorStageError::Invalid(
            "centered axisymmetric native inertia and aligned revolute-z root required",
        ));
    }
    forward_kinematics(model, base(), &[angle_rad], &[angular_velocity_rad_s])?;
    poll(cancelled)?;
    Ok(link.parent_from_child_zero().translation())
}

/// Admit bounded coupling policy before any field trial. Original joint effort
/// limits remain the caller's one-time `model.validate_applied_effort` check.
pub(crate) fn admit_controls(
    duration_s: f64,
    shaft_torque_nm: f64,
    controls: RotorControls,
) -> Result<()> {
    if !duration_s.is_finite()
        || duration_s <= 0.0
        || !shaft_torque_nm.is_finite()
        || !(1..=256).contains(&controls.maximum_iterations)
        || !controls.relaxation.is_finite()
        || controls.relaxation <= 0.0
        || controls.relaxation > 1.0
        || [
            controls.angular_velocity_tolerance_rad_s,
            controls.energy_tolerance_j,
            controls.angular_impulse_tolerance_nms,
            controls.linear_impulse_tolerance_ns,
        ]
        .into_iter()
        .any(|x| !x.is_finite() || x <= 0.0)
    {
        return Err(RotorStageError::Invalid(
            "positive finite coupling controls and interval required",
        ));
    }
    Ok(())
}

/// Evaluate exactly one mechanical stage from the original accepted speed.
/// A nonaccepted endpoint is never fed into another trial by this helper.
/// The owner uses the returned mechanical midpoint only as its next bounded
/// field-velocity proposal and retains its accepted histories until convergence.
#[allow(clippy::too_many_arguments)]
pub(crate) fn stage(
    model: &ArticulatedModel,
    angle_rad: f64,
    old_speed_rad_s: f64,
    field_speed_rad_s: f64,
    duration_s: f64,
    shaft_torque_nm: f64,
    loads: &SpatialLorentzLoads,
    field_mechanical_work_j: f64,
    field_energy_residual_j: f64,
    iteration: usize,
    controls: RotorControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<RotorStage> {
    poll(cancelled)?;
    admit_controls(duration_s, shaft_torque_nm, controls)?;
    if iteration == 0
        || iteration > controls.maximum_iterations
        || model.link_count() != 1
        || model.dof_count() != 1
    {
        return Err(RotorStageError::Invalid(
            "one bounded native rotor stage required",
        ));
    }
    let origin = model.links()[0].parent_from_child_zero().translation();
    let world = loads.wrench(&Se3::identity(), [origin.x, origin.y, origin.z])?;
    let midpoint_angle = finite(angle_rad + 0.5 * duration_s * field_speed_rad_s)?;
    let kinematics = forward_kinematics(model, base(), &[midpoint_angle], &[field_speed_rad_s])?;
    let body = loads.wrench(&kinematics.world_from_link[0].inverse()?, [0.0; 3])?;
    let response = apply_prescribed_base_impulses(
        model,
        base(),
        &[midpoint_angle],
        &[old_speed_rad_s],
        &[finite(duration_s * shaft_torque_nm)?],
        &[Wrench::new(
            body.torque.scale(duration_s),
            body.force.scale(duration_s),
        )],
    )?;
    poll(cancelled)?;
    let endpoint_speed = response.joint_velocity_after[0];
    let actual_midpoint = f64::midpoint(old_speed_rad_s, endpoint_speed);
    let velocity_residual_rad_s = finite(actual_midpoint - field_speed_rad_s)?.abs();
    let shaft_work_j = finite(duration_s * shaft_torque_nm * actual_midpoint)?;
    let load_work_residual_j = finite(
        (Dd::from_f64(field_mechanical_work_j)
            - Dd::from_f64(duration_s)
                * Dd::from_f64(field_speed_rad_s)
                * Dd::from_f64(world.torque.z))
        .to_f64(),
    )?;
    let mechanical_work_residual_j = finite(
        (Dd::from_f64(response.kinetic_energy_change)
            - Dd::from_f64(field_mechanical_work_j)
            - Dd::from_f64(shaft_work_j))
        .to_f64(),
    )?;
    let total_energy_residual_j = finite(
        (Dd::from_f64(field_energy_residual_j) + Dd::from_f64(mechanical_work_residual_j)).to_f64(),
    )?;
    if load_work_residual_j.abs() > controls.energy_tolerance_j {
        return Err(RotorStageError::Balance(
            "native Lorentz torque and native field work disagree",
        ));
    }
    let endpoint = if velocity_residual_rad_s <= controls.angular_velocity_tolerance_rad_s
        && mechanical_work_residual_j.abs() <= controls.energy_tolerance_j
        && total_energy_residual_j.abs() <= controls.energy_tolerance_j
    {
        if response.angular_momentum_residual_world.norm() > controls.angular_impulse_tolerance_nms
            || response.linear_momentum_residual_world.norm() > controls.linear_impulse_tolerance_ns
            || response.kinetic_work_residual.abs() > controls.energy_tolerance_j
        {
            return Err(RotorStageError::Balance(
                "native mechanical impulse/work residual exceeds coupling budget",
            ));
        }
        let angle = finite(angle_rad + duration_s * actual_midpoint)?;
        if duration_s * actual_midpoint != 0.0 && angle == angle_rad {
            return Err(RotorStageError::Invalid(
                "rotor angle increment disappeared in endpoint",
            ));
        }
        forward_kinematics(model, base(), &[angle], &[endpoint_speed])?;
        Some(RotorEndpoint {
            angle_rad: angle,
            angular_velocity_rad_s: endpoint_speed,
        })
    } else {
        None
    };
    let balance = RotorBalance {
        iterations: iteration,
        field_angular_velocity_rad_s: field_speed_rad_s,
        mechanical_midpoint_velocity_rad_s: actual_midpoint,
        velocity_residual_rad_s,
        electromagnetic_torque_nm: world.torque.z,
        shaft_work_j,
        load_work_residual_j,
        mechanical_work_residual_j,
        total_energy_residual_j,
    };
    poll(cancelled)?;
    Ok(RotorStage {
        mechanics: response,
        balance,
        endpoint,
    })
}
