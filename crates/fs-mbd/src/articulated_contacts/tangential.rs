//! Stateful tangential friction and elastic shear on articulated trees.
//!
//! One fixed tangent has the exact projected native mobility W = J M^-1 J^T.
//! Over a held-force interval its midpoint displacement response is
//! c = dt^2 W / 2. The original constitutive law solves
//! du = -dt v_relative_before - c F; its work-conjugate average effort applies
//! the impulse dt F to the SAME ABA owner. Actual endpoint velocities then
//! reconstruct du and close kinetic energy, reversible contact storage,
//! physical heat and translating-support work together.
//!
//! This is a frozen-configuration shear substep, not a normal contact solver
//! or a pose integrator. A prescribed normal load does not enforce normal
//! kinematics. The caller owns normal evolution, geometric integration, and
//! the unchanged physical tangent/history association between accepted steps.
//! A dry sticking branch constrains average slip over the interval, not the
//! endpoint velocity; resolving transitions and intrastep reversals requires
//! the enclosing integrator's time refinement.

use super::{
    ArticulatedContactError, ArticulatedContactSystem, ContactCounterpart, EventResponse,
    EventWorkspace, LinkContactPoint, Prepared, finite, magnitude, poll, prepare_point,
    unit_normal, vector,
};
use crate::articulated::FreeFloatingImpulseResponse;
use fs_ga::Vec3;
use fs_tribo::{TangentialLaw, TangentialSolveConfig, TangentialStep, TangentialStepInput};

/// One scalar shear coordinate; its positive effort acts on `point`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ArticulatedTangentialContact {
    /// Original articulated link material point.
    pub point: LinkContactPoint,
    /// Opposite dynamic point or prescribed translating support. For a support,
    /// `point_world` must identify the actual common interface point here.
    pub other: ContactCounterpart,
    /// Nonzero physical tangent, normalized before use. Constitutive axis 0
    /// retains this direction; transverse response is not silently discarded.
    pub tangent_world: Vec3,
}

/// Held normal/environment inputs and an explicit physical interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ArticulatedTangentialDrive {
    /// Positive physical duration [s]; body and joint poses remain fixed.
    pub dt_s: f64,
    /// Prescribed compressive load [N]. The law owns opening/load restrictions.
    pub normal_force_n: f64,
    /// Explicit contact boundary [K], passed unchanged to the original law.
    pub contact_temperature_k: Option<f64>,
    /// Explicit nominal pressure [Pa], never inferred from normal force.
    pub nominal_pressure_pa: Option<f64>,
}

/// Independent geometry, physical work and numerical acceptance controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ArticulatedTangentialConfig {
    /// Original constitutive root/force controls, with no implicit law defaults.
    pub solve: TangentialSolveConfig,
    /// Largest caller-admitted frozen-configuration interval [s].
    pub maximum_step_s: f64,
    /// Bound on actual midpoint relative speed [m/s].
    pub maximum_slip_speed_m_s: f64,
    /// Common-point geometry tolerance [m]; points are never projected.
    pub geometry_tolerance_m: f64,
    /// Actual displacement versus constitutive-solve allowance [m].
    pub displacement_tolerance_m: f64,
    /// Actual velocity versus native mobility prediction allowance [m/s].
    pub velocity_tolerance_m_s: f64,
    /// Absolute work/storage closure allowance [J]. Native mechanics also
    /// retains its existing scale-aware roundoff allowance.
    pub energy_tolerance_j: f64,
}

/// Complete private mechanical/constitutive candidate for one accepted substep.
/// The caller publishes the checkpoint and all body responses together.
#[derive(Clone, Debug)]
pub struct ArticulatedTangentialResult<S> {
    /// Original ABA state updates, in input-system order.
    pub responses: Vec<FreeFloatingImpulseResponse>,
    /// Original law result, including retained memory, average effort and heat.
    pub contact: TangentialStep<S>,
    /// Actual scalar mobility [1/kg]; never a lumped link mass.
    pub mobility_kg_inverse: f64,
    /// Held-force midpoint displacement compliance [m/N].
    pub compliance_m_per_n: f64,
    /// Positive-body minus counterpart endpoint speed [m/s].
    pub relative_speed_after_m_s: f64,
    /// Counterpart-minus-body midpoint relative increment [m].
    pub slip_increment_m: f64,
    /// Actual increment minus the constitutive returned increment [m].
    pub compatibility_residual_m: f64,
    /// Actual reversible contact energy before/after [J].
    pub contact_energy_before_j: f64,
    pub contact_energy_after_j: f64,
    /// Actual body kinetic-energy change [J], which may be positive on unloading.
    pub kinetic_energy_change_j: f64,
    /// Signed work supplied by prescribed support translation [J].
    pub support_work_j: f64,
    /// Actual contact storage change + heat - average traction * actual slip [J].
    pub contact_work_residual_j: f64,
    /// Actual kinetic + contact storage change + heat - support work [J].
    pub energy_balance_residual_j: f64,
    /// Actual world momentum changes after external support impulse [N s, N m s].
    pub linear_momentum_residual: Vec3,
    pub angular_momentum_residual: Vec3,
    /// Actual versus projected native point velocity discrepancy [m/s].
    pub endpoint_velocity_residual_m_s: f64,
}

/// Advance one stateful scalar shear port over at most 128 articulated trees.
///
/// Input mechanics, law and checkpoint are borrowed. Any refusal or cancellation
/// returns no candidate; retrying the same originals is deterministic. Elastic
/// unloading can return stored contact energy to the bodies, while exported
/// physical heat must remain nonnegative. No impact loss is relabeled as heat.
///
/// Only one tangent is assembled, so this does not approximate a general 2-D
/// mobility by a scalar. Multiple coupled shear ports, finite-patch spin/rolling,
/// normal-load evolution, pose advancement and thermal allocation are outside
/// this operation. A larger integrator must account for its own normal/external
/// forces and geometric splitting error separately.
pub fn advance_articulated_tangential<L: TangentialLaw>(
    systems: &[ArticulatedContactSystem<'_>],
    contact: ArticulatedTangentialContact,
    law: &L,
    state: &L::State,
    drive: ArticulatedTangentialDrive,
    config: ArticulatedTangentialConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedTangentialResult<L::State>, ArticulatedContactError> {
    poll(&mut cancelled)?;
    admit(systems.len(), drive, config)?;
    let before_energy = finite(law.stored_energy_j(state))?;
    if before_energy < 0.0 {
        return Err(ArticulatedContactError::Invalid("negative contact storage"));
    }
    let mut workspace = EventWorkspace::new(systems, &mut cancelled)?;
    let axis = unit_normal(contact.tangent_world)?;
    let (first, world_point) = prepare_point(contact.point, axis, &workspace.before)?;
    let (second, other_point, plane_speed) = match contact.other {
        ContactCounterpart::Body(point) => {
            if point.system_index == contact.point.system_index
                && point.link_index == contact.point.link_index
            {
                return Err(ArticulatedContactError::Invalid(
                    "shear requires distinct links",
                ));
            }
            let (term, position) = prepare_point(point, axis.scale(-1.0), &workspace.before)?;
            (Some(term), position, 0.0)
        }
        ContactCounterpart::PrescribedPlane {
            point_world,
            velocity_world,
        } => {
            vector(point_world)?;
            vector(velocity_world)?;
            (None, point_world, finite(axis.dot(velocity_world))?)
        }
    };
    if magnitude(world_point - other_point) > config.geometry_tolerance_m {
        return Err(ArticulatedContactError::Invalid(
            "shear material points do not coincide",
        ));
    }
    let row = Prepared {
        first,
        second,
        axis_world: axis,
        world_point,
        plane_speed,
        minimum_speed: None,
    };
    let before_speed = finite(row.velocity(&workspace.before))?;
    let rows = [row];
    let mobility = workspace.delassus(systems, &rows, &mut cancelled)?[0];
    if !mobility.is_finite() || mobility < 0.0 {
        return Err(ArticulatedContactError::Invalid(
            "nonphysical shear mobility",
        ));
    }
    let compliance = finite(0.5 * drive.dt_s * drive.dt_s * mobility)?;
    let step = law.advance(
        state,
        &TangentialStepInput {
            normal_n: drive.normal_force_n,
            slip_increment_m: [finite(-drive.dt_s * before_speed)?, 0.0],
            dt_s: drive.dt_s,
            contact_temperature_k: drive.contact_temperature_k,
            nominal_pressure_pa: drive.nominal_pressure_pa,
            compliance_m_per_n: compliance,
            solve: config.solve,
        },
    )?;
    poll(&mut cancelled)?;
    if step.traction_n[1] != 0.0 || step.slip_increment_m[1] != 0.0 {
        return Err(ArticulatedContactError::Invalid(
            "scalar shear refuses transverse response",
        ));
    }
    let force = finite(step.traction_n[0])?;
    if force.abs() > config.solve.maximum_force_n || step.iterations > config.solve.max_iterations {
        return Err(ArticulatedContactError::Unresolved(
            "constitutive effort/work cap exceeded",
        ));
    }
    for residual in step.compatibility_residual_m {
        if finite(residual)?.abs() > config.displacement_tolerance_m {
            return Err(ArticulatedContactError::Unresolved(
                "constitutive shear compatibility",
            ));
        }
    }
    let after_energy = finite(law.stored_energy_j(&step.state))?;
    let stored_change = finite(after_energy - before_energy)?;
    let heat = finite(step.dissipated_j)?;
    if after_energy < 0.0 || heat < 0.0 {
        return Err(ArticulatedContactError::Invalid(
            "negative contact storage or heat",
        ));
    }
    for residual in [
        finite(stored_change - step.stored_change_j)?,
        finite(step.closure_residual_j)?,
    ] {
        if residual.abs() > config.energy_tolerance_j {
            return Err(ArticulatedContactError::Unresolved(
                "constitutive contact energy",
            ));
        }
    }
    let impulse = finite(drive.dt_s * force)?;
    let predicted_speed = finite(before_speed + mobility * impulse)?;
    let mechanics = workspace.apply(
        systems,
        &rows,
        EventResponse {
            impulses: &[impulse],
            predicted_speeds: &[predicted_speed],
            stored_energy_change: stored_change,
            dissipated_energy: heat,
            velocity_tolerance: config.velocity_tolerance_m_s,
            energy_tolerance: config.energy_tolerance_j,
        },
        &mut cancelled,
    )?;
    let after_speed = mechanics.speeds_after[0];
    let slip = finite(-drive.dt_s * (0.5 * before_speed + 0.5 * after_speed))?;
    let compatibility = finite(slip - step.slip_increment_m[0])?;
    let work_residual = finite(stored_change + heat - force * slip)?;
    if compatibility.abs() > config.displacement_tolerance_m
        || finite(slip / drive.dt_s)?.abs() > config.maximum_slip_speed_m_s
        || work_residual.abs() > config.energy_tolerance_j
    {
        return Err(ArticulatedContactError::Unresolved(
            "actual shear displacement/work",
        ));
    }
    poll(&mut cancelled)?;
    Ok(ArticulatedTangentialResult {
        responses: mechanics.responses,
        contact: step,
        mobility_kg_inverse: mobility,
        compliance_m_per_n: compliance,
        relative_speed_after_m_s: after_speed,
        slip_increment_m: slip,
        compatibility_residual_m: compatibility,
        contact_energy_before_j: before_energy,
        contact_energy_after_j: after_energy,
        kinetic_energy_change_j: mechanics.kinetic_energy_change,
        support_work_j: mechanics.plane_work,
        contact_work_residual_j: work_residual,
        energy_balance_residual_j: mechanics.energy_balance_residual,
        linear_momentum_residual: mechanics.linear_momentum_residual,
        angular_momentum_residual: mechanics.angular_momentum_residual,
        endpoint_velocity_residual_m_s: mechanics.endpoint_velocity_residual,
    })
}

fn admit(
    count: usize,
    drive: ArticulatedTangentialDrive,
    config: ArticulatedTangentialConfig,
) -> Result<(), ArticulatedContactError> {
    if count > 128
        || count == 0
        || !(1..=128).contains(&config.solve.max_iterations)
        || !drive.dt_s.is_finite()
        || drive.dt_s <= 0.0
        || drive.dt_s > config.maximum_step_s
        || !drive.normal_force_n.is_finite()
        || drive.normal_force_n < 0.0
        || drive
            .contact_temperature_k
            .is_some_and(|t| !t.is_finite() || t <= 0.0)
        || drive
            .nominal_pressure_pa
            .is_some_and(|p| !p.is_finite() || p < 0.0)
        || [
            config.solve.maximum_force_n,
            config.solve.force_absolute_tolerance_n,
            config.solve.force_relative_tolerance,
            config.maximum_step_s,
            config.maximum_slip_speed_m_s,
            config.geometry_tolerance_m,
            config.displacement_tolerance_m,
            config.velocity_tolerance_m_s,
            config.energy_tolerance_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        || config.solve.force_relative_tolerance >= 1.0
    {
        return Err(ArticulatedContactError::Invalid(
            "invalid articulated shear controls",
        ));
    }
    Ok(())
}
