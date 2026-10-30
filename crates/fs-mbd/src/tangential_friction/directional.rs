//! Directional Coulomb friction through actual two-body tangential mobility.
//!
//! This frozen-geometry, held-load midpoint impulse phase retains full body
//! inertia tensors and off-centre recoil. Normal impulses, pose integration,
//! changing load, thermal stores and event localization belong to their owners.
//! A coarse phase that assists endpoint slip or turns moving midpoint stiction
//! into a mirror velocity refuses; the trajectory owner must refine or locate
//! arrest. Exact stopping is accepted when the held phase resolves it.

/// Continuous held-load directional motion through resolved finite arrest.
pub mod trajectory;

use super::{RigidContactBinding, components, magnitude, vector};
use crate::{DynamicsError, MassProperties, PointKinematics, RigidBodyState, Vec3};
use fs_contact::interface_binding::thermal_friction::elliptic::EllipticFrictionMode;
use fs_exec::CancelGate;
use fs_math::det;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::tangential::elliptic::{
    EllipticCoulomb, EllipticCoulombState, SymmetricTangentialCompliance,
};
use fs_tribo::{
    FrictionRegime, TangentialSolveConfig, TangentialStep, TangentialStepInput, TriboError,
};

/// Original coincident contact points, fixed tangent gauge and held normal data.
#[derive(Clone, Copy, Debug)]
pub struct RigidDirectionalContact {
    pub binding: RigidContactBinding,
    pub frame: TangentFrame,
    pub normal_load_n: f64,
    pub contact_temperature_k: Option<f64>,
    pub nominal_pressure_pa: Option<f64>,
}

#[derive(Clone, Copy, Debug)]
pub struct RigidDirectionalControls {
    pub maximum_bodies: usize,
    pub geometry_tolerance_m: f64,
    pub velocity_tolerance_m_s: f64,
    pub linear_momentum_tolerance_ns: f64,
    pub angular_momentum_tolerance_nms: f64,
    pub energy_tolerance_j: f64,
    pub contact: TangentialSolveConfig,
}

/// Independently reconstructed native balances, in physical units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidDirectionalAudit {
    pub body_linear_momentum_residual_ns: [Vec3; 2],
    pub body_angular_momentum_residual_nms: [Vec3; 2],
    pub linear_momentum_change_ns: Vec3,
    pub angular_momentum_change_nms: Vec3,
    pub endpoint_velocity_residual_m_s: f64,
    pub work_energy_residual_j: f64,
    pub total_energy_residual_j: f64,
    pub contact_work_residual_j: f64,
    pub maximum_assisting_work_j: f64,
}

#[derive(Clone, Debug)]
pub struct RigidDirectionalReport<S> {
    pub contact: TangentialStep<S>,
    pub mobility_per_kg: [[f64; 2]; 2],
    pub compliance_m_per_n: [[f64; 2]; 2],
    pub slip_velocity_before_m_s: [f64; 2],
    pub slip_velocity_after_m_s: [f64; 2],
    /// Actual normal recoil is exposed for the surrounding normal owner.
    pub normal_velocity_before_m_s: f64,
    pub normal_velocity_after_m_s: f64,
    pub actual_midpoint_displacement_m: [f64; 2],
    pub force_impulse_on_a_world_ns: Vec3,
    /// Evaluated from momentum increments and mean velocities, without
    /// subtracting two potentially much larger absolute kinetic energies.
    pub kinetic_energy_change_j: f64,
    pub midpoint_impulse_work_j: f64,
    pub audit: RigidDirectionalAudit,
}

#[derive(Debug)]
pub enum RigidDirectionalError<E = TriboError> {
    Invalid(&'static str),
    Dynamics(DynamicsError),
    Compliance(TriboError),
    Source(E),
    Cancelled,
    UnresolvedSlipEvent { assisting_work_j: f64 },
    Audit(RigidDirectionalAudit),
}
impl<E> From<DynamicsError> for RigidDirectionalError<E> {
    fn from(e: DynamicsError) -> Self {
        Self::Dynamics(e)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for RigidDirectionalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid directional friction: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for RigidDirectionalError<E> {}

/// Original native directional law with exactly the retained material history.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_directional(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidDirectionalContact,
    law: &EllipticCoulomb,
    history: &mut EllipticCoulombState,
    dt_s: f64,
    controls: RigidDirectionalControls,
    gate: &CancelGate,
) -> Result<RigidDirectionalReport<EllipticCoulombState>, RigidDirectionalError> {
    apply_rigid_directional_with(
        bodies,
        properties,
        contact,
        history,
        dt_s,
        controls,
        |state, input, compliance, mode| match mode {
            EllipticFrictionMode::StaticOrSliding => {
                law.advance_with_compliance(state, input, compliance)
            }
            EllipticFrictionMode::Kinetic => {
                law.advance_kinetic_with_compliance(state, input, compliance)
            }
        },
        gate,
    )
}

/// Resolve the same mechanical phase with a retained source-law owner.
///
/// The callback must use the supplied original history, free displacement and
/// matrix compliance. A law owning speed/path coordinates must close its query
/// against the actual solved motion inside the callback; sampling free motion
/// as physical slip is invalid. The callback owns its source/domain/work budget.
/// Exactly two bodies and the accepted history are staged and published together.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn apply_rigid_directional_with<S: Clone, E>(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidDirectionalContact,
    history: &mut S,
    dt_s: f64,
    controls: RigidDirectionalControls,
    resolve: impl FnOnce(
        &S,
        &TangentialStepInput,
        SymmetricTangentialCompliance,
        EllipticFrictionMode,
    ) -> Result<TangentialStep<S>, E>,
    gate: &CancelGate,
) -> Result<RigidDirectionalReport<S>, RigidDirectionalError<E>> {
    poll(gate)?;
    let binding = contact.binding;
    if bodies.len() != properties.len()
        || bodies.len() > controls.maximum_bodies
        || binding.body_a >= bodies.len()
        || binding.body_b >= bodies.len()
        || binding.body_a == binding.body_b
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || !contact.normal_load_n.is_finite()
        || contact.normal_load_n < 0.0
        || [
            controls.geometry_tolerance_m,
            controls.velocity_tolerance_m_s,
            controls.linear_momentum_tolerance_ns,
            controls.angular_momentum_tolerance_nms,
            controls.energy_tolerance_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err(RigidDirectionalError::Invalid(
            "bounded original bodies, distinct contact pair and finite physical gates required",
        ));
    }
    let before = [bodies[binding.body_a], bodies[binding.body_b]];
    let properties = [properties[binding.body_a], properties[binding.body_b]];
    let arms = [binding.arm_a_body, binding.arm_b_body];
    let axes = [
        vector(contact.frame.longitudinal()),
        vector(contact.frame.lateral()),
        vector(contact.frame.normal()),
    ];
    let old = points(before, properties, arms)?;
    if magnitude(old[0].point_world.sub(old[1].point_world)) > controls.geometry_tolerance_m {
        return Err(RigidDirectionalError::Invalid(
            "original contact points do not coincide",
        ));
    }
    let mobility = response(before, properties, arms, axes)?;
    let scale = finite(0.5 * dt_s * dt_s)?;
    if scale == 0.0 {
        return Err(RigidDirectionalError::Invalid(
            "unrepresentable midpoint compliance",
        ));
    }
    let matrix = mobility.map(|row| row.map(|x| x * scale));
    let compliance =
        SymmetricTangentialCompliance::new(matrix).map_err(RigidDirectionalError::Compliance)?;
    let relative_before = old[0].point_velocity_world.sub(old[1].point_velocity_world);
    let initial = [relative_before.dot(axes[0]), relative_before.dot(axes[1])];
    let input = TangentialStepInput {
        normal_n: contact.normal_load_n,
        slip_increment_m: [finite(dt_s * initial[0])?, finite(dt_s * initial[1])?],
        dt_s,
        compliance_m_per_n: 0.0,
        contact_temperature_k: contact.contact_temperature_k,
        nominal_pressure_pa: contact.nominal_pressure_pa,
        solve: controls.contact,
    };
    let mode = if initial == [0.0; 2] {
        EllipticFrictionMode::StaticOrSliding
    } else {
        EllipticFrictionMode::Kinetic
    };
    let resolved =
        resolve(history, &input, compliance, mode).map_err(RigidDirectionalError::Source)?;
    poll(gate)?;
    // Native traction is positive work INTO the contact, hence acts negatively on A.
    let force = axes[0]
        .scale(-resolved.traction_n[0])
        .add(axes[1].scale(-resolved.traction_n[1]));
    let impulse = force.scale(dt_s);
    let mut after = before;
    for (i, sign) in [(0, 1.0), (1, -1.0)] {
        after[i] = before[i]
            .apply_impulse_at_body_point(properties[i], arms[i], impulse.scale(sign))?
            .state_after;
    }
    let new = points(after, properties, arms)?;
    let relative_after = new[0].point_velocity_world.sub(new[1].point_velocity_world);
    let final_rate = [relative_after.dot(axes[0]), relative_after.dot(axes[1])];
    let midpoint = relative_before.scale(0.5).add(relative_after.scale(0.5));
    let displacement = [
        finite(dt_s * midpoint.dot(axes[0]))?,
        finite(dt_s * midpoint.dot(axes[1]))?,
    ];
    let midpoint_work = finite(impulse.dot(midpoint))?;
    let mut kinetic = 0.0;
    let mut linear = Vec3::ZERO;
    let mut angular = Vec3::ZERO;
    let mut body_linear = [Vec3::ZERO; 2];
    let mut body_angular = [Vec3::ZERO; 2];
    for (i, sign) in [(0, 1.0), (1, -1.0)] {
        let dp = after[i]
            .linear_momentum_world()
            .sub(before[i].linear_momentum_world());
        let dl = after[i]
            .angular_momentum_body()
            .sub(before[i].angular_momentum_body());
        let mean_v = old[i]
            .center_of_mass_velocity_world
            .scale(0.5)
            .add(new[i].center_of_mass_velocity_world.scale(0.5));
        let mean_omega = old[i]
            .angular_velocity_body
            .scale(0.5)
            .add(new[i].angular_velocity_body.scale(0.5));
        kinetic = finite(kinetic + finite(dp.dot(mean_v))? + finite(dl.dot(mean_omega))?)?;
        let expected = impulse.scale(sign);
        body_linear[i] = dp.sub(expected);
        body_angular[i] = dl.sub(
            arms[i].cross(
                before[i]
                    .pose()
                    .orientation()
                    .rotate_world_to_body(expected),
            ),
        );
        linear = linear.add(dp);
        angular = angular
            .add(before[i].pose().position_world().cross(dp))
            .add(before[i].pose().orientation().rotate_body_to_world(dl));
    }
    let velocity_residual = finite(det::hypot(
        final_rate[0] - (2.0 * resolved.slip_increment_m[0] / dt_s - initial[0]),
        final_rate[1] - (2.0 * resolved.slip_increment_m[1] / dt_s - initial[1]),
    ))?;
    let assisting = finite(
        dt_s * force
            .dot(relative_before)
            .max(force.dot(relative_after))
            .max(0.0),
    )?;
    let audit = RigidDirectionalAudit {
        body_linear_momentum_residual_ns: body_linear,
        body_angular_momentum_residual_nms: body_angular,
        linear_momentum_change_ns: linear,
        angular_momentum_change_nms: angular,
        endpoint_velocity_residual_m_s: velocity_residual,
        work_energy_residual_j: finite(kinetic - midpoint_work)?,
        total_energy_residual_j: finite(
            kinetic + resolved.stored_change_j + resolved.dissipated_j,
        )?,
        contact_work_residual_j: finite(
            midpoint_work + resolved.stored_change_j + resolved.dissipated_j,
        )?,
        maximum_assisting_work_j: assisting,
    };
    if [
        audit.work_energy_residual_j,
        audit.total_energy_residual_j,
        audit.contact_work_residual_j,
    ]
    .iter()
    .any(|x| x.abs() > controls.energy_tolerance_j)
        || !resolved.dissipated_j.is_finite()
        || resolved.dissipated_j < 0.0
        || velocity_residual > controls.velocity_tolerance_m_s
        || body_linear
            .into_iter()
            .chain([linear])
            .any(|x| !x.is_finite() || magnitude(x) > controls.linear_momentum_tolerance_ns)
        || body_angular
            .into_iter()
            .chain([angular])
            .any(|x| !x.is_finite() || magnitude(x) > controls.angular_momentum_tolerance_nms)
    {
        return Err(RigidDirectionalError::Audit(audit));
    }
    if assisting > controls.energy_tolerance_j
        || initial != [0.0; 2] && resolved.regime == FrictionRegime::Sticking
    {
        return Err(RigidDirectionalError::UnresolvedSlipEvent {
            assisting_work_j: assisting,
        });
    }
    let report = RigidDirectionalReport {
        contact: resolved,
        mobility_per_kg: mobility,
        compliance_m_per_n: matrix,
        slip_velocity_before_m_s: initial,
        slip_velocity_after_m_s: final_rate,
        normal_velocity_before_m_s: relative_before.dot(axes[2]),
        normal_velocity_after_m_s: relative_after.dot(axes[2]),
        actual_midpoint_displacement_m: displacement,
        force_impulse_on_a_world_ns: impulse,
        kinetic_energy_change_j: kinetic,
        midpoint_impulse_work_j: midpoint_work,
        audit,
    };
    let next_history = report.contact.state.clone();
    poll(gate)?;
    bodies[binding.body_a] = after[0];
    bodies[binding.body_b] = after[1];
    *history = next_history;
    Ok(report)
}

fn poll<E>(gate: &CancelGate) -> Result<(), RigidDirectionalError<E>> {
    if gate.is_requested() {
        Err(RigidDirectionalError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite<E>(x: f64) -> Result<f64, RigidDirectionalError<E>> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(RigidDirectionalError::Invalid(
            "unrepresentable directional rigid-body response",
        ))
    }
}
fn points(
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    arms: [Vec3; 2],
) -> Result<[PointKinematics; 2], DynamicsError> {
    Ok([
        bodies[0].point_kinematics(properties[0], arms[0])?,
        bodies[1].point_kinematics(properties[1], arms[1])?,
    ])
}
fn response<E>(
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    arms: [Vec3; 2],
    axes: [Vec3; 3],
) -> Result<[[f64; 2]; 2], RigidDirectionalError<E>> {
    let mut factor = [[0.0; 12]; 2];
    for (body, sign) in [(0, 1.0), (1, -1.0)] {
        let inverse_mass_root = det::sqrt(properties[body].mass()).recip();
        for axis in 0..2 {
            let direction = bodies[body]
                .pose()
                .orientation()
                .rotate_world_to_body(axes[axis]);
            let translation = components(axes[axis].scale(sign * inverse_mass_root));
            let rotation = components(
                properties[body].angular_mobility_factor(arms[body].cross(direction).scale(sign)),
            );
            for j in 0..3 {
                factor[axis][6 * body + j] = translation[j];
                factor[axis][6 * body + 3 + j] = rotation[j];
            }
        }
    }
    let mut matrix = [[0.0; 2]; 2];
    for i in 0..2 {
        for j in 0..=i {
            matrix[i][j] = finite(factor[i].iter().zip(&factor[j]).map(|(a, b)| a * b).sum())?;
            matrix[j][i] = matrix[i][j];
        }
    }
    Ok(matrix)
}
