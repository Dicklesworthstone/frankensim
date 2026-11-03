//! Shared original rigid force kick and post-contact kinetic drift.
//!
//! Both elastic and retained-history solid trajectories use these same native
//! rigid states and phase accounts. Callers poll cancellation per body and own
//! phase tolerances, interval refinement and complete-candidate publication.

use super::super::{
    free::{rotation_difference, rotational_change},
    norm,
};
use super::FiniteSolidSustainedRigidStep;
use fs_mbd::{
    DynamicsError, Gravity, MassProperties, MomentumImpulseReceipt, RigidBodyIntegrator,
    RigidBodyState, Vec3, Wrench,
};

pub(crate) fn rigid_force_kick(
    old: RigidBodyState,
    property: MassProperties,
    wrench: Wrench,
    gravity: Gravity,
    duration_s: f64,
) -> Result<MomentumImpulseReceipt, DynamicsError> {
    let force = wrench
        .force_world
        .add(gravity.acceleration_world().scale(property.mass()));
    old.apply_momentum_impulse(
        property,
        force.scale(duration_s),
        wrench.torque_body.scale(duration_s),
    )
}

/// All quantities use the actual kick, contacted state and native drift. The
/// body-frame torque kick retains its original discrete work convention.
#[allow(clippy::too_many_arguments)]
pub(crate) fn rigid_contact_drift(
    kick: MomentumImpulseReceipt,
    contacted: RigidBodyState,
    property: MassProperties,
    wrench: Wrench,
    gravity: Gravity,
    duration_s: f64,
    origin: Vec3,
) -> Result<FiniteSolidSustainedRigidStep, DynamicsError> {
    let drift = RigidBodyIntegrator::new(Gravity::ZERO).step(
        contacted,
        property,
        Wrench::ZERO,
        duration_s,
    )?;
    let next = drift.state_after;
    let dp = next
        .linear_momentum_world()
        .sub(contacted.linear_momentum_world());
    let p_average = contacted
        .linear_momentum_world()
        .scale(0.5)
        .add(next.linear_momentum_world().scale(0.5));
    let dx = next
        .pose()
        .position_world()
        .sub(contacted.pose().position_world());
    let r_average = contacted
        .pose()
        .position_world()
        .sub(origin)
        .scale(0.5)
        .add(next.pose().position_world().sub(origin).scale(0.5));
    let dl = next
        .angular_momentum_body()
        .sub(contacted.angular_momentum_body());
    let drift_orbital = r_average.cross(dp).add(dx.cross(p_average));
    let drift_spin =
        contacted
            .pose()
            .orientation()
            .rotate_body_to_world(dl)
            .add(rotation_difference(
                contacted.pose().orientation(),
                next.pose().orientation(),
                next.angular_momentum_body(),
            ));
    let drift_angular = drift_orbital.add(drift_spin);
    let kinetic_drift_defect = p_average.scale(property.mass().recip()).dot(dp)
        + rotational_change(
            property,
            contacted.angular_momentum_body(),
            next.angular_momentum_body(),
        )?;
    let force = wrench
        .force_world
        .add(gravity.acceleration_world().scale(property.mass()));
    let force_displacement_work = force.dot(
        next.pose()
            .position_world()
            .sub(kick.state_before.pose().position_world()),
    );
    let velocity_average = kick
        .state_before
        .center_of_mass_velocity_world(property)?
        .scale(0.5)
        .add(
            kick.state_after
                .center_of_mass_velocity_world(property)?
                .scale(0.5),
        );
    let omega_average = property
        .angular_velocity_body_checked(kick.state_before.angular_momentum_body())?
        .scale(0.5)
        .add(
            property
                .angular_velocity_body_checked(kick.state_after.angular_momentum_body())?
                .scale(0.5),
        );
    let force_kick_work = kick.linear_impulse_world.dot(velocity_average);
    let torque_kick_work = kick.angular_impulse_body.dot(omega_average);
    let force_splitting = force_kick_work - force_displacement_work;
    let kick_orbital_residual = kick
        .state_before
        .pose()
        .position_world()
        .sub(origin)
        .cross(kick.linear_residual_world);
    let kick_spin_residual = kick
        .state_before
        .pose()
        .orientation()
        .rotate_body_to_world(kick.angular_residual_body);
    let kick_angular_residual = kick_orbital_residual.add(kick_spin_residual);
    let absolute_angular_residual = norm(kick_orbital_residual)
        + norm(kick_spin_residual)
        + norm(drift_orbital)
        + norm(drift_spin);
    Ok(FiniteSolidSustainedRigidStep {
        kick,
        drift,
        force_impulse_work_j: force_kick_work,
        torque_impulse_work_j: torque_kick_work,
        force_displacement_work_j: force_displacement_work,
        force_splitting_work_j: force_splitting,
        kinetic_drift_defect_j: kinetic_drift_defect,
        drift_linear_increment_kg_m_s: dp,
        drift_angular_increment_kg_m2_s: drift_angular,
        kick_angular_residual_kg_m2_s: kick_angular_residual,
        absolute_angular_residual_kg_m2_s: absolute_angular_residual,
    })
}

impl FiniteSolidSustainedRigidStep {
    /// Original independent phase metrics, before any whole-system reduction.
    pub(crate) fn phase_residuals(&self) -> (f64, f64, f64) {
        (
            norm(self.kick.linear_residual_world) + norm(self.drift_linear_increment_kg_m_s),
            self.absolute_angular_residual_kg_m2_s,
            self.kick.work.work_energy_residual_j.abs()
                + (self.kick.work.impulse_work_j
                    - (self.force_impulse_work_j + self.torque_impulse_work_j))
                    .abs(),
        )
    }
}
