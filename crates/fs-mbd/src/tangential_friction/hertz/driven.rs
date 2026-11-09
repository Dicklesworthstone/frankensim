//! Powered Hertz shear/spin on the original finite rigid-body pair.
//!
//! Held world forces and body torques drive the original full contact mobility.
//! The pressure-resolved law owns static admission, breakaway and arrest. Its
//! integrated impulses and the external impulses update each original momentum
//! once. Source work includes common-body motion as well as relative slip/spin;
//! only the law's physical traction work becomes frictional heat.
//!
//! Poses, arms, inertia, normal load and patch radius stay fixed over the phase.
//! A whole-path normal-rate bound refuses untracked opening. Curved trajectories
//! have the timestep accuracy declared by the original driven-law report.

use super::trajectory::PreparedRigidHertzTrajectory;
use super::transient::normal_response_row;
use super::*;
use crate::{RigidKineticEnergyChange, Wrench};
use fs_tribo::sliding_patch::contact::driven::{
    HertzDrivenError, HertzDrivenInput, HertzDrivenReport,
};

pub use fs_tribo::sliding_patch::contact::driven::{HertzDrivenControls, HertzDrivenMode};

/// Independent checks of the actual original bodies and applied work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidHertzDrivenAudit {
    pub body_linear_momentum_residual_ns: [Vec3; 2],
    pub body_angular_momentum_residual_nms: [Vec3; 2],
    pub linear_momentum_residual_ns: Vec3,
    pub angular_momentum_residual_nms: Vec3,
    pub endpoint_velocity_residual_m_s: f64,
    pub endpoint_spin_residual_rad_s: f64,
    pub driving_velocity_residual_m_s: f64,
    /// Conservative bound throughout the entire held-geometry phase.
    pub normal_velocity_bound_m_s: f64,
    pub relative_drive_work_residual_j: f64,
    pub relative_kinetic_residual_j: f64,
    pub work_energy_residual_j: f64,
    /// Actual full-body kinetic increment + physical heat - external work.
    pub total_energy_residual_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzDrivenReport {
    pub driven: HertzDrivenReport,
    pub scaled_mobility_per_kg: [[f64; 3]; 3],
    pub slip_velocity_m_s: [[f64; 2]; 2],
    pub spin_rad_s: [f64; 2],
    pub normal_velocity_m_s: [f64; 2],
    /// Contact resultants on A; B receives their opposites.
    pub force_impulse_world_ns: Vec3,
    pub spin_angular_impulse_world_nms: Vec3,
    pub external: [Wrench; 2],
    /// Actual original-body storage change, reduced before rounding.
    pub kinetic_change_j: f64,
    /// External force/torque work over the physical path, including common motion.
    pub external_work_j: f64,
    /// Summed resultant impulses paired with actual native endpoint mean rates.
    /// This is an algebraic work check, not an external-work quadrature.
    pub impulse_work_j: f64,
    pub dissipated_energy_j: f64,
    pub audit: RigidHertzDrivenAudit,
}

/// Staged body pair and retained physical mode for joint thermal publication.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzDrivenCandidate {
    pub bodies_after: [RigidBodyState; 2],
    pub report: RigidHertzDrivenReport,
}

#[derive(Clone, Debug)]
pub struct PreparedRigidHertzDriven {
    pair: PreparedRigidHertzTrajectory,
    external: [Wrench; 2],
    mode: HertzDrivenMode,
    drive: [f64; 3],
    normal_acceleration_m_s2: f64,
    normal_response: Vec3,
    driven: HertzDrivenControls,
}

#[derive(Debug)]
pub enum RigidHertzDrivenError {
    Invalid(&'static str),
    Mechanics(RigidHertzError),
    Driven(HertzDrivenError),
    Audit(RigidHertzDrivenAudit),
}
impl From<RigidHertzError> for RigidHertzDrivenError {
    fn from(error: RigidHertzError) -> Self {
        Self::Mechanics(error)
    }
}
impl From<DynamicsError> for RigidHertzDrivenError {
    fn from(error: DynamicsError) -> Self {
        Self::Mechanics(error.into())
    }
}
impl From<RigidFrictionError> for RigidHertzDrivenError {
    fn from(error: RigidFrictionError) -> Self {
        Self::Mechanics(error.into())
    }
}
impl From<HertzDrivenError> for RigidHertzDrivenError {
    fn from(error: HertzDrivenError) -> Self {
        Self::Driven(error)
    }
}
impl core::fmt::Display for RigidHertzDrivenError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid driven Hertz trajectory: {self:?}")
    }
}
impl std::error::Error for RigidHertzDrivenError {}

/// Publish both original bodies and the retained mode after all independent gates.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_hertz_driven(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidHertzContact<'_>,
    external: [Wrench; 2],
    mode: &mut HertzDrivenMode,
    dt_s: f64,
    controls: RigidHertzControls,
    driven: HertzDrivenControls,
    gate: &CancelGate,
) -> Result<RigidHertzDrivenReport, RigidHertzDrivenError> {
    let candidate = PreparedRigidHertzDriven::new(
        bodies, properties, contact, external, *mode, dt_s, controls, driven, gate,
    )?
    .advance(contact.law, gate)?;
    poll(gate)?;
    bodies[contact.binding.body_a] = candidate.bodies_after[0];
    bodies[contact.binding.body_b] = candidate.bodies_after[1];
    *mode = candidate.report.driven.mode;
    Ok(candidate.report)
}

impl PreparedRigidHertzDriven {
    /// Capture actual geometry, inertia and held external loads once per trial set.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contact: RigidHertzContact<'_>,
        external: [Wrench; 2],
        mode: HertzDrivenMode,
        dt_s: f64,
        controls: RigidHertzControls,
        driven: HertzDrivenControls,
        gate: &CancelGate,
    ) -> Result<Self, RigidHertzDrivenError> {
        let pair =
            PreparedRigidHertzTrajectory::new(bodies, properties, contact, dt_s, controls, gate)?;
        let mut linear_acceleration = [Vec3::ZERO; 2];
        let mut angular_acceleration = [Vec3::ZERO; 2];
        for i in 0..2 {
            poll(gate)?;
            external[i].validate()?;
            let orientation = pair.before[i].pose().orientation();
            angular_acceleration[i] = orientation.rotate_body_to_world(
                pair.properties[i].angular_velocity_body_checked(external[i].torque_body)?,
            );
            linear_acceleration[i] = external[i]
                .force_world
                .scale(pair.properties[i].mass().recip());
        }
        // Cancel common COM acceleration before adding either rotational point
        // contribution, as the original relative-velocity owner does for rates.
        let mut relative = [Dd::ZERO; 3];
        add_components(&mut relative, linear_acceleration[0]);
        add_components(&mut relative, linear_acceleration[1].scale(-1.0));
        for (i, sign) in [(0, 1.0), (1, -1.0)] {
            add_components(
                &mut relative,
                angular_acceleration[i]
                    .cross(pair.points_before[i].arm_world)
                    .scale(sign),
            );
        }
        let relative = rounded_components(relative);
        let spin = angular_acceleration[0]
            .sub(angular_acceleration[1])
            .dot(pair.axes[2]);
        let drive = [
            finite(relative.dot(pair.axes[0]))?,
            finite(relative.dot(pair.axes[1]))?,
            finite(pair.radius * spin)?,
        ];
        let normal_acceleration_m_s2 = finite(relative.dot(pair.axes[2]))?;
        let normal_response = normal_response_row(&pair)?;
        Ok(Self {
            pair,
            external,
            mode,
            drive,
            normal_acceleration_m_s2,
            normal_response,
            driven,
        })
    }

    #[must_use]
    pub fn contact_points_world(&self) -> [Vec3; 2] {
        self.pair.contact_points_world()
    }

    /// Retry original coefficients from the same bodies, mode and external drive.
    #[allow(clippy::too_many_lines)]
    pub fn advance(
        &self,
        law: &HertzWrenchContact,
        gate: &CancelGate,
    ) -> Result<RigidHertzDrivenCandidate, RigidHertzDrivenError> {
        poll(gate)?;
        let pair = &self.pair;
        pair.admit_law(law)?;
        let controls = pair.controls;
        let axes = pair.axes;
        let dt = pair.dt_s;
        let initial = relative_point(pair.points_before)?;
        let spin = relative_spin(pair.points_before).dot(axes[2]);
        let velocity = [initial.dot(axes[0]), initial.dot(axes[1])];
        let capacity = law.parameters().normal_force_n * law.parameters().static_mu;
        let wrench_bound = capacity.hypot(3.0 * core::f64::consts::PI / 16.0 * capacity);
        let normal_bound = initial.dot(axes[2]).abs()
            + dt * (self.normal_acceleration_m_s2.abs()
                + magnitude(self.normal_response) * wrench_bound);
        if !normal_bound.is_finite() || normal_bound > controls.velocity_tolerance_m_s {
            return Err(RigidHertzDrivenError::Invalid(
                "driven traction cannot preserve held closed normal motion",
            ));
        }
        let path = law.advance_driven(
            &HertzDrivenInput {
                velocity_m_s: velocity,
                spin_rate_rad_s: spin,
                scaled_mobility_per_kg: pair.mobility,
                driving_acceleration_m_s2: [self.drive[0], self.drive[1]],
                driving_spin_acceleration_rad_s2: self.drive[2] / pair.radius,
                dt_s: dt,
                mode: self.mode,
            },
            controls.contact,
            self.driven,
            || gate.is_requested(),
        )?;
        let impulse = axes[0]
            .scale(path.force_impulse_ns[0])
            .add(axes[1].scale(path.force_impulse_ns[1]));
        let couple = axes[2].scale(path.spin_angular_impulse_nms);
        let mut after = pair.before;
        let mut kinetic = RigidKineticEnergyChange::new();
        let mut free_work = Dd::ZERO;
        let mut total_linear = [Dd::ZERO; 3];
        let mut total_angular = [Dd::ZERO; 3];
        let mut body_linear = [Vec3::ZERO; 2];
        let mut body_angular = [Vec3::ZERO; 2];
        let mut requested = [(Vec3::ZERO, Vec3::ZERO); 2];
        for (i, sign) in [(0, 1.0), (1, -1.0)] {
            poll(gate)?;
            let before = pair.before[i];
            let property = pair.properties[i];
            let orientation = before.pose().orientation();
            let external_linear = self.external[i].force_world.scale(dt);
            let external_angular = self.external[i].torque_body.scale(dt);
            let contact_linear = impulse.scale(sign);
            let contact_angular = pair.arms[i]
                .cross(orientation.rotate_world_to_body(contact_linear))
                .add(orientation.rotate_world_to_body(couple.scale(sign)));
            let linear = external_linear.add(contact_linear);
            let angular = external_angular.add(contact_angular);
            requested[i] = (linear, angular);
            after[i] = before
                .apply_momentum_impulse(property, linear, angular)?
                .state_after;
            kinetic.add_body(property, before, after[i])?;
            // Independently integrate the free external power. Do not recover
            // it from a rounded free endpoint: an unrepresentable source kick
            // must remain source work, not disappear from the energy account.
            let free_linear_increment = external_linear.scale(property.mass().recip());
            let free_angular_increment =
                property.angular_velocity_body_checked(external_angular)?;
            free_work = free_work
                + dot_wide(
                    components(external_linear),
                    components(pair.points_before[i].center_of_mass_velocity_world),
                )
                + dot_wide(
                    components(external_linear),
                    components(free_linear_increment),
                ) * Dd::from_f64(0.5)
                + dot_wide(
                    components(external_angular),
                    components(pair.points_before[i].angular_velocity_body),
                )
                + dot_wide(
                    components(external_angular),
                    components(free_angular_increment),
                ) * Dd::from_f64(0.5);
            let linear_change = after[i]
                .linear_momentum_world()
                .sub(before.linear_momentum_world());
            let angular_change = after[i]
                .angular_momentum_body()
                .sub(before.angular_momentum_body());
            body_linear[i] = linear_change.sub(linear);
            body_angular[i] = angular_change.sub(angular);
            let remaining_linear = linear_change.sub(external_linear);
            let remaining_angular = angular_change.sub(external_angular);
            add_components(&mut total_linear, remaining_linear);
            add_components(
                &mut total_angular,
                before.pose().position_world().cross(remaining_linear),
            );
            add_components(
                &mut total_angular,
                orientation.rotate_body_to_world(remaining_angular),
            );
        }
        let kinetic = kinetic.finish_retained()?;
        let new = points(after, pair.properties, pair.arms)?;
        let endpoint = relative_point(new)?;
        let endpoint_spin = relative_spin(new).dot(axes[2]);
        let q0 = [velocity[0], velocity[1], pair.radius * spin];
        let travel = [
            path.slip_increment_m[0],
            path.slip_increment_m[1],
            pair.radius * path.spin_increment_rad,
        ];
        let free_travel = std::array::from_fn(|i| q0[i] * dt + 0.5 * self.drive[i] * dt * dt);
        let drive_work = dot_wide(path.driving_scaled_wrench_n, travel);
        let complementary_work = free_work - dot_wide(path.driving_scaled_wrench_n, free_travel);
        let external_work = complementary_work + drive_work;
        let mut impulse_work = Dd::ZERO;
        for i in 0..2 {
            let mean_linear = pair.points_before[i]
                .center_of_mass_velocity_world
                .scale(0.5)
                .add(new[i].center_of_mass_velocity_world.scale(0.5));
            let mean_angular = pair.points_before[i]
                .angular_velocity_body
                .scale(0.5)
                .add(new[i].angular_velocity_body.scale(0.5));
            impulse_work = impulse_work
                + dot_wide(components(requested[i].0), components(mean_linear))
                + dot_wide(components(requested[i].1), components(mean_angular));
        }
        let driving_residual: [f64; 3] = std::array::from_fn(|i| {
            (dot_wide(pair.mobility[i], path.driving_scaled_wrench_n) - Dd::from_f64(self.drive[i]))
                .to_f64()
                * dt
        });
        let audit = RigidHertzDrivenAudit {
            body_linear_momentum_residual_ns: body_linear,
            body_angular_momentum_residual_nms: body_angular,
            linear_momentum_residual_ns: rounded_components(total_linear),
            angular_momentum_residual_nms: rounded_components(total_angular),
            endpoint_velocity_residual_m_s: (endpoint.dot(axes[0]) - path.velocity_after_m_s[0])
                .hypot(endpoint.dot(axes[1]) - path.velocity_after_m_s[1]),
            endpoint_spin_residual_rad_s: (endpoint_spin - path.spin_rate_after_rad_s).abs(),
            driving_velocity_residual_m_s: magnitude(vector(driving_residual)),
            normal_velocity_bound_m_s: normal_bound.max(endpoint.dot(axes[2]).abs()),
            relative_drive_work_residual_j: (drive_work - Dd::from_f64(path.drive_work_j)).to_f64(),
            relative_kinetic_residual_j: (kinetic
                - complementary_work
                - Dd::from_f64(path.kinetic_energy_change_j))
            .to_f64(),
            work_energy_residual_j: (kinetic - impulse_work).to_f64(),
            total_energy_residual_j: (kinetic + Dd::from_f64(path.dissipated_energy_j)
                - external_work)
                .to_f64(),
        };
        if [
            (
                audit.endpoint_velocity_residual_m_s,
                controls.velocity_tolerance_m_s,
            ),
            (
                audit.endpoint_spin_residual_rad_s,
                controls.angular_velocity_tolerance_rad_s,
            ),
            (
                audit.driving_velocity_residual_m_s,
                controls.velocity_tolerance_m_s,
            ),
            (
                audit.normal_velocity_bound_m_s,
                controls.velocity_tolerance_m_s,
            ),
        ]
        .into_iter()
        .any(|(v, tolerance)| !v.is_finite() || v > tolerance)
            || [
                audit.relative_drive_work_residual_j,
                audit.relative_kinetic_residual_j,
                audit.work_energy_residual_j,
                audit.total_energy_residual_j,
            ]
            .into_iter()
            .any(|v| !v.is_finite() || v.abs() > controls.energy_tolerance_j)
            || body_linear
                .into_iter()
                .chain([audit.linear_momentum_residual_ns])
                .any(|v| !v.is_finite() || magnitude(v) > controls.linear_momentum_tolerance_ns)
            || body_angular
                .into_iter()
                .chain([audit.angular_momentum_residual_nms])
                .any(|v| !v.is_finite() || magnitude(v) > controls.angular_momentum_tolerance_nms)
            || !external_work.to_f64().is_finite()
        {
            return Err(RigidHertzDrivenError::Audit(audit));
        }
        let report = RigidHertzDrivenReport {
            scaled_mobility_per_kg: pair.mobility,
            slip_velocity_m_s: [velocity, [endpoint.dot(axes[0]), endpoint.dot(axes[1])]],
            spin_rad_s: [spin, endpoint_spin],
            normal_velocity_m_s: [initial.dot(axes[2]), endpoint.dot(axes[2])],
            force_impulse_world_ns: impulse,
            spin_angular_impulse_world_nms: couple,
            external: self.external,
            kinetic_change_j: kinetic.to_f64(),
            external_work_j: external_work.to_f64(),
            impulse_work_j: impulse_work.to_f64(),
            dissipated_energy_j: path.dissipated_energy_j,
            audit,
            driven: path,
        };
        poll(gate)?;
        Ok(RigidHertzDrivenCandidate {
            bodies_after: after,
            report,
        })
    }
}

fn dot_wide(a: [f64; 3], b: [f64; 3]) -> Dd {
    (0..3).fold(Dd::ZERO, |sum, i| {
        sum + Dd::from_f64(a[i]) * Dd::from_f64(b[i])
    })
}
