//! Closed, deterministic reduced Euler-disc trajectory runner.
//!
//! This runner is deliberately a bounded numerical model, not a validation
//! claim.  It advances the production `fs_mbd` rigid-body state while retaining
//! separately evaluated contact, rolling, base, and gas wrenches.  The new
//! Euler adapter modules are exported by the crate, but this runner does not
//! yet consume their receipts: its channel laws remain deliberately reduced.

use core::fmt;

use crate::specimen::ResolvedDiscProfile;
use crate::{
    ContactDiscGeometry, ContactGeometry, contact_geometry, profile_contact_geometry,
    profile_state_at_ground_contact, state_at_ground_contact,
};
use fs_exec::Cx;
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

/// A factor set supplied by a campaign case, in SI units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoupledFactors {
    pub mass_kg: f64,
    pub radius_m: f64,
    pub thickness_m: f64,
    pub transverse_inertia_kg_m2: f64,
    pub axial_inertia_kg_m2: f64,
    pub gravity_m_per_s2: f64,
    /// Saturated gross-sliding Coulomb coefficient; no static-stick solve is performed.
    pub sliding_friction_coefficient: f64,
    pub rolling_resistance_m: f64,
    /// Kelvin-Voigt normal-contact stiffness [N/m].
    pub contact_stiffness_n_per_m: f64,
    /// Kelvin-Voigt normal-contact damping [N s/m].
    pub contact_damping_n_s_per_m: f64,
    /// Effective vertical inertia of the base mode [kg].
    pub base_effective_mass_kg: f64,
    pub base_stiffness_n_per_m: f64,
    pub base_damping_n_s_per_m: f64,
    pub gas_rotational_damping_n_m_s: f64,
    pub gas_translation_damping_n_s_per_m: f64,
}

/// Geometry-independent material and environment factors, in SI units.
///
/// Profile trajectories derive mass, center of mass, principal inertia, outer
/// radius, and thickness from one [`ResolvedDiscProfile`]. Keeping those
/// quantities out of this input prevents a caller from pairing one contact
/// shape with hand-entered inertia from another shape.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoupledChannelFactors {
    pub gravity_m_per_s2: f64,
    pub sliding_friction_coefficient: f64,
    pub rolling_resistance_m: f64,
    pub contact_stiffness_n_per_m: f64,
    pub contact_damping_n_s_per_m: f64,
    pub base_effective_mass_kg: f64,
    pub base_stiffness_n_per_m: f64,
    pub base_damping_n_s_per_m: f64,
    pub gas_rotational_damping_n_m_s: f64,
    pub gas_translation_damping_n_s_per_m: f64,
}

impl CoupledFactors {
    pub fn channel_factors(self) -> CoupledChannelFactors {
        CoupledChannelFactors {
            gravity_m_per_s2: self.gravity_m_per_s2,
            sliding_friction_coefficient: self.sliding_friction_coefficient,
            rolling_resistance_m: self.rolling_resistance_m,
            contact_stiffness_n_per_m: self.contact_stiffness_n_per_m,
            contact_damping_n_s_per_m: self.contact_damping_n_s_per_m,
            base_effective_mass_kg: self.base_effective_mass_kg,
            base_stiffness_n_per_m: self.base_stiffness_n_per_m,
            base_damping_n_s_per_m: self.base_damping_n_s_per_m,
            gas_rotational_damping_n_m_s: self.gas_rotational_damping_n_m_s,
            gas_translation_damping_n_s_per_m: self.gas_translation_damping_n_s_per_m,
        }
    }
}

/// Fixed deterministic run controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoupledControls {
    pub timestep_s: f64,
    pub maximum_steps: u32,
    pub terminal_inclination_rad: f64,
    pub reimpact_limit: u32,
}

/// Initial pose/twist.  Inclination is the angle of the disc plane above the base.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoupledInitialState {
    pub inclination_rad: f64,
    pub precession_rad_per_s: f64,
    pub spin_rad_per_s: f64,
}

/// An independently evaluated channel wrench and its body work over one step.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ChannelWrench {
    pub force_world_n: Vec3,
    pub torque_world_nm: Vec3,
    pub work_j: f64,
}

/// Channel ownership is explicit: no channel may borrow another's work.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ChannelOwnership {
    pub gravity: ChannelWrench,
    pub contact: ChannelWrench,
    pub rolling: ChannelWrench,
    pub base: ChannelWrench,
    pub gas: ChannelWrench,
}

/// One retained, time-evolving sample.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledSample {
    pub time_s: f64,
    pub inclination_rad: f64,
    pub precession_rad_per_s: f64,
    pub spin_rad_per_s: f64,
    /// Finite-difference `d(precession)/dt` [rad/s²].
    pub precession_acceleration_rad_per_s2: f64,
    /// Gap evaluated at the beginning of the interval ending at `time_s` [m].
    pub interval_start_gap_m: f64,
    /// Unilateral normal force evaluated over the interval ending at `time_s` [N].
    pub interval_normal_force_n: f64,
    pub contact_active: bool,
    /// Profile feature selected by the analytic support query at this
    /// sample's retained post-step pose. Cylinder-only compatibility runs do
    /// not expose a feature index.
    pub support_source_feature: Option<usize>,
    pub reimpact_count: u32,
    pub channels: ChannelOwnership,
    pub mechanical_energy_j: f64,
    pub energy_defect_j: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CoupledTerminal {
    TerminalInclination,
    HorizonReached,
    NumericalRefusal,
}

/// Deterministic, restartable reduced state.  It contains only actual state,
/// never an encoded duration or target outcome.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledCheckpoint {
    pub state: RigidBodyState,
    pub time_s: f64,
    pub base_deflection_m: f64,
    pub base_velocity_m_per_s: f64,
    pub accumulated_channel_work_j: [f64; 5],
    pub accumulated_energy_defect_j: f64,
    /// Total body-plus-base-plus-contact energy at the start of this trajectory [J].
    pub initial_total_energy_j: f64,
    /// Stable identity of factors, restart-relevant controls, initial twist, and energy.
    pub configuration_fingerprint: u64,
    pub reimpact_count: u32,
    pub was_in_contact: bool,
}

#[derive(Clone, Debug, PartialEq)]
pub struct CoupledRun {
    pub samples: Vec<CoupledSample>,
    pub checkpoint: CoupledCheckpoint,
    pub terminal: CoupledTerminal,
    pub applicability: &'static str,
    pub model_disagreement: &'static str,
}

#[derive(Clone, Debug, PartialEq)]
pub enum CoupledError {
    InvalidInput(&'static str),
    CheckpointMismatch,
    Dynamics(String),
}

impl fmt::Display for CoupledError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}
impl std::error::Error for CoupledError {}

enum CoupledGeometry<'profile, 'cx> {
    Cylinder(ContactDiscGeometry),
    Profile {
        profile: &'profile ResolvedDiscProfile,
        cx: &'profile Cx<'cx>,
    },
}

impl CoupledGeometry<'_, '_> {
    fn contact(&self, pose: Pose) -> Result<(ContactGeometry, Option<usize>), CoupledError> {
        match self {
            Self::Cylinder(geometry) => contact_geometry(*geometry, pose)
                .map(|contact| (contact, None))
                .map_err(|error| CoupledError::Dynamics(error.to_string())),
            Self::Profile { profile, cx } => {
                profile_contact_geometry(&profile.chart, profile.mass_properties, pose, cx)
                    .map(|geometry| (geometry.contact, Some(geometry.support_source_feature)))
                    .map_err(|error| CoupledError::Dynamics(error.to_string()))
            }
        }
    }

    fn state_at_ground_contact(
        &self,
        orientation: UnitQuaternion,
        linear_momentum_world: Vec3,
        angular_momentum_body: Vec3,
    ) -> Result<RigidBodyState, CoupledError> {
        match self {
            Self::Cylinder(geometry) => state_at_ground_contact(
                *geometry,
                orientation,
                linear_momentum_world,
                angular_momentum_body,
            )
            .map_err(|error| CoupledError::Dynamics(error.to_string())),
            Self::Profile { profile, cx } => profile_state_at_ground_contact(
                &profile.chart,
                profile.mass_properties,
                orientation,
                linear_momentum_world,
                angular_momentum_body,
                cx,
            )
            .map_err(|error| CoupledError::Dynamics(error.to_string())),
        }
    }

    fn identity_word(&self) -> u64 {
        match self {
            Self::Cylinder(_) => 0,
            Self::Profile { profile, .. } => profile.identity.0,
        }
    }
}

/// Starts or resumes a closed reduced trajectory.  Contact is unilateral and
/// may separate/reimpact; rolling/base/gas are recomputed from each new state.
pub fn run_closed_reduced(
    factors: CoupledFactors,
    controls: CoupledControls,
    initial: CoupledInitialState,
    restart: Option<CoupledCheckpoint>,
) -> Result<CoupledRun, CoupledError> {
    let geometry = CoupledGeometry::Cylinder(ContactDiscGeometry {
        radius_m: factors.radius_m,
        thickness_m: factors.thickness_m,
        mass_kg: factors.mass_kg,
    });
    run_closed_with_geometry(factors, controls, initial, restart, geometry)
}

/// Runs the same reduced channel model against a true resolved profile.
///
/// Geometry, support, mass, center of mass, and inertia all come from
/// `profile`. This removes the former cone-inertia/cylinder-contact surrogate;
/// it does not upgrade the reduced point-contact or loss laws to experimental
/// validation.
pub fn run_closed_profile_reduced(
    profile: &ResolvedDiscProfile,
    channels: CoupledChannelFactors,
    controls: CoupledControls,
    initial: CoupledInitialState,
    restart: Option<CoupledCheckpoint>,
    cx: &Cx<'_>,
) -> Result<CoupledRun, CoupledError> {
    let factors = CoupledFactors {
        mass_kg: profile.mass_properties.mass,
        radius_m: profile.dimensions.outer_radius_m,
        thickness_m: profile.dimensions.thickness_m,
        transverse_inertia_kg_m2: profile.mass_properties.principal_inertia.transverse,
        axial_inertia_kg_m2: profile.mass_properties.principal_inertia.axial,
        gravity_m_per_s2: channels.gravity_m_per_s2,
        sliding_friction_coefficient: channels.sliding_friction_coefficient,
        rolling_resistance_m: channels.rolling_resistance_m,
        contact_stiffness_n_per_m: channels.contact_stiffness_n_per_m,
        contact_damping_n_s_per_m: channels.contact_damping_n_s_per_m,
        base_effective_mass_kg: channels.base_effective_mass_kg,
        base_stiffness_n_per_m: channels.base_stiffness_n_per_m,
        base_damping_n_s_per_m: channels.base_damping_n_s_per_m,
        gas_rotational_damping_n_m_s: channels.gas_rotational_damping_n_m_s,
        gas_translation_damping_n_s_per_m: channels.gas_translation_damping_n_s_per_m,
    };
    run_closed_with_geometry(
        factors,
        controls,
        initial,
        restart,
        CoupledGeometry::Profile { profile, cx },
    )
}

fn run_closed_with_geometry(
    factors: CoupledFactors,
    controls: CoupledControls,
    initial: CoupledInitialState,
    restart: Option<CoupledCheckpoint>,
    geometry_model: CoupledGeometry<'_, '_>,
) -> Result<CoupledRun, CoupledError> {
    validate(factors, controls, initial)?;
    let mass = MassProperties::new(
        factors.mass_kg,
        Vec3::ZERO,
        Vec3::new(
            factors.transverse_inertia_kg_m2,
            factors.transverse_inertia_kg_m2,
            factors.axial_inertia_kg_m2,
        ),
    )
    .map_err(|e| CoupledError::Dynamics(e.to_string()))?;
    let gravity = Gravity::new(Vec3::new(0.0, 0.0, -factors.gravity_m_per_s2))
        .map_err(|e| CoupledError::Dynamics(e.to_string()))?;
    let integrator = RigidBodyIntegrator::new(gravity);
    let mut configured_checkpoint = initial_checkpoint(factors, initial, mass, &geometry_model)?;
    let configured_initial_energy_j = total_energy(
        integrator
            .diagnostics(configured_checkpoint.state, mass)
            .map_err(|e| CoupledError::Dynamics(e.to_string()))?
            .mechanical_energy,
        factors,
        configured_checkpoint.base_deflection_m,
        configured_checkpoint.base_velocity_m_per_s,
        initial_penetration_m(&geometry_model, &configured_checkpoint)?,
    );
    let configured_fingerprint = configuration_fingerprint(
        factors,
        controls,
        initial,
        configured_initial_energy_j,
        geometry_model.identity_word(),
    );
    let mut checkpoint = match restart {
        Some(checkpoint)
            if checkpoint.initial_total_energy_j.to_bits()
                == configured_initial_energy_j.to_bits()
                && checkpoint.configuration_fingerprint == configured_fingerprint =>
        {
            checkpoint
        }
        Some(_) => return Err(CoupledError::CheckpointMismatch),
        None => {
            configured_checkpoint.initial_total_energy_j = configured_initial_energy_j;
            configured_checkpoint.configuration_fingerprint = configured_fingerprint;
            configured_checkpoint
        }
    };
    let mut samples = Vec::with_capacity(controls.maximum_steps as usize);
    let mut previous_precession = qois(checkpoint.state, mass)?.1;

    for _ in 0..controls.maximum_steps {
        let state = checkpoint.state;
        let pose = state.pose();
        let normal = pose
            .orientation()
            .rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
        let inclination = normal.z.clamp(-1.0, 1.0).acos();
        if inclination <= controls.terminal_inclination_rad {
            return Ok(CoupledRun {
                samples,
                checkpoint,
                terminal: CoupledTerminal::TerminalInclination,
                applicability: applicability(),
                model_disagreement: disagreement(),
            });
        }
        let velocity = state
            .center_of_mass_velocity_world(mass)
            .map_err(|e| CoupledError::Dynamics(e.to_string()))?;
        let omega_body = mass
            .angular_velocity_body_checked(state.angular_momentum_body())
            .map_err(|e| CoupledError::Dynamics(e.to_string()))?;
        let omega_world = pose.orientation().rotate_body_to_world(omega_body);
        let (geometry, _) = geometry_model.contact(pose)?;
        let contact_arm = geometry.radius_world_m;
        let support_height = geometry.gap_m - checkpoint.base_deflection_m;
        let contact_velocity = velocity.add(omega_world.cross(contact_arm));
        let relative_normal_speed = contact_velocity.z - checkpoint.base_velocity_m_per_s;
        let penetration = (-support_height).max(0.0);
        let normal_force = unilateral_normal_force(
            penetration,
            relative_normal_speed,
            factors.contact_stiffness_n_per_m,
            factors.contact_damping_n_s_per_m,
        );
        let contact_active = normal_force > 0.0;
        if contact_active && !checkpoint.was_in_contact {
            checkpoint.reimpact_count += 1;
        }
        if checkpoint.reimpact_count > controls.reimpact_limit {
            return Ok(CoupledRun {
                samples,
                checkpoint,
                terminal: CoupledTerminal::NumericalRefusal,
                applicability: applicability(),
                model_disagreement: disagreement(),
            });
        }

        let tangent_velocity = Vec3::new(contact_velocity.x, contact_velocity.y, 0.0);
        let tangent_speed = tangent_velocity.norm_squared().sqrt();
        let friction_force = if contact_active && tangent_speed > 1.0e-14 {
            tangent_velocity
                .scale(-factors.sliding_friction_coefficient * normal_force / tangent_speed)
        } else {
            Vec3::ZERO
        };
        let contact_force = Vec3::new(friction_force.x, friction_force.y, normal_force);
        let contact_torque = contact_arm.cross(contact_force);
        let gas_force = velocity.scale(-factors.gas_translation_damping_n_s_per_m);
        let gas_torque = omega_world.scale(-factors.gas_rotational_damping_n_m_s);
        let provisional = integrator
            .step(
                state,
                mass,
                Wrench {
                    force_world: contact_force.add(gas_force),
                    torque_body: pose
                        .orientation()
                        .rotate_world_to_body(contact_torque.add(gas_torque)),
                },
                controls.timestep_s,
            )
            .map_err(|e| CoupledError::Dynamics(e.to_string()))?;
        let provisional_omega_body = mass
            .angular_velocity_body_checked(provisional.state_after.angular_momentum_body())
            .map_err(|e| CoupledError::Dynamics(e.to_string()))?;
        let provisional_omega = provisional
            .state_after
            .pose()
            .orientation()
            .rotate_body_to_world(provisional_omega_body);
        let midpoint_omega_for_rolling = omega_world.add(provisional_omega).scale(0.5);
        let midpoint_state_for_rolling = provisional.state_after;
        let midpoint_normal = midpoint_state_for_rolling
            .pose()
            .orientation()
            .rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
        let rolling_spin_speed = midpoint_omega_for_rolling.dot(midpoint_normal);
        let midpoint_velocity_for_rolling = velocity
            .add(
                midpoint_state_for_rolling
                    .center_of_mass_velocity_world(mass)
                    .map_err(|e| CoupledError::Dynamics(e.to_string()))?,
            )
            .scale(0.5);
        let (midpoint_contact, _) = geometry_model.contact(midpoint_state_for_rolling.pose())?;
        let midpoint_contact_velocity = midpoint_velocity_for_rolling
            .add(midpoint_omega_for_rolling.cross(midpoint_contact.radius_world_m));
        let rolling_contour_speed = Vec3::new(
            midpoint_contact_velocity.x,
            midpoint_contact_velocity.y,
            0.0,
        )
        .norm_squared()
        .sqrt();
        let declared_rolling_power_w = factors.rolling_resistance_m
            * normal_force
            * (rolling_spin_speed.abs() + rolling_contour_speed / factors.radius_m);
        let omega_squared = omega_world.norm_squared();
        let rolling_power_w = if omega_squared > 1.0e-28 {
            declared_rolling_power_w
        } else {
            0.0
        };
        let rolling_torque = if rolling_power_w > 0.0 {
            omega_world.scale(-rolling_power_w / omega_squared)
        } else {
            Vec3::ZERO
        };
        let total_force = contact_force.add(gas_force);
        let total_torque_world = contact_torque.add(rolling_torque).add(gas_torque);
        let total_torque_body = pose.orientation().rotate_world_to_body(total_torque_world);
        let receipt = integrator
            .step(
                state,
                mass,
                Wrench {
                    force_world: total_force,
                    torque_body: total_torque_body,
                },
                controls.timestep_s,
            )
            .map_err(|e| CoupledError::Dynamics(e.to_string()))?;
        let midpoint_velocity = velocity
            .add(
                receipt
                    .state_after
                    .center_of_mass_velocity_world(mass)
                    .map_err(|e| CoupledError::Dynamics(e.to_string()))?,
            )
            .scale(0.5);
        let midpoint_omega = omega_world
            .add(
                receipt
                    .state_after
                    .pose()
                    .orientation()
                    .rotate_body_to_world(
                        mass.angular_velocity_body_checked(
                            receipt.state_after.angular_momentum_body(),
                        )
                        .map_err(|e| CoupledError::Dynamics(e.to_string()))?,
                    ),
            )
            .scale(0.5);
        let work = |force: Vec3, torque: Vec3| {
            (force.dot(midpoint_velocity) + torque.dot(midpoint_omega)) * controls.timestep_s
        };
        let channels = ChannelOwnership {
            gravity: ChannelWrench {
                force_world_n: gravity.acceleration_world().scale(factors.mass_kg),
                torque_world_nm: Vec3::ZERO,
                work_j: work(
                    gravity.acceleration_world().scale(factors.mass_kg),
                    Vec3::ZERO,
                ),
            },
            contact: ChannelWrench {
                force_world_n: contact_force,
                torque_world_nm: contact_torque,
                work_j: work(contact_force, contact_torque),
            },
            rolling: ChannelWrench {
                force_world_n: Vec3::ZERO,
                torque_world_nm: rolling_torque,
                // The reduced generalized loss is non-negative by construction,
                // and this torque satisfies tau dot omega = -P_roll.
                work_j: -rolling_power_w * controls.timestep_s,
            },
            base: ChannelWrench {
                force_world_n: Vec3::ZERO,
                torque_world_nm: Vec3::ZERO,
                work_j: -factors.base_damping_n_s_per_m
                    * checkpoint.base_velocity_m_per_s.powi(2)
                    * controls.timestep_s,
            },
            gas: ChannelWrench {
                force_world_n: gas_force,
                torque_world_nm: gas_torque,
                work_j: work(gas_force, gas_torque),
            },
        };
        let new_base_acceleration = (-normal_force
            - factors.base_stiffness_n_per_m * checkpoint.base_deflection_m
            - factors.base_damping_n_s_per_m * checkpoint.base_velocity_m_per_s)
            / factors.base_effective_mass_kg;
        checkpoint.base_velocity_m_per_s += new_base_acceleration * controls.timestep_s;
        checkpoint.base_deflection_m += checkpoint.base_velocity_m_per_s * controls.timestep_s;
        checkpoint.state = receipt.state_after;
        checkpoint.time_s += controls.timestep_s;
        checkpoint.was_in_contact = contact_active;
        let tangential_work = work(friction_force, contact_arm.cross(friction_force));
        let contact_damping_work = if penetration > 0.0 {
            -factors.contact_damping_n_s_per_m
                * relative_normal_speed.min(0.0).powi(2)
                * controls.timestep_s
        } else {
            0.0
        };
        checkpoint.accumulated_channel_work_j[0] += channels.gravity.work_j;
        checkpoint.accumulated_channel_work_j[1] += tangential_work + contact_damping_work;
        checkpoint.accumulated_channel_work_j[2] += channels.rolling.work_j;
        checkpoint.accumulated_channel_work_j[3] += channels.base.work_j;
        checkpoint.accumulated_channel_work_j[4] += channels.gas.work_j;
        let (post_geometry, post_support_source_feature) =
            geometry_model.contact(checkpoint.state.pose())?;
        let post_penetration = (checkpoint.base_deflection_m - post_geometry.gap_m).max(0.0);
        let total_energy = total_energy(
            receipt.diagnostics_after.mechanical_energy,
            factors,
            checkpoint.base_deflection_m,
            checkpoint.base_velocity_m_per_s,
            post_penetration,
        );
        let defect = (total_energy - checkpoint.initial_total_energy_j)
            - (checkpoint.accumulated_channel_work_j[1]
                + checkpoint.accumulated_channel_work_j[2]
                + checkpoint.accumulated_channel_work_j[3]
                + checkpoint.accumulated_channel_work_j[4]);
        checkpoint.accumulated_energy_defect_j = defect;
        let (sample_inclination, precession, spin) = qois(checkpoint.state, mass)?;
        let precession_acceleration = (precession - previous_precession) / controls.timestep_s;
        previous_precession = precession;
        samples.push(CoupledSample {
            time_s: checkpoint.time_s,
            inclination_rad: sample_inclination,
            precession_rad_per_s: precession,
            spin_rad_per_s: spin,
            precession_acceleration_rad_per_s2: precession_acceleration,
            interval_start_gap_m: support_height,
            interval_normal_force_n: normal_force,
            contact_active,
            support_source_feature: post_support_source_feature,
            reimpact_count: checkpoint.reimpact_count,
            channels,
            mechanical_energy_j: total_energy,
            energy_defect_j: defect,
        });
        if !defect.is_finite() || !checkpoint.base_deflection_m.is_finite() {
            return Ok(CoupledRun {
                samples,
                checkpoint,
                terminal: CoupledTerminal::NumericalRefusal,
                applicability: applicability(),
                model_disagreement: disagreement(),
            });
        }
        // Classify an event reached on the final allowed step as physical,
        // rather than incorrectly turning it into a horizon censor. The event
        // time is still the first committed fixed-step boundary at or below
        // the threshold; this reduced runner does not claim substep event
        // interpolation.
        if sample_inclination <= controls.terminal_inclination_rad {
            return Ok(CoupledRun {
                samples,
                checkpoint,
                terminal: CoupledTerminal::TerminalInclination,
                applicability: applicability(),
                model_disagreement: disagreement(),
            });
        }
    }
    Ok(CoupledRun {
        samples,
        checkpoint,
        terminal: CoupledTerminal::HorizonReached,
        applicability: applicability(),
        model_disagreement: disagreement(),
    })
}

fn initial_checkpoint(
    factors: CoupledFactors,
    initial: CoupledInitialState,
    mass: MassProperties,
    geometry: &CoupledGeometry<'_, '_>,
) -> Result<CoupledCheckpoint, CoupledError> {
    let orientation =
        UnitQuaternion::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), initial.inclination_rad)
            .map_err(|e| CoupledError::Dynamics(e.to_string()))?;
    let normal = orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
    let omega_world =
        Vec3::new(0.0, 0.0, initial.precession_rad_per_s).add(normal.scale(initial.spin_rad_per_s));
    let omega_body = orientation.rotate_world_to_body(omega_world);
    let angular = Vec3::new(
        omega_body.x * mass.principal_inertia_body().x,
        omega_body.y * mass.principal_inertia_body().y,
        omega_body.z * mass.principal_inertia_body().z,
    );
    let ground_state = geometry.state_at_ground_contact(orientation, Vec3::ZERO, angular)?;
    let contact_arm = geometry.contact(ground_state.pose())?.0.radius_world_m;
    // This is the declared rolling initial condition: the current analytic
    // profile support point has zero world velocity before the preload offset.
    let center_of_mass_velocity = omega_world.cross(contact_arm).scale(-1.0);
    let preload_force_n = factors.mass_kg * factors.gravity_m_per_s2;
    let (base_deflection_m, contact_penetration_m) =
        if factors.contact_stiffness_n_per_m > 0.0 && factors.base_stiffness_n_per_m > 0.0 {
            (
                -preload_force_n / factors.base_stiffness_n_per_m,
                preload_force_n / factors.contact_stiffness_n_per_m,
            )
        } else {
            (0.0, 0.0)
        };
    let grounded_state = geometry.state_at_ground_contact(
        orientation,
        center_of_mass_velocity.scale(factors.mass_kg),
        angular,
    )?;
    let shifted_position = grounded_state.pose().position_world().add(Vec3::new(
        0.0,
        0.0,
        base_deflection_m - contact_penetration_m,
    ));
    let shifted_pose = Pose::new(shifted_position, orientation)
        .map_err(|error| CoupledError::Dynamics(error.to_string()))?;
    Ok(CoupledCheckpoint {
        state: RigidBodyState::new(
            shifted_pose,
            center_of_mass_velocity.scale(factors.mass_kg),
            angular,
        )
        .map_err(|error| CoupledError::Dynamics(error.to_string()))?,
        time_s: 0.0,
        base_deflection_m,
        base_velocity_m_per_s: 0.0,
        accumulated_channel_work_j: [0.0; 5],
        accumulated_energy_defect_j: 0.0,
        initial_total_energy_j: 0.0,
        configuration_fingerprint: 0,
        reimpact_count: 0,
        was_in_contact: contact_penetration_m > 0.0,
    })
}

/// Decomposes the symmetry-axis motion without treating an arbitrary horizontal
/// angular-velocity projection as precession.
pub fn qois(state: RigidBodyState, mass: MassProperties) -> Result<(f64, f64, f64), CoupledError> {
    let normal = state
        .pose()
        .orientation()
        .rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
    let inclination = normal.z.clamp(-1.0, 1.0).acos();
    let omega_body = mass
        .angular_velocity_body_checked(state.angular_momentum_body())
        .map_err(|e| CoupledError::Dynamics(e.to_string()))?;
    let omega = state.pose().orientation().rotate_body_to_world(omega_body);
    let normal_dot = omega.cross(normal);
    let horizontal_norm_squared = normal.x.mul_add(normal.x, normal.y * normal.y);
    if horizontal_norm_squared <= 1.0e-14 {
        return Err(CoupledError::InvalidInput(
            "precession undefined at symmetry-axis pole",
        ));
    }
    let precession = (normal.x * normal_dot.y - normal.y * normal_dot.x) / horizontal_norm_squared;
    // Intrinsic spin is the residual axial rate after removing the declared
    // world-vertical precession component.
    Ok((
        inclination,
        precession,
        omega.dot(normal) - precession * normal.z,
    ))
}

/// Returns the QoIs of the declared initial twist before any force is applied.
pub fn initial_qois(
    factors: CoupledFactors,
    initial: CoupledInitialState,
) -> Result<(f64, f64, f64), CoupledError> {
    validate(
        factors,
        CoupledControls {
            timestep_s: 1.0,
            maximum_steps: 1,
            terminal_inclination_rad: 0.001,
            reimpact_limit: 0,
        },
        initial,
    )?;
    let mass = MassProperties::new(
        factors.mass_kg,
        Vec3::ZERO,
        Vec3::new(
            factors.transverse_inertia_kg_m2,
            factors.transverse_inertia_kg_m2,
            factors.axial_inertia_kg_m2,
        ),
    )
    .map_err(|e| CoupledError::Dynamics(e.to_string()))?;
    let geometry = ContactDiscGeometry {
        radius_m: factors.radius_m,
        thickness_m: factors.thickness_m,
        mass_kg: factors.mass_kg,
    };
    qois(
        initial_checkpoint(factors, initial, mass, &CoupledGeometry::Cylinder(geometry))?.state,
        mass,
    )
}

/// Returns the finite-cylinder rim velocity at the initialized contact point.
/// The preload changes only position, so this verifies the rolling twist before
/// the preload offset as well.
pub fn initial_contact_point_velocity(
    factors: CoupledFactors,
    initial: CoupledInitialState,
) -> Result<Vec3, CoupledError> {
    validate(
        factors,
        CoupledControls {
            timestep_s: 1.0,
            maximum_steps: 1,
            terminal_inclination_rad: 0.001,
            reimpact_limit: 0,
        },
        initial,
    )?;
    let mass = MassProperties::new(
        factors.mass_kg,
        Vec3::ZERO,
        Vec3::new(
            factors.transverse_inertia_kg_m2,
            factors.transverse_inertia_kg_m2,
            factors.axial_inertia_kg_m2,
        ),
    )
    .map_err(|e| CoupledError::Dynamics(e.to_string()))?;
    let geometry = ContactDiscGeometry {
        radius_m: factors.radius_m,
        thickness_m: factors.thickness_m,
        mass_kg: factors.mass_kg,
    };
    let geometry = CoupledGeometry::Cylinder(geometry);
    let checkpoint = initial_checkpoint(factors, initial, mass, &geometry)?;
    let contact = geometry.contact(checkpoint.state.pose())?.0;
    let velocity = checkpoint
        .state
        .center_of_mass_velocity_world(mass)
        .map_err(|e| CoupledError::Dynamics(e.to_string()))?;
    let omega_body = mass
        .angular_velocity_body_checked(checkpoint.state.angular_momentum_body())
        .map_err(|e| CoupledError::Dynamics(e.to_string()))?;
    let omega_world = checkpoint
        .state
        .pose()
        .orientation()
        .rotate_body_to_world(omega_body);
    Ok(velocity.add(omega_world.cross(contact.radius_world_m)))
}

fn initial_penetration_m(
    geometry: &CoupledGeometry<'_, '_>,
    checkpoint: &CoupledCheckpoint,
) -> Result<f64, CoupledError> {
    let contact = geometry.contact(checkpoint.state.pose())?.0;
    Ok((checkpoint.base_deflection_m - contact.gap_m).max(0.0))
}

fn total_energy(
    body_mechanical_energy_j: f64,
    factors: CoupledFactors,
    base_deflection_m: f64,
    base_velocity_m_per_s: f64,
    penetration_m: f64,
) -> f64 {
    body_mechanical_energy_j
        + 0.5 * factors.base_effective_mass_kg * base_velocity_m_per_s.powi(2)
        + 0.5 * factors.base_stiffness_n_per_m * base_deflection_m.powi(2)
        + 0.5 * factors.contact_stiffness_n_per_m * penetration_m.powi(2)
}

fn unilateral_normal_force(
    penetration_m: f64,
    relative_normal_speed_m_per_s: f64,
    stiffness_n_per_m: f64,
    damping_n_s_per_m: f64,
) -> f64 {
    // Kelvin-Voigt damping belongs to the closed contact branch. Applying the
    // dashpot while the gap is open would create a non-physical attractive
    // pre-contact force and could turn harmless approach/separation into
    // artificial reimpact chatter.
    if penetration_m <= 0.0 {
        return 0.0;
    }
    (stiffness_n_per_m * penetration_m
        + damping_n_s_per_m * (-relative_normal_speed_m_per_s).max(0.0))
    .max(0.0)
}

#[cfg(test)]
mod tests {
    use super::unilateral_normal_force;

    #[test]
    fn kelvin_voigt_dashpot_is_inactive_across_an_open_gap() {
        assert_eq!(unilateral_normal_force(0.0, -3.0, 80_000.0, 3.0), 0.0);
        assert_eq!(unilateral_normal_force(-1.0e-6, -3.0, 80_000.0, 3.0), 0.0);
        assert!(unilateral_normal_force(1.0e-6, -3.0, 80_000.0, 3.0) > 0.0);
    }
}

/// Public integration plan for replacing each reduced channel law with the
/// corresponding exported receipt adapter after focused compilation proves the
/// adapter boundary: `mechanics` supplies contact/base, `rolling_contact`
/// supplies rolling, and `air` supplies gas.  Until then, this runner makes no
/// adapter-consumption claim.
pub const ADAPTER_INTEGRATION_PLAN: &str = "after focused adapter compilation: mechanics=>contact/base;rolling_contact=>rolling;air=>gas;retain independent work ownership and energy closure";

fn applicability() -> &'static str {
    "reduced-rigid-Kelvin-Voigt-contact-and-one-mode-base;profile path uses one analytic axisymmetric chart for support and mass properties;gross-sliding-Coulomb-without-static-stick-resolution;not-finite-patch-or-resolved-gas-film"
}
fn disagreement() -> &'static str {
    "profile geometry is physical line/arc support but contact remains a reduced point Kelvin-Voigt law;rolling/base/gas are reduced channel laws and the runner does not yet consume the higher-fidelity transactional adapters"
}

fn configuration_fingerprint(
    factors: CoupledFactors,
    controls: CoupledControls,
    initial: CoupledInitialState,
    initial_energy_j: f64,
    geometry_identity: u64,
) -> u64 {
    // A fixed FNV-1a mix of IEEE-754 encodings is portable and deterministic;
    // this is an identity binding, not a cryptographic digest.
    let mut fingerprint = 0xcbf2_9ce4_8422_2325_u64;
    for bits in [
        factors.mass_kg.to_bits(),
        factors.radius_m.to_bits(),
        factors.thickness_m.to_bits(),
        factors.transverse_inertia_kg_m2.to_bits(),
        factors.axial_inertia_kg_m2.to_bits(),
        factors.gravity_m_per_s2.to_bits(),
        factors.sliding_friction_coefficient.to_bits(),
        factors.rolling_resistance_m.to_bits(),
        factors.contact_stiffness_n_per_m.to_bits(),
        factors.contact_damping_n_s_per_m.to_bits(),
        factors.base_effective_mass_kg.to_bits(),
        factors.base_stiffness_n_per_m.to_bits(),
        factors.base_damping_n_s_per_m.to_bits(),
        factors.gas_rotational_damping_n_m_s.to_bits(),
        factors.gas_translation_damping_n_s_per_m.to_bits(),
        controls.timestep_s.to_bits(),
        controls.terminal_inclination_rad.to_bits(),
        u64::from(controls.reimpact_limit),
        initial.inclination_rad.to_bits(),
        initial.precession_rad_per_s.to_bits(),
        initial.spin_rad_per_s.to_bits(),
        initial_energy_j.to_bits(),
        geometry_identity,
    ] {
        fingerprint ^= bits;
        fingerprint = fingerprint.wrapping_mul(0x0000_0100_0000_01b3);
    }
    fingerprint
}
fn validate(
    f: CoupledFactors,
    c: CoupledControls,
    i: CoupledInitialState,
) -> Result<(), CoupledError> {
    for x in [
        f.mass_kg,
        f.radius_m,
        f.thickness_m,
        f.transverse_inertia_kg_m2,
        f.axial_inertia_kg_m2,
        f.gravity_m_per_s2,
        f.sliding_friction_coefficient,
        f.rolling_resistance_m,
        f.contact_stiffness_n_per_m,
        f.contact_damping_n_s_per_m,
        f.base_effective_mass_kg,
        f.base_stiffness_n_per_m,
        f.base_damping_n_s_per_m,
        f.gas_rotational_damping_n_m_s,
        f.gas_translation_damping_n_s_per_m,
        c.timestep_s,
        c.terminal_inclination_rad,
        i.inclination_rad,
        i.precession_rad_per_s,
        i.spin_rad_per_s,
    ] {
        if !x.is_finite() {
            return Err(CoupledError::InvalidInput("non-finite factor"));
        }
    }
    if f.mass_kg <= 0.0
        || f.radius_m <= 0.0
        || f.thickness_m <= 0.0
        || f.transverse_inertia_kg_m2 <= 0.0
        || f.axial_inertia_kg_m2 <= 0.0
        || f.gravity_m_per_s2 <= 0.0
        || f.sliding_friction_coefficient < 0.0
        || f.rolling_resistance_m < 0.0
        || f.contact_stiffness_n_per_m < 0.0
        || f.contact_damping_n_s_per_m < 0.0
        || f.base_effective_mass_kg <= 0.0
        || f.base_stiffness_n_per_m < 0.0
        || f.base_damping_n_s_per_m < 0.0
        || f.gas_rotational_damping_n_m_s < 0.0
        || f.gas_translation_damping_n_s_per_m < 0.0
        || c.timestep_s <= 0.0
        || c.maximum_steps == 0
        || c.terminal_inclination_rad <= 0.0
        || i.inclination_rad <= c.terminal_inclination_rad
    {
        return Err(CoupledError::InvalidInput("out-of-domain factor"));
    }
    Ok(())
}
