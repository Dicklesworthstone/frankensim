//! Material-bound anisotropic Stokes resistance on the native rigid-body owner.
//!
//! For an unbounded, no-slip triaxial ellipsoid, let
//! Delta(t)=sqrt((a1²+t)(a2²+t)(a3²+t)), chi=int_0^infinity dt/Delta,
//! alpha_i=int_0^infinity dt/((ai²+t)Delta). The resistance eigenvalues are
//! C_i=16*pi*mu/(chi+ai²*alpha_i) and
//! D_i=(16*pi*mu/3)*(aj²+ak²)/(aj²*alpha_j+ak²*alpha_k).
//! Translation/rotation cross-coupling vanishes at the ellipsoid centre.
//! See Anand & Narsimhan, JFM 983 (2024), Appendix B, equations B1--B4:
//! <https://doi.org/10.1017/jfm.2024.136>. Their integral definitions include
//! a1*a2*a3; that common factor is cancelled in the formulas above.
//!
//! Geometry and the original fluid determine all six resistances. A matching
//! uniform material region supplies mass and full inertia. Authored ellipsoid
//! axes are explicitly transformed into the native principal frame; their order
//! is never assumed to match the inertia eigensolver's sorted axes.
//!
//! The existing rigid-body integrator owns momentum, gyroscopic rotation,
//! position and attitude. A bounded load iteration closes rotary resistance
//! against its actual midpoint angular velocity. Translation is eliminated
//! algebraically in the ellipsoid axes at the native midpoint attitude. A
//! prescribed rigid ambient flow includes the local COM velocity and half-curl;
//! its spatial variation participates in the native midpoint position solve. All
//! work is paired with actual endpoint-average velocities, with an independent
//! storage/heat/drive/maintained-flow check before returning a candidate.
//!
//! Scope: frozen uniform Newtonian material in a steady translating/rotating
//! ambient flow, uniform rigid solid, full immersion,
//! no slip, low Reynolds number and fast fluid diffusion relative to body
//! relaxation. The longest semiaxis sets conservative Reynolds/diffusion scales.
//! Finite-Re corrections, fluid memory/added mass, ambient strain,
//! flow, walls, free surfaces, rarefaction and interaction between particles
//! need different models. Do not add an overlapping exterior or gap drag law.

use crate::creeping_flow::{
    CreepingFlowError, CreepingFlowLimits, CreepingLoad, CreepingStep, CreepingStepControls,
};
use crate::principal::PrincipalBody;
use crate::{MaterialInertiaError, VolumeMoments};
use fs_ga::{Mat3, Se3, Vec3 as GeometryVector};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_math::det;
use fs_mbd::{
    DynamicsError, Gravity, MassProperties, RigidBodyIntegrator, RigidBodyState, UnitQuaternion,
    Vec3, Wrench,
};

mod ambient;
pub use ambient::RigidAmbientFlow;
mod resistance;
pub use resistance::{EllipsoidQuadrature, EllipsoidResistance};

/// Source/shape, regime, numerical work or native mechanical admission failure.
#[derive(Clone, Debug, PartialEq)]
pub enum EllipsoidFlowError {
    Invalid(&'static str),
    Material(MaterialError),
    Geometry(MaterialInertiaError),
    Dynamics(DynamicsError),
    /// The shared creeping-flow Reynolds, diffusion or energy condition failed.
    Regime(CreepingFlowError),
    /// Actual midpoint force/torque failed the declared closure tolerance.
    Coupling {
        relative_residual: f64,
        tolerance: f64,
        iterations: usize,
    },
    Cancelled,
}
impl core::fmt::Display for EllipsoidFlowError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material creeping ellipsoid: {self:?}")
    }
}
impl std::error::Error for EllipsoidFlowError {}
impl From<MaterialError> for EllipsoidFlowError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<MaterialInertiaError> for EllipsoidFlowError {
    fn from(error: MaterialInertiaError) -> Self {
        Self::Geometry(error)
    }
}
impl From<DynamicsError> for EllipsoidFlowError {
    fn from(error: DynamicsError) -> Self {
        Self::Dynamics(error)
    }
}
type Error = EllipsoidFlowError;

/// Explicit actual exterior, independent of material or an inertia-fit inference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EllipsoidBoundary {
    /// Positive physical semiaxes [m], in the declared ellipsoid frame.
    pub semiaxes_m: GeometryVector,
    /// Places the centred ellipsoid in original material body's authored coordinates.
    pub body_from_ellipsoid: Se3,
}
impl EllipsoidBoundary {
    /// Feed the same analytic exterior to the material inertia owner. Matching
    /// moments is necessary, not proof of an arbitrary body's actual boundary.
    pub fn volume_moments(self) -> Result<VolumeMoments, MaterialInertiaError> {
        VolumeMoments::ellipsoid(self.semiaxes_m)?.transformed(self.body_from_ellipsoid)
    }
}

/// Numerical controls for the shared native midpoint step and load closure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EllipsoidStepControls {
    pub mechanical: CreepingStepControls,
    /// Bounded native rotational trial evaluations, in 1..=64.
    pub maximum_coupling_iterations: usize,
    /// Scaled force/torque closure residual in [64*epsilon, 1e-8]. Work is gated separately.
    pub relative_wrench_tolerance: f64,
    /// Maximum dt*|omega_mid| and dt*|Omega_ambient|, in (0,0.5] rad.
    /// This is not a phase-error certificate.
    pub maximum_rotation_increment_rad: f64,
}

/// Original native motion/work receipt and actual load-closure diagnostics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EllipsoidStep {
    /// Same load/work fields as creeping sphere. Reynolds uses longest semiaxis.
    pub flow: CreepingStep,
    pub rotational_iterations: usize,
    pub relative_wrench_residual: f64,
    pub rotation_increment_rad: f64,
    pub ambient_rotation_increment_rad: f64,
}

/// An immutable exterior-fluid source and one declared ellipsoidal boundary.
#[derive(Debug)]
pub struct CreepingEllipsoid<'a, F: FluidTransport + ?Sized> {
    source: &'a F,
    transport: FluidTransportState,
    boundary: EllipsoidBoundary,
    geometry: VolumeMoments,
    resistance: EllipsoidResistance,
    longest_semiaxis_m: f64,
    fluid_mass_kg: f64,
    diffusion_time_s: f64,
    gravity: Gravity,
    ambient: RigidAmbientFlow,
    limits: CreepingFlowLimits,
}
struct BodyFrame {
    properties: MassProperties,
    ellipsoid_from_principal: Mat3,
    principal_from_ellipsoid: Mat3,
    rotary_principal: Mat3,
}

impl<'a, F: FluidTransport + ?Sized> CreepingEllipsoid<'a, F> {
    /// Resolve source and six resistance coefficients once; later steps use only
    /// fixed-size tensors and native trial states, without quadrature or allocation.
    pub fn new(
        source: &'a F,
        boundary: EllipsoidBoundary,
        gravity: Gravity,
        ambient_velocity_world: Vec3,
        limits: CreepingFlowLimits,
        quadrature: EllipsoidQuadrature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        let geometry = boundary.volume_moments()?;
        if !ambient_velocity_world.is_finite()
            || !small_positive(limits.maximum_reynolds, 0.1)
            || !small_positive(limits.maximum_diffusion_ratio, 0.1)
        {
            return Err(Error::Invalid(
                "finite ambient flow and creeping limits in (0,0.1] required",
            ));
        }
        let transport = source.transport_state()?;
        let axes = [
            boundary.semiaxes_m.x,
            boundary.semiaxes_m.y,
            boundary.semiaxes_m.z,
        ];
        let resistance = resistance::prepare(
            axes,
            transport.dynamic_viscosity_pa_s(),
            quadrature,
            &mut cancelled,
        )?;
        let longest_semiaxis_m = positive(axes[0].max(axes[1]).max(axes[2]))?;
        let fluid_mass_kg = positive(transport.density_kg_m3() * geometry.volume_m3())?;
        let diffusion_time_s = positive(
            longest_semiaxis_m * (longest_semiaxis_m / transport.kinematic_viscosity_m2_s()),
        )?;
        poll(&mut cancelled)?;
        Ok(Self {
            source,
            transport,
            boundary,
            geometry,
            resistance,
            longest_semiaxis_m,
            fluid_mass_kg,
            diffusion_time_s,
            gravity,
            ambient: RigidAmbientFlow::new(Vec3::ZERO, ambient_velocity_world, Vec3::ZERO)?,
            limits,
        })
    }
    pub const fn source(&self) -> &'a F {
        self.source
    }
    pub const fn transport(&self) -> FluidTransportState {
        self.transport
    }
    pub const fn boundary(&self) -> EllipsoidBoundary {
        self.boundary
    }
    pub const fn resistance(&self) -> &EllipsoidResistance {
        &self.resistance
    }

    /// Select a steady rigid ambient flow while retaining the original source,
    /// geometry and prepared coefficients. The field must remain steady over
    /// the trajectory; abrupt source/flow changes need a fluid-memory model.
    /// Ambient strain and a finite rotating vessel are not represented.
    pub fn with_rigid_ambient(mut self, ambient: RigidAmbientFlow) -> Result<Self, Error> {
        let rotation = finite(ambient.angular_speed_rad_s() * self.diffusion_time_s)?;
        if rotation > self.limits.maximum_reynolds {
            return Err(Error::Regime(CreepingFlowError::Reynolds {
                translation: 0.0,
                rotation,
                maximum: self.limits.maximum_reynolds,
            }));
        }
        if rotation > self.limits.maximum_diffusion_ratio {
            return Err(Error::Regime(CreepingFlowError::Unsteady {
                diffusion_ratio: rotation,
                maximum: self.limits.maximum_diffusion_ratio,
            }));
        }
        self.ambient = ambient;
        Ok(self)
    }
    pub const fn ambient_flow(&self) -> RigidAmbientFlow {
        self.ambient
    }

    fn body_frame(&self, body: &PrincipalBody<'_>) -> Result<BodyFrame, Error> {
        let regions = body.source().regions();
        if regions.len() != 1 {
            return Err(Error::Invalid(
                "ellipsoid flow requires one uniform matching material region",
            ));
        }
        let actual = regions[0].geometry();
        let tensor_scale = self.geometry.inertia_m5().max_abs();
        let centroid_error = actual.centroid_m() - self.geometry.centroid_m();
        if (actual.volume_m3() / self.geometry.volume_m3() - 1.0).abs() > 1e-10
            || centroid_error
                .x
                .abs()
                .max(centroid_error.y.abs())
                .max(centroid_error.z.abs())
                > 1e-10 * self.longest_semiaxis_m
            || actual
                .inertia_m5()
                .m
                .iter()
                .zip(self.geometry.inertia_m5().m)
                .any(|(a, b)| (*a - b).abs() > 1e-10 * tensor_scale)
        {
            return Err(Error::Invalid(
                "declared ellipsoid and original material volume/centroid/tensor disagree",
            ));
        }
        let ellipsoid_from_principal = self
            .boundary
            .body_from_ellipsoid
            .rotation()
            .matrix()
            .transpose()
            .compose(body.body_from_principal().rotation().matrix());
        let principal_from_ellipsoid = ellipsoid_from_principal.transpose();
        let d = self.resistance.rotation_kg_m2_s;
        let rotary_principal = principal_from_ellipsoid
            .compose(Mat3 {
                m: [d[0], 0.0, 0.0, 0.0, d[1], 0.0, 0.0, 0.0, d[2]],
            })
            .compose(ellipsoid_from_principal);
        Ok(BodyFrame {
            properties: body.properties(),
            ellipsoid_from_principal,
            principal_from_ellipsoid,
            rotary_principal,
        })
    }
    /// Actual force/torque and powers; body gravity remains in native dynamics.
    /// Maintained-flow power is `F_drag·u(COM) + torque·Omega` in matching frames.
    /// The reported rotational Reynolds number uses the sum of relative-spin
    /// and ambient angular speeds, scaled by longest-semiaxis squared / nu.
    pub fn load(
        &self,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CreepingLoad, Error> {
        poll(&mut cancelled)?;
        let frame = self.body_frame(body)?;
        let load = self.state_load(&frame, state)?;
        poll(&mut cancelled)?;
        Ok(load)
    }
    fn state_load(&self, frame: &BodyFrame, state: RigidBodyState) -> Result<CreepingLoad, Error> {
        self.load_at(
            frame,
            state.pose().orientation(),
            state.pose().position_world(),
            state.center_of_mass_velocity_world(frame.properties)?,
            frame
                .properties
                .angular_velocity_body_checked(state.angular_momentum_body())?,
        )
    }
    fn load_at(
        &self,
        frame: &BodyFrame,
        attitude: UnitQuaternion,
        position: Vec3,
        velocity: Vec3,
        omega: Vec3,
    ) -> Result<CreepingLoad, Error> {
        let ambient_velocity = self.ambient.velocity_at(position)?;
        let relative = vector(velocity.sub(ambient_velocity))?;
        let ambient_omega = self.ambient.angular_velocity_world_rad_s();
        let (ambient_omega_body, relative_omega) = if ambient_omega == Vec3::ZERO {
            (Vec3::ZERO, omega)
        } else {
            let background = vector(attitude.rotate_world_to_body(ambient_omega))?;
            (background, vector(omega.sub(background))?)
        };
        let translation_reynolds = self
            .transport
            .reynolds(norm(relative)?, 2.0 * self.longest_semiaxis_m)?;
        // Co-rotation must not bypass the ambient-flow inertia restriction.
        let rotation_reynolds = finite(
            (norm(relative_omega)? + self.ambient.angular_speed_rad_s()) * self.diffusion_time_s,
        )?;
        if translation_reynolds > self.limits.maximum_reynolds
            || rotation_reynolds > self.limits.maximum_reynolds
        {
            return Err(Error::Regime(CreepingFlowError::Reynolds {
                translation: translation_reynolds,
                rotation: rotation_reynolds,
                maximum: self.limits.maximum_reynolds,
            }));
        }
        let relative_axes = apply(
            frame.ellipsoid_from_principal,
            attitude.rotate_world_to_body(relative),
        )?;
        let omega_axes = apply(frame.ellipsoid_from_principal, relative_omega)?;
        let force_axes = diagonal(self.resistance.translation_kg_s, relative_axes)?.scale(-1.0);
        let torque_axes = diagonal(self.resistance.rotation_kg_m2_s, omega_axes)?.scale(-1.0);
        let drag = vector(
            attitude.rotate_body_to_world(apply(frame.principal_from_ellipsoid, force_axes)?),
        )?;
        let rotary = apply(frame.principal_from_ellipsoid, torque_axes)?;
        let buoyancy = vector(self.gravity.acceleration_world().scale(-self.fluid_mass_kg))?;
        let dissipated_power_w =
            finite(-force_axes.dot(relative_axes) - torque_axes.dot(omega_axes))?;
        if dissipated_power_w < 0.0
            || (relative != Vec3::ZERO || relative_omega != Vec3::ZERO) && dissipated_power_w == 0.0
        {
            return Err(Error::Invalid(
                "nonzero ellipsoid slip/spin requires representable positive viscous power",
            ));
        }
        Ok(CreepingLoad {
            wrench: Wrench {
                force_world: vector(drag.add(buoyancy))?,
                torque_body: rotary,
            },
            drag_force_world_n: drag,
            rotary_torque_body_nm: rotary,
            buoyancy_force_world_n: buoyancy,
            dissipated_power_w,
            background_power_w: if ambient_omega == Vec3::ZERO {
                finite(drag.dot(ambient_velocity))?
            } else {
                finite(drag.dot(ambient_velocity) + rotary.dot(ambient_omega_body))?
            },
            translation_reynolds,
            rotation_reynolds,
        })
    }
    /// Native motion without an additional applied wrench.
    pub fn step(
        &self,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        controls: EllipsoidStepControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<EllipsoidStep, Error> {
        self.step_with_wrench(body, state, Wrench::ZERO, controls, cancelled)
    }
    /// Compose held world-frame COM force and principal-body free couple with
    /// anisotropic drag. Weight, buoyancy and drag must not be duplicated.
    /// Offset-point forces/world-held couples need their own stage conversion.
    /// Every trial starts from the same accepted state. Late failure or
    /// cancellation never advances accepted body state.
    pub fn step_with_wrench(
        &self,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        applied: Wrench,
        controls: EllipsoidStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EllipsoidStep, Error> {
        poll(&mut cancelled)?;
        validate_controls(controls)?;
        vector(applied.force_world)?;
        vector(applied.torque_body)?;
        let frame = self.body_frame(body)?;
        let p = frame.properties;
        let (shortest_time, diffusion_ratio) = self.timescales(&frame)?;
        let dt = controls.mechanical.duration_s;
        let ambient_omega = self.ambient.angular_velocity_world_rad_s();
        let ambient_rotation_increment_rad = finite(dt * self.ambient.angular_speed_rad_s())?;
        if ambient_rotation_increment_rad > controls.maximum_rotation_increment_rad {
            return Err(Error::Invalid(
                "ambient flow turns beyond the declared midpoint increment; reduce dt",
            ));
        }
        if finite(dt / shortest_time)? > controls.mechanical.maximum_relaxation_fraction {
            return Err(Error::Invalid(
                "ellipsoid step exceeds its body-relaxation fraction; reduce dt",
            ));
        }
        let load_before = self.state_load(&frame, state)?;
        let initial_velocity = state.center_of_mass_velocity_world(p)?;
        let initial_omega = p.angular_velocity_body_checked(state.angular_momentum_body())?;
        let buoyant_drive = vector(load_before.buoyancy_force_world_n.add(applied.force_world))?;
        let integrator = RigidBodyIntegrator::new(self.gravity);
        let mut torque = load_before.rotary_torque_body_nm;
        let mut angular_result = None;
        let mut residual = f64::INFINITY;
        for iteration in 1..=controls.maximum_coupling_iterations {
            poll(&mut cancelled)?;
            let trial = integrator.step(
                state,
                p,
                Wrench {
                    force_world: buoyant_drive,
                    torque_body: vector(torque.add(applied.torque_body))?,
                },
                dt,
            )?;
            let angular_mid = p.angular_velocity_body_checked(
                state
                    .angular_momentum_body()
                    .scale(0.5)
                    .add(trial.state_after.angular_momentum_body().scale(0.5)),
            )?;
            let increment = finite(dt * norm(angular_mid)?)?;
            if increment > controls.maximum_rotation_increment_rad {
                return Err(Error::Invalid(
                    "ellipsoid midpoint rotation exceeds its declared increment; reduce dt",
                ));
            }
            let relative_angular_mid = if ambient_omega == Vec3::ZERO {
                angular_mid
            } else {
                let trial_attitude_mid = midpoint_attitude(
                    state.pose().orientation(),
                    trial.state_after.pose().orientation(),
                )?;
                vector(angular_mid.sub(trial_attitude_mid.rotate_world_to_body(ambient_omega)))?
            };
            let desired = apply(frame.rotary_principal, relative_angular_mid)?.scale(-1.0);
            residual = relative_residual(
                desired,
                torque,
                norm(applied.torque_body)?.max(norm(load_before.rotary_torque_body_nm)?),
            )?;
            if residual <= controls.relative_wrench_tolerance {
                angular_result = Some((trial, angular_mid, iteration, increment));
                break;
            }
            torque = desired;
        }
        let (angular_trial, angular_mid, rotational_iterations, rotation_increment_rad) =
            angular_result.ok_or(Error::Coupling {
                relative_residual: residual,
                tolerance: controls.relative_wrench_tolerance,
                iterations: controls.maximum_coupling_iterations,
            })?;
        let attitude_mid = midpoint_attitude(
            state.pose().orientation(),
            angular_trial.state_after.pose().orientation(),
        )?;
        let half_dt = positive(0.5 * dt)?;
        let acceleration = vector(
            self.gravity
                .acceleration_world()
                .scale(1.0 - self.fluid_mass_kg / p.mass())
                .add(applied.force_world.scale(1.0 / p.mass())),
        )?;
        let initial_position = state.pose().position_world();
        let ambient_before = self.ambient.velocity_at(initial_position)?;
        let c = self.resistance.translation_kg_s;
        let velocity_mid = if ambient_omega == Vec3::ZERO {
            // Retain the original uniform-flow elimination and operation order.
            let right = vector(
                initial_velocity
                    .sub(ambient_before)
                    .add(acceleration.scale(half_dt)),
            )?;
            let right_axes = apply(
                frame.ellipsoid_from_principal,
                attitude_mid.rotate_world_to_body(right),
            )?;
            let axes_mid = vector(Vec3::new(
                right_axes.x / positive(1.0 + half_dt * (c[0] / p.mass()))?,
                right_axes.y / positive(1.0 + half_dt * (c[1] / p.mass()))?,
                right_axes.z / positive(1.0 + half_dt * (c[2] / p.mass()))?,
            ))?;
            let relative_mid = vector(
                attitude_mid.rotate_body_to_world(apply(frame.principal_from_ellipsoid, axes_mid)?),
            )?;
            vector(relative_mid.add(ambient_before))?
        } else {
            let to_axes = |value| {
                apply(
                    frame.ellipsoid_from_principal,
                    attitude_mid.rotate_world_to_body(value),
                )
            };
            let axes_mid = ambient::midpoint_velocity_axes(
                c.map(|coefficient| half_dt * (coefficient / p.mass())),
                to_axes(ambient_omega)?.scale(half_dt),
                to_axes(vector(initial_velocity.add(acceleration.scale(half_dt)))?)?,
                to_axes(ambient_before)?,
            )?;
            vector(
                attitude_mid.rotate_body_to_world(apply(frame.principal_from_ellipsoid, axes_mid)?),
            )?
        };
        let predicted_position_mid = vector(initial_position.add(velocity_mid.scale(half_dt)))?;
        let predicted_load = self.load_at(
            &frame,
            attitude_mid,
            predicted_position_mid,
            velocity_mid,
            angular_mid,
        )?;
        let dynamics = integrator.step(
            state,
            p,
            Wrench {
                force_world: vector(predicted_load.wrench.force_world.add(applied.force_world))?,
                torque_body: vector(torque.add(applied.torque_body))?,
            },
            dt,
        )?;
        // Use actual endpoint averages from the native owner for work and heat.
        let actual_velocity_mid = vector(
            initial_velocity.scale(0.5).add(
                dynamics
                    .state_after
                    .center_of_mass_velocity_world(p)?
                    .scale(0.5),
            ),
        )?;
        let actual_omega_mid = vector(
            initial_omega.scale(0.5).add(
                p.angular_velocity_body_checked(dynamics.state_after.angular_momentum_body())?
                    .scale(0.5),
            ),
        )?;
        let actual_position_mid = vector(
            initial_position.add(
                dynamics
                    .state_after
                    .pose()
                    .position_world()
                    .sub(initial_position)
                    .scale(0.5),
            ),
        )?;
        let load_midpoint = self.load_at(
            &frame,
            attitude_mid,
            actual_position_mid,
            actual_velocity_mid,
            actual_omega_mid,
        )?;
        let force_residual = relative_residual(
            load_midpoint.drag_force_world_n,
            predicted_load.drag_force_world_n,
            norm(applied.force_world)?.max(norm(load_before.drag_force_world_n)?),
        )?;
        residual = residual.max(force_residual).max(relative_residual(
            load_midpoint.rotary_torque_body_nm,
            torque,
            norm(applied.torque_body)?.max(norm(load_before.rotary_torque_body_nm)?),
        )?);
        if residual > controls.relative_wrench_tolerance {
            return Err(Error::Coupling {
                relative_residual: residual,
                tolerance: controls.relative_wrench_tolerance,
                iterations: rotational_iterations,
            });
        }
        let load_after = self.state_load(&frame, dynamics.state_after)?;
        let viscous_heat_j = finite(dt * load_midpoint.dissipated_power_w)?;
        let background_work_j = finite(dt * load_midpoint.background_power_w)?;
        let applied_power = finite(
            applied.force_world.dot(actual_velocity_mid)
                + applied.torque_body.dot(actual_omega_mid),
        )?;
        let applied_work_j = finite(dt * applied_power)?;
        if load_midpoint.dissipated_power_w > 0.0 && viscous_heat_j == 0.0
            || load_midpoint.background_power_w != 0.0 && background_work_j == 0.0
            || applied_power != 0.0 && applied_work_j == 0.0
        {
            return Err(Error::Invalid(
                "ellipsoid work or heat vanished in conversion to energy",
            ));
        }
        let dx = dynamics
            .state_after
            .pose()
            .position_world()
            .sub(state.pose().position_world());
        let before = dynamics.diagnostics_before;
        let after = dynamics.diagnostics_after;
        let stored_energy_change_j = finite(
            (after.translational_kinetic_energy - before.translational_kinetic_energy)
                + (after.rotational_kinetic_energy - before.rotational_kinetic_energy)
                - (p.mass() - self.fluid_mass_kg) * self.gravity.acceleration_world().dot(dx),
        )?;
        let energy_defect_j =
            finite(stored_energy_change_j + viscous_heat_j - background_work_j - applied_work_j)?;
        if energy_defect_j.abs() > controls.mechanical.maximum_energy_defect_j {
            return Err(Error::Regime(CreepingFlowError::EnergyDefect {
                defect_j: energy_defect_j,
                maximum_j: controls.mechanical.maximum_energy_defect_j,
            }));
        }
        poll(&mut cancelled)?;
        Ok(EllipsoidStep {
            flow: CreepingStep {
                dynamics,
                load_before,
                load_midpoint,
                load_after,
                viscous_heat_j,
                background_work_j,
                applied_work_j,
                stored_energy_change_j,
                energy_defect_j,
                diffusion_ratio,
            },
            rotational_iterations,
            relative_wrench_residual: residual,
            rotation_increment_rad,
            ambient_rotation_increment_rad,
        })
    }
    fn timescales(&self, frame: &BodyFrame) -> Result<(f64, f64), Error> {
        let p = frame.properties;
        let moments = p.principal_inertia_body();
        let roots = [
            det::sqrt(moments.x),
            det::sqrt(moments.y),
            det::sqrt(moments.z),
        ];
        let mut fastest = 0.0_f64;
        for c in self.resistance.translation_kg_s {
            fastest = fastest.max(positive(c / p.mass())?);
        }
        // Gershgorin bound of I^-1/2 D I^-1/2, conservative under frame roundoff.
        for i in 0..3 {
            let mut rate = 0.0;
            for j in 0..3 {
                rate += (frame.rotary_principal.m[3 * i + j] / roots[i] / roots[j]).abs();
            }
            fastest = fastest.max(positive(rate)?);
        }
        // The imposed fluid field must diffuse across the body faster than it
        // turns, even when the body initially co-moves with it and has no slip.
        fastest = fastest.max(self.ambient.angular_speed_rad_s());
        let shortest = positive(1.0 / fastest)?;
        let ratio = positive(self.diffusion_time_s * fastest)?;
        if ratio > self.limits.maximum_diffusion_ratio {
            return Err(Error::Regime(CreepingFlowError::Unsteady {
                diffusion_ratio: ratio,
                maximum: self.limits.maximum_diffusion_ratio,
            }));
        }
        Ok((shortest, ratio))
    }
}
fn validate_controls(controls: EllipsoidStepControls) -> Result<(), Error> {
    let step = controls.mechanical;
    if !small_positive(step.maximum_relaxation_fraction, 1.0)
        || !step.duration_s.is_finite()
        || step.duration_s <= 0.0
        || !step.maximum_energy_defect_j.is_finite()
        || step.maximum_energy_defect_j <= 0.0
        || !(1..=64).contains(&controls.maximum_coupling_iterations)
        || !controls.relative_wrench_tolerance.is_finite()
        || !(64.0 * f64::EPSILON..=1e-8).contains(&controls.relative_wrench_tolerance)
        || !small_positive(controls.maximum_rotation_increment_rad, 0.5)
    {
        return Err(Error::Invalid(
            "invalid ellipsoid native-step, coupling or rotation controls",
        ));
    }
    Ok(())
}
fn midpoint_attitude(
    before: UnitQuaternion,
    after: UnitQuaternion,
) -> Result<UnitQuaternion, Error> {
    let a = before.components();
    let b = after.components();
    let dot: f64 = a.iter().zip(b).map(|(a, b)| *a * b).sum();
    let sign = if dot < 0.0 { -1.0 } else { 1.0 };
    Ok(UnitQuaternion::new(
        a[0] + sign * b[0],
        a[1] + sign * b[1],
        a[2] + sign * b[2],
        a[3] + sign * b[3],
    )?)
}
fn relative_residual(actual: Vec3, expected: Vec3, source_scale: f64) -> Result<f64, Error> {
    let scale = source_scale.max(norm(actual)?).max(norm(expected)?);
    let difference = norm(actual.sub(expected))?;
    if scale == 0.0 {
        return Ok(0.0);
    }
    finite(difference / scale)
}
fn apply(matrix: Mat3, v: Vec3) -> Result<Vec3, Error> {
    let result = matrix.apply(GeometryVector::new(v.x, v.y, v.z));
    vector(Vec3::new(result.x, result.y, result.z))
}
fn diagonal(diagonal: [f64; 3], v: Vec3) -> Result<Vec3, Error> {
    vector(Vec3::new(
        diagonal[0] * v.x,
        diagonal[1] * v.y,
        diagonal[2] * v.z,
    ))
}
fn small_positive(value: f64, maximum: f64) -> bool {
    value.is_finite() && value > 0.0 && value <= maximum
}
fn finite(value: f64) -> Result<f64, Error> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::Invalid("unrepresentable ellipsoid arithmetic"))
    }
}
fn positive(value: f64) -> Result<f64, Error> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(Error::Invalid("unrepresentable positive ellipsoid scale"))
    }
}
fn vector(value: Vec3) -> Result<Vec3, Error> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::Invalid("unrepresentable ellipsoid vector"))
    }
}
fn norm(value: Vec3) -> Result<f64, Error> {
    finite(value.x.hypot(value.y).hypot(value.z))
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
