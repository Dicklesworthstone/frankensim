//! Conservative spatial gravity through the original articulated motion owner.
//!
//! Each field sample uses the actual link COM pose and full central inertia.
//! COM forces and tidal torques are shifted to the original link origins for
//! ABA; the same loads supply whole-tree momentum and work. The fixed source
//! does not recoil, and finite-step energy defects remain measured quantities.

use super::{
    ArticulatedError, ArticulatedModel, EvaluatedStage, FreeFloatingMotionError,
    FreeFloatingMotionLoads, FreeFloatingMotionState, FreeFloatingMotionStep, Kinematics,
    MAX_FREE_FLOATING_STEP_LINKS, components, poll, step_evaluated, sum_products, validate_slice,
    validate_vec3,
};
use crate::gravity::{GravityError, RigidGravityField};
use crate::{MassProperties, Pose, UnitQuaternion};
use fs_ga::{GaError, Se3, Vec3, Wrench};
use fs_math::dd::Dd;
use std::borrow::Cow;

/// Atomic refusal by the original motion or conservative-field owner.
#[derive(Clone, Debug, PartialEq)]
pub enum ArticulatedGravityError {
    /// Model, state, load, duration, arithmetic, or cancellation refused.
    Motion(FreeFloatingMotionError),
    /// Native full mass properties, source frame, or field domain refused.
    Field(GravityError),
}

impl core::fmt::Display for ArticulatedGravityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Motion(error) => write!(f, "articulated gravity: {error}"),
            Self::Field(error) => write!(f, "articulated {error}"),
        }
    }
}

impl std::error::Error for ArticulatedGravityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(match self {
            Self::Motion(error) => error,
            Self::Field(error) => error,
        })
    }
}

impl From<FreeFloatingMotionError> for ArticulatedGravityError {
    fn from(error: FreeFloatingMotionError) -> Self {
        Self::Motion(error)
    }
}
impl From<ArticulatedError> for ArticulatedGravityError {
    fn from(error: ArticulatedError) -> Self {
        Self::Motion(error.into())
    }
}
impl From<GaError> for ArticulatedGravityError {
    fn from(error: GaError) -> Self {
        Self::Motion(error.into())
    }
}
impl From<GravityError> for ArticulatedGravityError {
    fn from(error: GravityError) -> Self {
        Self::Field(error)
    }
}

/// Per-link field loads at one configuration, suitable for native ABA/RNEA.
#[derive(Clone, Debug, PartialEq)]
pub struct ArticulatedGravityLoad {
    /// Force and moment about each original link origin, in that link's axes.
    pub body_wrenches: Vec<Wrench>,
    /// Sum of the same field's per-link conservative potentials, in joules.
    pub potential_energy_j: f64,
}

/// Native trajectory and a distinct conservative gravitational energy account.
#[derive(Clone, Debug, PartialEq)]
pub struct ArticulatedGravityStep {
    /// Original motion receipt. Its external work already includes gravity work.
    pub motion: FreeFloatingMotionStep,
    /// Field potential of the supplied physical link COM poses, in joules.
    pub potential_energy_before_j: f64,
    /// Field potential of the accepted physical link COM poses, in joules.
    pub potential_energy_after_j: f64,
    /// Sum of per-link endpoint changes before rounding the total, in joules.
    pub potential_energy_change_j: f64,
    /// Signed midpoint field work, already included in `motion.external_work`.
    pub gravity_work_j: f64,
    /// Held external-wrench work sampled independently of the gravity source.
    pub other_external_work_j: f64,
    /// Potential change plus gravity work: a quadrature defect, not dissipation.
    pub gravity_work_residual_j: f64,
    /// Actual kinetic plus potential change minus nongravitational work, in J.
    /// This measures integration error; it does not assert discrete conservation.
    pub mechanical_work_residual_j: f64,
}

/// Full native mass properties prepared once for an immutable original tree.
///
/// Each COM remains at its original material offset. The native rigid inertia
/// admission retains all tensor components and runs only during preparation;
/// stage evaluation neither rebuilds principal frames nor changes tree inertia.
#[derive(Debug)]
pub struct PreparedArticulatedGravity<'a> {
    model: &'a ArticulatedModel,
    properties: Vec<MassProperties>,
}

impl<'a> PreparedArticulatedGravity<'a> {
    /// Bind the original tree, retaining each mass, COM, and full central inertia.
    ///
    /// Refuses more than the native motion link bound, nonphysical rigid-body
    /// principal moments (including triangle inequalities), or cancellation.
    pub fn new(
        model: &'a ArticulatedModel,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ArticulatedGravityError> {
        poll(&mut cancelled)?;
        if model.link_count() > MAX_FREE_FLOATING_STEP_LINKS {
            return Err(FreeFloatingMotionError::TooManyLinks {
                got: model.link_count(),
                maximum: MAX_FREE_FLOATING_STEP_LINKS,
            }
            .into());
        }
        let mut properties = Vec::with_capacity(model.link_count());
        for (i, link) in model.links().iter().enumerate() {
            if i % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let inertia = link.inertia();
            properties.push(
                MassProperties::from_inertia_tensor(
                    inertia.mass(),
                    crate::Vec3::ZERO,
                    inertia.inertia_com(),
                )
                .map_err(GravityError::from)?,
            );
        }
        poll(&mut cancelled)?;
        Ok(Self { model, properties })
    }

    /// Sample one stationary field at the supplied original link poses.
    ///
    /// The kinematics must describe this bound model in original link order.
    /// The receipt can drive prescribed or free native dynamics with their
    /// separate uniform-gravity input set to zero. A frame ID identifies the
    /// same inertial coordinates used by the field; no frame transform is inferred.
    pub fn load(
        &self,
        field: &impl RigidGravityField,
        frame: u64,
        kinematics: &Kinematics,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ArticulatedGravityLoad, ArticulatedGravityError> {
        poll(&mut cancelled)?;
        super::super::require_len(
            "articulated_gravity.world_from_link",
            kinematics.world_from_link.len(),
            self.model.link_count(),
        )?;
        let mut body_wrenches = Vec::with_capacity(self.model.link_count());
        let mut potential = Dd::ZERO;
        for (i, pose) in kinematics.world_from_link.iter().enumerate() {
            if i % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let load = self.sample_link(field, frame, i, *pose)?;
            body_wrenches.push(load.at_origin);
            potential = potential + Dd::from_f64(load.potential);
        }
        validate_slice(&[potential.to_f64()], "articulated_gravity.potential")?;
        poll(&mut cancelled)?;
        Ok(ArticulatedGravityLoad {
            body_wrenches,
            potential_energy_j: potential.to_f64(),
        })
    }

    /// Advance with field loads sampled at the native initial and midpoint poses.
    ///
    /// This reuses the original O(links) Lie-midpoint/ABA momentum owner and its
    /// actual endpoint kinematics. Held wrenches and actuator efforts compose
    /// normally. `loads.gravity_world` must be zero so gravity has one owner.
    /// The field is stationary and conservative as required by `RigidGravityField`;
    /// source recoil, mutual gravity, contacts, and constraints are not inferred.
    ///
    /// The endpoint is also sampled for physical field admission and potential.
    /// Any field, dynamics, work, or cancellation failure returns no candidate.
    /// Custom field callbacks own their side effects and must bound their work.
    pub fn step(
        &self,
        field: &impl RigidGravityField,
        frame: u64,
        state: &FreeFloatingMotionState,
        loads: FreeFloatingMotionLoads<'_>,
        duration: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ArticulatedGravityStep, ArticulatedGravityError> {
        poll(&mut cancelled)?;
        if loads.gravity_world != Vec3::default() {
            return Err(FreeFloatingMotionError::InvalidLoad(
                "field gravity requires zero separate uniform gravity",
            )
            .into());
        }
        let mut potential_before = Vec::with_capacity(self.model.link_count());
        let mut field_error = None;
        let mut gravity_power = Dd::ZERO;
        let mut other_power = Dd::ZERO;
        let evaluated = step_evaluated(
            self.model,
            state,
            loads,
            duration,
            &mut cancelled,
            |map, base, cancelled| {
                let result = (|| -> Result<_, ArticulatedGravityError> {
                    let mut accumulated = map.loads_accumulated(base, loads, cancelled)?;
                    other_power = accumulated.body_wrench_power;
                    gravity_power = Dd::ZERO;
                    let record_start = potential_before.is_empty();
                    let mut wrenches = loads.external_body_wrench.to_vec();
                    for (i, wrench) in wrenches.iter_mut().enumerate() {
                        if i % 64 == 0 {
                            poll(cancelled)?;
                        }
                        let relative = map.shape.world_from_link[i];
                        let pose = base.world_from_base.compose(relative)?;
                        let load = self.sample_link(field, frame, i, pose)?;
                        let center = self.model.links()[i].inertia().center_of_mass();
                        let arm = base
                            .world_from_base
                            .rotation()
                            .rotate(relative.transform_point(center)? - map.center_base)?;
                        // Reduce about the actual system COM before translating
                        // into the world, retaining the material COM lever arm.
                        accumulated.add_force(load.force_world);
                        accumulated.add_torque(pose.rotation().rotate(load.torque_com_body)?);
                        accumulated.add_torque_cross(arm, load.force_world);
                        let velocity = relative
                            .inverse()?
                            .adjoint()
                            .apply_twist(base.twist_body)
                            .plus(map.shape.body_twist[i]);
                        let com_velocity = velocity.linear + velocity.angular.cross(center);
                        gravity_power = gravity_power
                            + sum_products(
                                components(load.at_origin.force),
                                components(com_velocity),
                            )
                            + sum_products(
                                components(load.torque_com_body),
                                components(velocity.angular),
                            );
                        *wrench = Wrench::new(
                            wrench.torque + load.at_origin.torque,
                            wrench.force + load.at_origin.force,
                        );
                        accumulated.add_storage(load.potential);
                        if record_start {
                            potential_before.push(load.potential);
                        }
                    }
                    accumulated.body_wrench_power = accumulated.body_wrench_power + gravity_power;
                    let load = accumulated.finish();
                    validate_vec3(load.force_world, "articulated_gravity.stage_force")?;
                    validate_vec3(load.torque_com_world, "articulated_gravity.stage_torque")?;
                    validate_slice(
                        &[load.body_wrench_power, load.stored_energy],
                        "articulated_gravity.stage_energy",
                    )?;
                    Ok(EvaluatedStage {
                        body_wrench: Cow::Owned(wrenches),
                        load,
                    })
                })();
                result.map_err(|error| match error {
                    ArticulatedGravityError::Motion(error) => error,
                    ArticulatedGravityError::Field(error) => {
                        field_error = Some(error);
                        FreeFloatingMotionError::InvalidLoad("conservative gravity field")
                    }
                })
            },
            None,
        );
        let evaluated = evaluated.map_err(|error| match field_error {
            Some(error) => ArticulatedGravityError::Field(error),
            None => ArticulatedGravityError::Motion(error),
        })?;
        let mut potential_after = Dd::ZERO;
        let mut potential_change = Dd::ZERO;
        for (i, (pose, before)) in evaluated
            .kinematics_after
            .world_from_link
            .iter()
            .zip(potential_before)
            .enumerate()
        {
            if i % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let after = Dd::from_f64(self.sample_link(field, frame, i, *pose)?.potential);
            potential_after = potential_after + after;
            potential_change = potential_change + (after - Dd::from_f64(before));
        }
        let dt = Dd::from_f64(duration);
        let gravity_work = dt * gravity_power;
        let other_work = dt * other_power;
        let motion = evaluated.step;
        let gravity_residual = potential_change + gravity_work;
        let mechanical_residual = Dd::from_f64(motion.kinetic_energy_change) + potential_change
            - other_work
            - Dd::from_f64(motion.actuator_work);
        validate_slice(
            &[
                potential_after.to_f64(),
                potential_change.to_f64(),
                gravity_work.to_f64(),
                other_work.to_f64(),
                gravity_residual.to_f64(),
                mechanical_residual.to_f64(),
            ],
            "articulated_gravity.endpoint_energy",
        )?;
        poll(&mut cancelled)?;
        Ok(ArticulatedGravityStep {
            motion,
            potential_energy_before_j: evaluated.stored_energy_before,
            potential_energy_after_j: potential_after.to_f64(),
            potential_energy_change_j: potential_change.to_f64(),
            gravity_work_j: gravity_work.to_f64(),
            other_external_work_j: other_work.to_f64(),
            gravity_work_residual_j: gravity_residual.to_f64(),
            mechanical_work_residual_j: mechanical_residual.to_f64(),
        })
    }

    fn sample_link(
        &self,
        field: &impl RigidGravityField,
        frame: u64,
        link: usize,
        pose: Se3,
    ) -> Result<LinkGravity, ArticulatedGravityError> {
        let center = self.model.links()[link].inertia().center_of_mass();
        let position = pose.transform_point(center)?;
        let rotation = pose.rotation();
        let q = rotation.as_quat();
        let orientation = UnitQuaternion::from_canonical_components([q.w, q.x, q.y, q.z])
            .map_err(GravityError::from)?;
        let com_pose = Pose::new(
            crate::Vec3::new(position.x, position.y, position.z),
            orientation,
        )
        .map_err(GravityError::from)?;
        let load = field.load(frame, com_pose, self.properties[link])?;
        load.wrench.validate().map_err(GravityError::from)?;
        if !load.potential_energy_j.is_finite() {
            return Err(GravityError::Invalid("nonfinite articulated gravity potential").into());
        }
        let force_world = geometric(load.wrench.force_world);
        let force = rotation.inverse().rotate(force_world)?;
        let torque_com_body = geometric(load.wrench.torque_body);
        let at_origin = Wrench::new(torque_com_body + center.cross(force), force);
        super::super::validate_wrench(at_origin, "articulated_gravity.link_wrench")?;
        Ok(LinkGravity {
            at_origin,
            force_world,
            torque_com_body,
            potential: load.potential_energy_j,
        })
    }
}

struct LinkGravity {
    at_origin: Wrench,
    force_world: Vec3,
    torque_com_body: Vec3,
    potential: f64,
}

fn geometric(value: crate::Vec3) -> Vec3 {
    Vec3::new(value.x, value.y, value.z)
}
