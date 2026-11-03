//! Velocity-controlled separation of finite rigid components.
//!
//! A prescribed relative twist fixes every component's final motion relative
//! to one reference component's FINAL rigid velocity field. All components
//! therefore share translation and rotation recoil. Their full intrinsic
//! inertias participate alongside orbital angular momentum.
//!
//! This owner derives equivalent internal force/couple impulses and delegates
//! all state construction and work accounting to the original assembly release.
//! It does not discover contacts, locate real actuators, model fracture, or
//! supply the required mechanical energy.

use super::capture::relative_momentum;
use super::{
    AssemblyControls, AssemblyError, AssemblyPart, AssemblyRelease, RigidAssembly,
    SeparationImpulse, VectorSum, maximum, poll, vector,
};
use crate::{MassProperties, RigidBodyState, Vec3};

/// One component's final motion relative to a chosen reference component.
///
/// Both vectors use the ORIGINAL assembly axes, not either component's axes.
/// The reference component's entry must be exactly ZERO.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RelativeSeparationTwist {
    /// COM velocity minus the reference's final rigid velocity evaluated at
    /// this same COM [m/s]. This is not merely a difference of COM velocities.
    pub relative_velocity_assembly_m_s: Vec3,
    /// Component angular velocity minus reference angular velocity [rad/s].
    pub relative_angular_velocity_assembly_rad_s: Vec3,
}
impl RelativeSeparationTwist {
    /// Follow the reference's entire final rigid velocity field.
    pub const ZERO: Self = Self {
        relative_velocity_assembly_m_s: Vec3::ZERO,
        relative_angular_velocity_assembly_rad_s: Vec3::ZERO,
    };
}

/// Native event budgets, actual relative-motion gates, and available energy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RelativeSeparationControls {
    /// Original part/impulse caps and momentum/energy balance tolerances.
    pub assembly: AssemblyControls,
    /// Maximum component of actual relative COM-field velocity error [m/s].
    pub relative_velocity_tolerance_m_s: f64,
    /// Maximum component of actual relative angular-velocity error [rad/s].
    pub relative_angular_velocity_tolerance_rad_s: f64,
    /// Available separation work [J], finite and nonnegative. The composing
    /// energy source must debit the accepted native supplied work.
    pub maximum_supplied_work_j: f64,
}

/// One complete finite-component event, with the original native endpoints.
#[derive(Clone, Debug, PartialEq)]
pub struct RelativeSeparationResult {
    /// All native component states, joint separation work, and actual balances.
    pub release: AssemblyRelease,
    /// Equivalent nonzero star impulses against the reference component.
    /// Force acts at the other component's COM; its intrinsic spin impulse is
    /// a free couple. These are a wrench decomposition, not actuator geometry.
    pub impulses: Vec<SeparationImpulse>,
    /// Largest actual relative-velocity residual component [m/s].
    pub maximum_relative_velocity_residual_m_s: f64,
    /// Largest actual relative-angular-velocity residual component [rad/s].
    pub maximum_relative_angular_velocity_residual_rad_s: f64,
}

/// Native refusal, an unrealized velocity boundary condition, or insufficient
/// source energy. No variant publishes partially separated components.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RelativeSeparationError {
    /// Original admission, arithmetic, budget, balance, or cancellation refusal.
    Assembly(AssemblyError),
    /// Actual canonical endpoints do not realize the requested linear motion.
    RelativeVelocity {
        /// Original component index, or retained group index for a partition.
        part_index: usize,
        /// Actual minus requested velocity in original assembly axes [m/s].
        residual_assembly_m_s: Vec3,
        /// Declared maximum absolute component [m/s].
        tolerance_m_s: f64,
    },
    /// Actual canonical endpoints do not realize the requested spin difference.
    RelativeAngularVelocity {
        /// Original component index, or retained group index for a partition.
        part_index: usize,
        /// Actual minus requested angular velocity in assembly axes [rad/s].
        residual_assembly_rad_s: Vec3,
        /// Declared maximum absolute component [rad/s].
        tolerance_rad_s: f64,
    },
    /// Native joint impulse work exceeds its available mechanical source.
    SourceWork {
        /// Required native separation work [J].
        required_j: f64,
        /// Declared available work [J].
        maximum_j: f64,
    },
}
impl From<AssemblyError> for RelativeSeparationError {
    fn from(error: AssemblyError) -> Self {
        Self::Assembly(error)
    }
}
impl core::fmt::Display for RelativeSeparationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "relative rigid separation: {self:?}")
    }
}
impl std::error::Error for RelativeSeparationError {}

impl RigidAssembly {
    /// Release finite components at specified final relative linear/spin rates.
    ///
    /// There must be one target per original component; reference_part selects
    /// an existing component with a ZERO target. For centered COM arms r_i,
    /// intrinsic inertias J_i in assembly axes, total mass M and full composite
    /// inertia I, uncentered conservation gives
    /// dv = -sum(m_i*c_i)/M,
    /// dw = -I^-1 sum(r_i cross m_i*c_i + J_i*eta_i).
    /// Component impulses are m_i*(dv+dw cross r_i+c_i) and
    /// J_i*(dw+eta_i). A single original inertia action suffices; there is no
    /// dense 6N mobility matrix or sequence of independent pair solves.
    ///
    /// Computation uses an equivalent common-twist subtraction before these
    /// reductions. The largest-mass component anchors linear data and the
    /// largest intrinsic principal moment anchors angular data, retaining small
    /// recoil when a component dominates the total inventory or spin inertia.
    /// Input order breaks pivot ties and fixes compensated reduction order.
    ///
    /// Nonzero solved force/couple pairs consume the native impulse budget.
    /// A zero target can still require an impulse to follow the common recoil.
    /// Native release includes all shared-reference work cross terms. Its
    /// momentum/work gates, available source work, actual endpoint twist gates
    /// and cancellation must all pass before any endpoints are returned.
    ///
    /// Targets describe an instantaneous ideal internal actuation event. The
    /// caller owns its mechanical feasibility, contact/connectivity decision,
    /// actuator locations, finite-time force law and source-energy debit.
    pub fn release_relative_twists(
        &self,
        parent: RigidBodyState,
        reference_part: usize,
        targets: &[RelativeSeparationTwist],
        controls: RelativeSeparationControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RelativeSeparationResult, RelativeSeparationError> {
        poll(&mut cancelled)?;
        controls.assembly.validate(self.parts.len(), 0)?;
        if targets.len() != self.parts.len() || reference_part >= self.parts.len() {
            return Err(AssemblyError::Invalid(
                "one relative twist per part and an existing reference required",
            )
            .into());
        }
        if [
            controls.relative_velocity_tolerance_m_s,
            controls.relative_angular_velocity_tolerance_rad_s,
            controls.maximum_supplied_work_j,
        ]
        .into_iter()
        .any(|v| !v.is_finite() || v < 0.0)
        {
            return Err(AssemblyError::Invalid(
                "finite nonnegative relative-motion tolerances and source work required",
            )
            .into());
        }
        let mut mass_anchor = 0;
        let mut spin_anchor = 0;
        for (index, (part, target)) in self.parts.iter().zip(targets).enumerate() {
            poll(&mut cancelled)?;
            vector(target.relative_velocity_assembly_m_s)?;
            vector(target.relative_angular_velocity_assembly_rad_s)?;
            if part.properties.mass() > self.parts[mass_anchor].properties.mass() {
                mass_anchor = index;
            }
            if maximum(part.properties.principal_inertia_body())
                > maximum(self.parts[spin_anchor].properties.principal_inertia_body())
            {
                spin_anchor = index;
            }
        }
        if targets[reference_part] != RelativeSeparationTwist::ZERO {
            return Err(AssemblyError::Invalid(
                "reference component must have a zero relative twist",
            )
            .into());
        }

        let anchor_position = self.parts[mass_anchor].pose_in_reference.position_world();
        let anchor_velocity = targets[mass_anchor].relative_velocity_assembly_m_s;
        let anchor_spin = targets[spin_anchor].relative_angular_velocity_assembly_rad_s;
        let centered = |part: &AssemblyPart, target: &RelativeSeparationTwist| {
            // b_i = (c_i-c_anchor) - eta_anchor cross (r_i-r_anchor).
            // d_i = eta_i-eta_anchor. This is the same physical constraint
            // chart, with one common rigid twist removed before reduction.
            let arm = vector(part.pose_in_reference.position_world().sub(anchor_position))?;
            let velocity = vector(
                vector(target.relative_velocity_assembly_m_s.sub(anchor_velocity))?
                    .sub(vector(anchor_spin.cross(arm))?),
            )?;
            let spin = vector(
                target
                    .relative_angular_velocity_assembly_rad_s
                    .sub(anchor_spin),
            )?;
            Ok::<_, AssemblyError>((velocity, spin))
        };
        let mut translation = VectorSum::default();
        let mut angular = VectorSum::default();
        for (part, target) in self.parts.iter().zip(targets) {
            poll(&mut cancelled)?;
            let (relative, spin) = centered(part, target)?;
            let mass = part.properties.mass();
            translation.add(vector(relative.scale(-mass / self.properties.mass()))?);
            angular.add(vector(
                part.pose_in_reference
                    .position_world()
                    .cross(vector(relative.scale(mass))?),
            )?);
            angular.add(spin_in_assembly(part, spin)?);
        }
        let translation = vector(translation.value())?;
        let rotation = vector(
            self.properties
                .angular_velocity_body_checked(vector(angular.value())?)
                .map_err(AssemblyError::from)?
                .scale(-1.0),
        )?;
        let mut impulses = Vec::with_capacity(self.parts.len() - 1);
        for (index, (part, target)) in self.parts.iter().zip(targets).enumerate() {
            poll(&mut cancelled)?;
            if index == reference_part {
                continue;
            }
            let (relative, spin) = centered(part, target)?;
            let mut dv = VectorSum::default();
            dv.add(relative);
            dv.add(translation);
            dv.add(vector(
                rotation.cross(part.pose_in_reference.position_world()),
            )?);
            let impulse = vector(vector(dv.value())?.scale(part.properties.mass()))?;
            let couple = spin_in_assembly(part, vector(rotation.add(spin))?)?;
            if impulse == Vec3::ZERO && couple == Vec3::ZERO {
                continue;
            }
            controls
                .assembly
                .validate(self.parts.len(), impulses.len() + 1)?;
            impulses.push(SeparationImpulse {
                part_a: index,
                part_b: reference_part,
                point_assembly_m: part.pose_in_reference.position_world(),
                impulse_on_a_assembly_ns: impulse,
                couple_on_a_assembly_nms: couple,
            });
        }
        let release = self.release(parent, &impulses, controls.assembly, &mut cancelled)?;
        if release.supplied_work_j > controls.maximum_supplied_work_j {
            return Err(RelativeSeparationError::SourceWork {
                required_j: release.supplied_work_j,
                maximum_j: controls.maximum_supplied_work_j,
            });
        }

        // Audit actual stored momenta and actual represented component poses.
        // Removing one common boost before division cannot hide an unrepresented
        // velocity jump behind two separately rounded laboratory velocities.
        let boost = parent
            .center_of_mass_velocity_world(self.properties)
            .map_err(AssemblyError::from)?;
        let reference = release.bodies[reference_part];
        let reference_properties = self.parts[reference_part].properties;
        let reference_velocity = velocity_in_chart(reference, reference_properties, boost)?;
        let reference_spin = angular_velocity_world(reference, reference_properties)?;
        let orientation = parent.pose().orientation();
        let mut maximum_linear = 0.0_f64;
        let mut maximum_angular = 0.0_f64;
        for (index, ((part, target), &body)) in self
            .parts
            .iter()
            .zip(targets)
            .zip(&release.bodies)
            .enumerate()
        {
            poll(&mut cancelled)?;
            if index == reference_part {
                continue;
            }
            let arm = vector(
                body.pose()
                    .position_world()
                    .sub(reference.pose().position_world()),
            )?;
            let relative = vector(
                velocity_in_chart(body, part.properties, boost)?
                    .sub(reference_velocity)
                    .sub(vector(reference_spin.cross(arm))?),
            )?;
            let linear = vector(
                orientation
                    .rotate_world_to_body(relative)
                    .sub(target.relative_velocity_assembly_m_s),
            )?;
            let angular = vector(
                orientation
                    .rotate_world_to_body(vector(
                        angular_velocity_world(body, part.properties)?.sub(reference_spin),
                    )?)
                    .sub(target.relative_angular_velocity_assembly_rad_s),
            )?;
            maximum_linear = maximum_linear.max(maximum(linear));
            maximum_angular = maximum_angular.max(maximum(angular));
            if maximum(linear) > controls.relative_velocity_tolerance_m_s {
                return Err(RelativeSeparationError::RelativeVelocity {
                    part_index: index,
                    residual_assembly_m_s: linear,
                    tolerance_m_s: controls.relative_velocity_tolerance_m_s,
                });
            }
            if maximum(angular) > controls.relative_angular_velocity_tolerance_rad_s {
                return Err(RelativeSeparationError::RelativeAngularVelocity {
                    part_index: index,
                    residual_assembly_rad_s: angular,
                    tolerance_rad_s: controls.relative_angular_velocity_tolerance_rad_s,
                });
            }
        }
        poll(&mut cancelled)?;
        Ok(RelativeSeparationResult {
            release,
            impulses,
            maximum_relative_velocity_residual_m_s: maximum_linear,
            maximum_relative_angular_velocity_residual_rad_s: maximum_angular,
        })
    }
}

fn spin_in_assembly(part: &AssemblyPart, spin: Vec3) -> Result<Vec3, AssemblyError> {
    let axes = part.pose_in_reference.orientation();
    let local = vector(axes.rotate_world_to_body(spin))?;
    let momentum = part.properties.angular_momentum_body_checked(local)?;
    vector(axes.rotate_body_to_world(momentum))
}
fn velocity_in_chart(
    body: RigidBodyState,
    properties: MassProperties,
    boost: Vec3,
) -> Result<Vec3, AssemblyError> {
    let mass = properties.mass();
    let momentum = relative_momentum(body.linear_momentum_world(), mass, boost)?;
    vector(Vec3::new(
        momentum.x / mass,
        momentum.y / mass,
        momentum.z / mass,
    ))
}
fn angular_velocity_world(
    body: RigidBodyState,
    properties: MassProperties,
) -> Result<Vec3, AssemblyError> {
    let spin = properties.angular_velocity_body_checked(body.angular_momentum_body())?;
    vector(body.pose().orientation().rotate_body_to_world(spin))
}
