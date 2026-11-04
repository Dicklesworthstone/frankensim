//! One selected axial, frictionless rigid-body/piston impact.
//!
//! The original rigid point mobility and guided wall mass supply the existing
//! Newton contact law. Gas pressure is finite and supplies no zero-time
//! impulse. Positions, source gas, thermal history and prior port windows stay
//! unchanged. This is not contact search, friction, a shock/heat model, or an
//! unconstrained six-degree-of-freedom model of the translating wall.

use super::{PistonImpulseControls, PistonImpulseError, PistonImpulseReceipt, poll};
use crate::cavity_phs::piston::{GasPiston, PistonGas, PistonSupport};
use fs_contact::impact::{
    NormalImpactError, NormalImpactInput, NormalImpactResult, solve_normal_impact,
};
use fs_math::dd::Dd;
use fs_mbd::{DynamicsError, ImpulseReceipt, MassProperties, RigidBodyState, Vec3};

/// Frozen material-point witness on the piston axis, on either face.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidPistonImpactContact {
    /// Native rigid witness relative to its COM, in body coordinates [m].
    pub rigid_arm_body_m: Vec3,
    /// World origin of this face when the native coordinate V/A is zero [m].
    pub piston_axis_origin_world_m: Vec3,
    /// Nonzero world direction of increasing native gas volume; normalized here.
    pub expansion_axis_world: Vec3,
    /// Exactly +1 or -1; normal from the selected piston face to the rigid body.
    pub normal_sign: f64,
    /// Native normal-law proximity activation [m]; no penetration repair.
    pub activation_tolerance_m: f64,
    /// Explicit Newton coefficient in [0, 1].
    pub restitution: f64,
}

/// Event limits; the same momentum precision also bounds restitution closure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidPistonImpactControls {
    /// Original wall momentum, velocity and energy admission.
    pub piston: PistonImpulseControls,
    /// Maximum lateral separation of the witnesses from the declared axis [m].
    pub geometry_tolerance_m: f64,
    /// Maximum accepted native rigid kinetic energy [J].
    pub maximum_rigid_kinetic_energy_j: f64,
}
impl Default for RigidPistonImpactControls {
    fn default() -> Self {
        let piston = PistonImpulseControls::default();
        Self {
            maximum_rigid_kinetic_energy_j: piston.maximum_kinetic_energy_j,
            piston,
            geometry_tolerance_m: 1e-9,
        }
    }
}

/// Both native momentum candidates and independently audited event work.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidPistonImpactReceipt {
    /// Original contact law; loss is unallocated, never automatically gas heat.
    pub contact: NormalImpactResult,
    /// Guided wall impulse and any fixed-support reaction.
    pub piston: PistonImpulseReceipt,
    /// Original rigid point impulse, including full-inertia rotational recoil.
    pub rigid: ImpulseReceipt,
    /// Actual unchanged face point, from the original V/A coordinate [m].
    pub piston_face_world_m: Vec3,
    /// Admitted unit world normal, pointing from piston to rigid witness.
    pub normal_world: Vec3,
    /// Signed selected-witness gap [m].
    pub gap_m: f64,
    /// Positive values separate the original witnesses [m/s].
    pub relative_normal_speed_before_m_s: f64,
    /// Relative speed measured again on both actual candidates [m/s].
    pub relative_normal_speed_after_m_s: f64,
    /// Complete rigid-point plus guided-wall inverse effective mass [1/kg].
    pub inverse_effective_mass_kg_inv: f64,
    /// Polarized actual kinetic change of both native receivers [J].
    pub kinetic_energy_change_j: f64,
    /// Actual kinetic change plus the native contact loss [J].
    pub energy_residual_j: f64,
    /// Actual endpoint relative speed minus the Newton target [m/s].
    pub restitution_residual_m_s: f64,
    /// Largest represented/requested momentum defect in either receiver.
    pub maximum_relative_momentum_residual: f64,
}

/// A refusal publishes neither receiver and preserves every gas/port history.
#[derive(Debug)]
pub enum RigidPistonImpactError {
    /// Original wall admission, geometry or event energy refusal.
    Piston(PistonImpulseError),
    /// Original full-inertia rigid owner refusal.
    Rigid(DynamicsError),
    /// Original scalar contact law refusal.
    Contact(NormalImpactError),
    /// Unrepresented rigid momentum; coordinates 0..3 linear, 3..6 body angular.
    Momentum {
        coordinate: usize,
        relative: f64,
        tolerance: f64,
    },
    /// Actual velocities do not close the original restitution law.
    Restitution {
        residual_m_s: f64,
        tolerance_m_s: f64,
    },
}
impl From<PistonImpulseError> for RigidPistonImpactError {
    fn from(error: PistonImpulseError) -> Self {
        Self::Piston(error)
    }
}
impl From<DynamicsError> for RigidPistonImpactError {
    fn from(error: DynamicsError) -> Self {
        Self::Rigid(error)
    }
}
impl From<NormalImpactError> for RigidPistonImpactError {
    fn from(error: NormalImpactError) -> Self {
        Self::Contact(error)
    }
}
impl core::fmt::Display for RigidPistonImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native rigid/piston impact refused: {self:?}")
    }
}
impl std::error::Error for RigidPistonImpactError {}

impl<G: PistonGas> GasPiston<G> {
    /// Resolve one selected normal impact against the original guided wall.
    ///
    /// The witness must lie on the declared face axis. Off-centre rigid arms
    /// retain full rotational response. The piston remains a one-coordinate
    /// constrained wall, not an inferred rotational body. Both native states
    /// publish together after momentum, energy, restitution and cancellation
    /// admission. No material source is queried and no gas heat is credited.
    pub fn apply_rigid_normal_impact(
        &mut self,
        rigid: &mut RigidBodyState,
        properties: MassProperties,
        contact: RigidPistonImpactContact,
        controls: RigidPistonImpactControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RigidPistonImpactReceipt, RigidPistonImpactError> {
        poll(&mut cancelled)?;
        if !contact.piston_axis_origin_world_m.is_finite()
            || (contact.normal_sign != -1.0 && contact.normal_sign != 1.0)
            || !(controls.geometry_tolerance_m > 0.0 && controls.geometry_tolerance_m.is_finite())
            || !(controls.maximum_rigid_kinetic_energy_j > 0.0
                && controls.maximum_rigid_kinetic_energy_j.is_finite())
        {
            return Err(invalid(
                "finite axial witness and positive rigid-impact limits required",
            ));
        }
        let mobility = rigid.directional_effective_mass(
            properties,
            contact.rigid_arm_body_m,
            contact.expansion_axis_world.scale(contact.normal_sign),
        )?;
        let normal = mobility.unit_direction_world;
        let axis = normal.scale(contact.normal_sign);
        let face = contact
            .piston_axis_origin_world_m
            .add(axis.scale(self.position_m()));
        let before_point = rigid.point_kinematics(properties, contact.rigid_arm_body_m)?;
        // Difference native relative coordinates before rounding world points.
        // A common large translation must not erase the actual gap/lever arm.
        let origin = components(contact.piston_axis_origin_world_m);
        let center = components(rigid.pose().position_world());
        let arm = components(before_point.arm_world);
        let direction = components(axis);
        let separation: [Dd; 3] = core::array::from_fn(|i| {
            Dd::from_f64(center[i]) - Dd::from_f64(origin[i]) + Dd::from_f64(arm[i])
                - Dd::from_f64(direction[i]) * Dd::from_f64(self.position_m())
        });
        if !face.is_finite() {
            return Err(invalid("piston contact position is unrepresentable"));
        }
        let gap_dd = separation
            .into_iter()
            .zip(components(normal))
            .fold(Dd::ZERO, |sum, (distance, normal)| {
                sum + distance * Dd::from_f64(normal)
            });
        let gap = checked(gap_dd)?;
        let lateral: [f64; 3] = [
            checked(separation[0] - Dd::from_f64(normal.x) * gap_dd)?,
            checked(separation[1] - Dd::from_f64(normal.y) * gap_dd)?,
            checked(separation[2] - Dd::from_f64(normal.z) * gap_dd)?,
        ];
        if fs_math::det::hypot(fs_math::det::hypot(lateral[0], lateral[1]), lateral[2])
            > controls.geometry_tolerance_m
        {
            return Err(invalid("rigid witness is outside the selected piston axis"));
        }
        let wall_inverse_mass = if self.support() == PistonSupport::Fixed {
            0.0
        } else {
            self.mass_kg().recip()
        };
        let inverse_mass = checked(
            Dd::from_f64(mobility.inverse_mass_kg_inverse) + Dd::from_f64(wall_inverse_mass),
        )?;
        let before = dot(before_point.point_velocity_world, normal)
            - Dd::from_f64(contact.normal_sign) * Dd::from_f64(self.velocity_m_s());
        let before_value = checked(before)?;
        let law = solve_normal_impact(NormalImpactInput {
            gap,
            activation_tolerance: contact.activation_tolerance_m,
            relative_normal_speed: before_value,
            inverse_effective_mass: inverse_mass,
            restitution: contact.restitution,
        })?;
        let impulse_world = normal.scale(law.normal_impulse);
        for (direction, applied) in components(normal)
            .into_iter()
            .zip(components(impulse_world))
        {
            if direction != 0.0 && law.normal_impulse != 0.0 && applied == 0.0 {
                return Err(invalid("rigid impulse projection loses nonzero momentum"));
            }
        }
        let piston =
            self.stage_impulse(-contact.normal_sign * law.normal_impulse, controls.piston)?;
        let rigid_jump = rigid.apply_impulse_at_body_point(
            properties,
            contact.rigid_arm_body_m,
            impulse_world,
        )?;
        if rigid_jump.work.kinetic_energy_after_j > controls.maximum_rigid_kinetic_energy_j {
            return Err(PistonImpulseError::Budget {
                what: "rigid kinetic energy",
                value: rigid_jump.work.kinetic_energy_after_j,
                limit: controls.maximum_rigid_kinetic_energy_j,
            }
            .into());
        }
        let mut worst = piston.relative_momentum_residual;
        for (group, (p0, p1, requested)) in [
            (
                rigid_jump.state_before.linear_momentum_world(),
                rigid_jump.state_after.linear_momentum_world(),
                rigid_jump.impulse_world,
            ),
            (
                rigid_jump.state_before.angular_momentum_body(),
                rigid_jump.state_after.angular_momentum_body(),
                rigid_jump.angular_impulse_body,
            ),
        ]
        .into_iter()
        .enumerate()
        {
            for (coordinate, ((old, next), applied)) in components(p0)
                .into_iter()
                .zip(components(p1))
                .zip(components(requested))
                .enumerate()
            {
                let change = Dd::from_f64(next) - Dd::from_f64(old);
                let scale = checked(change)?.abs().max(applied.abs());
                let defect = checked(change - Dd::from_f64(applied))?;
                let relative = if scale == 0.0 {
                    0.0
                } else {
                    defect.abs() / scale
                };
                if relative > controls.piston.momentum_relative_tolerance {
                    return Err(RigidPistonImpactError::Momentum {
                        coordinate: 3 * group + coordinate,
                        relative,
                        tolerance: controls.piston.momentum_relative_tolerance,
                    });
                }
                worst = worst.max(relative);
            }
        }
        let after = dot(
            rigid_jump.point_kinematics_after.point_velocity_world,
            normal,
        ) - Dd::from_f64(contact.normal_sign) * Dd::from_f64(piston.velocity_m_s[1]);
        let after_value = checked(after)?;
        let restitution_residual = checked(after - Dd::from_f64(law.relative_normal_speed_after))?;
        let velocity_tolerance = controls.piston.momentum_relative_tolerance
            * before_value
                .abs()
                .max(after_value.abs())
                .max(law.relative_normal_speed_after.abs());
        if restitution_residual.abs() > velocity_tolerance {
            return Err(RigidPistonImpactError::Restitution {
                residual_m_s: restitution_residual,
                tolerance_m_s: velocity_tolerance,
            });
        }
        let wall_change = Dd::from_f64(self.mass_kg())
            * (Dd::from_f64(piston.velocity_m_s[1]) - Dd::from_f64(piston.velocity_m_s[0]))
            * (Dd::from_f64(0.5 * piston.velocity_m_s[0])
                + Dd::from_f64(0.5 * piston.velocity_m_s[1]));
        let kinetic = wall_change + rigid_kinetic_change(&rigid_jump);
        let kinetic_value = checked(kinetic)?;
        let energy_residual = checked(kinetic + Dd::from_f64(law.dissipated_energy))?;
        let energy_tolerance = checked(
            Dd::from_f64(controls.piston.energy_absolute_tolerance_j)
                + Dd::from_f64(controls.piston.energy_relative_tolerance)
                    * Dd::from_f64(kinetic_value.abs().max(law.dissipated_energy)),
        )?;
        if energy_residual.abs() > energy_tolerance {
            return Err(PistonImpulseError::Energy {
                residual_j: energy_residual,
                tolerance_j: energy_tolerance,
            }
            .into());
        }
        poll(&mut cancelled)?;
        if law.impacted {
            self.publish_impulse(&piston);
            *rigid = rigid_jump.state_after;
        }
        Ok(RigidPistonImpactReceipt {
            contact: law,
            piston,
            rigid: rigid_jump,
            piston_face_world_m: face,
            normal_world: normal,
            gap_m: gap,
            relative_normal_speed_before_m_s: before_value,
            relative_normal_speed_after_m_s: after_value,
            inverse_effective_mass_kg_inv: inverse_mass,
            kinetic_energy_change_j: kinetic_value,
            energy_residual_j: energy_residual,
            restitution_residual_m_s: restitution_residual,
            maximum_relative_momentum_residual: worst,
        })
    }
}

fn components(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn dot(a: Vec3, b: Vec3) -> Dd {
    components(a)
        .into_iter()
        .zip(components(b))
        .fold(Dd::ZERO, |sum, (a, b)| {
            sum + Dd::from_f64(a) * Dd::from_f64(b)
        })
}
fn invalid(what: &'static str) -> RigidPistonImpactError {
    PistonImpulseError::Invalid(what).into()
}
fn checked(value: Dd) -> Result<f64, RigidPistonImpactError> {
    if value.hi.is_finite() && value.lo.is_finite() && value.to_f64().is_finite() {
        Ok(value.to_f64())
    } else {
        Err(invalid("rigid/piston impulse account became nonfinite"))
    }
}
fn pairing(p0: Vec3, p1: Vec3, v0: Vec3, v1: Vec3) -> Dd {
    let [p0, p1, v0, v1] = [p0, p1, v0, v1].map(components);
    (0..3).fold(Dd::ZERO, |sum, i| {
        sum + (Dd::from_f64(p1[i]) - Dd::from_f64(p0[i]))
            * (Dd::from_f64(0.5 * v0[i]) + Dd::from_f64(0.5 * v1[i]))
    })
}
fn rigid_kinetic_change(jump: &ImpulseReceipt) -> Dd {
    let before = jump.point_kinematics_before;
    let after = jump.point_kinematics_after;
    pairing(
        jump.state_before.linear_momentum_world(),
        jump.state_after.linear_momentum_world(),
        before.center_of_mass_velocity_world,
        after.center_of_mass_velocity_world,
    ) + pairing(
        jump.state_before.angular_momentum_body(),
        jump.state_after.angular_momentum_body(),
        before.angular_velocity_body,
        after.angular_velocity_body,
    )
}
