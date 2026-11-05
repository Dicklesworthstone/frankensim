//! Instantaneous finite-patch shear and spin on original rigid-body momenta.
//!
//! The supplied normal impulse is the Hertz friction budget, not another
//! impulse applied here. The normal owner supplies that impulse and the fixed
//! pressure-impulse shape/radius. Off-centre shear and spin can change normal
//! velocity; both endpoints and the complete normal cross response are exposed.
//! This conditional phase does not impose restitution, integrate a pose, infer
//! a contact duration, discard elastic history or allocate mechanical loss as
//! heat.

use super::{
    ImpulseAccumulator, RigidImpactPoint, add_cross, add_vector, components, influence, magnitude,
    material_separation, rounded_vector,
};
use crate::{
    DynamicsError, MassProperties, MomentumImpulseReceipt, MomentumWork, PointKinematics,
    RigidBodyState, Vec3, momentum_difference,
};
use fs_math::dd::Dd;
use fs_tribo::sliding_patch::contact::impulse::{
    HertzImpulseControls, HertzImpulseError, HertzImpulseInput, HertzImpulsePatch,
    HertzImpulseReport,
};

/// The other physical patch centre and its instantaneous material motion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RigidHertzImpulseCounterpart {
    /// A second finite body. Both native bodies recoil.
    Body(RigidImpactPoint),
    /// Prescribed motion at the actual patch centre, not an arbitrary plane
    /// anchor. Angular velocity supplies the normal-axis spin work. The caller
    /// owns any relation between this local linear velocity and support rotation.
    PrescribedPatch {
        /// Original counterpart centre, which must meet the geometry gate [m].
        point_world: Vec3,
        /// Material velocity at that centre [m/s].
        velocity_world: Vec3,
        /// Instantaneous support angular velocity in world coordinates [rad/s].
        angular_velocity_world: Vec3,
    },
}

/// One caller-owned finite Hertz patch at the current original body poses.
#[derive(Clone, Copy, Debug)]
pub struct RigidHertzImpulseContact<'a> {
    /// Original COM-relative material point on the positive-impulse body.
    pub point: RigidImpactPoint,
    /// Opposite original material point or prescribed patch motion.
    pub other: RigidHertzImpulseCounterpart,
    /// Normal directed from the counterpart toward the positive-impulse body.
    pub normal_world: Vec3,
    /// Tangent gauge; a normal-parallel or nonfinite hint refuses.
    pub tangent_world: Vec3,
    /// Original radius, pressure-impulse support and one declared coefficient.
    pub law: &'a HertzImpulsePatch,
}

/// Native body acceptance gates in addition to the original impulse-law limits.
#[derive(Clone, Copy, Debug)]
pub struct RigidHertzImpulseControls {
    pub maximum_bodies: usize,
    pub geometry_tolerance_m: f64,
    /// Tangential and normal compatibility against actual native velocities.
    pub velocity_tolerance_m_s: f64,
    pub spin_tolerance_rad_s: f64,
    /// Bounds each native linear residual and the complete momentum balance.
    pub linear_momentum_tolerance_ns: f64,
    /// Bounds each native angular residual and the complete spatial balance.
    pub angular_momentum_tolerance_nms: f64,
    /// Bounds native work, total energy and constitutive-work mismatch.
    pub energy_tolerance_j: f64,
    pub contact: HertzImpulseControls,
}

/// Complete candidate mechanics. Preparation/solve alone publishes no body.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct RigidHertzImpulseReport {
    pub contact: HertzImpulseReport,
    /// [tangent one, tangent two, normal], normalized from the authored frame.
    pub axes_world: [Vec3; 3],
    /// Maps [p_x,p_y,L/a] to [dv_x,dv_y,a*domega], retaining every cross term.
    pub scaled_mobility_per_kg: [[f64; 3]; 3],
    /// Maps those same impulses to normal-velocity change. This phase does
    /// not impose any normal target; a surrounding normal solve may use it.
    pub normal_cross_mobility_per_kg: [f64; 3],
    pub slip_velocity_before_m_s: [f64; 2],
    pub slip_velocity_after_m_s: [f64; 2],
    pub spin_before_rad_s: f64,
    pub spin_after_rad_s: f64,
    pub normal_velocity_before_m_s: f64,
    pub normal_velocity_after_m_s: f64,
    pub force_impulse_world_ns: Vec3,
    pub spin_impulse_world_nms: Vec3,
    /// One original resultant update combines the point moment and spin couple.
    pub point_impulse: MomentumImpulseReceipt,
    /// Present exactly when the counterpart is another finite body.
    pub counterpart_impulse: Option<MomentumImpulseReceipt>,
    /// Absolute kinetic diagnostics for the touched finite bodies only [J].
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    /// Polarized original momentum increments, reduced across both bodies.
    pub kinetic_change_j: f64,
    /// Signed force/couple work supplied by the prescribed patch, or zero.
    pub support_work_j: f64,
    pub relative_impulse_work_j: f64,
    /// Actual native relative midpoint work, with no thermal partition [J].
    pub dissipated_energy_j: f64,
    pub work_energy_residual_j: f64,
    pub total_energy_residual_j: f64,
    pub contact_work_residual_j: f64,
    /// Actual world impulse supplied by a prescribed patch, otherwise zero.
    pub support_impulse_ns: Vec3,
    /// Its spatial moment about the world origin at the declared support point.
    pub support_angular_impulse_nms: Vec3,
    pub linear_momentum_residual_ns: Vec3,
    pub angular_momentum_residual_nms: Vec3,
    pub maximum_native_linear_residual_ns: f64,
    pub maximum_native_angular_residual_nms: f64,
    pub maximum_native_work_residual_j: f64,
    pub endpoint_tangent_residual_m_s: f64,
    pub endpoint_spin_residual_rad_s: f64,
    pub normal_response_residual_m_s: f64,
}

/// Every refusal preserves the caller's body states and reusable preparation.
#[derive(Debug)]
pub enum RigidHertzImpulseError {
    Invalid(&'static str),
    Dynamics(DynamicsError),
    Contact(HertzImpulseError),
    Cancelled,
    Audit(Box<RigidHertzImpulseReport>),
}

impl From<DynamicsError> for RigidHertzImpulseError {
    fn from(error: DynamicsError) -> Self {
        Self::Dynamics(error)
    }
}
impl From<HertzImpulseError> for RigidHertzImpulseError {
    fn from(error: HertzImpulseError) -> Self {
        match error {
            HertzImpulseError::Cancelled => Self::Cancelled,
            _ => Self::Contact(error),
        }
    }
}
impl core::fmt::Display for RigidHertzImpulseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native Hertz patch impulse: {self:?}")
    }
}
impl std::error::Error for RigidHertzImpulseError {}

#[derive(Clone, Copy, Debug)]
struct NativePoint {
    state: RigidBodyState,
    properties: MassProperties,
    arm_body: Vec3,
    kinematics: PointKinematics,
}

#[derive(Clone, Copy, Debug)]
struct Rates {
    relative_velocity: Vec3,
    relative_angular_velocity: Vec3,
    tangent: [f64; 2],
    spin: f64,
    normal: f64,
}

/// A fixed two-body/one-support snapshot with the original mobility cached.
///
/// Repeated normal-impulse or friction-coefficient trials start from this same
/// native state; the radius and physical geometry cannot change. Successful
/// preparation and trials use fixed stack storage and no heap allocation.
#[derive(Clone, Debug)]
pub struct PreparedRigidHertzImpulse {
    point: RigidImpactPoint,
    other: RigidHertzImpulseCounterpart,
    first: NativePoint,
    second: Option<NativePoint>,
    axes: [Vec3; 3],
    radius_m: f64,
    mobility: [[f64; 3]; 3],
    normal_cross: [f64; 3],
    before: Rates,
    controls: RigidHertzImpulseControls,
}

impl PreparedRigidHertzImpulse {
    /// Admit the original witnesses and cache full native shear/spin response.
    pub fn new(
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contact: RigidHertzImpulseContact<'_>,
        controls: RigidHertzImpulseControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, RigidHertzImpulseError> {
        poll(&mut cancelled)?;
        controls.contact.validate()?;
        if bodies.len() != properties.len()
            || controls.maximum_bodies == 0
            || bodies.len() > controls.maximum_bodies
            || contact.point.body >= bodies.len()
            || [
                controls.geometry_tolerance_m,
                controls.velocity_tolerance_m_s,
                controls.spin_tolerance_rad_s,
                controls.linear_momentum_tolerance_ns,
                controls.angular_momentum_tolerance_nms,
                controls.energy_tolerance_j,
            ]
            .into_iter()
            .any(|x| !x.is_finite() || x < 0.0)
        {
            return Err(RigidHertzImpulseError::Invalid(
                "matching bounded native bodies and nonnegative finite gates required",
            ));
        }
        let normal = contact.normal_world.normalized("Hertz normal")?;
        let hint = contact.tangent_world.normalized("Hertz tangent hint")?;
        let tangent = normal
            .cross(hint.cross(normal))
            .normalized("Hertz tangent projection")?;
        let axes = [
            tangent,
            normal.cross(tangent).normalized("Hertz second tangent")?,
            normal,
        ];
        let first = native_point(bodies, properties, contact.point)?;
        poll(&mut cancelled)?;
        let second = match contact.other {
            RigidHertzImpulseCounterpart::Body(point) => {
                if point.body >= bodies.len() || point.body == contact.point.body {
                    return Err(RigidHertzImpulseError::Invalid(
                        "counterpart must name a distinct finite body",
                    ));
                }
                Some(native_point(bodies, properties, point)?)
            }
            RigidHertzImpulseCounterpart::PrescribedPatch {
                point_world,
                velocity_world,
                angular_velocity_world,
            } => {
                if !point_world.is_finite()
                    || !velocity_world.is_finite()
                    || !angular_velocity_world.is_finite()
                {
                    return Err(RigidHertzImpulseError::Invalid(
                        "finite prescribed patch position and motion required",
                    ));
                }
                None
            }
        };
        let (other_centre, other_arm) = match second {
            Some(point) => (
                point.state.pose().position_world(),
                point.kinematics.arm_world,
            ),
            None => match contact.other {
                RigidHertzImpulseCounterpart::PrescribedPatch { point_world, .. } => {
                    (point_world, Vec3::ZERO)
                }
                RigidHertzImpulseCounterpart::Body(_) => unreachable!(),
            },
        };
        let separation = material_separation(
            first.state.pose().position_world(),
            first.kinematics.arm_world,
            other_centre,
            other_arm,
        );
        if !separation.is_finite() || magnitude(separation) > controls.geometry_tolerance_m {
            return Err(RigidHertzImpulseError::Invalid(
                "original patch centres do not coincide within the geometry gate",
            ));
        }
        let radius_m = contact.law.radius_m();
        let mut mobility = [[Dd::ZERO; 3]; 3];
        let mut normal_cross = [Dd::ZERO; 3];
        for point in [
            Some(contact.point),
            match contact.other {
                RigidHertzImpulseCounterpart::Body(point) => Some(point),
                RigidHertzImpulseCounterpart::PrescribedPatch { .. } => None,
            },
        ]
        .into_iter()
        .flatten()
        {
            poll(&mut cancelled)?;
            // The original force rows retain full inertia and original arms.
            // Replace only the third point-force row with the radius-scaled
            // normal-axis free couple. Its removed row supplies normal coupling.
            let mut factor = influence(bodies, properties, point, axes).factor;
            let normal_factor = factor[2];
            let rotation = properties[point.body].angular_mobility_factor(
                bodies[point.body]
                    .pose()
                    .orientation()
                    .rotate_world_to_body(normal.scale(radius_m)),
            );
            factor[2] = [0.0, 0.0, 0.0, rotation.x, rotation.y, rotation.z];
            for i in 0..3 {
                for j in 0..=i {
                    mobility[i][j] = mobility[i][j] + dot6(factor[i], factor[j]);
                    mobility[j][i] = mobility[i][j];
                }
                normal_cross[i] = normal_cross[i] + dot6(normal_factor, factor[i]);
            }
        }
        let mobility = mobility.map(|row| row.map(Dd::to_f64));
        let normal_cross = normal_cross.map(Dd::to_f64);
        if mobility
            .iter()
            .flatten()
            .chain(&normal_cross)
            .any(|x| !x.is_finite())
        {
            return Err(RigidHertzImpulseError::Invalid(
                "unrepresentable original Hertz mobility",
            ));
        }
        let before = rates(
            first.kinematics,
            second.map(|point| point.kinematics),
            contact.other,
            axes,
        )?;
        poll(&mut cancelled)?;
        Ok(Self {
            point: contact.point,
            other: contact.other,
            first,
            second,
            axes,
            radius_m,
            mobility,
            normal_cross,
            before,
            controls,
        })
    }

    #[must_use]
    pub const fn point(&self) -> RigidImpactPoint {
        self.point
    }
    #[must_use]
    pub const fn counterpart(&self) -> RigidHertzImpulseCounterpart {
        self.other
    }
    #[must_use]
    pub const fn scaled_mobility_per_kg(&self) -> &[[f64; 3]; 3] {
        &self.mobility
    }
    #[must_use]
    pub const fn normal_cross_mobility_per_kg(&self) -> &[f64; 3] {
        &self.normal_cross
    }

    /// Solve one candidate normal-budget/coefficient trial on the fixed snapshot.
    ///
    /// Normal impulse is only the original constitutive budget. This method
    /// applies no normal impulse and leaves every caller-owned state unchanged.
    /// The returned native receipts contain the complete candidate endpoints.
    #[allow(clippy::too_many_lines)]
    pub fn solve(
        &self,
        law: &HertzImpulsePatch,
        normal_impulse_ns: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RigidHertzImpulseReport, RigidHertzImpulseError> {
        poll(&mut cancelled)?;
        if law.radius_m().to_bits() != self.radius_m.to_bits() {
            return Err(RigidHertzImpulseError::Invalid(
                "prepared Hertz radius must remain unchanged",
            ));
        }
        let contact = law.solve(
            &HertzImpulseInput {
                normal_impulse_ns,
                velocity_before_m_s: self.before.tangent,
                spin_before_rad_s: self.before.spin,
                scaled_mobility_per_kg: self.mobility,
            },
            self.controls.contact,
            &mut cancelled,
        )?;
        let impulse = self.axes[0]
            .scale(contact.impulse_ns[0])
            .add(self.axes[1].scale(contact.impulse_ns[1]));
        let couple = self.axes[2].scale(contact.spin_impulse_nms);
        if !impulse.is_finite()
            || !couple.is_finite()
            || magnitude(impulse) > self.controls.contact.maximum_impulse_ns
            || magnitude(couple) > self.controls.contact.maximum_spin_impulse_nms
        {
            return Err(HertzImpulseError::ImpulseBudget {
                impulse_ns: magnitude(impulse),
                spin_impulse_nms: magnitude(couple),
            }
            .into());
        }
        poll(&mut cancelled)?;
        let point_impulse = apply(self.first, impulse, couple)?;
        poll(&mut cancelled)?;
        let counterpart_impulse = self
            .second
            .map(|point| apply(point, impulse.scale(-1.0), couple.scale(-1.0)))
            .transpose()?;
        let first_after = point_impulse
            .state_after
            .point_kinematics(self.first.properties, self.first.arm_body)?;
        let second_after = self
            .second
            .zip(counterpart_impulse)
            .map(|(point, receipt)| {
                receipt
                    .state_after
                    .point_kinematics(point.properties, point.arm_body)
            })
            .transpose()?;
        let after = rates(first_after, second_after, self.other, self.axes)?;
        let mut account = Account::default();
        for (point, receipt) in [
            Some((self.first, point_impulse)),
            self.second.zip(counterpart_impulse),
        ]
        .into_iter()
        .flatten()
        {
            poll(&mut cancelled)?;
            account.add(point, receipt)?;
        }
        let mut relative_work = MomentumWork::default();
        relative_work.add(
            self.before.relative_velocity,
            after.relative_velocity,
            components(impulse).map(Dd::from_f64),
        );
        relative_work.add(
            self.before.relative_angular_velocity,
            after.relative_angular_velocity,
            components(couple).map(Dd::from_f64),
        );
        let relative_work = relative_work.total();
        let mut support_work = Dd::ZERO;
        let mut support_linear = [Dd::ZERO; 3];
        let mut support_angular = [Dd::ZERO; 3];
        if let RigidHertzImpulseCounterpart::PrescribedPatch {
            point_world,
            velocity_world,
            angular_velocity_world,
        } = self.other
        {
            support_work = dot(impulse, velocity_world) + dot(couple, angular_velocity_world);
            add_vector(&mut support_linear, impulse);
            add_cross(&mut support_angular, point_world, impulse);
            add_vector(&mut support_angular, couple);
        }
        let kinetic = account.work.total();
        let normal_change = self
            .normal_cross
            .into_iter()
            .zip([
                contact.impulse_ns[0],
                contact.impulse_ns[1],
                contact.spin_impulse_nms / self.radius_m,
            ])
            .fold(Dd::ZERO, |sum, (w, p)| {
                sum + Dd::from_f64(w) * Dd::from_f64(p)
            });
        let endpoint_tangent_residual_m_s = (0..2)
            .map(|i| (after.tangent[i] - contact.velocity_after_m_s[i]).abs())
            .fold(0.0_f64, f64::max);
        let endpoint_spin_residual_rad_s = (after.spin - contact.spin_after_rad_s).abs();
        let normal_response_residual_m_s =
            (Dd::from_f64(after.normal) - Dd::from_f64(self.before.normal) - normal_change)
                .to_f64()
                .abs();
        let report = RigidHertzImpulseReport {
            axes_world: self.axes,
            scaled_mobility_per_kg: self.mobility,
            normal_cross_mobility_per_kg: self.normal_cross,
            slip_velocity_before_m_s: self.before.tangent,
            slip_velocity_after_m_s: after.tangent,
            spin_before_rad_s: self.before.spin,
            spin_after_rad_s: after.spin,
            normal_velocity_before_m_s: self.before.normal,
            normal_velocity_after_m_s: after.normal,
            force_impulse_world_ns: impulse,
            spin_impulse_world_nms: couple,
            point_impulse,
            counterpart_impulse,
            kinetic_before_j: account.before.to_f64(),
            kinetic_after_j: account.after.to_f64(),
            kinetic_change_j: kinetic.to_f64(),
            support_work_j: support_work.to_f64(),
            relative_impulse_work_j: relative_work.to_f64(),
            dissipated_energy_j: -relative_work.to_f64(),
            work_energy_residual_j: (kinetic - support_work - relative_work).to_f64(),
            total_energy_residual_j: (kinetic - support_work + Dd::from_f64(contact.dissipated_j))
                .to_f64(),
            contact_work_residual_j: (relative_work + Dd::from_f64(contact.dissipated_j)).to_f64(),
            support_impulse_ns: rounded_vector(support_linear),
            support_angular_impulse_nms: rounded_vector(support_angular),
            linear_momentum_residual_ns: rounded_vector(core::array::from_fn(|i| {
                account.linear[i] - support_linear[i]
            })),
            angular_momentum_residual_nms: rounded_vector(core::array::from_fn(|i| {
                account.angular[i] - support_angular[i]
            })),
            maximum_native_linear_residual_ns: account.linear_residual,
            maximum_native_angular_residual_nms: account.angular_residual,
            maximum_native_work_residual_j: account.work_residual,
            endpoint_tangent_residual_m_s,
            endpoint_spin_residual_rad_s,
            normal_response_residual_m_s,
            contact,
        };
        let c = self.controls;
        if [
            report.kinetic_before_j,
            report.kinetic_after_j,
            report.kinetic_change_j,
            report.support_work_j,
            report.relative_impulse_work_j,
            report.dissipated_energy_j,
            report.work_energy_residual_j,
            report.total_energy_residual_j,
            report.contact_work_residual_j,
            report.maximum_native_linear_residual_ns,
            report.maximum_native_angular_residual_nms,
            report.maximum_native_work_residual_j,
            endpoint_tangent_residual_m_s,
            endpoint_spin_residual_rad_s,
            normal_response_residual_m_s,
            magnitude(report.linear_momentum_residual_ns),
            magnitude(report.angular_momentum_residual_nms),
        ]
        .into_iter()
        .any(|x| !x.is_finite())
            || report.dissipated_energy_j < 0.0
            || report.work_energy_residual_j.abs() > c.energy_tolerance_j
            || report.total_energy_residual_j.abs() > c.energy_tolerance_j
            || report.contact_work_residual_j.abs() > c.energy_tolerance_j
            || report.maximum_native_work_residual_j > c.energy_tolerance_j
            || report.maximum_native_linear_residual_ns > c.linear_momentum_tolerance_ns
            || report.maximum_native_angular_residual_nms > c.angular_momentum_tolerance_nms
            || magnitude(report.linear_momentum_residual_ns) > c.linear_momentum_tolerance_ns
            || magnitude(report.angular_momentum_residual_nms) > c.angular_momentum_tolerance_nms
            || endpoint_tangent_residual_m_s > c.velocity_tolerance_m_s
            || normal_response_residual_m_s > c.velocity_tolerance_m_s
            || endpoint_spin_residual_rad_s > c.spin_tolerance_rad_s
        {
            return Err(RigidHertzImpulseError::Audit(Box::new(report)));
        }
        poll(&mut cancelled)?;
        Ok(report)
    }
}

/// Apply one conditional Hertz tangent/spin impulse phase atomically.
///
/// Only the original one or two touched bodies are staged. Each receives one
/// combined native linear/central-angular impulse; the normal budget is not
/// applied again. The caller owns normal contact, fixed-patch applicability,
/// source coefficients and subsequent pose integration.
pub fn apply_rigid_hertz_impulse(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidHertzImpulseContact<'_>,
    normal_impulse_ns: f64,
    controls: RigidHertzImpulseControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RigidHertzImpulseReport, RigidHertzImpulseError> {
    let prepared =
        PreparedRigidHertzImpulse::new(bodies, properties, contact, controls, &mut cancelled)?;
    let report = prepared.solve(contact.law, normal_impulse_ns, &mut cancelled)?;
    poll(&mut cancelled)?;
    bodies[contact.point.body] = report.point_impulse.state_after;
    if let (RigidHertzImpulseCounterpart::Body(point), Some(receipt)) =
        (contact.other, report.counterpart_impulse)
    {
        bodies[point.body] = receipt.state_after;
    }
    Ok(report)
}

fn native_point(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    point: RigidImpactPoint,
) -> Result<NativePoint, DynamicsError> {
    Ok(NativePoint {
        state: bodies[point.body],
        properties: properties[point.body],
        arm_body: point.arm_body,
        kinematics: bodies[point.body].point_kinematics(properties[point.body], point.arm_body)?,
    })
}

fn rates(
    first: PointKinematics,
    second: Option<PointKinematics>,
    other: RigidHertzImpulseCounterpart,
    axes: [Vec3; 3],
) -> Result<Rates, RigidHertzImpulseError> {
    let (relative_velocity, other_omega) = match (second, other) {
        (Some(second), RigidHertzImpulseCounterpart::Body(_)) => (
            first.relative_velocity_world(second)?,
            second.angular_velocity_world,
        ),
        (
            None,
            RigidHertzImpulseCounterpart::PrescribedPatch {
                velocity_world,
                angular_velocity_world,
                ..
            },
        ) => (
            first.velocity_relative_to_world(velocity_world)?,
            angular_velocity_world,
        ),
        _ => unreachable!(),
    };
    let relative_angular_velocity = rounded_vector(momentum_difference(
        other_omega,
        first.angular_velocity_world,
    ));
    let tangent = [
        dot(relative_velocity, axes[0]).to_f64(),
        dot(relative_velocity, axes[1]).to_f64(),
    ];
    let spin = dot(relative_angular_velocity, axes[2]).to_f64();
    let normal = dot(relative_velocity, axes[2]).to_f64();
    if !relative_angular_velocity.is_finite()
        || [tangent[0], tangent[1], spin, normal]
            .into_iter()
            .any(|x| !x.is_finite())
    {
        return Err(RigidHertzImpulseError::Invalid(
            "unrepresentable native patch rates",
        ));
    }
    Ok(Rates {
        relative_velocity,
        relative_angular_velocity,
        tangent,
        spin,
        normal,
    })
}

fn apply(
    point: NativePoint,
    impulse: Vec3,
    couple: Vec3,
) -> Result<MomentumImpulseReceipt, DynamicsError> {
    let mut resultant = ImpulseAccumulator::default();
    resultant.add(point.state, point.arm_body, impulse);
    add_vector(
        &mut resultant.angular_body,
        point
            .state
            .pose()
            .orientation()
            .rotate_world_to_body(couple),
    );
    point.state.apply_momentum_impulse(
        point.properties,
        rounded_vector(resultant.linear_world),
        rounded_vector(resultant.angular_body),
    )
}

#[derive(Default)]
struct Account {
    before: Dd,
    after: Dd,
    work: MomentumWork,
    linear: [Dd; 3],
    angular: [Dd; 3],
    linear_residual: f64,
    angular_residual: f64,
    work_residual: f64,
}
impl Account {
    fn add(
        &mut self,
        point: NativePoint,
        receipt: MomentumImpulseReceipt,
    ) -> Result<(), DynamicsError> {
        let old = point.state;
        let new = receipt.state_after;
        let dp = momentum_difference(old.linear_momentum_world(), new.linear_momentum_world());
        let dl = momentum_difference(old.angular_momentum_body(), new.angular_momentum_body());
        self.work.add(
            old.center_of_mass_velocity_world(point.properties)?,
            new.center_of_mass_velocity_world(point.properties)?,
            dp,
        );
        self.work.add(
            point
                .properties
                .angular_velocity_body_checked(old.angular_momentum_body())?,
            point
                .properties
                .angular_velocity_body_checked(new.angular_momentum_body())?,
            dl,
        );
        self.before = self.before + Dd::from_f64(receipt.work.kinetic_energy_before_j);
        self.after = self.after + Dd::from_f64(receipt.work.kinetic_energy_after_j);
        let position = components(old.pose().position_world()).map(Dd::from_f64);
        for d in 0..3 {
            self.linear[d] = self.linear[d] + dp[d];
            let j = (d + 1) % 3;
            let k = (d + 2) % 3;
            self.angular[d] = self.angular[d] + position[j] * dp[k] - position[k] * dp[j];
        }
        for part in [
            Vec3::new(dl[0].hi, dl[1].hi, dl[2].hi),
            Vec3::new(dl[0].lo, dl[1].lo, dl[2].lo),
        ] {
            add_vector(
                &mut self.angular,
                old.pose().orientation().rotate_body_to_world(part),
            );
        }
        self.linear_residual = self
            .linear_residual
            .max(magnitude(receipt.linear_residual_world));
        self.angular_residual = self
            .angular_residual
            .max(magnitude(receipt.angular_residual_body));
        self.work_residual = self
            .work_residual
            .max(receipt.work.work_energy_residual_j.abs());
        Ok(())
    }
}

fn dot(a: Vec3, b: Vec3) -> Dd {
    components(a)
        .into_iter()
        .zip(components(b))
        .fold(Dd::ZERO, |sum, (a, b)| {
            sum + Dd::from_f64(a) * Dd::from_f64(b)
        })
}
fn dot6(a: [f64; 6], b: [f64; 6]) -> Dd {
    a.into_iter().zip(b).fold(Dd::ZERO, |sum, (a, b)| {
        sum + Dd::from_f64(a) * Dd::from_f64(b)
    })
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), RigidHertzImpulseError> {
    if cancelled() {
        Err(RigidHertzImpulseError::Cancelled)
    } else {
        Ok(())
    }
}
