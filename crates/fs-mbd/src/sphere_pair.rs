//! Complete two-sphere flight and frictional impact under common acceleration.
//!
//! Both spheres have finite mass, isotropic inertia and arbitrary spin. The
//! common constant acceleration cancels from relative motion, allowing one
//! swept-sphere query, at most one impact, and exact flight on each side.
//! Normal and Coulomb response belong to fs-contact; this adapter supplies
//! geometry, the two-body effective masses and physical momentum application.
//!
//! Impulses act equally and oppositely at the real facing surface points.
//! Translation and torque-free attitude are analytic, up to floating arithmetic.
//! Global linear/angular momentum and force-work/energy balances must close
//! before a complete result is returned. No position or momentum is repaired.
//! This is a floating analytic sphere query, not interval-certified general CCD.
//! It does not model unequal accelerations, external torques, anisotropic inertia,
//! deformation, adhesion, rolling resistance, or a many-body contact network.

use crate::{DynamicsError, MassProperties, Pose, RigidBodyState, Vec3};
use fs_contact::impact::{NormalImpactError, NormalImpactInput, solve_normal_impact};
use fs_contact::swept_sphere::{SweptSphereError, SweptSphereInput, first_sphere_contact};
use fs_contact::tangential::impact::{
    CoulombImpulseError, CoulombImpulseInput, FrictionLaw, solve_isotropic_coulomb_impulse,
};

/// A finite-mass sphere with its centre of mass at the geometric centre.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpherePairBody {
    /// Existing validated rigid-body mass and principal inertia.
    pub properties: MassProperties,
    /// Positive finite enclosing radius, in metres.
    pub radius: f64,
}

/// Two spheres and the laws governing their free flight and possible impact.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpherePairScene {
    /// Sphere zero and sphere one; their masses and radii may differ.
    pub bodies: [SpherePairBody; 2],
    /// Common constant acceleration in world coordinates, in m/s².
    pub acceleration_world: Vec3,
    /// Newton normal restitution in the closed interval [0, 1].
    pub restitution: f64,
    /// Nonnegative ideal Coulomb coefficient, equal for static and kinetic slip.
    pub friction_coefficient: f64,
}

/// The single detected contact, including zero-impulse interior grazing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpherePairImpact {
    /// Contact time measured from the start of this advance, in seconds.
    pub time_s: f64,
    /// Unit normal from sphere zero towards sphere one.
    pub normal_world: Vec3,
    /// Facing surface points of sphere zero and sphere one at the event.
    ///
    /// Their separation may contain the reported admission's floating residual;
    /// neither point is moved to force exact coincidence.
    pub points_world: [Vec3; 2],
    /// Nonnegative repulsive normal impulse magnitude, in kg m/s.
    pub normal_impulse: f64,
    /// Tangential impulse on sphere one; sphere zero receives its negative.
    pub tangent_impulse_world: Vec3,
    /// Relative centre normal velocity (one minus zero) before impact, in m/s.
    pub relative_normal_speed_before: f64,
    /// Measured relative centre normal velocity after impact, in m/s.
    pub relative_normal_speed_after: f64,
}

/// Complete accepted advance. No input state is changed on success or refusal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpherePairAdvance {
    /// Both momentum-form states at the end of the requested duration.
    pub states_after: [RigidBodyState; 2],
    /// Earliest detected contact, if one occurred within this duration.
    pub impact: Option<SpherePairImpact>,
    /// Work of both common-acceleration forces along the actual trajectories.
    pub external_force_work: f64,
    /// Mechanical energy lost to normal restitution and ideal Coulomb friction.
    ///
    /// This is not an automatic heat assignment to either sphere.
    pub dissipated_energy: f64,
    /// Final kinetic energy minus initial energy minus force work plus loss.
    pub energy_balance_residual: f64,
    /// Total final minus initial linear momentum minus applied force impulses.
    pub linear_momentum_residual: Vec3,
    /// Total orbital-plus-spin change minus applied angular impulses about origin.
    pub angular_momentum_residual: Vec3,
    /// Actual final centre distance minus the sum of the radii, in metres.
    ///
    /// A detected event may leave a bounded signed floating-point residual.
    pub gap_after: f64,
}

/// Domain, detection, response or conservation refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum SpherePairError {
    /// Invalid shape, material law, duration or nonfinite derived arithmetic.
    Invalid(&'static str),
    /// Strictly negative initial geometric gap; no repair is attempted.
    Penetration {
        /// Measured initial signed gap.
        gap: f64,
    },
    /// Existing rigid-body admission or momentum operation refused.
    Dynamics(DynamicsError),
    /// Shared scalar normal law refused.
    Contact(NormalImpactError),
    /// Shared isotropic Coulomb law refused.
    Friction(CoulombImpulseError),
    /// Shared swept-sphere detection refused.
    Detection(SweptSphereError),
    /// Cancellation was observed before complete publication.
    Cancelled,
    /// Event geometry, endpoint response or a conservation balance did not close.
    Unresolved(&'static str),
}

impl From<DynamicsError> for SpherePairError {
    fn from(error: DynamicsError) -> Self {
        Self::Dynamics(error)
    }
}
impl From<NormalImpactError> for SpherePairError {
    fn from(error: NormalImpactError) -> Self {
        Self::Contact(error)
    }
}
impl From<CoulombImpulseError> for SpherePairError {
    fn from(error: CoulombImpulseError) -> Self {
        Self::Friction(error)
    }
}
impl From<SweptSphereError> for SpherePairError {
    fn from(error: SweptSphereError) -> Self {
        Self::Detection(error)
    }
}
impl std::fmt::Display for SpherePairError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "sphere pair advance refused: {self:?}")
    }
}
impl std::error::Error for SpherePairError {}

#[derive(Clone, Copy)]
struct Mechanics {
    shapes: [SpherePairBody; 2],
    acceleration: Vec3,
    forces: [Vec3; 2],
    radii_sum: f64,
    inverse_normal_mass: f64,
    inverse_tangent_mass: f64,
}
impl Mechanics {
    fn new(scene: &SpherePairScene) -> Result<Self, SpherePairError> {
        if !(scene.restitution.is_finite() && (0.0..=1.0).contains(&scene.restitution))
            || !(scene.friction_coefficient.is_finite() && scene.friction_coefficient >= 0.0)
        {
            return Err(SpherePairError::Invalid(
                "restitution or friction coefficient",
            ));
        }
        vector(scene.acceleration_world)?;
        let mut inverse_normal_mass = 0.0;
        let mut inverse_tangent_mass = 0.0;
        let mut forces = [Vec3::ZERO; 2];
        for (i, shape) in scene.bodies.iter().enumerate() {
            shape.properties.validate()?;
            if !(shape.radius.is_finite() && shape.radius > 0.0) {
                return Err(SpherePairError::Invalid(
                    "positive finite sphere radius required",
                ));
            }
            let moments = shape.properties.principal_inertia_body();
            let bound =
                finite((shape.properties.mass() * shape.radius) * shape.radius * (2.0 / 3.0))?;
            if moments.x != moments.y || moments.x != moments.z || moments.x > bound {
                return Err(SpherePairError::Invalid(
                    "isotropic inertia must fit within each sphere",
                ));
            }
            let inverse_mass = finite(shape.properties.mass().recip())?;
            inverse_normal_mass = finite(inverse_normal_mass + inverse_mass)?;
            inverse_tangent_mass = finite(
                inverse_tangent_mass + inverse_mass + (shape.radius / moments.x) * shape.radius,
            )?;
            forces[i] = scene.acceleration_world.scale(shape.properties.mass());
            vector(forces[i])?;
        }
        if inverse_normal_mass <= 0.0 || inverse_tangent_mass <= 0.0 {
            return Err(SpherePairError::Invalid(
                "positive representable effective mass required",
            ));
        }
        Ok(Self {
            shapes: scene.bodies,
            acceleration: scene.acceleration_world,
            forces,
            radii_sum: finite(scene.bodies[0].radius + scene.bodies[1].radius)?,
            inverse_normal_mass,
            inverse_tangent_mass,
        })
    }

    fn relative_position(self, bodies: [RigidBodyState; 2]) -> Result<Vec3, SpherePairError> {
        let value = bodies[1]
            .pose()
            .position_world()
            .sub(bodies[0].pose().position_world());
        vector(value)?;
        Ok(value)
    }
    fn velocities(self, bodies: [RigidBodyState; 2]) -> Result<[Vec3; 2], SpherePairError> {
        Ok([
            bodies[0].center_of_mass_velocity_world(self.shapes[0].properties)?,
            bodies[1].center_of_mass_velocity_world(self.shapes[1].properties)?,
        ])
    }
    fn gap(self, bodies: [RigidBodyState; 2]) -> Result<f64, SpherePairError> {
        finite(norm(self.relative_position(bodies)?)? - self.radii_sum)
    }
    fn geometry_allowance(self, bodies: [RigidBodyState; 2]) -> Result<f64, SpherePairError> {
        finite(
            512.0
                * f64::EPSILON
                * (norm(bodies[0].pose().position_world())?
                    + norm(bodies[1].pose().position_world())?
                    + self.radii_sum),
        )
    }
    fn kinetic(self, bodies: [RigidBodyState; 2]) -> Result<f64, SpherePairError> {
        finite(
            bodies[0].kinetic_energy(self.shapes[0].properties)?
                + bodies[1].kinetic_energy(self.shapes[1].properties)?,
        )
    }

    fn flight(
        self,
        bodies: [RigidBodyState; 2],
        h: f64,
        ledger: &mut Ledger,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[RigidBodyState; 2], SpherePairError> {
        if h == 0.0 {
            return Ok(bodies);
        }
        let mut result = bodies;
        for i in 0..2 {
            poll(cancelled)?;
            let body = bodies[i];
            let position = body.pose().position_world();
            let velocity = body.center_of_mass_velocity_world(self.shapes[i].properties)?;
            let displacement = velocity.add(self.acceleration.scale(0.5 * h)).scale(h);
            let spin = world_spin(body);
            let inertia = self.shapes[i].properties.principal_inertia_body().x;
            let orientation = body
                .pose()
                .orientation()
                .right_exp(body.angular_momentum_body().scale(h / inertia))?;
            result[i] = RigidBodyState::new(
                Pose::new(position.add(displacement), orientation)?,
                body.linear_momentum_world().add(self.forces[i].scale(h)),
                orientation.rotate_world_to_body(spin),
            )?;
            // The h³ term vanishes: each applied force is parallel to acceleration.
            let angular = position
                .cross(self.forces[i])
                .scale(h)
                .add(velocity.cross(self.forces[i]).scale(0.5 * h).scale(h));
            vector(angular)?;
            ledger.force_work = finite(ledger.force_work + self.forces[i].dot(displacement))?;
            ledger.force_angular = ledger.force_angular.add(angular);
            ledger.angular_scale = finite(ledger.angular_scale + norm(angular)?)?;
        }
        vector(ledger.force_angular)?;
        Ok(result)
    }
}
#[derive(Default)]
struct Ledger {
    force_work: f64,
    force_angular: Vec3,
    angular_scale: f64,
    loss: f64,
    impulse_scale: f64,
}

/// Advance both finite-mass spheres through their first detected impact.
///
/// Both forces are mass times the same constant acceleration. Relative free
/// flight is therefore linear. Once the accepted impact has a nonclosing normal
/// response, that linear relative trajectory cannot return to contact, so the
/// complete requested duration requires at most two flights and one impact.
/// No timestep, contact iteration or attitude microstep is hidden in this API.
///
/// Newly supplied negative gaps refuse strictly. At a detected event the actual
/// geometry and response must agree within floating-point roundoff; positions
/// are never snapped to contact. An unresolved endpoint or global balance
/// refuses atomically. Cancellation is polled between the bounded phases.
pub fn advance_sphere_pair(
    scene: &SpherePairScene,
    states: [RigidBodyState; 2],
    duration_s: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SpherePairAdvance, SpherePairError> {
    poll(&mut cancelled)?;
    if !(duration_s.is_finite() && duration_s > 0.0) {
        return Err(SpherePairError::Invalid(
            "positive finite duration required",
        ));
    }
    let mechanics = Mechanics::new(scene)?;
    for state in states {
        state.validate()?;
    }
    let initial_gap = mechanics.gap(states)?;
    if initial_gap < 0.0 {
        return Err(SpherePairError::Penetration { gap: initial_gap });
    }
    let initial_energy = mechanics.kinetic(states)?;
    let velocities = mechanics.velocities(states)?;
    let relative_velocity = velocities[1].sub(velocities[0]);
    vector(relative_velocity)?;
    let detection = first_sphere_contact(SweptSphereInput {
        relative_position_world: components(mechanics.relative_position(states)?),
        relative_velocity_world: components(relative_velocity),
        radius_sum: mechanics.radii_sum,
        duration_s,
    })?;
    poll(&mut cancelled)?;
    let mut ledger = Ledger::default();
    let (states_after, impact) = if let Some(hit) = detection {
        if !(hit.time_s.is_finite() && hit.time_s >= 0.0 && hit.time_s <= duration_s) {
            return Err(SpherePairError::Unresolved(
                "contact time outside the requested interval",
            ));
        }
        let before = mechanics.flight(states, hit.time_s, &mut ledger, &mut cancelled)?;
        let gap = mechanics.gap(before)?;
        if gap.abs() > mechanics.geometry_allowance(before)? {
            return Err(SpherePairError::Unresolved("detected event geometry"));
        }
        let normal = mechanics
            .relative_position(before)?
            .normalized("sphere contact normal")?;
        let predicted_normal = Vec3::new(
            hit.normal_world[0],
            hit.normal_world[1],
            hit.normal_world[2],
        );
        if norm(normal.sub(predicted_normal))? > 4096.0 * f64::EPSILON {
            return Err(SpherePairError::Unresolved("detected event normal"));
        }
        let axis = if normal.x.abs() < 0.8 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        let tangent0 = normal.cross(axis).normalized("sphere contact tangent")?;
        let tangent1 = normal
            .cross(tangent0)
            .normalized("sphere contact bitangent")?;
        let arms = [
            normal.scale(mechanics.shapes[0].radius),
            normal.scale(-mechanics.shapes[1].radius),
        ];
        let arms_body = [
            before[0].pose().orientation().rotate_world_to_body(arms[0]),
            before[1].pose().orientation().rotate_world_to_body(arms[1]),
        ];
        let points_world = [
            before[0].pose().position_world().add(arms[0]),
            before[1].pose().position_world().add(arms[1]),
        ];
        vector(points_world[0])?;
        vector(points_world[1])?;
        let point_before = [
            before[0].point_kinematics(mechanics.shapes[0].properties, arms_body[0])?,
            before[1].point_kinematics(mechanics.shapes[1].properties, arms_body[1])?,
        ];
        let center_relative = point_before[1]
            .center_of_mass_velocity_world
            .sub(point_before[0].center_of_mass_velocity_world);
        let relative_before = center_relative
            .add(
                point_before[1]
                    .angular_velocity_world
                    .cross(point_before[1].arm_world),
            )
            .sub(
                point_before[0]
                    .angular_velocity_world
                    .cross(point_before[0].arm_world),
            );
        vector(relative_before)?;
        let speed_before = finite(center_relative.dot(normal))?;
        let mut velocity_scale = norm(center_relative)?.max(norm(relative_before)?);
        for point in point_before {
            velocity_scale = velocity_scale
                .max(norm(point.center_of_mass_velocity_world)?)
                .max(norm(point.angular_velocity_world)? * norm(point.arm_world)?);
        }
        let velocity_allowance = finite(8192.0 * f64::EPSILON * velocity_scale)?;
        if (speed_before - hit.relative_normal_speed).abs() > velocity_allowance {
            return Err(SpherePairError::Unresolved("detected event velocity"));
        }
        let normal_result = solve_normal_impact(NormalImpactInput {
            gap: 0.0,
            activation_tolerance: 0.0,
            relative_normal_speed: speed_before,
            inverse_effective_mass: mechanics.inverse_normal_mass,
            restitution: scene.restitution,
        })?;
        let friction = solve_isotropic_coulomb_impulse(
            CoulombImpulseInput {
                normal_impulse: normal_result.normal_impulse,
                relative_tangent_velocity: [
                    relative_before.dot(tangent0),
                    relative_before.dot(tangent1),
                ],
                inverse_tangent_mass: mechanics.inverse_tangent_mass,
            },
            &FrictionLaw::Coulomb {
                static_mu: scene.friction_coefficient,
                kinetic_mu: scene.friction_coefficient,
            },
        )?;
        let tangent = tangent0
            .scale(friction.tangent_impulse[0])
            .add(tangent1.scale(friction.tangent_impulse[1]));
        let impulse = normal.scale(normal_result.normal_impulse).add(tangent);
        vector(impulse)?;
        poll(&mut cancelled)?;
        let applied0 = before[0].apply_impulse_at_body_point(
            mechanics.shapes[0].properties,
            arms_body[0],
            impulse.scale(-1.0),
        )?;
        let applied1 = before[1].apply_impulse_at_body_point(
            mechanics.shapes[1].properties,
            arms_body[1],
            impulse,
        )?;
        let after = [applied0.state_after, applied1.state_after];
        let center_after = mechanics.velocities(after)?;
        let speed_after = finite(center_after[1].sub(center_after[0]).dot(normal))?;
        let relative_after = center_after[1]
            .sub(center_after[0])
            .add(
                applied1
                    .point_kinematics_after
                    .angular_velocity_world
                    .cross(applied1.point_kinematics_after.arm_world),
            )
            .sub(
                applied0
                    .point_kinematics_after
                    .angular_velocity_world
                    .cross(applied0.point_kinematics_after.arm_world),
            );
        vector(relative_after)?;
        if speed_after < -velocity_allowance
            || (speed_after - normal_result.relative_normal_speed_after).abs() > velocity_allowance
            || (relative_after.dot(tangent0) - friction.relative_tangent_velocity_after[0]).abs()
                > velocity_allowance
            || (relative_after.dot(tangent1) - friction.relative_tangent_velocity_after[1]).abs()
                > velocity_allowance
        {
            return Err(SpherePairError::Unresolved(
                "physical impact endpoint velocity",
            ));
        }
        ledger.loss = finite(normal_result.dissipated_energy + friction.dissipated_energy)?;
        ledger.impulse_scale = norm(impulse)?;
        let relative_work = finite(impulse.dot(relative_before.add(relative_after).scale(0.5)))?;
        let work_scale = finite(
            norm(impulse)? * (norm(relative_before)? + norm(relative_after)?) + ledger.loss,
        )?;
        if (relative_work + ledger.loss).abs() > finite(8192.0 * f64::EPSILON * work_scale)? {
            return Err(SpherePairError::Unresolved("physical impact work balance"));
        }
        let impact = SpherePairImpact {
            time_s: hit.time_s,
            normal_world: normal,
            points_world,
            normal_impulse: normal_result.normal_impulse,
            tangent_impulse_world: tangent,
            relative_normal_speed_before: speed_before,
            relative_normal_speed_after: speed_after,
        };
        (
            mechanics.flight(after, duration_s - hit.time_s, &mut ledger, &mut cancelled)?,
            Some(impact),
        )
    } else {
        (
            mechanics.flight(states, duration_s, &mut ledger, &mut cancelled)?,
            None,
        )
    };
    poll(&mut cancelled)?;
    let gap_after = mechanics.gap(states_after)?;
    if gap_after < -mechanics.geometry_allowance(states_after)? {
        return Err(SpherePairError::Unresolved("penetrating final geometry"));
    }
    let final_energy = mechanics.kinetic(states_after)?;
    let energy_balance_residual =
        finite((final_energy - initial_energy) - ledger.force_work + ledger.loss)?;
    let external_impulse = mechanics.forces[0]
        .scale(duration_s)
        .add(mechanics.forces[1].scale(duration_s));
    let linear_momentum_residual = total_linear(states_after)
        .sub(total_linear(states))
        .sub(external_impulse);
    let angular_momentum_residual = total_angular(states_after)
        .sub(total_angular(states))
        .sub(ledger.force_angular);
    vector(linear_momentum_residual)?;
    vector(angular_momentum_residual)?;
    let energy_scale =
        finite(initial_energy.abs() + final_energy.abs() + ledger.force_work.abs() + ledger.loss)?;
    let mut linear_scale = norm(external_impulse)? + 2.0 * ledger.impulse_scale;
    let mut angular_scale = ledger.angular_scale;
    for i in 0..2 {
        linear_scale += norm(states[i].linear_momentum_world())?
            + norm(states_after[i].linear_momentum_world())?;
        angular_scale += angular_magnitude(states[i])? + angular_magnitude(states_after[i])?;
    }
    let allowance = 16384.0 * f64::EPSILON;
    if energy_balance_residual.abs() > finite(allowance * energy_scale)?
        || norm(linear_momentum_residual)? > finite(allowance * linear_scale)?
        || norm(angular_momentum_residual)? > finite(allowance * angular_scale)?
    {
        return Err(SpherePairError::Unresolved(
            "global momentum or energy balance",
        ));
    }
    Ok(SpherePairAdvance {
        states_after,
        impact,
        external_force_work: ledger.force_work,
        dissipated_energy: ledger.loss,
        energy_balance_residual,
        linear_momentum_residual,
        angular_momentum_residual,
        gap_after,
    })
}

fn world_spin(body: RigidBodyState) -> Vec3 {
    body.pose()
        .orientation()
        .rotate_body_to_world(body.angular_momentum_body())
}
fn total_linear(bodies: [RigidBodyState; 2]) -> Vec3 {
    bodies[0]
        .linear_momentum_world()
        .add(bodies[1].linear_momentum_world())
}
fn total_angular(bodies: [RigidBodyState; 2]) -> Vec3 {
    bodies[0]
        .pose()
        .position_world()
        .cross(bodies[0].linear_momentum_world())
        .add(world_spin(bodies[0]))
        .add(
            bodies[1]
                .pose()
                .position_world()
                .cross(bodies[1].linear_momentum_world()),
        )
        .add(world_spin(bodies[1]))
}
fn angular_magnitude(body: RigidBodyState) -> Result<f64, SpherePairError> {
    finite(
        norm(
            body.pose()
                .position_world()
                .cross(body.linear_momentum_world()),
        )? + norm(world_spin(body))?,
    )
}
fn components(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}
fn norm(value: Vec3) -> Result<f64, SpherePairError> {
    Ok(value.stable_norm("sphere vector magnitude")?)
}
fn vector(value: Vec3) -> Result<(), SpherePairError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(SpherePairError::Invalid("nonfinite sphere arithmetic"))
    }
}
fn finite(value: f64) -> Result<f64, SpherePairError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(SpherePairError::Invalid("nonfinite sphere arithmetic"))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), SpherePairError> {
    if cancelled() {
        Err(SpherePairError::Cancelled)
    } else {
        Ok(())
    }
}
