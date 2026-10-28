//! Event-resolved isotropic sphere motion against one prescribed plane.
//!
//! The plane normal, velocity and applied centre force are constant. Ideal
//! Coulomb friction has equal static/kinetic coefficients. Persistent contact
//! admits driven rolling and collinear sliding with exact slip-arrest events.
//! Changing-direction sliding uses explicitly bounded implicit Coulomb steps,
//! with first-order continuum convergence and timestep-resolved slip arrest.
//! Ballistic contact timing and both impulse laws remain in fs-contact.
//! Centre trajectories and world momenta use analytic constant-force segments;
//! their curved-sliding force is a time-discrete approximation of Coulomb flow.
//! Attitude is exact when spin and torque commute, including torque-free
//! flight and sliding from rest. Otherwise bounded second-order Lie midpoint
//! steps resolve the changing spin axis. There is no penetration projection.

use crate::{DynamicsError, MassProperties, Pose, RigidBodyState, Vec3};
use core::fmt;
use fs_contact::ballistic::{
    BallisticContactError, BallisticContactInput, first_ballistic_contact,
};
use fs_contact::impact::{NormalImpactError, NormalImpactInput, solve_normal_impact};
use fs_contact::tangential::impact::{
    CoulombImpulseError, CoulombImpulseInput, CoulombImpulseResult, FrictionLaw,
    solve_isotropic_coulomb_impulse,
};

/// Physical parameters of the bounded sphere/plane trajectory.
#[derive(Clone, Copy, Debug)]
pub struct SpherePlane {
    /// Positive mass and equal principal moments about the sphere centre.
    pub properties: MassProperties,
    /// Bounding sphere radius in metres; I <= 2 m r^2 / 3 is required.
    pub radius: f64,
    /// A world point on the plane at absolute time zero.
    pub point_world: Vec3,
    /// Nonzero fixed outward normal, normalized internally.
    pub normal_world: Vec3,
    /// Constant prescribed material velocity of the plane, in m/s.
    pub velocity_world: Vec3,
    /// Constant applied centre force in newtons. Resolved transverse driving
    /// uses the declared curved-sliding timestep and work bounds.
    pub force_world: Vec3,
    /// Newton normal restitution in [0,1].
    pub restitution: f64,
    /// Common nonnegative static and kinetic Coulomb coefficient.
    pub friction_coefficient: f64,
}

/// Contact mode; rolling means zero relative point slip to numerical resolution.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpherePlaneMode {
    /// Unconstrained motion, including separation after impact.
    Flight,
    /// Supported contact with tangential slip.
    Sliding,
    /// Supported contact without resolved tangential slip.
    Rolling,
}

/// Momentum state and private event-continuation information.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpherePlaneState {
    body: RigidBodyState,
    time_s: f64,
    mode: SpherePlaneMode,
    contact_known: bool,
    zero_normal: bool,
    velocity_scale: f64,
}
impl SpherePlaneState {
    /// Start a new trajectory. Geometry is checked strictly by advance:
    /// even a tiny negative initial gap refuses.
    pub fn new(body: RigidBodyState, time_s: f64) -> Result<Self, SpherePlaneError> {
        finite(time_s)?;
        Ok(Self {
            body,
            time_s,
            mode: SpherePlaneMode::Flight,
            contact_known: false,
            zero_normal: false,
            velocity_scale: 0.0,
        })
    }
    /// Physical state, using the existing rigid-body momentum owner.
    pub const fn body(self) -> RigidBodyState {
        self.body
    }
    /// Absolute time in seconds.
    pub const fn time_s(self) -> f64 {
        self.time_s
    }
    /// Last resolved contact mode.
    pub const fn mode(self) -> SpherePlaneMode {
        self.mode
    }
}

/// Event, curved-sliding and attitude accuracy/work limits.
#[derive(Clone, Copy, Debug)]
pub struct SpherePlaneConfig {
    /// Maximum published impacts and contact-mode events, at most 4096.
    pub max_events: usize,
    /// Maximum Lie-midpoint substep for noncommuting spin and contact torque.
    /// Commuting spin uses one exact exponential. Translation and momenta
    /// do not use these substeps.
    pub max_attitude_step_s: f64,
    /// Maximum attitude substeps per advance, at most one million.
    pub max_attitude_substeps: usize,
    /// Maximum step for changing-direction Coulomb sliding. This is a
    /// first-order velocity-implicit discretization, separate from the exact
    /// collinear segments and their exact slip-arrest times. Refine this value
    /// to resolve the curved trajectory and its timestep-boundary arrest.
    pub max_sliding_step_s: f64,
    /// Maximum curved-sliding trials per advance, at most one million.
    /// Rejected direction-reversal trials count toward this work bound.
    pub max_sliding_substeps: usize,
}
impl Default for SpherePlaneConfig {
    fn default() -> Self {
        Self {
            max_events: 128,
            max_attitude_step_s: 0.01,
            max_attitude_substeps: 65_536,
            max_sliding_step_s: 0.01,
            max_sliding_substeps: 65_536,
        }
    }
}

/// Discrete events in chronological order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpherePlaneEventKind {
    /// Normal/tangential instantaneous contact response.
    Impact,
    /// Entry into persistent normal support.
    Support,
    /// Sliding has reached zero relative point slip; driving force may then
    /// start opposite slip if the static traction limit is exceeded.
    Rolling,
    /// Applied tangent force exceeds static traction and starts sliding.
    Sliding,
    /// A supported sphere departs from the plane.
    LiftOff,
}

/// One event; continuous support is accounted in the aggregate impulses.
#[derive(Clone, Copy, Debug)]
pub struct SpherePlaneEvent {
    /// Absolute event time, in seconds.
    pub time_s: f64,
    /// Event kind.
    pub kind: SpherePlaneEventKind,
    /// Normal impact impulse; zero for mode-only events.
    pub normal_impulse: f64,
    /// Tangential impact impulse; zero for mode-only events.
    pub tangent_impulse_world: Vec3,
}

/// Complete advance, published only after momentum and work admission.
#[derive(Clone, Debug)]
pub struct SpherePlaneAdvance {
    /// State at exactly the requested end time.
    pub state_after: SpherePlaneState,
    /// Bounded chronological event sequence.
    pub events: Vec<SpherePlaneEvent>,
    /// Total impact plus sustained-contact impulse, in N s.
    pub contact_impulse_world: Vec3,
    /// Total contact angular impulse about the world origin, in N m s.
    pub contact_angular_impulse_world: Vec3,
    /// Work of the declared constant applied force, in J.
    pub external_force_work: f64,
    /// Signed work supplied by the moving plane, in J.
    pub plane_work: f64,
    /// Normal impact and actual tangent-force work loss, in J; not a heat
    /// allocation. Curved-sliding loss belongs to the declared time-discrete
    /// trajectory and converges with its sliding timestep.
    pub dissipated_energy: f64,
    /// Actual kinetic-energy change minus supplied work plus loss, in J.
    pub energy_balance_residual: f64,
    /// Actual momentum change minus applied force and contact impulses.
    pub linear_momentum_residual: Vec3,
    /// Actual orbital-plus-spin change minus all world-origin torque integrals.
    pub angular_momentum_residual: Vec3,
    /// Measured final geometric gap. Known contact may have only the recorded
    /// floating-point geometry residual; position is never projected.
    pub gap_after: f64,
    /// Accepted plus rejected curved-sliding trials; zero for analytic paths.
    pub sliding_substeps: usize,
}

/// Atomic refusal; no partial trajectory or repaired state is returned.
#[derive(Debug)]
pub enum SpherePlaneError {
    /// Unsupported physical parameters or invalid work limits.
    Invalid(&'static str),
    /// Newly supplied or unresolved penetrating geometry.
    Penetration {
        /// Signed sphere/plane gap in metres.
        gap: f64,
    },
    /// Existing rigid-body owner refused.
    Dynamics(DynamicsError),
    /// Existing normal impulse owner refused.
    Contact(NormalImpactError),
    /// Existing Coulomb owner refused.
    Friction(CoulombImpulseError),
    /// Ballistic event timing refused.
    Detection(BallisticContactError),
    /// Impact/mode event cap, including unresolved accumulating bounces.
    EventBudget,
    /// Attitude accuracy request exceeds its explicit work cap.
    AttitudeBudget,
    /// Curved-sliding accuracy request exceeds its explicit work cap.
    SlidingBudget,
    /// A positive event interval cannot advance the floating-point clock.
    NoProgress,
    /// Cancellation before publication.
    Cancelled,
    /// Physical geometry, velocity, momentum, or work did not resolve.
    Unresolved(&'static str),
}
impl fmt::Display for SpherePlaneError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "sphere/plane trajectory: {self:?}")
    }
}
impl std::error::Error for SpherePlaneError {}
impl From<DynamicsError> for SpherePlaneError {
    fn from(e: DynamicsError) -> Self {
        Self::Dynamics(e)
    }
}
impl From<NormalImpactError> for SpherePlaneError {
    fn from(e: NormalImpactError) -> Self {
        Self::Contact(e)
    }
}
impl From<CoulombImpulseError> for SpherePlaneError {
    fn from(e: CoulombImpulseError) -> Self {
        Self::Friction(e)
    }
}
impl From<BallisticContactError> for SpherePlaneError {
    fn from(e: BallisticContactError) -> Self {
        Self::Detection(e)
    }
}

#[derive(Default)]
struct Ledger {
    contact: Vec3,
    contact_angular: Vec3,
    external_angular: Vec3,
    force_work: f64,
    plane_work: f64,
    loss: f64,
    impulse_scale: f64,
    angular_scale: f64,
}
struct Engine {
    scene: SpherePlane,
    config: SpherePlaneConfig,
    n: Vec3,
    tangent: [Vec3; 2],
    mass: f64,
    inertia: f64,
    inverse_tangent_mass: f64,
    acceleration_normal: f64,
    force_tangent_world: Vec3,
    acceleration_tangent: [f64; 2],
    arm: Vec3,
    ledger: Ledger,
    events: Vec<SpherePlaneEvent>,
    attitude_steps: usize,
    sliding_steps: usize,
    velocity_scale: f64,
}
impl Engine {
    fn new(scene: SpherePlane, config: SpherePlaneConfig) -> Result<Self, SpherePlaneError> {
        scene.properties.validate()?;
        if !(scene.radius.is_finite() && scene.radius > 0.0)
            || !(scene.restitution.is_finite() && (0.0..=1.0).contains(&scene.restitution))
            || !(scene.friction_coefficient.is_finite() && scene.friction_coefficient >= 0.0)
            || config.max_events == 0
            || config.max_events > 4096
            || !(config.max_attitude_step_s.is_finite() && config.max_attitude_step_s > 0.0)
            || config.max_attitude_substeps == 0
            || config.max_attitude_substeps > 1_000_000
            || !(config.max_sliding_step_s.is_finite() && config.max_sliding_step_s > 0.0)
            || config.max_sliding_substeps == 0
            || config.max_sliding_substeps > 1_000_000
        {
            return Err(SpherePlaneError::Invalid(
                "sphere, restitution, friction, or work limits",
            ));
        }
        vector(scene.point_world)?;
        vector(scene.velocity_world)?;
        vector(scene.force_world)?;
        let n = scene.normal_world.normalized("plane normal")?;
        let mass = scene.properties.mass();
        let moments = scene.properties.principal_inertia_body();
        let inertia_bound = finite((mass * scene.radius) * scene.radius * (2.0 / 3.0))?;
        if moments.x != moments.y || moments.x != moments.z || moments.x > inertia_bound {
            return Err(SpherePlaneError::Invalid(
                "isotropic inertia must fit within the sphere",
            ));
        }
        let force_n = finite(scene.force_world.dot(n))?;
        let force_tangent_world = scene.force_world.sub(n.scale(force_n));
        vector(force_tangent_world)?;
        let axis = if n.x.abs() < 0.8 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        let t0 = n.cross(axis).normalized("plane tangent")?;
        let t1 = n.cross(t0).normalized("plane bitangent")?;
        let inverse_tangent_mass =
            finite(mass.recip() + (scene.radius / moments.x) * scene.radius)?;
        if inverse_tangent_mass <= 0.0 {
            return Err(SpherePlaneError::Invalid("tangent inverse mass"));
        }
        Ok(Self {
            scene,
            config,
            n,
            tangent: [t0, t1],
            mass,
            inertia: moments.x,
            inverse_tangent_mass,
            acceleration_normal: finite(force_n / mass)?,
            force_tangent_world,
            acceleration_tangent: [
                finite(force_tangent_world.dot(t0) / mass)?,
                finite(force_tangent_world.dot(t1) / mass)?,
            ],
            arm: n.scale(-scene.radius),
            ledger: Ledger::default(),
            events: Vec::new(),
            attitude_steps: 0,
            sliding_steps: 0,
            velocity_scale: 0.0,
        })
    }
    fn plane_point(&self, time: f64) -> Result<Vec3, SpherePlaneError> {
        let p = self
            .scene
            .point_world
            .add(self.scene.velocity_world.scale(time));
        vector(p)?;
        Ok(p)
    }
    fn gap(&self, body: RigidBodyState, time: f64) -> Result<f64, SpherePlaneError> {
        finite(
            body.pose()
                .position_world()
                .sub(self.plane_point(time)?)
                .dot(self.n)
                - self.scene.radius,
        )
    }
    fn geometry_allowance(&self, body: RigidBodyState, time: f64) -> Result<f64, SpherePlaneError> {
        finite(
            512.0
                * f64::EPSILON
                * (norm(body.pose().position_world())?
                    + norm(self.plane_point(time)?)?
                    + self.scene.radius),
        )
    }
    fn kinematics(&mut self, body: RigidBodyState) -> Result<(f64, [f64; 2]), SpherePlaneError> {
        let velocity = body.center_of_mass_velocity_world(self.scene.properties)?;
        let omega = world_spin(body).scale(self.inertia.recip());
        let relative = velocity
            .sub(self.scene.velocity_world)
            .add(omega.cross(self.arm));
        vector(relative)?;
        self.velocity_scale = self
            .velocity_scale
            .max(norm(velocity)?)
            .max(norm(self.scene.velocity_world)?)
            .max(norm(omega)? * self.scene.radius);
        Ok((
            finite(velocity.sub(self.scene.velocity_world).dot(self.n))?,
            [
                finite(relative.dot(self.tangent[0]))?,
                finite(relative.dot(self.tangent[1]))?,
            ],
        ))
    }
    fn tangent_vector(&self, components: [f64; 2]) -> Vec3 {
        self.tangent[0]
            .scale(components[0])
            .add(self.tangent[1].scale(components[1]))
    }
    fn friction(
        &self,
        slip: [f64; 2],
        normal_impulse: f64,
    ) -> Result<CoulombImpulseResult, SpherePlaneError> {
        Ok(solve_isotropic_coulomb_impulse(
            CoulombImpulseInput {
                normal_impulse,
                relative_tangent_velocity: slip,
                inverse_tangent_mass: self.inverse_tangent_mass,
            },
            &FrictionLaw::Coulomb {
                static_mu: self.scene.friction_coefficient,
                kinetic_mu: self.scene.friction_coefficient,
            },
        )?)
    }
    fn event(
        &mut self,
        time: f64,
        kind: SpherePlaneEventKind,
        normal: f64,
        tangent: Vec3,
    ) -> Result<(), SpherePlaneError> {
        if self.events.len() >= self.config.max_events {
            return Err(SpherePlaneError::EventBudget);
        }
        self.events.push(SpherePlaneEvent {
            time_s: time,
            kind,
            normal_impulse: normal,
            tangent_impulse_world: tangent,
        });
        Ok(())
    }
    fn impact(
        &mut self,
        body: RigidBodyState,
        time: f64,
    ) -> Result<(RigidBodyState, bool), SpherePlaneError> {
        let (normal_speed, slip) = self.kinematics(body)?;
        let normal = solve_normal_impact(NormalImpactInput {
            gap: 0.0,
            activation_tolerance: 0.0,
            relative_normal_speed: normal_speed,
            inverse_effective_mass: self.mass.recip(),
            restitution: self.scene.restitution,
        })?;
        let friction = self.friction(slip, normal.normal_impulse)?;
        let tangent = self.tangent_vector(friction.tangent_impulse);
        let impulse = self.n.scale(normal.normal_impulse).add(tangent);
        let arm_body = body.pose().orientation().rotate_world_to_body(self.arm);
        let result = body.apply_impulse_at_body_point(self.scene.properties, arm_body, impulse)?;
        let angular = body.pose().position_world().add(self.arm).cross(impulse);
        self.ledger.contact = self.ledger.contact.add(impulse);
        self.ledger.contact_angular = self.ledger.contact_angular.add(angular);
        self.ledger.plane_work =
            finite(self.ledger.plane_work + impulse.dot(self.scene.velocity_world))?;
        self.ledger.loss =
            finite(self.ledger.loss + normal.dissipated_energy + friction.dissipated_energy)?;
        self.ledger.impulse_scale = finite(self.ledger.impulse_scale + norm(impulse)?)?;
        self.ledger.angular_scale = finite(self.ledger.angular_scale + norm(angular)?)?;
        let (actual_normal, actual_slip) = self.kinematics(result.state_after)?;
        let allowance = 4096.0 * f64::EPSILON * self.velocity_scale;
        if (actual_normal - normal.relative_normal_speed_after).abs() > allowance
            || (actual_slip[0] - friction.relative_tangent_velocity_after[0]).abs() > allowance
            || (actual_slip[1] - friction.relative_tangent_velocity_after[1]).abs() > allowance
        {
            return Err(SpherePlaneError::Unresolved("impact endpoint velocity"));
        }
        self.event(
            time,
            SpherePlaneEventKind::Impact,
            normal.normal_impulse,
            tangent,
        )?;
        Ok((
            result.state_after,
            normal.relative_normal_speed_after == 0.0,
        ))
    }
    fn segment(
        &mut self,
        body: RigidBodyState,
        h: f64,
        normal_force: f64,
        tangent_force: Vec3,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<RigidBodyState, SpherePlaneError> {
        poll(cancelled)?;
        // Support cancels only the normal applied force. Tangential driving
        // remains in both the analytic trajectory and the external-work ledger.
        let supported = normal_force > 0.0;
        let reaction = self.n.scale(normal_force);
        let contact = reaction.add(tangent_force);
        let total = if supported {
            self.force_tangent_world.add(tangent_force)
        } else {
            self.scene.force_world.add(contact)
        };
        let acceleration = total.scale(self.mass.recip());
        let x0 = body.pose().position_world();
        let v0 = body.center_of_mass_velocity_world(self.scene.properties)?;
        let displacement = v0.add(acceleration.scale(0.5 * h)).scale(h);
        let x1 = x0.add(displacement);
        let p1 = body.linear_momentum_world().add(total.scale(h));
        let l0 = world_spin(body);
        // The normal support line passes through the centre exactly.
        let torque = self.arm.cross(tangent_force);
        let l1 = l0.add(torque.scale(h));
        let mut orientation = body.pose().orientation();
        if commuting_directions(l0, torque) {
            // Commuting world angular velocities have this exact integral,
            // including zero initial spin and torque-free flight.
            let rotation = l0.add(torque.scale(0.5 * h)).scale(h / self.inertia);
            orientation = orientation.right_exp(orientation.rotate_world_to_body(rotation))?;
        } else {
            let steps_float = (h / self.config.max_attitude_step_s).ceil().max(1.0);
            if !steps_float.is_finite()
                || steps_float > (self.config.max_attitude_substeps - self.attitude_steps) as f64
            {
                return Err(SpherePlaneError::AttitudeBudget);
            }
            let steps = steps_float as usize;
            self.attitude_steps += steps;
            let substep = h / steps as f64;
            for step in 0..steps {
                poll(cancelled)?;
                let midpoint = (step as f64 + 0.5) * substep;
                let world_rotation = l0.add(torque.scale(midpoint)).scale(substep / self.inertia);
                orientation =
                    orientation.right_exp(orientation.rotate_world_to_body(world_rotation))?;
            }
        }
        let next = RigidBodyState::new(
            Pose::new(x1, orientation)?,
            p1,
            orientation.rotate_world_to_body(l1),
        )?;
        let contact_impulse = contact.scale(h);
        let external_angular = torque_integral(x0, v0, acceleration, self.scene.force_world, h);
        let contact_angular = torque_integral(x0.add(self.arm), v0, acceleration, contact, h);
        vector(external_angular)?;
        vector(contact_angular)?;
        self.ledger.contact = self.ledger.contact.add(contact_impulse);
        self.ledger.contact_angular = self.ledger.contact_angular.add(contact_angular);
        self.ledger.external_angular = self.ledger.external_angular.add(external_angular);
        self.ledger.force_work =
            finite(self.ledger.force_work + self.scene.force_world.dot(displacement))?;
        self.ledger.plane_work =
            finite(self.ledger.plane_work + contact_impulse.dot(self.scene.velocity_world))?;
        self.ledger.impulse_scale = finite(self.ledger.impulse_scale + norm(contact_impulse)?)?;
        self.ledger.angular_scale =
            finite(self.ledger.angular_scale + norm(external_angular)? + norm(contact_angular)?)?;
        self.kinematics(next)?;
        Ok(next)
    }

    /// Velocity-implicit Coulomb time discretization for resolved curved slip.
    /// The disk owner chooses the full tangent impulse, including spin inertia.
    /// The native segment uses it as a constant force over this step. The
    /// endpoint satisfies Coulomb; continuum direction and arrest need refinement.
    fn curved_sliding(
        &mut self,
        body: RigidBodyState,
        time: f64,
        remaining: f64,
        normal_force: f64,
        slip: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(RigidBodyState, f64, bool), SpherePlaneError> {
        let slip_speed = finite(slip[0].hypot(slip[1]))?;
        let acceleration = self.acceleration_tangent;
        let acceleration_speed = finite(acceleration[0].hypot(acceleration[1]))?;
        let capacity = finite(self.scene.friction_coefficient * normal_force)?;
        let can_hold = self.inverse_tangent_mass > 1.0
            && capacity > f64::MAX / self.inverse_tangent_mass
            || acceleration_speed <= self.inverse_tangent_mass * capacity;
        let mut h = remaining.min(self.config.max_sliding_step_s);
        let (friction, arrested) = loop {
            poll(cancelled)?;
            if self.sliding_steps >= self.config.max_sliding_substeps {
                return Err(SpherePlaneError::SlidingBudget);
            }
            self.sliding_steps += 1;
            progress(time, h)?;
            let free_slip = [
                finite(acceleration[0].mul_add(h, slip[0]))?,
                finite(acceleration[1].mul_add(h, slip[1]))?,
            ];
            let friction = self.friction(free_slip, finite(normal_force * h)?)?;
            let after = friction.relative_tangent_velocity_after;
            let after_speed = finite(after[0].hypot(after[1]))?;
            let arrested = after_speed <= 4096.0 * f64::EPSILON * self.velocity_scale;
            // Resolve a turning trajectory without allowing the discrete force
            // to do positive work on an initial part of its linear slip path.
            // A supercritical transverse drive cannot arrest at a finite time;
            // a coarse disk projection to zero must therefore be refined.
            // Check arrested trials too: an opposing drive can carry the free
            // velocity through zero and make a coarse stopping impulse assist
            // the initial slip. Endpoint arrest alone does not make that passive.
            let impulse_speed = friction.tangent_impulse[0].hypot(friction.tangent_impulse[1]);
            let positive_initial_work = impulse_speed > 0.0
                && (friction.tangent_impulse[0] / impulse_speed).mul_add(
                    slip[0] / slip_speed,
                    (friction.tangent_impulse[1] / impulse_speed) * (slip[1] / slip_speed),
                ) > 0.0;
            if positive_initial_work || (arrested && !can_hold) {
                h *= 0.5;
                continue;
            }
            break (friction, arrested);
        };
        let impulse = self.tangent_vector(friction.tangent_impulse);
        let tangent_force = impulse.scale(h.recip());
        vector(tangent_force)?;
        let next = self.segment(body, h, normal_force, tangent_force, cancelled)?;
        let (_, after) = self.kinematics(next)?;
        let allowance = 8192.0 * f64::EPSILON * self.velocity_scale;
        if (after[0] - friction.relative_tangent_velocity_after[0]).abs() > allowance
            || (after[1] - friction.relative_tangent_velocity_after[1]).abs() > allowance
            || (arrested && after[0].hypot(after[1]) > allowance)
        {
            return Err(SpherePlaneError::Unresolved(
                "curved sliding endpoint velocity",
            ));
        }
        let midpoint = [
            0.5 * slip[0] + 0.5 * after[0],
            0.5 * slip[1] + 0.5 * after[1],
        ];
        let loss = finite(-impulse.dot(self.tangent_vector(midpoint)))?;
        if loss < 0.0 {
            return Err(SpherePlaneError::Unresolved("curved sliding friction work"));
        }
        self.ledger.loss = finite(self.ledger.loss + loss)?;

        // Normal support is exact in the declared equations. Also bound the
        // complete native quadratic path's floating geometry residual, not
        // merely its endpoint, without changing the pose or momentum.
        let gap0 = self.gap(body, time)?;
        let gap1 = self.gap(next, finite(time + h)?)?;
        let normal_velocity = body
            .center_of_mass_velocity_world(self.scene.properties)?
            .sub(self.scene.velocity_world)
            .dot(self.n);
        let normal_acceleration =
            self.force_tangent_world.add(tangent_force).dot(self.n) / self.mass;
        let extremum = -normal_velocity / normal_acceleration;
        let mut maximum_gap = gap0.abs().max(gap1.abs());
        if extremum > 0.0 && extremum < h {
            maximum_gap = maximum_gap.max(
                finite(gap0 + extremum * (normal_velocity + 0.5 * normal_acceleration * extremum))?
                    .abs(),
            );
        }
        let geometry_allowance = self
            .geometry_allowance(body, time)?
            .max(self.geometry_allowance(next, finite(time + h)?)?);
        if maximum_gap > geometry_allowance {
            return Err(SpherePlaneError::Unresolved(
                "curved sliding supported path",
            ));
        }
        Ok((next, h, arrested))
    }
}

/// Advance through detected impacts, normal support, Coulomb sliding and rolling.
///
/// Newly constructed penetrating states always refuse. Only a returned state's
/// private known-contact mode may carry a bounded floating-point geometry/zero-
/// velocity residual through continuation; no position or momentum is projected.
/// Rebound accumulation or an unresolvable event time refuses at the event/work
/// boundary, rather than silently changing restitution to stick.
///
/// Arbitrary spin, tangential velocity and constant centre force are admitted
/// in flight. Persistent contact admits arbitrary driven rolling, and sliding
/// when the tangent applied force is collinear with slip (including opposing
/// force and resolved reversal). Changing-direction sliding uses explicitly
/// bounded, first-order velocity-implicit Coulomb steps and resolves arrest at
/// a timestep boundary. Refinement is required for its continuum trajectory;
/// no exact curved slip-arrest time is claimed. Static traction keeps
/// zero slip when its required magnitude fits the Coulomb limit. Sliding arrest
/// and force-driven slip onset are resolved explicitly. Loss is measured from
/// actual midpoint slip work; curved-step force-direction error is part of the
/// declared discretization and is not assigned to thermal cells. Attitude is
/// exact for commuting spin/torque; the general
/// changing-axis case uses the declared second-order attitude substeps.
/// Spin about the normal has no torsional/rolling friction.
/// A constant-velocity plane supplies explicitly signed work. This is floating
/// analytic event detection, not interval-certified general-purpose CCD.
pub fn advance_sphere_plane(
    scene: &SpherePlane,
    state: SpherePlaneState,
    duration_s: f64,
    config: SpherePlaneConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SpherePlaneAdvance, SpherePlaneError> {
    poll(&mut cancelled)?;
    if !(duration_s.is_finite() && duration_s > 0.0) {
        return Err(SpherePlaneError::Invalid(
            "positive finite duration required",
        ));
    }
    let end_time = finite(state.time_s + duration_s)?;
    if end_time == state.time_s {
        return Err(SpherePlaneError::NoProgress);
    }
    let mut engine = Engine::new(*scene, config)?;
    engine.velocity_scale = state.velocity_scale;
    let initial_energy = state.body.kinetic_energy(scene.properties)?;
    let mut body = state.body;
    let mut mode = state.mode;
    let raw_gap = engine.gap(body, state.time_s)?;
    let geometry_allowance = engine.geometry_allowance(body, state.time_s)?;
    if raw_gap < 0.0 && (!state.contact_known || -raw_gap > geometry_allowance) {
        return Err(SpherePlaneError::Penetration { gap: raw_gap });
    }
    let mut touching =
        raw_gap == 0.0 || (state.contact_known && raw_gap.abs() <= geometry_allowance);
    let mut zero_normal = touching && state.zero_normal;
    if !touching {
        mode = SpherePlaneMode::Flight;
    }
    let mut elapsed = 0.0;
    while elapsed < duration_s {
        poll(&mut cancelled)?;
        let remaining = duration_s - elapsed;
        let time = finite(state.time_s + elapsed)?;
        let (normal_raw, slip) = engine.kinematics(body)?;
        let velocity_allowance = 4096.0 * f64::EPSILON * engine.velocity_scale;
        let normal_speed = if zero_normal && normal_raw.abs() <= velocity_allowance {
            0.0
        } else {
            normal_raw
        };
        let slip_speed = finite(slip[0].hypot(slip[1]))?;
        if touching && normal_speed < 0.0 {
            (body, zero_normal) = engine.impact(body, time)?;
            mode = SpherePlaneMode::Flight;
            if zero_normal && engine.acceleration_normal <= 0.0 {
                let (_, after_slip) = engine.kinematics(body)?;
                mode = if after_slip[0].hypot(after_slip[1])
                    <= 4096.0 * f64::EPSILON * engine.velocity_scale
                {
                    SpherePlaneMode::Rolling
                } else {
                    SpherePlaneMode::Sliding
                };
            }
            continue;
        }
        if touching && normal_speed == 0.0 && engine.acceleration_normal <= 0.0 {
            zero_normal = true;
            if mode == SpherePlaneMode::Flight {
                engine.event(time, SpherePlaneEventKind::Support, 0.0, Vec3::ZERO)?;
            }
            let normal_force = finite(-scene.force_world.dot(engine.n))?;
            let rolling = slip_speed == 0.0
                || (mode == SpherePlaneMode::Rolling && slip_speed <= velocity_allowance);
            let friction_force = finite(scene.friction_coefficient * normal_force)?;
            let acceleration = engine.acceleration_tangent;
            let acceleration_speed = finite(acceleration[0].hypot(acceleration[1]))?;
            if friction_force == 0.0 {
                // No traction means no restriction on the direction of driving.
                // Integrate the actual force, even when starting with zero slip.
                if rolling && acceleration_speed > 0.0 {
                    engine.event(time, SpherePlaneEventKind::Sliding, 0.0, Vec3::ZERO)?;
                }
                progress(time, remaining)?;
                body = engine.segment(body, remaining, normal_force, Vec3::ZERO, &mut cancelled)?;
                let (_, after_slip) = engine.kinematics(body)?;
                mode = if after_slip[0] == 0.0 && after_slip[1] == 0.0
                    || (rolling
                        && acceleration_speed == 0.0
                        && after_slip[0].hypot(after_slip[1])
                            <= 8192.0 * f64::EPSILON * engine.velocity_scale)
                {
                    SpherePlaneMode::Rolling
                } else {
                    SpherePlaneMode::Sliding
                };
                elapsed = duration_s;
                continue;
            }
            // Compare in acceleration units without overflowing a large finite
            // traction capacity. Only the finite demand is used for rolling.
            let holds_rolling = rolling
                && (engine.inverse_tangent_mass > 1.0
                    && friction_force > f64::MAX / engine.inverse_tangent_mass
                    || acceleration_speed <= engine.inverse_tangent_mass * friction_force);
            let mut stopping_time = f64::INFINITY;
            let mut transverse_acceleration = 0.0;
            if !holds_rolling {
                let direction = if rolling {
                    if acceleration_speed == 0.0 {
                        return Err(SpherePlaneError::Unresolved("sliding onset direction"));
                    }
                    engine.event(time, SpherePlaneEventKind::Sliding, 0.0, Vec3::ZERO)?;
                    [
                        acceleration[0] / acceleration_speed,
                        acceleration[1] / acceleration_speed,
                    ]
                } else {
                    [slip[0] / slip_speed, slip[1] / slip_speed]
                };
                if acceleration_speed > 0.0 && !rolling {
                    // Dimensionless collinearity admits rotation roundoff only;
                    // resolved transverse drive uses the bounded curved path.
                    let turn = direction[0].mul_add(
                        acceleration[1] / acceleration_speed,
                        -direction[1] * (acceleration[0] / acceleration_speed),
                    );
                    if turn.abs() > 256.0 * f64::EPSILON {
                        transverse_acceleration = finite(turn.abs() * acceleration_speed)?;
                    }
                }
                let along =
                    finite(direction[0].mul_add(acceleration[0], direction[1] * acceleration[1]))?;
                let slip_rate =
                    finite((-engine.inverse_tangent_mass).mul_add(friction_force, along))?;
                if rolling && slip_rate <= 0.0 {
                    return Err(SpherePlaneError::Unresolved("sliding onset acceleration"));
                }
                if !rolling && slip_rate < 0.0 {
                    // An infinite stop time is beyond every finite horizon.
                    stopping_time = slip_speed / -slip_rate;
                }
                mode = SpherePlaneMode::Sliding;
            } else {
                mode = SpherePlaneMode::Rolling;
            }
            let h = remaining.min(stopping_time);
            progress(elapsed, h)?;
            progress(time, h)?;
            // A nominally normal force may leave an arbitrary tiny tangent
            // direction after rotation. Admit it only below the velocity
            // resolution of this segment, never as a resolved curved path.
            if transverse_acceleration > (256.0 * f64::EPSILON * slip_speed) / h {
                let (next, step, arrested) = engine.curved_sliding(
                    body,
                    time,
                    remaining,
                    normal_force,
                    slip,
                    &mut cancelled,
                )?;
                progress(elapsed, step)?;
                body = next;
                elapsed = if step == remaining {
                    duration_s
                } else {
                    elapsed + step
                };
                if arrested {
                    mode = SpherePlaneMode::Rolling;
                    engine.event(
                        finite(state.time_s + elapsed)?,
                        SpherePlaneEventKind::Rolling,
                        0.0,
                        Vec3::ZERO,
                    )?;
                }
                continue;
            }
            // The existing disk owner selects the integrated traction from
            // the force-driven free slip. Event splitting ensures this impulse
            // corresponds to one constant-direction continuous-force segment.
            let admitted_slip = if rolling { [0.0; 2] } else { slip };
            let free_slip = [
                finite(acceleration[0].mul_add(h, admitted_slip[0]))?,
                finite(acceleration[1].mul_add(h, admitted_slip[1]))?,
            ];
            let friction = engine.friction(free_slip, finite(normal_force * h)?)?;
            let tangent_force = engine
                .tangent_vector(friction.tangent_impulse)
                .scale(h.recip());
            vector(tangent_force)?;
            body = engine.segment(body, h, normal_force, tangent_force, &mut cancelled)?;
            let (_, after_slip) = engine.kinematics(body)?;
            if !holds_rolling {
                // Impulse-law loss at free_slip includes fictitious splitting
                // work. Physical continuous loss uses actual endpoint slip.
                let impulse = tangent_force.scale(h);
                let midpoint = [
                    0.5 * slip[0] + 0.5 * after_slip[0],
                    0.5 * slip[1] + 0.5 * after_slip[1],
                ];
                let loss = finite(-impulse.dot(engine.tangent_vector(midpoint)))?;
                if loss < 0.0 {
                    return Err(SpherePlaneError::Unresolved("continuous friction work"));
                }
                engine.ledger.loss = finite(engine.ledger.loss + loss)?;
            }
            elapsed = if h == remaining {
                duration_s
            } else {
                elapsed + h
            };
            if holds_rolling || stopping_time <= remaining {
                if after_slip[0].hypot(after_slip[1])
                    > 8192.0 * f64::EPSILON * engine.velocity_scale
                {
                    return Err(SpherePlaneError::Unresolved("supported rolling velocity"));
                }
                mode = SpherePlaneMode::Rolling;
                if !holds_rolling {
                    engine.event(
                        finite(state.time_s + elapsed)?,
                        SpherePlaneEventKind::Rolling,
                        0.0,
                        Vec3::ZERO,
                    )?;
                }
            }
            continue;
        }
        if mode != SpherePlaneMode::Flight {
            engine.event(time, SpherePlaneEventKind::LiftOff, 0.0, Vec3::ZERO)?;
        }
        mode = SpherePlaneMode::Flight;
        let gap = if touching {
            0.0
        } else {
            engine.gap(body, time)?
        };
        let event = first_ballistic_contact(BallisticContactInput {
            gap,
            relative_normal_speed: normal_speed,
            relative_normal_acceleration: engine.acceleration_normal,
            duration_s: remaining,
        })?;
        let Some(event) = event else {
            progress(time, remaining)?;
            body = engine.segment(body, remaining, 0.0, Vec3::ZERO, &mut cancelled)?;
            elapsed = duration_s;
            touching = false;
            zero_normal = false;
            continue;
        };
        if event.time_s > 0.0 {
            progress(elapsed, event.time_s)?;
            progress(time, event.time_s)?;
            body = engine.segment(body, event.time_s, 0.0, Vec3::ZERO, &mut cancelled)?;
            elapsed = if event.time_s == remaining {
                duration_s
            } else {
                elapsed + event.time_s
            };
        } else if !touching && gap > 0.0 {
            return Err(SpherePlaneError::NoProgress);
        }
        let event_time = finite(state.time_s + elapsed)?;
        if engine.gap(body, event_time)?.abs() > engine.geometry_allowance(body, event_time)? {
            return Err(SpherePlaneError::Unresolved("detected contact geometry"));
        }
        touching = true;
        zero_normal = event.relative_normal_speed == 0.0;
        if event.relative_normal_speed < 0.0 {
            (body, zero_normal) = engine.impact(body, finite(state.time_s + elapsed)?)?;
            if zero_normal && engine.acceleration_normal <= 0.0 {
                let (_, after_slip) = engine.kinematics(body)?;
                mode = if after_slip[0].hypot(after_slip[1])
                    <= 4096.0 * f64::EPSILON * engine.velocity_scale
                {
                    SpherePlaneMode::Rolling
                } else {
                    SpherePlaneMode::Sliding
                };
            }
        } else if engine.acceleration_normal < 0.0 {
            engine.event(
                finite(state.time_s + elapsed)?,
                SpherePlaneEventKind::Support,
                0.0,
                Vec3::ZERO,
            )?;
            mode = if slip_speed == 0.0 {
                SpherePlaneMode::Rolling
            } else {
                SpherePlaneMode::Sliding
            };
        }
        // Grazing with outward acceleration returns to flight; the next
        // query starts with gap=0, normal speed=0 and therefore cannot loop.
    }
    poll(&mut cancelled)?;
    let gap_after = engine.gap(body, end_time)?;
    let geometry_allowance = engine.geometry_allowance(body, end_time)?;
    if (touching && gap_after.abs() > geometry_allowance) || (!touching && gap_after < 0.0) {
        return Err(SpherePlaneError::Unresolved("trajectory geometric gap"));
    }
    let final_energy = body.kinetic_energy(scene.properties)?;
    let ledger = engine.ledger;
    let energy_balance_residual = finite(
        (final_energy - initial_energy) - ledger.force_work - ledger.plane_work + ledger.loss,
    )?;
    let linear_momentum_residual = body
        .linear_momentum_world()
        .sub(state.body.linear_momentum_world())
        .sub(scene.force_world.scale(duration_s))
        .sub(ledger.contact);
    let angular_momentum_residual = total_angular(body)
        .sub(total_angular(state.body))
        .sub(ledger.external_angular)
        .sub(ledger.contact_angular);
    vector(linear_momentum_residual)?;
    vector(angular_momentum_residual)?;
    let tolerance = 8192.0 * f64::EPSILON * (engine.events.len() + 1) as f64;
    let energy_scale = initial_energy.abs()
        + final_energy.abs()
        + ledger.force_work.abs()
        + ledger.plane_work.abs()
        + ledger.loss;
    let linear_scale = norm(state.body.linear_momentum_world())?
        + norm(body.linear_momentum_world())?
        + norm(scene.force_world)? * duration_s
        + ledger.impulse_scale;
    let angular_scale =
        angular_magnitude(state.body)? + angular_magnitude(body)? + ledger.angular_scale;
    if energy_balance_residual.abs() > finite(tolerance * energy_scale)?
        || norm(linear_momentum_residual)? > finite(tolerance * linear_scale)?
        || norm(angular_momentum_residual)? > finite(tolerance * angular_scale)?
    {
        return Err(SpherePlaneError::Unresolved(
            "global momentum or energy balance",
        ));
    }
    Ok(SpherePlaneAdvance {
        state_after: SpherePlaneState {
            body,
            time_s: end_time,
            mode,
            contact_known: touching,
            zero_normal,
            velocity_scale: engine.velocity_scale,
        },
        events: engine.events,
        contact_impulse_world: ledger.contact,
        contact_angular_impulse_world: ledger.contact_angular,
        external_force_work: ledger.force_work,
        plane_work: ledger.plane_work,
        dissipated_energy: ledger.loss,
        energy_balance_residual,
        linear_momentum_residual,
        angular_momentum_residual,
        gap_after,
        sliding_substeps: engine.sliding_steps,
    })
}
// Compare directions, not dimensionful products: a common length or mass
// rescaling must not turn noncommuting angular motion into a zero cross product.
fn commuting_directions(spin: Vec3, torque: Vec3) -> bool {
    if spin == Vec3::ZERO || torque == Vec3::ZERO {
        return true;
    }
    let spin_scale = spin.x.abs().max(spin.y.abs()).max(spin.z.abs());
    let torque_scale = torque.x.abs().max(torque.y.abs()).max(torque.z.abs());
    let a = Vec3::new(
        spin.x / spin_scale,
        spin.y / spin_scale,
        spin.z / spin_scale,
    );
    let b = Vec3::new(
        torque.x / torque_scale,
        torque.y / torque_scale,
        torque.z / torque_scale,
    );
    a.cross(b) == Vec3::ZERO
}
fn world_spin(body: RigidBodyState) -> Vec3 {
    body.pose()
        .orientation()
        .rotate_body_to_world(body.angular_momentum_body())
}
fn total_angular(body: RigidBodyState) -> Vec3 {
    body.pose()
        .position_world()
        .cross(body.linear_momentum_world())
        .add(world_spin(body))
}
fn angular_magnitude(body: RigidBodyState) -> Result<f64, SpherePlaneError> {
    finite(
        norm(
            body.pose()
                .position_world()
                .cross(body.linear_momentum_world()),
        )? + norm(world_spin(body))?,
    )
}
fn torque_integral(x: Vec3, v: Vec3, a: Vec3, force: Vec3, h: f64) -> Vec3 {
    x.cross(force)
        .scale(h)
        .add(v.cross(force).scale(0.5 * h).scale(h))
        .add(a.cross(force).scale(h / 6.0).scale(h).scale(h))
}
fn norm(value: Vec3) -> Result<f64, SpherePlaneError> {
    Ok(value.stable_norm("sphere vector magnitude")?)
}
fn vector(value: Vec3) -> Result<(), SpherePlaneError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(SpherePlaneError::Invalid("nonfinite sphere arithmetic"))
    }
}
fn finite(value: f64) -> Result<f64, SpherePlaneError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(SpherePlaneError::Invalid("nonfinite sphere arithmetic"))
    }
}
fn progress(elapsed: f64, h: f64) -> Result<(), SpherePlaneError> {
    if h > 0.0 && elapsed + h > elapsed {
        Ok(())
    } else {
        Err(SpherePlaneError::NoProgress)
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), SpherePlaneError> {
    if cancelled() {
        Err(SpherePlaneError::Cancelled)
    } else {
        Ok(())
    }
}
