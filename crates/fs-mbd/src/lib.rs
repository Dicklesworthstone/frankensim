//! Deterministic, safe-Rust foundations for unconstrained rigid-body dynamics.
//!
//! This crate deliberately owns only smooth, unconstrained rigid-body motion.
//! Contact, joints, holonomic/nonholonomic constraints, impacts, and their
//! receipts require the `fs-contact`, `fs-kinematics`, and `fs-solver` lanes;
//! they are not approximated here. Orientations are canonical unit quaternions
//! mapping body vectors into the world frame. The canonical sign convention
//! removes the quaternion double-cover ambiguity deterministically.
//!
//! The fixed-step integrator uses a deterministic midpoint update: translation
//! under a piecewise-constant world force and a Lie-group attitude update at
//! midpoint body angular momentum. It is a production foundation for smooth
//! rigid bodies, not a claim of variational, contact, or constraint preservation.

#![deny(unsafe_code)]
#![warn(missing_docs)]

use core::fmt;

/// A three-dimensional vector in an explicitly documented frame.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec3 {
    /// x component.
    pub x: f64,
    /// y component.
    pub y: f64,
    /// z component.
    pub z: f64,
}

impl Vec3 {
    /// The zero vector.
    pub const ZERO: Self = Self::new(0.0, 0.0, 0.0);

    /// Creates a vector.
    #[must_use]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Returns whether every component is finite.
    #[must_use]
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    /// Dot product.
    #[must_use]
    pub fn dot(self, other: Self) -> f64 {
        self.x
            .mul_add(other.x, self.y.mul_add(other.y, self.z * other.z))
    }

    /// Cross product, preserving the operand order.
    #[must_use]
    pub fn cross(self, other: Self) -> Self {
        Self::new(
            self.y.mul_add(other.z, -(self.z * other.y)),
            self.z.mul_add(other.x, -(self.x * other.z)),
            self.x.mul_add(other.y, -(self.y * other.x)),
        )
    }

    /// Squared Euclidean norm.
    #[must_use]
    pub fn norm_squared(self) -> f64 {
        self.dot(self)
    }

    fn stable_norm(self, field: &'static str) -> Result<f64, DynamicsError> {
        if !self.is_finite() {
            return Err(DynamicsError::NonFinite(field));
        }
        let scale = self.x.abs().max(self.y.abs()).max(self.z.abs());
        if scale == 0.0 {
            return Ok(0.0);
        }
        let scaled = Self::new(self.x / scale, self.y / scale, self.z / scale);
        let magnitude = scale * scaled.dot(scaled).sqrt();
        if !magnitude.is_finite() {
            return Err(DynamicsError::UnrepresentableMagnitude(field));
        }
        Ok(magnitude)
    }

    fn normalized(self, field: &'static str) -> Result<Self, DynamicsError> {
        if !self.is_finite() {
            return Err(DynamicsError::NonFinite(field));
        }
        let scale = self.x.abs().max(self.y.abs()).max(self.z.abs());
        if scale == 0.0 {
            return Err(DynamicsError::InvalidOrientation);
        }
        let scaled = Self::new(self.x / scale, self.y / scale, self.z / scale);
        let inverse_norm = scaled.dot(scaled).sqrt().recip();
        Ok(scaled.scale(inverse_norm))
    }

    /// Scalar multiplication.
    #[must_use]
    pub fn scale(self, scalar: f64) -> Self {
        Self::new(scalar * self.x, scalar * self.y, scalar * self.z)
    }

    /// Componentwise sum.
    #[must_use]
    pub fn add(self, other: Self) -> Self {
        Self::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }

    /// Componentwise difference.
    #[must_use]
    pub fn sub(self, other: Self) -> Self {
        Self::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

/// An input or state error that prevents an admitted dynamics step.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DynamicsError {
    /// A named scalar or vector contained a NaN or infinity.
    NonFinite(&'static str),
    /// Mass must be finite and strictly positive.
    InvalidMass,
    /// A principal moment must be finite and strictly positive.
    InvalidPrincipalInertia,
    /// Principal moments fail a rigid-body triangle inequality.
    InconsistentPrincipalInertia,
    /// A quaternion with zero or non-finite norm cannot define an orientation.
    InvalidOrientation,
    /// A finite vector's mathematical magnitude cannot be represented in `f64`.
    UnrepresentableMagnitude(&'static str),
    /// This core admits center-of-mass reference frames only.
    UnsupportedReferenceOffset,
    /// A step duration must be finite and strictly positive.
    InvalidStepDuration,
}

impl fmt::Display for DynamicsError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFinite(field) => write!(formatter, "{field} must be finite"),
            Self::InvalidMass => formatter.write_str("mass must be finite and positive"),
            Self::InvalidPrincipalInertia => {
                formatter.write_str("each principal moment must be finite and positive")
            }
            Self::InconsistentPrincipalInertia => formatter
                .write_str("principal moments must satisfy rigid-body triangle inequalities"),
            Self::InvalidOrientation => {
                formatter.write_str("orientation must have a finite nonzero norm")
            }
            Self::UnrepresentableMagnitude(field) => {
                write!(formatter, "{field} magnitude is not representable as f64")
            }
            Self::UnsupportedReferenceOffset => formatter
                .write_str("this rigid-body core requires a center-of-mass reference point"),
            Self::InvalidStepDuration => {
                formatter.write_str("step duration must be finite and positive")
            }
        }
    }
}

impl std::error::Error for DynamicsError {}

/// A normalized quaternion mapping body-frame vectors into the world frame.
///
/// Its fields are private so an invalid orientation cannot enter a body state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UnitQuaternion {
    w: f64,
    x: f64,
    y: f64,
    z: f64,
}

impl UnitQuaternion {
    /// The identity orientation.
    pub const IDENTITY: Self = Self {
        w: 1.0,
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    /// Normalizes and canonically signs a quaternion.
    ///
    /// The first nonzero component in `(w, x, y, z)` is always positive,
    /// providing a deterministic representative of the double cover.
    pub fn new(w: f64, x: f64, y: f64, z: f64) -> Result<Self, DynamicsError> {
        if !w.is_finite() || !x.is_finite() || !y.is_finite() || !z.is_finite() {
            return Err(DynamicsError::InvalidOrientation);
        }
        let scale = w.abs().max(x.abs()).max(y.abs()).max(z.abs());
        if scale == 0.0 {
            return Err(DynamicsError::InvalidOrientation);
        }
        let scaled_w = w / scale;
        let scaled_x = x / scale;
        let scaled_y = y / scale;
        let scaled_z = z / scale;
        let inverse_norm = scaled_w
            .mul_add(
                scaled_w,
                scaled_x.mul_add(scaled_x, scaled_y.mul_add(scaled_y, scaled_z * scaled_z)),
            )
            .sqrt()
            .recip();
        Ok(Self::canonical(
            scaled_w * inverse_norm,
            scaled_x * inverse_norm,
            scaled_y * inverse_norm,
            scaled_z * inverse_norm,
        ))
    }

    /// Builds an orientation from a body-frame axis and angle in radians.
    pub fn from_axis_angle(axis_body: Vec3, angle_radians: f64) -> Result<Self, DynamicsError> {
        if !axis_body.is_finite() || !angle_radians.is_finite() {
            return Err(DynamicsError::NonFinite("axis or angle"));
        }
        let axis = axis_body.normalized("axis_body")?;
        let half_angle = 0.5 * angle_radians;
        Self::new(
            half_angle.cos(),
            half_angle.sin() * axis.x,
            half_angle.sin() * axis.y,
            half_angle.sin() * axis.z,
        )
    }

    /// Returns quaternion components in `(w, x, y, z)` order.
    #[must_use]
    pub const fn components(self) -> [f64; 4] {
        [self.w, self.x, self.y, self.z]
    }

    /// Rotates a body-frame vector into the world frame.
    #[must_use]
    pub fn rotate_body_to_world(self, vector_body: Vec3) -> Vec3 {
        let vector_quaternion = Self {
            w: 0.0,
            x: vector_body.x,
            y: vector_body.y,
            z: vector_body.z,
        };
        let rotated = self.multiply(vector_quaternion).multiply(self.conjugate());
        Vec3::new(rotated.x, rotated.y, rotated.z)
    }

    /// Right-composes a body-frame rotation vector through the exponential map.
    #[must_use]
    pub fn right_exp(self, rotation_vector_body: Vec3) -> Result<Self, DynamicsError> {
        if !rotation_vector_body.is_finite() {
            return Err(DynamicsError::NonFinite("rotation_vector_body"));
        }
        let half = rotation_vector_body.scale(0.5);
        let theta = half.stable_norm("rotation_vector_body")?;
        let (cosine, sine_over_theta) = if theta < 1e-6 {
            let theta_squared = theta * theta;
            (
                1.0 - 0.5 * theta_squared + theta_squared * theta_squared / 24.0,
                1.0 - theta_squared / 6.0 + theta_squared * theta_squared / 120.0,
            )
        } else {
            (theta.cos(), theta.sin() / theta)
        };
        let step = Self {
            w: cosine,
            x: sine_over_theta * half.x,
            y: sine_over_theta * half.y,
            z: sine_over_theta * half.z,
        };
        let product = self.multiply(step);
        Self::new(product.w, product.x, product.y, product.z)
    }

    fn canonical(w: f64, x: f64, y: f64, z: f64) -> Self {
        let sign = if w < 0.0
            || (w == 0.0 && (x < 0.0 || (x == 0.0 && (y < 0.0 || (y == 0.0 && z < 0.0)))))
        {
            -1.0
        } else {
            1.0
        };
        Self {
            w: sign * w,
            x: sign * x,
            y: sign * y,
            z: sign * z,
        }
    }

    fn conjugate(self) -> Self {
        Self {
            w: self.w,
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }

    fn multiply(self, other: Self) -> Self {
        Self {
            w: self.w.mul_add(
                other.w,
                -(self.x * other.x) - self.y * other.y - self.z * other.z,
            ),
            x: self.w.mul_add(
                other.x,
                self.x
                    .mul_add(other.w, self.y.mul_add(other.z, -(self.z * other.y))),
            ),
            y: self.w.mul_add(
                other.y,
                self.y
                    .mul_add(other.w, self.z.mul_add(other.x, -(self.x * other.z))),
            ),
            z: self.w.mul_add(
                other.z,
                self.z
                    .mul_add(other.w, self.x.mul_add(other.y, -(self.y * other.x))),
            ),
        }
    }
}

/// The world position and body-to-world orientation of a body reference point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Pose {
    /// Reference-point location in the world frame.
    position_world: Vec3,
    /// Body-to-world orientation.
    orientation: UnitQuaternion,
}

impl Pose {
    /// Constructs a center-of-mass pose from a finite world position and a
    /// previously validated body-to-world orientation.
    pub fn new(position_world: Vec3, orientation: UnitQuaternion) -> Result<Self, DynamicsError> {
        if !position_world.is_finite() {
            return Err(DynamicsError::NonFinite("pose.position_world"));
        }
        Ok(Self {
            position_world,
            orientation,
        })
    }

    /// The world-origin pose.
    #[must_use]
    pub const fn identity() -> Self {
        Self {
            position_world: Vec3::ZERO,
            orientation: UnitQuaternion::IDENTITY,
        }
    }

    /// Returns the center-of-mass position in the world frame.
    #[must_use]
    pub const fn position_world(self) -> Vec3 {
        self.position_world
    }

    /// Returns the body-to-world orientation.
    #[must_use]
    pub const fn orientation(self) -> UnitQuaternion {
        self.orientation
    }

    fn validate(self) -> Result<(), DynamicsError> {
        if !self.position_world.is_finite() {
            return Err(DynamicsError::NonFinite("pose.position_world"));
        }
        Ok(())
    }
}

/// Validated mass and diagonal inertia tensor in its principal body frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MassProperties {
    mass: f64,
    center_of_mass_body: Vec3,
    principal_inertia_body: Vec3,
}

impl MassProperties {
    /// Validates mass, center of mass, and principal moments.
    ///
    /// This foundation admits only a center-of-mass reference point, so
    /// `center_of_mass_body` must be exactly zero. Offset-reference spatial
    /// inertia is intentionally deferred rather than silently mishandled.
    pub fn new(
        mass: f64,
        center_of_mass_body: Vec3,
        principal_inertia_body: Vec3,
    ) -> Result<Self, DynamicsError> {
        if !mass.is_finite() || mass <= 0.0 {
            return Err(DynamicsError::InvalidMass);
        }
        if !center_of_mass_body.is_finite() {
            return Err(DynamicsError::NonFinite("center_of_mass_body"));
        }
        if center_of_mass_body != Vec3::ZERO {
            return Err(DynamicsError::UnsupportedReferenceOffset);
        }
        if !principal_inertia_body.is_finite()
            || principal_inertia_body.x <= 0.0
            || principal_inertia_body.y <= 0.0
            || principal_inertia_body.z <= 0.0
        {
            return Err(DynamicsError::InvalidPrincipalInertia);
        }
        let moments = principal_inertia_body;
        if moments.x > moments.y + moments.z
            || moments.y > moments.x + moments.z
            || moments.z > moments.x + moments.y
        {
            return Err(DynamicsError::InconsistentPrincipalInertia);
        }
        Ok(Self {
            mass,
            center_of_mass_body,
            principal_inertia_body,
        })
    }

    /// Returns mass in kilograms.
    #[must_use]
    pub const fn mass(self) -> f64 {
        self.mass
    }

    /// Returns the declared center of mass in the principal body frame.
    #[must_use]
    pub const fn center_of_mass_body(self) -> Vec3 {
        self.center_of_mass_body
    }

    /// Returns principal moments `(Ixx, Iyy, Izz)` in kg m².
    #[must_use]
    pub const fn principal_inertia_body(self) -> Vec3 {
        self.principal_inertia_body
    }

    /// Converts body angular momentum into body angular velocity.
    #[must_use]
    pub fn angular_velocity_body(self, angular_momentum_body: Vec3) -> Vec3 {
        Vec3::new(
            angular_momentum_body.x / self.principal_inertia_body.x,
            angular_momentum_body.y / self.principal_inertia_body.y,
            angular_momentum_body.z / self.principal_inertia_body.z,
        )
    }

    fn validate(self) -> Result<(), DynamicsError> {
        Self::new(
            self.mass,
            self.center_of_mass_body,
            self.principal_inertia_body,
        )
        .map(|_| ())
    }
}

/// Momentum-form state of one unconstrained rigid body.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidBodyState {
    /// Center-of-mass pose.
    pose: Pose,
    /// Linear momentum in kg m/s, expressed in the world frame.
    linear_momentum_world: Vec3,
    /// Angular momentum in kg m²/s, expressed in the principal body frame.
    angular_momentum_body: Vec3,
}

impl RigidBodyState {
    /// Validates a state with an already-valid pose.
    pub fn new(
        pose: Pose,
        linear_momentum_world: Vec3,
        angular_momentum_body: Vec3,
    ) -> Result<Self, DynamicsError> {
        pose.validate()?;
        if !linear_momentum_world.is_finite() {
            return Err(DynamicsError::NonFinite("linear_momentum_world"));
        }
        if !angular_momentum_body.is_finite() {
            return Err(DynamicsError::NonFinite("angular_momentum_body"));
        }
        Ok(Self {
            pose,
            linear_momentum_world,
            angular_momentum_body,
        })
    }

    /// Returns the center-of-mass pose.
    #[must_use]
    pub const fn pose(self) -> Pose {
        self.pose
    }

    /// Returns world-frame linear momentum.
    #[must_use]
    pub const fn linear_momentum_world(self) -> Vec3 {
        self.linear_momentum_world
    }

    /// Returns principal-body-frame angular momentum.
    #[must_use]
    pub const fn angular_momentum_body(self) -> Vec3 {
        self.angular_momentum_body
    }

    fn validate(self) -> Result<(), DynamicsError> {
        Self::new(
            self.pose,
            self.linear_momentum_world,
            self.angular_momentum_body,
        )
        .map(|_| ())
    }
}

/// Constant external wrench over one step, in declared world/body frames.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Wrench {
    /// Applied force in newtons, expressed in the world frame.
    pub force_world: Vec3,
    /// Applied torque in newton metres, expressed in the principal body frame.
    pub torque_body: Vec3,
}

impl Wrench {
    /// The zero wrench.
    pub const ZERO: Self = Self {
        force_world: Vec3::ZERO,
        torque_body: Vec3::ZERO,
    };

    fn validate(self) -> Result<(), DynamicsError> {
        if !self.force_world.is_finite() {
            return Err(DynamicsError::NonFinite("wrench.force_world"));
        }
        if !self.torque_body.is_finite() {
            return Err(DynamicsError::NonFinite("wrench.torque_body"));
        }
        Ok(())
    }
}

/// A uniform gravitational acceleration in the world frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Gravity {
    /// Acceleration in m/s², expressed in the world frame.
    acceleration_world: Vec3,
}

impl Gravity {
    /// No gravitational acceleration.
    pub const ZERO: Self = Self {
        acceleration_world: Vec3::ZERO,
    };

    /// Validates a gravitational acceleration.
    pub fn new(acceleration_world: Vec3) -> Result<Self, DynamicsError> {
        if !acceleration_world.is_finite() {
            return Err(DynamicsError::NonFinite("gravity.acceleration_world"));
        }
        Ok(Self { acceleration_world })
    }

    /// Returns world-frame gravitational acceleration.
    #[must_use]
    pub const fn acceleration_world(self) -> Vec3 {
        self.acceleration_world
    }

    fn validate(self) -> Result<(), DynamicsError> {
        Self::new(self.acceleration_world).map(|_| ())
    }
}

/// Conservation and balance diagnostics evaluated at one state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DynamicsDiagnostics {
    /// Translational kinetic energy in joules.
    pub translational_kinetic_energy: f64,
    /// Rotational kinetic energy in joules.
    pub rotational_kinetic_energy: f64,
    /// Potential energy of uniform gravity, using zero at the world origin.
    pub gravitational_potential_energy: f64,
    /// Total mechanical energy in joules.
    pub mechanical_energy: f64,
    /// Linear momentum in the world frame.
    pub linear_momentum_world: Vec3,
    /// Total angular momentum about the world origin, in the world frame.
    pub angular_momentum_world: Vec3,
}

/// A completed fixed step and its before/after diagnostics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StepReceipt {
    /// Fixed duration in seconds.
    pub duration_seconds: f64,
    /// State before the step.
    pub state_before: RigidBodyState,
    /// State after the step.
    pub state_after: RigidBodyState,
    /// Diagnostics before the step.
    pub diagnostics_before: DynamicsDiagnostics,
    /// Diagnostics after the step.
    pub diagnostics_after: DynamicsDiagnostics,
}

/// Result of advancing a sequence of complete fixed steps.
#[derive(Clone, Debug, PartialEq)]
pub enum AdvanceOutcome {
    /// Every requested step was completed.
    Completed {
        /// Number of completed steps.
        completed_steps: usize,
        /// Diagnostics after the final state.
        final_diagnostics: DynamicsDiagnostics,
    },
    /// Cancellation was observed before the next whole step began.
    Cancelled {
        /// Number of complete steps committed before cancellation.
        completed_steps: usize,
        /// Diagnostics of the last fully committed state.
        final_diagnostics: DynamicsDiagnostics,
    },
}

/// Deterministic fixed-step integrator for smooth, unconstrained rigid bodies.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidBodyIntegrator {
    gravity: Gravity,
}

impl RigidBodyIntegrator {
    /// Creates an integrator with the specified uniform gravity.
    #[must_use]
    pub const fn new(gravity: Gravity) -> Self {
        Self { gravity }
    }

    /// Returns the configured gravity.
    #[must_use]
    pub const fn gravity(self) -> Gravity {
        self.gravity
    }

    /// Evaluates energy and world-frame momentum diagnostics.
    #[must_use]
    pub fn diagnostics(
        self,
        state: RigidBodyState,
        properties: MassProperties,
    ) -> Result<DynamicsDiagnostics, DynamicsError> {
        self.gravity.validate()?;
        state.validate()?;
        properties.validate()?;
        let translational_kinetic_energy = finite_derived(
            state.linear_momentum_world.norm_squared() / (2.0 * properties.mass),
            "diagnostics.translational_kinetic_energy",
        )?;
        let angular_velocity_body = properties.angular_velocity_body(state.angular_momentum_body);
        if !angular_velocity_body.is_finite() {
            return Err(DynamicsError::NonFinite(
                "diagnostics.angular_velocity_body",
            ));
        }
        let rotational_kinetic_energy = finite_derived(
            0.5 * state.angular_momentum_body.dot(angular_velocity_body),
            "diagnostics.rotational_kinetic_energy",
        )?;
        let gravitational_potential_energy = finite_derived(
            -properties.mass
                * self
                    .gravity
                    .acceleration_world
                    .dot(state.pose.position_world),
            "diagnostics.gravitational_potential_energy",
        )?;
        let linear_angular_momentum_world = state
            .pose
            .orientation
            .rotate_body_to_world(state.angular_momentum_body);
        let orbital_angular_momentum_world =
            state.pose.position_world.cross(state.linear_momentum_world);
        let angular_momentum_world =
            orbital_angular_momentum_world.add(linear_angular_momentum_world);
        if !angular_momentum_world.is_finite() {
            return Err(DynamicsError::NonFinite(
                "diagnostics.angular_momentum_world",
            ));
        }
        Ok(DynamicsDiagnostics {
            translational_kinetic_energy,
            rotational_kinetic_energy,
            gravitational_potential_energy,
            mechanical_energy: finite_derived(
                translational_kinetic_energy
                    + rotational_kinetic_energy
                    + gravitational_potential_energy,
                "diagnostics.mechanical_energy",
            )?,
            linear_momentum_world: state.linear_momentum_world,
            angular_momentum_world,
        })
    }

    /// Performs one deterministic midpoint Lie-group step.
    ///
    /// Force and torque are constant over this step. No contact, constraint, or
    /// collision query is performed. A failed validation leaves no mutable state
    /// because this method returns a new state in its receipt.
    pub fn step(
        self,
        state: RigidBodyState,
        properties: MassProperties,
        wrench: Wrench,
        duration_seconds: f64,
    ) -> Result<StepReceipt, DynamicsError> {
        if !duration_seconds.is_finite() || duration_seconds <= 0.0 {
            return Err(DynamicsError::InvalidStepDuration);
        }
        self.gravity.validate()?;
        state.validate()?;
        properties.validate()?;
        wrench.validate()?;
        let diagnostics_before = self.diagnostics(state, properties)?;
        let total_force_world = wrench
            .force_world
            .add(self.gravity.acceleration_world().scale(properties.mass));
        let linear_momentum_mid = state
            .linear_momentum_world
            .add(total_force_world.scale(0.5 * duration_seconds));
        let position_after = state
            .pose
            .position_world
            .add(linear_momentum_mid.scale(duration_seconds / properties.mass));
        let linear_momentum_after = state
            .linear_momentum_world
            .add(total_force_world.scale(duration_seconds));

        let angular_rate_before =
            angular_momentum_rate(state.angular_momentum_body, properties, wrench.torque_body);
        let angular_momentum_mid = state
            .angular_momentum_body
            .add(angular_rate_before.scale(0.5 * duration_seconds));
        let angular_rate_mid =
            angular_momentum_rate(angular_momentum_mid, properties, wrench.torque_body);
        let angular_momentum_after = state
            .angular_momentum_body
            .add(angular_rate_mid.scale(duration_seconds));
        let angular_velocity_mid = properties.angular_velocity_body(angular_momentum_mid);
        let orientation_after = state
            .pose
            .orientation
            .right_exp(angular_velocity_mid.scale(duration_seconds))?;

        let state_after = RigidBodyState::new(
            Pose::new(position_after, orientation_after)?,
            linear_momentum_after,
            angular_momentum_after,
        )?;
        let diagnostics_after = self.diagnostics(state_after, properties)?;
        Ok(StepReceipt {
            duration_seconds,
            state_before: state,
            state_after,
            diagnostics_before,
            diagnostics_after,
        })
    }

    /// Advances at most `step_count` fixed steps and observes cancellation only
    /// at whole-step boundaries. A cancellation result therefore never contains
    /// a partially applied force, torque, orientation, or diagnostics update.
    pub fn advance<F>(
        self,
        state: &mut RigidBodyState,
        properties: MassProperties,
        wrench: Wrench,
        duration_seconds: f64,
        step_count: usize,
        mut is_cancelled: F,
    ) -> Result<AdvanceOutcome, DynamicsError>
    where
        F: FnMut(usize) -> bool,
    {
        for completed_steps in 0..step_count {
            if is_cancelled(completed_steps) {
                return Ok(AdvanceOutcome::Cancelled {
                    completed_steps,
                    final_diagnostics: self.diagnostics(*state, properties)?,
                });
            }
            *state = self
                .step(*state, properties, wrench, duration_seconds)?
                .state_after;
        }
        Ok(AdvanceOutcome::Completed {
            completed_steps: step_count,
            final_diagnostics: self.diagnostics(*state, properties)?,
        })
    }
}

fn finite_derived(value: f64, field: &'static str) -> Result<f64, DynamicsError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DynamicsError::NonFinite(field))
    }
}

fn angular_momentum_rate(
    angular_momentum_body: Vec3,
    properties: MassProperties,
    torque_body: Vec3,
) -> Vec3 {
    // Euler equation in the rotating body frame: Ldot = L × omega + tau.
    angular_momentum_body
        .cross(properties.angular_velocity_body(angular_momentum_body))
        .add(torque_body)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPSILON: f64 = 1e-10;

    fn assert_close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "expected {expected:.16e}, got {actual:.16e}, tolerance {tolerance:.16e}"
        );
    }

    fn disc_properties() -> MassProperties {
        // A generic axisymmetric disc in its center-of-mass principal frame.
        MassProperties::new(2.0, Vec3::ZERO, Vec3::new(0.25, 0.25, 0.5)).unwrap()
    }

    fn state() -> RigidBodyState {
        RigidBodyState::new(Pose::identity(), Vec3::ZERO, Vec3::ZERO).unwrap()
    }

    #[test]
    fn rejects_nonphysical_mass_and_inertia() {
        assert_eq!(
            MassProperties::new(0.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)),
            Err(DynamicsError::InvalidMass)
        );
        assert_eq!(
            MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, -1.0)),
            Err(DynamicsError::InvalidPrincipalInertia)
        );
        assert_eq!(
            MassProperties::new(1.0, Vec3::ZERO, Vec3::new(3.0, 1.0, 1.0)),
            Err(DynamicsError::InconsistentPrincipalInertia)
        );
        assert_eq!(
            MassProperties::new(1.0, Vec3::new(0.1, 0.0, 0.0), Vec3::new(1.0, 1.0, 1.0)),
            Err(DynamicsError::UnsupportedReferenceOffset)
        );
    }

    #[test]
    fn quaternion_sign_is_canonical_and_rotation_is_correct() {
        let positive = UnitQuaternion::new(0.0, 0.0, 0.0, 1.0).unwrap();
        let negative = UnitQuaternion::new(0.0, 0.0, 0.0, -2.0).unwrap();
        assert_eq!(positive, negative);
        let rotated = positive.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
        assert_close(rotated.x, -1.0, EPSILON);
        assert_close(rotated.y, 0.0, EPSILON);
    }

    #[test]
    fn right_handed_z_rotation_and_composition_match_known_oracles() {
        let quarter_turn =
            UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), core::f64::consts::FRAC_PI_2)
                .unwrap();
        let rotated = quarter_turn.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
        assert_close(rotated.x, 0.0, EPSILON);
        assert_close(rotated.y, 1.0, EPSILON);
        assert_close(rotated.z, 0.0, EPSILON);

        let half_turn = quarter_turn
            .right_exp(Vec3::new(0.0, 0.0, core::f64::consts::FRAC_PI_2))
            .unwrap();
        let composed = half_turn.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
        assert_close(composed.x, -1.0, EPSILON);
        assert_close(composed.y, 0.0, EPSILON);
        assert_close(composed.z, 0.0, EPSILON);
    }

    #[test]
    fn scaled_normalization_admits_extreme_finite_quaternion_axis_and_rotation_inputs() {
        let huge = UnitQuaternion::new(f64::MAX, f64::MAX, -f64::MAX, f64::MAX).unwrap();
        for component in huge.components() {
            assert!(component.is_finite());
        }
        let tiny = f64::from_bits(1);
        let tiny_quaternion = UnitQuaternion::new(tiny, -tiny, tiny, -tiny).unwrap();
        for component in tiny_quaternion.components() {
            assert!(component.is_finite());
        }
        let huge_axis = UnitQuaternion::from_axis_angle(
            Vec3::new(f64::MAX, 0.5 * f64::MAX, -0.25 * f64::MAX),
            core::f64::consts::FRAC_PI_2,
        )
        .unwrap();
        assert!(
            huge_axis
                .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0))
                .is_finite()
        );
        let huge_rotation = UnitQuaternion::IDENTITY
            .right_exp(Vec3::new(f64::MAX, f64::MAX, f64::MAX))
            .unwrap();
        for component in huge_rotation.components() {
            assert!(component.is_finite());
        }
    }

    #[test]
    fn gravity_matches_the_constant_acceleration_solution() {
        let gravity = Gravity::new(Vec3::new(0.0, 0.0, -9.81)).unwrap();
        let integrator = RigidBodyIntegrator::new(gravity);
        let receipt = integrator
            .step(state(), disc_properties(), Wrench::ZERO, 0.25)
            .unwrap();
        assert_close(
            receipt.state_after.linear_momentum_world().z,
            -4.905,
            EPSILON,
        );
        assert_close(
            receipt.state_after.pose().position_world().z,
            -0.3065625,
            EPSILON,
        );
        assert_close(
            receipt.diagnostics_before.mechanical_energy,
            receipt.diagnostics_after.mechanical_energy,
            EPSILON,
        );
    }

    #[test]
    fn constant_force_and_torque_update_the_declared_momentum_frames() {
        let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
        let wrench = Wrench {
            force_world: Vec3::new(4.0, 0.0, 0.0),
            torque_body: Vec3::new(0.0, 0.0, 3.0),
        };
        let receipt = integrator
            .step(state(), disc_properties(), wrench, 0.5)
            .unwrap();
        assert_close(receipt.state_after.linear_momentum_world().x, 2.0, EPSILON);
        assert_close(receipt.state_after.pose().position_world().x, 0.25, EPSILON);
        assert_close(receipt.state_after.angular_momentum_body().z, 1.5, EPSILON);
    }

    #[test]
    fn torque_free_spherical_body_preserves_energy_and_angular_momentum() {
        let properties = MassProperties::new(3.0, Vec3::ZERO, Vec3::new(2.0, 2.0, 2.0)).unwrap();
        let initial = RigidBodyState::new(
            Pose::identity(),
            Vec3::new(1.0, -2.0, 3.0),
            Vec3::new(2.0, -3.0, 5.0),
        )
        .unwrap();
        let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
        let initial_diagnostics = integrator.diagnostics(initial, properties).unwrap();
        let mut current = initial;
        for _ in 0..1_000 {
            current = integrator
                .step(current, properties, Wrench::ZERO, 0.01)
                .unwrap()
                .state_after;
        }
        let final_diagnostics = integrator.diagnostics(current, properties).unwrap();
        assert_close(
            final_diagnostics.mechanical_energy,
            initial_diagnostics.mechanical_energy,
            2e-12,
        );
        assert_close(
            final_diagnostics.angular_momentum_world.x,
            initial_diagnostics.angular_momentum_world.x,
            2e-12,
        );
        assert_close(
            final_diagnostics.angular_momentum_world.y,
            initial_diagnostics.angular_momentum_world.y,
            2e-12,
        );
        assert_close(
            final_diagnostics.angular_momentum_world.z,
            initial_diagnostics.angular_momentum_world.z,
            2e-12,
        );
    }

    #[test]
    fn axisymmetric_euler_top_phase_has_the_right_sign_and_refines() {
        let properties = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(2.0, 2.0, 3.5)).unwrap();
        let initial_omega = Vec3::new(1.1, 0.0, 0.8);
        let initial = RigidBodyState::new(
            Pose::identity(),
            Vec3::ZERO,
            Vec3::new(
                2.0 * initial_omega.x,
                2.0 * initial_omega.y,
                3.5 * initial_omega.z,
            ),
        )
        .unwrap();
        let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
        let phase_rate = (3.5 - 2.0) * initial_omega.z / 2.0;
        let duration = 1.0;
        let expected = Vec3::new(
            initial_omega.x * (phase_rate * duration).cos(),
            initial_omega.x * (phase_rate * duration).sin(),
            initial_omega.z,
        );

        let mut coarse = initial;
        for _ in 0..100 {
            coarse = integrator
                .step(coarse, properties, Wrench::ZERO, 0.01)
                .unwrap()
                .state_after;
        }
        let mut fine = initial;
        for _ in 0..200 {
            fine = integrator
                .step(fine, properties, Wrench::ZERO, 0.005)
                .unwrap()
                .state_after;
        }
        let coarse_error = properties
            .angular_velocity_body(coarse.angular_momentum_body())
            .sub(expected)
            .stable_norm("coarse_error")
            .unwrap();
        let fine_omega = properties.angular_velocity_body(fine.angular_momentum_body());
        let fine_error = fine_omega.sub(expected).stable_norm("fine_error").unwrap();
        assert!(
            fine_omega.y > 0.0,
            "positive axial spin must advance the phase toward +y"
        );
        assert!(
            fine_error < coarse_error,
            "fine {fine_error:e} must improve coarse {coarse_error:e}"
        );
    }

    #[test]
    fn torque_free_refinement_reduces_energy_defect_for_an_asymmetric_body() {
        let properties = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.5, 2.0)).unwrap();
        let initial =
            RigidBodyState::new(Pose::identity(), Vec3::ZERO, Vec3::new(0.4, 0.7, 1.1)).unwrap();
        let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
        let initial_energy = integrator
            .diagnostics(initial, properties)
            .unwrap()
            .mechanical_energy;
        let mut coarse = initial;
        let mut fine = initial;
        for _ in 0..100 {
            coarse = integrator
                .step(coarse, properties, Wrench::ZERO, 0.01)
                .unwrap()
                .state_after;
        }
        for _ in 0..200 {
            fine = integrator
                .step(fine, properties, Wrench::ZERO, 0.005)
                .unwrap()
                .state_after;
        }
        let coarse_defect = (integrator
            .diagnostics(coarse, properties)
            .unwrap()
            .mechanical_energy
            - initial_energy)
            .abs();
        let fine_defect = (integrator
            .diagnostics(fine, properties)
            .unwrap()
            .mechanical_energy
            - initial_energy)
            .abs();
        assert!(
            fine_defect < coarse_defect,
            "fine {fine_defect:e} must improve coarse {coarse_defect:e}"
        );
    }

    #[test]
    fn cancellation_occurs_before_a_complete_step_and_keeps_prior_state() {
        let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
        let mut current = state();
        let outcome = integrator
            .advance(
                &mut current,
                disc_properties(),
                Wrench {
                    force_world: Vec3::new(2.0, 0.0, 0.0),
                    torque_body: Vec3::ZERO,
                },
                0.5,
                5,
                |completed_steps| completed_steps == 2,
            )
            .unwrap();
        assert!(matches!(
            outcome,
            AdvanceOutcome::Cancelled {
                completed_steps: 2,
                ..
            }
        ));
        assert_close(current.linear_momentum_world().x, 2.0, EPSILON);
        assert_close(current.pose().position_world().x, 0.5, EPSILON);
    }

    #[test]
    fn overflowing_diagnostics_and_cancelled_advance_refuse_without_a_false_receipt() {
        let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
        let state =
            RigidBodyState::new(Pose::identity(), Vec3::new(f64::MAX, 0.0, 0.0), Vec3::ZERO)
                .unwrap();
        assert_eq!(
            integrator.diagnostics(state, disc_properties()),
            Err(DynamicsError::NonFinite(
                "diagnostics.translational_kinetic_energy"
            ))
        );
        let mut cancelled_state = state;
        assert_eq!(
            integrator.advance(
                &mut cancelled_state,
                disc_properties(),
                Wrench::ZERO,
                0.1,
                1,
                |_| true,
            ),
            Err(DynamicsError::NonFinite(
                "diagnostics.translational_kinetic_energy"
            ))
        );
        assert_eq!(cancelled_state, state);
    }
}
