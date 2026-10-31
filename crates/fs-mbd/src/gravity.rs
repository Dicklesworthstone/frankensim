//! Conservative gravity, sampled at both ends of a native rigid-body step.
//!
//! Uniform acceleration and a fixed Newtonian source share force/potential
//! ownership. A separately selected quadrupole approximation adds both the
//! gravity-gradient torque AND its matching translational force. The latter
//! matters: adding torque alone violates angular momentum and energy balance.
//! All quantities use SI and one caller-identified inertial frame.

use crate::{
    DynamicsDiagnostics, DynamicsError, Gravity, MassProperties, Pose, RigidBodyIntegrator,
    RigidBodyState, StepReceipt, Vec3, Wrench,
};

/// Invalid field/frame, inadmissible finite-body approximation, or native failure.
#[derive(Clone, Debug, PartialEq)]
pub enum GravityError {
    /// Invalid or unrepresentable physical input/output.
    Invalid(&'static str),
    /// Query and source do not share the declared inertial frame.
    FrameMismatch,
    /// The caller's finite-body expansion limit was exceeded.
    OutsideExpansionDomain,
    /// The native rigid-body owner refused a state or arithmetic result.
    Dynamics(DynamicsError),
    /// No candidate state was published.
    Cancelled,
}

impl core::fmt::Display for GravityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "gravity: {self:?}")
    }
}
impl std::error::Error for GravityError {}
impl From<DynamicsError> for GravityError {
    fn from(error: DynamicsError) -> Self {
        Self::Dynamics(error)
    }
}

/// One conservative field evaluation on a body's actual mass and pose.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GravityLoad {
    /// Force in the world frame and torque about COM in the body frame.
    pub wrench: Wrench,
    /// Potential in joules, with the field's documented reference zero.
    pub potential_energy_j: f64,
}

impl GravityLoad {
    fn validate(self) -> Result<Self, GravityError> {
        self.wrench.validate()?;
        finite(self.potential_energy_j)?;
        Ok(self)
    }
}

/// A stationary conservative field in the identified inertial frame.
///
/// Implementations must return force/torque equal to the negative translational/
/// rotational derivatives of the SAME potential. This semantic obligation is
/// not established merely by implementing the trait. Velocity-dependent and
/// time-dependent forces require a different work owner.
pub trait RigidGravityField {
    /// Evaluate the actual COM pose and the original native mass properties.
    fn load(
        &self,
        frame: u64,
        pose: Pose,
        properties: MassProperties,
    ) -> Result<GravityLoad, GravityError>;
}

/// The existing uniform gravity value bound to a frame and potential datum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UniformGravity {
    frame: u64,
    gravity: Gravity,
    zero_potential_at_world_m: Vec3,
}

impl UniformGravity {
    /// Bind explicit acceleration and the point where specific potential is zero.
    pub fn new(
        frame: u64,
        gravity: Gravity,
        zero_potential_at_world_m: Vec3,
    ) -> Result<Self, GravityError> {
        gravity.validate()?;
        if !zero_potential_at_world_m.is_finite() {
            return Err(GravityError::Invalid(
                "uniform gravity datum must be finite",
            ));
        }
        Ok(Self {
            frame,
            gravity,
            zero_potential_at_world_m,
        })
    }
}

impl RigidGravityField for UniformGravity {
    fn load(
        &self,
        frame: u64,
        pose: Pose,
        properties: MassProperties,
    ) -> Result<GravityLoad, GravityError> {
        same_frame(frame, self.frame)?;
        properties.validate()?;
        let displacement = pose.position_world().sub(self.zero_potential_at_world_m);
        let acceleration = self.gravity.acceleration_world();
        GravityLoad {
            wrench: Wrench {
                force_world: acceleration.scale(properties.mass()),
                torque_body: Vec3::ZERO,
            },
            potential_energy_j: -properties.mass() * acceleration.dot(displacement),
        }
        .validate()
    }
}

/// A fixed Newtonian point source, with no inferred planet or softening length.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CentralGravity {
    frame: u64,
    center_world_m: Vec3,
    gravitational_parameter_m3_s2: f64,
}

impl CentralGravity {
    /// Bind an inertial source position and positive gravitational parameter G M.
    pub fn new(
        frame: u64,
        center_world_m: Vec3,
        gravitational_parameter_m3_s2: f64,
    ) -> Result<Self, GravityError> {
        if !center_world_m.is_finite()
            || !gravitational_parameter_m3_s2.is_finite()
            || gravitational_parameter_m3_s2 <= 0.0
        {
            return Err(GravityError::Invalid(
                "central gravity needs finite center and positive G M",
            ));
        }
        Ok(Self {
            frame,
            center_world_m,
            gravitational_parameter_m3_s2,
        })
    }

    /// Acceleration [m/s²] and specific potential [J/kg], zero at infinity.
    /// No radius powers or reciprocal direction normalization can overflow.
    pub fn sample(&self, frame: u64, position_world_m: Vec3) -> Result<(Vec3, f64), GravityError> {
        let (direction, radius, potential) = self.radial(frame, position_world_m)?;
        let acceleration = direction.scale(finite(potential / radius)?);
        if !acceleration.is_finite() {
            return Err(GravityError::Invalid(
                "unrepresentable central acceleration",
            ));
        }
        Ok((acceleration, potential))
    }

    fn radial(&self, frame: u64, position: Vec3) -> Result<(Vec3, f64, f64), GravityError> {
        same_frame(frame, self.frame)?;
        let r = position.sub(self.center_world_m);
        let radius = r.stable_norm("gravity radius")?;
        if radius <= 0.0 {
            return Err(GravityError::Invalid(
                "central gravity is singular at the source",
            ));
        }
        let direction = Vec3::new(r.x / radius, r.y / radius, r.z / radius);
        let potential = finite(-self.gravitational_parameter_m3_s2 / radius)?;
        Ok((direction, radius, potential))
    }

    /// Select a finite-body quadrupole model with a caller-supplied COM bounding
    /// radius and maximum size/distance ratio, strictly between zero and one.
    /// This explicit domain is not a certified multipole truncation-error bound.
    pub fn quadrupole(
        self,
        body_radius_m: f64,
        maximum_radius_ratio: f64,
    ) -> Result<QuadrupoleGravity, GravityError> {
        if !body_radius_m.is_finite()
            || body_radius_m <= 0.0
            || !maximum_radius_ratio.is_finite()
            || maximum_radius_ratio <= 0.0
            || maximum_radius_ratio >= 1.0
        {
            return Err(GravityError::Invalid(
                "quadrupole needs positive body radius and ratio in (0,1)",
            ));
        }
        Ok(QuadrupoleGravity {
            source: self,
            body_radius_m,
            maximum_radius_ratio,
        })
    }
}

impl RigidGravityField for CentralGravity {
    fn load(
        &self,
        frame: u64,
        pose: Pose,
        properties: MassProperties,
    ) -> Result<GravityLoad, GravityError> {
        properties.validate()?;
        let (acceleration, potential) = self.sample(frame, pose.position_world())?;
        GravityLoad {
            wrench: Wrench {
                force_world: acceleration.scale(properties.mass()),
                torque_body: Vec3::ZERO,
            },
            potential_energy_j: properties.mass() * potential,
        }
        .validate()
    }
}

/// Monopole plus second mass moment of an extended rigid body.
///
/// Uses the actual full native inertia tensor, including rotated principal axes.
/// With n from source to COM and I in world coordinates,
/// U = -mu m/r + mu(3 n.I.n - tr I)/(2 r³),
/// tau = 3 mu/r³ n cross I n. Force differentiates the complete U.
/// See NASA-CR-188243 (1993), gravity partials/gravity-gradient torque:
/// <https://ntrs.nasa.gov/citations/19940025085>.
/// No source recoil, self-gravity, J2, tides in a deformable body or relativity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuadrupoleGravity {
    source: CentralGravity,
    body_radius_m: f64,
    maximum_radius_ratio: f64,
}

impl RigidGravityField for QuadrupoleGravity {
    fn load(
        &self,
        frame: u64,
        pose: Pose,
        properties: MassProperties,
    ) -> Result<GravityLoad, GravityError> {
        properties.validate()?;
        let (n, radius, specific_potential) = self.source.radial(frame, pose.position_world())?;
        if self.body_radius_m / radius > self.maximum_radius_ratio {
            return Err(GravityError::OutsideExpansionDomain);
        }
        let principal = properties.principal_inertia_body();
        let trace = finite(principal.x + principal.y + principal.z)?;
        // Necessary (not sufficient) consistency of declared COM bounding sphere:
        // tr(I)/(2m) = mean squared mass radius <= bounding radius squared.
        if (0.5 * trace / properties.mass()).sqrt()
            > self.body_radius_m * (1.0 + 32.0 * f64::EPSILON)
        {
            return Err(GravityError::Invalid(
                "inertia exceeds declared body radius",
            ));
        }
        let n_body = pose.orientation().rotate_world_to_body(n);
        // This is the existing full inertia action; no per-query diagonalization.
        let inertia_n_body = properties.angular_momentum_body_checked(n_body)?;
        let n_inertia_n = n_body.dot(inertia_n_body);
        let mu_over_r3 = finite(-specific_potential / radius / radius)?;
        let correction_potential = finite(0.5 * mu_over_r3 * (3.0 * n_inertia_n - trace))?;
        let inertia_n_world = pose.orientation().rotate_body_to_world(inertia_n_body);
        let correction_force = n
            .scale(5.0 * n_inertia_n - trace)
            .sub(inertia_n_world.scale(2.0))
            .scale(1.5 * mu_over_r3 / radius);
        GravityLoad {
            wrench: Wrench {
                force_world: n
                    .scale((specific_potential / radius) * properties.mass())
                    .add(correction_force),
                torque_body: n_body.cross(inertia_n_body).scale(3.0 * mu_over_r3),
            },
            potential_energy_j: properties.mass() * specific_potential + correction_potential,
        }
        .validate()
    }
}

/// Second-order kick/drift/kick under a spatial gravity field.
///
/// Uses the original native `step_verlet_with` map and mass properties.
/// Field force AND torque are resampled at the actual endpoint pose. Translation
/// is velocity Verlet; potential replaces the uniform-origin diagnostic. With
/// no other wrench, conservative fields and fixed steps give a symplectic map with a
/// timestep-dependent energy defect, not exact energy conservation. Bounded
/// long-time energy behavior requires resolved, smooth, bounded dynamics.
/// Other force/torque is held over this interval and excluded from stored energy.
/// Cancellation/error returns no partly advanced state. Source callbacks must
/// be repeatable and side-effect-free; a frame ID checks identity, not inertiality.
/// Endpoint admission does not certify the intervening path: close approaches,
/// source crossings and finite-body domain crossings require step refinement
/// and a caller-owned event/trajectory admission policy.
pub fn step_in_field(
    field: &(impl RigidGravityField + ?Sized),
    frame: u64,
    state: RigidBodyState,
    properties: MassProperties,
    other_wrench: Wrench,
    duration_seconds: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<StepReceipt, GravityError> {
    if !duration_seconds.is_finite() || duration_seconds <= 0.0 {
        return Err(DynamicsError::InvalidStepDuration.into());
    }
    state.validate()?;
    properties.validate()?;
    other_wrench.validate()?;
    poll(&mut cancelled)?;
    let before = field.load(frame, state.pose(), properties)?.validate()?;
    poll(&mut cancelled)?;
    let (mut receipt, potential_after) = RigidBodyIntegrator::new(Gravity::ZERO)
        .step_verlet_with::<_, GravityError>(
            state,
            properties,
            combine(before.wrench, other_wrench),
            duration_seconds,
            |drifted| {
                poll(&mut cancelled)?;
                let after = field.load(frame, drifted.pose(), properties)?.validate()?;
                Ok((combine(after.wrench, other_wrench), after.potential_energy_j))
            },
        )?;
    receipt.diagnostics_before =
        with_potential(receipt.diagnostics_before, before.potential_energy_j)?;
    receipt.diagnostics_after = with_potential(receipt.diagnostics_after, potential_after)?;
    poll(&mut cancelled)?;
    Ok(receipt)
}

fn combine(gravity: Wrench, other: Wrench) -> Wrench {
    Wrench {
        force_world: gravity.force_world.add(other.force_world),
        torque_body: gravity.torque_body.add(other.torque_body),
    }
}

fn with_potential(
    mut result: DynamicsDiagnostics,
    potential: f64,
) -> Result<DynamicsDiagnostics, GravityError> {
    result.gravitational_potential_energy = finite(potential)?;
    result.mechanical_energy =
        finite(result.translational_kinetic_energy + result.rotational_kinetic_energy + potential)?;
    Ok(result)
}

fn same_frame(actual: u64, expected: u64) -> Result<(), GravityError> {
    if actual == expected {
        Ok(())
    } else {
        Err(GravityError::FrameMismatch)
    }
}
fn finite(value: f64) -> Result<f64, GravityError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(GravityError::Invalid("nonfinite gravity arithmetic"))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), GravityError> {
    if cancelled() {
        Err(GravityError::Cancelled)
    } else {
        Ok(())
    }
}
