//! Conservative gravity, sampled at both ends of a native rigid-body step.
//!
//! Uniform acceleration and a fixed Newtonian source share force/potential
//! ownership. A separately selected quadrupole approximation adds both the
//! gravity-gradient torque AND its matching translational force. The latter
//! matters: adding torque alone violates angular momentum and energy balance.
//! All quantities use SI and one caller-identified inertial frame.

/// Conservative point-field values and analytic spatial derivatives.
pub mod field;
pub use field::{
    PointGravityField, PointGravitySample, QuadrupoleSourceExtent, QuadrupoleSourceGravity,
    UniformSphereGravity,
};

use crate::{
    DynamicsDiagnostics, DynamicsError, Gravity, MassProperties, Pose, RigidBodyIntegrator,
    RigidBodyState, StepReceipt, Vec3, Wrench,
};
use fs_math::dd::Dd;

/// Closed, recoiling systems with explicitly bounded direct Newtonian interactions.
pub mod mutual;

mod superposition;
pub use superposition::SuperposedGravity;

mod adaptive;
mod adaptive_field;
pub use adaptive::{
    GravityAdvanceError, GravityAdvanceOptions, GravityAdvanceReceipt, GravityAdvanceStats,
};
pub use adaptive_field::advance_in_field;

/// Invalid field/frame, inadmissible finite-body approximation, or native failure.
#[derive(Clone, Debug, PartialEq)]
pub enum GravityError {
    /// Invalid or unrepresentable physical input/output.
    Invalid(&'static str),
    /// Query and source do not share the declared inertial frame.
    FrameMismatch,
    /// The caller's finite-body expansion limit was exceeded.
    OutsideExpansionDomain,
    /// A pair endpoint or straight split drift reached its declared cutoff.
    /// Adaptive owners may refine a trial only after admitting the initial state.
    SeparationCutoff {
        /// Rejected endpoint or closest straight-drift separation [m].
        separation_m: f64,
        /// Required strict lower bound on separation [m].
        minimum_separation_m: f64,
    },
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

    /// Potential at `after` minus that at `before`, in joules.
    ///
    /// A direct difference can retain work smaller than the rounding of either
    /// absolute endpoint energy. It must use the same conservative potential,
    /// mass and frame as [`Self::load`]. The default leaves callers' retained
    /// endpoint loads as their fallback without sampling the field again.
    fn potential_energy_increment(
        &self,
        _frame: u64,
        _before: Pose,
        _after: Pose,
        _properties: MassProperties,
    ) -> Result<Option<f64>, GravityError> {
        Ok(None)
    }
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

    fn specific_potential_change(
        &self,
        frame: u64,
        before: Vec3,
        after: Vec3,
    ) -> Result<Dd, GravityError> {
        same_frame(frame, self.frame)?;
        if !before.is_finite() || !after.is_finite() {
            return Err(GravityError::Invalid("nonfinite gravity endpoint"));
        }
        let acceleration = self.gravity.acceleration_world();
        let mut change = Dd::ZERO;
        for (g, a, b) in [
            (acceleration.x, after.x, before.x),
            (acceleration.y, after.y, before.y),
            (acceleration.z, after.z, before.z),
        ] {
            // Cancel the datum symbolically, before either absolute potential
            // or the endpoint subtraction can round away the physical work.
            change = change - Dd::from_f64(g) * (Dd::from_f64(a) - Dd::from_f64(b));
        }
        Ok(change)
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

    fn potential_energy_increment(
        &self,
        frame: u64,
        before: Pose,
        after: Pose,
        properties: MassProperties,
    ) -> Result<Option<f64>, GravityError> {
        properties.validate()?;
        let change =
            self.specific_potential_change(frame, before.position_world(), after.position_world())?;
        Ok(Some(finite(
            (Dd::from_f64(properties.mass()) * change).to_f64(),
        )?))
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

    fn specific_potential_change(
        &self,
        frame: u64,
        before: Vec3,
        after: Vec3,
    ) -> Result<Dd, GravityError> {
        self.specific_potential_change_to(self, frame, before, after)
    }

    /// Specific interaction-potential change between two source snapshots [J/kg].
    ///
    /// Both the source COM and sampled material point may move. The snapshots
    /// must retain the same G M and inertial frame; this is not a variable-mass
    /// source law. The original radial domains admit their own endpoints.
    /// Direct relative coordinate increments retain source recoil and work below
    /// the rounding of either absolute potential. No force-work quadrature,
    /// source advance, reaction force or path-clearance guarantee is supplied.
    pub fn specific_potential_increment_to(
        &self,
        source_after: &Self,
        frame: u64,
        before_world_m: Vec3,
        after_world_m: Vec3,
    ) -> Result<f64, GravityError> {
        finite(
            self.specific_potential_change_to(source_after, frame, before_world_m, after_world_m)?
                .to_f64(),
        )
    }

    fn specific_potential_change_to(
        &self,
        source_after: &Self,
        frame: u64,
        before: Vec3,
        after: Vec3,
    ) -> Result<Dd, GravityError> {
        if self.gravitational_parameter_m3_s2 != source_after.gravitational_parameter_m3_s2 {
            return Err(GravityError::Invalid(
                "moving central source must retain the original G M",
            ));
        }
        let (_, radius_before, potential_before) = self.radial(frame, before)?;
        let (_, radius_after, potential_after) = source_after.radial(frame, after)?;
        if radius_before < 0.5 * radius_after || radius_after < 0.5 * radius_before {
            // Well-separated radii have no cancellation. This branch also
            // avoids an unrepresentable ratio between two admitted radii.
            return Ok(Dd::from_f64(potential_after) - Dd::from_f64(potential_before));
        }
        Ok(inverse_radius_potential_increment(
            [self.center_world_m, before],
            [source_after.center_world_m, after],
            [radius_before, radius_after],
            -potential_before,
        ))
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

// Comparable admitted radii only: rationalize U_after - U_before for U = -C/r.
// Both original positions participate before rounding the relative displacement,
// so this same owner retains fixed-source work and finite-source recoil.
fn inverse_radius_potential_increment(
    before: [Vec3; 2],
    after: [Vec3; 2],
    radii: [f64; 2],
    magnitude_before: f64,
) -> Dd {
    let scale = radii[0].max(radii[1]);
    let divisor = Dd::from_f64(scale);
    let mut radius_squared_change = Dd::ZERO;
    for (a0, a1, b0, b1) in [
        (after[0].x, after[1].x, before[0].x, before[1].x),
        (after[0].y, after[1].y, before[0].y, before[1].y),
        (after[0].z, after[1].z, before[0].z, before[1].z),
    ] {
        let [a0, a1, b0, b1] = [a0, a1, b0, b1].map(Dd::from_f64);
        let a_scaled = (a1 - a0) / divisor;
        let b_scaled = (b1 - b0) / divisor;
        let difference = (a1 - b1) - (a0 - b0);
        let difference = if difference.hi.is_finite() {
            difference / divisor
        } else {
            // Common translation can overflow each body's unscaled travel
            // even while both relative endpoint separations remain admitted.
            a_scaled - b_scaled
        };
        radius_squared_change = radius_squared_change + difference * (a_scaled + b_scaled);
    }
    // Use coordinate increments, not the difference of rounded norms. A tiny
    // transverse displacement can change potential while those norms are equal.
    let before_scaled = Dd::from_f64(radii[0] / scale);
    let after_scaled = Dd::from_f64(radii[1] / scale);
    let relative_change = radius_squared_change / ((before_scaled + after_scaled) * after_scaled);
    Dd::from_f64(magnitude_before) * relative_change
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

    fn potential_energy_increment(
        &self,
        frame: u64,
        before: Pose,
        after: Pose,
        properties: MassProperties,
    ) -> Result<Option<f64>, GravityError> {
        properties.validate()?;
        let change =
            self.specific_potential_change(frame, before.position_world(), after.position_world())?;
        Ok(Some(finite(
            (Dd::from_f64(properties.mass()) * change).to_f64(),
        )?))
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
        let trace = quadrupole_inertia_trace(properties, self.body_radius_m)?;
        let mu_over_r3 = finite(-specific_potential / radius / radius)?;
        let correction = quadrupole_correction(pose, properties, n, radius, mu_over_r3, trace)?;
        GravityLoad {
            wrench: Wrench {
                force_world: n
                    .scale((specific_potential / radius) * properties.mass())
                    .add(correction.wrench.force_world),
                torque_body: correction.wrench.torque_body,
            },
            potential_energy_j: properties.mass() * specific_potential
                + correction.potential_energy_j,
        }
        .validate()
    }

    fn potential_energy_increment(
        &self,
        frame: u64,
        before: Pose,
        after: Pose,
        properties: MassProperties,
    ) -> Result<Option<f64>, GravityError> {
        properties.validate()?;
        let trace = quadrupole_inertia_trace(properties, self.body_radius_m)?;
        let mut correction = [0.0; 2];
        for (value, pose) in correction.iter_mut().zip([before, after]) {
            let (n, radius, potential) = self.source.radial(frame, pose.position_world())?;
            if self.body_radius_m / radius > self.maximum_radius_ratio {
                return Err(GravityError::OutsideExpansionDomain);
            }
            let mu_over_r3 = finite(-potential / radius / radius)?;
            *value = quadrupole_correction(pose, properties, n, radius, mu_over_r3, trace)?
                .potential_energy_j;
        }
        let central = self.source.specific_potential_change(
            frame,
            before.position_world(),
            after.position_world(),
        )?;
        // Do not add a tiny tidal potential to the much larger monopole before
        // taking the difference. Rotation-only tidal work can then survive
        // even when both complete endpoint energies round to the same value.
        let change = Dd::from_f64(properties.mass()) * central
            + (Dd::from_f64(correction[1]) - Dd::from_f64(correction[0]));
        Ok(Some(finite(change.to_f64())?))
    }
}

// Shared second-mass-moment owner for fixed-source and recoiling pair gravity.
// The COM radius is a caller declaration; this necessary inertia check does not
// reconstruct a body's shape or certify its complete support.
fn quadrupole_inertia_trace(
    properties: MassProperties,
    body_radius_m: f64,
) -> Result<f64, GravityError> {
    let principal = properties.principal_inertia_body();
    let trace = finite(principal.x + principal.y + principal.z)?;
    if (0.5 * trace / properties.mass()).sqrt() > body_radius_m * (1.0 + 32.0 * f64::EPSILON) {
        return Err(GravityError::Invalid(
            "inertia exceeds declared body radius",
        ));
    }
    Ok(trace)
}

// n points from the other COM/source to this body. This returns only the
// quadrupole contribution, so a mutual pair counts its monopole exactly once.
fn quadrupole_correction(
    pose: Pose,
    properties: MassProperties,
    n: Vec3,
    radius: f64,
    mu_over_r3: f64,
    trace: f64,
) -> Result<GravityLoad, GravityError> {
    let n_body = pose.orientation().rotate_world_to_body(n);
    let inertia_n_body = properties.angular_momentum_body_checked(n_body)?;
    let n_inertia_n = n_body.dot(inertia_n_body);
    let correction_potential = finite(0.5 * mu_over_r3 * (3.0 * n_inertia_n - trace))?;
    let inertia_n_world = pose.orientation().rotate_body_to_world(inertia_n_body);
    let correction_force = n
        .scale(5.0 * n_inertia_n - trace)
        .sub(inertia_n_world.scale(2.0))
        .scale(1.5 * mu_over_r3 / radius);
    GravityLoad {
        wrench: Wrench {
            force_world: correction_force,
            torque_body: n_body.cross(inertia_n_body).scale(3.0 * mu_over_r3),
        },
        potential_energy_j: correction_potential,
    }
    .validate()
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
    cancelled: impl FnMut() -> bool,
) -> Result<StepReceipt, GravityError> {
    Ok(step_in_field_with_loads(
        field,
        frame,
        state,
        properties,
        other_wrench,
        duration_seconds,
        cancelled,
    )?
    .step)
}

/// The original native step and the actual conservative samples that drove it.
///
/// Coupled consumers can account for external impulses, potential decrease and
/// numerical work defects without evaluating the field again or reconstructing
/// an intermediate momentum from a rounded endpoint. No sample certifies the
/// intervening trajectory's field domain.
#[derive(Clone, Copy, Debug)]
pub struct GravityStepReceipt {
    /// Completed native state and diagnostics, including this field's potential.
    pub step: StepReceipt,
    /// Field load at the initial pose, before the first half-kick.
    pub load_before: GravityLoad,
    /// Field load at the actual drift endpoint, before the final half-kick.
    pub load_after: GravityLoad,
    /// Field-owned endpoint potential change, in joules. Built-in uniform and
    /// central fields retain small changes hidden by absolute-energy rounding;
    /// custom fields without a direct increment use the retained load values.
    pub potential_energy_change_j: f64,
    /// Native state after free drift and before the endpoint half-kick.
    pub state_after_drift: RigidBodyState,
}

/// Advance the same map as [`step_in_field`], retaining its original stages.
///
/// The returned loads exclude `other_wrench`; the native step includes both.
/// Field, frame, full-inertia, cancellation and endpoint-domain semantics are
/// exactly those of [`step_in_field`].
pub fn step_in_field_with_loads(
    field: &(impl RigidGravityField + ?Sized),
    frame: u64,
    state: RigidBodyState,
    properties: MassProperties,
    other_wrench: Wrench,
    duration_seconds: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GravityStepReceipt, GravityError> {
    if !duration_seconds.is_finite() || duration_seconds <= 0.0 {
        return Err(DynamicsError::InvalidStepDuration.into());
    }
    state.validate()?;
    properties.validate()?;
    other_wrench.validate()?;
    poll(&mut cancelled)?;
    let before = field.load(frame, state.pose(), properties)?.validate()?;
    poll(&mut cancelled)?;
    let (mut receipt, (after, state_after_drift)) = RigidBodyIntegrator::new(Gravity::ZERO)
        .step_verlet_with::<_, GravityError>(
        state,
        properties,
        combine(before.wrench, other_wrench),
        duration_seconds,
        |drifted| {
            poll(&mut cancelled)?;
            let after = field.load(frame, drifted.pose(), properties)?.validate()?;
            Ok((combine(after.wrench, other_wrench), (after, drifted)))
        },
    )?;
    receipt.diagnostics_before =
        with_potential(receipt.diagnostics_before, before.potential_energy_j)?;
    receipt.diagnostics_after =
        with_potential(receipt.diagnostics_after, after.potential_energy_j)?;
    let potential_energy_change_j = finite(
        field
            .potential_energy_increment(frame, state.pose(), state_after_drift.pose(), properties)?
            .unwrap_or_else(|| {
                (Dd::from_f64(after.potential_energy_j) - Dd::from_f64(before.potential_energy_j))
                    .to_f64()
            }),
    )?;
    poll(&mut cancelled)?;
    Ok(GravityStepReceipt {
        step: receipt,
        load_before: before,
        load_after: after,
        potential_energy_change_j,
        state_after_drift,
    })
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
