//! Exact held isotropic linear drag on the original rigid-body state.
//!
//! Translation solves `m v' = F + m g - c (v - U)` and spherical free spin
//! solves `I omega' = -c_rot omega`. The fixed spin axis makes its integrated
//! rotation an exact group increment. Displacement, relative-slip heat and
//! maintaining-flow work are integrated along the same exponential path;
//! they are not endpoint trapezoids or energy-balance remainders.
//!
//! This is a native mechanics owner, not a fluid/material constitutive law.
//! Its caller supplies admitted resistances and the physical regime. The force
//! and uniform bath velocity may be held or explicitly affine in time;
//! properties and resistances stay fixed over the interval. An accelerating
//! fluid's pressure-gradient, added-mass and history forces are not inferred.
//! Anisotropic inertia, applied couples, rotating baths, lift, contact and
//! changing sources require their existing more general owners.

use crate::{
    DynamicsError, Gravity, MassProperties, Pose, RigidBodyState, RigidKineticEnergyChange, Vec3,
    Wrench,
};
use fs_math::{dd::Dd, decay};

mod affine;

/// An exact-flow input or arithmetic refusal; no state has been modified.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LinearDragError {
    /// An unsupported physical assumption or unrepresentable derived value.
    Invalid(&'static str),
    /// The original state, inertia or rotation owner refused an input.
    Dynamics(DynamicsError),
}

impl core::fmt::Display for LinearDragError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(message) => write!(f, "isotropic linear drag: {message}"),
            Self::Dynamics(error) => write!(f, "isotropic linear drag: {error}"),
        }
    }
}
impl std::error::Error for LinearDragError {}
impl From<DynamicsError> for LinearDragError {
    fn from(error: DynamicsError) -> Self {
        Self::Dynamics(error)
    }
}

/// A uniform damping bath acting on one body with exactly isotropic inertia.
/// Its prepared velocity is the starting endpoint for an affine-flow step.
///
/// Preparation is allocation-free and performs no eigensolve. Nonnegative
/// coefficients include the zero-drag limit; no small positive drag is erased.
/// Equal principal moments are required exactly, without inertia projection.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsotropicLinearDrag {
    properties: MassProperties,
    gravity: Gravity,
    ambient_velocity_world: Vec3,
    translation_resistance_kg_s: f64,
    rotary_resistance_kg_m2_s: f64,
}

/// Exact-model motion and independent source/storage accounts, rounded to f64.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearDragStep {
    /// Physical interval [s].
    pub duration_s: f64,
    /// Original authoritative state at the beginning of the interval.
    pub state_before: RigidBodyState,
    /// New authoritative pose and momenta; no velocity-only state is substituted.
    pub state_after: RigidBodyState,
    /// Integrated COM path [m], retained separately from rounded coordinates.
    pub displacement_world: Vec3,
    /// Integrated body angular velocity [rad], about its fixed spin axis.
    pub rotation_vector_body: Vec3,
    /// Integrated translational drag force on the body [N s].
    pub drag_impulse_world: Vec3,
    /// Integrated rotary drag couple on the body in the world frame [N m s].
    pub drag_angular_impulse_world: Vec3,
    /// Actual linear momentum increment minus applied, gravity and drag impulses.
    pub linear_momentum_residual_world: Vec3,
    /// Actual COM angular momentum increment minus rotary drag impulse [N m s].
    pub angular_momentum_residual_world: Vec3,
    /// Integral of `c |v-U|^2` [J], independent of the state-energy difference.
    pub translation_heat_j: f64,
    /// Integral of `c_rot |omega|^2` [J], independent of the energy residual.
    pub rotation_heat_j: f64,
    /// Total nonnegative physical dissipation [J]; no thermal destination assumed.
    pub dissipated_energy_j: f64,
    /// Signed work of the actuator maintaining the bath: integral of
    /// `F_drag . U(t)` [J]. For held U this is `drag_impulse . U`.
    pub background_work_j: f64,
    /// Integral of the prescribed applied force paired with COM velocity [J].
    pub applied_work_j: f64,
    /// Gravity paired with the integrated COM path [J].
    pub gravity_work_j: f64,
    /// Actual represented momentum-storage change, evaluated by the shared owner [J].
    pub kinetic_energy_change_j: f64,
    /// Actual kinetic plus gravity-potential change between represented states [J].
    pub mechanical_energy_change_j: f64,
    /// `delta(K+V) + heat - background_work - applied_work` [J].
    /// Includes final-state rounding and lost coordinate increments; never heat.
    pub work_energy_residual_j: f64,
}

impl IsotropicLinearDrag {
    /// Admit a body and fixed isotropic translation/rotation resistance.
    pub fn new(
        properties: MassProperties,
        gravity: Gravity,
        ambient_velocity_world: Vec3,
        translation_resistance_kg_s: f64,
        rotary_resistance_kg_m2_s: f64,
    ) -> Result<Self, LinearDragError> {
        properties.validate()?;
        gravity.validate()?;
        let inertia = properties.principal_inertia_body();
        if inertia.x != inertia.y || inertia.x != inertia.z {
            return Err(LinearDragError::Invalid(
                "exact free-spin flow requires isotropic inertia",
            ));
        }
        if !ambient_velocity_world.is_finite()
            || !translation_resistance_kg_s.is_finite()
            || translation_resistance_kg_s < 0.0
            || !rotary_resistance_kg_m2_s.is_finite()
            || rotary_resistance_kg_m2_s < 0.0
        {
            return Err(LinearDragError::Invalid(
                "finite bath velocity and nonnegative finite resistances required",
            ));
        }
        Ok(Self {
            properties,
            gravity,
            ambient_velocity_world,
            translation_resistance_kg_s,
            rotary_resistance_kg_m2_s,
        })
    }

    /// Advance the exact held model without a relaxation-number step restriction.
    ///
    /// The applied force is held in world coordinates and excludes gravity and
    /// drag. Applied torque must be zero: arbitrary forced spin changes its
    /// axis, so its attitude cannot use this exact single-axis exponential.
    /// This bounded allocation-free step returns a new state only on success.
    /// The caller still resolves source changes, contacts and validity limits.
    pub fn step(
        self,
        state: RigidBodyState,
        applied: Wrench,
        duration_s: f64,
    ) -> Result<LinearDragStep, LinearDragError> {
        state.validate()?;
        applied.validate()?;
        if applied.torque_body != Vec3::ZERO {
            return Err(LinearDragError::Invalid(
                "exact free-spin flow does not admit an applied torque",
            ));
        }
        if !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(DynamicsError::InvalidStepDuration.into());
        }
        let mass = dd(self.properties.mass());
        let dt = dd(duration_s);
        let translation = Moments::new(
            self.translation_resistance_kg_s,
            self.properties.mass(),
            duration_s,
        )?;
        let before_p = components(state.linear_momentum_world()).map(dd);
        let force = components(applied.force_world).map(dd);
        let weight = components(self.gravity.acceleration_world()).map(|g| mass * dd(g));
        let total_force: [Dd; 3] = core::array::from_fn(|i| force[i] + weight[i]);
        let ambient = components(self.ambient_velocity_world).map(dd);
        let mut after_p = [Dd::ZERO; 3];
        let mut displacement = [Dd::ZERO; 3];
        let mut drag_impulse = [Dd::ZERO; 3];
        let mut translation_heat = Dd::ZERO;
        for i in 0..3 {
            let velocity = before_p[i] / mass;
            let acceleration = total_force[i] / mass;
            if self.translation_resistance_kg_s == 0.0 {
                after_p[i] = before_p[i] + total_force[i] * dt;
                displacement[i] = dt * (velocity + acceleration * dt * dd(0.5));
                continue;
            }
            let relative_momentum = before_p[i] - mass * ambient[i];
            let relative_velocity = relative_momentum / mass;
            let velocity_increment =
                acceleration * translation.integral_s - relative_velocity * dd(translation.lost);
            let relative_mean =
                relative_velocity * translation.phi1 + acceleration * translation.drive_integral_s;
            // For short intervals preserve the initial absolute velocity,
            // rather than cancelling a large bath velocity against -U.
            let absolute_mean = if translation.exposure < 0.25 {
                velocity + acceleration * translation.drive_integral_s
                    - relative_velocity * dd(translation.exposure) * dd(translation.phi2)
            } else {
                ambient[i] + relative_mean
            };
            displacement[i] = dt * absolute_mean;
            after_p[i] = if translation.exposure < 0.25 {
                before_p[i] + total_force[i] * translation.integral_s
                    - relative_momentum * dd(translation.lost)
            } else {
                mass * ambient[i]
                    + attenuated(relative_momentum, translation.exposure)
                    + total_force[i] * translation.integral_s
            };
            drag_impulse[i] = -dd(self.translation_resistance_kg_s) * dt * relative_mean;
            // The mean-square + variance decomposition is nonnegative even
            // when drive opposes the initial slip. It integrates the path,
            // independently of rounded endpoint momenta or source work.
            translation_heat = translation_heat
                + dd(self.translation_resistance_kg_s)
                    * dt
                    * (relative_mean * relative_mean
                        + velocity_increment * velocity_increment * dd(translation.variance));
        }
        self.finish(
            state,
            duration_s,
            TranslationPath {
                after_p,
                displacement,
                drag_impulse,
                translation_heat,
                applied_work: dot(force, displacement),
                total_impulse: total_force.map(|value| value * dt),
            },
        )
    }

    // All prescribed source histories use the original native free-spin flow,
    // represented-state storage and retained impulse/work accounting below.
    fn finish(
        self,
        state: RigidBodyState,
        duration_s: f64,
        path: TranslationPath,
    ) -> Result<LinearDragStep, LinearDragError> {
        let ambient = components(self.ambient_velocity_world).map(dd);
        let background_work = dot(path.drag_impulse, ambient);
        self.finish_with_background_work(state, duration_s, path, background_work)
    }

    fn finish_with_background_work(
        self,
        state: RigidBodyState,
        duration_s: f64,
        path: TranslationPath,
        background_work: Dd,
    ) -> Result<LinearDragStep, LinearDragError> {
        let TranslationPath {
            after_p,
            displacement,
            drag_impulse,
            translation_heat,
            applied_work,
            total_impulse,
        } = path;
        let mass = dd(self.properties.mass());
        let inertia = dd(self.properties.principal_inertia_body().x);
        let rotation = Moments::new(self.rotary_resistance_kg_m2_s, inertia.hi, duration_s)?;
        let before_p = components(state.linear_momentum_world()).map(dd);
        let before_l = components(state.angular_momentum_body()).map(dd);
        let weight = components(self.gravity.acceleration_world()).map(|g| mass * dd(g));
        let after_l = before_l.map(|value| attenuated(value, rotation.exposure));
        let rotation_vector_body =
            vector(before_l.map(|value| value / inertia * rotation.integral_s))?;
        let orientation = state.pose().orientation().right_exp(rotation_vector_body)?;
        let position = components(state.pose().position_world()).map(dd);
        let state_after = RigidBodyState::new(
            Pose::new(
                vector(core::array::from_fn(|i| position[i] + displacement[i]))?,
                orientation,
            )?,
            vector(after_p)?,
            vector(after_l)?,
        )?;
        let rotation_heat = if self.rotary_resistance_kg_m2_s == 0.0 {
            Dd::ZERO
        } else {
            dot(before_l, before_l) / inertia
                * dd(0.5)
                * dd(decay::decayed_fraction(2.0 * rotation.exposure))
        };
        let heat = translation_heat + rotation_heat;
        if scalar(translation_heat)? < 0.0 || scalar(rotation_heat)? < 0.0 {
            return Err(LinearDragError::Invalid("negative integrated drag heat"));
        }
        let gravity_work = dot(weight, displacement);
        let represented_displacement = core::array::from_fn(|i| {
            dd(components(state_after.pose().position_world())[i]) - position[i]
        });
        let mut kinetic = RigidKineticEnergyChange::new();
        kinetic.add_body(self.properties, state, state_after)?;
        let kinetic_change = kinetic.finish_retained()?;
        let mechanical_change = kinetic_change - dot(weight, represented_displacement);
        let linear_residual = core::array::from_fn(|i| {
            dd(components(state_after.linear_momentum_world())[i])
                - before_p[i]
                - total_impulse[i]
                - drag_impulse[i]
        });
        let drag_angular_impulse_world = state
            .pose()
            .orientation()
            .rotate_body_to_world(vector(before_l.map(|value| -value * dd(rotation.lost)))?);
        let initial_world_l = state
            .pose()
            .orientation()
            .rotate_body_to_world(state.angular_momentum_body());
        let final_world_l = orientation.rotate_body_to_world(state_after.angular_momentum_body());
        let angular_residual = core::array::from_fn(|i| {
            dd(components(final_world_l)[i])
                - dd(components(initial_world_l)[i])
                - dd(components(drag_angular_impulse_world)[i])
        });
        Ok(LinearDragStep {
            duration_s,
            state_before: state,
            state_after,
            displacement_world: vector(displacement)?,
            rotation_vector_body,
            drag_impulse_world: vector(drag_impulse)?,
            drag_angular_impulse_world,
            linear_momentum_residual_world: vector(linear_residual)?,
            angular_momentum_residual_world: vector(angular_residual)?,
            translation_heat_j: scalar(translation_heat)?,
            rotation_heat_j: scalar(rotation_heat)?,
            dissipated_energy_j: scalar(heat)?,
            background_work_j: scalar(background_work)?,
            applied_work_j: scalar(applied_work)?,
            gravity_work_j: scalar(gravity_work)?,
            kinetic_energy_change_j: scalar(kinetic_change)?,
            mechanical_energy_change_j: scalar(mechanical_change)?,
            work_energy_residual_j: scalar(
                mechanical_change + heat - background_work - applied_work,
            )?,
        })
    }
}

struct TranslationPath {
    after_p: [Dd; 3],
    displacement: [Dd; 3],
    drag_impulse: [Dd; 3],
    translation_heat: Dd,
    applied_work: Dd,
    total_impulse: [Dd; 3],
}

// phi1 = integral_0^1 exp(-x*s) ds; phi2 = (1-phi1)/x.
// variance is that of (1-exp(-x*s))/(1-exp(-x)) on [0,1].
// Its small-x series tends to 1/12, the variance of s, without subtracting
// nearly equal moments. Transcendentals remain in the shared decay owner.
struct Moments {
    exposure: f64,
    lost: f64,
    phi1: Dd,
    phi2: f64,
    integral_s: Dd,
    drive_integral_s: Dd,
    variance: f64,
}

impl Moments {
    fn new(resistance: f64, inertia: f64, duration: f64) -> Result<Self, LinearDragError> {
        let exposure = relaxation_exposure(resistance, inertia, duration)?;
        let lost = decay::decayed_fraction(exposure);
        let (phi1, phi2, variance) = if exposure < 0.25 {
            let mut phi2 = 0.5;
            let mut term = phi2;
            for order in 1..=16 {
                term *= -exposure / f64::from(order + 2);
                phi2 += term;
            }
            let x2 = exposure * exposure;
            let variance = 1.0 / 12.0
                + x2 * (-1.0 / 720.0
                    + x2 * (1.0 / 30_240.0
                        + x2 * (-1.0 / 1_209_600.0
                            + x2 * (1.0 / 47_900_160.0 - x2 * 691.0 / 1_307_674_368_000.0))));
            (dd(1.0) - dd(exposure) * dd(phi2), phi2, variance)
        } else {
            let phi1 = lost / exposure;
            (
                dd(phi1),
                (1.0 - phi1) / exposure,
                (0.5 * (1.0 + decay::remaining_fraction(exposure)) / lost - 1.0 / exposure)
                    / exposure,
            )
        };
        Ok(Self {
            exposure,
            lost,
            phi1,
            phi2,
            integral_s: dd(duration) * phi1,
            drive_integral_s: dd(duration) * dd(phi2),
            variance,
        })
    }
}

fn relaxation_exposure(
    resistance: f64,
    inertia: f64,
    duration: f64,
) -> Result<f64, LinearDragError> {
    let rate = resistance / inertia;
    let exposure = rate * duration;
    if !rate.is_finite()
        || !exposure.is_finite()
        || resistance > 0.0 && (rate == 0.0 || exposure == 0.0)
    {
        return Err(LinearDragError::Invalid(
            "drag relaxation scale is unrepresentable",
        ));
    }
    Ok(exposure)
}

fn dd(value: f64) -> Dd {
    Dd::from_f64(value)
}
fn components(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}
fn attenuated(value: Dd, exposure: f64) -> Dd {
    dd(decay::attenuate(value.hi, exposure)) + dd(decay::attenuate(value.lo, exposure))
}
fn dot(a: [Dd; 3], b: [Dd; 3]) -> Dd {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn scalar(value: Dd) -> Result<f64, LinearDragError> {
    let rounded = value.to_f64();
    if value.hi.is_finite() && value.lo.is_finite() && rounded.is_finite() {
        Ok(rounded)
    } else {
        Err(LinearDragError::Invalid("unrepresentable flow arithmetic"))
    }
}
fn vector(value: [Dd; 3]) -> Result<Vec3, LinearDragError> {
    Ok(Vec3::new(
        scalar(value[0])?,
        scalar(value[1])?,
        scalar(value[2])?,
    ))
}
