//! Variable-mass rocket attitude on the original native rigid-body algebra.
//!
//! The declared mass distribution has a fixed common centre of mass and fixed
//! principal axes, with I = I_dry + m_propellant * K_propellant. The angular
//! momentum equation is
//! dH_body/dt + omega cross H_body = torque_body - mdot*(K_exit*omega + j_swirl).
//! The exit co-rotation moment is an explicit input, independent of tank inertia:
//! retaining angular momentum, retaining angular velocity, and specifying a
//! nozzle jet-damping moment are different physical models.
//!
//! Symmetric torque/outflow kicks surround the original momentum-preserving
//! free-rotation map at midpoint inertia. Each outflow kick integrates the
//! linear-in-propellant inertia exactly. Translation reuses the original
//! analytic rocket step along the midpoint body attitude. The smooth forced
//! composition is second order; timestep refinement is still required.
//! Work is bounded by an explicit substep limit, with burnout as an exact event.
//! No heap allocation, new orientation representation, or atmospheric law is
//! introduced. Internal slosh/circulation, shifting COM/principal axes and
//! non-diagonal outflow moments require a more general open-system model.
//! See Eke, NASA/CR-1998-208246, https://ntrs.nasa.gov/citations/19980210404.

use super::{RocketError, RocketForcing, RocketState};
use crate::{DynamicsError, MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

/// Instantaneous two-body separation with conserved mass and momenta.
#[path = "rigid/staging.rs"]
pub mod staging;

/// Declared concentric dry and propellant mass distributions in common axes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidRocketInertia {
    /// Positive dry principal moments about the common COM [kg m²].
    pub dry_principal_kg_m2: Vec3,
    /// Propellant principal moment per unit remaining mass [m²].
    /// Nonnegative, including a concentrated mass at the COM, is admitted.
    pub propellant_specific_principal_m2: Vec3,
}
impl RigidRocketInertia {
    /// Derive the original validated native mass properties from actual masses.
    pub fn mass_properties(
        self,
        translation: RocketState,
    ) -> Result<MassProperties, RigidRocketError> {
        if !translation.dry_mass_kg.is_finite()
            || translation.dry_mass_kg <= 0.0
            || !translation.propellant_mass_kg.is_finite()
            || translation.propellant_mass_kg < 0.0
        {
            return Err(RigidRocketError::Invalid(
                "positive dry and nonnegative propellant masses required",
            ));
        }
        MassProperties::new(
            translation.dry_mass_kg,
            Vec3::ZERO,
            self.dry_principal_kg_m2,
        )?;
        physical_diagonal(
            self.propellant_specific_principal_m2,
            "propellant specific inertia",
        )?;
        Ok(MassProperties::new(
            translation.dry_mass_kg + translation.propellant_mass_kg,
            Vec3::ZERO,
            self.dry_principal_kg_m2.add(
                self.propellant_specific_principal_m2
                    .scale(translation.propellant_mass_kg),
            ),
        )?)
    }
}

/// Native translation plus native body-to-world attitude and body spin momentum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidRocketState {
    translation: RocketState,
    orientation: UnitQuaternion,
    angular_momentum_body_kg_m2_s: Vec3,
    inertia: RigidRocketInertia,
}

/// Explicit forcing, fixed in its stated frame during a requested interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidRocketForcing {
    /// Propellant consumption rate during powered flight [kg/s].
    pub mass_flow_kg_s: f64,
    /// Total thrust divided by mass flow, pointing thrustward in the body [m/s].
    /// Includes pressure thrust once. Zero is allowed for spin-only venting.
    pub effective_exhaust_velocity_body_m_s: Vec3,
    /// Thrust application point relative to the common COM, in body axes [m].
    pub thrust_application_arm_body_m: Vec3,
    /// Diagonal outflow mass-weighted co-rotation moment per mass [m²].
    /// The exported spin is mdot * K_exit * omega; this is not tank inertia.
    /// For a uniform circular exit normal to x, centred at x=-l, it is
    /// (R²/2, l²+R²/4, l²+R²/4). Supply zero only for zero such spin flux.
    /// This diagonal convention does not describe an arbitrary off-axis outlet.
    pub co_rotating_exhaust_specific_inertia_m2: Vec3,
    /// Additional exported body angular momentum per mass [m²/s], e.g. swirl.
    /// Excludes co-rotation and relative jet momentum already in thrust.
    pub exhaust_swirl_specific_angular_momentum_body_m2_s: Vec3,
    /// Uniform world-frame gravitational acceleration [m/s²].
    pub gravity_world_m_s2: Vec3,
    /// Other world-frame force at the COM [N], excluding gravity and propulsion.
    pub external_force_world_n: Vec3,
    /// Other body-frame torque about the COM [N m], including any aerodynamic arm.
    pub external_torque_body_n_m: Vec3,
}

/// Explicit temporal resolution and work budget for smooth attitude/forcing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidRocketSolve {
    /// Maximum physical substep [s]. Burnout can introduce a shorter step.
    pub maximum_substep_s: f64,
    /// Maximum accepted substeps, counting a split at burnout.
    pub maximum_substeps: usize,
}

/// One complete trajectory candidate, with actual native-state accounting.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidRocketStep {
    /// Complete new state; never returned on a refusal or cancellation.
    pub state_after: RigidRocketState,
    /// Sum of the actual native translation steps' consumed propellant [kg].
    pub burned_propellant_mass_kg: f64,
    /// Initial minus final inventory minus accepted consumed mass [kg].
    pub propellant_balance_residual_kg: f64,
    /// Powered duration, excluding any dry coast [s].
    pub powered_duration_s: f64,
    /// Remaining requested duration with engine/outflow exactly off [s].
    pub coast_duration_s: f64,
    /// Integral of the body thrust magnitude [N s].
    /// A rotating vehicle world impulse norm can be smaller than this scalar.
    pub axial_thrust_impulse_n_s: f64,
    /// Integrated directed thrust in the world frame [N s].
    pub thrust_impulse_world_n_s: Vec3,
    /// Applied lever-arm and other torque impulses in the discrete world map [N m s].
    pub applied_angular_impulse_world_kg_m2_s: Vec3,
    /// Exported co-rotation and swirl angular momentum in that same map [kg m²/s].
    pub exported_angular_momentum_world_kg_m2_s: Vec3,
    /// Actual world spin change minus applied torque plus exported spin [kg m²/s].
    /// This audits discrete kicks/free maps, not the unknown exact trajectory.
    pub angular_momentum_residual_world_kg_m2_s: Vec3,
    /// Number of complete physical substeps, including a burnout split.
    pub substeps: usize,
}

/// Input, native-owner, arithmetic, work-budget, or cancellation refusal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RigidRocketError {
    /// A declared model/input or derived quantity is outside the admitted domain.
    Invalid(&'static str),
    /// The original translational owner refused.
    Translation(RocketError),
    /// The original mass/orientation/rotation owner refused.
    Dynamics(DynamicsError),
    /// The requested physical interval exceeded the declared work budget.
    SubstepBudget {
        /// Complete tentative substeps before refusal.
        completed: usize,
        /// Caller-declared cap.
        maximum: usize,
    },
    /// Cancellation was observed before publishing the complete state.
    Cancelled,
}
impl core::fmt::Display for RigidRocketError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid variable-mass rocket: {self:?}")
    }
}
impl std::error::Error for RigidRocketError {}
impl From<RocketError> for RigidRocketError {
    fn from(value: RocketError) -> Self {
        Self::Translation(value)
    }
}
impl From<DynamicsError> for RigidRocketError {
    fn from(value: DynamicsError) -> Self {
        Self::Dynamics(value)
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), RigidRocketError> {
    if cancelled() {
        Err(RigidRocketError::Cancelled)
    } else {
        Ok(())
    }
}
fn physical_diagonal(value: Vec3, field: &'static str) -> Result<(), RigidRocketError> {
    if !value.is_finite()
        || value.x < 0.0
        || value.y < 0.0
        || value.z < 0.0
        || value.x > value.y + value.z
        || value.y > value.x + value.z
        || value.z > value.x + value.y
    {
        return Err(RigidRocketError::Invalid(field));
    }
    Ok(())
}
fn zero_forcing() -> RocketForcing {
    RocketForcing {
        mass_flow_kg_s: 0.0,
        effective_exhaust_velocity_m_s: Vec3::ZERO,
        gravity_m_s2: Vec3::ZERO,
        external_force_n: Vec3::ZERO,
    }
}
#[derive(Default)]
struct VectorSum {
    value: Vec3,
    compensation: Vec3,
}
impl VectorSum {
    fn add(&mut self, increment: Vec3) {
        let adjusted = increment.sub(self.compensation);
        let next = self.value.add(adjusted);
        self.compensation = next.sub(self.value).sub(adjusted);
        self.value = next;
    }
}
#[derive(Default)]
struct ScalarSum {
    value: f64,
    compensation: f64,
}
impl ScalarSum {
    fn add(&mut self, increment: f64) {
        let adjusted = increment - self.compensation;
        let next = self.value + adjusted;
        self.compensation = (next - self.value) - adjusted;
        self.value = next;
    }
}

impl RigidRocketState {
    /// Assemble validated native owners without an independent propellant inventory.
    pub fn new(
        translation: RocketState,
        orientation: UnitQuaternion,
        angular_momentum_body_kg_m2_s: Vec3,
        inertia: RigidRocketInertia,
    ) -> Result<Self, RigidRocketError> {
        translation.advance(zero_forcing(), 0.0, || false)?;
        let result = Self {
            translation,
            orientation,
            angular_momentum_body_kg_m2_s,
            inertia,
        };
        result
            .mass_properties()?
            .angular_velocity_body_checked(angular_momentum_body_kg_m2_s)?;
        result.rigid_body()?;
        Ok(result)
    }
    /// Original variable-mass translational state, in the declared common world frame.
    #[must_use]
    pub const fn translation(&self) -> &RocketState {
        &self.translation
    }
    /// Original native body-to-world orientation.
    #[must_use]
    pub const fn orientation(&self) -> UnitQuaternion {
        self.orientation
    }
    /// Native principal-body-frame spin angular momentum [kg m²/s].
    #[must_use]
    pub const fn angular_momentum_body_kg_m2_s(&self) -> Vec3 {
        self.angular_momentum_body_kg_m2_s
    }
    /// Declared fixed-COM, common-principal-axis inertia model.
    #[must_use]
    pub const fn inertia(&self) -> RigidRocketInertia {
        self.inertia
    }
    /// Actual mass and inertia at the remaining propellant inventory.
    pub fn mass_properties(&self) -> Result<MassProperties, RigidRocketError> {
        self.inertia.mass_properties(self.translation)
    }
    /// Expose the original momentum-form rigid-body owner for existing consumers.
    pub fn rigid_body(&self) -> Result<RigidBodyState, RigidRocketError> {
        let mass = self.mass_properties()?.mass();
        Ok(RigidBodyState::new(
            Pose::new(self.translation.position_m, self.orientation)?,
            self.translation.velocity_m_s.scale(mass),
            self.angular_momentum_body_kg_m2_s,
        )?)
    }

    /// Advance body-directed propulsion, inertia and spin in one bounded transaction.
    ///
    /// Inputs remain fixed in their stated frames. Environment/drag/changing gimbal
    /// consumers re-evaluate them between accepted calls. Fixed-COM mass depletion
    /// uses the native mass owner, translation uses its closed-form mass integrals,
    /// and rotation/force-direction coupling is second order for smooth inputs.
    /// Cancellation or exhausted work returns no partially accepted state.
    pub fn advance(
        &self,
        forcing: RigidRocketForcing,
        duration_s: f64,
        controls: RigidRocketSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RigidRocketStep, RigidRocketError> {
        poll(&mut cancelled)?;
        if !duration_s.is_finite()
            || duration_s < 0.0
            || !controls.maximum_substep_s.is_finite()
            || controls.maximum_substep_s <= 0.0
            || controls.maximum_substeps == 0
        {
            return Err(RigidRocketError::Invalid(
                "finite nonnegative duration, positive substep and nonzero work budget required",
            ));
        }
        self.translation.advance(
            RocketForcing {
                mass_flow_kg_s: forcing.mass_flow_kg_s,
                effective_exhaust_velocity_m_s: forcing.effective_exhaust_velocity_body_m_s,
                gravity_m_s2: forcing.gravity_world_m_s2,
                external_force_n: forcing.external_force_world_n,
            },
            0.0,
            &mut cancelled,
        )?;
        physical_diagonal(
            forcing.co_rotating_exhaust_specific_inertia_m2,
            "co-rotating exhaust specific inertia",
        )?;
        for vector in [
            forcing.thrust_application_arm_body_m,
            forcing.exhaust_swirl_specific_angular_momentum_body_m2_s,
            forcing.external_torque_body_n_m,
        ] {
            if !vector.is_finite() {
                return Err(RigidRocketError::Invalid(
                    "finite arm, swirl and external torque required",
                ));
            }
        }
        let exhaust_speed = forcing
            .effective_exhaust_velocity_body_m_s
            .stable_norm("rocket effective exhaust speed")?;
        let initial_spin = self
            .orientation
            .rotate_body_to_world(self.angular_momentum_body_kg_m2_s);
        let mut current = *self;
        let (mut remaining, mut substeps) = (duration_s, 0);
        let (mut elapsed, mut powered, mut coast, mut consumed) = (
            ScalarSum::default(),
            ScalarSum::default(),
            ScalarSum::default(),
            ScalarSum::default(),
        );
        let (mut thrust, mut applied, mut exported) = (
            VectorSum::default(),
            VectorSum::default(),
            VectorSum::default(),
        );
        while remaining > 0.0 {
            poll(&mut cancelled)?;
            if substeps == controls.maximum_substeps {
                return Err(RigidRocketError::SubstepBudget {
                    completed: substeps,
                    maximum: controls.maximum_substeps,
                });
            }
            let mut h = remaining.min(controls.maximum_substep_s);
            let burning =
                forcing.mass_flow_kg_s > 0.0 && current.translation.propellant_mass_kg > 0.0;
            if burning {
                let burnout = current.translation.propellant_mass_kg / forcing.mass_flow_kg_s;
                if burnout < h {
                    h = burnout;
                }
            }
            if !h.is_finite() || h <= 0.0 || (remaining - h == remaining) {
                return Err(RigidRocketError::Invalid(
                    "physical substep is not representable",
                ));
            }
            let mut probe = zero_forcing();
            probe.mass_flow_kg_s = forcing.mass_flow_kg_s;
            let mass_step = current.translation.advance(probe, h, &mut cancelled)?;
            let burned = mass_step.burned_propellant_mass_kg;
            if mass_step.coast_duration_s != 0.0 && burned != 0.0 {
                return Err(RigidRocketError::Invalid(
                    "burnout split is not representable",
                ));
            }
            let rate = if burned > 0.0 { burned / h } else { 0.0 };
            let midpoint_propellant = current.translation.propellant_mass_kg - 0.5 * burned;
            let midpoint_translation = RocketState {
                propellant_mass_kg: midpoint_propellant,
                ..current.translation
            };
            let before_properties = current.mass_properties()?;
            let midpoint_properties = current.inertia.mass_properties(midpoint_translation)?;
            let after_properties = current.inertia.mass_properties(mass_step.state_after)?;
            let thrust_body = forcing.effective_exhaust_velocity_body_m_s.scale(rate);
            let torque = forcing
                .external_torque_body_n_m
                .add(forcing.thrust_application_arm_body_m.cross(thrust_body));
            let swirl = forcing
                .exhaust_swirl_specific_angular_momentum_body_m2_s
                .scale(0.5 * burned);
            let half_torque = torque.scale(0.5 * h);
            if !half_torque.is_finite() || !swirl.is_finite() {
                return Err(RigidRocketError::Invalid(
                    "torque or swirl impulse overflow",
                ));
            }
            let momentum = current
                .angular_momentum_body_kg_m2_s
                .add(half_torque)
                .sub(swirl);
            applied.add(current.orientation.rotate_body_to_world(half_torque));
            exported.add(current.orientation.rotate_body_to_world(swirl));
            let (momentum, loss) = exhaust_spin_step(
                momentum,
                before_properties.principal_inertia_body(),
                midpoint_properties.principal_inertia_body(),
                current.inertia.propellant_specific_principal_m2,
                forcing.co_rotating_exhaust_specific_inertia_m2,
                0.5 * burned,
            )?;
            exported.add(current.orientation.rotate_body_to_world(loss));

            let (midpoint_orientation, midpoint_momentum) = crate::free_rotation_step(
                current.orientation,
                momentum,
                midpoint_properties,
                0.5 * h,
            )?;
            let (orientation, momentum) = crate::free_rotation_step(
                midpoint_orientation,
                midpoint_momentum,
                midpoint_properties,
                0.5 * h,
            )?;
            let (momentum, loss) = exhaust_spin_step(
                momentum,
                midpoint_properties.principal_inertia_body(),
                after_properties.principal_inertia_body(),
                current.inertia.propellant_specific_principal_m2,
                forcing.co_rotating_exhaust_specific_inertia_m2,
                0.5 * burned,
            )?;
            let momentum = momentum.add(half_torque).sub(swirl);
            exported.add(orientation.rotate_body_to_world(loss.add(swirl)));
            applied.add(orientation.rotate_body_to_world(half_torque));

            let motion = current.translation.advance(
                RocketForcing {
                    mass_flow_kg_s: forcing.mass_flow_kg_s,
                    effective_exhaust_velocity_m_s: midpoint_orientation
                        .rotate_body_to_world(forcing.effective_exhaust_velocity_body_m_s),
                    gravity_m_s2: forcing.gravity_world_m_s2,
                    external_force_n: forcing.external_force_world_n,
                },
                h,
                &mut cancelled,
            )?;
            if motion.burned_propellant_mass_kg != burned {
                return Err(RigidRocketError::Invalid(
                    "native propellant inventory changed between mass and motion",
                ));
            }
            thrust.add(motion.thrust_impulse_n_s);
            powered.add(motion.powered_duration_s);
            coast.add(motion.coast_duration_s);
            consumed.add(motion.burned_propellant_mass_kg);
            current = Self::new(motion.state_after, orientation, momentum, current.inertia)?;
            if h == remaining {
                remaining = 0.0;
            } else {
                elapsed.add(h);
                remaining = duration_s - elapsed.value;
            }
            substeps += 1;
        }
        let burned = consumed.value;
        let propellant_residual =
            self.translation.propellant_mass_kg - current.translation.propellant_mass_kg - burned;
        let axial_impulse = exhaust_speed * burned;
        let residual = current
            .orientation
            .rotate_body_to_world(current.angular_momentum_body_kg_m2_s)
            .sub(initial_spin)
            .sub(applied.value)
            .add(exported.value);
        if !axial_impulse.is_finite()
            || !thrust.value.is_finite()
            || !applied.value.is_finite()
            || !exported.value.is_finite()
            || !residual.is_finite()
            || !propellant_residual.is_finite()
        {
            return Err(RigidRocketError::Invalid(
                "integrated rocket impulse is not representable",
            ));
        }
        poll(&mut cancelled)?;
        Ok(RigidRocketStep {
            state_after: current,
            burned_propellant_mass_kg: burned,
            propellant_balance_residual_kg: propellant_residual,
            powered_duration_s: powered.value,
            coast_duration_s: coast.value,
            axial_thrust_impulse_n_s: axial_impulse,
            thrust_impulse_world_n_s: thrust.value,
            applied_angular_impulse_world_kg_m2_s: applied.value,
            exported_angular_momentum_world_kg_m2_s: exported.value,
            angular_momentum_residual_world_kg_m2_s: residual,
            substeps,
        })
    }
}

// Solve dH_i/dt = -mdot * K_exit_i / I_i(t) * H_i with linear I_i(m).
// Compute retained spin with exp and exported spin with expm1 separately:
// subtraction must not erase a tiny retained spin or a tiny exported spin.
fn exhaust_spin_step(
    momentum: Vec3,
    inertia_before: Vec3,
    inertia_after: Vec3,
    propellant_moment: Vec3,
    exit_moment: Vec3,
    expelled: f64,
) -> Result<(Vec3, Vec3), RigidRocketError> {
    let inputs = [
        (
            momentum.x,
            inertia_before.x,
            inertia_after.x,
            propellant_moment.x,
            exit_moment.x,
        ),
        (
            momentum.y,
            inertia_before.y,
            inertia_after.y,
            propellant_moment.y,
            exit_moment.y,
        ),
        (
            momentum.z,
            inertia_before.z,
            inertia_after.z,
            propellant_moment.z,
            exit_moment.z,
        ),
    ];
    let mut retained = [momentum.x, momentum.y, momentum.z];
    let mut lost = [0.0; 3];
    for (i, (spin, before, after, bulk, exit)) in inputs.into_iter().enumerate() {
        if expelled == 0.0 || exit == 0.0 {
            continue;
        }
        let fraction = (expelled / before) * bulk;
        let multiplier = if fraction < 0.125 {
            super::mass_integrals(fraction, before, after).1
        } else {
            (before.ln() - after.ln()) / fraction
        };
        let exponent = ((expelled / before) * exit) * multiplier;
        if !fraction.is_finite()
            || fraction < 0.0
            || !multiplier.is_finite()
            || multiplier <= 0.0
            || exponent.is_nan()
            || exponent < 0.0
        {
            return Err(RigidRocketError::Invalid(
                "co-rotation exhaust integral is not representable",
            ));
        }
        retained[i] = spin * (-exponent).exp();
        lost[i] = spin * -(-exponent).exp_m1();
    }
    let after = Vec3::new(retained[0], retained[1], retained[2]);
    let loss = Vec3::new(lost[0], lost[1], lost[2]);
    if !after.is_finite() || !loss.is_finite() {
        return Err(RigidRocketError::Invalid(
            "co-rotation exhaust angular momentum overflow",
        ));
    }
    Ok((after, loss))
}

#[cfg(test)]
#[path = "rigid/tests.rs"]
mod tests;
