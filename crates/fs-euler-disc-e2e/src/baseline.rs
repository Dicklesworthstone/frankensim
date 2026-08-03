//! Executable ideal-conservative Euler-disc baseline.
//!
//! This is deliberately a narrow, contact-constrained reference rung: a
//! parameterized squat rigid disc follows a steady, ideal no-slip rolling
//! reduction on a horizontal support. The support geometry is the lowest rim
//! of the finite-radius disc; gravity supplies the normal load, and the
//! requested path terminates when the static-friction capacity is insufficient.
//! The attitude update composes the generic `fs-time` Lie-group stepper. This
//! does not establish dissipation, impacts, finite stop time, or any
//! physical-video correspondence.

use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
};
use fs_time::{quat_exp_step, quat_rotate};
use fs_tribo::{
    ContactFrame, FrictionLaw, InputAuthority, InterfaceMedium, InterfaceSystemRef, TangentialSlip,
};

/// Gravitational acceleration used by the baseline, in m/s².
pub const STANDARD_GRAVITY_MPS2: f64 = 9.806_65;
/// Execution bound that keeps one invocation's retained trajectory finite.
pub const MAX_BASELINE_STEPS: u32 = 100_000;

/// Scope assigned to the executed reduced path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaselineDynamicsClass {
    /// Thin homogeneous disc, small inclination, steady no-slip rolling oracle.
    ThinHomogeneousSmallAngleSteadyRolling,
    /// Geometrically supported no-slip path whose rates are prescribed, not dynamics.
    PrescribedKinematicPath,
}

impl BaselineDynamicsClass {
    const fn model_id(self) -> &'static str {
        match self {
            Self::ThinHomogeneousSmallAngleSteadyRolling => {
                "ideal_conservative_thin_homogeneous_steady_rolling_disc"
            }
            Self::PrescribedKinematicPath => "prescribed_no_slip_rolling_disc_kinematics",
        }
    }
}

/// Parameterized physical input for a squat, axisymmetric rigid disc.
#[derive(Debug, Clone, PartialEq)]
pub struct SquatDiscInput {
    /// Declared scope: dynamic thin-disc oracle or explicitly kinematic path.
    pub dynamics_class: BaselineDynamicsClass,
    /// Outer radius, in metres.
    pub radius_m: f64,
    /// Axial thickness, in metres.
    pub thickness_m: f64,
    /// Total mass, in kilograms.
    pub mass_kg: f64,
    /// Principal moments `(transverse, transverse, axial)`, in kg·m².
    pub inertia_kg_m2: [f64; 3],
    /// Constant inclination of the symmetry axis from the support normal, in rad.
    pub inclination_from_vertical_rad: f64,
    /// Constant precession about the support normal, in rad/s.
    pub precession_rad_s: f64,
    /// Constant spin about the disc symmetry axis, in rad/s.
    pub spin_rad_s: f64,
    /// Available static-friction coefficient at the ideal dry support.
    pub static_friction_coefficient: f64,
    /// Ordered disc/support interface identity consumed by `fs-tribo`.
    pub interface_system_id: String,
    /// Declared history identity consumed by `fs-tribo`.
    pub interface_history_id: String,
    /// Caller-declared source identity retained by `fs-tribo`; not a material receipt.
    pub interface_source_id: String,
    /// Authority carried with the dry-contact inputs; it is retained, not upgraded.
    pub interface_authority: InputAuthority,
    /// Initial physical state at t = 0.
    pub initial_state: BaselineState,
    /// Fixed integration increment, in seconds.
    pub step_seconds: f64,
    /// Number of production integration steps to execute.
    pub steps: u32,
}

impl SquatDiscInput {
    /// A constructible, short ideal-rolling fixture; callers may parameterize
    /// all physical fields before running the production baseline.
    #[must_use]
    pub fn nominal() -> Self {
        let radius_m: f64 = 0.038;
        let thickness_m: f64 = 0.000_01;
        let mass_kg: f64 = 0.120;
        let transverse = mass_kg * (3.0 * radius_m * radius_m + thickness_m * thickness_m) / 12.0;
        let axial = 0.5 * mass_kg * radius_m * radius_m;
        let inclination_from_vertical_rad: f64 = 0.08;
        let precession_rad_s =
            (4.0 * STANDARD_GRAVITY_MPS2 / (radius_m * inclination_from_vertical_rad.sin())).sqrt();
        let orientation_body_to_world = [
            (0.5 * inclination_from_vertical_rad).cos(),
            0.0,
            (0.5 * inclination_from_vertical_rad).sin(),
            0.0,
        ];
        let mut input = Self {
            dynamics_class: BaselineDynamicsClass::ThinHomogeneousSmallAngleSteadyRolling,
            radius_m,
            thickness_m,
            mass_kg,
            inertia_kg_m2: [transverse, transverse, axial],
            inclination_from_vertical_rad,
            precession_rad_s,
            spin_rad_s: -precession_rad_s * inclination_from_vertical_rad.cos(),
            static_friction_coefficient: 1.0,
            interface_system_id: "baseline/squat-disc->planar-support".to_owned(),
            interface_history_id: "baseline/ideal-no-slip-v1".to_owned(),
            interface_source_id: "baseline/caller-declared-thin-disc-oracle".to_owned(),
            interface_authority: InputAuthority::CallerDeclared,
            initial_state: BaselineState {
                time_seconds: 0.0,
                position_m: [0.0; 3],
                linear_velocity_mps: [0.0; 3],
                orientation_body_to_world,
                angular_velocity_body_rad_s: [0.0; 3],
            },
            step_seconds: 0.000_05,
            steps: 100,
        };
        input.initial_state = supported_state(&input, orientation_body_to_world, 0.0, [0.0; 3]);
        input
    }

    /// Rebuild the t = 0 state on the declared rim-contact/no-slip manifold
    /// after changing the prescribed kinematic rates or geometry.
    pub fn reset_supported_initial_state(&mut self) {
        self.initial_state = supported_state(
            self,
            self.initial_state.orientation_body_to_world,
            0.0,
            self.initial_state.position_m,
        );
    }
}

/// Physical state retained at one deterministic sample point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BaselineState {
    /// Simulation time, in seconds.
    pub time_seconds: f64,
    /// Centre-of-mass position in a world frame, in metres.
    pub position_m: [f64; 3],
    /// Centre-of-mass velocity in the world frame, in m/s.
    pub linear_velocity_mps: [f64; 3],
    /// Unit quaternion `(w, x, y, z)` mapping body vectors into the world frame.
    pub orientation_body_to_world: [f64; 4],
    /// Angular velocity in the body frame, in rad/s.
    pub angular_velocity_body_rad_s: [f64; 3],
}

/// Energy accounting at one sample point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BaselineEnergyLedger {
    /// Translational kinetic energy, in J.
    pub translational_kinetic_j: f64,
    /// Rotational kinetic energy, in J.
    pub rotational_kinetic_j: f64,
    /// Gravitational potential energy relative to z = 0, in J.
    pub gravitational_potential_j: f64,
    /// Work done by the ideal support constraint, in J. It is identically zero.
    pub ideal_support_work_j: f64,
    /// Sum of all energy channels, in J.
    pub total_j: f64,
    /// Difference from the initial total energy, in J.
    pub residual_from_initial_j: f64,
}

/// Geometric support/contact diagnostics for one sample.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BaselineSupportDiagnostic {
    /// Rim contact point in the world frame, in metres.
    pub contact_point_m: [f64; 3],
    /// Signed plane gap at the contact point, in metres; zero means supported.
    pub plane_gap_m: f64,
    /// Kinematic no-slip residual at the contact point, in m/s.
    pub no_slip_residual_mps: f64,
    /// Normal reaction balancing gravity, in N.
    pub normal_force_n: f64,
    /// Required horizontal static-contact force, in N.
    pub required_tangential_force_n: f64,
    /// Available static-contact force, in N.
    pub available_static_friction_n: f64,
    /// Ideal support work, in J.
    pub support_work_j: f64,
}

/// Dynamic-oracle audit retained with each trajectory.
///
/// Its thin-disc scope uses `theta` as the angle between the disc symmetry axis
/// and vertical and `Omega` as world-vertical precession. In that convention,
/// the admitted steady no-slip small-angle oracle checks
/// `Omega² sin(theta) = 4 g / R` and `E ≈ 3/2 m g R theta`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BaselineEquilibriumReceipt {
    /// Whether this output is admitted as the narrow dynamic oracle or only kinematics.
    pub dynamics_class: BaselineDynamicsClass,
    /// `Omega² sin(theta) - 4g/R`, in s⁻².
    pub precession_balance_residual_s_inv2: f64,
    /// `spin + Omega cos(theta)`, in rad/s; zero keeps the CM stationary ideally.
    pub spin_closure_residual_rad_s: f64,
    /// Thin-disc small-angle energy approximation, in J.
    pub small_angle_energy_oracle_j: f64,
    /// Initial MBD mechanical energy minus the small-angle oracle, in J.
    pub small_angle_energy_residual_j: f64,
}

impl BaselineEquilibriumReceipt {
    fn is_admitted_dynamic_oracle(self, input: &SquatDiscInput) -> bool {
        if self.dynamics_class != BaselineDynamicsClass::ThinHomogeneousSmallAngleSteadyRolling {
            return false;
        }
        let theta = input.inclination_from_vertical_rad;
        let expected_transverse = 0.25 * input.mass_kg * input.radius_m * input.radius_m;
        let expected_axial = 0.5 * input.mass_kg * input.radius_m * input.radius_m;
        let thin_and_small_angle = input.thickness_m / input.radius_m <= 0.01 && theta <= 0.20;
        let homogeneous_inertia = relative_error(input.inertia_kg_m2[0], expected_transverse)
            <= 0.02
            && relative_error(input.inertia_kg_m2[1], expected_transverse) <= 0.02
            && relative_error(input.inertia_kg_m2[2], expected_axial) <= 0.02;
        let precession_scale = 4.0 * STANDARD_GRAVITY_MPS2 / input.radius_m;
        let spin_scale = input.precession_rad_s.abs().max(1.0);
        let energy_scale = self.small_angle_energy_oracle_j.abs().max(1.0e-12);
        thin_and_small_angle
            && homogeneous_inertia
            && self.precession_balance_residual_s_inv2.abs() <= 1.0e-10 * precession_scale
            && self.spin_closure_residual_rad_s.abs() <= 1.0e-10 * spin_scale
            && self.small_angle_energy_residual_j.abs() <= 0.03 * energy_scale
    }
}

/// One state plus its contemporaneous energy accounting.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BaselineSample {
    /// Zero-based logical integration step.
    pub step: u32,
    /// Physical state.
    pub state: BaselineState,
    /// Energy ledger for `state`.
    pub energy: BaselineEnergyLedger,
    /// Geometric and force-capacity diagnostics for the support constraint.
    pub support: BaselineSupportDiagnostic,
}

/// Terminal state for the ideal rolling baseline.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BaselineTerminal {
    /// The requested bounded time horizon was reached while static support held.
    TimeHorizonReached { completed_steps: u32 },
    /// The ideal no-slip support can no longer supply the required tangent force.
    StaticFrictionCapacityExceeded {
        /// First step outside the declared support capacity.
        step: u32,
        /// Tangential force demanded by the reduced path, in N.
        required_tangential_force_n: f64,
        /// Declared capacity, in N.
        available_static_friction_n: f64,
    },
}

impl BaselineTerminal {
    const fn code(self) -> &'static str {
        match self {
            Self::TimeHorizonReached { .. } => "time_horizon_reached",
            Self::StaticFrictionCapacityExceeded { .. } => "static_friction_capacity_exceeded",
        }
    }

    const fn disposition(self) -> &'static str {
        match self {
            Self::TimeHorizonReached { .. } => "completed",
            Self::StaticFrictionCapacityExceeded { .. } => "terminated",
        }
    }
}

/// Deterministic input-admission refusal categories.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaselineRefusalReason {
    /// A physical or integration input contained a non-finite value.
    NonFiniteInput,
    /// A positive geometric extent, mass, inertia, or timestep was required.
    NonPositivePhysicalParameter,
    /// The supplied thickness is not bounded by the supplied radius.
    NotASquatDisc,
    /// The requested production trajectory exceeds this rung's retention budget.
    StepBudgetExceeded,
    /// The initial orientation does not have a positive finite norm.
    InvalidOrientation,
    /// The input state does not lie on the declared support/no-slip manifold.
    InitialStateViolatesSupportConstraint,
    /// The committed dry-contact production API refused the declared interface.
    TribologyInputRejected,
    /// The committed rigid-body production API refused the declared state.
    MultibodyInputRejected,
    /// A dynamic-oracle invocation failed its thin homogeneous steady-rolling checks.
    DynamicOracleViolation,
}

impl BaselineRefusalReason {
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::NonFiniteInput => "non_finite_input",
            Self::NonPositivePhysicalParameter => "non_positive_physical_parameter",
            Self::NotASquatDisc => "not_a_squat_disc",
            Self::StepBudgetExceeded => "step_budget_exceeded",
            Self::InvalidOrientation => "invalid_orientation",
            Self::InitialStateViolatesSupportConstraint => {
                "initial_state_violates_support_constraint"
            }
            Self::TribologyInputRejected => "tribology_input_rejected",
            Self::MultibodyInputRejected => "multibody_input_rejected",
            Self::DynamicOracleViolation => "dynamic_oracle_violation",
        }
    }
}

/// Structured refusal, retaining no partial trajectory as evidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BaselineRefusal {
    /// Stable machine-readable reason.
    pub reason: BaselineRefusalReason,
}

/// Completed, deterministic production trajectory.
#[derive(Debug, Clone, PartialEq)]
pub struct BaselineTrajectory {
    /// Explicit ideal no-slip rolling model identifier.
    pub model_id: &'static str,
    /// Input used to produce the trajectory.
    pub input: SquatDiscInput,
    /// Ordered samples, including initial state and terminal state.
    pub samples: Vec<BaselineSample>,
    /// Dynamic-oracle receipt or explicit kinematic classification.
    pub equilibrium: BaselineEquilibriumReceipt,
    /// Completion disposition.
    pub terminal: BaselineTerminal,
}

/// Structured output for both accepted and refused production invocations.
#[derive(Debug, Clone, PartialEq)]
pub enum BaselineRunOutput {
    /// A real integration reached its requested bounded horizon.
    Completed(BaselineTrajectory),
    /// Admission refused before integration, with no partial trajectory.
    Refused(BaselineRefusal),
}

impl BaselineRunOutput {
    /// Deterministic one-line JSON suitable for a runner or retained log.
    #[must_use]
    pub fn structured_output(&self) -> String {
        match self {
            Self::Completed(trajectory) => {
                let final_sample = trajectory
                    .samples
                    .last()
                    .expect("completed trajectory has initial sample");
                format!(
                    concat!(
                        "{{\"model\":\"{}\",\"disposition\":\"{}\",",
                        "\"terminal\":\"{}\",\"samples\":{},",
                        "\"time_seconds\":{:.17e},\"energy_residual_j\":{:.17e},",
                        "\"precession_balance_residual_s_inv2\":{:.17e}}}"
                    ),
                    trajectory.model_id,
                    trajectory.terminal.disposition(),
                    trajectory.terminal.code(),
                    trajectory.samples.len(),
                    final_sample.state.time_seconds,
                    final_sample.energy.residual_from_initial_j,
                    trajectory.equilibrium.precession_balance_residual_s_inv2,
                )
            }
            Self::Refused(refusal) => format!(
                "{{\"model\":\"ideal_conservative_no_slip_rolling_disc\",\"disposition\":\"refused\",\"reason\":\"{}\"}}",
                refusal.reason.code()
            ),
        }
    }
}

/// Run the ideal conservative no-slip rolling-disc production baseline.
///
/// This integrates the reduced, constant-inclination steady rolling path. The
/// finite-radius rim is constrained to z = 0, translational velocity obeys
/// `v_cm + omega × r_contact = 0`, gravity sets normal load, and the path ends
/// if the inferred horizontal support demand exceeds static-friction capacity.
#[must_use]
pub fn run_ideal_conservative_baseline(input: SquatDiscInput) -> BaselineRunOutput {
    if let Some(reason) = validate_input(&input) {
        return BaselineRunOutput::Refused(BaselineRefusal { reason });
    }

    let mut state = input.initial_state;
    let initial_total_j = energy(&input, state, 0.0).total_j;
    let equilibrium = equilibrium_receipt(&input, initial_total_j);
    let mut samples = Vec::with_capacity(usize::try_from(input.steps).unwrap_or(0) + 1);
    samples.push(BaselineSample {
        step: 0,
        state,
        energy: energy(&input, state, initial_total_j),
        support: support_diagnostic(&input, state, [0.0; 3]),
    });

    let mut terminal = BaselineTerminal::TimeHorizonReached {
        completed_steps: input.steps,
    };
    for step in 1..=input.steps {
        let omega_world = world_angular_velocity(&input, state.orientation_body_to_world);
        let omega_body = world_to_body(state.orientation_body_to_world, omega_world);
        let orientation_body_to_world = quat_exp_step(
            state.orientation_body_to_world,
            omega_body,
            input.step_seconds,
        );
        let next_kinematic = supported_state(
            &input,
            orientation_body_to_world,
            f64::from(step) * input.step_seconds,
            state.position_m,
        );
        let position_m = add3(
            state.position_m,
            scale3(
                add3(
                    state.linear_velocity_mps,
                    next_kinematic.linear_velocity_mps,
                ),
                0.5 * input.step_seconds,
            ),
        );
        state = BaselineState {
            position_m,
            ..next_kinematic
        };
        let acceleration_mps2 = scale3(
            sub3(
                state.linear_velocity_mps,
                samples
                    .last()
                    .expect("initial sample")
                    .state
                    .linear_velocity_mps,
            ),
            1.0 / input.step_seconds,
        );
        let support = support_diagnostic(&input, state, acceleration_mps2);
        let capacity_exceeded =
            support.required_tangential_force_n > support.available_static_friction_n + 1.0e-12;
        samples.push(BaselineSample {
            step,
            state,
            energy: energy(&input, state, initial_total_j),
            support,
        });
        if capacity_exceeded {
            terminal = BaselineTerminal::StaticFrictionCapacityExceeded {
                step,
                required_tangential_force_n: support.required_tangential_force_n,
                available_static_friction_n: support.available_static_friction_n,
            };
            break;
        }
    }

    BaselineRunOutput::Completed(BaselineTrajectory {
        model_id: input.dynamics_class.model_id(),
        input,
        samples,
        equilibrium,
        terminal,
    })
}

fn validate_input(input: &SquatDiscInput) -> Option<BaselineRefusalReason> {
    if !input.radius_m.is_finite()
        || !input.thickness_m.is_finite()
        || !input.mass_kg.is_finite()
        || !input.step_seconds.is_finite()
        || !input.inclination_from_vertical_rad.is_finite()
        || !input.precession_rad_s.is_finite()
        || !input.spin_rad_s.is_finite()
        || !input.static_friction_coefficient.is_finite()
        || !all_finite(&input.inertia_kg_m2)
        || !state_is_finite(input.initial_state)
    {
        return Some(BaselineRefusalReason::NonFiniteInput);
    }
    if input.radius_m <= 0.0
        || input.thickness_m <= 0.0
        || input.mass_kg <= 0.0
        || input.step_seconds <= 0.0
        || input.static_friction_coefficient < 0.0
        || input.inertia_kg_m2.iter().any(|value| *value <= 0.0)
    {
        return Some(BaselineRefusalReason::NonPositivePhysicalParameter);
    }
    if input.thickness_m > input.radius_m {
        return Some(BaselineRefusalReason::NotASquatDisc);
    }
    if input.inclination_from_vertical_rad <= 1.0e-6
        || input.inclination_from_vertical_rad >= core::f64::consts::FRAC_PI_2 - 1.0e-6
    {
        return Some(BaselineRefusalReason::NotASquatDisc);
    }
    if input.steps > MAX_BASELINE_STEPS {
        return Some(BaselineRefusalReason::StepBudgetExceeded);
    }
    let orientation_norm_squared = dot4(input.initial_state.orientation_body_to_world);
    if !orientation_norm_squared.is_finite() || (orientation_norm_squared - 1.0).abs() > 1.0e-12 {
        return Some(BaselineRefusalReason::InvalidOrientation);
    }
    if mbd_diagnostics(input, input.initial_state).is_none() {
        return Some(BaselineRefusalReason::MultibodyInputRejected);
    }
    if static_friction_capacity(input).is_none() {
        return Some(BaselineRefusalReason::TribologyInputRejected);
    }
    if input.dynamics_class == BaselineDynamicsClass::ThinHomogeneousSmallAngleSteadyRolling
        && !equilibrium_receipt(
            input,
            mbd_diagnostics(input, input.initial_state)?.mechanical_energy,
        )
        .is_admitted_dynamic_oracle(input)
    {
        return Some(BaselineRefusalReason::DynamicOracleViolation);
    }
    let normal = disc_normal(input.initial_state.orientation_body_to_world);
    if (normal.z - input.inclination_from_vertical_rad.cos()).abs() > 1.0e-10 {
        return Some(BaselineRefusalReason::InitialStateViolatesSupportConstraint);
    }
    let expected = supported_state(
        input,
        input.initial_state.orientation_body_to_world,
        0.0,
        input.initial_state.position_m,
    );
    if (input.initial_state.position_m[2] - expected.position_m[2]).abs() > 1.0e-10
        || norm3(sub3(
            input.initial_state.linear_velocity_mps,
            expected.linear_velocity_mps,
        )) > 1.0e-10
        || norm3(sub3(
            input.initial_state.angular_velocity_body_rad_s,
            expected.angular_velocity_body_rad_s,
        )) > 1.0e-10
    {
        return Some(BaselineRefusalReason::InitialStateViolatesSupportConstraint);
    }
    None
}

fn energy(
    input: &SquatDiscInput,
    state: BaselineState,
    initial_total_j: f64,
) -> BaselineEnergyLedger {
    let diagnostics =
        mbd_diagnostics(input, state).expect("validated baseline state reaches fs-mbd");
    let total_j = diagnostics.mechanical_energy;
    BaselineEnergyLedger {
        translational_kinetic_j: diagnostics.translational_kinetic_energy,
        rotational_kinetic_j: diagnostics.rotational_kinetic_energy,
        gravitational_potential_j: diagnostics.gravitational_potential_energy,
        ideal_support_work_j: 0.0,
        total_j,
        residual_from_initial_j: total_j - initial_total_j,
    }
}

fn equilibrium_receipt(
    input: &SquatDiscInput,
    initial_mechanical_energy_j: f64,
) -> BaselineEquilibriumReceipt {
    let theta = input.inclination_from_vertical_rad;
    let small_angle_energy_oracle_j =
        1.5 * input.mass_kg * STANDARD_GRAVITY_MPS2 * input.radius_m * theta;
    BaselineEquilibriumReceipt {
        dynamics_class: input.dynamics_class,
        precession_balance_residual_s_inv2: input.precession_rad_s.mul_add(
            input.precession_rad_s * theta.sin(),
            -4.0 * STANDARD_GRAVITY_MPS2 / input.radius_m,
        ),
        spin_closure_residual_rad_s: input.spin_rad_s + input.precession_rad_s * theta.cos(),
        small_angle_energy_oracle_j,
        small_angle_energy_residual_j: initial_mechanical_energy_j - small_angle_energy_oracle_j,
    }
}

fn mbd_diagnostics(
    input: &SquatDiscInput,
    state: BaselineState,
) -> Option<fs_mbd::DynamicsDiagnostics> {
    let orientation = UnitQuaternion::new(
        state.orientation_body_to_world[0],
        state.orientation_body_to_world[1],
        state.orientation_body_to_world[2],
        state.orientation_body_to_world[3],
    )
    .ok()?;
    let properties = MassProperties::new(
        input.mass_kg,
        Vec3::ZERO,
        Vec3::new(
            input.inertia_kg_m2[0],
            input.inertia_kg_m2[1],
            input.inertia_kg_m2[2],
        ),
    )
    .ok()?;
    let pose = Pose::new(
        Vec3::new(
            state.position_m[0],
            state.position_m[1],
            state.position_m[2],
        ),
        orientation,
    )
    .ok()?;
    let state = RigidBodyState::new(
        pose,
        Vec3::new(
            input.mass_kg * state.linear_velocity_mps[0],
            input.mass_kg * state.linear_velocity_mps[1],
            input.mass_kg * state.linear_velocity_mps[2],
        ),
        Vec3::new(
            input.inertia_kg_m2[0] * state.angular_velocity_body_rad_s[0],
            input.inertia_kg_m2[1] * state.angular_velocity_body_rad_s[1],
            input.inertia_kg_m2[2] * state.angular_velocity_body_rad_s[2],
        ),
    )
    .ok()?;
    let gravity = Gravity::new(Vec3::new(0.0, 0.0, -STANDARD_GRAVITY_MPS2)).ok()?;
    RigidBodyIntegrator::new(gravity)
        .diagnostics(state, properties)
        .ok()
}

fn interface_reference(input: &SquatDiscInput) -> Option<InterfaceSystemRef> {
    InterfaceSystemRef::new(
        input.interface_system_id.clone(),
        input.interface_history_id.clone(),
        input.interface_source_id.clone(),
        input.interface_authority,
        InterfaceMedium::Dry,
    )
    .ok()
}

fn static_friction_capacity(input: &SquatDiscInput) -> Option<f64> {
    let law = FrictionLaw::Coulomb {
        static_mu: input.static_friction_coefficient,
        kinetic_mu: input.static_friction_coefficient,
    };
    let frame = ContactFrame::new([0.0, 0.0, 1.0]).ok()?;
    let zero_slip = TangentialSlip::new(&frame, [0.0; 3]).ok()?;
    law.evaluate(
        &interface_reference(input)?,
        input.mass_kg * STANDARD_GRAVITY_MPS2,
        zero_slip,
    )
    .ok()
    .map(|response| response.static_limit)
}

fn supported_state(
    input: &SquatDiscInput,
    orientation_body_to_world: [f64; 4],
    time_seconds: f64,
    previous_position_m: [f64; 3],
) -> BaselineState {
    let offset = rim_contact_offset(
        input.radius_m,
        input.thickness_m,
        disc_normal(orientation_body_to_world),
    );
    let omega_world = world_angular_velocity(input, orientation_body_to_world);
    BaselineState {
        time_seconds,
        position_m: [previous_position_m[0], previous_position_m[1], -offset.z],
        linear_velocity_mps: vec3_array(omega_world.cross(offset).scale(-1.0)),
        orientation_body_to_world,
        angular_velocity_body_rad_s: world_to_body(orientation_body_to_world, omega_world),
    }
}

fn support_diagnostic(
    input: &SquatDiscInput,
    state: BaselineState,
    acceleration_mps2: [f64; 3],
) -> BaselineSupportDiagnostic {
    let normal = disc_normal(state.orientation_body_to_world);
    let contact_offset = rim_contact_offset(input.radius_m, input.thickness_m, normal);
    let contact_point_m = add3(state.position_m, vec3_array(contact_offset));
    let omega_world = world_angular_velocity(input, state.orientation_body_to_world);
    let no_slip_residual_mps = norm3(add3(
        state.linear_velocity_mps,
        vec3_array(omega_world.cross(contact_offset)),
    ));
    let normal_force_n = input.mass_kg * STANDARD_GRAVITY_MPS2;
    let required_tangential_force_n = input.mass_kg
        * (acceleration_mps2[0].mul_add(
            acceleration_mps2[0],
            acceleration_mps2[1] * acceleration_mps2[1],
        ))
        .sqrt();
    BaselineSupportDiagnostic {
        contact_point_m,
        plane_gap_m: contact_point_m[2],
        no_slip_residual_mps,
        normal_force_n,
        required_tangential_force_n,
        available_static_friction_n: static_friction_capacity(input)
            .expect("validated baseline reaches fs-tribo"),
        support_work_j: 0.0,
    }
}

fn world_angular_velocity(input: &SquatDiscInput, orientation_body_to_world: [f64; 4]) -> Vec3 {
    let normal = disc_normal(orientation_body_to_world);
    Vec3::new(0.0, 0.0, input.precession_rad_s).add(normal.scale(input.spin_rad_s))
}

fn disc_normal(orientation_body_to_world: [f64; 4]) -> Vec3 {
    let vector = quat_rotate(orientation_body_to_world, [0.0, 0.0, 1.0]);
    Vec3::new(vector[0], vector[1], vector[2])
}

fn world_to_body(orientation_body_to_world: [f64; 4], vector_world: Vec3) -> [f64; 3] {
    let conjugate = [
        orientation_body_to_world[0],
        -orientation_body_to_world[1],
        -orientation_body_to_world[2],
        -orientation_body_to_world[3],
    ];
    quat_rotate(conjugate, [vector_world.x, vector_world.y, vector_world.z])
}

fn rim_contact_offset(radius_m: f64, thickness_m: f64, normal: Vec3) -> Vec3 {
    let in_plane_vertical = Vec3::new(
        -normal.z * normal.x,
        -normal.z * normal.y,
        1.0 - normal.z * normal.z,
    );
    let inverse_norm = in_plane_vertical.norm_squared().sqrt().recip();
    in_plane_vertical
        .scale(-radius_m * inverse_norm)
        .sub(normal.scale(0.5 * thickness_m))
}

fn state_is_finite(state: BaselineState) -> bool {
    state.time_seconds.is_finite()
        && all_finite(&state.position_m)
        && all_finite(&state.linear_velocity_mps)
        && all_finite(&state.orientation_body_to_world)
        && all_finite(&state.angular_velocity_body_rad_s)
}

fn all_finite<const N: usize>(values: &[f64; N]) -> bool {
    values.iter().all(|value| value.is_finite())
}

fn dot3(values: [f64; 3]) -> f64 {
    values[0].mul_add(
        values[0],
        values[1].mul_add(values[1], values[2] * values[2]),
    )
}

fn dot4(values: [f64; 4]) -> f64 {
    values[0].mul_add(
        values[0],
        values[1].mul_add(
            values[1],
            values[2].mul_add(values[2], values[3] * values[3]),
        ),
    )
}

fn vec3_array(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}

fn add3(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [left[0] + right[0], left[1] + right[1], left[2] + right[2]]
}

fn sub3(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
}

fn scale3(values: [f64; 3], scalar: f64) -> [f64; 3] {
    [values[0] * scalar, values[1] * scalar, values[2] * scalar]
}

fn norm3(values: [f64; 3]) -> f64 {
    dot3(values).sqrt()
}

fn relative_error(actual: f64, expected: f64) -> f64 {
    (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE)
}
