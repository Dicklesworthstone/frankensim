//! A bounded, dynamic unilateral-contact rung for a finite-thickness disc.
//!
//! The support is the fixed plane `z = 0`.  Each fixed step evolves full
//! momentum-form rigid-body state under gravity, solves a point normal impulse,
//! and admits only a Coulomb-static (sticking) tangential impulse.  The contact
//! point is the lowest rim point of a homogeneous finite cylinder, derived from
//! the current orientation rather than from a prescribed inclination curve.
//!
//! This is intentionally not a general contact engine.  It does not model a
//! finite contact patch, compliance, restitutional impact law, sliding, rolling
//! resistance, aerodynamic drag, or an asymptotic Euler-disc decay law. A failure of the
//! sticking cone or a separating contact terminates the run; it is not silently
//! replaced with a kinematic trajectory.

#![deny(unsafe_code)]

use core::fmt;

use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::{
    ContactFrame, FrictionLaw, FrictionRegime, InputAuthority, InterfaceSystemRef, TangentialSlip,
};

const GROUND_NORMAL: Vec3 = Vec3::new(0.0, 0.0, 1.0);
const MAX_STEPS: u32 = 20_000;
const GEOMETRY_EPSILON: f64 = 128.0 * f64::EPSILON;

/// The explicit no-claim boundary of this mechanics rung.
pub const NO_CLAIM_BOUNDARY: &str = "Point contact with a rigid horizontal plane; sticking only. \
    No sliding, restitutional-impact law, finite-patch, aerodynamic, rolling-resistance, asymptotic Euler-disc, \
    experimental-validation, or convergence-order claim is made.";

/// Refusal from the bounded unilateral-contact model.
#[derive(Debug, Clone, PartialEq)]
pub enum ContactDynamicsError {
    /// A named scalar or vector is non-finite or outside its documented domain.
    InvalidInput { field: &'static str },
    /// A finite intermediate result cannot be represented as an `f64`.
    NonFiniteDerived { field: &'static str },
    /// The cylinder axis is too close to the plane normal for a unique rim point.
    UnsupportedFaceContact,
    /// A horizontal cylinder has a lowest supporting line rather than a unique point.
    UnsupportedLineContact,
    /// The point-contact effective mass is singular or numerically indefinite.
    SingularContactMass,
    /// The caller supplied an initial state with material penetration beyond the
    /// declared numerical admission tolerance.
    InitialPenetrationExceeded { gap_m: f64, tolerance_m: f64 },
    /// A step's geometric drift exceeded the bounded position-projection tolerance.
    GeometricProjectionExceeded { gap_m: f64, tolerance_m: f64 },
    /// The coupled sticking solve met neither the declared normal nor tangent residual bound.
    ConstraintResidualExceeded {
        residual_m_per_s: f64,
        tolerance_m_per_s: f64,
    },
    /// The retained fixed-step bound was exceeded.
    StepBudgetExceeded,
    /// A lower-level checked rigid-body operation refused the state.
    RigidBodyRefusal { detail: String },
    /// The declared dry friction interface or law refused a query.
    DryLawRefusal { detail: String },
    /// The coarse and refined runs did not retain the same terminal class.
    IncomparableRefinement,
}

impl fmt::Display for ContactDynamicsError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for ContactDynamicsError {}

/// Homogeneous finite-cylinder geometry and mass.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DiscGeometry {
    /// Outer rim radius in metres.
    pub radius_m: f64,
    /// Axial thickness in metres.
    pub thickness_m: f64,
    /// Total mass in kilograms.
    pub mass_kg: f64,
}

impl DiscGeometry {
    /// Builds the center-of-mass inertia of a homogeneous solid cylinder.
    pub fn mass_properties(self) -> Result<MassProperties, ContactDynamicsError> {
        positive_finite(self.radius_m, "geometry.radius_m")?;
        positive_finite(self.thickness_m, "geometry.thickness_m")?;
        positive_finite(self.mass_kg, "geometry.mass_kg")?;
        let radius_squared = checked_mul(self.radius_m, self.radius_m, "radius_squared")?;
        let thickness_squared =
            checked_mul(self.thickness_m, self.thickness_m, "thickness_squared")?;
        let transverse = checked_mul(
            self.mass_kg,
            (3.0 * radius_squared + thickness_squared) / 12.0,
            "transverse_inertia",
        )?;
        let axial = checked_mul(self.mass_kg, radius_squared / 2.0, "axial_inertia")?;
        MassProperties::new(
            self.mass_kg,
            Vec3::ZERO,
            Vec3::new(transverse, transverse, axial),
        )
        .map_err(|error| ContactDynamicsError::RigidBodyRefusal {
            detail: error.to_string(),
        })
    }
}

/// One declared horizontal-plane contact setup.
#[derive(Debug, Clone, PartialEq)]
pub struct ContactDynamicsInput {
    /// Homogeneous finite-cylinder geometry used to derive mass and inertia.
    pub geometry: DiscGeometry,
    /// Full center-of-mass pose and world/body momenta at the initial instant.
    pub initial_state: RigidBodyState,
    /// Positive magnitude of downward gravitational acceleration in m/s².
    pub gravity_m_per_s2: f64,
    /// Coulomb static friction coefficient, retained as caller input.
    pub static_friction_coefficient: f64,
    /// Checked dry-interface provenance retained by `fs-tribo`.
    pub interface: InterfaceSystemRef,
    /// Fixed deterministic step duration in seconds.
    pub timestep_s: f64,
    /// Number of fixed steps to attempt, up to `MAX_STEPS`.
    pub maximum_steps: u32,
    /// Largest per-step geometric drift corrected by the reported normal
    /// position projection, in metres.
    pub contact_tolerance_m: f64,
    /// Largest initially admitted negative contact gap in metres.  Larger
    /// penetration is refused instead of being Baumgarte-corrected.
    pub maximum_initial_penetration_m: f64,
    /// Separating normal-speed threshold in m/s.
    pub release_speed_tolerance_m_per_s: f64,
}

impl ContactDynamicsInput {
    /// Validates every caller-controlled scalar and reconstructs mass properties.
    pub fn validate(&self) -> Result<MassProperties, ContactDynamicsError> {
        let mass_properties = self.geometry.mass_properties()?;
        positive_finite(self.gravity_m_per_s2, "gravity_m_per_s2")?;
        nonnegative_finite(
            self.static_friction_coefficient,
            "static_friction_coefficient",
        )?;
        positive_finite(self.timestep_s, "timestep_s")?;
        if self.maximum_steps == 0 || self.maximum_steps > MAX_STEPS {
            return Err(ContactDynamicsError::StepBudgetExceeded);
        }
        positive_finite(self.contact_tolerance_m, "contact_tolerance_m")?;
        nonnegative_finite(
            self.maximum_initial_penetration_m,
            "maximum_initial_penetration_m",
        )?;
        nonnegative_finite(
            self.release_speed_tolerance_m_per_s,
            "release_speed_tolerance_m_per_s",
        )?;
        if self.interface.ordered_system_id().trim().is_empty()
            || self.interface.history_id().trim().is_empty()
            || self.interface.provenance().source_id().trim().is_empty()
        {
            return Err(ContactDynamicsError::InvalidInput {
                field: "interface_identity",
            });
        }
        // `RigidBodyState` has private fields and a checked constructor.  Its
        // accessors are still checked here before any derived arithmetic.
        finite_vec(
            self.initial_state.pose().position_world(),
            "initial_state.position_world",
        )?;
        finite_vec(
            self.initial_state.linear_momentum_world(),
            "initial_state.linear_momentum_world",
        )?;
        finite_vec(
            self.initial_state.angular_momentum_body(),
            "initial_state.angular_momentum_body",
        )?;
        Ok(mass_properties)
    }
}

/// Geometry-derived unilateral contact data at one pose.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContactGeometry {
    /// Signed gap from the rim point to the plane, positive when separated [m].
    pub gap_m: f64,
    /// Center-of-mass-to-contact vector in world coordinates [m].
    pub radius_world_m: Vec3,
    /// Contact point in world coordinates [m].
    pub point_world_m: Vec3,
    /// Unit disc symmetry axis in world coordinates.
    pub symmetry_axis_world: Vec3,
}

/// Static-friction feasibility evaluated from the solved contact impulse.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StickFeasibility {
    /// Solved non-negative normal impulse [N s].
    pub normal_impulse_ns: f64,
    /// Required tangent impulse magnitude [N s].
    pub required_tangential_impulse_ns: f64,
    /// Coulomb static impulse capacity [N s].
    pub static_capacity_impulse_ns: f64,
    /// `static_capacity_impulse_ns - required_tangential_impulse_ns` [N s].
    pub friction_cone_margin_ns: f64,
    /// Whether the demanded tangent impulse lies inside the static cone.
    pub feasible: bool,
    /// Authority retained by the dry-law response; never upgraded here.
    pub input_authority: InputAuthority,
}

/// A measured, explicit mechanical-energy ledger for one contact step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnergyLedger {
    /// Translational plus rotational plus gravitational potential energy before [J].
    pub mechanical_energy_before_j: f64,
    /// Translational plus rotational plus gravitational potential energy after [J].
    pub mechanical_energy_after_j: f64,
    /// `after - before` [J].
    pub mechanical_energy_delta_j: f64,
    /// Gravity force dotted with the center-of-mass displacement [J].
    pub gravity_work_j: f64,
    /// Midpoint estimate of the solved contact-impulse work [J].
    pub contact_impulse_work_estimate_j: f64,
    /// Mechanical-energy change from the bounded numerical position projection [J].
    ///
    /// This is reported separately because it is a geometric correction, not a
    /// physical contact-work claim.
    pub geometric_projection_work_j: f64,
    /// Normal-position projection displacement used to remove geometric drift [m].
    pub normal_position_projection_m: f64,
    /// Potential-energy change caused by that reported position projection [J].
    pub projection_potential_shift_j: f64,
    /// `mechanical_delta - contact_impulse_work_estimate - geometric_projection_work` [J].
    ///
    /// This exposes the discrete defect; it is not a conservation certificate.
    pub mechanical_balance_residual_j: f64,
}

/// One completed sticking-contact time step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContactStepReceipt {
    /// State before gravity, gyroscopic drift, and contact impulses.
    pub state_before: RigidBodyState,
    /// State after contact and reported normal position projection.
    pub state_after: RigidBodyState,
    /// Geometry at the state before the step.
    pub contact_before: ContactGeometry,
    /// Geometry after orientation and normal position projection.
    pub contact_after: ContactGeometry,
    /// Solved non-negative normal impulse [N s].
    pub normal_impulse_ns: f64,
    /// Solved tangent impulse in world coordinates [N s].
    pub tangential_impulse_world_ns: Vec3,
    /// Static-friction feasibility, which was true for a completed step.
    pub stick: StickFeasibility,
    /// Mechanical-energy quantities measured rather than assumed conserved.
    pub energy: EnergyLedger,
    /// Contact-point velocity after the static impulse before pose advancement [m/s].
    pub post_impulse_contact_velocity_world_m_per_s: Vec3,
    /// Difference between the achieved and requested contact velocity [m/s].
    pub post_impulse_contact_velocity_residual_world_m_per_s: Vec3,
}

/// A terminal event returned instead of fabricating a continuation.
#[derive(Debug, Clone, PartialEq)]
pub enum ContactTermination {
    /// The requested fixed horizon was completed while sticking remained feasible.
    HorizonReached,
    /// The state was already separated or the free normal motion was releasing.
    ContactLost {
        /// Zero-based attempted step index.
        step_index: u32,
        /// Geometry-derived plane gap [m].
        gap_m: f64,
        /// Free/present normal contact speed [m/s], positive away from the plane.
        normal_velocity_m_per_s: f64,
    },
    /// The coupled sticking impulse demanded a negative normal reaction while
    /// the configuration was not already separated.  The rung refuses to
    /// relabel that infeasible constrained solve as a contact-loss trajectory.
    UnilateralReactionInfeasible {
        /// Zero-based attempted step index.
        step_index: u32,
        /// Coupled sticking normal impulse [N s].
        required_normal_impulse_ns: f64,
    },
    /// The required tangent impulse exceeded Coulomb static capacity.
    StickInfeasible {
        /// Zero-based attempted step index.
        step_index: u32,
        /// Checked static-friction diagnostics.
        stick: StickFeasibility,
    },
}

/// A finite deterministic run of the sticking-only rung.
#[derive(Debug, Clone, PartialEq)]
pub struct ContactDynamicsRun {
    /// Every admitted sticking step, in deterministic order.
    pub steps: Vec<ContactStepReceipt>,
    /// The state at the final admitted step, or the initial state on immediate termination.
    pub final_state: RigidBodyState,
    /// Why no further step was taken.
    pub termination: ContactTermination,
}

/// Comparison of a fixed horizon with the same horizon at half the timestep.
#[derive(Debug, Clone, PartialEq)]
pub struct TimestepRefinement {
    /// Result using the caller timestep.
    pub coarse: ContactDynamicsRun,
    /// Result using half the caller timestep and twice the requested steps.
    pub fine: ContactDynamicsRun,
    /// Same horizon at one quarter of the caller timestep, used as a deterministic reference.
    pub reference: ContactDynamicsRun,
    /// Euclidean center-of-mass endpoint difference [m].
    pub final_position_difference_m: f64,
    /// Euclidean world linear-momentum endpoint difference [kg m/s].
    pub final_linear_momentum_difference_kg_m_per_s: f64,
    /// Euclidean body angular-momentum endpoint difference [kg m²/s].
    pub final_angular_momentum_difference_kg_m2_per_s: f64,
    /// Fine-minus-coarse final mechanical energy [J].
    pub final_mechanical_energy_difference_j: f64,
    /// Coarse center-of-mass endpoint distance to the quarter-step reference [m].
    pub coarse_reference_position_difference_m: f64,
    /// Half-step center-of-mass endpoint distance to the quarter-step reference [m].
    pub fine_reference_position_difference_m: f64,
    /// Whether the half-step endpoint was no farther from the retained reference.
    pub refinement_improved: bool,
}

/// Returns the unique lowest rim contact of a tilted finite cylinder.
pub fn contact_geometry(
    geometry: DiscGeometry,
    pose: Pose,
) -> Result<ContactGeometry, ContactDynamicsError> {
    geometry.mass_properties()?;
    let axis = pose
        .orientation()
        .rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
    finite_vec(axis, "symmetry_axis_world")?;
    let axis_normal_component = checked_dot(axis, GROUND_NORMAL, "axis_normal_component")?;
    let plane_projection = checked_sub(
        GROUND_NORMAL,
        checked_scale(axis, axis_normal_component, "axis_projection")?,
        "axis_plane_projection",
    )?;
    let projection_norm = stable_norm(plane_projection, "axis_plane_projection")?;
    if projection_norm <= GEOMETRY_EPSILON {
        return Err(ContactDynamicsError::UnsupportedFaceContact);
    }
    if axis_normal_component.abs() <= GEOMETRY_EPSILON {
        return Err(ContactDynamicsError::UnsupportedLineContact);
    }
    let rim_down = checked_scale(
        plane_projection,
        -geometry.radius_m / projection_norm,
        "rim_down",
    )?;
    let lower_face_sign = if axis_normal_component >= 0.0 {
        -1.0
    } else {
        1.0
    };
    let lower_face = checked_scale(
        axis,
        lower_face_sign * 0.5 * geometry.thickness_m,
        "lower_face",
    )?;
    let radius_world_m = checked_add(lower_face, rim_down, "radius_world_m")?;
    let point_world_m = checked_add(pose.position_world(), radius_world_m, "point_world_m")?;
    let gap_m = point_world_m.z;
    finite_scalar(gap_m, "gap_m")?;
    Ok(ContactGeometry {
        gap_m,
        radius_world_m,
        point_world_m,
        symmetry_axis_world: axis,
    })
}

/// Places a supplied orientation and momentum state exactly on the plane.
///
/// This helper only resolves finite-cylinder geometry; it does not impose a
/// rolling velocity, precession law, or contact force.
pub fn state_at_ground_contact(
    geometry: DiscGeometry,
    orientation: UnitQuaternion,
    linear_momentum_world: Vec3,
    angular_momentum_body: Vec3,
) -> Result<RigidBodyState, ContactDynamicsError> {
    let provisional = Pose::new(Vec3::ZERO, orientation).map_err(rigid_refusal)?;
    let provisional_geometry = contact_geometry(geometry, provisional)?;
    let position = Vec3::new(0.0, 0.0, -provisional_geometry.radius_world_m.z);
    let pose = Pose::new(position, orientation).map_err(rigid_refusal)?;
    RigidBodyState::new(pose, linear_momentum_world, angular_momentum_body).map_err(rigid_refusal)
}

/// Runs the fixed-step unilateral contact solver until its horizon or an honest terminal event.
pub fn run_contact_dynamics(
    input: &ContactDynamicsInput,
) -> Result<ContactDynamicsRun, ContactDynamicsError> {
    let mass_properties = input.validate()?;
    let mut state = input.initial_state;
    let mut steps = Vec::with_capacity(input.maximum_steps as usize);

    for step_index in 0..input.maximum_steps {
        let geometry = contact_geometry(input.geometry, state.pose())?;
        if step_index == 0 && geometry.gap_m < -input.maximum_initial_penetration_m {
            return Err(ContactDynamicsError::InitialPenetrationExceeded {
                gap_m: geometry.gap_m,
                tolerance_m: input.maximum_initial_penetration_m,
            });
        }
        if geometry.gap_m > input.contact_tolerance_m {
            return Ok(ContactDynamicsRun {
                steps,
                final_state: state,
                termination: ContactTermination::ContactLost {
                    step_index,
                    gap_m: geometry.gap_m,
                    normal_velocity_m_per_s: 0.0,
                },
            });
        }

        let attempted = sticking_step(input, mass_properties, state, geometry)?;
        match attempted {
            AttemptedStep::Completed(receipt) => {
                state = receipt.state_after;
                steps.push(receipt);
            }
            AttemptedStep::ContactLost {
                gap_m,
                normal_velocity_m_per_s,
            } => {
                return Ok(ContactDynamicsRun {
                    steps,
                    final_state: state,
                    termination: ContactTermination::ContactLost {
                        step_index,
                        gap_m,
                        normal_velocity_m_per_s,
                    },
                });
            }
            AttemptedStep::StickInfeasible(stick) => {
                return Ok(ContactDynamicsRun {
                    steps,
                    final_state: state,
                    termination: ContactTermination::StickInfeasible { step_index, stick },
                });
            }
            AttemptedStep::UnilateralReactionInfeasible {
                required_normal_impulse_ns,
            } => {
                return Ok(ContactDynamicsRun {
                    steps,
                    final_state: state,
                    termination: ContactTermination::UnilateralReactionInfeasible {
                        step_index,
                        required_normal_impulse_ns,
                    },
                });
            }
        }
    }

    Ok(ContactDynamicsRun {
        steps,
        final_state: state,
        termination: ContactTermination::HorizonReached,
    })
}

/// Repeats a retained horizon at half and quarter timesteps and compares endpoints.
///
/// This is deterministic refinement evidence, not a claim of observed order.
pub fn refine_timestep_by_two(
    input: &ContactDynamicsInput,
) -> Result<TimestepRefinement, ContactDynamicsError> {
    input.validate()?;
    let fine_steps = input
        .maximum_steps
        .checked_mul(2)
        .ok_or(ContactDynamicsError::StepBudgetExceeded)?;
    let reference_steps = fine_steps
        .checked_mul(2)
        .ok_or(ContactDynamicsError::StepBudgetExceeded)?;
    if reference_steps > MAX_STEPS {
        return Err(ContactDynamicsError::StepBudgetExceeded);
    }
    let mut fine_input = input.clone();
    fine_input.timestep_s *= 0.5;
    fine_input.maximum_steps = fine_steps;
    let mut reference_input = input.clone();
    reference_input.timestep_s *= 0.25;
    reference_input.maximum_steps = reference_steps;
    let coarse = run_contact_dynamics(input)?;
    let fine = run_contact_dynamics(&fine_input)?;
    let reference = run_contact_dynamics(&reference_input)?;
    if terminal_class(&coarse.termination) != terminal_class(&fine.termination)
        || terminal_class(&coarse.termination) != terminal_class(&reference.termination)
    {
        return Err(ContactDynamicsError::IncomparableRefinement);
    }
    let final_position_difference_m = stable_norm(
        checked_sub(
            coarse.final_state.pose().position_world(),
            fine.final_state.pose().position_world(),
            "refinement.position_difference",
        )?,
        "refinement.position_difference",
    )?;
    let final_linear_momentum_difference_kg_m_per_s = stable_norm(
        checked_sub(
            coarse.final_state.linear_momentum_world(),
            fine.final_state.linear_momentum_world(),
            "refinement.linear_momentum_difference",
        )?,
        "refinement.linear_momentum_difference",
    )?;
    let final_angular_momentum_difference_kg_m2_per_s = stable_norm(
        checked_sub(
            coarse.final_state.angular_momentum_body(),
            fine.final_state.angular_momentum_body(),
            "refinement.angular_momentum_difference",
        )?,
        "refinement.angular_momentum_difference",
    )?;
    let coarse_energy = mechanical_energy(
        input.gravity_m_per_s2,
        input.geometry.mass_kg,
        input.geometry.mass_properties()?,
        coarse.final_state,
    )?;
    let fine_energy = mechanical_energy(
        input.gravity_m_per_s2,
        input.geometry.mass_kg,
        input.geometry.mass_properties()?,
        fine.final_state,
    )?;
    let final_mechanical_energy_difference_j =
        checked_scalar_sub(fine_energy, coarse_energy, "refinement.energy_difference")?;
    let coarse_reference_position_difference_m = endpoint_position_difference(
        coarse.final_state,
        reference.final_state,
        "refinement.coarse_reference_position_difference",
    )?;
    let fine_reference_position_difference_m = endpoint_position_difference(
        fine.final_state,
        reference.final_state,
        "refinement.fine_reference_position_difference",
    )?;
    Ok(TimestepRefinement {
        coarse,
        fine,
        reference,
        final_position_difference_m,
        final_linear_momentum_difference_kg_m_per_s,
        final_angular_momentum_difference_kg_m2_per_s,
        final_mechanical_energy_difference_j,
        coarse_reference_position_difference_m,
        fine_reference_position_difference_m,
        refinement_improved: fine_reference_position_difference_m
            <= coarse_reference_position_difference_m,
    })
}

enum AttemptedStep {
    Completed(ContactStepReceipt),
    ContactLost {
        gap_m: f64,
        normal_velocity_m_per_s: f64,
    },
    StickInfeasible(StickFeasibility),
    UnilateralReactionInfeasible {
        required_normal_impulse_ns: f64,
    },
}

fn sticking_step(
    input: &ContactDynamicsInput,
    mass_properties: MassProperties,
    state: RigidBodyState,
    contact_before: ContactGeometry,
) -> Result<AttemptedStep, ContactDynamicsError> {
    let duration = input.timestep_s;
    let mass = mass_properties.mass();
    let momentum = state.linear_momentum_world();
    let angular_momentum = state.angular_momentum_body();
    let angular_velocity = mass_properties.angular_velocity_body(angular_momentum);
    finite_vec(angular_velocity, "angular_velocity_body")?;
    let gravity_impulse = Vec3::new(0.0, 0.0, -mass * input.gravity_m_per_s2 * duration);
    finite_vec(gravity_impulse, "gravity_impulse")?;
    let free_momentum = checked_add(momentum, gravity_impulse, "free_momentum")?;
    // Euler's equation in principal body coordinates: dL/dt = L x omega.
    let free_angular_momentum = checked_add(
        angular_momentum,
        checked_scale(
            checked_cross(angular_momentum, angular_velocity, "gyroscopic_torque")?,
            duration,
            "gyroscopic_impulse",
        )?,
        "free_angular_momentum",
    )?;
    let free_state = RigidBodyState::new(state.pose(), free_momentum, free_angular_momentum)
        .map_err(rigid_refusal)?;
    let free_contact_velocity =
        contact_velocity(mass_properties, free_state, contact_before.radius_world_m)?;
    let target_normal_velocity = (-contact_before.gap_m * 0.2 / duration).max(0.0);
    finite_scalar(target_normal_velocity, "target_normal_velocity")?;
    if free_contact_velocity.z > input.release_speed_tolerance_m_per_s
        && contact_before.gap_m >= -input.maximum_initial_penetration_m
    {
        return Ok(AttemptedStep::ContactLost {
            gap_m: contact_before.gap_m,
            normal_velocity_m_per_s: free_contact_velocity.z,
        });
    }
    let target_contact_velocity = Vec3::new(0.0, 0.0, target_normal_velocity);
    let total_impulse = coupled_sticking_impulse(
        mass_properties,
        state.pose().orientation(),
        contact_before.radius_world_m,
        free_contact_velocity,
        target_contact_velocity,
    )?;
    let normal_impulse_ns = total_impulse.z;
    if normal_impulse_ns < 0.0 {
        return Ok(AttemptedStep::UnilateralReactionInfeasible {
            required_normal_impulse_ns: normal_impulse_ns,
        });
    }
    let tangent_impulse = Vec3::new(total_impulse.x, total_impulse.y, 0.0);
    let required_tangential_impulse_ns = stable_norm(tangent_impulse, "tangent_impulse")?;
    let normal_reaction_n = checked_div(normal_impulse_ns, duration, "normal_reaction_n")?;
    let frame = ContactFrame::new([0.0, 0.0, 1.0]).map_err(tribo_refusal)?;
    let zero_slip = TangentialSlip::new(&frame, [0.0, 0.0, 0.0]).map_err(tribo_refusal)?;
    let friction = FrictionLaw::Coulomb {
        static_mu: input.static_friction_coefficient,
        kinetic_mu: input.static_friction_coefficient,
    }
    .evaluate(&input.interface, normal_reaction_n, zero_slip)
    .map_err(tribo_refusal)?;
    if friction.regime != FrictionRegime::Sticking {
        return Err(ContactDynamicsError::DryLawRefusal {
            detail: "zero-slip Coulomb query did not report sticking".to_owned(),
        });
    }
    let static_capacity_impulse_ns = checked_mul(
        friction.static_limit,
        duration,
        "static_capacity_impulse_ns",
    )?;
    let friction_cone_margin_ns = checked_scalar_sub(
        static_capacity_impulse_ns,
        required_tangential_impulse_ns,
        "friction_cone_margin_ns",
    )?;
    let feasible = friction_cone_margin_ns
        >= -contact_impulse_tolerance(static_capacity_impulse_ns, required_tangential_impulse_ns);
    let stick = StickFeasibility {
        normal_impulse_ns,
        required_tangential_impulse_ns,
        static_capacity_impulse_ns,
        friction_cone_margin_ns,
        feasible,
        input_authority: friction.provenance().authority(),
    };
    if !feasible {
        return Ok(AttemptedStep::StickInfeasible(stick));
    }
    let post_impulse = apply_impulse(
        mass_properties,
        free_state,
        contact_before.radius_world_m,
        total_impulse,
    )?;
    let post_impulse_contact_velocity =
        contact_velocity(mass_properties, post_impulse, contact_before.radius_world_m)?;
    let post_impulse_contact_velocity_residual = checked_sub(
        post_impulse_contact_velocity,
        target_contact_velocity,
        "post_impulse_contact_velocity_residual",
    )?;
    let residual_norm = stable_norm(
        post_impulse_contact_velocity_residual,
        "post_impulse_contact_velocity_residual",
    )?;
    let residual_tolerance = contact_velocity_tolerance(target_contact_velocity)?;
    if residual_norm > residual_tolerance {
        return Err(ContactDynamicsError::ConstraintResidualExceeded {
            residual_m_per_s: residual_norm,
            tolerance_m_per_s: residual_tolerance,
        });
    }
    let post_angular_velocity =
        mass_properties.angular_velocity_body(post_impulse.angular_momentum_body());
    finite_vec(post_angular_velocity, "post_angular_velocity_body")?;
    let orientation = state
        .pose()
        .orientation()
        .right_exp(checked_scale(
            post_angular_velocity,
            duration,
            "rotation_vector_body",
        )?)
        .map_err(rigid_refusal)?;
    let raw_position = checked_add(
        state.pose().position_world(),
        checked_scale(
            post_impulse.linear_momentum_world(),
            duration / mass,
            "center_of_mass_displacement",
        )?,
        "raw_position",
    )?;
    let raw_pose = Pose::new(raw_position, orientation).map_err(rigid_refusal)?;
    let raw_contact = contact_geometry(input.geometry, raw_pose)?;
    if raw_contact.gap_m.abs() > input.contact_tolerance_m {
        return Err(ContactDynamicsError::GeometricProjectionExceeded {
            gap_m: raw_contact.gap_m,
            tolerance_m: input.contact_tolerance_m,
        });
    }
    let corrected_position = checked_sub(
        raw_position,
        Vec3::new(0.0, 0.0, raw_contact.gap_m),
        "normal_position_projection",
    )?;
    let position_projection_m = -raw_contact.gap_m;
    let pose_after = Pose::new(corrected_position, orientation).map_err(rigid_refusal)?;
    let state_after = RigidBodyState::new(
        pose_after,
        post_impulse.linear_momentum_world(),
        post_impulse.angular_momentum_body(),
    )
    .map_err(rigid_refusal)?;
    let raw_state = RigidBodyState::new(
        raw_pose,
        post_impulse.linear_momentum_world(),
        post_impulse.angular_momentum_body(),
    )
    .map_err(rigid_refusal)?;
    let contact_after = contact_geometry(input.geometry, pose_after)?;
    if contact_after.gap_m.abs() > input.contact_tolerance_m {
        return Err(ContactDynamicsError::NonFiniteDerived {
            field: "projected_contact_gap_m",
        });
    }
    let energy = energy_ledger(
        input.gravity_m_per_s2,
        mass_properties,
        state,
        raw_state,
        state_after,
        total_impulse,
        free_contact_velocity,
        post_impulse_contact_velocity,
        position_projection_m,
    )?;
    Ok(AttemptedStep::Completed(ContactStepReceipt {
        state_before: state,
        state_after,
        contact_before,
        contact_after,
        normal_impulse_ns,
        tangential_impulse_world_ns: tangent_impulse,
        stick,
        energy,
        post_impulse_contact_velocity_world_m_per_s: post_impulse_contact_velocity,
        post_impulse_contact_velocity_residual_world_m_per_s:
            post_impulse_contact_velocity_residual,
    }))
}

fn coupled_sticking_impulse(
    mass_properties: MassProperties,
    orientation: UnitQuaternion,
    radius_world_m: Vec3,
    free_contact_velocity: Vec3,
    target_contact_velocity: Vec3,
) -> Result<Vec3, ContactDynamicsError> {
    let directions = [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        GROUND_NORMAL,
    ];
    let responses = [
        contact_velocity_delta(mass_properties, orientation, radius_world_m, directions[0])?,
        contact_velocity_delta(mass_properties, orientation, radius_world_m, directions[1])?,
        contact_velocity_delta(mass_properties, orientation, radius_world_m, directions[2])?,
    ];
    let mut system = [[0.0; 4]; 3];
    for row in 0..3 {
        for column in 0..3 {
            system[row][column] =
                checked_dot(directions[row], responses[column], "coupled_contact_mass")?;
        }
    }
    let target_delta = checked_sub(
        target_contact_velocity,
        free_contact_velocity,
        "coupled_contact_velocity_delta",
    )?;
    system[0][3] = target_delta.x;
    system[1][3] = target_delta.y;
    system[2][3] = target_delta.z;
    let solution = solve_3x3(system)?;
    Ok(Vec3::new(solution[0], solution[1], solution[2]))
}

fn solve_3x3(mut system: [[f64; 4]; 3]) -> Result<[f64; 3], ContactDynamicsError> {
    let mut scale: f64 = 0.0;
    for row in system {
        for value in row.into_iter().take(3) {
            finite_scalar(value, "coupled_contact_mass")?;
            scale = scale.max(value.abs());
        }
    }
    let pivot_tolerance = 256.0 * f64::EPSILON * scale.max(1.0);
    for column in 0..3 {
        let mut pivot_row = column;
        for candidate in (column + 1)..3 {
            if system[candidate][column].abs() > system[pivot_row][column].abs() {
                pivot_row = candidate;
            }
        }
        if system[pivot_row][column].abs() <= pivot_tolerance {
            return Err(ContactDynamicsError::SingularContactMass);
        }
        if pivot_row != column {
            system.swap(pivot_row, column);
        }
        let pivot = system[column][column];
        for row in (column + 1)..3 {
            let factor = checked_div(system[row][column], pivot, "coupled_contact_elimination")?;
            for entry in column..4 {
                system[row][entry] = checked_scalar_sub(
                    system[row][entry],
                    checked_mul(factor, system[column][entry], "coupled_contact_elimination")?,
                    "coupled_contact_elimination",
                )?;
            }
        }
    }
    let z = checked_div(system[2][3], system[2][2], "coupled_contact_solution")?;
    let y = checked_div(
        checked_scalar_sub(
            system[1][3],
            checked_mul(system[1][2], z, "coupled_contact_solution")?,
            "coupled_contact_solution",
        )?,
        system[1][1],
        "coupled_contact_solution",
    )?;
    let x = checked_div(
        checked_scalar_sub(
            checked_scalar_sub(
                system[0][3],
                checked_mul(system[0][1], y, "coupled_contact_solution")?,
                "coupled_contact_solution",
            )?,
            checked_mul(system[0][2], z, "coupled_contact_solution")?,
            "coupled_contact_solution",
        )?,
        system[0][0],
        "coupled_contact_solution",
    )?;
    Ok([x, y, z])
}

fn apply_impulse(
    mass_properties: MassProperties,
    state: RigidBodyState,
    radius_world_m: Vec3,
    impulse_world_ns: Vec3,
) -> Result<RigidBodyState, ContactDynamicsError> {
    let linear = checked_add(
        state.linear_momentum_world(),
        impulse_world_ns,
        "post_impulse_linear_momentum",
    )?;
    let radius_body = world_to_body(state.pose().orientation(), radius_world_m)?;
    let impulse_body = world_to_body(state.pose().orientation(), impulse_world_ns)?;
    let angular = checked_add(
        state.angular_momentum_body(),
        checked_cross(radius_body, impulse_body, "contact_impulse_torque")?,
        "post_impulse_angular_momentum",
    )?;
    RigidBodyState::new(state.pose(), linear, angular).map_err(rigid_refusal)
}

fn contact_velocity_delta(
    mass_properties: MassProperties,
    orientation: UnitQuaternion,
    radius_world_m: Vec3,
    impulse_world: Vec3,
) -> Result<Vec3, ContactDynamicsError> {
    let impulse_body = world_to_body(orientation, impulse_world)?;
    let radius_body = world_to_body(orientation, radius_world_m)?;
    let angular_impulse_body = checked_cross(radius_body, impulse_body, "angular_impulse_body")?;
    let delta_angular_velocity_body = mass_properties.angular_velocity_body(angular_impulse_body);
    finite_vec(delta_angular_velocity_body, "delta_angular_velocity_body")?;
    let delta_angular_velocity_world =
        orientation.rotate_body_to_world(delta_angular_velocity_body);
    let rotational = checked_cross(
        delta_angular_velocity_world,
        radius_world_m,
        "contact_rotational_delta",
    )?;
    let translational = checked_scale(
        impulse_world,
        1.0 / mass_properties.mass(),
        "contact_translational_delta",
    )?;
    checked_add(translational, rotational, "contact_velocity_delta")
}

fn contact_velocity(
    mass_properties: MassProperties,
    state: RigidBodyState,
    radius_world_m: Vec3,
) -> Result<Vec3, ContactDynamicsError> {
    let linear = checked_scale(
        state.linear_momentum_world(),
        1.0 / mass_properties.mass(),
        "center_of_mass_velocity",
    )?;
    let angular_body = mass_properties.angular_velocity_body(state.angular_momentum_body());
    finite_vec(angular_body, "angular_velocity_body")?;
    let angular_world = state
        .pose()
        .orientation()
        .rotate_body_to_world(angular_body);
    checked_add(
        linear,
        checked_cross(angular_world, radius_world_m, "contact_angular_velocity")?,
        "contact_velocity",
    )
}

fn mechanical_energy(
    gravity_m_per_s2: f64,
    mass_kg: f64,
    mass_properties: MassProperties,
    state: RigidBodyState,
) -> Result<f64, ContactDynamicsError> {
    let momentum_norm = stable_norm(state.linear_momentum_world(), "linear_momentum")?;
    let translational = checked_div(
        checked_mul(momentum_norm, momentum_norm, "translational_energy")?,
        2.0 * mass_kg,
        "translational_energy",
    )?;
    let angular_velocity = mass_properties.angular_velocity_body(state.angular_momentum_body());
    let rotational = 0.5
        * checked_dot(
            state.angular_momentum_body(),
            angular_velocity,
            "rotational_energy",
        )?;
    let potential = checked_mul(
        mass_kg * gravity_m_per_s2,
        state.pose().position_world().z,
        "gravitational_potential",
    )?;
    checked_scalar_add(
        checked_scalar_add(translational, rotational, "mechanical_energy")?,
        potential,
        "mechanical_energy",
    )
}

#[allow(clippy::too_many_arguments)]
fn energy_ledger(
    gravity_m_per_s2: f64,
    mass_properties: MassProperties,
    state_before: RigidBodyState,
    state_before_projection: RigidBodyState,
    state_after: RigidBodyState,
    total_impulse: Vec3,
    contact_velocity_before: Vec3,
    contact_velocity_after: Vec3,
    normal_position_projection_m: f64,
) -> Result<EnergyLedger, ContactDynamicsError> {
    let mass = mass_properties.mass();
    let mechanical_energy_before_j =
        mechanical_energy(gravity_m_per_s2, mass, mass_properties, state_before)?;
    let mechanical_energy_before_projection_j = mechanical_energy(
        gravity_m_per_s2,
        mass,
        mass_properties,
        state_before_projection,
    )?;
    let mechanical_energy_after_j =
        mechanical_energy(gravity_m_per_s2, mass, mass_properties, state_after)?;
    let mechanical_energy_delta_j = checked_scalar_sub(
        mechanical_energy_after_j,
        mechanical_energy_before_j,
        "mechanical_energy_delta",
    )?;
    let displacement = checked_sub(
        state_after.pose().position_world(),
        state_before.pose().position_world(),
        "gravity_displacement",
    )?;
    let gravity_work_j = checked_mul(-mass * gravity_m_per_s2, displacement.z, "gravity_work")?;
    let midpoint_contact_velocity = checked_scale(
        checked_add(
            contact_velocity_before,
            contact_velocity_after,
            "midpoint_contact_velocity",
        )?,
        0.5,
        "midpoint_contact_velocity",
    )?;
    let contact_impulse_work_estimate_j = checked_dot(
        total_impulse,
        midpoint_contact_velocity,
        "contact_impulse_work",
    )?;
    let projection_potential_shift_j = checked_mul(
        mass * gravity_m_per_s2,
        normal_position_projection_m,
        "projection_potential_shift",
    )?;
    let geometric_projection_work_j = checked_scalar_sub(
        mechanical_energy_after_j,
        mechanical_energy_before_projection_j,
        "geometric_projection_work",
    )?;
    let mechanical_balance_residual_j = checked_scalar_sub(
        checked_scalar_sub(
            mechanical_energy_delta_j,
            contact_impulse_work_estimate_j,
            "mechanical_balance_residual",
        )?,
        geometric_projection_work_j,
        "mechanical_balance_residual",
    )?;
    Ok(EnergyLedger {
        mechanical_energy_before_j,
        mechanical_energy_after_j,
        mechanical_energy_delta_j,
        gravity_work_j,
        contact_impulse_work_estimate_j,
        geometric_projection_work_j,
        normal_position_projection_m,
        projection_potential_shift_j,
        mechanical_balance_residual_j,
    })
}

fn world_to_body(
    orientation: UnitQuaternion,
    vector_world: Vec3,
) -> Result<Vec3, ContactDynamicsError> {
    finite_vec(vector_world, "vector_world")?;
    let [w, x, y, z] = orientation.components();
    let result = Vec3::new(
        (1.0 - 2.0 * (y * y + z * z)).mul_add(
            vector_world.x,
            (2.0 * (x * y + w * z)).mul_add(vector_world.y, 2.0 * (x * z - w * y) * vector_world.z),
        ),
        (2.0 * (x * y - w * z)).mul_add(
            vector_world.x,
            (1.0 - 2.0 * (x * x + z * z))
                .mul_add(vector_world.y, 2.0 * (y * z + w * x) * vector_world.z),
        ),
        (2.0 * (x * z + w * y)).mul_add(
            vector_world.x,
            (2.0 * (y * z - w * x)).mul_add(
                vector_world.y,
                (1.0 - 2.0 * (x * x + y * y)) * vector_world.z,
            ),
        ),
    );
    finite_vec(result, "vector_body")?;
    Ok(result)
}

fn terminal_class(termination: &ContactTermination) -> u8 {
    match termination {
        ContactTermination::HorizonReached => 0,
        ContactTermination::ContactLost { .. } => 1,
        ContactTermination::StickInfeasible { .. } => 2,
        ContactTermination::UnilateralReactionInfeasible { .. } => 3,
    }
}

fn endpoint_position_difference(
    left: RigidBodyState,
    right: RigidBodyState,
    field: &'static str,
) -> Result<f64, ContactDynamicsError> {
    stable_norm(
        checked_sub(
            left.pose().position_world(),
            right.pose().position_world(),
            field,
        )?,
        field,
    )
}

fn contact_impulse_tolerance(capacity: f64, required: f64) -> f64 {
    256.0 * f64::EPSILON * capacity.abs().max(required.abs())
}

fn contact_velocity_tolerance(target: Vec3) -> Result<f64, ContactDynamicsError> {
    checked_mul(
        1024.0 * f64::EPSILON,
        stable_norm(target, "target_contact_velocity")?.max(1.0),
        "contact_velocity_tolerance",
    )
}

fn positive_finite(value: f64, field: &'static str) -> Result<(), ContactDynamicsError> {
    if !(value.is_finite() && value > 0.0) {
        return Err(ContactDynamicsError::InvalidInput { field });
    }
    Ok(())
}

fn nonnegative_finite(value: f64, field: &'static str) -> Result<(), ContactDynamicsError> {
    if !(value.is_finite() && value >= 0.0) {
        return Err(ContactDynamicsError::InvalidInput { field });
    }
    Ok(())
}

fn finite_scalar(value: f64, field: &'static str) -> Result<(), ContactDynamicsError> {
    if !value.is_finite() {
        return Err(ContactDynamicsError::NonFiniteDerived { field });
    }
    Ok(())
}

fn finite_vec(value: Vec3, field: &'static str) -> Result<(), ContactDynamicsError> {
    if !value.is_finite() {
        return Err(ContactDynamicsError::NonFiniteDerived { field });
    }
    Ok(())
}

fn stable_norm(value: Vec3, field: &'static str) -> Result<f64, ContactDynamicsError> {
    finite_vec(value, field)?;
    let scale = value.x.abs().max(value.y.abs()).max(value.z.abs());
    if scale == 0.0 {
        return Ok(0.0);
    }
    let scaled = Vec3::new(value.x / scale, value.y / scale, value.z / scale);
    let result = scale * scaled.dot(scaled).sqrt();
    finite_scalar(result, field)?;
    Ok(result)
}

fn checked_add(left: Vec3, right: Vec3, field: &'static str) -> Result<Vec3, ContactDynamicsError> {
    let value = left.add(right);
    finite_vec(value, field)?;
    Ok(value)
}

fn checked_sub(left: Vec3, right: Vec3, field: &'static str) -> Result<Vec3, ContactDynamicsError> {
    let value = left.sub(right);
    finite_vec(value, field)?;
    Ok(value)
}

fn checked_scalar_sub(
    left: f64,
    right: f64,
    field: &'static str,
) -> Result<f64, ContactDynamicsError> {
    let value = left - right;
    finite_scalar(value, field)?;
    Ok(value)
}

fn checked_scalar_add(
    left: f64,
    right: f64,
    field: &'static str,
) -> Result<f64, ContactDynamicsError> {
    let value = left + right;
    finite_scalar(value, field)?;
    Ok(value)
}

fn checked_scale(
    value: Vec3,
    scalar: f64,
    field: &'static str,
) -> Result<Vec3, ContactDynamicsError> {
    finite_scalar(scalar, field)?;
    let result = value.scale(scalar);
    finite_vec(result, field)?;
    Ok(result)
}

fn checked_cross(
    left: Vec3,
    right: Vec3,
    field: &'static str,
) -> Result<Vec3, ContactDynamicsError> {
    let value = left.cross(right);
    finite_vec(value, field)?;
    Ok(value)
}

fn checked_dot(left: Vec3, right: Vec3, field: &'static str) -> Result<f64, ContactDynamicsError> {
    let value = left.dot(right);
    finite_scalar(value, field)?;
    Ok(value)
}

fn checked_mul(left: f64, right: f64, field: &'static str) -> Result<f64, ContactDynamicsError> {
    let value = left * right;
    finite_scalar(value, field)?;
    Ok(value)
}

fn checked_div(left: f64, right: f64, field: &'static str) -> Result<f64, ContactDynamicsError> {
    if right == 0.0 || !right.is_finite() {
        return Err(ContactDynamicsError::NonFiniteDerived { field });
    }
    let value = left / right;
    finite_scalar(value, field)?;
    Ok(value)
}

fn rigid_refusal(error: fs_mbd::DynamicsError) -> ContactDynamicsError {
    ContactDynamicsError::RigidBodyRefusal {
        detail: error.to_string(),
    }
}

fn tribo_refusal(error: fs_tribo::TriboError) -> ContactDynamicsError {
    ContactDynamicsError::DryLawRefusal {
        detail: error.to_string(),
    }
}
