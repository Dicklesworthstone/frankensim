//! First-order force/contact/drift for loaded free articulated trees.
//!
//! Loads kick the original coupled inertia, one simultaneous Coulomb solve
//! supplies contact reactions, and the native unforced Lie midpoint advances
//! every base and joint. Geometry must admit the complete candidate path; the
//! independent endpoint gap checks are not continuous collision detection.

use super::friction::{
    ArticulatedFrictionalContact, ArticulatedFrictionalContactResult,
    resolve_articulated_frictional_contacts,
};
use super::trajectory::{ArticulatedTrajectoryMomentum, ArticulatedTrajectorySystem};
use super::{
    ArticulatedContactError, ArticulatedContactSystem, ContactCounterpart, LinkContactPoint,
};
use crate::articulated::motion::{
    FreeFloatingMotionError, FreeFloatingMotionLoads, FreeFloatingMotionState,
    FreeFloatingMotionStep, MAX_FREE_FLOATING_STEP_LINKS, step_free_floating,
};
use crate::articulated::{
    ArticulatedError, ArticulatedModel, ArticulatedMomentum, BaseState,
    FreeFloatingImpulseResponse, Kinematics, apply_free_floating_impulses, forward_kinematics,
    free_floating_momentum, validate_external_wrenches, validate_forward_dynamics_header,
};
use fs_contact::frictional_multicontact::FrictionalContactConfig;
use fs_ga::{Se3, Twist, Vec3, Wrench};

/// Independent physical gates for the composed step; SI units throughout.
#[derive(Debug, Clone, Copy)]
pub struct ArticulatedContactStepConfig {
    /// Finite positive duration, in seconds.
    pub duration: f64,
    /// Fixed world point about which aggregate angular momenta are reported.
    pub angular_origin_world: Vec3,
    /// Original simultaneous contact solver controls.
    pub contact: FrictionalContactConfig,
    /// Maximum absolute per-tree and aggregate unforced drift energy change, in J.
    pub drift_energy_tolerance: f64,
    /// Maximum complete-step work/energy residual, in joules.
    pub energy_tolerance: f64,
    /// Maximum norm of complete-step linear impulse residual, in N s.
    pub linear_momentum_tolerance: f64,
    /// Maximum norm of complete-step angular impulse residual, in N m s.
    pub angular_momentum_tolerance: f64,
}

/// Actual signed gap at an original material witness, in stable-ID order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ArticulatedContactGap {
    /// Original stable contact identity.
    pub id: u64,
    /// Normal separation in metres; negative values always refuse.
    pub signed_gap: f64,
}

/// Original candidate numerical path supplied to the required geometry owner.
///
/// The geometry owner must establish feature coverage and applicability of the
/// supplied frozen normals/witnesses throughout the interval. A finite set of
/// samples alone does not establish that claim. The callback must not publish
/// external history before the complete step succeeds.
pub struct ArticulatedContactPath<'a> {
    /// Original borrowed models, pre-step states, and held loads.
    pub systems: &'a [ArticulatedTrajectorySystem<'a>],
    /// Coupled states immediately after all force and contact impulses.
    pub post_contact: &'a [FreeFloatingMotionState],
    /// Actual candidate native endpoints.
    pub endpoint: &'a [FreeFloatingMotionState],
    /// Material contact features at the step start; planes translate once.
    pub contacts: &'a [ArticulatedFrictionalContact],
    /// Duration of the original unforced drift, in seconds.
    pub duration: f64,
}
impl ArticulatedContactPath<'_> {
    /// Sample the original zero-load Lie drift at a fraction in [0,1].
    /// This is a numerical path query, not a no-missed-collision certificate.
    pub fn sample(
        &self,
        fraction: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<FreeFloatingMotionState>, FreeFloatingMotionError> {
        if cancelled() {
            return Err(FreeFloatingMotionError::Cancelled);
        }
        if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
            return Err(FreeFloatingMotionError::InvalidDuration);
        }
        if fraction == 0.0 {
            return Ok(self.post_contact.to_vec());
        }
        if fraction == 1.0 {
            return Ok(self.endpoint.to_vec());
        }
        self.systems
            .iter()
            .zip(self.post_contact)
            .map(|(system, state)| {
                Ok(drift(
                    system.model,
                    state,
                    fraction * self.duration,
                    &mut cancelled,
                )?
                .state_after)
            })
            .collect()
    }
}

/// Whole-step refusal; all caller states remain borrowed and unchanged.
#[derive(Debug)]
pub enum ArticulatedContactStepError<E> {
    /// Invalid bounds, configuration payload, or finite-value admission.
    Invalid(&'static str),
    /// Original tree-state, force-limit, or impulse response refused.
    Dynamics {
        system_index: usize,
        error: ArticulatedError,
    },
    /// Original unforced drift refused.
    Drift {
        system_index: usize,
        error: FreeFloatingMotionError,
    },
    /// Original simultaneous contact owner refused.
    Contact(ArticulatedContactError),
    /// Required whole-path geometry admission refused.
    Path(E),
    /// A supplied material contact penetrates at an actual endpoint.
    Penetration {
        contact_id: u64,
        endpoint: bool,
        signed_gap: f64,
    },
    /// A physical balance exceeded its independent caller-owned allowance.
    Balance(&'static str),
    /// Cancellation before publication.
    Cancelled,
}
impl<E: core::fmt::Debug> core::fmt::Display for ArticulatedContactStepError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "articulated contact step: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for ArticulatedContactStepError<E> {}

/// Accepted physical state, native responses, and separated work/impulse terms.
pub struct ArticulatedContactStepResult<P> {
    /// Reusable physical endpoints in input-system order.
    pub states_after: Vec<FreeFloatingMotionState>,
    /// Native force/gravity/joint impulse responses before contact.
    pub load_responses: Vec<FreeFloatingImpulseResponse>,
    /// Native simultaneous contact response, retaining its determinacy scope.
    /// Its loss includes removal of the free load kick by support reactions;
    /// in sustained contact that splitting loss is not continuum friction heat.
    pub contact: ArticulatedFrictionalContactResult,
    /// Caller geometry receipt, published only with all accepted states.
    pub path_admission: P,
    /// Actual initial witness gaps.
    pub initial_gaps: Vec<ArticulatedContactGap>,
    /// Actual endpoint gaps, with prescribed planes advanced once.
    pub endpoint_gaps: Vec<ArticulatedContactGap>,
    /// Independent physical initial momentum and energy.
    pub momentum_before: ArticulatedTrajectoryMomentum,
    /// Independent physical final momentum and energy.
    pub momentum_after: ArticulatedTrajectoryMomentum,
    /// Total load and support impulse; dynamic contact pairs cancel globally.
    pub external_linear_impulse_world: Vec3,
    /// Total applied angular impulse about the configured fixed origin.
    pub external_angular_impulse_about_origin_world: Vec3,
    /// Exact joint impulse work at average pre/post-kick joint rates, in J.
    pub actuator_impulse_work: f64,
    /// Exact gravity impulse work at average pre/post-kick link velocities, in J.
    pub gravity_impulse_work: f64,
    /// Exact link-frame wrench impulse work, in J. Its continuous path work is
    /// not inferred from endpoints when the link frame rotates.
    pub body_wrench_impulse_work: f64,
    /// Exact work of the held generalized effort along actual joint displacement.
    pub actuator_displacement_work: f64,
    /// Exact uniform-gravity work along every actual link COM displacement.
    pub gravity_displacement_work: f64,
    /// Joint kick work minus actual displacement work: a first-order splitting
    /// defect, never physical contact heat.
    pub actuator_splitting_work: f64,
    /// Gravity kick work minus actual displacement work, also numerical splitting.
    pub gravity_splitting_work: f64,
    /// Initial gravity-kick angular impulse minus the exact gravity angular
    /// integral over the accepted constant-total-momentum COM drift, in N m s.
    pub gravity_angular_splitting_impulse_world: Vec3,
    /// Actual kinetic change during the zero-load native drift, in joules.
    pub drift_energy_defect: f64,
    /// Drift angular-momentum defect about the configured origin, in N m s.
    pub drift_angular_momentum_defect: Vec3,
    /// Actual total kinetic change minus kick/support work plus native impact loss.
    /// This includes the measured drift error and does not declare a heat partition.
    pub kinetic_work_residual: f64,
    /// Actual linear change minus all applied impulses, in N s.
    pub linear_momentum_residual: Vec3,
    /// Actual angular change minus all applied angular impulses, in N m s.
    pub angular_momentum_residual: Vec3,
}

/// Advance loaded articulated systems through one sustained contact step.
///
/// Internal joint forces retain equal-and-opposite base recoil. Gravity acts at
/// every original link COM. All contact rows are solved simultaneously through
/// original ABA mobility; links are never replaced by independent lumped masses.
/// Resting support ordinarily uses zero restitution/target. Declared positive
/// restitution remains an impact law and is never silently changed.
///
/// Every supplied initial/endpoint gap must be nonnegative. The mandatory path
/// callback must admit omitted-feature coverage, intermediate separation, and
/// continued applicability of the frozen witnesses and normals. No configuration
/// or momentum projection repairs a failure. All native states stay private until
/// the callback, physical gates, and final cancellation check succeed.
///
/// This first-order force/contact/drift composition is not general collision
/// discovery, a complete Moreau--Jean method, or an exact continuous-force flow.
/// Refine the step for force/drift splitting and nonlinear articulated motion.
/// At most 128 systems, 4096 total links, and 21 contacts are admitted. Native
/// tree passes are O(links); only the bounded contact-space matrix is dense.
#[allow(clippy::too_many_lines)]
pub fn step_articulated_frictional_contacts<P, E>(
    systems: &[ArticulatedTrajectorySystem<'_>],
    contacts: &[ArticulatedFrictionalContact],
    config: ArticulatedContactStepConfig,
    admit_path: impl FnOnce(&ArticulatedContactPath<'_>) -> Result<P, E>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedContactStepResult<P>, ArticulatedContactStepError<E>> {
    poll(&mut cancelled)?;
    admission(systems, contacts, config)?;
    let initial: Vec<_> = systems.iter().map(|system| system.state.clone()).collect();
    let before = kinematics(systems, &initial, &mut cancelled)?;
    let initial_gaps = gaps(&before, contacts, 0.0, &mut cancelled)?;
    let momentum_before = totals(
        systems,
        &initial,
        config.angular_origin_world,
        &mut cancelled,
    )?;
    let mut post_contact = initial.clone();
    let mut load_responses = Vec::with_capacity(systems.len());
    let mut actuator_work = 0.0;
    let mut gravity_work = 0.0;
    let mut body_work = 0.0;
    let mut linear_impulse = Vec3::default();
    let mut angular_impulse = Vec3::default();
    for (system_index, system) in systems.iter().enumerate() {
        poll(&mut cancelled)?;
        let loads = system.loads;
        let state = system.state;
        validate_forward_dynamics_header(
            system.model,
            state.base.twist_body,
            &state.joint_position,
            &state.joint_velocity,
            loads.generalized_force,
            loads.gravity_world,
            loads.external_body_wrench,
        )
        .map_err(|error| dynamics(system_index, error))?;
        validate_external_wrenches(loads.external_body_wrench)
            .map_err(|error| dynamics(system_index, error))?;
        let joint_impulses: Vec<_> = loads
            .generalized_force
            .iter()
            .map(|force| config.duration * force)
            .collect();
        let mut body_impulses = Vec::with_capacity(system.model.link_count());
        let mut gravity_impulses = Vec::with_capacity(system.model.link_count());
        for (link_index, link) in system.model.links().iter().enumerate() {
            if link_index % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let gravity_force = before[system_index].world_from_link[link_index]
                .rotation()
                .inverse()
                .rotate(loads.gravity_world)
                .map_err(|error| dynamics(system_index, error.into()))?
                .scale(link.inertia().mass() * config.duration);
            let gravity = Wrench::new(
                link.inertia().center_of_mass().cross(gravity_force),
                gravity_force,
            );
            let applied = loads.external_body_wrench[link_index];
            body_impulses.push(Wrench::new(
                gravity.torque + applied.torque.scale(config.duration),
                gravity.force + applied.force.scale(config.duration),
            ));
            gravity_impulses.push(gravity);
        }
        let response = apply_free_floating_impulses(
            system.model,
            state.base,
            &state.joint_position,
            &state.joint_velocity,
            &joint_impulses,
            &body_impulses,
        )
        .map_err(|error| dynamics(system_index, error))?;
        post_contact[system_index].base = response.base_after;
        post_contact[system_index]
            .joint_velocity
            .clone_from(&response.joint_velocity_after);
        let after_kick = forward_kinematics(
            system.model,
            BaseState::prescribed(
                response.base_after.world_from_base,
                response.base_after.twist_body,
                Twist::zero(),
            ),
            &state.joint_position,
            &response.joint_velocity_after,
        )
        .map_err(|error| dynamics(system_index, error))?;
        for (link_index, gravity) in gravity_impulses.iter().enumerate() {
            if link_index % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let midpoint = before[system_index].body_twist[link_index]
                .scale(0.5)
                .plus(after_kick.body_twist[link_index].scale(0.5));
            gravity_work += gravity.pairing(midpoint);
            body_work += config.duration * loads.external_body_wrench[link_index].pairing(midpoint);
        }
        for (index, impulse) in joint_impulses.iter().enumerate() {
            actuator_work += impulse
                * (0.5 * state.joint_velocity[index] + 0.5 * response.joint_velocity_after[index]);
        }
        linear_impulse = linear_impulse + response.external_linear_impulse_world;
        angular_impulse = angular_impulse + response.external_angular_impulse_world
            - config
                .angular_origin_world
                .cross(response.external_linear_impulse_world);
        load_responses.push(response);
    }
    let inputs: Vec<_> = systems
        .iter()
        .zip(&post_contact)
        .map(|(system, state)| ArticulatedContactSystem {
            model: system.model,
            base: state.base,
            joint_position: &state.joint_position,
            joint_velocity: &state.joint_velocity,
        })
        .collect();
    let contact =
        resolve_articulated_frictional_contacts(&inputs, contacts, config.contact, &mut cancelled)
            .map_err(|error| {
                if matches!(error, ArticulatedContactError::Cancelled) {
                    ArticulatedContactStepError::Cancelled
                } else {
                    ArticulatedContactStepError::Contact(error)
                }
            })?;
    for (state, response) in post_contact.iter_mut().zip(&contact.responses) {
        state.base = response.base_after;
        state
            .joint_velocity
            .clone_from(&response.joint_velocity_after);
        linear_impulse = linear_impulse + response.external_linear_impulse_world;
        angular_impulse = angular_impulse + response.external_angular_impulse_world
            - config
                .angular_origin_world
                .cross(response.external_linear_impulse_world);
    }
    let before_drift = totals(
        systems,
        &post_contact,
        config.angular_origin_world,
        &mut cancelled,
    )?;
    let mut states_after = Vec::with_capacity(systems.len());
    let mut drift_energy_defect = 0.0;
    for (index, (system, state)) in systems.iter().zip(&post_contact).enumerate() {
        let result =
            drift(system.model, state, config.duration, &mut cancelled).map_err(|error| {
                if matches!(error, FreeFloatingMotionError::Cancelled) {
                    ArticulatedContactStepError::Cancelled
                } else {
                    ArticulatedContactStepError::Drift {
                        system_index: index,
                        error,
                    }
                }
            })?;
        if result.kinetic_work_residual.abs() > config.drift_energy_tolerance {
            return Err(ArticulatedContactStepError::Balance(
                "one tree exceeds the unforced drift energy gate",
            ));
        }
        drift_energy_defect += result.kinetic_work_residual;
        states_after.push(result.state_after);
    }
    let after = kinematics(systems, &states_after, &mut cancelled)?;
    let endpoint_gaps = gaps(&after, contacts, config.duration, &mut cancelled)?;
    let momentum_after = totals(
        systems,
        &states_after,
        config.angular_origin_world,
        &mut cancelled,
    )?;
    let mut actuator_displacement_work = 0.0;
    let mut gravity_displacement_work = 0.0;
    let mut gravity_angular_splitting = Vec3::default();
    for (system_index, system) in systems.iter().enumerate() {
        for (dof, &effort) in system.loads.generalized_force.iter().enumerate() {
            actuator_displacement_work += effort
                * (states_after[system_index].joint_position[dof]
                    - system.state.joint_position[dof]);
        }
        for (link_index, link) in system.model.links().iter().enumerate() {
            if link_index % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let local_com = link.inertia().center_of_mass();
            let first = before[system_index].world_from_link[link_index]
                .transform_point(local_com)
                .map_err(|error| dynamics(system_index, error.into()))?;
            let last = after[system_index].world_from_link[link_index]
                .transform_point(local_com)
                .map_err(|error| dynamics(system_index, error.into()))?;
            let force = system.loads.gravity_world.scale(link.inertia().mass());
            gravity_displacement_work += force.dot(last - first);
            // Individual link COM paths need not be linear. Their mass-weighted
            // sum is the native unforced tree COM and IS linear, so this sum
            // gives the exact total uniform-gravity angular splitting term.
            gravity_angular_splitting = gravity_angular_splitting
                + (first - last).cross(force).scale(0.5 * config.duration);
        }
    }
    let work_residual = momentum_after.kinetic_energy
        - momentum_before.kinetic_energy
        - actuator_work
        - gravity_work
        - body_work
        - contact.plane_work
        + contact.contact.dissipated_energy;
    let linear_residual =
        momentum_after.linear_world - momentum_before.linear_world - linear_impulse;
    let angular_residual = momentum_after.angular_about_origin_world
        - momentum_before.angular_about_origin_world
        - angular_impulse;
    let drift_angular =
        momentum_after.angular_about_origin_world - before_drift.angular_about_origin_world;
    let actuator_splitting_work = actuator_work - actuator_displacement_work;
    let gravity_splitting_work = gravity_work - gravity_displacement_work;
    if ![
        actuator_work,
        gravity_work,
        body_work,
        actuator_displacement_work,
        gravity_displacement_work,
        work_residual,
        drift_energy_defect,
        actuator_splitting_work,
        gravity_splitting_work,
    ]
    .iter()
    .all(|x| x.is_finite())
        || ![
            linear_impulse,
            angular_impulse,
            linear_residual,
            angular_residual,
            drift_angular,
            gravity_angular_splitting,
        ]
        .iter()
        .all(|x| finite_vector(*x))
    {
        return Err(ArticulatedContactStepError::Invalid(
            "nonfinite step balance",
        ));
    }
    if drift_energy_defect.abs() > config.drift_energy_tolerance
        || work_residual.abs() > config.energy_tolerance
        || linear_residual.norm() > config.linear_momentum_tolerance
        || angular_residual.norm() > config.angular_momentum_tolerance
    {
        return Err(ArticulatedContactStepError::Balance(
            "native step exceeds work or momentum gate",
        ));
    }
    poll(&mut cancelled)?;
    let path_admission = admit_path(&ArticulatedContactPath {
        systems,
        post_contact: &post_contact,
        endpoint: &states_after,
        contacts,
        duration: config.duration,
    })
    .map_err(ArticulatedContactStepError::Path)?;
    poll(&mut cancelled)?;
    Ok(ArticulatedContactStepResult {
        states_after,
        load_responses,
        contact,
        path_admission,
        initial_gaps,
        endpoint_gaps,
        momentum_before,
        momentum_after,
        external_linear_impulse_world: linear_impulse,
        external_angular_impulse_about_origin_world: angular_impulse,
        actuator_impulse_work: actuator_work,
        gravity_impulse_work: gravity_work,
        body_wrench_impulse_work: body_work,
        actuator_displacement_work,
        gravity_displacement_work,
        actuator_splitting_work,
        gravity_splitting_work,
        gravity_angular_splitting_impulse_world: gravity_angular_splitting,
        drift_energy_defect,
        drift_angular_momentum_defect: drift_angular,
        kinetic_work_residual: work_residual,
        linear_momentum_residual: linear_residual,
        angular_momentum_residual: angular_residual,
    })
}

fn admission<E>(
    systems: &[ArticulatedTrajectorySystem<'_>],
    contacts: &[ArticulatedFrictionalContact],
    config: ArticulatedContactStepConfig,
) -> Result<(), ArticulatedContactStepError<E>> {
    if systems.is_empty()
        || systems.len() > 128
        || contacts.len() > 21
        || contacts.len() > config.contact.maximum_contacts
        || !finite_vector(config.angular_origin_world)
        || [
            config.duration,
            config.drift_energy_tolerance,
            config.energy_tolerance,
            config.linear_momentum_tolerance,
            config.angular_momentum_tolerance,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err(ArticulatedContactStepError::Invalid(
            "bounded systems/contacts and positive finite gates required",
        ));
    }
    let mut links = 0;
    for system in systems {
        links += system.model.link_count();
        if links > MAX_FREE_FLOATING_STEP_LINKS {
            return Err(ArticulatedContactStepError::Invalid(
                "total link cap exceeded",
            ));
        }
    }
    Ok(())
}
fn drift(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    duration: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<FreeFloatingMotionStep, FreeFloatingMotionError> {
    let efforts = vec![0.0; model.dof_count()];
    let wrenches = vec![Wrench::default(); model.link_count()];
    step_free_floating(
        model,
        state,
        FreeFloatingMotionLoads {
            generalized_force: &efforts,
            external_body_wrench: &wrenches,
            gravity_world: Vec3::default(),
        },
        duration,
        cancelled,
    )
}
fn kinematics<E>(
    systems: &[ArticulatedTrajectorySystem<'_>],
    states: &[FreeFloatingMotionState],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<Kinematics>, ArticulatedContactStepError<E>> {
    systems
        .iter()
        .zip(states)
        .enumerate()
        .map(|(index, (system, state))| {
            poll(cancelled)?;
            forward_kinematics(
                system.model,
                BaseState::prescribed(
                    state.base.world_from_base,
                    state.base.twist_body,
                    Twist::zero(),
                ),
                &state.joint_position,
                &state.joint_velocity,
            )
            .map_err(|error| dynamics(index, error))
        })
        .collect()
}
fn totals<E>(
    systems: &[ArticulatedTrajectorySystem<'_>],
    states: &[FreeFloatingMotionState],
    origin: Vec3,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ArticulatedTrajectoryMomentum, ArticulatedContactStepError<E>> {
    let mut total = ArticulatedTrajectoryMomentum {
        linear_world: Vec3::default(),
        angular_about_origin_world: Vec3::default(),
        kinetic_energy: 0.0,
    };
    for (index, (system, state)) in systems.iter().zip(states).enumerate() {
        poll(cancelled)?;
        let mut base = state.base;
        base.world_from_base = Se3::from_parts(
            base.world_from_base.rotation(),
            base.world_from_base.translation() - origin,
        )
        .map_err(|error| dynamics(index, error.into()))?;
        let value: ArticulatedMomentum = free_floating_momentum(
            system.model,
            base,
            &state.joint_position,
            &state.joint_velocity,
        )
        .map_err(|error| dynamics(index, error))?;
        total.linear_world = total.linear_world + value.linear_momentum_world;
        total.angular_about_origin_world =
            total.angular_about_origin_world + value.angular_momentum_world;
        total.kinetic_energy += value.kinetic_energy;
    }
    Ok(total)
}
fn gaps<E>(
    states: &[Kinematics],
    contacts: &[ArticulatedFrictionalContact],
    time: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<ArticulatedContactGap>, ArticulatedContactStepError<E>> {
    let mut result = Vec::with_capacity(contacts.len());
    for contact in contacts {
        poll(cancelled)?;
        let point = position(states, contact.point)?;
        let other = match contact.other {
            ContactCounterpart::Body(other) => position(states, other)?,
            ContactCounterpart::PrescribedPlane {
                point_world,
                velocity_world,
            } => point_world + velocity_world.scale(time),
        };
        let normal = super::unit_normal(contact.normal_world)
            .map_err(ArticulatedContactStepError::Contact)?;
        let gap = normal.dot(point - other);
        if !gap.is_finite() {
            return Err(ArticulatedContactStepError::Invalid(
                "nonfinite witness gap",
            ));
        }
        if gap < 0.0 {
            return Err(ArticulatedContactStepError::Penetration {
                contact_id: contact.id,
                endpoint: time != 0.0,
                signed_gap: gap,
            });
        }
        result.push(ArticulatedContactGap {
            id: contact.id,
            signed_gap: gap,
        });
    }
    result.sort_unstable_by_key(|gap| gap.id);
    Ok(result)
}
fn position<E>(
    states: &[Kinematics],
    point: LinkContactPoint,
) -> Result<Vec3, ArticulatedContactStepError<E>> {
    let pose = states
        .get(point.system_index)
        .and_then(|state| state.world_from_link.get(point.link_index))
        .ok_or(ArticulatedContactStepError::Invalid(
            "contact link index outside supplied systems",
        ))?;
    pose.transform_point(point.point_body)
        .map_err(|error| dynamics(point.system_index, error.into()))
}
fn finite_vector(value: Vec3) -> bool {
    value.x.is_finite() && value.y.is_finite() && value.z.is_finite()
}
fn dynamics<E>(system_index: usize, error: ArticulatedError) -> ArticulatedContactStepError<E> {
    ArticulatedContactStepError::Dynamics {
        system_index,
        error,
    }
}
fn poll<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), ArticulatedContactStepError<E>> {
    if cancelled() {
        Err(ArticulatedContactStepError::Cancelled)
    } else {
        Ok(())
    }
}
