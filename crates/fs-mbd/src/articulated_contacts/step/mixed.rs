//! Loaded sustained contact for original free and prescribed articulated trees.
//!
//! A held-load impulse, one simultaneous native Coulomb event, and an unloaded
//! native smooth step compose a first-order force/contact/drift interval. A
//! prescribed boundary remains a physical source during all three stages.

use super::{
    ArticulatedContactGap, ArticulatedContactStepConfig, ArticulatedContactStepError,
    com_displacement, drift, dynamics, finite_vector, gaps, impulse_pairing, poll,
    relative_kinematics,
};
use crate::articulated::motion::prescribed::{
    PrescribedBaseMotionError, PrescribedBaseMotionState, PrescribedBaseMotionStep,
    PrescribedBasePath, step_prescribed_base,
};
use crate::articulated::motion::{
    FreeFloatingMotionError, FreeFloatingMotionLoads, FreeFloatingMotionStep,
    MAX_FREE_FLOATING_STEP_LINKS,
};
use crate::articulated::{
    ArticulatedModel, BaseState, Kinematics, forward_kinematics, momentum_from_kinematics,
    validate_external_wrenches, validate_forward_dynamics_header,
};
use crate::articulated_contacts::ArticulatedContactError;
use crate::articulated_contacts::friction::{
    ArticulatedFrictionalContact, ArticulatedFrictionalContactResult,
};
use crate::articulated_contacts::mixed::{
    ArticulatedContactBoundary, ContactResponse, ContactSystem, MixedArticulatedContactSystem,
    MixedArticulatedImpulseResponse, resolve_mixed_articulated_frictional_contacts,
};
use crate::articulated_contacts::trajectory::mixed::{
    MixedArticulatedTrajectoryState, MixedArticulatedTrajectorySystem,
};
use crate::articulated_contacts::trajectory::{ArticulatedTrajectoryMomentum, Impulses};
use fs_contact::frictional_multicontact::FrictionalContactResult;
use fs_ga::{Se3, Vec3, Wrench};
use fs_math::dd::Dd;

/// Actual native unloaded drift, retaining its independent source and defect terms.
#[derive(Debug, Clone)]
pub enum MixedArticulatedContactDrift {
    /// Original free Lie midpoint drift; no physical external load is supplied.
    Free(FreeFloatingMotionStep),
    /// Original prescribed-base drift. A driven support can still supply work
    /// and spatial impulse even though gravity, body loads and efforts are zero.
    Prescribed(PrescribedBaseMotionStep),
}

/// Complete candidate supplied to the mandatory whole-path geometry owner.
///
/// The owner must establish intermediate separation, omitted-feature coverage,
/// and continued applicability of the frozen contact witnesses and normals for
/// the original native drift and the caller's actual smooth support paths.
/// Midpoint/end samples alone do not establish that claim. No interpolation of
/// a prescribed support path is invented here. Callback-owned history must be
/// staged until the entire returned step is accepted.
pub struct MixedArticulatedContactPath<'a> {
    /// Original borrowed models, accepted states and held physical loads.
    pub systems: &'a [MixedArticulatedTrajectorySystem<'a>],
    /// Coupled states immediately after the load kick and simultaneous contact.
    pub post_contact: &'a [MixedArticulatedTrajectoryState],
    /// Actual candidate native endpoints, including prescribed support clocks.
    pub endpoint: &'a [MixedArticulatedTrajectoryState],
    /// Actual caller-supplied samples, one per system; free systems have None.
    pub support_paths: &'a [Option<PrescribedBasePath>],
    /// Original material contact features; planes translate once over the step.
    pub contacts: &'a [ArticulatedFrictionalContact],
    /// Common initial physical clock, in seconds.
    pub time_before_s: f64,
    /// Duration of the original native unloaded drift, in seconds.
    pub duration: f64,
}

/// Accepted states and native evidence for one mixed sustained-contact step.
pub struct MixedArticulatedContactStepResult<P> {
    /// Common initial physical clock, in seconds.
    pub time_before_s: f64,
    /// Common accepted physical clock, in seconds.
    pub time_after_s: f64,
    /// Reusable native endpoints in original system order.
    pub states_after: Vec<MixedArticulatedTrajectoryState>,
    /// Original load kick responses, including full prescribed support reactions.
    pub load_responses: Vec<MixedArticulatedImpulseResponse>,
    /// Original simultaneous law, including work determinacy and support sources.
    /// Its loss includes removal of the free load kick by support reactions; that
    /// force/contact splitting loss is numerical, not continuum friction heat.
    pub contact: ArticulatedFrictionalContactResult<
        FrictionalContactResult,
        MixedArticulatedImpulseResponse,
    >,
    /// Original native drift receipts, including independent support quadrature.
    pub drift_responses: Vec<MixedArticulatedContactDrift>,
    /// Actual supplied support-path samples; free systems have None.
    pub support_paths: Vec<Option<PrescribedBasePath>>,
    /// Geometry owner's admission of the complete candidate path.
    pub path_admission: P,
    /// Actual signed gaps at the original witnesses before the step.
    pub initial_gaps: Vec<ArticulatedContactGap>,
    /// Actual endpoint gaps, with translating planes advanced exactly once.
    pub endpoint_gaps: Vec<ArticulatedContactGap>,
    /// Physical native initial totals, retained as independent diagnostics.
    pub momentum_before: ArticulatedTrajectoryMomentum,
    /// Physical native endpoint totals, retained as independent diagnostics.
    pub momentum_after: ArticulatedTrajectoryMomentum,
    /// Actual per-stage native kinetic changes, reduced before rounding, J.
    pub kinetic_energy_change: f64,
    /// Actual per-stage native linear changes, in world axes, N s.
    pub linear_momentum_change_world: Vec3,
    /// Actual angular changes about the configured fixed origin, N m s.
    pub angular_momentum_change_about_origin_world: Vec3,
    /// Body/gravity, plane and prescribed support impulses ON all retained trees.
    /// Native equal-and-opposite body contact impulses cancel in the reduction.
    pub external_linear_impulse_world: Vec3,
    /// Corresponding angular impulse about the configured fixed origin, N m s.
    pub external_angular_impulse_about_origin_world: Vec3,
    /// Relative joint impulses paired with actual average pre/post-kick rates, J.
    pub actuator_impulse_work: f64,
    /// Gravity impulses paired with actual average pre/post-kick link twists, J.
    pub gravity_impulse_work: f64,
    /// Link-frame wrench impulse work, J. Continuous work of a rotating body
    /// wrench cannot generally be inferred from its two endpoint poses.
    pub body_wrench_impulse_work: f64,
    /// Native prescribed-base support work during the held-load kick, J.
    pub load_support_work: f64,
    /// Native prescribed-base support work during the unloaded drift, J.
    pub drift_support_work: f64,
    /// Held relative effort paired with the actual joint displacement, J.
    pub actuator_displacement_work: f64,
    /// Uniform-gravity work along actual link COM displacements, J.
    pub gravity_displacement_work: f64,
    /// Joint kick work minus displacement work; numerical splitting, never heat.
    pub actuator_splitting_work: f64,
    /// Gravity kick work minus displacement work; numerical splitting, never heat.
    pub gravity_splitting_work: f64,
    /// Native drift kinetic change minus driven support work, J. A support's
    /// physical work is not an integration defect or dissipation.
    pub drift_energy_defect: f64,
    /// Native drift actual linear change minus independent support impulse, N s.
    pub drift_linear_momentum_defect: Vec3,
    /// Corresponding angular defect about the configured fixed origin, N m s.
    pub drift_angular_momentum_defect: Vec3,
    /// Actual kinetic change minus kick/boundary/drift work plus contact loss, J.
    pub kinetic_work_residual: f64,
    /// Actual linear change minus every original external/support impulse, N s.
    pub linear_momentum_residual: Vec3,
    /// Corresponding angular impulse residual about the fixed origin, N m s.
    pub angular_momentum_residual: Vec3,
}

/// Advance supported and free mechanisms through one loaded contact interval.
///
/// The original free/prescribed ABA inertia operators apply held forces once as
/// impulses over the duration. One simultaneous native Coulomb solve then uses
/// their actual coupled mobility. Original native smooth owners drift with zero
/// held forces; prescribed support paths still produce their physical reaction
/// impulse and work. Root-joint actuator transmission and all base recoil remain
/// native, without fictitious support masses or linkwise mass lumping.
///
/// `support_path` is called once per prescribed system, with its unchanged
/// accepted base/clock and post-contact joint rates. It must supply midpoint and
/// endpoint samples of the same smooth path as that initial base. The complete
/// path must be admitted by `admit_path`; derivative consistency and geometry
/// cannot be inferred from these two samples. Both callbacks must stage external
/// history until this function succeeds. Inputs remain borrowed on any refusal.
///
/// Initial and endpoint gaps must be nonnegative; no pose or momentum projection
/// repairs a failure. The existing absolute work/momentum gates apply, and the
/// drift gate limits kinetic change minus driven-support work. Per-stage actual
/// changes and source terms are reduced in double-double precision. Absolute
/// endpoint totals are diagnostics, not the basis of the balance subtraction.
///
/// This first-order force/contact/drift composition is not a complete Moreau--
/// Jean scheme or contact discovery. Resting support ordinarily uses a zero
/// normal target/restitution; a caller's impact law is never silently changed.
/// No continuous angular gravity integral is inferred for curved driven paths,
/// and no continuous rotating-body-wrench work or thermal partition is invented.
/// Native stages are O(links); only the original bounded contact matrix is dense.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn step_mixed_articulated_frictional_contacts<P, E>(
    systems: &[MixedArticulatedTrajectorySystem<'_>],
    contacts: &[ArticulatedFrictionalContact],
    start_time_s: f64,
    config: ArticulatedContactStepConfig,
    mut support_path: impl FnMut(
        usize,
        &PrescribedBaseMotionState,
        f64,
    ) -> Result<PrescribedBasePath, PrescribedBaseMotionError>,
    admit_path: impl FnOnce(&MixedArticulatedContactPath<'_>) -> Result<P, E>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MixedArticulatedContactStepResult<P>, ArticulatedContactStepError<E>> {
    poll(&mut cancelled)?;
    let end_time_s = admission(systems, contacts, start_time_s, config)?;
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
    let mut changes = Impulses::default();
    let mut sources = Impulses::default();
    let mut kinetic_change = Dd::ZERO;
    let mut actuator_work = Dd::ZERO;
    let mut gravity_work = Dd::ZERO;
    let mut body_work = Dd::ZERO;
    let mut load_support_work = Dd::ZERO;
    for (system_index, system) in systems.iter().enumerate() {
        poll(&mut cancelled)?;
        let input = native_input(system.model, system.state);
        let view = input.system_view();
        let loads = system.loads;
        validate_forward_dynamics_header(
            system.model,
            view.base.twist_body,
            view.joint_position,
            view.joint_velocity,
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
        let response = input
            .apply(&joint_impulses, &body_impulses)
            .map_err(|error| dynamics(system_index, error))?;
        let after_kick = response.impulse_view();
        let after_kinematics = forward_kinematics(
            system.model,
            after_kick.base_after,
            view.joint_position,
            after_kick.joint_velocity_after,
        )
        .map_err(|error| dynamics(system_index, error))?;
        for (link_index, gravity) in gravity_impulses.iter().enumerate() {
            if link_index % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let first = before[system_index].body_twist[link_index];
            let last = after_kinematics.body_twist[link_index];
            gravity_work = gravity_work + impulse_pairing(*gravity, first, last);
            let applied = loads.external_body_wrench[link_index];
            body_work = body_work
                + impulse_pairing(
                    Wrench::new(
                        applied.torque.scale(config.duration),
                        applied.force.scale(config.duration),
                    ),
                    first,
                    last,
                );
        }
        for (dof, &impulse) in joint_impulses.iter().enumerate() {
            actuator_work = actuator_work
                + Dd::from_f64(impulse)
                    * (Dd::from_f64(0.5 * view.joint_velocity[dof])
                        + Dd::from_f64(0.5 * after_kick.joint_velocity_after[dof]));
        }
        if let Some((_, work)) = after_kick.support {
            load_support_work = load_support_work + Dd::from_f64(work);
        }
        add_response(
            &response,
            &mut changes,
            &mut sources,
            &mut kinetic_change,
            config.angular_origin_world,
        );
        apply_response(&mut post_contact[system_index], &response)?;
        load_responses.push(response);
    }
    let inputs: Vec<_> = systems
        .iter()
        .zip(&post_contact)
        .map(|(system, state)| native_input(system.model, state))
        .collect();
    let contact = resolve_mixed_articulated_frictional_contacts(
        &inputs,
        contacts,
        config.contact,
        &mut cancelled,
    )
    .map_err(|error| {
        if matches!(error, ArticulatedContactError::Cancelled) {
            ArticulatedContactStepError::Cancelled
        } else {
            ArticulatedContactStepError::Contact(error)
        }
    })?;
    for (state, response) in post_contact.iter_mut().zip(&contact.responses) {
        add_response(
            response,
            &mut changes,
            &mut sources,
            &mut kinetic_change,
            config.angular_origin_world,
        );
        apply_response(state, response)?;
    }
    let mut states_after = Vec::with_capacity(systems.len());
    let mut drift_responses = Vec::with_capacity(systems.len());
    let mut support_paths = Vec::with_capacity(systems.len());
    let mut drift_support_work = Dd::ZERO;
    let mut drift_energy_defect = Dd::ZERO;
    let mut drift_momentum_defect = Impulses::default();
    for (system_index, (system, state)) in systems.iter().zip(&post_contact).enumerate() {
        poll(&mut cancelled)?;
        let (receipt, defect) = match state {
            MixedArticulatedTrajectoryState::Free(state) => {
                let step = drift(system.model, state, config.duration, &mut cancelled).map_err(
                    |error| {
                        if matches!(error, FreeFloatingMotionError::Cancelled) {
                            ArticulatedContactStepError::Cancelled
                        } else {
                            ArticulatedContactStepError::Drift {
                                system_index,
                                error,
                            }
                        }
                    },
                )?;
                sources.add(
                    step.external_linear_impulse_world,
                    step.external_angular_impulse_world,
                    config.angular_origin_world,
                );
                changes.add(
                    step.linear_momentum_change_world,
                    step.angular_momentum_change_world,
                    config.angular_origin_world,
                );
                drift_momentum_defect.add(
                    step.linear_momentum_residual_world,
                    step.angular_momentum_residual_world,
                    config.angular_origin_world,
                );
                kinetic_change = kinetic_change + Dd::from_f64(step.kinetic_energy_change);
                states_after.push(MixedArticulatedTrajectoryState::Free(
                    step.state_after.clone(),
                ));
                support_paths.push(None);
                let defect = step.kinetic_work_residual;
                (MixedArticulatedContactDrift::Free(step), defect)
            }
            MixedArticulatedTrajectoryState::Prescribed(state) => {
                let path = support_path(system_index, state, config.duration)
                    .map_err(|error| prescribed_error(system_index, error))?;
                poll(&mut cancelled)?;
                let efforts = vec![0.0; system.model.dof_count()];
                let wrenches = vec![Wrench::default(); system.model.link_count()];
                let step = step_prescribed_base(
                    system.model,
                    state,
                    path,
                    FreeFloatingMotionLoads {
                        generalized_force: &efforts,
                        external_body_wrench: &wrenches,
                        gravity_world: Vec3::default(),
                    },
                    config.duration,
                    &mut cancelled,
                )
                .map_err(|error| prescribed_error(system_index, error))?;
                if step.state_after.time_s != end_time_s {
                    return Err(ArticulatedContactStepError::Invalid(
                        "prescribed drift changed the common clock",
                    ));
                }
                sources.add(
                    step.external_linear_impulse_world,
                    step.external_angular_impulse_world,
                    config.angular_origin_world,
                );
                sources.add(
                    step.support_linear_impulse_world,
                    step.support_angular_impulse_world,
                    config.angular_origin_world,
                );
                changes.add(
                    step.linear_momentum_change_world,
                    step.angular_momentum_change_world,
                    config.angular_origin_world,
                );
                drift_momentum_defect.add(
                    step.linear_momentum_residual_world,
                    step.angular_momentum_residual_world,
                    config.angular_origin_world,
                );
                kinetic_change = kinetic_change + Dd::from_f64(step.kinetic_energy_change);
                drift_support_work = drift_support_work + Dd::from_f64(step.support_work);
                states_after.push(MixedArticulatedTrajectoryState::Prescribed(
                    step.state_after.clone(),
                ));
                support_paths.push(Some(path));
                let defect = step.kinetic_work_residual;
                (MixedArticulatedContactDrift::Prescribed(step), defect)
            }
        };
        if !defect.is_finite() || defect.abs() > config.drift_energy_tolerance {
            return Err(ArticulatedContactStepError::Balance(
                "one tree exceeds the unloaded drift work gate",
            ));
        }
        drift_energy_defect = drift_energy_defect + Dd::from_f64(defect);
        drift_responses.push(receipt);
    }
    let after = kinematics(systems, &states_after, &mut cancelled)?;
    let endpoint_gaps = gaps(&after, contacts, config.duration, &mut cancelled)?;
    let momentum_after = totals(
        systems,
        &states_after,
        config.angular_origin_world,
        &mut cancelled,
    )?;
    let mut actuator_displacement_work = Dd::ZERO;
    let mut gravity_displacement_work = Dd::ZERO;
    for (system_index, system) in systems.iter().enumerate() {
        let first = native_input(system.model, system.state);
        let last = native_input(system.model, &states_after[system_index]);
        for (dof, &effort) in system.loads.generalized_force.iter().enumerate() {
            actuator_displacement_work = actuator_displacement_work
                + Dd::from_f64(effort)
                    * (Dd::from_f64(last.joint_position[dof])
                        - Dd::from_f64(first.joint_position[dof]));
        }
        if system.loads.gravity_world == Vec3::default() {
            continue;
        }
        poll(&mut cancelled)?;
        let first_shape = relative_kinematics(system.model, first.joint_position)
            .map_err(|error| dynamics(system_index, error))?;
        poll(&mut cancelled)?;
        let last_shape = relative_kinematics(system.model, last.joint_position)
            .map_err(|error| dynamics(system_index, error))?;
        let first_base = first.system_view().base.world_from_base;
        let last_base = last.system_view().base.world_from_base;
        for (link_index, link) in system.model.links().iter().enumerate() {
            if link_index % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let local_com = link.inertia().center_of_mass();
            let displacement = com_displacement(
                first_base,
                last_base,
                first_shape.world_from_link[link_index],
                last_shape.world_from_link[link_index],
                local_com,
            )
            .map_err(|error| dynamics(system_index, error.into()))?;
            for (displacement, gravity) in displacement.into_iter().zip([
                system.loads.gravity_world.x,
                system.loads.gravity_world.y,
                system.loads.gravity_world.z,
            ]) {
                gravity_displacement_work = gravity_displacement_work
                    + Dd::from_f64(link.inertia().mass()) * Dd::from_f64(gravity) * displacement;
            }
        }
    }
    let kinetic_work_residual = (kinetic_change
        - actuator_work
        - gravity_work
        - body_work
        - load_support_work
        - Dd::from_f64(contact.boundary_work)
        - drift_support_work
        + Dd::from_f64(contact.contact.dissipated_energy))
    .to_f64();
    let (linear_momentum_residual, angular_momentum_residual) = changes.difference(&sources);
    let actuator_splitting_work = (actuator_work - actuator_displacement_work).to_f64();
    let gravity_splitting_work = (gravity_work - gravity_displacement_work).to_f64();
    let [
        kinetic_energy_change,
        actuator_impulse_work,
        gravity_impulse_work,
        body_wrench_impulse_work,
        load_support_work,
        drift_support_work,
        actuator_displacement_work,
        gravity_displacement_work,
        drift_energy_defect,
    ] = [
        kinetic_change,
        actuator_work,
        gravity_work,
        body_work,
        load_support_work,
        drift_support_work,
        actuator_displacement_work,
        gravity_displacement_work,
        drift_energy_defect,
    ]
    .map(Dd::to_f64);
    if ![
        kinetic_energy_change,
        actuator_impulse_work,
        gravity_impulse_work,
        body_wrench_impulse_work,
        load_support_work,
        drift_support_work,
        actuator_displacement_work,
        gravity_displacement_work,
        drift_energy_defect,
        actuator_splitting_work,
        gravity_splitting_work,
        kinetic_work_residual,
    ]
    .iter()
    .all(|x| x.is_finite())
        || ![
            changes.linear,
            changes.angular,
            sources.linear,
            sources.angular,
            drift_momentum_defect.linear,
            drift_momentum_defect.angular,
            linear_momentum_residual,
            angular_momentum_residual,
        ]
        .iter()
        .all(|x| finite_vector(*x))
    {
        return Err(ArticulatedContactStepError::Invalid(
            "nonfinite mixed step balance",
        ));
    }
    if drift_energy_defect.abs() > config.drift_energy_tolerance
        || kinetic_work_residual.abs() > config.energy_tolerance
        || linear_momentum_residual.norm() > config.linear_momentum_tolerance
        || angular_momentum_residual.norm() > config.angular_momentum_tolerance
    {
        return Err(ArticulatedContactStepError::Balance(
            "native mixed step exceeds work or momentum gate",
        ));
    }
    poll(&mut cancelled)?;
    let path_admission = admit_path(&MixedArticulatedContactPath {
        systems,
        post_contact: &post_contact,
        endpoint: &states_after,
        support_paths: &support_paths,
        contacts,
        time_before_s: start_time_s,
        duration: config.duration,
    })
    .map_err(ArticulatedContactStepError::Path)?;
    poll(&mut cancelled)?;
    Ok(MixedArticulatedContactStepResult {
        time_before_s: start_time_s,
        time_after_s: end_time_s,
        states_after,
        load_responses,
        contact,
        drift_responses,
        support_paths,
        path_admission,
        initial_gaps,
        endpoint_gaps,
        momentum_before,
        momentum_after,
        kinetic_energy_change,
        linear_momentum_change_world: changes.linear,
        angular_momentum_change_about_origin_world: changes.angular,
        external_linear_impulse_world: sources.linear,
        external_angular_impulse_about_origin_world: sources.angular,
        actuator_impulse_work,
        gravity_impulse_work,
        body_wrench_impulse_work,
        load_support_work,
        drift_support_work,
        actuator_displacement_work,
        gravity_displacement_work,
        actuator_splitting_work,
        gravity_splitting_work,
        drift_energy_defect,
        drift_linear_momentum_defect: drift_momentum_defect.linear,
        drift_angular_momentum_defect: drift_momentum_defect.angular,
        kinetic_work_residual,
        linear_momentum_residual,
        angular_momentum_residual,
    })
}

fn admission<E>(
    systems: &[MixedArticulatedTrajectorySystem<'_>],
    contacts: &[ArticulatedFrictionalContact],
    time: f64,
    config: ArticulatedContactStepConfig,
) -> Result<f64, ArticulatedContactStepError<E>> {
    let end_time = time + config.duration;
    if systems.is_empty()
        || systems.len() > 128
        || contacts.len() > 21
        || contacts.len() > config.contact.maximum_contacts
        || !finite_vector(config.angular_origin_world)
        || !time.is_finite()
        || !end_time.is_finite()
        || end_time <= time
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
            "bounded systems/contacts, advancing clock and positive finite gates required",
        ));
    }
    let mut links = 0usize;
    for system in systems {
        links = links.checked_add(system.model.link_count()).ok_or(
            ArticulatedContactStepError::Invalid("total link count overflow"),
        )?;
        if links > MAX_FREE_FLOATING_STEP_LINKS {
            return Err(ArticulatedContactStepError::Invalid(
                "total link cap exceeded",
            ));
        }
        if let MixedArticulatedTrajectoryState::Prescribed(state) = system.state {
            if state.time_s != time {
                return Err(ArticulatedContactStepError::Invalid(
                    "prescribed clock differs from the common step clock",
                ));
            }
        }
    }
    Ok(end_time)
}

pub(super) fn native_input<'a>(
    model: &'a ArticulatedModel,
    state: &'a MixedArticulatedTrajectoryState,
) -> MixedArticulatedContactSystem<'a> {
    let (boundary, joint_position, joint_velocity) = match state {
        MixedArticulatedTrajectoryState::Free(state) => (
            ArticulatedContactBoundary::Free(state.base),
            &state.joint_position,
            &state.joint_velocity,
        ),
        MixedArticulatedTrajectoryState::Prescribed(state) => (
            ArticulatedContactBoundary::Prescribed(state.base),
            &state.joint_position,
            &state.joint_velocity,
        ),
    };
    MixedArticulatedContactSystem {
        model,
        boundary,
        joint_position,
        joint_velocity,
    }
}

fn apply_response<E>(
    state: &mut MixedArticulatedTrajectoryState,
    response: &MixedArticulatedImpulseResponse,
) -> Result<(), ArticulatedContactStepError<E>> {
    match (state, response) {
        (
            MixedArticulatedTrajectoryState::Free(state),
            MixedArticulatedImpulseResponse::Free(response),
        ) => {
            state.base = response.base_after;
            state
                .joint_velocity
                .clone_from(&response.joint_velocity_after);
        }
        (
            MixedArticulatedTrajectoryState::Prescribed(state),
            MixedArticulatedImpulseResponse::Prescribed(response),
        ) => {
            state.base = response.base_after;
            state
                .joint_velocity
                .clone_from(&response.joint_velocity_after);
        }
        _ => {
            return Err(ArticulatedContactStepError::Invalid(
                "native impulse changed boundary kind",
            ));
        }
    }
    Ok(())
}

fn add_response(
    response: &MixedArticulatedImpulseResponse,
    changes: &mut Impulses,
    sources: &mut Impulses,
    kinetic: &mut Dd,
    origin: Vec3,
) {
    let view = response.impulse_view();
    changes.add(
        view.linear_momentum_change_world,
        view.angular_momentum_change_world,
        origin,
    );
    *kinetic = *kinetic + Dd::from_f64(view.kinetic_energy_change);
    match response {
        MixedArticulatedImpulseResponse::Free(response) => sources.add(
            response.external_linear_impulse_world,
            response.external_angular_impulse_world,
            origin,
        ),
        MixedArticulatedImpulseResponse::Prescribed(response) => {
            sources.add(
                response.external_impulse_world.force,
                response.external_impulse_world.torque,
                origin,
            );
            sources.add(
                response.support_impulse_world.force,
                response.support_impulse_world.torque,
                origin,
            );
        }
    }
}

fn kinematics<E>(
    systems: &[MixedArticulatedTrajectorySystem<'_>],
    states: &[MixedArticulatedTrajectoryState],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<Kinematics>, ArticulatedContactStepError<E>> {
    systems
        .iter()
        .zip(states)
        .enumerate()
        .map(|(index, (system, state))| {
            poll(cancelled)?;
            let input = native_input(system.model, state);
            let view = input.system_view();
            forward_kinematics(
                system.model,
                view.base,
                view.joint_position,
                view.joint_velocity,
            )
            .map_err(|error| dynamics(index, error))
        })
        .collect()
}

pub(super) fn totals<E>(
    systems: &[MixedArticulatedTrajectorySystem<'_>],
    states: &[MixedArticulatedTrajectoryState],
    origin: Vec3,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ArticulatedTrajectoryMomentum, ArticulatedContactStepError<E>> {
    let mut spatial = Impulses::default();
    let mut kinetic = Dd::ZERO;
    for (index, (system, state)) in systems.iter().zip(states).enumerate() {
        poll(cancelled)?;
        let input = native_input(system.model, state);
        let view = input.system_view();
        let base = BaseState {
            world_from_base: Se3::from_parts(
                view.base.world_from_base.rotation(),
                view.base.world_from_base.translation() - origin,
            )
            .map_err(|error| dynamics(index, error.into()))?,
            ..view.base
        };
        let kinematics =
            forward_kinematics(system.model, base, view.joint_position, view.joint_velocity)
                .map_err(|error| dynamics(index, error))?;
        let momentum = momentum_from_kinematics(system.model, &kinematics)
            .map_err(|error| dynamics(index, error))?;
        spatial.add(
            momentum.linear_momentum_world,
            momentum.angular_momentum_world,
            Vec3::default(),
        );
        kinetic = kinetic + Dd::from_f64(momentum.kinetic_energy);
    }
    let total = ArticulatedTrajectoryMomentum {
        linear_world: spatial.linear,
        angular_about_origin_world: spatial.angular,
        kinetic_energy: kinetic.to_f64(),
    };
    if !finite_vector(total.linear_world)
        || !finite_vector(total.angular_about_origin_world)
        || !total.kinetic_energy.is_finite()
    {
        return Err(ArticulatedContactStepError::Invalid(
            "nonfinite mixed total momentum",
        ));
    }
    Ok(total)
}

fn prescribed_error<E>(
    system_index: usize,
    error: PrescribedBaseMotionError,
) -> ArticulatedContactStepError<E> {
    if matches!(error, PrescribedBaseMotionError::Cancelled) {
        ArticulatedContactStepError::Cancelled
    } else {
        ArticulatedContactStepError::PrescribedDrift {
            system_index,
            error,
        }
    }
}
