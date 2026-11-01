//! Sustained Coulomb contact through the original dry/fluid Lie midpoint.
//!
//! Normal and tangent reactions act on the same complete effective mass used
//! for fluid bias and implicit joint impedance. This is a bounded fixed-feature
//! sustained mode, not an impact reset or an automatic collision detector.

use super::{
    ArticulatedLoads, ArticulatedState, ArticulatedStep, ArticulatedStepConfig, JointImpedance,
    VelocityResponseContext, admit_state, cancelled, pose_at, step_with_velocity_response,
};
use crate::Refusal;
use crate::addedmass::joint::native_index;
use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactConfig, FrictionalContactDeterminacy,
    FrictionalContactError, FrictionalContactResiduals, FrictionalContactResult, NormalVelocityLaw,
    solve_frictional_contacts,
};
use fs_ga::{Twist, Vec3, Wrench};
use fs_mbd::articulated::contact_geometry::{
    ContactCounterpart, ContactFrameGeometry, prepare_contact_frame,
};
use fs_mbd::articulated::generalized::GeneralizedKinematics;
use fs_mbd::articulated::{ArticulatedModel, BaseState, Kinematics, forward_kinematics};
pub use fs_mbd::articulated_contacts::friction::ArticulatedFrictionalContact;

/// The original smooth and contact controls retain their distinct units.
#[derive(Clone, Copy, Debug)]
pub struct ContactStepConfig<'a> {
    pub smooth: ArticulatedStepConfig<'a>,
    pub contact: FrictionalContactConfig,
}

/// Candidate numerical configuration path for the REQUIRED geometry owner.
///
/// The callback must admit complete feature coverage, nonpenetration and the
/// continued applicability of these material witnesses and frozen world frames.
/// A finite set of samples is not that proof. Geometry/history stays private
/// until this callback and final cancellation admission both succeed.
pub struct ContactPath<'a> {
    pub model: &'a ArticulatedModel,
    pub initial: &'a ArticulatedState,
    pub midpoint: &'a ArticulatedState,
    pub endpoint: &'a ArticulatedState,
    /// Ascending stable identity order; plane positions are at the window start.
    pub contacts: &'a [ArticulatedFrictionalContact],
    pub start_time_s: f64,
    pub duration_s: f64,
}
impl ContactPath<'_> {
    /// Query the original constant-midpoint-twist / linear-joint interpolant.
    /// Returned twists differentiate this CONFIGURATION path, not a re-integrated
    /// force trajectory. This floating query is not certified continuous collision
    /// detection or an enclosure of the exact dynamical trajectory.
    pub fn sample(
        &self,
        fraction: f64,
        mut is_cancelled: impl FnMut() -> bool,
    ) -> Result<Kinematics, Refusal> {
        poll(&mut is_cancelled)?;
        if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
            return Err(invalid("path fraction must be finite and in [0,1]"));
        }
        let mut motion = self.midpoint.motion.clone();
        if fraction == 0.0 {
            motion.base.world_from_base = self.initial.motion.base.world_from_base;
            motion
                .joint_position
                .clone_from(&self.initial.motion.joint_position);
        } else if fraction == 1.0 {
            motion.base.world_from_base = self.endpoint.motion.base.world_from_base;
            motion
                .joint_position
                .clone_from(&self.endpoint.motion.joint_position);
        } else {
            let time = fraction * self.duration_s;
            motion.base.world_from_base = pose_at(
                self.initial.motion.base.world_from_base,
                self.midpoint.motion.base.twist_body,
                time,
            )?;
            for (i, q) in motion.joint_position.iter_mut().enumerate() {
                *q = self.initial.motion.joint_position[i]
                    + time * self.midpoint.motion.joint_velocity[i];
            }
        }
        let result = forward_kinematics(
            self.model,
            BaseState::prescribed(
                motion.base.world_from_base,
                motion.base.twist_body,
                Twist::zero(),
            ),
            &motion.joint_position,
            &motion.joint_velocity,
        )
        .map_err(invalid)?;
        poll(&mut is_cancelled)?;
        Ok(result)
    }
}

/// One accepted sustained-contact step, with actual midpoint reaction work.
#[derive(Clone, Debug)]
pub struct ContactStep<P> {
    pub step: ArticulatedStep,
    /// Original auxiliary half-step solve, aligned with the velocity increment
    /// actually accepted by the smooth owner. Its `dissipated_energy` describes
    /// that auxiliary projection; it is NOT physical friction work or heat.
    pub reaction_solve: FrictionalContactResult,
    pub contact_ids: Vec<u64>,
    /// Actual native material-point velocities, in normal/tangent/tangent order.
    pub midpoint_relative_velocities: Vec<[f64; 3]>,
    pub endpoint_relative_velocities: Vec<[f64; 3]>,
    /// Actual material gaps, in the same stable identity order.
    pub initial_gaps: Vec<f64>,
    pub endpoint_gaps: Vec<f64>,
    /// Original normal speeds are retained, including arithmetic-sized negative
    /// values. No input state or gap is projected to erase them.
    pub initial_normal_speeds: Vec<f64>,
    /// Per-contact arithmetic allowance at the accepted midpoint, m/s. It scales
    /// uncancelled original material/support motion and the ORIGINAL free velocity
    /// increment before contact. User solver tolerances do not enter this value.
    pub normal_speed_roundoff_mps: Vec<f64>,
    /// Residuals reconstructed from accepted midpoint velocities and the aligned
    /// reaction impulses, independently of the auxiliary solver's receipt.
    pub midpoint_law_residuals: FrictionalContactResiduals,
    /// Actual midpoint velocity minus its aligned auxiliary prediction, m/s.
    pub midpoint_velocity_residual_mps: f64,
    /// Reaction impulse paired with actual midpoint material velocities, J.
    pub reaction_work_j: f64,
    /// Normal reaction paired with actual relative midpoint normal speed, J.
    pub normal_work_j: f64,
    /// SIGNED tangential impulse work at the actual midpoint, J; dissipative
    /// within the declared residual budget. No heat/deformation split is inferred.
    pub friction_work_j: f64,
    /// Work delivered by prescribed translating supports, J.
    pub support_work_j: f64,
    /// External contact impulse from prescribed supports; internal pairs cancel.
    pub external_linear_impulse_world: Vec3,
    pub external_angular_impulse_world: Vec3,
    pub path_admission: P,
}

struct AcceptedReaction {
    solve: FrictionalContactResult,
    midpoint: ArticulatedState,
    geometry: Vec<ContactFrameGeometry>,
    residuals: FrictionalContactResiduals,
    velocity_residual: f64,
    normal_speed_roundoff: Vec<f64>,
}

fn refusal(code: &'static str, message: impl ToString) -> Refusal {
    Refusal {
        code,
        message: message.to_string(),
        ranked_repairs: vec![
            "resolve impacts and contact mode transitions, reduce the smooth step, or inspect the declared geometry and residual budgets".into(),
        ],
    }
}
fn invalid(message: impl ToString) -> Refusal {
    refusal("articulated-contact-invalid", message)
}
fn unresolved(message: impl ToString) -> Refusal {
    refusal("articulated-contact-unresolved", message)
}
fn finite(value: f64) -> Result<f64, Refusal> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("nonfinite contact arithmetic"))
    }
}
fn vector(value: Vec3) -> Result<Vec3, Refusal> {
    finite(value.x)?;
    finite(value.y)?;
    finite(value.z)?;
    Ok(value)
}
fn poll(is_cancelled: &mut dyn FnMut() -> bool) -> Result<(), Refusal> {
    if is_cancelled() {
        Err(cancelled())
    } else {
        Ok(())
    }
}
fn velocity(state: &ArticulatedState) -> Vec<f64> {
    let u = state.motion.base.twist_body.to_array();
    let mut result: Vec<_> = (0..6).map(|i| u[native_index(i)]).collect();
    result.extend_from_slice(&state.motion.joint_velocity);
    result
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn kinematics(model: &ArticulatedModel, state: &ArticulatedState) -> Result<Kinematics, Refusal> {
    forward_kinematics(
        model,
        BaseState::prescribed(
            state.motion.base.world_from_base,
            state.motion.base.twist_body,
            Twist::zero(),
        ),
        &state.motion.joint_position,
        &state.motion.joint_velocity,
    )
    .map_err(invalid)
}
fn frames(
    native: &Kinematics,
    contacts: &[ArticulatedFrictionalContact],
    elapsed: f64,
    is_cancelled: &mut dyn FnMut() -> bool,
) -> Result<Vec<ContactFrameGeometry>, Refusal> {
    contacts
        .iter()
        .enumerate()
        .map(|(index, c)| {
            poll(is_cancelled)?;
            let other = match c.other {
                ContactCounterpart::PrescribedPlane {
                    point_world,
                    velocity_world,
                } => ContactCounterpart::PrescribedPlane {
                    point_world: point_world + velocity_world.scale(elapsed),
                    velocity_world,
                },
                body => body,
            };
            prepare_contact_frame(
                c.point,
                other,
                c.normal_world,
                c.tangent_world,
                index,
                &[native],
            )
            .map_err(invalid)
        })
        .collect()
}
fn nonnegative_gaps(frames: &[ContactFrameGeometry]) -> Result<Vec<f64>, Refusal> {
    frames
        .iter()
        .map(|frame| {
            let gap = finite(frame.rows[0].gap)?;
            if gap < 0.0 {
                Err(unresolved(
                    "an actual material witness penetrates its counterpart",
                ))
            } else {
                Ok(gap)
            }
        })
        .collect()
}

// A persistent zero normal speed can retain floating error from rotating
// frames or from cancelling free load/contact responses. Bound those actual
// uncancelled motion terms, never a fixed speed floor or a user solve tolerance.
// In unforced pure approach this allowance scales the approach itself, so even
// an arbitrarily small resolved closing speed still requires an impact reset.
fn normal_speed_admission(
    context: &VelocityResponseContext<'_>,
    free_increment: &[f64],
    contacts: &[ArticulatedFrictionalContact],
    before: &[ContactFrameGeometry],
    initial_native: &Kinematics,
) -> Result<Vec<f64>, Refusal> {
    let mut allowances = Vec::with_capacity(contacts.len());
    for (c, frame) in contacts.iter().zip(before) {
        let mut scale = match c.other {
            ContactCounterpart::PrescribedPlane { velocity_world, .. } => {
                finite(velocity_world.norm())?
            }
            ContactCounterpart::Body(_) => 0.0,
        };
        for term in core::iter::once(frame.rows[0].first).chain(frame.rows[0].second) {
            let arm = finite(term.point.point_body.norm())?;
            let input = initial_native.body_twist[term.point.link_index];
            scale = finite(scale + input.linear.norm() + input.angular.norm() * arm)?;
            let columns = context.kinematics.body_jacobian(term.point.link_index);
            for (i, du) in free_increment.iter().enumerate() {
                let column = columns[native_index(i)];
                scale = finite(
                    scale + du.abs() * (column.linear.norm() + column.angular.norm() * arm),
                )?;
            }
        }
        let allowance = finite(256.0 * f64::EPSILON * scale)?;
        let row = frame.rows[0];
        if row.gap <= c.activation_tolerance && row.relative_normal_speed < -allowance {
            return Err(refusal(
                "articulated-contact-impact-required",
                format!(
                    "active contact {} approaches at {:e} m/s beyond arithmetic allowance {:e}; resolve its impact before sustained evolution",
                    c.id, row.relative_normal_speed, allowance,
                ),
            ));
        }
        allowances.push(allowance);
    }
    Ok(allowances)
}

// Every row uses the original point wrench and the SAME native stage map.
fn solve_midpoint(
    context: &VelocityResponseContext<'_>,
    increment: &mut [f64],
    initial_velocity: &[f64],
    contacts: &[ArticulatedFrictionalContact],
    initial_gaps: &[f64],
    geometry: &[ContactFrameGeometry],
    config: FrictionalContactConfig,
    is_cancelled: &mut dyn FnMut() -> bool,
) -> Result<FrictionalContactResult, Refusal> {
    let n = increment.len();
    let free: Vec<_> = initial_velocity
        .iter()
        .zip(increment.iter())
        .map(|(u, du)| u + 0.5 * du)
        .collect();
    let mut rows = Vec::with_capacity(3 * contacts.len());
    let mut responses = Vec::with_capacity(3 * contacts.len());
    let mut laws = Vec::with_capacity(contacts.len());
    for (i, frame) in geometry.iter().enumerate() {
        poll(is_cancelled)?;
        let mut speeds = [0.0; 3];
        for (a, g) in frame.rows.iter().enumerate() {
            let mut row = vec![0.0; n];
            for term in core::iter::once(g.first).chain(g.second) {
                let columns = context.kinematics.body_jacobian(term.point.link_index);
                for (j, value) in row.iter_mut().enumerate() {
                    *value = finite(*value + term.direction.pairing(columns[native_index(j)]))?;
                }
            }
            speeds[a] = finite(dot(&row, &free) - g.plane_speed)?;
            responses.push(context.mass.solve(&row)?);
            rows.push(row);
        }
        laws.push(FrictionalContact {
            id: contacts[i].id,
            // Activation belongs to the admitted incoming configuration. Signed
            // provisional midpoint geometry is NOT published or gap-repaired.
            gap: initial_gaps[i],
            activation_tolerance: contacts[i].activation_tolerance,
            relative_velocity: speeds,
            friction_coefficient: contacts[i].friction_coefficient,
            normal_law: NormalVelocityLaw::Target { speed: 0.0 },
        });
    }
    let count = rows.len();
    let mut matrix = vec![0.0; count * count];
    for i in 0..count {
        poll(is_cancelled)?;
        for j in 0..=i {
            // One half of the physical velocity response enters the midpoint.
            let value =
                finite(0.25 * (dot(&rows[i], &responses[j]) + dot(&rows[j], &responses[i])))?;
            matrix[i * count + j] = value;
            matrix[j * count + i] = value;
        }
    }
    let solved =
        solve_frictional_contacts(&laws, &matrix, config, &mut *is_cancelled).map_err(|e| {
            if matches!(e, FrictionalContactError::Cancelled) {
                cancelled()
            } else {
                unresolved(e)
            }
        })?;
    let mut impulse = vec![0.0; n];
    for (row, p) in rows.iter().zip(solved.impulses.iter().flatten()) {
        for (total, j) in impulse.iter_mut().zip(row) {
            *total = finite(*total + p * j)?;
        }
    }
    let response = context.mass.solve(&impulse)?;
    for (du, change) in increment.iter_mut().zip(response) {
        *du = finite(*du + change)?;
    }
    Ok(solved)
}

fn actual_law(
    solved: &FrictionalContactResult,
    frames: &[ContactFrameGeometry],
    contacts: &[ArticulatedFrictionalContact],
    config: FrictionalContactConfig,
) -> Result<(FrictionalContactResiduals, f64, bool), Refusal> {
    let mut r = FrictionalContactResiduals::default();
    let mut velocity_residual = 0.0_f64;
    for (i, frame) in frames.iter().enumerate() {
        for a in 0..3 {
            velocity_residual = velocity_residual.max(
                finite(
                    frame.rows[a].relative_normal_speed - solved.relative_velocities_after[i][a],
                )?
                .abs(),
            );
        }
    }
    for &i in &solved.active_indices {
        let p = solved.impulses[i];
        let v = frames[i].rows.map(|r| r.relative_normal_speed);
        let radius = finite(contacts[i].friction_coefficient * p[0])?;
        let speed = finite(v[1].hypot(v[2]))?;
        r.normal_impulse_violation = r.normal_impulse_violation.max((-p[0]).max(0.0));
        r.normal_velocity_violation = r.normal_velocity_violation.max((-v[0]).max(0.0));
        r.normal_complementarity = r.normal_complementarity.max(finite(p[0] * v[0])?.abs());
        r.cone_violation = r
            .cone_violation
            .max(finite(p[1].hypot(p[2]) - radius)?.max(0.0));
        r.tangential_duality_gap = r
            .tangential_duality_gap
            .max(finite(radius * speed + p[1] * v[1] + p[2] * v[2])?.abs());
        if speed > config.velocity_tolerance {
            r.sliding_impulse_residual = r.sliding_impulse_residual.max(finite(
                (p[1] + radius * v[1] / speed).hypot(p[2] + radius * v[2] / speed),
            )?);
        }
    }
    let ready = velocity_residual <= config.velocity_tolerance
        && r.normal_impulse_violation <= config.impulse_tolerance
        && r.normal_velocity_violation <= config.velocity_tolerance
        && r.normal_complementarity <= config.complementarity_tolerance
        && r.cone_violation <= config.impulse_tolerance
        && r.sliding_impulse_residual <= config.impulse_tolerance
        && r.tangential_duality_gap <= config.complementarity_tolerance;
    Ok((r, velocity_residual, ready))
}

// Contact impulses must remain physically internal independently of the smooth
// integrator's finite-step momentum defect. Never hide a finite-gap couple.
fn spatial_impulses(
    frames: &[ContactFrameGeometry],
    native: &Kinematics,
    impulses: &[[f64; 3]],
) -> Result<(Vec3, Vec3), Refusal> {
    let (mut total_p, mut total_l, mut external_p, mut external_l) = (
        Vec3::default(),
        Vec3::default(),
        Vec3::default(),
        Vec3::default(),
    );
    let (mut p_scale, mut l_scale) = (0.0, 0.0);
    for (frame, impulse) in frames.iter().zip(impulses) {
        for (row, &p) in frame.rows.iter().zip(impulse) {
            for term in core::iter::once(row.first).chain(row.second) {
                let pose = native.world_from_link[term.point.link_index];
                let local = Wrench::new(
                    term.direction.torque.scale(p),
                    term.direction.force.scale(p),
                );
                let force = vector(pose.rotation().rotate(local.force).map_err(invalid)?)?;
                let spin = vector(pose.rotation().rotate(local.torque).map_err(invalid)?)?;
                let orbital = vector(pose.translation().cross(force))?;
                p_scale = finite(p_scale + force.norm())?;
                l_scale = finite(l_scale + spin.norm() + orbital.norm())?;
                total_p = vector(total_p + force)?;
                total_l = vector(total_l + spin + orbital)?;
            }
            if row.second.is_none() {
                let force = vector(row.axis_world.scale(p))?;
                external_p = vector(external_p + force)?;
                external_l = vector(external_l + row.world_point.cross(force))?;
            }
        }
    }
    let p_residual = finite((total_p - external_p).norm())?;
    let l_residual = finite((total_l - external_l).norm())?;
    if p_residual > 512.0 * f64::EPSILON * p_scale || l_residual > 512.0 * f64::EPSILON * l_scale {
        return Err(unresolved(
            "internal point reactions create a nonphysical spatial impulse or finite-gap couple",
        ));
    }
    Ok((external_p, external_l))
}

/// Advance one sustained, fixed-feature contact step with full dry/fluid inertia.
///
/// All witnesses belong to this one tree (`system_index == 0`) or prescribed
/// translating planes. Only zero normal targets/restitution are admitted. A
/// resolved approaching active normal requires the separate physical impact
/// reset. Arithmetic-sized negative input speeds remain unprojected and are
/// reported with allowances from uncancelled point motion and free load response;
/// no user solver tolerance or absolute speed floor authorizes this admission.
/// An unresolved friction stopping/reversal or release/recontact refuses;
/// shrinking the step never silently replaces an impact with midpoint support.
///
/// The original midpoint evaluates fluid bias, articulated geometry and loads,
/// and eliminates implicit joint impedance. One simultaneous Coulomb solve uses
/// half of THAT effective mobility. The reaction receipt is aligned with the
/// previous response which generated the accepted increment, then independently
/// checked against actual accepted midpoint contact speeds. It is not the next
/// fixed-point iterate's unapplied impulse.
///
/// Every initial, actual midpoint and endpoint gap is nonnegative. The required
/// geometry callback admits the complete configuration path and omitted-feature
/// coverage; endpoint checks alone are not collision detection. Reactions at a
/// positive gap are resolved only within the declared activation band. Nonlinear
/// smooth energy/momentum defects retain the original refinement-based scope;
/// exact finite-step global conservation and physical heat allocation are absent.
#[allow(clippy::too_many_arguments)]
pub fn step_with_contacts<F, P>(
    model: &ArticulatedModel,
    impedance: &[JointImpedance],
    state: &ArticulatedState,
    memory_inputs: &[f64],
    t_s: f64,
    dt_s: f64,
    contacts: &[ArticulatedFrictionalContact],
    config: ContactStepConfig<'_>,
    evaluate: F,
    admit_path: impl FnOnce(&ContactPath<'_>) -> Result<P, Refusal>,
    mut is_cancelled: impl FnMut() -> bool,
) -> Result<ContactStep<P>, Refusal>
where
    F: Fn(f64, &ArticulatedState, &GeneralizedKinematics) -> Result<ArticulatedLoads, Refusal>,
{
    poll(&mut is_cancelled)?;
    if model.link_count() > super::MAX_ARTICULATED_LINKS
        || model.dof_count() > super::MAX_CONTROL_COORDS
        || contacts.len() > 21
        || contacts.len() > config.contact.maximum_contacts
        || !t_s.is_finite()
        || !dt_s.is_finite()
        || dt_s <= 0.0
    {
        return Err(invalid("invalid bounded model, contacts or step duration"));
    }
    let mut ordered = contacts.to_vec();
    ordered.sort_unstable_by_key(|c| c.id);
    if ordered.windows(2).any(|pair| pair[0].id == pair[1].id) {
        return Err(invalid("duplicate contact identity"));
    }
    for c in &ordered {
        if c.point.system_index != 0
            || matches!(c.other,ContactCounterpart::Body(p) if p.system_index!=0)
            || !c.activation_tolerance.is_finite()
            || c.activation_tolerance < 0.0
            || !c.friction_coefficient.is_finite()
            || c.friction_coefficient < 0.0
            || !matches!(
                c.normal_law,
                NormalVelocityLaw::Target { speed: 0.0 }
                    | NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 }
            )
        {
            return Err(invalid(
                "sustained contacts require this tree, finite friction and zero normal target/restitution",
            ));
        }
    }
    admit_state(model, &state.motion)?;
    let initial_native = kinematics(model, state)?;
    let before = frames(&initial_native, &ordered, 0.0, &mut is_cancelled)?;
    let initial_gaps = nonnegative_gaps(&before)?;
    let initial_velocity = velocity(state);
    let mut previous = None;
    let mut accepted = None;
    let step = step_with_velocity_response(
        model,
        impedance,
        state,
        memory_inputs,
        t_s,
        dt_s,
        config.smooth,
        evaluate,
        |context, increment, cancel| {
            let normal_speed_roundoff =
                normal_speed_admission(&context, increment, &ordered, &before, &initial_native)?;
            let geometry = frames(
                context.kinematics.native_kinematics(),
                &ordered,
                0.5 * dt_s,
                cancel,
            )?;
            let next = solve_midpoint(
                &context,
                increment,
                &initial_velocity,
                &ordered,
                &initial_gaps,
                &geometry,
                config.contact,
                cancel,
            )?;
            // The smooth owner publishes its evaluated candidate, whose delta
            // came from the PREVIOUS response. Keep the reaction aligned with it.
            let aligned = previous.replace(next);
            if let Some(solve) = aligned {
                let (residuals, velocity_residual, ready) =
                    actual_law(&solve, &geometry, &ordered, config.contact)?;
                accepted = ready.then(|| AcceptedReaction {
                    solve,
                    midpoint: context.midpoint.clone(),
                    geometry,
                    residuals,
                    velocity_residual,
                    normal_speed_roundoff,
                });
                Ok(ready)
            } else {
                Ok(false)
            }
        },
        &mut is_cancelled,
    )?;
    let accepted = accepted.ok_or_else(|| unresolved("no aligned accepted reaction"))?;
    nonnegative_gaps(&accepted.geometry)?;
    let after = frames(
        &kinematics(model, &step.state)?,
        &ordered,
        dt_s,
        &mut is_cancelled,
    )?;
    let endpoint_gaps = nonnegative_gaps(&after)?;
    for &i in &accepted.solve.active_indices {
        let p = accepted.solve.impulses[i];
        let first = before[i].rows.map(|r| r.relative_normal_speed);
        let last = after[i].rows.map(|r| r.relative_normal_speed);
        let c = ordered[i];
        if p[0] > config.contact.impulse_tolerance
            && (accepted.geometry[i].rows[0].gap > c.activation_tolerance
                || endpoint_gaps[i] > c.activation_tolerance
                || first[0] > config.contact.velocity_tolerance
                || last[0] < -config.contact.velocity_tolerance)
        {
            return Err(unresolved(
                "resolve separation or recontact before applying a sustained normal reaction",
            ));
        }
        let first_slip = finite(first[1].hypot(first[2]))?;
        if p[1].hypot(p[2]) > config.contact.impulse_tolerance
            && first_slip > config.contact.velocity_tolerance
            && finite(first[1] * last[1] + first[2] * last[2])?
                < -config.contact.velocity_tolerance * first_slip
        {
            return Err(unresolved(
                "the step crosses an unresolved friction stopping or turning event",
            ));
        }
    }
    if matches!(
        accepted.solve.determinacy,
        FrictionalContactDeterminacy::UniqueVelocity(_)
    ) {
        let mut common = None;
        for &i in &accepted.solve.active_indices {
            if let ContactCounterpart::PrescribedPlane { velocity_world, .. } = ordered[i].other {
                if common.is_some_and(|v| v != velocity_world) {
                    return Err(unresolved(
                        "redundant reactions do not establish work for independently moving supports",
                    ));
                }
                common = Some(velocity_world);
            }
        }
    }
    let midpoint_native = kinematics(model, &accepted.midpoint)?;
    let (external_linear_impulse_world, external_angular_impulse_world) = spatial_impulses(
        &accepted.geometry,
        &midpoint_native,
        &accepted.solve.impulses,
    )?;
    let (mut reaction_work, mut normal_work, mut friction_work, mut support_work) =
        (0.0, 0.0, 0.0, 0.0);
    for (g, p) in accepted.geometry.iter().zip(&accepted.solve.impulses) {
        normal_work = finite(normal_work + p[0] * g.rows[0].relative_normal_speed)?;
        for a in 0..3 {
            let relative = finite(p[a] * g.rows[a].relative_normal_speed)?;
            let support = finite(p[a] * g.rows[a].plane_speed)?;
            reaction_work = finite(reaction_work + relative + support)?;
            support_work = finite(support_work + support)?;
            if a > 0 {
                friction_work = finite(friction_work + relative)?;
            }
        }
    }
    if friction_work > config.contact.complementarity_tolerance
        || normal_work.abs()
            > config.contact.complementarity_tolerance * (ordered.len().max(1) as f64)
    {
        return Err(unresolved(
            "accepted contact work exceeds the declared physical residual budget",
        ));
    }
    poll(&mut is_cancelled)?;
    let path_admission = admit_path(&ContactPath {
        model,
        initial: state,
        midpoint: &accepted.midpoint,
        endpoint: &step.state,
        contacts: &ordered,
        start_time_s: t_s,
        duration_s: dt_s,
    })?;
    poll(&mut is_cancelled)?;
    Ok(ContactStep {
        step,
        contact_ids: ordered.iter().map(|c| c.id).collect(),
        midpoint_relative_velocities: accepted
            .geometry
            .iter()
            .map(|g| g.rows.map(|r| r.relative_normal_speed))
            .collect(),
        endpoint_relative_velocities: after
            .iter()
            .map(|g| g.rows.map(|r| r.relative_normal_speed))
            .collect(),
        initial_gaps,
        endpoint_gaps,
        initial_normal_speeds: before
            .iter()
            .map(|g| g.rows[0].relative_normal_speed)
            .collect(),
        normal_speed_roundoff_mps: accepted.normal_speed_roundoff,
        midpoint_law_residuals: accepted.residuals,
        midpoint_velocity_residual_mps: accepted.velocity_residual,
        reaction_work_j: reaction_work,
        normal_work_j: normal_work,
        friction_work_j: friction_work,
        support_work_j: support_work,
        external_linear_impulse_world,
        external_angular_impulse_world,
        reaction_solve: accepted.solve,
        path_admission,
    })
}
