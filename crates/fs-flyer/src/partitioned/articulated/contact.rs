//! Sustained Coulomb contact through the original dry/fluid Lie midpoint.
//!
//! Normal and tangent reactions act on the same complete effective mass used
//! for fluid bias and implicit joint impedance. This is a bounded fixed-feature
//! sustained mode, not an impact reset or an automatic collision detector.

use super::{
    ArticulatedLoads, ArticulatedStage, ArticulatedState, ArticulatedStep, ArticulatedStepConfig,
    JointImpedance, StepSystem, VelocityResponseContext, admit_state, cancelled, pose_at,
    step_systems_with_velocity_response,
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

/// Simultaneous sustained contact among independently assembled fluid-bearing trees.
#[path = "contact/systems.rs"]
pub mod systems;

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
    midpoints: Vec<ArticulatedState>,
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
    native: &[&Kinematics],
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
                native,
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
    contexts: &[VelocityResponseContext<'_>],
    free_increments: &[Vec<f64>],
    contacts: &[ArticulatedFrictionalContact],
    before: &[ContactFrameGeometry],
    initial_native: &[Kinematics],
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
            let system = term.point.system_index;
            let arm = finite(term.point.point_body.norm())?;
            let input = initial_native[system].body_twist[term.point.link_index];
            scale = finite(scale + input.linear.norm() + input.angular.norm() * arm)?;
            let columns = contexts[system]
                .kinematics
                .body_jacobian(term.point.link_index);
            for (i, du) in free_increments[system].iter().enumerate() {
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

// A row touches at most two independently prepared trees. The same-tree case
// combines both endpoint wrenches BEFORE solving, preserving the original row.
struct RowPart {
    system: usize,
    jacobian: Vec<f64>,
    response: Vec<f64>,
}

// Every row uses the original point wrench and the SAME native stage map.
fn solve_midpoint(
    contexts: &[VelocityResponseContext<'_>],
    increments: &mut [Vec<f64>],
    initial_velocities: &[Vec<f64>],
    contacts: &[ArticulatedFrictionalContact],
    initial_gaps: &[f64],
    geometry: &[ContactFrameGeometry],
    config: FrictionalContactConfig,
    is_cancelled: &mut dyn FnMut() -> bool,
) -> Result<FrictionalContactResult, Refusal> {
    let free: Vec<Vec<_>> = initial_velocities
        .iter()
        .zip(increments.iter())
        .map(|(initial, increment)| {
            initial
                .iter()
                .zip(increment)
                .map(|(u, du)| u + 0.5 * du)
                .collect()
        })
        .collect();
    let mut rows: Vec<Vec<RowPart>> = Vec::with_capacity(3 * contacts.len());
    let mut laws = Vec::with_capacity(contacts.len());
    for (i, frame) in geometry.iter().enumerate() {
        poll(is_cancelled)?;
        let mut speeds = [0.0; 3];
        for (a, g) in frame.rows.iter().enumerate() {
            let mut parts: Vec<RowPart> = Vec::with_capacity(2);
            for term in core::iter::once(g.first).chain(g.second) {
                let system = term.point.system_index;
                let index = if let Some(index) = parts.iter().position(|p| p.system == system) {
                    index
                } else {
                    parts.push(RowPart {
                        system,
                        jacobian: vec![0.0; increments[system].len()],
                        response: Vec::new(),
                    });
                    parts.len() - 1
                };
                let columns = contexts[system]
                    .kinematics
                    .body_jacobian(term.point.link_index);
                for (j, value) in parts[index].jacobian.iter_mut().enumerate() {
                    *value = finite(*value + term.direction.pairing(columns[native_index(j)]))?;
                }
            }
            let mut speed = 0.0;
            for part in &mut parts {
                speed = finite(speed + dot(&part.jacobian, &free[part.system]))?;
                part.response = if part.jacobian.iter().any(|j| *j != 0.0) {
                    contexts[part.system].mass.solve(&part.jacobian)?
                } else {
                    vec![0.0; part.jacobian.len()]
                };
            }
            speeds[a] = finite(speed - g.plane_speed)?;
            rows.push(parts);
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
            let (mut ij, mut ji) = (0.0, 0.0);
            for first in &rows[i] {
                for second in &rows[j] {
                    if first.system == second.system {
                        ij = finite(ij + dot(&first.jacobian, &second.response))?;
                        ji = finite(ji + dot(&second.jacobian, &first.response))?;
                    }
                }
            }
            // One half of the physical velocity response enters the midpoint.
            let value = finite(0.25 * (ij + ji))?;
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
    let mut impulses: Vec<Option<Vec<f64>>> = vec![None; contexts.len()];
    for (parts, p) in rows.iter().zip(solved.impulses.iter().flatten()) {
        if *p == 0.0 {
            continue;
        }
        for part in parts {
            let impulse =
                impulses[part.system].get_or_insert_with(|| vec![0.0; part.jacobian.len()]);
            for (total, j) in impulse.iter_mut().zip(&part.jacobian) {
                *total = finite(*total + p * j)?;
            }
        }
    }
    // Preserve the original single summed-impulse solve per affected tree.
    // Unconnected trees and exactly cancelling zero blocks require no solve.
    for (system, impulse) in impulses.iter().enumerate() {
        if let Some(impulse) = impulse.as_ref().filter(|p| p.iter().any(|x| *x != 0.0)) {
            let response = contexts[system].mass.solve(impulse)?;
            for (du, change) in increments[system].iter_mut().zip(response) {
                *du = finite(*du + change)?;
            }
        }
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
    native: &[&Kinematics],
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
                let pose = native[term.point.system_index].world_from_link[term.point.link_index];
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
    let inputs = [systems::ContactSystem {
        model,
        impedance,
        state,
        memory_inputs,
        smooth: config.smooth,
    }];
    let result = systems::step_with_contacts(
        &inputs,
        t_s,
        dt_s,
        contacts,
        config.contact,
        |time, stages| Ok(vec![evaluate(time, stages[0].state, stages[0].kinematics)?]),
        |path| {
            admit_path(&ContactPath {
                model,
                initial: state,
                midpoint: &path.midpoints[0],
                endpoint: path.endpoints[0],
                contacts: path.contacts,
                start_time_s: path.start_time_s,
                duration_s: path.duration_s,
            })
        },
        &mut is_cancelled,
    )?;
    let step = result
        .steps
        .into_iter()
        .next()
        .ok_or_else(|| invalid("missing accepted tree"))?;
    Ok(ContactStep {
        step,
        reaction_solve: result.reaction_solve,
        contact_ids: result.contact_ids,
        midpoint_relative_velocities: result.midpoint_relative_velocities,
        endpoint_relative_velocities: result.endpoint_relative_velocities,
        initial_gaps: result.initial_gaps,
        endpoint_gaps: result.endpoint_gaps,
        initial_normal_speeds: result.initial_normal_speeds,
        normal_speed_roundoff_mps: result.normal_speed_roundoff_mps,
        midpoint_law_residuals: result.midpoint_law_residuals,
        midpoint_velocity_residual_mps: result.midpoint_velocity_residual_mps,
        reaction_work_j: result.reaction_work_j,
        normal_work_j: result.normal_work_j,
        friction_work_j: result.friction_work_j,
        support_work_j: result.support_work_j,
        external_linear_impulse_world: result.external_linear_impulse_world,
        external_angular_impulse_world: result.external_angular_impulse_world,
        path_admission: result.path_admission,
    })
}
