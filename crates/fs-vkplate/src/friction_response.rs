//! Joint normal, Coulomb and retained shear contact on the same body inverse.

use super::{Body, Contact, ContactQueue, Friction, contact_response};
use crate::contact_island::{
    Controls as NormalControls, SolveError as ContactError, SolveReport as NormalReport,
};
use crate::friction_island::{
    Controls, SolveError, SolveReport, Tangential, Workspace as IslandWorkspace, check_friction,
};
use crate::mindlin_contact::{Candidate as MindlinCandidate, MindlinContact};

#[derive(Clone, Default)]
pub(crate) struct Workspace {
    core: IslandWorkspace,
    active: Vec<usize>,
    tangents: Vec<Tangential>,
    normal_slot: Vec<usize>,
    predicted: Vec<f64>,
    previous_force: Vec<f64>,
    bound: Vec<f64>,
    paths: Vec<f64>,
    nonlocal: Vec<usize>,
    nodes: Vec<(usize, usize)>,
    seen: Vec<bool>,
    stuck: Vec<bool>,
    // Signed applied reactions, indexed like Instrument::frictions. These
    // accepted values also warm-start the next solve; trials never change them.
    applied: Vec<[f64; 2]>,
    mindlin_candidates: Vec<MindlinCandidate>,
}

impl Workspace {
    pub(crate) fn tangential_forces(&self) -> &[[f64; 2]] {
        &self.applied
    }
}

fn touch(work: &mut Workspace, bodies: &[Body], queue: &mut ContactQueue, b: usize, i: usize) {
    if bodies[b].nonlocal_response() {
        if !work.seen[b] {
            work.nonlocal.push(b);
            work.seen[b] = true;
        }
    } else if queue.mark[b][i] != queue.batch {
        queue.mark[b][i] = queue.batch;
        work.nodes.push((b, i));
    }
}

/// q contains positive compressive normal forces and opposing tangential
/// forces. The physical correction is -A^-1 B q for normal and all tangential
/// footprint columns, including the retained scalar shear ports.
fn apply_response(
    bodies: &mut [Body],
    contacts: &[Contact],
    frictions: &[Friction],
    mindlin: &[MindlinContact],
    normal_indices: &[usize],
    active: &[usize],
    nonlocal: &[usize],
    nodes: &[(usize, usize)],
    k: f64,
    force: &[f64],
    response: &mut [f64],
) {
    for &b in nonlocal {
        bodies[b].tmp.fill(0.0);
    }
    for &(b, i) in nodes {
        bodies[b].tmp[i] = 0.0;
    }
    for (slot, &ci) in normal_indices.iter().enumerate() {
        let c = &contacts[ci];
        for (terms, sign) in [(&c.a, 1.0), (&c.b, -1.0)] {
            for &(b, i, w) in terms {
                if w != 0.0 {
                    bodies[b].tmp[i] += sign * w * force[slot];
                }
            }
        }
    }
    let n = normal_indices.len();
    for (slot, &fi) in active.iter().enumerate() {
        let fr = &frictions[fi];
        for (axis, terms) in [&fr.sx, &fr.sy].into_iter().enumerate() {
            for &(b, i, w) in terms {
                if w != 0.0 {
                    bodies[b].tmp[i] += w * force[n + 2 * slot + axis];
                }
            }
        }
    }
    let retained_start = n + 2 * active.len();
    for (slot, contact) in mindlin.iter().enumerate() {
        for &(b, i, w) in contact.footprint() {
            if w != 0.0 {
                bodies[b].tmp[i] += w * force[retained_start + slot];
            }
        }
    }
    for &b in nonlocal {
        bodies[b].solve_force_response(k);
    }
    for &(b, i) in nodes {
        bodies[b].tmp[i] *= bodies[b].compliance(i, k);
    }
    for (slot, &ci) in normal_indices.iter().enumerate() {
        let c = &contacts[ci];
        response[slot] = 0.0;
        for (terms, sign) in [(&c.a, 1.0), (&c.b, -1.0)] {
            for &(b, i, w) in terms {
                if w != 0.0 {
                    response[slot] += sign * w * bodies[b].tmp[i];
                }
            }
        }
    }
    for (slot, &fi) in active.iter().enumerate() {
        let fr = &frictions[fi];
        for (axis, terms) in [&fr.sx, &fr.sy].into_iter().enumerate() {
            response[n + 2 * slot + axis] = terms
                .iter()
                .filter(|t| t.2 != 0.0)
                .map(|&(b, i, w)| w * bodies[b].tmp[i])
                .sum();
        }
    }
    for (slot, contact) in mindlin.iter().enumerate() {
        response[retained_start + slot] = contact
            .footprint()
            .iter()
            .filter(|t| t.2 != 0.0)
            .map(|&(b, i, w)| w * bodies[b].tmp[i])
            .sum();
    }
}

/// Validate every normal and tangential candidate at one common rounded body
/// pose, then publish the whole stage. Earlier body prediction is caller-owned.
/// Normal-only samples retain the existing independent/component fast paths.
pub(crate) fn solve_contacts_with_friction_checked(
    bodies: &mut [Body],
    contacts: &mut [Contact],
    frictions: &mut [Friction],
    mindlin: &mut [MindlinContact],
    queue: &mut ContactQueue,
    normal: &mut contact_response::Workspace,
    work: &mut Workspace,
    k: f64,
    normal_controls: NormalControls,
    controls: Controls,
) -> Result<(NormalReport, SolveReport), ContactError> {
    let fail = ContactError::Friction;
    let mut has_friction = false;
    for fr in frictions.iter() {
        let c = contacts.get(fr.contact).ok_or_else(|| {
            fail(SolveError::InvalidInput(
                "friction contact index is out of range",
            ))
        })?;
        if !fr.mu.is_finite() || fr.mu < 0.0 {
            return Err(fail(SolveError::InvalidInput(
                "friction coefficient must be finite and nonnegative",
            )));
        }
        has_friction |= c.active && fr.mu > 0.0;
    }
    if !has_friction && mindlin.is_empty() {
        let report = contact_response::solve_contacts_checked(
            bodies,
            contacts,
            queue,
            normal,
            k,
            normal_controls,
        )?;
        for fr in frictions.iter_mut() {
            fr.force = 0.0;
            fr.stuck = false;
        }
        work.applied.resize(frictions.len(), [0.0; 2]);
        work.applied.fill([0.0; 2]);
        return Ok((report, SolveReport::default()));
    }
    contact_response::prepare(bodies, contacts, queue, normal, k)?;
    let n = normal.indices.len();
    work.active.clear();
    work.tangents.clear();
    work.normal_slot.resize(contacts.len(), usize::MAX);
    work.normal_slot.fill(usize::MAX);
    for (slot, &ci) in normal.indices.iter().enumerate() {
        work.normal_slot[ci] = slot;
    }
    work.predicted.clear();
    work.predicted.extend_from_slice(&normal.predicted);
    work.previous_force.clear();
    work.previous_force
        .extend_from_slice(&normal.previous_force);
    work.bound.clear();
    work.bound.extend_from_slice(&normal.diagonal_bound);
    work.nonlocal.clear();
    work.nodes.clear();
    work.seen.resize(bodies.len(), false);
    work.seen.fill(false);
    queue.new_batch();
    for &ci in &normal.indices {
        for &(b, i, w) in contacts[ci].a.iter().chain(&contacts[ci].b) {
            if w != 0.0 {
                touch(work, bodies, queue, b, i);
            }
        }
    }
    for (fi, fr) in frictions.iter().enumerate() {
        let slot = work.normal_slot[fr.contact];
        if slot == usize::MAX || fr.mu == 0.0 {
            continue;
        }
        work.active.push(fi);
        work.tangents.push(Tangential {
            normal: slot,
            mu: fr.mu,
        });
        let previous = work.applied.get(fi).copied().unwrap_or([0.0; 2]);
        for (axis, terms) in [&fr.sx, &fr.sy].into_iter().enumerate() {
            let (mut slip, mut norm_bound) = (0.0, 0.0);
            for &(b, i, w) in terms {
                let body = bodies.get(b).ok_or_else(|| {
                    fail(SolveError::InvalidInput(
                        "friction body index is out of range",
                    ))
                })?;
                if i >= body.n || !w.is_finite() {
                    return Err(fail(SolveError::InvalidInput(
                        "friction node or footprint weight is invalid",
                    )));
                }
                let diagonal = body.compliance(i, k);
                if !diagonal.is_finite() || diagonal <= 0.0 {
                    return Err(fail(SolveError::InvalidResponse));
                }
                if w != 0.0 {
                    slip += w * (body.wn[i] - body.wp[i]);
                    norm_bound += w.abs() * diagonal.sqrt();
                    touch(work, bodies, queue, b, i);
                }
            }
            let bound = norm_bound * norm_bound;
            if !slip.is_finite() || !bound.is_finite() {
                return Err(fail(SolveError::InvalidResponse));
            }
            work.predicted.push(slip);
            work.previous_force.push(-previous[axis]);
            work.bound.push(bound);
        }
    }
    for contact in mindlin.iter() {
        let (mut slip, mut norm_bound) = (0.0, 0.0);
        for &(b, i, w) in contact.footprint() {
            let body = bodies.get(b).ok_or_else(|| {
                fail(SolveError::InvalidInput(
                    "Mindlin body index is out of range",
                ))
            })?;
            if i >= body.n || !w.is_finite() {
                return Err(fail(SolveError::InvalidInput(
                    "Mindlin node or footprint weight is invalid",
                )));
            }
            let diagonal = body.compliance(i, k);
            if !diagonal.is_finite() || diagonal <= 0.0 {
                return Err(fail(SolveError::InvalidResponse));
            }
            if w != 0.0 {
                slip += w * (body.wn[i] - body.wp[i]);
                norm_bound += w.abs() * diagonal.sqrt();
                touch(work, bodies, queue, b, i);
            }
        }
        let bound = norm_bound * norm_bound;
        if !slip.is_finite() || !bound.is_finite() {
            return Err(fail(SolveError::InvalidResponse));
        }
        work.predicted.push(slip);
        work.previous_force.push(contact.report().traction_n);
        work.bound.push(bound);
    }
    work.paths.resize(work.predicted.len(), 0.0);
    work.stuck.resize(work.active.len(), false);
    work.mindlin_candidates.clear();
    let Workspace {
        core,
        active,
        tangents,
        predicted,
        previous_force,
        bound,
        paths,
        nonlocal,
        nodes,
        stuck,
        applied,
        mindlin_candidates,
        ..
    } = work;
    let mut report = core
        .solve(
            predicted,
            previous_force,
            bound,
            n,
            tangents,
            mindlin.len(),
            controls,
            |force, response| {
                apply_response(
                    bodies,
                    contacts,
                    frictions,
                    mindlin,
                    &normal.indices,
                    active,
                    nonlocal,
                    nodes,
                    k,
                    force,
                    response,
                )
            },
            |slot, d| {
                contact_response::law(
                    &contacts[normal.indices[slot]],
                    normal.previous[slot],
                    normal.crush[slot].as_ref(),
                    d,
                    k,
                )
            },
            |slot, path| {
                mindlin[slot]
                    .trial(path, k)
                    .map_err(|error| SolveError::Mindlin {
                        contact: slot,
                        error,
                    })
            },
        )
        .map_err(|error| {
            fail(match error {
                SolveError::InvalidLaw { contact } => SolveError::InvalidLaw {
                    contact: normal.indices[contact],
                },
                SolveError::Coulomb { friction, error } => SolveError::Coulomb {
                    friction: active[friction],
                    error,
                },
                other => other,
            })
        })?;
    let forces = core.forces();
    // The last Krylov action is not the accepted force. Apply the entire
    // accepted aggregate once, including zero forces, before publication.
    apply_response(
        bodies,
        contacts,
        frictions,
        mindlin,
        &normal.indices,
        active,
        nonlocal,
        nodes,
        k,
        forces,
        paths,
    );
    report.compliance_applications += 1;
    for &b in nonlocal.iter() {
        if bodies[b]
            .wn
            .iter()
            .zip(&bodies[b].tmp)
            .any(|(&x, &dx)| !(x - dx).is_finite())
        {
            return Err(fail(SolveError::InvalidResponse));
        }
    }
    for &(b, i) in nodes.iter() {
        if !(bodies[b].wn[i] - bodies[b].tmp[i]).is_finite() {
            return Err(fail(SolveError::InvalidResponse));
        }
    }
    let mut displacement_residual = 0.0_f64;
    for (slot, &ci) in normal.indices.iter().enumerate() {
        let c = &contacts[ci];
        let mut d = -c.gap;
        for (terms, sign) in [(&c.a, 1.0), (&c.b, -1.0)] {
            for &(b, i, w) in terms {
                if w != 0.0 {
                    d += sign * w * (bodies[b].wn[i] - bodies[b].tmp[i]);
                }
            }
        }
        paths[slot] = d;
        let (expected, slope) =
            contact_response::law(c, normal.previous[slot], normal.crush[slot].as_ref(), d, k);
        let force = forces[slot];
        if !d.is_finite()
            || !force.is_finite()
            || force < 0.0
            || !expected.is_finite()
            || expected < 0.0
            || !slope.is_finite()
            || slope < 0.0
        {
            return Err(fail(SolveError::InvalidResponse));
        }
        let scale = controls.force_tolerance + controls.relative_tolerance * force.max(expected);
        if !scale.is_finite() {
            return Err(fail(SolveError::InvalidResponse));
        }
        report.normal_residual = report.normal_residual.max((force - expected).abs() / scale);
    }
    report.friction_work_j = 0.0;
    report.dissipated_j = 0.0;
    report.stored_change_j = 0.0;
    for (slot, &fi) in active.iter().enumerate() {
        let fr = &frictions[fi];
        let at = n + 2 * slot;
        let q = [forces[at], forces[at + 1]];
        let mut s = [0.0; 2];
        for (axis, terms) in [&fr.sx, &fr.sy].into_iter().enumerate() {
            for &(b, i, w) in terms {
                if w != 0.0 {
                    s[axis] += w * ((bodies[b].wn[i] - bodies[b].tmp[i]) - bodies[b].wp[i]);
                }
            }
            paths[at + axis] = s[axis];
        }
        let normal_force = forces[tangents[slot].normal];
        let slip = s[0].hypot(s[1]);
        let slip_scale = controls.displacement_tolerance
            + controls.relative_tolerance * slip.max(predicted[at].hypot(predicted[at + 1]));
        let check = check_friction(
            normal_force,
            q,
            s,
            [predicted[at], predicted[at + 1]],
            fr.mu,
            [bound[at], bound[at + 1]],
            controls,
        )
        .map_err(|error| {
            fail(SolveError::Coulomb {
                friction: fi,
                error,
            })
        })?;
        report.friction_residual = report.friction_residual.max(check.residual);
        // Signed mechanical work stays observable even at the declared
        // roundoff tolerance; it is never silently relabelled as heat.
        report.friction_work_j += check.friction_work_j;
        report.dissipated_j += (-check.friction_work_j).max(0.0);
        stuck[slot] = normal_force > 0.0 && slip <= slip_scale;
    }
    let retained_start = n + 2 * active.len();
    for (slot, contact) in mindlin.iter().enumerate() {
        let at = retained_start + slot;
        let path: f64 = contact
            .footprint()
            .iter()
            .filter(|t| t.2 != 0.0)
            .map(|&(b, i, w)| w * ((bodies[b].wn[i] - bodies[b].tmp[i]) - bodies[b].wp[i]))
            .sum();
        paths[at] = path;
        let candidate = contact
            .candidate(0.5 * path, forces[at], k, controls)
            .map_err(|error| {
                fail(SolveError::Mindlin {
                    contact: slot,
                    error,
                })
            })?;
        report.retained_residual = report.retained_residual.max(candidate.residual);
        report.friction_work_j += candidate.report.mechanical_work_j;
        report.stored_change_j += candidate.report.stored_change_j;
        report.dissipated_j += candidate.report.dissipated_j;
        mindlin_candidates.push(candidate);
    }
    for i in 0..paths.len() {
        let scale = controls.displacement_tolerance
            + controls.relative_tolerance
                * paths[i]
                    .abs()
                    .max(predicted[i].abs())
                    .max(core.paths()[i].abs());
        if !paths[i].is_finite() || !scale.is_finite() {
            return Err(fail(SolveError::InvalidResponse));
        }
        displacement_residual =
            displacement_residual.max((paths[i] - core.paths()[i]).abs() / scale);
    }
    if !report.normal_residual.is_finite()
        || !report.friction_residual.is_finite()
        || !report.retained_residual.is_finite()
        || !report.friction_work_j.is_finite()
        || !report.dissipated_j.is_finite()
        || !report.stored_change_j.is_finite()
        || report.normal_residual > 1.0
        || report.friction_residual > 1.0
        || report.retained_residual > 1.0
        || displacement_residual > 1.0
    {
        return Err(fail(SolveError::NoConvergence {
            newton_iterations: report.newton_iterations,
            normal_residual: report.normal_residual.max(displacement_residual),
            friction_residual: report.friction_residual,
            retained_residual: report.retained_residual,
            diagnosis: None,
        }));
    }
    // Every check precedes every physical write, including diagnostics and
    // the accepted signed reaction used by the next warm start.
    for &b in nonlocal.iter() {
        let body = &mut bodies[b];
        for (next, &delta) in body.wn.iter_mut().zip(&body.tmp) {
            *next -= delta;
        }
    }
    for &(b, i) in nodes.iter() {
        bodies[b].wn[i] -= bodies[b].tmp[i];
    }
    for (slot, &ci) in normal.indices.iter().enumerate() {
        let c = &mut contacts[ci];
        c.force = forces[slot];
        c.peak_force = c.peak_force.max(c.force);
        c.d_now = paths[slot];
        c.crush = normal.crush[slot];
    }
    for fr in frictions.iter_mut() {
        fr.force = 0.0;
        fr.stuck = false;
    }
    applied.resize(frictions.len(), [0.0; 2]);
    applied.fill([0.0; 2]);
    for (slot, &fi) in active.iter().enumerate() {
        let at = n + 2 * slot;
        applied[fi] = [-forces[at], -forces[at + 1]];
        frictions[fi].force = forces[at].hypot(forces[at + 1]);
        frictions[fi].stuck = stuck[slot];
    }
    for (contact, candidate) in mindlin.iter_mut().zip(mindlin_candidates.drain(..)) {
        contact.publish(candidate);
    }
    Ok((
        NormalReport {
            newton_iterations: report.newton_iterations,
            cg_iterations: 0,
            compliance_applications: report.compliance_applications,
            displacement_residual,
            force_residual: report.normal_residual,
        },
        report,
    ))
}
