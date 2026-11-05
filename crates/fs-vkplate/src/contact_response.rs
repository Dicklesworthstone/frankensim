//! Adapt the normal-contact solve to the same frozen inverse as body prediction.

use super::{
    Body, Contact, ContactQueue, Crush, Law, contact_force, contact_phi, crush_force, crush_phi,
    viscoelastic_force, viscoelastic_phi,
};
use crate::contact_island::{Controls, SolveError, SolveReport, Workspace as IslandWorkspace};

#[path = "contact_components.rs"]
mod components;

#[derive(Clone, Default)]
pub(crate) struct Workspace {
    core: IslandWorkspace,
    pub(super) indices: Vec<usize>,
    pub(super) predicted: Vec<f64>,
    pub(super) previous: Vec<f64>,
    pub(super) previous_force: Vec<f64>,
    pub(super) diagonal_bound: Vec<f64>,
    pub(super) crush: Vec<Option<Crush>>,
    pub(super) touched: Vec<usize>,
    seen: Vec<bool>,
    components: components::Components,
    diagonal: Vec<f64>,
    scalar_force: Vec<f64>,
    scalar_pose: Vec<f64>,
    actual_pose: Vec<f64>,
}

/// Apply N^T A^-1 N without forming a contact matrix. Scatter all contact
/// forces before solving each touched body once. Repeated and opposite-signed
/// footprint entries therefore combine before the body inverse is applied.
fn apply_compliance(
    bodies: &mut [Body],
    contacts: &[Contact],
    indices: &[usize],
    touched: &[usize],
    k: f64,
    force: &[f64],
    response: &mut [f64],
) {
    for &bi in touched {
        bodies[bi].tmp.fill(0.0);
    }
    for (&ci, &f) in indices.iter().zip(force) {
        let c = &contacts[ci];
        for &(b, i, w) in &c.a {
            bodies[b].tmp[i] += w * f;
        }
        for &(b, i, w) in &c.b {
            bodies[b].tmp[i] -= w * f;
        }
    }
    for &bi in touched {
        bodies[bi].solve_force_response(k);
    }
    for (&ci, value) in indices.iter().zip(response) {
        let c = &contacts[ci];
        *value =
            c.a.iter()
                .map(|&(b, i, w)| w * bodies[b].tmp[i])
                .sum::<f64>()
                - c.b
                    .iter()
                    .map(|&(b, i, w)| w * bodies[b].tmp[i])
                    .sum::<f64>();
    }
}

pub(super) fn law(c: &Contact, previous: f64, crush: Option<&Crush>, d: f64, k: f64) -> (f64, f64) {
    let Law::Power {
        k: stiffness,
        alpha,
        mu,
    } = c.law;
    match crush {
        Some(crush) => crush_phi(d, previous, stiffness, crush, mu, k),
        None => match c.retardation_s {
            Some(retardation) => viscoelastic_phi(d, previous, stiffness, retardation, k),
            None => contact_phi(d, previous, stiffness, alpha, mu, k),
        },
    }
}

/// Pure scalar candidates: reuse the queue's interleaved Hertz arithmetic,
/// without its sequential publication. Independence must already be proved.
fn independent_forces(
    contacts: &[Contact],
    indices: &[usize],
    predicted: &[f64],
    previous: &[f64],
    diagonal: &[f64],
    crush: &[Option<Crush>],
    queue: &mut ContactQueue,
    k: f64,
    force: &mut [f64],
    pose: &mut [f64],
) {
    queue.j.clear();
    for v in [
        &mut queue.dp,
        &mut queue.cm,
        &mut queue.x0,
        &mut queue.hi,
        &mut queue.dn,
        &mut queue.kk,
        &mut queue.mu,
        &mut queue.f,
    ] {
        v.clear();
    }
    queue.fx.clear();
    for (slot, &ci) in indices.iter().enumerate() {
        let c = &contacts[ci];
        let Law::Power {
            k: stiffness,
            alpha,
            mu,
        } = c.law;
        if alpha == 1.5
            && crush[slot].is_none()
            && c.retardation_s.is_none()
            && !(predicted[slot] <= 0.0 && previous[slot] <= 0.0)
        {
            let hi = contact_phi(predicted[slot], previous[slot], stiffness, alpha, mu, k).0;
            queue.j.push(slot);
            queue.dp.push(predicted[slot]);
            queue.cm.push(diagonal[slot]);
            queue.x0.push(if c.force > 0.0 && c.force < hi {
                c.force
            } else {
                0.5 * hi
            });
            queue.hi.push(hi);
            queue.dn.push(previous[slot]);
            queue.kk.push(stiffness);
            queue.mu.push(mu);
            queue.fx.push(c.fixed_point_exit);
        } else {
            force[slot] = match crush[slot].as_ref() {
                Some(crush) => crush_force(
                    predicted[slot],
                    diagonal[slot],
                    c.force,
                    previous[slot],
                    c.law,
                    crush,
                    c.fixed_point_exit,
                    k,
                ),
                None => match c.retardation_s {
                    Some(retardation) => viscoelastic_force(
                        predicted[slot],
                        diagonal[slot],
                        c.force,
                        previous[slot],
                        stiffness,
                        retardation,
                        c.fixed_point_exit,
                        k,
                    ),
                    None => contact_force(
                        predicted[slot],
                        diagonal[slot],
                        c.force,
                        previous[slot],
                        c.law,
                        c.fixed_point_exit,
                        k,
                    ),
                },
            };
        }
    }
    queue.solve(k);
    for (&slot, &f) in queue.j.iter().zip(&queue.f) {
        force[slot] = f;
    }
    queue.j.clear(); // No pending sequential publication may escape this helper.
    for slot in 0..indices.len() {
        pose[slot] = predicted[slot] - diagonal[slot] * force[slot];
    }
}

pub(super) struct Prepared {
    independent: bool,
    nonlocal: bool,
    partitioned: bool,
}

/// Read and validate the frozen normal-contact frame without advancing any
/// physical state. The normal-only and joint friction adapters share this
/// exact footprint convention, previous pose and crushing-history selection.
pub(super) fn prepare(
    bodies: &[Body],
    contacts: &[Contact],
    queue: &mut ContactQueue,
    work: &mut Workspace,
    k: f64,
) -> Result<Prepared, SolveError> {
    if !k.is_finite() || k <= 0.0 {
        return Err(SolveError::InvalidInput(
            "time step must be positive and finite",
        ));
    }
    work.indices.clear();
    work.predicted.clear();
    work.previous.clear();
    work.previous_force.clear();
    work.diagonal_bound.clear();
    work.diagonal.clear();
    work.crush.clear();
    work.touched.clear();
    work.seen.resize(bodies.len(), false);
    work.seen.fill(false);
    work.components.reset(bodies, contacts.len());
    queue.mark.resize_with(bodies.len(), Vec::new);
    for (mark, body) in queue.mark.iter_mut().zip(bodies.iter()) {
        mark.resize(body.n, 0);
    }
    queue.new_batch();
    let mut independent = true;
    let mut nonlocal = false;
    for (ci, c) in contacts.iter().enumerate().filter(|(_, c)| c.active) {
        let slot = work.indices.len();
        let Law::Power {
            k: stiffness,
            alpha,
            mu,
        } = c.law;
        if !stiffness.is_finite()
            || stiffness < 0.0
            || !alpha.is_finite()
            || alpha <= 0.0
            || !mu.is_finite()
            || mu < 0.0
            || !c.gap.is_finite()
            || !c.force.is_finite()
            || c.force < 0.0
        {
            return Err(SolveError::InvalidInput(
                "contact law, gap or previous force is invalid",
            ));
        }
        if c.retardation_s
            .is_some_and(|a| !a.is_finite() || a < 0.0 || alpha != 1.5 || c.crush.is_some())
        {
            return Err(SolveError::InvalidInput(
                "retardation requires a nonnegative Hertz law without crushing",
            ));
        }
        if c.crush.is_some_and(|a| {
            alpha != 1.5
                || stiffness <= 0.0
                || !a.area.is_finite()
                || a.area <= 0.0
                || !a.pressure.is_finite()
                || a.pressure <= 0.0
                || !a.deepest.is_finite()
                || a.deepest < 0.0
        }) {
            return Err(SolveError::InvalidInput(
                "crushing needs positive Hertz stiffness, area and pressure",
            ));
        }
        let (mut dp, mut previous, mut norm_bound, mut exact_diagonal) = (-c.gap, -c.gap, 0.0, 0.0);
        for (terms, sign) in [(&c.a, 1.0), (&c.b, -1.0)] {
            for &(b, i, w) in terms {
                let body = bodies.get(b).ok_or(SolveError::InvalidInput(
                    "contact body index is out of range",
                ))?;
                if i >= body.n || !w.is_finite() {
                    return Err(SolveError::InvalidInput(
                        "contact node or footprint weight is invalid",
                    ));
                }
                let diagonal = body.compliance(i, k);
                if !diagonal.is_finite() || diagonal <= 0.0 {
                    return Err(SolveError::InvalidResponse);
                }
                dp += sign * w * body.wn[i];
                previous += sign * w * body.wp[i];
                // A >= D, so C_ii <= ||D^-1/2 N_i||^2. The triangle
                // bound also covers repeated signed nodes without sorting or
                // allocating. It is only a preconditioner, never the physics.
                norm_bound += w.abs() * diagonal.sqrt();
                exact_diagonal += w * w * diagonal;
                if w != 0.0 {
                    independent &=
                        work.components
                            .claim(slot, b, i, body.nonlocal_response(), queue);
                    nonlocal |= body.nonlocal_response();
                }
                if !work.seen[b] {
                    work.seen[b] = true;
                    work.touched.push(b);
                }
            }
        }
        let bound = norm_bound * norm_bound;
        if !dp.is_finite() || !previous.is_finite() || !bound.is_finite() {
            return Err(SolveError::InvalidResponse);
        }
        let crush = c.crush.map(|mut a| {
            a.deepest = a.deepest.max(previous);
            a
        });
        work.indices.push(ci);
        work.predicted.push(dp);
        work.previous.push(previous);
        work.previous_force.push(c.force);
        work.diagonal_bound.push(bound);
        work.diagonal.push(exact_diagonal);
        work.crush.push(crush);
    }
    let count = work.indices.len();
    let partitioned = !independent && work.components.finish(count) > 1;
    work.scalar_force.resize(count, 0.0);
    work.scalar_pose.resize(count, 0.0);
    work.actual_pose.resize(count, 0.0);
    Ok(Prepared {
        independent,
        nonlocal,
        partitioned,
    })
}

/// Publish all normal forces, crush memories and body displacements together.
/// A refusal changes only private numerical scratch. A caller needing friction
/// must use the joint adapter so its final normal laws remain consistent.
pub(crate) fn solve_contacts_checked(
    bodies: &mut [Body],
    contacts: &mut [Contact],
    queue: &mut ContactQueue,
    work: &mut Workspace,
    k: f64,
    controls: Controls,
) -> Result<SolveReport, SolveError> {
    controls.validate()?;
    let Prepared {
        independent,
        nonlocal,
        partitioned,
    } = prepare(bodies, contacts, queue, work, k)?;
    let count = work.indices.len();
    let Workspace {
        core,
        indices,
        predicted,
        previous,
        previous_force,
        diagonal_bound,
        crush,
        touched,
        diagonal,
        scalar_force,
        scalar_pose,
        actual_pose,
        components,
        ..
    } = work;
    // Reduced work budgets retain the checked coupled driver. The inherited
    // scalar lanes have a separate fixed 60-iteration per-contact safeguard.
    let mut scalar =
        independent && controls.max_newton_iterations >= Controls::default().max_newton_iterations;
    let mut report = SolveReport::default();
    if scalar {
        if nonlocal && count != 0 {
            // With independence proved, C * 1 is its exact diagonal. Apply
            // each touched body's inverse once for the whole scalar batch.
            scalar_force.fill(1.0);
            apply_compliance(
                bodies,
                contacts,
                indices,
                touched,
                k,
                scalar_force,
                diagonal,
            );
            report.compliance_applications += 1;
        }
        if diagonal.iter().all(|a| a.is_finite() && *a >= 0.0) {
            independent_forces(
                contacts,
                indices,
                predicted,
                previous,
                diagonal,
                crush,
                queue,
                k,
                scalar_force,
                scalar_pose,
            );
        } else {
            scalar = false;
        }
    }
    let frame = components::Frame {
        indices,
        predicted,
        previous,
        previous_force,
        diagonal_bound,
        diagonal,
        crush,
    };
    if partitioned {
        report = components.solve(
            bodies,
            contacts,
            queue,
            core,
            frame,
            k,
            controls,
            scalar_force,
            scalar_pose,
        )?;
    }
    let mut retried_scalar_components = false;
    loop {
        if !scalar && !partitioned {
            let spent = report.compliance_applications;
            let result = core.solve(
                predicted,
                previous_force,
                diagonal_bound,
                controls,
                |force, response| {
                    apply_compliance(bodies, contacts, indices, touched, k, force, response)
                },
                |slot, d| {
                    law(
                        &contacts[indices[slot]],
                        previous[slot],
                        crush[slot].as_ref(),
                        d,
                        k,
                    )
                },
            );
            report = result.map_err(|error| match error {
                SolveError::InvalidLaw { contact } => SolveError::InvalidLaw {
                    contact: indices[contact],
                },
                other => other,
            })?;
            report.compliance_applications += spent;
        }
        let (forces, penetrations) = if scalar || partitioned {
            (&scalar_force[..], &scalar_pose[..])
        } else {
            (core.forces(), core.penetrations())
        };
        let sparse = scalar && !nonlocal;
        let candidate_report = report;
        let accepted = (|| {
            if !forces.iter().all(|f| f.is_finite() && *f >= 0.0)
                || !penetrations.iter().all(|d| d.is_finite())
            {
                return Err(SolveError::InvalidResponse);
            }
            // Do not rely on whichever Krylov/trial application happened last. The
            // accepted aggregate response is computed once more for publication.
            if sparse {
                // Every nonzero node occurs exactly once. Only those diagonal
                // responses can move; leave unrelated body scratch and nodes alone.
                for (&ci, &f) in indices.iter().zip(forces) {
                    let c = &contacts[ci];
                    for (terms, sign) in [(&c.a, 1.0), (&c.b, -1.0)] {
                        for &(b, i, w) in terms {
                            if w != 0.0 {
                                bodies[b].tmp[i] = sign * w * bodies[b].compliance(i, k) * f;
                            }
                        }
                    }
                }
            } else {
                apply_compliance(bodies, contacts, indices, touched, k, forces, actual_pose);
            }
            if !indices.is_empty() {
                report.compliance_applications += 1;
            }
            if sparse {
                for &ci in indices.iter() {
                    let c = &contacts[ci];
                    for &(b, i, w) in c.a.iter().chain(&c.b) {
                        if w != 0.0 && !(bodies[b].wn[i] - bodies[b].tmp[i]).is_finite() {
                            return Err(SolveError::InvalidResponse);
                        }
                    }
                }
            } else {
                for &bi in touched.iter() {
                    let body = &bodies[bi];
                    if body
                        .wn
                        .iter()
                        .zip(&body.tmp)
                        .any(|(&next, &delta)| !(next - delta).is_finite())
                    {
                        return Err(SolveError::InvalidResponse);
                    }
                }
            }
            // Re-evaluate the actual rounded nodal pose, rather than only dp - C f.
            // Large common offsets must not hide a constitutive force mismatch.
            for (slot, &ci) in indices.iter().enumerate() {
                let c = &contacts[ci];
                let mut d = -c.gap;
                for (terms, sign) in [(&c.a, 1.0), (&c.b, -1.0)] {
                    for &(b, i, w) in terms {
                        if w != 0.0 {
                            d += sign * w * (bodies[b].wn[i] - bodies[b].tmp[i]);
                        }
                    }
                }
                let (actual_force, slope) = law(c, previous[slot], crush[slot].as_ref(), d, k);
                if !d.is_finite()
                    || !actual_force.is_finite()
                    || actual_force < 0.0
                    || !slope.is_finite()
                    || slope < 0.0
                {
                    return Err(SolveError::InvalidResponse);
                }
                let f = forces[slot];
                let force_scale = controls.force_tolerance
                    + controls.relative_tolerance * f.abs().max(actual_force.abs());
                let displacement_scale = controls.displacement_tolerance
                    + controls.relative_tolerance
                        * d.abs()
                            .max(penetrations[slot].abs())
                            .max(predicted[slot].abs());
                if !force_scale.is_finite() || !displacement_scale.is_finite() {
                    return Err(SolveError::InvalidResponse);
                }
                report.force_residual = report
                    .force_residual
                    .max((f - actual_force).abs() / force_scale);
                report.displacement_residual = report
                    .displacement_residual
                    .max((d - penetrations[slot]).abs() / displacement_scale);
                actual_pose[slot] = d;
            }
            if !report.force_residual.is_finite()
                || !report.displacement_residual.is_finite()
                || report.force_residual > 1.0
                || report.displacement_residual > 1.0
            {
                return Err(SolveError::NoConvergence {
                    newton_iterations: report.newton_iterations,
                    displacement_residual: report.displacement_residual,
                    force_residual: report.force_residual,
                });
            }
            Ok(())
        })();
        if let Err(error) = accepted {
            if scalar {
                // A scalar stopping rule is only a candidate. Preserve the
                // original warm forces and retry the simultaneous owner before
                // refusing a physically valid but more demanding tolerance.
                scalar = false;
                continue;
            }
            if partitioned
                && !retried_scalar_components
                && components.has_scalar_candidates()
                && controls.max_newton_iterations >= Controls::default().max_newton_iterations
            {
                // Discard the failed final-gate residual maxima, retaining
                // accepted component residuals and every application spent.
                let spent = report.compliance_applications;
                report = candidate_report;
                report.compliance_applications = spent;
                let retry = components.retry_scalar(
                    bodies,
                    contacts,
                    queue,
                    core,
                    frame,
                    k,
                    controls,
                    scalar_force,
                    scalar_pose,
                )?;
                components::add_report(&mut report, retry);
                retried_scalar_components = true;
                continue;
            }
            return Err(error);
        }
        if sparse {
            for &ci in indices.iter() {
                let c = &contacts[ci];
                for &(b, i, w) in c.a.iter().chain(&c.b) {
                    if w != 0.0 {
                        bodies[b].wn[i] -= bodies[b].tmp[i];
                    }
                }
            }
        } else {
            for &bi in touched.iter() {
                let body = &mut bodies[bi];
                for (next, &delta) in body.wn.iter_mut().zip(&body.tmp) {
                    *next -= delta;
                }
            }
        }
        for (slot, &ci) in indices.iter().enumerate() {
            let c = &mut contacts[ci];
            c.force = forces[slot];
            c.peak_force = c.peak_force.max(c.force);
            c.d_now = actual_pose[slot];
            c.crush = crush[slot];
        }
        return Ok(report);
    }
}

/// The exported sequential stepper is equivalent only for independent
/// footprints and diagonal body response (or one contact per nonlocal body).
/// Include inactive contacts too: a scheduled strike can activate them later.
pub(crate) fn requires_coupled_contacts(bodies: &[Body], contacts: &[Contact]) -> bool {
    let mut nodes: Vec<Vec<bool>> = bodies.iter().map(|b| vec![false; b.n]).collect();
    let mut owner = vec![usize::MAX; bodies.len()];
    for (ci, c) in contacts.iter().enumerate() {
        for &(b, i, w) in c.a.iter().chain(&c.b) {
            if w == 0.0 {
                continue;
            }
            if b >= bodies.len() || i >= bodies[b].n || !w.is_finite() {
                return true;
            }
            if nodes[b][i]
                || (bodies[b].nonlocal_response() && owner[b] != usize::MAX && owner[b] != ci)
            {
                return true;
            }
            nodes[b][i] = true;
            owner[b] = ci;
        }
    }
    false
}
