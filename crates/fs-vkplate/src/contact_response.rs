//! Adapt the normal-contact solve to the same frozen inverse as body prediction.

use super::{Body, Contact, ContactQueue, Crush, Law, contact_phi, crush_phi, viscoelastic_phi};
use crate::contact_island::{Controls, SolveError, SolveReport, Workspace as IslandWorkspace};

#[derive(Clone, Default)]
pub(crate) struct Workspace {
    core: IslandWorkspace,
    indices: Vec<usize>,
    predicted: Vec<f64>,
    previous: Vec<f64>,
    previous_force: Vec<f64>,
    diagonal_bound: Vec<f64>,
    crush: Vec<Option<Crush>>,
    touched: Vec<usize>,
    seen: Vec<bool>,
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

fn law(c: &Contact, previous: f64, crush: Option<&Crush>, d: f64, k: f64) -> (f64, f64) {
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

/// Publish all normal forces, crush memories and body displacements together.
/// A refusal changes only private numerical scratch. Friction is a subsequent,
/// separate stage and is not included in this normal-contact consistency claim.
pub(crate) fn solve_contacts_checked(
    bodies: &mut [Body],
    contacts: &mut [Contact],
    _queue: &mut ContactQueue,
    work: &mut Workspace,
    k: f64,
    controls: Controls,
) -> Result<SolveReport, SolveError> {
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
    work.crush.clear();
    work.touched.clear();
    work.seen.resize(bodies.len(), false);
    work.seen.fill(false);
    for (ci, c) in contacts.iter().enumerate().filter(|(_, c)| c.active) {
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
        let (mut dp, mut previous, mut norm_bound) = (-c.gap, -c.gap, 0.0);
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
        work.crush.push(crush);
    }
    let Workspace {
        core,
        indices,
        predicted,
        previous,
        previous_force,
        diagonal_bound,
        crush,
        touched,
        ..
    } = work;
    let result = core.solve(
        predicted,
        previous_force,
        diagonal_bound,
        controls,
        |force, response| apply_compliance(bodies, contacts, indices, touched, k, force, response),
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
    let mut report = result.map_err(|error| match error {
        SolveError::InvalidLaw { contact } => SolveError::InvalidLaw {
            contact: indices[contact],
        },
        other => other,
    })?;
    // Do not rely on whichever Krylov/trial application happened last. The
    // accepted aggregate response is computed once more for publication.
    apply_compliance(
        bodies,
        contacts,
        indices,
        touched,
        k,
        core.forces(),
        previous_force,
    );
    if !indices.is_empty() {
        report.compliance_applications += 1;
    }
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
    // Re-evaluate the actual rounded nodal pose, rather than only dp - C f.
    // Large common offsets must not hide a constitutive force mismatch.
    for (slot, &ci) in indices.iter().enumerate() {
        let c = &contacts[ci];
        let mut d = -c.gap;
        for (terms, sign) in [(&c.a, 1.0), (&c.b, -1.0)] {
            for &(b, i, w) in terms {
                d += sign * w * (bodies[b].wn[i] - bodies[b].tmp[i]);
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
        let f = core.forces()[slot];
        let force_scale = controls.force_tolerance
            + controls.relative_tolerance * f.abs().max(actual_force.abs());
        let displacement_scale = controls.displacement_tolerance
            + controls.relative_tolerance
                * d.abs()
                    .max(core.penetrations()[slot].abs())
                    .max(predicted[slot].abs());
        report.force_residual = report
            .force_residual
            .max((f - actual_force).abs() / force_scale);
        report.displacement_residual = report
            .displacement_residual
            .max((d - core.penetrations()[slot]).abs() / displacement_scale);
        previous_force[slot] = d;
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
    for &bi in touched.iter() {
        let body = &mut bodies[bi];
        for (next, &delta) in body.wn.iter_mut().zip(&body.tmp) {
            *next -= delta;
        }
    }
    for (slot, &ci) in indices.iter().enumerate() {
        let c = &mut contacts[ci];
        c.force = core.forces()[slot];
        c.peak_force = c.peak_force.max(c.force);
        c.d_now = previous_force[slot];
        c.crush = crush[slot];
    }
    Ok(report)
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
