//! Stable response-graph components and staged contact candidates.

use super::{Body, Contact, ContactQueue, Crush, independent_forces, law};
use crate::contact_island::{Controls, SolveError, SolveReport, Workspace as IslandWorkspace};

#[derive(Clone, Copy)]
pub(super) struct Frame<'a> {
    pub indices: &'a [usize],
    pub predicted: &'a [f64],
    pub previous: &'a [f64],
    pub previous_force: &'a [f64],
    pub diagonal_bound: &'a [f64],
    pub diagonal: &'a [f64],
    pub crush: &'a [Option<Crush>],
}

#[derive(Clone, Default)]
pub(super) struct Components {
    parent: Vec<usize>,
    size: Vec<usize>,
    repeated: Vec<bool>,
    node_owner: Vec<Vec<usize>>,
    body_owner: Vec<usize>,
    group_of: Vec<usize>,
    root_group: Vec<usize>,
    counts: Vec<usize>,
    offsets: Vec<usize>,
    slots: Vec<usize>,
    scalar_slots: Vec<usize>,
    group: Group,
}

fn root(parent: &mut [usize], mut i: usize) -> usize {
    while parent[i] != i {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    i
}

impl Components {
    pub(super) fn reset(&mut self, bodies: &[Body], count: usize) {
        self.parent.resize(count, 0);
        for (i, parent) in self.parent.iter_mut().enumerate() {
            *parent = i;
        }
        self.size.resize(count, 1);
        self.size.fill(1);
        self.repeated.resize(count, false);
        self.repeated.fill(false);
        self.body_owner.resize(bodies.len(), usize::MAX);
        self.body_owner.fill(usize::MAX);
        self.node_owner.resize_with(bodies.len(), Vec::new);
        for (owners, body) in self.node_owner.iter_mut().zip(bodies) {
            owners.resize(body.n, 0);
        }
        self.scalar_slots.clear();
    }

    fn join(&mut self, a: usize, b: usize) {
        let (mut a, mut b) = (root(&mut self.parent, a), root(&mut self.parent, b));
        if a == b {
            return;
        }
        if self.size[a] < self.size[b] || (self.size[a] == self.size[b] && a > b) {
            std::mem::swap(&mut a, &mut b);
        }
        self.parent[b] = a;
        self.size[a] += self.size[b];
    }

    /// Claim one validated nonzero footprint term. False also covers repeated
    /// terms within one contact, which must retain aggregate response physics.
    pub(super) fn claim(
        &mut self,
        slot: usize,
        b: usize,
        i: usize,
        nonlocal: bool,
        queue: &mut ContactQueue,
    ) -> bool {
        let mut independent = true;
        if queue.mark[b][i] == queue.batch {
            let other = self.node_owner[b][i];
            self.repeated[slot] |= other == slot;
            self.join(slot, other);
            independent = false;
        } else {
            queue.mark[b][i] = queue.batch;
            self.node_owner[b][i] = slot;
        }
        if nonlocal && self.body_owner[b] != usize::MAX && self.body_owner[b] != slot {
            self.join(slot, self.body_owner[b]);
            independent = false;
        }
        self.body_owner[b] = slot;
        independent
    }

    /// Component order follows first contact occurrence; contacts inside each
    /// component retain input order, independent of union-tree shape.
    pub(super) fn finish(&mut self, count: usize) -> usize {
        self.group_of.resize(count, 0);
        self.root_group.resize(count, usize::MAX);
        self.root_group.fill(usize::MAX);
        self.counts.clear();
        for slot in 0..count {
            let r = root(&mut self.parent, slot);
            if self.root_group[r] == usize::MAX {
                self.root_group[r] = self.counts.len();
                self.counts.push(0);
            }
            let group = self.root_group[r];
            self.group_of[slot] = group;
            self.counts[group] += 1;
        }
        self.offsets.clear();
        self.offsets.push(0);
        for &count in &self.counts {
            self.offsets.push(self.offsets.last().unwrap() + count);
        }
        self.counts.fill(0);
        self.slots.resize(count, 0);
        for slot in 0..count {
            let group = self.group_of[slot];
            self.slots[self.offsets[group] + self.counts[group]] = slot;
            self.counts[group] += 1;
        }
        self.counts.len()
    }

    pub(super) fn has_scalar_candidates(&self) -> bool {
        !self.scalar_slots.is_empty()
    }

    pub(super) fn solve(
        &mut self,
        bodies: &mut [Body],
        contacts: &[Contact],
        queue: &mut ContactQueue,
        core: &mut IslandWorkspace,
        frame: Frame<'_>,
        k: f64,
        controls: Controls,
        forces: &mut [f64],
        poses: &mut [f64],
    ) -> Result<SolveReport, SolveError> {
        self.scalar_slots.clear();
        for group in 0..self.counts.len() {
            let slots = &self.slots[self.offsets[group]..self.offsets[group + 1]];
            if slots.len() == 1 && !self.repeated[slots[0]] {
                self.scalar_slots.push(slots[0]);
            }
        }
        let mut report = SolveReport::default();
        if !self.scalar_slots.is_empty() {
            add_report(
                &mut report,
                self.group.solve(
                    &self.scalar_slots,
                    frame,
                    bodies,
                    contacts,
                    queue,
                    core,
                    k,
                    controls,
                    true,
                    forces,
                    poses,
                )?,
            );
        }
        for group in 0..self.counts.len() {
            let slots = &self.slots[self.offsets[group]..self.offsets[group + 1]];
            if slots.len() == 1 && !self.repeated[slots[0]] {
                continue;
            }
            add_report(
                &mut report,
                self.group.solve(
                    slots, frame, bodies, contacts, queue, core, k, controls, false, forces, poses,
                )?,
            );
        }
        Ok(report)
    }

    pub(super) fn retry_scalar(
        &mut self,
        bodies: &mut [Body],
        contacts: &[Contact],
        queue: &mut ContactQueue,
        core: &mut IslandWorkspace,
        frame: Frame<'_>,
        k: f64,
        controls: Controls,
        forces: &mut [f64],
        poses: &mut [f64],
    ) -> Result<SolveReport, SolveError> {
        self.group.solve(
            &self.scalar_slots,
            frame,
            bodies,
            contacts,
            queue,
            core,
            k,
            controls,
            false,
            forces,
            poses,
        )
    }
}

pub(super) fn add_report(total: &mut SolveReport, part: SolveReport) {
    total.newton_iterations += part.newton_iterations;
    total.cg_iterations += part.cg_iterations;
    total.compliance_applications += part.compliance_applications;
    total.displacement_residual = total.displacement_residual.max(part.displacement_residual);
    total.force_residual = total.force_residual.max(part.force_residual);
}

#[derive(Clone, Default)]
struct Group {
    indices: Vec<usize>,
    predicted: Vec<f64>,
    previous: Vec<f64>,
    previous_force: Vec<f64>,
    bound: Vec<f64>,
    diagonal: Vec<f64>,
    crush: Vec<Option<Crush>>,
    force: Vec<f64>,
    pose: Vec<f64>,
    nonlocal: Vec<usize>,
    nodes: Vec<(usize, usize)>,
    seen: Vec<bool>,
}

impl Group {
    fn prepare(
        &mut self,
        slots: &[usize],
        frame: Frame<'_>,
        bodies: &[Body],
        contacts: &[Contact],
        queue: &mut ContactQueue,
    ) {
        self.indices.clear();
        for v in [
            &mut self.predicted,
            &mut self.previous,
            &mut self.previous_force,
            &mut self.bound,
            &mut self.diagonal,
        ] {
            v.clear();
        }
        self.crush.clear();
        self.nonlocal.clear();
        self.nodes.clear();
        self.seen.resize(bodies.len(), false);
        self.seen.fill(false);
        queue.new_batch();
        for &slot in slots {
            let ci = frame.indices[slot];
            self.indices.push(ci);
            self.predicted.push(frame.predicted[slot]);
            self.previous.push(frame.previous[slot]);
            self.previous_force.push(frame.previous_force[slot]);
            self.bound.push(frame.diagonal_bound[slot]);
            self.diagonal.push(frame.diagonal[slot]);
            self.crush.push(frame.crush[slot]);
            let c = &contacts[ci];
            for &(b, i, w) in c.a.iter().chain(&c.b) {
                if w == 0.0 {
                    continue;
                }
                if bodies[b].nonlocal_response() {
                    if !self.seen[b] {
                        self.nonlocal.push(b);
                        self.seen[b] = true;
                    }
                } else if queue.mark[b][i] != queue.batch {
                    queue.mark[b][i] = queue.batch;
                    self.nodes.push((b, i));
                }
            }
        }
        self.force.resize(slots.len(), 0.0);
        self.pose.resize(slots.len(), 0.0);
    }

    fn solve(
        &mut self,
        slots: &[usize],
        frame: Frame<'_>,
        bodies: &mut [Body],
        contacts: &[Contact],
        queue: &mut ContactQueue,
        core: &mut IslandWorkspace,
        k: f64,
        controls: Controls,
        mut scalar: bool,
        forces: &mut [f64],
        poses: &mut [f64],
    ) -> Result<SolveReport, SolveError> {
        self.prepare(slots, frame, bodies, contacts, queue);
        scalar &= controls.max_newton_iterations >= Controls::default().max_newton_iterations;
        let mut report = SolveReport::default();
        if scalar {
            if !self.nonlocal.is_empty() {
                self.force.fill(1.0);
                apply(
                    bodies,
                    contacts,
                    &self.indices,
                    &self.nonlocal,
                    &self.nodes,
                    k,
                    &self.force,
                    &mut self.diagonal,
                );
                report.compliance_applications += 1;
            }
            if self.diagonal.iter().all(|a| a.is_finite() && *a >= 0.0) {
                independent_forces(
                    contacts,
                    &self.indices,
                    &self.predicted,
                    &self.previous,
                    &self.diagonal,
                    &self.crush,
                    queue,
                    k,
                    &mut self.force,
                    &mut self.pose,
                );
            } else {
                scalar = false;
            }
        }
        if !scalar {
            let spent = report.compliance_applications;
            report = core
                .solve(
                    &self.predicted,
                    &self.previous_force,
                    &self.bound,
                    controls,
                    |force, response| {
                        apply(
                            bodies,
                            contacts,
                            &self.indices,
                            &self.nonlocal,
                            &self.nodes,
                            k,
                            force,
                            response,
                        )
                    },
                    |local, d| {
                        law(
                            &contacts[self.indices[local]],
                            self.previous[local],
                            self.crush[local].as_ref(),
                            d,
                            k,
                        )
                    },
                )
                .map_err(|error| match error {
                    SolveError::InvalidLaw { contact } => SolveError::InvalidLaw {
                        contact: self.indices[contact],
                    },
                    other => other,
                })?;
            report.compliance_applications += spent;
            self.force.copy_from_slice(core.forces());
            self.pose.copy_from_slice(core.penetrations());
        }
        for (local, &slot) in slots.iter().enumerate() {
            forces[slot] = self.force[local];
            poses[slot] = self.pose[local];
        }
        Ok(report)
    }
}

/// Component-local response: sparse diagonal supports may share a body with
/// another component, so never scan or multiply the entire diagonal body.
fn apply(
    bodies: &mut [Body],
    contacts: &[Contact],
    indices: &[usize],
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
    for (&ci, &f) in indices.iter().zip(force) {
        let c = &contacts[ci];
        for (terms, sign) in [(&c.a, 1.0), (&c.b, -1.0)] {
            for &(b, i, w) in terms {
                if w != 0.0 {
                    bodies[b].tmp[i] += sign * w * f;
                }
            }
        }
    }
    for &b in nonlocal {
        bodies[b].solve_force_response(k);
    }
    for &(b, i) in nodes {
        bodies[b].tmp[i] *= bodies[b].compliance(i, k);
    }
    for (&ci, value) in indices.iter().zip(response) {
        let c = &contacts[ci];
        *value = 0.0;
        for (terms, sign) in [(&c.a, 1.0), (&c.b, -1.0)] {
            for &(b, i, w) in terms {
                if w != 0.0 {
                    *value += sign * w * bodies[b].tmp[i];
                }
            }
        }
    }
}
