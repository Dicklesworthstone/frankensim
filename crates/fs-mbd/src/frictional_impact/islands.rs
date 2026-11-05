//! Independent simultaneous contacts without a dense matrix for the whole scene.
//!
//! Only shared dynamic bodies connect islands. A prescribed support supplies
//! velocity and work but has no reciprocal state, so a common floor or belt
//! does not couple otherwise independent bodies. Each island uses the original
//! native impulse owner, including its complete constitutive and body audits.
//! All islands remain private until their aggregate balances also pass.

use super::{
    ImpactResponse, PhysicalImpulseView, RigidFrictionalContact, RigidFrictionalImpactError,
    RigidFrictionalImpactReport, RigidImpactControls, RigidImpactCounterpart,
    apply_rigid_frictional_impacts,
};
use crate::{MassProperties, RigidBodyState};
use fs_contact::frictional_multicontact::elliptic::{
    EllipticFrictionalContactError, EllipticFrictionalContactResult,
};
use fs_contact::frictional_multicontact::{FrictionalContactError, FrictionalContactResult};
use fs_exec::CancelGate;
use fs_math::dd::Dd;

/// One original native solve and its mapping back to the supplied scene.
#[derive(Clone, Debug)]
pub struct ContactIsland<R = FrictionalContactResult> {
    /// Original body indices, in increasing order. Every body belongs to at
    /// most one island; bodies with no supplied contact are absent.
    pub body_indices: Vec<usize>,
    /// Original contact indices in the native report's stable contact-ID order.
    pub contact_indices: Vec<usize>,
    /// Constitutive determinacy, work and native balances for this island.
    /// A velocity-only result never becomes a unique local impulse or heat
    /// allocation merely because it is part of an independent scene.
    pub impact: RigidFrictionalImpactReport<R>,
}

/// Physical scene response with every original island witness retained.
#[derive(Clone, Debug)]
pub struct IslandContactResult<R = FrictionalContactResult> {
    /// Ordered by each island's minimum stable contact ID.
    pub islands: Vec<ContactIsland<R>>,
    /// Canonical global contact-ID order, matching the outer report's IDs.
    pub impulses: Vec<[f64; 3]>,
    /// Actual endpoint contact velocities, in the same global order.
    pub relative_velocities_after: Vec<[f64; 3]>,
    /// Sum of admitted relative impulse work, J with SI inputs.
    pub relative_work: f64,
    /// Total mechanical loss, with no inferred heat partition.
    pub dissipated_energy: f64,
}

impl<R> ImpactResponse for IslandContactResult<R> {
    fn physical(&self) -> PhysicalImpulseView<'_> {
        PhysicalImpulseView {
            impulses: &self.impulses,
            relative_velocities_after: &self.relative_velocities_after,
            relative_work: self.relative_work,
            dissipated_energy: self.dissipated_energy,
        }
    }
}

/// Failure of the bounded island decomposition or one original native solve.
#[derive(Debug)]
pub enum IslandContactError<R = FrictionalContactResult, E = FrictionalContactError> {
    /// A connected set cannot be split without discarding physical coupling.
    ContactLimit {
        /// Smallest stable contact ID in the connected set.
        first_contact_id: u64,
        /// Supplied contacts in that set, including inactive witnesses.
        contacts: usize,
        /// Existing constitutive solver's per-island limit.
        maximum: usize,
    },
    /// `None` identifies the empty native configuration admission; otherwise
    /// the ID names the failed island. No earlier island has been published.
    Solve {
        /// Smallest stable contact ID, absent for configuration admission.
        first_contact_id: Option<u64>,
        /// Unchanged native refusal, including the original physical receipt.
        source: Box<RigidFrictionalImpactError<R, E>>,
    },
}

impl<R: core::fmt::Debug, E: core::fmt::Debug> core::fmt::Display for IslandContactError<R, E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid contact islands: {self:?}")
    }
}
impl<R: core::fmt::Debug, E: core::fmt::Debug> std::error::Error for IslandContactError<R, E> {}

/// Whole-scene native audit retaining each independent constitutive result.
pub type IslandImpactReport<R = FrictionalContactResult> =
    RigidFrictionalImpactReport<IslandContactResult<R>>;
/// Scene admission, constituent native failure, cancellation or global audit.
pub type IslandImpactError<R = FrictionalContactResult, E = FrictionalContactError> =
    RigidFrictionalImpactError<IslandContactResult<R>, IslandContactError<R, E>>;
/// Native island audit for directional Coulomb friction.
pub type EllipticIslandImpactReport = IslandImpactReport<EllipticFrictionalContactResult>;
/// Directional island refusal with the original ellipse solve retained.
pub type EllipticIslandImpactError =
    IslandImpactError<EllipticFrictionalContactResult, EllipticFrictionalContactError>;

/// Apply independent native normal/Coulomb events as one atomic scene update.
///
/// `controls.maximum_contacts` bounds the whole scene, while
/// `controls.contact.maximum_contacts` bounds each connected island (21 by
/// default). Explicitly larger limits also require the constitutive owner's
/// scalar matrix budget. More than 63 active scalar rows require its
/// positive-definite determinacy admission; they do not broaden its bounded
/// redundant-contact path. The latter cap is not
/// a total scene limit. Other constitutive work limits also apply per island;
/// the scene contact cap bounds their total number. Assembly needs O(B+C)
/// graph storage plus the sum of the original island receipts, and builds
/// only the individual dense matrices, never a global `(3C)^2` matrix.
///
/// All supplied witnesses connect their named dynamic bodies, including
/// geometrically inactive witnesses whose endpoint velocities can change.
/// This conservative partition preserves the native candidate/geometry checks;
/// it does not discover contacts or silently drop inactive input constraints.
/// No connected island is approximated by sequential pairwise collisions.
///
/// Islands are solved in minimum-contact-ID order. Every native balance gate
/// and the same global tolerance apply before final publication, so numerical
/// allowances do not grow with the number of independent impacts. Cancellation
/// or any late failure leaves every supplied body unchanged.
pub fn apply_rigid_contact_islands(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contacts: &[RigidFrictionalContact],
    controls: RigidImpactControls,
    gate: &CancelGate,
) -> Result<IslandImpactReport, IslandImpactError> {
    apply_islands_with(
        bodies,
        properties,
        contacts,
        controls,
        gate,
        apply_rigid_frictional_impacts,
    )
}

/// Directional friction uses the same exact dynamic-body decomposition.
pub fn apply_rigid_elliptic_contact_islands(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contacts: &[super::elliptic::RigidEllipticContact],
    controls: RigidImpactControls,
    gate: &CancelGate,
) -> Result<EllipticIslandImpactReport, EllipticIslandImpactError> {
    apply_islands_with(
        bodies,
        properties,
        contacts,
        controls,
        gate,
        super::elliptic::apply_rigid_elliptic_impacts,
    )
}

struct Component {
    bodies: Vec<usize>,
    contacts: Vec<usize>,
}

#[allow(clippy::too_many_lines)]
fn apply_islands_with<F: Copy, R: ImpactResponse, E>(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contacts: &[RigidFrictionalContact<F>],
    controls: RigidImpactControls,
    gate: &CancelGate,
    mut solve: impl FnMut(
        &mut [RigidBodyState],
        &[MassProperties],
        &[RigidFrictionalContact<F>],
        RigidImpactControls,
        &CancelGate,
    )
        -> Result<RigidFrictionalImpactReport<R>, RigidFrictionalImpactError<R, E>>,
) -> Result<IslandImpactReport<R>, IslandImpactError<R, E>> {
    let invalid = RigidFrictionalImpactError::Invalid;
    super::poll_with(gate)?;
    if bodies.len() != properties.len()
        || bodies.len() > controls.maximum_bodies
        || contacts.len() > controls.maximum_contacts
    {
        return Err(invalid(
            "matching bounded scene bodies and contacts required",
        ));
    }
    // Reuse native control and constitutive admission, including empty-scene
    // behavior, before allocating the graph. This builds only an empty matrix.
    solve(&mut [], &[], &[], controls, gate).map_err(|source| {
        RigidFrictionalImpactError::Contact(IslandContactError::Solve {
            first_contact_id: None,
            source: Box::new(source),
        })
    })?;
    let mut order: Vec<_> = (0..contacts.len()).collect();
    order.sort_unstable_by_key(|&i| contacts[i].id);
    let mut parent: Vec<_> = (0..bodies.len()).collect();
    let mut rank = vec![0_u8; bodies.len()];
    for (slot, &index) in order.iter().enumerate() {
        super::poll_with(gate)?;
        let c = contacts[index];
        if c.point.body >= bodies.len() || (slot > 0 && c.id == contacts[order[slot - 1]].id) {
            return Err(invalid(
                "valid body indices and distinct contact identities required",
            ));
        }
        if let RigidImpactCounterpart::Body(other) = c.other {
            if other.body >= bodies.len() || other.body == c.point.body {
                return Err(invalid("counterpart must name another valid dynamic body"));
            }
            join(&mut parent, &mut rank, c.point.body, other.body);
        }
    }
    let mut component_of_root = vec![usize::MAX; bodies.len()];
    let mut components: Vec<Component> = Vec::new();
    for &index in &order {
        super::poll_with(gate)?;
        let root = find(&mut parent, contacts[index].point.body);
        let slot = &mut component_of_root[root];
        if *slot == usize::MAX {
            *slot = components.len();
            components.push(Component {
                bodies: Vec::new(),
                contacts: Vec::new(),
            });
        }
        components[*slot].contacts.push(index);
    }
    let mut untouched_energy = Dd::ZERO;
    for (index, (&body, &property)) in bodies.iter().zip(properties).enumerate() {
        super::poll_with(gate)?;
        let slot = component_of_root[find(&mut parent, index)];
        if slot != usize::MAX {
            components[slot].bodies.push(index);
        } else {
            // Contacted bodies are checked by their native solve. Admit each
            // unchanged body once, without re-evaluating every island's energies.
            untouched_energy = untouched_energy + Dd::from_f64(body.kinetic_energy(property)?);
        }
    }
    for component in &components {
        super::poll_with(gate)?;
        if component.contacts.len() > controls.contact.maximum_contacts {
            return Err(RigidFrictionalImpactError::Contact(
                IslandContactError::ContactLimit {
                    first_contact_id: contacts[component.contacts[0]].id,
                    contacts: component.contacts.len(),
                    maximum: controls.contact.maximum_contacts,
                },
            ));
        }
    }
    let mut canonical_index = vec![0; contacts.len()];
    for (slot, &index) in order.iter().enumerate() {
        canonical_index[index] = slot;
    }
    let mut local_index = vec![0; bodies.len()];
    let mut candidate = bodies.to_vec();
    let mut islands = Vec::with_capacity(components.len());
    let mut impulses = vec![[0.0; 3]; contacts.len()];
    let mut velocities = vec![[0.0; 3]; contacts.len()];
    let mut sums = Sums::default();
    for component in components {
        super::poll_with(gate)?;
        let mut local_bodies = Vec::with_capacity(component.bodies.len());
        let mut local_properties = Vec::with_capacity(component.bodies.len());
        for (local, &original) in component.bodies.iter().enumerate() {
            local_index[original] = local;
            local_bodies.push(bodies[original]);
            local_properties.push(properties[original]);
        }
        let local_contacts: Vec<_> = component
            .contacts
            .iter()
            .map(|&original| {
                let mut c = contacts[original];
                c.point.body = local_index[c.point.body];
                if let RigidImpactCounterpart::Body(ref mut other) = c.other {
                    other.body = local_index[other.body];
                }
                c
            })
            .collect();
        let impact = solve(
            &mut local_bodies,
            &local_properties,
            &local_contacts,
            controls,
            gate,
        )
        .map_err(|source| {
            RigidFrictionalImpactError::Contact(IslandContactError::Solve {
                first_contact_id: Some(local_contacts[0].id),
                source: Box::new(source),
            })
        })?;
        super::poll_with(gate)?;
        let physical = impact.contact.physical();
        for (local, &original) in component.contacts.iter().enumerate() {
            let slot = canonical_index[original];
            impulses[slot] = physical.impulses[local];
            velocities[slot] = physical.relative_velocities_after[local];
        }
        for (&original, &body) in component.bodies.iter().zip(&local_bodies) {
            candidate[original] = body;
        }
        sums.add(&impact);
        islands.push(ContactIsland {
            body_indices: component.bodies,
            contact_indices: component.contacts,
            impact,
        });
    }
    let report = RigidFrictionalImpactReport {
        contact_ids: order.iter().map(|&i| contacts[i].id).collect(),
        contact: IslandContactResult {
            islands,
            impulses,
            relative_velocities_after: velocities,
            relative_work: sums.relative_work.to_f64(),
            dissipated_energy: sums.loss.to_f64(),
        },
        kinetic_before_j: (untouched_energy + sums.before).to_f64(),
        kinetic_after_j: (untouched_energy + sums.after).to_f64(),
        kinetic_change_j: sums.kinetic_change.to_f64(),
        support_work_j: sums.support_work.to_f64(),
        relative_impulse_work_j: sums.native_work.to_f64(),
        work_energy_residual_j: (sums.kinetic_change - sums.support_work - sums.native_work)
            .to_f64(),
        total_energy_residual_j: (sums.kinetic_change - sums.support_work + sums.loss).to_f64(),
        support_impulse_ns: super::rounded_vector(sums.support_impulse),
        support_angular_impulse_nms: super::rounded_vector(sums.support_angular),
        linear_momentum_residual_ns: super::rounded_vector(sums.linear_residual),
        angular_momentum_residual_nms: super::rounded_vector(sums.angular_residual),
        endpoint_velocity_residual_mps: sums.endpoint_error,
    };
    if [
        report.kinetic_before_j,
        report.kinetic_after_j,
        report.kinetic_change_j,
        report.support_work_j,
        report.relative_impulse_work_j,
        report.contact.relative_work,
        report.contact.dissipated_energy,
        report.work_energy_residual_j,
        report.total_energy_residual_j,
        super::magnitude(report.support_impulse_ns),
        super::magnitude(report.support_angular_impulse_nms),
        super::magnitude(report.linear_momentum_residual_ns),
        super::magnitude(report.angular_momentum_residual_nms),
    ]
    .iter()
    .any(|value| !value.is_finite())
        || report.work_energy_residual_j.abs() > controls.energy_tolerance_j
        || report.total_energy_residual_j.abs() > controls.energy_tolerance_j
        || (sums.native_work - sums.relative_work).to_f64().abs() > controls.energy_tolerance_j
        || super::magnitude(report.linear_momentum_residual_ns)
            > controls.linear_momentum_tolerance_ns
        || super::magnitude(report.angular_momentum_residual_nms)
            > controls.angular_momentum_tolerance_nms
        || report.endpoint_velocity_residual_mps > controls.velocity_tolerance_mps
    {
        return Err(RigidFrictionalImpactError::Audit(Box::new(report)));
    }
    super::poll_with(gate)?;
    bodies.copy_from_slice(&candidate);
    Ok(report)
}

#[derive(Default)]
struct Sums {
    before: Dd,
    after: Dd,
    kinetic_change: Dd,
    support_work: Dd,
    native_work: Dd,
    relative_work: Dd,
    loss: Dd,
    support_impulse: [Dd; 3],
    support_angular: [Dd; 3],
    linear_residual: [Dd; 3],
    angular_residual: [Dd; 3],
    endpoint_error: f64,
}

impl Sums {
    fn add<R: ImpactResponse>(&mut self, impact: &RigidFrictionalImpactReport<R>) {
        self.before = self.before + Dd::from_f64(impact.kinetic_before_j);
        self.after = self.after + Dd::from_f64(impact.kinetic_after_j);
        self.kinetic_change = self.kinetic_change + Dd::from_f64(impact.kinetic_change_j);
        self.support_work = self.support_work + Dd::from_f64(impact.support_work_j);
        self.native_work = self.native_work + Dd::from_f64(impact.relative_impulse_work_j);
        self.relative_work =
            self.relative_work + Dd::from_f64(impact.contact.physical().relative_work);
        self.loss = self.loss + Dd::from_f64(impact.contact.physical().dissipated_energy);
        super::add_vector(&mut self.support_impulse, impact.support_impulse_ns);
        super::add_vector(
            &mut self.support_angular,
            impact.support_angular_impulse_nms,
        );
        super::add_vector(
            &mut self.linear_residual,
            impact.linear_momentum_residual_ns,
        );
        super::add_vector(
            &mut self.angular_residual,
            impact.angular_momentum_residual_nms,
        );
        self.endpoint_error = self
            .endpoint_error
            .max(impact.endpoint_velocity_residual_mps);
    }
}

fn find(parent: &mut [usize], mut node: usize) -> usize {
    while parent[node] != node {
        parent[node] = parent[parent[node]];
        node = parent[node];
    }
    node
}

fn join(parent: &mut [usize], rank: &mut [u8], a: usize, b: usize) {
    let mut a = find(parent, a);
    let mut b = find(parent, b);
    if a == b {
        return;
    }
    if rank[a] < rank[b] || (rank[a] == rank[b] && a > b) {
        core::mem::swap(&mut a, &mut b);
    }
    parent[b] = a;
    if rank[a] == rank[b] {
        rank[a] += 1;
    }
}
