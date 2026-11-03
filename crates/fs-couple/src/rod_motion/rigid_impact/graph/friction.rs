//! Coupled Coulomb impulses on the original mixed rod/rigid contact graph.
//!
//! Endpoint geometry, signed mobility and native momentum publication belong to
//! the common graph owner. This adapter only selects the original three-axis
//! contact law; normal and tangent impulses share the same graph response.

use super::super::friction::{RodRigidFrictionalControls, RodRigidFrictionalError, contact_frame};
use super::super::{RodRigidImpactError, project_velocity};
use super::{
    RodRigidGraphBudget, RodRigidGraphContact, RodRigidGraphTransfer, RodRigidParticipant,
    validate_graph_budget,
};
use crate::rod_motion::{RodMotion, RodMotionInertia, checkpoint};
use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactResult, NormalVelocityLaw, solve_frictional_contacts,
};
use fs_exec::CancelGate;
use fs_mbd::RigidBodyState;

/// One oriented graph edge and the coefficient resolved for this event.
#[derive(Clone, Copy, Debug)]
pub struct RodRigidFrictionalGraphContact {
    pub geometry: RodRigidGraphContact,
    /// World direction projected into the admitted contact plane.
    pub tangent_world: [f64; 3],
    /// Finite, nonnegative, event-local Coulomb coefficient.
    pub friction_coefficient: f64,
}

/// The original contact solution and its actual native graph response.
/// A redundant impulse representative does not authorize unique local loads or
/// per-edge heat allocations; the original determinacy/work witness is retained.
#[derive(Clone, Debug)]
pub struct RodRigidFrictionalGraphReport {
    pub contact: FrictionalContactResult,
    /// Complete row-major 3N-by-3N mobility in original contact order [kg^-1].
    /// Each edge uses its normal, tangent-1 and tangent-2 coordinates.
    pub delassus_kg_inverse: Vec<f64>,
    pub gaps_m: Vec<f64>,
    pub frames_world: Vec<[[f64; 3]; 3]>,
    pub relative_velocities_before_m_s: Vec<[f64; 3]>,
    /// Actual endpoint velocities, native receipts and aggregate physical work.
    pub transfer: RodRigidGraphTransfer<3>,
}

impl<I: RodMotionInertia> RodMotion<I> {
    /// Resolve simultaneous normal and tangent impulses on a mixed graph.
    ///
    /// Edges can join rod material attachments to finite rigid bodies, two
    /// distinct rigid bodies, or two distinct rod nodes. Every shared native
    /// participant contributes its full signed cross mobility, including all
    /// normal/tangent terms, original central tensors, eccentric arms and
    /// admitted stationary supports. The same world impulse acts positively at
    /// the edge's positive endpoint and negatively at its negative endpoint.
    ///
    /// The original non-associated Signorini--Coulomb owner solves the entire
    /// graph and retains its passivity, determinacy and total-work refusals.
    /// The common graph owner applies one original rod batch and one resultant
    /// update per affected body, then checks every actual endpoint triple,
    /// support-inclusive momentum and total kinetic work before returning any
    /// participant. Tangential impulses require actual witness collocation.
    /// Unchanged participants retain their native states and omit receipts.
    ///
    /// Contact/participant budgets are checked before graph allocation. The
    /// original contact limit is 21 edges / 63 scalar rows. Matrix assembly is
    /// O(edges^2) with O(edges) original inverse-inertia actions; native staging
    /// is O(rod nodes + bodies + edges), with no per-edge trajectory clones.
    /// This frozen-configuration event preserves clocks and material history.
    /// It does not discover collisions, evolve persistent friction, or allocate
    /// mechanical loss among thermal or other material destinations.
    pub fn apply_rigid_frictional_contact_graph(
        &self,
        bodies: &[RodRigidParticipant],
        witnesses: &[RodRigidFrictionalGraphContact],
        controls: RodRigidFrictionalControls,
        budget: RodRigidGraphBudget,
        gate: &CancelGate,
    ) -> Result<(Self, Vec<RigidBodyState>, RodRigidFrictionalGraphReport), RodRigidFrictionalError>
    {
        checkpoint(gate)?;
        let count = witnesses.len();
        if count > 21
            || count > controls.contact.maximum_contacts
            || 3 * count > controls.contact.matrix.max_variables
        {
            return Err(RodRigidImpactError::Invalid("frictional graph contact budget").into());
        }
        validate_graph_budget(bodies.len(), budget)?;
        let edges: Vec<_> = witnesses.iter().map(|witness| witness.geometry).collect();
        let geometry = self.prepare_contact_graph(bodies, &edges, controls.impact, budget, gate)?;
        let mut frames = Vec::with_capacity(count);
        let mut contacts = Vec::with_capacity(count);
        for (i, witness) in witnesses.iter().enumerate() {
            checkpoint(gate)?;
            let point = &geometry.contacts[i];
            let frame = contact_frame(point.normal_world, witness.tangent_world)?;
            let mut before = [0.0; 3];
            for a in 0..3 {
                before[a] = project_velocity(frame[a], point.relative_velocity_world_m_s)?;
            }
            frames.push(frame);
            contacts.push(FrictionalContact {
                id: i as u64,
                gap: point.gap_m,
                activation_tolerance: witness.geometry.activation_tolerance_m,
                relative_velocity: before,
                friction_coefficient: witness.friction_coefficient,
                normal_law: NormalVelocityLaw::NewtonRestitution {
                    coefficient: witness.geometry.restitution,
                },
            });
        }
        let delassus = self.contact_graph_mobility(bodies, &geometry, &frames, gate)?;
        let contact = solve_frictional_contacts(&contacts, &delassus, controls.contact, || {
            gate.is_requested()
        })?;
        let (rod, states, transfer) = self.stage_contact_graph(
            bodies,
            &geometry,
            &frames,
            &contact.impulses,
            &contact.relative_velocities_after,
            contact.dissipated_energy,
            contact.relative_work,
            controls.impact,
            gate,
        )?;
        let report = RodRigidFrictionalGraphReport {
            contact,
            delassus_kg_inverse: delassus,
            gaps_m: contacts.iter().map(|point| point.gap).collect(),
            frames_world: frames,
            relative_velocities_before_m_s: contacts
                .iter()
                .map(|point| point.relative_velocity)
                .collect(),
            transfer,
        };
        checkpoint(gate)?;
        Ok((rod, states, report))
    }
}
