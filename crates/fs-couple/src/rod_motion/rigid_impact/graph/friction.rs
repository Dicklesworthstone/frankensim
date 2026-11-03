//! Coupled Coulomb impulses on the original mixed rod/rigid contact graph.
//!
//! Endpoint geometry, signed mobility and native momentum publication belong to
//! the common graph owner. This adapter only selects the original three-axis
//! contact law; normal and tangent impulses share the same graph response.

use super::super::RodRigidImpactError;
use super::super::friction::{RodRigidFrictionalControls, RodRigidFrictionalError};
use super::network::friction::{
    RodRigidFrictionalNetworkContact, apply_rigid_frictional_contact_network,
};
use super::network::{RodRigidNetworkBudget, validate_network_budget};
use super::{
    RodRigidGraphBudget, RodRigidGraphContact, RodRigidGraphTransfer, RodRigidParticipant,
};
use crate::rod_motion::{RodMotion, RodMotionInertia, checkpoint};
use fs_contact::frictional_multicontact::FrictionalContactResult;
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
        let budget = RodRigidNetworkBudget {
            maximum_rods: 1,
            maximum_bodies: budget.maximum_bodies,
            maximum_rod_nodes: self.model().inertia().len(),
        };
        validate_network_budget(core::slice::from_ref(self), bodies.len(), budget)?;
        let contacts: Vec<_> = witnesses
            .iter()
            .map(|witness| RodRigidFrictionalNetworkContact {
                geometry: witness.geometry.into(),
                tangent_world: witness.tangent_world,
                friction_coefficient: witness.friction_coefficient,
            })
            .collect();
        let (mut rods, states, report) = apply_rigid_frictional_contact_network(
            core::slice::from_ref(self),
            bodies,
            &contacts,
            controls,
            budget,
            gate,
        )?;
        let rod = rods.pop().ok_or(RodRigidImpactError::Invalid(
            "single-rod contact result missing its participant",
        ))?;
        let report = RodRigidFrictionalGraphReport {
            contact: report.contact,
            delassus_kg_inverse: report.delassus_kg_inverse,
            gaps_m: report.gaps_m,
            frames_world: report.frames_world,
            relative_velocities_before_m_s: report.relative_velocities_before_m_s,
            transfer: report.transfer.into_single_rod(),
        };
        checkpoint(gate)?;
        Ok((rod, states, report))
    }
}
