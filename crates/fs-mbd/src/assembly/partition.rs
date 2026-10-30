//! Partial separation: release stages without disassembling their payloads.
//!
//! A prepared partition caches each retained group's complete inertia and
//! original component geometry. Impulses identify original components, but act
//! on the full rigid group containing them. This retains internal rigidity and
//! includes all shared-group impulse-work cross terms. Releasing every leaf
//! first and then relatching groups would introduce a fictitious inelastic
//! event and the wrong separation-work account.
//!
//! Preparation and storage are O(parts + groups); each event is O(groups +
//! impulses). Expanding evolved groups back to their original material bodies
//! is O(parts), using two reusable largest-group workspaces. No dense multibody
//! matrix or additional physical solver is introduced.

use super::{
    AssemblyControls, AssemblyError, AssemblyPart, AssemblyRelease, ImpulseTotals, RigidAssembly,
    SeparationImpulse, poll, release_components_into,
};
use crate::{Pose, RigidBodyState, UnitQuaternion, Vec3};

/// A deterministic grouping of every original component into one retained
/// rigid subassembly. Group IDs are contiguous from zero; their geometry and
/// mass properties are prepared once and reused for separation events.
#[derive(Clone, Debug, PartialEq)]
pub struct AssemblyPartition {
    grouped: RigidAssembly,
    groups: Vec<RigidAssembly>,
    membership: Vec<usize>,
    offsets: Vec<usize>,
    original_indices: Vec<usize>,
}

impl AssemblyPartition {
    /// Prepare partial separation of `assembly`. `membership[i]` selects the
    /// retained group of original part `i`; every ID from zero to the largest
    /// supplied ID must occur. Original order is retained within each group.
    ///
    /// The caller owns connectivity and geometric admissibility, as for the
    /// original assembly owner. This operation neither moves bodies nor admits
    /// a contact, fracture, elastic-storage, or attachment-strength model.
    pub fn new(
        assembly: &RigidAssembly,
        membership: &[usize],
        controls: AssemblyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AssemblyError> {
        poll(&mut cancelled)?;
        let count = assembly.parts.len();
        controls.validate(count, 0)?;
        if membership.len() != count {
            return Err(AssemblyError::Invalid(
                "one group ID per original part required",
            ));
        }
        let mut sizes = vec![0; count];
        let mut group_count = 0;
        for &group in membership {
            poll(&mut cancelled)?;
            if group >= count {
                return Err(AssemblyError::Invalid(
                    "group ID exceeds original part count",
                ));
            }
            sizes[group] += 1;
            group_count = group_count.max(group + 1);
        }
        let mut offsets = Vec::with_capacity(group_count + 1);
        offsets.push(0);
        for &size in &sizes[..group_count] {
            poll(&mut cancelled)?;
            if size == 0 {
                return Err(AssemblyError::Invalid(
                    "group IDs must be contiguous and nonempty",
                ));
            }
            offsets.push(offsets.last().copied().unwrap_or(0) + size);
        }
        let mut cursors = offsets[..group_count].to_vec();
        let mut original_indices = vec![0; count];
        for (index, &group) in membership.iter().enumerate() {
            poll(&mut cancelled)?;
            original_indices[cursors[group]] = index;
            cursors[group] += 1;
        }
        let mut groups = Vec::with_capacity(group_count);
        let mut group_parts = Vec::with_capacity(group_count);
        let largest = sizes[..group_count].iter().copied().max().unwrap_or(0);
        let mut parts = Vec::with_capacity(largest);
        for group in 0..group_count {
            poll(&mut cancelled)?;
            parts.clear();
            for &index in &original_indices[offsets[group]..offsets[group + 1]] {
                poll(&mut cancelled)?;
                parts.push(assembly.parts[index]);
            }
            let child = RigidAssembly::new(&parts, controls, &mut cancelled)?;
            // Every group retains the parent's reference axes. Each original
            // component keeps its own orientation within the child assembly.
            group_parts.push(AssemblyPart {
                properties: child.properties(),
                pose_in_reference: Pose::new(
                    child.center_of_mass_in_reference(),
                    UnitQuaternion::IDENTITY,
                )?,
            });
            groups.push(child);
        }
        poll(&mut cancelled)?;
        Ok(Self {
            grouped: RigidAssembly {
                parts: group_parts,
                properties: assembly.properties,
                center_in_reference: Vec3::ZERO,
            },
            groups,
            membership: membership.to_vec(),
            offsets,
            original_indices,
        })
    }

    /// Cached child assemblies in group-ID order. Their COM-relative geometry
    /// can itself be partitioned for later staging or released normally.
    pub fn groups(&self) -> &[RigidAssembly] {
        &self.groups
    }

    /// Group ID for each original component, in original component order.
    pub fn membership(&self) -> &[usize] {
        &self.membership
    }

    /// Original component indices corresponding to a child's local part order.
    /// Returns `None` when the group does not exist.
    pub fn original_part_indices(&self, group: usize) -> Option<&[usize]> {
        let &start = self.offsets.get(group)?;
        let &end = self.offsets.get(group.checked_add(1)?)?;
        Some(&self.original_indices[start..end])
    }

    /// Release the prepared groups with the original parent's COM state.
    /// Impulse indices refer to ORIGINAL components, and their physical points
    /// remain relative to the PARENT assembly COM in its axes. Each impulse
    /// must cross a released interface: endpoints within one retained group
    /// refuse, since a rigid model cannot describe their internal deformation.
    ///
    /// Returned states are in group-ID order. Each carries the entire group's
    /// inherited spin/orbital motion and net impulse. The shared native release
    /// owner checks actual world momentum and supplied-work balances.
    pub fn release(
        &self,
        parent: RigidBodyState,
        impulses: &[SeparationImpulse],
        controls: AssemblyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AssemblyRelease, AssemblyError> {
        poll(&mut cancelled)?;
        controls.validate(self.membership.len(), impulses.len())?;
        let mut grouped_impulses = Vec::with_capacity(impulses.len());
        for impulse in impulses {
            poll(&mut cancelled)?;
            let Some((&part_a, &part_b)) = self
                .membership
                .get(impulse.part_a)
                .zip(self.membership.get(impulse.part_b))
            else {
                return Err(AssemblyError::Invalid(
                    "separation endpoint is not an original part",
                ));
            };
            if part_a == part_b {
                return Err(AssemblyError::Invalid(
                    "separation impulse must cross retained groups",
                ));
            }
            grouped_impulses.push(SeparationImpulse {
                part_a,
                part_b,
                ..*impulse
            });
        }
        self.grouped
            .release(parent, &grouped_impulses, controls, cancelled)
    }

    /// Recover all original material-body states from the groups' current
    /// states, in original component order. Groups may have translated and
    /// rotated independently since separation. Every body inherits its group's
    /// full rigid velocity field; there is no new impulse or inelastic loss.
    ///
    /// Each actual child expansion must pass the original momentum/energy
    /// gates. All staging is private until every group succeeds, including late
    /// cancellation. Caller states and this prepared partition remain unchanged.
    pub fn expand(
        &self,
        group_states: &[RigidBodyState],
        controls: AssemblyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<RigidBodyState>, AssemblyError> {
        poll(&mut cancelled)?;
        controls.validate(self.membership.len(), 0)?;
        if group_states.len() != self.groups.len() {
            return Err(AssemblyError::Invalid(
                "one state per retained group required",
            ));
        }
        let first = group_states[0];
        let largest = self.groups.iter().map(|g| g.parts.len()).max().unwrap_or(0);
        let mut bodies = vec![first; self.membership.len()];
        let mut staged = vec![first; largest];
        let mut deltas = vec![ImpulseTotals::default(); largest];
        for (group, (assembly, &state)) in self.groups.iter().zip(group_states).enumerate() {
            poll(&mut cancelled)?;
            let count = assembly.parts.len();
            let receipt = release_components_into(
                &assembly.parts,
                assembly.properties,
                state,
                &[],
                &mut deltas[..count],
                &mut staged[..count],
                &mut cancelled,
            )?;
            controls.admit(receipt.balance)?;
            for (&index, &body) in self.original_indices
                [self.offsets[group]..self.offsets[group + 1]]
                .iter()
                .zip(&staged[..count])
            {
                poll(&mut cancelled)?;
                bodies[index] = body;
            }
        }
        poll(&mut cancelled)?;
        Ok(bodies)
    }
}
