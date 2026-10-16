//! STABLE PERSISTENT ENTITY IDENTITY (the R3 amendment, bead lmp4.10):
//! patches, interfaces, and load cases carry first-class IDs that
//! SURVIVE EDITS — assigned at creation, transformed EXPLICITLY by
//! ledgered operations, never reconstructed heuristically. This is the
//! structural mitigation of the classic topological-naming problem:
//! FrankenSim controls its own kernel, so identity is an invariant the
//! kernel maintains rather than a guess a differ recovers.

use std::collections::BTreeMap;

/// A stable persistent entity id (patches, interfaces, load cases…).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EntityId(pub u64);

/// How one ledgered edit transformed entity identity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IdTransform {
    /// The entity survives unchanged (identity preserved).
    Preserved(EntityId),
    /// The entity was replaced (e.g. re-fit): old → new.
    Replaced(EntityId, EntityId),
    /// The entity split (e.g. a boolean cut): old → parts.
    Split(EntityId, Vec<EntityId>),
    /// Entities merged: parts → new.
    Merged(Vec<EntityId>, EntityId),
    /// A brand-new entity appeared.
    Created(EntityId),
    /// The entity was deleted.
    Deleted(EntityId),
}

/// The identity map of one edit history: op id → its transforms. Edits
/// are LEDGERED operations, so the map is replayable provenance, not a
/// heuristic reconstruction.
#[derive(Debug, Clone, Default)]
pub struct IdentityMap {
    /// Transforms per op, in op order (BTreeMap: deterministic).
    pub edits: BTreeMap<i64, Vec<IdTransform>>,
}

impl IdentityMap {
    /// An empty map.
    #[must_use]
    pub fn new() -> Self {
        IdentityMap::default()
    }

    /// Record one op's transforms.
    pub fn record(&mut self, op: i64, transforms: Vec<IdTransform>) {
        self.edits.entry(op).or_default().extend(transforms);
    }

    /// The ops (in order) that TOUCHED an entity — directly or through
    /// its replacement/split/merge ancestry. This is the attribution
    /// walk the semantic diff uses.
    #[must_use]
    pub fn ops_touching(&self, entity: EntityId) -> Vec<i64> {
        // Ancestry closure: ids that flow INTO `entity` over history.
        let mut lineage: std::collections::BTreeSet<EntityId> = [entity].into();
        // Walk edits newest→oldest growing the lineage set, then collect
        // touching ops oldest→newest.
        for (_, transforms) in self.edits.iter().rev() {
            for t in transforms {
                match t {
                    IdTransform::Replaced(old, new) if lineage.contains(new) => {
                        lineage.insert(*old);
                    }
                    IdTransform::Split(old, parts) if parts.iter().any(|p| lineage.contains(p)) => {
                        lineage.insert(*old);
                    }
                    IdTransform::Merged(parts, new) if lineage.contains(new) => {
                        lineage.extend(parts.iter().copied());
                    }
                    _ => {}
                }
            }
        }
        let touches = |t: &IdTransform| -> bool {
            match t {
                IdTransform::Preserved(id)
                | IdTransform::Created(id)
                | IdTransform::Deleted(id) => lineage.contains(id),
                IdTransform::Replaced(old, new) => lineage.contains(old) || lineage.contains(new),
                IdTransform::Split(old, parts) => {
                    lineage.contains(old) || parts.iter().any(|p| lineage.contains(p))
                }
                IdTransform::Merged(parts, new) => {
                    lineage.contains(new) || parts.iter().any(|p| lineage.contains(p))
                }
            }
        };
        self.edits
            .iter()
            .filter(|(_, ts)| ts.iter().any(touches))
            .map(|(op, _)| *op)
            .collect()
    }
}
