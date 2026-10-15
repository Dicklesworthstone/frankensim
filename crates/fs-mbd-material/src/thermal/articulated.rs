//! Apply complete thermal configurations to the existing material mechanism owner.
use std::collections::BTreeMap;
use fs_mbd::articulated::{ArticulatedModel, Link};
use crate::articulated::ResolvedArticulated;
use crate::{MaterialInertiaError, poll};
use super::{ExpandedBody, ThermalExpansionBinding, expand_body};

/// Complete region temperature paths for one named reference link.
#[derive(Clone, Copy, Debug)]
pub struct LinkThermalExpansion<'a> {
    /// Exact existing link name. All links, including fixed links, are required.
    pub link: &'a str,
    /// Exactly one source-bound expansion per region of this link.
    pub regions: &'a [ThermalExpansionBinding<'a>],
}

/// Thermally configured mechanism with original material inventories retained.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedArticulated<'a> {
    reference: &'a ResolvedArticulated,
    model: ArticulatedModel,
    bodies: Vec<ExpandedBody<'a>>,
}
impl<'a> ExpandedArticulated<'a> {
    /// Unchanged reference mechanism and its complete material-source bundles.
    pub const fn reference(&self) -> &'a ResolvedArticulated { self.reference }
    /// Feed directly to the existing kinematics and articulated dynamics solvers.
    pub const fn model(&self) -> &ArticulatedModel { &self.model }
    /// Thermal geometry and receipts in the original link-tree order.
    pub fn bodies(&self) -> &[ExpandedBody<'a>] { &self.bodies }
}

/// Reconfigure all link inertias from source-integrated free thermal expansion.
///
/// Exactly one complete binding per reference link is required. Binding order
/// cannot change the tree, generalized-coordinate order, joints, limits or
/// body frames. Every link's retained mass inventory is conserved and its
/// thermal-source reference is checked by `expand_body`; no stale hand-entered
/// inertia survives a missing assignment. The whole result is transactional.
///
/// Joint placements remain PRESCRIBED. This does not infer thermally displaced
/// joints, solve bonded-region compatibility, or perform a momentum/energy
/// conserving live-state transition. The output is an immutable configured
/// model for the existing solver, with its geometry and source receipts.
/// Allocation and work are bounded by the reference mechanism's admitted global
/// region count. Cancellation is checked before work, per link and before return.
pub fn expand_articulated<'a>(
    reference: &'a ResolvedArticulated,
    bindings: &[LinkThermalExpansion<'_>],
    mut cancelled: impl FnMut() -> bool,
) -> Result<ExpandedArticulated<'a>, MaterialInertiaError> {
    poll(&mut cancelled)?;
    if bindings.len() != reference.model().link_count() {
        return Err(MaterialInertiaError::Invalid("thermal assignments must cover every reference link exactly once"));
    }
    let mut by_name = BTreeMap::new();
    for binding in bindings {
        poll(&mut cancelled)?;
        if binding.link.trim().is_empty() || binding.link.len() > 128
            || by_name.insert(binding.link, binding.regions).is_some() {
            return Err(MaterialInertiaError::Invalid("thermal link names must be bounded, nonblank and unique"));
        }
    }
    let mut numerical = Vec::with_capacity(reference.model().link_count());
    let mut bodies = Vec::with_capacity(reference.bodies().len());
    for (link, original_body) in reference.model().links().iter().zip(reference.bodies()) {
        poll(&mut cancelled)?;
        let regions = by_name.remove(link.name())
            .ok_or(MaterialInertiaError::Invalid("missing or foreign thermal link assignment"))?;
        let body = expand_body(original_body, regions, &mut cancelled)?;
        numerical.push(Link::new(link.name(), link.parent(), link.parent_from_child_zero(),
            link.joint(), body.inertia()));
        bodies.push(body);
    }
    let model = ArticulatedModel::new(numerical)?;
    poll(&mut cancelled)?;
    Ok(ExpandedArticulated { reference, model, bodies })
}
