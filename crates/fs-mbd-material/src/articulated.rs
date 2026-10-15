//! Material/geometry descriptions lowered to the existing articulated-body owner.
use std::collections::BTreeMap;

use fs_ga::Se3;
use fs_mbd::articulated::{ArticulatedModel, JointModel, Link};

use crate::{MAX_BODY_REGIONS, MaterialInertiaError, RegionBinding, ResolvedBody, poll, resolve_body};

/// One material-defined link, in root-first tree order.
/// No dummy mass or hand-entered inertia is required to build a mechanism.
#[derive(Clone, Debug)]
pub struct MaterialLink<'a> {
    /// Unique link name; preserved verbatim by the existing model owner.
    pub name: &'a str,
    /// Index of an earlier parent, or None for the sole root at index zero.
    pub parent: Option<usize>,
    /// Existing zero-coordinate child-to-parent placement.
    pub parent_from_child_zero: Se3,
    /// Existing joint type, axis and physical limits; no new joint schema.
    pub joint: JointModel,
    /// Non-overlapping rigid mass regions in THIS link's body frame.
    pub regions: &'a [RegionBinding<'a>],
}

/// An existing articulated model with all material/geometry source bundles retained.
/// Bodies and links have exactly the same order. Topology and generalized state
/// coordinates stay owned by `fs-mbd`, not by this compilation boundary.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedArticulated {
    model: ArticulatedModel,
    bodies: Vec<ResolvedBody>,
}
impl ResolvedArticulated {
    /// Feed directly to existing forward, inverse, kinematic and free-flight solvers.
    pub const fn model(&self) -> &ArticulatedModel { &self.model }
    /// Complete per-link source-bound bodies in original tree order.
    pub fn bodies(&self) -> &[ResolvedBody] { &self.bodies }
}

/// Compile a complete mechanism from geometry and material assignments.
///
/// `max_regions` is a GLOBAL budget across all links, not a per-link multiplier.
/// Every link requires at least one region. Joint, frame, name and parent
/// semantics are checked by the existing articulated model constructor.
/// Source selection, state-domain refusal and cancellation are inherited from
/// `resolve_body`; any failure leaves no partly updated mechanism.
///
/// This constructs a model, not an integrator. Applied efforts, gravity,
/// external wrenches and q/qdot remain explicit inputs to existing solvers.
/// Creating another material model does not account for the mass, momentum or
/// energy of physically exchanging material in an already moving mechanism.
pub fn resolve_articulated(
    links: &[MaterialLink<'_>],
    max_regions: usize,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ResolvedArticulated, MaterialInertiaError> {
    poll(&mut cancelled)?;
    if max_regions == 0 || max_regions > MAX_BODY_REGIONS
        || links.is_empty() || links.len() > max_regions {
        return Err(MaterialInertiaError::Invalid("links need a positive global region budget within the 4096 cap"));
    }
    let mut count = 0usize;
    for link in links {
        if link.name.trim().is_empty() || link.name.len() > 128 || link.regions.is_empty() {
            return Err(MaterialInertiaError::Invalid("each link needs a bounded nonblank name and at least one material region"));
        }
        count = count.checked_add(link.regions.len())
            .filter(|&n| n <= max_regions)
            .ok_or(MaterialInertiaError::Invalid("total mechanism regions exceed the global budget"))?;
    }
    let mut bodies = Vec::with_capacity(links.len());
    let mut numerical = Vec::with_capacity(links.len());
    for link in links {
        poll(&mut cancelled)?;
        let body = resolve_body(link.regions, link.regions.len(), &mut cancelled)?;
        numerical.push(Link::new(link.name, link.parent, link.parent_from_child_zero,
            link.joint, body.inertia()));
        bodies.push(body);
    }
    let model = ArticulatedModel::new(numerical)?;
    poll(&mut cancelled)?;
    Ok(ResolvedArticulated { model, bodies })
}

/// Explicit material replacement for one named link of an existing mechanism.
#[derive(Clone, Debug)]
pub struct LinkMaterials<'a> {
    /// Exact existing link name, not an inferred material/chemistry name.
    pub link: &'a str,
    /// Complete new mass partition in that link's original body frame.
    pub regions: &'a [RegionBinding<'a>],
}

/// Recompile ALL link inertias while preserving the existing mechanism topology.
///
/// Exactly one binding per existing link is required. Missing, repeated and
/// foreign names refuse rather than silently retaining a stale authored mass.
/// Binding input order is immaterial; original link order, names, parents,
/// frames, joint axes and limits are retained. Numerical dynamics and source
/// bundles are published together only after complete successful compilation.
/// The immutable template is never modified, including on late cancellation.
pub fn rebind_articulated(
    template: &ArticulatedModel,
    bindings: &[LinkMaterials<'_>],
    max_regions: usize,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ResolvedArticulated, MaterialInertiaError> {
    poll(&mut cancelled)?;
    if max_regions == 0 || max_regions > MAX_BODY_REGIONS
        || bindings.len() != template.link_count() || bindings.len() > max_regions {
        return Err(MaterialInertiaError::Invalid("exactly one material binding per link and a bounded global region budget required"));
    }
    let mut by_name = BTreeMap::new();
    for binding in bindings {
        if binding.link.trim().is_empty() || binding.link.len() > 128
            || by_name.insert(binding.link, binding.regions).is_some() {
            return Err(MaterialInertiaError::Invalid("link material names must be bounded, nonblank and unique"));
        }
    }
    let mut links = Vec::with_capacity(template.link_count());
    for link in template.links() {
        let regions = by_name.remove(link.name())
            .ok_or(MaterialInertiaError::Invalid("missing link material assignment or foreign link name"))?;
        links.push(MaterialLink { name: link.name(), parent: link.parent(),
            parent_from_child_zero: link.parent_from_child_zero(), joint: link.joint(), regions });
    }
    resolve_articulated(&links, max_regions, &mut cancelled)
}
