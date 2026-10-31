//! Apply complete thermal configurations to the existing material mechanism owner.
use super::{ExpandedBody, ThermalExpansionBinding, expand_body};
use crate::articulated::ResolvedArticulated;
use crate::{MaterialInertiaError, poll};
use fs_ga::{Se3, Vec3};
use fs_mbd::articulated::{ArticulatedError, ArticulatedModel, Link};
use std::collections::BTreeMap;

/// Complete region temperature paths for one named reference link.
#[derive(Clone, Copy, Debug)]
pub struct LinkThermalExpansion<'a> {
    /// Exact existing link name. All links, including fixed links, are required.
    pub link: &'a str,
    /// Exactly one source-bound expansion per region of this link.
    pub regions: &'a [ThermalExpansionBinding<'a>],
}

/// A joint origin attached to an explicit material region on each side.
///
/// The original zero-placement translation is the reference attachment point
/// in the parent frame; the child-side point is its body-frame origin. The
/// parent point follows its region's source-integrated expansion. The selected
/// child region must expand about its origin, so the native joint axis and
/// motion subspace retain their meaning without an invented frame shift.
/// Region membership and attachment compatibility are caller geometry inputs.
/// Unlisted joints retain their explicitly prescribed reference placement.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct JointMaterialAttachment {
    /// Exact existing non-root link name.
    pub link: String,
    /// Region key in this link's immediate parent.
    pub parent_region: String,
    /// Region key in this link, anchored at its body-frame origin.
    pub child_region: String,
}

/// Thermally configured mechanism with original material inventories retained.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedArticulated<'a> {
    reference: &'a ResolvedArticulated,
    model: ArticulatedModel,
    bodies: Vec<ExpandedBody<'a>>,
    joint_attachments: Vec<JointMaterialAttachment>,
}
impl<'a> ExpandedArticulated<'a> {
    /// Unchanged reference mechanism and its complete material-source bundles.
    pub const fn reference(&self) -> &'a ResolvedArticulated { self.reference }
    /// Feed directly to the existing kinematics and articulated dynamics solvers.
    pub const fn model(&self) -> &ArticulatedModel { &self.model }
    /// Thermal geometry and receipts in the original link-tree order.
    pub fn bodies(&self) -> &[ExpandedBody<'a>] { &self.bodies }
    /// Material-driven joints in original link order; others remain prescribed.
    pub fn joint_attachments(&self) -> &[JointMaterialAttachment] {
        &self.joint_attachments
    }

    /// Replace the optional joint bindings with the SAME material deformation
    /// as link inertia; unlisted placements revert to their reference values.
    /// Always rebuilds from original placements, never a previously expanded
    /// translation. This is an immutable configuration, not a live-state remap;
    /// a continuing motion owner must account for changed momentum and work.
    pub fn with_joint_attachments(
        &self,
        attachments: &[JointMaterialAttachment],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MaterialInertiaError> {
        poll(&mut cancelled)?;
        Self::from_bodies_with_joint_attachments(
            self.reference,
            self.bodies.clone(),
            attachments,
            &mut cancelled,
        )
    }

    // Shared assembly for prescribed temperatures and actual heat-transport
    // endpoints. The original inventory, body order and topology stay bound.
    pub(super) fn from_bodies(
        reference: &'a ResolvedArticulated,
        bodies: Vec<ExpandedBody<'a>>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MaterialInertiaError> {
        Self::from_bodies_with_joint_attachments(reference, bodies, &[], cancelled)
    }

    pub(super) fn from_bodies_with_joint_attachments(
        reference: &'a ResolvedArticulated,
        bodies: Vec<ExpandedBody<'a>>,
        attachments: &[JointMaterialAttachment],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MaterialInertiaError> {
        poll(&mut cancelled)?;
        if bodies.len() != reference.bodies().len() {
            return Err(MaterialInertiaError::Invalid("configured bodies must cover the reference mechanism"));
        }
        if attachments.len() >= bodies.len() {
            return Err(MaterialInertiaError::Invalid(
                "joint attachments may only name non-root links",
            ));
        }
        let mut by_name = BTreeMap::new();
        for attachment in attachments {
            poll(&mut cancelled)?;
            if [
                &attachment.link,
                &attachment.parent_region,
                &attachment.child_region,
            ]
            .into_iter()
            .any(|s| s.trim().is_empty() || s.len() > 128)
                || by_name
                    .insert(attachment.link.as_str(), attachment)
                    .is_some()
            {
                return Err(MaterialInertiaError::Invalid(
                    "joint attachment names must be bounded, nonblank and unique",
                ));
            }
        }
        let mut numerical = Vec::with_capacity(bodies.len());
        let mut joint_attachments = Vec::with_capacity(attachments.len());
        for ((link, original), body) in reference.model().links().iter()
            .zip(reference.bodies()).zip(&bodies) {
            poll(&mut cancelled)?;
            if !core::ptr::eq(original, body.reference()) {
                return Err(MaterialInertiaError::Invalid("configured body belongs to a different reference link"));
            }
            let original_placement = link.parent_from_child_zero();
            let placement = if let Some(attachment) = by_name.remove(link.name()) {
                let parent = link.parent().ok_or(MaterialInertiaError::Invalid(
                    "the root placement belongs to the external base, not a parent material"))?;
                let parent_region = region(&bodies[parent], &attachment.parent_region)?;
                let child_region = region(body, &attachment.child_region)?;
                if child_region.fixed_point_m() != Vec3::default() {
                    return Err(MaterialInertiaError::Invalid(
                        "child joint material must expand about its body-frame origin"));
                }
                let point = parent_region.transport_reference_point(original_placement.translation())?;
                joint_attachments.push(attachment.clone());
                if point == original_placement.translation() {
                    original_placement
                } else {
                    Se3::from_parts(original_placement.rotation(), point)
                        .map_err(ArticulatedError::from)?
                }
            } else {
                original_placement
            };
            numerical.push(Link::new(link.name(), link.parent(), placement,
                link.joint(), body.inertia()));
        }
        if !by_name.is_empty() {
            return Err(MaterialInertiaError::Invalid(
                "joint attachment names an unknown link",
            ));
        }
        let model = ArticulatedModel::new(numerical)?;
        poll(&mut cancelled)?;
        Ok(Self {
            reference,
            model,
            bodies,
            joint_attachments,
        })
    }
}

fn region<'a, 'b>(
    body: &'b ExpandedBody<'a>,
    key: &str,
) -> Result<&'b super::ExpandedRegion, MaterialInertiaError> {
    let index = body
        .reference()
        .regions()
        .iter()
        .position(|r| r.key() == key)
        .ok_or(MaterialInertiaError::Invalid(
            "joint attachment names an unknown material region",
        ))?;
    Ok(&body.regions()[index])
}

/// Reconfigure all link inertias from source-integrated free thermal expansion.
///
/// Exactly one complete binding per reference link is required. Binding order
/// cannot change the tree, generalized-coordinate order, joints, limits or
/// body frames. Every link's retained mass inventory is conserved and its
/// thermal-source reference is checked by `expand_body`; no stale hand-entered
/// inertia survives a missing assignment. The whole result is transactional.
///
/// Joint placements remain prescribed unless `with_joint_attachments` explicitly
/// binds them to material geometry. This does not solve bonded-region compatibility
/// or perform a momentum/energy conserving live-state transition. The output
/// is an immutable configured model, with its geometry and source receipts.
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
    let mut bodies = Vec::with_capacity(reference.bodies().len());
    for (link, original_body) in reference.model().links().iter().zip(reference.bodies()) {
        poll(&mut cancelled)?;
        let regions = by_name.remove(link.name())
            .ok_or(MaterialInertiaError::Invalid("missing or foreign thermal link assignment"))?;
        let body = expand_body(original_body, regions, &mut cancelled)?;
        bodies.push(body);
    }
    ExpandedArticulated::from_bodies(reference, bodies, cancelled)
}
