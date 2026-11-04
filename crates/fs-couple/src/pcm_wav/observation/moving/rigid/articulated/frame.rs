//! Atomic accepted frames for multiple articulated acoustic attachments.

use super::{
    ArticulatedAcousticError, BodyFrameCompactSourceStrength, PreparedLinkPointMotion,
    RigidAcousticHistory, RigidAcousticSample,
};
use super::super::{invalid, poll};
use super::super::relative::RelativeAcousticMotion;

/// One source/history binding for an accepted mechanical event.
///
/// Multiple entries may borrow the SAME prepared link snapshot. Distinct
/// mutable history borrows prevent recording one history twice in a frame.
/// A microphone may use zero strength while retaining its original motion.
pub struct LinkAcousticUpdate<'a> {
    /// Original preallocated history, retaining its fixed link-relative anchor.
    pub history: &'a mut RigidAcousticHistory,
    /// Original native receipt for this attachment's link at the supplied time.
    pub motion: &'a PreparedLinkPointMotion,
    /// Explicit force ON FLUID and volume injection in this link's axes.
    pub strength: BodyFrameCompactSourceStrength,
}

/// Publish one common physical timestamp across all supplied fixed attachments.
///
/// Stage every original point/force sample in caller-owned scratch, then
/// publish ALL paired histories only after the final cancellation checkpoint.
/// Any history/source/motion refusal changes none of the histories or mechanics.
/// Scratch is disposable on failure and is overwritten completely on a retry;
/// its incoming contents can never authorize publication.
///
/// Provide exactly one scratch slot per update. There is no heap allocation,
/// history clone, repeated dynamics solve, or hidden sample clock. Work and
/// scratch are O(attachments); the final non-interruptible publication is one
/// paired push per attachment into already reserved capacity. Callers bound
/// commit latency by the size of this explicit frame, not by weakening atomicity.
/// Mechanical identity, world frame and physical source provenance remain with
/// the caller, exactly as for `RigidAcousticHistory::append_link`.
pub fn append_link_frame(
    time_s: f64,
    updates: &mut [LinkAcousticUpdate<'_>],
    scratch: &mut [Option<RigidAcousticSample>],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ArticulatedAcousticError> {
    append_frame(time_s, updates, None, scratch, cancelled)
}

/// Atomically record fixed and migrating sources/microphones in the same frame.
///
/// Supply one local displacement/velocity/acceleration record per update;
/// a zero record uses the unchanged fixed-attachment path. All records describe
/// this same accepted physical event. Source derivatives must include any
/// variation along the moving local path, in the carrier's coordinate axes.
/// Shape, motion, source, clock and cancellation refusals publish nothing.
/// This reuses the original frame staging, storage and publication owner.
pub fn append_link_frame_with_relative_motion(
    time_s: f64,
    updates: &mut [LinkAcousticUpdate<'_>],
    relative: &[RelativeAcousticMotion],
    scratch: &mut [Option<RigidAcousticSample>],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ArticulatedAcousticError> {
    append_frame(time_s, updates, Some(relative), scratch, cancelled)
}

fn append_frame(
    time_s: f64,
    updates: &mut [LinkAcousticUpdate<'_>],
    relative: Option<&[RelativeAcousticMotion]>,
    scratch: &mut [Option<RigidAcousticSample>],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ArticulatedAcousticError> {
    poll(cancelled)?;
    if updates.is_empty() || scratch.len() != updates.len()
        || relative.is_some_and(|motion| motion.len() != updates.len())
    {
        return Err(invalid("articulated acoustic frame needs matching scratch and motion counts").into());
    }
    for (i, (update, slot)) in updates.iter().zip(scratch.iter_mut()).enumerate() {
        poll(cancelled)?;
        update.history.admit_append(time_s)?;
        *slot = Some(match relative {
            Some(motion) => update.history.attachment.sample_link_relative(
                time_s, update.motion, motion[i], update.strength,
            )?,
            None => update.history.attachment.sample_link(
                time_s, update.motion, update.strength,
            )?,
        });
    }
    poll(cancelled)?;
    // Unique mutable history borrows and successful preflight guarantee all
    // capacities and clocks remain valid. No fallible work follows this point.
    for (update, slot) in updates.iter_mut().zip(scratch.iter_mut()) {
        update.history.publish(slot.take().expect("every frame sample was staged"));
    }
    Ok(())
}
