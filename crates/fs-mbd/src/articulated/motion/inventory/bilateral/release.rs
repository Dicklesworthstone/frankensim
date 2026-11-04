//! Release and ejection through the original mixed rigid-constraint response.

use super::*;
use crate::articulated::motion::inventory::release::{
    PayloadBatchRelease, PayloadRelease, PayloadReleaseInput, SeparationImpulse,
};

/// Original two-body release with the retained tree's simultaneous reaction.
#[derive(Debug)]
pub struct BilateralPayloadRelease {
    /// Actual surviving states, inventory identity and common-frame accounting.
    pub release: PayloadRelease,
    /// Retained-tree event, absent for exact legacy empty-loop delegation.
    /// Its source work concerns the retained tree alone in the lab frame.
    /// The release receipt pairs both sides of separation in its common frame.
    pub reaction: Option<BilateralImpulseResponse>,
}

/// One simultaneous release retaining all material and ideal joint constraints.
#[derive(Debug)]
pub struct BilateralPayloadBatchRelease {
    /// All actual surviving bodies and their complete common-frame accounting.
    pub release: PayloadBatchRelease,
    /// Auxiliary retained-tree response to the summed separation source and
    /// simultaneous row reactions, absent for exact empty-loop delegation.
    pub reaction: Option<BilateralImpulseResponse>,
}

/// Detach several payloads through one closed mechanism at one physical instant.
///
/// All retained inertias are restored before a single mixed native impulse
/// response transmits the accumulated point impulses and couples. Every
/// detached body inherits the original loaded velocity field, then responds
/// through its own full inertia. Pair work uses the same final retained state,
/// preserving recoil cross terms between different separation sources.
///
/// On repeated targets, inputs must consume handles in reverse capture order;
/// this orders inventory validation only. The physical event is simultaneous.
/// Passive or exactly cancelling retained-side sources preserve state bits
/// after geometry, tangent, rank and condition admission. The complete signed
/// separation work and actual constraint work keep their independent limits.
/// Empty material loops delegate to the original batch release exactly.
///
/// # Errors
/// Empty/stale/misordered inventory, unresolved constraints, native endpoint,
/// momentum/work/energy budgets or cancellation refuse atomically.
#[allow(clippy::too_many_arguments)]
pub fn release_payloads(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    inputs: &[PayloadReleaseInput<'_>],
    loops: &[BilateralLoop],
    rows: &[JointConstraint],
    controls: PayloadControls<'_>,
    config: BilateralImpulseConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<BilateralPayloadBatchRelease, PayloadError> {
    if loops.is_empty() {
        return super::super::release::release_payloads(
            model, state, inputs, rows, controls, cancelled,
        )
        .map(|release| BilateralPayloadBatchRelease {
            release,
            reaction: None,
        });
    }
    super::super::release::release_payloads_evaluated(
        model,
        state,
        inputs,
        rows,
        controls,
        Some((loops, config)),
        cancelled,
    )
    .map(|(release, reaction)| BilateralPayloadBatchRelease { release, reaction })
}

/// Detach retained inventory while material closures and gears transmit recoil.
///
/// The original attachment restores the exact retained inertia. One native
/// event applies its separation impulse and all body/joint reactions together;
/// the detached payload responds through its original full spatial inertia.
/// Physical endpoint speed limits apply to that combined response.
///
/// Passive and zero-impulse release admit the retained mechanism's geometry,
/// tangent tolerances, combined rank and condition policy, then preserve its
/// state bits. Finite admitted defects do not activate a projection. An active
/// separation retains actual signed constraint work separately from the
/// equal/opposite source pair, both in the existing common translating frame.
///
/// The existing absolute source-work and constraint-work caps remain
/// independent. This does not infer fracture, engagement or a contact law.
/// Empty material loops reproduce the original release path exactly.
///
/// # Errors
/// Changed inventory, invalid state or supplied impulse, unresolved rows,
/// physical limits, momentum/work/energy budgets, or cancellation refuse
/// without changing the input model, state, attachment or caller's budgets.
#[allow(clippy::too_many_arguments)]
pub fn release_payload(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    attachment: &PayloadAttachment,
    impulse: Option<SeparationImpulse>,
    loops: &[BilateralLoop],
    rows: &[JointConstraint],
    controls: PayloadControls<'_>,
    config: BilateralImpulseConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<BilateralPayloadRelease, PayloadError> {
    if loops.is_empty() {
        return super::super::release::release_payload(
            model, state, attachment, impulse, rows, controls, cancelled,
        )
        .map(|release| BilateralPayloadRelease {
            release,
            reaction: None,
        });
    }
    super::super::release::release_evaluated(
        model,
        state,
        attachment,
        impulse,
        rows,
        controls,
        Some((loops, config)),
        cancelled,
    )
    .map(|(release, reaction)| BilateralPayloadRelease { release, reaction })
}
